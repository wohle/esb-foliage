lm_indiv_calc <- function(code_tree_species_input, code_leaves_type_input,
                          dat, params_oi){
  
  # Input parameters to function:
  
  # code_tree_species: int of tree species code according to dict (https://icp-forests.org/documentation/Dictionaries/d_tree_spec_fo.html)
  
  # code_leaves_type: int values of 0 - 4 representing needle age classes;
  # 0: current season, 1: one year-old needles, 2: two year-old needles, etc.;
  # age class of all broadleaved tree species is 0
  
  # dat: data frame containing complete data set
  
  # params_oi: data frame that contains one variable (column) called "params";
  # the row values of params comprise character values of variable (column) names 
  # that are part of the data frame dat
  
  ######################################################################################################|
  ## Prepare data set ####
  
  # Filter for tree species and foliage age class
  dat_species_age <- dat %>% filter(code_tree_species == code_tree_species_input,
                                    code_leaves_type == code_leaves_type_input) %>%
    mutate(key_param_species_ageclass_plot = 
             paste(param, code_tree_species, 
                   code_leaves_type, key_plot, sep = "_"))
  
  # Identify plots with the time series >= 6 years
  n_series_datapoints <- dat_species_age %>% 
    filter(!is.na(param_value_avg)) %>%
    group_by(param, code_tree_species, code_leaves_type, key_plot) %>%
    summarize(n_years_collected = n_distinct(survey_year), 
              year_min = min(survey_year), year_max = max(survey_year)) %>%
    filter(n_years_collected >= 6) %>% #random: filter time series >= 6 years
    mutate(n_years_potential = year_max - year_min + 1) %>%
    mutate(n_years_gaps = n_years_potential - n_years_collected) %>%
    mutate(key_param_species_ageclass_plot = 
             paste(param, code_tree_species, 
                   code_leaves_type, key_plot, sep = "_"))
  
  # Filter for plots with long time series
  dat_species_age <- dat_species_age %>%
    filter(key_param_species_ageclass_plot %in% 
             n_series_datapoints$key_param_species_ageclass_plot)
  
  # Filter out unknown parameters "al" and "b"
  dat_species_age <- dat_species_age %>% filter(param != "al", param != "b")
  
  ######################################################################################################|
  ## Define functions to apply to nested data ####
  
    # Function for centering and scaling
  scaling_fun <- function(df) {
    df %>% mutate(year_centered = survey_year - 2012,
                  var_scaled = 
                    as.numeric(scale(param_value_avg, center = T, scale = T)))
  }
  
  # Linear model of scaled variable values vs. centered year
  param_plot_lm <- function(df) {
    lm(var_scaled ~ year_centered, data = df)
  }
  
  # # Mann-Kendall trend test with scaled variables does not work
  # # -> error message: first mandatory vector should be positive: year or year+fraction
  # # second mandatory vector should be numerical: measured data (ties allowed)
  # # third optional vector should be positive integer: season, month, site or a unique code for season and site
  # # fourth optional vector should be numerical: covariable
  # mk_output <- function(df) {
  #   rkt(as.data.frame(df)$year_centered, as.data.frame(df)$var_scaled)
  # }
  
  # Mann-Kendall trend test
  mk_output <- function(df) {
    rkt(as.data.frame(df)$survey_year, as.data.frame(df)$param_value_avg)
  }
  
  # Shapiro-Wilk Test for testing normality
  sw_test <- function(df) {
    shapiro.test(df$var_scaled)[2] #p value of Shapiro-Wilk test
  }
  
  ######################################################################################################|
  ## Model execution ####
  
  # Group data for parameters and forest plots, nest, scale, model
  dat_model <- dat_species_age %>% group_by(param, key_plot) %>%
    mutate(avg_value_per_plot = mean(param_value_avg, na.rm = T),
           sd_value_per_plot = sd(param_value_avg, na.rm = T)) %>%
    nest() %>%
    mutate(data_scaled = map(data, scaling_fun)) %>%
    mutate(lm_results = map(data_scaled, param_plot_lm)) %>%
    mutate(mk_results = map(data, mk_output)) %>%
    mutate(sw_test_pvalue = map(data_scaled, sw_test))
 
  ######################################################################################################|
  ## Add lm coefficients to data frame ####
  
  # initialize
  dat_model <- dat_model %>%
    mutate(slope_lm = NA, sd_error_slope_lm = NA, intercept_lm = NA, rel_change = NA, 
           abs_change = NA, p_value = NA, R2 = NA, sens_slope = NA, mk_p = NA)
  
  # get coefficients
  for (i in 1:length(dat_model$lm_results)) {
    
    dat_model$slope_lm[i] <- summary(dat_model$lm_results[[i]])$coefficients[2,1]
    
    dat_model$sd_error_slope_lm[i] <- summary(dat_model$lm_results[[i]])$coefficients[2,2]
    
    dat_model$intercept_lm[i] <- summary(dat_model$lm_results[[i]])$coefficients[1,1]
    
    dat_model$rel_change[i] <- #percent (*100) change over 10 years (*10) rel. to y at centered year (backscaled intercept y value)
      10*dat_model$slope_lm[i]*dat_model$data_scaled[[i]]$sd_value_per_plot[1]*100/
      (((summary(dat_model$lm_results[[i]])$coefficients[1,1]*dat_model$data_scaled[[i]]$sd_value_per_plot[1])) + 
         dat_model$data_scaled[[i]]$avg_value_per_plot[1])
    # y_data_scaled = (y_data - mean(y_data))/sd(y_data) -> y_data = (y_data_scaled * sd(y_data)) + mean(y_data)
    
    # Explanation: (dy/dx)*(1/sd(y)) = slope_scaled -> (dy/dx) = slope_scaled*sd(y)
    # multiply with 10 to calculate change over 10 years
    dat_model$abs_change[i] <- 10*dat_model$slope_lm[i]*dat_model$data_scaled[[i]]$sd_value_per_plot[1]
    
    # p-value is insensitive to scaling
    dat_model$p_value[i] <- summary(dat_model$lm_results[[i]])$coefficients[2,4]
    
    dat_model$R2[i] <- summary(dat_model$lm_results[[i]])$r.squared
    
    dat_model$sens_slope[i] <- dat_model$mk_results[[i]]$B
    
    dat_model$mk_p[i] <- dat_model$mk_results[[i]]$sl
    
  }
  
  # transform p-values of Shapiro-Wilk test to numeric
  dat_model$sw_test_pvalue <- as.numeric(unlist(dat_model$sw_test_pvalue))
  
  dat_model$data <- NULL
  
  return(dat_model)
  
}