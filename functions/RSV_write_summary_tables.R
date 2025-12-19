############################################################################ #
# This file is part of the RSV modelling project.
# 
# => WRITE TABLE WITH SUMMARY STATISTICS
#
#  Copyright 2025, CHERMID, UNIVERSITY OF ANTWERP
############################################################################ #

write_output_tables <- function(sim_output_filename)
{
  
  print(paste('Create tables with all results for', sim_output_filename))
  
  # load data
  sim_output <- readRDS(paste0(sim_output_filename, '.rds'))
  
  # select column names
  colnames_burden_program <- names(sim_output)[grepl('_program$', names(sim_output))]
  colnames_burden         <- names(sim_output)[grepl('_reference$', names(sim_output))]
  
  # realisation id
  num_sim <- unique(sim_output$num_sim)
  sim_output$sim_id <- 1:num_sim

  # convert file names into levels => semi-numeric
  col_file_name <- which(grepl('file', names(sim_output)))
  for(i_col in col_file_name){
    sim_output[, i_col] <- factor(sim_output[, i_col])
  }
  
  if(any(is.na(sim_output$efficacy_maternal))){
    sim_output$efficacy_maternal <- 0
  }
 
  # aggregate: sum 
  system.time(sim_data_global <- aggregate(. ~ config_tag + scenario + perspective +  population_code + intervention + season + sim_id +scenario_id, data = sim_output, sum, na.rm = TRUE, na.action = NULL))
  
  # proportion deaths averted
  sim_data_global$rsv_deaths_averted_proportion <- sim_data_global$rsv_deaths_averted / sim_data_global$rsv_deaths_reference
  sim_data_global$rsv_deaths_averted_proportion[sim_data_global$rsv_deaths_reference == 0] <- 0
  
  # FVP
  sim_data_global$FVP <- (sim_data_global$FVP_maternal + sim_data_global$FVP_mAb)
  
  # Cases averted per 100K fully vaccinated mothers
  sim_data_global$total_cases_averted_per_100K_FVP <-  sim_data_global$total_cases_averted / (sim_data_global$FVP / 100000)
  
  # Deaths averted per 100K fully vaccinated mothers
  sim_data_global$rsv_deaths_averted_per_100K_FVP <-  sim_data_global$rsv_deaths_averted / (sim_data_global$FVP / 100000)
  
  # aggregate 
  # note: included scenario_id to sort the results according the scenarios in the configuration file
  sim_data_mean  <- aggregate(. ~ config_tag + scenario + perspective + population_code + intervention + season + scenario_id, data = sim_data_global, mean, na.rm = TRUE, na.action = NULL) 
  sim_data_CI_LR <- aggregate(. ~ config_tag + scenario + perspective + population_code + intervention + season + scenario_id, data = sim_data_global, quantile, 0.025, na.rm = TRUE, na.action = NULL)
  sim_data_CI_HR <- aggregate(. ~ config_tag + scenario + perspective + population_code + intervention + season + scenario_id, data = sim_data_global, quantile, 0.975, na.rm = TRUE, na.action = NULL)
  
  # FVP per death averted
  sim_data_mean$FVP_per_average_death_averted  <- sim_data_mean$FVP / sim_data_mean$rsv_deaths_averted
  sim_data_CI_LR$FVP_per_average_death_averted <- sim_data_CI_LR$FVP / sim_data_CI_LR$rsv_deaths_averted
  sim_data_CI_HR$FVP_per_average_death_averted <- sim_data_CI_HR$FVP / sim_data_CI_HR$rsv_deaths_averted
  
  # round => mean
  sim_data_mean  <- round_sim_data(sim_data_mean)
  
  # write table with mean and CI
  sim_data_all <- generate_mean_CI_matrix(sim_data_mean,
                                          sim_data_CI_LR,
                                          sim_data_CI_HR,
                                          f_digits = 2)
  
  burden_base_flag        <- which(sim_data_mean$season == "allyear")[1]
  global_data_burden      <- sim_data_mean[burden_base_flag, colnames_burden]
  global_data_burden_ci   <- sim_data_all[burden_base_flag, colnames_burden]
  
  # remove "_reference" tag from column names
  names(global_data_burden)      <- gsub('_reference', '', names(global_data_burden))
  names(global_data_burden_ci)   <- gsub('_reference', '', names(global_data_burden_ci))
  
  intervention_scenario_opt             <- unique(sim_data_mean[,c('intervention', 'scenario', 'perspective', 'population_code')])
  row.names(intervention_scenario_opt)  <- apply(intervention_scenario_opt, 1, paste, collapse = '_')
  
  global_data_mean  <- NULL
  global_data_all   <- NULL 
  global_data_all_k <- NULL
  
  for(i in 1:nrow(intervention_scenario_opt)){
    intervention_scenario_flag <- sim_data_mean$intervention == intervention_scenario_opt$intervention[i] & 
                                  sim_data_mean$scenario == intervention_scenario_opt$scenario[i] &
                                  sim_data_mean$perspective == intervention_scenario_opt$perspective[i] &
                                  sim_data_mean$population_code == intervention_scenario_opt$population_code[i]
    global_data_mean <- rbind(global_data_mean,
                              sim_data_mean[intervention_scenario_flag, colnames_burden_program])
    global_data_all  <- rbind(global_data_all,
                               sim_data_all[intervention_scenario_flag, colnames_burden_program])
  }
  row.names(global_data_mean)  <- apply(intervention_scenario_opt, 1, paste, collapse = '_')
  row.names(global_data_all)   <- apply(intervention_scenario_opt, 1, paste, collapse = '_')
    
 names(global_data_mean) <- gsub('_program', '', colnames_burden_program) # to enable rbind()
 global_data_mean <- rbind(no_intervention = global_data_burden, global_data_mean)
 
 names(global_data_all) <- gsub('_program', '', colnames_burden_program) # to enable rbind()
 global_data_all <- rbind(no_intervention = global_data_burden_ci, global_data_all)
 
 # add column with parameter names
 table_mean  <- data.frame(output = names(global_data_mean), t(global_data_mean), stringsAsFactors = FALSE)
 table_all   <- data.frame(output = names(global_data_all), t(global_data_all), stringsAsFactors = FALSE)
 
 # tables with all and mean values ----
 write.table(table_mean,  file.path(dirname(sim_output_filename), paste0(run_tag, '_model_output_mean.csv')),  sep=',', row.names=F)
 write.table(table_all,   file.path(dirname(sim_output_filename), paste0(run_tag, '_model_output_all.csv')),   sep=',', row.names=F)

 }


################################################################################# #
## SUMMARY TABLES
################################################################################# #
write_summary_tables <- function(sim_output_filename)
{
  
  print(paste('Create summary statistics table for',sim_output_filename))
  
  # load data
  sim_output <- readRDS(paste0(sim_output_filename, '.rds'))
  
  # calculate total averted treatment cost (hospital and non-hospital cases)
  sim_output$treatment_cost_disc_averted   <- sim_output$total_medical_cost_disc_averted 
  
  # express costs in thousands =
  sim_output$treatment_cost_disc_averted_k       <- sim_output$treatment_cost_disc_averted / 1000
  sim_output$intervention_cost_program_k         <- sim_output$intervention_cost_program / 1000
  sim_output$incremental_medical_cost_program_k  <- sim_output$incremental_medical_cost_program / 1000
  sim_output$incremental_cost_disc_program_k     <- sim_output$incremental_cost_disc_program / 1000
  
  # select output columns names
  col_selection <- c('intervention',
                     'scenario',
                     'perspective',
                     'scenario_id',
                     'year',
                     'population_code',
                     'stochastic_proc_id',
                     'total_cases_disc_averted', 
                     'hosp_cases_disc_averted', 
                     'non_ma_cases_disc_averted',
                     'rsv_deaths_disc_averted',
                     'total_QALY_disc_averted',
                     'cost_rsv_non_ma_disc_averted',
                     'cost_rsv_outpatient_disc_averted',
                     'cost_rsv_hosp_disc_averted',
                     'cost_rsv_icu_disc_averted',
                     'treatment_cost_disc_averted',
                     'treatment_cost_disc_averted_k', 
                     'intervention_cost_program',
                     'intervention_cost_program_k',
                     'incremental_medical_cost_program',
                     'incremental_medical_cost_program_k',
                     'total_cost_productivity_loss_disc_averted',
                     'total_cost_productivity_loss_disc_averted',
                     'incremental_cost_disc_program',
                     'incremental_cost_disc_program_k')
  
  # take a subset of the sim_output
  country_data <- sim_output[,col_selection]
  
  # if different years present, aggregate per scenario and country
  if(length(unique(country_data$year)) > 1){
    country_data <- aggregate(. ~ intervention + scenario + perspective + population_code + stochastic_proc_id + scenario_id, data=country_data, sum, na.rm = TRUE, na.action = NULL)
  }
  
  # to compare mAb base case vs. maternal base case
  condition_col = c('total_cases_disc_averted', 
                    'hosp_cases_disc_averted', 
                    'non_ma_cases_disc_averted',
                    'rsv_deaths_disc_averted',
                    'total_QALY_disc_averted',
                    'treatment_cost_disc_averted',
                    'treatment_cost_disc_averted_k', 
                    'intervention_cost_program',
                    'intervention_cost_program_k',
                    'incremental_cost_disc_program',
                    'incremental_cost_disc_program_k')
  
  # remove the column "year" and "stochastic_proc_id"
  country_data$year               <- NULL
  country_data$stochastic_proc_id <- NULL
  
  # calculate mean and confidence intervals over all stochastic realisations
  country_data_mean      <- aggregate(. ~ perspective + population_code + intervention + scenario + scenario_id, data = country_data, mean, na.rm = TRUE, na.action = NULL) 
  country_data_CI_LR     <- aggregate(. ~ perspective + population_code + intervention + scenario + scenario_id, data = country_data, quantile, 0.025, na.rm = TRUE, na.action = NULL)
  country_data_CI_HR     <- aggregate(. ~ perspective + population_code + intervention + scenario + scenario_id, data = country_data, quantile, 0.975, na.rm = TRUE, na.action = NULL)
  
  # round mean output
  country_data_mean  <- round_sim_data(country_data_mean)
  
  # combine mean with confidence intervals (CI)
  country_data_all <- generate_mean_CI_matrix(country_data_mean,
                                              country_data_CI_LR,
                                              country_data_CI_HR)
  # write table with mean values
  write.table(country_data_mean, file.path(dirname(sim_output_filename), paste0(run_tag, '_summary_statistics_mean.csv')),
              sep = ',', row.names = FALSE)
  
  # write table with all statistics
  write.table(country_data_all, file.path(dirname(sim_output_filename), paste0(run_tag, '_summary_statistics_all.csv')),
              sep = ',', row.names = FALSE)
}

## HELP FUNCTION
round_sim_data <- function(sim_data_x,digits_x=2){
  
  # by default, use 2 digits
  flag_numeric              <- unlist(lapply(sim_data_x[1,], is.numeric))
  sim_data_x[,flag_numeric] <- round(sim_data_x[,flag_numeric], digits=2)
  
  # get the columns with absolute values greater or equal to 1 ("ge1")
  flag_less1     <- sim_data_x[,flag_numeric] != 0 & abs(sim_data_x[,flag_numeric]) < 1
  flag_col_less1 <- colSums(flag_less1) == 0 
  col_ge1        <- names(flag_col_less1)[flag_col_less1]
  
  # round with given number of digits for columns with absolute values ≥1
  sim_data_x[,col_ge1] <- round(sim_data_x[,col_ge1], digits = digits_x)
  
  return(sim_data_x)
}

round_sim_data_scale <- function(sim_data_x, digits_x = 2, scale = 1){
  
  if(!is.null(dim(sim_data_x))){
    flag_numeric              <- unlist(lapply(sim_data_x[1,],is.numeric))

    flag_a <- !is.na(sim_data_x[1,]) & flag_numeric & !grepl('_id', names(sim_data_x))
    sim_data_x[,flag_a] <- round(sim_data_x[,flag_a] / scale, digits=digits_x)  
  }
  
  return(sim_data_x)
}

# combine mean with confidence intervals (CI)
generate_mean_CI_matrix <- function(f_mean,f_CI_LR,f_CI_HR,f_digits=2,f_scale = 1){
  
  # round
  f_mean  <- round_sim_data_scale(f_mean, f_digits, f_scale)
  f_CI_LR <- round_sim_data_scale(f_CI_LR, f_digits, f_scale)
  f_CI_HR <- round_sim_data_scale(f_CI_HR, f_digits, f_scale)
  
  # get output columns with numeric content
  flag_numeric     <- unlist(lapply(f_mean[1,], is.numeric)) & (!grepl('_id', names(f_mean)))
  colnames_output  <- colnames(f_mean)[flag_numeric]
  
  # start from the mean...
  f_data_all <- f_mean  
  
  # help function to format the numbers (number of decimals, but no leading spaces)
  specify_decimal <- function(x, k) {
    
    # rounds the values to the specified number of significant digits
    if(all(abs(x[x!=0]) > 10^(k-1),na.rm = TRUE)){
      x <- round(x)
      k <- 0
    } else {
      x <- round(x, k-1)
      k <- k-1
    }
    
    # return
    trimws(format(x, nsmall=k, big.mark=","))
    }
  
  # for each output variable => add median and CI
  for(col_out in colnames_output){
    f_data_all[,col_out] <- specify_decimal(f_mean[,col_out],f_digits)
    
    if(col_out %in% names(f_CI_LR) &&
       col_out %in% names(f_CI_HR) &&
       !all(is.na(c(f_CI_LR[, col_out], f_CI_HR[, col_out])))){
      sel_row <- !(is.na(f_CI_LR[,col_out]) | is.na(f_CI_HR[, col_out]))
      f_data_all[sel_row, col_out] <- paste0(f_data_all[sel_row, col_out],
                                     '  [', specify_decimal(f_CI_LR[sel_row, col_out], f_digits),
                                     ' ; ', specify_decimal(f_CI_HR[sel_row, col_out], f_digits), ']')
    }
  }
  
  # return result
  return(f_data_all)
}


