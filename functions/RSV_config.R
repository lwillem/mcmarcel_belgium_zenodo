############################################################################ #
# This file is part of the McMarcel framework.
# 
# => DEFAULT CONFIGURATION FOR THE RSV ANALYSIS
#
#  Copyright 2025, CHERMID, UNIVERSITY OF ANTWERP
############################################################################ #

# Main function to retrieve model settings
# configList <- sim_config_matrix[i_scen,]
get_rsv_ce_config <- function(configList)
{
  
  # Initialize config variable
  config <- list()
  
  # Add directory names for input and output 
  config$inputFileDir  <- "./input/"
  config$outputFileDir <- configList$outputFileDir

  # time horizon ####
  config$nYearsOfAges  <- 5
  config$monthsInYear  <- 12
  
  # add default config parameters ####
  config$num_sim                     <- configList$num_sim
  config$rng_seed                    <- NA
  config$year                        <- configList$year
  config$population_code             <- configList$population_code
  config$scenario                    <- NA
  config$coverage_maternal           <- NA
  config$coverage_mAb                <- NA
  config$intervention_start          <- NA 
  config$intervention_end            <- NA
  config$coverage_catch_up           <- NA
  config$scenario_id                 <- NA
  
  # QALY parameters for outpatient and hospital cases ####
  # from Zhuxin Mao et al adapted for Belgium
  config$QALYloss_hosp_mean = 0.00982  
  config$QALYloss_hosp_ll   = 0.0085
  config$QALYloss_hosp_ul   = 0.0114
  config$QALYloss_outpatient_mean = 0.00589
  config$QALYloss_outpatient_ll   = 0.0052
  config$QALYloss_outpatient_ul   = 0.0067
  config$QALYloss_hosp_se         = get_se_based_on_ci(config$QALYloss_hosp_ll,config$QALYloss_hosp_ul)
  config$QALYloss_outpatient_se   = get_se_based_on_ci(config$QALYloss_outpatient_ll,config$QALYloss_outpatient_ul)
  
  config$QALY_non_ma_mean = 0.00337 # this is a new value, using the Belgium value set
  config$QALY_non_ma_ll   = 0.0028
  config$QALY_non_ma_ul   = 0.004
  config$QALY_non_ma_se   = get_se_based_on_ci(config$QALY_non_ma_ll,config$QALY_non_ma_ul)
  
  # QALY parameters for caregivers ####
  config$boolean_caregiver_QALYloss   = FALSE    # not included by default
  config$caregiver_QALYloss_hosp_mean = 0.00289
  config$caregiver_QALYloss_hosp_ll   = 0.0016
  config$caregiver_QALYloss_hosp_ul   = 0.0045
  config$caregiver_QALYloss_outpatient_mean = 0.000514
  config$caregiver_QALYloss_outpatient_ll   = 0.0002
  config$caregiver_QALYloss_outpatient_ul   = 0.0009
  config$caregiver_QALY_non_ma_mean = 0 #-0.000284
  config$caregiver_QALY_non_ma_ll   = 0 # -0.0007
  config$caregiver_QALY_non_ma_ul   = 0.0001
  config$caregiver_QALYloss_hosp_se        = get_se_based_on_ci(config$caregiver_QALYloss_hosp_ll,config$caregiver_QALYloss_hosp_ul)
  config$caregiver_QALYloss_outpatient_se  = get_se_based_on_ci(config$caregiver_QALYloss_outpatient_ll,config$caregiver_QALYloss_outpatient_ul)
  config$caregiver_QALY_non_ma_se          = get_se_based_on_ci(config$caregiver_QALY_non_ma_ll,config$caregiver_QALY_non_ma_ul)
  
  # Births ####
  # derive target population, still birth rate and incomplete maternal transfer for 2020-2030, including country name and iso3
  all_country_data <- rbind(read.table('./input/20250106_BE_country_details.csv',sep=',',header=T,stringsAsFactors = F)) 
  
  # check if required country (and year) is present in the data base  
  flag  <- all_country_data$population_code == configList$population_code & all_country_data$year == configList$year 
  if(!any(flag))
  { 
    flag_population         <- all_country_data$population_code == configList$population_code
    country_years           <- all_country_data$year[flag_population]
    country_year_selection  <- country_years[which.min(abs(country_years - as.numeric(configList$year)))]
    flag_year               <- all_country_data$year == country_year_selection
    flag                    <- flag_population & flag_year
    
    if(!exists('warning_population')){
      print(paste('WARNING: NO TARGET POPULATION AND BIRTH RATE DATA AVAILABLE FOR ',configList$population_code,configList$year))
      print(paste('=======> USE DEMOGRAPHIC DATA FROM ',country_year_selection))
      warning_population <<- TRUE
    }
  }
  config$target_population  <- all_country_data$target_population[flag]
  config$pre_term_rate      <- all_country_data$incomplete_maternal_transfer_rate[flag]   
  config$stillbirth_rate    <- all_country_data$stillbirth_rate[flag]
  config$country_iso        <- all_country_data$country_iso[flag]
  
  # Set discounting values 
  config$disc_rate_cost      <- 0.030
  config$disc_rate_effect    <- 0.015 
  
  # RSV hospital incidence ####
  config$filename_hospdata     <- './input/BELX_totalnrhospadmissions_byage_bycalendermonth_perc_S2018_2019_4seasons.csv' 
  config$hospdata_birth_cohort <- 108680 # No adjustment for birth cohort, so assume the 2024 birth cohort 
  config$hospdata_num_seasons  <- 1 
  
  # RSV ICU incidence ####
  config$filename_icudata      <- './input/BELX_totalnrICUadmissions_byage_bycalendermonth_perc_S2018_2019_4seasons.csv' 
  config$icudata_birth_cohort  <- 108680
  config$icudata_num_seasons   <- 1 
  
  # RSV mortality ####
  config$filename_mortality    <- './input/BELX_hCFR_byAge.csv' 
  
  # RSV outpatient visits ####
  config$filename_outpatient_data <- './input/ESPX_outpatient_incidence_1000_infant_months.csv'# hosp outpatient incidence per 1000 infant-months, by age (rows) for different samples (columns)
  
  # Wheezing and asthma ####
  config$wheezing_horizon   <- 0 # if 0, do not account for recurrent wheeze and asthma due to hospital admission in the first year of life

  config$wheezing_under1_shape1 <- 13
  config$wheezing_under1_shape2 <- 42-13
  config$wheezing_2_shape1 <- 135
  config$wheezing_2_shape2 <- 504-135
  config$wheezing_3_shape1 <- 87 
  config$wheezing_3_shape2 <- 504-87  
  config$wheezing_4_shape1 <- 80 
  config$wheezing_4_shape2 <- 504-80  
  config$wheezing_5_shape1 <- 50  
  config$wheezing_5_shape2 <- 504-50 
  
  # QALY loss due to recurrent wheezing and asthma per year for children
  config$QALY_wheezing_mean <- 0.0183   
  config$QALY_wheezing_se   <- 0.00001

  # Non-medically attended cases ####
  config$filename_proportion_non_ma <- './input/RESCEUX_nonMA_cases_0-59m.csv' 
  
  # Cost parameters ####
  # includes hospital, ICU, nonMA, wheezing and asthma
  config$filename_cost <- './input/BELX_cost_perspectives_regions.csv'
  config$perspective   <- 'payer' # health care payer
  
  # administration cost
  config$admin_cost_maternal     <- 0 # switched off for now 
  config$admin_cost_mAb          <- 0 # switched off for now
  config$price_dose_maternal_min <- 200 # if equal to max, always the same price     
  config$price_dose_maternal_max <- 200 # if equal to min, always the same price   
  config$price_dose_mAb_min      <- 200 # if equal to max, always the same price 
  config$price_dose_mAb_max      <- 200 # if equal to min, always the same price 
  
  # dose wastage
  config$wastage_maternal <- 0 
  config$wastage_mAb      <- 0 
  
  # yearly fixed implementation costs
  config$implementation_cost_mAb      <- 0 # switched off for now
  config$implementation_cost_maternal <- 0 # switched off for now

  # productivity loss estimates for societal perspective 
  config$productivity_cost_outpatient_3to11m_mean  <- 636.5
  config$productivity_cost_outpatient_3to11m_se    <- 134.9
  config$productivity_cost_outpatient_12to47m_mean <- 1296.2
  config$productivity_cost_outpatient_12to47m_se   <- 281.0
  config$productivity_cost_inpatient_3to11m_mean   <- 1318.1
  config$productivity_cost_inpatient_3to11m_se     <- 207.4
  config$productivity_cost_inpatient_12to47m_mean  <- 2107.0
  config$productivity_cost_inpatient_12to47m_se    <- 930.1
   
  # Efficacy parameters ####
  # differential efficacy for outpatient, hospitalization and mortality
  config$efficacy_maternal_nonMA_mean       <- NA # NA means to use file specified in csv
  config$efficacy_maternal_nonMA_stdev      <- NA 
  config$efficacy_maternal_outpatient_mean  <- NA 
  config$efficacy_maternal_outpatient_stdev <- NA
  config$efficacy_maternal_hosp_mean        <- NA
  config$efficacy_maternal_hosp_stdev       <- NA
  config$efficacy_maternal_icu_mean         <- NA
  config$efficacy_maternal_icu_stdev        <- NA
  config$efficacy_maternal_mortality_mean   <- NA
  config$efficacy_maternal_mortality_stdev  <- NA
  config$dur_protection_maternal            <- NA
  
  config$efficacy_mAb_nonMA_mean       <- NA 
  config$efficacy_mAb_nonMA_stdev      <- NA
  config$efficacy_mAb_outpatient_mean  <- NA 
  config$efficacy_mAb_outpatient_stdev <- NA
  config$efficacy_mAb_hosp_mean       <- NA
  config$efficacy_mAb_hosp_stdev      <- NA
  config$efficacy_mAb_icu_mean        <- NA
  config$efficacy_mAb_icu_stdev       <- NA
  config$efficacy_mAb_mortality_mean  <- NA
  config$efficacy_mAb_mortality_stdev <- NA
  config$dur_protection_mAb           <- NA 
  
  # external efficacy file names (in the 'input' directory)
  config$efficacy_mAb_nonMA_file      <- NA
  config$efficacy_mAb_outpatient_file <- NA
  config$efficacy_mAb_hosp_file       <- NA
  config$efficacy_mAb_icu_file        <- NA 
  config$efficacy_mAb_mortality_file  <- NA 
  config$efficacy_maternal_nonMA_file <- NA
  config$efficacy_maternal_outpatient_file <- NA
  config$efficacy_maternal_hosp_file  <- NA
  config$efficacy_maternal_icu_file   <- NA 
  config$efficacy_maternal_mortality_file  <- NA 
  
  # option to save the efficacy values to file
  config$bool_save_protection_files <- NA
  
  # Transfer provided configuration details ####
  ##################################################### #
  # Check if the provided configuration, has only valid parameters
  # valid if all column names are present in the default configuration
  if(any(!names(configList) %in% names(config))){
    stop(paste(c('!! STOP EXECUTION => UNUSED CONFIG PARAMETER(S): ',
                paste(names(configList)[!names(configList) %in% names(config)],collapse=', '))))
  }
  
  # transfer the function parameter to 'config'
  config[names(configList)] <- configList
  
  # check file names and complete with input dir if needed
  sel_filename          <- names(config)[grepl('filename', names(config))]
  sel_filename_complete <- dirname(as.character(config[sel_filename])) != './input'
  sel_filename          <- sel_filename[sel_filename_complete]
  config[sel_filename]  <- paste0('./input/', config[sel_filename], '.csv')
  
  # Select cost data ####
  ################################# #
  # national or regional? && perspective?
  db_cost         <- read.csv(config$filename_cost,header=T)
  opt_perspective <- c(db_cost$perspective,'societal')
  
  if(!(config$perspective %in% opt_perspective)){
    stop(paste('Perspective in model config does not match the cost data:',paste(opt_perspective,collapse=', ')))
  }
  
  # Account for the perspective
  # if societal: use payer perspective, and use productivity loss information
  # if not: use given perspective and do not account for productivity loss
  if(config$perspective == 'societal'){
    db_cost <- db_cost[db_cost$perspective == 'payer',]
    config$bool_productivity_cost <- TRUE
  } else{
    db_cost <- db_cost[db_cost$perspective == config$perspective,]
    config$bool_productivity_cost <- FALSE
  }
  
  # identify population/regional code
  opt_pop_codes  <- unlist(strsplit(config$population_code, '_'))
  population_set <- opt_pop_codes[length(opt_pop_codes)]
  
  # select national/regional costs
  sel_col <- grepl(population_set,names(db_cost)) | # to capture regional values
              grepl('hosp', names(db_cost)) |  # include national hospital estimates
              grepl('icu', names(db_cost)) |   # include national icu estimates
              grepl('non_ma', names(db_cost))  # include national non-MA estimates
  db_cost <- db_cost[,sel_col]
  
  # remove national/regional tag
  names(db_cost) <- gsub('_...$', '', names(db_cost))
  
  # estimate sd for outpatient costs 
  db_cost$cost_outpatient_under1y_sd = get_se_based_on_ci(db_cost$cost_outpatient_under1y_uci, db_cost$cost_outpatient_under1y_lci)
  db_cost$cost_outpatient_1y_sd = get_se_based_on_ci(db_cost$cost_outpatient_1y_uci, db_cost$cost_outpatient_1y_lci)
  db_cost$cost_outpatient_2y_sd = get_se_based_on_ci(db_cost$cost_outpatient_2y_uci, db_cost$cost_outpatient_2y_lci)
  db_cost$cost_outpatient_3y_sd = get_se_based_on_ci(db_cost$cost_outpatient_3y_uci, db_cost$cost_outpatient_3y_lci)
  db_cost$cost_outpatient_4y_sd = get_se_based_on_ci(db_cost$cost_outpatient_4y_uci, db_cost$cost_outpatient_4y_lci)
  
  # add cost db to config list
  config$db_cost <- db_cost
  
  # cost for recurrent wheezing and asthma, per year
  config$wheezing_cost_peryear <- config$db_cost$cost_wheezing_asthma
  
  # cost for non-medically attended case
  config$non_ma_cost <- config$db_cost$cost_non_ma
  
  # Consistency checks ####
  ################################# #
  
  flag <- (grepl('efficacy', names(config)) |
           grepl('prob', names(config)) |
           grepl('rate', names(config)) |
           grepl('coverage', names(config))) &
           !grepl('file', names(config))
  
  if(any(unlist(config[flag])>1,na.rm=T)){
    smd_print('INCONSISTENT CONFIGURATION: EFFICACY, RATE, PROBABILITY OR COVERAGE > 1 WITH', configList$config_tag,WARNING = TRUE, FORCED = TRUE)
  }
  
  # check whether the requested efficacy files exist
  file_names_efficacy <- unlist(config[grepl('efficacy_.*_file', names(config))])
  file_names_efficacy <- file_names_efficacy[!is.na(file_names_efficacy)]
  if(length(file_names_efficacy) > 0){
    # update file names with full path
    file_names_efficacy_full <- paste0('input/', file_names_efficacy, '.csv')
    config[names(file_names_efficacy)] <- file_names_efficacy_full
    # check availability
    efficacy_file_present <- file.exists(file_names_efficacy_full)
    if((any(!efficacy_file_present))){
      print(paste('REQUESTED EFFICACY FILE(S) NOT AVAILABLE IN THE "input" DIRECTORY: ', paste(file_names_efficacy[!efficacy_file_present], collapse = ', ')))
    }
  }

  # if given, option to include full path for the non-medically attendance file
  if(!is.na(config$filename_proportion_non_ma) & !grepl('.csv', config$filename_proportion_non_ma)){
    config$filename_proportion_non_ma <- paste0('input/', config$filename_proportion_non_ma,'.csv')
  }
  
  # check whether intervention start and end are valid months, i.e. part of [1:12]
  intervention_months       <- unlist(config[c('intervention_start','intervention_end')])
  intervention_months_valid <- intervention_months %in% 1:12
  if(any(!intervention_months_valid)){
    print(paste('INCONSISTENT INTERVENTION TIMING(S): ', paste(names(intervention_months)[!intervention_months_valid], collapse = ', ')))
  }
  
  # check whether given booleans, are 0 or 1
  boolean_values <- unlist(config[grepl('boolean_.*',names(config))])
  boolean_values_valid <- (boolean_values %in% 0:1)
  if(any(!boolean_values_valid)){
    print(paste('INCONSISTENT BOOLEAN VALUE(S): ', paste(names(boolean_values)[!boolean_values_valid],collapse = ', ')))
  } else{
    config[names(boolean_values)] <- config[names(boolean_values)] == 1
  }

  # Include intervention details and labels ----
  config$intervention_months    <- elapsed_months(config$intervention_start, config$intervention_end) + 1
  config$proportion_target_pop  <- ifelse(config$coverage_catch_up > 0,
                                          1,
                                          config$intervention_months / 12)
  if(config$coverage_mAb > 0){
    config$intervention <- paste(c(config$intervention, 'mAb'), collapse = '_')
  }
  if(config$coverage_maternal > 0){
    config$intervention <- paste(c(config$intervention, 'maternal'), collapse = '_')
  }
  
  if(config$intervention_months == 12){
    config$season <- 'allyear'
  } else{
    config$season <- paste0(month_id2name(config$intervention_start),
                            month_id2name(config$intervention_end),
                            ifelse(config$coverage_catch_up > 0,
                                   '_CB',
                                   ''))
    config$intervention <- paste(config$intervention, config$season, sep = '_')
  }
  
  # Include config tag ----
  config$config_tag <- paste(config$intervention, config$scenario, sep = '_')
  
  # return config list
  return(config)
}

## HELP FUNCTIONS 

# get standard error based on confidence interval
get_se_based_on_ci <- function(ll, ul){
  return((ul - ll) / (2 * 1.96))
  }

# sample from Beta distribution
sample_beta_dist <- function(num_samples,distr_mean,distr_sd){
  # sample from Beta distribution
    beta.a     <- (distr_mean^2 * (1 - distr_mean) / distr_sd^2) - distr_mean
    beta.b     <- beta.a * (1 - distr_mean) / distr_mean
    sample_all <- rbeta(num_samples,shape1=beta.a,shape2=beta.b)
    
  # truncate negative values
  sample_all[sample_all < 0] <- 0 
  
  #return samples
  return(sample_all)
}
