############################################################################ #
# This file is part of the RSV modelling project.
# 
# => GET EFFECTIVE PROTECTION OVER TIME 
# 
#  Copyright 2025, CHERMID, UNIVERSITY OF ANTWERP
############################################################################ #


# Get protection over time given the provided parameters
# note: currently, sampled protection from a distribution is constant over time
get_protection <-function(config,
                          intervention){
  
  # if both efficacy files and no mean and stdev values are given, STOP
  file_names_efficacy    <- unlist(config[grepl('efficacy_.*_file', names(config)) & grepl(intervention, names(config))])
  sample_values_efficacy <- unlist(config[(grepl('efficacy_.*_mean', names(config)) |
                                           grepl('efficacy_.*_stdev', names(config)) |
                                           grepl('dur_protection', names(config))) & 
                                           grepl(intervention, names(config)) ])
  

  # if both efficacy files and mean and stdev values are given, STOP
  file_names_efficacy    <- file_names_efficacy[!is.na(file_names_efficacy)]
  sample_values_efficacy <- sample_values_efficacy[!is.na(sample_values_efficacy)]
  if(length(file_names_efficacy) > 0 & length(sample_values_efficacy) > 0){
    stop('Conflicting information for the efficacy, both external files as mean and stdev values are provided. The model cannot choose, so we stop the calculation.')
  }
  
  ######################## DISTRIBUTION ######################################################################################### #
  # if mean and stdev values are given and NO file names: sample from distribution
  if(length(sample_values_efficacy) >= 5){
    
    # remove intervention label
    names(sample_values_efficacy) <- names(sample_values_efficacy) <- gsub(paste0('_',intervention),'', names(sample_values_efficacy))
    
    # if ICU VE is missing, use hospital VE
    if(!any(grepl('icu', names(file_names_efficacy)))){
      sample_values_efficacy['efficacy_icu_mean']  <- sample_values_efficacy['efficacy_hosp_mean']
      sample_values_efficacy['efficacy_icu_stdev'] <- sample_values_efficacy['efficacy_hosp_stdev']
    }
    
    # if mortality VE is missing, use ICU VE
    if(!any(grepl('mortality', names(file_names_efficacy)))){
      sample_values_efficacy['efficacy_mortality_mean']  <- sample_values_efficacy['efficacy_icu_mean'] 
      sample_values_efficacy['efficacy_mortality_stdev'] <- sample_values_efficacy['efficacy_icu_stdev'] 
    }
    
    # if nonMA VE is missing, use outpatient VE
    if(!any(grepl('nonMA', names(file_names_efficacy)))){
      sample_values_efficacy['efficacy_nonMA_mean']  <- sample_values_efficacy['efficacy_outpatient_mean']
      sample_values_efficacy['efficacy_nonMA_stdev'] <- sample_values_efficacy['efficacy_outpatient_stdev']
    }
    
    protection_data <- sample_effective_protection_constant(config=config,
                                                              dur_protection       = sample_values_efficacy['dur_protection'],
                                                              efficacy_nonMA_mean  = sample_values_efficacy['efficacy_nonMA_mean'],
                                                              efficacy_nonMA_stdev = sample_values_efficacy['efficacy_nonMA_stdev'],
                                                              efficacy_outpatient_mean  = sample_values_efficacy['efficacy_outpatient_mean'],
                                                              efficacy_outpatient_stdev = sample_values_efficacy['efficacy_outpatient_stdev'],
                                                              efficacy_hosp_mean  = sample_values_efficacy['efficacy_hosp_mean'],
                                                              efficacy_hosp_stdev = sample_values_efficacy['efficacy_hosp_stdev'],
                                                              efficacy_icu_mean   = sample_values_efficacy['efficacy_icu_mean'],
                                                              efficacy_icu_stdev  = sample_values_efficacy['efficacy_icu_stdev'])
    
    # assume protection against mortality equal to protection against ICU admission
    protection_data$mortality <- protection_data$icu 
    
  } else if(length(sample_values_efficacy) > 0){
    stop('Insufficient data provided on the protection of ',intervention,': the model needs at least mean efficacy and stdev for the outpatient and hospital cases + duration or protection. The model cannot continue, so we stop the calculation.')
  }
  
  ######################## FILE ######################################################################################### #
  # if efficacy files are given and NO mean and stdev values: use the data in the files
   if(length(file_names_efficacy) > 1){

     # remove intervention label
     names(file_names_efficacy) <- gsub(paste0('_',intervention),'', names(file_names_efficacy))
         
     # if ICU VE is missing, use hospital VE
     if(!any(grepl('icu', names(file_names_efficacy)))){
       file_names_efficacy['efficacy_icu_file'] <- file_names_efficacy['efficacy_hosp_file'] 
     }
     
     # if mortality VE is missing, use ICU VE
     if(!any(grepl('mortality', names(file_names_efficacy)))){
       file_names_efficacy['efficacy_mortality_file'] <- file_names_efficacy['efficacy_icu_file']
     }
     
     # if nonMA VE is missing, use outpatient VE
     if(!any(grepl('nonMA', names(file_names_efficacy)))){
       file_names_efficacy['efficacy_nonMA_file'] <- file_names_efficacy['efficacy_outpatient_file']
     }
     
    # initiate list
    protection_data <- list()
    
    # specify disease outcomes
    db_severity <- data.frame(category = gsub('efficacy_','', gsub('_file','', names(file_names_efficacy))),
                              file = file_names_efficacy)
    
    # fix for hosp
    db_severity$category <- gsub('hosp$', 'hospital', db_severity$category)
    
    # load values
    i_severity <- 1
    for(i_severity in 1:nrow(db_severity)){
      
      age_protection    <- load_efficacy_from_file(config = config, db_severity$file[i_severity])
      age_protection_0m <- age_protection[1,]
      protection_data[[db_severity$category[i_severity]]] <- age_protection
      
      # option to save the efficacy values to file
      if(!is.na(config$bool_save_protection_files) & config$bool_save_protection_files){
        write.table(t(protection_data$outpatient),
                    file = paste0(config$outputFileDir, '/VE_', config$config_tag, '_outpatient.csv'), sep=',', row.names = FALSE, col.names = FALSE)
      }
    }
      
   } else if(length(file_names_efficacy) > 0){
     stop('Insufficient data provided on the protection of ',intervention,': the model needs three efficacy files. The model cannot continue, so we stop the calculation.')
  }
  
  return(protection_data)
}

# sample effective protection, constant over time, given the provided mean and stdev
load_efficacy_from_file <-function(config,
                                   efficacy_file){
  
  # make sure the ages are in the rows, and sample in the columns
  age_protection  <- t(read.table(efficacy_file,sep=',',header = FALSE))

  # account for the number of simulations
  if(config$num_sim <= ncol(age_protection)){
    age_protection  <- age_protection[,1:config$num_sim]
  } else{
    print('WARNING: more efficacy values required than provided, (re)sample from provided efficacy estimates')
    sample_col      <- sample(ncol(age_protection),config$num_sim,replace = TRUE)
    age_protection  <- age_protection[,sample_col]
  }
  
  # make sure you have config$nMonthsOfAges rows
  if(nrow(age_protection) < config$nMonthsOfAges){
    age_protection <- rbind(age_protection,
                           matrix(data = 0,
                                  nrow = config$nMonthsOfAges-nrow(age_protection),
                                  ncol = ncol(age_protection)))
  }
  
  # return values
  return(age_protection)
}

# sample effective protection, constant over time, given the provided mean and stdev
sample_effective_protection_constant <-function(config,
                                                dur_protection,
                                                efficacy_nonMA_mean, efficacy_nonMA_stdev,
                                                efficacy_outpatient_mean, efficacy_outpatient_stdev,
                                                efficacy_hosp_mean, efficacy_hosp_stdev,
                                                efficacy_icu_mean, efficacy_icu_stdev){
  
  # convert duration of protection from years to months
  dur_prot <- round(config$monthsInYear * dur_protection) 
  
  # Calculate the effective protection by outcome 
  # init matrices (i.e. age (in months)*sim empty matrices)
  ageEffectiveProtection_outpatient   <- matrix(0, ncol=config$num_sim, nrow=config$nMonthsOfAges)
  ageEffectiveProtection_nonMA        <- ageEffectiveProtection_outpatient
  ageEffectiveProtection_hospital     <- ageEffectiveProtection_outpatient
  ageEffectiveProtection_icu          <- ageEffectiveProtection_outpatient
  
  # Option to use differential efficacy for outpatient and hospitalisation
  efficacy_outpatient   <- sample_beta_dist(config$num_sim, efficacy_outpatient_mean,     
                                                    efficacy_outpatient_stdev)
  # re-use the all-cause efficacy if the provided mean and stdev are equal
  if(efficacy_outpatient_mean == efficacy_hosp_mean &&
     efficacy_outpatient_stdev == efficacy_hosp_stdev){
     efficacy_hospital <- efficacy_outpatient
  } else{ # or sample new efficacy values
    efficacy_hospital  <- sample_beta_dist(config$num_sim, efficacy_hosp_mean, efficacy_hosp_stdev)
  }
  
  # nonMA: sample or re-use the outpatient efficacy if the provided mean and stdev are equal
  if( efficacy_nonMA_mean == efficacy_outpatient_mean  &&
      efficacy_nonMA_stdev == efficacy_outpatient_stdev){
    efficacy_nonMA <- efficacy_outpatient
  } else {
    efficacy_nonMA  <- sample_beta_dist(config$num_sim, efficacy_nonMA_mean,     
                                        efficacy_nonMA_stdev)
  }
  
  # ICU efficacy
  # re-use the hospital efficacy if the provided mean and stdev are equal
  if(efficacy_icu_mean == efficacy_hosp_mean &&
     efficacy_icu_stdev == efficacy_hosp_stdev){
    efficacy_icu <- efficacy_hospital
  } else{ # or sample new efficacy values
    efficacy_icu  <- sample_beta_dist(config$num_sim, efficacy_icu_mean,
                                           efficacy_icu_stdev)
  }

  # Fill matrices: without uncertainty for duration of protection
  i_eff <- 1
  for(i_eff in 1:config$num_sim){
      ageEffectiveProtection_nonMA[1:dur_prot,i_eff]      <- efficacy_nonMA[i_eff] 
      ageEffectiveProtection_outpatient[1:dur_prot,i_eff] <- efficacy_outpatient[i_eff] 
      ageEffectiveProtection_hospital[1:dur_prot,i_eff]   <- efficacy_hospital[i_eff]
      ageEffectiveProtection_icu[1:dur_prot,i_eff]        <- efficacy_icu[i_eff]
  }
  
  # return values
  return(list(nonMA      = ageEffectiveProtection_nonMA,
              outpatient = ageEffectiveProtection_outpatient,
              hospital   = ageEffectiveProtection_hospital,
              icu        = ageEffectiveProtection_icu))
  
}
