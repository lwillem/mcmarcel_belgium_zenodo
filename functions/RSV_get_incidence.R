############################################################################ #
# This file is part of the RSV modelling project.
# 
# => GET INCIDENCE DATA 
# 
#  Copyright 2025, CHERMID, UNIVERSITY OF ANTWERP
############################################################################ #

# Get disease incidence (prevented) over time given the provided config$season and effective protection
get_incidence_data <-function(config){
  
  # RSV hospitalisations ####
  RSV_hospdata <- read.table(config$filename_hospdata, sep = ',', header = T, stringsAsFactors = F) # contains total nr of hospitalisations (sum over X seasons) by age (in months, up to age 5, rows) and by calendar month (columns)
  RSV_hosp_age_calendarmonth_average <- RSV_hospdata[,-1] / config$hospdata_num_seasons # 60*12 matrix; remove first column (containing age groups)
  RSV_hosp_age_calendarmonth         <- as.matrix(RSV_hosp_age_calendarmonth_average) / config$hospdata_birth_cohort * config$target_population
  
  # get age-specific distribution of hospital admissions by calendar month
  RSV_prop_hosp_age_calendarmonth <- multiply_matrix_by_age(RSV_hosp_age_calendarmonth,1/rowSums(RSV_hosp_age_calendarmonth))
  
  # RSV ICU ####
  RSV_icudata <- read.table(config$filename_icudata, sep = ',', header = T, stringsAsFactors = F) # contains total nr of ICU admissions (sum over X seasons) by age (in months, up to age 5, rows) and by calendar month (columns)
  RSV_icu_age_calendarmonth_average <- RSV_icudata[,-1] / config$icudata_num_seasons # 60*12 matrix; remove first column (containing age groups)
  RSV_icu_age_calendarmonth         <- as.matrix(RSV_icu_age_calendarmonth_average) / config$icudata_birth_cohort * config$target_population

  # RSV outpatient cases ####
  RSV_outpatient <- read.table(config$filename_outpatient_data, sep = ',', header = T, stringsAsFactors = F)
  RSV_outpatient_age_sim <- as.matrix(RSV_outpatient) # contains hosp outpatient incidence by age (in months, up to age 5, rows) and by calendar month (columns)
  RSV_outpatient_age_sim <- RSV_outpatient_age_sim[,-1] # remove age column
  RSV_outpatient_age_sim <- RSV_outpatient_age_sim / 1000 * config$target_population

  # apply the distribution by calendar month from the hospital incidence to the outpatient incidence
  RSV_outpatient_age_calendarmonth_sim <- array(data = NA, dim = c(nrow(RSV_prop_hosp_age_calendarmonth), ncol(RSV_prop_hosp_age_calendarmonth), config$num_sim))
  RSV_outpatient_age_cohort_sim        <- RSV_outpatient_age_calendarmonth_sim
  i_sim <- 1
  for(i_sim in 1:config$num_sim){
    # sample random row index if i_sim > number of samples
    if(i_sim > ncol(RSV_outpatient_age_sim)){ 
      i_sim <- sample(ncol(RSV_outpatient_age_sim),1) 
    }
    RSV_outpatient_age_calendarmonth_sim[,,i_sim] <- RSV_prop_hosp_age_calendarmonth * RSV_outpatient_age_sim[,rep(i_sim,12)]
    RSV_outpatient_age_cohort_sim[,,i_sim]        <- time_age_to_cohort_age(RSV_outpatient_age_calendarmonth_sim[,,i_sim])
  }
  
  # RSV deaths #### 
  RSV_hCFR_age   <- data.frame(read.table(config$filename_mortality, sep = ',', header = T, stringsAsFactors = F)) # contains hosp outpatient incidence by age (in months, up to age 5, rows) and by calendar month (columns)
  RSV_hCFR_age   <- RSV_hCFR_age[,-1] # remove age column

  # make sure that hCFR is strictly positive
  RSV_hCFR_age[RSV_hCFR_age < 0] <- 0
  
  RSV_deaths_hosp_age <- multiply_matrix_by_age(RSV_hosp_age_calendarmonth, RSV_hCFR_age$hCFR_hosp)
  RSV_deaths_icu_age  <- multiply_matrix_by_age(RSV_icu_age_calendarmonth, RSV_hCFR_age$hCFR_icu)
  RSV_deaths_age      <- RSV_deaths_hosp_age + RSV_deaths_icu_age
  
  # Wheezing ####
  # Probability recurrent wheezing and asthma after RSV hospital admission in the first year of life
  # for these wheezing/asthma cases, consequences (QALY loss and costs) in the first 3 years of life (wheezing) or up to age 13 (asthma) are accounted for
  prob_wheezing_under1 = rbeta(config$num_sim, config$wheezing_under1_shape1, config$wheezing_under1_shape2) # probability recurrent wheezing in <1 year olds
  
  # probability recurrent wheezing in the second, third, fourth and fifth year of life, if hospitalised for RSV in the first year of life
  # these probabilities are used to calculate total costs and QALY loss due to recurrent wheezing and asthma,
  # i.e. with age, fewer children (who initially developed recurrent wheeze/asthma) keep suffering from recurrent wheezing and asthma
  prob_wheezing_2 <- rbeta(config$num_sim, config$wheezing_2_shape1, config$wheezing_2_shape2)
  prob_wheezing_3 <- rbeta(config$num_sim, config$wheezing_3_shape1, config$wheezing_3_shape2)
  prob_wheezing_4 <- rbeta(config$num_sim, config$wheezing_4_shape1, config$wheezing_4_shape2)
  prob_wheezing_5 <- rbeta(config$num_sim, config$wheezing_5_shape1, config$wheezing_5_shape2)
  
  # make sure "config$wheezing_horizon" is an integer that can be used as index
  wheezing_horizon <- floor(config$wheezing_horizon)
  
  # 12 months x number of simulations
  hosp_cases_0to1y <- matrix(1, nrow=12, ncol=config$num_sim)
  wheezing_incidence_time  <- abind(hosp_cases_0to1y * rep(prob_wheezing_under1, each=12),
                                hosp_cases_0to1y * rep(prob_wheezing_2, each=12),
                                hosp_cases_0to1y * rep(prob_wheezing_3, each=12),
                                hosp_cases_0to1y * rep(prob_wheezing_4, each=12),
                                hosp_cases_0to1y * rep(prob_wheezing_5, each=12), 
                                along=3) 
  
  # extend the age dimension if more years are needed
  wheezing_age_included <- dim(wheezing_incidence_time)[3]
  if(wheezing_horizon > wheezing_age_included){
    sel_age_index <- c(1:wheezing_age_included, rep(wheezing_age_included, wheezing_horizon - wheezing_age_included)) 
    wheezing_incidence_time <- wheezing_incidence_time[,,sel_age_index]
  }
  # set all counts for (year > wheezing_horizon) to 0
  sel_age_index <- which(!seq(dim(wheezing_incidence_time)[3]) %in%  (0:wheezing_horizon))
  wheezing_incidence_time[,,sel_age_index] <- wheezing_incidence_time[,,sel_age_index] * 0
  
  # # Non-medically attended cases ####
  # load proportion non-MA / MA
  if(!is.na(config$filename_proportion_non_ma)){ 
    proportion_non_ma=read.table(config$filename_proportion_non_ma, sep = ',', header = T, stringsAsFactors = F)
    # adjust to "time*sim" and if needed, extend up to 5 years of age
    proportion_non_ma <- t(proportion_non_ma)
    if(nrow(proportion_non_ma)==12){
      proportion_non_ma <- rbind(proportion_non_ma,matrix(0, nrow=60-12, ncol = ncol(proportion_non_ma)))
    }
  } else{
    proportion_non_ma <- matrix(0, nrow = nrow(RSV_hosp_age_calendarmonth), ncol = config$num_sim)
  }

  if(ncol(proportion_non_ma)>=config$num_sim){
    proportion_non_ma <- proportion_non_ma[,1:config$num_sim]
  } else{
    proportion_non_ma <- proportion_non_ma[, sample(ncol(proportion_non_ma), config$num_sim,replace = T)]
  }

  # translate proportion into incidence based on the medically attended cases
  # proportion = nonMA / MA
  RSV_non_ma_age_cohort_sim <- array(data=NA, dim=c(nrow(RSV_hosp_age_calendarmonth), ncol(RSV_hosp_age_calendarmonth), config$num_sim))
  i_sim <- 1
  for(i_sim in 1:config$num_sim){
    RSV_ma_age_calendarmonth           <- RSV_hosp_age_calendarmonth + RSV_outpatient_age_calendarmonth_sim[,,i_sim] 
    RSV_non_ma_age_calendarmonth       <- multiply_matrix_by_age(RSV_ma_age_calendarmonth, proportion_non_ma[,i_sim])
    RSV_non_ma_age_cohort_sim[,,i_sim] <- time_age_to_cohort_age(RSV_non_ma_age_calendarmonth)
  }
  
  # define summary metrics for EVPPI
  proportion_non_ma_mean    = colMeans(proportion_non_ma) #ideally, should be proportion_non_ma[1,], ie uncertainty of one age group but degree of uncertainty (slope of fitted regression) is assumed to be fixed

  # Burden of disease (by cohort) ####
  RSV_hosp_age_cohort           <- time_age_to_cohort_age(RSV_hosp_age_calendarmonth)
  RSV_icu_age_cohort            <- time_age_to_cohort_age(RSV_icu_age_calendarmonth)
  RSV_deaths_age_cohort         <- time_age_to_cohort_age(RSV_deaths_age)
  
  # aggregate by age
  RSV_hosp_age           <- rowSums(RSV_hosp_age_cohort)
  RSV_icu_age            <- rowSums(RSV_icu_age_cohort)
  RSV_outpatient_age     <- apply(RSV_outpatient_age_cohort_sim, 3, rowSums)
  RSV_non_ma_age         <- apply(RSV_non_ma_age_cohort_sim, 3, rowSums)
  RSV_deaths_age         <- rowSums(RSV_deaths_age_cohort)
  
  # duplicate to get one column per simulation
  RSV_hosp_age    <- matrix(nrow=length(RSV_hosp_age), ncol = config$num_sim, RSV_hosp_age) # allow for uncertainty
  RSV_icu_age     <- matrix(nrow=length(RSV_icu_age), ncol = config$num_sim, RSV_icu_age) # allow for uncertainty
  RSV_deaths_age  <- matrix(nrow=length(RSV_deaths_age), ncol = config$num_sim, RSV_deaths_age) # allow for uncertainty
  
  # Wheezing ----
  wheezing_cases_time <- wheezing_incidence_time * 0
  hosp_cases_0to1y    <- RSV_hosp_age[1:12,]
  
  for(i_year in 1:dim(wheezing_cases_time)[3]){
    wheezing_cases_time[,,i_year]         <- hosp_cases_0to1y * wheezing_incidence_time[,,i_year]
  }
  
  # hospital admissions after first year of age do not account for wheezing or asthma
  wheezing_cases_complement <- matrix(0, nrow = 60-12, ncol = config$num_sim)
  wheezing_cases            <- rbind(wheezing_cases_time[,,1], wheezing_cases_complement)
  
  # fix for when number of simulation is one
  if(config$num_sim == 1)
  {
    RSV_outpatient_age  <- as.matrix(RSV_outpatient_age, ncol=1)
    RSV_hosp_age        <- as.matrix(RSV_hosp_age, ncol=1)
    RSV_icu_age         <- as.matrix(RSV_icu_age, ncol=1)
    RSV_deaths_age      <- as.matrix(RSV_deaths_age, ncol=1)
    RSV_non_ma_age      <- as.matrix(RSV_non_ma_age, ncol=1) 
    wheezing_cases      <- as.matrix(wheezing_cases, ncol=1) 
  }
  
  # total cases
  RSV_total_cases <- RSV_hosp_age + RSV_icu_age + RSV_non_ma_age + RSV_outpatient_age  
  
  # combine incidence rates and probabilities
  df_country  <- list(hosp_cases       = RSV_hosp_age,
                      icu_cases        = RSV_icu_age,
                      outpatient_cases = RSV_outpatient_age,
                      rsv_deaths       = RSV_deaths_age,
                      non_ma_cases     = RSV_non_ma_age,
                      wheezing_cases   = wheezing_cases,
                      wheezing_cases_time = wheezing_cases_time,
                      total_cases      = RSV_total_cases,
                      
                      outpatient_agedist_under6months  = RSV_outpatient_age[1:6,], #matrix(0, 6, config$num_sim),  
                      prob_wheezing_under1      = prob_wheezing_under1,
                      proportion_non_ma_mean    = proportion_non_ma_mean,
                      
                      # export age-cohort-sim specific arrays
                      RSV_hosp_age_cohort           = RSV_hosp_age_cohort, 
                      RSV_icu_age_cohort            = RSV_icu_age_cohort,
                      RSV_outpatient_age_cohort_sim = RSV_outpatient_age_cohort_sim, 
                      RSV_non_ma_age_cohort_sim     = RSV_non_ma_age_cohort_sim, 
                      RSV_deaths_age_cohort         = RSV_deaths_age_cohort,
                      wheezing_incidence_time       = wheezing_incidence_time
                      )
                      
  # load data and return
  return(df_country) 

}


# Get disease incidence (prevented) over time given the provided config$season and effective protection
get_prevented_incidence <-function(config,
                                   df_country,
                                   ageEfficacy,
                                   coverage,
                                   coverage_catch_up){
  
  # get incidence
  RSV_hosp_age_cohort           <- df_country$RSV_hosp_age_cohort
  RSV_icu_age_cohort            <- df_country$RSV_icu_age_cohort
  RSV_outpatient_age_cohort_sim <- df_country$RSV_outpatient_age_cohort_sim 
  RSV_non_ma_age_cohort_sim     <- df_country$RSV_non_ma_age_cohort_sim
  RSV_deaths_age_cohort         <- df_country$RSV_deaths_age_cohort
  wheezing_incidence_time       <- df_country$wheezing_incidence_time
  
  # Initiate averted Cases 
  # duplicate to get one column per simulation
  outpatient_cases_averted <- matrix(NA, nrow = nrow(RSV_outpatient_age_cohort_sim), ncol = config$num_sim) 
  non_ma_cases_averted     <- matrix(NA, nrow = nrow(RSV_non_ma_age_cohort_sim), ncol = config$num_sim) 
  hosp_cases_averted       <- matrix(NA, nrow = nrow(RSV_hosp_age_cohort), ncol = config$num_sim) 
  icu_cases_averted        <- matrix(NA, nrow = nrow(RSV_icu_age_cohort), ncol = config$num_sim) 
  rsv_deaths_averted       <- matrix(NA, nrow = nrow(RSV_deaths_age_cohort), ncol = config$num_sim)
  
  i_sim <- 1
  for(i_sim in 1:config$num_sim){
    outpatient_protection_by_cohort <- get_protection_by_cohort(ageEfficacy$outpatient[,i_sim],
                                                                config$intervention_start,
                                                                config$intervention_end,
                                                                coverage,
                                                                coverage_catch_up)
    outpatient_cases_averted[,i_sim] <- rowSums(RSV_outpatient_age_cohort_sim[,,i_sim] * outpatient_protection_by_cohort)
    
    nonMA_protection_by_cohort <- get_protection_by_cohort(ageEfficacy$nonMA[,i_sim],
                                                           config$intervention_start,
                                                           config$intervention_end,
                                                           coverage,
                                                           coverage_catch_up)
    non_ma_cases_averted[,i_sim] <- rowSums(RSV_non_ma_age_cohort_sim[,,i_sim] * nonMA_protection_by_cohort)
    
    hosp_protection_by_cohort <- get_protection_by_cohort(ageEfficacy$hospital[,i_sim],
                                                          config$intervention_start,
                                                          config$intervention_end,
                                                          coverage,
                                                          coverage_catch_up)
    hosp_cases_averted[,i_sim] <- rowSums(RSV_hosp_age_cohort * hosp_protection_by_cohort)
    
    icu_protection_by_cohort <- get_protection_by_cohort(ageEfficacy$icu[,i_sim],
                                                         config$intervention_start,
                                                         config$intervention_end,
                                                         coverage,
                                                         coverage_catch_up)
    icu_cases_averted[,i_sim] <- rowSums(RSV_icu_age_cohort * icu_protection_by_cohort)
    
    mort_protection_by_cohort <- get_protection_by_cohort(ageEfficacy$mortality[,i_sim],
                                                          config$intervention_start,
                                                          config$intervention_end,
                                                          coverage,
                                                          coverage_catch_up)
    rsv_deaths_averted[,i_sim] <- rowSums(RSV_deaths_age_cohort * mort_protection_by_cohort)
  }
  
  # Wheezing
  wheezing_cases_time_averted  <- wheezing_incidence_time * 0
  hosp_cases_0to1y_averted     <- hosp_cases_averted[1:12,]

  for(i_year in 1:dim(wheezing_cases_time_averted)[3]){
    wheezing_cases_time_averted[,,i_year] <- hosp_cases_0to1y_averted * wheezing_incidence_time[,,i_year]
  }

  # hospital admissions after first year of age do not account for wheezing or asthma
  wheezing_cases_complement <- matrix(0, nrow=60-12, ncol=config$num_sim)
  wheezing_cases_averted    <- rbind(wheezing_cases_time_averted[,,1], wheezing_cases_complement)

  # total cases averted
  total_cases_averted       <- outpatient_cases_averted + non_ma_cases_averted + hosp_cases_averted + icu_cases_averted
  
  # total uptake
  # note: use protection function (with efficacy 1) to obtain a cohort-specific coverage matrix
  uptake_cohort_age <- get_protection_by_cohort(c(1, rep(0, 59)),
                                                config$intervention_start,
                                                config$intervention_end,
                                                coverage,
                                                coverage_catch_up) * config$target_population / 12
  # combine incidence rates and probabilities
  df_country <- list(outpatient_cases_averted    = outpatient_cases_averted,
                     hosp_cases_averted          = hosp_cases_averted,
                     icu_cases_averted           = icu_cases_averted,
                     rsv_deaths_averted          = rsv_deaths_averted,
                     non_ma_cases_averted        = non_ma_cases_averted,
                     wheezing_cases_averted      = wheezing_cases_averted,
                     wheezing_cases_time_averted = wheezing_cases_time_averted,
                     total_cases_averted         = total_cases_averted,
                     uptake_cohort_age           = uptake_cohort_age
  )
}

# start from time_age matrix, convert to cohort matrix and define preventable burden
# throw warning if preventable burden of disease spans more than 1 calendar year
get_preventable_burden_catch_up <- function(f_time_age,
                                            intervention_start,
                                            intervention_end,
                                            dur_protection = NA){
  
  # transform in cohort-specific matrix
  full_burden <- time_age_to_cohort_age(f_time_age)
  
  # make sure the "dur_protection" parameter is valid 
  if((any(is.na(dur_protection))) || any(dur_protection>nrow(full_burden))){
    dur_protection <- nrow(full_burden)
  } else{
    dur_protection <- floor(dur_protection)[1]
  }
  
  # define preventable burden
  preventable_burden <- full_burden
  if(intervention_start > intervention_end){
    months_intervention <- unique(c(intervention_start:12, 1:intervention_end))
  } else{
    months_intervention <- unique(intervention_start:intervention_end)
  }
  for(i_cohort in 1:ncol(preventable_burden)){
    if(i_cohort %in% months_intervention){
      preventable_burden[-(1:(dur_protection)),i_cohort] <- 0
    } else{
      preventable_burden[-(1:dur_protection) - elapsed_months(i_cohort,intervention_start),i_cohort] <- 0
    }
  }
  
  # check time horizon
  if(sum(preventable_burden[-(1:12),])>0){
    print("WARNING: Preventable burden of disease with catch-up program covers ages above 1 year! 
              This interferes with the model implementation on discounting, which assumes that the
              modelling time horizon aligns with the age of the cohort. This is OK with a catch up 
              program if all ages remains below 1 year of age. In this case, the burden of disease
              above 1 year of age for the catch up cohorts is underestimated.")
  }
  
  return(preventable_burden)
}

get_preventable_burden_season <- function(f_time_age, intervention_start, 
                                          intervention_end, dur_protection = NA){
  
  # transform in cohort specific matrix
  full_burden <- time_age_to_cohort_age(f_time_age)
  
  # define preventable burden
  preventable_burden <- full_burden
  if(intervention_start > intervention_end){
    preventable_burden[,(intervention_end+1):(intervention_start-1)] <- 0
  } else{
    preventable_burden[,-(intervention_start:intervention_end)] <- 0
  }
  
  # account for the duration of protection
  if(any(!is.na(dur_protection))){
    preventable_burden[-(1:dur_protection),] <- 0
  }
  
  return(preventable_burden)
}

time_age_to_cohort_age <- function(f_time_age){
  
  # define maximum age in years
  max_age_yr <- floor(nrow(f_time_age)/12)
  
  # init new matrix
  cohort_age <- f_time_age*0
  
  # fill new matrix
  for(i_yr in 1:max_age_yr){
    month_index <- (1:12) + 12 * (i_yr-1)
    cohort_age[month_index,1] <- diag(f_time_age[month_index,])
    for(i in 2:12){
      cohort_age[month_index,i] = diag(f_time_age[month_index, c(i:12,1:(i-1))])
    }
  }
  # add column names
  colnames(cohort_age) <- colnames(f_time_age)
  
  # return result
  return(cohort_age)
}


# convert month id into date object equal to the 1st of the month in 2022
month_to_date <- function(month_id){
  return(as.POSIXlt(as.Date(paste0('2022-',month_id,'-01'))))
}

month_id2name <- function(month_id){
  # return(format(month_to_date(month_id),'%b'))
  month_names <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  return(month_names[month_id])
}

# months between the 1e of "month start" and the 1e of "month_end" (i.e. excluding month_end)
elapsed_months <- function(month_start,month_end) {
  date_end   <- month_to_date(month_end)
  date_start <- month_to_date(month_start)
  if(month_end < month_start){
    date_end$year <- date_end$year+1
  }
  abs(12 * (date_end$year - date_start$year) + (date_end$mon - date_start$mon))
}

get_protection_by_cohort <- function(f_time,
                                     intervention_start,
                                     intervention_end,
                                     coverage = 1,
                                     coverage_catch_up = 1){
  
  # default: year round and seasonal programs
  protection_by_cohort <- matrix(nrow=length(f_time), ncol=12,f_time)

  if(intervention_start > intervention_end){
    months_intervention <- unique(c(intervention_start:12,1:intervention_end))
  } else{
    months_intervention <- unique(intervention_start:intervention_end)
  }
  
  for(i_cohort in 1:ncol(protection_by_cohort)){
    if(!i_cohort %in% months_intervention){
      protection_by_cohort[,i_cohort]  <- 0
      if(coverage_catch_up > 0){
        time_shift <- elapsed_months(i_cohort, intervention_start)
        protection_by_cohort[-(1:time_shift),i_cohort] <- f_time[1:(length(f_time)-time_shift)] * coverage_catch_up
      }
    } else{
      protection_by_cohort[,i_cohort] <- protection_by_cohort[,i_cohort] * coverage
    }
  }
  
  # check time horizon
  if(sum(protection_by_cohort[-(1:12),])>0 & !exists('bool_warning_printed')){
    print("WARNING: Preventable burden of disease with catch-up program covers ages above 1 year! 
              This interferes with the model implementation on discounting, which assumes that the
              modelling time horizon aligns with the age of the cohort. This is OK with a catch up 
              program if all ages remains below 1 year of age. In this case, the burden of disease
              above 1 year of age for the catch up cohorts is underestimated.")
    bool_warning_printed <<- TRUE
  }
  
  return(protection_by_cohort)
}

