############################################################################ #
# This file is part of the RSV modelling project.
# 
# => BURDEN FUNCTION FOR THE RSV ANALYSIS
#
#  Copyright 2025, CHERMID, UNIVERSITY OF ANTWERP
############################################################################ #

# code to debug the function, is never executed by default
if(0==1)
{
  # for debugging starting from RSV_main.R, run these 2 lines, and then run everything in function starting from set.seed 
    i_scen <- 5  #intervention option considered
    configList <- sim_config_matrix[i_scen,] 
}
  
# function to calculate the burden
get_burden <- function (configList) {

  # set RNG seed
  set.seed(configList$rng_seed)
  
  ################################ #
  # Configuration               ####
  ################################ #
  
  # get model parameters (based on country, year, num_sim, etc... as specified in configList)
  config <- get_rsv_ce_config(configList) 

	################################ #
	# Pre-processing              ####
	################################ #
	
  # get 'mAb' target population (one value)
  config$target_population_lifebirth <- config$target_population * (1 - config$stillbirth_rate)

  # use the nYearsOfAges to construct age-dependent vectors
  config$nMonthsOfAges   <- config$nYearsOfAges*config$monthsInYear
  
	# births => cohort size
	# get the life table for the given country and period 
  # by age in months (for 100 months) and includes lx, le, nMx, lx_rate, le_disc
	#life_table <- get_life_table(config$country_iso,config$year,config$disc_rate_effect) # generic
	life_table <- get_life_table_BEL2023(config$disc_rate_effect)                         # for Belgium
  
	# select the target ages [0-59 months]  
	life_table <- life_table[1:config$nMonthsOfAges,] 
	
	# transfer mortality-related QALY losses
  config$hosp_mortality_QALYloss      <- life_table$LE
  config$hosp_mortality_QALYloss_disc <- life_table$LE_disc
  config$life_table <- life_table
  
  ############################## #
  # Intervention: mAb  ####
  ############################## #
  
  ageEfficacy_mAb <- get_protection(config              = config,
                                    intervention        = 'mAb')
  
  # keep track of protection values
  config$efficacy_mAb_nonMA          <- ageEfficacy_mAb$nonMA[1,]
  config$efficacy_mAb_outpatient     <- ageEfficacy_mAb$outpatient[1,]
  config$efficacy_mAb_hospital       <- ageEfficacy_mAb$hospital[1,]
  config$efficacy_mAb_icu            <- ageEfficacy_mAb$icu[1,]
  config$efficacy_mAb_mortality      <- ageEfficacy_mAb$mortality[1,]
  
  ############################################# #
  # Intervention: Maternal immunisation ####
  ############################################# #
  
  ageEfficacy_maternal <- get_protection(config       = config,
                                         intervention = 'maternal')
  
  # keep track of protection values
  config$efficacy_maternal_nonMA          <- ageEfficacy_maternal$nonMA[1,]
  config$efficacy_maternal_outpatient     <- ageEfficacy_maternal$outpatient[1,]
  config$efficacy_maternal_hospital       <- ageEfficacy_maternal$hospital[1,]
  config$efficacy_maternal_icu            <- ageEfficacy_maternal$icu[1,]
  config$efficacy_maternal_mortality      <- ageEfficacy_maternal$mortality[1,]
  
  ###################################### #
	# Disease incidence  #### 
	###################################### #
  
  # Health care use 
  df_country <- get_incidence_data(config)
  
  ###################################### #
  # Prevented Disease incidence  #### 
  ###################################### #
  
  # mAb-related prevented burden 
  df_prevented_mAb <- get_prevented_incidence(config,
                                              df_country,
                                              ageEfficacy_mAb,
                                              coverage = config$coverage_mAb,
                                              coverage_catch_up = config$coverage_catch_up) 
  
  # maternal-related prevented burden 
  df_prevented_maternal <- get_prevented_incidence(config,
                                                df_country,
                                                ageEfficacy_maternal,
                                                coverage = config$coverage_maternal * (1 - config$pre_term_rate),
                                                coverage_catch_up = 0) 
  
  # get total prevented incidence
  df_prevented <- sum_lists_by_element(df_prevented_mAb, df_prevented_maternal)

  ###################################### #
  # Burden of disease  #### 
  ###################################### #
  
  # sample QALY loss due to (non) severe infection
  #note: severe estim. counts for both hospital and ICU cases 
	config$QALYloss_hospital    <- sample_rgamma(config$QALYloss_hosp_mean, config$QALYloss_hosp_se,config$num_sim)
  config$QALYloss_outpatient  <- sample_rgamma(config$QALYloss_outpatient_mean, config$QALYloss_outpatient_se,config$num_sim)
  
  #QALY loss due to recurrent wheezing and asthma
  config$wheezing_QALYloss_peryear <- (1-sample_rgamma(config$QALY_wheezing_mean, config$QALY_wheezing_se,config$num_sim))
  
  #QALY loss due to non-medically attended cases
  config$non_ma_QALYloss <- sample_rgamma(config$QALY_non_ma_mean, config$QALY_non_ma_se, config$num_sim)

  # option to include QALY loss for caregiver(s)
  if(config$boolean_caregiver_QALYloss){
    
    config$caregiver_QALYloss_hospital   <- sample_rgamma(config$caregiver_QALYloss_hosp_mean, config$caregiver_QALYloss_hosp_se, config$num_sim)
    config$caregiver_QALYloss_outpatient <- sample_rgamma(config$caregiver_QALYloss_outpatient_mean, config$caregiver_QALYloss_outpatient_se, config$num_sim)
    config$caregiver_non_ma_QALYloss     <- sample_rgamma(config$caregiver_QALY_non_ma_mean, config$caregiver_QALY_non_ma_se, config$num_sim)
    
    config$QALYloss_hospital   <- config$QALYloss_hospital + config$caregiver_QALYloss_hospital
    config$QALYloss_outpatient <- config$QALYloss_outpatient + config$caregiver_QALYloss_outpatient
    config$non_ma_QALYloss     <- config$non_ma_QALYloss + config$caregiver_non_ma_QALYloss
  }
  
	################################### #
	# Costs   ####
	################################### #
	
  # extrapolate hospital cost data and reshape into an [age, sim] structure
  config$hosp_cost <- get_rgamma_age_sim(matrix_mean_sd = rbind(c(config$db_cost$cost_hosp_under1y_mean, config$db_cost$cost_hosp_under1y_sd),
                                                                c(config$db_cost$cost_hosp_1to4y_mean, config$db_cost$cost_hosp_1to4y_sd)),
                                         age_breaks     = c(0,11,59),
                                         num_sim        = config$num_sim)
	
	# ICU cost
	config$icu_cost <- get_rgamma_age_sim(matrix_mean_sd = rbind(c(config$db_cost$cost_icu_under1y_mean, config$db_cost$cost_icu_under1y_sd),
	                                                             c(config$db_cost$cost_icu_1to4y_mean, config$db_cost$cost_icu_1to4y_sd)),
	                                      age_breaks     = c(0,11,59),
	                                      num_sim        = config$num_sim)
	
	
	# outpatient
	config$outpatient_cost <- get_rgamma_age_sim(matrix_mean_sd = rbind(c(config$db_cost$cost_outpatient_under1y_mean, config$db_cost$cost_outpatient_under1y_sd), 
	                                                                    c(config$db_cost$cost_outpatient_1y_mean, config$db_cost$cost_outpatient_1y_sd),
	                                                                    c(config$db_cost$cost_outpatient_2y_mean, config$db_cost$cost_outpatient_2y_sd),
	                                                                    c(config$db_cost$cost_outpatient_3y_mean, config$db_cost$cost_outpatient_3y_sd),
	                                                                    c(config$db_cost$cost_outpatient_4y_mean, config$db_cost$cost_outpatient_4y_sd)),
	                                             age_breaks     = c(0,(1:5)*12-1),
	                                             num_sim        = config$num_sim)
	
	### Intervention costs: maternal
	config$price_dose_maternal   <- runif(config$num_sim, config$price_dose_maternal_min, config$price_dose_maternal_max)
	uptake_maternal_cohort_age   <- df_prevented_maternal$uptake_cohort_age  * (1 + config$stillbirth_rate) # pregnancies = newborns + stillbirths
	config$unit_price_maternal   <- config$admin_cost_maternal + (config$price_dose_maternal * (1 + config$wastage_maternal))
	if(sum(uptake_maternal_cohort_age) > 0){
	  config$unit_price_maternal <- config$price_dose_maternal  + config$implementation_cost_maternal	/ sum(df_prevented_maternal$uptake_cohort_age)
	}
	cost_maternal <- multiply_age_by_sim(rowSums(uptake_maternal_cohort_age), config$unit_price_maternal)
	cost_maternal_delivery <- sum(uptake_maternal_cohort_age) * config$admin_cost_maternal
	
	### Intervention costs: mAb
	config$price_dose_mAb <- runif(config$num_sim, config$price_dose_mAb_min, config$price_dose_mAb_max)
	uptake_mAb_cohort_age <- df_prevented_mAb$uptake_cohort_age  # target pop = newborns
	config$unit_price_mAb <- config$admin_cost_mAb + (config$price_dose_mAb * (1 + config$wastage_mAb))
	if(sum(uptake_mAb_cohort_age) > 0){
	  config$unit_price_mAb      <- config$unit_price_mAb  + config$implementation_cost_mAb	/ sum(uptake_mAb_cohort_age)
	 }
	cost_mAb <- multiply_age_by_sim(rowSums(uptake_mAb_cohort_age), config$unit_price_mAb)
	cost_mAb_delivery <- sum(uptake_mAb_cohort_age) * config$admin_cost_mAb
	
	# total intervention cost
	config$intervention_cost <- (cost_maternal + cost_mAb)
	
	# implementation cost
	config$implementation_cost_unc  <- rep(any(uptake_mAb_cohort_age > 0) * config$implementation_cost_mAb + 
	                                       any(uptake_maternal_cohort_age > 0) * config$implementation_cost_maternal, 
	                                       config$num_sim)  
	
	##  Productivity costs: societal perspective
	if(config$bool_productivity_cost){
	  config$unit_cost_productivity_loss_outpatient <- get_rgamma_age_sim(matrix_mean_sd = rbind(c(0, 0), 
            	                                                                                 c(config$productivity_cost_outpatient_3to11m_mean, config$productivity_cost_outpatient_3to11m_se),
            	                                                                                 c(config$productivity_cost_outpatient_12to47m_mean, config$productivity_cost_outpatient_12to47m_se),
            	                                                                                 c(0,0)),
            	                                                          age_breaks     = c(0, 3, 11, 47, 59),
            	                                                          num_sim        = config$num_sim)
	  config$unit_cost_productivity_loss_inpatient  <- get_rgamma_age_sim(matrix_mean_sd = rbind(c(0, 0), 
            	                                                                                 c(config$productivity_cost_inpatient_3to11m_mean, config$productivity_cost_inpatient_3to11m_se),
            	                                                                                 c(config$productivity_cost_inpatient_12to47m_mean, config$productivity_cost_inpatient_12to47m_se),
            	                                                                                 c(0, 0)),
            	                                                          age_breaks     = c(0, 3, 11, 47, 59),
            	                                                          num_sim        = config$num_sim)
	} else {
	  config$unit_cost_productivity_loss_outpatient <- 0
	  config$unit_cost_productivity_loss_inpatient  <- 0
	}
	

	############################### #
	# Discounting                ####
	############################### #
	
	# set-up discounting-over-time vector (in months, i.e. same value for 12 months of each year)
	years_after_vaccin       <- rep(0:(config$nYearsOfAges - 1), each = config$monthsInYear)
	config$disc_time_effect  <- 1/((1 + config$disc_rate_effect)^years_after_vaccin)
	config$disc_time_cost    <- 1/((1 + config$disc_rate_cost)^years_after_vaccin)
	
	# set-up discounting-over-time vector (in years, for discounting long-term burden asthma and wheezing)
	lifeyears_vector              <- 0:99
	config$disc_time_effect_yr    <- 1/((1 + config$disc_rate_effect)^lifeyears_vector)
	config$disc_time_cost_yr      <- 1/((1 + config$disc_rate_cost)^lifeyears_vector) 
	
	config$QALYloss_outpatient_disc <- get_discounted_life_years(config$QALYloss_outpatient, config$disc_rate_effect)
  config$QALYloss_hospital_disc   <- get_discounted_life_years(config$QALYloss_hospital, config$disc_rate_effect)
  config$non_ma_QALYloss_disc     <- get_discounted_life_years(config$non_ma_QALYloss, config$disc_rate_effect)
	
  ############################### #
  # Burden (prevented)    ####
  ############################### #
  
  # get burden estimates for default burden
  df_country_all <- calculate_burden_estimates(df_burden = df_country, config = config)
  
  # get averted burden estimates
  df_averted_all <- calculate_burden_estimates(df_burden = df_prevented, config = config)
  
	# obtain estimated burden with intervention in place
	df_program_all <- sum_lists_by_element(df_country_all, df_averted_all, bool_negation = TRUE)
	
	# incremental cost = intervention cost - averted costs over time
	df_program_all$intervention_cost        <- config$intervention_cost
	df_program_all$incremental_medical_cost <- config$intervention_cost - df_averted_all$total_medical_cost
	df_program_all$incremental_cost         <- config$intervention_cost - df_averted_all$total_cost_societal
	
	# discounted incremental costs
	# note: discounting has no impact on intervention cost if it takes place in the 1st year
	df_program_all$intervention_cost_disc           <- config$intervention_cost * config$disc_time_cost
	df_program_all$incremental_medical_cost_disc    <- df_program_all$intervention_cost_disc - df_averted_all$total_medical_cost_disc
	df_program_all$incremental_cost_disc            <- df_program_all$intervention_cost_disc - df_averted_all$total_medical_cost_disc - df_averted_all$total_cost_productivity_loss_disc
	
	################################# #
	## AGGREGATE OUTPUT                      ####
	################################# #
	## FUNCTION TO AGGREGATE BURDEN OUTPUT
	# note: use function to get age-specific results, without duplicating code
	#age_function=select_0to5y;age_tag=''
	
	aggregate_output <- function(age_function,age_tag){
	  
	  # start with empty data.frame, only specify number of rows
	  output_all <- data.frame(row.names = 1:config$num_sim)
	  
	  # add aggregated reference age-group specific estimates
	  i_elem <- names(df_country_all)[2]
	  for(i_elem in names(df_country_all)){
	    output_all[paste0(i_elem,'_reference')] <- colSums(age_function(df_country_all[[i_elem]]))
	    output_all[paste0(i_elem,'_averted')]   <- colSums(age_function(df_averted_all[[i_elem]]))
	  }
	  
	  # add aggregated program age-group specific estimates
	  for(i_elem in names(df_program_all)){
	    output_all[paste0(i_elem,'_program')]  <- colSums(age_function(df_program_all[[i_elem]]))
	  }
	  
	  # include also reference incremental and intervention cost, equal to 0
	  colnames_incr <- names(df_program_all)[grepl('incremental', names(df_program_all)) | grepl('intervention', names(df_program_all))]
	  output_all[paste0(colnames_incr, '_reference')]  <- output_all[, paste0(colnames_incr, '_program')] * 0
	  
	  # adjust column names
	  if(nchar(age_tag)>0){
	    
	    # get base names (e.g. total_cases from "total_cases_disc_averted" or "total_cases")
	    b_names <- gsub('_disc', '', names(output_all))
	    b_names <- unique(gsub('_averted', '', b_names))
	    b_names <- unique(gsub('_reference', '', b_names))
	    b_names <- unique(gsub('_program', '', b_names))
	    
	    # add age category (e.g. "total_cases_0to1y")
	    b_names_age <- paste(b_names, age_tag, sep = '_')
	    
	    # replace base name by age-specific name
	    for(i_name in 1:length(b_names)){
	      names(output_all) <- gsub(b_names[i_name], b_names_age[i_name], names(output_all))
	    }
	  }
	  
	  # return
	  return(output_all)
	} 
	
	output_all <- data.frame(aggregate_output(select_0to5y,''), 
	                         aggregate_output(select_0to1y,'0to1y'),
	                         aggregate_output(select_1to5y,'1to5y'),
	                         aggregate_output(select_0to2mo,'0to2mo'),
	                         aggregate_output(select_3to5mo,'3to5mo'), 
	                         aggregate_output(select_6to11mo,'6to11mo'))  
	
	# add tags
	output_all$intervention    <- config$intervention 
	output_all$season          <- config$season 
	output_all$config_tag      <- config$config_tag 
	
	# add input parameters
	output_all$perspective        <- config$perspective  
	output_all$outpatient_cost    <- config$outpatient_cost[1,]  
	output_all$non_ma_cost        <- config$non_ma_cost 
	output_all$hosp_cost          <- config$hosp_cost[1,]    # (NOTE! cost is NOT the same for all age groups)
	output_all$wheezing_cost      <- config$wheezing_cost_peryear  
	output_all$price_dose_maternal_min<- config$price_dose_maternal_min  
	output_all$price_dose_mAb_min     <- config$price_dose_mAb_min 
	output_all$price_dose_maternal_max<- config$price_dose_maternal_max 
	output_all$price_dose_mAb_max     <- config$price_dose_mAb_max 
	output_all$price_dose_maternal<- config$price_dose_maternal  
	output_all$price_dose_mAb     <- config$price_dose_mAb  
	output_all$admin_cost_maternal<- config$admin_cost_maternal 
	output_all$admin_cost_mAb     <- config$admin_cost_mAb  
	output_all$implementation_cost<- config$implementation_cost_unc
	output_all$dur_prot_maternal  <- config$dur_prot_maternal
	output_all$dur_prot_mAb       <- config$dur_prot_mAb
	output_all$rsv_outpatient_agedist_0months    <- df_country$outpatient_agedist_under6months[1,]
	output_all$rsv_outpatient_agedist_1months    <- df_country$outpatient_agedist_under6months[2,]
	output_all$rsv_outpatient_agedist_2months    <- df_country$outpatient_agedist_under6months[3,]
	output_all$rsv_outpatient_agedist_3months    <- df_country$outpatient_agedist_under6months[4,]
	output_all$rsv_outpatient_agedist_4months    <- df_country$outpatient_agedist_under6months[5,]
	output_all$rsv_outpatient_agedist_5months    <- df_country$outpatient_agedist_under6months[6,]
	output_all$rsv_wheezing_prob_under1  <- df_country$prob_wheezing_under1
	output_all$proportion_non_ma_mean    <- df_country$proportion_non_ma_mean
	
	output_all$QALYloss_hospital         <- config$QALYloss_hospital 
	output_all$QALYloss_outpatient       <- config$QALYloss_outpatient  
	output_all$QALYloss_wheezing         <- config$wheezing_QALYloss_peryear
	output_all$QALYloss_non_ma           <- config$non_ma_QALYloss
	output_all$FVP_maternal              <- config$target_population * config$coverage_maternal * config$proportion_target_pop 
	output_all$FVP_mAb                   <- config$target_population_lifebirth * config$coverage_mAb * config$proportion_target_pop 
	output_all$FVP_total                 <- output_all$FVP_maternal + output_all$FVP_mAb
	output_all$population_lifebirth      <- round(config$target_population_lifebirth) 
	output_all$target_population         <- config$target_population 
	
	# add cost # only used under 1 year, as the intervention only protection for 1st RSV season 
	output_all$cost_hosp_under1y         <- colMeans(config$hosp_cost[1:12,])
	output_all$cost_icu_under1y          <- colMeans(config$icu_cost[1:12,])
	output_all$cost_outpatient_under1y   <- colMeans(config$outpatient_cost[1:12,]) 
	
	output_all$cost_maternal_delivery    <- cost_maternal_delivery
  output_all$cost_mAb_delivery         <- cost_mAb_delivery 
  
  output_all$efficacy_maternal_outpatient   <- config$efficacy_maternal_outpatient 
  output_all$efficacy_maternal_hospital     <- config$efficacy_maternal_hospital 
  output_all$efficacy_maternal_icu          <- config$efficacy_maternal_icu 
  output_all$efficacy_maternal_mortality    <- config$efficacy_maternal_mortality
  output_all$efficacy_mAb_outpatient        <- config$efficacy_mAb_outpatient 
  output_all$efficacy_mAb_hospital          <- config$efficacy_mAb_hospital 
  output_all$efficacy_mAb_icu               <- config$efficacy_mAb_icu 
  output_all$efficacy_mAb_mortality         <- config$efficacy_mAb_mortality 
  output_all$efficacy_total_num_na          <- config$total_efficacy_num_na
  
  output_all$stochastic_proc_id  <- 1:config$num_sim

	# return results
	return(output_all)
}

calculate_burden_estimates <- function(df_burden, config){
  
  # adjust column names if needed
  names(df_burden) <- gsub('_averted','',names(df_burden))
  
  # start from estimated cases and deaths only
  df_burden <- df_burden[grepl('cases',names(df_burden)) | grepl('deaths$',names(df_burden))]
  names(df_burden)
  
  # order names
  df_burden <- df_burden[order(names(df_burden))]
  
  # initiate dummy to extend time horizon to 5 years
  zero_cases_over1 <- matrix(0,nrow=60-12,ncol=config$num_sim)
  
  # store wheezing_cases_time outside df_burden
  # note: to prevent different dimensions in output
  wheezing_cases_time <- df_burden$wheezing_cases_time
  df_burden$wheezing_cases_time <- NULL
  
  ############################### #
  # Costs       
  ###############################	#
  
  # medical costs: baseline
  df_burden$cost_rsv_outpatient    <- df_burden$outpatient_cases    * config$outpatient_cost 
  df_burden$cost_rsv_non_ma        <- df_burden$non_ma_cases        * config$non_ma_cost     
  df_burden$cost_rsv_hosp          <- df_burden$hosp_cases          * config$hosp_cost  
  df_burden$cost_rsv_icu           <- df_burden$icu_cases           * config$icu_cost  
  cost_rsv_wheezing_under1         <- wheezing_cases_time           * config$wheezing_cost_peryear
  df_burden$cost_rsv_wheezing      <- rbind(aggregate_cohort_by_year(cost_rsv_wheezing_under1),zero_cases_over1)
  
  # time loss related costs
  df_burden$cost_productivity_loss_outpatient <- df_burden$outpatient_cases * config$unit_cost_productivity_loss_outpatient 
  df_burden$cost_productivity_loss_hospital   <- df_burden$hosp_cases * config$unit_cost_productivity_loss_inpatient
  df_burden$cost_productivity_loss_icu        <- df_burden$icu_cases * config$unit_cost_productivity_loss_inpatient
  
  ####################### #
  # DISCOUNTING   
  ####################### #
  
  for(i_elem in names(df_burden)){
      if(grepl('cost',i_elem)){
        elem_disc <- df_burden[[i_elem]] * config$disc_time_cost
      } else{
        elem_disc <- df_burden[[i_elem]] * config$disc_time_effect
      }
      df_burden[[paste0(i_elem,'_disc')]] = elem_disc
  }
  
  # special case for wheezing
  cost_rsv_wheezing_under1_disc        <- multiply_cohort_year(cost_rsv_wheezing_under1, config$disc_time_effect_yr)
  df_burden$cost_rsv_wheezing_disc     <- rbind(aggregate_cohort_by_year(cost_rsv_wheezing_under1_disc), zero_cases_over1)
  
  ####################### #
  # QALY   
  ####################### #
  
  # QALY
  # YLD = years lived with disability
  df_burden$outpatient_YLD        <- multiply_age_sim(df_burden$outpatient_cases, config$QALYloss_outpatient)
  df_burden$non_ma_YLD            <- multiply_age_sim(df_burden$non_ma_cases, config$non_ma_QALYloss) 
  df_burden$hosp_YLD              <- multiply_age_sim(df_burden$hosp_cases, config$QALYloss_hospital)
  df_burden$icu_YLD               <- multiply_age_sim(df_burden$icu_cases, config$QALYloss_hospital) 
  
  # Wheezing and asthma 
  wheezing_YLD_under1             <- multiply_cohort_sim(wheezing_cases_time, config$wheezing_QALYloss_peryear)
  df_burden$wheezing_YLD          <- rbind(aggregate_cohort_by_year(wheezing_YLD_under1), zero_cases_over1)
  
  # Discounted burden using discounted QALY estimates and time preference
  df_burden$outpatient_YLD_disc   <- multiply_age_sim(df_burden$outpatient_cases, config$QALYloss_outpatient_disc) * config$disc_time_effect
  df_burden$non_ma_YLD_disc       <- multiply_age_sim(df_burden$non_ma_cases, config$non_ma_QALYloss_disc)         * config$disc_time_effect 
  df_burden$hosp_YLD_disc         <- multiply_age_sim(df_burden$hosp_cases, config$QALYloss_hospital_disc)         * config$disc_time_effect
  df_burden$icu_YLD_disc          <- multiply_age_sim(df_burden$icu_cases, config$QALYloss_hospital_disc)          * config$disc_time_effect 
  
  # Discounted wheezing and asthma 
  wheezing_YLD_under1_disc        <- multiply_cohort_year(wheezing_YLD_under1, config$disc_time_effect_yr)
  df_burden$wheezing_YLD_disc     <- rbind(aggregate_cohort_by_year(wheezing_YLD_under1_disc), zero_cases_over1)
  
  ####################### #
  # TOTALS   
  ####################### #
  
  # aggregate (discounted) medical cost, productivity cost, and YLD
  out_list <- c('cost_rsv','cost_productivity_loss','YLD')
  for(i_out in out_list){
    for(i_disc in c(FALSE,TRUE)){
      sel_elem <- grepl(i_out, names(df_burden)) & 
                    !grepl('total_', names(df_burden)) &
                    (grepl('_disc', names(df_burden)) == i_disc)
      name_out <- paste0('total_', i_out, ifelse(i_disc,'_disc',''))
      df_burden[[name_out]] <- sum_elements(df_burden[sel_elem])
    }
  }

  # update medical cost variable name
  names(df_burden) <- gsub('total_cost_rsv','total_medical_cost',names(df_burden))
  
  # total societal cost
  df_burden$total_cost_societal                  <- df_burden$total_medical_cost + df_burden$total_cost_productivity_loss 
  df_burden$total_cost_societal_disc             <- df_burden$total_medical_cost_disc + df_burden$total_cost_productivity_loss_disc 
  
  # Life years lost => YLL 
  df_burden$total_YLL      <- df_burden$rsv_deaths * config$life_table$LE 
  df_burden$total_YLL_disc <- df_burden$rsv_deaths * config$life_table$LE_disc * config$disc_time_effect 
  
  df_burden$total_HrQA_YLL       <- df_burden$rsv_deaths * config$life_table$HrQALE 
  df_burden$total_HrQA_YLL_disc  <- df_burden$rsv_deaths * config$life_table$HrQALE_disc * config$disc_time_effect 
  
  # QALY loss  
  df_burden$total_QALY      <- df_burden$total_YLD + df_burden$total_HrQA_YLL
  df_burden$total_QALY_disc <- df_burden$total_YLD_disc + df_burden$total_HrQA_YLL_disc	
  
  ####################### #
  # RETURN   
  ####################### #
  
  # return
  return(df_burden)
}


# Sample QALY values
sample_rgamma <- function(mean,stdev,num_sim){
  alpha <- (mean^2) / (stdev^2)
  beta  <- (stdev^2) / (mean)
  sample <- rgamma(num_sim, shape=alpha, rate = 1/beta)
}

## FUNCTION TO SELECT ALL BURDEN OUTPUT (dummy)
select_0to5y <- function(x){
  return(x)
}

## FUNCTION TO SELECT BURDEN OUTPUT FOR (1,5] YEARS OF AGE
select_1to5y <- function(x){
  x_1to5y         <- x*0
  x_1to5y[13:60,] <- x[13:60,]
  return(x_1to5y)
}

## FUNCTION TO SELECT BURDEN OUTPUT FOR [0,1] YEARS OF AGE
select_0to1y <- function(x){
  x_0to1y        <- x*0
  x_0to1y[1:12,] <- x[1:12,]
  return(x_0to1y)
}

## FUNCTION TO SELECT BURDEN OUTPUT FOR [0,2] MONTHS OF AGE
select_0to2mo <- function(x){
  x_0to2mo        <- x*0
  x_0to2mo[1:3,] <- x[1:3,]
  return(x_0to2mo)
} 

## FUNCTION TO SELECT BURDEN OUTPUT FOR (3,5] MONTHS OF AGE
select_3to5mo <- function(x){
  x_3to5mo        <- x*0
  x_3to5mo[4:6,] <- x[4:6,]
  return(x_3to5mo)
}

## FUNCTION TO SELECT BURDEN OUTPUT FOR (6,11] MONTHS OF AGE
select_6to11mo <- function(x){
  x_6to11mo        <- x*0
  x_6to11mo[7:12,] <- x[7:12,]
  return(x_6to11mo)
} 

# multiply age-specific vector with simulation-specific vector
multiply_age_by_sim <- function(vector_age, vector_sim){
  matrix1_age_sim <- matrix(rep(vector_age, length(vector_sim)), ncol = length(vector_sim), byrow = FALSE)
  matrix2_age_sim <- matrix(rep(vector_sim, length(vector_age)), ncol = length(vector_sim), byrow = TRUE)
  return(matrix1_age_sim * matrix2_age_sim)
}

# multiply age-sim matrix with simulation-specific vector (i.e. by row)
multiply_age_sim <- function(matrix_age_sim, vector_sim){
  return(matrix_age_sim * rep(vector_sim[1:ncol(matrix_age_sim)],
                              each = nrow(matrix_age_sim)))
}

# multiply age-sim (or age-cohort) matrix with age-specific vector (i.e. by column)
multiply_matrix_by_age <- function(matrix_age_sim,vector_age){
  return(matrix_age_sim * rep(vector_age,
                              times = ncol(matrix_age_sim)))
}

# multiply age-sim-year array with simulation-specific vector
multiply_cohort_sim <- function(array_age_sim_year,vector_sim){
  return(array_age_sim_year * rep(vector_sim[1:dim(array_age_sim_year)[2]],
                                  times = dim(array_age_sim_year)[3],
                                  each  = dim(array_age_sim_year)[1]))
}

# multi 3D age-sim-year array with year specific vector
multiply_cohort_year <- function(array_age_sim_year,vector_year){
  return(array_age_sim_year * rep(vector_year[1:dim(array_age_sim_year)[3]],
                                  each = dim(array_age_sim_year)[1] * dim(array_age_sim_year)[2]))
}

# sum the 3D age-sim-year array by the cohort over the years included
aggregate_cohort_by_year <- function(array_age_sim_year){
  return(apply(array_age_sim_year,1:2,sum))
}

# function to obtain sampled age-specific [age,sim] matrix from age-group-specific summary statistics
get_rgamma_age_sim <- function(matrix_mean_sd,   # e.g. rbind(c(100,10),c(500,20))
                               age_breaks,       # e.g. c(0,11,59)
                               num_sim){         # e.g. 10
  # define [age_group - sim] matrix
  matrix_age_sim <- matrix(NA,ncol = num_sim, nrow = nrow(matrix_mean_sd))
  i_age <- 1
  for(i_age in 1:nrow(matrix_age_sim)){
    
    if(all(matrix_mean_sd[i_age,] == 0)){
      matrix_age_sim[i_age,] <- 0
    } else{
      matrix_age_sim[i_age,] <- sample_rgamma(matrix_mean_sd[i_age,1], matrix_mean_sd[i_age,2], num_sim)
    }
  }
  
  # extrapolate ages
  data_age_sim <- approx_age_by_sim(matrix_age_sim, age_breaks)
  
  # return
  return(data_age_sim)
}

# function to obtain an age-specific [age,sim] matrix from age-group-specific data
approx_age_by_sim <- function(matrix_age_sim,   # e.g. rbind(c(42,56,30),c(30,40,50))
                              age_breaks){      # e.g. c(0,11,59)

  # define input and output age breaks
  x_in  <- age_breaks[-length(age_breaks)]
  x_out <- seq(min(age_breaks), max(age_breaks))
  
  # define [age, sim] structure
  data_age_sim <- matrix(NA, nrow = length(x_out), ncol = ncol(matrix_age_sim))
  
  # extrapolate given values, by simulation
  for(i_sim in 1:ncol(matrix_age_sim)){
    data_age_sim[,i_sim] <- approx(x = x_in,
                                   y = matrix_age_sim[,i_sim],
                                   xout = x_out,
                                   rule = 2,
                                   method = 'constant')$y
  }
  
  # return
  return(data_age_sim)
}

# note: option to in
sum_lists_by_element <- function(list1, list2, bool_negation = FALSE){
  
  # check names and dimensions
  if(!all(names(list1) == names(list1)) && length(unlist(list1)) == length(unlist(list2))){
    stop("list_sum_by_element error: lists are not compatible")
  }
  
  # element wise sum
  for(i_elem in names(list1)){
    list1[[i_elem]] <- list1[[i_elem]] + list2[[i_elem]] * ifelse(bool_negation, -1, 1)
  }
  
  # return result
  return(list1)
}

sum_elements <- function(x_list){
  
  # check dimensions
  if(length(unique(lapply(x_list,dim))) > 1){
    stop("sum_elements error: elements are not compatible")
  }

  # sum elements
  for(i_elem in 1:length(x_list)){
    if(i_elem == 1){
      x_out <- x_list[[i_elem]]
    } else{
      x_out <- x_out + x_list[[i_elem]]
    }
  }
  
  # return result
  return(x_out)
}
