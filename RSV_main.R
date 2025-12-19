############################################################################ #
# This file is part of the McMarcel framework
# 
# => MAIN SCRIPT TO RUN RSV ANALYSES FOR BELGIUM
#
#  Copyright 2025, CHERMID, UNIVERSITY OF ANTWERP
############################################################################ #

# clear workspace
rm(list=ls())

## set working directory (or open RStudio with this script)
# setwd("C:\User\path\to\the\rcode\folder") ## WINDOWS
# setwd("/Users/path/to/the/rcode/folder") ## MAC

# (re)load packages and functions
source('functions/RSV_load_all.R')

####################### #
## SETTINGS          ####
####################### #

# Number of stochastic simulations (should be ≥2)
num_sim <- 20

# RUN TAG
# => to find the configuration file at ./config/<run_tag>.csv 
run_tag <- 'RSV_BELX_basecase'           # base case settings 

# random number generator seed
# note: does not necessarily gives the exact same result when used on another system (e.g. other computer)
rng_seed                  <- 20190118

####################### #
## MODEL SETUP ----
####################### #

print("****** START MC MARCEL ******")
print(paste("WORK DIR:", getwd()))

# output directory postfix
output_dir_postfix <- paste0(run_tag, '_n', num_sim)

# add time stamp (month/day/hour/min/sec) to output directory name
output_dir <- paste0('output/', format(Sys.time(), '%m%d%H%M%S_'), output_dir_postfix) 

# config file name
config_filename <- paste0('./config/', run_tag, '.xlsx')

# check number of simulations
if(num_sim<=1){
  stop("Number of simulations ('num_sim') should be >1")
}
  
####################### #
## LOAD CONFIG ----
####################### #

# load config file in csv format
sim_config_matrix <- read.xlsx(config_filename, sheet = 'scenarios', na.strings = "-")

# output file name
sim_output_filename <- file.path(output_dir, gsub('/', '_', run_tag))

# add simulation details
sim_config_matrix$num_sim       <- num_sim
sim_config_matrix$scenario_id   <- 1:dim(sim_config_matrix)[1]
if(!'rng_seed' %in% names(sim_config_matrix)) {
  sim_config_matrix$rng_seed    <- rng_seed  
}

# check if output dir exists, and create one if not
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

# Count number of scenarios (i.e. number of interventions compared, excluding comparator ('current practice'))
num_scen <- length(sim_config_matrix$scenario_id)

####################################### #
## BURDEN ----
####################################### #

# An advantage of using ‘foreach’ instead of ‘for’ is that it has an automatic option to combine the result.
# The method is specified via the ".combine” function parameter.

print('PROCESSING BURDEN [FOREACH]')
i_scen <- 1
sim_output <- foreach(i_scen = 1:num_scen,
                      .combine = 'rbind') %do%
{
  # print progress
  print(paste0(i_scen,'/', num_scen))
  
  # run get_burden function
  run_output_long <- get_burden(sim_config_matrix[i_scen,])
  
  # add scenario id
  run_output_long$scenario_id <- i_scen
  
  # return results
  return(run_output_long)
}

# print progress to the console
print('GET BURDEN COMPLETE')

################################ #
## POST-PROCESSING ----
################################ #

# add config details to output
print('MERGE INPUT & OUTPUT')
sim_output_num_row <- nrow(sim_output)
sim_output         <- merge(sim_config_matrix, sim_output, all=FALSE)
if(sim_output_num_row != nrow(sim_output)){
  print("ISSUE WHEN MERGING SIMULATION OUTPUT WITH CONFIGURATION FILE")
}

# sort output on scenario_id
sim_output <- sim_output[order(sim_output$scenario_id),]

# save full output as rds file
saveRDS(sim_output, file = paste0(sim_output_filename, '.rds'))

# print progress to the console
print('POST-PROCESSING COMPLETE')

############################### #
## CONSISTENCY CHECK ----
############################### #

if(any(sim_output$target_population     < sim_output$FVP_total) ||
   any(sim_output$total_cases_averted   < sim_output$hosp_cases_averted, na.rm=T) ||
   any(sim_output$hosp_cases_averted    > sim_output$hosp_cases_reference, na.rm=T) || 
   any(sim_output$rsv_deaths_averted    > sim_output$rsv_deaths_reference, na.rm=T)){
  print('INCONSISTENT OUTPUT: POPULATION ~ CASES ~ AVERTED CASES')
}

####################### #
## RESULTS ----
####################### #

# CEA results 
plot_CEA_results(sim_output_filename)

# aggregated global and country tables
write_output_tables(sim_output_filename)
write_summary_tables(sim_output_filename)

# print progress to the console
print('SCRIPT COMPLETE!')


  
  ## CHECK
  
  out_files <- dir(dirname(sim_output_filename),full.names = TRUE)
  
  ref_files <- dir('output/export_reference/',full.names = TRUE)
  if(any(basename(out_files) != basename(ref_files))){
    stop('output file names changed')
  }
  
  
  xx <- readRDS(out_files[grepl('.rds',out_files)])
  yy <- readRDS(ref_files[grepl('.rds',ref_files)])
  
  yy$outputFileDir <- NULL
  
  if(any(xx != yy)){
    stop('output changed')
  }
  
  print('TEST OK')
