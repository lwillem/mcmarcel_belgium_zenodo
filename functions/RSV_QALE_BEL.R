############################################################################ #
# This file is part of the McMarcel framework.
# 
# => CODE TO OBTAIN HEALTH RELATED QUALITY ADJUSTED LIFE EXPECTATION
#
#  Copyright 2025, CHERMID, UNIVERSITY OF ANTWERP
############################################################################ #

popnormINDEX_BE <- function () 
{
  
  ## code to prepare the input file
  # library(devtools)
  # devtools::install_github("brechtdv/EQ5D.be")
  # library(EQ5D.be) 
  # be_qol <- popnormINDEX(age = 15:100, sex = 'B', region = 'BE', year = 2018)
  # saveRDS(be_qol, 'input/eq5d/popnormINDEX.rds')
  
  # load file
  be_qol <- readRDS('input/eq5d/popnormINDEX.rds')
  
  # extrapolate 15y old to  0:14 year olds
  be_qol <- be_qol[c(rep(1,15), 1:nrow(be_qol)), ]
  be_qol$age[1:15] <- 0:14
  
  # remove row names
  row.names(be_qol) <- NULL
  
  return(be_qol)
}

# calculate age-specific health-related (Hr) Quality-of-Life (QoL)
get_HrQALE <- function(life_expectancy_year, 
                       f_disc_rate_qaly,
                       qol_age_year = NA){
  
  # make sure the QoL vector is set
  if(any(is.na(qol_age_year)))
  {
    qol_age_year <- rep(1, length(life_expectancy_year))
  }
  
  # make sure the QoL vector has the correct dimensions
  # duplicate last values if needed
  if(length(life_expectancy_year) > length(qol_age_year))
  {
    qol_age_year <- c(qol_age_year,
                      rep(qol_age_year[length(qol_age_year)], length(life_expectancy_year) - length(qol_age_year)))
    #warning("extended the age-specific Quality-of-Life vector")
  }
  
  # initialize the output vector
  life_expectancy_HrQoL <- vector(length = length(life_expectancy_year))
  
  # loop over the ages, and apply discounting
  i <- 1
  for(i in 1:length(life_expectancy_year)){
    
    tmp_years_full      <- 1:floor(life_expectancy_year[i])
    tmp_year_incomplete <- ceiling(life_expectancy_year[i])
    tmp_year_remaining  <- life_expectancy_year[i] - floor(life_expectancy_year[i])
    
    # QoL-adjusted and discounted life expectancy
    life_expectancy_HrQoL_floor     <- sum(qol_age_year[(i - 1) + tmp_years_full] / ((1 + f_disc_rate_qaly)^tmp_years_full))
    life_expectancy_HrQoL_remaining <- sum(qol_age_year[tmp_year_incomplete] / ((1 + f_disc_rate_qaly)^tmp_year_incomplete)) * tmp_year_remaining
    life_expectancy_HrQoL[i]        <- life_expectancy_HrQoL_floor + life_expectancy_HrQoL_remaining
  }
  
  return(life_expectancy_HrQoL)
}

get_LE_year2months <- function(age, LE_year, max_months){
  num_months_year       <- 12
  life_expectancy_month <- approx(age*num_months_year, LE_year, xout = seq(0, max_months, 1))
  life_table            <- data.frame(life_expectancy_month)
  names(life_table)     <- c('age','LE')
  return(life_table)
}

get_life_table_BEL2023 <- function(f_disc_rate_qaly = 0, max_age_months = 100){
  
  # load life table
  life_table_BE   <- read.table('input/life_table_Belgium2023.csv', sep = ',', header = TRUE)
  life_table_year <- data.frame(age = life_table_BE$X,
                                LE = life_table_BE$EX)
  
  # load health-related Quality of Life
  qol_age_year <- popnormINDEX_BE()$mean
  
  # calculate discounted LE
  life_table_year$LE_disc <- get_HrQALE(life_table_year$LE,
                                       f_disc_rate_qaly = f_disc_rate_qaly,
                                       qol_age_year = NA)
  # calculate HrQoL adjusted LE
  life_table_year$HrQALE <- get_HrQALE(life_table_year$LE,
                                          f_disc_rate_qaly = 0,
                                           qol_age_year)
  # calculate discounted HrQoL adjusted LE
  life_table_year$HrQALE_disc <- get_HrQALE(life_table_year$LE,
                                                f_disc_rate_qaly = f_disc_rate_qaly,
                                                qol_age_year)
  life_table_year
  
  # switch from years to months
  # linear approximation for life_expectancy to age in MONTHS (for 100 months)
  life_table             <- get_LE_year2months(life_table_year$age, life_table_year$LE, max_age_months)
  life_table$LE_disc     <- get_LE_year2months(life_table_year$age, life_table_year$LE_disc, max_age_months)$LE
  life_table$HrQALE      <- get_LE_year2months(life_table_year$age, life_table_year$HrQALE, max_age_months)$LE
  life_table$HrQALE_disc <- get_LE_year2months(life_table_year$age, life_table_year$HrQALE_disc, max_age_months)$LE
  
  #return
  return(life_table)
}

