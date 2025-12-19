############################################################################ #
# This file is part of the RSV modelling project.
# 
# => GET LIFE TABLE
#
#  Copyright 2025, CHERMID, UNIVERSITY OF ANTWERP
############################################################################ #


# Get the life table with monthly life expectancy, based on mortality rates from the UN 
# f_country         country name
# f_year            number => is converted into string with 5 year interval, such as 2015-2020
# f_disc_rate_qaly  discount rate for QALY

# f_country_iso='GBR'
# f_year=2020
# f_disc_rate_qaly=0.04
get_life_table <- function(f_country_iso,f_year,f_disc_rate_qaly,max_age_months = 100)
{

  # get period
  period <- get_year_category(f_year)[1]
  
  # convert the country iso3 code into the country name (corresponding the wpp2019 package)
  f_country <- get_UN_country_name(f_country_iso)
  
  # load population data
  data(mxM, package = "wpp2019")
  data(mxF, package = "wpp2019")
  data(sexRatio, package = "wpp2019")

  if(!f_country %in% sexRatio$name){
    print(paste("Unknown country: ", f_country))
    print("Make sure that the country name starts with a capital letter")
    print("All country names are provided by the function 'get_country_list()'")
    return(NA)
  } 
  
  # get column for year interval
  sel_col     <- names(mxM) == sub("_","-",period)
  
  # age: A character string representing an age interval (given by the starting age of the interval).
  # note: Data for ages 85-100+ are not the official UN data. While the published UN mortality datasets 
  # contain data only up to 85+, data for ages 85-100+ in this dataset were derived from UN published 
  # life table quantities.
  wpp_data <- data.frame(age_group_min = mxM$age[mxM$name %in% f_country],
                         mxM = mxM[mxM$name %in% f_country,sel_col],
                         mxF = mxF[mxF$name %in% f_country,sel_col],
                         sexRatio = sexRatio[sexRatio$name %in% f_country,sel_col])
  
  # Adjust ages in years
  wpp_data$age_group_mean <- wpp_data$age_group_min + diff(c(wpp_data$age_group_min, max(wpp_data$age_group_min+5))) / 2
  
  # Account for gender balance
  # Estimates and projections of the sex ratio at birth derived as the number of males divided by the number of females
  wpp_data$prop_female <- 1/(1 + wpp_data$sexRatio) 
  wpp_data$prop_male   <- 1 - wpp_data$prop_female
  
  # Get population nMx: mortality rate for 5-year age groups
  wpp_data$nMx <- wpp_data$mxM * wpp_data$prop_male + wpp_data$mxF * wpp_data$prop_female
  
  # Approximation under linear model: get mortality rate by age in years
  life_table_year<- approx(wpp_data$age_group_min,wpp_data$nMx,xout=seq(0,100,1),rule = 2:1) 
  names(life_table_year) <- c('age','nMx')
  life_table_year<- data.frame(life_table_year) # store as data.frame

  # complete life table in years
  life_table_year<- create_full_life_table(life_table_year)
  
  # calculate discounted LE
  life_table_year$LE_disc <- get_HrQALE(life_table_year$LE,
                                        f_disc_rate_qaly = f_disc_rate_qaly,
                                        qol_age_year = NA)
  
  zz <- get_discounted_life_years(life_table_year$LE, f_disc_rate_qaly)
  head(zz)
  head(life_table_year$LE_disc)
  # load health-related Quality of Life
  if(f_country_iso != 'BEL') {
    warning("Use Belgian data on age-specific health-releated Quality of Life")
  }
  qol_age_year <- popnormINDEX_BE()$mean
  
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
  
  return(life_table)
}

get_country_list <- function()
{
  data(sexRatio, package = "wpp2019", envir = environment())
  return(levels(sexRatio$name))
}

get_UN_country_name <- function(country_iso3){
  
  # generated with create_UN_country_database
  # computational expensive procedure! => best to store local copy
  UNcountries <- readRDS('input/dictionary_UN_countrynames.rds')
  
  # look for the given iso3 code in the country database
  flag <- country_iso3 == UNcountries$iso3 
  
  # return the country name
  return(UNcountries$name[flag])
}

get_country_iso3 <- function(country_name){
  
  # get UN country name - ISO database 
  # generated with create_UN_country_database
  # computational expensive procedure! => best to store local copy
  UNcountries <- readRDS('input/dictionary_UN_countrynames.rds')
  
  # look for the given iso3 code in the country database
  flag <- country_name == UNcountries$name 
  
  # return the country name
  return(UNcountries$iso3[flag])
}

# computational expensive procedure! => best to store local copy
create_UN_country_database <- function(){
  
  # file name for the [country name - iso3 code] database
  dictionary_UN_countrynames_file <- file.path('output','dictionary_UN_countrynames.rds')
  
  smd_print('Create UN country database')
  data("UNlocations", package = 'wpp2019') # does not contain iso code
  UNcountries      <- data.frame(name = UNlocations$name) # vector with all country names
  UNcountries$iso3 <- countrycode(UNcountries$name,'country.name', 'iso3c', warn = FALSE) # add iso codes
  UNcountries$iso3[is.na(UNcountries$iso3)] <- 0
  
  saveRDS(UNcountries, file=dictionary_UN_countrynames_file)
  smd_print('UN country database completed')
}



## FUNCTION TO COMPLETE LIFE TABLE BASED ON MORTALITY RATE
#f_life_table=life_table_year
create_full_life_table <- function(f_life_table)
{

  # Copy life table
  # age - years
  # nMx - age-specific death rate between ages x and x+n
  life_table_out <- f_life_table
  
  # Convert rate to probability
  # nqx - probability of dying between ages x and x+n
  life_table_out$nqx <- 1 - exp(-life_table_out$nMx)
  
  # Calculate:
  # - lx - number of people left alive at age x
  # - ndx - number of people dying between ages x and x+n
  # - nLx - person-years lived between ages x and x+n  
  life_table_out$lx  <- NA # init columns, otherwise we cannot use row indices
  life_table_out$ndx <- NA # init columns, otherwise we cannot use row indices
  life_table_out$nLx <- NA # init columns, otherwise we cannot use row indices
  for(i_age in 1:dim(life_table_out)[1]){
     if(i_age == 1) {
      # lx starts with 100 000 and every year it is reduced by the number of death in ndx
      life_table_out$lx[1] <- 100000 
    } else {
      life_table_out$lx[i_age] <- life_table_out$lx[i_age-1] - life_table_out$ndx[i_age-1]  
    }
    life_table_out$ndx[i_age] <- life_table_out$lx[i_age] * life_table_out$nqx[i_age]
    life_table_out$nLx[i_age] <- life_table_out$ndx[i_age] / life_table_out$nMx[i_age] 
  }
  
  # calculate Tx - person-years lived above age x
  i_age <- 1
  num_ages <- dim(life_table_out)[1]
  for(i_age in 1:num_ages){
    life_table_out$Tx[i_age] <- sum(life_table_out$nLx[i_age:num_ages])
  }
  
  # Calculate life expectation of life at age x  (by Ex=Tx/Lx)
  life_table_out$LE <- life_table_out$Tx / life_table_out$lx
  
  return(life_table_out)
}

get_year_category <- function(f_year)
{
  year_cutoff <- (0:10)*5 + 2000
  year_index  <- max(which(f_year >= year_cutoff))
  period      <- paste(year_cutoff[c(year_index, year_index+1)], collapse = '_')
  
  return(c(period,year_cutoff[c(year_index)]))
}

# function to get discounted life years (element-wise operation for vectors)
get_discounted_life_years <- function(num_years_vector, disc_rate){
  
  foreach(num_years = num_years_vector,
          .combine='c')%do%{
    
            if(num_years<1){
              return(num_years)
            } else{
              vector_years <- 1:floor(num_years)
              num_years_disc_floor <- sum(1 / ((1 + disc_rate)^vector_years))
              
              remainder_years <- ceiling(num_years)
              num_years_disc_remaining <- sum(1 / ((1 + disc_rate)^remainder_years)) * (num_years - floor(num_years))
              
              # return
              return(num_years_disc_floor + num_years_disc_remaining)               
            }

  } -> num_years_disc 
  
  # return result
  return(num_years_disc)
}
