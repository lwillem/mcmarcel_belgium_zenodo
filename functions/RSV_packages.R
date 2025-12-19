#############################################################################
# This file is part of the RSV modelling project.
# 
# => LOAD REQUIRED R PACKAGES AND INSTALL THEM IF NOT PRESENT
#
#  Copyright 2025, CHERMID, UNIVERSITY OF ANTWERP
#############################################################################

# list all required CRAN packages for this project
#
# wpp2019       to obtain UN population data 
# countrycode   to convert country names into iso3 codes
# openxlsx      general features for xlsx files (configuration files)
# abind         to combine multidimensional arrays
# scales        to use transparent colours in the CE planes 
# ggplot2       to make (fancy) plots
# dplyr         to aggregate simulation results (without removing entire rows containing an NA)
# RColorBrewer  to increase the number of colours for the EVPPI plot
# foreach       to run foreach loops 
# dampack       to calculate ICERS and identify dominated strategies

all_packages <- c('wpp2019', 'countrycode', 'openxlsx', 'abind', 'scales', 
                  'ggplot2', 'dplyr', 'RColorBrewer', 'foreach',
                  'dampack') 

# load all packages (and install them if not present yet)
for (i_package in all_packages) {
  if (!i_package %in% rownames(installed.packages())) {
    install.packages(i_package)
  }
  library(i_package, character.only = TRUE, quietly = TRUE, 
          warn.conflicts = FALSE, verbose = FALSE)
}


