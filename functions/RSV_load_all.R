#############################################################################
# This file is part of the RSV modelling project.
# 
# => SCRIPT TO LOAD ALL PACKAGES AND FUNCTIONS OF MCMARCEL
#
#  Copyright 2025, CHERMID, UNIVERSITY OF ANTWERP
#############################################################################

# make sure the time and date are defined in ENG
Sys.setlocale(category = "LC_ALL",locale = "en_US.UTF-8")

# load packages
source('functions/RSV_packages.R')

# load RSV-related functions
source('functions/RSV_config.R')
source('functions/RSV_get_burden.R')
source('functions/RSV_get_incidence.R') 
source('functions/RSV_get_life_table.R') 
source('functions/RSV_get_protection.R')
source('functions/RSV_plot_CEA_results.R')
source('functions/RSV_QALE_BEL.R') # specific for Belgium 
source('functions/RSV_write_summary_tables.R')

