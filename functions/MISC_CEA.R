############################################################################ #
# This file is part of the RSV modelling project.
# 
# => MISCELLANEOUS HELP FUNCTION FOR COST-EFFECTIVENESS ANALYSES
#
#  Copyright 2025, CHERMID, UNIVERSITY OF ANTWERP
############################################################################ #


# calculate the net benefit, given a price per DALY 'wtp_level'
get_net_benefit <- function(daly_averted, incr_cost_disc, wtp_level)
{
  NB <- wtp_level * daly_averted - incr_cost_disc
  return(NB)
}


# calculate the EVPPI for continuous uncertain parameters
evppi_gam <- function(NB, model_parameter_values){
  
  # note: current should have NB == 0
  D_opt  <- which(!colSums(NB) == 0)
  D <- ncol(NB)
  N <- nrow(NB)
  
  g.hat_new <- matrix(0,nrow=N,ncol=D)
  for(d in D_opt)
  {
    model <- gam(NB[,d] ~ model_parameter_values) 
    g.hat_new[,d] <- model$fitted
  }  
  #plot(model_parameter_values,NB[,28])
  #points(model_parameter_values,g.hat_new[,28],pch=20)
  #points(model_parameter_values,g.hat_new[,20],col='green')
  #points(model_parameter_values,g.hat_new[,22],col='blue')
  
  perfect.info  <- mean(apply(g.hat_new,1,max))
  baseline      <- max(colSums(g.hat_new)/N)
  partial.evpi  <- round(perfect.info - baseline, digits=4) ## estimate EVPPI 
  
  return(partial.evpi)
}

# get EVPPI for discrete uncertain parameters
evppi_gam_discr <- function(NB, model_parameter_values){
  
  # note: current should have NB == 0
  D_opt  <- which(!colSums(NB) == 0)
  D <- ncol(NB)
  N <- nrow(NB)
  
  g.hat_new <- matrix(0,nrow=N,ncol=D)
  for(d in D_opt)
  {
    model <- gam(NB[,d] ~ model_parameter_values) 
    g.hat_new[,d] <- model$fitted
  }   
  
  perfect.info  <- mean(apply(g.hat_new,1,max))
  baseline      <- max(colSums(g.hat_new)/N)
  partial.evpi  <- round(perfect.info - baseline, digits=4) ## estimate EVPPI 
  
  return(partial.evpi)
}

# get EVPPI for age-specific uncertain parameters with 6 age classes
evppi_gam_agedist <- function(NB, model_parameter_values){

  # note: current should have NB == 0
  D_opt  <- which(!colSums(NB) == 0)
  D <- ncol(NB)
  N <- nrow(NB)
  
  g.hat_new <- matrix(0,nrow=N,ncol=D)
  for(d in D_opt)
  {
  model <- gam(NB[,d] ~ model_parameter_values[,1] + model_parameter_values[,2] + model_parameter_values[,3]
               + model_parameter_values[,4] + model_parameter_values[,5] + model_parameter_values[,6]
              ) 
   g.hat_new[,d] <- model$fitted
  }   
  
  perfect.info  <- mean(apply(g.hat_new,1,max))
  baseline      <- max(colSums(g.hat_new)/N)
  partial.evpi  <- round(perfect.info - baseline, digits=4) ## estimate EVPPI 
  
  return(partial.evpi)
}
