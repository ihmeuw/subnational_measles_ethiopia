

## clear environment
rm(list=ls())
rm(list = ls(all.names = TRUE))

## Set repo location and indicator group
user               <- Sys.info()['user']
core_repo          <- sprintf('FILEPATH/',user)
vax_repo           <- sprintf('FILEPATH/',user)
measles_repo       <- sprintf('FILEPATH/',user)
remote             <- 'origin'
branch             <- 'develop'
pullgit            <- FALSE

## sort some directory stuff
commondir      <- sprintf('FILEPATH')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

library(MASS)

library(rgeos)
library(raster)
library(INLA)
library(seegSDM)
library(seegMBG)
library(gbm)
library(foreign)
library(parallel)
library(doParallel)
library(grid)
library(gridExtra)
library(pacman)
library(data.table)
library(gtools)
library(glmnet)
library(ggplot2)
library(RMySQL)
library(plyr)
library(tictoc)
##library(dplyr)
library(magrittr)
library(data.table)
source(paste0('FILEPATH/misc_vaccine_functions.R'))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0('FILEPATH/setup.R'))
library(kernlab, lib.loc='FILEPATH')
library(data.table)
library(ggplot2)

iso3 <- 'ETH'

library(argparse)
# 

parser <- argparse::ArgumentParser()
parser$add_argument('--boot_i')
# parser$add_argument('--outer_number')
parser$add_argument('--run_date')


args <- parser$parse_args()
boot_i <- args$boot_i
run_date <- args$run_date




dir.create(paste0('FILEPATH/',run_date))

bootstrap_lookup <- matrix(NA, nrow=2120, ncol=79)
bootstrap_lookup_temp <- matrix(NA, nrow=371, ncol=79)
for(d in 1:79){
  message(d)
  sample_boot <- sample(1:371, size = 93, replace=F)  
  
  for(i in 1:371){
    if(i %in% sample_boot){
      bootstrap_lookup_temp[i,d] = 1
    }else{
      bootstrap_lookup_temp[i,d] = 0
    }
  }
}

bootstrap_lookup[1750:2120,] <- bootstrap_lookup_temp
save(bootstrap_lookup, file = paste0('FILEPATH/',run_date,'/bootstrap_',boot_i,'_bootstrap_lookup.RData') )

############################################################################################################
## load processed data
############################################################################################################
load(paste0('FILEPATH/', iso3,'_processed_births_by_adm2.RData'))
load(paste0('FILEPATH/', iso3,'_processed_population_by_epiweek_measles_model_age_bins_by_adm2.RData'))
load(paste0('FILEPATH/', iso3_,'_processed_mobility_theta99.RData'))
M <- M_theta99

load(paste0('FILEPATH/', iso3,'_processed_contact_patterns_array_district.RData'))
load(paste0('FILEPATH/', iso3,'_processed_population.RData'))
load(paste0('FILEPATH/', iso3,'_processed_cases.RData')) # agg nationally
load(paste0('FILEPATH/', iso3,'_processed_cases_loess.RData')) # subnational
load(paste0('FILEPATH/', iso3,'_processed_national_cases.RData'))
load(paste0('FILEPATH/', iso3,'_processed_coverage.RData'))
load(paste0('FILEPATH/', iso3,'_processed_contact_patterns.RData'))
load(paste0('FILEPATH/', iso3,'_processed_mortality.RData'))
load(paste0('FILEPATH/', iso3,'_processed_migration.RData'))
load('FILEPATH/W_M_composite.RData')

i_ticker <- rep(1:53, length(1980:2019))
i_ticker_1 <- i_ticker / 53
i_ticker_1 <- i_ticker_1[1:2120]

all_strat_cases = cases_age_bin0_4[1] + cases_age_bin5_9[1] + cases_age_bin10_14[1] + cases_age_bin15_19[1] + cases_age_bin20_24[1] + cases_age_bin25_29[1] +
  cases_age_bin30_34[1] + cases_age_bin35_39[1] + cases_age_bin40_44[1] + cases_age_bin45_49[1] + cases_age_bin50_54[1] + cases_age_bin55_59[1] +
  cases_age_bin60_64[1] + cases_age_bin65_69[1] + cases_age_bin70_74[1] + cases_age_bin75_79[1] + cases_age_bin80_84[1]

age_cov_prevalence=ri_sia_time_series_mcv1
age_cov_prevalence_counts=age_cov_prevalence * pop_array
age_cov_prevalence_counts[,107,] <- 0
age_cov_prevalence_counts_red <- array(NA, dim=c(2120, 24,79))
age_cov_prevalence_counts_red[,1:16,] = age_cov_prevalence_counts[1:2120,1:16,1:79]

age_cov2_prevalence=ri_sia_time_series_mcv2
age_cov2_prevalence_counts=age_cov2_prevalence * pop_array
age_cov2_prevalence_counts[,107,] <- 0
age_cov2_prevalence_counts_red <- array(NA, dim=c(2120, 24,79))
age_cov2_prevalence_counts_red[,1:16,] = age_cov2_prevalence_counts[1:2120,1:16,1:79]

pop_array_red <- array(NA, dim=c(2120, 24,79))
pop_array_red[,1:16,] = pop_array[1:2120,1:16,1:79]

for(i in 1:2120){
  for(d in 1:79){
    age_cov_prevalence_counts_red[i,17,d] = sum(age_cov_prevalence_counts[i,17:21, d])
    age_cov_prevalence_counts_red[i,18,d] = sum(age_cov_prevalence_counts[i,22:26, d])
    age_cov_prevalence_counts_red[i,19,d] = sum(age_cov_prevalence_counts[i,27:36, d])
    age_cov_prevalence_counts_red[i,20,d] = sum(age_cov_prevalence_counts[i,37:46, d])
    age_cov_prevalence_counts_red[i,21,d] = sum(age_cov_prevalence_counts[i,47:56, d])
    age_cov_prevalence_counts_red[i,22,d] = sum(age_cov_prevalence_counts[i,57:66, d])
    age_cov_prevalence_counts_red[i,23,d] = sum(age_cov_prevalence_counts[i,67:76, d])
    age_cov_prevalence_counts_red[i,24,d] = sum(age_cov_prevalence_counts[i,77:107, d])
    
    age_cov2_prevalence_counts_red[i,17,d] = sum(age_cov2_prevalence_counts[i,17:21, d])
    age_cov2_prevalence_counts_red[i,18,d] = sum(age_cov2_prevalence_counts[i,22:26, d])
    age_cov2_prevalence_counts_red[i,19,d] = sum(age_cov2_prevalence_counts[i,27:36, d])
    age_cov2_prevalence_counts_red[i,20,d] = sum(age_cov2_prevalence_counts[i,37:46, d])
    age_cov2_prevalence_counts_red[i,21,d] = sum(age_cov2_prevalence_counts[i,47:56, d])
    age_cov2_prevalence_counts_red[i,22,d] = sum(age_cov2_prevalence_counts[i,57:66, d])
    age_cov2_prevalence_counts_red[i,23,d] = sum(age_cov2_prevalence_counts[i,67:76, d])
    age_cov2_prevalence_counts_red[i,24,d] = sum(age_cov2_prevalence_counts[i,77:107, d])
    
    pop_array_red[i,17,d] = sum(pop_array[i,17:21, d])
    pop_array_red[i,18,d] = sum(pop_array[i,22:26, d])
    pop_array_red[i,19,d] = sum(pop_array[i,27:36, d])
    pop_array_red[i,20,d] = sum(pop_array[i,37:46, d])
    pop_array_red[i,21,d] = sum(pop_array[i,47:56, d])
    pop_array_red[i,22,d] = sum(pop_array[i,57:66, d])
    pop_array_red[i,23,d] = sum(pop_array[i,67:76, d])
    pop_array_red[i,24,d] = sum(pop_array[i,77:107, d])
  }
}

age_cov_prevalence_counts_vector <- as.vector(age_cov_prevalence_counts_red)
age_cov2_prevalence_counts_vector <- as.vector(age_cov2_prevalence_counts_red)
pop_districts <- matrix(NA, 2120, 79)

for(d in 1:79){
  pop_districts[,d] <- rowSums(pop_array_red[,,d])
}
pop_array_vector_district <- as.vector(pop_districts)
pop_array_vector <- as.vector(pop_array_red)

# first set -- under 1s
W[,1:16,,] <- W[,1:16,,] / 2
W_vector <- as.vector(W)

national_cases[6] <- (national_cases[5] + national_cases[7]) / 2
zero_matrix <- matrix(0, nrow = 24, ncol=79)

cases_age_bin0_4_matrix[cases_age_bin0_4_matrix < 0] <- 0
cases_age_bin5_9_matrix[cases_age_bin5_9_matrix < 0] <- 0
cases_age_bin10_14_matrix[cases_age_bin10_14_matrix < 0] <- 0
cases_age_bin15_24_matrix[cases_age_bin15_24_matrix < 0] <- 0
cases_age_bin25_34_matrix[cases_age_bin25_34_matrix < 0] <- 0
cases_age_bin35_44_matrix[cases_age_bin35_44_matrix < 0] <- 0
cases_age_bin45_54_matrix[cases_age_bin45_54_matrix < 0] <- 0
cases_age_bin55_64_matrix[cases_age_bin55_64_matrix < 0] <- 0
cases_age_bin65_plus_matrix[cases_age_bin65_plus_matrix < 0] <- 0

autocorrelation <- matrix(1, nrow=79,ncol=9)
########################################################################################################################
library(tmvtnorm)
library(Rcpp)
library(dfoptim)
library(arm)
library(tictoc)

setwd('FILEPATH')
sourceCpp("/FILEPATH/Rcpp_rolling_average_estimated_no_eff_subnat_only_down_adjust_regional_reporting_vax_eff_vary2_loess_seasonal_serial_interval_incidence_autocorrelation_nbinom5_single_reporting_bootstrap_district.cpp")
sourceCpp("FILEPATH/Rcpp_rolling_average_estimated_no_eff_subnat_only_down_adjust_regional_reporting_vax_eff_vary2_smoothed_seasonal_with_predict_serial_incidence2.cpp")
sourceCpp('FILEPATH/Rcpp_BCD_rho_fitting_bootstrap_district_single_reporting.cpp')

births = as.matrix(births_matrix)
mortality_matrix = as.matrix(mortality_matrix)
migration_matrix = as.matrix(migration_matrix)

##################################################################################################################################
#### set vaccine effectiveness paramters
logit_vax_eff <- 0

##################################################################################################################################
#### define internal functions

Likelihoodcpp_measles_transmission2_omp_short <- function(pars){
  
  liklihood_init <- Likelihoodcpp_measles_transmission2_omp(pop=pop, pop_array_vector=pop_array_vector,
                                                            age_cov_prevalence_counts_vector = age_cov_prevalence_counts_vector,
                                                            age_cov2_prevalence_counts_vector = age_cov2_prevalence_counts_vector,
                                                            W_vector = W_vector,
                                                            pop_array_vector_district=pop_array_vector_district,
                                                            
                                                            national_cases = national_cases, 
                                                            all_strat_cases = all_strat_cases,
                                                            
                                                            zero_matrix = zero_matrix,
                                                            
                                                            births=births_matrix,
                                                            
                                                            mortality_matrix = as.matrix(mortality_matrix), 
                                                            migration_matrix=as.matrix(migration_matrix),
                                                            
                                                            cases_age_bin0_4 = cases_age_bin0_4, cases_age_bin5_9 = cases_age_bin5_9,
                                                            cases_age_bin10_14 = cases_age_bin10_14, cases_age_bin15_19 = cases_age_bin15_19,
                                                            cases_age_bin20_24 = cases_age_bin20_24, cases_age_bin25_29 = cases_age_bin25_29,
                                                            cases_age_bin30_34 = cases_age_bin30_34, cases_age_bin35_39 = cases_age_bin35_39,
                                                            cases_age_bin40_44 = cases_age_bin40_44, cases_age_bin45_49 = cases_age_bin45_49,
                                                            cases_age_bin50_54 = cases_age_bin50_54, cases_age_bin55_59 = cases_age_bin55_59,
                                                            cases_age_bin60_64 = cases_age_bin60_64, cases_age_bin65_69 = cases_age_bin65_69,
                                                            cases_age_bin70_74 = cases_age_bin70_74, cases_age_bin75_79 = cases_age_bin75_79,
                                                            cases_age_bin80_84 = cases_age_bin80_84,
                                                            
                                                            cases_age_bin0_4_matrix = as.matrix(cases_age_bin0_4_matrix), 
                                                            cases_age_bin5_9_matrix = as.matrix(cases_age_bin5_9_matrix),
                                                            cases_age_bin10_14_matrix = as.matrix(cases_age_bin10_14_matrix), 
                                                            
                                                            cases_age_bin15_24_matrix = as.matrix(cases_age_bin15_24_matrix),
                                                            cases_age_bin25_34_matrix = as.matrix(cases_age_bin25_34_matrix), 
                                                            cases_age_bin35_44_matrix = as.matrix(cases_age_bin35_44_matrix), 
                                                            cases_age_bin45_54_matrix = as.matrix(cases_age_bin45_54_matrix), 
                                                            cases_age_bin55_64_matrix = as.matrix(cases_age_bin55_64_matrix), 
                                                            cases_age_bin65_plus_matrix = as.matrix(cases_age_bin65_plus_matrix), 
                                                            
                                                            i_ticker_1 = i_ticker_1,
                                                            M = M,
                                                            autocorrelation = as.matrix(autocorrelation),
                                                            
                                                            max_beta = pars[1], 
                                                            min_beta = pars[2], 
                                                            rho=invlogit(logit_rho_vals[i+1]),
                                                            logit_vax_efficacy = logit_vax_eff,
                                                            bootstrap_lookup = bootstrap_lookup)
  return(liklihood_init)

}


Likelihood_reporting_short <- function(pars){
  
  -1 * likelihood_reporting_fit(rho = invlogit(pars[1]), 
                                I_annual_by_district_by_week_0_4 = I_annual_by_district_by_week_0_4, 
                                I_annual_by_district_by_week_5_9 = I_annual_by_district_by_week_5_9, 
                                I_annual_by_district_by_week_10_14 = I_annual_by_district_by_week_10_14, 
                                I_annual_by_district_by_week_15_24 = I_annual_by_district_by_week_15_24, 
                                I_annual_by_district_by_week_25_34 = I_annual_by_district_by_week_15_24, 
                                I_annual_by_district_by_week_35_44 = I_annual_by_district_by_week_35_44, 
                                I_annual_by_district_by_week_45_54 = I_annual_by_district_by_week_45_54, 
                                I_annual_by_district_by_week_55_64 = I_annual_by_district_by_week_55_64, 
                                I_annual_by_district_by_week_65_plus = I_annual_by_district_by_week_65_plus, 
                                cases_age_bin0_4_matrix = as.matrix(cases_age_bin0_4_matrix), 
                                cases_age_bin5_9_matrix = as.matrix(cases_age_bin5_9_matrix), 
                                cases_age_bin10_14_matrix = as.matrix(cases_age_bin10_14_matrix),
                                cases_age_bin15_24_matrix = as.matrix(cases_age_bin15_24_matrix), 
                                cases_age_bin25_34_matrix = as.matrix(cases_age_bin25_34_matrix), 
                                cases_age_bin35_44_matrix = as.matrix(cases_age_bin35_44_matrix), 
                                cases_age_bin45_54_matrix = as.matrix(cases_age_bin45_54_matrix), 
                                cases_age_bin55_64_matrix = as.matrix(cases_age_bin55_64_matrix), 
                                cases_age_bin65_plus_matrix = as.matrix(cases_age_bin65_plus_matrix),
                                autocorrelation = as.matrix(autocorrelation),
                                bootstrap_lookup = bootstrap_lookup)

}

##################################################################################################################################
#### start BCD

beta_max_vals <- rep(NA, 11)
beta_min_vals <- rep(NA, 11)
logit_rho_vals <- rep(NA, 11)

beta_max_vals[1] <- 0.25
beta_min_vals[1] <- 0.12
logit_rho_vals[1] <- -5

beta_val_ll <- rep(NA, 11)
rho_val_ll <- rep(NA, 11)


for(i in 1:10){
  
  
  message("...... predicting and processing")
  
  prediction_measles_cases <- predict_Likelihoodcpp_measles_transmission2_omp(i_ticker_1 = i_ticker_1,
                                                                              national_cases = national_cases, all_strat_cases = all_strat_cases,
                                                                              
                                                                              pop=pop,
                                                                              
                                                                              cases_age_bin0_4 = cases_age_bin0_4, cases_age_bin5_9 = cases_age_bin5_9,
                                                                              cases_age_bin10_14 = cases_age_bin10_14, cases_age_bin15_19 = cases_age_bin15_19,
                                                                              cases_age_bin20_24 = cases_age_bin20_24, cases_age_bin25_29 = cases_age_bin25_29,
                                                                              cases_age_bin30_34 = cases_age_bin30_34, cases_age_bin35_39 = cases_age_bin35_39,
                                                                              cases_age_bin40_44 = cases_age_bin40_44, cases_age_bin45_49 = cases_age_bin45_49,
                                                                              cases_age_bin50_54 = cases_age_bin50_54, cases_age_bin55_59 = cases_age_bin55_59,
                                                                              cases_age_bin60_64 = cases_age_bin60_64, cases_age_bin65_69 = cases_age_bin65_69,
                                                                              cases_age_bin70_74 = cases_age_bin70_74, cases_age_bin75_79 = cases_age_bin75_79,
                                                                              cases_age_bin80_84 = cases_age_bin80_84,
                                                                              
                                                                              
                                                                              cases_age_bin0_4_matrix = as.matrix(cases_age_bin0_4_matrix), 
                                                                              cases_age_bin5_9_matrix = as.matrix(cases_age_bin5_9_matrix),
                                                                              cases_age_bin10_14_matrix = as.matrix(cases_age_bin10_14_matrix), 
                                                                              
                                                                              cases_age_bin15_24_matrix = as.matrix(cases_age_bin15_24_matrix),
                                                                              cases_age_bin25_34_matrix = as.matrix(cases_age_bin25_34_matrix), 
                                                                              cases_age_bin35_44_matrix = as.matrix(cases_age_bin35_44_matrix), 
                                                                              cases_age_bin45_54_matrix = as.matrix(cases_age_bin45_54_matrix), 
                                                                              cases_age_bin55_64_matrix = as.matrix(cases_age_bin55_64_matrix), 
                                                                              cases_age_bin65_plus_matrix = as.matrix(cases_age_bin65_plus_matrix), 
                                                                              
                                                                              zero_matrix = zero_matrix,
                                                                              
                                                                              births = as.matrix(births_matrix),
                                                                              
                                                                              age_cov_prevalence_counts_vector = age_cov_prevalence_counts_vector, age_cov2_prevalence_counts_vector = age_cov2_prevalence_counts_vector,
                                                                              
                                                                              pop_array_vector = pop_array_vector,
                                                                              pop_array_vector_district = pop_array_vector_district,
                                                                              
                                                                              mortality_matrix = as.matrix(mortality_matrix), migration_matrix = as.matrix(migration_matrix),
                                                                              W_vector = W_vector, M = M,
                                                                              
                                                                              max_beta = beta_max_vals[i], 
                                                                              min_beta = beta_min_vals[i], 
                                                                              logit_rho1 = logit_rho_vals[i], 
                                                                              logit_rho2 = logit_rho_vals[i], 
                                                                              logit_rho3 = logit_rho_vals[i], 
                                                                              logit_rho4 = logit_rho_vals[i], 
                                                                              logit_rho5 = logit_rho_vals[i], 
                                                                              logit_rho6 = logit_rho_vals[i], 
                                                                              logit_rho7 = logit_rho_vals[i], 
                                                                              logit_rho8 = logit_rho_vals[i], 
                                                                              logit_rho9 = logit_rho_vals[i], 
                                                                              logit_rho10 = logit_rho_vals[i], 
                                                                              logit_rho11 = logit_rho_vals[i], 
                                                                              
                                                                              logit_vax_efficacy = logit_vax_eff)
  
  prediction_measles_cases <- as.data.frame(prediction_measles_cases)
  
  I_predicted <- array(NA, dim=c(2120,24,79))
  
  for(w in 1750:2120){
    for(a in 1:24){
      for(d in 1:79){
        I_predicted[w,a,d] = prediction_measles_cases[w + 2120*(a-1) + 2120*24*(d-1),5];
      }
    }
  }
  #############################################################################################################
  I_final_counts <- I_predicted #* pop_array_red
  
  I_annual_by_district_by_week_0_4 <- matrix(0, 2120, 79)
  I_annual_by_district_by_week_5_9 <- matrix(0, 2120, 79)
  I_annual_by_district_by_week_10_14 <- matrix(0, 2120, 79)
  I_annual_by_district_by_week_15_24 <- matrix(0, 2120, 79)
  I_annual_by_district_by_week_25_34 <- matrix(0, 2120, 79)
  I_annual_by_district_by_week_35_44 <- matrix(0, 2120, 79)
  I_annual_by_district_by_week_45_54 <- matrix(0, 2120, 79)
  I_annual_by_district_by_week_55_64 <- matrix(0, 2120, 79)
  I_annual_by_district_by_week_65_plus <- matrix(0, 2120, 79)
  
  for (w in 1750:2120){
    for(d in 1:79){
      for(a in 1:16){
        I_annual_by_district_by_week_0_4[w,d] = I_annual_by_district_by_week_0_4[w, d] + I_final_counts[w, a, d]
      }
      I_annual_by_district_by_week_5_9[w,d] = I_annual_by_district_by_week_5_9[w, d] + I_final_counts[w, 17, d]
      I_annual_by_district_by_week_10_14[w,d] = I_annual_by_district_by_week_10_14[w, d] + I_final_counts[w, 18, d]
      I_annual_by_district_by_week_15_24[w,d] = I_annual_by_district_by_week_15_24[w, d] + I_final_counts[w, 19, d]
      I_annual_by_district_by_week_25_34[w,d] = I_annual_by_district_by_week_25_34[w, d] + I_final_counts[w, 20, d]
      I_annual_by_district_by_week_35_44[w,d] = I_annual_by_district_by_week_35_44[w, d] + I_final_counts[w, 21, d]
      I_annual_by_district_by_week_45_54[w,d] = I_annual_by_district_by_week_45_54[w, d] + I_final_counts[w, 22, d]
      I_annual_by_district_by_week_55_64[w,d] = I_annual_by_district_by_week_55_64[w, d] + I_final_counts[w, 23, d]
      I_annual_by_district_by_week_65_plus[w,d] = I_annual_by_district_by_week_65_plus[w, d] + I_final_counts[w, 24, d]
      
      if(I_annual_by_district_by_week_0_4[w,d] < 0.001){
        I_annual_by_district_by_week_0_4[w,d] = 0.001
      }
      if(I_annual_by_district_by_week_5_9[w,d] < 0.001){
        I_annual_by_district_by_week_5_9[w,d] = 0.001
      }
      if(I_annual_by_district_by_week_10_14[w,d] < 0.001){
        I_annual_by_district_by_week_10_14[w,d] = 0.001
      }
      if(I_annual_by_district_by_week_15_24[w,d] < 0.001){
        I_annual_by_district_by_week_15_24[w,d] = 0.001
      }
      if(I_annual_by_district_by_week_25_34[w,d] < 0.001){
        I_annual_by_district_by_week_25_34[w,d] = 0.001
      }
      if(I_annual_by_district_by_week_35_44[w,d] < 0.001){
        I_annual_by_district_by_week_35_44[w,d] = 0.001
      }
      if(I_annual_by_district_by_week_45_54[w,d] < 0.001){
        I_annual_by_district_by_week_45_54[w,d] = 0.001
      }
      if(I_annual_by_district_by_week_55_64[w,d] < 0.001){
        I_annual_by_district_by_week_55_64[w,d] = 0.001
      }
      if(I_annual_by_district_by_week_65_plus[w,d] < 0.001){
        I_annual_by_district_by_week_65_plus[w,d] = 0.001
      }
    }
  }
  
  message("......fitting rho likelihood")
  
  model_fit_object_new_rho <- nlm(Likelihood_reporting_short, c(logit_rho_vals[i]))
  
  logit_rho_vals[i+1] <- model_fit_object_new_rho$estimate
  
  rho_val_ll[i+1] <- -1 * model_fit_object_new_rho$minimum
  
  
  message("...Now starting iteration: ", i)
  message("......fitting beta likelihood")
  tic()
  model_fit_object_betas <- nmkb(fn= Likelihoodcpp_measles_transmission2_omp_short, 
                           par=c(beta_max_vals[i], beta_min_vals[i]), 
                           lower = c(0.1, 0), 
                           upper = c(1, 0.99),  
                           control = list(trace = T, maxfeval = 20000, maximize = T))
  
  toc()
  
  beta_max_vals[i+1] <- model_fit_object_betas$par[1]
  beta_min_vals[i+1] <- model_fit_object_betas$par[2]
  
  beta_val_ll[i+1] <- model_fit_object_betas$value
  
  if(i == 10){
    
    prediction_measles_cases <- predict_Likelihoodcpp_measles_transmission2_omp(i_ticker_1 = i_ticker_1,
                                                                                national_cases = national_cases, all_strat_cases = all_strat_cases,
                                                                                
                                                                                pop=pop,
                                                                                
                                                                                cases_age_bin0_4 = cases_age_bin0_4, cases_age_bin5_9 = cases_age_bin5_9,
                                                                                cases_age_bin10_14 = cases_age_bin10_14, cases_age_bin15_19 = cases_age_bin15_19,
                                                                                cases_age_bin20_24 = cases_age_bin20_24, cases_age_bin25_29 = cases_age_bin25_29,
                                                                                cases_age_bin30_34 = cases_age_bin30_34, cases_age_bin35_39 = cases_age_bin35_39,
                                                                                cases_age_bin40_44 = cases_age_bin40_44, cases_age_bin45_49 = cases_age_bin45_49,
                                                                                cases_age_bin50_54 = cases_age_bin50_54, cases_age_bin55_59 = cases_age_bin55_59,
                                                                                cases_age_bin60_64 = cases_age_bin60_64, cases_age_bin65_69 = cases_age_bin65_69,
                                                                                cases_age_bin70_74 = cases_age_bin70_74, cases_age_bin75_79 = cases_age_bin75_79,
                                                                                cases_age_bin80_84 = cases_age_bin80_84,
                                                                                
                                                                                
                                                                                cases_age_bin0_4_matrix = as.matrix(cases_age_bin0_4_matrix), 
                                                                                cases_age_bin5_9_matrix = as.matrix(cases_age_bin5_9_matrix),
                                                                                cases_age_bin10_14_matrix = as.matrix(cases_age_bin10_14_matrix), 
                                                                                
                                                                                cases_age_bin15_24_matrix = as.matrix(cases_age_bin15_24_matrix),
                                                                                cases_age_bin25_34_matrix = as.matrix(cases_age_bin25_34_matrix), 
                                                                                cases_age_bin35_44_matrix = as.matrix(cases_age_bin35_44_matrix), 
                                                                                cases_age_bin45_54_matrix = as.matrix(cases_age_bin45_54_matrix), 
                                                                                cases_age_bin55_64_matrix = as.matrix(cases_age_bin55_64_matrix), 
                                                                                cases_age_bin65_plus_matrix = as.matrix(cases_age_bin65_plus_matrix), 
                                                                                
                                                                                zero_matrix = zero_matrix,
                                                                                
                                                                                births = as.matrix(births_matrix),
                                                                                
                                                                                age_cov_prevalence_counts_vector = age_cov_prevalence_counts_vector, age_cov2_prevalence_counts_vector = age_cov2_prevalence_counts_vector,
                                                                                
                                                                                pop_array_vector = pop_array_vector,
                                                                                pop_array_vector_district = pop_array_vector_district,
                                                                                
                                                                                mortality_matrix = as.matrix(mortality_matrix), migration_matrix = as.matrix(migration_matrix),
                                                                                W_vector = W_vector, M = M,
                                                                                
                                                                                max_beta = beta_max_vals[i+1], 
                                                                                min_beta = beta_min_vals[i+1], 
                                                                                logit_rho1 = logit_rho_vals[i+1], 
                                                                                logit_rho2 = logit_rho_vals[i+1], 
                                                                                logit_rho3 = logit_rho_vals[i+1], 
                                                                                logit_rho4 = logit_rho_vals[i+1], 
                                                                                logit_rho5 = logit_rho_vals[i+1], 
                                                                                logit_rho6 = logit_rho_vals[i+1], 
                                                                                logit_rho7 = logit_rho_vals[i+1], 
                                                                                logit_rho8 = logit_rho_vals[i+1], 
                                                                                logit_rho9 = logit_rho_vals[i+1], 
                                                                                logit_rho10 = logit_rho_vals[i+1], 
                                                                                logit_rho11 = logit_rho_vals[i+1], 
                                                                                
                                                                                logit_vax_efficacy = logit_vax_eff)
    
    prediction_measles_cases <- as.data.frame(prediction_measles_cases)
    
    I_predicted <- array(NA, dim=c(2120,24,79))
    
    for(w in 1750:2120){
      for(a in 1:24){
        for(d in 1:79){
          I_predicted[w,a,d] = prediction_measles_cases[w + 2120*(a-1) + 2120*24*(d-1),5];
        }
      }
    }
  }
  
  save(beta_max_vals, beta_min_vals, logit_rho_vals, beta_val_ll, rho_val_ll, I_predicted, file = paste0('FILEPATH/',run_date,'/bootstrap_',boot_i,'_BCD_iteration_', i, '.RData') )
  
  
}
          





