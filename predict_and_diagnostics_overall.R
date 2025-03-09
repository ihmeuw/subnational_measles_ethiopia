
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
# 
dir.create(paste0('FILEPATH/',run_date,'/'))

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

cases_age_bin0_4_matrix = as.matrix(cases_age_bin0_4_matrix)
cases_age_bin5_9_matrix = as.matrix(cases_age_bin5_9_matrix)
cases_age_bin10_14_matrix = as.matrix(cases_age_bin10_14_matrix)

cases_age_bin15_24_matrix = as.matrix(cases_age_bin15_19_matrix) + as.matrix(cases_age_bin20_24_matrix)
cases_age_bin15_24_matrix = as.matrix(cases_age_bin15_24_matrix)

cases_age_bin25_34_matrix = as.matrix(cases_age_bin25_29_matrix) + as.matrix(cases_age_bin30_34_matrix)
cases_age_bin25_34_matrix = as.matrix(cases_age_bin25_34_matrix)

cases_age_bin35_44_matrix = as.matrix(cases_age_bin35_39_matrix) + as.matrix(cases_age_bin40_44_matrix)
cases_age_bin35_44_matrix = as.matrix(cases_age_bin35_44_matrix)

cases_age_bin45_54_matrix = as.matrix(cases_age_bin45_49_matrix) + as.matrix(cases_age_bin50_54_matrix)
cases_age_bin45_54_matrix = as.matrix(cases_age_bin45_54_matrix)

cases_age_bin55_64_matrix = as.matrix(cases_age_bin55_59_matrix) + as.matrix(cases_age_bin60_64_matrix)
cases_age_bin55_64_matrix = as.matrix(cases_age_bin55_64_matrix)

cases_age_bin65_plus_matrix = as.matrix(cases_age_bin65_69_matrix) + as.matrix(cases_age_bin70_74_matrix) +
  as.matrix(cases_age_bin75_79_matrix) + as.matrix(cases_age_bin80_84_matrix)
cases_age_bin65_plus_matrix = as.matrix(cases_age_bin65_plus_matrix)


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
cases_age_bin0_4_matrix[cases_age_bin5_9_matrix < 0] <- 0
cases_age_bin10_14_matrix[cases_age_bin10_14_matrix < 0] <- 0
cases_age_bin15_24_matrix[cases_age_bin15_24_matrix < 0] <- 0
cases_age_bin25_34_matrix[cases_age_bin25_34_matrix < 0] <- 0
cases_age_bin35_44_matrix[cases_age_bin35_44_matrix < 0] <- 0
cases_age_bin45_54_matrix[cases_age_bin45_54_matrix < 0] <- 0
cases_age_bin55_64_matrix[cases_age_bin55_64_matrix < 0] <- 0
cases_age_bin65_plus_matrix[cases_age_bin65_plus_matrix < 0] <- 0

########################################################################################################################
########################################################################################################################

library(tmvtnorm)

births = as.matrix(births_matrix)
mortality_matrix = as.matrix(mortality_matrix)
migration_matrix = as.matrix(migration_matrix)

########################################################################################################################
########################################################################################################################
library(Rcpp)
#.libPaths("~/myRlib")
setwd('FILEPATH')
sourceCpp("FILEPATH/Rcpp_rolling_average_estimated_no_eff_subnat_only_down_adjust_regional_reporting_vax_eff_vary2_smoothed_seasonal_with_predict_serial_incidence2.cpp")


###################################################################################################
##################################################################################################

run_date <- 'RUNDATE'
load(paste0('FILEPATH/',run_date,'/model_fit_object.RData'))

#val_i <- 1471


max_beta_vec <- rep(NA, 100)
min_beta_vec <- rep(NA, 100)
logit_rho1_vec <- rep(NA, 100)
logit_rho2_vec <- rep(NA, 100)
logit_rho3_vec <- rep(NA, 100)
logit_rho4_vec <- rep(NA, 100)
logit_rho5_vec <- rep(NA, 100)
logit_rho6_vec <- rep(NA, 100)
logit_rho7_vec <- rep(NA, 100)
logit_rho8_vec <- rep(NA, 100)
logit_rho9_vec <- rep(NA, 100)
logit_rho10_vec <- rep(NA, 100)
logit_rho11_vec <- rep(NA, 100)




for(i in 1:100){
  
  
  load(paste0('FILEPATH/',run_date,'/bootstrap_',i,'_BCD_iteration_10.RData'))
  
  max_beta_vec[i] <- beta_max_vals[11]
  min_beta_vec[i] <- beta_min_vals[11]
  logit_rho1_vec[i] <- logit_rho1_vals[11]
  logit_rho2_vec[i] <- logit_rho2_vals[11]
  logit_rho3_vec[i] <- logit_rho3_vals[11]
  logit_rho4_vec[i] <- logit_rho4_vals[11]
  logit_rho5_vec[i] <- logit_rho5_vals[11]
  logit_rho6_vec[i] <- logit_rho6_vals[11]
  logit_rho7_vec[i] <- logit_rho7_vals[11]
  logit_rho8_vec[i] <- logit_rho8_vals[11]
  logit_rho9_vec[i] <- logit_rho9_vals[11]
  logit_rho10_vec[i] <- logit_rho10_vals[11]
  logit_rho11_vec[i] <- logit_rho11_vals[11]
  
  
  
}

max_beta_median <- median(max_beta_vec) #model_fit_object$par[1]#median(out2[2000:(i), 2], na.rm=T)
min_beta_median <- median(min_beta_vec) #median(out2[2000:(i), 3], na.rm=T)
logit_rho1_median <- median(logit_rho1_vec) #median(out2[2000:(i),4], na.rm=T)
logit_rho2_median <- median(logit_rho2_vec) #median(out2[2000:(i),5], na.rm=T)
logit_rho3_median <- median(logit_rho3_vec) #median(out2[2000:(i),6], na.rm=T)
logit_rho4_median <- median(logit_rho4_vec) #median(out2[2000:(i),7], na.rm=T)
logit_rho5_median <- median(logit_rho5_vec) #median(out2[2000:(i),8], na.rm=T)
logit_rho6_median <- median(logit_rho6_vec) #median(out2[2000:(i),9], na.rm=T)
logit_rho7_median <- median(logit_rho7_vec) #median(out2[2000:(i),10], na.rm=T)
logit_rho8_median <- median(logit_rho8_vec) #median(out2[2000:(i),11], na.rm=T)
logit_rho9_median <- median(logit_rho9_vec) #median(out2[2000:(i),12], na.rm=T)
logit_rho10_median <- median(logit_rho10_vec) #median(out2[2000:(i),13], na.rm=T)
logit_rho11_median <- median(logit_rho11_vec) #median(out2[2000:(i),14], na.rm=T)
logit_rho_t_median <- 9999 #median(out2[500:(i*10),14], na.rm=T)
logit_vax_efficacy_median <- 0 #model_fit_object$par[14]#model_fit_object$par[1]#median(out2[2000:(i),15], na.rm=T)
alpha_median <- 0.99#median(out2[500:(i*10),7], na.rm=T)
theta_median <- 0.99#median(out2[500:(i*10),8], na.rm=T)

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
                                                                            
                                                                            max_beta = max_beta_median, 
                                                                            min_beta = min_beta_median, 
                                                                            logit_rho1 = logit_rho1_median, 
                                                                            logit_rho2 = logit_rho2_median, 
                                                                            logit_rho3 = logit_rho3_median, 
                                                                            logit_rho4 = logit_rho4_median, 
                                                                            logit_rho5 = logit_rho5_median, 
                                                                            logit_rho6 = logit_rho6_median, 
                                                                            logit_rho7 = logit_rho7_median, 
                                                                            logit_rho8 = logit_rho8_median, 
                                                                            logit_rho9 = logit_rho9_median, 
                                                                            logit_rho10 = logit_rho10_median, 
                                                                            logit_rho11 = logit_rho11_median, 
                                                                            
                                                                            logit_vax_efficacy = logit_vax_efficacy_median)


prediction_measles_cases <- as.data.frame(prediction_measles_cases)

M_predicted <- array(NA, dim=c(2120,24,79))
S_predicted <- array(NA, dim=c(2120,24,79))
prevalence_predicted <- array(NA, dim=c(2120,24,79))
R_predicted <- array(NA, dim=c(2120,24,79))
I_predicted <- array(NA, dim=c(2120,24,79))

foi_estimated <- array(NA, dim=c(2120,24,79))

for(i in 1:2120){
  message(i)
  for(a in 1:24){
    for(d in 1:79){
      M_predicted[i,a,d] = prediction_measles_cases[i + 2120*(a-1) + 2120*24*(d-1),1];
      S_predicted[i,a,d] = prediction_measles_cases[i + 2120*(a-1) + 2120*24*(d-1),2];
      prevalence_predicted[i,a,d] = prediction_measles_cases[i + 2120*(a-1) + 2120*24*(d-1),3];
      R_predicted[i,a,d] = prediction_measles_cases[i + 2120*(a-1) + 2120*24*(d-1),4];
      I_predicted[i,a,d] = prediction_measles_cases[i + 2120*(a-1) + 2120*24*(d-1),5];
      
      foi_estimated[i,a,d] = prediction_measles_cases[i + 2120*(a-1) + 2120*24*(d-1),6];
    }
  }
}
mod_name <-  run_date #<- '2023_04_07_ETH_loess_serial_interval_incidence_vax_eff075_HIGHESTdraw1471_NEW_justI3'


#############################################################################################################
M_final_counts <- M_predicted #* pop_array_red
S_final_counts <- S_predicted #* pop_array_red
I_final_counts <- I_predicted #* pop_array_red
R_final_counts <- R_predicted #* pop_array_red
prevalence_final_counts <- prevalence_predicted
# 

dir.create(paste0("FILEPATH/",run_date,"/predict/"), recursive = T)

save(M_final_counts, S_final_counts, I_final_counts, R_final_counts, prevalence_final_counts, file = paste0("FILEPATH/",run_date,"/predict/count_compartments.RData"))

# R_vax_final_counts <- R_vax_predicted * pop_array_red
# 
# R_vax2_final_counts <- R_vax2_predicted * pop_array_red
library(arm)
#############################################################################################################
year_list <- 1980:2019
I_annual <- rep(0, 40)
adj_cases_annual <- rep(NA, 40)
for(y in 2:40){
  I_annual[y] = 0
  for (i in ((y-1)*53 + 1):(y*53)){
    for(a in 1:24){
      for(d in 1:79){
        I_annual[y] = I_annual[y] + I_final_counts[i, a, d]
      }
    }
  }
  adj_cases_annual[y] = national_cases[y] / invlogit(-6)
}

I_annual_for_sim_data <- (I_annual)  
results <- data.table(I_annual, year_list)

### make plots
I_summary_mean <- results
I_summary_mean$group <- "TSIR"
colnames(I_summary_mean) <- c("measles", "time", "group")

I_summary_cases <- data.table(cbind(year_list, adj_cases_annual))
I_summary_cases$group <- "(adjusted) case notifications -- fake reporting"
colnames(I_summary_cases) <- c("time", "measles", "group")
# 
# I_summary_cases2 <- data.table(cbind(year_list, national_cases / 0.0025))
# I_summary_cases2$group <- "(FAKE adjusted) case notifications"
# colnames(I_summary_cases2) <- c("time", "measles", "group")

I_summary_cases_raw <- data.table(cbind(year_list, national_cases))
I_summary_cases_raw$group <- "(raw) case notifications"
colnames(I_summary_cases_raw) <- c("time", "measles", "group")

I_summary <- rbind(I_summary_mean, I_summary_cases, I_summary_cases_raw)
I_summary_all <- rbind(I_summary)

### make plot
p2 <- ggplot(I_summary_all, aes(x=time, y=measles, group=group, color=group)) +
  geom_line(aes(linetype=group))+ theme_bw() + xlab("Year") + ylab("Cases") +
  geom_point(aes(shape=group)) + ggtitle("Cases")

run_date <- mod_name

dir.create(paste0("FILEPATH/",run_date,"/predict/"), recursive=T)
png(file = paste0("FILEPATH/",run_date,"/predict/",iso3,"p1_case_time_series.png"),
    width = 10,
    height = 7,
    units = "in",
    res = 400,
    type = "cairo")
print(p2)
dev.off()

##### weekly at end plots
I_summary_mean <- results
I_summary_mean$group <- "TSIR"
colnames(I_summary_mean) <- c("measles", "time", "group")

subnat_cases_all_age_yr <- rep(0, 7)
subnat_cases_all_age_yr_raw <- rep(0, 7)

weekly_cases <- cases_age_bin0_4_matrix + cases_age_bin5_9_matrix + cases_age_bin10_14_matrix + 
  cases_age_bin15_24_matrix + 
  cases_age_bin25_34_matrix + 
  cases_age_bin35_44_matrix + 
  cases_age_bin45_54_matrix + 
  cases_age_bin55_64_matrix + 
  cases_age_bin65_plus_matrix 


weekly_cases_rows_raw <- rowSums(weekly_cases)
weekly_cases[,1] <- weekly_cases[,1] / invlogit(logit_rho1_median )
weekly_cases[,c(2, 12, 20, 28, 34)] <- weekly_cases[,c(2, 12, 20, 28, 34)] / invlogit(logit_rho2_median)
weekly_cases[,c(3, 13, 21, 29, 35, 40, 44, 48, 52, 56, 59, 62)] <- weekly_cases[,c(3, 13, 21, 29, 35, 40, 44, 48, 52, 56, 59, 62)] / invlogit(logit_rho3_median)
weekly_cases[,c(4, 14, 22)] <- weekly_cases[,c(4, 14, 22)] / invlogit(logit_rho4_median)
weekly_cases[,c(5)] <- weekly_cases[,c(5)] / invlogit(logit_rho5_median)
weekly_cases[,c(6, 15, 23)] <- weekly_cases[,c(6, 15, 23)] / invlogit(logit_rho6_median)
weekly_cases[,c(7)] <- weekly_cases[,c(7)] / invlogit(logit_rho7_median)
weekly_cases[,c(8, 16, 24, 30, 36, 41, 45, 49, 53, 57, 60, 63, 65, 67, 69, 71, 73)] <- weekly_cases[,c(8, 16, 24, 30, 36, 41, 45, 49, 53, 57, 60, 63, 65, 67, 69, 71, 73)] / invlogit(logit_rho8_median)
weekly_cases[,c(9, 17, 25, 31, 37, 42, 46, 50, 54)] <- weekly_cases[,c(9, 17, 25, 31, 37, 42, 46, 50, 54)] / invlogit(logit_rho9_median)
weekly_cases[,c(10, 18, 26, 32, 38, 43, 47, 51, 55, 58, 61, 64, 66, 68, 70, 72, 74, 75, 76, 77, 78, 79)] <- weekly_cases[,c(10, 18, 26, 32, 38, 43, 47, 51, 55, 58, 61, 64, 66, 68, 70, 72, 74, 75, 76, 77, 78, 79)] / invlogit(logit_rho10_median)
weekly_cases[,c(11, 19, 27, 33, 39)] <- weekly_cases[,c(11, 19, 27, 33, 39)] / invlogit(logit_rho11_median)


weekly_cases_rows <- rowSums(weekly_cases)


for(i in 1:7){
  for(w in (((i-1)*53)+1) :(i*53) ){
    subnat_cases_all_age_yr[i] = subnat_cases_all_age_yr[i] + weekly_cases_rows[w+1749]
    subnat_cases_all_age_yr_raw[i] = subnat_cases_all_age_yr_raw[i] + weekly_cases_rows_raw[w+1749]
    
  }
}

national_cases_subnat_end <- c(rep(NA,33), subnat_cases_all_age_yr)
national_cases_subnat_end_raw <- c(rep(NA,33), subnat_cases_all_age_yr_raw)

I_summary_cases <- data.table(cbind(year_list, national_cases_subnat_end ))
I_summary_cases$group <- "(adjusted) case notifications"
colnames(I_summary_cases) <- c("time", "measles", "group")

I_summary_cases_raw <- data.table(cbind(year_list, national_cases_subnat_end_raw))
I_summary_cases_raw$group <- "(raw) case notifications"
colnames(I_summary_cases_raw) <- c("time", "measles", "group")

I_summary <- rbind(I_summary_mean, I_summary_cases, I_summary_cases_raw)
I_summary_all <- rbind(I_summary)


### make plot
p2 <- ggplot(I_summary_all, aes(x=time, y=measles, group=group, color=group)) +  geom_rect(aes(xmin = 2013, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1),fill = "lightgrey") + 
  geom_line(aes(linetype=group))+ theme_bw() + xlab("Year") + ylab("Cases") +
  geom_point(aes(shape=group)) + ggtitle("Cases") 

png(file = paste0("FILEPATH/",run_date,"/predict/",iso3,"p1_case_time_series_with_subnat.png"),
    width = 10,
    height = 7,
    units = "in",
    res = 400,
    type = "cairo")

print(p2)
dev.off()



####### Make alt national level plots with implied scalar from models


eth_ad2 <- sf::st_read(paste0("FILEPATH/","2020_05_21","/lbd_standard_admin_2.shp"))
eth_ad2 <- subset(eth_ad2, NAME_0 == "Ethiopia")
district_labels <- colnames(births_matrix)


############################################################################################################
##### age pattern over time from estimated cases

I_annual_by_district_by_week_0_4 <- matrix(0, 2120, 79)
I_annual_by_district_by_week_5_9 <- matrix(0, 2120, 79)
I_annual_by_district_by_week_10_14 <- matrix(0, 2120, 79)
for (i in 1:2120){
  for(d in 1:79){
    for(a in 1:16){
      I_annual_by_district_by_week_0_4[i,d] = I_annual_by_district_by_week_0_4[i, d] + I_final_counts[i, a, d]
    }
    I_annual_by_district_by_week_5_9[i,d] = I_annual_by_district_by_week_5_9[i, d] + I_final_counts[i, 17, d]
    I_annual_by_district_by_week_10_14[i,d] = I_annual_by_district_by_week_10_14[i, d] + I_final_counts[i, 18, d]
  }
}

I_annual_by_district_0_4 <- matrix(0, 40, 79)
I_annual_by_district_5_9 <- matrix(0, 40, 79)
I_annual_by_district_10_14 <- matrix(0, 40, 79)

for(y in 2:40){
  for (i in ((y-1)*53 + 1):(y*53)){
    for(d in 1:79){
      for(a in 1:24){
        I_annual_by_district_0_4[y,d] = I_annual_by_district_0_4[y, d] + I_annual_by_district_by_week_0_4[i, d]
        I_annual_by_district_5_9[y,d] = I_annual_by_district_5_9[y, d] + I_annual_by_district_by_week_5_9[i, d]
        I_annual_by_district_10_14[y,d] = I_annual_by_district_10_14[y, d] + I_annual_by_district_by_week_10_14[i, d]
      }
    }
  }
}

cases_annual_by_district_0_4 <- matrix(0, 40, 79)
cases_annual_by_district_5_9 <- matrix(0, 40, 79)
cases_annual_by_district_10_14 <- matrix(0, 40, 79)

for(y in 2:40){
  for (i in ((y-1)*53 + 1):(y*53)){
    for(d in 1:79){
      for(a in 1:24){
        cases_annual_by_district_0_4[y,d] = cases_annual_by_district_0_4[y, d] + cases_age_bin0_4_matrix[i, d]
        cases_annual_by_district_5_9[y,d] = cases_annual_by_district_5_9[y, d] + cases_age_bin5_9_matrix[i, d]
        cases_annual_by_district_10_14[y,d] = cases_annual_by_district_10_14[y, d] + cases_age_bin10_14_matrix[i, d]
      }
    }
  }
}


dir.create(paste0("FILEPATH/",run_date,"/predict/age_patterns/reported_cases"), recursive=T)
dir.create(paste0("FILEPATH/",run_date,"/predict/age_patterns/estimated_cases"), recursive=T)

for(d in 1:79){
  message(d)
  district_name <- subset(eth_ad2, ADM2_CODE == district_labels[d])$ADM2_NAME
  region_name <- subset(eth_ad2, ADM2_CODE == district_labels[d])$ADM1_NAME
  
  cases_0_4_df <- data.table(I_annual_by_district_0_4[,d])
  cases_0_4_df$age_bin <- '0_4'
  cases_0_4_df$year <- 1980:2019
  colnames(cases_0_4_df) <- c("cases", "age_bin", "year")
  
  cases_5_9_df <- data.table(I_annual_by_district_5_9[,d])
  cases_5_9_df$age_bin <- '5_9'
  cases_5_9_df$year <- 1980:2019
  colnames(cases_5_9_df) <- c("cases", "age_bin", "year")
  
  cases_10_14_df <- data.table(I_annual_by_district_10_14[,d])
  cases_10_14_df$age_bin <- '10_14'
  cases_10_14_df$year <- 1980:2019
  colnames(cases_10_14_df) <- c("cases", "age_bin", "year")
  
  cases_df <- rbind(cases_0_4_df, cases_5_9_df, cases_10_14_df)
  
  cases_df$age_bin <- factor(cases_df$age_bin, levels = c("0_4", "5_9", "10_14"))
  
  cases_df <- subset(cases_df, year > 1980)
  
  cases_df$cases_adjusted <- cases_df$cases 
  
  gg1 <- ggplot(data = subset(cases_df, year > 1989), aes(x=age_bin, y=cases_adjusted)) + 
    ggtitle(paste0("Estimated cases -- ", district_name, " -- ", region_name)) + 
    geom_bar( stat='identity', aes(fill = age_bin)) + facet_wrap(~year, scales = 'free') +
    theme_bw() + xlab("Age group") + ylab("Estimated cases")
  
  rep_cases_0_4_df <- data.table(cases_annual_by_district_0_4[,d])
  rep_cases_0_4_df$age_bin <- '0_4'
  rep_cases_0_4_df$year <- 1980:2019
  colnames(rep_cases_0_4_df) <- c("cases", "age_bin", "year")
  
  rep_cases_5_9_df <- data.table(cases_annual_by_district_5_9[,d])
  rep_cases_5_9_df$age_bin <- '5_9'
  rep_cases_5_9_df$year <- 1980:2019
  colnames(rep_cases_5_9_df) <- c("cases", "age_bin", "year")
  
  
  rep_cases_10_14_df <- data.table(cases_annual_by_district_10_14[,d])
  rep_cases_10_14_df$age_bin <- '10_14'
  rep_cases_10_14_df$year <- 1980:2019
  colnames(rep_cases_10_14_df) <- c("cases", "age_bin", "year")
  
  rep_cases_df <- rbind(rep_cases_0_4_df, rep_cases_5_9_df, rep_cases_10_14_df)
  
  rep_cases_df$age_bin <- factor(rep_cases_df$age_bin, levels = c("0_4", "5_9", "10_14"))
  
  rep_cases_df <- subset(rep_cases_df, year > 2012)
  
  gg2 <- ggplot(data = rep_cases_df, aes(x=age_bin, y=cases )) + ggtitle(paste0("Reported cases -- ", district_name, " -- ", region_name)) + 
    geom_bar( stat='identity', aes(fill = age_bin)) + facet_wrap(~year) +
    theme_bw() + xlab("Age group") + ylab("Reported Cases")
  
  
  png(file = paste0("FILEPATH/",run_date,"/predict/age_patterns/estimated_cases/",iso3,"_",district_name,".png"),
      width = 10,
      height = 7,
      units = "in",
      res = 400,
      type = "cairo")
  print(gg1)
  dev.off()
  
  png(file = paste0("FILEPATH/",run_date,"/predict/age_patterns/reported_cases/",iso3,"_",district_name,".png"),
      width = 10,
      height = 7,
      units = "in",
      res = 400,
      type = "cairo")
  print(gg2)
  dev.off()
  
}





############################################################################################################
#### let's do by district -- weekly plots AGE

I_annual_by_district_by_week_0_4 <- matrix(0, 2120, 79)
I_annual_by_district_by_week_5_9 <- matrix(0, 2120, 79)
I_annual_by_district_by_week_10_14 <- matrix(0, 2120, 79)
for (i in 1:2120){
  for(d in 1:79){
    for(a in 1:16){
      I_annual_by_district_by_week_0_4[i,d] = I_annual_by_district_by_week_0_4[i, d] + I_final_counts[i, a, d]
    }
    
    I_annual_by_district_by_week_5_9[i,d] = I_annual_by_district_by_week_5_9[i, d] + I_final_counts[i, 17, d]
    
    
    I_annual_by_district_by_week_10_14[i,d] = I_annual_by_district_by_week_10_14[i, d] + I_final_counts[i, 18, d]
    
  }
}






dir.create(paste0("FILEPATH/",run_date,"/predict/district_weekly_cases_subnational_age_0_4_downadjusted/"), recursive=T)
dir.create(paste0("FILEPATH/",run_date,"/predict/district_weekly_cases_subnational_age_5_9_downadjusted/"), recursive=T)
dir.create(paste0("FILEPATH/",run_date,"/predict/district_weekly_cases_subnational_age_10_14_downadjusted/"), recursive=T)



for(d in 1:79){
  message(d)
  
  
  district_name <- subset(eth_ad2, ADM2_CODE == district_labels[d])$ADM2_NAME
  region_name <- subset(eth_ad2, ADM2_CODE == district_labels[d])$ADM1_NAME
  
  ## 0 - 4
  I_annual_for_sim_data <- (I_annual_by_district_by_week_0_4[,d]) 
  results <- data.table(c(rep(NA, 5), I_annual_for_sim_data[6:2120]), 1:2120)
  
  ### make plots
  I_summary_mean <- results
  I_summary_mean$group <- "TSIR"
  colnames(I_summary_mean) <- c("measles", "time", "group")
  if(d %in% c(1)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho1_median )
  }
  if(d %in% c(2, 12, 20, 28, 34)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho2_median )
  }
  if(d %in% c(3, 13, 21, 29, 35, 40, 44, 48, 52, 56, 59, 62)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho3_median )
  }
  if(d %in% c(4, 14, 22)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho4_median )
  }
  if(d %in% c(5)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho5_median )
  }
  if(d %in% c(6, 15, 23)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho6_median )
  }
  if(d %in% c(7)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho7_median )
  }
  if(d %in% c(8, 16, 24, 30, 36, 41, 45, 49, 53, 57, 60, 63, 65, 67, 69, 71, 73)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho8_median )
  }
  if(d %in% c(9, 17, 25, 31, 37, 42, 46, 50, 54)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho9_median )
  }
  if(d %in% c(10, 18, 26, 32, 38, 43, 47, 51, 55, 58, 61, 64, 66, 68, 70, 72, 74, 75, 76, 77, 78, 79)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho10_median )
  }
  if(d %in% c(11, 19, 27, 33, 39)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho11_median )
  }
  
  
  cases_age_bin0_4_matrix2 <- matrix(NA, 2120,79)
  cases_age_bin0_4_matrix2[1:2098,] <- cases_age_bin0_4_matrix[1:2098,] 
  
  cases_age_bin0_4_matrix <- as.matrix(cases_age_bin0_4_matrix)
  I_summary_cases_raw <- data.table(cbind(1:2120, c(rep(NA, 1749), cases_age_bin0_4_matrix[1750:2120,d] )))
  I_summary_cases_raw$group <- "(raw) case notifications"
  colnames(I_summary_cases_raw) <- c("time", "measles", "group")
  
  I_summary <- rbind(I_summary_mean, I_summary_cases_raw)
  I_summary_all <- rbind(I_summary)
  
  ### make plot
  p2 <- ggplot() +  
    #geom_rect(aes(xmin = 0, xmax = 1749, ymin = -Inf, ymax = Inf), alpha = 0.1, color='lightgrey', fill = "lightgrey") + 
    theme_bw() + xlab("Year") + ylab("Cases") + ggtitle(paste0(district_name, " -- ", region_name, " -- 0-4 yos")  ) + 
    geom_line(data = subset(I_summary_all, time > 1749 & group != "TSIR"), mapping = aes(x=(time / 53) + 1980, y=measles, group=group, color=group)) + 
    geom_line(data = subset(I_summary_all, time > 1749 & group == "TSIR"), lwd = 1, mapping = aes(x=(time / 53) + 1980, y=measles, group=group, color=group))
  
  png(file = paste0("FILEPATH/",run_date,"/predict/district_weekly_cases_subnational_age_0_4_downadjusted/",iso3,"p1_weekly_case_time_series_with_subnat_",district_name,".png"),
      width = 10,
      height = 7,
      units = "in",
      res = 400,
      type = "cairo")
  print(p2)
  dev.off()
  
  
  
  ## 5 - 9
  I_annual_for_sim_data <- (I_annual_by_district_by_week_5_9[,d])
  results <- data.table(c(rep(NA, 5), I_annual_for_sim_data[6:2120]), 1:2120)
  
  ### make plots
  I_summary_mean <- results
  I_summary_mean$group <- "TSIR"
  colnames(I_summary_mean) <- c("measles", "time", "group")
  
  if(d %in% c(1)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho1_median)
  }
  if(d %in% c(2, 12, 20, 28, 34)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho2_median )
  }
  if(d %in% c(3, 13, 21, 29, 35, 40, 44, 48, 52, 56, 59, 62)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho3_median )
  }
  if(d %in% c(4, 14, 22)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho4_median )
  }
  if(d %in% c(5)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho5_median )
  }
  if(d %in% c(6, 15, 23)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho6_median )
  }
  if(d %in% c(7)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho7_median )
  }
  if(d %in% c(8, 16, 24, 30, 36, 41, 45, 49, 53, 57, 60, 63, 65, 67, 69, 71, 73)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho8_median )
  }
  if(d %in% c(9, 17, 25, 31, 37, 42, 46, 50, 54)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho9_median )
  }
  if(d %in% c(10, 18, 26, 32, 38, 43, 47, 51, 55, 58, 61, 64, 66, 68, 70, 72, 74, 75, 76, 77, 78, 79)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho10_median )
  }
  if(d %in% c(11, 19, 27, 33, 39)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho11_median )
  }
  
  # I_summary_cases <- data.table(cbind(1:2120, c(rep(NA, 1749), cases_age_bin5_9_matrix[1750:2120,d] / invlogit(logit_rho_median))))
  # I_summary_cases$group <- "(adjusted) case notifications"
  # colnames(I_summary_cases) <- c("time", "measles", "group")
  
  cases_age_bin5_9_matrix2 <- matrix(NA, 2120,79)
  cases_age_bin5_9_matrix2[1:2098,] <- cases_age_bin5_9_matrix[1:2098,]
  
  cases_age_bin5_9_matrix <- as.matrix(cases_age_bin5_9_matrix)
  I_summary_cases_raw <- data.table(cbind(1:2120, c(rep(NA, 1749), cases_age_bin5_9_matrix[1750:2120,d] )))
  I_summary_cases_raw$group <- "(raw) case notifications"
  colnames(I_summary_cases_raw) <- c("time", "measles", "group")
  
  I_summary <- rbind(I_summary_mean, I_summary_cases_raw)
  I_summary_all <- rbind(I_summary)
  
  ### make plot
  p2 <- ggplot() +
    #geom_rect(aes(xmin = 0, xmax = 1749, ymin = -Inf, ymax = Inf), alpha = 0.1, color='lightgrey', fill = "lightgrey") +
    theme_bw() + xlab("Year") + ylab("Cases") + ggtitle(paste0(district_name, " -- ", region_name, " -- 5-9 yos") ) +
    geom_line(data = subset(I_summary_all, time > 1749 & group != "TSIR"), mapping = aes(x=(time / 53) + 1980, y=measles, group=group, color=group)) +
    geom_line(data = subset(I_summary_all, time > 1749 & group == "TSIR"), lwd = 1, mapping = aes(x=(time / 53) + 1980, y=measles, group=group, color=group))
  
  png(file = paste0("FILEPATH/",run_date,"/predict/district_weekly_cases_subnational_age_5_9_downadjusted/",iso3,"p1_weekly_case_time_series_with_subnat_",district_name,".png"),
      width = 10,
      height = 7,
      units = "in",
      res = 400,
      type = "cairo")
  print(p2)
  dev.off()
  
  
  ## 10 - 14
  I_annual_for_sim_data <- (I_annual_by_district_by_week_10_14[,d])
  results <- data.table(c(rep(NA, 5), I_annual_for_sim_data[6:2120]), 1:2120)
  
  ### make plots
  I_summary_mean <- results
  I_summary_mean$group <- "TSIR"
  colnames(I_summary_mean) <- c("measles", "time", "group")
  
  if(d %in% c(1)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho1_median )
  }
  if(d %in% c(2, 12, 20, 28, 34)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho2_median )
  }
  if(d %in% c(3, 13, 21, 29, 35, 40, 44, 48, 52, 56, 59, 62)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho3_median )
  }
  if(d %in% c(4, 14, 22)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho4_median )
  }
  if(d %in% c(5)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho5_median )
  }
  if(d %in% c(6, 15, 23)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho6_median )
  }
  if(d %in% c(7)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho7_median )
  }
  if(d %in% c(8, 16, 24, 30, 36, 41, 45, 49, 53, 57, 60, 63, 65, 67, 69, 71, 73)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho8_median )
  }
  if(d %in% c(9, 17, 25, 31, 37, 42, 46, 50, 54)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho9_median )
  }
  if(d %in% c(10, 18, 26, 32, 38, 43, 47, 51, 55, 58, 61, 64, 66, 68, 70, 72, 74, 75, 76, 77, 78, 79)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho10_median )
  }
  if(d %in% c(11, 19, 27, 33, 39)){
    I_summary_mean$measles <- I_summary_mean$measles * invlogit(logit_rho11_median )
  }
  
  # I_summary_cases <- data.table(cbind(1:2120, c(rep(NA, 1749), cases_age_bin10_14_matrix[1750:2120,d] / invlogit(logit_rho_median))))
  # I_summary_cases$group <- "(adjusted) case notifications"
  # colnames(I_summary_cases) <- c("time", "measles", "group")
  
  
  cases_age_bin10_14_matrix2 <- matrix(NA, 2120,79)
  cases_age_bin10_14_matrix2[1:2098,] <- cases_age_bin10_14_matrix[1:2098,]
  #
  cases_age_bin10_14_matrix <- as.matrix(cases_age_bin10_14_matrix)
  I_summary_cases_raw <- data.table(cbind(1:2120, c(rep(NA, 1749), cases_age_bin10_14_matrix[1750:2120,d] )))
  I_summary_cases_raw$group <- "(raw) case notifications"
  colnames(I_summary_cases_raw) <- c("time", "measles", "group")
  
  I_summary <- rbind(I_summary_mean, I_summary_cases_raw)
  I_summary_all <- rbind(I_summary)
  
  ### make plot
  p2 <- ggplot() +
    #geom_rect(aes(xmin = 0, xmax = 1749, ymin = -Inf, ymax = Inf), alpha = 0.1, color='lightgrey', fill = "lightgrey") +
    theme_bw() + xlab("Year") + ylab("Cases") + ggtitle(paste0(district_name, " -- ", region_name, " -- 10-14 yos")  ) +
    geom_line(data = subset(I_summary_all, time > 1749 & group != "TSIR"), mapping = aes(x=(time / 53) + 1980, y=measles, group=group, color=group)) +
    geom_line(data = subset(I_summary_all, time > 1749 & group == "TSIR"), lwd = 1, mapping = aes(x=(time / 53) + 1980, y=measles, group=group, color=group))
  
  png(file = paste0("FILEPATH/",run_date,"/predict/district_weekly_cases_subnational_age_10_14_downadjusted/",iso3,"p1_weekly_case_time_series_with_subnat_",district_name,".png"),
      width = 10,
      height = 7,
      units = "in",
      res = 400,
      type = "cairo")
  print(p2)
  dev.off()
  
  
  
  
}



############################################################################################################    
### make MSIR plots 

make_MSIR_plots <- T

if(make_MSIR_plots){
  
  cols <- c("0_5_mo" = "#9e0142", "6_8_mo" = "#d53e4f", "9_11_mo" = "#f46d43", "1_yo" = "#fdae61",
            "2_yo" = "#fee08b", "3_yo" = "#ffffbf", "4_yo" = "#e6f598", "5_14_yo" = "#abdda4",
            "15_44_yo" = "#66c2a5", "45_64_yo" = "#3288bd", "65+_yo" = "#5e4fa2")
  #########################################################################
  pop_by_age_week <- matrix(0,2120,24)
  M_by_age_week <- matrix(0,2120,24)
  S_by_age_week <- matrix(0,2120,24)
  I_by_age_week <- matrix(0,2120,24)
  R_by_age_week <- matrix(0,2120,24)
  R_vax_by_age_week <- matrix(0,2120,24)
  R_vax2_by_age_week <- matrix(0,2120,24)
  for(i in 1:2120){
    for(a in 1:24){
      for(d in 1:79){
        pop_by_age_week[i,a] = pop_by_age_week[i,a] + pop_array_red[i,a,d]
        M_by_age_week[i,a]= M_by_age_week[i,a] + M_final_counts[i,a,d]
        S_by_age_week[i,a]= S_by_age_week[i,a] + S_final_counts[i,a,d]
        I_by_age_week[i,a]= I_by_age_week[i,a] + I_final_counts[i,a,d]
        R_by_age_week[i,a]= R_by_age_week[i,a] + R_final_counts[i,a,d]
        # R_vax_by_age_week[i,a]= R_vax_by_age_week[i,a] + R_vax_final_counts[i,a,d]
        # R_vax2_by_age_week[i,a]= R_vax2_by_age_week[i,a] + R_vax2_final_counts[i,a,d]
        
      }  
    }
  }
  
  data <- as.data.frame(pop_by_age_week)
  
  ##### format populations
  data <- as.data.frame(pop_by_age_week)
  data$id <- 1:nrow(data) 
  plot_data <- melt(data,id.var=c("id"))
  
  plot_data$age_bin <- c(rep("0_5_mo", 6*2120), rep("6_8_mo", 3*2120), rep("9_11_mo", 3*2120), rep("1_yo", 1*2120),
                         rep("2_yo", 1*2120), rep("3_yo", 1*2120), rep("4_yo", 1*2120), rep("5_14_yo", 2*2120),
                         rep("15_44_yo", 3*2120),rep("45_64_yo", 2*2120), rep("65+_yo", 1*2120))
  
  population_data_v2 <- aggregate(plot_data$value, by=list(plot_data$id, plot_data$age_bin), FUN= "sum")
  colnames(population_data_v2) <- c("week", "age_bin", "population")
  
  ##### susceptibles
  data <- as.data.frame(S_by_age_week)
  data$id <- 1:nrow(data)
  plot_data <- melt(data,id.var=c("id"))
  
  plot_data$age_bin <- c(rep("0_5_mo", 6*2120), rep("6_8_mo", 3*2120), rep("9_11_mo", 3*2120), rep("1_yo", 1*2120),
                         rep("2_yo", 1*2120), rep("3_yo", 1*2120), rep("4_yo", 1*2120), rep("5_14_yo", 2*2120),
                         rep("15_44_yo", 3*2120),rep("45_64_yo", 2*2120), rep("65+_yo", 1*2120))
  
  plot_data_v2 <- aggregate(plot_data$value, by=list(plot_data$id, plot_data$age_bin), FUN= "sum")
  colnames(plot_data_v2) <- c("week", "age_bin", "value")
  
  plot_data_v3 <- merge(plot_data_v2, population_data_v2, by=c("week", "age_bin") )
  
  plot_data_v3$proportion <- plot_data_v3$value / plot_data_v3$population
  
  plot_data_v3$age_bin <- factor(plot_data_v3$age_bin, levels=c("0_5_mo", "6_8_mo", "9_11_mo", "1_yo", "2_yo", "3_yo", "4_yo",
                                                                "5_14_yo", "15_44_yo", "45_64_yo", "65+_yo"))
  
  p2_S_by_age <- ggplot(plot_data_v3, aes(x=week,y=proportion,group=age_bin,colour=age_bin)) +
    geom_line() + ylab("Proportion") + xlab("Week") + scale_colour_manual(values = cols) +
    ggtitle("Susceptible")
  
  ##### maternally immune
  data <- as.data.frame(M_by_age_week)
  data$id <- 1:nrow(data)
  plot_data <- melt(data,id.var=c("id"))
  
  plot_data$age_bin <- c(rep("0_5_mo", 6*2120), rep("6_8_mo", 3*2120), rep("9_11_mo", 3*2120), rep("1_yo", 1*2120),
                         rep("2_yo", 1*2120), rep("3_yo", 1*2120), rep("4_yo", 1*2120), rep("5_14_yo", 2*2120),
                         rep("15_44_yo", 3*2120),rep("45_64_yo", 2*2120), rep("65+_yo", 1*2120))
  
  plot_data_v2 <- aggregate(plot_data$value, by=list(plot_data$id, plot_data$age_bin), FUN= "sum")
  colnames(plot_data_v2) <- c("week", "age_bin", "value")
  
  plot_data_v3 <- merge(plot_data_v2, population_data_v2, by=c("week", "age_bin") )
  
  plot_data_v3$proportion <- plot_data_v3$value / plot_data_v3$population
  
  plot_data_v3$age_bin <- factor(plot_data_v3$age_bin, levels=c("0_5_mo", "6_8_mo", "9_11_mo", "1_yo", "2_yo", "3_yo", "4_yo",
                                                                "5_14_yo", "15_44_yo", "45_64_yo", "65+_yo"))
  
  p2_M_by_age <- ggplot(plot_data_v3, aes(x=week,y=proportion,group=age_bin,colour=age_bin)) +
    geom_line() + ylab("Proportion") + xlab("Week") + scale_colour_manual(values = cols) +
    ggtitle("Maternally immune")
  
  
  ##### infected
  data <- as.data.frame(I_by_age_week)
  data$id <- 1:nrow(data) 
  plot_data <- melt(data,id.var=c("id"))
  
  plot_data$age_bin <- c(rep("0_5_mo", 6*2120), rep("6_8_mo", 3*2120), rep("9_11_mo", 3*2120), rep("1_yo", 1*2120),
                         rep("2_yo", 1*2120), rep("3_yo", 1*2120), rep("4_yo", 1*2120), rep("5_14_yo", 2*2120),
                         rep("15_44_yo", 3*2120),rep("45_64_yo", 2*2120), rep("65+_yo", 1*2120))
  
  plot_data_v2 <- aggregate(plot_data$value, by=list(plot_data$id, plot_data$age_bin), FUN= "sum")
  colnames(plot_data_v2) <- c("week", "age_bin", "value")
  
  plot_data_v3 <- merge(plot_data_v2, population_data_v2, by=c("week", "age_bin") )
  
  plot_data_v3$proportion <- plot_data_v3$value / plot_data_v3$population
  
  plot_data_v3$age_bin <- factor(plot_data_v3$age_bin, levels=c("0_5_mo", "6_8_mo", "9_11_mo", "1_yo", "2_yo", "3_yo", "4_yo",
                                                                "5_14_yo", "15_44_yo", "45_64_yo", "65+_yo"))
  
  p2_I_by_age <- ggplot(subset(plot_data_v3, week > 0), aes(x=week,y=proportion,group=age_bin,colour=age_bin)) +
    geom_line() + ylab("Proportion") + xlab("Week") + scale_colour_manual(values = cols) + 
    ggtitle("Infected")
  
  ##### recovered -- all
  data <- as.data.frame(R_by_age_week)
  data$id <- 1:nrow(data)
  plot_data <- melt(data,id.var=c("id"))
  
  plot_data$age_bin <- c(rep("0_5_mo", 6*2120), rep("6_8_mo", 3*2120), rep("9_11_mo", 3*2120), rep("1_yo", 1*2120),
                         rep("2_yo", 1*2120), rep("3_yo", 1*2120), rep("4_yo", 1*2120), rep("5_14_yo", 2*2120),
                         rep("15_44_yo", 3*2120),rep("45_64_yo", 2*2120), rep("65+_yo", 1*2120))
  
  plot_data_v2 <- aggregate(plot_data$value, by=list(plot_data$id, plot_data$age_bin), FUN= "sum")
  colnames(plot_data_v2) <- c("week", "age_bin", "value")
  
  plot_data_v3 <- merge(plot_data_v2, population_data_v2, by=c("week", "age_bin") )
  
  plot_data_v3$proportion <- plot_data_v3$value / plot_data_v3$population
  
  plot_data_v3$age_bin <- factor(plot_data_v3$age_bin, levels=c("0_5_mo", "6_8_mo", "9_11_mo", "1_yo", "2_yo", "3_yo", "4_yo",
                                                                "5_14_yo", "15_44_yo", "45_64_yo", "65+_yo"))
  
  p2_R_by_age <- ggplot(plot_data_v3, aes(x=week,y=proportion,group=age_bin,colour=age_bin)) +
    geom_line() + ylab("Proportion") + xlab("Week") + scale_colour_manual(values = cols) +
    ggtitle("Recovered")
  # 
  # 
  ##### make plots
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  the_legend <- g_legend(p2_R_by_age)
  
  p2_M_by_age <- p2_M_by_age + theme(legend.position = "none")
  p2_S_by_age <- p2_S_by_age + theme(legend.position = "none")
  p2_I_by_age <- p2_I_by_age + theme(legend.position = "none")
  p2_R_by_age <- p2_R_by_age + theme(legend.position = "none")
  # p2_R_vax_by_age <- p2_R_vax_by_age + theme(legend.position = "none")
  # p2_R_vax2_by_age <- p2_R_vax2_by_age + theme(legend.position = "none")
  
  lay <- rbind(c(1,1,1,1,2,2,2,2,5),
               c(1,1,1,1,2,2,2,2,5),
               c(1,1,1,1,2,2,2,2,5),
               c(3,3,3,3,4,4,4,4,5),
               c(3,3,3,3,4,4,4,4,5),
               c(3,3,3,3,4,4,4,4,5))
  
  plot_all <- arrangeGrob(p2_M_by_age, p2_S_by_age, p2_I_by_age, p2_R_by_age, the_legend,
                          layout_matrix = lay,
                          heights = c(1,1,1,1,1,1))
  
  png(file = paste0("FILEPATH/",run_date,"/predict/",iso3,"p2_MSIR_by_age_time.png"),
      width = 14,
      height = 9,
      units = "in",
      res = 400,
      type = "cairo")
  grid.draw(plot_all)
  dev.off()
  
} # end if(make_MSIR_plots)

############################################################################################################    
### make plots -- compartments by district 

make_MSIR_plots_by_district <- T


eth_ad2 <- sf::st_read(paste0("FILEPATH/","2020_05_21","/lbd_standard_admin_2.shp"))
eth_ad2 <- subset(eth_ad2, NAME_0 == "Ethiopia")
district_labels <- colnames(births_matrix)

if(make_MSIR_plots_by_district){
  dir.create(paste0("FILEPATH/",run_date,"/predict/district_compartments/"))
  for(d in 1:79){
    
    message(d)
    ############################################################################################################    
    ### make plots -- SUSCEPTIBLE
    cols <- c("0_5_mo" = "#9e0142", "6_8_mo" = "#d53e4f", "9_11_mo" = "#f46d43", "1_yo" = "#fdae61",
              "2_yo" = "#fee08b", "3_yo" = "#ffffbf", "4_yo" = "#e6f598", "5_14_yo" = "#abdda4",
              "15_44_yo" = "#66c2a5", "45_64_yo" = "#3288bd", "65+_yo" = "#5e4fa2")
    
    pop_by_age_week <- matrix(0,2120,24)
    M_by_age_week <- matrix(0,2120,24)
    S_by_age_week <- matrix(0,2120,24)
    I_by_age_week <- matrix(0,2120,24)
    R_by_age_week <- matrix(0,2120,24)
    for(i in 1:2120){
      for(a in 1:24){
        
        pop_by_age_week[i,a] =  pop_array_red[i,a,d]
        M_by_age_week[i,a]=  M_final_counts[i,a,d]
        S_by_age_week[i,a]= S_final_counts[i,a,d]
        I_by_age_week[i,a]= I_final_counts[i,a,d]
        R_by_age_week[i,a]= R_final_counts[i,a,d]
        
      }  
    }
    
    data <- as.data.frame(pop_by_age_week)
    data$id <- 1:nrow(data) 
    plot_data <- melt(data,id.var=c("id"))
    
    plot_data$age_bin <- c(rep("0_5_mo", 6*2120), rep("6_8_mo", 3*2120), rep("9_11_mo", 3*2120), rep("1_yo", 1*2120),
                           rep("2_yo", 1*2120), rep("3_yo", 1*2120), rep("4_yo", 1*2120), rep("5_14_yo", 2*2120),
                           rep("15_44_yo", 3*2120),rep("45_64_yo", 2*2120), rep("65+_yo", 1*2120))
    
    
    population_data_v2 <- aggregate(plot_data$value, by=list(plot_data$id, plot_data$age_bin), FUN= "sum")
    colnames(population_data_v2) <- c("week", "age_bin", "population")
    
    
    ### susceptibles
    data <- as.data.frame(S_by_age_week)
    data$id <- 1:nrow(data)
    plot_data <- melt(data,id.var=c("id"))
    
    plot_data$age_bin <- c(rep("0_5_mo", 6*2120), rep("6_8_mo", 3*2120), rep("9_11_mo", 3*2120), rep("1_yo", 1*2120),
                           rep("2_yo", 1*2120), rep("3_yo", 1*2120), rep("4_yo", 1*2120), rep("5_14_yo", 2*2120),
                           rep("15_44_yo", 3*2120),rep("45_64_yo", 2*2120), rep("65+_yo", 1*2120))
    
    plot_data_v2 <- aggregate(plot_data$value, by=list(plot_data$id, plot_data$age_bin), FUN= "sum")
    colnames(plot_data_v2) <- c("week", "age_bin", "value")
    
    plot_data_v3 <- merge(plot_data_v2, population_data_v2, by=c("week", "age_bin") )
    
    plot_data_v3$proportion <- plot_data_v3$value / plot_data_v3$population
    
    plot_data_v3$age_bin <- factor(plot_data_v3$age_bin, levels=c("0_5_mo", "6_8_mo", "9_11_mo", "1_yo", "2_yo", "3_yo", "4_yo",
                                                                  "5_14_yo", "15_44_yo", "45_64_yo", "65+_yo"))
    
    p2_S_by_age <- ggplot(plot_data_v3, aes(x=week,y=proportion,group=age_bin,colour=age_bin)) +
      geom_line() + ylab("Proportion") + xlab("Week") + scale_colour_manual(values = cols) +
      ggtitle("Susceptible")
    
    ### maternally immune
    data <- as.data.frame(M_by_age_week)
    data$id <- 1:nrow(data)
    plot_data <- melt(data,id.var=c("id"))
    
    plot_data$age_bin <- c(rep("0_5_mo", 6*2120), rep("6_8_mo", 3*2120), rep("9_11_mo", 3*2120), rep("1_yo", 1*2120),
                           rep("2_yo", 1*2120), rep("3_yo", 1*2120), rep("4_yo", 1*2120), rep("5_14_yo", 2*2120),
                           rep("15_44_yo", 3*2120),rep("45_64_yo", 2*2120), rep("65+_yo", 1*2120))
    
    plot_data_v2 <- aggregate(plot_data$value, by=list(plot_data$id, plot_data$age_bin), FUN= "sum")
    colnames(plot_data_v2) <- c("week", "age_bin", "value")
    
    plot_data_v3 <- merge(plot_data_v2, population_data_v2, by=c("week", "age_bin") )
    
    plot_data_v3$proportion <- plot_data_v3$value / plot_data_v3$population
    
    plot_data_v3$age_bin <- factor(plot_data_v3$age_bin, levels=c("0_5_mo", "6_8_mo", "9_11_mo", "1_yo", "2_yo", "3_yo", "4_yo",
                                                                  "5_14_yo", "15_44_yo", "45_64_yo", "65+_yo"))
    
    p2_M_by_age <- ggplot(plot_data_v3, aes(x=week,y=proportion,group=age_bin,colour=age_bin)) +
      geom_line() + ylab("Proportion") + xlab("Week") + scale_colour_manual(values = cols) +
      ggtitle("Maternally immune")
    
    
    ### infected
    data <- as.data.frame(I_by_age_week)
    data$id <- 1:nrow(data) 
    plot_data <- melt(data,id.var=c("id"))
    
    
    plot_data$age_bin <- c(rep("0_5_mo", 6*2120), rep("6_8_mo", 3*2120), rep("9_11_mo", 3*2120), rep("1_yo", 1*2120),
                           rep("2_yo", 1*2120), rep("3_yo", 1*2120), rep("4_yo", 1*2120), rep("5_14_yo", 2*2120),
                           rep("15_44_yo", 3*2120),rep("45_64_yo", 2*2120), rep("65+_yo", 1*2120))
    
    plot_data_v2 <- aggregate(plot_data$value, by=list(plot_data$id, plot_data$age_bin), FUN= "sum")
    colnames(plot_data_v2) <- c("week", "age_bin", "value")
    
    plot_data_v3 <- merge(plot_data_v2, population_data_v2, by=c("week", "age_bin") )
    
    plot_data_v3$proportion <- plot_data_v3$value / plot_data_v3$population
    
    plot_data_v3$age_bin <- factor(plot_data_v3$age_bin, levels=c("0_5_mo", "6_8_mo", "9_11_mo", "1_yo", "2_yo", "3_yo", "4_yo",
                                                                  "5_14_yo", "15_44_yo", "45_64_yo", "65+_yo"))
    
    p2_I_by_age <- ggplot(subset(plot_data_v3, week > 2), aes(x=week,y=proportion,group=age_bin,colour=age_bin)) +
      geom_line() + ylab("Proportion") + xlab("Week") + scale_colour_manual(values = cols) + 
      ggtitle("Infected")
    
    
    ### recovered -- all
    data <- as.data.frame(R_by_age_week)
    data$id <- 1:nrow(data)
    plot_data <- melt(data,id.var=c("id"))
    
    plot_data$age_bin <- c(rep("0_5_mo", 6*2120), rep("6_8_mo", 3*2120), rep("9_11_mo", 3*2120), rep("1_yo", 1*2120),
                           rep("2_yo", 1*2120), rep("3_yo", 1*2120), rep("4_yo", 1*2120), rep("5_14_yo", 2*2120),
                           rep("15_44_yo", 3*2120),rep("45_64_yo", 2*2120), rep("65+_yo", 1*2120))
    
    plot_data_v2 <- aggregate(plot_data$value, by=list(plot_data$id, plot_data$age_bin), FUN= "sum")
    colnames(plot_data_v2) <- c("week", "age_bin", "value")
    
    plot_data_v3 <- merge(plot_data_v2, population_data_v2, by=c("week", "age_bin") )
    
    plot_data_v3$proportion <- plot_data_v3$value / plot_data_v3$population
    
    plot_data_v3$age_bin <- factor(plot_data_v3$age_bin, levels=c("0_5_mo", "6_8_mo", "9_11_mo", "1_yo", "2_yo", "3_yo", "4_yo",
                                                                  "5_14_yo", "15_44_yo", "45_64_yo", "65+_yo"))
    
    p2_R_by_age <- ggplot(plot_data_v3, aes(x=week,y=proportion,group=age_bin,colour=age_bin)) +
      geom_line() + ylab("Proportion") + xlab("Week") + scale_colour_manual(values = cols) +
      ggtitle("Recovered")
    
    district_name <- subset(eth_ad2, ADM2_CODE == district_labels[d])$ADM2_NAME
    
    
    g_legend<-function(a.gplot){
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)}
    
    the_legend <- g_legend(p2_R_by_age)
    
    titleobject <- textGrob(paste0(district_name),gp=gpar(fontsize=20))
    
    p2_M_by_age <- p2_M_by_age + theme(legend.position = "none")
    p2_S_by_age <- p2_S_by_age + theme(legend.position = "none")
    p2_I_by_age <- p2_I_by_age + theme(legend.position = "none")
    p2_R_by_age <- p2_R_by_age + theme(legend.position = "none")
    
    lay <- rbind(c(1,1,1,1,1,1,1,1,1),
                 c(2,2,2,2,3,3,3,3,6),
                 c(2,2,2,2,3,3,3,3,6),
                 c(2,2,2,2,3,3,3,3,6),
                 c(4,4,4,4,5,5,5,5,6),
                 c(4,4,4,4,5,5,5,5,6),
                 c(4,4,4,4,5,5,5,5,6))
    
    plot_all <- arrangeGrob(titleobject, p2_M_by_age, p2_S_by_age, p2_I_by_age, p2_R_by_age, the_legend,
                            layout_matrix = lay,
                            heights = c(1,1,1,1,1,1, 1))
    
    png(file = paste0("FILEPATH/",run_date,"/predict/district_compartments/",iso3,"p2_MSIR_by_age_time_",district_name,".png"),
        width = 14,
        height = 9,
        units = "in",
        res = 400,
        type = "cairo")
    
    grid.draw(plot_all)
    
    dev.off()
  }
  
  
}





