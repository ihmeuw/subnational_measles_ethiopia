
## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- 'USERNAME'
core_repo          <- sprintf('FILEPATH',user)
vax_repo           <- sprintf('FILEPATH',user)
measles_repo       <- sprintf('FILEPATH',user)
remote             <- 'origin'
branch             <- 'develop'
pullgit            <- FALSE

## sort some directory stuff
commondir      <- sprintf('FILEPATH')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

# Load MBG packages and functions
# message('Loading in required R packages and MBG functions')
# source(paste0(core_repo, '/mbg_central/setup.R'))
# mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0('FILEPATH/misc_vaccine_functions.R'))
library(kernlab)
library(data.table)
library(ggplot2)

library(TMB)
# library(tmbstan)

library(numDeriv)

# Custom load indicator-specific functions
source(paste0(vax_repo,'FILEPATH/misc_vaccine_functions.R'))
library(readxl)


###############################################################################################################   
### set up run date, logdir, and paramters to launch models / predictions
iso3 <- 'ETH'
year_list <- 1980:2021
processing_date <- '2023_01_03'

############################################################################################################
## load processed data
############################################################################################################

# ######################################################

contact <- fread(paste0('FILEPATH/contact_matrix_',iso3,'.csv'))

contact_ori <- as.matrix(contact)#contact_all[[icty]]


###################################################################################################### 
################# lets pull / compute population by epiweek
###################################################################################################### 
source("FILEPATH/get_population.R")

if(iso3 == "ETH") iso3_loc_id <- 179

year_list <- '2019'

population_age_bin0 <- get_population(age_group_id = 28, single_year=T, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin0 <- population_age_bin0$population

population_age_bin1 <- get_population(age_group_id = 238, single_year=T, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin1 <- population_age_bin1$population

population_age_bin2 <- get_population(age_group_id = 50, single_year=T, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin2 <- population_age_bin2$population

population_age_bin3 <- get_population(age_group_id = 51, single_year=T, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin3 <- population_age_bin3$population

population_age_bin4 <- get_population(age_group_id = 52, single_year=T, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin4 <- population_age_bin4$population

population_age_bin5 <- get_population(age_group_id = 53, single_year=T, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin5 <- population_age_bin5$population

population_age_bin6 <- get_population(age_group_id = 54, single_year=T, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin6 <- population_age_bin6$population

population_age_bin7 <- get_population(age_group_id = 55, single_year=T, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin7 <- population_age_bin7$population

population_age_bin8 <- get_population(age_group_id = 56, single_year=T, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin8 <- population_age_bin8$population

population_age_bin9 <- get_population(age_group_id = 57, single_year=T, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin9 <- population_age_bin9$population

population_age_bin10 <- get_population(age_group_id = 58, single_year=T, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin10 <- population_age_bin10$population

population_age_bin11 <- get_population(age_group_id = 59, single_year=T, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin11 <- population_age_bin11$population

population_age_bin12 <- get_population(age_group_id = 60, single_year=T, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin12 <- population_age_bin12$population

population_age_bin13 <- get_population(age_group_id = 61, single_year=T, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin13 <- population_age_bin13$population

population_age_bin14 <- get_population(age_group_id = 62, single_year=T, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin14 <- population_age_bin14$population
# 

#### 15-19 age bin 8
population_age_bin15_19 <- get_population(age_group_id = 8, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin15_19 <- population_age_bin15_19$population

#### 20-24 age bin 9
population_age_bin20_24 <- get_population(age_group_id = 9, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin20_24 <- population_age_bin20_24$population

#### 25-29 age bin 10
population_age_bin25_29 <- get_population(age_group_id = 10, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin25_29 <- population_age_bin25_29$population

#### 30-34 age bin 11
population_age_bin30_34 <- get_population(age_group_id = 11, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin30_34 <- population_age_bin30_34$population

#### 35-39 age bin 12
population_age_bin35_39 <- get_population(age_group_id = 12, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin35_39 <- population_age_bin35_39$population

#### 40-44 age bin 13
population_age_bin40_44 <- get_population(age_group_id = 13, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin40_44 <- population_age_bin40_44$population

#### 45-49 age bin 14
population_age_bin45_49 <- get_population(age_group_id = 14, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin45_49 <- population_age_bin45_49$population

#### 50-54 age bin 15
population_age_bin50_54 <- get_population(age_group_id = 15, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin50_54 <- population_age_bin50_54$population

#### 55-59 age bin 16
population_age_bin55_59 <- get_population(age_group_id = 16, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin55_59 <- population_age_bin55_59$population

#### 60-64 age bin 17
population_age_bin60_64 <- get_population(age_group_id = 17, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin60_64 <- population_age_bin60_64$population

#### 65-69 age bin 18
population_age_bin65_69 <- get_population(age_group_id = 18, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin65_69 <- population_age_bin65_69$population

#### 70-74 age bin 19
population_age_bin70_74 <- get_population(age_group_id = 19, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin70_74 <- population_age_bin70_74$population

#### 75-79 age bin 20
population_age_bin75_79 <- get_population(age_group_id = 20, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin75_79 <- population_age_bin75_79$population

#### 80-84 age bin ---- 30
population_age_bin80_84 <- get_population(age_group_id = 30, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin80_84 <- population_age_bin80_84$population

#### 85-89 age bin ---- 31
population_age_bin85_89 <- get_population(age_group_id = 31, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin85_89 <- population_age_bin85_89$population

#### 90-94 age bin ---- 32
population_age_bin90_94 <- get_population(age_group_id = 32, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin90_94 <- population_age_bin90_94$population

#### 95+ age bin 24
population_age_bin95_plus <- get_population(age_group_id = 235, location_id = iso3_loc_id, sex_id = 3, gbd_round_id = 7, decomp_step="step3",  year_id = year_list)
population_age_bin95_plus <- population_age_bin95_plus$population



population_age_bin5_9 <- population_age_bin5 + population_age_bin6 + population_age_bin7 + population_age_bin8 + 
  population_age_bin9 
population_age_bin10_14 <- population_age_bin10 + population_age_bin11 + population_age_bin12 + population_age_bin13 + 
  population_age_bin14 

population_age_bin15_24 <- population_age_bin15_19 + population_age_bin20_24 
population_age_bin25_34 <- population_age_bin25_29 + population_age_bin30_34 
population_age_bin35_44 <- population_age_bin35_39 + population_age_bin40_44 
population_age_bin45_54 <- population_age_bin45_49 + population_age_bin50_54 
population_age_bin55_64 <- population_age_bin55_59 + population_age_bin60_64 
population_age_bin65_plus <- population_age_bin65_69 + population_age_bin70_74 + population_age_bin75_79 + population_age_bin80_84 +
  population_age_bin85_89 + population_age_bin90_94 + population_age_bin95_plus 

pop_data <- c(population_age_bin0, population_age_bin1, population_age_bin2, population_age_bin3, population_age_bin4,
              population_age_bin5_9, population_age_bin10_14, 
              population_age_bin15_24, population_age_bin25_34, population_age_bin35_44, 
              population_age_bin45_54, population_age_bin55_64, population_age_bin65_plus)


pop_cty  <- pop_data

# expand contactees (break down in to single-year age groups by proportion)
contact_col <- matrix (0, ncol = 13, nrow = dim(contact_ori)[1])

# 0-4 years olds
for (icol in 1){
  pop_prp <- pop_cty [5*(icol-1)+(1:5)] / sum(pop_cty[5*(icol-1)+(1:5)])
  contactees <- rep(contact_ori[, icol], 5)
  contact_col [, 5*(icol-1)+(1:5)] <- sweep (matrix (contactees, ncol = 5), 2, pop_prp, "*")
}

# 5-9 and 10-14 year olds
for (icol in 2:3){
  print(icol)
  contact_col[, icol + 4] <- contact_ori[, icol]
}

# 15-24 year olds
contact_col[, 8] <- contact_ori[, 4] + contact_ori[, 5]

# 25-34 year olds
contact_col[, 9] <- contact_ori[, 6] + contact_ori[, 7]

# 35-44 year olds
contact_col[, 10] <- contact_ori[, 8] + contact_ori[, 9]

# 45-54 year olds
contact_col[, 11] <- contact_ori[, 10] + contact_ori[, 11]

# 55-64 year olds
contact_col[, 12] <- contact_ori[, 12] + contact_ori[, 13]

# 65+ year olds
contact_col[, 13] <- contact_ori[, 14] + contact_ori[, 15] + contact_ori[, 16]

# expand contactors (apply the average to different groups)
contactor_0to4   <- matrix (rep(contact_col[1,], each = 5), nrow = 5)
contactor_5to14   <- matrix (rep(contact_col[2:3,], each = 1), nrow = 2)
contactor_15to24   <- matrix (rep((contact_col[4,] + contact_col[5,]), each = 1), nrow = 1)
contactor_25to34   <- matrix (rep((contact_col[6,] + contact_col[7,]), each = 1), nrow = 1)
contactor_35to44   <- matrix (rep((contact_col[8,] + contact_col[9,]), each = 1), nrow = 1)
contactor_45to54   <- matrix (rep((contact_col[10,] + contact_col[11,]), each = 1), nrow = 1)
contactor_55to64   <- matrix (rep((contact_col[12,] + contact_col[13,]), each = 1), nrow = 1)
contactor_65plus   <- matrix (rep((contact_col[14,] + contact_col[15,] + contact_col[16,]), each = 1), nrow = 1)

contact_13 <- rbind (contactor_0to4, contactor_5to14, contactor_15to24, contactor_25to34, contactor_35to44,
                     contactor_45to54, contactor_55to64, contactor_65plus)

c_contact <- as.matrix(contact_13)

# expand 0-2 years old to weekly age strata
s         <- 12 # number of finer stages within an age band (weekly ages, so 52)
jt        <- 1  # how many ages to expand to s (or to weeks)
beta_full <- matrix (0,
                     ncol = 24,
                     nrow = 24)

# ------------------------------------------------------------------------------
expandMatrix <- function (A,
                          expand_rows  = 1,
                          expand_cols  = 1,
                          rescale_rows = F,
                          rescale_cols = F) {
  
  if(!is.matrix(A)){
    stop("A is not a matrix")
  }
  
  matvals <- numeric(0)
  rows <- nrow(A)
  cols <- ncol(A)
  
  for(c in 1:cols) {
    matvals <- c(
      matvals,
      rep(
        A[,c],
        expand_cols,
        each = expand_rows
      )
    )
  }
  
  B <- matrix (matvals,
               nrow = rows * expand_rows,
               ncol = cols * expand_cols)
  
  if(rescale_rows & rescale_cols){
    B <- B/(expand_rows*expand_cols)
  } else if(rescale_rows){
    B <- B/expand_rows
  } else if(rescale_cols){
    B <- B/expand_cols
  }
  
  return (B)
  
} # end of function -- expandMatrix
# ------------------------------------------------------------------------------

beta_full[(1:(s*jt)), (1:(s*jt))] <- expandMatrix (
  A = as.matrix(c_contact [1:jt, 1:jt]/s),  # needs to be divided by 52 so that the mean total number of contacts stays the same
  expand_rows =  s, expand_cols =  s,
  rescale_rows = FALSE, rescale_cols = FALSE)

beta_full[1:(s*jt),((s*jt)+1):(ncol(beta_full))] <- expandMatrix(
  A = as.matrix(c_contact [1:jt,(jt+1):ncol(c_contact)]),
  expand_rows = s, expand_cols = 1,
  rescale_rows = F, rescale_cols = F)

beta_full[((s*jt)+1):(nrow(beta_full)), 1:(s*jt)] <- expandMatrix(
  A = as.matrix(c_contact [(jt+1):nrow(c_contact),1:jt]/s),  # adjust to ensure the mean total number of contacts stays the same
  expand_rows = 1, expand_cols = s,
  rescale_rows = F, rescale_cols = F)

beta_full[((s*jt)+1):(nrow(beta_full)), ((s*jt)+1):(ncol(beta_full))] <-
  c_contact [(jt+1):nrow(c_contact),(jt+1):ncol(c_contact)]

contact_matrix <- beta_full
save(contact_matrix, file=paste0('FILEPATH', iso3, '/',processing_date,'/', iso3,'_processed_contact_patterns.RData'))


###########################################################################
#### save array version after accounting for population size 

load(paste0('FILEPATH', iso3,'/',processing_date,'/', iso3,'_processed_population_by_epiweek_measles_model_age_bins_by_adm2.RData'))
pop_array <- pop_array_red

W <- array(NA, dim=c(24, 24, dim(pop_array)[1], dim(pop_array)[3]))

for (i in 2:dim(pop_array)[1]){
  message(i)
  ###########################################################################
  ### make matricies reciprical
  for(a in 1:24){
    for(b in 1:24){
      for(d in 1:dim(pop_array)[3]){
        W[a,b,i,d] <- (1 / pop_array[i,a,d])*(0.5) * ((contact_matrix[a,b]*pop_array[i,a,d]) + (contact_matrix[b,a]* pop_array[i,b,d]))
      }
    }
  }
  
}

save(W, file=paste0('FILEPATH', iso3, '/',processing_date,'/', iso3,'_processed_contact_patterns_array_district.RData'))

