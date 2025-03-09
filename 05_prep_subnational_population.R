
## clear environment
rm(list=ls())

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

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))

source('FILEPATH/prep_functions.R')
source('FILEPATH/shapefile_functions.R')
source('FILEPATH/misc_functions.R')
source('FILEPATHlocation_metadata_functions.R')
mbg_central/util/location_metadata_functions.R
# mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0('FILEPATH/misc_vaccine_functions.R'))
# library(kernlab)
library(data.table)
library(ggplot2)


# Custom load indicator-specific functions

source(paste0('FILEPATH/misc_vaccine_functions.R'))
library(readxl)

year_list <- 1980:2020
#processing_date <- make_time_stamp(time_stamp=T)
processing_date <- '2023_01_03'
iso3 <- 'ETH'


###################################################################################################### 
################# grab important age and location codes
###################################################################################################### 
source("FILEPATH/get_population.R")

#Table of regional locations and names
code <- get_adm0_codes(iso3)

loc <- get_location_code_mapping(shapefile_version='2020_05_21')
loc <- loc[ADM_CODE==code, list(country=ihme_lc_id, ADM0_CODE=ADM_CODE,gbd_id=loc_id)]

#Grab the region name
if(iso3 == 'ETH') region <- 'vax_essa' 

sample<-fread(paste0('FILEPATH/age_bin',1,'/',region,'_admin_2.csv'))
sample <- sample[ADM0_CODE==code & !is.na(pop)] #Drop admin 2s where we have NAs in the RI data (i.e. some wonky places in indonesia)
adm2_names <-sort(unique(sample$ADM2_CODE))
nadm2s = length(adm2_names)

#All the ages we will loop over
ages <- c(42, 388, 2:7, 53:142,235)

#Get age lookup table to correct age IDs
age_lookup <- fread('FILEPATH/vax_age_lookup.csv')

#Connect ages that span months/years to individual age-weeks
week_age_association <- as.data.table(cbind(c(1:5036), c(rep(42, 4), 
                                                         rep(388, 21), 
                                                         rep(2, 14), 
                                                         rep(3, 14), 
                                                         rep(c(4:7,53:142), each=53), 
                                                         235)))
names(week_age_association) <- c('age_week','age_group')

#using all gbd age names
week_age_association2 <- as.data.table(cbind(c(1:5036), c(rep(42, 4), rep(388, 21), rep(2, 14), rep(433, 14), rep(49:142, each=53), 235)))
names(week_age_association2) <- c('age_week','age_group')


#using all gbd age names
week_age_association_measles <- as.data.table(cbind(c(1:107), c(rep(42, 1), # 0 month olds
                                                                rep(388, 5), # 1 - 5 month olds
                                                                rep(2, 3), # 6 - 8 month olds
                                                                rep(3, 3), # 9 - 11 month olds
                                                                rep(c(4:7, 53:142)), # 1 year to 94 years old
                                                                235))) # 95+ year olds
names(week_age_association_measles) <- c('age_bin','age_group')


###################################################################################################### 
################# lets get / compute population by epiweek
###################################################################################################### 
message('working on population')
#Arrays to which all the weekly data will be assigned--age by weeks and regular ages
pop_array <- array(dim = c(2173,107,nadm2s),dimnames = list(1:2173 , 1:107,c(paste0('adm2_',adm2_names))))
# 
# pop_array_age_agg <- array(dim = c(length(ages),2120,nadm2s),dimnames = list(ages,1:2120 , c(paste0('adm2_',adm2_names))))


#loop over ages
for(age in ages){
  print(age)
  

    if(age %in% 2:7){ #Grab files from alyssa's custom age bins
      pop <- fread(paste0('FILEPATH/age_bin',age,'/',region,'_admin_2.csv'))
    }else{ #grab files with standard gbd values
      pop <- fread(paste0('FILEPATH/gbd_age_', age, '/',region,'_admin_2.csv'))
    }

  message("for age ", age, " <0 is ", length(which(pop$pop < 0)))
  
  
  pop[,V1:=NULL]
  pop <- pop[ADM0_CODE==loc$ADM0_CODE]
  pop <- pop[year %in% year_list]
  
  for(ad in adm2_names){
    #subset to population for one admin2
    pop_ad <- pop[ADM2_CODE==ad]$pop
    
    #Divvy up changes in population that occur over a year to occur gradually over the weeks
    diff_vec <- c()
    year_all <- c()
    year_part <- c()
    year_all_pop <- c()
    year_part_pop <- c()
    
    for(i in 1:length(year_list)){
      diff_vec[i] <- (pop_ad[i+1] - pop_ad[i]) / 53
      if(is.na(diff_vec[i])) diff_vec[i] <- diff_vec[i - 1]   #!#!#!#!#!hot fix to compensate for missing 2020 population data, used to calculate pop changes across year 2019
      
      year_part <- rep(diff_vec[i], 53)
      year_all <- c(year_all, year_part)
      year_part_pop <- rep(pop_ad[i], 53)
      year_all_pop <- c(year_all_pop, year_part_pop)
    }
    for(j in 2:(length(year_list) * 53)){
      year_all_pop[j] <- year_all_pop[j-1] + year_all[j]
    }
    pop_ad <- year_all_pop  
    
    
    
    #grab the age weeks for the relevant age year/month
    age_weeks <- week_age_association_measles[age_group==age]$age_bin
    
    #Apply to relevant subset of the population arrays  
    pop_array[1:(length(year_list)*53),age_weeks,paste0('adm2_', ad)] <- rep(pop_ad, each=length(age_weeks)) /length(age_weeks)  #!#!#!#Pop ad might NOT?? need to be divided by 53 or number of weeks in age group? reassess once set up
    
    #pop_array_age_agg[as.character(age),,paste0('adm2_', ad)] <- pop_ad
  }
  
}


pop_array_red <- array(NA, dim=c(2279, 24,nadm2s))

pop_array_red[1:2173,1:16,] = pop_array[1:2173,1:16,]
  
for(i in 1:2173){
  for(d in 1:nadm2s){
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
  pop_array_red[2174:2279,,] <- pop_array_red[2173,,]
}


#### save population data in a RData set
save(pop_array, pop_array_red, 
     file = paste0('FILEPATH', iso3, '/',processing_date,'/', iso3,'_processed_population_by_epiweek_measles_model_age_bins_by_adm2.RData'))


