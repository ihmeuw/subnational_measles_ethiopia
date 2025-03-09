
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
source(paste0('FILEPATH/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0('FILEPATH/misc_vaccine_functions.R'))
source(paste0('FILEPATH/gravity_functions.R'))


## load iso3 and shapefile version
iso3 <- 'ETH'
processing_date <- '2023_01_03'

library(plyr)


full_mat <- array(NA, dim=c(2279, 107, nadmin_units))
#full_mat_too_add <- array(NA, dim=c(53, 107, 79))

for(i in 1980:2019)  {
  message(i)

  eth <- readRDS(paste0('FILEPATH/cohort_counts_for_year_',i,'_',iso3,'.RDS'))

  for(a in unique(eth$age_group_id)){
    #message(a)
    eth_sub <- subset(eth, val == "n_ri_sia" & age_group_id == a & dose == "mcv1")
    
    for(ad in unique(eth_sub$ADM2_CODE)){
      #message(ad)
      eth_sub2 <- subset(eth_sub, ADM2_CODE == ad)
      
      index_ad <- which( unique(eth_sub$ADM2_CODE) == ad)
    
      full_mat[((i-1980 )*53+ 1):((i+1-1980 )*53), which(eth$age_group_id == a)[1],index_ad] <- eth_sub2$cov
    }
  }    
}


ri_sia_time_series_mcv1 <- full_mat


full_mat <- array(NA, dim=c(2279, 107, nadmin_units))

for(i in 1980:2019)  {
  message(i)
  eth <- readRDS(paste0('FILEPATH/cohort_counts_for_year_',i,'_',iso3,'.RDS'))
  
  for(a in unique(eth$age_group_id)){
    #message(a)
    eth_sub <- subset(eth, val == "n_ri_sia" & age_group_id == a & dose == "mcv2")
    
    for(ad in unique(eth_sub$ADM2_CODE)){
      #message(ad)
      eth_sub2 <- subset(eth_sub, ADM2_CODE == ad)
      
      index_ad <- which( unique(eth_sub$ADM2_CODE) == ad)

      full_mat[((i-1980 )*53+ 1):((i+1-1980 )*53), which(eth$age_group_id == a)[1],index_ad] <- eth_sub2$cov
    }
  }    
}


ri_sia_time_series_mcv2 <- full_mat

save(ri_sia_time_series_mcv1, ri_sia_time_series_mcv2, file = paste0('FILEPATH/', iso3,'/',processing_date,'/', iso3,'_processed_coverage.RData'))
