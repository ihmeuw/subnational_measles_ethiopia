
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

source('FILEPATH/prep_functions.R')
source('FILEPATH/shapefile_functions.R')
source('FILEPATH/misc_functions.R')
source('FILEPATH/location_metadata_functions.R')
mbg_central/util/location_metadata_functions.R
# mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0('FILEPATH/misc_vaccine_functions.R'))
# library(kernlab)
library(data.table)
library(ggplot2)


source(paste0('FILEPATH/misc_vaccine_functions.R'))
library(readxl)

year_list <- 1980:2019
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
  


###################################################################################################### 
################# lets births by epiweek
###################################################################################################### 
births <- fread(paste0('FILEPATH/gbd_age_', 164, '/',region,'_admin_2.csv'))
births <- subset(births, ADM0_CODE == code)

births$births <- births$pop
births<-births[year %in% year_list]
births$births <- (births$births) / 53

births <- as.data.table(cbind(rep(births[,get("ADM2_CODE")], each = 53),cbind(rep(births$year, each = 53),rep(births$births, each = 53))))
colnames(births) <- c("ADM2_CODE","year", "births")
births[,births:=as.numeric(births)]

births_matrix <- matrix(0, 2279, nadm2s)
rownames(births_matrix) <- 1:2279
colnames(births_matrix) <- adm2_names

for(a in adm2_names){
  message(a)
  births_matrix[1:2120,as.character(a)] <- births[get('ADM2_CODE')==a]$births
}

births_matrix[2121:2279,] <- births_matrix[2120,]

#### save births data in a RData set
save(births_matrix,
     file = paste0('FILEPATH', iso3,'/', processing_date, '/', iso3,'_processed_births_by_adm2.RData'))


