

## clear environment
rm(list=ls())

library(readxl)

processing_date <- '2023_01_03'

iso3 <- 'ETH'

dir.create(paste0("FILEPATH", iso3) , recursive=T)
year_list <- 1980:2021

dir.create(paste0('FILEPATH',iso3,'/', processing_date))

###################################################################################################### 
################# pre stratified cases
###################################################################################################### 

data <- data.table(read_excel('FILEPATH/incidence_series_2020.xls', sheet='Measles'))
data <- subset(data, ISO_code == iso3)

t_new <- t(data)
t_new <- t_new[-c(1:4),]
new <- as.numeric(rev(as.vector(t_new)))
national_cases <- new

#### save national pre-stratified case data in a RData set
save(national_cases,
     file = paste0('FILEPATH',iso3,'/', processing_date,'/', iso3,'_processed_national_cases.RData'))


