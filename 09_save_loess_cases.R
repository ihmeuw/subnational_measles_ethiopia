


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
load(paste0('FILEPATH/', '2022_02_01','/', iso3,'_processed_cases.RData')) # subnational


cases_age_bin0_4_matrix

dim(cases_age_bin0_4_matrix)

loess(cases_age_bin0_4_matrix[,11], sm)

smoothed_cases_0_4 <- data.table()
smoothed_cases_5_9 <- data.table()
smoothed_cases_10_14 <- data.table()
smoothed_cases_15_24 <- data.table()
smoothed_cases_25_34 <- data.table()
smoothed_cases_35_44 <- data.table()
smoothed_cases_45_54 <- data.table()
smoothed_cases_55_64 <- data.table()
smoothed_cases_65_plus <- data.table()
for(d in 1:79){
  
  message(d)
  cases_df <- data.table()
  cases_df$cases <- cases_age_bin0_4_matrix[1750:2120,d]
  cases_df$week <- 1:371
  loessMod10 <- loess(cases ~ week, data=cases_df, span=0.20)
  smoothed_cases_0_4 <- cbind(smoothed_cases_0_4, c(predict(loessMod10)))
  
  
  cases_df <- data.table()
  cases_df$cases <- cases_age_bin5_9_matrix[1750:2120,d]
  cases_df$week <- 1:371
  loessMod10 <- loess(cases ~ week, data=cases_df, span=0.20)
  smoothed_cases_5_9 <- cbind(smoothed_cases_5_9, c(predict(loessMod10)))
  
  
  cases_df <- data.table()
  cases_df$cases <- cases_age_bin10_14_matrix[1750:2120,d]
  cases_df$week <- 1:371
  loessMod10 <- loess(cases ~ week, data=cases_df, span=0.20)
  smoothed_cases_10_14 <- cbind(smoothed_cases_10_14, c(predict(loessMod10)))
  
  
  cases_df <- data.table()
  cases_df$cases <- cases_age_bin15_19_matrix[1750:2120,d] + cases_age_bin20_24_matrix[1750:2120,d]
  cases_df$week <- 1:371
  loessMod10 <- loess(cases ~ week, data=cases_df, span=0.20)
  smoothed_cases_15_24 <- cbind(smoothed_cases_15_24, c(predict(loessMod10)))
  
  
  cases_df <- data.table()
  cases_df$cases <- cases_age_bin25_29_matrix[1750:2120,d] + cases_age_bin30_34_matrix[1750:2120,d]
  cases_df$week <- 1:371
  loessMod10 <- loess(cases ~ week, data=cases_df, span=0.20)
  smoothed_cases_25_34 <- cbind(smoothed_cases_25_34, c(predict(loessMod10)))
  
  
  cases_df <- data.table()
  cases_df$cases <- cases_age_bin35_39_matrix[1750:2120,d] + cases_age_bin40_44_matrix[1750:2120,d]
  cases_df$week <- 1:371
  loessMod10 <- loess(cases ~ week, data=cases_df, span=0.20)
  smoothed_cases_35_44 <- cbind(smoothed_cases_35_44, c(predict(loessMod10)))
  
  
  cases_df <- data.table()
  cases_df$cases <- cases_age_bin45_49_matrix[1750:2120,d] + cases_age_bin50_54_matrix[1750:2120,d]
  cases_df$week <- 1:371
  loessMod10 <- loess(cases ~ week, data=cases_df, span=0.20)
  smoothed_cases_45_54 <- cbind(smoothed_cases_45_54, c(predict(loessMod10)))
  
  
  cases_df <- data.table()
  cases_df$cases <- cases_age_bin55_59_matrix[1750:2120,d] + cases_age_bin60_64_matrix[1750:2120,d]
  cases_df$week <- 1:371
  loessMod10 <- loess(cases ~ week, data=cases_df, span=0.20)
  smoothed_cases_55_64 <- cbind(smoothed_cases_55_64, c(predict(loessMod10)))
  
  
  cases_df <- data.table()
  cases_df$cases <- cases_age_bin65_69_matrix[1750:2120,d] + cases_age_bin70_74_matrix[1750:2120,d] + 
                      cases_age_bin75_79_matrix[1750:2120,d] + cases_age_bin80_84_matrix[1750:2120,d]
  cases_df$week <- 1:371
  loessMod10 <- loess(cases ~ week, data=cases_df, span=0.20)
  smoothed_cases_65_plus <- cbind(smoothed_cases_65_plus, c(predict(loessMod10)))

  
}

filler <- matrix(NA, nrow = 1749, ncol = 79)

cases_age_bin0_4_matrix <- rbind(filler, smoothed_cases_0_4, use.names = F)
cases_age_bin5_9_matrix <- rbind(filler, smoothed_cases_5_9, use.names = F)
cases_age_bin10_14_matrix <- rbind(filler, smoothed_cases_10_14, use.names = F)
cases_age_bin15_24_matrix <- rbind(filler, smoothed_cases_15_24, use.names = F)
cases_age_bin25_34_matrix <- rbind(filler, smoothed_cases_25_34, use.names = F)
cases_age_bin35_44_matrix <- rbind(filler, smoothed_cases_35_44, use.names = F)
cases_age_bin45_54_matrix <- rbind(filler, smoothed_cases_45_54, use.names = F)
cases_age_bin55_64_matrix <- rbind(filler, smoothed_cases_55_64, use.names = F)
cases_age_bin65_plus_matrix <- rbind(filler, smoothed_cases_65_plus, use.names = F)


save(cases_age_bin0_4_matrix, cases_age_bin5_9_matrix, cases_age_bin10_14_matrix,
     cases_age_bin15_24_matrix, cases_age_bin25_34_matrix, cases_age_bin35_44_matrix,
     cases_age_bin45_54_matrix, cases_age_bin55_64_matrix, cases_age_bin65_plus_matrix,
     file = paste0('FILEPATH/', '2022_02_01','/', iso3,'_processed_cases_loess.RData')
     )





iso3 <- 'ETH'



load(
  file = paste0('FILEPATH/', '2022_02_01','/', iso3,'_processed_cases_loess.RData')
)

cases_age_bin0_4_matrix[cases_age_bin0_4_matrix < 0] <- 0
cases_age_bin0_4_matrix_smooth <- cases_age_bin0_4_matrix

load(
  file = paste0('FILEPATH/', '2022_02_01','/', iso3,'_processed_cases.RData')
)



matplot(cases_age_bin0_4_matrix[1751:2120, 1], type='l')

test <- melt(cases_age_bin0_4_matrix[1751:2120,])
colnames(test) <- c('week', 'district', 'value')
test$class <- 'reported'

cases_age_bin0_4_matrix_smooth <- as.matrix(cases_age_bin0_4_matrix_smooth)
colnames(cases_age_bin0_4_matrix_smooth) <- 1:79
test2 <- melt(as.matrix(cases_age_bin0_4_matrix_smooth[1751:2120,]))
colnames(test2) <- c('week', 'district', 'value')
test2$class <- 'smoothed'


df <- rbind(test, test2)


ggplot(data=subset(df, district %in% c(1,3,4,8))) + geom_line(aes(x=week, y=value, color=as.factor(district))) + facet_wrap(~class+district, nrow=2, scales = 'free')
