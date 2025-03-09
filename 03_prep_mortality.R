

## clear environment
rm(list=ls())

library(readxl)

processing_date <- '2023_01_03'

iso3 <- 'ETH'

dir.create(paste0("FILEPATH", iso3) , recursive=T)
year_list <- 1980:2022

dir.create(paste0('/FILEPATH',iso3,'/', processing_date))


if(iso3 == "ETH") ihme_loc_id <- 179

###################################################################################################### 
################# pull mortality estimates
###################################################################################################### 
source("FILEPATH/get_envelope.R")

mortality_master <- get_envelope(year_id=year_list, sex_id=3, location_id = ihme_loc_id, 
                                 gbd_round_id=7, decomp_step='step3', with_hiv=1, with_shock=1, rates=1, age_group_id = -1)

mortality_age_bin0 <- subset(mortality_master, age_group_id == 28)
mortality_age_bin0 <- mortality_age_bin0$mean

mortality_age_bin1 <- subset(mortality_master, age_group_id == 238)
mortality_age_bin1 <- mortality_age_bin1$mean

mortality_age_bin2 <- subset(mortality_master, age_group_id == 5)
mortality_age_bin2 <- mortality_age_bin2$mean

mortality_age_bin3 <- subset(mortality_master, age_group_id == 5)
mortality_age_bin3 <- mortality_age_bin3$mean

mortality_age_bin4 <- subset(mortality_master, age_group_id == 5)
mortality_age_bin4 <- mortality_age_bin4$mean

mortality_age_bin5 <- subset(mortality_master, age_group_id == 6)
mortality_age_bin5 <- mortality_age_bin5$mean

mortality_age_bin6 <- subset(mortality_master, age_group_id == 6)
mortality_age_bin6 <- mortality_age_bin6$mean

mortality_age_bin7 <- subset(mortality_master, age_group_id == 6)
mortality_age_bin7 <- mortality_age_bin7$mean

mortality_age_bin8 <- subset(mortality_master, age_group_id == 6)
mortality_age_bin8 <- mortality_age_bin8$mean

mortality_age_bin9 <- subset(mortality_master, age_group_id == 6)
mortality_age_bin9 <- mortality_age_bin9$mean

mortality_age_bin10 <- subset(mortality_master, age_group_id == 7)
mortality_age_bin10 <- mortality_age_bin10$mean

mortality_age_bin11 <- subset(mortality_master, age_group_id == 7)
mortality_age_bin11 <- mortality_age_bin11$mean

mortality_age_bin12 <- subset(mortality_master, age_group_id == 7)
mortality_age_bin12 <- mortality_age_bin12$mean

mortality_age_bin13 <- subset(mortality_master, age_group_id == 7)
mortality_age_bin13 <- mortality_age_bin13$mean

mortality_age_bin14 <- subset(mortality_master, age_group_id == 7)
mortality_age_bin14 <- mortality_age_bin14$mean

mortality_age_bin15 <- subset(mortality_master, age_group_id == 8)
mortality_age_bin15 <- mortality_age_bin15$mean

mortality_age_bin16 <- subset(mortality_master, age_group_id == 8)
mortality_age_bin16 <- mortality_age_bin16$mean

mortality_age_bin17 <- subset(mortality_master, age_group_id == 8)
mortality_age_bin17 <- mortality_age_bin17$mean

mortality_age_bin18 <- subset(mortality_master, age_group_id == 8)
mortality_age_bin18 <- mortality_age_bin18$mean

mortality_age_bin19 <- subset(mortality_master, age_group_id == 8)
mortality_age_bin19 <- mortality_age_bin19$mean

mortality_age_bin20 <- subset(mortality_master, age_group_id == 9)
mortality_age_bin20 <- mortality_age_bin20$mean

mortality_age_bin21 <- subset(mortality_master, age_group_id == 9)
mortality_age_bin21 <- mortality_age_bin21$mean

mortality_age_bin22 <- subset(mortality_master, age_group_id == 9)
mortality_age_bin22 <- mortality_age_bin22$mean

mortality_age_bin23 <- subset(mortality_master, age_group_id == 9)
mortality_age_bin23 <- mortality_age_bin23$mean

mortality_age_bin24 <- subset(mortality_master, age_group_id == 9)
mortality_age_bin24 <- mortality_age_bin24$mean

mortality_age_bin25 <- subset(mortality_master, age_group_id == 10)
mortality_age_bin25 <- mortality_age_bin25$mean

mortality_age_bin26 <- subset(mortality_master, age_group_id == 10)
mortality_age_bin26 <- mortality_age_bin26$mean

mortality_age_bin27 <- subset(mortality_master, age_group_id == 10)
mortality_age_bin27 <- mortality_age_bin27$mean

mortality_age_bin28 <- subset(mortality_master, age_group_id == 10)
mortality_age_bin28 <- mortality_age_bin28$mean

mortality_age_bin29 <- subset(mortality_master, age_group_id == 10)
mortality_age_bin29 <- mortality_age_bin29$mean

mortality_age_bin30 <- subset(mortality_master, age_group_id == 11)
mortality_age_bin30 <- mortality_age_bin30$mean

mortality_age_bin31 <- subset(mortality_master, age_group_id == 11)
mortality_age_bin31 <- mortality_age_bin31$mean

mortality_age_bin32 <- subset(mortality_master, age_group_id == 11)
mortality_age_bin32 <- mortality_age_bin32$mean

mortality_age_bin33 <- subset(mortality_master, age_group_id == 11)
mortality_age_bin33 <- mortality_age_bin33$mean

mortality_age_bin34 <- subset(mortality_master, age_group_id == 11)
mortality_age_bin34 <- mortality_age_bin34$mean

mortality_age_bin35 <- subset(mortality_master, age_group_id == 12)
mortality_age_bin35 <- mortality_age_bin35$mean

mortality_age_bin36 <- subset(mortality_master, age_group_id == 12)
mortality_age_bin36 <- mortality_age_bin36$mean

mortality_age_bin37 <- subset(mortality_master, age_group_id == 12)
mortality_age_bin37 <- mortality_age_bin37$mean

mortality_age_bin38 <- subset(mortality_master, age_group_id == 12)
mortality_age_bin38 <- mortality_age_bin38$mean

mortality_age_bin39 <- subset(mortality_master, age_group_id == 12)
mortality_age_bin39 <- mortality_age_bin39$mean

mortality_age_bin40 <- subset(mortality_master, age_group_id == 13)
mortality_age_bin40 <- mortality_age_bin40$mean

mortality_age_bin41 <- subset(mortality_master, age_group_id == 13)
mortality_age_bin41 <- mortality_age_bin41$mean

mortality_age_bin42 <- subset(mortality_master, age_group_id == 13)
mortality_age_bin42 <- mortality_age_bin42$mean

mortality_age_bin43 <- subset(mortality_master, age_group_id == 13)
mortality_age_bin43 <- mortality_age_bin43$mean

mortality_age_bin44 <- subset(mortality_master, age_group_id == 13)
mortality_age_bin44 <- mortality_age_bin44$mean

mortality_age_bin45 <- subset(mortality_master, age_group_id == 14)
mortality_age_bin45 <- mortality_age_bin45$mean

mortality_age_bin46 <- subset(mortality_master, age_group_id == 14)
mortality_age_bin46 <- mortality_age_bin46$mean

mortality_age_bin47 <- subset(mortality_master, age_group_id == 14)
mortality_age_bin47 <- mortality_age_bin47$mean

mortality_age_bin48 <- subset(mortality_master, age_group_id == 14)
mortality_age_bin48 <- mortality_age_bin48$mean

mortality_age_bin49 <- subset(mortality_master, age_group_id == 14)
mortality_age_bin49 <- mortality_age_bin49$mean

mortality_age_bin50 <- subset(mortality_master, age_group_id == 15)
mortality_age_bin50 <- mortality_age_bin50$mean

mortality_age_bin51 <- subset(mortality_master, age_group_id == 15)
mortality_age_bin51 <- mortality_age_bin51$mean

mortality_age_bin52 <- subset(mortality_master, age_group_id == 15)
mortality_age_bin52 <- mortality_age_bin52$mean

mortality_age_bin53 <- subset(mortality_master, age_group_id == 15)
mortality_age_bin53 <- mortality_age_bin53$mean

mortality_age_bin54 <- subset(mortality_master, age_group_id == 15)
mortality_age_bin54 <- mortality_age_bin54$mean

mortality_age_bin55 <- subset(mortality_master, age_group_id == 16)
mortality_age_bin55 <- mortality_age_bin55$mean

mortality_age_bin56 <- subset(mortality_master, age_group_id == 16)
mortality_age_bin56 <- mortality_age_bin56$mean

mortality_age_bin57 <- subset(mortality_master, age_group_id == 16)
mortality_age_bin57 <- mortality_age_bin57$mean

mortality_age_bin58 <- subset(mortality_master, age_group_id == 16)
mortality_age_bin58 <- mortality_age_bin58$mean

mortality_age_bin59 <- subset(mortality_master, age_group_id == 16)
mortality_age_bin59 <- mortality_age_bin59$mean

mortality_age_bin60 <- subset(mortality_master, age_group_id == 17)
mortality_age_bin60 <- mortality_age_bin60$mean

mortality_age_bin61 <- subset(mortality_master, age_group_id == 17)
mortality_age_bin61 <- mortality_age_bin61$mean

mortality_age_bin62 <- subset(mortality_master, age_group_id == 17)
mortality_age_bin62 <- mortality_age_bin62$mean

mortality_age_bin63 <- subset(mortality_master, age_group_id == 17)
mortality_age_bin63 <- mortality_age_bin63$mean

mortality_age_bin64 <- subset(mortality_master, age_group_id == 17)
mortality_age_bin64 <- mortality_age_bin64$mean

mortality_age_bin65 <- subset(mortality_master, age_group_id == 18)
mortality_age_bin65 <- mortality_age_bin65$mean

mortality_age_bin66 <- subset(mortality_master, age_group_id == 18)
mortality_age_bin66 <- mortality_age_bin66$mean

mortality_age_bin67 <- subset(mortality_master, age_group_id == 18)
mortality_age_bin67 <- mortality_age_bin67$mean

mortality_age_bin68 <- subset(mortality_master, age_group_id == 18)
mortality_age_bin68 <- mortality_age_bin68$mean

mortality_age_bin69 <- subset(mortality_master, age_group_id == 18)
mortality_age_bin69 <- mortality_age_bin69$mean

mortality_age_bin70 <- subset(mortality_master, age_group_id == 19)
mortality_age_bin70 <- mortality_age_bin70$mean

mortality_age_bin71 <- subset(mortality_master, age_group_id == 19)
mortality_age_bin71 <- mortality_age_bin71$mean

mortality_age_bin72 <- subset(mortality_master, age_group_id == 19)
mortality_age_bin72 <- mortality_age_bin72$mean

mortality_age_bin73 <- subset(mortality_master, age_group_id == 19)
mortality_age_bin73 <- mortality_age_bin73$mean

mortality_age_bin74 <- subset(mortality_master, age_group_id == 19)
mortality_age_bin74 <- mortality_age_bin74$mean

mortality_age_bin75 <- subset(mortality_master, age_group_id == 20)
mortality_age_bin75 <- mortality_age_bin75$mean

mortality_age_bin76 <- subset(mortality_master, age_group_id == 20)
mortality_age_bin76 <- mortality_age_bin76$mean

mortality_age_bin77 <- subset(mortality_master, age_group_id == 20)
mortality_age_bin77 <- mortality_age_bin77$mean

mortality_age_bin78 <- subset(mortality_master, age_group_id == 20)
mortality_age_bin78 <- mortality_age_bin78$mean

mortality_age_bin79 <- subset(mortality_master, age_group_id == 20)
mortality_age_bin79 <- mortality_age_bin79$mean

mortality_age_bin80 <- subset(mortality_master, age_group_id == 30)
mortality_age_bin80 <- mortality_age_bin80$mean

mortality_age_bin81 <- subset(mortality_master, age_group_id == 30)
mortality_age_bin81 <- mortality_age_bin81$mean

mortality_age_bin82 <- subset(mortality_master, age_group_id == 30)
mortality_age_bin82 <- mortality_age_bin82$mean

mortality_age_bin83 <- subset(mortality_master, age_group_id == 30)
mortality_age_bin83 <- mortality_age_bin83$mean

mortality_age_bin84 <- subset(mortality_master, age_group_id == 30)
mortality_age_bin84 <- mortality_age_bin84$mean

mortality_age_bin85 <- subset(mortality_master, age_group_id == 31)
mortality_age_bin85 <- mortality_age_bin85$mean

mortality_age_bin86 <- subset(mortality_master, age_group_id == 31)
mortality_age_bin86 <- mortality_age_bin86$mean

mortality_age_bin87 <- subset(mortality_master, age_group_id == 31)
mortality_age_bin87 <- mortality_age_bin87$mean

mortality_age_bin88 <- subset(mortality_master, age_group_id == 31)
mortality_age_bin88 <- mortality_age_bin88$mean

mortality_age_bin89 <- subset(mortality_master, age_group_id == 31)
mortality_age_bin89 <- mortality_age_bin89$mean

mortality_age_bin90 <- subset(mortality_master, age_group_id == 32)
mortality_age_bin90 <- mortality_age_bin90$mean

mortality_age_bin91 <- subset(mortality_master, age_group_id == 32)
mortality_age_bin91 <- mortality_age_bin91$mean

mortality_age_bin92 <- subset(mortality_master, age_group_id == 32)
mortality_age_bin92 <- mortality_age_bin92$mean

mortality_age_bin93 <- subset(mortality_master, age_group_id == 32)
mortality_age_bin93 <- mortality_age_bin93$mean

mortality_age_bin94 <- subset(mortality_master, age_group_id == 32)
mortality_age_bin94 <- mortality_age_bin94$mean

mortality_age_bin95 <- subset(mortality_master, age_group_id == 235 )
mortality_age_bin95 <- mortality_age_bin95$mean



mortality_age_bin_month_0 <- subset(mortality_master, age_group_id == 42)
mortality_age_bin_month_0 <- mortality_age_bin_month_0$mean / 12

mortality_age_bin_month_1 <- subset(mortality_master, age_group_id == 388)
mortality_age_bin_month_1 <- mortality_age_bin_month_1$mean / 12 

mortality_age_bin_month_2 <- subset(mortality_master, age_group_id == 388)
mortality_age_bin_month_2 <- mortality_age_bin_month_2$mean / 12

mortality_age_bin_month_3 <- subset(mortality_master, age_group_id == 388)
mortality_age_bin_month_3 <- mortality_age_bin_month_3$mean / 12

mortality_age_bin_month_4 <- subset(mortality_master, age_group_id == 388)
mortality_age_bin_month_4 <- mortality_age_bin_month_4$mean / 12

mortality_age_bin_month_5 <- subset(mortality_master, age_group_id == 388)
mortality_age_bin_month_5 <- mortality_age_bin_month_5$mean / 12

mortality_age_bin_month_6 <- subset(mortality_master, age_group_id == 389)
mortality_age_bin_month_6 <- mortality_age_bin_month_6$mean / 12

mortality_age_bin_month_7 <- subset(mortality_master, age_group_id == 389)
mortality_age_bin_month_7 <- mortality_age_bin_month_7$mean/ 12

mortality_age_bin_month_8 <- subset(mortality_master, age_group_id == 389)
mortality_age_bin_month_8 <- mortality_age_bin_month_8$mean/ 12

mortality_age_bin_month_9 <- subset(mortality_master, age_group_id == 389)
mortality_age_bin_month_9 <- mortality_age_bin_month_9$mean / 12

mortality_age_bin_month_10 <- subset(mortality_master, age_group_id == 389)
mortality_age_bin_month_10 <- mortality_age_bin_month_10$mean / 12

mortality_age_bin_month_11 <- subset(mortality_master, age_group_id == 389)
mortality_age_bin_month_11 <- mortality_age_bin_month_11$mean / 12





mortality_matrix <- cbind(mortality_age_bin_month_0, mortality_age_bin_month_1, mortality_age_bin_month_2, mortality_age_bin_month_3, 
                          mortality_age_bin_month_4, mortality_age_bin_month_5, mortality_age_bin_month_6, mortality_age_bin_month_7, 
                          mortality_age_bin_month_8, mortality_age_bin_month_9, mortality_age_bin_month_10, mortality_age_bin_month_11, 
                          
                          mortality_age_bin1, mortality_age_bin2, mortality_age_bin3, mortality_age_bin4,
                          
                          mortality_age_bin5, mortality_age_bin6, mortality_age_bin7, mortality_age_bin8, mortality_age_bin9, 
                          mortality_age_bin10, mortality_age_bin11, mortality_age_bin12, mortality_age_bin13, mortality_age_bin14,
                          mortality_age_bin15, mortality_age_bin16, mortality_age_bin17, mortality_age_bin18, mortality_age_bin19, 
                          mortality_age_bin20, mortality_age_bin21, mortality_age_bin22, mortality_age_bin23, mortality_age_bin24,
                          mortality_age_bin25, mortality_age_bin26, mortality_age_bin27, mortality_age_bin28, mortality_age_bin29, 
                          mortality_age_bin30, mortality_age_bin31, mortality_age_bin32, mortality_age_bin33, mortality_age_bin34,
                          mortality_age_bin35, mortality_age_bin36, mortality_age_bin37, mortality_age_bin38, mortality_age_bin39, 
                          mortality_age_bin40, mortality_age_bin41, mortality_age_bin42, mortality_age_bin43, mortality_age_bin44,
                          mortality_age_bin45, mortality_age_bin46, mortality_age_bin47, mortality_age_bin48, mortality_age_bin49, 
                          mortality_age_bin50, mortality_age_bin51, mortality_age_bin52, mortality_age_bin53, mortality_age_bin54,
                          mortality_age_bin55, mortality_age_bin56, mortality_age_bin57, mortality_age_bin58, mortality_age_bin59, 
                          mortality_age_bin60, mortality_age_bin61, mortality_age_bin62, mortality_age_bin63, mortality_age_bin64,
                          mortality_age_bin65, mortality_age_bin66, mortality_age_bin67, mortality_age_bin68, mortality_age_bin69, 
                          mortality_age_bin70, mortality_age_bin71, mortality_age_bin72, mortality_age_bin73, mortality_age_bin74,
                          mortality_age_bin75, mortality_age_bin76, mortality_age_bin77, mortality_age_bin78, mortality_age_bin79, 
                          mortality_age_bin80, mortality_age_bin81, mortality_age_bin82, mortality_age_bin83, mortality_age_bin84,
                          mortality_age_bin85, mortality_age_bin86, mortality_age_bin87, mortality_age_bin88, mortality_age_bin89, 
                          mortality_age_bin90, mortality_age_bin91, mortality_age_bin92, mortality_age_bin93, mortality_age_bin94,
                          mortality_age_bin95)



rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
full_mat <- data.table()
for(i in 1:43){
  temp <- rep.row(mortality_matrix[i,],53)
  full_mat <- rbind(full_mat, temp)
}

mortality_matrix <- copy(full_mat)

save(mortality_matrix ,
     file = paste0('FILEPATH',iso3,'/', processing_date,'/', iso3,'_processed_mortality.RData'))





mortality_matrix_weekly <- matrix(mortality_age_bin_month_0, 43, 4)
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin_month_1, 43, 4))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin_month_2, 43, 4))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin_month_3, 43, 4))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin_month_4, 43, 4))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin_month_5, 43, 5))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin_month_6, 43, 4))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin_month_7, 43, 5))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin_month_8, 43, 5))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin_month_9, 43, 4))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin_month_10, 43, 4))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin_month_11, 43, 5))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin1, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin2, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin3, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin4, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin5, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin6, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin7, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin8, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin9, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin10, 43, 53))

mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin11, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin12, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin13, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin14, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin15, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin16, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin17, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin18, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin19, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin20, 43, 53))

mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin21, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin22, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin23, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin24, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin25, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin26, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin27, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin28, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin29, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin30, 43, 53))

mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin31, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin32, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin33, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin34, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin35, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin36, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin37, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin38, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin39, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin40, 43, 53))

mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin41, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin43, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin43, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin44, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin45, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin46, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin47, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin48, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin49, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin50, 43, 53))

mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin51, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin52, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin53, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin54, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin55, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin56, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin57, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin58, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin59, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin60, 43, 53))

mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin61, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin62, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin63, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin64, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin65, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin66, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin67, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin68, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin69, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin70, 43, 53))

mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin71, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin72, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin73, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin74, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin75, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin76, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin77, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin78, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin79, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin80, 43, 53))

mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin81, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin82, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin83, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin84, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin85, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin86, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin87, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin88, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin89, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin90, 43, 53))

mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin91, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin92, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin93, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin94, 43, 53))
mortality_matrix_weekly <- cbind(mortality_matrix_weekly, matrix(mortality_age_bin95, 43, 1))

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
full_mat <- data.table()
for(i in 1:43){
  temp <- rep.row(mortality_matrix_weekly[i,],53)
  full_mat <- rbind(full_mat, temp)
}

mortality_matrix_weekly <- copy(full_mat)


save(mortality_matrix_weekly ,
     file = paste0('FILEPATH', iso3,'/', processing_date, '/', iso3,'_processed_mortality.RData'))




source("FILEPATH/get_envelope.R")

mortality_master <- get_envelope(year_id=year_list, sex_id=3, location_id = ihme_loc_id,
                                 gbd_round_id=7, decomp_step='step3', with_hiv=1, rates=1)

mortality_aggregated_all_age <- mortality_master$mean
save(mortality_aggregated_all_age,
     file = paste0('FILEPATH', iso3,'/',processing_date,'/', iso3,'_processed_mortality_all_age.RData'))


