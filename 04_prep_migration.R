

## clear environment
rm(list=ls())

library(readxl)

processing_date <- '2023_01_03'

iso3 <- 'ETH'

dir.create(paste0("FILEPATH", iso3) , recursive=T)
year_list <- 1980:2022

dir.create(paste0('FILEPATH/',iso3,'/', processing_date))

if(iso3 == "ETH") ihme_loc_id <- 179

############################################################################################################
## migration
############################################################################################################

migration <- fread('FILEPATH/net_migration_single_year_264.csv')
migration <- subset(migration, sex_id == 3)
migration <- subset(migration, location_id == ihme_loc_id)
migration <- subset(migration, year_id %in% year_list)

migration_rates <- subset(migration, measure_id == 18)

migration_age_bin0 <- subset(migration_rates, age_group_id == 28)
migration_age_bin0 <- migration_age_bin0$mean
migration_age_bin0[43] <- migration_age_bin0[42]

migration_age_bin_month_0 <- migration_age_bin_month_1 <- migration_age_bin_month_2 <- migration_age_bin_month_3 <- 
  migration_age_bin_month_4 <- migration_age_bin_month_5 <- migration_age_bin_month_6 <- migration_age_bin_month_7 <- 
  migration_age_bin_month_8 <- migration_age_bin_month_9 <- migration_age_bin_month_10 <- migration_age_bin_month_11 <- migration_age_bin0

migration_age_bin1 <- subset(migration_rates, age_group_id == 238)
migration_age_bin1 <- migration_age_bin1$mean
migration_age_bin1[43] <- migration_age_bin1[42]

migration_age_bin2 <- subset(migration_rates, age_group_id == 50)
migration_age_bin2 <- migration_age_bin2$mean
migration_age_bin2[43] <- migration_age_bin2[42]

migration_age_bin3 <- subset(migration_rates, age_group_id == 51)
migration_age_bin3 <- migration_age_bin3$mean
migration_age_bin3[43] <- migration_age_bin3[42]

migration_age_bin4 <- subset(migration_rates, age_group_id == 52)
migration_age_bin4 <- migration_age_bin4$mean
migration_age_bin4[43] <- migration_age_bin4[42]

migration_age_bin5 <- subset(migration_rates, age_group_id == 53)
migration_age_bin5 <- migration_age_bin5$mean
migration_age_bin5[43] <- migration_age_bin5[42]

migration_age_bin6 <- subset(migration_rates, age_group_id == 54)
migration_age_bin6 <- migration_age_bin6$mean
migration_age_bin6[43] <- migration_age_bin6[42]

migration_age_bin7 <- subset(migration_rates, age_group_id == 55)
migration_age_bin7 <- migration_age_bin7$mean
migration_age_bin7[43] <- migration_age_bin7[42]

migration_age_bin8 <- subset(migration_rates, age_group_id == 56)
migration_age_bin8 <- migration_age_bin8$mean
migration_age_bin8[43] <- migration_age_bin8[42]

migration_age_bin9 <- subset(migration_rates, age_group_id == 57)
migration_age_bin9 <- migration_age_bin9$mean
migration_age_bin9[43] <- migration_age_bin9[42]

migration_age_bin10 <- subset(migration_rates, age_group_id == 58)
migration_age_bin10 <- migration_age_bin10$mean
migration_age_bin10[43] <- migration_age_bin10[42]

migration_age_bin11 <- subset(migration_rates, age_group_id == 59)
migration_age_bin11 <- migration_age_bin11$mean
migration_age_bin11[43] <- migration_age_bin11[42]

migration_age_bin12 <- subset(migration_rates, age_group_id == 60)
migration_age_bin12 <- migration_age_bin12$mean
migration_age_bin12[43] <- migration_age_bin12[42]

migration_age_bin13 <- subset(migration_rates, age_group_id == 61)
migration_age_bin13 <- migration_age_bin13$mean
migration_age_bin13[43] <- migration_age_bin13[42]

migration_age_bin14 <- subset(migration_rates, age_group_id == 62)
migration_age_bin14 <- migration_age_bin14$mean
migration_age_bin14[43] <- migration_age_bin14[42]

migration_age_bin15 <- subset(migration_rates, age_group_id == 63)
migration_age_bin15 <- migration_age_bin15$mean
migration_age_bin15[43] <- migration_age_bin15[42]

migration_age_bin16 <- subset(migration_rates, age_group_id == 64)
migration_age_bin16 <- migration_age_bin16$mean
migration_age_bin16[43] <- migration_age_bin16[42]

migration_age_bin17 <- subset(migration_rates, age_group_id == 65)
migration_age_bin17 <- migration_age_bin17$mean
migration_age_bin17[43] <- migration_age_bin17[42]

migration_age_bin18 <- subset(migration_rates, age_group_id == 66)
migration_age_bin18 <- migration_age_bin18$mean
migration_age_bin18[43] <- migration_age_bin18[42]

migration_age_bin19 <- subset(migration_rates, age_group_id == 67)
migration_age_bin19 <- migration_age_bin19$mean
migration_age_bin19[43] <- migration_age_bin19[42]

migration_age_bin20 <- subset(migration_rates, age_group_id == 68)
migration_age_bin20 <- migration_age_bin20$mean
migration_age_bin20[43] <- migration_age_bin20[42]

migration_age_bin21 <- subset(migration_rates, age_group_id == 69)
migration_age_bin21 <- migration_age_bin21$mean
migration_age_bin21[43] <- migration_age_bin21[42]

migration_age_bin22 <- subset(migration_rates, age_group_id == 70)
migration_age_bin22 <- migration_age_bin22$mean
migration_age_bin22[43] <- migration_age_bin22[42]

migration_age_bin23 <- subset(migration_rates, age_group_id == 71)
migration_age_bin23 <- migration_age_bin23$mean
migration_age_bin23[43] <- migration_age_bin23[42]

migration_age_bin24 <- subset(migration_rates, age_group_id == 72)
migration_age_bin24 <- migration_age_bin24$mean
migration_age_bin24[43] <- migration_age_bin24[42]

migration_age_bin25 <- subset(migration_rates, age_group_id == 73)
migration_age_bin25 <- migration_age_bin25$mean
migration_age_bin25[43] <- migration_age_bin25[42]

migration_age_bin26 <- subset(migration_rates, age_group_id == 74)
migration_age_bin26 <- migration_age_bin26$mean
migration_age_bin26[43] <- migration_age_bin26[42]

migration_age_bin27 <- subset(migration_rates, age_group_id == 75)
migration_age_bin27 <- migration_age_bin27$mean
migration_age_bin27[43] <- migration_age_bin27[42]

migration_age_bin28 <- subset(migration_rates, age_group_id == 76)
migration_age_bin28 <- migration_age_bin28$mean
migration_age_bin28[43] <- migration_age_bin28[42]

migration_age_bin29 <- subset(migration_rates, age_group_id == 77)
migration_age_bin29 <- migration_age_bin29$mean
migration_age_bin29[43] <- migration_age_bin29[42]

migration_age_bin30 <- subset(migration_rates, age_group_id == 78)
migration_age_bin30 <- migration_age_bin30$mean
migration_age_bin30[43] <- migration_age_bin30[42]

migration_age_bin31 <- subset(migration_rates, age_group_id == 79)
migration_age_bin31 <- migration_age_bin31$mean
migration_age_bin31[43] <- migration_age_bin31[42]

migration_age_bin32 <- subset(migration_rates, age_group_id == 80)
migration_age_bin32 <- migration_age_bin32$mean
migration_age_bin32[43] <- migration_age_bin32[42]

migration_age_bin33 <- subset(migration_rates, age_group_id == 81)
migration_age_bin33 <- migration_age_bin33$mean
migration_age_bin33[43] <- migration_age_bin33[42]

migration_age_bin34 <- subset(migration_rates, age_group_id == 82)
migration_age_bin34 <- migration_age_bin34$mean
migration_age_bin34[43] <- migration_age_bin34[42]

migration_age_bin35 <- subset(migration_rates, age_group_id == 83)
migration_age_bin35 <- migration_age_bin35$mean
migration_age_bin35[43] <- migration_age_bin35[42]

migration_age_bin36 <- subset(migration_rates, age_group_id == 84)
migration_age_bin36 <- migration_age_bin36$mean
migration_age_bin36[43] <- migration_age_bin36[42]

migration_age_bin37 <- subset(migration_rates, age_group_id == 85)
migration_age_bin37 <- migration_age_bin37$mean
migration_age_bin37[43] <- migration_age_bin37[42]

migration_age_bin38 <- subset(migration_rates, age_group_id == 86)
migration_age_bin38 <- migration_age_bin38$mean
migration_age_bin38[43] <- migration_age_bin38[42]

migration_age_bin39 <- subset(migration_rates, age_group_id == 87)
migration_age_bin39 <- migration_age_bin39$mean
migration_age_bin39[43] <- migration_age_bin39[42]

migration_age_bin40 <- subset(migration_rates, age_group_id == 88)
migration_age_bin40 <- migration_age_bin40$mean
migration_age_bin40[43] <- migration_age_bin40[42]

migration_age_bin41 <- subset(migration_rates, age_group_id == 89)
migration_age_bin41 <- migration_age_bin41$mean
migration_age_bin41[43] <- migration_age_bin41[42]

migration_age_bin42 <- subset(migration_rates, age_group_id == 90)
migration_age_bin42 <- migration_age_bin42$mean
migration_age_bin42[43] <- migration_age_bin42[42]

migration_age_bin43 <- subset(migration_rates, age_group_id == 91)
migration_age_bin43 <- migration_age_bin43$mean
migration_age_bin43[43] <- migration_age_bin43[42]

migration_age_bin44 <- subset(migration_rates, age_group_id == 92)
migration_age_bin44 <- migration_age_bin44$mean
migration_age_bin44[43] <- migration_age_bin44[42]

migration_age_bin45 <- subset(migration_rates, age_group_id == 93)
migration_age_bin45 <- migration_age_bin45$mean
migration_age_bin45[43] <- migration_age_bin45[42]

migration_age_bin46 <- subset(migration_rates, age_group_id == 94)
migration_age_bin46 <- migration_age_bin46$mean
migration_age_bin46[43] <- migration_age_bin46[42]

migration_age_bin47 <- subset(migration_rates, age_group_id == 95)
migration_age_bin47 <- migration_age_bin47$mean
migration_age_bin47[43] <- migration_age_bin47[42]

migration_age_bin48 <- subset(migration_rates, age_group_id == 96)
migration_age_bin48 <- migration_age_bin48$mean
migration_age_bin48[43] <- migration_age_bin48[42]

migration_age_bin49 <- subset(migration_rates, age_group_id == 97)
migration_age_bin49 <- migration_age_bin49$mean
migration_age_bin49[43] <- migration_age_bin49[42]

migration_age_bin50 <- subset(migration_rates, age_group_id == 98)
migration_age_bin50 <- migration_age_bin50$mean
migration_age_bin50[43] <- migration_age_bin50[42]

migration_age_bin51 <- subset(migration_rates, age_group_id == 99)
migration_age_bin51 <- migration_age_bin51$mean
migration_age_bin51[43] <- migration_age_bin51[42]

migration_age_bin52 <- subset(migration_rates, age_group_id == 100)
migration_age_bin52 <- migration_age_bin52$mean
migration_age_bin52[43] <- migration_age_bin52[42]

migration_age_bin53 <- subset(migration_rates, age_group_id == 101)
migration_age_bin53 <- migration_age_bin53$mean
migration_age_bin53[43] <- migration_age_bin53[42]

migration_age_bin54 <- subset(migration_rates, age_group_id == 102)
migration_age_bin54 <- migration_age_bin54$mean
migration_age_bin54[43] <- migration_age_bin54[42]

migration_age_bin55 <- subset(migration_rates, age_group_id == 103)
migration_age_bin55 <- migration_age_bin55$mean
migration_age_bin55[43] <- migration_age_bin55[42]

migration_age_bin56 <- subset(migration_rates, age_group_id == 104)
migration_age_bin56 <- migration_age_bin56$mean
migration_age_bin56[43] <- migration_age_bin56[42]

migration_age_bin57 <- subset(migration_rates, age_group_id == 105)
migration_age_bin57 <- migration_age_bin57$mean
migration_age_bin57[43] <- migration_age_bin57[42]

migration_age_bin58 <- subset(migration_rates, age_group_id == 106)
migration_age_bin58 <- migration_age_bin58$mean
migration_age_bin58[43] <- migration_age_bin58[42]

migration_age_bin59 <- subset(migration_rates, age_group_id == 107)
migration_age_bin59 <- migration_age_bin59$mean
migration_age_bin59[43] <- migration_age_bin59[42]

migration_age_bin60 <- subset(migration_rates, age_group_id == 108)
migration_age_bin60 <- migration_age_bin60$mean
migration_age_bin60[43] <- migration_age_bin60[42]

migration_age_bin61 <- subset(migration_rates, age_group_id == 109)
migration_age_bin61 <- migration_age_bin61$mean
migration_age_bin61[43] <- migration_age_bin61[42]

migration_age_bin62 <- subset(migration_rates, age_group_id == 110)
migration_age_bin62 <- migration_age_bin62$mean
migration_age_bin62[43] <- migration_age_bin62[42]

migration_age_bin63 <- subset(migration_rates, age_group_id == 111)
migration_age_bin63 <- migration_age_bin63$mean
migration_age_bin63[43] <- migration_age_bin63[42]

migration_age_bin64 <- subset(migration_rates, age_group_id == 112)
migration_age_bin64 <- migration_age_bin64$mean
migration_age_bin64[43] <- migration_age_bin64[42]

migration_age_bin65 <- subset(migration_rates, age_group_id == 113)
migration_age_bin65 <- migration_age_bin65$mean
migration_age_bin65[43] <- migration_age_bin65[42]

migration_age_bin66 <- subset(migration_rates, age_group_id == 114)
migration_age_bin66 <- migration_age_bin66$mean
migration_age_bin66[43] <- migration_age_bin66[42]

migration_age_bin67 <- subset(migration_rates, age_group_id == 115)
migration_age_bin67 <- migration_age_bin67$mean
migration_age_bin67[43] <- migration_age_bin67[42]

migration_age_bin68 <- subset(migration_rates, age_group_id == 116)
migration_age_bin68 <- migration_age_bin68$mean
migration_age_bin68[43] <- migration_age_bin68[42]

migration_age_bin69 <- subset(migration_rates, age_group_id == 117)
migration_age_bin69 <- migration_age_bin69$mean
migration_age_bin69[43] <- migration_age_bin69[42]

migration_age_bin70 <- subset(migration_rates, age_group_id == 118)
migration_age_bin70 <- migration_age_bin70$mean
migration_age_bin70[43] <- migration_age_bin70[42]

migration_age_bin71 <- subset(migration_rates, age_group_id == 119)
migration_age_bin71 <- migration_age_bin71$mean
migration_age_bin71[43] <- migration_age_bin71[42]

migration_age_bin72 <- subset(migration_rates, age_group_id == 120)
migration_age_bin72 <- migration_age_bin72$mean
migration_age_bin72[43] <- migration_age_bin72[42]

migration_age_bin73 <- subset(migration_rates, age_group_id == 121)
migration_age_bin73 <- migration_age_bin73$mean
migration_age_bin73[43] <- migration_age_bin73[42]

migration_age_bin74 <- subset(migration_rates, age_group_id == 122)
migration_age_bin74 <- migration_age_bin74$mean
migration_age_bin74[43] <- migration_age_bin74[42]

migration_age_bin75 <- subset(migration_rates, age_group_id == 123)
migration_age_bin75 <- migration_age_bin75$mean
migration_age_bin75[43] <- migration_age_bin75[42]

migration_age_bin76 <- subset(migration_rates, age_group_id == 124)
migration_age_bin76 <- migration_age_bin76$mean
migration_age_bin76[43] <- migration_age_bin76[42]

migration_age_bin77 <- subset(migration_rates, age_group_id == 125)
migration_age_bin77 <- migration_age_bin77$mean
migration_age_bin77[43] <- migration_age_bin77[42]

migration_age_bin78 <- subset(migration_rates, age_group_id == 126)
migration_age_bin78 <- migration_age_bin78$mean
migration_age_bin78[43] <- migration_age_bin78[42]

migration_age_bin79 <- subset(migration_rates, age_group_id == 127)
migration_age_bin79 <- migration_age_bin79$mean
migration_age_bin79[43] <- migration_age_bin79[42]

migration_age_bin80 <- subset(migration_rates, age_group_id == 128)
migration_age_bin80 <- migration_age_bin80$mean
migration_age_bin80[43] <- migration_age_bin80[42]

migration_age_bin81 <- subset(migration_rates, age_group_id == 129)
migration_age_bin81 <- migration_age_bin81$mean
migration_age_bin81[43] <- migration_age_bin81[42]

migration_age_bin82 <- subset(migration_rates, age_group_id == 130)
migration_age_bin82 <- migration_age_bin82$mean
migration_age_bin82[43] <- migration_age_bin82[42]

migration_age_bin83 <- subset(migration_rates, age_group_id == 131)
migration_age_bin83 <- migration_age_bin83$mean
migration_age_bin83[43] <- migration_age_bin83[42]

migration_age_bin84 <- subset(migration_rates, age_group_id == 132)
migration_age_bin84 <- migration_age_bin84$mean
migration_age_bin84[43] <- migration_age_bin84[42]

migration_age_bin85 <- subset(migration_rates, age_group_id == 133)
migration_age_bin85 <- migration_age_bin85$mean
migration_age_bin85[43] <- migration_age_bin85[42]

migration_age_bin86 <- subset(migration_rates, age_group_id == 134)
migration_age_bin86 <- migration_age_bin86$mean
migration_age_bin86[43] <- migration_age_bin86[42]

migration_age_bin87 <- subset(migration_rates, age_group_id == 135)
migration_age_bin87 <- migration_age_bin87$mean
migration_age_bin87[43] <- migration_age_bin87[42]

migration_age_bin88 <- subset(migration_rates, age_group_id == 136)
migration_age_bin88 <- migration_age_bin88$mean
migration_age_bin88[43] <- migration_age_bin88[42]

migration_age_bin89 <- subset(migration_rates, age_group_id == 137)
migration_age_bin89 <- migration_age_bin89$mean
migration_age_bin89[43] <- migration_age_bin89[42]

migration_age_bin90 <- subset(migration_rates, age_group_id == 138)
migration_age_bin90 <- migration_age_bin90$mean
migration_age_bin90[43] <- migration_age_bin90[42]

migration_age_bin91 <- subset(migration_rates, age_group_id == 139)
migration_age_bin91 <- migration_age_bin91$mean
migration_age_bin91[43] <- migration_age_bin91[42]

migration_age_bin92 <- subset(migration_rates, age_group_id == 140)
migration_age_bin92 <- migration_age_bin92$mean
migration_age_bin92[43] <- migration_age_bin92[42]

migration_age_bin93 <- subset(migration_rates, age_group_id == 141)
migration_age_bin93 <- migration_age_bin93$mean
migration_age_bin93[43] <- migration_age_bin93[42]

migration_age_bin94 <- subset(migration_rates, age_group_id == 142)
migration_age_bin94 <- migration_age_bin94$mean
migration_age_bin94[43] <- migration_age_bin94[42]

migration_age_bin95 <- subset(migration_rates, age_group_id == 235 )
migration_age_bin95 <- migration_age_bin95$mean
migration_age_bin95[43] <- migration_age_bin95[42]


migration_matrix1 <- cbind(migration_age_bin_month_0, migration_age_bin_month_1, migration_age_bin_month_2, migration_age_bin_month_3,
                           migration_age_bin_month_4,  migration_age_bin_month_5, migration_age_bin_month_6, migration_age_bin_month_7, 
                           migration_age_bin_month_8,  migration_age_bin_month_9, migration_age_bin_month_10, migration_age_bin_month_11,
                           migration_age_bin1, migration_age_bin2, migration_age_bin3, migration_age_bin4,
                           migration_age_bin5, migration_age_bin6, migration_age_bin7, migration_age_bin8, migration_age_bin9, 
                           migration_age_bin10, migration_age_bin11, migration_age_bin12, migration_age_bin13, migration_age_bin14,
                           migration_age_bin15, migration_age_bin16, migration_age_bin17, migration_age_bin18, migration_age_bin19, 
                           migration_age_bin20, migration_age_bin21, migration_age_bin22, migration_age_bin23, migration_age_bin24,
                           migration_age_bin25, migration_age_bin26, migration_age_bin27, migration_age_bin28, migration_age_bin29, 
                           migration_age_bin30, migration_age_bin31, migration_age_bin32, migration_age_bin33, migration_age_bin34,
                           migration_age_bin35, migration_age_bin36, migration_age_bin37, migration_age_bin38, migration_age_bin39, 
                           migration_age_bin40, migration_age_bin41, migration_age_bin42, migration_age_bin43, migration_age_bin44,
                           migration_age_bin45, migration_age_bin46, migration_age_bin47, migration_age_bin48, migration_age_bin49, 
                           migration_age_bin50, migration_age_bin51, migration_age_bin52, migration_age_bin53, migration_age_bin54,
                           migration_age_bin55, migration_age_bin56, migration_age_bin57, migration_age_bin58, migration_age_bin59, 
                           migration_age_bin60, migration_age_bin61, migration_age_bin62, migration_age_bin63, migration_age_bin64,
                           migration_age_bin65, migration_age_bin66, migration_age_bin67, migration_age_bin68, migration_age_bin69, 
                           migration_age_bin70, migration_age_bin71, migration_age_bin72, migration_age_bin73, migration_age_bin74,
                           migration_age_bin75, migration_age_bin76, migration_age_bin77, migration_age_bin78, migration_age_bin79, 
                           migration_age_bin80, migration_age_bin81, migration_age_bin82, migration_age_bin83, migration_age_bin84,
                           migration_age_bin85, migration_age_bin86, migration_age_bin87, migration_age_bin88, migration_age_bin89, 
                           migration_age_bin90, migration_age_bin91, migration_age_bin92, migration_age_bin93, migration_age_bin94,
                           migration_age_bin95)

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
full_mat <- data.table()
for(i in 1:43){
  temp <- rep.row(migration_matrix1[i,],53)
  full_mat <- rbind(full_mat, temp)
}

migration_matrix <- copy(full_mat)


save(migration_matrix, 
     file = paste0('FILEPATH',iso3,'/', processing_date,'/', iso3,'_processed_migration.RData'))








migration_matrix_weekly <- matrix(migration_age_bin_month_0, 43, 4)
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin_month_1, 43, 4))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin_month_2, 43, 4))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin_month_3, 43, 4))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin_month_4, 43, 4))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin_month_5, 43, 5))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin_month_6, 43, 4))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin_month_7, 43, 5))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin_month_8, 43, 5))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin_month_9, 43, 4))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin_month_10, 43, 4))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin_month_11, 43, 5))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin1, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin2, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin3, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin4, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin5, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin6, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin7, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin8, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin9, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin10, 43, 53))

migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin11, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin12, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin13, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin14, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin15, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin16, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin17, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin18, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin19, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin20, 43, 53))

migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin21, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin22, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin23, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin24, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin25, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin26, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin27, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin28, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin29, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin30, 43, 53))

migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin31, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin32, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin33, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin34, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin35, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin36, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin37, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin38, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin39, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin40, 43, 53))

migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin41, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin42, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin43, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin44, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin45, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin46, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin47, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin48, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin49, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin50, 43, 53))

migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin51, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin52, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin53, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin54, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin55, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin56, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin57, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin58, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin59, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin60, 43, 53))

migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin61, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin62, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin63, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin64, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin65, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin66, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin67, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin68, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin69, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin70, 43, 53))

migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin71, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin72, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin73, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin74, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin75, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin76, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin77, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin78, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin79, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin80, 43, 53))

migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin81, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin82, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin83, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin84, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin85, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin86, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin87, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin88, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin89, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin90, 43, 53))

migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin91, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin92, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin93, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin94, 43, 53))
migration_matrix_weekly <- cbind(migration_matrix_weekly, matrix(migration_age_bin95, 43, 1))




dim(migration_matrix_weekly)





rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
full_mat <- data.table()
for(i in 1:43){
  temp <- rep.row(migration_matrix_weekly[i,],53)
  full_mat <- rbind(full_mat, temp)
}

migration_matrix_weekly <- copy(full_mat)





save(migration_matrix_weekly ,
     file = paste0('FILEPATH', iso3,'/',processing_date,'/', iso3,'_processed_migration.RData'))




