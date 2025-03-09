rm(list=ls())

pred_list_0_4 <- array(data = NA, c(100,371,79))
pred_list_5_9 <- array(data = NA, c(100,371,79))
pred_list_10_14 <- array(data = NA, c(100,371,79))
beta_max_list <- array(data = NA, c(100,11))
beta_min_list <- array(data = NA, c(100,11))
rho_list <- array(data = NA, c(100,11))
final_ll <- array(data = NA, c(100,11))

run_date <- 'RUNDATE'

for(i in 1:100){
  message(i)
  load(paste0('FILEPATH/',run_date,'/bootstrap_',i,'_BCD_iteration_', 10, '.RData') )
  #pred_list[[i]] <- I_predicted
  
  
  I_annual_by_district_by_week_0_4 <- matrix(0, 2120, 79)
  I_annual_by_district_by_week_5_9 <- matrix(0, 2120, 79)
  I_annual_by_district_by_week_10_14 <- matrix(0, 2120, 79)
  for (w in 1750:2120){
    for(d in 1:79){
      for(a in 1:16){
        I_annual_by_district_by_week_0_4[w,d] = I_annual_by_district_by_week_0_4[w, d] + (I_predicted[w, a, d] * invlogit(logit_rho_vals[11]))
      }
      I_annual_by_district_by_week_5_9[w,d] = I_annual_by_district_by_week_5_9[w, d] + (I_predicted[w, 17, d]* invlogit(logit_rho_vals[11]))
      I_annual_by_district_by_week_10_14[w,d] = I_annual_by_district_by_week_10_14[w, d] + (I_predicted[w, 18, d]* invlogit(logit_rho_vals[11]))
    }
  }
  pred_list_0_4[i,,] <- I_annual_by_district_by_week_0_4[1750:2120,1:79]
  pred_list_5_9[i,,] <- I_annual_by_district_by_week_5_9[1750:2120,1:79]
  pred_list_10_14[i,,] <- I_annual_by_district_by_week_10_14[1750:2120,1:79]
  
  beta_max_list[i,] <- (beta_max_vals)
  beta_min_list[i,] <- (beta_min_vals)
  rho_list[i,] <- (logit_rho_vals)
  final_ll[i,] <- rho_val_ll
}



rho_long <- reshape2::melt(rho_list)
colnames(rho_long) <- c("bootstrap","iteration","value")


beta_max_long <- reshape2::melt(beta_max_list)
colnames(beta_max_long) <- c("bootstrap","iteration","value")


beta_min_long <- reshape2::melt(beta_min_list)
colnames(beta_min_long) <- c("bootstrap","iteration","value")


rho_val_ll_long <- reshape2::melt(final_ll)
colnames(rho_val_ll_long) <- c("bootstrap","iteration","value")



load(paste0('FILEPATH/', 'ETH','_processed_cases_loess.RData')) # subnational


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
cases_age_bin0_4_matrix[cases_age_bin0_4_matrix < 0] <- 0
cases_age_bin5_9_matrix[cases_age_bin5_9_matrix < 0] <- 0
cases_age_bin10_14_matrix[cases_age_bin10_14_matrix < 0] <- 0
cases_age_bin15_24_matrix[cases_age_bin15_24_matrix < 0] <- 0
cases_age_bin25_34_matrix[cases_age_bin25_34_matrix < 0] <- 0
cases_age_bin35_44_matrix[cases_age_bin35_44_matrix < 0] <- 0
cases_age_bin45_54_matrix[cases_age_bin45_54_matrix < 0] <- 0
cases_age_bin55_64_matrix[cases_age_bin55_64_matrix < 0] <- 0
cases_age_bin65_plus_matrix[cases_age_bin65_plus_matrix < 0] <- 0


ln_rho <-ggplot(rho_long, aes(x=iteration, y=value, group = bootstrap)) +
  geom_line(alpha = 0.2) + theme_bw() #+ facet_wrap(~bootstrap)
ln_rho


ln_beta_min<-ggplot(beta_min_long, aes(x=iteration, y=value, group = bootstrap)) +
  geom_line(alpha = 0.2) + theme_bw() #+ facet_wrap(~bootstrap)
ln_beta_min


ln_beta_max <-ggplot(beta_max_long, aes(x=iteration, y=value, group = bootstrap)) +
  geom_line(alpha = 0.2) + theme_bw() #+ facet_wrap(~bootstrap)
ln_beta_max

ln_lls <-ggplot(rho_val_ll_long, aes(x=iteration, y=value, group = bootstrap)) +
  geom_line(alpha = 0.2) + theme_bw() #+ facet_wrap(~bootstrap)
ln_lls



pred_0_4_long <- reshape2::melt(pred_list_0_4)
colnames(pred_0_4_long) <- c("bootstrap","week","district","value")


pred_5_9_long <- reshape2::melt(pred_list_5_9)
colnames(pred_5_9_long) <- c("bootstrap","week","district","value")


pred_10_14_long <- reshape2::melt(pred_list_10_14)
colnames(pred_10_14_long) <- c("bootstrap","week","district","value")

colnames(cases_age_bin0_4_matrix) <- 1:79
cases_age_bin0_4_matrix_long <- reshape2::melt(cases_age_bin0_4_matrix)
colnames(cases_age_bin0_4_matrix_long) <- c('week', 'district', 'value')
cases_age_bin0_4_matrix_long <- subset(cases_age_bin0_4_matrix_long, week > 1749)

colnames(cases_age_bin5_9_matrix) <- 1:79
cases_age_bin5_9_matrix_long <- reshape2::melt(cases_age_bin5_9_matrix)
colnames(cases_age_bin5_9_matrix_long) <- c('week', 'district', 'value')
cases_age_bin5_9_matrix_long <- subset(cases_age_bin5_9_matrix_long, week > 1749)


colnames(cases_age_bin10_14_matrix) <- 1:79
cases_age_bin10_14_matrix_long <- reshape2::melt(cases_age_bin10_14_matrix)
colnames(cases_age_bin10_14_matrix_long) <- c('week', 'district', 'value')
cases_age_bin10_14_matrix_long <- subset(cases_age_bin10_14_matrix_long, week > 1749)


ln_0_4 <-ggplot(pred_0_4_long, aes(x=week, y=value)) +
  geom_line(alpha = 0.2, color='dodgerblue4', aes(group = bootstrap)) + theme_bw() + 
  facet_wrap(~district, scales = 'free') + xlab("Week") + ylab("Cases") + 
  geom_line(data =cases_age_bin0_4_matrix_long, aes(x=week-1749, y=value), color='black', lwd=1) 
ln_0_4


ln_5_9 <-ggplot(pred_5_9_long, aes(x=week, y=value)) +
  geom_line(alpha = 0.2, color='darkgreen', aes(group = bootstrap)) + theme_bw() + 
  facet_wrap(~district, scales = 'free') +xlab("Week") + ylab("Cases") + 
  geom_line(data =cases_age_bin5_9_matrix_long, aes(x=week-1749, y=value), color='black', lwd=1) 
ln_5_9


ln_10_14 <-ggplot(pred_10_14_long, aes(x=week, y=value)) +
  geom_line(alpha = 0.2, color='orange', aes(group = bootstrap)) + theme_bw() + 
  facet_wrap(~district, scales = 'free') +xlab("Week") + ylab("Cases") + 
  geom_line(data =cases_age_bin10_14_matrix_long, aes(x=week-1749, y=value), color='black', lwd=1) 
ln_10_14


rho_long_final <- subset(rho_long, iteration ==11)
quantile(rho_long_final$value, probs = c(0.5, 0.025, 0.975), na.rm=T)


beta_max_long_final <- subset(beta_max_long, iteration ==11)
quantile(beta_max_long_final$value, probs = c(0.5, 0.025, 0.975), na.rm=T)


beta_min_long_final <- subset(beta_min_long, iteration ==11)
quantile(beta_min_long_final$value, probs = c(0.5, 0.025, 0.975), na.rm=T)


ll_long_final <- subset(rho_val_ll_long, iteration ==11)
quantile(ll_long_final$value, probs = c(0.5, 0.025, 0.975), na.rm=T)


beta_min_long_final$variable <- 'beta_min'
beta_max_long_final$variable <- 'beta_max'
rho_long_final$variable <- 'logit_rho'

all_pars <- rbind(beta_max_long_final, beta_min_long_final, rho_long_final)

gg_density <- ggplot(data = all_pars) + geom_density(aes(x = value, fill=variable)) + facet_wrap(~variable, scales = 'free') + theme_bw()



dir.create(paste0("FILEPATH/",run_date,"/diagnostics/"), recursive = T)



png(file = paste0("FILEPATH/",run_date,"/diagnostics/beta_max_iterations.png"),
    width = 14,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
plot(ln_beta_max)
dev.off()


png(file = paste0("FILEPATH/",run_date,"/diagnostics/beta_min_iterations.png"),
    width = 14,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
plot(ln_beta_min)
dev.off()


png(file = paste0("FILEPATH/",run_date,"/diagnostics/rho_iterations.png"),
    width = 14,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
plot(ln_rho)
dev.off()



png(file = paste0("FILEPATH/",run_date,"/diagnostics/ll_iterations.png"),
    width = 14,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
plot(ln_lls)
dev.off()




png(file = paste0("FILEPATH/",run_date,"/diagnostics/cases_0_4.png"),
    width = 14,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
plot(ln_0_4)
dev.off()


png(file = paste0("FILEPATH/",run_date,"/diagnostics/cases_5_9.png"),
    width = 14,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
plot(ln_5_9)
dev.off()





png(file = paste0("FILEPATH/",run_date,"/diagnostics/cases_10_14.png"),
    width = 14,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
plot(ln_10_14)
dev.off()



png(file = paste0("FILEPATH/",run_date,"/diagnostics/density_plot.png"),
    width = 14,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
plot(gg_density)
dev.off()


ln_0_4_same_scale <-ggplot(pred_0_4_long, aes(x=week, y=value)) +
  geom_line(alpha = 0.2, color='dodgerblue4', aes(group = bootstrap)) + theme_bw() + 
  facet_wrap(~district) + xlab("Week") + ylab("Cases") + 
  geom_line(data =cases_age_bin0_4_matrix_long, aes(x=week-1749, y=value), color='black', lwd=1) 
ln_0_4_same_scale


ln_5_9_same_scale <-ggplot(pred_5_9_long, aes(x=week, y=value)) +
  geom_line(alpha = 0.2, color='darkgreen', aes(group = bootstrap)) + theme_bw() + 
  facet_wrap(~district) +xlab("Week") + ylab("Cases") + 
  geom_line(data =cases_age_bin5_9_matrix_long, aes(x=week-1749, y=value), color='black', lwd=1) 
ln_5_9_same_scale


ln_10_14_same_scale <-ggplot(pred_10_14_long, aes(x=week, y=value)) +
  geom_line(alpha = 0.2, color='orange', aes(group = bootstrap)) + theme_bw() + 
  facet_wrap(~district) +xlab("Week") + ylab("Cases") + 
  geom_line(data =cases_age_bin0_4_matrix_long, aes(x=week-1749, y=value), color='black', lwd=1) 
ln_10_14_same_scale





png(file = paste0("FILEPATH/",run_date,"/diagnostics/cases_0_4_same_scale.png"),
    width = 14,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
plot(ln_0_4_same_scale)
dev.off()


png(file = paste0("FILEPATH/",run_date,"/diagnostics/cases_5_9_same_scale.png"),
    width = 14,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
plot(ln_5_9_same_scale)
dev.off()





png(file = paste0("FILEPATH/",run_date,"/diagnostics/cases_10_14_same_scale.png"),
    width = 14,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
plot(ln_10_14_same_scale)
dev.off()


load(paste0('FILEPATH/ETH_processed_population_by_epiweek_measles_model_age_bins_by_adm2.RData'))


pop_array_red <- array(NA, dim=c(2120, 24,79))
pop_array_red[,1:16,] = pop_array[1:2120,1:16,1:79]

for(i in 1750:2120){
  for(d in 1:79){
    
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


pop_array_red_0_4 <- array(NA, dim=c(2120, 79))
pop_array_red_5_9 <- (pop_array_red[,17,])
pop_array_red_10_14 <- pop_array_red[,18,]
for(i in 1750:2120){
  for(d in 1:79){
    pop_array_red_0_4[i,d] = sum(pop_array[i,1:16, d])
  }
}



pop_array_red_0_4_long <- reshape2::melt(pop_array_red_0_4)
pop_array_red_5_9_long <- reshape2::melt(pop_array_red_5_9)
pop_array_red_10_14_long <- reshape2::melt(pop_array_red_10_14)

colnames(pop_array_red_0_4_long) <- c("week", "district", "pop")
colnames(pop_array_red_5_9_long) <- c("week", "district", "pop")
colnames(pop_array_red_10_14_long) <- c("week", "district", "pop")


pop_array_red_0_4_long <- subset(pop_array_red_0_4_long, week > 1749)
pop_array_red_5_9_long <- subset(pop_array_red_5_9_long, week > 1749)
pop_array_red_10_14_long <- subset(pop_array_red_10_14_long, week > 1749)


pop_array_red_0_4_long$week <- pop_array_red_0_4_long$week - 1749
pop_array_red_5_9_long$week <- pop_array_red_5_9_long$week - 1749
pop_array_red_10_14_long$week <- pop_array_red_10_14_long$week - 1749


pred_pop_0_4 <- merge(pred_0_4_long, pop_array_red_0_4_long, by=c('week', 'district'))
pred_pop_5_9 <- merge(pred_5_9_long, pop_array_red_5_9_long, by=c('week', 'district'))
pred_pop_10_14 <- merge(pred_10_14_long, pop_array_red_10_14_long, by=c('week', 'district'))


cases_age_bin0_4_matrix_long$week <- cases_age_bin0_4_matrix_long$week-1749
cases_age_bin5_9_matrix_long$week <- cases_age_bin5_9_matrix_long$week-1749
cases_age_bin10_14_matrix_long$week <- cases_age_bin10_14_matrix_long$week-1749

cases_pop_0_4 <- merge(cases_age_bin0_4_matrix_long, pop_array_red_0_4_long, by=c('week', 'district'))
cases_pop_5_9 <- merge(cases_age_bin5_9_matrix_long, pop_array_red_5_9_long, by=c('week', 'district'))
cases_pop_10_14 <- merge(cases_age_bin10_14_matrix_long, pop_array_red_10_14_long, by=c('week', 'district'))




ln_0_4_same_scale_incidence <-ggplot(pred_pop_0_4, aes(x=week, y=value/pop)) +
  geom_line(alpha = 0.2, color='lightblue', aes(group = bootstrap)) + theme_bw() + 
  facet_wrap(~district) + xlab("Week") + ylab("Incidence") + 
  geom_line(data =cases_pop_0_4, aes(x=week, y=value/pop), color='black', lwd=1) 
ln_0_4_same_scale_incidence



ln_0_4_incidence <-ggplot(pred_pop_0_4, aes(x=week, y=value/pop)) +
  geom_line(alpha = 0.2, color='lightblue', aes(group = bootstrap)) + theme_bw() + 
  facet_wrap(~district, scales = 'free') + xlab("Week") + ylab("Incidence") + 
  geom_line(data =cases_pop_0_4, aes(x=week, y=value/pop), color='black', lwd=1) 
ln_0_4_incidence



ln_5_9_same_scale_incidence <-ggplot(pred_pop_5_9, aes(x=week, y=value/pop)) +
  geom_line(alpha = 0.2, color='darkgreen', aes(group = bootstrap)) + theme_bw() + 
  facet_wrap(~district) + xlab("Week") + ylab("Incidence") + 
  geom_line(data =cases_pop_5_9, aes(x=week, y=value/pop), color='black', lwd=1) 
ln_5_9_same_scale_incidence



ln_5_9_incidence <-ggplot(pred_pop_5_9, aes(x=week, y=value/pop)) +
  geom_line(alpha = 0.2, color='darkgreen', aes(group = bootstrap)) + theme_bw() + 
  facet_wrap(~district, scales = 'free') + xlab("Week") + ylab("Incidence") + 
  geom_line(data =cases_pop_5_9, aes(x=week, y=value/pop), color='black', lwd=1) 
ln_5_9_incidence



ln_10_14_same_scale_incidence <-ggplot(pred_pop_10_14, aes(x=week, y=value/pop)) +
  geom_line(alpha = 0.2, color='orange', aes(group = bootstrap)) + theme_bw() + 
  facet_wrap(~district) + xlab("Week") + ylab("Incidence") + 
  geom_line(data =cases_pop_10_14, aes(x=week, y=value/pop), color='black', lwd=1) 
ln_10_14_same_scale_incidence



ln_10_14_incidence <-ggplot(pred_pop_10_14, aes(x=week, y=value/pop)) +
  geom_line(alpha = 0.2, color='orange', aes(group = bootstrap)) + theme_bw() + 
  facet_wrap(~district, scales = 'free') + xlab("Week") + ylab("Incidence") + 
  geom_line(data =cases_pop_10_14, aes(x=week, y=value/pop), color='black', lwd=1) 
ln_10_14_incidence




png(file = paste0("FILEPATH/",run_date,"/diagnostics/cases_0_4_same_scale_incidence.png"),
    width = 14,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
plot(ln_0_4_same_scale_incidence)
dev.off()


png(file = paste0("FILEPATH/",run_date,"/diagnostics/cases_5_9_same_scale_incidence.png"),
    width = 14,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
plot(ln_5_9_same_scale_incidence)
dev.off()





png(file = paste0("FILEPATH/",run_date,"/diagnostics/cases_10_14_same_scale_incidence.png"),
    width = 14,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
plot(ln_10_14_same_scale_incidence)
dev.off()


png(file = paste0("FILEPATH/",run_date,"/diagnostics/cases_0_4_incidence.png"),
    width = 14,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
plot(ln_0_4_incidence)
dev.off()


png(file = paste0("FILEPATH/",run_date,"/diagnostics/cases_5_9_incidence.png"),
    width = 14,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
plot(ln_5_9_incidence)
dev.off()





png(file = paste0("FILEPATH/",run_date,"/diagnostics/cases_10_14_incidence.png"),
    width = 14,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
plot(ln_10_14_incidence)
dev.off()

