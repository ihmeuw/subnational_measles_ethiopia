//#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

#include <omp.h>
// [[Rcpp::plugins(openmp)]]


// [[Rcpp::export]]
double dpois_nonint(double x, double lambda){
  
  double logresult = 0;
  
  logresult = -lambda + x*log(lambda) - lgamma(x+1);
  return(logresult);
}

// [[Rcpp::export]]
double dnbinom_nonint(double x, double mu, double size){
  
  double logresult = 0;
  double p =0;
  
  p = size / (size + mu);
  logresult = lgamma(x+size) - lgamma(size) - lgamma(x+1) + size*log(p)+x*log(1-p);
  return(logresult);
}
//#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

#include <omp.h>
// [[Rcpp::plugins(openmp)]]


// [[Rcpp::export]]
double Priorcpp_measles(double max_beta, 
                        double min_beta, 
                        
                        double logit_rho1, 
                        double logit_vax_efficacy) {
  //return 1/sqrt(2*M_PI)*exp(-pow(beta0,2)/2);
  return R::dnorm(max_beta, 3, 0.75, true) +  
    R::dnorm(min_beta, 1, 0.75, true) +  
    R::dnorm(logit_rho1, -5.9, 0.5, true) + 
    R::dnorm(logit_vax_efficacy, 1.1, 0.05, true);
  
}

//   
// // [[Rcpp::export]]
// NumericMatrix mobility_thetas(NumericMatrix G, double theta){
//   
//   double one_minus_theta = 0;
//   NumericVector G_denoms(79);
//   NumericMatrix M(79,79);
//   
//   one_minus_theta = 1 - theta;
//   G_denoms = Cpp_colSums(G);
//   
//   for (int i = 0; i < 79; i++) {
//     M(_,i) = (G(_,i) / G_denoms(i)) * one_minus_theta;
//     M(i,i) = theta;
//   }
//   
//   return(M);
// }

// [[Rcpp::export]]
double Likelihoodcpp_measles_transmission2_omp(NumericVector i_ticker_1, 
                                           NumericVector W_vector, NumericMatrix M, 
                                           NumericMatrix mortality_matrix, NumericMatrix migration_matrix, 
                                           NumericVector national_cases, 
                                           NumericVector pop,
                                           double all_strat_cases,
                                           NumericVector cases_age_bin0_4, NumericVector cases_age_bin5_9, 
                                           NumericVector cases_age_bin10_14, NumericVector cases_age_bin15_19, 
                                           NumericVector cases_age_bin20_24, NumericVector cases_age_bin25_29, 
                                           NumericVector cases_age_bin30_34, NumericVector cases_age_bin35_39, 
                                           NumericVector cases_age_bin40_44, NumericVector cases_age_bin45_49, 
                                           NumericVector cases_age_bin50_54, NumericVector cases_age_bin55_59, 
                                           NumericVector cases_age_bin60_64, NumericVector cases_age_bin65_69, 
                                           NumericVector cases_age_bin70_74, NumericVector cases_age_bin75_79, 
                                           NumericVector cases_age_bin80_84, 
                                           
                                           NumericMatrix cases_age_bin0_4_matrix, 
                                           NumericMatrix cases_age_bin5_9_matrix, 
                                           NumericMatrix cases_age_bin10_14_matrix, 
                                           NumericMatrix cases_age_bin15_24_matrix, 
                                           NumericMatrix cases_age_bin25_34_matrix, 
                                           NumericMatrix cases_age_bin35_44_matrix, 
                                           NumericMatrix cases_age_bin45_54_matrix, 
                                           NumericMatrix cases_age_bin55_64_matrix, 
                                           NumericMatrix cases_age_bin65_plus_matrix, 
                                           
                                           NumericMatrix zero_matrix,
                                           
                                           NumericMatrix births,
                                           
                                           NumericVector age_cov_prevalence_counts_vector, NumericVector age_cov2_prevalence_counts_vector,
                                           
                                           NumericVector pop_array_vector,
                                           
                                           NumericVector pop_array_vector_district,
                                           
                                           NumericMatrix autocorrelation,
                                           
                                           NumericVector max_beta,
                                           NumericVector min_beta,
                                           
                                           double rho1, 
                                           double rho2, 
                                           double rho3, 
                                           double rho4, 
                                           double rho5, 
                                           double rho6, 
                                           double rho7, 
                                           double rho8, 
                                           double rho9, 
                                           double rho10, 
                                           double rho11,
                                           
                                           NumericVector logit_vax_efficacy,
                                           NumericMatrix bootstrap_lookup) {
  double likelihood = 0;
  double loglike= 0;
  double foi_beta0= 0;
  double foi_beta1= 0;
  double foi_cos_inside1a= 0;
  double foi_cos_inside1b= 0; 
  double foi_cos_inside2= 0;
  double foi_cos= 0;
  double composite_beta_part1= 0;
  double composite_beta_part2= 0;
  double composite_beta= 0;
  double amplitude= 0; 
  double vertical_shift= 0; 
  NumericVector mort_mig_val(24);
  NumericVector to_export(40);
  
  double rho;
  
  double vax_efficacy;
  double composite_dose1_efficacy;
  double composite_dose1_efficacy_1minus;
  double composite_dose2_efficacy;
  double composite_dose2_efficacy_1minus;
  
  int nobs_stratified = 2120;
  int nadmin_units = 79;
  int nage_groups = 24;
  
  NumericMatrix foi(nage_groups, nadmin_units);
  NumericMatrix foi2(nage_groups, nadmin_units);
  
  // NumericMatrix S_unvax_premort_preage_prematernal(nage_groups,nadmin_units);
  // NumericMatrix S_vax_premort_preage_prematernal(nage_groups,nadmin_units);
  // NumericMatrix S_vax2_premort_preage_prematernal(nage_groups,nadmin_units);
  
  // NumericMatrix I_unvax_premort_preage(nage_groups,nadmin_units);
  // NumericMatrix I_vax_premort_preage(nage_groups,nadmin_units);
  // NumericMatrix I_vax2_premort_preage(nage_groups,nadmin_units);
  
  NumericVector numerator_temp(18);
  
  // NumericMatrix R_unvax_premort_preage(nage_groups,nadmin_units);
  // NumericMatrix R_vax_premort_preage(nage_groups,nadmin_units);
  // NumericMatrix R_vax2_premort_preage(nage_groups,nadmin_units);
  
  // NumericMatrix R_unvax_preage(nage_groups,nadmin_units);
  // NumericMatrix R_vax_preage(nage_groups,nadmin_units);
  // NumericMatrix R_vax2_preage(nage_groups,nadmin_units);
  
  NumericMatrix M_unvax_preage_prematernal(nage_groups,nadmin_units);
  NumericMatrix M_vax_preage_prematernal(nage_groups,nadmin_units);
  NumericMatrix M_vax2_preage_prematernal(nage_groups,nadmin_units);
  
  NumericMatrix S_unvax_preage_prematernal(nage_groups,nadmin_units);
  NumericMatrix S_vax_preage_prematernal(nage_groups,nadmin_units);
  NumericMatrix S_vax2_preage_prematernal(nage_groups,nadmin_units);
  
  NumericMatrix I_unvax_preage(nage_groups,nadmin_units);
  NumericMatrix I_vax_preage(nage_groups,nadmin_units);
  NumericMatrix I_vax2_preage(nage_groups,nadmin_units);
  
  NumericMatrix I_unvax_preage_INCIDENCE(nage_groups,nadmin_units);
  NumericMatrix I_vax_preage_INCIDENCE(nage_groups,nadmin_units);
  NumericMatrix I_vax2_preage_INCIDENCE(nage_groups,nadmin_units);
  
  NumericMatrix maternal_unvax_calc(nage_groups,nadmin_units);
  NumericMatrix maternal_vax_calc(nage_groups,nadmin_units);
  NumericMatrix maternal_vax2_calc(nage_groups,nadmin_units);
  
  NumericMatrix M_unvax_preage(nage_groups,nadmin_units);
  NumericMatrix M_vax_preage(nage_groups,nadmin_units);
  NumericMatrix M_vax2_preage(nage_groups,nadmin_units);
  
  NumericMatrix S_unvax_preage(nage_groups,nadmin_units);
  NumericMatrix S_vax_preage(nage_groups,nadmin_units);
  NumericMatrix S_vax2_preage(nage_groups,nadmin_units);
  
  NumericMatrix R_unvax_precalib(nage_groups,nadmin_units);
  NumericMatrix R_vax_precalib(nage_groups,nadmin_units);
  NumericMatrix R_vax2_precalib(nage_groups,nadmin_units);
  
  NumericMatrix M_unvax_precalib(nage_groups,nadmin_units);
  NumericMatrix M_vax_precalib(nage_groups,nadmin_units);
  NumericMatrix M_vax2_precalib(nage_groups,nadmin_units);
  
  NumericMatrix S_unvax_precalib(nage_groups,nadmin_units);
  NumericMatrix S_vax_precalib(nage_groups,nadmin_units);
  NumericMatrix S_vax2_precalib(nage_groups,nadmin_units);
  
  NumericMatrix I_unvax_precalib(nage_groups,nadmin_units);
  NumericMatrix I_vax_precalib(nage_groups,nadmin_units);
  NumericMatrix I_vax2_precalib(nage_groups,nadmin_units);
  
  NumericMatrix I_unvax_precalib_INCIDENCE(nage_groups,nadmin_units);
  NumericMatrix I_vax_precalib_INCIDENCE(nage_groups,nadmin_units);
  NumericMatrix I_vax2_precalib_INCIDENCE(nage_groups,nadmin_units);
  
  NumericMatrix prop_denominator(nage_groups,nadmin_units);
  
  NumericMatrix M_unvax_prevax(nage_groups,nadmin_units);
  NumericMatrix M_vax_prevax(nage_groups,nadmin_units);
  NumericMatrix M_vax2_prevax(nage_groups,nadmin_units);
  
  NumericMatrix S_unvax_prevax(nage_groups,nadmin_units);
  NumericMatrix S_vax_prevax(nage_groups,nadmin_units);
  NumericMatrix S_vax2_prevax(nage_groups,nadmin_units);
  
  NumericMatrix I_unvax_prevax(nage_groups,nadmin_units);
  NumericMatrix I_vax_prevax(nage_groups,nadmin_units);
  NumericMatrix I_vax2_prevax(nage_groups,nadmin_units);
  
  NumericMatrix I_unvax_prevax_INCIDENCE(nage_groups,nadmin_units);
  NumericMatrix I_vax_prevax_INCIDENCE(nage_groups,nadmin_units);
  NumericMatrix I_vax2_prevax_INCIDENCE(nage_groups,nadmin_units);
  
  NumericMatrix R_unvax_prevax(nage_groups,nadmin_units);
  NumericMatrix R_vax_prevax(nage_groups,nadmin_units);
  NumericMatrix R_vax2_prevax(nage_groups,nadmin_units);
  
  NumericMatrix M_unvax_prevax_counts(nage_groups,nadmin_units);
  NumericMatrix M_vax_prevax_counts(nage_groups,nadmin_units);
  NumericMatrix M_vax2_prevax_counts(nage_groups,nadmin_units);
  
  NumericMatrix S_unvax_prevax_counts(nage_groups,nadmin_units);
  NumericMatrix S_vax_prevax_counts(nage_groups,nadmin_units);
  NumericMatrix S_vax2_prevax_counts(nage_groups,nadmin_units);
  
  NumericMatrix I_unvax_prevax_counts(nage_groups,nadmin_units);
  NumericMatrix I_vax_prevax_counts(nage_groups,nadmin_units);
  NumericMatrix I_vax2_prevax_counts(nage_groups,nadmin_units);
  
  NumericMatrix I_unvax_prevax_counts_INCIDENCE(nage_groups,nadmin_units);
  NumericMatrix I_vax_prevax_counts_INCIDENCE(nage_groups,nadmin_units);
  NumericMatrix I_vax2_prevax_counts_INCIDENCE(nage_groups,nadmin_units);
  
  NumericMatrix R_unvax_prevax_counts(nage_groups,nadmin_units);
  NumericMatrix R_vax_prevax_counts(nage_groups,nadmin_units);
  NumericMatrix R_vax2_prevax_counts(nage_groups,nadmin_units);
  
  NumericMatrix unvax_denom(nage_groups,nadmin_units);
  NumericMatrix temp_vax_calc(nage_groups,nadmin_units);
  NumericMatrix vax2_calc(nage_groups,nadmin_units);
  NumericMatrix vax_denom(nage_groups,nadmin_units);
  
  NumericMatrix M_prop(nage_groups,nadmin_units);
  NumericMatrix S_prop(nage_groups,nadmin_units);
  NumericMatrix I_prop(nage_groups,nadmin_units);
  NumericMatrix R_prop(nage_groups,nadmin_units);
  
  NumericMatrix M_unvax_postvax_counts(nage_groups,nadmin_units);
  NumericMatrix S_unvax_postvax_counts(nage_groups,nadmin_units);
  NumericMatrix I_unvax_postvax_counts(nage_groups,nadmin_units);
  NumericMatrix I_unvax_postvax_counts_INCIDENCE(nage_groups,nadmin_units);
  NumericMatrix R_unvax_postvax_counts(nage_groups,nadmin_units);
  
  NumericMatrix M_vax_postvax_counts(nage_groups,nadmin_units);
  NumericMatrix S_vax_postvax_counts(nage_groups,nadmin_units);
  NumericMatrix I_vax_postvax_counts(nage_groups,nadmin_units);
  NumericMatrix I_vax_postvax_counts_INCIDENCE(nage_groups,nadmin_units);
  NumericMatrix R_vax_postvax_counts(nage_groups,nadmin_units);
  
  NumericMatrix M_prop2(nage_groups,nadmin_units);
  NumericMatrix S_prop2(nage_groups,nadmin_units);
  NumericMatrix I_prop2(nage_groups,nadmin_units);
  NumericMatrix R_prop2(nage_groups,nadmin_units);
  
  NumericMatrix M_vax_postvax2_counts(nage_groups,nadmin_units);
  NumericMatrix S_vax_postvax2_counts(nage_groups,nadmin_units);
  NumericMatrix I_vax_postvax2_counts(nage_groups,nadmin_units);
  NumericMatrix I_vax_postvax2_counts_INCIDENCE(nage_groups,nadmin_units);
  NumericMatrix R_vax_postvax2_counts(nage_groups,nadmin_units);
  
  NumericMatrix M_vax2_postvax2_counts(nage_groups,nadmin_units);
  NumericMatrix S_vax2_postvax2_counts(nage_groups,nadmin_units);
  NumericMatrix I_vax2_postvax2_counts(nage_groups,nadmin_units);
  NumericMatrix I_vax2_postvax2_counts_INCIDENCE(nage_groups,nadmin_units);
  NumericMatrix R_vax2_postvax2_counts(nage_groups,nadmin_units);
  
  NumericMatrix temp_calc_efficacy_1minus(nage_groups,nadmin_units);
  NumericMatrix temp_calc_efficacy(nage_groups,nadmin_units);
  NumericMatrix temp_calc2_efficacy_1minus(nage_groups,nadmin_units);
  NumericMatrix temp_calc2_efficacy(nage_groups,nadmin_units);
  
  NumericMatrix M_unvax_final(nage_groups, nadmin_units);
  NumericMatrix M_vax_final(nage_groups, nadmin_units);
  NumericMatrix M_vax2_final(nage_groups, nadmin_units);
  NumericMatrix M_final(nage_groups, nadmin_units);
  
  NumericMatrix S_unvax_final(nage_groups, nadmin_units);
  NumericMatrix S_vax_final(nage_groups, nadmin_units);
  NumericMatrix S_vax2_final(nage_groups, nadmin_units);
  NumericMatrix S_final(nage_groups, nadmin_units);
  
  NumericMatrix I_unvax_final(nage_groups, nadmin_units);
  NumericMatrix I_vax_final(nage_groups, nadmin_units);
  NumericMatrix I_vax2_final(nage_groups, nadmin_units);
  NumericMatrix I_final(nage_groups, nadmin_units);
  
  NumericMatrix I_unvax_final_INCIDENCE(nage_groups, nadmin_units);
  NumericMatrix I_vax_final_INCIDENCE(nage_groups, nadmin_units);
  NumericMatrix I_vax2_final_INCIDENCE(nage_groups, nadmin_units);
  NumericMatrix I_final_INCIDENCE(nage_groups, nadmin_units);
  
  NumericMatrix I_final_counts(nage_groups, nadmin_units);
  NumericMatrix I_final_district_counts(nobs_stratified, nadmin_units);
  
  NumericMatrix I_final_counts_INCIDENCE(nage_groups, nadmin_units);
  NumericMatrix I_final_district_counts_INCIDENCE(nobs_stratified, nadmin_units);
  
  NumericMatrix R_unvax_final(nage_groups, nadmin_units);
  NumericMatrix R_vax_final(nage_groups, nadmin_units);
  NumericMatrix R_vax2_final(nage_groups, nadmin_units);
  NumericMatrix R_final(nage_groups, nadmin_units);
  
  NumericMatrix age_cov_prevalence_counts(nage_groups, nadmin_units);
  NumericMatrix age_cov2_prevalence_counts(nage_groups, nadmin_units);
  
  NumericVector I_annual(40);
  NumericVector adj_cases_annual(40);
  
  
  NumericMatrix estimated_persons_last_week(nage_groups, nadmin_units);
  NumericMatrix this_week_scalar(nage_groups, nadmin_units);
  
  NumericMatrix M_unvax_preage_counts(nage_groups, nadmin_units);
  NumericMatrix M_vax_preage_counts(nage_groups, nadmin_units);
  NumericMatrix M_vax2_preage_counts(nage_groups, nadmin_units);
  
  NumericMatrix S_unvax_preage_counts(nage_groups, nadmin_units);
  NumericMatrix S_vax_preage_counts(nage_groups, nadmin_units);
  NumericMatrix S_vax2_preage_counts(nage_groups, nadmin_units);
  
  NumericMatrix I_unvax_preage_counts(nage_groups, nadmin_units);
  NumericMatrix I_vax_preage_counts(nage_groups, nadmin_units);
  NumericMatrix I_vax2_preage_counts(nage_groups, nadmin_units);
  
  NumericMatrix I_unvax_preage_counts_INCIDENCE(nage_groups, nadmin_units);
  NumericMatrix I_vax_preage_counts_INCIDENCE(nage_groups, nadmin_units);
  NumericMatrix I_vax2_preage_counts_INCIDENCE(nage_groups, nadmin_units);
  
  NumericMatrix R_unvax_preage_counts(nage_groups, nadmin_units);
  NumericMatrix R_vax_preage_counts(nage_groups, nadmin_units);
  NumericMatrix R_vax2_preage_counts(nage_groups, nadmin_units);
  
  NumericMatrix M_unvax_precalib_counts(nage_groups, nadmin_units);
  NumericMatrix M_vax_precalib_counts(nage_groups, nadmin_units);
  NumericMatrix M_vax2_precalib_counts(nage_groups, nadmin_units);
  
  NumericMatrix S_unvax_precalib_counts(nage_groups, nadmin_units);
  NumericMatrix S_vax_precalib_counts(nage_groups, nadmin_units);
  NumericMatrix S_vax2_precalib_counts(nage_groups, nadmin_units);
  
  NumericMatrix I_unvax_precalib_counts(nage_groups, nadmin_units);
  NumericMatrix I_vax_precalib_counts(nage_groups, nadmin_units);
  NumericMatrix I_vax2_precalib_counts(nage_groups, nadmin_units);
  
  NumericMatrix I_unvax_precalib_counts_INCIDENCE(nage_groups, nadmin_units);
  NumericMatrix I_vax_precalib_counts_INCIDENCE(nage_groups, nadmin_units);
  NumericMatrix I_vax2_precalib_counts_INCIDENCE(nage_groups, nadmin_units);
  
  NumericMatrix R_unvax_precalib_counts(nage_groups, nadmin_units);
  NumericMatrix R_vax_precalib_counts(nage_groups, nadmin_units);
  NumericMatrix R_vax2_precalib_counts(nage_groups, nadmin_units);
  
  NumericMatrix M_final_counts(nage_groups,nadmin_units);
  NumericMatrix M_unvax_final_counts(nage_groups,nadmin_units);
  NumericMatrix M_vax_final_counts(nage_groups,nadmin_units);
  NumericMatrix M_vax2_final_counts(nage_groups,nadmin_units);
  
  NumericMatrix S_final_counts(nage_groups,nadmin_units);
  NumericMatrix S_unvax_final_counts(nage_groups,nadmin_units);
  NumericMatrix S_vax_final_counts(nage_groups,nadmin_units);
  NumericMatrix S_vax2_final_counts(nage_groups,nadmin_units);
  
  NumericMatrix I_unvax_final_counts(nage_groups,nadmin_units);
  NumericMatrix I_vax_final_counts(nage_groups,nadmin_units);
  NumericMatrix I_vax2_final_counts(nage_groups,nadmin_units);
  
  NumericMatrix I_unvax_final_counts_INCIDENCE(nage_groups,nadmin_units);
  NumericMatrix I_vax_final_counts_INCIDENCE(nage_groups,nadmin_units);
  NumericMatrix I_vax2_final_counts_INCIDENCE(nage_groups,nadmin_units);
  
  NumericMatrix R_final_counts(nage_groups,nadmin_units);
  NumericMatrix R_unvax_final_counts(nage_groups,nadmin_units);
  NumericMatrix R_vax_final_counts(nage_groups,nadmin_units);
  NumericMatrix R_vax2_final_counts(nage_groups,nadmin_units);
  
  NumericMatrix overall_I_in(nadmin_units,nage_groups);
  NumericMatrix overall_I_in_alpha(nadmin_units,nage_groups);
  
  NumericMatrix I_final_district_counts_0_4(nobs_stratified,nadmin_units);
  NumericMatrix I_final_district_counts_5_9(nobs_stratified,nadmin_units);
  NumericMatrix I_final_district_counts_10_14(nobs_stratified,nadmin_units);
  NumericMatrix I_final_district_counts_15_24(nobs_stratified,nadmin_units);
  NumericMatrix I_final_district_counts_25_34(nobs_stratified,nadmin_units);
  NumericMatrix I_final_district_counts_35_44(nobs_stratified,nadmin_units);
  NumericMatrix I_final_district_counts_45_54(nobs_stratified,nadmin_units);
  NumericMatrix I_final_district_counts_55_64(nobs_stratified,nadmin_units);
  NumericMatrix I_final_district_counts_65_plus(nobs_stratified,nadmin_units);
  
  NumericMatrix I_final_district_counts_0_4_rolling(nobs_stratified,nadmin_units);
  NumericMatrix I_final_district_counts_5_9_rolling(nobs_stratified,nadmin_units);
  NumericMatrix I_final_district_counts_10_14_rolling(nobs_stratified,nadmin_units);
  NumericMatrix I_final_district_counts_15_24_rolling(nobs_stratified,nadmin_units);
  NumericMatrix I_final_district_counts_25_34_rolling(nobs_stratified,nadmin_units);
  NumericMatrix I_final_district_counts_35_44_rolling(nobs_stratified,nadmin_units);
  NumericMatrix I_final_district_counts_45_54_rolling(nobs_stratified,nadmin_units);
  NumericMatrix I_final_district_counts_55_64_rolling(nobs_stratified,nadmin_units);
  NumericMatrix I_final_district_counts_65_plus_rolling(nobs_stratified,nadmin_units);

  NumericMatrix I_quarterly_district_0_4(40,nadmin_units);
  NumericMatrix I_quarterly_district_5_9(40,nadmin_units);
  NumericMatrix I_quarterly_district_10_14(40,nadmin_units);
  NumericMatrix I_quarterly_district_15_24(40,nadmin_units);
  NumericMatrix I_quarterly_district_25_34(40,nadmin_units);
  NumericMatrix I_quarterly_district_35_44(40,nadmin_units);
  NumericMatrix I_quarterly_district_45_54(40,nadmin_units);
  NumericMatrix I_quarterly_district_55_64(40,nadmin_units);
  NumericMatrix I_quarterly_district_65_plus(40,nadmin_units);
  
  NumericMatrix cases_quarterly_district_0_4(40,nadmin_units);
  NumericMatrix cases_quarterly_district_5_9(40,nadmin_units);
  NumericMatrix cases_quarterly_district_10_14(40,nadmin_units);
  NumericMatrix cases_quarterly_district_15_24(40,nadmin_units);
  NumericMatrix cases_quarterly_district_25_34(40,nadmin_units);
  NumericMatrix cases_quarterly_district_35_44(40,nadmin_units);
  NumericMatrix cases_quarterly_district_45_54(40,nadmin_units);
  NumericMatrix cases_quarterly_district_55_64(40,nadmin_units);
  NumericMatrix cases_quarterly_district_65_plus(40,nadmin_units);
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //convert various inputs from vectors to arrays
  vax_efficacy = 0;
  vax_efficacy = exp(logit_vax_efficacy[0]) / (1 + exp(logit_vax_efficacy[0]));
  composite_dose1_efficacy = 0.93 * vax_efficacy;
  composite_dose1_efficacy_1minus = 1 - composite_dose1_efficacy;
  
  composite_dose2_efficacy_1minus = composite_dose1_efficacy_1minus * composite_dose1_efficacy_1minus;
  composite_dose2_efficacy = 1 - composite_dose2_efficacy_1minus;
  
  //rho = exp(logit_rho1[0]) / (1 + exp(logit_rho1[0]));
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ////// Set up starting state
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  
  numerator_temp(0) = 0;
  numerator_temp(0)  = (((national_cases(0) / 53) * (cases_age_bin0_4(0) / all_strat_cases)) / 60) / 0.04742587;
  numerator_temp(1)  = (((national_cases(0) / 53) * (cases_age_bin0_4(0) / all_strat_cases)) / 5) / 0.04742587;
  numerator_temp(2)  = (((national_cases(0) / 53) * ((cases_age_bin5_9(0) )/ all_strat_cases)) ) / 0.04742587;
  numerator_temp(3)  = (((national_cases(0) / 53) * ((cases_age_bin10_14(0) )/ all_strat_cases)) ) / 0.04742587;
  numerator_temp(4)  = (((national_cases(0) / 53) * ((cases_age_bin15_19(0) + cases_age_bin20_24(0) )/ all_strat_cases)) ) / 0.04742587;
  numerator_temp(5)  = (((national_cases(0) / 53) * ((cases_age_bin25_29(0) + cases_age_bin30_34(0) )/ all_strat_cases)) ) / 0.04742587;
  numerator_temp(6)  = (((national_cases(0) / 53) * ((cases_age_bin35_39(0) + cases_age_bin40_44(0) )/ all_strat_cases)) ) / 0.04742587;
  numerator_temp(7)  = (((national_cases(0) / 53) * ((cases_age_bin45_49(0) + cases_age_bin50_54(0) )/ all_strat_cases)) ) / 0.04742587;
  numerator_temp(8)  = (((national_cases(0) / 53) * ((cases_age_bin55_59(0) + cases_age_bin60_64(0) )/ all_strat_cases)) ) / 0.04742587;
  numerator_temp(9)  = (((national_cases(0) / 53) * ((cases_age_bin65_69(0) + cases_age_bin70_74(0) +
                                                       cases_age_bin75_79(0) + cases_age_bin80_84(0) )/ all_strat_cases)) ) / 0.04742587;

  for(int d=0; d<nadmin_units; d++){
    
    for(int a=0; a<12; a++){
      I_unvax_final_counts(a,d) = (numerator_temp(0) / 79);
    }
    for(int a=12; a<16; a++){
      I_unvax_final_counts(a,d) = (numerator_temp(1) / 79);
    }
    I_unvax_final_counts(16,d) = (numerator_temp(2) / 79);
    I_unvax_final_counts(17,d) = (numerator_temp(3) / 79);
    I_unvax_final_counts(18,d) = (numerator_temp(4) / 79);
    I_unvax_final_counts(19,d) = (numerator_temp(5) / 79);
    I_unvax_final_counts(20,d) = (numerator_temp(6) / 79);
    I_unvax_final_counts(21,d) = (numerator_temp(7) / 79);    
    I_unvax_final_counts(22,d) = (numerator_temp(8) / 79);
    I_unvax_final_counts(23,d) = (numerator_temp(9) / 79);
    
    
    for(int a=0; a<6; a++){
      R_unvax_final_counts(a,d) = 0.25 * pop_array_vector[0 + nobs_stratified*a + nobs_stratified*nage_groups*d]; 
    }
    for(int a=6; a<10; a++){
      R_unvax_final_counts(a,d) = 0.60 * pop_array_vector[0 + nobs_stratified*a + nobs_stratified*nage_groups*d];
    }
    for(int a=10; a<nage_groups; a++){
      R_unvax_final_counts(a,d) = 0.98 * pop_array_vector[0 + nobs_stratified*a + nobs_stratified*nage_groups*d];
    }
    
    for(int a=0; a<6; a++){
      M_unvax_final_counts(a,d) =  0.4 * pop_array_vector[0 + nobs_stratified*a + nobs_stratified*nage_groups*d];
    }
    for(int a=6; a<10; a++){
      M_unvax_final_counts(a,d) = 0.1 * pop_array_vector[0 + nobs_stratified*a + nobs_stratified*nage_groups*d];
    }
    for(int a=10; a<nage_groups; a++){
      M_unvax_final_counts(a,d) = 0.000000000000000000001;
    }
    
    
    for(int a = 0; a < nage_groups; a++){
      
      M_vax_final_counts(a,d) = M_vax2_final_counts(a,d) = 0.000000000000000000001;
      I_vax_final_counts(a,d) = I_vax2_final_counts(a,d) = 0.000000000000000000001;
      R_vax_final_counts(a,d) = R_vax2_final_counts(a,d) = 0.000000000000000000001;
      
      M_final_counts(a,d) = M_unvax_final_counts(a,d) + M_vax_final_counts(a,d) + M_vax2_final_counts(a,d);
      I_final_counts(a,d) = I_unvax_final_counts(a,d) + I_vax_final_counts(a,d) + I_vax2_final_counts(a,d);
      R_final_counts(a,d) = R_unvax_final_counts(a,d) + R_vax_final_counts(a,d) + R_vax2_final_counts(a,d);
      
      S_unvax_final_counts(a,d) = pop_array_vector[0 + nobs_stratified*a + nobs_stratified*nage_groups*d] - M_final_counts(a,d)  - I_final_counts(a,d) - R_final_counts(a,d) - 0.000000000000000000001 - 0.000000000000000000001;
      S_vax_final_counts(a,d) = 0.000000000000000000001;
      S_vax2_final_counts(a,d) = 0.000000000000000000001;
      
      S_final_counts(a,d) = S_unvax_final_counts(a,d) + S_vax_final_counts(a,d) + S_vax2_final_counts(a,d);
      
    }
  }
  
  amplitude = (max_beta[0] - min_beta[0]) / 2;
  vertical_shift = (max_beta[0] + min_beta[0]) / 2; 

 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  //// TSIR modelling
  for (int i=1;i<2120; i++){

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //// use beta and combine to make composite beta (with seasonality)
    
    composite_beta = (amplitude * (sin((2*3.14159 * i_ticker_1(i)) + 2))) + vertical_shift;
     
    int t_id;
    int a_size;
    int c_size;
    int d_size;
    int z_size;
    omp_set_num_threads(30);
      #pragma omp parallel
    int t_id = omp_get_thread_num();
    a_size = 24 / omp_get_num_threads();
    d_size = 79 / omp_get_num_threads();
    c_size = 24 / omp_get_num_threads();
    z_size = 79 / omp_get_num_threads();

    
    for (int a = t_id * a_size; a < (t_id+1) * a_size; a ++){
      for (int d = t_id * d_size; d < (t_id+1) * d_size; d ++){
        
        foi(a,d) = 0;
        for (int c = t_id * c_size; c < (t_id+1) * c_size; c ++){
          overall_I_in(d,c) = 0;
          overall_I_in_alpha(d,c) = 0;
          for (int z = t_id * z_size; z < (t_id+1) * z_size; z ++){
            overall_I_in(d,c) +=  M(z,d) * I_final_counts(c,z);
          }
          overall_I_in_alpha(d,c) = pow(overall_I_in(d,c), 0.99);
          foi(a,d) += composite_beta * W_vector[a + nage_groups*c + nage_groups*nage_groups*i + nage_groups*nage_groups*nobs_stratified*d] * (overall_I_in_alpha(d,c) / pop_array_vector[i-1 + nobs_stratified*c + nobs_stratified*nage_groups*d]) ;
        }
        
        // NOTE : W[a][c][i][d] * (overall_I_in[d][i][c] / pop_array[i-1][c][d]) is number of contacts for age a who are infected.
        // NOTE : beta is probability those result in infection
        if(foi(a,d) > 1){
          foi(a,d) = 1;
        }
        if(foi(a,d) < 0){
          foi(a,d) = 0;
        }
        foi2(a,d) = 1 - foi(a,d);
        
      }}
    
      
      
      
      
    int t_id2;
    int a_size2;
    int d_size2;
    //omp_set_num_threads(62);
#pragma omp parallel
    int t_id2 = omp_get_thread_num();
    a_size2 = 24 / omp_get_num_threads();
    d_size2 = 79 / omp_get_num_threads();
  
    for(int a1 = t_id2 * a_size2; a1 < (t_id2+1) * a_size2; a1 ++){
      for(int d1 = t_id2 * d_size2; d1 < (t_id2+1) * d_size2; d1 ++){

        ///////////////////////////////////////////////////////////////////////////////////////////////////////
        //// infect people
        S_unvax_preage_prematernal(a1,d1) = S_unvax_final_counts(a1,d1) * foi2(a1,d1);
        S_vax_preage_prematernal(a1,d1) = S_vax_final_counts(a1,d1) * foi2(a1,d1);
        S_vax2_preage_prematernal(a1,d1) = S_vax2_final_counts(a1,d1) * foi2(a1,d1);

        I_unvax_preage_counts_INCIDENCE(a1,d1) = (S_unvax_final_counts(a1,d1) * (foi(a1,d1))) ;
        I_vax_preage_counts_INCIDENCE(a1,d1) = (S_vax_final_counts(a1,d1) * (foi(a1,d1))) ;
        I_vax2_preage_counts_INCIDENCE(a1,d1) = (S_vax2_final_counts(a1,d1) * (foi(a1,d1))) ;
        
        I_unvax_preage_counts(a1,d1) = (S_unvax_final_counts(a1,d1) * (foi(a1,d1))) + (I_unvax_final_counts(a1,d1) / 2);
        I_vax_preage_counts(a1,d1) = (S_vax_final_counts(a1,d1) * (foi(a1,d1))) + (I_vax_final_counts(a1,d1) / 2);
        I_vax2_preage_counts(a1,d1) = (S_vax2_final_counts(a1,d1) * (foi(a1,d1))) + (I_vax2_final_counts(a1,d1) / 2);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////
        //// allow all people in previously infected state to recover
        R_unvax_preage_counts(a1,d1) = R_unvax_final_counts(a1,d1) + (I_unvax_final_counts(a1,d1) / 2);
        R_vax_preage_counts(a1,d1) = R_vax_final_counts(a1,d1) + (I_vax_final_counts(a1,d1) / 2);
        R_vax2_preage_counts(a1,d1) = R_vax2_final_counts(a1,d1) + (I_vax2_final_counts(a1,d1) / 2);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////
        //// carry over maternal class
        M_unvax_preage_prematernal(a1,d1) = M_unvax_final_counts(a1,d1) ;
        M_vax_preage_prematernal(a1,d1) = M_vax_final_counts(a1,d1) ;
        M_vax2_preage_prematernal(a1,d1) = M_vax2_final_counts(a1,d1) ;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////
        ///// loss of maternal immunity -- on average, takes 3.3 months or about 13 weeks
        maternal_unvax_calc(a1,d1) = M_unvax_preage_prematernal(a1,d1)/13;
        maternal_vax_calc(a1,d1) = M_vax_preage_prematernal(a1,d1)/13;
        maternal_vax2_calc(a1,d1) = M_vax2_preage_prematernal(a1,d1)/13;

        M_unvax_preage_counts(a1,d1) = M_unvax_preage_prematernal(a1,d1) - maternal_unvax_calc(a1,d1);
        M_vax_preage_counts(a1,d1) = M_vax_preage_prematernal(a1,d1) - maternal_vax_calc(a1,d1);
        M_vax2_preage_counts(a1,d1) = M_vax2_preage_prematernal(a1,d1) - maternal_vax2_calc(a1,d1);

        S_unvax_preage_counts(a1,d1) = S_unvax_preage_prematernal(a1,d1) + maternal_unvax_calc(a1,d1);
        S_vax_preage_counts(a1,d1) = S_vax_preage_prematernal(a1,d1) + maternal_vax_calc(a1,d1);
        S_vax2_preage_counts(a1,d1) = S_vax2_preage_prematernal(a1,d1) + maternal_vax2_calc(a1,d1);
      }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //// age everyone

    //// for youngest group, we also need to add the births
    M_unvax_precalib_counts(0,_) = M_unvax_preage_counts(0,_) + births(i,_) - (M_unvax_preage_counts(0,_)/4);
    M_vax_precalib_counts(0,_) = M_vax_preage_counts(0,_)  - (M_vax_preage_counts(0,_)/4);
    M_vax2_precalib_counts(0,_) = M_vax2_preage_counts(0,_)  - (M_vax2_preage_counts(0,_)/4);

    S_unvax_precalib_counts(0,_) = S_unvax_preage_counts(0,_)  - (S_unvax_preage_counts(0,_)/4);
    S_vax_precalib_counts(0,_) = S_vax_preage_counts(0,_)  - (S_vax_preage_counts(0,_)/4);
    S_vax2_precalib_counts(0,_) = S_vax2_preage_counts(0,_)  - (S_vax2_preage_counts(0,_)/4);

    I_unvax_precalib_counts(0,_) = I_unvax_preage_counts(0,_)  - (I_unvax_preage_counts(0,_)/4);
    I_vax_precalib_counts(0,_) = I_vax_preage_counts(0,_)  - (I_vax_preage_counts(0,_)/4);
    I_vax2_precalib_counts(0,_) = I_vax2_preage_counts(0,_)  - (I_vax2_preage_counts(0,_)/4);
    
    I_unvax_precalib_counts_INCIDENCE(0,_) = I_unvax_preage_counts_INCIDENCE(0,_)  - (I_unvax_preage_counts_INCIDENCE(0,_)/4);
    I_vax_precalib_counts_INCIDENCE(0,_) = I_vax_preage_counts_INCIDENCE(0,_)  - (I_vax_preage_counts_INCIDENCE(0,_)/4);
    I_vax2_precalib_counts_INCIDENCE(0,_) = I_vax2_preage_counts_INCIDENCE(0,_)  - (I_vax2_preage_counts_INCIDENCE(0,_)/4);

    R_unvax_precalib_counts(0,_) = R_unvax_preage_counts(0,_)  - (R_unvax_preage_counts(0,_)/4);
    R_vax_precalib_counts(0,_) = R_vax_preage_counts(0,_)  - (R_vax_preage_counts(0,_)/4);
    R_vax2_precalib_counts(0,_) = R_vax2_preage_counts(0,_)  - (R_vax2_preage_counts(0,_)/4);


    for (int a=1; a < 12; a++){

      M_unvax_precalib_counts(a,_) = M_unvax_preage_counts(a,_) + (M_unvax_preage_counts(a-1,_)/4) - (M_unvax_preage_counts(a,_)/4);
      M_vax_precalib_counts(a,_) = M_vax_preage_counts(a,_) + (M_vax_preage_counts(a-1,_)/4) - (M_vax_preage_counts(a,_)/4);
      M_vax2_precalib_counts(a,_) = M_vax2_preage_counts(a,_) + (M_vax2_preage_counts(a-1,_)/4) - (M_vax2_preage_counts(a,_)/4);

      S_unvax_precalib_counts(a,_) = S_unvax_preage_counts(a,_) + (S_unvax_preage_counts(a-1,_)/4) - (S_unvax_preage_counts(a,_)/4);
      S_vax_precalib_counts(a,_) = S_vax_preage_counts(a,_) + (S_vax_preage_counts(a-1,_)/4) - (S_vax_preage_counts(a,_)/4);
      S_vax2_precalib_counts(a,_) = S_vax2_preage_counts(a,_) + (S_vax2_preage_counts(a-1,_)/4) - (S_vax2_preage_counts(a,_)/4);

      I_unvax_precalib_counts(a,_) = I_unvax_preage_counts(a,_) + (I_unvax_preage_counts(a-1,_)/4) - (I_unvax_preage_counts(a,_)/4);
      I_vax_precalib_counts(a,_) = I_vax_preage_counts(a,_) + (I_vax_preage_counts(a-1,_)/4) - (I_vax_preage_counts(a,_)/4);
      I_vax2_precalib_counts(a,_) = I_vax2_preage_counts(a,_) + (I_vax2_preage_counts(a-1,_)/4) - (I_vax2_preage_counts(a,_)/4);

      I_unvax_precalib_counts_INCIDENCE(a,_) = I_unvax_preage_counts_INCIDENCE(a,_) + (I_unvax_preage_counts_INCIDENCE(a-1,_)/4) - (I_unvax_preage_counts_INCIDENCE(a,_)/4);
      I_vax_precalib_counts_INCIDENCE(a,_) = I_vax_preage_counts_INCIDENCE(a,_) + (I_vax_preage_counts_INCIDENCE(a-1,_)/4) - (I_vax_preage_counts_INCIDENCE(a,_)/4);
      I_vax2_precalib_counts_INCIDENCE(a,_) = I_vax2_preage_counts_INCIDENCE(a,_) + (I_vax2_preage_counts_INCIDENCE(a-1,_)/4) - (I_vax2_preage_counts_INCIDENCE(a,_)/4);
      
      R_unvax_precalib_counts(a,_) = R_unvax_preage_counts(a,_) + (R_unvax_preage_counts(a-1,_)/4) - (R_unvax_preage_counts(a,_)/4);
      R_vax_precalib_counts(a,_) = R_vax_preage_counts(a,_) + (R_vax_preage_counts(a-1,_)/4) - (R_vax_preage_counts(a,_)/4);
      R_vax2_precalib_counts(a,_) = R_vax2_preage_counts(a,_) + (R_vax2_preage_counts(a-1,_)/4) - (R_vax2_preage_counts(a,_)/4);

    }

    M_unvax_precalib_counts(12,_) = M_unvax_preage_counts(12,_) + (M_unvax_preage_counts(11,_)/4) - (M_unvax_preage_counts(12,_)/53);
    M_vax_precalib_counts(12,_) = M_vax_preage_counts(12,_) + (M_vax_preage_counts(11,_)/4) - (M_vax_preage_counts(12,_)/53);
    M_vax2_precalib_counts(12,_) = M_vax2_preage_counts(12,_) + (M_vax2_preage_counts(11,_)/4) - (M_vax2_preage_counts(12,_)/53);

    S_unvax_precalib_counts(12,_) = S_unvax_preage_counts(12,_) + (S_unvax_preage_counts(11,_)/4) - (S_unvax_preage_counts(12,_)/53);
    S_vax_precalib_counts(12,_) = S_vax_preage_counts(12,_) + (S_vax_preage_counts(11,_)/4) - (S_vax_preage_counts(12,_)/53);
    S_vax2_precalib_counts(12,_) = S_vax2_preage_counts(12,_) + (S_vax2_preage_counts(11,_)/4) - (S_vax2_preage_counts(12,_)/53);

    I_unvax_precalib_counts(12,_) = I_unvax_preage_counts(12,_) + (I_unvax_preage_counts(11,_)/4) - (I_unvax_preage_counts(12,_)/53);
    I_vax_precalib_counts(12,_) = I_vax_preage_counts(12,_) + (I_vax_preage_counts(11,_)/4) - (I_vax_preage_counts(12,_)/53);
    I_vax2_precalib_counts(12,_) = I_vax2_preage_counts(12,_) + (I_vax2_preage_counts(11,_)/4) - (I_vax2_preage_counts(12,_)/53);

    I_unvax_precalib_counts_INCIDENCE(12,_) = I_unvax_preage_counts_INCIDENCE(12,_) + (I_unvax_preage_counts_INCIDENCE(11,_)/4) - (I_unvax_preage_counts_INCIDENCE(12,_)/53);
    I_vax_precalib_counts_INCIDENCE(12,_) = I_vax_preage_counts_INCIDENCE(12,_) + (I_vax_preage_counts_INCIDENCE(11,_)/4) - (I_vax_preage_counts_INCIDENCE(12,_)/53);
    I_vax2_precalib_counts_INCIDENCE(12,_) = I_vax2_preage_counts_INCIDENCE(12,_) + (I_vax2_preage_counts_INCIDENCE(11,_)/4) - (I_vax2_preage_counts_INCIDENCE(12,_)/53);
    
    R_unvax_precalib_counts(12,_) = R_unvax_preage_counts(12,_) + (R_unvax_preage_counts(11,_)/4) - (R_unvax_preage_counts(12,_)/53);
    R_vax_precalib_counts(12,_) = R_vax_preage_counts(12,_) + (R_vax_preage_counts(11,_)/4) - (R_vax_preage_counts(12,_)/53);
    R_vax2_precalib_counts(12,_) = R_vax2_preage_counts(12,_) + (R_vax2_preage_counts(11,_)/4) - (R_vax2_preage_counts(12,_)/53);

    for (int a = 13; a < 16; a++){

      M_unvax_precalib_counts(a,_) = M_unvax_preage_counts(a,_) + (M_unvax_preage_counts(a-1,_)/53) - (M_unvax_preage_counts(a,_)/53);
      M_vax_precalib_counts(a,_) = M_vax_preage_counts(a,_) + (M_vax_preage_counts(a-1,_)/53) - (M_vax_preage_counts(a,_)/53);
      M_vax2_precalib_counts(a,_) = M_vax2_preage_counts(a,_) + (M_vax2_preage_counts(a-1,_)/53) - (M_vax2_preage_counts(a,_)/53);

      S_unvax_precalib_counts(a,_) = S_unvax_preage_counts(a,_) + (S_unvax_preage_counts(a-1,_)/53) - (S_unvax_preage_counts(a,_)/53);
      S_vax_precalib_counts(a,_) = S_vax_preage_counts(a,_) + (S_vax_preage_counts(a-1,_)/53) - (S_vax_preage_counts(a,_)/53);
      S_vax2_precalib_counts(a,_) = S_vax2_preage_counts(a,_) + (S_vax2_preage_counts(a-1,_)/53) - (S_vax2_preage_counts(a,_)/53);

      I_unvax_precalib_counts(a,_) = I_unvax_preage_counts(a,_) + (I_unvax_preage_counts(a-1,_)/53) - (I_unvax_preage_counts(a,_)/53);
      I_vax_precalib_counts(a,_) = I_vax_preage_counts(a,_) + (I_vax_preage_counts(a-1,_)/53) - (I_vax_preage_counts(a,_)/53);
      I_vax2_precalib_counts(a,_) = I_vax2_preage_counts(a,_) + (I_vax2_preage_counts(a-1,_)/53) - (I_vax2_preage_counts(a,_)/53);

      I_unvax_precalib_counts_INCIDENCE(a,_) = I_unvax_preage_counts_INCIDENCE(a,_) + (I_unvax_preage_counts_INCIDENCE(a-1,_)/53) - (I_unvax_preage_counts_INCIDENCE(a,_)/53);
      I_vax_precalib_counts_INCIDENCE(a,_) = I_vax_preage_counts_INCIDENCE(a,_) + (I_vax_preage_counts_INCIDENCE(a-1,_)/53) - (I_vax_preage_counts_INCIDENCE(a,_)/53);
      I_vax2_precalib_counts_INCIDENCE(a,_) = I_vax2_preage_counts_INCIDENCE(a,_) + (I_vax2_preage_counts_INCIDENCE(a-1,_)/53) - (I_vax2_preage_counts_INCIDENCE(a,_)/53);
      
      R_unvax_precalib_counts(a,_) = R_unvax_preage_counts(a,_) + (R_unvax_preage_counts(a-1,_)/53) - (R_unvax_preage_counts(a,_)/53);
      R_vax_precalib_counts(a,_) = R_vax_preage_counts(a,_) + (R_vax_preage_counts(a-1,_)/53) - (R_vax_preage_counts(a,_)/53);
      R_vax2_precalib_counts(a,_) = R_vax2_preage_counts(a,_) + (R_vax2_preage_counts(a-1,_)/53) - (R_vax2_preage_counts(a,_)/53);

    }
    
    M_unvax_precalib_counts(16,_) = M_unvax_preage_counts(16,_) + (M_unvax_preage_counts(15,_)/53) - (M_unvax_preage_counts(16,_)/265);
    M_vax_precalib_counts(16,_) = M_vax_preage_counts(16,_) + (M_vax_preage_counts(15,_)/53) - (M_vax_preage_counts(16,_)/265);
    M_vax2_precalib_counts(16,_) = M_vax2_preage_counts(16,_) + (M_vax2_preage_counts(15,_)/53) - (M_vax2_preage_counts(16,_)/265);
    
    S_unvax_precalib_counts(16,_) = S_unvax_preage_counts(16,_) + (S_unvax_preage_counts(15,_)/53) - (S_unvax_preage_counts(16,_)/265);
    S_vax_precalib_counts(16,_) = S_vax_preage_counts(16,_) + (S_vax_preage_counts(15,_)/53) - (S_vax_preage_counts(16,_)/265);
    S_vax2_precalib_counts(16,_) = S_vax2_preage_counts(16,_) + (S_vax2_preage_counts(15,_)/53) - (S_vax2_preage_counts(16,_)/265);
    
    I_unvax_precalib_counts(16,_) = I_unvax_preage_counts(16,_) + (I_unvax_preage_counts(15,_)/53) - (I_unvax_preage_counts(16,_)/265);
    I_vax_precalib_counts(16,_) = I_vax_preage_counts(16,_) + (I_vax_preage_counts(15,_)/53) - (I_vax_preage_counts(16,_)/265);
    I_vax2_precalib_counts(16,_) = I_vax2_preage_counts(16,_) + (I_vax2_preage_counts(15,_)/53) - (I_vax2_preage_counts(16,_)/265);

    I_unvax_precalib_counts_INCIDENCE(16,_) = I_unvax_preage_counts_INCIDENCE(16,_) + (I_unvax_preage_counts_INCIDENCE(15,_)/53) - (I_unvax_preage_counts_INCIDENCE(16,_)/265);
    I_vax_precalib_counts_INCIDENCE(16,_) = I_vax_preage_counts_INCIDENCE(16,_) + (I_vax_preage_counts_INCIDENCE(15,_)/53) - (I_vax_preage_counts_INCIDENCE(16,_)/265);
    I_vax2_precalib_counts_INCIDENCE(16,_) = I_vax2_preage_counts_INCIDENCE(16,_) + (I_vax2_preage_counts_INCIDENCE(15,_)/53) - (I_vax2_preage_counts_INCIDENCE(16,_)/265);
        
    R_unvax_precalib_counts(16,_) = R_unvax_preage_counts(16,_) + (R_unvax_preage_counts(15,_)/53) - (R_unvax_preage_counts(16,_)/265);
    R_vax_precalib_counts(16,_) = R_vax_preage_counts(16,_) + (R_vax_preage_counts(15,_)/53) - (R_vax_preage_counts(16,_)/265);
    R_vax2_precalib_counts(16,_) = R_vax2_preage_counts(16,_) + (R_vax2_preage_counts(15,_)/53) - (R_vax2_preage_counts(16,_)/265);
    
    
    M_unvax_precalib_counts(17,_) = M_unvax_preage_counts(17,_) + (M_unvax_preage_counts(16,_)/265) - (M_unvax_preage_counts(17,_)/265);
    M_vax_precalib_counts(17,_) = M_vax_preage_counts(17,_) + (M_vax_preage_counts(16,_)/265) - (M_vax_preage_counts(17,_)/265);
    M_vax2_precalib_counts(17,_) = M_vax2_preage_counts(17,_) + (M_vax2_preage_counts(16,_)/265) - (M_vax2_preage_counts(17,_)/265);
    
    S_unvax_precalib_counts(17,_) = S_unvax_preage_counts(17,_) + (S_unvax_preage_counts(16,_)/265) - (S_unvax_preage_counts(17,_)/265);
    S_vax_precalib_counts(17,_) = S_vax_preage_counts(17,_) + (S_vax_preage_counts(16,_)/265) - (S_vax_preage_counts(17,_)/265);
    S_vax2_precalib_counts(17,_) = S_vax2_preage_counts(17,_) + (S_vax2_preage_counts(16,_)/265) - (S_vax2_preage_counts(17,_)/265);
    
    I_unvax_precalib_counts(17,_) = I_unvax_preage_counts(17,_) + (I_unvax_preage_counts(16,_)/265) - (I_unvax_preage_counts(17,_)/265);
    I_vax_precalib_counts(17,_) = I_vax_preage_counts(17,_) + (I_vax_preage_counts(16,_)/265) - (I_vax_preage_counts(17,_)/265);
    I_vax2_precalib_counts(17,_) = I_vax2_preage_counts(17,_) + (I_vax2_preage_counts(16,_)/265) - (I_vax2_preage_counts(17,_)/265);

    I_unvax_precalib_counts_INCIDENCE(17,_) = I_unvax_preage_counts_INCIDENCE(17,_) + (I_unvax_preage_counts_INCIDENCE(16,_)/265) - (I_unvax_preage_counts_INCIDENCE(17,_)/265);
    I_vax_precalib_counts_INCIDENCE(17,_) = I_vax_preage_counts_INCIDENCE(17,_) + (I_vax_preage_counts_INCIDENCE(16,_)/265) - (I_vax_preage_counts_INCIDENCE(17,_)/265);
    I_vax2_precalib_counts_INCIDENCE(17,_) = I_vax2_preage_counts_INCIDENCE(17,_) + (I_vax2_preage_counts_INCIDENCE(16,_)/265) - (I_vax2_preage_counts_INCIDENCE(17,_)/265);
    
    R_unvax_precalib_counts(17,_) = R_unvax_preage_counts(17,_) + (R_unvax_preage_counts(16,_)/265) - (R_unvax_preage_counts(17,_)/265);
    R_vax_precalib_counts(17,_) = R_vax_preage_counts(17,_) + (R_vax_preage_counts(16,_)/265) - (R_vax_preage_counts(17,_)/265);
    R_vax2_precalib_counts(17,_) = R_vax2_preage_counts(17,_) + (R_vax2_preage_counts(16,_)/265) - (R_vax2_preage_counts(17,_)/265);
    
    
    M_unvax_precalib_counts(18,_) = M_unvax_preage_counts(18,_) + (M_unvax_preage_counts(17,_)/265) - (M_unvax_preage_counts(18,_)/530);
    M_vax_precalib_counts(18,_) = M_vax_preage_counts(18,_) + (M_vax_preage_counts(17,_)/265) - (M_vax_preage_counts(18,_)/530);
    M_vax2_precalib_counts(18,_) = M_vax2_preage_counts(18,_) + (M_vax2_preage_counts(17,_)/265) - (M_vax2_preage_counts(18,_)/530);
    
    S_unvax_precalib_counts(18,_) = S_unvax_preage_counts(18,_) + (S_unvax_preage_counts(17,_)/265) - (S_unvax_preage_counts(18,_)/530);
    S_vax_precalib_counts(18,_) = S_vax_preage_counts(18,_) + (S_vax_preage_counts(17,_)/265) - (S_vax_preage_counts(18,_)/530);
    S_vax2_precalib_counts(18,_) = S_vax2_preage_counts(18,_) + (S_vax2_preage_counts(17,_)/265) - (S_vax2_preage_counts(18,_)/530);
    
    I_unvax_precalib_counts(18,_) = I_unvax_preage_counts(18,_) + (I_unvax_preage_counts(17,_)/265) - (I_unvax_preage_counts(18,_)/530);
    I_vax_precalib_counts(18,_) = I_vax_preage_counts(18,_) + (I_vax_preage_counts(17,_)/265) - (I_vax_preage_counts(18,_)/530);
    I_vax2_precalib_counts(18,_) = I_vax2_preage_counts(18,_) + (I_vax2_preage_counts(17,_)/265) - (I_vax2_preage_counts(18,_)/530);

    I_unvax_precalib_counts_INCIDENCE(18,_) = I_unvax_preage_counts_INCIDENCE(18,_) + (I_unvax_preage_counts_INCIDENCE(17,_)/265) - (I_unvax_preage_counts_INCIDENCE(18,_)/530);
    I_vax_precalib_counts_INCIDENCE(18,_) = I_vax_preage_counts_INCIDENCE(18,_) + (I_vax_preage_counts_INCIDENCE(17,_)/265) - (I_vax_preage_counts_INCIDENCE(18,_)/530);
    I_vax2_precalib_counts_INCIDENCE(18,_) = I_vax2_preage_counts_INCIDENCE(18,_) + (I_vax2_preage_counts_INCIDENCE(17,_)/265) - (I_vax2_preage_counts_INCIDENCE(18,_)/530);
    
    R_unvax_precalib_counts(18,_) = R_unvax_preage_counts(18,_) + (R_unvax_preage_counts(17,_)/265) - (R_unvax_preage_counts(18,_)/530);
    R_vax_precalib_counts(18,_) = R_vax_preage_counts(18,_) + (R_vax_preage_counts(17,_)/265) - (R_vax_preage_counts(18,_)/530);
    R_vax2_precalib_counts(18,_) = R_vax2_preage_counts(18,_) + (R_vax2_preage_counts(17,_)/265) - (R_vax2_preage_counts(18,_)/530);

    for (int a = 19; a < 23; a++){

      M_unvax_precalib_counts(a,_) = M_unvax_preage_counts(a,_) + (M_unvax_preage_counts(a-1,_)/530) - (M_unvax_preage_counts(a,_)/530);
      M_vax_precalib_counts(a,_) = M_vax_preage_counts(a,_) + (M_vax_preage_counts(a-1,_)/530) - (M_vax_preage_counts(a,_)/530);
      M_vax2_precalib_counts(a,_) = M_vax2_preage_counts(a,_) + (M_vax2_preage_counts(a-1,_)/530) - (M_vax2_preage_counts(a,_)/530);

      S_unvax_precalib_counts(a,_) = S_unvax_preage_counts(a,_) + (S_unvax_preage_counts(a-1,_)/530) - (S_unvax_preage_counts(a,_)/530);
      S_vax_precalib_counts(a,_) = S_vax_preage_counts(a,_) + (S_vax_preage_counts(a-1,_)/530) - (S_vax_preage_counts(a,_)/530);
      S_vax2_precalib_counts(a,_) = S_vax2_preage_counts(a,_) + (S_vax2_preage_counts(a-1,_)/530) - (S_vax2_preage_counts(a,_)/530);

      I_unvax_precalib_counts(a,_) = I_unvax_preage_counts(a,_) + (I_unvax_preage_counts(a-1,_)/530) - (I_unvax_preage_counts(a,_)/530);
      I_vax_precalib_counts(a,_) = I_vax_preage_counts(a,_) + (I_vax_preage_counts(a-1,_)/530) - (I_vax_preage_counts(a,_)/530);
      I_vax2_precalib_counts(a,_) = I_vax2_preage_counts(a,_) + (I_vax2_preage_counts(a-1,_)/530) - (I_vax2_preage_counts(a,_)/530);

      I_unvax_precalib_counts_INCIDENCE(a,_) = I_unvax_preage_counts_INCIDENCE(a,_) + (I_unvax_preage_counts_INCIDENCE(a-1,_)/530) - (I_unvax_preage_counts_INCIDENCE(a,_)/530);
      I_vax_precalib_counts_INCIDENCE(a,_) = I_vax_preage_counts_INCIDENCE(a,_) + (I_vax_preage_counts_INCIDENCE(a-1,_)/530) - (I_vax_preage_counts_INCIDENCE(a,_)/530);
      I_vax2_precalib_counts_INCIDENCE(a,_) = I_vax2_preage_counts_INCIDENCE(a,_) + (I_vax2_preage_counts_INCIDENCE(a-1,_)/530) - (I_vax2_preage_counts_INCIDENCE(a,_)/530);
      
      R_unvax_precalib_counts(a,_) = R_unvax_preage_counts(a,_) + (R_unvax_preage_counts(a-1,_)/530) - (R_unvax_preage_counts(a,_)/530);
      R_vax_precalib_counts(a,_) = R_vax_preage_counts(a,_) + (R_vax_preage_counts(a-1,_)/530) - (R_vax_preage_counts(a,_)/530);
      R_vax2_precalib_counts(a,_) = R_vax2_preage_counts(a,_) + (R_vax2_preage_counts(a-1,_)/530) - (R_vax2_preage_counts(a,_)/530);

    }

    M_unvax_precalib_counts(23,_) = M_unvax_preage_counts(23,_) + (M_unvax_preage_counts(22,_)/530);
    M_vax_precalib_counts(23,_) = M_vax_preage_counts(23,_) + (M_vax_preage_counts(22,_)/530);
    M_vax2_precalib_counts(23,_) = M_vax2_preage_counts(23,_) + (M_vax2_preage_counts(22,_)/530);

    S_unvax_precalib_counts(23,_) = S_unvax_preage_counts(23,_) + (S_unvax_preage_counts(22,_)/530);
    S_vax_precalib_counts(23,_) = S_vax_preage_counts(23,_) + (S_vax_preage_counts(22,_)/530);
    S_vax2_precalib_counts(23,_) = S_vax2_preage_counts(23,_) + (S_vax2_preage_counts(22,_)/530);

    I_unvax_precalib_counts(23,_) = I_unvax_preage_counts(23,_) + (I_unvax_preage_counts(22,_)/530);
    I_vax_precalib_counts(23,_) = I_vax_preage_counts(23,_) + (I_vax_preage_counts(22,_)/530);
    I_vax2_precalib_counts(23,_) = I_vax2_preage_counts(23,_) + (I_vax2_preage_counts(22,_)/530);

    I_unvax_precalib_counts_INCIDENCE(23,_) = I_unvax_preage_counts_INCIDENCE(23,_) + (I_unvax_preage_counts_INCIDENCE(22,_)/530);
    I_vax_precalib_counts_INCIDENCE(23,_) = I_vax_preage_counts_INCIDENCE(23,_) + (I_vax_preage_counts_INCIDENCE(22,_)/530);
    I_vax2_precalib_counts_INCIDENCE(23,_) = I_vax2_preage_counts_INCIDENCE(23,_) + (I_vax2_preage_counts_INCIDENCE(22,_)/530);
    
    R_unvax_precalib_counts(23,_) = R_unvax_preage_counts(23,_) + (R_unvax_preage_counts(22,_)/530);
    R_vax_precalib_counts(23,_) = R_vax_preage_counts(23,_) + (R_vax_preage_counts(22,_)/530);
    R_vax2_precalib_counts(23,_) = R_vax2_preage_counts(23,_) + (R_vax2_preage_counts(22,_)/530);


    
    int t_id3;
    int a_size3;
    int d_size3;
    //omp_set_num_threads(62);
    #pragma omp parallel
    int t_id3 = omp_get_thread_num();
    a_size3 = 24 / omp_get_num_threads();
    d_size3 = 79 / omp_get_num_threads();

    for(int d2 = t_id3 * d_size3; d2 < (t_id3+1) * d_size3; d2 ++){    
      for(int a2 = t_id3 * a_size3; a2 < (t_id3+1) * a_size3; a2 ++){

    // for(int a2=0; a2 < nage_groups; a2++){
    //   for(int d2=0; d2< nadmin_units; d2++){
    // 
    //     ///////////////////////////////////////////////////////////////////////////////////////////////////////
        //// scale populations to this week's size
        estimated_persons_last_week(a2,d2) =
          M_unvax_precalib_counts(a2,d2) + S_unvax_precalib_counts(a2,d2) + I_unvax_precalib_counts(a2,d2) + R_unvax_precalib_counts(a2,d2) +
          M_vax_precalib_counts(a2,d2) + S_vax_precalib_counts(a2,d2) + I_vax_precalib_counts(a2,d2) + R_vax_precalib_counts(a2,d2) +
          M_vax2_precalib_counts(a2,d2) + S_vax2_precalib_counts(a2,d2) + I_vax2_precalib_counts(a2,d2) + R_vax2_precalib_counts(a2,d2);

        this_week_scalar(a2,d2) = pop_array_vector[i + nobs_stratified*a2 + nobs_stratified*nage_groups*d2] / estimated_persons_last_week(a2,d2);

        M_unvax_prevax_counts(a2,d2) = M_unvax_precalib_counts(a2,d2) * this_week_scalar(a2,d2);
        M_vax_prevax_counts(a2,d2) = M_vax_precalib_counts(a2,d2) * this_week_scalar(a2,d2);
        M_vax2_prevax_counts(a2,d2) = M_vax2_precalib_counts(a2,d2) * this_week_scalar(a2,d2);

        S_unvax_prevax_counts(a2,d2) = S_unvax_precalib_counts(a2,d2) * this_week_scalar(a2,d2);
        S_vax_prevax_counts(a2,d2) = S_vax_precalib_counts(a2,d2) * this_week_scalar(a2,d2);
        S_vax2_prevax_counts(a2,d2) = S_vax2_precalib_counts(a2,d2) * this_week_scalar(a2,d2);

        I_unvax_prevax_counts(a2,d2) = I_unvax_precalib_counts(a2,d2) * this_week_scalar(a2,d2);
        I_vax_prevax_counts(a2,d2) = I_vax_precalib_counts(a2,d2) * this_week_scalar(a2,d2);
        I_vax2_prevax_counts(a2,d2) = I_vax2_precalib_counts(a2,d2) * this_week_scalar(a2,d2);
        
        I_unvax_prevax_counts_INCIDENCE(a2,d2) = I_unvax_precalib_counts_INCIDENCE(a2,d2) * this_week_scalar(a2,d2);
        I_vax_prevax_counts_INCIDENCE(a2,d2) = I_vax_precalib_counts_INCIDENCE(a2,d2) * this_week_scalar(a2,d2);
        I_vax2_prevax_counts_INCIDENCE(a2,d2) = I_vax2_precalib_counts_INCIDENCE(a2,d2) * this_week_scalar(a2,d2);

        R_unvax_prevax_counts(a2,d2) = R_unvax_precalib_counts(a2,d2) * this_week_scalar(a2,d2);
        R_vax_prevax_counts(a2,d2) = R_vax_precalib_counts(a2,d2) * this_week_scalar(a2,d2);
        R_vax2_prevax_counts(a2,d2) = R_vax2_precalib_counts(a2,d2) * this_week_scalar(a2,d2);

        /////////////////////////////////////////////////////////////////////////////////////////////////////
        //// vaccinate persons
        unvax_denom(a2,d2) = (M_unvax_prevax_counts(a2,d2) + S_unvax_prevax_counts(a2,d2) + I_unvax_prevax_counts(a2,d2) +  R_unvax_prevax_counts(a2,d2));

        temp_vax_calc(a2,d2) = (age_cov_prevalence_counts_vector[i + nobs_stratified*a2 + nobs_stratified*nage_groups*d2] -
          M_vax_prevax_counts(a2,d2) - M_vax2_prevax_counts(a2,d2) -
          S_vax_prevax_counts(a2,d2) - S_vax2_prevax_counts(a2,d2) -
          I_vax_prevax_counts(a2,d2) - I_vax2_prevax_counts(a2,d2) -
          R_vax_prevax_counts(a2,d2) - R_vax2_prevax_counts(a2,d2));

        vax2_calc(a2,d2) = (age_cov2_prevalence_counts_vector[i + nobs_stratified*a2 + nobs_stratified*nage_groups*d2] -
          M_vax2_prevax_counts(a2,d2) -
          S_vax2_prevax_counts(a2,d2) -
          I_vax2_prevax_counts(a2,d2) -
          R_vax2_prevax_counts(a2,d2));

        if(temp_vax_calc(a2,d2) > 0){

          if(unvax_denom(a2,d2) > 0){

            M_prop(a2,d2) = M_unvax_prevax_counts(a2,d2) / unvax_denom(a2,d2);
            S_prop(a2,d2) = S_unvax_prevax_counts(a2,d2) / unvax_denom(a2,d2);
            I_prop(a2,d2) = I_unvax_prevax_counts(a2,d2) / unvax_denom(a2,d2);
            R_prop(a2,d2) = R_unvax_prevax_counts(a2,d2) / unvax_denom(a2,d2);

          }else{
            M_prop(a2,d2) = 0;
            S_prop(a2,d2) = 0;
            I_prop(a2,d2) = 0;
            R_prop(a2,d2) = 0;
          }

          temp_calc_efficacy_1minus(a2,d2) = temp_vax_calc(a2,d2) * composite_dose1_efficacy_1minus;
          temp_calc_efficacy(a2,d2) = temp_vax_calc(a2,d2) * composite_dose1_efficacy;

          M_unvax_postvax_counts(a2,d2) = M_unvax_prevax_counts(a2,d2) - (M_prop(a2,d2) * temp_vax_calc(a2,d2));
          S_unvax_postvax_counts(a2,d2) = S_unvax_prevax_counts(a2,d2) - (S_prop(a2,d2) * temp_vax_calc(a2,d2));
          I_unvax_postvax_counts(a2,d2) = I_unvax_prevax_counts(a2,d2) - (I_prop(a2,d2) * temp_vax_calc(a2,d2));
          R_unvax_postvax_counts(a2,d2) = R_unvax_prevax_counts(a2,d2) - (R_prop(a2,d2) * temp_vax_calc(a2,d2));

          M_vax_postvax_counts(a2,d2) = M_vax_prevax_counts(a2,d2) + (M_prop(a2,d2) * temp_calc_efficacy_1minus(a2,d2));
          S_vax_postvax_counts(a2,d2) = S_vax_prevax_counts(a2,d2) + (S_prop(a2,d2) * temp_calc_efficacy_1minus(a2,d2));
          I_vax_postvax_counts(a2,d2) = I_vax_prevax_counts(a2,d2) + (I_prop(a2,d2) * temp_calc_efficacy_1minus(a2,d2));
          R_vax_postvax_counts(a2,d2) = R_vax_prevax_counts(a2,d2) + (R_prop(a2,d2) * temp_vax_calc(a2,d2)) +
            (M_prop(a2,d2) * temp_calc_efficacy(a2,d2)) +
            (S_prop(a2,d2) * temp_calc_efficacy(a2,d2)) +
            (I_prop(a2,d2) * temp_calc_efficacy(a2,d2)) ;

        }else{
          M_unvax_postvax_counts(a2,d2) = M_unvax_prevax_counts(a2,d2);
          S_unvax_postvax_counts(a2,d2) = S_unvax_prevax_counts(a2,d2);
          I_unvax_postvax_counts(a2,d2) = I_unvax_prevax_counts(a2,d2);
          R_unvax_postvax_counts(a2,d2) = R_unvax_prevax_counts(a2,d2);

          M_vax_postvax_counts(a2,d2) = M_vax_prevax_counts(a2,d2);
          S_vax_postvax_counts(a2,d2) = S_vax_prevax_counts(a2,d2);
          I_vax_postvax_counts(a2,d2) = I_vax_prevax_counts(a2,d2);
          R_vax_postvax_counts(a2,d2) = R_vax_prevax_counts(a2,d2);
        }

          if(vax2_calc(a2,d2) > 0){

            vax_denom(a2,d2) = M_vax_postvax_counts(a2,d2) + S_vax_postvax_counts(a2,d2) + I_vax_postvax_counts(a2,d2) + R_vax_postvax_counts(a2,d2);

            if(vax_denom(a2,d2) > 0){
              M_prop2(a2,d2) = M_vax_postvax_counts(a2,d2) / vax_denom(a2,d2);
              S_prop2(a2,d2) = S_vax_postvax_counts(a2,d2) / vax_denom(a2,d2);
              I_prop2(a2,d2) = I_vax_postvax_counts(a2,d2) / vax_denom(a2,d2);
              R_prop2(a2,d2) = R_vax_postvax_counts(a2,d2) / vax_denom(a2,d2);
            }else{
              M_prop2(a2,d2) = 0;
              S_prop2(a2,d2) = 0;
              I_prop2(a2,d2) = 0;
              R_prop2(a2,d2) = 0;
            }

            temp_calc2_efficacy_1minus(a2,d2) = vax2_calc(a2,d2) * composite_dose1_efficacy_1minus;
            temp_calc2_efficacy(a2,d2) = vax2_calc(a2,d2) * composite_dose1_efficacy;

            // taking away all the ones that need vaccinating with second dose, proportional to first dose coverage by health state
            M_vax_postvax2_counts(a2,d2) = M_vax_postvax_counts(a2,d2) - (M_prop2(a2,d2) * vax2_calc(a2,d2));
            S_vax_postvax2_counts(a2,d2) = S_vax_postvax_counts(a2,d2) - (S_prop2(a2,d2) * vax2_calc(a2,d2));
            I_vax_postvax2_counts(a2,d2) = I_vax_postvax_counts(a2,d2) - (I_prop2(a2,d2) * vax2_calc(a2,d2));
            R_vax_postvax2_counts(a2,d2) = R_vax_postvax_counts(a2,d2) - (R_prop2(a2,d2) * vax2_calc(a2,d2));

            // adding second dose coverage
            M_vax2_postvax2_counts(a2,d2) = M_vax2_prevax_counts(a2,d2) + (M_prop2(a2,d2) * temp_calc2_efficacy_1minus(a2,d2));
            S_vax2_postvax2_counts(a2,d2) = S_vax2_prevax_counts(a2,d2) + (S_prop2(a2,d2) * temp_calc2_efficacy_1minus(a2,d2));
            I_vax2_postvax2_counts(a2,d2) = I_vax2_prevax_counts(a2,d2) + (I_prop2(a2,d2) * temp_calc2_efficacy_1minus(a2,d2));
            R_vax2_postvax2_counts(a2,d2) = R_vax2_prevax_counts(a2,d2) +
              (R_prop2(a2,d2) * temp_calc2_efficacy_1minus(a2,d2)) +
              (R_prop2(a2,d2) * temp_calc2_efficacy(a2,d2)) +
              (M_prop2(a2,d2) * temp_calc2_efficacy(a2,d2)) +
              (S_prop2(a2,d2) * temp_calc2_efficacy(a2,d2)) +
              (I_prop2(a2,d2) * temp_calc2_efficacy(a2,d2)) ;

          }else{
            M_vax_postvax2_counts(a2,d2) = M_vax_postvax_counts(a2,d2);
            S_vax_postvax2_counts(a2,d2) = S_vax_postvax_counts(a2,d2);
            I_vax_postvax2_counts(a2,d2) = I_vax_postvax_counts(a2,d2);
            R_vax_postvax2_counts(a2,d2) = R_vax_postvax_counts(a2,d2);

            M_vax2_postvax2_counts(a2,d2) = M_vax2_prevax_counts(a2,d2);
            S_vax2_postvax2_counts(a2,d2) = S_vax2_prevax_counts(a2,d2);
            I_vax2_postvax2_counts(a2,d2) = I_vax2_prevax_counts(a2,d2);
            R_vax2_postvax2_counts(a2,d2) = R_vax2_prevax_counts(a2,d2);
          }

          M_unvax_final_counts(a2,d2) = M_unvax_postvax_counts(a2,d2);
          M_vax_final_counts(a2,d2) = M_vax_postvax2_counts(a2,d2) ;
          M_vax2_final_counts(a2,d2)  =  M_vax2_postvax2_counts(a2,d2) ;

          S_unvax_final_counts(a2,d2)  = S_unvax_postvax_counts(a2,d2)  ;
          S_vax_final_counts(a2,d2)  = S_vax_postvax2_counts(a2,d2)  ;
          S_vax2_final_counts(a2,d2)  = S_vax2_postvax2_counts(a2,d2) ;

          I_unvax_final_counts(a2,d2)  = I_unvax_postvax_counts(a2,d2)  ;
          I_vax_final_counts(a2,d2)  = I_vax_postvax2_counts(a2,d2)  ;
          I_vax2_final_counts(a2,d2)  = I_vax2_postvax2_counts(a2,d2) ;
          
          I_unvax_final_counts_INCIDENCE(a2,d2)  = I_unvax_prevax_counts_INCIDENCE(a2,d2)  ;
          I_vax_final_counts_INCIDENCE(a2,d2)  = I_vax_prevax_counts_INCIDENCE(a2,d2)  ;
          I_vax2_final_counts_INCIDENCE(a2,d2)  = I_vax2_prevax_counts_INCIDENCE(a2,d2) ;

          R_unvax_final_counts(a2,d2) = R_unvax_postvax_counts(a2,d2);
          R_vax_final_counts(a2,d2) = R_vax_postvax2_counts(a2,d2) ;
          R_vax2_final_counts(a2,d2) = R_vax2_postvax2_counts(a2,d2) ;

          M_final_counts(a2,d2) = M_unvax_final_counts(a2,d2) + M_vax_final_counts(a2,d2) + M_vax2_final_counts(a2,d2);
          S_final_counts(a2,d2) = S_unvax_final_counts(a2,d2) + S_vax_final_counts(a2,d2) + S_vax2_final_counts(a2,d2);
          I_final_counts(a2,d2) = I_unvax_final_counts(a2,d2) + I_vax_final_counts(a2,d2) + I_vax2_final_counts(a2,d2);
          I_final_counts_INCIDENCE(a2,d2) = I_unvax_final_counts_INCIDENCE(a2,d2) + I_vax_final_counts_INCIDENCE(a2,d2) + I_vax2_final_counts_INCIDENCE(a2,d2);
          R_final_counts(a2,d2) = R_unvax_final_counts(a2,d2) + R_vax_final_counts(a2,d2) + R_vax2_final_counts(a2,d2);

          I_final_district_counts(i,d2) += I_final_counts(a2,d2);

        } // stop a2 loop0
        for (int a3=1; a3 < 16; a3++){
         I_final_district_counts_0_4(i,d2) += I_final_counts_INCIDENCE(a3,d2);
        }
        I_final_district_counts_5_9(i,d2) += I_final_counts_INCIDENCE(16,d2);
        I_final_district_counts_10_14(i,d2) += I_final_counts_INCIDENCE(17,d2);
        I_final_district_counts_15_24(i,d2) += I_final_counts_INCIDENCE(18,d2);
        I_final_district_counts_25_34(i,d2) += I_final_counts_INCIDENCE(19,d2);
        I_final_district_counts_35_44(i,d2) += I_final_counts_INCIDENCE(20,d2);
        I_final_district_counts_45_54(i,d2) += I_final_counts_INCIDENCE(21,d2);
        I_final_district_counts_55_64(i,d2) += I_final_counts_INCIDENCE(22,d2);
        I_final_district_counts_65_plus(i,d2) += I_final_counts_INCIDENCE(23,d2);
        
        if(I_final_district_counts_0_4(i,d2) < 0.001){
          I_final_district_counts_0_4(i,d2) = 0.001;
        }
        if(I_final_district_counts_5_9(i,d2) < 0.001){
          I_final_district_counts_5_9(i,d2) = 0.001;
        }
        if(I_final_district_counts_10_14(i,d2) < 0.001){
          I_final_district_counts_10_14(i,d2) = 0.001;
        }
        if(I_final_district_counts_15_24(i,d2) < 0.001){
          I_final_district_counts_15_24(i,d2) = 0.001;
        }
        if(I_final_district_counts_25_34(i,d2) < 0.001){
          I_final_district_counts_25_34(i,d2) = 0.001;
        }
        if(I_final_district_counts_35_44(i,d2) < 0.001){
          I_final_district_counts_35_44(i,d2) = 0.001;
        }
        if(I_final_district_counts_45_54(i,d2) < 0.001){
          I_final_district_counts_45_54(i,d2) = 0.001;
        }
        if(I_final_district_counts_55_64(i,d2) < 0.001){
          I_final_district_counts_55_64(i,d2) = 0.001;
        }
        if(I_final_district_counts_65_plus(i,d2) < 0.001){
          I_final_district_counts_65_plus(i,d2) = 0.001;
        }

      } // stop d2 loop
    } // stop i loop


  NumericVector v1{0};
  NumericVector v2{1, 11, 19, 27, 33};
  NumericVector v3{2, 12, 20, 28, 34, 39, 43, 47, 51, 55, 58, 61};
  NumericVector v4{3, 13, 21};
  NumericVector v5{4};
  NumericVector v6{5, 14, 22};
  NumericVector v7{6};
  NumericVector v8{7, 15, 23, 29, 35, 40, 44, 48, 52, 56, 59, 62, 64, 66, 68, 70, 72};
  NumericVector v9{8, 16, 24, 30, 36, 41, 45, 49, 53};
  NumericVector v10{9, 17, 25, 31, 37, 42, 46, 50, 54, 57, 60, 63, 65, 67, 69, 71, 73, 74, 75, 76, 77, 78};
  NumericVector v11{10, 18, 26, 32, 38};
  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  for(int q8 = 1749; q8 < 2120 ; q8 ++){
    
    // Addis Abeba 
    for(int VALUES = 0; VALUES < 79 ; VALUES ++){
      
      if(bootstrap_lookup(q8,VALUES) == 1){
        
        //for(VALUES : v1){
        
        if( std::find(v1.begin(), v1.end(), VALUES) != v1.end()){
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_final_district_counts_0_4(q8,VALUES) * rho1, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_final_district_counts_5_9(q8,VALUES) * rho1, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_final_district_counts_10_14(q8,VALUES) * rho1, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_final_district_counts_15_24(q8,VALUES) * rho1, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_final_district_counts_25_34(q8,VALUES) * rho1, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_final_district_counts_35_44(q8,VALUES) * rho1, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_final_district_counts_45_54(q8,VALUES) * rho1, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_final_district_counts_55_64(q8,VALUES) * rho1, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_final_district_counts_65_plus(q8,VALUES) * rho1, 5) * autocorrelation(VALUES,8));
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if( std::find(v2.begin(), v2.end(), VALUES) != v2.end()){ // Afar
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_final_district_counts_0_4(q8,VALUES) * rho2, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_final_district_counts_5_9(q8,VALUES) * rho2, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_final_district_counts_10_14(q8,VALUES) * rho2, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_final_district_counts_15_24(q8,VALUES) * rho2, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_final_district_counts_25_34(q8,VALUES) * rho2, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_final_district_counts_35_44(q8,VALUES) * rho2, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_final_district_counts_45_54(q8,VALUES) * rho2, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_final_district_counts_55_64(q8,VALUES) * rho2, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_final_district_counts_65_plus(q8,VALUES) * rho2, 5) * autocorrelation(VALUES,8));
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if( std::find(v3.begin(), v3.end(), VALUES) != v3.end()){// Amhara
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_final_district_counts_0_4(q8,VALUES) * rho3, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_final_district_counts_5_9(q8,VALUES) * rho3, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_final_district_counts_10_14(q8,VALUES) * rho3, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_final_district_counts_15_24(q8,VALUES) * rho3, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_final_district_counts_25_34(q8,VALUES) * rho3, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_final_district_counts_35_44(q8,VALUES) * rho3, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_final_district_counts_45_54(q8,VALUES) * rho3, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_final_district_counts_55_64(q8,VALUES) * rho3, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_final_district_counts_65_plus(q8,VALUES) * rho3, 5) * autocorrelation(VALUES,8));
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if( std::find(v4.begin(), v4.end(), VALUES) != v4.end()){ // Benshangul-Gumaz
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_final_district_counts_0_4(q8,VALUES) * rho4, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_final_district_counts_5_9(q8,VALUES) * rho4, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_final_district_counts_10_14(q8,VALUES) * rho4, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_final_district_counts_15_24(q8,VALUES) * rho4, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_final_district_counts_25_34(q8,VALUES) * rho4, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_final_district_counts_35_44(q8,VALUES) * rho4, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_final_district_counts_45_54(q8,VALUES) * rho4, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_final_district_counts_55_64(q8,VALUES) * rho4, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_final_district_counts_65_plus(q8,VALUES) * rho4, 5) * autocorrelation(VALUES,8));
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if( std::find(v5.begin(), v5.end(), VALUES) != v5.end()){ /// Dire Dawa
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_final_district_counts_0_4(q8,VALUES) * rho5, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_final_district_counts_5_9(q8,VALUES) * rho5, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_final_district_counts_10_14(q8,VALUES) * rho5, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_final_district_counts_15_24(q8,VALUES) * rho5, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_final_district_counts_25_34(q8,VALUES) * rho5, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_final_district_counts_35_44(q8,VALUES) * rho5, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_final_district_counts_45_54(q8,VALUES) * rho5, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_final_district_counts_55_64(q8,VALUES) * rho5, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_final_district_counts_65_plus(q8,VALUES) * rho5, 5) * autocorrelation(VALUES,8));
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if( std::find(v6.begin(), v6.end(), VALUES) != v6.end()){ // Gambela Peoples
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_final_district_counts_0_4(q8,VALUES) * rho6, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_final_district_counts_5_9(q8,VALUES) * rho6, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_final_district_counts_10_14(q8,VALUES) * rho6, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_final_district_counts_15_24(q8,VALUES) * rho6, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_final_district_counts_25_34(q8,VALUES) * rho6, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_final_district_counts_35_44(q8,VALUES) * rho6, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_final_district_counts_45_54(q8,VALUES) * rho6, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_final_district_counts_55_64(q8,VALUES) * rho6, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_final_district_counts_65_plus(q8,VALUES) * rho6, 5) * autocorrelation(VALUES,8));
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if( std::find(v7.begin(), v7.end(), VALUES) != v7.end()){ // Harari People
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_final_district_counts_0_4(q8,VALUES) * rho7, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_final_district_counts_5_9(q8,VALUES) * rho7, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_final_district_counts_10_14(q8,VALUES) * rho7, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_final_district_counts_15_24(q8,VALUES) * rho7, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_final_district_counts_25_34(q8,VALUES) * rho7, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_final_district_counts_35_44(q8,VALUES) * rho7, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_final_district_counts_45_54(q8,VALUES) * rho7, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_final_district_counts_55_64(q8,VALUES) * rho7, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_final_district_counts_65_plus(q8,VALUES) * rho7, 5) * autocorrelation(VALUES,8));
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if( std::find(v8.begin(), v8.end(), VALUES) != v8.end()){ // Oromia
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_final_district_counts_0_4(q8,VALUES) * rho8, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_final_district_counts_5_9(q8,VALUES) * rho8, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_final_district_counts_10_14(q8,VALUES) * rho8, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_final_district_counts_15_24(q8,VALUES) * rho8, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_final_district_counts_25_34(q8,VALUES) * rho8, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_final_district_counts_35_44(q8,VALUES) * rho8, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_final_district_counts_45_54(q8,VALUES) * rho8, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_final_district_counts_55_64(q8,VALUES) * rho8, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_final_district_counts_65_plus(q8,VALUES) * rho8, 5) * autocorrelation(VALUES,8));
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if( std::find(v9.begin(), v9.end(), VALUES) != v9.end()){ // Somali
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_final_district_counts_0_4(q8,VALUES) * rho9, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_final_district_counts_5_9(q8,VALUES) * rho9, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_final_district_counts_10_14(q8,VALUES) * rho9, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_final_district_counts_15_24(q8,VALUES) * rho9, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_final_district_counts_25_34(q8,VALUES) * rho9, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_final_district_counts_35_44(q8,VALUES) * rho9, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_final_district_counts_45_54(q8,VALUES) * rho9, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_final_district_counts_55_64(q8,VALUES) * rho9, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_final_district_counts_65_plus(q8,VALUES) * rho9, 5) * autocorrelation(VALUES,8));
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if( std::find(v10.begin(), v10.end(), VALUES) != v10.end()){ // SNNP
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_final_district_counts_0_4(q8,VALUES) * rho10, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_final_district_counts_5_9(q8,VALUES) * rho10, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_final_district_counts_10_14(q8,VALUES) * rho10, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_final_district_counts_15_24(q8,VALUES) * rho10, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_final_district_counts_25_34(q8,VALUES) * rho10, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_final_district_counts_35_44(q8,VALUES) * rho10, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_final_district_counts_45_54(q8,VALUES) * rho10, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_final_district_counts_55_64(q8,VALUES) * rho10, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_final_district_counts_65_plus(q8,VALUES) * rho10, 5) * autocorrelation(VALUES,8));
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if( std::find(v11.begin(), v11.end(), VALUES) != v11.end()){ // Tigray
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_final_district_counts_0_4(q8,VALUES) * rho11, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_final_district_counts_5_9(q8,VALUES) * rho11, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_final_district_counts_10_14(q8,VALUES) * rho11, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_final_district_counts_15_24(q8,VALUES) * rho11, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_final_district_counts_25_34(q8,VALUES) * rho11, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_final_district_counts_35_44(q8,VALUES) * rho11, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_final_district_counts_45_54(q8,VALUES) * rho11, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_final_district_counts_55_64(q8,VALUES) * rho11, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_final_district_counts_65_plus(q8,VALUES) * rho11, 5) * autocorrelation(VALUES,8));
        }
        
      }  
      }
      
} // end q8 for week!


  return loglike; 
}





