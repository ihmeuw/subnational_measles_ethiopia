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

double logit(double x) {
  double result; 
    result = log( x / (1.0 - x) );
  return result;
}

double invlogit(double x) {
  double result;
    result = 1.0 / (1.0 + exp (-1.0 * x));
  return result;
}

// [[Rcpp::export]]
double dnbinom_nonint(double x, double mu, double size){
  
  double logresult = 0;
  double p =0;
  
  p = size / (size + mu);
  logresult = lgamma(x+size) - lgamma(size) - lgamma(x+1) + size*log(p)+x*log(1-p);
  return(logresult);
}




// [[Rcpp::export]]
NumericVector stl_sort(NumericVector x) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
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
                        double logit_rho2, 
                        double logit_rho3, 
                        double logit_rho4, 
                        double logit_rho5, 
                        double logit_rho6, 
                        double logit_rho7, 
                        double logit_rho8, 
                        double logit_rho9, 
                        double logit_rho10, 
                        double logit_rho11, 
                        double logit_rho_a, 
                        double logit_vax_efficacy) {
  //return 1/sqrt(2*M_PI)*exp(-pow(beta0,2)/2);
  return R::dnorm(max_beta, 3, 0.75, true) +  
    R::dnorm(min_beta, 1, 0.75, true) +  
    R::dnorm(logit_rho1, -5.9, 0.5, true) + 
    R::dnorm(logit_rho2, -5.9, 0.5, true) + 
    R::dnorm(logit_rho3, -5.9, 0.5, true) + 
    R::dnorm(logit_rho4, -5.9, 0.5, true) + 
    R::dnorm(logit_rho5, -5.9, 0.5, true) + 
    R::dnorm(logit_rho6, -5.9, 0.5, true) + 
    R::dnorm(logit_rho7, -5.9, 0.5, true) + 
    R::dnorm(logit_rho8, -5.9, 0.5, true) + 
    R::dnorm(logit_rho9, -5.9, 0.5, true) + 
    R::dnorm(logit_rho10, -5.9, 0.5, true) + 
    R::dnorm(logit_rho11, -5.9, 0.5, true) + 
    R::dnorm(logit_rho_a, 0, 0.5, true) + 
    R::dnorm(logit_vax_efficacy, 1.1, 0.05, true);
  
}



// [[Rcpp::export]]
NumericVector Rcpp_sort(NumericVector x, NumericVector y) {
  // Order the elements of x by sorting y
  // First create a vector of indices
  IntegerVector idx = seq_along(x) - 1;
  // Then sort that vector by the values of y
  std::sort(idx.begin(), idx.end(), [&](int i, int j){return y[i] < y[j];});
  // And return x in that order
  return x[idx];
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
double likelihood_reporting_fit(double rho, 
                                
                                NumericMatrix I_annual_by_district_by_week_0_4, 
                                NumericMatrix I_annual_by_district_by_week_5_9, 
                                NumericMatrix I_annual_by_district_by_week_10_14, 
                                NumericMatrix I_annual_by_district_by_week_15_24, 
                                NumericMatrix I_annual_by_district_by_week_25_34, 
                                NumericMatrix I_annual_by_district_by_week_35_44, 
                                NumericMatrix I_annual_by_district_by_week_45_54, 
                                NumericMatrix I_annual_by_district_by_week_55_64, 
                                NumericMatrix I_annual_by_district_by_week_65_plus, 
                                           
                                           NumericMatrix cases_age_bin0_4_matrix, 
                                           NumericMatrix cases_age_bin5_9_matrix, 
                                           NumericMatrix cases_age_bin10_14_matrix, 
                                           NumericMatrix cases_age_bin15_24_matrix, 
                                           NumericMatrix cases_age_bin25_34_matrix, 
                                           NumericMatrix cases_age_bin35_44_matrix, 
                                           NumericMatrix cases_age_bin45_54_matrix, 
                                           NumericMatrix cases_age_bin55_64_matrix, 
                                           NumericMatrix cases_age_bin65_plus_matrix, 
                                           
                                           
                                           NumericMatrix autocorrelation,
                                           NumericMatrix bootstrap_lookup) {

  double loglike= 0;
  
  
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
    // Addis Abeba 
    for(int VALUES = 0; VALUES < 79 ; VALUES ++){
      
      if(bootstrap_lookup(q8,VALUES) == 1){
        
        //for(VALUES : v1){
        
        if( std::find(v1.begin(), v1.end(), VALUES) != v1.end()){
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_annual_by_district_by_week_0_4(q8,VALUES) * rho, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_annual_by_district_by_week_5_9(q8,VALUES) * rho, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_annual_by_district_by_week_10_14(q8,VALUES) * rho, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_annual_by_district_by_week_15_24(q8,VALUES) * rho, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_annual_by_district_by_week_25_34(q8,VALUES) * rho, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_annual_by_district_by_week_35_44(q8,VALUES) * rho, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_annual_by_district_by_week_45_54(q8,VALUES) * rho, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_annual_by_district_by_week_55_64(q8,VALUES) * rho, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_annual_by_district_by_week_65_plus(q8,VALUES) * rho, 5) * autocorrelation(VALUES,8));
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if( std::find(v2.begin(), v2.end(), VALUES) != v2.end()){ // Afar
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_annual_by_district_by_week_0_4(q8,VALUES) * rho, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_annual_by_district_by_week_5_9(q8,VALUES) * rho, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_annual_by_district_by_week_10_14(q8,VALUES) * rho, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_annual_by_district_by_week_15_24(q8,VALUES) * rho, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_annual_by_district_by_week_25_34(q8,VALUES) * rho, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_annual_by_district_by_week_35_44(q8,VALUES) * rho, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_annual_by_district_by_week_45_54(q8,VALUES) * rho, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_annual_by_district_by_week_55_64(q8,VALUES) * rho, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_annual_by_district_by_week_65_plus(q8,VALUES) * rho, 5) * autocorrelation(VALUES,8));
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if( std::find(v3.begin(), v3.end(), VALUES) != v3.end()){// Amhara
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_annual_by_district_by_week_0_4(q8,VALUES) * rho, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_annual_by_district_by_week_5_9(q8,VALUES) * rho, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_annual_by_district_by_week_10_14(q8,VALUES) * rho, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_annual_by_district_by_week_15_24(q8,VALUES) * rho, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_annual_by_district_by_week_25_34(q8,VALUES) * rho, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_annual_by_district_by_week_35_44(q8,VALUES) * rho, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_annual_by_district_by_week_45_54(q8,VALUES) * rho, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_annual_by_district_by_week_55_64(q8,VALUES) * rho, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_annual_by_district_by_week_65_plus(q8,VALUES) * rho, 5) * autocorrelation(VALUES,8));
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if( std::find(v4.begin(), v4.end(), VALUES) != v4.end()){ // Benshangul-Gumaz
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_annual_by_district_by_week_0_4(q8,VALUES) * rho, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_annual_by_district_by_week_5_9(q8,VALUES) * rho, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_annual_by_district_by_week_10_14(q8,VALUES) * rho, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_annual_by_district_by_week_15_24(q8,VALUES) * rho, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_annual_by_district_by_week_25_34(q8,VALUES) * rho, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_annual_by_district_by_week_35_44(q8,VALUES) * rho, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_annual_by_district_by_week_45_54(q8,VALUES) * rho, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_annual_by_district_by_week_55_64(q8,VALUES) * rho, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_annual_by_district_by_week_65_plus(q8,VALUES) * rho, 5) * autocorrelation(VALUES,8));
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if( std::find(v5.begin(), v5.end(), VALUES) != v5.end()){ /// Dire Dawa
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_annual_by_district_by_week_0_4(q8,VALUES) * rho, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_annual_by_district_by_week_5_9(q8,VALUES) * rho, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_annual_by_district_by_week_10_14(q8,VALUES) * rho, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_annual_by_district_by_week_15_24(q8,VALUES) * rho, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_annual_by_district_by_week_25_34(q8,VALUES) * rho, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_annual_by_district_by_week_35_44(q8,VALUES) * rho, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_annual_by_district_by_week_45_54(q8,VALUES) * rho, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_annual_by_district_by_week_55_64(q8,VALUES) * rho, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_annual_by_district_by_week_65_plus(q8,VALUES) * rho, 5) * autocorrelation(VALUES,8));
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if( std::find(v6.begin(), v6.end(), VALUES) != v6.end()){ // Gambela Peoples
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_annual_by_district_by_week_0_4(q8,VALUES) * rho, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_annual_by_district_by_week_5_9(q8,VALUES) * rho, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_annual_by_district_by_week_10_14(q8,VALUES) * rho, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_annual_by_district_by_week_15_24(q8,VALUES) * rho, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_annual_by_district_by_week_25_34(q8,VALUES) * rho, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_annual_by_district_by_week_35_44(q8,VALUES) * rho, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_annual_by_district_by_week_45_54(q8,VALUES) * rho, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_annual_by_district_by_week_55_64(q8,VALUES) * rho, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_annual_by_district_by_week_65_plus(q8,VALUES) * rho, 5) * autocorrelation(VALUES,8));
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if( std::find(v7.begin(), v7.end(), VALUES) != v7.end()){ // Harari People
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_annual_by_district_by_week_0_4(q8,VALUES) * rho, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_annual_by_district_by_week_5_9(q8,VALUES) * rho, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_annual_by_district_by_week_10_14(q8,VALUES) * rho, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_annual_by_district_by_week_15_24(q8,VALUES) * rho, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_annual_by_district_by_week_25_34(q8,VALUES) * rho, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_annual_by_district_by_week_35_44(q8,VALUES) * rho, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_annual_by_district_by_week_45_54(q8,VALUES) * rho, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_annual_by_district_by_week_55_64(q8,VALUES) * rho, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_annual_by_district_by_week_65_plus(q8,VALUES) * rho, 5) * autocorrelation(VALUES,8));
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if( std::find(v8.begin(), v8.end(), VALUES) != v8.end()){ // Oromia
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_annual_by_district_by_week_0_4(q8,VALUES) * rho, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_annual_by_district_by_week_5_9(q8,VALUES) * rho, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_annual_by_district_by_week_10_14(q8,VALUES) * rho, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_annual_by_district_by_week_15_24(q8,VALUES) * rho, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_annual_by_district_by_week_25_34(q8,VALUES) * rho, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_annual_by_district_by_week_35_44(q8,VALUES) * rho, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_annual_by_district_by_week_45_54(q8,VALUES) * rho, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_annual_by_district_by_week_55_64(q8,VALUES) * rho, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_annual_by_district_by_week_65_plus(q8,VALUES) * rho, 5) * autocorrelation(VALUES,8));
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if( std::find(v9.begin(), v9.end(), VALUES) != v9.end()){ // Somali
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_annual_by_district_by_week_0_4(q8,VALUES) * rho, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_annual_by_district_by_week_5_9(q8,VALUES) * rho, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_annual_by_district_by_week_10_14(q8,VALUES) * rho, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_annual_by_district_by_week_15_24(q8,VALUES) * rho, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_annual_by_district_by_week_25_34(q8,VALUES) * rho, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_annual_by_district_by_week_35_44(q8,VALUES) * rho, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_annual_by_district_by_week_45_54(q8,VALUES) * rho, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_annual_by_district_by_week_55_64(q8,VALUES) * rho, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_annual_by_district_by_week_65_plus(q8,VALUES) * rho, 5) * autocorrelation(VALUES,8));
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if( std::find(v10.begin(), v10.end(), VALUES) != v10.end()){ // SNNP
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_annual_by_district_by_week_0_4(q8,VALUES) * rho, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_annual_by_district_by_week_5_9(q8,VALUES) * rho, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_annual_by_district_by_week_10_14(q8,VALUES) * rho, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_annual_by_district_by_week_15_24(q8,VALUES) * rho, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_annual_by_district_by_week_25_34(q8,VALUES) * rho, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_annual_by_district_by_week_35_44(q8,VALUES) * rho, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_annual_by_district_by_week_45_54(q8,VALUES) * rho, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_annual_by_district_by_week_55_64(q8,VALUES) * rho, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_annual_by_district_by_week_65_plus(q8,VALUES) * rho, 5) * autocorrelation(VALUES,8));
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        if( std::find(v11.begin(), v11.end(), VALUES) != v11.end()){ // Tigray
          loglike += (dnbinom_nonint(cases_age_bin0_4_matrix(q8,VALUES), I_annual_by_district_by_week_0_4(q8,VALUES) * rho, 5) * autocorrelation(VALUES,0));
          loglike += (dnbinom_nonint(cases_age_bin5_9_matrix(q8,VALUES), I_annual_by_district_by_week_5_9(q8,VALUES) * rho, 5) * autocorrelation(VALUES,1));
          loglike += (dnbinom_nonint(cases_age_bin10_14_matrix(q8,VALUES), I_annual_by_district_by_week_10_14(q8,VALUES) * rho, 5) * autocorrelation(VALUES,2));
          loglike += (dnbinom_nonint(cases_age_bin15_24_matrix(q8,VALUES), I_annual_by_district_by_week_15_24(q8,VALUES) * rho, 5) * autocorrelation(VALUES,3));
          loglike += (dnbinom_nonint(cases_age_bin25_34_matrix(q8,VALUES), I_annual_by_district_by_week_25_34(q8,VALUES) * rho, 5) * autocorrelation(VALUES,4));
          loglike += (dnbinom_nonint(cases_age_bin35_44_matrix(q8,VALUES), I_annual_by_district_by_week_35_44(q8,VALUES) * rho, 5) * autocorrelation(VALUES,5));
          loglike += (dnbinom_nonint(cases_age_bin45_54_matrix(q8,VALUES), I_annual_by_district_by_week_45_54(q8,VALUES) * rho, 5) * autocorrelation(VALUES,6));
          loglike += (dnbinom_nonint(cases_age_bin55_64_matrix(q8,VALUES), I_annual_by_district_by_week_55_64(q8,VALUES) * rho, 5) * autocorrelation(VALUES,7));
          loglike += (dnbinom_nonint(cases_age_bin65_plus_matrix(q8,VALUES), I_annual_by_district_by_week_65_plus(q8,VALUES) * rho, 5) * autocorrelation(VALUES,8));
        }
        
      }  
    }
} // end q8 for week!


  return loglike; 
}

 
 
 


