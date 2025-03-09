Exploring the utility of subnational case notifications in fitting dynamic measles models in Ethiopia

**For data processing scripts, see the following files:**
- 01_prep_subnational_case_notifications.R
- 02_prep_subnational_births.R
- 03_prep_mortality.R
- 04_prep_migration.R
- 05_prep_subnational_population.R
- 06_prep_contact_matrices.R
- 07_prep_gravity_matrices.R
- 08_prep_coverage.R
- 09_save_loess_cases.R

To run models with single reporting rates, use script BCD_bootstrap_single_reporting.R, which calls the following Rcpp function scripts:
- Rcpp_BCD_rho_fitting_bootstrap_district_single_age.cpp
- Rcpp_rolling_average_estimated_no_eff_subnat_only_down_adjust_regional_reporting_vax_eff_vary2_loess_seasonal_serial_interval_incidence_autocorrelation_nbinom5_single_reporting_bootstrap_district.cpp

To run models with regional reporting rates, use script BCD_bootstrap_regional_reporting.R, which calls the following Rcpp function scripts:
- Rcpp_BCD_rho_fitting_bootstrap_district_regional_age.cpp
- Rcpp_rolling_average_estimated_no_eff_subnat_only_down_adjust_regional_reporting_vax_eff_vary2_loess_seasonal_serial_interval_incidence_autocorrelation_nbinom5_regional_reporting_bootstrap_district.cpp


The remaining files are used for diagnostics for boostrapped samples (diagnostics_bootstrap.R) and making predictions and overall diganostic plots (predict_and_diagnostics_overall.R).

