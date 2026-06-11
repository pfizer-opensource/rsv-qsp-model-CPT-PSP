Code repository for simulating the RSV viral dynamics model and key results used for the publication Utilizing virtual clinical trials to inform target coverage that drives RSV antiviral efficacy.

Key files:
- rsv_model.m: RSV viral load ODE model equations.
- SolveBalances_rsv.m: ODE simulation function.
- PK_check.m: One-compartment PK analytical helper function for simulating out sisunatovir, zelicapavir pharmacokinetics.
- terminal_slope_ls.m: Helper function to calculate post-peak slopes of viral load.
- treatment_sims.m: Runs zelicapavir treatment simulations for the healthy adult virtual population (merged VPops from sisunatovir and zelicapavir part 1 & 2 challenge studies) and the pediatric virtual population for varying treatment intervention timing.
- treatment_fig_plots.m: Plots out the results from treatment_sims.m. Also, plots out bar plot for the study meta-analysis performed on challenge, RCT, and observational studies.
- bootstrap_success_trials.m: Performs the supplemental, basic power analysis on day 3 VL reduction for varying effect sizes using bootstrapping. Tests that the 95% confidence limit is above a given effect size.
- rsv_qsp_unit_test.m: Simulate and plot the supplementary unit tests for a standard, increased, and impaired immune response. 

MAT files:
- vpop_X.mat: Virtual populations matched to healthy adults (3 challenge studies, for sisunatovir and zelicapavir part 1 & 2) and pediatrics (zelicapavir Ph2).
- study_analysis.mat: Results (e.g., peak viral load,  slope post-peak) for the meta-analysis performed on challenge, RCT, and observational studies.
- fit_EDP938_params.mat: PK parameters for the one-compartment PK model for zelicapavir.
- best_fit_Finhib.mat: Best fitting parameter set for sisunatovir VPop, used to perform unit tests.
