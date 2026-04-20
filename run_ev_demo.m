%% RUN_EV_DEMO

% Author(s): Carlos Formoso, Ousmane Diagne

clear;
clc;
close all;

params = get_ev_user_inputs();

results = ev_vaccine_calculator(params);

disp('EV THERAPY DESIGN SUMMARY');
disp(results.summary_table);


% Example parameter set for a simple EV vaccine formulation scenario.

%{
params.cell_count = 2.0e8;
params.media_volume_ml = 200;
params.ev_per_cell = 1500;
params.isolation_efficiency = 0.18;
params.drug_loading_efficiency = 0.25;
params.drug_input_ug = 600;
params.vaccine_target_ug_per_dose = 40;
params.response_threshold_ug = 35;
params.mc_samples = 5000;
params.max_loadable_drug_ug = 220;
params.loading_half_efficiency = 0.15;
params.mean_uptake_factor = 1.0;
params.mean_immune_factor = 1.0;
params.uptake_noise = 0.12;
params.immune_noise = 0.10;
params.response_slope = 0.06;
%}