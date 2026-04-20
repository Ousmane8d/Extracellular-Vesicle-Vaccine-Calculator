%% GET_EV_USER_INPUTS

function params = get_ev_user_inputs()

% Author(s): Carlos Formoso, Ousmane Diagne
%
% Purpose:
% This script demonstrates how to run the EV vaccine calculator with a
% single set of example formulation parameters.
%
% Inputs:
% No external inputs are required. The script defines a parameter
% structure named "params" directly in the file.
%
% Outputs:
% 1) MATLAB figures showing Monte Carlo and sensitivity plots
% 2) A summary table printed in the Command Window
%
% Important variables:
% params.cell_count               - EV-producing cell count
% params.ev_per_cell              - estimated EV yield per cell
% params.isolation_efficiency     - EV recovery fraction
% params.drug_loading_efficiency  - nominal cargo loading efficiency
% params.drug_input_ug            - initial drug or antigen amount
% params.vaccine_target_ug_per_dose - target cargo per vaccine dose
% params.response_threshold_ug    - midpoint of response curve
% params.mc_samples               - number of Monte Carlo simulations


fprintf('\n============================\n');
fprintf('EV THERAPY INPUT INTERFACE\n');
fprintf('============================\n\n');

params.cell_count = input('Enter producing cell count (e.g. 2e8): ');
params.media_volume_ml = input('Enter media volume (mL): ');
params.ev_per_cell = input('Enter EV yield per cell (e.g. 1500): ');
params.isolation_efficiency = input('Enter isolation efficiency (0–1): ');
params.drug_loading_efficiency = input('Enter drug loading efficiency (0–1): ');
params.drug_input_ug = input('Enter total drug/antigen input (ug): ');

params.vaccine_target_ug_per_dose = input('Enter target ug per dose: ');

params.mean_uptake_factor = input('Enter mean uptake factor (default 1.0): ');
params.mean_immune_factor = input('Enter mean immune factor (default 1.0): ');
params.uptake_noise = input('Enter uptake noise (e.g. 0.12): ');
params.immune_noise = input('Enter immune noise (e.g. 0.10): ');

params.response_threshold_ug = input('Enter response threshold (ug): ');
params.mc_samples = input('Enter Monte Carlo sample size (e.g. 5000): ');

fprintf('\nInputs successfully loaded.\n\n');

end


