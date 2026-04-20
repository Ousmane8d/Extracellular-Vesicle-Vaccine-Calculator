%% EV_VACCINE_CALCULATOR

function results = ev_vaccine_calculator(params)

% Author(s): Carlos Formoso, Ousmane Diagne
%
% Purpose:
% This function estimates extracellular vesicle (EV) yield, loaded cargo,
% number of available vaccine doses, and an approximate response score.
% The model is designed as a simple MATLAB project for EV-based cancer
% vaccine formulation and includes uncertainty analysis through Monte Carlo
% simulation.
%
% Inputs:
% params - structure containing model inputs. Required fields:
%   cell_count                - number of EV-producing cells
%   media_volume_ml           - culture media volume in mL
%   ev_per_cell               - estimated EVs produced per cell
%   isolation_efficiency      - fraction of EVs recovered after isolation
%   drug_loading_efficiency   - nominal loading efficiency into EVs
%   drug_input_ug             - initial drug or antigen amount in micrograms
%   vaccine_target_ug_per_dose- target micrograms per administered dose
%   response_threshold_ug     - payload value where response begins rising
%   mc_samples                - number of Monte Carlo trials
%
% Optional fields:
%   max_loadable_drug_ug      - saturation limit for loaded cargo
%   loading_half_efficiency   - efficiency value giving half-max loading
%   mean_uptake_factor        - mean biological uptake multiplier
%   mean_immune_factor        - mean immune-response multiplier
%   uptake_noise              - variability in uptake
%   immune_noise              - variability in immune response
%   response_slope            - steepness of the sigmoid response curve
%   response_max              - maximum response score
%
% Outputs:
% results - structure containing:
%   summary_table             - deterministic summary of key outputs
%   monte_carlo               - Monte Carlo distributions for model outputs
%   sensitivity               - loading-efficiency sensitivity analysis
%
% Important variables:
% total_ev_produced           - estimated total EVs released by all cells
% recovered_ev                - EVs recovered after isolation losses
% loaded_drug_ug              - loaded therapeutic cargo in micrograms
% effective_payload_ug        - loaded cargo after uptake/immune modifiers
% available_doses             - number of doses that can be prepared
% predicted_response_score    - nonlinear response estimate from payload
%
% Model notes:
% This version keeps the model simple but adds:
% 1) saturating EV loading
% 2) nonlinear immune response
% 3) modest biological variability in uptake and response

required_fields = { ...
    'cell_count', ...
    'media_volume_ml', ...
    'ev_per_cell', ...
    'isolation_efficiency', ...
    'drug_loading_efficiency', ...
    'drug_input_ug', ...
    'vaccine_target_ug_per_dose', ...
    'response_threshold_ug', ...
    'mc_samples'};

for k = 1:numel(required_fields)
    if ~isfield(params, required_fields{k})
        error('Missing required field: %s', required_fields{k});
    end
end

params = set_default_params(params);

% Compute baseline deterministic outputs from the user-provided parameters.
total_ev_produced = params.cell_count * params.ev_per_cell;
recovered_ev = total_ev_produced * params.isolation_efficiency;
loaded_drug_ug = saturating_loading(params.drug_loading_efficiency, params);
available_doses = loaded_drug_ug / params.vaccine_target_ug_per_dose;
effective_payload_ug = loaded_drug_ug * params.mean_uptake_factor * params.mean_immune_factor;
predicted_response_score = sigmoid_response(effective_payload_ug, params);

summary_table = table( ...
    total_ev_produced, ...
    recovered_ev, ...
    loaded_drug_ug, ...
    effective_payload_ug, ...
    available_doses, ...
    predicted_response_score, ...
    'VariableNames', { ...
    'TotalEVProduced', ...
    'RecoveredEV', ...
    'LoadedDrug_ug', ...
    'EffectivePayload_ug', ...
    'AvailableDoses', ...
    'PredictedResponseScore'});

mc = monte_carlo_ev(params);
sense = sensitivity_scan(params);

results.summary_table = summary_table;
results.monte_carlo = mc;
results.sensitivity = sense;

make_ev_plots(params, mc, sense);
end

function mc = monte_carlo_ev(params)
% MONTE_CARLO_EV
% Runs repeated trials with random variation in isolation efficiency,
% loading efficiency, EV yield, uptake, and immune response.

n = params.mc_samples;

iso = clamp01(params.isolation_efficiency + 0.04 * randn(n, 1));
load_eff = clamp01(params.drug_loading_efficiency + 0.05 * randn(n, 1));
ev_per_cell = max(0, params.ev_per_cell + 0.15 * params.ev_per_cell * randn(n, 1));
uptake_factor = max(0, params.mean_uptake_factor + params.uptake_noise * randn(n, 1));
immune_factor = max(0, params.mean_immune_factor + params.immune_noise * randn(n, 1));

recovered_ev = params.cell_count .* ev_per_cell .* iso;
loaded_drug_ug = saturating_loading(load_eff, params);
available_doses = loaded_drug_ug ./ params.vaccine_target_ug_per_dose;
effective_payload_ug = loaded_drug_ug .* uptake_factor .* immune_factor;
response_prob = sigmoid_response(effective_payload_ug, params);

mc.recovered_ev = recovered_ev;
mc.loaded_drug_ug = loaded_drug_ug;
mc.effective_payload_ug = effective_payload_ug;
mc.available_doses = available_doses;
mc.response_prob = response_prob;
end

function sense = sensitivity_scan(params)
% SENSITIVITY_SCAN
% Evaluates how loaded drug, dose count, and response change as loading
% efficiency varies across a reasonable range.

eff_range = linspace(0.05, 0.50, 60);

loaded_drug = saturating_loading(eff_range, params);
doses = loaded_drug ./ params.vaccine_target_ug_per_dose;
effective_payload = loaded_drug .* params.mean_uptake_factor .* params.mean_immune_factor;
response_score = sigmoid_response(effective_payload, params);

sense.loading_efficiency = eff_range;
sense.loaded_drug_ug = loaded_drug;
sense.available_doses = doses;
sense.response_score = response_score;
end

function make_ev_plots(params, mc, sense)
% MAKE_EV_PLOTS
% Creates summary figures for Monte Carlo output and sensitivity trends.

figure('Name', 'EV Therapy Design', 'Color', 'w');

subplot(2,2,1);
histogram(mc.loaded_drug_ug, 30, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'none');
xlabel('Loaded drug (\mug)');
ylabel('Count');
title('Monte Carlo Loaded Drug');
grid on;

subplot(2,2,2);
histogram(mc.available_doses, 30, 'FaceColor', [0.1 0.7 0.5], 'EdgeColor', 'none');
xlabel('Available doses');
ylabel('Count');
title('Monte Carlo Dose Count');
grid on;

subplot(2,2,3);
plot(sense.loading_efficiency, sense.loaded_drug_ug, 'LineWidth', 2, 'Color', [0.85 0.33 0.1]);
xlabel('Drug loading efficiency');
ylabel('Loaded drug (\mug)');
title('Saturating Loading Curve');
grid on;

subplot(2,2,4);
yyaxis left;
plot(sense.loading_efficiency, sense.available_doses, 'LineWidth', 2, 'Color', [0.49 0.18 0.56]);
ylabel('Available doses');
yyaxis right;
plot(sense.loading_efficiency, sense.response_score, '--', 'LineWidth', 2, 'Color', [0.93 0.69 0.13]);
ylabel('Response score');
xlabel('Drug loading efficiency');
title('Dose and Response Trend');
grid on;
end

function loaded_drug_ug = saturating_loading(load_eff, params)
% SATURATING_LOADING
% Uses a simple saturation curve so loading does not increase linearly
% forever at high efficiencies.

loaded_drug_ug = params.max_loadable_drug_ug .* load_eff ./ ...
    (params.loading_half_efficiency + load_eff + eps);
loaded_drug_ug = min(loaded_drug_ug, params.drug_input_ug);
end

function response = sigmoid_response(effective_payload_ug, params)
% SIGMOID_RESPONSE
% Converts effective payload into a bounded biological response score.

response = params.response_max ./ ...
    (1 + exp(-params.response_slope .* (effective_payload_ug - params.response_threshold_ug)));
end

function params = set_default_params(params)
% SET_DEFAULT_PARAMS
% Adds default optional parameters when they are not provided by the user.

defaults.max_loadable_drug_ug = 220;
defaults.loading_half_efficiency = 0.15;
defaults.mean_uptake_factor = 1.0;
defaults.mean_immune_factor = 1.0;
defaults.uptake_noise = 0.12;
defaults.immune_noise = 0.10;
defaults.response_slope = 0.06;
defaults.response_max = 1.0;

fields = fieldnames(defaults);
for k = 1:numel(fields)
    if ~isfield(params, fields{k})
        params.(fields{k}) = defaults.(fields{k});
    end
end
end

function x = clamp01(x)
% CLAMP01
% Restricts a value or vector to the interval from 0 to 1.

x(x < 0) = 0;
x(x > 1) = 1;
end