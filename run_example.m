%% RUN_EXAMPLE
% End-to-end walkthrough of DetritalThermoFilter using synthetic data.
%
% Run this script from the repo root (one level above the example/ folder)
% to verify your MATLAB setup before using your own data.
%
% Expected runtime: < 1 minute on any modern machine.

%% -------------------------------------------------------------------------
% Step 0: Generate synthetic catchment data
% -------------------------------------------------------------------------
% This creates a Catchments_example/ folder with three fictional catchments
% (CatchmentA, CatchmentB, CatchmentC), each containing the four required
% CSV files. Grain ages are drawn from Gaussian mixture distributions with
% known pulse parameters so you can check the pipeline recovers them.

fprintf("=== Generating synthetic data ===\n");
example_dir = fullfile(fileparts(mfilename('fullpath')));
addpath(example_dir);                   % so generate_synthetic_data is found
addpath(fileparts(example_dir));        % so pipeline functions are found

generate_synthetic_data("Catchments_example");

%% -------------------------------------------------------------------------
% Step 1: Run the full pipeline
% -------------------------------------------------------------------------
% This is the same call you would use on your own data — just swap
% "Catchments_example" for your actual Catchments folder path.

fprintf("\n=== Running pipeline ===\n");
run_detrital_pipeline( ...
    "Catchments_example", ...   % input folder
    "Output_example", ...       % output folder
    "Kmax",      5, ...         % test K = 2..5 (3 is plenty for synthetic data)
    "NSigma",    1.0, ...       % pulse window = mu +/- 1*sigma (default)
    "P_thresh",  0.65, ...      % probability threshold for exclusion/flagging (default)
    "Delta",     8);            % magmatic lag window in Ma (default)

%% -------------------------------------------------------------------------
% Step 2: Inspect the sensitivity table
% -------------------------------------------------------------------------
fprintf("\n=== Sensitivity table ===\n");
sens = readtable("Output_example/_pipeline_summary/sensitivity_by_system.csv");
disp(sens(:, {'Catchment','System','N_keep_standard','N_keep_conservative','Sensitivity'}));

%% -------------------------------------------------------------------------
% Step 3: Inspect the pipeline summary (pulse windows + model start times)
% -------------------------------------------------------------------------
fprintf("\n=== Pipeline summary ===\n");
psum = readtable("Output_example/_pipeline_summary/pipeline_summary.csv");
disp(psum(:, {'Catchment','mu_young','sigma_young','model_start_Ma','model_start_err_Ma','Tyoung_lo','Tyoung_hi'}));

%% -------------------------------------------------------------------------
% Step 4: What to check
% -------------------------------------------------------------------------
% After running, open Output_example/ and inspect:
%
%   README.txt
%       - Auto-generated guide to all output files and reason codes.
%         Read this first if anything is unclear.
%
%   _pipeline_summary/pipeline_summary.csv
%       - mu_young should be close to the known pulse means:
%           CatchmentA: ~90 Ma
%           CatchmentB: ~75 Ma
%           CatchmentC: ~100 Ma
%       - Tyoung_lo / Tyoung_hi = mu +/- NSigma*sigma (with NSigma=1.0):
%           CatchmentA: ~85–95 Ma  (mu=90, sigma~5)
%           CatchmentB: ~69–81 Ma  (mu=75, sigma~6)
%           CatchmentC: ~93–107 Ma (mu=100, sigma~7)
%       - model_start_Ma = mu + 1*sigma; use as QTQt oldest constraint
%       - model_start_err_Ma = sigma_young; use as +/- for QTQt bounding box
%       - K_used: for clean 3-population synthetic data, BIC typically selects K=3
%
%   _pipeline_summary/sensitivity_by_system.csv
%       - Hb_ArAr will likely show HIGH sensitivity (expected — hornblende
%         ages straddle the pulse window boundary by design)
%       - Ap_He+Ap_UPb and Zrn_He+Zrn_UPb should show low-moderate sensitivity
%
%   CatchmentA/ZPb_QA/ZrnPb_youngest_pulse_QA.png
%       - BIC panel: should show a clear minimum or plateau at K=3
%       - Age distribution panel: shaded window covers ~85–95 Ma peak;
%         red vertical line marks model_start_Ma (~95 Ma)
%
%   CatchmentA/standard/summary_counts.csv
%       - Ap_He+Ap_UPb: expect ~80 keep_exhumation, ~10 exclude_legacy,
%         ~10 flag_magmatic (the tight 1-5 Myr synthetic lags fall well
%         within the 8 Ma Delta window and will reliably be flagged)
%       - Hb_ArAr: expect ~55-60 keep_exhumation, ~40-45 exclude_legacy
%         (hornblende ages cluster near the pulse window edge by design)
%       - Zrn_He+Zrn_UPb: expect ~65 keep_exhumation, ~13 flag_magmatic,
%         ~8 exclude_legacy
%
%   CatchmentA/standard/all_data_classified.csv
%       - reason_code column: KE = kept, EL1/EL2 = legacy, FM1/FM2 = magmatic
%       - code_ApPb / class_ApPb: independent U-Pb classification for apatites
%       - confidence column (0-1): higher = more confidently classified as KE

fprintf("\nExample run complete. Outputs written to: Output_example/\n");
fprintf("See comments in this script (Step 4) for what to look for in the outputs.\n");
fprintf("Key file for thermal modeling: Output_example/CatchmentA/standard/kept_strict.csv\n");
fprintf("Key file for QTQt setup:       Output_example/_pipeline_summary/pipeline_summary.csv\n");
