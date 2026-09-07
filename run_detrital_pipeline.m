function run_detrital_pipeline(catchments_root, outdir_root, opts)
% RUN_DETRITAL_PIPELINE
% Discovers all catchment subfolders under catchments_root, runs
% infer_youngest_pulse_from_ZPb then filter_detrital_thermo for each,
% and writes cross-catchment summary and sensitivity CSVs.
% Public outputs use observation-based screening terminology. Geological
% interpretations remain separate from the computed classifications.
%
% Expected file layout (names fixed per catchment):
%   <catchments_root>/<CatchmentName>/ZrnPb.csv
%   <catchments_root>/<CatchmentName>/HblAr.csv
%   <catchments_root>/<CatchmentName>/ApHeApPb.csv
%   <catchments_root>/<CatchmentName>/ZrnHeZrnPb.csv
%
% -----------------------------------------------------------------------
% OUTPUT FOLDER STRUCTURE
% -----------------------------------------------------------------------
%   <outdir_root>/
%     README.txt
%     pipeline_summary.csv          — target-component window + GMM statistics
%     output_summary.csv            — counts by catchment, chronometer, action
%     filter_code_lookup.csv        — code definitions written once per run
%     <CatchmentName>/
%       youngest_zircon_component/  — GMM plot and component summary
%       filter_output/              — older-bound reference (default)
%         filter_results_full.csv   — complete one-date-per-row review table
%         filter_results_coded.csv  — compact table using numeric lookup IDs
%         model_input_ages.csv      — eligible dates for downstream models
%         excluded_ages.csv         — only dates recommended for exclusion
%         review_flags.csv          — dates carrying a separate review flag
%         output_summary.csv        — catchment-level counts
%       sensitivity/                — only when run_sensitivity=true
%         reference_boundary_comparison.csv
%     reference_boundary_sensitivity_summary.csv — when run_sensitivity=true
%
% Only older-than-reference results are recommended for exclusion. Review
% flags do not remove rows from model_input_ages.csv.
%
% -----------------------------------------------------------------------
% USAGE
% -----------------------------------------------------------------------
%   run_detrital_pipeline("Catchments", "Output")
%   run_detrital_pipeline("Catchments", "Output", P_thresh=0.50, Delta=12)
%
% -----------------------------------------------------------------------
% OPTIONS  (name-value)
% -----------------------------------------------------------------------
%   Kmax            max K to test in BIC selection (default 6)
%   K_override      if >0, use this K for ALL catchments (skip BIC)
%   K_override_map  containers.Map of catchment name -> K override
%                   e.g. containers.Map({"TC","RC"},{3,4})
%   Delta           short crystallization-to-cooling interval threshold
%   P_thresh        screening decision probability (default 0.65)
%   Nmc             Monte Carlo draws per grain for ZPb GMM (default 50)
%   TargetComponentAgeRange  allowed range for the selected GMM component
%                   mean (default [-Inf Inf]; e.g. [70 300] Ma). The GMM
%                   still fits the full distribution. PulseAgeRange is a
%                   deprecated compatibility alias.
%   run_sensitivity false (default): writes the primary older-bound result.
%                   true: also evaluates younger-bound and midpoint choices,
%                   but condenses them into one comparison table rather than
%                   creating two additional result-folder trees.
%
% -----------------------------------------------------------------------
% FILTER KNOBS — quick reference
% -----------------------------------------------------------------------
%   Parameter        Default   Tighter        Effect of tightening
%   P_thresh         0.65      0.50           Fewer borderline grains pass
%   Delta            8 Ma      12-15 Ma       More short intervals identified
%   reference mode   older     younger         More grains older than reference

arguments
    catchments_root  (1,1) string  = "Catchments"
    outdir_root      (1,1) string  = "Output"
    opts.Kmax             (1,1) double  = 6
    opts.K_override       (1,1) double  = 0
    opts.K_override_map                = []
    opts.BoundsMethod     (1,1) string  = "gmm_sigma_window"
    opts.NSigma           (1,1) double  = 1.0        % sigma multiplier for window + model start
    opts.Delta            (1,1) double  = 8
    opts.P_thresh         (1,1) double  = 0.65
    opts.Nmc              (1,1) double  = 50
    opts.TargetComponentAgeRange (1,2) double = [NaN NaN]
    opts.PulseAgeRange    (1,2) double  = [-Inf Inf] % deprecated alias
    opts.run_sensitivity  (1,1) logical = false
end

if ~isfolder(catchments_root)
    error("Catchments root folder not found: %s", catchments_root);
end
target_age_range = resolve_target_age_range( ...
    opts.TargetComponentAgeRange, opts.PulseAgeRange);
if ~isfolder(outdir_root), mkdir(outdir_root); end

% Write README on every run so it stays current with the parameters used.
write_readme(outdir_root, opts);

% ---- Discover catchment subfolders ----
entries = dir(catchments_root);
is_dot  = strcmp({entries.name}, '.') | strcmp({entries.name}, '..');
catchment_names = string({entries([entries.isdir] & ~is_dot).name});

if isempty(catchment_names)
    error("No subfolders found in %s", catchments_root);
end
fprintf("Found %d catchment(s): %s\n", numel(catchment_names), ...
    strjoin(catchment_names, ", "));

% ---- Modes to run ----
if opts.run_sensitivity
    modes = ["older_bound", "younger_bound", "midpoint"];
else
    modes = "older_bound";
end

% ---- Storage ----
summary_rows = cell(numel(catchment_names), 1);
% results_store{i, mi} = one-date-per-row results for catchment i, mode mi
results_store = cell(numel(catchment_names), numel(modes));
output_summary_rows = cell(numel(catchment_names), 1);
code_lookup = table();

% ---- Per-catchment loop ----
for i = 1:numel(catchment_names)
    cname = catchment_names(i);
    cdir  = fullfile(catchments_root, cname);

    fprintf("\n=== Processing catchment: %s ===\n", cname);

    % -- File paths --
    f_zpb = fullfile(cdir, "ZrnPb.csv");
    f_ar  = fullfile(cdir, "HblAr.csv");
    f_ap  = fullfile(cdir, "ApHeApPb.csv");
    f_zrn = fullfile(cdir, "ZrnHeZrnPb.csv");

    required = [f_zpb, f_ar, f_ap, f_zrn];
    labels   = ["ZrnPb.csv","HblAr.csv","ApHeApPb.csv","ZrnHeZrnPb.csv"];
    skip = false;
    for k = 1:numel(required)
        if ~isfile(required(k))
            warning("Missing file for catchment %s: %s — skipping.", cname, labels(k));
            summary_rows{i} = make_error_row(cname);
            skip = true; break
        end
    end
    if skip, continue; end

    % -- Determine K for this catchment --
    K_use = opts.K_override;
    if ~isempty(opts.K_override_map) && isKey(opts.K_override_map, char(cname))
        K_use = opts.K_override_map(char(cname));
        fprintf("  Using manual K override = %d for %s\n", K_use, cname);
    end

    % -- Step 1: infer youngest pulse (once per catchment, shared across modes) --
    try
        pulse_opts = {"Kmax", opts.Kmax, "Nmc", opts.Nmc, ...
                      "BoundsMethod", opts.BoundsMethod, "NSigma", opts.NSigma, ...
                      "TargetComponentAgeRange", target_age_range};
        if K_use > 0
            pulse_opts = [pulse_opts, {"K_override", K_use}]; %#ok<AGROW>
        end
        pulse = infer_youngest_pulse_from_ZPb(f_zpb, ...
            fullfile(outdir_root, cname, "youngest_zircon_component"), pulse_opts{:});
        fprintf("  Target-component age window: %.1f - %.1f Ma  (K=%d, selection=%s)\n", ...
            pulse.Tyoung(1), pulse.Tyoung(2), pulse.K_used, pulse.K_selection_method);
        fprintf("  Candidate model-start age: %.1f Ma  (component mu=%.1f + %.0f*sigma; sigma=%.1f Ma)\n", ...
            pulse.model_start_Ma, ...
            pulse.mu_young, pulse.NSigma_used, pulse.sigma_young);
    catch ME
        warning("Pulse inference failed for %s: %s", cname, ME.message);
        summary_rows{i} = make_error_row(cname);
        continue
    end

    % -- Step 2: write the primary result and calculate optional sensitivity --
    mode_results = cell(1, numel(modes));
    for mi = 1:numel(modes)
        mode = modes(mi);
        is_primary = mi == 1;
        cout = fullfile(outdir_root, cname, "filter_output");
        try
            mode_results{mi} = filter_detrital_thermo( ...
                f_ar, f_ap, f_zrn, cout, pulse.Tyoung, ...
                "Delta",         opts.Delta, ...
                "P_thresh",      opts.P_thresh, ...
                "ReferenceMode", mode, ...
                "WriteOutputs", is_primary, ...
                "WriteCodeLookup", false);

            results_store{i, mi} = mode_results{mi}.filter_results;
            if is_primary
                catchment_summary = addvars( ...
                    mode_results{mi}.output_summary, ...
                    repmat(cname, height(mode_results{mi}.output_summary), 1), ...
                    'Before', 1, 'NewVariableNames', 'Catchment');
                output_summary_rows{i} = catchment_summary;
                if isempty(code_lookup)
                    code_lookup = mode_results{mi}.code_lookup;
                end
            end
        catch ME
            warning("Filtering failed for %s (mode=%s): %s", cname, mode, ME.message);
        end
    end

    if opts.run_sensitivity && all(~cellfun(@isempty, mode_results))
        sensitivity_dir = fullfile(outdir_root, cname, "sensitivity");
        if ~isfolder(sensitivity_dir), mkdir(sensitivity_dir); end
        boundary_comparison = build_reference_boundary_comparison( ...
            mode_results{1}.filter_results, ...
            mode_results{2}.filter_results, ...
            mode_results{3}.filter_results);
        writetable(boundary_comparison, ...
            fullfile(sensitivity_dir, "reference_boundary_comparison.csv"));
    end

    % -- Collect pulse summary row --
    summary_rows{i} = table(cname, pulse.Tyoung(1), pulse.Tyoung(2), ...
        pulse.mu_young, pulse.sigma_young, pulse.weight_young, ...
        pulse.model_start_Ma, pulse.NSigma_used, ...
        string(pulse.bounds_method), pulse.pulse_age_range(1), ...
        pulse.pulse_age_range(2), pulse.K_used, pulse.K_bic, ...
        string(pulse.K_selection_method), pulse.K_override_applied, ...
        pulse.Nages, pulse.Nselected, ...
        'VariableNames', {'Catchment','target_window_lo_Ma','target_window_hi_Ma', ...
        'target_component_mean_Ma','target_component_sigma_Ma','target_component_weight', ...
        'candidate_model_start_Ma','NSigma','BoundsMethod', ...
        'target_component_search_lo_Ma','target_component_search_hi_Ma', ...
        'K_used','K_bic_selected','K_selection_method','K_override_applied', ...
        'N_ZPb_ages','N_ZPb_assigned'});

    fprintf("  Done.\n");
end

% ---- Write pipeline summary ----
valid_rows = summary_rows(~cellfun(@isempty, summary_rows));
if ~isempty(valid_rows)
    summary_tbl = vertcat(valid_rows{:});
    writetable(summary_tbl, fullfile(outdir_root, "pipeline_summary.csv"));
    fprintf("\nPipeline summary written to: %s\n", ...
        fullfile(outdir_root, "pipeline_summary.csv"));
else
    warning("No catchments processed successfully.");
    return
end

valid_output_rows = output_summary_rows(~cellfun(@isempty, output_summary_rows));
if ~isempty(valid_output_rows)
    writetable(vertcat(valid_output_rows{:}), ...
        fullfile(outdir_root, "output_summary.csv"));
end
if ~isempty(code_lookup)
    writetable(code_lookup, fullfile(outdir_root, "filter_code_lookup.csv"));
end

% ---- Build and write sensitivity table ----
if opts.run_sensitivity && numel(modes) > 1
    sens_tbl = build_sensitivity_table(catchment_names, modes, results_store);
    if ~isempty(sens_tbl)
        sens_path = fullfile(outdir_root, "reference_boundary_sensitivity_summary.csv");
        writetable(sens_tbl, sens_path);
        fprintf("Sensitivity table written to: %s\n", sens_path);
        print_sensitivity_summary(sens_tbl);
    end
end

end % main function

% =======================================================================
% REFERENCE-BOUNDARY COMPARISON
% =======================================================================
function C = build_reference_boundary_comparison(older, younger, midpoint)
% Condense all three reference choices into one row per dated analysis.
key_vars = {'GrainID','PairID','AnalysisID','System','Chronometer','PairRole','Age_Ma'};
assert(height(older) == height(younger) && height(older) == height(midpoint), ...
    "Reference-boundary result tables have different row counts.");
for v = key_vars
    a = older.(v{1});
    b = younger.(v{1});
    c = midpoint.(v{1});
    if isnumeric(a)
        same_b = all((a == b) | (isnan(a) & isnan(b)));
        same_c = all((a == c) | (isnan(a) & isnan(c)));
    else
        same_b = isequal(string(a), string(b));
        same_c = isequal(string(a), string(c));
    end
    assert(same_b && same_c, ...
        "Reference-boundary result rows do not align for variable %s.", v{1});
end

C = older(:, {'GrainID','PairID','AnalysisID','System','Mineral', ...
    'Chronometer','PairRole','PairStatus','Age_Ma','Age_1sigma_Ma', ...
    'ReviewRecommended','ReviewCode'});

C.PrimaryReferenceAge_Ma = older.ReferenceAge_Ma;
C.PrimaryP_OlderThanReference = older.P_OlderThanReference;
C.PrimaryReferenceClass = older.ReferenceClass;
C.PrimaryAction = older.Action;
C.PrimaryModelInclude = older.ModelInclude;

C.MidpointReferenceAge_Ma = midpoint.ReferenceAge_Ma;
C.MidpointP_OlderThanReference = midpoint.P_OlderThanReference;
C.MidpointReferenceClass = midpoint.ReferenceClass;
C.MidpointAction = midpoint.Action;
C.MidpointModelInclude = midpoint.ModelInclude;

C.YoungerBoundaryReferenceAge_Ma = younger.ReferenceAge_Ma;
C.YoungerBoundaryP_OlderThanReference = younger.P_OlderThanReference;
C.YoungerBoundaryReferenceClass = younger.ReferenceClass;
C.YoungerBoundaryAction = younger.Action;
C.YoungerBoundaryModelInclude = younger.ModelInclude;
end

% =======================================================================
% SENSITIVITY TABLE BUILDER
% =======================================================================
function tbl = build_sensitivity_table(catchment_names, modes, results_store)
% For each catchment x chronometer, report the actual boundary ages and the
% resulting number of model-input dates. No qualitative sensitivity label
% is assigned; users see the direct numerical effect of each choice.

rows = {};

for i = 1:numel(catchment_names)
    cname = catchment_names(i);

    primary = results_store{i, 1};
    if isempty(primary), continue; end
    chronometers = unique(string(primary.Chronometer), 'stable');
    chronometers(chronometers == "ZrnUPb") = []; % reference context only

    for si = 1:numel(chronometers)
        chrono = chronometers(si);
        n_vals = nan(1, numel(modes));
        boundary_vals = nan(1, numel(modes));
        n_total = nnz(string(primary.Chronometer) == chrono);

        for mi = 1:numel(modes)
            L = results_store{i, mi};
            if isempty(L), continue; end
            mask = string(L.Chronometer) == chrono;
            if any(mask)
                n_vals(mi) = nnz(logical(L.ModelInclude(mask)));
                boundary_vals(mi) = L.ReferenceAge_Ma(find(mask, 1, 'first'));
            end
        end

        n_primary = n_vals(1);      % older edge; primary result
        n_younger = n_vals(2);      % younger edge
        n_midpoint = n_vals(3);     % midpoint
        primary_minus_younger = n_primary - n_younger;
        if n_total > 0
            primary_minus_younger_pct_total = 100 * primary_minus_younger / n_total;
        else
            primary_minus_younger_pct_total = NaN;
        end

        rows{end+1} = table( ...
            cname, chrono, n_total, ...
            boundary_vals(1), n_primary, ...
            boundary_vals(3), n_midpoint, ...
            boundary_vals(2), n_younger, ...
            primary_minus_younger, primary_minus_younger_pct_total, ...
            'VariableNames', { ...
                'Catchment', 'Chronometer', 'TotalDatedAnalyses', ...
                'PrimaryOlderBoundary_Ma', 'ModelInputs_PrimaryOlderBoundary', ...
                'MidpointBoundary_Ma', 'ModelInputs_Midpoint', ...
                'YoungerBoundary_Ma', 'ModelInputs_YoungerBoundary', ...
                'PrimaryMinusYounger_Count', 'PrimaryMinusYounger_PercentOfTotal'}); %#ok<AGROW>
    end
end

if isempty(rows)
    tbl = table();
    return
end

tbl = sortrows(vertcat(rows{:}), {'Catchment','Chronometer'});

end

% =======================================================================
% CONSOLE SENSITIVITY SUMMARY
% =======================================================================
function print_sensitivity_summary(tbl)

if isempty(tbl), return; end

fprintf("\n--- Reference-boundary comparison (model-input counts) ---\n");
fprintf("  %-10s  %-12s  %7s  %9s  %9s  %9s  %9s\n", ...
    "Catchment", "Chronometer", "N_total", "Primary", "Midpoint", "Younger", "Difference");
fprintf("  %s\n", repmat('-', 1, 78));

for r = 1:height(tbl)
    fprintf("  %-10s  %-12s  %7d  %9d  %9d  %9d  %+9d\n", ...
        tbl.Catchment(r), tbl.Chronometer(r), tbl.TotalDatedAnalyses(r), ...
        tbl.ModelInputs_PrimaryOlderBoundary(r), tbl.ModelInputs_Midpoint(r), ...
        tbl.ModelInputs_YoungerBoundary(r), tbl.PrimaryMinusYounger_Count(r));
end
fprintf("\n");

end

% =======================================================================
% HELPERS
% =======================================================================
function row = make_error_row(cname)
row = table(cname, NaN, NaN, NaN, NaN, NaN, NaN, NaN, "failed", ...
    NaN, NaN, NaN, NaN, "failed", false, NaN, NaN, ...
    'VariableNames', {'Catchment','target_window_lo_Ma','target_window_hi_Ma', ...
    'target_component_mean_Ma','target_component_sigma_Ma','target_component_weight', ...
    'candidate_model_start_Ma','NSigma','BoundsMethod', ...
    'target_component_search_lo_Ma','target_component_search_hi_Ma', ...
    'K_used','K_bic_selected','K_selection_method','K_override_applied', ...
    'N_ZPb_ages','N_ZPb_assigned'});
end

% =======================================================================
% README WRITER
% =======================================================================
function write_readme(outdir_root, opts)
% Write a run-specific guide using mechanism-neutral terminology.
target_age_range = resolve_target_age_range( ...
    opts.TargetComponentAgeRange, opts.PulseAgeRange);

fid = fopen(fullfile(outdir_root, "README.txt"), "w");
fprintf(fid, "DetritalChronFilter — Detrital Thermochronology Screening Output\n");
fprintf(fid, "================================================================\n\n");

fprintf(fid, "Generated by run_detrital_pipeline.m\n");
fprintf(fid, "Date: %s\n", datestr(now, "yyyy-mm-dd HH:MM:SS")); %#ok<TNOW1,DATST>
fprintf(fid, "Terminology version: neutral-v5-named-full-and-coded\n\n");

fprintf(fid, "INTERPRETATION POLICY\n");
fprintf(fid, "---------------------\n");
fprintf(fid, "Only a cooling age that meets the older-than-reference probability\n");
fprintf(fid, "threshold receives an exclusion recommendation. Paired-age order and\n");
fprintf(fid, "short-interval patterns are review flags only. They may be informative,\n");
fprintf(fid, "but do not independently establish a thermal mechanism or bad analysis.\n\n");

fprintf(fid, "WHICH FILES TO USE\n");
fprintf(fid, "------------------\n");
fprintf(fid, "Each catchment has one filter_output folder with six tables:\n");
fprintf(fid, "  filter_results_full.csv\n");
fprintf(fid, "    All dated analyses, one date per row. Paired dates share GrainID\n");
fprintf(fid, "    and PairID. This version includes full explanations for review.\n");
fprintf(fid, "  filter_results_coded.csv\n");
fprintf(fid, "    Compact publication version. ReferenceResultID and ReviewFlagID\n");
fprintf(fid, "    map to filter_code_lookup.csv; ReviewFlagID 0 means no flag.\n");
fprintf(fid, "  model_input_ages.csv\n");
fprintf(fid, "    Dates eligible for downstream modeling, including review-flagged\n");
fprintf(fid, "    dates that were not recommended for exclusion.\n");
fprintf(fid, "  excluded_ages.csv\n");
fprintf(fid, "    Only dates classified older_than_reference. The filename refers\n");
fprintf(fid, "    to individual dates because paired dates can have different results.\n");
fprintf(fid, "  review_flags.csv\n");
fprintf(fid, "    All dates carrying a review flag. Related paired ages are repeated\n");
fprintf(fid, "    beside them. A flag never causes exclusion, although a flagged date\n");
fprintf(fid, "    can be excluded independently by the older-than-reference rule.\n");
fprintf(fid, "  output_summary.csv\n");
fprintf(fid, "    Counts by chronometer and action.\n\n");
fprintf(fid, "Run-level tables are written once at the output root:\n");
fprintf(fid, "  pipeline_summary.csv, output_summary.csv, filter_code_lookup.csv\n\n");
if opts.run_sensitivity
    fprintf(fid, "Optional reference-boundary sensitivity:\n");
    fprintf(fid, "  <CatchmentName>/sensitivity/reference_boundary_comparison.csv\n");
    fprintf(fid, "  reference_boundary_sensitivity_summary.csv\n");
    fprintf(fid, "    Actual boundary ages and model-input counts for the primary older\n");
    fprintf(fid, "    edge, midpoint, and younger edge are placed side by side;\n");
    fprintf(fid, "    no duplicate per-mode result folders are created.\n\n");
end
fprintf(fid, "Target-component QA:\n");
fprintf(fid, "  <CatchmentName>/youngest_zircon_component/\n");
fprintf(fid, "    BIC curve and GMM fit to ZrnPb ages. Inspect to verify K\n");
fprintf(fid, "    selection and target-component interpretation before use.\n\n");

fprintf(fid, "REFERENCE AND REVIEW CODES\n");
fprintf(fid, "--------------------------\n");
fprintf(fid, "  RT   eligible_after_reference_screen; retain unless review flag applies\n");
fprintf(fid, "  OR1  older_than_reference; probability >= 0.90\n");
fprintf(fid, "  OR2  older_than_reference; decision threshold <= probability < 0.90\n");
fprintf(fid, "       OR1 and OR2 are the only codes that recommend exclusion.\n");
fprintf(fid, "  SC1  review: short crystallization-to-cooling interval; probability >= 0.90\n");
fprintf(fid, "  SC2  review: short interval; moderate probability\n");
fprintf(fid, "  AOI  review: cooling age is older than paired U-Pb\n");
fprintf(fid, "       age beyond combined 2-sigma uncertainty\n");
fprintf(fid, "  AOU  review: nominal age order overlaps within 2-sigma\n");
fprintf(fid, "  II   review: required uncertainty is absent/invalid\n");
fprintf(fid, "  NA   paired-age assessment is not applicable\n");
fprintf(fid, "  RC   zircon U-Pb target-component reference context; not screened as a model-input age\n\n");
fprintf(fid, "The CodeID-to-code-to-definition mapping is also written as\n");
fprintf(fid, "filter_code_lookup.csv. Numeric IDs are nominal labels only.\n");
fprintf(fid, "Full values are reported in P_older_than_reference, P_short_interval,\n");
fprintf(fid, "crystallization_to_cooling_interval_Ma, and interval_1sigma_Ma.\n\n");

fprintf(fid, "SCREENING PARAMETERS USED IN THIS RUN\n");
fprintf(fid, "-------------------------------------\n");
fprintf(fid, "  Delta       = %.0f Ma   (short crystallization-to-cooling interval threshold)\n", opts.Delta);
fprintf(fid, "  P_thresh    = %.2f     (screening decision probability)\n", opts.P_thresh);
fprintf(fid, "  BoundsMethod= %s\n", opts.BoundsMethod);
fprintf(fid, "  NSigma      = %.1f     (component-sigma window multiplier)\n", opts.NSigma);
fprintf(fid, "  Kmax        = %d        (max GMM components tested by BIC)\n", opts.Kmax);
fprintf(fid, "  Nmc         = %d        (Monte Carlo draws per grain for ZPb GMM)\n", opts.Nmc);
fprintf(fid, "  Component search = %.1f to %.1f Ma (allowed selected-component mean)\n", ...
    target_age_range(1), target_age_range(2));
fprintf(fid, "                   Ages outside this range remain in the GMM fit but\n");
fprintf(fid, "                   their components are ineligible for target selection.\n");
if opts.run_sensitivity
    fprintf(fid, "  Sensitivity = true      (one condensed three-boundary comparison)\n\n");
else
    fprintf(fid, "  Sensitivity = false     (primary older-bound reference only; default)\n\n");
end
fprintf(fid, "The primary result uses the older edge of the target-component window.\n");
fprintf(fid, "Boundary labels state the numerical choice and do not imply a preferred\n");
fprintf(fid, "geological interpretation.\n\n");

fclose(fid);
fprintf("  README written to: %s\n", fullfile(outdir_root, "README.txt"));
end

function range = resolve_target_age_range(primary, legacy)
% Prefer the neutral public name while preserving the former API.
if all(isnan(primary))
    range = legacy;
elseif any(isnan(primary))
    error("TargetComponentAgeRange must contain two numeric bounds.");
else
    if ~isequal(legacy, [-Inf Inf]) && ~isequal(primary, legacy)
        error("Specify either TargetComponentAgeRange or PulseAgeRange, not conflicting values for both.");
    end
    range = primary;
end
assert(range(1) < range(2), ...
    "TargetComponentAgeRange must be an increasing [younger older] interval.");
end
