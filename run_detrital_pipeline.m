function run_detrital_pipeline(samples_root, outdir_root, opts)
% RUN_DETRITAL_PIPELINE
% Processes one sample folder directly or discovers multiple sample
% subfolders under samples_root. Runs infer_target_component followed by
% filter_detrital_thermo and writes cross-sample summaries.
% Public outputs use observation-based screening terminology. Geological
% interpretations remain separate from the computed classifications.
%
% Expected files in each sample folder:
%   ReferenceDistribution.csv  complete distribution used for GMM fitting
%   ChronometerData.csv        any number of chronometer systems, long format
%
% -----------------------------------------------------------------------
% OUTPUT FOLDER STRUCTURE
% -----------------------------------------------------------------------
%   <outdir_root>/
%     README.txt
%     pipeline_summary.csv          — target-component window + GMM statistics
%     output_summary.csv            — age ranges and counts by sample and chronometer
%     filter_code_lookup.csv        — code definitions written once per run
%     <SampleName>/
%       target_component/           — GMM plot and component summary
%       filter_output/              — older-bound reference (default)
%         filter_results_full.csv   — complete one-date-per-row review table
%         filter_results_coded.csv  — compact table using numeric lookup IDs
%         model_input_ages.csv      — eligible dates for downstream models
%         excluded_ages.csv         — only dates recommended for exclusion
%         review_flags.csv          — dates carrying a separate review flag
%         output_summary.csv        — chronometer age ranges and counts
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
%   run_detrital_pipeline("Samples", "Output")
%   run_detrital_pipeline("Samples/SampleA", "Output")
%   run_detrital_pipeline("Samples", "Output", P_thresh=0.50, Delta=12)
%
% -----------------------------------------------------------------------
% OPTIONS  (name-value)
% -----------------------------------------------------------------------
%   Kmax            max K to test in BIC selection (default 6)
%   K_override      if >0, use this K for ALL samples (skip BIC)
%   K_override_map  containers.Map of sample name -> K override
%                   e.g. containers.Map({"TC","RC"},{3,4})
%   Delta           short paired-age interval threshold
%   P_thresh        screening decision probability (default 0.65)
%   Nmc             Monte Carlo draws per reference age for GMM (default 50)
%   TargetComponentAgeRange  allowed range for the selected GMM component
%                   mean (default [-Inf Inf]; e.g. [50 200] Ma). The GMM
%                   still fits the full distribution.
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
    samples_root     (1,1) string  = "Samples"
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
    opts.run_sensitivity  (1,1) logical = false
end

if ~isfolder(samples_root)
    error("Samples root folder not found: %s", samples_root);
end
target_age_range = resolve_target_age_range(opts.TargetComponentAgeRange);

% ---- Discover one direct sample or multiple sample subfolders ----
[sample_names, sample_dirs, input_layout] = discover_samples(samples_root);
fprintf("Found %d sample(s) using %s layout: %s\n", numel(sample_names), ...
    input_layout, strjoin(sample_names, ", "));

if ~isfolder(outdir_root), mkdir(outdir_root); end
% Write README on every valid run so it stays current with the parameters used.
write_readme(outdir_root, opts);

% ---- Modes to run ----
if opts.run_sensitivity
    modes = ["older_bound", "younger_bound", "midpoint"];
else
    modes = "older_bound";
end

% ---- Storage ----
summary_rows = cell(numel(sample_names), 1);
% results_store{i, mi} = one-date-per-row results for sample i, mode mi
results_store = cell(numel(sample_names), numel(modes));
output_summary_rows = cell(numel(sample_names), 1);
code_lookup = table();

% ---- Per-sample loop ----
for i = 1:numel(sample_names)
    sample_name = sample_names(i);
    sample_dir = sample_dirs(i);

    fprintf("\n=== Processing sample: %s ===\n", sample_name);

    % -- File paths --
    f_reference = fullfile(sample_dir, "ReferenceDistribution.csv");
    f_chronometers = fullfile(sample_dir, "ChronometerData.csv");

    % -- Determine K for this sample --
    K_use = opts.K_override;
    if ~isempty(opts.K_override_map) && isKey(opts.K_override_map, char(sample_name))
        K_use = opts.K_override_map(char(sample_name));
        fprintf("  Using manual K override = %d for %s\n", K_use, sample_name);
    end

    % -- Step 1: infer target component (once per sample) --
    try
        target_opts = {"Kmax", opts.Kmax, "Nmc", opts.Nmc, ...
                      "BoundsMethod", opts.BoundsMethod, "NSigma", opts.NSigma, ...
                      "TargetComponentAgeRange", target_age_range};
        if K_use > 0
            target_opts = [target_opts, {"K_override", K_use}]; %#ok<AGROW>
        end
        target = infer_target_component(f_reference, ...
            fullfile(outdir_root, sample_name, "target_component"), target_opts{:});
        fprintf("  Target-component age window: %.1f - %.1f Ma  (K=%d, selection=%s)\n", ...
            target.target_window_Ma(1), target.target_window_Ma(2), ...
            target.K_used, target.K_selection_method);
        fprintf("  Candidate model-start age: %.1f Ma  (component mu=%.1f + %.0f*sigma; sigma=%.1f Ma)\n", ...
            target.candidate_model_start_Ma, ...
            target.target_component_mean_Ma, target.NSigma_used, ...
            target.target_component_sigma_Ma);
    catch ME
        warning("Target-component inference failed for %s: %s", sample_name, ME.message);
        summary_rows{i} = make_error_row(sample_name);
        continue
    end

    % -- Step 2: write the primary result and calculate optional sensitivity --
    mode_results = cell(1, numel(modes));
    for mi = 1:numel(modes)
        mode = modes(mi);
        is_primary = mi == 1;
        cout = fullfile(outdir_root, sample_name, "filter_output");
        try
            mode_results{mi} = filter_detrital_thermo( ...
                f_chronometers, cout, target.target_window_Ma, ...
                "Delta",         opts.Delta, ...
                "P_thresh",      opts.P_thresh, ...
                "ReferenceMode", mode, ...
                "WriteOutputs", is_primary, ...
                "WriteCodeLookup", false);

            results_store{i, mi} = mode_results{mi}.filter_results;
            if is_primary
                sample_summary = addvars( ...
                    mode_results{mi}.output_summary, ...
                    repmat(sample_name, height(mode_results{mi}.output_summary), 1), ...
                    'Before', 1, 'NewVariableNames', 'Sample');
                output_summary_rows{i} = sample_summary;
                if isempty(code_lookup)
                    code_lookup = mode_results{mi}.code_lookup;
                end
            end
        catch ME
            warning("Filtering failed for %s (mode=%s): %s", sample_name, mode, ME.message);
        end
    end

    if opts.run_sensitivity && all(~cellfun(@isempty, mode_results))
        sensitivity_dir = fullfile(outdir_root, sample_name, "sensitivity");
        if ~isfolder(sensitivity_dir), mkdir(sensitivity_dir); end
        boundary_comparison = build_reference_boundary_comparison( ...
            mode_results{1}.filter_results, ...
            mode_results{2}.filter_results, ...
            mode_results{3}.filter_results);
        writetable(boundary_comparison, ...
            fullfile(sensitivity_dir, "reference_boundary_comparison.csv"));
    end

    % -- Collect target-component summary row --
    summary_rows{i} = table(sample_name, target.reference_system, ...
        target.target_window_Ma(1), target.target_window_Ma(2), ...
        target.target_component_mean_Ma, target.target_component_sigma_Ma, ...
        target.target_component_weight, target.candidate_model_start_Ma, ...
        target.NSigma_used, string(target.bounds_method), ...
        target.target_component_age_range(1), ...
        target.target_component_age_range(2), target.K_used, target.K_bic, ...
        string(target.K_selection_method), target.K_override_applied, ...
        target.Nages, target.Nselected, ...
        'VariableNames', {'Sample','ReferenceSystem', ...
        'target_window_lo_Ma','target_window_hi_Ma', ...
        'target_component_mean_Ma','target_component_sigma_Ma','target_component_weight', ...
        'candidate_model_start_Ma','NSigma','BoundsMethod', ...
        'target_component_search_lo_Ma','target_component_search_hi_Ma', ...
        'K_used','K_bic_selected','K_selection_method','K_override_applied', ...
        'N_reference_ages','N_reference_ages_assigned'});

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
    warning("No samples processed successfully.");
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
    sens_tbl = build_sensitivity_table(sample_names, modes, results_store);
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
key_vars = {'GrainID','PairID','AnalysisID','Chronometer','PairRole', ...
    'PairStatus','UseForModel','Age_Ma'};
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

C = older(:, {'GrainID','PairID','AnalysisID','Chronometer','PairRole', ...
    'PairStatus','UseForModel','Age_Ma','Age_1sigma_Ma', ...
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
function tbl = build_sensitivity_table(sample_names, modes, results_store)
% For each sample x chronometer, report the actual boundary ages and the
% resulting number of model-input dates. No qualitative sensitivity label
% is assigned; users see the direct numerical effect of each choice.

rows = {};

for i = 1:numel(sample_names)
    sample_name = sample_names(i);

    primary = results_store{i, 1};
    if isempty(primary), continue; end
    chronometers = unique(string(primary.Chronometer(primary.UseForModel)), 'stable');

    for si = 1:numel(chronometers)
        chrono = chronometers(si);
        n_vals = nan(1, numel(modes));
        boundary_vals = nan(1, numel(modes));
        n_total = nnz(string(primary.Chronometer) == chrono & primary.UseForModel);

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
            sample_name, chrono, n_total, ...
            boundary_vals(1), n_primary, ...
            boundary_vals(3), n_midpoint, ...
            boundary_vals(2), n_younger, ...
            primary_minus_younger, primary_minus_younger_pct_total, ...
            'VariableNames', { ...
                'Sample', 'Chronometer', 'TotalDatedAnalyses', ...
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

tbl = sortrows(vertcat(rows{:}), {'Sample','Chronometer'});

end

% =======================================================================
% CONSOLE SENSITIVITY SUMMARY
% =======================================================================
function print_sensitivity_summary(tbl)

if isempty(tbl), return; end

fprintf("\n--- Reference-boundary comparison (model-input counts) ---\n");
fprintf("  %-10s  %-12s  %7s  %9s  %9s  %9s  %9s\n", ...
    "Sample", "Chronometer", "N_total", "Primary", "Midpoint", "Younger", "Difference");
fprintf("  %s\n", repmat('-', 1, 78));

for r = 1:height(tbl)
    fprintf("  %-10s  %-12s  %7d  %9d  %9d  %9d  %+9d\n", ...
        tbl.Sample(r), tbl.Chronometer(r), tbl.TotalDatedAnalyses(r), ...
        tbl.ModelInputs_PrimaryOlderBoundary(r), tbl.ModelInputs_Midpoint(r), ...
        tbl.ModelInputs_YoungerBoundary(r), tbl.PrimaryMinusYounger_Count(r));
end
fprintf("\n");

end

% =======================================================================
% HELPERS
% =======================================================================
function [sample_names, sample_dirs, layout] = discover_samples(samples_root)
% Accept either one sample directly or multiple samples in subfolders.
reference_name = "ReferenceDistribution.csv";
chronometer_name = "ChronometerData.csv";
root_has_reference = isfile(fullfile(samples_root, reference_name));
root_has_chronometers = isfile(fullfile(samples_root, chronometer_name));

entries = dir(samples_root);
is_dot = strcmp({entries.name}, '.') | strcmp({entries.name}, '..');
subfolders = entries([entries.isdir] & ~is_dot);
subfolder_names = string({subfolders.name});
subfolder_dirs = fullfile(samples_root, subfolder_names);
sub_has_reference = false(size(subfolder_dirs));
sub_has_chronometers = false(size(subfolder_dirs));
for i = 1:numel(subfolder_dirs)
    sub_has_reference(i) = isfile(fullfile(subfolder_dirs(i), reference_name));
    sub_has_chronometers(i) = isfile(fullfile(subfolder_dirs(i), chronometer_name));
end
sub_has_any_input = sub_has_reference | sub_has_chronometers;

if root_has_reference || root_has_chronometers
    if ~(root_has_reference && root_has_chronometers)
        error("DetritalChronFilter:IncompleteSample", ...
            "Single-sample folder %s must contain both %s and %s.", ...
            samples_root, reference_name, chronometer_name);
    end
    if any(sub_has_any_input)
        error("DetritalChronFilter:AmbiguousInputLayout", ...
            "Ambiguous input layout in %s: recognized files occur both " + ...
            "directly and in sample subfolders. Use one layout per run.", ...
            samples_root);
    end
    normalized_root = char(samples_root);
    while numel(normalized_root) > 1 && normalized_root(end) == filesep
        normalized_root(end) = [];
    end
    [~, root_name] = fileparts(normalized_root);
    if isempty(root_name)
        root_name = "Sample";
    end
    sample_names = string(root_name);
    sample_dirs = samples_root;
    layout = "single-sample";
    return
end

incomplete = xor(sub_has_reference, sub_has_chronometers);
if any(incomplete)
    bad_name = subfolder_names(find(incomplete, 1));
    error("DetritalChronFilter:IncompleteSample", ...
        "Sample folder %s is incomplete. Each sample must contain both %s and %s.", ...
        bad_name, reference_name, chronometer_name);
end

valid = sub_has_reference & sub_has_chronometers;
sample_names = subfolder_names(valid);
sample_dirs = subfolder_dirs(valid);
if isempty(sample_names)
    error("DetritalChronFilter:NoSampleInputs", ...
        "No valid sample inputs found in %s. Supply both recognized CSVs " + ...
        "directly or place them together in each sample subfolder.", ...
        samples_root);
end
layout = "multi-sample";
end

function row = make_error_row(sample_name)
row = table(sample_name, "failed", NaN, NaN, NaN, NaN, NaN, NaN, NaN, "failed", ...
    NaN, NaN, NaN, NaN, "failed", false, NaN, NaN, ...
    'VariableNames', {'Sample','ReferenceSystem', ...
    'target_window_lo_Ma','target_window_hi_Ma', ...
    'target_component_mean_Ma','target_component_sigma_Ma','target_component_weight', ...
    'candidate_model_start_Ma','NSigma','BoundsMethod', ...
    'target_component_search_lo_Ma','target_component_search_hi_Ma', ...
    'K_used','K_bic_selected','K_selection_method','K_override_applied', ...
    'N_reference_ages','N_reference_ages_assigned'});
end

% =======================================================================
% README WRITER
% =======================================================================
function write_readme(outdir_root, opts)
% Write a run-specific guide using mechanism-neutral terminology.
target_age_range = resolve_target_age_range(opts.TargetComponentAgeRange);

fid = fopen(fullfile(outdir_root, "README.txt"), "w");
fprintf(fid, "DetritalChronFilter — Detrital Thermochronology Screening Output\n");
fprintf(fid, "================================================================\n\n");

fprintf(fid, "Generated by run_detrital_pipeline.m\n");
fprintf(fid, "Date: %s\n", datestr(now, "yyyy-mm-dd HH:MM:SS")); %#ok<TNOW1,DATST>
fprintf(fid, "Terminology version: sample-generic-v1\n\n");

fprintf(fid, "INTERPRETATION POLICY\n");
fprintf(fid, "---------------------\n");
fprintf(fid, "Only a model-candidate age that meets the older-than-reference\n");
fprintf(fid, "threshold receives an exclusion recommendation. Paired-age order and\n");
fprintf(fid, "short-interval patterns are review flags only. They may be informative,\n");
fprintf(fid, "but do not independently establish a thermal mechanism or bad analysis.\n\n");

fprintf(fid, "WHICH FILES TO USE\n");
fprintf(fid, "------------------\n");
fprintf(fid, "Each sample has one filter_output folder with six tables:\n");
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
fprintf(fid, "    One row per chronometer with observed age range, model-input\n");
fprintf(fid, "    age range and median, and excluded/review/context-only counts.\n");
fprintf(fid, "    Range comparisons are descriptive and do not apply a\n");
fprintf(fid, "    closure-temperature ordering rule.\n\n");
fprintf(fid, "Run-level tables are written once at the output root:\n");
fprintf(fid, "  pipeline_summary.csv, output_summary.csv, filter_code_lookup.csv\n\n");
if opts.run_sensitivity
    fprintf(fid, "Optional reference-boundary sensitivity:\n");
    fprintf(fid, "  <SampleName>/sensitivity/reference_boundary_comparison.csv\n");
    fprintf(fid, "  reference_boundary_sensitivity_summary.csv\n");
    fprintf(fid, "    Actual boundary ages and model-input counts for the primary older\n");
    fprintf(fid, "    edge, midpoint, and younger edge are placed side by side;\n");
    fprintf(fid, "    no duplicate per-mode result folders are created.\n\n");
end
fprintf(fid, "Target-component QA:\n");
fprintf(fid, "  <SampleName>/target_component/\n");
fprintf(fid, "    BIC curve and GMM fit to ReferenceDistribution.csv. Inspect K\n");
fprintf(fid, "    selection and target-component interpretation before use.\n\n");

fprintf(fid, "REFERENCE AND REVIEW CODES\n");
fprintf(fid, "--------------------------\n");
fprintf(fid, "  RT   eligible_after_reference_screen; retain unless review flag applies\n");
fprintf(fid, "  OR1  older_than_reference; probability >= 0.90\n");
fprintf(fid, "  OR2  older_than_reference; decision threshold <= probability < 0.90\n");
fprintf(fid, "       OR1 and OR2 are the only codes that recommend exclusion.\n");
fprintf(fid, "  SI1  review: short paired-age interval; probability >= 0.90\n");
fprintf(fid, "  SI2  review: short paired-age interval; moderate probability\n");
fprintf(fid, "  AOI  review: expected-younger age is older than the paired\n");
fprintf(fid, "       expected-older age beyond combined 2-sigma uncertainty\n");
fprintf(fid, "  AOU  review: expected age order overlaps within 2-sigma\n");
fprintf(fid, "  II   review: required uncertainty is absent/invalid\n");
fprintf(fid, "  CX   UseForModel=false; retained as context, not screened\n\n");
fprintf(fid, "The CodeID-to-code-to-definition mapping is also written as\n");
fprintf(fid, "filter_code_lookup.csv. Numeric IDs are nominal labels only.\n");
fprintf(fid, "Full values are reported in P_OlderThanReference, P_ShortInterval,\n");
fprintf(fid, "PairInterval_Ma, and PairInterval_1sigma_Ma.\n\n");

fprintf(fid, "SCREENING PARAMETERS USED IN THIS RUN\n");
fprintf(fid, "-------------------------------------\n");
fprintf(fid, "  Delta       = %.0f Ma   (short paired-age interval threshold)\n", opts.Delta);
fprintf(fid, "  P_thresh    = %.2f     (screening decision probability)\n", opts.P_thresh);
fprintf(fid, "  BoundsMethod= %s\n", opts.BoundsMethod);
fprintf(fid, "  NSigma      = %.1f     (component-sigma window multiplier)\n", opts.NSigma);
fprintf(fid, "  Kmax        = %d        (max GMM components tested by BIC)\n", opts.Kmax);
fprintf(fid, "  Nmc         = %d        (Monte Carlo draws per reference age)\n", opts.Nmc);
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

function range = resolve_target_age_range(primary)
if all(isnan(primary))
    range = [-Inf Inf];
elseif any(isnan(primary))
    error("TargetComponentAgeRange must contain two numeric bounds.");
else
    range = primary;
end
assert(range(1) < range(2), ...
    "TargetComponentAgeRange must be an increasing [younger older] interval.");
end
