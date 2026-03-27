function run_detrital_pipeline(catchments_root, outdir_root, opts)
% RUN_DETRITAL_PIPELINE
% Discovers all catchment subfolders under catchments_root, runs
% infer_youngest_pulse_from_ZPb then filter_detrital_thermo for each,
% and writes cross-catchment summary and sensitivity CSVs.
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
%     README.txt                    — guide to all outputs and reason codes
%     _pipeline_summary/
%       pipeline_summary.csv        — pulse window + GMM stats per catchment
%       sensitivity_by_system.csv   — keep_exhumation N under standard vs
%                                     conservative vs midpoint (if run_sensitivity)
%     <CatchmentName>/
%       ZPb_QA/                     — GMM QA plots for pulse inference
%       standard/                   — filtered outputs using standard mode
%         kept_strict.csv           *** PRIMARY USE FILE ***
%         kept_plus_flagged.csv
%         all_data_classified.csv
%         excluded_legacy.csv
%         discordant.csv
%         flag_discordant.csv
%         indeterminate.csv
%         summary_counts.csv
%       conservative/               — same outputs, conservative mode
%       midpoint/                   — same outputs, midpoint mode
%
% The standard/kept_strict.csv file is the recommended input to TSF
% and Pecube workflows. The conservative/ and midpoint/ folders exist
% for sensitivity testing only and are not used in primary analysis.
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
%   Delta           magmatic lag window in Ma (default 8)
%                   Increase to flag more grains as magmatic
%   P_thresh        probability threshold for exclusion/flagging (default 0.65)
%                   Lower = stricter filter; P_thresh > 0.75 triggers a warning
%   Nmc             Monte Carlo draws per grain for ZPb GMM (default 50)
%   run_sensitivity true (default): runs standard + conservative + midpoint
%                   modes and writes sensitivity_by_system.csv.
%                   false: runs standard mode only (faster).
%
% -----------------------------------------------------------------------
% FILTER KNOBS — quick reference
% -----------------------------------------------------------------------
%   Parameter        Default   Tighter        Effect of tightening
%   P_thresh         0.65      0.50           Fewer borderline grains pass
%   Delta            8 Ma      12-15 Ma       More grains flagged magmatic
%   P_legacy_mode    standard  conservative   More grains excluded as legacy

arguments
    catchments_root  (1,1) string  = "Catchments"
    outdir_root      (1,1) string  = "Output"
    opts.Kmax             (1,1) double  = 6
    opts.K_override       (1,1) double  = 0
    opts.K_override_map                = []
    opts.BoundsMethod     (1,1) string  = "gmm_ci"   % "gmm_ci" | "quantile"
    opts.NSigma           (1,1) double  = 1.0        % sigma multiplier for window + model start
    opts.Delta            (1,1) double  = 8
    opts.P_thresh         (1,1) double  = 0.65
    opts.Nmc              (1,1) double  = 50
    opts.run_sensitivity  (1,1) logical = true
end

if ~isfolder(catchments_root)
    error("Catchments root folder not found: %s", catchments_root);
end
if ~isfolder(outdir_root), mkdir(outdir_root); end

% Summary outputs go in a dedicated subfolder so they are clearly
% separated from per-catchment filtered data.
summary_dir = fullfile(outdir_root, "_pipeline_summary");
if ~isfolder(summary_dir), mkdir(summary_dir); end

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
    modes = ["standard", "conservative", "midpoint"];
else
    modes = "standard";
end

% ---- Storage ----
summary_rows = cell(numel(catchment_names), 1);
% counts_store{i, mi} = summary_counts table for catchment i, mode mi
counts_store = cell(numel(catchment_names), numel(modes));

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
                      "BoundsMethod", opts.BoundsMethod, "NSigma", opts.NSigma};
        if K_use > 0
            pulse_opts = [pulse_opts, {"K_override", K_use}]; %#ok<AGROW>
        end
        pulse = infer_youngest_pulse_from_ZPb(f_zpb, ...
            fullfile(outdir_root, cname, "ZPb_QA"), pulse_opts{:});
        fprintf("  Pulse window:     %.1f - %.1f Ma  (K=%d, BIC=%d)\n", ...
            pulse.Tyoung(1), pulse.Tyoung(2), pulse.K_used, pulse.K_bic);
        fprintf("  Model start time: %.1f +/- %.1f Ma  (mu=%.1f + %.0f*sigma=%.1f)\n", ...
            pulse.model_start_Ma, pulse.model_start_err_Ma, ...
            pulse.mu_young, pulse.NSigma_used, pulse.sigma_young);
    catch ME
        warning("Pulse inference failed for %s: %s", cname, ME.message);
        summary_rows{i} = make_error_row(cname);
        continue
    end

    % -- Step 2: filter for each mode --
    for mi = 1:numel(modes)
        mode = modes(mi);
        cout = fullfile(outdir_root, cname, mode);
        try
            filter_detrital_thermo(f_ar, f_ap, f_zrn, cout, pulse.Tyoung, ...
                "Delta",         opts.Delta, ...
                "P_thresh",      opts.P_thresh, ...
                "P_legacy_mode", mode);

            % Read back summary_counts.csv for the sensitivity table
            sc_path = fullfile(cout, "summary_counts.csv");
            if isfile(sc_path)
                counts_store{i, mi} = readtable(sc_path, "VariableNamingRule","preserve");
            end
        catch ME
            warning("Filtering failed for %s (mode=%s): %s", cname, mode, ME.message);
        end
    end

    % -- Collect pulse summary row --
    summary_rows{i} = table(cname, pulse.Tyoung(1), pulse.Tyoung(2), ...
        pulse.mu_young, pulse.sigma_young, pulse.weight_young, ...
        pulse.model_start_Ma, pulse.model_start_err_Ma, pulse.NSigma_used, ...
        string(pulse.bounds_method), pulse.K_used, pulse.K_bic, ...
        pulse.Nages, pulse.Nselected, ...
        'VariableNames', {'Catchment','Tyoung_lo','Tyoung_hi', ...
        'mu_young','sigma_young','weight_young', ...
        'model_start_Ma','model_start_err_Ma','NSigma','BoundsMethod', ...
        'K_used','K_bic_selected','N_ZPb_ages','N_ZPb_assigned'});

    fprintf("  Done.\n");
end

% ---- Write pipeline summary ----
valid_rows = summary_rows(~cellfun(@isempty, summary_rows));
if ~isempty(valid_rows)
    summary_tbl = vertcat(valid_rows{:});
    writetable(summary_tbl, fullfile(summary_dir, "pipeline_summary.csv"));
    fprintf("\nPipeline summary written to: %s\n", ...
        fullfile(summary_dir, "pipeline_summary.csv"));
else
    warning("No catchments processed successfully.");
    return
end

% ---- Build and write sensitivity table ----
if opts.run_sensitivity && numel(modes) > 1
    sens_tbl = build_sensitivity_table(catchment_names, modes, counts_store);
    if ~isempty(sens_tbl)
        sens_path = fullfile(summary_dir, "sensitivity_by_system.csv");
        writetable(sens_tbl, sens_path);
        fprintf("Sensitivity table written to: %s\n", sens_path);
        print_sensitivity_summary(sens_tbl);
    end
end

end % main function

% =======================================================================
% SENSITIVITY TABLE BUILDER
% =======================================================================
function tbl = build_sensitivity_table(catchment_names, modes, counts_store)
% For each catchment × system, extract keep_exhumation N under each mode,
% compute absolute and % differences between standard and conservative.

rows = {};

for i = 1:numel(catchment_names)
    cname = catchment_names(i);

    % Gather all system names seen across any mode for this catchment
    all_systems = string([]);
    for mi = 1:numel(modes)
        sc = counts_store{i, mi};
        if ~isempty(sc)
            all_systems = union(all_systems, string(sc.System));
        end
    end

    for si = 1:numel(all_systems)
        sys = all_systems(si);

        % Extract keep_exhumation N for each mode
        n_vals = nan(1, numel(modes));
        n_total = NaN;
        for mi = 1:numel(modes)
            sc = counts_store{i, mi};
            if isempty(sc), continue; end
            sc_sys = sc(string(sc.System) == sys, :);
            % Total N = sum across all classes for this system (from any mode)
            if mi == 1 && ~isempty(sc_sys)
                n_total = sum(sc_sys.N);
            end
            sc_keep = sc_sys(string(sc_sys.class) == "keep_exhumation", :);
            if ~isempty(sc_keep)
                n_vals(mi) = sc_keep.N(1);
            else
                % System present but zero kept grains — record as 0, not NaN
                if ~isempty(sc_sys)
                    n_vals(mi) = 0;
                end
            end
        end

        n_std = n_vals(1);   % standard
        n_con = n_vals(2);   % conservative
        n_mid = n_vals(3);   % midpoint

        abs_diff = n_std - n_con;
        if isfinite(n_std) && n_std > 0
            pct_diff = 100 * abs_diff / n_std;
        elseif isfinite(n_std) && n_std == 0
            pct_diff = 0;
        else
            pct_diff = NaN;
        end

        % Sensitivity flag based on % difference standard vs conservative
        if ~isfinite(pct_diff)
            sens_flag = "indeterminate";
        elseif abs(pct_diff) > 20
            sens_flag = "HIGH";
        elseif abs(pct_diff) > 5
            sens_flag = "moderate";
        else
            sens_flag = "low";
        end

        rows{end+1} = table( ...
            cname, sys, n_total, ...
            n_std, n_con, n_mid, ...
            abs_diff, pct_diff, sens_flag, ...
            'VariableNames', { ...
                'Catchment', 'System', 'N_total_grains', ...
                'N_keep_standard', 'N_keep_conservative', 'N_keep_midpoint', ...
                'Diff_std_minus_con', 'PctDiff_std_minus_con', ...
                'Sensitivity'}); %#ok<AGROW>
    end
end

if isempty(rows)
    tbl = table();
    return
end

tbl = vertcat(rows{:});

% Sort: HIGH sensitivity first, then by catchment and system
% Use a numeric priority key to avoid string-length sorting artifacts
sens_priority = zeros(height(tbl), 1);
sens_priority(tbl.Sensitivity == "HIGH")          = 1;
sens_priority(tbl.Sensitivity == "moderate")      = 2;
sens_priority(tbl.Sensitivity == "low")           = 3;
sens_priority(tbl.Sensitivity == "indeterminate") = 4;
[~, ord] = sortrows([sens_priority, double(categorical(tbl.Catchment)), ...
                     double(categorical(tbl.System))]);
tbl = tbl(ord, :);

end

% =======================================================================
% CONSOLE SENSITIVITY SUMMARY
% =======================================================================
function print_sensitivity_summary(tbl)

if isempty(tbl), return; end

fprintf("\n--- Sensitivity summary (standard vs conservative) ---\n");
fprintf("  %-10s  %-20s  %7s  %7s  %7s  %9s  %8s  %s\n", ...
    "Catchment", "System", "N_total", "N_std", "N_con", "N_mid", "PctDiff", "Sensitivity");
fprintf("  %s\n", repmat("-", 1, 88));

for r = 1:height(tbl)
    n_mid_str = "-";
    if isfinite(tbl.N_keep_midpoint(r))
        n_mid_str = sprintf("%d", tbl.N_keep_midpoint(r));
    end
    n_tot_str = "-";
    if isfinite(tbl.N_total_grains(r))
        n_tot_str = sprintf("%d", tbl.N_total_grains(r));
    end
    pct_str = "-";
    if isfinite(tbl.PctDiff_std_minus_con(r))
        pct_str = sprintf("%+.1f%%", tbl.PctDiff_std_minus_con(r));
    end
    fprintf("  %-10s  %-20s  %7s  %7d  %7d  %9s  %8s  %s\n", ...
        tbl.Catchment(r), tbl.System(r), n_tot_str, ...
        tbl.N_keep_standard(r), tbl.N_keep_conservative(r), ...
        n_mid_str, pct_str, tbl.Sensitivity(r));
end
fprintf("\n");

end

% =======================================================================
% HELPERS
% =======================================================================
function row = make_error_row(cname)
row = table(cname, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, "failed", NaN, NaN, NaN, NaN, ...
    'VariableNames', {'Catchment','Tyoung_lo','Tyoung_hi', ...
    'mu_young','sigma_young','weight_young', ...
    'model_start_Ma','model_start_err_Ma','NSigma','BoundsMethod', ...
    'K_used','K_bic_selected','N_ZPb_ages','N_ZPb_assigned'});
end

% =======================================================================
% README WRITER
% =======================================================================
function write_readme(outdir_root, opts)
% Writes README.txt to outdir_root explaining the folder structure,
% which files to use, all reason codes, and the filter parameters used
% in this run. Overwrites any previous README so it stays current.

fid = fopen(fullfile(outdir_root, "README.txt"), "w");
fprintf(fid, "MultichronFitTSF — Detrital Thermochronometry Filter Output\n");
fprintf(fid, "============================================================\n\n");

fprintf(fid, "Generated by run_detrital_pipeline.m\n");
fprintf(fid, "Date: %s\n\n", datestr(now, "yyyy-mm-dd HH:MM:SS")); %#ok<TNOW1,DATST>

fprintf(fid, "---------------------------------------------------------------\n");
fprintf(fid, "WHICH FILES TO USE\n");
fprintf(fid, "---------------------------------------------------------------\n");
fprintf(fid, "Primary analysis:\n");
fprintf(fid, "  <CatchmentName>/standard/kept_strict.csv\n");
fprintf(fid, "    Contains keep_exhumation grains only. This is the file\n");
fprintf(fid, "    to use as input to TSF and Pecube workflows.\n\n");
fprintf(fid, "Sensitivity testing (do not use for primary analysis):\n");
fprintf(fid, "  <CatchmentName>/conservative/kept_strict.csv\n");
fprintf(fid, "  <CatchmentName>/midpoint/kept_strict.csv\n");
fprintf(fid, "    Same filter applied with different legacy reference ages.\n");
fprintf(fid, "    Compare grain counts against standard/ to assess sensitivity.\n\n");
fprintf(fid, "Cross-catchment summaries:\n");
fprintf(fid, "  _pipeline_summary/pipeline_summary.csv\n");
fprintf(fid, "    Pulse window and GMM parameters for each catchment.\n");
fprintf(fid, "  _pipeline_summary/sensitivity_by_system.csv\n");
fprintf(fid, "    keep_exhumation N under standard vs conservative vs midpoint,\n");
fprintf(fid, "    with sensitivity flags per catchment and thermochronologic system.\n\n");
fprintf(fid, "QA plots:\n");
fprintf(fid, "  <CatchmentName>/ZPb_QA/\n");
fprintf(fid, "    BIC curve and GMM fit to ZrnPb ages. Inspect to verify K\n");
fprintf(fid, "    selection and pulse window before using results.\n\n");

fprintf(fid, "---------------------------------------------------------------\n");
fprintf(fid, "REASON CODES  (columns: reason_code, code_ApPb)\n");
fprintf(fid, "---------------------------------------------------------------\n");
fprintf(fid, "Each grain in all_data_classified.csv carries a reason_code\n");
fprintf(fid, "indicating why it was assigned its class. The 'reason' column\n");
fprintf(fid, "gives a short phrase; full probability values are in P_legacy\n");
fprintf(fid, "and P_magmatic. The confidence threshold split is fixed at 0.90.\n\n");
fprintf(fid, "  Code  Class              Condition\n");
fprintf(fid, "  ----  -----------------  ----------------------------------------\n");
fprintf(fid, "  KE    keep_exhumation    passes all filters\n");
fprintf(fid, "  IN    indeterminate      missing or zero analytical uncertainty\n");
fprintf(fid, "  DC    discordant         He older than U-Pb age by >2-sigma;\n");
fprintf(fid, "                           likely analytical artifact\n");
fprintf(fid, "  FD    flag_discordant    He nominally older than U-Pb within 2-sigma;\n");
fprintf(fid, "                           inspect: could be noise, implantation, or reset\n");
fprintf(fid, "  EL1   exclude_legacy     P_legacy >= 0.90  (high-confidence legacy)\n");
fprintf(fid, "  EL2   exclude_legacy     P_thresh <= P_legacy < 0.90  (moderate)\n");
fprintf(fid, "  FM1   flag_magmatic      P_magmatic >= 0.90  (high-confidence magmatic)\n");
fprintf(fid, "  FM2   flag_magmatic      P_thresh <= P_magmatic < 0.90  (moderate)\n\n");
fprintf(fid, "Apatite U-Pb codes (code_ApPb / class_ApPb) — independent of He classification:\n");
fprintf(fid, "  code_ApPb is the short code; class_ApPb is the full class name.\n");
fprintf(fid, "  KE  / keep_exhumation  U-Pb post-pulse; usable as mid-T constraint (~500 C)\n");
fprintf(fid, "  EL1 / exclude_legacy   U-Pb likely pre-pulse, high confidence (P >= 0.90)\n");
fprintf(fid, "  EL2 / exclude_legacy   U-Pb likely pre-pulse, moderate confidence\n");
fprintf(fid, "  IN  / indeterminate    missing/zero U-Pb uncertainty\n");
fprintf(fid, "  (blank)                non-apatite system — not applicable\n\n");
fprintf(fid, "Interpreting He + ApPb code combinations:\n");
fprintf(fid, "  KE + KE   both ages usable; brackets ~500 C to ~70 C cooling path\n");
fprintf(fid, "  KE + EL   He age retained; U-Pb predates pulse, not a valid mid-T\n");
fprintf(fid, "            constraint (do not include ApPb in thermal model)\n");
fprintf(fid, "  EL + KE   He excluded; U-Pb post-pulse (rare — note in methods)\n\n");

fprintf(fid, "---------------------------------------------------------------\n");
fprintf(fid, "FILTER PARAMETERS USED IN THIS RUN\n");
fprintf(fid, "---------------------------------------------------------------\n");
fprintf(fid, "  Delta       = %.0f Ma   (magmatic lag window)\n", opts.Delta);
fprintf(fid, "  P_thresh    = %.2f     (exclusion/flag probability threshold)\n", opts.P_thresh);
fprintf(fid, "  BoundsMethod= %s\n", opts.BoundsMethod);
fprintf(fid, "  NSigma      = %.1f     (sigma multiplier for pulse window and model start)\n", opts.NSigma);
fprintf(fid, "  Kmax        = %d        (max GMM components tested by BIC)\n", opts.Kmax);
fprintf(fid, "  Nmc         = %d        (Monte Carlo draws per grain for ZPb GMM)\n", opts.Nmc);
if opts.run_sensitivity
    fprintf(fid, "  Sensitivity = true      (standard + conservative + midpoint modes run)\n\n");
else
    fprintf(fid, "  Sensitivity = false     (standard mode only)\n\n");
end
fprintf(fid, "P_legacy_mode is always 'standard' for the standard/ folder,\n");
fprintf(fid, "'conservative' for conservative/, and 'midpoint' for midpoint/.\n\n");

fprintf(fid, "---------------------------------------------------------------\n");
fprintf(fid, "FILTER KNOBS — adjusting strictness\n");
fprintf(fid, "---------------------------------------------------------------\n");
fprintf(fid, "  Parameter       Default   Tighter   Effect of tightening\n");
fprintf(fid, "  P_thresh        0.65      0.50      Excludes/flags more borderline grains\n");
fprintf(fid, "  Delta           8 Ma      12+ Ma    Flags more grains as magmatic\n");
fprintf(fid, "  P_legacy_mode   standard  conserv.  Removes more legacy grains\n\n");
fprintf(fid, "A warning is printed if P_thresh > 0.75 (relaxed setting).\n\n");

fprintf(fid, "---------------------------------------------------------------\n");
fprintf(fid, "REFERENCE\n");
fprintf(fid, "---------------------------------------------------------------\n");
fprintf(fid, "Thermochronometric Scaling Function (TSF) methodology:\n");
fprintf(fid, "  Gallagher & Parra (2020), EPSL\n");
fprintf(fid, "This multi-chronometer extension:\n");
fprintf(fid, "  Giblin et al. (in prep)\n");

fclose(fid);
fprintf("  README written to: %s\n", fullfile(outdir_root, "README.txt"));
end
