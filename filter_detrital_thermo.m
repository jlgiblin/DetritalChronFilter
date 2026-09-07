function result = filter_detrital_thermo(file_arar, file_ap, file_zrn, outdir, Tyoung, opts)
% FILTER_DETRITAL_THERMO
% Probabilistic screening of detrital thermochronologic data relative to a
% catchment-specific target zircon U-Pb age-component window.
%
% Reads three files:
%   (1) Hornblende Ar/Ar      — cooling age only
%   (2) Apatite double-dated  — (U-Th)/He + U-Pb
%   (3) Zircon double-dated   — (U-Th)/He + U-Pb
%
% -----------------------------------------------------------------------
% REFERENCE SCREEN AND REVIEW FLAGS
% The older-than-reference probability is the only exclusion rule.
%   OR1/OR2  older_than_reference       — recommended_action="exclude"
%   RT       eligible_after_reference_screen
%   II       insufficient_information   — recommended_action="review"
%
% Paired-age diagnostics never recommend exclusion:
%   SC1/SC2  short crystallization-to-cooling interval
%   AOI      cooling age older than paired U-Pb beyond combined 2-sigma
%   AOU      nominally older cooling age, overlapping within combined 2-sigma
%
% -----------------------------------------------------------------------
% CODES
%   RT    eligible_after_reference_screen
%   OR1   older_than_reference, probability >= 0.90
%   OR2   older_than_reference, moderate probability
%   SC1   short_crystallization_cooling_interval, probability >= 0.90
%   SC2   short_crystallization_cooling_interval, moderate probability
%   AOI   age_order_inconsistent
%   AOU   age_order_unresolved
%   II    insufficient_information
%
% -----------------------------------------------------------------------
% APATITE U-Pb INDEPENDENT REFERENCE ASSESSMENT
% -----------------------------------------------------------------------
% For apatite rows, the paired U-Pb age is assessed against the reference
% independently from the He result. A missing paired age is "not_applicable",
% not failed information. This assessment never overrides the He result.
%
% -----------------------------------------------------------------------
% STREAMLINED OUTPUTS written to outdir
% -----------------------------------------------------------------------
%   filter_results_full.csv — complete review table, one date per row
%   filter_results_coded.csv— compact publication table using lookup IDs
%   model_input_ages.csv    — dates eligible for downstream modeling
%   excluded_ages.csv       — only dates recommended for exclusion
%   review_flags.csv        — all dates carrying a separate review flag
%   output_summary.csv      — counts by chronometer and action
%   filter_code_lookup.csv  — ID definitions (unless WriteCodeLookup=false)
%
% -----------------------------------------------------------------------
% USAGE
% -----------------------------------------------------------------------
%   filter_detrital_thermo(f_ar, f_ap, f_zrn, outdir, [80 100])
%   filter_detrital_thermo(f_ar, f_ap, f_zrn, outdir, [80 100], Delta=8)
%   filter_detrital_thermo(f_ar, f_ap, f_zrn, outdir, [80 100], ...
%       ReferenceMode="younger_bound", P_thresh=0.65)
%
% -----------------------------------------------------------------------
% PARAMETERS  (name-value, all optional)
% -----------------------------------------------------------------------
%   Tyoung         [lo hi] Ma — target-component window from
%                  infer_youngest_pulse_from_ZPb
%
%   Delta          Ma — upper bound of the short crystallization-to-cooling
%                  interval screen (default 8). This is a sensitivity
%                  parameter, not a universal thermal-mechanism boundary.
%
%   P_thresh       Probability threshold for older-than-reference exclusion
%                  and for reporting short-interval review flags (default
%                  0.65). Lower values trigger more results.
%                  See FILTER KNOBS section below for guidance.
%
%   ReferenceMode  Reference-boundary choice:
%
%     "older_bound" (default)
%         Reference = Tyoung_hi (older edge of target-component window).
%         P_older_than_reference = P(t_cool > Tyoung_hi)
%
%     "younger_bound"
%         Reference = Tyoung_lo (younger edge of target-component window).
%         P_older_than_reference = P(t_cool > Tyoung_lo)
%
%     "midpoint"
%         Reference = mean(Tyoung_lo, Tyoung_hi).
%         P_older_than_reference = P(t_cool > midpoint)
%
% -----------------------------------------------------------------------
% FILTER KNOBS — how to tighten or relax the filter
% -----------------------------------------------------------------------
%   Parameter        Default        Alternative      Effect
%   P_thresh         0.65           0.50             More dates meet a
%                                                     decision threshold
%   Delta            8 Ma           12-15 Ma         More paired dates meet
%                                                     the short-interval flag
%   ReferenceMode    older_bound    younger_bound    Uses the younger edge;
%                                                     more dates may be older
%                                                     than the reference
%
%   Note: P_thresh > 0.75 is considered relaxed. A warning is printed
%   if this threshold is exceeded, prompting sensitivity testing.
%
%   Note on high-uncertainty dates: dates where 1-sigma uncertainty
%   is a large fraction of the measured age (rule of thumb: sig/age > 0.30)
%   may have P_older_than_reference near 0.5 and remain below the decision
%   threshold. This is the expected probabilistic outcome, but it means
%   high-uncertainty dates are not automatically excluded. Recommended
%   practice: pre-filter input data on relative uncertainty before running
%   the pipeline when appropriate, and inspect the reported probabilities.

arguments
    file_arar (1,1) string
    file_ap   (1,1) string
    file_zrn  (1,1) string
    outdir    (1,1) string
    Tyoung    (1,2) double          % [lo hi] Ma
    opts.Delta          (1,1) double = 8     % Ma; short-interval window
    opts.P_thresh       (1,1) double = 0.65  % probability threshold
    opts.ReferenceMode  (1,1) string {mustBeMember(opts.ReferenceMode, ...
                            ["older_bound","younger_bound","midpoint"])} = "older_bound"
    opts.WriteOutputs       (1,1) logical = true
    opts.WriteCodeLookup    (1,1) logical = true
end

% ---- Parameter warnings ----
if opts.P_thresh > 0.75
    fprintf(['  NOTE: P_thresh = %.2f is above 0.75 (relaxed). Fewer grains\n' ...
             '        will meet a screening condition. Consider P_thresh <= 0.65\n' ...
             '        and report a sensitivity comparison.\n'], ...
        opts.P_thresh);
end

if opts.WriteOutputs && ~isfolder(outdir), mkdir(outdir); end

% Store params for use in sub-functions
params.Tyoung_lo        = Tyoung(1);
params.Tyoung_hi        = Tyoung(2);
params.P_exclude_threshold = opts.P_thresh;
params.P_review_threshold  = opts.P_thresh;
params.Delta            = opts.Delta;
reference_mode = opts.ReferenceMode;

% Compute the reference age used by the older-than-reference calculation.
switch reference_mode
    case "older_bound"
        params.reference_age = Tyoung(2);        % upper (older) bound
        params.reference_mode = "older_bound";
    case "younger_bound"
        params.reference_age = Tyoung(1);        % lower (younger) bound
        params.reference_mode = "younger_bound";
    case "midpoint"
        params.reference_age = mean(Tyoung);     % midpoint
        params.reference_mode = "midpoint";
end

fprintf("  Reference mode = '%s'  (reference age = %.1f Ma)\n", ...
    params.reference_mode, params.reference_age);

% ---- Load & validate inputs ----
T_arar = load_and_check(file_arar, ["HblGrain","HblArDate"]);
T_ap   = load_and_check(file_ap,   ["ApGrain","ApHeDate","ApPbDate"]);
T_zrn  = load_and_check(file_zrn,  ["ZrnGrain","ZrnHeDate","ZrnPbDate"]);

% ---- Standardize into a single table ----

% 1) Hornblende Ar/Ar — cooling age only; 1σ already reported
T_arar_std = table();
T_arar_std.Grain      = string(T_arar.HblGrain);
T_arar_std.System     = repmat("Hb_ArAr", height(T_arar), 1);
T_arar_std.t_th_Ma    = T_arar.HblArDate;
T_arar_std.sig_th_Ma  = read_1sigma(T_arar, "HblAr1sigerr", "HblAr2sigerr", file_arar);
T_arar_std.t_c_Ma     = nan(height(T_arar), 1);
T_arar_std.sig_c_Ma   = nan(height(T_arar), 1);
T_arar_std.is_apatite = false(height(T_arar), 1);

% 2) Apatite double-dated — canonical public inputs are 1σ
T_ap_std = table();
T_ap_std.Grain      = string(T_ap.ApGrain);
T_ap_std.System     = repmat("Ap_He+Ap_UPb", height(T_ap), 1);
T_ap_std.t_th_Ma    = T_ap.ApHeDate;
T_ap_std.sig_th_Ma  = read_1sigma(T_ap, "ApHe1sigerr", "ApHe2sigerr", file_ap);
T_ap_std.t_c_Ma     = T_ap.ApPbDate;
T_ap_std.sig_c_Ma   = read_1sigma(T_ap, "ApPb1sigerr", "ApPb2sigerr", file_ap);
T_ap_std.is_apatite = true(height(T_ap), 1);

% 3) Zircon double-dated — canonical public inputs are 1σ
T_zrn_std = table();
T_zrn_std.Grain      = string(T_zrn.ZrnGrain);
T_zrn_std.System     = repmat("Zrn_He+Zrn_UPb", height(T_zrn), 1);
T_zrn_std.t_th_Ma    = T_zrn.ZrnHeDate;
T_zrn_std.sig_th_Ma  = read_1sigma(T_zrn, "ZrnHe1sigerr", "ZrnHe2sigerr", file_zrn);
T_zrn_std.t_c_Ma     = T_zrn.ZrnPbDate;
T_zrn_std.sig_c_Ma   = read_1sigma(T_zrn, "ZrnPb1sigerr", "ZrnPb2sigerr", file_zrn);
T_zrn_std.is_apatite = false(height(T_zrn), 1);

T = [T_arar_std; T_ap_std; T_zrn_std];

% ---- Classify all grains ----
T = classify_all(T, params);

% ---- Write neutral public outputs ----
T_public = T;

% Reader-facing table: one dated analysis per row. Cooling and paired U-Pb
% dates retain the same GrainID and PairID so their relationship is visible
% without placing two chronometers in one wide record.
T_long = make_long_results_table(T_public, params);
T_coded = make_coded_results_table(T_long);
T_model = T_long(T_long.ModelInclude, :);
T_excluded = T_long(T_long.Action == "exclude", :);
T_review = make_review_flags_table(T_long);

[G_long, chrono_long, action_long] = findgroups(T_long.Chronometer, T_long.Action);
N_long = splitapply(@numel, T_long.GrainID, G_long);
long_summary = table(chrono_long, action_long, N_long, ...
    'VariableNames', {'Chronometer','Action','N'});

% Reference-screen summary counts. Review flags are intentionally separate
% because they do not determine exclusion.
[G, sys, cls] = findgroups(T.System, T.reference_screen_class);
counts = splitapply(@numel, T.Grain, G);
summary_tbl = table(sys, cls, counts, ...
    'VariableNames', {'System','reference_screen_class','N'});

[G_review, sys_review, review_code] = findgroups(T.System, T.review_codes);
review_counts = splitapply(@numel, T.Grain, G_review);
review_summary = table(sys_review, review_code, review_counts, ...
    'VariableNames', {'System','review_codes','N'});
review_summary(review_summary.review_codes == "", :) = [];

% Return the computed tables so the pipeline can build sensitivity results
% without writing temporary per-mode folders.
result = struct();
result.filter_results = T_long;
result.coded_results = T_coded;
result.model_input_ages = T_model;
result.excluded_ages = T_excluded;
result.review_flags = T_review;
result.output_summary = long_summary;
result.system_screening_counts = summary_tbl;
result.review_counts = review_summary;
result.wide_results = T_public;
result.code_lookup = screening_code_lookup();

if opts.WriteOutputs
    writetable(T_long, fullfile(outdir, "filter_results_full.csv"));
    writetable(T_coded, fullfile(outdir, "filter_results_coded.csv"));
    writetable(T_model, fullfile(outdir, "model_input_ages.csv"));
    writetable(T_excluded, fullfile(outdir, "excluded_ages.csv"));
    writetable(T_review, fullfile(outdir, "review_flags.csv"));
    writetable(long_summary, fullfile(outdir, "output_summary.csv"));
    if opts.WriteCodeLookup
        writetable(screening_code_lookup(), ...
            fullfile(outdir, "filter_code_lookup.csv"));
    end

end

if opts.WriteOutputs
    fprintf("  Screening complete — streamlined outputs in: %s\n", outdir);
else
    fprintf("  Screening calculations complete — no per-mode files written.\n");
end
print_summary(summary_tbl);

% Print neutral apatite U-Pb screening summary separately
T_ap_only = T(T.is_apatite, :);
if ~isempty(T_ap_only) && any(T_ap_only.paired_age_reference_code ~= "")
    fprintf("  Apatite U-Pb (independent reference screening):\n");
    apb_codes = unique(T_ap_only.paired_age_reference_code(T_ap_only.paired_age_reference_code ~= ""));
    for c = apb_codes'
        n_c = sum(T_ap_only.paired_age_reference_code == c);
        fprintf("    %-6s  %d\n", c, n_c);
    end
    T_eligible = T_ap_only(T_ap_only.eligible_after_reference_screen, :);
    if ~isempty(T_eligible)
        n_candidate = sum(T_eligible.paired_age_reference_code == "RT");
        n_older = sum(T_eligible.paired_age_reference_code == "OR1" | ...
            T_eligible.paired_age_reference_code == "OR2");
        fprintf("    Of %d eligible apatite He ages: %d paired U-Pb ages pass the reference screen; %d are older than reference\n", ...
            height(T_eligible), n_candidate, n_older);
    end
end

end

% =======================================================================
% CLASSIFICATION
% =======================================================================
function T = classify_all(T, params)
n = height(T);

% Calculate only observation-based quantities used by the public policy.
T.P_older_than_reference = nan(n, 1);
T.P_short_interval = nan(n, 1);
T.crystallization_to_cooling_interval_Ma = nan(n, 1);
T.interval_1sigma_Ma = nan(n, 1);

for i = 1:n
    t_th  = T.t_th_Ma(i);
    s_th  = T.sig_th_Ma(i);
    t_c   = T.t_c_Ma(i);
    s_c   = T.sig_c_Ma(i);

    has_cooling = isfinite(t_th) && isfinite(s_th) && s_th > 0;
    has_cryst   = isfinite(t_c)  && isfinite(s_c)  && s_c  > 0;

    if has_cooling
        T.P_older_than_reference(i) = ...
            1 - normcdf(params.reference_age, t_th, s_th);
    end

    if has_cryst
        interval = t_c - t_th;
        interval_sigma = sqrt(s_c^2 + s_th^2);
        T.crystallization_to_cooling_interval_Ma(i) = interval;
        T.interval_1sigma_Ma(i) = interval_sigma;

        % The short-interval probability is only evaluated when the nominal
        % age order is nonnegative. Negative intervals receive AOI/AOU review
        % flags below and are never assigned a causal interpretation.
        if interval >= 0
            T.P_short_interval(i) = ...
                normcdf(params.Delta, interval, interval_sigma) ...
                - normcdf(0, interval, interval_sigma);
        end
    end
end

% Public policy: only the zircon-reference comparison recommends exclusion.
% Paired-age patterns are independent review flags and never exclude a row.
T = apply_public_policy(T, params);

end

function T = apply_public_policy(T, params)
% Separate the one exclusion screen from non-excluding diagnostic flags.
n = height(T);

T.reference_screen_class = repmat("", n, 1);
T.reference_screen_code = repmat("", n, 1);
T.reference_screen_basis = repmat("", n, 1);
T.eligible_after_reference_screen = false(n, 1);
T.review_recommended = false(n, 1);
T.review_codes = repmat("", n, 1);
T.review_basis = repmat("", n, 1);
T.recommended_action = repmat("", n, 1);

for i = 1:n
    has_cooling = isfinite(T.t_th_Ma(i)) && isfinite(T.sig_th_Ma(i)) && T.sig_th_Ma(i) > 0;
    has_pair = isfinite(T.t_c_Ma(i)) && isfinite(T.sig_c_Ma(i)) && T.sig_c_Ma(i) > 0;

    if ~has_cooling
        T.reference_screen_class(i) = "insufficient_information";
        T.reference_screen_code(i) = "II";
        T.reference_screen_basis(i) = "Cooling age or uncertainty is missing or invalid; reference screen not evaluated";
        T.review_recommended(i) = true;
        T.review_codes(i) = "II";
        T.review_basis(i) = T.reference_screen_basis(i);
        T.recommended_action(i) = "review";
        continue
    end

    if T.P_older_than_reference(i) >= params.P_exclude_threshold
        T.reference_screen_class(i) = "older_than_reference";
        if T.P_older_than_reference(i) >= 0.90
            T.reference_screen_code(i) = "OR1";
            T.reference_screen_basis(i) = "High probability that cooling age is older than the selected reference age";
        else
            T.reference_screen_code(i) = "OR2";
            T.reference_screen_basis(i) = "Moderate probability that cooling age is older than the selected reference age";
        end
        T.recommended_action(i) = "exclude";
    else
        T.reference_screen_class(i) = "eligible_after_reference_screen";
        T.reference_screen_code(i) = "RT";
        T.reference_screen_basis(i) = "Older-than-reference probability is below the exclusion threshold";
        T.eligible_after_reference_screen(i) = true;
        T.recommended_action(i) = "retain";
    end

    % Paired-age diagnostics add review information without changing the
    % reference-screen eligibility or an existing exclusion recommendation.
    if has_pair
        interval = T.crystallization_to_cooling_interval_Ma(i);
        interval_sigma = T.interval_1sigma_Ma(i);
        if interval < -2 * interval_sigma
            T.review_codes(i) = "AOI";
            T.review_basis(i) = "Cooling age is older than paired U-Pb age beyond combined 2-sigma uncertainty";
        elseif interval < 0
            T.review_codes(i) = "AOU";
            T.review_basis(i) = "Cooling age is nominally older than paired U-Pb age but overlaps within combined 2-sigma uncertainty";
        elseif isfinite(T.P_short_interval(i)) && T.P_short_interval(i) >= params.P_review_threshold
            if T.P_short_interval(i) >= 0.90
                T.review_codes(i) = "SC1";
                T.review_basis(i) = "High probability of a short crystallization-to-cooling interval";
            else
                T.review_codes(i) = "SC2";
                T.review_basis(i) = "Moderate probability of a short crystallization-to-cooling interval";
            end
        end
    end

    if T.review_codes(i) ~= ""
        T.review_recommended(i) = true;
        if T.recommended_action(i) ~= "exclude"
            T.recommended_action(i) = "review";
        end
    end
end

% Assess apatite U-Pb independently from the paired He result. This pass
% executes even if the He result has missing information or an age-order flag.
T.P_ApPb_older_than_reference = nan(n, 1);
T.paired_age_status = repmat("", n, 1);
T.paired_age_reference_class = repmat("", n, 1);
T.paired_age_reference_code = repmat("", n, 1);
T.paired_age_reference_basis = repmat("", n, 1);
T.paired_age_eligible_after_reference_screen = nan(n, 1);

for i = find(T.is_apatite)'
    has_age = isfinite(T.t_c_Ma(i)) && T.t_c_Ma(i) > 0;
    has_uncertainty = isfinite(T.sig_c_Ma(i)) && T.sig_c_Ma(i) > 0;
    if ~has_age
        T.paired_age_status(i) = "not_provided";
        T.paired_age_reference_class(i) = "not_applicable";
        T.paired_age_reference_code(i) = "NA";
        T.paired_age_reference_basis(i) = "No paired apatite U-Pb age was provided";
    elseif ~has_uncertainty
        T.paired_age_status(i) = "insufficient_information";
        T.paired_age_reference_class(i) = "insufficient_information";
        T.paired_age_reference_code(i) = "II";
        T.paired_age_reference_basis(i) = "Apatite U-Pb uncertainty is missing or invalid";
    else
        p_apb = 1 - normcdf(params.reference_age, T.t_c_Ma(i), T.sig_c_Ma(i));
        T.P_ApPb_older_than_reference(i) = p_apb;
        T.paired_age_status(i) = "evaluated";
        if p_apb >= params.P_exclude_threshold
            T.paired_age_reference_class(i) = "older_than_reference";
            T.paired_age_eligible_after_reference_screen(i) = 0;
            if p_apb >= 0.90
                T.paired_age_reference_code(i) = "OR1";
                T.paired_age_reference_basis(i) = "High probability that apatite U-Pb age is older than the selected reference age";
            else
                T.paired_age_reference_code(i) = "OR2";
                T.paired_age_reference_basis(i) = "Moderate probability that apatite U-Pb age is older than the selected reference age";
            end
        else
            T.paired_age_reference_class(i) = "eligible_after_reference_screen";
            T.paired_age_reference_code(i) = "RT";
            T.paired_age_reference_basis(i) = "Apatite U-Pb older-than-reference probability is below the exclusion threshold";
            T.paired_age_eligible_after_reference_screen(i) = 1;
        end
    end
end

end

% =======================================================================
% HELPERS
% =======================================================================
function T = load_and_check(filepath, required_cols)
% Load CSV and validate that required columns exist.
assert(isfile(filepath), "File not found: %s", filepath);
T = readtable(filepath, "Delimiter",",", "VariableNamingRule","preserve");
vnames = string(T.Properties.VariableNames);
missing = required_cols(~ismember(required_cols, vnames));
if ~isempty(missing)
    error("In file '%s', missing required column(s): %s\nFound columns: %s", ...
        filepath, strjoin(missing, ", "), strjoin(vnames, ", "));
end
end

function sigma1 = read_1sigma(T, canonical_name, legacy_2sigma_name, filepath)
% Read a 1σ uncertainty. Exactly one supported column must be present.
vnames = string(T.Properties.VariableNames);
has_1sigma = any(vnames == canonical_name);
has_2sigma = any(vnames == legacy_2sigma_name);

if has_1sigma && has_2sigma
    error("File '%s' contains both %s and %s. Keep only one uncertainty convention.", ...
        filepath, canonical_name, legacy_2sigma_name);
elseif has_1sigma
    sigma1 = T.(canonical_name);
elseif has_2sigma
    sigma1 = T.(legacy_2sigma_name) ./ 2;
    warning("Deprecated 2-sigma input '%s' in %s was converted to 1-sigma. Prefer '%s'.", ...
        legacy_2sigma_name, filepath, canonical_name);
else
    error("File '%s' must contain the 1-sigma uncertainty column '%s'. " + ...
        "The deprecated 2-sigma alternative '%s' is also accepted.", ...
        filepath, canonical_name, legacy_2sigma_name);
end
end

function L = make_long_results_table(T, params)
% One dated analysis per row. Paired measurements share GrainID/PairID,
% allowing a cooling age and its U-Pb age to carry separate decisions.
n = height(T);
source_row = (1:n)';
role_order = ones(n, 1);

GrainID = string(T.Grain);
System = string(T.System);
Mineral = repmat("", n, 1);
Mineral(System == "Hb_ArAr") = "hornblende";
Mineral(System == "Ap_He+Ap_UPb") = "apatite";
Mineral(System == "Zrn_He+Zrn_UPb") = "zircon";
Chronometer = repmat("", n, 1);
Chronometer(System == "Hb_ArAr") = "HblArAr";
Chronometer(System == "Ap_He+Ap_UPb") = "ApHe";
Chronometer(System == "Zrn_He+Zrn_UPb") = "ZrnHe";
PairID = repmat("not_paired", n, 1);
PairID(System == "Ap_He+Ap_UPb") = "Ap:" + GrainID(System == "Ap_He+Ap_UPb");
PairID(System == "Zrn_He+Zrn_UPb") = "Zrn:" + GrainID(System == "Zrn_He+Zrn_UPb");
AnalysisID = Chronometer + ":" + GrainID;
PairRole = repmat("cooling_age", n, 1);
PairRole(System == "Hb_ArAr") = "single_age";
has_pair_age = isfinite(T.t_c_Ma) & T.t_c_Ma > 0;
PairStatus = repmat("unpaired", n, 1);
PairStatus(has_pair_age) = "paired";
Age_Ma = T.t_th_Ma;
Age_1sigma_Ma = T.sig_th_Ma;
ReferenceMode = repmat(string(params.reference_mode), n, 1);
ReferenceAge_Ma = repmat(params.reference_age, n, 1);
P_OlderThanReference = T.P_older_than_reference;
ReferenceClass = string(T.reference_screen_class);
ReferenceCode = string(T.reference_screen_code);
ReviewRecommended = logical(T.review_recommended);
ReviewCode = string(T.review_codes);
ReviewCode(~ReviewRecommended) = "NF";
P_ShortInterval = T.P_short_interval;
PairInterval_Ma = T.crystallization_to_cooling_interval_Ma;
PairInterval_1sigma_Ma = T.interval_1sigma_Ma;
Action = string(T.recommended_action);
ModelInclude = logical(T.eligible_after_reference_screen);
ScreeningBasis = string(T.reference_screen_basis);
ReviewBasis = string(T.review_basis);
ReviewBasis(~ReviewRecommended) = "No separate review flag assigned";

C = table(source_row, role_order, GrainID, PairID, AnalysisID, System, ...
    Mineral, Chronometer, PairRole, PairStatus, Age_Ma, Age_1sigma_Ma, ...
    ReferenceMode, ReferenceAge_Ma, P_OlderThanReference, ReferenceClass, ...
    ReferenceCode, ReviewRecommended, ReviewCode, P_ShortInterval, ...
    PairInterval_Ma, PairInterval_1sigma_Ma, Action, ModelInclude, ...
    ScreeningBasis, ReviewBasis);

% Add one paired U-Pb row wherever an age is present. Apatite U-Pb uses
% its independent reference assessment. Zircon U-Pb is labeled as reference
% context and is not screened as a model-input age.
paired_rows = find(has_pair_age);
np = numel(paired_rows);
P = T(paired_rows, :);
source_row = paired_rows;
role_order = 2 * ones(np, 1);
GrainID = string(P.Grain);
System = string(P.System);
Mineral = repmat("zircon", np, 1);
Mineral(P.is_apatite) = "apatite";
Chronometer = repmat("ZrnUPb", np, 1);
Chronometer(P.is_apatite) = "ApUPb";
PairID = repmat("Zrn:", np, 1) + GrainID;
PairID(P.is_apatite) = "Ap:" + GrainID(P.is_apatite);
AnalysisID = Chronometer + ":" + GrainID;
PairRole = repmat("paired_u_pb_age", np, 1);
PairStatus = repmat("paired", np, 1);
Age_Ma = P.t_c_Ma;
Age_1sigma_Ma = P.sig_c_Ma;
ReferenceMode = repmat(string(params.reference_mode), np, 1);
ReferenceAge_Ma = repmat(params.reference_age, np, 1);
P_OlderThanReference = nan(np, 1);
ReferenceClass = repmat("reference_context_not_screened", np, 1);
ReferenceCode = repmat("RC", np, 1);
ReviewRecommended = false(np, 1);
ReviewCode = repmat("NF", np, 1);
P_ShortInterval = nan(np, 1);
PairInterval_Ma = nan(np, 1);
PairInterval_1sigma_Ma = nan(np, 1);
Action = repmat("reference_only", np, 1);
ModelInclude = false(np, 1);
ScreeningBasis = repmat( ...
    "Zircon U-Pb age provides target-component reference context and is not screened as a model-input age", ...
    np, 1);
ReviewBasis = repmat("No separate review flag assigned", np, 1);

is_ap = logical(P.is_apatite);
P_OlderThanReference(is_ap) = P.P_ApPb_older_than_reference(is_ap);
ReferenceClass(is_ap) = string(P.paired_age_reference_class(is_ap));
ReferenceCode(is_ap) = string(P.paired_age_reference_code(is_ap));
ScreeningBasis(is_ap) = string(P.paired_age_reference_basis(is_ap));
ap_eligible = is_ap & P.paired_age_eligible_after_reference_screen == 1;
ap_older = is_ap & string(P.paired_age_reference_class) == "older_than_reference";
ap_review = is_ap & string(P.paired_age_reference_class) == "insufficient_information";
Action(is_ap) = "not_applicable";
Action(ap_eligible) = "retain";
Action(ap_older) = "exclude";
Action(ap_review) = "review";
ModelInclude(ap_eligible) = true;
ReviewRecommended(ap_review) = true;
ReviewCode(ap_review) = "II";
ReviewBasis(ap_review) = string(P.paired_age_reference_basis(ap_review));

U = table(source_row, role_order, GrainID, PairID, AnalysisID, System, ...
    Mineral, Chronometer, PairRole, PairStatus, Age_Ma, Age_1sigma_Ma, ...
    ReferenceMode, ReferenceAge_Ma, P_OlderThanReference, ReferenceClass, ...
    ReferenceCode, ReviewRecommended, ReviewCode, P_ShortInterval, ...
    PairInterval_Ma, PairInterval_1sigma_Ma, Action, ModelInclude, ...
    ScreeningBasis, ReviewBasis);

L = [C; U];
L = sortrows(L, {'source_row','role_order'});
L = removevars(L, {'source_row','role_order'});
L = add_action_explanations(L, params);
end

function L = add_action_explanations(L, params)
% Put the numerical rule and its consequence in one plain-language field.
% This complements, rather than replaces, the shorter machine-readable
% ScreeningBasis and ReviewBasis fields.
n = height(L);
P_Threshold = repmat(params.P_exclude_threshold, n, 1);
ShortIntervalThreshold_Ma = nan(n, 1);
ShortIntervalThreshold_Ma(isfinite(L.P_ShortInterval)) = params.Delta;
ActionReason = strings(n, 1);

for i = 1:n
    action = string(L.Action(i));
    p_old = L.P_OlderThanReference(i);
    ref = L.ReferenceAge_Ma(i);
    review_basis = string(L.ReviewBasis(i));
    screening_basis = string(L.ScreeningBasis(i));

    if action == "exclude"
        ActionReason(i) = sprintf( ...
            "Exclude: P(age older than %.2f Ma) = %.3f meets the %.2f decision threshold.", ...
            ref, p_old, params.P_exclude_threshold);
        if L.ReviewRecommended(i) && strlength(review_basis) > 0
            ActionReason(i) = ActionReason(i) + ...
                " Separate non-excluding review flag: " + review_basis + ".";
        end
    elseif action == "retain"
        ActionReason(i) = sprintf( ...
            "Retain: P(age older than %.2f Ma) = %.3f is below the %.2f decision threshold.", ...
            ref, p_old, params.P_exclude_threshold);
    elseif action == "review" && L.ModelInclude(i)
        if isfinite(p_old)
            prefix = string(sprintf( ...
                "Include for modeling: P(age older than %.2f Ma) = %.3f is below the %.2f decision threshold. Review flag: ", ...
                ref, p_old, params.P_exclude_threshold));
        else
            prefix = "Include for modeling, with review flag: ";
        end
        if strlength(review_basis) > 0
            ActionReason(i) = prefix + review_basis + ".";
        else
            ActionReason(i) = prefix + screening_basis + ".";
        end
    elseif action == "review"
        ActionReason(i) = "Review: " + screening_basis + ".";
    elseif action == "reference_only"
        ActionReason(i) = "Reference context only; this zircon U-Pb date is not a model-input decision.";
        P_Threshold(i) = NaN;
    elseif action == "not_applicable"
        ActionReason(i) = "No screening action applies to this row.";
        P_Threshold(i) = NaN;
    else
        ActionReason(i) = screening_basis;
    end
end

L = addvars(L, P_Threshold, ShortIntervalThreshold_Ma, ActionReason, ...
    'After', 'ReviewBasis');
end

function R = make_review_flags_table(L)
% A self-contained subset for users who want to inspect only flagged dates.
% For paired analyses, include the related date directly beside the flag.
R = L(L.ReviewRecommended, :);
n = height(R);
RelatedChronometer = repmat("not_paired", n, 1);
RelatedAge_Ma = nan(n, 1);
RelatedAge_1sigma_Ma = nan(n, 1);

for i = 1:n
    if R.PairStatus(i) == "unpaired"
        continue
    end
    related = find(L.PairID == R.PairID(i) & L.AnalysisID ~= R.AnalysisID(i));
    if numel(related) == 1
        RelatedChronometer(i) = L.Chronometer(related);
        RelatedAge_Ma(i) = L.Age_Ma(related);
        RelatedAge_1sigma_Ma(i) = L.Age_1sigma_Ma(related);
    end
end

R = addvars(R, RelatedChronometer, RelatedAge_Ma, RelatedAge_1sigma_Ma, ...
    'After', 'Age_1sigma_Ma');
end

function C = make_coded_results_table(L)
% Publication-oriented one-date-per-row table. Long definitions and prose
% are replaced by nominal numeric IDs defined in filter_code_lookup.csv.
lookup = screening_code_lookup();
ReferenceResultID = nan(height(L), 1);
ReviewFlagID = zeros(height(L), 1); % 0 = no review flag
for i = 1:height(lookup)
    ReferenceResultID(L.ReferenceCode == lookup.Code(i)) = lookup.CodeID(i);
    ReviewFlagID(L.ReviewCode == lookup.Code(i)) = lookup.CodeID(i);
end

C = table(L.GrainID, L.PairID, L.AnalysisID, L.Chronometer, L.PairRole, ...
    L.Age_Ma, L.Age_1sigma_Ma, L.ReferenceAge_Ma, ...
    L.P_OlderThanReference, L.P_ShortInterval, L.PairInterval_Ma, ...
    ReferenceResultID, ReviewFlagID, L.Action, L.ModelInclude, ...
    'VariableNames', {'GrainID','PairID','AnalysisID','Chronometer','PairRole', ...
    'Age_Ma','Age_1sigma_Ma','ReferenceAge_Ma', ...
    'P_OlderThanReference','P_ShortInterval','PairInterval_Ma', ...
    'ReferenceResultID','ReviewFlagID','Action','ModelInclude'});
end

function L = screening_code_lookup()
% Machine-readable and publication-ready definitions for compact outputs.
CodeID = (0:10)';
Code = ["NF";"RT";"OR1";"OR2";"SC1";"SC2";"AOI";"AOU";"II";"NA";"RC"];
Role = ["review";"reference";"reference";"reference";"review";"review"; ...
    "review";"review";"reference_or_review";"paired_age_status";"reference_context"];
Meaning = [ ...
    "no_review_flag";
    "eligible_after_reference_screen";
    "older_than_reference";
    "older_than_reference";
    "short_crystallization_cooling_interval";
    "short_crystallization_cooling_interval";
    "age_order_inconsistent";
    "age_order_unresolved";
    "insufficient_information";
    "not_applicable";
    "reference_context_not_screened"];
Definition = [ ...
    "No separate review flag was assigned";
    "Older-than-reference probability is below the exclusion threshold";
    "High probability that cooling age is older than the selected reference age";
    "Moderate probability that cooling age is older than the selected reference age";
    "High probability of a short crystallization-to-cooling interval";
    "Moderate probability of a short crystallization-to-cooling interval";
    "Cooling age is older than the paired age beyond combined 2-sigma uncertainty";
    "Cooling age is nominally older than the paired age but overlaps within combined 2-sigma uncertainty";
    "Required age or uncertainty is missing or nonpositive";
    "Paired age was not provided, so paired-age assessment is not applicable";
    "Zircon U-Pb age provides target-component reference context and is not screened as a model-input age"];
DefaultAction = ["none";"retain";"exclude";"exclude";"review";"review"; ...
    "review";"review";"review";"not_applicable";"reference_only"];
IdentifierNote = repmat("CodeID is a nominal identifier, not a rank or continuous quantity", 11, 1);
L = table(CodeID, Code, Role, Meaning, Definition, DefaultAction, IdentifierNote);
end

function print_summary(summary_tbl)
% Print a readable summary table to the console.
if ismember("reference_screen_class", string(summary_tbl.Properties.VariableNames))
    class_values = summary_tbl.reference_screen_class;
elseif ismember("screening_class", string(summary_tbl.Properties.VariableNames))
    class_values = summary_tbl.screening_class;
else
    class_values = summary_tbl.class;
end
classes = unique(class_values);
systems = unique(summary_tbl.System);
fprintf("  %-22s", "");
for c = classes'
    fprintf("  %-16s", c);
end
fprintf("\n");
for s = systems'
    fprintf("  %-22s", s);
    for c = classes'
        idx = summary_tbl.System == s & class_values == c;
        if any(idx)
            fprintf("  %-16d", summary_tbl.N(idx));
        else
            fprintf("  %-16s", "—");
        end
    end
    fprintf("\n");
end
end
