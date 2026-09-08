function result = filter_detrital_thermo(chronometer_csv, outdir, target_window, opts)
% FILTER_DETRITAL_THERMO
% Probabilistically screens a flexible, long-format chronometer table
% relative to a user-selected reference-component age window.
%
% REQUIRED INPUT COLUMNS
%   Chronometer, GrainID, Age_Ma, Age_1sigma_Ma
%
% OPTIONAL INPUT COLUMNS
%   PairID       Shared identifier for two analyses of the same grain.
%   PairRole     Blank, "expected_younger", or "expected_older".
%   UseForModel  true/false; defaults to true when omitted or blank.
%
% Chronometer is a grouping label only. The program never infers closure
% temperature, mineral, kinetics, or expected age order from its text.
% Pair-order diagnostics are evaluated only when a PairID occurs exactly
% twice with one expected_younger and one expected_older role.
%
% REFERENCE SCREEN
%   OR1/OR2  older_than_reference       action="exclude"
%   RT       eligible_after_reference_screen
%   II       insufficient_information   action="review"
%   CX       context_only_not_screened  action="context_only"
%
% OPTIONAL PAIR REVIEW FLAGS
%   SI1/SI2  short_pair_interval
%   AOI      age_order_inconsistent
%   AOU      age_order_unresolved
%
% Pair flags never recommend exclusion. Only OR1 and OR2 do so.
%
% USAGE
%   filter_detrital_thermo("ChronometerData.csv", "filter_output", [80 100])
%   filter_detrital_thermo("ChronometerData.csv", "filter_output", ...
%       [80 100], Delta=8, P_thresh=0.65)

arguments
    chronometer_csv (1,1) string
    outdir          (1,1) string
    target_window   (1,2) double
    opts.Delta          (1,1) double = 8
    opts.P_thresh       (1,1) double = 0.65
    opts.ReferenceMode  (1,1) string {mustBeMember(opts.ReferenceMode, ...
                            ["older_bound","younger_bound","midpoint"])} = "older_bound"
    opts.WriteOutputs       (1,1) logical = true
    opts.WriteCodeLookup    (1,1) logical = true
end

if ~(all(isfinite(target_window)) && target_window(1) < target_window(2))
    error("target_window must be an increasing [younger older] interval.");
end
if opts.P_thresh > 0.75
    fprintf(['  NOTE: P_thresh = %.2f is above 0.75. Fewer dates will meet\n' ...
        '        a screening condition; consider reporting a sensitivity test.\n'], ...
        opts.P_thresh);
end
if opts.WriteOutputs && ~isfolder(outdir), mkdir(outdir); end

params.target_window_lo = target_window(1);
params.target_window_hi = target_window(2);
params.P_exclude_threshold = opts.P_thresh;
params.P_review_threshold = opts.P_thresh;
params.Delta = opts.Delta;
params.reference_mode = opts.ReferenceMode;
switch opts.ReferenceMode
    case "older_bound"
        params.reference_age = target_window(2);
    case "younger_bound"
        params.reference_age = target_window(1);
    case "midpoint"
        params.reference_age = mean(target_window);
end

fprintf("  Reference mode = '%s'  (reference age = %.1f Ma)\n", ...
    params.reference_mode, params.reference_age);

T = load_chronometer_data(chronometer_csv);
L = classify_observations(T, params);
L = add_action_explanations(L, params);

T_coded = make_coded_results_table(L);
T_model = L(L.ModelInclude, :);
T_excluded = L(L.Action == "exclude", :);
T_review = make_review_flags_table(L);
output_summary = make_output_summary(L);

[G, chronometer, reference_class] = findgroups(L.Chronometer, L.ReferenceClass);
N = splitapply(@numel, L.GrainID, G);
screening_counts = table(chronometer, reference_class, N, ...
    'VariableNames', {'Chronometer','ReferenceClass','N'});

review_rows = L.ReviewRecommended;
if any(review_rows)
    [G_review, review_chronometer, review_code] = findgroups( ...
        L.Chronometer(review_rows), L.ReviewCode(review_rows));
    review_N = splitapply(@numel, L.GrainID(review_rows), G_review);
    review_counts = table(review_chronometer, review_code, review_N, ...
        'VariableNames', {'Chronometer','ReviewCode','N'});
else
    review_counts = table(strings(0,1), strings(0,1), zeros(0,1), ...
        'VariableNames', {'Chronometer','ReviewCode','N'});
end

result = struct();
result.filter_results = L;
result.coded_results = T_coded;
result.model_input_ages = T_model;
result.excluded_ages = T_excluded;
result.review_flags = T_review;
result.output_summary = output_summary;
result.system_screening_counts = screening_counts;
result.review_counts = review_counts;
result.code_lookup = screening_code_lookup();

if opts.WriteOutputs
    writetable(L, fullfile(outdir, "filter_results_full.csv"));
    writetable(T_coded, fullfile(outdir, "filter_results_coded.csv"));
    writetable(T_model, fullfile(outdir, "model_input_ages.csv"));
    writetable(T_excluded, fullfile(outdir, "excluded_ages.csv"));
    writetable(T_review, fullfile(outdir, "review_flags.csv"));
    writetable(output_summary, fullfile(outdir, "output_summary.csv"));
    if opts.WriteCodeLookup
        writetable(screening_code_lookup(), ...
            fullfile(outdir, "filter_code_lookup.csv"));
    end
end

if opts.WriteOutputs
    fprintf("  Screening complete — outputs in: %s\n", outdir);
else
    fprintf("  Screening calculations complete — no per-mode files written.\n");
end
print_summary(screening_counts);
end

% ======================================================================
% INPUT
% ======================================================================
function T = load_chronometer_data(filepath)
assert(isfile(filepath), "File not found: %s", filepath);
T_in = readtable(filepath, "Delimiter",",", ...
    "VariableNamingRule","preserve", "TextType","string");
vnames = string(T_in.Properties.VariableNames);
required = ["Chronometer","GrainID","Age_Ma","Age_1sigma_Ma"];
missing = required(~ismember(required, vnames));
if ~isempty(missing)
    error("File '%s' is missing required column(s): %s. Found: %s", ...
        filepath, strjoin(missing, ", "), strjoin(vnames, ", "));
end
if any(endsWith(vnames, "2sigerr", "IgnoreCase", true))
    error("File '%s' contains a 2-sigma uncertainty column. Convert all uncertainties to 1-sigma and use Age_1sigma_Ma.", filepath);
end

n = height(T_in);
Chronometer = normalize_text(T_in.Chronometer);
GrainID = normalize_text(T_in.GrainID);
Age_Ma = numeric_column(T_in.Age_Ma, "Age_Ma", filepath);
Age_1sigma_Ma = numeric_column(T_in.Age_1sigma_Ma, "Age_1sigma_Ma", filepath);

if any(strlength(Chronometer) == 0)
    error("Chronometer must be nonblank in every row of %s.", filepath);
end
if any(strlength(GrainID) == 0)
    error("GrainID must be nonblank in every row of %s.", filepath);
end

if ismember("PairID", vnames)
    PairID = normalize_text(T_in.PairID);
else
    PairID = repmat("", n, 1);
end
if ismember("PairRole", vnames)
    PairRole = lower(normalize_text(T_in.PairRole));
else
    PairRole = repmat("", n, 1);
end
if ismember("UseForModel", vnames)
    UseForModel = parse_model_flag(T_in.UseForModel, filepath);
else
    UseForModel = true(n, 1);
end

allowed_roles = ["","expected_younger","expected_older"];
bad_role = ~ismember(PairRole, allowed_roles);
if any(bad_role)
    error("Unsupported PairRole '%s' in %s. Use expected_younger, expected_older, or blank.", ...
        PairRole(find(bad_role, 1)), filepath);
end
if any(strlength(PairRole) > 0 & strlength(PairID) == 0)
    error("PairRole requires a nonblank PairID in %s.", filepath);
end

AnalysisID = Chronometer + ":" + GrainID;
if numel(unique(AnalysisID)) ~= n
    [unique_ids, ~, group_id] = unique(AnalysisID);
    id_counts = accumarray(group_id, 1);
    duplicates = unique_ids(id_counts > 1);
    error("Duplicate Chronometer + GrainID identifier in %s: %s", ...
        filepath, duplicates(1));
end

PairStatus = repmat("unpaired", n, 1);
paired_ids = unique(PairID(strlength(PairID) > 0), "stable");
for pair_id = paired_ids'
    rows = find(PairID == pair_id);
    if numel(rows) ~= 2
        error("PairID '%s' must occur exactly twice in %s; found %d rows.", ...
            pair_id, filepath, numel(rows));
    end
    roles = PairRole(rows);
    if all(strlength(roles) == 0)
        PairStatus(rows) = "paired_without_order";
    elseif nnz(roles == "expected_younger") == 1 && ...
            nnz(roles == "expected_older") == 1
        PairStatus(rows) = "paired_with_order";
    else
        error("PairID '%s' must have one expected_younger and one expected_older PairRole, or two blank roles.", pair_id);
    end
end

PairID(strlength(PairID) == 0) = "not_paired";
T = table(GrainID, PairID, AnalysisID, Chronometer, PairRole, PairStatus, ...
    UseForModel, Age_Ma, Age_1sigma_Ma);
end

function values = normalize_text(values)
values = strtrim(string(values));
values(ismissing(values)) = "";
end

function values = numeric_column(values, column_name, filepath)
if isnumeric(values)
    values = double(values(:));
else
    values = str2double(string(values(:)));
end
if numel(values) == 0
    error("Column %s is empty in %s.", column_name, filepath);
end
end

function flags = parse_model_flag(values, filepath)
if islogical(values)
    flags = values(:);
    return
elseif isnumeric(values)
    values = values(:);
    bad = ~(values == 0 | values == 1 | isnan(values));
    if any(bad)
        error("UseForModel must contain true/false, yes/no, 1/0, or blank in %s.", filepath);
    end
    flags = true(size(values));
    flags(values == 0) = false;
    return
end

text_values = lower(strtrim(string(values(:))));
text_values(ismissing(text_values)) = "";
flags = true(size(text_values));
flags(ismember(text_values, ["false","no","0"])) = false;
valid = ismember(text_values, ["","true","false","yes","no","1","0"]);
if any(~valid)
    error("UseForModel contains unsupported value '%s' in %s.", ...
        text_values(find(~valid, 1)), filepath);
end
end

% ======================================================================
% CLASSIFICATION
% ======================================================================
function L = classify_observations(T, params)
L = T;
n = height(L);
L.ReferenceMode = repmat(string(params.reference_mode), n, 1);
L.ReferenceAge_Ma = repmat(params.reference_age, n, 1);
L.P_OlderThanReference = nan(n, 1);
L.ReferenceClass = repmat("", n, 1);
L.ReferenceCode = repmat("", n, 1);
L.ReviewRecommended = false(n, 1);
L.ReviewCode = repmat("NF", n, 1);
L.P_ShortInterval = nan(n, 1);
L.PairInterval_Ma = nan(n, 1);
L.PairInterval_1sigma_Ma = nan(n, 1);
L.Action = repmat("", n, 1);
L.ModelInclude = false(n, 1);
L.ScreeningBasis = repmat("", n, 1);
L.ReviewBasis = repmat("No separate review flag assigned", n, 1);

valid_age = isfinite(L.Age_Ma) & L.Age_Ma > 0;
valid_uncertainty = isfinite(L.Age_1sigma_Ma) & L.Age_1sigma_Ma > 0;
valid_measurement = valid_age & valid_uncertainty;

for i = 1:n
    if ~L.UseForModel(i)
        L.ReferenceClass(i) = "context_only_not_screened";
        L.ReferenceCode(i) = "CX";
        L.Action(i) = "context_only";
        L.ScreeningBasis(i) = "UseForModel is false; the age is retained only as pair or contextual information";
        if ~valid_measurement(i)
            L.ReviewRecommended(i) = true;
            L.ReviewCode(i) = "II";
            L.ReviewBasis(i) = "Context age or 1-sigma uncertainty is missing or nonpositive";
        end
        continue
    end

    if ~valid_measurement(i)
        L.ReferenceClass(i) = "insufficient_information";
        L.ReferenceCode(i) = "II";
        L.ReviewRecommended(i) = true;
        L.ReviewCode(i) = "II";
        L.Action(i) = "review";
        L.ScreeningBasis(i) = "Age or 1-sigma uncertainty is missing or nonpositive; reference screen not evaluated";
        L.ReviewBasis(i) = L.ScreeningBasis(i);
        continue
    end

    p_older = 1 - normcdf(params.reference_age, ...
        L.Age_Ma(i), L.Age_1sigma_Ma(i));
    L.P_OlderThanReference(i) = p_older;
    if p_older >= params.P_exclude_threshold
        L.ReferenceClass(i) = "older_than_reference";
        L.Action(i) = "exclude";
        if p_older >= 0.90
            L.ReferenceCode(i) = "OR1";
            L.ScreeningBasis(i) = "High probability that the age is older than the selected reference age";
        else
            L.ReferenceCode(i) = "OR2";
            L.ScreeningBasis(i) = "Moderate probability that the age is older than the selected reference age";
        end
    else
        L.ReferenceClass(i) = "eligible_after_reference_screen";
        L.ReferenceCode(i) = "RT";
        L.Action(i) = "retain";
        L.ModelInclude(i) = true;
        L.ScreeningBasis(i) = "Older-than-reference probability is below the exclusion threshold";
    end
end

ordered_pairs = unique(L.PairID(L.PairStatus == "paired_with_order"), "stable");
for pair_id = ordered_pairs'
    younger = find(L.PairID == pair_id & L.PairRole == "expected_younger");
    older = find(L.PairID == pair_id & L.PairRole == "expected_older");
    pair_rows = [younger; older];

    if numel(younger) ~= 1 || numel(older) ~= 1
        error("Internal pair validation failed for PairID '%s'.", pair_id);
    end
    if ~all(valid_measurement(pair_rows))
        if valid_measurement(younger) && L.ReviewCode(younger) == "NF"
            L.ReviewRecommended(younger) = true;
            L.ReviewCode(younger) = "II";
            L.ReviewBasis(younger) = "Paired age or 1-sigma uncertainty is missing or nonpositive";
            if L.Action(younger) == "retain"
                L.Action(younger) = "review";
            end
        end
        continue
    end

    interval = L.Age_Ma(older) - L.Age_Ma(younger);
    interval_sigma = sqrt(L.Age_1sigma_Ma(older)^2 + ...
        L.Age_1sigma_Ma(younger)^2);
    L.PairInterval_Ma(pair_rows) = interval;
    L.PairInterval_1sigma_Ma(pair_rows) = interval_sigma;
    if interval >= 0
        p_short = normcdf(params.Delta, interval, interval_sigma) - ...
            normcdf(0, interval, interval_sigma);
        L.P_ShortInterval(pair_rows) = p_short;
    else
        p_short = NaN;
    end

    review_code = "";
    review_basis = "";
    if interval < -2 * interval_sigma
        review_code = "AOI";
        review_basis = "Expected-younger age is older than expected-older age beyond combined 2-sigma uncertainty";
    elseif interval < 0
        review_code = "AOU";
        review_basis = "Expected-younger age is nominally older than expected-older age but overlaps within combined 2-sigma uncertainty";
    elseif isfinite(p_short) && p_short >= params.P_review_threshold
        if p_short >= 0.90
            review_code = "SI1";
            review_basis = "High probability of a short interval between the paired ages";
        else
            review_code = "SI2";
            review_basis = "Moderate probability of a short interval between the paired ages";
        end
    end

    if strlength(review_code) > 0
        L.ReviewRecommended(younger) = true;
        L.ReviewCode(younger) = review_code;
        L.ReviewBasis(younger) = review_basis;
        if L.Action(younger) == "retain"
            L.Action(younger) = "review";
        end
    end
end
end

% ======================================================================
% OUTPUT TABLES
% ======================================================================
function S = make_output_summary(L)
chronometers = unique(string(L.Chronometer), 'stable');
n = numel(chronometers);
N_ReportedRows = zeros(n, 1);
N_ValidAges = zeros(n, 1);
ObservedMinAge_Ma = nan(n, 1);
ObservedMaxAge_Ma = nan(n, 1);
N_ModelInput = zeros(n, 1);
ModelInputMinAge_Ma = nan(n, 1);
ModelInputMedianAge_Ma = nan(n, 1);
ModelInputMaxAge_Ma = nan(n, 1);
N_Excluded = zeros(n, 1);
N_ReviewFlagged = zeros(n, 1);
N_ContextOnly = zeros(n, 1);

for i = 1:n
    rows = string(L.Chronometer) == chronometers(i);
    valid = rows & isfinite(L.Age_Ma) & L.Age_Ma > 0;
    model = valid & logical(L.ModelInclude);
    observed_ages = L.Age_Ma(valid);
    model_ages = L.Age_Ma(model);
    N_ReportedRows(i) = nnz(rows);
    N_ValidAges(i) = nnz(valid);
    if ~isempty(observed_ages)
        ObservedMinAge_Ma(i) = min(observed_ages);
        ObservedMaxAge_Ma(i) = max(observed_ages);
    end
    N_ModelInput(i) = nnz(model);
    if ~isempty(model_ages)
        ModelInputMinAge_Ma(i) = min(model_ages);
        ModelInputMedianAge_Ma(i) = median(model_ages);
        ModelInputMaxAge_Ma(i) = max(model_ages);
    end
    N_Excluded(i) = nnz(rows & L.Action == "exclude");
    N_ReviewFlagged(i) = nnz(rows & L.ReviewRecommended);
    N_ContextOnly(i) = nnz(rows & L.Action == "context_only");
end

S = table(chronometers, N_ReportedRows, N_ValidAges, ...
    ObservedMinAge_Ma, ObservedMaxAge_Ma, N_ModelInput, ...
    ModelInputMinAge_Ma, ModelInputMedianAge_Ma, ModelInputMaxAge_Ma, ...
    N_Excluded, N_ReviewFlagged, N_ContextOnly, ...
    'VariableNames', {'Chronometer','N_ReportedRows','N_ValidAges', ...
    'ObservedMinAge_Ma','ObservedMaxAge_Ma','N_ModelInput', ...
    'ModelInputMinAge_Ma','ModelInputMedianAge_Ma','ModelInputMaxAge_Ma', ...
    'N_Excluded','N_ReviewFlagged','N_ContextOnly'});
end

function L = add_action_explanations(L, params)
n = height(L);
P_Threshold = repmat(params.P_exclude_threshold, n, 1);
ShortIntervalThreshold_Ma = nan(n, 1);
ShortIntervalThreshold_Ma(isfinite(L.P_ShortInterval)) = params.Delta;
ActionReason = strings(n, 1);

for i = 1:n
    action = L.Action(i);
    p_old = L.P_OlderThanReference(i);
    ref = L.ReferenceAge_Ma(i);
    if action == "exclude"
        ActionReason(i) = sprintf( ...
            "Exclude: P(age older than %.2f Ma) = %.3f meets the %.2f decision threshold.", ...
            ref, p_old, params.P_exclude_threshold);
        if L.ReviewRecommended(i)
            ActionReason(i) = ActionReason(i) + ...
                " Separate non-excluding review flag: " + L.ReviewBasis(i) + ".";
        end
    elseif action == "retain"
        ActionReason(i) = sprintf( ...
            "Retain: P(age older than %.2f Ma) = %.3f is below the %.2f decision threshold.", ...
            ref, p_old, params.P_exclude_threshold);
    elseif action == "review" && L.ModelInclude(i)
        ActionReason(i) = sprintf( ...
            "Include for modeling: P(age older than %.2f Ma) = %.3f is below the %.2f decision threshold. Review flag: %s.", ...
            ref, p_old, params.P_exclude_threshold, L.ReviewBasis(i));
    elseif action == "review"
        ActionReason(i) = "Review: " + L.ReviewBasis(i) + ".";
    elseif action == "context_only"
        ActionReason(i) = "Context only: UseForModel is false; this age is not a model-input decision.";
        if L.ReviewRecommended(i)
            ActionReason(i) = ActionReason(i) + " Review flag: " + L.ReviewBasis(i) + ".";
        end
        P_Threshold(i) = NaN;
    else
        ActionReason(i) = L.ScreeningBasis(i);
    end
end

L = addvars(L, P_Threshold, ShortIntervalThreshold_Ma, ActionReason, ...
    'After', 'ReviewBasis');
end

function R = make_review_flags_table(L)
R = L(L.ReviewRecommended, :);
n = height(R);
RelatedChronometer = repmat("not_paired", n, 1);
RelatedAge_Ma = nan(n, 1);
RelatedAge_1sigma_Ma = nan(n, 1);
for i = 1:n
    if R.PairID(i) == "not_paired"
        continue
    end
    related = find(L.PairID == R.PairID(i) & ...
        L.AnalysisID ~= R.AnalysisID(i));
    if isscalar(related)
        RelatedChronometer(i) = L.Chronometer(related);
        RelatedAge_Ma(i) = L.Age_Ma(related);
        RelatedAge_1sigma_Ma(i) = L.Age_1sigma_Ma(related);
    end
end
R = addvars(R, RelatedChronometer, RelatedAge_Ma, RelatedAge_1sigma_Ma, ...
    'After', 'Age_1sigma_Ma');
end

function C = make_coded_results_table(L)
lookup = screening_code_lookup();
ReferenceResultID = nan(height(L), 1);
ReviewFlagID = zeros(height(L), 1);
for i = 1:height(lookup)
    ReferenceResultID(L.ReferenceCode == lookup.Code(i)) = lookup.CodeID(i);
    ReviewFlagID(L.ReviewCode == lookup.Code(i)) = lookup.CodeID(i);
end
C = table(L.GrainID, L.PairID, L.AnalysisID, L.Chronometer, ...
    L.PairRole, L.PairStatus, L.UseForModel, L.Age_Ma, L.Age_1sigma_Ma, ...
    L.ReferenceAge_Ma, L.P_OlderThanReference, L.P_ShortInterval, ...
    L.PairInterval_Ma, ReferenceResultID, ReviewFlagID, L.Action, ...
    L.ModelInclude, ...
    'VariableNames', {'GrainID','PairID','AnalysisID','Chronometer', ...
    'PairRole','PairStatus','UseForModel','Age_Ma','Age_1sigma_Ma', ...
    'ReferenceAge_Ma','P_OlderThanReference','P_ShortInterval', ...
    'PairInterval_Ma','ReferenceResultID','ReviewFlagID','Action', ...
    'ModelInclude'});
end

function L = screening_code_lookup()
CodeID = (0:9)';
Code = ["NF";"RT";"OR1";"OR2";"SI1";"SI2";"AOI";"AOU";"II";"CX"];
Role = ["review";"reference";"reference";"reference";"review"; ...
    "review";"review";"review";"reference_or_review";"context"];
Meaning = [ ...
    "no_review_flag";
    "eligible_after_reference_screen";
    "older_than_reference";
    "older_than_reference";
    "short_pair_interval";
    "short_pair_interval";
    "age_order_inconsistent";
    "age_order_unresolved";
    "insufficient_information";
    "context_only_not_screened"];
Definition = [ ...
    "No separate review flag was assigned";
    "Older-than-reference probability is below the exclusion threshold";
    "High probability that the age is older than the selected reference age";
    "Moderate probability that the age is older than the selected reference age";
    "High probability of a short interval between paired ages";
    "Moderate probability of a short interval between paired ages";
    "Expected-younger age is older than expected-older age beyond combined 2-sigma uncertainty";
    "Expected-younger age is nominally older than expected-older age but overlaps within combined 2-sigma uncertainty";
    "Required age or 1-sigma uncertainty is missing or nonpositive";
    "UseForModel is false; age is retained only as pair or contextual information"];
DefaultAction = ["none";"retain";"exclude";"exclude";"review"; ...
    "review";"review";"review";"review";"context_only"];
IdentifierNote = repmat( ...
    "CodeID is a nominal identifier, not a rank or continuous quantity", 10, 1);
L = table(CodeID, Code, Role, Meaning, Definition, DefaultAction, IdentifierNote);
end

function print_summary(summary_tbl)
classes = unique(summary_tbl.ReferenceClass, "stable");
chronometers = unique(summary_tbl.Chronometer, "stable");
fprintf("  %-22s", "");
for c = classes'
    fprintf("  %-16s", c);
end
fprintf("\n");
for chronometer = chronometers'
    fprintf("  %-22s", chronometer);
    for c = classes'
        idx = summary_tbl.Chronometer == chronometer & ...
            summary_tbl.ReferenceClass == c;
        if any(idx)
            fprintf("  %-16d", summary_tbl.N(idx));
        else
            fprintf("  %-16s", "—");
        end
    end
    fprintf("\n");
end
end
