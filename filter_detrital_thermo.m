function filter_detrital_thermo(file_arar, file_ap, file_zrn, outdir, Tyoung, opts)
% FILTER_DETRITAL_THERMO
% Probabilistic classification of detrital thermochronologic data relative
% to a catchment-specific youngest magmatic pulse window.
%
% Reads three files:
%   (1) Hornblende Ar/Ar      — cooling age only
%   (2) Apatite double-dated  — (U-Th)/He + U-Pb
%   (3) Zircon double-dated   — (U-Th)/He + U-Pb
%
% -----------------------------------------------------------------------
% CLASSIFICATION CLASSES  (column: class)
% Mutually exclusive, evaluated in priority order:
%   indeterminate    — missing or zero analytical uncertainty; cannot classify
%   discordant       — He age older than U-Pb age by >2-sigma; likely
%                      analytical artifact; excluded
%   flag_discordant  — He age nominally older than U-Pb age but within
%                      2-sigma; inspect before deciding
%   exclude_legacy   — P(t_cool > reference) >= P_thresh
%   flag_magmatic    — P(0 < lag < Delta) >= P_thresh
%   keep_exhumation  — passes all filters
%
% -----------------------------------------------------------------------
% REASON CODES  (column: reason_code)
% Short codes for quick table scanning. Full definitions in README.txt.
%   KE    keep_exhumation   — passes all filters
%   IN    indeterminate     — missing/zero uncertainty
%   DC    discordant        — He older than U-Pb by >2-sigma
%   FD    flag_discordant   — He nominally older than U-Pb, within 2-sigma
%   EL1   exclude_legacy    — P_legacy >= 0.90  (high-confidence)
%   EL2   exclude_legacy    — P_thresh <= P_legacy < 0.90  (moderate)
%   FM1   flag_magmatic     — P_magmatic >= 0.90  (high-confidence)
%   FM2   flag_magmatic     — P_thresh <= P_magmatic < 0.90  (moderate)
%
% The 'reason' column gives a short matching phrase. Full probability
% values are in columns P_legacy and P_magmatic.
%
% -----------------------------------------------------------------------
% APATITE U-Pb INDEPENDENT CLASSIFICATION  (columns: code_ApPb, reason_ApPb)
% -----------------------------------------------------------------------
% For apatite double-dated grains only, the U-Pb age (~450-550 C closure)
% is assessed independently as a mid-temperature cooling constraint.
% This never overrides the He classification. Codes used: KE, EL1, EL2, IN.
% Possible He + ApPb combinations and their meaning:
%   KE + KE  — both ages usable; brackets ~500 C to ~70 C cooling path
%   KE + EL  — use He age only; U-Pb predates pulse, not a valid mid-T constraint
%   EL + KE  — He excluded but U-Pb is post-pulse (rare; note in methods)
% code_ApPb and reason_ApPb are blank for HblAr and ZrnHe/ZrnPb rows.
%
% -----------------------------------------------------------------------
% OUTPUTS written to outdir
% -----------------------------------------------------------------------
%   all_data_classified.csv  — every grain with all probability, class,
%                              reason_code, reason, lag, and ApPb columns:
%                              P_legacy_ApPb, code_ApPb, class_ApPb,
%                              reason_ApPb (apatite rows only; blank elsewhere)
%   kept_strict.csv          — keep_exhumation only  [PRIMARY USE FILE]
%   kept_plus_flagged.csv    — keep_exhumation + flag_magmatic
%   excluded_legacy.csv      — exclude_legacy grains
%   discordant.csv           — hard discordant (DC)
%   flag_discordant.csv      — soft discordant (FD) — inspect before use
%   indeterminate.csv        — grains with missing uncertainties
%   summary_counts.csv       — N per system x class
%
% -----------------------------------------------------------------------
% USAGE
% -----------------------------------------------------------------------
%   filter_detrital_thermo(f_ar, f_ap, f_zrn, outdir, [80 100])
%   filter_detrital_thermo(f_ar, f_ap, f_zrn, outdir, [80 100], Delta=8)
%   filter_detrital_thermo(f_ar, f_ap, f_zrn, outdir, [80 100], ...
%       P_legacy_mode="conservative", P_thresh=0.65)
%
% -----------------------------------------------------------------------
% PARAMETERS  (name-value, all optional)
% -----------------------------------------------------------------------
%   Tyoung         [lo hi] Ma — youngest pulse window from
%                  infer_youngest_pulse_from_ZPb
%
%   Delta          Ma — lag window for magmatic flag (default 8).
%                  Grains whose cooling-to-crystallization lag is
%                  probably shorter than Delta are flagged as magmatic.
%                  Geologically: for Sierra Nevada plutons cooling at
%                  ~50-100 C/Myr, 8 Ma spans ~400 C to ~70 C.
%                  Increase to flag more grains; decrease to pass more.
%
%   P_thresh       Probability threshold for exclusion and flagging
%                  (default 0.65). A grain is excluded/flagged only if
%                  the relevant probability exceeds this value.
%                  Lower = stricter filter (fewer grains retained).
%                  See FILTER KNOBS section below for guidance.
%
%   P_legacy_mode  Which pulse boundary to use as the legacy reference:
%
%     "standard" (default)
%         Reference = Tyoung_hi (older/upper bound of pulse window).
%         Least aggressive — excludes only grains that very likely
%         predate even the oldest edge of the pulse.
%         P_legacy = P(t_cool > Tyoung_hi)
%
%     "conservative"
%         Reference = Tyoung_lo (younger/lower bound of pulse window).
%         Most aggressive — excludes grains that likely predate the
%         youngest edge of the pulse.
%         P_legacy = P(t_cool > Tyoung_lo)
%
%     "midpoint"
%         Reference = mean(Tyoung_lo, Tyoung_hi). Intermediate.
%         P_legacy = P(t_cool > midpoint)
%
%   Recommended workflow: run with "standard" first, then "conservative"
%   and compare kept_strict.csv counts. If insensitive to the choice,
%   either is defensible. If sensitive, report both in supplementary.
%
% -----------------------------------------------------------------------
% FILTER KNOBS — how to tighten or relax the filter
% -----------------------------------------------------------------------
%   Parameter        Default   Tighter value   Effect of tightening
%   P_thresh         0.65      0.50            Excludes/flags grains at
%                                              lower probability — fewer
%                                              borderline grains pass
%   Delta            8 Ma      12-15 Ma        Widens magmatic lag window;
%                                              more grains flagged magmatic
%   P_legacy_mode    standard  conservative    Shifts legacy reference to
%                                              younger pulse edge; more
%                                              grains excluded as legacy
%
%   Note: P_thresh > 0.75 is considered relaxed. A warning is printed
%   if this threshold is exceeded, prompting sensitivity testing.
%
%   Note on high-uncertainty grains: grains where 1-sigma uncertainty
%   is a large fraction of the measured age (rule of thumb: sig/age > 0.30)
%   will have P_legacy near 0.5 regardless of their measured age, and will
%   pass through as keep_exhumation with low confidence. This is the correct
%   probabilistic outcome — the code cannot exclude what it cannot constrain —
%   but it means high-uncertainty grains are not filtered out. Recommended
%   practice: pre-filter input data on relative uncertainty before running
%   the pipeline, particularly for apatite U-Pb where large uncertainties
%   are common. The confidence column in all_data_classified.csv flags
%   retained grains with low classification certainty (< 0.4 warrants
%   inspection).

arguments
    file_arar (1,1) string
    file_ap   (1,1) string
    file_zrn  (1,1) string
    outdir    (1,1) string
    Tyoung    (1,2) double          % [lo hi] Ma
    opts.Delta          (1,1) double = 8     % Ma; magmatic lag window
    opts.P_thresh       (1,1) double = 0.65  % probability threshold
    opts.P_legacy_mode  (1,1) string {mustBeMember(opts.P_legacy_mode, ...
                            ["standard","conservative","midpoint"])} = "standard"
end

% ---- Parameter warnings ----
if opts.P_thresh > 0.75
    fprintf(['  NOTE: P_thresh = %.2f is above 0.75 (relaxed). Fewer grains\n' ...
             '        will be excluded/flagged. Consider P_thresh <= 0.65 for\n' ...
             '        a stricter filter. Run with run_sensitivity=true to compare.\n'], ...
        opts.P_thresh);
end

if ~isfolder(outdir), mkdir(outdir); end

% Store params for use in sub-functions
params.Tyoung_lo        = Tyoung(1);
params.Tyoung_hi        = Tyoung(2);
params.P_legacy_exclude = opts.P_thresh;
params.P_magmatic_flag  = opts.P_thresh;
params.Delta            = opts.Delta;
params.P_legacy_mode    = opts.P_legacy_mode;

% Compute the reference age used for P_legacy based on chosen mode
switch opts.P_legacy_mode
    case "standard"
        params.legacy_ref = Tyoung(2);           % upper (older) bound
    case "conservative"
        params.legacy_ref = Tyoung(1);           % lower (younger) bound
    case "midpoint"
        params.legacy_ref = mean(Tyoung);        % midpoint
end

fprintf("  P_legacy_mode = '%s'  (reference age = %.1f Ma)\n", ...
    opts.P_legacy_mode, params.legacy_ref);

% ---- Load & validate inputs ----
T_arar = load_and_check(file_arar, ["HblGrain","HblArDate","HblAr1sigerr"]);
T_ap   = load_and_check(file_ap,   ["ApGrain","ApHeDate","ApHe2sigerr","ApPbDate","ApPb2sigerr"]);
T_zrn  = load_and_check(file_zrn,  ["ZrnGrain","ZrnHeDate","ZrnHe2sigerr","ZrnPbDate","ZrnPb2sigerr"]);

% ---- Standardize into a single table ----

% 1) Hornblende Ar/Ar — cooling age only; 1σ already reported
T_arar_std = table();
T_arar_std.Grain      = string(T_arar.HblGrain);
T_arar_std.System     = repmat("Hb_ArAr", height(T_arar), 1);
T_arar_std.t_th_Ma    = T_arar.HblArDate;
T_arar_std.sig_th_Ma  = T_arar.HblAr1sigerr;  % already 1σ
T_arar_std.t_c_Ma     = nan(height(T_arar), 1);
T_arar_std.sig_c_Ma   = nan(height(T_arar), 1);
T_arar_std.is_apatite = false(height(T_arar), 1);

% 2) Apatite double-dated — 2σ reported, convert to 1σ
T_ap_std = table();
T_ap_std.Grain      = string(T_ap.ApGrain);
T_ap_std.System     = repmat("Ap_He+Ap_UPb", height(T_ap), 1);
T_ap_std.t_th_Ma    = T_ap.ApHeDate;
T_ap_std.sig_th_Ma  = T_ap.ApHe2sigerr ./ 2;
T_ap_std.t_c_Ma     = T_ap.ApPbDate;
T_ap_std.sig_c_Ma   = T_ap.ApPb2sigerr ./ 2;
T_ap_std.is_apatite = true(height(T_ap), 1);

% 3) Zircon double-dated — 2σ reported, convert to 1σ
T_zrn_std = table();
T_zrn_std.Grain      = string(T_zrn.ZrnGrain);
T_zrn_std.System     = repmat("Zrn_He+Zrn_UPb", height(T_zrn), 1);
T_zrn_std.t_th_Ma    = T_zrn.ZrnHeDate;
T_zrn_std.sig_th_Ma  = T_zrn.ZrnHe2sigerr ./ 2;
T_zrn_std.t_c_Ma     = T_zrn.ZrnPbDate;
T_zrn_std.sig_c_Ma   = T_zrn.ZrnPb2sigerr ./ 2;
T_zrn_std.is_apatite = false(height(T_zrn), 1);

T = [T_arar_std; T_ap_std; T_zrn_std];

% ---- Classify all grains ----
T = classify_all(T, params);

% ---- Write outputs ----
writetable(T, fullfile(outdir, "all_data_classified.csv"));

write_subset(T, T.class == "keep_exhumation",                                        outdir, "kept_strict.csv");
write_subset(T, T.class == "keep_exhumation" | T.class == "flag_magmatic",           outdir, "kept_plus_flagged.csv");
write_subset(T, T.class == "exclude_legacy",                                         outdir, "excluded_legacy.csv");
write_subset(T, T.class == "discordant",                                             outdir, "discordant.csv");
write_subset(T, T.class == "flag_discordant",                                        outdir, "flag_discordant.csv");
write_subset(T, T.class == "indeterminate",                                          outdir, "indeterminate.csv");

% Summary counts
[G, sys, cls] = findgroups(T.System, T.class);
counts = splitapply(@numel, T.Grain, G);
summary_tbl = table(sys, cls, counts, 'VariableNames', {'System','class','N'});
writetable(summary_tbl, fullfile(outdir, "summary_counts.csv"));

fprintf("  Filtering complete — outputs in: %s\n", outdir);
print_summary(summary_tbl);

% Print apatite U-Pb classification summary separately
T_ap_only = T(T.is_apatite, :);
if ~isempty(T_ap_only) && any(T_ap_only.code_ApPb ~= "")
    fprintf("  Apatite U-Pb (mid-T constraint assessment):\n");
    apb_codes = unique(T_ap_only.code_ApPb(T_ap_only.code_ApPb ~= ""));
    for c = apb_codes'
        n_c = sum(T_ap_only.code_ApPb == c);
        fprintf("    %-6s  %d\n", c, n_c);
    end
    T_kept = T_ap_only(T_ap_only.class == "keep_exhumation", :);
    if ~isempty(T_kept)
        n_both = sum(T_kept.code_ApPb == "KE");
        n_excl = sum(T_kept.code_ApPb == "EL1" | T_kept.code_ApPb == "EL2");
        fprintf("    Of %d He-kept apatites: %d have usable ApPb (KE+KE), %d have legacy ApPb (KE+EL)\n", ...
            height(T_kept), n_both, n_excl);
    end
end

end

% =======================================================================
% CLASSIFICATION
% =======================================================================
function T = classify_all(T, params)

n = height(T);

% Pre-allocate output columns
T.P_legacy   = nan(n,1);
T.P_magmatic = nan(n,1);
T.lag_Ma     = nan(n,1);
T.lag_1sig   = nan(n,1);
T.class       = repmat("", n, 1);
T.reason_code = repmat("", n, 1);   % short code — see header and README
T.reason      = repmat("", n, 1);   % short phrase matching the code

% Confidence = 1 - max(P_legacy, P_magmatic) for kept grains; NaN otherwise.
T.confidence = nan(n,1);

% Apatite U-Pb independent classification — only populated for apatite rows.
T.P_legacy_ApPb = nan(n,1);
T.code_ApPb     = repmat("", n, 1);   % KE | EL1 | EL2 | IN
T.class_ApPb    = repmat("", n, 1);   % keep_exhumation | exclude_legacy | indeterminate
T.reason_ApPb   = repmat("", n, 1);

for i = 1:n
    t_th  = T.t_th_Ma(i);
    s_th  = T.sig_th_Ma(i);
    t_c   = T.t_c_Ma(i);
    s_c   = T.sig_c_Ma(i);

    has_cooling = isfinite(t_th) && isfinite(s_th) && s_th > 0;
    has_cryst   = isfinite(t_c)  && isfinite(s_c)  && s_c  > 0;

    % ---- Indeterminate: missing or zero uncertainty ----
    if ~has_cooling
        T.class(i)       = "indeterminate";
        T.reason_code(i) = "IN";
        T.reason(i)      = "missing/zero cooling age uncertainty";
        continue
    end

    % ---- Legacy probability ----
    % P(t_cool > legacy_ref): probability the cooling age predates the
    % reference point defined by P_legacy_mode.
    %   "standard"     → ref = Tyoung_hi (least aggressive exclusion)
    %   "conservative" → ref = Tyoung_lo (most aggressive exclusion)
    %   "midpoint"     → ref = mean(Tyoung_lo, Tyoung_hi)
    T.P_legacy(i) = 1 - normcdf(params.legacy_ref, t_th, s_th);

    % ---- Lag + magmatic probability ----
    if has_cryst
        lag_mu  = t_c - t_th;    % positive = crystallization older than cooling (normal)
        lag_sig = sqrt(s_c^2 + s_th^2);

        T.lag_Ma(i)   = lag_mu;
        T.lag_1sig(i) = lag_sig;

        % ---- Discordance check ----
        % Physically, t_c >= t_th must hold (crystallization must predate cooling).
        % We distinguish two cases:
        %
        %   Hard discordance: lag_mu < -2*lag_sig
        %     He age is older than U-Pb age beyond 2σ combined uncertainty.
        %     Almost certainly an analytical artifact. Excluded.
        %
        %   Soft discordance: -2*lag_sig <= lag_mu < 0
        %     He age is nominally older than U-Pb age but within 2σ.
        %     Could be analytical noise, He implantation, or partial reset.
        %     Flagged for inspection but not automatically excluded.
        if lag_mu < -2 * lag_sig
            T.class(i)       = "discordant";
            T.reason_code(i) = "DC";
            T.reason(i)      = sprintf("He > U-Pb by >2-sigma (lag %.1f +/- %.1f Ma)", ...
                lag_mu, lag_sig);
            continue
        elseif lag_mu < 0
            T.class(i)       = "flag_discordant";
            T.reason_code(i) = "FD";
            T.reason(i)      = sprintf("He nominally > U-Pb within 2-sigma (lag %.1f +/- %.1f Ma)", ...
                lag_mu, lag_sig);
            continue
        end

        % P(0 < lag < Delta): probability cooling occurred within Delta Ma
        % of crystallization, suggesting magmatic thermal relaxation rather
        % than long-term exhumational cooling.
        T.P_magmatic(i) = normcdf(params.Delta, lag_mu, lag_sig) ...
                        - normcdf(0,             lag_mu, lag_sig);
    end

    % ---- Classification (priority order) ----
    if T.P_legacy(i) >= params.P_legacy_exclude
        T.class(i) = "exclude_legacy";
        if T.P_legacy(i) >= 0.90
            T.reason_code(i) = "EL1";
            T.reason(i) = sprintf("legacy; P=%.2f [%s, ref %.1f Ma]", ...
                T.P_legacy(i), params.P_legacy_mode, params.legacy_ref);
        else
            T.reason_code(i) = "EL2";
            T.reason(i) = sprintf("legacy, moderate confidence; P=%.2f [%s, ref %.1f Ma]", ...
                T.P_legacy(i), params.P_legacy_mode, params.legacy_ref);
        end

    elseif isfinite(T.P_magmatic(i)) && T.P_magmatic(i) >= params.P_magmatic_flag
        T.class(i) = "flag_magmatic";
        if T.P_magmatic(i) >= 0.90
            T.reason_code(i) = "FM1";
            T.reason(i) = sprintf("magmatic lag; P=%.2f [Delta=%.0f Ma]", ...
                T.P_magmatic(i), params.Delta);
        else
            T.reason_code(i) = "FM2";
            T.reason(i) = sprintf("magmatic lag, moderate confidence; P=%.2f [Delta=%.0f Ma]", ...
                T.P_magmatic(i), params.Delta);
        end

    else
        T.class(i)       = "keep_exhumation";
        T.reason_code(i) = "KE";
        T.reason(i)      = "passes all filters";
        p_leg = T.P_legacy(i);
        p_mag = T.P_magmatic(i);
        if isfinite(p_mag)
            T.confidence(i) = 1 - max(p_leg, p_mag);
        else
            T.confidence(i) = 1 - p_leg;
        end
    end

    % ---- Apatite U-Pb independent classification ----
    if T.is_apatite(i) && has_cryst
        p_apb = 1 - normcdf(params.legacy_ref, t_c, s_c);
        T.P_legacy_ApPb(i) = p_apb;

        if ~isfinite(s_c) || s_c <= 0
            T.code_ApPb(i)   = "IN";
            T.reason_ApPb(i) = "missing/zero U-Pb uncertainty";
        elseif p_apb >= params.P_legacy_exclude
            if p_apb >= 0.90
                T.code_ApPb(i)   = "EL1";
                T.reason_ApPb(i) = sprintf("ApPb legacy; P=%.2f — not a valid mid-T constraint", p_apb);
            else
                T.code_ApPb(i)   = "EL2";
                T.reason_ApPb(i) = sprintf("ApPb legacy, moderate confidence; P=%.2f", p_apb);
            end
        else
            T.code_ApPb(i)   = "KE";
            T.reason_ApPb(i) = sprintf("ApPb post-pulse; P=%.2f — usable mid-T constraint", p_apb);
        end
    elseif T.is_apatite(i) && ~has_cryst
        T.code_ApPb(i)   = "IN";
        T.reason_ApPb(i) = "no ApPb age available";
    end
    % Non-apatite rows: code_ApPb and reason_ApPb stay blank
end

T.class       = string(T.class);
T.reason_code = string(T.reason_code);
T.reason      = string(T.reason);
T.code_ApPb   = string(T.code_ApPb);
T.reason_ApPb = string(T.reason_ApPb);

% Derive class_ApPb from code_ApPb so the ApPb assessment uses the same
% class vocabulary as the main 'class' column (keep_exhumation,
% exclude_legacy, indeterminate). Blank for non-apatite rows.
T.class_ApPb = repmat("", height(T), 1);
T.class_ApPb(T.code_ApPb == "KE")                              = "keep_exhumation";
T.class_ApPb(T.code_ApPb == "EL1" | T.code_ApPb == "EL2")     = "exclude_legacy";
T.class_ApPb(T.code_ApPb == "IN")                              = "indeterminate";
T.class_ApPb = string(T.class_ApPb);

end

% =======================================================================
% HELPERS
% =======================================================================
function T = load_and_check(filepath, required_cols)
% Load CSV and validate that required columns exist.
assert(isfile(filepath), "File not found: %s", filepath);
T = readtable(filepath, "VariableNamingRule","preserve");
vnames = string(T.Properties.VariableNames);
missing = required_cols(~ismember(required_cols, vnames));
if ~isempty(missing)
    error("In file '%s', missing required column(s): %s\nFound columns: %s", ...
        filepath, strjoin(missing, ", "), strjoin(vnames, ", "));
end
end

function write_subset(T, mask, outdir, filename)
% Write a subset of T to CSV; if empty write header-only file.
T_sub = T(mask, :);
writetable(T_sub, fullfile(outdir, filename));
end

function print_summary(summary_tbl)
% Print a readable summary table to the console.
classes = unique(summary_tbl.class);
systems = unique(summary_tbl.System);
fprintf("  %-22s", "");
for c = classes'
    fprintf("  %-16s", c);
end
fprintf("\n");
for s = systems'
    fprintf("  %-22s", s);
    for c = classes'
        idx = summary_tbl.System == s & summary_tbl.class == c;
        if any(idx)
            fprintf("  %-16d", summary_tbl.N(idx));
        else
            fprintf("  %-16s", "—");
        end
    end
    fprintf("\n");
end
end
