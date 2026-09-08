function out = infer_target_component(reference_csv, outdir, opts)
% INFER_TARGET_COMPONENT
% Selects the youngest eligible component from a Gaussian mixture model
% (GMM) fit to a user-defined reference age distribution.
%
% Its mean (mu) and component standard deviation (sigma) define a target
% age window and a candidate model-start age. Interpreting this statistical
% component geologically requires independent evidence.
%
% Uses Monte Carlo jittering to propagate the required 1-sigma analytical
% uncertainties in Age_1sigma_Ma into GMM fitting.
%
% K selection: BIC over K = 2..Kmax by default. Use K_override to fix K
% after inspecting QA plots.
%
% -----------------------------------------------------------------------
% BOUNDS METHODS  (BoundsMethod parameter)
% -----------------------------------------------------------------------
%   "gmm_sigma_window" (default)
%       Window = [mu - NSigma*sigma,  mu + NSigma*sigma]
%       Uses the GMM component's own mean and std to define the window.
%       Robust to outlier grains at the tails of the assigned population.
%       NSigma=1.0 (default) gives ~68% of the component distribution;
%       NSigma=2.0 gives ~95%. The upper bound (mu + NSigma*sigma) also
%       is also reported as a candidate model-start age for QTQt/Pecube.
%       Its suitability must be evaluated for the specific model.
%       Minimum window floor: 2 Ma (numerical guard only; the GMM sigma
%       already encodes real component width, so the floor is rarely active).
%
%   "quantile"
%       Window = quantile(assigned_ages, BoundsQuantiles) + Buffer.
%       Original method. More sensitive to outlier grains at the tails.
%       Use if the GMM sigma is poorly constrained (small N in component).
%
% -----------------------------------------------------------------------
% REQUIRED COLUMNS
% -----------------------------------------------------------------------
%   ReferenceSystem, GrainID, Age_Ma, Age_1sigma_Ma
%
% TargetComponentAgeRange limits which fitted GMM component can be selected
% as the target. The GMM is still fit to the complete age distribution, so
% components outside the range can be represented but cannot be selected.
%
% -----------------------------------------------------------------------
% OUTPUT STRUCT
% -----------------------------------------------------------------------
%   out.reference_system  label supplied in ReferenceSystem
%   out.target_window_Ma  [lo hi] Ma target window passed to screening
%   out.target_component_mean_Ma   GMM component mean (Ma)
%   out.target_component_sigma_Ma  GMM component standard deviation (Ma)
%   out.target_component_weight    GMM component mixture weight
%   out.candidate_model_start_Ma   component mean + NSigma * component sigma
%   out.NSigma_used       NSigma value used to compute model_start_Ma
%   out.bounds_method     BoundsMethod used
%   out.K_used            K used for final GMM
%   out.K_bic             K selected by BIC; NaN when K is overridden
%   out.K_selection_method  "bic" or "manual_override"
%   out.K_override_applied  whether K_override controlled the final fit
%   out.Nages             number of valid input ages
%   out.Nselected         ages assigned to the selected component
%   out.QAplot            path to saved QA figure
%   out.used_errors       logical — were analytical uncertainties incorporated?
%   out.Nmc_used          Monte Carlo draws used (0 if no errors)
%   out.target_component_age_range  eligible component-mean range [lo hi] Ma
%
% -----------------------------------------------------------------------
% USAGE
% -----------------------------------------------------------------------
%   out = infer_target_component("SampleA/ReferenceDistribution.csv", "SampleA/QA")
%   out = infer_target_component("SampleA/ReferenceDistribution.csv", "SampleA/QA", NSigma=1.5)
%   out = infer_target_component("SampleA/ReferenceDistribution.csv", "SampleA/QA", K_override=3)
%   out = infer_target_component("SampleA/ReferenceDistribution.csv", "SampleA/QA", BoundsMethod="quantile")
%   out = infer_target_component("SampleA/ReferenceDistribution.csv", "SampleA/QA", ...
%       TargetComponentAgeRange=[50 200])

arguments
    reference_csv  (1,1) string
    outdir   (1,1) string  = "target_component"
    opts.Kmax             (1,1) double = 6
    opts.K_override       (1,1) double = 0
    opts.MinWeight        (1,1) double = 0.05
    opts.MemberProbThresh (1,1) double = 0.7
    opts.TargetComponentAgeRange (1,2) double = [NaN NaN]
    opts.BoundsMethod     (1,1) string ...
                          {mustBeMember(opts.BoundsMethod, ...
                          ["gmm_sigma_window","quantile"])} = "gmm_sigma_window"
    opts.NSigma           (1,1) double = 1.0
    opts.BoundsQuantiles  (1,2) double = [0.02 0.98]  % used by quantile method only
    opts.MinimumWindowWidth (1,1) double = NaN
    opts.Buffer           (1,1) double = 2     % Ma; added to each bound (quantile only)
    opts.Bandwidth        (1,1) double = NaN   % KDE bandwidth; NaN = auto
    opts.Nmc              (1,1) double = 50
    opts.MaxTotalSamples  (1,1) double = 200000
end

target_age_range = resolve_target_age_range(opts.TargetComponentAgeRange);

% Resolve minimum window width: sigma window uses 2 Ma (numerical floor; the
% sigma already encodes component width). quantile uses 10 Ma (guards against
% absurdly narrow empirical windows when few grains are assigned).
if isnan(opts.MinimumWindowWidth)
    if opts.BoundsMethod == "gmm_sigma_window"
        minWidth = 2;
    else
        minWidth = 10;
    end
else
    minWidth = opts.MinimumWindowWidth;
end

if ~isfolder(outdir), mkdir(outdir); end

% ---- Load and validate ----
T = readtable(reference_csv, "Delimiter",",", "VariableNamingRule","preserve", ...
    "TextType","string");
vnames = string(T.Properties.VariableNames);

required = ["ReferenceSystem","GrainID","Age_Ma","Age_1sigma_Ma"];
missing = required(~ismember(required, vnames));
if ~isempty(missing)
    error("File '%s' is missing required column(s): %s. Found: %s", ...
        reference_csv, strjoin(missing, ", "), strjoin(vnames, ", "));
end
if any(endsWith(vnames, "2sigerr", "IgnoreCase", true))
    error("File '%s' contains a 2-sigma uncertainty column. Convert all uncertainties to 1-sigma and use Age_1sigma_Ma.", reference_csv);
end

reference_labels = unique(strtrim(string(T.ReferenceSystem)), "stable");
reference_labels(ismissing(reference_labels) | strlength(reference_labels) == 0) = [];
if numel(reference_labels) ~= 1
    error("ReferenceSystem must contain one nonblank system label in %s.", reference_csv);
end
reference_system = reference_labels(1);

grain_ids = strtrim(string(T.GrainID));
grain_ids(ismissing(grain_ids)) = "";
if any(strlength(grain_ids) == 0)
    error("GrainID must be nonblank in every row of %s.", reference_csv);
end
if numel(unique(grain_ids)) ~= numel(grain_ids)
    error("GrainID values must be unique in %s.", reference_csv);
end

ages = to_numeric(T.Age_Ma);
sig_input = to_numeric(T.Age_1sigma_Ma);

% Basic cleaning
valid = isfinite(ages) & ages > 0;
ages  = ages(valid);
sig_input = sig_input(valid);

assert(~isempty(ages), "No valid positive ages found in %s", reference_csv);
if any(~isfinite(sig_input) | sig_input <= 0)
    error("Every valid reference age in %s must have a positive Age_1sigma_Ma value.", reference_csv);
end
N = numel(ages);

% ---- Monte Carlo jittering to propagate uncertainties into GMM ----
Nmc = max(1, round(opts.Nmc));
total = N * Nmc;
if total > opts.MaxTotalSamples
    Nmc = max(1, floor(opts.MaxTotalSamples / N));
end
rng(0);
A = repmat(ages, 1, Nmc) + repmat(sig_input, 1, Nmc) .* randn(N, Nmc);
z = A(:);

z = z(isfinite(z) & z > 0);

% ---- BIC-based K selection ----
Kmin = 2;
Kmax = max(Kmin, opts.Kmax);

if opts.K_override > 0
    % Manual override: fit only the requested K
    K_bic  = NaN;
    K_use  = opts.K_override;
    K_selection_method = "manual_override";
    K_override_applied = true;
    gm     = fit_gmm(z, K_use);
    bic_vals = NaN(Kmax - Kmin + 1, 1);  % empty for QA plot
    K_range  = (Kmin:Kmax)';
else
    % BIC sweep
    K_range  = (Kmin:Kmax)';
    bic_vals = NaN(numel(K_range), 1);
    gm_store = cell(numel(K_range), 1);

    for ki = 1:numel(K_range)
        try
            gm_store{ki} = fit_gmm(z, K_range(ki));
            bic_vals(ki) = gm_store{ki}.BIC;
        catch
            bic_vals(ki) = Inf;
        end
    end

    [~, best_ki] = min(bic_vals);
    K_bic = K_range(best_ki);
    K_use = K_bic;
    K_selection_method = "bic";
    K_override_applied = false;
    gm    = gm_store{best_ki};
end

% ---- Extract youngest eligible component within the requested age range ----
mu  = gm.mu(:);
sig = sqrt(squeeze(gm.Sigma));
w   = gm.ComponentProportion(:);

inAgeRange = mu >= target_age_range(1) & mu <= target_age_range(2);
assert(any(inAgeRange), ...
    "No fitted GMM component mean falls within TargetComponentAgeRange [%.1f, %.1f] Ma.", ...
    target_age_range(1), target_age_range(2));

validComp   = w >= opts.MinWeight & inAgeRange;
mu_masked   = mu;
mu_masked(~validComp) = Inf;
[mu_y, idx_y] = min(mu_masked);
if isinf(mu_y)
    mu_masked = mu;
    mu_masked(~inAgeRange) = Inf;
    [mu_y, idx_y] = min(mu_masked);  % fallback: in-range component below MinWeight
    warning("All GMM components within TargetComponentAgeRange are below MinWeight; " + ...
        "using the youngest in-range component.");
end
sig_y = sig(idx_y);
w_y   = w(idx_y);

if all(isfinite(target_age_range))
    fprintf("  Target-component eligibility range: %.1f - %.1f Ma\n", ...
        target_age_range(1), target_age_range(2));
end

% Posterior membership on original ages — used for Nselected reporting
% regardless of BoundsMethod, and for the quantile method's grain selection.
post = posterior(gm, ages(:));
p_y  = post(:, idx_y);

% ---- Define target-component window ----
% model_start_Ma is always mu + NSigma*sigma regardless of BoundsMethod.
% This candidate age is the selected-component mean plus NSigma times the
% component standard deviation. The component sigma is dispersion, not an
% uncertainty on the candidate age.
model_start_Ma = mu_y + opts.NSigma * sig_y;

switch opts.BoundsMethod

    case "gmm_sigma_window"
        % Window defined symmetrically from GMM component parameters.
        % Not sensitive to outlier grains at the tails.
        lo = mu_y - opts.NSigma * sig_y;
        hi = mu_y + opts.NSigma * sig_y;
        fprintf("  BoundsMethod = gmm_sigma_window | mu=%.1f, sigma=%.1f, NSigma=%.1f\n", ...
            mu_y, sig_y, opts.NSigma);

    case "quantile"
        % Original method: empirical quantiles of posterior-assigned grains.
        % Requires enough grains assigned to the selected component.
        sel = p_y >= opts.MemberProbThresh;

        if nnz(sel) >= 10
            ages_y = ages(sel);
        else
            ages_y = ages(ages >= (mu_y - 2*sig_y) & ages <= (mu_y + 2*sig_y));
            if numel(ages_y) < 10
                ages_y = ages;
                warning("infer_target_component: few grains assigned to the selected " + ...
                    "component — target bounds may be poorly constrained. " + ...
                    "Consider BoundsMethod=gmm_sigma_window or K_override.", reference_csv);
            end
        end

        q  = quantile(ages_y, opts.BoundsQuantiles);
        lo = q(1) - opts.Buffer;
        hi = q(2) + opts.Buffer;
        fprintf("  BoundsMethod = quantile | quantiles=[%.2f %.2f], buffer=%.1f Ma\n", ...
            opts.BoundsQuantiles(1), opts.BoundsQuantiles(2), opts.Buffer);
end

% Enforce minimum window width
% GMM sigma window: floor is 2 Ma (numerical guard only)
% quantile: floor is 10 Ma (guards against sparse-assignment edge cases)
if (hi - lo) < minWidth
    mid = (hi + lo) / 2;
    lo  = mid - minWidth / 2;
    hi  = mid + minWidth / 2;
    fprintf("  NOTE: window expanded to minimum width of %.0f Ma\n", minWidth);
end

Ty = [lo, hi];

fprintf("  Target-component window: %.1f - %.1f Ma | candidate model-start age = %.1f Ma (mu+%.0fsigma; component sigma=%.1f Ma)\n", ...
    Ty(1), Ty(2), model_start_Ma, opts.NSigma, sig_y);

% ---- KDE for QA (original ages only, not jittered) ----
age_min = max(0, prctile(ages, 0.5) - 50);
age_max = prctile(ages, 99.5) + 50;
x_kde   = linspace(age_min, age_max, 1200);

if isnan(opts.Bandwidth)
    [f_kde, ~] = ksdensity(ages, x_kde);
else
    [f_kde, ~] = ksdensity(ages, x_kde, "Bandwidth", opts.Bandwidth);
end

gmm_pdf    = pdf(gm, x_kde(:));
gmm_pdf    = gmm_pdf    / trapz(x_kde, gmm_pdf);
f_kde_norm = f_kde      / trapz(x_kde, f_kde);

% ---- QA figure: two panels (BIC curve + age distributions) ----
fig = figure("Visible","off","Position",[100 100 900 650]);

% Panel 1: BIC curve (skip if K was manually overridden)
ax1 = subplot(2,1,1);
if opts.K_override > 0
    text(0.5, 0.5, sprintf("K manually set to %d (BIC not run)", K_use), ...
        "HorizontalAlignment","center","Units","normalized","FontSize",11);
    axis off
else
    plot(ax1, K_range, bic_vals, "o-k", "LineWidth", 1.2, "MarkerFaceColor","k");
    hold(ax1,"on");
    plot(ax1, K_use, bic_vals(K_range==K_use), "ro", ...
        "MarkerSize", 9, "LineWidth", 1.5);
    xlabel(ax1, "Number of components K");
    ylabel(ax1, "BIC");
    title(ax1, sprintf("BIC model selection — best K = %d", K_use));
    grid(ax1, "on");
    legend(ax1, "BIC", "Selected K", "Location","best");
end

% Panel 2: KDE + GMM pdf + target window + model start time
ax2 = subplot(2,1,2);
plot(ax2, x_kde, f_kde_norm, "LineWidth", 1.4); hold(ax2,"on");
plot(ax2, x_kde, gmm_pdf,    "--",  "LineWidth", 1.4);
yl = ylim(ax2);
patch(ax2, [Ty(1) Ty(2) Ty(2) Ty(1)], [yl(1) yl(1) yl(2) yl(2)], ...
    0.85*[1 1 1], "FaceAlpha", 0.35, "EdgeColor","none");
plot(ax2, [mu_y mu_y], yl, ":", "LineWidth", 1.3, "Color", [0.2 0.2 0.2]);
plot(ax2, [model_start_Ma model_start_Ma], yl, "-", ...
    "LineWidth", 1.5, "Color", [0.8 0.2 0.2]);
xlabel(ax2, sprintf("%s age (Ma)", reference_system));
ylabel(ax2, "Normalized density");
override_note = "";
if opts.K_override > 0
    override_note = " [K fixed]";
end
title(ax2, sprintf("Target window: %.1f-%.1f Ma | mu=%.1f, sig=%.1f, N-sigma=%.1f | candidate start=%.1f Ma | K=%d%s", ...
    Ty(1), Ty(2), mu_y, sig_y, opts.NSigma, model_start_Ma, K_use, override_note));
legend(ax2, "KDE (raw ages)", ...
    sprintf("GMM (K=%d%s)", K_use, override_note), ...
    sprintf("Target-component window [%s]", normalize_bounds_name(opts.BoundsMethod)), ...
    "Component mean (mu)", ...
    sprintf("Candidate model-start age (mu+%.0fsig = %.1f Ma)", opts.NSigma, model_start_Ma), ...
    "Location","best");
grid(ax2, "on");

[~, base] = fileparts(reference_csv);
pngpath = fullfile(outdir, "target_component_plot.png");
saveas(fig, pngpath);
close(fig);

% ---- Output struct ----
out = struct();
out.reference_system = reference_system;
out.NSigma_used         = opts.NSigma;
out.bounds_method       = normalize_bounds_name(opts.BoundsMethod);
out.K_used          = K_use;
out.K_bic           = K_bic;
out.K_selection_method = K_selection_method;
out.K_override_applied = K_override_applied;
out.Nages           = N;
out.Nselected       = nnz(p_y >= opts.MemberProbThresh);
out.QAplot          = pngpath;
out.used_errors     = true;
out.Nmc_used        = Nmc;
out.target_component_age_range = target_age_range;
out.target_window_Ma = Ty;
out.target_component_mean_Ma = mu_y;
out.target_component_sigma_Ma = sig_y;
out.target_component_weight = w_y;
out.candidate_model_start_Ma = model_start_Ma;

% Neutral public summary. The candidate start age is reported without
% treating the component sigma as its uncertainty.
out_tbl = table(string(base), reference_system, Ty(1), Ty(2), mu_y, sig_y, w_y, ...
    model_start_Ma, opts.NSigma, string(normalize_bounds_name(opts.BoundsMethod)), ...
    target_age_range(1), target_age_range(2), ...
    K_use, K_bic, K_selection_method, K_override_applied, ...
    N, out.Nselected, true, Nmc, ...
    'VariableNames', {'Dataset','ReferenceSystem','target_window_lo_Ma','target_window_hi_Ma', ...
    'target_component_mean_Ma','target_component_sigma_Ma','target_component_weight', ...
    'candidate_model_start_Ma','NSigma','BoundsMethod', ...
    'target_component_search_lo_Ma','target_component_search_hi_Ma', ...
    'K_used','K_bic','K_selection_method','K_override_applied', ...
    'Nages','Nselected','used_errors','Nmc_per_grain'});
writetable(out_tbl, fullfile(outdir, "target_component_summary.csv"));

end

% -----------------------------------------------------------------------
function name = normalize_bounds_name(name)
name = string(name);
end

% -----------------------------------------------------------------------
function values = to_numeric(values)
if isnumeric(values)
    values = double(values(:));
else
    values = str2double(string(values(:)));
end
end

% -----------------------------------------------------------------------
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

% -----------------------------------------------------------------------
function gm = fit_gmm(z, K)
% Fit a K-component GMM with regularization and multiple random starts.
% Use a K-specific seed so the same K gives the same fit whether it is
% reached through the BIC sweep or requested as a manual override.
rng(1000 + round(K), "twister");
gm = fitgmdist(z, K, ...
    "RegularizationValue", 1e-3, ...
    "Replicates",          15, ...
    "Options",             statset("MaxIter", 500));
end
