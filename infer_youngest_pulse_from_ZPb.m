function out = infer_youngest_pulse_from_ZPb(zpb_csv, outdir, opts)
% INFER_YOUNGEST_PULSE_FROM_ZPB
% Infers youngest magmatic pulse window from detrital zircon U-Pb ages
% using a Gaussian Mixture Model (GMM) fit to the age distribution.
%
% The youngest GMM component (by mean) is identified as the youngest
% magmatic pulse. Its mean (mu) and standard deviation (sigma) are used
% to define the pulse window and a recommended model start time.
%
% If ZrnPb2sigerr column exists, uses Monte Carlo jittering to propagate
% 2-sigma analytical uncertainties into GMM fitting. Posterior assignment
% uses original measured ages for interpretability.
%
% K selection: BIC over K = 2..Kmax by default. Use K_override to fix K
% after inspecting QA plots.
%
% -----------------------------------------------------------------------
% BOUNDS METHODS  (BoundsMethod parameter)
% -----------------------------------------------------------------------
%   "gmm_ci" (default)
%       Window = [mu - NSigma*sigma,  mu + NSigma*sigma]
%       Uses the GMM component's own mean and std to define the window.
%       Robust to outlier grains at the tails of the assigned population.
%       NSigma=1.0 (default) gives ~68% of the component distribution;
%       NSigma=2.0 gives ~95%. The upper bound (mu + NSigma*sigma) also
%       serves as a principled model start time for QTQt/Pecube — it
%       represents the mean crystallization age of the youngest pulse
%       plus a tolerance for within-pluton age variability.
%       Minimum window floor: 2 Ma (numerical guard only; the GMM sigma
%       already encodes real pulse width, so the floor is rarely active).
%
%   "quantile"
%       Window = quantile(assigned_ages, BoundsQuantiles) + Buffer.
%       Original method. More sensitive to outlier grains at the tails.
%       Use if the GMM sigma is poorly constrained (small N in component).
%
% -----------------------------------------------------------------------
% REQUIRED / OPTIONAL COLUMNS
% -----------------------------------------------------------------------
%   Required : ZrnPbDate
%   Optional : ZrnPb2sigerr  (2-sigma errors; used for MC jittering)
%
% -----------------------------------------------------------------------
% OUTPUT STRUCT
% -----------------------------------------------------------------------
%   out.Tyoung            [lo hi] Ma — pulse window passed to filter
%   out.mu_young          GMM component mean (Ma)
%   out.sigma_young       GMM component std (Ma)
%   out.weight_young      GMM component mixture weight
%   out.model_start_Ma      mu + NSigma*sigma — recommended QTQt/Pecube
%                           start time
%   out.model_start_err_Ma  sigma_young — use as +/- for QTQt bounding box;
%                           represents real within-pluton age spread, not
%                           just uncertainty on mu
%   out.NSigma_used       NSigma value used to compute model_start_Ma
%   out.bounds_method     BoundsMethod used ("gmm_ci" or "quantile")
%   out.K_used            K used for final GMM
%   out.K_bic             K selected by BIC (= K_used unless overridden)
%   out.Nages             number of valid input ages
%   out.Nselected         grains assigned to youngest component
%   out.QAplot            path to saved QA figure
%   out.used_errors       logical — were 2-sigma errors incorporated?
%   out.Nmc_used          Monte Carlo draws used (0 if no errors)
%
% -----------------------------------------------------------------------
% USAGE
% -----------------------------------------------------------------------
%   out = infer_youngest_pulse_from_ZPb("TC/ZrnPb.csv", "TC/QA")
%   out = infer_youngest_pulse_from_ZPb("TC/ZrnPb.csv", "TC/QA", NSigma=1.5)
%   out = infer_youngest_pulse_from_ZPb("TC/ZrnPb.csv", "TC/QA", K_override=3)
%   out = infer_youngest_pulse_from_ZPb("TC/ZrnPb.csv", "TC/QA", BoundsMethod="quantile")

arguments
    zpb_csv  (1,1) string
    outdir   (1,1) string  = "ZPb_pulse_QA"
    opts.Kmax             (1,1) double = 6
    opts.K_override       (1,1) double = 0
    opts.MinWeight        (1,1) double = 0.05
    opts.MemberProbThresh (1,1) double = 0.7
    opts.BoundsMethod     (1,1) string ...
                          {mustBeMember(opts.BoundsMethod, ...
                          ["gmm_ci","quantile"])} = "gmm_ci"
    opts.NSigma           (1,1) double = 1.0   % used by gmm_ci and for model_start_Ma
    opts.BoundsQuantiles  (1,2) double = [0.02 0.98]  % used by quantile method only
    opts.MinPulseWidth    (1,1) double = NaN   % Ma; NaN = auto (2 Ma for gmm_ci, 10 Ma for quantile)
    opts.Buffer           (1,1) double = 2     % Ma; added to each bound (quantile only)
    opts.Bandwidth        (1,1) double = NaN   % KDE bandwidth; NaN = auto
    opts.Nmc              (1,1) double = 50
    opts.MaxTotalSamples  (1,1) double = 200000
end

% Resolve MinPulseWidth: gmm_ci uses 2 Ma (just a numerical floor — the
% sigma already encodes pulse width). quantile uses 10 Ma (guards against
% absurdly narrow empirical windows when few grains are assigned).
if isnan(opts.MinPulseWidth)
    if opts.BoundsMethod == "gmm_ci"
        minWidth = 2;
    else
        minWidth = 10;
    end
else
    minWidth = opts.MinPulseWidth;  % user override
end

if ~isfolder(outdir), mkdir(outdir); end

% ---- Load and validate ----
T = readtable(zpb_csv, "VariableNamingRule","preserve");
vnames = string(T.Properties.VariableNames);

assert(any(vnames == "ZrnPbDate"), ...
    "Column 'ZrnPbDate' not found in %s.\nFound: %s", ...
    zpb_csv, strjoin(vnames, ", "));

ages = T.("ZrnPbDate")(:);
hasErr = any(vnames == "ZrnPb2sigerr");
e2 = hasErr * ones(size(ages));   % placeholder
if hasErr, e2 = T.("ZrnPb2sigerr")(:); end

% Basic cleaning
valid = isfinite(ages) & ages > 0;
ages  = ages(valid);
e2    = e2(valid);

assert(~isempty(ages), "No valid ZrnPbDate ages found in %s", zpb_csv);
N = numel(ages);

% ---- Convert 2σ → 1σ; handle missing/zero errors ----
sig1 = nan(N, 1);
if hasErr
    sig1 = e2 ./ 2;
    sig1(~isfinite(sig1) | sig1 <= 0) = NaN;
end

useErrors = hasErr && any(isfinite(sig1));

% Fallback sigma: trimmed mean of valid sigmas (more robust than median
% when a few grains have anomalously large reported uncertainties)
if useErrors
    s_valid = sig1(isfinite(sig1));
    trim_lo = prctile(s_valid, 10);
    trim_hi = prctile(s_valid, 90);
    sig_fallback = mean(s_valid(s_valid >= trim_lo & s_valid <= trim_hi));
    if ~isfinite(sig_fallback) || sig_fallback <= 0
        sig_fallback = 2.0; % Ma hard floor
    end
    sig_use = sig1;
    sig_use(~isfinite(sig_use)) = sig_fallback;
end

% ---- Monte Carlo jittering to propagate uncertainties into GMM ----
if useErrors
    Nmc = max(1, round(opts.Nmc));
    total = N * Nmc;
    if total > opts.MaxTotalSamples
        Nmc = max(1, floor(opts.MaxTotalSamples / N));
    end
    rng(0);
    A = repmat(ages, 1, Nmc) + repmat(sig_use, 1, Nmc) .* randn(N, Nmc);
    z = A(:);
else
    Nmc = 0;
    z = ages(:);
end

z = z(isfinite(z) & z > 0);

% ---- BIC-based K selection ----
Kmin = 2;
Kmax = max(Kmin, opts.Kmax);

if opts.K_override > 0
    % Manual override: fit only the requested K
    K_bic  = opts.K_override;   % record that BIC wasn't used
    K_use  = opts.K_override;
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
    gm    = gm_store{best_ki};
end

% ---- Extract youngest valid component ----
mu  = gm.mu(:);
sig = sqrt(squeeze(gm.Sigma));
w   = gm.ComponentProportion(:);

validComp   = w >= opts.MinWeight;
mu_masked   = mu;
mu_masked(~validComp) = Inf;
[mu_y, idx_y] = min(mu_masked);
if isinf(mu_y)
    [mu_y, idx_y] = min(mu);   % fallback: all components below MinWeight
end
sig_y = sig(idx_y);
w_y   = w(idx_y);

% Posterior membership on original ages — used for Nselected reporting
% regardless of BoundsMethod, and for the quantile method's grain selection.
post = posterior(gm, ages(:));
p_y  = post(:, idx_y);

% ---- Define pulse window ----
% model_start_Ma is always mu + NSigma*sigma regardless of BoundsMethod.
% This is the recommended start time for QTQt/Pecube: it represents the
% GMM-estimated mean crystallization age of the youngest pulse, plus a
% tolerance for within-pluton age variability. NSigma=1 means ~1 standard
% deviation above the mean — i.e. the upper edge of the core of the pulse
% distribution, not the far tail.
model_start_Ma = mu_y + opts.NSigma * sig_y;

switch opts.BoundsMethod

    case "gmm_ci"
        % Window defined symmetrically from GMM component parameters.
        % Not sensitive to outlier grains at the tails.
        lo = mu_y - opts.NSigma * sig_y;
        hi = mu_y + opts.NSigma * sig_y;
        fprintf("  BoundsMethod = gmm_ci | mu=%.1f, sigma=%.1f, NSigma=%.1f\n", ...
            mu_y, sig_y, opts.NSigma);

    case "quantile"
        % Original method: empirical quantiles of posterior-assigned grains.
        % Requires enough grains assigned to youngest component.
        sel = p_y >= opts.MemberProbThresh;

        if nnz(sel) >= 10
            ages_y = ages(sel);
        else
            ages_y = ages(ages >= (mu_y - 2*sig_y) & ages <= (mu_y + 2*sig_y));
            if numel(ages_y) < 10
                ages_y = ages;
                warning("infer_youngest_pulse: few grains assigned to youngest " + ...
                    "component — pulse bounds may be poorly constrained. " + ...
                    "Consider BoundsMethod=gmm_ci or K_override.", zpb_csv);
            end
        end

        q  = quantile(ages_y, opts.BoundsQuantiles);
        lo = q(1) - opts.Buffer;
        hi = q(2) + opts.Buffer;
        fprintf("  BoundsMethod = quantile | quantiles=[%.2f %.2f], buffer=%.1f Ma\n", ...
            opts.BoundsQuantiles(1), opts.BoundsQuantiles(2), opts.Buffer);
end

% Enforce minimum window width
% gmm_ci: floor is 2 Ma (numerical guard only; sigma encodes real pulse width)
% quantile: floor is 10 Ma (guards against sparse-assignment edge cases)
if (hi - lo) < minWidth
    mid = (hi + lo) / 2;
    lo  = mid - minWidth / 2;
    hi  = mid + minWidth / 2;
    fprintf("  NOTE: window expanded to minimum width of %.0f Ma\n", minWidth);
end

Ty = [lo, hi];

fprintf("  Pulse window: %.1f - %.1f Ma | model_start = %.1f +/- %.1f Ma (mu+%.0fsig, err=sigma)\n", ...
    Ty(1), Ty(2), model_start_Ma, sig_y, opts.NSigma);

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

% Panel 2: KDE + GMM pdf + pulse window + model start time
ax2 = subplot(2,1,2);
plot(ax2, x_kde, f_kde_norm, "LineWidth", 1.4); hold(ax2,"on");
plot(ax2, x_kde, gmm_pdf,    "--",  "LineWidth", 1.4);
yl = ylim(ax2);
patch(ax2, [Ty(1) Ty(2) Ty(2) Ty(1)], [yl(1) yl(1) yl(2) yl(2)], ...
    0.85*[1 1 1], "FaceAlpha", 0.35, "EdgeColor","none");
plot(ax2, [mu_y mu_y], yl, ":", "LineWidth", 1.3, "Color", [0.2 0.2 0.2]);
plot(ax2, [model_start_Ma model_start_Ma], yl, "-", ...
    "LineWidth", 1.5, "Color", [0.8 0.2 0.2]);
xlabel(ax2, "Zircon U-Pb age (Ma)");
ylabel(ax2, "Normalized density");
override_note = "";
if opts.K_override > 0
    override_note = " [K fixed]";
end
title(ax2, sprintf("Window: %.1f-%.1f Ma | mu=%.1f, sig=%.1f, N-sigma=%.1f | model start=%.1f Ma | K=%d%s", ...
    Ty(1), Ty(2), mu_y, sig_y, opts.NSigma, model_start_Ma, K_use, override_note));
legend(ax2, "KDE (raw ages)", ...
    sprintf("GMM (K=%d%s)", K_use, override_note), ...
    sprintf("Pulse window [%s]", opts.BoundsMethod), ...
    "Component mean (mu)", ...
    sprintf("Model start (mu+%.0fsig = %.1f Ma)", opts.NSigma, model_start_Ma), ...
    "Location","best");
grid(ax2, "on");

[~, base] = fileparts(zpb_csv);
pngpath = fullfile(outdir, base + "_youngest_pulse_QA.png");
saveas(fig, pngpath);
close(fig);

% ---- Output struct ----
out = struct();
out.Tyoung          = Ty;
out.mu_young        = mu_y;
out.sigma_young     = sig_y;
out.weight_young    = w_y;
out.model_start_Ma      = model_start_Ma;   % mu + NSigma*sigma — use as QTQt/Pecube start
out.model_start_err_Ma  = sig_y;            % sigma_young — use as +/- for QTQt bounding box
out.NSigma_used         = opts.NSigma;
out.bounds_method       = opts.BoundsMethod;
out.K_used          = K_use;
out.K_bic           = K_bic;
out.Nages           = N;
out.Nselected       = nnz(p_y >= opts.MemberProbThresh);
out.QAplot          = pngpath;
out.used_errors     = useErrors;
out.Nmc_used        = Nmc;

% Summary CSV — includes model_start_Ma so every catchment run produces
% a ready-to-use recommended start time alongside the pulse window.
out_tbl = table(string(base), Ty(1), Ty(2), mu_y, sig_y, w_y, ...
    model_start_Ma, sig_y, opts.NSigma, string(opts.BoundsMethod), ...
    K_use, K_bic, N, out.Nselected, useErrors, Nmc, ...
    'VariableNames', {'Dataset','Tyoung_lo','Tyoung_hi', ...
    'mu_young','sigma_young','weight_young', ...
    'model_start_Ma','model_start_err_Ma','NSigma','BoundsMethod', ...
    'K_used','K_bic','Nages','Nselected','used_errors','Nmc_per_grain'});
writetable(out_tbl, fullfile(outdir, base + "_youngest_pulse_summary.csv"));

end

% -----------------------------------------------------------------------
function gm = fit_gmm(z, K)
% Fit a K-component GMM with regularization and multiple random starts.
gm = fitgmdist(z, K, ...
    "RegularizationValue", 1e-3, ...
    "Replicates",          15, ...
    "Options",             statset("MaxIter", 500));
end
