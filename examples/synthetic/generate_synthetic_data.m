function generate_synthetic_data(outdir)
% GENERATE_SYNTHETIC_DATA
% Creates three synthetic catchment datasets for testing DetritalChronFilter.
%
% Each catchment has a known youngest target component (defined by mu and
% sigma) and a mix of grains designed to exercise all screening paths:
%   - eligible_after_reference_screen
%   - older_than_reference
%   - short_crystallization_cooling_interval
%   - age_order_unresolved
%
% With the default gmm_sigma_window BoundsMethod and NSigma=1.0, the expected
% target-component windows (mu +/- 1*sigma) are:
%   CatchmentA: ~85–95 Ma  (mu=90, sigma=5)
%   CatchmentB: ~69–81 Ma  (mu=75, sigma=6)
%   CatchmentC: ~93–107 Ma (mu=100, sigma=7)
%
% Usage:
%   generate_synthetic_data              % writes to inputs_generated/
%   generate_synthetic_data("my_dir")    % writes to my_dir/

if nargin < 1 || isempty(outdir)
    example_dir = string(fileparts(mfilename("fullpath")));
    outdir = fullfile(example_dir, "inputs_generated");
end

rng(42);  % reproducible

% ---- Define three synthetic catchments ----
catchments = struct();

% Catchment A: clear youngest pulse ~90 Ma, well-separated from older populations
catchments(1).name         = "CatchmentA";
catchments(1).pulse_mu     = 90;    % Ma — youngest pulse mean
catchments(1).pulse_sig    = 5;     % Ma — pulse width
catchments(1).older_mus    = [140, 200];  % older zircon U-Pb components
catchments(1).older_sigs   = [10, 15];
catchments(1).older_weights = [0.35, 0.25];  % fraction of ZPb grains
catchments(1).young_weight  = 0.40;

% Catchment B: youngest pulse ~75 Ma, closer to older population
catchments(2).name         = "CatchmentB";
catchments(2).pulse_mu     = 75;
catchments(2).pulse_sig    = 6;
catchments(2).older_mus    = [110, 170];
catchments(2).older_sigs   = [8, 20];
catchments(2).older_weights = [0.40, 0.20];
catchments(2).young_weight  = 0.40;

% Catchment C: youngest pulse ~100 Ma, larger uncertainties
catchments(3).name         = "CatchmentC";
catchments(3).pulse_mu     = 100;
catchments(3).pulse_sig    = 7;
catchments(3).older_mus    = [150, 220];
catchments(3).older_sigs   = [12, 18];
catchments(3).older_weights = [0.30, 0.30];
catchments(3).young_weight  = 0.40;

% ---- Generate data for each catchment ----
for ci = 1:numel(catchments)
    c    = catchments(ci);
    cdir = fullfile(outdir, c.name);
    if ~isfolder(cdir), mkdir(cdir); end

    Ty_lo = c.pulse_mu - 2*c.pulse_sig;
    Ty_hi = c.pulse_mu + 2*c.pulse_sig;

    % -- ZrnPb.csv: detrital zircon U-Pb (used for pulse inference) --
    N_zpb = 120;
    weights = [c.young_weight, c.older_weights];
    ns      = round(weights / sum(weights) * N_zpb);
    mus     = [c.pulse_mu, c.older_mus];
    sigs    = [c.pulse_sig, c.older_sigs];

    zpb_ages = [];
    zpb_errs = [];
    for k = 1:numel(mus)
        ages_k = mus(k) + sigs(k) .* randn(ns(k), 1);
        errs_k = max(1, sigs(k)*0.3 + 0.5*randn(ns(k),1));  % generated as 2σ for continuity
        zpb_ages = [zpb_ages; ages_k]; %#ok<AGROW>
        zpb_errs = [zpb_errs; errs_k]; %#ok<AGROW>
    end
    % Ensure positive ages
    zpb_ages = max(10, zpb_ages);
    zpb_errs = max(0.5, zpb_errs);

    T_zpb = table(zpb_ages, zpb_errs ./ 2, ...
        'VariableNames', {'ZrnPbDate','ZrnPb1sigerr'});
    writetable(T_zpb, fullfile(cdir, "ZrnPb.csv"));

    % -- Helper: make cooling ages relative to a crystallization age --
    % longer-interval grains: cool 20-60 Myr after the selected component
    % older-than-reference grains: cooling ages older than Ty_hi
    % short-interval grains: cool within about 5 Myr of crystallization

    N_ap  = 100;
    N_hbl = 100;
    N_zrn = 90;

    % ---- ApHeApPb.csv ----
    % ~80% longer interval, ~10% older than reference, ~10% short interval
    n_keep = round(0.80 * N_ap);
    n_leg  = round(0.10 * N_ap);
    n_mag  = N_ap - n_keep - n_leg;

    % Crystallization ages: drawn from pulse
    tc_keep = c.pulse_mu + c.pulse_sig .* randn(n_keep, 1);
    tc_leg  = c.pulse_mu + c.pulse_sig .* randn(n_leg,  1);
    tc_mag  = c.pulse_mu + c.pulse_sig .* randn(n_mag,  1);

    % Cooling ages
    th_keep = tc_keep - (20 + 40*rand(n_keep,1));   % 20-60 Myr lag
    th_leg  = Ty_hi + 5 + 20*rand(n_leg,  1);       % older than pulse window
    th_mag  = tc_mag  - (2 + 10*rand(n_mag, 1));   % 2-12 Myr lag — straddles Delta=8

    th_all = [th_keep; th_leg; th_mag];
    tc_all = [tc_keep; tc_leg; tc_mag];

    % Add a few grains with nominally reversed age order that overlap within
    % combined 2-sigma uncertainty so AOU appears in the example outputs.
    % These are created by making the He age slightly older than the U-Pb age
    % with realistic uncertainties so the discordance is within ~1.5 sigma.
    n_fd   = 4;
    tc_fd  = c.pulse_mu + c.pulse_sig .* randn(n_fd, 1);
    th_fd  = tc_fd + (0.5 + 1.5*rand(n_fd,1));  % He slightly OLDER than U-Pb (reversed)
    th_all = [th_all; th_fd];
    tc_all = [tc_all; tc_fd];
    n_tot  = numel(th_all);

    ap_grains   = arrayfun(@(i) sprintf("Ap%03d", i), (1:n_tot)', 'UniformOutput', false);
    ap_he_err   = max(0.5, 0.08*th_all + randn(n_tot,1)) ./ 2;
    ap_pb_err   = max(1.0, 0.05*tc_all + randn(n_tot,1)) ./ 2;

    T_ap = table(string(ap_grains), th_all, ap_he_err, tc_all, ap_pb_err, ...
        'VariableNames', {'ApGrain','ApHeDate','ApHe1sigerr','ApPbDate','ApPb1sigerr'});
    writetable(T_ap, fullfile(cdir, "ApHeApPb.csv"));

    % ---- HblAr.csv ----
    % Hornblende ages span the reference: about 60% younger than the older
    % boundary and about 40% older than it.
    n_keep_hbl = round(0.60 * N_hbl);
    n_leg_hbl  = N_hbl - n_keep_hbl;

    % Hbl cooling ages include values near or younger than the upper bound.
    th_hbl_keep = (Ty_hi - 5) - 15*rand(n_keep_hbl, 1);   % just below Ty_hi
    th_hbl_leg  = Ty_hi + 2  + 30*rand(n_leg_hbl, 1);      % above Ty_hi

    th_hbl = [th_hbl_keep; th_hbl_leg];
    n_hbl  = numel(th_hbl);

    hbl_grains  = arrayfun(@(i) sprintf("Hbl%03d", i), (1:n_hbl)', 'UniformOutput', false);
    hbl_err_1sig = max(0.5, 0.03*th_hbl + 0.3*randn(n_hbl,1));  % ~3% 1sig

    T_hbl = table(string(hbl_grains), th_hbl, hbl_err_1sig, ...
        'VariableNames', {'HblGrain','HblArDate','HblAr1sigerr'});
    writetable(T_hbl, fullfile(cdir, "HblAr.csv"));

    % ---- ZrnHeZrnPb.csv ----
    % ~75% longer interval, ~15% short interval, ~10% older than reference
    n_keep_z = round(0.75 * N_zrn);
    n_mag_z  = round(0.15 * N_zrn);
    n_leg_z  = N_zrn - n_keep_z - n_mag_z;

    tc_zkeep = c.pulse_mu + c.pulse_sig .* randn(n_keep_z, 1);
    tc_zmag  = c.pulse_mu + c.pulse_sig .* randn(n_mag_z,  1);
    tc_zleg  = c.pulse_mu + c.pulse_sig .* randn(n_leg_z,  1);

    th_zkeep = tc_zkeep - (15 + 35*rand(n_keep_z, 1));   % 15-50 Myr lag
    th_zmag  = tc_zmag  - (2 + 10*rand(n_mag_z, 1));  % 2-12 Myr (straddles Delta=8)
    th_zleg  = Ty_hi + 5 + 20*rand(n_leg_z, 1);          % older than reference

    th_z = [th_zkeep; th_zmag; th_zleg];
    tc_z = [tc_zkeep; tc_zmag; tc_zleg];
    n_z  = numel(th_z);

    zrn_grains  = arrayfun(@(i) sprintf("Zrn%03d", i), (1:n_z)', 'UniformOutput', false);
    zhe_err     = max(0.5, 0.07*th_z + randn(n_z,1)) ./ 2;
    zpb_err_dd  = max(1.0, 0.04*tc_z + randn(n_z,1)) ./ 2;

    T_zrn = table(string(zrn_grains), th_z, zhe_err, tc_z, zpb_err_dd, ...
        'VariableNames', {'ZrnGrain','ZrnHeDate','ZrnHe1sigerr','ZrnPbDate','ZrnPb1sigerr'});
    writetable(T_zrn, fullfile(cdir, "ZrnHeZrnPb.csv"));

    fprintf("Generated synthetic data for %s  (component mu=%.0f Ma, sigma=%.0f Ma, sigma window: ~%.0f-%.0f Ma)\n", ...
        c.name, c.pulse_mu, c.pulse_sig, c.pulse_mu - c.pulse_sig, c.pulse_mu + c.pulse_sig);
end

fprintf("\nSynthetic catchment data written to: %s\n", outdir);
fprintf("Run  run_example  to process with the full pipeline.\n");

end
