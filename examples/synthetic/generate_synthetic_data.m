function generate_synthetic_data(outdir)
% GENERATE_SYNTHETIC_DATA
% Creates three synthetic sample datasets for testing DetritalChronFilter.
%
% Each sample has a known target component (defined by mu and
% sigma) and a mix of analyses designed to exercise all screening paths:
%   - eligible_after_reference_screen
%   - older_than_reference
%   - short_pair_interval
%   - age_order_unresolved
%
% With the default gmm_sigma_window BoundsMethod and NSigma=1.0, the expected
% target-component windows (mu +/- 1*sigma) are:
%   SampleA: ~85–95 Ma  (mu=90, sigma=5)
%   SampleB: ~69–81 Ma  (mu=75, sigma=6)
%   SampleC: ~93–107 Ma (mu=100, sigma=7)
%
% Usage:
%   generate_synthetic_data              % writes to inputs_generated/
%   generate_synthetic_data("my_dir")    % writes to my_dir/

if nargin < 1 || isempty(outdir)
    example_dir = string(fileparts(mfilename("fullpath")));
    outdir = fullfile(example_dir, "inputs_generated");
end

rng(42);  % reproducible

% ---- Define three synthetic samples ----
samples = struct();

% Sample A: clear target component ~90 Ma, well-separated from older components
samples(1).name         = "SampleA";
samples(1).target_mu    = 90;    % Ma — target-component mean
samples(1).target_sig   = 5;     % Ma — component width
samples(1).older_mus    = [140, 200];  % older zircon U-Pb components
samples(1).older_sigs   = [10, 15];
samples(1).older_weights = [0.35, 0.25];  % fraction of reference ages
samples(1).young_weight  = 0.40;

% Sample B: target component ~75 Ma, closer to an older component
samples(2).name         = "SampleB";
samples(2).target_mu    = 75;
samples(2).target_sig   = 6;
samples(2).older_mus    = [110, 170];
samples(2).older_sigs   = [8, 20];
samples(2).older_weights = [0.40, 0.20];
samples(2).young_weight  = 0.40;

% Sample C: target component ~100 Ma, larger uncertainties
samples(3).name         = "SampleC";
samples(3).target_mu    = 100;
samples(3).target_sig   = 7;
samples(3).older_mus    = [150, 220];
samples(3).older_sigs   = [12, 18];
samples(3).older_weights = [0.30, 0.30];
samples(3).young_weight  = 0.40;

% ---- Generate data for each sample ----
for ci = 1:numel(samples)
    c    = samples(ci);
    cdir = fullfile(outdir, c.name);
    if ~isfolder(cdir), mkdir(cdir); end

    Ty_hi = c.target_mu + 2*c.target_sig;

    % -- ReferenceDistribution.csv: full distribution for GMM inference --
    n_reference = 120;
    weights = [c.young_weight, c.older_weights];
    ns      = round(weights / sum(weights) * n_reference);
    mus     = [c.target_mu, c.older_mus];
    sigs    = [c.target_sig, c.older_sigs];

    reference_ages = [];
    reference_errors = [];
    for k = 1:numel(mus)
        ages_k = mus(k) + sigs(k) .* randn(ns(k), 1);
        errs_k = max(0.5, sigs(k)*0.15 + 0.25*randn(ns(k),1));
        reference_ages = [reference_ages; ages_k]; %#ok<AGROW>
        reference_errors = [reference_errors; errs_k]; %#ok<AGROW>
    end
    % Ensure positive ages
    reference_ages = max(10, reference_ages);
    reference_errors = max(0.5, reference_errors);

    reference_ids = compose("Ref%03d", (1:numel(reference_ages))');
    T_reference = table(repmat("ZrnUPb", numel(reference_ages), 1), ...
        reference_ids, reference_ages, reference_errors, ...
        'VariableNames', {'ReferenceSystem','GrainID','Age_Ma','Age_1sigma_Ma'});
    writetable(T_reference, fullfile(cdir, "ReferenceDistribution.csv"));

    % -- Helper: make cooling ages relative to a crystallization age --
    % longer-interval grains: cool 20-60 Myr after the selected component
    % older-than-reference grains: cooling ages older than Ty_hi
    % short-interval grains: cool within about 5 Myr of crystallization

    N_ap  = 100;
    N_hbl = 100;
    N_zrn = 90;

    % ---- Paired apatite example rows ----
    % ~80% longer interval, ~10% older than reference, ~10% short interval
    n_keep = round(0.80 * N_ap);
    n_leg  = round(0.10 * N_ap);
    n_mag  = N_ap - n_keep - n_leg;

    % Crystallization ages: drawn from the target component
    tc_keep = c.target_mu + c.target_sig .* randn(n_keep, 1);
    tc_leg  = c.target_mu + c.target_sig .* randn(n_leg,  1);
    tc_mag  = c.target_mu + c.target_sig .* randn(n_mag,  1);

    % Cooling ages
    th_keep = tc_keep - (20 + 40*rand(n_keep,1));   % 20-60 Myr lag
    th_leg  = Ty_hi + 5 + 20*rand(n_leg,  1);       % older than target window
    th_mag  = tc_mag  - (2 + 10*rand(n_mag, 1));   % 2-12 Myr lag — straddles Delta=8

    th_all = [th_keep; th_leg; th_mag];
    tc_all = [tc_keep; tc_leg; tc_mag];

    % Add a few grains with nominally reversed age order that overlap within
    % combined 2-sigma uncertainty so AOU appears in the example outputs.
    % These are created by making the He age slightly older than the U-Pb age
    % with realistic uncertainties so the discordance is within ~1.5 sigma.
    n_fd   = 4;
    tc_fd  = c.target_mu + c.target_sig .* randn(n_fd, 1);
    th_fd  = tc_fd + (0.5 + 1.5*rand(n_fd,1));  % He slightly OLDER than U-Pb (reversed)
    th_all = [th_all; th_fd];
    tc_all = [tc_all; tc_fd];
    n_tot  = numel(th_all);

    ap_grains   = arrayfun(@(i) sprintf("Ap%03d", i), (1:n_tot)', 'UniformOutput', false);
    ap_he_err   = max(0.5, 0.08*th_all + randn(n_tot,1)) ./ 2;
    ap_pb_err   = max(1.0, 0.05*tc_all + randn(n_tot,1)) ./ 2;

    T_ap = table(string(ap_grains), th_all, ap_he_err, tc_all, ap_pb_err, ...
        'VariableNames', {'ApGrain','ApHeDate','ApHe1sigerr','ApPbDate','ApPb1sigerr'});

    % ---- Unpaired hornblende example rows ----
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

    % ---- Paired zircon example rows ----
    % ~75% longer interval, ~15% short interval, ~10% older than reference
    n_keep_z = round(0.75 * N_zrn);
    n_mag_z  = round(0.15 * N_zrn);
    n_leg_z  = N_zrn - n_keep_z - n_mag_z;

    tc_zkeep = c.target_mu + c.target_sig .* randn(n_keep_z, 1);
    tc_zmag  = c.target_mu + c.target_sig .* randn(n_mag_z,  1);
    tc_zleg  = c.target_mu + c.target_sig .* randn(n_leg_z,  1);

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
    % ---- Flexible long-format ChronometerData.csv ----
    ap_pair_ids = "ApPair:" + string(T_ap.ApGrain);
    zrn_pair_ids = "ZrnPair:" + string(T_zrn.ZrnGrain);
    T_chronometers = [ ...
        make_rows("HblAr", T_hbl.HblGrain, T_hbl.HblArDate, ...
            T_hbl.HblAr1sigerr, repmat("", height(T_hbl), 1), ...
            repmat("", height(T_hbl), 1), true(height(T_hbl), 1));
        make_rows("ApHe", T_ap.ApGrain, T_ap.ApHeDate, ...
            T_ap.ApHe1sigerr, ap_pair_ids, ...
            repmat("expected_younger", height(T_ap), 1), true(height(T_ap), 1));
        make_rows("ApUPb", T_ap.ApGrain, T_ap.ApPbDate, ...
            T_ap.ApPb1sigerr, ap_pair_ids, ...
            repmat("expected_older", height(T_ap), 1), true(height(T_ap), 1));
        make_rows("ZrnHe", T_zrn.ZrnGrain, T_zrn.ZrnHeDate, ...
            T_zrn.ZrnHe1sigerr, zrn_pair_ids, ...
            repmat("expected_younger", height(T_zrn), 1), true(height(T_zrn), 1));
        make_rows("ZrnUPb", T_zrn.ZrnGrain, T_zrn.ZrnPbDate, ...
            T_zrn.ZrnPb1sigerr, zrn_pair_ids, ...
            repmat("expected_older", height(T_zrn), 1), false(height(T_zrn), 1))];
    writetable(T_chronometers, fullfile(cdir, "ChronometerData.csv"));

    fprintf("Generated synthetic data for %s  (component mu=%.0f Ma, sigma=%.0f Ma, sigma window: ~%.0f-%.0f Ma)\n", ...
        c.name, c.target_mu, c.target_sig, c.target_mu - c.target_sig, c.target_mu + c.target_sig);
end

fprintf("\nSynthetic sample data written to: %s\n", outdir);
fprintf("Run  run_example  to process with the full pipeline.\n");

end

function T = make_rows(chronometer, grain_id, age, uncertainty, ...
        pair_id, pair_role, use_for_model)
n = numel(age);
T = table(repmat(string(chronometer), n, 1), string(grain_id), ...
    age(:), uncertainty(:), string(pair_id), string(pair_role), ...
    logical(use_for_model), ...
    'VariableNames', {'Chronometer','GrainID','Age_Ma','Age_1sigma_Ma', ...
    'PairID','PairRole','UseForModel'});
end
