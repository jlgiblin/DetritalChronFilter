# DetritalChronFilter

A MATLAB workflow for probabilistic filtering of detrital thermochronologic data to isolate post-magmatic exhumation signals from catchments with complex, multi-pulse magmatic histories.

---

## Overview

In source regions with prolonged or episodic magmatism, detrital thermochronologic age distributions reflect cooling signals from multiple magmatic pulses — not all of which record the exhumation episode of interest. This toolbox identifies and removes grains whose cooling ages most likely predate the youngest magmatic pulse, and flags grains whose cooling ages are indistinguishable from crystallization (potentially recording magmatic thermal relaxation rather than exhumation). The result is a filtered dataset suitable for inverse thermal history modeling targeting the youngest pulse of magmatism.

The workflow:
1. Infers a **catchment-specific youngest magmatic pulse window** from detrital zircon U–Pb age distributions using Gaussian mixture modeling (GMM) with BIC-based component selection
2. Derives a **recommended model start time** (`mu + NSigma*sigma` of the youngest GMM component) for use as a QTQt/Pecube boundary condition
3. Classifies each thermochronologic grain probabilistically relative to that window
4. Produces filtered output CSVs, per-grain reason codes, and a cross-catchment sensitivity table

Supported thermochronologic systems:
- Hornblende ⁴⁰Ar/³⁹Ar
- Apatite (U–Th)/He + U–Pb (double-dated)
- Zircon (U–Th)/He + U–Pb (double-dated)

---

## Requirements

- MATLAB R2020b or later
- Statistics and Machine Learning Toolbox (for `fitgmdist`, `ksdensity`, `normcdf`)

No additional toolboxes or third-party dependencies are required.

---

## Repository structure

```
DetritalChronFilter/
├── run_detrital_pipeline.m          % Main entry point — run this
├── infer_youngest_pulse_from_ZPb.m  % GMM pulse window inference
├── filter_detrital_thermo.m         % Probabilistic grain classification
├── pipeline_cheatsheet.md           % Copy-pasteable command examples
├── example/
│   ├── generate_synthetic_data.m    % Generates toy catchment data
│   └── run_example.m                % End-to-end example run
├── LICENSE
└── README.md
```

---

## Quick start

### 1. Organise your data

Place your input CSVs in the following folder structure. File names must match exactly:

```
Catchments/
├── CatchmentA/
│   ├── ZrnPb.csv        % Detrital zircon U–Pb (used to infer pulse window)
│   ├── HblAr.csv        % Hornblende Ar/Ar
│   ├── ApHeApPb.csv     % Apatite (U–Th)/He + U–Pb
│   └── ZrnHeZrnPb.csv   % Zircon (U–Th)/He + U–Pb
├── CatchmentB/
│   └── ...
```

Each catchment is a subfolder. The four file names are fixed. See [Input format](#input-format) for required column names.

### 2. Run the pipeline

```matlab
run_detrital_pipeline("Catchments", "Output")
```

That's it. The pipeline will:
- Auto-discover all catchment subfolders
- Infer the youngest magmatic pulse window for each catchment using GMM
- Report a recommended model start time (`mu + NSigma*sigma`) for QTQt/Pecube setup
- Filter thermochronologic data under three `P_legacy_mode` options (`standard`, `conservative`, `midpoint`) for sensitivity testing
- Write filtered CSVs, a cross-catchment sensitivity table, and a README explaining all outputs

### 3. Check outputs

```
Output/
├── README.txt                        % Auto-generated guide to outputs and reason codes
├── _pipeline_summary/
│   ├── pipeline_summary.csv          % Pulse window + model start time per catchment
│   └── sensitivity_by_system.csv     % keep_exhumation N across modes — HIGH sensitivity flagged
├── CatchmentA/
│   ├── ZPb_QA/                       % BIC curve + GMM diagnostic plots — inspect before use
│   ├── standard/                     % Filtered outputs using standard P_legacy_mode
│   │   ├── kept_strict.csv           *** USE THIS for thermal history modeling ***
│   │   ├── kept_plus_flagged.csv
│   │   ├── all_data_classified.csv
│   │   ├── excluded_legacy.csv
│   │   ├── discordant.csv
│   │   ├── flag_discordant.csv
│   │   ├── indeterminate.csv
│   │   └── summary_counts.csv
│   ├── conservative/                 % Sensitivity testing only
│   └── midpoint/                     % Sensitivity testing only
└── CatchmentB/
    └── ...
```

**For thermal history modeling, use `kept_strict.csv`** from the `standard/` subfolder. The `conservative/` and `midpoint/` folders exist for sensitivity testing only — do not use them for primary analysis.

**For QTQt/Pecube model setup**, read `_pipeline_summary/pipeline_summary.csv`. The `model_start_Ma` column gives the recommended start time (`mu + NSigma*sigma`) and `model_start_err_Ma` gives the uncertainty (= `sigma_young`, representing within-pluton age spread) for use as the ± of a QTQt bounding box.

---

## Input format

All files are CSV with a header row. Column names must match exactly (case-sensitive).

### ZrnPb.csv — Detrital zircon U–Pb
| Column | Description |
|--------|-------------|
| `ZrnPbDate` | Zircon U–Pb age (Ma) — **required** |
| `ZrnPb2sigerr` | 2σ analytical uncertainty (Ma) — optional but recommended |

### HblAr.csv — Hornblende ⁴⁰Ar/³⁹Ar
| Column | Description |
|--------|-------------|
| `HblGrain` | Grain identifier |
| `HblArDate` | Hornblende Ar/Ar cooling age (Ma) |
| `HblAr1sigerr` | **1σ** analytical uncertainty (Ma) |

### ApHeApPb.csv — Apatite double-dated
| Column | Description |
|--------|-------------|
| `ApGrain` | Grain identifier |
| `ApHeDate` | Apatite (U–Th)/He cooling age (Ma) |
| `ApHe2sigerr` | **2σ** analytical uncertainty (Ma) |
| `ApPbDate` | Apatite U–Pb age (Ma) |
| `ApPb2sigerr` | **2σ** analytical uncertainty (Ma) |

### ZrnHeZrnPb.csv — Zircon double-dated
| Column | Description |
|--------|-------------|
| `ZrnGrain` | Grain identifier |
| `ZrnHeDate` | Zircon (U–Th)/He cooling age (Ma) |
| `ZrnHe2sigerr` | **2σ** analytical uncertainty (Ma) |
| `ZrnPbDate` | Zircon U–Pb crystallization age (Ma) |
| `ZrnPb2sigerr` | **2σ** analytical uncertainty (Ma) |

> **Note on uncertainty conventions:** Hornblende Ar/Ar uncertainties are expected as **1σ**. All He and U–Pb uncertainties are expected as **2σ** and are converted internally to 1σ before use.

---

## Options and parameters

All options are passed as name-value arguments to `run_detrital_pipeline`. See `pipeline_cheatsheet.md` for copy-pasteable examples.

```matlab
run_detrital_pipeline("Catchments", "Output", ...
    NSigma=1.0, ...          % Sigma multiplier for pulse window and model start time (default 1.0)
    BoundsMethod="gmm_ci",...% Pulse window method: "gmm_ci" (default) or "quantile"
    Delta=8, ...             % Magmatic lag window in Ma (default 8)
    P_thresh=0.65, ...       % Probability threshold for exclusion/flagging (default 0.65)
    Kmax=6, ...              % Max GMM components to test via BIC (default 6)
    K_override=0, ...        % Fix K for all catchments; 0 = use BIC (default)
    Nmc=50, ...              % Monte Carlo draws per grain for GMM fitting (default 50)
    run_sensitivity=true)    % Run all 3 P_legacy_mode options (default true)
```

### Pulse window method — `BoundsMethod`

| Method | Window definition | When to use |
|--------|------------------|-------------|
| `"gmm_ci"` (default) | `[mu − NSigma*sigma, mu + NSigma*sigma]` | Recommended. Robust to outlier grains; uses the GMM's own estimate of pulse spread |
| `"quantile"` | Empirical quantiles of posterior-assigned grains + buffer | Fallback if GMM sigma is poorly constrained |

### Model start time — `NSigma`

`NSigma` controls both the pulse window width and the recommended model start time:

```
model_start_Ma     = mu + NSigma * sigma_young
model_start_err_Ma = sigma_young   (use as ± for QTQt bounding box)
```

- `NSigma=1.0` — window spans ~68% of the GMM component; start time = 1 sigma above mean
- `NSigma=1.5` — wider; retains more grains near the upper tail of the pulse
- `NSigma=2.0` — window spans ~95% of the GMM component

**QTQt setup example:** if `mu_young = 85.5` and `sigma_young = 2.2`, with `NSigma=1.0`:
→ `model_start_Ma = 87.7 Ma`, `model_start_err_Ma = 2.2 Ma`
→ Set your oldest QTQt constraint at **87.7 ± 2.2 Ma**

### P_legacy_mode — sensitivity testing

The three output subfolders differ only in which edge of the pulse window is used as the legacy reference age:

| Mode | Reference age | Aggressiveness |
|------|--------------|----------------|
| `standard` | `Tyoung_hi` = `mu + NSigma*sigma` | Least — most grains retained |
| `midpoint` | Midpoint of window | Intermediate |
| `conservative` | `Tyoung_lo` = `mu − NSigma*sigma` | Most — fewest grains retained |

Check `sensitivity_by_system.csv` after a run. Systems flagged `HIGH` (>20% difference between standard and conservative) warrant discussion in your methods. Apatite He typically shows low sensitivity; hornblende Ar/Ar typically shows high sensitivity.

### Filter strictness knobs

| Parameter | Default | Tighter | Effect of tightening |
|-----------|---------|---------|----------------------|
| `P_thresh` | 0.65 | 0.50 | Excludes/flags more borderline grains |
| `Delta` | 8 Ma | 12–15 Ma | Flags more grains as magmatic |
| `P_legacy_mode` | standard | conservative | Removes more legacy grains |

A warning is printed if `P_thresh > 0.75`.

To override K for specific catchments after inspecting QA plots:

```matlab
K_map = containers.Map({"TC", "WP"}, {3, 4});
run_detrital_pipeline("Catchments", "Output", K_override_map=K_map)
```

---

## Output classifications

Each grain in `all_data_classified.csv` is assigned one of six classes:

| Class | Reason code | Meaning | In `kept_strict.csv`? |
|-------|-------------|---------|----------------------|
| `keep_exhumation` | `KE` | Post-magmatic cooling; passes all filters | ✅ Yes |
| `flag_magmatic` | `FM1` / `FM2` | Cooling closely follows crystallization; possible magmatic relaxation | ❌ No (in `kept_plus_flagged.csv`) |
| `exclude_legacy` | `EL1` / `EL2` | High probability of pre-pulse cooling | ❌ No |
| `flag_discordant` | `FD` | He nominally older than U-Pb but within 2σ — inspect before deciding | ❌ No (in `flag_discordant.csv`) |
| `discordant` | `DC` | Crystallization age younger than cooling age beyond 2σ | ❌ No |
| `indeterminate` | `IN` | Missing or zero analytical uncertainty | ❌ No |

`1` and `2` suffixes on reason codes indicate confidence level: `1` = P ≥ 0.90 (high confidence), `2` = P_thresh ≤ P < 0.90 (moderate confidence).

The `confidence` column in `kept_strict.csv` gives a continuous measure of classification certainty for retained grains: `1 − max(P_legacy, P_magmatic)`. Higher values indicate grains less likely to be misclassified.

### Apatite U–Pb independent classification

For apatite double-dated grains, the U–Pb age is assessed independently as a potential mid-temperature cooling constraint (~450–550°C), in addition to the He classification. This never overrides the He class. Three columns are added for apatite rows:

| Column | Description |
|--------|-------------|
| `P_legacy_ApPb` | Probability that the U–Pb age predates the pulse |
| `code_ApPb` | Short code: `KE`, `EL1`, `EL2`, `IN` |
| `class_ApPb` | Full class name: `keep_exhumation`, `exclude_legacy`, `indeterminate` |
| `reason_ApPb` | Short phrase explaining the ApPb classification |

Key combinations:
- `KE` + `KE` — both He and U–Pb are post-pulse; brackets cooling path from ~500°C to ~70°C
- `KE` + `EL` — He retained for thermal modeling; U–Pb predates pulse, do not use as mid-T constraint

---

## Worked example

A self-contained synthetic example is provided in `example/`. It generates toy data for three fictional catchments with known properties and runs the full pipeline.

```matlab
cd example
generate_synthetic_data   % writes Catchments_example/ folder
run_example               % runs pipeline, writes Output_example/
```

Expected runtime: < 1 minute.

---

## Known limitations

**Single-system grains without crystallization ages** (hornblende Ar/Ar) cannot be assessed for magmatic lag and therefore cannot be assigned a `flag_magmatic` classification. In catchments where hornblende ages cluster within the pulse window, sensitivity to `P_legacy_mode` will be high — this is expected behaviour reflecting genuine geological ambiguity rather than a methodological flaw.

**The workflow assumes normally distributed analytical uncertainties.** Asymmetric or non-Gaussian uncertainties are approximated as symmetric Gaussian. If your dataset has strongly asymmetric uncertainties, consider symmetrising them before input.

**GMM component selection via BIC is not guaranteed to recover geologically meaningful components.** Always inspect the QA plots in `ZPb_QA/` for each catchment. If BIC selects a K that looks geologically unreasonable, use `K_override` or `K_override_map` to set K manually.

**The youngest pulse window is inferred from detrital zircon U–Pb only.** If your zircon U–Pb dataset is small (< ~30 grains assigned to the youngest component), the window bounds may be poorly constrained. The code issues a warning in this case.

**No inter-system age concordance is enforced.** Grains are classified independently by system.

**Grains with large analytical uncertainties relative to their age can evade probabilistic exclusion.** When 1-sigma uncertainty is a large fraction of the age itself, the probability distribution for that grain is so wide that P_legacy converges toward ~0.5 regardless of the measured age — meaning the grain will neither be excluded nor flagged, and will pass through as `keep_exhumation` with low confidence. This is the honest probabilistic outcome (we genuinely cannot classify grains we know little about), but it means that grains with poor age precision do not get filtered out. A practical rule of thumb: grains where 1-sigma uncertainty exceeds ~30% of the measured age (i.e. `sig_1sigma / age > 0.30`) are likely to behave this way. For apatite U–Pb in particular, where large uncertainties are common due to low U concentrations, it is recommended to pre-filter your input data on relative uncertainty before running the pipeline. The `confidence` column in `all_data_classified.csv` provides a useful post-hoc check — retained grains with low confidence scores (< 0.4) should be inspected, as they may include high-uncertainty grains that passed through on the basis of an imprecise age rather than a geologically meaningful one.

---

## Citation

If you use this workflow in published work, please cite the associated paper (see `CITATION.cff`) and include the following in your methods:

> Detrital thermochronologic data were filtered using the DetritalChronFilter workflow (Giblin et al., in prep), which classifies grains probabilistically relative to a catchment-specific youngest magmatic pulse window inferred from detrital zircon U–Pb Gaussian mixture modeling.

---

## License

MIT License — see `LICENSE` for details.

---

## Contact

Questions, bug reports, and suggestions are welcome via the [GitHub Issues](../../issues) page.
