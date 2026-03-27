# MultichronFitTSF — Pipeline Command Cheat Sheet

## Basic run (all defaults)
```matlab
run_detrital_pipeline("Catchments", "Output")
```

---

## Key knobs and what they do

### Pulse window width — `NSigma`
Controls how wide the filtering window is around the youngest GMM pulse mean.
Window = `[mu - NSigma*sigma,  mu + NSigma*sigma]`
Model start time = `mu + NSigma*sigma`

```matlab
% Default (1 sigma — recommended starting point)
run_detrital_pipeline("Catchments", "Output", NSigma=1.0)

% Wider window (more permissive — retains grains further from pulse mean)
run_detrital_pipeline("Catchments", "Output_1p5sig", NSigma=1.5)

% Even wider
run_detrital_pipeline("Catchments", "Output_2sig", NSigma=2.0)
```

**QTQt bounding box:** use `model_start_Ma` as the center and `sigma_young`
as the +/- uncertainty. Both are reported in `_pipeline_summary/pipeline_summary.csv`.
Pecube does not use an uncertainty on start time — use `model_start_Ma` directly.

---

### Pulse window method — `BoundsMethod`
Controls how the window boundaries are calculated.

```matlab
% gmm_ci (default) — window from GMM component mean ± NSigma*sigma
% Robust to outlier grains at the tails. Recommended.
run_detrital_pipeline("Catchments", "Output", BoundsMethod="gmm_ci")

% quantile — window from empirical quantiles of assigned grains + buffer
% Original method. More sensitive to tail outliers.
% NSigma still controls model_start_Ma even in quantile mode.
run_detrital_pipeline("Catchments", "Output", BoundsMethod="quantile")
```

---

### Legacy filter strictness — `P_thresh`
Probability threshold for excluding/flagging grains.
Lower = stricter (fewer grains pass). A warning prints if P_thresh > 0.75.

```matlab
% Default
run_detrital_pipeline("Catchments", "Output", P_thresh=0.65)

% Stricter
run_detrital_pipeline("Catchments", "Output_strict", P_thresh=0.50)

% Relaxed (warning will print)
run_detrital_pipeline("Catchments", "Output_relaxed", P_thresh=0.80)
```

---

### Legacy filter reference age — `P_legacy_mode`
This is set automatically per output folder when `run_sensitivity=true` (default).
Standard mode uses `Tyoung_hi` (mu + NSigma*sigma) as the reference.
Conservative mode uses `Tyoung_lo` (mu - NSigma*sigma).

Use the sensitivity output to decide which result to report:
- If `N_keep_standard` ≈ `N_keep_conservative` → results are robust; use standard
- If they differ substantially → report both in supplementary material

---

### Magmatic lag window — `Delta`
Grains whose cooling age falls within `Delta` Ma of their crystallization age
are flagged as potentially magmatic rather than exhumation-driven.
Geologically motivated: ~8 Ma spans ~400°C to ~70°C at Sierra Nevada cooling rates.

```matlab
% Default (8 Ma)
run_detrital_pipeline("Catchments", "Output", Delta=8)

% Wider — flags more grains as potentially magmatic
run_detrital_pipeline("Catchments", "Output_delta12", Delta=12)
```

---

### Fix K for a specific catchment (after inspecting QA plots)
```matlab
% Override K for all catchments
run_detrital_pipeline("Catchments", "Output", K_override=4)

% Override K per catchment
K_map = containers.Map({"TC", "WP"}, {3, 4});
run_detrital_pipeline("Catchments", "Output", K_override_map=K_map)
```

---

### Run standard mode only (faster, no sensitivity CSVs)
```matlab
run_detrital_pipeline("Catchments", "Output", run_sensitivity=false)
```

---

## Combining knobs
Parameters can be combined freely:
```matlab
run_detrital_pipeline("Catchments", "Output_tight", ...
    NSigma=1.5, ...
    P_thresh=0.55, ...
    Delta=10, ...
    BoundsMethod="gmm_ci")
```

---

## Output files — what to use

| File | Use for |
|---|---|
| `<Catchment>/standard/kept_strict.csv` | **Primary input to TSF and Pecube** |
| `<Catchment>/standard/all_data_classified.csv` | Full grain-level inspection |
| `_pipeline_summary/pipeline_summary.csv` | Model start times, pulse windows across catchments |
| `_pipeline_summary/sensitivity_by_system.csv` | Standard vs conservative grain count comparison |
| `<Catchment>/ZPb_QA/*.png` | Visual check of GMM fit and pulse window — inspect before use |
| `<Catchment>/conservative/kept_strict.csv` | Sensitivity comparison only |

---

## Reading `pipeline_summary.csv` for QTQt setup

| Column | Meaning |
|---|---|
| `mu_young` | GMM mean of youngest pulse (Ma) |
| `sigma_young` | GMM std of youngest pulse (Ma) — use as ± for QTQt bounding box |
| `model_start_Ma` | `mu + NSigma*sigma` — recommended model start time |
| `Tyoung_lo / hi` | Full pulse window passed to the filter |
| `NSigma` | NSigma value used in this run |

**QTQt example for TC:**
If `mu_young = 85.5`, `sigma_young = 2.2`, `model_start_Ma = 87.7`:
→ Set QTQt oldest constraint at **87.7 ± 2.2 Ma**

