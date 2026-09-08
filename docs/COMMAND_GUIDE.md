# DetritalChronFilter command guide

All age-uncertainty inputs use absolute 1σ values in Ma. Use input headings ending in `1sigerr`.

See [`USER_MANUAL.md`](USER_MANUAL.md) for complete input, QA, troubleshooting,
output, and reproducibility guidance.

## Default streamlined run

```matlab
run_detrital_pipeline("Catchments", "Output")
```

This writes one `filter_output` folder per catchment. Optional reference-boundary comparisons and deprecated compatibility files are off by default.

## Restrict target-component selection

```matlab
run_detrital_pipeline("Catchments", "Output_50_200Ma", ...
    TargetComponentAgeRange=[50 200])
```

The GMM still uses the complete valid zircon U–Pb distribution. The range only controls which component means can be selected; the youngest eligible component is chosen.

## Adjust the component window

```matlab
run_detrital_pipeline("Catchments", "Output", NSigma=1.0)
run_detrital_pipeline("Catchments", "Output_1p5sigma", NSigma=1.5)
```

With `BoundsMethod="gmm_sigma_window"`, the window is the selected component mean ± `NSigma` × component standard deviation.

## Compare reference boundaries

```matlab
run_detrital_pipeline("Catchments", "Output_sensitivity", ...
    run_sensitivity=true)
```

This adds one `reference_boundary_comparison.csv` per catchment and one root-level `reference_boundary_sensitivity_summary.csv`. The root table gives the actual boundary ages and model-input counts. It does not create three full result folders.

## Change screening settings

```matlab
run_detrital_pipeline("Catchments", "Output_p50", P_thresh=0.50)
run_detrital_pipeline("Catchments", "Output_delta12", Delta=12)
```

`P_thresh` controls when a probability becomes a screening result or review flag. `Delta` defines the short paired-interval calculation. Both should be reported and justified.

## Override mixture-component count

After inspecting the QA figures:

```matlab
K_map = containers.Map({"WP"}, {3});
run_detrital_pipeline("Catchments", "Output", K_override_map=K_map)
```

## Files to use first

| File | Use |
|---|---|
| `<Catchment>/filter_output/filter_results_full.csv` | Complete one-date-per-row results with full explanations; paired rows share `GrainID` and `PairID` |
| `<Catchment>/filter_output/filter_results_coded.csv` | Compact publication table using IDs from `filter_code_lookup.csv` |
| `<Catchment>/filter_output/model_input_ages.csv` | Dates eligible for downstream modeling, including dates with non-excluding review flags |
| `<Catchment>/filter_output/excluded_ages.csv` | Only dates meeting the older-than-reference exclusion rule |
| `<Catchment>/filter_output/review_flags.csv` | All review flags, with related paired ages included; any exclusion is independently due to the older-than-reference rule |
| `<Catchment>/filter_output/output_summary.csv` | Counts by chronometer and action |
| `pipeline_summary.csv` | Selected-component statistics and candidate model-start ages |
| `output_summary.csv` | Combined filtering counts across catchments |
| `filter_code_lookup.csv` | Code definitions and nominal numeric IDs |
| `<Catchment>/youngest_zircon_component/` | Mixture-model plot and selected-component summary |

Only `OR1` and `OR2` recommend exclusion. `SC1`, `SC2`, `AOI`, `AOU`, and `II` request review. `RC` denotes zircon U–Pb reference context and is not a model-input decision.
