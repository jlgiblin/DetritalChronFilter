# DetritalChronFilter command guide

Each catchment folder needs `ReferenceDistribution.csv` and
`ChronometerData.csv`. All ages and absolute 1σ uncertainties use Ma.

## Default run

```matlab
addpath(pwd)
run_detrital_pipeline("Catchments", "Output")
```

## Restrict target-component selection

```matlab
run_detrital_pipeline("Catchments", "Output", ...
    TargetComponentAgeRange=[50 200])
```

The full reference distribution remains in the GMM fit.

## Catchment-specific K override

```matlab
K_map = containers.Map({'CatchmentC'}, {3});
run_detrital_pipeline("Catchments", "Output", K_override_map=K_map)
```

Use an override only after inspecting the target-component QA output.

## Adjust target-window width

```matlab
run_detrital_pipeline("Catchments", "Output", NSigma=1.5)
```

## Compare reference boundaries

```matlab
run_detrital_pipeline("Catchments", "Output", run_sensitivity=true)
```

This adds condensed older-edge, midpoint, and younger-edge comparisons.

## Run the synthetic example

```matlab
run("examples/synthetic/run_example.m")
```

## Run the tests

```matlab
run_release_tests
```

## Main outputs

| Path | Purpose |
|---|---|
| `<Catchment>/target_component/` | Mixture-model plot and selected-component summary |
| `<Catchment>/filter_output/filter_results_full.csv` | Full one-row-per-analysis results |
| `<Catchment>/filter_output/filter_results_coded.csv` | Compact publication table |
| `<Catchment>/filter_output/model_input_ages.csv` | Retained model-input candidates |
| `<Catchment>/filter_output/excluded_ages.csv` | OR1/OR2 exclusions only |
| `<Catchment>/filter_output/review_flags.csv` | Pair and information review flags |
| `<Catchment>/filter_output/output_summary.csv` | Observed and retained ranges by chronometer |
| `pipeline_summary.csv` | Target-component and K-selection record |
| `output_summary.csv` | Combined chronometer summaries |
| `filter_code_lookup.csv` | Code and numeric-ID definitions |

See `USER_MANUAL.md` for input-column definitions and interpretation guidance.
