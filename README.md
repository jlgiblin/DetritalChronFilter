# DetritalChronFilter

A MATLAB workflow for purpose-specific, probabilistic screening of detrital thermochronologic data using a selected statistical component of the zircon U–Pb age distribution as the reference.

Current release candidate: **v0.1.0-rc1**. This version intentionally supports the documented four-file input layout. Generic chronometer imports and additional component-selection methods are planned, but are not part of v0.1.0.

## Documentation

- [`docs/USER_MANUAL.md`](docs/USER_MANUAL.md) provides the complete setup,
  operation, output, quality-control, troubleshooting, and reporting guidance.
- [`docs/COMMAND_GUIDE.md`](docs/COMMAND_GUIDE.md) is a short command reference.
- [`input_templates/INPUT_DATA_CHECKLIST.md`](input_templates/INPUT_DATA_CHECKLIST.md)
  summarizes input preparation.
- [`CHANGELOG.md`](CHANGELOG.md) records release changes.

## Terminology and decision policy

The public results describe calculated age relationships rather than assigning a geological cause to an individual analysis.

- `older_than_reference` means that the probability of an age being older than the selected reference reaches the user-set threshold. This is the only result that recommends exclusion.
- `eligible_after_reference_screen` means that an age did not reach that exclusion threshold.
- `age_order_inconsistent`, `age_order_unresolved`, and `short_crystallization_cooling_interval` are review flags. They do not automatically exclude an age.
- `reference_context_not_screened` identifies paired zircon U–Pb dates that provide reference context but are not candidate model-input ages.

Users retain responsibility for interpreting the geological cause of these patterns and for deciding whether review-flagged ages suit their research question.

## Workflow

1. Fit Gaussian mixture models to each catchment's complete valid detrital zircon U–Pb distribution.
2. Select the youngest statistical component whose mean falls within `TargetComponentAgeRange`. Components outside this search range remain in the fit but cannot be selected.
3. Define the target-component window and use its older edge as the primary reference age.
4. Calculate each eligible date's probability of being older than that reference.
5. Calculate paired-age order and crystallization-to-cooling interval diagnostics as non-excluding review information.
6. Write a small set of analysis-level, model-input, excluded-age, summary, and QA outputs.

Supported systems are hornblende 40Ar/39Ar, apatite (U–Th)/He with paired apatite U–Pb, and zircon (U–Th)/He with paired zircon U–Pb.

## Requirements

- MATLAB R2021a or later
- Statistics and Machine Learning Toolbox

## Input layout

Place input CSVs in one folder per catchment:

```text
Catchments/
├── CatchmentA/
│   ├── ZrnPb.csv
│   ├── HblAr.csv
│   ├── ApHeApPb.csv
│   └── ZrnHeZrnPb.csv
└── CatchmentB/
    └── ...
```

Required columns are:

| File | Columns |
|---|---|
| `ZrnPb.csv` | `ZrnPbDate`, `ZrnPb1sigerr` |
| `HblAr.csv` | `HblGrain`, `HblArDate`, `HblAr1sigerr` |
| `ApHeApPb.csv` | `ApGrain`, `ApHeDate`, `ApHe1sigerr`, `ApPbDate`, `ApPb1sigerr` |
| `ZrnHeZrnPb.csv` | `ZrnGrain`, `ZrnHeDate`, `ZrnHe1sigerr`, `ZrnPbDate`, `ZrnPb1sigerr` |

All datasets must include absolute 1σ uncertainties in Ma, using the required column names ending in `1sigerr`. Other uncertainty conventions are not accepted.

## Quick start

Download or clone the repository, open MATLAB in the repository folder, and add the code to the MATLAB path:

```matlab
addpath(pwd)
```

```matlab
run_detrital_pipeline("Catchments", "Output", ...
    TargetComponentAgeRange=[50 200])
```

The default run creates one primary result set per catchment:

```text
Output/
├── README.txt
├── pipeline_summary.csv
├── output_summary.csv
├── filter_code_lookup.csv
└── CatchmentA/
    ├── youngest_zircon_component/
    └── filter_output/
        ├── filter_results_full.csv
        ├── filter_results_coded.csv
        ├── model_input_ages.csv
        ├── excluded_ages.csv
        ├── review_flags.csv
        └── output_summary.csv
```

The six catchment tables have distinct roles:

| File | Purpose |
|---|---|
| `filter_results_full.csv` | Complete review table, with one dated analysis per row and full explanations |
| `filter_results_coded.csv` | Compact publication table using numeric result and review IDs |
| `model_input_ages.csv` | Dates eligible for downstream modeling, including non-excluded dates carrying review flags |
| `excluded_ages.csv` | Only dates that meet the older-than-reference exclusion rule |
| `review_flags.csv` | All dates carrying a review flag, with the related paired date repeated beside them |
| `output_summary.csv` | One row per chronometer with observed and model-input age ranges, model-input median, and exclusion/review counts |

Paired measurements occupy separate rows in both filter-results tables and share `GrainID` and `PairID`. This preserves the relationship while allowing, for example, apatite He and apatite U–Pb dates to have separate filtering decisions. The file is named `excluded_ages.csv`, rather than “excluded grains,” because two dates from one paired grain may receive different results.

The root-level `pipeline_summary.csv` records target-component and model-start information for all catchments. The root-level `output_summary.csv` combines the chronometer summaries across catchments, and `filter_code_lookup.csv` defines the numeric IDs once for the entire run.

The observed minimum and maximum include every valid reported date. The model-input minimum, median, and maximum describe only dates retained after the reference screen, including non-excluding review flags. Comparing these ranges can reveal patterns worth reviewing, but the summary does not assign closure temperatures or label a cross-chronometer ordering problem.

In the full table, `ReviewCode="NF"` means that no separate review flag was assigned. The corresponding coded-table value is `ReviewFlagID=0`. All numeric IDs map directly to `filter_code_lookup.csv`; they are identifiers, not ranked scores. Unpaired rows use `PairID="not_paired"` so spreadsheet and MATLAB imports preserve the column as text.

A review flag never causes exclusion. A row in `review_flags.csv` can nevertheless have `Action="exclude"` if that same date independently meets the older-than-reference rule; `ActionReason` makes that distinction explicit.

For reliable programmatic import in MATLAB, specify the comma delimiter explicitly:

```matlab
T = readtable("filter_results_full.csv", ...
    Delimiter=",", VariableNamingRule="preserve", TextType="string");
```

Explicit delimiter selection prevents MATLAB from mistaking prose-containing CSV files for whitespace-delimited text.

## Optional reference-boundary comparison

The default is `run_sensitivity=false`. To compare three reference choices without producing three duplicate folder trees:

```matlab
run_detrital_pipeline("Catchments", "Output", ...
    TargetComponentAgeRange=[50 200], ...
    run_sensitivity=true)
```

This adds:

```text
Output/
├── reference_boundary_sensitivity_summary.csv
└── CatchmentA/
    └── sensitivity/
        └── reference_boundary_comparison.csv
```

The root summary reports the actual boundary ages and resulting model-input counts for the primary older edge, midpoint, and younger edge. The catchment comparison gives analysis-level results side by side. This is a stress test of the boundary choice—not a separate “conservative” interpretation—and the code no longer assigns vague high/moderate/low sensitivity labels.

## Included examples

Run the complete synthetic example from the repository root:

```matlab
run("examples/synthetic/run_example.m")
```

The example writes generated results to `examples/synthetic/output/`, which is ignored by Git.

## Reference and review codes

| Code | Role | Meaning | Action |
|---|---|---|---|
| `RT` | Reference screen | Probability is below the exclusion threshold | Retain unless a review flag applies |
| `OR1` | Reference screen | High probability that the date is older than the reference | Exclude |
| `OR2` | Reference screen | Moderate probability that the date is older than the reference | Exclude |
| `SC1` | Review | High probability of a short crystallization-to-cooling interval | Review; do not exclude |
| `SC2` | Review | Moderate probability of a short interval | Review; do not exclude |
| `AOI` | Review | Cooling age is older than its paired U–Pb age beyond combined 2σ uncertainty | Review; do not exclude |
| `AOU` | Review | Nominal age order is unexpected but overlaps within combined 2σ uncertainty | Review; do not exclude |
| `II` | Reference/review | A required date or uncertainty is missing or nonpositive | Review |
| `NA` | Paired-age status | A paired age was not provided | Not applicable |
| `RC` | Reference context | Paired zircon U–Pb date provides target-component context and is not screened as model input | Reference only |

`CodeID` values in the lookup table are nominal identifiers only. They are not ranks or probability values.

## Important columns in `filter_results_full.csv`

| Column | Meaning |
|---|---|
| `GrainID`, `PairID` | Identifiers that keep paired dates linked |
| `Chronometer` | Date type for the individual row |
| `PairRole` | `cooling_age`, `single_age`, or `paired_u_pb_age` |
| `Age_Ma`, `Age_1sigma_Ma` | Date and absolute 1σ uncertainty |
| `ReferenceAge_Ma` | Reference boundary used for the calculation |
| `P_OlderThanReference` | Probability that the row's date is older than the reference |
| `ReferenceClass`, `ReferenceCode` | Reference-screen result in full and coded form |
| `ReviewRecommended`, `ReviewCode` | Non-excluding review status |
| `P_ShortInterval` | Probability that the paired interval is positive and shorter than `Delta` |
| `PairInterval_Ma` | Paired U–Pb date minus cooling age |
| `Action` | `retain`, `review`, `exclude`, or `reference_only` |
| `ModelInclude` | Whether the date appears in `model_input_ages.csv` |
| `P_Threshold` | Probability threshold used for the decision |
| `ShortIntervalThreshold_Ma` | Interval threshold used for SC1/SC2, where applicable |
| `ActionReason` | Plain-language statement combining the numerical result, threshold, and action |

## Main options

```matlab
run_detrital_pipeline("Catchments", "Output", ...
    NSigma=1.0, ...
    BoundsMethod="gmm_sigma_window", ...
    Delta=8, ...
    P_thresh=0.65, ...
    Kmax=6, ...
    K_override=0, ...
    Nmc=50, ...
    TargetComponentAgeRange=[50 200], ...
    run_sensitivity=false)
```

- `TargetComponentAgeRange` controls which component means are eligible for selection; it does not remove dates before fitting.
- `BoundsMethod="gmm_sigma_window"` defines the component window as its mean ± `NSigma` × component standard deviation.
- `Delta` defines the upper bound of the short-interval review calculation. It is a user-set sensitivity parameter, not a universal process boundary.
- `P_thresh` is the probability required to assign an older-than-reference result or short-interval review flag.
- `K_override` or `K_override_map` can be used after inspecting the mixture-model QA output.

## Candidate model-start age

The pipeline reports:

```text
candidate_model_start_Ma = target_component_mean_Ma + NSigma × target_component_sigma_Ma
```

This is a workflow-derived candidate boundary. `target_component_sigma_Ma` describes dispersion within the selected statistical component; it is not uncertainty on the candidate model-start age.

## Scientific cautions

- Mixture components remain statistical components until geological meaning is evaluated using independent evidence.
- A short paired interval or unexpected nominal age order does not uniquely determine its cause.
- An age passing the reference screen is eligible for the stated workflow; it is not thereby validated for every possible use.
- Single-system ages without paired U–Pb dates cannot be evaluated for paired-age order or interval.
- Large analytical uncertainty can keep a probability below the decision threshold. Inspect the probabilities and uncertainties, not only the action label.
- The automatic model-input table is purpose-specific. Other research questions may require different inclusion decisions.

## Testing

Run the release test from the repository root:

```matlab
run_release_tests
```

The tests generate temporary synthetic data, run automatic and overridden component selection, verify the output layout and public result-field names, check the coded/full table mapping, and confirm that review flags do not cause exclusion.

## Scope and roadmap

The supported v0.1.0 interface is the four-file catchment layout documented above. Proposed generic chronometer input, unpaired-system support, and alternative component estimators are described in [`docs/ROADMAP.md`](docs/ROADMAP.md) and should not be cited as implemented features.

## Citation status

Citation metadata are provided in `CITATION.cff`. An archived-release DOI can be added after the public v0.1.0 release is connected to an archival service such as Zenodo.
