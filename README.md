# DetritalChronFilter

A MATLAB workflow for purpose-specific, probabilistic screening of detrital
chronometer ages against a selected component of a user-defined reference age
distribution.

Current release candidate: **v0.1.0-rc2**.

## What the program does

For each sample, DetritalChronFilter:

1. fits Gaussian mixture models to the complete reference age distribution;
2. selects the youngest eligible component, using BIC or a documented manual K;
3. defines a target-component window and reference boundary;
4. calculates the probability that every model-candidate age is older than that
   reference;
5. calculates optional paired-age order and interval diagnostics when the user
   supplies pair metadata; and
6. writes full, coded, model-input, excluded-age, review, summary, and QA outputs.

Only the older-than-reference result recommends exclusion. Pair diagnostics are
review information and never independently exclude an age.

Chronometer and reference-system names are metadata. The program does not infer
mineral type, closure temperature, kinetics, or expected age order from a label.

## Requirements

- MATLAB R2021a or later
- Statistics and Machine Learning Toolbox
- Ages and absolute 1σ analytical uncertainties reported in Ma

## Input layout

For multiple samples, each sample folder contains two CSV files:

```text
Samples/
├── SampleA/
│   ├── ReferenceDistribution.csv
│   └── ChronometerData.csv
└── SampleB/
    ├── ReferenceDistribution.csv
    └── ChronometerData.csv
```

For one sample, users may either retain that structure or pass the sample folder
itself. The program accepts the two files directly when both occur in the
selected folder. It does not mix direct files and sample subfolders in one run.

### `ReferenceDistribution.csv`

```csv
ReferenceSystem,GrainID,Age_Ma,Age_1sigma_Ma
ZrnUPb,R001,92.4,1.3
ZrnUPb,R002,95.1,1.1
ZrnUPb,R003,181.6,2.2
```

| Column | Requirement | Meaning |
|---|---|---|
| `ReferenceSystem` | Required | One system label repeated for the complete distribution; metadata only |
| `GrainID` | Required | Unique identifier for each reference age |
| `Age_Ma` | Required | Measured age in Ma |
| `Age_1sigma_Ma` | Required | Absolute 1σ analytical uncertainty in Ma |

This file must contain the complete distribution supplied to the mixture model,
not only ages believed to belong to the target component. The system may be
zircon U–Pb or another independently justified reference chronometer.

### `ChronometerData.csv`

```csv
Chronometer,GrainID,Age_Ma,Age_1sigma_Ma,PairID,PairRole,UseForModel
ZrnHe,Z001,48.2,2.1,Zrn-Z001,expected_younger,true
ZrnUPb,Z001,92.4,1.3,Zrn-Z001,expected_older,false
ApHe,A001,35.8,1.8,Ap-A001,expected_younger,true
ApUPb,A001,65.1,1.1,Ap-A001,expected_older,true
HblAr,H001,78.2,1.5,,,true
```

| Column | Requirement | Meaning |
|---|---|---|
| `Chronometer` | Required | User-defined grouping label; not scientifically interpreted by the code |
| `GrainID` | Required | Identifier unique within a chronometer |
| `Age_Ma` | Required | Measured age in Ma |
| `Age_1sigma_Ma` | Required | Absolute 1σ analytical uncertainty in Ma |
| `PairID` | Optional | Shared identifier for two analyses of one grain |
| `PairRole` | Optional | `expected_younger`, `expected_older`, or blank |
| `UseForModel` | Optional | `true` or `false`; blank or omitted defaults to `true` |

The optional column headings are included in the template. Users without paired
analyses leave `PairID` and `PairRole` blank.

When `PairRole` is used, a `PairID` must occur exactly twice, with one
`expected_younger` and one `expected_older` row. Two blank roles preserve the
pair relationship without requesting an age-order calculation.

Set `UseForModel=false` when an age is included only to provide pair or contextual
information. Such a row is retained in the full results but is not screened as a
model-input candidate.

`ChronometerData.csv` may contain one chronometer or any number of chronometers.

## Quick start

Open MATLAB in the repository folder and add it to the path:

```matlab
addpath(pwd)
run_detrital_pipeline("Samples", "Output")
```

For a single sample whose two CSVs are directly inside `SampleA/`:

```matlab
run_detrital_pipeline("SampleA", "Output")
```

To limit which fitted component means can be selected:

```matlab
run_detrital_pipeline("Samples", "Output", ...
    TargetComponentAgeRange=[50 200])
```

The mixture model still uses every valid age in `ReferenceDistribution.csv`.
The range affects component selection only.

## Outputs

```text
Output/
├── README.txt
├── pipeline_summary.csv
├── output_summary.csv
├── filter_code_lookup.csv
└── SampleA/
    ├── target_component/
    │   ├── target_component_plot.png
    │   └── target_component_summary.csv
    └── filter_output/
        ├── filter_results_full.csv
        ├── filter_results_coded.csv
        ├── model_input_ages.csv
        ├── excluded_ages.csv
        ├── review_flags.csv
        └── output_summary.csv
```

| File | Purpose |
|---|---|
| `filter_results_full.csv` | Complete one-row-per-analysis results with probabilities and explanations |
| `filter_results_coded.csv` | Compact table using numeric IDs from `filter_code_lookup.csv` |
| `model_input_ages.csv` | Ages eligible for downstream modeling, including non-excluding review flags |
| `excluded_ages.csv` | Only model-candidate ages assigned OR1 or OR2 |
| `review_flags.csv` | All review-flagged rows with the related paired age beside them when available |
| Sample `output_summary.csv` | Observed and model-input age ranges and counts for each chronometer |
| Root `pipeline_summary.csv` | Reference system, target-component parameters, K selection, and candidate start age |
| Root `output_summary.csv` | Combined chronometer summaries across samples |

The summary reports observed minimum and maximum ages plus the minimum, median,
and maximum of the retained model inputs. These descriptive ranges help users
spot cross-chronometer patterns, but the code does not apply a closure-temperature
ordering rule.

## Reference and review codes

| Code | Role | Meaning | Default action |
|---|---|---|---|
| `RT` | Reference screen | Older-than-reference probability is below the threshold | Retain |
| `OR1` | Reference screen | High probability that the age is older than the reference | Exclude |
| `OR2` | Reference screen | Moderate probability that the age is older than the reference | Exclude |
| `SI1` | Pair review | High probability of a short paired-age interval | Review only |
| `SI2` | Pair review | Moderate probability of a short paired-age interval | Review only |
| `AOI` | Pair review | Expected-younger age is older beyond combined 2σ uncertainty | Review only |
| `AOU` | Pair review | Expected order is reversed nominally but unresolved within combined 2σ | Review only |
| `II` | Reference/review | Required age or uncertainty is missing or nonpositive | Review |
| `CX` | Context | `UseForModel=false`; row is not screened as model input | Context only |
| `NF` | Review | No separate review flag | None |

All input uncertainty columns are 1σ. AOI/AOU use two times the combined 1σ
uncertainty as an internal overlap criterion; this is not a 2σ input convention.

## Main options

| Option | Default | Meaning |
|---|---:|---|
| `Kmax` | `6` | Largest component count evaluated by BIC |
| `K_override` | `0` | Fixed K for all samples; zero uses BIC |
| `K_override_map` | empty | Sample-specific K values in a `containers.Map` |
| `TargetComponentAgeRange` | unrestricted | Eligible component-mean range |
| `BoundsMethod` | `"gmm_sigma_window"` | Target-window method |
| `NSigma` | `1.0` | Component-standard-deviation multiplier |
| `P_thresh` | `0.65` | OR2 and SI2 probability threshold |
| `Delta` | `8` Ma | Upper bound for the optional short-pair-interval review |
| `Nmc` | `50` | Monte Carlo draws per reference age |
| `run_sensitivity` | `false` | Compare older edge, midpoint, and younger edge |

## Examples and testing

Run the complete synthetic example:

```matlab
run("examples/synthetic/run_example.m")
```

Run the release tests:

```matlab
run_release_tests
```

See `docs/USER_MANUAL.md` for detailed preparation, QA, troubleshooting, and
reporting guidance. Citation metadata are in `CITATION.cff`; the software is
distributed under the MIT License.
