# DetritalChronFilter user manual

DetritalChronFilter screens ages from one or more chronometer systems against a
selected component of a user-defined reference age distribution. It calculates
age relationships and probabilities; it does not assign a geological mechanism
to an individual analysis.

This manual describes the `0.1.0-rc2` two-file workflow.

## 1. Requirements

- MATLAB R2021a or later
- Statistics and Machine Learning Toolbox
- Ages in Ma
- Absolute analytical uncertainties reported at 1σ

The program does not convert 2σ inputs.

## 2. Package contents

| File or folder | Purpose |
|---|---|
| `run_detrital_pipeline.m` | Recommended entry point for one or more samples |
| `convert_fixed_inputs_to_two_file.m` | Converts original four-file inputs without overwriting them |
| `infer_target_component.m` | Fits reference-distribution mixture models and selects the target component |
| `filter_detrital_thermo.m` | Applies the reference screen and optional pair diagnostics |
| `input_templates/` | Header-only templates and preparation checklist |
| `examples/synthetic/` | Complete two-file example |
| `docs/COMMAND_GUIDE.md` | Short command reference |
| `run_release_tests.m` | Automated test entry point |

## 3. Folder structure

Create one folder per sample. Each folder must contain the two recognized
input filenames; unrelated files are ignored.

```text
Samples/
├── SampleA/
│   ├── ReferenceDistribution.csv
│   └── ChronometerData.csv
└── SampleB/
    ├── ReferenceDistribution.csv
    └── ChronometerData.csv
```

For one sample, the selected input folder may itself contain the two CSVs. For
multiple samples, select their parent folder. The program does not mix these two
layouts in one run, and it reports incomplete sample folders as input errors.

## 4. Reference distribution

`ReferenceDistribution.csv` supplies the complete distribution used for Gaussian
mixture modeling.

```csv
ReferenceSystem,GrainID,Age_Ma,Age_1sigma_Ma
ZrnUPb,R001,92.4,1.3
ZrnUPb,R002,95.1,1.1
ZrnUPb,R003,181.6,2.2
```

All four columns are required:

- `ReferenceSystem`: one nonblank label repeated throughout the file;
- `GrainID`: a unique nonblank identifier;
- `Age_Ma`: measured age in Ma; and
- `Age_1sigma_Ma`: positive absolute 1σ uncertainty in Ma.

The system label is reported in outputs but is not used to select a model,
closure temperature, or interpretation. A reference system other than zircon
U–Pb may be used when scientifically justified.

Do not pre-trim the file to the desired component. Components outside
`TargetComponentAgeRange` remain necessary for representing the complete
distribution even though they cannot be selected as the target.

## 5. Chronometer data

`ChronometerData.csv` is a long-format table containing one dated analysis per
row. It may contain one chronometer or any number of chronometers.

```csv
Chronometer,GrainID,Age_Ma,Age_1sigma_Ma,PairID,PairRole,UseForModel
ZrnHe,Z001,48.2,2.1,Zrn-Z001,expected_younger,true
ZrnUPb,Z001,92.4,1.3,Zrn-Z001,expected_older,false
ApHe,A001,35.8,1.8,Ap-A001,expected_younger,true
ApUPb,A001,65.1,1.1,Ap-A001,expected_older,true
HblAr,H001,78.2,1.5,,,true
```

### Required columns

| Column | Rule |
|---|---|
| `Chronometer` | Nonblank user-defined label; groups rows in outputs only |
| `GrainID` | Nonblank and unique within each `Chronometer` label |
| `Age_Ma` | Measured age in Ma |
| `Age_1sigma_Ma` | Absolute 1σ uncertainty in Ma |

`Chronometer` text is never mapped to a mineral, closure temperature, or kinetic
model. For example, `ZrnHe`, `ZHe`, and `CustomSystem` are simply distinct labels.
Use one spelling consistently if the rows should remain in the same group.

### Optional pairing columns

The template includes the optional columns so users can leave cells blank rather
than changing the header.

- `PairID` links exactly two rows representing analyses of the same grain.
- `PairRole` may be `expected_younger`, `expected_older`, or blank.
- `UseForModel` may be `true` or `false`; blank or an absent column defaults to
  `true`.

Supported `UseForModel` entries are `true/false`, `yes/no`, or `1/0`.

For an ordered pair, the same `PairID` must occur exactly twice, with one row of
each role. The role communicates the user's expected age relation; it is not
inferred from `Chronometer`. If both roles are blank, the pair remains linked in
the output but no order or interval calculation is performed.

For an unpaired analysis, leave both `PairID` and `PairRole` blank. Do not create
a one-row PairID.

Set `UseForModel=false` for an age that is supplied only as paired or contextual
information. The row receives code `CX`, remains available in the full results
and as related pair information when applicable, and is omitted from
`model_input_ages.csv`.

## 6. Running the program

Open MATLAB in the repository folder:

```matlab
addpath(pwd)
run_detrital_pipeline("Samples", "Output")
```

To restrict the eligible target-component means:

```matlab
run_detrital_pipeline("Samples", "Output", ...
    TargetComponentAgeRange=[50 200])
```

To apply a reviewed component-count override to one sample:

```matlab
K_map = containers.Map({'SampleC'}, {3});
run_detrital_pipeline("Samples", "Output", ...
    K_override_map=K_map)
```

The pipeline records whether K was selected by BIC or supplied manually.

### Original four-file inputs

Users updating from the original public version can convert a copy of their
inputs into the current layout:

```matlab
convert_fixed_inputs_to_two_file("OriginalInputs", "ConvertedInputs")
```

The destination must differ from the source. Existing converted CSVs are not
overwritten. The converter retains apatite U–Pb as a model candidate and paired
zircon U–Pb as context, matching the prior public workflow.

## 7. Target-component selection

For each sample, the program:

1. reads the complete valid reference distribution;
2. propagates 1σ analytical uncertainties through Monte Carlo sampling;
3. fits K = 2 through `Kmax` and selects the lowest-BIC fit unless K is overridden;
4. identifies components whose means fall inside `TargetComponentAgeRange`;
5. selects the youngest eligible component; and
6. defines the target window from the selected component.

With `BoundsMethod="gmm_sigma_window"`, the window is:

```text
component mean ± NSigma × component standard deviation
```

The default `NSigma=1` represents the central approximately 68% of a Gaussian
component. Component standard deviation describes age dispersion; it is not the
uncertainty on the component mean.

Inspect `target_component_plot.png` and `target_component_summary.csv` before
accepting the component selection.

## 8. Reference screen

The primary reference is the older edge of the target-component window. For each
valid row with `UseForModel=true`, the program calculates:

```text
P(age > reference age)
```

With default `P_thresh=0.65`:

- `OR1`: probability at least 0.90; recommend exclusion;
- `OR2`: probability at least 0.65 but below 0.90; recommend exclusion; and
- `RT`: probability below 0.65; retain after this reference screen.

OR1 and OR2 are the only automatic exclusion recommendations.

## 9. Optional pair diagnostics

For a validated ordered pair, the program calculates:

```text
pair interval = expected-older age − expected-younger age
combined 1σ uncertainty = sqrt(σolder² + σyounger²)
```

The following flag is attached to the `expected_younger` row:

- `AOI`: expected order is reversed beyond two times combined 1σ uncertainty;
- `AOU`: nominal order is reversed but unresolved within that overlap criterion;
- `SI1`: at least 0.90 probability that the positive interval is below `Delta`;
- `SI2`: probability at least `P_thresh` but below 0.90 for that interval; or
- `II`: pair information required for the calculation is missing or invalid.

These are review flags only. A flagged row remains in the model-input table when
it passes the reference screen.

The program does not compare age distributions according to closure-temperature
ordering. The observed and retained age ranges in `output_summary.csv` allow such
patterns to be inspected manually or handled in a downstream model.

## 10. Outputs

Each sample receives:

```text
SampleA/
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

The root output directory contains a combined `pipeline_summary.csv`, combined
`output_summary.csv`, shared `filter_code_lookup.csv`, and a run-specific
`README.txt`.

`output_summary.csv` gives one row per chronometer with:

- reported and valid age counts;
- complete observed minimum and maximum;
- retained model-input count, minimum, median, and maximum;
- excluded count;
- review-flag count; and
- context-only count.

The coded table uses nominal identifiers, not ranks. Always distribute
`filter_code_lookup.csv` with a coded results table.

## 11. Reference-boundary sensitivity

Set `run_sensitivity=true` to compare the older edge, midpoint, and younger edge:

```matlab
run_detrital_pipeline("Samples", "Output", run_sensitivity=true)
```

This adds one analysis-level comparison per sample and one combined summary.
It does not create three duplicate output-folder trees.

## 12. Quality-control checklist

Before using the model-input table:

1. Confirm that all ages and uncertainties are in Ma and 1σ.
2. Confirm that `ReferenceSystem` contains one consistent label.
3. Inspect the BIC curve and fitted reference distribution.
4. Confirm that the selected component is appropriate for the stated objective.
5. Record and justify any K override or component search range.
6. Review OR1/OR2 ages and all pair flags.
7. Compare observed and retained age ranges among chronometers.
8. Preserve the run README, summaries, code lookup, and full result table.

## 13. Common input errors

| Message | Likely cause | Correction |
|---|---|---|
| Missing recognized file | One of the two filenames is absent | Use the supplied filenames exactly |
| Missing required column | Header differs from the template | Copy the template header |
| Duplicate identifier | Same chronometer and GrainID occur more than once | Make GrainID unique within that system |
| PairID must occur twice | Pair identifier is missing its companion row | Add the paired row or clear PairID |
| Invalid PairRole | Pair roles are incomplete or inconsistent | Use one expected-younger and one expected-older role |
| No eligible component | Search range excludes every fitted component | Inspect QA and revise the range if justified |

## 14. Reproducibility

Record the software version, input tables, K settings, target-component range,
window method, `NSigma`, `P_thresh`, `Delta`, sensitivity setting, and all manual
decisions. Preserve `README.txt` from each run with its output tables.

Run the automated tests from the repository root:

```matlab
run_release_tests
```
