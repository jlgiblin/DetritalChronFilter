# DetritalChronFilter user manual

DetritalChronFilter is a MATLAB workflow for screening detrital thermochronologic
ages relative to a selected component of a catchment's detrital zircon U-Pb age
distribution. It calculates age relationships and probabilities. It does not
assign a geological mechanism to an individual analysis.

This manual describes the supported `0.1.0` workflow. The current interface uses
four CSV files per catchment. Proposed generic chronometer imports and additional
component-selection methods are described separately in `ROADMAP.md` and are not
implemented in this release.

## Requirements

- MATLAB R2021a or later
- Statistics and Machine Learning Toolbox
- Input ages and absolute 1 sigma analytical uncertainties reported in Ma

The automated release tests were run with MATLAB R2025b.

## Files in the package

| File or folder | Purpose |
|---|---|
| `run_detrital_pipeline.m` | Recommended entry point for one or more catchments |
| `infer_youngest_pulse_from_ZPb.m` | Fits zircon U-Pb mixture models and selects the target component |
| `filter_detrital_thermo.m` | Applies the reference screen and paired-age review calculations |
| `input_templates/` | Header-only CSV templates and an input checklist |
| `examples/synthetic/` | Complete synthetic example with three catchments |
| `examples/dissertation_chapter1/` | Settings used for the Chapter 1 analysis, without dissertation data |
| `docs/COMMAND_GUIDE.md` | Short command reference |
| `run_release_tests.m` | Automated test entry point |

## Prepare the input folders

Create one subfolder for each catchment. Every catchment folder must contain the
four supported filenames.

```text
Catchments/
├── CatchmentA/
│   ├── ZrnPb.csv
│   ├── HblAr.csv
│   ├── ApHeApPb.csv
│   └── ZrnHeZrnPb.csv
└── CatchmentB/
    ├── ZrnPb.csv
    ├── HblAr.csv
    ├── ApHeApPb.csv
    └── ZrnHeZrnPb.csv
```

The required columns are:

| File | Required columns |
|---|---|
| `ZrnPb.csv` | `ZrnPbDate`, `ZrnPb1sigerr` |
| `HblAr.csv` | `HblGrain`, `HblArDate`, `HblAr1sigerr` |
| `ApHeApPb.csv` | `ApGrain`, `ApHeDate`, `ApHe1sigerr`, `ApPbDate`, `ApPb1sigerr` |
| `ZrnHeZrnPb.csv` | `ZrnGrain`, `ZrnHeDate`, `ZrnHe1sigerr`, `ZrnPbDate`, `ZrnPb1sigerr` |

Use the templates in `input_templates/` rather than constructing headers from
memory. The following rules prevent the most common input problems:

- Report ages and absolute uncertainties in Ma.
- Supply analytical uncertainties as 1 sigma values in columns ending in
  `1sigerr`.
- Keep every age paired with the uncertainty from the same analysis and reporting
  convention.
- Use a unique grain or analysis identifier within each thermochronometer file.
- Leave paired U-Pb cells blank when an apatite or zircon He analysis is unpaired.
  Do not enter zero as a missing-value placeholder.
- Investigate zero, negative, missing, or nonnumeric ages and uncertainties before
  running the program.

Older `2sigerr` columns remain accepted so earlier analyses can be reproduced.
The code converts those values to 1 sigma and prints a warning. Do not provide
both 1 sigma and 2 sigma columns for the same age; the program will stop rather
than guess which convention is correct.

The fixed four-file layout does not require every He analysis to have a paired
U-Pb date. An unpaired analysis can be included by leaving its paired U-Pb age and
uncertainty blank. The reference screen will still be applied to the He age, but
paired-age order and interval calculations will be reported as not applicable.

## Run the program

Open MATLAB in the repository folder and add that folder to the path:

```matlab
addpath(pwd)
```

For a first test, run the included synthetic example:

```matlab
run("examples/synthetic/run_example.m")
```

To process your own catchment folders with default settings:

```matlab
run_detrital_pipeline("Catchments", "Output")
```

If the scientific question targets a known age interval, restrict which fitted
component means are eligible for selection:

```matlab
run_detrital_pipeline("Catchments", "Output", ...
    TargetComponentAgeRange=[70 300])
```

This range does not remove ages from the zircon U-Pb distribution. The mixture
model is fit to all valid zircon U-Pb dates, and the range is applied only when
choosing the target component.

## How the zircon reference component is selected

For each catchment, the program performs the following operations:

1. It reads the complete valid age distribution from `ZrnPb.csv`.
2. If analytical uncertainties are available, it propagates them through Monte
   Carlo resampling using the supplied 1 sigma values.
3. It fits Gaussian mixture models with component counts from `K=2` through
   `Kmax=6` by default.
4. It chooses the component count with the lowest Bayesian Information Criterion
   unless the user supplies an override.
5. It identifies the youngest component that satisfies the minimum component
   weight and whose mean falls within `TargetComponentAgeRange`.
6. With the default `gmm_sigma_window` method, it defines the component window as
   the selected mean plus or minus `NSigma` times the selected component standard
   deviation.

The default `NSigma=1` window therefore spans the selected component mean plus or
minus one component standard deviation. The older edge of this window is the
primary reference age used by the filter.

The component is a statistical feature until its geological significance is
evaluated using independent evidence. The program reports it as the selected or
target zircon component rather than automatically calling it a magmatic pulse.

### Manual component-count overrides

Use a component-count override only after inspecting the mixture-model QA figure.
For one K value applied to every catchment:

```matlab
run_detrital_pipeline("Catchments", "Output", K_override=3)
```

For selected catchments:

```matlab
K_map = containers.Map({"WP"}, {3});
run_detrital_pipeline("Catchments", "Output", K_override_map=K_map)
```

The pipeline summary records the selected K, the BIC-preferred K, whether an
override was applied, and the selection method. Report and justify any override
in the associated methods or data-processing description.

### Candidate model-start age

The pipeline also reports a candidate model-start age:

```text
candidate model-start age = component mean + NSigma x component standard deviation
```

This is the older edge of the selected component window. The component standard
deviation describes dispersion within the fitted component; it is not uncertainty
on the component mean or on the candidate model-start age.

## How the thermochronologic dates are screened

Every eligible thermochronologic date is compared with the selected reference
age. Given the measured date and its 1 sigma uncertainty, the program calculates
`P_OlderThanReference`, the probability that the date is older than the reference.

With the default `P_thresh=0.65`:

- `OR1` is assigned when `P_OlderThanReference` is at least 0.90.
- `OR2` is assigned when the probability is at least `P_thresh` but below 0.90.
- `RT` is assigned when the probability is below `P_thresh`.

`OR1` and `OR2` are the only results that recommend exclusion. `RT` means that a
date remains eligible for the stated reference-screening purpose; it does not
validate the analysis for every possible use.

### Paired-age review information

For paired apatite and zircon dates, the program also calculates the U-Pb age
minus the He age and its combined analytical uncertainty. These calculations
provide review information without assigning a mechanism:

- `AOI` identifies a cooling age older than its paired U-Pb age beyond combined
  2 sigma uncertainty.
- `AOU` identifies a nominally older cooling age that overlaps the paired U-Pb
  age within combined 2 sigma uncertainty.
- `SC1` identifies at least 0.90 probability that the crystallization-to-cooling
  interval lies between zero and `Delta`.
- `SC2` identifies a probability between `P_thresh` and 0.90 for that short
  interval.
- `II` identifies missing or nonpositive information required for a calculation.

All input uncertainty columns remain 1 sigma. The AOI and AOU checks use two times
the combined 1 sigma uncertainty when evaluating paired-age overlap.

Review flags never cause exclusion by themselves. If a row has both a review flag
and `Action="exclude"`, the exclusion is independently caused by an OR1 or OR2
reference result. `ActionReason` records both facts.

Paired apatite U-Pb dates are written as separate rows and screened independently
against the zircon reference. Paired zircon U-Pb dates are written as reference
context (`RC`) rather than as candidate model-input ages.

## Main pipeline options

| Option | Default | Meaning |
|---|---:|---|
| `Kmax` | `6` | Largest GMM component count tested during BIC selection |
| `K_override` | `0` | Fixed K for all catchments; zero retains automatic selection |
| `K_override_map` | empty | Catchment-specific K values supplied with a `containers.Map` |
| `BoundsMethod` | `"gmm_sigma_window"` | Method used to define the selected component window |
| `NSigma` | `1.0` | Component standard-deviation multiplier used for the window |
| `Delta` | `8` Ma | Upper bound of the short paired-interval review calculation |
| `P_thresh` | `0.65` | Probability required for OR2 or SC2 assignment |
| `Nmc` | `50` | Monte Carlo draws per zircon U-Pb date when uncertainties are used |
| `TargetComponentAgeRange` | unrestricted | Eligible range for the selected component mean |
| `run_sensitivity` | `false` | Whether to compare the three available reference boundaries |

Changing a parameter changes what the calculation asks. Record nondefault values
and explain why they suit the research question. In particular:

- Lowering `P_thresh` causes more dates to meet the older-than-reference exclusion
  rule and more paired dates to receive short-interval review flags.
- Increasing `Delta` causes more paired dates to meet the short-interval review
  condition. It does not exclude them.
- Increasing `NSigma` moves the primary older-edge reference to an older age and
  widens the selected component window.
- Restricting `TargetComponentAgeRange` changes which fitted component can be
  selected but does not change the dates used to fit the mixture model.

## Output folders and tables

A default run produces:

```text
Output/
├── README.txt
├── pipeline_summary.csv
├── output_summary.csv
├── filter_code_lookup.csv
└── CatchmentA/
    ├── youngest_zircon_component/
    │   ├── youngest_zircon_component_plot.png
    │   └── youngest_zircon_component_summary.csv
    └── filter_output/
        ├── filter_results_full.csv
        ├── filter_results_coded.csv
        ├── model_input_ages.csv
        ├── excluded_ages.csv
        ├── review_flags.csv
        └── output_summary.csv
```

Use the tables as follows:

| Table | Use |
|---|---|
| `filter_results_full.csv` | Review every dated analysis with full descriptions and probabilities |
| `filter_results_coded.csv` | Prepare compact publication or supplementary tables using numeric IDs |
| `model_input_ages.csv` | Obtain dates eligible for downstream modeling under the stated rule |
| `excluded_ages.csv` | Review only dates assigned OR1 or OR2 |
| `review_flags.csv` | Review all non-excluding paired-age and information flags |
| Catchment `output_summary.csv` | Check counts by chronometer and action |
| Root `pipeline_summary.csv` | Review component parameters, K selection, search range, and candidate start ages |
| Root `output_summary.csv` | Compare filtering counts across catchments |
| Root `filter_code_lookup.csv` | Translate codes and numeric identifiers |

The two filter-results tables contain one dated analysis per row. Paired dates
share `GrainID` and `PairID` but occupy separate rows so each chronometer receives
its own result. `ReviewCode="NF"` and `ReviewFlagID=0` mean that no separate review
flag was assigned. Numeric IDs are lookup identifiers, not ranks or scores.

For reliable MATLAB imports of prose-containing output tables, specify the comma
delimiter:

```matlab
T = readtable("filter_results_full.csv", ...
    Delimiter=",", VariableNamingRule="preserve", TextType="string");
```

## Reference-boundary comparison

The optional sensitivity run compares the older edge, midpoint, and younger edge
of the selected component window:

```matlab
run_detrital_pipeline("Catchments", "Output_sensitivity", ...
    run_sensitivity=true)
```

This adds one analysis-level comparison table per catchment and one root-level
summary. It does not create three duplicate output trees. Use this comparison to
show how the selected boundary affects model-input counts; do not treat the three
boundaries as predefined geological interpretations.

## Quality-control procedure

Complete these checks before using `model_input_ages.csv` in a thermal-history
model or publication:

1. Confirm that every input uncertainty uses the intended 1 sigma convention.
2. Inspect each `youngest_zircon_component_plot.png` rather than accepting the
   fitted component automatically.
3. Confirm that the selected component lies within the scientifically intended
   age interval and is supported by a meaningful part of the distribution.
4. Compare `K_used`, `K_bic_selected`, `K_selection_method`, and
   `K_override_applied` in `pipeline_summary.csv`.
5. Record the reason for every manual K override or target-component age range.
6. Inspect `excluded_ages.csv` and `review_flags.csv`, including the measured ages,
   uncertainties, probabilities, and paired dates.
7. Investigate `II`, `AOI`, and `AOU` rows against the original analytical data.
8. If the reference boundary materially affects model-input counts, retain and
   report the sensitivity comparison.
9. Confirm that the final modeling table matches the intended scientific question
   rather than treating the automated recommendation as a universal decision.

## Troubleshooting

### The selected component is an implausibly young minor population

First confirm that the component is real in the QA figure. If the research
question targets a known interval, set `TargetComponentAgeRange` so components
outside that interval remain in the fit but cannot become the target. If the
component structure itself is poor, inspect the BIC curve and fitted densities
before considering a K override. Do not use an override solely to obtain a
preferred filtering outcome.

### No component falls within the target age range

The program stops rather than selecting an out-of-range component. Check the age
units and requested bounds, inspect the zircon U-Pb distribution, and revise the
range only if the scientific target was specified incorrectly or too narrowly.

### The selected component window looks too narrow or too wide

Inspect the fitted component and its standard deviation. Adjust `NSigma` if the
scientific question requires a different fraction of the fitted component. Record
the chosen value. Avoid treating `NSigma` as an uncertainty confidence level on
the component mean.

### The component fit looks unstable

Small or weakly separated age populations can produce unstable mixture models.
Treat the result as provisional, retain the QA figure, inspect neighboring K
solutions, and avoid presenting the automatic selection as definitive. Increasing
`Nmc` can reduce Monte Carlo sampling noise but cannot create information absent
from the measured age distribution.

### Nearly all dates are excluded

Check the 1 sigma uncertainty convention, age units, selected component, reference
age, and `P_thresh`. Run the reference-boundary comparison to determine whether
the result is dominated by the boundary choice. Do not change the threshold only
to reach a desired retained count.

### Many paired dates have short-interval flags

Confirm the pair identifiers and paired ages, then review the measured interval,
combined uncertainty, `Delta`, and `P_ShortInterval`. A short-interval flag is
descriptive review information and does not establish a geological mechanism or
recommend exclusion.

### A catchment is skipped

Confirm that its folder contains all four required filenames with the exact
capitalization shown in this manual. The pipeline reports missing files in the
MATLAB command window and continues to the next catchment.

## Reproducibility record

For each reported run, retain or report:

- DetritalChronFilter version and Git commit
- MATLAB release and toolbox availability
- Input uncertainty convention and any conversion from older files
- `TargetComponentAgeRange`
- `Kmax`, the BIC-selected K, the K used, and any override rationale
- `BoundsMethod` and `NSigma`
- `P_thresh` and `Delta`
- Whether reference-boundary sensitivity was run
- The QA figures, pipeline summary, code lookup, and complete full-results table
- Any manual decisions made after reviewing flagged or excluded dates

The output `README.txt` records the principal settings used in each run. Keep it
with the result tables rather than relying on folder names or memory.

## Dissertation Chapter 1 settings

The included Chapter 1 configuration applies automatic BIC selection to OP and TC
and a reviewed `K=3` override to WP. It also limits eligible target-component
means to 70-300 Ma and enables reference-boundary sensitivity:

```matlab
addpath("examples/dissertation_chapter1")
run_chapter1_configuration("Ch1_Input", "Ch1_Output")
```

The repository contains only the settings. Dissertation data are not distributed.

## Testing

Run the automated release suite from the repository root:

```matlab
run_release_tests
```

The tests use temporary synthetic inputs, exercise automatic and overridden K
selection, validate the full and coded tables, verify output counts, and confirm
that review flags do not cause exclusion.

## Citation and license

Citation metadata are provided in `CITATION.cff`. The software is distributed
under the MIT License. Add an archived-release DOI to the citation metadata after
one is available.
