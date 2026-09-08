# Changelog

All notable changes to DetritalChronFilter are documented here.

## 0.1.0-rc2 — 2026-09-07

- Replaces four fixed chronometer files with `ReferenceDistribution.csv`
  and a flexible long-format `ChronometerData.csv`.
- Allows any number of user-labelled chronometer systems; labels are used
  for grouping and display, not scientific inference.
- Makes paired relationships optional through `PairID` and `PairRole`.
- Adds `UseForModel` so paired or contextual observations can be retained
  without entering the screened model-input dataset.
- Renames zircon-specific target-component outputs to generic
  `target_component` terminology.
- Adds descriptive post-filter minimum, median, and maximum ages for each
  chronometer to support manual review without imposing temperature order.
- Requires 1-sigma uncertainties directly and removes public 2-sigma
  conversion options.
- Updates templates, synthetic examples, tests, and documentation for the
  generic two-file workflow.

## 0.1.0-rc1 — 2026-09-06

Initial release candidate for the revised public workflow.

- Uses observation-based, mechanism-neutral result terminology.
- Makes older-than-reference probability the only automatic exclusion rule.
- Retains paired-age order and short-interval results as non-excluding review flags.
- Requires absolute 1σ analytical uncertainties for every input dataset.
- Supports automatic BIC selection and explicit per-catchment K overrides.
- Records component search range, component window, K-selection method, and override status.
- Writes aligned full and coded one-date-per-row result tables.
- Separates model inputs, excluded dates, review flags, and output summaries.
- Condenses optional reference-boundary sensitivity into comparison tables.
- Uses `filter_output` and zircon-specific target-component folder names.
- Uses `filter_results` and `output_summary` for the corresponding returned
  MATLAB result fields as well as the written output concepts.
- Adds a complete user manual with input preparation, component-selection,
  filtering, QA, troubleshooting, output, and reproducibility guidance.
