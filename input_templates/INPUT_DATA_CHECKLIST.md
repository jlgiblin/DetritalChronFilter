# Input-data checklist

Use the two header-only CSV templates in this folder for each sample.

## Reference distribution

- Name the file `ReferenceDistribution.csv`.
- Include the complete distribution used for component fitting.
- Use one consistent, nonblank `ReferenceSystem` label.
- Give every reference age a unique, nonblank `GrainID`.
- Report `Age_Ma` and a positive absolute `Age_1sigma_Ma`.

## Chronometer data

- Name the file `ChronometerData.csv`.
- Put one dated analysis on each row.
- Use a consistent `Chronometer` label for rows that should be grouped.
- Make the combination of `Chronometer` and `GrainID` unique.
- Report `Age_Ma` and absolute `Age_1sigma_Ma`.
- Leave `PairID` and `PairRole` blank for unpaired analyses.
- For an ordered pair, repeat one `PairID` exactly twice and use one
  `expected_younger` and one `expected_older` role.
- Use two blank roles when retaining the pair link without an order diagnostic.
- Set `UseForModel=false` only for a row supplied as pair/context information.
- Leave `UseForModel` blank or omit the column to default all rows to `true`.

The code does not convert 2σ uncertainties and does not infer scientific
meaning from the system labels.
