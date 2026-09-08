# Input-data checklist

Use the four header-only CSV templates in this folder for each catchment.

- Report all ages in Ma.
- Report every analytical uncertainty as an absolute 1σ value in Ma.
- Use column names ending in `1sigerr`. Other uncertainty conventions are not accepted.
- Keep the age and its uncertainty from the same analysis and reporting convention.
- Use a unique grain/analysis identifier within each input file.
- Leave paired-age cells blank for genuinely unpaired observations; do not insert zero as a placeholder.
- Investigate zero, negative, missing, or nonnumeric ages and uncertainties before running the program.
- Preserve the four filenames expected by the current pipeline: `ZrnPb.csv`, `HblAr.csv`, `ApHeApPb.csv`, and `ZrnHeZrnPb.csv`.

Convert uncertainties to 1σ before preparing these files; the code does not convert them automatically.
