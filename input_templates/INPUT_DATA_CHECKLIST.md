# Chapter 1 input-data checklist

Use the four header-only CSV templates in this folder for each catchment.

- Report all ages in Ma.
- Report every analytical uncertainty as an absolute 1σ value in Ma.
- Use column names ending in `1sigerr`; do not include a corresponding `2sigerr` column.
- Keep the age and its uncertainty from the same analysis and reporting convention.
- Use a unique grain/analysis identifier within each input file.
- Leave paired-age cells blank for genuinely unpaired observations; do not insert zero as a placeholder.
- Investigate zero, negative, missing, or nonnumeric ages and uncertainties before running the program.
- Preserve the four filenames expected by the current pipeline: `ZrnPb.csv`, `HblAr.csv`, `ApHeApPb.csv`, and `ZrnHeZrnPb.csv`.

The code will stop if both 1σ and 2σ columns are supplied for the same age. Deprecated 2σ columns remain accepted only to reproduce older runs.
