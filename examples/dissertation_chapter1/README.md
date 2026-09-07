# Dissertation Chapter 1 configuration

This directory records the settings used for the dissertation Chapter 1 analysis without distributing the underlying data.

- OP and TC use automatic BIC component-count selection.
- WP uses a reviewed `K=3` override.
- The selected component mean must fall between 70 and 300 Ma.
- Component bounds use the mean ± 1 component standard deviation.
- Input analytical uncertainties are absolute 1σ values.
- The older-than-reference and short-interval probability threshold is 0.65.
- The short-interval window is 8 Ma.
- Reference-boundary sensitivity output is enabled.

From the repository root:

```matlab
addpath(pwd)
addpath("examples/dissertation_chapter1")
run_chapter1_configuration("path/to/Ch1_Input", "path/to/Ch1_Output")
```

The input directory must contain catchment subdirectories named `OP`, `TC`, and `WP`, using the standard four-file input layout.
