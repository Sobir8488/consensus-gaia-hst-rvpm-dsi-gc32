# MATLAB code

This directory contains MATLAB scripts for the GC32 dynamical-state-index analysis.

## Files

- `run_gc32_dsi_analysis_pipeline.m` - reads `data/input`, recomputes DSI components, uncertainty estimates, null diagnostics, and machine-readable tables.
- `make_gc32_dsi_manuscript_figures.m` - regenerates manuscript-style figures from archived or recomputed tables.
- `load_gc32_dsi_input_tables.m` - lightweight input-table loader and consistency check.

## Recompute the analysis

From the package root in MATLAB:

```matlab
addpath(genpath('code/matlab'));
run('code/matlab/run_gc32_dsi_analysis_pipeline.m');
```

The recomputed files are written to `outputs/recomputed`.

## Regenerate figures from archived tables

```matlab
addpath(genpath('code/matlab'));
make_gc32_dsi_manuscript_figures('outputs/tables_csv', ...
    'data/input/13_gc32_structural_validation.csv', ...
    'outputs/recomputed/figures', false);
```

The archived figures used for submission are in `figures`.
