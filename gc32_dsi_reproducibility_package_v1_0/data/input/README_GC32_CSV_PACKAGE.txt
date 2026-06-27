DSI_GC32_CSV_PACKAGE

Strict sample: 32 Galactic globular clusters.
Definition: Gaia EDR3 profile ∩ HST PM dispersion profile ∩ combined RV/PM velocity-dispersion profile.

Important:
- All core kinematic CSV files contain observed rows filtered from the supplied files.
- APOGEE VAC data exist for 23/32 clusters in the supplied data; the remaining 9 clusters are not artificially filled.
- fiducial DSI is intentionally not computed in the CSV package. It should be computed in MATLAB using robust z-score, bootstrap uncertainty and leave-one-observable-out validation.

Main MATLAB file:
load_GC32_DSI_csv_package.m
