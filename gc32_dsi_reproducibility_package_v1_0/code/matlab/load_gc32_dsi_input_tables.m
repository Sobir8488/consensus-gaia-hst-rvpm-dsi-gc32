%% load_gc32_dsi_input_tables.m
% Load the GC32 DSI input CSV tables from the package data/input directory.
clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
packageRoot = fileparts(fileparts(scriptDir));
pkgDir = fullfile(packageRoot,'data','input');
if ~isfolder(pkgDir)
    pkgDir = pwd;
end

T_sample     = readtable(fullfile(pkgDir,'01_gc32_sample_list.csv'));
T_coverage   = readtable(fullfile(pkgDir,'02_gc32_data_coverage_matrix.csv'));
T_hstpm      = readtable(fullfile(pkgDir,'03_gc32_hst_pm_dispersion_profiles_masyr.csv'));
T_gaiaProf   = readtable(fullfile(pkgDir,'04_gc32_gaia_edr3_rotation_dispersion_profiles_masyr.csv'));
T_gaiaMem    = readtable(fullfile(pkgDir,'05_gc32_gaia_member_catalogue_summary.csv'));
T_vel        = readtable(fullfile(pkgDir,'06_gc32_combined_rv_pm_velocity_dispersion_profiles_kms.csv'));
T_hugs       = readtable(fullfile(pkgDir,'07_gc32_hst_hugs_exposure_metadata.csv'));
T_apogeePar  = readtable(fullfile(pkgDir,'08_gc32_apogee_gc_parameters_available23.csv'));
T_apogeeChem = readtable(fullfile(pkgDir,'09_gc32_apogee_abundance_summary_available23.csv'));
T_apogeeStar = readtable(fullfile(pkgDir,'10_gc32_apogee_star_members_available23.csv'));
T_features   = readtable(fullfile(pkgDir,'11_gc32_profile_feature_matrix_derived.csv'));
T_main       = readtable(fullfile(pkgDir,'12_gc32_main_dsi_input_matrix_matlab.csv'));
T_audit      = readtable(fullfile(pkgDir,'13_gc32_missing_data_audit.csv'));
T_structural = readtable(fullfile(pkgDir,'13_gc32_structural_validation.csv'));
T_weights    = readtable(fullfile(pkgDir,'16_gc32_dsi_weights_template.csv'));

assert(height(T_sample)==32, 'Strict sample must contain exactly 32 clusters.');
assert(all(T_coverage.core_dynamics_complete==1), 'Core kinematic coverage is incomplete.');

fprintf('GC32 DSI input package loaded successfully.\n');
fprintf('Clusters: %d\n', height(T_sample));
fprintf('HST PM rows: %d\n', height(T_hstpm));
fprintf('Gaia profile rows: %d\n', height(T_gaiaProf));
fprintf('Velocity profile rows: %d\n', height(T_vel));
fprintf('APOGEE chemistry clusters: %d\n', height(T_apogeeChem));
