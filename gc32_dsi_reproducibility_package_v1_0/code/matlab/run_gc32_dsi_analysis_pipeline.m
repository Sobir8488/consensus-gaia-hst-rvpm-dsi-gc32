%% run_gc32_dsi_analysis_pipeline.m
% =========================================================================
% MATLAB analysis pipeline for the GC32 dynamical-state-index package.
%
% Article:
%   "A Gaia--HST--Radial Velocity Based Dynamical-State Index
%    for Galactic Globular Clusters"
%
% Purpose
% -------
% This script reads the 32-cluster CSV input tables, computes the dynamical
% components, constructs fixed-weight, learned-weight and latent-factor DSI
% diagnostics, propagates Monte-Carlo uncertainties, evaluates null models,
% and writes machine-readable result tables.
%
% Input directory
% ---------------
% By default, the script auto-detects the package directory data/input.
% Set dataDir manually below to use another CSV input directory.
%
% Main output directory
% ---------------------
% outputs/recomputed/tables
% outputs/recomputed/figures
% outputs/recomputed/logs
% outputs/recomputed/mat
% =========================================================================

clear; clc; close all;
rng(20260529, 'twister');

%% ============================ USER SETTINGS =============================
dataDir = '';        % empty = auto-detect package data/input directory
runMode = 'publication';

switch lower(runMode)
    case 'fast'
        Nmc = 1000;
        Nnull = 5000;
        Nweight = 2500;
    case 'publication'
        Nmc = 5000;
        Nnull = 50000;
        Nweight = 20000;
    otherwise
        error('Unknown runMode: %s', runMode);
end

% Methodological settings.
% These values are intentionally moderate to keep the 32-GC run feasible in MATLAB Online.
if strcmpi(runMode,'publication')
    Nbayes = 8000; NbayesBurn = 3000; NbayesThin = 5;
    NphysicalNull = min(Nnull,30000);
    NprofileNull = min(3000,Nnull);
else
    Nbayes = 2500; NbayesBurn = 1000; NbayesThin = 5;
    NphysicalNull = min(Nnull,5000);
    NprofileNull = min(1000,Nnull);
end

useParallel = true;       % parfor when available

% Graphics-safe defaults for MATLAB Online.
% The heavy publication run should first write tables/results; figures can be
% generated in a second light run by setting doFigures=true.
doFigures = false;
doPerClusterProfiles = false;
graphicsSafeMode = true;
maxLabelledPoints = 10;

MAD_TO_SIGMA = 1.4826;
zCap = 3.5;
minBinsSlope = 3;
minBinsCurvature = 5;
% stability-weight safeguards
wMax = 0.35;                  % cap any learned component weight to avoid one-component collapse
balancePenaltyStrength = 1.00; % stronger balance penalty in stability-weight learning
% consensus weights are set automatically after adding classical and Bayesian latent diagnostics


% Baseline transparent fixed weights retained for traceability.
wFixed = [0.22, 0.16, 0.14, 0.16, 0.16, 0.16];
wFixed = wFixed ./ sum(wFixed);

componentNames = ["GradientStrength","Anisotropy","RotationSupport", ...
                  "ProfileCurvature","CentralEnhancement","OuterDisturbance"];

warning('off','MATLAB:rankDeficientMatrix');
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');

%% ============================ PATHS =====================================
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir); scriptDir = pwd; end
if isempty(dataDir); dataDir = autoDetectDataDir(scriptDir); end
assert(isfolder(dataDir), 'Data directory not found: %s', dataDir);

packageRoot = fileparts(fileparts(scriptDir));
outDir = fullfile(packageRoot, 'outputs', 'recomputed');
tabDir = fullfile(outDir, 'tables');
figDir = fullfile(outDir, 'figures');
logDir = fullfile(outDir, 'logs');
matDir = fullfile(outDir, 'mat');
mkdirIfNeeded(outDir); mkdirIfNeeded(tabDir); mkdirIfNeeded(figDir);
mkdirIfNeeded(logDir); mkdirIfNeeded(matDir);

diary(fullfile(logDir, 'gc32_dsi_analysis_log.txt'));
fprintf('============================================================\n');
fprintf('GC32 DSI analysis pipeline started\n');
fprintf('Date/time: %s\n', datestr(now));
fprintf('Data dir : %s\n', dataDir);
fprintf('Out dir  : %s\n', outDir);
fprintf('runMode  : %s | Nmc=%d | Nnull=%d | Nweight=%d\n', runMode, Nmc, Nnull, Nweight);
fprintf('============================================================\n');

%% =============================== LOAD DATA ===============================
F = struct();
F.sample   = fullfile(dataDir,'01_gc32_sample_list.csv');
F.coverage = fullfile(dataDir,'02_gc32_data_coverage_matrix.csv');
F.hst      = fullfile(dataDir,'03_gc32_hst_pm_dispersion_profiles_masyr.csv');
F.gaia     = fullfile(dataDir,'04_gc32_gaia_edr3_rotation_dispersion_profiles_masyr.csv');
F.gmem     = fullfile(dataDir,'05_gc32_gaia_member_catalogue_summary.csv');
F.vel      = fullfile(dataDir,'06_gc32_combined_rv_pm_velocity_dispersion_profiles_kms.csv');
F.hugs     = fullfile(dataDir,'07_gc32_hst_hugs_exposure_metadata.csv');
F.apPar    = fullfile(dataDir,'08_gc32_apogee_gc_parameters_available23.csv');
F.apChem   = fullfile(dataDir,'09_gc32_apogee_abundance_summary_available23.csv');
F.apStar   = fullfile(dataDir,'10_gc32_apogee_star_members_available23.csv');
F.main     = fullfile(dataDir,'12_gc32_main_dsi_input_matrix_matlab.csv');
checkFiles(F);

% Optional but recommended external-validation file. It is intentionally NOT
% passed to checkFiles(), because the core DSI can still run without it.
F.structural = fullfile(dataDir,'13_gc32_structural_validation.csv');

T_sample = readtable(F.sample,'PreserveVariableNames',true);
T_cov    = readtable(F.coverage,'PreserveVariableNames',true);
T_hst    = readtable(F.hst,'PreserveVariableNames',true);
T_gaia   = readtable(F.gaia,'PreserveVariableNames',true);
T_gmem   = readtable(F.gmem,'PreserveVariableNames',true);
T_vel    = readtable(F.vel,'PreserveVariableNames',true);
T_hugs   = readtable(F.hugs,'PreserveVariableNames',true);
T_apPar  = readtable(F.apPar,'PreserveVariableNames',true);
T_apChem = readtable(F.apChem,'PreserveVariableNames',true);
T_apStar = readtable(F.apStar,'PreserveVariableNames',true);
T_main   = readtable(F.main,'PreserveVariableNames',true);

if isfile(F.structural)
    T_structural = readtable(F.structural,'PreserveVariableNames',true);
    if ismember("cluster_id", string(T_structural.Properties.VariableNames))
        T_structural.cluster_id = string(T_structural.cluster_id);
    end
else
    T_structural = table();
    warning('Structural validation file not found: %s. Structural external validation will be skipped.', F.structural);
end

tables = {'T_sample','T_cov','T_hst','T_gaia','T_gmem','T_vel','T_hugs','T_apPar','T_apChem','T_apStar','T_main'};
for ii=1:numel(tables)
    eval(sprintf('if ismember("cluster_id", string(%s.Properties.VariableNames)); %s.cluster_id = string(%s.cluster_id); end', tables{ii}, tables{ii}, tables{ii}));
end

clusterID = string(T_sample.cluster_id);
nC = numel(clusterID);
assert(nC==32,'Strict sample must contain 32 clusters.');
assert(all(T_cov.core_dynamics_complete==1),'Core dynamic coverage incomplete.');

fprintf('\nLoaded:\n');
fprintf('  clusters       : %d\n', nC);
fprintf('  HST rows       : %d\n', height(T_hst));
fprintf('  Gaia rows      : %d\n', height(T_gaia));
fprintf('  RV/PM rows     : %d\n', height(T_vel));
fprintf('  APOGEE chem GCs: %d\n', height(T_apChem));

%% ======================== BASELINE FEATURES ==============================
fprintf('\nComputing observed-data dynamic features and physical profile fits.\n');

[Tfeat, Xraw] = computeAllFeatures(clusterID,T_hst,T_gaia,T_vel,minBinsSlope,minBinsCurvature);
chem = computeChemComplexity(clusterID,T_apChem);
Tfeat.chemical_complexity = chem;
Qdata = computeQualityScore(Tfeat);
Tfeat.data_quality_score = Qdata;
TfitQuality = classifyProfileFitQuality(Tfeat);


% EM-PCA imputation instead of basic median.
[Ximp, imputeInfo] = emPCAImpute(Xraw, 1, 100, 1e-6);
if any(~isfinite(Ximp(:)))
    [Ximp, imputeInfoFallback] = medianImpute(Xraw);
    imputeInfo.method = "median_fallback";
    imputeInfo.fallback = imputeInfoFallback;
end

[Zuncap, zLoc, zScale] = robustZ(Ximp, MAD_TO_SIGMA);
Z = winsorizeMatrix(Zuncap, zCap);

%% ======================== LEARNED WEIGHTS ================================
fprintf('\nLearning weights: PCA, entropy and stability optimization.\n');

[wPCA, pcaLoadings, pcaExplained] = pcaWeights(Z);
wEntropy = entropyWeights(Z);

% Monte Carlo components first; stability weights use MC component ensemble.
fprintf('\nRunning covariance-informed profile Monte-Carlo.\n');
useParActual = canUseParallel(useParallel);
[MC_Z, MC_X] = monteCarloComponents(clusterID,T_hst,T_gaia,T_vel,zLoc,zScale,zCap,Nmc, ...
    minBinsSlope,minBinsCurvature,useParActual);

[DSI_fixed, DSI_fixed_MC, DSI_fixed_BC, DSI_fixed_err, DSI_fixed_bias] = ...
    biasCorrectedIndex(Z,MC_Z,wFixed);

[DSI_pca, DSI_pca_MC, DSI_pca_BC, DSI_pca_err, DSI_pca_bias] = ...
    biasCorrectedIndex(Z,MC_Z,wPCA);

[DSI_entropy, DSI_entropy_MC, DSI_entropy_BC, DSI_entropy_err, DSI_entropy_bias] = ...
    biasCorrectedIndex(Z,MC_Z,wEntropy);

[wStability, stabilityInfo] = learnStabilityWeights(Z,MC_Z,Nweight,wMax,balancePenaltyStrength);
[DSI_stability, DSI_stability_MC, DSI_stability_BC, DSI_stability_err, DSI_stability_bias] = ...
    biasCorrectedIndex(Z,MC_Z,wStability);

% Quality-aware conservative score
Qshrink = sqrt(max(0,min(1,Qdata)));
DSI_stability_QA = DSI_stability_BC .* Qshrink;

%% ======================== LATENT VARIABLE MODELS ==========================
fprintf('\nEstimating classical and Bayesian latent dynamical complexity.\n');

% Classical one-factor score retained for backward comparability.
[theta_latent, lambda_latent, latentInfo] = latentOneFactor(Z);
thetaMC = nan(nC,Nmc);
for b=1:Nmc
    thetaMC(:,b) = MC_Z(:,:,b) * lambda_latent(:);
end
thetaMC_med = median(thetaMC,2,'omitnan');
thetaMC_p16 = prctileLocal(thetaMC,16,2);
thetaMC_p84 = prctileLocal(thetaMC,84,2);
thetaMC_err = (thetaMC_p84-thetaMC_p16)./2;
theta_bias = thetaMC_med - theta_latent;
theta_latent_BC = theta_latent - theta_bias;
theta_latent_QA = theta_latent_BC .* Qshrink;

% v5.0: hierarchical Bayesian latent factor model with measurement scatter.
% Model: z_ij ~ Normal(alpha_j + lambda_j*theta_i, sigma_int_j^2 + u_ij),
% where u_ij is estimated from the covariance-informed MC component ensemble.
[theta_bayes_BC, theta_bayes_err, lambda_bayes, lambda_bayes_err, ...
    alpha_bayes, sigma_int_bayes, bayesInfo] = bayesianLatentDSI(Z,MC_Z,Nbayes,NbayesBurn,NbayesThin);
theta_bayes_QA = theta_bayes_BC .* Qshrink;
wBayes = abs(lambda_bayes(:));
if sum(wBayes)<=0 || any(~isfinite(wBayes)), wBayes = ones(nK,1)./nK; else, wBayes = wBayes./sum(wBayes); end
[DSI_bayesLoad, DSI_bayesLoad_MC, DSI_bayesLoad_BC, DSI_bayesLoad_err, DSI_bayesLoad_bias] = ...
    biasCorrectedIndex(Z,MC_Z,wBayes);

% v5 consensus index: agreement across fixed, PCA, entropy, stability,
% classical latent and Bayesian latent diagnostics.
ConsensusMatrix = [ ...
    standardizeVector(DSI_fixed_BC), ...
    standardizeVector(DSI_pca_BC), ...
    standardizeVector(DSI_entropy_BC), ...
    standardizeVector(DSI_stability_BC), ...
    standardizeVector(theta_latent_BC), ...
    standardizeVector(theta_bayes_BC)];
consensusWeights = ones(1,size(ConsensusMatrix,2))./size(ConsensusMatrix,2);
DSI_consensus_BC = ConsensusMatrix * consensusWeights(:);
DSI_consensus_QA = DSI_consensus_BC .* Qshrink;

%% ======================== NULL AND PHYSICAL-NULL MODELS ===================
fprintf('\nRunning empirical, physical-proxy, rotation-free and profile-shuffle null simulations.\n');

% Empirical incoherent-component null: independently permutes component values
% across clusters. This is no longer described as a literal physical relaxed
% cluster model; it is a coherence-destroying empirical null.
[nullDSI, pNull] = relaxedNullPermutation(Z,wStability,DSI_stability_BC,Nnull);
qNull = bhFDR(pNull);
nullSignificance = classifyNullSignificance(pNull,qNull);

[nullTheta, pNullTheta] = relaxedNullPermutation(Z,lambda_latent(:).'/sum(abs(lambda_latent)),theta_latent_BC,Nnull);
qNullTheta = bhFDR(pNullTheta);

% v5.0 Bayesian-latent null using empirical component shuffling.
[nullThetaBayes, pNullThetaBayes] = relaxedNullPermutation(Z,lambda_bayes(:).'/sum(abs(lambda_bayes)),theta_bayes_BC,Nnull);
qNullThetaBayes = bhFDR(pNullThetaBayes);

% v5.0 consensus null: combines stability DSI, classical latent and Bayesian latent analogues.
nullConsensus = nan(nC, Nnull);
for bb = 1:Nnull
    nullConsensus(:,bb) = mean([ ...
        standardizeVector(nullDSI(:,bb)), ...
        standardizeVector(nullTheta(:,bb)), ...
        standardizeVector(nullThetaBayes(:,bb))], 2, 'omitnan');
end
pNullConsensus = nan(nC,1);
for ii = 1:nC
    pNullConsensus(ii) = (1 + sum(nullConsensus(ii,:) >= DSI_consensus_BC(ii))) / (Nnull + 1);
end
qNullConsensus = bhFDR(pNullConsensus);

% v5.0 physical relaxed proxy null in component space. This is not an N-body
% model; it is a conservative low-rotation/low-curvature/regular-profile null.
[physNullScores, pNullPhysical, physNullInfo] = physicalRelaxedNullProxy(Z,wStability,DSI_stability_BC,NphysicalNull);
qNullPhysical = bhFDR(pNullPhysical);

% v5.0 rotation-free null. Tests whether high DSI survives when the rotation
% support component is drawn from the low-rotation empirical regime.
[rotFreeScores, pNullRotationFree, rotationFreeStress] = rotationFreeNull(Z,wStability,DSI_stability_BC,NphysicalNull);
qNullRotationFree = bhFDR(pNullRotationFree);

% v5.0 profile-shuffle null. Destroys radial ordering inside each observed
% profile while preserving the measured values and error scale.
[profileShuffleScores, pNullProfileShuffle] = profileShuffleNull(clusterID,T_hst,T_gaia,T_vel, ...
    zLoc,zScale,zCap,wStability,DSI_stability_BC,NprofileNull,minBinsSlope,minBinsCurvature,useParActual);
qNullProfileShuffle = bhFDR(pNullProfileShuffle);

% Conservative combined null: a cluster is treated as robust only if it is
% unusual relative to all null families. This avoids overclaiming from any one null.
pNullCombinedConservative = max([pNull, pNullPhysical, pNullRotationFree, pNullProfileShuffle],[],2,'omitnan');
qNullCombinedConservative = bhFDR(pNullCombinedConservative);

%% ======================== MIXTURE CLASSIFICATION =========================
fprintf('\nFitting 3-state Gaussian mixture to primary DSI and latent score.\n');

Ymix = [standardizeVector(DSI_stability_BC), standardizeVector(theta_latent_BC)];
TgmmSelect = gmmModelSelection(Ymix, 2:4, 300, 1e-7);
GMM = fitGMM_EM(Ymix,3,300,1e-7);
mixClass = classifyMixture(GMM,Ymix);
mixProb = GMM.gamma;
mixEntropy = -sum(mixProb.*log(max(mixProb,eps)),2);

%% ======================== LOO + WEIGHT ROBUSTNESS ========================
fprintf('\nLeave-one-component-out and robustness diagnostics.\n');

nK = numel(componentNames);
DSI_loo = nan(nC,nK);
for k=1:nK
    wk = wStability; wk(k)=0; wk=wk./sum(wk);
    DSI_loo(:,k) = Z*wk(:);
end
looDelta = DSI_loo - DSI_stability;
looMaxAbsDelta = max(abs(looDelta),[],2,'omitnan');

rankMC = nan(nC,Nmc);
for b=1:Nmc
    rankMC(:,b) = tiedRanking(-DSI_stability_MC(:,b));
end
rankMC_sigma = std(rankMC,0,2,'omitnan');
rankPrimary = tiedRanking(-DSI_stability_BC);

%% ======================== EXTERNAL VALIDATION ============================
fprintf('\nExternal validation and correlations.\n');

Tjoin = joinValidationTables(clusterID,T_main,T_gmem,T_apPar,T_apChem,Tfeat);

% Merge optional Harris/Baumgardt structural-validation table when present.
% These quantities are external validation variables only; they are NOT used
% in DSI construction, weight learning, latent modelling or null models.
if ~isempty(T_structural) && height(T_structural)>0 && ismember("cluster_id", string(T_structural.Properties.VariableNames))
    T_structural = unique(T_structural,'rows');
    Tjoin = outerjoin(Tjoin,T_structural, ...
        'Keys','cluster_id', ...
        'MergeKeys',true, ...
        'Type','left');
end

[Tjoin, TstructMap] = addStructuralDynamicalValidationVariables(Tjoin, clusterID);
candidateVars = ["log10_mass","mass","Rgc","Fe_H_med","Fe_H_sig","dist", ...
    "RV_sig","PM_sig","N_member_p50","N_member_p90", ...
    "pmra_member_sigma_masyr","pmdec_member_sigma_masyr", ...
    "qflag_pass_fraction","N_apogee_total","FE_H_sigma","N_FE_sigma", ...
    "AL_FE_sigma","MG_FE_sigma","NA_FE_sigma","chemical_complexity", ...
    "hst_eta","gaia_eta","vel_eta","lambda_R_gaia", ...
    "concentration_c","log10_trc_yr","log10_trh_yr", ...
    "trc_yr","trh_yr","dyn_age_rc","dyn_age_rh", ...
    "baumgardt_mass_msun","baumgardt_Rgc_kpc","baumgardt_dist_kpc"];

Tcorr_DSI = computeCorrelations(Tjoin,clusterID,DSI_stability_BC,candidateVars);
Tcorr_theta = computeCorrelations(Tjoin,clusterID,theta_latent_BC,candidateVars);
if height(Tcorr_DSI)>0
    Tcorr_DSI.BH_q_spearman = bhFDR(Tcorr_DSI.spearman_p_approx);
    Tcorr_DSI = sortrows(Tcorr_DSI,'abs_spearman_rho','descend');
end
if height(Tcorr_theta)>0
    Tcorr_theta.BH_q_spearman = bhFDR(Tcorr_theta.spearman_p_approx);
    Tcorr_theta = sortrows(Tcorr_theta,'abs_spearman_rho','descend');
end

% Focused validation: DSI versus relaxation time, concentration and
% core-collapse status. These tests are written even if the corresponding
% structural variables are partially unavailable.
structuralVars = ["concentration_c","log10_trc_yr","log10_trh_yr", ...
    "trc_yr","trh_yr","dyn_age_rc","dyn_age_rh"];

Tcorr_struct_DSI = computeCorrelations(Tjoin,clusterID,DSI_consensus_BC,structuralVars);
if height(Tcorr_struct_DSI)>0
    Tcorr_struct_DSI.BH_q_spearman = bhFDR(Tcorr_struct_DSI.spearman_p_approx);
    Tcorr_struct_DSI = sortrows(Tcorr_struct_DSI,'abs_spearman_rho','descend');
end

Tcorr_struct_stability = computeCorrelations(Tjoin,clusterID,DSI_stability_BC,structuralVars);
if height(Tcorr_struct_stability)>0
    Tcorr_struct_stability.BH_q_spearman = bhFDR(Tcorr_struct_stability.spearman_p_approx);
    Tcorr_struct_stability = sortrows(Tcorr_struct_stability,'abs_spearman_rho','descend');
end

Tcc_consensus = compareCoreCollapseGroups(Tjoin,clusterID,DSI_consensus_BC, ...
    "DSI_consensus_BC");
Tcc_stability = compareCoreCollapseGroups(Tjoin,clusterID,DSI_stability_BC, ...
    "DSI_stability_BC");

writetable(TstructMap, fullfile(tabDir,'GC32_v5_1_structural_validation_column_map.csv'));
writetable(Tcorr_struct_DSI, fullfile(tabDir,'GC32_v5_1_DSI_vs_relaxation_concentration.csv'));
writetable(Tcorr_struct_stability, fullfile(tabDir,'GC32_v5_1_stabilityDSI_vs_relaxation_concentration.csv'));
writetable(Tcc_consensus, fullfile(tabDir,'GC32_v5_1_core_collapse_group_test_consensusDSI.csv'));
writetable(Tcc_stability, fullfile(tabDir,'GC32_v5_1_core_collapse_group_test_stabilityDSI.csv'));

%% ======================== PUBLICATION READINESS ==========================
readiness = computeReadiness(Qdata,DSI_stability_err,DSI_stability_bias,rankMC_sigma,looMaxAbsDelta,pNull);
readinessClass = classifyReadiness(readiness);

interpretFlag = strings(nC,1);
for i=1:nC
    if Qdata(i)<0.70
        interpretFlag(i)="caution_low_coverage";
    elseif abs(DSI_stability_bias(i))>max(0.35,DSI_stability_err(i))
        interpretFlag(i)="caution_MC_bias";
    elseif rankMC_sigma(i)>4
        interpretFlag(i)="caution_rank_instability";
    elseif looMaxAbsDelta(i)>0.75
        interpretFlag(i)="caution_component_sensitive";
    elseif qNull(i)<0.10
        interpretFlag(i)="robust_null_significant";
    else
        interpretFlag(i)="robust";
    end
end

%% ======================== OUTPUT TABLES =================================
fprintf('\nWriting publication tables.\n');

Tresults = table();
Tresults.cluster_id = clusterID;
Tresults.DSI_fixed_BC = DSI_fixed_BC;
Tresults.DSI_pca_BC = DSI_pca_BC;
Tresults.DSI_entropy_BC = DSI_entropy_BC;
Tresults.DSI_stability_BC = DSI_stability_BC;
Tresults.DSI_stability_QA = DSI_stability_QA;
Tresults.DSI_stability_err = DSI_stability_err;
Tresults.DSI_stability_bias = DSI_stability_bias;
Tresults.theta_latent_BC = theta_latent_BC;
Tresults.theta_latent_QA = theta_latent_QA;
Tresults.theta_bayes_BC = theta_bayes_BC;
Tresults.theta_bayes_QA = theta_bayes_QA;
Tresults.theta_bayes_err = theta_bayes_err;
Tresults.DSI_bayesLoad_BC = DSI_bayesLoad_BC;
Tresults.DSI_bayesLoad_err = DSI_bayesLoad_err;
Tresults.DSI_consensus_BC = DSI_consensus_BC;
Tresults.DSI_consensus_QA = DSI_consensus_QA;
Tresults.theta_latent_err = thetaMC_err;
Tresults.p_null_DSI = pNull;
Tresults.q_null_DSI = qNull;
Tresults.p_null_theta = pNullTheta;
Tresults.q_null_theta = qNullTheta;
Tresults.p_null_theta_bayes = pNullThetaBayes;
Tresults.q_null_theta_bayes = qNullThetaBayes;
Tresults.p_null_consensus = pNullConsensus;
Tresults.q_null_consensus = qNullConsensus;
Tresults.p_null_physical_proxy = pNullPhysical;
Tresults.q_null_physical_proxy = qNullPhysical;
Tresults.p_null_rotation_free = pNullRotationFree;
Tresults.q_null_rotation_free = qNullRotationFree;
Tresults.p_null_profile_shuffle = pNullProfileShuffle;
Tresults.q_null_profile_shuffle = qNullProfileShuffle;
Tresults.p_null_combined_conservative = pNullCombinedConservative;
Tresults.q_null_combined_conservative = qNullCombinedConservative;
Tresults.null_significance = nullSignificance;
Tresults.mixture_class = mixClass;
Tresults.mix_prob_relaxed = mixProb(:,1);
Tresults.mix_prob_intermediate = mixProb(:,2);
Tresults.mix_prob_disturbed = mixProb(:,3);
Tresults.mix_entropy = mixEntropy;
Tresults.rank_primary = rankPrimary;
Tresults.rank_MC_sigma = rankMC_sigma;
Tresults.loo_max_abs_delta = looMaxAbsDelta;
Tresults.data_quality_score = Qdata;
Tresults.publication_readiness = readiness;
Tresults.publication_readiness_class = readinessClass;
Tresults.interpretation_flag = interpretFlag;
Tresults.hst_n_bins = Tfeat.hst_n_bins;
Tresults.gaia_n_bins = Tfeat.gaia_n_bins;
Tresults.vel_n_bins = Tfeat.vel_n_bins;
Tresults.hst_fit_flag = TfitQuality.hst_fit_flag;
Tresults.gaia_fit_flag = TfitQuality.gaia_fit_flag;
Tresults.vel_fit_flag = TfitQuality.vel_fit_flag;
Tresults.n_problematic_profile_fits = TfitQuality.n_problematic_profile_fits;
Tresults.lambda_R_gaia = Tfeat.lambda_R_gaia;
Tresults.rotation_support_robust = Tfeat.rotation_support_robust;
Tresults.chemistry_available = isfinite(chem);

% Add dominant component from stability contributions
contrib = bsxfun(@times, Z, wStability(:)');
[~,domIdx] = max(abs(contrib),[],2);
dom = strings(nC,1);
for i=1:nC; dom(i)=componentNames(domIdx(i)); end
Tresults.dominant_component = dom;

T_rank_primary = sortrows(Tresults,'DSI_stability_BC','descend');
T_rank_QA = sortrows(Tresults,'DSI_stability_QA','descend');
T_rank_theta = sortrows(Tresults,'theta_latent_BC','descend');
T_rank_consensus = sortrows(Tresults,'DSI_consensus_BC','descend');
T_rank_consensus_QA = sortrows(Tresults,'DSI_consensus_QA','descend');

% Component tables
Traw = table();
Traw.cluster_id = clusterID;
for k=1:nK; Traw.(componentNames(k)) = Xraw(:,k); end
Traw.chemical_complexity = chem;
Traw.data_quality_score = Qdata;

Tz = array2table(Z,'VariableNames',componentNames);
Tz.cluster_id=clusterID; Tz=movevars(Tz,'cluster_id','Before',1);

TzUncap = array2table(Zuncap,'VariableNames',componentNames);
TzUncap.cluster_id=clusterID; TzUncap=movevars(TzUncap,'cluster_id','Before',1);

Tweights = table(componentNames(:),wFixed(:),wPCA(:),wEntropy(:),wStability(:),wBayes(:), ...
    lambda_latent(:),lambda_bayes(:),lambda_bayes_err(:), ...
    'VariableNames',{'component','w_fixed','w_PCA','w_entropy','w_stability','w_bayesian_loading','latent_loading_classical','latent_loading_bayes','latent_loading_bayes_err'});

Tprofile = Tfeat;

Tlatent = table();
Tlatent.component = componentNames(:);
Tlatent.loading_classical = lambda_latent(:);
Tlatent.loading_bayes = lambda_bayes(:);
Tlatent.loading_bayes_err = lambda_bayes_err(:);
Tlatent.alpha_bayes = alpha_bayes(:);
Tlatent.sigma_int_bayes = sigma_int_bayes(:);
Tlatent.pca_loading_PC1 = pcaLoadings(:);
Tlatent.pca_explained_PC1_percent = repmat(pcaExplained(1),nK,1);

Tmc = table(clusterID,DSI_stability,DSI_stability_BC,DSI_stability_err,DSI_stability_bias, ...
    theta_latent,theta_latent_BC,thetaMC_err,theta_bias,theta_bayes_BC,theta_bayes_err, ...
    DSI_bayesLoad_BC,DSI_bayesLoad_err, ...
    'VariableNames',{'cluster_id','DSI_stability_raw','DSI_stability_BC','DSI_stability_err','DSI_stability_bias','theta_raw','theta_BC','theta_err','theta_bias','theta_bayes_BC','theta_bayes_err','DSI_bayesLoad_BC','DSI_bayesLoad_err'});

Tnull = table(clusterID,pNull,qNull,pNullTheta,qNullTheta,pNullThetaBayes,qNullThetaBayes, ...
    pNullConsensus,qNullConsensus,pNullPhysical,qNullPhysical,pNullRotationFree,qNullRotationFree, ...
    pNullProfileShuffle,qNullProfileShuffle,pNullCombinedConservative,qNullCombinedConservative,nullSignificance, ...
    'VariableNames',{'cluster_id','p_null_DSI','q_null_DSI','p_null_theta','q_null_theta','p_null_theta_bayes','q_null_theta_bayes','p_null_consensus','q_null_consensus','p_null_physical_proxy','q_null_physical_proxy','p_null_rotation_free','q_null_rotation_free','p_null_profile_shuffle','q_null_profile_shuffle','p_null_combined_conservative','q_null_combined_conservative','empirical_null_significance'});

Tmix = table(clusterID,mixClass,mixProb(:,1),mixProb(:,2),mixProb(:,3),mixEntropy, ...
    'VariableNames',{'cluster_id','mixture_class','P_relaxed','P_intermediate','P_disturbed','mixture_entropy'});

Tbayes = table();
Tbayes.cluster_id = clusterID;
Tbayes.theta_bayes_BC = theta_bayes_BC;
Tbayes.theta_bayes_err = theta_bayes_err;
Tbayes.theta_bayes_QA = theta_bayes_QA;
Tbayes.bayes_latent_acceptance_proxy = repmat(bayesInfo.nPosteriorSamples,nC,1);

TphysicalNull = table(clusterID,pNullPhysical,qNullPhysical,pNullRotationFree,qNullRotationFree, ...
    pNullProfileShuffle,qNullProfileShuffle,pNullCombinedConservative,qNullCombinedConservative,rotationFreeStress, ...
    'VariableNames',{'cluster_id','p_physical_proxy','q_physical_proxy','p_rotation_free','q_rotation_free','p_profile_shuffle','q_profile_shuffle','p_combined_conservative','q_combined_conservative','rotation_free_DSI_drop'});

writetable(T_rank_primary, fullfile(tabDir,'gc32_dsi_primary_ranked_dsi_stability_bc.csv'));
writetable(T_rank_QA, fullfile(tabDir,'gc32_dsi_quality_aware_ranked_dsi_stability_qa.csv'));
writetable(T_rank_theta, fullfile(tabDir,'gc32_dsi_classical_latent_ranked_theta_bc.csv'));
writetable(T_rank_consensus, fullfile(tabDir,'gc32_dsi_consensus_ranked_bc.csv'));
writetable(T_rank_consensus_QA, fullfile(tabDir,'gc32_dsi_consensus_ranked_qa.csv'));
writetable(Tresults, fullfile(tabDir,'gc32_dsi_all_results_original_order.csv'));
writetable(Traw, fullfile(tabDir,'gc32_dsi_raw_components.csv'));
writetable(Tz, fullfile(tabDir,'gc32_dsi_robust_z_components_capped.csv'));
writetable(TzUncap, fullfile(tabDir,'gc32_dsi_robust_z_components_uncapped.csv'));
writetable(Tweights, fullfile(tabDir,'gc32_dsi_learned_weights_and_latent_loadings.csv'));
writetable(Tprofile, fullfile(tabDir,'gc32_dsi_profile_physical_fit_diagnostics.csv'));
writetable(TfitQuality, fullfile(tabDir,'gc32_dsi_profile_fit_quality_flags.csv'));
writetable(Tlatent, fullfile(tabDir,'gc32_dsi_latent_factor_model.csv'));
writetable(Tmc, fullfile(tabDir,'gc32_dsi_mc_bias_uncertainty_summary.csv'));
writetable(Tnull, fullfile(tabDir,'gc32_dsi_null_family_significance.csv'));
writetable(Tmix, fullfile(tabDir,'gc32_dsi_mixture_model_classes.csv'));
writetable(Tbayes, fullfile(tabDir,'gc32_dsi_bayesian_latent_scores.csv'));
writetable(TphysicalNull, fullfile(tabDir,'gc32_dsi_component_space_relaxed_proxy_null.csv'));
writetable(TgmmSelect, fullfile(tabDir,'gc32_dsi_gmm_model_selection.csv'));
writetable(Tcorr_DSI, fullfile(tabDir,'gc32_dsi_external_validation_correlations_dsi.csv'));
writetable(Tcorr_theta, fullfile(tabDir,'gc32_dsi_external_validation_correlations_theta.csv'));

% Short tables for manuscript
writetable(T_rank_primary(1:min(12,nC),:), fullfile(tabDir,'gc32_dsi_top12_primary.csv'));
writetable(T_rank_consensus(1:min(12,nC),:), fullfile(tabDir,'gc32_dsi_top12_consensus.csv'));
writetable(T_rank_QA(1:min(12,nC),:), fullfile(tabDir,'gc32_dsi_top12_quality_aware.csv'));
writetable(T_rank_primary(max(1,nC-11):nC,:), fullfile(tabDir,'gc32_dsi_bottom12_primary.csv'));

%% ======================== FIGURES =======================================
if doFigures
    fprintf('\nGenerating figures.\n');
    if graphicsSafeMode
        set(0,'DefaultFigureVisible','off');
    end
    order = orderBy(DSI_stability_BC,'descend');

    fig=figure('Color','w','Position',[100 100 1350 650]);
    bar(1:nC,DSI_stability_BC(order),'FaceAlpha',0.85); hold on;
    errorbar(1:nC,DSI_stability_BC(order),DSI_stability_err(order),'k.','LineWidth',1);
    xticks(1:nC); xticklabels(clusterID(order)); xtickangle(60);
    ylabel('DSI_{stability,BC}');
    title('Primary learned-weight bias-corrected DSI ranking');
    grid on; box on; saveFig(fig,figDir,'fig01_primary_dsi_stability_bc_ranked'); close(fig);

    fig=figure('Color','w','Position',[100 100 900 650]);
    scatter(DSI_stability_BC,theta_latent_BC,90,-log10(max(pNull,1e-12)),'filled'); hold on;
    addLabelsTopN(DSI_stability_BC,theta_latent_BC,clusterID,DSI_stability_BC,maxLabelledPoints);
    xlabel('DSI_{stability,BC}'); ylabel('\theta_{latent,BC}');
    title('Composite DSI versus latent dynamical complexity');
    cb=colorbar; cb.Label.String='-log_{10}(p_{null})';
    grid on; box on; saveFig(fig,figDir,'fig02_dsi_vs_latent_theta'); close(fig);

    fig=figure('Color','w','Position',[100 100 1050 740]);
    imagesc(Z(order,:)); colormap(parula); cb=colorbar; cb.Label.String='robust capped z';
    yticks(1:nC); yticklabels(clusterID(order));
    xticks(1:nK); xticklabels(componentNames); xtickangle(45);
    title('DSI component matrix sorted by primary DSI');
    saveFig(fig,figDir,'fig03_component_heatmap'); close(fig);

    fig=figure('Color','w','Position',[100 100 1050 740]);
    imagesc(contrib(order,:)); colormap(parula); cb=colorbar; cb.Label.String='weighted contribution';
    yticks(1:nC); yticklabels(clusterID(order));
    xticks(1:nK); xticklabels(componentNames); xtickangle(45);
    title('Stability-weighted component contributions');
    saveFig(fig,figDir,'fig04_component_contributions'); close(fig);

    fig=figure('Color','w','Position',[100 100 900 650]);
    scatter(Qdata,DSI_stability_BC,90,readiness,'filled'); hold on;
    addLabelsTopN(Qdata,DSI_stability_BC,clusterID,DSI_stability_BC,maxLabelledPoints);
    xlabel('Data quality score'); ylabel('DSI_{stability,BC}');
    title('DSI versus data quality and publication readiness');
    cb=colorbar; cb.Label.String='readiness score';
    grid on; box on; saveFig(fig,figDir,'fig05_quality_readiness'); close(fig);

    fig=figure('Color','w','Position',[100 100 900 650]);
    scatter(Tfeat.lambda_R_gaia,DSI_stability_BC,90,Tfeat.rotation_support_robust,'filled'); hold on;
    addLabelsTopN(Tfeat.lambda_R_gaia,DSI_stability_BC,clusterID,DSI_stability_BC,maxLabelledPoints);
    xlabel('\lambda_R-like Gaia spin proxy'); ylabel('DSI_{stability,BC}');
    title('Spin proxy and dynamical complexity');
    cb=colorbar; cb.Label.String='P90(|V|)/median(\sigma)';
    grid on; box on; saveFig(fig,figDir,'fig06_lambdaR_rotation'); close(fig);

    fig=figure('Color','w','Position',[100 100 900 650]);
    scatter(DSI_stability_BC,-log10(max(pNull,1e-12)),90,Qdata,'filled'); hold on;
    addLabelsTopN(DSI_stability_BC,-log10(max(pNull,1e-12)),clusterID,DSI_stability_BC,maxLabelledPoints);
    xlabel('DSI_{stability,BC}'); ylabel('-log_{10}(p_{null})');
    title('Relaxed-null significance of DSI');
    cb=colorbar; cb.Label.String='Q_{data}';
    grid on; box on; saveFig(fig,figDir,'fig07_null_significance'); close(fig);

    fig=figure('Color','w','Position',[100 100 900 650]);
    hold on;
    clsU = unique(mixClass, 'stable');
    markerList = {'o','s','^','d','v','>'};
    for kk = 1:numel(clsU)
        idxCls = mixClass == clsU(kk);
        scatter(DSI_stability_BC(idxCls), theta_latent_BC(idxCls), 85, markerList{min(kk,numel(markerList))}, 'filled', ...
            'DisplayName', char(clsU(kk)));
    end
    addLabelsTopN(DSI_stability_BC,theta_latent_BC,clusterID,DSI_stability_BC,maxLabelledPoints);
    xlabel('DSI_{stability,BC}'); ylabel('\theta_{latent,BC}');
    title('GMM-inferred dynamical-state classes');
    legend('Location','best');
    grid on; box on; saveFig(fig,figDir,'fig08_mixture_classes'); close(fig);

    fig=figure('Color','w','Position',[100 100 900 650]);
    scatter(Tfeat.hst_eta, Tfeat.vel_eta, 90, DSI_stability_BC, 'filled'); hold on;
    addLabelsTopN(Tfeat.hst_eta,Tfeat.vel_eta,clusterID,DSI_stability_BC,maxLabelledPoints);
    xlabel('\eta_{HST} profile'); ylabel('\eta_{RV/PM} profile');
    title('Physical profile-shape diagnostics');
    cb=colorbar; cb.Label.String='DSI_{stability,BC}';
    grid on; box on; saveFig(fig,figDir,'fig09_physical_profile_eta'); close(fig);

    scatterValidation(Tjoin,clusterID,DSI_stability_BC,"Rgc","R_{GC}",figDir,"fig10_dsi_vs_rgc");
    scatterValidation(Tjoin,clusterID,DSI_stability_BC,"Fe_H_med","[Fe/H]",figDir,"fig11_dsi_vs_feh");
    scatterValidation(Tjoin,clusterID,DSI_stability_BC,"log10_mass","log_{10}(M/M_\odot)",figDir,"fig12_dsi_vs_mass");
    scatterValidation(Tjoin,clusterID,DSI_stability_BC,"dist","Distance",figDir,"fig13_dsi_vs_distance");

    % Structural/evolutionary external-validation figures.
    scatterValidation(Tjoin,clusterID,DSI_consensus_BC,"log10_trh_yr", ...
        "log_{10}(t_{rh}/{\rm yr})",figDir,"fig14_v5_1_DSI_vs_half_mass_relaxation_time");
    scatterValidation(Tjoin,clusterID,DSI_consensus_BC,"log10_trc_yr", ...
        "log_{10}(t_{rc}/{\rm yr})",figDir,"fig15_v5_1_DSI_vs_core_relaxation_time");
    scatterValidation(Tjoin,clusterID,DSI_consensus_BC,"concentration_c", ...
        "Concentration c",figDir,"fig16_v5_1_DSI_vs_concentration");
    plotCoreCollapseBox(Tjoin,clusterID,DSI_consensus_BC,figDir, ...
        "fig17_v5_1_DSI_by_core_collapse_status");

    if doPerClusterProfiles
        % Select top objects from primary and latent rankings:
        ordTheta = orderBy(theta_latent_BC,'descend');
        selected = unique([clusterID(order(1:min(6,nC))); clusterID(ordTheta(1:min(6,nC)))], 'stable');
        for i=1:numel(selected)
            fig=plotProfilePanel(selected(i),T_hst,T_gaia,T_vel);
            saveFig(fig,figDir,"profile_"+selected(i));
            close(fig);
        end
    end
end

%% ======================== SAVE MAT + SUMMARY =============================
saveFile = fullfile(matDir,'gc32_dsi_workspace.mat');
try
    save(saveFile, 'Tresults','T_rank_primary','T_rank_QA','T_rank_theta','T_rank_consensus', ...
        'Traw','Tz','Tweights','Tprofile','TfitQuality','Tlatent','Tmc','Tnull','Tmix', ...
        'TgmmSelect','Tcorr_DSI','Tcorr_theta','TstructMap','Tcorr_struct_DSI','Tcorr_struct_stability','Tcc_consensus','Tcc_stability','Z','Zuncap','Xraw','Ximp', ...
        'DSI_fixed_BC','DSI_pca_BC','DSI_entropy_BC','DSI_stability_BC','DSI_consensus_BC', ...
        'theta_latent_BC','theta_bayes_BC','theta_bayes_err','lambda_bayes','lambda_bayes_err','MC_Z','wFixed','wPCA','wEntropy','wStability','wBayes','lambda_latent', ...
        'Qdata','pNull','pNullTheta','pNullThetaBayes','pNullConsensus','pNullPhysical','pNullRotationFree','pNullProfileShuffle','pNullCombinedConservative','componentNames','bayesInfo','physNullInfo','-v7.3');
catch saveErr
    warning('v4.2:saveWarning','Large MAT save failed: %s. Saving compact workspace instead.', saveErr.message);
    save(saveFile, 'Tresults','T_rank_primary','T_rank_consensus','Tweights','TfitQuality', ...
        'Tcorr_DSI','Tcorr_theta','DSI_stability_BC','DSI_consensus_BC','theta_latent_BC','-v7');
end
fprintf('MAT workspace saved: %s\n', saveFile);

summaryFile = fullfile(logDir,'gc32_dsi_analysis_summary.txt');
fid=fopen(summaryFile,'w');
fprintf(fid,'GC32 DSI analysis summary\n');
fprintf(fid,'Generated: %s\n',datestr(now));
fprintf(fid,'runMode=%s, Nmc=%d, Nnull=%d, Nweight=%d\n\n',runMode,Nmc,Nnull,Nweight);
fprintf(fid,'Implemented safe analysis upgrades:\n');
fprintf(fid,'  learned weights: fixed/PCA/entropy/stability\n');
fprintf(fid,'  latent one-factor dynamical-complexity score\n');
fprintf(fid,'  empirical incoherent-component null significance model\n');
fprintf(fid,'  physical relaxed-proxy, rotation-free and radial profile-shuffle nulls\n');
fprintf(fid,'  hierarchical Bayesian latent factor model with posterior theta/lambda uncertainty\n');
fprintf(fid,'  mixture-model classification\n');
fprintf(fid,'  covariance-informed MC with shared systematic/radial perturbations\n');
fprintf(fid,'  physical profile-shape fitting and lambda_R spin proxy\n');
fprintf(fid,'  EM-PCA imputation and FDR-corrected external validation\n\n');

fprintf(fid,'Top 12 primary DSI_stability_BC clusters:\n');
for i=1:min(12,height(T_rank_primary))
    fprintf(fid,'%2d %-8s DSI=% .4f err=%.4f theta=% .4f pNull=%.4g mix=%s flag=%s dom=%s\n', ...
        i,T_rank_primary.cluster_id(i),T_rank_primary.DSI_stability_BC(i),T_rank_primary.DSI_stability_err(i), ...
        T_rank_primary.theta_latent_BC(i),T_rank_primary.p_null_DSI(i),T_rank_primary.mixture_class(i), ...
        T_rank_primary.interpretation_flag(i),T_rank_primary.dominant_component(i));
end

fprintf(fid,'\nRecommended manuscript primary table: gc32_dsi_primary_ranked_dsi_stability_bc.csv\n');
fprintf(fid,'Recommended robustness table: gc32_dsi_quality_aware_ranked_dsi_stability_qa.csv\n');
fprintf(fid,'Recommended latent table: gc32_dsi_classical_latent_ranked_theta_bc.csv\n');
fprintf(fid,'Recommended consensus table: gc32_dsi_consensus_ranked_bc.csv\n');
fprintf(fid,'Caution: DSI is a complexity indicator, not direct IMBH evidence.\n');
fclose(fid);

fprintf('\n============================================================\n');
fprintf('GC32 DSI analysis pipeline completed\n');
fprintf('Primary results: %s\n', fullfile(tabDir,'gc32_dsi_primary_ranked_dsi_stability_bc.csv'));
fprintf('Summary: %s\n', summaryFile);
fprintf('============================================================\n');
diary off;

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function d=autoDetectDataDir(scriptDir)
    packageRoot = fileparts(fileparts(scriptDir));
    candidates = { ...
        fullfile(packageRoot,'data','input'), ...
        fullfile(fileparts(scriptDir),'data','input'), ...
        fullfile(scriptDir,'data','input'), ...
        fullfile(pwd,'data','input'), ...
        pwd, scriptDir};
    target = '12_gc32_main_dsi_input_matrix_matlab.csv';
    d = '';
    for i = 1:numel(candidates)
        if isfolder(candidates{i}) && isfile(fullfile(candidates{i},target))
            d = candidates{i};
            return;
        end
    end
    error('Could not auto-detect the input CSV directory; set dataDir manually.');
end

function mkdirIfNeeded(d), if ~isfolder(d), mkdir(d); end, end

function checkFiles(F)
    n=fieldnames(F);
    for i=1:numel(n)
        if ~isfile(F.(n{i})), error('Missing file: %s',F.(n{i})); end
    end
end

function [Tfeat,X]=computeAllFeatures(clusterID,HST,GAIA,VEL,minSlope,minCurv)
    n=numel(clusterID);
    Tfeat=table(); Tfeat.cluster_id=clusterID;
    names=["hst_n_bins","gaia_n_bins","vel_n_bins","gradient_strength","anisotropy_abs","rotation_support_robust", ...
        "profile_curvature","central_enhance","outer_disturbance","lambda_R_gaia", ...
        "hst_eta","hst_sigma0","hst_rc","gaia_eta","gaia_sigma0","gaia_rc","vel_eta","vel_sigma0","vel_rc"];
    for nm=names, Tfeat.(nm)=nan(n,1); end

    for i=1:n
        cid=clusterID(i);
        H=HST(HST.cluster_id==cid,:);
        G=GAIA(GAIA.cluster_id==cid & GAIA.radius_arcsec>0,:);
        V=VEL(VEL.cluster_id==cid,:);
        hf=hstFeatures(H,minSlope,minCurv);
        gf=gaiaFeatures(G,minSlope,minCurv);
        vf=velFeatures(V,minSlope,minCurv);

        Tfeat.hst_n_bins(i)=height(H); Tfeat.gaia_n_bins(i)=height(G); Tfeat.vel_n_bins(i)=height(V);
        Tfeat.gradient_strength(i)=-meanFinite([hf.slope,gf.slope,vf.slope]);
        Tfeat.anisotropy_abs(i)=hf.anisotropy_abs;
        Tfeat.rotation_support_robust(i)=gf.rotation_support_robust;
        Tfeat.profile_curvature(i)=meanFinite([hf.curvature,gf.curvature,vf.curvature]);
        Tfeat.central_enhance(i)=safeDiv(vf.center_sigma,vf.median_sigma);
        Tfeat.outer_disturbance(i)=meanFinite([hf.outer_inner,gf.outer_inner,vf.outer_inner]);
        Tfeat.lambda_R_gaia(i)=gf.lambda_R;

        Tfeat.hst_eta(i)=hf.eta; Tfeat.hst_sigma0(i)=hf.sigma0; Tfeat.hst_rc(i)=hf.rc;
        Tfeat.gaia_eta(i)=gf.eta; Tfeat.gaia_sigma0(i)=gf.sigma0; Tfeat.gaia_rc(i)=gf.rc;
        Tfeat.vel_eta(i)=vf.eta; Tfeat.vel_sigma0(i)=vf.sigma0; Tfeat.vel_rc(i)=vf.rc;
    end
    X=[Tfeat.gradient_strength,Tfeat.anisotropy_abs,Tfeat.rotation_support_robust, ...
       Tfeat.profile_curvature,Tfeat.central_enhance,Tfeat.outer_disturbance];
end

function f=hstFeatures(T,minSlope,minCurv)
    f=initF(); if isempty(T), return; end
    r=T.radius_arcsec; s=T.pm_sigma_total_masyr; e=T.pm_sigma_total_err_masyr;
    sr=T.pm_sigma_radial_masyr; st=T.pm_sigma_tangential_masyr;
    f.slope=logSlope(r,s,e,minSlope); f.curvature=logCurvature(r,s,e,minCurv);
    f.outer_inner=outerInner(r,s); f.center_sigma=firstByRadius(r,s); f.median_sigma=medianFinitePos(s);
    beta=1-(st.^2./sr.^2); f.anisotropy_abs=median(betaAbs(beta),'omitnan');
    [f.sigma0,f.rc,f.eta]=fitPhysicalProfile(r,s,e);
end

function f=gaiaFeatures(T,minSlope,minCurv)
    f=initF(); if isempty(T), return; end
    r=T.radius_arcsec; s=T.pm_disp_p50_masyr; rot=T.pm_rot_p50_masyr;
    e=0.5.*(T.pm_disp_p84_masyr-T.pm_disp_p16_masyr);
    f.slope=logSlope(r,s,e,minSlope); f.curvature=logCurvature(r,s,e,minCurv);
    f.outer_inner=outerInner(r,s); f.center_sigma=firstByRadius(r,s); f.median_sigma=medianFinitePos(s);
    f.rotation_support_robust=safeDiv(percentile1d(abs(rot(isfinite(rot))),90),f.median_sigma);
    f.lambda_R=lambdaR(r,rot,s);
    [f.sigma0,f.rc,f.eta]=fitPhysicalProfile(r,s,e);
end

function f=velFeatures(T,minSlope,minCurv)
    f=initF(); if isempty(T), return; end
    r=T.radius_arcsec; s=T.sigma_kms; e=mean([T.sigma_err_up_kms,T.sigma_err_low_kms],2,'omitnan');
    f.slope=logSlope(r,s,e,minSlope); f.curvature=logCurvature(r,s,e,minCurv);
    f.outer_inner=outerInner(r,s); f.center_sigma=firstByRadius(r,s); f.median_sigma=medianFinitePos(s);
    [f.sigma0,f.rc,f.eta]=fitPhysicalProfile(r,s,e);
end

function f=initF()
    f.slope=NaN; f.curvature=NaN; f.outer_inner=NaN; f.center_sigma=NaN; f.median_sigma=NaN;
    f.anisotropy_abs=NaN; f.rotation_support_robust=NaN; f.lambda_R=NaN;
    f.sigma0=NaN; f.rc=NaN; f.eta=NaN;
end

function b=betaAbs(beta), b=abs(beta(isfinite(beta))); end
function m=medianFinitePos(x), x=x(isfinite(x)&x>0); if isempty(x), m=NaN; else, m=median(x,'omitnan'); end, end

function val=lambdaR(r,v,s)
    r=double(r(:)); v=double(v(:)); s=double(s(:));
    ok=isfinite(r)&isfinite(v)&isfinite(s)&(r>0)&(s>0);
    r=r(ok); v=v(ok); s=s(ok);
    if isempty(r), val=NaN; return; end
    val=sum(r.*abs(v),'omitnan')/sum(r.*sqrt(v.^2+s.^2),'omitnan');
end

function [sigma0,rc,eta]=fitPhysicalProfile(r,s,e)
    r=double(r(:)); s=double(s(:)); e=double(e(:));
    ok=isfinite(r)&isfinite(s)&r>0&s>0;
    if numel(e)==numel(s), ok=ok&isfinite(e)&e>0; end
    r=r(ok); s=s(ok);
    if numel(e)==numel(ok), e=e(ok); else, e=ones(size(s)); end
    if numel(s)<4 || numel(unique(round(log10(r),12)))<3
        sigma0=NaN; rc=NaN; eta=NaN; return;
    end
    rmed=median(r,'omitnan'); smax=max(s); eta0=max(0.1,-logSlope(r,s,e,3));
    p0=log([smax,max(rmed,eps),max(eta0,0.05)]);
    obj=@(p) sum(((s - exp(p(1)).*(1+(r./exp(p(2))).^2).^(-exp(p(3))/2))./max(e,median(e,'omitnan')*0.2)).^2,'omitnan');
    old=warning('off','all');
    try
        p=fminsearch(obj,p0,optimset('Display','off','MaxIter',1000,'MaxFunEvals',3000));
        sigma0=exp(p(1)); rc=exp(p(2)); eta=exp(p(3));
        if eta>10 || rc<=0 || sigma0<=0, sigma0=NaN; rc=NaN; eta=NaN; end
    catch
        sigma0=NaN; rc=NaN; eta=NaN;
    end
    warning(old);
end

function s=logSlope(r,y,e,minN)
    [x,yy,w]=prepLog(r,y,e);
    if numel(x)<minN || numel(unique(round(x,12)))<2, s=NaN; return; end
    p=weightedPolyfit(x,yy,1,w); s=p(1); if ~isfinite(s), s=NaN; end
end

function c=logCurvature(r,y,e,minN)
    [x,yy,w]=prepLog(r,y,e);
    if numel(x)<minN || numel(unique(round(x,12)))<3, c=NaN; return; end
    p=weightedPolyfit(x,yy,2,w); c=abs(2*p(1)); if ~isfinite(c), c=NaN; end
end

function [x,yy,w]=prepLog(r,y,e)
    r=double(r(:)); y=double(y(:)); e=double(e(:));
    ok=isfinite(r)&isfinite(y)&r>0&y>0;
    if numel(e)==numel(y), ok=ok&isfinite(e)&e>0; end
    r=r(ok); y=y(ok);
    x=log10(r); yy=log10(y);
    if numel(e)==numel(ok)
        e=e(ok); w=1./max(e,median(e,'omitnan')*0.05);
    else
        w=ones(size(y));
    end
end

function p=weightedPolyfit(x,y,deg,w)
    x=x(:); y=y(:); w=w(:);
    ok=isfinite(x)&isfinite(y)&isfinite(w)&w>0;
    x=x(ok); y=y(ok); w=w(ok);
    if numel(x)<deg+1 || numel(unique(round(x,12)))<deg+1
        p=nan(1,deg+1); return;
    end
    V=zeros(numel(x),deg+1);
    for j=0:deg, V(:,deg+1-j)=x.^j; end
    A=V.*sqrt(w); b=y.*sqrt(w);
    if rank(A)<deg+1, p=nan(1,deg+1); return; end
    old=warning('off','all'); c=A\b; warning(old); p=c(:).';
end

function val=outerInner(r,y)
    r=double(r(:)); y=double(y(:));
    ok=isfinite(r)&isfinite(y)&y>0; r=r(ok); y=y(ok);
    if numel(y)<3, val=NaN; return; end
    [~,ord]=sort(r); y=y(ord); k=max(1,ceil(0.2*numel(y)));
    val=safeDiv(median(y(end-k+1:end),'omitnan'),median(y(1:k),'omitnan'));
end

function val=firstByRadius(r,y)
    r=double(r(:)); y=double(y(:)); ok=isfinite(r)&isfinite(y)&r>0&y>0;
    r=r(ok); y=y(ok); if isempty(y), val=NaN; return; end
    [~,i]=min(r); val=y(i);
end

function v=meanFinite(x), x=x(isfinite(x)); if isempty(x), v=NaN; else, v=mean(x); end, end
function z=safeDiv(a,b), if ~isfinite(a)||~isfinite(b)||b==0, z=NaN; else, z=a/b; end, end

function chem=computeChemComplexity(clusterID,T)
    chem=nan(numel(clusterID),1);
    if isempty(T)||~ismember("cluster_id",string(T.Properties.VariableNames)), return; end
    T.cluster_id=string(T.cluster_id); vars=string(T.Properties.VariableNames);
    comps=["FE_H_sigma","N_FE_sigma","AL_FE_sigma","MG_FE_sigma"];
    for i=1:numel(clusterID)
        idx=find(T.cluster_id==clusterID(i),1); if isempty(idx), continue; end
        ss=0; n=0;
        for c=comps
            if ismember(c,vars)
                val=T.(c)(idx); if isfinite(val), ss=ss+val^2; n=n+1; end
            end
        end
        if n>0, chem(i)=sqrt(ss); end
    end
end

function Q=computeQualityScore(T)
    qH=min(1,T.hst_n_bins./8); qG=min(1,T.gaia_n_bins./80); qV=min(1,T.vel_n_bins./12);
    comp=[T.gradient_strength,T.anisotropy_abs,T.rotation_support_robust,T.profile_curvature,T.central_enhance,T.outer_disturbance];
    qC=1-0.08*sum(~isfinite(comp),2);
    Q=(qH.*qG.*qV.*max(0.4,qC)).^(1/4); Q(~isfinite(Q))=0; Q=max(0,min(1,Q));
end

function [Ximp,info]=emPCAImpute(X,k,maxIter,tol)
    Ximp=X; miss=~isfinite(X);
    info.method="EM_PCA"; info.converged=false; info.iterations=0;
    [Ximp,~]=medianImpute(Ximp);
    prev=Ximp;
    for it=1:maxIter
        [Z,loc,scale]=robustZ(Ximp,1); %#ok<ASGLU>
        try
            [coeff,score]=pcaFallback(Z);
            Zhat=score(:,1:k)*coeff(:,1:k)';
            Xhat=Zhat.*scale+loc;
            Ximp(miss)=Xhat(miss);
        catch
            info.method="EM_PCA_failed";
            return;
        end
        delta=norm(Ximp(:)-prev(:))/max(1,norm(prev(:)));
        prev=Ximp;
        if delta<tol
            info.converged=true; info.iterations=it; return;
        end
    end
    info.iterations=maxIter;
end

function [Ximp,info]=medianImpute(X)
    Ximp=X; miss=~isfinite(X); info.method="median";
    for j=1:size(X,2)
        m=median(X(:,j),'omitnan'); if ~isfinite(m), m=0; end
        Ximp(miss(:,j),j)=m;
    end
end

function [Z,loc,scale]=robustZ(X,c)
    loc=median(X,1,'omitnan'); scale=nan(1,size(X,2));
    for j=1:size(X,2)
        scale(j)=c*median(abs(X(:,j)-loc(j)),'omitnan');
        if ~isfinite(scale(j))||scale(j)<=0, scale(j)=std(X(:,j),0,'omitnan'); end
        if ~isfinite(scale(j))||scale(j)<=0, scale(j)=1; end
    end
    Z=(X-loc)./scale;
end

function Z=winsorizeMatrix(Z,cap), Z=min(max(Z,-cap),cap); end

function [w,load1,expl]=pcaWeights(Z)
    [coeff,~,expl]=pcaFallback(Z);
    load1=coeff(:,1);
    w=abs(load1(:)); if sum(w)==0, w=ones(size(w)); end
    w=w./sum(w);
end

function [coeff,score,explained]=pcaFallback(Z)
    Z=Z-mean(Z,1,'omitnan');
    try
        [coeff,score,~,~,explained]=pca(Z);
    catch
        [U,S,V]=svd(Z,'econ'); score=U*S; coeff=V;
        latent=diag(S).^2/(size(Z,1)-1); explained=100*latent/sum(latent);
    end
end

function w=entropyWeights(Z)
    X=abs(Z); X=X-min(X,[],1); X=X+eps;
    P=X./sum(X,1);
    n=size(Z,1);
    E=-sum(P.*log(P),1)./log(n);
    d=1-E; if sum(d)<=0, d=ones(size(d)); end
    w=d(:)./sum(d);
end

function ok=canUseParallel(useParallel)
    ok=false;
    if ~useParallel, return; end
    try
        v=ver; hasPCT=any(strcmp({v.Name},'Parallel Computing Toolbox'));
        if hasPCT
            p=gcp('nocreate'); if isempty(p), parpool('threads'); end
            ok=true;
        end
    catch
        ok=false;
    end
end

function [MC_Z,MC_X]=monteCarloComponents(clusterID,HST,GAIA,VEL,loc,scale,zCap,Nmc,minSlope,minCurv,usePar)
    n=numel(clusterID); K=numel(loc);
    MC_Z=nan(n,K,Nmc); MC_X=nan(n,K,Nmc);
    if usePar
        parfor b=1:Nmc
            [Zb,Xb]=mcOne(clusterID,HST,GAIA,VEL,loc,scale,zCap,minSlope,minCurv);
            MC_Z(:,:,b)=Zb; MC_X(:,:,b)=Xb;
        end
    else
        for b=1:Nmc
            if mod(b,max(1,round(Nmc/10)))==0, fprintf('  MC %d/%d\n',b,Nmc); end
            [Zb,Xb]=mcOne(clusterID,HST,GAIA,VEL,loc,scale,zCap,minSlope,minCurv);
            MC_Z(:,:,b)=Zb; MC_X(:,:,b)=Xb;
        end
    end
end

function [Zb,Xb]=mcOne(clusterID,HST,GAIA,VEL,loc,scale,zCap,minSlope,minCurv)
    n=numel(clusterID); K=numel(loc);
    Xb=nan(n,K); Zb=nan(n,K);
    % Shared systematic terms approximate covariance within each data source.
    sysH=0.03*randn(n,1); sysG=0.03*randn(n,1); sysV=0.03*randn(n,1);
    for i=1:n
        cid=clusterID(i);
        H=perturbHST(HST(HST.cluster_id==cid,:),sysH(i));
        G=perturbGaia(GAIA(GAIA.cluster_id==cid & GAIA.radius_arcsec>0,:),sysG(i));
        V=perturbVel(VEL(VEL.cluster_id==cid,:),sysV(i));
        hf=hstFeatures(H,minSlope,minCurv); gf=gaiaFeatures(G,minSlope,minCurv); vf=velFeatures(V,minSlope,minCurv);
        x=[-meanFinite([hf.slope,gf.slope,vf.slope]),hf.anisotropy_abs,gf.rotation_support_robust, ...
            meanFinite([hf.curvature,gf.curvature,vf.curvature]),safeDiv(vf.center_sigma,vf.median_sigma), ...
            meanFinite([hf.outer_inner,gf.outer_inner,vf.outer_inner])];
        for k=1:K, if ~isfinite(x(k)), x(k)=loc(k); end, end
        z=(x-loc)./scale; z=min(max(z,-zCap),zCap);
        Xb(i,:)=x; Zb(i,:)=z;
    end
end

function H=perturbHST(H,sys)
    if isempty(H), return; end
    r = H.radius_arcsec;
    cols=["pm_sigma_total_masyr","pm_sigma_radial_masyr","pm_sigma_tangential_masyr"];
    errs=["pm_sigma_total_err_masyr","pm_sigma_radial_err_masyr","pm_sigma_tangential_err_masyr"];
    for j=1:numel(cols)
        if ismember(cols(j),string(H.Properties.VariableNames)) && ismember(errs(j),string(H.Properties.VariableNames))
            e=H.(errs(j));
            noise = drawRadialCorrelatedNoise(r,e,0.55);
            H.(cols(j))=H.(cols(j)).*(1+sys)+noise;
            H.(cols(j))(H.(cols(j))<=0)=NaN;
        end
    end
end

function G=perturbGaia(G,sys)
    if isempty(G), return; end
    r = G.radius_arcsec;
    e=0.5*(G.pm_disp_p84_masyr-G.pm_disp_p16_masyr);
    G.pm_disp_p50_masyr=G.pm_disp_p50_masyr.*(1+sys)+drawRadialCorrelatedNoise(r,e,0.75);
    G.pm_disp_p50_masyr(G.pm_disp_p50_masyr<=0)=NaN;
    er=0.5*(G.pm_rot_p84_masyr-G.pm_rot_p16_masyr);
    G.pm_rot_p50_masyr=G.pm_rot_p50_masyr.*(1+sys)+drawRadialCorrelatedNoise(r,er,0.75);
end

function V=perturbVel(V,sys)
    if isempty(V), return; end
    r = V.radius_arcsec;
    e=mean([V.sigma_err_up_kms,V.sigma_err_low_kms],2,'omitnan');
    V.sigma_kms=V.sigma_kms.*(1+sys)+drawRadialCorrelatedNoise(r,e,0.65);
    V.sigma_kms(V.sigma_kms<=0)=NaN;
end

function noise = drawRadialCorrelatedNoise(r,e,ell)
    % Gaussian-kernel radial covariance approximation in log-radius space.
    % This is a covariance-informed perturbation, not a claim of full observed covariance.
    r=double(r(:)); e=double(e(:));
    n=numel(r);
    noise=zeros(n,1);
    ok=isfinite(r)&r>0&isfinite(e)&e>0;
    if sum(ok)==0, return; end
    rr=log10(r(ok)); ee=e(ok);
    if sum(ok)==1
        noise(ok)=ee.*randn;
        return;
    end
    D=rr-rr';
    C=(ee*ee').*exp(-0.5*(D./max(ell,eps)).^2);
    C=C+diag((0.03*ee).^2+1e-10);
    [L,p]=chol(C,'lower');
    if p~=0
        C=diag(max(ee,eps).^2);
        L=chol(C,'lower');
    end
    noise(ok)=L*randn(sum(ok),1);
end

function [D, Dmc, Dbc, Derr, Dbias]=biasCorrectedIndex(Z,MC_Z,w)
    D=Z*w(:);
    n=size(Z,1); B=size(MC_Z,3);
    Dmc=nan(n,B);
    for b=1:B, Dmc(:,b)=MC_Z(:,:,b)*w(:); end
    med=median(Dmc,2,'omitnan'); p16=prctileLocal(Dmc,16,2); p84=prctileLocal(Dmc,84,2);
    Dbias=med-D; Derr=(p84-p16)/2; Dbc=D-Dbias;
end

function [w,info]=learnStabilityWeights(Z,MC_Z,N,wMax,balancePenaltyStrength)
    % v4.2.1 PATCH:
    % Accepts wMax and balancePenaltyStrength. This prevents one-component
    % collapse and fixes the "Too many input arguments" error.
    K=size(Z,2); B=size(MC_Z,3);
    bestScore=Inf;
    equalW=ones(1,K)/K;
    bestW=projectCappedSimplex(equalW,wMax);

    % Deterministic candidates.
    candidateW=zeros(4,K);
    candidateW(1,:)=equalW;
    [wpca,~,~]=pcaWeights(Z);
    candidateW(2,:)=projectCappedSimplex(wpca(:)',wMax);
    candidateW(3,:)=projectCappedSimplex(entropyWeights(Z)',wMax);
    invScale=1./max(std(Z,0,1,'omitnan'),eps);
    candidateW(4,:)=projectCappedSimplex(invScale./sum(invScale),wMax);

    maxB=min(B,300);
    for t=1:N
        if t<=size(candidateW,1)
            wt=candidateW(t,:);
        else
            wt=projectCappedSimplex(rand(1,K),wMax);
        end

        D=Z*wt(:);
        r0=tiedRanking(-D);
        rankSig=0;
        for b=1:maxB
            Db=MC_Z(:,:,b)*wt(:);
            rb=tiedRanking(-Db);
            rankSig=rankSig+mean(abs(rb-r0),'omitnan');
        end
        rankSig=rankSig/maxB;

        balancePenalty=balancePenaltyStrength*sum((wt-equalW).^2);
        effectiveN=1/sum(wt.^2);
        effectivePenalty=0.10*(K/effectiveN-1); % discourages concentration

        score=rankSig+balancePenalty+effectivePenalty;
        if score<bestScore
            bestScore=score;
            bestW=wt;
        end
    end

    w=projectCappedSimplex(bestW,wMax)';
    info.bestScore=bestScore;
    info.wMax=wMax;
    info.maxWeight=max(w);
    info.effectiveN=1/sum(w.^2);
end

function w=projectCappedSimplex(w,wMax)
    % Simple cap-and-redistribute projection onto nonnegative simplex.
    w=max(w(:)',0);
    K=numel(w);
    if sum(w)==0
        w=ones(1,K)/K;
        return;
    end
    w=w/sum(w);
    if wMax<1/K
        wMax=1/K;
    end
    for it=1:200
        over=w>wMax+1e-12;
        if ~any(over)
            break;
        end
        excess=sum(w(over)-wMax);
        w(over)=wMax;
        under=~over;
        if any(under)
            denom=sum(w(under));
            if denom>0
                w(under)=w(under)+excess*w(under)/denom;
            else
                w=ones(1,K)/K;
            end
        end
        w=max(w,0);
        w=w/sum(w);
    end
    w=min(w,wMax);
    w=w/sum(w);
end


function [thetaMed,thetaErr,lambdaMed,lambdaErr,alphaMed,sigmaIntMed,info]=bayesianLatentDSI(Z,MC_Z,nIter,burn,thin)
    % Lightweight hierarchical Bayesian one-factor latent model.
    % z_ij ~ Normal(alpha_j + lambda_j*theta_i, sigma_int_j^2 + u_ij)
    % u_ij is estimated from MC propagated component scatter.
    Z=double(Z); [n,K]=size(Z);
    U=var(MC_Z,0,3,'omitnan');
    U(~isfinite(U))=0;
    uFloor=0.03^2;
    U=max(U,uFloor);

    [theta0,lambda0,~]=latentOneFactor(Z);
    theta=standardizeVector(theta0);
    lambda=lambda0(:)';
    alpha=zeros(1,K);
    sig2=max(0.05^2, 0.25*var(Z,0,1,'omitnan'));

    keepIdx = burn+1:thin:nIter;
    nKeep=numel(keepIdx);
    TH=nan(n,nKeep); LA=nan(K,nKeep); AL=nan(K,nKeep); SG=nan(K,nKeep);
    kk=0;

    for it=1:nIter
        V=U + repmat(sig2,n,1);
        % theta update
        for i=1:n
            prec=1 + sum((lambda.^2)./V(i,:), 'omitnan');
            mu=sum(lambda.*(Z(i,:)-alpha)./V(i,:), 'omitnan')/max(prec,eps);
            theta(i)=mu + sqrt(1/max(prec,eps))*randn;
        end
        % centre/scale theta for identifiability
        muT=mean(theta,'omitnan'); sdT=std(theta,0,'omitnan'); if ~isfinite(sdT)||sdT<=0, sdT=1; end
        alpha=alpha + lambda*muT;
        lambda=lambda*sdT;
        theta=(theta-muT)/sdT;

        % alpha and lambda updates by component
        V=U + repmat(sig2,n,1);
        for j=1:K
            vj=V(:,j);
            % alpha_j | rest
            precA=1/25 + sum(1./vj,'omitnan');
            muA=sum((Z(:,j)-lambda(j)*theta)./vj,'omitnan')/max(precA,eps);
            alpha(j)=muA + sqrt(1/max(precA,eps))*randn;
            % lambda_j | rest
            precL=1/4 + sum((theta.^2)./vj,'omitnan');
            muL=sum(theta.*(Z(:,j)-alpha(j))./vj,'omitnan')/max(precL,eps);
            lambda(j)=muL + sqrt(1/max(precL,eps))*randn;
            % intrinsic scatter approximate inverse-gamma update after removing known measurement variance
            resid=Z(:,j)-alpha(j)-lambda(j)*theta;
            eff=max(resid.^2-U(:,j),0);
            a=2 + 0.5*n;
            b=0.05 + 0.5*sum(eff,'omitnan');
            g=randg(a)/max(b,eps); % Gamma(shape=a, rate=b)
            sig2(j)=max(1/max(g,eps),0.005^2);
        end

        % orient high theta with high total complexity
        if corrLocal(theta,sum(Z,2))<0
            theta=-theta; lambda=-lambda;
        end

        if it>burn && mod(it-burn,thin)==0
            kk=kk+1; TH(:,kk)=theta(:); LA(:,kk)=lambda(:); AL(:,kk)=alpha(:); SG(:,kk)=sqrt(sig2(:));
        end
    end

    thetaMed=median(TH,2,'omitnan');
    thetaErr=(prctileLocal(TH,84,2)-prctileLocal(TH,16,2))./2;
    lambdaMed=median(LA,2,'omitnan');
    lambdaErr=(prctileLocal(LA,84,2)-prctileLocal(LA,16,2))./2;
    alphaMed=median(AL,2,'omitnan');
    sigmaIntMed=median(SG,2,'omitnan');
    % final orientation
    if corrLocal(thetaMed,sum(Z,2))<0
        thetaMed=-thetaMed; lambdaMed=-lambdaMed;
    end
    info.nIter=nIter; info.burn=burn; info.thin=thin; info.nPosteriorSamples=nKeep;
    info.meanThetaErr=mean(thetaErr,'omitnan');
    info.meanAbsLambda=sum(abs(lambdaMed),'omitnan');
end

function [scores,p,info]=physicalRelaxedNullProxy(Z,w,obs,N)
    % Component-space relaxed proxy. Components are centred on the empirically
    % regular part of the sample and shrunk toward low rotation/anisotropy/curvature.
    [n,K]=size(Z); %#ok<ASGLU>
    D=Z*w(:);
    baseIdx = D <= percentile1d(D,40);
    if sum(baseIdx)<5, baseIdx = D <= median(D,'omitnan'); end
    mu=median(Z(baseIdx,:),1,'omitnan');
    % force physically regular tendencies for the most diagnostic components
    % columns: G, A, R, C, E, T in robust-z space.
    if K>=2, mu(2)=min(mu(2), percentile1d(Z(:,2),35)); end      % low anisotropy
    if K>=3, mu(3)=min(mu(3), percentile1d(Z(:,3),25)); end      % low rotation
    if K>=4, mu(4)=min(mu(4), percentile1d(Z(:,4),35)); end      % low curvature
    if K>=6, mu(6)=min(mu(6), percentile1d(abs(Z(:,6)),45)); end % regular outer structure
    C=cov(Z(baseIdx,:),'omitrows');
    if any(~isfinite(C(:))) || rank(C)<K
        C=cov(Z,'omitrows');
    end
    if any(~isfinite(C(:))), C=eye(K)*0.25; end
    C=0.45*C + 0.04*eye(K); % shrink: conservative relaxed proxy
    L=safeChol(C);
    scores=nan(N,1);
    for b=1:N
        zsim=mu + randn(1,K)*L';
        zsim=min(max(zsim,-3.5),3.5);
        scores(b)=zsim*w(:);
    end
    p=nan(numel(obs),1);
    for i=1:numel(obs)
        p(i)=(1+sum(scores>=obs(i)))/(N+1);
    end
    info.mu=mu; info.cov=C; info.N=N; info.description="component_space_physical_relaxed_proxy";
end

function [scores,p,drop]=rotationFreeNull(Z,w,obs,N)
    [n,K]=size(Z); scores=nan(n,N);
    lowRot = Z(:,3) <= percentile1d(Z(:,3),25);
    rotPool = Z(lowRot,3); if isempty(rotPool), rotPool=median(Z(:,3),'omitnan'); end
    for b=1:N
        Zp=nan(size(Z));
        for k=1:K
            if k==3
                Zp(:,k)=rotPool(randi(numel(rotPool),n,1));
            else
                Zp(:,k)=Z(randperm(n),k);
            end
        end
        scores(:,b)=Zp*w(:);
    end
    p=nan(n,1);
    for i=1:n
        p(i)=(1+sum(scores(i,:)>=obs(i)))/(N+1);
    end
    noRotScore=median(scores,2,'omitnan');
    drop=obs-noRotScore;
end

function [scores,p]=profileShuffleNull(clusterID,HST,GAIA,VEL,loc,scale,zCap,w,obs,N,minSlope,minCurv,usePar)
    n=numel(clusterID); scores=nan(n,N);
    if usePar
        parfor b=1:N
            scores(:,b)=profileShuffleOne(clusterID,HST,GAIA,VEL,loc,scale,zCap,w,minSlope,minCurv);
        end
    else
        for b=1:N
            if mod(b,max(1,round(N/5)))==0, fprintf('  profile-shuffle null %d/%d\n',b,N); end
            scores(:,b)=profileShuffleOne(clusterID,HST,GAIA,VEL,loc,scale,zCap,w,minSlope,minCurv);
        end
    end
    p=nan(n,1);
    for i=1:n
        p(i)=(1+sum(scores(i,:)>=obs(i)))/(N+1);
    end
end

function score=profileShuffleOne(clusterID,HST,GAIA,VEL,loc,scale,zCap,w,minSlope,minCurv)
    n=numel(clusterID); K=numel(loc); score=nan(n,1);
    for i=1:n
        cid=clusterID(i);
        H=shuffleProfileRadius(HST(HST.cluster_id==cid,:));
        G=shuffleProfileRadius(GAIA(GAIA.cluster_id==cid & GAIA.radius_arcsec>0,:));
        V=shuffleProfileRadius(VEL(VEL.cluster_id==cid,:));
        hf=hstFeatures(H,minSlope,minCurv); gf=gaiaFeatures(G,minSlope,minCurv); vf=velFeatures(V,minSlope,minCurv);
        x=[-meanFinite([hf.slope,gf.slope,vf.slope]),hf.anisotropy_abs,gf.rotation_support_robust, ...
            meanFinite([hf.curvature,gf.curvature,vf.curvature]),safeDiv(vf.center_sigma,vf.median_sigma), ...
            meanFinite([hf.outer_inner,gf.outer_inner,vf.outer_inner])];
        for k=1:K, if ~isfinite(x(k)), x(k)=loc(k); end, end
        z=(x-loc)./scale; z=min(max(z,-zCap),zCap);
        score(i)=z*w(:);
    end
end

function T=shuffleProfileRadius(T)
    if isempty(T) || ~ismember("radius_arcsec",string(T.Properties.VariableNames)), return; end
    r=T.radius_arcsec;
    if numel(r)>2, T.radius_arcsec=r(randperm(numel(r))); end
end

function L=safeChol(C)
    C=(C+C')/2;
    [L,p]=chol(C,'lower');
    if p~=0
        [V,D]=eig(C); d=max(diag(D),1e-5); L=V*diag(sqrt(d));
    end
end

function [theta,lambda,info]=latentOneFactor(Z)
    [coeff,score,explained]=pcaFallback(Z);
    lambda=coeff(:,1);
    % Orient so high theta correlates positively with row-sum complexity.
    if corrLocal(score(:,1),sum(Z,2))<0, lambda=-lambda; score(:,1)=-score(:,1); end
    lambda=lambda./sum(abs(lambda));
    theta=Z*lambda;
    info.explainedPC1=explained(1);
end

function [nullScores,p]=relaxedNullPermutation(Z,w,obs,N)
    n=size(Z,1); K=size(Z,2);
    nullScores=nan(n,N);
    for b=1:N
        Zp=nan(size(Z));
        for k=1:K
            Zp(:,k)=Z(randperm(n),k);
        end
        nullScores(:,b)=Zp*w(:);
    end
    p=nan(n,1);
    for i=1:n
        p(i)=(1+sum(nullScores(i,:)>=obs(i)))/(N+1);
    end
end

function sig=classifyNullSignificance(p,q)
    sig=strings(numel(p),1);
    for i=1:numel(p)
        if q(i)<0.05, sig(i)="FDR_significant";
        elseif p(i)<0.05, sig(i)="nominal_significant";
        elseif p(i)<0.10, sig(i)="marginal";
        else, sig(i)="not_significant";
        end
    end
end

function G=fitGMM_EM(Y,K,maxIter,tol)
    Y=double(Y); [n,d]=size(Y);
    [~,ord]=sort(Y(:,1)); cuts=round(linspace(1,n,K+1));
    mu=nan(K,d);
    for k=1:K, mu(k,:)=mean(Y(ord(cuts(k):cuts(k+1)),:),1,'omitnan'); end
    Sigma=repmat(eye(d),1,1,K); piK=ones(1,K)/K; llOld=-Inf;
    gamma=ones(n,K)/K;
    for it=1:maxIter
        dens=zeros(n,K);
        for k=1:K, dens(:,k)=piK(k)*mvnpdfLocal(Y,mu(k,:),Sigma(:,:,k)); end
        denom=sum(dens,2)+eps; gamma=dens./denom;
        Nk=sum(gamma,1);
        for k=1:K
            piK(k)=Nk(k)/n;
            mu(k,:)=sum(Y.*gamma(:,k),1)/max(Nk(k),eps);
            Xc=Y-mu(k,:);
            Sigma(:,:,k)=(Xc'*(Xc.*gamma(:,k)))/max(Nk(k),eps)+1e-4*eye(d);
        end
        ll=sum(log(denom),'omitnan');
        if abs(ll-llOld)<tol, break; end
        llOld=ll;
    end
    % Sort states by mean DSI dimension.
    [~,idx]=sort(mu(:,1),'ascend');
    G.mu=mu(idx,:); G.Sigma=Sigma(:,:,idx); G.pi=piK(idx); G.gamma=gamma(:,idx); G.logLik=llOld; G.K=K;
end

function p=mvnpdfLocal(X,mu,S)
    d=size(X,2); Xc=X-mu;
    [R,flag]=chol(S);
    if flag~=0, S=S+1e-3*eye(d); R=chol(S); end
    q=sum((Xc/R).^2,2);
    p=exp(-0.5*q)/((2*pi)^(d/2)*prod(diag(R)));
end

function cls=classifyMixture(G,Y)
    gamma=G.gamma; [~,idx]=max(gamma,[],2);
    names=["GMM_relaxed","GMM_intermediate","GMM_disturbed"];
    cls=names(idx)';
end


function Tq = classifyProfileFitQuality(Tfeat)
    Tq = table();
    Tq.cluster_id = Tfeat.cluster_id;
    Tq.hst_fit_flag = fitFlag(Tfeat.hst_n_bins,Tfeat.hst_sigma0,Tfeat.hst_rc,Tfeat.hst_eta);
    Tq.gaia_fit_flag = fitFlag(Tfeat.gaia_n_bins,Tfeat.gaia_sigma0,Tfeat.gaia_rc,Tfeat.gaia_eta);
    Tq.vel_fit_flag = fitFlag(Tfeat.vel_n_bins,Tfeat.vel_sigma0,Tfeat.vel_rc,Tfeat.vel_eta);
    Tq.n_problematic_profile_fits = double(Tq.hst_fit_flag~="valid_fit") + ...
        double(Tq.gaia_fit_flag~="valid_fit") + double(Tq.vel_fit_flag~="valid_fit");
end

function flag = fitFlag(nBins,sigma0,rc,eta)
    flag = strings(numel(nBins),1);
    for i=1:numel(nBins)
        if nBins(i) < 4 || ~isfinite(sigma0(i)) || ~isfinite(rc(i)) || ~isfinite(eta(i))
            flag(i)="insufficient_or_failed";
        elseif sigma0(i) <= 0 || rc(i) <= 0 || eta(i) <= 0
            flag(i)="nonphysical_fit";
        elseif eta(i) > 5 || rc(i) < 1e-3 || rc(i) > 1e5
            flag(i)="boundary_or_extreme_fit";
        else
            flag(i)="valid_fit";
        end
    end
end

function Tsel = gmmModelSelection(Y,Klist,maxIter,tol)
    rows = {};
    n = size(Y,1); d = size(Y,2);
    for K = Klist
        G = fitGMM_EM(Y,K,maxIter,tol);
        logL = G.logLik;
        nParams = (K-1) + K*d + K*d*(d+1)/2;
        BIC = -2*logL + nParams*log(n);
        AIC = -2*logL + 2*nParams;
        rows(end+1,:) = {K,logL,nParams,AIC,BIC}; %#ok<AGROW>
    end
    Tsel = cell2table(rows,'VariableNames',{'K','logLik','nParams','AIC','BIC'});
    Tsel = sortrows(Tsel,'BIC','ascend');
end

function T=joinValidationTables(clusterID,T_main,T_gmem,T_apPar,T_apChem,Tfeat)
    T=table(); T.cluster_id=clusterID;
    src={T_main,T_gmem,T_apPar,T_apChem,Tfeat};
    for s=1:numel(src)
        A=src{s}; if isempty(A)||~ismember("cluster_id",string(A.Properties.VariableNames)), continue; end
        A.cluster_id=string(A.cluster_id); vars=string(A.Properties.VariableNames);
        for v=vars
            if v=="cluster_id"||ismember(v,string(T.Properties.VariableNames)), continue; end
            vals=nan(numel(clusterID),1);
            for i=1:numel(clusterID)
                idx=find(A.cluster_id==clusterID(i),1);
                if ~isempty(idx)
                    val=A.(v)(idx);
                    if isnumeric(val)||islogical(val), vals(i)=double(val); else, vals(i)=str2double(string(val)); end
                end
            end
            if any(isfinite(vals)), T.(v)=vals; end
        end
    end
end

function Tcorr=computeCorrelations(T,clusterID,y,candidates)
    rows={};
    for v=candidates
        if ~ismember(v,string(T.Properties.VariableNames)), continue; end
        x=T.(v); ok=isfinite(x)&isfinite(y); N=sum(ok);
        if N<5, continue; end
        [rp,pp]=pearsonLocal(x(ok),y(ok)); [rs,ps]=spearmanLocal(x(ok),y(ok));
        rows(end+1,:)={v,N,rp,pp,rs,ps,abs(rs)}; %#ok<AGROW>
    end
    if isempty(rows), Tcorr=table(); else
        Tcorr=cell2table(rows,'VariableNames',{'variable','N','pearson_r','pearson_p_approx','spearman_rho','spearman_p_approx','abs_spearman_rho'});
    end
end

function [r,p]=pearsonLocal(x,y)
    r=corrLocal(x,y); n=numel(x);
    if n>2 && abs(r)<1
        t=r*sqrt((n-2)/(1-r^2)); p=2*(1-studentTCDF(abs(t),n-2));
    else, p=NaN; end
end

function [rho,p]=spearmanLocal(x,y), [rho,p]=pearsonLocal(tiedRanking(x),tiedRanking(y)); end
function r=corrLocal(x,y)
    x=x(:); y=y(:); x=x-mean(x,'omitnan'); y=y-mean(y,'omitnan');
    r=sum(x.*y,'omitnan')/sqrt(sum(x.^2,'omitnan')*sum(y.^2,'omitnan'));
end

function F=studentTCDF(t,nu)
    try, F=1-0.5*betainc(nu/(nu+t^2),nu/2,0.5);
    catch, F=0.5*(1+erf(t/sqrt(2))); end
end

function q=bhFDR(p)
    p=p(:); q=nan(size(p)); ok=isfinite(p); ps=p(ok); [sp,ord]=sort(ps); m=numel(ps);
    if m==0, return; end
    qs=sp.*m./(1:m)'; for i=m-1:-1:1, qs(i)=min(qs(i),qs(i+1)); end
    qo=nan(size(ps)); qo(ord)=min(qs,1); q(ok)=qo;
end

function readiness=computeReadiness(Q,err,bias,rankSig,loo,pnull)
    eS=exp(-max(0,err)); bS=exp(-abs(bias)/0.6); rS=exp(-rankSig/5); lS=exp(-loo/0.75);
    nS=min(1,-log10(max(pnull,1e-12))/2);
    readiness=100*(0.30*Q+0.20*eS+0.15*bS+0.15*rS+0.10*lS+0.10*nS);
end

function cls=classifyReadiness(x)
    cls=strings(numel(x),1);
    for i=1:numel(x)
        if x(i)>=80, cls(i)="A_publication_ready";
        elseif x(i)>=65, cls(i)="B_usable_with_caution";
        elseif x(i)>=50, cls(i)="C_candidate_only";
        else, cls(i)="D_low_confidence"; end
    end
end

function z=standardizeVector(x), z=(x-median(x,'omitnan'))/(1.4826*median(abs(x-median(x,'omitnan')),'omitnan')); end
function r=tiedRanking(x)
    x=x(:); [xs,ord]=sort(x,'ascend','MissingPlacement','last'); r=nan(size(x)); i=1;
    while i<=numel(x)
        if isnan(xs(i)), r(ord(i:end))=NaN; break; end
        j=i; while j<numel(x)&&xs(j+1)==xs(i), j=j+1; end
        r(ord(i:j))=(i+j)/2; i=j+1;
    end
end

function p=prctileLocal(X,q,dim)
    if dim==2
        p=nan(size(X,1),1); for i=1:size(X,1), p(i)=percentile1d(X(i,:),q); end
    else
        p=nan(1,size(X,2)); for j=1:size(X,2), p(j)=percentile1d(X(:,j),q); end
    end
end

function p=percentile1d(x,q)
    x=sort(x(isfinite(x))); if isempty(x), p=NaN; return; end
    if numel(x)==1, p=x; return; end
    pos=1+(q/100)*(numel(x)-1); lo=floor(pos); hi=ceil(pos);
    if lo==hi, p=x(lo); else, p=x(lo)+(pos-lo)*(x(hi)-x(lo)); end
end

function ord=orderBy(x,mode)
    if strcmpi(mode,'descend'), [~,ord]=sort(x,'descend'); else, [~,ord]=sort(x,'ascend'); end
end

function saveFig(fig,dir,name)
    name=char(regexprep(string(name),'[^\w\-]','_'));
    print(fig,fullfile(dir,[name '.png']),'-dpng','-r300');
    try, exportgraphics(fig,fullfile(dir,[name '.pdf']),'ContentType','vector');
    catch, print(fig,fullfile(dir,[name '.pdf']),'-dpdf','-painters'); end
end

function addLabels(x,y,lab)
    for i=1:numel(x)
        if isfinite(x(i))&&isfinite(y(i)), text(x(i),y(i)," "+lab(i),'FontSize',8,'Interpreter','none'); end
    end
end


function addLabelsTopN(x,y,lab,score,N)
    x = double(x(:)); y = double(y(:)); score = double(score(:));
    ok = isfinite(x) & isfinite(y) & isfinite(score);
    if sum(ok)==0, return; end
    idxAll = find(ok);
    [~,ord] = sort(abs(score(ok)),'descend');
    idx = idxAll(ord(1:min(N,numel(ord))));
    for ii = 1:numel(idx)
        i = idx(ii);
        text(x(i),y(i)," "+lab(i),'FontSize',7,'Interpreter','none', ...
            'HorizontalAlignment','left','VerticalAlignment','bottom');
    end
end

function scatterValidation(T,clusterID,y,var,xlab,figDir,name)
    if ~ismember(var,string(T.Properties.VariableNames)), return; end
    x=T.(var); ok=isfinite(x)&isfinite(y); if sum(ok)<5, return; end
    fig=figure('Color','w','Position',[100 100 850 650]);
    scatter(x(ok),y(ok),80,y(ok),'filled'); hold on; addLabelsTopN(x(ok),y(ok),clusterID(ok),y(ok),8);
    xlabel(xlab); ylabel('DSI_{stability,BC}'); title("External validation: "+var);
    cb=colorbar; cb.Label.String='DSI'; grid on; box on; saveFig(fig,figDir,name); close(fig);
end

function fig=plotProfilePanel(cid,HST,GAIA,VEL)
    fig=figure('Color','w','Position',[100 100 1100 760]); tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
    H=HST(HST.cluster_id==cid,:); G=GAIA(GAIA.cluster_id==cid & GAIA.radius_arcsec>0,:); V=VEL(VEL.cluster_id==cid,:);
    nexttile; if ~isempty(H), errorbar(H.radius_arcsec,H.pm_sigma_total_masyr,H.pm_sigma_total_err_masyr,'o-'); set(gca,'XScale','log'); grid on; xlabel('R'); ylabel('HST \sigma_{PM}'); title(cid+" HST total"); end
    nexttile; if ~isempty(H), plot(H.radius_arcsec,H.pm_sigma_radial_masyr,'o-'); hold on; plot(H.radius_arcsec,H.pm_sigma_tangential_masyr,'s-'); set(gca,'XScale','log'); grid on; legend({'radial','tangential'}); title(cid+" HST anisotropy"); end
    nexttile; if ~isempty(G), plot(G.radius_arcsec,G.pm_disp_p50_masyr,'-'); hold on; plot(G.radius_arcsec,abs(G.pm_rot_p50_masyr),'--'); set(gca,'XScale','log'); grid on; legend({'dispersion','|rotation|'}); title(cid+" Gaia"); end
    nexttile; if ~isempty(V), types=unique(string(V.data_type)); hold on; for k=1:numel(types), Vk=V(string(V.data_type)==types(k),:); errorbar(Vk.radius_arcsec,Vk.sigma_kms,Vk.sigma_err_up_kms,'o','DisplayName',types(k)); end; set(gca,'XScale','log'); grid on; legend('Location','best'); title(cid+" RV/PM"); end
    sgtitle("v5 profile diagnostics: "+cid,'Interpreter','none');
end


function [T, Tmap] = addStructuralDynamicalValidationVariables(T, clusterID)
    % Adds standardized validation variables from Harris/Baumgardt columns.
    % These are external validation variables only.
    n = numel(clusterID);
    cAlias = ["concentration_c","concentration","c","C","c_harris", ...
              "harris_c","King_c","king_c","conc"];
    trcAlias = ["log10_trc_yr","log_trc","log10_trc","log_t_rc", ...
                "log10_t_rc","log_rc_relax","log_core_relaxation_time", ...
                "trc_yr","t_rc_yr","core_relaxation_time_yr"];
    trhAlias = ["log10_trh_yr","log_trh","log10_trh","log_t_rh", ...
                "log10_t_rh","log_rh_relax","log_half_mass_relaxation_time", ...
                "trh_yr","t_rh_yr","half_mass_relaxation_time_yr"];
    ageAlias = ["age_gyr","Age_Gyr","age","cluster_age_gyr","Age"];
    ccAlias = ["core_collapse_flag","core_collapsed","core_collapse", ...
               "is_core_collapsed","post_core_collapse","PCC","cc","CC", ...
               "core_status","collapse_status","classification","core_collapse_status"];

    [cVal, cSrc] = getNumericByAliases(T,cAlias);
    [trcRaw, trcSrc] = getNumericByAliases(T,trcAlias);
    [trhRaw, trhSrc] = getNumericByAliases(T,trhAlias);
    [ageGyr, ageSrc] = getNumericByAliases(T,ageAlias);
    [ccFlag, ccSrc] = getCoreCollapseFlagByAliases(T,ccAlias);

    if any(isfinite(cVal)), T.concentration_c = cVal; else, T.concentration_c = nan(n,1); end
    [T.log10_trc_yr, T.trc_yr] = standardizeRelaxationTime(trcRaw, trcSrc);
    [T.log10_trh_yr, T.trh_yr] = standardizeRelaxationTime(trhRaw, trhSrc);

    if any(isfinite(ageGyr))
        ageYr = ageGyr .* 1e9;
        T.dyn_age_rc = ageYr ./ T.trc_yr;
        T.dyn_age_rh = ageYr ./ T.trh_yr;
    else
        T.dyn_age_rc = nan(n,1);
        T.dyn_age_rh = nan(n,1);
    end

    if any(isfinite(ccFlag))
        T.core_collapse_flag = ccFlag;
        ccMode = "explicit_column";
    elseif any(isfinite(T.concentration_c))
        T.core_collapse_flag = double(T.concentration_c >= 2.49);
        ccMode = "derived_from_high_concentration_c_ge_2p49";
    else
        T.core_collapse_flag = nan(n,1);
        ccMode = "not_available";
    end

    Tmap = table( ...
        ["concentration_c";"log10_trc_yr";"log10_trh_yr";"age_gyr";"core_collapse_flag"], ...
        [cSrc;trcSrc;trhSrc;ageSrc;ccSrc], ...
        ["numeric";"numeric/log_standardized";"numeric/log_standardized";"numeric";ccMode], ...
        'VariableNames',{'standard_variable','source_column','mapping_note'});

    fprintf('\nStructural validation variables:\n');
    fprintf('  concentration_c source : %s\n', cSrc);
    fprintf('  log10_trc_yr source    : %s\n', trcSrc);
    fprintf('  log10_trh_yr source    : %s\n', trhSrc);
    fprintf('  core-collapse source   : %s (%s)\n', ccSrc, ccMode);
end

function [x, src] = getNumericByAliases(T, aliases)
    n = height(T); x = nan(n,1); src = "not_found";
    vars = string(T.Properties.VariableNames);
    for a = aliases
        idx = find(strcmpi(vars,a),1);
        if ~isempty(idx)
            vname = vars(idx); raw = T.(vname);
            if isnumeric(raw) || islogical(raw), x = double(raw);
            else, x = str2double(string(raw)); end
            src = vname; return;
        end
    end
end

function [flag, src] = getCoreCollapseFlagByAliases(T, aliases)
    n = height(T); flag = nan(n,1); src = "not_found";
    vars = string(T.Properties.VariableNames);
    for a = aliases
        idx = find(strcmpi(vars,a),1);
        if isempty(idx), continue; end
        vname = vars(idx); raw = T.(vname); src = vname;
        if isnumeric(raw) || islogical(raw)
            flag = double(raw); flag(flag~=0 & flag~=1) = nan; return;
        end
        s = lower(strtrim(string(raw))); f = nan(n,1);
        yesPat = ["cc","core collapsed","core-collapsed","collapsed", ...
                  "pcc","post core collapse","post-core-collapse","yes","true","1"];
        noPat = ["non core collapsed","non-core-collapsed","not collapsed", ...
                 "normal","king","no","false","0"];
        for i=1:n
            if any(contains(s(i),yesPat)), f(i) = 1;
            elseif any(contains(s(i),noPat)), f(i) = 0; end
        end
        flag = f; return;
    end
end

function [logt, t] = standardizeRelaxationTime(x, src)
    x = double(x(:)); logt = nan(size(x)); t = nan(size(x));
    if ~any(isfinite(x)), return; end
    srcLow = lower(string(src));
    if contains(srcLow,"log")
        logt = x; t = 10.^logt; return;
    end
    medx = median(x,'omitnan');
    if medx > 4 && medx < 12
        logt = x; t = 10.^logt;
    elseif medx > 1e4
        t = x; logt = log10(t);
    elseif medx > 0 && medx < 1e4
        t = x .* 1e6; logt = log10(t);
    end
end

function Ttest = compareCoreCollapseGroups(T, clusterID, y, yName)
    %#ok<INUSD>
    if ~ismember("core_collapse_flag", string(T.Properties.VariableNames))
        Ttest = table(string(yName),0,0,nan,nan,nan,nan,nan,"core_collapse_flag_not_available", ...
            'VariableNames',{'metric','N_nonCC','N_CC','median_nonCC','median_CC', ...
            'delta_median_CC_minus_nonCC','p_perm','cliffs_delta_CC_vs_nonCC','note'});
        return;
    end
    cc = T.core_collapse_flag;
    ok = isfinite(y) & isfinite(cc) & (cc==0 | cc==1);
    y0 = y(ok & cc==0); y1 = y(ok & cc==1);
    if numel(y0)<3 || numel(y1)<2
        Ttest = table(string(yName),numel(y0),numel(y1),median(y0,'omitnan'),median(y1,'omitnan'),nan,nan,nan, ...
            "insufficient_group_size", ...
            'VariableNames',{'metric','N_nonCC','N_CC','median_nonCC','median_CC', ...
            'delta_median_CC_minus_nonCC','p_perm','cliffs_delta_CC_vs_nonCC','note'});
        return;
    end
    deltaObs = median(y1,'omitnan') - median(y0,'omitnan');
    pPerm = permutationGroupP(y0,y1,20000);
    cd = cliffsDelta(y1,y0);
    Ttest = table(string(yName),numel(y0),numel(y1),median(y0,'omitnan'),median(y1,'omitnan'),deltaObs,pPerm,cd, ...
        "CC_vs_nonCC_permutation_test", ...
        'VariableNames',{'metric','N_nonCC','N_CC','median_nonCC','median_CC', ...
        'delta_median_CC_minus_nonCC','p_perm','cliffs_delta_CC_vs_nonCC','note'});
end

function p = permutationGroupP(y0,y1,Nperm)
    y0 = y0(:); y1 = y1(:); y = [y0;y1]; n1 = numel(y1);
    deltaObs = median(y1,'omitnan') - median(y0,'omitnan'); cnt = 0;
    for b=1:Nperm
        idx = randperm(numel(y)); y1b = y(idx(1:n1)); y0b = y(idx(n1+1:end));
        deltaB = median(y1b,'omitnan') - median(y0b,'omitnan');
        if abs(deltaB) >= abs(deltaObs), cnt = cnt + 1; end
    end
    p = (1+cnt)/(Nperm+1);
end

function d = cliffsDelta(y1,y0)
    y1 = y1(:); y0 = y0(:);
    if isempty(y1) || isempty(y0), d = NaN; return; end
    gt = 0; lt = 0;
    for i=1:numel(y1)
        gt = gt + sum(y1(i) > y0); lt = lt + sum(y1(i) < y0);
    end
    d = (gt - lt) / (numel(y1)*numel(y0));
end

function plotCoreCollapseBox(T,clusterID,y,figDir,outName)
    if ~ismember("core_collapse_flag", string(T.Properties.VariableNames)), return; end
    cc = T.core_collapse_flag; ok = isfinite(y) & isfinite(cc) & (cc==0 | cc==1);
    if sum(ok)<5, return; end
    fig = figure('Color','w','Position',[100 100 850 650]);
    boxplot(y(ok), categorical(cc(ok),[0 1],{'non-CC','CC/PCC'})); hold on;
    xcat = double(categorical(cc(ok),[0 1],{'non-CC','CC/PCC'}));
    scatter(xcat, y(ok), 70, 'filled', 'MarkerFaceAlpha',0.70);
    xlabel('Core-collapse status'); ylabel('DSI_{consensus,BC}');
    title('Consensus DSI by core-collapse status');
    addLabelsTopN(xcat,y(ok),clusterID(ok),y(ok),10);
    grid on; box on; saveFig(fig,figDir,outName); close(fig);
end
