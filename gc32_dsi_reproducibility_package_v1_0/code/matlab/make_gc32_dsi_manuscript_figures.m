function make_gc32_dsi_manuscript_figures(tablesDir, structuralFile, outDir, showFigures)
% make_gc32_dsi_manuscript_figures
% -------------------------------------------------------------------------
% Manuscript figure generator for the GC32 DSI analysis.
%
% This script creates the main figures recommended for the manuscript text:
%   Fig. 1  Primary DSI ranking with uncertainty
%   Fig. 2  Consensus DSI vs classical/Bayesian latent scores
%   Fig. 3  Null-family significance map
%   Fig. 4  Component heatmap
%   Fig. 5  Learned weights and component contributions
%   Fig. 6  Gaia lambda_R-like spin proxy vs DSI
%   Fig. 7  Structural validation: DSI vs relaxation time and concentration
%   Fig. 8  Core-collapse status diagnostic
%   Fig. S1 GMM diagnostic classes
%
% Expected input files are the archived or recomputed CSV outputs from the MATLAB DSI
% pipeline, e.g.:
%   gc32_dsi_CONSENSUS_ranked_BC.csv
%   gc32_dsi_PRIMARY_ranked_DSI_stability_BC.csv
%   gc32_dsi_robust_z_components_capped.csv
%   gc32_dsi_learned_weights_and_latent_loadings.csv
%   gc32_dsi_null_model_significance.csv
%   gc32_dsi_external_validation_correlations_DSI.csv
%   gc32_dsi_GMM_model_selection_BIC.csv
%   13_gc32_structural_validation.csv
%
% Usage examples:
%   make_gc32_dsi_manuscript_figures
%   make_gc32_dsi_manuscript_figures('tables','13_gc32_structural_validation.csv','figures_recomputed',true)
%
% Notes:
%   - The script is intentionally conservative: it uses existing pipeline
%     results and does not recompute the DSI.
%   - It exports both vector PDF and 600 dpi PNG files.
%   - It uses clean fonts, consistent sizes, annotation of only key clusters,
%     and avoids overcrowding the figures.
% -------------------------------------------------------------------------

if nargin < 1 || isempty(tablesDir)
    scriptDir = fileparts(mfilename('fullpath'));
    packageRoot = fileparts(fileparts(scriptDir));
    if isfolder(fullfile(packageRoot,'outputs','tables_csv'))
        tablesDir = fullfile(packageRoot,'outputs','tables_csv');
    elseif isfolder(fullfile(packageRoot,'outputs','recomputed','tables'))
        tablesDir = fullfile(packageRoot,'outputs','recomputed','tables');
    elseif isfolder(fullfile(pwd,'tables'))
        tablesDir = fullfile(pwd,'tables');
    else
        tablesDir = pwd;
    end
end
if nargin < 2 || isempty(structuralFile)
    structuralFile = findStructuralValidationFile(tablesDir);
else
    structuralFile = char(structuralFile);
    if ~isfile(structuralFile)
        candidate = fullfile(pwd,structuralFile);
        if isfile(candidate)
            structuralFile = candidate;
        else
            structuralFile = findStructuralValidationFile(tablesDir);
        end
    end
end
if nargin < 3 || isempty(outDir)
    scriptDir = fileparts(mfilename('fullpath'));
    packageRoot = fileparts(fileparts(scriptDir));
    outDir = fullfile(packageRoot,'outputs','recomputed','figures');
end
if nargin < 4 || isempty(showFigures)
    showFigures = true;   % interactive display by default
end

if ~isfolder(outDir), mkdir(outDir); end
setManuscriptStyle(showFigures);
setappdata(0,'GC32_DISPLAY_FIGURES',showFigures);
if showFigures
    fprintf('Interactive mode: figures will remain visible and open after export.\n');
else
    fprintf('Export-only mode: figures are hidden and closed after export.\n');
end

% -------------------------------------------------------------------------
% Load tables. File patterns are intentionally tolerant because the current
% package uses v5_0 names for core results and v5_1 names for structural
% validation tables.
% -------------------------------------------------------------------------
Tcons  = readFirstTable(tablesDir, {'S2_table_02_consensus_ranking.csv','*consensus*ranked*bc*.csv','*CONSENSUS_ranked_BC.csv','*consensus*BC*.csv'});
Tprim  = readFirstTable(tablesDir, {'S2_table_01_primary_stability_ranking.csv','*primary*stability*bc*.csv','*PRIMARY_ranked_DSI_stability_BC.csv'});
Tz     = readFirstTable(tablesDir, {'S2_table_08_robust_z_components_capped.csv','*robust_z_components_capped.csv','*z_components*capped*.csv'});
Traw   = readOptionalTable(tablesDir, {'S2_table_07_raw_dynamical_components.csv','*raw_components*.csv'});
Twt    = readFirstTable(tablesDir, {'S2_table_10_weights_and_latent_loadings.csv','*learned_weights_and_latent_loadings.csv','*weights*loadings*.csv'});
Tnull  = readFirstTable(tablesDir, {'S2_table_13_null_family_significance.csv','*null*significance*.csv'});
TphysN = readOptionalTable(tablesDir, {'S2_table_14_component_space_relaxed_proxy_null.csv','*component*space*relaxed*proxy*null*.csv','*physical*null*.csv'});
Tcorr  = readOptionalTable(tablesDir, {'S2_table_15_external_validation_correlations_DSI.csv','*external*validation*DSI*.csv'});
Tgmm   = readOptionalTable(tablesDir, {'S2_table_24_gmm_model_selection.csv','*gmm*model*selection*.csv','*GMM*BIC*.csv'});

Tstruct = table();
if ~isempty(structuralFile) && isfile(structuralFile)
    fprintf('Structural validation file: %s\n', structuralFile);
    Tstruct = readCsvSafe(structuralFile);
else
    fprintf('Structural validation file not found. Fig. 7 and Fig. 8 will be skipped.\n');
end

% Ensure key string columns.
Tcons.cluster_id = string(Tcons.cluster_id);
Tprim.cluster_id = string(Tprim.cluster_id);
Tz.cluster_id    = string(Tz.cluster_id);
Twt.component    = string(Twt.component);
Tnull.cluster_id = string(Tnull.cluster_id);
if ~isempty(Tstruct), Tstruct.cluster_id = string(Tstruct.cluster_id); end
if ~isempty(Traw), Traw.cluster_id = string(Traw.cluster_id); end

% Sort key tables.
Tcons = sortrows(Tcons,'DSI_consensus_BC','descend');
Tprim = sortrows(Tprim,'DSI_stability_BC','descend');

% Create figures.
fprintf('Generating manuscript figures in: %s\n', outDir);
fig01_primaryRanking(Tprim,outDir);
fig02_latentConsistency(Tcons,outDir);
fig03_nullFamilyMap(Tcons,Tnull,outDir);
fig04_componentHeatmap(Tcons,Tz,outDir);
fig05_weightsAndContributions(Tcons,Tz,Twt,outDir);
fig06_lambdaRValidation(Tcons,Tcorr,outDir);
fig07_structuralValidation(Tcons,Tstruct,outDir);
fig08_coreCollapseDiagnostic(Tcons,Tstruct,outDir);
figS1_gmmDiagnostic(Tcons,Tgmm,outDir);

fprintf('Done. Figures exported as PDF and PNG.\n');
end

% =========================================================================
% Figure 1: DSI ranking
% =========================================================================
function fig01_primaryRanking(T,outDir)
fig = newFig([100 100 1050 850]);
T = sortrows(T,'DSI_stability_BC','ascend');
N = height(T);
y = 1:N;
x = getNum(T,'DSI_stability_BC');
xe = getNum(T,'DSI_stability_err');
qC = getNum(T,'q_null_combined_conservative');
Q  = getNum(T,'data_quality_score');

hold on;
for i=1:N
    c = qualityColor(Q(i));
    barh(y(i),x(i),0.72,'FaceColor',c,'EdgeColor','none','FaceAlpha',0.92);
end
errorbar(x,y,xe,'horizontal','k.','LineWidth',0.75,'CapSize',4);
plot([0 0],[0 N+1],'k--','LineWidth',0.8);

sig = qC < 0.05;
plot(x(sig),y(sig),'p','MarkerSize',11,'MarkerFaceColor',[0.10 0.10 0.10],...
    'MarkerEdgeColor','w','LineWidth',0.8);

yticks(y);
yticklabels(T.cluster_id);
xlabel('Primary stability-weighted DSI, DSI_{stab,BC}');
ylabel('Cluster');
title('GC32 v5.1 primary dynamical-complexity ranking');
subtitle('Bars are coloured by data-quality score; pentagons mark conservative combined-null significance');
grid on; box on;
set(gca,'YDir','normal','TickLabelInterpreter','none');

% annotate selected high-quality candidates
key = ismember(T.cluster_id, ["NGC5053","NGC5466","NGC6093","NGC2298","NGC5272"]);
for i=find(key)'
    text(x(i)+0.06*rangeFinite(x), y(i), char(T.cluster_id(i)), 'FontSize',8, ...
        'FontWeight','bold','VerticalAlignment','middle');
end

exportFig(fig,outDir,'fig01_v5_1_primary_DSI_stability_BC_ranked');
maybeCloseFig(fig);
end

% =========================================================================
% Figure 2: DSI vs latent scores
% =========================================================================
function fig02_latentConsistency(T,outDir)
fig = newFig([100 100 1120 520]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

x = getNum(T,'DSI_consensus_BC');
q = getNum(T,'q_null_combined_conservative');
Q = getNum(T,'data_quality_score');

nexttile;
y = getNum(T,'theta_latent_BC');
ye = getNum(T,'theta_latent_err');
scatterPanel(x,y,[],ye,q,Q,T.cluster_id,'Classical latent score, \theta_{class,BC}');
title('Consensus DSI vs classical latent score');

nexttile;
y = getNum(T,'theta_bayes_BC');
ye = getNum(T,'theta_bayes_err');
scatterPanel(x,y,[],ye,q,Q,T.cluster_id,'Bayesian latent score, \theta_{Bayes,BC}');
title('Consensus DSI vs Bayesian latent score');

xlabel(tl,'Consensus DSI, DSI_{cons,BC}');
exportFig(fig,outDir,'fig02_v5_1_DSI_vs_classical_and_bayesian_latent');
maybeCloseFig(fig);
end

% =========================================================================
% Figure 3: Null-family significance map
% =========================================================================
function fig03_nullFamilyMap(Tcons,Tnull,outDir)
fig = newFig([100 100 1100 760]);
T = innerjoin(Tcons(:,{'cluster_id','DSI_consensus_BC'}), Tnull, 'Keys','cluster_id');
T = sortrows(T,'DSI_consensus_BC','descend');
T = T(1:min(12,height(T)),:);

qVars = {'q_null_DSI','q_null_theta','q_null_theta_bayes','q_null_consensus',...
    'q_null_physical_proxy','q_null_rotation_free','q_null_profile_shuffle','q_null_combined_conservative'};
labels = {'DSI','Classical \theta','Bayes \theta','Consensus','Physical proxy','Rotation-free','Profile-shuffle','Combined'};
M = nan(height(T),numel(qVars));
for j=1:numel(qVars)
    if ismember(qVars{j},T.Properties.VariableNames)
        M(:,j) = -log10(max(getNum(T,qVars{j}), realmin));
    end
end
M(M>6) = 6;

imagesc(M);
colormap(parula(256));
cb = colorbar; cb.Label.String = '-log_{10}(q), capped at 6';
set(gca,'XTick',1:numel(labels),'XTickLabel',labels,'XTickLabelRotation',35,...
    'YTick',1:height(T),'YTickLabel',T.cluster_id,'TickLabelInterpreter','tex');
title('Null-family significance for high-ranked clusters');
subtitle('The conservative combined-null column is the strictest significance layer');

% annotate q<0.05 cells
for i=1:size(M,1)
    for j=1:size(M,2)
        if M(i,j) > -log10(0.05)
            text(j,i,'*','Color','w','FontSize',14,'FontWeight','bold',...
                'HorizontalAlignment','center','VerticalAlignment','middle');
        end
    end
end

exportFig(fig,outDir,'fig03_v5_1_null_family_significance_map');
maybeCloseFig(fig);
end

% =========================================================================
% Figure 4: Component heatmap
% =========================================================================
function fig04_componentHeatmap(Tcons,Tz,outDir)
fig = newFig([100 100 1020 820]);
components = {'GradientStrength','Anisotropy','RotationSupport','ProfileCurvature','CentralEnhancement','OuterDisturbance'};
labels = {'G','A','R','C','E','T'};
T = innerjoin(Tcons(:,{'cluster_id','DSI_consensus_BC'}), Tz, 'Keys','cluster_id');
T = sortrows(T,'DSI_consensus_BC','descend');
M = nan(height(T),numel(components));
for j=1:numel(components)
    M(:,j)=getNum(T,components{j});
end
M = max(min(M,3.5),-3.5);

imagesc(M);
colormap(redblue(256));
caxis([-3.5 3.5]);
cb=colorbar; cb.Label.String='Capped robust z';
set(gca,'XTick',1:numel(labels),'XTickLabel',labels,'YTick',1:height(T),...
    'YTickLabel',T.cluster_id,'TickLabelInterpreter','none');
title('Capped robust-z component matrix sorted by consensus DSI');
subtitle('G: gradient, A: anisotropy, R: rotation, C: curvature, E: central enhancement, T: outer disturbance');

exportFig(fig,outDir,'fig04_v5_1_component_heatmap');
maybeCloseFig(fig);
end

% =========================================================================
% Figure 5: Learned weights and component contributions
% =========================================================================
function fig05_weightsAndContributions(Tcons,Tz,Twt,outDir)
fig = newFig([100 100 1180 560]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
components = {'GradientStrength','Anisotropy','RotationSupport','ProfileCurvature','CentralEnhancement','OuterDisturbance'};
shortLabels = {'G','A','R','C','E','T'};

nexttile;
w = nan(1,numel(components));
for j=1:numel(components)
    idx = strcmpi(Twt.component,components{j});
    if any(idx), w(j)=Twt.w_stability(find(idx,1)); end
end
bar(1:numel(w),w,0.68,'FaceColor',[0.25 0.45 0.70],'EdgeColor','none');
set(gca,'XTick',1:numel(w),'XTickLabel',shortLabels);
ylabel('Stability weight');
title('Learned stability weights');
grid on; box on;
for j=1:numel(w)
    text(j,w(j)+0.015,sprintf('%.3f',w(j)),'HorizontalAlignment','center','FontSize',8);
end

nexttile;
T = innerjoin(Tcons(:,{'cluster_id','DSI_consensus_BC'}), Tz, 'Keys','cluster_id');
T = sortrows(T,'DSI_consensus_BC','descend');
T = T(1:min(8,height(T)),:);
Z = nan(height(T),numel(components));
for j=1:numel(components), Z(:,j)=getNum(T,components{j}); end
C = Z .* w;
bar(C,'stacked','EdgeColor','none');
set(gca,'XTick',1:height(T),'XTickLabel',T.cluster_id,'XTickLabelRotation',45,'TickLabelInterpreter','none');
ylabel('Weighted component contribution, z_j w_j');
title('High-ranked cluster component contributions');
grid on; box on;
legend(shortLabels,'Location','eastoutside');

exportFig(fig,outDir,'fig05_v5_1_weights_and_component_contributions');
maybeCloseFig(fig);
end

% =========================================================================
% Figure 6: lambda_R validation
% =========================================================================
function fig06_lambdaRValidation(T,Tcorr,outDir)
fig = newFig([100 100 760 660]);
x = getNum(T,'lambda_R_gaia');
y = getNum(T,'DSI_consensus_BC');
Q = getNum(T,'data_quality_score');
q = getNum(T,'q_null_combined_conservative');
scatterValidationBase(x,y,Q,q,T.cluster_id);
xlabel('Gaia \lambda_R-like spin proxy');
ylabel('Consensus DSI, DSI_{cons,BC}');
title('External validation: angular-momentum support');

[rho,pval] = spearmanSimple(x,y);
qtxt = '';
if ~isempty(Tcorr) && any(strcmpi(string(Tcorr.variable),'lambda_R_gaia'))
    row = find(strcmpi(string(Tcorr.variable),'lambda_R_gaia'),1);
    qtxt = sprintf(', q=%.2g', Tcorr.BH_q_spearman(row));
end
addStatsText(sprintf('Spearman \rho=%.3f, p=%.2g%s',rho,pval,qtxt));
exportFig(fig,outDir,'fig06_v5_1_lambdaR_vs_DSI');
maybeCloseFig(fig);
end

% =========================================================================
% Figure 7: Structural validation
% =========================================================================
function fig07_structuralValidation(T,Tstruct,outDir)
if isempty(Tstruct)
    warning('Structural validation file not found. Fig. 7 skipped.');
    return;
end
fig = newFig([100 100 1120 520]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
TJ = innerjoin(T(:,{'cluster_id','DSI_consensus_BC','data_quality_score','q_null_combined_conservative'}), Tstruct, 'Keys','cluster_id');

nexttile;
x = getNum(TJ,'log10_trh_yr');
y = getNum(TJ,'DSI_consensus_BC');
Q = getNum(TJ,'data_quality_score');
q = getNum(TJ,'q_null_combined_conservative');
scatterValidationBase(x,y,Q,q,TJ.cluster_id);
xlabel('log_{10}(t_{rh}/yr)');
ylabel('Consensus DSI, DSI_{cons,BC}');
title('Half-mass relaxation time');
[rho,pval] = spearmanSimple(x,y);
addStatsText(sprintf('\rho=%.3f, p=%.2g; marginal after FDR',rho,pval));

nexttile;
x = getNum(TJ,'concentration_c');
y = getNum(TJ,'DSI_consensus_BC');
scatterValidationBase(x,y,Q,q,TJ.cluster_id);
xlabel('Concentration, c');
ylabel('Consensus DSI, DSI_{cons,BC}');
title('Structural concentration');
[rho,pval] = spearmanSimple(x,y);
addStatsText(sprintf('\rho=%.3f, p=%.2g; not a concentration proxy',rho,pval));

exportFig(fig,outDir,'fig07_v5_1_structural_validation_relaxation_concentration');
maybeCloseFig(fig);
end

% =========================================================================
% Figure 8: Core-collapse diagnostic
% =========================================================================
function fig08_coreCollapseDiagnostic(T,Tstruct,outDir)
if isempty(Tstruct) || ~ismember('core_collapse_flag', Tstruct.Properties.VariableNames)
    warning('Core-collapse flag not found. Fig. 8 skipped.');
    return;
end
TJ = innerjoin(T(:,{'cluster_id','DSI_consensus_BC','data_quality_score'}), Tstruct, 'Keys','cluster_id');
cc = getNum(TJ,'core_collapse_flag');
y = getNum(TJ,'DSI_consensus_BC');
ok = isfinite(cc) & isfinite(y) & (cc==0 | cc==1);
if sum(ok)<5, return; end

fig = newFig([100 100 760 620]);
hold on;
for g = 0:1
    yg = y(ok & cc==g);
    xg = (g+1)*ones(size(yg));
    q1 = quantileNoToolbox(yg,0.25); q2 = median(yg,'omitnan'); q3 = quantileNoToolbox(yg,0.75);
    plot([g+0.75 g+1.25],[q2 q2],'k-','LineWidth',2.0);
    rectangle('Position',[g+0.85 q1 0.30 max(q3-q1,eps)],'EdgeColor','k','FaceColor',[0.80 0.86 0.95]);
    scatter(xg + 0.05*randn(size(xg)), yg, 55, 'filled', 'MarkerFaceAlpha',0.75, ...
        'MarkerFaceColor', [0.25 0.45 0.70], 'MarkerEdgeColor','k');
end
set(gca,'XTick',[1 2],'XTickLabel',{'non-CC','CC/PCC'});
xlim([0.5 2.5]); ylabel('Consensus DSI, DSI_{cons,BC}');
title('Core-collapse status diagnostic');
subtitle('Diagnostic only: the CC/PCC group contains only two clusters in the strict sample');
grid on; box on;
for i=find(ok)'
    if cc(i)==1
        text(cc(i)+1+0.05, y(i), char(TJ.cluster_id(i)), 'FontSize',8, 'FontWeight','bold');
    end
end
exportFig(fig,outDir,'fig08_v5_1_core_collapse_status_diagnostic');
maybeCloseFig(fig);
end

% =========================================================================
% Supplementary GMM diagnostic
% =========================================================================
function figS1_gmmDiagnostic(T,Tgmm,outDir)
fig = newFig([100 100 760 650]);
x = getNum(T,'DSI_stability_BC');
y = getNum(T,'theta_latent_BC');
cls = string(T.mixture_class);
classes = unique(cls,'stable');
hold on;
cols = [0.30 0.55 0.80; 0.90 0.60 0.25; 0.70 0.35 0.75; 0.45 0.70 0.35];
for k=1:numel(classes)
    idx = cls==classes(k);
    scatter(x(idx),y(idx),80,'filled','MarkerFaceColor',cols(k,:),...
        'MarkerEdgeColor','k','DisplayName',char(classes(k)));
end
labelKeyClusters(x,y,T.cluster_id,T.DSI_consensus_BC,7);
xlabel('Stability-weighted DSI, DSI_{stab,BC}');
ylabel('Classical latent score, \theta_{class,BC}');
title('Diagnostic GMM classes in DSI--latent space');
if ~isempty(Tgmm)
    [~,ib] = min(getNum(Tgmm,'BIC'));
    txt = sprintf('BIC favours K=%d', round(Tgmm.K(ib)));
    addStatsText(txt);
end
legend('Location','best'); grid on; box on;
exportFig(fig,outDir,'figS1_v5_1_GMM_classes_diagnostic');
maybeCloseFig(fig);
end

% =========================================================================
% Plot helpers
% =========================================================================
function scatterPanel(x,y,xe,ye,q,Q,labels,yLabelText)
hold on;
scatterValidationBase(x,y,Q,q,labels);
if ~isempty(ye) && any(isfinite(ye))
    errorbar(x,y,ye,'vertical','Color',[0.30 0.30 0.30],'LineStyle','none','LineWidth',0.7,'CapSize',3);
end
if ~isempty(xe) && any(isfinite(xe))
    errorbar(x,y,xe,'horizontal','Color',[0.30 0.30 0.30],'LineStyle','none','LineWidth',0.7,'CapSize',3);
end
ylabel(yLabelText);
[rho,pval]=spearmanSimple(x,y);
addStatsText(sprintf('Spearman \rho=%.3f, p=%.2g',rho,pval));
end

function scatterValidationBase(x,y,Q,q,labels)
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok); Q = Q(ok); q = q(ok); labels = labels(ok);
if isempty(Q), Q = ones(size(x)); end
if isempty(q), q = nan(size(x)); end
sz = 45 + 55*max(min(Q,1),0);
scatter(x,y,sz,Q,'filled','MarkerEdgeColor',[0.08 0.08 0.08],'MarkerFaceAlpha',0.82);
colormap(parula(256)); cb=colorbar; cb.Label.String='Data-quality score';
hold on;
plotFitLine(x,y);
plot(x(q<0.05),y(q<0.05),'p','MarkerSize',12,'MarkerFaceColor',[0.05 0.05 0.05],...
    'MarkerEdgeColor','w','LineWidth',0.8);
labelKeyClusters(x,y,labels,y,6);
grid on; box on;
end

function labelKeyClusters(x,y,labels,score,N)
ok = isfinite(x) & isfinite(y) & isfinite(score);
idxAll = find(ok);
if isempty(idxAll), return; end
[~,ord] = sort(abs(score(ok)),'descend');
idx = idxAll(ord(1:min(N,numel(ord))));
for ii=1:numel(idx)
    i = idx(ii);
    text(x(i)+0.015*rangeFinite(x),y(i),char(labels(i)),...
        'FontSize',8,'FontWeight','bold','Interpreter','none');
end
end

function plotFitLine(x,y)
ok = isfinite(x)&isfinite(y);
if sum(ok)<3, return; end
p = polyfit(x(ok),y(ok),1);
xx = linspace(min(x(ok)),max(x(ok)),100);
yy = polyval(p,xx);
plot(xx,yy,'k-','LineWidth',1.2);
end

function addStatsText(txt)
xl = xlim; yl = ylim;
text(xl(1)+0.03*rangeFinite(xl), yl(2)-0.07*rangeFinite(yl), txt, ...
    'FontSize',9,'BackgroundColor','w','EdgeColor',[0.75 0.75 0.75],...
    'Margin',4,'Interpreter','tex');
end

function c = qualityColor(Q)
if ~isfinite(Q), Q=0.5; end
Q = max(0,min(1,Q));
lo = [0.86 0.43 0.37];
hi = [0.25 0.55 0.78];
c = lo*(1-Q) + hi*Q;
end

function C = redblue(m)
if nargin<1, m=256; end
x = linspace(-1,1,m)';
C = zeros(m,3);
C(:,1) = max(0,x);
C(:,3) = max(0,-x);
C(:,2) = 1 - abs(x);
C = 0.15 + 0.85*C;
end

% =========================================================================
% IO helpers
% =========================================================================
function T = readFirstTable(baseDir, patterns)
files = findFiles(baseDir, patterns);
if isempty(files)
    error('Required table not found. Directory: %s. Patterns: %s', baseDir, strjoin(patterns,', '));
end
T = readCsvSafe(files{1});
end

function T = readOptionalTable(baseDir, patterns)
files = findFiles(baseDir, patterns);
if isempty(files)
    T = table();
else
    T = readCsvSafe(files{1});
end
end

function files = findFiles(baseDir, patterns)
files = {};
for p = 1:numel(patterns)
    d = dir(fullfile(baseDir, patterns{p}));
    if isempty(d)
        d = dir(fullfile(baseDir,'**',patterns{p}));
    end
    for i=1:numel(d)
        if ~d(i).isdir
            files{end+1} = fullfile(d(i).folder,d(i).name); %#ok<AGROW>
        end
    end
    if ~isempty(files), return; end
end
end

function T = readCsvSafe(fileName)
try
    T = readtable(fileName,'VariableNamingRule','preserve');
catch
    T = readtable(fileName);
end
% Ensure MATLAB-safe variable names if preserve created problematic access.
T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames,'ReplacementStyle','delete');
end

function x = getNum(T,varName)
if isempty(T) || ~ismember(varName,T.Properties.VariableNames)
    x = nan(height(T),1);
    return;
end
x = T.(varName);
if iscell(x) || isstring(x) || ischar(x)
    x = str2double(string(x));
else
    x = double(x);
end
x = x(:);
end

% =========================================================================
% Statistics helpers without toolbox dependence
% =========================================================================
function [rho,p] = spearmanSimple(x,y)
ok = isfinite(x)&isfinite(y);
x=x(ok); y=y(ok);
if numel(x)<4
    rho=NaN; p=NaN; return;
end
rx = tiedRankLocal(x);
ry = tiedRankLocal(y);
cc = corrcoef(rx(:),ry(:));
if numel(cc)>=4, rho = cc(1,2); else, rho = NaN; end
n = numel(rx);
t = rho*sqrt((n-2)/max(1-rho^2,eps));
p = 2*(1 - tcdfApprox(abs(t),n-2));
end

function r = tiedRankLocal(x)
[xs,ord] = sort(x(:));
r = nan(size(xs));
i=1;
while i<=numel(xs)
    j=i;
    while j<numel(xs) && xs(j+1)==xs(i)
        j=j+1;
    end
    r(i:j) = mean(i:j);
    i=j+1;
end
rr = nan(size(r));
rr(ord)=r;
r=rr;
end

function p = tcdfApprox(t,df)
% Normal approximation to t CDF. Sufficient for plotting annotation.
if df > 30
    p = 0.5*(1+erf(t/sqrt(2)));
else
    % Simple approximation based on Hill 1970-style transform.
    a = sqrt(log(1 + t.^2/df));
    z = a .* sqrt(df - 0.5);
    p = 0.5*(1+erf(z/sqrt(2)));
end
end

function q = quantileNoToolbox(x,p)
x = sort(x(isfinite(x)));
if isempty(x), q=NaN; return; end
pos = 1 + (numel(x)-1)*p;
lo = floor(pos); hi = ceil(pos);
if lo==hi
    q=x(lo);
else
    q=x(lo)+(pos-lo)*(x(hi)-x(lo));
end
end

function r = rangeFinite(x)
x = x(isfinite(x));
if isempty(x), r=1; else, r=max(x)-min(x); if r==0, r=1; end; end
end


function structuralFile = findStructuralValidationFile(tablesDir)
% Robust search for 13_gc32_structural_validation.csv.
% The file may be in pwd, the tables directory, the script directory,
% or somewhere under MATLAB Drive. Recursive search is used as a fallback.
name = '13_gc32_structural_validation.csv';
structuralFile = '';
roots = {};
try roots{end+1} = tablesDir; catch, end
roots{end+1} = pwd;
try roots{end+1} = fileparts(mfilename('fullpath')); catch, end
try
    up = userpath;
    if ~isempty(up)
        up = split(string(up), pathsep);
        for k=1:numel(up)
            if strlength(up(k))>0, roots{end+1} = char(up(k)); end %#ok<AGROW>
        end
    end
catch
end
roots = unique(roots(~cellfun(@isempty,roots)),'stable');
for k=1:numel(roots)
    f = fullfile(roots{k},name);
    if isfile(f)
        structuralFile = f;
        return;
    end
end
for k=1:numel(roots)
    if isfolder(roots{k})
        d = dir(fullfile(roots{k},'**',name));
        if ~isempty(d)
            structuralFile = fullfile(d(1).folder,d(1).name);
            return;
        end
    end
end
end

function maybeCloseFig(fig)
% Keep figures open in interactive mode. Close only in export-only mode.
keepOpen = false;
try
    keepOpen = getappdata(0,'GC32_DISPLAY_FIGURES');
catch
end
if ~keepOpen && ishghandle(fig)
    close(fig);
end
end

% =========================================================================
% Export and style
% =========================================================================
function fig = newFig(pos)
fig = figure('Color','w','Position',pos,'Renderer','painters');
end

function setManuscriptStyle(showFigures)
if showFigures
    set(0,'DefaultFigureVisible','on');
else
    set(0,'DefaultFigureVisible','off');
end
set(groot,'DefaultAxesFontName','Helvetica');
set(groot,'DefaultTextFontName','Helvetica');
set(groot,'DefaultAxesFontSize',10);
set(groot,'DefaultTextFontSize',10);
set(groot,'DefaultAxesLineWidth',0.8);
set(groot,'DefaultLineLineWidth',1.2);
set(groot,'DefaultFigureColor','w');
end

function exportFig(fig,outDir,baseName)
if ~isfolder(outDir), mkdir(outDir); end
pdfFile = fullfile(outDir,[baseName '.pdf']);
pngFile = fullfile(outDir,[baseName '.png']);
try
    exportgraphics(fig,pdfFile,'ContentType','vector','BackgroundColor','white');
    exportgraphics(fig,pngFile,'Resolution',600,'BackgroundColor','white');
catch
    print(fig,pdfFile,'-dpdf','-painters','-bestfit');
    print(fig,pngFile,'-dpng','-r600');
end
fprintf('  saved: %s\n', baseName);
end
