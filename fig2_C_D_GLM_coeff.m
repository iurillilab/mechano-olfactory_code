clear; clc; close all;
%%

dataFilesDirList = {['M1', filesep, '102721'], ...
                    ['M2', filesep, '102821'],...
                    ['M3', filesep, '102921'],...
                    ['M5', filesep, '120221']};
mouseNameList = {'M1_102721', 'M2_102821','M3_102921','M5_120221'};
dataPathPrefix = '/Users/galileo/Dropbox (Personal)/SniffOdor_project_processedData/sniffOdorProject/Conc_Series';
%  dataPathPrefix = 'D:\sniffOdorProject\Conc_Series';
%%
constantFactor= [];
alpha = [];
beta = [];
gamma = [];
alphaPval = [];
betaPval = [];
gammaPval = [];
%%
for mn = 1 : length(mouseNameList)
    %%
    mouseName = mouseNameList{mn};
    mouseDirName = fullfile(dataPathPrefix,'processedDataStorage',...
        mouseName);
    load(fullfile(mouseDirName,...
            sprintf('%s_glmModel_with_baseline.mat', mouseName)));
   %%
   numOfneurons = size(glmModel.regModelCoeff_results_o1_shuffled, 2);
   numOfReasmpele = size(glmModel.regModelCoeff_results_o1_shuffled, 3);
   umOfShuffle = size(glmModel.regModelCoeff_results_o1_shuffled, 4);
   
   constantFactor = [constantFactor;...
       [nanmean(glmModel.regModelCoeff_results_o1(1,:,:),3)',...
       nanmean(glmModel.regModelCoeff_results_o2(1,:,:),3)']];
       
       
   alpha = [alpha;...
       [nanmean(glmModel.regModelCoeff_results_o1(2,:,:),3)',...
       nanmean(glmModel.regModelCoeff_results_o2(2,:,:),3)']];
       
   beta = [beta;...
       [nanmean(glmModel.regModelCoeff_results_o1(3,:,:),3)',...
       nanmean(glmModel.regModelCoeff_results_o2(3,:,:),3)']];
   
   gamma = [gamma;...
       [nanmean(glmModel.regModelCoeff_results_o1(4,:,:),3)',...
        nanmean(glmModel.regModelCoeff_results_o2(4,:,:),3)']];
%    mouse_alpha_p_vala = nan(2, numOfneurons, numOfReasmpele);
%    mouse_betta_p_vala = nan(2, numOfneurons, numOfReasmpele);
%    mouse_gamma_p_vala = nan(2, numOfneurons, numOfReasmpele);
%    
%     for ni = 1 : numOfneurons
%         for i = 1: numOfReasmpele
%             p_alpha_odor1 = Efron_and_Tibshirani_hypothesis_testing(...
%                 squeeze(glmModel.regModelCoeff_results_o1(2,ni,i)),...
%                 squeeze(glmModel.regModelCoeff_results_o1_shuffled(2,ni,i,:))); 
%             p_alpha_odor2 = Efron_and_Tibshirani_hypothesis_testing(...
%                 squeeze(glmModel.regModelCoeff_results_o2(2,ni,i)),...
%                 squeeze(glmModel.regModelCoeff_results_o2_shuffled(2,ni,i,:)));
%             
%             p_betta_odor1 = Efron_and_Tibshirani_hypothesis_testing(...
%                 squeeze(glmModel.regModelCoeff_results_o1(3,ni,i)),...
%                 squeeze(glmModel.regModelCoeff_results_o1_shuffled(3,ni,i,:))); 
%             p_betta_odor2 = Efron_and_Tibshirani_hypothesis_testing(...
%                 squeeze(glmModel.regModelCoeff_results_o2(3,ni,i)),...
%                 squeeze(glmModel.regModelCoeff_results_o2_shuffled(3,ni,i,:)));
%             
%             p_gamma_odor1 = Efron_and_Tibshirani_hypothesis_testing(...
%                 squeeze(glmModel.regModelCoeff_results_o1(4,ni,i)),...
%                 squeeze(glmModel.regModelCoeff_results_o1_shuffled(4,ni,i,:))); 
%             p_gamma_odor2 = Efron_and_Tibshirani_hypothesis_testing(...
%                 squeeze(glmModel.regModelCoeff_results_o2(4,ni,i)),...
%                 squeeze(glmModel.regModelCoeff_results_o2_shuffled(4,ni,i,:)));
%             
%             mouse_alpha_p_vala(1, ni, i) = p_alpha_odor1;
%             mouse_alpha_p_vala(2, ni, i) = p_alpha_odor2;
%             
%             mouse_betta_p_vala(1, ni, i) = p_betta_odor1;
%             mouse_betta_p_vala(2, ni, i) = p_betta_odor2;
%             
%             mouse_gamma_p_vala(1, ni, i) = p_gamma_odor1;
%             mouse_gamma_p_vala(2, ni, i) = p_gamma_odor2;
%         end
%     end
%    alphaPval = [alphaPval,...
%        mean(mouse_alpha_p_vala,3)];
%    betaPval = [betaPval,...
%        mean(mouse_betta_p_vala, 3)];
%    gammaPval = [gammaPval,...
%        mean(mouse_gamma_p_vala, 3)];
end
%%
pVal = .01;
% sigGamma = (gammaPval <= pVal)';
% sigAlpha = (alphaPval <= pVal)';
% sigBeta = (betaPval <= pVal)';
%%
xLim = [-2,2];
yLim = xLim;
xLimRainCloud = [-.3, 1.5];
markerSize = 25;
transparancyAlfa = .5;
dotFaceColor = .4*[1,1,1];
dotEdgeColor = .1*[1,1,1];
%
% [colorSet, ~] = set_plot_seting(15, 8);
pltDir = fullfile(final_figs_path('paola'),...
    'fig4_C_D_GLM_with_baseline_2_wind_notABS', 'single-odors');
if ~exist(pltDir , 'dir')
   mkdir(pltDir)
end
%%
for oId = 1 : 2
   close all; figure, hold on
   plot([0,0], xLim, 'k:')
   plot( xLim, [0,0], 'k:')
   mdl = fitlm(alpha(:, oId),...
        beta(:, oId));
   scatter(alpha(:, oId),...
        beta(:, oId),...
        markerSize,...
        'MarkerEdgeColor', dotEdgeColor,...
        'MarkerFaceColor', dotFaceColor, 'MarkerFaceAlpha', transparancyAlfa);
  
    plot_lr_lineAnd_inerval(mdl, xLim')
    R = corrcoef(alpha(:, oId),...
        beta(:, oId)); 
    ylim(yLim)
    xlim(xLim)
    xlabel('inhalation type regressor')
    ylabel('concentration regressor')    
    title({sprintf('Odor%d - GLM model',...
        oId),...
        sprintf('Corr. coeff.= %.3f',...
        R(1,2)), ''})
    grid on
    axis square
    set(gcf,'Position',[200,200,400,600]);
    print_it(pltDir, sprintf('od%d_alphaVsBetta', oId),...
        'pooled')
    %%
    figure, hold on
    plot([0,0], xLim, 'k--')
    plot( xLim, [0,0], 'k--')
    mdl = fitlm(beta(:, oId),...
        gamma(:, oId));
    scatter(beta(:, oId),...
        gamma(:, oId),...
        markerSize,...
        'MarkerEdgeColor', dotEdgeColor,...
        'MarkerFaceColor', dotFaceColor, 'MarkerFaceAlpha', transparancyAlfa);
    plot_lr_lineAnd_inerval(mdl, xLim')
    ylim(yLim)
    xlim(xLim)
    R = corrcoef(beta(:, oId),...
        gamma(:, oId)); 
    ylabel('interaction regressor')
    xlabel('concentration regressor')

    title({sprintf('Odor%d - GLM model ',...
        oId),...
        sprintf('Corr. coeff.= %.3f',...
        R(1,2)), ''})
    grid on
    axis square
    set(gcf,'Position',[200,200,400,600]);
   print_it(pltDir, sprintf('od%d_BetaVsGamma', oId),...
        'pooled')
    %%
    figure, hold on
    plot([0,0], xLim, 'k--')
    plot( xLim, [0,0], 'k--')
    mdl = fitlm(alpha(:, oId),...
        gamma(:, oId));
    scatter(alpha(:, oId),...
        gamma(:, oId),...
        markerSize,...
        'MarkerEdgeColor', dotEdgeColor,...
        'MarkerFaceColor', dotFaceColor, 'MarkerFaceAlpha', transparancyAlfa);
    plot_lr_lineAnd_inerval(mdl, xLim')
    ylim(yLim)
    xlim(xLim)
    R = corrcoef(alpha(:, oId),...
        gamma(:, oId)); 
    ylabel('interaction regressor')
    xlabel('inhalation type regressor')

    title({sprintf('Odor%d - GLM model ',...
        oId),...
        sprintf('Corr. coeff.= %.3f',...
        R(1,2)), ''})
    grid on
    axis square
    set(gcf,'Position',[200,200,400,600]);
    print_it(pltDir, sprintf('od%d_alohaVsGamma', oId),...
        'pooled')
    %%
    close all; figure;
    raincloud_plot((alpha(:, oId)),...
        'color',[.5,.5,.5], 'box_on', 1, 'alpha', .5,...
        'lwr_bnd', 1, 'line_width', .5, 'band_width', .005,...
        'dot_size', 3);
    xlim(xLimRainClouth)
    xlabel('inhalation type regressor')
    title({sprintf('Odor%d - GLM model Alpha',...
        oId)})
    box off
    set(gcf,'Position',[350,350,550,300]);
    print_it(pltDir, sprintf('od%d_alpha', oId),...
         'pooled')
    %%
    close all; figure;
    raincloud_plot((beta(:, oId)),...
        'color',[.5,.5,.5], 'box_on', 1, 'alpha', .5,...
        'lwr_bnd', 1, 'line_width',.5, 'band_width', .005,...
        'dot_size', 3);
    xlim(xLimRainClouth)
    xlabel('concentration regressor')
    title({sprintf('Odor%d - GLM model Beta',...
        oId)})
    box off
    set(gcf,'Position',[350,350,550,300]);
    print_it(pltDir, sprintf('od%d_beta', oId),...
         'pooled')
    %%
    close all; figure;
    raincloud_plot((gamma(:, oId)),...
        'color',[.5,.5,.5], 'box_on', 1, 'alpha', .5,...
        'lwr_bnd', 1, 'line_width', .5, 'band_width', .005,...
        'dot_size', 3);
    xlim(xLimRainClouth)
    xlabel('intraction regressor')
    title({sprintf('Odor%d - GLM model Gamma',...
        oId)})
    box off
    set(gcf,'Position',[350,350,550,300]);
    print_it(pltDir, sprintf('od%d_gamma', oId),...
         'pooled')
    %%
end
%%
% [colorSet, ~] = set_plot_seting(15, 8);
pltDir = fullfile(final_figs_path('paola'),...
    'fig4_C_D_GLM_with_baseline_2_wind_notABS', 'pooled-odors');
if ~exist(pltDir , 'dir')
   mkdir(pltDir)
end
%%
close all; figure, hold on
plot([0,0], xLim, 'k:')
plot( xLim, [0,0], 'k:')
mdl = fitlm(alpha(:),...
    beta(:));
scatter(alpha(:),...
    beta(:),...
    markerSize,...
    'MarkerEdgeColor', dotEdgeColor,...
    'MarkerFaceColor', dotFaceColor,...
    'MarkerFaceAlpha', transparancyAlfa);

plot_lr_lineAnd_inerval(mdl, xLim')
R = corrcoef(alpha(:),...
    beta(:)); 
ylim(yLim)
xlim(xLim)
xlabel('inhalation type regressor')
ylabel('concentration regressor')    
title({sprintf('pooled odors - GLM model'),...
    sprintf('Corr. coeff.= %.3f',...
    R(1,2)), ''})
grid on
axis square
set(gcf,'Position',[200,200,400,600]);
print_it(pltDir, sprintf('pooled_odors_alphaVsBetta'),...
    'pooled')
%%
figure, hold on
plot([0,0], xLim, 'k--')
plot( xLim, [0,0], 'k--')
% mdl = fitlm(beta(:),...
%     gamma(:));
scatter(beta(:),...
    gamma(:),...
    markerSize,...
    'MarkerEdgeColor', dotEdgeColor,...
    'MarkerFaceColor', dotFaceColor,...
    'MarkerFaceAlpha', transparancyAlfa);
%  plot_lr_lineAnd_inerval(mdl, xLim')
ylim(yLim)
xlim(xLim)
R = corrcoef(beta(:),...
    gamma(:)); 

ylabel('interaction regressor')
xlabel('concentration regressor')

title({sprintf('pooled odors - GLM model'),...
    sprintf('Corr. coeff.= %.3f',...
    R(1,2)), ''})
grid on
axis square
set(gcf,'Position',[200,200,400,600]);
print_it(pltDir, sprintf('pooled_odors_BetaVsGamma'),...
    'pooled')
%%
figure, hold on
plot([0,0], xLim, 'k--')
plot( xLim, [0,0], 'k--')
% mdl = fitlm(alpha(:, :),...
%         gamma(:, :));
scatter(alpha(:),...
    gamma(:),...
    markerSize,...
    'MarkerEdgeColor', dotEdgeColor,...
    'MarkerFaceColor', dotFaceColor,...
    'MarkerFaceAlpha', transparancyAlfa);
% plot_lr_lineAnd_inerval(mdl, xLim')
ylim(yLim)
xlim(xLim)
R = corrcoef(alpha(:),...
    gamma(:)); 
ylabel('interaction regressor')
xlabel('inhalation type regressor')

title({sprintf('pooled odors - GLM model '),...
    sprintf('Corr. coeff.= %.3f',...
    R(1,2)), ''})
grid on
axis square
set(gcf,'Position',[200,200,400,600]);
print_it(pltDir, sprintf('pooled_odors_alphaVsGamma'),...
    'pooled')
%%
DOTSIZE = 10; %doesnt do anything
xLimRainCloud = [-2, 2];
close all; figure;
raincloud_plot((alpha(:)),...
    'color',[.5,.5,.5], 'box_on', 1, 'alpha', .5,...
    'lwr_bnd', 1, 'line_width', .5, 'band_width', .005,...
    'dot_size', DOTSIZE);
xlim(xLimRainCloud)
xlabel('inhalation type regressor')
title({sprintf('pooled odors - GLM model Alpha')})
box off
set(gcf,'Position',[350,350,550,300]);
print_it(pltDir, sprintf('pooled_odors_alpha'),...
     'pooled')
%%
close all; figure;
raincloud_plot((beta(:)),...
    'color',[.5,.5,.5], 'box_on', 1, 'alpha', .5,...
    'lwr_bnd', 1, 'line_width', .5, 'band_width', .005,...
    'dot_size', DOTSIZE);
xlim(xLimRainCloud)
xlabel('concentration regressor')
title({sprintf('pooled odors - GLM model Beta')})
box off
set(gcf,'Position',[350,350,550,300]);
print_it(pltDir, sprintf('pooled_odors_beta'),...
     'pooled')
%%
close all; figure;
raincloud_plot((gamma(:)),...
    'color',[.5,.5,.5], 'box_on', 1, 'alpha', .5,...
    'lwr_bnd', 1, 'line_width', .5, 'band_width', .005,...
    'dot_size', DOTSIZE);
xlim(xLimRainCloud)
xlabel('intraction regressor')
title({sprintf('pooled odors - GLM model Gamma')})
box off
set(gcf,'Position',[350,350,550,300]);
print_it(pltDir, sprintf('pooled_odors_gamma'),...
     'pooled')

 
 
 
 
 
%% scattered data boxplot (not absolute values)

boxeddata = cat(1, alpha(:), beta(:), gamma(:));
boxedlabels = cat(1, ones(size(alpha(:))), 2*ones(size(beta(:))), 3*ones(size(gamma(:))));
figure; hold on

% add single data points
figure; hold on
for ei = 1:3
    X = boxeddata(boxedlabels==ei);
    X(X>2) = 2;
    Y = ei + (rand(sum(boxedlabels==ei),1)-0.5)*0.4;
    scatter(X,Y, 100, 0.5*[1 1 1], 'filled')   
end
xlim([-2 2])
box off
xticks(-2:1:2)

print_it(pltDir, sprintf('pooled_odors_ALL'),'pooled')


%% boxplot only (absolute values)
figure; hold on
boxplot(abs(boxeddata), boxedlabels, 'Colors', 'k', 'Jitter',1, 'Symbol','', 'Orientation', 'horizontal', 'Whisker', 3)
ax = gca;
ax.YDir = 'reverse';
ax.YAxis.Visible = 'off';
ax.XTickLabel = [];
ax.XAxis.LineWidth = 0.3;
axx = {ax.Children.Children};
for i = 16:21
    axx{1}(i).LineStyle = '-';
end
% whiskers are adapted to the 95% lenght post hoc (vectorially) with the
% correct length reported in the figure:
text(0.6,1, sprintf('%1.3f',prctile(abs(reshape(alpha(rs,:),[],1)),95)))
text(0.6,2, sprintf('%1.3f',prctile(abs(reshape(beta(rs,:),[],1)),95)))
text(0.6,3, sprintf('%1.3f',prctile(abs(reshape(gamma(rs,:),[],1)),95)))



%% separate in fast spiking and regular spiking for supplementary figure:
load('/Users/galileo/Dropbox (Personal)/SniffOdor_project_processedData/sniffOdorProject/Bulbectomy/processedDataStorage/odorResp_RSFS_mechano_tags_piriform.mat')
fs = fastspikelabel(1:464)==1;
rs = fastspikelabel(1:464)==2;

%% boxplot only (absolute values)
figure; hold on

%% rs:
subplot(3,1,1)
boxeddata = cat(2, alpha(rs,:), beta(rs,:), gamma(rs,:));
boxeddata = reshape(boxeddata, [], 3);
boxedlabels = cat(1, ones(numel(alpha(rs,:)),1), 2*ones(numel(alpha(rs,:)),1), 3*ones(numel(alpha(rs,:)),1));
boxplot(abs(boxeddata), boxedlabels, 'Colors', 'k', 'Jitter',1, 'Symbol','', 'Orientation', 'horizontal', 'Whisker', 3)
ax = gca;
ax.YDir = 'reverse';
ax.YAxis.Visible = 'off';
ax.XTickLabel = [];
ax.XAxis.LineWidth = 0.3;
axx = {ax.Children.Children};
for i = 16:21
    axx{1}(i).LineStyle = '-';
end
xlim([0 1])
box off
xticks(0:0.1:1)
%
text(0.6,1, sprintf('%1.3f',prctile(abs(reshape(alpha(rs,:),[],1)),95)))
text(0.6,2, sprintf('%1.3f',prctile(abs(reshape(beta(rs,:),[],1)),95)))
text(0.6,3, sprintf('%1.3f',prctile(abs(reshape(gamma(rs,:),[],1)),95)))

%% fs:
subplot(3,1,2)
boxeddata = cat(2, alpha(fs,:), beta(fs,:), gamma(fs,:));
boxeddata = reshape(boxeddata, [], 3);
boxedlabels = cat(1, ones(numel(alpha(fs,:)),1), 2*ones(numel(alpha(fs,:)),1), 3*ones(numel(alpha(fs,:)),1));
boxplot(abs(boxeddata), boxedlabels, 'Colors', 'k', 'Jitter',1, 'Symbol','', 'Orientation', 'horizontal', 'Whisker', 3)
ax = gca;
ax.YDir = 'reverse';
ax.YAxis.Visible = 'off';
ax.XTickLabel = [];
ax.XAxis.LineWidth = 0.3;
axx = {ax.Children.Children};
for i = 16:21
    axx{1}(i).LineStyle = '-';
end
xlim([0 1])
box off
xticks(0:0.1:1)
%
text(0.6,1, sprintf('%1.3f',prctile(abs(reshape(alpha(fs,:),[],1)),95)))
text(0.6,2, sprintf('%1.3f',prctile(abs(reshape(beta(fs,:),[],1)),95)))
text(0.6,3, sprintf('%1.3f',prctile(abs(reshape(gamma(fs,:),[],1)),95)))

savefig(fullfile(pltDir, 'temp_subplots_1_2_RsFs.fig'))

%% add 70ms in subplot(3)
% check dedicated script

