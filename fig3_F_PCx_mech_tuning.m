clear; clc; close all;
[dataFilesDirList, mouseNameList, dataPathPrefixList] =....
    get_all_healthy_mice('uniTn1');
add_ndt_paths_and_init_rand_generator
%%
rasterMapPooled = [];
rasterMap_sem_Pooled = [];
slopeVecPooled = [];
slopePvalVecPooled = [];
unitsIdsPooled = [];
mouseIdPooled = [];
numOfQuantileSplit = 9;

    
for mn = 1  : length(mouseNameList)    
    %% loading data 
    mouseName = mouseNameList{mn};
    fprintf('Analysis %s data... \n', mouseName);
    mouseNameToPrint = mouseName(4:end);
    dataPathPrefix = dataPathPrefixList{mn};
    mouseDirName = fullfile(dataPathPrefix,'processedDataStorage',...
        mouseName);
    load(fullfile(mouseDirName,...
        sprintf('\\%s_respData.mat', mouseName)));
    load(fullfile(mouseDirName,...
        sprintf('\\%s_respLabelsData.mat', mouseName)));
    load(fullfile(mouseDirName,...
        sprintf('\\%s_inh_NR.mat', mouseName)));
    load(fullfile(mouseDirName,...
        sprintf('\\%s_recordingData.mat', mouseName)));
    load(fullfile(mouseDirName,...
        sprintf('\\%s_respTrace.mat', mouseName)));
    %% get OCtc units
    goodAndSafe_unitsFlag =...
        get_goodAndSafe_units(dataPathPrefix, mouseName,...
        6, .5);
    goodAndSafe_unitsId = find(goodAndSafe_unitsFlag);
    %% cuting airflow trace to the quantiles
    baselineRespMaxAirFlow =...
        abs(resp(respData.xPeaks(respData.xpeakCondLabel == 0)));
    inhQuantiles = quantile(baselineRespMaxAirFlow, numOfQuantileSplit);
    inhQuantileBin = nan(size(baselineRespMaxAirFlow));
    for qIdx = 1: length(inhQuantiles)
        if qIdx == 1
            inhQuantileBin(baselineRespMaxAirFlow < inhQuantiles(qIdx)) =...
                qIdx;
        elseif qIdx == length(inhQuantiles)
            inhQuantileBin(...
                (baselineRespMaxAirFlow >= inhQuantiles(qIdx-1))&...
                (baselineRespMaxAirFlow < inhQuantiles(qIdx))) = qIdx;
            inhQuantileBin(baselineRespMaxAirFlow >= inhQuantiles(qIdx)) =...
                qIdx+1;
        else
            inhQuantileBin(...
                (baselineRespMaxAirFlow >= inhQuantiles(qIdx-1))&...
                (baselineRespMaxAirFlow < inhQuantiles(qIdx))) = qIdx;
        end       
    end 
    %%
    avg_fr = inh_NR.after_sc/...
        (inh_NR.cacheInfo.windSizeSec);
    baseline_inhs_avg_fr =...
        avg_fr(respData.xpeakCondLabel == 0, :);
    baseline_avg_fr_normalized =...
        zscore(baseline_inhs_avg_fr);
%     baseline_avg_fr_normalized(:,...
%         mean(baseline_inhs_avg_fr) == 0) = 0;
    
   %% geting tuning slope and geting rastermap
   rasterMap = nan(sum(goodAndSafe_unitsFlag),...
       max(inhQuantileBin)); 
   rasterMap_sem = nan(sum(goodAndSafe_unitsFlag),...
       max(inhQuantileBin));
   
   slopeVec = nan(sum(goodAndSafe_unitsFlag), 1);
   slopePvalVec = nan(sum(goodAndSafe_unitsFlag), 1);
   for nn = 1 : length(goodAndSafe_unitsId)
       ni = goodAndSafe_unitsId(nn);
       %% caculate the slope    
       neuronResp = baseline_avg_fr_normalized(:, ni);
       
       mdl_fited = fitlm(inhQuantileBin,...
           neuronResp);
       slopeVec(nn) = mdl_fited.Coefficients.Estimate(2);
       slopePvalVec(nn) = mdl_fited.Coefficients.pValue(2);
       %%--
       for bi = 1 : max(inhQuantileBin)
           rasterMap(nn, bi) =...
               mean(neuronResp(...
               inhQuantileBin== bi));
           rasterMap_sem(nn, bi) =...
               std(neuronResp(...
               inhQuantileBin == bi))/...
               sqrt(sum(inhQuantileBin== bi));
       end
   end
   %% 
    rasterMapPooled = [rasterMapPooled;...
        rasterMap];
    rasterMap_sem_Pooled = [rasterMap_sem_Pooled;...
        rasterMap_sem];
    slopeVecPooled = [slopeVecPooled;...
        slopeVec];
    slopePvalVecPooled = [slopePvalVecPooled;...
        slopePvalVec];
    unitsIdsPooled = [unitsIdsPooled; goodAndSafe_unitsId];
    mouseIdPooled = [mouseIdPooled;...
        ones(size(goodAndSafe_unitsId))*mn];
end
%% sorting slope
pval_cutOff = .01;
[sortedSlope, sortedSlopeIdx] = sort(slopeVecPooled);
rasterMap_sorted = rasterMapPooled(sortedSlopeIdx, :);
rasterMap_sem_sorted = rasterMap_sem_Pooled(sortedSlopeIdx, :);
slopePvalVecPooledSorted = slopePvalVecPooled(sortedSlopeIdx);
sigSlopeVecSorted = slopePvalVecPooledSorted <= pval_cutOff;
unitsIdsPooled_sotred = unitsIdsPooled(sortedSlopeIdx);
mouseIdPooled_sotred = mouseIdPooled(sortedSlopeIdx);

rasterMap_sorted_significant =....
    rasterMap_sorted(sigSlopeVecSorted, :);
rasterMap_sem_sorted_significant =....
    rasterMap_sem_sorted(sigSlopeVecSorted, :);
unitsIdsPooled_sotred_sigOnes = unitsIdsPooled_sotred(sigSlopeVecSorted);
mouseIdPooled_sotred_sigOnes = mouseIdPooled_sotred(sigSlopeVecSorted);
%% Spearman distrubution
[spearman_RHO,spearman_Pval] =...
    corr(rasterMap_sorted',(1:numOfQuantileSplit+1)',...
    'Type','Spearman');
sigSpearmanCoeff = spearman_Pval <= pval_cutOff;
%%
colorSet = set_plot_seting(15);
plotDir = fullfile(final_figs_path('uniTn1'),...
    'fig_F_H');
if ~exist(plotDir , 'dir')
   mkdir(plotDir)
end
%% Spearman dist plot
close all;
figure;
hold on
histogram((spearman_RHO), -1:.05:1,...
    'EdgeAlpha', 0, 'FaceColor', [.3,.3,.3]);
histogram((spearman_RHO(sigSpearmanCoeff)), -1:.05:1,...
    'EdgeAlpha', 0, 'FaceColor', [0,.5,.1]);
ylabel('frequency')
xlabel('Spearmann coefficient')
ylim([0, 95])
legend({'not significant', 'significant'}, 'Location', 'bestoutside')
legend boxoff
box off
set(gcf,'Position',[100 100 700 250]);
print_it(plotDir, 'tuning_curves_SpearmanCoeff_distribution', 'pooled')
%% raster map
num_of_bins = numOfQuantileSplit+1;
xTiksLabels =...
    round((1/num_of_bins : 1/num_of_bins:1)-...
    (1/(2*num_of_bins)), 2);
xtikesSet = 1:2:num_of_bins;
cRange = [-.8 1.5];
close all; figure;
imagesc(rasterMap_sorted_significant)
colorbar
box off
xtickangle(45)
caxis(cRange)
xticks(xtikesSet)
xticklabels(xTiksLabels(xtikesSet));
xtickangle(45)
yticks([])
xlabel('inh. max airflow quantile')
ylabel('Unins')
set(gcf,'Position',[100 100 350 800]);
colormap(bluewhitered)
print_it(plotDir, 'rastermap_maxAirFlow_tuning', 'pooled')
%%
plot_bin_size = .009;
close all;
figure;
histogram(sortedSlope, -.25:plot_bin_size:.25,...
    'EdgeAlpha', 0, 'FaceColor', [.3,.3,.3]);
hold on
histogram(sortedSlope(sigSlopeVecSorted),...
    -.25:plot_bin_size:.25,...
    'EdgeAlpha', 0, 'FaceColor', [0,.5,.1]);
% legend({'', 'significant slops'})
ylabel('frequency')
xlabel('linear regression slope')
ylim([0, 95])
box off
legend({'not significant', 'significant'}, 'Location', 'bestoutside')
legend boxoff
set(gcf,'Position',[100 100 700 250]);
print_it(plotDir, 'lin_reg_slope_distribution', 'pooled')
%%
fileID = fopen(fullfile(plotDir,...
            'linear_regression_slope_data.txt'),'w');
fprintf(fileID, '%d of %d niuts have significant lin. reg. slpoe. (p-val cutoff = %.5f) \n',...
    sum(sigSlopeVecSorted), length(sortedSlope), pval_cutOff);
fprintf(fileID, 'Number of units with a significant positive lin. reg. slope = %d. \n',...
    sum(sortedSlope(sigSlopeVecSorted)>0));
fprintf(fileID, 'Number of units with a significant negetive lin. reg. slope = %d. \n',...
    sum(sortedSlope(sigSlopeVecSorted)<0));
fclose(fileID);

fileID = fopen(fullfile(plotDir,...
            'spearman_coeff_data.txt'),'w');
fprintf(fileID, '%d of %d niuts have significant spearman coeff. (p-val cutoff = %.5f) \n',...
    sum(sigSpearmanCoeff), length(sigSpearmanCoeff), pval_cutOff);
fprintf(fileID, 'ads(spearman coeff): mean +- SEM: %.5f +- %.5f. \n',...
    mean(abs(spearman_RHO)),...
    std(spearman_RHO)/sqrt(length(spearman_RHO)) );
fclose(fileID);
%%
% brforeInhOnsetWind = .6;
% afterInhOnsetWind = .6;
% PITH_time = -1*brforeInhOnsetWind:...
%             .001 : afterInhOnsetWind +.001; 
% smoothingCurnel = 10;
% ylimRange = [-13,40];
% xlimRange = [-.3, .4];
% rastersDotSize = 2;
% lineWidth = 2.5;
% colors = cbrewer2('seq', 'Blues', numOfQuantileSplit+1);
% 
% c = 1;
% for i = [743, 739, 697, 675, 669, 638, 602]
%     %%
%     plotDirExampels =fullfile(plotDir, 'examples_of_F',...
%         sprintf('example%d', c));
%     if ~exist(plotDirExampels , 'dir')
%        mkdir(plotDirExampels)
%     end
%     c = c+1;
%     %%
%     mn = mouseIdPooled_sotred_sigOnes(i);
%     ni = unitsIdsPooled_sotred_sigOnes(i);
%     %%
%     mouseName = mouseNameList{mn};
%     dataPathPrefix = dataPathPrefixList{mn};
%     mouseDirName = fullfile(dataPathPrefix,'processedDataStorage',...
%         mouseName);
%     respRasterDir = fullfile(mouseDirName, 'RstersFiles\inhRasters');
%     load(fullfile(mouseDirName,...
%         sprintf('\\%s_respData.mat', mouseName)));
%     load(fullfile(mouseDirName,...
%         sprintf('\\%s_respLabelsData.mat', mouseName)));
%     load(fullfile(mouseDirName,...
%         sprintf('\\%s_inh_NR.mat', mouseName)));
%     load(fullfile(mouseDirName,...
%         sprintf('\\%s_unitsBaselineFR.mat', mouseName)));
%     load(fullfile(mouseDirName,...
%         sprintf('\\%s_recordingData.mat', mouseName)));
%     load(fullfile(mouseDirName,...
%             sprintf('\\%s_respTrace.mat', mouseName)));
%     %%----
%     unitRasterData = load([respRasterDir,...
%             sprintf('\\%s_NPX_unit%d_resps_raster_data.mat',...
%             mouseName, ni)]);
%     %%---
%     neuronInhSpikeRaster_baseline = unitRasterData.raster_data(...
%         unitRasterData.raster_labels.xpeakCondLabel == 0, :);
% 
%     baselineRespMaxAirFlow =...
%         abs(resp(respData.xPeaks(respData.xpeakCondLabel == 0)));
% 
%     inhQuantiles = quantile(baselineRespMaxAirFlow, numOfQuantileSplit);
%     inhQuantileBin = nan(size(baselineRespMaxAirFlow));
%     for qIdx = 1: length(inhQuantiles)
%         if qIdx == 1
%             inhQuantileBin(baselineRespMaxAirFlow < inhQuantiles(qIdx)) =...
%                 qIdx;
%         elseif qIdx == length(inhQuantiles)
%             inhQuantileBin(...
%                 (baselineRespMaxAirFlow >= inhQuantiles(qIdx-1))&...
%                 (baselineRespMaxAirFlow < inhQuantiles(qIdx))) = qIdx;
%             inhQuantileBin(baselineRespMaxAirFlow >= inhQuantiles(qIdx)) =...
%                 qIdx+1;
%         else
%             inhQuantileBin(...
%                 (baselineRespMaxAirFlow >= inhQuantiles(qIdx-1))&...
%                 (baselineRespMaxAirFlow < inhQuantiles(qIdx))) = qIdx;
%         end       
%     end 
%     %%
%     neuon_mech_unning_curve =...
%         nan(max(inhQuantileBin), 1);
%     neuon_mech_unning_curve_sem =...
%         nan(max(inhQuantileBin), 1);
%     
%     avg_fr = inh_NR.after_sc/...
%         (inh_NR.cacheInfo.windSizeSec);
%     baseline_inhs_avg_fr =...
%         avg_fr(respData.xpeakCondLabel == 0, :);
%     
%     neuronResp = baseline_inhs_avg_fr(:, ni);
%     %%--
%     for bi = 1 : max(inhQuantileBin)
%        neuon_mech_unning_curve(bi) =...
%            mean(neuronResp(...
%            inhQuantileBin== bi));
%        neuon_mech_unning_curve_sem( bi) =...
%            std(neuronResp(...
%            inhQuantileBin == bi))/...
%            sqrt(sum(inhQuantileBin== bi));
%     end
%     neuon_mech_unning_curve =...
%         neuon_mech_unning_curve - unitsBaselineFR(ni);
%     %%
%     close all
%     errorbar(1:10, neuon_mech_unning_curve,...
%         neuon_mech_unning_curve_sem,...
%         'Marker', 'o', 'Color', colorSet(2,:),...
%         'CapSize', 3, 'LineWidth', lineWidth, 'MarkerSize', 4)
%     box off
%     xtickangle(45)
%     caxis(cRange)
%     xticks(xtikesSet)
%     xticklabels(xTiksLabels(xtikesSet));
%     xtickangle(45)
%     xlabel('inh. max airflow quantile')
%     ylabel('\Delta firing rate (sp./s)')
%     title(sprintf('Id = %d', i))
%     set(gcf,'Position',[100 100 500 300]);
%     print_it(plotDirExampels,'tuning_curve', num2str(i))
%     %%
%     neuronBinRasterCell = cell(numOfQuantileSplit+1, 1);
%     neuronsBinPITHCell = cell(numOfQuantileSplit+1, 1);
% 
%     for bi = 1 : max(inhQuantileBin)
%        neuronBinRaster = neuronInhSpikeRaster_baseline(...
%             inhQuantileBin == bi, :);
%        neuronsRtype_PITH =...
%            gen_fx_gsmooth((mean(neuronBinRaster,1)*1000)-...
%            unitsBaselineFR(ni),...
%            smoothingCurnel);
%        neuronBinRasterCell{bi} = neuronBinRaster;
%        neuronsBinPITHCell{bi} = neuronsRtype_PITH;
%     end
% 
%     %% plot PITH
%     close all
%     figure;
%     hold on
%     for bi = 1 :  max(inhQuantileBin)
%         plot(PITH_time,...
%             neuronsBinPITHCell{bi},'Color', colors(bi,:),...
%             'LineWidth', lineWidth)
%     end
%     cb = colorbar('FontSize', 12 , 'Box',...
%     'off', 'Location', 'eastoutside');
%     colormap(colors)
%     p  = patchline([0,0],ylim,'linestyle','--','edgecolor','k',...
%     'linewidth',2,'edgealpha',0.7);
%     xlim(xlimRange)
%     xticks([-.2,0,.2,.4])
%     xlabel('time(sec)');
%     ylabel('\Delta firing rate (sp./s)')
%     box off
%     set(gcf,'Position',[100 100 400 400]);
%     print_it(plotDirExampels,sprintf('PITH_id%d_ni%d', i, ni), mouseName)
%     
%         %%
% end
% 
% % plot(rasterMap_sorted_significant(760,:))