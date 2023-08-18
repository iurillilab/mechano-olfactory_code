clear; clc; close all;
savingDir =...
    'D:\Reza\sniffOdorProject\AI_PCtx';
% savingDir =...
%     'G:\My Drive\sinff-odor-project\raw_data';
mouseNameList = { 'M3_221117'...
    'M4_221118', 'M5_221119'};
%%
pooled_rasterMaps = []; 
pooled_slope = [];
pooled_pVal = [];
pooled_OCtx_single_good_Ids = [];
%%
stimuliDuration = .15;
frThreshold = .1;
pVal_cuttOff = 0.01;
beforeOnsetWind = .5;
afterOnsetWind = .5;
smoothingKurnel = 10;
num_of_vac_inh = 200;
afRAnge = [50, 100, 150, 200, 250];
cRange = [-.8, 1.5];
%%
[colorSet, SandR_colors] = set_plot_seting(15, 4);
pltDir = fullfile(final_figs_path('uniTn1'),...
    'fig5_J_K');
if ~exist(pltDir , 'dir')
   mkdir(pltDir)
end
%%
for mn = 1 : length(mouseNameList)
    mouseName = mouseNameList{mn};
    %% load data 
    load(fullfile(savingDir, mouseName,...
        sprintf('%s_eventsInfo.mat', mouseName)));
    load(fullfile(savingDir, mouseName,...
        sprintf('%s_spikeData.mat', mouseName)));
    load(fullfile(savingDir, mouseName,...
        sprintf('%s_expParams.mat', mouseName)));
    load(fullfile(savingDir, mouseName,...
        sprintf('%s_anatomyReconstructionTable.mat', mouseName)));
   %% get single units spikes time
    [singleUnits, singleUnitsIds, singleUnitsFr] =...
        get_singleUnits(spikeStruct);
    %% define container

    tuningTestResultMean = [];
    slopeVecMean = [];
    slopePvalVecMean = [];
    rasterMapsMean =...
        nan(length(singleUnits),...
        length(afRAnge));

    unitAirflowTuningMatMean_mean = nan(length(singleUnits),...
        length(afRAnge));
    unitAirflowTuningMatMean_sem = nan(length(singleUnits),...
        length(afRAnge));

    %% set af_vec
    af_vec = [];
    af_idx = 1;
    for af = afRAnge
       af_vec = [af_vec;...
            repmat(af_idx, [sum(eventsInfo(:, 3) == af), 1])];
       af_idx = af_idx +1;
    end
    %% get neurons with tuning and raster map and PSTH mat
    
    for sui = 1 : length(singleUnits)
        unitST = singleUnits{sui}.unitSt;
        [raster_mat, unitBaseLine_fr, evoke_raster_mat] =...
            get_raster_mat(unitST, eventsInfo(:, 1), beforeOnsetWind,...
            afterOnsetWind, beforeOnsetWind);
        %%
        evokeRespToAirflowMat = nan(num_of_vac_inh,...
            length(afRAnge));
        evokeRespToAirflowVec = [];
        afi = 1; 
        for af = afRAnge
            if sum(eventsInfo(:, 3) == af) > 0
                af_trials_deltafr = trials_fr_from_raster(...
                    raster_mat(eventsInfo(:, 3) == af, :),...
                    round(beforeOnsetWind*1000) + 1,...
                    round((beforeOnsetWind+afterOnsetWind)*1000))-...
                    trials_fr_from_raster(...
                    raster_mat(eventsInfo(:, 3) == af, :),...
                    1, round(beforeOnsetWind*1000));
                af_trials_baseline = trials_fr_from_raster(...
                    raster_mat(eventsInfo(:, 3) == af, :),...
                    1, round(beforeOnsetWind*1000));
                evokeRespToAirflowMat(:, afi)=...
                    af_trials_deltafr;
                evokeRespToAirflowVec = [evokeRespToAirflowVec;...
                    af_trials_deltafr];
            end
            afi = afi+1;
        end   
        %--- mean responce
        unitAirflowTuningMatMean_mean(sui, :) =...
            nanmean(evokeRespToAirflowMat);
        unitAirflowTuningMatMean_sem(sui, :) =...
            nanstd(evokeRespToAirflowMat)/...
            sqrt(size(evokeRespToAirflowMat, 1));
     
        rasterMapsMean(sui, :) =...
            (nanmean(evokeRespToAirflowMat)-...
             nanmean(evokeRespToAirflowMat, [1,2]))/...
             nanstd(evokeRespToAirflowMat, [], [1,2]);
        
        mdl_fited = fitlm(af_vec, zscore(evokeRespToAirflowVec));
        slopeVecMean(sui) = mdl_fited.Coefficients.Estimate(2);
        slopePvalVecMean(sui) = mdl_fited.Coefficients.pValue(2);       
    end
    %% sorting olfactory units
    margingToRemove = 50;
    OCtc_unitsFlag = false(size(unitsAreaInfoTable.K));
    Olf_stractures = {'OLF', 'EPd', 'PIR'};
    for i = 1: length(OCtc_unitsFlag)
         OCtx_flag = ismember(unitsAreaInfoTable.acronym{i},...
            Olf_stractures);
         if isempty(unitsAreaInfoTable.upperRegion{i})
             OCtx_flag = OCtx_flag &....
                ((abs(unitsAreaInfoTable.estimatedDepth(i) -...
                unitsAreaInfoTable.upperBorder(i)) > 50));
         else
             OCtx_flag = OCtx_flag &....
                ((abs(unitsAreaInfoTable.estimatedDepth(i) -...
                unitsAreaInfoTable.upperBorder(i)) > 50) |...
                ismember(unitsAreaInfoTable.upperRegion{i},...
                Olf_stractures) );
         end

         if isempty(unitsAreaInfoTable.lowerRegion{i})
             OCtx_flag = OCtx_flag &....
                ((abs(unitsAreaInfoTable.estimatedDepth(i) -...
                unitsAreaInfoTable.lowerBorder(i)) > 50));
         else
             OCtx_flag = OCtx_flag &....
                ((abs(unitsAreaInfoTable.estimatedDepth(i) -...
                unitsAreaInfoTable.lowerBorder(i)) > 50) |...
                ismember(unitsAreaInfoTable.lowerRegion{i},...
                Olf_stractures) );
         end


         OCtc_unitsFlag(i) = OCtx_flag;

    end
    OCtc_unitsId = unitsAreaInfoTable.K(...
        OCtc_unitsFlag);
    %%
    all_Octx_UnitsId = find((singleUnitsFr > frThreshold) &...
        ismember(singleUnitsIds, OCtc_unitsId));
    %%
    pooled_rasterMaps = [pooled_rasterMaps;...
        rasterMapsMean(...
        all_Octx_UnitsId, :)]; 
    pooled_slope = [pooled_slope;...
        slopeVecMean(all_Octx_UnitsId)'];
    pooled_pVal = [pooled_pVal;...
        slopePvalVecMean(all_Octx_UnitsId)'];
    pooled_OCtx_single_good_Ids = [pooled_OCtx_single_good_Ids;...
        [singleUnitsIds(all_Octx_UnitsId),...
        ones(length(all_Octx_UnitsId), 1)*mn]];
    %% sorting units
    Octx_unitsWithSigSlope =...
        (slopePvalVecMean' <= pVal_cuttOff) &...
        (singleUnitsFr > frThreshold) &...
        ismember(singleUnitsIds, OCtc_unitsId);
    Octx_unitsWithSigSlopeId = find(Octx_unitsWithSigSlope);
    %%

    slopeVec_Octx_unitsWithSigSlope =...
        slopeVecMean(Octx_unitsWithSigSlope);
    rasterMap_Octx_unitsWithSigSlope =...
        rasterMapsMean(Octx_unitsWithSigSlope, :);
    [lr_slopeSorted, lr_slopeSorted_id]=...
        sort(slopeVec_Octx_unitsWithSigSlope, 'descend');
    rasterMap_Octx_unitsWithSigSlope_sorted =...
        rasterMap_Octx_unitsWithSigSlope(...
        lr_slopeSorted_id, :);
    %----
    imAlpha=ones(size(rasterMap_Octx_unitsWithSigSlope_sorted));
    imAlpha(isnan(rasterMap_Octx_unitsWithSigSlope_sorted))=0;
    close all; figure;
    imagesc(rasterMap_Octx_unitsWithSigSlope_sorted,...
        'AlphaData',imAlpha)
    caxis(cRange)
    colorbar
    colormap(bluewhitered)
    box off
    xticks(1:length(afRAnge))
    yticks([])
    xticklabels(afRAnge)
    xtickangle(45)
    mouseManeForTitle = split(mouseName, '_');
    mouseManeForTitle = sprintf('%s(%s)',...
         mouseManeForTitle{1},  mouseManeForTitle{2});     
    title({mouseManeForTitle,...
        ''})
    xlabel('airflow (ml/s)')
    ylabel(sprintf('Unit Ids (%d units)', length(Octx_unitsWithSigSlopeId)));
    set(gcf,'Position',[100 100 300 650]);
    print_it(pltDir, 'rasterMap_AI_exp', mouseName)
end
%%
Octx_unitsWithSigSlope =...
    pooled_pVal <= pVal_cuttOff;

slopeVec_Octx_unitsWithSigSlope =...
    pooled_slope(Octx_unitsWithSigSlope);
rasterMap_Octx_unitsWithSigSlope =...
    pooled_rasterMaps(Octx_unitsWithSigSlope, :);
[lr_slopeSorted, lr_slopeSorted_id]=...
    sort(slopeVec_Octx_unitsWithSigSlope, 'descend');
rasterMap_Octx_unitsWithSigSlope_sorted =...
    rasterMap_Octx_unitsWithSigSlope(...
    lr_slopeSorted_id, :);
%----

close all; figure;
imAlpha=ones(size(rasterMap_Octx_unitsWithSigSlope_sorted));
imAlpha(isnan(rasterMap_Octx_unitsWithSigSlope_sorted))=0;
imagesc(rasterMap_Octx_unitsWithSigSlope_sorted, 'AlphaData',imAlpha)
caxis(cRange)
colorbar
colormap(bluewhitered)
box off
xticks(1:length(afRAnge))
yticks([])
xtickangle(45)
mouseManeForTitle = split(mouseName, '_');
mouseManeForTitle = sprintf('%s(%s)',...
     mouseManeForTitle{1},  mouseManeForTitle{2});     
title({'Pooled all neurons',...
    ''})
xlabel('airflow')
ylabel(sprintf('Unit Ids (%d units)',...
    sum(Octx_unitsWithSigSlope)));
set(gcf,'Position',[100 100 300 650]);
print_it(pltDir, 'rasterMap_AI_exp', 'pooeld')

%%
plot_bin_size = .03;
close all;
figure;
histogram(pooled_slope, -.4:plot_bin_size:.7,...
    'EdgeAlpha', 0, 'FaceColor', [.3,.3,.3]);
hold on
histogram(pooled_slope(Octx_unitsWithSigSlope),...
     -.4:plot_bin_size:.7,...
    'EdgeAlpha', 0, 'FaceColor', [0,.5,.1]);
% legend({'', 'significant slops'})
ylabel('frequency')
xlabel('linear regression slope')
ylim([0, 50])
box off
legend({'not significant', 'significant'}, 'Location', 'bestoutside')
legend boxoff
set(gcf,'Position',[100 100 700 250]);
print_it(pltDir, 'lin_reg_slope_distribution_AI_exp', 'pooled')
%%
fileID = fopen(fullfile(pltDir,...
            'linear_regression_slope_data.txt'),'w');
fprintf(fileID, '%d of %d niuts have significant lin. reg. slpoe. (p-val cutoff = %.5f) \n',...
    sum(Octx_unitsWithSigSlope), length(pooled_slope), pVal_cuttOff);
fprintf(fileID, 'Number of units with a significant positive lin. reg. slope = %d. \n',...
    sum(pooled_slope(Octx_unitsWithSigSlope)>0));
fprintf(fileID, 'Number of units with a significant negetive lin. reg. slope = %d. \n',...
    sum(pooled_slope(Octx_unitsWithSigSlope)<0));
fclose(fileID);
%%
[colorSet, SandR_colors] = set_plot_seting(15, 4);
pltDir = fullfile(final_figs_path('uniTn1'),...
    'fig5_I');
if ~exist(pltDir , 'dir')
   mkdir(pltDir)
end
%%
mn = 3;
sui = 23;
mouseName = mouseNameList{mn};
%% load data 
load(fullfile(savingDir, mouseName,...
    sprintf('%s_eventsInfo.mat', mouseName)));
load(fullfile(savingDir, mouseName,...
    sprintf('%s_spikeData.mat', mouseName)));
load(fullfile(savingDir, mouseName,...
    sprintf('%s_expParams.mat', mouseName)));
%% get single units spikes time
[singleUnits, singleUnitsIds, singleUnitsFr] =...
    get_singleUnits(spikeStruct);

%% get PSTH
afterOnsetWind_PITH = 1;
unitST = singleUnits{sui}.unitSt;
[raster_mat, unitBaseLine_fr, evoke_raster_mat] =...
    get_raster_mat(unitST, eventsInfo(:, 1), beforeOnsetWind,...
    afterOnsetWind_PITH, beforeOnsetWind);
%%--
PITH_time = -beforeOnsetWind:.001:afterOnsetWind_PITH;
unitPITHs = nan(length(PITH_time),...
    length(afRAnge));
af_idx = 1;
evokeRespToAirflowMat = [];
for af = afRAnge
    unitPITHs( :, af_idx) =...
        gen_fx_gsmooth(mean(evoke_raster_mat(...
            eventsInfo(:, 3) == af, :))*1000,...
            smoothingKurnel);
    af_idx = af_idx+1;
    %%
    af_trials_deltafr = trials_fr_from_raster(...
        raster_mat(eventsInfo(:, 3) == af, :),...
        round(beforeOnsetWind*1000) + 1,...
        round((beforeOnsetWind+afterOnsetWind)*1000))-...
        trials_fr_from_raster(...
        raster_mat(eventsInfo(:, 3) == af, :),...
        1, round(beforeOnsetWind*1000));
    evokeRespToAirflowMat = [evokeRespToAirflowMat,...
        af_trials_deltafr];
end
%% ---
close all
colors = cbrewer2('seq', 'Blues', 10);
colors = colors(5:end,:);
numOfRows = 0;
afi = 1;
figure; hold on;
for af = afRAnge
    af_ids = find(eventsInfo(:, 3) == af);
    af_ids = af_ids(randperm(length(af_ids)));
    [rowRaster,colRaster] =...
        find(raster_mat(af_ids, :)==1);
    plot(colRaster, rowRaster+numOfRows,'.',...
        'Color', colors(afi, :),'MarkerSize', 2);
    numOfRows = numOfRows + sum(eventsInfo(:, 3) == af);
    plot(xlim, [numOfRows, numOfRows],...
        'LineStyle', ':','Color', 'k', 'LineWidth', 1);
    afi = afi + 1;
end
yl_PITH = ylim;
plot([500, 500], yl_PITH, 'LineStyle', '--','Color', 'k',  'LineWidth', 1.5);
fill([500,500, 500+1000*stimuliDuration, 500+1000*stimuliDuration],...
    [yl_PITH, yl_PITH(2), yl_PITH(1) ],...
    [.5, .5, .5],...
    'EdgeAlpha', 0, 'FaceAlpha', .3)
yticks([])
xticks([1, 501, 1001, 1501])
ylim(ylim-1)
xticklabels(PITH_time([1, 501, 1001, 1501]))
xlabel('time (sec)')
set(gcf,'Position',[100 100 400 400]);
print_it(pltDir,sprintf('raster_SUI%d', sui), mouseName)
%%
lineWidth = 2.5;
close all
figure;
hold on
for bi = 1 :  expParams.numOFAirflows
    plot(PITH_time,...
        unitPITHs(:, bi),'Color', colors(bi,:),...
        'LineWidth', lineWidth)
end
cb = colorbar('FontSize', 12 , 'Box',...
'off', 'Location', 'eastoutside');
colormap(colors)
p  = patchline([0,0],ylim,'linestyle','--','edgecolor','k',...
'linewidth',2,'edgealpha',0.7);
xticks([-.5, 0, .5, 1])
xlabel('time(sec)');
ylabel('\Delta firing rate (sp./s)')
box off
set(gcf,'Position',[100 100 400 400]);
print_it(pltDir,sprintf('PITH_SUI%d_AI_exp', sui), mouseName)
%%

evokeRespToAirflowMat_zs=...
    (evokeRespToAirflowMat-...
     nanmean(evokeRespToAirflowMat, [1,2]))/...
     nanstd(evokeRespToAirflowMat, [], [1,2]);

neuon_mech_unning_curve = mean(evokeRespToAirflowMat_zs);
neuon_mech_unning_curve_sem = std(evokeRespToAirflowMat_zs)/...
    sqrt(size(evokeRespToAirflowMat_zs, 1));

close all
errorbar(afRAnge, neuon_mech_unning_curve,...
    neuon_mech_unning_curve_sem,...
    'Marker', 'o', 'Color', colorSet(2,:),...
    'CapSize', 3, 'LineWidth', lineWidth, 'MarkerSize', 4)
box off
xtickangle(45)
xlim([25, 300])
xlabel('inh. airflow (ml/s)')
ylabel('z-score firing rate (un)')
title(sprintf('SUI = %d', sui))
set(gcf,'Position',[100 100 500 300]);
print_it(pltDir,'tuning_curve_AI_exp', num2str(i))