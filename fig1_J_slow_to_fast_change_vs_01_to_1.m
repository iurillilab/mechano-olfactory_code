clear; clc; close all;
%% 
%This script may not work with the storage data. So then you should fixt
%the real trial id index to reduc ambiguty and inrease readebilty and then
%run it

clear; clc; close all;
[~, mouseNameList, dataPathPrefixList] =....
    get_all_conc_mice('paola');
%%
[tBeginning, tEnd, trialsToRemove]=...
    get_removal_variables();
pVal_threshold = 0.05;
lowConc = 1;
highConc = 2;
conc_set = [lowConc, highConc];
inhInducedChange = [];
concInducedChange = [];
diffRespTestResults = [];

count = 0;

%%
for mn = 1 : length(mouseNameList)
    %%
    mouseName = mouseNameList{mn};
    dataPathPrefix = dataPathPrefixList{mn};
    mouseDirName = fullfile(dataPathPrefix,'processedDataStorage',...
        mouseName);
    load(fullfile(mouseDirName,...
        sprintf('%s_respData.mat', mouseName)));
    load(fullfile(mouseDirName,...
        sprintf('%s_respLabelsData.mat', mouseName)));
    load(fullfile(mouseDirName,...
        sprintf('%s_inh_NR.mat', mouseName)));
    load(fullfile(mouseDirName,...
        sprintf('%s_unitsBaselineFR.mat', mouseName)));
    load(fullfile(mouseDirName,...
        sprintf('%s_recordingData.mat', mouseName)));
    load(fullfile(mouseDirName,...
        sprintf('%s_odorStatistics_wholeStimulus.mat', mouseNameList{mn})))
    %%
    if isfield(respData,'rawXpeakPuffFree')
        xpeaksTrialId = respData.xpeakRealTrialId;
    else
        xpeaksTrialId = respData.xpeakTrialId;
    end
    %%
    goodAndSafe_unitsFlag =...
        get_goodAndSafe_units(dataPathPrefix, mouseName,...
        6, .5);
    %%   
    avg_fr = inh_NR.after_sc/inh_NR.cacheInfo.windSizeSec;
    %%
    pvals = odorStatistics_wholeStimulus.graterTestPVal(:,:,conc_set); %not limiting the number of odors included
    pvals = reshape(pvals, size(pvals,1), [], 1);
    pvals = pvals';
    clear c_pvals
    for i = 1:size(pvals,2)
        [c_pvals(:,i), ~, ~, ~] = fdr_BH(pvals(:,i), pVal_threshold);
    end
    c_results = sum(c_pvals < pVal_threshold) > 0;
    c_results = c_results(:);
    
    % revert c_pvals to original form
    c_pvals = c_pvals';
    c_pvals = reshape(c_pvals, size(c_pvals,1), [], length(conc_set));
    
%     odor_responsive_flag = c_results(goodAndSafe_unitsFlag); use later to concatenate for ct indexing
    
    odorResponsiveUnits = goodAndSafe_unitsFlag & c_results;
    sum(odorResponsiveUnits)
    
    %%
    for ni = find(odorResponsiveUnits)'
        for oId = 1 : max(respData.xpeakOdorId)
%             c_pvals(ni, oId, conc_set)
            if sum(c_pvals(ni, oId, conc_set)<pVal_threshold) >= 1
%             if((sum(odorStatistics_wholeStimulus.graterTestResult(...
%                 ni, oId, conc_set), 3) >= 1)&...
%                 goodAndSafe_unitsFlag(ni))
            count = count+1;
                lowConcSlow_FR = avg_fr(xpeaksTrialId >...
                    trialsToRemove &...
                    respData.xPeaksTimeLog > tBeginning &...
                    respData.xPeaksTimeLog <= tEnd &...
                    respLabelsData.respLabels == 2 &...
                    respData.xpeakConcId == lowConc &...
                    respData.xpeakOdorId == oId, ni);
                highConcSlow_FR = avg_fr(xpeaksTrialId >...
                    trialsToRemove &...
                    respData.xPeaksTimeLog > tBeginning &...
                    respData.xPeaksTimeLog <= tEnd &...
                    respLabelsData.respLabels == 2 &...
                    respData.xpeakConcId == highConc &...
                    respData.xpeakOdorId == oId, ni);
                
                lowConcFast_FR = avg_fr(xpeaksTrialId >...
                    trialsToRemove &...
                    respData.xPeaksTimeLog > tBeginning &...
                    respData.xPeaksTimeLog <= tEnd &...
                    respLabelsData.respLabels == 1 &...
                    respData.xpeakConcId == lowConc &...
                    respData.xpeakOdorId == oId, ni);

                concInducedChange_cell = mean(highConcSlow_FR)-...
                    mean(lowConcSlow_FR);
                inhInducedChange_cell = mean(lowConcFast_FR)-...
                    mean(lowConcSlow_FR);
                [pVal_diff] =...
                    ranksum(lowConcSlow_FR, highConcSlow_FR,...
                    'tail', 'both');
                %%
                inhInducedChange = [inhInducedChange;...
                    inhInducedChange_cell];
                concInducedChange = [concInducedChange;...
                    concInducedChange_cell]; 
                diffRespTestResults =[diffRespTestResults;...
                    pVal_diff <= pVal_threshold];
            end    
        end
    end
end
disp(count)
%%
addpath('/Users/galileo/GitHub/REZA/sniff-odor-project/projectOwnToolBox')
addpath(genpath('/Users/galileo/Dropbox (Personal)/SniffOdor_project_processedData/sniffOdorProject/toolBox/Tools-master/Tools-master/plotting'))


[colorSet, SandR_colors] =...
    set_plot_seting(15, 4);
pltDir = fullfile(final_figs_path('paola'),...
    'fig1_K');
if ~exist(pltDir , 'dir')
   mkdir(pltDir)
end
fileId = fopen(fullfile(pltDir,...
    'fig1_K_num_and_stat.txt'),...
    'w');
%%
markerSize = 15;
axisLim = [-.5, 25];
transparancyAlfa = .5;
dotFaceColor = .4*[1,1,1];
dotEdgeColor = .1*[1,1,1];
pltTiks = 0:10:25;

close all; figure;
hold on;
mdl = fitlm(abs(inhInducedChange),...
        abs(concInducedChange));  
scatter(abs(inhInducedChange(diffRespTestResults ~= 1)),...
    abs(concInducedChange(diffRespTestResults ~= 1)),markerSize,...
    'MarkerEdgeColor',dotEdgeColor, 'MarkerFaceColor',[1,1,1],...
   'MarkerFaceAlpha',transparancyAlfa,'MarkerEdgeAlpha',transparancyAlfa)

mdl_justSig = fitlm(abs(inhInducedChange(diffRespTestResults == 1)),...
        abs(concInducedChange(diffRespTestResults == 1)));  
    
scatter(abs(inhInducedChange(diffRespTestResults == 1)),...
    abs(concInducedChange(diffRespTestResults == 1)),markerSize,...
    'MarkerEdgeColor',dotEdgeColor, 'MarkerFaceColor',dotFaceColor,...
   'MarkerFaceAlpha',transparancyAlfa,'MarkerEdgeAlpha',transparancyAlfa)

plot_lr_lineAnd_inerval(mdl, axisLim')
plot_lr_lineAnd_inerval(mdl_justSig, axisLim', [],...
    colorSet(2,:), colorSet(2,:), .2)
xlabel('delta firing rate, fast to Slow')
ylabel('delta firing rate, .1% to 1%')
xlim(axisLim); ylim(axisLim);
xticks(pltTiks); yticks(pltTiks)
axis square
set(gcf,'Position',[350 350 400 410]);
print_it(pltDir, 'fig1_K_6_mice_with_2_conc','pooled')
%%
axisLimZoom = [-.2, 10];
xlim(axisLimZoom); ylim(axisLimZoom);
xticks(0:3:10); yticks(0:3:10)
print_it(pltDir, 'fig1_K_6_mice_with_2_conc_zoom_in','pooled')
%---
fprintf(fileId,'all units linear reg: r2 = %.5f, p_val = %.5f. \n',...
    mdl.Rsquared.Ordinary, mdl.Coefficients.pValue(2));
fprintf(fileId, 'just sig. units linear reg: r2 = %.5f, p_val = %.5f. \n',...
    mdl_justSig.Rsquared.Ordinary, mdl_justSig.Coefficients.pValue(2));
%%
axisLim = [-3, 35];
yLim = [-.23, .23];
close all
figure
raincloud_plot(abs(inhInducedChange),...
    'color', .5*[1,1,1], 'box_on', 1, 'alpha', .5,...
    'lwr_bnd', 1, 'line_width', 1.5);
box off
xlabel('abs delta firing rate, fast to slow')
title('slow to fast change')
xlim(axisLim);
ylim(yLim)
set(gcf,'Position',[350 350 500 350])
print_it(pltDir, 'fig1_K_slow2fast_raincloud_6_mice_with_2_conc','pooled')
%%--
close all
raincloud_plot(abs(concInducedChange),...
    'color', .5*[1,1,1], 'box_on', 1, 'alpha', .5,...
    'lwr_bnd', 1, 'line_width', 1.5);
box off
xlabel('abs delta firing rate, .1% to 1%')
title('.1% to 1% change')
xlim(axisLim);
ylim(yLim)
set(gcf,'Position',[350 350 500 350])
print_it(pltDir, 'fig1_K_01to1_raincloud_6_mice_with_2_conc','pooled')
% print_it(pltDir, 'fig1_K_01to1_raincloud_6_mice_with_2_conc_combined','pooled')

%% --
close all
raincloud_plot(abs(inhInducedChange(diffRespTestResults == 1)),...
    'color', colorSet(2, :), 'box_on', 1, 'alpha', .5,...
    'lwr_bnd', 1, 'line_width', 1.5);
box off
xlabel('abs delta firing rate, fast to slow')
title({'slow to fast change',...
    'significant difference response conc. change',''})
xlim(axisLim);
ylim(yLim)
set(gcf,'Position',[350 350 500 400])
print_it(pltDir, 'fig1_K_slow2fast_sig_diff_raincloud_6_mice_with_2_conc','pooled')
%%--
close all
raincloud_plot(abs(concInducedChange(diffRespTestResults == 1)),...
    'color', colorSet(2, :), 'box_on', 1, 'alpha', .5,...
    'lwr_bnd', 1, 'line_width', 1.5);
box off
xlabel('abs delta firing rate, .1% to 1%')
title({'.1% to 1% change',...
    'significant difference response conc. change',''})
xlim(axisLim);
ylim(yLim)
set(gcf,'Position',[350 350 500 400])
print_it(pltDir, 'fig1_K_01to1_sig_diff_raincloud_6_mice_with_2_conc','pooled')
%%
[pVal_all_units] =...
    ranksum(abs(inhInducedChange),...
        abs(concInducedChange),...
    'tail', 'both');
[pVal_just_sig] =...
    ranksum(abs(inhInducedChange(diffRespTestResults == 1)),...
        abs(concInducedChange(diffRespTestResults == 1)),...
    'tail', 'both');

[h, p] = ttest(abs(concInducedChange),abs(inhInducedChange))
  pVal_all_units 
  pVal_just_sig
%%
fprintf(fileId,'all units rank sum: p_val = %.5f. \n',...
    pVal_all_units);
fprintf(fileId, 'just sig. rank sum: p_val = %.5f. \n',...
    pVal_just_sig);
fclose(fileId);