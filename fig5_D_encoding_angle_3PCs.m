clear; clc; close all;
mouseNameList = {'M1_102721', 'M2_102821','M3_102921',...
    'M5_120221'};
% dataPathPrefix = 'E:\Ephys\Conc_Series';
dataPathPrefix = 'D:\Reza\sniffOdorProject\Conc_Series';
%%
[tBeginning, tEnd, trialsToRemove]=...
    get_removal_variables();
add_ndt_paths_and_init_rand_generator
%%
conc_set = [1,2,4];
numberOfBalanceSampling = 100;
concVec_Stype_vs_samplingVec =...
    nan(numberOfBalanceSampling, 2);
concVec_Rtype_vs_samplingVec =...
    nan(numberOfBalanceSampling, 2);
concVec_Stype_vs_concVec_Rtype =...
    nan(numberOfBalanceSampling, 2);
%%
filesRoot = 'C:\Users\a.abolghasemi\Documents\MATLAB';
addpath(genpath(fullfile(filesRoot, 'drtoolbox')))
for perm = 1 : numberOfBalanceSampling
        %%
        [fSpaceCell, fSpaceLbelsCell] =...
            get_aSample(dataPathPrefix, mouseNameList,....
                        conc_set, 1:2, tBeginning, tEnd,trialsToRemove,...
                        1, 0, 6, .5);
         %%
        pseudoPopulation =...
            [fSpaceCell{1}, fSpaceCell{2},fSpaceCell{3}, fSpaceCell{4}];
        [score,pcaMap] = pca(zscore(pseudoPopulation),3);
        fSpaceLbels = fSpaceLbelsCell{1};
        %%
        for oId = 1 : 2
            %%
            labels_Stype =...
                (fSpaceLbels.OdorId == oId &...
                fSpaceLbels.RespLabel == 1);
            
            label_Rtype_flag =...
                (fSpaceLbels.OdorId == oId &...
                fSpaceLbels.RespLabel == 1);
            %%
            [~, mapping] = lda(score(...
                fSpaceLbels.OdorId == oId &...
                fSpaceLbels.RespLabel == 1, :),...
                fSpaceLbels.ConcId(...
                fSpaceLbels.OdorId == oId &...
                fSpaceLbels.RespLabel == 1), 1);
            
            concVec_Stype = mapping.M(:, 1);

            [~, mapping] = lda(score(...
                fSpaceLbels.OdorId == oId &...
                fSpaceLbels.RespLabel == 2, :),...
                fSpaceLbels.ConcId(...
                fSpaceLbels.OdorId == oId &...
                fSpaceLbels.RespLabel == 2), 1);
            
            concVec_Rtype = mapping.M(:, 1);
            
            [~, mapping] = lda(score(...
                fSpaceLbels.OdorId == oId, :),...
                fSpaceLbels.RespLabel(...
                fSpaceLbels.OdorId == oId),1);
            samplingVec = mapping.M(:, 1);
            %%
            concVec_Stype_vs_samplingVec(perm, oId) =...
                real(acos(dot(concVec_Stype, samplingVec)/...
                (norm(concVec_Stype)* norm(samplingVec)) ));
            
            concVec_Rtype_vs_samplingVec(perm, oId) =...
                real(acos( dot(concVec_Rtype, samplingVec)/...
                (norm(concVec_Rtype)* norm(samplingVec)) ));
            
            concVec_Stype_vs_concVec_Rtype(perm, oId) =...
                real(acos( dot(concVec_Rtype, concVec_Stype)/...
                (norm(concVec_Rtype)* norm(concVec_Stype)) ));      
         end
end
rmpath(genpath(fullfile(filesRoot, 'drtoolbox')))
%%
pltDir = fullfile(final_figs_path('uniTn1'),...
    'Fig5_D');
if ~exist(pltDir , 'dir')
   mkdir(pltDir)
end
[colorSet, SandR_colors] =...
    set_plot_seting(15, 8);
%%
conditions = {'odor 1', 'odor 2', 'pooled-two-odor'};
for i = 1: 3
    %% calculate transferred angle
    switch i 
        case 1
            concEncodingsAngleTransferred =...
                abs(90- rad2deg(concVec_Stype_vs_concVec_Rtype(:,1)));
            concVsSampelingAngleTransferred =...
                abs(90- rad2deg([concVec_Stype_vs_samplingVec(:,1);...
                concVec_Rtype_vs_samplingVec(:,1)]));
        case 2
            concEncodingsAngleTransferred =...
                abs(90- rad2deg(concVec_Stype_vs_concVec_Rtype(:,2)));
            concVsSampelingAngleTransferred =...
                abs(90- rad2deg([concVec_Stype_vs_samplingVec(:,2);...
                concVec_Rtype_vs_samplingVec(:,2)]));
        case 3
            concEncodingsAngleTransferred =...
                abs(90- rad2deg(concVec_Stype_vs_concVec_Rtype(:)));
            concVsSampelingAngleTransferred =...
                abs(90- rad2deg([concVec_Stype_vs_samplingVec(:);...
                concVec_Rtype_vs_samplingVec(:)]));
    end

    axisLim = [-5, 95];
    yLim = [-.1, .1];
    close all
    raincloud_plot(concEncodingsAngleTransferred,...
        'color', colorSet(1,:), 'box_on', 1, 'alpha', .5,...
        'lwr_bnd', 1, 'line_width', 1.5);
    box off
    xlabel('|90 - angle|')
    title({'slow vs. fast conc. encoding angle',...
        conditions{i}})
    xlim(axisLim);
    ylim(yLim)
    set(gcf,'Position',[350 350 500 350])
    print_it(pltDir,...
        sprintf('slow_and_fast_conc_encoding_angle_%s',...
        conditions{i}),'pooled')

    close all
    raincloud_plot(concVsSampelingAngleTransferred,...
        'color', colorSet(2,:), 'box_on', 1, 'alpha', .5,...
        'lwr_bnd', 1, 'line_width', 1.5);
    box off
    xlabel('|90 - angle|')
    title({'conc. encoding vs inh. speed encoding angle',...
        conditions{i}})
    xlim(axisLim);
    ylim(yLim)
    set(gcf,'Position',[350 350 500 350])
    print_it(pltDir,...
        sprintf('conc_encoding_vs_inh_speed_encoding_angle_%s',...
        conditions{i}),'pooled')

    %%
    fileID = fopen(fullfile(pltDir,...
            sprintf('consDir_samplingDir_angle_%s.txt',...
            conditions{i})),'w');

    fprintf(fileID, 'R-type vs. S-type conc encoding transferd angle: mean: %.3f, std: %.5f; median: %.5f.\n',....
        (mean(concEncodingsAngleTransferred)),...
        (std(concEncodingsAngleTransferred)),...
        (median(concEncodingsAngleTransferred)));

    fprintf(fileID, 'conc. encoding vs inh. speed encoding transferd angle: mean: %.3f, std: %.5f; median: %.5f.\n',....
         (mean(concVsSampelingAngleTransferred)),...
        (std(concVsSampelingAngleTransferred)),...
        (median(concVsSampelingAngleTransferred)));

    [testResults, ~] = ranksum(concEncodingsAngleTransferred,...
        concVsSampelingAngleTransferred);
    fprintf(fileID, 'rank sum p_Val of transferd angle: %.6f.\n', testResults);
    fclose(fileID);
end