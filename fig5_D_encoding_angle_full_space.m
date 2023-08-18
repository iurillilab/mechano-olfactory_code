clear; clc; close all;

dataDir = fullfile('G:\My Drive',...
    'sinff-odor-project\Figuers_data_Marh2023');
load(fullfile(dataDir, 'fig5_encoding_angle_full_space_data',...
    'encoding_angle_full_space_data'));
%%
pltDir = fullfile(final_figs_path('uniTn1'),...
    'Fig5_D_full_space');
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
                90 - concVec_Stype_vs_concVec_Rtype(:,1);
            concVsSampelingAngleTransferred =...
                90 - [concVec_Stype_vs_samplingVec(:,1);...
                concVec_Rtype_vs_samplingVec(:,1)];
        case 2
            concEncodingsAngleTransferred =...
                90 - concVec_Stype_vs_concVec_Rtype(:,2);
            concVsSampelingAngleTransferred =...
                90 - [concVec_Stype_vs_samplingVec(:,2);...
                concVec_Rtype_vs_samplingVec(:,2)];
        case 3
            concEncodingsAngleTransferred =...
                90 - concVec_Stype_vs_concVec_Rtype(:);
            concVsSampelingAngleTransferred =...
                90 - [concVec_Stype_vs_samplingVec(:);...
                concVec_Rtype_vs_samplingVec(:)];
    end

    axisLim = [0, 95];
    close all
    raincloud_plot(concEncodingsAngleTransferred,...
        'color', colorSet(1,:), 'box_on', 1, 'alpha', .5,...
        'lwr_bnd', 1, 'line_width', 1.5);
    box off
    xlabel('|90 - angle|')
    xticks([0,45,90])
    title({'slow vs. fast conc. encoding angle',...
        conditions{i}})
    xlim(axisLim);
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
    xticks([0,45,90])
    title({'conc. encoding vs inh. speed encoding angle',...
        conditions{i}})
    xlim(axisLim);
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