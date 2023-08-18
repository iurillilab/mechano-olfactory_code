clear; clc; close all;
%% 
%This script may not work with the storage data. So then you should fixt
%the real trial id index to reduc ambiguty and inrease readebilty and then
%run it

clear; clc; close all;
[dataFilesDirList, mouseNameList, dataPathPrefixList] =....
    get_all_conc_mice('rezaExternal');

%%
[tBeginning, tEnd, trialsToRemove]=...
    get_removal_variables();
conc_set = [1,2];
num_of_odors = 2;
%%
brforeInhOnsetWind = .6;
afterInhOnsetWind = .6;
PITH_time = -1*brforeInhOnsetWind:...
     .001 : afterInhOnsetWind +.001;
smoothingCurnel = 10;

%%
[colorSet, SandR_colors] = set_plot_seting(15, 4);
pltDir = fullfile(final_figs_path('uniTn1'),...
    'fig1_I_and_J', 'example3');
if ~exist(pltDir , 'dir')
   mkdir(pltDir)
end

%% selection 

ni = 79;
oId = 3;
mn = 5;
ylimRange = [-10,60];
%% the best ones
% mn = 3;
% ni = 67;
% oId = 1;

% mn = 5;
% ni = 76;
% oId = 4;

% mn = 5;
% ni = 79;
% oId = 3;
%%
mouseName = mouseNameList{mn};
dataPathPrefix = dataPathPrefixList{mn};
mouseDirName = fullfile(dataPathPrefix,'processedDataStorage',...
    mouseName);
respRasterDir = fullfile(mouseDirName, 'RstersFiles\inhRasters');
TrialsDir = fullfile(mouseDirName,'RstersFiles\allTrialsRasters');
load(fullfile(mouseDirName,...
    sprintf('\\%s_respData.mat', mouseName)));
load(fullfile(mouseDirName,...
    sprintf('\\%s_respLabelsData.mat', mouseName)));
load(fullfile(mouseDirName,...
    sprintf('\\%s_inh_NR.mat', mouseName)));
load(fullfile(mouseDirName,...
    sprintf('\\%s_unitsBaselineFR.mat', mouseName)));
load(fullfile(mouseDirName,...
    sprintf('\\%s_recordingData.mat', mouseName)));
load(fullfile(mouseDirName,...
    sprintf('\\%s_odorStatistics_wholeStimulus.mat', mouseNameList{mn})))
%%    
unitRasterData = load([respRasterDir,...
    sprintf('\\%s_NPX_unit%d_resps_raster_data.mat',...
    mouseName, ni)]);
unitTrialRasterData = load([TrialsDir,...
    sprintf('\\%s_NPX_unit%d_allTrials_raster_data.mat',...
    mouseNameList{mn},ni)]);
cell_conc_fast_PITH = cell(4,1);
cell_conc_slow_PITH = cell(4,1);
cell_conc_fast_raster = cell(4,1);
cell_conc_slow_raster = cell(4,1);

 for ci = conc_set

    slowInhRasters_conc =...
        unitRasterData.raster_data(...
        unitRasterData.raster_labels.xpeakTrialId >...
        trialsToRemove &...
        unitRasterData.raster_labels.xPeaksTimeLog >...
        tBeginning &...
        unitRasterData.raster_labels.xPeaksTimeLog <=...
        tEnd &...
        unitRasterData.raster_labels.labels == 2 &...
        unitRasterData.raster_labels.xpeakConcId == ci &...
        unitRasterData.raster_labels.xpeakOdorId == oId, :);
    fastInhRasters_conc =...
        unitRasterData.raster_data(...
        unitRasterData.raster_labels.xpeakTrialId >...
        trialsToRemove &...
        unitRasterData.raster_labels.xPeaksTimeLog >...
        tBeginning &...
        unitRasterData.raster_labels.xPeaksTimeLog <=...
        tEnd &...
        unitRasterData.raster_labels.labels == 1 &...
        unitRasterData.raster_labels.xpeakConcId == ci &...
        unitRasterData.raster_labels.xpeakOdorId == oId, :);
   % test cell-odor prefrence
   
   slow_PITH_conc = gen_fx_gsmooth(...
       mean(slowInhRasters_conc*1000 - unitsBaselineFR(ni)),...
        smoothingCurnel);

   fast_PITH_conc = gen_fx_gsmooth(...
       mean(fastInhRasters_conc*1000 - unitsBaselineFR(ni)),...
        smoothingCurnel);

   cell_conc_fast_PITH{ci} =...
     fast_PITH_conc;
   cell_conc_slow_PITH{ci} =...
     slow_PITH_conc; 
   cell_conc_fast_raster{ci} = fastInhRasters_conc;
   cell_conc_slow_raster{ci} = slowInhRasters_conc;
 end
%%
%%

xlimRange = [-.3, .4];
lineWidth = 2.5;
rstar_dot_size = 4;
%%--
close all
hold on
rasterHigh = 0;
rasteMeadel = [];
for ci = conc_set
    [row_slow,col_slow] =...
          find(cell_conc_slow_raster{ci} == 1);
    plot(col_slow, row_slow+rasterHigh,'.', 'Color',...
        SandR_colors.r_shade(ci,:), 'MarkerSize', rstar_dot_size);
    
    rasteMeadel = [rasteMeadel;...
        rasterHigh+...
        (size(cell_conc_slow_raster{ci}, 1)/2)];
    rasterHigh = rasterHigh+...
        size(cell_conc_slow_raster{ci}, 1);
    plot([0,1200], [rasterHigh, rasterHigh], 'k:',...
        'LineWidth', 1)
end
xlim(xlimRange*1000+600)
ylim([-1,rasterHigh-1])
yticks(rasteMeadel)
xticks([-.2,0,.2,.4]*1000+600)
xtickangle(45)
ytickangle(45)
yticklabels({'.01%','1%'})
xticklabels({'-0.2', '0', '0.2', '0.4'})
box off
set(gcf,'Position',[200 640 300 150]);
print_it(pltDir, sprintf('01_to_1_conc_slow_raster_ni%d_oId%d', ni, oId),...
    mouseName)
%%--
figure
hold on
for ci = conc_set
    plot(PITH_time, cell_conc_slow_PITH{ci},...
        'Color',SandR_colors.r_shade(ci, :),...
        'linewidth',lineWidth);
end

p  = patchline([0,0],[-50,150],'linestyle','--','edgecolor','k',...
'linewidth',2,'edgealpha',0.7);
ylim(ylimRange)
xlim(xlimRange)
xticks([-.2,0,.2,.4])
xlabel('Time (sec)');
ylabel('Firing rate (Hz)')
box off
set(gcf,'Position',[200 200 300 350]);
print_it(pltDir, sprintf('01_to_1_conc_slow_PITH_ni%d_oId%d', ni, oId),...
    mouseName)
%% ---
figure;
hold on
plot(PITH_time,...
    cell_conc_slow_PITH{1},'Color',...
    SandR_colors.r_shade(1,:), 'LineWidth', lineWidth)
plot(PITH_time,...
    cell_conc_fast_PITH{1},'Color',...
    SandR_colors.s_shade(1,:), 'LineWidth', lineWidth)
p  = patchline([0,0],[-50,150],'linestyle','--','edgecolor','k',...
'linewidth',2,'edgealpha',0.7);
ylim(ylimRange)
xlim(xlimRange)
xticks([-.2,0,.2,.4])
xlabel('time(sec)');
ylabel('\Delta firing rate (sp./s)')
box off

set(gcf,'Position',[820 200 300 350]);
print_it(pltDir, sprintf('slow_to_fast_PITH_ni%d_01Conc', ni),...
    mouseName)
%---
figure
hold on

hold on
[rowRtype,colRtype] = find(...
    cell_conc_slow_raster{1}==1);
[rowStype,colStype] = find(...
    cell_conc_fast_raster{1}==1);

plot(colRtype, rowRtype+max(rowStype),'.', 'Color',...
    SandR_colors.r_shade(1,:), 'MarkerSize', rstar_dot_size);
plot(colStype, rowStype,'.',...
    'Color', SandR_colors.s_shade(1,:),'MarkerSize', rstar_dot_size);

p  = patchline([600,600],...
    [0,max(rowStype)+max(rowRtype)],...
    'linestyle','--','edgecolor','k',...
'linewidth',2,'edgealpha',0.7);
xticks([100,300,500,700]+300)
xticklabels({'-0.2', '0', '0.2', '0.4', '0.6'})
yticks([max(rowStype)/2,...
    max(rowStype)+max(rowRtype)/2])
yticklabels({'Fast', 'Slow'})
ytickangle(45)
xlim([300, 1000])
ylim([1, max(rowStype)+max(rowRtype)]) 
xlabel('time (s)')
set(gcf,'Position',[800 640 300 150]);
print_it(pltDir, sprintf('slow_to_fast_Raster_ni%d_no_odor', ni),...
    mouseName)