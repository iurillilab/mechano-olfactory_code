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
numOfCons = 4;
num_of_odors = 2;
%%
brforeInhOnsetWind = .6;
afterInhOnsetWind = .6;
PITH_time = -1*brforeInhOnsetWind:...
     .001 : afterInhOnsetWind +.001;
smoothingCurnel = 10;
conc_set = [1,2];
%%
[colorSet, SandR_colors] = set_plot_seting(15, 4);
pltDir = fullfile(final_figs_path('uniTn1'),...
    'fig1_E','example5');
if ~exist(pltDir , 'dir')
   mkdir(pltDir)
end

%% selection 
% mn = 1;
% ni = 10;
% oId = 1;

% mn = 3;
% ni = 86;
% oId = 2;

% mn = 3;
% ni = 88;
% oId = 1;

% mn = 4; %very good
% ni = 99;
% oId = 1;

% mn = 2;%very good
% ni = 48;
% oId = 2;
% mn = 5;%very good
% ni = 75;
% oId = 2;
%%
% mn = 5;
% ni = 52;
% oId = 4;
%%
% ni = 74; selected
% oId = 1;
% mn = 3;
%%
% ni = 99;  %selected
% oId = 1;
% mn = 4;
%%
%  mn = 2;%very good
%  ni = 48;
%  oId = 2;
%%
 mn = 5;
 ni = 52;
 oId = 4;
ylimRange = [-20,100];
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
if isfield(respData,'rawXpeakPuffFree')
    xpeaksTrialId = respData.xpeakRealTrialId;
else
    xpeaksTrialId = respData.xpeakTrialId;
end
%%
unitRasterData = load([respRasterDir,...
    sprintf('\\%s_NPX_unit%d_resps_raster_data.mat',...
    mouseName, ni)]);
cell_conc_fast_PITH = cell(4,1);
cell_conc_slow_PITH = cell(4,1);
cell_conc_fast_raster = cell(4,1);
cell_conc_slow_raster = cell(4,1);

 for ci = 1 : numOfCons

    slowInhRasters =...
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
    fastInhRasters =...
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
   
   slow_PITH = gen_fx_gsmooth(...
       mean(slowInhRasters*1000 - unitsBaselineFR(ni)),...
        smoothingCurnel);

   fast_PITH = gen_fx_gsmooth(...
       mean(fastInhRasters*1000 - unitsBaselineFR(ni)),...
        smoothingCurnel);

   cell_conc_fast_PITH{ci} =...
     fast_PITH;
   cell_conc_slow_PITH{ci} =...
     slow_PITH; 
   cell_conc_fast_raster{ci} = fastInhRasters;
   cell_conc_slow_raster{ci} = slowInhRasters;

 end
%%
avg_fr = inh_NR.after_sc/inh_NR.cacheInfo.windSizeSec;
slowInhFr = avg_fr(xpeaksTrialId >...
    trialsToRemove &...
    respData.xPeaksTimeLog > tBeginning &...
    respData.xPeaksTimeLog <= tEnd &...
    respLabelsData.respLabels == 2 &...
    respData.xpeakConcId == 1 &...
    respData.xpeakOdorId == oId, ni);
fastInhFr = avg_fr(xpeaksTrialId >...
    trialsToRemove &...
    respData.xPeaksTimeLog > tBeginning &...
    respData.xPeaksTimeLog <= tEnd &...
    respLabelsData.respLabels == 1 &...
    respData.xpeakConcId == 1 &...
    respData.xpeakOdorId == oId, ni);
                
avrage_responce = [mean(slowInhFr);...
    mean(fastInhFr)];
%%
neuronSlowRaster_noOdor = unitRasterData.raster_data(...
    unitRasterData.raster_labels.labels == 2 &...
    unitRasterData.raster_labels.xpeakCondLabel == 0, :);

neuronFastRaster_noOdor = unitRasterData.raster_data(...
    unitRasterData.raster_labels.labels == 1 &...
    unitRasterData.raster_labels.xpeakCondLabel == 0, :);

neuronFastRaster_Mix = unitRasterData.raster_data(...
    unitRasterData.raster_labels.xpeakCondLabel == 0, :);

neuronsSlow_PITH_noOdor = gen_fx_gsmooth(...
       mean(neuronSlowRaster_noOdor*1000 - unitsBaselineFR(ni)),...
        smoothingCurnel);
neuronsFast_PITH_noOdor = gen_fx_gsmooth(...
    mean(neuronFastRaster_noOdor*1000 - unitsBaselineFR(ni)),...
        smoothingCurnel);

neuronsMix_PITH_noOdor = gen_fx_gsmooth(...
    mean(neuronFastRaster_Mix*1000 - unitsBaselineFR(ni)),...
        smoothingCurnel);

%%
xlimRange = [-.3, .4];
singleMicePITH_alpa = .25;
lineWidth = 2.5;
alfa_set = .15:.28:1;
rstar_dot_size = 3;

%% .01 conc
ci = 1;
figure
hold on
rasterHigh = 0;
rasteMeadel = [];
[row_slow,col_slow] =...
      find(cell_conc_slow_raster{ci} == 1);
scatter(col_slow, row_slow+rasterHigh, rstar_dot_size,...
    'MarkerEdgeColor', SandR_colors.r_shade(ci+1,:),...
    'MarkerFaceColor', SandR_colors.r_shade(ci+1,:))

rasteMeadel = [rasteMeadel;...
    rasterHigh+...
    (size(cell_conc_slow_raster{ci}, 1)/2)];
rasterHigh = rasterHigh+...
    size(cell_conc_slow_raster{ci}, 1);
plot([0,1200], [rasterHigh, rasterHigh], 'k:',...
    'LineWidth', 1)
    
[row_fast,col_fast] =...
      find(cell_conc_fast_raster{ci} == 1);
scatter(col_fast, row_fast+rasterHigh, rstar_dot_size,...
    'MarkerEdgeColor', SandR_colors.s_shade(ci+1,:),...
    'MarkerFaceColor', SandR_colors.s_shade(ci+1,:)); 
rasteMeadel = [rasteMeadel;...
    rasterHigh+...
    (size(cell_conc_fast_raster{ci}, 1)/2)];
rasterHigh = rasterHigh+...
    size(cell_conc_fast_raster{ci}, 1);

plot([0,1200], [rasterHigh, rasterHigh], 'k:',...
    'LineWidth', 1)

xlim(xlimRange*1000+600)
ylim([0,rasterHigh-1])
yticks(rasteMeadel)
xticks([-.2,0,.2,.4]*1000+600)
xtickangle(45)
ytickangle(45)
yticklabels({'slow','fast'})
xticklabels({'-0.2', '0', '0.2', '0.4'})
title('0.1% conc.')
box off
set(gcf,'Position',[200 700 300 250]);
print_it(pltDir, sprintf('01_conc_raster_ni%d_oId%d', ni, oId),...
    mouseName)
%%--
figure
hold on 
plot(PITH_time, cell_conc_slow_PITH{ci},...
    'Color',SandR_colors.r_shade(ci+1, :),...
    'linewidth',lineWidth);
plot(PITH_time, cell_conc_fast_PITH{ci},...
    'Color',SandR_colors.s_shade(ci+1, :),...
    'linewidth',lineWidth);


p  = patchline([0,0],[-100,100],'linestyle','--','edgecolor','k',...
'linewidth',2,'edgealpha',0.7);
ylim(ylimRange)
xlim(xlimRange)
xticks([-.2,0,.2,.4])
xlabel('Time (sec)');
ylabel('\Delta firing rate (Hz)')
title('0.1% conc.')
box off
set(gcf,'Position',[200 200 300 350]);
print_it(pltDir, sprintf('01_conc_PITH_ni%d_oId%d', ni, oId),...
    mouseName)
%%
markerSize = 25;
transparancyAlfa = 1;
scaterAxisLim = [-2, 60];
dotFaceColor = colorSet(1, :);
dotEdgeColor = colorSet(1, :);
pltTiks = 0:30:60;

close all; 
figure;


scatter(avrage_responce(1), avrage_responce(2),markerSize,...
    'MarkerEdgeColor',dotEdgeColor, 'MarkerFaceColor',dotFaceColor,...
   'MarkerFaceAlpha',transparancyAlfa,'MarkerEdgeAlpha',1)

p  = patchline([.5,100],[.5,100],'linestyle','--','edgecolor','k',...
    'linewidth',2,'edgealpha', 1);
    
ylabel('firing rate (sp./s), S-type')
xlabel('firing rate (sp./s), R-type')
title({'6 mice with two conc',''})
xlim(scaterAxisLim); ylim(scaterAxisLim);
xticks(pltTiks); yticks(pltTiks)
axis square
set(gcf,'Position',[350 350 400 430]);
print_it(pltDir,...
        sprintf('FR_scater_plots_single_point'), 'pooled')
%%
figure;
hold on
plot(PITH_time,...
    neuronsSlow_PITH_noOdor,...
    'Color',  SandR_colors.r, 'LineWidth', lineWidth)
plot(PITH_time,...
    neuronsFast_PITH_noOdor,...
    'Color',  SandR_colors.s, 'LineWidth', lineWidth)
p  = patchline([0,0],[-100,100],'linestyle','--','edgecolor','k',...
'linewidth',2,'edgealpha',0.7);
ylim(ylimRange)
xlim(xlimRange)
xticks([-.2,0,.2,.4])
xlabel('time(sec)');
ylabel('\Delta firing rate (sp./s)')
box off

set(gcf,'Position',[820 200 300 350]);
print_it(pltDir, sprintf('slow_and_fast_PITH_ni%d_no_odor', ni),...
    mouseName)
%---
figure
hold on

rastersDotSize = 2;
[~, Rtype_sortOrder] =...
    sort(respData.inhLen(...
    unitRasterData.raster_labels.labels == 2 &...
    unitRasterData.raster_labels.xpeakCondLabel == 0));
[~, Stype_sortOrder] =...
    sort(respData.inhLen(...
    unitRasterData.raster_labels.labels == 1 &...
    unitRasterData.raster_labels.xpeakCondLabel == 0));

hold on
[rowRtype,colRtype] = find(...
    neuronSlowRaster_noOdor(Rtype_sortOrder,:)==1);
[rowStype,colStype] = find(...
    neuronFastRaster_noOdor(Stype_sortOrder, :)==1);

plot(colRtype, rowRtype+length(Stype_sortOrder),'.', 'Color',...
    SandR_colors.r, 'MarkerSize', rastersDotSize);
plot(colStype, rowStype,'.',...
    'Color', SandR_colors.s,'MarkerSize', rastersDotSize);

p  = patchline([600,600],...
    [0,length(Rtype_sortOrder)+length(Stype_sortOrder)],...
    'linestyle','--','edgecolor','k',...
'linewidth',2,'edgealpha',0.7);
xticks([100,300,500,700]+300)
xticklabels({'-0.2', '0', '0.2', '0.4', '0.6'})
yticks([length(Stype_sortOrder)/2,...
    length(Stype_sortOrder)+length(Rtype_sortOrder)/2])
yticklabels({'S type', 'R type'})
ytickangle(45)
xlim([300, 1000])
ylim([1, length(Rtype_sortOrder)+length(Stype_sortOrder)]) 
xlabel('time (s)')
set(gcf,'Position',[820 640 300 240]);
print_it(pltDir, sprintf('slow_to_fast_Raster_ni%d_no_odor', ni),...
    mouseName)
%% ----
figure;
hold on
plot(PITH_time,...
    neuronsMix_PITH_noOdor,...
    'Color',  [.5,.5,.5], 'LineWidth', lineWidth)
p  = patchline([0,0],[-100,100],'linestyle','--','edgecolor','k',...
'linewidth',2,'edgealpha',0.7);
ylim(ylimRange)
xlim(xlimRange)
xticks([-.2,0,.2,.4])
xlabel('time(sec)');
ylabel('\Delta firing rate (sp./s)')
box off

set(gcf,'Position',[1220 200 300 350]);
print_it(pltDir, sprintf('mix_PITH_ni%d_no_odor', ni),...
    mouseName)
% ---
figure
hold on

rastersDotSize = 2;
[~, sortOrder] =...
    sort(respData.inhLen(...
    unitRasterData.raster_labels.xpeakCondLabel == 0));

hold on
[row,col] = find(...
    neuronFastRaster_Mix(sortOrder,:)==1);

plot(col, row,'.', 'Color',...
    [.5,.5,.5], 'MarkerSize', rastersDotSize);
p  = patchline([600,600],...
    [0,length(sortOrder)],...
    'linestyle','--','edgecolor','k',...
    'linewidth',2,'edgealpha',0.7); 

xticks([100,300,500,700]+300)
xticklabels({'-0.2', '0', '0.2', '0.4', '0.6'})
ytickangle(45)
xlim([300, 1000])
ylim([1, length(sortOrder)])
yticks([])
xlabel('time (s)')
set(gcf,'Position',[1220 640 300 240]);
print_it(pltDir, sprintf('mix_Raster_ni%d_no_odor', ni),...
    mouseName)
%%
