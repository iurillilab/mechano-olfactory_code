clearvars -except T B R tests M
folderDataSave = '/Users/galileo/Dropbox (Personal)/SniffOdor_project_processedData/sniffOdorProject/Bulbectomy/processedDataStorage';
generalOutputDir = fullfile('supp_bulbectomized_tuning_newCode', 'firstSniff', 'across_comparisons');
generalPlotDir = fullfile(final_figs_path('paola'), generalOutputDir, 'Aug14_0to180vsSLOW_ttest_RightTail_p001', 'rightTail_test3');
mkdir(generalPlotDir)
cd(generalPlotDir)
testNuse = 2;

%% load data
if ~exist('B', 'var')
    load('/Users/galileo/Dropbox (Personal)/SniffOdor_project_processedData/sniffOdorProject/Bulbectomy/processedDataStorage/analysisTables.mat', 'T', 'B', 'R');
end
if exist(fullfile(folderDataSave, 'firstSniffTest_metadata.mat'), 'file')
    tests = load(fullfile(folderDataSave, 'firstSniffTest_metadata.mat'));
else
    tests = [];
end
load('/Users/galileo/Dropbox (Personal)/SniffOdor_project_processedData/sniffOdorProject/Bulbectomy/processedDataStorage/sigOdorResponse_bulbect.mat')

fnames = {'fs', 'lr', 'slr', 'sniffs', 'rRr'};
brforeInhOnsetWind = .6;
afterInhOnsetWind = .6;
PITH_time = -1*brforeInhOnsetWind:.001 : afterInhOnsetWind +.001; 

% plot settings
YLIM = [-4 20];
XLIM = [-0.3 0.4];
fastCol = [136 87 164]./255;
regCol = [0 166 81]./255;
alphashade = 0.2;
set(0,'DefaultFigureWindowStyle', 'docked')
close all
figure(101)


%% subset and plot timecourses
YLIM = [-5 18];
for ei = 1:3
    expName = T.exp{find(T.expID==ei,1)};
    test = tests.(expName);


usefield_p = test(testNuse).usefield_p;         
usefield_resp = test(testNuse).usefield_resp;   

useT = B(:, 1:2);
% useT.p = table2array(B(:, ismember(B.Properties.VariableNames, usefield_p)));
useT.dR = table2array(B(:, ismember(B.Properties.VariableNames, usefield_resp)));
useT.sig  = test(testNuse).cellsIncluded;


subT = T(useT.sig==1 & useT.expID == ei, :);
subT.dr = useT.dR(useT.sig==1 & useT.expID == ei);
subT = sortrows(subT, {'dr'}, 'descend');

fs = cat(1,subT.firstsniff_PITH{:});
alls = cat(1,subT.sniff_PITH{:});
rRr = cat(1,subT.rRr_PITH{:});
slr = cat(1,subT.seclastreg_PITH{:});
lr = cat(1,subT.lastreg_PITH{:});

% % plot grand-average ~
PCTLremove = 1;

%all cells (minus some too noisy to be zscored)
% zscore over baseline
preW = PITH_time <= -0.2;
baseline_mean = mean(fs(:,preW), 2);
baseline_std = std(fs(:,preW), [], 2);
keepcells = baseline_std>=prctile(baseline_std,PCTLremove); %just remove extreme noise
zs_fs = bsxfun(@minus, fs(keepcells,:), baseline_mean(keepcells));
zs_fs = bsxfun(@rdivide, zs_fs, baseline_std(keepcells));
sum(keepcells)

baseline_mean = mean(rRr(:,preW), 2);
baseline_std = std(rRr(:,preW), [], 2);
keepcells = baseline_std>=prctile(baseline_std,PCTLremove); %just remove extreme noise
zs_slow = bsxfun(@minus, rRr(keepcells,:), baseline_mean(keepcells));
zs_slow = bsxfun(@rdivide, zs_slow, baseline_std(keepcells));


figure; hold on
clear h  %first sniff
stdshade_modified(zs_fs, alphashade, fastCol ,PITH_time,[],[], 2)   
stdshade_modified(zs_slow, alphashade, regCol ,PITH_time,[],[], 2)  
% plot([0 0], YLIM, 'LineWidth', 0.1, 'Color', [0.5 0.5 0.5])
h(2) = plot(PITH_time, mean(zs_slow), 'LineWidth',2, 'Color', regCol);
h(1) = plot(PITH_time, mean(zs_fs), 'LineWidth',2, 'Color', fastCol);
ylabel(sprintf('zscored on baseline (N = %d/%d)', sum(keepcells), length(keepcells)))
title(sprintf('%s: all sig %s cells', expName, usefield_p), 'Interpreter', 'none')
legend(h, 'first sniff', 'rRr')
legend boxoff
box off
ylim(YLIM)
xlim(XLIM)

print_it(generalPlotDir, sprintf('grandAverage_zscONBas_firstSniff_rRs_', usefield_p), expName)
% close

end


%% plot distribution of amplitude of responses d(fast_post - slow_post)
cab(101)
dGrater = [];
groupGrater = [];
g = [];
for ei = 1:3
    expName = T.exp{find(T.expID==ei,1)};
    test = tests.(expName);
    usefield_p = test(testNuse).usefield_p;         %'genRT_fs' ;%'rfpR'; %R is right tail
    usefield_resp = test(testNuse).usefield_resp;   %'gendResp_fs';%'rfresp';

    useT = R(:, 1:3);
    useT.dR = table2array(B(:, ismember(B.Properties.VariableNames, usefield_resp)));
    
    cellsIncluded = test(testNuse).cellsIncluded;
    subT = useT(cellsIncluded,:);
    %     subT = sortrows(subT, {'dr'}, 'descend');
    
    dGrater = cat(1, dGrater, subT.dR);
    groupGrater = cat(1, groupGrater, subT.exp);
    g = cat(1, g, ei*ones(sum(cellsIncluded),1));
end

% box plot + ANOVA - right tail
figure; hold on
boxplot(dGrater, groupGrater, 'Colors', 'k')
ylim([0, 60])
title(sprintf('%s:\namplitude of fast-slow in 0-180 ms\nfor the first sniff in responsive neurons',usefield_p))
box off
ax = gca;
delete(ax.Children.Children(1:3))

% add single cell dResps
for ei = 1:3
    expName = T.exp{find(T.expID==ei,1)};
    test = tests.(expName);
    usefield_p = test(testNuse).usefield_p;         %'genRT_fs' ;%'rfpR'; %R is right tail
    usefield_resp = test(testNuse).usefield_resp;   %'gendResp_fs';%'rfresp';
    
    resp2plot = table2array(B(test(testNuse).cellsIncluded, ismember(B.Properties.VariableNames, usefield_resp)));
    
    X = ei + (rand(sum(test(testNuse).cellsIncluded),1)-0.5)*0.4;
    scatter(X,resp2plot, 50, 'k', 'filled')   
end


print_it(generalPlotDir, sprintf('boxplot_amplitude_firstSniff_%s', usefield_p), 'rightTail')


diary boxplot_amplitude_rightTailsComparison.txt
diary on
disp('one-way ANOVA on fast-preferring units across regions:')
[p,tbl,stats] = anova1(dGrater, groupGrater);
disp('stats:')
disp(stats)
disp('group names:')
disp(stats.gnames)
disp('ANOVA:')
disp(tbl)
figure;
disp('multiple comparison corrected post-hoc test (tukey-kramer):')
c = multcompare(stats, 'CType','tukey-kramer');
disp(c)

diary off
diary boxplot_amplitude_rightTailsComparison.txt
diary on
disp('two-sample ttest betw pirif and bulbect only:')
apirif = dGrater(g==1);
abulb = dGrater(g==3);
[h,p] = vartest2(apirif, abulb);
if h
    [h,p] = ttest2(apirif, abulb, 'Vartype','unequal');
    disp('two-sided test, unequal variances:')
else
    [h,p] = ttest2(apirif, abulb, 'Vartype','equal');
    disp('two-sided test, equal variances:')
end
diary off
%% NEW boxplot of amplitude with single cell lines
respW = PITH_time > 0 & PITH_time <= 0.180;
all2plot = [];
group2plot = [];
for ei = 1:3
    expName = T.exp{find(T.expID==ei,1)};
    test = tests.(expName);
    usefield_p = test(testNuse).usefield_p;         %'genRT_fs' ;%'rfpR'; %R is right tail
    cellsIncluded = test(testNuse).cellsIncluded;

   
    pithS = cat(1, T.seclastreg_PITH{cellsIncluded});
    pithS = mean(pithS(:,respW),2);
    
    pithF = cat(1, T.firstsniff_PITH{cellsIncluded});
    pithF = mean(pithF(:,respW),2);
    
    
    all2plot = cat(1, all2plot, (cat(2, pithS, pithF)));
    group2plot = cat(1, group2plot, ei*ones(sum(cellsIncluded),1));
end

% figure; 
hold on
for ei = 1:3 % add single cell lines
    expName = T.exp{find(T.expID==ei,1)};
    xs = ei+[-0.35,0.35];
    plot(xs, all2plot(group2plot==ei,:)', 'Color', 0.4*ones(1,3), 'LineWidth', 0.25)
end
ylim([-20 60])
print_it(generalPlotDir, sprintf('boxplot_amplitude_firstSniff_%s_withCellLines', usefield_p), 'rightTail')

% %% repeat mixed model anova as previously done (do with and without motor cortex)
% respW = PITH_time > 0 & PITH_time <= 0.180;
% condUseTags = {'seclastreg_PITH', 'firstsniff_PITH'};
% tagsused = {};
% all2plot = [];
% groupNames = {};
% g = [];
% count = 0;
% for ei = [1,3]
%     expName = T.exp{find(T.expID==ei,1)};
%     test = tests.(expName);
%     usefield_p = test(testNuse).usefield_p;         %'genRT_fs' ;%'rfpR'; %R is right tail
%     cellsIncluded = test(testNuse).cellsIncluded;
%     groupNames = cat(1, groupNames, T.exp(cellsIncluded));  %single repet (one per cell)
%     for r = 1: length(condUseTags)  % in here two repets per cell
%         count = count+1;
%         pith = cat(1, T.(condUseTags{r}){cellsIncluded});
%         pith = mean(pith(:,respW),2);
%         
%         all2plot = cat(1, all2plot, pith);  
%         g = cat(1, g, count*ones(sum(cellsIncluded),1));
%     end
%     tagsused = cat(1,tagsused, expName);
% end
% 
% figure; 
% hold on
% for ei = [1,3] % add single cell lines
%     expName = T.exp{find(T.expID==ei,1)};
%     xs = ei*2+[-1,0];
%     plot(xs, cat(2, all2plot(g==xs(1)), all2plot(g==xs(2)))', 'Color', 0.4*ones(1,3), 'LineWidth', 0.25)
% end
% ylim([-20 60])
% boxplot(all2plot, g, 'Colors', 'k')
% a = gca;
% a.XTick = (2:2:6)-0.5;
% a.XTickLabel = tagsused;
% 
% print_it(generalPlotDir, sprintf('boxplot_amplitude_firstSniff_SLOWFAST_SLPIT_%s_withCellLines', usefield_p), 'rightTail')
% 
% % ANOVA between fast amplitudes
% %
% clc
% diary boxplot_fastSniffAmplitude_comparison.txt
% diary on
% disp('one-way ANOVA on amplitude of first sniff in sig cells across regions:')
% [p,tbl,stats] = anova1(all2plot(ismember(g,2:2:6)), groupNames);
% stats
% stats.gnames
% tbl
% figure;
% disp('multiple comparison corrected post-hoc test (tukey-kramer):')
% c = multcompare(stats, 'CType','tukey-kramer')
% 
% diary off



% % % mixed model anova - DO NOT USE BC within subject significance is how
% these cells were defined in the first place. no need to recalculate ttest

% betweenSbjFactor = nan(size(g));
% betweenSbjFactor(g==1 | g==2) = 1;
% betweenSbjFactor(g==3 | g==4) = 2;
% betweenSbjFactor(g==5 | g==6) = 3;
% 
% withinSbjFactor = nan(size(g));
% withinSbjFactor(mod(g,2)==1) = 1;
% withinSbjFactor(mod(g,2)==0) = 2;
% % withinSbjFactor = g;
% 
% subjectcodes = [];
% subjectcodes = cat(1, subjectcodes, repmat((1:sum(g==1))',2,1));
% subjectcodes = cat(1, subjectcodes, repmat(sum(g==1)+(1:sum(g==3))',2,1));
% subjectcodes = cat(1, subjectcodes, repmat(sum(g==1)+sum(g==3)+(1:sum(g==5))',2,1));
% 
% clear X
% X(:,1) = all2plot;
% X(:,2) = betweenSbjFactor;
% X(:,3) = withinSbjFactor;
% X(:,4) = subjectcodes;
% 
% [SSQs, DFs, MSQs, Fs, Ps]=mixed_between_within_anova(X);

%% box plot of distribution of responsive cells by region and animals
%find proportion of responsive cells, mouse by mouse
groupedProportions = [];
groupNames = [];
groupID = [];
for ei = 1:3
    expName = T.exp{find(T.expID==ei,1)};
    useT = R(:, 1:3);    
    useT.cellsIncluded = tests.(expName)(testNuse).cellsIncluded;
    useT = useT(useT.expID==ei,:);

    sigCellMice = useT.mouseOrdinal.*useT.cellsIncluded;
    sigCellMice(sigCellMice==0) = [];
    exp(ei).numer = [];
    exp(ei).denom = [];
    for m = 1:max(useT.mouseOrdinal)
        exp(ei).numer(m) = sum(sigCellMice==m);
        exp(ei).denom(m) = sum(useT.mouseOrdinal==m);
        pr = exp(ei).numer(m)/exp(ei).denom(m);
        groupedProportions = cat(1, groupedProportions, pr);
        groupNames = cat(1, groupNames, {expName});
        groupID = cat(1,groupID, ei);
    end
    figure
    histogram(sigCellMice)
    title(expName)
    ylabel(sprintf('significantly modulated to test %d', testNuse))
    xlabel('mouse ordinal')
    box off
    print_it(generalPlotDir, sprintf('boxplot_firstSniffResponsiveCells_%s', usefield_p), expName)

end


figure; hold on
boxplot(groupedProportions, groupNames, 'Colors', 'k')
box off
% add single mice data
for ei = 1:3
    Y = groupedProportions(groupID==ei);
    X = ei + (rand(sum(groupID==ei),1)-0.5)*0.3;
    scatter(X,Y, 100, 'k', 'filled')   
end

title(sprintf('%s:\nfraction of neurons responsive to first sniff',usefield_p))
ax = gca;
print_it(generalPlotDir, sprintf('boxplot_firstSniffResponsiveFraction_%s', usefield_p), 'rightTail')


%
clc
diary boxplot_firstSniffResponsiveFraction_comparison.txt
diary on
disp('one-way ANOVA on fast-preferring units across regions:')
[p,tbl,stats] = anova1(groupedProportions, groupNames);
stats
stats.gnames
tbl
figure;
disp('multiple comparison corrected post-hoc test (tukey-kramer):')
c = multcompare(stats, 'CType','tukey-kramer')

diary off

close all


%% calculate coefficient of variation of all responsive neurons
% get the raster, 
% take the postW average (no subtractions) for each inhalation
% std/mean
usefield_p = tests.bulbect(testNuse).usefield_p;
usefield_resp = tests.bulbect(testNuse).usefield_resp;

postW = PITH_time >= 0 & PITH_time < 0.18;

useSig = false(size(T,1),1);
for ei = 1:3
    expName = T.exp{find(T.expID==ei,1)};
    useSig = useSig | tests.(expName)(testNuse).cellsIncluded;
end
useT = R(:, 1:3);  
% get all cells included:
% p = table2array(B(:, ismember(B.Properties.VariableNames, usefield_p)));
% dR = table2array(B(:, ismember(B.Properties.VariableNames, usefield_resp)));
useT.sig  = useSig;
useT.firstsniff_Raster = R.firstsniff_Raster;
useT = useT(useT.sig==1,:);

trialNumbers = [];
for i = 1:size(useT,1)
    trialNumbers(i) = size(useT.firstsniff_Raster{i},1);
end
useT.trialNumbers = trialNumbers(:);
groupedTrialNumber = [];
for ei = 1:3
    mIDS = unique(useT.mouseOrdinal(useT.expID==ei), 'stable')';
    length(mIDS)
    for m = 1:length(mIDS)
        groupedTrialNumber = cat(1, groupedTrialNumber, unique(useT.trialNumbers(useT.expID==ei & useT.mouseOrdinal==mIDS(m))) );
    end
end

CV = [];
for i = 1:size(useT,1)
    raster = useT.firstsniff_Raster{i};
    raster = raster(:,postW);
    meanresp = mean(raster,2);
    rmean = mean(meanresp);
    rstd = std(meanresp);
    cv = rstd/rmean;
    CV = cat(1, CV, cv);
end
useT.CV = CV(:);

figure; hold on
boxplot(useT.CV, useT.exp, 'Colors', 'k')
title('coefficient of variation, resp to first sniff')
ax = gca;
delete(ax.Children.Children(1:3))

% add single neuron data, grouped by mouse
colors = distinguishable_colors(10);
colors(2,3) = 0.1;
mousegroup = useT.expID*10 + useT.mouseOrdinal; %not used
for ei = 1:3
    mIDS = unique(useT.mouseOrdinal(useT.expID==ei), 'stable')';
    for m = 1:length(mIDS)
        Y = useT.CV(useT.expID==ei & useT.mouseOrdinal==mIDS(m));
        X = ei + (rand(length(Y),1)-0.5)*0.4;
        scatter(X,Y, 50, 'k', 'filled') 
    end
end
box off

print_it(generalPlotDir, sprintf('boxplot_CV_%s', usefield_p), 'rightTail')

%
clc
diary boxplot_CVcomparison.txt
diary on
disp('one-way ANOVA on fast-preferring units across regions:')
[p,tbl,stats] = anova1(useT.CV, useT.exp);
stats
stats.gnames
tbl
figure;
disp('multiple comparison corrected post-hoc test (tukey-kramer):')
c = multcompare(stats, 'CType','tukey-kramer')

diary off



diary boxplot_CVcomparison.txt
diary on
disp('two-sample ttest betw pirif and bulbect CVs only:')
apirif = CV(useT.expID==1);
abulb = CV(useT.expID==3);
[h,p] = vartest2(apirif, abulb);
if h
    [h,p] = ttest2(apirif, abulb, 'Vartype','unequal');
    disp('two-sided test, unequal variances:')
else
    [h,p] = ttest2(apirif, abulb, 'Vartype','equal');
    disp('two-sided test, equal variances:')
end
h
p
diary off

%% latence vs rRr
binSize  = 10; %ms

t = PITH_time(1:end-2);
tb = reshape(t,binSize,[]);
tcenters = mean(tb(floor(binSize/2):ceil(binSize/2),:), 1);

p = nan(length(tcenters),3);
for ei = 1:3
    binnedF = [];
    binnedS = [];
    expName = T.exp{find(T.expID==ei,1)};
    useT = R(:, 1:3);    
    useT.cellsIncluded = tests.(expName)(testNuse).cellsIncluded;
    useT.firstsniff_PITH = T.firstsniff_PITH;                                %T.lastreg_PITH;%to compare
    useT.rRr_PITH = T.rRr_PITH;    
    useT.cellsIncluded = tests.(expName)(testNuse).cellsIncluded;
    cellsIncluded = useT.cellsIncluded;
    
    useT = useT(cellsIncluded ,:);
    disp(size(useT))
    
    
    for i = 1:size(useT,1)
        F = useT.firstsniff_PITH{i};
        S = useT.rRr_PITH{i};
        F = F(:,1:end-2);
        S = S(:,1:end-2);
        
        F = mean(reshape(F,binSize,[]), 1);
        S = mean(reshape(S,binSize,[]), 1);
        binnedF = cat(1, binnedF, F);
        binnedS = cat(1, binnedS, S);
    end
    
    for nb = 1:length(tcenters)
        [~, p(nb,ei)] = ttest(binnedF(:,nb), binnedS(:,nb), 'tail', 'right');
    end
    
end
clear c_pvalues c_alpha h
for i = 1:3
    [c_pvalues(:,i), c_alpha(:,i), h(i,:), ~] = fdr_BH(p(:,i), 0.05);
end

figure; hold on
plot(tcenters,c_pvalues)
ax= gca;
ax.YAxis.Scale = 'log';
legend('control', 'motor', 'bulbect')

% plot(tcenters,p(:,3))
% ax= gca;
% ax.YAxis.Scale = 'log';
% legend('bulbect')

hold on
for i = 1:3
    plot(tcenters, 10^(-20-i)*(c_pvalues(:,i)<0.05), 'LineWidth',2);
end
xlim([-0.3 0.4])


legend boxoff
title(sprintf('firstsniff vs rRr - binsize: %d - ttest',binSize))
xlabel('time')
ylabel('delta (first Sniff - rRr) FDR-BH-corrected p')
grid on

plot([tcenters(1), tcenters(end)], [0.05 0.05], '-k')
print_it(generalPlotDir, sprintf('p_overTime_firstSniff_vs_rRr_%s', usefield_p), 'rightTail')



