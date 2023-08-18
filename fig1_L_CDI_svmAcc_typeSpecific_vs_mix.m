clear; clc; close all;
[~, mouseNameList, dataPathPrefixList] =....
    get_all_conc_mice('paola');
%%
dataDir = fullfile('G:\My Drive\sinff-odor-project\Figuers_data_Marh2023',...
    'fig1_L_data_auROC_svmAcc_typeSpesefic_vs_mix');

[colorSet, SandR_colors] = set_plot_seting(15, 4);
plotDir = fullfile(final_figs_path('uniTn1'),...
    'fig1_L');
if ~exist(plotDir , 'dir')
   mkdir(plotDir)
end
%%
pooled_Rtype_auROC = [];
pooled_Stype_auROC = [];
pooled_mixType_auROC = [];

pooled_Rtype_SVM_Acc = [];
pooled_Stype_SVM_Acc = [];
pooled_mixType_SVM_Acc = [];
%%
for mn = 1 : length(mouseNameList)
    %%
    mouseName = mouseNameList{mn};
    dataPathPrefix = dataPathPrefixList{mn};
    load(fullfile(dataDir, mouseName))
    
   pooled_Rtype_auROC = [pooled_Rtype_auROC,...
       neuron_discriminability_index_Rtytpe];
   pooled_Stype_auROC = [pooled_Stype_auROC,...
       neuron_discriminability_index_Stytpe];
   pooled_mixType_auROC = [pooled_mixType_auROC,...
       neuron_discriminability_index_MixTytpe];

   pooled_Rtype_SVM_Acc = [pooled_Rtype_SVM_Acc,...
       neuron_svm_acc_index_Rtytpe];
   pooled_Stype_SVM_Acc = [pooled_Stype_SVM_Acc,...
       neuron_svm_acc_index_Stytpe];
   pooled_mixType_SVM_Acc = [pooled_mixType_SVM_Acc,...
       neuron_svm_acc_index_MixTytpe];
end
%%
close all
subplot(1,3,1)
hold on
h = cdfplot(pooled_Rtype_auROC);
set(h, 'Color', SandR_colors.r);
h = cdfplot(pooled_mixType_auROC);
set(h, 'Color', colorSet(1, :));
ylabel('CDF', 'FontAngle', 'italic')
box off
title([])

subplot(1,3,2)
hold on
h = cdfplot(pooled_Stype_auROC);
set(h, 'Color', SandR_colors.s);
h = cdfplot(pooled_mixType_auROC);
set(h, 'Color', colorSet(1, :));
xlabel({'CDI'})
box off
title([])

subplot(1,3,3)
hold on
h = cdfplot(pooled_Rtype_auROC);
set(h, 'Color', SandR_colors.r);
h = cdfplot(pooled_Stype_auROC);
set(h, 'Color', SandR_colors.s);
h = cdfplot(pooled_mixType_auROC);
set(h, 'Color', colorSet(1, :));
legend('R type inh. type','S type inh. type', 'mix inh. type',...
    'Location', 'best')
box off
legend boxoff
title([])

sgtitle({'10-fold conc discrimination, all cell-odor pairs'},...
    'fontsize', 20)
set(gcf,'Position',[350 350 1000 350]);
print_it(plotDir, 'discriminability_index_CDF', 'pooled')
%%
%close all;
figure
hold on
bh = boxplot([pooled_mixType_auROC',...
    pooled_Rtype_auROC',...
    pooled_Stype_auROC'],...
    'Labels',{'mix type', 'R type', 'S type'},...
    'Colors',[colorSet(1,:); SandR_colors.r; SandR_colors.s]);
set(bh,'LineWidth', 2);
ylabel({'CDI', ''})
xtickangle(50)
xlim([.5, 3.5])
ylim([-.05,.8])
box off
title({'10-fold conc' 'discrimination'})
set(gcf,'Position',[350 350 330 350]);
print_it(plotDir, 'discriminability_index_boxPlots', 'pooled')
%%
pRankSumRtype =...
    ranksum(pooled_Rtype_auROC ,...
    pooled_mixType_auROC );
pRankSumStype =...
    ranksum(pooled_Stype_auROC ,...
    pooled_mixType_auROC );
[~,pKs_mix_R] = kstest2(pooled_Rtype_auROC ,...
    pooled_mixType_auROC );
[~,pKs_mix_S] = kstest2(pooled_Stype_auROC ,...
    pooled_mixType_auROC );
%% 
[~,~,stats] = kruskalwallis(...
    [pooled_mixType_auROC';...
    pooled_Rtype_auROC';...
    pooled_Stype_auROC'],...
    [ones(length(pooled_mixType_auROC), 1);...
    2*ones(length(pooled_Rtype_auROC), 1);...
    3*ones(length(pooled_Stype_auROC), 1)]);
c = multcompare(stats);
%%
fileId = fopen(fullfile(plotDir, ...
    'CDI_typeSpesefic_vs_mix_num_stat.txt'), 'w');
fprintf(fileId, 'mix v.s. R-type CDF unequality Kolmogorov-Smirnov test p val = %.5f \n',...
    pKs_mix_R);
fprintf(fileId, 'mix v.s. R-type CDF unequality Kolmogorov-Smirnov test p val = %.5f \n',...
    pKs_mix_S);

fprintf(fileId, 'mix v.s. R-type ranksum test for discriminability index difference = %.5f. \n',...
    pRankSumRtype);
fprintf(fileId, 'mix v.s. S-type ranksum test for discriminability index difference = %.5f. \n',...
    pRankSumStype);

fprintf(fileId, 'Discriminability index R type = %.2f +- %.2f. \n',...
    mean(pooled_Rtype_auROC ),...
    std(pooled_Rtype_auROC )/...
    sqrt(length(pooled_Rtype_auROC )));
fprintf(fileId, 'Discriminability index S type = %.2f +- %.2f. \n',...
    mean(pooled_Stype_auROC ),...
    std(pooled_Stype_auROC )/...
    sqrt(length(pooled_Stype_auROC )));
fprintf(fileId, 'Discriminability index mix type = %.2f +- %.2f. \n',...
    mean(pooled_mixType_auROC ),...
    std(pooled_mixType_auROC )/...
    sqrt(length(pooled_mixType_auROC )));
fclose(fileId);
%%
typeSpecificAverage_auROC = ...
    mean([pooled_Rtype_auROC; pooled_Stype_auROC]);

close all
hold on
h = cdfplot(typeSpecificAverage_auROC);
set(h, 'Color', colorSet(2, :));
h = cdfplot(pooled_mixType_auROC);
set(h, 'Color', colorSet(1, :));

legend('R-S type inh. avrage', 'mix inh. type',...
    'Location', 'bestoutside')
ylabel('CDF', 'FontAngle', 'italic')
xlabel({'CDI'})
title({'10-fold conc discrimination','all cell-odor pairs'})
box off
legend boxoff
set(gcf,'Position',[350 350 600 300]);
print_it(plotDir, 'discriminability_index_CDF_avrageTypes_vs_mix', 'pooled')
%%
%close all;
figure
hold on
bh = boxplot([pooled_mixType_auROC',...
    typeSpecificAverage_auROC'],...
    'Labels',{'mix type', 'R-S type avg.'},...
    'Colors',[colorSet(1,:); colorSet(2,:)]);
set(bh,'LineWidth', 2);
ylabel({'CDI', ''})
xtickangle(50)
xlim([.5, 2.5])
ylim([-.05,.8])
box off
title({'10-fold conc' 'discrimination'})
set(gcf,'Position',[350 350 300 350]);
print_it(plotDir, 'discriminability_index_boxPlots_avrageTypes_vs_mix', 'pooled')
%%
axisLim = [-.08,.8];
yLim = [-4,4];


close all
raincloud_plot(pooled_mixType_auROC,...
    'color', colorSet(2, :), 'box_on', 1, 'alpha', .5,...
    'lwr_bnd', 1, 'line_width', 1.5);
box off
xlabel('CDI')
title({'mix type',''})
xlim(axisLim);
ylim(yLim)
set(gcf,'Position',[350 350 500 400])
print_it(plotDir, 'CDI_raincloud_vs_mix', 'pooled')
%%---
close all
raincloud_plot(typeSpecificAverage_auROC,...
    'color', colorSet(1, :), 'box_on', 1, 'alpha', .5,...
    'lwr_bnd', 1, 'line_width', 1.5);
box off
xlabel('CDI')
title({'R-S type avg.',''})
xlim(axisLim);
ylim(yLim)
set(gcf,'Position',[350 350 500 400])
print_it(plotDir, 'CDI_raincloud_vs_avrageTypes', 'pooled')
%%
pRanksum_avg_mix =...
    ranksum(typeSpecificAverage_auROC ,...
    pooled_mixType_auROC );
[~,pKs_avg_mix] = kstest2(typeSpecificAverage_auROC ,...
    pooled_mixType_auROC );

fileId = fopen(fullfile(plotDir, ...
    'auROC_typesAvrage_vs_mix_num_stat.txt'), 'w');
fprintf(fileId, 'mix v.s. R-S-type avrage CDF unequality Kolmogorov-Smirnov test p val = %.5f \n',...
    pRanksum_avg_mix);

fprintf(fileId, 'mix v.s. R-S-type avrage ranksum test for discriminability index difference = %.5f. \n',...
    pKs_avg_mix);

fprintf(fileId, 'Discriminability index R-S-type avrage type = %.2f +- %.2f. \n',...
    mean(typeSpecificAverage_auROC ),...
    std(typeSpecificAverage_auROC )/...
    sqrt(length(typeSpecificAverage_auROC )));
fprintf(fileId, 'Discriminability index mix type = %.2f +- %.2f. \n',...
    mean(pooled_mixType_auROC ),...
    std(pooled_mixType_auROC )/...
    sqrt(length(pooled_mixType_auROC )));
fclose(fileId);
%% CVM resulta

close all
subplot(1,3,1)
hold on
h = cdfplot(pooled_Rtype_SVM_Acc);
set(h, 'Color', SandR_colors.r);
h = cdfplot(pooled_mixType_SVM_Acc);
set(h, 'Color', colorSet(1, :));
ylabel('CDF', 'FontAngle', 'italic')
box off
title([])

subplot(1,3,2)
hold on
h = cdfplot(pooled_Stype_SVM_Acc);
set(h, 'Color', SandR_colors.s);
h = cdfplot(pooled_mixType_SVM_Acc);
set(h, 'Color', colorSet(1, :));
xlabel({'accuracy (%)'})
box off
title([])

subplot(1,3,3)
hold on
h = cdfplot(pooled_Rtype_SVM_Acc);
set(h, 'Color', SandR_colors.r);
h = cdfplot(pooled_Stype_SVM_Acc);
set(h, 'Color', SandR_colors.s);
h = cdfplot(pooled_mixType_SVM_Acc);
set(h, 'Color', colorSet(1, :));
legend('R type inh. type','S type inh. type', 'mix inh. type',...
    'Location', 'best')
box off
legend boxoff
title([])

sgtitle({'10-fold conc decoding, all cell-odor pairs'},...
    'fontsize', 20)
set(gcf,'Position',[350 350 1000 350]);
print_it(plotDir, 'SVM_accuracy_CDF', 'pooled')
%%
%close all;
figure
hold on
plot([0,4], [50,50], 'k--')
bh = boxplot([pooled_mixType_SVM_Acc',...
    pooled_Rtype_SVM_Acc',...
    pooled_Stype_SVM_Acc'],...
    'Labels',{'mix type', 'R type', 'S type'},...
    'Colors',[colorSet(1,:); SandR_colors.r; SandR_colors.s]);

set(bh,'LineWidth', 2);
ylabel({'accuracy (%)', ''})
xtickangle(50)
xlim([.5, 3.5])
ylim([47,70])
box off
title({'10-fold conc' 'decoding'})
set(gcf,'Position',[350 350 330 400]);
print_it(plotDir, 'SVM_accuracy_boxPlots', 'pooled')
%%
pRankSumRtype =...
    ranksum(pooled_Rtype_SVM_Acc ,...
    pooled_mixType_SVM_Acc );
pRankSumStype =...
    ranksum(pooled_Stype_SVM_Acc ,...
    pooled_mixType_SVM_Acc );
[~,pKs_mix_R] = kstest2(pooled_Rtype_SVM_Acc ,...
    pooled_mixType_SVM_Acc );
[~,pKs_mix_S] = kstest2(pooled_Stype_SVM_Acc ,...
    pooled_mixType_SVM_Acc );
%% 
[p,tbl,stats] = kruskalwallis(...
    [pooled_mixType_SVM_Acc';...
    pooled_Rtype_SVM_Acc';...
    pooled_Stype_SVM_Acc'],...
    [ones(length(pooled_mixType_SVM_Acc), 1);...
    2*ones(length(pooled_Rtype_SVM_Acc), 1);...
    3*ones(length(pooled_Stype_SVM_Acc), 1)]);
c = multcompare(stats);
%%
fileId = fopen(fullfile(plotDir, ...
    'SVM_accuracy_typeSpesefic_vs_mix_num_stat.txt'), 'w');
fprintf(fileId, 'mix v.s. R-type SVM accuracy CDF unequality Kolmogorov-Smirnov test p val = %.5f \n',...
    pKs_mix_R);
fprintf(fileId, 'mix v.s. R-type SVM accuracy CDF unequality Kolmogorov-Smirnov test p val = %.5f \n',...
    pKs_mix_S);

fprintf(fileId, 'mix v.s. R-type ranksum test for SVM accuracy difference = %.5f. \n',...
    pRankSumRtype);
fprintf(fileId, 'mix v.s. S-type ranksum test for SVM accuracy difference = %.5f. \n',...
    pRankSumStype);

fprintf(fileId, 'SVM accuracy R type = %.2f +- %.2f. \n',...
    mean(pooled_Rtype_SVM_Acc ),...
    std(pooled_Rtype_SVM_Acc )/...
    sqrt(length(pooled_Rtype_SVM_Acc )));
fprintf(fileId, 'SVM accuracy S type = %.2f +- %.2f. \n',...
    mean(pooled_Stype_SVM_Acc ),...
    std(pooled_Stype_SVM_Acc )/...
    sqrt(length(pooled_Stype_SVM_Acc )));
fprintf(fileId, 'SVM accuracy mix type = %.2f +- %.2f. \n',...
    mean(pooled_mixType_SVM_Acc ),...
    std(pooled_mixType_SVM_Acc )/...
    sqrt(length(pooled_mixType_SVM_Acc )));
fclose(fileId);
%%
typeSpecificAverage_svm = ...
    mean([pooled_Rtype_SVM_Acc; pooled_Stype_SVM_Acc]);

close all
hold on
h = cdfplot(typeSpecificAverage_svm);
set(h, 'Color', colorSet(2, :));
h = cdfplot(pooled_mixType_SVM_Acc);
set(h, 'Color', colorSet(1, :));

legend('R-S type inh. avrage', 'mix inh. type',...
    'Location', 'bestoutside')
ylabel('CDF', 'FontAngle', 'italic')
xlabel({'accuracy (%)'})
title({'10-fold conc decoding','all cell-odor pairs'})
box off
legend boxoff
set(gcf,'Position',[350 350 600 300]);
print_it(plotDir, 'SVM_accuracy_CDF_avrageTypes_vs_mix', 'pooled')
%%
%close all;
figure
hold on
plot([0,4], [50,50], 'k--')
bh = boxplot([pooled_mixType_SVM_Acc',...
    typeSpecificAverage_svm'],...
    'Labels',{'mix type', 'R-S type avg.'},...
    'Colors',[colorSet(1,:); colorSet(2,:)]);
set(bh,'LineWidth', 2);
ylabel({'accuracy (%)', ''})
xtickangle(50)
xlim([.5, 2.5])
ylim([47,70])
box off
title({'10-fold conc' 'decoding'})
set(gcf,'Position',[350 350 300 400]);
print_it(plotDir, 'SVM_accuracy_boxPlots_avrageTypes_vs_mix', 'pooled')
%% rainclouad plots 
axisLim = [45,70];
yLim = [-.2,.2];
%%--
close all
raincloud_plot(pooled_mixType_SVM_Acc,...
    'color', colorSet(2, :), 'box_on', 1, 'alpha', .5,...
    'lwr_bnd', 1, 'line_width', 1.5);
box off
xlabel('accuracy (%)')
title({'mix type',''})
xlim(axisLim);
ylim(yLim)
set(gcf,'Position',[350 350 500 400])
print_it(plotDir, 'SVM_accuracy_raincloud_vs_mix', 'pooled')
%%--
close all
raincloud_plot(typeSpecificAverage_svm,...
    'color', colorSet(1, :), 'box_on', 1, 'alpha', .5,...
    'lwr_bnd', 1, 'line_width', 1.5);
box off
xlabel('accuracy (%)')
title({'R-S type avg.',''})
xlim(axisLim);
ylim(yLim)
set(gcf,'Position',[350 350 500 400])
print_it(plotDir, 'SVM_accuracy_raincloud_vs_avrageTypes', 'pooled')
%%
pRanksum_avg_mix =...
    ranksum(typeSpecificAverage_svm ,...
    pooled_mixType_SVM_Acc );
[~,pKs_avg_mix] = kstest2(typeSpecificAverage_svm ,...
    pooled_mixType_SVM_Acc );

fileId = fopen(fullfile(plotDir, ...
    'SVM_accuracy_typesAvrage_vs_mix_num_stat.txt'), 'w');
fprintf(fileId, 'mix v.s. R-S-type SVM accuracy avrage CDF unequality Kolmogorov-Smirnov test p val = %.5f \n',...
    pRanksum_avg_mix);

fprintf(fileId, 'mix v.s. R-S-type avrage ranksum test for SVM accuracy difference = %.5f. \n',...
    pKs_avg_mix);

fprintf(fileId, 'SVM accuracy R-S-type avrage type = %.2f +- %.2f. \n',...
    mean(typeSpecificAverage_auROC ),...
    std(typeSpecificAverage_auROC )/...
    sqrt(length(typeSpecificAverage_auROC )));
fprintf(fileId, 'SVM accuracy mix type = %.2f +- %.2f. \n',...
    mean(pooled_mixType_SVM_Acc ),...
    std(pooled_mixType_SVM_Acc )/...
    sqrt(length(pooled_mixType_SVM_Acc )));
fclose(fileId);