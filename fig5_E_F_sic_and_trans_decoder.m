clear; clc; close all;

dataDir = fullfile('C:\Users\a.abolghasemi',...
    'OneDrive - Fondazione Istituto Italiano Tecnologia',...
    'sniff-odor-project_results',...
    'Figuers_data_Marh2023');

[colorSet, SandR_colors] = set_plot_seting(15, 4);
pltDir = fullfile(final_figs_path('uniTn1'),...
    'fig5_E_F');
if ~exist(pltDir , 'dir')
   mkdir(pltDir)
end
fileDir = pltDir;
load(fullfile(dataDir, 'Fig5_cic_and_trans_decoding_of_conc_data',...
    'Fig5_cic_and_trans_decoding_of_conc_data'));
%%
fileID = fopen(fullfile(pltDir,...
        'fig5_E_F_sic_tans_decoding_avrage.txt'),'w');

linseSizePLT = 2.5;
close all
f1 = figure;
hold on
p  = plot([0 6],[100/3,100/3],'linestyle','-.','Color','k');
errorbar(1, nanmean([fullSpaceDecodingAcc_R_o1(:);...
            fullSpaceDecodingAcc_R_o2(:)]),...
            nanstd([fullSpaceDecodingAcc_R_o1(:);...
            fullSpaceDecodingAcc_R_o2(:)]),...
            'Marker', 'o', 'Color', [.5,.5,.5],...
            'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)
errorbar(2, nanmean([fullSpaceDecodingAcc_S_o1(:);...
             fullSpaceDecodingAcc_S_o2(:)]),...
            nanstd([fullSpaceDecodingAcc_S_o1(:);...
            fullSpaceDecodingAcc_S_o2(:)]),...
            'Marker', 'o', 'Color', [.5,.5,.5],...
            'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)
errorbar(3, nanmean([fullSpaceDecodingAcc_Strain_Rtest_o1(:);...
            fullSpaceDecodingAcc_Strain_Rtest_o2(:)]),...
            nanstd([fullSpaceDecodingAcc_Strain_Rtest_o1(:);...
            fullSpaceDecodingAcc_Strain_Rtest_o2(:)]),...
            'Marker', 'o', 'Color', [.5,.5,.5],...
            'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)
errorbar(4, nanmean([fullSpaceDecodingAcc_Rtrain_Stest_o1(:);...
            fullSpaceDecodingAcc_Rtrain_Stest_o2(:)]),...
            nanstd([fullSpaceDecodingAcc_Rtrain_Stest_o1(:);...
            fullSpaceDecodingAcc_Rtrain_Stest_o2(:)]),...
            'Marker', 'o', 'Color', [.5,.5,.5],...
            'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)

yticks([25, 50,75,100])
xticks(1:4);
xticklabels({'R-type', 'S-type',...
    'S-train/R-test', 'R-train/S-test'})
xtickangle(45)
box off
xlim([0,4.5])
ylim([20, 102])

title({'Conc. decoding', 'full space - avrage',''})
ylabel('accuracy(%)')

set(gcf,'Position',[120,120,270,350]);
print_it(pltDir, 'full_space', 'pseudoPopulation')
%%--
fprintf(fileID, 'full space conc decoding in S type acc: %.2f +- %.2f.\n',....
    nanmean([fullSpaceDecodingAcc_S_o1(:);...
            fullSpaceDecodingAcc_S_o2(:)]),...
    nanstd([fullSpaceDecodingAcc_S_o1(:);...
            fullSpaceDecodingAcc_S_o2(:)]));
fprintf(fileID, 'full space conc decoding in R type acc: %.2f +- %.2f.\n',....
    nanmean([fullSpaceDecodingAcc_R_o1(:);...
            fullSpaceDecodingAcc_R_o2(:)]),...
    nanstd([fullSpaceDecodingAcc_R_o1(:);...
            fullSpaceDecodingAcc_R_o2(:)]));
fprintf(fileID, 'full space conc decoding in S train R test  acc: %.2f +- %.2f.\n',....
    nanmean([fullSpaceDecodingAcc_Strain_Rtest_o1(:);...
            fullSpaceDecodingAcc_Strain_Rtest_o2(:)]),...
    nanstd([fullSpaceDecodingAcc_Strain_Rtest_o1(:);...
            fullSpaceDecodingAcc_Strain_Rtest_o2(:)]));
fprintf(fileID, 'full space conc decoding in R train S test  acc: %.2f +- %.2f.\n\n',....
    nanmean([fullSpaceDecodingAcc_Rtrain_Stest_o1(:);...
            fullSpaceDecodingAcc_Rtrain_Stest_o2(:)]),...
    nanstd([fullSpaceDecodingAcc_Rtrain_Stest_o1(:);...
            fullSpaceDecodingAcc_Rtrain_Stest_o2(:)]));
%%--
linseSizePLT = 2.5;
close all
f1 = figure;
hold on
p  = plot([0 30],[100/3,100/3],'linestyle','-.','Color','k');
errorbar(1, nanmean([PC_2_and_3_DecodingAcc_R_o1(:);...
      PC_2_and_3_DecodingAcc_R_o2(:)]),...
      nanstd([PC_2_and_3_DecodingAcc_R_o1(:);...
      PC_2_and_3_DecodingAcc_R_o2(:)]),...
      'Marker', 'o', 'Color', [.5,.5,.5],...
      'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)
errorbar(2, nanmean([PC_2_and_3_DecodingAcc_S_o1(:);...
       PC_2_and_3_DecodingAcc_S_o2(:)]),...
       nanstd([PC_2_and_3_DecodingAcc_S_o1(:);...
       PC_2_and_3_DecodingAcc_S_o2(:)]),...
       'Marker', 'o', 'Color', [.5,.5,.5],...
       'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)
errorbar(3, nanmean([PC_2_and_3_DecodingAcc_Strain_Rtest_o1(:);...
        PC_2_and_3_DecodingAcc_Strain_Rtest_o2(:)]),...
        nanstd([PC_2_and_3_DecodingAcc_Strain_Rtest_o1(:);...
        PC_2_and_3_DecodingAcc_Strain_Rtest_o2(:)]),...
        'Marker', 'o', 'Color', [.5,.5,.5],...
        'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)
errorbar(4, nanmean([PC_2_and_3_DecodingAcc_Rtrain_Stest_o1(:);...
        PC_2_and_3_DecodingAcc_Rtrain_Stest_o2(:)]),...
        nanstd([PC_2_and_3_DecodingAcc_Rtrain_Stest_o1(:);...
        PC_2_and_3_DecodingAcc_Rtrain_Stest_o2(:)]),...
        'Marker', 'o', 'Color', [.5,.5,.5],...
        'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)

yticks([25, 50,75,100])
xticks(1:4);
xticklabels({'R-type', 'S-type',...
    'S-train/R-test', 'R-train/S-test'})
xtickangle(45)
box off
xlim([0,4.5])
ylim([20, 102])

title({'Conc. decoding', 'PC 2 and 3 - avrage',''})
ylabel('accuracy(%)')

set(gcf,'Position',[120,120,270,350]);
print_it(pltDir, 'PC_2_and_3', 'pseudoPopulation')
%%--
fprintf(fileID, 'PC 2 and 3 conc decoding in S type acc: %.2f +- %.2f.\n',....
        nanmean([PC_2_and_3_DecodingAcc_S_o1(:);...
                PC_2_and_3_DecodingAcc_S_o2(:)]),...
        nanstd([PC_2_and_3_DecodingAcc_S_o1(:);...
                PC_2_and_3_DecodingAcc_S_o2(:)]));
fprintf(fileID, 'PC 2 and 3 conc decoding in R type acc: %.2f +- %.2f.\n',....
        nanmean([PC_2_and_3_DecodingAcc_R_o1(:);...
                PC_2_and_3_DecodingAcc_R_o2(:)]),...
        nanstd([PC_2_and_3_DecodingAcc_R_o1(:);...
                PC_2_and_3_DecodingAcc_R_o2(:)]));
fprintf(fileID, 'PC 2 and 3 conc decoding in S train R test  acc: %.2f +- %.2f.\n',....
        nanmean([PC_2_and_3_DecodingAcc_Strain_Rtest_o1(:);...
                PC_2_and_3_DecodingAcc_Strain_Rtest_o2(:)]),...
        nanstd([PC_2_and_3_DecodingAcc_Strain_Rtest_o1(:);...
                PC_2_and_3_DecodingAcc_Strain_Rtest_o2(:)]));
fprintf(fileID, 'PC 2 and 3 conc decoding in R train S test  acc: %.2f +- %.2f.\n\n',....
        nanmean([PC_2_and_3_DecodingAcc_Rtrain_Stest_o1(:);...
                PC_2_and_3_DecodingAcc_Rtrain_Stest_o2(:)]),...
        nanstd([PC_2_and_3_DecodingAcc_Rtrain_Stest_o1(:);...
                PC_2_and_3_DecodingAcc_Rtrain_Stest_o2(:)]));
%
fclose(fileID);
%%
fileID = fopen(fullfile(pltDir,...
        'fig5_E_F_sic_tans_decoding_odor1.txt'),'w');
    
linseSizePLT = 2.5;
close all
f1 = figure;
hold on
p  = plot([0 6],[100/3,100/3],'linestyle','-.','Color','k');
errorbar(1, nanmean(fullSpaceDecodingAcc_R_o1(:)),...
            nanstd(fullSpaceDecodingAcc_R_o1(:)),...
            'Marker', 'o', 'Color', [.5,.5,.5],...
            'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)
errorbar(2, nanmean(fullSpaceDecodingAcc_S_o1(:)),...
            nanstd(fullSpaceDecodingAcc_S_o1(:)),...
            'Marker', 'o', 'Color', [.5,.5,.5],...
            'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)
errorbar(3, nanmean(fullSpaceDecodingAcc_Strain_Rtest_o1(:)),...
            nanstd(fullSpaceDecodingAcc_Strain_Rtest_o1(:)),...
            'Marker', 'o', 'Color', [.5,.5,.5],...
            'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)
errorbar(4, nanmean(fullSpaceDecodingAcc_Rtrain_Stest_o1(:)),...
            nanstd(fullSpaceDecodingAcc_Rtrain_Stest_o1(:)),...
            'Marker', 'o', 'Color', [.5,.5,.5],...
            'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)

yticks([25, 50,75,100])
xticks(1:4);
xticklabels({'R-type', 'S-type',...
    'S-train/R-test', 'R-train/S-test'})
xtickangle(45)
box off
xlim([0,4.5])
ylim([20, 102])

title({'Conc. decoding', 'full space - odor 1',''})
ylabel('accuracy(%)')

set(gcf,'Position',[120,120,270,350]);
print_it(pltDir, 'full_space_odor1', 'pseudoPopulation')
%%--
fprintf(fileID, 'full space conc decoding in S type acc: %.2f +- %.2f.\n',....
    nanmean(fullSpaceDecodingAcc_S_o1(:)),...
    nanstd(fullSpaceDecodingAcc_S_o1(:)));
fprintf(fileID, 'full space conc decoding in R type acc: %.2f +- %.2f.\n',....
    nanmean(fullSpaceDecodingAcc_R_o1(:)),...
    nanstd(fullSpaceDecodingAcc_R_o1(:)));
fprintf(fileID, 'full space conc decoding in S train R test  acc: %.2f +- %.2f.\n',....
    nanmean(fullSpaceDecodingAcc_Strain_Rtest_o1(:)),...
    nanstd(fullSpaceDecodingAcc_Strain_Rtest_o1(:)));
fprintf(fileID, 'full space conc decoding in R train S test  acc: %.2f +- %.2f.\n\n',....
    nanmean(fullSpaceDecodingAcc_Rtrain_Stest_o1(:)),...
    nanstd(fullSpaceDecodingAcc_Rtrain_Stest_o1(:)));
%%--
linseSizePLT = 2.5;
close all
f1 = figure;
hold on
p  = plot([0 30],[100/3,100/3],'linestyle','-.','Color','k');
errorbar(1, nanmean(PC_2_and_3_DecodingAcc_R_o1(:)),...
      nanstd(PC_2_and_3_DecodingAcc_R_o1(:)),...
      'Marker', 'o', 'Color', [.5,.5,.5],...
      'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)
errorbar(2, nanmean(PC_2_and_3_DecodingAcc_S_o1(:)),...
       nanstd(PC_2_and_3_DecodingAcc_S_o1(:)),...
       'Marker', 'o', 'Color', [.5,.5,.5],...
       'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)
errorbar(3, nanmean(PC_2_and_3_DecodingAcc_Strain_Rtest_o1(:)),...
        nanstd(PC_2_and_3_DecodingAcc_Strain_Rtest_o1(:)),...
        'Marker', 'o', 'Color', [.5,.5,.5],...
        'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)
errorbar(4, nanmean(PC_2_and_3_DecodingAcc_Rtrain_Stest_o1(:)),...
        nanstd(PC_2_and_3_DecodingAcc_Rtrain_Stest_o1(:)),...
        'Marker', 'o', 'Color', [.5,.5,.5],...
        'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)

yticks([25, 50,75,100])
xticks(1:4);
xticklabels({'R-type', 'S-type',...
    'S-train/R-test', 'R-train/S-test'})
xtickangle(45)
box off
xlim([0,4.5])
ylim([20, 102])

title({'Conc. decoding', 'PC 2 and 3 - odor 1',''})
ylabel('accuracy(%)')

set(gcf,'Position',[120,120,270,350]);
print_it(pltDir, 'PC_2_and_3_odor1', 'pseudoPopulation')
%%--
fprintf(fileID, 'PC 2 and 3 conc decoding in S type acc: %.2f +- %.2f.\n',....
        nanmean(PC_2_and_3_DecodingAcc_S_o1(:)),...
        nanstd(PC_2_and_3_DecodingAcc_S_o1(:)));
fprintf(fileID, 'PC 2 and 3 conc decoding in R type acc: %.2f +- %.2f.\n',....
        nanmean(PC_2_and_3_DecodingAcc_R_o1(:)),...
        nanstd(PC_2_and_3_DecodingAcc_R_o1(:)));
fprintf(fileID, 'PC 2 and 3 conc decoding in S train R test  acc: %.2f +- %.2f.\n',....
        nanmean(PC_2_and_3_DecodingAcc_Strain_Rtest_o1(:)),...
        nanstd(PC_2_and_3_DecodingAcc_Strain_Rtest_o1(:)));
fprintf(fileID, 'PC 2 and 3 conc decoding in R train S test  acc: %.2f +- %.2f.\n\n',....
        nanmean(PC_2_and_3_DecodingAcc_Rtrain_Stest_o1(:)),...
        nanstd(PC_2_and_3_DecodingAcc_Rtrain_Stest_o1(:)));
%
fclose(fileID);
%%
fileID = fopen(fullfile(pltDir,...
        'fig5_E_F_sic_tans_decoding_odor2.txt'),'w');
    
linseSizePLT = 2.5;
close all
f1 = figure;
hold on
p  = plot([0 6],[100/3,100/3],'linestyle','-.','Color','k');
errorbar(1, nanmean(fullSpaceDecodingAcc_R_o2(:)),...
            nanstd(fullSpaceDecodingAcc_R_o2(:)),...
            'Marker', 'o', 'Color', [.5,.5,.5],...
            'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)
errorbar(2, nanmean(fullSpaceDecodingAcc_S_o2(:)),...
            nanstd(fullSpaceDecodingAcc_S_o2(:)),...
            'Marker', 'o', 'Color', [.5,.5,.5],...
            'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)
errorbar(3, nanmean(fullSpaceDecodingAcc_Strain_Rtest_o2(:)),...
            nanstd(fullSpaceDecodingAcc_Strain_Rtest_o2(:)),...
            'Marker', 'o', 'Color', [.5,.5,.5],...
            'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)
errorbar(4, nanmean(fullSpaceDecodingAcc_Rtrain_Stest_o2(:)),...
            nanstd(fullSpaceDecodingAcc_Rtrain_Stest_o2(:)),...
            'Marker', 'o', 'Color', [.5,.5,.5],...
            'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)

yticks([25, 50,75,100])
xticks(1:4);
xticklabels({'R-type', 'S-type',...
    'S-train/R-test', 'R-train/S-test'})
xtickangle(45)
box off
xlim([0,4.5])
ylim([20, 102])

title({'Conc. decoding', 'full space - odor 2',''})
ylabel('accuracy(%)')

set(gcf,'Position',[120,120,270,350]);
print_it(pltDir, 'full_space_odor2', 'pseudoPopulation')
%%--
fprintf(fileID, 'full space conc decoding in S type acc: %.2f +- %.2f.\n',....
    nanmean(fullSpaceDecodingAcc_S_o2(:)),...
    nanstd(fullSpaceDecodingAcc_S_o2(:)));
fprintf(fileID, 'full space conc decoding in R type acc: %.2f +- %.2f.\n',....
    nanmean(fullSpaceDecodingAcc_R_o2(:)),...
    nanstd(fullSpaceDecodingAcc_R_o2(:)));
fprintf(fileID, 'full space conc decoding in S train R test  acc: %.2f +- %.2f.\n',....
    nanmean(fullSpaceDecodingAcc_Strain_Rtest_o2(:)),...
    nanstd(fullSpaceDecodingAcc_Strain_Rtest_o2(:)));
fprintf(fileID, 'full space conc decoding in R train S test  acc: %.2f +- %.2f.\n\n',....
    nanmean(fullSpaceDecodingAcc_Rtrain_Stest_o2(:)),...
    nanstd(fullSpaceDecodingAcc_Rtrain_Stest_o2(:)));
%%--
linseSizePLT = 2.5;
close all
f1 = figure;
hold on
p  = plot([0 30],[100/3,100/3],'linestyle','-.','Color','k');
errorbar(1, nanmean(PC_2_and_3_DecodingAcc_R_o2(:)),...
      nanstd(PC_2_and_3_DecodingAcc_R_o2(:)),...
      'Marker', 'o', 'Color', [.5,.5,.5],...
      'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)
errorbar(2, nanmean(PC_2_and_3_DecodingAcc_S_o2(:)),...
       nanstd(PC_2_and_3_DecodingAcc_S_o2(:)),...
       'Marker', 'o', 'Color', [.5,.5,.5],...
       'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)
errorbar(3, nanmean(PC_2_and_3_DecodingAcc_Strain_Rtest_o2(:)),...
        nanstd(PC_2_and_3_DecodingAcc_Strain_Rtest_o2(:)),...
        'Marker', 'o', 'Color', [.5,.5,.5],...
        'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)
errorbar(4, nanmean(PC_2_and_3_DecodingAcc_Rtrain_Stest_o2(:)),...
        nanstd(PC_2_and_3_DecodingAcc_Rtrain_Stest_o2(:)),...
        'Marker', 'o', 'Color', [.5,.5,.5],...
        'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 4)

yticks([25, 50,75,100])
xticks(1:4);
xticklabels({'R-type', 'S-type',...
    'S-train/R-test', 'R-train/S-test'})
xtickangle(45)
box off
xlim([0,4.5])
ylim([20, 102])

title({'Conc. decoding', 'PC 2 and 3 - odor 2',''})
ylabel('accuracy(%)')

set(gcf,'Position',[120,120,270,350]);
print_it(pltDir, 'PC_2_and_3_odor2', 'pseudoPopulation')
%%--
fprintf(fileID, 'PC 2 and 3 conc decoding in S type acc: %.2f +- %.2f.\n',....
        nanmean(PC_2_and_3_DecodingAcc_S_o2(:)),...
        nanstd(PC_2_and_3_DecodingAcc_S_o2(:)));
fprintf(fileID, 'PC 2 and 3 conc decoding in R type acc: %.2f +- %.2f.\n',....
        nanmean(PC_2_and_3_DecodingAcc_R_o2(:)),...
        nanstd(PC_2_and_3_DecodingAcc_R_o2(:)));
fprintf(fileID, 'PC 2 and 3 conc decoding in S train R test  acc: %.2f +- %.2f.\n',....
        nanmean(PC_2_and_3_DecodingAcc_Strain_Rtest_o2(:)),...
        nanstd(PC_2_and_3_DecodingAcc_Strain_Rtest_o2(:)));
fprintf(fileID, 'PC 2 and 3 conc decoding in R train S test  acc: %.2f +- %.2f.\n\n',....
        nanmean(PC_2_and_3_DecodingAcc_Rtrain_Stest_o2(:)),...
        nanstd(PC_2_and_3_DecodingAcc_Rtrain_Stest_o2(:)));
%
fclose(fileID);
