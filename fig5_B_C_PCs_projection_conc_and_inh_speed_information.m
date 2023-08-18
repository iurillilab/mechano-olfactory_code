clear; clc; close all;

dataDir = fullfile('C:\Users\Alireza',...
    'OneDrive - Fondazione Istituto Italiano Tecnologia',...
    'sniff-odor-project_results',...
    'Figuers_data_Marh2023');

[colorSet, SandR_colors] = set_plot_seting(15, 4);
pltDir = fullfile(final_figs_path('uniTn1'),...
    'fig5_B_C');
if ~exist(pltDir , 'dir')
   mkdir(pltDir)
end
fileDir = pltDir;
load(fullfile(dataDir, 'Fig5_Decoding_of_conc_and_type_in_each_PC_data',...
    'Fig5_Decoding_of_conc_and_type_in_each_PC_data'));
%%
conc_decoding_in_each_PC =...
   cat(1, singlePCs_DecodingAcc_concO1,...
   singlePCs_DecodingAcc_concO2);
singlePCs_DecodingAcc_type = ...
    cat(1, singlePCs_DecodingAcc_typeO1,...
   singlePCs_DecodingAcc_typeO2);

mean_conc_decoding_acc_in_each_PC = ...
    mean(conc_decoding_in_each_PC, [1,3]);
std_conc_decoding_acc_in_each_PC = ...
    std(conc_decoding_in_each_PC, [],[1,3]);
%%
mean_inh_type_decoding_acc_in_each_PC = ...
    mean(singlePCs_DecodingAcc_type, [1,3]);
std_inh_type_decoding_acc_in_each_PC = ...
    std(singlePCs_DecodingAcc_type, [],[1,3]);
%%
linseSizePLT = 2.5;
close all
hold on
num_of_pc_to_plot = 15;
for pcId = 1: num_of_pc_to_plot
    errorbar(pcId,mean_conc_decoding_acc_in_each_PC(pcId),...
        std_conc_decoding_acc_in_each_PC(pcId),...
        'Marker', 'o', 'Color', [.5,.5,.5],...
        'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 5);
end
plot([0,num_of_pc_to_plot],[100,100]/3,'linestyle','--','Color','k',...
    'linewidth',2);
ylim([20, 102])
yticks([25, 50,75,100])
ylabel('accuracy (%)')
xlabel('#PC')
title('conc. decoding')
box off
set(gcf,'Position',[200,200,525,320]);
print_it(pltDir, 'conc_decoding_single_PC_projection',...
    'pseudoPopulation');
%%
linseSizePLT = 2.5;
close all
hold on
num_of_pc_to_plot = 15;
for pcId = 1: num_of_pc_to_plot
    errorbar(pcId,mean_inh_type_decoding_acc_in_each_PC(pcId),...
        std_inh_type_decoding_acc_in_each_PC(pcId),...
        'Marker', 'o', 'Color', [.5,.5,.5],...
        'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 5);
end
plot([0,num_of_pc_to_plot],[100,100]/2,'linestyle','--','Color','k',...
    'linewidth',2);
ylim([20, 102])
yticks([25, 50,75,100])
ylabel('accuracy (%)')
xlabel('#PC')
title('inh. speed decoding')
box off
set(gcf,'Position',[200,200,525,320]);
print_it(pltDir, 'inh_speed_decoding_single_PC_projection',...
    'pseudoPopulation');
%%
file_id = fopen(fullfile(pltDir, 'decoding_in_each_PCx.txt'),'w');
for PCi = 1 : 15
    fprintf(file_id, 'PC %d: conc. %.2f +- %.2f %%; inh. %.2f +%.2f %%. \n',...
        PCi, mean_conc_decoding_acc_in_each_PC(PCi),...
        std_conc_decoding_acc_in_each_PC(PCi),...
        mean_inh_type_decoding_acc_in_each_PC(PCi),...
        std_inh_type_decoding_acc_in_each_PC(PCi));
end
fclose(file_id);
%%
for oId = 1 : 2
    switch oId
        case 1
            mean_conc_decoding_acc_in_each_PC = ...
                mean(singlePCs_DecodingAcc_concO1, [1,3]);
            std_conc_decoding_acc_in_each_PC = ...
                std(singlePCs_DecodingAcc_concO1, [],[1,3]);

            mean_inh_type_decoding_acc_in_each_PC = ...
                mean(singlePCs_DecodingAcc_typeO1  , [1,3]);
            std_inh_type_decoding_acc_in_each_PC = ...
                std(singlePCs_DecodingAcc_typeO1, [],[1,3]);
        case 2
            mean_conc_decoding_acc_in_each_PC = ...
                mean(singlePCs_DecodingAcc_concO2, [1,3]);
            std_conc_decoding_acc_in_each_PC = ...
                std(singlePCs_DecodingAcc_concO2, [],[1,3]);

            mean_inh_type_decoding_acc_in_each_PC = ...
                mean(singlePCs_DecodingAcc_typeO2  , [1,3]);
            std_inh_type_decoding_acc_in_each_PC = ...
                std(singlePCs_DecodingAcc_typeO2, [],[1,3]);
    end
            
    close all
    hold on
    num_of_pc_to_plot = 15;
    for pcId = 1: num_of_pc_to_plot
           errorbar(pcId,mean_conc_decoding_acc_in_each_PC(pcId),...
            std_conc_decoding_acc_in_each_PC(pcId),...
            'Marker', 'o', 'Color', [.5,.5,.5],...
            'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 5);
    end
    plot([0,num_of_pc_to_plot],[100,100]/3,'linestyle','--','Color','k',...
        'linewidth',2);
    ylim([20, 102])
    yticks([25, 50,75,100])
    ylabel('accuracy (%)')
    xlabel('#PC')
    title(sprintf('conc. decoding - odor %d', oId))
    box off
    set(gcf,'Position',[200,200,525,320]);
    print_it(pltDir, sprintf('conc_decoding_single_PC_projection_odor%d', oId),...
        'pseudoPopulation');
    %%
    linseSizePLT = 2.5;
    close all
    hold on
    num_of_pc_to_plot = 15;
    for pcId = 1: num_of_pc_to_plot
        errorbar(pcId,mean_inh_type_decoding_acc_in_each_PC(pcId),...
            std_inh_type_decoding_acc_in_each_PC(pcId),...
            'Marker', 'o', 'Color', [.5,.5,.5],...
            'CapSize', 3, 'LineWidth', linseSizePLT, 'MarkerSize', 5);
    end
    plot([0,num_of_pc_to_plot],[100,100]/2,'linestyle','--','Color','k',...
        'linewidth',2);
    ylim([20, 102])
    yticks([25, 50,75,100])
    ylabel('accuracy (%)')
    xlabel('#PC')
    title(sprintf('inh. speed decoding - odor %d', oId))
    box off
    set(gcf,'Position',[200,200,525,320]);
    print_it(pltDir, sprintf('inh_speed_decoding_single_PC_projection_odor%d', oId),...
        'pseudoPopulation');
    %%
    file_id = fopen(fullfile(pltDir,...
        sprintf('decoding_in_each_PCx_%d.txt',oId)),'w');
    for PCi = 1 : 15
        fprintf(file_id, 'PC %d: conc. %.2f +- %.2f %%; inh. %.2f +%.2f %%. \n',...
            PCi, mean_conc_decoding_acc_in_each_PC(PCi),...
            std_conc_decoding_acc_in_each_PC(PCi),...
            mean_inh_type_decoding_acc_in_each_PC(PCi),...
            std_inh_type_decoding_acc_in_each_PC(PCi));
    end
    fclose(file_id);

end