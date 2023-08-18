clear; clc; close all;
[dataFilesDirList, mouseNameList, dataPathPrefixList] =....
    get_all_healthy_mice('beast1');
%% dfine countainer
add_ndt_paths_and_init_rand_generator
%%
numberOfBlanceSampling = 100;
num_cv_splits = 10;
numOfResample = 100;
the_feature_preprocessors{1} =...
    zscore_normalize_FP;
the_classifier = libsvm_CL;
blanceSmapleSize = 400;
%%
each_mice_population_inh_decodind_cell = cell(size(mouseNameList));
%%
for mn = 1 : length(mouseNameList)    
    %%
    decodingResultsMousePapulation =... 
        nan(numberOfBlanceSampling, numOfResample); 
    %%
    parfor perm = 1 : numberOfBlanceSampling
        [pseudoPopulation, labelsSet] =...
            get_baiselineSample(mouseNameList(mn),...
            dataPathPrefixList(mn),...
            6, blanceSmapleSize);
        %%
        totalSpace = pseudoPopulation;
        binedData = cell(1,size(totalSpace,2));
        binedLabe = cell(1,size(totalSpace,2));
        for ni = 1 : size(totalSpace,2)
            binedData{ni} =...
                totalSpace(:,ni);
            binedLabe{ni} =...
                labelsSet;
        end
        regionsBasicDs =...
            basic_DS(binedData, binedLabe, num_cv_splits);
        %%
        inhTypeDs = regionsBasicDs;
        the_cross_validator_inhType =...
            standard_resample_CV(inhTypeDs, the_classifier, the_feature_preprocessors);
        the_cross_validator_inhType.num_resample_runs = numOfResample;
        the_cross_validator_inhType.display_progress.resample_run_time = 0;
        the_cross_validator_inhType.display_progress.zero_one_loss = 0;
        the_cross_validator_inhType.test_only_at_training_times = 1;
        DECODING_RESULTS_inhTypeDs_pooled =...
            the_cross_validator_inhType.run_cv_decoding;
        decodingResultsMousePapulation(perm, :) =...
            mean(DECODING_RESULTS_inhTypeDs_pooled.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
    end
    each_mice_population_inh_decodind_cell{mn} = ...
        decodingResultsMousePapulation;
end
%%
pltDir = fullfile(final_figs_path('beast1'),...
    'fig3_decoding_of_inh_type_data_single_mice');
if ~exist(pltDir , 'dir')
   mkdir(pltDir)
end

fileDir = pltDir;
%%
save(fullfile(pltDir,'fig3_decoding_of_inh_type_data_single_mice.mat'),...
    'each_mice_population_inh_decodind_cell', '-v7.3');
