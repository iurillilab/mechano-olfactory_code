clear; clc; close all;
mouseNameList = {'M1_122121', 'M4_011522', 'M5_021022', 'M6_021122'};
dataPathPrefixList = {'E:\Ephys\Bulbectomy', 'E:\Ephys\Bulbectomy',...
    'E:\Ephys\Bulbectomy', 'E:\Ephys\Bulbectomy'};
%% dfine countainer
add_ndt_paths_and_init_rand_generator
%%
regionsAcroNames = {'Octx', 'MCtx'};
regionsNames = {'priform','motor'};
regionsMappedId = [6,1];
%%
numberOfBlanceSampling = 100;
num_cv_splits = 10;
numOfResample = 100;
the_feature_preprocessors{1} =...
    zscore_normalize_FP;
the_classifier = libsvm_CL;
blanceSmapleSize = 400;
%%
maxNumberOfNeurons = 175;
num_of_neurons_set = 5:10:maxNumberOfNeurons;
%%
for ri = 1 : length(regionsMappedId)    
    %%
    rgName = regionsNames{ri};
    regionDecodingResults =... 
            nan(numberOfBlanceSampling, numOfResample, maxNumberOfNeurons);
    regionDecodingResultsMaxUnits =... 
        nan(numberOfBlanceSampling, numOfResample); 
    %%
    parfor perm = 1 : numberOfBlanceSampling
        [pseudoPopulation, labelsSet] =...
            get_baiselineSample(mouseNameList, dataPathPrefixList,...
            regionsMappedId(ri), blanceSmapleSize);
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
        regionDecodingResultsOnePerm =... 
            nan(numOfResample, maxNumberOfNeurons);
        for nn = num_of_neurons_set
            if size(pseudoPopulation,2) > nn
                inhTypeDs = regionsBasicDs;
                inhTypeDs.sample_sites_with_replacement = 1;
                inhTypeDs.num_resample_sites = nn;
                the_cross_validator_inhType =...
                    standard_resample_CV(inhTypeDs, the_classifier, the_feature_preprocessors);
                the_cross_validator_inhType.num_resample_runs = numOfResample;
                the_cross_validator_inhType.display_progress.resample_run_time = 0;
                the_cross_validator_inhType.display_progress.zero_one_loss = 0;
                the_cross_validator_inhType.test_only_at_training_times = 1;
                DECODING_RESULTS_inhTypeDs_pooled = the_cross_validator_inhType.run_cv_decoding;
                regionDecodingResultsOnePerm(:,nn) =...
                    mean(DECODING_RESULTS_inhTypeDs_pooled.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
            end
        end
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
        regionDecodingResultsMaxUnits(perm, :) =...
            mean(DECODING_RESULTS_inhTypeDs_pooled.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
        %%
        regionDecodingResults(perm, :, :) =...
            regionDecodingResultsOnePerm;
    end
    decodingResults.(rgName) = regionDecodingResults;
    decodingResultsMaxUnits.(rgName) = regionDecodingResultsMaxUnits;
end
%%
[colorSet, SandR_colors] = set_plot_seting(20, 8);
pltDir = fullfile(final_figs_path('beast1'),...
    'fig3_decoding_of_inh_type_data','bulbectomized');
if ~exist(pltDir , 'dir')
   mkdir(pltDir)
end
%%
save(fullfile(pltDir,...
    'fig3_decoding_of_inh_type_data_bulbectomized.mat'),...
    'decodingResults', 'decodingResultsMaxUnits',...
    'num_of_neurons_set','-v7.3');