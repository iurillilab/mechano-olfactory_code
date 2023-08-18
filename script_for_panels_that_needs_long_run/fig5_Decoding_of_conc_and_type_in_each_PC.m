clear; clc; close all;

dataFilesDirList = {'M1\102721', 'M2\102821','M3\102921',...
    'M5\120221'};
mouseNameList = {'M1_102721', 'M2_102821','M3_102921',...
    'M5_120221'};
dataPathPrefix = 'E:\Ephys\Conc_Series';

%%
[tBeginning, tEnd, trialsToRemove]=...
    get_removal_variables();
add_ndt_paths_and_init_rand_generator
%%
numberOfBalanceSampling = 100;
num_cv_splits = 10;
numOfResample = 100;
the_feature_preprocessors{1} =...
    zscore_normalize_FP;
the_classifier = libsvm_CL;
%%
fullSpaceDecodingAcc_concO1 =...
    nan(numOfResample, numberOfBalanceSampling);
fullSpaceDecodingAcc_concO2 =...
    nan(numOfResample, numberOfBalanceSampling);
fullSpaceDecodingAcc_type =...
    nan(numOfResample, numberOfBalanceSampling);
%%
numberOfSinglePC = 20;
singlePCs_DecodingAcc_concO1 =...
    nan(numOfResample, numberOfSinglePC, numberOfBalanceSampling);
singlePCs_DecodingAcc_concO2 =...
    nan( numOfResample, numberOfSinglePC, numberOfBalanceSampling);
singlePCs_DecodingAcc_typeO1 =...
    nan( numOfResample, numberOfSinglePC, numberOfBalanceSampling);
singlePCs_DecodingAcc_typeO2 =...
    nan( numOfResample, numberOfSinglePC, numberOfBalanceSampling);
conc_set = [1,2,4];
%%
parfor perm = 1 : numberOfBalanceSampling
    [fSpaceCell, fSpaceLbelsCell] =...
        get_aSample(dataPathPrefix, mouseNameList,....
                    conc_set, 1:2, tBeginning, tEnd,trialsToRemove,...
                    1, 0, 6, .5);
    %%
    pseudoPopulation =...
        [fSpaceCell{1}, fSpaceCell{2},fSpaceCell{3}, fSpaceCell{4}];
    [~,score,latent] = pca(zscore(pseudoPopulation));
    fSpaceLbels = fSpaceLbelsCell{1};
    %%
    onePermSinglePCs_DecodingAcc_concO1 =...
        nan(numOfResample, numberOfSinglePC);
    onePermSinglePCs_DecodingAcc_concO2 =...
        nan(numOfResample, numberOfSinglePC);
    onePermSinglePCs_DecodingAcc_typeO1 =...
        nan(numOfResample, numberOfSinglePC);
    onePermSinglePCs_DecodingAcc_typeO2 =...
        nan(numOfResample, numberOfSinglePC);
    %%
    for pcNum = 1 : numberOfSinglePC
        %
        totalSpace = score(:, pcNum);
        %%
        binedDataConcO1 = cell(1,size(totalSpace,2));
        binedLabeConcO1 = cell(1,size(totalSpace,2));
        for ni = 1 : size(totalSpace,2)
            binedDataConcO1{ni} =...
                totalSpace(fSpaceLbels.OdorId == 1, ni);
            binedLabeConcO1{ni} =...
                fSpaceLbels.ConcId(fSpaceLbels.OdorId == 1);
        end
        DS_concO1 = basic_DS(binedDataConcO1,...
            binedLabeConcO1, num_cv_splits);
        %%
        binedDataConcO2 = cell(1,size(totalSpace,2));
        binedLabeConcO2 = cell(1,size(totalSpace,2));
        for ni = 1 : size(totalSpace,2)
            binedDataConcO2{ni} =...
                totalSpace(fSpaceLbels.OdorId == 2, ni);
            binedLabeConcO2{ni} =...
                fSpaceLbels.ConcId(fSpaceLbels.OdorId == 2);
        end
        DS_concO2 = basic_DS(binedDataConcO2,...
            binedLabeConcO2, num_cv_splits);
        %%
        the_cross_validator_ConcO1 =...
            standard_resample_CV(DS_concO1, the_classifier, the_feature_preprocessors);
        the_cross_validator_ConcO1.num_resample_runs = numOfResample;
        the_cross_validator_ConcO1.display_progress.resample_run_time = 0;
        the_cross_validator_ConcO1.display_progress.zero_one_loss = 0;
        the_cross_validator_ConcO1.test_only_at_training_times = 1;
        DECODING_RESULTS_ConcO1 = the_cross_validator_ConcO1.run_cv_decoding;
        onePermSinglePCs_DecodingAcc_concO1(:,pcNum) =...
            mean(DECODING_RESULTS_ConcO1.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
        %%
        the_cross_validator_ConcO2 =...
            standard_resample_CV(DS_concO2, the_classifier, the_feature_preprocessors);
        the_cross_validator_ConcO2.num_resample_runs = numOfResample;
        the_cross_validator_ConcO2.display_progress.resample_run_time = 0;
        the_cross_validator_ConcO2.display_progress.zero_one_loss = 0;
        the_cross_validator_ConcO2.test_only_at_training_times = 1;
        DECODING_RESULTS_ConcO2 = the_cross_validator_ConcO2.run_cv_decoding;
        onePermSinglePCs_DecodingAcc_concO2(:,pcNum) =...
            mean(DECODING_RESULTS_ConcO2.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
        %%
        binedDataTypeO1 = cell(1,size(totalSpace,2));
        binedLabeTypeO1 = cell(1,size(totalSpace,2));
        for ni = 1 : size(totalSpace,2)
            binedDataTypeO1{ni} =...
                totalSpace(fSpaceLbels.OdorId == 1,ni);
            binedLabeTypeO1{ni} =...
                fSpaceLbels.RespLabel(fSpaceLbels.OdorId == 1);
        end
        DS_Type_O1 = basic_DS(binedDataTypeO1,...
            binedLabeTypeO1, num_cv_splits);
        %%
        binedDataTypeO2 = cell(1,size(totalSpace,2));
        binedLabeTypeO2 = cell(1,size(totalSpace,2));
        for ni = 1 : size(totalSpace,2)
            binedDataTypeO2{ni} =...
                totalSpace(fSpaceLbels.OdorId == 2,ni);
            binedLabeTypeO2{ni} =...
                fSpaceLbels.RespLabel(fSpaceLbels.OdorId == 2);
        end
        DS_Type_O2 = basic_DS(binedDataTypeO2,...
            binedLabeTypeO2, num_cv_splits);
        %%
        the_cross_validator_typeO1 =...
            standard_resample_CV(DS_Type_O1, the_classifier, the_feature_preprocessors);
        the_cross_validator_typeO1.num_resample_runs = numOfResample;
        the_cross_validator_typeO1.display_progress.resample_run_time = 0;
        the_cross_validator_typeO1.display_progress.zero_one_loss = 0;
        the_cross_validator_typeO1.test_only_at_training_times = 1;
        DECODING_RESULTS_typeO1 = the_cross_validator_typeO1.run_cv_decoding;
        onePermSinglePCs_DecodingAcc_typeO1(:, pcNum) =...
            mean(DECODING_RESULTS_typeO1.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
        %%
        the_cross_validator_typeO2 =...
            standard_resample_CV(DS_Type_O2, the_classifier, the_feature_preprocessors);
        the_cross_validator_typeO2.num_resample_runs = numOfResample;
        the_cross_validator_typeO2.display_progress.resample_run_time = 0;
        the_cross_validator_typeO2.display_progress.zero_one_loss = 0;
        the_cross_validator_typeO2.test_only_at_training_times = 1;
        DECODING_RESULTS_typeO2 = the_cross_validator_typeO2.run_cv_decoding;
        onePermSinglePCs_DecodingAcc_typeO2(:, pcNum) =...
            mean(DECODING_RESULTS_typeO2.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
    end
    %%
    singlePCs_DecodingAcc_concO1(:,:,perm) =...
        onePermSinglePCs_DecodingAcc_concO1;
    singlePCs_DecodingAcc_concO2(:,:,perm) =...
        onePermSinglePCs_DecodingAcc_concO2;
    singlePCs_DecodingAcc_typeO1(:,:,perm) =...
        onePermSinglePCs_DecodingAcc_typeO1;
    singlePCs_DecodingAcc_typeO2(:,:,perm) =...
        onePermSinglePCs_DecodingAcc_typeO2;
end
%%
pltDir = fullfile(final_figs_path('beast1'),...
    'Fig5_Decoding_of_conc_and_type_in_each_PC_data');
if ~exist(pltDir , 'dir')
   mkdir(pltDir)
end
%%
save(fullfile(pltDir,'Fig5_Decoding_of_conc_and_type_in_each_PC_data.mat'),...
    'singlePCs_DecodingAcc_concO1',...
    'singlePCs_DecodingAcc_concO2',...
    'singlePCs_DecodingAcc_typeO1',...
    "singlePCs_DecodingAcc_typeO2",...
    '-v7.3');
