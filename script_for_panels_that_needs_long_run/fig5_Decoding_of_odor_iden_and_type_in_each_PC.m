clear; clc; close all;

mouseNameList = {'M1_050622', 'M2_051022'};
% dataPathPrefix = 'E:\Ephys\sniffOdorProject\5Od_2Conc_PCx';
dataPathPrefix = 'D:\sniffOdorProject\5Od_2Conc_PCx';
%%
[tBeginning, tEnd, trialsToRemove]=...
    get_removal_variables();
add_ndt_paths_and_init_rand_generator
%%
numberOfBalanceSampling = 50
num_cv_splits = 10;
numOfResample = 50;
the_feature_preprocessors{1} =...
    zscore_normalize_FP;
the_classifier = libsvm_CL;
%%
fullSpaceDecodingAcc_idenC1 =...
    nan(numOfResample, numberOfBalanceSampling);
fullSpaceDecodingAcc_idenC2 =...
    nan(numOfResample, numberOfBalanceSampling);
fullSpaceDecodingAcc_type =...
    nan(numOfResample, numberOfBalanceSampling);
%%
numberOfSinglePC = 20;
singlePCs_DecodingAcc_idenC1 =...
    nan(numOfResample, numberOfSinglePC, numberOfBalanceSampling);
singlePCs_DecodingAcc_idenC2 =...
    nan( numOfResample, numberOfSinglePC, numberOfBalanceSampling);
singlePCs_DecodingAcc_typeC1 =...
    nan( numOfResample, numberOfSinglePC, numberOfBalanceSampling);
singlePCs_DecodingAcc_typeC2 =...
    nan( numOfResample, numberOfSinglePC, numberOfBalanceSampling);
conc_set = [1, 2];
%%
for perm = 1 : numberOfBalanceSampling
    [fSpaceCell, fSpaceLbelsCell] =...
        get_aSample(dataPathPrefix, mouseNameList,....
                    conc_set, 1:5, tBeginning, tEnd,trialsToRemove,...
                    1, 0, 6, .5);
    %%
    pseudoPopulation =...
        [fSpaceCell{1}, fSpaceCell{2}];
    [~,score,latent] = pca(zscore(pseudoPopulation));
    fSpaceLbels = fSpaceLbelsCell{1};
    %%
    onePermSinglePCs_DecodingAcc_idenC1 =...
        nan(numOfResample, numberOfSinglePC);
    onePermSinglePCs_DecodingAcc_idenC2 =...
        nan(numOfResample, numberOfSinglePC);
    onePermSinglePCs_DecodingAcc_typeC1 =...
        nan(numOfResample, numberOfSinglePC);
    onePermSinglePCs_DecodingAcc_typeC2 =...
        nan(numOfResample, numberOfSinglePC);
    %%
    for pcNum = 1 : numberOfSinglePC
        %
        totalSpace = score(:, pcNum);
        %%
        binedDataIdenC1 = cell(1,size(totalSpace,2));
        binedLabeIdenC1 = cell(1,size(totalSpace,2));
        for ni = 1 : size(totalSpace,2)
            binedDataIdenC1{ni} =...
                totalSpace(fSpaceLbels.ConcId == 1, ni);
            binedLabeIdenC1{ni} =...
                fSpaceLbels.OdorId(fSpaceLbels.ConcId == 1);
        end
        DS_idenC1 = basic_DS(binedDataIdenC1,...
            binedLabeIdenC1, num_cv_splits);
        %%
        binedDataIdenC2 = cell(1,size(totalSpace,2));
        binedLabeIdenC2 = cell(1,size(totalSpace,2));
        for ni = 1 : size(totalSpace,2)
            binedDataIdenC2{ni} =...
                totalSpace(fSpaceLbels.ConcId == 2, ni);
            binedLabeIdenC2{ni} =...
                fSpaceLbels.OdorId(fSpaceLbels.ConcId == 2);
        end
        DS_idenC2 = basic_DS(binedDataIdenC2,...
            binedLabeIdenC2, num_cv_splits);
        %%
        the_cross_validator_IdenC1 =...
            standard_resample_CV(DS_idenC1, the_classifier, the_feature_preprocessors);
        the_cross_validator_IdenC1.num_resample_runs = numOfResample;
        the_cross_validator_IdenC1.display_progress.resample_run_time = 0;
        the_cross_validator_IdenC1.display_progress.zero_one_loss = 0;
        the_cross_validator_IdenC1.test_only_at_training_times = 1;
        DECODING_RESULTS_IdenC1 = the_cross_validator_IdenC1.run_cv_decoding;
        onePermSinglePCs_DecodingAcc_idenC1(:,pcNum) =...
            mean(DECODING_RESULTS_IdenC1.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
        %%
        the_cross_validator_IdenC2 =...
            standard_resample_CV(DS_idenC2, the_classifier, the_feature_preprocessors);
        the_cross_validator_IdenC2.num_resample_runs = numOfResample;
        the_cross_validator_IdenC2.display_progress.resample_run_time = 0;
        the_cross_validator_IdenC2.display_progress.zero_one_loss = 0;
        the_cross_validator_IdenC2.test_only_at_training_times = 1;
        DECODING_RESULTS_IdenC2 = the_cross_validator_IdenC2.run_cv_decoding;
        onePermSinglePCs_DecodingAcc_idenC2(:,pcNum) =...
            mean(DECODING_RESULTS_IdenC2.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
        %%
        binedDataTypeO1 = cell(1,size(totalSpace,2));
        binedLabeTypeO1 = cell(1,size(totalSpace,2));
        for ni = 1 : size(totalSpace,2)
            binedDataTypeO1{ni} =...
                totalSpace(fSpaceLbels.ConcId == 1,ni);
            binedLabeTypeO1{ni} =...
                fSpaceLbels.RespLabel(fSpaceLbels.ConcId == 1);
        end
        DS_Type_O1 = basic_DS(binedDataTypeO1,...
            binedLabeTypeO1, num_cv_splits);
        %%
        binedDataTypeO2 = cell(1,size(totalSpace,2));
        binedLabeTypeO2 = cell(1,size(totalSpace,2));
        for ni = 1 : size(totalSpace,2)
            binedDataTypeO2{ni} =...
                totalSpace(fSpaceLbels.ConcId == 2,ni);
            binedLabeTypeO2{ni} =...
                fSpaceLbels.RespLabel(fSpaceLbels.ConcId == 2);
        end
        DS_Type_O2 = basic_DS(binedDataTypeO2,...
            binedLabeTypeO2, num_cv_splits);
        %%
        the_cross_validator_typeC1 =...
            standard_resample_CV(DS_Type_O1, the_classifier, the_feature_preprocessors);
        the_cross_validator_typeC1.num_resample_runs = numOfResample;
        the_cross_validator_typeC1.display_progress.resample_run_time = 0;
        the_cross_validator_typeC1.display_progress.zero_one_loss = 0;
        the_cross_validator_typeC1.test_only_at_training_times = 1;
        DECODING_RESULTS_typeC1 = the_cross_validator_typeC1.run_cv_decoding;
        onePermSinglePCs_DecodingAcc_typeC1(:, pcNum) =...
            mean(DECODING_RESULTS_typeC1.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
        %%
        the_cross_validator_typeC2 =...
            standard_resample_CV(DS_Type_O2, the_classifier, the_feature_preprocessors);
        the_cross_validator_typeC2.num_resample_runs = numOfResample;
        the_cross_validator_typeC2.display_progress.resample_run_time = 0;
        the_cross_validator_typeC2.display_progress.zero_one_loss = 0;
        the_cross_validator_typeC2.test_only_at_training_times = 1;
        DECODING_RESULTS_typeC2 = the_cross_validator_typeC2.run_cv_decoding;
        onePermSinglePCs_DecodingAcc_typeC2(:, pcNum) =...
            mean(DECODING_RESULTS_typeC2.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
    end
    %%
    singlePCs_DecodingAcc_idenC1(:,:,perm) =...
        onePermSinglePCs_DecodingAcc_idenC1;
    singlePCs_DecodingAcc_idenC2(:,:,perm) =...
        onePermSinglePCs_DecodingAcc_idenC2;
    singlePCs_DecodingAcc_typeC1(:,:,perm) =...
        onePermSinglePCs_DecodingAcc_typeC1;
    singlePCs_DecodingAcc_typeC2(:,:,perm) =...
        onePermSinglePCs_DecodingAcc_typeC2;
end
%%
pltDir = fullfile(final_figs_path('uniTn1'),...
    'Fig5_Decoding_of_odor_iden_and_type_in_each_PC');
if ~exist(pltDir , 'dir')
   mkdir(pltDir)
end
%%
save(fullfile(pltDir,'Fig5_Decoding_of_odor_iden_and_type_in_each_PC.mat'),...
    'singlePCs_DecodingAcc_idenC1',...
    'singlePCs_DecodingAcc_idenC2',...
    'singlePCs_DecodingAcc_typeC1',...
    "singlePCs_DecodingAcc_typeC2",...
    '-v7.3');
