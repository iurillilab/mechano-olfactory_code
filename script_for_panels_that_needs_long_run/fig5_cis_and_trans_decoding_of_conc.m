clear; clc; close all;

dataFilesDirList = {'M1\102721', 'M2\102821','M3\102921',...
    'M5\120221'};
mouseNameList = {'M1_102721', 'M2_102821','M3_102921',...
    'M5_120221'};
 dataPathPrefix = 'E:\Ephys\Conc_Series';
% dataPathPrefix = 'D:\Reza\sniffOdorProject\Conc_Series';

%%
[tBeginning, tEnd, trialsToRemove]=...
    get_removal_variables();
add_ndt_paths_and_init_rand_generator
conc_set = [1,2,4];

%%
numberOfBalanceSampling = 100;
num_cv_splits = 10;
numOfResample = 100;
the_feature_preprocessors{1} =...
    zscore_normalize_FP;
the_classifier = libsvm_CL;
%%
fullSpaceDecodingAcc_R_o1 =...
    nan(numOfResample, numberOfBalanceSampling);
fullSpaceDecodingAcc_S_o1 =...
    nan(numOfResample, numberOfBalanceSampling);
fullSpaceDecodingAcc_Rtrain_Stest_o1 =...
    nan(numOfResample, numberOfBalanceSampling);
fullSpaceDecodingAcc_Strain_Rtest_o1 =...
    nan(numOfResample, numberOfBalanceSampling);

fullSpaceDecodingAcc_R_o2 =...
    nan(numOfResample, numberOfBalanceSampling);
fullSpaceDecodingAcc_S_o2 =...
    nan(numOfResample, numberOfBalanceSampling);
fullSpaceDecodingAcc_Rtrain_Stest_o2 =...
    nan(numOfResample, numberOfBalanceSampling);
fullSpaceDecodingAcc_Strain_Rtest_o2 =...
    nan(numOfResample, numberOfBalanceSampling);
%%
PC_2_and_3_DecodingAcc_R_o1 =...
    nan(numOfResample, numberOfBalanceSampling);
PC_2_and_3_DecodingAcc_S_o1 =...
    nan(numOfResample, numberOfBalanceSampling);
PC_2_and_3_DecodingAcc_Rtrain_Stest_o1 =...
    nan(numOfResample, numberOfBalanceSampling);
PC_2_and_3_DecodingAcc_Strain_Rtest_o1 =...
    nan(numOfResample, numberOfBalanceSampling);

PC_2_and_3_DecodingAcc_R_o2 =...
    nan(numOfResample, numberOfBalanceSampling);
PC_2_and_3_DecodingAcc_S_o2 =...
    nan(numOfResample, numberOfBalanceSampling);
PC_2_and_3_DecodingAcc_Rtrain_Stest_o2 =...
    nan(numOfResample, numberOfBalanceSampling);
PC_2_and_3_DecodingAcc_Strain_Rtest_o2 =...
    nan(numOfResample, numberOfBalanceSampling);
%%
for perm = 1 : numberOfBalanceSampling
    %%
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
    totalSpace = pseudoPopulation;
    %%
    binedDataStype_o1 = cell(1,size(totalSpace,2));
    binedLabeStype_o1 = cell(1,size(totalSpace,2));
    for ni = 1 : size(totalSpace,2)
        binedDataStype_o1{ni} =...
            totalSpace(fSpaceLbels.OdorId == 1 &...
                       fSpaceLbels.RespLabel == 1, ni);
        binedLabeStype_o1{ni} =...
            fSpaceLbels.ConcId(fSpaceLbels.OdorId == 1 &...
                               fSpaceLbels.RespLabel == 1);
    end
    DS_S_O1 = basic_DS(binedDataStype_o1,...
        binedLabeStype_o1, num_cv_splits);
    %%
    binedDataStype_o2 = cell(1,size(totalSpace,2));
    binedLabeStype_o2 = cell(1,size(totalSpace,2));
    for ni = 1 : size(totalSpace,2)
        binedDataStype_o2{ni} =...
            totalSpace(fSpaceLbels.OdorId == 2 &...
                       fSpaceLbels.RespLabel == 1, ni);
        binedLabeStype_o2{ni} =...
            fSpaceLbels.ConcId(fSpaceLbels.OdorId == 2 &...
                               fSpaceLbels.RespLabel == 1);
    end
    DS_S_O2 = basic_DS(binedDataStype_o2,...
        binedLabeStype_o2, num_cv_splits);
    %%
    binedDataRtype_o1 = cell(1,size(totalSpace,2));
    binedLabeRtype_o1 = cell(1,size(totalSpace,2));
    for ni = 1 : size(totalSpace,2)
        binedDataRtype_o1{ni} =...
            totalSpace(fSpaceLbels.OdorId == 1 &...
                       fSpaceLbels.RespLabel == 2, ni);
        binedLabeRtype_o1{ni} =...
            fSpaceLbels.ConcId(fSpaceLbels.OdorId == 1 &...
                               fSpaceLbels.RespLabel == 2);
    end
    DS_R_O1 = basic_DS(binedDataRtype_o1,...
        binedLabeRtype_o1, num_cv_splits);
    %%
    binedDataRtype_o2 = cell(1,size(totalSpace,2));
    binedLabeRtype_o2 = cell(1,size(totalSpace,2));
    for ni = 1 : size(totalSpace,2)
        binedDataRtype_o2{ni} =...
            totalSpace(fSpaceLbels.OdorId == 2 &...
                       fSpaceLbels.RespLabel == 2, ni);
        binedLabeRtype_o2{ni} =...
            fSpaceLbels.ConcId(fSpaceLbels.OdorId == 2 &...
                               fSpaceLbels.RespLabel == 2);
    end
    DS_R_O2 = basic_DS(binedDataRtype_o2,...
        binedLabeRtype_o2, num_cv_splits);
    %%
    the_cross_validator_S_o1 =...
        standard_resample_CV(DS_S_O1, the_classifier, the_feature_preprocessors);
    the_cross_validator_S_o1.num_resample_runs = numOfResample;
    the_cross_validator_S_o1.display_progress.resample_run_time = 0;
    the_cross_validator_S_o1.display_progress.zero_one_loss = 0;
    the_cross_validator_S_o1.test_only_at_training_times = 1;
    DECODING_RESULTS_s_o1 = the_cross_validator_S_o1.run_cv_decoding;
    fullSpaceDecodingAcc_S_o1(:,perm) =...
        mean(DECODING_RESULTS_s_o1.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
    %%
    the_cross_validator_S_o2 =...
        standard_resample_CV(DS_S_O2, the_classifier, the_feature_preprocessors);
    the_cross_validator_S_o2.num_resample_runs = numOfResample;
    the_cross_validator_S_o2.display_progress.resample_run_time = 0;
    the_cross_validator_S_o2.display_progress.zero_one_loss = 0;
    the_cross_validator_S_o2.test_only_at_training_times = 1;
    DECODING_RESULTS_S_o2 = the_cross_validator_S_o2.run_cv_decoding;
    fullSpaceDecodingAcc_S_o2(:,perm) =...
        mean(DECODING_RESULTS_S_o2.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
    %%
    the_cross_validator_R_o1 =...
        standard_resample_CV(DS_R_O1, the_classifier, the_feature_preprocessors);
    the_cross_validator_R_o1.num_resample_runs = numOfResample;
    the_cross_validator_R_o1.display_progress.resample_run_time = 0;
    the_cross_validator_R_o1.display_progress.zero_one_loss = 0;
    the_cross_validator_R_o1.test_only_at_training_times = 1;
    DECODING_RESULTS_R_o1 = the_cross_validator_R_o1.run_cv_decoding;
    fullSpaceDecodingAcc_R_o1(:,perm) =...
        mean(DECODING_RESULTS_R_o1.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
    %%
    the_cross_validator_R_o2 =...
        standard_resample_CV(DS_R_O2, the_classifier, the_feature_preprocessors);
    the_cross_validator_R_o2.num_resample_runs = numOfResample;
    the_cross_validator_R_o2.display_progress.resample_run_time = 0;
    the_cross_validator_R_o2.display_progress.zero_one_loss = 0;
    the_cross_validator_R_o2.test_only_at_training_times = 1;
    DECODING_RESULTS_R_o2 = the_cross_validator_R_o2.run_cv_decoding;
    fullSpaceDecodingAcc_R_o2(:,perm) =...
        mean(DECODING_RESULTS_R_o2.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
    %%
    the_cross_validator_Rtrain_Stest_O1 =...
        standard_resample_CV(DS_R_O1, the_classifier, the_feature_preprocessors);
    the_cross_validator_Rtrain_Stest_O1.num_resample_runs = numOfResample;
    the_cross_validator_Rtrain_Stest_O1.use_different_train_and_test_set = 1;
    the_cross_validator_Rtrain_Stest_O1.datasource2 = DS_S_O1;
    the_cross_validator_Rtrain_Stest_O1.display_progress.resample_run_time = 0;
    the_cross_validator_Rtrain_Stest_O1.display_progress.zero_one_loss = 0;
    the_cross_validator_Rtrain_Stest_O1.test_only_at_training_times = 1;
    DECODING_RESULTS_Rtrain_Stest_o1 = the_cross_validator_Rtrain_Stest_O1.run_cv_decoding;        
    fullSpaceDecodingAcc_Rtrain_Stest_o1(:,perm) =...
        mean(DECODING_RESULTS_Rtrain_Stest_o1.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;

    %----%%----
    the_cross_validator_Rtrain_Stest_O2 =...
        standard_resample_CV(DS_R_O2, the_classifier, the_feature_preprocessors);
    the_cross_validator_Rtrain_Stest_O2.num_resample_runs = numOfResample;
    the_cross_validator_Rtrain_Stest_O2.use_different_train_and_test_set = 1;
    the_cross_validator_Rtrain_Stest_O2.datasource2 = DS_S_O2;
    the_cross_validator_Rtrain_Stest_O2.display_progress.resample_run_time = 0;
    the_cross_validator_Rtrain_Stest_O2.display_progress.zero_one_loss = 0;
    the_cross_validator_Rtrain_Stest_O2.test_only_at_training_times = 1;
    DECODING_RESULTS_Rtrain_Stest_O2 = the_cross_validator_Rtrain_Stest_O2.run_cv_decoding;          
    fullSpaceDecodingAcc_Rtrain_Stest_o2(:,perm) =...
        mean(DECODING_RESULTS_Rtrain_Stest_O2.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
    %%
    the_cross_validator_Strain_Rtest_o1 =...
        standard_resample_CV(DS_S_O1, the_classifier, the_feature_preprocessors);
    the_cross_validator_Strain_Rtest_o1.num_resample_runs = numOfResample;
    the_cross_validator_Strain_Rtest_o1.use_different_train_and_test_set = 1;
    the_cross_validator_Strain_Rtest_o1.datasource2 = DS_R_O1;
    the_cross_validator_Strain_Rtest_o1.display_progress.resample_run_time = 0;
    the_cross_validator_Strain_Rtest_o1.display_progress.zero_one_loss = 0;
    the_cross_validator_Strain_Rtest_o1.test_only_at_training_times = 1;
    DECODING_RESULTS_Strain_Rtest_o1 = the_cross_validator_Strain_Rtest_o1.run_cv_decoding;        
    fullSpaceDecodingAcc_Strain_Rtest_o1(:,perm) =...
        mean(DECODING_RESULTS_Strain_Rtest_o1.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;

    %----%%----
    the_cross_validator_Strain_Rtest_o2 =...
        standard_resample_CV(DS_S_O2, the_classifier, the_feature_preprocessors);
    the_cross_validator_Strain_Rtest_o2.num_resample_runs = numOfResample;
    the_cross_validator_Strain_Rtest_o2.use_different_train_and_test_set = 1;
    the_cross_validator_Strain_Rtest_o2.datasource2 = DS_R_O2;
    the_cross_validator_Strain_Rtest_o2.display_progress.resample_run_time = 0;
    the_cross_validator_Strain_Rtest_o2.display_progress.zero_one_loss = 0;
    the_cross_validator_Strain_Rtest_o2.test_only_at_training_times = 1;
    DECODING_RESULTS_Strain_Rtest_o2 = the_cross_validator_Strain_Rtest_o2.run_cv_decoding;          
    fullSpaceDecodingAcc_Strain_Rtest_o2(:,perm) =....
        mean(DECODING_RESULTS_Strain_Rtest_o2.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
    %% in PC 2
    totalSpace = score(:,2:3);
    %%
    binedDataStype_o1 = cell(1,size(totalSpace,2));
    binedLabeStype_o1 = cell(1,size(totalSpace,2));
    for ni = 1 : size(totalSpace,2)
        binedDataStype_o1{ni} =...
            totalSpace(fSpaceLbels.OdorId == 1 &...
                       fSpaceLbels.RespLabel == 1, ni);
        binedLabeStype_o1{ni} =...
            fSpaceLbels.ConcId(fSpaceLbels.OdorId == 1 &...
                               fSpaceLbels.RespLabel == 1);
    end
    DS_S_O1 = basic_DS(binedDataStype_o1,...
        binedLabeStype_o1, num_cv_splits);
    %%
    binedDataStype_o2 = cell(1,size(totalSpace,2));
    binedLabeStype_o2 = cell(1,size(totalSpace,2));
    for ni = 1 : size(totalSpace,2)
        binedDataStype_o2{ni} =...
            totalSpace(fSpaceLbels.OdorId == 2 &...
                       fSpaceLbels.RespLabel == 1, ni);
        binedLabeStype_o2{ni} =...
            fSpaceLbels.ConcId(fSpaceLbels.OdorId == 2 &...
                               fSpaceLbels.RespLabel == 1);
    end
    DS_S_O2 = basic_DS(binedDataStype_o2,...
        binedLabeStype_o2, num_cv_splits);
    %%
    binedDataRtype_o1 = cell(1,size(totalSpace,2));
    binedLabeRtype_o1 = cell(1,size(totalSpace,2));
    for ni = 1 : size(totalSpace,2)
        binedDataRtype_o1{ni} =...
            totalSpace(fSpaceLbels.OdorId == 1 &...
                       fSpaceLbels.RespLabel == 2, ni);
        binedLabeRtype_o1{ni} =...
            fSpaceLbels.ConcId(fSpaceLbels.OdorId == 1 &...
                               fSpaceLbels.RespLabel == 2);
    end
    DS_R_O1 = basic_DS(binedDataRtype_o1,...
        binedLabeRtype_o1, num_cv_splits);
    %%
    binedDataRtype_o2 = cell(1,size(totalSpace,2));
    binedLabeRtype_o2 = cell(1,size(totalSpace,2));
    for ni = 1 : size(totalSpace,2)
        binedDataRtype_o2{ni} =...
            totalSpace(fSpaceLbels.OdorId == 2 &...
                       fSpaceLbels.RespLabel == 2, ni);
        binedLabeRtype_o2{ni} =...
            fSpaceLbels.ConcId(fSpaceLbels.OdorId == 2 &...
                               fSpaceLbels.RespLabel == 2);
    end
    DS_R_O2 = basic_DS(binedDataRtype_o2,...
        binedLabeRtype_o2, num_cv_splits);
    %%
    the_cross_validator_S_o1 =...
        standard_resample_CV(DS_S_O1, the_classifier, the_feature_preprocessors);
    the_cross_validator_S_o1.num_resample_runs = numOfResample;
    the_cross_validator_S_o1.display_progress.resample_run_time = 0;
    the_cross_validator_S_o1.display_progress.zero_one_loss = 0;
    the_cross_validator_S_o1.test_only_at_training_times = 1;
    DECODING_RESULTS_s_o1 = the_cross_validator_S_o1.run_cv_decoding;
    PC_2_and_3_DecodingAcc_S_o1(:,perm) =...
        mean(DECODING_RESULTS_s_o1.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
    %%
    the_cross_validator_S_o2 =...
        standard_resample_CV(DS_S_O2, the_classifier, the_feature_preprocessors);
    the_cross_validator_S_o2.num_resample_runs = numOfResample;
    the_cross_validator_S_o2.display_progress.resample_run_time = 0;
    the_cross_validator_S_o2.display_progress.zero_one_loss = 0;
    the_cross_validator_S_o2.test_only_at_training_times = 1;
    DECODING_RESULTS_S_o2 = the_cross_validator_S_o2.run_cv_decoding;
    PC_2_and_3_DecodingAcc_S_o2(:,perm) =...
        mean(DECODING_RESULTS_S_o2.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
    %%
    the_cross_validator_R_o1 =...
        standard_resample_CV(DS_R_O1, the_classifier, the_feature_preprocessors);
    the_cross_validator_R_o1.num_resample_runs = numOfResample;
    the_cross_validator_R_o1.display_progress.resample_run_time = 0;
    the_cross_validator_R_o1.display_progress.zero_one_loss = 0;
    the_cross_validator_R_o1.test_only_at_training_times = 1;
    DECODING_RESULTS_R_o1 = the_cross_validator_R_o1.run_cv_decoding;
    PC_2_and_3_DecodingAcc_R_o1(:,perm) =...
        mean(DECODING_RESULTS_R_o1.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
    %%
    the_cross_validator_R_o2 =...
        standard_resample_CV(DS_R_O2, the_classifier, the_feature_preprocessors);
    the_cross_validator_R_o2.num_resample_runs = numOfResample;
    the_cross_validator_R_o2.display_progress.resample_run_time = 0;
    the_cross_validator_R_o2.display_progress.zero_one_loss = 0;
    the_cross_validator_R_o2.test_only_at_training_times = 1;
    DECODING_RESULTS_R_o2 = the_cross_validator_R_o2.run_cv_decoding;
    PC_2_and_3_DecodingAcc_R_o2(:,perm) =...
        mean(DECODING_RESULTS_R_o2.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
    %%
    the_cross_validator_Rtrain_Stest_O1 =...
        standard_resample_CV(DS_R_O1, the_classifier, the_feature_preprocessors);
    the_cross_validator_Rtrain_Stest_O1.num_resample_runs = numOfResample;
    the_cross_validator_Rtrain_Stest_O1.use_different_train_and_test_set = 1;
    the_cross_validator_Rtrain_Stest_O1.datasource2 = DS_S_O1;
    the_cross_validator_Rtrain_Stest_O1.display_progress.resample_run_time = 0;
    the_cross_validator_Rtrain_Stest_O1.display_progress.zero_one_loss = 0;
    the_cross_validator_Rtrain_Stest_O1.test_only_at_training_times = 1;
    DECODING_RESULTS_Rtrain_Stest_o1 = the_cross_validator_Rtrain_Stest_O1.run_cv_decoding;        
    PC_2_and_3_DecodingAcc_Rtrain_Stest_o1(:,perm) =...
        mean(DECODING_RESULTS_Rtrain_Stest_o1.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;

    %----%%----
    the_cross_validator_Rtrain_Stest_O2 =...
        standard_resample_CV(DS_R_O2, the_classifier, the_feature_preprocessors);
    the_cross_validator_Rtrain_Stest_O2.num_resample_runs = numOfResample;
    the_cross_validator_Rtrain_Stest_O2.use_different_train_and_test_set = 1;
    the_cross_validator_Rtrain_Stest_O2.datasource2 = DS_S_O2;
    the_cross_validator_Rtrain_Stest_O2.display_progress.resample_run_time = 0;
    the_cross_validator_Rtrain_Stest_O2.display_progress.zero_one_loss = 0;
    the_cross_validator_Rtrain_Stest_O2.test_only_at_training_times = 1;
    DECODING_RESULTS_Rtrain_Stest_O2 = the_cross_validator_Rtrain_Stest_O2.run_cv_decoding;          
    PC_2_and_3_DecodingAcc_Rtrain_Stest_o2(:,perm) =...
        mean(DECODING_RESULTS_Rtrain_Stest_O2.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
    %%
    the_cross_validator_Strain_Rtest_o1 =...
        standard_resample_CV(DS_S_O1, the_classifier, the_feature_preprocessors);
    the_cross_validator_Strain_Rtest_o1.num_resample_runs = numOfResample;
    the_cross_validator_Strain_Rtest_o1.use_different_train_and_test_set = 1;
    the_cross_validator_Strain_Rtest_o1.datasource2 = DS_R_O1;
    the_cross_validator_Strain_Rtest_o1.display_progress.resample_run_time = 0;
    the_cross_validator_Strain_Rtest_o1.display_progress.zero_one_loss = 0;
    the_cross_validator_Strain_Rtest_o1.test_only_at_training_times = 1;
    DECODING_RESULTS_Strain_Rtest_o1 = the_cross_validator_Strain_Rtest_o1.run_cv_decoding;        
    PC_2_and_3_DecodingAcc_Strain_Rtest_o1(:,perm) =...
        mean(DECODING_RESULTS_Strain_Rtest_o1.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;

    %----%%----
    the_cross_validator_Strain_Rtest_o2 =...
        standard_resample_CV(DS_S_O2, the_classifier, the_feature_preprocessors);
    the_cross_validator_Strain_Rtest_o2.num_resample_runs = numOfResample;
    the_cross_validator_Strain_Rtest_o2.use_different_train_and_test_set = 1;
    the_cross_validator_Strain_Rtest_o2.datasource2 = DS_R_O2;
    the_cross_validator_Strain_Rtest_o2.display_progress.resample_run_time = 0;
    the_cross_validator_Strain_Rtest_o2.display_progress.zero_one_loss = 0;
    the_cross_validator_Strain_Rtest_o2.test_only_at_training_times = 1;
    DECODING_RESULTS_Strain_Rtest_o2 = the_cross_validator_Strain_Rtest_o2.run_cv_decoding;          
    PC_2_and_3_DecodingAcc_Strain_Rtest_o2(:,perm) =....
        mean(DECODING_RESULTS_Strain_Rtest_o2.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
end
%%
saveDir = fullfile(final_figs_path('beast1'),...
    'Fig5_cic_and_trans_decoding_of_conc_data','pooled');
if ~exist(saveDir , 'dir')
   mkdir(saveDir)
end
%%
save(fullfile(saveDir,'Fig5_cic_and_trans_decoding_of_conc_data.mat'),...
    'fullSpaceDecodingAcc_R_o1',...
    'fullSpaceDecodingAcc_S_o1',...
    'fullSpaceDecodingAcc_Rtrain_Stest_o1',...
    'fullSpaceDecodingAcc_Strain_Rtest_o1',...
    'fullSpaceDecodingAcc_R_o2',...
    'fullSpaceDecodingAcc_S_o2',...
    'fullSpaceDecodingAcc_Rtrain_Stest_o2',...
    'fullSpaceDecodingAcc_Strain_Rtest_o2',...
    'PC_2_and_3_DecodingAcc_R_o1',...
    'PC_2_and_3_DecodingAcc_S_o1',...
    'PC_2_and_3_DecodingAcc_Rtrain_Stest_o1',...
    'PC_2_and_3_DecodingAcc_Strain_Rtest_o1',...
    'PC_2_and_3_DecodingAcc_R_o2',...
    'PC_2_and_3_DecodingAcc_S_o2',...
    'PC_2_and_3_DecodingAcc_Rtrain_Stest_o2',...
    'PC_2_and_3_DecodingAcc_Strain_Rtest_o2',...
    '-v7.3');
