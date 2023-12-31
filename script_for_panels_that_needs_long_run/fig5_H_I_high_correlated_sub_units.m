clear; clc; close all;
%%
dataFilesDirList = {'M1\102721', 'M2\102821','M3\102921',...
    'M5\120221'};
mouseNameList = {'M1_102721', 'M2_102821','M3_102921',...
    'M5_120221'};
  dataPathPrefix = 'E:\Ephys\sniffOdorProject\Conc_Series';
%  dataPathPrefix = 'D:\Reza\sniffOdorProject\Conc_Series';
%dataPathPrefix = 'D:\sniffOdorProject\Conc_Series';
%%
linearity_Coeff_poold = [];
linearity_pVal_poold = [];
nonlinearity_Coeff_poold = [];
nonlinearity_pVal_poold = [];
%%
for mn = 1 : length(mouseNameList)
    %%
    mouseName = mouseNameList{mn};
    mouseDirName = fullfile(dataPathPrefix,'processedDataStorage',...
        mouseName);
    load(fullfile(mouseDirName,...
            sprintf('\\%s_lrModel_zscored_no_baseline.mat', mouseName)));
    goodAndSafe_unitsFlag =...
            get_goodAndSafe_units(dataPathPrefix, mouseName,...
            6, .5);
   %%
   linearity_Coeff_poold =...
       [linearity_Coeff_poold,...
       nanmean(lrModel.linearlity_model_coeff(...
       :, goodAndSafe_unitsFlag,:, :), 4)];
   linearity_pVal_poold =...
       [linearity_pVal_poold,...
       nanmean(lrModel.linearlity_model_pVal(...
       :, goodAndSafe_unitsFlag,:, :), 4)];
   
   nonlinearity_Coeff_poold = ...
       [nonlinearity_Coeff_poold,...
       nanmean(lrModel.nonlinearlity_model_coeff(...
       :, goodAndSafe_unitsFlag,:, :), 4)];
   nonlinearity_pVal_poold =...
       [nonlinearity_pVal_poold,...
       nanmean(lrModel.nonlinearlity_model_pVal(...
       :, goodAndSafe_unitsFlag,:, :), 4)];
   
end
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
angleRange = deg2rad(30);
pVal = .01;
% sigGamma = squeeze(nonlinearity_pVal_poold(1, :, :)) <= pVal;
% sigAlpha = squeeze(linearity_pVal_poold(2, :, :)) <= pVal;
% sigBeta = squeeze(linearity_pVal_poold(3, :, :)) <= pVal;
alpha = squeeze(linearity_Coeff_poold(2,:, :));
beta = squeeze(linearity_Coeff_poold(3,:, :));
gamma = squeeze(nonlinearity_Coeff_poold(1,:, :));

%%
numberOfOdors = 2;
numOfSpace = 1;
selectionsCondition = 3;

decodingAcc_Mix =...
    nan(numOfResample, numOfSpace,...
    numberOfOdors, selectionsCondition, numberOfBalanceSampling);
decodingAcc_R =...
    nan(numOfResample, numOfSpace, ...
    numberOfOdors, selectionsCondition, numberOfBalanceSampling);
decodingAcc_S =...
    nan(numOfResample, numOfSpace, ...
    numberOfOdors, selectionsCondition, numberOfBalanceSampling);
decodingAcc_Rtrain_Stest =...
    nan(numOfResample, numOfSpace, ...
    numberOfOdors, selectionsCondition, numberOfBalanceSampling);
decodingAcc_Strain_Rtest =...
    nan(numOfResample, numOfSpace, ...
    numberOfOdors, selectionsCondition, numberOfBalanceSampling);
%%
concSet = [1,2,4];
parfor perm = 1 : numberOfBalanceSampling
    %%
    [fSpaceCell, fSpaceLbelsCell, eventLabelsCell] =...
        get_aSample(dataPathPrefix, mouseNameList,....
                    concSet, 1:2, tBeginning, tEnd,trialsToRemove,...
                    1, 0, 6, .5);
    %%
    pseudoPopulation =...
        [fSpaceCell{1}, fSpaceCell{2},fSpaceCell{3}, fSpaceCell{4}];
    fSpaceLbels = fSpaceLbelsCell{1};
    %%
    unitIdsCell = cell(selectionsCondition,numberOfOdors);
    for oId = 1 : numberOfOdors
        
        controlSelection =...
            find((abs(beta(:, oId)) <= abs(alpha(:, oId)*tan(pi/4+angleRange))) &...
            (abs(beta(:, oId)) >= abs(alpha(:, oId)*(pi/4-angleRange))));
        
        selectedUnits_neg_corr =...
            find((abs(beta(:, oId)) <= abs(alpha(:, oId)*tan(pi/4+angleRange))) &...
            (abs(beta(:, oId)) >= abs(alpha(:, oId)*(pi/4-angleRange))) &...
            (sign(alpha(:, oId)) ~= sign(beta(:, oId))));

        selectedUnits_poss_corr =...
            find((abs(beta(:, oId)) <= abs(alpha(:, oId)*tan(pi/4+angleRange))) &...
            (abs(beta(:, oId)) >= abs(alpha(:, oId)*(pi/4-angleRange))) &...
            (sign(alpha(:, oId)) == sign(beta(:, oId))));
        
        min_unit_size = min([length(selectedUnits_neg_corr),...
            length(selectedUnits_poss_corr)]);
        
        
        unitIdsCell{1, oId} = sort(selectedUnits_neg_corr(...
            randperm(size(selectedUnits_neg_corr,1),...
            min_unit_size)));                  
        unitIdsCell{2, oId} = sort(selectedUnits_poss_corr(...
            randperm(size(selectedUnits_poss_corr,1),...
            min_unit_size)));       
        unitIdsCell{3, oId} = sort(controlSelection(...
            randperm(size(controlSelection,1),...
            min_unit_size)));
    end
    %%
    decodingAcc_Mix_l3 =...
        nan(numOfResample, numOfSpace,...
        numberOfOdors, selectionsCondition);
    decodingAcc_R_l3 =...
        nan(numOfResample, numOfSpace, ...
        numberOfOdors, selectionsCondition);
    decodingAcc_S_l3 =...
        nan(numOfResample, numOfSpace, ...
        numberOfOdors, selectionsCondition);
    decodingAcc_Rtrain_Stest_l3 =...
        nan(numOfResample, numOfSpace, ...
        numberOfOdors, selectionsCondition);
    decodingAcc_Strain_Rtest_l3 =...
        nan(numOfResample, numOfSpace, ...
        numberOfOdors, selectionsCondition);
    %===
    for sc = 1: selectionsCondition
        decodingAcc_Mix_l2 =...
            nan(numOfResample, numOfSpace, numberOfOdors);
        decodingAcc_R_l2 =...
            nan(numOfResample, numOfSpace, numberOfOdors);
        decodingAcc_S_l2 =...
            nan(numOfResample, numOfSpace, numberOfOdors);
        decodingAcc_Rtrain_Stest_l2 =...
            nan(numOfResample, numOfSpace, numberOfOdors);
        decodingAcc_Strain_Rtest_l2 =...
            nan(numOfResample, numOfSpace, numberOfOdors);
        %----
        for oId = 1 : numberOfOdors
            decodingAcc_Mix_l1 =...
                nan(numOfResample, numOfSpace)
            decodingAcc_R_l1 =...
                nan(numOfResample, numOfSpace)
            decodingAcc_S_l1 =...
                nan(numOfResample, numOfSpace)
            decodingAcc_Rtrain_Stest_l1 =...
                nan(numOfResample, numOfSpace)
            decodingAcc_Strain_Rtest_l1 =...
                nan(numOfResample, numOfSpace)
            for sn = 1 : numOfSpace
                %----
                if sn == 1
                    totalSpace =...
                        pseudoPopulation(:, unitIdsCell{sc, oId});
                elseif sn == 2
                    totalSpace =...
                        pseudoPopulation(:, unitIdsCell{sc, oId});
                    [~,score,~] = pca(zscore(totalSpace));
                    if rem(sc, 2) == 1
                        totalSpace = score(:,1:2);
                    else 
                        totalSpace = score(:,2:3);
                    end
                end 
                %---------
%                 sample_size = fix(sum(...
%                     fSpaceLbels.OdorId == oId &...
%                     fSpaceLbels.RespLabel == 2 &...
%                     fSpaceLbels.ConcId == 1)/2);
%                 mix_sample = [];
%                 for ii = 1 : 2
%                     for cii = 1 :length(concSet)
%                         cc = concSet(cii);
%                         evets_sample =...
%                             find(fSpaceLbels.OdorId == oId &...
%                             fSpaceLbels.RespLabel == ii &...
%                             fSpaceLbels.ConcId == cc);
%                         mix_sample = [mix_sample;...
%                             datasample(evets_sample,...
%                             sample_size, 'Replace', false)];
% 
%                     end
%                 end
%                 mix_sample = sort(mix_sample);

%                 DS_mix = get_ndt_basicDs(...
%                     totalSpace(mix_sample, :),...
%                     fSpaceLbels.ConcId(mix_sample),...
%                     num_cv_splits)
                DS_S = get_ndt_basicDs(...
                    totalSpace(fSpaceLbels.OdorId == oId &...
                    fSpaceLbels.RespLabel == 1, :),...
                    fSpaceLbels.ConcId(fSpaceLbels.OdorId == oId &...
                    fSpaceLbels.RespLabel == 1), num_cv_splits)
                DS_R = get_ndt_basicDs(...
                    totalSpace(fSpaceLbels.OdorId == oId &...
                    fSpaceLbels.RespLabel == 2, :),...
                    fSpaceLbels.ConcId(fSpaceLbels.OdorId == oId &...
                    fSpaceLbels.RespLabel == 2), num_cv_splits)
                %--------
%                 the_cross_validator_Mix =...
%                 standard_resample_CV(DS_mix, the_classifier, the_feature_preprocessors);
%                 the_cross_validator_Mix.num_resample_runs = numOfResample;
%                 the_cross_validator_Mix.display_progress.resample_run_time = 0;
%                 the_cross_validator_Mix.display_progress.zero_one_loss = 0;
%                 the_cross_validator_Mix.test_only_at_training_times = 1;
%                 DECODING_RESULTS_Mix = the_cross_validator_Mix.run_cv_decoding;
%                 decodingAcc_Mix_l1(:,sn) =...
%                     mean(DECODING_RESULTS_Mix.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
%                 
%                 the_cross_validator_S =...
%                 standard_resample_CV(DS_S, the_classifier, the_feature_preprocessors);
%                 the_cross_validator_S.num_resample_runs = numOfResample;
%                 the_cross_validator_S.display_progress.resample_run_time = 0;
%                 the_cross_validator_S.display_progress.zero_one_loss = 0;
%                 the_cross_validator_S.test_only_at_training_times = 1;
%                 DECODING_RESULTS_s = the_cross_validator_S.run_cv_decoding;
%                 decodingAcc_S_l1(:,sn) =...
%                     mean(DECODING_RESULTS_s.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
% 
%                 the_cross_validator_R =...
%                     standard_resample_CV(DS_R, the_classifier, the_feature_preprocessors);
%                 the_cross_validator_R.num_resample_runs = numOfResample;
%                 the_cross_validator_R.display_progress.resample_run_time = 0;
%                 the_cross_validator_R.display_progress.zero_one_loss = 0;
%                 the_cross_validator_R.test_only_at_training_times = 1;
%                 DECODING_RESULTS_R = the_cross_validator_R.run_cv_decoding;
%                 decodingAcc_R_l1(:,sn) =...
%                     mean(DECODING_RESULTS_R.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
                
%                 the_cross_validator_Rtrain_Stest =...
%                     standard_resample_CV(DS_R, the_classifier, the_feature_preprocessors);
%                 the_cross_validator_Rtrain_Stest.num_resample_runs = numOfResample;
%                 the_cross_validator_Rtrain_Stest.use_different_train_and_test_set = 1;
%                 the_cross_validator_Rtrain_Stest.datasource2 = DS_S;
%                 the_cross_validator_Rtrain_Stest.display_progress.resample_run_time = 0;
%                 the_cross_validator_Rtrain_Stest.display_progress.zero_one_loss = 0;
%                 the_cross_validator_Rtrain_Stest.test_only_at_training_times = 1;
%                 DECODING_RESULTS_Rtrain_Stest = the_cross_validator_Rtrain_Stest.run_cv_decoding;        
%                 decodingAcc_Rtrain_Stest_l1(:,sn) =...
%                     mean(DECODING_RESULTS_Rtrain_Stest.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;

                %----%%----
                the_cross_validator_Strain_Rtest =...
                    standard_resample_CV(DS_S, the_classifier, the_feature_preprocessors);
                the_cross_validator_Strain_Rtest.num_resample_runs = numOfResample;
                the_cross_validator_Strain_Rtest.use_different_train_and_test_set = 1;
                the_cross_validator_Strain_Rtest.datasource2 = DS_R;
                the_cross_validator_Strain_Rtest.display_progress.resample_run_time = 0;
                the_cross_validator_Strain_Rtest.display_progress.zero_one_loss = 0;
                the_cross_validator_Strain_Rtest.test_only_at_training_times = 1;
                DECODING_RESULTS_Strain_Rtest = the_cross_validator_Strain_Rtest.run_cv_decoding;        
                decodingAcc_Strain_Rtest_l1(:,sn) =...
                    mean(DECODING_RESULTS_Strain_Rtest.ZERO_ONE_LOSS_RESULTS.decoding_results,2)*100;
                %-------
            end
            decodingAcc_Mix_l2(:,:,oId) =...
                decodingAcc_Mix_l1;
            decodingAcc_R_l2(:,:,oId) =...
                decodingAcc_R_l1;
            decodingAcc_S_l2(:,:,oId) =...
                decodingAcc_S_l1;
            decodingAcc_Rtrain_Stest_l2(:,:,oId) =...
                decodingAcc_Rtrain_Stest_l1;
            decodingAcc_Strain_Rtest_l2(:,:,oId) =...
                decodingAcc_Strain_Rtest_l1;
        end
        decodingAcc_Mix_l3(:,:,:,sc) =...
            decodingAcc_Mix_l2;
        decodingAcc_R_l3(:,:,:,sc) =...
            decodingAcc_R_l2;
        decodingAcc_S_l3(:,:,:,sc) =...
            decodingAcc_S_l2;
        decodingAcc_Rtrain_Stest_l3(:,:,:,sc) =...
            decodingAcc_Rtrain_Stest_l2;
        decodingAcc_Strain_Rtest_l3(:,:,:,sc) =...
            decodingAcc_Strain_Rtest_l2;
    end
    decodingAcc_Mix(:,:,:,:, perm) =...
            decodingAcc_Mix_l3;
    decodingAcc_R(:,:,:,:, perm) =...
        decodingAcc_R_l3;
    decodingAcc_S(:,:,:,:, perm) =...
        decodingAcc_S_l3;
    decodingAcc_Rtrain_Stest(:,:,:,:, perm) =...
        decodingAcc_Rtrain_Stest_l3;
    decodingAcc_Strain_Rtest(:,:,:,:, perm) =...
        decodingAcc_Strain_Rtest_l3;
end
%%
[~, SandR_colors] = set_plot_seting(17, 8);
save_dir = fullfile(final_figs_path('beast1'),...
    'fig5_H_I_high_correlated_sub_units_data_same_size_fixation');
if ~exist(save_dir , 'dir')
   mkdir(save_dir)
end
%%
save(fullfile(save_dir, 'fig5_H_I_high_correlated_sub_units_data.mat'),...
    'decodingAcc_Mix',...
    'decodingAcc_R',...
    'decodingAcc_S', ...
    'decodingAcc_Rtrain_Stest',...
    'decodingAcc_Strain_Rtest',...
    '-v7.3')