clear; clc; close all;

dataFilesDirList = {'M1\102721', 'M2\102821','M3\102921',...
    'M5\120221'};
mouseNameList = {'M1_102721', 'M2_102821','M3_102921',...
    'M5_120221'};
dataPathPrefix = 'E:\Ephys\sniffOdorProject\Conc_Series';
% dataPathPrefix = 'D:\Reza\sniffOdorProject\Conc_Series';
%%
tBeginning = .4;
tEnd = 5;
trialsToRemove = 2;
add_ndt_paths_and_init_rand_generator
numberOfBalanceSampling = 100;
%%

ptionsGlmnet= glmnetSet;
ptionsGlmnet.alpha = .95; 

optsGlmfit = statset('glmfit');
optsGlmfit.MaxIter = 1000;
numberOFShuffleing = 100;
%%
conc_set = [1,2,4];
baselineCond = 1;
%%
for mn = 1 : length(mouseNameList)
    %%
    mn
    mouseName = mouseNameList{mn};
    %%
    [fSpaceCell, fSpaceLbelsCell] =...
        get_aSample(dataPathPrefix, {mouseName},....
                    conc_set, 1:2, tBeginning, tEnd,trialsToRemove,...
                    1, baselineCond, 6, .5);
    totalSpace = fSpaceCell{1};
    %%
    regModelCoeff_results_o1 = nan(4, size(totalSpace,2),...
        numberOfBalanceSampling);
    regModelCoeff_results_o1_shuffled = nan(4, size(totalSpace,2),...
        numberOFShuffleing, numberOfBalanceSampling);
    noRegModelCoeff_results_o1 = nan(4, size(totalSpace,2),...
        numberOfBalanceSampling);
    
    regModelCoeff_results_o2 = nan(4, size(totalSpace,2),...
        numberOfBalanceSampling);
    regModelCoeff_results_o2_shuffled = nan(4, size(totalSpace,2),...
        numberOFShuffleing, numberOfBalanceSampling);
    noRegModelCoeff_results_o2 = nan(4, size(totalSpace,2),...
        numberOfBalanceSampling);
    %%
    [~, minSize] = ...
        get_minSize_mat(dataPathPrefix, {mouseName},....
            conc_set, 1:2,  tBeginning, tEnd, trialsToRemove);
    nFold = min([minSize-3, 10]);
    %%
    parfor perm = 1 : numberOfBalanceSampling 
        %%
        [fSpaceCell, fSpaceLbelsCell] =...
            get_aSample(dataPathPrefix, {mouseName},....
                        conc_set, 1:2, tBeginning, tEnd,trialsToRemove,...
                        1, baselineCond, 6, .5);
        %%
        totalSpace = fSpaceCell{1};
        fSpaceLbels = fSpaceLbelsCell{1};
        %%
        designMatO1 = zeros(length(fSpaceLbels.OdorId),3);
        designMatO1(fSpaceLbels.RespLabel==1, 1) = 1;
        designMatO1(fSpaceLbels.ConcId == 1, 2) = .33;
        designMatO1(fSpaceLbels.ConcId == 2, 2) = .66;
        designMatO1(fSpaceLbels.ConcId == 4, 2) = 1;
        designMatO1(fSpaceLbels.OdorId == 2, :) = []; 
        designMatO1(:, 3) = designMatO1(:, 1).* designMatO1(:, 2);
        %%
        regModelCoeffMAt = nan(4, size(totalSpace,2));
        noRegModelCoeffMAt = nan(4, size(totalSpace,2));
        regModelShuffle = nan(4, size(totalSpace,2),numberOFShuffleing);
        %%
        for ni = 1 : size(totalSpace,2)
            %%
            try
                respVecO1 = totalSpace(fSpaceLbels.OdorId ~= 2,ni);
   
                mdl_WI_reg = cvglmnet(designMatO1,respVecO1,....
                    'poisson',ptionsGlmnet,[],nFold);       
%                 mdl_WI_noReg =  fitglm(designMatO1,respVecO1,...
%                     'linear','Distribution','poisson');

                regModelCoeff = cvglmnetCoef(mdl_WI_reg, 'lambda_1se');
%                 noRegModelCoeff = mdl_WI_noReg.Coefficients.Estimate;

                regModelCoeffMAt(:,ni) = regModelCoeff;
%                 noRegModelCoeffMAt(:, ni) = noRegModelCoeff;
                
                for si = 1 : numberOFShuffleing
                    respVecO1Shuffle = respVecO1(randperm(length(respVecO1)));
                    mdl_WI_reg_shuufle = cvglmnet(designMatO1,respVecO1Shuffle,....
                        'poisson',ptionsGlmnet,[],nFold);    
                    regModelCoeffShuffle = cvglmnetCoef(mdl_WI_reg_shuufle, 'lambda_1se');
                    regModelShuffle(:, ni,si) = regModelCoeffShuffle;
                end
            end
            %%
        end
        regModelCoeff_results_o1(:,:,perm) = regModelCoeffMAt;
        noRegModelCoeff_results_o1(:,:,perm) = noRegModelCoeffMAt
        regModelCoeff_results_o1_shuffled(:,:,:,perm) =...
            regModelShuffle;
        %%
        %%
        designMatO2 = zeros(length(fSpaceLbels.OdorId),3);
        designMatO2(fSpaceLbels.RespLabel==1, 1) = 1;
        designMatO2(fSpaceLbels.ConcId == 1, 2) = .33;
        designMatO2(fSpaceLbels.ConcId == 2, 2) = .66;
        designMatO2(fSpaceLbels.ConcId == 4, 2) = 1;
        designMatO2(fSpaceLbels.OdorId == 1, :) = []; 
        designMatO2(:, 3) = designMatO2(:, 1).* designMatO2(:, 2);
        %%
        regModelCoeffMAt = nan(4, size(totalSpace,2));
        noRegModelCoeffMAt = nan(4, size(totalSpace,2));
        regModelShuffle = nan(4, size(totalSpace,2),...
                                numberOFShuffleing);
        %%
        
        %%
        for ni = 1 : size(totalSpace,2)
            %%
            try
                respVecO2 = totalSpace(fSpaceLbels.OdorId ~= 1,ni);
                
                mdl_WI_reg = cvglmnet(designMatO2,respVecO2,....
                    'poisson',ptionsGlmnet,[],nFold);      
                
%                 mdl_WI_noReg =  fitglm(designMatO2,respVecO2,...
%                     'linear','Distribution','poisson');

                regModelCoeff = cvglmnetCoef(mdl_WI_reg, 'lambda_1se');
%                 noRegModelCoeff = mdl_WI_noReg.Coefficients.Estimate;

                regModelCoeffMAt(:,ni) = regModelCoeff;
%                 noRegModelCoeffMAt(:, ni) = noRegModelCoeff;
                
                for si = 1 : numberOFShuffleing
                    respVecO2Shuffle = respVecO2(randperm(length(respVecO2)));
                    mdl_WI_reg_shuufle = cvglmnet(designMatO2,respVecO2Shuffle,....
                        'poisson',ptionsGlmnet,[],nFold);
                    regModelCoeffShuffle = cvglmnetCoef(mdl_WI_reg_shuufle, 'lambda_1se');
                    regModelShuffle(:, ni,si) = regModelCoeffShuffle;
                end
            end
            %%
        end
        regModelCoeff_results_o2(:,:,perm) = regModelCoeffMAt;
        noRegModelCoeff_results_o2(:,:,perm) = noRegModelCoeffMAt
        regModelCoeff_results_o2_shuffled(:,:,:,perm) = regModelShuffle;
    end
    glmModel.regModelCoeff_results_o2 = regModelCoeff_results_o2;
    glmModel.noRegModelCoeff_results_o2 = noRegModelCoeff_results_o2;
    glmModel.regModelCoeff_results_o2_shuffled =...
        regModelCoeff_results_o2_shuffled;
    
    glmModel.regModelCoeff_results_o1 = regModelCoeff_results_o1;
    glmModel.noRegModelCoeff_results_o1 = noRegModelCoeff_results_o1;
    glmModel.regModelCoeff_results_o1_shuffled =...
        regModelCoeff_results_o1_shuffled;
    glmModel.numberOFShuffleing =...
        numberOFShuffleing;
    glmModel.numberOfBalanceSampling =...
        numberOfBalanceSampling;
    %%
    mouseDirName = fullfile(dataPathPrefix,'processedDataStorage',...
        mouseName);
    save(fullfile(mouseDirName,...
        sprintf('\\%s_glmModel_with_baseline.mat', mouseName)),...
        'glmModel', '-v7.3');
end