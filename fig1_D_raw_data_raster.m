clear; clc; close all;
[dataFilesDirList, mouseNameList, dataPathPrefixList] =....
    get_all_conc_mice('uniTn1');
%%
mn = 3;
oi = 51;
%for mn = 1: length(dataFilesDirList)
     %%
    mouseName = mouseNameList{mn};
    fprintf('Analysis %s data... \n', mouseName);
    mouseNameToPrint = mouseName(4:end);
    dataPathPrefix = dataPathPrefixList{mn};
    mouseDirName = fullfile(dataPathPrefix,'processedDataStorage',...
        mouseName);
    load(fullfile(mouseDirName,...
        sprintf('\\%s_respData.mat', mouseName)));
    load(fullfile(mouseDirName,...
        sprintf('\\%s_respLabelsData.mat', mouseName)));
    load(fullfile(mouseDirName,...
        sprintf('\\%s_inh_NR.mat', mouseName)));
    load(fullfile(mouseDirName,...
        sprintf('\\%s_recordingData.mat', mouseName)));
    load(fullfile(mouseDirName,...
        sprintf('\\%s_respTrace.mat', mouseName)));
    load(fullfile(mouseDirName,...
        sprintf('\\%s_spikesStract_sync.mat', mouseName)));
    load(fullfile(mouseDirName,...
        sprintf('\\%s_ADC_sync.mat', mouseName)));
    load(fullfile(mouseDirName,...
        sprintf('\\%s_eventOnsets.mat', mouseName)));
    %% get OCtc units
    [colorSet,  SandR_color]= set_plot_seting(15);
    plotDir = fullfile(final_figs_path('uniTn1'),...
        'fig1_B');
    if ~exist(plotDir , 'dir')
       mkdir(plotDir)
    end
    %%
    goodAndSafe_unitsFlag =...
        get_goodAndSafe_units(dataPathPrefix, mouseName,...
        6, .5);
    goodAndSafe_unitsId = find(goodAndSafe_unitsFlag);
    %% get units spikes
    selected_units = cell(length(goodAndSafe_unitsId), 1);
    singleUnitsIds = nan(length(goodAndSafe_unitsId), 1);
    singleUnitsFr = nan(length(goodAndSafe_unitsId), 1);
    for uIdx = 1 : length(goodAndSafe_unitsId)
        unitData.ui = spikesStract_sync.cids(uIdx);
        unitData.unitSt = spikesStract_sync.st(...
            spikesStract_sync.clu == unitData.ui);
        unitData.unitFr = length(unitData.unitSt)/...
            (unitData.unitSt(end) - unitData.unitSt(1));
        selected_units{uIdx} = unitData;
        clearvars unitData
    end
    
    %%
    if isfield(eventOnsets,'realTrialId')
        presentaionTrialIds = eventOnsets.realTrialId;
    else
        presentaionTrialIds = eventOnsets.trialId;
    end
    %% get raster
    before_onset_time = 10;
    after_onset_time = 10;
    odorwindow = 5;
%     for oi = 1: length(eventOnsets.odorOnset)
        if presentaionTrialIds(oi) > 2
            odor_onst = eventOnsets.odorOnset(oi);
            units_reaste =...
                zeros(length(selected_units),...
                (before_onset_time+....
                after_onset_time)*1000+1);
            for uIdx = 1 : length(selected_units)
                unitSt = selected_units{uIdx}.unitSt;
                spikes_times_in_rster_window =unitSt(...
                    (unitSt > (odor_onst-before_onset_time)) &...
                    (unitSt < (odor_onst+after_onset_time)))-...
                    (odor_onst-before_onset_time);
                units_reaste(uIdx,...
                    round(spikes_times_in_rster_window*1000)+1) = 1;
            end
            %%
            trimed_resp = resp(...
                (TimeStamps_ADC_sync > (odor_onst-before_onset_time))&...
                (TimeStamps_ADC_sync < (odor_onst+after_onset_time)));
            
            trimed_wind_xPeaks_flag =...
                (respData.rawXpeakPuffFree > (odor_onst-before_onset_time)*1000)&...
                (respData.rawXpeakPuffFree < (odor_onst+after_onset_time)*1000);
            trimed_resp_peaks =...
                respData.rawXpeakPuffFree(trimed_wind_xPeaks_flag)-...
                (odor_onst-before_onset_time)*1000;
            trimed_resp_labels =...
                respLabelsData.respLabelsRaw(trimed_wind_xPeaks_flag);
            %% --
            % [~, row_sorting_order] = sort(slopeVec);
            % ops.nCall = size(units_reaste);
            bin_size = 50;
            units_binned_reaste = nan(length(selected_units),...
                fix(size(units_reaste, 2)/bin_size));
            c = 1;
            for bt = 1 : bin_size :size(units_reaste, 2)-bin_size
                units_binned_reaste(:, c) =...
                    sum(units_reaste(:, bt: bt+50-1), 2);
                c = c +1;
            end
            %
            ops.iPC = 1:100;
            ops.upsamp = 100;
            [isort1, isort2, Sm] = mapTmap(units_binned_reaste, ops);
            %% --
            rstar_dot_size = .45;
            rstar_dot_colot = .4*[1,1,1];
            [row,col] =...
                  find(units_reaste(isort1,:) == 1);

            close all;
            figure
            subplot(5,1,1)
            hold on
            plot(trimed_resp, 'Color', .7*[1,1,1],...
                'LineWidth', 1)
            for i = 1 : length(trimed_resp_peaks)
                inhOnSet = max(find(trimed_resp(1:...
                    trimed_resp_peaks(i))>=0));
                if isempty(inhOnSet)
                    inhOnSet = 1;
                end
                inhOffSet = min(find(trimed_resp(...
                    trimed_resp_peaks(i): end)>=0)) + trimed_resp_peaks(i)-1;
                if isempty(inhOffSet)
                    inhOffSet = length(trimed_resp);
                end
                if trimed_resp_labels(i) == 1
                    plot(inhOnSet:inhOffSet,...
                        trimed_resp(inhOnSet: inhOffSet),...
                        'Color', SandR_color.s,...
                        'LineWidth', 1.1)
                elseif trimed_resp_labels(i) == 2
                    plot(inhOnSet:inhOffSet,...
                        trimed_resp(inhOnSet: inhOffSet),...
                        'Color', SandR_color.r,...
                        'LineWidth', 1.1)
                end
            end
            plot_y_lim = ylim;
            fill([before_onset_time,before_onset_time,...
                before_onset_time + odorwindow,...
                before_onset_time + odorwindow]*1000,...
                [plot_y_lim, plot_y_lim(2), plot_y_lim(1) ],...
                [.5, .5, .5],...
                'EdgeAlpha', 0, 'FaceAlpha', .15,...
                'FaceColor', colorSet(1,:))
            xlim([0, (before_onset_time+after_onset_time)*1000])
            yticks([])
            xticks([])

            box off

            subplot(5,1,[2:4])
            hold on
            scatter(col, row,rstar_dot_size,...
                'MarkerEdgeColor', rstar_dot_colot,...
                'MarkerFaceColor', rstar_dot_colot)
            plot_y_lim = ylim;
            xlim([0, (before_onset_time+after_onset_time)*1000])
            xticklabels(xticks/1000-before_onset_time)
            xlabel('time to odor onst (s)')
            yticks([])
            box off
            set(gcf,'Position',[100 100 900 400]);
            print_it(plotDir, sprintf('raw_data_raster_%d',oi), mouseName)
        end
%     end
%end