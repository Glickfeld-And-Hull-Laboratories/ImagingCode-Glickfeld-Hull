function [ret_npSub_tc, ret_distance, resp_by_stim, ret_dfof_trial, trialIndSourceUsed] = ...
    retinotopy_for_singleday(day_id, expt, mouse, inputStructure, referenceFOV, mask_cell, timingSource, rc, validation_choice)

expDate   = expt(day_id).date;
ImgFolder = expt(day_id).ret_run;
time      = expt(day_id).ret_time;

retDataPath = fullfile(rc.data, mouse, expDate, ImgFolder);
load(fullfile(retDataPath, [ImgFolder '_000_000.mat']));
nframes = info.config.frames;
ret_data_temp = squeeze(sbxread(fullfile(retDataPath, [ImgFolder '_000_000']), 0, nframes));
fprintf('Loaded %d frames\n', nframes);

fName = ['Z:\Behavior\Data\data-' mouse '-' expDate '-' time '.mat'];
loadedData = load(fName);
ret_inputStructure = loadedData.input;
clear loadedData

retAvrg = mean(ret_data_temp, 3);
[~, retAvrg_registered] = stackRegGPU(retAvrg, referenceFOV);
[~, ret_data_registered] = stackRegGPU(ret_data_temp, retAvrg_registered);

figure; imagesc(mean(ret_data_registered, 3)); colormap gray;
hold on
bound = cell2mat(bwboundaries(mask_cell > 0));
plot(bound(:,2), bound(:,1), '.', 'color', 'b', 'MarkerSize', .1);
drawnow

ret_data_tc = stackGetTimeCourses(ret_data_registered, mask_cell);
mask_np     = imCellNeuropil(mask_cell, 3, 5);

nCells = size(ret_data_tc, 2);
sz     = size(ret_data_registered);
down   = 5;
data_reg_down = stackGroupProject(ret_data_registered, down);
np_tc      = zeros(sz(3), nCells);
np_tc_down = zeros(floor(sz(3)/down), nCells);

for i = 1:nCells
    np_tc(:,i)      = stackGetTimeCourses(ret_data_registered, mask_np(:,:,i));
    np_tc_down(:,i) = stackGetTimeCourses(data_reg_down, mask_np(:,:,i));
    fprintf('     Cell #%d\n', i);
end

data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
ii = 0.01:0.01:1;
x  = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down - tcRemoveDC(np_tc_down * ii(i)));
end
[~, ind] = max(x, [], 1);
np_w = 0.01 * ind;

ret_npSub_tc = ret_data_tc - bsxfun(@times, tcRemoveDC(np_tc), np_w);

switch timingSource
    case 'PD'
        [ret_stimOns, ~] = photoFrameFinder_Sanworks(info.frame);
        trialIndSourceUsed = 'PD';
    case 'MW'
        input_correct = counterValCorrect_noPhotodiode(ret_inputStructure);
        ret_stimOns = cell2mat(input_correct.cStimOn);
        trialIndSourceUsed = 'MW';
        clear input_correct
    case 'cS'
        ret_stimOns = cell2mat(ret_inputStructure.cStimOn);
        trialIndSourceUsed = 'cS';
    otherwise
        error('Unrecognised timing source: %s. Expected PD, MW, or cS.', timingSource);
end

nOn     = ret_inputStructure(1).nScansOn;
nOff    = ret_inputStructure(1).nScansOff;
nTrials = length(ret_stimOns);

ret_data_trial = nan(nOn + nOff, nCells, nTrials);
for iTrial = 1:nTrials
    if ~isnan(ret_stimOns(iTrial)) && ...
       (ret_stimOns(iTrial) - nOff/2) >= 1 && ...
       (ret_stimOns(iTrial) + nOn + nOff/2) <= size(ret_npSub_tc, 1)
        ret_data_trial(:,:,iTrial) = ret_npSub_tc(ret_stimOns(iTrial) - nOff/2 : ret_stimOns(iTrial) - 1 + nOn + nOff/2, :);
    end
end

baselineFrames = (nOff/4 + 1):nOff/2;
stimFrames     = (nOff/2 + 1):(nOff/2 + nOn);

F0 = mean(ret_data_trial(baselineFrames,:,:), 1);
ret_dfof_trial = (ret_data_trial - F0) ./ F0;

ret_mean_resp = squeeze(mean(ret_dfof_trial(stimFrames,:,:), 1));

trialAz = celleqel2mat_padded(ret_inputStructure.tGratingAzimuthDeg);
trialEl = celleqel2mat_padded(ret_inputStructure.tGratingElevationDeg);
azs = unique(trialAz);
els = flipud(unique(trialEl));

resp_by_stim = nan(length(els), length(azs), nCells);
for i_el = 1:length(els)
    this_el_trials = find(trialEl == els(i_el));
    for i_az = 1:length(azs)
        these_trials = intersect(this_el_trials, find(trialAz == azs(i_az)));
        resp_by_stim(i_el, i_az, :) = nanmean(ret_mean_resp(:, these_trials), 2);
    end
end

[nElev, nAzim, ~] = size(resp_by_stim);

% Keep max for validation plot highlighting
resp_reshaped = reshape(resp_by_stim, [], nCells);
[~, maxIdx]   = max(resp_reshaped, [], 1);
[maxElev, maxAzim] = ind2sub([nElev, nAzim], maxIdx);

% Weighted average RF position: weights = positive responses only
[azGrid, elGrid] = meshgrid(azs, els);
prefAzimDeg = nan(1, nCells);
prefElevDeg = nan(1, nCells);
for iCell = 1:nCells
    w = max(resp_by_stim(:,:,iCell), 0);
    total_w = sum(w(:));
    if total_w > 0
        prefAzimDeg(iCell) = sum(w(:) .* azGrid(:)) / total_w;
        prefElevDeg(iCell) = sum(w(:) .* elGrid(:)) / total_w;
    end
end
nNaN = sum(isnan(prefAzimDeg));
fprintf('%d/%d cells have NaN preferred position (all responses <= 0)\n', nNaN, nCells);

if validation_choice && nCells >= 1
    rng('shuffle');
    selected_cells = randperm(nCells, min(5, nCells));

    for i_cell = 1:length(selected_cells)
        this_cell  = selected_cells(i_cell);
        all_traces = [];
        for i_el = 1:length(els)
            for i_az = 1:length(azs)
                these_trials = intersect(find(trialEl == els(i_el)), find(trialAz == azs(i_az)));
                if ~isempty(these_trials)
                    mean_trace = squeeze(nanmean(ret_dfof_trial(:, this_cell, these_trials), 3));
                    all_traces = [all_traces; mean_trace(:)];
                end
            end
        end
        y_limits = [min(all_traces) max(all_traces)];

        figure('Name', sprintf('Cell %d', this_cell));
        for i_el = 1:length(els)
            for i_az = 1:length(azs)
                subplot(length(els), length(azs), (i_el-1)*length(azs) + i_az);
                these_trials = intersect(find(trialEl == els(i_el)), find(trialAz == azs(i_az)));
                is_max = (i_el == maxElev(this_cell)) && (i_az == maxAzim(this_cell));
                if ~isempty(these_trials)
                    mean_trace = squeeze(nanmean(ret_dfof_trial(:, this_cell, these_trials), 3));
                    if is_max
                        plot(mean_trace, 'r', 'LineWidth', 2);
                    else
                        plot(mean_trace, 'k', 'LineWidth', 1);
                    end
                    hold on
                    xline(nOff/2, 'r--');
                    xline(nOff/2 + nOn, 'r--');
                end
                ylim(y_limits);
                if is_max
                    title(sprintf('Az:%d El:%d *MAX*', azs(i_az), els(i_el)), 'Color', 'r', 'FontWeight', 'bold');
                else
                    title(sprintf('Az:%d El:%d', azs(i_az), els(i_el)));
                end
                set(gca, 'TickDir', 'out'); grid off; box off;
                if i_el == length(els),  xlabel('Frame'); end
                if i_az == 1,            ylabel('dF/F');  end
            end
        end
        if isnan(prefAzimDeg(this_cell))
            sgtitle(sprintf('Cell %d | Weighted pref: NaN (no positive responses)', this_cell));
        else
            sgtitle(sprintf('Cell %d | Weighted pref: Az=%.1f deg, El=%.1f deg', this_cell, prefAzimDeg(this_cell), prefElevDeg(this_cell)));
        end
        drawnow
    end
end

finalAzim = double(inputStructure.gratingAzimuthDeg);
finalElev = double(inputStructure.gratingElevationDeg);

ret_distance = sqrt((prefAzimDeg - finalAzim).^2 + (prefElevDeg - finalElev).^2);

end