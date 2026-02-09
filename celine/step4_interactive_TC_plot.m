% Convert raw calcium timecourses to trial-structured dF/F data
data_dfof_trial_match = cell(1, nd);
fractTimeActive_match = cell(1, nd);
cellstd_match = cell(1, nd);
for id = 1:nd
    cStimOnTemp = stimOns{id};
    nTrials(id) = length(cStimOnTemp);
    [nFrames, nCells] = size(cellTCs_match{id});
    data_trial_match = nan(nOn + nOff, nTrials(id), nCells);
    for iTrial = 1:nTrials(id)
        if ~isnan(cStimOnTemp(iTrial)) && (cStimOnTemp(iTrial) + nOn + nOff/2) <= nFrames && (cStimOnTemp(iTrial) - nOff/2) >= 1
            data_trial_match(:, iTrial, :) = cellTCs_match{id}(cStimOnTemp(iTrial) - nOff/2:cStimOnTemp(iTrial) - 1 + nOn + nOff/2, :);
        end
    end
    fractTimeActive_match{id} = zeros(1, nCells);
    data_f_match = mean(data_trial_match(1:(nOff/2), :, :), 1);
    data_dfof_trial_match{id} = bsxfun(@rdivide, bsxfun(@minus, data_trial_match, data_f_match), data_f_match);
    
    % Interactive plot
    figure;
    h = plot(squeeze(mean(data_dfof_trial_match{id}(:,:,:),2, 'omitmissing')));
    set(gca, 'TickDir', 'out'); grid off; box off;
    dcm = datacursormode(gcf);
    datacursormode on
    set(dcm, 'UpdateFcn', @(~,event) sprintf('Cell %d\nFrame: %d\ndF/F: %.3f', ...
        event.Target.SeriesIndex, event.Position(1), event.Position(2)));
    
    meansub_match = cellTCs_match{id} - nanmean(cellTCs_match{id}, 1);
    cellstd = nanstd(meansub_match, [], 1);
    cellstd_match{id} = cellstd;
    for iCell = 1:nCells
        fractTimeActive_match{id}(:, iCell) = length(find(meansub_match(:, iCell) > 3.*cellstd(1, iCell))) ./ nFrames;
    end
    clear data_trial_match data_f_match cellstd cStimOnTemp
end