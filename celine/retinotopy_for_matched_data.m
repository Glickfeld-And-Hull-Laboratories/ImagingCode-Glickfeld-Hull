function [ret_npSub_tc_matched, ret_distance_matched, resp_by_stim_matched, ret_dfof_trial_matched, ...
          trialIndSourceUsed, lbub_fits_matched, goodfit_ind_matched, r2_vec_matched, ...
          fitAzimDeg_matched, fitElevDeg_matched, prefAzimDeg_matched, prefElevDeg_matched, ...
          Azs_matched, Els_matched] = ...
    retinotopy_for_matched_data(nd, allDays, expt, mouse, fov_avg, masks, fitGeoTAf, instructions, inputStructure, match_ind, validation_choice, fnOut)

ret_npSub_tc_matched   = cell(1,nd);
ret_distance_matched   = cell(1,nd);
resp_by_stim_matched   = cell(1,nd);
ret_dfof_trial_matched = cell(1,nd);
lbub_fits_matched      = cell(1,nd);
goodfit_ind_matched    = cell(1,nd);
r2_vec_matched         = cell(1,nd);
fitAzimDeg_matched     = cell(1,nd);
fitElevDeg_matched     = cell(1,nd);
prefAzimDeg_matched    = cell(1,nd);
prefElevDeg_matched    = cell(1,nd);
Azs_matched            = cell(1,nd);
Els_matched            = cell(1,nd);

for id = 1:nd
    sess = allDays(id);
    date = expt(sess).date;
    subnum = expt(sess).mouse;
    ImgFolder = expt(sess).ret_run;
    imgMatFile = [ImgFolder '_000_000.mat'];
    time = expt(sess).ret_time;

    retDataPath = fullfile('Z:\home\ACh\Data\2p_data',mouse,date,ImgFolder);
    load(fullfile(retDataPath,imgMatFile));
    nframes = info.config.frames;
    ret_data_temp = squeeze(sbxread(fullfile(retDataPath,[ImgFolder '_000_000']),0,nframes));
    fprintf(['Loaded ' num2str(nframes) ' frames \r\n'])

    fName = ['Z:\Behavior\Data\data-' subnum '-' date '-' time '.mat'];
    loadedData = load(fName);
    ret_inputStructure = loadedData.input;
    clear input

    referenceFOV = fov_avg{id};
    retAvrg = mean(ret_data_temp,3);
    [~, retAvrg_registered] = stackRegGPU(retAvrg,referenceFOV);
    [~, ret_data_registered] = stackRegGPU(ret_data_temp,retAvrg_registered);

    if id == 2
        ret_data_registered = imwarp(ret_data_registered,fitGeoTAf, 'OutputView', imref2d(size(ret_data_registered)));
    end
    ret_data_registered_FOV = mean(ret_data_registered,3);
    referenceMasks = masks{id};

    figure; imagesc(ret_data_registered_FOV), colormap gray; caxis([200 4000]);
    hold on
    bound = cell2mat(bwboundaries(referenceMasks(:,:,1)));
    plot(bound(:,2),bound(:,1),'.','color','b','MarkerSize',.1); hold on;
    drawnow

    ret_data_tc = stackGetTimeCourses(ret_data_registered, referenceMasks);
    mask_np = imCellNeuropil(referenceMasks, 3, 5);

    nCells = size(ret_data_tc,2);
    data_tc_down = stackGetTimeCourses(stackGroupProject(ret_data_registered,5), referenceMasks);
    sz = size(ret_data_registered);
    down = 5;
    data_reg_down = stackGroupProject(ret_data_registered,down);
    np_tc = zeros(sz(3),nCells);
    np_tc_down = zeros(floor(sz(3)./down), nCells);
    for i = 1:nCells
        np_tc(:,i) = stackGetTimeCourses(ret_data_registered,mask_np(:,:,i));
        np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
        fprintf(['     Cell #' num2str(i) '%s/n'])
    end

    ii = 0.01:0.01:1;
    x = zeros(length(ii), nCells);
    for i = 1:100
        x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
    end
    [~, ind] = max(x,[],1);
    np_w = 0.01*ind;

    ret_npSub_tc = ret_data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);

    if id==1
        ret_npSub_tc=ret_npSub_tc(:,match_ind);
    end
    nMatch = length(match_ind);

    nOn = ret_inputStructure(1).nScansOn;
    nOff = ret_inputStructure(1).nScansOff;

    switch instructions.tIdxSource
        case 'PD'
            [ret_stimOns_temp, ret_stimOffs] = photoFrameFinder_Sanworks(info.frame);
            ret_inputStructure.stimOns_photodiode = ret_stimOns_temp;
            ret_inputStructure.stimOffs_photodiode = ret_stimOffs;
            ret_inputStructure.stimOns_mwCounter = [];
            ret_inputStructure.stimTimingSource = 'PD';
            trialIndSourceUsed = 'PD';
        case 'MW'
            input_correct = counterValCorrect_noPhotodiode(input);
            ret_inputStructure.stimOns_mwCounter = cell2mat(input_correct.cStimOn);
            ret_inputStructure.counterValues = input_correct.counterValues;
            ret_inputStructure.counterTimesUs = input_correct.counterTimesUs;
            ret_inputStructure.stimOns_photodiode = [];
            ret_inputStructure.stimOffs_photodiode = [];
            ret_inputStructure.stimTimingSource = 'MW';
            ret_stimOns_temp = ret_inputStructure.stimOns_mwCounter;
            trialIndSourceUsed = 'MW';
            clear input_correct
        case 'cS'
            ret_inputStructure.stimTimingSource = 'cS';
            ret_stimOns_temp = cell2mat(ret_inputStructure.cStimOn);
            trialIndSourceUsed = 'cS';
        otherwise
            error('No valid trial indexing source specificed in instr file.');
    end

    nTrials = length(ret_stimOns_temp);
    ret_data_trial = nan(nOn+nOff,nMatch,nTrials);
    for itrial = 1:nTrials
        if ~isnan(ret_stimOns_temp(itrial)) & (ret_stimOns_temp(itrial)+nOn+nOff/2)<size(ret_npSub_tc,1)
            ret_data_trial(:,:,itrial) = ret_npSub_tc((ret_stimOns_temp(itrial)-nOff/2:ret_stimOns_temp(itrial)-1+nOn+nOff/2),:);
        end
    end

    baselineFrames = (nOff/4+1):nOff/2;
    stimFrames = (nOff/2+1):(nOff/2+nOn);

    F0 = mean(ret_data_trial(baselineFrames,:,:), 1);
    ret_dfof_trial = (ret_data_trial - F0) ./ F0;

    ret_mean_resp = squeeze(mean(ret_dfof_trial(stimFrames,:,:), 1));  % nMatch x nTrials

    trialAz = celleqel2mat_padded(ret_inputStructure.tGratingAzimuthDeg);
    azs = unique(trialAz);
    trialEl = celleqel2mat_padded(ret_inputStructure.tGratingElevationDeg);
    els = unique(trialEl);
    if min(els) < 0, els = fliplr(els); end

    resp_by_stim = nan(length(els),length(azs),nMatch);
    for i_el = 1:length(els)
        this_el_trials = find(trialEl == els(i_el));
        for i_az = 1:length(azs)
            this_az_trials = find(trialAz==azs(i_az));
            these_trials = intersect(this_el_trials, this_az_trials);
            resp_by_stim(i_el, i_az,:) = nanmean(ret_mean_resp(:,these_trials),2);
        end
    end

    [nElev, nAzim, nMatch] = size(resp_by_stim);

    % Keep max for validation plot highlighting
    resp_reshaped = reshape(resp_by_stim, [], nMatch);
    [~, maxIdx] = max(resp_reshaped, [], 1);
    [maxElev, maxAzim] = ind2sub([nElev, nAzim], maxIdx);

    % Weighted average RF position
    [azGrid, elGrid] = meshgrid(azs, els);
    prefAzimDeg = nan(1, nMatch);
    prefElevDeg = nan(1, nMatch);
    for iCell = 1:nMatch
        w = max(resp_by_stim(:,:,iCell), 0);
        total_w = sum(w(:));
        if total_w > 0
            prefAzimDeg(iCell) = sum(w(:) .* azGrid(:)) / total_w;
            prefElevDeg(iCell) = sum(w(:) .* elGrid(:)) / total_w;
        end
    end
    nNaN = sum(isnan(prefAzimDeg));
    fprintf('Day %d: %d/%d cells have NaN preferred position (all responses <= 0)\n', id, nNaN, nMatch);

    % Validation trace plots
    if validation_choice && nMatch >= 5
        rng('shuffle');
        selected_cells = randperm(nMatch, min(5, nMatch));
        for i_cell = 1:length(selected_cells)
            this_cell = selected_cells(i_cell);
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
            figure('Name', sprintf('Day %d - Cell %d', id, this_cell));
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
                        hold on;
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
                    if i_el == length(els), xlabel('Frame'); end
                    if i_az == 1, ylabel('dF/F'); end
                end
            end
            drawnow
        end
    end

    %% 2D Gaussian fitting with shuffling
    nStim = length(azs) * length(els);

    % Build Stims and Ind_struct
    Stims = zeros(nStim, 2);
    idx = 1;
    for iEl = 1:length(els)
        for iAz = 1:length(azs)
            Stims(idx,:) = [els(iEl) azs(iAz)];
            idx = idx + 1;
        end
    end

    Ind_struct = [];
    for iStim = 1:nStim
        Ind_struct(iStim).all_trials = intersect(find(trialEl == Stims(iStim,1)), find(trialAz == Stims(iStim,2)));
    end

    % Grid for fitting scripts
    grid2.AzAz = azGrid;
    grid2.ElEl = elGrid;
    dAz = median(diff(azs));
    dEl = median(diff(els));
    Az_vec00 = azs(1):(dAz/10):azs(end);
    El_vec00 = els(1):(dEl/10):els(end);
    [AzAz00, ElEl00] = meshgrid(Az_vec00, El_vec00);
    grid2.AzAz00 = AzAz00;
    grid2.ElEl00 = ElEl00;

    Nshuf = 100;
    Fit_struct = struct('True', cell(1, nMatch), 'Shuf', cell(1, nMatch));
    resp_dFoverF = ret_mean_resp;  % nMatch x nTrials

    fprintf('Day %d: fitting retinotopy (Nshuf=%d)...\n', id, Nshuf)
   
for count_shuf = 1:Nshuf
    for iCell = 1:nMatch
        if ~isempty(Fit_struct(iCell).Shuf) && length(Fit_struct(iCell).Shuf) >= count_shuf
            shuf_entry = Fit_struct(iCell).Shuf(count_shuf);
            if isstruct(shuf_entry) && isfield(shuf_entry, 's_') && isstruct(shuf_entry.s_) && isfield(shuf_entry.s_, 'x')
                shuf_s = shuf_entry.s_;
                fit_shuf_vec(iCell,:,count_shuf) = [shuf_s.x, ...
                    shuf_s.Elhicut_50, shuf_s.Azhicut_50, ...
                    shuf_s.Elhicut_10, shuf_s.Azhicut_10];
            end
        end
    end
end
    % Extract fit parameter vectors
    fit_true_vec = NaN(nMatch, 10);
    for iCell = 1:nMatch
        if ~isempty(Fit_struct(iCell).True)
            fit_true_vec(iCell,:) = [Fit_struct(iCell).True.s_.x, ...
                Fit_struct(iCell).True.s_.Elhicut_50, Fit_struct(iCell).True.s_.Azhicut_50, ...
                Fit_struct(iCell).True.s_.Elhicut_10, Fit_struct(iCell).True.s_.Azhicut_10];
        end
    end

    fit_shuf_vec = NaN(nMatch, 10, Nshuf);
    for count_shuf = 1:Nshuf
        for iCell = 1:nMatch
            if ~isempty(Fit_struct(iCell).Shuf) && length(Fit_struct(iCell).Shuf) >= count_shuf && isfield(Fit_struct(iCell).Shuf(count_shuf), 's_')
                fit_shuf_vec(iCell,:,count_shuf) = [Fit_struct(iCell).Shuf(count_shuf).s_.x, ...
                    Fit_struct(iCell).Shuf(count_shuf).s_.Elhicut_50, Fit_struct(iCell).Shuf(count_shuf).s_.Azhicut_50, ...
                    Fit_struct(iCell).Shuf(count_shuf).s_.Elhicut_10, Fit_struct(iCell).Shuf(count_shuf).s_.Azhicut_10];
            end
        end
    end

    % lbub confidence intervals
    Npars       = size(fit_shuf_vec, 2);
    lbub_fits   = NaN(nMatch, Npars, 5);
    alpha_bound = 0.025;
    ind_shuf_lb = ceil(Nshuf * alpha_bound);
    ind_shuf_ub = ceil(Nshuf * (1-alpha_bound));

    for iCell = 1:nMatch
        for count2 = 1:Npars
            i_sorted = sort(squeeze(fit_shuf_vec(iCell,count2,:)));
            lbub_fits(iCell,count2,1) = i_sorted(ind_shuf_lb);
            lbub_fits(iCell,count2,2) = i_sorted(ind_shuf_ub);
            lbub_fits(iCell,count2,3) = mean(i_sorted);
            lbub_fits(iCell,count2,5) = std(i_sorted);
        end
        lbub_fits(iCell,:,4) = fit_true_vec(iCell,:);
    end

    lbub_diff = lbub_fits(:,:,2) - lbub_fits(:,:,1);

    % lbub goodfit: Az0 and El0 CI width < 2x step
    goodfit_ind = [];
    azStep = ret_inputStructure(1).gratingAzimuthStepDeg;
    for iCell = 1:nMatch
        if lbub_diff(iCell,4) < azStep*2 && lbub_diff(iCell,5) < azStep*2
            goodfit_ind = [goodfit_ind iCell];
        end
    end

    % Remove RFs at retinotopy perimeter
    goodfit_ind2 = zeros(size(goodfit_ind));
    for i = 1:length(goodfit_ind)
        if sum(round(lbub_fits(goodfit_ind(i),4,4)) == [min(azs) max(azs)]) || ...
           sum(round(lbub_fits(goodfit_ind(i),5,4)) == [min(els) max(els)])
            continue
        end
        goodfit_ind2(i) = goodfit_ind(i);
    end
    goodfit_ind = goodfit_ind2(goodfit_ind2 ~= 0);
    fprintf('Day %d: %d good-fit cells (lbub)\n', id, length(goodfit_ind))

    % R� goodness of fit
    r2_threshold = 0.65;
    r2_vec = nan(1, nMatch);
    for iCell = 1:nMatch
        if isempty(Fit_struct(iCell).True), continue; end
        x_fit = Fit_struct(iCell).True.s_.x;
        Az0 = x_fit(4); El0 = x_fit(5); xi = x_fit(6);
        sig_az = x_fit(2); sig_el = x_fit(3);
        az_rot =  (azGrid - Az0) .* cos(xi) + (elGrid - El0) .* sin(xi);
        el_rot = -(azGrid - Az0) .* sin(xi) + (elGrid - El0) .* cos(xi);
        fit_surf = x_fit(1) .* exp(-(az_rot.^2./(2*sig_az^2) + el_rot.^2./(2*sig_el^2)));
        actual = resp_by_stim(:,:,iCell);
        ss_tot = sum((actual(:) - mean(actual(:))).^2);
        if ss_tot > 0
            r2_vec(iCell) = 1 - sum((actual(:) - fit_surf(:)).^2) / ss_tot;
        end
    end
    goodfit_ind_r2 = find(r2_vec >= r2_threshold);
    fprintf('Day %d: %d good-fit cells (R2 >= %.2f)\n', id, length(goodfit_ind_r2), r2_threshold)

    fitAzimDeg = lbub_fits(:,4,4)';
    fitElevDeg = lbub_fits(:,5,4)';

    finalAzim = double(inputStructure(id).gratingAzimuthDeg);
    finalElev  = double(inputStructure(id).gratingElevationDeg);

    ret_distance_matched{id}   = sqrt((prefAzimDeg - finalAzim).^2 + (prefElevDeg - finalElev).^2);
    ret_npSub_tc_matched{id}   = ret_npSub_tc;
    resp_by_stim_matched{id}   = resp_by_stim;
    ret_dfof_trial_matched{id} = ret_dfof_trial;
    lbub_fits_matched{id}      = lbub_fits;
    goodfit_ind_matched{id}    = goodfit_ind;
    r2_vec_matched{id}         = r2_vec;
    fitAzimDeg_matched{id}     = fitAzimDeg;
    fitElevDeg_matched{id}     = fitElevDeg;
    prefAzimDeg_matched{id}    = prefAzimDeg;
    prefElevDeg_matched{id}    = prefElevDeg;
    Azs_matched{id}            = azs;
    Els_matched{id}            = els;
end

save(fullfile(fnOut, [mouse '_ret_matched.mat']), ...
    'ret_npSub_tc_matched', 'ret_distance_matched', 'resp_by_stim_matched', 'ret_dfof_trial_matched', ...
    'trialIndSourceUsed', 'lbub_fits_matched', 'goodfit_ind_matched', 'r2_vec_matched', ...
    'fitAzimDeg_matched', 'fitElevDeg_matched', 'prefAzimDeg_matched', 'prefElevDeg_matched', ...
    'Azs_matched', 'Els_matched', '-v7.3');
fprintf('Saved to %s\n', fnOut)

end
