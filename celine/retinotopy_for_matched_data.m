function [ret_npSub_tc_matched, ret_distance_matched, resp_by_stim_matched, ret_dfof_trial_matched, ...
          trialIndSourceUsed, lbub_fits_matched, goodfit_ind_matched, r2_vec_matched, ...
          fitAzimDeg_matched, fitElevDeg_matched, prefAzimDeg_matched, prefElevDeg_matched, ...
          Azs_matched, Els_matched, distMap_matched, dist_vec_matched] = ...
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
distMap_matched        = cell(1,nd);
dist_vec_matched       = cell(1,nd);

if computer == 'GLNXA64'
    isilonName = '/home/cc735@dhe.duke.edu/GlickfeldLabShare/All_staff';
else
    isilonName = 'Z:';
end

for id = 1:nd
    sess      = allDays(id);
    date      = expt(sess).date;
    subnum    = expt(sess).mouse;
    ImgFolder = expt(sess).ret_run;
    time      = expt(sess).ret_time;

    % Load imaging data
    retDataPath = fullfile(isilonName, '/home/ACh/Data/2p_data', mouse, date, ImgFolder);
    imgMatFile  = [ImgFolder '_000_000.mat'];
    load(fullfile(retDataPath, imgMatFile));
    nframes   = info.config.frames;
    fprintf('Day %d: reading %d frames...\n', id, nframes)
    data_temp = squeeze(sbxread(fullfile(retDataPath, [ImgFolder '_000_000']), 0, nframes));

    % Load behavior data
    fName = fullfile(isilonName, '/Behavior/Data', ['data-' subnum '-' date '-' time '.mat']);
    loadedData = load(fName);
    ret_inputStructure = loadedData.input;
    nOn  = ret_inputStructure.nScansOn;
    nOff = ret_inputStructure.nScansOff;

    % Find stimulus onsets
    switch instructions.tIdxSource
        case 'PD'
            [stimOns, ~] = photoFrameFinder_Sanworks(info.frame);
            trialIndSourceUsed = 'PD';
        case 'MW'
            input_correct = counterValCorrect_noPhotodiode(ret_inputStructure);
            stimOns = cell2mat(input_correct.cStimOn);
            trialIndSourceUsed = 'MW';
        case 'cS'
            stimOns = cell2mat(ret_inputStructure.cStimOn);
            trialIndSourceUsed = 'cS';
        otherwise
            error('No valid trial indexing source specified in instr file.');
    end
    nTrials = length(stimOns);
    fprintf('Day %d: %d trials detected\n', id, nTrials)

    % Register ret data to reference FOV
    referenceFOV   = fov_avg{id};
    referenceMasks = masks{id};
    nCells         = max(referenceMasks(:));

    fprintf('Day %d: registering to reference FOV...\n', id)
    retAvrg = mean(data_temp, 3);
    [~, retAvrg_reg] = stackRegGPU(retAvrg, referenceFOV);
    [~, data_reg]    = stackRegGPU(data_temp, retAvrg_reg);
    clear data_temp retAvrg retAvrg_reg

    % Apply geometric transform for matched day
    if id == 2
        data_reg = imwarp(data_reg, fitGeoTAf, 'OutputView', imref2d(size(data_reg)));
    end
    sz = size(data_reg);

    % Extract timecourses with neuropil subtraction
    fprintf('Day %d: extracting timecourses...\n', id)
    down          = 5;
    data_reg_down = stackGroupProject(data_reg, down);
    data_tc       = stackGetTimeCourses(data_reg, referenceMasks);
    data_tc_down  = stackGetTimeCourses(data_reg_down, referenceMasks);

    mask_np    = imCellNeuropil(referenceMasks, 3, 5);
    np_tc      = zeros(sz(3), nCells);
    np_tc_down = zeros(floor(sz(3)/down), nCells);
    for i = 1:nCells
        np_tc(:,i)      = stackGetTimeCourses(data_reg, mask_np(:,:,i));
        np_tc_down(:,i) = stackGetTimeCourses(data_reg_down, mask_np(:,:,i));
    end
    clear data_reg data_reg_down

    ii = 0.01:0.01:1;
    x  = zeros(100, nCells);
    for i = 1:100
        x(i,:) = skewness(data_tc_down - tcRemoveDC(np_tc_down * ii(i)));
    end
    [~, ind] = max(x, [], 1);
    np_w     = 0.01 * ind;
    npSub_tc = data_tc - tcRemoveDC(np_tc) .* np_w;
    clear data_tc data_tc_down np_tc np_tc_down

    % Subset to matched cells (day 1 only; day 2 already in matched space)
    if id == 1
        npSub_tc_match = npSub_tc(:, match_ind);
    else
        npSub_tc_match = npSub_tc;
    end
    nMatch = size(npSub_tc_match, 2);

    % Build trial matrix and compute dF/F
    tc_mat = nan(nOn+nOff, nMatch, nTrials);
    for itrial = 1:nTrials
        on = stimOns(itrial);
        if ~isnan(on) && (on + nOn + nOff/2) <= size(npSub_tc_match, 1)
            tc_mat(:,:,itrial) = npSub_tc_match((on - nOff/2):(on - 1 + nOn + nOff/2), :);
        end
    end

    baselineFrames = (nOff/4+1):nOff/2;
    stimFrames     = (nOff/2+1):(nOff/2+nOn);

    tc_f    = mean(tc_mat(baselineFrames, :, :), 1);
    tc_dfof = (tc_mat - tc_f) ./ tc_f;
    clear tc_mat tc_f

    % Stimulus info
    Az  = celleqel2mat_padded(ret_inputStructure.tGratingAzimuthDeg);
    El  = celleqel2mat_padded(ret_inputStructure.tGratingElevationDeg);


    Azs = unique(Az);
    Els = unique(El);
    if min(Els) < 0, Els = fliplr(Els); end

    nStim = length(Azs) * length(Els);
    Stims = zeros(nStim, 2);
    idx = 1;
    for iEl = 1:length(Els)
        for iAz = 1:length(Azs)
            Stims(idx,:) = [Els(iEl) Azs(iAz)];
            idx = idx + 1;
        end
    end

    % Response per Az/El combination
    resp_dFoverF = squeeze(mean(tc_dfof(stimFrames, :, :), 1));  % nMatch x nTrials

    resp_by_stim = nan(length(Els), length(Azs), nMatch);
    for i_el = 1:length(Els)
        indE = find(El == Els(i_el));
        for i_az = 1:length(Azs)
            indA = find(Az == Azs(i_az));
            resp_by_stim(i_el, i_az, :) = nanmean(resp_dFoverF(:, intersect(indE, indA)), 2);
        end
    end

    % Build Ind_struct
    Ind_struct = [];
    for iStim = 1:nStim
        Ind_struct(iStim).all_trials = intersect(find(El == Stims(iStim,1)), find(Az == Stims(iStim,2)));
    end

    % Build grid
    [AzAz, ElEl] = meshgrid(Azs, Els);
    grid2.AzAz   = AzAz;
    grid2.ElEl   = ElEl;
    dAz = median(diff(Azs));
    dEl = median(diff(Els));
    [AzAz00, ElEl00] = meshgrid(Azs(1):(dAz/10):Azs(end), Els(1):(dEl/10):Els(end));
    grid2.AzAz00 = AzAz00;
    grid2.ElEl00 = ElEl00;

    % Fit retinotopy with shuffling
    Nshuf       = 100;
    Fit_struct  = [];
    PLOTIT_FIT  = 0;
    SAVEALLDATA = 0;

    fprintf('Day %d: fitting retinotopy (Nshuf=%d)...\n', id, Nshuf)
    for count_shuf = 0:Nshuf
        fprintf('  count_shuf: %d/%d\n', count_shuf, Nshuf)
        Im_mat_USE = zeros(nMatch, nStim);
        for iCond = 1:nStim
            ind_all = Ind_struct(iCond).all_trials;
            if count_shuf > 0
                ind_all = ind_all(randsample(length(ind_all), length(ind_all), 1));
            end
            Im_mat_USE(:, iCond) = mean(resp_dFoverF(:, ind_all), 2);
        end

        for iCell = 1:nMatch
            a = Im_mat_USE(iCell,:);
            if max(a,[],2) > 0
                b    = reshape(a', length(Azs), length(Els));
                data = b';
                PLOTIT_FIT  = 0;
                SAVEALLDATA = (count_shuf == 0);
                Fit_2Dellipse_ret_lbub
                if count_shuf == 0
                    Fit_struct(iCell).True.s_ = s;
                else
                    Fit_struct(iCell).Shuf(count_shuf).s_ = s;
                end
            end
        end
    end
    fprintf('Day %d: fitting done\n', id)

    % Goodness-of-fit assessment (lbub)
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
            if isfield(Fit_struct(iCell), 'Shuf') && ~isempty(Fit_struct(iCell).Shuf) && length(Fit_struct(iCell).Shuf) >= count_shuf && isfield(Fit_struct(iCell).Shuf(count_shuf), 's_') && ~isempty(Fit_struct(iCell).Shuf(count_shuf).s_)
                fit_shuf_vec(iCell,:,count_shuf) = [Fit_struct(iCell).Shuf(count_shuf).s_.x, ...
                    Fit_struct(iCell).Shuf(count_shuf).s_.Elhicut_50, Fit_struct(iCell).Shuf(count_shuf).s_.Azhicut_50, ...
                    Fit_struct(iCell).Shuf(count_shuf).s_.Elhicut_10, Fit_struct(iCell).Shuf(count_shuf).s_.Azhicut_10];
            end
        end
    end

    Npars       = size(fit_shuf_vec, 2);
    lbub_fits   = NaN(nMatch, Npars, 5);
    alpha_bound = 0.025;
    ind_shuf_lb = ceil(Nshuf * alpha_bound);
    ind_shuf_ub = ceil(Nshuf * (1 - alpha_bound));

    for iCell = 1:nMatch
        for count2 = 1:Npars
            i_sorted = sort(squeeze(fit_shuf_vec(iCell, count2, :)));
            lbub_fits(iCell, count2, 1) = i_sorted(ind_shuf_lb);
            lbub_fits(iCell, count2, 2) = i_sorted(ind_shuf_ub);
            lbub_fits(iCell, count2, 3) = mean(i_sorted);
            lbub_fits(iCell, count2, 5) = std(i_sorted);
        end
        lbub_fits(iCell,:,4) = fit_true_vec(iCell,:);
    end

    lbub_diff = lbub_fits(:,:,2) - lbub_fits(:,:,1);

    goodfit_ind  = [];
    azStep = ret_inputStructure.gratingAzimuthStepDeg;
    for iCell = 1:nMatch
        if lbub_diff(iCell,4) < azStep*2 && lbub_diff(iCell,5) < azStep*2
            goodfit_ind = [goodfit_ind iCell];
        end
    end

    % Remove RFs at retinotopy perimeter
    goodfit_ind2 = zeros(size(goodfit_ind));
    for i = 1:length(goodfit_ind)
        if sum(round(lbub_fits(goodfit_ind(i),4,4)) == [min(Azs) max(Azs)]) || ...
           sum(round(lbub_fits(goodfit_ind(i),5,4)) == [min(Els) max(Els)])
            continue
        end
        goodfit_ind2(i) = goodfit_ind(i);
    end
    goodfit_ind = goodfit_ind2(goodfit_ind2 ~= 0);
    fprintf('Day %d: %d good-fit cells (lbub)\n', id, length(goodfit_ind))

    fitAzimDeg = lbub_fits(:,4,4)';
    fitElevDeg = lbub_fits(:,5,4)';

    % RF distance from stimulus center (goodfit cells only)
    finalAzim = double(inputStructure(id).gratingAzimuthDeg);
    finalElev  = double(inputStructure(id).gratingElevationDeg);
    distAzim = nan(1, nMatch);
    distElev = nan(1, nMatch);
    distAzim(goodfit_ind) = fitAzimDeg(goodfit_ind);
    distElev(goodfit_ind) = fitElevDeg(goodfit_ind);
    ret_distance = sqrt((distAzim - finalAzim).^2 + (distElev - finalElev).^2);

    % Weighted average preferred position
    prefAzimDeg = nan(1, nMatch);
    prefElevDeg = nan(1, nMatch);
    for iCell = 1:nMatch
        w = max(resp_by_stim(:,:,iCell), 0);
        total_w = sum(w(:));
        if total_w > 0
            prefAzimDeg(iCell) = sum(w(:) .* AzAz(:)) / total_w;
            prefElevDeg(iCell) = sum(w(:) .* ElEl(:)) / total_w;
        end
    end

    % RF distance map via plotRFdistanceMap
    % For day 1, remap mask so pixel values index into matched cell space
    if id == 1
        mask_plot = zeros(size(referenceMasks));
        for c = 1:nMatch
            mask_plot(referenceMasks == match_ind(c)) = c;
        end
    else
        mask_plot = referenceMasks;
    end
    [distMap_matched{id}, dist_vec_matched{id}] = plotRFdistanceMap(...
        lbub_fits(:,:,4), goodfit_ind, mask_plot, finalAzim, finalElev);

    ret_npSub_tc_matched{id}   = npSub_tc_match;
    ret_distance_matched{id}   = ret_distance;
    resp_by_stim_matched{id}   = resp_by_stim;
    ret_dfof_trial_matched{id} = tc_dfof;
    lbub_fits_matched{id}      = lbub_fits;
    goodfit_ind_matched{id}    = goodfit_ind;
    r2_vec_matched{id}         = nan(1, nMatch);
    fitAzimDeg_matched{id}     = fitAzimDeg;
    fitElevDeg_matched{id}     = fitElevDeg;
    prefAzimDeg_matched{id}    = prefAzimDeg;
    prefElevDeg_matched{id}    = prefElevDeg;
    Azs_matched{id}            = Azs;
    Els_matched{id}            = Els;
end

save(fullfile(fnOut, [mouse '_ret_matched.mat']), ...
    'ret_npSub_tc_matched', 'ret_distance_matched', 'resp_by_stim_matched', 'ret_dfof_trial_matched', ...
    'trialIndSourceUsed', 'lbub_fits_matched', 'goodfit_ind_matched', 'r2_vec_matched', ...
    'fitAzimDeg_matched', 'fitElevDeg_matched', 'prefAzimDeg_matched', 'prefElevDeg_matched', ...
    'Azs_matched', 'Els_matched', 'distMap_matched', 'dist_vec_matched', '-v7.3');
fprintf('Saved to %s\n', fnOut)

end