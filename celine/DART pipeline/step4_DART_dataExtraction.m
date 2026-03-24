% Multi-day matched dataset trial extraction and response analysis.
% Loads step3 outputs (matched cell timecourses, alignment, input struct),
% segments data into trials, identifies visually responsive cells, and computes
% trial-averaged timecourses at each cell's preferred direction for every
% contrast x size combination.
clear all
if ~exist('instr', 'var')
    instr = input('Enter name of instructions file: ', 's');
end
run(instr);

ds     = instructions.ds;
run(ds);

rc = behavConstsDART;

day_id = str2double(instructions.session);

if day_id > length(expt)
    error('day_id %d not valid for this dataset', day_id);
end

match_day        = expt(day_id).multiday_matchdays;
nd               = 2;
mouse            = expt(day_id).mouse;
ExperimentFolder = expt(day_id).exptType;

if expt(day_id).multiday_timesincedrug_hours > 0
    dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end
fn_multi = fullfile(rc.analysis, ExperimentFolder, mouse, ['multiday_' dart_str]);

switch instructions.refDay
    case '1'
        pre = 1; post = 2;
        allDays = [match_day, day_id];
        fprintf('Baseline used as reference\n');
    case '2'
        pre = 2; post = 1;
        allDays = [day_id, match_day];
        fprintf('Post-DART used as reference\n');
end
fprintf('Analyzing sessions: %s\n', num2str(allDays));
clear instr

cd(fn_multi);
load(fullfile(fn_multi, 'timecourses.mat'))
load(fullfile(fn_multi, 'multiday_alignment.mat'))
load(fullfile(fn_multi, 'input.mat'))
inputStructure = input;
clear input

% Pupil data
if ~exist('includePupil', 'var')
    includePupil = input('Include pupil data? Eye analysis must already be done. (y/n): ', 's');
end
if includePupil == 'y'
    pupil = cell(1, nd);
    for id = 1:nd
        dayPath = fullfile(rc.analysis, ExperimentFolder, expt(allDays(id)).mouse, ...
            expt(allDays(id)).date, expt(allDays(id)).contrastxori_runs{1});
        pupil{id} = load(fullfile(dayPath, 'pupil.mat'));
    end
else
    fprintf('Not including pupil data\n');
end

% Stimulus parameters
nOn  = inputStructure(1).nScansOn;
nOff = inputStructure(1).nScansOff;

tCon_match  = cell(1, nd);
tDir_match  = cell(1, nd);
tOri_match  = cell(1, nd);
tSize_match = cell(1, nd);
nTrials     = zeros(1, nd);

for id = 1:nd
    nTrials(id) = size(cellTCs_match{id}, 1) / (nOn + nOff);
    tCon_match{id}  = celleqel2mat_padded(inputStructure(id).tGratingContrast(1:nTrials(id)));
    tDir_match{id}  = celleqel2mat_padded(inputStructure(id).tGratingDirectionDeg(1:nTrials(id)));
    tOri_match{id}  = tDir_match{id};
    tOri_match{id}(tDir_match{id} >= 180) = tDir_match{id}(tDir_match{id} >= 180) - 180;
    tSize_match{id} = celleqel2mat_padded(inputStructure(id).tGratingDiameterDeg(1:nTrials(id)));
end

cons  = unique(tCon_match{1});  nCon  = length(cons);
dirs  = unique(tDir_match{1});  nDir  = length(dirs);
oris  = unique(tOri_match{1});  nOri  = length(oris);
sizes = unique(tSize_match{1}); nSize = length(sizes);

% Stimulus onset times
stimOns = cell(1, nd);

for id = 1:nd
    if isfield(inputStructure(id), 'stimTimingSource') && ~isempty(inputStructure(id).stimTimingSource)
        timingSource = inputStructure(id).stimTimingSource;
        fprintf('Day %d: using previously saved timing source: %s\n', id, timingSource);
    elseif isfield(inputStructure(id), 'stimOns_photodiode') && ~isempty(inputStructure(id).stimOns_photodiode)
        timingSource = 'PD';
        fprintf('Day %d: no stimTimingSource field - using stimOns_photodiode\n', id);
    elseif isfield(inputStructure(id), 'stimOns_mwCounter') && ~isempty(inputStructure(id).stimOns_mwCounter)
        timingSource = 'MW';
        fprintf('Day %d: no stimTimingSource field - using stimOns_mwCounter\n', id);
    else
        fprintf('Day %d: could not determine timing source automatically.\n', id);
        timingSource = input(sprintf('Day %d - enter timing source (PD/MW/cS): ', id), 's');
        inputStructure(id).stimTimingSource = timingSource;
    end

    switch timingSource
        case 'MW'
            if isfield(inputStructure(id), 'stimOns_mwCounter') && ~isempty(inputStructure(id).stimOns_mwCounter)
                stimOns{id} = inputStructure(id).stimOns_mwCounter;
            else
                input_correct = counterValCorrect_noPhotodiode(inputStructure(id));
                stimOns{id} = cell2mat(input_correct.cStimOn);
            end
        case 'PD'
            if isfield(inputStructure(id), 'stimOns_photodiode') && ~isempty(inputStructure(id).stimOns_photodiode)
                stimOns{id} = inputStructure(id).stimOns_photodiode;
            else
                imgFolder  = expt(allDays(id)).contrastxori_runs{1};
                rawMatFile = fullfile(rc.data, expt(allDays(id)).mouse, expt(allDays(id)).date, ...
                    imgFolder, [imgFolder '_000_000.mat']);
                rawData = load(rawMatFile, 'info');
                [stimOns{id}, ~] = photoFrameFinder_Sanworks(rawData.info.frame);
            end
        case 'cS'
            stimOns{id} = cell2mat(inputStructure(id).cStimOn);
        otherwise
            error('Unrecognised timing source for day %d: %s. Expected PD, MW, or cS.', id, timingSource);
    end
    nTrials(id) = length(stimOns{id});
end

input = inputStructure; %#ok<NASGU>
save(fullfile(fn_multi, 'input.mat'), 'input');
clear input

% Trial dropping
if instructions.tDropBool
    oldinput = inputStructure;
    for id = 1:nd
        if id == 1
            dropTrials = instructions.tDropRefDay;
        else
            dropTrials = instructions.tDropMatchDay;
        end
        fprintf('Dropped %d trials from day %d\n', length(dropTrials), id);
        inputStructure(id) = trialDropper(oldinput(id), dropTrials, instructions.tDropAction);
        switch inputStructure(id).stimTimingSource
            case 'PD'; stimOns{id} = inputStructure(id).stimOns_photodiode;
            case 'MW'; stimOns{id} = cell2mat(inputStructure(id).stimOns_mwCounter);
            case 'cS'; stimOns{id} = cell2mat(inputStructure(id).cStimOn);
        end
        nTrials(id) = length(stimOns{id});
    end
    input = inputStructure; %#ok<NASGU>
    save(fullfile(fn_multi, 'input.mat'), 'input');
    clear input oldinput dropTrials
end

% Split into trials and compute dF/F
data_dfof_trial_match = cell(1, nd);

for id = 1:nd
    [nFrames, nCells] = size(cellTCs_match{id});
    data_trial = nan(nOn + nOff, nTrials(id), nCells);
    for iTrial = 1:nTrials(id)
        if ~isnan(stimOns{id}(iTrial)) && ...
           (stimOns{id}(iTrial) - nOff/2) >= 1 && ...
           (stimOns{id}(iTrial) - 1 + nOn + nOff/2) <= nFrames
            data_trial(:, iTrial, :) = cellTCs_match{id}(stimOns{id}(iTrial) - nOff/2 : stimOns{id}(iTrial) - 1 + nOn + nOff/2, :);
        end
    end
    data_f = mean(data_trial(1:(nOff/2), :, :), 1);
    data_dfof_trial_match{id} = bsxfun(@rdivide, bsxfun(@minus, data_trial, data_f), data_f);
end
clear data_trial data_f

% Analysis windows
stimStart = nOff/2;
stimEnd   = stimStart + nOn;
resp_win  = (stimStart + 1):(stimEnd + 1);
base_win  = 1:(stimStart - 1);

% Behavioral state: running vs stationary
wheel_tc          = cell(1, nd);
wheel_tc_raw      = cell(1, nd);
wheel_trial_avg   = cell(1, nd);
wheel_trial_avg_raw = cell(1, nd);
RIx = cell(1, nd);

for id = 1:nd
    [nFrames, ~] = size(cellTCs_match{id});
    ws = wheelSpeedCalc(inputStructure(id), 32, expt(allDays(1)).wheelColor);
    ws_clean = ws;
    ws_clean(abs(ws_clean) < 4.9) = 0;
    wheel_tc{id}     = nan(nOn + nOff, nTrials(id));
    wheel_tc_raw{id} = nan(nOn + nOff, nTrials(id));
    for iTrial = 1:nTrials(id)
        if ~isnan(stimOns{id}(iTrial)) && ...
           (stimOns{id}(iTrial) - nOff/2) >= 1 && ...
           (stimOns{id}(iTrial) - 1 + nOn + nOff/2) <= nFrames
            wheel_tc{id}(:, iTrial)     = ws_clean(stimOns{id}(iTrial) - nOff/2 : stimOns{id}(iTrial) - 1 + nOn + nOff/2);
            wheel_tc_raw{id}(:, iTrial) = abs(ws(stimOns{id}(iTrial) - nOff/2 : stimOns{id}(iTrial) - 1 + nOn + nOff/2));
        end
    end
    wheel_trial_avg{id}     = mean(wheel_tc{id}(nOff/2+1:nOn+nOff/2, :), 1, 'omitnan');
    wheel_trial_avg_raw{id} = mean(wheel_tc_raw{id}(nOff/2+1:nOn+nOff/2, :), 1, 'omitnan');
    RIx{id} = wheel_trial_avg{id} > 2;
    fprintf('Day %d: %d/%d running trials (%.1f%%)\n', id, sum(RIx{id}), nTrials(id), 100*sum(RIx{id})/nTrials(id));
end
clear ws ws_clean

% Behavioral state: pupil size (arousal)
pupilMeans   = nan(nd, 3);
motorByPupil = nan(nd, 2);
PIx_stat     = cell(2, nd);

if includePupil == 'y'
    statPupilBothDays = horzcat(pupil{pre}.rad.stim(~RIx{pre}), pupil{post}.rad.stim(~RIx{post}));
    statPupilThreshold = prctile(statPupilBothDays, 50);
    for id = 1:nd
        PIx_temp = pupil{id}.rad.stim > statPupilThreshold;
        PIx_stat{1, id} = logical(PIx_temp .* ~RIx{id});
        PIx_stat{2, id} = logical(~PIx_temp .* ~RIx{id});
        pupilMeans(id, 1) = mean(pupil{id}.rad.stim(PIx_stat{1, id}), 'omitmissing');
        pupilMeans(id, 2) = mean(pupil{id}.rad.stim(PIx_stat{2, id}), 'omitmissing');
        pupilMeans(id, 3) = mean(pupil{id}.rad.stim(RIx{id}),         'omitmissing');
        motorByPupil(id, 1) = mean(wheel_trial_avg_raw{id}(PIx_stat{1, id}), 'omitmissing');
        motorByPupil(id, 2) = mean(wheel_trial_avg_raw{id}(PIx_stat{2, id}), 'omitmissing');
        fprintf('Day %d pupil: Large=%.2f, Small=%.2f, Running=%.2f\n', ...
            id, pupilMeans(id,1), pupilMeans(id,2), pupilMeans(id,3));
    end
    clear statPupilBothDays statPupilThreshold PIx_temp
else
    fprintf('No pupil data - all stationary trials used without arousal split\n');
    for id = 1:nd
        PIx_stat{1, id} = false(1, nTrials(id));
        PIx_stat{2, id} = logical(~RIx{id});
    end
end

save(fullfile(fn_multi, 'behavioral_state.mat'), 'RIx', 'wheel_tc', 'wheel_trial_avg', 'PIx_stat', 'pupilMeans', 'motorByPupil');

% Responsive cells: significant response in at least one stimulus condition
h_match = cell(1, nd);
p_match = cell(1, nd);
responsiveCellsMatch = cell(1, nd);

for id = 1:nd
    [~, nCells] = size(cellTCs_match{id});
    h_id = zeros(nCells, nDir, nCon, nSize);
    p_id = ones(nCells, nDir, nCon, nSize);
    tCon = tCon_match{id}(:, 1:nTrials(id));
    tDir = tDir_match{id}(:, 1:nTrials(id));
    tSize = tSize_match{id}(:, 1:nTrials(id));
    for iDir = 1:nDir
        ind_dir = find(tDir == dirs(iDir));
        for iCon = 1:nCon
            ind_con = find(tCon == cons(iCon));
            for iSize = 1:nSize
                ind = intersect(intersect(ind_dir, ind_con), find(tSize == sizes(iSize)));
                if length(ind) >= 3
                    [h_id(:, iDir, iCon, iSize), p_id(:, iDir, iCon, iSize)] = ttest( ...
                        nanmean(data_dfof_trial_match{id}(resp_win, ind, :), 1), ...
                        nanmean(data_dfof_trial_match{id}(base_win, ind, :), 1), ...
                        'dim', 2, 'tail', 'right', ...
                        'alpha', 0.05 / (nDir * nCon * nSize - 1));
                end
            end
        end
    end
    h_match{id} = h_id;
    p_match{id} = p_id;
    responsiveCellsMatch{id} = squeeze(any(any(any(h_id, 2), 3), 4));
end
clear h_id p_id tCon tDir tSize ind_dir ind_con ind

resp_either_day = responsiveCellsMatch{1} | responsiveCellsMatch{2};
keep_cells_temp = find(resp_either_day);
fprintf('%d cells responsive on at least one day (%d red, %d green)\n', ...
    length(keep_cells_temp), sum(red_ind_match(keep_cells_temp)), sum(~red_ind_match(keep_cells_temp)));

% Outlier removal
outlier_method = 3;
%outlier_method = input('Choose outlier removal: (1) Manual, (2) STD threshold, (3) None: ');
outliers_all = [];
if outlier_method == 1
    outliers_all = input('Enter outlier cell numbers (e.g., [1 5 10]): ');
    keep_cells = setdiff(keep_cells_temp, outliers_all);
elseif outlier_method == 2
    std_threshold = input('Standard deviation threshold (e.g., 3): ');
    for id = 1:nd
        resp_max_temp = squeeze(max(nanmean(data_dfof_trial_match{id}(resp_win, :, keep_cells_temp), 1), [], 2));
        thresh = nanmean(resp_max_temp) + std_threshold * std(resp_max_temp);
        outliers_all = union(outliers_all, keep_cells_temp(resp_max_temp > thresh));
    end
    fprintf('Removing %d outlier cells (>%.1f SD from mean)\n', length(outliers_all), std_threshold);
    keep_cells = setdiff(keep_cells_temp, outliers_all);
else
    fprintf('Skipping outlier removal\n');
    keep_cells = keep_cells_temp;
end

nKeep            = length(keep_cells);
red_cells_keep   = logical(red_ind_match(keep_cells));
green_cells_keep = logical(~red_ind_match(keep_cells));

cell_summary = table(...
    {'Total'; 'Red (HTP+)'; 'Green (HTP-)'}, ...
    [length(red_ind_match); sum(red_ind_match); sum(~red_ind_match)], ...
    [length(keep_cells_temp); sum(red_ind_match(keep_cells_temp)); sum(~red_ind_match(keep_cells_temp))], ...
    [length(outliers_all); sum(red_ind_match(outliers_all)); sum(~red_ind_match(outliers_all))], ...
    [nKeep; sum(red_cells_keep); sum(green_cells_keep)], ...
    'VariableNames', {'CellType', 'Total_Cells', 'Responsive_Cells', 'Outliers_Removed', 'Final_Kept'});
disp(cell_summary);

save(fullfile(fn_multi, 'cell_analysis.mat'), 'keep_cells', 'red_cells_keep', 'green_cells_keep', ...
    'h_match', 'p_match', 'responsiveCellsMatch', 'cell_summary');

clear keep_cells_temp resp_either_day outliers_all resp_max_temp thresh

% Subset data to keep cells and compute response matrices for all conditions
data_dfof_trial_keep = cell(1, nd);
fullTC_keep          = cell(1, nd);
data_resp_keep       = cell(1, nd);

for id = 1:nd
    data_dfof_trial_keep{id} = data_dfof_trial_match{id}(:, :, keep_cells);
    fullTC_keep{id}          = cellTCs_match{id}(:, keep_cells);
    tCon  = tCon_match{id}(:, 1:nTrials(id));
    tDir  = tDir_match{id}(:, 1:nTrials(id));
    tSize = tSize_match{id}(:, 1:nTrials(id));
    data_resp = zeros(nKeep, nDir, nCon, nSize, 2);
    for iDir = 1:nDir
        ind_dir = find(tDir == dirs(iDir));
        for iCon = 1:nCon
            ind_con = find(tCon == cons(iCon));
            for iSize = 1:nSize
                ind = intersect(intersect(ind_dir, ind_con), find(tSize == sizes(iSize)));
                data_resp(:, iDir, iCon, iSize, 1) = squeeze(nanmean(nanmean(data_dfof_trial_keep{id}(resp_win, ind, :), 1), 2));
                data_resp(:, iDir, iCon, iSize, 2) = squeeze(std(nanmean(data_dfof_trial_keep{id}(resp_win, ind, :), 1), [], 2) ./ sqrt(length(ind)));
            end
        end
    end
    data_resp_keep{id} = data_resp;
end
clear data_resp tCon tDir tSize ind_dir ind_con ind

% Preferred direction per cell
prefDir_keep = cell(1, nd);
for id = 1:nd
    resp_dir_avg = squeeze(mean(mean(data_resp_keep{id}(:, :, :, :, 1), 4), 3));
    [~, prefDir_keep{id}] = max(resp_dir_avg, [], 2);
end
clear resp_dir_avg

% Response analysis: trial-averaged timecourses and mean responses,
% split by behavioral state, at each cell's preferred direction

pref_responses_stat            = cell(1, nd);
pref_responses_loc             = cell(1, nd);
pref_responses_stat_largePupil = cell(1, nd);
pref_responses_stat_smallPupil = cell(1, nd);
pref_allTrials_stat            = cell(nCon, nSize, nd);
pref_allTrials_loc             = cell(nCon, nSize, nd);
pref_allTrials_largePupil      = cell(nCon, nSize, nd);
pref_allTrials_smallPupil      = cell(nCon, nSize, nd);
h_keep                         = cell(1, nd);
h_largeVsPeak_keep             = cell(1, nd);
p_largeVsPeak_keep             = cell(1, nd);
tc_trial_avrg_stat             = cell(1, nd);
tc_trial_avrg_loc              = cell(1, nd);
tc_trial_avrg_stat_largePupil  = cell(1, nd);
tc_trial_avrg_stat_smallPupil  = cell(1, nd);
stat_resp_keep                 = cell(1, nd);
conBySize_resp_stat_keep       = cell(1, nd);
conBySize_resp_loc_keep        = cell(1, nd);
trialCounts = cell(2, nd);

for id = 1:nd
    tCon = tCon_match{id}(:, 1:nTrials(id));
    tDir = tDir_match{id}(:, 1:nTrials(id));
    tSize = tSize_match{id}(:, 1:nTrials(id));

    stat_inds           = find(~RIx{id});
    loc_inds            = find(RIx{id});
    ind_stat_largePupil = find(PIx_stat{1, id});
    ind_stat_smallPupil = find(PIx_stat{2, id});

    stat_resp     = zeros(nKeep, nDir, nCon, nSize);
    h             = zeros(nKeep, nDir, nCon, nSize);
    temp_pref_stat = zeros(nKeep, nCon, nSize);
    temp_pref_loc  = zeros(nKeep, nCon, nSize);
    temp_pref_lp   = zeros(nKeep, nCon, nSize);
    temp_pref_sp   = zeros(nKeep, nCon, nSize);
    temp_tc_stat   = nan(nOn + nOff, nKeep, nCon, nSize);
    temp_tc_loc    = nan(nOn + nOff, nKeep, nCon, nSize);
    temp_tc_lp     = nan(nOn + nOff, nKeep, nCon, nSize);
    temp_tc_sp     = nan(nOn + nOff, nKeep, nCon, nSize);
    trialCounts{1, id} = [];
    trialCounts{2, id} = [];

    % Significance test and stat_resp across all directions
    for iDir = 1:nDir
        ind_dir = find(tDir == dirs(iDir));
        for iCon = 1:nCon
            ind_con = find(tCon == cons(iCon));
            for iSize = 1:nSize
                ind      = intersect(intersect(ind_dir, ind_con), find(tSize == sizes(iSize)));
                ind_stat = intersect(ind, stat_inds);
                if length(ind) >= 3
                    [h(:, iDir, iCon, iSize), ~] = ttest( ...
                        nanmean(data_dfof_trial_keep{id}(resp_win, ind, :), 1), ...
                        nanmean(data_dfof_trial_keep{id}(base_win, ind, :), 1), ...
                        'dim', 2, 'tail', 'right', 'alpha', 0.05 / (nDir * nCon * nSize - 1));
                end
                stat_resp(:, iDir, iCon, iSize) = squeeze(nanmean(nanmean(data_dfof_trial_keep{id}(resp_win, ind_stat, :), 1), 2));
            end
        end
    end

    % Preferred direction timecourses and scalar responses
    for iCell = 1:nKeep
        dir_val  = dirs(prefDir_keep{id}(iCell));
        ind_dir  = find(tDir == dir_val);
        for iCon = 1:nCon
            ind_con = find(tCon == cons(iCon));
            for iSize = 1:nSize
                ind_pref = intersect(intersect(ind_dir, ind_con), find(tSize == sizes(iSize)));
                ind_s    = intersect(ind_pref, stat_inds);
                ind_l    = intersect(ind_pref, loc_inds);
                ind_lp   = intersect(ind_pref, ind_stat_largePupil);
                ind_sp   = intersect(ind_pref, ind_stat_smallPupil);

                temp_tc_stat(:, iCell, iCon, iSize) = nanmean(data_dfof_trial_keep{id}(:, ind_s,  iCell), 2);
                temp_tc_loc(:, iCell, iCon, iSize)  = nanmean(data_dfof_trial_keep{id}(:, ind_l,  iCell), 2);
                temp_tc_lp(:, iCell, iCon, iSize)   = nanmean(data_dfof_trial_keep{id}(:, ind_lp, iCell), 2);
                temp_tc_sp(:, iCell, iCon, iSize)   = nanmean(data_dfof_trial_keep{id}(:, ind_sp, iCell), 2);

                all_s  = nanmean(data_dfof_trial_keep{id}(resp_win, ind_s,  iCell), 1);
                all_l  = nanmean(data_dfof_trial_keep{id}(resp_win, ind_l,  iCell), 1);
                all_lp = nanmean(data_dfof_trial_keep{id}(resp_win, ind_lp, iCell), 1);
                all_sp = nanmean(data_dfof_trial_keep{id}(resp_win, ind_sp, iCell), 1);

                temp_pref_stat(iCell, iCon, iSize) = nanmean(all_s);
                temp_pref_loc(iCell, iCon, iSize)  = nanmean(all_l);
                temp_pref_lp(iCell, iCon, iSize)   = nanmean(all_lp);
                temp_pref_sp(iCell, iCon, iSize)   = nanmean(all_sp);

                pref_allTrials_stat{iCon, iSize, id}{iCell}       = all_s;
                pref_allTrials_loc{iCon, iSize, id}{iCell}        = all_l;
                pref_allTrials_largePupil{iCon, iSize, id}{iCell} = all_lp;
                pref_allTrials_smallPupil{iCon, iSize, id}{iCell} = all_sp;

                trialCounts{1, id} = [trialCounts{1, id}, length(ind_s)];
                trialCounts{2, id} = [trialCounts{2, id}, length(ind_l)];
            end
        end
    end

    % Surround suppression: largest vs peak size response
    h_largeVsPeak = zeros(nKeep, 1);
    p_largeVsPeak = zeros(nKeep, 1);
    for iCell = 1:nKeep
        resp_by_size_cell = squeeze(stat_resp(iCell, prefDir_keep{id}(iCell), end, :));
        [~, peakSize] = max(resp_by_size_cell);
        if peakSize == nSize
            p_largeVsPeak(iCell) = 1;
        else
            ind_base  = intersect(find(tDir == dirs(prefDir_keep{id}(iCell))), find(tCon == cons(end)));
            ind_large = intersect(intersect(ind_base, find(tSize == sizes(end))),    stat_inds);
            ind_peak  = intersect(intersect(ind_base, find(tSize == sizes(peakSize))), stat_inds);
            if length(ind_large) >= 3 && length(ind_peak) >= 3
                resp_large = squeeze(nanmean(data_dfof_trial_keep{id}(resp_win, ind_large, iCell), 1));
                resp_peak  = squeeze(nanmean(data_dfof_trial_keep{id}(resp_win, ind_peak,  iCell), 1));
                [h_largeVsPeak(iCell), p_largeVsPeak(iCell)] = ttest2(resp_large, resp_peak, 'tail', 'left', 'alpha', 0.05);
            end
        end
    end

    stat_resp_keep{id}                 = stat_resp;
    h_keep{id}                         = h;
    h_largeVsPeak_keep{id}             = h_largeVsPeak;
    p_largeVsPeak_keep{id}             = p_largeVsPeak;
    pref_responses_stat{id}            = temp_pref_stat;
    pref_responses_loc{id}             = temp_pref_loc;
    pref_responses_stat_largePupil{id} = temp_pref_lp;
    pref_responses_stat_smallPupil{id} = temp_pref_sp;
    tc_trial_avrg_stat{id}             = temp_tc_stat;
    tc_trial_avrg_loc{id}              = temp_tc_loc;
    tc_trial_avrg_stat_largePupil{id}  = temp_tc_lp;
    tc_trial_avrg_stat_smallPupil{id}  = temp_tc_sp;
    conBySize_resp_stat_keep{id}       = temp_pref_stat;
    conBySize_resp_loc_keep{id}        = temp_pref_loc;
end

fprintf('Surround suppression: %d/%d cells pass large vs peak test (%.1f%%)\n', ...
    sum(cellfun(@sum, h_largeVsPeak_keep)), sum(cellfun(@length, h_largeVsPeak_keep)), ...
    100 * sum(cellfun(@sum, h_largeVsPeak_keep)) / sum(cellfun(@length, h_largeVsPeak_keep)));

clear tCon tDir tSize stat_resp h h_largeVsPeak p_largeVsPeak resp_by_size_cell
clear stat_inds loc_inds ind_stat_largePupil ind_stat_smallPupil
clear temp_tc_stat temp_tc_loc temp_tc_lp temp_tc_sp
clear temp_pref_stat temp_pref_loc temp_pref_lp temp_pref_sp
clear all_s all_l all_lp all_sp ind_dir ind_con ind_pref ind_s ind_l ind_lp ind_sp
clear ind_base ind_large ind_peak ind ind_stat resp_large resp_peak

% Normalized direction tuning
order = 1:nDir;
norm_dir_resp = cell(1, nd);
for id = 1:nd
    normDir_temp = nan(nKeep, nDir, nCon, nSize);
    for iCell = 1:nKeep
        newOrder = circshift(order, (prefDir_keep{id}(iCell) - 1) * -1);
        normDir_temp(iCell, :, :, :) = stat_resp_keep{id}(iCell, newOrder, :, :);
    end
    norm_dir_resp{id} = normDir_temp;
end
clear order normDir_temp newOrder

% Normalized difference
[norm_diff, bsln_std] = calculateNormalizedDifference(pref_allTrials_stat, ...
    pref_allTrials_loc, pre, post, nCon, nKeep, nSize);

% Save
save(fullfile(fn_multi, 'tc_keep.mat'), 'tc_trial_avrg_stat', 'tc_trial_avrg_loc', ...
    'tc_trial_avrg_stat_largePupil', 'tc_trial_avrg_stat_smallPupil', ...
    'fullTC_keep', 'data_dfof_trial_keep', 'prefDir_keep');

save(fullfile(fn_multi, 'resp_keep.mat'), 'stat_resp_keep', 'pref_responses_stat', 'pref_responses_loc', ...
    'pref_responses_stat_largePupil', 'pref_responses_stat_smallPupil', 'trialCounts', ...
    'conBySize_resp_stat_keep', 'conBySize_resp_loc_keep', 'h_keep', ...
    'h_largeVsPeak_keep', 'p_largeVsPeak_keep', 'norm_dir_resp', 'data_resp_keep', 'norm_diff', 'bsln_std');

fprintf('Saved to %s\n', fn_multi);

% Correlation analysis
trialResp      = cell(1, nd);
subTrialResp   = cell(1, nd);
conditionMeans = cell(1, nd);

for id = 1:nd
    tCon  = tCon_match{id}(:, 1:nTrials(id));
    tDir  = tDir_match{id}(:, 1:nTrials(id));
    tSize = tSize_match{id}(:, 1:nTrials(id));
    trialResp{id}      = squeeze(mean(data_dfof_trial_keep{id}(stimStart:(stimStart + nOn - 1), :, :), 1, 'omitnan'));
    subTrialResp{id}   = nan(size(trialResp{id}));
    conditionMeans{id} = nan(nDir, nCon, nSize, nKeep);
    for iDir = 1:nDir
        ind_dir = find(tDir == dirs(iDir));
        for iCon = 1:nCon
            ind_con = find(tCon == cons(iCon));
            for iSize = 1:nSize
                ind_cond = intersect(intersect(intersect(ind_dir, ind_con), find(tSize == sizes(iSize))), find(~RIx{id}));
                if ~isempty(ind_cond)
                    tempData  = trialResp{id}(ind_cond, :);
                    cellMeans = mean(tempData, 1, 'omitnan');
                    subTrialResp{id}(ind_cond, :)          = tempData - cellMeans;
                    conditionMeans{id}(iDir, iCon, iSize, :) = cellMeans;
                end
            end
        end
    end
end
clear tCon tDir tSize ind_dir ind_con ind_cond tempData cellMeans

noiseCorr = cell(1, nd);
sigCorr   = cell(1, nd);
for id = 1:nd
    noiseCorr{id}  = nan(2, nKeep);
    sigCorr{id}    = nan(2, nKeep);
    condMeansReshaped = reshape(conditionMeans{id}, nCon * nDir * nSize, nKeep);
    stat_trials = find(~RIx{id});
    for iCell = 1:nKeep
        if red_cells_keep(iCell)
            otherCells = find(green_cells_keep);
        else
            otherCells = setdiff(find(green_cells_keep), iCell);
        end
        if ~isempty(otherCells)
            [R_noise, p_noise] = calculateCorrelation(subTrialResp{id}(stat_trials, iCell), ...
                                                      subTrialResp{id}(stat_trials, otherCells));
            noiseCorr{id}(:, iCell) = [R_noise; p_noise];
            [R_signal, p_signal]    = calculateCorrelation(condMeansReshaped(:, iCell), ...
                                                           condMeansReshaped(:, otherCells));
            sigCorr{id}(:, iCell)   = [R_signal; p_signal];
        end
    end
end
clear condMeansReshaped stat_trials R_noise p_noise R_signal p_signal otherCells

save(fullfile(fn_multi, 'HT_pyr_relationship.mat'), 'conditionMeans', 'sigCorr', 'noiseCorr', ...
    'trialResp', 'subTrialResp');
fprintf('Correlation analysis saved to HT_pyr_relationship.mat\n');

% Optional retinotopy alignment
if ~exist('doRetino', 'var')
    doRetino = strcmpi(input('Complete retinotopy alignment? (y/n): ', 's'), 'y');
end

if doRetino
    [ret_npSub_tc_matched, ret_distance_matched, resp_by_stim_matched, ret_dfof_trial_matched, trialIndSourceUsed] = ...
        retinotopy_for_matched_data(nd, allDays, expt, mouse, fov_avg, masks, fitGeoTAf, ...
            instructions, inputStructure, match_ind, false);

    ret_npSub_tc_keep   = cell(1, nd);
    ret_distance_keep   = cell(1, nd);
    resp_by_stim_keep   = cell(1, nd);
    ret_dfof_trial_keep = cell(1, nd);
    for id = 1:nd
        ret_npSub_tc_keep{id}   = ret_npSub_tc_matched{id}(:, keep_cells);
        ret_distance_keep{id}   = ret_distance_matched{id}(:, keep_cells);
        resp_by_stim_keep{id}   = resp_by_stim_matched{id}(:, :, keep_cells);
        ret_dfof_trial_keep{id} = ret_dfof_trial_matched{id}(:, keep_cells, :);
    end
    save(fullfile(fn_multi, 'retino_aligned.mat'), ...
        'ret_npSub_tc_matched', 'ret_distance_matched', 'resp_by_stim_matched', 'ret_dfof_trial_matched', ...
        'ret_npSub_tc_keep', 'ret_distance_keep', 'resp_by_stim_keep', 'ret_dfof_trial_keep', 'trialIndSourceUsed');
    fprintf('Retinotopy alignment saved to %s\n', fn_multi);
else
    fprintf('Skipping retinotopy alignment\n');
end

fprintf('\n=== ANALYSIS COMPLETE ===\n');