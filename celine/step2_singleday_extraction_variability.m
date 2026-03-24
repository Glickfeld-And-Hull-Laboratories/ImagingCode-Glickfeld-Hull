% Single-day trial extraction and response analysis.
% Loads step1 outputs (neuropil-subtracted timecourses, masks, input struct),
% segments data into trials, identifies visually responsive cells, and computes
% trial-averaged Median abs devation timecourses at each cell's preferred direction for every
% contrast x size combination.


% ds and day_id can be injected by the batch runner, or entered interactively.
if ~exist('ds', 'var')
    ds = input('Enter name of datasheet file: ', 's');
end
if ~exist('day_id', 'var')
    day_id = input('Enter session number: ');
end

run(ds);

rc = behavConstsDART;
dataStructLabels = {'contrastxori'};

if day_id > length(expt)
    error('day_id %d not valid for this dataset', day_id);
end

mouse            = expt(day_id).mouse;
expDate          = expt(day_id).date;
ExperimentFolder = expt(day_id).exptType;
runs             = eval(['expt(day_id).' cell2mat(dataStructLabels) '_runs']);
runFolder        = runs{1};

fnout = fullfile(rc.analysis, ExperimentFolder, mouse, expDate, runFolder);

% Load step1 outputs:
%   npSub_tc     - neuropil-subtracted timecourses (nFrames x nCells)
%   mask_cell    - labeled cell mask
%   mask_np      - neuropil masks
%   mask_label   - logical array indicating red (interneuron) cells
%   data_dfof    - dfof images from step1 segmentation
%   input        - MWorks behavioral/stimulus structure
load(fullfile(fnout, 'TCs.mat'))
load(fullfile(fnout, 'mask_cell.mat'))
load(fullfile(fnout, 'input.mat'))
load(fullfile(fnout, 'regOuts&Img.mat'))


% Rename to avoid conflict with MATLAB's built-in input()
inputStructure = input;
clear input

%% Stimulus parameters

nOn  = inputStructure.nScansOn;   % number of imaging frames during stimulus
nOff = inputStructure.nScansOff;  % number of imaging frames during pre/post-stimulus period

%% Stimulus onset times

% Determine timing source in order of preference:
%   1) input.stimTimingSource field (set on a previous run of this script)
%   2) whichever of stimOns_photodiode / stimOns_mwCounter is non-empty
%   3) prompt the user, then save the choice back to input.mat for future runs

if isfield(inputStructure, 'stimTimingSource') && ~isempty(inputStructure.stimTimingSource)
    timingSource = inputStructure.stimTimingSource;
    fprintf('Using previously saved timing source: %s\n', timingSource);

elseif isfield(inputStructure, 'stimOns_photodiode') && ~isempty(inputStructure.stimOns_photodiode)
    timingSource = 'PD';
    fprintf('No stimTimingSource field found - using stimOns_photodiode\n');

elseif isfield(inputStructure, 'stimOns_mwCounter') && ~isempty(inputStructure.stimOns_mwCounter)
    timingSource = 'MW';
    fprintf('No stimTimingSource field found - using stimOns_mwCounter\n');

else
    fprintf('Could not determine timing source automatically.\n');
    timingSource = input('Enter timing source to use (PD = photodiode, MW = mWorks counter, cS = native cStimOn): ', 's');
    inputStructure.stimTimingSource = timingSource;
    input = inputStructure; %#ok<NASGU>
    save(fullfile(fnout, 'input.mat'), 'input');
    clear input
    fprintf('Timing source saved to input.mat\n');
end

switch timingSource
    case 'MW'
        if isfield(inputStructure, 'stimOns_mwCounter') && ~isempty(inputStructure.stimOns_mwCounter)
            stimOns = inputStructure.stimOns_mwCounter;
        else
            fprintf('stimOns_mwCounter empty - calculating from counterValCorrect_noPhotodiode\n');
            input_correct = counterValCorrect_noPhotodiode(inputStructure);
            stimOns = cell2mat(input_correct.cStimOn);
        end
    case 'PD'
        if isfield(inputStructure, 'stimOns_photodiode') && ~isempty(inputStructure.stimOns_photodiode)
            stimOns = inputStructure.stimOns_photodiode;
        else
            fprintf('stimOns_photodiode empty - calculating from photoFrameFinder_Sanworks\n');
            rawMatFile = fullfile(rc.data, mouse, expDate, runFolder, [runFolder '_000_000.mat']);
            rawData = load(rawMatFile, 'info');
            if ~isfield(rawData.info, 'frame')
                error('No photodiode data found in %s', rawMatFile);
            end
            [stimOns, ~] = photoFrameFinder_Sanworks(rawData.info.frame);
        end
    case 'cS'
        stimOns = cell2mat(inputStructure.cStimOn);
    otherwise
        error('Unrecognised timing source: %s. Expected PD, MW, or cS.', timingSource);
end

nTrials           = length(stimOns);
[nFrames, nCells] = size(npSub_tc);

tCon  = celleqel2mat_padded(inputStructure.tGratingContrast(1:nTrials));
tDir  = celleqel2mat_padded(inputStructure.tGratingDirectionDeg(1:nTrials));
tOri  = tDir;
tOri(tDir >= 180) = tDir(tDir >= 180) - 180;
tSize = celleqel2mat_padded(inputStructure.tGratingDiameterDeg(1:nTrials));

cons  = unique(tCon);  nCon  = length(cons);
dirs  = unique(tDir);  nDir  = length(dirs);
oris  = unique(tOri);  nOri  = length(oris);
sizes = unique(tSize); nSize = length(sizes);

%% Split into trials and compute dF/F

data_trial = nan(nOn + nOff, nTrials, nCells);

for iTrial = 1:nTrials
    if ~isnan(stimOns(iTrial)) && ...
       (stimOns(iTrial) - nOff/2) >= 1 && ...
       (stimOns(iTrial) - 1 + nOn + nOff/2) <= nFrames
        data_trial(:, iTrial, :) = npSub_tc(stimOns(iTrial) - nOff/2 : stimOns(iTrial) - 1 + nOn + nOff/2, :);
    end
end

data_f          = mean(data_trial(1:(nOff/2), :, :), 1);
data_dfof_trial = bsxfun(@rdivide, bsxfun(@minus, data_trial, data_f), data_f);
clear data_trial data_f

%% Analysis windows

stimStart = nOff/2;
stimEnd   = stimStart + nOn;

resp_win = (stimStart + 1):(stimEnd + 1);%changed this from a 3-frame shifts to a 1-frame shifts on 3/3/26
base_win = 1:(stimStart - 1);

%% Behavioral state classification: running vs stationary

wheel_speed = wheelSpeedCalc(inputStructure, 32, expt(day_id).wheelColor);
wheel_speed_clean = wheel_speed;
wheel_speed_clean(abs(wheel_speed_clean) < 4.9) = 0;

wheel_tc = nan(nOn + nOff, nTrials);
for iTrial = 1:nTrials
    if ~isnan(stimOns(iTrial)) && ...
       (stimOns(iTrial) - nOff/2) >= 1 && ...
       (stimOns(iTrial) - 1 + nOn + nOff/2) <= nFrames
        wheel_tc(:, iTrial) = wheel_speed_clean(stimOns(iTrial) - nOff/2 : stimOns(iTrial) - 1 + nOn + nOff/2);
    end
end

wheel_trial_avg = mean(wheel_tc(nOff/2+1 : nOn+nOff/2, :), 1, 'omitnan');
RIx = wheel_trial_avg > 2;
fprintf('%d/%d running trials (%.1f%%)\n', sum(RIx), nTrials, 100*sum(RIx)/nTrials);

%% Behavioral state classification: pupil size (arousal)

if ~exist('includePupil', 'var')
    includePupil = input('Include pupil data? Eye analysis must already be done. (y/n): ', 's');
end

if includePupil == 'y'
    pupilFile = load(fullfile(fnout, 'pupil.mat'));
    pupilVars = fieldnames(pupilFile);

    pupil = [];
    for iVar = 1:length(pupilVars)
        candidate = pupilFile.(pupilVars{iVar});
        if isstruct(candidate) && isfield(candidate, 'rad') && isfield(candidate.rad, 'stim')
            pupil = candidate;
            break
        end
    end

    if isempty(pupil)
        fprintf('NOTE: pupil.mat found but does not contain expected .rad.stim structure - skipping pupil analysis\n');
        fprintf('      Variables in file: %s\n', strjoin(pupilVars, ', '));
        includePupil = 'n';
    end
end

if includePupil == 'y'
    statPupilThreshold = prctile(pupil.rad.stim(~RIx), 50);
    PIx_large = logical((pupil.rad.stim > statPupilThreshold) .* ~RIx);
    PIx_small = logical((pupil.rad.stim <= statPupilThreshold) .* ~RIx);
    fprintf('Pupil threshold: %.2f | Large stationary: %d trials | Small stationary: %d trials\n', ...
        statPupilThreshold, sum(PIx_large), sum(PIx_small));
else
    fprintf('No pupil data - all stationary trials used without arousal split\n');
    PIx_large = false(1, nTrials);
    PIx_small = ~RIx;
end

%% Responsive cells: significant response in at least one stimulus condition

h_resp = zeros(nCells, nDir, nCon, nSize);
for iDir = 1:nDir
    ind_dir = find(tDir == dirs(iDir));
    for iCon = 1:nCon
        ind_con = find(tCon == cons(iCon));
        for iSize = 1:nSize
            ind_size = find(tSize == sizes(iSize));
            ind      = intersect(intersect(ind_dir, ind_con), ind_size);
            if length(ind) >= 3
                [h_resp(:, iDir, iCon, iSize), ~] = ttest( ...
                    nanmean(data_dfof_trial(resp_win, ind, :), 1), ...
                    nanmean(data_dfof_trial(base_win, ind, :), 1), ...
                    'dim', 2, 'tail', 'right', ...
                    'alpha', 0.05 / (nDir * nCon * nSize - 1));
            end
        end
    end
end

% Keep any cell that passes significance in at least one condition
resp_cells = squeeze(any(any(any(h_resp, 2), 3), 4));
keep_single = find(resp_cells);
nKeep      = length(keep_single);
red_cells  = logical(mask_label(keep_single));
fprintf('%d/%d cells responsive (%d red, %d green)\n', nKeep, nCells, sum(red_cells), sum(~red_cells));

% Subset trial data and significance results to responsive cells only
data_dfof_trial_keep = data_dfof_trial(:, :, keep_single);
h_keep               = h_resp(keep_single, :, :, :);

% Per-size responsiveness mask: nKeep x nSize
% True if cell is significant in at least one direction x contrast at that size
resp_by_size = squeeze(any(any(h_keep, 2), 3));  % nKeep x nSize
fprintf('Cells responsive by size:\n');
for iSize = 1:nSize
    fprintf('  %g deg: %d cells\n', sizes(iSize), sum(resp_by_size(:, iSize)));
end

%% Preferred direction per cell


data_resp = zeros(nKeep, nDir, nCon, nSize);
for iDir = 1:nDir
    ind_dir = find(tDir == dirs(iDir));
    for iCon = 1:nCon
        ind_con = find(tCon == cons(iCon));
        for iSize = 1:nSize
            ind_size = find(tSize == sizes(iSize));
            ind      = intersect(intersect(ind_dir, ind_con), ind_size);
            data_resp(:, iDir, iCon, iSize) = squeeze(nanmean(nanmean(data_dfof_trial_keep(resp_win, ind, :), 1), 2));
        end
    end
end


resp_dir_avg     = squeeze(mean(mean(data_resp, 4), 3));
[~, prefDir_idx] = max(resp_dir_avg, [], 2);

%% Trial-averaged timecourses at preferred direction, split by behavioral state

stat_inds = find(~RIx);
loc_inds  = find(RIx);
ind_large = find(PIx_large);
ind_small = find(PIx_small);

tc_trial_MAD_stat       = nan(nOn + nOff, nKeep, nCon, nSize);
tc_trial_MAD_loc        = nan(nOn + nOff, nKeep, nCon, nSize);
tc_trial_MAD_largePupil = nan(nOn + nOff, nKeep, nCon, nSize);
tc_trial_MAD_smallPupil = nan(nOn + nOff, nKeep, nCon, nSize);

conBySize_MAD_stat       = zeros(nKeep, nCon, nSize);
conBySize_MAD_loc        = zeros(nKeep, nCon, nSize);
conBySize_MAD_largePupil = zeros(nKeep, nCon, nSize);
conBySize_MAD_smallPupil = zeros(nKeep, nCon, nSize);

for iCell = 1:nKeep
    pref_dir_val = dirs(prefDir_idx(iCell));
    ind_dir      = find(tDir == pref_dir_val);
    for iCon = 1:nCon
        ind_con = find(tCon == cons(iCon));
        for iSize = 1:nSize
            ind_size  = find(tSize == sizes(iSize));
            ind_stim  = intersect(intersect(ind_dir, ind_con), ind_size);

            ind_s  = intersect(ind_stim, stat_inds);
            ind_l  = intersect(ind_stim, loc_inds);
            ind_lp = intersect(ind_stim, ind_large);
            ind_sp = intersect(ind_stim, ind_small);

            if ~isempty(ind_s)
                tc_trial_MAD_stat(:, iCell, iCon, iSize)       = nanmean(data_dfof_trial_keep(:, ind_s, iCell), 2);
                conBySize_MAD_stat(iCell, iCon, iSize)          = nanmean(nanmean(data_dfof_trial_keep(resp_win, ind_s, iCell), 1), 2);
            end
            if ~isempty(ind_l)
                tc_trial_MAD_loc(:, iCell, iCon, iSize)         = nanmean(data_dfof_trial_keep(:, ind_l, iCell), 2);
                conBySize_MAD_loc(iCell, iCon, iSize)           = nanmean(nanmean(data_dfof_trial_keep(resp_win, ind_l, iCell), 1), 2);
            end
            if ~isempty(ind_lp)
                tc_trial_MAD_largePupil(:, iCell, iCon, iSize)  = nanmean(data_dfof_trial_keep(:, ind_lp, iCell), 2);
                conBySize_MAD_largePupil(iCell, iCon, iSize)     = nanmean(nanmean(data_dfof_trial_keep(resp_win, ind_lp, iCell), 1), 2);
            end
            if ~isempty(ind_sp)
                tc_trial_MAD_smallPupil(:, iCell, iCon, iSize)  = nanmean(data_dfof_trial_keep(:, ind_sp, iCell), 2);
                conBySize_MAD_smallPupil(iCell, iCon, iSize)     = nanmean(nanmean(data_dfof_trial_keep(resp_win, ind_sp, iCell), 1), 2);
            end
        end
    end
end
figure;plot(mean(tc_trial_MAD_stat(:,:,3,3),2));title('grand mean across cells, con 3 size 3, keepSingle')
%% Save

save(fullfile(fnout, 'singleday_extraction.mat'), ...
    'tc_trial_MAD_stat', 'tc_trial_MAD_loc', ...
    'tc_trial_MAD_largePupil', 'tc_trial_MAD_smallPupil', ...
    'conBySize_MAD_stat', 'conBySize_MAD_loc', ...
    'conBySize_MAD_largePupil', 'conBySize_MAD_smallPupil', ...
    'data_dfof_trial_keep', 'keep_single', 'red_cells', 'prefDir_idx', ...
    'RIx', 'PIx_large', 'PIx_small', 'wheel_tc', 'wheel_trial_avg', ...
    'h_keep', 'resp_by_size', ...
    'cons', 'sizes', 'dirs', 'oris', ...
    'nOn', 'nOff', 'resp_win', 'base_win', 'tCon', 'tDir', 'tOri', 'tSize');

fprintf('Saved to %s\n', fullfile(fnout, 'singleday_extraction.mat'));
%% Optional retinotopy alignment

if isfield(expt(day_id), 'ret_run') && ~isempty(expt(day_id).ret_run)
    if ~exist('doRetino', 'var')
        response = input('Complete retinotopy alignment? (y/n): ', 's');
        doRetino = strcmpi(response, 'y') || strcmpi(response, 'yes');
    end
else
    doRetino = false;
    fprintf('No ret_run defined for this session - skipping retinotopy alignment\n');
end

if doRetino
    referenceFOV = data_avg;

    if ~exist('validation_choice', 'var')
        validation_choice = strcmpi(input('Plot validation traces for random cells? (y/n): ', 's'), 'y');
    end

    [ret_npSub_tc, ret_distance, resp_by_stim, ret_dfof_trial, trialIndSourceUsed] = ...
        retinotopy_for_singleday(day_id, expt, mouse, inputStructure, ...
            referenceFOV, mask_cell, timingSource, rc, validation_choice);

    ret_npSub_tc_keep   = ret_npSub_tc(:, keep_single);
    ret_distance_keep   = ret_distance(keep_single);
    resp_by_stim_keep   = resp_by_stim(:, :, keep_single);
    ret_dfof_trial_keep = ret_dfof_trial(:, keep_single, :);

    save(fullfile(fnout, 'retino_aligned.mat'), ...
        'ret_npSub_tc', 'ret_distance', 'resp_by_stim', 'ret_dfof_trial', ...
        'ret_npSub_tc_keep', 'ret_distance_keep', 'resp_by_stim_keep', 'ret_dfof_trial_keep', ...
        'trialIndSourceUsed');

    fprintf('Retinotopy alignment saved to %s\n', fullfile(fnout, 'retino_aligned.mat'));
end