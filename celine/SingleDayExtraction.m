% Single-day trial extraction and response analysis.
% Loads step1 outputs (neuropil-subtracted timecourses, masks, input struct),
% segments data into trials, identifies visually responsive cells, and computes
% trial-averaged timecourses at each cell's preferred direction for every
% contrast x size combination.

clear all; clear global; close all
clc

prompt = 'Enter name of instructions file: ';
instr = input(prompt, 's');
clear prompt
run(instr);

ds = instructions.ds;
run(ds);

rc = behavConstsDART;
dataStructLabels = {'contrastxori'};

day_id = str2double(instructions.session);
if day_id > length(expt)
    error('Day_id %d not valid for this dataset', day_id);
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

% Rename to avoid conflict with MATLAB's built-in input()
inputStructure = input;
clear input

%% Stimulus parameters

nOn  = inputStructure.nScansOn;   % number of imaging frames during stimulus
nOff = inputStructure.nScansOff;  % number of imaging frames during pre/post-stimulus period

% Estimate trial count from total recording length as a consistency check;
% will be overwritten below using the actual number of stimulus onsets
nTrials = size(npSub_tc, 1) / (nOn + nOff);

% Extract per-trial stimulus parameters
tCon  = celleqel2mat_padded(inputStructure.tGratingContrast(1:nTrials));
tDir  = celleqel2mat_padded(inputStructure.tGratingDirectionDeg(1:nTrials));
tOri  = tDir;
tOri(tDir >= 180) = tDir(tDir >= 180) - 180;  % collapse directions to orientations (0-179 deg)
tSize = celleqel2mat_padded(inputStructure.tGratingDiameterDeg(1:nTrials));

% Unique stimulus values and counts
cons  = unique(tCon);  nCon  = length(cons);
dirs  = unique(tDir);  nDir  = length(dirs);
oris  = unique(tOri);  nOri  = length(oris);
sizes = unique(tSize); nSize = length(sizes);

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
    % No timing source can be determined automatically; ask the user
    fprintf('Could not determine timing source automatically.\n');
    timingSource = input('Enter timing source to use (PD = photodiode, MW = mWorks counter, cS = native cStimOn): ', 's');

    % Save the choice into the input structure so this prompt is skipped next time
    inputStructure.stimTimingSource = timingSource;
    input = inputStructure; %#ok<NASGU>
    save(fullfile(fnout, 'input.mat'), 'input');
    clear input
    fprintf('Timing source saved to input.mat\n');
end

switch timingSource
    case 'MW'
        stimOns = inputStructure.stimOns_mwCounter;
    case 'PD'
        stimOns = inputStructure.stimOns_photodiode;
    case 'cS'
        stimOns = cell2mat(inputStructure.cStimOn);
    otherwise
        error('Unrecognised timing source: %s. Expected PD, MW, or cS.', timingSource);
end

nTrials           = length(stimOns);       % redefine using actual number of detected onsets
[nFrames, nCells] = size(npSub_tc);

%% Split into trials and compute dF/F

% Each trial window runs from nOff/2 frames before stimulus onset to
% nOff/2 frames after stimulus offset, giving a (nOn+nOff) x nTrials x nCells array.
% Trials that would fall outside the recording boundaries are left as NaN.
data_trial = nan(nOn + nOff, nTrials, nCells);

for iTrial = 1:nTrials
    if ~isnan(stimOns(iTrial)) && ...
       (stimOns(iTrial) - nOff/2) >= 1 && ...
       (stimOns(iTrial) - 1 + nOn + nOff/2) <= nFrames
        data_trial(:, iTrial, :) = npSub_tc(stimOns(iTrial) - nOff/2 : stimOns(iTrial) - 1 + nOn + nOff/2, :);
    end
end

% Baseline F = mean fluorescence over the first half of the pre-stimulus window
data_f          = mean(data_trial(1:(nOff/2), :, :), 1);
data_dfof_trial = bsxfun(@rdivide, bsxfun(@minus, data_trial, data_f), data_f);
clear data_trial data_f

%% Analysis windows

stimStart = nOff/2;           % frame index of stimulus onset within each trial
stimEnd   = stimStart + nOn;  % frame index of stimulus offset

% Response window is shifted 3 frames (~200 ms at 15 Hz) to account for
% calcium indicator rise time
resp_win = (stimStart + 3):(stimEnd + 3);
base_win = 1:(stimStart - 1);  % pre-stimulus baseline window

%% Responsive cells: significant response in at least one stimulus condition

% For each direction x contrast x size combination, test whether the mean
% response during the stimulus window is significantly greater than baseline.
% Uses a Bonferroni correction across all conditions. Requires >= 3 trials.
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
resp_cells = squeeze(any(any(any(h_resp, 2), 3), 4));  % nCells x 1 logical
keep_cells = find(resp_cells);
nKeep      = length(keep_cells);
red_cells  = logical(mask_label(keep_cells));  % red = interneuron (HTP+)
fprintf('%d/%d cells responsive (%d red, %d green)\n', nKeep, nCells, sum(red_cells), sum(~red_cells));

% Subset trial data and significance results to responsive cells only
data_dfof_trial_keep = data_dfof_trial(:, :, keep_cells);
h_keep               = h_resp(keep_cells, :, :, :);

%% Preferred direction per cell

% Compute mean response in the response window for each direction x contrast x size,
% then define preferred direction as the direction with the highest mean response
% averaged across all contrasts and sizes
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

resp_dir_avg     = squeeze(mean(mean(data_resp, 4), 3));  % nKeep x nDir, averaged across cons and sizes
[~, prefDir_idx] = max(resp_dir_avg, [], 2);              % index into dirs for each cell

%% Trial-averaged timecourses at preferred direction, for each contrast x size

% tc_trial_avrg:  (nOn+nOff) x nKeep x nCon x nSize
%   Full trial timecourse averaged over trials at each cell's preferred direction
% conBySize_resp: nKeep x nCon x nSize
%   Mean response amplitude in the response window at preferred direction

tc_trial_avrg  = nan(nOn + nOff, nKeep, nCon, nSize);
conBySize_resp = zeros(nKeep, nCon, nSize);

for iCell = 1:nKeep
    pref_dir_val = dirs(prefDir_idx(iCell));
    ind_dir      = find(tDir == pref_dir_val);
    for iCon = 1:nCon
        ind_con = find(tCon == cons(iCon));
        for iSize = 1:nSize
            ind_size = find(tSize == sizes(iSize));
            ind      = intersect(intersect(ind_dir, ind_con), ind_size);
            if ~isempty(ind)
                tc_trial_avrg(:, iCell, iCon, iSize) = nanmean(data_dfof_trial_keep(:, ind, iCell), 2);
                conBySize_resp(iCell, iCon, iSize)   = nanmean(nanmean(data_dfof_trial_keep(resp_win, ind, iCell), 1), 2);
            end
        end
    end
end

%% Save

save(fullfile(fnout, 'singleday_extraction.mat'), ...
    'tc_trial_avrg', 'conBySize_resp', 'data_dfof_trial_keep', ...
    'keep_cells', 'red_cells', 'prefDir_idx', ...
    'h_keep', 'cons', 'sizes', 'dirs', 'oris', ...
    'nOn', 'nOff', 'resp_win', 'base_win', 'tCon', 'tDir', 'tOri', 'tSize');

fprintf('Saved to %s\n', fullfile(fnout, 'singleday_extraction.mat'));
%% load in data for all the days
%concatenate it
%deal with the missing sizes by adding NaNs
%% find responsive cells




