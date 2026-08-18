% Single-day trial extraction and response analysis.
% Loads step1 outputs (neuropil-subtracted timecourses, masks, input struct),
% segments data into trials, and computes population-level MAD and dF/F
% timecourses for HTP- and HTP+ cells separately across contrast x size x direction.

clear all
close all
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
    input = inputStructure;
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

%% Analysis windows

stimStart = nOff/2;
stimEnd   = stimStart + nOn;

resp_win = (stimStart + 1):(stimEnd + 1);
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

%% Population-level MAD and dF/F: HTP- and HTP+ cells, stationary trials, per con x size x dir
% MAD is baseline-normalized (percent change from baseline), matching Terlau et al. 2026.
% _tc  = full timecourse (nFrames x nCon x nSize x nDir)
% _resp = mean over response window (nCon x nSize x nDir)

stat_inds     = find(~RIx);
green_idx_all = find(~mask_label);
red_idx_all   = find(mask_label);

pop_MAD_tc_HTPminus   = nan(nOn + nOff, nCon, nSize, nDir);
pop_MAD_resp_HTPminus = nan(nCon, nSize, nDir);
pop_MAD_tc_HTPplus    = nan(nOn + nOff, nCon, nSize, nDir);
pop_MAD_resp_HTPplus  = nan(nCon, nSize, nDir);

pop_dfof_tc_HTPminus   = nan(nOn + nOff, nCon, nSize, nDir);
pop_dfof_resp_HTPminus = nan(nCon, nSize, nDir);
pop_dfof_tc_HTPplus    = nan(nOn + nOff, nCon, nSize, nDir);
pop_dfof_resp_HTPplus  = nan(nCon, nSize, nDir);

for iCon = 1:nCon
    ind_con = find(tCon == cons(iCon));
    for iSize = 1:nSize
        ind_size = find(tSize == sizes(iSize));
        for iDir = 1:nDir
            ind_dir = find(tDir == dirs(iDir));
            ind_s   = intersect(intersect(intersect(ind_con, ind_size), ind_dir), stat_inds);
            if ~isempty(ind_s)
                if ~isempty(green_idx_all)
                    pop_tc = mean(data_trial(:, ind_s, green_idx_all), 3);
                    dev    = abs(pop_tc - nanmedian(pop_tc, 2));
                    mad_tc = nanmedian(dev ./ pop_tc, 2);
                    mad_tc(nanmedian(pop_tc, 2) <= 0) = NaN;
                    baseline_mad = nanmean(mad_tc(base_win));
                    mad_tc = (mad_tc - baseline_mad) / baseline_mad;
                    pop_MAD_tc_HTPminus(:, iCon, iSize, iDir)   = mad_tc;
                    pop_MAD_resp_HTPminus(iCon, iSize, iDir)     = nanmean(mad_tc(resp_win));

                    dfof_tc = mean(mean(data_dfof_trial(:, ind_s, green_idx_all), 3), 2);
                    pop_dfof_tc_HTPminus(:, iCon, iSize, iDir)   = dfof_tc;
                    pop_dfof_resp_HTPminus(iCon, iSize, iDir)     = nanmean(dfof_tc(resp_win));
                end
                if ~isempty(red_idx_all)
                    pop_tc = mean(data_trial(:, ind_s, red_idx_all), 3);
                    dev    = abs(pop_tc - nanmedian(pop_tc, 2));
                    mad_tc = nanmedian(dev ./ pop_tc, 2);
                    mad_tc(nanmedian(pop_tc, 2) <= 0) = NaN;
                    baseline_mad = nanmean(mad_tc(base_win));
                    mad_tc = (mad_tc - baseline_mad) / baseline_mad;
                    pop_MAD_tc_HTPplus(:, iCon, iSize, iDir)     = mad_tc;
                    pop_MAD_resp_HTPplus(iCon, iSize, iDir)       = nanmean(mad_tc(resp_win));

                    dfof_tc = mean(mean(data_dfof_trial(:, ind_s, red_idx_all), 3), 2);
                    pop_dfof_tc_HTPplus(:, iCon, iSize, iDir)     = dfof_tc;
                    pop_dfof_resp_HTPplus(iCon, iSize, iDir)       = nanmean(dfof_tc(resp_win));
                end
            end
        end
    end
end

%% Figure: direction-averaged MAD and dF/F, nSize rows x 2 cols, contrast as darkness

t_frames   = (1:(nOn+nOff)) - nOff/2 - 1;
graylevels = linspace(0.75, 0, nCon);
ylbls      = {'MAD (frac. change)', '\DeltaF/F'};
type_vars  = {{pop_MAD_tc_HTPminus, pop_dfof_tc_HTPminus}, ...
              {pop_MAD_tc_HTPplus,  pop_dfof_tc_HTPplus}};
type_lbls  = {'HTP-', 'HTP+'};

for iType = 1:2
    metrics = type_vars{iType};
    figure('Position', [50 50 500 200*nSize]);
    axH = gobjects(nSize, 2);
    for iSize = 1:nSize
        for iCol = 1:2
            axH(iSize, iCol) = subplot(nSize, 2, (iSize-1)*2 + iCol);
            hold on
            for iCon = 1:nCon
                tc = nanmean(metrics{iCol}(:, iCon, iSize, :), 4);
                plot(t_frames, tc, 'Color', graylevels(iCon)*[1 1 1], 'LineWidth', 1.2, ...
                    'DisplayName', sprintf('%g%%', cons(iCon)*100));
            end
            yl = ylim;
            plot([0 nOn], [yl(1) yl(1)], 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
            set(gca, 'TickDir', 'out', 'Box', 'off');
            xlabel('Frame'); ylabel(ylbls{iCol});
            if iSize == 1, title(ylbls{iCol}); end
            if iCol == 1, text(-nOff/2, mean(yl), sprintf('%g deg', sizes(iSize)), ...
                'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'); end
        end
    end
    for iCol = 1:2
        yl = cell2mat(arrayfun(@(ax) ax.YLim, axH(:,iCol), 'UniformOutput', false));
        ylShared = [min(yl(:,1)) max(yl(:,2))];
        set(axH(:,iCol), 'YLim', ylShared);
        for iSize = 1:nSize
            plot(axH(iSize,iCol), [0 nOn], [ylShared(1) ylShared(1)], 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
        end
    end
    sgtitle(sprintf('Population MAD and dF/F – %s cells, stationary', type_lbls{iType}));
end

%% Save

save(fullfile(fnout, 'singleday_extraction_MAD.mat'), ...
    'pop_MAD_tc_HTPminus',   'pop_MAD_resp_HTPminus',   'pop_dfof_tc_HTPminus',  'pop_dfof_resp_HTPminus', ...
    'pop_MAD_tc_HTPplus',    'pop_MAD_resp_HTPplus',     'pop_dfof_tc_HTPplus',   'pop_dfof_resp_HTPplus', ...
    'cons', 'sizes', 'dirs', 'nOn', 'nOff', 'resp_win');

fprintf('Saved to %s\n', fullfile(fnout, 'singleday_extraction_MAD.mat'));