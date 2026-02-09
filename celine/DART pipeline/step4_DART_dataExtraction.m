%% ===== EXPERIMENTAL SETUP AND DATA LOADING =====
clear all; clear global; clc
prompt = 'Enter name of instructions file: ';
instr = input(prompt, 's');
clear prompt
run(instr);

ds=instructions.ds;
run(ds); 

dataStructLabels = {'contrastxori'};
rc = behavConstsDART; % directories

% Input validation
if ~exist('expt', 'var')
    error('Dataset %s not found', ds);
end

day_id = str2double(instructions.session);


if length(expt) < day_id
    error('Day_id %d not valid for this dataset', day_id);
else
    match_day = expt(day_id).multiday_matchdays;

end

nd = 2; % hardcoding the number of days for now
mouse = expt(day_id).mouse;
experimentFolder = expt(day_id).exptType;
fnout = fullfile(rc.analysis, mouse);

% Set up file naming based on drug condition
if expt(day_id).multiday_timesincedrug_hours > 0
    dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end
fn_multi = fullfile(rc.analysis, experimentFolder, mouse, ['multiday_' dart_str]);

% Determine which session was used as reference for cell matching
x = instructions.refDay;
switch x
    case '1'
        pre = 1;  % baseline session, used as reference, is in the 1st position
        post = 2;
        fprintf('Baseline used as reference\n');
        allDays = [match_day,day_id];
        fprintf('Analyzing sessions: %s\n', num2str(allDays));
    case '2'
        pre = 2;
        post = 1;  % post-DART session, used as reference, is in the 1st position
        fprintf('Post-DART used as reference\n');
        allDays = [day_id, match_day];
        fprintf('Analyzing sessions: %s\n', num2str(allDays));
end
clear x instr

% Load the matched cell data
cd(fn_multi)
load(fullfile(fn_multi, 'timecourses.mat'))
load(fullfile(fn_multi, 'multiday_alignment.mat'))
load(fullfile(fn_multi, 'input.mat'))

% Rename input to inputStructure to avoid conflict with MATLAB's input function
inputStructure = input;
clear input

% Load pupil data - optional because some experiments don't have this
prompt = 'Include pupil data? Eye analysis must already be done for this experiment. y/n: ';
includePupil = input(prompt, 's');
if includePupil == 'y'
    pupil = cell(1, nd);
    for id = 1:nd
        % Get the appropriate folder path for each day
        mouse_temp = expt(allDays(id)).mouse;
        date = expt(allDays(id)).date;
        imgFolder = expt(allDays(id)).contrastxori_runs{1};
        dayPath = fullfile(rc.analysis, experimentFolder, mouse_temp, date, imgFolder);
        % Load pupil data from the correct day's folder
        pupil{id} = load(fullfile(dayPath, 'pupil.mat'));
    end
    clear mouse_temp date imgFolder dayPath
else
    fprintf('Not including pupil data\n');
end
clear prompt

%% ===== STIMULUS PARAMETER EXTRACTION =====
% Extract contrast, direction, orientation, and size for each trial
% Handle cases where actual trial count differs from planned trials

nOn = inputStructure(1).nScansOn;
nOff = inputStructure(1).nScansOff;

% Initialize trial condition arrays
tCon_match = cell(1, nd);
tDir_match = cell(1, nd);
tOri_match = cell(1, nd);
tSize_match = cell(1, nd);

% Extract stimulus parameters for each day
for id = 1:nd
    % Calculate actual number of trials based on recorded frames
    % (accounts for disruptions before full trial set completion)
    nTrials(id) = size(cellTCs_match{id}, 1) / (nOn + nOff);
    
    tCon_match{id} = celleqel2mat_padded(inputStructure(id).tGratingContrast(1:nTrials(id)));
    tDir_match{id} = celleqel2mat_padded(inputStructure(id).tGratingDirectionDeg(1:nTrials(id)));
    tOri_match{id} = tDir_match{id};
    % Convert directions to orientations (0-179d)
    tOri_match{id}(tDir_match{id} >= 180) = tDir_match{id}(tDir_match{id} >= 180) - 180;
    tSize_match{id} = celleqel2mat_padded(inputStructure(id).tGratingDiameterDeg(1:nTrials(id)));
end

% Get unique stimulus parameters
oris = unique(tOri_match{1});
dirs = unique(tDir_match{1});
cons = unique(tCon_match{1});
sizes = unique(tSize_match{1});
nOri = length(oris);
nCon = length(cons);
nDir = length(dirs);
nSize = length(sizes);

%% ===== TRIAL DEFINITION USING STIMULUS ONSET =====
% Extract cStimOn with variable method - automatically uses photodiode info if available,
% otherwise uses corrected input structure timing

stimOns =cell(1,nd); %this is the set of trial start times that will be used throughout the rest of this script
% correctedInputStructure = cell(1, nd);


for id = 1:nd
    mouse_temp = expt(allDays(id)).mouse;
    date = expt(allDays(id)).date;
    imgFolder = expt(allDays(id)).contrastxori_runs{1};
    imgMatFile = [imgFolder '_000_000.mat'];
    dataPath = fullfile(rc.data, mouse_temp, date, imgFolder);
    load(fullfile(dataPath, imgMatFile));
    if isfield(inputStructure(id),'tIdxSource')
        if ~isempty(inputStructure(id).tIdxSource)
            fprintf('Found previously corrected stim on timings for day %i.\n',id);
            switch inputStructure(id).tIdxSource
                case 'MW'
                    % stimOns{id} = cell2mat(inputStructure(id).stimOns_mwCounter);
                    stimOns{id} = inputStructure(id).stimOns_mwCounter;
                    disp('Using corrected mWorks counter.');
                case 'PD'
                    stimOns{id} = inputStructure(id).stimOns_photodiode;
                    disp('Using photodiode onsets.');
                case 'cS'
                    stimOns{id} = cell2mat(inputStructure(id).cStimOn);
                    disp('Using native cStimOn.');
            end
        else
            fprintf('No assigned stim on timing source for day %i,\n',id);
            sourceSel = input('Make selection here (PD/MW/cS): ','s');
            switch sourceSel
                case 'MW'
                    input_correct = counterValCorrect_noPhotodiode(inputStructure(id));
                    stimOns{id} = cell2mat(input_correct.cStimOn);
                    disp('Using newly calculated mWorks counter.');
                case 'PD'
                    if isfield(info,'frame')
                        [stimOns{id},~] = photoFrameFinder_Sanworks(info.frame);
                        disp('Using newly calculated photodiode onsets.');
                    else
                        error('No photodiode data for selected day.')
                    end
                case 'cS'
                    stimOns{id} = cell2mat(input_correct.cStimOn);
                    disp('Using native cStimOn.');
            end
        end
    else
        fprintf('Input struct was not previously corrected for day %i.\n',id);
        if isfield(info,'frame')
            [stimOns{id},~] = photoFrameFinder_Sanworks(info.frame);
            disp('Calculating stimOns from photodiode data.\n');
        else
            input_correct = counterValCorrect_noPhotodiode(inputStructure(id));
            stimOns{id} = cell2mat(input_correct.cStimOn);
            disp('Calculating stimOns from MWCounter.\n');
        end
    end
end

% save(fullfile(fn_multi,'correctedInputStructure.mat'),'correctedInputStructure')
% clear correctedInputStructure
clear mouse_temp date imgFolder imgMatFile dataPath info

%% drop trials
if instructions.tDropBool == true
    oldinput = inputStructure; % save input struct before modification
    clear inputStructure

    
    for id = 1:nd
        if id == 1
            dropTrials = instructions.tDropRefDay;
            fprintf('dropped %s\n trials from ref day', num2str(length(instructions.tDropRefDay)))
        elseif id == 2
            dropTrials = instructions.tDropMatchDay;
            fprintf('dropped %s\n trials from matched day', num2str(length(instructions.tDropMatchDay)))
        end
        newinput = trialDropper(oldinput(id),dropTrials,instructions.tDropAction);
        
        inputStructure(id) = newinput;
        clear newinput
        switch inputStructure(id).tIdxSource
            case 'PD'
                stimOns{id} = inputStructure(id).stimOns_photodiode;
            case 'MW'
                stimOns{id} = cell2mat(inputStructure(id).stimOns_mwCounter);
            case 'cS'
                stimOns{id} = cell2mat(inputStructure(id).cStimOn);
        end
    end
end
input=inputStructure;
save(fullfile(fn_multi,'input.mat'),'input')
clear input

%% Convert raw calcium timecourses to trial-structured dF/F data
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
    
    % Calculate activity statistics
    fractTimeActive_match{id} = zeros(1, nCells);
    
    % Calculate dF/F using baseline period
    data_f_match = mean(data_trial_match(1:(nOff/2), :, :), 1);
    data_dfof_trial_match{id} = bsxfun(@rdivide, bsxfun(@minus, data_trial_match, data_f_match), data_f_match);
    figure;plot(squeeze(mean(data_dfof_trial_match{id}(:,:,:),2, 'omitmissing')))
    % Calculate cell standard deviations for activity thresholding
    meansub_match = cellTCs_match{id} - nanmean(cellTCs_match{id}, 1);
    cellstd = nanstd(meansub_match, [], 1);
    cellstd_match{id} = cellstd;
    
    for iCell = 1:nCells
        fractTimeActive_match{id}(:, iCell) = length(find(meansub_match(:, iCell) > 3.*cellstd(1, iCell))) ./ nFrames;
    end
    
    clear data_trial_match data_f_match cellstd cStimOnTemp
end
clear meansub_match

% Define analysis windows
stimStart = nOff/2;
stimEnd = stimStart + nOn;
resp_win = (stimStart + 3):(stimEnd + 3); % at 15 Hz, 3 frames = ~200 ms
base_win = 1:(stimStart - 1);

%% ===== CELL RESPONSIVENESS ANALYSIS =====
% Identify cells that respond significantly to visual stimuli
[h_match, p_match, responsiveCellsMatch] = findResponsiveCells(data_dfof_trial_match, ...
    tCon_match, tDir_match, tSize_match, stimStart, stimEnd, nTrials, ...
    nCells, nDir, nCon, nSize, dirs, cons, sizes);

% Find all cells responsive on at least one day
resp_either_day = responsiveCellsMatch{1} | responsiveCellsMatch{2};
keep_cells_temp = find(resp_either_day);

% Count initial responsive cells by type
initial_red_resp = sum(red_ind_match(keep_cells_temp));
initial_green_resp = sum(~red_ind_match(keep_cells_temp));
initial_total_resp = length(keep_cells_temp);

% Total cells by type
total_red = sum(red_ind_match);
total_green = sum(~red_ind_match);
total_cells = length(red_ind_match);

% Choose outlier removal method
outlier_method = input('Choose outlier removal: (1) Manual, (2) STD threshold, (3) None: ');

outliers_all = [];
if outlier_method == 1
    % Manual outlier entry
    outliers_all = input('Enter outlier cell numbers (e.g., [1 5 10]): ');
    fprintf('Removing %d manually specified outlier cells\n', length(outliers_all));
    keep_cells = setdiff(keep_cells_temp, outliers_all);
    
elseif outlier_method == 2
    % STD threshold method
    std_threshold = input('Standard deviation threshold for outlier removal (e.g., 3): ');
    for id = 1:nd
        data_temp = data_dfof_trial_match{id};
        resp_max_temp = squeeze(max(nanmean(data_temp(resp_win, :, keep_cells_temp), 1), [], 2));
        thresh = nanmean(resp_max_temp) + std_threshold * std(resp_max_temp);
        outliers_day = keep_cells_temp(resp_max_temp > thresh);
        outliers_all = union(outliers_all, outliers_day);
    end
    fprintf('Removing %d outlier cells (>%.1f SD from mean)\n', length(outliers_all), std_threshold);
    keep_cells = setdiff(keep_cells_temp, outliers_all);
    
else
    % No outlier removal
    fprintf('Skipping outlier removal\n');
    keep_cells = keep_cells_temp;
end

% Count final keep cells by type
final_red_keep = sum(red_ind_match(keep_cells));
final_green_keep = sum(~red_ind_match(keep_cells));
nKeep = length(keep_cells);

% Within the keep cells, identify which ones are red vs green
red_cells_keep = logical(red_ind_match(keep_cells));
green_cells_keep = logical(~red_ind_match(keep_cells));

% Create and display summary table
cell_summary = table(...
    {'Total'; 'Red (HTP+)'; 'Green (HTP-)'}, ...
    [total_cells; total_red; total_green], ...
    [initial_total_resp; initial_red_resp; initial_green_resp], ...
    [length(outliers_all); sum(red_ind_match(outliers_all)); sum(~red_ind_match(outliers_all))], ...
    [nKeep; final_red_keep; final_green_keep], ...
    'VariableNames', {'CellType', 'Total_Cells', 'Responsive_Cells', 'Outliers_Removed', 'Final_Kept'});
fprintf('\n--- Cell Filtering Summary ---\n');
disp(cell_summary);

% Save cell analysis results
save('cell_analysis.mat', 'keep_cells', 'red_cells_keep', 'green_cells_keep', ...
    'h_match', 'p_match', 'responsiveCellsMatch', 'cell_summary');
save('cell_filtering_summary.mat', 'cell_summary');
writetable(cell_summary, 'cell_filtering_summary.csv');
fprintf('\nCell analysis saved to: cell_analysis.mat\n');
fprintf('Summary saved to: cell_filtering_summary.mat and cell_filtering_summary.csv\n');

clear initial_red_resp initial_green_resp initial_total_resp
clear final_red_keep final_green_keep total_red total_green total_cells
clear outliers_all keep_cells_temp resp_either_day data_temp resp_max_temp thresh outliers_day

%% ===== INITIAL RESPONSE MATRIX CALCULATION =====
% Create response matrices for all stimulus conditions using keep cells only

data_dfof_trial_keep = cell(1, nd);
data_resp_keep = cell(1, nd);
fullTC_keep = cell(1, nd);

for id = 1:nd
    % Subset data to keep cells only
    data_dfof_trial_keep{id} = data_dfof_trial_match{id}(:, :, keep_cells);
    fullTC_keep{id} = cellTCs_match{id}(:, keep_cells);
    
    % Create response matrix for all sizes/contrasts/directions
    data_resp = zeros(nKeep, nDir, nCon, nSize, 2);
    tCon = tCon_match{id}(:, 1:nTrials(id));
    tSize = tSize_match{id}(:, 1:nTrials(id));
    tDir = tDir_match{id}(:, 1:nTrials(id));
    
    for iDir = 1:nDir
        ind_dir = find(tDir == dirs(iDir));
        for iCon = 1:nCon
            ind_con = find(tCon == cons(iCon));
            for iSize = 1:nSize
                ind_size = find(tSize == sizes(iSize));
                ind_temp = intersect(ind_dir, ind_con);
                ind = intersect(ind_temp, ind_size);
                
                % Mean response and SEM
                data_resp(:, iDir, iCon, iSize, 1) = squeeze(nanmean(nanmean(data_dfof_trial_keep{id}(resp_win, ind, :), 1), 2));
                data_resp(:, iDir, iCon, iSize, 2) = squeeze(std(nanmean(data_dfof_trial_keep{id}(resp_win, ind, :), 1), [], 2) ./ sqrt(length(ind)));
            end
        end
    end
    data_resp_keep{id} = data_resp;
end

% Find preferred direction for each keep cell
prefDir_keep = cell(1, nd);
for id = 1:nd
    prefDir_keep{id} = findPreferredDirection(data_resp_keep{id});
end

clear tCon tSize tDir ind_dir ind_con ind_size ind_temp ind data_resp

%% ===== BEHAVIORAL STATE CLASSIFICATION =====
% Classify trials as running vs stationary based on wheel speed
% Classify trials by pupil size (arousal state) if available

fprintf('Analyzing wheel speed and behavioral states...\n');

% Import and process wheel data
wheel_speed = cell(1, nd);
for id = 1:nd
    wheel_speed{id} = wheelSpeedCalc(inputStructure(id), 32, expt(allDays(1)).wheelColor);
    fprintf('Day %d mean wheel speed: %.2f\n', id, nanmean(wheel_speed{id}));
end

% Remove wheel encoder "ticks" (small spurious movements)
wheel_speed_clean = cell(1, nd);
for id = 1:nd
    wheel_speed_clean{id} = wheel_speed{id};
    wheel_speed_clean{id}(abs(wheel_speed_clean{id}) < 4.9) = 0;
end

% Convert to trial-by-trial wheel timecourse and classify running trials
wheel_tc = cell(1, nd);
wheel_trial_avg = cell(1, nd);
wheel_tc_raw = cell(1, nd);
wheel_trial_avg_raw = cell(1, nd);
RIx = cell(1, nd); % Running index (logical)

for id = 1:nd
    cStimOnTemp = stimOns{id};
    wheel_tc{id} = nan(nOn + nOff, nTrials(id));
    wheel_tc_raw{id} = nan(nOn + nOff, nTrials(id));
    
    for iTrial = 1:nTrials(id)
        if ~isnan(cStimOnTemp(iTrial)) && (cStimOnTemp(iTrial) + nOn + nOff/2) <= nFrames && (cStimOnTemp(iTrial) - nOff/2) >= 1
            wheel_tc{id}(:, iTrial) = wheel_speed_clean{id}(cStimOnTemp(iTrial) - nOff/2:cStimOnTemp(iTrial) - 1 + nOn + nOff/2);
            wheel_tc_raw{id}(:, iTrial) = abs(wheel_speed{id}(cStimOnTemp(iTrial) - nOff/2:cStimOnTemp(iTrial) - 1 + nOn + nOff/2));
        end
    end
    
    % Calculate trial-averaged wheel speed
    wheel_trial_avg{id} = mean(wheel_tc{id}(nOff/2+1:nOn+nOff/2, :), 1, 'omitnan');
    wheel_trial_avg_raw{id} = mean(wheel_tc_raw{id}(nOff/2+1:nOn+nOff/2, :), 1, 'omitnan');
    
    % Classify running trials (threshold = 2 units above noise level)
    RIx{id} = wheel_trial_avg{id} > 2;
    
    fprintf('Day %d: %d/%d running trials (%.1f%%)\n', id, sum(RIx{id}), length(RIx{id}), 100*sum(RIx{id})/length(RIx{id}));
    clear cStimOnTemp
end

% Pupil analysis - classify large vs small pupil trials during stationary periods
pupilMeans = nan(nd, 3);
PIx_stat = cell(2, nd); % pupil index: {1} = large pupil stationary, {2} = small pupil stationary
motorByPupil = nan(nd, 2);

if includePupil == 'y'
    fprintf('Analyzing pupil size for arousal state classification...\n');
    
    % Combine pupil data from both days to set threshold
    statPupilBothDays = horzcat(pupil{pre}.rad.stim(~RIx{pre}), pupil{post}.rad.stim(~RIx{post}));
    statPupilThreshold = prctile(statPupilBothDays, 50);
    
    for id = 1:nd
        PIx_temp = pupil{id}.rad.stim > statPupilThreshold;
        PIx_stat{1, id} = logical(PIx_temp .* ~RIx{id}); % large pupil AND stationary
        PIx_stat{2, id} = logical(~PIx_temp .* ~RIx{id}); % small pupil AND stationary
        
        pupilMeans(id, 1) = mean(pupil{id}.rad.stim(PIx_stat{1, id}), 'omitmissing'); % large pupil stationary
        pupilMeans(id, 2) = mean(pupil{id}.rad.stim(PIx_stat{2, id}), 'omitmissing'); % small pupil stationary  
        pupilMeans(id, 3) = mean(pupil{id}.rad.stim(RIx{id}), 'omitmissing'); % running (any pupil)
        
        motorByPupil(id, 1) = mean(wheel_trial_avg_raw{id}(PIx_stat{1, id}), 'omitmissing');
        motorByPupil(id, 2) = mean(wheel_trial_avg_raw{id}(PIx_stat{2, id}), 'omitmissing');
        
        fprintf('Day %d pupil: Large=%.2f, Small=%.2f, Running=%.2f\n', ...
                id, pupilMeans(id,1), pupilMeans(id,2), pupilMeans(id,3));
    end
    clear statPupilBothDays statPupilThreshold PIx_temp
else
    fprintf('No pupil data - all stationary trials classified as small pupil\n');
    for id = 1:nd
        PIx_stat{2, id} = logical(~RIx{id}); % all stationary trials = small pupil
        PIx_stat{1, id} = false(size(RIx{id})); % no large pupil trials
    end
end

% Save behavioral state data
save('behavioral_state.mat', 'RIx', 'wheel_tc', 'wheel_trial_avg', 'PIx_stat', 'pupilMeans', 'motorByPupil');
clear wheel_speed wheel_speed_clean

%% ===== RESPONSE ANALYSIS =====
% Calculate responses at preferred direction for different behavioral states
% Test for surround suppression effects (large vs peak size responses)

fprintf('Calculating responses by behavioral state...\n');

% Initialize response arrays
pref_responses_stat = cell(1, nd);
pref_responses_loc = cell(1, nd); 
pref_responses_stat_largePupil = cell(1, nd);
pref_responses_stat_smallPupil = cell(1, nd);

% Individual trial response arrays
pref_allTrials_stat = cell(nCon, nSize, nd);
pref_allTrials_loc = cell(nCon, nSize, nd);
pref_allTrials_largePupil = cell(nCon, nSize, nd);
pref_allTrials_smallPupil = cell(nCon, nSize, nd);

% Significance testing arrays
h_keep = cell(1, nd);
h_largeVsPeak_keep = cell(1, nd);
p_largeVsPeak_keep = cell(1, nd);

% Time course arrays
tc_trial_avrg_stat = cell(1, nd);
tc_trial_avrg_loc = cell(1, nd);
tc_trial_avrg_stat_largePupil = cell(1, nd);
tc_trial_avrg_stat_smallPupil = cell(1, nd);

% Additional response matrices
stat_resp_keep = cell(1, nd);
conBySize_resp_stat_keep = cell(1, nd);
conBySize_resp_loc_keep = cell(1, nd);

% Trial count tracking
trialCounts = cell(2, nd);

for id = 1:nd
    fprintf('Processing day %d...\n', id);
    
    % Trial info for this day
    tCon = tCon_match{id}(:, 1:nTrials(id));
    tSize = tSize_match{id}(:, 1:nTrials(id));
    tDir = tDir_match{id}(:, 1:nTrials(id));
    
    % Behavioral state indices
    stat_inds = find(~RIx{id});
    loc_inds = find(RIx{id});
    ind_stat_largePupil = intersect(stat_inds, find(PIx_stat{1, id}));
    ind_stat_smallPupil = intersect(stat_inds, find(PIx_stat{2, id}));
    
    % Initialize response arrays for this day
    stat_resp = zeros(nKeep, nDir, nCon, nSize);
    temp_pref_responses_stat = zeros(nKeep, nCon, nSize);
    temp_pref_responses_loc = zeros(nKeep, nCon, nSize);
    temp_pref_responses_stat_largePupil = zeros(nKeep, nCon, nSize);
    temp_pref_responses_stat_smallPupil = zeros(nKeep, nCon, nSize);
    
    temp_tc_stat = nan((nOn + nOff), nKeep, nCon, nSize);
    temp_tc_loc = nan((nOn + nOff), nKeep, nCon, nSize);
    temp_tc_stat_largePupil = nan((nOn + nOff), nKeep, nCon, nSize);
    temp_tc_stat_smallPupil = nan((nOn + nOff), nKeep, nCon, nSize);
    
    % Significance testing array
    h = zeros(nKeep, nDir, nCon, nSize);
    trialCounts{1, id} = [];
    trialCounts{2, id} = [];
    
    % Calculate significance testing for all conditions
    for iDir = 1:nDir
        ind_dir = find(tDir == dirs(iDir));
        for iCon = 1:nCon
            ind_con = find(tCon == cons(iCon));
            for iSize = 1:nSize
                ind_size = find(tSize == sizes(iSize));
                ind_temp = intersect(ind_dir, ind_con);
                ind = intersect(ind_temp, ind_size);
                
                % Statistical test: response vs baseline (Bonferroni corrected)
                [h(:, iDir, iCon, iSize), ~] = ttest(nanmean(data_dfof_trial_keep{id}(resp_win, ind, :), 1), ...
                                                     nanmean(data_dfof_trial_keep{id}(base_win, ind, :), 1), ...
                                                     'dim', 2, 'tail', 'right', ...
                                                     'alpha', 0.05./(nDir*nCon*nSize-1));
            end
        end
    end
    
    % Calculate responses for all conditions (stationary trials only for stat_resp)
    for iDir = 1:nDir
        ind_dir = find(tDir == dirs(iDir));
        for iCon = 1:nCon
            ind_con = find(tCon == cons(iCon));
            for iSize = 1:nSize
                ind_size = find(tSize == sizes(iSize));
                ind_temp = intersect(ind_dir, ind_con);
                ind = intersect(ind_temp, ind_size);
                ind_stat = intersect(ind, stat_inds);
                
                % Calculate stationary responses for all conditions
                stat_resp(:, iDir, iCon, iSize) = squeeze(nanmean(nanmean(data_dfof_trial_keep{id}(resp_win, ind_stat, :), 1), 2));
            end
        end
    end
    
    % Calculate preferred direction responses using existing preferred directions
    for iCon = 1:nCon
        for iSize = 1:nSize
            % Initialize individual trial arrays for this contrast/size
            temp_all_stat = cell(1, nKeep);
            temp_all_loc = cell(1, nKeep);
            temp_all_largePupil = cell(1, nKeep);
            temp_all_smallPupil = cell(1, nKeep);
            
            for i = 1:nKeep
                % Get preferred direction for this cell
                temp_dir = dirs(prefDir_keep{id}(i));
                
                % Find trials matching preferred direction, contrast, and size
                dir_inds = find(tDir == temp_dir);
                con_inds = find(tCon == cons(iCon));
                size_inds = find(tSize == sizes(iSize));
                temp_trials = intersect(intersect(dir_inds, con_inds), size_inds);
                
                % Split by behavioral condition
                temp_trials_stat = intersect(temp_trials, stat_inds);
                temp_trials_loc = intersect(temp_trials, loc_inds);
                temp_trials_stat_largePupil = intersect(temp_trials, ind_stat_largePupil);
                temp_trials_stat_smallPupil = intersect(temp_trials, ind_stat_smallPupil);
                
                % Calculate time courses
                temp_tc_stat(:, i, iCon, iSize) = nanmean(data_dfof_trial_keep{id}(:, temp_trials_stat, i), 2);
                temp_tc_loc(:, i, iCon, iSize) = nanmean(data_dfof_trial_keep{id}(:, temp_trials_loc, i), 2);
                temp_tc_stat_largePupil(:, i, iCon, iSize) = nanmean(data_dfof_trial_keep{id}(:, temp_trials_stat_largePupil, i), 2);
                temp_tc_stat_smallPupil(:, i, iCon, iSize) = nanmean(data_dfof_trial_keep{id}(:, temp_trials_stat_smallPupil, i), 2);
                
                % Calculate individual trial responses (for pref_allTrials)
                temp_all_stat{i} = nanmean(data_dfof_trial_keep{id}(resp_win, temp_trials_stat, i), 1);
                temp_all_loc{i} = nanmean(data_dfof_trial_keep{id}(resp_win, temp_trials_loc, i), 1);
                temp_all_largePupil{i} = nanmean(data_dfof_trial_keep{id}(resp_win, temp_trials_stat_largePupil, i), 1);
                temp_all_smallPupil{i} = nanmean(data_dfof_trial_keep{id}(resp_win, temp_trials_stat_smallPupil, i), 1);
                
                % Calculate preferred responses (mean of response window)
                temp_pref_responses_stat(i, iCon, iSize) = nanmean(temp_all_stat{i}, 2);
                temp_pref_responses_loc(i, iCon, iSize) = nanmean(temp_all_loc{i}, 2);
                temp_pref_responses_stat_largePupil(i, iCon, iSize) = nanmean(temp_all_largePupil{i}, 2);
                temp_pref_responses_stat_smallPupil(i, iCon, iSize) = nanmean(temp_all_smallPupil{i}, 2);
                
                % Track trial counts
                trialCounts{1, id} = [trialCounts{1, id}, length(temp_trials_stat)];
                trialCounts{2, id} = [trialCounts{2, id}, length(temp_trials_loc)];
            end
            
            % Store individual trial responses for this contrast/size
            pref_allTrials_stat{iCon, iSize, id} = temp_all_stat;
            pref_allTrials_loc{iCon, iSize, id} = temp_all_loc;
            pref_allTrials_largePupil{iCon, iSize, id} = temp_all_largePupil;
            pref_allTrials_smallPupil{iCon, iSize, id} = temp_all_smallPupil;
        end
    end
    
    % Test for surround suppression: largest size response vs peak size response
    h_largeVsPeak = zeros(nKeep, 1);
    p_largeVsPeak = zeros(nKeep, 1);
    
    for iCell = 1:nKeep
        % Find peak size at preferred direction and highest contrast
        resp_by_size = squeeze(stat_resp(iCell, prefDir_keep{id}(iCell), end, :));
        [~, peakSize] = max(resp_by_size);
        
        % If peak size equals largest size, no suppression to test
        if peakSize == nSize
            h_largeVsPeak(iCell) = 0;
            p_largeVsPeak(iCell) = 1;
        else
            % Get trial indices for preferred direction and highest contrast
            ind_dir = find(tDir == dirs(prefDir_keep{id}(iCell)));
            ind_con = find(tCon == cons(end));
            ind_temp = intersect(ind_dir, ind_con);
            
            % Get trial indices for largest size and peak size (stationary only)
            ind_large = intersect(intersect(ind_temp, find(tSize == sizes(end))), stat_inds);
            ind_peak = intersect(intersect(ind_temp, find(tSize == sizes(peakSize))), stat_inds);
            
            % Perform statistical test if sufficient trials available
            if length(ind_large) >= 3 && length(ind_peak) >= 3
                resp_large = squeeze(nanmean(data_dfof_trial_keep{id}(resp_win, ind_large, iCell), 1));
                resp_peak = squeeze(nanmean(data_dfof_trial_keep{id}(resp_win, ind_peak, iCell), 1));
                [h_largeVsPeak(iCell), p_largeVsPeak(iCell)] = ttest2(resp_large, resp_peak, 'tail', 'left', 'alpha', 0.05);
            end
        end
    end
    
    % Store results for this day
    stat_resp_keep{id} = stat_resp;
    h_keep{id} = h;
    h_largeVsPeak_keep{id} = h_largeVsPeak;
    p_largeVsPeak_keep{id} = p_largeVsPeak;
    
    pref_responses_stat{id} = temp_pref_responses_stat;
    pref_responses_loc{id} = temp_pref_responses_loc;
    pref_responses_stat_largePupil{id} = temp_pref_responses_stat_largePupil;
    pref_responses_stat_smallPupil{id} = temp_pref_responses_stat_smallPupil;
    
    tc_trial_avrg_stat{id} = temp_tc_stat;
    tc_trial_avrg_loc{id} = temp_tc_loc;
    tc_trial_avrg_stat_largePupil{id} = temp_tc_stat_largePupil;
    tc_trial_avrg_stat_smallPupil{id} = temp_tc_stat_smallPupil;
end

% Summary of surround suppression test results
total_pass = sum(cellfun(@(x) sum(x == 1), h_largeVsPeak_keep));
total_cells = sum(cellfun(@length, h_largeVsPeak_keep));
fprintf('\nSurround suppression analysis: %d/%d cells pass large vs peak test (%.1f%%)\n', ...
        total_pass, total_cells, 100*total_pass/total_cells);

% Offer to plot surround suppression results
response = input('Plot surround suppression sanity check? (y/n): ', 's');
if strcmpi(response, 'y')
    sanityCheckPeakVsLarge(1, h_largeVsPeak_keep, prefDir_keep, stat_resp_keep, ...
                          data_dfof_trial_keep, tCon_match, tSize_match, tDir_match, ...
                          RIx, nTrials, dirs, cons, sizes, resp_win);
end

% Clean up temporary variables
clear temp_* ind_* dir_inds con_inds size_inds temp_trials* temp_dir
clear stat_inds loc_inds ind_stat_largePupil ind_stat_smallPupil
clear tCon tSize tDir stat_resp h h_largeVsPeak p_largeVsPeak

% Calculate normalized differences
[norm_diff, bsln_std] = calculateNormalizedDifference(pref_allTrials_stat, ...
    pref_allTrials_loc, pre, post, nCon, nKeep,nSize);
% in norm_diff, the first dimension is behavioral state. 1=stationary
% 2=running.
%% ===== NORMALIZED DIRECTION TUNING ANALYSIS =====
% Normalize direction tuning so preferred direction = 0d for each cell
% Uses stationary trials only to avoid locomotion effects on tuning

fprintf('Calculating normalized direction tuning...\n');

order = 1:nDir; % direction indices
norm_dir_resp = cell(1, nd);

for id = 1:nd
    normDir_temp = nan(nKeep, nDir, nCon, nSize);
    for iCell = 1:nKeep
        % Get preferred direction index for this cell
        pref_idx = prefDir_keep{id}(iCell);
        % Circularly shift directions so preferred = position 1
        newOrder = circshift(order, (pref_idx - 1) * -1);
        normDir_temp(iCell, :, :, :) = stat_resp_keep{id}(iCell, newOrder, :, :);
    end
    norm_dir_resp{id} = normDir_temp;
end
clear order normDir_temp pref_idx newOrder

%% ===== SAVE RESPONSE DATA =====
% Save all calculated response matrices and analysis results

save('tc_keep.mat', 'tc_trial_avrg_stat', 'tc_trial_avrg_loc', ...
     'tc_trial_avrg_stat_largePupil', 'tc_trial_avrg_stat_smallPupil', ...
     'fullTC_keep', 'data_dfof_trial_keep', 'prefDir_keep');

save('resp_keep.mat', 'stat_resp_keep', 'pref_responses_stat', 'pref_responses_loc', ...
     'pref_responses_stat_largePupil', 'pref_responses_stat_smallPupil', 'trialCounts', ...
     'conBySize_resp_stat_keep', 'conBySize_resp_loc_keep', 'h_keep', ...
     'h_largeVsPeak_keep', 'p_largeVsPeak_keep', 'norm_dir_resp', 'data_resp_keep','norm_diff','bsln_std');

fprintf('Response data saved to: tc_keep.mat and resp_keep.mat\n');

%% ===== CORRELATION ANALYSIS (INTERNEURON-PYRAMIDAL RELATIONSHIPS) =====
% Calculate noise correlations (trial-to-trial covariability) and 
% signal correlations (stimulus tuning similarity) between cell types

fprintf('Beginning correlation analysis...\n');

% Initialize data structures for correlation analysis
trialResp = cell(1, nd);        % Raw trial responses for each cell
subTrialResp = cell(1, nd);     % Mean-subtracted trial responses (for noise correlation)
conditionMeans = cell(1, nd);   % Mean response for each stimulus condition

% Calculate trial responses and condition means
for id = 1:nd
    fprintf('Calculating trial responses for day %d...\n', id);
    
    % Extract trial responses during stimulus period
    trialResp{id} = squeeze(mean(data_dfof_trial_keep{id}(stimStart:(stimStart + nOn - 1), :, :), 1, 'omitnan'));
    subTrialResp{id} = nan(size(trialResp{id}));
    conditionMeans{id} = nan(nDir, nCon, nSize, nKeep);
    
    % Get trial parameters
    tCon = tCon_match{id}(:, 1:nTrials(id));
    tDir = tDir_match{id}(:, 1:nTrials(id));
    tSize = tSize_match{id}(:, 1:nTrials(id));
    
    % Calculate condition means and subtract from trial responses (for noise correlation)
    for iDir = 1:nDir
        ind_dir = find(tDir == dirs(iDir));
        
        for iCon = 1:nCon
            ind_con = find(tCon == cons(iCon));
            ind_dir_con = intersect(ind_dir, ind_con);
            
            for iSize = 1:nSize
                ind_size = find(tSize == sizes(iSize));
                % Use stationary trials only to avoid locomotion artifacts
                ind_condition = intersect(intersect(ind_dir_con, ind_size), find(~RIx{id}));
                
                if ~isempty(ind_condition)
                    % Calculate condition mean for each cell
                    tempData = trialResp{id}(ind_condition, :);
                    cellMeans = mean(tempData, 1, 'omitnan');
                    
                    % Subtract condition mean from individual trials (removes stimulus-driven covariation)
                    subTrialResp{id}(ind_condition, :) = tempData - cellMeans;
                    conditionMeans{id}(iDir, iCon, iSize, :) = cellMeans;
                end
            end
        end
    end
end

% Setup correlation plotting preferences
response = input('Plot correlation figures as a sanity check? (y/n): ', 's');
doPlot = strcmpi(response, 'y') || strcmpi(response, 'yes');

if doPlot
    numToPlot = min(5, nKeep);
    cellsToPlot = randperm(nKeep, numToPlot);
    fprintf('Plotting enabled for %d randomly selected cells.\n', numToPlot);
else
    fprintf('Plotting disabled.\n');
    cellsToPlot = [];
end

% Calculate noise and signal correlations
noiseCorr = cell(1, nd);    % Noise correlations (trial-to-trial variability)
sigCorr = cell(1, nd);      % Signal correlations (stimulus tuning similarity)

for id = 1:nd
    fprintf('Calculating correlations for day %d: cell ', id);
    
    noiseCorr{id} = nan(2, nKeep);  % [correlation; p-value]
    sigCorr{id} = nan(2, nKeep);
    
    % Reshape condition means for signal correlation analysis
    condMeansReshaped = reshape(conditionMeans{id}, nCon * nDir * nSize, nKeep);
    
    for iCell = 1:nKeep
        if mod(iCell, 10) == 0, fprintf('%d ', iCell); end
        
        % Determine comparison cells based on cell type
        if ismember(iCell, red_cells_keep)
            % Red cell (interneuron): compare to all green cells (pyramidal)
            otherCells = find(green_cells_keep);
        else
            % Green cell (pyramidal): compare to other green cells (exclude self)
            otherCells = setdiff(find(green_cells_keep), iCell);
        end
        
        if ~isempty(otherCells)
            % Noise correlation: correlation of trial-to-trial fluctuations
            [R_noise, p_noise] = calculateCorrelation(subTrialResp{id}(find(~RIx{id}), iCell), ...
                                                      subTrialResp{id}(find(~RIx{id}), otherCells));
            noiseCorr{id}(:, iCell) = [R_noise; p_noise];
            
            % Optional plotting for noise correlation
            if doPlot && ismember(iCell, cellsToPlot)
                plotCorrelation(subTrialResp{id}(:, otherCells), subTrialResp{id}(:, iCell), ...
                               R_noise, iCell, 'NoiseCorr', id, pre);
            end
            
            % Signal correlation: correlation of stimulus tuning profiles
            [R_signal, p_signal] = calculateCorrelation(condMeansReshaped(:, iCell), ...
                                                        condMeansReshaped(:, otherCells));
            sigCorr{id}(:, iCell) = [R_signal; p_signal];
            

        end
    end
    fprintf('\n');
end

% Save correlation results
save(fullfile(fn_multi, 'HT_pyr_relationship.mat'), 'conditionMeans', 'sigCorr', 'noiseCorr', ...
     'trialResp', 'subTrialResp');

fprintf('Correlation analysis complete. Results saved to: HT_pyr_relationship.mat\n');

%% optional retinotopy alignment
response = input('Complete retinotopy alignement? (y/n): ', 's');
doRetino = strcmpi(response, 'y') || strcmpi(response, 'yes');

if doRetino
    [ret_npSub_tc_matched, ret_distance_matched,resp_by_stim_matched,ret_dfof_trial_matched,trialIndSourceUsed] = retinotopy_for_matched_data(nd, ...
        allDays, expt, mouse, fov_avg, masks, fitGeoTAf, ...
        instructions, inputStructure,match_ind,false);

   %make "keep" subsets for the matched retino data 
    ret_npSub_tc_keep = cell(1,nd);
    ret_distance_keep = cell(1,nd);
    resp_by_stim_keep=cell(1,nd);
    ret_dfof_trial_keep=cell(1,nd);

    for id = 1:nd
        ret_npSub_tc_keep{id}=ret_npSub_tc_matched{id}(:,keep_cells);
        ret_distance_keep{id}=ret_distance_matched{id}(:,keep_cells);
        resp_by_stim_keep{id}=resp_by_stim_matched{id}(:,:,keep_cells);
        ret_dfof_trial_keep{id}=ret_dfof_trial_matched{id}(:,keep_cells,:);
    end

save(fullfile(fn_multi, 'retino_aligned.mat'), 'ret_npSub_tc_matched', ...
    'ret_distance_matched','resp_by_stim_matched','ret_dfof_trial_matched',...
    'ret_npSub_tc_keep','ret_distance_keep','resp_by_stim_keep','ret_dfof_trial_keep','trialIndSourceUsed');

else
    fprint('Not including retinotopy aligment')
end

%%
% Clean up correlation analysis variables
clear condMeansReshaped tempData cellMeans numToPlot cellsToPlot
clear R_noise p_noise R_signal p_signal otherCells tCon tDir tSize
clear ind_dir ind_con ind_dir_con ind_size ind_condition

fprintf('\n=== ANALYSIS COMPLETE ===\n');
