%% set paths and read in data
clear all; clear global; clc

prompt = 'Enter ds (e.g., DART_V1_YM90K_Celine): ';
ds = input(prompt, 's');
clear x prompt

dataStructLabels = {'contrastxori'};
rc =  behavConstsDART; %directories
eval(ds);
day_id = input('Enter reference day id ');
pre_day = expt(day_id).multiday_matchdays;
allDays = [day_id,pre_day];

nd=2; %hardcoding the number of days for now
mouse = expt(day_id).mouse;
experimentFolder = expt(day_id).exptType;
fnout = fullfile(rc.achAnalysis,mouse);

if expt(day_id).multiday_timesincedrug_hours>0
    dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end
fn_multi = fullfile(rc.achAnalysis,experimentFolder,mouse,['multiday_' dart_str]);

%check which day is pre and which is post, set variables accordingly
prompt = 'Which sesson was used as reference for matching: 0- baseline, 1- post-DART';
            x = input(prompt);
            switch x
                case 0
                    pre=1; %baeline session, used as reference, is in the 1st position
                    post=2;
                    "baseline used as reference"
                case 1
                  pre=2;
                  post=1; %post-DART session, used as reference, is in the 1st position  
                  "post-DART used as reference"
            end
clear x prompt

% load the data
cd(fn_multi)
load(fullfile(fn_multi,'timecourses.mat'))
load(fullfile(fn_multi,'multiday_alignment.mat'))
load(fullfile(fn_multi,'input.mat'))
% Rename input to inputStructure to avoid conflict with MATLAB's input function
inputStructure = input;
clear input

% load pupil data - optional because some experiments don't have this
prompt = 'Include pupil data? Eye analysis must be done for this experiment. 1 = yes, 0 = no'
    includePupil=input(prompt);
    if includePupil
            
            pupil=cell(1,nd);
            for id=1:nd
                % Get the appropriate folder path for each day
                mouse = expt(allDays(id)).mouse;
                date = expt(allDays(id)).date;
                imgFolder = expt(allDays(id)).contrastxori_runs{1};
                dayPath = fullfile(rc.achAnalysis,experimentFolder, mouse, date, imgFolder);
                % Load pupil data from the correct day's folder
                pupil{id}=load(fullfile(dayPath,'pupil.mat'));
            end
    else
            fprint('Not including pupil data')

    end


%% BEGIN DATA EXTRACTION
%% Get stimulus parameters
nOn = inputStructure(1).nScansOn;
nOff = inputStructure(1).nScansOff;


%tells the contrast, direction and orientation for each trial each day
tCon_match = cell(1,nd);
tDir_match = cell(1,nd);
tOri_match = cell(1,nd);
tSize_match = cell(1,nd);

%find the contrasts, directions and orientations for each day
%in case of instances where the number of trails actually collected was not
%consistent with the number mWorks thinks occured, only take 1:nTrials as
%dictated above based on the number of frames recorded
for id = 1:nd
    nTrials(id) = size(cellTCs_match{id},1)/(nOn+nOff); %to account for times 
%when there is a disruption before the full set of trials is collcted, I'm 
%determining the number of trials each day by how many frames of data I 
%have divided by the number of frames per trial
    tCon_match{id} = celleqel2mat_padded(inputStructure(id).tGratingContrast(1:nTrials(id)));
    tDir_match{id} = celleqel2mat_padded(inputStructure(id).tGratingDirectionDeg(1:nTrials(id)));
    tOri_match{id} = tDir_match{id};
    tOri_match{id}(find(tDir_match{id}>=180)) = tDir_match{id}(find(tDir_match{id}>=180))-180;
    tSize_match{id} = celleqel2mat_padded(inputStructure(id).tGratingDiameterDeg(1:nTrials(id)));
end
oris = unique(tOri_match{1});
dirs = unique(tDir_match{1});
cons = unique(tCon_match{1});
sizes = unique(tSize_match{1});
nOri = length(oris);
nCon = length(cons);
nDir = length(dirs);
nSize = length(sizes);
%% Extract cStimOn with variable method
% this will automatically use the photodiode info.frames if that field
% exists, and use the corrected input structure if info.frames does not
% exist. With either method, a cStimOn field for each day will
% be generated and added to the input strucutre, then used for trial definition




for id = 1:nd %cycle through days
    mouse = expt(allDays(id)).mouse;
    date = expt(allDays(id)).date;
    imgFolder = expt(allDays(id)).contrastxori_runs{1};
    imgMatFile = [imgFolder '_000_000.mat'];
    dataPath = fullfile(rc.achData, mouse, date, imgFolder);
    load(fullfile(dataPath,imgMatFile));
    if isfield(inputStructure,'cStimOn')
        fprintf('cStimOn already exists in the input structure\n');
    elseif isfield(info, 'frame')
        fprintf('No cStim on, but field "frame" exists in info structure\n and will be used to find cStimOn');
        [cStimOnTemp stimOffsTemp] = photoFrameFinder_Sanworks(info.frame);
        inputStructure(id).cStimOn = cStimOnTemp;
        inputStructure(id).stimOffs=stimOffsTemp;
        clear  cStimOnTemp stimOffsTemp
    else 
        fprintf('No cStim on and field "frame" does not exist in info structure,\n will run counterValCorrect_noPhotodiode');
        correctedInputStructure=NaN(1,nd);
        correctedInputStructure(:,id)=counterValCorrect_noPhotodiode(inputStructure(:,id));
        inputStructure=correctedInputStructure; %replace the input structure with the new version  
    end
end
    
%% Use cStimOn to define trials
data_dfof_trial_match = cell(1,nd); %make an empty array that is 1 by however many days there are (1X2 usually)
fractTimeActive_match = cell(1,nd);
cellstd_match = cell(1,nd);

for id = 1:nd
    cStimOnTemp=cell2mat(inputStructure(id).cStimOn);
    nTrials(id) = length(cStimOnTemp);
    [nFrames nCells] = size(cellTCs_match{id});
    
    data_trial_match = nan(nOn+nOff,nTrials(id),nCells);
    
    for iTrial = 1:nTrials(id)
      if ~isnan(cStimOnTemp(iTrial)) & (cStimOnTemp(iTrial)+nOn+nOff/2)<nFrames
        data_trial_match(:,iTrial,:) = cellTCs_match{id}(cStimOnTemp(iTrial)-nOff/2:cStimOnTemp(iTrial)-1+nOn+nOff/2,:);
        %cellTCs_match is the raw F timecourse for each cell for the entire
        %session
      end
    end
    
    
    fractTimeActive_match{id} = zeros(1,nCells);
    
    
    data_f_match = mean(data_trial_match(1: (nOff/2),:,:),1);
    data_dfof_trial_match{id} = bsxfun(@rdivide,bsxfun(@minus,data_trial_match,data_f_match),data_f_match);
    meansub_match = cellTCs_match{id}-nanmean(cellTCs_match{id},1);
    cellstd = nanstd(meansub_match,[],1);
    cellstd_match{id}=cellstd;
    for iCell = 1:nCells
        fractTimeActive_match{id}(:,iCell) = length(find(meansub_match(:,iCell)>3.*cellstd(1,iCell)))./nFrames; %not currently saving this or doing anything with it - do we want to?
    end
    clear data_trial_match data_f_match cellstd
end

stimStart = nOff/2;
stimEnd = stimStart+nOn;

resp_win = (stimStart+3):(stimEnd+3); %at 15 hz, 3 frames = ~200 ms.
base_win = 1: stimStart-1;


%% Narrow down cells and make response matrix
% Find cells that were responsive on either day ("keep") and identify red
%  cells within that subset

% This is a function that finds the stimulus conditions for which each cell was
% significantly responsive and makes a list of which cells are responsive
[h_match, p_match, responsiveCellsMatch] = findResponsiveCells(data_dfof_trial_match, ...
    tCon_match, tDir_match, tSize_match, stimStart, stimEnd, nTrials, ...
    nCells, nDir, nCon, nSize, dirs, cons, sizes);

% User parameters
remove_outliers = input('Remove outliers? (1=yes, 0=no): ');
if remove_outliers
    std_threshold = input('Standard deviation threshold for outlier removal (e.g., 3): ');
end

% Find all cells responsive on at least one day (regardless of red/green)
resp_either_day = responsiveCellsMatch{1} | responsiveCellsMatch{2};
keep_cells_temp = find(resp_either_day);

% Count initial responsive cells by type
initial_red_resp = sum(red_ind_match(keep_cells_temp));
initial_green_resp = sum(~red_ind_match(keep_cells_temp));
initial_total_resp = length(keep_cells_temp);

% Total cells by type (for reference)
total_red = sum(red_ind_match);
total_green = sum(~red_ind_match);
total_cells = length(red_ind_match);

% Optional outlier removal
if remove_outliers
    outliers_all = [];
    for id = 1:nd
        data_temp = responsiveCellsMatch{id}(keep_cells_temp,:,:,:,1);
        resp_max_temp = max(squeeze(max(data_temp(:,:,:,1),[],3)),[],2);
        thresh = nanmean(resp_max_temp) + std_threshold*std(resp_max_temp);
        outliers_all = union(outliers_all, keep_cells_temp(resp_max_temp > thresh));
    end
    
    fprintf('Removing %d outlier cells (>%.1f SD from mean)\n', length(outliers_all), std_threshold);
    keep_cells = setdiff(keep_cells_temp, outliers_all);
else
    fprintf('Skipping outlier removal\n');
    keep_cells = keep_cells_temp;
    outliers_all = [];
end

% Count final keep cells by type
final_red_keep = sum(red_ind_match(keep_cells));
final_green_keep = sum(~red_ind_match(keep_cells));
nKeep = length(keep_cells);

% Within the keep cells, identify which ones are red
red_cells_keep = red_ind_match(keep_cells);
green_cells_keep=~red_ind_match(keep_cells);

% Create summary table
cell_summary = table(...
    {'Total'; 'Red (HTP+)'; 'Green (HTP-)'}, ...
    [total_cells; total_red; total_green], ...
    [initial_total_resp; initial_red_resp; initial_green_resp], ...
    [length(outliers_all); sum(red_ind_match(outliers_all)); sum(~red_ind_match(outliers_all))], ...
    [nKeep; final_red_keep; final_green_keep], ...
    'VariableNames', {'CellType', 'Total_Cells', 'Responsive_Cells', 'Outliers_Removed','Final_Kept'});

% Display summary
fprintf('\n--- Cell Filtering Summary ---\n');
disp(cell_summary);

% Save summary table
save('cell_filtering_summary.mat', 'cell_summary');
writetable(cell_summary, 'cell_filtering_summary.csv');

fprintf('\nSummary saved to: cell_filtering_summary.mat and cell_filtering_summary.csv\n');

% Clean up temporary variables
clear initial_red_resp initial_green_resp initial_total_resp
clear final_red_keep final_green_keep 
clear total_red total_green total_cells outliers_all
%% Reduce the data to only the keep cells and make response matrices
data_dfof_trial_keep=cell(1, nd);
data_resp_keep=cell(1,nd);
fullTC_keep=cell(1,nd);
for id = 1:nd
    %make a set of data for the keep cells only and use this moving forward
    data_dfof_trial_keep{id}=data_dfof_trial_match{id}(:,:,keep_cells);
    fullTC_keep{id}=cellTCs_match{id}(:,keep_cells);
    % restrucutre this as a response matrix for all
    % sizes/contrasts/orientations
    data_resp = zeros(nKeep, nDir, nCon,nSize,2);
    tCon = tCon_match{id}(:,1:nTrials(id));
    tSize = tSize_match{id}(:,1:nTrials(id));  
    tDir = tDir_match{id}(:,1:nTrials(id));
        for iDir = 1:nDir
            ind_dir = find(tDir == dirs(iDir));
            for iCon = 1:nCon
                ind_con = find(tCon == cons(iCon));
                for iSize = 1:nSize                
                    ind_size = find(tSize == sizes(iSize));
                    ind_temp = intersect(ind_dir,ind_con); %for every orientation and then every contrast, find trials with that con/dir/size combination
                    ind=intersect(ind_temp,ind_size);
                    data_resp(:,iDir,iCon,iSize,1) = squeeze(nanmean(nanmean(data_dfof_trial_keep{id}(resp_win,ind,:),1),2));
                    data_resp(:,iDir,iCon,iSize,2) = squeeze(std(nanmean(data_dfof_trial_keep{id}(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
                end
            end
        end
data_resp_keep{id} = data_resp;
end

%% For keep cells only, find preferred direction
%for the keep cells subset, find the preferred direction for each cell,
%then average over all trials at that direction for all sizes and contrasts
%to find each cell's trial-averaged timecourse ad each cell's
%trial-averaged mean response (ie, the mean over the respons window, so
%dropping the time dimension) for each size and contrast.
for id  = 1:nd
    prefDir_keep{id}=findPreferredDirection(data_resp_keep{id});
end
%% Get wheel speed, running vs. stationary trials
wheel_speed = cell(1,nd);
%import the wheel data
for id = 1:nd
    wheel_speed{id} = wheelSpeedCalc(inputStructure(id),32,expt(allDays(1)).wheelColor); 
    nanmean(wheel_speed{id})
end

%make a version where the wheel "ticks" are removed
wheel_speed_clean = cell(1,nd);
for id = 1:nd
    wheel_speed_clean{id}=wheel_speed{id};
    wheel_speed_clean{id}(abs(wheel_speed_clean{id})<4.9)=0;
end

%convert to trial-by-trial wheel timecourse and mean wheel speed per trial
%and find a set of ind for running trials (RIx, logical)
wheel_tc = cell(1,nd);
wheel_trial_avg= cell(1,nd);
wheel_tc_raw = cell(1,nd);
wheel_trial_avg_raw= cell(1,nd);
RIx = cell(1,nd);


for id = 1:nd
    cStimOnTemp=cell2mat(inputStructure(id).cStimOn);
    wheel_tc{id}=nan(nOn+nOff, nTrials(id));
    wheel_tc_raw{id}=nan(nOn+nOff, nTrials(id));
    for iTrial = 1:nTrials(id)
        if ~isnan(cStimOnTemp(iTrial)) & (cStimOnTemp(iTrial)+nOn+nOff/2)<nFrames
            wheel_tc{id}(:,iTrial) = wheel_speed_clean{id}(cStimOnTemp(iTrial)-nOff/2:cStimOnTemp(iTrial)-1+nOn+nOff/2);
            wheel_tc_raw{id}(:,iTrial) = abs(wheel_speed{id}(cStimOnTemp(iTrial)-nOff/2:cStimOnTemp(iTrial)-1+nOn+nOff/2));
        end
    end
    wheel_trial_avg{id} = mean(wheel_tc{id}(nOff/2:nOn+nOff/2,:),1);
    wheel_trial_avg_raw{id} = mean(wheel_tc_raw{id}(nOff/2:nOn+nOff/2,:),1);
    RIx{id} = wheel_trial_avg{id}>2; %~5 is the noise level in the wheel movement, but this should get averaged out over a trial
    clear cStimOnTemp
end



%% Find large vs small pupil trials - 
% %this will only run if the user chose Y to including pupil data above
%Make matrices of NaNs that will get populated if there is pupil data or
%serve as dummies if there is not
pupilMeans = nan(nd,3);
PIx_stat = cell(2,nd); %pupil index for each day, first cell is inds for 
% stationary large pupil, second cell is inds for stationary small pupil
motorByPupil = nan(nd,2);

if includePupil=='Y'
    statPupilBothDays =horzcat(pupil{pre}.rad.stim(~RIx{pre}),pupil{post}.rad.stim(~RIx{post})); %combine all the pupil values for the two days
    statPupilThreshold=prctile(statPupilBothDays,50);
    %plot(statPupilBothDays); hline(statPupilThreshold); xlabel('trials, both days');ylabel('pupil diam'); hold off
    %print(fullfile(fn_multi,'pupilTraceWThreshold.pdf'),'-dpdf');
    
    pupilMeans = nan(nd,3);
    % for each day, the first column is the mean pupil size for stat trials
    % below threshold, the second column is the mean pupil size for stat trials
    % above threshold, and the third column is the mean pupil size for all
    % running trials
    PIx_stat = cell(2,nd); %pupil index for each day, first cell is inds for 
    % stationary large pupil, second cell is inds for stationary small pupil
    motorByPupil = nan(nd,2);
    for id = 1:nd
        PIx_temp=pupil{id}.rad.stim > statPupilThreshold;
        PIx_stat{1,id}= logical(PIx_temp.*~RIx{id});
        PIx_stat{2,id}= logical(~PIx_temp.*~RIx{id});
        pupilMeans(id,1)=mean(pupil{id}.rad.stim(PIx_stat{1,id}), 'omitmissing'); %passes pupil threshold but isn't running
        pupilMeans(id,2)=mean(pupil{id}.rad.stim(PIx_stat{2,id}), 'omitmissing'); %doesn't pass pupil threshold AND isn't running
        pupilMeans(id,3)=mean(pupil{id}.rad.stim(RIx{id}), 'omitmissing'); %is running, regardless of pupil size
        motorByPupil(id,1)=mean(wheel_trial_avg_raw{id}(PIx_stat{1,id}),'omitmissing');
        motorByPupil(id,2)=mean(wheel_trial_avg_raw{id}(PIx_stat{2,id}),'omitmissing');
    end
else 
    for id = 1:nd
        PIx_stat{2,id}= logical(~RIx{id}); %calling all stationary trials small pupil if there is no real pupil data
    end
end
    save(fullfile(fn_multi,'pupilMeans.mat'),'pupilMeans','motorByPupil','PIx_stat');

%% Create timecourse matrices and response matrices (averaged over time) at preferred direction, based on behavioral state

% Initialize cell arrays for keep data
data_dfof_trial_keep = cell(1, nd);
data_resp_keep = cell(1, nd);
fullTC_keep = cell(1, nd);

% Preferred responses
pref_responses_stat = cell(1, nd);
pref_responses_loc = cell(1, nd);
pref_responses_stat_largePupil = cell(1, nd);
pref_responses_stat_smallPupil = cell(1, nd);

% ConBySize matrices (responses at pref direction for each contrast/size)
conBySize_resp_stat_keep = cell(1, nd);
conBySize_resp_loc_keep = cell(1, nd);

% Significance testing
h_keep = cell(1, nd);
h_largeVsPeak_keep = cell(1, nd);
p_largeVsPeak_keep = cell(1, nd);

% Time courses
tc_trial_avrg_stat = cell(1, nd);
tc_trial_avrg_loc = cell(1, nd);
tc_trial_avrg_stat_largePupil = cell(1, nd);
tc_trial_avrg_stat_smallPupil = cell(1, nd);

% Peak responses
pref_peak_stat = cell(1, nd);
pref_peak_loc = cell(1, nd);

% All trials data for statistics
pref_allTrials_stat = cell(nCon, nSize, nd);
pref_allTrials_loc = cell(nCon, nSize, nd);
pref_allTrials_largePupil = cell(nCon, nSize, nd);
pref_allTrials_smallPupil = cell(nCon, nSize, nd);

% Added from old code
stat_resp_keep = cell(1, nd);

% Trial counts
trialCounts = cell(2, nd);

for id = 1:nd
    % Subset data to keep cells only
    data_dfof_trial_keep{id} = data_dfof_trial_match{id}(:, :, keep_cells);
    fullTC_keep{id} = cellTCs_match{id}(:, keep_cells);
    
    % Trial info for this day
    tCon = tCon_match{id}(:, 1:nTrials(id));
    tSize = tSize_match{id}(:, 1:nTrials(id));
    tDir = tDir_match{id}(:, 1:nTrials(id));
    
    % Condition indices
    stat_inds = find(~RIx{id});
    loc_inds = find(RIx{id});
    ind_stat_largePupil = intersect(stat_inds, find(PIx_stat{1, id}));
    ind_stat_smallPupil = intersect(stat_inds, find(PIx_stat{2, id}));
    
    % Initialize arrays for this day
    data_resp = zeros(nKeep, nDir, nCon, nSize, 2);
    stat_resp = zeros(nKeep, nDir, nCon, nSize);  % Added from old code
    temp_pref_responses_stat = zeros(nKeep, nCon, nSize);
    temp_pref_responses_loc = zeros(nKeep, nCon, nSize);
    temp_pref_responses_stat_largePupil = zeros(nKeep, nCon, nSize);
    temp_pref_responses_stat_smallPupil = zeros(nKeep, nCon, nSize);
    
    temp_tc_stat = nan((nOn+nOff), nKeep, nCon, nSize);
    temp_tc_loc = nan((nOn+nOff), nKeep, nCon, nSize);
    temp_tc_stat_largePupil = nan((nOn+nOff), nKeep, nCon, nSize);
    temp_tc_stat_smallPupil = nan((nOn+nOff), nKeep, nCon, nSize);
    
    % Initialize significance testing array
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
                
                % Significance test: response vs baseline
                [h(:, iDir, iCon, iSize), ~] = ttest(nanmean(data_dfof_trial_keep{id}(resp_win, ind, :), 1), ...
                                                     nanmean(data_dfof_trial_keep{id}(base_win, ind, :), 1), ...
                                                     'dim', 2, 'tail', 'right', ...
                                                     'alpha', 0.01./(nDir*nCon*nSize-1));
            end
        end
    end
    
    % Calculate responses for all conditions
    for iDir = 1:nDir
        ind_dir = find(tDir == dirs(iDir));
        for iCon = 1:nCon
            ind_con = find(tCon == cons(iCon));
            for iSize = 1:nSize
                ind_size = find(tSize == sizes(iSize));
                ind_temp = intersect(ind_dir, ind_con);
                ind = intersect(ind_temp, ind_size);
                ind_stat = intersect(ind, stat_inds);  % Added stationary trials only
                
                % Calculate mean response and SEM for all cells
                data_resp(:, iDir, iCon, iSize, 1) = squeeze(nanmean(nanmean(data_dfof_trial_keep{id}(resp_win, ind, :), 1), 2));
                data_resp(:, iDir, iCon, iSize, 2) = squeeze(std(nanmean(data_dfof_trial_keep{id}(resp_win, ind, :), 1), [], 2) ./ sqrt(length(ind)));
                
                % Calculate stationary responses for all conditions (added from old code)
                stat_resp(:, iDir, iCon, iSize) = squeeze(nanmean(nanmean(data_dfof_trial_keep{id}(resp_win, ind_stat, :), 1), 2));
            end
        end
    end
    
    % Calculate preferred responses using existing preferred directions
    for iCon = 1:nCon
        for iSize = 1:nSize
            for i = 1:nKeep
                % Get preferred direction for this cell (already calculated)
                temp_dir = dirs(prefDir_keep{id}(i));
                
                % Find trials matching preferred direction, contrast, and size
                dir_inds = find(tDir == temp_dir);
                con_inds = find(tCon == cons(iCon));
                size_inds = find(tSize == sizes(iSize));
                temp_trials = intersect(intersect(dir_inds, con_inds), size_inds);
                
                % Split by condition
                temp_trials_stat = intersect(temp_trials, stat_inds);
                temp_trials_loc = intersect(temp_trials, loc_inds);
                temp_trials_stat_largePupil = intersect(temp_trials, ind_stat_largePupil);
                temp_trials_stat_smallPupil = intersect(temp_trials, ind_stat_smallPupil);
                
                % Calculate time courses
                temp_tc_stat(:, i, iCon, iSize) = nanmean(data_dfof_trial_keep{id}(:, temp_trials_stat, i), 2);
                temp_tc_loc(:, i, iCon, iSize) = nanmean(data_dfof_trial_keep{id}(:, temp_trials_loc, i), 2);
                temp_tc_stat_largePupil(:, i, iCon, iSize) = nanmean(data_dfof_trial_keep{id}(:, temp_trials_stat_largePupil, i), 2);
                temp_tc_stat_smallPupil(:, i, iCon, iSize) = nanmean(data_dfof_trial_keep{id}(:, temp_trials_stat_smallPupil, i), 2);
                
                % Calculate preferred responses (mean of response window)
                temp_pref_responses_stat(i, iCon, iSize) = nanmean(nanmean(data_dfof_trial_keep{id}(resp_win, temp_trials_stat, i), 1));
                temp_pref_responses_loc(i, iCon, iSize) = nanmean(nanmean(data_dfof_trial_keep{id}(resp_win, temp_trials_loc, i), 1));
                temp_pref_responses_stat_largePupil(i, iCon, iSize) = nanmean(nanmean(data_dfof_trial_keep{id}(resp_win, temp_trials_stat_largePupil, i), 1));
                temp_pref_responses_stat_smallPupil(i, iCon, iSize) = nanmean(nanmean(data_dfof_trial_keep{id}(resp_win, temp_trials_stat_smallPupil, i), 1));
                
                % Track trial counts
                trialCounts{1, id} = [trialCounts{1, id}, length(temp_trials_stat)];
                trialCounts{2, id} = [trialCounts{2, id}, length(temp_trials_loc)];
            end
        end
    end
    
    % Test if largest size response > peak size response (adapted from old code)
    h_largeVsPeak = zeros(nKeep, 1);
    p_largeVsPeak = zeros(nKeep, 1);
    
    for iCell = 1:nKeep
        % Find peak size at preferred direction and highest contrast
        resp_by_size = squeeze(stat_resp(iCell, prefDir_keep{id}(iCell), end, :));
        [~, peakSize] = max(resp_by_size);
        
        % If peak size = largest size, no significant difference
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
            
            % Perform statistical test if we have sufficient trials
            if length(ind_large) >= 3 && length(ind_peak) >= 3
                resp_large = squeeze(nanmean(data_dfof_trial_keep{id}(resp_win, ind_large, iCell), 1));
                resp_peak = squeeze(nanmean(data_dfof_trial_keep{id}(resp_win, ind_peak, iCell), 1));
                [h_largeVsPeak(iCell), p_largeVsPeak(iCell)] = ttest2(resp_large, resp_peak, 'tail', 'right', 'alpha', 0.05);
            end
        end
    end
    
    % Store results for this day
    data_resp_keep{id} = data_resp;
    stat_resp_keep{id} = stat_resp;  % Added from old code
    h_keep{id} = h;
    h_largeVsPeak_keep{id} = h_largeVsPeak;  % Added from old code
    p_largeVsPeak_keep{id} = p_largeVsPeak;  % Added from old code
    
    pref_responses_stat{id} = temp_pref_responses_stat;
    pref_responses_loc{id} = temp_pref_responses_loc;
    pref_responses_stat_largePupil{id} = temp_pref_responses_stat_largePupil;
    pref_responses_stat_smallPupil{id} = temp_pref_responses_stat_smallPupil;
    
    tc_trial_avrg_stat{id} = temp_tc_stat;
    tc_trial_avrg_loc{id} = temp_tc_loc;
    tc_trial_avrg_stat_largePupil{id} = temp_tc_stat_largePupil;
    tc_trial_avrg_stat_smallPupil{id} = temp_tc_stat_smallPupil;
end

% Summary of t-test results
total_pass = sum(cellfun(@(x) sum(x == 1), h_largeVsPeak_keep));
total_cells = sum(cellfun(@length, h_largeVsPeak_keep));
fprintf('%d/%d cells pass large vs peak test (%.1f%%)\n', total_pass, total_cells, 100*total_pass/total_cells);

% Ask to plot
response = input('Plot boxplots? (y/n): ', 's');
if strcmpi(response, 'y')
    plotPeakVsLargeBoxplot(1, h_largeVsPeak_keep, prefDir_keep, stat_resp_keep, ...
                          data_dfof_trial_keep, tCon_match, tSize_match, tDir_match, ...
                          RIx, nTrials, dirs, cons, sizes, resp_win);
end
% Save timecourse data
save('tc_keep.mat', 'tc_trial_avrg_stat', 'tc_trial_avrg_loc', ...
     'tc_trial_avrg_stat_largePupil', 'tc_trial_avrg_stat_smallPupil');

clear temp_* ind_* dir_inds con_inds size_inds temp_trials* temp_dir


%% Normalized direction tuning
%this will normalize direction tuning so that the preferred direction is
%referred to as 0 for each cell. This only uses stationary trials.
% Note that there are still nDir directions, they are just
%rearranged. 
order = [1:nDir]; %make a list of indices for the directions
norm_dir_resp = cell(1,nd);

for id = 1:nd
    normDir_temp = nan(nKeep,nDir,nCon,nSize); %make an empty matrix to hold the data for this day
    for iCell = 1:nKeep
        index=prefDir_keep{id}(iCell);
        newOrder=circshift(order,((index-1)*-1));
        normDir_temp(iCell,:,:,:)=stat_resp_keep{id}(iCell,newOrder,:,:);
    end
    norm_dir_resp{id}=normDir_temp;

end


% Save response data
save('resp_keep.mat', 'stat_resp_keep', 'pref_responses_stat', 'pref_responses_loc', ...
     'pref_responses_stat_largePupil', 'pref_responses_stat_smallPupil', 'trialCounts', ...
     'conBySize_resp_stat_keep', 'conBySize_resp_loc_keep', 'h_keep', ...
     'h_largeVsPeak_keep', 'p_largeVsPeak_keep','norm_dir_resp');

%% Interneuron/Pyramidal Cell Relationship Analysis
% Initialize data structures
trialResp = cell(1, nd);        % Raw trial responses for each cell 
subTrialResp = cell(1, nd);     % Mean-subtracted trial responses (for noise correlation)
conditionMeans = cell(1, nd);   % Mean response for each stimulus condition

% Calculate trial responses and condition means
for id = 1:nd
    % Extract trial responses during stimulus period
    trialResp{id} = squeeze(mean(data_dfof_trial_keep{id}(stimStart:(stimStart+nOn-1), :, :), 1, 'omitmissing'));
    subTrialResp{id} = nan(size(trialResp{id}));
    conditionMeans{id} = nan(nDir, nCon, nSize, nKeep);
    
    % Get trial parameters
    tCon = tCon_match{id}(:, 1:nTrials(id));
    tDir = tDir_match{id}(:, 1:nTrials(id));
    
    % Calculate condition means and subtract from trial responses
    for iDir = 1:nDir
        ind_dir = find(tDir == dirs(iDir));
        
        for iCon = 1:nCon
            ind_con = find(tCon == cons(iCon));
            ind_dir_con = intersect(ind_dir, ind_con);
            
            for iSize = 1:nSize
                ind_size = find(tSize == sizes(iSize));
                ind_condition = intersect(intersect(ind_dir_con, ind_size), find(~RIx{id})); % Stationary trials only
                
                if ~isempty(ind_condition)
                    % Calculate condition mean for each cell
                    tempData = trialResp{id}(ind_condition, :);
                    cellMeans = mean(tempData, 1, 'omitmissing');
                    
                    % Subtract condition mean from individual trials
                    subTrialResp{id}(ind_condition, :) = tempData - cellMeans;
                    conditionMeans{id}(iDir, iCon, iSize, :) = cellMeans;
                end
            end
        end
    end
end

% Setup plotting preferences
response = input('Plot correlation figures as a sanity check? (y/n): ', 's');
doPlot = strcmpi(response, 'y') || strcmpi(response, 'yes');

if doPlot 
    numToPlot = min(5, nKeep);
    cellsToPlot = randperm(nKeep, numToPlot);
    disp(['Plotting enabled for ' num2str(numToPlot) ' randomly selected cells.']);
else
    disp('Plotting disabled.');
    cellsToPlot = [];
end

% Calculate correlations
noiseCorr = cell(1, nd);    % Noise correlations (trial-to-trial variability)
sigCorr = cell(1, nd);      % Signal correlations (stimulus tuning similarity)

for id = 1:nd
    noiseCorr{id} = nan(2, nKeep);  % [correlation; p-value]
    sigCorr{id} = nan(2, nKeep);
    
    % Reshape condition means for signal correlation analysis
    condMeansReshaped = reshape(conditionMeans{id}, nCon*nDir*nSize, nKeep);
    
    for iCell = 1:nKeep
        % Determine comparison cells based on cell type
        if ismember(iCell, red_cells_keep)
            % Red cell: compare to all green cells
            otherCells = green_cells_keep;
        else
            % Green cell: compare to other green cells (exclude self)
            otherCells = setdiff(green_cells_keep, iCell);
        end
        
        if ~isempty(otherCells)
            % Noise correlation calculation
            [R_noise, p_noise] = calculateCorrelation(subTrialResp{id}(:, iCell), ...
                                                      subTrialResp{id}(:, otherCells));
            noiseCorr{id}(:, iCell) = [R_noise; p_noise];
            
            % Optional plotting for noise correlation
            if doPlot && ismember(iCell, cellsToPlot)
                plotCorrelation(subTrialResp{id}(:, otherCells), subTrialResp{id}(:, iCell), ...
                               R_noise, iCell, 'NoiseCorr', id, pre);
            end
            
            % Signal correlation calculation
            [R_signal, p_signal] = calculateCorrelation(condMeansReshaped(:, iCell), ...
                                                        condMeansReshaped(:, otherCells));
            sigCorr{id}(:, iCell) = [R_signal; p_signal];
            
            % Optional plotting for signal correlation
            if doPlot && ismember(iCell, cellsToPlot)
                plotCorrelation(condMeansReshaped(:, otherCells), condMeansReshaped(:, iCell), ...
                               R_signal, iCell, 'SigCorr', id, pre);
            end
        end
    end
end

% Save results
save(fullfile(fn_multi, 'HT_pyr_relationship.mat'), 'conditionMeans', 'sigCorr', 'noiseCorr');



