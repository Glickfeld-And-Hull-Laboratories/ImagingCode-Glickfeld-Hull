clear all; clear global; 
close all
clc
ds = 'DART_expt_info'; %dataset info
dataStructLabels = {'contrastxori'};
rc =  behavConstsDART; %directories
eval(ds);
includeAllGCaMP = 0; 

if includeAllGCaMP == 1
    session_id = [19 21 23 25 27 30 33 36 39 41 43];
else % only include 8f
    session_id = [27 30 33 36 39 41 43];
end
experimentFolder = 'PV_YM90K';
fnroot = fullfile(rc.achAnalysis,experimentFolder);
fnout = fullfile(fnroot,'summary_analyses','baselineResponse_forYM');
% mkdir(fullfile(fnroot,'summary_analyses','baselineResponse_forYM'));
idx = 1;


%% 

% a table with the mean and std deviation of PV and Pyr responses at each contrast and each behavioral state
% only use cells that have running at all contrasts; ONLY large size

for idx = 1:length(session_id)
    this_day = session_id(idx);
    this_mouse = expt(this_day).mouse;
    this_date = expt(this_day).date;
    run_folder = expt(this_day).contrastxori_runs;
    this_time = expt(this_day).contrastxori_time;
    % pull first day data
    folder_path = fullfile(fnroot,this_mouse,this_date,run_folder);
    TC_path = fullfile(folder_path{1},"TCs.mat");
    load(TC_path,'npSub_tc');
    load(fullfile(folder_path{1},'input.mat'));
    load(fullfile(folder_path{1},'mask_cell.mat'),'mask_label');
    % Behav_Name = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' this_mouse '-' this_date '-' this_time{1} '.mat'];
    % calculate dfof from npSub_tc
    nOn = input.nScansOn;
    nOff = input.nScansOff;
    ntrials = size(input.tGratingDirectionDeg,2);
    nCells = size(npSub_tc,2);
    tc_trial = reshape(npSub_tc,nOn+nOff,ntrials,nCells);
    offFramesIdx = nOff/2+1:nOff;
    meanBase_trial = mean(tc_trial(offFramesIdx,:,:),1);
    dfof_trial = (tc_trial-meanBase_trial)./meanBase_trial; % nFrames x nTrials x nCells
    
    % find running
    loc_dat = input.wheelSpeedValues;
    wheelspd = zeros(length(loc_dat),1);
    for iTrial = 1:length(wheelspd)
      wheelspd(iTrial) = mean(loc_dat{iTrial});
    end
    runidx = find(wheelspd>2);
    runTrialBoolean = wheelspd > 2;
    % find contrast 
    tCon = cell2mat(input.tGratingContrast);
    contrasts = unique(tCon);
    nCon = length(contrasts);
    % align and calculate NOTES: input.tGratingDiameterDeg is size
    tSize = celleqel2mat_padded(input.tGratingDiameterDeg);
    nSize = length(unique(tSize));
    sizeBoolean = tSize == 1000;
    
    
    
    fprintf(['Finished extracting ' this_mouse '\n']);
end
