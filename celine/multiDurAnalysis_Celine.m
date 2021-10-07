% to analyze individual day's data
clear all; clear global; close all

%identifying animal and run
date = '211004';
imgFolder = '002';
time = '1442';
mouse = 'i2015';

%setting my paths
fnIn_base = 'Z:\home\Celine\Analysis\2p_analysis\';
fnOut_base = 'Z:\home\Celine\Analysis\2p_analysis\';

fnIn = fullfile(fnIn_base,mouse,date,imgFolder);
fnOut = fullfile(fnOut_base,mouse,date,imgFolder);
mkdir(fnOut);
cd(fnIn);

run_str = catRunName(imgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];

%% load daata
%data_tc =load(fullfile([datemouserun '_TCs.mat']));
load(fullfile([datemouserun '_trial_TCs.mat']));
load('input.mat');

stimOneTime = celleqel2mat_padded(input.tStimOneGratingOnTimeMs); %duration for each trial
stimOneTimes = unique(stimOneTime); %list of durations used
nTime = length(stimOneTimes); %how many different durations were used
nTrials = length(cStimOne);
ntrials = size(input.tThisTrialStartTimeMs,2);