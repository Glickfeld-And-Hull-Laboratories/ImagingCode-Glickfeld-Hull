% cells are already registered, segmented
% read in timecourse data
% read in behavioral data

% do this for each run I want to compare

clear all; clear global; close all

%identifying animal and run
date = '210501';
imgFolder = '001';
time = '1201';
mouse = 'CC05';
frame_rate = 60; %enter the frame rate, or I can edit this to enter the stimulus duration


%setting my paths
fn_base = 'Z:\home\celine\Analysis\2p_analysis\';
fn = fullfile(fn_base,mouse,date,imgFolder);
cd(fn)

run_str = catRunName(imgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];
TC_file = [datemouserun '_TCs.mat'];
load(TC_file);

beh_prefix = strcat('Z:\Behavior\Data\data-i''');
beh_file = [beh_prefix mouse '''-' date '-' time '.mat'];
load(beh_file); %load the mworks behavioral file
behData = input;
clear beh_prefix beh_file input TC_file

nCells = size(npSub_tc,2);
nOn = behData.nScansOn;
nOff = behData.nScansOff;
ntrials = size(behData.tGratingDirectionDeg,2); %this is a cell array with one value per trial, so length = ntrialstDir = celleqel2mat_padded(behData.tGratingDirectionDeg); %transforms cell array into matrix (1 x ntrials)
tDir = celleqel2mat_padded(behData.tGratingDirectionDeg); %transforms cell array into matrix (1 x ntrials)
Dirs = unique(tDir);
nDirs = length(Dirs);

data_tc_trial = reshape(npSub_tc, [nOn+nOff,ntrials,nCells]);
data_f_trial = mean(data_tc_trial(nOff/2:nOff,:,:),1);
data_dfof_trial = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial, data_f_trial), data_f_trial);

%% plot time courses
tc_cell_avrg = mean(data_dfof_trial,3);%average pver cells, one row per trial
tc_trial_avrg = squeeze(mean(data_dfof_trial,2));%average over trials, one row per cell
tc_cell_trial_avrg = mean(tc_cell_avrg,2);%average over trials and cells

figure;
plot(tc_trial_avrg, 'LineWidth',.005,'color',[.25 .25 .25]);
hold on;
plot(tc_cell_trial_avrg, 'LineWidth',2, 'color','k');
hold on;
vline(nOff,'g')
title('Timecourses after np subtraction');
hold off

%% remove outliers based on max response? 

%find the max trial-averaged dfof for each cell
max_values = max(tc_trial_avrg);
mean_max = mean(max_values);
sd_max = std(max_values);
outliers_max = find(max_values > mean_max + (3*sd_max));

tc_trial_avrg(:,outliers_max)=[];
dfof_cleaned = data_dfof_trial;
dfof_cleaned(:,:,outliers_max)=[]; %un-averaged time courses with outliers removed


figure;
plot(tc_trial_avrg, 'LineWidth',.005,'color',[.25 .25 .25]);
hold on;
plot(tc_cell_trial_avrg, 'LineWidth',2, 'color','k');
hold on;
vline(nOff,'g')
title('Timecourses after removing outliers');
hold off
%% looking only at time right before and during stim
%find how many frames is half a second before the stimulus turns on
%half_sec = (frame_rate * .5)
%start_time = nOff - half_sec;
%make a matrix of cleaned dfof time courses for the period half a second
%before the stim through the stim
dfof_stim_win = dfof_cleaned((nOff+1):(nOff+nOn),:,:);
%same as above but activity of each cell averaged over the trials, to get
%started
avrg_dfof_stim_win = tc_trial_avrg((nOff+1):(nOff+nOn),:,:);

figure;
t=(1:size(avrg_dfof_stim_win,1))-half_sec;
plot(avrg_dfof_stim_win, 'LineWidth',.005,'color',[.25 .25 .25]);
%vline(half_sec,'g')


% What frame does the response peak at? I is the index where the peak
% occures
[M I] = max(avrg_dfof_stim_win);
mean(I-half_sec) %the average number of frames after stim onset that the TC peaks

% What frame does the response hit its lowest point?
[Mmin Imin] = min(avrg_dfof_stim_win);
mean(Imin-half_sec) 
%negative value indicates that cells usually do not fall below baseline in the sitm on period


deriv = diff(avrg_dfof_stim_win);
figure;
plot(deriv, 'LineWidth',.005,'color',[.25 .25 .25]);
% was not too informative, too noisy

% maybe look at how the dfof at 
