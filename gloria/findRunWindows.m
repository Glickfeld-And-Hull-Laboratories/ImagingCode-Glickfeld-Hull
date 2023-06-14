%This script is to analyze the RAW behavior data of movingDots stimulus/running behavior during two photon imaging
%The analysis will calculate average running speed, average duration of
%running periods, and average duration of stationary periods.
%This analysis will also graph speed across time, and make histograms of durations of each behavioral state.
%will save speed, cells during different behavioral states, median& std& ave of duration of each behavioral state.

%% Section I: set paths and create analysis folders for each session
%define the directory and files

clear;
close all;

fn_base = 'Z:\\All_Staff';
mworks_fn = fullfile(fn_base, 'Behavior\Data');
Aging_pn = fullfile(fn_base, 'home\ACh\Aging\Gloria');
tc_fn = fullfile(Aging_pn, '\analysis\2p'); 

CD = [Aging_pn '\eyeTCs'];
cd(CD);
eyeParameters = uigetfile('*.mat');
load(eyeParameters);

fprintf(['Loaded ' eyeParameters '\r\n'])

%% Section II : calculate for each frame,find relative frames for each behavioral state and plot speed.

mworks = fullfile(mworks_fn, ['data-' mouse '-' date '-' time '.mat']);
load(mworks);

timepoints = datetime({DOB; record_date}, 'InputFormat','yyyy-MM-dd');
age = caldiff(timepoints, 'weeks'); % calculate age at time of recording in months
age = split(age, 'weeks');

smooth_pupil_rad = smooth(pupil_rad);

    lenframe = 65;
    speed = wheel_speed; %wheelSpeedCalc does the same thing as calculate_speed_2P
    % this gives you the average for each frame (ave speed in number of units during that 65ms)
    speed = double(speed);
    speed(abs(speed)<4.9) = 0; %remove all ticks below 4.9 that is jitter/noise
    speed = max(speed,0); %set any remaining negative values to 0
    speed = smooth(speed);
  
    % find relative behavioral states and save to behavior analysis file
    smallestspd = ceil(1/lenframe*1000);%smallestspd in unit/second, quadrature taken every 1ms
    frm_maxGap = 5; % if the animal is still for less than 300ms during running, the running before and after the short still should still be counted as one part
    [frames,frames_stay_cell, frames_bf_cell, frames_run_cell, ...
        frames_move_cell] = findFrames_behavStates_2P(speed,smallestspd,frm_maxGap);

    % plot speed and save figure
    fig_speedtc = figure;
    plot(frames, speed); hold on;
    axis([0 54000 -40 60])
    title (['average speed every frame (65ms)' '  - ' mouse ' (' num2str(age) 'm' ')']);
    xlabel ('frames');
    ylabel ('speed (pulses/s)');
    %saveas(fig_speedtc,[Aging_pn '/' 'runOnset' '/' datemouserun '_speed.tif']);

%% Plot average TCs of speed and pupil during stationary/running transition

period = 30; %# of frames right before and after running--- probably not useful
befoRunStay = 30; %1s, # of frames that speed=0 before running
befoRun = 30; %# of frames to include before running starts
totalT1 = 45; %# total T = of frames before running onset + # of frames following running onset/ # of frames before running offset + # of frames after running offset
totalT2 = 45;
aftRunStay = 30; %2s,# of frames that speed = 0 after running
aftRun = 15; %# of frames to include after running ends

[frames_befo_run_cell,frames_aft_run_cell,frames_runTrigger_mat,frms_runoff_mat,frames_runoff_include, ...
    frames_run_mat] = findFrames_runWindows_2P(speed,frames_run_cell,period,befoRunStay,totalT1,totalT2,aftRunStay,befoRun,aftRun);
%% 
nWindows = size(frames_runTrigger_mat,2);
window_ext = zeros(15, nWindows);
runTrigger_speed = zeros(totalT1+15, nWindows);
runTrigger_pupil = zeros(totalT1+15, nWindows);

avg_runTrigger_speed = zeros(totalT1+15,1); %averaged across all run windows into a single average TC
avg_runTrigger_pupil = zeros(totalT1+15,1); %averaged across all run windows into a single average TC
std_err_speed = zeros(totalT1+15,1);
std_err_pupil = zeros(totalT1+15,1);

for i = 1:nWindows
    window_ext(:,i) = frames_runTrigger_mat(end,i)+1:frames_runTrigger_mat(end,i)+15; %we set the criteria as running for 1s, but still want to look at data 1s past that cutoff
    frames_runTrigger_mat_ext = cat(1,frames_runTrigger_mat,window_ext);
    runTrigger_speed(:,i) = speed(frames_runTrigger_mat_ext(:,i)); %access the avg speeds for all frames in each run-triggered window
    runTrigger_pupil(:,i) = smooth_pupil_rad(frames_runTrigger_mat_ext(:,i)); %access the pupil radius for all frames in each run-triggered window
end

baseline_pupil = mean(runTrigger_pupil(1:befoRun/2,:)); %create baseline pupil matrix by averaging from first half of each window's stationary period
runTrigger_pupil_dfof = (runTrigger_pupil-baseline_pupil)./baseline_pupil; %normalize pupil sizes to baseline to get dfof

avg_runTrigger_speed(:,1) = mean(runTrigger_speed(:,:),2); % averages speed across all windows for each frame in running window
avg_runTrigger_pupil(:,1) = mean(runTrigger_pupil_dfof(:,:),2,'omitnan'); % averages speed across all windows for each frame in running window

std_err_speed(:,1) = (std(runTrigger_speed,0,2))/(sqrt(nWindows)); % standard error across all running windows for each frame in window
std_err_pupil(:,1) = (std(runTrigger_pupil_dfof,0,2,'omitnan'))/(sqrt(nWindows)); % standard error across all running windows for each frame in window

save(fullfile(Aging_pn, "runOnset", [datemouserun '_runOnset.mat']), 'mouse', 'age', 'datemouserun', 'frames_runTrigger_mat_ext', 'avg_runTrigger_pupil', 'avg_runTrigger_speed', 'std_err_speed', 'std_err_pupil', 'nWindows', 'runTrigger_speed', 'runTrigger_pupil_dfof')

%% 
%create vector for time axis
t = [1:size(avg_runTrigger_speed,1)];
t = (t-30)./15; %re-centers start of running (frame 30) as 0 and converts frames to seconds

left_color = [0 0 0];
right_color = [0 0 0];
set(figure,'defaultAxesColorOrder',[left_color; right_color]);
shadedErrorBar(t,avg_runTrigger_speed(:),std_err_speed(:),'lineProps','-r'); % plot average time course with error bars of speed across all running windows
axis([t(1) t(end) -2 15])
vline(0,'black'); % indicate when running starts
hold on
ylabel("Speed (pulses/s)")
xlabel('Time (s)')
title(['Avg time course for running onset' '  - ' mouse ' (' num2str(age) 'wks' ')'], "FontSize", 11)

yyaxis right
axis([t(1) t(end) -0.1 0.8])
ylabel("Pupil Fractional Change (from baseline)")
shadedErrorBar(t,avg_runTrigger_pupil(:),std_err_pupil(:), 'lineProps','-b'); % plot average time course with error bars of speed across all running windows

%% Subsampling windows to determine nWindow cutoff for inclusion

iterations = floor(nWindows/5);
std_subsample_speed = zeros(60,iterations);
std_subsample_pupil = zeros(60,iterations);
sub_nWindows = zeros(1,iterations);

[row col] = subplotn(iterations); 
left_color = [0 0 0];
right_color = [0 0 0];

set(figure,'defaultAxesColorOrder',[left_color; right_color]);
for i = 1:iterations
    sub_nWindows(1,i) = i*5;
    subsample = randsample(nWindows,sub_nWindows(1,i));
    subsample_speed = runTrigger_speed(:,subsample);
    subsample_pupil = runTrigger_pupil_dfof(:,subsample);

    avg_subsample_speed(:,1) = mean(subsample_speed(:,:),2);
    avg_subsample_pupil(:,1) = mean(subsample_pupil(:,:),2, 'omitnan');

    std_subsample_speed(:,i) = std(subsample_speed,0,2);
    std_subsample_pupil(:,i) = std(subsample_pupil,0,2,'omitnan');

    std_err_subsample_speed(:,1) = std_subsample_speed(:,i)/(sqrt(sub_nWindows(1,i)));
    std_err_subsample_pupil(:,1) = std_subsample_pupil(:,i)/(sqrt(sub_nWindows(1,i)));

    subplot(row,col,i)
    shadedErrorBar(t,avg_subsample_speed(:),std_err_subsample_speed(:),'lineProps','-r'); % plot average time course with error bars of speed across all running windows
    axis([t(1) t(end) -2 20])
    vline(0,'black'); % indicate running onset
    hold on
    ylabel("Speed (pulses/s)")
    xlabel('Time (s)')
    title([mouse ' (' num2str(age) 'm' ')' ' - ' num2str(sub_nWindows(1,i)) ' windows'], "FontSize", 10)
 
    yyaxis right
    axis([t(1) t(end) -0.2 0.8])
    ylabel("Pupil Fractional Change", 'FontSize', 9)
    shadedErrorBar(t,avg_subsample_pupil(:),std_err_subsample_pupil(:), 'lineProps','-b'); % plot average time course with error bars of speed across all running windows
    
    sgtitle('Average time course for subsampling');
end

avg_std_subsample_speed = zeros(1,iterations);
avg_std_subsample_pupil = zeros(1,iterations);

avg_std_subsample_speed(:) = mean(std_subsample_speed(:,:),1);
avg_std_subsample_pupil(:) = mean(std_subsample_pupil(:,:),1);

print(fullfile(Aging_pn, "Figures", [datemouserun '_subsample.pdf']),'-dpdf','-fillpage')

figure;
subplot(1,2,1)
    scatter(sub_nWindows,avg_std_subsample_speed)
    hold on
    axis([0 nWindows 0 4.0])
    xlabel('# Windows Subsampled')
    ylabel("STD")
    title(['STD of Speed across Subsamples'], "FontSize", 9)
subplot(1,2,2)
    scatter(sub_nWindows,avg_std_subsample_pupil)
    hold on 
    axis([0 nWindows 0 0.3])
    xlabel('# Windows Subsampled')
    ylabel("STD")
    title(['STD of Pupil Change across Subsamples'], "FontSize", 9)

print(fullfile(Aging_pn, "Figures", [datemouserun '_subsample_STD.pdf']),'-dpdf','-fillpage')

