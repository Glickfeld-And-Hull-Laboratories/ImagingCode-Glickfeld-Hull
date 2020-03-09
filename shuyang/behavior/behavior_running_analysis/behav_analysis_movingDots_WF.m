%This script is to analyze the RAW behavior data of movingDots stimulus/running behavior during wide field imaging
%The analysis will calculate average running speed, average duration of
%running periods, and average duration of stationary periods.
%This analysis will also graph speed across time, and make histograms of durations of each behavioral state.
%will save speed, cells during different behavioral states, median& std& ave of duration of each behavioral state.

%% Section I: set paths and create analysis folders for each session
%define the directory and files
clear;
folder = 'Z:\Data\Behv_MovingDots\behavior_raw_WF\';
sessionID = '1029-190618';% this  variable name is confusing, this session ID is just tha date and the subject#, 
%there might be more than 1 sessions on a single subject on the same day
filename = dir([folder 'data-i' '*' sessionID  '*' ]);
for i = 1: size(filename,1)
    behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' sessionID '_' num2str(i)];
    if ~exist(behav_dest)
        mkdir(behav_dest);
    end
end

%% Section II : calculate for each frame,find relative frames for each behaviroal state and plot speed.
for i = 1:size(filename,1)
    file = [folder, filename(i).name];
    load(file);
    % later analysis will need to know if this session is a reverse session, and if yes, the frames of the reverse stimuli. 
    cReverse = input.cReverse;
    cReverse_vec = [cReverse{:}];
    speed = calculate_speed(input);
    % find relative behavioral states and save to behavior analysis file
    [frames,frames_stay_cell, frames_bf_cell, frames_run_cell, frames_move_cell] = findFrames_behavStates(speed);
    save([behav_dest '\' sessionID '_' num2str(i) '_behavAnalysis.mat' ],...
        'cReverse_vec','speed','frames','frames_stay_cell','frames_bf_cell',...
        'frames_run_cell','frames_move_cell');
    % plot speed and save figure
    fig_speedtc = figure;
    plot(frames, speed); hold on;
    if isempty(cReverse_vec) == 0
    vline(cReverse_vec, 'r');
    end
    title (['average speed every frame(100ms)', '  ', sessionID '-' num2str(i)]);
    xlabel ('frames');
    ylabel ('speed(pulses/s)');
    saveas(fig_speedtc,[behav_dest '\' sessionID '_' num2str(i) '_speed']);
end
% check if found frames are what I want to find
figure;plot(speed);
hold on;
plot(cell2mat(frames_bf_cell), 10*ones(1,length(cell2mat(frames_bf_cell))),'r.');

%% draw distribution of running duration
%num_frames_run_all = {};
num_frames_run = [];
for i =  1:size(filename,1)
    for n = 1:size(frames_run_cell,2)
    num_frames_run = [num_frames_run size(frames_run_cell{n},2)];
    end
    aveDuration_runFrames = mean(num_frames_run);
    medianDura_runFrames = median(num_frames_run);
    stdDura_runFrames = std(num_frames_run);
    save([behav_dest '\' sessionID '_' num2str(i) '_behavAnalysis.mat' ],...
        'aveDuration_runFrames', 'medianDura_runFrames','stdDura_runFrames', '-append');
    % plot distribution of all running windows' duration
    fig_runDuradist = figure;
    edges = (0:2:max(num_frames_run));
    histogram(num_frames_run,edges);
    title (['distribution of running duration', '  ',filename(i).name]);
    xlabel ('time length (number of frames)');
    ylabel ('frequency');
    vline(medianDura_runFrames, 'r');
    saveas(fig_runDuradist,[behav_dest '\' sessionID '_' num2str(i) '_runDuraDist']);
end


%% Section IV: Calculate average stationary duration:
num_frames_stay = [];
for i = 1:size(filename,1)
    for k = 1: size(frames_stay_cell,2)
        num_frames_stay = [num_frames_stay size(frames_stay_cell{k},2)];
    end
    aveDuration_stayFrames = mean(num_frames_stay);
    medianDura_stayFrames = median(num_frames_stay);
    stdDura_stayFrames = std(num_frames_stay);
    save([behav_dest '\' sessionID '_' num2str(i) '_behavAnalysis.mat' ],...
        'aveDuration_stayFrames', 'medianDura_stayFrames','stdDura_stayFrames', '-append');
    %plot distribution of number of frames in each stationary state
    fig_staydist = figure;
    edges = (0:5:max(num_frames_stay));
    histogram(num_frames_stay,edges);
    title (['distribution of stationary duration', '  ',filename(i).name]);
    xlabel ('time length (number of frames)');
    xlim([1 500]);
    ylabel ('frequency');
    vline(medianDura_stayFrames, 'r');
    saveas(fig_staydist,[behav_dest '\' sessionID '_' num2str(i) '_stayDistribution' ]);
end


%% Section V: calculate average speed and ave&median duration of running above the runThreshold(mean)
% calculate the average speed to get running threshold
for i = 1:size(filename,1)
    frames_run = frames(cell2mat(frames_run_cell));
    speed_run = speed(frames_run);
    aveRunSpeed = mean(speed_run);
    %threshold of reverse stimuli.
    runThreshold = roundn(aveRunSpeed,1); %round averun to the nearest multiple of 10.Need to do this because all speed values are multiple of 10, so need a multiple of 10 for setting the
   
    % find above threshold running durations
    dura_runAbvThre = [];
    for p = 1: size(frames_run_cell,2)
        dura_runAbvThre(p) = sum((speed(frames_run_cell{p}) >= runThreshold));
    end
    aveDura_runAbvThre = mean(dura_runAbvThre);
    medianDura_runAbvThre = median(dura_runAbvThre);
    stdDura_runAbvThre = std(dura_runAbvThre);
    save([behav_dest '\' sessionID '_' num2str(i) '_behavAnalysis.mat' ],...
        'runThreshold', 'aveDura_runAbvThre','medianDura_runAbvThre','stdDura_runAbvThre', '-append');
    %plot distribution of number of frames in each stationary state
    fig_ABTspeeddist = figure;
    edges = (0:2:max(dura_runAbvThre));
    histogram(dura_runAbvThre,edges);
    title (['distribution of abvThreRun duration', '  ',filename(i).name]);
    xlabel ('time length (number of frames)');
    ylabel ('frequency');
    vline(runThreshold, 'r');
    saveas(fig_ABTspeeddist,[behav_dest '\' sessionID '_' num2str(i) '_ATR Distribution' ]);
end

