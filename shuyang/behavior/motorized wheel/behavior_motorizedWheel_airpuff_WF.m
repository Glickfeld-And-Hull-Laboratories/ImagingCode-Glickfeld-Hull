%%This script is to extract frame index of different speed modes.
% during motorized treadmill experiment in 2P imaging
% runoff_fast_mat, runoff_slow_mat,runtrig_fast_mat, runtrig_slow_mat,
% run_off_mat, runtrig_mat are maxtrix containing
% frame index for running triggered average and runing offset average.
% each row is a trial 
% matrix for airpuff trigger during different behavioral states:
% airpuff_fast_mat,airpuff_slow_mat,airpuff_stay_mat,airpuff_run_mat
% lastly, stay_woAipuff includes all stationary frames excluding the first second 
% of stationary trial and 1s after airpuff delivery.

%%
clear;
%folder = 'Z:\home\shuyang\Data\behavior\motorized_WF\';
folder = 'Z:\home\shuyang\Data\behavior\motorized_WF_drug\';
days = '1911-220126-1629';% if there're more than 1 sessions on a single subject on the same day, put the time in here after the date
%1728
filename = dir([folder 'data-i' '*' days  '*' ]);
for i = 1: size(filename,1)
    %behav_dest{i} = ['Z:\home\shuyang\Analysis\motorizedWheel_Analysis\WF\behavioral_analysis\' days '_' num2str(i)];
    behav_dest{i} = ['Z:\home\shuyang\Analysis\motorizedWheel_Analysis\WF_drug\behavioral_analysis\' days '_' num2str(i)];
    if ~exist(behav_dest{i})
        mkdir(behav_dest{i});
    end
end

%% matrix for running trig and running offset, speed change, and airpuff trigger
%--------------------------------------running-----------------------------------------------------------
for i = 1: size(filename,1)
    % load data
    file = [folder, filename(i).name];
    load(file);
    countTime = cell2mat(input.counterTimesUs);
    countValue = input.counterValues;
    nstart = find(countValue{1} == 1, 1, 'last');
    countValue{1} = countValue{1}(nstart:end);
    %countValue = countValue(1:10);
    
    airpuffsOn = input.cAirpuffsOn;
    airpuffsOn{1} = airpuffsOn{1}(2:end); %there's always a zero at the beginning of the experiment, so delete this
    %airpuffsOn = airpuffsOn(1:10);
    airpuffsOn_vec = cell2mat(airpuffsOn);
    airpuffsOn_vec = double(airpuffsOn_vec);
    commandVolt = cell2mat(input.commandVoltage);
    %commandVolt = commandVolt(1:10);
    
    % if all running trials and stationary trials are the same length:
    ntrials = length(countValue);
    % when commnadVolt == 0, speed == 0
    % when commandVolt == 2, speed == 9.6-12cm/s
    % when commandVolt == 3, speed == 14.4-18cm/s
    % count how many trials there are for each condition
    nslowTrials = sum(commandVolt == 2);
    nfastTrials = sum(commandVolt == 3);
    nstayTrials = sum(commandVolt == 0);
    
    t1 = 0;
    t2 = 0;
    t3 = 0;
    t4 = 0;
    t5 = 0;
    t6 = 0;
    
    framerateHz = 10;
    secondsafter = 3;
    secondsbefore = 1;
    nframesafter = secondsafter*framerateHz;
    nframesbefore = secondsbefore*framerateHz-1;
    len = nframesafter+nframesbefore+1;
    runtrig_slow_mat = zeros(nslowTrials,len); % each row is a trial
    runoff_slow_mat = zeros(nslowTrials,len);
    runtrig_fast_mat = zeros(nfastTrials,len);
    runoff_fast_mat = zeros(nfastTrials,len);
    spd_decrease_mat = zeros(nfastTrials,len);
    spd_increase_mat = zeros(nfastTrials,len);
    %run trig: stationary ---> running trial
    %run off: running trial ---> stationary
    %speed change: increase & decrease
    
    for t = 1:ntrials-1
        if     commandVolt(t) == 0 && commandVolt(t+1) == 2 % runtrig slow
            t1 = t1 + 1;
            runtrig_slow_mat(t1,:) = countValue{t}(end)-nframesbefore:countValue{t+1}(nframesafter); % each row is a trial
        elseif commandVolt(t) == 0 && commandVolt(t+1) == 3 % runtrig fast
            t2 = t2 + 1;
            runtrig_fast_mat(t2,:) = countValue{t}(end)-nframesbefore:countValue{t+1}(nframesafter);
        elseif commandVolt(t) == 2 && commandVolt(t+1) == 0 % runoff slow
            t3 = t3 + 1;
            runoff_slow_mat(t3,:) = countValue{t}(end)-nframesbefore:countValue{t+1}(nframesafter);
        elseif commandVolt(t) == 3 && commandVolt(t+1) == 0 % runoff fast
            t4 = t4 + 1;
            runoff_fast_mat(t4,:) = countValue{t}(end)-nframesbefore:countValue{t+1}(nframesafter);
        elseif commandVolt(t) == 3 && commandVolt(t+1) == 2 % speed decrease
            t5 = t5 + 1;
            spd_decrease_mat(t5,:) = countValue{t}(end)-nframesbefore:countValue{t+1}(nframesafter);
        elseif commandVolt(t) == 2 && commandVolt(t+1) == 3 % speed increase
            t6 = t6 + 1;
            spd_increase_mat(t6,:) = countValue{t}(end)-nframesbefore:countValue{t+1}(nframesafter);
        end
    end
    % delete rows with zeros because not all fast/slow speed trials will meet the criteria in the for loop
    runtrig_slow_mat = runtrig_slow_mat(any(runtrig_slow_mat,2),:);
    runtrig_fast_mat = runtrig_fast_mat(any(runtrig_fast_mat,2),:);
    runoff_slow_mat = runoff_slow_mat(any(runoff_slow_mat,2),:);
    runoff_fast_mat = runoff_fast_mat(any(runoff_fast_mat,2),:);
    spd_decrease_mat = spd_decrease_mat(any(spd_decrease_mat,2),:);
    spd_increase_mat = spd_increase_mat(any(spd_increase_mat,2),:);
    
    runtrig_mat = cat(1,runtrig_slow_mat,runtrig_fast_mat);
    runoff_mat = cat(1,runoff_slow_mat,runoff_fast_mat);
    
    %--------------------------------------------airpuff-----------------------------------------------------
    t1 = 0;
    t2 = 0;
    t3 = 0;
    airpuffsOn_stay = cell(1,nstayTrials);
    airpuffsOn_slow = cell(1,nslowTrials);
    airpuffsOn_fast = cell(1,nfastTrials);
    run_slow_vec = [];
    run_fast_vec = [];
    stay_vec = [];
    for t = 1: ntrials
        if     commandVolt(t) == 0
            t1 = t1 + 1;
            airpuffsOn_stay{t1} = airpuffsOn{t}; %each trial has 2 airpuf stimuli
            stay_vec = cat(2,stay_vec,countValue{t}(21:end));% get rid of the transition time in between speed modes, so just delete the first 2 seconds
        elseif commandVolt(t) == 2
            t2 = t2 + 1;
            airpuffsOn_slow{t2} = airpuffsOn{t};
            run_slow_vec = cat(2,run_slow_vec,countValue{t}(21:end));
        elseif commandVolt(t) == 3
            t3 = t3 + 1;
            airpuffsOn_fast{t3} = airpuffsOn{t};
            run_fast_vec = cat(2,run_fast_vec,countValue{t}(21:end));
        end
    end
    airpuffsOn_stay = cell2mat(airpuffsOn_stay);
    airpuffsOn_slow = cell2mat(airpuffsOn_slow);
    airpuffsOn_fast = cell2mat(airpuffsOn_fast);
    
    
    seconds_around_air = 1;
    len2 = framerateHz*seconds_around_air*2; % taking equal amount of time before and after airpuff delivery
    nframes_around_air = framerateHz*seconds_around_air;
    airpuff_slow_mat = zeros(2*nslowTrials,len2); % each row is an airpuff trial (each running trial has 2 airpuffs)
    airpuff_fast_mat = zeros(2*nfastTrials,len2);
    airpuff_stay_mat = zeros(2*nstayTrials,len2);
    puffresp_slow_vec = [];
    puffresp_fast_vec = [];
    puffresp_stay_vec = [];
    for a = 1:length(airpuffsOn_slow)
        airpuff_slow_mat(a,:) =  airpuffsOn_slow(a)-nframes_around_air+1:airpuffsOn_slow(a)+nframes_around_air; % the 30th column is the frame when airpuff happens
        puffresp_slow_vec = cat(2,puffresp_slow_vec,airpuffsOn_slow(a):airpuffsOn_slow(a)+3);
    end
    for a = 1:length(airpuffsOn_fast)
        airpuff_fast_mat(a,:) =  airpuffsOn_fast(a)-nframes_around_air+1:airpuffsOn_fast(a)+nframes_around_air; 
        puffresp_fast_vec = cat(2,puffresp_fast_vec,airpuffsOn_fast(a):airpuffsOn_fast(a)+3);
    end
    for a = 1:length(airpuffsOn_stay)
        airpuff_stay_mat(a,:) =  airpuffsOn_stay(a)-nframes_around_air+1:airpuffsOn_stay(a)+nframes_around_air; 
        puffresp_stay_vec = cat(2,puffresp_stay_vec,airpuffsOn_stay(a):airpuffsOn_stay(a)+3);
    end
    airpuff_run_mat = cat(1,airpuff_slow_mat,airpuff_fast_mat);
    runslow_nopuff = setdiff(run_slow_vec,puffresp_slow_vec);
    runfast_nopuff = setdiff(run_fast_vec,puffresp_fast_vec);
    stay_nopuff = setdiff(stay_vec,puffresp_stay_vec);
    
    save([behav_dest{i} '\' days '_' num2str(i) '_behavAnalysis.mat' ],...
        'airpuffsOn','airpuffsOn_vec', 'runtrig_slow_mat','runtrig_fast_mat','runoff_slow_mat',...
        'runoff_fast_mat','spd_decrease_mat','spd_increase_mat','runtrig_mat',...
        'runoff_mat','airpuffsOn_stay','airpuffsOn_slow','airpuffsOn_fast',...
        'airpuff_slow_mat','airpuff_fast_mat','airpuff_stay_mat','airpuff_run_mat',...
        'puffresp_slow_vec','puffresp_fast_vec','puffresp_stay_vec','runslow_nopuff',...
        'runfast_nopuff','stay_nopuff','stay_vec','run_slow_vec','run_fast_vec');
end
%% write a vector that has commandvoltages and is the same length as frames 
for i = 1: size(filename,1)
    % need a for loop
    speed_inVolt  = [];
    for v = 1: length(commandVolt)
        speed_inVolt(countValue{v}) = commandVolt(v);
    end
    
    figure; plot(speed_inVolt);
    ylim([-0.5 3.5]);ylabel('Volts');xlabel('frames');
    vline(airpuffsOn_vec,'red');
    
    savefig([behav_dest{i} '\' days '_' num2str(i) '_speed.fig']);
    
    save([behav_dest{i} '\' days '_' num2str(i) '_behavAnalysis.mat' ],...
        'speed_inVolt','-append');
end

% %% frames during stationary without airpuff ---- useful for deconvolution, sort out bad cells
% for i = 1: size(filename,1)
%     stay_all_frames = [];
%     for c = 1: length(commandVolt)
%         if commandVolt(c) == 0
%             stay_all_frames = cat(2,stay_all_frames,countValue{c}(31:end)); % delete the first second because the mice might be transitioning from the last state
%         end
%     end
%     airpuff_trig_stay = [];
%     for a = 1: length(airpuffsOn_stay)
%         airpuff_trig_stay = cat(2,airpuff_trig_stay,airpuffsOn_stay(a)-1:airpuffsOn_stay(a)+30);
%     end
%     stay_woAirpuff = setdiff(stay_all_frames,airpuff_trig_stay);
%     stay_woAirpuff = double(stay_woAirpuff);
%     
%     save([behav_dest{i} '\' sessionID '_' num2str(i) '_behavAnalysis.mat' ],...
%         'stay_woAirpuff','-append');
% end
% 


