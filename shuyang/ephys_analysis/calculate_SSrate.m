
%{
%% get and save master 8 times, will later use this to align the ephys data and behavior data
% load data, drag the file into matlab
channel = '191210_master8'; % used spikes2 to sort only master 8 pulses in this channel, so the times are the times for each pulse, start of the ephys recording is 0.
spk2_out_all = converted_ns6_191210_004_SJ_Ch37;
times_m8= spk2_out_all.times; %in seconds
time_diff = diff(times_m8);
logical(sum(round(time_diff,1)==0.1)==length(times_m8)-1); %check if all pulses are 0.1s apart from the beginning to the end
analysis_dest = 'Y:\home\nathan\DATA\S-probe\Awake\useful_ephys\analysis\';
save([analysis_dest channel '.mat'],'times_m8','');
%}

%% load ephys neuronal data and calculate # of simple spikes per 100ms
%channels = {'191206_Ch16','191210_D004_Chan13','191213_D002_Ch15',...
 %   '191213_D002_Ch18','191216_Ch6_MM'};
analysis_dest = 'Y:\home\nathan\DATA\S-probe\Awake\useful_ephys\analysis\';
m8channel = '191216_master8';
master8 = load([analysis_dest m8channel '.mat']);
times_m8 = master8.times_m8;
channel = '191216_Ch6_MM';
PC = load([analysis_dest channel '.mat']);
allspk_times = PC.spk_time; %in seconds
% find how many spikes there are between each pulse, this gives you # of
% spikes per 100ms. this firing rate is having the same resolution as
% speed, and they start at the same time. 
%(during recording, ephys is started first, then MWroks, then master8. MWorks is waiting for the first master8 pulse to start, just like it's waiting for the pulse from the camera during imaging. So when counter ==1 in behavior, it's also the time of the first pulse in the master8 channel.)
nspk_100ms = zeros(1,length(times_m8));
for t = 1:length(times_m8)
    if t == 1
        nspk_100ms(t) = length(find(allspk_times>times_m8(t)-0.1 & allspk_times <= times_m8(t)));
    else
        nspk_100ms(t) = length(find(allspk_times>times_m8(t-1) & allspk_times <= times_m8(t)));
    end
end
save([analysis_dest channel '.mat'],'nspk_100ms','-append');


%% calculate firing rate per second and see where the FR drops during recording
%channels = {'191206_Ch16','191210_D004_Chan13','191213_D002_Ch15',...
 %   '191213_D002_Ch18','191216_Ch6_MM'};
channel = '191216_Ch6_MM';
analysis_dest = 'Y:\home\nathan\DATA\S-probe\Awake\useful_ephys\analysis\';
PC = load([analysis_dest channel '.mat']);
nspk_100ms = PC.nspk_100ms;
nspk_s = zeros(1,floor(length(nspk_100ms)/10));
for i = 1:floor(length(nspk_100ms)/10)
    nspk_s(i) = sum(nspk_100ms((i-1)*10+1:10*i));
end

SSrate = figure;
plot(nspk_s);
title(channel);
xlabel('time(s)');
ylabel('simple spike rates(Hz)');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
SSrate.Units = 'centimeters';
SSrate.Position = [1 1 15 5];
fig_name = ['SSrate_everysec_' channel];
path = [analysis_dest 'figures\']; 
%orient(dfOvF_behavStates,'landscape');
print(SSrate,[path,fig_name],'-r600','-dpdf');
%saveas(SSrate,[path,fig_name]);


% cut nspk_100ms based on the firing rate plot above and save new nspk_100ms_cut
tcut = 1437; % the second where the FR drop down
nspk_100ms_cutted = nspk_100ms(1:tcut*10);
save([analysis_dest channel '.mat'],'nspk_100ms_cutted','tcut','nspk_s','-append');


%% load behavior data and calculate firing rate during running vs. stationary
%channels = {'191206_Ch16','191210_D004_Chan13','191213_D002_Ch15',...
 %   '191213_D002_Ch18','191216_Ch6_MM'};
%sessions = {'1050-191206','1045-191210','1045-191213','1046-191216'};
behav_dest = 'Y:\home\nathan\DATA\S-probe\Awake\useful_ephys\analysis\behavior';
sessionID = '1046-191216';
behav_output = load([behav_dest '\' sessionID '_behavAnalysis.mat' ]);
frames_run = cell2mat(behav_output.frames_run_cell);
frames_stay = cell2mat(behav_output.frames_stay_cell);
channel = '191216_Ch6_MM';
analysis_dest = 'Y:\home\nathan\DATA\S-probe\Awake\useful_ephys\analysis\';
PC = load([analysis_dest channel '.mat']);
nspk_100ms_cutted = PC.nspk_100ms_cutted;
tcut = PC.tcut;
frames_run_cut = frames_run(frames_run<tcut*10);
frames_stay_cut = frames_stay(frames_stay<tcut*10);
nspk_100ms_stay = nspk_100ms_cutted(frames_stay_cut);
nspk_100ms_run = nspk_100ms_cutted(frames_run_cut);
%firing rate = total number of spikes/time length
SS_stay = sum(nspk_100ms_stay)*10/length(frames_stay_cut);
SS_run = sum(nspk_100ms_run)*10/length(frames_run_cut);
save([analysis_dest channel '.mat'],'SS_stay','SS_run','-append');
save([behav_dest '\' sessionID '_behavAnalysis.mat'],'frames_run_cut','frames_stay_cut','-append');

%% SS aligned to running onset and offset
%channels = {'191206_Ch16','191210_D004_Chan13','191213_D002_Ch15',...
 %   '191213_D002_Ch18','191216_Ch6_MM'};
%sessions = {'1050-191206','1045-191210','1045-191213','1046-191216'};

behav_dest = 'Y:\home\nathan\DATA\S-probe\Awake\useful_ephys\analysis\behavior';
sessionID = '1046-191216';
behav_output = load([behav_dest '\' sessionID '_behavAnalysis.mat' ]);
frames_runTrigger_mat = behav_output.frames_runTrigger_mat;
frames_runoff_mat = behav_output.frames_runoff_mat;
speed = behav_output.speed;
speed = double(speed);
channel = '191216_Ch6_MM';
analysis_dest = 'Y:\home\nathan\DATA\S-probe\Awake\useful_ephys\analysis\';
PC = load([analysis_dest channel '.mat']);
nspk_100ms_cutted = PC.nspk_100ms_cutted;
tcut = PC.tcut;

% get rid of the trials that are later then ephys FR drop off time point
[inx1,inx2] = find(frames_runTrigger_mat>tcut*10); %inx2 tells you which column (i.e. trial) is later then tcut
cut_cols = unique(inx2);
frames_runTrigger_mat_cut = frames_runTrigger_mat;
frames_runTrigger_mat_cut(:,cut_cols)=[];

[inx1,inx2] = find(frames_runoff_mat>tcut*10); %inx2 tells you which column (i.e. trial) is later then tcut
cut_cols = unique(inx2);
frames_runoff_mat_cut = frames_runoff_mat;
frames_runoff_mat_cut(:,cut_cols)=[];

runTrigSpeed = speed(frames_runTrigger_mat_cut);
if size(runTrigSpeed,2)&&size(runTrigSpeed,1)>1
    ave_speed_runTrigger = mean(runTrigSpeed,2);
    ste_speed_runTrigger = std(runTrigSpeed,0,2)/sqrt(size(runTrigSpeed,2));
else
    ave_speed_runTrigger = runTrigSpeed';
    ste_speed_runTrigger = zeros(26,1);
end

runoffSpeed = speed(frames_runoff_mat_cut);
if size(runoffSpeed,2)&&size(runoffSpeed,1)>1
    ave_speed_runoff = mean(runoffSpeed,2);
    ste_speed_runoff = std(runoffSpeed,0,2)/sqrt(size(runoffSpeed,2));
else
    ave_speed_runoff = runoffSpeed';
    ste_speed_runoff = zeros(26,1);
end
save([behav_dest '\' sessionID '_behavAnalysis.mat' ],...
    'runTrigSpeed','runoffSpeed','ave_speed_runTrigger','ave_speed_runoff',...
    'ste_speed_runTrigger','ste_speed_runoff','frames_runTrigger_mat_cut',...
    'frames_runoff_mat_cut','-append');


SStrig = nspk_100ms_cutted(frames_runTrigger_mat_cut);
if size(SStrig,2)>1 && size(SStrig,1)>1
    ave_SStrig = mean(SStrig,2)*10;
    ste_SStrig = std(SStrig,0,2)*10/sqrt(size(SStrig,2));
else
    ave_SStrig = (SStrig*10)';
    ste_SStrig = zeros(26,1);
end

SSoff = nspk_100ms_cutted(frames_runoff_mat_cut);
if size(SSoff,2)>1 && size(SSoff,1) >1
    ave_SSoff = mean(SSoff,2)*10;
    ste_SSoff = std(SSoff,0,2)*10/sqrt(size(SSoff,2));
else
    ave_SSoff = (SSoff*10)';
    ste_SSoff = zeros(26,1);
end
save([analysis_dest channel '.mat'],'ave_SStrig','ave_SSoff','ste_SStrig','ste_SSoff','-append');


