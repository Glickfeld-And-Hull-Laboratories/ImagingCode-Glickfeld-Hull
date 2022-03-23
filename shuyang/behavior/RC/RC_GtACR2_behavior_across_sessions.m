%% peak lick rate on each trial
clear;
% naive day
%sessions = {'210909_img1098','210921_img1901'};
% CS1+GtACR2 after training day before imaging
%sessions = {'210920_img1098','210927_img1901'};

% CS1+GtACR2 after training
% sessions = {'210924_img1098','210929_img1901'};
% CS1 training without IO inhibition
sessions = {'211006_img1098','211008_img1901'}; 

run = {'000','000','000','000','000','000','000'};
analysis_out    = 'Z:\home\shuyang\2P_analysis\';
color_code = 'b';
% fig_name1 = 'across_session_lickpeaktime_CS1_naive';
% fig_name2 = 'across_session_lickpeakrate_CS1_naive';
% fig_name1 = 'across_session_lickpeaktime_CS1GtACR2_trained';
% fig_name2 = 'across_session_lickpeakrate_CS1GtACR2_trained';
fig_name1 = 'across_session_lickpeaktime_plus7dayCS1';
fig_name2 = 'across_session_lickpeakrate_plus7dayCS1';

max_tt_all = [];
max_hz_all = [];

for ii = 1:length(sessions)
    %load things
    lickAlign_file = dir([analysis_out sessions{ii} '\' sessions{ii} '_' run{ii} '*' 'cueAlignLick.mat']);
    assert(length(lickAlign_file)) = 1;
    lickAlign_output = load([analysis_out sessions{ii} '\' lickAlign_file.name]);
    lickAlign_thisSession = lickAlign_output.lickCueAlign; % frame*trial
    tt = lickAlign_output.tt;
    lickAlign_thisSession = lickAlign_output.lickCueAlign; % frame*trial
    lickAlign_thisSession_smooth = zeros(size(lickAlign_thisSession));
    for i = 1:size(lickAlign_thisSession,2)
        lickAlign_thisSession_smooth(:,i) = smooth(lickAlign_thisSession(:,i),10); %averaging every 10 frames
    end   
    tpluscue = 200; % 200ms threshold after cue, anything before then is a false alarm
    [tt_min tt_start] = min(abs(tt-tpluscue)); %find the index in tt when it's 200
    [max_val max_ind] = max(lickAlign_thisSession_smooth(tt_start:end,:),[],1);
    max_ind(find(max_val==0)) = size(lickAlign_thisSession,1)-tt_start; % if max licking rate is 0, set the index of max value to the end of the trial (-tt_start b/c you add it later)
    %max_ind(find(isnan(max_ind))) = size(lickAlign_thisSession,1)-tt_start;
    %max_val(find(isnan(max_ind))) = 0; % if no max value, set it to 0.
    max_tt = tt(max_ind+tt_start-1);
    frameRate = min(diff(tt));
    max_hz = max_val.*(1000/frameRate);
    
    max_tt_all = [max_tt_all max_tt];
    max_hz_all = [max_hz_all max_hz];
end

maxttcum_fig = figure;
maxttcum_fig.Units = 'inches';
maxttcum_fig.Position = [1 1 4 2.5];
cdfplot(max_tt_all);
vline(0,color_code);vline(700,'k');
ylabel('Fraction of max lick times');
xlim([-500 3500]);
ylim([0 1]);
%title('Pre-learning, CS1');
%title('Post CS1+IO inhibition before imaging day');
%title('Post CS1+IO inhibition training');
title('+ 8 DAYS CS1 training');
xlabel('Time from cue (ms)');
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'FontSize',7);
%xlim([-1.1 1.6]);ylim([0 3]);
box off;
path = 'Z:\home\shuyang\2P_analysis\forRO1\';
print(maxttcum_fig,[path,fig_name1],'-r600','-dpdf');
print(maxttcum_fig,[path,fig_name1],'-djpeg');
savefig([path fig_name1]);

maxhzcum_fig = figure;
maxhzcum_fig.Units = 'inches';
maxhzcum_fig.Position = [1 1 4 2.5];
cdfplot(max_hz_all);
%vline(0,color_code);vline(700,'k');
ylabel('Fraction of max lick rates');
%xlim([-500 3000]);
ylim([0 1]);
%title('Pre-learning, CS1');
%title('Post CS1+IO inhibition before imaging day');
%title('Post CS1+IO inhibition training');
title('+ 8 DAYS CS1 training');
xlabel('lick rates (Hz)');
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'FontSize',7);
%xlim([-1.1 1.6]);ylim([0 3]);
box off;
path = 'Z:\home\shuyang\2P_analysis\forRO1\';
print(maxhzcum_fig,[path,fig_name2],'-r600','-dpdf');
print(maxhzcum_fig,[path,fig_name2],'-djpeg');
savefig([path fig_name2]);


%% plot cumulative distribution function plot for RT
% assign document paths and experimental sessions
clear;
% naive day
%sessions = {'210909_img1098','210921_img1901'};
% CS1+GtACR2 after training day before imaging
%sessions = {'210920_img1098','210927_img1901'};

% CS1+GtACR2 after training
 sessions = {'210924_img1098','210929_img1901'};
% CS1 training without IO inhibition
%sessions = {'211006_img1098','211008_img1901'}; 

behav_analysis_out = 'Z:\home\shuyang\behavior_analysis\RC\';
color_code = 'b';
fig_name = 'across_session_RTCDF_CS1GtACR2_trained';
%fig_name = 'across_session_RTCDF_CS1_trained';
ntrials_total = 0;
RT = [];
timeaxis = 1:1:4700; % 0 is time of cue, 4700 is 4000ms after reward
for ii = 1:length(sessions)    
    %load(['Z:\home\shuyang\behavior_analysis\RC\RT_data\201113_img1078-1308.mat'])
    rawRT = load((fullfile([behav_analysis_out '\RT_data\' sessions{ii} '.mat'])));
    RT_this_session = rawRT.RT_this_session_raw;
    ntrials_total = ntrials_total + length(RT_this_session);
    %RT = [RT timeaxis(find(RT_this_session))]; 
    RT = [RT RT_this_session];
end
RT(RT<=0)=[];
% save(fullfile(['Z:\2P_analysis\across_sessions_RC\',data_file]),...
%     'ave_lick_targetAlign_acrossessions','ste_lick_targetAlign_acrossessions','tt',...
%     'ntrials_total','target_align_licks_all_session','-append');
% 

% plot the above section
close all;
RTcum_fig = figure; 
RTcum_fig.Units = 'inches';
RTcum_fig.Position = [1 1 4 2.5];
cdfplot(RT);
vline(0,color_code);vline(700,'k');
xlim([-500 5000]);
ylabel('Fraction of RTs');
ylim([0 1]);
%title('Pre-learning, CS1');
%title('Post CS1+IO inhibition before imaging day');
title('Post CS1+IO inhibition training');
%title('+ 8 DAYS CS1 training');
xlabel('Time from cue (ms)');
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'FontSize',7);
%xlim([-1.1 1.6]);ylim([0 3]);
box off;
path = 'Z:\home\shuyang\2P_analysis\forRO1\';
print(RTcum_fig,[path,fig_name],'-r600','-dpdf');
print(RTcum_fig,[path,fig_name],'-djpeg');
savefig([path fig_name]);


%% plot cumulative distribution function plot for licking data
% assign document paths and experimental sessions
clear;
% naive day
%sessions = {'210909_img1098','210921_img1901'};
% CS1+GtACR2 after training day before imaging
% sessions = {'210922_img1098','210928_img1901'};

% CS1+GtACR2 after training
% sessions = {'210924_img1098','210929_img1901'};
% CS1 training without IO inhibition
 sessions = {'211006_img1098','211008_img1901'}; 
run = {'000','000','000','000','000','000','000'}; 
analysis_out    = 'Z:\home\shuyang\2P_analysis\';
color_code = 'b';
% fig_name = 'across_session_targetalign_lickCDF_CS1GtACR2_daybeforeimg';
% fig_name = 'across_session_targetalign_lickCDF_CS1GtACR2';
 fig_name = 'across_session_targetalign_lickCDF_plus7daysCS1';

target_align_licks_all_session = [];
ntrials_total = 0;
lick_times = [];
for ii = 1:length(sessions)
    %load things
    lickAlign_file = dir([analysis_out sessions{ii} '\' sessions{ii} '_' run{ii} '*' 'cueAlignLick.mat']);
    assert(length(lickAlign_file)) = 1;
    lickAlign_output = load([analysis_out sessions{ii} '\' lickAlign_file.name]);
    lickAlign_thisSession = lickAlign_output.lickCueAlign; % frame*trial
    ntrials_total = ntrials_total + size(lickAlign_thisSession,2);
    
    for itrial = 1:size(lickAlign_thisSession,2)
        lick_times = [lick_times lickAlign_output.tt(find(lickAlign_thisSession(:,itrial)))];
    end
    
end

tt = lickAlign_output.tt;

% save(fullfile(['Z:\2P_analysis\across_sessions_RC\',data_file]),...
%     'ave_lick_targetAlign_acrossessions','ste_lick_targetAlign_acrossessions','tt',...
%     'ntrials_total','target_align_licks_all_session','-append');
% 

% plot the above section
lickcum_fig = figure; 
lickcum_fig.Units = 'inches';
lickcum_fig.Position = [1 1 4 2.5];
cdfplot(lick_times);
vline(0,color_code);vline(700,'k');
ylabel('Fraction of licks');
ylim([0 1]);
%title('Pre-learning, CS1');
%title('Post CS1+IO inhibition before imaging day');
%title('Post CS1+IO inhibition training');
title('+ 8 DAYS CS1 training');
xlabel('Time from cue (ms)');
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'FontSize',7);
%xlim([-1.1 1.6]);ylim([0 3]);
box off;
path = 'Z:\home\shuyang\2P_analysis\across_sessions_RC_GtACR2\';
print(lickcum_fig,[path,fig_name],'-r600','-dpdf');
print(lickcum_fig,[path,fig_name],'-djpeg');
savefig([path fig_name]);



%% licking aligned to cue, histogram
clear;
% naive day
%sessions = {'210909_img1098','210921_img1901'};
% CS1+GtACR2 after training
sessions = {'210922_img1098','210928_img1901'};

% CS1+GtACR2 after training
% sessions = {'210924_img1098','210929_img1901'};
% CS1 training without IO inhibition
%sessions = {'211006_img1098','211008_img1901'}; 


analysis_out    = 'Z:\home\shuyang\2P_analysis\';
behav_out  = 'Z:\home\shuyang\behavior_analysis\RC\hist_data_across_animals\';
% fig_name = 'across_session_firstLickAlign_spike_naive';
% data_file = 'across_sessions_RC_data_naive';
fig_name = 'across_session_firstLickAlign_spike_CS1_GtACR2dayBeforeImg';
%data_file = 'across_sessions_RC_data_CS1_GtACR2dayBeforeImg.mat';
% fig_name = 'across_session_firstLickAlign_spike_CS1_GtACR2trained';
% data_file = 'across_sessions_RC_data_CS1_GtACR2trained.mat';
% fig_name = 'across_session_firstLickAlign_spike_plus7days_CS1trained';
% data_file = 'across_sessions_RC_data_plus7days_CS1trained.mat';

% lick align average
licks_hist_all_session = [];
for ii = 1:length(sessions) 
    %load things  
    lickAlign_hist = load((fullfile([behav_out sessions{ii} 'rew_hist.mat'])));
    lickhist_thisSession = lickAlign_hist.full_trial_licks_rewarded_bin; %1*241 time bins (100ms per bin)
    licks_hist_all_session = cat(1,licks_hist_all_session,lickhist_thisSession);
end

licks_hist_all_avg = mean(licks_hist_all_session,1);
licks_hist_all_sem = std(licks_hist_all_session,0,1)/sqrt(size(licks_hist_all_session,1));
% save(fullfile(['Z:\home\shuyang\2P_analysis\across_sessions_RC_GtACR2\',data_file]),...
%     'licks_hist_all_avg','licks_hist_all_sem','-append');


num_animals = size(licks_hist_all_session, 1);
x_axis = (([1:length(licks_hist_all_avg)]-181)*100 + 50)/1000;

lickhist_fig = figure; 
lickhist_fig.Units = 'centimeters';
lickhist_fig.Position = [1 1 18 8];
bar(x_axis, licks_hist_all_avg,'FaceColor',[0.6 0.6 0.6],'EdgeColor',[0.6 0.6 0.6]); hold on;
errorbar(x_axis, licks_hist_all_avg, licks_hist_all_sem,'Color',[0.5 0.5 0.5], 'LineStyle', 'none', 'CapSize', 1);
%title('Pre-learning, CS1');
title('Post CS1+IO inhibition before imaging day');
%title('Post CS1+IO inhibition training');
%title('+ 8 DAYS CS1 training');
xlabel('time relative to cue (s)');
ylabel('lick rate (Hz)');
ylim([0 10]);
vline(0.7, 'k');
vline(0, 'r')
xlim([-1 4]);
text(3.4, 9.4,['n = ', num2str(num_animals)],'FontSize', 20,'FontName','Arial');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',20,'FontName','Arial');
box off;
path = 'Z:\home\shuyang\2P_analysis\across_sessions_RC_GtACR2\';
print(lickhist_fig,[path,fig_name],'-r600','-dpdf');
print(lickhist_fig,[path,fig_name],'-djpeg');
savefig([path fig_name]);

