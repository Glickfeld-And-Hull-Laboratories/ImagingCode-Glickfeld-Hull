%% SECTION - assign pathnames and datasets to be analyzed/written. 
clear;
sessions = {'190429_img1021','190430_img1023','190507_img1024','190603_img1025'};
days = {'1021-190429_1','1023-190430_1','1024-190507_1','1025-190603_1'};

%sessions = {'190430_img1023'};
%days = {'1023-190430_1'};
%there might be more than 1 sessions on a single subject on the same day
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder

% behavior analysis results 
color_code = {'r','g','m','y','b'};

%% 1. when does the first peak happen?
runoff_bi_all = [];
for i = 1:length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{i}, '\'];
    dfOvF_output = load([image_analysis_dest sessions{i} '_dfOvF.mat']);
    runoff_bi = dfOvF_output.runoff_bi; %frames*trials
    runoff_bi_all = cat(2,runoff_bi_all,runoff_bi);
end
fig_dest = 'Z:\Analysis\2P_MovingDots_Analysis\across_sessions\';

x = (0.1:0.1:2.5);
y = [0 size(runoff_bi_all,2)];
figure; imagesc(x,y,runoff_bi_all');
ylabel('trials');
xlabel('time(s)');
title('across sessions');
savefig([fig_dest 'across_sessions_runoff_1peakTime']);

%% when does all big peaks happen?
big_pks_all = [];
for i = 1:length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{i}, '\'];
    dfOvF_output = load([image_analysis_dest sessions{i} '_dfOvF.mat']);
    big_pks = dfOvF_output.big_pks; %frames*trials
    big_pks_all = cat(2,big_pks_all,big_pks);
end
fig_dest = 'Z:\Analysis\2P_MovingDots_Analysis\across_sessions\';

x = (0.1:0.1:2.5);
y = [0 size(big_pks_all,2)];
figure; imagesc(x,y,big_pks_all');
ylabel('trials');
xlabel('time(s)');
title('across sessions, time of all big peaks');
savefig([fig_dest 'across_sessions_runoff_allpeaksTime']);

%% is trial length and average speed correlated with when the first peak happens/how many peaks are there?
peak_loc_all = [];
sumpks_all = [];
trial_lens_all = [];
avespd_runofftrials_all = [];

for i = 1:length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{i}, '\'];
    dfOvF_output = load([image_analysis_dest sessions{i} '_dfOvF.mat']);
    peak_loc = dfOvF_output.peak_loc;
    peak_loc_all = cat(2,peak_loc_all,peak_loc);
    sumpks = dfOvF_output.sumpks;
    sumpks_all = cat(2,sumpks_all,sumpks);
    avespd_runofftrials = dfOvF_output.avespd_runofftrials;
    avespd_runofftrials_all = cat(2,avespd_runofftrials_all,avespd_runofftrials);
    trial_lens = dfOvF_output.trial_lens;
    trial_lens_all = cat(2,trial_lens_all,trial_lens);
end

fig_dest = 'Z:\Analysis\2P_MovingDots_Analysis\across_sessions\';
%sumpks and ave speed/trial lengths

figure; scatter(avespd_runofftrials_all*2*3.1415926*7.5/128,sumpks_all,'filled','MarkerFacecolor',[0.9373 0.396 0.2824]);
ylim([0 max(sumpks_all)]); xlim([0 max(avespd_runofftrials_all*2*3.1415926*7.5/128)]);
ylabel('total number of df/f peaks');
xlabel('average speed (cm/s)');
title('across sessions');
savefig([fig_dest 'across_sessions_runoff_sum_pks_Vs_speed']);

figure; scatter(trial_lens_all./30,sumpks_all,'filled','MarkerFacecolor',[0.9373 0.396 0.2824]);
ylim([0 max(sumpks_all)]);xlim([0 max(trial_lens_all./30)]);
ylabel('total number of df/f peaks');
xlabel('running trial duration (s)');
title('across sessions');
savefig([fig_dest 'across_sessions_sum_pks_Vs_duration']);

figure; scatter(avespd_runofftrials_all*2*3.1415926*7.5/128,peak_loc_all./30,'filled','MarkerFacecolor',[0.9373 0.396 0.2824]);
ylim([0 max(peak_loc_all/30)]); xlim([0 max(avespd_runofftrials_all*2*3.1415926*7.5/128)]);
ylabel('time of first peak after running offset(s)');
xlabel('average speed (cm/s)');
title('across sessions');
savefig([fig_dest 'across_sessions_1pk_time_Vs_speed']);

figure; scatter(trial_lens_all./30,peak_loc_all./30,'filled','MarkerFacecolor',[0.9373 0.396 0.2824]);
ylim([0 max(peak_loc_all/30)]);xlim([0 max(trial_lens_all./30)]);
ylabel('time of first peak after running offset(s)');
xlabel('running trial duration (s)');
title('across sessions');
axis square;
savefig([fig_dest 'across_sessions_1pk_time_Vs_duration']);



