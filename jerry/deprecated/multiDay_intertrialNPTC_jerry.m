clear all; clear global; close all;
clc

ds = 'DART_expt_info'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories
eval(ds);
nd=2; 

day_id = [42]; % concat days, enter the post-DART
mouse = expt(day_id).mouse;

fnout = fullfile(rc.achAnalysis,'PV_YM90K',mouse);
% if expt(day_id).multiday_timesincedrug_hours>0
%     dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
% else
%     dart_str = 'control';
% end
% 
% fn_multi = fullfile(rc.achAnalysis,'PV_YM90K',mouse,['multiday_' dart_str]);

matched_day = expt(day_id).multiday_matchdays;
matched_day_run = expt(matched_day).contrastxori_runs;
fn_single_d1 = fullfile(rc.achAnalysis,'PV_YM90K',mouse,expt(matched_day).date,matched_day_run);
fn_single_d2 = fullfile(rc.achAnalysis,'PV_YM90K',mouse,expt(day_id).date,expt(day_id).contrastxori_runs);
fn_across_tp = fullfile(rc.achAnalysis,'PV_YM90K',mouse,'across_tps');
mkdir(fn_across_tp);

nOn = 30;
nOff = 60;
nTot = nOn+nOff;
ntrials = 960;
%% read matched cell time courses
% load(fullfile(fn_multi,'timecourses.mat'));
% d1TCs = cellTCs_match{1};
% d2TCs = cellTCs_match{2};
% twtc_d1 = reshape(d1TCs,[],ntrials,size(d1TCs,2));
% twtc_d2 = reshape(d2TCs,[],ntrials,size(d2TCs,2));
% 
% twtc_off = mean(squeeze(mean(twtc_d1(nOff/2:nOff,:,:),2)),2);
% 
% all_frames1  = mean(d1TCs,2);
% all_frames2  = mean(d2TCs,2);
% 
% plot(all_frames1);
% figure
% plot(all_frames2);

%% read neuropil timecourse from each day, NOT MATCHED
% day1 
load(fullfile(fn_single_d1{1},'TCs.mat'));

avg_np_tc1 = mean(np_tc,2);
avg_np_tc1_trial = reshape(avg_np_tc1,nTot,[]);
avg_np_tc1_f0 = avg_np_tc1_trial(nOff/2+1:nOff,:);
avg_np_tc1_off = reshape(avg_np_tc1_f0,[],1);

mean_sub_tc1 = avg_np_tc1_off - mean(avg_np_tc1_off);

figure
plot(avg_np_tc1_off);
ylim([400 700]);
saveas(gcf,fullfile(fn_across_tp,['OffFrames_' expt(matched_day).multiday_timesincedrug_hours '_raw.pdf']));

figure
plot(mean_sub_tc1);
ylim([-100 200]);
saveas(gcf,fullfile(fn_across_tp,['OffFrames_' expt(matched_day).multiday_timesincedrug_hours '_ms.pdf']));

clear np_tc npSub_tc

% day2
load(fullfile(fn_single_d2{1},'TCs.mat'));

avg_np_tc2 = mean(np_tc,2);
avg_np_tc2_trial = reshape(avg_np_tc2,nTot,[]);
avg_np_tc2_f0 = avg_np_tc2_trial(nOff/2+1:nOff,:);
avg_np_tc2_off = reshape(avg_np_tc2_f0,[],1);

mean_sub_tc2 = avg_np_tc2_off - mean(avg_np_tc2_off);

figure
plot(avg_np_tc2_off);
ylim([400 700]);
saveas(gcf,fullfile(fn_across_tp,['OffFrames_' expt(day_id).multiday_timesincedrug_hours '_raw.pdf']));

figure
plot(mean_sub_tc2);
ylim([-100 200]);
saveas(gcf,fullfile(fn_across_tp,['OffFrames_' expt(day_id).multiday_timesincedrug_hours '_ms.pdf']));

clear np_tc npSub_tc

