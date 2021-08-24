% Reward Conditioning Task: average neural align, lickalign, and piezoalign across sessions 
%% align neural data with cue across sessions 
% assign document paths and experimental sessions
clear;
% CS1 naive
% sessions = {'201113_img1078','201114_img1079','210112_img1083',...
% '210405_img1087','210416_img1088','210412_img1090'};
% CS1 trained
% sessions = {'201125_img1078','201125_img1079','210119_img1083','210416_img1087',...
% '210428_img1088','210427_img1090'};
% CS1+CS2 test day after training with CS1+CS2
% sessions = {'201208_img1078','201209_img1079','210202_img1083','210428_img1087',...
%     '210507_img1088','210506_img1090'};
% CS1 test day after training with CS1+CS2 %!!!!!!!!!!!!!!!!!201209_img1079 the run is 001
sessions = {'201211_img1078','201209_img1079','210203_img1083','210429_img1087',...
 '210508_img1088','210507_img1090'};
% CS2 test day after training with CS1+CS2
% sessions = {'201210_img1078','201210_img1079','210205_img1083','210430_img1087',...
% '210509_img1088','210511_img1090'};

% control group- 1901 and 1092
% naive day
%  sessions = {'210616_img1091','210616_img1092'};
% trained day 9
% sessions = {'210628_img1091','210628_img1092'};

run = {'000','001','000','000','000','000'}; %!!!!!!!!!!!!!!!!! 201209_img1079 the run is 001 for CS1 test day after training with CS1+2
analysis_out    = 'Z:\2P_analysis\';
color_code = 'b';
% fig_name = 'across_session_targetalign_spike_CS1_naive_control';
% data_file = 'across_sessions_RC_data_CS1_naive_control.mat';
% fig_name = 'across_session_targetalign_spike_CS1_naive';
% data_file = 'across_sessions_RC_data_CS1_naive.mat';
% fig_name = 'across_session_targetalign_spike_CS1_trained';
% data_file = 'across_sessions_RC_data_CS1_trained.mat';
% fig_name = 'across_session_targetalign_spike_CS1_trained_control';
% data_file = 'across_sessions_RC_data_CS1_trained_control.mat';
% fig_name = 'across_session_targetalign_spike_CS1CS2_testCS1CS2';
% data_file = 'across_sessions_RC_data_CS1CS2_testCS1CS2.mat';
fig_name = 'across_session_targetalign_spike_CS1CS2_testCS1';
data_file = 'across_sessions_RC_data_CS1CS2_testCS1.mat';
% fig_name = 'across_session_targetalign_spike_CS1CS2_testCS2';
% data_file = 'across_sessions_RC_data_CS1CS2_testCS2.mat';


% target align average
target_align_events_all_session = [];
nPCs_total = 0;
ntrials_total = 0;
nPCs = [];
for ii = 1:length(sessions) 
    %load things  
    targetAlign_output = load((fullfile([analysis_out sessions{ii} '\' sessions{ii} '_' run{ii} '_targetAlign.mat'])));
    targetAlign_events_thisSession = targetAlign_output.targetAlign_events; % frame*cell*trial
    nPCs_total = nPCs_total + size(targetAlign_events_thisSession,2);
    nPCs = cat(2,nPCs,size(targetAlign_events_thisSession,2));
    ntrials_total = ntrials_total + size(targetAlign_events_thisSession,3);
    
    %for spike data, average first. average across cells and trials
    mean_targetAlign_events = squeeze(nanmean(targetAlign_events_thisSession,3));    
    target_align_events_all_session = cat(2,target_align_events_all_session,mean_targetAlign_events);
end

tt = targetAlign_output.tt;

ave_spike_targetAlign_acrossessions = mean(target_align_events_all_session,2)*30;
ste_spike_targetAlign_acrossessions = std(target_align_events_all_session,0,2)*30/sqrt(size(target_align_events_all_session,2));

save(fullfile(['Z:\2P_analysis\across_sessions_RC\',data_file]),...
    'ave_spike_targetAlign_acrossessions','target_align_events_all_session','tt',...
    'nPCs_total','nPCs','ntrials_total','target_align_events_all_session','-append');


% plot the above section
targetAlign_fig = figure; 
targetAlign_fig.Units = 'inches';
targetAlign_fig.Position = [1 1 4 2.5];
%targetAlign_fig.PaperOrientation = 'landscape';
shadedErrorBar(tt, ave_spike_targetAlign_acrossessions, ste_spike_targetAlign_acrossessions);
ylim([0 3]);vline(0,color_code);vline(700,'k');
ylabel('Spike rate (Hz)');
%title('CS1 naive day, across sessions');
%title('CS1 trained day, across sessions');
%title('CS1+CS2 trained test w/ CS1+CS2, across sessions');
title('CS1+CS2 trained test w/ CS1, across sessions');
%title('CS1+CS2 trained test w/ CS2, across sessions');
xlabel('Time from cue (ms)');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',10);
%xlim([-1.1 1.6]);ylim([0 3]);
box off;
path = 'Z:\2P_analysis\across_sessions_RC\';
print(targetAlign_fig,[path,fig_name],'-r600','-dpdf');
print(targetAlign_fig,[path,fig_name],'-djpeg');
savefig([path fig_name]);


%% align neural data with cue onset on learned days, late lick trials vs. early lick trials
% results for this section come from RC_lickAlign_SJ
clear
% CS1 naive
% sessions = {'201113_img1078','201114_img1079','210112_img1083',...
% '210405_img1087','210412_img1090','210416_img1088'};
% CS1 trained
% sessions = {'201125_img1078','201125_img1079','210119_img1083','210416_img1087',...
% '210428_img1088','210427_img1090'};

% control group- 1901 and 1092
% naive day
% sessions = {'210616_img1091','210616_img1092'};
% trained day 9
sessions = {'210628_img1091','210628_img1092'};
run = {'000','000','000','000','000','000'};%,'000'};
analysis_out    = 'Z:\2P_analysis\';
color_code = 'b';
% fig_name = 'across_session_CS1_naive_byLickTime';
% data_file = 'across_sessions_RC_data_CS1_naive.mat';
% fig_name = 'across_session_CS1_naive_control_byLickTime';
% data_file = 'across_sessions_RC_data_CS1_naive_control.mat';
% fig_name = 'across_session_CS1_trained_byLickTime';
% data_file = 'across_sessions_RC_data_CS1_trained.mat';
fig_name = 'across_session_CS1_trained_control_byLickTime';
data_file = 'across_sessions_RC_data_CS1_trained_control.mat';
targetAlign_events_earlybst = [];
targetAlign_events_latebst = [];
nPCs_total = 0;
earlyBst_sessions = [];
lateBst_sessions = [];
for ii = 1:length(sessions) 
    %load things  
    targetAlign_output = load((fullfile([analysis_out sessions{ii} '\' sessions{ii} '_' run{ii} '_targetAlign.mat'])));
    targetAlign_events_thisSession = targetAlign_output.targetAlign_events; % frame*cell*trial
    prewin_frames = targetAlign_output.prewin_frames;
    frameRateHz = targetAlign_output.frameRateHz;
    nPCs_total = nPCs_total + size(targetAlign_events_thisSession,2);
    lickAlign_output = load((fullfile([analysis_out sessions{ii} '\' sessions{ii} '_' run{ii} '_cueAlignLick.mat'])));
    ind_early_bst_thisSession = lickAlign_output.ind_early_bst; % index of the early and late burst trials
    ind_late_bst_thisSession = lickAlign_output.ind_late_bst;
    %need to average across trials before cat, otherwise # of cells and
    %trials are both different across sessions, can't cat
    targetAlign_events_earlybst_thisSession = targetAlign_events_thisSession(:,:,ind_early_bst_thisSession);
    targetAlign_events_earlybst_cells = nanmean(targetAlign_events_earlybst_thisSession,3);
    targetAlign_events_earlybst = cat(2,targetAlign_events_earlybst,targetAlign_events_earlybst_cells);
    
    targetAlign_events_latebst_thisSession = targetAlign_events_thisSession(:,:,ind_late_bst_thisSession);
    targetAlign_events_latebst_cells = nanmean(targetAlign_events_latebst_thisSession,3);
    targetAlign_events_latebst = cat(2,targetAlign_events_latebst,targetAlign_events_latebst_cells);
    
    lickBurstStart_thisSession = lickAlign_output.lickBurstStart; % when does burst start in each trial
    earlyBst_thisSession = (lickBurstStart_thisSession(:,ind_early_bst_thisSession)-prewin_frames)*1000/frameRateHz;%how many ms after cue do lick busts happen in these early trials of this session
    earlyBst_sessions = cat(2,earlyBst_sessions,earlyBst_thisSession);
    lateBst_thisSession = (lickBurstStart_thisSession(:,ind_late_bst_thisSession)-prewin_frames)*1000/frameRateHz;%how many ms after cue do lick busts happen in these early trials of this session
    lateBst_sessions = cat(2,lateBst_sessions,lateBst_thisSession);
    
%     earlyBst_sessions(ii) = (mean(lickBurstStart_thisSession(:,ind_early_bst_thisSession))-prewin_frames)*1000/frameRateHz;%on average how many ms lick busts happen in these early trials of this session
%     lateBst_sessions(ii) = (mean(lickBurstStart_thisSession(:,ind_late_bst_thisSession))-prewin_frames)*1000/frameRateHz;%on average how many ms lick busts happen in these late trials of this session
end
tt = lickAlign_output.tt;

neural_by_lick_time = figure; %plot neural data of trials of early vs. late bursts
neural_by_lick_time.Units = 'inches';
neural_by_lick_time.Position = [1 1 4 2.5];
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',10);
shadedErrorBar(tt,nanmean(targetAlign_events_earlybst,2).*(1000./frameRateHz), (nanstd(targetAlign_events_earlybst,0,2)).*(1000./frameRateHz)./sqrt(nPCs_total), 'k');
hold on;
shadedErrorBar(tt,nanmean(targetAlign_events_latebst,2).*(1000./frameRateHz), (nanstd(targetAlign_events_latebst,0,2)).*(1000./frameRateHz)./sqrt(nPCs_total), 'b');
errorbar(mean(earlyBst_sessions),2,std(earlyBst_sessions),'horizontal','-s','LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k','Color','k');
errorbar(mean(lateBst_sessions),2,std(lateBst_sessions),'horizontal','-s','LineStyle','none','MarkerFaceColor','b','MarkerEdgeColor','b','Color','b');
xlabel('Time from cue (ms)');
ylabel('Spike rate (Hz)');
ylim([0 3]);
vline(0,color_code);vline(700,'k');
%title(['CS1 naive, - lick bursts: early = ' num2str(mean(earlyBst_sessions)) 'ms; late = ' num2str(mean(lateBst_sessions)) 'ms']);
title(['CS1 trained, - lick bursts: early = ' num2str(mean(earlyBst_sessions)) 'ms; late = ' num2str(mean(lateBst_sessions)) 'ms']);
hold off;
box off;
path = 'Z:\2P_analysis\across_sessions_RC\';
print(neural_by_lick_time,[path,fig_name],'-r600','-dpdf');
print(neural_by_lick_time,[path,fig_name],'-djpeg');
savefig([path fig_name]);

save(fullfile(['Z:\2P_analysis\across_sessions_RC\',data_file]),...
    'targetAlign_events_earlybst','targetAlign_events_latebst',...
    'earlyBst_sessions','lateBst_sessions',...
    'nPCs_total','-append');


%% align neural data with lickonset across sessions 
% assign document paths and experimental sessions
clear;
% CS1 naive
% sessions = {'201113_img1078','201114_img1079','210112_img1083','210405_img1087','210416_img1088','210412_img1090'};
% CS1 trained
% sessions = {'201125_img1078','201125_img1079','210119_img1083','210416_img1087',...
% '210427_img1090','210428_img1088'};

% control group- 1901 and 1092
% naive day
% sessions = {'210616_img1091','210616_img1092'};
% trained day 9
 sessions = {'210628_img1091','210628_img1092'};

%CS2 test day after training with CS1+CS2
% sessions = {'201210_img1078','201210_img1079','210205_img1083','210430_img1087',...
%  '210509_img1088','210511_img1090'};
run = {'000','000','000','000','000','000'};
analysis_out    = 'Z:\2P_analysis\';
color_code = 'k';
%data_file = 'across_sessions_RC_data_CS1_naive_lickAlign.mat';
%data_file = 'across_sessions_RC_data_CS1_naive_control_lickAlign.mat';
%data_file = 'across_sessions_RC_data_CS1_trained_lickAlign.mat';
data_file = 'across_sessions_RC_data_CS1_trained_control_lickAlign.mat';
%data_file = 'across_sessions_RC_data_CS1CS2_testCS2_lickAlign.mat';

% target align average
firstlick_align_events_all_session = [];% frames*total number of sessions 
lickBurstAlign_events_all_session = [];
firstLickAlign_all_session = [];
lickBurstAlign_postRew_all_session = [];

nPCs_total = 0;
ntrials_total = 0;
nPCs = [];
for ii = 1:length(sessions) 
    %load things  
    lickAlign_output = load((fullfile([analysis_out sessions{ii} '\' sessions{ii} '_' run{ii} '_cueAlignLick.mat'])));
    
    firstLickAlign_events_thisSession = lickAlign_output.firstPostRew_lickAlignEvents_1500_3000ms_scale; % frame*cell*trial
    firstLickAlign_thisSession = lickAlign_output.firstPostRew_lickAlign_1500_3000ms_scale; % frame*trial
    
    lickburstAlign_postRew_events_thisSession = lickAlign_output.postRew_lickAlignEvents_1500_3000ms_scale;
    lickbursAlign_postRew_thisSession = lickAlign_output.postRew_lickAlign_1500_3000ms_scale;
    
    nPCs_total = nPCs_total + size(firstLickAlign_events_thisSession,2);
    nPCs = cat(2,nPCs,size(firstLickAlign_events_thisSession,2));
    ntrials_total = ntrials_total + size(firstLickAlign_events_thisSession,3);
    
    %for spike data, average first. average across cells and trials
    mean_firstLickAlign_events = squeeze(nanmean(firstLickAlign_events_thisSession,3));
    firstlick_align_events_all_session = cat(2,firstlick_align_events_all_session,mean_firstLickAlign_events);
    
    mean_lickBurstAlign_events = squeeze(nanmean(lickburstAlign_postRew_events_thisSession,3));
    lickBurstAlign_events_all_session = cat(2,lickBurstAlign_events_all_session,mean_lickBurstAlign_events);
    
    % average licking data
    firstLickAlign_all_session = cat(2,firstLickAlign_all_session,firstLickAlign_thisSession);
    lickBurstAlign_postRew_all_session = cat(2,lickBurstAlign_postRew_all_session,lickbursAlign_postRew_thisSession);
    
    
end

tt = lickAlign_output.tt;

ave_spike_firstLickAlign_acrossessions = nanmean(firstlick_align_events_all_session,2)*30;
ste_spike_firstLickAlign_acrossessions = nanstd(firstlick_align_events_all_session,0,2)*30/sqrt(size(firstlick_align_events_all_session,2));

ave_spike_lickBurstAlign_acrossessions = nanmean(lickBurstAlign_events_all_session,2)*30;
ste_spike_lickBurstAlign_acrossessions = nanstd(lickBurstAlign_events_all_session,0,2)*30/sqrt(size(lickBurstAlign_events_all_session,2));

ave_firstLickAlign_acrossessions = nanmean(firstLickAlign_all_session,2)*30;
ste_firstLickAlign_acrossessions = nanstd(firstLickAlign_all_session,0,2)*30/sqrt(size(firstLickAlign_all_session,2));

ave_lickBurstAlign_acrossessions = nanmean(lickBurstAlign_postRew_all_session,2)*30;
ste_lickBurstAlign_acrossessions = nanstd(lickBurstAlign_postRew_all_session,0,2)*30/sqrt(size(lickBurstAlign_postRew_all_session,2));


save(fullfile(['Z:\2P_analysis\across_sessions_RC\',data_file]),...
    'ave_spike_firstLickAlign_acrossessions','ste_spike_firstLickAlign_acrossessions',...
    'ave_spike_lickBurstAlign_acrossessions','ste_spike_lickBurstAlign_acrossessions',...
    'ave_firstLickAlign_acrossessions','ste_firstLickAlign_acrossessions',...
    'ave_lickBurstAlign_acrossessions','ste_lickBurstAlign_acrossessions',...
    'firstlick_align_events_all_session','lickBurstAlign_events_all_session',...
    'firstLickAlign_all_session','lickBurstAlign_postRew_all_session',...
    'tt','nPCs_total','ntrials_total');


% plot the above section
firstLickAlign_fig = figure; 
firstLickAlign_fig.Units = 'inches';
firstLickAlign_fig.Position = [1 1 3.5 3];
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',10);
subplot(2,1,1);
shadedErrorBar(tt, ave_spike_firstLickAlign_acrossessions, ste_spike_firstLickAlign_acrossessions);
ylim([0 3.5]); vline(0,color_code);ylabel('Spike rate (Hz)');
%title('CS1 naive,  align to first lick post reward, across sessions');
title('CS1 trained,  align to first lick post reward, across sessions');
%title('CS1+CS2 trained, test w/CS2,  align to first lick post reward, across sessions');
box off;
subplot(2,1,2);
shadedErrorBar(tt, ave_firstLickAlign_acrossessions, ste_firstLickAlign_acrossessions);
ylim([0 35]); vline(0,color_code);ylabel('Lick rate (Hz)');
xlabel('Time from first lick after reward (ms)');
box off;
%xlim([-1.1 1.6]);ylim([0 3]);
path = 'Z:\2P_analysis\across_sessions_RC\';
%fig_name = 'across_session_firstLickAlign_CS1_naive';
%fig_name = 'across_session_firstLickAlign_CS1_naive_control';
%fig_name = 'across_session_firstLickAlign_CS1_trained';
fig_name = 'across_session_firstLickAlign_CS1_trained_control';
%fig_name = 'across_session_firstLickAlign_CS1CS2_testCS2';
print(firstLickAlign_fig,[path,fig_name],'-r600','-dpdf');
print(firstLickAlign_fig,[path,fig_name],'-djpeg');
savefig([path fig_name]);

lickBurstAlign_fig = figure; 
lickBurstAlign_fig.Units = 'inches';
lickBurstAlign_fig.Position = [1 1 3.5 3];
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',10);
subplot(2,1,1);
shadedErrorBar(tt, ave_spike_lickBurstAlign_acrossessions, ste_spike_lickBurstAlign_acrossessions);
ylim([0 3.5]); vline(0,color_code);ylabel('Spike rate (Hz)');
%title('CS1 naive,  align to lick burst post reward, across sessions');
title('CS1 trained,  align to lick burst post reward, across sessions');
%title('CS1+CS2 trained, test w/CS2, align to lick burst post reward, across sessions');
box off;
subplot(2,1,2);
shadedErrorBar(tt, ave_lickBurstAlign_acrossessions, ste_lickBurstAlign_acrossessions);
ylim([0 35]); vline(0,color_code);ylabel('Lick rate (Hz)'); 
xlabel('Time from lick burst after reward (ms)');
box off;
path = 'Z:\2P_analysis\across_sessions_RC\';
%fig_name = 'across_session_lickBurstAlign_CS1_naive';
%fig_name = 'across_session_lickBurstAlign_CS1_naive_control';
%fig_name = 'across_session_lickBurstAlign_CS1_trained';
fig_name = 'across_session_lickBurstAlign_CS1_trained_control';
%fig_name = 'across_session_lickBurstAlign_CS1CS2_testCS2';
print(lickBurstAlign_fig,[path,fig_name],'-r600','-dpdf');
print(lickBurstAlign_fig,[path,fig_name],'-djpeg');
savefig([path fig_name]);


%% align neural data with cue across sessions, interleaved day 
% assign document paths and experimental sessions
clear;
% interleaved day
sessions = {'201217_img1078','201211_img1079','210209_img1083','210503_img1087',...
'210511_img1088','210509_img1090'};
%no response to CS2
run = {'000','000','000','000','000','000'};
analysis_out    = 'Z:\2P_analysis\';
fig_name = 'across_session_targetalign_spike_CS1CS2_interleave';
data_file = 'across_sessions_RC_data_CS1CS2_interleave.mat';

% target align average
nPCs = [];
target_align_events_audi_all_session = [];
target_align_events_vis_all_session = [];
target_align_events_av_all_session = [];
nPCs_total = 0;
ntrials_total = 0;
for ii = 1:length(sessions) 
    %load things  
    targetAlign_output = load((fullfile([analysis_out sessions{ii} '\' sessions{ii} '_' run{ii} '_targetAlign.mat'])));
    targetAlign_events_thisSession = targetAlign_output.targetAlign_events; % frame*cell*trial
    nPCs_total = nPCs_total + size(targetAlign_events_thisSession,2);
    ntrials_total = ntrials_total + size(targetAlign_events_thisSession,3);
    nPCs = cat(2,nPCs,size(targetAlign_events_thisSession,2));
    trailTypes = load (['Z:\behavior_analysis\RC\' sessions{ii} '_IL_trial_inx.mat']);
    auditrials = trailTypes.auditrials;
    vistrials = trailTypes.vistrials;
    avtrials = trailTypes.avtrials;
    %for spike data, average first. average across cells and trials
    targetAlign_events_audi_thisSession  = targetAlign_events_thisSession(:,:,auditrials);
    targetAlign_events_vis_thisSession  = targetAlign_events_thisSession(:,:,vistrials);
    targetAlign_events_av_thisSession  = targetAlign_events_thisSession(:,:,avtrials);
    
    mean_targetAlign_audi_events = squeeze(nanmean(targetAlign_events_audi_thisSession,3));
    mean_targetAlign_vis_events = squeeze(nanmean(targetAlign_events_vis_thisSession,3));
    mean_targetAlign_av_events = squeeze(nanmean(targetAlign_events_av_thisSession,3));
    
    target_align_events_audi_all_session = cat(2,target_align_events_audi_all_session,mean_targetAlign_audi_events);
    target_align_events_vis_all_session = cat(2,target_align_events_vis_all_session,mean_targetAlign_vis_events);
    target_align_events_av_all_session = cat(2,target_align_events_av_all_session,mean_targetAlign_av_events);
end

tt = targetAlign_output.tt;

ave_spike_targetAlign_audi_acrossessions = mean(target_align_events_audi_all_session,2)*30;
ste_spike_targetAlign_audi_acrossessions = std(target_align_events_audi_all_session,0,2)*30/sqrt(size(target_align_events_audi_all_session,2));
ave_spike_targetAlign_vis_acrossessions = mean(target_align_events_vis_all_session,2)*30;
ste_spike_targetAlign_vis_acrossessions = std(target_align_events_vis_all_session,0,2)*30/sqrt(size(target_align_events_vis_all_session,2));
ave_spike_targetAlign_av_acrossessions = mean(target_align_events_av_all_session,2)*30;
ste_spike_targetAlign_av_acrossessions = std(target_align_events_av_all_session,0,2)*30/sqrt(size(target_align_events_av_all_session,2));

save(fullfile(['Z:\2P_analysis\across_sessions_RC\',data_file]),...
    'target_align_events_audi_all_session','target_align_events_vis_all_session',...
    'target_align_events_av_all_session','tt',...
    'nPCs_total','ntrials_total');


targetAlign_interleave_fig = figure; 
targetAlign_interleave_fig.Units = 'inches';
targetAlign_interleave_fig.Position = [1 1 3.5 3.5];
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',10);
subplot(3,1,1);
shadedErrorBar(tt, ave_spike_targetAlign_audi_acrossessions, ste_spike_targetAlign_audi_acrossessions);
ylim([0 2.5]);vline(0,'b');vline(700,'k');ylabel('Spike rate (Hz)');
title('CS1+CS2 trained test w/ interleaved, across sessions');
box off;
subplot(3,1,2);
shadedErrorBar(tt, ave_spike_targetAlign_vis_acrossessions, ste_spike_targetAlign_vis_acrossessions);
ylim([0 2.5]);vline(0,'g');vline(700,'k');ylabel('Spike rate (Hz)');
box off;
subplot(3,1,3);
shadedErrorBar(tt, ave_spike_targetAlign_av_acrossessions, ste_spike_targetAlign_av_acrossessions);
ylim([0 2.5]);vline(0,'r');vline(700,'k');ylabel('Spike rate (Hz)');
box off;
xlabel('Time from cue (ms)');
path = 'Z:\2P_analysis\across_sessions_RC\';
print(targetAlign_interleave_fig,[path,fig_name],'-r600','-dpdf');
print(targetAlign_interleave_fig,[path,fig_name],'-djpeg');
savefig([path fig_name]);


%% align neural data with lickonset across sessions - interleaved session
% assign document paths and experimental sessions
clear;
sessions = {'201211_img1079','201217_img1078','210209_img1083','210503_img1087',...
'210511_img1088','210509_img1090'};
run = {'000','000','000','000','000','000'};
analysis_out    = 'Z:\2P_analysis\';
color_code = 'k';
%data_file = 'across_sessions_RC_data_CS1_naive_lickAlign.mat';
data_file = 'across_sessions_RC_data_interleaved_lickAlign.mat';

% target align average
firstlick_align_events_audi_all_session = [];% frames*total number of sessions 
firstlick_align_events_vis_all_session = [];
firstlick_align_events_av_all_session = [];
lickBurstAlign_events_audi_all_session = [];
lickBurstAlign_events_vis_all_session = [];
lickBurstAlign_events_av_all_session = [];

firstLickAlign_audi_all_session = [];
firstLickAlign_vis_all_session = [];
firstLickAlign_av_all_session = [];

lickBurstAlign_postRew_audi_all_session = [];
lickBurstAlign_postRew_vis_all_session = [];
lickBurstAlign_postRew_av_all_session = [];

nPCs_total = 0;
ntrials_total = 0;
for ii = 1:length(sessions) 
    %load things  
    lickAlign_output = load((fullfile([analysis_out sessions{ii} '\' sessions{ii} '_' run{ii} '_cueAlignLick.mat'])));
    
    firstLickAlign_events_thisSession = lickAlign_output.firstPostRew_lickAlignEvents_1500_3000ms_scale; % frame*cell*trial
    firstLickAlign_thisSession = lickAlign_output.firstPostRew_lickAlign_1500_3000ms_scale; % frame*trial
    
    lickburstAlign_postRew_events_thisSession = lickAlign_output.postRew_lickAlignEvents_1500_3000ms_scale;
    lickbursAlign_postRew_thisSession = lickAlign_output.postRew_lickAlign_1500_3000ms_scale;
    
    trailTypes = load (['Z:\behavior_analysis\RC\' sessions{ii} '_IL_trial_inx.mat']);
    auditrials = trailTypes.auditrials;
    vistrials = trailTypes.vistrials;
    avtrials = trailTypes.avtrials;
    
    nPCs_total = nPCs_total + size(firstLickAlign_events_thisSession,2);
    ntrials_total = ntrials_total + size(firstLickAlign_events_thisSession,3);
    
    %for spike data, average across trials before cat - need to do this for different trial types
    mean_firstLickAlign_events_audi = squeeze(nanmean(firstLickAlign_events_thisSession(:,:,auditrials),3));
    firstlick_align_events_audi_all_session = cat(2,firstlick_align_events_audi_all_session,mean_firstLickAlign_events_audi);
    
    mean_firstLickAlign_events_vis = squeeze(nanmean(firstLickAlign_events_thisSession(:,:,vistrials),3));
    firstlick_align_events_vis_all_session = cat(2,firstlick_align_events_vis_all_session,mean_firstLickAlign_events_vis);
    
    mean_firstLickAlign_events_av = squeeze(nanmean(firstLickAlign_events_thisSession(:,:,avtrials),3));
    firstlick_align_events_av_all_session = cat(2,firstlick_align_events_av_all_session,mean_firstLickAlign_events_av);
    
    mean_lickBurstAlign_events_audi = squeeze(nanmean(lickburstAlign_postRew_events_thisSession(:,:,auditrials),3));
    lickBurstAlign_events_audi_all_session = cat(2,lickBurstAlign_events_audi_all_session,mean_lickBurstAlign_events_audi);
    
    mean_lickBurstAlign_events_vis = squeeze(nanmean(lickburstAlign_postRew_events_thisSession(:,:,vistrials),3));
    lickBurstAlign_events_vis_all_session = cat(2,lickBurstAlign_events_vis_all_session,mean_lickBurstAlign_events_vis);
    
    mean_lickBurstAlign_events_av = squeeze(nanmean(lickburstAlign_postRew_events_thisSession(:,:,avtrials),3));
    lickBurstAlign_events_av_all_session = cat(2,lickBurstAlign_events_av_all_session,mean_lickBurstAlign_events_av);
    
    % average licking data
    firstLickAlign_audi_all_session = cat(2,firstLickAlign_audi_all_session,firstLickAlign_thisSession(:,auditrials));
    firstLickAlign_vis_all_session = cat(2,firstLickAlign_vis_all_session,firstLickAlign_thisSession(:,vistrials));
    firstLickAlign_av_all_session = cat(2,firstLickAlign_av_all_session,firstLickAlign_thisSession(:,avtrials));
    
    lickBurstAlign_postRew_audi_all_session = cat(2,lickBurstAlign_postRew_audi_all_session,lickbursAlign_postRew_thisSession(:,auditrials));
    lickBurstAlign_postRew_vis_all_session = cat(2,lickBurstAlign_postRew_vis_all_session,lickbursAlign_postRew_thisSession(:,vistrials));
    lickBurstAlign_postRew_av_all_session = cat(2,lickBurstAlign_postRew_av_all_session,lickbursAlign_postRew_thisSession(:,avtrials));    
    
end

tt = lickAlign_output.tt;

ave_spike_firstLickAlign_audi_acrossessions = nanmean(firstlick_align_events_audi_all_session,2)*30;
ste_spike_firstLickAlign_audi_acrossessions = nanstd(firstlick_align_events_audi_all_session,0,2)*30/sqrt(size(firstlick_align_events_audi_all_session,2));

ave_spike_firstLickAlign_vis_acrossessions = nanmean(firstlick_align_events_vis_all_session,2)*30;
ste_spike_firstLickAlign_vis_acrossessions = nanstd(firstlick_align_events_vis_all_session,0,2)*30/sqrt(size(firstlick_align_events_vis_all_session,2));

ave_spike_firstLickAlign_av_acrossessions = nanmean(firstlick_align_events_av_all_session,2)*30;
ste_spike_firstLickAlign_av_acrossessions = nanstd(firstlick_align_events_av_all_session,0,2)*30/sqrt(size(firstlick_align_events_av_all_session,2));

ave_spike_lickBurstAlign_audi_acrossessions = nanmean(lickBurstAlign_events_audi_all_session,2)*30;
ste_spike_lickBurstAlign_audi_acrossessions = nanstd(lickBurstAlign_events_audi_all_session,0,2)*30/sqrt(size(lickBurstAlign_events_audi_all_session,2));

ave_spike_lickBurstAlign_vis_acrossessions = nanmean(lickBurstAlign_events_vis_all_session,2)*30;
ste_spike_lickBurstAlign_vis_acrossessions = nanstd(lickBurstAlign_events_vis_all_session,0,2)*30/sqrt(size(lickBurstAlign_events_vis_all_session,2));

ave_spike_lickBurstAlign_av_acrossessions = nanmean(lickBurstAlign_events_av_all_session,2)*30;
ste_spike_lickBurstAlign_av_acrossessions = nanstd(lickBurstAlign_events_av_all_session,0,2)*30/sqrt(size(lickBurstAlign_events_av_all_session,2));

ave_firstLickAlign_audi_acrossessions = nanmean(firstLickAlign_audi_all_session,2)*30;
ste_firstLickAlign_audi_acrossessions = nanstd(firstLickAlign_audi_all_session,0,2)*30/sqrt(size(firstLickAlign_audi_all_session,2));

ave_firstLickAlign_vis_acrossessions = nanmean(firstLickAlign_vis_all_session,2)*30;
ste_firstLickAlign_vis_acrossessions = nanstd(firstLickAlign_vis_all_session,0,2)*30/sqrt(size(firstLickAlign_vis_all_session,2));

ave_firstLickAlign_av_acrossessions = nanmean(firstLickAlign_av_all_session,2)*30;
ste_firstLickAlign_av_acrossessions = nanstd(firstLickAlign_av_all_session,0,2)*30/sqrt(size(firstLickAlign_av_all_session,2));

ave_lickBurstAlign_audi_acrossessions = nanmean(lickBurstAlign_postRew_audi_all_session,2)*30;
ste_lickBurstAlign_audi_acrossessions = nanstd(lickBurstAlign_postRew_audi_all_session,0,2)*30/sqrt(size(lickBurstAlign_postRew_audi_all_session,2));

ave_lickBurstAlign_vis_acrossessions = nanmean(lickBurstAlign_postRew_vis_all_session,2)*30;
ste_lickBurstAlign_vis_acrossessions = nanstd(lickBurstAlign_postRew_vis_all_session,0,2)*30/sqrt(size(lickBurstAlign_postRew_vis_all_session,2));

ave_lickBurstAlign_av_acrossessions = nanmean(lickBurstAlign_postRew_av_all_session,2)*30;
ste_lickBurstAlign_av_acrossessions = nanstd(lickBurstAlign_postRew_av_all_session,0,2)*30/sqrt(size(lickBurstAlign_postRew_av_all_session,2));


save(fullfile(['Z:\2P_analysis\across_sessions_RC\',data_file]),...
    'ave_spike_firstLickAlign_audi_acrossessions','ste_spike_firstLickAlign_audi_acrossessions',...
    'ave_spike_firstLickAlign_vis_acrossessions','ste_spike_firstLickAlign_vis_acrossessions',...
    'ave_spike_firstLickAlign_av_acrossessions','ste_spike_firstLickAlign_av_acrossessions',...
    'ave_spike_lickBurstAlign_audi_acrossessions','ste_spike_lickBurstAlign_audi_acrossessions',...
    'ave_spike_lickBurstAlign_vis_acrossessions','ste_spike_lickBurstAlign_vis_acrossessions',...
    'ave_spike_lickBurstAlign_av_acrossessions','ste_spike_lickBurstAlign_av_acrossessions',...
    'ave_firstLickAlign_audi_acrossessions','ste_firstLickAlign_audi_acrossessions',...
    'ave_firstLickAlign_vis_acrossessions','ste_firstLickAlign_vis_acrossessions',...
    'ave_firstLickAlign_av_acrossessions','ste_firstLickAlign_av_acrossessions',...
    'ave_lickBurstAlign_audi_acrossessions','ste_lickBurstAlign_audi_acrossessions',...
    'ave_lickBurstAlign_vis_acrossessions','ste_lickBurstAlign_vis_acrossessions',...
    'ave_lickBurstAlign_av_acrossessions','ste_lickBurstAlign_av_acrossessions',...
    'tt','nPCs_total','ntrials_total');


% plot the above section
firstLickAlign_fig = figure; 
firstLickAlign_fig.Units = 'inches';
firstLickAlign_fig.Position = [1 1 9 5];
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',10);
subplot(2,3,1);
shadedErrorBar(tt, ave_spike_firstLickAlign_audi_acrossessions, ste_spike_firstLickAlign_audi_acrossessions,'k');
ylim([0 3.5]);vline(0,'k');ylabel('spike rate (Hz)');
title('CS1');box off;
subplot(2,3,4);
shadedErrorBar(tt, ave_firstLickAlign_audi_acrossessions, ste_firstLickAlign_audi_acrossessions,'k');
ylim([0 40]);vline(0,'k');ylabel('lick rate (Hz)');xlabel('time from first lick (ms)');
box off;
subplot(2,3,2);
shadedErrorBar(tt, ave_spike_firstLickAlign_av_acrossessions, ste_spike_firstLickAlign_av_acrossessions,'k');
ylim([0 3.5]);vline(0,'k');ylabel('spike rate (Hz)');
title('CS1+CS2');box off;
subplot(2,3,5);
shadedErrorBar(tt, ave_firstLickAlign_av_acrossessions, ste_firstLickAlign_av_acrossessions,'k');
ylim([0 40]);vline(0,'k');ylabel('lick rate (Hz)');xlabel('time from first lick (ms)');
box off;
subplot(2,3,3);
shadedErrorBar(tt, ave_spike_firstLickAlign_vis_acrossessions, ste_spike_firstLickAlign_vis_acrossessions,'k');
ylim([0 3.5]);vline(0,'k');ylabel('spike rate (Hz)');
title('CS2');box off;
subplot(2,3,6);
shadedErrorBar(tt, ave_firstLickAlign_vis_acrossessions, ste_firstLickAlign_vis_acrossessions,'k');
ylim([0 40]);vline(0,'k');ylabel('lick rate (Hz)');xlabel('time from first lick (ms)');
box off;
supertitle('align to post reward first lick interleaved day');
path = 'Z:\2P_analysis\across_sessions_RC\';
fig_name = 'across_session_firstLickAlign_CS1CS2_interleave';
print(firstLickAlign_fig,[path,fig_name],'-r600','-dpdf');
print(firstLickAlign_fig,[path,fig_name],'-djpeg');
savefig([path fig_name]);

lickBurstAlign_fig = figure; 
lickBurstAlign_fig.Units = 'inches';
lickBurstAlign_fig.Position = [1 1 9 5];
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',10);
subplot(2,3,1);
shadedErrorBar(tt, ave_spike_lickBurstAlign_audi_acrossessions, ste_spike_lickBurstAlign_audi_acrossessions,'k');
ylim([0 3.5]);vline(0,'k');ylabel('spike rate (Hz)');
title('CS1');box off;
subplot(2,3,4);
shadedErrorBar(tt, ave_lickBurstAlign_audi_acrossessions, ste_lickBurstAlign_audi_acrossessions,'k');
ylim([0 40]);vline(0,'k');ylabel('lick rate (Hz)');xlabel('time from first lick (ms)');
box off;
subplot(2,3,2);
shadedErrorBar(tt, ave_spike_lickBurstAlign_av_acrossessions, ste_spike_lickBurstAlign_av_acrossessions,'k');
ylim([0 3.5]);vline(0,'k');ylabel('spike rate (Hz)');
title('CS1+CS2');box off;
subplot(2,3,5);
shadedErrorBar(tt, ave_lickBurstAlign_av_acrossessions, ste_lickBurstAlign_av_acrossessions,'k');
ylim([0 40]);vline(0,'k');ylabel('lick rate (Hz)');xlabel('time from first lick (ms)');
box off;
subplot(2,3,3);
shadedErrorBar(tt, ave_spike_lickBurstAlign_vis_acrossessions, ste_spike_lickBurstAlign_vis_acrossessions,'k');
ylim([0 3.5]);vline(0,'k');ylabel('spike rate (Hz)');
title('CS2');box off;
subplot(2,3,6);
shadedErrorBar(tt, ave_lickBurstAlign_vis_acrossessions, ste_lickBurstAlign_vis_acrossessions,'k');
ylim([0 40]);vline(0,'k');ylabel('lick rate (Hz)');xlabel('time from first lick (ms)');
box off;
supertitle('align to post reward lick burst interleaved day');
path = 'Z:\2P_analysis\across_sessions_RC\';
fig_name = 'across_session_lickBurstAlign_CS1CS2_interleave';
print(firstLickAlign_fig,[path,fig_name],'-r600','-dpdf');
print(firstLickAlign_fig,[path,fig_name],'-djpeg');
savefig([path fig_name]);


%% align licking data with cue across sessions
% assign document paths and experimental sessions
clear;
% CS1 naive
% sessions = {'201113_img1078','201114_img1079','210112_img1083','210113_img1082',...
% '210405_img1087','210412_img1090','210416_img1088'};
% CS1 trained
% sessions = {'201125_img1078','201125_img1079','210119_img1083','210416_img1087',...
% '210427_img1090','210428_img1088'};
% CS1 test day after training with CS1+CS2 %!!!!!!!!!!!!!!!!!201209_img1079 the run is 001
% sessions = {'201211_img1078','201209_img1079','210203_img1083','210429_img1087',...
%     '210508_img1088','210507_img1090'};
% CS2 test day after training with CS1+CS2
% sessions = {'201210_img1078','201210_img1079','210205_img1083','210430_img1087',...
% '210509_img1088','210511_img1090'};
% CS1+CS2 test day after training with CS1+CS2
sessions = {'201208_img1078','201209_img1079','210202_img1083','210428_img1087',...
    '210507_img1088','210506_img1090'};

% control group- 1901 and 1092
% naive day
% sessions = {'210616_img1091','210616_img1092'};
% trained day 9
% sessions = {'210628_img1091','210628_img1092'};

run = {'000','000','000','000','000','000','000'}; 
% run = {'000','001','000','000','000','000','000'}; %!!!!!!!!!!!!!!!!!201209_img1079 the run is 001
analysis_out    = 'Z:\2P_analysis\';
color_code = 'b';
% fig_name = 'across_session_targetalign_lick_CS1_naive';
% data_file = 'across_sessions_RC_data_CS1_naive.mat';
% fig_name = 'across_session_targetalign_lick_CS1_naive_control';
% data_file = 'across_sessions_RC_data_CS1_naive_control.mat';
% fig_name = 'across_session_targetalign_lick_CS1_trained';
% data_file = 'across_sessions_RC_data_CS1_trained.mat';
% fig_name = 'across_session_targetalign_lick_CS1_trained_control';
% data_file = 'across_sessions_RC_data_CS1_trained_control.mat';
% fig_name = 'across_session_targetalign_lick_CS1CS2_testCS1';
% data_file = 'across_sessions_RC_data_CS1CS2_testCS1.mat';
% fig_name = 'across_session_targetalign_lick_CS1CS2_testCS2';
% data_file = 'across_sessions_RC_data_CS1CS2_testCS2.mat';
fig_name = 'across_session_targetalign_lick_CS1CS2_testCS1CS2';
data_file = 'across_sessions_RC_data_CS1CS2_testCS1CS2.mat';
% lick align average
target_align_licks_all_session = [];
ntrials_total = 0;
for ii = 1:length(sessions) 
    %load things  
    lickAlign_output = load((fullfile([analysis_out sessions{ii} '\' sessions{ii} '_' run{ii} '_cueAlignLick.mat'])));
    lickAlign_thisSession = lickAlign_output.lickCueAlign; % frame*trial
    ntrials_total = ntrials_total + size(lickAlign_thisSession,2);
    
    % average across trials
    mean_targetAlign_licks = squeeze(nanmean(lickAlign_thisSession,2));
    target_align_licks_all_session = cat(2,target_align_licks_all_session,mean_targetAlign_licks);
end

tt = lickAlign_output.tt;

ave_lick_targetAlign_acrossessions = mean(target_align_licks_all_session,2)*30;
ste_lick_targetAlign_acrossessions = std(target_align_licks_all_session,0,2)*30/sqrt(size(target_align_licks_all_session,2));

save(fullfile(['Z:\2P_analysis\across_sessions_RC\',data_file]),...
    'ave_lick_targetAlign_acrossessions','ste_lick_targetAlign_acrossessions','tt',...
    'ntrials_total','target_align_licks_all_session','-append');


% plot the above section
lickAlign_fig = figure; 
lickAlign_fig.Units = 'inches';
lickAlign_fig.Position = [1 1 4 2.5];
shadedErrorBar(tt, ave_lick_targetAlign_acrossessions, ste_lick_targetAlign_acrossessions);
ylim([0 10]);vline(0,color_code);vline(700,'k');
ylabel('Lick rate (Hz)');
%title('CS1 naive day, across sessions');
%title('CS1 trained day, across sessions');
%title('CS1+CS2 trained test w/ CS1, across sessions');
%title('CS1+CS2 trained test w/ CS2, across sessions');
title('CS1+CS2 trained test w/ CS1+CS2, across sessions');
xlabel('Time from cue (ms)');
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'FontSize',7);
%xlim([-1.1 1.6]);ylim([0 3]);
box off;
path = 'Z:\2P_analysis\across_sessions_RC\';
print(lickAlign_fig,[path,fig_name],'-r600','-dpdf');
print(lickAlign_fig,[path,fig_name],'-djpeg');
savefig([path fig_name]);


%% align neural data with cue onset on learned days, high lick rate trials vs. low lick rate trials 
% results for this section come from RC_lickAlign_SJ
clear;
% CS1 naive
% sessions = {'201113_img1078','201114_img1079','210112_img1083',...
% '210405_img1087','210412_img1090','210416_img1088'};
% CS1 trained
% sessions = {'201125_img1078','201125_img1079','210119_img1083','210416_img1087',...
% '210427_img1090','210428_img1088'};

% control group- 1901 and 1092
% naive day
% sessions = {'210616_img1091','210616_img1092'};
% trained day 9
 sessions = {'210628_img1091','210628_img1092'};

run = {'000','000','000','000','000','000','000'};
analysis_out    = 'Z:\2P_analysis\';
color_code = 'b';
% fig_name = 'across_session_CS1_naive_byLickRate';
% data_file = 'across_sessions_RC_data_CS1_naive.mat';
% fig_name = 'across_session_CS1_naive_control_byLickRate';
% data_file = 'across_sessions_RC_data_CS1_naive_control.mat';
% fig_name = 'across_session_CS1_trained_byLickRate';
% data_file = 'across_sessions_RC_data_CS1_trained.mat';
fig_name = 'across_session_CS1_trained_control_byLickRate';
data_file = 'across_sessions_RC_data_CS1_trained_control.mat';
targetAlign_events_highlick_postrew = [];
targetAlign_events_lowlick_postrew = [];
nPCs_total = 0;
highlick_postrew_sessions = [];
lowlick_postrew_sessions = [];
    
for ii = 1:length(sessions) 
    %load things  
    targetAlign_output = load((fullfile([analysis_out sessions{ii} '\' sessions{ii} '_' run{ii} '_targetAlign.mat'])));
    targetAlign_events_thisSession = targetAlign_output.targetAlign_events; % frame*cell*trial
    prewin_frames = targetAlign_output.prewin_frames;
    frameRateHz = targetAlign_output.frameRateHz;
    nPCs_total = nPCs_total + size(targetAlign_events_thisSession,2);
    lickAlign_output = load((fullfile([analysis_out sessions{ii} '\' sessions{ii} '_' run{ii} '_cueAlignLick.mat'])));
    ind_low_postrew_thisSession = lickAlign_output.ind_low_postrew; % index of the early and late burst trials
    ind_high_postrew_thisSession = lickAlign_output.ind_high_postrew;
    %need to average across trials before cat, otherwise # of cells and
    %trials are both different across sessions, can't cat
    targetAlign_events_lowlick_postrew_thisSession = targetAlign_events_thisSession(:,:,ind_low_postrew_thisSession);
    targetAlign_events_lowlick_postrew_cells = nanmean(targetAlign_events_lowlick_postrew_thisSession,3);
    targetAlign_events_lowlick_postrew = cat(2,targetAlign_events_lowlick_postrew,targetAlign_events_lowlick_postrew_cells);
    
    targetAlign_events_highlick_postrew_thisSession = targetAlign_events_thisSession(:,:,ind_high_postrew_thisSession);
    targetAlign_events_highlick_postrew_cells = nanmean(targetAlign_events_highlick_postrew_thisSession,3);
    targetAlign_events_highlick_postrew = cat(2,targetAlign_events_highlick_postrew,targetAlign_events_highlick_postrew_cells);
    
    postRew_lickBurstHz_thisSession = lickAlign_output.postRew_lickBurstHz; % lick burst rate after reward delivery in each trial
    lowlick_postrew_sessions = cat(2,lowlick_postrew_sessions,postRew_lickBurstHz_thisSession(ind_low_postrew_thisSession));
    highlick_postrew_sessions = cat(2,highlick_postrew_sessions,postRew_lickBurstHz_thisSession(ind_high_postrew_thisSession));    
    
end
tt = lickAlign_output.tt;
ave_low_postrew_BstHz = mean(lowlick_postrew_sessions);
ave_high_postrew_BstHz = mean(highlick_postrew_sessions);

neural_by_lick_rate = figure; %plot neural data of trials of high vs. low lick burst rate
neural_by_lick_rate.Units = 'inches';
neural_by_lick_rate.Position = [1 1 4 2.5];
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',10);
shadedErrorBar(tt,nanmean(targetAlign_events_lowlick_postrew,2).*(1000./frameRateHz), (nanstd(targetAlign_events_lowlick_postrew,0,2)).*(1000./frameRateHz)./sqrt(nPCs_total), 'b');
hold on;
shadedErrorBar(tt,nanmean(targetAlign_events_highlick_postrew,2).*(1000./frameRateHz), (nanstd(targetAlign_events_highlick_postrew,0,2)).*(1000./frameRateHz)./sqrt(nPCs_total), 'k');
xlabel('Time from cue (ms)');
ylabel('Spike rate (Hz)');
ylim([0 2.5]);
vline(0,color_code);vline(700,'k');
%title(['CS1 naive, - lick bursts: high rate = ' num2str(ave_high_postrew_BstHz) 'Hz; low rate = ' num2str(ave_low_postrew_BstHz) 'Hz']);
title(['CS1 trained, - lick bursts: high rate = ' num2str(ave_high_postrew_BstHz) 'Hz; low rate = ' num2str(ave_low_postrew_BstHz) 'Hz']);
hold off;
box off;
path = 'Z:\2P_analysis\across_sessions_RC\';
print(neural_by_lick_rate,[path,fig_name],'-r600','-dpdf');
print(neural_by_lick_rate,[path,fig_name],'-djpeg');
savefig([path fig_name]);

save(fullfile(['Z:\2P_analysis\across_sessions_RC\',data_file]),...
    'targetAlign_events_lowlick_postrew','targetAlign_events_highlick_postrew',...
    'lowlick_postrew_sessions','highlick_postrew_sessions',...
    'nPCs_total','-append');

%% align neural data and piezo data with cue onset, most 10% movement vs. least 10% movement during the pre-reward window 
clear;

% CS1 naive
% sessions = {'210112_img1083','210113_img1082','210405_img1087','210412_img1090','210416_img1088'};
% CS1 trained
% sessions = {'201125_img1078','201125_img1079','210119_img1083','210416_img1087',...
% '210427_img1090','210428_img1088'};

% control group- 1901 and 1092
% naive day
% sessions = {'210616_img1091','210616_img1092'};
% trained day 9
 sessions = {'210628_img1091','210628_img1092'};

run = {'000','000','000','000','000','000'};
analysis_out    = 'Z:\2P_analysis\';
color_code = 'b';
% fig_name = 'across_session_CS1_naive_bypiezo10';
% data_file = 'across_sessions_RC_data_CS1_naive.mat';
% fig_name = 'across_session_CS1_naive_control_bypiezo10';
% data_file = 'across_sessions_RC_data_CS1_naive_control.mat';
% fig_name = 'across_session_CS1_trained_bypiezo10';
% data_file = 'across_sessions_RC_data_CS1_trained.mat';
fig_name = 'across_session_CS1_trained_control_bypiezo10';
data_file = 'across_sessions_RC_data_CS1_trained_control.mat';
targetAlign_events_mostmove_prerew = [];
targetAlign_events_lessmove_prerew = [];
nPCs_total = 0;
mostmove_prerew_sessions = [];
lessmove_prerew_sessions = [];
    
for ii = 1:length(sessions) 
    %load things  
    targetAlign_output = load((fullfile([analysis_out sessions{ii} '\' sessions{ii} '_' run{ii} '_targetAlign.mat'])));
    targetAlign_events_thisSession = targetAlign_output.targetAlign_events; % frame*cell*trial
    prewin_frames = targetAlign_output.prewin_frames;
    frameRateHz = targetAlign_output.frameRateHz;
    nPCs_total = nPCs_total + size(targetAlign_events_thisSession,2);
    piezoAlign_output = load((fullfile([analysis_out sessions{ii} '\' sessions{ii} '_' run{ii} '_cueAlignPiezo.mat'])));
    ind_low10piezo_prerew_thisSession = piezoAlign_output.ind_low10piezo_prerew;% index of the trials that have most and least movement
    ind_high10piezo_prerew_thisSession = piezoAlign_output.ind_high10piezo_prerew;
    %need to average across trials before cat, otherwise # of cells and
    %trials are both different across sessions, can't cat
    targetAlign_events_lessmove_prerew_thisSession = targetAlign_events_thisSession(:,:,ind_low10piezo_prerew_thisSession);
    targetAlign_events_lessmove_prerew_cells = nanmean(targetAlign_events_lessmove_prerew_thisSession,3);
    targetAlign_events_lessmove_prerew = cat(2,targetAlign_events_lessmove_prerew,targetAlign_events_lessmove_prerew_cells);
    
    targetAlign_events_mostmove_prerew_thisSession = targetAlign_events_thisSession(:,:,ind_high10piezo_prerew_thisSession);
    targetAlign_events_mostmove_prerew_cells = nanmean(targetAlign_events_mostmove_prerew_thisSession,3);
    targetAlign_events_mostmove_prerew = cat(2,targetAlign_events_mostmove_prerew,targetAlign_events_mostmove_prerew_cells);
    
    targetAlign_piezo_thisSession = piezoAlign_output.targetAlign_piezo; % lick burst rate after reward delivery in each trial
    lessmove_prerew_sessions = cat(2,lessmove_prerew_sessions,targetAlign_piezo_thisSession(:,ind_low10piezo_prerew_thisSession));
    mostmove_prerew_sessions = cat(2,mostmove_prerew_sessions,targetAlign_piezo_thisSession(:,ind_high10piezo_prerew_thisSession));    
    
end
tt = targetAlign_output.tt;

neural_by_piezo = figure; %plot neural data of trials of early vs. late bursts
neural_by_piezo.Units = 'inches';
neural_by_piezo.Position = [1 1 3.5 3];
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',10);
subplot(2,1,1);
shadedErrorBar(tt,nanmean(lessmove_prerew_sessions,2).*(1000./frameRateHz), (nanstd(lessmove_prerew_sessions,0,2)).*(1000./frameRateHz)./sqrt(size(lessmove_prerew_sessions,2)), 'b');
hold on;
shadedErrorBar(tt,nanmean(mostmove_prerew_sessions,2).*(1000./frameRateHz), (nanstd(mostmove_prerew_sessions,0,2)).*(1000./frameRateHz)./sqrt(size(lessmove_prerew_sessions,2)), 'k');
ylabel('Piezo voltage (V)');
ylim([0 10]);
vline(0,color_code);vline(700,'k');
%title('CS1 naive, - black: most 10% movement, blue: least 10% movement' );
title('CS1 trained, - black: most 10% movement, blue: least 10% movement' );
box off;
subplot(2,1,2);
shadedErrorBar(tt,nanmean(targetAlign_events_lessmove_prerew,2).*(1000./frameRateHz), (nanstd(targetAlign_events_lessmove_prerew,0,2)).*(1000./frameRateHz)./sqrt(nPCs_total), 'b');
hold on;
shadedErrorBar(tt,nanmean(targetAlign_events_mostmove_prerew,2).*(1000./frameRateHz), (nanstd(targetAlign_events_mostmove_prerew,0,2)).*(1000./frameRateHz)./sqrt(nPCs_total), 'k');
xlabel('Time from cue (ms)');
ylabel('Spike rate (Hz)');
ylim([0 2.5]);
vline(0,color_code);vline(700,'k');
hold off;
box off;
path = 'Z:\2P_analysis\across_sessions_RC\';
print(neural_by_piezo,[path,fig_name],'-r600','-dpdf');
print(neural_by_piezo,[path,fig_name],'-djpeg');
savefig([path fig_name]);

save(fullfile(['Z:\2P_analysis\across_sessions_RC\',data_file]),...
    'targetAlign_events_lessmove_prerew','targetAlign_events_mostmove_prerew',...
    'lessmove_prerew_sessions','mostmove_prerew_sessions',...
    'nPCs_total','-append');


