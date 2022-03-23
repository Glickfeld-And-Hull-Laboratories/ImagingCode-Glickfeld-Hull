%% plot for RT, all trials, across animals
clear;
behav_analysis_out = 'Z:\home\shuyang\behavior_analysis\RC\';
% % sessions = {'201113_img1078-1308','201113_img1078-1341','201114_img1078','201115_img1078',...
%      '201116_img1078','201117_img1078','201118_img1078','201119_img1078','201120_img1078',...
%      '201123_img1078','201124_img1078','201125_img1078','201127_img1078','201130_img1078',...
%      '201201_img1078','201202_img1078','201203_img1078','201204_img1078','201206_img1078',...
%      '201208_img1078','201210_img1078-1134','201211_img1078-1251','201211_img1078-1326',...
%      '201216_img1078','201217_img1078'};
% testDay_mouse = [1078,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,2];
% sessions = {'201114_img1079','201115_img1079','201116_img1079','201117_img1079',...
%     '201118_img1079','201119_img1079','201120_img1079','201123_img1079','201124_img1079',...
%     '201125_img1079','201127_img1079','201129_img1079','201130_img1079','201201_img1079',...
%     '201202_img1079','201206_img1079','201207_img1079','201209_img1079-1251',...
%     '201209_img1079-1325','201210_img1079-1252','201216_img1079','201211_img1079'};
% sessions = {'210112_img1083','210113_img1083','210114_img1083','210115_img1083',...
%     '210117_img1083','210118_img1083','210119_img1083','210121_img1083','210122_img1083',...
%     '210124_img1083','210125_img1083','210126_img1083','210128_img1083','210129_img1083',...
%     '210131_img1083','210201_img1083','210202_img1083','210203_img1083','210204_img1083-1141',...
%     '210205_img1083-1140','210209_img1083'};
% sessions = {'210405_img1087','210406_img1087','210407_img1087','210408_img1087',...
%     '210409_img1087','210411_img1087','210412_img1087','210413_img1087','210414_img1087',...
%     '210415_img1087','210416_img1087','210419_img1087','210420_img1087','210421_img1087',...
%     '210422_img1087','210423_img1087','210426_img1087','210427_img1087','210428_img1087',...
%     '210429_img1087','210430_img1087','210503_img1087'};
% testDay_mouse = [1087,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,2];
% sessions = {'210416_img1088','210419_img1088','210420_img1088','210421_img1088',...
%     '210422_img1088','210423_img1088','210426_img1088','210427_img1088','210428_img1088',...
%     '210429_img1088','210430_img1088','210501_img1088','210502_img1088','210503_img1088',...
%     '210504_img1088','210505_img1088','210506_img1088','210507_img1088','210508_img1088',...
%     '210509_img1088','210511_img1088'};
% sessions = {'210412_img1090','210413_img1090','210414_img1090','210415_img1090',...
%     '210416_img1090','210419_img1090','210420_img1090','210422_img1090',...
%     '210423_img1090','210426_img1090','210427_img1090','210428_img1090',...
%     '210429_img1090','210430_img1090','210501_img1090','210502_img1090',...
%     '210503_img1090','210504_img1090','210505_img1090','210506_img1090',...
%     '210507_img1090','210508_img1090','210509_img1090','210510_img1090',...
%     '210511_img1090'};
% from naive day to CS1 imaging day: 
ID1078 = {'201113_img1078-1308','201113_img1078-1341','201114_img1078','201115_img1078',...
      '201116_img1078','201117_img1078','201118_img1078','201119_img1078','201120_img1078',...
      '201123_img1078','201124_img1078','201125_img1078'};
ID1079 = {'201114_img1079','201115_img1079','201116_img1079','201117_img1079',...
     '201118_img1079','201119_img1079','201120_img1079','201123_img1079','201124_img1079',...
     '201125_img1079'};
ID1083 = {'210112_img1083','210113_img1083','210114_img1083','210115_img1083',...
     '210117_img1083','210118_img1083','210119_img1083'};
ID1087 = {'210405_img1087','210406_img1087','210407_img1087','210408_img1087',...
     '210409_img1087','210411_img1087','210412_img1087','210413_img1087','210414_img1087',...
     '210415_img1087','210416_img1087'};
ID1088 = {'210416_img1088','210419_img1088','210420_img1088','210421_img1088',...
     '210422_img1088','210423_img1088','210426_img1088','210427_img1088','210428_img1088'};
ID1090 = {'210412_img1090','210413_img1090','210414_img1090','210415_img1090',...
     '210416_img1090','210419_img1090','210420_img1090','210422_img1090',...
     '210423_img1090','210426_img1090','210427_img1090'};
IDs = {ID1078,ID1079,ID1083,ID1087,ID1088,ID1090};
RT_all_animals = {};
RT_all_lengths = zeros(1,length(IDs));
for i = 1:length(IDs)
    % stack all RTs from each animal into a vector first
    animal = IDs{i};
    RT_all_animals{i} = [];
    for s = 1:length(animal)
        rawRT = load((fullfile([behav_analysis_out '\RT_data\' animal{s} '.mat'])));
        RT_this_session = rawRT.RT_this_session_raw;
        RT_this_session(RT_this_session<0) = 0;
        RT_all_animals{i} = [RT_all_animals{i} RT_this_session];
    end
    RT_all_lengths(i) = length(RT_all_animals{i});
end
RT_all_animals_mat = nan(length(IDs),max(RT_all_lengths));
for i = 1:length(IDs)
    RT_all_animals_mat(i,1:length(RT_all_animals{i})) = RT_all_animals{i};
end

x = 1:1:max(RT_all_lengths);
colorcode = brewermap(6,'Dark2');
RTall_fig = figure;
RTall_fig.Units = 'inches';
RTall_fig.Position = [1 1 4 2.5];
% for i = 1:length(IDs)
%     scatter(x,RT_all_animals_mat(i,:),8,'filled','MarkerEdgeColor',colorcode(i,:),'MarkerFaceColor',colorcode(i,:));hold on;
% end
errorbar(x,nanmean(RT_all_animals_mat),nanstd(RT_all_animals_mat,0,1)/sqrt(length(IDs)),'Color',[0.5 0.5 0.5], 'LineStyle', 'none', 'CapSize', 1);
xlabel('trial number');
ylabel('response time (ms)');
title('CS1 training, across animals');
box off;
ylim([0 2500]);
% path = 'Z:\home\shuyang\2P_analysis\forRO1\';
% print(RTall_fig,[path,fig_name1],'-r600','-dpdf');
% print(RTall_fig,[path,fig_name1],'-djpeg');
% savefig([path fig_name1]);


%% peak lick rate on each trial
clear;
% CS1 naive
sessions = {'201113_img1078','201114_img1079','210112_img1083',...
'210405_img1087','210416_img1088','210412_img1090'};
% CS1 trained
% sessions = {'201125_img1078','201125_img1079','210119_img1083','210416_img1087',...
%     '210427_img1090','210428_img1088'};
run = {'000','000','000','000','000','000','000'};
% run = {'000','001','000','000','000','000','000'}; %!!!!!!!!!!!!!!!!!201209_img1079 the run is 001
analysis_out    = 'Z:\home\shuyang\2P_analysis\';
color_code = 'b';
fig_name1 = 'across_session_lickpeaktime_CS1_naive';
fig_name2 = 'across_session_lickpeakrate_CS1_naive';
% fig_name1 = 'across_session_lickpeaktime_CS1_trained';
% fig_name2 = 'across_session_lickpeakrate_CS1_trained';

max_tt_all = [];
max_hz_all = [];

for ii = 1:length(sessions)
    %load things
    lickAlign_output = load((fullfile([analysis_out sessions{ii} '\' sessions{ii} '_' run{ii} '_cueAlignLick.mat'])));    lickAlign_thisSession = lickAlign_output.lickCueAlign; % frame*trial
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
title('CS1 naive day, across sessions');
%title('CS1 trained day, across sessions');
%title('CS1+CS2 trained test w/ CS1, across sessions');
%title('CS1+CS2 trained test w/ CS2, across sessions');
%title('CS1+CS2 trained test w/ CS1+CS2, across sessions');
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
title('CS1 naive day, across sessions');
%title('CS1 trained day, across sessions');
%title('CS1+CS2 trained test w/ CS1, across sessions');
%title('CS1+CS2 trained test w/ CS2, across sessions');
%title('CS1+CS2 trained test w/ CS1+CS2, across sessions');
xlabel('lick rates (Hz)');
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'FontSize',7);
%xlim([-1.1 1.6]);ylim([0 3]);
box off;
path = 'Z:\home\shuyang\2P_analysis\forRO1\';
print(maxhzcum_fig,[path,fig_name2],'-r600','-dpdf');
print(maxhzcum_fig,[path,fig_name2],'-djpeg');
savefig([path fig_name2]);

%% plot cumulative distribution function plot for licking data
% assign document paths and experimental sessions
clear;
% CS1 naive
% sessions = {'201113_img1078-1308','201114_img1079','210112_img1083','210405_img1087','210412_img1090',...
%     '210416_img1088','201107_img1073','201107_img1074','210616_img1091','210616_img1092','210909_img1098',...
%     '210921_img1901'};
% CS1 trained
sessions = {'201125_img1078','201125_img1079','210119_img1083','210416_img1087',...
'210427_img1090','210428_img1088'};
% CS1 test day after training with CS1+CS2 %!!!!!!!!!!!!!!!!!201209_img1079 the run is 001
% sessions = {'201211_img1078','201209_img1079','210203_img1083','210429_img1087',...
%     '210508_img1088','210507_img1090'};
% CS2 test day after training with CS1+CS2
% sessions = {'201210_img1078','201210_img1079','210205_img1083','210430_img1087',...
% '210509_img1088','210511_img1090'};
% CS1+CS2 test day after training with CS1+CS2
% sessions = {'201208_img1078','201209_img1079','210202_img1083','210428_img1087',...
%     '210507_img1088','210506_img1090'};

% control group- 1901 and 1092
% naive day
% sessions = {'210616_img1091','210616_img1092'};
% trained day 9
% sessions = {'210628_img1091','210628_img1092'};

behav_analysis_out = 'Z:\home\shuyang\behavior_analysis\RC\';
color_code = 'b';
%fig_name = 'across_session_RTCDF_CS1_naive';
fig_name = 'across_session_RTCDF_CS1_trained';
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
RTcum_fig = figure; 
RTcum_fig.Units = 'inches';
RTcum_fig.Position = [1 1 4 2.5];
cdfplot(RT);
vline(0,color_code);vline(700,'k');
xlim([-500 5000]);
ylabel('Fraction of RTs');
ylim([0 1]);
%title('CS1 naive day, across sessions');
title('CS1 trained day, across sessions');
%title('CS1+CS2 trained test w/ CS1, across sessions');
%title('CS1+CS2 trained test w/ CS2, across sessions');
%title('CS1+CS2 trained test w/ CS1+CS2, across sessions');
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
% CS1 naive
% sessions = {'201113_img1078','201114_img1079','210112_img1083','210113_img1082',...
% '210405_img1087','210412_img1090','210416_img1088'};
% CS1 trained
sessions = {'201125_img1078','201125_img1079','210119_img1083','210416_img1087',...
'210427_img1090','210428_img1088'};
% CS1 test day after training with CS1+CS2 %!!!!!!!!!!!!!!!!!201209_img1079 the run is 001
% sessions = {'201211_img1078','201209_img1079','210203_img1083','210429_img1087',...
%     '210508_img1088','210507_img1090'};
% CS2 test day after training with CS1+CS2
% sessions = {'201210_img1078','201210_img1079','210205_img1083','210430_img1087',...
% '210509_img1088','210511_img1090'};
% CS1+CS2 test day after training with CS1+CS2
% sessions = {'201208_img1078','201209_img1079','210202_img1083','210428_img1087',...
%     '210507_img1088','210506_img1090'};

% control group- 1901 and 1092
% naive day
% sessions = {'210616_img1091','210616_img1092'};
% trained day 9
% sessions = {'210628_img1091','210628_img1092'};

run = {'000','000','000','000','000','000','000'}; 
% run = {'000','001','000','000','000','000','000'}; %!!!!!!!!!!!!!!!!!201209_img1079 the run is 001
analysis_out    = 'Z:\home\shuyang\2P_analysis\';
color_code = 'b';
fig_name = 'across_session_targetalign_lickCDF_CS1_trained';
target_align_licks_all_session = [];
ntrials_total = 0;
lick_times = [];
for ii = 1:length(sessions)
    %load things
    lickAlign_output = load((fullfile([analysis_out sessions{ii} '\' sessions{ii} '_' run{ii} '_cueAlignLick.mat'])));    lickAlign_thisSession = lickAlign_output.lickCueAlign; % frame*trial
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
%title('CS1 naive day, across sessions');
title('CS1 trained day, across sessions');
%title('CS1+CS2 trained test w/ CS1, across sessions');
%title('CS1+CS2 trained test w/ CS2, across sessions');
%title('CS1+CS2 trained test w/ CS1+CS2, across sessions');
xlabel('Time from cue (ms)');
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'FontSize',7);
%xlim([-1.1 1.6]);ylim([0 3]);
box off;
path = 'Z:\home\shuyang\2P_analysis\across_sessions_RC\';
print(lickcum_fig,[path,fig_name],'-r600','-dpdf');
print(lickcum_fig,[path,fig_name],'-djpeg');
savefig([path fig_name]);

%% licking aligned to cue, histogram
clear;
% CS1 naive
% sessions = {'201113_img1078-1308','201114_img1079','210112_img1083','210405_img1087','210412_img1090',...
%     '210416_img1088','201107_img1073','201107_img1074','210616_img1091','210616_img1092','210909_img1098',...
%     '210921_img1901'};
% CS1 trained
% sessions = {'201125_img1078','201125_img1079','210119_img1083','210416_img1087',...
% '210427_img1090','210428_img1088','201111_img1074'};
% CS1 test day after training with CS1+CS2 %!!!!!!!!!!!!!!!!!201209_img1079 the run is 001
% sessions = {'201211_img1078-1251','201209_img1079-1325','210203_img1083','210429_img1087',...
%     '210508_img1088','210507_img1090'};
% CS2 test day after training with CS1+CS2
% sessions = {'201210_img1078-1134','201210_img1079-1252','210205_img1083-1140','210430_img1087',...
% '210509_img1088','210511_img1090'};
% CS1+CS2 test day after training with CS1+CS2
sessions = {'201208_img1078','201209_img1079-1251','210202_img1083','210428_img1087',...
    '210507_img1088','210506_img1090'};

% control group- 1901 and 1092
% naive day
% sessions = {'210616_img1091','210616_img1092'};
% trained day 9
% sessions = {'210628_img1091','210628_img1092'};

analysis_out    = 'Z:\home\shuyang\2P_analysis\';
behav_out  = 'Z:\home\shuyang\behavior_analysis\RC\hist_data_across_animals\';
color_code = 'r';
% fig_name = 'across_session_targetalign_lick_CS1_naive_hist';
% data_file = 'across_sessions_RC_data_CS1_naive.mat';
% fig_name = 'across_session_targetalign_lick_CS1_trained_hist';
% data_file = 'across_sessions_RC_data_CS1_trained.mat';
fig_name = 'across_session_targetalign_lick_CS1CS2_testCS1CS2_hist';
data_file = 'across_sessions_RC_data_CS1CS2_testCS1CS2.mat';
% fig_name = 'across_session_targetalign_lick_CS1CS2_testCS1_hist';
% data_file = 'across_sessions_RC_data_CS1CS2_testCS1.mat';
% fig_name = 'across_session_targetalign_lick_CS1CS2_testCS2_hist';
% data_file = 'across_sessions_RC_data_CS1CS2_testCS2.mat';

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
save(fullfile(['Z:\home\shuyang\2P_analysis\across_sessions_RC_HPfiltered\',data_file]),...
    'licks_hist_all_avg','licks_hist_all_sem','-append');


num_animals = size(licks_hist_all_session, 1);
x_axis = (([1:length(licks_hist_all_avg)]-181)*100 + 50)/1000;

lickhist_fig = figure; 
lickhist_fig.Units = 'centimeters';
lickhist_fig.Position = [1 1 18 8];
bar(x_axis, licks_hist_all_avg,'FaceColor',[0.6 0.6 0.6],'EdgeColor',[0.6 0.6 0.6]); hold on;
errorbar(x_axis, licks_hist_all_avg, licks_hist_all_sem,'Color',[0.5 0.5 0.5], 'LineStyle', 'none', 'CapSize', 1);
%title(['Pre-learning, CS1']);
%title(['Post-learning, CS1 trained']);
title('+ 8 days CS1+CS2 training, CS1+CS2');
%title('+ 8 days CS1+CS2 training, CS1');
%title('+ 8 days CS1+CS2 training, CS2');
xlabel('time relative to cue (s)');
ylabel('lick rate (Hz)');
ylim([0 10]);
vline(0.7, color_code);
vline(0, 'k')
xlim([-1 4]);
text(3.4, 9.4,['n = ', num2str(num_animals)],'FontSize', 20,'FontName','Arial');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',20,'FontName','Arial');
box off;
path = 'Z:\home\shuyang\2P_analysis\across_sessions_RC_HPfiltered\2021retreat\';
print(lickhist_fig,[path,fig_name],'-r600','-dpdf');
print(lickhist_fig,[path,fig_name],'-djpeg');
savefig([path fig_name]);


