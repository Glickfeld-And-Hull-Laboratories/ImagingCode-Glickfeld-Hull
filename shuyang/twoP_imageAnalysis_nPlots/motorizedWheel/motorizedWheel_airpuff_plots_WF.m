% average df/f during slow speed vs. high speed vs. stationary
% same thins for spikes
%dfOvF_fast: trial*frame*cell (dfOvF of all cells and frames during all fast running trials, first second of each trial is thrown away)
%dfOvF_fast_cells: dfOvF average across trials and frames, each cell has one value
%FR_fast_cells: firing rate of each cell average across trials
%runtrig ave and runoff ave for df/f and spikes


%% assign document paths and experimental sessions
clear;
sessions = '220126_img1911_1mMGabazine1uMAlexa1.5ul_1'; 
%image_analysis_base = 'Z:\home\shuyang\Analysis\motorizedWheel_Analysis\WF\imaging_analysis\'; 
image_analysis_base = 'Z:\home\shuyang\Analysis\motorizedWheel_Analysis\WF_drug\imaging_analysis\'; 
image_analysis_dest = [image_analysis_base, sessions, '\'];

% behavior analysis results 
days = '1911-220126-1629_1'; 
%behav_dest = ['Z:\home\shuyang\Analysis\motorizedWheel_Analysis\WF\behavioral_analysis\' days '\'];
behav_dest = ['Z:\home\shuyang\Analysis\motorizedWheel_Analysis\WF_drug\behavioral_analysis\' days '\'];
color_code = {'b','r','k','c'};

%% load data
behav_output = load([behav_dest days '_behavAnalysis.mat']);
airpuff_fast_mat = behav_output.airpuff_fast_mat;
airpuff_stay_mat = behav_output.airpuff_stay_mat;
airpuff_slow_mat = behav_output.airpuff_slow_mat;
runtrig_slow_mat = behav_output.runtrig_slow_mat;
runtrig_fast_mat = behav_output.runtrig_fast_mat;
runoff_slow_mat = behav_output.runoff_slow_mat;
runoff_fast_mat = behav_output.runoff_fast_mat;
spd_decrease_mat = behav_output.spd_decrease_mat;
spd_increase_mat = behav_output.spd_increase_mat;
runslow_nopuff = behav_output.runslow_nopuff;
runfast_nopuff = behav_output.runfast_nopuff;
stay_nopuff = behav_output.stay_nopuff;

dfOvF = load([image_analysis_dest sessions '_dfOvF.mat']);
%dfOvF = load([image_analysis_dest sessions '_dfOvF_ctrlbase.mat']);
dfOvF = dfOvF.dfOvF_btmbase; % ROI*frames
%dfOvF = dfOvF.dfOvF_ctrlbase;


%% calculate the amount of suppression during running
dfOvF_fast = dfOvF(:,runfast_nopuff); % using the vector should be the same thing as using the matrix above
dfOvF_slow = dfOvF(:,runslow_nopuff);
dfOvF_stay = dfOvF(:,stay_nopuff);

dfOvF_fast_ROIs = mean(dfOvF_fast,2); % average across frames
%dfOvF_fast_all = mean(dfOvF_fast_ROIs);
%dfOvF_fast_ste = std(dfOvF_fast_ROIs)/sqrt(size(dfOvF,1));

dfOvF_slow_ROIs = mean(dfOvF_slow,2);
%dfOvF_slow_all = mean(dfOvF_slow_ROIs);
%dfOvF_slow_ste = std(dfOvF_slow_ROIs)/sqrt(size(dfOvF,1));

dfOvF_stay_ROIs = mean(dfOvF_stay,2);
%dfOvF_stay_all = mean(dfOvF_stay_ROIs);
%dfOvF_stay_ste = std(dfOvF_stay_ROIs)/sqrt(size(dfOvF,1));

suppress_fastrun_ROIs = dfOvF_stay_ROIs - dfOvF_fast_ROIs;
suppress_slowrun_ROIs = dfOvF_stay_ROIs - dfOvF_slow_ROIs;

save([image_analysis_dest sessions '_dfOvF.mat'],'dfOvF_fast_ROIs','dfOvF_slow_ROIs','dfOvF_stay_ROIs',...
    'suppress_fastrun_ROIs','suppress_slowrun_ROIs','-append');

%% df/f runtrig ave and runoffset ave and speed transit
dfOvF_runoff_fast = zeros(size(runoff_fast_mat,1),size(runoff_fast_mat,2),size(dfOvF,1)); %trial*frame*ROI
dfOvF_runoff_slow = zeros(size(runoff_slow_mat,1),size(runoff_slow_mat,2),size(dfOvF,1));
dfOvF_runtrig_fast = zeros(size(runtrig_fast_mat,1),size(runtrig_fast_mat,2),size(dfOvF,1));
dfOvF_runtrig_slow = zeros(size(runtrig_slow_mat,1),size(runtrig_slow_mat,2),size(dfOvF,1));
dfOvF_spd_decrease = zeros(size(spd_decrease_mat,1),size(spd_decrease_mat,2),size(dfOvF,1));
dfOvF_spd_increase = zeros(size(spd_increase_mat,1),size(spd_increase_mat,2),size(dfOvF,1));
for r = 1:size(dfOvF,1)
    for f = 1:size(runoff_fast_mat,1)
        dfOvF_runoff_fast(f,:,r) = dfOvF(r,runoff_fast_mat(f,:));
    end
    for g = 1:size(runoff_slow_mat,1)
        dfOvF_runoff_slow(g,:,r) = dfOvF(r,runoff_slow_mat(g,:));
    end
    for h = 1:size(runtrig_fast_mat,1)
        dfOvF_runtrig_fast(h,:,r) = dfOvF(r,runtrig_fast_mat(h,:));
    end
    for i = 1:size(runtrig_slow_mat,1)
        dfOvF_runtrig_slow(i,:,r) = dfOvF(r,runtrig_slow_mat(i,:));
    end
    for j = 1:size(spd_decrease_mat,1)
        dfOvF_spd_decrease(j,:,r) = dfOvF(r,spd_decrease_mat(j,:));
    end
    for k = 1:size(spd_increase_mat,1)
        dfOvF_spd_increase(k,:,r) = dfOvF(r,spd_increase_mat(k,:));
    end
end

%average across trials
dfOvF_runoff_fast_ROIs = squeeze(mean(dfOvF_runoff_fast,1));
dfOvF_runoff_slow_ROIs = squeeze(mean(dfOvF_runoff_slow,1));
dfOvF_runtrig_fast_ROIs = squeeze(mean(dfOvF_runtrig_fast,1));
dfOvF_runtrig_slow_ROIs = squeeze(mean(dfOvF_runtrig_slow,1));
dfOvF_spd_decrease_ROIs = squeeze(mean(dfOvF_spd_decrease,1));
dfOvF_spd_increase_ROIs = squeeze(mean(dfOvF_spd_increase,1));

ste_runoff_fast = squeeze(std(dfOvF_runoff_fast,0,1)/sqrt(size(dfOvF_runoff_fast,1)));
ste_runoff_slow = squeeze(std(dfOvF_runoff_slow,0,1)/sqrt(size(dfOvF_runoff_slow,1)));
ste_runtrig_fast = squeeze(std(dfOvF_runtrig_fast,0,1)/sqrt(size(dfOvF_runtrig_fast,1)));
ste_runtrig_slow = squeeze(std(dfOvF_runtrig_slow,0,1)/sqrt(size(dfOvF_runtrig_slow,1)));
ste_spd_decrease = squeeze(std(dfOvF_spd_decrease,0,1)/sqrt(size(dfOvF_spd_decrease,1)));
ste_spd_increase = squeeze(std(dfOvF_spd_increase,0,1)/sqrt(size(dfOvF_spd_increase,1)));
% if size(dfOvF,1)>1 % if more than 1 ROI
%     dfOvF_runoff_fast_ave = squeeze(mean(dfOvF_runoff_fast,2)); % average across ROIs
%     dfOvF_runoff_slow_ave = squeeze(mean(dfOvF_runoff_slow,2));
%     dfOvF_runtrig_fast_ave = squeeze(mean(dfOvF_runtrig_fast,2));
%     dfOvF_runtrig_slow_ave = squeeze(mean(dfOvF_runtrig_slow,2));
%     dfOvF_spd_decrease_ave = squeeze(mean(dfOvF_spd_decrease,2));
%     dfOvF_spd_increase_ave = squeeze(mean(dfOvF_spd_increase,2));
% else
%     dfOvF_runoff_fast_ave = dfOvF_runoff_fast_ROIs; 
%     dfOvF_runoff_slow_ave = dfOvF_runoff_slow_ROIs;
%     dfOvF_runtrig_fast_ave = dfOvF_runtrig_fast_ROIs;
%     dfOvF_runtrig_slow_ave = dfOvF_runtrig_slow_ROIs;
%     dfOvF_spd_decrease_ave = dfOvF_spd_decrease_ROIs;
%     dfOvF_spd_increase_ave = dfOvF_spd_increase_ROIs;
% end

%save([image_analysis_dest sessions '_dfOvF_ctrlbase.mat'],'dfOvF_runoff_fast_ROIs',...
save([image_analysis_dest sessions '_dfOvF.mat'],'dfOvF_runoff_fast_ROIs',...
    'dfOvF_runoff_slow_ROIs','dfOvF_runtrig_fast_ROIs','dfOvF_runtrig_slow_ROIs',...
    'dfOvF_spd_decrease_ROIs','dfOvF_spd_increase_ROIs','-append');

colorcode = {[0.8941 0.1020 0.1098],[0.2157 0.4941 0.7216],[0.3020 0.6863 0.2902],[0.5961 0.3059 0.6392],[1.0000 0.4980 0]};
dfOvF_trig = figure;
x = (1:size(dfOvF_runoff_fast_ROIs,1))/10 - 1; %secondsbefore in behavioral analysis = 1
subplot(2,3,1);
for r = 1:size(dfOvF,1)
shadedErrorBar(x,dfOvF_runtrig_slow_ROIs(:,r),ste_runtrig_slow(:,r),{'color',colorcode{r}}); hold on;
end
hold off;
ylabel('df/f'); xlim([-1 3]); ylim([0 0.5]); 
vline(0,'k');
title('stay-slow speed');box off;
subplot(2,3,4);
for r = 1:size(dfOvF,1)
shadedErrorBar(x,dfOvF_runoff_slow_ROIs(:,r),ste_runoff_slow(:,r),{'color',colorcode{r}});hold on;
end
hold off;
ylabel('df/f'); xlim([-1 3]);ylim([0 0.5]); 
vline(0,'k');
title('slow speed-stay');box off;
subplot(2,3,2);
for r = 1:size(dfOvF,1)
shadedErrorBar(x,dfOvF_runtrig_fast_ROIs(:,r),ste_runtrig_fast(:,r),{'color',colorcode{r}});hold on;
end
hold off;
xlim([-1 3]); ylim([0 0.5]); vline(0,'k'); xlabel('time(s)');
title('stay-fast speed'); box off;
subplot(2,3,5);
for r = 1:size(dfOvF,1)
shadedErrorBar(x,dfOvF_runoff_fast_ROIs(:,r),ste_runoff_fast(:,r),{'color',colorcode{r}});hold on;
end
hold off;
xlim([-1 3]); ylim([0 0.5]); vline(0,'k'); xlabel('time(s)');
title('fast speed-stay');box off;
subplot(2,3,3);
for r = 1:size(dfOvF,1)
shadedErrorBar(x,dfOvF_spd_decrease_ROIs(:,r),ste_spd_decrease(:,r),{'color',colorcode{r}});hold on;
end
hold off;
xlim([-1 3]);ylim([0 0.5]);  %xlim([0 2.5]);
vline(0,'k');
title('speed decrease');box off;
subplot(2,3,6);
for r = 1:size(dfOvF,1)
shadedErrorBar(x,dfOvF_spd_increase_ROIs(:,r),ste_spd_increase(:,r),{'color',colorcode{r}});hold on;
end
hold off;
xlim([-1 3]);ylim([0 0.5]); 
vline(0,'k');
title('speed increase');box off;
sgtitle(['dfOvF' sessions]);
savefig([image_analysis_dest sessions '_dfOvF_state_transit.fig']);
%savefig([image_analysis_dest sessions '_dfOvFctrl_state_transit.fig']);

%% df/f airpuff triggered average in different behavioral states
dfOvF_airpuff_fast = zeros(size(airpuff_fast_mat,1),size(airpuff_fast_mat,2),size(dfOvF,1)); %trial*frame*ROI
dfOvF_airpuff_slow = zeros(size(airpuff_slow_mat,1),size(airpuff_slow_mat,2),size(dfOvF,1));
dfOvF_airpuff_stay = zeros(size(airpuff_stay_mat,1),size(airpuff_stay_mat,2),size(dfOvF,1));

for r = 1:size(dfOvF,1)
    for f = 1:size(airpuff_fast_mat,1)
        dfOvF_airpuff_fast(f,:,r) = dfOvF(r,airpuff_fast_mat(f,:));
    end
    for g = 1:size(airpuff_slow_mat,1)
        dfOvF_airpuff_slow(g,:,r) = dfOvF(r,airpuff_slow_mat(g,:));
    end
    for h = 1:size(airpuff_stay_mat,1)
        dfOvF_airpuff_stay(h,:,r) = dfOvF(r,airpuff_stay_mat(h,:));
    end
end

%average across trials
dfOvF_airpuff_fast_ROIs = squeeze(mean(dfOvF_airpuff_fast,1)); %frame*cell
dfOvF_airpuff_slow_ROIs = squeeze(mean(dfOvF_airpuff_slow,1));
dfOvF_airpuff_stay_ROIs = squeeze(mean(dfOvF_airpuff_stay,1));

%combine matrix of 2 running speeds
dfOvF_airpuff_run = cat(1,dfOvF_airpuff_fast,dfOvF_airpuff_slow);

%save([image_analysis_dest sessions '_dfOvF_ctrlbase.mat'],'dfOvF_airpuff_fast_ROIs',...
save([image_analysis_dest sessions '_dfOvF.mat'],'dfOvF_airpuff_fast_ROIs',...
    'dfOvF_airpuff_slow_ROIs','dfOvF_airpuff_stay_ROIs',...
    'dfOvF_airpuff_run','dfOvF_airpuff_stay','-append');

ste_airpuff_fast = squeeze(std(dfOvF_airpuff_fast,0,1)/sqrt(size(dfOvF_airpuff_fast,1)));
ste_airpuff_slow = squeeze(std(dfOvF_airpuff_slow,0,1)/sqrt(size(dfOvF_airpuff_slow,1)));
ste_airpuff_stay = squeeze(std(dfOvF_airpuff_stay,0,1)/sqrt(size(dfOvF_airpuff_stay,1)));

colorcode = {[0.8941 0.1020 0.1098],[0.2157 0.4941 0.7216],[0.3020 0.6863 0.2902],[0.5961 0.3059 0.6392],[1.0000 0.4980 0]};
dfOvFairpuff_trig = figure;
x = (1:length(dfOvF_airpuff_fast_ROIs))/10 - 1; %seconds_around_air = 1 in behavioral analysis
subplot(1,3,1);
for r = 1:size(dfOvF,1)
    shadedErrorBar(x,dfOvF_airpuff_fast_ROIs(:,r),ste_airpuff_fast(:,r),{'color',colorcode{r}}); hold on;
end
hold off;
ylabel('df/f'); ylim([0 0.5]); xlim([-1 1]);
vline(0,'k');
title('fast speed');box off;
%xlim([0,1]);
subplot(1,3,2);
for r = 1:size(dfOvF,1)
    shadedErrorBar(x,dfOvF_airpuff_slow_ROIs(:,r),ste_airpuff_slow(:,r),{'color',colorcode{r}}); hold on;
end
hold off;
ylim([0 0.5]);xlim([-1 1]);
xlabel('time(s)');
vline(0,'k');
title('slow speed');box off;
subplot(1,3,3);
for r = 1:size(dfOvF,1)
    shadedErrorBar(x,dfOvF_airpuff_stay_ROIs(:,r),ste_airpuff_stay(:,r),{'color',colorcode{r}}); hold on;
end
hold off;
ylim([0 0.5]); 
vline(0,'k'); xlim([-1 1]);
title('stationary');box off;
sgtitle(['tactile stimulus' sessions]);
savefig([image_analysis_dest sessions '_dfOvF_airpuff.fig']);
%savefig([image_analysis_dest sessions '_dfOvFctrl_airpuff.fig']);

