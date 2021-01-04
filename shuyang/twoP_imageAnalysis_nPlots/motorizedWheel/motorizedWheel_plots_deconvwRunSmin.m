% analysis for motorized whell experiment, plot both df/f and FR.
% average df/f during slow speed vs. high speed vs. stationary
% same thins for spikes
%dfOvF_fast: trial*frame*cell (dfOvF of all cells and frames during all fast running trials, first second of each trial is thrown away)
%dfOvF_fast_cells: dfOvF average across trials and frames, each cell has one value
%FR_fast_cells: firing rate of each cell average across trials
%runtrig ave and runoff ave for df/f and spikes


%% assign document paths and experimental sessions
clear;
%sessions = {'200116_img1041','200217_img1061','200225_img1049',...
%    '200319_img1064','200319_img1064_2'};
%days = {'1041-200116_1','1061-200217_1','1049-200225_1',...
%    '1064-200319_1','1064-200319_2'};
sessions = '200319_img1064_2'; 
image_analysis_base = 'Z:\Analysis\motorizedWheel_Analysis\running\imaging_analysis\'; 
image_analysis_dest = [image_analysis_base, sessions, '\'];

% behavior analysis results 
days = '1064-200319_2';
behav_dest = ['Z:\Analysis\motorizedWheel_Analysis\running\behavioral_analysis\' days '\'];
color_code = {'b','r','k','c'};

%% load data
behav_output = load([behav_dest days '_behavAnalysis.mat']);
run_fast_mat = behav_output.run_fast_mat;
run_slow_mat = behav_output.run_slow_mat;
stationary_mat = behav_output.staytionary_mat;
runtrig_slow_mat = behav_output.runtrig_slow_mat;
runtrig_fast_mat = behav_output.runtrig_fast_mat;
runoff_slow_mat = behav_output.runoff_slow_mat;
runoff_fast_mat = behav_output.runoff_fast_mat;
spd_decrease_mat = behav_output.spd_decrease_mat;
spd_increase_mat = behav_output.spd_increase_mat;
runtrig_mat = behav_output.runtrig_mat;
runoff_mat = behav_output.runoff_mat;
run_fast_vec = behav_output.run_fast_vec;
run_slow_vec = behav_output.run_slow_vec;
stay_vec = behav_output.stay_vec;

spk_deconv_output = load([image_analysis_dest,'deconv_wfastRunSmin\' sessions, '_spk_deconvolve_staynrun_seperate.mat']);
dfOvF = load([image_analysis_dest sessions '_dfOvF.mat']);
dfOvF = dfOvF.dfOvF_btm_cl;

%% df/F scatter plot during slow, fast, and stationary
%{
% dfOvF_fast = zeros(size(run_fast_mat,1),size(run_fast_mat,2),size(dfOvF,2)); %trial*frame*cell
% dfOvF_slow = zeros(size(run_slow_mat,1),size(run_slow_mat,2),size(dfOvF,2));
% dfOvF_stationary = zeros(size(stationary_mat,1),size(stationary_mat,2),size(dfOvF,2));
% for c = 1:size(dfOvF,2) %for each cell
%     for f = 1:size(run_fast_mat,1)
%         dfOvF_fast(f,:,c) = dfOvF(run_fast_mat(f,:),c);
%     end
%     for s = 1:size(run_slow_mat,1)
%         dfOvF_slow(s,:,c) = dfOvF(run_slow_mat(s,:),c);
%     end
%     for t = 1:size(stationary_mat,1)
%         dfOvF_stay(t,:,c) = dfOvF(stationary_mat(t,:),c);
%     end
% end
% 
% % save df/F matrix: trial*frame*cell and frame*cell
% save([image_analysis_dest sessions '_dfOvF.mat'],'dfOvF_fast','dfOvF_slow','dfOvF_stationary','-append');

dfOvF_fast = dfOvF(run_fast_vec,:); % using the vector should be the same thing as using the matrix above
dfOvF_slow = dfOvF(run_slow_vec,:);
dfOvF_stay = dfOvF(stay_vec,:);

dfOvF_fast_cells = mean(dfOvF_fast,1); % average across frames
dfOvF_fast_all = mean(dfOvF_fast_cells);
dfOvF_fast_ste = std(dfOvF_fast_cells)/sqrt(size(dfOvF,2));

dfOvF_slow_cells = mean(dfOvF_slow,1);
dfOvF_slow_all = mean(dfOvF_slow_cells);
dfOvF_slow_ste = std(dfOvF_slow_cells)/sqrt(size(dfOvF,2));

dfOvF_stay_cells = mean(dfOvF_stay,1);
dfOvF_stay_all = mean(dfOvF_stay_cells);
dfOvF_stay_ste = std(dfOvF_stay_cells)/sqrt(size(dfOvF,2));

save([image_analysis_dest sessions '_dfOvF.mat'],'dfOvF_fast_cells','dfOvF_slow_cells','dfOvF_stay_cells','-append');

dfOvF_fig = figure;
x = [1,2,3];
dfOvF_plot = [dfOvF_fast_all,dfOvF_slow_all,dfOvF_stay_all];
dfOvF_ste_plot = [dfOvF_fast_ste,dfOvF_slow_ste,dfOvF_stay_ste];
% x_plot = repmat(x,size(dfOvF,1),1);
errorbar(x,dfOvF_plot,dfOvF_ste_plot,'.','LineStyle','none','MarkerSize',20);

xlim([0.5 3.5]);
x1= [1,2,3];
set(gca,'XTick',x1,'XTicklabel',{'fast running','slow running','stationary'});
ylabel('df/f');
title(sessions);
savefig([image_analysis_dest sessions '_dfOvF_scatter.fig']);
%}

%% FR scatter plot
spk_logic = spk_deconv_output.spk_logic;

FR_fast_cells = zeros(1,size(spk_logic,2));
FR_slow_cells = zeros(1,size(spk_logic,2));
FR_stay_cells = zeros(1,size(spk_logic,2));

for c = 1:size(spk_logic,2) %for each cell
    FR_fast_cells(c) = sum(spk_logic(run_fast_vec,c))/(length(run_fast_vec)/30); %using the vector this way should be the same as using the matirx
    FR_slow_cells(c) = sum(spk_logic(run_slow_vec,c))/(length(run_fast_vec)/30);
    FR_stay_cells(c) = sum(spk_logic(stay_vec,c))/(length(stay_vec)/30);
%     for f = 1:size(run_fast_mat,1)
%         %spk_fast(f,:,c) = spk_logic(run_fast_mat(f,:),c);
%         %calculate firing rate: for each cell each trial
%         FR_fast(f,c) = sum(spk_fast(f,:,c)==1)/(size(run_fast_mat,2)/30);
%     end
%     for s = 1:size(run_slow_mat,1)
%         spk_slow(s,:,c) = spk_logic(run_slow_mat(s,:),c);
%         FR_slow(s,c) = sum(spk_slow(s,:,c)==1)/(size(run_slow_mat,2)/30);
%     end
%     for t = 1:size(stationary_mat,1)
%         spk_stationary(t,:,c) = spk_logic(stationary_mat(t,:),c);
%         FR_stay(t,c) = sum(spk_stationary(t,:,c)==1)/(size(stationary_mat,2)/30);
%         
%     end
end
% FR_fast_cells = mean(FR_fast);
% FR_slow_cells = mean(FR_slow);
% FR_stay_cells = mean(FR_stay);

FR_fast_all = mean(FR_fast_cells);
FR_slow_all = mean(FR_slow_cells);
FR_stay_all = mean(FR_stay_cells);

save([image_analysis_dest,'deconv_wfastRunSmin\' sessions, '_spk_deconvolve_staynrun_seperate.mat'],...
    'FR_fast_cells','FR_slow_cells','FR_stay_cells','-append');

FR_fast_ste = std(FR_fast_cells)/sqrt(size(spk_logic,2));
FR_slow_ste = std(FR_slow_cells)/sqrt(size(spk_logic,2));
FR_stay_ste = std(FR_stay_cells)/sqrt(size(spk_logic,2));

FR_fig = figure;
x = [1,2,3];
FR_plot = [FR_fast_all,FR_slow_all,FR_stay_all];
FR_ste_plot = [FR_fast_ste,FR_slow_ste,FR_stay_ste];
% x_plot = repmat(x,size(dfOvF,1),1);
errorbar(x,FR_plot,FR_ste_plot,'.','LineStyle','none','MarkerSize',20);
xlim([0.5 3.5]);
x1= [1,2,3];
set(gca,'XTick',x1,'XTicklabel',{'fast running','slow running','stationary'});
ylabel('firing rate');
title(sessions);
savefig([image_analysis_dest,'deconv_wfastRunSmin\' sessions '_FR_scatter.fig']);

%% df/f runtrig ave and runoffset ave and speed transit
%{
dfOvF_runoff_fast = zeros(size(runoff_fast_mat,1),size(runoff_fast_mat,2),size(dfOvF,2)); %trial*frame*cell
dfOvF_runoff_slow = zeros(size(runoff_slow_mat,1),size(runoff_slow_mat,2),size(dfOvF,2));
dfOvF_runtrig_fast = zeros(size(runtrig_fast_mat,1),size(runtrig_fast_mat,2),size(dfOvF,2));
dfOvF_runtrig_slow = zeros(size(runtrig_slow_mat,1),size(runtrig_slow_mat,2),size(dfOvF,2));
dfOvF_spd_decrease = zeros(size(spd_decrease_mat,1),size(spd_decrease_mat,2),size(dfOvF,2));
dfOvF_spd_increase = zeros(size(spd_increase_mat,1),size(spd_increase_mat,2),size(dfOvF,2));

for c = 1:size(dfOvF,2) %for each cell
    for f = 1:size(runoff_fast_mat,1)
        dfOvF_runoff_fast(f,:,c) = dfOvF(runoff_fast_mat(f,:),c);
    end
    for g = 1:size(runoff_slow_mat,1)
        dfOvF_runoff_slow(g,:,c) = dfOvF(runoff_slow_mat(g,:),c);
    end
    for h = 1:size(runtrig_fast_mat,1)
        dfOvF_runtrig_fast(h,:,c) = dfOvF(runtrig_fast_mat(h,:),c);
    end
    for i = 1:size(runtrig_slow_mat,1)
        dfOvF_runtrig_slow(i,:,c) = dfOvF(runtrig_slow_mat(i,:),c);
    end
    for j = 1:size(spd_decrease_mat,1)
        dfOvF_spd_decrease(j,:,c) = dfOvF(spd_decrease_mat(j,:),c);
    end
    for k = 1:size(spd_increase_mat,1)
        dfOvF_spd_increase(k,:,c) = dfOvF(spd_increase_mat(k,:),c);
    end
end

%average across trials
dfOvF_runoff_fast_cells = squeeze(mean(dfOvF_runoff_fast,1)); %frame*cell
dfOvF_runoff_slow_cells = squeeze(mean(dfOvF_runoff_slow,1));
dfOvF_runtrig_fast_cells = squeeze(mean(dfOvF_runtrig_fast,1));
dfOvF_runtrig_slow_cells = squeeze(mean(dfOvF_runtrig_slow,1));
dfOvF_spd_decrease_cells = squeeze(mean(dfOvF_spd_decrease,1));
dfOvF_spd_increase_cells = squeeze(mean(dfOvF_spd_increase,1));

save([image_analysis_dest sessions '_dfOvF.mat'],'dfOvF_runoff_fast_cells',...
    'dfOvF_runoff_slow_cells','dfOvF_runtrig_fast_cells','dfOvF_runtrig_slow_cells',...
    'dfOvF_spd_decrease_cells','dfOvF_spd_increase_cells','-append');

%average across cells
dfOvF_runoff_fast_session = mean(dfOvF_runoff_fast_cells,2);
dfOvF_runoff_slow_session = mean(dfOvF_runoff_slow_cells,2);
dfOvF_runtrig_fast_session = mean(dfOvF_runtrig_fast_cells,2);
dfOvF_runtrig_slow_session = mean(dfOvF_runtrig_slow_cells,2);
dfOvF_spd_decrease_session = mean(dfOvF_spd_decrease_cells,2);
dfOvF_spd_increase_session = mean(dfOvF_spd_increase_cells,2);

ncells = size(dfOvF,2);
ste_runoff_fast_cells = std(dfOvF_runoff_fast_cells,0,2)/sqrt(ncells);
ste_runoff_slow_cells = std(dfOvF_runoff_slow_cells,0,2)/sqrt(ncells);
ste_runtrig_fast_cells = std(dfOvF_runtrig_fast_cells,0,2)/sqrt(ncells);
ste_runtrig_slow_cells = std(dfOvF_runtrig_slow_cells,0,2)/sqrt(ncells);
ste_spd_decrease_cells = std(dfOvF_spd_decrease_cells,0,2)/sqrt(ncells);
ste_spd_increase_cells = std(dfOvF_spd_increase_cells,0,2)/sqrt(ncells);

dfOvF_trig = figure;
x = (1:length(dfOvF_runoff_fast_session))/30 - 1;
subplot(2,3,1);
shadedErrorBar(x,dfOvF_runtrig_slow_session,ste_runtrig_slow_cells,{'color',[0.1922 0.6392 0.3294]});
ylabel('df/f'); ylim([0 1]); xlim([-1 1.5]);
vline(0,'k');
title('stay --> slow run');
%xlim([0,1]);
subplot(2,3,4);
shadedErrorBar(x,dfOvF_runoff_slow_session,ste_runoff_slow_cells,{'color',[0.1922 0.6392 0.3294]});
ylabel('df/f'); ylim([0 1]);xlim([-1 1.5]);
vline(0,'k');
title('slow run --> stay');
subplot(2,3,2);
shadedErrorBar(x,dfOvF_runtrig_fast_session,ste_runtrig_fast_cells,{'color',[0.1922 0.6392 0.3294]});
ylim([0 1]); vline(0,'k'); xlabel('time from behavioral state change (s)');xlim([-1 1.5]);
title('stay --> fast run');
subplot(2,3,5);
shadedErrorBar(x,dfOvF_runoff_fast_session,ste_runoff_fast_cells,{'color',[0.1922 0.6392 0.3294]});
ylim([0 1]); vline(0,'k'); xlabel('time from behavioral state change (s)');xlim([-1 1.5]);
title('fast run --> stay');
subplot(2,3,3);
shadedErrorBar(x,dfOvF_spd_decrease_session,ste_spd_decrease_cells,{'color',[0.1922 0.6392 0.3294]});
ylim([0 1]); xlim([-1 1.5]);vline(0,'k');
title('speed decrease');
subplot(2,3,6);
shadedErrorBar(x,dfOvF_spd_increase_session,ste_spd_increase_cells,{'color',[0.1922 0.6392 0.3294]});
xlim([-1 1.5]); ylim([0 1]); vline(0,'k');
title('speed increase');
supertitle(['dfOvF' sessions]);

savefig([image_analysis_dest sessions '_dfOvF_state_transit.fig']);
%}

%% FR run trig ave and runoff ave
spk_logic = spk_deconv_output.spk_logic;
%spk_logic = spk_deconv_output.spk_logic;
spk_runoff_fast = zeros(size(runoff_fast_mat,1),size(runoff_fast_mat,2),size(spk_logic,2)); %trial*frame*cell
spk_runoff_slow = zeros(size(runoff_slow_mat,1),size(runoff_slow_mat,2),size(spk_logic,2));
spk_runtrig_fast = zeros(size(runtrig_fast_mat,1),size(runtrig_fast_mat,2),size(spk_logic,2));
spk_runtrig_slow = zeros(size(runtrig_slow_mat,1),size(runtrig_slow_mat,2),size(spk_logic,2));
spk_spd_decrease = zeros(size(spd_decrease_mat,1),size(spd_decrease_mat,2),size(spk_logic,2));
spk_spd_increase = zeros(size(spd_increase_mat,1),size(spd_increase_mat,2),size(spk_logic,2));

for c = 1:size(spk_logic,2) %for each cell
    for f = 1:size(runoff_fast_mat,1)
        spk_runoff_fast(f,:,c) = spk_logic(runoff_fast_mat(f,:),c);
    end
    for g = 1:size(runoff_slow_mat,1)
        spk_runoff_slow(g,:,c) = spk_logic(runoff_slow_mat(g,:),c);
    end
    for h = 1:size(runtrig_fast_mat,1)
        spk_runtrig_fast(h,:,c) = spk_logic(runtrig_fast_mat(h,:),c);
    end
    for i = 1:size(runtrig_slow_mat,1)
        spk_runtrig_slow(i,:,c) = spk_logic(runtrig_slow_mat(i,:),c);
    end
    for j = 1:size(spd_decrease_mat,1)
        spk_spd_decrease(j,:,c) = spk_logic(spd_decrease_mat(j,:),c);
    end
    for k = 1:size(spd_increase_mat,1)
        spk_spd_increase(k,:,c) = spk_logic(spd_increase_mat(k,:),c);
    end
end

%average across trials and get average firing rate for each cell
FR_runoff_fast_cells = squeeze(mean(spk_runoff_fast,1))*30;
FR_runoff_slow_cells = squeeze(mean(spk_runoff_slow,1))*30;
FR_runtrig_fast_cells = squeeze(mean(spk_runtrig_fast,1))*30;
FR_runtrig_slow_cells = squeeze(mean(spk_runtrig_slow,1))*30;
FR_spd_decrease_cells = squeeze(mean(spk_spd_decrease,1))*30;
FR_spd_increase_cells = squeeze(mean(spk_spd_increase,1))*30;

save([image_analysis_dest,'deconv_wfastRunSmin\' sessions, '_spk_deconvolve_staynrun_seperate.mat'],...
    'FR_runoff_fast_cells','FR_runoff_slow_cells','FR_runtrig_fast_cells',...
    'FR_runtrig_slow_cells','FR_spd_decrease_cells','FR_spd_increase_cells','-append');

%average across cells
FR_runoff_fast_session = mean(FR_runoff_fast_cells,2);
FR_runoff_slow_session = mean(FR_runoff_slow_cells,2);
FR_runtrig_fast_session = mean(FR_runtrig_fast_cells,2);
FR_runtrig_slow_session = mean(FR_runtrig_slow_cells,2);
FR_spd_decrease_session = mean(FR_spd_decrease_cells,2);
FR_spd_increase_session = mean(FR_spd_increase_cells,2);

ncells = size(spk_logic,2);
steFR_runoff_fast_cells = std(FR_runoff_fast_cells,0,2)/sqrt(ncells);
steFR_runoff_slow_cells = std(FR_runoff_slow_cells,0,2)/sqrt(ncells);
steFR_runtrig_fast_cells = std(FR_runtrig_fast_cells,0,2)/sqrt(ncells);
steFR_runtrig_slow_cells = std(FR_runtrig_slow_cells,0,2)/sqrt(ncells);
steFR_spd_decrease_cells = std(FR_spd_decrease_cells,0,2)/sqrt(ncells);
steFR_spd_increase_cells = std(FR_spd_increase_cells,0,2)/sqrt(ncells);

FR_trig = figure;
x = (1:length(FR_runoff_fast_session))/30;
subplot(2,3,1);
shadedErrorBar(x,FR_runtrig_slow_session,steFR_runtrig_slow_cells,{'color',[0.1373 0.5451 0.2706]});
ylabel('firing rate'); 
xlim([0 2.5]); ylim([-0.05 4.5]);vline(1,'k');
title('stay-slow speed');
%xlim([0,1]);
subplot(2,3,4);
shadedErrorBar(x,FR_runoff_slow_session,steFR_runoff_slow_cells,{'color',[0.1373 0.5451 0.2706]});
xlim([0 2.5]); ylim([-0.05 4.5]);ylabel('firing rate'); 
vline(1,'k');
title('slow speed-stay');
subplot(2,3,2);
shadedErrorBar(x,FR_runtrig_fast_session,steFR_runtrig_fast_cells,{'color',[0.1373 0.5451 0.2706]});
xlabel('time(s)');xlim([0 2.5]); ylim([-0.05 4.5]); vline(1,'k');
title('stay-fast speed');
subplot(2,3,5);
shadedErrorBar(x,FR_runoff_fast_session,steFR_runoff_fast_cells,{'color',[0.1373 0.5451 0.2706]});
xlabel('time(s)');xlim([0 2.5]); ylim([-0.05 4.5]); vline(1,'k');
title('fast speed-stay');
subplot(2,3,3);
shadedErrorBar(x,FR_spd_decrease_session,steFR_spd_decrease_cells,{'color',[0.1373 0.5451 0.2706]});
xlim([0 2.5]); ylim([-0.05 4.5]); vline(1,'k');
title('speed decrease');
subplot(2,3,6);
shadedErrorBar(x,FR_spd_increase_session,steFR_spd_increase_cells,{'color',[0.1373 0.5451 0.2706]});
xlim([0 2.5]); ylim([-0.05 4.5]); vline(1,'k');
title('speed increase');
supertitle(['Firing rate' sessions]);

savefig([image_analysis_dest,'deconv_wfastRunSmin\' sessions '_FR_state_transit.fig']);

