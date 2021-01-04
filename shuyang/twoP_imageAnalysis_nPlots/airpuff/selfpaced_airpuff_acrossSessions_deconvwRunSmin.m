% airpuff trigger response during running and stationary, acorss all
% self-paced sessions

%% assign document paths and experimental sessions
clear;
sessions = {'191114_img1040','191115_img1039','191115_img1041','191115_img1042'};%,'200316_img1064_airpuff_2'};
image_analysis_base  = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
% behavior analysis results
days = {'1040-191114_1','1039-191115_1','1041-191115_1','1042-191115_1'};%,'1064-200316_2'};
color_code = {'b','r','k','c'};

%% df/f and spk
ave_dfOvF_runairpuff_allsession = [];% frames*total number of cells in all sessions
ave_dfOvF_stayairpuff_allsession = [];
aveSpd_run_all_session_100ms = [];% frames*total number of sessions that have runtrig
aveSpd_stay_all_session_100ms = [];
avespk_runairpuff_allsession_100ms = [];
avespk_stayairpuff_allsession_100ms = [];
nPCs_total = 0;
% ntrials_run_total = 0;
% ntrials_stay_total = 0;

% 'spd_airstimAve_allrun','spd_airstimAve_stay','spk_airpuff_allrun_cells','spk_airpuff_stay_cells'
% 'dfOvFbtm_airpuff_allrun_cells','dfOvFbtm_airpuff_stay_cells'
for ii = 1:length(sessions) 
    %load things
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\']; 
    dfOvF_strct = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    dfOvFbtm_airpuff_allrun_cells = dfOvF_strct.dfOvFbtm_airpuff_allrun_cells;%frame*cells
    dfOvFbtm_airpuff_stay_cells = dfOvF_strct.dfOvFbtm_airpuff_stay_cells;
    behav_dest = ['Z:\Analysis\Airpuff_analysis\behavioral_analysis\' days{ii}];
    behav_struct = load([behav_dest '\' days{ii} '_behavAnalysis.mat']);
    spd_airstim_allrun_100ms = behav_struct.spd_airstim_allrun_100ms; %trials * frames
    spd_airstim_stay_100ms = behav_struct.spd_airstim_stay_100ms;
  
    spk_deconv_output = load([image_analysis_base, sessions{ii}, '\deconv_wRunSmin\',sessions{ii} '_spk_deconvolve_staynrun_seperate.mat']);
    spk_airpuff_allrun_cells = spk_deconv_output.spk_airpuff_allrun_cells;
    spk_airpuff_stay_cells = spk_deconv_output.spk_airpuff_stay_cells;
    
    nPCs_total = nPCs_total + size(spk_airpuff_stay_cells,2);
    %ntrials_run_total = ntrials_run_total + size(,2);
    %ntrials_stay_total = ntrials_stay_total + size(,2);
    
    ave_dfOvF_runairpuff_allsession = cat(2,ave_dfOvF_runairpuff_allsession,dfOvFbtm_airpuff_allrun_cells);%when cat matrix, the matrix with NaN will be cat. And you have matrix with NaN because for some sessions there's no running trials fullfill the runtrig/runoff criteria
    ave_dfOvF_stayairpuff_allsession = cat(2,ave_dfOvF_stayairpuff_allsession,dfOvFbtm_airpuff_stay_cells);
    aveSpd_run_all_session_100ms = cat(1,aveSpd_run_all_session_100ms,spd_airstim_allrun_100ms);
    aveSpd_stay_all_session_100ms = cat(1,aveSpd_stay_all_session_100ms,spd_airstim_stay_100ms);
    avespk_runairpuff_allsession_100ms = cat(2,avespk_runairpuff_allsession_100ms,spk_airpuff_allrun_cells);
    avespk_stayairpuff_allsession_100ms = cat(2,avespk_stayairpuff_allsession_100ms,spk_airpuff_stay_cells);
end

ave_dfOvF_runairpuff_acrossessions = mean(ave_dfOvF_runairpuff_allsession,2); % average across cells
ste_dfOvF_runairpuff_acrossessions = std(ave_dfOvF_runairpuff_allsession,0,2)/sqrt(size(ave_dfOvF_runairpuff_allsession,2));

ave_dfOvF_stayairpuff_acrossessions = mean(ave_dfOvF_stayairpuff_allsession,2);
ste_dfOvF_stayairpuff_acrossessions = std(ave_dfOvF_stayairpuff_allsession,0,2)/sqrt(size(ave_dfOvF_stayairpuff_allsession,2));

aveFR_runairpuff_acrosssessions = mean(avespk_runairpuff_allsession_100ms,2)*30; % average across cells
steFR_runairpuff_acrosssessions = std(avespk_runairpuff_allsession_100ms,0,2)/sqrt(size(avespk_runairpuff_allsession_100ms,2))*30;

aveFR_stayairpuff_acrosssessions = mean(avespk_stayairpuff_allsession_100ms,2)*30;
steFR_stayairpuff_acrosssessions = std(avespk_stayairpuff_allsession_100ms,0,2)/sqrt(size(avespk_stayairpuff_allsession_100ms,2))*30;

aveSpd_run_acrossessions = mean(aveSpd_run_all_session_100ms,1);
steSpd_run_acrossessions = std(aveSpd_run_all_session_100ms,0,1)/sqrt(size(aveSpd_run_all_session_100ms,1));
aveSpd_stay_acrossessions = mean(aveSpd_stay_all_session_100ms,1);
steSpd_stay_acrossessions = std(aveSpd_stay_all_session_100ms,0,1)/sqrt(size(aveSpd_stay_all_session_100ms,1));


%%
%plot
airpuff_across_fig = figure;

x = (1: size(ave_dfOvF_runairpuff_acrossessions,1))/30- 0.5-(1/30);
subplot(3,2,1);
shadedErrorBar(x,ave_dfOvF_runairpuff_acrossessions,ste_dfOvF_runairpuff_acrossessions,{'color',[0.8431 0.0980 0.1098]});
%title('airpuff stim running'); 
ylabel('df/f');
xlim([-0.5,0.5]);ylim([0,0.6]);
vline(16/30- 0.5-(1/30),'k');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',8);
subplot(3,2,3);
shadedErrorBar(x,aveFR_runairpuff_acrosssessions,steFR_runairpuff_acrosssessions,{'color',[0.1922 0.6392 0.3294]});
ylabel('firing rate');
xlim([-0.5,0.5]); ylim([0 5]);
vline(16/30- 0.5-(1/30),'k');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',8);
subplot(3,2,5);
shadedErrorBar(x,aveSpd_run_acrossessions*2*3.1415926*7.5/128,steSpd_run_acrossessions*2*3.1415926*7.5/128,{'color','k'});
xlabel('time from airpuff (s)'); ylabel('speed(cm/s)'); 
xlim([-0.5,0.5]); ylim([-0.5 15]);
vline(16/30- 0.5-(1/30),'k');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',8);

subplot(3,2,2);
shadedErrorBar(x,ave_dfOvF_stayairpuff_acrossessions,ste_dfOvF_stayairpuff_acrossessions,{'color',[0.8431 0.0980 0.1098]});
%title('airpuff stim stay'); %ylabel('df/f');
xlim([-0.5,0.5]); ylim([0,0.6]);
vline(16/30- 0.5-(1/30),'k');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',8);
subplot(3,2,4);
shadedErrorBar(x,aveFR_stayairpuff_acrosssessions,steFR_stayairpuff_acrosssessions,{'color',[0.1922 0.6392 0.3294]});
%ylabel('firing rate');
xlim([-0.5,0.5]); ylim([0 5]);
vline(16/30- 0.5-(1/30),'k');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',8);
subplot(3,2,6);
shadedErrorBar(x,aveSpd_stay_acrossessions*2*3.1415926*7.5/128,steSpd_stay_acrossessions*2*3.1415926*7.5/128,{'color','k'});
xlabel('time from airpuff(s)'); %ylabel('speed(cm/s)');
xlim([-0.5,0.5]);ylim([-0.5 15]); 
vline(16/30- 0.5-(1/30),'k');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',8);

airpuff_across_fig.Units = 'centimeters';
airpuff_across_fig.Position = [3 3 11 10];
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',8);
fig_name = 'across_airpufftrig';
path = 'Z:\Analysis\Airpuff_analysis\across_sessions\deconv_wRunSmin\';
print(airpuff_across_fig,[path,fig_name],'-r600','-depsc');

% title('airpuff response across sessions');
% savefig(['Z:\Analysis\Airpuff_analysis\across_sessions\' 'across_sessions_puffResponse.fig']);


