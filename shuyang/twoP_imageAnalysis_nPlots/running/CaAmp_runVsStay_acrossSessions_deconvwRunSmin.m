% use this script after "CaAmp_runVsStayAirpuff"
% this script is comparing amplitude of well isolated Ca transients of all
% neurons in all sessions
% does the amplitude of Ca transient change under different conditions?
% condition 1 : stationary (not including 1s after running offset and 1s before running onset)
% condition 2: running (not including 1s after running onset and 1s? before running offset)
% condition 3: airpuff response during running (within 200ms of airpuff onset)
% condition 4: airpuff response during stationary (within 200ms of airpuff onset)
% events need to be isolated: no other events within 600ms before and/or after
% frm_stay_midpart: all stationary frames without beginning and end of stationary
% isolated_inx_neurons: index of isolated event peaks 1*#of cells, each cell element is a vector
% spk_iso_stay_neurons: index of isolated event peaks during stationary, 1*# of neurons, each cell element is a vector
% isoevent_stay_neurons: index of whole isolated event(1s) during stationary, 1*#of neurons, each cell element is a matrix(event*frames)
% dfOvF_btm_stay_isoevents: %df/f of whole isolated event(1s) during stationary, 1*#of neurons, each cell element is a matrix(event*frames)
% dfOvF_btm_stay_isoevents_neurons: average df/f of isolated event of each cell, averaged across events, neurons*frames
% ave_dfOvF_btm_stay_iso_session: average df/f of isolated event across cells, use this to plot

%% Section I: set paths and create analysis folders 
%define the directory and files
clear;
sessions = {'190429_img1021','190430_img1023','190507_img1024','190603_img1025'};
days = {'1021-190429_1','1023-190430_1','1024-190507_1','1025-190603_1'};
%sessions = {'190117_img1016'};
%days = {'1016-190117_1'};

image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; 
color_code = {'c','r','y','g'};

%% 
ave_dfOvF_btm_stay_iso_Allsessions = zeros(length(sessions),30);
ste_dfOvF_btm_stay_iso_Allsessions = zeros(length(sessions),30);
ave_dfOvF_btm_run_iso_Allsessions = zeros(length(sessions),30);
ste_dfOvF_btm_run_iso_Allsessions = zeros(length(sessions),30);
nevents_stay_Allsessions = zeros(1,length(sessions));
nevents_run_Allsessions = zeros(1,length(sessions));
nPCs_Allsessions = zeros(1,length(sessions));
nPCs_run = zeros(1,length(sessions));%cells included when comparing amplitude during running. only cells fire during early, middle, abd late running are included because we need to draw the same number of events from every condition
dfOvF_btm_run_isoevents_neurons_all = [];
dfOvF_btm_stay_isoevents_neurons_all = [];
dfOvF_spk_run_begin_neurons_all = [];
dfOvF_spk_run_middle_neurons_all = [];
dfOvF_spk_run_late_neurons_all = [];
isospk_run_begin_all = 0;
isospk_run_middle_all = 0;
isospk_run_late_all = 0;
for ii = 1: length(sessions)
    image_analysis_dest_deconv_staywrunSmin = [image_analysis_base, sessions{ii}, '\deconvolution_staywRunSmin\'];
    isoCa_output = load([image_analysis_dest_deconv_staywrunSmin '_isoCaEvent.mat']);
    ave_dfOvF_btm_stay_iso_Allsessions(ii,:) = isoCa_output.ave_dfOvF_btm_stay_iso_session;
    ave_dfOvF_btm_run_iso_Allsessions(ii,:) = isoCa_output.ave_dfOvF_btm_run_iso_session;
    dfOvF_btm_run_isoevents_neurons_all = cat(1,dfOvF_btm_run_isoevents_neurons_all,isoCa_output.dfOvF_btm_run_isoevents_neurons);
    dfOvF_btm_stay_isoevents_neurons_all = cat(1,dfOvF_btm_stay_isoevents_neurons_all,isoCa_output.dfOvF_btm_stay_isoevents_neurons);
    dfOvF_spk_run_begin_neurons_all = cat(1,dfOvF_spk_run_begin_neurons_all,isoCa_output.dfOvF_spk_run_begin_neurons);
    dfOvF_spk_run_middle_neurons_all = cat(1,dfOvF_spk_run_middle_neurons_all,isoCa_output.dfOvF_spk_run_middle_neurons);
    dfOvF_spk_run_late_neurons_all = cat(1,dfOvF_spk_run_late_neurons_all,isoCa_output.dfOvF_spk_run_late_neurons);
    nPCs_run (ii) = size(isoCa_output.dfOvF_spk_run_begin_neurons,1);
    isospk_run_begin_all = isospk_run_begin_all + isoCa_output.isospk_run_begin;
    isospk_run_middle_all = isospk_run_middle_all + isoCa_output.isospk_run_middle;
    isospk_run_late_all = isospk_run_late_all + isoCa_output.isospk_run_late;

    % total number of events from all neurons in each session
%     nevents_stay_Allsessions(ii) = isoCa_output.nevents_stay;
%     nevents_run_Allsessions(ii) = isoCa_output.nevents_run;
    
    % total number of neurons in each session
    nPCs_Allsessions(ii) = size(isoCa_output.dfOvF_btm_stay_isoevents_neurons,1);  
   
end
% average across sessions
ave_dfOvF_btm_stay_iso_across = mean(ave_dfOvF_btm_stay_iso_Allsessions);
ave_dfOvF_btm_run_iso_across = mean(ave_dfOvF_btm_run_iso_Allsessions);

% 
ste_dfOvF_btm_stay_iso_across = std(dfOvF_btm_stay_isoevents_neurons_all)/sqrt(size(dfOvF_btm_stay_isoevents_neurons_all,1));
ste_dfOvF_btm_run_iso_across = std(dfOvF_btm_run_isoevents_neurons_all)/sqrt(size(dfOvF_btm_run_isoevents_neurons_all,1));

% ntotalevents_stay = sum(nevents_stay_Allsessions);
% ntotalevents_run = sum(nevents_run_Allsessions);
ntotalPCs = sum(nPCs_Allsessions);
x = (1:30)/30;
f1 = figure;
errorbar(x,ave_dfOvF_btm_stay_iso_across,ste_dfOvF_btm_stay_iso_across,'.','LineStyle','-','linewidth', 1,'MarkerSize',4,'Capsize',3); hold on;
errorbar(x,ave_dfOvF_btm_run_iso_across,ste_dfOvF_btm_run_iso_across,'.','LineStyle','-','linewidth', 1,'MarkerSize',4,'Capsize',3); hold on;
%legend(['stationary nevents=' num2str(ntotalevents_stay)],['running nevents=' num2str(ntotalevents_run)]);
legend('statinary','running'); %legend('boxoff');
%text(0.75,0.45, ['n total PCs=' num2str(ntotalPCs)]);
ylim([0 0.6]);
xlabel('time(s)'); ylabel('df/F');
f1.Units = 'centimeters';
f1.Position = [3 3 5 5];
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
path = 'Z:\Analysis\2P_MovingDots_Analysis\across_sessions\deconvolve_staywRunSmin\';
fig_name = 'across_session_CaAmp';
print(f1,[path,fig_name],'-r600','-depsc');
savefig([path fig_name '.fig']);

%plot: make the y start at the same place
diff = ave_dfOvF_btm_stay_iso_across(1)-ave_dfOvF_btm_run_iso_across(1);
ave_dfOvF_btm_stay_iso_across_plot = ave_dfOvF_btm_stay_iso_across - diff;
x = (1:30)/30;
CaAmp_dfOvF = figure;
%haven't figured out how to use RGB colors in errorbar
errorbar(x,ave_dfOvF_btm_stay_iso_across_plot,ste_dfOvF_btm_stay_iso_across,'.','LineStyle','-','linewidth', 1,'MarkerSize',4,'Capsize',3); hold on;
errorbar(x,ave_dfOvF_btm_run_iso_across,ste_dfOvF_btm_run_iso_across,'.','LineStyle','-','linewidth', 1,'MarkerSize',4,'Capsize',3); hold on;
legend('stationary','running');%legend('boxoff');
axis square;
%text(0.75,0.45, ['n total PCs=' num2str(ntotalPCs)]);
ylim([0 0.6]);xlim([0 1.1]);
xlabel('time(s)'); ylabel('df/F');
CaAmp_dfOvF.Units = 'centimeters';
CaAmp_dfOvF.Position = [3 3 5 5];
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
fig_name = 'across_session_CaAmp_samey';
orient(CaAmp_dfOvF,'landscape')
print(CaAmp_dfOvF,[path,fig_name],'-r600','-depsc');
savefig([path fig_name '.fig']);

%ttest
[h,p] = ttest(ave_dfOvF_btm_stay_iso_across_plot,ave_dfOvF_btm_run_iso_across);

% when does isolated event happen during running?
total_event = isospk_run_begin_all+isospk_run_middle_all+isospk_run_late_all;
p_isospk_run_begin_all = isospk_run_begin_all*100/total_event;
p_isospk_run_middle_all = isospk_run_middle_all*100/total_event;
p_isospk_run_late_all = isospk_run_late_all*100/total_event;

figure; 
bar([p_isospk_run_begin_all,p_isospk_run_middle_all,p_isospk_run_late_all]);
set(gca,'XTickLabel',{'0-0.5s','0.5-1s','later than 1s'});
ylabel('percentage (%)');
title('time of isolated peaks during running windows');
savefig(['Z:\Analysis\2P_MovingDots_Analysis\across_sessions\' 'across_sessions_isoCaEvent_run_time.fig']);

% how does time of event during running influecne amplitude
figure;
fast_errbar(x,dfOvF_spk_run_begin_neurons_all,1,'color',[0.1373 0.5451 0.2706]);hold on;
fast_errbar(x,dfOvF_spk_run_middle_neurons_all,1,'color',[0.2549 0.6706 0.3647]); hold on;
fast_errbar(x,dfOvF_spk_run_late_neurons_all,1,'color',[0.7294 0.8941 0.7020]); hold on;
legend('0-0.5s','0.5-1s','later than 1s');
text(0.75,0.3,['n total PCs = ' num2str(sum(nPCs_run))]);
title('rawF of isolated events during different time period of running trials');
ylabel('rawF');
xlabel('time(s)');
savefig([path 'across_sessions_CaAmpRun3conditions.fig']);



