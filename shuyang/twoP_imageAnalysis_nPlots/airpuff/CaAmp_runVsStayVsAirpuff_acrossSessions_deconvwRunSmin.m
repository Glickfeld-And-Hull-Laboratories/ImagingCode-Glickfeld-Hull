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
sessions = {'191114_img1040','191115_img1039','191115_img1041','191115_img1042'};
days = {'1040-191114_1','1039-191115_1','1041-191115_1','1042-191115_1'};
image_analysis_base = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\'; color_code = {'c','r','y','g'};

%% 
ave_dfOvF_btm_stay_iso_Allsessions = zeros(length(sessions),30);
ave_dfOvF_btm_run_iso_Allsessions = zeros(length(sessions),30);
ave_dfOvF_airpuff_stay_iso_Allsessions = zeros(length(sessions),30);
ave_dfOvF_airpuff_run_iso_Allsessions = zeros(length(sessions),30);
nPCsrun_Allsessions = zeros(1,length(sessions));
nPCsstay_Allsessions = zeros(1,length(sessions));
dfOvF_btm_run_isoevents_neurons_all = [];
dfOvF_btm_stay_isoevents_neurons_all = [];
dfOvF_btm_airpuff_stay_isoevents_neurons_all = [];
dfOvF_btm_airpuff_run_isoevents_neurons_all = [];

for ii = 1: length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    isoCa_output = load([image_analysis_dest '\deconv_wRunSmin\' sessions{ii} '_isoCaEvent.mat']);
    ave_dfOvF_btm_stay_iso_Allsessions(ii,:) = isoCa_output.ave_dfOvF_btm_stay_iso_session;
    ave_dfOvF_btm_run_iso_Allsessions(ii,:) = isoCa_output.ave_dfOvF_btm_run_iso_session;
    ave_dfOvF_airpuff_stay_iso_Allsessions(ii,:) = isoCa_output.ave_dfOvF_btm_airpuff_stay_iso_session;
    ave_dfOvF_airpuff_run_iso_Allsessions(ii,:) = isoCa_output.ave_dfOvF_btm_airpuff_run_iso_session;
    dfOvF_btm_run_isoevents_neurons_all = cat(1,dfOvF_btm_run_isoevents_neurons_all,isoCa_output.dfOvF_btm_run_isoevents_neurons);
    dfOvF_btm_stay_isoevents_neurons_all = cat(1,dfOvF_btm_stay_isoevents_neurons_all,isoCa_output.dfOvF_btm_stay_isoevents_neurons);
    dfOvF_btm_airpuff_stay_isoevents_neurons_all = cat(1,dfOvF_btm_airpuff_stay_isoevents_neurons_all, isoCa_output.dfOvF_btm_airpuff_stay_isoevents_neurons);
    dfOvF_btm_airpuff_run_isoevents_neurons_all = cat(1,dfOvF_btm_airpuff_run_isoevents_neurons_all, isoCa_output.dfOvF_btm_airpuff_run_isoevents_neurons);
    % total number of events from all neurons in each session
%     nevents_stay_Allsessions(ii) = isoCa_output.nevents_stay;
%     nevents_run_Allsessions(ii) = isoCa_output.nevents_run;
    
    % total number of neurons in each session
    nPCsrun_Allsessions(ii) = size(isoCa_output.dfOvF_btm_run_isoevents_neurons,1);
    nPCsstay_Allsessions(ii) = size(isoCa_output.dfOvF_btm_stay_isoevents_neurons,1);
    % you will have less cells here than the airpuff plot across sessions,
    % because we are drawing some # of events during running and
    % stationary, thus if a cell doesn't have well-isolated events during
    % running, then that cell is being sorted out
    
end
% average across sessions
ave_dfOvF_btm_stay_iso_across = mean(ave_dfOvF_btm_stay_iso_Allsessions);
ave_dfOvF_btm_run_iso_across = mean(ave_dfOvF_btm_run_iso_Allsessions);
ave_dfOvF_btm_airpuff_stay_iso_across = mean(ave_dfOvF_airpuff_stay_iso_Allsessions);
ave_dfOvF_btm_airpuff_run_iso_across = mean(ave_dfOvF_airpuff_run_iso_Allsessions);

ste_dfOvF_btm_stay_iso_across = std(dfOvF_btm_stay_isoevents_neurons_all)/sqrt(size(dfOvF_btm_stay_isoevents_neurons_all,1));
ste_dfOvF_btm_run_iso_across = std(dfOvF_btm_run_isoevents_neurons_all)/sqrt(size(dfOvF_btm_run_isoevents_neurons_all,1));
ste_dfOvF_btm_airpuff_stay_iso_across = std(dfOvF_btm_airpuff_stay_isoevents_neurons_all)/sqrt(size(dfOvF_btm_airpuff_stay_isoevents_neurons_all,1));
ste_dfOvF_btm_airpuff_run_iso_across = std(dfOvF_btm_airpuff_run_isoevents_neurons_all)/sqrt(size(dfOvF_btm_airpuff_run_isoevents_neurons_all,1));

% ntotalevents_stay = sum(nevents_stay_Allsessions);
% ntotalevents_run = sum(nevents_run_Allsessions);
ntotalPCs_run = sum(nPCsrun_Allsessions);
ntotalPCs_stay = sum(nPCsstay_Allsessions);

x = (1:30)/30;
% figure;
% errorbar(x,ave_dfOvF_btm_stay_iso_across,ste_dfOvF_btm_stay_iso_across,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
% errorbar(x,ave_dfOvF_btm_run_iso_across,ste_dfOvF_btm_run_iso_across,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); 
% errorbar(x,ave_dfOvF_btm_airpuff_stay_iso_across,ste_dfOvF_btm_airpuff_stay_iso_across,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); 
% errorbar(x,ave_dfOvF_btm_airpuff_run_iso_across,ste_dfOvF_btm_airpuff_run_iso_across,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); 

% %legend(['stationary nevents=' num2str(ntotalevents_stay)],['running nevents=' num2str(ntotalevents_run)]);
% legend('stationary', 'running','airpuff response_stationary');
% text(0.75,0.45, ['n total PCs=' num2str(ntotalPCs)]);
% %ylim([0 0.6]);
% xlabel('time(s)'); ylabel('df/f');
% savefig(['Z:\Analysis\Airpuff_analysis\across_sessions\' 'across_sessions_CaAmpRunvsStayVsAirpuff.fig']);
%%
% start with same y
diff1 = ave_dfOvF_btm_stay_iso_across(1) - ave_dfOvF_btm_airpuff_stay_iso_across(1);
diff2 = ave_dfOvF_btm_run_iso_across(1) - ave_dfOvF_btm_airpuff_stay_iso_across(1);
diff3 = ave_dfOvF_btm_airpuff_run_iso_across(1) - ave_dfOvF_btm_airpuff_stay_iso_across(1);
ave_dfOvF_btm_stay_iso_across_plot = ave_dfOvF_btm_stay_iso_across - diff1;
ave_dfOvF_btm_run_iso_across_plot = ave_dfOvF_btm_run_iso_across - diff2;
ave_dfOvF_btm_airpuff_run_iso_across_plot = ave_dfOvF_btm_airpuff_run_iso_across - diff3;
CaAmp_dfOvF = figure;hold on;
errorbar(x,ave_dfOvF_btm_stay_iso_across_plot,ste_dfOvF_btm_stay_iso_across,...
    '.','LineStyle','-','linewidth', 1,'MarkerSize',4,'color',[0, 0.4470, 0.7410],'Capsize',3); % matlab default blue color
errorbar(x,ave_dfOvF_btm_run_iso_across_plot,ste_dfOvF_btm_run_iso_across,...
    '.','LineStyle','-','linewidth', 1,'MarkerSize',4,'color',[0.8500, 0.3250, 0.0980],'Capsize',3); % matlab default red color
errorbar(x,ave_dfOvF_btm_airpuff_stay_iso_across,ste_dfOvF_btm_airpuff_stay_iso_across,...
    '.','LineStyle','-','linewidth', 1,'MarkerSize',4,'color','k','Capsize',3);
errorbar(x,ave_dfOvF_btm_airpuff_run_iso_across_plot,ste_dfOvF_btm_airpuff_run_iso_across,...
    '.','LineStyle','-','linewidth', 1,'MarkerSize',4,'color','g','Capsize',3);
%legend(['stationary nevents=' num2str(ntotalevents_stay)],['running nevents=' num2str(ntotalevents_run)]);
legend('stationary', 'running','airpuff response stay','airpuff response run'); %legend boxoff;
%text(0.75,0.45, ['n total PCs=' num2str(ntotalPCs)]);
ylim([0 0.8]);xlim([0 1.1])
xlabel('time(s)'); ylabel('df/f');
CaAmp_dfOvF.Units = 'centimeters';
CaAmp_dfOvF.Position = [3 3 5 5];
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
fig_name = 'across_session_CaAmpRunvsStayVsrunAirpuffVsstayAirpuff_samey';
path = 'Z:\Analysis\Airpuff_analysis\across_sessions\deconv_wRunSmin\';
print(CaAmp_dfOvF,[path,fig_name],'-r600','-depsc');
savefig(['Z:\Analysis\Airpuff_analysis\across_sessions\deconv_wRunSmin\' 'across_sessions_CaAmpRunvsStayVsrunAirpuffVsstayAirpuff_samey.fig']);




