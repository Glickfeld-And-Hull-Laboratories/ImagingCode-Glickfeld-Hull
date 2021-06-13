% use this script after "CaAmp_runVsStayAirpuff"
% this script is comparing amplitude of well isolated Ca transients of all
% neurons in all sessions
% does the amplitude of Ca transient change under different conditions?
% condition 1 : stationary (not including 1s after running offset and 1s before running onset)
% condition 2: fast running (not including 1s after running onset and 1s? before running offset)
% condition 3: slow running

%% Section I: set paths and create analysis folders 
%define the directory and files
clear;
sessions = {'200116_img1041','200217_img1061','200225_img1049','200319_img1064','200319_img1064_2'};
days = {'1041-200116_1','1061-200217_1','1049-200225_1','1064-200319_1','1064-200319_2'};
image_analysis_base = 'Z:\Analysis\motorizedWheel_Analysis\running\imaging_analysis\'; 
color_code = {'c','r','y','g'};

%% rawF 
ave_rawF_stay_iso_Allsessions = zeros(length(sessions),30);
ave_rawF_runfast_iso_Allsessions = zeros(length(sessions),30);
ave_rawF_runslow_iso_Allsessions = zeros(length(sessions),30);
nPCs_Allsessions = zeros(1,length(sessions));
rawF_runfast_isoevents_neurons_all = [];
rawF_stay_isoevents_neurons_all = [];
rawF_runslow_isoevents_neurons_all = [];

for ii = 1: length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    isoCa_output = load([image_analysis_dest sessions{ii} '_isoCaEvent.mat']);
    ave_rawF_stay_iso_Allsessions(ii,:) = isoCa_output.ave_rawF_stay_iso_session;
    ave_rawF_runfast_iso_Allsessions(ii,:) = isoCa_output.ave_rawF_runfast_iso_session;
    ave_rawF_runslow_iso_Allsessions(ii,:) = isoCa_output.ave_rawF_runslow_iso_session;
    rawF_runfast_isoevents_neurons_all = cat(1,rawF_runfast_isoevents_neurons_all,isoCa_output.rawF_runfast_isoevents_neurons);
    rawF_stay_isoevents_neurons_all = cat(1,rawF_stay_isoevents_neurons_all,isoCa_output.rawF_stay_isoevents_neurons);
    rawF_runslow_isoevents_neurons_all = cat(1,rawF_runslow_isoevents_neurons_all, isoCa_output.rawF_runslow_isoevents_neurons);
    % total number of events from all neurons in each session
%     nevents_stay_Allsessions(ii) = isoCa_output.nevents_stay;
%     nevents_run_Allsessions(ii) = isoCa_output.nevents_run;
    
    % total number of neurons in each session
    nPCs_Allsessions(ii) = size(isoCa_output.rawF_stay_isoevents_neurons,1);  
   
end
% average across sessions
ave_rawF_stay_iso_across = mean(ave_rawF_stay_iso_Allsessions);
ave_rawF_runfast_iso_across = mean(ave_rawF_runfast_iso_Allsessions);
ave_rawF_runslow_iso_across = mean(ave_rawF_runslow_iso_Allsessions);

ste_rawF_stay_iso_across = std(rawF_stay_isoevents_neurons_all)/sqrt(size(rawF_stay_isoevents_neurons_all,1));
ste_rawF_runfast_iso_across = std(rawF_runfast_isoevents_neurons_all)/sqrt(size(rawF_runfast_isoevents_neurons_all,1));
ste_rawF_runslow_iso_across = std(rawF_runslow_isoevents_neurons_all)/sqrt(size(rawF_runslow_isoevents_neurons_all,1));
% ntotalevents_stay = sum(nevents_stay_Allsessions);
% ntotalevents_run = sum(nevents_run_Allsessions);
ntotalPCs = sum(nPCs_Allsessions);
x = (1:30)/30;
figure;
%haven't figured out how to use RGB colors in errorbar
errorbar(x,ave_rawF_stay_iso_across,ste_rawF_stay_iso_across,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
errorbar(x,ave_rawF_runfast_iso_across,ste_rawF_runfast_iso_across,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
errorbar(x,ave_rawF_runslow_iso_across,ste_rawF_runslow_iso_across,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
%legend(['stationary nevents=' num2str(ntotalevents_stay)],['running nevents=' num2str(ntotalevents_run)]);
legend('stationary', 'fast running','slow running');
text(0.75,6000, ['n total PCs=' num2str(ntotalPCs)]);
%ylim([0 0.6]);
xlabel('time(s)'); ylabel('rawF');
savefig(['Z:\Analysis\motorizedWheel_Analysis\running\across_sessions\' 'across_sessions_CaAmpRunvsStay_rawF.fig']);

%% df/F
ave_dfOvF_stay_iso_Allsessions = zeros(length(sessions),30);
ave_dfOvF_runfast_iso_Allsessions = zeros(length(sessions),30);
ave_dfOvF_runslow_iso_Allsessions = zeros(length(sessions),30);
nPCs_Allsessions = zeros(1,length(sessions));
dfOvF_runfast_isoevents_neurons_all = [];
dfOvF_stay_isoevents_neurons_all = [];
dfOvF_runslow_isoevents_neurons_all = [];

for ii = 1: length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    isoCa_output = load([image_analysis_dest,'deconv_wfastRunSmin\' sessions{ii} '_isoCaEvent.mat']);
    ave_dfOvF_stay_iso_Allsessions(ii,:) = isoCa_output.ave_dfOvF_btm_stay_iso_session;
    ave_dfOvF_runfast_iso_Allsessions(ii,:) = isoCa_output.ave_dfOvF_btm_runfast_iso_session;
    ave_dfOvF_runslow_iso_Allsessions(ii,:) = isoCa_output.ave_dfOvF_btm_runslow_iso_session;
    dfOvF_runfast_isoevents_neurons_all = cat(1,dfOvF_runfast_isoevents_neurons_all,isoCa_output.dfOvF_btm_runfast_isoevents_neurons);
    dfOvF_stay_isoevents_neurons_all = cat(1,dfOvF_stay_isoevents_neurons_all,isoCa_output.dfOvF_btm_stay_isoevents_neurons);
    dfOvF_runslow_isoevents_neurons_all = cat(1,dfOvF_runslow_isoevents_neurons_all, isoCa_output.dfOvF_btm_runslow_isoevents_neurons);
    % total number of events from all neurons in each session
%     nevents_stay_Allsessions(ii) = isoCa_output.nevents_stay;
%     nevents_run_Allsessions(ii) = isoCa_output.nevents_run;
    
    % total number of neurons in each session
    nPCs_Allsessions(ii) = size(isoCa_output.dfOvF_btm_stay_isoevents_neurons,1);  
   
end
% average across sessions
ave_dfOvF_stay_iso_across = mean(ave_dfOvF_stay_iso_Allsessions);
ave_dfOvF_runfast_iso_across = mean(ave_dfOvF_runfast_iso_Allsessions);
ave_dfOvF_runslow_iso_across = mean(ave_dfOvF_runslow_iso_Allsessions);

ste_dfOvF_stay_iso_across = std(dfOvF_stay_isoevents_neurons_all)/sqrt(size(dfOvF_stay_isoevents_neurons_all,1));
ste_dfOvF_runfast_iso_across = std(dfOvF_runfast_isoevents_neurons_all)/sqrt(size(dfOvF_runfast_isoevents_neurons_all,1));
ste_dfOvF_runslow_iso_across = std(dfOvF_runslow_isoevents_neurons_all)/sqrt(size(dfOvF_runslow_isoevents_neurons_all,1));
% ntotalevents_stay = sum(nevents_stay_Allsessions);
% ntotalevents_run = sum(nevents_run_Allsessions);
ntotalPCs = sum(nPCs_Allsessions);
x = (1:30)/30;
% figure;
% %haven't figured out how to use RGB colors in errorbar
% errorbar(x,ave_dfOvF_stay_iso_across,ste_dfOvF_stay_iso_across,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
% errorbar(x,ave_dfOvF_runfast_iso_across,ste_dfOvF_runfast_iso_across,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
% errorbar(x,ave_dfOvF_runslow_iso_across,ste_dfOvF_runslow_iso_across,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
% %legend(['stationary nevents=' num2str(ntotalevents_stay)],['running nevents=' num2str(ntotalevents_run)]);
% legend('stationary', 'fast running','slow running');
% text(0.75,6000, ['n total PCs=' num2str(ntotalPCs)]);
% %ylim([0 0.6]);
% xlabel('time(s)'); ylabel('df/f');
% savefig(['Z:\Analysis\motorizedWheel_Analysis\running\across_sessions\' 'across_sessions_CaAmpRunvsStay.fig']);

% plot start from the same y
diff1 = ave_dfOvF_stay_iso_across(1) - ave_dfOvF_runfast_iso_across(1);
diff2 = ave_dfOvF_runslow_iso_across(1) - ave_dfOvF_runfast_iso_across(1);
ave_dfOvF_stay_iso_across_plot = ave_dfOvF_stay_iso_across - diff1;
ave_dfOvF_runslow_iso_across_plot = ave_dfOvF_runslow_iso_across - diff2;
CaAmp_dfOvF = figure;
errorbar(x,ave_dfOvF_stay_iso_across_plot,ste_dfOvF_stay_iso_across,'.','LineStyle','-','linewidth', 1,'MarkerSize',4,'color','k','Capsize',3); hold on;
errorbar(x,ave_dfOvF_runfast_iso_across,ste_dfOvF_runfast_iso_across,'.',...
    'LineStyle','-','linewidth', 1,'MarkerSize',4,'color',[0, 0.4470, 0.7410],'Capsize',3); % matlab default color
errorbar(x,ave_dfOvF_runslow_iso_across_plot,ste_dfOvF_runslow_iso_across,'.',...
    'LineStyle','-','linewidth', 1,'MarkerSize',4,'color',[0.8500, 0.3250, 0.0980],'Capsize',3); %matlab default color

%legend(['stationary nevents=' num2str(ntotalevents_stay)],['running nevents=' num2str(ntotalevents_run)]);
legend('stationary', 'fast','slow'); %legend boxoff;
%text(0.75,0.55, ['n total PCs=' num2str(ntotalPCs)]);
ylim([0 0.85]); xlim([0 1.1]);
xlabel('time(s)'); ylabel('df/F');
CaAmp_dfOvF.Units = 'centimeters';
CaAmp_dfOvF.Position = [3 3 5 5];
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
fig_name = 'across_session_CaAmp_motorized_running';
path = 'Z:\Analysis\motorizedWheel_Analysis\running\across_sessions\deconv_wRunSmin\';
%orient(CaAmp_dfOvF,'landscape')
print(CaAmp_dfOvF,[path,fig_name],'-r600','-depsc');
%savefig(['Z:\Analysis\motorizedWheel_Analysis\running\across_sessions\' 'across_sessions_CaAmpRunvsStay_samey.fig']);

