%% Section I: set paths and create analysis folders for each session
%define the directory and files
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = {'191114_img1040','191115_img1039','191115_img1041','191115_img1042'};
days = {'1040-191114_1','1039-191115_1','1041-191115_1','1042-191115_1'};
image_analysis_base  = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder

%image_dest_base    = ['Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the movingDots analysis folder
% behavior analysis results 
color_code = {'c','r','y','g'};

%% plot PSTH for airpuff 
FR_airstimRun_cells_across = [];
FR_airstimStay_cells_across = [];

for ii = 1: length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\deconv_wRunSmin\'];
    spk_deconv_output = load([image_analysis_dest sessions{ii},'_spk_deconvolve_staynrun_seperate.mat']);
    FR_airstimAve_stay = spk_deconv_output.FR_airstimAve_stay;
    FR_airstimAve_run = spk_deconv_output.FR_airstimAve_allrun;
    FR_airstimRun_cells_across = cat(2,FR_airstimRun_cells_across,FR_airstimAve_run);
    FR_airstimStay_cells_across = cat(2,FR_airstimStay_cells_across,FR_airstimAve_stay);
end
FR_airstimRun_acrossSession = mean(FR_airstimRun_cells_across,2);
FR_airstimStay_acrossSession = mean(FR_airstimStay_cells_across,2);

% bin every 3 frames -> 0.1s, calculate average FR during every 0.1s
n = 3;% number of frames want to bin together
vecEnd = n*floor(length(FR_airstimAve_run)/n);% vector length can't be divided by 3, floor to the nearest number that can be divided by 3.
FRbin_airstim_across_stay1 = reshape(FR_airstimStay_acrossSession(1:vecEnd),[n,30/n]);
FRbin_airstim_across_stay = mean(FRbin_airstim_across_stay1);
FRbin_airstim_across_run1 = reshape(FR_airstimRun_acrossSession(1:vecEnd),[n,30/n]);
FRbin_airstim_across_run = mean(FRbin_airstim_across_run1);

PSTH = figure;
x = (0:0.1:1)-0.5;
subplot(1,2,1);
histogram('BinEdges',x, 'BinCounts', FRbin_airstim_across_run,'Facecolor',[0.9569 0.6471 0.5098]);
ylabel('firing rate (Hz)'); xlabel('time(s)'); 
%title('airpuff stim run');
ylim([0 3]);vline(0,'k');
%saveas(PSTH, [image_analysis_dest '\' days{ii} '_airpuff_PSTH']);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',8);

subplot(1,2,2)
histogram('BinEdges',x, 'BinCounts', FRbin_airstim_across_stay,'Facecolor',[0.9569 0.6471 0.5098]);
%title(['airpuff stim stay ' sessions{ii}]);
xlabel('time(s)');
ylim([0 3]);vline(0,'k');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',8);

PSTH.Units = 'centimeters';
PSTH.Position = [1 3 5.5 5];
fig_name = 'airpuff_PSTH_across';
path = 'Z:\Analysis\Airpuff_analysis\across_sessions\deconv_wRunSmin\';
print(PSTH,[path,fig_name],'-r600','-dpdf');
