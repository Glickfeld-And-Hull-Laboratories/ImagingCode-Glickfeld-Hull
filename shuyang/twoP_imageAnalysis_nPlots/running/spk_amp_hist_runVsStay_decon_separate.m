sessions = '190603_img1025'; 
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
image_analysis_dest = [image_analysis_base, sessions, '\'];
%load data, use data of only good cells
rawF_output = load([image_analysis_dest sessions, '_deconvolution_thresh-4_TCave_cl.mat']);
rawF = rawF_output.TCave_cl;

spk = load([image_analysis_base, sessions, '\deconvolution_sep\' sessions '_spk_deconvolve_staynrun_seperate.mat' ]);
spk_inx_stay = spk.spk_inx_stay;
spk_inx_run = spk.spk_inx_run;
peaks_run = [];
peaks_stay = [];
for c = 1: size(spk_inx_stay,2)
    peaks_stay = cat(1,peaks_stay,rawF(spk_inx_stay{c},c));
    peaks_run = cat(1,peaks_run,rawF(spk_inx_run{c},c));
end
hist_fig = figure; 
hist_fig.Units = 'centimeters';
hist_fig.Position = [1 1 7.5 6];
hold on;subplot(1,2,1);
h1 = histogram(peaks_run,'BinWidth',200);
h1.FaceColor = 'b'; %This way the color is the most different ,don't know how to change the colors respectively
h1.EdgeColor = 'b';
xlim([0 max(peaks_run)]);
xlabel('rawF of peaks');
ylabel('number of events');
%title(['rawF of spike peaks running(' num2str(threshold_run) 'std)']);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
subplot(1,2,2);
h2 = histogram(peaks_stay,'BinWidth',200);
h2.FaceColor = 'r';
h2.EdgeColor = 'r';
xlim([0 max(peaks_stay)]);
%xlabel('rawF of peaks');
title('stay');
%saveas(hist_fig,[image_analysis_base, sessions, '\deconvolution_sep\' sessions '_SpkeventAmp_hist_sep']);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
hold off;
fig_name = [sessions '_SpkeventAmp_hist_sep'];
path = [image_analysis_dest 'deconvolution_sep\'];
orient(hist_fig,'landscape');
print(hist_fig,[path,fig_name],'-r600','-dpdf');

