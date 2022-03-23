% spike identification using deconvolution. (function deconvolve Ca got from OASIS on git hub)
% to use the deconvolution function, need to add OASIS function to search
% path: >>oasis_setup
% function deconvolve Ca: denoise raw fluorescence trace and identifies spike events. 
% right now using FOOPSI method and the model is ar1 (also tried ar2: see commented bottom part),ar1 works better than ar2.
% then SECTION II draws the GUI showing raw traces with identified events and deconvolved trace.
% 
%% assign document paths and experimental sessions
%sessions = {'190429_img1021','190430_img1023','190507_img1024','190603_img1025'};
%days = {'1021-190429_1','1023-190430_1','1024-190507_1','1025-190603_1'};

clear;
sessions = '190603_img1025'; 
days = '1025-190603_1';

image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
image_analysis_dest = [image_analysis_base, sessions, '\'];
% behavior analysis results
behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days '\'];
%load data, use data of only good cells
rawF_output = load([image_analysis_dest sessions, '_deconvolution_thresh-4_TCave_cl.mat']);
TCave = rawF_output.TCave_cl;
rawF = rawF_output.TCave_cl;
% dfOvF_strct = load([image_analysis_dest sessions '_dfOvF.mat']);
% dfOvF_btm = dfOvF_strct.dfOvF_btm_cl;

behav_output = load([behav_dest days '_behavAnalysis.mat']);
frm_stay_cell = behav_output.frames_stay_cell;
frm_stay = cell2mat(frm_stay_cell);
frm_run_cell = behav_output.frames_run_cell;
frm_run = cell2mat(frm_run_cell);

%% SECTION I deconvolution
% input: fluorescence trace: T*1, threshold for spike size (how big the spikes are supposed to be: normal set as -3(3 sigma, 3 times bigger than noise levels)).
% kernal: denoised trace, T*1 vector
% spikes: all of the values in a spike(thus, you will get more than 1 values if a spike is longer than a frame. e.g.: if a spike is 120ms long, you will get 4 values back to back which all belong to a single spike)
    
frames = 1:1:size(TCave,1);
threshold_run = -4; %changing the threhold basically changes the identification of spikes when the peak amplitude is small (those small peaks), doesn't change anything with the bigger jittered ones. -3 or -3.5 gives a FR close to 1 during stationary

% deconvolution
kernel_run = zeros(length(frm_run),size(TCave,2));
spk_run = zeros(length(frm_run),size(TCave,2));
smin_run = zeros(1,size(TCave,2)); % smin is the threshold (in fluorescence values) of smallest spikes being detected
% spk_peak_run = {};
% spk_inx_run = {};
% spk_logic_run = zeros(length(frm_run),size(TCave,2));
% num_spks_cell_run = zeros(1,size(TCave,2));
% FRrun_cell = zeros(1,size(TCave,2));
for c = 1: size(TCave,2)
    [kernel_run(:,c), spk_run(:,c), options] = deconvolveCa(TCave(frm_run,c), 'optimize_pars', true, ...
        'optimize_b', true, 'method','foopsi', 'smin', threshold_run);
    smin_run(c) = options.smin;
%     % get only the peaks of each spike
%     [spk_peak_run{c},spk_inx_run{c}] = findpeaks(spk_run(:,c));
%     spk_inx_run{c} = frm_run(spk_inx_run{c});
%     % spike logic
%     spk_logic_run(:,c) = (ismember(frm_run,spk_inx_run{c}))';
%     num_spks_cell_run(c) = sum(spk_logic_run(:,c)==1);
%     FRrun_cell(c)= num_spks_cell_run(c)/length(frm_run)*30; % firing rate = # of spikes/duration(s)
end
%aveFR_run = mean(FRrun_cell);

% use the smin from running to do deconvolution for the whole experiment
kernel = zeros(size(TCave,1),size(TCave,2));
spk = zeros(size(TCave,1),size(TCave,2));
spk_peak = {};
spk_inx = {};
spk_logic = zeros(size(TCave,1),size(TCave,2));
num_spks_cell = zeros(1,size(TCave,2));
FRstay_cell = zeros(1,size(TCave,2));
for c = 1: size(TCave,2)
    [kernel(:,c), spk(:,c), options] = deconvolveCa_stay(TCave(:,c), 'optimize_pars', true, ...
        'optimize_b', true, 'method','foopsi', 'smin', smin_run(c)); % deconvolveCa_stay is the function I modified that applies the smin of running to the whole dataset (only changed for FOOPSI model OR1) I think I just assigned smin_run(c) to gmax
    % get only the peaks of each spike
    [spk_peak{c},spk_inx{c}] = findpeaks(spk(:,c)); % the index you get here is the index in frm_stay
    % spike logic
    spk_logic(:,c) = (ismember(frames,spk_inx{c}))';
    num_spks_cell(c) = sum(spk_logic(frm_stay,c)==1);
    FRstay_cell(c)= num_spks_cell(c)/length(frm_stay)*30; % firing rate = # of spikes/duration(s)
end
aveFR_stay = mean(FRstay_cell);

peaks_run = [];
peaks_stay = [];
for c = 1: size(spk_inx,2)
    temp_inx = spk_inx{c};
    spk_inx_stay = intersect(frm_stay,temp_inx);
    spk_inx_run = intersect(frm_run,temp_inx);
    peaks_stay = cat(1,peaks_stay,rawF(spk_inx_stay,c));
    peaks_run = cat(1,peaks_run,rawF(spk_inx_run,c));
end

hist_fig = figure; 
hist_fig.Units = 'centimeters';
hist_fig.Position = [1 1 7.5 6];
hold on;subplot(1,2,1);
h1 = histogram(peaks_run,'BinWidth',200);
h1.FaceColor = 'b'; %This way the color is the most different ,don't know how to change the colors respectively
h1.EdgeColor = 'b';
xlim([-0.1 max(peaks_run)]);
xlabel('rawF of peaks');
ylabel('number of events');
%title(['rawF of spike peaks running(' num2str(threshold_run) 'std)']);
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'FontSize',7);
subplot(1,2,2);
h2 = histogram(peaks_stay,'BinWidth',200);
h2.FaceColor = 'r';
h2.EdgeColor = 'r';
xlim([-0.1 max(peaks_stay)]);
%xlabel('rawF of peaks');
title('stay using running smin');
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'FontSize',7);
hold off;
fig_name = [sessions '_SpkeventAmp_hist_sep'];
image_analysis_dest_deconv_staywRunSmin = [image_analysis_base, sessions, '\deconvolution_staywRunSmin\'];
path = image_analysis_dest_deconv_staywRunSmin;
orient(hist_fig,'landscape');
print(hist_fig,[path,fig_name],'-r600','-dpdf');
% if ~exist(image_analysis_dest_deconv_staywRunSmin)
%     mkdir(image_analysis_dest_deconv_staywRunSmin);
% end
saveas(hist_fig,[image_analysis_dest_deconv_staywRunSmin fig_name]);

save([image_analysis_dest_deconv_staywRunSmin sessions '_spk_deconvolve_staynrun_seperate.mat' ],...
    'threshold_run','smin_run','FRstay_cell', 'aveFR_stay','options',...
    'spk_logic','spk','kernel','spk_peak','spk_inx');

%%
%GUI
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
image_analysis_dest_deconv_staywRunSmin = [image_analysis_base, sessions, '\deconvolution_staywRunSmin\'];
deconvolve_output = load([image_analysis_dest_deconv_staywRunSmin sessions, '_spk_deconvolve_staynrun_seperate.mat']);
threshold_run = deconvolve_output.threshold_run;
kernel = deconvolve_output.kernel;
spk = deconvolve_output.spk;

% [fig_deconvolve_run] = GUI_rawTrace_denoiseNdeconv(TCave(frm_run,:),kernel_run,spk_run,sessions,threshold_run);
% savefig([image_analysis_dest_deconv_staywRunSmin sessions '_GUI_TCave_deconvolution_running_threshold' num2str(threshold_run) '.fig']);

[fig_deconvolve_stay] = GUI_rawTrace_denoiseNdeconv(TCave,kernel,spk,sessions,threshold_run);
savefig([image_analysis_dest_deconv_staywRunSmin sessions '_GUI_TCave_deconvolution_whole_exper' '.fig']);
