% spike identification using deconvolution. (function deconvolve Ca got from OASIS on git hub)
% to use the deconvolution function, need to add OASIS function to search
% path: >>oasis_setup
% function deconvolve Ca: denoise raw fluorescence trace and identifies spike events. 
% right now using FOOPSI method and the model is ar1 (also tried ar2: see commented bottom part),ar1 works better than ar2.
% then SECTION II draws the GUI showing raw traces with identified events and deconvolved trace.
% 
% in this script, deconvolution during running and stationary are done separately
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
    
threshold_stay = -4; %changing the threhold basically changes the identification of spikes when the peak amplitude is small (those small peaks), doesn't change anything with the bigger jittered ones. -3 or -3.5 gives a FR close to 1 during stationary

% deconvolution
kernel_stay = zeros(length(frm_stay),size(TCave,2));
spk_stay = zeros(length(frm_stay),size(TCave,2));
spk_peak_stay = {};
spk_inx_stay = {};
spk_logic_stay = zeros(length(frm_stay),size(TCave,2));
num_spks_cell_stay = zeros(1,size(TCave,2));
FRstay_cell = zeros(1,size(TCave,2));
for c = 1: size(TCave,2)
    [kernel_stay(:,c), spk_stay(:,c), options] = deconvolveCa(TCave(frm_stay,c), 'optimize_pars', true, ...
        'optimize_b', true, 'method','foopsi', 'smin', threshold_stay);
    % get only the peaks of each spike
    [spk_peak_stay{c},spk_inx_stay{c}] = findpeaks(spk_stay(:,c)); % the index you get here is the index in frm_stay
    spk_inx_stay{c} = frm_stay(spk_inx_stay{c}); %now the index is back to frames from 1:end of whole experiment
    % spike logic
    spk_logic_stay(:,c) = (ismember(frm_stay,spk_inx_stay{c}))';
    num_spks_cell_stay(c) = sum(spk_logic_stay(:,c)==1);
    FRstay_cell(c)= num_spks_cell_stay(c)/length(frm_stay)*30; % firing rate = # of spikes/duration(s)
end
aveFR_stay = mean(FRstay_cell);

threshold_run = -4; %changing the threhold basically changes the identification of spikes when the peak amplitude is small (those small peaks), doesn't change anything with the bigger jittered ones. -3 or -3.5 gives a FR close to 1 during stationary

% deconvolution
kernel_run = zeros(length(frm_run),size(TCave,2));
spk_run = zeros(length(frm_run),size(TCave,2));
spk_peak_run = {};
spk_inx_run = {};
spk_logic_run = zeros(length(frm_run),size(TCave,2));
num_spks_cell_run = zeros(1,size(TCave,2));
FRrun_cell = zeros(1,size(TCave,2));
for c = 1: size(TCave,2)
    [kernel_run(:,c), spk_run(:,c), options] = deconvolveCa(TCave(frm_run,c), 'optimize_pars', true, ...
        'optimize_b', true, 'method','foopsi', 'smin', threshold_run);
    % get only the peaks of each spike
    [spk_peak_run{c},spk_inx_run{c}] = findpeaks(spk_run(:,c));
    spk_inx_run{c} = frm_run(spk_inx_run{c});
    % spike logic
    spk_logic_run(:,c) = (ismember(frm_run,spk_inx_run{c}))';
    num_spks_cell_run(c) = sum(spk_logic_run(:,c)==1);
    FRrun_cell(c)= num_spks_cell_run(c)/length(frm_run)*30; % firing rate = # of spikes/duration(s)
end
aveFR_run = mean(FRrun_cell);


peaks_run = [];
peaks_stay = [];
for c = 1: size(spk_inx_stay,2)
    peaks_stay = cat(1,peaks_stay,rawF(spk_inx_stay{c},c));
    peaks_run = cat(1,peaks_run,rawF(spk_inx_run{c},c));
end

hist = figure; hold on;subplot(1,2,1);
h1 = histogram(peaks_run,'BinWidth',200);
h1.FaceColor = 'b'; %This way the color is the most different ,don't know how to change the colors respectively
h1.EdgeColor = 'b';
xlim([-0.1 max(peaks_run)]);
xlabel('rawF of peaks');
ylabel('number of events');
title(['rawF of spike peaks running(' num2str(threshold_run) 'std)']);
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'FontSize',18);
subplot(1,2,2);
h2 = histogram(peaks_stay,'BinWidth',200);
h2.FaceColor = 'r';
h2.EdgeColor = 'r';
xlim([-0.1 max(peaks_stay)]);
xlabel('rawF of peaks');
title(['stay(' num2str(threshold_stay) 'std)']);
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'FontSize',18);
hold off;

image_analysis_dest_deconv_sep = [image_analysis_base, sessions, '\deconvolution_sep\'];
if ~exist(image_analysis_dest_deconv_sep)
    mkdir(image_analysis_dest_deconv_sep);
end
saveas(hist,[image_analysis_dest_deconv_sep sessions '_SpkeventAmp_hist_sep']);

save([image_analysis_dest_deconv_sep sessions '_spk_deconvolve_staynrun_seperate.mat' ],...
    'threshold_stay','threshold_run','FRstay_cell', 'aveFR_stay','aveFR_run','num_spks_cell_stay',...
    'num_spks_cell_run','options','spk_logic_stay','spk_logic_run','spk_stay',...
    'spk_run','kernel_stay','kernel_run','spk_peak_stay','spk_peak_run','spk_inx_stay','spk_inx_run');

%%
%GUI
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
image_analysis_dest_deconv_sep = [image_analysis_base, sessions, '\deconvolution_sep\'];
deconvolve_output = load([image_analysis_dest_deconv_sep sessions, '_spk_deconvolve_staynrun_seperate.mat']);
kernel_run = deconvolve_output.kernel_run;
spk_run = deconvolve_output.spk_run;
threshold_run = deconvolve_output.threshold_run;
kernel_stay = deconvolve_output.kernel_stay;
spk_stay = deconvolve_output.spk_stay;
threshold_stay = deconvolve_output.threshold_stay;

[fig_deconvolve_run] = GUI_rawTrace_denoiseNdeconv(TCave(frm_run,:),kernel_run,spk_run,sessions,threshold_run);
savefig([image_analysis_dest_deconv_sep sessions '_GUI_TCave_deconvolution_running_threshold' num2str(threshold_run) '.fig']);

[fig_deconvolve_stay] = GUI_rawTrace_denoiseNdeconv(TCave(frm_stay,:),kernel_stay,spk_stay,sessions,threshold_stay);
savefig([image_analysis_dest_deconv_sep sessions '_GUI_TCave_deconvolution_stay_threshold' num2str(threshold_stay) '.fig']);

%%
% plot top 10% of spikes, 10%-20%, 20%-30%... and see what amplitude of
% peaks gives you the spikes that look like real spike

image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
image_analysis_dest_deconv_sep = [image_analysis_base, sessions, '\deconvolution_sep\'];

deconvolve_output = load([image_analysis_dest_deconv_sep sessions, '_spk_deconvolve_staynrun_seperate.mat']);
spk_inx_stay = deconvolve_output.spk_inx_stay;
spk_inx_run = deconvolve_output.spk_inx_run;

rawF_perc_peak_stay_cells = zeros(size(spk_inx_stay,2),10,30);%cell*percentage categories*frame
for c = 1: size(spk_inx_stay,2)
    temp_cell_peaks = rawF(spk_inx_stay{c},c);
    spk_sort = sort(temp_cell_peaks,1,'descend'); %sort, biggest first
    sum_sort = floor(length(spk_sort)/10)*10; %round this list to the nearest number that can be divided by 10
    len_p = sum_sort/10; %length of each percentage catogery(how many spikes for top 10%, 10-20%.....)
    spk_inx_stay_cell = spk_inx_stay{c};
    for p = 1:10 % take every 10%
        inx = (1+(p-1)*len_p: len_p*p);
        spk_p = spk_sort(inx);
        spk_p_logic = ismember(temp_cell_peaks,spk_p);
        spk_p_inx_stay = spk_inx_stay_cell(spk_p_logic==1);
        whole_peak_stay = zeros(length(spk_p_inx_stay),30); %whole peak indices of bottom 10% spikes in terms of amplitude
        rawF_peak_stay_p = zeros(length(spk_p_inx_stay),30); %peaks*frames
        for s = 1:length(spk_p_inx_stay)
            if spk_p_inx_stay(s) > 11 % if the spike happens before the 11 frame after experiment starts, there's no way to pull out 11 frames before then
                whole_peak_stay(s,:) = spk_p_inx_stay(s)-11:spk_p_inx_stay(s)+18;
                rawF_peak_stay_p(s,:) = rawF(whole_peak_stay(s,:),c);
            end
        end
        rawF_perc_peak_stay_cells(c,p,:) = mean(rawF_peak_stay_p); %average across peaks
    end
end
rawF_perc_peak_stay = squeeze(mean(rawF_perc_peak_stay_cells)); %average across cells
ste_rawF_perc_peak_stay = squeeze(std(rawF_perc_peak_stay_cells)/sqrt(size(rawF_perc_peak_stay_cells,1)));
%calculate the peak size (peak value-baseline)


rawF_perc_peak_run_cells = zeros(size(spk_inx_run,2),10,30);%cell*percentage categories*frame
for c = 1: size(spk_inx_run,2)
    temp_cell_peaks = rawF(spk_inx_run{c},c);
    spk_sort = sort(temp_cell_peaks,1,'descend'); %sort, biggest first
    sum_sort = floor(length(spk_sort)/10)*10; %round this list to the nearest number that can be divided by 10
    len_p = sum_sort/10; %length of each percentage catogery(how many spikes for top 10%, 10-20%.....), if a cell has too few spikes (less than 10 spikes during running),len_p is gonna be zero
    spk_inx_run_cell = spk_inx_run{c};
    for p = 1:10 % take every 10%
        inx = (1+(p-1)*len_p: len_p*p);
        spk_p = spk_sort(inx);
        spk_p_logic = ismember(temp_cell_peaks,spk_p);
        spk_p_inx_run = spk_inx_run_cell(spk_p_logic==1);
        whole_peak_run = zeros(length(spk_p_inx_run),30); %whole peak indices of bottom 10% spikes in terms of amplitude
        rawF_peak_run_p = zeros(length(spk_p_inx_run),30); %peaks*frames
        for s = 1:length(spk_p_inx_run)
            if spk_p_inx_run(s) > 11 % if the spike happens before the 11 frame after experiment starts, there's no way to pull out 11 frames before then
                whole_peak_run(s,:) = spk_p_inx_run(s)-11:spk_p_inx_run(s)+18;
                rawF_peak_run_p(s,:) = rawF(whole_peak_run(s,:),c);
            end
        end
        rawF_perc_peak_run_cells(c,p,:) = mean(rawF_peak_run_p); %average across peaks
    end
end
% if a cell has less than 10 spikes during running, len_p will be 0, should
% delete that cell before averaging
rawF_perc_peak_run_cells(all(all(isnan(rawF_perc_peak_run_cells),3),2),:,:)=[];
rawF_perc_peak_run = squeeze(mean(rawF_perc_peak_run_cells)); %average across cells
ste_rawF_perc_peak_run = squeeze(std(rawF_perc_peak_run_cells)/sqrt(size(rawF_perc_peak_run_cells,1)));

[n1 n2] = subplotn(10);
Caevent = figure;
x = (1:30)/30;
peak_amp_stay = zeros(1,10);
peak_amp_run = zeros(1,10);
for i = 1:10
    peak_amp_stay(i) = max(rawF_perc_peak_stay(i,:))-mean(rawF_perc_peak_stay(i,1:10));
    peak_amp_run(i) = max(rawF_perc_peak_run(i,:))-mean(rawF_perc_peak_run(i,1:10));
    subplot(n1,n2,i);
    errorbar(x,rawF_perc_peak_stay(i,:),ste_rawF_perc_peak_stay(i,:),'.','LineStyle','-','linewidth', 1.25,'MarkerSize',12); hold on;
    %text(0.35,2200, ['peakAmpStay=' num2str(peak_amp_stay(i))]);
    errorbar(x,rawF_perc_peak_run(i,:),ste_rawF_perc_peak_run(i,:),'.','LineStyle','-','linewidth', 1.25,'MarkerSize',12); hold on;
    %text(0.35,2800, ['peakAmpRun=' num2str(peak_amp_run(i))]);
    title(['top' num2str((i-1)*10) '-' num2str(10*i) '%']);
end
supertitle(sessions);
saveas(Caevent,[image_analysis_dest_deconv_sep sessions '_peak_Amp_every_10perc_rawF']);
save([image_analysis_dest_deconv_sep sessions '_spk_deconvolve_staynrun_seperate.mat' ],...
    'rawF_perc_peak_stay_cells','rawF_perc_peak_stay','ste_rawF_perc_peak_stay', ...
    'rawF_perc_peak_run_cells','rawF_perc_peak_run','ste_rawF_perc_peak_run',...
    'peak_amp_stay','peak_amp_run','-append');

