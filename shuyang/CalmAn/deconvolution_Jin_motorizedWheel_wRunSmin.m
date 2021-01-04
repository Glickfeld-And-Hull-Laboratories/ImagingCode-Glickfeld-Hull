% spike identification using deconvolution. (function deconvolve Ca got from OASIS on git hub)
% to use the deconvolution function, need to add OASIS function to search
% path: >>oasis_setup
% function deconvolve Ca: denoise raw fluorescence trace and identifies spike events. 
% right now using FOOPSI method and the model is ar1 (also tried ar2: see commented bottom part),ar1 works better than ar2.
% then SECTION II draws the GUI showing raw traces with identified events and deconvolved trace.

%% paths
clear;
% sessions = {'200116_img1041','200217_img1061','200225_img1049','200319_img1064','200319_img1064_2'};
sessions = '200319_img1064_2'; 
%days = {'1041-200116_1','1061-200217_1','1049-200225_1','1064-200319_1','1064-200319_2'};
days = '1064-200319_2';
image_analysis_base = 'Z:\Analysis\motorizedWheel_Analysis\running\imaging_analysis\'; 
color_code = {'b','r','k','c'};

%% deconvolution and plot
% paths and load data-------------------------------------------------------------------------------------------------------
image_analysis_dest = [image_analysis_base, sessions, '\'];
deconvolutionDest = [image_analysis_dest,'deconv_wfastRunSmin\'];
if ~exist(deconvolutionDest)
    mkdir(deconvolutionDest);
end
% behavior analysis results
behav_dest = ['Z:\Analysis\motorizedWheel_Analysis\running\behavioral_analysis\' days '\'];
%load data
filename = dir([image_analysis_dest 'getTC\' '*' '_TCave.mat']);
rawF_output = load([image_analysis_dest sessions, '_deconvolution_thresh-4_TCave_cl.mat']);
TCave = rawF_output.TCave_cl;
behav_output = load([behav_dest days '_behavAnalysis.mat']);
frm_stay = behav_output.stay_vec;
frm_fastrun = behav_output.run_fast_vec;

%Deconvolution------------------------------------------------------------------------------------------------------------------
% input: fluorescence trace: T*1, threshold for spike size (how big the spikes are supposed to be: normal set as -3(3 sigma, 3 times bigger than noise levels)).
% kernal: denoised trace, T*1 vector
% spikes: all of the values in a spike(thus, you will get more than 1 values if a spike is longer than a frame. e.g.: if a spike is 120ms long, you will get 4 values back to back which all belong to a single spike)

% deconvolution on fastest speed and get smin
frames = 1:1:size(TCave,1);
threshold_run = -4; %changing the threhold basically changes the identification of spikes when the peak amplitude is small (those small peaks), doesn't change anything with the bigger jittered ones. -3 or -3.5 gives a FR close to 1 during stationary
kernel_run = zeros(length(frm_fastrun),size(TCave,2));
spk_run = zeros(length(frm_fastrun),size(TCave,2));
smin_run = zeros(1,size(TCave,2)); % smin is the threshold (in fluorescence values) of smallest spikes being detected
for c = 1: size(TCave,2)
    [kernel_run(:,c), spk_run(:,c), options] = deconvolveCa(TCave(frm_fastrun,c), 'optimize_pars', true, ...
        'optimize_b', true, 'method','foopsi', 'smin', threshold_run);
    smin_run(c) = options.smin;
end

% deconvolution for whole experiment
kernel = zeros(size(TCave,1),size(TCave,2));
spk = zeros(size(TCave,1),size(TCave,2));
spk_peak = {};
spk_inx = {};
spk_logic = zeros(size(TCave,1),size(TCave,2));
num_spks_cell = zeros(1,size(TCave,2));
FRstay_cell = zeros(1,size(TCave,2));
for c = 1: size(TCave,2)
    [kernel(:,c), spk(:,c), options] = deconvolveCa_stay(TCave(:,c), 'optimize_pars', true, ...
        'optimize_b', true, 'method','foopsi', 'smin', smin_run(c)); % deconvolveCa_stay is the function I modified that applies the smin of running to the whole dataset (only changed for FOOPSI model OR1)
    % get only the peaks of each spike
    [spk_peak{c},spk_inx{c}] = findpeaks(spk(:,c));
    % spike logic
    spk_logic(:,c) = (ismember(frames,spk_inx{c}))';
    num_spks_cell(c) = sum(spk_logic(frm_stay,c)==1);
    FRstay_cell(c)= num_spks_cell(c)/length(frm_stay)*30; % firing rate = # of spikes/duration(s)
end
aveFR = mean(FRstay_cell);
% hist_FR = figure;
% histogram(FRstay_cell,'BinWidth',0.1);
% title([sessions 'deconvolution-firing rate-stationary-threshold' num2str(threshold)]);
% saveas(hist_FR,[deconvolutionFigureDest sessions '_histFR_deconvolution_threshold' num2str(threshold) '.fig']);

%save variables
save([deconvolutionDest sessions '_spk_deconvolve_staynrun_seperate.mat'],'threshold_run',...
   'smin_run', 'FRstay_cell', 'aveFR','options','spk_logic','spk','kernel','spk_peak','spk_inx');


%% 
%{
% plots---------------------------------------------------------------------------------------------------------------------------
% write a GUI with subplots: TCave with red dots and kernel
[fig_deconvolve] = GUI_rawTrace_nDeconvolve(TCave,kernel,spk_logic,sessions,threshold);
savefig([deconvolutionFigureDest sessions '_GUI_TCave_deconvolution_threshold' num2str(threshold) '.fig']);

% spike logic
% plot spike logic
figure;
colormap gray;
imagesc(imcomplement(spk_logic'));
xlabel('Frame #');
ylabel('Cell #');
%xlim([25000 28992]);
title ([sessions ' deconvolution']);
%print([image_analysis_dest sessions{i} '_scatter_spikeLogic_deconvolve_threshold' num2str(threshold)],'-dpdf','-fillpage');
% this figure is before deleting the bad cells!
savefig([deconvolutionFigureDest sessions '_scatter_spikeLogic_deconvolve_threshold' num2str(threshold) '.fig']);

% plot spike rates for each cell during stationary
figure;
scatter(1:length(FRstay_cell), FRstay_cell,70, 'ko');
ylabel('Spike Rate (Hz)');
xlabel('Cell #'); axis square;
ylim([0 1.4]);
title ([sessions ' deconvolution']);
savefig([deconvolutionFigureDest sessions '_scatter_spikeRate_deconvolve_threshold' num2str(threshold) '.fig']);


%% delete the bad cells
nobase_maxperiod = zeros(1,size(spk_peak,2));
for c = 1:size(spk_peak,2)
    nobase_maxperiod(c) = max(diff(find(kernel(:,c)<100)));
    %kernel values smaller than 100 is considered as basline. diff(find) gives you how many continous frames the denoised signal doesn't return to baseline
end
figure; plot(nobase_maxperiod);
title('not return to baseline period');
savefig([deconvolutionFigureDest sessions '_Caevent_nobase' num2str(threshold) '.fig']);
%if the maximum not return to baseline period is bigger than n # of frames,
%through out this cell
nobase_thres = 1000;
nobase_cell = find(nobase_maxperiod > nobase_thres);

baseline_btm = zeros(1, size(TCave,2));
for n = 1:size(TCave,2)                     % for each cell
        TCave_sort = sort(TCave(:,n),1,'ascend');
        btm = TCave_sort(1:floor(length(TCave_sort)*0.1));
        %        % test if btm always belongs to running
        %         btm_inx = frames(ismember(TCave_cl(:,n),btm));
        %         NbtmInRun = sum((ismember(frm_run,btm_inx))==1);
        baseline_btm(n) = mean(btm);
end
figure; plot(baseline_btm);
title('baseline fluorescence');
savefig([deconvolutionFigureDest sessions '_baselineF' num2str(threshold) '.fig']);
baseline_thres = 2000;
%baseline_thres = [];
low_Fcell = find(baseline_btm < baseline_thres);
%low_Fcell = [];

badPCs = union(nobase_cell,low_Fcell);
TCave_cl = TCave;TCave_cl(:,badPCs) = [];
FRstay_cell_cl = FRstay_cell;FRstay_cell_cl(badPCs) = [];
aveFR_cl = mean(FRstay_cell_cl);
spk_logic_cl = spk_logic; spk_logic_cl(:,badPCs) = [];
spk_cl = spk; spk_cl(:,badPCs) = [];
kernel_cl = kernel; kernel_cl(:,badPCs) = [];
spk_peak_cl = spk_peak; spk_peak_cl(badPCs) = [];
spk_inx_cl = spk_inx; spk_inx_cl(badPCs) = [];


% make sure you save everything you want to save
save([image_analysis_dest sessions '_deconvolution_thresh', num2str(threshold), '_TCave_cl.mat'],...
    'TCave_cl','threshold','badPCs','nobase_cell','low_Fcell','nobase_thres','baseline_thres');
save([image_analysis_dest sessions '_spk_deconvolve_threshold' num2str(threshold) '.mat' ],'FRstay_cell_cl', ...
    'aveFR_cl','badPCs','spk_logic_cl','spk_cl','kernel_cl','spk_peak_cl',...
    'spk_inx_cl','threshold','badPCs','nobase_cell','low_Fcell','nobase_thres',...
    'baseline_thres','-append');

%% make a new mask figure
% load mask3D final
filename2 = dir([image_analysis_dest 'getTC\' '*' 'thresh97_coor0.8_mask3D.mat']);
mask3D = load([image_analysis_dest 'getTC\' filename2.name]);
mask3D = mask3D.mask3D;
allcells = 1:size(TCave,2);
goodcells = setdiff(allcells,badPCs);
figure; imshow([image_analysis_dest 'getTC\' 'AVG_' sessions '_000_rgstr_tiff_1_29999_50_ref50_' 'jpeg.jpg']); hold on;
for m  = 1:length(goodcells)
    bound = cell2mat(bwboundaries(mask3D(:,:,goodcells(m))));
    randcolor = rand(1,4);
    plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;
    text(mean(bound(:,2)),mean(bound(:,1)), ...
        num2str(goodcells(m)), 'color', 'y', 'FontSize', 8);
end
hold on;
title([sessions ' masks only good cells']);
savefig([image_analysis_dest sessions '_mask_wdendrites_goodcells.fig']);



%%
% frames = 1:1:60000;
% % what pengcheng told me:
% [c1, s1, options1] = deconvolveCa(TCave(:,1), 'optimize_pars', true, ...
%     'optimize_b', true, 'method','foopsi', 'smin', -3.5);
% % get only the peaks of each spike 
% [peaks1,locs1] = findpeaks(s1);
% % spike logic
% s_logic1 = ismember(frames,locs1);
% s_logic1 = s_logic1';
% num_spks1 = sum(s_logic1==1);
% 
% % don't know why ar2 doesn't work well, spike numbers are not reasonable,
% % 13 times bigger than using ar1.
% [c2, s2, options2] = deconvolveCa(TCave(:,1), 'foopsi','ar2',...
%     'optimize_pars', true, 'optimize_b', true,'smin', -4);
% [peaks2,locs2] = findpeaks(s2);
% s_logic2 = ismember(frames,locs2);
% s_logic2 = s_logic2';
% num_spks2 = sum(s_logic2==1);
% 
% [c3, s3, options3] = deconvolveCa(TCave(:,1), 'optimize_pars', true, ...
%     'optimize_b', true, 'method','foopsi');
% s_logic3 = logical(s3);
% num_spks3 = sum(s_logic3==1);
% 
% %%
% plotred = TCave(:,1).*s_logic2;
% plotred(plotred==0) = NaN;
% figure; 
% subplot(2,1,1)
% plot(TCave(:,1)); hold on;
% plot(plotred,'ro','MarkerSize',8); ylabel('TCave');
% xlim([200,3500]);
% subplot(2,1,2)
% plot(c2);
% xlim([200,3500]);ylabel('c');xlabel('frame');

%}