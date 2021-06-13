% spike identification using deconvolution. (function deconvolve Ca got from OASIS on git hub)
% to use the deconvolution function, need to add OASIS function to search
% path: >>oasis_setup
% function deconvolve Ca: denoise raw fluorescence trace and identifies spike events. 
% right now using FOOPSI method and the model is ar1 (also tried ar2: see commented bottom part),ar1 works better than ar2.
% then SECTION II draws the GUI showing raw traces with identified events and deconvolved trace.

%% paths
clear;
% sessions = {'200305_img1049','200319_img1064_airpuff','200319_img1064_airpuff_2'};
% days = {'1049-200305_1','1064-200319_1','1064-200319_2'};
sessions = '200305_img1049'; 
days = '1049-200305_1';
image_analysis_base = 'Z:\Analysis\motorizedWheel_Analysis\airpuff\imaging_analysis\'; 
color_code = {'b','r','k','c'};

%% deconvolution and plot
% paths and load data-------------------------------------------------------------------------------------------------------
image_analysis_dest = [image_analysis_base, sessions, '\'];
deconvolutionDest = [image_analysis_dest,'deconv_wfastRunSmin\'];
if ~exist(deconvolutionDest)
    mkdir(deconvolutionDest);
end
% behavior analysis results
behav_dest = ['Z:\Analysis\motorizedWheel_Analysis\airpuff\behavioral_analysis\' days '\'];
%load data
filename = dir([image_analysis_dest 'getTC\' '*' '_TCave.mat']);
TCave = load([image_analysis_dest 'getTC\' filename.name]);
TCave = TCave.tc_avg;
behav_output = load([behav_dest days '_behavAnalysis.mat']);
stay_nopuff = behav_output.stay_nopuff;
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
% saveas(hist_FR,[deconvolutionDest sessions '_histFR_deconvolution_threshold' num2str(threshold) '.fig']);

%save variables
save([deconvolutionDest sessions '_spk_deconvolve_staynrun_seperate.mat' ],'threshold_run',...
   'smin_run', 'FRstay_cell', 'aveFR','options','spk_logic','spk','kernel','spk_peak','spk_inx');

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
nobase_thres = 1500;
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
baseline_thres = 2500;
low_Fcell = find(baseline_btm<baseline_thres);
%---------------------------------------------------------------------------------------------------------------------------------
% find the ones that has super low FR (fake cells), delete those in TCave, spk_logic,spk,kernel,spk_peak,spk_inx
%FRthres_low = 0.1; %Lisberger paper showed the lowest firing rate is 0.1
%badPCs1 = find(FRstay_cell<FRthres_low);
%I think I should use other features to filter out the bad cells instead of
%simply looking at the firing rate. if using the features below, should
%filter out the cells that are firing too much or too little (max no fire,
%max no baseline, max fire...). These features looks through the neural
%activity throughtout the whole time, whereas FR is just an average and a
%lot of detials can be missing

% use other criteria to sort out bad cells:
% peak_std = zeros(1,size(spk_peak,2));
% peak_mean = zeros(1,size(spk_peak,2));
% for c = 1:size(spk_peak,2)                  % for each cell
%     peak_std(c) = std(spk_peak{c});
%     peak_mean(c) = mean(spk_peak{c});
% end
% peak_variation = peak_std./peak_mean; %if peak variation is too big, this cell might not be good.
% figure; plot(peak_variation);
% title('peak variation');
% savefig([deconvolutionFigureDest sessions '_CaeventPeak_variation' num2str(threshold) '.fig']);
% peak_vartoomuch = find(peak_variation>0.64);
%peak_vartoomuch = [];
% cannot totally through out all of the cells in here, need to go back and
% look at each cell that has a high variation and decide.
% but this at least gives me a good pool to look at
% cell # 18,9,7,4,1 in session 191206_img1038 can all be sorted out using this way
% ----------------------------------------------------------------------------------------------
% another feature is that the deconvolved signal doesn't go back to
% baseline for a long time: cell#88,53,91,and 101

% %----------------------------------------------------------------------------------------------------
% % doesn't fire for a long time
% max_nofire = zeros(1,size(spk_peak,2));
% for c = 1:size(spk_peak,2)
%     max_nofire(c) = max(diff(spk_inx{c}));
% end
% figure; plot(max_nofire);
% title('not firing');
% savefig([deconvolutionFigureDest sessions '_Caevent_nofire' num2str(threshold) '.fig']);
% max_nofire_thres = 7000;
% nofire_cell = find(max_nofire > max_nofire_thres);
% %----------------------------------------------------------------------------------------------------
% %a third feature is that the cells are firing a lot 
% %calculate the firing rate of a period of time and if the FR is high for a
% %long time then through the cell out first calculate # of spikes per second
% maxfire = zeros(1,size(spk_peak,2));
% binwidth = 900; % number of frames
% bins = floor(size(spk_logic,1)/binwidth);
% for c = 1:size(spk_peak,2)
%     nfire = []; t = 1;
%     while t < size(spk_logic,1)-binwidth % use a band to scan through the whole session, count the number of total spikes in time length binwidth                               
%         nfire = [nfire sum(spk_logic(t:t+binwidth-1,c)==1)]; 
%         t = t+1;
%     end
%     maxfire(c) = max(nfire); % for each cell get the highest firing rate during any bindwidth
% end
% figure; plot(maxfire); title('firing too much');
% savefig([deconvolutionFigureDest sessions '_Caevent_maxfire' num2str(threshold) '.fig']);
% figure; hist(maxfire); title('histogram # of events in 900 frames');
% savefig([deconvolutionFigureDest sessions '_Caevent_maxfire_hist' num2str(threshold) '.fig']);
% FR_toohigh = find(maxfire>65);
%FR_toohigh = [];

% badcells2 = union(nofire_cell,nobase_cell);
% badPCs2 = union(badcells2,FR_toohigh);

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
    'aveFR_cl','badPCs','spk_logic_cl','spk_cl','kernel_cl','spk_peak_cl','spk_inx_cl',...
    'threshold','badPCs','nobase_cell','low_Fcell','nobase_thres','baseline_thres','-append');

%% make a new mask figure
% load mask3D final
filename2 = dir([image_analysis_dest 'getTC\' '*' 'thresh97_coor0.8_mask3D.mat']);
mask3D = load([image_analysis_dest 'getTC\' filename2.name]);
mask3D = mask3D.mask3D;
allcells = 1:size(TCave,2);
goodcells = setdiff(allcells,badPCs);
figure; imshow([image_analysis_dest 'getTC\' 'AVG_' sessions '_000_rgstr_tiff_1_29999_50_ref43_' 'jpeg.jpg']); hold on;
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