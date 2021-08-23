clear;
analysis_out = 'Z:\2P_analysis\';
bdata_source = 'Z:\Data\behavior\RC\';

ttl =   1;
noreg = 0;
imgreglaseron = 0;

mouse = '1079';
date = '201114';
run = '000';

fprintf([date ' ' mouse '\n'])
img_fn = [date '_img' mouse '\getTC_' run '\'];
filename = dir([analysis_out,img_fn '*' '_TCave.mat']);
load([analysis_out,img_fn, filename.name]);
nIC = size(tc_avg,2);
figure; plot(mean(tc_avg,2));

filter_out = [date '_img' mouse '\HPfilter_' run '\'];
if ~exist(fullfile(analysis_out,filter_out))
    mkdir(fullfile(analysis_out, filter_out))
end

decon_out = [date '_img' mouse '\deconvolution_' run '\'];
if ~exist(fullfile(analysis_out,decon_out))
    mkdir(fullfile(analysis_out, decon_out))
end

%% high pass filtering
% butterworth filter: a type of signal processing filter designed to have a
% frequency response as flat as possible in the passband, also called
% maximally flat magnitude filter. 
FrameRate = 30;
cutoffF =[0.005, 0.01, 0.02, 0.035, 0.05];
filter_output = {};
for i = 1:length(cutoffF)
    [b,a]=butter(3,(cutoffF(i)*2/FrameRate),'high');%get transfer function coefficients. first input is filter order, second input is cutoff frequency, formula is on matlab page
    tc_avg_filtered = filtfilt(b,a,tc_avg); %zero-phase filtering, in both forward and reverse directions. this function operates along the first array dimension of x with size greater than 1
    
    figure; plot(mean(tc_avg,2)); hold on; plot(mean(tc_avg_filtered,2));
    title ([date ' img' mouse ' ' run ' ' num2str(cutoffF(i))]); ylabel('rawF(blue) filtered F(red)');
    xlabel ('frame');
    box off;
    savefig([analysis_out filter_out date '_img' mouse '_' run '_TCaveVsHPfilter' num2str(cutoffF(i)) '.fig']);
    
    % for each cell, find peaks in raw F and filtered F. if the number of peaks
    % are identical, then high pass filtering didn't get rid of any peaks
    npks_raw = zeros(1,size(tc_avg,2));
    npks_filtered = zeros(1,size(tc_avg,2));
    for c = 1:size(tc_avg,2)
        [pks1,locs1] = findpeaks(tc_avg(:,c));
        npks_raw(c) = length(locs1);
        
        [pks2,locs2] = findpeaks(tc_avg_filtered(:,c));
        npks_filtered(c) = length(locs2);
    end
    
    diff_pks = npks_filtered - npks_raw;
    figure; bar(diff_pks);
    ylim([-10 8]);xlabel('cell #');
    box off;
    title([date ' img' mouse ' ' run ' ' num2str(cutoffF(i)) ' npksHP - npksRaw']);
    savefig([analysis_out filter_out date '_img' mouse '_' run '_diff_pks' num2str(cutoffF(i)) '.fig']);
    filter_output{i}.cutoff = cutoffF(i);
    filter_output{i}.tc_avg_filtered = tc_avg_filtered;
    filter_output{i}.diff_pks = diff_pks;
end
save([analysis_out filter_out date '_img' mouse '_' run '_HPfilter.mat' ],'filter_output');
close all;
    
%% plot on a cell by cell basis, with different filtered frequencies
filename = dir([analysis_out filter_out '*' '_HPfilter.mat']);
load([analysis_out,filter_out, filename.name]);
[fig_filtered] = GUI_rawTrace_nHPfiltered(tc_avg,filter_output,[date '_img' mouse]);
savefig([analysis_out filter_out date '_img' mouse '_' run '_GUI_TCave_HPfilter.fig']);
%[fig_filtered_peaks] = GUI_rawTrace_nHPfiltered_peaks(tc_avg,filter_output,[date '_img' mouse]);

%% Deconvolution
% find the baseline states: 4s-120frames before cTargetOn
filename = dir([analysis_out filter_out '*' '_HPfilter.mat']);
load([analysis_out,filter_out, filename.name]);

baseline = zeros(1,120*length(cTargetOn_cutted));
a = 1;
b = 60;
for i = 1:length(cTargetOn_cutted)
    if cTargetOn_cutted(i)-59>0
        baseline(a:b) = cTargetOn_cutted(i)-59:cTargetOn_cutted(i);
        a = a+60;
        b = b+60;
    end
end
baseline(baseline==0) = []; % this should only happen when you manually cut the neural data

% 0.005, 0.01, 0.02, 0.035, 0.05
tc_avg_filtered = filter_output{4}.tc_avg_filtered;
cutoffF = filter_output{4}.cutoff;
TCave = tc_avg_filtered+abs(min(min(tc_avg_filtered)))+10; %bring all the filtered trace to above 0 
frames = 1:1:size(TCave,1);
threshold = -3; %changing the threhold basically changes the identification of spikes when the peak amplitude is small (those small peaks), doesn't change anything with the bigger jittered ones. -3 or -3.5 gives a FR close to 1 during stationary
% deconvolution
kernel = zeros(size(TCave,1),size(TCave,2));
spk = zeros(size(TCave,1),size(TCave,2));
spk_peak = {};
spk_inx = {};
spk_logic = zeros(size(TCave,1),size(TCave,2));
num_spks_cell = zeros(1,size(TCave,2));
FRstay_cell = zeros(1,size(TCave,2));
for c = 1: size(TCave,2)
    [kernel(:,c), spk(:,c), options] = deconvolveCa(TCave(:,c), 'optimize_pars', true, ...
        'optimize_b', true, 'method','foopsi', 'smin', threshold); %since optimize baseline is true, the actuatual deconvolution takes raw data - baseline (denoised trace)
    % get only the peaks of each spike
    [spk_peak{c},spk_inx{c}] = findpeaks(spk(:,c));
    % spike logic
    spk_logic(:,c) = (ismember(frames,spk_inx{c}))';
    num_spks_cell(c) = sum(spk_logic(baseline,c)==1);
    FRstay_cell(c)= num_spks_cell(c)/length(baseline)*30; % firing rate = # of spikes/duration(s)
end
aveFR = mean(FRstay_cell);
%save variables
save([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_spk_cutoff_' num2str(cutoffF) '_deconvolve_threshold' num2str(threshold) '.mat' ],'threshold',...
    'cutoffF','FRstay_cell', 'aveFR','num_spks_cell','options','spk_logic','spk','kernel','spk_peak','spk_inx');


% plots
% write a GUI with subplots: TCave with red dots and kernel
hist_FR = figure;
histogram(FRstay_cell,'BinWidth',0.1);
title([date ' img' mouse 'deconvolution-firing rate-stationary-threshold' num2str(threshold) 'cutoffF ' num2str(cutoffF)]);
saveas(hist_FR,[analysis_out decon_out date '_img' mouse '_histFR_deconvolution_threshold' num2str(threshold) '_cutoffF ' num2str(cutoffF) '.fig']);

[fig_deconvolve] = GUI_rawTrace_nDeconvolve(TCave,kernel,spk_logic,[date '_img' mouse],threshold);
savefig([analysis_out decon_out date '_img' mouse '_GUI_TCave_deconvolution_threshold' num2str(threshold) '_cutoffF ' num2str(cutoffF) '.fig']);

% spike logic
% plot spike logic
figure;
colormap gray;
imagesc(imcomplement(spk_logic'));
xlabel('Frame #');
ylabel('Cell #');
%xlim([25000 28992]);
title ([date ' img' mouse ' deconvolution']);
savefig([analysis_out decon_out date '_img' mouse '_scatter_spikeLogic_deconvolve_threshold' num2str(threshold) '_cutoffF ' num2str(cutoffF) '.fig']);

% plot spike rates for each cell during stationary
figure;
scatter(1:length(FRstay_cell), FRstay_cell,70, 'ko');
ylabel('Spike Rate (Hz)');
xlabel('Cell #'); axis square;
%ylim([0 1.4]);
title ([date ' img' mouse ' deconvolution']);
savefig([analysis_out decon_out date '_img' mouse '_scatter_spikeRate_deconvolve_threshold' num2str(threshold) '_cutoffF ' num2str(cutoffF) '.fig']);


%% delete the bad cells
nobase_maxperiod = zeros(1,size(spk_peak,2));
for c = 1:size(spk_peak,2)
    if sum(kernel(:,c)<100)>0
        nobase_maxperiod(c) = max(diff(find(kernel(:,c)<100)));% got a first cell that doesn't have a value lower than 100 after 3 years
    else 
        nobase_maxperiod(c) = 5000;% the cell doesn't have a value lower than 100 
    end
    %kernel values smaller than 100 is considered as basline. diff(find) gives you how many continous frames the denoised signal doesn't return to baseline
end
figure; plot(nobase_maxperiod);
title('not return to baseline period');
savefig([analysis_out decon_out date '_img' mouse '_Caevent_nobase' num2str(threshold) '.fig']);
%if the maximum not return to baseline period is bigger than n # of frames,
%throw out this cell
nobase_thres = 2000; % determine this number after looking at the plot. sort out the outliers, I typically use a number in between 2000-4000. 
nobase_cell = find(nobase_maxperiod > nobase_thres);

baseline_btm = zeros(1, size(TCave,2));
for n = 1:size(TCave,2)                     % for each cell
    TCave_sort = sort(TCave(:,n),1,'ascend');
    btm = TCave_sort(1:floor(length(TCave_sort)*0.1));
    baseline_btm(n) = mean(btm);
end
figure; plot(baseline_btm);
title('baseline fluorescence');
savefig([analysis_out decon_out date '_img' mouse '_baselineF' num2str(threshold) '.fig']);
baseline_thres = 2000; % same as above, if all of the cells look bright to you, can also set this to empty. 
%baseline_thres = [];
low_Fcell = find(baseline_btm < baseline_thres);
%low_Fcell = [];

badPCs = union(nobase_cell,low_Fcell);
TCave_cl = TCave;TCave_cl(:,badPCs) = [];
all_TCave_cl = tc_avg_all; all_TCave_cl(:,badPCs) = [];
FRstay_cell_cl = FRstay_cell;FRstay_cell_cl(badPCs) = [];
aveFR_cl = mean(FRstay_cell_cl);
spk_logic_cl = spk_logic; spk_logic_cl(:,badPCs) = [];
spk_cl = spk; spk_cl(:,badPCs) = [];
kernel_cl = kernel; kernel_cl(:,badPCs) = [];
spk_peak_cl = spk_peak; spk_peak_cl(badPCs) = [];
spk_inx_cl = spk_inx; spk_inx_cl(badPCs) = [];

% make sure you save everything you want to save
save([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_deconvolution_thresh', num2str(threshold), '_TCave_cl.mat'],...
    'TCave_cl','threshold','badPCs','nobase_cell','low_Fcell','nobase_thres','baseline_thres','all_TCave_cl');
save([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_spk_deconvolve_threshold' num2str(threshold) '.mat' ],'FRstay_cell_cl', ...
    'aveFR_cl','badPCs','spk_logic_cl','spk_cl','kernel_cl','spk_peak_cl',...
    'spk_inx_cl','threshold','badPCs','nobase_cell','low_Fcell','nobase_thres',...
    'baseline_thres','-append');

%% make a new mask figure
% load mask3D final
filename2 = dir([analysis_out,img_fn '*' 'thresh97_coor0.8_mask3D.mat']);
mask3D = load([analysis_out img_fn filename2.name]);
mask3D = mask3D.mask3D;
allcells = 1:size(TCave,2);
goodcells = setdiff(allcells,badPCs);
figure;
%needs to modify the line below due to different file names
imshow([analysis_out,img_fn 'AVG_' date '_img' mouse '_' run '_rgstr_tiff_every50_ref1_jpeg.jpg']); hold on;
for m  = 1:length(goodcells)
    bound = cell2mat(bwboundaries(mask3D(:,:,goodcells(m))));
    randcolor = rand(1,4);
    plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;
    text(mean(bound(:,2)),mean(bound(:,1)), ...
        num2str(goodcells(m)), 'color', 'y', 'FontSize', 8);
end
hold on;
title([date ' img' mouse ' masks only good cells']);
savefig([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_mask_wdendrites_goodcells.fig']);



