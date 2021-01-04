clear;
analysis_out = 'Z:\2P_analysis\';
bdata_source = 'Z:\Data\behavior\RC\';

ttl =   1;
noreg = 0;
imgreglaseron = 0;

sessions = '201216_img1078';% for behavior data
mouse = '1078';
date = '201216';
run = '000';

fprintf([date ' ' mouse '\n'])
img_fn = [date '_img' mouse '\getTC_' run '\'];
filename = dir([analysis_out,img_fn '*' '_TCave.mat']);
load([analysis_out,img_fn, filename.name]);
nIC = size(tc_avg,2);

decon_out = [date '_img' mouse '\deconvolution_' run '\'];
if ~exist(fullfile(analysis_out,decon_out))
    mkdir(fullfile(analysis_out, decon_out))
end

%% Deconvolution
% find the baseline states: 4s-120frames before cTargetOn
baseline = zeros(1,120*length(cTargetOn_cutted));
a = 1;
b = 120;
for i = 1:length(cTargetOn_cutted)
    baseline(a:b) = cTargetOn_cutted(i)-119:cTargetOn_cutted(i);
    a = a+120;
    b = b+120;
end

TCave = tc_avg;
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
save([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_spk_deconvolve_threshold' num2str(threshold) '.mat' ],'threshold',...
    'FRstay_cell', 'aveFR','num_spks_cell','options','spk_logic','spk','kernel','spk_peak','spk_inx');

%% plots
% write a GUI with subplots: TCave with red dots and kernel
hist_FR = figure;
histogram(FRstay_cell,'BinWidth',0.1);
title([date ' img' mouse 'deconvolution-firing rate-stationary-threshold' num2str(threshold)]);
saveas(hist_FR,[analysis_out decon_out date '_img' mouse '_histFR_deconvolution_threshold' num2str(threshold) '.fig']);

[fig_deconvolve] = GUI_rawTrace_nDeconvolve(TCave,kernel,spk_logic,[date '_img' mouse],threshold);
savefig([analysis_out decon_out date '_img' mouse '_GUI_TCave_deconvolution_threshold' num2str(threshold) '.fig']);

% spike logic
% plot spike logic
figure;
colormap gray;
imagesc(imcomplement(spk_logic'));
xlabel('Frame #');
ylabel('Cell #');
%xlim([25000 28992]);
title ([date ' img' mouse ' deconvolution']);
savefig([analysis_out decon_out date '_img' mouse '_scatter_spikeLogic_deconvolve_threshold' num2str(threshold) '.fig']);

% plot spike rates for each cell during stationary
figure;
scatter(1:length(FRstay_cell), FRstay_cell,70, 'ko');
ylabel('Spike Rate (Hz)');
xlabel('Cell #'); axis square;
%ylim([0 1.4]);
title ([date ' img' mouse ' deconvolution']);
savefig([analysis_out decon_out date '_img' mouse '_scatter_spikeRate_deconvolve_threshold' num2str(threshold) '.fig']);


%% delete the bad cells
nobase_maxperiod = zeros(1,size(spk_peak,2));
for c = 1:size(spk_peak,2)
    nobase_maxperiod(c) = max(diff(find(kernel(:,c)<100)));
    %kernel values smaller than 100 is considered as basline. diff(find) gives you how many continous frames the denoised signal doesn't return to baseline
end
figure; plot(nobase_maxperiod);
title('not return to baseline period');
savefig([analysis_out decon_out date '_img' mouse '_Caevent_nobase' num2str(threshold) '.fig']);
%if the maximum not return to baseline period is bigger than n # of frames,
%throw out this cell
nobase_thres = 8000;
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
baseline_thres = 3000;
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
imshow([analysis_out,img_fn 'AVG_' date '_img' mouse '_' run '_rgstr_tiff_every50_ref2_jpeg.jpg']); hold on;
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


