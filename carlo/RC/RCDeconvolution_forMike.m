clear; close all;
analysis_out = 'A:\home\carlo\mikeAnalysis\2P\';
bdata_source = 'A:\home\mike\Data\Behavior\';

ttl =   1;
noreg = 0;
imgreglaseron = 0;

RCExptListMike_Inter %ORDER FOR EXPT LIST:: 1507,1511,1510,1512,1513,1516,1520
pairings = loadRCList_Mike;
mouse = strtrim(expt(pairings{1,1}).mouse);
date = expt(pairings{1,1}).date;
run = expt(pairings{1,1}).run;
sessions = [date '_img' mouse];
region = '_';
fprintf([date ' ' mouse '\n'])
img_fn = [date '_img' mouse '\getTC_' run '\'];
filename = dir([analysis_out img_fn '*' '_TCave.mat']);
load([analysis_out,img_fn, filename.name]);
nIC = size(tc_avg,2);

decon_out = [date '_img' mouse '\deconvolution_' run '\'];
if ~exist(fullfile(analysis_out,decon_out))
    mkdir(fullfile(analysis_out, decon_out))
end

%% Deconvolution
% find the baseline states: 4s-120frames before cTargetOn
baseline = zeros(1,120*length(cTargetOn_cutted(1,2:end)));
a = 1;
b = 120;
for i = 1:length(cTargetOn_cutted(1,2:end))
    baseline(a:b) = cTargetOn_cutted(i+1)-119:cTargetOn_cutted(i+1);
    a = a+120;
    b = b+120;
end

preFilterFst200Frms = nanmean(raw_tc_avg(1:200,:),1); %take the first 200 frames from raw trace
TCAvg = tc_avg+preFilterFst200Frms; %add that value to the filtered TC to return to its original mean (in an attempt to make the filtered data AS comparable to the raw trace as possible, since in an ideal world they would be identical)
TCAvg_All=NaN(length(pockel_tc),size(raw_tc_avg,2));
TCAvg_All(laseron,:) = tc_avg_all(laseron,:)+preFilterFst200Frms; %add that value to the filtered TC to return to its original mean (in an attempt to make the filtered data AS comparable to the raw trace as possible, since in an ideal world they would be identical)
%plot TCAvg to make sure the only change was a shift returning to where the raw imaging data began
figure; plot(nanmean(tc_avg,2),'r'); hold on; plot(nanmean(TCAvg_All,2),'m'); plot(nanmean(raw_tc_avg,2),'k'); plot(nanmean(TCAvg,2),'g');
%should NOT expect TCAvg_All trace to align with the 3 others along the
%x-axis, since this TC contains NaN for every ITI. Even if you zoom in,
%this trace will not align with the others, it is here solely to check and
%ensure it's y-axis is where it should be after the processing above (i.e.
%in line with the first 200 frames of raw_tc_avg)

TCave = TCAvg;
frames = 1:1:size(TCave,1);
thresholdDeco = -4; %don't change this value across animals
%%this threshold simply sets the degree to which spikes are selected from
%%the TC, lower threshold means more spikes will be found throughout a TC
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
        'optimize_b', true, 'method','foopsi', 'smin', thresholdDeco); %since optimize baseline is true, the actuatual deconvolution takes raw data - baseline (denoised trace)
    % get only the peaks of each spike
    [spk_peak{c},spk_inx{c}] = findpeaks(spk(:,c));
    % spike logic
    spk_logic(:,c) = (ismember(frames,spk_inx{c}))';
    num_spks_cell(c) = sum(spk_logic(baseline,c)==1);
    FRstay_cell(c)= num_spks_cell(c)/length(baseline)*30; % firing rate = # of spikes/duration(s)
end
aveFR = mean(FRstay_cell);
%save variables
save([analysis_out date '_img' mouse '\' date '_img' mouse '_spk_deconvolve_threshold' num2str(thresholdDeco) '.mat' ],'thresholdDeco',...
    'FRstay_cell', 'aveFR','num_spks_cell','options','spk_logic','spk','kernel','spk_peak','spk_inx');

%% plots
% write a GUI with subplots: TCave with red dots and kernel
hist_FR = figure;
histogram(FRstay_cell,'BinWidth',0.1);
title([date ' img' mouse 'deconvolution-firing rate-stationary-threshold' num2str(thresholdDeco)]);
saveas(hist_FR,[analysis_out decon_out date '_img' mouse '_' region '_histFR_deconvolution_threshold' num2str(thresholdDeco) '.fig']);

%raw traces for each neuron
[fig_deconvolve] = GUI_rawTrace_nDeconvolve(TCave,kernel,spk_logic,[date '_img' mouse],thresholdDeco);
savefig([analysis_out decon_out date '_img' mouse '_' region '_GUI_TCave_deconvolution_threshold' num2str(thresholdDeco) '.fig']);

% spike logic
% plot spike logic - raster-like plot [cell by nFrames]
figure;
colormap gray;
imagesc(imcomplement(spk_logic'));
xlabel('Frame #');
ylabel('Cell #');
%xlim([25000 28992]);
title ([date ' img' mouse ' deconvolution']);
savefig([analysis_out decon_out date '_img' mouse '_' region '_scatter_spikeLogic_deconvolve_threshold' num2str(thresholdDeco) '.fig']);

% plot spike rates for each cell during stationary - baseline spike frequency
figure;
scatter(1:length(FRstay_cell), FRstay_cell,70, 'ko');
ylabel('Spike Rate (Hz)');
xlabel('Cell #'); axis square;
%ylim([0 1.4]);
title ([date ' img' mouse ' deconvolution']);
savefig([analysis_out decon_out date '_img' mouse '_' region '_scatter_spikeRate_deconvolve_threshold' num2str(thresholdDeco) '.fig']);


%% delete the bad cells
nobase_maxperiod = zeros(1,size(spk_peak,2));
for c = 1:size(spk_peak,2)
    nobase_maxperiod(c) = max(diff(find(kernel(:,c)<100)));
    %kernel values smaller than 100 is considered as basline. diff(find) gives you how many continous frames the denoised signal doesn't return to baseline
end
figure; plot(nobase_maxperiod);
title('not return to baseline period');
savefig([analysis_out decon_out date '_img' mouse '_' region '_Caevent_nobase' num2str(thresholdDeco) '.fig']);
%if the maximum not return to baseline period is bigger than n # of frames,
%throw out this cell
nobase_thres = 1000; %33 second threshold for return to baseline
nobase_cell = find(nobase_maxperiod > nobase_thres);

baseline_btm = zeros(1, size(TCave,2));
for n = 1:size(TCave,2)                     % for each cell
    TCave_sort = sort(TCave(:,n),1,'ascend');
    btm = TCave_sort(1:floor(length(TCave_sort)*0.1));
    baseline_btm(n) = mean(btm);
end
figure; plot(baseline_btm);
title('baseline fluorescence');
savefig([analysis_out decon_out date '_img' mouse '_' region '_baselineF' num2str(thresholdDeco) '.fig']);
baseline_thres = input('baseline: ');
%baseline_thres = [];
low_Fcell = find(baseline_btm < baseline_thres);
%low_Fcell = [];

lowFR = FRstay_cell<0.075;
for f = 1:size(TCave,2)
    if lowFR(1,f) == 1
        lowFRcell{1,f} = f;
    elseif sum(lowFR)==0
        lowFRcell=[];
    end
end
lowFRcell = cell2mat(lowFRcell); % 59.97.
    
badPCs = union(low_Fcell,lowFRcell);
TCave_cl = TCave; TCave_cl(:,badPCs) = [];
all_TCave_cl = TCAvg_All; all_TCave_cl(:,badPCs) = [];
FRstay_cell_cl = FRstay_cell;FRstay_cell_cl(badPCs) = [];
aveFR_cl = mean(FRstay_cell_cl);
spk_logic_cl = spk_logic; spk_logic_cl(:,badPCs) = [];
spk_cl = spk; spk_cl(:,badPCs) = [];
kernel_cl = kernel; kernel_cl(:,badPCs) = [];
spk_peak_cl = spk_peak; spk_peak_cl(badPCs) = [];
spk_inx_cl = spk_inx; spk_inx_cl(badPCs) = [];
allcells = 1:size(TCave,2);
goodcells = setdiff(allcells,badPCs);

% make sure you save everything you want to save
load([analysis_out date '_img' mouse '\getTC_' run '\' date '_img' mouse '_' run  'saveOutputs.mat']);
load([analysis_out date '_img' mouse '\getTC_' run '\' date '_img' mouse '_' run '_nPCA' num2str(nPCA) '_mu' num2str(mu) '_nIC' num2str(nIC) '_thresh' num2str(cluster_threshold) '_coor' num2str(threshold) run 'mask3D_final.mat'], 'mask3D');
mask3D_final = mask3D; mask3D_final(:,:,badPCs) = [];
save([analysis_out date '_img' mouse '\getTC_' run '\' date '_img' mouse '_' run '_nPCA' num2str(nPCA) '_mu' num2str(mu) '_nIC' num2str(nIC) '_thresh' num2str(cluster_threshold) '_coor' num2str(threshold) run 'mask3D_final.mat'], 'mask3D_final','-append');
save([analysis_out date '_img' mouse '\' date '_img' mouse '_deconvolution_thresh', num2str(thresholdDeco), '_TCave_cl.mat'],...
    'TCave_cl','thresholdDeco','badPCs','nobase_cell','low_Fcell','nobase_thres','baseline_thres','all_TCave_cl','goodcells');
save([analysis_out date '_img' mouse '\' date '_img' mouse '_spk_deconvolve_threshold' num2str(thresholdDeco) '.mat' ],'FRstay_cell_cl', ...
    'aveFR_cl','badPCs','spk_logic_cl','spk_cl','kernel_cl','spk_peak_cl',...
    'spk_inx_cl','badPCs','nobase_cell','low_Fcell','nobase_thres',...
    'baseline_thres','-append');

%% make a new mask figure
% load mask3D final
filename2 = dir([analysis_out,img_fn '*' 'thresh' num2str(cluster_threshold) '_coor' num2str(threshold) '_mask3D.mat']);
mask3D = load([analysis_out img_fn filename2.name]);
mask3D = mask3D.mask3D;

figure;
%needs to modify the line below due to different file names
imshow([analysis_out,img_fn 'AVG_' date '_img' mouse '_' run '_rgstr_tiff_every50_ref' num2str(ref) '_jpeg.jpg']); hold on;
for m  = 1:length(goodcells)
    bound = cell2mat(bwboundaries(mask3D(:,:,goodcells(m))));
    randcolor = rand(1,4);
    plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;
    text(mean(bound(:,2)),mean(bound(:,1)), ...
        num2str(goodcells(m)), 'color', 'y', 'FontSize', 8);
end
hold on;
title([date ' img' mouse ' masks only good cells']);
savefig([analysis_out date '_img' mouse '\' date '_img' mouse '_mask_wdendrites_goodcells.fig']);

figure;
%needs to modify the line below due to different file names
imshow([analysis_out,img_fn 'AVG_' date '_img' mouse '_' run '_rgstr_tiff_every50_ref' num2str(ref) '_jpeg.jpg']); hold on;
for m  = 1:length(badPCs)
    bound = cell2mat(bwboundaries(mask3D(:,:,badPCs(m))));
    randcolor = rand(1,4);
    plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;
    text(mean(bound(:,2)),mean(bound(:,1)), ...
        num2str(badPCs(m)), 'color', 'y', 'FontSize', 8);
end
hold on;
title([date ' img' mouse ' masks only bad cells']);
savefig([analysis_out date '_img' mouse '\' date '_img' mouse '_mask_wdendrites_badcells.fig']);



