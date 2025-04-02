clear all; close all; clc

data_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\josh\2p\';
analysis_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\court\2P_Analysis\';
bdata_source = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\josh\Behavior\crp_data\';

exp = 1;
session = 2;

RCExptListJosh

mouse = expt(exp).mouse;
date = expt(exp).date{session};
date_b = [date(5:6), date(1:2), date(3:4)];
time = expt(exp).time{session};
run = expt(exp).run{session};

sessions = [date '_img' mouse];

thresholdDeco = -2; %don't change this value across animals
%%this threshold simply sets the degree to which spikes are selected from
%%the TC, lower threshold means more spikes will be found throughout a TC
if isempty(thresholdDeco)
img_fn = [date '_img' mouse '\getTC\'];
else
img_fn = [date '_img' mouse '\getTC_' num2str(thresholdDeco) '\'];    
end
fprintf([date ' ' mouse '\n'])
filename = dir([analysis_out,date '_img' mouse '\getTC\' '*' '_TCave.mat']);
load([analysis_out,date '_img' mouse '\getTC\', filename.name]);
% load cTargetOn;
b_data = load(fullfile(bdata_source,['data-i' mouse '-' date_b '-' time '.mat']));
mworks = b_data.input;
%
nIC = size(tc_avg,2);

decon_out = [date '_img' mouse '\deconvolution\'];
if ~exist(fullfile(analysis_out,decon_out))
    mkdir(fullfile(analysis_out, decon_out))
end

%% Deconvolution
% find the baseline states: 4s-120frames before cTargetOn
baseline = zeros(1,50*length(cTargetOn_cutted(1,2:end)));
a = 1;
b = 50;
for i = 1:length(cTargetOn_cutted)-1
    baseline(a:b) = cTargetOn_cutted(i+1)-49:cTargetOn_cutted(i+1);
    a = a+50;
    b = b+50;
end

start=1; finish=50; baselineTC = cell(1,size(raw_tc_avg,2));
for i = 1:length(cTargetOn_cutted)-1
    for ii = 1:size(raw_tc_avg,2)
baselineTC{ii} = [baselineTC{ii} raw_tc_avg(baseline(1,start:finish),ii)];
    end
    start = start+50; finish = finish+50;
end

for i = 1:size(raw_tc_avg,2)
    baselineTC{i} = mean(baselineTC{i},2,'omitnan');
end
baselineTC= cell2mat(baselineTC);
singleValue_baselineTC = mean(baselineTC,1,'omitnan');

tempFig = setFigure;
shadedErrorBar(1:size(baselineTC,1),mean(baselineTC,2,'omitnan'),std(baselineTC,[],2,'omitnan')./sqrt(size(raw_tc_avg,2)));


preFilterFst200Frms = mean(baselineTC,1); %take the first 50 frames from each trial  raw trace
TCAvg = tc_avg+preFilterFst200Frms; %add that value to the filtered TC to return to its original mean (in an attempt to make the filtered data AS comparable to the raw trace as possible, since in an ideal world they would be identical)
TCAvg_All=NaN(length(pockel_tc),size(tc_avg,2));
TCAvg_All(laseron,:) = tc_avg_all(laseron,:)+preFilterFst200Frms; %add that value to the filtered TC to return to its original mean (in an attempt to make the filtered data AS comparable to the raw trace as possible, since in an ideal world they would be identical)
%plot TCAvg to make sure the only change was a shift returning to where the raw imaging data began
tempFig=setFigure; hold on; plot(nanmean(tc_avg,2),'r'); hold on; plot(nanmean(TCAvg_All,2),'m'); plot(nanmean(raw_tc_avg,2),'k'); plot(nanmean(TCAvg,2),'g');
%should NOT expect TCAvg_All trace to align with the 3 others along the
%x-axis, since this TC contains NaN for every ITI. Even if you zoom in,
%this trace will not align with the others, it is here solely to check and
%ensure it's y-axis is where it should be after the processing above (i.e.
%in line with the first 200 frames of raw_tc_avg)

%% look at individual baseline activity for each cell, if certain cell's baseline are 2o above the population mean, remove them as corrupt 
%
tempFig = setFigure; hold on;
for i = 1:size(tc_avg,2)
baselineCell(i) = mean(tc_avg(baseline,i),1,'omitnan');

yline(baselineCell(i),'Color',rand(3,1));
hold on;
end

baselineCell_avg = mean(baselineCell,2);
baselineCell_std = std(baselineCell,[],2);
bleachCell_ind = baselineCell(baselineCell > (baselineCell_avg+(baselineCell_std*2)));
%}
%%
%TCAvg(:,bleachCell_ind)=[];

TCave = TCAvg;
frames = 1:1:size(TCave,1);
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
        'optimize_b', true, 'method','foopsi', 'smin', thresholdDeco); %since optimize baseline is true, the actual deconvolution takes raw data - baseline (denoised trace)
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
hist_FR = setFigure;
histogram(FRstay_cell,'BinWidth',0.1);
title([date ' img' mouse 'deconvolution-firing rate-stationary-threshold' num2str(thresholdDeco)]);
saveas(hist_FR,[analysis_out decon_out date '_img' mouse '_histFR_deconvolution_threshold' num2str(thresholdDeco) '.fig']);

%raw traces for each neuron
[fig_deconvolve] = GUI_rawTrace_nDeconvolve_Carlo(TCave,kernel,spk_logic,[date '_img' mouse],thresholdDeco,singleValue_baselineTC,zeros(size(singleValue_baselineTC)));
savefig([analysis_out decon_out date '_img' mouse '_GUI_TCave_deconvolution_threshold' num2str(thresholdDeco) '.fig']);

% spike logic
% plot spike logic - raster-like plot [cell by nFrames]
tempFig=setFigure; hold on;
colormap gray;
imagesc(imcomplement(spk_logic'));
xlabel('Frame #');
ylabel('Cell #');
%xlim([25000 28992]);
title ([date ' img' mouse ' deconvolution']);
savefig([analysis_out decon_out date '_img' mouse '_scatter_spikeLogic_deconvolve_threshold' num2str(thresholdDeco) '.fig']);

% plot spike rates for each cell during stationary - baseline spike frequency
tempFig=setFigure; hold on;
scatter(1:length(FRstay_cell), FRstay_cell,70, 'ko');
ylabel('Spike Rate (Hz)');
xlabel('Cell #'); axis square;
%ylim([0 1.4]);
title ([date ' img' mouse ' deconvolution']);
savefig([analysis_out decon_out date '_img' mouse '_scatter_spikeRate_deconvolve_threshold' num2str(thresholdDeco) '.fig']);

% finds baseline kernel values for first 50 frames of every trial for every cell
start=1; finish=50;baselineKernel = cell(1,size(raw_tc_avg,2));
for i = 1:length(cTargetOn_cutted)-1
    for ii = 1:size(raw_tc_avg,2)
baselineKernel{ii} = [baselineKernel{ii}; kernel(baseline(1,start:finish),ii)];
    end
    start = start+50; finish = finish+50;
end
for i = 1:size(raw_tc_avg,2)
    baselineKernel{i} = mean(baselineKernel{i},1,'omitnan');
end
baselineKernel= cell2mat(baselineKernel);


%% delete the bad cells
nobase_maxperiod = zeros(1,size(spk_peak,2));
for c = 1:size(spk_peak,2)
    if sum(find(kernel(:,c)<baselineKernel(1,c))) > 1
    nobase_maxperiod(c) = max(diff(find(kernel(:,c)<baselineKernel(1,c))));
    %added by Carlo - 11/08/21 - to correct error where the cell does not return to baseline [values lower than 100]
    elseif sum(find(kernel(:,c)<baselineKernel(1,c))) < 1
    nobase_maxperiod(c) = length(kernel(:,c));
    end
    %kernel values smaller than 100 is considered as basline. diff(find) gives you how many continous frames the denoised signal doesn't return to baseline
end
tempFig=setFigure; hold on; plot(nobase_maxperiod); hold on; hline(1000,'r');
title('not return to baseline period');
savefig([analysis_out decon_out date '_img' mouse '_Caevent_nobase' num2str(thresholdDeco) '.fig']);
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
tempFig=setFigure; hold on; plot(baseline_btm);
title('baseline fluorescence');
savefig([analysis_out decon_out date '_img' mouse '_baselineF' num2str(thresholdDeco) '.fig']);
baseline_thres = 1000;%str2double(input('Set the baseline: '));
%baseline_thres = [];
low_Fcell = find(baseline_btm < baseline_thres);
%low_Fcell = [];
low_Fcell = [low_Fcell nobase_cell];

lowFR = FRstay_cell<aveFR-(3*std(FRstay_cell,[],2));
for f = 1:size(TCave,2)
    if lowFR(1,f) == 1
        lowFRcell{1,f} = f;
    elseif sum(lowFR)==0
        lowFRcell=[];
    end
end
lowFRcell = cell2mat(lowFRcell); % 59.97.
lowFRcell = [lowFRcell find(FRstay_cell<0.05)];    
    
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
load([analysis_out date '_img' mouse '\getTC\' date '_img' mouse '_saveOutputs.mat']);
load([analysis_out date '_img' mouse '\getTC\' date ...
    '_img' mouse '_nPCA' num2str(nPCA) '_mu' num2str(mu) ...
    '_nIC' num2str(nIC) '_thresh' num2str(cluster_threshold) ...
    '_coor' num2str(threshold) '_mask3D.mat'], 'mask3D');
mask3D_final = mask3D; mask3D_final(:,:,badPCs) = [];
save([analysis_out date '_img' mouse '\getTC\' date '_img' mouse '_nPCA' num2str(nPCA) '_mu' num2str(mu) '_nIC' num2str(nIC) '_thresh' num2str(cluster_threshold) '_coor' num2str(threshold) '_mask3D.mat'], 'mask3D_final','-append');

save([analysis_out date '_img' mouse '\' date '_img' mouse '_deconvolution_thresh', num2str(thresholdDeco), '_TCave_cl.mat'],...
    'TCave_cl','TCave','thresholdDeco','badPCs','nobase_cell','low_Fcell','nobase_thres','baseline_thres','all_TCave_cl','goodcells');
save([analysis_out date '_img' mouse '\' date '_img' mouse '_spk_deconvolve_threshold' num2str(thresholdDeco) '.mat' ],'FRstay_cell_cl', ...
    'aveFR_cl','badPCs','spk_logic_cl','spk_logic','aveFR','spk','kernel','spk_peak','spk_inx','spk_cl','kernel_cl','spk_peak_cl',...
    'spk_inx_cl','badPCs','nobase_cell','low_Fcell','nobase_thres',...
    'baseline_thres','-append');

%% make a new mask figure
% load mask3D final
filename2 = dir([analysis_out,date '_img' mouse '\getTC\' '*' 'thresh' num2str(cluster_threshold) '_coor' num2str(threshold) '_mask3D.mat']);
mask3D = load([analysis_out date '_img' mouse '\getTC\' filename2.name]);
mask3D = mask3D.mask3D;
tempFig=setFigure; hold on;
%needs to modify the line below due to different file names
imshow([analysis_out,date '_img' mouse '\getTC\' 'AVG_' date '_img' mouse '_rgstr_tiff_every50_ref' num2str(ref) '.tif']); hold on;
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

tempFig=setFigure; hold on;
%needs to modify the line below due to different file names
imshow([analysis_out,date '_img' mouse '\getTC\' 'AVG_' date '_img' mouse '_rgstr_tiff_every50_ref' num2str(ref) '.tif']); hold on;
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


%close all;
% 