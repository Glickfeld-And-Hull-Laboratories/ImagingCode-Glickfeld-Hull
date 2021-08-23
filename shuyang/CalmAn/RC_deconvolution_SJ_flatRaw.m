clear;
analysis_out = 'Z:\2P_analysis\';
bdata_source = 'Z:\Data\behavior\RC\';

ttl =   1;
noreg = 0;
imgreglaseron = 0;

mouse = '1083';
date = '210112';
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


%%
% CellDataFlat01 
%   flat the cell data by fitting with polynomial curve fit
data_length = size(tc_avg,1);
fit_order = 4;
fit_tc_avg = zeros(size(tc_avg,1),size(tc_avg,2));
for c = 1:size(tc_avg,2)
    pfit = polyfit((1:1:data_length),(tc_avg(:,c))',fit_order);
    fit_tc_avg(:,c) = polyval(pfit,(1:1:data_length)); 
end

% flat data
tc_avg_flat = tc_avg - fit_tc_avg;

% ------- fit flat data
fit_tc_avg_flat = zeros(size(tc_avg,1),size(tc_avg,2));
for c = 1:size(tc_avg,2)
    pfit2 = polyfit((1:1:data_length),(tc_avg_flat(:,c))',0);
    fit_tc_avg_flat(:,c) = polyval(pfit2,(1:1:data_length));
end

for c = 1:size(tc_avg,2)
    figure;
    plot((1:1:data_length),tc_avg(:,c),'b-');hold on;
    plot((1:1:data_length),fit_tc_avg(:,c),'r-','linewidth',3);
    plot((1:1:data_length),tc_avg_flat(:,c),'k-');
    plot((1:1:data_length),fit_tc_avg_flat(:,c),'r--','linewidth',3);
    hold off;
end

% % ------- parameters 
% test_cell_number = 10; % pick up a cell for testing 
% fit_order = 4; 
% 
% % ------- generate test data 
% test_cell_data = tc_avg(:,test_cell_number)';
% test_cell_data_length = length(test_cell_data);
% test_cell_data_frame = 1:1:test_cell_data_length;
% 
% % ------- fit data,polynomial curve fitting 
% p_fit = polyfit(test_cell_data_frame,test_cell_data,fit_order);
% test_cell_data_fit = polyval(p_fit, test_cell_data_frame);

% % ------- flat data 
% test_cell_data_flat = test_cell_data-test_cell_data_fit;

% % ------- fit flat data 
% p2_fit = polyfit(test_cell_data_frame,test_cell_data_flat,0);
% test_cell_data_flat_fit = polyval(p2_fit, test_cell_data_frame);

% ------- plot 
% f1 = figure;
% set(gcf,'unit','centimeters','position',[2 2 25 15]);
% hold on 
% plot(test_cell_data_frame,test_cell_data,'b-')
% plot(test_cell_data_frame,test_cell_data_fit,'r-','linewidth',3)
% plot(test_cell_data_frame,test_cell_data_flat,'k-')
% plot(test_cell_data_frame,test_cell_data_flat_fit,'r--','linewidth',3)
% hold off 
% legend('raw data','fit data','flat data','flat data fit','location','bestoutside')
% legend boxoff 
% title(['flattened for cell #' num2str(test_cell_number) ' by polynomial fit (power of ' num2str(fit_order) ')'])
% set(gca,'FontSize',18);
% set(gca,'FontName','Times New Roman');
% 
% % ------- figure save 
% fig_name = ['CellDataFlat_Cell' num2str(test_cell_number) '_PloyFit_PowerOf' num2str(fit_order)];
% print(f1,fig_name,'-r600','-djpeg')
% savefig(f1,fig_name)



%% Deconvolution
% find the baseline states: 4s-120frames before cTargetOn
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

% deconvolution on flattened data
TCave = tc_avg_flat;
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
save([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_spk_deconvolve_flattened_threshold' num2str(threshold) '.mat' ],'threshold',...
    'FRstay_cell', 'aveFR','num_spks_cell','options','spk_logic','spk','kernel','spk_peak','spk_inx');

%% plots
% write a GUI with subplots: TCave with red dots and kernel
hist_FR = figure;
histogram(FRstay_cell,'BinWidth',0.1);
title([date ' img' mouse 'deconvolution-flattened-firing rate-stationary-threshold' num2str(threshold)]);
saveas(hist_FR,[analysis_out decon_out date '_img' mouse '_histFR_deconvolution_flattened_threshold' num2str(threshold) '.fig']);

[fig_deconvolve] = GUI_rawTrace_nDeconvolve(TCave,kernel,spk_logic,[date '_img' mouse],threshold);
savefig([analysis_out decon_out date '_img' mouse '_GUI_TCave_deconvolution_flattened_threshold' num2str(threshold) '.fig']);

% spike logic
% plot spike logic
figure;
colormap gray;
imagesc(imcomplement(spk_logic'));
xlabel('Frame #');
ylabel('Cell #');
%xlim([25000 28992]);
title ([date ' img' mouse ' deconvolution-flattened']);
savefig([analysis_out decon_out date '_img' mouse '_scatter_spikeLogic_deconvolve_flattened_threshold' num2str(threshold) '.fig']);

% plot spike rates for each cell during stationary
figure;
scatter(1:length(FRstay_cell), FRstay_cell,70, 'ko');
ylabel('Spike Rate (Hz)');
xlabel('Cell #'); axis square;
%ylim([0 1.4]);
title ([date ' img' mouse ' deconvolution-flattened']);
savefig([analysis_out decon_out date '_img' mouse '_scatter_spikeRate_deconvolve_flattened_threshold' num2str(threshold) '.fig']);


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
nobase_thres = 10000; % determine this number after looking at the plot. sort out the outliers, I typically use a number in between 2000-4000. 
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
imshow([analysis_out,img_fn 'AVG_' date '_img' mouse '_' run '_rgstr_tiff_every50_ref4_jpeg.jpg']); hold on;
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


