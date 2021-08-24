clear; close all;
% ref=1;
analysis_out = 'A:\home\carlo\analysis\2P\';
bdata_source = 'A:\home\carlo\rawData\behavioral\';
write_tiff = true; do_jake_pockel = true;
ttl =   1;
noreg = 0;
imgreglaseron = 0;
% sessions = '2105170_img1705';% for behavior data
% mouse = '1705';
% date = '2105170';
% run = '000';

%RCExptListCarlo_Inter;
%RCExptListCarlo_Block;
expType=input('monomodal(1) or dimodal(2) session: ');
seshType=input('interleaved (1) or blocked (2) session: ');
if seshType==1 && expType == 1
RCExptListCarlo_Inter;
elseif seshType==2 && expType == 1
RCExptListCarlo_Block;
elseif seshType==1 && expType == 2
RCExptListCarlo_VA_Inter;
end

pairings = loadRCList_Carlo(expt,seshType);
iexp = pairings{1,1}; isesh = pairings{3,1};
mouse = strtrim(expt{1,iexp}.mouse{1,isesh});
date = expt{1,iexp}.date{1,isesh};
run = expt{1,iexp}.run{1,isesh};
sessions = [date '_img' mouse];

%cut will be filled with frames with poor illumination, noise, etc.
cut = [];
%time = expt(id).time(iexp,:);
fprintf([date ' ' mouse ' ' run '\n'])
img_fn = [date '_img' mouse '\getTC_' run '\'];
if ~exist(fullfile(analysis_out,img_fn))
    mkdir(fullfile(analysis_out, img_fn))
end

%% load and remove laser off period
cd(fullfile('A:\home\carlo\rawData\2P\', [date '_img' mouse]));
load([mouse '_' date '_' run '.mat']);

%identify behavior file
if strlength(sessions) < 16
    bfile = dir([bdata_source 'data-i' mouse '-' date '*' ]);
else
    bfile = dir([bdata_source 'data-i' mouse '-' date sessions(end-4:end) '.mat']);
end
%load behavior file
behav_dest = [bdata_source bfile.name];
assert(length(bfile)) = 1;
b_data = load(behav_dest);
mworks = b_data.input;

clear ttl_log
nf= mworks.counterValues{end}(end);
if mod(nf,10000)==2
    nf = nf-2;
elseif mod(nf,10000)==1
    nf = nf-1;
end

if ttl == 1
    load([mouse '_' date '_' run '_realtime.mat'])
    disp(['Loading ' num2str(nf) ' frames'])
    data = sbxread([mouse '_' date '_' run],0,nf);
    data = squeeze(data);
    if do_jake_pockel && length(unique(ttl_log)) == 1 || sum(ttl_log) < 0.05*(length(ttl_log))  %if ttl_log is faulty then identify pockel transitions via changes in F
        pockel_tc = find_pockel_tc_SJ(data, 0);
        if pockel_tc(1)==1
            pockel_tc(1,1:2) = 0; %delete the first two frames because the image is partial when experiment first starts. don't need to do this if laser if off from the beginning
            data_cut = data(:,:,(pockel_tc==1)); % get rid of laser off period
        else
            data_cut = data(:,:,(pockel_tc==1));
        end
    else
        ttl_trans = find(abs(diff(ttl_log)));
        n = length(ttl_trans);
        for i = 1:n
            ttl_log(ttl_trans(i)-4:ttl_trans(i)+4,:) = 0;
        end
        ttl_ind = find(ttl_log(1:nf,:) ==0);
        data_cut = data;
        data_cut(:,:,ttl_ind) = [];
    end
else
    if nf > 60000
        nf = 60000;
    end
    fprintf('Loading data \n')
    data = sbxread(['img' mouse '_000_' run],0,nf);
    data = squeeze(data);
end

laseron = find(pockel_tc==1);
disp(['Data is now ' num2str(size(data_cut,3)) ' frames long after removing laser-off periods' ])
if write_tiff && ~exist([analysis_out, img_fn, date '_img' mouse '_0k.tif'])  
for i = 1:floor(length(data_cut)/5000)
    saveastiff(data_cut(:,:,0+((i-1)*5000)+1:1000+((i-1)*5000)),fullfile(analysis_out,img_fn, [date '_img' mouse '_' num2str((i-1)*5000) 'k.tif']));
end %automated incrimental tiff save for cut data (will progress from 1:1000-100001:101000 with no deviation [if necessary]
end

%   saveastiff(data_cut(:,:,1:1000),fullfile(analysis_out,img_fn, [date '_img' mouse '_1k.tif']));


%% register
if exist(fullfile(analysis_out, img_fn, [date '_img' mouse '_reg.mat']))
    load(fullfile(analysis_out, img_fn, [date '_img' mouse '_reg.mat']))
else
    figure;
    nplot = floor(size(data_cut,3)./2000);
    [n n2] = subplotn(nplot);
    for i = 1:nplot
        subplot(n,n2,i)
        imagesc(mean(data_cut(:,:,1+(i-1).*2000:500+(i-1).*2000),3))
        title(num2str(i));
    end
    prompt = 'Choose a registration frame \n';
    ref = input(prompt);
    img_ref = mean(data_cut(:,:,1+(ref-1).*2000:300+(ref-1).*2000),3);
    [npw, nph, nt] = size(data_cut);
    save(fullfile([analysis_out, img_fn, date '_img' mouse  '_reg.mat']), 'ref', 'img_ref', 'npw', 'nph', 'nt');
end

if size(data_cut,3) >= 60000
    data_cut = data_cut(:,:,1:2:end);
end
fprintf(['Registering ' num2str(size(data_cut,3)) ' frames'])
[rgs_out,img_rgs]= stackRegister(data_cut, img_ref);

gap = 50;
idx = (1:gap:nt);
writetiff(img_rgs(:,:,idx),[analysis_out,img_fn date '_img' mouse '_' run...
    '_rgstr_tiff_' 'every',num2str(gap) '_ref' num2str(ref)]);
prompt = 'look at the tiff, any shifts? Yes there is shift/No did not see any: \n';
response = input(prompt,'s');
switch response
    case 'No'; disp('NICE!! PLEASE SAVE average image as jpeg and go to PCA');
        %dont actually save anything??
        clear data;
    case 'Yes'; disp('go back and re-register :(');
    otherwise; disp('Please enter "Y" for yes and "N" for no.')
end

%% this image should tell you if your cells get bleached over time during imaging, if the overall F goes up gradually, there's nothing you can do but you know there's bleach. And that might influence deconvolution
figure; plot(squeeze(mean(mean(img_rgs,2),1)));
savefig([analysis_out,img_fn date '_img' mouse '_' run '_raw_F_wholeview']);

%% do this as needed: remove part of the data because of z shift
% cutting the data shouldn't influence the registration as long as you
% didn't take the ref frame from the part you need to cut. But it will
% influence PCA ICA, so do it before then.
%!!!!!! if you do this, you also need to rewrite cTargetOn. Otherwise you'll
%be referencing some trials with no neural data in the later analysis.
%cut = 25000:size(data_cut,3);

% cut1 = 1;
% cut2 = 2432;
% cut = cut1:1:cut2;
% img_rgs(:,:,cut) = [];
% laseron(cut1:cut2) = [];

%% generate cTargetOn_cutted (index of cTargetOn in the cutted neural data)
cTargetOn = double(cell2mat(mworks.cTargetOn));

%cTargetOn = cTargetOn(1,2:end);
cTargetOn(1) = NaN; % get rid of the first trial because juice = 0

cTargetOn_cutted = [];
for i = 1:length(cTargetOn)
    if ~isnan(cTargetOn(i))
        cTargetOn_cutted = [cTargetOn_cutted find(laseron==cTargetOn(i))];  % find the index of cTargetOn in cutted data
    end
end

%deconvolution needs 120 frames before cTargetOn, so need to remove trials
%that begin with frame numbers <120
for i = 1:length(cTargetOn_cutted)
    if cTargetOn_cutted(1,i) < 120
        cTargetOn_cutted(1,i) = NaN;
    end
end


if length(cTargetOn_cutted) <= length(cTargetOn) % if there're cue presentations when there's no nice imaging data (this happens if you cut the some part of the imaging data due to z shift or other reasons)
    %cTargetOn_wneural = [];
    prewin_frames = round(1500./mworks.frameRateHz);
    for i = 1:length(cTargetOn)
        if ~isnan(cTargetOn(i))
            if find(laseron==cTargetOn(i))  % cue was during laser on
                ind = find(laseron==cTargetOn(i));
            else
                cTargetOn(i) = nan;
%                 if sum(find(diff(laseron(ind-prewin_frames: ind))>1))==0 %and the time before cue is long enough
%                     %cTargetOn_wneural = [cTargetOn_wneural,cTargetOn(i)];
%                     continue;
%                 else
%                     cTargetOn(i) = nan;
%                 end
            end
        end
    end
    % overwrite cTargetOn with only the ones that has neural data
    mworks.cTargetOn = cTargetOn;
    input = mworks;
    save(fullfile(behav_dest), 'input'); % overwrite cTargetOn
    clearvars input
end

%% PCA
nPCA = 196; %100 for old datasets, 500 for newer 
imgsize = size (img_rgs,3);
[mixedsig_PCA, mixedfilters_PCA, CovEvals_PCA, ~, ~, ~] = CellsortPCA_2P_Jin(img_rgs,[1 imgsize], nPCA,[], []);
figure; plot(CovEvals_PCA); % looks like can take 2-60th PCA or sth
savefig([analysis_out,img_fn date '_img' mouse '_' run '_PCA', num2str(nPCA) '_CoEvals']);

figure;
for i = 1:196%(floor(sqrt(nPCA)))^2
    subplot(14,14,i); imagesc(mixedfilters_PCA(:,:,i));
    set(gca,'xticklabel',[],'yticklabel',[]);
    text(.8,.1,num2str(i),'fontsize',12,'color','w','fontweight','bold','unit','norm');
    colormap gray
end
% for i = 197:392%(floor(sqrt(nPCA)))^2
%     subplot(14,14,i-196); imagesc(mixedfilters_PCA(:,:,i));
%     set(gca,'xticklabel',[],'yticklabel',[]);
%     text(.8,.1,num2str(i),'fontsize',12,'color','w','fontweight','bold','unit','norm');
%     colormap gray
% end
savefig([analysis_out,img_fn date '_img' mouse '_' run '_PCA', num2str(nPCA)]);
save([analysis_out,img_fn date '_img' mouse '_' run '_PCA_variables_', num2str(nPCA),'.mat'], ...
    'mixedsig_PCA', 'mixedfilters_PCA', 'CovEvals_PCA', 'nPCA');


%% ICA: seperates independent spatial and temporal components
%PCuse =       1:125;
PCuse =       1:size(mixedfilters_PCA,3);
mu =          0.3; % weight of temporal info in spatio-teporal ICA
nIC =         input('Look at PCA output, how many IC should be included in ICA: '); % cannot be bigger than nPCA. If CoEvals doesn't change in later ICs, it will not converge!
ica_A_guess = []; %If this is empty than matlab will randomdize it and you can get different results, can see the random number generator in CellsortICA2P
termtol =      1e-6;
maxrounds =   2000;
npw =         264;
nph =         796;
[ica_sig, mixedfilters_ICA, ica_A, numiter] = CellsortICA_2P(mixedsig_PCA,...
    mixedfilters_PCA,...
    CovEvals_PCA, ...
    PCuse,...
    mu,...
    nIC,...
    [],...
    termtol,...
    maxrounds);

icasig = permute(mixedfilters_ICA,[2,3,1]);% make it the same dimension as the raw data (pixel by pixel by time)
ICA_variables.mu = mu; ICA_variables.nIC = nIC;  ICA_variables.termtol = termtol; ICA_variables.naxrounds = maxrounds;
ICA_variables.npw = npw;  ICA_variables.nph = nph;
save([analysis_out,img_fn date '_img' mouse '_' run '_nIC' num2str(nIC) '_variables_for_nPCA_', ...
    num2str(nPCA), '.mat'], 'ica_sig', 'icasig', 'ICA_variables');
writetiff(icasig, [analysis_out,img_fn date '_img' mouse '_' run '_icasig_for_nPCA', num2str(nPCA) '_nIC' num2str(nIC)]);
figure; imagesc(sum(icasig,3));
savefig([analysis_out,img_fn date '_img' mouse '_' run '_icasig_sum_', num2str(nPCA) '_nIC' num2str(nIC)]);


%% smooth and threshold the mask
%stack filter - acts as a low pass spatial filter, reduces noise. It will apply the Gaussian filter, and filter out low spatial frequency noise such as small blips and bloops in the ICs which do not belong to dendrites.
icasig_filt = stackFilter(icasig);

%set threshold a threshold for which pixels to include in a given dendrite's mask.
nIC = size(icasig_filt, 3);
cluster_threshold = 97; % this is using the top 3 percent of the fluorescence values, so brightest 3% is yes (1), and the rest is no (0)
%tried lower values for the threshold and turns out to have some wierd
%masks
mask_cell = zeros(size(icasig_filt));

for ic = 1:nIC
    %initialize sm_logical
    sm_logical = zeros(npw,nph);
    %convert to a binary mask (0 and 1)
    icasig_filt(:,:,ic) = imclearborder(icasig_filt(:,:,ic));
    sm_logical((icasig_filt(:,:,ic)> mean([max(prctile(icasig_filt(:,:,ic),cluster_threshold,1)) max(prctile(icasig_filt(:,:,ic),cluster_threshold,2))])))=1;
    sm_logical((icasig_filt(:,:,ic)<=mean([max(prctile(icasig_filt(:,:,ic),cluster_threshold,1)) max(prctile(icasig_filt(:,:,ic),cluster_threshold,2))])))=0;
    sm_logical = logical(sm_logical);
    %bwlabel identifies unconnected objects within a single mask. So 0=background 1=Object#1 2=Object#2
    if sum(sum(sm_logical,2),1) <51
        sm_logical = 0;
    end
    mask_cell(:,:,ic) = bwlabel(sm_logical);
end
save([analysis_out,img_fn date '_img' mouse '_' run '_nPCA', ...
    num2str(nPCA), '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', ...
    num2str(cluster_threshold), '_mask.mat'], 'mask_cell');

%visualize the masks
figure('rend', 'painters', 'pos', [50 150 (796*1.5) (264*1.5)]); imagesc(sum(mask_cell,3));...
    title([date '_img' mouse '_' run '_', ' nPCA ', num2str(nPCA), ' mu ', num2str(mu), ' nIC ', num2str(nIC)]);
savefig([analysis_out,img_fn date '_img' mouse '_' run '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), '_mask_cell_sum.fig']);

figure;
for i = 1:100
    subplot(10,10,i); imagesc(mask_cell(:,:,i));
    colormap gray
end
savefig([analysis_out,img_fn date '_img' mouse '_' run, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), '_mask_cell_indi.fig']);

%% separate overlapping masks and combine the ones that have a high correlation
%split individual masks, remove small masks, deal with overlapping
mask_final = processMask(mask_cell);
mask_raw = reshape(mask_final, npw, nph);
figure; imagesc(mask_raw); truesize;
savefig([analysis_out,img_fn date '_img' mouse '_' run, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), '_maskS_process.fig']);

% combine highly correlated ICs to one
threshold = 0.99; %got the same thing when threshold = 0.8 and 0.9, tried different values, doesn't seem to change things
[ ~, mask3D, ~] = finalMask_Jin(img_rgs, mask_final, threshold);
%figure; imshow([image_analysis_dest 'AVG_' sessions '_' order '_rgstr_tiff_0_' num2str(nframes) '_50_jpeg.jpg']); hold on;
figure; imshow([analysis_out,img_fn 'AVG_' date '_img' mouse '_' run '_rgstr_tiff_every' num2str(gap) '_ref' num2str(ref) '_jpeg.jpg']); hold on;
for i  = 1:size(mask3D,3)
    bound = cell2mat(bwboundaries(mask3D(:,:,i)));
    randcolor = rand(1,4);
    plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;
    text(mean(bound(:,2)),mean(bound(:,1)), ...
        num2str(i), 'color', 'y', 'FontSize', 8);
end
%figure; imagesc(sum(mask3D,3)); truesize; % got the same thing as the figure above
savefig([analysis_out,img_fn date '_img' mouse '_' run, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), ...
    '_coor' num2str(threshold),'_mask_wdendrites.fig']);
save([analysis_out,img_fn date '_img' mouse '_' run, '_nPCA', ...
    num2str(nPCA), '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', ...
    num2str(cluster_threshold), '_coor' num2str(threshold) '_mask3D.mat'], 'mask*');


%% look at the masks, mananully delete the ones that don't look like dendrites
startRemove = input('remove any ROIs? (1-yes | 2-no): '); allROItoRemove = [];
if startRemove == 1
doneRemove = 2;
while doneRemove == 2
removeCell = input('remove which ROI: ');
if removeCell == 0
    doneRemove=1;
elseif removeCell ~= 0
allROItoRemove = [allROItoRemove removeCell];
doneRemove = input('finished removing ROI? (1-yes | 2-no): ');
end

end
end

mask3D(:,:,[allROItoRemove]) = [];
%figure; imshow([image_analysis_dest 'AVG_' sessions '_' order '_rgstr_tiff_0_' num2str(nframes) '_50_jpeg.jpg']); hold on;
figure; imshow([analysis_out,img_fn 'AVG_' date '_img' mouse '_' run '_rgstr_tiff_every' num2str(gap) '_ref' num2str(ref) '_jpeg.jpg']); hold on;
for i  = 1:size(mask3D,3)
    bound = cell2mat(bwboundaries(mask3D(:,:,i)));
    randcolor = rand(1,4);
    plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;
    text(mean(bound(:,2)),mean(bound(:,1)), ...
        num2str(i), 'color', 'y', 'FontSize', 8);
end
hold on;
title([date '_img' mouse '_masks_final']);
savefig([analysis_out,img_fn date '_img' mouse '_' run, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), ...
    '_coor' num2str(threshold),'_mask_wdendrites_final.fig']);
save([analysis_out,img_fn date '_img' mouse '_' run, '_nPCA', ...
    num2str(nPCA), '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', ...
    num2str(cluster_threshold), '_coor' num2str(threshold) '_mask3D_final.mat'], 'mask3D');


%% get TCs
nmask = size(mask3D,3);
FrameRate = 30;
raw_tc_avg = getTC(img_rgs, mask3D, nmask);


%%
figure;

cutoffF = .035;
%high-pass filter the time course before replacing laser off periods with NaN
[b,a]=butter(3,(cutoffF*2/FrameRate),'high');
tc_avg = filtfilt(b,a,raw_tc_avg); %zero-phase filtering using parameters from butterworth high-pass filter to output filtered TC
%to test the alignment by removing the session mean from the unfiltered TC
baseOrigTC = nanmean(nanmean(raw_tc_avg(1:1000,:),2),1); %find average of entire TC
origTC = nanmean(raw_tc_avg,2);  %call TC collapsed across neurons
normTC = tc_avg+baseOrigTC; %remove mean from collapsed TC
%plot raw and filtered TC
plot(origTC,'r');
hold on;
plot(nanmean(normTC,2),'k'); %collapse filtered TC across neurons to plot
title(['cutoff frequency at ' num2str(cutoffF) ' Hz']);
legend('raw TC','filtered TC','Location','northwest');
savefig('rawVSfilteredTC_0075.fig');


rgstr_sum = sum(img_rgs,3);
plotTC_Jin(tc_avg, mask3D, rgstr_sum, 1:size(tc_avg,2), FrameRate);
% if plotting all the neurons make all lines look flat, plot some of them
% plotTC_Jin(tc_avg,~,~,1:40,~);
savefig([analysis_out,img_fn date '_img' mouse '_' run, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), '_TC.fig']);
mask_flat = plotMask_Jin(mask3D,1); %second input is either 0 or 1, if 1, makes a figure.
savefig([analysis_out,img_fn date '_img' mouse '_' run, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), '_mask_plotMask.fig']);

data_corr = corrcoef(tc_avg);
figure; fig = imagesc(data_corr);
savefig([analysis_out,img_fn date '_img' mouse '_' run, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), '_corrcoef.fig']);

save([analysis_out,img_fn date '_img' mouse '_' run, '_nPCA', ...
    num2str(nPCA), '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', ...
    num2str(cluster_threshold), '_TCave.mat'], 'tc_avg','raw_tc_avg','pockel_tc','cutoffF');

%% make movie and revert tc_avg to inital counters by adding Nans

tc_avg_all = NaN(length(pockel_tc),nmask);
tc_avg_all(laseron,:) = tc_avg;

fprintf('Making movie \n')
sz = size(img_rgs);
nTrials = length(cTargetOn_cutted)-1;
prewin_frames = round(1500./mworks.frameRateHz); %
postwin_frames = round(3000./mworks.frameRateHz);
rgs_align = nan(sz(1),sz(2),prewin_frames+postwin_frames,nTrials);
for itrial = 1:nTrials
    if cTargetOn_cutted(itrial)+postwin_frames<size(img_rgs,3) 
        rgs_align(:,:,:,itrial) = img_rgs(:,:,cTargetOn_cutted(itrial)-prewin_frames:cTargetOn_cutted(itrial)+postwin_frames-1);
    end
end
rgs_align_avg = mean(rgs_align,4); % your rgs_align shouldn't have nan in it
writetiff(rgs_align_avg, fullfile([analysis_out,img_fn date '_img' mouse '_' run,'_cueAlign.tif']));

%% BUTTERWORTH HIGH-PASS FILTER USED ON TC_AVG

save([analysis_out,img_fn date '_img' mouse '_' run, '_nPCA', ...
    num2str(nPCA), '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', ...
    num2str(cluster_threshold), '_TCave.mat'],'cut', 'tc_avg', 'tc_avg_all','cTargetOn_cutted','laseron','-append');

save([analysis_out,img_fn date '_img' mouse '_' run, 'saveOutputs.mat'], 'nPCA','mu','nIC','cluster_threshold','threshold','ref');