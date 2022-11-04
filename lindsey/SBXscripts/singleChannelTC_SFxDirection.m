%% Load, register, segment and neuropil correct 2P data
close all
clear all global
clc

%Path names
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
lg_fn = fullfile(fn_base, 'home\lindsey');
data_fn = fullfile(lg_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data');
fnout = fullfile(lg_fn, 'Analysis\2P');

%Specific experiment information
date = '220905';
ImgFolder = '001';
time = '1125';
mouse = 'i1377';

ds = 'SFxDir_ExptList';
eval(ds)

iexp = 11;

mouse = expt(iexp).mouse;
date = expt(iexp).date;
time = cell2mat(expt(iexp).sfTime);
ImgFolder = cell2mat(expt(iexp).sfFolder);

frame_rate = 15;
run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];

%Load 2P data
%Load mworks data- this has information about experiment (e.g. the visual stimuli presented and synchronization with the microscope)
fName = fullfile(mworks_fn, ['data-' mouse '-' date '-' time '.mat']);
load(fName);
%Load 2P metadata- this has information about the information about the imaging session (e.g. frame count, zoom)
CD = fullfile(data_fn, mouse, date, ImgFolder);
cd(CD);
imgMatFile = [ImgFolder '_000_000.mat'];
load(imgMatFile);
%Load 2P images
totframes = input.counterValues{end}(end); %this is from the mworks structure- finds the last value clocked for frame count
fprintf(['Reading ' num2str(totframes) ' frames \r\n'])
data = sbxread([ImgFolder '_000_000'],0,totframes); %loads the .sbx files with imaging data (path, nframes to skip, nframes to load)
%Data is nPMT x nYpix x nXpix x nframes. 
fprintf(['Data is ' num2str(size(data))])
%When imaging a single channel, nPMT = 1, so squeeze:
data = squeeze(data);

%% Register 2P data
%Goal here is to remove X-Y movement artifacts
%1. Find a stable target
%    a. Plot average of 500 frames throughout stack
if ~exist(fullfile(fnout, datemouse, datemouserun, [datemouserun '_reg_shifts.mat']))
nframes = 500; %nframes to average for target
nep = 9;
nskip = floor(totframes./nep); %nframes to skip for each average
[n, n2] = subplotn(nep); 
figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:nep 
    subplot(n,n2,i); 
    imagesc(mean(data(:,:,1+((i-1)*nskip):nframes+((i-1)*nskip)),3)); 
    title([num2str(1+((i-1)*nskip)) '-' num2str(nframes+((i-1)*nskip))]); 
end

%    b. GUI to select target image- choose one that is sharp and close to center of stack
f=gcf;
w = waitforbuttonpress; %click on subplot
if w == 0
    axesClicked = gca;
    allAxes = flipud(findobj(f.Children,'Type','axes'));
    numClicked = find(axesClicked==allAxes);
    close all
end
fprintf(['Selected subplot ' num2str(numClicked)])
%    c. Create target image
data_avg = mean(data(:,:,1+((numClicked-1)*nskip):nframes+((numClicked-1)*nskip)),3); %average 500 frames to make target
%2. stackRegister minimizes the difference of each frame from the target
[out, data_reg] = stackRegister(data,data_avg);
%New average image after registration
data_reg_avg = mean(data_reg,3);
%Save registration shifts and target, and mworks data
mkdir(fullfile(fnout, datemouse, datemouserun))
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_reg_shifts.mat']), 'data_reg_avg', 'out', 'data_avg')
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_input.mat']), 'input')
%Test registration
%    a. Make sure first and last images are sharp and similar. Focus on vasculature and nuclei- not cells since F changes
ind = [1 nep];
for i = 1:length(ind) 
    subplot(2,1,i); 
    ix = ind(i);
    imagesc(mean(data_reg(:,:,1+((ix-1)*nskip):nframes+((ix-1)*nskip)),3)); 
    title([num2str(1+((ix-1)*nskip)) '-' num2str(nframes+((ix-1)*nskip))]); 
end
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_FOV_first&last.pdf']), '-dpdf')
%    b. Average of all frames should be sharp
imagesq(data_reg_avg); 
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_FOV_avg.pdf']), '-dpdf')
else
    load(fullfile(fnout, datemouse, datemouserun, [datemouserun '_reg_shifts.mat']))
    [out, data_reg] = stackRegister_MA(data,[],[],out);
end

clear data

%% Segment 2P data
%Goal here is to create a cell mask to extract fluorescence timecourses
%1. Find activated cells
%    a. Reshape data stack to segregate trials 
%    First need to know how many frames per trial and how many trials
nOn = input.nScansOn;
nOff = input.nScansOff;
ntrials = size(input.tGratingDirectionDeg,2); %this is a cell array with one value per trial, so length = ntrials
sz = size(data_reg);
%        Each trial has nOff frames followed by nOn frames, so can reshape stack so nYpix x nXpix x nFrames/Trial (nOn+nOff) x nTrials
data_tr = reshape(data_reg,[sz(1), sz(2), nOn+nOff, ntrials]);
data_tr_inv = mean(mean(data_tr(:,:,nOff/2:nOff,:),3),4)-mean(mean(data_tr(:,:,nOff+1:nOn+nOff,:),3),4);

fprintf(['Size of data_tr is ' num2str(size(data_tr)) '\n'])
%    b. Find baseline F from last half of off period- avoids decay of previous on trial
fprintf('Making data_f \n')
data_f = mean(data_tr(:,:,nOff/2:nOff,:),3); 
%    c. Find dF/F for each trial
fprintf('Making data_df \n')
data_df = bsxfun(@minus, double(data_tr), data_f);
fprintf('Making data_dfof \n')
data_dfof = bsxfun(@rdivide,data_df, data_f); 
fprintf('Clearing data \n')
clear data_f data_df data_tr
%    d. Find average dF/F for each stimulus condition (this is for an experiment with changing grating direction)
Dir = celleqel2mat_padded(input.tGratingDirectionDeg); %transforms cell array into matrix (1 x ntrials)
if find(Dir>+360)
    Dir(find(Dir>=360)) = Dir(find(Dir>=360))-360;
end
Dirs = unique(Dir);
nDirs = length(Dirs);

SF_mat = celleqel2mat_padded(input.tGratingSpatialFreqCPD); %transforms cell array into matrix (1 x ntrials)
SFs = unique(SF_mat);
nSF = length(SFs);
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_stimData.mat']), 'Dir', 'Dirs', 'nDirs', 'SF_mat', 'SFs', 'nSF', ...
    'nOn', 'nOff', 'ntrials')
nStim = nDirs+nSF;
data_dfof_avg = zeros(sz(1),sz(2),nStim); %create empty matrix with FOV for each direction: nYpix x nXPix x nDir

fprintf('Averaging Dirs: \n')
figure;
[n, n2] = subplotn(nDirs); %function to optimize subplot number/dimensions
for idir = 1:nDirs
    fprintf(num2str(idir))
    ind = find(Dir == Dirs(idir)); %find all trials with each direction
    data_dfof_avg(:,:,idir) = mean(mean(data_dfof(:,:,nOff+1:nOn+nOff,ind),3),4); %average all On frames and all trials
    subplot(n,n2,idir)
    imagesc(data_dfof_avg(:,:,idir))
    title(['Dir- ' num2str(Dirs(idir))])
end
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_allDirs.pdf']), '-dpdf')

figure;
[n, n2] = subplotn(nSF); %function to optimize subplot number/dimensions
fprintf('\n Averaging SFs: \n')
for iSF = 1:nSF
    fprintf(num2str(iSF))
    ind = find(SF_mat == SFs(iSF)); %find all trials with each SF
    data_dfof_avg(:,:,iSF+nDirs) = mean(mean(data_dfof(:,:,nOff+1:nOn+nOff,ind),3),4); %average all On frames and all trials
    subplot(n,n2,iSF)
    imagesc(data_dfof_avg(:,:,iSF+nDirs))
    title(['SF- ' num2str(SFs(iSF))])
end
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_allSF.pdf']), '-dpdf')

data_dfof_avg = cat(3,data_reg_avg,data_dfof_avg);
data_dfof_avg = cat(3,data_dfof_avg,data_tr_inv);
clear data_dfof
%        Filtering data helps make cells more visible for selection
myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_avg_all = imfilter(data_dfof_avg,myfilter);
data_dfof_max = max(data_dfof_avg_all,[],3); %finds all active cells by taking max projection

save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg_all', 'nStim')
%    2. Create cell masks from active cells
%        Concatenate images of max and all directions
data_dfof = cat(3,data_dfof_max, data_dfof_avg_all);
sz = size(data_dfof);
%        Set up empty matrices for segmenting cells
mask_data = data_dfof;
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
%        a. Click on bright cells in mask to create mask- blue shading should fill bright region of single cell evenly.  
%            If cell not filled well, press 'z' or shift+click to remove. Try again by 1) clicking on different spot- brighter spot will make fill smaller- if there are two cells also helps to fill brighter cell first; 2) change threshold- higher thresh will make fill smaller; 3) wait for different image where cell might be more distinct from background
%            Press 'return' when done with each image. 
%            Next images will have black spots where cells were previously selected- don't let new masks encroach on old cells or will fuse. imCellBuffer helps with this.
for iStim = 1:size(data_dfof,3)    
    mask_data_temp = mask_data(:,:,iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0; %blacks out old cells
    bwout = imCellEditInteractiveLG(mask_data_temp); %selection GUI
    mask_all = mask_all+bwout; %adds new cells to old cells
    mask_exp = imCellBuffer(mask_all,3)+mask_all; %creates buffer around cells to avoid fusing
    close all
end
%        b. Number independent regions for mask
mask_cell = bwlabel(mask_all); %turns logical into numbered cells
figure;
imagesc(mask_cell)
%        c. Create neuropil masks (these are regions around each cell without overlap from neighboring cells)
nMaskPix = 5; %thickness of neuropil ring in pixels
nBuffPix = 3; %thickness of buffer between cell and ring
mask_np = imCellNeuropil(mask_cell,nBuffPix,nMaskPix);

save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_np')
clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_2 data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 

%% Neuropil subtraction
%Goal is to remove contamination from out-of-focus fluorescence
%Extract cell timecourses
data_tc = stackGetTimeCourses(data_reg, mask_cell); %applies mask to stack (averages all pixels in each frame for each cell) to get timecourses
            %Timecourses are nFrames x nCells
fprintf(['data_tc is ' num2str(size(data_tc))]) 
[nFrames, nCells] = size(data_tc);
%        1.  Downsampled timecourses for neuropil subtraction- averageing decreases noise
down = 5; %number of frames to average
data_reg_down = stackGroupProject(data_reg,down); %averages every 5 frames in stack  
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);    
%        2.  Extract neuropil timecourses (full and downsampled)
np_tc = zeros(nFrames,nCells);
np_tc_down = zeros(floor(nFrames./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '%s\n']) 
end
%        3. Find best neuropil weights by maximizing skew on downsampled subtractions
%            Assumes calcium signals are 1) sparse and 2) positive.  Too little subtraction will decrease sparseness and too much will make signals negative. 
%            Skewness decribes the shape of a distribution- a gaussian has a skew of 0; long tail to the right is a positive skew; long tail to the left if a negative skew.  
%            Thus, the best neuropil subtraction should maximize sparseness and minimize negative values.  
%            The best subtraction will therefore yield the highest skew for the distribution of fluorescence values for each cell.
%        a. Measure skew for all weights:
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
%        b. Find index with highest skew for each cell
[max_skew, ind] =  max(x,[],1); 
%        c. convert to weight
np_w = 0.01*ind; 
%        4. Subtract weighted neuropil response from full timecourses
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
clear data_reg data_reg_down

npSub_tc_crop = npSub_tc(nOff/2+1:end-nOn-(nOff/2),:);
data_trial = reshape(npSub_tc_crop, [nOn+nOff ntrials-1 nCells]);
data_f = mean(data_trial(nOff/4:nOff/2,:,:),1);
data_dfof_tc = bsxfun(@rdivide,bsxfun(@minus,data_trial,data_f),data_f);
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc','data_dfof_tc')


%% make tuning curves
base_win = nOff/2:nOff;
resp_win = nOff+5:nOff+nOn;
data_trial = reshape(npSub_tc, [nOn+nOff ntrials nCells]);
data_f = mean(data_trial(base_win,:,:),1);
data_dfof = bsxfun(@rdivide,bsxfun(@minus,data_trial,data_f),data_f);
data_dfof_resp = squeeze(mean(data_dfof(resp_win,:,:),1));
data_dfof_base = squeeze(mean(data_dfof(base_win,:,:),1));

resp_cell = cell(nDirs,nSF);
base_cell = cell(nDirs,nSF);
data_dfof_dir_sf = zeros(nCells,nDirs,nSF,2);
data_dfof_dir = zeros(nCells,nDirs,2);
data_dfof_ori = zeros(nCells,nDirs/2,2);
data_dfof_sf = zeros(nCells,nSF,2);
h_dir_sf = zeros(nDirs,nSF,nCells);
p_dir_sf = zeros(nDirs,nSF,nCells);
for iDir = 1:nDirs
    ind1 = find(Dir==Dirs(iDir));
    for iSF = 1:nSF
        ind2 = find(SF_mat==SFs(iSF));
        ind_use = intersect(ind1,ind2);
        resp_cell{iDir,iSF} = squeeze(mean(data_dfof(resp_win,ind_use,:),1));
        base_cell{iDir,iSF} = squeeze(mean(data_dfof(base_win,ind_use,:),1));
        [h_dir_sf(iDir,iSF,:), p_dir_sf(iDir,iSF,:)] = ttest(resp_cell{iDir,iSF},base_cell{iDir,iSF},'tail','right','alpha',0.05./(nDirs-1));
        data_dfof_dir_sf(:,iDir,iSF,1) = squeeze(mean(mean(data_dfof(resp_win,ind_use,:),1),2));
        data_dfof_dir_sf(:,iDir,iSF,2) = squeeze(std(mean(data_dfof(resp_win,ind_use,:),1),[],2)./sqrt(length(ind_use)));
        data_dfof_sf(:,iSF,1) = squeeze(mean(mean(data_dfof(resp_win,ind2,:),1),2));
        data_dfof_sf(:,iSF,2) = squeeze(std(mean(data_dfof(resp_win,ind2,:),1),[],2)./sqrt(length(ind2)));
    end
    data_dfof_dir(:,iDir,1) = squeeze(mean(mean(data_dfof(resp_win,ind1,:),1),2));
    data_dfof_dir(:,iDir,2) = squeeze(std(mean(data_dfof(resp_win,ind1,:),1),[],2)./sqrt(length(ind1)));
    if iDir<= nDirs/2
        ind_ori = [ind1 find(Dir==Dirs(iDir+nDirs/2))];
        data_dfof_ori(:,iDir,1) = squeeze(mean(mean(data_dfof(resp_win,ind_ori,:),1),2));
        data_dfof_ori(:,iDir,2) = squeeze(std(mean(data_dfof(resp_win,ind_ori,:),1),[],2)./sqrt(length(ind_ori)));
    end
end

% figure;
% [n1 n2] = subplotn(nCells);
% for iC = 1:nCells
%     subplot(n1,n2,iC)
%     imagesc(squeeze(data_dfof_dir_sf(iC,:,:,1)))
%     title([num2str(chop(max(max(data_dfof_dir_sf(iC,:,:,1),[],2),[],3),2)) ' ' num2str(sum(sum(h_dir_sf(:,:,iC))))])
%     set(gca,'XTickLabel',SFs,'YTickLabel',Dirs(2:2:end))
% end
% figure;
% [n1 n2] = subplotn(nCells);
% for iC = 1:nCells
%     subplot(n1,n2,iC)
%     for iDir = 1:nDirs
%         errorbar(SFs,squeeze(data_dfof_dir_sf(iC,iDir,:,1)),squeeze(data_dfof_dir_sf(iC,iDir,:,2)))
%         hold on
%     end
% end
[val prefSF] = max(data_dfof_sf(:,:,1),[],2);
[prefVal prefDir] = max(data_dfof_dir(:,:,1),[],2);
oppDir = prefDir+nDirs/2;
oppDir(find(oppDir>nDirs)) = oppDir(find(oppDir>nDirs))-nDirs;
oppVal = indOnly(data_dfof_dir(:,:,1),oppDir);
prefVal(find(prefVal<0)) = 0;
oppVal(find(oppVal<0)) = 0;
DSI = (prefVal-oppVal)./(prefVal+oppVal);

% figure;
% [n1 n2] = subplotn(nCells);
% for iC = 1:nCells
%     subplot(n1,n2,iC)
%     errorbar(SFs,squeeze(data_dfof_sf(iC,:,1)),squeeze(data_dfof_sf(iC,:,2)))
%     hold on
%     errorbar(SFs,squeeze(data_dfof_dir_sf(iC,prefDir(iC),:,1)),squeeze(data_dfof_dir_sf(iC,prefDir(iC),:,2)))
%     title(num2str(SFs(prefSF(iC))))
% end
% 
% figure;
% [n1 n2] = subplotn(nCells);
% for iC = 1:nCells
%     subplot(n1,n2,iC)
%     for iSF = 1:nSF
%         errorbar(Dirs,data_dfof_dir_sf(iC,:,iSF,1),data_dfof_dir_sf(iC,:,iSF,2))
%         hold on
%     end
% end
% figure;
% [n1 n2] = subplotn(nCells);
% for iC = 1:nCells
%     subplot(n1,n2,iC)
%     errorbar(Dirs,squeeze(data_dfof_dir(iC,:,1)),squeeze(data_dfof_dir(iC,:,2)))
%     hold on
%     errorbar(Dirs,squeeze(data_dfof_dir_sf(iC,:,prefSF(iC),1)),squeeze(data_dfof_dir_sf(iC,:,prefSF(iC),2)))
% end
% 
% figure;
% [n1 n2] = subplotn(nCells);
% for iC = 1:nCells
%     subplot(n1,n2,iC)
%     errorbar(Dirs(1:nDirs/2),squeeze(data_dfof_ori(iC,:,1)),squeeze(data_dfof_ori(iC,:,2)))
% end
[prefVal prefOri] = max(data_dfof_ori(:,:,1),[],2);
nOri = nDirs/2;
orthOri = prefOri+nOri/2;
orthOri(find(orthOri>nOri)) = orthOri(find(orthOri>nOri))-nOri;
orthVal = indOnly(data_dfof_ori(:,:,1),orthOri);
prefVal(find(prefVal<0)) = 0;
orthVal(find(orthVal<0)) = 0;
OSI = (prefVal-orthVal)./(prefVal+orthVal);
figure;
subplot(2,2,1)
histogram(DSI,[0:0.1:1])
xlim([0 1])
xlabel('DSI')
subplot(2,2,2)
histogram(OSI,[0:0.1:1])
xlim([0 1])
xlabel('OSI')
subplot(2,2,3)
histogram(prefSF)
set(gca,'XtickLabel',SFs)
xlabel('Pref SF')
sgtitle([mouse ' ' date])
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_DSI_OSI_prefSF.pdf']),'-dpdf','-bestfit')

data_dfof_stim_dir = zeros(nDirs,nSF,nCells);
h_dir = zeros(nDirs,nSF,nCells);
p_dir = zeros(nDirs,nSF,nCells);
trialInd = cell(nDirs,nSF);
for iDir = 1:nDirs
    ind_dir = find(Dir == Dirs(iDir));
    for iSF = 1:nSF
        ind_SF = find(SF_mat == SFs(iSF));
        trialInd{iDir,iSF} = intersect(ind_dir,ind_SF);
        data_dfof_stim_dir(iDir, iSF, :) = mean(data_dfof_resp(trialInd{iDir,iSF},:),1);
        [h_dir(iDir,iSF,:) p_dir(iDir,iSF,:)] = ttest(data_dfof_resp(trialInd{iDir,iSF},:), data_dfof_base(trialInd{iDir,iSF},:),'dim',1, 'tail', 'right', 'alpha', 0.05./((nDirs.*nSF)-1));
    end
end

h_dir_all = find(sum(sum(h_dir,1),2));

[max_dir_val max_dir_ind] = max(squeeze(mean(data_dfof_stim_dir,2)),[],1);
max_sf_val = zeros(1,nCells);
max_sf_ind = zeros(1,nCells);
fit_out = cell(1,nCells);
g_fit = cell(1,nCells);
prefSF = zeros(1,nCells);
figure; movegui('center')
start = 1;
if nCells<49
    [n n2] = subplotn(nCells);
else
    n = 7;
    n2 = 7;
end
for iCell = 1:nCells
    if start >49
        figure;
        start = 1;
    end
    [max_sf_val(1,iCell) max_sf_ind(1,iCell)] = max(data_dfof_stim_dir(max_dir_ind(iCell),:,iCell),[],2);
    data_use = mean(data_dfof_stim_dir(find(h_dir(:,max_sf_ind(iCell),iCell)),:,iCell),1);
    if ~isnan(sum(data_use))
        [s i] = sort(data_use);
        if abs(i(nSF)-i(nSF-1))<3  || s(nSF-1)./s(nSF)<0.4
            x = log2([0.01 SFs 0.64]);
            options = fitoptions('gauss1', 'Lower', [s(end) x(i(nSF)) 0], 'Upper', [s(end).*1.5 x(i(nSF)+2) Inf]);
            [fit_out{iCell} g_fit{iCell}] =fit(log2(SFs)',data_use','gauss1',options);
        else
            fit_out{iCell}.b1 = NaN;
            g_fit{iCell}.rsquare = NaN;
        end
        subplot(n,n2,start)
        plot(log2(SFs),data_use,'o')
        hold on
        if ~isnan(fit_out{iCell}.b1)
            plot(fit_out{iCell})
        end
        prefSF(1,iCell)= fit_out{iCell}.b1;
        RsqSF(1,iCell)= g_fit{iCell}.rsquare;
        if find(h_dir_all == iCell)
            title('Sig')
        end
        legend('off')
        if rem(iCell, 10) == 0
            fprintf([num2str(iCell) '\n'])
        end
    else
        fit_out{iCell}.b1 = nan;
        prefSF(1,iCell)= nan;
        RsqSF(1,iCell)= nan;
    end
    start = start+1;
end

figure;
start = 1;
for iCell = 1:nCells
    if start >49
        figure;
        start = 1;
    end
    subplot(n,n2,start)
    imagesc(data_dfof_stim_dir(:,:,iCell))
    start = start+1;
end

prefSF_cut = prefSF;
prefSF_cut(find(prefSF>max(log2(SFs)))) = max(log2(SFs));
prefSF_cut(find(prefSF<min(log2(SFs)))) = min(log2(SFs));
figure; movegui('center');
RsqSF(find(RsqSF<0)) = 0;
subplot(2,1,1); histogram(RsqSF(h_dir_all),[0:.1:1.1]); vline(0.5)
ind_sf = intersect(find(RsqSF>0.5),h_dir_all); xlabel('Rsq')
subplot(2,1,2); histogram(prefSF_cut(ind_sf))
set(gca, 'XTick', log2(SFs), 'XTickLabels', SFs)
vline(mean(prefSF_cut(ind_sf),2))
xlabel('SF (cpd)')
suptitle([date ' ' mouse])
print(fullfile(fnout, datemouse, datemouserun, [datemouserun  '_SFfitDist.pdf']), '-dpdf')


save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_respData.mat']),'resp_win','base_win','resp_cell','base_cell','data_dfof_dir_sf','data_dfof_dir','data_dfof_sf','data_dfof_ori','h_dir_sf','prefSF','DSI','OSI')
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_SFfits.mat']), 'fit_out','g_fit', 'prefSF', 'RsqSF','data_dfof_resp','data_dfof_base','h_dir_all', 'h_dir', 'trialInd','max_sf_ind','max_dir_ind')

data_dfof_stim_dir_shift = zeros(size(data_dfof_stim_dir));
for iCell = 1:nCells
    data_dfof_stim_dir_shift(:,:,iCell) = circshift(data_dfof_stim_dir(:,:,iCell),nDirs/2-max_dir_ind(iCell),1);
end
figure; 
subplot(2,1,1)
imagesc(mean(data_dfof_stim_dir_shift,3))
set(gca,'YTickLabels',Dirs(2:2:end), 'XTick', 1:nSF, 'XTickLabels', SFs)
subplot(2,1,2)
imagesc(mean(data_dfof_stim_dir_shift./max(max(data_dfof_stim_dir_shift,1),2),3))
set(gca,'YTickLabels',Dirs(2:2:end), 'XTick', 1:nSF, 'XTickLabels', SFs)
title('Normalized')
suptitle([date ' ' mouse])
print(fullfile(fnout, datemouse, datemouserun, [datemouserun  '_avgSFxDir.pdf']), '-dpdf')
movegui('center')
%% wheel data
if ~exist(fullfile(fnout,datemouse,datemouserun,[datemouserun '_wheelSpeed.mat']))
    wheel_speed = wheelSpeedCalc(input,32,'orange');
    nC = size(npSub_tc,2);
    r_wheel = zeros(nC,1);
    for iC = 1:nC
        r_wheel(iC,1) = triu2vec(corrcoef(smooth(npSub_tc(:,iC)',10), smooth(wheel_speed,10)));
    end
    save(fullfile(fnout,datemouse,datemouserun,[datemouserun '_wheelSpeed.mat']),'wheel_speed','r_wheel')
else
    load(fullfile(fnout,datemouse,datemouserun,[datemouserun '_wheelSpeed.mat']))
    nC = size(npSub_tc,2);
    r_wheel = zeros(nC,1);
    for iC = 1:nC
        r_wheel(iC,1) = triu2vec(corrcoef(smooth(npSub_tc(:,iC)',10), smooth(wheel_speed,10)));
    end
    save(fullfile(fnout,datemouse,datemouserun,[datemouserun '_wheelSpeed.mat']),'wheel_speed','r_wheel')
end

%% pupil data
if ~exist(fullfile(fnout,datemouse,datemouserun,[datemouserun '_pupil.mat']))
    data = load(fullfile(data_fn,mouse,date,ImgFolder,[ImgFolder '_000_000_eye.mat']));
    data=squeeze(data.data);
    [data_crop rect] = cropEyeData(data); 
    eyeData = extractEyeData(data_crop,[4 25]);
    pupil_tc = sqrt(eyeData.Area/pi);
    pupil_tc = pupil_tc(1:size(npSub_tc,1));
else
    load(fullfile(fnout,datemouse,datemouserun,[datemouserun '_pupil.mat']))
end

pupil_sm = smooth(pupil_tc,10);
nC = size(npSub_tc,2);
r_pupil = zeros(nC,1);
for iC = 1:nC
    r_pupil(iC,1) = triu2vec(corrcoef(smooth(npSub_tc(:,iC)',10), pupil_sm));
end

pos_tc = eyeData.Centroid;
pos_diff = diff(pos_tc(:,1));
figure; plot(pos_diff)
std_pos = std(pos_diff,[],1,'omitnan');
ind = find(abs(pos_diff)>2.*std_pos);
ind_use = ind;
ind_use(find(diff(ind)<3)+1) = [];
hold on; 
plot(ind_use',pos_diff(ind_use,1),'or')
hold on; plot((pos_tc(:,1)-mean(pos_tc(:,1)))./max(pos_tc(:,1)))

s_pos = [];
s_neg = [];
s_pos_vis = [];
s_neg_vis = [];
nframes = size(pos_diff,1);
ind_use(find(ind_use<30))= [];
ind_use(find(ind_use>nframes-30))= [];
ind_mat = reshape(1:nframes,[nOn+nOff ntrials]);
ind_on = reshape(ind_mat(nOff+1:end,:),[nOn.*ntrials 1]);
ind_off = reshape(ind_mat(1:nOff,:),[nOff.*ntrials 1]);
for i = 1:length(ind_use)
    if pos_diff(ind_use(i))>0
        if find(ind_on==ind_use(i))
            s_pos_vis = [s_pos_vis ind_use(i)];
        else
            s_pos = [s_pos ind_use(i)];
        end
    else
        if find(ind_on==ind_use(i))
            s_neg_vis = [s_neg_vis ind_use(i)];
        else
            s_neg = [s_neg ind_use(i)];
        end
    end
end
figure; 
nCells = size(npSub_tc,2);
npSub_f = prctile(npSub_tc,10,1);
npSub_dfof = (npSub_tc-npSub_f)./npSub_f;
pos_tc_base = median(pos_tc(:,1));
pos_tc_delta = (pos_tc(:,1)-pos_tc_base);
s_pos_eyetc = zeros(60,length(s_pos));
s_pos_celltc = zeros(60,nCells,length(s_pos));
for i = 1:length(s_pos)
    s_pos_eyetc(:,i) = pos_tc_delta(s_pos(i)-29:s_pos(i)+30,1);
    s_pos_celltc(:,:,i) = npSub_dfof(s_pos(i)-29:s_pos(i)+30,:);
end
tt = (-29:30).*1./frame_rate;
subplot(2,2,1)
%plot(tt,mean(s_pos_eyetc,2)-mean(mean(s_pos_eyetc,2),1)); 
hold on; 
plot(tt,mean(s_pos_celltc,3)-mean(mean(s_pos_celltc,1),3))
vline(0)
title('Pos- off')
ylabel('dF/F')
xlabel('Time from eye movement (s)')

s_pos_vis_eyetc = zeros(60,length(s_pos_vis));
s_pos_vis_celltc = zeros(60,nCells,length(s_pos_vis));
for i = 1:length(s_pos_vis)
    s_pos_vis_eyetc(:,i) = pos_tc_delta(s_pos_vis(i)-29:s_pos_vis(i)+30,1);
    s_pos_vis_celltc(:,:,i) = npSub_dfof(s_pos_vis(i)-29:s_pos_vis(i)+30,:);
end
subplot(2,2,2)
%plot(tt,mean(s_pos_vis_eyetc,2)-mean(mean(s_pos_vis_eyetc,2),1)); 
hold on; 
plot(tt,mean(s_pos_vis_celltc,3)-mean(mean(s_pos_vis_celltc,1),3))
title('Pos- on')
ylabel('dF/F')
xlabel('Time from eye movement (s)')

s_neg_eyetc = zeros(60,length(s_neg));
s_neg_celltc = zeros(60,nCells,length(s_neg));
for i = 1:length(s_neg)
    s_neg_eyetc(:,i) = pos_tc_delta(s_neg(i)-29:s_neg(i)+30,1);
    s_neg_celltc(:,:,i) = npSub_dfof(s_neg(i)-29:s_neg(i)+30,:);
end
subplot(2,2,3)
%plot(tt,mean(s_neg_eyetc,2)-mean(mean(s_neg_eyetc,2),1)); 
hold on; 
plot(tt,mean(s_neg_celltc,3)-mean(mean(s_neg_celltc,1),3))
vline(0)
title('Neg- off')
ylabel('dF/F')
xlabel('Time from eye movement (s)')

s_neg_vis_eyetc = zeros(60,length(s_neg_vis));
s_neg_vis_celltc = zeros(60,nCells,length(s_neg_vis));
for i = 1:length(s_neg_vis)
    s_neg_vis_eyetc(:,i) = pos_tc_delta(s_neg_vis(i)-29:s_neg_vis(i)+30,1);
    s_neg_vis_celltc(:,:,i) = npSub_dfof(s_neg_vis(i)-29:s_neg_vis(i)+30,:);
end
subplot(2,2,4)
%plot(tt,mean(s_neg_vis_eyetc,2)-mean(mean(s_neg_vis_eyetc,2),1)); 
hold on; 
plot(tt,mean(s_neg_vis_celltc,3)-mean(mean(s_neg_vis_celltc,1),3))
vline(0)
title('Neg- on')
ylabel('dF/F')
xlabel('Time from eye movement (s)')
sgtitle([mouse ' ' date])
print(fullfile(fnout,datemouse,datemouserun,[datemouserun '_eyeMvmt.pdf']),'-dpdf','-bestfit')

save(fullfile(fnout,datemouse,datemouserun,[datemouserun '_pupil.mat']),'eyeData','r_pupil','rect','pupil_tc','s_neg_vis_celltc','s_neg_vis_eyetc','s_pos_vis_celltc','s_pos_vis_eyetc','s_neg_celltc','s_neg_eyetc','s_pos_celltc','s_pos_eyetc','ind_use')