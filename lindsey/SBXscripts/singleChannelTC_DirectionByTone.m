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
date = '221104';
ImgFolder = '001';
time = '1230';
mouse = 'i1377';
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
fprintf(['Data is ' num2str(size(data)) '\n'])
%When imaging a single channel, nPMT = 1, so squeeze:
data = squeeze(data);

%% Register 2P data

if exist(fullfile(fnout, datemouse, datemouserun, [datemouserun '_reg_shifts.mat']))
    load(fullfile(fnout, datemouse, datemouserun, [datemouserun '_reg_shifts.mat']))
    [out, data_reg] = stackRegister_MA(data,[],[],out);
    %New average image after registration
    data_reg_avg = mean(data_reg,3);
else
%Goal here is to remove X-Y movement artifacts
%1. Find a stable target
%    a. Plot average of 500 frames throughout stack
nframes = 500; %nframes to average for target
nskip = 5000; %nframes to skip for each average

nep = floor(size(data,3)./nskip);
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
    allAxes = findobj(f.Children,'Type','axes');
    numClicked = find(axesClicked==flipud(allAxes));
    close all
end
fprintf(['Selected subplot ' num2str(numClicked) '\n'])
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
figure;
for i = 1:length(ind) 
    subplot(2,1,i); 
    ix = ind(i);
    imagesc(mean(data_reg(:,:,1+((ix-1)*nskip):nframes+((ix-1)*nskip)),3)); 
    title([num2str(1+((ix-1)*nskip)) '-' num2str(nframes+((ix-1)*nskip))]); 
end
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_FOV_first&last.pdf']), '-dpdf')
%    b. Average of all frames should be sharp
figure;
imagesq(data_reg_avg); 
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_FOV_avg.pdf']), '-dpdf')
end
clear data
%% Segment 2P data
%Goal here is to create a cell mask to extract fluorescence timecourses
%1. Find activated cells
%    a. Reshape data stack to segregate trials 
%    First need to know how many frames per trial and how many trials

nOn = unique(celleqel2mat_padded(input.nStimOneFramesOn));
nOff = unique(celleqel2mat_padded(input.tItiWaitFrames));
ntrials = size(input.tStimOneGratingDirectionDeg,2); %this is a cell array with one value per trial, so length = ntrials
cStimOneOn = celleqel2mat_padded(input.cStimOneOn);
sz = size(data_reg);
%        Each trial has nOff frames followed by nOn frames, so can reshape stack so nYpix x nXpix x nFrames/Trial (nOn+nOff) x nTrials
data_tr = zeros(sz(1),sz(2),nOn+nOff,ntrials);
for itrial = 1:ntrials
    data_tr(:,:,:,itrial) = data_reg(:,:,cStimOneOn(itrial)-nOff+1:cStimOneOn(itrial)+nOn);
end

data_tr_inv = mean(mean(data_tr(:,:,nOff/2:nOff,:),3),4)-mean(mean(data_tr(:,:,nOff+1:nOn+nOff,:),3),4);

fprintf(['Size of data_tr is ' num2str(size(data_tr))])
%    b. Find baseline F from last half of off period- avoids decay of previous on trial
data_f = mean(data_tr(:,:,nOff/2:nOff,:),3); 
%    c. Find dF/F for each trial
data_df = bsxfun(@minus, double(data_tr), data_f); 
clear data_tr
data_dfof = bsxfun(@rdivide,data_df, data_f); 
clear data_f data_df 
%    d. Find average dF/F for each stimulus condition (this is for an experiment with changing grating direction)
tDir = celleqel2mat_padded(input.tStimOneGratingDirectionDeg); %transforms cell array into matrix (1 x ntrials)
Dirs = unique(tDir);
nDirs = length(Dirs);

tCon = celleqel2mat_padded(input.tStimOneGratingContrast);
tAud = ~celleqel2mat_padded(input.tBlock2TrialNumber);

data_dfof_avg = zeros(sz(1),sz(2),nDirs+1); %create empty matrix with FOV for each direction: nYpix x nXPix x nDir

nStim = nDirs;
figure; movegui('center')
[n, n2] = subplotn(nDirs); %function to optimize subplot number/dimensions
for idir = 1:nDirs
    ind = intersect(find(tCon),find(tDir == Dirs(idir))); %find all trials with each direction
    data_dfof_avg(:,:,idir) = mean(mean(data_dfof(:,:,nOff+1:nOn+nOff,ind),3),4); %average all On frames and all trials
    subplot(n,n2,idir)
    imagesc(data_dfof_avg(:,:,idir))
end
ind = intersect(find(tCon==0),find(tAud)); %find all trials with tone only
data_dfof_avg(:,:,idir+1) = mean(mean(data_dfof(:,:,nOff+1:nOn+nOff,ind),3),4); %average all On frames and all trials
figure; imagesc(data_dfof_avg(:,:,idir+1))
title('Tone Resp')

%        Filtering data helps make cells more visible for selection
myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_avg_all = imfilter(data_dfof_avg,myfilter);
data_dfof_max = max(data_dfof_avg_all,[],3); %finds all active cells by taking max projection
figure; movegui('center'); imagesc(data_dfof_max);

data_dfof = cat(3,data_reg_avg,cat(3,data_dfof_max,cat(3,data_dfof_avg_all,data_tr_inv)));

save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg_all', 'nStim')
%    2. Create cell masks from active cells
%        Concatenate images of max and all directions
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

save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

%% make tuning curves
base_win = nOff/2:nOff;
resp_win = nOff+5:nOff+nOn;

data_trial = zeros(nOn+nOff,nCells,ntrials);
for itrial = 1:ntrials
    data_trial(:,:,itrial) = npSub_tc(cStimOneOn(itrial)-nOff+1:cStimOneOn(itrial)+nOn,:);
end
data_trial = permute(data_trial, [1 3 2]);
data_f = mean(data_trial(base_win,:,:),1);
data_dfof = bsxfun(@rdivide,bsxfun(@minus,data_trial,data_f),data_f);

resp_cell_dir = cell(nDirs+1,2);
base_cell_dir = cell(nDirs+1,2);
ind = cell(nDirs+1,2);
indn = zeros(nDirs+1,2);
data_dfof_dir = zeros(nCells,nDirs+1,2,2);
h_dir = zeros(nDirs+1,2,nCells);
p_dir = zeros(nDirs+1,2,nCells);
for iA = 1:2
    for iDir = 1:nDirs
        ind{iDir,iA} = intersect(find(tAud == iA-1),intersect(find(tCon),find(tDir==Dirs(iDir))));
        indn(iDir,iA) = length(ind{iDir,iA});
        resp_cell_dir{iDir,iA} = squeeze(mean(data_dfof(resp_win,ind{iDir,iA},:),1));
        base_cell_dir{iDir,iA} = squeeze(mean(data_dfof(base_win,ind{iDir,iA},:),1));
        [h_dir(iDir,iA,:), p_dir(iDir,iA,:)] = ttest(resp_cell_dir{iDir},base_cell_dir{iDir},'tail','right','alpha',0.05./(nDirs-1));
        data_dfof_dir(:,iDir,iA,1) = squeeze(mean(mean(data_dfof(resp_win,ind{iDir,iA},:),1),2));
        data_dfof_dir(:,iDir,iA,2) = squeeze(std(mean(data_dfof(resp_win,ind{iDir,iA},:),1),[],2)./sqrt(length(ind{iDir,iA})));
    end
    ind{iDir+1,iA} = intersect(find(tAud == iA-1),find(tCon==0));
    indn(iDir+1,iA) = length(ind{iDir+1,iA});
    resp_cell_dir{iDir+1,iA} = squeeze(mean(data_dfof(resp_win,ind{iDir+1,iA},:),1));
    base_cell_dir{iDir+1,iA} = squeeze(mean(data_dfof(base_win,ind{iDir+1,iA},:),1));
    [h_dir(iDir+1,iA,:), p_dir(iDir+1,iA,:)] = ttest(resp_cell_dir{iDir},base_cell_dir{iDir},'tail','right','alpha',0.05./(nDirs-1));
    data_dfof_dir(:,iDir+1,iA,1) = squeeze(mean(mean(data_dfof(resp_win,ind{iDir+1,iA},:),1),2));
    data_dfof_dir(:,iDir+1,iA,2) = squeeze(std(mean(data_dfof(resp_win,ind{iDir+1,iA},:),1),[],2)./sqrt(length(ind{iDir+1,iA})));
end

h_all_dir = sum(squeeze(h_dir(:,1,:)>0),1);
resp_ind = find(h_all_dir);
start=1;
n = 1;
figure;
movegui('center')
for iCell = 1:nCells
    if start>25
        print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_cellTuningDir' num2str(n) '.pdf']),'-dpdf','-bestfit')
        figure;movegui('center');
        start = 1;
        n = n+1;
    end
    subplot(5,5,start)
    errorbar(Dirs, data_dfof_dir(iCell,1:nDirs,1,1), data_dfof_dir(iCell,1:nDirs,1,2), '-o')
    hold on 
    errorbar(Dirs, data_dfof_dir(iCell,1:nDirs,2,1), data_dfof_dir(iCell,1:nDirs,2,2), '-o')
    errorbar(360, data_dfof_dir(iCell,end,1,1), data_dfof_dir(iCell,end,1,2), '-o')
    errorbar(360, data_dfof_dir(iCell,end,2,1), data_dfof_dir(iCell,end,2,2), '-o')
    title(['R = ' num2str(h_all_dir(iCell))])
    start = start +1;
end
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_cellTuningDir' num2str(n) '.pdf']),'-dpdf','-bestfit')

tOri = tDir;
tOri(find(tDir>=180)) = tOri(find(tDir>=180))-180;
Oris = unique(tOri);
nOri = length(Oris);

resp_cell_ori = cell(nOri,2);
base_cell_ori = cell(nOri,2);
data_dfof_ori = zeros(nCells,nOri,2,2);
h_ori = zeros(nOri,2,nCells);
p_ori = zeros(nOri,2,nCells);
for iA = 1:2
    for iOri = 1:nOri
        ind = intersect(find(tAud == iA-1),intersect(find(tCon),find(tOri==Oris(iOri))));
        resp_cell_ori{iOri,iA} = squeeze(mean(data_dfof(resp_win,ind,:),1));
        base_cell_ori{iOri,iA} = squeeze(mean(data_dfof(base_win,ind,:),1));
        [h_ori(iOri,iA,:) p_ori(iOri,iA,:)] = ttest(resp_cell_ori{iOri},base_cell_ori{iOri},'tail','right','alpha',0.05./(nOri-1));
        data_dfof_ori(:,iOri,iA,1) = squeeze(mean(mean(data_dfof(resp_win,ind,:),1),2));
        data_dfof_ori(:,iOri,iA,2) = squeeze(std(mean(data_dfof(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
    end
end

h_all_ori = sum(squeeze(h_ori(:,1,:)>0),1);

start=1;
n = 1;
figure;
movegui('center')
for iCell = 1:nCells
    if start>25
        print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_cellTuningOri' num2str(n) '.pdf']),'-dpdf','-bestfit')
        figure;movegui('center');
        start = 1;
        n = n+1;
    end
    subplot(5,5,start)
    errorbar(Oris, data_dfof_ori(iCell,:,1,1), data_dfof_ori(iCell,:,1,2), '-o')
    hold on
    errorbar(Oris, data_dfof_ori(iCell,:,2,1), data_dfof_ori(iCell,:,2,2), '-o')
    title(['R = ' num2str(h_all_ori(iCell))])
    start = start +1;
end
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_cellTuningOri' num2str(n) '.pdf']),'-dpdf','-bestfit')

 b_ori = nan(2,nCells);
    k1_ori = nan(2,nCells);
    R1_ori = nan(2,nCells);
    u1_ori = nan(2,nCells);
    R_square_ori = nan(2,nCells);
    sse_ori = nan(2,nCells);
    stim_DSI = nan(2,nCells);
    stim_OSI = nan(2,nCells);
    
    Oris = Dirs(1:nDirs/2);
    nOris = length(Oris);
    data_dfof_ori = squeeze(mean(reshape(data_dfof_dir(:,1:nDirs,:,1),[nCells nOris 2 2]),3));
    for iCell = 1:nCells
        for i = 1
            data = [data_dfof_ori(iCell,:,i) data_dfof_ori(iCell,1,i)];
            theta = [deg2rad(Oris) pi];
            [b_ori(i,iCell),k1_ori(i,iCell),R1_ori(i,iCell),u1_ori(i,iCell),sse_ori(i,iCell),R_square_ori(i,iCell)] ...
                = miaovonmisesfit_ori(theta,data);
            [max_val max_ind] = max(data_dfof_ori(iCell,:,i),[],2);
            null_ind = max_ind+(nOris./2);
            null_ind(find(null_ind>nOris)) = null_ind(find(null_ind>nOris))-nOris;
            min_val = data_dfof_ori(iCell,null_ind,i);
            if min_val<0
                min_val = 0;
            end
            stim_OSI(i,iCell) = (max_val-min_val)./(max_val+min_val);
            [max_val max_ind] = max(data_dfof_dir(iCell,1:nDirs,i,1),[],2);
            null_ind = max_ind+(nDirs./2);
            null_ind(find(null_ind>nDirs)) = null_ind(find(null_ind>nDirs))-nDirs;
            min_val = data_dfof_dir(iCell,null_ind,i,1);
            if min_val<0
                min_val = 0;
            end
            stim_DSI(i,iCell) = (max_val-min_val)./(max_val+min_val);
        end
    end
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_oriResp.mat']), 'data_dfof_dir', 'data_dfof_ori','base_win','resp_win','h_dir','k1_ori','b_ori','R1_ori','u1_ori','stim_DSI','stim_OSI','R_square_ori')

[max_val_ori max_ori] = max(mean(data_dfof_ori,3),[],2);
[max_val_dir max_dir] = max(mean(data_dfof_dir(:,1:nDirs,:,1),3),[],2);
data_dfof_ori_shift = data_dfof_ori;
data_dfof_dir_shift = data_dfof_dir(:,1:nDirs,:,1);
for iCell = 1:nCells
    data_dfof_ori_shift(iCell,:,:) = circshift(data_dfof_ori(iCell,:,:),-max_ori(iCell)+4,2);
    data_dfof_dir_shift(iCell,:,:) = circshift(data_dfof_dir(iCell,1:nDirs,:,1),-max_dir(iCell)+8,2);
end

data_dfof_ori_shift_norm = data_dfof_ori_shift./max_val_ori;
data_dfof_dir_shift_norm = data_dfof_dir_shift./max_val_dir;

if length(resp_ind)>1
figure; 
subplot(2,2,1)
errorbar(Oris,mean(data_dfof_ori_shift_norm(resp_ind,:,1),1),std(data_dfof_ori_shift_norm(resp_ind,:,1),1)./sqrt(length(resp_ind)))
hold on
errorbar(Oris,mean(data_dfof_ori_shift_norm(resp_ind,:,2),1),std(data_dfof_ori_shift_norm(resp_ind,:,2),1)./sqrt(length(resp_ind)))
subplot(2,2,2)
errorbar(Dirs,mean(data_dfof_dir_shift_norm(resp_ind,:,1),1),std(data_dfof_dir_shift_norm(resp_ind,:,1),1)./sqrt(length(resp_ind)))
hold on
errorbar(Dirs,mean(data_dfof_dir_shift_norm(resp_ind,:,2),1),std(data_dfof_dir_shift_norm(resp_ind,:,2),1)./sqrt(length(resp_ind)))
subplot(2,2,3)
errorbar(1,mean(data_dfof_dir(resp_ind,end,1),1),std(data_dfof_dir(resp_ind,end,1),1)./sqrt(length(resp_ind)))
hold on
errorbar(1,mean(data_dfof_dir(resp_ind,end,2),1),std(data_dfof_dir(resp_ind,end,2),1)./sqrt(length(resp_ind)))
end

%%
data = load(fullfile(data_fn,mouse,date,ImgFolder,[ImgFolder '_000_000_eye.mat']));
data=squeeze(data.data);
[data_crop rect] = cropEyeData(data); 
eyeData = extractEyeData(data_crop,[3 25]);
pupil_tc = sqrt(eyeData.Area/pi);
pupil_tc = pupil_tc(1:size(npSub_tc,1));

pupil_trial = zeros(nOn+nOff,ntrials);
for itrial = 1:ntrials
    pupil_trial(:,itrial) = pupil_tc(cStimOneOn(itrial)-nOff+1:cStimOneOn(itrial)+nOn,:);
end
pupil_resp = zeros(nOn+nOff,2);
pupil_resp(:,1) = mean(pupil_trial(:,find(~tAud)),2,'omitnan');
pupil_resp(:,2) = mean(pupil_trial(:,find(tAud)),2,'omitnan');
save(fullfile(fnout,datemouse,datemouserun,[datemouserun '_pupil.mat']),'eyeData','rect','pupil_trial','pupil_resp')

%%
wheel_speed = wheelSpeedCalc(input,32,'orange');
figure; plot(wheel_speed)
wheel_trial = zeros(nOn+nOff,ntrials);
for itrial = 1:ntrials
   wheel_trial(:,itrial) = wheel_speed(:,cStimOneOn(itrial)-nOff+1:cStimOneOn(itrial)+nOn);
end
figure; plot(mean(wheel_trial,1))
save(fullfile(fnout,datemouse,datemouserun,[datemouserun '_wheel.mat']),'wheel_speed','wheel_trial')