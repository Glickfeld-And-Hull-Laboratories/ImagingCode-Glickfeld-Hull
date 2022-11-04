%% Load, register, segment and neuropil correct 2P data
close all
clear all global

%Path names
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
cam_fn = fullfile(fn_base, 'home\camaron');
lg_fn = fullfile(fn_base, 'home\lindsey');
data_fn = fullfile(lg_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data');
fnout = fullfile(lg_fn, 'Analysis\2P');

%Specific experiment information
date = '220904';
ImgFolder = '002';
time = '1343';
mouse = 'i1372';
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
    numClicked = find(axesClicked==allAxes);
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
fprintf(['Size of data_tr is ' num2str(size(data_tr))])
%    b. Find baseline F from last half of off period- avoids decay of previous on trial
data_f = mean(data_tr(:,:,nOff/2:nOff,:),3); 
%    c. Find dF/F for each trial
data_df = bsxfun(@minus, double(data_tr), data_f); 
data_dfof = bsxfun(@rdivide,data_df, data_f); 
clear data_f data_df data_tr
%    d. Find average dF/F for each stimulus condition (this is for an experiment with changing grating direction)
tDir = celleqel2mat_padded(input.tDotDirectionDeg); %transforms cell array into matrix (1 x ntrials)
Dirs = unique(tDir);
nDirs = length(Dirs);
data_dfof_avg = zeros(sz(1),sz(2),nDirs+1); %create empty matrix with FOV for each direction: nYpix x nXPix x nDir
tCoh = celleqel2mat_padded(input.tDotCoherence);
nStim = nDirs;
figure; movegui('center')
[n, n2] = subplotn(nDirs+1); %function to optimize subplot number/dimensions
for idir = 1:nDirs
    ind = intersect(find(tCoh == 1), find(tDir == Dirs(idir))); %find all trials with each direction
    data_dfof_avg(:,:,idir) = mean(mean(data_dfof(:,:,nOff+1:nOn+nOff,ind),3),4); %average all On frames and all trials
    subplot(n,n2,idir)
    imagesc(data_dfof_avg(:,:,idir))
end
ind = find(tCoh < 1); %find all trials with each direction
data_dfof_avg(:,:,idir+1) = mean(mean(data_dfof(:,:,nOff+1:nOn+nOff,ind),3),4); %average all On frames and all trials
subplot(n,n2,idir+1)
imagesc(data_dfof_avg(:,:,idir+1))

clear data_dfof
%        Filtering data helps make cells more visible for selection
myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_avg_all = imfilter(data_dfof_avg,myfilter);
data_dfof_max = max(data_dfof_avg_all,[],3); %finds all active cells by taking max projection
figure; movegui('center'); imagesc(data_dfof_max);
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

save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

%% make tuning curves
base_win = nOff/2:nOff;
resp_win = nOff+5:nOff+nOn;
data_trial = reshape(npSub_tc, [nOn+nOff ntrials nCells]);
data_f = mean(data_trial(base_win,:,:),1);
data_dfof = bsxfun(@rdivide,bsxfun(@minus,data_trial,data_f),data_f);

resp_cell_dir = cell(1,nDirs+1);
base_cell_dir = cell(1,nDirs+1);
ind = cell(1,nDirs+1);
data_dfof_dir = zeros(nCells,nDirs+1,2);
h_dir = zeros(nDirs+1,nCells);
p_dir = zeros(nDirs+1,nCells);
for iDir = 1:nDirs
    ind{iDir} = intersect(find(tCoh == 1),find(tDir==Dirs(iDir)));
    resp_cell_dir{iDir} = squeeze(mean(data_dfof(resp_win,ind{iDir},:),1));
    base_cell_dir{iDir} = squeeze(mean(data_dfof(base_win,ind{iDir},:),1));
    [h_dir(iDir,:), p_dir(iDir,:)] = ttest(resp_cell_dir{iDir},base_cell_dir{iDir},'tail','right','alpha',0.05./(nDirs-1));
    data_dfof_dir(:,iDir,1) = squeeze(mean(mean(data_dfof(resp_win,ind{iDir},:),1),2));
    data_dfof_dir(:,iDir,2) = squeeze(std(mean(data_dfof(resp_win,ind{iDir},:),1),[],2)./sqrt(length(ind{iDir})));
end
iDir = iDir+1;
ind{iDir} = find(tCoh<1);
resp_cell_dir{iDir} = squeeze(mean(data_dfof(resp_win,ind{iDir},:),1));
base_cell_dir{iDir} = squeeze(mean(data_dfof(base_win,ind{iDir},:),1));
[h_dir(iDir,:), p_dir(iDir,:)] = ttest(resp_cell_dir{iDir},base_cell_dir{iDir},'tail','right','alpha',0.05./(nDirs-1));
data_dfof_dir(:,iDir,1) = squeeze(mean(mean(data_dfof(resp_win,ind{iDir},:),1),2));
data_dfof_dir(:,iDir,2) = squeeze(std(mean(data_dfof(resp_win,ind{iDir},:),1),[],2)./sqrt(length(ind{iDir})));

group_str = {'Dir=0','Dir=180','Coh=0'};

h_all_dir = sum(h_dir,1);
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
    errorbar(1:3, data_dfof_dir(iCell,:,1), data_dfof_dir(iCell,:,2), '-o')
    set(gca, 'XTick', 1:3, 'XTickLabel', group_str)
    xlim([0 4])
    ylim([0 inf])
    title(['R = ' num2str(h_all_dir(iCell))])
    start = start +1;
end
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_cellTuningDir' num2str(n) '.pdf']),'-dpdf','-bestfit')

%%
wheel_speed = wheelSpeedCalc(input,32,'purple');
wheel_tr_tc = reshape(wheel_speed,[nOn+nOff, ntrials]);
wheel_tr = mean(wheel_tr_tc(nOff+1:end,:),1);
figure; plot(wheel_tr)
RIx = wheel_tr>1;


ind_norun_1 = find(ismember(ind{1},find(~RIx)));
ind_norun_2 = find(ismember(ind{2},find(~RIx)));
[h_ds_0, p_ds_0] = ttest2(resp_cell_dir{1}(ind_norun_1,:),resp_cell_dir{2}(ind_norun_2,:),'tail','right');
[h_ds_180, p_ds_180] = ttest2(resp_cell_dir{1}(ind_norun_1,:),resp_cell_dir{2}(ind_norun_2,:),'tail','left');

data_dfof_dir_run = zeros(nCells,nDirs+1,2);
data_dfof_dir_norun = zeros(nCells,nDirs+1,2);
for iDir = 1:nDirs+1
    run_ind = intersect(ind{iDir},find(RIx));
    length(run_ind)
    norun_ind = intersect(ind{iDir},find(~RIx));
    data_dfof_dir_run(:,iDir,1) = squeeze(mean(mean(data_dfof(resp_win,run_ind,:),1),2));
    data_dfof_dir_run(:,iDir,2) = squeeze(std(mean(data_dfof(resp_win,run_ind,:),1),[],2)./sqrt(length(run_ind)));
    data_dfof_dir_norun(:,iDir,1) = squeeze(mean(mean(data_dfof(resp_win,norun_ind,:),1),2));
    data_dfof_dir_norun(:,iDir,2) = squeeze(std(mean(data_dfof(resp_win,norun_ind,:),1),[],2)./sqrt(length(norun_ind)));
end

figure;
subplot(2,2,1)
errorbar(1:3, mean(data_dfof_dir_run(:,:,1),1), std(data_dfof_dir_run(:,:,1),[],1)./sqrt(nCells), '-o')
hold on
errorbar(1:3, mean(data_dfof_dir_norun(:,:,1),1), std(data_dfof_dir_norun(:,:,1),[],1)./sqrt(nCells), '-o')
set(gca, 'XTick', 1:3, 'XTickLabel', group_str)
xlim([0 4])
ylim([0 inf])
ylabel('dF/F')
legend({'Running','Stationary'},'location','southeast')
subplot(2,2,2)
[max_val_run max_ind_run] = max(data_dfof_dir_run(:,:,1),[],2);
data_dfof_dir_run_norm = data_dfof_dir_run(:,:,1)./max_val_run;
[max_val max_ind_norun] = max(data_dfof_dir_norun(:,:,1),[],2);
data_dfof_dir_norun_norm = data_dfof_dir_norun(:,:,1)./max_val_run;
errorbar(1:3, mean(data_dfof_dir_run_norm,1), std(data_dfof_dir_run_norm,[],1)./sqrt(nCells), '-o')
hold on
errorbar(1:3, mean(data_dfof_dir_norun_norm,1), std(data_dfof_dir_norun_norm,[],1)./sqrt(nCells), '-o')
set(gca, 'XTick', 1:3, 'XTickLabel', group_str)
xlim([0 4])
ylim([0 1])
ylabel('Normalized dF/F')
subplot(2,2,3)
histogram(max_ind_run)
title('Max resp- Running')
set(gca, 'XTick', 1:3, 'XTickLabel', group_str)
subplot(2,2,4)
histogram(max_ind_norun)
title('Max resp- Stationary')
set(gca, 'XTick', 1:3, 'XTickLabel', group_str)
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_coherenceSummary.pdf']),'-dpdf','-bestfit')

%%
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
    errorbar(1:3, data_dfof_dir_run(iCell,:,1), data_dfof_dir_run(iCell,:,2), '-o')
    hold on
    errorbar(1:3, data_dfof_dir_norun(iCell,:,1), data_dfof_dir_norun(iCell,:,2), '-o')
    set(gca, 'XTick', 1:3, 'XTickLabel', group_str)
    xlim([0 4])
    ylim([0 inf])
    title(num2str(iCell))
    start = start +1;
end
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_cellTuningDir' num2str(n) '.pdf']),'-dpdf','-bestfit')

figure;

for i = 1:3
    subplot(3,1,i)
    ind_cell = find(max_ind_norun==i);
    errorbar(1:3, mean(data_dfof_dir_run_norm(ind_cell,:),1), std(data_dfof_dir_run_norm(ind_cell,:),[],1)./sqrt(length(ind_cell)), '-o')
    hold on
    errorbar(1:3, mean(data_dfof_dir_norun_norm(ind_cell,:),1), std(data_dfof_dir_norun_norm(ind_cell,:),[],1)./sqrt(length(ind_cell)), '-o')
    set(gca, 'XTick', 1:3, 'XTickLabel', group_str)
    xlim([0 4])
    ylim([0 inf])
    title([group_str{i} '- n = ' num2str(length(ind_cell))])
    ylabel('Normalized dF/F')
    legend({'Running','Stationary'},'location','southeast')
end