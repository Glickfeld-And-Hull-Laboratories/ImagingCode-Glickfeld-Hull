%% This code allows us to identify visually responsive cells and create a 
% timecourse of how those cells change fluorescence across all frames

clear all
clear all global
%%  Load, register, segment, and neuropil correct 2P data
% Path names

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff'; % base file name
td_fn = fullfile(fn_base, 'home\Tierney'); % Tierney's personal file folder
data_fn = fullfile(td_fn, '\2p_Imaging\2p_Tutorial'); % This is where the data is 
mworks_fn = fullfile(fn_base, 'Behavior\Data'); % This is where your mworks data is (from 2p stimulus)
fnout = fullfile(td_fn, 'Analysis\2p_Analysis\2p_Tutorial_Analysis'); % This is where we designate our ouput folder
%% Specific experiment information

date = '211020'; % Specifies which imaging data folder you want to grab 
ImgFolder = '002'; % Imaging that you're interested in 
time = '0917'; % Time that the imaging was started. Add 0 before three digit numbers
mouse = 'tj_092021'; % Specifies mouse
frame_rate = 15.5; % Specifies frame rate 
run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' mouse]; % Creates a variable of date and mouse together
datemouserun = [date '_' mouse '_' run_str]; % Creates a variable of date, mouse, and run together
%% Load 2P data
% Load mworks data (stimulus and behavior information)

fName = fullfile(mworks_fn, ['data-i' mouse '-' date '-' time '.mat']); % This just adds data-i before name so my file matches the mworks output; also loads input variable 
load(fName);
%% Load 2P metadata

CD = fullfile(data_fn, datemouse, ImgFolder); % Sets current directory in folder path
cd(CD); % Changes current directory
imgMatFile = [ImgFolder '_000_000.mat']; % Finds the matlab file for whichever run you want to use
load(imgMatFile);

%% Load 2P images

totframes = input.counterValues{end}(end); %this is from the mworks structure: finds the last value clocked for frame count
fprintf(['Reading ' num2str(totframes) ' frames \r\n']) % This will spit out "Reading xxx frames" in the command window
data = sbxread([ImgFolder '_000_000'],0,totframes); % This reads the scanbox data (input, how many frames we want to skip, total frames we want)
%% Data is nPMT x nYpix x nXpix x nframes. 

fprintf(['Data is ' num2str(size(data))]) % Tells me dimensions of my dataset
%% When imaging a single channel, nPMT = 1, so squeeze:

data = squeeze(data); % This removes PMT from the dimensions set above because we are only imaging one channel (PMT)
%% Register 2P data. Register can correct for motion in x-y plane, but not in z
% plane.
% Goal here is to remove X-Y movement artifacts. We will pick one average
% of 500 frames to match the rest of the frames to.
% 
% 1. Find a stable target
% 
% a. Plot average of 500 frames throughout stack

nframes = 500; %nframes to average
nskip = 1500; %nframes to skip for each average; this may change depending on how many frames we collect total

nep = floor(size(data,3)./nskip); % This determines how many subplots will be in the figure to choose from
[n n2] = subplotn(nep); 
figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:nep; % Makes subplots of 500 frame average for subplots 1 through nep (9)
    subplot(n,n2,i); 
    imagesc(mean(data(:,:,1+((i-1)*nskip):nframes+((i-1)*nskip)),3)); 
    title([num2str(1+((i-1)*nskip)) '-' num2str(nframes+((i-1)*nskip))]); 
end

%% 
% b. GUI to select target image: choose one that is sharp and close to center 
% of stack.

f=gcf;
w = waitforbuttonpress; %click on subplot
if w == 0
    axesClicked = gca;
    allAxes = flipud(findobj(f.Children,'Type','axes'));
    numClicked = find(axesClicked==allAxes);
    close all
end
fprintf(['Selected subplot ' num2str(numClicked)])
%% 
% c. Create target image
% This designates the frame that I picked in the last section to be used to make adjustments
% to all other frames

data_avg = mean(data(:,:,1+((numClicked-1)*nskip):nframes+((numClicked-1)*nskip)),3); 
%% 
% 2. stackRegister minimizes the difference of each frame from the target

% This section takes a long time to run. If I've already Registered once
% and have created the folder, the "if" statement will pull in that folder
% so I don't have to redo this whole step

if exist(fullfile(fnout, [date '_' mouse '_' run_str]))
    load(fullfile(fnout, [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(fnout, [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(data(:,:,:),[],[],out);
else
    [out, data_reg] = stackRegister(data,data_avg);
    % two outputs from stackRegister: out is a matrix that tells you how
    % much you need to shift to match your target image and
    % data_reg is the registered stack; version of data that is now fixed
    data_reg_avg = mean(data_reg,3);
    reg = data_reg_avg;
    mkdir(fullfile(fnout,[date '_' mouse '_' run_str]))
    save(fullfile(fnout, [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'data_reg_avg', 'out', 'data_avg')
    save(fullfile(fnout, [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
% %% 
% % New average image after registration
% 
% data_reg_avg = mean(data_reg,3);
% %% 
% % Save registration shifts and target, and mworks data
% 
% mkdir(fullfile(fnout, datemouse, datemouserun))
% save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_reg_shifts.mat']), 'data_reg_avg', 'out', 'data_avg')
% save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_input.mat']), 'input')
%% 
% Test registration
% 
% a. Make sure first and last images are sharp and similar. Focus on vasculature 
% and nuclei- not cells since F changes. If images are not clear/similar,
% try redoing Register using a different reference image. If this doesn't
% work, data might not be good for analysis.

ind = [1 nep];
for i = 1:length(ind) 
    subplot(2,1,i); 
    ix = ind(i);
    imagesc(mean(data_reg(:,:,1+((ix-1)*nskip):nframes+((ix-1)*nskip)),3)); 
    title([num2str(1+((ix-1)*nskip)) '-' num2str(nframes+((ix-1)*nskip))]); 
end
print(fullfile(fnout, [date '_' mouse '_' run_str], [datemouserun '_FOV_first&last.pdf']), '-dpdf')

%% 
% b. Average of all frames should be sharp

figure;
imagesq(data_reg_avg); 
print(fullfile(fnout, [date '_' mouse '_' run_str], [datemouserun '_FOV_avg.pdf']), '-dpdf')
clear data
%% Segment 2P data
% Goal here is to create a cell mask to extract fluorescence timecourses
% 
% 1. Find activated cells
% 
% a. Reshape data stack to segregate trials 
% 
% First need to know how many frames per trial and how many trials

nOn = input.nScansOn; % number of On frames per trial
nOff = input.nScansOff; % number of Off frames per trial
ntrials = size(input.tGratingDirectionDeg,2); %this is a cell array with one value per trial, so length = ntrials
sz = size(data_reg);
%% 
% Each trial has nOff frames followed by nOn frames, so can reshape stack so 
% nYpix x nXpix x nFrames/Trial (nOn+nOff) x nTrials

data_tr = reshape(data_reg,[sz(1), sz(2), nOn+nOff, ntrials]);
fprintf(['Size of data_tr is ' num2str(size(data_tr))])
%% 
% b. Find baseline F from last half of off period- avoids decay of previous 
% on trial

data_f = mean(data_tr(:,:,nOff/2:nOff,:),3); 
%% 
% c. Find dF/F for each trial

data_df = bsxfun(@minus, double(data_tr), data_f); 
data_dfof = bsxfun(@rdivide,data_df, data_f); 
clear data_f data_df data_tr
%% 
% d. Find average dF/F for each stimulus condition (this is for an experiment 
% with changing grating direction)

Dir = celleqel2mat_padded(input.tGratingDirectionDeg); %transforms cell array into matrix (1 x ntrials)
Dirs = unique(Dir);
nDirs = length(Dirs);
data_dfof_avg = zeros(sz(1),sz(2),nDirs); %create empty matrix with FOV for each direction: nYpix x nXPix x nDir

nStim = nDirs;
[n n2] = subplotn(nDirs); %function to optimize subplot number/dimensions
for idir = 1:nDirs
    ind = find(Dir == Dirs(idir)); %find all trials with each direction
    data_dfof_avg(:,:,idir) = mean(mean(data_dfof(:,:,nOff+1:nOn+nOff,ind),3),4); %average all On frames and all trials
    subplot(n,n2,idir)
    imagesc(data_dfof_avg(:,:,idir))
end
clear data_dfof
%% 
% Filtering data helps make cells more visible for selection

myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_avg_all = imfilter(data_dfof_avg,myfilter);
data_dfof_max = max(data_dfof_avg_all,[],3); %finds all active cells by taking max projection

save(fullfile(fnout, [date '_' mouse '_' run_str], [datemouserun '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg_all', 'nStim')
%% 
% 2. Create cell masks from active cells
% 
% Concatenate images of max and all directions

data_dfof = cat(3,data_dfof_max, data_dfof_avg_all);
%% 
% Set up empty matrices for segmenting cells

% This section is creating empty matrices for the masks with specific
% dimensions (dimensions are the number of x and y pixels).
% Mask_all is the cell, mask_exp (expanded) is the cell plus a buffer
% region around the cell. The buffer region prevents two cells from being
% fused together if they are in close proximity. Can look at teh difference
% using the following command:
% figure(3); imagesc(mask_all)
% figure(4); imagesc(mask_exp)

mask_data = data_dfof;
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
%% 
% a. Click on bright cells in mask to create mask- blue shading should fill 
% bright region of single cell evenly.  
% 
% If cell not filled well, press 'z' or shift+click to remove. Try again by 
% 1) clicking on different spot- brighter spot will make fill smaller- if there 
% are two cells also helps to fill brighter cell first; 2) change threshold- higher 
% thresh will make fill smaller; 3) wait for different image where cell might 
% be more distinct from background
% 
% Press 'return' when done with each image. 
% 
% Next images will have black spots where cells were previously selected- don't 
% let new masks encroach on old cells or will fuse. imCellBuffer helps with this.

for iStim = 1:size(data_dfof,3)    
    mask_data_temp = mask_data(:,:,iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0; %blacks out old cells
    bwout = imCellEditInteractiveLG(mask_data_temp); %selection GUI
    mask_all = mask_all+bwout; %adds new cells to old cells
    mask_exp = imCellBuffer(mask_all,3)+mask_all; %creates buffer around cells to avoid fusing
    close all
end
%% 
% b. Number independent regions for mask

mask_cell = bwlabel(mask_all); %turns logical into numbered cells
figure;
imagesc(mask_cell)
%% 
% c. Create neuropil masks (these are regions around each cell without overlap 
% from neighboring cells)

nMaskPix = 5; %thickness of neuropil ring in pixels
nBuffPix = 3; %thickness of buffer between cell and ring
mask_np = imCellNeuropil(mask_cell,nBuffPix,nMaskPix);

save(fullfile(fnout, [date '_' mouse '_' run_str], [datemouserun '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_np')
clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_2 data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 
%% Neuropil subtraction
% Goal is to remove contamination from out-of-focus fluorescence
%% 
% Extract cell timecourses

data_tc = stackGetTimeCourses(data_reg, mask_cell); %applies mask to stack to get timecourses
%% 
% Timecourses are nFrames x nCells

fprintf(['data_tc is ' num2str(size(data_tc))]) 
nCells = size(data_tc,2);
%% 
% Downsampled timecourses for neuropil subtraction- averageing decreases noise

down = 5; %number of frames to average
data_reg_down = stackGroupProject(data_reg,down); %averages every 5 frames in stack  
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);    
%% 
% 2.  Extract neuropil timecourses (full and downsampled)

sz = size(data_reg);
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '%s\n']) 
end
%% 
% 3. Find best neuropil weights by maximizing skew on downsampled subtractions
% 
% Assumes calcium signals are 1) sparse and 2) positive.  Too little subtraction 
% will decrease sparseness and too much will make signals negative. Thus, the 
% best neuropil subtraction should yield the highest skew.

% a. Measure skew for all weights:

ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
%% 
% b. Find index with highest skew for each cell

[max_skew ind] =  max(x,[],1); 
%% 
% c. convert to weight

np_w = 0.01*ind; 
%% 
% 4. Subtract weighted neuropil response from full timecourses

npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
clear data_reg data_reg_down

save(fullfile(fnout, [date '_' mouse '_' run_str], [datemouserun '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

%% OriAnalysis

%% 1. Transform the time courses into a trial structure and calculate dF/F. 

% dF/F = (F-F0)/F0). F0 is the F in the baseline condition. 

% Defines matrix for fluorescent response.
data_f = reshape(npSub_tc, [nOff + nOn, ntrials, nCells]);

% Plot the average timecourse of fluorescence across all trials and cells,
% This allows you to tailor your response window according to the delay in visual response.
data_f_tc = mean(mean(data_f, 2), 3);
plot(data_f_tc);

% Define windows that reference the frames on each trial for baseline and response.. 
base_win = [nOff/2:nOff];
resp_win = [nOff+1:nOff+nOn];

% Create a matrix that is in dF/F instead of F. 
data_f0 = mean(data_f(base_win,:,:),1); 
data_df = bsxfun(@minus, data_f, data_f0); 
data_dfof = bsxfun(@rdivide, data_df, data_f0); 

%% 2. Identify cells that are significantly responsive to at least one orientation. 

%Translate directions to orientations
dir_180 = find(Dir > 179);
norm_dir = Dir(dir_180)-180;
Ori = Dir;
Ori(dir_180) = norm_dir;
Oris = unique(Ori);
nOris = length(Oris);

% Create variables for output of ttest.
h_ori = zeros(nOris, nCells);
p_ori = zeros(nOris, nCells);
resp_avg = zeros(nCells, nOris);
resp_sem = zeros(nCells, nOris);

% Use a ttest to compare baseline and stimulus response across trials.
for iOri = 1:nOris
    ind = find(Ori==Oris(iOri)); %find all trials with each orientation
    base_ori = squeeze(mean(data_dfof(base_win,ind,:),1)); 
    resp_ori = squeeze(mean(data_dfof(resp_win,ind,:),1));
    [h_ori(iOri,:), p_ori(iOri,:)] = ttest(resp_ori,base_ori,'tail','right','alpha',0.05./(nOris-1));
    resp_avg(:,iOri) = squeeze(mean(mean(data_dfof(resp_win,ind,:),1),2));
    resp_sem(:,iOri) = squeeze(std(mean(data_dfof(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
end

% How many cells are significant for at least 1 orientation?
h_ori_sum = sum(h_ori, 1);
h_ori_sig = find(h_ori_sum > 0); %or find(h_ori_sum) because default is > 0
nCells_sig = length(h_ori_sig)

%% 3. Plot average tuning curve for each cell 

nep = nCells; % This determines how many subplots will be in the figure to choose from
[n n2] = subplotn(nep);
for i = 1:nep;
    subplot(n,n2,i); 
    errorbar(Oris, resp_avg(i,:), resp_sem(i,:), '-o')
    % PLOT SIG ORI in diff color here; will use indexing
    if i == h_ori_sum
        hold on
        errorbar(Oris, resp_avg(i,:), resp_sem(i,:), '-o');
    end
    xlabel('Orientation');
    ylabel('Average dF/F')
    title(['Sig ori = ' num2str(h_ori_sum(:,i))])
end

nep = nCells; % This determines how many subplots will be in the figure to choose from
[n n2] = subplotn(nep);
for i = 1:nep;
    subplot(n,n2,i); 
    errorbar(Oris, resp_avg(i,:), resp_sem(i,:), '-o')
    % PLOT SIG ORI in diff color here; will use indexing
    hold on
    errorbar(Oris, resp_avg(,:), resp_sem(i,:); '-o');
    xlabel('Orientation');
    ylabel('Average dF/F')
    title(['Sig ori = ' num2str(h_ori_sum(:,i))])
end

%% 4. Fit data with a von Mises function 

theta = deg2rad(Oris)

% PREALLOCATE MATRICES FOR OUTPUT OF VON MISES
b = zeros(nCells,1) % baseline response
k1 = zeros(nCells,1) % sharpness of orientation tuning
R = zeros(nCells,1) % peak response
u1 = zeros(nCells,1) % preferred orientation
sse = zeros(nCells,1) % sum of squared error
R_square = zeros(nCells,1)

for i = 1:nep;
    for iOri = 1:nOris;
        ind = find(Ori==Oris(iOri)); %find all trials with each orientation
        base_avg(i,iOri) = squeeze(mean(mean(data_dfof(base_win,ind,i),1),2)); 
        resp_avg(i,iOri) = squeeze(mean(mean(data_dfof(resp_win,ind,i),1),2));
        resp_norm = resp_avg(i,iOri) - base_avg(i,iOri);
    end
    [b(i,:), k1(i,:), R(i,:), u1(i,:), sse(i,:), R_square(i,:)] = miaovonmisesfit_ori(theta, resp_norm);
end