%[text] #  Load, register, segment and neuropil correct 2P data
%%
%[text] ## Path names
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff';
lg_fn = fullfile(fn_base, 'home\lindsey');
data_fn = fullfile(lg_fn, 'Tutorials');
mworks_fn = fullfile(fn_base, 'Behavior\Data');
fnout = fullfile(fn_base, 'home\David', 'Analysis');
%%
%[text] ## Specific experiment information
date = '221129';
ImgFolder = '003';
time = '1526';
mouse = 'i2080';
frame_rate = 15.5;
run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];
%%
%[text] ## Load 2P data
%[text] Load mworks data (stimulus and behavior information)
fName = fullfile(mworks_fn, ['data-' mouse '-' date '-' time '.mat']);
load(fName);
%[text] Load 2P metadata
CD = fullfile(data_fn, '2p_data', mouse, date, ImgFolder);
cd(CD);
imgMatFile = [ImgFolder '_000_000.mat'];
load(imgMatFile);
%[text] Load 2P images
totframes = input.counterValues{end}(end); %this is from the mworks structure- finds the last value clocked for frame count
fprintf(['Reading ' num2str(totframes) ' frames \r\n']) %[output:6bbd9536]
data = sbxread([ImgFolder '_000_000'],0,totframes);
%[text] Data is nPMT x nYpix x nXpix x nframes. 
fprintf(['Data is ' num2str(size(data))]) %[output:8364a6b0]
%[text] When imaging a single channel, nPMT = 1, so squeeze:
data = squeeze(data);
sz = size(data);
%%
%[text] ## Register 2P data
%[text] Goal here is to remove X-Y movement artifacts
%[text] 1\. Find a stable target
%[text]     a. Plot average of 500 frames throughout stack
nframes = 500; %nframes to average
nskip = 1500; %nframes to skip for each average

nep = floor(size(data,3)./nskip);
[n n2] = subplotn(nep); 
figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:nep; 
    subplot(n,n2,i); 
    imagesc(mean(data(:,:,1+((i-1)*nskip):nframes+((i-1)*nskip)),3)); 
    title([num2str(1+((i-1)*nskip)) '-' num2str(nframes+((i-1)*nskip))]); 
end

%[text]     b. GUI to select target image- choose one that is sharp and close to center of stack
f=gcf;
w = waitforbuttonpress; %click on subplot
if w == 0
    axesClicked = gca;
    allAxes = flipud(findobj(f.Children,'Type','axes'));
    numClicked = find(axesClicked==allAxes);
    close all
end
fprintf(['Selected subplot ' num2str(numClicked)])
%[text]     c. Create target image
data_avg = mean(data(:,:,1+((numClicked-1)*nskip):nframes+((numClicked-1)*nskip)),3); 
%[text] 2\. stackRegister minimizes the difference of each frame from the target
[out, data_reg] = stackRegister(data,data_avg);
%[text] New average image after registration
data_reg_avg = mean(data_reg,3);
%[text] Save registration shifts and target, and mworks data
mkdir(fullfile(fnout, datemouse, datemouserun))
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_reg_shifts.mat']), 'data_reg_avg', 'out', 'data_avg')
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_input.mat']), 'input')
%[text] Test registration
%[text]      a. Make sure first and last images are sharp and similar. Focus on vasculature and nuclei- not cells since F changes
ind = [1 nep];
for i = 1:length(ind) 
    subplot(2,1,i); 
    ix = ind(i);
    imagesc(mean(data_reg(:,:,1+((ix-1)*nskip):nframes+((ix-1)*nskip)),3)); 
    title([num2str(1+((ix-1)*nskip)) '-' num2str(nframes+((ix-1)*nskip))]); 
end
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_FOV_first&last.pdf']), '-dpdf')
%[text]     b. Average of all frames should be sharp
imagesq(data_reg_avg); 
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_FOV_avg.pdf']), '-dpdf')
clear data
%%
%[text] ## Segment 2P data
%[text] Goal here is to create a cell mask to extract fluorescence timecourses
%[text] 1\. Find activated cells
%[text]     a. Reshape data stack to segregate trials 
%[text]     First need to know how many frames per trial and how many trials
nOn = input.nScansOn;
nOff = input.nScansOff;
ntrials = size(input.tGratingDirectionDeg,2); %this is a cell array with one value per trial, so length = ntrials
%[text]         Each trial has nOff frames followed by nOn frames, so can reshape stack so nYpix x nXpix x nFrames/Trial (nOn+nOff) x nTrials
data_tr = reshape(data_reg,[sz(1), sz(2), nOn+nOff, ntrials]);
fprintf(['Size of data_tr is ' num2str(size(data_tr))]) %[output:197654d8]
%[text]     b. Find baseline F from last half of off period- avoids decay of previous on trial
data_f = mean(data_tr(:,:,nOff/2:nOff,:),3); 
%[text]     c. Find dF/F for each trial
data_df = bsxfun(@minus, double(data_tr), data_f); 
data_dfof = bsxfun(@rdivide,data_df, data_f); 
clear data_f data_df data_tr
%[text]     d. Find average dF/F for each stimulus condition (this is for an experiment with changing grating direction)
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
%[text]         Filtering data helps make cells more visible for selection
myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_avg_all = imfilter(data_dfof_avg,myfilter);
data_dfof_max = max(data_dfof_avg_all,[],3); %finds all active cells by taking max projection

save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg_all', 'nStim')
%[text]     2\. Create cell masks from active cells
%[text]         Concatenate images of max and all directions
data_dfof = cat(3,data_dfof_max, data_dfof_avg_all);
%[text]         Set up empty matrices for segmenting cells
mask_data = data_dfof;
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
%[text]         a. Click on bright cells in mask to create mask- blue shading should fill bright region of single cell evenly.  
%[text]             If cell not filled well, press 'z' or shift+click to remove. Try again by 1) clicking on different spot- brighter spot will make fill smaller- if there are two cells also helps to fill brighter cell first; 2) change threshold- higher thresh will make fill smaller; 3) wait for different image where cell might be more distinct from background
%[text]             Press 'return' when done with each image. 
%[text]             Next images will have black spots where cells were previously selected- don't let new masks encroach on old cells or will fuse. imCellBuffer helps with this.
for iStim = 1:size(data_dfof,3)    
    mask_data_temp = mask_data(:,:,iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0; %blacks out old cells
    bwout = imCellEditInteractiveLG(mask_data_temp); %selection GUI
    mask_all = mask_all+bwout; %adds new cells to old cells
    mask_exp = imCellBuffer(mask_all,3)+mask_all; %creates buffer around cells to avoid fusing
    close all
end
%[text]         b. Number independent regions for mask
mask_cell = bwlabel(mask_all); %turns logical into numbered cells
figure;
imagesc(mask_cell)
%[text]         c. Create neuropil masks (these are regions around each cell without overlap from neighboring cells)
nMaskPix = 5; %thickness of neuropil ring in pixels
nBuffPix = 3; %thickness of buffer between cell and ring
mask_np = imCellNeuropil(mask_cell,nBuffPix,nMaskPix);

save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_np')
clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_2 data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 
%%
%[text] ## Neuropil subtraction
%[text] Goal is to remove contamination from out-of-focus fluorescence
%[text] 1. Extract cell timecourses \
data_tc = stackGetTimeCourses(data_reg, mask_cell); %applies mask to stack to get timecourses
%[text]             Timecourses are nFrames x nCells
fprintf(['data_tc is ' num2str(size(data_tc))])  %[output:983e04a3]
nCells = size(data_tc,2);
%[text]             Downsampled timecourses for neuropil subtraction- averageing decreases noise
down = 5; %number of frames to average
data_reg_down = stackGroupProject(data_reg,down); %averages every 5 frames in stack  
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);    
%[text]         2\.  Extract neuropil timecourses (full and downsampled)
sz = size(data_reg);
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells %[output:group:1e95ca64]
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '%s\n'])  %[output:7f517253]
end %[output:group:1e95ca64]
%[text]         3\. Find best neuropil weights by maximizing skew on downsampled subtractions
%[text]             Assumes calcium signals are 1) sparse and 2) positive.  Too little subtraction will decrease sparseness and too much will make signals negative. Thus, the best neuropil subtraction should yield the highest skew.
%[text]         a. Measure skew for all weights:
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
%[text]         b. Find index with highest skew for each cell
[max_skew ind] =  max(x,[],1); 
%[text]         c. convert to weight
np_w = 0.01*ind; 
%[text]         4\. Subtract weighted neuropil response from full timecourses
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
clear data_reg data_reg_down

save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":33.3}
%---
%[output:6bbd9536]
%   data: {"dataType":"text","outputData":{"text":"Reading 14400 frames \n","truncated":false}}
%---
%[output:8364a6b0]
%   data: {"dataType":"text","outputData":{"text":"Data is 1    529    796  14400","truncated":false}}
%---
%[output:197654d8]
%   data: {"dataType":"text","outputData":{"text":"Size of data_tr is 529  796   90  160","truncated":false}}
%---
%[output:983e04a3]
%   data: {"dataType":"text","outputData":{"text":"data_tc is 14400    148","truncated":false}}
%---
%[output:7f517253]
%   data: {"dataType":"text","outputData":{"text":"Cell #1Cell #2Cell #3Cell #4Cell #5Cell #6Cell #7Cell #8Cell #9Cell #10Cell #11Cell #12Cell #13Cell #14Cell #15Cell #16Cell #17Cell #18Cell #19Cell #20Cell #21Cell #22Cell #23Cell #24Cell #25Cell #26Cell #27Cell #28Cell #29Cell #30Cell #31Cell #32Cell #33Cell #34Cell #35Cell #36Cell #37Cell #38Cell #39Cell #40Cell #41Cell #42Cell #43Cell #44Cell #45Cell #46Cell #47Cell #48Cell #49Cell #50Cell #51Cell #52Cell #53Cell #54Cell #55Cell #56Cell #57Cell #58Cell #59Cell #60Cell #61Cell #62Cell #63Cell #64Cell #65Cell #66Cell #67Cell #68Cell #69Cell #70Cell #71Cell #72Cell #73Cell #74Cell #75Cell #76Cell #77Cell #78Cell #79Cell #80Cell #81Cell #82Cell #83Cell #84Cell #85Cell #86Cell #87Cell #88Cell #89Cell #90Cell #91Cell #92Cell #93Cell #94Cell #95Cell #96Cell #97Cell #98Cell #99Cell #100Cell #101Cell #102Cell #103Cell #104Cell #105Cell #106Cell #107Cell #108Cell #109Cell #110Cell #111Cell #112Cell #113Cell #114Cell #115Cell #116Cell #117Cell #118Cell #119Cell #120Cell #121Cell #122Cell #123Cell #124Cell #125Cell #126Cell #127Cell #128Cell #129Cell #130Cell #131Cell #132Cell #133Cell #134Cell #135Cell #136Cell #137Cell #138Cell #139Cell #140Cell #141Cell #142Cell #143Cell #144Cell #145Cell #146Cell #147Cell #148","truncated":false}}
%---
