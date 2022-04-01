%% Load, register, segment and neuropil correct 2P data
close all
clear all global
close all

% Pull in experiment data
ds = 'AllExperimentData_Multiday_TD'; %dataset info
dataStructLabels = {'runs'};
eval(ds)

day_id = 18;

%Path names
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
td_fn = fullfile(fn_base, 'home\Tierney');
data_fn = fullfile(td_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data');
fnout = fullfile(fn_base, 'home\Tierney\Analysis\2P');

%Specific experiment information
mouse = expt(day_id).mouse;
date = expt(day_id).date;
ImgFolder = expt(day_id).stimruns;
time = expt(day_id).time;
eye_str = expt(day_id).eye_str;

frame_rate = 15;
contra = strcmp(eye_str,'Contra'); %blocks for eye stimulation: 1 is when contra is open; 0 is contra closed
nrun = size(ImgFolder,1);
run_str = catRunName(ImgFolder, nrun);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];

%% Load 2P data

data = [];
clear temp
offset_frames = 0;
offset_time = 0;

for irun = 1:nrun
    %Load mworks data- this has information about experiment (e.g. the visual stimuli presented and synchronization with the microscope)
    fName = fullfile(mworks_fn, ['data-' mouse '-' date '-' time{irun} '.mat']);
    load(fName);
    ntrials = size(input.tGratingDirectionDeg,2);
    totframes = input.counterValues{end}(end); %this is from the mworks structure- finds the last value clocked for frame count
    input.contraTrialNumber = mat2cell(contra(irun).*ones(1,ntrials),1,ones(1,ntrials)); %create a field for tracking which eye stimulated
    temp(irun) = input;
    
    if irun>1
        for itrial = 1:ntrials
            temp(irun).counterValues{itrial} = bsxfun(@plus,temp(irun).counterValues{itrial},offset_frames);
            temp(irun).counterTimesUs{itrial} = bsxfun(@plus,temp(irun).counterTimesUs{itrial},offset_time-temp(irun).counterValues{1}(1)+1000);
            temp(irun).wheelSpeedTimesUs{itrial} = bsxfun(@plus,temp(irun).wheelSpeedTimesUs{itrial},offset_time-temp(irun).counterValues{1}(1)+1000);
        end
    end
    offset_frames = offset_frames+totframes;
    offset_time = offset_time+temp(irun).counterTimesUs{end}(end);
    
    %Load 2P metadata- this has information about the information about the imaging session (e.g. frame count, zoom)
    CD = fullfile(data_fn, mouse, date, ImgFolder(irun,:));
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    %Load 2P images
    
    fprintf(['Reading ' num2str(totframes) ' frames \r\n'])
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,totframes); %loads the .sbx files with imaging data (path, nframes to skip, nframes to load)
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
end

input = concatenateDataBlocks(temp);
clear data_temp
clear temp
%% Register 2P data
 
%Goal here is to remove X-Y movement artifacts
%1. Find a stable target
%    a. Plot average of 500 frames throughout stack

nep = 10; %numer of epochs
nframes = 500; %nframes to average for target
nskip = floor(size(data,3)./nep); %nframes to skip for each average

if exist(fullfile(fnout, datemouse, datemouserun)) 
    load(fullfile(fnout, datemouse, datemouserun, [datemouserun '_reg_shifts.mat']))
    [outs, data_reg]=stackRegister_MA(data(:,:,:),[],[],out);
    save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_input.mat']), 'input')
else
    [n, n2] = subplotn(nep); 
    figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:nep 
        subplot(n,n2,i); 
        imagesc(mean(data(:,:,1+((i-1)*nskip):nframes+((i-1)*nskip)),3)); 
        title([num2str(1+((i-1)*nskip)) '-' num2str(nframes+((i-1)*nskip))]); 
    end

    % b. GUI to select target image- choose one that is sharp and close to center of stack
    f=gcf;
    w = waitforbuttonpress; %click on subplot
    if  w == 0
        axesClicked = gca;
        allAxes = findobj(f.Children,'Type','axes');
        numClicked = find(axesClicked==allAxes);
        close all
    end
    fprintf(['Selected subplot ' num2str(numClicked) '\n'])

    % c. Create target image
    data_avg = mean(data(:,:,1+((numClicked-1)*nskip):nframes+((numClicked-1)*nskip)),3); %average 500 frames to make target
    %2. stackRegister minimizes the difference of each frame from the target
    [out, data_reg] = stackRegister(data,data_avg);
    %New average image after registration
    data_reg_avg = mean(data_reg,3);
    %Save registration shifts and target, and mworks data
    mkdir(fullfile(fnout, datemouse, datemouserun))
    save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_reg_shifts.mat']), 'data_reg_avg', 'out', 'data_avg')
    save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_input.mat']), 'input')
end

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
figure
imagesq(data_reg_avg); 
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_FOV_avg.pdf']), '-dpdf')

%% Segment 2P data
clear data
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
data_tr = data_tr(:,:,nOff/2:end,:);
fprintf(['Size of data_tr is ' num2str(size(data_tr))])
%    b. Find baseline F from last half of off period- avoids decay of previous on trial
data_f = single(mean(data_tr(:,:,1:nOff/2,:),3)); 
%    c. Find dF/F for each trial
data_df = single(data_tr)-data_f;
%data_df = bsxfun(@minus, double(data_tr), data_f);
clear data_tr
data_dfof = bsxfun(@rdivide,data_df, data_f); 
clear data_f data_df
%    d. Find average dF/F for each stimulus condition (this is for an experiment with changing grating direction)
tDir = celleqel2mat_padded(input.tGratingDirectionDeg); %transforms cell array into matrix (1 x ntrials)
Dirs = unique(tDir);
nDirs = length(Dirs);

tContra = celleqel2mat_padded(input.contraTrialNumber); %transforms cell array into matrix (1 x ntrials)
Eyes = unique(tContra);
nEye = length(Eyes);

data_dfof_avg = zeros(sz(1),sz(2),nDirs,nEye); %create empty matrix with FOV for each direction: nYpix x nXPix x nDir
data_dfof_avg_all = data_dfof_avg;
myfilter = fspecial('gaussian',[20 20], 0.5);
for ieye = 1:nEye
    ind_eye = find(tContra == Eyes(ieye)); %find all trials with each eye
    figure; movegui('center')
    [n, n2] = subplotn(nDirs); %function to optimize subplot number/dimensions
    for idir = 1:nDirs
        ind_dir = find(tDir == Dirs(idir)); %find all trials with each direction
        ind = intersect(ind_eye,ind_dir);
        data_dfof_avg(:,:,idir,ieye) = mean(mean(data_dfof(:,:,nOff/2+1:end,ind),3),4); %average all On frames and all trials
        subplot(n,n2,idir)
        imagesc(data_dfof_avg(:,:,idir,ieye))
    end
    data_dfof_avg_all(:,:,:,ieye) = imfilter(data_dfof_avg(:,:,:,ieye),myfilter);
end
clear data_dfof
nStim = nDirs;
%        Filtering data helps make cells more visible for selection

data_dfof_max = max(max(data_dfof_avg_all,[],4),[],3); %finds all active cells by taking max projection
figure; movegui('center'); imagesc(data_dfof_max);
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_stimData_Dir.mat']), 'nOn','nOff','ntrials','tDir','Dirs','nDirs','tContra','Eyes','nEye');
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_stimActFOV_Dir.mat']), 'data_dfof_max', 'data_dfof_avg_all', 'nStim')
%% 2. Create cell masks from active cells

% Get pixel correlation image, another method to identify cells
data_down = stackGroupProject(data_reg,100);
corrImg = getPixelCorrelationImage(data_down);
figure; imagesc(corrImg); movegui('center');title('pixel correlation');
clear data_down

%        Concatenate images of max and all directions
data_dfof = cat(3,data_dfof_max, max(data_dfof_avg_all,[],4),corrImg);

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

save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_mask_cell.mat']), 'corrImg', 'data_dfof_max', 'mask_cell', 'mask_np')
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

save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_TCs_Dir.mat']), 'data_tc', 'np_tc', 'npSub_tc', 'nCells')

%% looking at time courses

data_tc_trial = reshape(npSub_tc, [nOn+nOff,ntrials,nCells]);
data_f_trial = mean(data_tc_trial(nOff/2:nOff,:,:),1);
data_dfof_trial = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial, data_f_trial), data_f_trial);

%looking at data with np subtracted
tc_cell_avrg = mean(data_dfof_trial,3);%average pver cells, one row per trial
tc_trial_avrg = squeeze(mean(data_dfof_trial,2));%average over trials, one row per cell
tc_cell_trial_avrg = mean(tc_cell_avrg,2);%average over trials and cells

figure;
plot(tc_trial_avrg, 'LineWidth',.005,'color',[.25 .25 .25]);
hold on;
plot(tc_cell_trial_avrg, 'LineWidth',2, 'color','k');
hold on;
title('Timecourses with np subtracted');
hold off

%% make tuning curves
base_win = nOff/2:nOff;
resp_win = nOff+5:nOff+nOn;
data_trial = reshape(npSub_tc, [nOn+nOff ntrials nCells]);
data_f = mean(data_trial(base_win,:,:),1);
data_dfof = bsxfun(@rdivide,bsxfun(@minus,data_trial,data_f),data_f);

resp_cell_dir = cell(nEye,nDirs);
base_cell_dir = cell(nEye,nDirs);
data_dfof_dir = zeros(nCells,nDirs,nEye,2);
h_dir = zeros(nDirs,nCells,nEye);
p_dir = zeros(nDirs,nCells,nEye);
for iEye = 1:nEye
    ind_eye = find(tContra == Eyes(iEye));
    for iDir = 1:nDirs
        ind_dir = find(tDir==Dirs(iDir));
        ind = intersect(ind_eye,ind_dir);
        resp_cell_dir{iEye,iDir} = squeeze(mean(data_dfof(resp_win,ind,:),1));
        base_cell_dir{iEye,iDir} = squeeze(mean(data_dfof(base_win,ind,:),1));
        [h_dir(iDir,:,iEye), p_dir(iDir,:,iEye)] = ttest(resp_cell_dir{iEye,iDir},base_cell_dir{iEye,iDir},'tail','right','alpha',0.05./((nDirs*nEye)-1));
        data_dfof_dir(:,iDir,iEye,1) = squeeze(mean(mean(data_dfof(resp_win,ind,:),1),2));
        data_dfof_dir(:,iDir,iEye,2) = squeeze(std(mean(data_dfof(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
    end
end

h_all_dir = squeeze(sum(h_dir,1));
resp_ind = find(sum(h_all_dir,2));
ipsi_resp_ind = find(h_all_dir(:,1));
contra_resp_ind = find(h_all_dir(:,2)); 

save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_respData_Dir.mat']), 'h_dir', 'resp_ind', 'ipsi_resp_ind', 'contra_resp_ind', 'resp_cell_dir', 'base_cell_dir', 'data_dfof_dir','base_win','resp_win','data_dfof')

%plot histogram of number of significant directions for contra and ipsi
figure;
[n n2] = subplotn(nEye);
for iEye = 1:nEye
subplot(n,n2,iEye)
histogram(h_all_dir(:,iEye),[0:1:nDirs])
title(eye_str{find(contra==Eyes(iEye))})
xlabel('Significant directions')
ylabel('Cells')
ylim([0 50]);
end
sgtitle([mouse ' ' date])
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_sigDirsHist.pdf']),'-dpdf','-bestfit')

%plot all cells for all direction to both contra and ipsi
start=1;
n = 1;
figure;
movegui('center')
for iCell = 1:nCells
    if start>25
        sgtitle([mouse ' ' date])
        print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_cellTuningDir' num2str(n) '.pdf']),'-dpdf','-bestfit')
        figure;movegui('center');
        start = 1;
        n = n+1;
    end
    subplot(5,5,start)
    for iEye = 1:2
        errorbar(Dirs, data_dfof_dir(iCell,:,iEye,1), data_dfof_dir(iCell,:,iEye,2), '-o')
        hold on
    end
    title(['R = ' num2str(h_all_dir(iCell,:))])
    start = start +1;
end
sgtitle([mouse ' ' date])
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_cellTuningDir' num2str(n) '.pdf']),'-dpdf','-bestfit')

[contra_resp max_ind_contra] = max(data_dfof_dir(:,:,find(Eyes),1),[],2);
max_dir_contra = Dirs(max_ind_contra);
[ipsi_resp max_ind_ipsi] = max(data_dfof_dir(:,:,find(~Eyes),1),[],2); 
max_dir_ipsi = Dirs(max_ind_ipsi);

real_contra_resp = contra_resp;
real_contra_ind = find(real_contra_resp < 0);
real_contra_resp(real_contra_ind) = 0;

real_ipsi_resp = ipsi_resp;
real_ipsi_ind = find(real_ipsi_resp < 0);
real_ipsi_resp(real_ipsi_ind) = 0;

ODI = (real_contra_resp-real_ipsi_resp)./(real_contra_resp+real_ipsi_resp);

%scatter of max response to contra and ipsi
figure;
subplot(2,3,1)
contra_resp_any = max(data_dfof_dir(resp_ind,:,find(Eyes),1),[],2);
ipsi_resp_any = max(data_dfof_dir(resp_ind,:,find(~Eyes),1),[],2);
scatter(contra_resp_any,ipsi_resp_any)
axis square
refline(1)
xlabel('Max (Contra)')
ylabel('Max (Ipsi)')
xlim([0 1])
ylim([0 1])
title('Responsive (any)')

subplot(2,3,4)
ODI_any = (contra_resp_any-ipsi_resp_any)./(contra_resp_any+ipsi_resp_any);
histogram(ODI_any,[-1:0.1:1])
[temp_counts_any,~] = histcounts(ODI_any,[-1:0.1:1])
xlabel('ODI')
ylabel('Cells')
axis square
title(['Sig cells = ' num2str(length(resp_ind))])

subplot(2,3,2)
contra_resp_contraonly = max(data_dfof_dir(contra_resp_ind,:,find(Eyes),1),[],2); 
ipsi_resp_contraonly = max(data_dfof_dir(contra_resp_ind,:,find(~Eyes),1),[],2);
scatter(contra_resp_contraonly,ipsi_resp_contraonly)
axis square
refline(1)
xlabel('Max (Contra)')
ylabel('Max (Ipsi)')
xlim([0 1])
ylim([0 1])
title('Contra responsive')

subplot(2,3,5)
ODI_contraonly = (contra_resp_contraonly-ipsi_resp_contraonly)./(contra_resp_contraonly+ipsi_resp_contraonly);
histogram(ODI_contraonly,[-1:0.1:1])
[temp_counts_contraonly,~] = histcounts(ODI_contraonly,[-1:0.1:1])
xlabel('ODI')
ylabel('Cells')
axis square
title(['Sig cells = ' num2str(length(contra_resp_ind))])

subplot(2,3,3)
contra_resp_ipsionly = max(data_dfof_dir(ipsi_resp_ind,:,find(Eyes),1),[],2);
ipsi_resp_ipsionly = max(data_dfof_dir(ipsi_resp_ind,:,find(~Eyes),1),[],2);
scatter(contra_resp_ipsionly,ipsi_resp_ipsionly)
axis square
refline(1)
xlabel('Max (Contra)')
ylabel('Max (Ipsi)')
xlim([0 1])
ylim([0 1])
title('Ipsi responsive')

subplot(2,3,6)
ODI_ipsionly = (contra_resp_ipsionly-ipsi_resp_ipsionly)./(contra_resp_ipsionly+ipsi_resp_ipsionly);
histogram(ODI_ipsionly,[-1:0.1:1])
[temp_counts_ipsionly,~] = histcounts(ODI_ipsionly,[-1:0.1:1]);
xlabel('ODI')
ylabel('Cells')
axis square
title(['Sig cells = ' num2str(length(ipsi_resp_ind))])

max_histcounts = max([temp_counts_any temp_counts_contraonly temp_counts_ipsionly]);
subplot(2,3,4)
axis([-1 1 0 max_histcounts+1])
subplot(2,3,5)
axis([-1 1 0 max_histcounts+1])
subplot(2,3,6)
axis([-1 1 0 max_histcounts+1])

print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_OD_Dir.pdf']),'-dpdf','-bestfit')

save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_ODI_Dir.mat']), 'ODI', 'real_contra_resp', 'real_ipsi_resp', 'max_dir_contra', 'max_dir_ipsi')

%% von mises 

    b_ori = zeros(nEye,nCells);
    k1_ori = zeros(nEye,nCells);
    R1_ori = zeros(nEye,nCells);
    u1_ori = zeros(nEye,nCells);
    R_square_ori = zeros(nEye,nCells);
    sse_ori = zeros(nEye,nCells);
    stim_DSI = zeros(nEye,nCells);
    stim_OSI = zeros(nEye,nCells);
    theta_hires= deg2rad(0:180);
    y_fit = zeros(length(theta_hires),nEye,nCells);
    
%     Oris = Dirs(1:nDirs/2);
%     nOris = length(Oris);
    for iEye = 1:nEye
        data_dfof_ori(:,:,iEye)= mean(reshape(data_dfof_dir(:,:,iEye,1),[nCells nOris 2]),3);
        for iCell = 1:nCells
            data = [data_dfof_ori(iCell,:,iEye) data_dfof_ori(iCell,1,iEye)];
            theta = [deg2rad(Oris) pi];
            [b_ori(iEye,iCell),k1_ori(iEye,iCell),R1_ori(iEye,iCell),u1_ori(iEye,iCell),sse_ori(iEye,iCell),R_square_ori(iEye,iCell)] ...
                = miaovonmisesfit_ori(theta,data);
            [max_val max_ind] = max(data_dfof_ori(iCell,:,iEye),[],2);
            null_ind = max_ind+(nOris./2);
            null_ind(find(null_ind>nOris)) = null_ind(find(null_ind>nOris))-nOris;
            min_val = data_dfof_ori(iCell,null_ind,iEye);
            if min_val<0
                min_val = 0;
            end
            stim_OSI(1,iCell) = (max_val-min_val)./(max_val+min_val);
            [max_val max_ind] = max(data_dfof_dir(iCell,:,iEye,1),[],2); % SUPPOSED TO BE DIR OR ORI??????
            null_ind = max_ind+(nDirs./2);
            null_ind(find(null_ind>nDirs)) = null_ind(find(null_ind>nDirs))-nDirs;
            min_val = data_dfof_dir(iCell,null_ind,iEye,1);
            if min_val<0
                min_val = 0;
            end
            stim_DSI(1,iCell) = (max_val-min_val)./(max_val+min_val);
            y_fit(:,iEye,iCell) = b_ori(iEye,iCell) + R1_ori(iEye,iCell) .* exp(k1_ori(iEye,iCell).*(cos(2.*(theta_hires-u1_ori(iEye,iCell)))-1));
        end
    end
    
    [yfit_max, yfit_max_ind] = max(y_fit(:,:,:),[],1);
    prefOri_yfit = squeeze(theta_hires(yfit_max_ind));
    
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_oriResp_Dir.mat']), 'data_dfof_dir', 'data_dfof_ori','base_win','resp_win','h_dir', 'b_ori', 'k1_ori', 'R1_ori', 'u1_ori', 'R_square_ori', 'sse_ori','stim_DSI','stim_OSI', 'yfit_max', 'yfit_max_ind', 'prefOri_yfit') 