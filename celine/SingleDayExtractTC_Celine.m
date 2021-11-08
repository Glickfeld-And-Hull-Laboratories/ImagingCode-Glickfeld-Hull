% to analyze individual day's data
clear all; clear global; close all

%identifying animal and run
date = '211004';
imgFolder = '003';
time = '1457';
mouse = 'i2015';
frame_rate = 30; %enter the frame rate, or I can edit this to enter the stimulus duration

%setting my paths

fnIn_base = 'Z:\home\Celine\Data\2p_data\';
fnOut_base = 'Z:\home\Celine\Analysis\2p_analysis\';

%fnIn_base = 'Z:\home\ACh\Data\2p_data\';
%fnOut_base = 'Z:\home\ACh\Analysis\2p_analysis\';


fnIn = fullfile(fnIn_base,mouse,date,imgFolder);
fnOut = fullfile(fnOut_base,mouse,date,imgFolder);
mkdir(fnOut);
cd(fnIn);
%% load data
imgMatFile = [imgFolder '_000_000.mat'];
load(imgMatFile);

beh_prefix = strcat('Z:\Behavior\Data\data-');
beh_file = [beh_prefix mouse '-' date '-' time '.mat'];
load(beh_file); %load the mworks behavioral file

%get the number of frames entered by the 2p user and the number actually
%collected and use whichever is lower
nFrames = min([input.counterValues{end}(end) info.config.frames]);

%determine whether both PMTs were on, or only the green one, then reas in
%the data
fprintf(['Reading run ' imgMatFile '; ' num2str(nFrames) ' frames \r\n']);
if info.config.pmt1_gain > 0.5
    fprintf('Both green and red PMTS active');
    data_g = sbxread(imgMatFile(1,1:11),0,min(nFrames));
    data_r = squeeze(data_g(2,:,:,:));
    data_g = squeeze(data_g(1,:,:,:));
else
    fprintf('Only green PMT active');
    data_g = sbxread(imgMatFile(1,1:11),0,min(nFrames));
    data_g = squeeze(data_g(1,:,:,:));
end

fprintf(['Loaded ' num2str(nFrames) ' frames \r\n']);
   
behData = input;

clear beh_prefix input
%% register the green data
if exist(fullfile(fnOut,'regOuts&Img.mat')) %check if there is already registration info for this data
    load(fullfile(fnOut,'regOuts&Img.mat'))
    [~,data_g_reg] = stackRegister_MA(data_g,[],[],double(outs));
    data_avg = mean(data_g_reg,3);
    figure;imagesc(data_avg);colormap gray
    save(fullfile(fnOut,'regOuts&Img.mat'),'outs','regImg','data_avg')    
    save(fullfile(fnOut,'input.mat'),'behData')
    clear data_g input
    
else %if not, must register. Start by showing average for each of four 500-frame bins so user can choose one
    nep = floor(size(data_g,3)./10000);
    [n n2] = subplotn(nep);
    figure; 
    movegui('center')
    for i = 1:nep 
        subplot(n,n2,i); 
        imagesc(mean(data_g(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); 
        title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); 
        colormap gray; 
        %clim([0 3000]); 
    end
    
    regImgStartFrame = input('Enter Registration Image Start Frame:');
    regImg = mean(data_g(:,:,regImgStartFrame:(regImgStartFrame+499)),3);
    [outs,data_g_reg] = stackRegister(data_g,regImg);
    data_avg = mean(data_g_reg,3);
    figure;imagesc(data_avg);colormap gray; truesize; clim([0 3000]);
    print(fullfile(fnOut,'avgFOV.pdf'),'-dpdf','-bestfit')
    clear data_g
    save(fullfile(fnOut,'regOuts&Img.mat'),'outs','regImg','data_avg')
end


%reg red data 
%register the red data from the 920 nm run (same run used for green
%above)to the output of the green registration
if info.config.pmt1_gain > 0.5
    [~,data_r_reg] = stackRegister_MA(data_r,[],[],double(outs));
    redChImg = mean(data_r_reg,3);
    clear data_r clear data_r_reg
end

%% will add something to register the red cells

%% find activated cells
%find number of frames per trial and temporarily reshape data into trials
%overal goal here is to get green data in terms of df/f
nOn = behData.nScansOn;
nOff = behData.nScansOff;
ntrials = size(behData.tGratingDirectionDeg,2); %this is a cell array with one value per trial, so length = ntrials
%dimension of each frame in pixels, and number of frames
sz = size(data_g_reg);
data_tr = reshape(data_g_reg,[sz(1), sz(2), nOn+nOff, ntrials]);
%The trial data frame is in x and y pixels per frmae, then n frames per
%trial then ntrials
fprintf(['Size of data_tr is ' num2str(size(data_tr))])
%find baseline fluorescence, using the second half of the stim off period
%select the frames for the second half of the baseline, then average over  
%those frames
data_f = mean(data_tr(:,:,nOff/2:nOff,:),3); 
% subtract baseline fluorescence frrom raw
data_df = bsxfun(@minus, double(data_tr), data_f); 
%normalize to the baseline fluorescence
data_dfof = bsxfun(@rdivide,data_df, data_f); 
% clear the intermediate data frames

%all the dfof data will be from the green channel, so maybe it doesn't make
%
clear data_f data_df data_tr

%find the stimulus directions
tCons = celleqel2mat_padded(behData.tGratingContrast); %transforms cell array into matrix (1 x ntrials)
Cons = unique(tCons);
nCons = length(Cons);

tSize = celleqel2mat_padded(behData.tGratingDiameterDeg); %transforms cell array into matrix (1 x ntrials)
Sizes = unique(tSize);
nSizes = length(Sizes);
%% segmenting green cells
%create empty matrix with FOV for each direction: nYpix x nXPix x nDir
%we will find the average dfof for each of the directions
data_dfof_avg = zeros(sz(1),sz(2),nCons+nDirs); 
%images for segmentation will go through the different stimuli
figure; movegui('center')
[n, n2] = subplotn(nCons); %function to optimize subplot number/dimensions
for iCon = 1:nCons
    ind = find(tCons == Cons(iCon)); %find all trials with each direction
    %average all On frames and all trials
    data_dfof_avg(:,:,iCon) = mean(mean(data_dfof(:,:,nOff+1:nOn+nOff,ind),3),4);
    subplot(n,n2,iCon)
    imagesc(data_dfof_avg(:,:,iCon))
end

[n, n2] = subplotn(nSizes); %function to optimize subplot number/dimensions
for iSize = 1:nSizes
    ind = find(tSize == Sizes(iSize)); %find all trials with each direction
    %average all On frames and all trials
    data_dfof_avg(:,:,nCons+iSize) = mean(mean(data_dfof(:,:,nOff+1:nOn+nOff,ind),3),4);
    subplot(n,n2,iSize)
    imagesc(data_dfof_avg(:,:,iSize))
end

%make a gaussian filter an apply it to the data, which shoould smooth it
%not sure how to optimize this filter
myfilter = fspecial('gaussian',[20 20], 0.5); %maybe filter less?
data_dfof_avg_filtered = imfilter(data_dfof_avg,myfilter);
figure; movegui('center')
[n, n2] = subplotn(nCons); %function to optimize subplot number/dimensions
for iCon = 1:nCons
    ind = find(tCons == Cons(iCon)); %find all trials with each direction
    %average all On frames and all trials
    
    subplot(n,n2,iCon)
    imagesc(data_dfof_avg_filtered(:,:,iCon))
end
%the plots don't look very different, so I'm not sure what are good filter
%settings
%take the max projection - this give a single "frame" that is the maximum
%value for each pixel, which should amount to the maximum response for each
%cell to any stimulus condition
data_dfof_max = max(data_dfof_avg_filtered,[],3); 

data_reg_avg = mean(data_g_reg,3);
% -->I could also do a pixel correlation and use that for segmenting, too

%tack max projection onto the front of the set of averages
data_dfof_for_masks = cat(3,data_dfof_max, data_dfof_avg_filtered,data_reg_avg);
%find the size of this data frame in order to make mask data frames of the
%same dimensions
sz = size(data_dfof_for_masks);
%make empty 2D dataframes to hold mask data moving forward
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));

%loop through the max projection and then the average for each stim
%condition and show that image, using theimCellEditInteractiveLG function
%to select cells

for iStim = 1:size(data_dfof_for_masks,3)    
    mask_data_temp = data_dfof_for_masks(:,:,iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0; %blacks out old cells
    bwout = imCellEditInteractiveLG(mask_data_temp); %selection GUI
    mask_all = mask_all+bwout; %adds new cells to old cells
    mask_exp = imCellBuffer(mask_all,3)+mask_all; %creates buffer around cells to avoid fusing.
    %this buffer is carried over to the next iteration of the loop
    close all
end
%masks are 0/1, so mask all is a ypix by x pix array of 0/1 values.
mask_cell = bwlabel(mask_all); 
%turns logical into numbered cells, numbered from upper left to lower rihgt?
figure;
%shows all the cells, color graded by cell number
imagesc(mask_cell)


% extract timecourses before np subtracktion using stackGetTimeCourses
data_tc = stackGetTimeCourses(data_g_reg, mask_cell); %applies mask to stack (averages all pixels in each frame for each cell) to get timecourses

nCells =  max(max(mask_cell));

%make a ring-shaped mask around each cell, which is what we will designate
%as neuropil for that cell
nMaskPix = 5; %thickness of neuropil ring in pixels
nBuffPix = 3; %thickness of buffer between cell and ring
mask_np = imCellNeuropil(mask_cell,nBuffPix,nMaskPix);

% save mask data and clear variables that aren't needed anymore
run_str = catRunName(imgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];
save(fullfile(fnOut, [datemouserun '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_np')
clear mask_data mask_all mask_2 data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all  data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 

%% remove the neuropil

%for np subtraction the time courses will be downsampled to reduce
%noisiness. We will use the downsampled tc's to identify the optimal np
%subtraction, then apply that with the full data set
down = 5; %number of frames to average
data_reg_down = stackGroupProject(data_g_reg,down); %averages every 5 frames in stack  
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);    

%make an empty matrix the size of the full data
np_tc = zeros(nFrames,nCells);
%make an empty matrix the size of the downsampled data
np_tc_down = zeros(floor(nFrames./down), nCells);
%fill in those matrices: for each cell, create a corresponding timecourse
%for the np around that cell. Do this for the full and the downsampled data
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_g_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '%s\n']) 
end

%make a list of the weights to try for the np
ii= 0.01:0.01:1;
%make am empty matric to store the skewness of each np weight
x = zeros(length(ii), nCells);
%populate that matrix with the skewness by subtracting the np tc's,
%multiplied by each weight, from the corresponding cell tc's
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
%identify the maximum skew and the weight at which that skew occured for
%each cell
[max_skew, ind] =  max(x,[],1); 
% multiple those best weights by 0.01, but I'm not sure why
np_w = 0.01*ind; 
%use this to subtract the np timecourses at full sampling rate from the 
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
clear data_reg data_reg_down

save(fullfile(fnOut, [datemouserun '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

%% Now we have timecourses for each cell with the neuropil removed, can start analyzing them
% looking at time courses
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
vline(nOff,'g')
title('Timecourses');
hold off

print(fullfile(fnOut, [datemouserun 'timecourses']),'-dpdf');

