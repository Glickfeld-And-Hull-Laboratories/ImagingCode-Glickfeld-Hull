% to analyze individual day's data
clear all; clear global; close all

%identifying animal and run
date = '211004';
imgFolder = '002';
time = '1442';
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
   

clear beh_prefix 
%% register the green data
if exist(fullfile(fnOut,'regOuts&Img.mat')) %check if there is already registration info for this data
    load(fullfile(fnOut,'regOuts&Img.mat'))
    [~,data_g_reg] = stackRegister_MA(data_g,[],[],double(outs));
    data_avg = mean(data_g_reg,3);
    figure;imagesc(data_avg);colormap gray
    save(fullfile(fnOut,'regOuts&Img.mat'),'outs','regImg','data_avg')    
    save(fullfile(fnOut,'input.mat'),'input')
    clear data_g input
    
else %if not, must register. Start by showing average for each of four 500-frame bins so user can choose one
    nep = floor(size(data_g,3)./1000);
    [n n2] = subplotn(nep);
    figure; 
    movegui('center')
    for i = 1:nep 
        subplot(n,n2,i); 
        imagesc(mean(data_g(:,:,1+((i-1)*1000):500+((i-1)*1000)),3)); 
        title([num2str(1+((i-1)*1000)) '-' num2str(500+((i-1)*1000))]); 
        colormap gray; 
        %clim([0 3000]); 
    end
    beh_struct = input; clear input
    regImgStartFrame = input('Enter Registration Image Start Frame:');
    regImg = mean(data_g(:,:,regImgStartFrame:(regImgStartFrame+499)),3);
    [outs,data_g_reg] = stackRegister(data_g,regImg);
    data_avg = mean(data_g_reg,3);
    figure;imagesc(data_avg);colormap gray; truesize; clim([0 3000]);
    print(fullfile(fnOut,'avgFOV.pdf'),'-dpdf','-bestfit')
    clear data_g
    save(fullfile(fnOut,'regOuts&Img.mat'),'outs','regImg','data_avg')
    input = beh_struct;
end


%reg red data 
%register the red data from the 920 nm run (same run used for green
%above)to the output of the green registration
if info.config.pmt1_gain > 0.5
    [~,data_r_reg] = stackRegister_MA(data_r,[],[],double(outs));
    redChImg = mean(data_r_reg,3);
    clear data_r clear data_r_reg
end


%% find activated cells


%find the relevant parameters from the input structure
ntrials = size(input.tThisTrialStartTimeMs,2); %this is a cell array with one value per trial, so length = ntrials
%dimension of each frame in pixels, and number of frames
sz = size(data_g_reg);

%right now I'm only doing one contrast
stimOneContrast = celleqel2mat_padded(input.tStimOneGratingContrast);
stimOneCons = unique(stimOneContrast);

%right now I'm only doing one orientation
ori_mat = celleqel2mat_padded(input.tStimOneGratingDirectionDeg);
oris = unique(ori_mat);
nOri = length(oris);

stimOneTime = celleqel2mat_padded(input.tStimOneGratingOnTimeMs); %duration for each trial
stimOneTimes = unique(stimOneTime); %list of durations used
nTime = length(stimOneTimes); %how many different durations were used
nTrials = length(cStimOne);

%empty data frames that are the size of the FOV by the number of trials
data_f = nan(sz(1),sz(2),nTrials); 
data_one = nan(sz(1),sz(2),nTrials);
%loop through each trial
for itrial = 1:nTrials
    if (cStimOne(itrial)+30)<sz(3) %make sure trials isn't too close to the end
        data_f(:,:,itrial) = mean(data_g_reg(:,:,cStimOne(itrial)-20:cStimOne(itrial)-1),3); %baseline F is the 20 frames before the stim
        data_one(:,:,itrial) = mean(data_g_reg(:,:,cStimOne(itrial):cStimOne(itrial)+round(stimOneTime(itrial)*(double(input.frameRateHz)/1000))),3); %response F is the stimulus duration
    end
end
data_one_dfof = (data_one-data_f)./data_f;


if input.doRandStimOnTime
    data_dfof_stim = zeros(sz(1),sz(2),nTime+1);
    [n n2] = subplotn(nTime+1);
    figure;
    for it = 1:nTime %loop through the stim durations
        ind = find(stimOneTime == stimOneTimes(it)); %finr trials with that duration
        data_dfof_stim(:,:,it) = nanmean(data_one_dfof(:,:,ind),3); %average dfof for each duration
        subplot(n,n2,it)
        imagesc(data_dfof_stim(:,:,it))
        title(num2str(stimOneTimes(it)))
    end
end

data_dfof = cat(3,data_dfof_stim, mean(data_g_reg,3),double(max(data_g_reg,[],3)));


%% segmenting green cells

mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
mask_data = data_dfof;

for iStim = 1:size(data_dfof,3)
    mask_data_temp = mask_data(:,:,end+1-iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0;
    bwout = imCellEditInteractiveLG(mask_data_temp);
    mask_all = mask_all+bwout;
    mask_exp = imCellBuffer(mask_all,3)+mask_all;
    close all
end
mask_cell= bwlabel(mask_all);
figure; imagesc(mask_cell)

%% extract timecourses before np subtracktion using stackGetTimeCourses
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
save(fullfile(fnOut, [datemouserun '_mask_cell.mat']), 'mask_cell', 'mask_np')
clear mask_data mask_all

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

%split into trials using maximum length
nOnMax = max(stimOneTime)*(frame_rate/1000); %find the maximum number of frames for the on period
nFrameTrial = nOnMax+3*frame_rate; %to add one second on either side

data_tc_trial = nan(nFrameTrial,nTrials,nCells); %make empty datafrmae
for iTrial = 1:nTrials
    if cStimOne(iTrial)+nOnMax+frame_rate < nFrames %to remove trials too close to the end
        tempF =  mean(data_tc(cStimOne(iTrial)-20 :cStimOne(iTrial)-1,:),1);
        data_tc_trial_temp = data_tc(cStimOne(iTrial)-frame_rate+1 : cStimOne(iTrial)+nOnMax+(2*frame_rate),:);
        data_tc_trial(:,iTrial,:) = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial_temp, tempF), tempF);
     end
end

save(fullfile(fnOut, [datemouserun '_trial_TCs.mat']), 'data_tc_trial')



%% identifying visually responsive cells
%for now, I will use the number of frames corresponding to the min duration
nOnMin =  round(min(stimOneTimes)*frame_rate/1000);
resp_win = frame_rate+1:frame_rate+nOnMin;
base_win = frame_rate/2 : frame_rate;%half second before the stim

data_resp = zeros(nCells,nTime,2);
h = zeros(nCells, nTime);
p = zeros(nCells, nTime);

for  it = 1:nTime
   inds = find(stimOneTime == stimOneTimes(it));
   data_resp(:,it,1) =squeeze(nanmean(nanmean(data_tc_trial(resp_win,inds,:)),2));
   data_resp(:,it,2) =squeeze(std(nanmean(data_tc_trial(resp_win,inds,:)),[],2))./sqrt(length(inds));
   
   [h(:,it), p(:,it)] = ttest(squeeze(nanmean(data_tc_trial(resp_win,inds,:),1)), squeeze(nanmean(data_tc_trial(base_win,inds,:),1)),'dim',1,'tail','right','alpha',0.05./(nOri.*3-1));
    
end

h_all = sum(h,2);
resp=logical(h_all);

%% looking at time courses
t = 1:(size(data_tc_trial,1));
t=(t-(frame_rate))/frame_rate;

[n n2] = subplotn(nTime);
figure;
    for it = 1:nTime %loop through the stim durations
        inds = find(stimOneTime == stimOneTimes(it)); %find trials with that duration
        temp_trials = squeeze(nanmean(data_tc_trial(:,inds,resp),2));
        subplot(n,n2,it)
        plot(t,temp_trials)
        hold on
        plot(t,mean(squeeze(nanmean(data_tc_trial(:,inds,resp),2)),2),'k');
        title(num2str(stimOneTimes(it)))
    end
print(fullfile(fnOut, [datemouserun 'timecourses']),'-dpdf');


figure;
    for it = 1:nTime %loop through the stim durations
        inds = find(stimOneTime == stimOneTimes(it)); %find trials with that duration
        temp_mean = mean(squeeze(nanmean(data_tc_trial(:,inds,resp),2)),2);
        temp_se= std(squeeze(nanmean(data_tc_trial(:,inds,resp),2)),[],2)/sqrt(sum(resp));
        subplot(n,n2,it)
        shadedErrorBar(t,temp_mean,temp_se);
%         vline(frame_rate+1)
%         vline(frame_rate + 1+(stimOneTimes(it)*frame_rate/1000))
        title(num2str(stimOneTimes(it)))
    end
print(fullfile(fnOut, [datemouserun 'shadedEB_timecourses']),'-dpdf');