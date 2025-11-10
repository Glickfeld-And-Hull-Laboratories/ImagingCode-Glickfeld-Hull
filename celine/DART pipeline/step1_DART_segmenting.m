clear all; clear global;  close all
clc
prompt = 'Enter name of instructions file: '; % enter instr filename without .m
instr = input(prompt, 's');
clear prompt
run(instr);

ds=instructions.ds;
run(ds); 

dataStructLabels = {'contrastxori'}; %enter the variable in your datasheet that indicates the folder for the run where you showed stimuli

rc = behavConstsDART; %this sets directory pathes based on the user's netID and login server

doMWCmPD = true; % generate the MW counter - photodiode counter plot or not

day_id = str2double(instructions.session);
if length(expt) < day_id
    error('Day_id %d not valid for this dataset', num2str(day_id));
else
    fprintf('Analyzing sessions: %s\n', num2str(day_id));
end

%% load data for day

mouse = expt(day_id).mouse;
expDate = expt(day_id).date;
ExperimentFolder = expt(day_id).exptType;

fn = fullfile(rc.analysis,ExperimentFolder,mouse,expDate); %can make this flexible if folder structure is different
mkdir(fn)

runs = eval(['expt(day_id).' cell2mat(dataStructLabels) '_runs']);
times = eval(['expt(day_id).' cell2mat(dataStructLabels) '_time']);
% nruns = length(runs);
runFolder = [];
% for irun = 1:nruns % this is only needed to concat runs
imgFolder = runs{1};
% if nruns == 1
    runFolder = imgFolder;
    fnout = fullfile(fn,runFolder);
    mkdir(fullfile(fn,runFolder))
    fName = [imgFolder '_000_000'];
% elseif irun < nruns
%     runFolder = [runFolder '_' imgFolder];
%     fName = [imgFolder '_000_000'];
% else
%     runFolder = [runFolder '_' imgFolder];
%     fnout = fullfile(fn,runFolder);
%     mkdir(fullfile(fn,runFolder))
%     fName = [imgFolder '_000_000'];
% end

root = rc.data;
CD = fullfile(root, mouse, expDate, runFolder);
dat = 'data-';
cd(CD);

imgMatFile = [imgFolder '_000_000.mat'];
load(imgMatFile);
tHostname = lower(hostname);
[s,tUsername] = dos('ECHO %USERNAME%');
switch tHostname
    case {'nuke'}
         fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\' dat mouse '-' expDate '-' times{1} '.mat'];
    case{'nb-hubel'}
        fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\' dat mouse '-' expDate '-' times{1} '.mat'];     
end

load(fName); %load the mworks behavioral file

% temp(irun) = input; %load the data from the mworks file into temp

nframes = [input.counterValues{end}(end) info.config.frames];

% fprintf(['Reading run ' num2str(irun) '- ' num2str(min(nframes)) ' frames \r\n'])
fprintf(['Reading ' num2str(min(nframes)) ' frames \r\n'])

% if info.config.pmt1_gain > 0.5
%     data_temp_g = sbxread(imgMatFile(1,1:11),0,min(nframes));
%     data_temp_r = squeeze(data_temp_g(2,:,:,:));
%     data_temp_g = squeeze(data_temp_g(1,:,:,:));
% else
    data_temp_g = sbxread(imgMatFile(1,1:11),0,min(nframes));
    data_temp_g = squeeze(data_temp_g(1,:,:,:));
% end

fprintf(['Loaded ' num2str(min(nframes)) ' frames \r\n'])

% if nruns == 1 || irun == 1
    data_g = data_temp_g;
    clear data_temp_g
    % if info.config.pmt1_gain > 0.5
    %     data_r = data_temp_r;
    %     clear data_temp_r
    % end
% else
%     data_g = cat(3, data_g, data_temp_g);
%     clear data_temp_g
    % if info.config.pmt1_gain > 0.5
    %     data_r = cat(3, data_r, data_temp_r);
    %     clear data_temp_r
    % end
% end
% end


% register data for each day
%reg green data

if exist(fullfile(fnout,'regOuts&Img.mat')) %check if there is already registration info for this data
    msgbox("Frame shift info already exists for this session, registration is being done with saved info. Double check this is correct.","Caution","warn")
    load(fullfile(fnout,'regOuts&Img.mat'))
    [~,data_g_reg] = stackRegGPU(data_g,[],[],double(outs));
    data_avg = mean(data_g_reg,3);
    figure;imagesc(data_avg);colormap gray
    save(fullfile(fnout,'regOuts&Img.mat'),'outs','regImg','data_avg')
    % input = concatenateStructuresLG(temp);    
    save(fullfile(fnout,'input.mat'),'input')
    clear data_g
else %if not, must register. Start by showing average for each of four 500-frame bins so user can choose one
    nep = floor(size(data_g,3)./3000);
    [n n2] = subplotn(nep);
    figure; 
    movegui('center')
    for i = 1:nep 
        subplot(n,n2,i); 
        imagesc(mean(data_g(:,:,1+((i-1)*3000):500+((i-1)*3000)),3)); 
        title([num2str(1+((i-1)*3000)) '-' num2str(500+((i-1)*3000))]); 
        colormap gray; 
        clim([100 3000]); 
    end
    drawnow;

    %temporarily change name of the input structure
    inputTemp = input;
    clear input;

    regImgStartFrame = input('Enter Registration Image Start Frame, ENTER INTO DS:');
    input = inputTemp; %change the name back
    clear inputTemp

    regImg = mean(data_g(:,:,regImgStartFrame:(regImgStartFrame+499)),3);
    [outs,data_g_reg] = stackRegGPU(data_g,regImg);
    data_avg = mean(data_g_reg,3);
    figure;imagesc(data_avg);colormap gray; truesize; clim([100 3000]);
    title('Average FOV of Registered Stack')
    print(fullfile(fnout,'avgFOV.pdf'),'-dpdf','-bestfit')
    clear data_g 
    save(fullfile(fnout,'regOuts&Img.mat'),'outs','regImg','data_avg')
    % input = concatenateStructuresLG(temp);    
    save(fullfile(fnout,'input.mat'),'input')
end

    
%% find activated cells
%find number of frames per trial and temporarily reshape data into trials
%overal goal here is to get green data in terms of df/f
nOn = input.nScansOn;
nOff = input.nScansOff;
sz = size(data_g_reg);
ntrials = input.stopAfterNTrials;


% find coarse dfof to calculate active cells
data_g_trial = reshape(data_g_reg, [sz(1) sz(2) nOn+nOff ntrials]);
data_g_f = squeeze(mean(data_g_trial(:,:,nOff/2:nOff,:),3));
data_g_on = squeeze(mean(data_g_trial(:,:,nOff+2:nOff+nOn,:),3));
data_g_dfof = (data_g_on-data_g_f)./data_g_f;
clear data_g_trial data_g_on data_g_f

%find the different directions and sizes
tDir = celleqel2mat_padded(input.tGratingDirectionDeg(1:ntrials));
dirs = unique(tDir);
nDir = length(dirs);

tSize = celleqel2mat_padded(input.tGratingDiameterDeg(1:ntrials));
sizes = unique(tSize);
nSize = length(sizes);
ind_dir = find(tDir == max(sizes(:)));

data_g_size = zeros(sz(1),sz(2), nSize);
data_temp = zeros(sz(1),sz(2), nSize, nDir);
[n n2] = subplotn(nSize);
figure; movegui('center');

for iSize = 1:nSize %for every size
    data_g_dir = zeros(sz(1),sz(2), nDir); %make a data frame is the size of one imaging frame X the number of directions
    for iDir = 1:nDir %for every direction
        ind_dir = find(tDir == dirs(iDir)); %find the indices of trials with that direction
        ind_size = intersect(ind_dir, find(tSize == sizes(iSize)));%find ind of intersection between that direction and the size we're looking at
        data_g_dir(:,:,iDir) = mean(data_g_dfof(:,:,ind_size),3); %pull out all those trials and average over the trials
        data_temp(:,:,iSize,iDir,nDir) = mean(data_g_dfof(:,:,ind_size),3);
    end

    data_g_size(:,:,iSize) = max(data_g_dir,[],3); %find the direction with the max response for that size (max projection across dirs)
    subplot(n,n2,iSize)
    imagesc(data_g_size(:,:,iSize))
    title(num2str(sizes(iSize)))
end

%print image

data_size_max = max(data_g_size,[],3); % find the max projection across sizes of the max projection across dirs
data_dfof = cat(3, data_size_max,data_size_max,data_size_max,data_g_size);
figure; imagesc(data_size_max); movegui('center');title('data size max');
clear data_g_dfof


%get pixel correlation image, another method to identify cells
data_g_down = stackGroupProject(data_g_reg,100);
corrImg = getPixelCorrelationImage(data_g_down);
figure; imagesc(corrImg); movegui('center');title('pixel correlation');
data_dfof = cat(3, data_dfof, data_avg,data_avg,corrImg,corrImg);
clear data_g_down



%% load red cells
%this is where we use the 1040, 1000-frame run
if exist(fullfile(fnout,'redImage.mat'))
    load(fullfile(fnout,'redImage'))
elseif ~isempty(expt(day_id).redChannelRun) %if there IS a red channel run, find and load it
    redRun = expt(day_id).redChannelRun;
    imgMatFile = [redRun '_000_000.mat'];
        root = rc.data;
        cd(fullfile(root, mouse, expDate, redRun));

    load(imgMatFile);

    % fprintf(['Reading run ' num2str(irun) '- ' num2str(info.config.frames) ' frames \r\n'])
    fprintf(['Reading ' num2str(info.config.frames) ' frames \r\n'])
    data_temp = sbxread(imgMatFile(1,1:11),0,info.config.frames);
    if size(data_temp,1) == 2
        data_rg = squeeze(data_temp(1,:,:,:));
        data_rr = squeeze(data_temp(2,:,:,:));
    else
        data_rr = squeeze(data_temp(1,:,:,:));
    end
    clear data_temp

    if exist('redChImg')
        [out, data_rr_reg] = stackRegGPU(stackGroupProject(data_rr,100), redChImg);
        redChImg = mean(data_rr_reg,3);
        disp('used option 1');
    elseif info.config.pmt0_gain>0.5 %if there is a green channel in this run, it gets registered to the registration image from green channel from the 920 run
        %data_rr = padarray(data_rr,9,0,'pre');
        redAvg = mean(data_rr,3);
        [out, data_rr_reg] = stackRegGPU(data_rr,redAvg);
        [~, data_rg_reg] = stackRegGPU(data_rg,[],[],out);
        redChImgTemp = mean(data_rr_reg,3); 
        rg_avg = mean(data_rg_reg,3);
        [out2, ~] = stackRegGPU(rg_avg,data_avg);
        [~,redChImg]=stackRegGPU(redChImgTemp,[],[],out2);
        disp('used option 2');
%         [out, data_rg_reg] = stackRegGPU(data_rg,data_avg); %register the green channel from the 1040 run to the green channel from the 920 run
%         [~, data_rr_reg]=stackRegGPU(data_rr,[],[],out); %use those shifts to register the red 1040 run
%         redChImg = mean(data_rr_reg,3);

        
    else %if there is no green channel in this run
        redAvg = mean(data_rr,3);
        [out, data_rr_reg] = stackRegGPU(data_rr,redAvg);
        redChImgTemp = mean(data_rr_reg,3);
        [~,redChImg] = stackRegGPU(redChImgTemp,data_avg);
        disp('used option 3');
    end
    
    figure; colormap gray; imagesc(redChImg);  movegui('center');title('registration image for red channel');
   
    rgb = zeros(sz(1),sz(2),3);
    rgb(:,:,1) = redChImg./max(redChImg(:));
    rgb(:,:,2) = regImg./max(regImg(:));
    figure; image(rgb);  movegui('center')
    title('Green-920 + Red-1040')
    print(fullfile(fnout,'red_green_FOV.pdf'),'-dpdf','-bestfit')

    
    save(fullfile(fnout,'redImage'),'redChImg')
elseif isempty(expt(day_id).redChannelRun) %if there is NOT a red channel run make a dummy that is all zeros
    redChImg = zeros(size(regImg));
end


%create red image where any pixel value above a certain percentile of the max is set to 90%
%of the max - removing the highest 10% of pixel values to create a lower
%direction image for segmenting
threshPercentile = 99;

highValues = find(redChImg>prctile(redChImg,threshPercentile,'all'));
redThresh = redChImg;
redThresh(highValues)=prctile(redChImg,threshPercentile,'all');
figure; imagesc(redChImg);colormap gray;
figure; imagesc(redThresh);colormap gray;

%clear data_rr data_rg data_rg_reg data_rr_reg

% %% troubleshoot movie
% 
%  writerObj = VideoWriter('myVideo.avi');
%  writerObj.FrameRate = 15;
% 
%  % % set the seconds per image
%  % secsPerImage = [5 10 15];
%  % open the video writer
%  open(writerObj);
%  u8 = im2uint8(data_rr);
%  % write the frames to the video
%  for u=1:size(data_rr,3)
%      % convert the image to a frame
%      frame = im2frame(u8(:,:,u));
%      writeVideo(writerObj, frame);
%  end
%  % close the writer object
%  close(writerObj);

%% segment cells
close all

redForSegmenting = cat(3, redThresh,redThresh,redThresh); %make a dataframe that repeats the red channel image multiple times
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
%find and label the red cells - this is the first segmentation figure that
%comes up
if ~isempty(expt(day_id).redChannelRun)
    for iStim=1:size(redForSegmenting,3)
        mask_data_temp=redForSegmenting(:,:,iStim);
        mask_data_temp(find(mask_exp >= 1)) = 0;
        bwout = imCellEditInteractiveLG(mask_data_temp);
        mask_all = mask_all+bwout;
        mask_exp = imCellBuffer(mask_all,3)+mask_all;
        close all
    end
end

%this version does not pad the red cells


mask_cell_red = bwlabel(mask_all);
mask_data = data_dfof; %this is the registered data from the 920 run
%after making masks for all the red cells, go through different stimuli and
%identify cells that are visible for each

for iStim = 1:size(data_dfof,3)
    mask_data_temp = mask_data(:,:,iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0;
    bwout = imCellEditInteractiveLG(mask_data_temp);
    mask_all = mask_all+bwout;
    mask_exp = imCellBuffer(mask_all,3)+mask_all;
    close all
end
mask_cell = bwlabel(mask_all);
figure; imagesc(mask_cell) 

nCells = max(mask_cell(:));
%mask_label = ones(1,nCells); %this is for the EMX and similar lines only
mask_label = zeros(1,nCells);
for i = 1:nCells
    if mask_cell_red(find(mask_cell == i, 1))
        mask_label(1,i) = 1;
    end
end

mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile(fnout, 'mask_cell.mat'), 'data_dfof', 'mask_cell', 'mask_cell_red', 'mask_np','mask_label')


rgb = zeros(sz(1),sz(2),3);
    rgb(:,:,1) = redChImg./(max(redChImg(:))*.5);
    rgb(:,:,2) = regImg./max(regImg(:));
    figure; image(rgb);  movegui('center')

% hold on
% bound = cell2mat(bwboundaries(mask_cell_red));
% plot(bound(:,2),bound(:,1),'.','color','b','MarkerSize',2);
% hold off

%% extract timecourses


data_tc = stackGetTimeCourses(data_g_reg, mask_cell);
nCells = size(data_tc,2);
data_tc_down = stackGetTimeCourses(stackGroupProject(data_g_reg,5), mask_cell);
clear np_tc np_tc_down
sz = size(data_g_reg);
down = 5;
data_reg_down  = stackGroupProject(data_g_reg,down);
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_g_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['     Cell #' num2str(i) '%s/n']) 
end
%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
[max_skew ind] =  max(x,[],1);
np_w = 0.01*ind;
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);

save(fullfile(fnout, 'TCs.mat'), 'data_tc','np_tc','npSub_tc')

clear data_g_reg data_reg_down



%% find stim info from photodiode

switch instructions.tIdxSource
    case 'PD'
        cd(CD);
        load([runFolder '_000_000.mat']);
        [stimOns stimOffs] = photoFrameFinder_Sanworks(info.frame);
        % tCon = cell2mat(input.tGratingContrast);
        input.stimOns_photodiode = stimOns;
        input.stimOffs_photodiode = stimOffs;
        input.stimOns_mwCounter = [];
        input.stimTimingSource = 'PD';
    case 'MW'
        input_correct = counterValCorrect_noPhotodiode(input);
        input.stimOns_mwCounter = input_correct.cStimOn;
        input.counterValues = input_correct.counterValues;
        input.counterTimesUs = input_correct.counterTimesUs;
        input.stimOns_photodiode = [];
        input.stimOffs_photodiode = [];
        input.stimTimingSource = 'MW';
        stimOns = input.stimOns_mwCounter;
        clear input_correct
    otherwise
        error('No valid trial indexing source specificed in instr file. Use "PD" or "MW".');
end

[nFrames nCells] = size(npSub_tc);
nOn = input.nScansOn(1);
nOff = input.nScansOff(1);

% stimOns = stimOns_correct;
% stimOffs = stimOffs_correct;

save('input.mat','input') %resave the input structure named input that includes the photodiode-identified trial start times
nTrials = length(stimOns);
data_tc = nan(nOn+nOff,nCells,nTrials);
for itrial = 1:nTrials
  if ~isnan(stimOns(itrial)) & (stimOns(itrial)+nOn+nOff/2)<nFrames
    data_tc(:,:,itrial) = npSub_tc(stimOns(itrial)-nOff/2:stimOns(itrial)-1+nOn+nOff/2,:);
  end
end
%% dfof calculation
%data_tc = data_tc(:,:,1:end-1);
data_f_trial = nanmean(data_tc(1:nOff/2,:,:),1);
data_dfof_trial = bsxfun(@rdivide, bsxfun(@minus,data_tc, data_f_trial), data_f_trial);
data_dfof_trial = permute(data_dfof_trial,[1 3 2]); %nFrames x nTrials x nCells
% plot tc sanity check
% tavg = squeeze(nanmean(data_dfof_trial,3));
% figure;
% hold on
% for icell = 1:nCells
%     plot(tavg(:,icell), 'LineWidth',.005,'color',[.25 .25 .25]);
% end
% tc_cellavg = mean(tavg,2);
% plot(tc_cellavg);


%% calculate responsive cells

tCon = cell2mat(input.tGratingContrast);
contrasts = unique(tCon);
nCon = length(contrasts);
resp_win = nOff/2:nOff/2+nOn;
base_win = 1:nOff/2;

data_resp = zeros(nCells, nSize, nDir,2);
h = zeros(nCells, nSize, nDir);
p = zeros(nCells, nSize, nDir);
for iSize = 1:nSize
    ind_size = find(tSize == sizes(iSize));
    for iDir = 1:nDir
        ind_dir = find(tDir == dirs(iDir));
        ind = intersect(ind_size,ind_dir); %for every size and then every direction, find trials with that dir/size combination
        data_resp(:,iSize,iDir,1) = squeeze(mean(mean(data_dfof_trial(resp_win,ind,:),1),2));
        data_resp(:,iSize,iDir,2) = squeeze(std(mean(data_dfof_trial(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
        [h(:,iSize,iDir), p(:,iSize,iDir)] = ttest(mean(data_dfof_trial(resp_win,ind,:),1), mean(data_dfof_trial(base_win,ind,:),1),'dim',2,'tail','right','alpha',0.05./(nDir*nCon*nSize-1));
    end
end
%%
h_all = sum(sum(h,2),3);

resp=logical(h_all);
red=mask_label';
resp_red=logical(resp.*red);
sum(resp)
sum(resp_red)

if length(find(n))<36
    [n n2] = subplotn(length(find(h_all)));
    tot = n.*n2;
else
    n = 6;
    n2 = 6;
    tot = 36;
end

% plot size tuning at each direction
figure;
movegui('center')
start = 1;
pref_size = zeros(1,nCells);
for iCell = 1:nCells
    if start>tot
        figure; movegui('center')
        start = 1;
    end
    subplot(n,n2,start)
    %if find(find(h_all)==iCell)
        for iDir = 1:nDir
            errorbar(sizes, data_resp(iCell,:,iDir,1), data_resp(iCell,:,iDir,2),'-o')
            hold on
        end
        if find(find(mask_label)==iCell)
            title('R')
        end
        start= start+1;
        ylim([-0.1 inf])
    %end
    [max_val, pref_size(1,iCell)] = max(mean(data_resp(iCell,:,:,1),3),[],2);
end

% plots direction preference at preferred size
figure;
movegui('center')
start = 1;
for iCell = 1:nCells 
    if start>tot
        figure; movegui('center')
        start = 1;
    end
    subplot(n,n2,start)
    if find(find(h_all)==iCell)
        errorbar(dirs, squeeze(data_resp(iCell,pref_size(iCell),:,1)), squeeze(data_resp(iCell,pref_size(iCell),:,2)),'-o')
        if find(find(mask_label)==iCell)
            title('R')
        end
        ylim([-0.1 inf])
        start = start+1;
    end
end


%% looking at time courses
red_tcs = npSub_tc(:,find(mask_label));
green_inds = 1:nCells;
green_inds = setdiff(green_inds, find(mask_label));
green_tcs = npSub_tc(:,green_inds);

% data_tc_trial = reshape(npSub_tc, [nOn+nOff,nTrials,nCells]);

% data_f_trial = mean(data_tc_trial(nOff/2:nOff,:,:),1);
% data_dfof_trial = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial, data_f_trial), data_f_trial);

%looking at data with np subtracted
tc_cell_avrg = nanmean(data_dfof_trial(:,1:400,resp),3);%average pver cells, one row per trial
tc_trial_avrg = squeeze(nanmean(data_dfof_trial(:,1:400,resp),2));%average over trials, one row per cell
tc_cell_trial_avrg = nanmean(tc_cell_avrg,2);%average over trials and cells

figure;
plot(tc_trial_avrg, 'LineWidth',.005,'color',[.25 .25 .25]);
hold on;
plot(tc_cell_trial_avrg, 'LineWidth',2, 'color','k');
hold on;
vline(nOff/2,'g')
title('');
hold off
ylim([-.02 .18])


% addition 07/05/24 to save semi-raw tc and state # of cells
title(['Responsive cells (',num2str(sum(resp)), ' total, ', num2str(sum(resp_red)), ' red), out of ', num2str(size(data_tc,2)), ' total cells']);
print(fullfile(fnout,'rawTCs.pdf'),'-dpdf','-bestfit');

%% add diagnostic plots

%[stimOns stimOffs] = photoFrameFinder_Sanworks(info.frame);
nTrials = length(stimOns);
nCells = size(data_dfof_trial,3);
trialLength = nOn+nOff;
nFrames = trialLength*nTrials;
AllTrialsNFrames = nan(nTrials,1);
l1 = length(stimOns);
l2 = length(stimOffs);
if l1 > l2
    stimOffs(end+1) = nFrames;
end
%MAXnFrames
for itrial = 1:nTrials
    if itrial == 1
        startFrame = 0;
    else
        startFrame = stimOffs(itrial-1);
    end
    iOn_pd = stimOns(itrial);
    nOff_pd = iOn_pd - startFrame; 
    nOn_pd = stimOffs(itrial) - stimOns(itrial);
    thisTrialNFrames = nOn_pd+nOff_pd;
    AllTrialsNFrames(itrial) = thisTrialNFrames;
end
MAXnFrames = max(AllTrialsNFrames);
data_tc = nan(MAXnFrames,nCells,nTrials);
data_f_trial = nan(nCells,nTrials);
%stimOffs(end+1) = nFrames;
for itrial = 1:nTrials
    if itrial == 1
        firstFrame = 0;
        startFrame = 1; % for indexing later
    else
        firstFrame = stimOffs(itrial-1);
        startFrame = firstFrame;
    end
    iOn_pd = stimOns(itrial);
    nOff_pd = iOn_pd - startFrame; 
    nOn_pd = stimOffs(itrial) - stimOns(itrial);
    thisTrialNFrames = nOn_pd+nOff_pd;
    trialFrames = nan(MAXnFrames,nCells);
    trialFrames(1:nOn_pd+nOff_pd,:) = npSub_tc(startFrame:stimOffs(itrial)-1,:);
    %if ~isnan(stimOns(itrial)) & (stimOns(itrial)+nOn+nOff/2)<nFrames
        data_tc(:,:,itrial) = trialFrames;
    %end
end

elements = zeros(size(data_tc,3),1);
for m = 1:nTrials
    elements(m,1) = sum(~isnan(data_tc(:,1,m)));
end
% plot 
figure
[counts centers] = hist(elements);
histogram(elements);
set(gca,'YScale','log')
sgtitle('Log Trial Length Distribution')
ylim([0.1 max(n)+max(n)*1/10])
xlabel('nFrames')
ylabel('log number of trials')

print(fullfile(fnout,'LogTrialLengthDistribution.pdf'),'-dpdf','-bestfit');

figure
histogram(elements);
sgtitle('Trial Length Distribution')
ylim([0.1 max(n)+max(n)*1/10])
xlabel('nFrames')
ylabel('number of trials')

print(fullfile(fnout,'TrialLengthDistribution.pdf'),'-dpdf','-bestfit');

if doMWCmPD == true
    figure
    ptdcounter = stimOns;
    mwcounter = double([60:90:nFrames]);
    plot(ptdcounter-mwcounter);
    sgtitle('Counter Value Diff')
    xlabel('trial number')
    ylabel('Photodiode Counter - mWorks Counter')
    print(fullfile(fnout,'ptd-mw.pdf'),'-dpdf','-bestfit');
else
    fprintf('No MWorks Counter - Photodiode Counter plot generated');
end


%% look at running trials

loc_dat = input.wheelSpeedValues;
wheelspd = zeros(length(loc_dat),1);
for iTrial = 1:length(wheelspd)
    wheelspd(iTrial) = mean(loc_dat{iTrial});
end

figure
plot(wheelspd);
sgtitle('Trial Average Wheel Speed')
print(fullfile(fnout,'wheelspd.pdf'),'-dpdf','-bestfit');

runidx = find(wheelspd>2);
runTrialBoolean = wheelspd > 2;
tCon = cell2mat(input.tGratingContrast);
contrasts = unique(tCon);
nCon = length(contrasts);
tot_conds = nSize * nDir * nCon;
all_conds = cell(tot_conds,1);
RunTrialsN = nan(tot_conds,1);
cond_counter = 1;
%all_idx = [];

for iSize = 1:nSize
    ind_size = find(tSize == sizes(iSize));
    for iDir = 1:nDir
        ind_dir = find(tDir == dirs(iDir));
        for iCon = 1:nCon
            ind_con = find(tCon == contrasts(iCon));
            this_ind = intersect(intersect(ind_size,ind_dir,'stable'),ind_con,'stable'); 
            this_condition = ['Con-' num2str(contrasts(iCon)) '-Dir-' num2str(dirs(iDir)) '-Size-' num2str(sizes(iSize))];
            all_conds{cond_counter} = this_condition;
            haveRun = intersect(this_ind,runidx);
            RunTrialsN(cond_counter) = length(haveRun);
            cond_counter = cond_counter + 1;
            %all_idx = [all_idx length(this_ind)];
        end
    end
end

figure
bar(RunTrialsN)
sgtitle('Number of Running Trials in Each Stimulus Condition')
ylim([0 max(RunTrialsN)+1])
yticks([0:1:max(RunTrialsN)+1])
xticks([1:1:length(all_conds)])
xticklabels(all_conds)
ylabel('# of Running Trials')

ax = gca;
ax.FontSize = 8; 
print(fullfile(fnout,'nRunTrialsInCond.pdf'),'-dpdf','-bestfit');