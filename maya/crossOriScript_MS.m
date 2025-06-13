clc; clear all; close all; %clears everything in command windows, resets all variables
doRedChannel = 0;   %new variable doRedChannel
ds = 'CrossOriRandDirFourPhase_ExptList_MS'; % Pull experimental info from my experiment list
iexp = 2;  % experiment number, I chose a L4 experiment for you to practice on
doPhaseAfterDir = 0;  %new variables
doDirAfterPass = 0;
eval(ds) %evaluates the matlab code in this 

saveDir = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\maya\'; %creates a place to save the experiment info
 
frame_rate = params.frameRate; %pulls the frame rate from the parameters of the experimental info

%%
mouse = expt(iexp).mouse;   %mouse identification(i1404)
date = expt(iexp).date;     %date identification (241029)
area = expt(iexp).img_loc{1};    %area imaged(V1)
if doPhaseAfterDir            %if doPhaseAfterDir exists
    ImgFolder = expt(iexp).coFolder;  %retrieve the imaging folder number
    nrun = length(ImgFolder);      %and the number of runs based on the length of the ImgFolder
    ref_str = catRunName(cell2mat(ImgFolder), nrun);   %create a name string 'cell2mat(002)'
    ImgFolder = expt(iexp).copFolder; %then resaves ImgFolder
    time = expt(iexp).copTime; %and the time the run started
elseif doDirAfterPass      %same but checks
    ImgFolder = expt(iexp).passFolder;
    nrun = length(ImgFolder);
    ref_str = catRunName(cell2mat(ImgFolder), nrun);
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
else
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
end
nrun = length(ImgFolder);  
% run_str = catRunName(cell2mat(ImgFolder), nrun);
run_str = 'runs-002'; 

base = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\' expt(iexp).saveLoc]; %saves the base of the file folder 

fprintf(['2P imaging cross ori analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([ImgFolder{irun} ' - ' time{irun} '\n'])
end

%prints the selected data of the mouse, date, and experiment number and
%time
%% load and register
tic  %starts a timer of how long matlab takes to load and register data
data = [];  %creates an empty array/matrix 'data'
clear temp   %clears variable temp
trial_n = [];    %new matrix trial_n
offset = 0;   %new variable offset with value 0
for irun = 1:nrun     %for each run from 1 to nrun, assign variable base + Data\2P_images\mouse\date\folder(run)
    CD = [base '\Data\2P_images\' mouse '\' date '\' ImgFolder{irun}];
    %CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\AW68\two-photon imaging\' date '\' ImgFolder(irun,:)];
    %CD = [base '\Data\2P_images\' mouse '-KM' '\' date '_' mouse '\' ImgFolder(irun,:)];
    %CD = ['\\CRASH.dhe.duke.edu\data\home\kevin\Data\2P\' date '_' mouse '\' ImgFolder(irun,:)];
    cd(CD); 
    imgMatFile = [ImgFolder{irun} '_000_000.mat']; 
    load(imgMatFile); %load matlab file of run
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time{irun} '.mat'];
    load(fName); %load behavior data for the run
    
    temp(irun) = input;  %load each runs input data from the behavior file
    nframes = [temp(irun).counterValues{end}(end) info.config.frames]; %pulls number of frames collected vs. frames supposed to be collecte from input
    
    fprintf(['Reading run ' num2str(irun) '- ' num2str(min(nframes)) ' frames \r\n']) %prints 'reading run irun' and frames
    data_temp = sbxread(imgMatFile(1,1:11),0,min(nframes));  %data temporary variable is equal to the number of frames starting from zero and the consecutive frames after with N>1
    if size(data_temp,1)== 2 %if data temp spits out a variable that is only size 2, make it a larger matrix (4x)
        data_temp = data_temp(1,:,:,:);
    end
    
    if isfield(input, 'cStimOneOn')  %if one of the fieldss of input structure array is cStimOneOn
        if irun>1 %and if irun is greater than one
            ntrials = size(input.cStimOneOn,2); %make ntrials the value of the second column of cStimOneOn
            for itrial = 1:ntrials %for each trial, 1 through the number of trials
                temp(irun).cStimOneOn{itrial} = temp(irun).cStimOneOn{itrial}+offset; %input.cStimOneOn(each trial) = input.cStimOneOn+offset
                temp(irun).cStimOneOff{itrial} = temp(irun).cStimOneOff{itrial}+offset; %and cStimOneOff+offset as well
            end
        end
    end
    offset = offset+min(nframes); %offset is set to minimum number of frames collected
        
    data_temp = squeeze(data_temp); %get rid of empty dimension
    data = cat(3,data,data_temp); %concatenates data and data_temp to create a 3 dimension matrix
    trial_n = [trial_n nframes]; %transforms trial_n into trial_n + the actual nframes
end
input = concatenateStructuresLG(temp); %concatenates multiple mat files (temp)
clear data_temp
clear temp
toc %finishes counting time to load and register

% Choose register interval
regIntv = 5000; %registered interval
nep = floor(size(data,3)./regIntv);  %variable nep is assigned to the rounded value of the 3rd dimension of 'data' divided by 5000
[n n2] = subplotn(nep); %gives minimum subplot dimensions for nep
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*regIntv):500+((i-1)*regIntv)),3)); title([num2str(1+((i-1)*regIntv)) '-' num2str(500+((i-1)*regIntv))]); colormap gray; clim([0 3000]); end
movegui('center') %plots the registered intervals of 500 frames at a time (pulls each 5000'th frame)
%% Register data
data_avg = mean(data(:,:,85001:85500),3); %creates a matrix where data_avg is the mean of data of the 3rd dimesion on rows 85001 and 85500 (averages over time)
if doPhaseAfterDir || doDirAfterPass      %if either one of these is true 
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_reg_shifts.mat'])) %loads a file with the date, mouse, cell12mat002_reg_shifts.mat (basically the run data)
    [out, data_reg] = stackRegister(data,data_avg);  %returns the file data_reg with the registration of data and data_avg 
    mkdir(fullfile(saveDir, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str])) %creates a folder in my folder with this name
    save(fullfile(saveDir, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg') %saves this file into the created folder
    save(fullfile(saveDir, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input') %also saves this file into the folder
elseif exist(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat'])) %elseif the file with the 2P data exists...
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat'])) %load it
    [outs, data_reg] = stackRegister_MA(data,[],[],out); %and register it to data_reg but only register data with both green and red channels, and this time use the "outs" function so its a by 4 vector concatenated
else
    [out, data_reg] = stackRegister(data,data_avg); %otherwise just register the data normally 
    mkdir(fullfile(saveDir, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str])) %create a folder in my folder
    save(fullfile(saveDir, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg') %save the registered data into a file in the new folder
    save(fullfile(saveDir, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input') %save all the input parameters to a file in the folder
end

% test stability

%create a figure where for each of the 32 segmented stacks, its divided
%into n x n2  grid with axes at i, and displahys it as a colormap of the
%average values of data_reg at registered intervals, and name it 500 +
%whatever the interval is every 5000 frames (32 maps are created)
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*regIntv):500+((i-1)*regIntv)),3)); title([num2str(1+((i-1)*regIntv)) '-' num2str(500+((i-1)*regIntv))]); end
print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_byFrame.pdf']),'-dpdf', '-bestfit') % saves the figure to a certain file
movegui('center') %moves the figure to the center and displays it 
figure; imagesq(mean(data_reg(:,:,1:regIntv),3)); truesize; % just displays one average heat map plot for the entire run
print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit') %saves this as the average

i = 1; %var
sz = size(data_reg); %sz is the size of data reg [529,796,161987]
rg = zeros(sz(1),sz(2),3); % creates a matrix 529x796x[]
first = mean(data_reg(:,:,1+((i-1)*regIntv):500+((i-1)*regIntv)),3); %first frame of the 500 frames
rg(:,:,1) = first./max(first(:)); %image first row or red channel is equal to first 500 frames of each segment/max value of the first 500 frames (normalize)
i = nep; %i now equals 32 (nep)
last = mean(data_reg(:,:,1+((i-1)*regIntv):500+((i-1)*regIntv)),3); %last frame of the 500 frames
rg(:,:,2) = last./max(last(:)); %image second row or green channel is equal to first 500 frames of each segment/max value of the first 500 frames (normalize)
figure; image(rg) %make a figure
movegui('center') %center it 
%% if red channel data
if doRedChannel & (~doPhaseAfterDir & ~doDirAfterPass) %if doRedChannel (which equals 0) and dophaseafterrdir and dodirafterpass are not true
    ImgFolderRed = expt(iexp).redImg; %this variable equals the red image number in the experiment file
    CD = [base '\Data\2P_images\' mouse '\' date '\' ImgFolderRed{1}]; %CD is the base of the file and just pulls the data specifically for this run for the red 
    cd(CD); 
    imgMatFile = [ImgFolderRed{1} '_000_000.mat']; 
    load(imgMatFile); %loads specifically red channel data 

    nframes = info.config.frames; %number of frames
    fprintf(['Reading red data- ' num2str(nframes) ' frames \r\n']) %print that the red data is being read, with the number of frames completed
    data_red = sbxread(imgMatFile(1,1:11),0,nframes); %data_red is equal to reading imgMatFile(red data) the first column(green column) and second column (red column) to 11 from frame 0 to nframes
    data_red_g = squeeze(data_red(1,:,:,:)); %squeeze the first dimension 
    data_red_r = squeeze(data_red(2,:,:,:)); %squeeze the second dimension

    [rg_out rg_reg] = stackRegister(data_red_g,data_avg); %in rg_red, we store the data of data_red_g and data_avg registered to rg out 
    [rr_out rr_reg] = stackRegister_MA(data_red_r,[],[],rg_out); %rr_reg stores the output of coregistering data_red_r to rg_out (if you have a 3rd input, register frames based on input)

    
    data_red_avg = mean(rr_reg,3); %the red average data is the average of rr_reg of the 3rd dimension (average fov)
    figure; imagesc(data_red_avg); %heatmap figure of the average of the red data
    
    ImgFolderRed = expt(iexp).redImg; %resets variable ImgFolderRed to the experiment log for red channel
    CD = [base '\Data\2P_images\' mouse '\' date '\' ImgFolderRed{1}];        
    cd(CD);
    imgMatFile = [expt(iexp).redImg{1} '_000_000.mat'];
    load(imgMatFile); 
    nframes = info.config.frames; %all same as before (loading imgMatFile data for red channel)
    fprintf(['Reading run ' expt(iexp).redImg{1} '- ' num2str(min(nframes)) ' frames \r\n']) %printing number of frames
    data = sbxread(imgMatFile(1,1:11),0,nframes);  %data_red is equal to reading imgMatFile(red data) the first column and second column to 11 from frame 0 to nframes
    if size(data,1) == 2 %if size of first dimension of data = 2
        red_data = squeeze(data(2,:,:,:)); %squeeze the second dimension of data
        green_data = squeeze(data(1,:,:,:)); %squeeze the first dimension
        [out, green_data_reg] = stackRegister(green_data,data_avg); %register green data to the data average, save in green_data_reg
        [out2, red_data_reg] = stackRegister_MA(red_data,[],[],out); %registers red_data to a 4 dimension matrix corregistered to green_data_reg (in case you have a 3rd input)
        red_data_avg = mean(red_data_reg,3); %take the average of the 3rd dimension
        figure; imagesc(red_data_avg) %create a figure for red channel averages
        title('Red')
        green_data_avg = mean(green_data_reg,3);
        figure; imagesc(green_data_avg) %create a figure for green channel averages 
        title('Green')
        if size(expt(iexp).redImg,2) == 2 %if the size of the second dimension of the red image data is 2 
            CD = [base '\Data\2P_images\' mouse '\' date '\' ImgFolderRed{2}];                  
            cd(CD);
            imgMatFile = [expt(iexp).redImg{2} '_000_000.mat'];
            load(imgMatFile);
            nframes = info.config.frames; %same thing as before, store it in a imgMatFile and run registration
            fprintf(['Reading run ' expt(iexp).redImg{2} '- ' num2str(min(nframes)) ' frames \r\n'])
            data = sbxread(imgMatFile(1,1:11),0,nframes);
            red_data = squeeze(data(2,:,:,:));
            green_data = squeeze(data(1,:,:,:));
%             [out, green_data_reg] = stackRegister(mean(green_dat,3),green_data_avg);
%             [outs, red_data_reg] = stackRegister_MA(mean(red_data,3),[],[],out);
%             red_data_avg = red_data_reg;
            [out, green_data_reg] = stackRegister(green_data,green_data_avg);
            [outs, red_data_reg] = stackRegister_MA(red_data,[],[],out);
            red_data_avg = mean(red_data_reg,3);
            figure; imagesc(red_data_avg); movegui('center')
            title('Red') %create a figure for red channel data
        end
        save(fullfile([saveDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_redData.mat']), 'green_data_avg', 'red_data_avg') %save it in my folder
    end

    data_avg = mean(data_reg(:,:,size(data_reg,3)-10000:end),3); %take the average of data_reg in the third dimension from 151987 to the end
    figure; 
    subplot(2,2,1) %create a figure and subdivide into 2x2 figure with axes at 1
    sz = size(data_avg); %take the size of data_avg
    rgb = zeros(sz(1),sz(2),3); %create an empty matrix that is 529x796x3 
    rgb(:,:,1) = red_data_avg./max(red_data_avg(:)); %make the first dimension the red data average divided by the max valkue of the red data average
    imagesc(red_data_avg); %create that figure
    colormap gray %make it grey
    subplot(2,2,2) %new subplot of that data 2x2 with axes at 2
    rgb(:,:,2) = data_avg./max(data_avg(:)); %second dimension of rgb is data avg divided by max value of data avg
    imagesc(rgb); %create that figure
    if size(expt(iexp).redImg,2) == 2 %if the size of the second dimension of the red image data is 2, we have 2 different wavelengths for the PMIs
        title('Red at 1040; Green at 920')
    else
        title('Red at 920; Green at 920') %otherwise we only have one (920)
    end
    print(fullfile([base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit') %print the file name 

end

%% find activated cells
clear data out   %clear all data for outs (stack registration)
cStimOn = celleqel2mat_padded(input.cStimOneOn); %create a new matrix to store the time assigned to when stimulus is on
nTrials = length(cStimOn); %ntrials is how many times stimulus is on
sz = size(data_reg); %size of data registered

data_resp = nan(sz(1),sz(2),nTrials); %nan value array 529x796x2219
data_f = nan(sz(1),sz(2),nTrials); %nan value array 529x796x2219

for itrial = 1:nTrials % for each trial (1-2219)
    if cStimOn(itrial) + 20 < sz(3) % if the time for each trial + 20 is less than 161987...
        data_resp(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)+5:cStimOn(itrial)+20),3); % (RESPONSE) for the value assigned to the trial in data_reg, store the mean of that trial from i+5-i+20 in the 3rd dimension, on vs off windows, baseline vs response, for each trial average image over time
    end
    data_f(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)-15:cStimOn(itrial)),3); % (BASELINE) then in data f index of that trial, store the average of i-15 in the 3rd dimension
end

data_resp_dfof = (data_resp-data_f)./data_f; %CHANGE IN FLOURESCENE YIPPEE
%yeah so im not gonna go through all this (Maya starts losing her sanity
%pt. 1) but basically assign each aspect of the stimulus conditions to a
%variable 
stimCon_all = celleqel2mat_padded(input.tStimOneGratingContrast);
maskCon_all = celleqel2mat_padded(input.tMaskOneGratingContrast);
stimCons = unique(stimCon_all);
maskCons = unique(maskCon_all);
nStimCon = length(stimCons);
nMaskCon = length(maskCons);
maskPhas_all = celleqel2mat_padded(input.tMaskOneGratingPhaseDeg);
maskPhas = unique(maskPhas_all);
nMaskPhas = length(maskPhas);
stimDir_all = celleqel2mat_padded(input.tStimOneGratingDirectionDeg);
maskDir_all = rad2deg(wrapTo2Pi(deg2rad(celleqel2mat_padded(input.tMaskOneGratingDirectionDeg))));
maskDir_all(find(maskDir_all==360)) = 0;
stimDirs = unique(stimDir_all);
nStimDir = length(stimDirs);
maskDirs = unique(maskDir_all);
nMaskDir = length(maskDirs);
maskDiff_all = celleqel2mat_padded(input.tMaskOneGratingDirectionDeg) - celleqel2mat_padded(input.tStimOneGratingDirectionDeg);
maskDiffs = unique(maskDiff_all);
nMaskDiff = length(maskDiffs);
SF_all = celleqel2mat_padded(input.tStimOneGratingSpatialFreqCPD);
SFs = unique(SF_all);
nSF = length(SFs);
TF_all = celleqel2mat_padded(input.tMaskOneGratingTemporalFreqCPS);
TFs = unique(TF_all);
nTF = length(TFs);

if input.doTwoStimTogether %if input.doTwoStimTogether is true (non zero) (plaid)
    maskCon_all = celleqel2mat_padded(input.tStimTwoGratingContrast); %maskCon_all is equal to a matrix of the two stim grating contrasts
    maskCons = unique(maskCon_all); %no repetitions, in sorted order
    nMaskCon = length(maskCons); %lenght of that matrix 
    maskPhas_all = celleqel2mat_padded(input.tStimTwoGratingPhaseDeg); %degrees of two stims matrix
    maskPhas = unique(maskPhas_all); %sorted and unique
    nMaskPhas = length(maskPhas); %length
    maskDir_all = rad2deg(wrapTo2Pi(deg2rad(celleqel2mat_padded(input.tStimTwoGratingDirectionDeg)))); %degrees of direction matrix
    maskDir_all(find(maskDir_all==360)) = 0; %find all 360 degree values
    maskDirs = unique(maskDir_all); %sort and unique
    nMaskDir = length(maskDirs); %length
end

if nStimCon >= 2  || nMaskCon >=2  %if the number of stimulus conditions or mask conditions is greater than or equal to 2
    nStim1 = nStimCon*nMaskCon*nMaskPhas; %assign nStim1 to the product of the # of stimulus conditions x # of mask conditions x # of degrees
    data_dfof = nan(sz(1),sz(2),nStim1); %create new nan matrix that is 529x796x16

    start = 1;
    for is = 1:nStimCon %for is = 1:2 one to end of stimulus phase conditions for contrast (find all the times mark or stim was alone), finding gratings
        ind_stim = find(stimCon_all == stimCons(is)); %find where in the matrix of stimulus conditions, the index of stimulus conditions is in increasing unique order
        for im = 1:nMaskCon % for im = 1:2 and for each mask condition
            ind_mask = find(maskCon_all == maskCons(im)); %the index of that mask condition in unique order
            if im>1 & is>1 %if both stim and mask conditions have more than one instance
                for ip = 1:nMaskPhas % for ip = 1:4 check mask phase in order and find index of where we reach each mask phase (if both are more than one, both are visible, its a plaid, what phase did it start at?)
                    ind_p = find(maskPhas_all == maskPhas(ip));
                    ind_all = intersect(ind_p,intersect(ind_stim,ind_mask)); %where is the data common(where do they all overlap)
                    data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,ind_all),3); %each dimension of data_dfof, return mean of 3rd dimension nan values of change in flourescence data_resp_dfof
                    start = start+1; %add one to start and continue 
                end
            else
                ind_all = intersect(ind_stim,ind_mask); %else just find the intersection of the stimulus and mask
                data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,ind_all),3); %and for each of the conditions
                start = start+1;
            end
        end
    end
else 
    nStim1 = 0; %yeah and if none of that is true set nStim to 0 and move on
    data_dfof = []; %clear data_dfof
end

if nStimDir > 1 & ~input.doTwoStimTogether % if stimulus direction is greater than 1 and you're not doing two stimuli together
    nStim2 = nStimDir.*(nMaskDiff+1); %set nStim2 to the number of stimulus direction times the unique number of mask-stim directions +1
    data_dfof = cat(3, data_dfof, zeros(sz(1),sz(2), nStim2)); %data_dfof now equals the concatenation of data_dfof and zeroes(526x796x24)
    start = nStim1+1; %starting value is 0+1 (17)
    for is = 1:nStimDir %for is is 1:12, find all the times stim was alone and mask was alone
        %ind_stimalone is when stimulus conditions match when the mask
        %condition is 0, and when stimulus directions match the index of
        %the simulus direction (same for mask)
        ind_stimalone = intersect(intersect(find(stimCon_all),find(maskCon_all==0)),find(stimDir_all == stimDirs(is)));
        ind_maskalone = intersect(intersect(find(stimCon_all==0),find(maskCon_all)),find(maskDir_all == maskDirs(is)));
        ind_alone(is) = length([ind_stimalone ind_maskalone]); %concatenate both, how many grating trials alone
        data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,[ind_stimalone ind_maskalone]),3); %data_df_of filled starting with the 1st column of this stupid ahh matrix, get rid of nan values, find mean of data_resp_dfof starting with stim and mask for 3rd dim
        start = start+1; %repeat through all nStim+1 to nStimDir (24)
        for id = 1:nMaskDiff %for id = 1:1, plaid stimulus cross angle potentials, wiht index of stim mask differences
            ind_plaidstim = intersect(intersect(intersect(find(stimCon_all),find(maskCon_all)),find(stimDir_all == stimDirs(is))),find(maskDiff_all == maskDiffs(id)));
            ind_n(is,id) = length(ind_plaidstim); %how many plaids presented
            data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,ind_plaidstim),3); %mean without nans of data_resp_dfof along 3rd dimension, for the column eqaul to ind_plaidstim
            start = start+1;
        end
    end
    figure; %create a figure
    subplot(2,1,1); %2x1 with axes at 1
    imagesc(mean(data_dfof(:,:,nStim1+1:1+nMaskDiff:end),3)) %the image is a heatmap plot of mean of 3rd dimension of data_dfof from the first stimulus condition, to mask, stim difference
    title('Grating') 
    colormap gray
    % clim([0 1])
    subplot(2,1,2); %next part of figure
    imagesc(mean(data_dfof(:,:,nStim1+2:1+nMaskDiff:end),3))% instead with nstim1 + 2 which sets stimulus condition to plaid
    title('Plaid')
    colormap gray
    % clim([0 1])
    sgtitle([mouse ' ' date]) %title
    print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_GratingVsPlaid.pdf']),'-dpdf') %save this in a file
elseif nStimDir > 1 & input.doTwoStimTogether %if stimulus directions (12) > 1 and input.doTwoStimTogether is true
    nStim2 = nStimDir.*2; %nStim2 gets multiplied by two
    data_dfof = cat(3, data_dfof, zeros(sz(1),sz(2), nStim2)); %concatenate data_dfof and an empty matrix 529x796x24 along 3rd dimension
    start = nStim1+1; %start counting from nStim1 0 onward
    for is = 1:nStimDir %for each stimulus direction 1:12
        ind_same = intersect(find(stimDir_all == stimDirs(is)), find(maskDir_all == stimDirs(is))); %the intersection of where the mask direction equals the stimulus direction 1:12, and when the stimulus direction equals that same degree
        data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,ind_same),3); %each dimension of data_dfof filled with the mean of 3rd dimension nan values of the intersection values of mask and stimulus of the 3rd dimension
        ind_both = find(stimCon_all&maskCon_all); %where both stimCon_all and maskCon_all are true
        ind_diff = setdiff(find(stimDir_all == stimDirs(is)), find(maskDir_all == stimDirs(is))); %finds the indexed values of directions 1:12 that are not in maskDir_all but are in stimDir_all
        data_dfof(:,:,start+1) = nanmean(data_resp_dfof(:,:,intersect(ind_both,ind_diff)),3); %now data_dfof starting from the n+1 column is filled with the intersecting values of where both stimulus and mask are the same and when they don't align
        start = start+2; %now start is bumped by two
    end
else
    nStim2 = 0; %set nStim2 to zero
end

if nSF>1 %if number of different spatial frequencies is greater than 1
    nStim3 = nSF; %set nStim3 to  nSF
    data_dfof = cat(3, data_dfof, zeros(sz(1),sz(2), nStim3)); %set data_dfof to the concatenated matrix of data_dfof and an empty matrix of 529x796x1 along the 3rd dimension
    start = nStim1+nStim2+1; %and then set the starting point to 1
    for iSF = 1:nSF %for each different SF value 1:1
        ind_SF = find(SF_all == SFs(iSF)); %find all unique sf values
        data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,ind_SF),3); %and then starting from the next dimension of data_dfof, fill in every trial for the individual SF's per trial along the 3rd dimension
        start = start+1; %then increase starting point
    end
else
    nStim3 = 0; %then set our nStim3 to 0
end

if nTF>1 %next repeat the process for temporal frequency (1) 
    nStim3 = nTF;
    data_dfof = cat(3, data_dfof, zeros(sz(1),sz(2), nStim3));
    start = nStim1+nStim2+1;
    for iTF = 1:nTF
        ind_TF = find(TF_all == TFs(iTF));
        data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,ind_TF),3);
        start = start+1; %we end up with the same temporal frequency for all trials for all cells, which should not change data_dfof
    end
else
    nStim3 = 0; %set nStim3 to zero
end

data_dfof(:,:,isnan(mean(mean(data_dfof,1),2))) = []; %now set data_dfof in any column that is nan to the mean of second dimension of the mean of the first dimension of data_dfof, avg over first and second
myfilter = fspecial('gaussian',[20 20], 0.5); %create a filter for a figure
data_dfof_max = max(imfilter(data_dfof,myfilter),[],3); %puts data_dfof through a gaussian filter, and then returns the max values of every pixel the 3rd dimension
figure;
movegui('center') %create a figure in the center of the screen
imagesc(data_dfof_max) %heatmap of the values of data_dfof_max
print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_maxdfof.pdf']),'-dpdf') %print and save this data

data_dfof = cat(3, data_dfof, data_dfof_max); %return data_dfof to the filtered and maximum data values along the 3rd dimension
if doRedChannel & (~doPhaseAfterDir & ~doDirAfterPass)
    data_dfof = cat(3,data_dfof,red_data_avg); %then if the red channel is true, andnot the two other conditions, concatenate so the red channel data is added to data_dfof
end
if (strcmp(expt(iexp).driver,'SOM') || strcmp(expt(iexp).driver,'PV')) & ~doRedChannel %if these two are equal (SOM and PV) and the red channel condition is not true...
    data_dfof = cat(3,data_dfof,data_avg); %concatenate the two matricies
end

save(fullfile(saveDir, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']), 'cStimOn', 'maskCon_all', 'stimCon_all', 'stimCons', 'maskCons', 'nStimCon', 'nMaskCon', 'stimDir_all', 'stimDirs', 'nStimDir', 'maskDir_all', 'maskDirs', 'nMaskDir', 'maskPhas_all', 'maskPhas', 'nMaskPhas', 'maskDiff_all','maskDiffs','nMaskDiff','SF_all', 'SFs', 'nSF','TF_all', 'TFs', 'nTF', 'frame_rate', 'nTrials') %and save all this data into a file

%% cell segmentation 
if ~doPhaseAfterDir & ~doDirAfterPass %if two conditions are not true (like earlier)
    mask_exp = zeros(sz(1),sz(2)); %create these empty matricies
    mask_all = zeros(sz(1), sz(2));
    mask_data = data_dfof; %and copy data_dfof
    
    if strcmp(expt(iexp).driver{1},'SOM') || strcmp(expt(iexp).driver{1},'PV') %AND if the mice are SOMxPV
        if ~doRedChannel %AND if doRedChannel is not true
            bwout = imCellEditInteractiveLG(mean(data_reg,3)); %add cells to mask by clicking them (like in tutorial)
            mask_all = mask_all+bwout; %creates a mask layer with all the clicked cells
            mask_exp = imCellBuffer(mask_all,3)+mask_all; %buffer mask around selected cells
            close all
        end
    end
    
    for iStim = 1:size(data_dfof,3)   %for each stimulus from 1 to the size of the 3rd dimension of data_dfof
        mask_data_temp = mask_data(:,:,end+1-iStim); %create a temporary data mask that fills values of mask_data_temp from columns at the end minus each stimulus
        mask_data_temp(find(mask_exp >= 1)) = 0; %set every instance of mask_exp that is greater than or equal to one and set the equivalent index to zero in mask_data_temp
        bwout = imCellEditInteractiveLG(mask_data_temp); %now do the same thing to add cells to a mask layer by clicking them
        if doRedChannel & iStim==1 %if doRedChannel is true and iStim is equal to one
            red_mask = bwout; %new red_mask is equal to actual mask
        end
        mask_all = mask_all+bwout; %the mask of all cells is equal to the addition of all the masks together for every stimulus condition direction
        mask_exp = imCellBuffer(mask_all,3)+mask_all; %then the mask buffer is also equal to the entire mask along the third dimension 
        close all
    end
    
    mask_cell= bwlabel(mask_all); %returns the labels of the connected items in mask_all
    if doRedChannel %if red channel is true
        red_cells = unique(mask_cell(find(red_mask))); %find the unique instances of where mask_cell and red_mask correlate
    else
        red_cells = []; %clear red cells
    end
    figure; movegui('center') %create a figure and a heatmap
    imagesc(mask_cell)
end
    

clear data_adapt data_adapt_dfof data_test data_test_dfof data_test_avg data_resp data_resp_dfof bwout %clear all variables

%% neuropil subtraction
 
if ~doPhaseAfterDir & ~doDirAfterPass %if both conditions are not true
    mask_np = imCellNeuropil(mask_cell, 3, 5); %creates a neuropil mask of mask_cell
    save(fullfile(saveDir, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof', 'mask_cell', 'mask_np', 'red_cells') %saves the neuropil mask
else
    load(fullfile(saveDir, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_mask_cell.mat'])) %or just loads the actual cell data
end
clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_data_temp mask_exp data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 
 %clear all these variables

down = 5;
sz = size(data_reg); %sise of data_reg (registered data)

data_tc = stackGetTimeCourses(data_reg, mask_cell); %get the time courses for data_reg and mask_cell stacked
data_reg_down  = stackGroupProject(data_reg,down); %get the downsampled version of data_reg witht the ratio 5
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell); %then get the time courses for the downsampled data_reg and mask_cell
nCells = size(data_tc,2); %collect the number of cells sampled
np_tc = zeros(sz(3),nCells); %empty neuropil matrix in the size of 161987xnCells
np_tc_down = zeros(floor(sz(3)./down), nCells); %gets the downsized version of the size (divided by down, and rounded)
for i = 1:nCells %for every single cell...
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i)); %for each column in np_tc, get the individual timecourse from data_reg corresponding to the neuropil mask for each 3rd dimension
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i)); %do the same thing for the downsampled data
     fprintf(['Cell #' num2str(i) '%s/n'])  %print the number of cells
end
%get weights by maximizing skew
ii= 0.01:0.01:1; %weights
x = zeros(length(ii), nCells); %set x to a matrix of length ii and width of the cell number
for i = 1:100 %for 1-100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i))); %for matrix x, fill the first dimension with a skewness of the downsample minus the original data
end
[max_skew ind] =  max(x,[],1); %set a max skew matrix
np_w = 0.01*ind; %the weight of the skew
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w); %more weight stuff
clear data_reg data_reg_down

save(fullfile(saveDir, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc', 'nCells', 'sz') %save

clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell

%and clear
%BOOM IM DONE 
 