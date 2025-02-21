%%
% We are starting with data that are in nPix x nPix x nFrames for mice that
% have viewed 16dir grating stimuli; the ultimate goal is to be able to
% identify visually responsive cells (to the stimuli) that we can measure
% df/f (how much fluorescence changed with stimulus - normalized); with
% the df/f data, we will se how df/f varies with orientation/direction ->
% AKA which cells prefer which directions -> we can then ask how these
% cells change/maintain preference over multiple days?
%% PART 1: Locating and loading raw data files to be used in analysis %%

%clear workspace
clear all
clear global
close all
clc
%% get path names
date = '250106';
ImgFolder = strvcat('002'); 
time = strvcat('1530');
mouse = 'i1405';
run = strvcat('002'); 
nrun = size(ImgFolder,1);
frame_rate = 15;
run_str = catRunName(ImgFolder, nrun);
ref_str = catRunName(run, size(run,1)); 
doNewD1 = 0;
tj_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Data\2P_images';
%tj_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\2P_Imaging';
fnout = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P',[date '_' mouse], [date '_' mouse '_' run_str]);
if doNewD1
    fnout = fullfile(fnout,'newD1');
end
% fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\Day1_recycled';
%fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\reverse_match';
behav_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
%% load data
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun 
    %CD = fullfile(tj_fn, [mouse '\' date '_' mouse '\' ImgFolder(irun,:)]); %change current dierectory;
    CD = fullfile(tj_fn, [mouse '\' date '\' ImgFolder(irun,:)]); %change current dierectory;
    cd(CD); %set CD
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat']; %add the 0s to the imaging file
    load(imgMatFile); %**load this file from CD
%     fName = fullfile(behav_fn, ['data-i' '''' mouse '''' '-' date '-' time(irun,:) '.mat']); %find behavior data
    fName = fullfile(behav_fn, ['data-'  mouse  '-' date '-' time(irun,:) '.mat']); %find behavior data
    load(fName); %load behavior data
    
    %input is behavioral parameters
    %info is imaging parameters

    nframes = info.config.frames; %find nframes from info
    fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n']) %graphic display of frames reading
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes); %reads data from 0 to nframes from raw file
    %nPMT x nYpix x nXpix x nframes
    
    
    temp(irun) = input; 
    if isfield(input, 'nScansOn') %checks that nScansOn is in input var
        nOn = temp(irun).nScansOn; %find nOn in input
        nOff = temp(irun).nScansOff; %find nOff in input
        ntrials = size(temp(irun).tGratingDirectionDeg,2); %find ntrials based on input

        data_temp = squeeze(data_temp); %use only 1 PMT so squeeze data
        if nframes>ntrials*(nOn+nOff) %make sure nframes matches
            nframes = ntrials*(nOn+nOff); %if not make it match
            data_temp = data_temp(:,:,1:ntrials*(nOn+nOff)); %will restructure it to nYpix x nXpix x nframes
        elseif nframes<ntrials*(nOn+nOff) %similar to above
            temp(irun) = trialChopper(temp(irun),1:ceil(nframes./(nOn+nOff))); %rounds up to number of frames
        end
    end
    
    offset = offset+nframes; 

    data_temp = squeeze(data_temp); 
    data = cat(3,data,data_temp); %concatenate data and data_temp along 3rd dimension
    trial_n = [trial_n nframes]; 
end
input = concatenateDataBlocks(temp); %combines mat files; 
clear data_temp
%clear temp

%% PART 2: Registering data %%
%We want to remove X-Y artifacts by averaging several hundred frames
%together to find a suitable target image; then all frames will be shifted
%accordingly to match that target

 %% Choose register interval
t = 2000; %nframes to skip for each average; could add nframes to not hard code number to average
nep = floor(size(data,3)./t); %divides frames by skips and rounds down
[n n2] = subplotn(nep); %finds ideal number of subplots to make
figure; %makes figure
for i = 1:nep; %for the number of plots
    subplot(n,n2,i); %this subplot
    imagesc(mean(data(:,:,1+((i-1)*t):100+((i-1)*t)),3)); %scaled color image of mean of frames for specified range
    title([num2str(1+((i-1)*t)) '-' num2str(100+((i-1)*t))]); %titled based on frame numbers
end
%these figures are taking averages of each pixel value across certain sets
%of frames; ex: what is the avg pixel value of pixel 1,1 for these 500
%frames; what about pixel 1,2 etc.
%% Register data

data_avg = mean(data(:,:,12001:12100),3); %mean of pixel values over selected range of frames

if exist(fullfile(fnout)) %if this folder exists)
    load(fullfile(fnout, [date '_' mouse '_' run_str '_reg_shifts.mat'])) %load this mat file
    save(fullfile(fnout, [date '_' mouse '_' run_str '_input.mat']), 'input') %save input?
    [outs, data_reg]=stackRegister_MA(data(:,:,:),[],[],out); %using shifts data to move all frames to target image
%    clear out outs
else
    [out, data_reg] = stackRegister(data,data_avg); %stacks 3d frames data to 2d avg target
    data_reg_avg = mean(data_reg,3); %mean of all registered frames
    reg = data_reg_avg; %sets reg = to above
    mkdir(fullfile(fnout)) %make new directory and save
    save(fullfile(fnout, [date '_' mouse '_' run_str '_reg_shifts.mat']), 'data_reg_avg', 'out', 'data_avg')
    save(fullfile(fnout, [date '_' mouse '_' run_str '_input.mat']), 'input')
end
%clear data


%% test stability
figure; imagesq(data_reg_avg); truesize; % avg pixel value of all frames registered***
print(fullfile(fnout, [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf','-bestfit') %save as pdf that fits the page

%% PART 3: Find activated cells %%
%we want trial level data in df/f; f is the baseline fluorescence (noise)
%and df is the stimulus-specific fluorescence with baseline removed;
%dividing will standardize;
%we will create a cell mask to help with this process
%goal is to get cells that are significantly responsive to stimuli


%calculating trial level data, f, df, and df/f values
nOn = input.nScansOn; %same as above
nOff = input.nScansOff; %same as above
if isfield (info,'frame')
    [stimOns stimOffs] = photoFrameFinder_Sanworks(info.frame);
    nTrials = length(stimOns);
    sz = size(data_reg);
    data_f = nan(sz(1),sz(2),nTrials);
    data_tr = nan(sz(1),sz(2),nTrials);
    for itrial = 1:nTrials
        if ~isnan(stimOns(itrial)) & (stimOns(itrial)+nOn)<sz(3)
            data_f(:,:,itrial) = mean(data_reg(:,:,stimOns(itrial)-nOff/2:stimOns(itrial)-1),3);
            data_tr(:,:,itrial) = mean(data_reg(:,:,stimOns(itrial)+5:stimOns(itrial)+nOn),3);
        end
    end
    data_dfof = (data_tr-data_f)./data_f;
    clear data_tr data_f
else
    report('No Photodiode Data')
    nTrials = length(input.counterValues);
    vals = [];
    times = [];
    for i = 1:nTrials
        vals = [vals input.counterValues{i}];
        times = [times input.counterTimesUs{i}];
    end
    ind = find(vals == 1,1,'last');
    figure
    subplot(2,1,1)
    plot(diff(vals(ind:end)))
    subplot(2,1,2)
    plot(diff(times(ind:end)))
    stimOns = nOff+1:nOff+nOn:nTrials.*(nOff+nOn);
end

if input.doDirStim 
%     obtaining avg and max dfof images for each dir
    Dir = cell2mat_padded(input.tGratingDirectionDeg); %transforms cell array into matrix (1 x ntrials) with dir for each trial
    Dir = Dir(1:ntrials); 
    Dirs = unique(Dir); %what are all the dirs possible
    data_dfof_avg = zeros(sz(1),sz(2),length(Dirs));
    nDirs = length(Dirs); %how mmay directions
    [n n2] = subplotn(nDirs); %find optimal subplot number for ndirs
    figure;
    for idir = 1:length(Dirs) %for each direction
        ind = find(Dir(1:ntrials-1) == Dirs(idir));
        data_dfof_avg(:,:,idir) = mean(data_dfof(:,:,ind),3); %avg across ON frames and trials for each dir
        subplot(n,n2,idir)
        imagesc(data_dfof_avg(:,:,idir)) %plot image avg for each dir
    end
    %clear data_dfof
    myfilter = fspecial('gaussian',[20 20], 0.5); %making filter
    data_dfof_avg_all = imfilter(data_dfof_avg,myfilter); %applying filter to data; filters each pixel for clarity
    data_dfof_max = max(data_dfof_avg_all,[],3); %max of each pixel (within certain direction - not all frames?)
    
%     average dfof image for each direction
    figure; 
    Stims = Dirs; 
    nStim = length(Dirs);
    [n n2] = subplotn(nDirs);
    data_dfof_avg_ori = zeros(sz(1), sz(2), nDirs/2);
    for i = 1:nStim  %for each dir
        subplot(n,n2,i);  %find subplot
        imagesc(data_dfof_avg_all(:,:,i)); %make image of filtered image of particular dir
        clim([0 max(data_dfof_avg_all(:))]) %sets limit to max of dfofavg
        title(num2str(Dirs(i))) %titles each
        colormap(gray) %grayscale
        if i<(nDirs/2)+1
            data_dfof_avg_ori(:,:,i) = mean(data_dfof_avg_all(:,:,[i i+nDirs/2]),3); 
        end
    end
    print(fullfile(fnout, [date '_' mouse '_' run_str '_allStim.pdf']), '-dpdf')

%     average dfof image for each orientation
    if find(Dirs>180)
        figure;
        [n n2] = subplotn(nDirs/2);
        for i = 1:nStim/2
            subplot(n,n2,i)
            imagesc(data_dfof_avg_ori(:,:,i));
            clim([0 max(data_dfof_avg_ori(:))])
            title(num2str(Dirs(i)))
            axis off %this is plotting each ori
        end
        subplot(n,n2,i+1)
        imagesc(max(data_dfof_avg_ori,[],3)) %show max and add to subplot
        title('dfof Max')
        axis off
        data_dfof = cat(3,data_dfof_avg_ori,max(data_dfof_avg_ori,[],3)); %concatenate avg ori and max along 3rd dimension
        print(fullfile(fnout, [date '_' mouse '_' run_str '_activeCells.pdf']), '-dpdf')
    
        figure;
        imagesc(max(data_dfof_avg_ori,[],3)) %plot max only
        title('dfof Max')
        axis off
        print(fullfile(fnout, [date '_' mouse '_' run_str '_dfofMax.pdf']), '-dpdf')
    end
end

 save(fullfile(fnout, [date '_' mouse '_' run_str '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg_all')


%% cell segmentation - helping to identify cells at each ori and max
mask_exp = zeros(sz(1),sz(2)); 
mask_all = zeros(sz(1),sz(2)); 
mask_data = data_dfof_avg_all;

% click on cells in each direction's avg dfof image
for iStim = 1:size(data_dfof_avg_all,3) %for each dir
  mask_data_temp = mask_data(:,:,iStim);
  mask_data_temp(find(mask_exp >= 1)) = 0; %black out old cells
  bwout = imCellEditInteractiveLG(mask_data_temp); %selection GUI
  mask_all = mask_all+bwout; %add new cells to old
  mask_exp = imCellBuffer(mask_all,3)+mask_all; %creates buffer around cells to reduce fusing 2 together
  close all
end
mask_cell = bwlabel(mask_all); %turn logical into numbered cells -> rather than saying there is a cell here, what number is it
figure; imagesc(mask_cell) %colored image of each cell as different color


figure; 
[n n2] = subplotn(nStim);
for i = 1:nStim; 
    subplot(n,n2,i); 
    shade_img = imShade(data_dfof_avg_all(:,:,i), mask_all); %shade grayscale image with mask
    imagesc(shade_img)
    if input.doSizeStim 
    title([num2str(szs(i)) ' deg'])
    elseif input.doRetStim %?
        title([num2str(Stims(i,:))])
    end
    clim([0 max(data_dfof_avg_all(:))])
    colormap(gray)
end
print(fullfile(fnout, [date '_' mouse '_' run_str '_FOV_overlay.pdf']), '-dpdf')

%Create neuropil masks (these are regions around each cell without overlap 
% from neighboring cells); 'noise' from axons/dendrites
mask_np = imCellNeuropil(mask_cell, 3, 5); %creates neuropil mask to remove biological noise from surrounding cells
save(fullfile(fnout, [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_np')

% creating pixel correlation image from a smaller data reg stack -> ?***
data_reg_3hz = stackGroupProject(data_reg,5); %averaging every 5 frames for less noisy stack
pix = getPixelCorrelationImage(data_reg_3hz); %how much each pixel is correlated with surrounding?
pix(isnan(pix))=0;
save(fullfile(fnout, [date '_' mouse '_' run_str '_pixel.mat']),'pix')

%clear data_reg_3hz data_dfof data_dfof_avg max_dfof mask_data mask_all mask_2 data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 
%% neuropil subtraction
%data_tc will be a time course following each cell 
data_tc = stackGetTimeCourses(data_reg, mask_cell); %pixel intensities of stack based on mask
data_tc_down = stackGetTimeCourses(stackGroupProject(data_reg,5), mask_cell); %downsampled version of above
nCells = size(data_tc,2); %number of cells
sz = size(data_reg);
down = 5;
data_reg_down  = stackGroupProject(data_reg,down); %downsample of data reg
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells %for each cell
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i)); %get time course of NP mask in  full and downsample
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf([' Cell #' num2str(i) '%s/n']) 
end

%get weights by maximizing skew
% Find best neuropil weights by maximizing skew on downsampled subtractions
% 
% Assumes calcium signals are 1) sparse and 2) positive.  Too little subtraction 
% will decrease sparseness and too much will make signals negative. Thus, the 
% best neuropil subtraction should yield the highest skew.
ii= 0.01:0.01:1; %potential weights for skew
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i))); 
end
[max_skew ind] =  max(x,[],1); %find max skew and value
np_w = 0.01*ind; %convert to decimal
npSub_tc = data_tc(:,:)-bsxfun(@times,tcRemoveDC(np_tc(:,:)),np_w); %? -> subtract np weighted tc data from tc data



save(fullfile(fnout, [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
save(fullfile(fnout, [date '_' mouse '_' run_str '_input.mat']), 'input')

%% PART 4 - Transform TC into trial structure that is also in dF/F
%We have an np-subtracted TC of cells across frames; we now want to put data into a matrix 
%that is nframes/trial x ncells x ntrials so we can see which cells were
%responding to which stimuli; we also need F to be dF/F
ntrials = size(input.tGratingDirectionDeg,2); %these lines are same as above
Dir = cell2mat_padded(input.tGratingDirectionDeg);
Dir = Dir(1:ntrials);
    Dirs = unique(Dir);
    nDirs = length(Dirs);
if isfield(input, 'nScansOn')
    nOn = input.nScansOn;
    nOff = input.nScansOff;
    nCells = size(npSub_tc,2);
    data_mat = nan(nOn+nOff, nCells, ntrials);
    for itrial = 1:ntrials
        if stimOns(itrial)+nOn+nOff/2-1 <= size(npSub_tc,1)
            data_mat(:,:,itrial) = npSub_tc(stimOns(itrial)-nOff/2:stimOns(itrial)+nOn+nOff/2-1,:); 
        end
    end
    data_f = mean(data_mat(1:nOff/2,:,:),1);
    data_df = bsxfun(@minus, data_mat, data_f);
    data_dfof = bsxfun(@rdivide, data_df, data_f); %dfof in nframes x ncells x ntrials
    
    ndir = length(Dirs);
    [n, n2] = subplotn(nCells);
    h_dir = zeros(nCells, ndir);
    p_dir = zeros(nCells, ndir);
    base_win = 1:nOff/2;
    resp_win = nOff/2+5:nOff/2+nOn;
    base = squeeze(mean(data_dfof(base_win,:,:),1)); %averaging across baseline window
    resp = squeeze(mean(data_dfof(resp_win,:,:),1)); %averaging across response window
    resp_mat = resp-base;
    dir_resp = zeros(nCells,ndir);
    [x y] = ttest(resp', base', 'tail','right'); %ttest comparing base to resp for significance
    no_match = find(isnan(x)); %?***
    max_dir = zeros(nCells,1);
    figure;
    for i = 1:nCells
        if ~sum(no_match==i) 
        subplot(n, n2, i)
            for idir = 1:ndir %for each dir
                if nOn>29
                    ind = find(Dir == Dirs(idir)); %find trials of that dir
                else
                    ind = find(Dir(1:ntrials-1) == Dirs(idir));
                end
                [h_dir(i,idir), p_dir(i,idir)] = ttest(resp(i,ind), base(i,ind),'tail','right','alpha', 0.05/(ndir-1)); %ttest of each cell at each dir for base vs. resp
                if h_dir(i,idir) %if this cell/dir is sig make red, if not make black and plot
                    errorbar(Dirs(idir), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'or')
                else
                    errorbar(Dirs(idir), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'ok')
                end
                dir_resp(i,idir) = mean(resp(i,ind)-base(i,ind),2); %avg response at each dir for each cell
                hold on
            end
            if sum(h_dir(i,:),2)>0 %if cell has one sig dir
                temp_resp = dir_resp(i,:);
                temp_resp(find(h_dir(i,:)==0)) = NaN;
                [max_val max_ind] = max(temp_resp,[],2); %find max index and value
                max_dir(i,:) = max_ind; %which index/dir was highest?
            else
                [max_val max_ind] = max(dir_resp(i,:),[],2);
                max_dir(i,:) = max_ind;
            end
            title([num2str(Dirs(max_dir(i,:))) ' deg'])
        end
    end
        print(fullfile(fnout, [date '_' mouse '_' run_str '_dirTuning.pdf']),'-dpdf')
    %this section is similar to above but for ori rather than dir
    if find(Dirs>180)
    nori = length(Dirs)/2;
    Ori = Dir;
    for iori = 1:nori
        ind = find(Dir == Dirs(iori+nori));
        Ori(ind) = Dirs(iori);
    end
    Oris = unique(Ori);
    h_ori = zeros(nCells, nori);
    p_ori = zeros(nCells, nori);
    ori_resp = zeros(nCells,nori);
    max_ori = zeros(nCells,1);
    figure;
    for i = 1:nCells
        if ~sum(no_match==i)
            subplot(n, n2, i)
            for iori = 1:nori
                if nOn>29
                    ind = find(Ori == Oris(iori));
                else
                    ind = find(Ori(1:ntrials-1) == Oris(iori));
                end
                [h_ori(i,iori), p_ori(i,iori)] = ttest(resp(i,ind), base(i,ind),'tail','right','alpha', 0.05/(nori-1));
                if h_ori(i,iori)
                    errorbar(Oris(iori), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'or')
                else 
                    errorbar(Oris(iori), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'ok')
                end
                ori_resp(i,iori) = mean(resp(i,ind)-base(i,ind),2);
                hold on
            end
            if sum(h_ori(i,:),2)>0
                temp_resp = ori_resp(i,:);
                temp_resp(find(h_ori(i,:)==0)) = NaN;
                [max_val max_ind] = max(temp_resp,[],2);
                max_ori(i,:) = max_ind;
            else
                [max_val, max_ind] = max(ori_resp(i,:),[],2);
                max_ori(i,:) = max_ind;
            end
            title([num2str(Oris(max_ori(i,:))) ' deg'])
        end
    end
    
    good_ind = unique([find(x)'; find(sum(h_dir,2)>0); find(sum(h_ori,2)>0)]); %index of responsive cells
    print(fullfile(fnout, [date '_' mouse '_' run_str '_oriTuning.pdf']),'-dpdf')
    else
        Ori = Dir;
        nOri = ndir;
        Oris = Dirs;
        h_ori = h_dir;
        max_ori = max_dir;
        good_ind = unique([find(x)'; find(sum(h_dir,2)>0); find(sum(h_ori,2)>0)]); %index of responsive cells
    end
    save(fullfile(fnout, [date '_' mouse '_' run_str '_trialData.mat']),'data_dfof','max_dir','h_dir', 'h_ori', 'max_ori','good_ind','dir_resp','resp_win','base_win','resp','base','resp_mat')
end

%% ori fitting - getting a von Mises fit over our raw data to calculate theoretical pref ori
tOrientation = Ori;
tOrientation(tOrientation > 179) = tOrientation(tOrientation > 179) - 180;
[orientationInd, orientations] = findgroups(tOrientation);
nStim = length(orientations);

nBoot = 1000;
avgResponseEaOri = nan(nCells,nStim);
semResponseEaOri = nan(nCells,nStim);
tuningResamp = nan(nCells,nStim,nBoot);

for istim = 1:nStim
    ind = find(orientationInd == istim);
    avgResponseEaOri(:,istim) = mean(resp_mat(:,ind),2,'omitnan');
    semResponseEaOri(:,istim) = std(resp_mat(:,ind),[],2,'omitnan')./sqrt(length(ind));
    for iboot = 1:nBoot
        n = length(ind);
        randTrials = randsample(ind,n,1);
        tuningResamp(:,istim,iboot) = mean(resp_mat(:,randTrials),2,'omitnan');
    end
end
[vonMisesFitAllCellsAllBoots,~,fitReliability,R_square] = vonmisesReliableFit(avgResponseEaOri,...
    tuningResamp,orientations,nBoot);


save(fullfile(fnout, [date '_' mouse '_' run_str '_oriTuningAndFits.mat']),...
            'avgResponseEaOri','semResponseEaOri','vonMisesFitAllCellsAllBoots','fitReliability','R_square')

%% plotting the ori tuning curves - the von Mises fit is overlaying the raw data
dir_mat = celleqel2mat_padded(input.tGratingDirectionDeg);
ori_mat = dir_mat;
ori_mat(dir_mat>=180) = ori_mat(dir_mat>=180)-180;
oris = unique(ori_mat);
figure; 
if nCells<49 
    [n n2] = subplotn(nCells);
else
    [n, n2] = subplotn(49);
end
start = 1;
x = 0;
for ic = 1:nCells
    if start > 49
        suptitle([mouse ' ' date ' n = ' num2str(length(find(fitReliability<22.5))) '/' num2str(nCells) '- well-fit'])
        print(fullfile(fnout, [date '_' mouse '_' run_str '_oriTuningFits_cells' num2str(start-49) '-' num2str(start-1) '_newD1.pdf']),'-dpdf','-fillpage')
        start = 1;
        x = x+1;
        figure;
    end
    subplot(n,n2,ic-(x.*49))
    errorbar(oris,avgResponseEaOri(ic,:), semResponseEaOri(ic,:),'-o')
    hold on
    plot(0:180,vonMisesFitAllCellsAllBoots(:,1,ic));
    tit_str = num2str(chop(R_square(1,ic),2));
    if fitReliability(ic)<22.5
        tit_str = [tit_str '- R'];
    end
    title(tit_str)
    start = start+1;
end
suptitle([mouse ' ' date ' n = ' num2str(length(find(fitReliability<22.5))) '/' num2str(nCells) '- well-fit'])
print(fullfile(fnout, [date '_' mouse '_' run_str '_oriTuningFits' num2str(start-49) '-' num2str(start-1) '.pdf']),'-dpdf','-fillpage')

[max_resp prefOri] = max(vonMisesFitAllCellsAllBoots,[],1);
prefOri = squeeze(prefOri)-1;
prefOri_bootdiff = abs(prefOri(2:end,:)-prefOri(1,:));
prefOri_bootdiff(find(prefOri_bootdiff>90)) = 180-prefOri_bootdiff(find(prefOri_bootdiff>90));
ind_theta90 = find(prctile(prefOri_bootdiff,90,1)<22.5); %this is an index of cells w/ reliable fit over bootstraps
edges = [0 22.5:45:180]; 
[bin ind_bin] = histc(prefOri(1,:),edges); %how many pref oris are between 0-22.5? 22.5-67.5? etc.* - this is not really used
ind_bin(find(ind_bin==5)) = 1;
bin(1) = bin(1)+bin(5);
bin(5) = [];

tunedCells = cell(1,length(bin)); 
for i = 1:length(bin)
    tunedCells{i} = intersect(find(ind_bin==i),ind_theta90);
end

save(fullfile(fnout, [date '_' mouse '_' run_str '_oriTuningInfo.mat']),...
    'prefOri', 'prefOri_bootdiff', 'ind_theta90', 'tunedCells');

%% identify cells that are red - skip this if only imaging green cells
%If you imaged w/ 2 colors, you should have a snapshot for green, a snapshot for red, and the full
%run for green (which was analyzed above)

% rename variables from green imaging 
fov_avg{1} = data_reg_avg;
dfmax{1} = data_dfof_max;
corrmap{1} = pix;
masks{1} = mask_cell;
corrmap_norm{1} = uint8((corrmap{1}./max(corrmap{1}(:))).*255);
brightnessScaleFactor = 0.5;
fov_norm{1} = uint8((fov_avg{1}./max(fov_avg{1}(:))).*255);
fov_norm{1}(fov_norm{1} > (brightnessScaleFactor*255)) = brightnessScaleFactor*255;

% process the red channel from a 1000 frame run

%1040 run - make sure to set the right image folder
irun = 1;
WL = '1040'; 
ImgFolder = strvcat('002'); %
run = catRunName(ImgFolder, nrun);
imgMatFile = [ImgFolder '_000_000.mat'];
CD = fullfile(tj_fn, [mouse '\' date '_' mouse '\' ImgFolder(irun,:)]);
cd(CD);
load(imgMatFile);
data_temp = sbxread(imgMatFile(1,1:11),0,info.config.frames); 
mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run]));

%below is the best pipeline to register the red imaging snapshot with the green from the full run
%for maximal red/green cell matching
if size(data_temp,1)>1
data_1040_green = squeeze(data_temp(1,:,:,:)); %PMT 0 (green)
data_1040_red = squeeze(data_temp(2,:,:,:)); %PMT 1 (red)
[out_1040_red_regtoself data_1040_red_regtoself] = stackRegister(data_1040_red,mean(data_1040_red,3)); %register 1040 red channel to self
[out_1040_green_regtomaingreen data_1040_green_regtomaingreen] = stackRegister(data_1040_green,data_reg_avg); %register 1040 green to master green
[out_1040_red_regto1040green data_red_regto1040green] = stackRegister_MA(mean(data_1040_red_regtoself,3),[],[],out_1040_green_regtomaingreen); %red registered to 1040 green
data_red_regto1040green_avg = mean(data_red_regto1040green,3);
red = data_red_regto1040green_avg;
fov_red{1} = uint8((red./max(red(:))).*255);
end



%% select red cells - this will produce a new index of which cells were red (and green)
% size of cell box
clear input
close all
w=30;
h=30;
buf = 3;

% get cell centroids
redGreenCells = struct;

cellPosition = regionprops(masks{1});
nc = length(cellPosition);

xCenter = cellfun(@(a) round(a(1)),{cellPosition.Centroid});
yCenter = cellfun(@(a) round(a(2)),{cellPosition.Centroid});

% index cells NOT too close to edge and NOT in black part of transformation
[ypix,xpix] = size(fov_avg{1});
goodCells = xCenter>(w/2) & xCenter<xpix-(w/2) & yCenter>(h/2) & yCenter<ypix-(h/2);

goodCells = goodCells & ...
    arrayfun(@(x) sum(sum(fov_norm{1}(masks{1}==x)))>0,1:nc);

    % fine register each cell    
mask_exp = zeros(size(fov_avg{1}));
mask_all = zeros(size(fov_avg{1}));

for icell = 1:nc
    if goodCells(icell)
        % find best shift
        day1_cell_avg = fov_norm{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        
        corr = corrmap_norm{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        
        day1_cell_max = dfmax{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        
        day1_red_avg = fov_red{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        
        mask = masks{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);

            pass = true;
            figure;
            movegui('center')
            start = 1;
            subplot(1,4,start)
            imagesc(day1_cell_avg);axis image
            hold on
            bound = cell2mat(bwboundaries(mask(:,:,1)));
            plot(bound(:,2),bound(:,1),'-','color','r','MarkerSize',2);axis image
            title('avg')
            subplot(1,4,start+1)
            imagesc(corr);axis image
            hold on
            plot(bound(:,2),bound(:,1),'-','color','r','MarkerSize',2);axis image
            title('corr')
            subplot(1,4,start+2)
            imagesc(day1_cell_max);axis image
            hold on
            plot(bound(:,2),bound(:,1),'-','color','r','MarkerSize',2);
            title('dfof')
            subplot(1,4,start+3)
            imagesc(day1_red_avg);axis image
            hold on
            plot(bound(:,2),bound(:,1),'-','color','r','MarkerSize',2);
            title('Red')
            
            prompt = 'Choose one: 1- good, 2-okay, 3- none: ';
            drawnow
            x = input(prompt);
            switch x
                case 1
                    pass = true;
                    faint = false;
                case 2
                    pass = false;
                    faint = true;
                case 3
                    pass = false;
                    faint = false;
            end
    if pass
        redGreenCells(icell).pass = pass;
        redGreenCells(icell).faint = false;  
    elseif faint
        redGreenCells(icell).pass = false;
        redGreenCells(icell).faint = faint;  
    else 
        redGreenCells(icell).pass = false;
        redGreenCells(icell).faint = false;  
    end
    close all
    end  
end

goodCells = find([redGreenCells.pass]);
okayCells = find([redGreenCells.faint]);
redCells = sort([goodCells okayCells]);
redChImg = fov_red{1};
save(fullfile(fnout, [date '_' mouse '_' run_str '_multiday_alignment.mat']),'goodCells','okayCells','redCells','redChImg');
