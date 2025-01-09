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
ref_date = '250106';
ref_str = 'runs-002';
date = '250108';
ImgFolder = strvcat('003'); 
time = strvcat('1641');
mouse = 'i1405';
reg_run = strvcat('002'); 
nrun = size(ImgFolder,1);
frame_rate = 15;
run_str = catRunName(ImgFolder, nrun);
reg_str = catRunName(reg_run, size(reg_run,1)); 
doNewD1 = 0;
tj_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Data\2P_images';
%tj_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\2P_Imaging';
fnout = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P',[date '_' mouse], [date '_' mouse '_' run_str]);
fnout_reg = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P',[date '_' mouse], [date '_' mouse '_' reg_str]);
fnout_ref = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P',[ref_date '_' mouse], [ref_date '_' mouse '_' ref_str]);
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

%% Register data

if exist(fullfile(fnout)) %if this folder exists)
    load(fullfile(fnout, [date '_' mouse '_' run_str '_reg_shifts.mat'])) %load this mat file
    save(fullfile(fnout, [date '_' mouse '_' run_str '_input.mat']), 'input') %save input?
    [outs, data_reg]=stackRegister_MA(data(:,:,:),[],[],out); %using shifts data to move all frames to target image
%    clear out outs
else
    load(fullfile(fnout_reg, [date '_' mouse '_' reg_str '_reg_shifts.mat']), 'data_avg')
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

%% align, get masks, tcs, and neuropil subtraction
%data_tc will be a time course following each cell 
%load rotation from ref
load(fullfile(fnout_reg, [date '_' mouse '_' reg_str '_rematch_shifts.mat']))
data3 = imwarp(data_reg,fitGeoTAf, 'OutputView', imref2d(size(data_reg))); %displacing d2 data reg by the fits above
clear data_reg

% load masks from reg
load(fullfile(fnout_reg, [date '_' mouse '_' reg_str '_mask_cell.mat']))
TCs_D1 = load(fullfile(fnout_ref, [ref_date '_' mouse '_' ref_str '_TCs.mat']));
cellTCs_all{1} = TCs_D1.npSub_tc; %np-subtracted TC from d1
load(fullfile(fnout_reg, [date '_' mouse '_' reg_str '_multiday_alignment.mat']))
match = find([cellImageAlign.pass]); %finds matched cell indices
match_ind = match;
cellTCs_match{1} = cellTCs_all{1}(:,match_ind);

data_tc = stackGetTimeCourses(data3, mask_cell); %pixel intensities of stack based on mask
data_tc_down = stackGetTimeCourses(stackGroupProject(data3,5), mask_cell); %downsampled version of above
nCells = size(data_tc,2); %number of cells
sz = size(data3);
down = 5;
data_reg_down  = stackGroupProject(data3,down); %downsample of data reg
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells %for each cell
     np_tc(:,i) = stackGetTimeCourses(data3,mask_np(:,:,i)); %get time course of NP mask in  full and downsample
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
cellTCs_match{2} = npSub_tc(:,:);


save(fullfile(fnout, [date '_' mouse '_' run_str '_TCs.mat']),'cellTCs_match', 'cellTCs_all')
clear data3

%% PART 4 - Transform TC into trial structure that is also in dF/F
%We have an np-subtracted TC of cells across frames; we now want to put data into a matrix 
%that is nframes/trial x ncells x ntrials so we can see which cells were
%responding to which stimuli; we also need F to be dF/F
ntrials = size(input.tGratingDirectionDeg,2); %these lines are same as above
Dir = cell2mat_padded(input.tGratingDirectionDeg);
Dir = Dir(1:ntrials);
Dirs = unique(Dir);
nDirs = length(Dirs);
[stimOns stimOffs] = photoFrameFinder_Sanworks(info.frame);
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
end
base_win = 1:nOff/2;
resp_win = nOff/2+5:nOff/2+nOn;
[h_resp p_resp] = ttest(mean(data_dfof(resp_win,:,:),1),mean(data_dfof(base_win,:,:),1),'Dim',3,'tail','right');
data_resp = squeeze(mean(data_dfof(resp_win,:,:)));
save(fullfile(fnout, [date '_' mouse '_' run_str '_respData.mat']), 'data_dfof', 'data_resp', 'base_win', 'resp_win', 'h_resp', 'p_resp')
