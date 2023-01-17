%% Load, register, segment and neuropil correct 2P data
close all
clear all global
clc

ds = 'FlashStim_ExptList';
eval(ds);

iexp  = 4;

mouse = expt(iexp).mouse;
date = expt(iexp).date;
ImgFolder = expt(iexp).postFolder;
time = cell2mat(expt(iexp).postTime);
nrun = length(ImgFolder);
run_str = catRunName(ImgFolder, nrun);
pre_run_str = catRunName(expt(iexp).preFolder,length(expt(iexp).preFolder));
ImgFolder = cell2mat(ImgFolder);

%Path names
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';

if strcmp(expt(iexp).saveLoc,'tj')
    data_fn = fullfile(fn_base, 'home\tj\2p_Imaging\');
end
lg_fn = fullfile(fn_base, 'home\lindsey');
mworks_fn = fullfile(fn_base, 'Behavior\Data');
fnout = fullfile(lg_fn, 'Analysis\2P');

frame_rate = params.frameRate;
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];
datemouseprerun = [date '_' mouse '_' pre_run_str];

%Load 2P data
%Load mworks data- this has information about experiment (e.g. the visual stimuli presented and synchronization with the microscope)
fName = fullfile(mworks_fn, ['data-' mouse '-' date '-' time '.mat']);
load(fName);
%Load 2P metadata- this has information about the information about the imaging session (e.g. frame count, zoom)
CD = fullfile(data_fn, mouse, datemouse, ImgFolder);
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
 
load(fullfile(fnout, datemouse, datemouseprerun, [datemouseprerun '_reg_shifts.mat']))

mkdir(fullfile(fnout, datemouse, datemouserun));
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_input.mat']), 'input')

[out data_reg] = stackRegister(data,data_avg);
data_reg_avg = mean(data_reg,3);
save(fullfile(fnout, datemouse, datemouseprerun, [datemouseprerun '_reg_shifts.mat']),'out','data_avg','data_reg_avg')
%Test registration
%    b. Average of all frames should be sharp
imagesq(data_reg_avg); 
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_FOV_avg.pdf']), '-dpdf')

clear data

%% Neuropil subtraction
load(fullfile(fnout, datemouse, datemouseprerun, [datemouseprerun '_mask_cell.mat']))
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
nOff = input.nScansOff;
nOn = input.nScansOn;
tDir = celleqel2mat_padded(input.tGratingDirectionDeg);
Dirs = unique(tDir);
nDirs = length(Dirs);
ntrials = length(tDir);

base_win = nOff/2:nOff;
resp_win = nOff+5:nOff+nOn;
data_trial = reshape(npSub_tc, [nOn+nOff ntrials nCells]);
data_f = mean(data_trial(base_win,:,:),1);
data_dfof = bsxfun(@rdivide,bsxfun(@minus,data_trial,data_f),data_f);

resp_cell_dir = cell(1,nDirs);
base_cell_dir = cell(1,nDirs);
data_dfof_dir = zeros(nCells,nDirs,2);
h_dir = zeros(nDirs,nCells);
p_dir = zeros(nDirs,nCells);
for iDir = 1:nDirs
    ind = find(tDir==Dirs(iDir));
    resp_cell_dir{iDir} = squeeze(mean(data_dfof(resp_win,ind,:),1));
    base_cell_dir{iDir} = squeeze(mean(data_dfof(base_win,ind,:),1));
    [h_dir(iDir,:), p_dir(iDir,:)] = ttest(resp_cell_dir{iDir},base_cell_dir{iDir},'tail','right','alpha',0.05./(nDirs-1));
    data_dfof_dir(:,iDir,1) = squeeze(mean(mean(data_dfof(resp_win,ind,:),1),2));
    data_dfof_dir(:,iDir,2) = squeeze(std(mean(data_dfof(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
end

h_all_dir = sum(h_dir,1);
resp_ind = find(h_all_dir);

tOri = tDir;
tOri(find(tDir>=180)) = tOri(find(tDir>=180))-180;
Oris = unique(tOri);
nOri = length(Oris);

resp_cell_ori = cell(1,nOri);
base_cell_ori = cell(1,nOri);
data_dfof_ori = zeros(nCells,nOri,2);
h_ori = zeros(nOri,nCells);
p_ori = zeros(nOri,nCells);
for iOri = 1:nOri
    ind = find(tOri==Oris(iOri));
    resp_cell_ori{iOri} = squeeze(mean(data_dfof(resp_win,ind,:),1));
    base_cell_ori{iOri} = squeeze(mean(data_dfof(base_win,ind,:),1));
    [h_ori(iOri,:) p_ori(iOri,:)] = ttest(resp_cell_ori{iOri},base_cell_ori{iOri},'tail','right','alpha',0.05./(nOri-1));
    data_dfof_ori(:,iOri,1) = squeeze(mean(mean(data_dfof(resp_win,ind,:),1),2));
    data_dfof_ori(:,iOri,2) = squeeze(std(mean(data_dfof(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
end

h_all_ori = sum(h_ori,1);
b_ori = zeros(1,nCells);
k1_ori = zeros(1,nCells);
R1_ori = zeros(1,nCells);
u1_ori = zeros(1,nCells);
R_square_ori = zeros(1,nCells);
sse_ori = zeros(1,nCells);
stim_DSI = zeros(1,nCells);
stim_OSI = zeros(1,nCells);

Oris = Dirs(1:nDirs/2);
nOris = length(Oris);
start=1;
n = 1;
figure;
movegui('center')
theta_high = [0:179];
y_fit = zeros(180,nCells);
for iCell = 1:nCells
    if start>25
        print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_cellTuningOri' num2str(n) '.mat']),'-dpdf','-bestfit')
        figure;movegui('center');
        start = 1;
        n = n+1;
    end
    subplot(5,5,start)
    errorbar(Oris, data_dfof_ori(iCell,:,1), data_dfof_ori(iCell,:,2), '-o')
    title(['R = ' num2str(h_all_ori(iCell))])

    data = [data_dfof_ori(iCell,:,1) data_dfof_ori(iCell,1,1)];
    theta = [deg2rad(Oris) pi];
    [b_ori(:,iCell),k1_ori(:,iCell),R1_ori(:,iCell),u1_ori(:,iCell),sse_ori(:,iCell),R_square_ori(:,iCell)] ...
        = miaovonmisesfit_ori(theta,data);
    y_fit(:,iCell) = b_ori(:,iCell)+R1_ori(:,iCell).*exp(k1_ori(:,iCell).*(cos(2.*(deg2rad(theta_high)-u1_ori(:,iCell)))-1));
    hold on 
    plot(theta_high,y_fit(:,iCell))
    [max_val max_ind] = max(data_dfof_ori(iCell,:,1),[],2);
    null_ind = max_ind+(nOris./2);
    null_ind(find(null_ind>nOris)) = null_ind(find(null_ind>nOris))-nOris;
    min_val = data_dfof_ori(iCell,null_ind,1);
    if min_val<0
        min_val = 0;
    end
    stim_OSI(1,iCell) = (max_val-min_val)./(max_val+min_val);
    [max_val max_ind] = max(data_dfof_dir(iCell,:,1),[],2);
    null_ind = max_ind+(nDirs./2);
    null_ind(find(null_ind>nDirs)) = null_ind(find(null_ind>nDirs))-nDirs;
    min_val = data_dfof_dir(iCell,null_ind,1);
    if min_val<0
        min_val = 0;
    end
    stim_DSI(1,iCell) = (max_val-min_val)./(max_val+min_val);
    start = start +1;
end
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_cellTuningOri' num2str(n) '.mat']),'-dpdf','-bestfit')

save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_oriResp.mat']), 'data_dfof_dir', 'data_dfof_ori','base_win','resp_win','h_dir','h_ori','k1_ori','b_ori','R1_ori','u1_ori','stim_DSI','stim_OSI','R_square_ori','y_fit','Oris','Dirs')
