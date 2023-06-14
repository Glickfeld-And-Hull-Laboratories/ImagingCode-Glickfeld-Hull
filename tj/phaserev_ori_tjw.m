%clear everything
clear all
clear all global
clc
close all

%%
tj_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\2P_Imaging';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\SRP_segmented';
% fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P';
behav_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
dataset = 'exp_list_phaserev_tjw'; %experiment list to pick files from
eval(dataset); %load dataset
%%
expt_id = 1;
date = expt(expt_id).date;
ImgFolder = expt(expt_id).runs; %could we use char() instead here?
time = expt(expt_id).time_mat;
mouse = expt(expt_id).mouse;
run = ImgFolder; %multiple depths?***
nrun = size(ImgFolder,1); %what is this?***
run_str = catRunName(ImgFolder, nrun);
%%
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun 
    CD = fullfile(tj_fn, [mouse '\' date '_' mouse '\' ImgFolder(irun,:)]); %identify current dierectory;
    cd(CD); %set CD
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat']; %add the 0s to the imaging file
    load(imgMatFile); %**load this file from CD
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
    data = cat(3,data,data_temp); %concatenate data and data_temp along 3rd dimension;
    trial_n = [trial_n nframes]; 
end
input = concatenateDataBlocks(temp); %combines mat files;
clear data_temp
%clear temp

 %% Choose register interval
t = 2000; %nframes to skip for each average; could add nframes to not hard code number to average
nep = floor(size(data,3)./t); %divides frames by skips and rounds down
[n n2] = subplotn(nep); %finds ideal number of subplots to make
figure; %makes figure
for i = 1:nep; %for the number of plots
    subplot(n,n2,i); %this subplot
    imagesc(mean(data(:,:,1+((i-1)*t):500+((i-1)*t)),3)); %scaled color image of mean of frames for specified range
    title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]); %titled based on frame numbers
end
%these figures are taking averages of each pixel value across certain sets
%of frames; ex: what is the avg pixel value of pixel 1,1 for these 500
%frames; what about pixel 1,2 etc.

%% Register data

data_avg = mean(data(:,:,4001:4500),3); %mean of pixel values over selected range of frames

if exist(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str])) %if this folder exists)
    load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat'])) %load this mat file
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input') %save input?
    [outs, data_reg]=stackRegister_MA(data(:,:,:),[],[],out); %using shifts data to move all frames to target image
%    clear out outs
else
    [out, data_reg] = stackRegister(data,data_avg); %stacks 3d frames data to 2d avg target
    data_reg_avg = mean(data_reg,3); %mean of all registered frames
    reg = data_reg_avg; %sets reg = to above
    mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str])) %make new directory and save
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'data_reg_avg', 'out', 'data_avg')
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
%clear data


%% test stability
figure; imagesq(data_reg_avg); truesize; %why not imagesc?; avg pixel value of all frames registered***
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf','-bestfit') %save as pdf that fits the page

%%
nOn = double(input.nScansOn);
nOff = double(input.nScansOff);
phaseCyc = double(input.nScansPhaseCyc);
dir_mat = celleqel2mat_padded(input.tGratingDirectionDeg);
nTrials = length(dir_mat);
sz = size(data_reg);

data_resp = reshape(data_reg, [sz(1) sz(2) nOn+nOff nTrials]);
data_f = mean(data_resp(:,:,nOff-15:nOff,:),3);
data_resp_dfof = (double(data_resp)-data_f)./data_f; %Ypix x Xpix x frames x trials
data_inv_resp_dfof = (data_f-double(data_resp))./data_f; %Ypix x Xpix x frames x trials

clear data_resp data_f

dirs = unique(dir_mat);
nDir = length(dirs);

nStim = double(nDir); %number of stimulus combination selected randomly from list of 90 w/ replacement

% [n n2] = subplotn(nStim);
% data_dfof_stim = nan(sz(1),sz(2),nDir);
% start = 1;
% for iDir = 1:nDir
%     ind = find(dir_mat == dirs(iDir));
%     data_dfof_stim(:,:,iDir) = mean(mean(data_resp_dfof(:,:,nOff+1:nOff+nOn,ind),4),3);
%     subplot(n,n2,start)
%     imagesc(data_dfof_stim(:,:,iDir))
%     start = start+1;
% end

[n n2] = subplotn(nStim);
data_dfof_stim = nan(sz(1),sz(2),nDir);
start = 1;
for iDir = 1:nDir
    ind = find(dir_mat == dirs(iDir));
%     data_dfof_stim(:,:,iDir,1) = mean(mean(data_resp_dfof(:,:,nOff+1:nOff+input.nScansPhaseCyc,ind),4),3);
    data_dfof_stim(:,:,iDir,1) = mean(mean(data_resp_dfof(:,:,nOff+1:nOn+nOff,ind),3),4);
    data_inv_dfof_stim(:,:,iDir,1) = mean(mean(data_inv_resp_dfof(:,:,nOff+1:nOn+nOff,ind),3),4);
    subplot(n,n2,start)
    imagesc(data_dfof_stim(:,:,iDir,1))
    imagesc(data_inv_dfof_stim(:,:,iDir,1))
    start = start+1;
%     data_dfof_stim(:,:,iDir,2) = mean(mean(data_resp_dfof(:,:,nOff+1+input.nScansPhaseCyc:nOff+input.nScansPhaseCyc.*2,ind),4),3);
%     subplot(n,n2,start)
%     imagesc(data_dfof_stim(:,:,iDir,2))
%     start = start+1;
end


myfilter = fspecial('gaussian',[20 20], 0.5); %making filter

data_dfof_stim_all = reshape(data_dfof_stim, [sz(1) sz(2) nDir]);
data_dfof_stim_all = imfilter(data_dfof_stim_all, myfilter);
data_inv_dfof_stim_all = imfilter(data_inv_dfof_stim, myfilter);
data_dfof_max = max(data_dfof_stim_all,[],3);
clear data_dfof_stim data_resp_dfof data_inv_resp_dfof

figure; imagesc(data_dfof_max); movegui('center')

save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_segmentData.mat']), 'data_dfof_max', 'data_dfof_stim_all', 'nStim')

%%
%pixel correlation
data_reg_3hz = stackGroupProject(data_reg,5); %averaging every 5 frames for less noisy stack
pix = getPixelCorrelationImage(data_reg_3hz); %how much each pixel is correlated with surrounding? -> could save this figure**
pix(isnan(pix))=0;
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pixel.mat']),'pix')

data_dfof = cat(3,data_dfof_stim_all,data_inv_dfof_stim_all,pix);


%%
sz = size(data_dfof);
mask_data = data_dfof;
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
for iStim = 1:size(data_dfof,3)    
    mask_data_temp = mask_data(:,:,iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0; %blacks out old cells
    bwout = imCellEditInteractiveLG(mask_data_temp); %selection GUI
    mask_all = mask_all+bwout; %adds new cells to old cells
    mask_exp = imCellBuffer(mask_all,3)+mask_all; %creates buffer around cells to avoid fusing
    close all
end
mask_cell = bwlabel(mask_all); %turns logical into numbered cells
figure;
imagesc(mask_cell)
nMaskPix = 5; %thickness of neuropil ring in pixels
nBuffPix = 3; %thickness of buffer between cell and ring
mask_np = imCellNeuropil(mask_cell,nBuffPix,nMaskPix);

save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'mask_cell', 'mask_np')

clear data_dfof mask_data mask_all data_f data_dfof_max



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


clear data data_reg data_reg_down data_tc_down np_tc_down
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

%% 

nCells = size(npSub_tc,2);
cycPerTrial = floor(nOn/(phaseCyc*2)); %how many complete phase transitions happened during the nOn period - 2 in this case
data_trial = permute(reshape(npSub_tc,[nOn+nOff nTrials nCells]),[1 3 2]);

%%
data_f = mean(data_trial(nOff-15:nOff,:,:),1);
data_dfof = (data_trial-data_f)./data_f;

data_dfof_dir = zeros(nOn+nOff, nCells, nDir);

for iDir = 1:nDir
    ind_dir = find(dir_mat == dirs(iDir));
    data_dfof_dir(:,:,iDir) = squeeze(mean(data_dfof(:,:,ind_dir),3));
end

%%

ntrials = size(input.tGratingDirectionDeg,2); %these lines are same as above
Dir = cell2mat_padded(input.tGratingDirectionDeg);
Dir = Dir(1:ntrials);
    Dirs = unique(Dir);
    nDirs = length(Dirs);

%%
base_win = nOff-15:nOff;
resp_win = nOff+1:nOn;
base = squeeze(mean(data_dfof(base_win,:,:),1)); %avg across base
resp = squeeze(mean(data_dfof(resp_win,:,:),1));

ptrials = nTrials-1;
psig = double(0.05/ptrials);


for i = 1:nCells

                [h(i,:), p(i,:)] = ttest(resp(i,:), base(i,:),'tail','right');
                
end

sig_cell = zeros(1,nCells);

for i = 1:nCells
    sig_cell(i) = sum(h(i,:)); %sum of all direction significance levels (0 or 1) for each cell
end

n_sig_cell = sum(sig_cell>=1);



%%
%significantly responsive cells only
sig_index = find(sig_cell == 1);

data_dfof = data_dfof(:,sig_index,:);

avg_window_resp = squeeze(mean(data_dfof(resp_win,:,:),1));

avg_cell_resp = mean(mean(data_dfof(resp_win,:,:),1),3);
std_cell_resp = std(avg_window_resp,[],2);
se_cell_resp = std_cell_resp/sqrt(size(avg_window_resp,2));

avg_pop_resp = mean(avg_window_resp,1);
std_pop_resp = std(avg_window_resp,[],1);
se_pop_resp = std_pop_resp/sqrt(size(avg_window_resp,1));

figure; plot(mean(mean(data_dfof,2),3))

figure;
subplot(1,2,1);
errorbar(avg_cell_resp, se_cell_resp);
subplot(1,2,2);
errorbar(avg_pop_resp, se_pop_resp);

%%
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_keyvars.mat']), 'avg_cell_resp', 'se_cell_resp', ...
    'avg_pop_resp', 'se_pop_resp', 'sig_cell', 'sig_index', 'n_sig_cell', 'data_dfof')


%%

figure;
for i = 1:25
    subplot(5,5,i);
    plot(mean(data_dfof(:,i,:),3))
end

%%
% %%
% 
%  
%     ndir = length(Dirs);
%     [n, n2] = subplotn(nCells);
%     h_dir = zeros(nCells, ndir);
%     p_dir = zeros(nCells, ndir);
%     base_win = nOff-15:nOff;
%     resp_win = nOff+1:nOn;
%     base = squeeze(mean(data_dfof(base_win,:,:),1)); %averaging across baseline window
%     resp = squeeze(mean(data_dfof(resp_win,:,:),1)); %averaging across response window
%     dir_resp = zeros(nCells,ndir);
%     [x y] = ttest(resp', base', 'tail','right'); %ttest comparing base to resp for significance
%     no_match = find(isnan(x)); %?***
%     max_dir = zeros(nCells,1);
%     figure;
%     for i = 1:nCells
%         if ~sum(no_match==i) %not getting this?***
%         subplot(n, n2, i)
%             for idir = 1:ndir %for each dir
%                 if nOn>29
%                     ind = find(Dir == Dirs(idir)); %find trials of that dir
%                 else
%                     ind = find(Dir(1:ntrials-1) == Dirs(idir));
%                 end
%                 [h_dir(i,idir), p_dir(i,idir)] = ttest(resp(i,ind), base(i,ind),'tail','right','alpha', 0.05/(ndir-1)); %ttest of each cell at each dir for base vs. resp
%                 if h_dir(i,idir) %if this cell/dir is sig make red, if not make black and plot
%                     errorbar(Dirs(idir), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'or')
%                 else
%                     errorbar(Dirs(idir), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'ok')
%                 end
%                 dir_resp(i,idir) = mean(resp(i,ind)-base(i,ind),2); %avg response at each dir for each cell
%                 hold on
%             end
%             if sum(h_dir(i,:),2)>0 %if cell has one sig dir
%                 temp_resp = dir_resp(i,:);
%                 temp_resp(find(h_dir(i,:)==0)) = NaN;
%                 [max_val max_ind] = max(temp_resp,[],2); %find max index and value
%                 max_dir(i,:) = max_ind; %which index/dir was highest?
%             else
%                 [max_val max_ind] = max(dir_resp(i,:),[],2);
%                 max_dir(i,:) = max_ind;
%             end
%             title([num2str(Dirs(max_dir(i,:))) ' deg'])
%         end
%     end
% 
% %%
% nOri = nDir;
% Oris = dirs;
% nori = nOri;
% Ori = dir_mat;
% Oris = unique(Ori);
%     h_ori = zeros(nCells, nori);
%     p_ori = zeros(nCells, nori);
%     ori_resp = zeros(nCells,nori);
% 
%     max_ori = zeros(nCells,1);
%     figure;
%     for i = 1:nCells
%         if ~sum(no_match==i)
%             subplot(n, n2, i)
%             for iori = 1:nori
%                 if nOn>29
%                     ind = find(Ori == Oris(iori));
%                 else
%                     ind = find(Ori(1:ntrials-1) == Oris(iori));
%                 end
%                 [h_ori(i,iori), p_ori(i,iori)] = ttest(resp(i,ind), base(i,ind),'tail','right','alpha', 0.05/(nori-1));
%                 if h_ori(i,iori)
%                     errorbar(Oris(iori), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'or')
%                 else %should this be avg resp or avg resp-base?***
%                     errorbar(Oris(iori), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'ok')
%                 end
%                 ori_resp(i,iori) = mean(resp(i,ind)-base(i,ind),2);
%                 ori_stderror(i,iori) = std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind));
%                 hold on
%                
%             end
%             if sum(h_ori(i,:),2)>0
%                 temp_resp = ori_resp(i,:);
%                 temp_resp(find(h_ori(i,:)==0)) = NaN;
%                 [max_val max_ind] = max(temp_resp,[],2);
%                 max_ori(i,:) = max_ind;
%             else
%                 [max_val, max_ind] = max(ori_resp(i,:),[],2);
%                 max_ori(i,:) = max_ind;
%             end
%             title([num2str(Oris(max_ori(i,:))) ' deg'])
%         end
%     end
% 
% avg_resp_ori = ori_resp.';
% std_err_ori = ori_stderror.';
% 
% print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuning.pdf']),'-dpdf')
% save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_trialData.mat']),'data_dfof', 'h_ori', 'max_ori','avg_resp_ori',"std_err_ori")
% 
% 
% %%
% %replace stim ids with ori values
% 
%  [Aval, ~, indAval] = unique(max_ori);
% 
%  Avalnew = [unique(dirs)]; 
%  max_ori = Avalnew(indAval);
% pref_ori = max_ori;
% 
% save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_prefori.mat']),'pref_ori')
