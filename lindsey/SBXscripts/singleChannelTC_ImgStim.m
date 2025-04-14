%% get path names
close all 
clear all global
clc
date = '250307';
ImgFolder = {'002'};
time = strvcat('1108');
mouse = 'i1411';
doFromRef = 0;
ref = strvcat('002');
nrun = size(ImgFolder,2);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);


LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

%% load
tic
data = [];
CD = [LG_base '\Data\2P_images\' mouse '\' date '\' ImgFolder{1}];
%CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\AW68\two-photon imaging\' date '\' ImgFolder(irun,:)];
%CD = [LG_base '\Data\2P_images\' mouse '-KM' '\' date '_' mouse '\' ImgFolder(irun,:)];
%CD = ['\\CRASH.dhe.duke.edu\data\home\kevin\Data\2P\' date '_' mouse '\' ImgFolder(irun,:)];
cd(CD);
imgMatFile = [ImgFolder{1} '_000_000.mat'];
load(imgMatFile);
fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time '.mat'];
load(fName);
nframes = [input.counterValues{end}(end) info.config.frames];
fprintf(['Reading run ' ImgFolder{1} '- ' num2str(min(nframes)) ' frames \r\n'])
data = sbxread(imgMatFile(1,1:11),1,min(nframes));
if size(data,1)== 2
    data = data(1,:,:,:);
end
data = squeeze(data);
toc

%% Choose register interval
nep = floor(size(data,3)./10000);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end

%% Register data
data_avg = mean(data(:,:,100001:100500),3);
if exist(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(data,[],[],out);
    clear out outs
elseif doFromRef
    ref_str = ['runs-' ref];
    if size(ref,1)>1
        ref_str = [ref_str '-' ref(size(ref,1),:)];
    end
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_reg_shifts.mat']))
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    %load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_mask_cell.mat']))
    %load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_trialData.mat']))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
else
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
clear data out

%% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_byFrame.pdf']),'-dpdf', '-bestfit')

figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

figure; n = [1 nep]; for ii = 1:2; subplot(2,1,ii); i = n(ii); imagesc(mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end

%% find activated cells

[stimOns stimOffs] = photoFrameFinder_Sanworks(info.frame);
cStimOne = stimOns(1:2:end);
cStimTwo = stimOns(2:2:end);

nTrials = length(cStimOne);
sz = size(data_reg);
data_f = nan(sz(1),sz(2),nTrials);
data_one = nan(sz(1),sz(2),nTrials);
for itrial = 1:nTrials
    if ~isnan(cStimOne(itrial)) & (cStimOne(itrial)+30)<sz(3)
        data_f(:,:,itrial) = mean(data_reg(:,:,cStimOne(itrial)-20:cStimOne(itrial)-1),3);
        data_one(:,:,itrial) = mean(data_reg(:,:,cStimOne(itrial)+5:cStimOne(itrial)+25),3);
    end
end
data_one_dfof = (data_one-data_f)./data_f;
clear data_one


stim_mat = celleqel2mat_padded(input.tstimOne);
stims = unique(stim_mat);
nStim = length(stims);
data_dfof_stim = zeros(sz(1),sz(2),nStim);
figure;
start = 1;
for is = 1:nStim
    ind = find(stim_mat == stims(is));
    data_dfof_stim(:,:,start) = nanmean(data_one_dfof(:,:,ind),3);
    subplot(3,3,start)
    imagesc(data_dfof_stim(:,:,start))
    start = 1+start;
end

ind = 1:nTrials;
data_dfof_stim = cat(3,data_dfof_stim,nanmean(data_one_dfof(:,:,ind),3));

subplot(3,3,start)
imagesc(data_dfof_stim(:,:,end))
title('Avg all')

data_dfof = cat(3,data_dfof_stim,max(data_dfof_stim,[],3));
clear data_f
%% cell segmentation 
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
% bwout = imCellEditInteractive(data_dfof_max);
% mask_cell = bwlabel(bwout);


%% neuropil mask and subtraction
mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof', 'mask_cell', 'mask_np')

clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_data_temp mask_exp data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 

% neuropil subtraction
down = 5;
sz = size(data_reg);

data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_reg_down  = stackGroupProject(data_reg,down);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
nCells = size(data_tc,2);
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '%s/n']) 
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
clear data_reg data_reg_down

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell
%% Stim two analysis

nTrials = length(cStimOne);
[sz(3) nCells] = size(npSub_tc);

tc_one = nan(100,nCells,nTrials);
tc_two = nan(100,nCells,nTrials);
   
for itrial = 1:nTrials
    if ~isnan(cStimOne(itrial)) & ~isnan(cStimTwo(itrial)) & (cStimOne(itrial)+79)<sz(3) & (cStimTwo(itrial)+19)<sz(3)
        tc_one(:,:,itrial) = npSub_tc(cStimOne(itrial)-20:cStimOne(itrial)+79,:);
        tc_two(:,:,itrial) = npSub_tc(cStimTwo(itrial)-20:cStimTwo(itrial)+79,:);
    end
end
tc_one_f = mean(tc_one(1:20,:,:));
tc_two_f = mean(tc_two(1:20,:,:));
tc_one_dfof = (tc_one-tc_one_f)./tc_one_f;
tc_two_dfof = (tc_two-tc_one_f)./tc_one_f;

base_win = 20:22;
resp_win = 24:26;
figure;
subplot(2,1,1)
shadedErrorBar(1:100,squeeze(nanmean(nanmean(tc_one_dfof(:,:,:),3),2)),squeeze(nanstd(nanmean(tc_one_dfof(:,:,:),3),[],2))./sqrt(5));%-mean(tc_one_dfof_all(base_win,:,it),1),2)))
vline([base_win resp_win])
subplot(2,1,2)
shadedErrorBar(1:100,squeeze(nanmean(nanmean(tc_two_dfof(:,:,:),3),2)),squeeze(nanstd(nanmean(tc_two_dfof(:,:,:),3),[],2))./sqrt(5));%-mean(tc_one_dfof_all(base_win,:,it),1),2)))
vline([base_win resp_win])

resp_mat = zeros(nCells,nTrials,2);
resp_mat(:,:,1) = squeeze(mean(tc_one_dfof(resp_win,:,:),1)-mean(tc_one_dfof(base_win,:,:),1));
resp_mat(:,:,2) = squeeze(mean(tc_two_dfof(resp_win,:,:),1)-mean(tc_two_dfof(base_win,:,:),1));

h1_ori = zeros(nCells,2,nStim);
h2_ori = zeros(nCells,2,nStim);
p1_ori = zeros(nCells,2,nStim);
p2_ori = zeros(nCells,2,nStim);
resp_stim = zeros(nCells,2,nStim);
%cells stim1/stim2 img 

for is = 1:nStim
    ind = find(stim_mat == stims(is));
    [h1_ori(:,1,is) p1_ori(:,1,is)] = ttest(squeeze(mean(tc_one_dfof(resp_win,:,ind),1)),squeeze(mean(tc_one_dfof(base_win,:,ind),1)),'tail','right','dim',2,'alpha',0.05./(nStim));
    [h2_ori(:,2,is) p2_ori(:,2,is)] = ttest(squeeze(mean(tc_two_dfof(resp_win,:,ind),1)),squeeze(mean(tc_two_dfof(base_win,:,ind),1)),'tail','right','dim',2,'alpha',0.05./(nStim));
    resp_stim(:,:,is) = squeeze(mean(resp_mat(:,ind,:),2,'omitnan'));
end
   
[h1 p1] = ttest(squeeze(mean(tc_one_dfof(resp_win,:,:),1)),squeeze(mean(tc_one_dfof(base_win,:,:),1)),'tail','right','dim',2);
[h2 p2] = ttest(squeeze(mean(tc_two_dfof(resp_win,:,:),1)),squeeze(mean(tc_two_dfof(base_win,:,:),1)),'tail','right','dim',2);

good_ind = unique([find(h1); find(sum(sum(sum(h1_ori,3),4),5))]);



save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']), 'resp_mat', 'cStimOne','cStimTwo', 'stim_mat', 'stims', 'nStim','good_ind', 'base_win', 'resp_win', 'h1','h1_ori');
        
%% if do second run for tuning
ImgFolder = {'003'};
time = strvcat('1336');
singleStim = 0;
ref_str = run_str;
run_str = catRunName(ImgFolder, nrun);

data = [];
CD = [LG_base '\Data\2P_images\' mouse '\' date '\' ImgFolder{1}];
cd(CD);
imgMatFile = [ImgFolder{1} '_000_000.mat'];
load(imgMatFile);
fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time '.mat'];
load(fName);
nframes = [input.counterValues{end}(end) info.config.frames];
fprintf(['Reading run ' ImgFolder{1} '- ' num2str(min(nframes)) ' frames \r\n'])
data = sbxread(imgMatFile(1,1:11),1,min(nframes));
if size(data,1)== 2
    data = data(1,:,:,:);
end
data = squeeze(data);
toc

load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_reg_shifts.mat']), 'data_avg')
[out, data_reg] = stackRegister(data,data_avg);
mkdir(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')

[stimOns stimOffs] = photoFrameFinder_Sanworks(info.frame);
if singleStim
    cStimOne = stimOns;
else
    cStimOne = stimOns(1:2:end);
    cStimTwo = stimOns(2:2:end);
end

nTrials = length(cStimOne);
sz = size(data_reg);
data_f = nan(sz(1),sz(2),nTrials);
data_one = nan(sz(1),sz(2),nTrials);
for itrial = 1:nTrials
    if ~isnan(cStimOne(itrial)) & (cStimOne(itrial)+30)<sz(3)
        data_f(:,:,itrial) = mean(data_reg(:,:,cStimOne(itrial)-20:cStimOne(itrial)-1),3);
        data_one(:,:,itrial) = mean(data_reg(:,:,cStimOne(itrial)+5:cStimOne(itrial)+25),3);
    end
end
data_one_dfof = (data_one-data_f)./data_f;
clear data_one

stim_mat = celleqel2mat_padded(input.tstimOne);
stims = unique(stim_mat);
nStim = length(stims);

load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_mask_cell.mat']))
data_dfof_stim = zeros(sz(1),sz(2),nStim);
[n1 n2] = subplotn(nStim);
figure;
start = 1;
for is = 1:nStim
    ind = find(stim_mat == stims(is));
    data_dfof_stim(:,:,start) = nanmean(data_one_dfof(:,:,ind),3);
    shade_img = imShade(data_dfof_stim(:,:,start), mask_cell);
    subplot(n1,n2,start)
    imagesc(shade_img)
    start = 1+start;
end

% neuropil subtraction
down = 5;
sz = size(data_reg);

data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_reg_down  = stackGroupProject(data_reg,down);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
nCells = size(data_tc,2);
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '%s/n']) 
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
clear data_reg data_reg_down

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

nTrials = length(cStimOne);
[sz(3) nCells] = size(npSub_tc);

tc_one = nan(100,nCells,nTrials);
tc_two = nan(100,nCells,nTrials);
   
for itrial = 1:nTrials
    if ~isnan(cStimOne(itrial)) & (cStimOne(itrial)+79)<sz(3)
        tc_one(:,:,itrial) = npSub_tc(cStimOne(itrial)-20:cStimOne(itrial)+79,:);
        if ~singleStim
            if ~isnan(cStimTwo(itrial)) & (cStimTwo(itrial)+79)<sz(3)
                tc_two(:,:,itrial) = npSub_tc(cStimTwo(itrial)-20:cStimTwo(itrial)+79,:);
            end
        end
    end
end
tc_one_f = mean(tc_one(1:20,:,:));
tc_one_dfof = (tc_one-tc_one_f)./tc_one_f;
tc_two_dfof = (tc_two-tc_one_f)./tc_one_f;

%resp_win = 26:28;
figure;
shadedErrorBar(1:100,squeeze(nanmean(nanmean(tc_one_dfof(:,:,:),3),2)),squeeze(nanstd(nanmean(tc_one_dfof(:,:,:),3),[],2))./sqrt(5));%-mean(tc_one_dfof_all(base_win,:,it),1),2)))
vline([base_win resp_win])
    
resp_mat = zeros(nCells,nTrials,2);
resp_mat(:,:,1) = squeeze(mean(tc_one_dfof(resp_win,:,:),1)-mean(tc_one_dfof(base_win,:,:),1));
resp_mat(:,:,2) = squeeze(mean(tc_two_dfof(resp_win,:,:),1)-mean(tc_two_dfof(base_win,:,:),1));

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']), 'resp_mat', 'cStimOne','stim_mat', 'stims', 'nStim', 'base_win', 'resp_win');
