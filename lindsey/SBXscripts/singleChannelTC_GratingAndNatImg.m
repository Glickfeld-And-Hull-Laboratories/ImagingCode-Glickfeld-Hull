%% get path names
close all 
clear all global
clc
date = '251018';
ImgFolder = {'002'};
time = strvcat('1446');
mouse = 'i1426';
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
data_avg = mean(data(:,:,40001:40500),3);
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

%% find activated cells

cStimOne = cell2mat(input.cStimOneOn);
cStimTwo = cell2mat(input.cStimTwoOn);

nTrials = length(cStimOne);
sz = size(data_reg);
data_f = nan(sz(1),sz(2),nTrials);
data_one = nan(sz(1),sz(2),nTrials);
data_tc = nan(sz(1),sz(2),41,nTrials);
for itrial = 1:nTrials
    if ~isnan(cStimOne(itrial)) & (cStimOne(itrial)+15)<sz(3)
        data_f(:,:,itrial) = mean(data_reg(:,:,cStimOne(itrial)-20:cStimOne(itrial)-1),3);
        data_one(:,:,itrial) = mean(data_reg(:,:,cStimOne(itrial)+5:cStimOne(itrial)+15),3);
    end
end
data_one_dfof = (data_one-data_f)./data_f;
clear data_one

tGrating = celleqel2mat_padded(input.tDoGratingStim);
tStimOne = celleqel2mat_padded(input.tstimOne);
stimOne = unique(tStimOne);
nImage = length(stimOne);
tGratingSF = celleqel2mat_padded(input.tGratingSpatialFreqCPD);
gratingSFs = unique(tGratingSF(find(tGrating)));
nGrating = length(gratingSFs);

nStim = nImage + nGrating;

if nStim>50
    data_dfof_avg = zeros(sz(1),sz(2),20);
    for i = 1:20
        ind = randperm(nTrials,50);
        data_dfof_avg(:,:,i) = mean(data_one_dfof(:,:,ind),3);
    end
    data_dfof_max = max(data_dfof_avg,[],3);
    
    data_dfof_stim = zeros(sz(1),sz(2),21);
    b = 5;
    for i = 1:20
        ind = randperm(nTrials,50);
        corr_map = zeros(sz(1),sz(2));
        for ix = b:sz(2)-b
            for iy = b:sz(1)-b
                block = reshape(data_one_dfof(iy-1:iy+1,ix-1:ix+1,ind),[9 50]);
                TC = block(5,:);
                surround = mean(block([1:4 6:9],:),1);
                corr_map(iy,ix) = corr(TC',surround','rows','complete');
            end
        end
        data_dfof_stim(:,:,i) = corr_map;
    end
else
    data_dfof_stim = zeros(sz(1),sz(2),nStim);
    figure;
    [n n2] = subplotn(nImage);
    for i = 1:nImage
        ind = intersect(find(tGrating==0), find(tStimOne == stimOne(i)));
        data_dfof_stim(:,:,i) = mean(data_one_dfof(:,:,ind),3,'omitnan');
        subplot(n,n2,i)
        imagesc(data_dfof_stim(:,:,i))
    end
    figure;
    [n n2] = subplotn(nGrating);
    for i = 1:nGrating
        ind = intersect(find(tGrating==1), find(tGratingSF == gratingSFs(i)));
        data_dfof_stim(:,:,i+nImage) = mean(data_one_dfof(:,:,ind),3,'omitnan');
        subplot(n,n2,i)
        imagesc(data_dfof_stim(:,:,i+nImage))
    end
end
data_dfof_stim(:,:,nStim+1) = mean(data_one_dfof,3,'omitnan');
data_dfof_stim(:,:,nStim+2) = max(data_dfof_stim,[],3);
data_dfof = data_dfof_stim;

figure;
imagesc(data_dfof_stim(:,:,nStim+2))
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

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc', 'cStimOne')

clear data_tc data_tc_down np_tc np_tc_down mask_np
%% Stim analysis

tc_one = nan(50,nCells,nTrials);
tc_two = nan(50,nCells,nTrials);
   
for itrial = 1:nTrials
    if ~isnan(cStimOne(itrial)) & (cStimOne(itrial)+29)<sz(3)
        tc_one(:,:,itrial) = npSub_tc(cStimOne(itrial)-20:cStimOne(itrial)+29,:);
    end
    if ~isnan(cStimTwo(itrial)) & (cStimTwo(itrial)+29)<sz(3)
        tc_two(:,:,itrial) = npSub_tc(cStimTwo(itrial)-20:cStimTwo(itrial)+29,:);
    end
end
tc_one_f = mean(tc_one(1:20,:,:));
tc_one_dfof = (tc_one-tc_one_f)./tc_one_f;
tc_two_dfof = (tc_two-tc_one_f)./tc_one_f;

base_win = 20:22;
resp_win = 29:34;
figure;
subplot(2,1,1)
shadedErrorBar(1:50,squeeze(nanmean(nanmean(tc_one_dfof(:,:,:),3),2)),squeeze(nanstd(nanmean(tc_one_dfof(:,:,:),3),[],2))./sqrt(nCells));%-mean(tc_one_dfof_all(base_win,:,it),1),2)))
vline([base_win resp_win])
subplot(2,1,2)
shadedErrorBar(1:50,squeeze(nanmean(nanmean(tc_two_dfof(:,:,:),3),2)),squeeze(nanstd(nanmean(tc_two_dfof(:,:,:),3),[],2))./sqrt(nCells));%-mean(tc_one_dfof_all(base_win,:,it),1),2)))
vline([base_win resp_win])

%% trial responses
dfof_resp_one = squeeze(mean(tc_one_dfof(resp_win,:,:),1));
dfof_base_one = squeeze(mean(tc_one_dfof(base_win,:,:),1));
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']), 'dfof_resp_one', 'dfof_base_one', 'base_win', 'resp_win','tStimOne','tc_one_dfof','tGrating','stimOne','tStimOne','nImage','tGratingSF','gratingSFs','nGrating')


%% stim responses
dfof_resp_stim = zeros(nCells,nStim);
h_stim = zeros(nCells,nStim);
p_stim = zeros(nCells,nStim);

for i = 1:nImage
    ind = intersect(find(tGrating==0),find(tStimOne==stimOne(i)));
    [h_stim(:,i) p_stim(:,i)] = ttest(dfof_resp_one(:,ind),dfof_base_one(:,ind),'Dim', 2, 'Tail','right');
    dfof_resp_stim(:,i) = nanmean(dfof_resp_one(:,ind),2);
end
for i = 1:nGrating
    ind = intersect(find(tGrating==1),find(tGratingSF==gratingSFs(i)));
    [h_stim(:,i+nImage) p_stim(:,i+nImage)] = ttest(dfof_resp_one(:,ind),dfof_base_one(:,ind),'Dim', 2, 'Tail','right');
    dfof_resp_stim(:,i+nImage) = nanmean(dfof_resp_one(:,ind),2);
end

image_resp = find(sum(h_stim(:,1:nImage),2));
grating_resp = find(sum(h_stim(:,nImage+1:nImage+nGrating),2));

h_stim_thresh = h_stim & dfof_resp_stim>0.05;

mask_props = regionprops(mask_cell);
nResp = sum(h_stim,1);
nResp_thresh = sum(h_stim_thresh,1);

for i = 1:nStim
    ind_cell = find(h_stim(:,i));
    centroid_dist{i} = [];
    centroid_min_dist{i} = [];
    for iC1 = 1:nResp(i)
        all_dist = [];
        for iC2 = 1:nResp(i)
            if iC1<iC2
                centroid_diff = mask_props(ind_cell(iC1)).Centroid - mask_props(ind_cell(iC2)).Centroid;
                all_dist = [all_dist sqrt(centroid_diff(1).^2 + centroid_diff(2).^2)];
            end
        end
        centroid_dist{i} = [centroid_dist{i} all_dist];
        centroid_min_dist{i} = [centroid_min_dist{i} min(all_dist)];
    end
end
mean_dist = cellfun(@mean,centroid_dist);
mean_min_dist = cellfun(@mean,centroid_min_dist);
figure;
subplot(2,2,1)
scatter(mean_dist,nResp);
subplot(2,2,2)
scatter(mean_min_dist,nResp);
subplot(2,2,3)
scatter(1:nStim,nResp)
ylim([0 400])
subplot(2,2,4)
for i = 1:nStim
    errorbar(i,nanmean(dfof_resp_stim(find(h_stim(:,i)),i),1), nanstd(dfof_resp_stim(find(h_stim(:,i)),i),[],1)./sqrt(length(find(h_stim(:,i)))))
    hold on
end
ylim([0 0.5])

fractOverlap = nan(nStim,nStim);
for iStim1 = 1:nStim
    for iStim2 = 1:nStim
        if iStim1<iStim2
            ind1 = find(h_stim_thresh(:,iStim1));
            ind2 = find(h_stim_thresh(:,iStim2));
            tot_cells = length(unique([ind1; ind2]));
            same_cells = length(intersect(ind1,ind2));
            fractOverlap(iStim1,iStim2) = same_cells./tot_cells;
        end
    end
end
figure; imagesc(fractOverlap);
ax = gca; ax.CLim = [0 1];
%% running speed
wheel_speed = wheelSpeedCalc(input,32,'purple');
figure;
plot(wheel_speed)
wheel_trial = nan(1,nTrials);
for itrial = 1:nTrials
    if ~isnan(cStimOne(itrial)) & ~isnan(cStimTwo(itrial)) & (cStimOne(itrial)+14)<sz(3)
        wheel_trial(:,itrial) = mean(wheel_speed(:,cStimOne(itrial):cStimOne(itrial)+14),2);
    end
end
figure; 
plot(wheel_trial)
%% pupil data
fn = [cell2mat(ImgFolder) '_000_000_eye.mat'];

%load data
data_temp = load(fn);
data_temp = squeeze(data_temp.data);

%crop frames to match mworks data
nFrames = input.counterValues{end}(end);
data = data_temp(:,:,1:nFrames);      % the raw images...
[data_crop rect] = cropEyeData(data);

%% measure pupil position/diameter
rad_range = [3 25]; %adjust to expected range of pupil size (if low end is too small then may find noisy bright stuff)
Eye_data = extractEyeData(data_crop,rad_range);
[rad centroid] = alignEyeData(Eye_data,input);

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_eyeAndWheel.mat']), 'wheel_speed', 'wheel_trial', 'centroid', 'rad', 'rect', 'rad_range')
