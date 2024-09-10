%% get path names
close all 
clear all global
clc
date = '240625';
ImgFolder = {'002'};
time = strvcat('1507');
mouse = 'i1396';
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
data_avg = mean(data(:,:,80001:80500),3);
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

data_dfof_stim = zeros(sz(1),sz(2),20);
for i = 1:20
    corr_map = zeros(sz(1),sz(2));
    for ix = b:sz(2)-b
        for iy = b:sz(1)-b
            block = reshape(data_reg(iy-1:iy+1,ix-1:ix+1,1+5000.*(i-1):5000.*i),[9 5000]);
            TC = block(5,:);
            surround = nanmean(block([1:4 6:9],:),1);
            corr_map(iy,ix) = corr(double(TC'),surround','rows','complete');
        end
    end
    data_dfof_stim(:,:,i) = corr_map;
end
data_dfof_stim(:,:,i+1) = max(data_dfof_stim,[],3);

data_dfof = data_dfof_stim;

figure;
imagesc(data_dfof_stim(:,:,end))
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
cStimOne = cell2mat(input.cStimOneOn);

nTrials = length(cStimOne);
sz = size(npSub_tc);

tStimOne = celleqel2mat_padded(input.tstimOne);
stimOne = unique(tStimOne);
nStim = length(stimOne);
nImage = length(stimOne);

tc_one = nan(30,sz(2),nTrials);
   
for itrial = 1:nTrials
    if ~isnan(cStimOne(itrial)) & (cStimOne(itrial)+19)<sz(1)
        tc_one(:,:,itrial) = npSub_tc(cStimOne(itrial)-10:cStimOne(itrial)+19,:);
    end
end
tc_one_f = mean(tc_one(1:10,:,:));
tc_one_dfof = (tc_one-tc_one_f)./tc_one_f;

base_win = 10:12;
resp_win = 16:18;
figure;
shadedErrorBar(1:30,squeeze(nanmean(nanmean(tc_one_dfof(:,:,:),3),2)),squeeze(nanstd(nanmean(tc_one_dfof(:,:,:),3),[],2))./sqrt(5));%-mean(tc_one_dfof_all(base_win,:,it),1),2)))
vline([base_win resp_win])

%% trial responses
dfof_resp_one = squeeze(mean(tc_one_dfof(resp_win,:,:),1))-squeeze(mean(tc_one_dfof(base_win,:,:),1));
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']), 'dfof_resp_one', 'base_win', 'resp_win','tStimOne','tc_one_dfof','stimOne','tStimOne','nImage')



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
