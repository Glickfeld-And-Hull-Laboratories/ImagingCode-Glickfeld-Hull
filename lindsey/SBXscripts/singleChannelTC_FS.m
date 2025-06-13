%% get path names

close all
clear all
clc

date = '250207';
ImgFolder = [{'003'}];
time = strvcat('1529');
mouse = 'i2586';
doFromRef = 0;
ref = strvcat('002');
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
%LG_base = '\\CRASH.dhe.duke.edu\data\home\lindsey';

%% load and register
tic
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = [LG_base '\Data\2P_images\' mouse '\' date '\' ImgFolder{irun}];
    %CD = [LG_base '\Data\2P_images\' date '_' mouse '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder{irun} '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
    load(fName);
    
    temp(irun) = input;
    nframes = [temp(irun).counterValues{end}(end) info.config.frames];

    
    fprintf(['Reading run ' num2str(irun) '- ' num2str(min(nframes)) ' frames \r\n'])
    data_temp = sbxread(imgMatFile(1,1:11),0,min(nframes));
    if size(data_temp,1)== 2
        data_temp = data_temp(1,:,:,:);
    end
    
    if isfield(input, 'cLeverUp') 
        if irun>1
            ntrials = size(input.trialOutcomeCell,2);
            for itrial = 1:ntrials
                %temp(irun).counterValues{itrial} = bsxfun(@plus,temp(irun).counterValues{itrial},offset);
                temp(irun).cLeverDown{itrial} = temp(irun).cLeverDown{itrial}+offset;
                temp(irun).cFirstStim{itrial} = temp(irun).cFirstStim{itrial}+offset;
                temp(irun).cStimOn{itrial} = temp(irun).cStimOn{itrial}+offset;
                if ~isempty(temp(irun).cLeverUp{itrial})
                    temp(irun).cLeverUp{itrial} = temp(irun).cLeverUp{itrial}+offset;
                else
                    temp(irun).cLeverUp{itrial} = temp(irun).cLeverUp{itrial};
                end
                if ~isempty(temp(irun).cTargetOn{itrial})
                    temp(irun).cTargetOn{itrial} = temp(irun).cTargetOn{itrial}+offset;
                else
                    temp(irun).cTargetOn{itrial} = temp(irun).cTargetOn{itrial};
                end
            end
        end
    end
    offset = offset+min(nframes);
        
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    trial_n = [trial_n nframes];
end
input = concatenateDataBlocks(temp);
clear data_temp
clear temp
toc
%% Choose register interval
nep = floor(size(data,3)./10000);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*10000):300+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(300+((i-1)*10000))]); end

data_avg = mean(data(:,:,40001:40500),3);
%% Register data

if exist(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(double(data),[],[],out_bx);
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

tCyc = cell2mat(input.tCyclesOn);
cStart = cell2mat(input.cFirstStim);
cStim = cell2mat(input.cStimOn);
cTarget = celleqel2mat_padded(input.cTargetOn);

[cStimOns cStimOffs] = photoFrameFinder_Sanworks(info.frame);
cStart = cStimOns(1:2:end);
cStim = cStimOns(1:2:end);
cTarget = cStimOns(2:2:end);

nTrials = length(tCyc);
sz = size(data_reg);
data_f = zeros(sz(1),sz(2),nTrials);
data_base = zeros(sz(1),sz(2),nTrials);
data_base2 = zeros(sz(1),sz(2),nTrials);
data_targ = zeros(sz(1),sz(2),nTrials);
for itrial = 1:nTrials
    if ~isnan(cStart(itrial))
        data_f(:,:,itrial) = mean(data_reg(:,:,cStart(itrial)-20:cStart(itrial)-1),3);
        data_base(:,:,itrial) = mean(data_reg(:,:,cStart(itrial)+10:cStart(itrial)+20),3);
        if cStim(itrial) > cStart(itrial) 
            data_base2(:,:,itrial) = mean(data_reg(:,:,cStim(itrial)+9:cStim(itrial)+19),3);
        else
            data_base2(:,:,itrial) = nan(sz(1),sz(2));
        end
    else
        data_f(:,:,itrial) = nan(sz(1),sz(2));
        data_base(:,:,itrial) = nan(sz(1),sz(2));
        data_base2(:,:,itrial) = nan(sz(1),sz(2));
    end
    if ~isnan(cTarget(itrial))
        if cTarget(itrial)+19 < sz(3)
            data_targ(:,:,itrial) = mean(data_reg(:,:,cTarget(itrial)+5:cTarget(itrial)+10),3);
        else
            data_targ(:,:,itrial) = nan(sz(1),sz(2));
        end
    else
        data_targ(:,:,itrial) = nan(sz(1),sz(2));
    end
end
data_base_dfof = (data_base-data_f)./data_f;
data_base2_dfof = (data_base2-data_f)./data_f;
data_targ_dfof = (data_targ-data_f)./data_f;
targCon = celleqel2mat_padded(input.tGratingContrast);
if input.doRandCon
        baseCon = ones(size(targCon));
else
    baseCon = celleqel2mat_padded(input.tBaseGratingContrast);
end
ind_con = intersect(find(targCon == 1),find(baseCon == 0));
baseDir = celleqel2mat_padded(input.tBaseGratingDirectionDeg);
dirs = unique(baseDir);
ndir = length(dirs);
targetDelta = round(celleqel2mat_padded(input.tGratingDirectionDeg),0);
deltas = unique(targetDelta);
nDelta = length(deltas);
data_dfof_dir = zeros(sz(1),sz(2),ndir);
data_dfof2_dir = zeros(sz(1),sz(2),ndir);
[n n2] = subplotn(ndir);
figure;
for idir = 1:ndir
    ind = setdiff(find(baseDir == dirs(idir)),ind_con);
    data_dfof_dir(:,:,idir) = nanmean(data_base_dfof(:,:,ind),3);
    data_dfof2_dir(:,:,idir) = nanmean(data_base2_dfof(:,:,ind),3);
    subplot(n,n2,idir)
    imagesc(data_dfof_dir(:,:,idir))
    title(dirs(idir))
end
if sum(~isnan(data_dfof2_dir))>1
    data_dfof_dir_all = cat(3, data_dfof_dir, data_dfof2_dir);
else
    data_dfof_dir_all = data_dfof_dir;
end
data_dfof_targ = zeros(sz(1),sz(2),nDelta);
[n n2] = subplotn(nDelta);
figure;
for idir = 1:nDelta
    ind = find(targetDelta == deltas(idir));
    data_dfof_targ(:,:,idir) = nanmean(data_targ_dfof(:,:,ind),3);
    subplot(n,n2,idir)
    imagesc(data_dfof_targ(:,:,idir))
    title(deltas(idir))
end
data_dfof = cat(3,data_dfof_dir_all,data_dfof_targ);
myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_max = max(imfilter(data_dfof,myfilter),[],3);
figure;
imagesc(data_dfof_max)
data_dfof = cat(3,data_dfof,data_dfof_max);
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
%% FS cycle analysis

if iscell(input.nFramesOn)
    nOn = input.nFramesOn{1};
else
    nOn = input.nFramesOn;
end
prewin_frames = 30;
postwin_frames = 90;
tCyc = cell2mat(input.tCyclesOn);
% cStart = celleqel2mat_padded(input.cFirstStim);
% cTarget = celleqel2mat_padded(input.cTargetOn);
nTrials = length(tCyc);
nCells = size(npSub_tc,2);
maxCyc = max(tCyc,[],2);
data_trial = nan(prewin_frames+postwin_frames,nCells,maxCyc+1,nTrials);

tFramesOff = nan(nTrials,maxCyc);
SIx = strcmp(input.trialOutcomeCell, 'success');
MIx = strcmp(input.trialOutcomeCell, 'ignore');
FIx = strcmp(input.trialOutcomeCell, 'failure');
nCyc = tCyc;
nCyc([find(MIx) find(SIx)]) = tCyc([find(MIx) find(SIx)])+1;
for itrial = 1:nTrials
    if isfield(input, 'tFramesOff')
        if length(input.tFramesOff{itrial}>0)
            tempFramesOff = input.tFramesOff{itrial};
        else
            tempFramesOff = input.nFramesOff{itrial}.*(ones(1,tCyc(itrial)));
            input.tFramesOff{itrial} = tempFramesOff;
        end
    else
        if iscell(input.nFramesOff)
            tempFramesOff = input.nFramesOff{itrial}.*(ones(1,tCyc(itrial)));
        else
            tempFramesOff = input.nFramesOff.*(ones(1,tCyc(itrial)));
        end
    end

    tFramesOff(itrial,1:tCyc(itrial)) = tempFramesOff(1:tCyc(itrial));
    if ~isnan(cStart(itrial))
        for icyc = 1:nCyc(itrial)
            if icyc > 1
                cyc_add = ((icyc-1)*nOn)+sum(tempFramesOff(1:icyc-1));
            else
                cyc_add = 0;
            end
            if cStart(itrial)+postwin_frames-1+cyc_add <= size(npSub_tc,1)
                data_trial(:,:,icyc,itrial) = npSub_tc(cStart(itrial)-prewin_frames+cyc_add:cStart(itrial)+postwin_frames+cyc_add-1,:);
            else
                data_trial(:,:,icyc,itrial) = NaN(prewin_frames+postwin_frames,nCells);
            end 
        end
    else
        data_trial(:,:,icyc,itrial) = NaN(prewin_frames+postwin_frames,nCells);
    end
end
data_f = nanmean(data_trial(1:prewin_frames,:,1,:),1);
data_dfof = bsxfun(@rdivide,bsxfun(@minus,data_trial,data_f),data_f);

targCon = celleqel2mat_padded(input.tGratingContrast);
if isfield(input,'doRandCon') & input.doRandCon
	baseCon = nan(maxCyc,nTrials);
    for itrial = 1:nTrials
        baseCon(:,itrial) = input.tBaseGratingContrast{itrial}(1:tCyc(itrial));
    end
    ind_con = [];
else
    baseCon = celleqel2mat_padded(input.tBaseGratingContrast);
    ind_con = intersect(find(targCon == 1),find(baseCon == 0));
end
baseDir = celleqel2mat_padded(input.tBaseGratingDirectionDeg);
dirs = unique(baseDir);
ndir = length(dirs);
tGratingDir = round(double(celleqel2mat_padded(input.tGratingDirectionDeg)),0);
if sum(tGratingDir-baseDir) == 0
    targetDelta = tGratingDir-baseDir;
else
    targetDelta = tGratingDir;
end
deltas = unique(targetDelta);
nDelta = length(deltas);
offs = unique(tFramesOff(:,1));
noff = length(offs);
frameRateHz = input.frameRateHz;

base_win =31:33;
resp_win =36:38; 

% figure;
% if nCells<25
%     ii = nCells;
% else
%     ii = 25;
% end
% for i = 1:ii
%     subplot(5,5,i)
% if length(ind_con)>10
%     plot(squeeze(nanmean(mean(data_dfof(20:50,i,2,ind_con),2),4)))
% elseif noff>1
%     ind = find(tFramesOff(:,1) == offs(noff));
%     plot(squeeze(nanmean(mean(data_dfof(20:50,i,1,:),2),4)))
% else
%     plot(squeeze(nanmean(mean(data_dfof(20:50,i,1,:),2),4)))
% end
% vline(base_win-19)
% vline(resp_win-19)
% end

figure;
subplot(2,1,1)
plot(squeeze(nanmean(mean(data_dfof(:,:,1,:),2),4)));
vline(base_win,'k:')
vline(resp_win,'r:')
title('Baseline')
subplot(2,1,2)
sz = size(data_dfof);
data_targ = zeros(sz(1),sz(2),length([find(SIx)]));
for itrial = 1:sz(4);
    %if find([find(SIx)] == itrial)
        data_targ(:,:,itrial) = data_dfof(:,:,nCyc(itrial),itrial);
    %end
end
plot(squeeze(nanmean(mean(data_targ,2),3)));
title('Target')
vline(base_win,'k:')
vline(resp_win,'r:')

%%

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']), 'data_dfof', 'prewin_frames', 'postwin_frames', 'base_win','resp_win')
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']),'prewin_frames','baseDir', 'dirs', 'ndir', 'tFramesOff', 'offs', 'noff', 'baseCon', 'ind_con', 'tGratingDir', 'targetDelta', 'deltas', 'nDelta', 'tCyc', 'nCyc', 'maxCyc', 'nCells', 'frameRateHz', 'nTrials', 'SIx', 'MIx', 'FIx', 'cTarget', 'cStart', 'base_win','resp_win')
