%% get path names
close all 
clear all global
clc
date = '240129';
ImgFolder = {'002'};
time = strvcat('1421');
mouse = 'i1381';
doFromRef = 0;
ref = strvcat('002');
nrun = size(ImgFolder,2);
frame_rate = 15;
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
data_avg = mean(data(:,:,50001:50500),3);
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

cStimOn = cell2mat(input.cStimTwoOn);
nTrials = length(cStimOn);

sz = size(data_reg);
data_f = nan(sz(1),sz(2),nTrials);
data_one = nan(sz(1),sz(2),nTrials);
for itrial = 1:nTrials
    if ~isnan(cStimOn(itrial)) & (cStimOn(itrial)+30)<sz(3)
        data_f(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)-20:cStimOn(itrial)-1),3);
        data_one(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)+5:cStimOn(itrial)+25),3);
    end
end
data_resp_dfof = (data_one-data_f)./data_f;
clear data_one data_f

adaptStim = celleqel2mat_padded(input.tStimOneGratingContrast);
stimCon_all = celleqel2mat_padded(input.tStimTwoGratingContrast);
maskCon_all = celleqel2mat_padded(input.tMaskTwoGratingContrast);
stimCons = unique(stimCon_all);
maskCons = unique(maskCon_all);
nStimCon = length(stimCons);
nMaskCon = length(maskCons);

stimDir_all = celleqel2mat_padded(input.tStimTwoGratingDirectionDeg);
maskDir_all = rad2deg(wrapTo2Pi(deg2rad(celleqel2mat_padded(input.tMaskTwoGratingDirectionDeg))));
maskDir_all(find(maskDir_all==360)) = 0;
stimDirs = unique(stimDir_all);
nStimDir = length(stimDirs);
maskDirs = unique(maskDir_all);
nMaskDir = length(maskDirs);

ind_stimAlone = intersect(find(stimCon_all),find(maskCon_all==0));
ind_plaid = intersect(find(stimCon_all),find(maskCon_all));

adaptStim = celleqel2mat_padded(input.tStimOneGratingContrast);
ind_adapt = find(adaptStim);
ind_noadapt = find(adaptStim==0);

nStim = nStimDir*2;

if nStimDir > 1 & ~input.doTwoStimTogether
    data_dfof = zeros(sz(1),sz(2), nStim);
    start = 1;
    for is = 1:nStimDir
        ind_stim = intersect(intersect(ind_stimAlone,find(stimDir_all == stimDirs(is))),find(adaptStim==0));
        data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,ind_stim),3);
        figure(1)
        subplot(3,4,is)
        imagesc(data_dfof(:,:,start))
        start = start+1;
        ind_plaidstim = intersect(intersect(ind_plaid,find(stimDir_all == stimDirs(is))),find(adaptStim==0));
        data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,ind_plaidstim),3);        
        figure(2)
        subplot(3,4,is)
        imagesc(data_dfof(:,:,start))
        start = start+1;
    end
    figure; 
    subplot(2,1,1); 
    imagesc(mean(data_dfof(:,:,1:2:end),3,"omitnan"))
    title('Grating')
    colormap gray
    subplot(2,1,2); 
    imagesc(mean(data_dfof(:,:,2:2:end),3,"omitnan"))
    title('Plaid')
    colormap gray
    sgtitle([mouse ' ' date])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_GratingVsPlaid.pdf']),'-dpdf')
end

data_dfof = cat(3,data_dfof,max(data_dfof,[],3));

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']), 'cStimOn', 'maskCon_all', 'stimCon_all', 'stimCons', 'maskCons', 'nStimCon', 'nMaskCon', 'stimDir_all', 'stimDirs', 'nStimDir', 'maskDir_all', 'maskDirs', 'nMaskDir', 'adaptStim', 'ind_adapt', 'ind_noadapt', 'ind_stimAlone', 'ind_plaid', 'frame_rate', 'nTrials')

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

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell

%%
prewin_frames = frame_rate;
nFramesOn = unique(celleqel2mat_padded(input.nStimTwoFramesOn));
postwin_frames = frame_rate*2;
tt = (1-prewin_frames:postwin_frames).*(1/frame_rate);
data_resp = nan(prewin_frames+postwin_frames,nCells,nTrials);
data_f = nan(1,nCells,nTrials);

for itrial = 1:nTrials
    if cStimOn(itrial) + postwin_frames < sz(3)
        data_resp(:,:,itrial) = npSub_tc(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
    end
end

data_f = mean(data_resp(prewin_frames-ceil(prewin_frames/4):prewin_frames,:,:),1);
data_dfof_tc = (data_resp-data_f)./data_f;

resp_cell_noadapt = cell(nStimDir,2);
resp_cell_adapt = cell(nStimDir,2);
base_cell_noadapt = cell(nStimDir,2);
base_cell_adapt = cell(nStimDir,2);
trialInd_adapt = cell(nStimDir,2);
trialInd_noadapt = cell(nStimDir,2);
trialsperstim_noadapt = zeros(nStimDir,2);
trialsperstim_adapt = zeros(nStimDir,2);
h_resp =zeros(nCells,nStimDir,2);
p_resp =zeros(nCells,nStimDir,2);

base_win = prewin_frames-ceil(frame_rate/4):prewin_frames;
resp_win = prewin_frames+5:prewin_frames+nFramesOn;
data_dfof_dir_tc_avg_adapt = nan(prewin_frames+postwin_frames, nCells, nStimDir, 2);
data_dfof_dir_tc_avg_noadapt = nan(prewin_frames+postwin_frames, nCells, nStimDir, 2);

avg_resp_dir = zeros(nCells, nStimDir, 2, 2, 2); % 1- cells; 2- directions; 3- grating/plaid; 4- no adapt/sing adapt; 5- mean/sem

all_resp_dir = [];
all_resp_plaid = [];
all_dir = [];
all_plaid = [];
for iDir = 1:nStimDir
    ind_stimdir = find(stimDir_all == stimDirs(iDir));
    ind_diralone = intersect(ind_stimdir, ind_stimAlone);
    ind_dirplaid = intersect(ind_stimdir, ind_plaid);
    trialsperstim_noadapt(iDir,1) = length(intersect(ind_noadapt,ind_diralone));
    trialsperstim_noadapt(iDir,2) = length(intersect(ind_noadapt,ind_dirplaid));
    trialInd_noadapt{iDir,1} = intersect(ind_noadapt,ind_diralone);
    trialInd_noadapt{iDir,2} = intersect(ind_noadapt,ind_dirplaid);
    resp_cell_noadapt{iDir,1} = squeeze(mean(data_dfof_tc(resp_win,:,intersect(ind_noadapt,ind_diralone)),1));
    resp_cell_noadapt{iDir,2} = squeeze(mean(data_dfof_tc(resp_win,:,intersect(ind_noadapt,ind_dirplaid)),1));
    base_cell_noadapt{iDir,1} = squeeze(mean(data_dfof_tc(base_win,:,intersect(ind_noadapt,ind_diralone)),1));
    base_cell_noadapt{iDir,2} = squeeze(mean(data_dfof_tc(base_win,:,intersect(ind_noadapt,ind_dirplaid)),1));
    avg_resp_dir(:,iDir,1,1,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,intersect(ind_noadapt,ind_diralone)),1),3,'omitnan'));
    avg_resp_dir(:,iDir,1,1,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,intersect(ind_noadapt,ind_diralone)),1),[],3,'omitnan')./sqrt(length(intersect(ind_noadapt,ind_diralone))));
    avg_resp_dir(:,iDir,2,1,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,intersect(ind_noadapt,ind_dirplaid)),1),3,'omitnan'));
    avg_resp_dir(:,iDir,2,1,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,intersect(ind_noadapt,ind_dirplaid)),1),[],3,'omitnan')./sqrt(length(intersect(ind_noadapt,ind_dirplaid))));
    data_dfof_dir_tc_avg_noadapt(:,:,iDir,1) = squeeze(mean(data_dfof_tc(:,:,intersect(ind_noadapt,ind_diralone)),3,'omitnan'));
    data_dfof_dir_tc_avg_noadapt(:,:,iDir,2) = squeeze(mean(data_dfof_tc(:,:,intersect(ind_noadapt,ind_dirplaid)),3,'omitnan'));
    all_resp_dir = [all_resp_dir squeeze(mean(data_dfof_tc(resp_win,:,intersect(ind_noadapt,ind_diralone)),1))];
    all_resp_plaid = [all_resp_plaid squeeze(mean(data_dfof_tc(resp_win,:,intersect(ind_noadapt,ind_dirplaid)),1))];
    all_dir = [all_dir iDir.*ones(size(intersect(ind_noadapt,ind_diralone)))];
    all_plaid = [all_plaid iDir.*ones(size(intersect(ind_noadapt,ind_dirplaid)))];
    [h_resp(:,iDir,1), p_resp(:,iDir,1)] = ttest2(resp_cell_noadapt{iDir,1},base_cell_noadapt{iDir,1},'dim',2,'tail','right','alpha', 0.05./nStim);
    [h_resp(:,iDir,2), p_resp(:,iDir,2)] = ttest2(resp_cell_noadapt{iDir,2},base_cell_noadapt{iDir,2},'dim',2,'tail','right','alpha', 0.05./nStim);
    trialsperstim_adapt(iDir,1) = length(intersect(ind_adapt,ind_diralone));
    trialsperstim_adapt(iDir,2) = length(intersect(ind_adapt,ind_dirplaid));
    trialInd_adapt{iDir,1} = intersect(ind_adapt,ind_diralone);
    trialInd_adapt{iDir,2} = intersect(ind_adapt,ind_dirplaid);
    resp_cell_adapt{iDir,1} = squeeze(mean(data_dfof_tc(resp_win,:,intersect(ind_adapt,ind_diralone)),1));
    resp_cell_adapt{iDir,2} = squeeze(mean(data_dfof_tc(resp_win,:,intersect(ind_adapt,ind_dirplaid)),1));
    base_cell_adapt{iDir,1} = squeeze(mean(data_dfof_tc(base_win,:,intersect(ind_adapt,ind_diralone)),1));
    base_cell_adapt{iDir,2} = squeeze(mean(data_dfof_tc(base_win,:,intersect(ind_adapt,ind_dirplaid)),1));
    avg_resp_dir(:,iDir,1,2,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,intersect(ind_adapt,ind_diralone)),1),3,'omitnan'));
    avg_resp_dir(:,iDir,1,2,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,intersect(ind_adapt,ind_diralone)),1),[],3,'omitnan')./sqrt(length(intersect(ind_adapt,ind_diralone))));
    avg_resp_dir(:,iDir,2,2,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,intersect(ind_adapt,ind_dirplaid)),1),3,'omitnan'));
    avg_resp_dir(:,iDir,2,2,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,intersect(ind_adapt,ind_dirplaid)),1),[],3,'omitnan')./sqrt(length(intersect(ind_adapt,ind_dirplaid))));
    data_dfof_dir_tc_avg_adapt(:,:,iDir,1) = squeeze(mean(data_dfof_tc(:,:,intersect(ind_adapt,ind_diralone)),3,'omitnan'));
    data_dfof_dir_tc_avg_adapt(:,:,iDir,2) = squeeze(mean(data_dfof_tc(:,:,intersect(ind_adapt,ind_dirplaid)),3,'omitnan'));
end

%only measuring significance in no adapt condition
resp_ind = find(sum(sum(h_resp,2),3));
resp_ind_dir = find(sum(h_resp(:,:,1),2));
resp_ind_plaid = find(sum(h_resp(:,:,2),2));
p_anova_dir = zeros(1,nCells);
p_anova_plaid = zeros(1,nCells);
for iCell = 1:nCells
    p_anova_dir(iCell) = anova1(all_resp_dir(iCell,:), all_dir, 'off');
    p_anova_plaid(iCell) = anova1(all_resp_plaid(iCell,:), all_plaid, 'off');
end

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']), 'data_dfof_dir_tc_avg_adapt', 'resp_cell_adapt','base_cell_adapt','data_dfof_dir_tc_avg_noadapt', 'resp_cell_noadapt','base_cell_noadapt', 'data_dfof_tc', 'avg_resp_dir', 'h_resp','resp_ind', 'tt', 'frame_rate');
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']), 'prewin_frames', 'postwin_frames', 'resp_win', 'base_win', 'trialsperstim_adapt','trialsperstim_noadapt','trialInd_adapt','trialInd_noadapt');

%%
avg_resp_dir_shift = avg_resp_dir;
avg_resp_dir_shift(:,:,2,:,:) = circshift(avg_resp_dir_shift(:,:,2,:,:),2,2);
[maxVal maxDir] = max(avg_resp_dir(:,:,1,1,1),[],2);

for i = 1:5
    ind_orth = intersect(resp_ind_dir, find(maxDir==i));
    figure;
    suptitle(['Ori = ' num2str(stimDirs(i)) '; n = ' num2str(length(ind_orth))])
    subplot(2,2,1)
    plot(mean(data_dfof_dir_tc_avg_noadapt(:,ind_orth,1,1),2))
    hold on 
    plot(mean(data_dfof_dir_tc_avg_adapt(:,ind_orth,1,1),2))
    ylim([-0.2 0.5])
    title('0 deg grating')
    subplot(2,2,2)
    plot(mean(data_dfof_dir_tc_avg_noadapt(:,ind_orth,4,1),2))
    hold on 
    plot(mean(data_dfof_dir_tc_avg_adapt(:,ind_orth,4,1),2))
    ylim([-0.2 0.5])
    title('90 deg grating')
    subplot(2,2,3)
    plot(mean(data_dfof_dir_tc_avg_noadapt(:,ind_orth,1,2),2))
    hold on 
    plot(mean(data_dfof_dir_tc_avg_adapt(:,ind_orth,1,2),2))
    ylim([-0.2 0.5])
    title('60 deg plaid- 0 deg grating')
    subplot(2,2,4)
    plot(mean(data_dfof_dir_tc_avg_noadapt(:,ind_orth,4,2),2))
    hold on 
    plot(mean(data_dfof_dir_tc_avg_adapt(:,ind_orth,4,2),2))
    ylim([-0.2 0.5])
    title('150 deg plaid- 90 deg grating')
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respTCs_pref' num2str(stimDirs(i)) '.pdf']), '-dpdf', '-bestfit')
end

figure;
start = 0;
for i = 1:5
ind_orth = intersect(resp_ind_dir, find(maxDir==i));
subplot(5,2,1+start)
title('Gratings')
polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir(ind_orth,:,1,1,1),1) mean(avg_resp_dir(ind_orth,1,1,1,1),1)])
hold on
polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir(ind_orth,:,1,2,1),1) mean(avg_resp_dir(ind_orth,1,1,2,1),1)])
subplot(5,2,2+start)
title('Plaids')
polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind_orth,:,2,1,1),1) mean(avg_resp_dir_shift(ind_orth,1,2,1,1),1)])
hold on
polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind_orth,:,2,2,1),1) mean(avg_resp_dir_shift(ind_orth,1,2,2,1),1)])
start = start+2;
end
suptitle([mouse ' ' date])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgTuningCurves.pdf']), '-dpdf', '-bestfit')

figure;
ind_orth = intersect(resp_ind_dir, find(maxDir==5));
for i = 1:nStimDir
    subplot(3,4,i)
    plot(mean(data_dfof_dir_tc_avg_noadapt(:,ind_orth,i,1),2))
    hold on 
    plot(mean(data_dfof_dir_tc_avg_adapt(:,ind_orth,i,1),2))
    ylim([-0.2 0.5])
end

figure;
ind_orth = intersect(resp_ind_dir, find(maxDir==5));
for i = 1:nStimDir
    subplot(3,4,i)
    plot(mean(data_dfof_dir_tc_avg_noadapt(:,ind_orth,i,2),2))
    hold on 
    plot(mean(data_dfof_dir_tc_avg_adapt(:,ind_orth,i,2),2))
    ylim([-0.2 0.5])
end


%% pattern selectivity

maskDiff_all = celleqel2mat_padded(input.tMaskTwoGratingDirectionDeg) - celleqel2mat_padded(input.tStimTwoGratingDirectionDeg);
maskDiffs = unique(maskDiff_all);
nCells = size(avg_resp_dir,1);

int = unique(diff(stimDirs));
component = squeeze(avg_resp_dir(:,:,1,:,1)+circshift(avg_resp_dir(:,:,1,:,1),-maskDiffs./int,2));
pattern = squeeze(circshift(avg_resp_dir(:,:,1,:,1),-(maskDiffs/2)./int,2));

comp_corr = zeros(2,nCells);
patt_corr = zeros(2,nCells);
comp_patt_corr = zeros(2,nCells);

for iCell = 1:nCells
    for i = 1:2
        comp_corr(i,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,2,i,1),component(iCell,:,i)));
        patt_corr(i,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,2,i,1),pattern(iCell,:,i)));
        comp_patt_corr(i,iCell) = triu2vec(corrcoef(component(iCell,:,i),pattern(iCell,:,i)));
    end
end

Rp = ((patt_corr)-(comp_corr.*comp_patt_corr))./sqrt((1-comp_corr.^2).*(1-comp_patt_corr.^2));
Rc = ((comp_corr)-(patt_corr.*comp_patt_corr))./sqrt((1-patt_corr.^2).*(1-comp_patt_corr.^2));
Zp = (0.5.*log((1+Rp)./(1-Rp)))./sqrt(1./(nStimDir-3));
Zc = (0.5.*log((1+Rc)./(1-Rc)))./sqrt(1./(nStimDir-3));

Zp_use = intersect(resp_ind_dir, intersect(find(Zp(1,:)>1.28), find(Zp(1,:)-Zc(1,:)>1.28)));
Zc_use = intersect(resp_ind_dir, intersect(find(Zc(1,:)>1.28), find(Zc(1,:)-Zp(1,:)>1.28)));

figure;
subplot(3,2,1)
scatter(Zc(1,resp_ind_dir),Zp(1,resp_ind_dir))
hold on
plotZcZpBorders
xlim([-4 8])
ylim([-4 8])
xlabel('Zc')
ylabel('Zp')
title('control')
axis square
subplot(3,2,2)
scatter(Zc(2,resp_ind_dir),Zp(2,resp_ind_dir))
hold on
plotZcZpBorders
xlim([-4 8])
ylim([-4 8])
xlabel('Zc')
ylabel('Zp')
title('adapt')
axis square
subplot(3,2,3)
ind_orth = intersect(resp_ind_dir,find(maxDir==5));
ind_use = intersect(ind_orth,Zc_use);
polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind_use,:,1,1,1),1) mean(avg_resp_dir_shift(ind_use,1,1,1,1),1)])
hold on
polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind_use,:,1,2,1),1) mean(avg_resp_dir_shift(ind_use,1,1,2,1),1)])
title('Gratings- 90 deg- Zc cells')
subplot(3,2,4)
ind_use = intersect(ind_orth,Zp_use);
polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind_use,:,1,1,1),1) mean(avg_resp_dir_shift(ind_use,1,1,1,1),1)])
hold on
polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind_use,:,1,2,1),1) mean(avg_resp_dir_shift(ind_use,1,1,2,1),1)])
title('Gratings- 90 deg- Zp cells')
subplot(3,2,5)
ind_use = intersect(ind_orth,Zc_use);
polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind_use,:,2,1,1),1) mean(avg_resp_dir_shift(ind_use,1,2,1,1),1)])
hold on
polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind_use,:,2,2,1),1) mean(avg_resp_dir_shift(ind_use,1,2,2,1),1)])
title('Plaids- 90 deg- Zc cells')
subplot(3,2,6)
ind_use = intersect(ind_orth,Zp_use);
polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind_use,:,2,1,1),1) mean(avg_resp_dir_shift(ind_use,1,2,1,1),1)])
hold on
polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind_use,:,2,2,1),1) mean(avg_resp_dir_shift(ind_use,1,2,2,1),1)])
title('Plaids- 90 deg- Zp cells')

stim_use = [3 11];
stim_use_str = num2str(stimDirs(stim_use));

ind_use = intersect(resp_ind_dir, find(ismember(maxDir,stim_use)));
[h p_c] = ttest(Zc(1,ind_use),Zc(2,ind_use));
[h p_z] = ttest(Zp(1,ind_use),Zp(2,ind_use));

adapt_corr = zeros(2,nCells);
for iCell = 1:nCells
    adapt_corr(1,iCell) = triu2vec(corrcoef(avg_resp_dir_shift(iCell,:,1,1,1),avg_resp_dir_shift(iCell,:,1,2,1)));
    adapt_corr(2,iCell) = triu2vec(corrcoef(avg_resp_dir_shift(iCell,:,2,1,1),avg_resp_dir_shift(iCell,:,2,2,1)));
end
[h p_a] = ttest(adapt_corr(1,ind_use),adapt_corr(2,ind_use));

ZpZc = Zp-Zc;
[h p_zc] = ttest(ZpZc(1,ind_use),ZpZc(2,ind_use));

figure;
subplot(2,2,1)
scatter(Zc(1,ind_use),Zp(1,ind_use))
hold on
plotZcZpBorders
xlim([-4 8])
ylim([-4 8])
xlabel('Zc')
ylabel('Zp')
title('control')
axis square
subplot(2,2,2)
scatter(Zc(2,ind_use),Zp(2,ind_use))
hold on
plotZcZpBorders
xlim([-4 8])
ylim([-4 8])
xlabel('Zc')
ylabel('Zp')
title('adapt')
axis square
suptitle([mouse ' ' date ' n = ' num2str(length(ind_use)) ' - pref ' stim_use_str])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_ZcZpScatter_' stim_use_str '.pdf']), '-dpdf', '-bestfit')

figure;
subplot(2,2,1)
scatter(Zc(1,ind_use),Zc(2,ind_use));
xlabel('Control')
ylabel('Adapt')
xlim([-4 4])
ylim([-4 4])
axis square
refline(1)
title(['Zc- p=' num2str(chop(p_c,2))])
subplot(2,2,2)
scatter(Zp(1,ind_use),Zp(2,ind_use));
xlabel('Control')
ylabel('Adapt')
xlim([-4 4])
ylim([-4 4])
axis square
refline(1)
title(['Zp- p=' num2str(chop(p_z,2))])
subplot(2,2,3)
scatter(ZpZc(1,ind_use),ZpZc(2,ind_use))
xlabel('Control')
ylabel('Adapt')
xlim([-4 4])
ylim([-4 4])
axis square
refline(1)
title(['Zp-Zc- p=' num2str(chop(p_zc,2))])
subplot(2,2,4)
scatter(adapt_corr(1,ind_use),adapt_corr(2,ind_use),[],ZpZc(1,ind_use))
xlabel('Grating')
ylabel('Plaid')
title(['Correlation- p=' num2str(chop(p_a,2))])
xlim([-1 1])
ylim([-1 1])
axis square
refline(1)
colorbar
suptitle([mouse ' ' date ' n = ' num2str(length(ind_use)) ' - pref ' stim_use_str])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_ZcZpQuant_' stim_use_str '.pdf']), '-dpdf', '-bestfit')

figure;
for i = 1:nStimDir
    ind = intersect(resp_ind_dir, find(ismember(maxDir,i)));
    subplot(3,4,i)
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind,:,1,1,1),1) mean(avg_resp_dir_shift(ind,1,1,1,1),1)])
    hold on
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind,:,2,1,1),1) mean(avg_resp_dir_shift(ind,1,2,1,1),1)])
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind,:,1,2,1),1) mean(avg_resp_dir_shift(ind,1,1,2,1),1)])        
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind,:,2,2,1),1) mean(avg_resp_dir_shift(ind,1,2,2,1),1)])
    title(num2str(stimDirs(i)))
end
suptitle([mouse ' ' date])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_allPolarAllPref.pdf']), '-dpdf', '-bestfit')


figure; 
[n n2] = subplotn(length(ind_use));
for i = 1:length(ind_use)
    subplot(n,n2,i)
    ind = ind_use(i);
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind,:,1,1,1),1) mean(avg_resp_dir_shift(ind,1,1,1,1),1)])
    hold on
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind,:,2,1,1),1) mean(avg_resp_dir_shift(ind,1,2,1,1),1)])
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind,:,1,2,1),1) mean(avg_resp_dir_shift(ind,1,1,2,1),1)])        
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind,:,2,2,1),1) mean(avg_resp_dir_shift(ind,1,2,2,1),1)])        
    title([num2str(chop(Zp(1,ind)-Zc(1,ind),2)) ' ' num2str(chop(Zp(2,ind)-Zc(2,ind),2))])
end

figure; 
component_shift = circshift(component,2,2);
[n n2] = subplotn(length(ind_use));
for i = 1:length(ind_use)
    subplot(n,n2,i)
    ind = ind_use(i);
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(component_shift(ind,:,1),1) mean(component_shift(ind,1,1),1)])
    hold on
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind,:,2,1,1),1) mean(avg_resp_dir_shift(ind,1,2,1,1),1)])
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(component_shift(ind,:,2),1) mean(component_shift(ind,1,2),1)])        
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind,:,2,2,1),1) mean(avg_resp_dir_shift(ind,1,2,2,1),1)])        
    title([num2str(chop(Zp(1,ind)-Zc(1,ind),2)) ' ' num2str(chop(Zp(2,ind)-Zc(2,ind),2))])
end

%% cross ori interactions
ind_orth = intersect(find(h_resp(:,5,1)), find(avg_resp_dir(:,5,1,1,1)>avg_resp_dir(:,1,1,1,1)));
xori = squeeze(avg_resp_dir(:,1,2,:,1)-(avg_resp_dir(:,1,1,:,1)+avg_resp_dir(:,5,1,:,1))./(avg_resp_dir(:,1,2,:,1)+avg_resp_dir(:,1,1,:,1)+avg_resp_dir(:,5,1,:,1)));
figure;
scatter(xori(ind_orth,1),xori(ind_orth,2))
xlabel('Control')
ylabel('Adapt')
xlim([-1 1])
ylim([-1 1])
axis square
refline(1)
hline(0)
