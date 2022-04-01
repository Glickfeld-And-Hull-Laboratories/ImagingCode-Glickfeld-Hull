%% get path names
close all 
clear all global
clc

date = '211210';
ImgFolder = strvcat('002');
time = strvcat('1322');
mouse = 'i1367';
doFromRef = 0;
ref = strvcat('002');
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

%% load
tic
data = [];
CD = [LG_base '\Data\2P_images\' mouse '\' date '\' ImgFolder];
%CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\AW68\two-photon imaging\' date '\' ImgFolder(irun,:)];
%CD = [LG_base '\Data\2P_images\' mouse '-KM' '\' date '_' mouse '\' ImgFolder(irun,:)];
%CD = ['\\CRASH.dhe.duke.edu\data\home\kevin\Data\2P\' date '_' mouse '\' ImgFolder(irun,:)];
cd(CD);
imgMatFile = [ImgFolder '_000_000.mat'];
load(imgMatFile);
fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time '.mat'];
load(fName);
nframes = [input.counterValues{end}(end) info.config.frames];
fprintf(['Reading run ' ImgFolder '- ' num2str(min(nframes)) ' frames \r\n'])
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

data_avg = mean(data(:,:,40001:40500),3);
%% Register data

if exist(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(double(data),[],[],out);
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
stimOneContrast = celleqel2mat_padded(input.tStimOneGratingContrast);
stimOneCons = unique(stimOneContrast);

ori_mat = celleqel2mat_padded(input.tStimOneGratingDirectionDeg);
oris = unique(ori_mat);
nOri = length(oris);
stimOneTime = celleqel2mat_padded(input.tStimOneGratingOnTimeMs);
stimOneTimes = unique(stimOneTime);
nTime = length(stimOneTimes);
isiTime = celleqel2mat_padded(input.tISITimeMs);
isiTimes = unique(isiTime);
nISI = length(isiTimes);
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

zeroConInd = find(stimOneContrast==0);

if input.doRandStimOnTime
    data_dfof_stim = zeros(sz(1),sz(2),nTime+1);
    [n n2] = subplotn(nTime+1);
    figure;
    for it = 1:nTime
        if length(zeroConInd)>0 
            ind = setdiff(find(stimOneTime == stimOneTimes(it)),zeroConInd);
        else
            ind = find(stimOneTime == stimOneTimes(it));
        end
        data_dfof_stim(:,:,it) = nanmean(data_one_dfof(:,:,ind),3);
        subplot(n,n2,it)
        imagesc(data_dfof_stim(:,:,it))
        title(num2str(stimOneTimes(it)))
    end
elseif input.doRandDir
    data_dfof_stim = zeros(sz(1),sz(2),nOri+1);
    [n n2] = subplotn(nOri+1);
    figure;
    for it = 1:nOri
        if length(zeroConInd)>0 
            ind = setdiff(find(ori_mat == oris(it)),zeroConInd);
        else
            ind = find(ori_mat == oris(it));
        end
        data_dfof_stim(:,:,it) = nanmean(data_one_dfof(:,:,ind),3);
        subplot(n,n2,it)
        imagesc(data_dfof_stim(:,:,it))
        title(num2str(oris(it)))
    end
elseif input.doRandISITime
    data_dfof_stim = zeros(sz(1),sz(2),nISI+1);
    [n n2] = subplotn(nISI+1);
    figure;
    for it = 1:nISI
        if length(zeroConInd)>0 
            ind = setdiff(find(isiTime == isiTimes(it)),zeroConInd);
        else
            ind = find(isiTime == isiTimes(it));
        end
        data_dfof_stim(:,:,it) = nanmean(data_one_dfof(:,:,ind),3);
        subplot(n,n2,it)
        imagesc(data_dfof_stim(:,:,it))
        title(num2str(isiTimes(it)))
    end
    
end
subplot(n,n2,it+1)
if length(zeroConInd)>0 
    ind = setdiff(1:nTrials,zeroConInd);
else
    ind = 1:nTrials;
end
data_dfof_stim(:,:,it+1) = nanmean(data_one_dfof(:,:,ind),3);
imagesc(data_dfof_stim(:,:,it+1))
title('Avg all')

data_dfof = cat(3,data_dfof_stim,max(data_dfof_stim,[],3));

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
tISITime = celleqel2mat_padded(input.tISITimeMs);
ISIs = unique(tISITime);
nISI = length(ISIs);

tISIframes = cStimTwo-cStimOne;
ISIf = unique(tISIframes);


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

tt1 = [-20:79].*(1000/frame_rate);
tt2 = [-20:79].*(1000/frame_rate);
base_win = [1:20];
resp_win = [26:45];

for iOri = 1:nOri
    ind = find(ori_mat == oris(iOri));
    [h1_ori(:,iOri) p1_ori(:,iOri)] = ttest(squeeze(mean(tc_one_dfof(resp_win,:,ind),1)),squeeze(mean(tc_one_dfof(base_win,:,ind),1)),'tail','right','dim',2);
    [h2_ori(:,iOri) p2_ori(:,iOri)] = ttest(squeeze(mean(tc_two_dfof(resp_win,:,ind),1)),squeeze(mean(tc_two_dfof(base_win,:,ind),1)),'tail','right','dim',2);
end
[h1 p1] = ttest(squeeze(mean(tc_one_dfof(resp_win,:,:),1)),squeeze(mean(tc_one_dfof(base_win,:,:),1)),'tail','right','dim',2);
[h2 p2] = ttest(squeeze(mean(tc_two_dfof(resp_win,:,:),1)),squeeze(mean(tc_two_dfof(base_win,:,:),1)),'tail','right','dim',2);

figure;
tc_one_dfof_all = zeros(100,nCells,nISI);
tc_two_dfof_all = zeros(100,nCells,nISI);
subplot(2,2,1)
shadedErrorBar(tt1,squeeze(nanmean(nanmean(tc_one_dfof(:,:,:),3),2)),squeeze(nanstd(nanmean(tc_one_dfof(:,:,:),3),[],2))./sqrt(5));%-mean(tc_one_dfof_all(base_win,:,it),1),2)))
ylim([-0.01 0.3])
subplot(2,2,2)
for it = 1:nOri
    ind = setdiff(find(ori_mat == oris(it)),zeroConInd);
    plot(tt1,squeeze(nanmean(nanmean(tc_one_dfof(:,:,ind),3),2)));%-mean(tc_one_dfof_all(base_win,:,it),1),2)))
    hold on
end
ylim([-0.01 0.3])
legend(num2str(oris'),'location','southeast')
for it = 1:nISI
    ind = setdiff(find(tISITime == ISIs(it)),zeroConInd);
    tc_one_dfof_all(:,:,it) = nanmean(tc_one_dfof(:,:,ind),3);
    tc_two_dfof_all(:,:,it) = nanmean(tc_two_dfof(:,:,ind),3);
    subplot(2,2,3)
    plot(tt1,squeeze(mean(tc_one_dfof_all(:,:,it),2)));
    ylim([-0.01 0.3])
    hold on
    subplot(2,2,4)
    plot(tt2,squeeze(mean(tc_two_dfof_all(:,:,it),2)));%-mean(tc_two_dfof_all(base_win,:,it),1),2)))
    hold on
    ylim([-0.01 0.3])
end
legend(num2str(ISIs'),'location','southeast')
tc_two_dfof_all(:,:,it+1) = nanmean(tc_two_dfof(:,:,zeroConInd),3);
plot(tt2,squeeze(mean(tc_two_dfof_all(:,:,it+1),2)))
xlabel('Time (ms)')
ylabel('dF/F')
title('Stim two')
subplot(2,2,1)
xlabel('Time (ms)')
ylabel('dF/F')
title('Stim one')

figure;
[n n2] = subplotn(nCells);
for iC = 1:nCells
    subplot(n,n2,iC)
    for i = 1:nISI
        ind = setdiff(find(tISITime == ISIs(i)),zeroConInd);
        plot(tt1,nanmean(tc_one_dfof(:,iC,ind),3))
        hold on
    end
end
suptitle('ISI')

figure;
[n n2] = subplotn(nCells);
for iC = 1:nCells
    subplot(n,n2,iC)
    for i = 1:nOri
        ind = setdiff(find(ori_mat == oris(i)),zeroConInd);
        plot(tt1,nanmean(tc_one_dfof(:,iC,ind),3))
        hold on
    end
end
suptitle('Ori')

figure;
[n n2] = subplotn(nCells);
for iC = 1:10
    figure;
    for i = 1:nISI
        subplot(2,3,i)
        ind1 = find(tISITime == ISIs(i));
        for is = 1:nOri
            ind = intersect(ind1,setdiff(find(ori_mat == oris(is)),zeroConInd));
            plot(tt1,nanmean(tc_one_dfof(:,iC,ind),3))
            hold on
        end
        if i<nISI
        vline(ISIs(i))
        end
    end
end
suptitle('Ori')

figure;
for i = 1:nISI
    subplot(3,2,i)
    ind = setdiff(find(tISITime == ISIs(i)),zeroConInd);
    shadedErrorBar(tt1,nanmean(nanmean(tc_one_dfof(:,:,ind),3),2),nanstd(nanmean(tc_one_dfof(:,:,ind),3),[],2)./sqrt(nCells));
    ylim([-0.05 0.25])
    vline(ISIs(i))
end
subplot(3,2,i+1)
shadedErrorBar(tt1,nanmean(nanmean(tc_one_dfof(:,:,:),3),2),nanstd(nanmean(tc_one_dfof(:,:,:),3),[],2)./sqrt(nCells));
    ylim([-0.05 0.25])
    
for iCell = 1:nCells
    figure;
for i = 1:nISI
    subplot(3,2,i)
    ind = setdiff(find(tISITime == ISIs(i)),zeroConInd);
    shadedErrorBar(tt1,nanmean(nanmean(tc_one_dfof(:,iCell,ind),3),2),nanstd(nanmean(tc_one_dfof(:,iCell,ind),2),[],3)./sqrt(length(ind)));
    ylim([-0.05 0.25])
    vline(ISIs(i))
end
subplot(3,2,i+1)
shadedErrorBar(tt1,nanmean(nanmean(tc_one_dfof(:,iCell,:),3),2),nanstd(nanmean(tc_one_dfof(:,iCell,:),2),[],3)./sqrt(length(cStimOne)));
    ylim([-0.05 0.25])   
end

figure;
tc_one_dfof_resp = mean(tc_one_dfof(resp_win,:,:),1)-mean(tc_one_dfof(base_win,:,:),1);
tc_two_dfof_resp = mean(tc_two_dfof(resp_win,:,:),1)-mean(tc_two_dfof(base_win,:,:),1);
tc_two_dfof_resp_stim = zeros(nCells,nTime+1);
for it = 1:nTime
    ind = setdiff(find(stimOneTime == stimOneTimes(it)),zeroConInd);
    subplot(2,2,2)
    errorbar(stimOneTimes(it),squeeze(mean(nanmean(tc_one_dfof_resp(:,:,ind),3),2)), squeeze(std(nanmean(tc_one_dfof_resp(:,:,ind),3),[],2))./sqrt(nCells),'o')
    hold on
    subplot(2,2,4)
    errorbar(stimOneTimes(it),squeeze(mean(nanmean(tc_two_dfof_resp(:,:,ind),3),2)), squeeze(std(nanmean(tc_two_dfof_resp(:,:,ind),3),[],2))./sqrt(nCells),'o')
    hold on
    tc_two_dfof_resp_stim(:,it) = nanmean(tc_two_dfof_resp(:,:,ind),3);
end
errorbar(0,squeeze(mean(nanmean(tc_two_dfof_resp(:,:,zeroConInd),3),2)), squeeze(std(nanmean(tc_two_dfof_resp(:,:,zeroConInd),3),[],2))./sqrt(nCells),'o')
tc_two_dfof_resp_stim(:,it+1) = nanmean(tc_two_dfof_resp(:,:,zeroConInd),3);
xlabel('Stim one time (ms)')
ylabel('dF/F')
title('Stim two')
xlim([-100 4500])
ylim([0 0.15])
subplot(2,2,2)
xlabel('Stim one time (ms)')
ylabel('dF/F')
title('Stim one')
xlim([-100 4500])
ylim([0 0.15])
suptitle([date ' ' mouse ])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_adaptDuration.pdf']),'-dpdf','-fillpage')

tc_two_dfof_resp_norm = tc_two_dfof_resp_stim./tc_two_dfof_resp_stim(:,end);
figure;
for it = 1:nTime
    errorbar(stimOneTimes(it),squeeze(mean(tc_two_dfof_resp_norm(:,it),1)), squeeze(std(tc_two_dfof_resp_norm(:,it),[],1))./sqrt(nCells),'o')
    hold on
end
xlabel('Stim one time (ms)')
ylabel('Normalized dF/F')
ylim([0 1.2])
set(gca,'xscale','log')

%%
base_win = 21:23;
resp_win = 27:29;

resp_dfof_dir = zeros(nCells,nOri,nISI,2);
tr_ind = cell(nOri,nISI);
for i = 1:nOri
    ind = find(ori_mat == oris(i));
    for ii = 1:nISI
        ind_use = intersect(ind, find(tISITime == ISIs(ii)));
        tr_ind{i,ii} = ind_use;
        resp_temp1 = nanmean(tc_one_dfof(:,:,ind_use),3);
        resp_temp2 = nanmean(tc_two_dfof(:,:,ind_use),3);
        resp_dfof_dir(:,i,ii,1) = squeeze(mean(resp_temp1(resp_win,:),1)-mean(resp_temp1(base_win,:),1));
        resp_dfof_dir(:,i,ii,2) = squeeze(mean(resp_temp2(resp_win,:),1)-mean(resp_temp2(base_win,:),1));
    end
end
max_dir = zeros(1,nCells);
max_val = zeros(1,nCells);
h = zeros(1,nCells);
norm_df = zeros(nCells,nISI);
tc_dfof = zeros(size(tc_one_dfof,1),nCells,nISI);
tc_dfof_pref = zeros(size(tc_one_dfof,1),nCells,nISI,2);

for iC = 1:nCells
    [max_val(iC) max_dir(iC)] = max(mean(resp_dfof_dir(iC,:,:,1),3),[],2);
    ind = find(ori_mat == oris(max_dir(iC)));
    h(iC) = ttest(squeeze(mean(tc_one_dfof(resp_win,iC,ind),1)),squeeze(mean(tc_one_dfof(base_win,iC,ind),1)),'tail','right');
    norm_df(iC,:,1) = squeeze(resp_dfof_dir(iC,max_dir(iC),:,2)./resp_dfof_dir(iC,max_dir(iC),:,1));
    norm_df(iC,:,2) = squeeze(resp_dfof_dir(iC,max_dir(iC),:,2)./mean(resp_dfof_dir(iC,max_dir(iC),:,1),3));
    for i = 1:nISI
        tc_dfof(:,iC,i,1) = squeeze((nanmean(tc_two_dfof(:,iC,tr_ind{max_dir(iC),i}),3)-nanmean(mean(tc_two_dfof(base_win,iC,tr_ind{max_dir(iC),i}),1),3))./resp_dfof_dir(iC,max_dir(iC),i,1));
        tc_dfof(:,iC,i,2) = squeeze((nanmean(tc_two_dfof(:,iC,tr_ind{max_dir(iC),i}),3)-nanmean(mean(tc_two_dfof(base_win,iC,tr_ind{max_dir(iC),i}),1),3))./mean(resp_dfof_dir(iC,max_dir(iC),:,1),3));
        tc_dfof_pref(:,iC,i,1) = nanmean(tc_one_dfof(:,iC,tr_ind{max_dir(iC),i}),3);
        tc_dfof_pref(:,iC,i,2) = nanmean(tc_two_dfof(:,iC,tr_ind{max_dir(iC),i}),3)-nanmean(mean(tc_two_dfof(base_win,iC,tr_ind{max_dir(iC),i}),1),3);
    end
end
resp_ind= find(h & max_val>0.05);
figure; 
subplot(2,2,1)
for i = 1:nISI
    plot(tt2,mean(tc_dfof(:,resp_ind,i,1),2))
    hold on
end
ylim([-0.1 1.2])
xlabel('Time (ms)')
ylabel('Normalized dF/F')
title('Matched')
subplot(2,2,2)
for i = 1:nISI
    plot(tt2,mean(tc_dfof(:,resp_ind,i,2),2))
    hold on
end
ylim([-0.1 1.2])
xlabel('Time (ms)')
ylabel('Normalized dF/F')
title('All')
subplot(2,2,3)
errorbar(ISIs,mean(norm_df(resp_ind,:,1),1),std(norm_df(resp_ind,:,1),[],1)./sqrt(length(resp_ind)),'o');
xlabel('ISI (ms)')
ylabel('Normalized dF/F')
ylim([0 1.2])
subplot(2,2,4)
errorbar(ISIs,mean(norm_df(resp_ind,:,2),1),std(norm_df(resp_ind,:,2),[],1)./sqrt(length(resp_ind)),'o');
xlabel('ISI (ms)')
ylabel('Normalized dF/F')
ylim([0 1.2])

sgtitle([date ' ' mouse '- n=' num2str(length(resp_ind)) ' cells'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_adaptDuration.pdf']),'-dpdf','-fillpage')
    

figure;
for i= 1:nISI
    subplot(2,3,i)
    plot(tt2,mean(mean(tc_dfof_pref(:,resp_ind,i,1),3),2)); 
    hold on; 
    plot(tt2,mean(tc_dfof_pref(:,resp_ind,i,2),2));
    ylim([-0.02 .2])
end


%%

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']), 'tc_dfof','resp_dfof_dir', 'max_dir', 'norm_df', 'tc_dfof_pref',  'tc_one_dfof', 'tc_two_dfof', 'nCells','frame_rate', 'resp_ind')
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']),'cStimOne','cStimTwo','ori_mat','nOri','oris','tISITime','ISIs','nISI','stimOneTime', 'stimOneTimes', 'nTime', 'stimOneContrast', 'stimOneCons','base_win','resp_win')

%%

calib = 1/26.6; %mm per pixel

% Load and combine eye tracking data
nrun = length(ImgFolder);
data = [];
for irun =  1:nrun
    fn = [ImgFolder(irun,:) '_000_000_eye.mat'];

    data_temp = load(fn);          % should be a '*_eye.mat' file
    data_temp = squeeze(data_temp.data);

    nFrames = input.counterValues{end}(end);
    data = cat(3, data, data_temp(:,:,1:nFrames));      % the raw images...
end
figure;
data_avg = mean(data,3);
imagesc(data_avg);
movegui('center')
ax = gca;
rect = getrect(ax);
datat = data_avg(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));
figure;
imagesc(datat)
movegui('center')
while 1  % till broken out of

    % interactively get clicks
    [X Y selectionType] = getAPoint(gca);

    if isnan(X)
        key = lower(Y);
        switch key
          case char(13) % return
            break;  % out of loop, done
          case 'z' 
            imagesc(datat)
            rect = getrect(ax);
            datat = data(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));
            imagesc(datat)
        end
        continue
    end
end
close all
data = data(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3),:);

%%
rad_range = [4 15];
warning off;
A = cell(size(data,3),1);
B = cell(size(data,3),1);
C = cell(size(data,3),1);
D = cell(size(data,3),1);
for n = 1:size(data,3)
    A{n} = [0,0];
    B{n} = [0];
    C{n} = [0];
    D{n} = [0];
end
eye = struct('Centroid',A,'Area',B,'Val',C,'SNR',D);
radii = [];
for n = 1:size(data,3)
    [center,radii,metric] = imfindcircles(squeeze(data(:,:,n)),rad_range,'Sensitivity',0.95);
              % pick the circle with best score
    if(isempty(center))
        eye(n).Centroid = [NaN NaN];    % could not find anything...
        eye(n).Area = NaN;
        eye(n).Val = NaN;
        eye(n).SNR = NaN;
    else
        snr = zeros(1,size(center,1));
        for idx = 1:size(center,1)
            t = double(data(:,:,n));
            vector_of_y_values = (1:size(data,1)) - center(idx,2);
            vector_of_x_values = (1:size(data,2)) - center(idx,1);
            [Yg, Xg] = ndgrid(vector_of_y_values, vector_of_x_values);
            idx1 = find(Xg.^2 + Yg.^2 < (radii(idx)/2).^2);
            idx2 = find(Xg.^2 + Yg.^2 < (radii(idx).*2.5).^2 & Xg.^2 + Yg.^2 > (radii(idx).*1.5).^2);
            snr(idx) = mean(t(idx1))./mean(t(idx2));
        end
        [v,idx] = max(snr);
        val = metric(idx);
        t = double(data(:,:,n));
        vector_of_y_values = (1:size(data,1)) - center(idx,2);
        vector_of_x_values = (1:size(data,2)) - center(idx,1);
        [Yg, Xg] = ndgrid(vector_of_y_values, vector_of_x_values);
        idx1 = find(Xg.^2 + Yg.^2 < (radii(idx)/2).^2);
        idx2 = find(Xg.^2 + Yg.^2 < (radii(idx).*2.5).^2 & Xg.^2 + Yg.^2 > (radii(idx).*1.5).^2);
        snr = mean(t(idx1))./mean(t(idx2));
        eye(n).SNR = snr;
        eye(n).Val = val;
        eye(n).Centroid = center(idx,:);
        eye(n).Area = pi*radii(idx)^2;
    end
    if mod(n,100)==0
        fprintf('Frame %d/%d\n',n,size(data,3));
    end
end
Centroid = cell2mat({eye.Centroid}');
Area = cell2mat({eye.Area}');
Val = double(cell2mat({eye.Val}'));
SNR = double(cell2mat({eye.SNR}'));
Eye_data = data;
nanframes = length(find(isnan(Area)));

% no measurement frames
figure; 
subplot(2,2,1)
hist(sqrt(Area./pi));
xlabel('radius')
subplot(2,2,2)
hist(SNR);
xlabel('SNR')
subplot(2,2,3)
hist(Val);
xlabel('Metric')
movegui('center')

x1 = find(isnan(Area));
x2 = find(~isnan(Area));
x3 = unique([find(Val<0.1); find(Val<0.20 & SNR<1.7)]);

x = unique([x1; x3]);
if length(x)>25
    minx = 25;
else
    minx = length(x);
end

frames = sort(randsample(length(x),minx));
figure;
start = 1;
for i = 1:minx
    subplot(5,5,start);
    imagesq(data(:,:,x(frames(i)))); 
    hold on;
    scatter(Centroid(x(frames(i)),1), Centroid(x(frames(i)),2))
    title([num2str(chop(SNR(x(frames(i))),2)) ' ' num2str(chop(Val(x(frames(i))),2))])
    %title(num2str(x(frames(i))))
    start = start+1;
end
movegui('center')
sgtitle(['No pupil detected- ' num2str(length(x)) ' frames'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_noPupil2.pdf']),'-dpdf','-fillpage');

x = setdiff(x2,x3);
if length(x)>25
    minx = 25;
else
    minx = length(x);
end
frames = sort(randsample(length(x),minx));
figure;
start = 1;
for i = 1:minx
    subplot(5,5,start);
    imagesq(data(:,:,x(frames(i)))); 
    hold on;
    scatter(Centroid(x(frames(i)),1), Centroid(x(frames(i)),2))
    title([num2str(chop(SNR(x(frames(i))),2)) ' ' num2str(chop(Val(x(frames(i))),2))])
    %title(num2str(x(frames(i))))
    start = start+1;
end
movegui('center')
sgtitle('Pupil detected')
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Pupil.pdf']),'-dpdf','-fillpage');
%%
tc_pupil_one = nan(100,nTrials);
tc_pupil_two = nan(100,nTrials);
   
for itrial = 1:nTrials
    if ~isnan(cStimOne(itrial)) & ~isnan(cStimTwo(itrial)) & (cStimOne(itrial)+79)<sz(3) & (cStimTwo(itrial)+19)<sz(3)
        tc_pupil_one(:,itrial) = Area(cStimOne(itrial)-20:cStimOne(itrial)+79);
        tc_pupil_two(:,itrial) = Area(cStimTwo(itrial)-20:cStimTwo(itrial)+79);
    end
end

tc_pupil = zeros(size(tc_pupil_one,1),nISI);
prestim2_pupil = zeros(nISI,2);
for i = 1:nISI
    ind_isi = find(tISITime == ISIs(i));
    tc_pupil(:,i) = squeeze(nanmean(tc_pupil_one(:,ind_isi),2));
    prestim2_pupil(i,1) = nanmean(mean(tc_pupil_two(1:20,ind_isi),1),2);
    prestim2_pupil(i,2) = nanstd(mean(tc_pupil_two(1:20,ind_isi),1),[],2)./sqrt(length(ind_isi));
end

figure;
subplot(3,1,1)
plot(tt2,tc_pupil)
subplot(3,1,2)
shadedErrorBar(tt2, nanmean(tc_pupil_one,2),nanstd(tc_pupil_one,[],2)./sqrt(nTrials))
vline(ISIs(1:4))
subplot(3,1,3)
errorbar(ISIs,prestim2_pupil(:,1),prestim2_pupil(:,2),'ok')
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_PupilTC.pdf']),'-dpdf','-fillpage');

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']), 'rect', 'Area', 'Centroid', 'SNR', 'Val', 'frame_rate' , 'tc_pupil_one','tc_pupil_two');
