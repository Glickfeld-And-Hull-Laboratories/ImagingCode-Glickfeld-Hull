%% get path names
close all 
clear all global
clc
date = '220818';
ImgFolder = {'002'};
time = strvcat('1116');
mouse = 'i1365';
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
data_avg = mean(data(:,:,60001:60500),3);
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
sf_mat = celleqel2mat_padded(input.tStimOneGratingSpatialFreqCPD);
sfs = unique(sf_mat);
nSF = length(sfs);
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
elseif input.doRandDir & input.doRandSF
    data_dfof_stim = zeros(sz(1),sz(2),nOri*nSF);
    figure;
    it = 1;
    n= nOri;
    n2 = nSF;
    for io = 1:nOri
        ind = find(ori_mat == oris(io));
        for is = 1:nSF
            ind_use = intersect(ind,find(sf_mat == sfs(is)));
            subplot(n,n2,it)
            data_dfof_stim(:,:,it) = nanmean(data_one_dfof(:,:,ind_use),3);
            imagesc(data_dfof_stim(:,:,it))
            it=it+1;
        end
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
elseif input.doRandSF
    data_dfof_stim = zeros(sz(1),sz(2),nSF);
    figure;
    it = 1;
    [n n2] = subplotn(nSF+1);
    for it = 1:nSF
        ind_use = find(sf_mat == sfs(it));
        subplot(n,n2,it)
        data_dfof_stim(:,:,it) = nanmean(data_one_dfof(:,:,ind_use),3);
        imagesc(data_dfof_stim(:,:,it))
        title(num2str(chop(sfs(it),2)))
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

data_dfof = cat(3,mean(data_f,3,'omitnan'),cat(3,data_dfof_stim,max(data_dfof_stim,[],3)));
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

base_win = 21:23;
resp_win = 28:30;
figure;
subplot(2,1,1)
shadedErrorBar(1:100,squeeze(nanmean(nanmean(tc_one_dfof(:,:,:),3),2)),squeeze(nanstd(nanmean(tc_one_dfof(:,:,:),3),[],2))./sqrt(5));%-mean(tc_one_dfof_all(base_win,:,it),1),2)))
vline([base_win resp_win])
subplot(2,1,2)
shadedErrorBar(1:100,squeeze(nanmean(nanmean(tc_two_dfof(:,:,:),3),2)),squeeze(nanstd(nanmean(tc_two_dfof(:,:,:),3),[],2))./sqrt(5));%-mean(tc_one_dfof_all(base_win,:,it),1),2)))
vline([base_win resp_win])

sf_mat = celleqel2mat_padded(input.tStimOneGratingSpatialFreqCPD);
sfs = unique(sf_mat);
nSF = length(sfs);

h1_ori = zeros(nCells,nOri,nSF);
h2_ori = zeros(nCells,nOri,nSF);
p1_ori = zeros(nCells,nOri,nSF);
p2_ori = zeros(nCells,nOri,nSF);


for iOri = 1:nOri
    ind = find(ori_mat == oris(iOri));
    for iSF = 1:nSF
        ind_use = intersect(ind, find(sf_mat == sfs(iSF)));
        [h1_ori(:,iOri,iSF) p1_ori(:,iOri,iSF)] = ttest(squeeze(mean(tc_one_dfof(resp_win,:,ind_use),1)),squeeze(mean(tc_one_dfof(base_win,:,ind_use),1)),'tail','right','dim',2,'alpha',0.05./(nSF*nOri));
        [h2_ori(:,iOri,iSF) p2_ori(:,iOri,iSF)] = ttest(squeeze(mean(tc_two_dfof(resp_win,:,ind_use),1)),squeeze(mean(tc_two_dfof(base_win,:,ind_use),1)),'tail','right','dim',2,'alpha',0.05./(nSF*nOri));
    end
end
[h1 p1] = ttest(squeeze(mean(tc_one_dfof(resp_win,:,:),1)),squeeze(mean(tc_one_dfof(base_win,:,:),1)),'tail','right','dim',2);
[h2 p2] = ttest(squeeze(mean(tc_two_dfof(resp_win,:,:),1)),squeeze(mean(tc_two_dfof(base_win,:,:),1)),'tail','right','dim',2);

%%
tt1 = [-20:79].*1000/frame_rate;
tt2 = tt1;
figure;
tc_one_dfof_all = zeros(100,nCells,nISI);
tc_two_dfof_all = zeros(100,nCells,nISI);
tc_two_dfof_dir = zeros(100,nCells,nOri);
subplot(2,2,1)
shadedErrorBar(tt1,squeeze(nanmean(nanmean(tc_one_dfof(:,:,:),3),2)),squeeze(nanstd(nanmean(tc_one_dfof(:,:,:),3),[],2))./sqrt(5));%-mean(tc_one_dfof_all(base_win,:,it),1),2)))
ylim([-0.01 0.3])
subplot(2,2,2)
for it = 1:nOri
    ind = setdiff(find(ori_mat == oris(it)),zeroConInd);
    tc_one_dfof_dir(:,:,it) = nanmean(tc_one_dfof(:,:,ind),3);
    plot(tt1,squeeze(nanmean(nanmean(tc_one_dfof(:,:,ind),3),2)));%-mean(tc_one_dfof_all(base_win,:,it),1),2)))
    hold on
end
ylim([-0.01 0.3])
legend(num2str(oris'),'location','southeast')
for it = 1:nISI
    ind = setdiff(find(tISITime == ISIs(it)),zeroConInd);
    tc_one_dfof_all(:,:,it) = mean(tc_one_dfof(:,:,ind),3,'omitnan');
    tc_two_dfof_all(:,:,it) = mean(tc_two_dfof(:,:,ind),3,'omitnan')-(mean(mean(tc_two_dfof(base_win,:,ind),1,'omitnan'),3,'omitnan'));
    subplot(2,2,3)
    plot(tt1,squeeze(mean(tc_one_dfof_all(:,:,it),2)));
    ylim([-0.01 0.5])
    hold on
    subplot(2,2,4)
    plot(tt2,squeeze(mean(tc_two_dfof_all(:,:,it),2)));%-mean(tc_two_dfof_all(base_win,:,it),1),2)))
    hold on
    ylim([-0.01 0.5])
end
legend(num2str(ISIs'),'location','southeast')
tc_two_dfof_all(:,:,it+1) = nanmean(tc_two_dfof(:,:,zeroConInd),3);
plot(tt2,squeeze(mean(tc_two_dfof_all(:,:,it+1),2)-(mean(mean(tc_two_dfof_all(base_win,:,it+1),1,'omitnan'),2,'omitnan'))))
xlabel('Time (ms)')
ylabel('dF/F')
title('Stim two')
subplot(2,2,1)
xlabel('Time (ms)')
ylabel('dF/F')
title('Stim one')
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_adaptSummary.pdf']),'-dpdf','-fillpage')

tc_one_dfof_resp = mean(tc_one_dfof(resp_win,:,:),1)-mean(tc_one_dfof(base_win,:,:),1);
tc_two_dfof_resp = mean(tc_two_dfof(resp_win,:,:),1)-mean(tc_two_dfof(base_win,:,:),1);

if nISI>1
    resp_dfof_norm = zeros(nCells,nISI);
    figure;
    for it = 1:nISI
        ind = setdiff(find(tISITime == ISIs(it)),zeroConInd);
        tc_one_dfof_all(:,:,it) = mean(tc_one_dfof(:,:,ind),3,'omitnan');
        tc_two_dfof_all(:,:,it) = mean(tc_two_dfof(:,:,ind),3,'omitnan')-(mean(mean(tc_two_dfof(base_win,:,ind),1,'omitnan'),3,'omitnan'));
        subplot(2,2,1)
        plot(tt1,squeeze(mean(tc_one_dfof_all(:,:,it),2)));
        hold on
        subplot(2,2,2)
        plot(tt2,squeeze(mean(tc_two_dfof_all(:,:,it),2)));%-mean(tc_two_dfof_all(base_win,:,it),1),2)))
        hold on
        subplot(2,2,3)
        resp_dfof_norm(:,it) = squeeze(mean(tc_two_dfof_resp(:,:,ind),3,'omitnan')./mean(tc_one_dfof_resp(:,:,ind),3,'omitnan'));
        errorbar(ISIs(it), mean(resp_dfof_norm(:,it),1), std(resp_dfof_norm(:,it),[],1)./sqrt(nCells),'ok')
        ylim([0 1.5])
        hold on
    end
    sgtitle([mouse ' ' date])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_adaptISI.pdf']),'-dpdf','-fillpage')
end

if nTime>1
    figure;
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

end
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']), 'tc_one_dfof_all', 'tc_two_dfof_all', 'base_win', 'resp_win', 'tc_one_dfof_resp', 'tc_two_dfof_resp', 'resp_dfof_norm', 'tt1', 'tt2')
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']), 'tISITime', 'ISIs', 'nISI')

%%

wheel_speed = wheelSpeedCalc(input,32,'orange');

%%
if input.doRandDir
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

    figure; 
    resp_ind = find(sum(h1_ori,[2 3]));
    [n n2] = subplotn(length(resp_ind));
    for iC = 1:length(resp_ind)
        iCell = resp_ind(iC);
        subplot(n,n2,iC)
        imagesc(squeeze(resp_dfof_dir(iCell,:,1,1)))
        %title(num2str(sum(sum(h1_ori(iCell,:,:),2),3)>0))
        axis off
        clim([0 0.4])
    end
end

if input.doRandSF & ~input.doRandDir 
    resp_dfof_stim = zeros(nCells,nSF,2);
    resp_dfof_var = zeros(nCells,nSF,2);
    tr_ind = cell(1,nSF);
    for i = 1:nSF
        ind_use = find(sf_mat == sfs(i));
        tr_ind{1,i} = ind_use;
        resp_dfof_stim(:,i,1) = mean(mean(tc_one_dfof(resp_win,:,ind_use),1)-mean(tc_one_dfof(base_win,:,ind_use),1),3,'omitnan')';
        resp_dfof_stim(:,i,2) = mean(mean(tc_two_dfof(resp_win,:,ind_use),1)-mean(tc_two_dfof(base_win,:,ind_use),1),3,'omitnan')';
        resp_dfof_var(:,i,1) = var(mean(tc_one_dfof(resp_win,:,ind_use),1)-mean(tc_one_dfof(base_win,:,ind_use),1),[],3,'omitnan')';
        resp_dfof_var(:,i,2) = var(mean(tc_two_dfof(resp_win,:,ind_use),1)-mean(tc_two_dfof(base_win,:,ind_use),1),[],3,'omitnan')';
    end

    figure; 
    resp_ind = find(sum(h1_ori,[2 3]));
    [n n2] = subplotn(length(resp_ind));
    for iC = 1:length(resp_ind)
        iCell = resp_ind(iC);
        subplot(n,n2,iC)
        plot(sfs,resp_dfof_stim(iCell,:,1),'o')
    end
    resp_dfof_sfnan = resp_dfof_stim(:,:,1);
    resp_dfof_sfnan(squeeze(~h1_ori)) = nan;

    [pref_val pref_sf] = max(resp_dfof_sfnan(:,:,1),[],2);
    norm_dfof_stim = resp_dfof_stim(:,:,2)./resp_dfof_stim(:,:,1);
    norm_dfof_stim_nan = norm_dfof_stim;
    norm_dfof_stim_nan(squeeze(~h1_ori)) = nan;
    norm_dfof_stim_nan(norm_dfof_stim_nan<0.05) = nan;
    resp_dfof_pref = indOnly(resp_dfof_stim,pref_sf);
    ind_use = intersect(resp_ind,find(resp_dfof_pref>0.05));
    norm_dfof_stim_pref = indOnly(norm_dfof_stim,pref_sf);

    figure;
    subplot(1,2,1)
    scatter(repmat(sfs,[nCells 1]), norm_dfof_stim_nan,'ok')
    hold on
    errorbar(sfs,mean(norm_dfof_stim_nan,1,'omitnan'),std(norm_dfof_stim_nan,[],1,'omitnan')./sqrt(sum(~isnan(norm_dfof_stim_nan),1)),'or')
    xlabel('SF')
    ylabel('Norm resp')
    xlim([0 0.7])

    subplot(1,2,2)
    scatter(sfs(pref_sf(ind_use)),norm_dfof_stim_pref(ind_use),'ok')
    hold on
    for i = 1:nSF
        ind = intersect(ind_use,find(pref_sf==i));
        errorbar(sfs(i),mean(norm_dfof_stim_pref(ind),1,'omitnan'),std(norm_dfof_stim_pref(ind),[],1,'omitnan'),'or')
    end
    xlabel('Pref SF')
    ylabel('Norm resp')
    xlim([0 0.7])
    sgtitle([mouse ' ' date '- ' num2str(length(ind_use)) ' resp cells'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_adaptBySF.pdf']),'-dpdf','-fillpage')

    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']), 'norm_dfof_stim_nan', 'norm_dfof_stim_pref', 'resp_dfof_stim', 'pref_sf', 'norm_dfof_stim', 'tc_one_dfof', 'tc_two_dfof', 'nCells','frame_rate', 'h1_ori')
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']),'cStimOne','cStimTwo','ori_mat','nOri','oris','sf_mat', 'sfs', 'nSF', 'tISITime','ISIs','nISI','stimOneTime', 'stimOneTimes', 'nTime', 'stimOneContrast', 'stimOneCons','base_win','resp_win','tr_ind')
end

if input.doRandSF & input.doRandDir
    resp_dfof_stim = zeros(nCells,nOri,nSF,nISI,2);
    tr_ind = cell(nOri,nSF,nISI);
    tr_n = zeros(nOri,nSF,nISI);
    for i = 1:nOri
        ind = find(ori_mat == oris(i));
        for ii = 1:nSF
            ind2 = find(sf_mat == sfs(ii));
            for iii = 1:nISI
                ind_use = intersect(ind, intersect(ind2, find(tISITime == ISIs(iii))));
                tr_ind{i,ii,iii} = ind_use;
                tr_n(i,ii,iii) = length(ind_use);
                resp_temp1 = nanmean(tc_one_dfof(:,:,ind_use),3);
                resp_temp2 = nanmean(tc_two_dfof(:,:,ind_use),3);
                resp_dfof_stim(:,i,ii,iii,1) = squeeze(mean(resp_temp1(resp_win,:),1)-mean(resp_temp1(base_win,:),1));
                resp_dfof_stim(:,i,ii,iii,2) = squeeze(mean(resp_temp2(resp_win,:),1)-mean(resp_temp2(base_win,:),1));
            end
        end
    end
    
    figure; 
    resp_ind = find(sum(h1_ori,[2 3]));
    [n n2] = subplotn(length(resp_ind));
    for iC = 1:length(resp_ind)
        iCell = resp_ind(iC);
        subplot(n,n2,iC)
        imagesc(squeeze(resp_dfof_stim(iCell,:,:,1,1)))
        %title(num2str(sum(sum(h1_ori(iCell,:,:),2),3)>0))
        axis off
        clim([0 0.4])
    end

    [pref_val pref_sf] = max(mean(resp_dfof_stim(:,:,:,1,1),2),[],3);
    [pref_val pref_dir] = max(resp_dfof_stim(:,:,pref_sf,1,1),[],2);
    norm_dfof_stim = resp_dfof_stim(:,:,:,1,2)./resp_dfof_stim(:,:,:,1,1);
    norm_sf_avg = zeros(nSF,2);
    norm_sf_all = [];
    norm_prefsf_avg = zeros(nSF,2);
    norm_prefsf_all = [];
    norm_dir_avg = zeros(nSF,2);
    norm_dir_all = [];
    norm_prefdir_avg = zeros(nSF,2);
    norm_prefdir_all = [];
    resp_ind_sf = cell(1,nSF);
    resp_sf_n = zeros(1,nSF);
    for iSF = 1:nSF
        ind_sf = find(sum(sum(h1_ori(:,:,iSF),2),3));
        ind_sf_pref = intersect(ind_sf,find(pref_sf == iSF));
        resp_ind_sf{iSF} = ind_sf;
        resp_sf_n(iSF) = length(ind_sf); 
        norm_sf = [];
        norm_prefsf = [];
        norm_dir = [];
        norm_prefdir = [];
        for iOri = 1:nOri
            ind_dir = find(sum(sum(h1_ori(:,iOri,:),2),3));
            ind_dir_pref = find(pref_dir == iOri);
            ind_use = intersect(find(resp_dfof_stim(:,iOri,iSF,1,1)>0.05),intersect(ind_dir_pref,ind_sf));
            norm_sf = [norm_sf; [norm_dfof_stim(ind_use,iOri,iSF) sfs(iSF).*ones(length(ind_use),1) oris(iOri).*ones(length(ind_use),1) resp_dfof_stim(ind_use,iOri,iSF,1,1) ind_use]];
            ind_use = intersect(find(resp_dfof_stim(:,iOri,iSF,1,1)>0.05),intersect(ind_dir_pref,ind_sf_pref));
            norm_prefsf = [norm_prefsf; [norm_dfof_stim(ind_use,iOri,iSF) sfs(iSF).*ones(length(ind_use),1) oris(iOri).*ones(length(ind_use),1) resp_dfof_stim(ind_use,iOri,iSF,1,1) ind_use]];
            ind_use = intersect(find(resp_dfof_stim(:,iOri,iSF,1,1)>0.05),intersect(ind_dir,ind_sf_pref));
            norm_dir = [norm_dir; [norm_dfof_stim(ind_use,iOri,iSF) sfs(iSF).*ones(length(ind_use),1) oris(iOri).*ones(length(ind_use),1) resp_dfof_stim(ind_use,iOri,iSF,1,1) ind_use]];
            ind_use = intersect(find(resp_dfof_stim(:,iOri,iSF,1,1)>0.05),intersect(ind_dir_pref,ind_sf_pref));
            norm_prefdir = [norm_prefdir; [norm_dfof_stim(ind_use,iOri,iSF) sfs(iSF).*ones(length(ind_use),1) oris(iOri).*ones(length(ind_use),1) resp_dfof_stim(ind_use,iOri,iSF,1,1) ind_use]];
        end
        norm_sf_all = [norm_sf_all; norm_sf];
        norm_temp = norm_sf(:,1);
        norm_sf_avg(iSF,1) = mean(norm_temp);
        norm_sf_avg(iSF,2) = std(norm_temp)./sqrt(size(norm_temp,1));
    
        norm_prefsf_all = [norm_prefsf_all; norm_prefsf];
        norm_temp = norm_prefsf(:,1);
        norm_prefsf_avg(iSF,1) = mean(norm_temp);
        norm_prefsf_avg(iSF,2) = std(norm_temp)./sqrt(size(norm_temp,1));

        norm_dir_all = [norm_dir_all; norm_dir];
        norm_temp = norm_dir(:,1);
        norm_dir_avg(iSF,1) = mean(norm_temp);
        norm_dir_avg(iSF,2) = std(norm_temp)./sqrt(size(norm_temp,1));
    
        norm_prefdir_all = [norm_prefdir_all; norm_prefdir];
        norm_temp = norm_prefdir(:,1);
        norm_prefdir_avg(iSF,1) = mean(norm_temp);
        norm_prefdir_avg(iSF,2) = std(norm_temp)./sqrt(size(norm_temp,1));
    end
    
    figure;
    subplot(2,1,1)
    swarmchart(norm_sf_all(:,2), norm_sf_all(:,1),'k')
    hold on
    errorbar(sfs,norm_sf_avg(:,1),norm_sf_avg(:,2),'-or')
    title('Responsive')
    subplot(2,1,2)
    swarmchart(norm_prefsf_all(:,2), norm_prefsf_all(:,1),'k')
    hold on
    errorbar(sfs,norm_prefsf_avg(:,1),norm_prefsf_avg(:,2),'-or')
    title('Preferred')
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_AdaptBySF.pdf']),'-dpdf','-fillpage')
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']), 'resp_dfof_stim', 'pref_dir', 'pref_sf', 'norm_dfof_stim', 'norm_prefsf_all', 'norm_sf_all','norm_prefdir_all', 'norm_dir_all', 'tc_one_dfof', 'tc_two_dfof', 'nCells','frame_rate', 'h1_ori')
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']),'cStimOne','cStimTwo','ori_mat','nOri','oris','sf_mat', 'sfs', 'nSF', 'tISITime','ISIs','nISI','stimOneTime', 'stimOneTimes', 'nTime', 'stimOneContrast', 'stimOneCons','base_win','resp_win','tr_ind')
    
    figure; 
    scatter(norm_sf_all(:,1),norm_sf_all(:,3))
end
%%
max_dir = zeros(1,nCells);
neighbor_dir = zeros(2,nCells);
max_val = zeros(1,nCells);
snr = zeros(1,nCells);
snr_neighbor = zeros(1,nCells);
snr_resp = zeros(1,nCells);
h = zeros(1,nCells);
norm_df = zeros(nCells,nISI);
tc_dfof = zeros(size(tc_one_dfof,1),nCells,nISI);
tc_dfof_pref = zeros(size(tc_one_dfof,1),nCells,nISI,2);
resp_dfof_pref = zeros(nCells,nISI,2);
tc_dfof_neighbor = zeros(size(tc_one_dfof,1),nCells,nISI,2);
resp_dfof_neighbor = zeros(nCells,nISI,2);

for iC = 1:nCells
    [max_val(iC) max_dir(iC)] = max(mean(resp_dfof_dir(iC,:,:,1),3),[],2);
    neighbor_dir(:,iC) = [max_dir(iC)+1 max_dir(iC)-1];
    if find(neighbor_dir(:,iC) == 0)
        neighbor_dir(find(neighbor_dir(:,iC) == 0),iC) = nOri;
    end
    if find(neighbor_dir(:,iC) > nOri)
        neighbor_dir(find(neighbor_dir(:,iC) > nOri),iC) = 1;
    end
    ind = find(ori_mat == oris(max_dir(iC)));
    h(iC) = ttest(squeeze(mean(tc_one_dfof(resp_win,iC,ind),1)),squeeze(mean(tc_one_dfof(base_win,iC,ind),1)),'tail','right');
    norm_df(iC,:,1) = squeeze(resp_dfof_dir(iC,max_dir(iC),:,2)./resp_dfof_dir(iC,max_dir(iC),:,1));
    norm_df(iC,:,2) = squeeze(resp_dfof_dir(iC,max_dir(iC),:,2)./mean(resp_dfof_dir(iC,max_dir(iC),:,1),3));
    for i = 1:nISI
        tc_dfof(:,iC,i,1) = squeeze((nanmean(tc_two_dfof(:,iC,tr_ind{max_dir(iC),i}),3)-nanmean(mean(tc_two_dfof(base_win,iC,tr_ind{max_dir(iC),i}),1),3))./resp_dfof_dir(iC,max_dir(iC),i,1));
        tc_dfof(:,iC,i,2) = squeeze((nanmean(tc_two_dfof(:,iC,tr_ind{max_dir(iC),i}),3)-nanmean(mean(tc_two_dfof(base_win,iC,tr_ind{max_dir(iC),i}),1),3))./mean(resp_dfof_dir(iC,max_dir(iC),:,1),3));
        tc_dfof_pref(:,iC,i,1) = nanmean(tc_one_dfof(:,iC,tr_ind{max_dir(iC),i}),3);
        tc_dfof_pref(:,iC,i,2) = nanmean(tc_two_dfof(:,iC,tr_ind{max_dir(iC),i}),3)-nanmean(mean(tc_two_dfof(base_win,iC,tr_ind{max_dir(iC),i}),1),3);
        resp_dfof_pref(iC,i,1) = squeeze(mean(tc_dfof_pref(resp_win,iC,i,1),1)-mean(tc_dfof_pref(base_win,iC,i,1),1));
        resp_dfof_pref(iC,i,2) = squeeze(mean(tc_dfof_pref(resp_win,iC,i,2),1)-mean(tc_dfof_pref(base_win,iC,i,2),1));
        tc_dfof_neighbor(:,iC,i,1) = nanmean(tc_one_dfof(:,iC,[tr_ind{neighbor_dir(1,iC),i} tr_ind{neighbor_dir(2,iC),i}]),3);
        tc_dfof_neighbor(:,iC,i,2) = nanmean(tc_two_dfof(:,iC,[tr_ind{neighbor_dir(1,iC),i} tr_ind{neighbor_dir(2,iC),i}]),3)-nanmean(mean(tc_two_dfof(base_win,iC,[tr_ind{neighbor_dir(1,iC),i} tr_ind{neighbor_dir(2,iC),i}]),1),3);
        resp_dfof_neighbor(iC,i,1) = squeeze(mean(tc_dfof_neighbor(resp_win,iC,i,1),1)-mean(tc_dfof_neighbor(base_win,iC,i,1),1));
        resp_dfof_neighbor(iC,i,2) = squeeze(mean(tc_dfof_neighbor(resp_win,iC,i,2),1)-mean(tc_dfof_neighbor(base_win,iC,i,2),1));
    end
    snr(iC) = max_val(iC)./std(tc_dfof_pref(1:20,iC,1,1),[],1);
    snr_neighbor(iC) = resp_dfof_neighbor(iC,1,1)./std(tc_dfof_neighbor(1:20,iC,1,1),[],1);
    snr_resp(iC) = max_val(iC)./std(mean(tc_one_dfof(resp_win,iC,tr_ind{max_dir(iC),i}),1),[],3);
end
resp_ind= find(h & max_val>0.03);

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

figure;
subplot(2,2,1)
norm_resp_pref = resp_dfof_pref(:,:,2)./resp_dfof_pref(:,:,1);
% errorbar(mean(resp_dfof_pref(resp_ind,:,1),1), mean(norm_resp_pref,1), std(norm_resp_pref,[],1)./sqrt(length(resp_ind)), std(norm_resp_pref,[],1)./sqrt(length(resp_ind)), std(resp_dfof_pref(resp_ind,:,1),[],1)./sqrt(length(resp_ind)),std(resp_dfof_pref(resp_ind,:,1),[],1)./sqrt(length(resp_ind)),'o')
% hold on
% resp_ind_neighbor = intersect(resp_ind, find(resp_dfof_neighbor(:,:,1)>0.03));
norm_resp_neighbor = resp_dfof_neighbor(resp_ind,:,2)./resp_dfof_neighbor(resp_ind,:,1);
norm_resp_pref_n = resp_dfof_pref(resp_ind,:,2)./resp_dfof_pref(resp_ind,:,1);
% errorbar(mean(resp_dfof_neighbor(resp_ind_neighbor,:,1),1), mean(norm_resp_neighbor,1), std(norm_resp_neighbor,[],1)./sqrt(length(resp_ind_neighbor)), std(norm_resp_neighbor,[],1)./sqrt(length(resp_ind_neighbor)), std(resp_dfof_neighbor(resp_ind_neighbor,:,1),[],1)./sqrt(length(resp_ind_neighbor)),std(resp_dfof_neighbor(resp_ind_neighbor,:,1),[],1)./sqrt(length(resp_ind_neighbor)),'o')
% ylim([0 1])
% xlim([0 .2])
% ylabel('Normalized resp')
% xlabel('Response amplitude (dF/F)')

scatter(resp_dfof_pref(resp_ind,:,1),norm_resp_pref(resp_ind,:))
ylabel('Normalized resp')
xlabel('Response amplitude (dF/F)')
ylim([-.5 2])
xlim([0 .5])
subplot(2,2,2)
scatter(snr(resp_ind),norm_resp_pref(resp_ind,:))
ylabel('Normalized resp')
xlabel('SNR')
ylim([-.5 2])
xlim([0 70])
subplot(2,2,3)
scatter(snr_resp(resp_ind),norm_resp_pref(resp_ind,:))
ylabel('Normalized resp')
xlabel('SNR- resp')
ylim([-.5 2])
xlim([0 1.5])
subplot(2,2,3)

scatter(snr_neighbor(resp_ind),norm_resp_neighbor)
ylabel('Normalized resp- neighbor')
xlabel('SNR neighbor')
ylim([-.5 2])
xlim([0 70])
% scatter(norm_resp_pref_n, norm_resp_neighbor)
% xlabel('Normalized resp- pref')
% ylabel('Normalized resp- neighbor')
% xlim([-.5 2])
% ylim([-.5 2])
% refline(1)
title(num2str(chop(triu2vec(corrcoef(norm_resp_pref_n, norm_resp_neighbor)),2)))
select = resp_dfof_pref(resp_ind,:,1)./sum(resp_dfof_dir(resp_ind,:,:,1),2);
subplot(2,2,4)
scatter(norm_resp_pref,select);
xlim([-.5 2])
ylim([0 1])
title(num2str(chop(triu2vec(corrcoef(norm_resp_pref, select)),2)))
xlabel('Normalized resp- pref')
ylabel('Selectivity')

figure;
ind = intersect(resp_ind, find(resp_dfof_pref(:,:,1)<0.05));
for iC = 1:length(ind)
    subplot(6,6,iC)
    iCell = ind(iC);
    plot(tc_dfof_pref(:,iCell,i,1))
    hold on
    plot(tc_dfof_pref(:,iCell,i,2))
    title(num2str(chop(resp_dfof_pref(iCell,:,2)./resp_dfof_pref(iCell,:,1),2)))
    vline([base_win resp_win])
end

figure;
start = 1;
for iC = 7:12
    subplot(6,2,start)
    iCell = ind(iC);
    plot(tc_dfof_pref(:,iCell,i,1))
    hold on
    plot(tc_dfof_pref(:,iCell,i,2))
    title(num2str(chop(resp_dfof_pref(iCell,:,2)./resp_dfof_pref(iCell,:,1),2)))
    vline([base_win resp_win])
    subplot(6,2,start+1)
    iCell = ind(iC);
    plot(tc_dfof_neighbor(:,iCell,i,1))
    hold on
    plot(tc_dfof_neighbor(:,iCell,i,2))
    %title(num2str(chop(resp_dfof_neighbor(iCell,:,2)./resp_dfof_neighbor(iCell,:,1),2)))
    title(num2str(chop(select(iCell),2)))
    vline([base_win resp_win])
    start = start+2;
end

figure;
for iCell = 1:36
    iC= resp_ind(iCell);
    subplot(6,6,iCell)
    plot(1:100,squeeze(tc_one_dfof(:,iC,tr_ind{max_dir(iC),i})),'c')
    hold on
    plot(1:100,mean(tc_one_dfof(:,iC,tr_ind{max_dir(iC),i}),3,'omitnan'),'k')
end
figure;
for iCell = 1:36
    iC= resp_ind(iCell);
    subplot(6,6,iCell)
    n = max_dir(iC)+1;
    if n>nOri
        n = 1;
    end
    plot(1:100,squeeze(tc_one_dfof(:,iC,tr_ind{n,i})),'c')
    hold on
    plot(1:100,mean(tc_one_dfof(:,iC,tr_ind{n,i}),3,'omitnan'),'k')
end
    
    
%%

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']), 'tc_dfof','resp_dfof_dir', 'max_dir', 'norm_df', 'tc_dfof_pref',  'tc_one_dfof', 'tc_two_dfof', 'nCells','frame_rate', 'resp_ind', 'resp_dfof_pref', 'resp_dfof_neighbor','norm_resp_pref','snr','snr_resp')
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']),'cStimOne','cStimTwo','ori_mat','nOri','oris','tISITime','ISIs','nISI','stimOneTime', 'stimOneTimes', 'nTime', 'stimOneContrast', 'stimOneCons','base_win','resp_win','tr_ind')

    
%%
tc_one_boot = zeros(100,nCells,1000);
tc_two_boot = zeros(100,nCells,1000);
for iOri = 1:nOri
    cell_ind = find(max_dir==iOri);
    for i = 1:1000
        tr_sub = randperm(length(tr_ind{iOri,1}),floor(length(tr_ind{iOri,1})./2));
        tc_one_boot(:,cell_ind,i) = mean(tc_one_dfof(:,cell_ind,tr_ind{iOri,1}(tr_sub)),3,'omitnan');
        tc_two_boot(:,cell_ind,i) = mean(tc_two_dfof(:,cell_ind,tr_ind{iOri,1}(tr_sub)),3,'omitnan');
    end
end

resp_one_boot = mean(tc_one_boot(resp_win,:,:),1)-mean(tc_one_boot(base_win,:,:),1);
resp_two_boot = mean(tc_two_boot(resp_win,:,:),1)-mean(tc_two_boot(base_win,:,:),1);

norm_resp_boot = squeeze(resp_two_boot./resp_one_boot);

norm_resp_sort = sort(norm_resp_boot,2);

[norm_resp_med norm_resp_med_ind] = sort(norm_resp_sort(resp_ind,500),1);
norm_resp_90 = norm_resp_sort(resp_ind(norm_resp_med_ind),[100 900]);
figure; scatter(1:length(resp_ind),norm_resp_med-1);
hold on
scatter(1:length(resp_ind),norm_resp_90(:,1)-1);
scatter(1:length(resp_ind),norm_resp_90(:,2)-1);
scatter(1:length(resp_ind),norm_resp_pref(resp_ind(norm_resp_med_ind))-1);
hold on;
med_ad = median(abs(norm_resp_med-1));
hline([-med_ad med_ad])
less_ad1 = find(norm_resp_med-1>-med_ad,1,'first');
less_ad2 = find(norm_resp_med-1<med_ad,1,'last');
less_ad_all = [less_ad1:less_ad2];
more_ad_all = setdiff(1:length(resp_ind), less_ad_all);

ylabel('Adapt Ix')
xlabel('Cell sorted by Adapt Ix')
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_adaptCI.pdf']),'-dpdf','-bestfit')

figure; scatter(1:length(resp_ind),resp_dfof_pref(resp_ind(norm_resp_med_ind),:,1))
ylabel('R1 dF/F')
xlabel('Cell sorted by Adapt Ix')
vline([less_ad1 less_ad2])
figure; scatter(1:length(resp_ind),norm_resp_90(:,2)-norm_resp_90(:,1));
vline([less_ad1 less_ad2])
ylabel('10-90 Boot Diff of Adapt Ix')
xlabel('Cell sorted by Adapt Ix')
figure; scatter(resp_dfof_pref(resp_ind(norm_resp_med_ind),:,1),norm_resp_90(:,2)-norm_resp_90(:,1))
xlabel('R1 dF/F')
ylabel('10-90 Boot Diff of Adapt Ix')
figure; scatter(1:length(resp_ind),sum(norm_resp_sort(resp_ind(norm_resp_med_ind),:)<1,2));
xlabel('Cell sorted by Adapt Ix')
vline([less_ad1 less_ad2])
ylabel('Fraction of boots with adaptation <1')
figure; scatter(resp_dfof_pref(resp_ind(norm_resp_med_ind),:,1),sum(norm_resp_sort(resp_ind(norm_resp_med_ind),:)<1,2));

%[u1 s1 v1] = pca(squeeze(mean(tc_one_dfof(resp_win,resp_ind(norm_resp_med_ind(1:floor(length(resp_ind/2)))),:),1)));

%%

calib = 1/26.6; %mm per pixel

% Load and combine eye tracking data
nrun = length(ImgFolder);
data = [];
for irun =  1:nrun
    fn = [ImgFolder{irun} '_000_000_eye.mat'];

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
