clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriSingleStimRandDirAdapt_ExptList';

rc = behavConstsAV;
eval(ds)
nexp = size(expt,2);

frame_rate = 15;
for iexp  = 5
%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};
ImgFolder = expt(iexp).coFolder;
time = expt(iexp).coTime;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
%LG_base = '\\CRASH.dhe.duke.edu\data\home\lindsey';

fprintf([mouse ' ' date '\n'])

%% Test stim analysis
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))

%%
% if doRedChannel == 0
%     red_cells = [];
% end

nCells = size(npSub_tc,2);
prewin_frames = frame_rate;
nFramesOn = unique(celleqel2mat_padded(input.nStimTwoFramesOn));
postwin_frames = frame_rate*5;
tt = (1-prewin_frames:postwin_frames).*(1/frame_rate);
data_resp = nan(prewin_frames+postwin_frames,nCells,nTrials);
data_f = nan(1,nCells,nTrials);

for itrial = 1:nTrials
    if cStimOn(itrial) + postwin_frames < sz(3)
        data_resp(:,:,itrial) = npSub_tc(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
    end
end

data_f = mean(data_resp(ceil(prewin_frames/2):prewin_frames,:,:),1);
data_dfof_tc = (data_resp-data_f)./data_f;

ind_stimAlone = intersect(find(stimCon_all),find(maskCon_all==0));
ind_maskAlone = intersect(find(stimCon_all==0),find(maskCon_all));
ind_plaid = intersect(find(stimCon_all),find(maskCon_all));
ind_blank = intersect(find(stimCon_all==0),find(maskCon_all==0));

ind_noadapt = find(adaptStim == 0);
ind_singadapt = find(adaptStim);

resp_cell_noadapt = cell(nStimDir,2);
resp_cell_singadapt = cell(nStimDir,2);
base_cell_noadapt = cell(nStimDir,2);
base_cell_singadapt = cell(nStimDir,2);
trialInd_singadapt = cell(nStimDir,2);
trialInd_noadapt = cell(nStimDir,2);
trialsperstim_noadapt = zeros(nStimDir,2);
trialsperstim_singadapt = zeros(nStimDir,2);
h_resp =zeros(nCells,nStimDir,2);
p_resp =zeros(nCells,nStimDir,2);

base_win = prewin_frames-ceil(frame_rate/2):prewin_frames;
resp_win = prewin_frames+5:prewin_frames+nFramesOn;
data_dfof_dir_tc_avg_singadapt = nan(prewin_frames+postwin_frames, nCells, nStimDir, 2);
data_dfof_dir_tc_avg_noadapt = nan(prewin_frames+postwin_frames, nCells, nStimDir, 2);

avg_resp_dir = zeros(nCells, nStimDir, 2, 2, 2); % 1- cells; 2- directions; 3- grating/plaid; 4- no adapt/sing adapt; 5- mean/sem

all_resp_dir = [];
all_resp_plaid = [];
all_dir = [];
all_plaid = [];
nStim = nStimDir;
for iDir = 1:nStimDir
    ind_stimdir = find(stimDir_all == stimDirs(iDir));
    ind_maskdir = find(maskDir_all == maskDirs(iDir));
    ind_diralone = [intersect(ind_stimdir, ind_stimAlone), intersect(ind_maskdir, ind_maskAlone)];
    ind_dirplaid = [intersect(ind_stimdir, ind_plaid)];
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
    trialsperstim_singadapt(iDir,1) = length(intersect(ind_singadapt,ind_diralone));
    trialsperstim_singadapt(iDir,2) = length(intersect(ind_singadapt,ind_dirplaid));
    trialInd_singadapt{iDir,1} = intersect(ind_singadapt,ind_diralone);
    trialInd_singadapt{iDir,2} = intersect(ind_singadapt,ind_dirplaid);
    resp_cell_singadapt{iDir,1} = squeeze(mean(data_dfof_tc(resp_win,:,intersect(ind_singadapt,ind_diralone)),1));
    resp_cell_singadapt{iDir,2} = squeeze(mean(data_dfof_tc(resp_win,:,intersect(ind_singadapt,ind_dirplaid)),1));
    avg_resp_dir(:,iDir,1,2,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,intersect(ind_singadapt,ind_diralone)),1),3,'omitnan'));
    avg_resp_dir(:,iDir,1,2,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,intersect(ind_singadapt,ind_diralone)),1),[],3,'omitnan')./sqrt(length(intersect(ind_singadapt,ind_diralone))));
    avg_resp_dir(:,iDir,2,2,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,intersect(ind_singadapt,ind_dirplaid)),1),3,'omitnan'));
    avg_resp_dir(:,iDir,2,2,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,intersect(ind_singadapt,ind_dirplaid)),1),[],3,'omitnan')./sqrt(length(intersect(ind_singadapt,ind_dirplaid))));
    data_dfof_dir_tc_avg_singadapt(:,:,iDir,1) = squeeze(mean(data_dfof_tc(:,:,intersect(ind_singadapt,ind_diralone)),3,'omitnan'));
    data_dfof_dir_tc_avg_singadapt(:,:,iDir,2) = squeeze(mean(data_dfof_tc(:,:,intersect(ind_singadapt,ind_dirplaid)),3,'omitnan'));
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

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']), 'data_dfof_dir_tc_avg_singadapt', 'resp_cell_singadapt','base_cell_singadapt','data_dfof_dir_tc_avg_noadapt', 'resp_cell_noadapt','base_cell_noadapt', 'data_dfof_tc', 'avg_resp_dir', 'h_resp','resp_ind', 'tt', 'frame_rate');
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']), 'prewin_frames', 'postwin_frames', 'resp_win', 'base_win', 'trialsperstim_singadapt','trialsperstim_noadapt','trialInd_singadapt','trialInd_noadapt');

end
%%
maskDiff_all = celleqel2mat_padded(input.tMaskTwoGratingDirectionDeg) - celleqel2mat_padded(input.tStimTwoGratingDirectionDeg);
maskDiffs = unique(maskDiff_all);

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

Zc_use = intersect(resp_ind, find(Zc_all(1,:)>1.28 & Zc_all(1,:)-Zp_all(1,:)>1.28));
Zp_use = intersect(resp_ind, find(Zp_all(1,:)>1.28 & Zp_all(1,:)-Zc_all(1,:)>1.28));

[maxVal maxDir] = max(avg_resp_dir(:,:,1,1,1),[],2);
tests = [12 1 2];
orth = [3 4 5 9 10 11];
opp = [6 7 8];
tests_ind = [];
for i = 1:length(tests)
    tests_ind = [tests_ind; find(maxDir == tests(i))];
end
opp_ind = [];
for i = 1:length(opp)
    opp_ind = [opp_ind; find(maxDir == opp(i))];
end
orth_ind = [];
for i = 1:length(orth)
    orth_ind = [orth_ind; find(maxDir == orth(i))];
end

figure;
subplot(2,3,1)
Zc_use_temp = intersect(Zc_use,tests_ind);
plot(Zc(:,Zc_use_temp),'k')
set(gca,'XTick',[1 2],'XtickLabel',{'control','adapt'})
ylabel('Zc')
title('Test')
subplot(2,3,2)
Zc_use_temp = intersect(Zc_use,opp_ind);
plot(Zc(:,Zc_use_temp),'k')
set(gca,'XTick',[1 2],'XtickLabel',{'control','adapt'})
ylabel('Zc')
title('Opp')
subplot(2,3,3)
Zc_use_temp = intersect(Zc_use,orth_ind);
plot(Zc(:,Zc_use_temp),'k')
set(gca,'XTick',[1 2],'XtickLabel',{'control','adapt'})
ylabel('Zc')
title('Orth')

subplot(2,3,4)
Zp_use_temp = intersect(Zp_use,tests_ind);
plot(Zp(:,Zp_use_temp),'k')
set(gca,'XTick',[1 2],'XtickLabel',{'control','adapt'})
ylabel('Zp')
title('Test')
subplot(2,3,5)
Zp_use_temp = intersect(Zp_use,opp_ind);
plot(Zp(:,Zp_use_temp),'k')
set(gca,'XTick',[1 2],'XtickLabel',{'control','adapt'})
ylabel('Zp')
title('Opp')
subplot(2,3,6)
Zp_use_temp = intersect(Zp_use,orth_ind);
plot(Zp(:,Zp_use_temp),'k')
set(gca,'XTick',[1 2],'XtickLabel',{'control','adapt'})
ylabel('Zc')
title('Orth')

[h p] = ttest(Zc(1,Zc_use),Zc(2,Zc_use));
title(num2str(chop(p,3)))
subplot(2,3,2)
for i = 1:length(Zp_use)
    plot(Zp(:,Zp_use(i)),'k')
    hold on
end
set(gca,'XTick',[1 2],'XtickLabel',{'control','adapt'})
ylabel('Zp')
[h p] = ttest(Zp(1,Zp_use),Zp(2,Zp_use));
title(num2str(chop(p,3)))

%% direction and pattern selectivity
[maxVal maxDir] = max(avg_resp_dir(:,:,1,1,1),[],2);
oppDir = maxDir+nStimDir/2;
oppDir(find(oppDir>nStimDir)) = oppDir(find(oppDir>nStimDir))-nStimDir;
oppVal = indOnly(avg_resp_dir(:,:,1,1,1),oppDir);
oppVal(find(oppVal<0)) = 0;
DSI = (maxVal-oppVal)./(maxVal+oppVal);
gDSI = circ_var(avg_resp_dir(:,:,1,1,1),[],[],2);
maxSig = indOnly(h_resp(:,:,1),maxDir);

ind_DSI = find(DSI>0.5);
ind_use = intersect(find(p_anova_dir<0.05),resp_ind_dir);
ind_pref = intersect(resp_ind_dir, find(maxDir==1));
ind_orth = intersect(resp_ind_dir, find(maxDir==4));
ind_anova = find(p_anova_dir);

avg_resp_dir_shift = avg_resp_dir;
avg_resp_dir_shift(:,:,2,:,:) = circshift(avg_resp_dir_shift(:,:,2,:,:),3,2);
resp = intersect(resp_ind_dir, find(maxVal>0.1));
figure;
start = 1;
for i = 1:length(resp)
    if start > 36
        figure;
        start = 1;
    end
    subplot(6,6,start)
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir(resp(i),:,1,1,1),1) mean(avg_resp_dir(resp(i),1,1,1,1),1)])
    hold on
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(resp(i),:,2,1,1),1) mean(avg_resp_dir_shift(resp(i),1,2,1,1),1)])
    start = start+1;
end

figure;
start = 1;
for i = 1:36% length(ind_use)
    if start > 36
        figure;
        start = 1;
    end
    subplot(6,6,start)
    errorbar(stimDirs,avg_resp_dir_shift(ind_use(i),:,1,1,1),avg_resp_dir_shift(ind_use(i),:,1,1,2))
    hold on
    errorbar(stimDirs,avg_resp_dir_shift(ind_use(i),:,1,2,1),avg_resp_dir_shift(ind_use(i),:,1,2,2))
    title(num2str(ind_use(i)))
    start = start+1;
end

figure;
start = 1;
for i = 1:36 % length(ind_use)
    if start > 36
        figure;
        start = 1;
    end
    subplot(6,6,start)
    errorbar(stimDirs,avg_resp_dir_shift(ind_use(i),:,2,1,1),avg_resp_dir_shift(ind_use(i),:,2,1,2))
    hold on
    errorbar(stimDirs,avg_resp_dir_shift(ind_use(i),:,2,2,1),avg_resp_dir_shift(ind_use(i),:,2,2,2))
    title(num2str(ind_use(i)))
    start = start+1;
end


figure;
start = 1;
for i = 1:36 %length(ind_use)
    if start > 36
        figure;
        start = 1;
    end
    subplot(6,6,start)
    errorbar(stimDirs,avg_resp_dir_shift(i,:,2,1,1),avg_resp_dir_shift(i,:,2,1,2))
    hold on
    errorbar(stimDirs,avg_resp_dir_shift(i,:,2,2,1),avg_resp_dir_shift(i,:,2,2,2))
    title(num2str(i))
    start = start+1;
    vline([60 300])
end

figure;
subplot(2,1,1)
errorbar(stimDirs,mean(avg_resp_dir_shift(ind_pref,:,1,1,1),1),std(avg_resp_dir_shift(ind_pref,:,1,1,1),[],1)./sqrt(length(ind_pref)))
hold on
errorbar(stimDirs,mean(avg_resp_dir_shift(ind_pref,:,1,2,1),1),std(avg_resp_dir_shift(ind_pref,:,1,2,1),[],1)./sqrt(length(ind_pref)))
vline([0 180])
subplot(2,1,2)
errorbar(stimDirs,mean(avg_resp_dir_shift(ind_pref,:,2,1,1),1),std(avg_resp_dir_shift(ind_pref,:,2,1,1),[],1)./sqrt(length(ind_pref)))
hold on
errorbar(stimDirs,mean(avg_resp_dir_shift(ind_pref,:,2,2,1),1),std(avg_resp_dir_shift(ind_pref,:,2,2,1),[],1)./sqrt(length(ind_pref)))
vline([60 300])

figure;
subplot(2,1,1)
errorbar(stimDirs,mean(avg_resp_dir_shift(ind_orth,:,1,1,1),1),std(avg_resp_dir_shift(ind_orth,:,1,1,1),[],1)./sqrt(length(ind_orth)))
hold on
errorbar(stimDirs,mean(avg_resp_dir_shift(ind_orth,:,1,2,1),1),std(avg_resp_dir_shift(ind_orth,:,1,2,1),[],1)./sqrt(length(ind_orth)))
vline([0 180])
subplot(2,1,2)
errorbar(stimDirs,mean(avg_resp_dir_shift(ind_orth,:,2,1,1),1),std(avg_resp_dir_shift(ind_orth,:,2,1,1),[],1)./sqrt(length(ind_orth)))
hold on
errorbar(stimDirs,mean(avg_resp_dir_shift(ind_orth,:,2,2,1),1),std(avg_resp_dir_shift(ind_orth,:,2,2,1),[],1)./sqrt(length(ind_orth)))
vline([60 300])

ii=6;
figure;
for i = 1:12
    subplot(4,4,i)
    plot(data_dfof_dir_tc_avg_noadapt(:,ii,i,1))
    hold on
    plot(data_dfof_dir_tc_avg_singadapt(:,ii,i,1))
    ylim([-0.2 0.4])
end
subplot(4,4,i+1)
polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ii,:,1,1,1),1) mean(avg_resp_dir_shift(ii,1,1,1,1),1)])
hold on
polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ii,:,1,2,1),1) mean(avg_resp_dir_shift(ii,1,1,2,1),1)])

figure;
for i = 1:12
    subplot(4,4,i)
    plot(data_dfof_dir_tc_avg_noadapt(:,ii,i,2))
    hold on
    plot(data_dfof_dir_tc_avg_singadapt(:,ii,i,2))
    ylim([-0.2 0.4])
end
subplot(4,4,i+1)
polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ii,:,2,1,1),1) mean(avg_resp_dir_shift(ii,1,2,1,1),1)])
hold on
polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ii,:,2,2,1),1) mean(avg_resp_dir_shift(ii,1,2,2,1),1)])


for iCell = 1:10; %length(resp_ind_dir)
    iC = ind_use(iCell);
    subplot(5,2,iCell)
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(iC,:,1,1,1),1) mean(avg_resp_dir_shift(iC,1,1,1,1),1)])
    hold on
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(iC,:,2,1,1),1) mean(avg_resp_dir_shift(iC,1,2,1,1),1)])
end

figure; 
for i = 1:12
    subplot(3,4,i)
    shadedErrorBar(tt,mean(data_dfof_dir_tc_avg_noadapt(:,ind_pref,i,1),2),std(data_dfof_dir_tc_avg_noadapt(:,ind_pref,i,1),[],2)./sqrt(length(ind_pref)))
    hold on
    shadedErrorBar(tt,mean(data_dfof_dir_tc_avg_singadapt(:,ind_pref,i,1),2),std(data_dfof_dir_tc_avg_singadapt(:,ind_pref,i,1),[],2)./sqrt(length(ind_pref)),'lineProps','r')
    ylim([-0.1 1.5])
end
suptitle([mouse ' ' date '- Grating: control (black), adapt (red)- Pref 0'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_gratingTCs_prefAdapt.pdf']),'-dpdf','-bestfit')

figure; 
for i = 1:12
    subplot(3,4,i)
    shadedErrorBar(tt,mean(data_dfof_dir_tc_avg_noadapt(:,ind_pref,i,2),2),std(data_dfof_dir_tc_avg_noadapt(:,ind_pref,i,2),[],2)./sqrt(length(ind_pref)))
    hold on
    shadedErrorBar(tt,mean(data_dfof_dir_tc_avg_singadapt(:,ind_pref,i,2),2),std(data_dfof_dir_tc_avg_singadapt(:,ind_pref,i,2),[],2)./sqrt(length(ind_pref)),'lineProps','r')
    ylim([-0.1 1.5])
end
suptitle([mouse ' ' date '- Plaid: control (black), adapt (red) - Pref 0'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_plaidTCs_prefAdapt.pdf']),'-dpdf','-bestfit')

figure; 
for i = 1:12
    subplot(3,4,i)
    shadedErrorBar(tt,mean(data_dfof_dir_tc_avg_noadapt(:,ind_orth,i,1),2),std(data_dfof_dir_tc_avg_noadapt(:,ind_orth,i,1),[],2)./sqrt(length(ind_orth)))
    hold on
    shadedErrorBar(tt,mean(data_dfof_dir_tc_avg_singadapt(:,ind_orth,i,1),2),std(data_dfof_dir_tc_avg_singadapt(:,ind_orth,i,1),[],2)./sqrt(length(ind_orth)),'lineProps','r')
    ylim([-0.1 1.5])
end
suptitle([mouse ' ' date '- Grating: control (black), adapt (red)- Pref 90'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_gratingTCs_orthAdapt.pdf']),'-dpdf','-bestfit')

figure; 
for i = 1:12
    subplot(3,4,i)
    shadedErrorBar(tt,mean(data_dfof_dir_tc_avg_noadapt(:,ind_orth,i,2),2),std(data_dfof_dir_tc_avg_noadapt(:,ind_orth,i,2),[],2)./sqrt(length(ind_orth)))
    hold on
    shadedErrorBar(tt,mean(data_dfof_dir_tc_avg_singadapt(:,ind_orth,i,2),2),std(data_dfof_dir_tc_avg_singadapt(:,ind_orth,i,2),[],2)./sqrt(length(ind_orth)),'lineProps','r')
    ylim([-0.1 1.5])
end
suptitle([mouse ' ' date '- Plaid: control (black), adapt (red)- Pref 90'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_plaidTCs_orthAdapt.pdf']),'-dpdf','-bestfit')


figure; 
    subplot(1,2,1)
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir(ind_use,:,1,1,1),1) mean(avg_resp_dir(ind_use,1,1,1,1),1)])
    hold on
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir(ind_use,:,1,2,1),1) mean(avg_resp_dir(ind_use,1,1,2,1),1)])
    subplot(1,2,2)
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir(ind_use,:,2,1,1),1) mean(avg_resp_dir(ind_use,1,2,1,1),1)])
    hold on
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir(ind_use,:,2,2,1),1) mean(avg_resp_dir(ind_use,1,2,2,1),1)])
    

figure; 
    subplot(1,2,1)
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir(ind_pref,:,1,1,1),1) mean(avg_resp_dir(ind_pref,1,1,1,1),1)])
    hold on
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir(ind_pref,:,1,2,1),1) mean(avg_resp_dir(ind_pref,1,1,2,1),1)])
    subplot(1,2,2)
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir(ind_pref,:,2,1,1),1) mean(avg_resp_dir(ind_pref,1,2,1,1),1)])
    hold on
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir(ind_pref,:,2,2,1),1) mean(avg_resp_dir(ind_pref,1,2,2,1),1)])
sgtitle(num2str(length(ind_pref)))

figure; 
    subplot(1,2,1)
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir(ind_orth,:,1,1,1),1) mean(avg_resp_dir(ind_orth,1,1,1,1),1)])
    hold on
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir(ind_orth,:,1,2,1),1) mean(avg_resp_dir(ind_orth,1,1,2,1),1)])
    subplot(1,2,2)
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir(ind_orth,:,2,1,1),1) mean(avg_resp_dir(ind_orth,1,2,1,1),1)])
    hold on
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir(ind_orth,:,2,2,1),1) mean(avg_resp_dir(ind_orth,1,2,2,1),1)])
sgtitle(num2str(length(ind_orth)))

    ii = 6
figure; 
for i = 1:12
    subplot(3,4,i)
    plot(tt,mean(data_dfof_dir_tc_avg_noadapt(:,ii,i,1),2))
    hold on
    plot(tt,mean(data_dfof_dir_tc_avg_singadapt(:,ii,i,1),2))
    plot(tt,mean(data_dfof_dir_tc_avg_noadapt(:,ii,i,2),2))
    plot(tt,mean(data_dfof_dir_tc_avg_singadapt(:,ii,i,2),2))
    ylim([-0.1 1.5])
end   

figure; 
for i = 1:12
    subplot(3,4,i)
    plot(tt,mean(data_dfof_dir_tc_avg_noadapt(:,:,i,1),2))
    hold on
    plot(tt,mean(data_dfof_dir_tc_avg_singadapt(:,:,i,1),2))
    ylim([-0.1 1.5])
end 

figure; 
for i = 1:12
    subplot(3,4,i)
    plot(tt,mean(data_dfof_dir_tc_avg_noadapt(:,:,i,2),2))
    hold on
    plot(tt,mean(data_dfof_dir_tc_avg_singadapt(:,:,i,2),2))
    ylim([-0.1 1.5])
end  


figure;
subplot(2,2,1)
plot(mean(data_dfof_dir_tc_avg_noadapt(:,ind_pref,1,1),2))
hold on 
plot(mean(data_dfof_dir_tc_avg_singadapt(:,ind_pref,1,1),2))
ylim([-0.05 0.3])
subplot(2,2,2)
plot(mean(data_dfof_dir_tc_avg_noadapt(:,ind_pref,4,1),2))
hold on 
plot(mean(data_dfof_dir_tc_avg_singadapt(:,ind_pref,4,1),2))
ylim([-0.05 0.3])
subplot(2,2,3)
plot(mean(data_dfof_dir_tc_avg_noadapt(:,ind_pref,1,2),2))
hold on 
plot(mean(data_dfof_dir_tc_avg_singadapt(:,ind_pref,1,2),2))
ylim([-0.05 0.3])
subplot(2,2,4)
plot(mean(data_dfof_dir_tc_avg_noadapt(:,ind_pref,4,2),2))
hold on 
plot(mean(data_dfof_dir_tc_avg_singadapt(:,ind_pref,4,2),2))
ylim([-0.05 0.3])

for i = 1:5
ind_orth = intersect(resp_ind_dir, find(maxDir==i));
figure;
suptitle(['Ori = ' num2str(stimDirs(i))])
subplot(2,2,1)
plot(mean(data_dfof_dir_tc_avg_noadapt(:,ind_orth,1,1),2))
hold on 
plot(mean(data_dfof_dir_tc_avg_singadapt(:,ind_orth,1,1),2))
ylim([-0.05 0.5])
title('0 deg grating')
subplot(2,2,2)
plot(mean(data_dfof_dir_tc_avg_noadapt(:,ind_orth,4,1),2))
hold on 
plot(mean(data_dfof_dir_tc_avg_singadapt(:,ind_orth,4,1),2))
ylim([-0.05 0.5])
title('90 deg grating')
subplot(2,2,3)
plot(mean(data_dfof_dir_tc_avg_noadapt(:,ind_orth,1,2),2))
hold on 
plot(mean(data_dfof_dir_tc_avg_singadapt(:,ind_orth,1,2),2))
ylim([-0.05 0.5])
title('60 deg plaid- 0 deg grating')
subplot(2,2,4)
plot(mean(data_dfof_dir_tc_avg_noadapt(:,ind_orth,4,2),2))
hold on 
plot(mean(data_dfof_dir_tc_avg_singadapt(:,ind_orth,4,2),2))
ylim([-0.05 0.5])
title('150 deg plaid- 90 deg grating')
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

i = 1;
ind_orth = intersect(resp_ind_dir, find(maxDir==i));
figure;
[n n2] = subplotn(length(ind_orth));
for i = 1:length(ind_orth)
    subplot(n,n2,i)
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind_orth(i),:,1,1,1),1) mean(avg_resp_dir_shift(ind_orth(i),1,1,1,1),1)])
    hold on
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind_orth(i),:,1,2,1),1) mean(avg_resp_dir_shift(ind_orth(i),1,1,2,1),1)])
end
suptitle('Gratings')
figure;
for i = 1:length(ind_orth)
    subplot(n,n2,i)
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind_orth(i),:,2,1,1),1) mean(avg_resp_dir_shift(ind_orth(i),1,2,1,1),1)])
    hold on
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind_orth(i),:,2,2,1),1) mean(avg_resp_dir_shift(ind_orth(i),1,2,2,1),1)])
end
suptitle('Plaids')

i = 1;
ind_orth = intersect(resp_ind_dir, find(maxDir==i));
figure;
[n n2] = subplotn(length(ind_orth));
for i = 1:length(ind_orth)
    subplot(n,n2,i)
    errorbar(stimDirs, mean(avg_resp_dir_shift(ind_orth(i),:,1,1,1),1), mean(avg_resp_dir_shift(ind_orth(i),:,1,1,2),1))
    hold on
    errorbar(stimDirs, mean(avg_resp_dir_shift(ind_orth(i),:,2,1,1),1), mean(avg_resp_dir_shift(ind_orth(i),:,2,1,2),1))
end
suptitle('Gratings')
figure;
for i = 1:length(ind_orth)
    subplot(n,n2,i)
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind_orth(i),:,2,1,1),1) mean(avg_resp_dir_shift(ind_orth(i),1,2,1,1),1)])
    hold on
    polar(deg2rad([stimDirs stimDirs(1)]), [mean(avg_resp_dir_shift(ind_orth(i),:,2,2,1),1) mean(avg_resp_dir_shift(ind_orth(i),1,2,2,1),1)])
end
suptitle('Plaids')