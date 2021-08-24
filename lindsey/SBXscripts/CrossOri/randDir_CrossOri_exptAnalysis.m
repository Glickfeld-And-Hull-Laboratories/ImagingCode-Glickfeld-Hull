clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandDirRandPhase_ExptList';
eval(ds)
nexp = length(expt);

for iexp = 39

frame_rate = 15;

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc;
ImgFolder = expt(iexp).coFolder;
time = expt(iexp).coTime;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
%LG_base = '\\CRASH.dhe.duke.edu\data\home\lindsey';

fprintf([mouse ' ' date '\n'])

%% Pref direction analysis
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))

%%
if doRedChannel == 0
    red_cells = [];
end

prewin_frames = unique(celleqel2mat_padded(input.tItiWaitFrames))./3;
nFramesOn = unique(celleqel2mat_padded(input.nStimOneFramesOn));
postwin_frames = unique(celleqel2mat_padded(input.nStimOneFramesOn));
tt = (1-prewin_frames:postwin_frames).*(1/frame_rate);
data_resp = nan(prewin_frames+postwin_frames,nCells,nTrials-1);
data_f = nan(1,nCells,nTrials-1);

for itrial = 1:nTrials
    if cStimOn(itrial) + postwin_frames < sz(3)
        data_resp(:,:,itrial) = npSub_tc(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
    end
end

data_f = mean(data_resp(1:prewin_frames,:,:),1);
data_dfof_tc = (data_resp-data_f)./data_f;

ind_stimAlone = intersect(find(stimCon_all),find(maskCon_all==0));
ind_maskAlone = intersect(find(stimCon_all==0),find(maskCon_all));
ind_plaid = intersect(find(stimCon_all),find(maskCon_all));
ind_blank = intersect(find(stimCon_all==0),find(maskCon_all==0));

resp_cell = cell(nStimDir,2);
trialsperstim = zeros(nStimDir,2);
h_resp =zeros(nCells,nStimDir,2);
p_resp =zeros(nCells,nStimDir,2);
resp_win = prewin_frames+5:prewin_frames+nFramesOn;
resp_blank = squeeze(mean(data_dfof_tc(prewin_frames+frame_rate:prewin_frames+nFramesOn,:,ind_blank),1));
avg_resp_dir = zeros(nCells, nStimDir, 2, 2);
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
    trialsperstim(iDir,1) = length(ind_diralone);
    trialsperstim(iDir,2) = length(ind_dirplaid);
    resp_cell{iDir,1} = squeeze(mean(data_dfof_tc(resp_win,:,ind_diralone),1));
    resp_cell{iDir,2} = squeeze(mean(data_dfof_tc(resp_win,:,ind_dirplaid),1));
    [h_resp(:,iDir,1), p_resp(:,iDir,1)] = ttest2(resp_cell{iDir,1},resp_blank,'dim',2,'tail','right','alpha', 0.05./nStim);
    [h_resp(:,iDir,2), p_resp(:,iDir,2)] = ttest2(resp_cell{iDir,2},resp_blank,'dim',2,'tail','right','alpha', 0.05./nStim);
    avg_resp_dir(:,iDir,1,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind_diralone),1),3));
    avg_resp_dir(:,iDir,1,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,ind_diralone),1),[],3)./sqrt(length(ind_diralone)));
    avg_resp_dir(:,iDir,2,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind_dirplaid),1),3));
    avg_resp_dir(:,iDir,2,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,ind_dirplaid),1),[],3)./sqrt(length(ind_dirplaid)));
    all_resp_dir = [all_resp_dir squeeze(mean(data_dfof_tc(resp_win,:,ind_diralone),1))];
    all_resp_plaid = [all_resp_plaid squeeze(mean(data_dfof_tc(resp_win,:,ind_dirplaid),1))];
    all_dir = [all_dir iDir.*ones(size(ind_diralone))];
    all_plaid = [all_plaid iDir.*ones(size(ind_dirplaid))];
end

resp_ind = find(sum(sum(h_resp,2),3));
resp_ind_dir = find(sum(h_resp(:,:,1),2));
resp_ind_plaid = find(sum(h_resp(:,:,2),2));
p_anova_dir = zeros(1,nCells);
p_anova_plaid = zeros(1,nCells);
for iCell = 1:nCells
    p_anova_dir(iCell) = anova1(all_resp_dir(iCell,:), all_dir, 'off');
    p_anova_plaid(iCell) = anova1(all_resp_plaid(iCell,:), all_plaid, 'off');
end

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']), 'resp_cell', 'data_dfof_tc', 'resp_ind', 'frame_rate', 'h_resp', 'avg_resp_dir');
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']), 'prewin_frames', 'postwin_frames', 'resp_win', 'ind_stimAlone', 'ind_maskAlone', 'ind_plaid', 'ind_blank');


%%

b_dir = nan(1,nCells);
k1_dir = nan(1,nCells);
R1_dir = nan(1,nCells);
R2_dir = nan(1,nCells);
u1_dir = nan(1,nCells);
sse_dir = nan(1,nCells);
R_square_dir = nan(1,nCells);
range = 1:1:360;
y_dir_fit = nan(nCells,length(range));

figure; 
movegui('center')
start = 1;
n = 1;
for iC = 1:length(resp_ind)
    iCell = resp_ind(iC);
    if start>25
        sgtitle([date ' ' mouse ' Direction Tuning'])
        print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_plaidResp_dirTuning_' num2str(n) '.pdf']),'-dpdf', '-fillpage')
        n = n+1;
        figure; 
        movegui('center')
        start = 1;
    end
    subplot(5,5,start) 
    errorbar(stimDirs, avg_resp_dir(iCell,:,1,1), avg_resp_dir(iCell,:,1,2))
    hold on
    data = [avg_resp_dir(iCell,:,1,1) avg_resp_dir(iCell,1,1,1)];
    theta = [deg2rad(stimDirs) 2.*pi];
    [b_dir(:,iCell),k1_dir(:,iCell),R1_dir(:,iCell),R2_dir(:,iCell),u1_dir(:,iCell),u2_dir(:,iCell),sse_dir(:,iCell),R_square_dir(:,iCell)] ...
        = miaovonmisesfit_dir(theta,data);
    y_dir_fit(iCell,:) = b_dir(:,iCell)+R1_dir(:,iCell).*exp(k1_dir(:,iCell).*(cos(deg2rad(range)-u1_dir(:,iCell))-1))+R2_dir(:,iCell).*exp(k1_dir(:,iCell).*(cos(deg2rad(range)-u1_dir(:,iCell)-pi)-1));
    plot(range,y_dir_fit(iCell,:))
    errorbar(stimDirs, avg_resp_dir(iCell,:,2,1), avg_resp_dir(iCell,:,2,2))
    start = start+1;
    if p_anova_plaid(iCell)<0.05
        title('*')
    end
end
sgtitle([date ' ' mouse ' Direction Tuning'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_plaidResp_dirTuning_' num2str(n) '.pdf']),'-dpdf', '-fillpage')       

int = unique(diff(stimDirs));
component = avg_resp_dir(:,:,1,1)+circshift(avg_resp_dir(:,:,1,1),-90./int,2);
pattern = circshift(avg_resp_dir(:,:,1,1),-45./int,2);

comp_corr = zeros(1,nCells);
patt_corr = zeros(1,nCells);
comp_patt_corr = zeros(1,nCells);

for iCell = 1:nCells
    comp_corr(iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,2,1),component(iCell,:)));
    patt_corr(iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,2,1),pattern(iCell,:)));
    comp_patt_corr(iCell) = triu2vec(corrcoef(component(iCell,:),pattern(iCell,:)));
end
Rp = ((patt_corr)-(comp_corr.*comp_patt_corr))./sqrt((1-comp_corr.^2).*(1-comp_patt_corr.^2));
Rc = ((comp_corr)-(patt_corr.*comp_patt_corr))./sqrt((1-patt_corr.^2).*(1-comp_patt_corr.^2));
Zp = (0.5.*log((1+Rp)./(1-Rp)))./sqrt(1./(nStimDir-3));
Zc = (0.5.*log((1+Rc)./(1-Rc)))./sqrt(1./(nStimDir-3));

figure; 
movegui('center')
scatter(Zc(resp_ind), Zp(resp_ind))
xlabel('Zc')
ylabel('Zp')
ylim([-4 8])
xlim([-4 8])
hold on
plotZcZpBorders
sgtitle([date ' ' mouse])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_ZcZp.pdf']),'-dpdf', '-fillpage')       



ind = find(Zp>1);
figure; 
movegui('center')
start = 1;
[n n2] = subplotn(length(ind));
for iC = 1:length(ind)
    iCell = ind(iC);
    subplot(n,n2,start) 
    errorbar(stimDirs, avg_resp_dir(iCell,:,1,1), avg_resp_dir(iCell,:,1,2))
    hold on
    errorbar(stimDirs, avg_resp_dir(iCell,:,2,1), avg_resp_dir(iCell,:,2,2))
    start = start+1;
    title(['Zc = ' num2str(chop(Zc(iCell),2)) '; Zp= ' num2str(chop(Zp(iCell),2))])
end
sgtitle([date ' ' mouse ' Cells with Zp>1'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirTuning_Zpabove1.pdf']),'-dpdf', '-fillpage')       


ind = find(Zc>2);
figure; 
movegui('center')
start = 1;
[n n2] = subplotn(length(ind));
for iC = 1:length(ind)
    iCell = ind(iC);
    subplot(n,n2,start) 
    errorbar(stimDirs, avg_resp_dir(iCell,:,1,1), avg_resp_dir(iCell,:,1,2))
    hold on
    errorbar(stimDirs, avg_resp_dir(iCell,:,2,1), avg_resp_dir(iCell,:,2,2))
    start = start+1;
    title(['Zc = ' num2str(chop(Zc(iCell),2)) '; Zp= ' num2str(chop(Zp(iCell),2))])
end
sgtitle([date ' ' mouse ' Cells with Zc>2'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirTuning_Zcabove2.pdf']),'-dpdf', '-fillpage')       



stim_SI = zeros(1,nCells);
stim_OSI = zeros(1,nCells);
plaid_OSI = zeros(1,nCells);
stim_DSI = zeros(1,nCells);
plaid_DSI = zeros(1,nCells);
plaid_SI = zeros(1,nCells);
plaid_SI_45 = zeros(1,nCells);
h_plaid_SI = zeros(1,nCells);
stim_prefDir = zeros(1,nCells);
plaid_prefDir = zeros(1,nCells);
stim_prefOri = zeros(1,nCells);
plaid_prefOri = zeros(1,nCells);
resp_stim_prefDir = zeros(1,nCells);
resp_mask_prefDir = zeros(1,nCells);
resp_plaid_prefDir = zeros(1,nCells);
resp_stim_45Dir = zeros(1,nCells);
resp_mask_45Dir = zeros(1,nCells);
resp_plaid_45Dir = zeros(1,nCells);
avg_resp_ori = squeeze(mean(reshape(avg_resp_dir, [nCells nStimDir./2 2 2 2]),3));
avg_resp_ori_rect = avg_resp_ori;
avg_resp_ori_rect(find(avg_resp_ori<0)) = 0;
avg_resp_dir_rect = avg_resp_dir;
avg_resp_dir_rect(find(avg_resp_dir<0)) = 0;
nStimOri = size(avg_resp_ori,2);
for iCell = 1:nCells
    [max_val max_ind] = max(avg_resp_dir_rect(iCell,:,2,1));
    plaid_prefDir(iCell) = stimDirs(max_ind);
    null_ind = max_ind+(nStimDir./2);
    null_ind(find(null_ind>nStimDir)) = null_ind(find(null_ind>nStimDir))-nStimDir;
    min_val = avg_resp_dir_rect(iCell,null_ind,2,1);
    plaid_DSI(iCell) = (max_val-min_val)./(max_val+min_val);
    [max_val max_ind] = max(avg_resp_dir_rect(iCell,:,1,1));
    mask_ind = max_ind+(90./int);
    mask_ind(find(mask_ind>nStimDir)) = mask_ind(find(mask_ind>nStimDir))-nStimDir;
    plaid_val = avg_resp_dir_rect(iCell,max_ind,2,1);
    mask_val = avg_resp_dir_rect(iCell,mask_ind,1,1);
    plaid_SI(iCell) = (plaid_val-max_val-mask_val)./(plaid_val+max_val+mask_val);
    stim_SI(iCell) = (max_val-mask_val)./(max_val+mask_val);
    resp_stim_prefDir(iCell) = max_val;
    resp_mask_prefDir(iCell) = mask_val;
    resp_plaid_prefDir(iCell) = plaid_val;
    max_ind_45 = max_ind+2;
    max_ind_45(find(max_ind_45>nStimDir)) = max_ind_45(find(max_ind_45>nStimDir))-nStimDir;
    max_val_45 = avg_resp_dir_rect(iCell,max_ind_45,1,1);
    plaid_val_45 = avg_resp_dir_rect(iCell,max_ind_45,2,1);
    mask_ind_45 = mask_ind+2;
    mask_ind_45(find(mask_ind_45>nStimDir)) = mask_ind_45(find(mask_ind_45>nStimDir))-nStimDir;
    mask_val_45 = avg_resp_dir_rect(iCell,mask_ind_45,1,1);
    plaid_SI_45(iCell) = (plaid_val_45-max_val_45-mask_val_45)./(plaid_val_45+max_val_45+mask_val_45);
    resp_stim_45Dir(iCell) = max_val_45;
    resp_mask_45Dir(iCell) = mask_val_45;
    resp_plaid_45Dir(iCell) = plaid_val_45;
    h_plaid_SI(iCell) = ttest(resp_cell{max_ind,2}(iCell,:),max_val+mask_val);
    stim_prefDir(iCell) = stimDirs(max_ind);
    null_ind = max_ind+(nStimDir./2);
    null_ind(find(null_ind>nStimDir)) = null_ind(find(null_ind>nStimDir))-nStimDir;
    min_val = avg_resp_dir_rect(iCell,null_ind,1,1);
    stim_DSI(iCell) = (max_val-min_val)./(max_val+min_val);
    [max_val max_ind] = max(avg_resp_ori_rect(iCell,:,2,1));
    plaid_prefOri(iCell) = stimDirs(max_ind);
    null_ind = max_ind+(nStimOri./2);
    null_ind(find(null_ind>nStimOri)) = null_ind(find(null_ind>nStimDir/2))-nStimOri;
    min_val = avg_resp_ori_rect(iCell,null_ind,2,1);
    plaid_OSI(iCell) = (max_val-min_val)./(max_val+min_val);
    [max_val max_ind] = max(avg_resp_ori_rect(iCell,:,1,1));
    stim_prefOri(iCell) = stimDirs(max_ind);
    null_ind = max_ind+(nStimOri./2);
    null_ind(find(null_ind>nStimOri)) = null_ind(find(null_ind>nStimDir/2))-nStimOri;
    min_val = avg_resp_ori_rect(iCell,null_ind,1,1);
    stim_OSI(iCell) = (max_val-min_val)./(max_val+min_val);
end
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirAnalysis.mat']), 'component','pattern','Rp', 'Rc', 'Zp', 'Zc', 'stim_SI', 'stim_OSI', 'stim_DSI', 'plaid_OSI', 'plaid_DSI', 'plaid_SI', 'plaid_SI_45','h_plaid_SI', 'nCells','resp_stim_prefDir','resp_mask_prefDir','resp_plaid_prefDir','resp_stim_45Dir','resp_mask_45Dir','resp_plaid_45Dir', 'b_dir', 'k1_dir', 'R1_dir', 'R2_dir', 'u1_dir', 'sse_dir', 'R_square_dir', 'y_dir_fit');


figure; 
movegui('center')
subplot(3,3,1)
scatter(stim_DSI(resp_ind_plaid),plaid_SI(resp_ind_plaid))
xlabel('Stim DSI')
ylabel('SI')
subplot(3,3,2)
scatter(stim_DSI(resp_ind_plaid),Rc(resp_ind_plaid))
xlabel('Stim DSI')
ylabel('Rc')
ylim([-2 15])
subplot(3,3,3)
scatter(stim_DSI(resp_ind_plaid),Rp(resp_ind_plaid))
xlabel('Stim DSI')
ylabel('Rp')
ylim([-2 15])
subplot(3,3,4)
scatter(plaid_DSI(resp_ind_plaid),plaid_SI(resp_ind_plaid))
xlabel('Plaid DSI')
ylabel('SI')
subplot(3,3,5)
scatter(plaid_DSI(resp_ind_plaid),Rc(resp_ind_plaid))
xlabel('Plaid DSI')
ylabel('Rc')
ylim([-2 15])
subplot(3,3,6)
scatter(plaid_DSI(resp_ind_plaid),Rp(resp_ind_plaid))
xlabel('Plaid DSI')
ylabel('Rp')
ylim([-2 15])
subplot(3,3,7)
cdfplot(stim_DSI(resp_ind_plaid))
hold on
cdfplot(plaid_DSI(resp_ind_plaid))
legend({'Stim','Plaid'},'location','southeast')
title('')
xlabel('DSI')
subplot(3,3,8)
scatter(plaid_SI(resp_ind_plaid),Rc(resp_ind_plaid))
xlabel('SI')
ylabel('Rc')
ylim([-2 15])
subplot(3,3,9)
scatter(plaid_SI(resp_ind_plaid),Rp(resp_ind_plaid))
xlabel('SI')
ylabel('Rp')
ylim([-2 15])
sgtitle([date ' ' mouse])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_plaidResp_DSI.pdf']),'-dpdf', '-fillpage')       

figure; 
subplot(3,3,1)
movegui('center')
scatter(stim_OSI(resp_ind_plaid),plaid_SI(resp_ind_plaid))
xlabel('Stim OSI')
ylabel('SI')
subplot(3,3,2)
scatter(stim_OSI(resp_ind_plaid),Rc(resp_ind_plaid))
xlabel('Stim OSI')
ylabel('Rc')
ylim([-2 15])
subplot(3,3,3)
scatter(stim_OSI(resp_ind_plaid),Rp(resp_ind_plaid))
xlabel('Stim OSI')
ylabel('Rp')
ylim([-2 15])
subplot(3,3,4)
scatter(plaid_OSI(resp_ind_plaid),plaid_SI(resp_ind_plaid))
xlabel('Plaid OSI')
ylabel('SI')
xlim([0 1])
subplot(3,3,5)
scatter(plaid_OSI(resp_ind_plaid),Rc(resp_ind_plaid))
xlabel('Plaid OSI')
ylabel('Rc')
ylim([-2 15])
xlim([0 1])
subplot(3,3,6)
scatter(plaid_OSI(resp_ind_plaid),Rp(resp_ind_plaid))
xlabel('Plaid OSI')
ylabel('Rp')
ylim([-2 15])
xlim([0 1])
subplot(3,3,7)
cdfplot(stim_OSI(resp_ind_plaid))
hold on
cdfplot(plaid_OSI(resp_ind_plaid))
legend({'Stim','Plaid'})
xlabel('OSI')
title('')
sgtitle([date ' ' mouse])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_plaidResp_OSI.pdf']),'-dpdf', '-fillpage')       
end
