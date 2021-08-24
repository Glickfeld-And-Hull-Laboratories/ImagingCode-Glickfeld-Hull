clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandDirFF_ExptList';
rc = behavConstsAV;
eval(ds)
nexp = length(expt);

for iexp = 1:nexp

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
trialsperstim = zeros(nStimDir,nMaskDiff+1);
h_resp =zeros(nCells,nStimDir,nMaskDiff+1);
p_resp =zeros(nCells,nStimDir,nMaskDiff+1);
resp_win = prewin_frames+5:prewin_frames+nFramesOn;
resp_blank = squeeze(mean(data_dfof_tc(prewin_frames+frame_rate:prewin_frames+nFramesOn,:,ind_blank),1));
avg_resp_dir = zeros(nCells, nStimDir, nMaskDiff+1, 2);
all_resp_dir = [];
all_resp_plaid = cell(1,nMaskDiff);
all_dir = [];
all_plaid = cell(1,nMaskDiff);
nStim = nStimDir;
for iDir = 1:nStimDir
    ind_stimdir = find(stimDir_all == stimDirs(iDir));
    ind_maskdir = find(maskDir_all == maskDirs(iDir));
    ind_diralone = [intersect(ind_stimdir, ind_stimAlone), intersect(ind_maskdir, ind_maskAlone)];
    trialsperstim(iDir,1) = length(ind_diralone);
    resp_cell{iDir,1} = squeeze(mean(data_dfof_tc(resp_win,:,ind_diralone),1));
	[h_resp(:,iDir,1), p_resp(:,iDir,1)] = ttest2(resp_cell{iDir,1},resp_blank,'dim',2,'tail','right','alpha', 0.05./nStim);
    avg_resp_dir(:,iDir,1,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind_diralone),1),3));
    avg_resp_dir(:,iDir,1,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,ind_diralone),1),[],3)./sqrt(length(ind_diralone)));
    all_resp_dir = [all_resp_dir squeeze(mean(data_dfof_tc(resp_win,:,ind_diralone),1))];
    all_dir = [all_dir iDir.*ones(size(ind_diralone))];
    for id = 1:nMaskDiff
        ind_dirplaid = [intersect(find(maskDiff_all == maskDiffs(id)),intersect(ind_stimdir, ind_plaid))];
        trialsperstim(iDir,1+id) = length(ind_dirplaid);
        resp_cell{iDir,1+id} = squeeze(mean(data_dfof_tc(resp_win,:,ind_dirplaid),1));
        [h_resp(:,iDir,1+id), p_resp(:,iDir,1+id)] = ttest2(resp_cell{iDir,2},resp_blank,'dim',2,'tail','right','alpha', 0.05./nStim);
        avg_resp_dir(:,iDir,1+id,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind_dirplaid),1),3));
        avg_resp_dir(:,iDir,1+id,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,ind_dirplaid),1),[],3)./sqrt(length(ind_dirplaid)));
        all_resp_plaid{id} = [all_resp_plaid{id} squeeze(mean(data_dfof_tc(resp_win,:,ind_dirplaid),1))];
        all_plaid{id} = [all_plaid{id} iDir.*ones(size(ind_dirplaid))];
    end
end

resp_ind = find(sum(sum(h_resp,2),3));
resp_ind_dir = find(sum(h_resp(:,:,1),2));
resp_ind_plaid = find(sum(sum(h_resp(:,:,2:end),2),3));
p_anova_dir = zeros(1,nCells);
p_anova_plaid = cell(1,nMaskDiff);
for iCell = 1:nCells
    p_anova_dir(iCell) = anova1(all_resp_dir(iCell,:), all_dir, 'off');
    for id = 1:nMaskDiff
        p_anova_plaid{id} = zeros(1,nCells);
        p_anova_plaid{id}(iCell) = anova1(all_resp_plaid{id}(iCell,:), all_plaid{id}, 'off');
    end
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
        suptitle([date ' ' mouse ' Direction Tuning'])
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
    for id = 1:nMaskDiff
        errorbar(stimDirs, avg_resp_dir(iCell,:,1+id,1), avg_resp_dir(iCell,:,1+id,2))
    end
    start = start+1;
end
suptitle([date ' ' mouse ' Direction Tuning'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_plaidResp_dirTuning_' num2str(n) '.pdf']),'-dpdf', '-fillpage')       

int = 15;
component = zeros(length(0:15:360),nCells,nMaskDiff);
pattern = zeros(length(0:15:360),nCells,nMaskDiff);
comp_corr = zeros(nCells,nMaskDiff);
patt_corr = zeros(nCells,nMaskDiff);
comp_patt_corr = zeros(nCells,nMaskDiff);
Rc = zeros(nCells,nMaskDiff);
Rp = zeros(nCells,nMaskDiff);
Zc = zeros(nCells,nMaskDiff);
Zp = zeros(nCells,nMaskDiff);
avg_resp_dir_interp = zeros(length(0:15:360),nCells,nMaskDiff+1);
for i = 1:nMaskDiff+1
   avg_resp_dir_interp(:,:,i) = interp1([0:30:360]',[avg_resp_dir(:,:,i,1) avg_resp_dir(:,1,i,1)]',0:15:360);
end
figure(1)
for id = 1:nMaskDiff
    component(:,:,id) = avg_resp_dir_interp(:,:,1)+circshift(avg_resp_dir_interp(:,:,1),-maskDiffs(id)./int,2);
    pattern(:,:,id) = circshift(avg_resp_dir_interp(:,:,1),-maskDiffs(id)./2./int,2);
    for iCell = 1:nCells
        comp_corr(iCell,id) = triu2vec(corrcoef(avg_resp_dir_interp(:,iCell,id+1),component(:,iCell,id)));
        patt_corr(iCell,id) = triu2vec(corrcoef(avg_resp_dir_interp(:,iCell,id+1),pattern(:,iCell,id)));
        comp_patt_corr(iCell,id) = triu2vec(corrcoef(component(:,iCell,id),pattern(:,iCell,id)));
    end
    Rp(:,id) = ((patt_corr(:,id))-(comp_corr(:,id).*comp_patt_corr(:,id)))./sqrt((1-comp_corr(:,id).^2).*(1-comp_patt_corr(:,id).^2));
    Rc(:,id) = ((comp_corr(:,id))-(patt_corr(:,id).*comp_patt_corr(:,id)))./sqrt((1-patt_corr(:,id).^2).*(1-comp_patt_corr(:,id).^2));
    Zp(:,id) = (0.5.*log((1+Rp(:,id))./(1-Rp(:,id))))./sqrt(1./(nStimDir-3));
    Zc(:,id) = (0.5.*log((1+Rc(:,id))./(1-Rc(:,id))))./sqrt(1./(nStimDir-3));

    movegui('center')
    subplot(2,2,id)
    scatter(Zc(resp_ind,id), Zp(resp_ind,id))
    xlabel('Zc')
    ylabel('Zp')
    ylim([-4 8])
    xlim([-4 8])
    axis square
    title(num2str(maskDiffs(id)))
    hold on
    plotZcZpBorders
end
subplot(2,2,4)
 scatter(mean(Zc(resp_ind,:),2), mean(Zp(resp_ind,:),2))
xlabel('Zc')
ylabel('Zp')
ylim([-4 8])
xlim([-4 8])
axis square
title('All')
hold on
plotZcZpBorders
suptitle([date ' ' mouse])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_ZcZp.pdf']),'-dpdf', '-fillpage')       

Zp_ind = find(mean(Zp,2)>1.28 & mean(Zp,2)-mean(Zc,2)>1.28);
[n n2] = subplotn(length(Zp_ind));
figure;
for i = 1:length(Zp_ind)
    subplot(n,n2,i)
    plot(0:15:360,avg_resp_dir_interp(:,Zp_ind(i),1))
    hold on
    for id = 1:nMaskDiff
        plot(0:15:360,circshift(avg_resp_dir_interp(:,Zp_ind(i),id+1),maskDiffs(id)/2/15))
    end
    title(['Zc- ' num2str(chop(mean(Zc(Zp_ind(i),:),2),2)) '; Zp- ' num2str(chop(mean(Zp(Zp_ind(i),:),2),2))])
end

Zc_ind = find(mean(Zc,2)>1.28 & mean(Zc,2)-mean(Zp,2)>1.28);
[n n2] = subplotn(length(Zc_ind));
figure;
for i = 1:length(Zc_ind)
    subplot(n,n2,i)
    plot(0:15:360,avg_resp_dir_interp(:,Zc_ind(i),1))
    hold on
    for id = 1:nMaskDiff
        plot(0:15:360,circshift(avg_resp_dir_interp(:,Zc_ind(i),id+1),maskDiffs(id)/2/15))
    end
    title(['Zc- ' num2str(chop(mean(Zc(Zc_ind(i),:),2),2)) '; Zp- ' num2str(chop(mean(Zp(Zc_ind(i),:),2),2))])
end
    


%     ind = find(Zp>1);
%     figure; 
%     movegui('center')
%     start = 1;
%     [n n2] = subplotn(length(ind));
%     for iC = 1:length(ind)
%         iCell = ind(iC);
%         subplot(n,n2,start) 
%         errorbar(stimDirs, avg_resp_dir(iCell,:,1,1), avg_resp_dir(iCell,:,1,2))
%         hold on
%         errorbar(stimDirs, avg_resp_dir(iCell,:,2,1), avg_resp_dir(iCell,:,2,2))
%         start = start+1;
%         title(['Zc = ' num2str(chop(Zc(iCell),2)) '; Zp= ' num2str(chop(Zp(iCell),2))])
%     end
%     suptitle([date ' ' mouse ' Cells with Zp>1'])
%     print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirTuning_Zpabove1.pdf']),'-dpdf', '-fillpage')       
% 
% 
%     ind = find(Zc>2);
%     figure; 
%     movegui('center')
%     start = 1;
%     [n n2] = subplotn(length(ind));
%     for iC = 1:length(ind)
%         iCell = ind(iC);
%         subplot(n,n2,start) 
%         errorbar(stimDirs, avg_resp_dir(iCell,:,1,1), avg_resp_dir(iCell,:,1,2))
%         hold on
%         errorbar(stimDirs, avg_resp_dir(iCell,:,2,1), avg_resp_dir(iCell,:,2,2))
%         start = start+1;
%         title(['Zc = ' num2str(chop(Zc(iCell),2)) '; Zp= ' num2str(chop(Zp(iCell),2))])
%     end
%     suptitle([date ' ' mouse ' Cells with Zc>2'])
%     print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirTuning_Zcabove2.pdf']),'-dpdf', '-fillpage')       



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
suptitle([date ' ' mouse])
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
suptitle([date ' ' mouse])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_plaidResp_OSI.pdf']),'-dpdf', '-fillpage')       
end
