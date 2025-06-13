clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandDirTwoPhase_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = length(expt);

for iexp = 27
frame_rate = 15;

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc;
ImgFolder = expt(iexp).coFolder;
time = expt(iexp).coTime;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
%base = '\\CRASH.dhe.duke.edu\data\home\lindsey';

fprintf([mouse ' ' date '\n'])

%% Pref direction analysis
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))

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
ind_p = cell(1,nMaskPhas);
for ip = 1:nMaskPhas
    ind_p{ip} = find(maskPhas_all == maskPhas(ip));
end
resp_cell = cell(nStimDir,nMaskPhas,2);
trialsperstim = zeros(nStimDir,nMaskPhas,2);
h_resp =zeros(nCells,nStimDir,nMaskPhas,2);
p_resp =zeros(nCells,nStimDir,nMaskPhas,2);
resp_win = prewin_frames+5:prewin_frames+nFramesOn;
resp_blank = squeeze(mean(data_dfof_tc(prewin_frames+frame_rate:prewin_frames+nFramesOn,:,ind_blank),1));
avg_resp_dir = zeros(nCells, nStimDir,nMaskPhas, 2, 2);
all_resp_dir = [];
all_resp_plaid = cell(1,nMaskPhas);
all_dir = [];
all_plaid = cell(1,nMaskPhas);
nStim = nStimDir;
for iDir = 1:nStimDir
    ind_stimdir = find(stimDir_all == stimDirs(iDir));
    ind_maskdir = find(maskDir_all == stimDirs(iDir));
    ind_diralone = [intersect(ind_stimdir, ind_stimAlone), intersect(ind_maskdir, ind_maskAlone)];
    ind_dirplaid = [intersect(ind_stimdir, ind_plaid)];
    trialsperstim(iDir,1,1) = length(ind_diralone);
    resp_cell{iDir,1,1} = squeeze(mean(data_dfof_tc(resp_win,:,ind_diralone),1));
    [h_resp(:,iDir,1,1), p_resp(:,iDir,1,1)] = ttest2(resp_cell{iDir,1},resp_blank,'dim',2,'tail','right','alpha', 0.05./nStim);
    avg_resp_dir(:,iDir,1,1,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind_diralone),1),3));
    avg_resp_dir(:,iDir,1,1,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,ind_diralone),1),[],3)./sqrt(length(ind_diralone)));
    all_resp_dir = [all_resp_dir squeeze(mean(data_dfof_tc(resp_win,:,ind_diralone),1))];
    all_dir = [all_dir iDir.*ones(size(ind_diralone))];
    for ip = 1:nMaskPhas
        ind_dpplaid = intersect(ind_dirplaid,ind_p{ip});
        trialsperstim(iDir,ip,2) = length(ind_dpplaid);
        resp_cell{iDir,ip,2} = squeeze(mean(data_dfof_tc(resp_win,:,ind_dpplaid),1));
        [h_resp(:,iDir,ip,2), p_resp(:,iDir,ip,2)] = ttest2(resp_cell{iDir,ip,2},resp_blank,'dim',2,'tail','right','alpha', 0.05./nStim);
        avg_resp_dir(:,iDir,ip,2,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind_dpplaid),1),3));
        avg_resp_dir(:,iDir,ip,2,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,ind_dpplaid),1),[],3)./sqrt(length(ind_dpplaid)));
        if iDir == 1
            all_resp_plaid{ip} = [];
            all_plaid{ip} = [];
        end
        all_resp_plaid{ip} = [all_resp_plaid{ip} squeeze(mean(data_dfof_tc(resp_win,:,ind_dpplaid),1))];
        all_plaid{ip} = [all_plaid{ip} iDir.*ones(size(ind_dpplaid))];
    end
end

resp_ind = find(sum(sum(sum(h_resp,2),3),4));
resp_ind_dir = find(sum(h_resp(:,:,1,1),2));
resp_ind_plaid = find(sum(sum(h_resp(:,:,:,2),2),3));
p_anova_dir = zeros(1,nCells);
p_anova_plaid = zeros(2,nCells);
for iCell = 1:nCells
    p_anova_dir(iCell) = anova1(all_resp_dir(iCell,:), all_dir, 'off');
    for ip = 1:nMaskPhas
        p_anova_plaid(ip,iCell) = anova1(all_resp_plaid{ip}(iCell,:), all_plaid{ip}, 'off');
    end
end

save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']), 'resp_cell', 'data_dfof_tc', 'resp_ind', 'frame_rate', 'h_resp', 'avg_resp_dir','p_anova_dir','p_anova_plaid');
save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']), 'prewin_frames', 'postwin_frames', 'resp_win', 'ind_stimAlone', 'ind_maskAlone', 'ind_plaid', 'ind_blank');


%%
% 
% figure; 
% movegui('center')
% start = 1;
% n = 1;
% for iC = 1:length(resp_ind)
%     iCell = resp_ind(iC);
%     if start>25
%         suptitle([date ' ' mouse ' Direction Tuning'])
%         print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_plaidResp_dirTuning_' num2str(n) '.pdf']),'-dpdf', '-fillpage')
%         n = n+1;
%         figure; 
%         movegui('center')
%         start = 1;
%     end
%     subplot(5,5,start) 
%     errorbar(stimDirs, avg_resp_dir(iCell,:,1,1,1), avg_resp_dir(iCell,:,1,1,2))
%     hold on
%     errorbar(stimDirs, avg_resp_dir(iCell,:,1,2,1), avg_resp_dir(iCell,:,1,2,2))
%      errorbar(stimDirs, avg_resp_dir(iCell,:,2,2,1), avg_resp_dir(iCell,:,2,2,2))
%     start = start+1;
% end
% suptitle([date ' ' mouse ' Direction Tuning'])
% print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_plaidResp_dirTuning_' num2str(n) '.pdf']),'-dpdf', '-fillpage')       


%%

avg_resp_dir_rand = zeros(nCells,nStimDir,2);
for i = 1:nStimDir
    n = size(resp_cell{i,1,2},2);
    ind1 = randsample(1:n,ceil(n/2));
    ind2 = setdiff(1:n,ind1);
    avg_resp_dir_rand(:,i,1) = mean(resp_cell{i,1,2}(:,ind1),2);
    avg_resp_dir_rand(:,i,2) = mean(resp_cell{i,1,2}(:,ind2),2);
end
    
int = unique(diff(stimDirs));
component = avg_resp_dir(:,:,1,1,1)+circshift(avg_resp_dir(:,:,1,1,1),-90./int,2);
pattern = circshift(avg_resp_dir(:,:,1,1,1),-45./int,2);

comp_corr = zeros(nMaskPhas,nCells);
patt_corr = zeros(nMaskPhas,nCells);
comp_patt_corr = zeros(nMaskPhas,nCells);
plaid_corr = zeros(1,nCells);
plaid_corr_rand = zeros(1,nCells);

for iCell = 1:nCells
    for ip = 1:nMaskPhas
        comp_corr(ip,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,ip,2,1),component(iCell,:)));
        patt_corr(ip,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,ip,2,1),pattern(iCell,:)));
        comp_patt_corr(ip,iCell) = triu2vec(corrcoef(component(iCell,:),pattern(iCell,:)));
    end
    plaid_corr(1,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,1,2,1),avg_resp_dir(iCell,:,2,2,1)));
    plaid_corr_rand(1,iCell) = triu2vec(corrcoef(avg_resp_dir_rand(iCell,:,1),avg_resp_dir_rand(iCell,:,2)));
end
Rp = ((patt_corr)-(comp_corr.*comp_patt_corr))./sqrt((1-comp_corr.^2).*(1-comp_patt_corr.^2));
Rc = ((comp_corr)-(patt_corr.*comp_patt_corr))./sqrt((1-patt_corr.^2).*(1-comp_patt_corr.^2));
Zp = (0.5.*log((1+Rp)./(1-Rp)))./sqrt(1./(nStimDir-3));
Zc = (0.5.*log((1+Rc)./(1-Rc)))./sqrt(1./(nStimDir-3));

ZcZp_diff = Zc-Zp;
ind1 = intersect(find(Zp(1,:)>1.28),find(Zp(1,:)-Zc(1,:)>1.28));
ind2 = intersect(find(Zp(2,:)>1.28),find(Zp(2,:)-Zc(2,:)>1.28));
figure; 
movegui('center')
subplot(2,2,1)
scatter(Zc(1,resp_ind), Zp(1,resp_ind))
hold on
scatter(Zc(1,ind1),Zp(1,ind1));
scatter(Zc(1,ind2),Zp(1,ind2));
xlabel('Zc')
ylabel('Zp')
ylim([-4 8])
xlim([-4 8])
hold on
plotZcZpBorders
subplot(2,2,2)
scatter(Zc(2,resp_ind), Zp(2,resp_ind))
hold on
scatter(Zc(2,ind1),Zp(2,ind1));
scatter(Zc(2,ind2),Zp(2,ind2));
xlabel('Zc')
ylabel('Zp')
ylim([-4 8])
xlim([-4 8])
hold on
plotZcZpBorders
subplot(2,2,3)
scatter(ZcZp_diff(1,resp_ind), ZcZp_diff(2,resp_ind))
xlabel('Zc-Zp (0 deg)')
ylabel('Zc-Zp (90 deg)')
ylim([-4 8])
xlim([-4 8])
refline(1)
subplot(2,2,4)
histogram(plaid_corr)
hold on
histogram(plaid_corr_rand)
xlim([-1 1])
xlabel('Correlation of 0/90 deg plaids')
legend({'across','within'},'location','northwest')
subtitle([date ' ' mouse])
print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_ZcZp.pdf']),'-dpdf', '-fillpage')       

%%

plaid_SI = zeros(2,nCells);
plaid_SI_45 = zeros(2,nCells);
h_plaid_SI = zeros(2,nCells);
stim_prefDir = zeros(1,nCells);
stim_prefOri = zeros(1,nCells);
resp_stim_prefDir = zeros(1,nCells);
resp_mask_prefDir = zeros(1,nCells);
resp_plaid_prefDir = zeros(2,nCells);
resp_stim_45Dir = zeros(1,nCells);
resp_mask_45Dir = zeros(1,nCells);
resp_plaid_45Dir = zeros(2,nCells);
avg_resp_ori = squeeze(mean(reshape(avg_resp_dir, [nCells nStimDir./2 2 2 2 2]),3));
avg_resp_ori_rect = avg_resp_ori;
avg_resp_ori_rect(find(avg_resp_ori<0)) = 0;
avg_resp_dir_rect = avg_resp_dir;
avg_resp_dir_rect(find(avg_resp_dir<0)) = 0;
nStimOri = size(avg_resp_ori,2);
stim_DSI = zeros(1,nCells);
stim_OSI = zeros(1,nCells);

for iCell = 1:nCells
    [max_val max_ind] = max(avg_resp_dir_rect(iCell,:,1,1,1));
    opp_ind = max_ind+(nStimDir/2);
    opp_ind(find(opp_ind>nStimDir)) = opp_ind(find(opp_ind>nStimDir))-nStimDir;
    opp_val = avg_resp_dir_rect(iCell,opp_ind,1,1,1);
    stim_DSI(iCell) = (max_val-opp_val)./(max_val+opp_val);
    mask_ind = max_ind+(90./int);
    mask_ind(find(mask_ind>nStimDir)) = mask_ind(find(mask_ind>nStimDir))-nStimDir;
    plaid_val = squeeze(avg_resp_dir_rect(iCell,max_ind,:,2,1));
    mask_val = avg_resp_dir_rect(iCell,mask_ind,1,1,1);
    plaid_SI(:,iCell) = (plaid_val-max_val-mask_val)./(plaid_val+max_val+mask_val);
    resp_stim_prefDir(1,iCell) = max_val;
    resp_mask_prefDir(1,iCell) = mask_val;
    resp_plaid_prefDir(:,iCell) = plaid_val;
    max_ind_45 = max_ind+2;
    max_ind_45(find(max_ind_45>nStimDir)) = max_ind_45(find(max_ind_45>nStimDir))-nStimDir;
    max_val_45 = avg_resp_dir_rect(iCell,max_ind_45,1,1,1);
    plaid_val_45 = squeeze(avg_resp_dir_rect(iCell,max_ind_45,:,2,1));
    mask_ind_45 = mask_ind+2;
    mask_ind_45(find(mask_ind_45>nStimDir)) = mask_ind_45(find(mask_ind_45>nStimDir))-nStimDir;
    mask_val_45 = avg_resp_dir_rect(iCell,mask_ind_45,1,1,1);
    plaid_SI_45(:,iCell) = (plaid_val_45-max_val_45-mask_val_45)./(plaid_val_45+max_val_45+mask_val_45);
    resp_stim_45Dir(1,iCell) = max_val_45;
    resp_mask_45Dir(1,iCell) = mask_val_45;
    resp_plaid_45Dir(:,iCell) = plaid_val_45;
    h_plaid_SI(1,iCell) = ttest(resp_cell{max_ind,1,2}(iCell,:),max_val+mask_val);
    h_plaid_SI(2,iCell) = ttest(resp_cell{max_ind,2,2}(iCell,:),max_val+mask_val);
    [max_val max_ind] = max(avg_resp_ori_rect(iCell,:,1,1));
    stim_prefOri(iCell) = stimDirs(max_ind);
    null_ind = max_ind+(nStimOri./2);
    null_ind(find(null_ind>nStimOri)) = null_ind(find(null_ind>nStimDir/2))-nStimOri;
    min_val = avg_resp_ori_rect(iCell,null_ind,1,1);
    stim_OSI(iCell) = (max_val-min_val)./(max_val+min_val);
end

save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirAnalysis.mat']), 'component','pattern','Rp', 'Rc', 'Zp', 'Zc', 'plaid_SI', 'plaid_SI_45','h_plaid_SI', 'nCells','resp_stim_prefDir','resp_mask_prefDir','resp_plaid_prefDir','resp_stim_45Dir','resp_mask_45Dir','resp_plaid_45Dir','plaid_corr','plaid_corr_rand','stim_OSI','stim_DSI');

end