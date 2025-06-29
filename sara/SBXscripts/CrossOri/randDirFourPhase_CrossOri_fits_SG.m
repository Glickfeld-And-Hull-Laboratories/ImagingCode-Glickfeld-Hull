clc; clear all; close all;
doRedChannel = 1;
ds = 'CrossOriRandDirFourPhase_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = length(expt);
frame_rate = 15;
seed = rng;

max_dist = 10;

for iexp = [122] %63 64 107 109

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc;
ImgFolder = expt(iexp).coFolder;
time = expt(iexp).coTime;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
LGbase = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

fprintf([mouse ' ' date '\n'])

%% Pref direction analysis
load(fullfile(LGbase, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
load(fullfile(LGbase, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
load(fullfile(LGbase, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']))

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

ind_stimAlone = intersect(intersect(find(stimCon_all),find(maskCon_all==0)),find(centroid_dist<max_dist));
ind_maskAlone = intersect(intersect(find(stimCon_all==0),find(maskCon_all)),find(centroid_dist<max_dist));
ind_plaid = intersect(intersect(find(stimCon_all),find(maskCon_all)),find(centroid_dist<max_dist));
ind_blank = intersect(find(stimCon_all==0),find(maskCon_all==0));
ind_p = cell(1,nMaskPhas);
% maskPhas_shuf = maskPhas_all(randperm(length(maskPhas_all)));
% ind_p_shuf = cell(1,nMaskPhas);
for ip = 1:nMaskPhas
    ind_p{ip} = intersect(find(maskPhas_all == maskPhas(ip)),find(centroid_dist<max_dist));
%     ind_p_shuf{ip} = intersect(find(maskPhas_shuf == maskPhas(ip)),find(centroid_dist<max_dist));
end

resp_cell = cell(nStimDir,nMaskPhas,2);
trialsperstim = zeros(nStimDir,nMaskPhas,2);
trialInd = cell(nStimDir,nMaskPhas,2);
h_resp =zeros(nCells,nStimDir,nMaskPhas,2);
p_resp =zeros(nCells,nStimDir,nMaskPhas,2);
resp_win = prewin_frames+5:prewin_frames+nFramesOn;
resp_blank = squeeze(mean(data_dfof_tc(prewin_frames+frame_rate:prewin_frames+nFramesOn,:,ind_blank),1));
avg_resp_dir = zeros(nCells, nStimDir,nMaskPhas, 2, 2);
% avg_resp_dir_shuf = zeros(nCells, nStimDir,nMaskPhas, 2, 2);
all_resp_dir = [];
all_resp_plaid = cell(1,nMaskPhas);
all_dir = [];
all_plaid = cell(1,nMaskPhas);
nStim = nStimDir;
for iDir = 1:nStimDir
    ind_stimdir = intersect(find(stimDir_all == stimDirs(iDir)),find(centroid_dist<max_dist));
    ind_maskdir = intersect(find(maskDir_all == stimDirs(iDir)),find(centroid_dist<max_dist));
    ind_diralone = [intersect(ind_stimdir, ind_stimAlone), intersect(ind_maskdir, ind_maskAlone)];
    ind_dirplaid = [intersect(ind_stimdir, ind_plaid)];
    trialsperstim(iDir,1,1) = length(ind_diralone);
    trialInd{iDir,1,1} = ind_diralone;
    resp_cell{iDir,1,1} = squeeze(mean(data_dfof_tc(resp_win,:,ind_diralone),1));
    [h_resp(:,iDir,1,1), p_resp(:,iDir,1,1)] = ttest2(resp_cell{iDir,1,1},resp_blank,'dim',2,'tail','right','alpha', 0.05./nStim);
    avg_resp_dir(:,iDir,1,1,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind_diralone),1),3));
    avg_resp_dir(:,iDir,1,1,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,ind_diralone),1),[],3)./sqrt(length(ind_diralone)));
    all_resp_dir = [all_resp_dir squeeze(mean(data_dfof_tc(resp_win,:,ind_diralone),1))];
    all_dir = [all_dir iDir.*ones(size(ind_diralone))];
    for ip = 1:nMaskPhas
        ind_dpplaid = intersect(ind_dirplaid,ind_p{ip});
%         ind_dpplaid_shuf = intersect(ind_dirplaid,ind_p_shuf{ip});
        trialsperstim(iDir,ip,2) = length(ind_dpplaid);
        trialInd{iDir,ip,2} = ind_dpplaid;
        resp_cell{iDir,ip,2} = squeeze(mean(data_dfof_tc(resp_win,:,ind_dpplaid),1));
        [h_resp(:,iDir,ip,2), p_resp(:,iDir,ip,2)] = ttest2(resp_cell{iDir,ip,2},resp_blank,'dim',2,'tail','right','alpha', 0.05./nStim);
        avg_resp_dir(:,iDir,ip,2,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind_dpplaid),1),3));
        avg_resp_dir(:,iDir,ip,2,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,ind_dpplaid),1),[],3)./sqrt(length(ind_dpplaid)));
%         avg_resp_dir_shuf(:,iDir,ip,2,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind_dpplaid_shuf),1),3));
        if iDir == 1
            all_resp_plaid{ip} = [];
            all_plaid{ip} = [];
        end
        all_resp_plaid{ip} = [all_resp_plaid{ip} squeeze(mean(data_dfof_tc(resp_win,:,ind_dpplaid),1))];
        all_plaid{ip} = [all_plaid{ip} iDir.*ones(size(ind_dpplaid))];
    end
end

%shuffle trials across phase
ntrial = {};
resps = double.empty(nCells,0);
avg_resp_dir_shuf = zeros(nCells, nStimDir, nMaskPhas);
for id = 1:nStimDir
    for ip = 1:nMaskPhas
        ntrial{id,ip} = size(resp_cell{id,ip,2},2);
        resps = [resps, resp_cell{id,ip,2}];
    end
    resps_shuf = resps(:,randperm(size(resps,2)));
    start = 1;
    t=0;
    for ip = 1:nMaskPhas
        t = t + ntrial{id,ip};
        resp_cell_shuf{id,ip} = resps_shuf(:,start:t);
        avg_resp_dir_shuf(:,id,ip) =  mean(resp_cell_shuf{id,ip}(:,:),2);
    end
end



resp_ind = find(sum(sum(sum(h_resp,2),3),4));
resp_ind_dir = find(sum(h_resp(:,:,1,1),2)); %sig responsive to gratings
resp_ind_plaid = find(sum(sum(h_resp(:,:,:,2),2),3));
p_anova_dir = zeros(1,nCells);
p_anova_plaid = zeros(2,nCells);
for iCell = 1:nCells
    p_anova_dir(iCell) = anova1(all_resp_dir(iCell,:), all_dir, 'off'); %direction selective to gratings
    for ip = 1:nMaskPhas
        p_anova_plaid(ip,iCell) = anova1(all_resp_plaid{ip}(iCell,:), all_plaid{ip}, 'off'); %direction selective to plaids
    end
end


    p_dir = find(p_anova_dir<0.05);
    p_plaid1 = find(p_anova_plaid(1,:)<0.05);
    p_plaid2 = find(p_anova_plaid(2,:)<0.05);
    p_plaid3 = find(p_anova_plaid(3,:)<0.05);
    p_plaid4 = find(p_anova_plaid(4,:)<0.05);
    p_all = unique([p_dir,p_plaid1,p_plaid2,p_plaid3,p_plaid4]); %significantly responsive to a direction (anova) for gratings or any plaid set
    p_dir = find(p_anova_dir<0.05);
    p_plaid1 = find(p_anova_plaid(1,:)<0.05);
    p_all = unique([p_dir,p_plaid1]); %significantly responsive to a direction (anova) for gratings or any plaid set


for iCell = 1:nCells
    [max_val max_ind] = max(avg_resp_dir(iCell,:,1,1,1));
    null_ind = max_ind+(nStimDir./2);
    null_ind(find(null_ind>nStimDir)) = null_ind(find(null_ind>nStimDir))-nStimDir;
    min_val = avg_resp_dir(iCell,null_ind,1,1,1);
    if min_val < 0; min_val = 0; end
    DSI(iCell) = (max_val-min_val)./(max_val+min_val);
    DSI_maxInd(iCell) = max_ind; 
end

DSI_ind = find(DSI>0.5); %direction selective to gratings
OSI_ind = find(DSI<0.5);


%Calculate "reliability" of preferred grating direction
DSI_maxInd_boot = bootstrap_fourphase_randsamptrials(resp_cell,100)';

for iCell = 1:nCells
    idx = find(DSI_maxInd_boot(:,iCell)==DSI_maxInd(iCell));
    prefDir_resamp(iCell) = length(idx);
end

if ~exist(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' max_dist]),'dir')
    mkdir(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' num2str(max_dist)]))
end
save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_respData.mat']), 'resp_cell', 'data_dfof_tc', 'resp_ind', 'frame_rate', 'h_resp', 'avg_resp_dir','p_anova_dir','p_anova_plaid', 'trialInd');
save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_stimData.mat']), 'prewin_frames', 'postwin_frames', 'resp_win', 'ind_stimAlone', 'ind_maskAlone', 'ind_plaid', 'ind_blank','trialsperstim','DSI_ind','OSI_ind','resp_ind_dir','p_dir','prefDir_resamp');


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
component = avg_resp_dir(:,:,1,1,1)+circshift(avg_resp_dir(:,:,1,1,1),-120./int,2);
pattern = circshift(avg_resp_dir(:,:,1,1,1),-60./int,2);

%Calculate pattern and component prediction
comp_corr = zeros(nMaskPhas,nCells);
patt_corr = zeros(nMaskPhas,nCells);
comp_patt_corr = zeros(nMaskPhas,nCells);
plaid_corr = zeros(1,nCells);
plaid_corr1 = zeros(1,nCells);
plaid_corr2 = zeros(1,nCells);
plaid_corr3 = zeros(1,nCells);
plaid_corr4 = zeros(1,nCells);
plaid_corr5 = zeros(1,nCells);
plaid_corr6 = zeros(1,nCells);
plaid_corr_rand = zeros(1,nCells);

for iCell = 1:nCells
    for ip = 1:nMaskPhas
        comp_corr(ip,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,ip,2,1),component(iCell,:)));
        patt_corr(ip,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,ip,2,1),pattern(iCell,:)));
        comp_patt_corr(ip,iCell) = triu2vec(corrcoef(component(iCell,:),pattern(iCell,:)));
    end
    plaid_corr1(1,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,1,2,1),avg_resp_dir(iCell,:,2,2,1)));
    plaid_corr2(1,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,1,2,1),avg_resp_dir(iCell,:,3,2,1)));
    plaid_corr3(1,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,1,2,1),avg_resp_dir(iCell,:,4,2,1)));
    plaid_corr4(1,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,2,2,1),avg_resp_dir(iCell,:,3,2,1)));
    plaid_corr5(1,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,2,2,1),avg_resp_dir(iCell,:,4,2,1)));
    plaid_corr6(1,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,3,2,1),avg_resp_dir(iCell,:,4,2,1)));
    plaid_corr(1,iCell) = (plaid_corr1(1,iCell)+plaid_corr2(1,iCell)+plaid_corr3(1,iCell)+plaid_corr4(1,iCell)+plaid_corr5(1,iCell)+plaid_corr6(1,iCell))/6;
    plaid_corr_rand(1,iCell) = triu2vec(corrcoef(avg_resp_dir_rand(iCell,:,1),avg_resp_dir_rand(iCell,:,2)));
end

Rp = ((patt_corr)-(comp_corr.*comp_patt_corr))./sqrt((1-comp_corr.^2).*(1-comp_patt_corr.^2));
Rc = ((comp_corr)-(patt_corr.*comp_patt_corr))./sqrt((1-patt_corr.^2).*(1-comp_patt_corr.^2));
Zp = (0.5.*log((1+Rp)./(1-Rp)))./sqrt(1./(nStimDir-3));
Zc = (0.5.*log((1+Rc)./(1-Rc)))./sqrt(1./(nStimDir-3));

ZcZp_diff = Zc-Zp;
ind1 = intersect(find(Zp(1,:)>1.28),find(Zp(1,:)-Zc(1,:)>1.28));
ind2 = intersect(find(Zp(2,:)>1.28),find(Zp(2,:)-Zc(2,:)>1.28));
ind3 = intersect(find(Zp(3,:)>1.28),find(Zp(3,:)-Zc(3,:)>1.28));
ind4 = intersect(find(Zp(4,:)>1.28),find(Zp(4,:)-Zc(4,:)>1.28));
ind = {ind1,ind2,ind3,ind4};
ind_arr = [ind1,ind2,ind3,ind4];
ind_pds = unique(ind_arr);

save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_plaidCorr.mat']), 'plaid_corr', 'plaid_corr1', 'plaid_corr2', 'plaid_corr3', 'plaid_corr4', 'plaid_corr5', 'plaid_corr6');


% %Calculate SHUFFLED pattern and component prediction
% int = unique(diff(stimDirs));
% component_sh = avg_resp_dir(:,:,1,1,1)+circshift(avg_resp_dir(:,:,1,1,1),-120./int,2);
% pattern_sh = circshift(avg_resp_dir(:,:,1,1,1),-60./int,2);
% 
% comp_corr_sh = zeros(nMaskPhas,nCells);
% patt_corr_sh = zeros(nMaskPhas,nCells);
% comp_patt_corr_sh = zeros(nMaskPhas,nCells);
% 
% for iCell = 1:nCells
%     for ip = 1:nMaskPhas
%         comp_corr_sh(ip,iCell) = triu2vec(corrcoef(avg_resp_dir_shuf(iCell,:,ip),component_sh(iCell,:)));
%         patt_corr_sh(ip,iCell) = triu2vec(corrcoef(avg_resp_dir_shuf(iCell,:,ip),pattern_sh(iCell,:)));
%         comp_patt_corr_sh(ip,iCell) = triu2vec(corrcoef(component(iCell,:),pattern_sh(iCell,:)));
%     end
% end
% Rp_sh = ((patt_corr_sh)-(comp_corr_sh.*comp_patt_corr_sh))./sqrt((1-comp_corr_sh.^2).*(1-comp_patt_corr_sh.^2));
% Rc_sh = ((comp_corr_sh)-(patt_corr_sh.*comp_patt_corr_sh))./sqrt((1-patt_corr_sh.^2).*(1-comp_patt_corr_sh.^2));
% Zp_sh = (0.5.*log((1+Rp_sh)./(1-Rp_sh)))./sqrt(1./(nStimDir-3));
% Zc_sh = (0.5.*log((1+Rc_sh)./(1-Rc_sh)))./sqrt(1./(nStimDir-3));
% 
% ZcZp_diff_sh = Zc_sh-Zp_sh;


%% Calculate pairwise euclidean distance between ZpZc pairs

ZpZc = [];
for iCell = 1:nCells
    ZpZc(:,1,iCell) = Zp(:,iCell);
    ZpZc(:,2,iCell) = Zc(:,iCell);
end

ZpZcPWdist = double.empty(6,0);
for iCell = 1:nCells
    ZpZcPWdist(:,iCell) = pdist(ZpZc(:,:,iCell));
end

save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_ZpZc_pairwiseDist.mat']), 'ZpZcPWdist');


%% Trials per stimulus condition

figure;
subplot(2,1,1)
    bar(squeeze(trialsperstim(:,:,1))) 
    title('Trials per stim - GRATINGS')
    ylabel('# of trials')
    xlabel('direction')
subplot(2,1,2)
    bar(squeeze(trialsperstim(:,:,2))) 
    title('Trials per stim - PLAIDS')
    ylabel('# of trials')
    xlabel('direction')
    legend('0','90','180','270')
    
   
print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_TrialsPerStim.pdf']),'-dpdf', '-fillpage')       

    

%% ZcZp of population (resp_ind)

figure; 
movegui('center')
for i = 1:4
    subplot(4,4,i)
    scatter(Zc(i,resp_ind), Zp(i,resp_ind),'.')
    hold on
    scatter(Zc(i,ind{1}),Zp(i,ind{1}),'.');
    xlabel('Zc')
    ylabel('Zp')
    ylim([-4 8])
    xlim([-4 8])
    hold on
    if i==1; title('pattern cells at 0'); end
    plotZcZpBorders
end
for i = 1:4
    subplot(4,4,i+4)
    scatter(Zc(i,resp_ind), Zp(i,resp_ind),'.')
    hold on
    scatter(Zc(i,ind{2}),Zp(i,ind{2}),'.');
    xlabel('Zc')
    ylabel('Zp')
    ylim([-4 8])
    xlim([-4 8])
    hold on
    if i==1; title('pattern cells at 90'); end
    plotZcZpBorders
end
for i = 1:4
    subplot(4,4,i+8)
    scatter(Zc(i,resp_ind), Zp(i,resp_ind),'.')
    hold on
    scatter(Zc(i,ind{3}),Zp(i,ind{3}),'.');
    xlabel('Zc')
    ylabel('Zp')
    ylim([-4 8])
    xlim([-4 8])
    hold on
    if i==1; title('pattern cells at 180'); end
    plotZcZpBorders
end
for i = 1:4
    subplot(4,4,i+12)
    scatter(Zc(i,resp_ind), Zp(i,resp_ind),'.')
    hold on
    scatter(Zc(i,ind{4}),Zp(i,ind{4}),'.');
    xlabel('Zc')
    ylabel('Zp')
    ylim([-4 8])
    xlim([-4 8])
    hold on
    if i==1; title('pattern cells at 270'); end
    plotZcZpBorders
end
sgtitle('Pattern direction selective cells at four phases')
print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_ZcZp.pdf']),'-dpdf', '-fillpage')       

%% Set responsive index

% ind = 1:nCells;
ind = intersect(intersect(resp_ind_dir,DSI_ind),p_dir); %resp to 1 grating, DSI>0.5, resp to one direction of gratings
% ind = ind1;

%% PCI modulation

% PCI = (Zp-Zc)./(Zp+Zc);
PCI = (Zp-Zc);
% PCI_high = find(PCI>100);
% PCI(PCI_high) = NaN;

figure;
for i = 1:nMaskPhas
    cdfplot(PCI(i,:));
    hold on
end


%fit sinusoid
phase = [0 90 180 270];
phase_range = 0:1:359;
figure;
start=1;
n=1;
for iCell = 1:nCells
    subplot(5,4,start)
        scatter(phase,PCI(:,iCell),'LineWidth',1.25);
        hold on
        [b_hat_all(iCell,1), amp_hat_all(iCell,1), per_hat_all(iCell,1),pha_hat_all(iCell,1),sse_all(iCell,1),R_square_all(iCell,1)] = sinefit_PCI(deg2rad(phase),PCI(:,iCell));
        yfit_all(iCell,:,1) = b_hat_all(iCell,1)+amp_hat_all(iCell,1).*(sin(2*pi*deg2rad(phase_range)./per_hat_all(iCell,1) + 2.*pi/pha_hat_all(iCell,1)));
        idx = iCell==ind; %does iCell equal any value in the index?
        if any(idx)  %if yes, plot in red. if not, plot in black
            plot(phase_range, yfit_all(iCell,:,1),'k'); 
            subtitle(['cell ' num2str(iCell) ', Rsq ' num2str(R_square_all(iCell),'%.3f'), ', SSE ' num2str(sse_all(iCell),'%.2f')],'fontweight','bold')
        else
            plot(phase_range, yfit_all(iCell,:,1),'k:');
            subtitle(['cell ' num2str(iCell) ', Rsq ' num2str(R_square_all(iCell),'%.2f'), ', SSE ' num2str(sse_all(iCell),'%.2f')])
        end
        ylabel('Zp-Zc'); xlabel('Mask phase'); ylim([-7 7])
        xlim([0 360]); xticks([0 180 360]); set(gca,'TickDir','out'); axis square
    start = start+1;
    if start >20
        sgtitle([mouse ' ' date ' PCI modulation across mask phase by cell'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_PCImodulation_' num2str(n) '.pdf']),'-dpdf', '-fillpage')       
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle([mouse ' ' date ' PCI modulation across mask phase by cell'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_PCImodulation_' num2str(n) '.pdf']), '-dpdf','-fillpage')
    end        
end
    save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_PatternCorrelationIndexFits.mat']), 'ind', 'Zp', 'Zc', 'PCI', 'yfit_all', 'b_hat_all', 'amp_hat_all', 'per_hat_all', 'pha_hat_all', 'sse_all', 'R_square_all','DSI')
    close all
    
%% Shuffled modulation

for iCell = 1:nCells
    realorder = PCI(:,iCell);
    randorder = randperm(4);
    PCI_sh(:,iCell) = realorder(randorder);
end

%fit sinusoid
phase = [0 90 180 270];
phase_range = 0:1:359;

for iCell = 1:nCells
    [b_hat_all(iCell,1), amp_hat_all(iCell,1), per_hat_all(iCell,1),pha_hat_all(iCell,1),sse_all(iCell,1),R_square_all(iCell,1)] = sinefit_PCI(deg2rad(phase),PCI_sh(:,iCell));
    yfit_all(iCell,:,1) = b_hat_all(iCell,1)+amp_hat_all(iCell,1).*(sin(2*pi*deg2rad(phase_range)./per_hat_all(iCell,1) + 2.*pi/pha_hat_all(iCell,1)));
end
    save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_PatternCorrelationIndexFits_Shuffled.mat']), 'yfit_all', 'b_hat_all', 'amp_hat_all', 'per_hat_all', 'pha_hat_all', 'sse_all', 'R_square_all')

%% Direction tuning curve FIT

y_fits = zeros(360,nCells);
dirs = deg2rad(0:1:359);
for iCell = 1:nCells
    [b_hat_all(iCell,1), k1_hat_all(iCell,1), R1_hat_all(iCell,1), R2_hat_all(iCell,1), u1_hat_all(iCell,1), u2_hat_all(iCell,1), sse_all(iCell,1),R_square_all(iCell,1)] = miaovonmisesfit_dir(deg2rad(stimDirs),avg_resp_dir(iCell,:,1,1,1));
    dir_yfit_all(:,iCell) = b_hat_all(iCell,1)+R1_hat_all(iCell,1).*exp(k1_hat_all(iCell,1).*(cos(dirs-u1_hat_all(iCell,1))-1))+R2_hat_all(iCell,1).*exp(k1_hat_all(iCell,1).*(cos(dirs-u1_hat_all(iCell,1))-1));
end


% global DSI
    angs = 0:30:330;
    eps_val = 1e-3;  % Small epsilon to prevent division by zero
    for j = 1:nCells
        amps        = avg_resp_dir(j,:,1,1,1); %our responses  
        amps(amps < 0) = 0;

        total_response = sum(amps); % Normalize factor
        if total_response < eps_val
            g_dsi(j) = NaN; % Assign NaN if total response is too small
            g_osi(j) = NaN;
            ang(j) = NaN;
            ang_ori(j) = NaN;
        continue
        end

        % Direction Selectivity Index (gDSI)
            vec_x = sum(cos(deg2rad(angs)) .* amps);
            vec_y = sum(sin(deg2rad(angs)) .* amps);
            g_dsi(j) = sqrt(vec_x^2 + vec_y^2) / total_response;
    
        % Orientation Selectivity Index (gOSI)
            vec_x_ori = sum(cos(deg2rad(2 * angs)) .* amps);
            vec_y_ori = sum(sin(deg2rad(2 * angs)) .* amps);
            g_osi(j) = sqrt(vec_x_ori^2 + vec_y_ori^2) / total_response;
        
        % Preferred direction (in degrees)
            ang_dir(j) = mod(rad2deg(atan2(vec_y, vec_x)), 360);
        
        % Preferred orientation (in degrees, folded into [0, 180))
            ang_ori(j) = mod(0.5 * rad2deg(atan2(vec_y_ori, vec_x_ori)), 180);
    end

    save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_DirectionTuningFit.mat']), 'g_dsi', 'g_osi', 'ang_dir', 'ang_ori', 'dir_yfit_all', 'b_hat_all', 'k1_hat_all', 'R1_hat_all','R2_hat_all', 'u1_hat_all','u2_hat_all', 'sse_all', 'R_square_all')

    
%% Direction tuning curves plotted
respgrat = avg_resp_dir(:,:,1,1,1);
[prefresp, prefdir] = max(respgrat,[],2);
vecshift = zeros(size(respgrat,1),1)+6; 
vecshift = vecshift - prefdir;
respplaid = avg_resp_dir(:,:,:,2,1);

avg_resp_grat =[];
avg_resp_plaid = [];

for iCell = 1:size(respgrat,1)
    avg_resp_grat(iCell,:) = circshift(respgrat(iCell,:),vecshift(iCell),2);
    for i = 1:nMaskPhas
        avg_resp_plaid(iCell,:,i) = circshift(respplaid(iCell,:,i),vecshift(iCell)+(60./int),2);
    end
end


figure;
start = 1;
n = 1;
x=[-150:30:180];
for iCell = 1:nCells
    subplot(5,4,start)
        for im = 1:nMaskPhas
            plot(x, avg_resp_plaid(iCell,:,im))
            hold on
        end
        plot(x, avg_resp_grat(iCell,:),'k') 
        if iCell ==1; legend('0 deg','90 deg','180 deg', '270 deg', 'gratings'); end;
        xlabel('plaid direction')
        ylabel('df/f')
        idx = iCell==ind;
        if any(idx)
            subtitle(['cell ' num2str(iCell)],'fontweight','bold')
        else
            subtitle(['cell ' num2str(iCell)]); 
        end
    start = start+1;
    if start >20
        sgtitle([mouse ' ' date ' Direction tuning curves aligned to preferred grating direction'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_DirectionTuningGratingsAndPlaids_' num2str(n) '.pdf']),'-dpdf', '-fillpage')       
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle([mouse ' ' date ' Direction tuning curves aligned to preferred grating direction'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_DirectionTuningGratingsAndPlaids_' num2str(n) '.pdf']), '-dpdf','-fillpage')
    end
end
close all

%% Polar plots by individual cell

figure;
start = 1;
n = 1;

x=[-150:30:180];
x_rad = deg2rad(x);
for iCell =1:nCells
    subplot(5,4,start)
        for im = 1:nMaskPhas
            polarplot([x_rad x_rad(1)], [avg_resp_plaid(iCell,:,im) avg_resp_plaid(iCell,1,im)])
            hold on
        end
        polarplot([x_rad x_rad(1)], [avg_resp_grat(iCell,:) avg_resp_grat(iCell,1)],'k', 'LineWidth',2) 
        idx = iCell==ind;
        if any(idx)
            subtitle(['cell ' num2str(iCell)],'fontweight','bold')
        else
            subtitle(['cell ' num2str(iCell)])
        end
    start = start+1;    
    if start>20
        sgtitle([mouse ' ' date])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_PolarPlots_' num2str(n) '.pdf']), '-dpdf','-fillpage')
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle([mouse ' ' date])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_PolarPlots_' num2str(n) '.pdf']), '-dpdf','-fillpage')
    end
end     
close all

%% ZcZp by individual cell
figure;
start = 1;
n = 1;

for iCell = 1:nCells
    subplot(5,4,start)
        for im = 1:4
            scatter(Zc(im,iCell), Zp(im,iCell))
            hold on
        end
        ylabel('Zp'); ylim([-4 8]);
        xlabel('Zc'); xlim([-4 8]);
        if iCell ==1; legend('0 deg','90 deg','180 deg', '270 deg'); end;
        idx = iCell==ind;
        if any(idx)
            subtitle(['cell ' num2str(iCell)],'fontweight','bold')
        else
            subtitle(['cell ' num2str(iCell)])
        end
        plotZcZpBorders; set(gca,'TickDir','out'); axis square
    start = start+1;    
    if start>20
        sgtitle([mouse ' ' date ' - Zp Zc by cell'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_ZcZpByCell_' num2str(n) '.pdf']), '-dpdf','-fillpage')
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle([mouse ' ' date ' - Zp Zc by cell'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_ZcZpByCell_' num2str(n) '.pdf']), '-dpdf','-fillpage')
    end        
end
close all


stop
%% look at df/f across individual trials with comparison to statistical tests (ttest for grating, anova for plaid)

% c = linspace(0,round(max(centroid_dist)),round(max(centroid_dist)));


figure;
start=1;
n=1;
for iCell = 1:nCells
    subplot(5,4,start)
        for iDir = 1:nStimDir
%             cent = round(centroid_dist(trialInd{iDir,1,1}));
            scatter(iDir*ones(1,size(resp_cell{iDir,1,1},2)),resp_cell{iDir,1,1}(iCell,:));
            hold on
            max_df = max(resp_cell{iDir,1,1}(iCell,:),[],2); %find max df/f to use to print t-test result
            tt = h_resp(iCell,iDir,1,1);
            if tt == 1; text(iDir, max_df,'*'); end;
        end
        idx = iCell==ind;
        if any(idx)
            subtitle(['cell ' num2str(iCell) ', DSI ' num2str(DSI(iCell),'%.2f') ', p=' num2str(p_anova_dir(iCell),'%.3f')],'fontweight','bold')
        else
            subtitle(['cell ' num2str(iCell) ', DSI ' num2str(DSI(iCell),'%.2f') ', p=' num2str(p_anova_dir(iCell),'%.3f')])
        end
        ylabel('df/f'); xlabel('direction');
    start = start+1;    
    if start>20
        sgtitle([mouse ' ' date ' - Df/f response to gratings, t-test result (0/1), and anova'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_GratingResponses_' num2str(n) '.pdf']), '-dpdf','-fillpage')
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle([mouse ' ' date ' - Df/f response to gratings, t-test result (0/1), and anova'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_GratingResponses_' num2str(n) '.pdf']), '-dpdf','-fillpage')
    end     
end
    close all
   
   
figure;
start=1;
n=1;
for iCell = 1:nCells
    subplot(5,4,start)
        [p_min, p_ind] = min(p_anova_plaid(:,iCell),[],1);
        for iDir = 1:nStimDir
            scatter(iDir*ones(1,size(resp_cell{iDir,p_ind,2},2)),resp_cell{iDir,p_ind,2}(iCell,:));
            hold on
        end
        idx = iCell==ind;
        if any(idx)
            subtitle(['cell ' num2str(iCell) ' - phase ' num2str(p_ind) ', p=' num2str(p_min,'%.3f')],'fontweight','bold')
        else
            subtitle(['cell ' num2str(iCell) ' - phase ' num2str(p_ind) ', p=' num2str(p_min,'%.3f')])
        end
        ylabel('df/f'); xlabel('direction');
    start = start+1;    
    if start>20
        sgtitle([mouse ' ' date ' - Df/f response to plaids with the lowest p (anova)'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_PlaidResponses_' num2str(n) '.pdf']), '-dpdf','-fillpage')
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle([mouse ' ' date ' - Df/f response to plaids with the lowest p (anova)'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_PlaidResponses_' num2str(n) '.pdf']), '-dpdf','-fillpage')
    end        
end
  close all
 end
  %%
  
%   stop
% close all
% %   ind = intersect(intersect(resp_ind_dir,p_dir),OSI_ind);
%   ind = intersect(resp_ind_dir,p_dir);
%   figure;
%   subplot(3,2,1)
%     scatter(DSI,amp_hat_all);
%     xlabel('DSI')
%     ylabel('mod amp')
%   subplot(3,2,2)
%     cdfplot(amp_hat_all(find(DSI>0.5)))
%     hold on
%     cdfplot(amp_hat_all(find(DSI<0.5)))
%     xlabel('mod amp')
%     legend('DS (DSI>0.5)','OS (DSI<0.5)')
%   subplot(3,2,3)
%     cdfplot(mean(PCI(find(DSI>0.5)),1))
%     hold on
%     cdfplot(mean(PCI(find(DSI<0.5)),1))
%     xlabel('mean PCI')
%     legend('DS (DSI>0.5)','OS (DSI<0.5)')
%   subplot(3,2,4)
%     cdfplot(b_hat_all(find(DSI>0.5)))
%     hold on
%     cdfplot(b_hat_all(find(DSI<0.5)))
%     xlabel('mod baseline')
%     legend('DS (DSI>0.5)','OS (DSI<0.5)')  
%   subplot(3,2,5)
%     cdfplot(mean(PCI(:,(find(DSI>0.7))),1))
%     hold on
%     cdfplot(mean(PCI(:,(find(DSI<0.3))),1))
%     xlabel('mean PCI (Zp-Zc)')
%     legend('DS (DSI>0.7)','OS (DSI<0.3)')  
%   subplot(3,2,6)
%     cdfplot(mean(Zp(:,(find(DSI>0.6))),1))
%     hold on
%     cdfplot(mean(Zp(:,(find(DSI<0.4))),1))
%     xlabel('mean Zp')
%     legend('DS (DSI>0.6)','OS (DSI<0.4)')  