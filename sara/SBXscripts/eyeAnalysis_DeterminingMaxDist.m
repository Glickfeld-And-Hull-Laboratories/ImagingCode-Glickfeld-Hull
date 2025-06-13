clc; clear all; close all;
doRedChannel = 1;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';
ds = 'CrossOriRandDirFourPhase_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = length(expt);
frame_rate = 15;
seed = rng;

num=1;
figure;

for iexp = [9 10 13 26 40 51 53 63 64]

    clearvars -except doRedChannel base outDir svName ds rc nexp frame_rate seed num iexp
    eval(ds)

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc;
ImgFolder = expt(iexp).coFolder;
time = expt(iexp).coTime;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

fprintf([mouse ' ' date '\n'])

%% Pref direction analysis
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
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
 
for max_dist = [3 5 15]
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
    
    
    %% Set responsive index
    
    % ind = 1:nCells;
    ind = intersect(intersect(resp_ind_dir,DSI_ind),p_dir); %resp to 1 grating, DSI>0.5, resp to one direction of gratings
    % ind = ind1;
    
    %% PCI modulation
    
    PCI = (Zp-Zc);
    
    %fit sinusoid
    phase = [0 90 180 270];
    phase_range = 0:1:359;
    
    for iCell = 1:nCells
        [b_hat_all(iCell,1), amp_hat_all(iCell,1), per_hat_all(iCell,1),pha_hat_all(iCell,1),sse_all(iCell,1),R_square_all(iCell,1)] = sinefit_PCI(deg2rad(phase),PCI(:,iCell));
        yfit_all(iCell,:,1) = b_hat_all(iCell,1)+amp_hat_all(iCell,1).*(sin(2*pi*deg2rad(phase_range)./per_hat_all(iCell,1) + 2.*pi/pha_hat_all(iCell,1)));
    end

    subplot(3,4,num)
        cdfplot(amp_hat_all)
        hold on
        xlabel('fit amplitude')
        set(gca,'TickDir','out'); box off; axis square; grid off
        subtitle([date ' ' mouse])
        if num==1
            legend('3','5','15')
        end

end
num=num+1;

end

sgtitle('Fit amplitude distributions with different pupil centroid max distances')
print(fullfile(outDir, [svName '_eyeAnalysis_determiningMaxDist.pdf']),'-dpdf', '-fillpage') 

