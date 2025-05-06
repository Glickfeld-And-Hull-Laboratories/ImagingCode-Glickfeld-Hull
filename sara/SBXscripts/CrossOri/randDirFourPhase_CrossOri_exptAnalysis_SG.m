
clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandDirFourPhase_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = length(expt);
iexp = 113; 
frame_rate = 15;

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc;
ImgFolder = expt(iexp).coFolder;
time = expt(iexp).coTime;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

dataBase = (['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\' expt(iexp).saveLoc]);
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
%base = '\\CRASH.dhe.duke.edu\data\home\lindsey';

fprintf([mouse ' ' date '\n'])

%% Pref direction analysis
load(fullfile(dataBase, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
load(fullfile(dataBase, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
load(fullfile(dataBase, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))

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
trialInd = cell(nStimDir,nMaskPhas,2);
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
    ind_stimdir = ind_stimdir(ind_stimdir<2219);
    ind_maskdir = ind_maskdir(ind_maskdir<2219);
    ind_diralone = [intersect(ind_stimdir, ind_stimAlone), intersect(ind_maskdir, ind_maskAlone)];
    ind_dirplaid = [intersect(ind_stimdir, ind_plaid)];
    trialsperstim(iDir,1,1) = length(ind_diralone);
    resp_cell{iDir,1,1} = squeeze(mean(data_dfof_tc(resp_win,:,ind_diralone),1));
    trialInd{iDir,1,1} = ind_diralone;
    [h_resp(:,iDir,1,1), p_resp(:,iDir,1,1)] = ttest2(resp_cell{iDir,1,1},resp_blank,'dim',2,'tail','right','alpha', 0.05./nStim);
    avg_resp_dir(:,iDir,1,1,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind_diralone),1),3));
    avg_resp_dir(:,iDir,1,1,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,ind_diralone),1),[],3)./sqrt(length(ind_diralone)));
    all_resp_dir = [all_resp_dir squeeze(mean(data_dfof_tc(resp_win,:,ind_diralone),1))];
    all_dir = [all_dir iDir.*ones(size(ind_diralone))];
    for ip = 1:nMaskPhas
        ind_dpplaid = intersect(ind_dirplaid,ind_p{ip});
        trialsperstim(iDir,ip,2) = length(ind_dpplaid);
        resp_cell{iDir,ip,2} = squeeze(mean(data_dfof_tc(resp_win,:,ind_dpplaid),1));
        trialInd{iDir,ip,2} = ind_dpplaid;
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
if nMaskPhas > 1
    p_plaid2 = find(p_anova_plaid(2,:)<0.05);
    p_all = unique([p_dir,p_plaid1,p_plaid2]);
    if nMaskPhas > 2
        p_plaid3 = find(p_anova_plaid(3,:)<0.05);
        p_all = unique([p_dir,p_plaid1,p_plaid2,p_plaid3]);
        if nMaskPhas > 3
            p_plaid4 = find(p_anova_plaid(4,:)<0.05);
            p_all = unique([p_dir,p_plaid1,p_plaid2,p_plaid3,p_plaid4]); %significantly responsive to a direction (anova) for gratings or any plaid set
        end
    end 
end

% avg_resp_dir_rect(find(avg_resp_dir(:,:,1,1,1)<0)) = 0;
avg_resp_dir_rect = avg_resp_dir;
for iCell = 1:nCells
    [max_val max_ind] = max(avg_resp_dir_rect(iCell,:,1,1,1));
    null_ind = max_ind+(nStimDir./2);
    null_ind(find(null_ind>nStimDir)) = null_ind(find(null_ind>nStimDir))-nStimDir;
    min_val = avg_resp_dir_rect(iCell,null_ind,1,1,1);
    if min_val < 0; min_val = 0; end;
    DSI(iCell) = (max_val-min_val)./(max_val+min_val);
end

DSI_ind = find(DSI>0.5); %direction selective to gratings

if ~exist(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]), 'dir')
    mkdir(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]));
end
save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']), 'resp_cell', 'data_dfof_tc', 'resp_ind', 'frame_rate', 'h_resp', 'avg_resp_dir','p_anova_dir','p_anova_plaid');
save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']), 'prewin_frames', 'postwin_frames', 'resp_win', 'ind_stimAlone', 'ind_maskAlone', 'ind_plaid', 'ind_blank');
save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_varsforGLM.mat']), 'resp_cell','stimDir_all','trialInd','ind_diralone','ind_dirplaid','ind_dpplaid');

if nMaskPhas == 1
    stop
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
    
% Do all fits at once
    [DSIstruct, ZpZcStruct, plaid_corr, gratingFitStruct, ZpZcPWdist, phaseModStruct] = bigFits(avg_resp_dir);

    % Get direction selectivity
        DSI         = DSIstruct.DSI;
        DSI_ind     = DSIstruct.DS_ind;
        DSI_maxInd  = DSIstruct.prefDir;
    % Get direction tuning curve fit
        dir_b_hat_all           = gratingFitStruct.b;
        k1_hat_all          = gratingFitStruct.k1;
        R1_hat_all          = gratingFitStruct.R1;
        R2_hat_all          = gratingFitStruct.R2;
        u1_hat_all          = gratingFitStruct.u1;
        u2_hat_all          = gratingFitStruct.u2;
        dir_sse_all         = gratingFitStruct.sse;
        dir_R_square_all    = gratingFitStruct.Rsq;
        dir_yfit_all        = gratingFitStruct.yfit;
    % Get partial correlations
        Zp      = ZpZcStruct.Zp;
        Zc      = ZpZcStruct.Zc;
        Rp      = ZpZcStruct.Rp;
        Rc      = ZpZcStruct.Rc;
        ind     = ZpZcStruct.PDSind_byphase;
        ind_pds = ZpZcStruct.PDSind_all;
    % Get PCI fit, get amplitude and baseline
        PCI             = phaseModStruct.PCI;
        yfit_all        = phaseModStruct.yfit;
        amp_hat_all     = phaseModStruct.amp;
        b_hat_all       = phaseModStruct.b;
        sse_all         = phaseModStruct.sse;
        R_square_all    = phaseModStruct.rsq;



%% ZcZp of population (resp_ind)

figure; 
movegui('center')
for i = 1:4
    subplot(4,4,i)
    scatter(Zc(i,resp_ind), Zp(i,resp_ind))
    hold on
    scatter(Zc(i,ind{1}),Zp(i,ind{1}));
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
    scatter(Zc(i,resp_ind), Zp(i,resp_ind))
    hold on
    scatter(Zc(i,ind{2}),Zp(i,ind{2}));
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
    scatter(Zc(i,resp_ind), Zp(i,resp_ind))
    hold on
    scatter(Zc(i,ind{3}),Zp(i,ind{3}));
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
    scatter(Zc(i,resp_ind), Zp(i,resp_ind))
    hold on
    scatter(Zc(i,ind{4}),Zp(i,ind{4}));
    xlabel('Zc')
    ylabel('Zp')
    ylim([-4 8])
    xlim([-4 8])
    hold on
    if i==1; title('pattern cells at 270'); end
    plotZcZpBorders
end
sgtitle('Pattern direction selective cells at four phases')
print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_ZcZp.pdf']),'-dpdf', '-fillpage')       

%% Set responsive index

ind = intersect(intersect(resp_ind_dir,DSI_ind),p_dir); %resp to 1 grating, DSI>0.5, resp to one direction of gratings
% ind = ind1;

%% PCI modulation


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
        idx = iCell==ind; %does iCell equal any value in the index?
        if any(idx)  %if yes, plot in red. if not, plot in black
            plot(phase_range, yfit_all(iCell,:,1),'k'); 
            subtitle(['cell ' num2str(iCell) ', Rsq ' num2str(R_square_all(iCell),'%.3f'), ', SSE ' num2str(sse_all(iCell),'%.2f')],'fontweight','bold')
        else
            plot(phase_range, yfit_all(iCell,:,1),'k:');
            subtitle(['cell ' num2str(iCell) ', Rsq ' num2str(R_square_all(iCell),'%.2f'), ', SSE ' num2str(sse_all(iCell),'%.2f')])
        end
        ylabel('Zp-Zc'); xlabel('Mask phase'); ylim([-6 6])
    start = start+1;
    if start >20
        sgtitle([mouse ' ' date ' PCI modulation across mask phase by cell'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_PCImodulation_' num2str(n) '.pdf']),'-dpdf', '-fillpage')       
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle([mouse ' ' date ' PCI modulation across mask phase by cell'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_PCImodulation_' num2str(n) '.pdf']), '-dpdf','-fillpage')
    end        
end
    save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_PatternCorrelationIndexFits.mat']), 'ind', 'Zp', 'Zc', 'PCI', 'yfit_all', 'b_hat_all', 'amp_hat_all', 'sse_all', 'R_square_all')
    close all



%% direction tuning curves

[avg_resp_grat, avg_resp_plaid] = getAlignedGratPlaidTuning(avg_resp_dir);

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
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_DirectionTuningGratingsAndPlaids_' num2str(n) '.pdf']),'-dpdf', '-fillpage')       
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle([mouse ' ' date ' Direction tuning curves aligned to preferred grating direction'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_DirectionTuningGratingsAndPlaids_' num2str(n) '.pdf']), '-dpdf','-fillpage')
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
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_PolarPlots_' num2str(n) '.pdf']), '-dpdf','-fillpage')
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle([mouse ' ' date])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_PolarPlots_' num2str(n) '.pdf']), '-dpdf','-fillpage')
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
        plotZcZpBorders
    start = start+1;    
    if start>20
        sgtitle([mouse ' ' date ' - Zp Zc by cell'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_ZcZpByCell_' num2str(n) '.pdf']), '-dpdf','-fillpage')
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle([mouse ' ' date ' - Zp Zc by cell'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_ZcZpByCell_' num2str(n) '.pdf']), '-dpdf','-fillpage')
    end        
end
close all

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
    
print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TrialsPerStim.pdf']),'-dpdf', '-fillpage')       

    stop
%% look at df/f across individual trials with comparison to statistical tests (ttest for grating, anova for plaid)


figure;
start=1;
n=1;
for iCell = 1:nCells
    subplot(5,4,start)
        for iDir = 1:nStimDir
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
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_GratingResponses_' num2str(n) '.pdf']), '-dpdf','-fillpage')
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle([mouse ' ' date ' - Df/f response to gratings, t-test result (0/1), and anova'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_GratingResponses_' num2str(n) '.pdf']), '-dpdf','-fillpage')
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
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_PlaidResponses_' num2str(n) '.pdf']), '-dpdf','-fillpage')
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle([mouse ' ' date ' - Df/f response to plaids with the lowest p (anova)'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_PlaidResponses_' num2str(n) '.pdf']), '-dpdf','-fillpage')
    end        
end
  close all
   
