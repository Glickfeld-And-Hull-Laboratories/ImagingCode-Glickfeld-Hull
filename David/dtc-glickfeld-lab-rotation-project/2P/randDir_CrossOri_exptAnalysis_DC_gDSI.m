clc; clear all; close all; clear all global;
startup
doRedChannel = 0;
ds = 'CrossOriRandDir_ExptList_DC';
%rc = behavConstsAV;
eval(ds)
nexp = length(expt);
iexp = 7;
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
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\David';

% Add Sara's function library to path
addpath(genpath('Z:\All_Staff\home\David\repositories\dtc-glickfeld-lab-rotation-project\CrossOri_PatternV1_SG_DC'));

fprintf([mouse ' ' date '\n'])

%% Pref direction analysis
fprintf('Current stage: Loading data\n')
load(fullfile(dataBase, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
load(fullfile(dataBase, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
load(fullfile(dataBase, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))

%% Derive plaid offset from data
signed_offset = unique(maskDiffs(maskDiffs ~= 0));
if isempty(signed_offset)
    signed_offset = 45;
end
plaid_offset = abs(signed_offset(1));  % magnitude (for alignment functions)
offset_rad = deg2rad(signed_offset(1));  % signed (for VA/IOC)
fprintf('Plaid offset: %+.0f deg (mask - stim direction)\n', signed_offset(1))

% Component speeds for VA/IOC prediction lines
SF_val = expt(iexp).SF;
if isfield(expt, 'TF_stim')
    speed_stim = expt(iexp).TF_stim / SF_val;
    speed_mask = expt(iexp).TF_mask / SF_val;
else
    speed_stim = expt(iexp).TF / SF_val;
    speed_mask = speed_stim;
end
fprintf('VA/IOC params: speed_stim=%.1f, speed_mask=%.1f deg/s\n', speed_stim, speed_mask)

% VS/IOC angular shifts (radians) — stimulus property, computed once
alpha_VS = atan2(speed_mask * sin(offset_rad), speed_stim + speed_mask * cos(offset_rad));
alpha_IOC = atan2(speed_mask - speed_stim * cos(offset_rad), speed_stim * sin(offset_rad));
shift_VS_deg = -rad2deg(alpha_VS);
shift_IOC_deg = -rad2deg(alpha_IOC);

% C2 (mask-tracking) prediction
% If a cell tracks C2 (mask grating), delta_pref = +signed_offset(1)
c2_angle_deltapref = signed_offset(1);     % C2 prediction on delta_pref polar plot
shift_C2_deg = -signed_offset(1);          % for scatter plots: y = x + shift_C2
fprintf('C2 (mask-tracking) prediction: delta_pref = %+.1f deg, scatter shift = %+.1f deg\n', ...
    c2_angle_deltapref, shift_C2_deg)

%%
fprintf('Current stage: Response extraction\n')
if doRedChannel == 0
    red_cells = [];
end

prewin_frames = unique(celleqel2mat_padded(input.tItiWaitFrames))./3;
nFramesOn = unique(celleqel2mat_padded(input.nStimOneFramesOn));
postwin_frames = unique(celleqel2mat_padded(input.nStimOneFramesOn));
tt = (1-prewin_frames:postwin_frames).*(1/frame_rate);
data_resp = nan(prewin_frames+postwin_frames,nCells,nTrials);
data_f = nan(1,nCells,nTrials);

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
resp_blank = squeeze(mean(data_dfof_tc(prewin_frames+5:prewin_frames+nFramesOn,:,ind_blank),1));
avg_resp_dir = zeros(nCells, nStimDir,nMaskPhas, 2, 2);
all_resp_dir = [];
all_resp_plaid = cell(1,nMaskPhas);
all_dir = [];
all_plaid = cell(1,nMaskPhas);
nStim = nStimDir;
ind_diralone_all = cell(1,nStimDir);
ind_dirplaid_all = cell(1,nStimDir);
ind_dpplaid_all  = cell(nStimDir,nMaskPhas);
for iDir = 1:nStimDir
    ind_stimdir = find(stimDir_all == stimDirs(iDir));
    ind_maskdir = find(maskDir_all == stimDirs(iDir));
    ind_diralone = [intersect(ind_stimdir, ind_stimAlone), intersect(ind_maskdir, ind_maskAlone)];
    ind_dirplaid = [intersect(ind_stimdir, ind_plaid)];
    ind_diralone_all{iDir} = ind_diralone;
    ind_dirplaid_all{iDir} = ind_dirplaid;
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
        ind_dpplaid_all{iDir,ip} = ind_dpplaid;
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
p_anova_plaid = zeros(nMaskPhas,nCells);
for iCell = 1:nCells
    p_anova_dir(iCell) = anova1(all_resp_dir(iCell,:), all_dir, 'off'); %direction selective to gratings
    for ip = 1:nMaskPhas
        p_anova_plaid(ip,iCell) = anova1(all_resp_plaid{ip}(iCell,:), all_plaid{ip}, 'off'); %direction selective to plaids
    end
end

p_dir = find(p_anova_dir<0.05);
p_plaid1 = find(p_anova_plaid(1,:)<0.05);
p_all = unique([p_dir, p_plaid1]);
if nMaskPhas > 1
    p_plaid2 = find(p_anova_plaid(2,:)<0.05);
    p_all = unique([p_all, p_plaid2]);
end
if nMaskPhas > 2
    p_plaid3 = find(p_anova_plaid(3,:)<0.05);
    p_all = unique([p_all, p_plaid3]);
end
if nMaskPhas > 3
    p_plaid4 = find(p_anova_plaid(4,:)<0.05);
    p_all = unique([p_all, p_plaid4]); %significantly responsive to a direction (anova) for gratings or any plaid set
end

% --- Old peak DSI (commented out — replaced by gDSI) ---
% avg_resp_dir_rect = avg_resp_dir;
% for iCell = 1:nCells
%     [max_val max_ind] = max(avg_resp_dir_rect(iCell,:,1,1,1));
%     null_ind = max_ind+(nStimDir./2);
%     null_ind(find(null_ind>nStimDir)) = null_ind(find(null_ind>nStimDir))-nStimDir;
%     min_val = avg_resp_dir_rect(iCell,null_ind,1,1,1);
%     if min_val < 0; min_val = 0; end;
%     DSI(iCell) = (max_val-min_val)./(max_val+min_val);
% end
% DSI_ind = find(DSI>0.5);
%
% DSI_plaid = zeros(1,nCells);
% for iCell = 1:nCells
%     [max_val_plaid, max_ind_plaid] = max(avg_resp_dir_rect(iCell,:,1,2,1));
%     null_ind_plaid = max_ind_plaid + (nStimDir./2);
%     null_ind_plaid(find(null_ind_plaid > nStimDir)) = null_ind_plaid(find(null_ind_plaid > nStimDir)) - nStimDir;
%     min_val_plaid = avg_resp_dir_rect(iCell,null_ind_plaid,1,2,1);
%     if min_val_plaid < 0; min_val_plaid = 0; end;
%     DSI_plaid(iCell) = (max_val_plaid - min_val_plaid)./(max_val_plaid + min_val_plaid);
% end
% DSI_plaid_ind = find(DSI_plaid > 0.5);

outDir = fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], 'gDSI');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
save(fullfile(outDir,[date '_' mouse '_' run_str '_respData.mat']), 'resp_cell', 'data_dfof_tc', 'resp_ind', 'frame_rate', 'h_resp', 'avg_resp_dir','p_anova_dir','p_anova_plaid');
save(fullfile(outDir,[date '_' mouse '_' run_str '_stimData.mat']), 'prewin_frames', 'postwin_frames', 'resp_win', 'ind_stimAlone', 'ind_maskAlone', 'ind_plaid', 'ind_blank');
save(fullfile(outDir,[date '_' mouse '_' run_str '_varsforGLM.mat']), 'resp_cell','stimDir_all','trialInd','ind_diralone_all','ind_dirplaid_all','ind_dpplaid_all');

%%
fprintf('Current stage: Fitting (DSI, tuning curves)\n')

% --- Split-half random subsampling (unused, commented out) ---
% avg_resp_dir_rand = zeros(nCells,nStimDir,2);
% for i = 1:nStimDir
%     n = size(resp_cell{i,1,2},2);
%     ind1 = randsample(1:n,ceil(n/2));
%     ind2 = setdiff(1:n,ind1);
%     avg_resp_dir_rand(:,i,1) = mean(resp_cell{i,1,2}(:,ind1),2);
%     avg_resp_dir_rand(:,i,2) = mean(resp_cell{i,1,2}(:,ind2),2);
% end
    
% Core fits — work with any number of phases
% DSIstruct        = getDSIstruct(avg_resp_dir);  % replaced by getGDSI below
% ZpZcStruct       = getZpZcStruct(avg_resp_dir, 'alignedTestDir', plaid_offset);

% VS/IOC-resolved partial correlations
ZpZcVSIOC = getZpZcStruct_VSIOC(avg_resp_dir, alpha_VS, alpha_IOC, 'whole_cell', plaid_offset);
fprintf('  Zp/Zc VSIOC: degenerate=%d, VS-PDS=%d, IOC-PDS=%d, CDS=%d\n', ...
    ZpZcVSIOC.is_degenerate, length(ZpZcVSIOC.ind_VS_PDS), ...
    length(ZpZcVSIOC.ind_IOC_PDS), length(ZpZcVSIOC.ind_CDS))

% gDSI — Li et al. (2025, Curr Biol)
gDSI_grat_struct = getGDSI(avg_resp_dir(:,:,1,1,1), stimDirs);
gDSI_plaid_struct = getGDSI(avg_resp_dir(:,:,1,2,1), stimDirs);

gDSI = gDSI_grat_struct.gDSI;
gDSI_prefDir = gDSI_grat_struct.prefDir_deg;
gDSI_ind = gDSI_grat_struct.DS_ind;

gDSI_plaid = gDSI_plaid_struct.gDSI;
gDSI_plaid_prefDir = gDSI_plaid_struct.prefDir_deg;
gDSI_plaid_ind = gDSI_plaid_struct.DS_ind;

% --- Peak DSI (for comparison with gDSI) ---
DSI = zeros(1, nCells);
DSI_maxInd = zeros(1, nCells);
for iCell = 1:nCells
    [max_val, max_ind] = max(avg_resp_dir(iCell,:,1,1,1));
    null_ind = max_ind + (nStimDir/2);
    if null_ind > nStimDir; null_ind = null_ind - nStimDir; end
    min_val = avg_resp_dir(iCell, null_ind, 1, 1, 1);
    if min_val < 0; min_val = 0; end
    DSI(iCell) = (max_val - min_val) / (max_val + min_val);
    DSI_maxInd(iCell) = max_ind;
end

gratingFitStruct = getGratingTuningCurveFit(avg_resp_dir);

% Multi-phase fits — require 4 phases
if nMaskPhas >= 4
    plaid_corr     = getPlaidTuningCorrelations(avg_resp_dir);
    % ZpZcPWdist     = getZpZcPWdist(ZpZcStruct);
    % phaseModStruct = get4PhaseModulationFit(ZpZcStruct);
else
    plaid_corr     = [];
    % ZpZcPWdist     = [];
    % phaseModStruct = [];
end

    % Get direction selectivity (old peak DSI from getDSIstruct — replaced by gDSI)
        % DSI         = DSIstruct.DSI;
        % DSI_ind     = DSIstruct.DS_ind;
        % DSI_maxInd  = DSIstruct.prefDir;

% --- Weighted circular mean DSI (replaced by gDSI above) ---
% NOTE: The circular mean (angle(sum(R_theta * exp(i*theta)))) can flip 180 deg
% for cells with weak/noisy tuning where negative dF/F values dominate the
% vector sum. This causes negative DSI_w values and 180-deg-off VS/IOC
% predictions on the weighted polar plots. Li et al. (2025, Curr Biol) avoid
% this by using gDSI = |sum(R*exp(i*theta))| / sum(R) with a 0.2 threshold,
% which naturally filters out cells where the circular mean is unreliable.
% DSI_w = zeros(1, nCells);
% DSI_w_plaid = zeros(1, nCells);
% dirs_deg = stimDirs(:)';
% dirs_rad = deg2rad(dirs_deg);
% dirs_wrap = [dirs_deg, dirs_deg(1) + 360];
%
% for iCell = 1:nCells
%     % Grating
%     tc_g = avg_resp_dir(iCell,:,1,1,1);
%     tc_g_wrap = [tc_g, tc_g(1)];
%     pref_deg_g = mod(rad2deg(angle(sum(tc_g .* exp(1i * dirs_rad)))), 360);
%     null_deg_g = mod(pref_deg_g + 180, 360);
%     R_pref_g = interp1(dirs_wrap, tc_g_wrap, pref_deg_g, 'linear');
%     R_null_g = interp1(dirs_wrap, tc_g_wrap, null_deg_g, 'linear');
%     if R_null_g < 0; R_null_g = 0; end
%     DSI_w(iCell) = (R_pref_g - R_null_g) / (R_pref_g + R_null_g);
%
%     % Plaid (mask phase 1)
%     tc_p = avg_resp_dir(iCell,:,1,2,1);
%     tc_p_wrap = [tc_p, tc_p(1)];
%     pref_deg_p = mod(rad2deg(angle(sum(tc_p .* exp(1i * dirs_rad)))), 360);
%     null_deg_p = mod(pref_deg_p + 180, 360);
%     R_pref_p = interp1(dirs_wrap, tc_p_wrap, pref_deg_p, 'linear');
%     R_null_p = interp1(dirs_wrap, tc_p_wrap, null_deg_p, 'linear');
%     if R_null_p < 0; R_null_p = 0; end
%     DSI_w_plaid(iCell) = (R_pref_p - R_null_p) / (R_pref_p + R_null_p);
% end
%
% DSI_w_ind = find(DSI_w > 0.5);
% DSI_w_plaid_ind = find(DSI_w_plaid > 0.5);

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
    %     Zp      = ZpZcStruct.Zp;
    %     Zc      = ZpZcStruct.Zc;
    %     Rp      = ZpZcStruct.Rp;
    %     Rc      = ZpZcStruct.Rc;
    %     ind     = ZpZcStruct.PDSind_byphase;
    %     ind_pds = ZpZcStruct.PDSind_all;
    % Get PCI fit, get amplitude and baseline (requires 4 phases)
    % if nMaskPhas >= 4
    %     PCI             = phaseModStruct.PCI;
    %     yfit_all        = phaseModStruct.yfit;
    %     amp_hat_all     = phaseModStruct.amp;
    %     b_hat_all       = phaseModStruct.b;
    %     sse_all         = phaseModStruct.sse;
    %     R_square_all    = phaseModStruct.rsq;
    % end



% %% ZcZp of population (resp_ind)
%
% figure('Visible','off', 'Position', [0 0 1920 1080]);
% movegui('center')
% for ip = 1:nMaskPhas
%     subplot(1, nMaskPhas, ip)
%     scatter(Zc(ip, resp_ind), Zp(ip, resp_ind))
%     hold on
%     scatter(Zc(ip, ind{ip}), Zp(ip, ind{ip}));
%     xlabel('Zc'); ylabel('Zp')
%     ylim([-4 8]); xlim([-4 8])
%     title(['PDS cells at phase ' num2str(maskPhas(ip))])
%     plotZcZpBorders
% end
% sgtitle('Pattern direction selective cells by phase')
% print(fullfile(outDir,[date '_' mouse '_' run_str '_ZcZp.png']),'-dpng', '-r300')

%% Set responsive index

ind = intersect(intersect(resp_ind_dir, gDSI_ind), p_dir); %resp to 1 grating, gDSI>0.2, resp to one direction of gratings
ind_dsi_both = find(gDSI > 0.2 & gDSI_plaid > 0.2);
% ind_w = intersect(intersect(resp_ind_dir, DSI_w_ind), p_dir);  % replaced by gDSI
% ind_dsi_both_w = find(DSI_w > 0.5 & DSI_w_plaid > 0.5);  % replaced by gDSI
% ind = ind1;

%% Cell filtering summary
fprintf('Current stage: Plotting cell filtering summary\n')

funnel_labels = {'Total imaged', 'Resp (any)', 'Resp (grating)', ...
    'ANOVA dir-sel', 'Resp + ANOVA', 'Resp + ANOVA + gDSI', 'gDSI both'};
funnel_counts = [nCells, length(resp_ind), length(resp_ind_dir), ...
    length(p_dir), length(intersect(resp_ind_dir, p_dir)), ...
    length(ind), length(ind_dsi_both)];

figure('Visible','off', 'Position', [0 0 900 500]);
bh = barh(flip(funnel_counts), 'FaceColor', 'flat');
colors = [linspace(0.85, 0.1, 7); linspace(0.85, 0.3, 7); linspace(0.85, 0.7, 7)]';
bh.CData = flip(colors);
yticks(1:7)
yticklabels(flip(funnel_labels))
xlabel('Number of cells')
for i = 1:7
    text(funnel_counts(8-i) + max(funnel_counts)*0.02, i, ...
        sprintf('%d (%.0f%%)', funnel_counts(8-i), 100*funnel_counts(8-i)/nCells), ...
        'VerticalAlignment', 'middle', 'FontSize', 9)
end
title([mouse ' ' date ' — Cell Filtering Summary'], 'FontSize', 11)
xlim([0 max(funnel_counts)*1.25])

print(fullfile(outDir, ...
    [date '_' mouse '_' run_str '_CellFilteringSummary_gDSI.png']), '-dpng', '-r300')

% %% PCI modulation
%
% if nMaskPhas >= 4
%     figure('Visible','off', 'Position', [0 0 1920 1080]);
%     for i = 1:nMaskPhas
%         cdfplot(PCI(i,:));
%         hold on
%     end
%
%     %fit sinusoid
%     phase = [0 90 180 270];
%     phase_range = 0:1:359;
%     figure('Visible','off', 'Position', [0 0 1920 1080]);
%     start=1;
%     n=1;
%     for iCell = 1:nCells
%         subplot(5,4,start)
%             scatter(phase,PCI(:,iCell),'LineWidth',1.25);
%             hold on
%             idx = iCell==ind; %does iCell equal any value in the index?
%             if any(idx)  %if yes, plot in red. if not, plot in black
%                 plot(phase_range, yfit_all(iCell,:,1),'k');
%                 subtitle(['cell ' num2str(iCell) ', Rsq ' num2str(R_square_all(iCell),'%.3f'), ', SSE ' num2str(sse_all(iCell),'%.2f')],'fontweight','bold')
%             else
%                 plot(phase_range, yfit_all(iCell,:,1),'k:');
%                 subtitle(['cell ' num2str(iCell) ', Rsq ' num2str(R_square_all(iCell),'%.2f'), ', SSE ' num2str(sse_all(iCell),'%.2f')])
%             end
%             ylabel('Zp-Zc'); xlabel('Mask phase'); ylim([-6 6])
%         start = start+1;
%         if start >20
%             sgtitle([mouse ' ' date ' PCI modulation across mask phase by cell'])
%             print(fullfile(outDir,[date '_' mouse '_' run_str '_PCImodulation_' num2str(n) '.png']),'-dpng', '-r300')
%             figure('Visible','off', 'Position', [0 0 1920 1080]);
%     
%             start = 1;
%             n = n+1;
%         end
%         if iCell == nCells
%             sgtitle([mouse ' ' date ' PCI modulation across mask phase by cell'])
%             print(fullfile(outDir,[date '_' mouse '_' run_str '_PCImodulation_' num2str(n) '.png']), '-dpng', '-r300')
%         end
%     end
%     save(fullfile(outDir,[date '_' mouse '_' run_str '_PatternCorrelationIndexFits.mat']), 'ind', 'Zp', 'Zc', 'PCI', 'yfit_all', 'b_hat_all', 'amp_hat_all', 'sse_all', 'R_square_all')
%     %close all
% else
%     % Save what we have for 1-phase
%     save(fullfile(outDir,[date '_' mouse '_' run_str '_PatternCorrelationIndexFits.mat']), 'ind', 'Zp', 'Zc')
% end



%% Unaligned direction tuning curves (absolute directions, with von Mises fit)
fprintf('Current stage: Plotting unaligned tuning curves\n')

% --- Old unaligned tuning curves (peak DSI, commented out) ---
% figure('Visible','off', 'Position', [0 0 1920 1080]);
% start = 1;
% n = 1;
% for iCell = 1:nCells
%     subplot(5,4,start)
%         plot(stimDirs, avg_resp_dir(iCell,:,1,1,1), 'ko-', 'MarkerSize', 3)
%         hold on
%         plot(stimDirs, avg_resp_dir(iCell,:,1,2,1), 'bo-', 'MarkerSize', 3)
%         pref_dir = stimDirs(DSI_maxInd(iCell));
%         vs_pred = mod(pref_dir - rad2deg(alpha_VS), 360);
%         ioc_pred = mod(pref_dir - rad2deg(alpha_IOC), 360);
%         xline(vs_pred, 'g--', 'LineWidth', 1);
%         xline(ioc_pred, 'm--', 'LineWidth', 1);
%         if iCell == 1; legend({'grating','plaid','VS pred','IOC pred'},'FontSize',6); end
%         xlabel('direction (deg)')
%         ylabel('df/f')
%         idx = iCell==ind;
%         if any(idx)
%             subtitle(['cell ' num2str(iCell)],'fontweight','bold')
%         else
%             subtitle(['cell ' num2str(iCell)])
%         end
%     start = start+1;
%     if start > 20
%         sgtitle([mouse ' ' date ' Unaligned tuning curves: grating (black) + plaid (blue) vs absolute direction'], 'FontSize', 10)
%         print(fullfile(outDir,[date '_' mouse '_' run_str '_UnalignedTuning_' num2str(n) '.png']),'-dpng', '-r300')
%         figure('Visible','off', 'Position', [0 0 1920 1080]);
%         start = 1;
%         n = n+1;
%     end
%     if iCell == nCells
%         sgtitle([mouse ' ' date ' Unaligned tuning curves: grating (black) + plaid (blue) vs absolute direction'], 'FontSize', 10)
%         print(fullfile(outDir,[date '_' mouse '_' run_str '_UnalignedTuning_' num2str(n) '.png']),'-dpng', '-r300')
%     end
% end

% --- New unaligned tuning curves (gDSI pref dir) ---
figure('Visible','off', 'Position', [0 0 1920 1080]);
start = 1;
n = 1;
for iCell = 1:nCells
    subplot(5,4,start)
        plot(stimDirs, avg_resp_dir(iCell,:,1,1,1), 'ko-', 'MarkerSize', 3)
        hold on
        plot(stimDirs, avg_resp_dir(iCell,:,1,2,1), 'bo-', 'MarkerSize', 3)
        pref_dir = gDSI_prefDir(iCell);
        vs_pred = mod(pref_dir - rad2deg(alpha_VS), 360);
        ioc_pred = mod(pref_dir - rad2deg(alpha_IOC), 360);
        xline(vs_pred, 'g--', 'LineWidth', 1);
        xline(ioc_pred, 'm--', 'LineWidth', 1);
        if iCell == 1; legend({'grating','plaid','VS pred','IOC pred'},'FontSize',6); end
        xlabel('direction (deg)')
        ylabel('df/f')
        idx = iCell==ind;
        if any(idx)
            subtitle(['cell ' num2str(iCell)],'fontweight','bold')
        else
            subtitle(['cell ' num2str(iCell)])
        end
    start = start+1;
    if start > 20
        sgtitle([mouse ' ' date ' Unaligned tuning curves (gDSI pref dir): grating (black) + plaid (blue)'], 'FontSize', 10)
        print(fullfile(outDir,[date '_' mouse '_' run_str '_UnalignedTuning_gDSI_' num2str(n) '.png']),'-dpng', '-r300')
        figure('Visible','off', 'Position', [0 0 1920 1080]);
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle([mouse ' ' date ' Unaligned tuning curves (gDSI pref dir): grating (black) + plaid (blue)'], 'FontSize', 10)
        print(fullfile(outDir,[date '_' mouse '_' run_str '_UnalignedTuning_gDSI_' num2str(n) '.png']),'-dpng', '-r300')
    end
end
%close all

%% Polar plots (absolute direction, peak pref dir)
fprintf('Current stage: Plotting polar plots (absolute direction)\n')

% --- Old polar plots (all cells, peak DSI argmax, commented out) ---
% figure('Visible','off', 'Position', [0 0 1920 1080]);
% start = 1;
% n = 1;
% x_rad = deg2rad(stimDirs);
% for iCell = 1:nCells
%     subplot(5,4,start)
%         % Plaid tuning (colored, one per mask phase)
%         for im = 1:nMaskPhas
%             h = polarplot([x_rad x_rad(1)], [avg_resp_dir(iCell,:,im,2,1) avg_resp_dir(iCell,1,im,2,1)]);
%             if im == 1; h_plaid = h; end
%             hold on
%         end
%         % Grating tuning (black, thick)
%         h_grat = polarplot([x_rad x_rad(1)], [avg_resp_dir(iCell,:,1,1,1) avg_resp_dir(iCell,1,1,1,1)], 'k', 'LineWidth', 2);
%         % Peak preferred direction lines (argmax bin)
%         rmax = max(rlim);
%         theta_pref = deg2rad(stimDirs(DSI_maxInd(iCell)));
%         polarplot([theta_pref theta_pref], [0 rmax], 'k-', 'LineWidth', 1.5)
%         colors = get(gca, 'ColorOrder');
%         for im = 1:nMaskPhas
%             [~, plaid_max_ind] = max(avg_resp_dir(iCell,:,im,2,1));
%             theta_plaid_peak = deg2rad(stimDirs(plaid_max_ind));
%             polarplot([theta_plaid_peak theta_plaid_peak], [0 rmax], '-', 'Color', colors(im,:), 'LineWidth', 1.5)
%         end
%         % VS and IOC predicted plaid preferred direction
%         h_vs = polarplot([theta_pref - alpha_VS, theta_pref - alpha_VS], [0 rmax], 'g--', 'LineWidth', 1.5);
%         h_ioc = polarplot([theta_pref - alpha_IOC, theta_pref - alpha_IOC], [0 rmax], 'm--', 'LineWidth', 1.5);
%         if iCell == 1; legend([h_plaid, h_grat, h_vs, h_ioc], {'plaid','grating','VS','IOC'}, 'FontSize', 5, 'Location', 'best'); end
%         % Subtitle with DSI info
%         idx = iCell==ind;
%         if any(idx)
%             subtitle(['cell ' num2str(iCell) ', DSI ' num2str(DSI(iCell),'%.2f') ', pDSI ' num2str(DSI_plaid(iCell),'%.2f')], 'fontweight', 'bold')
%         else
%             subtitle(['cell ' num2str(iCell) ', DSI ' num2str(DSI(iCell),'%.2f') ', pDSI ' num2str(DSI_plaid(iCell),'%.2f')])
%         end
%     start = start+1;
%     if start > 20
%         sgtitle([mouse ' ' date ' Polar plots (peak pref dir): grating (black) + plaid (colored) + VS (green) + IOC (magenta)'], 'FontSize', 10)
%         print(fullfile(outDir,[date '_' mouse '_' run_str '_PolarPlots_peak_' num2str(n) '.png']), '-dpng', '-r300')
%         figure('Visible','off', 'Position', [0 0 1920 1080]);
%         start = 1;
%         n = n+1;
%     end
%     if iCell == nCells
%         sgtitle([mouse ' ' date ' Polar plots (peak pref dir): grating (black) + plaid (colored) + VS (green) + IOC (magenta)'], 'FontSize', 10)
%         print(fullfile(outDir,[date '_' mouse '_' run_str '_PolarPlots_peak_' num2str(n) '.png']), '-dpng', '-r300')
%     end
% end

% --- New polar plots (all cells, gDSI pref dir) ---
figure('Visible','off', 'Position', [0 0 1920 1080]);
start = 1;
n = 1;
x_rad = deg2rad(stimDirs);
for iCell = 1:nCells
    subplot(5,4,start)
        % Plaid tuning (colored, one per mask phase)
        for im = 1:nMaskPhas
            h = polarplot([x_rad x_rad(1)], [avg_resp_dir(iCell,:,im,2,1) avg_resp_dir(iCell,1,im,2,1)]);
            if im == 1; h_plaid = h; end
            hold on
        end
        % Grating tuning (black, thick)
        h_grat = polarplot([x_rad x_rad(1)], [avg_resp_dir(iCell,:,1,1,1) avg_resp_dir(iCell,1,1,1,1)], 'k', 'LineWidth', 2);
        % gDSI preferred direction lines (circular mean)
        rmax = max(rlim);
        theta_pref = deg2rad(gDSI_prefDir(iCell));
        polarplot([theta_pref theta_pref], [0 rmax], 'k-', 'LineWidth', 1.5)
        colors = get(gca, 'ColorOrder');
        for im = 1:nMaskPhas
            % Use plaid gDSI pref dir for phase 1, argmax for others
            if im == 1
                theta_plaid_pref = deg2rad(gDSI_plaid_prefDir(iCell));
            else
                [~, plaid_max_ind] = max(avg_resp_dir(iCell,:,im,2,1));
                theta_plaid_pref = deg2rad(stimDirs(plaid_max_ind));
            end
            polarplot([theta_plaid_pref theta_plaid_pref], [0 rmax], '-', 'Color', colors(im,:), 'LineWidth', 1.5)
        end
        % VS and IOC predicted plaid preferred direction
        h_vs = polarplot([theta_pref - alpha_VS, theta_pref - alpha_VS], [0 rmax], 'g--', 'LineWidth', 1.5);
        h_ioc = polarplot([theta_pref - alpha_IOC, theta_pref - alpha_IOC], [0 rmax], 'm--', 'LineWidth', 1.5);
        if iCell == 1; legend([h_plaid, h_grat, h_vs, h_ioc], {'plaid','grating','VS','IOC'}, 'FontSize', 5, 'Location', 'best'); end
        % Subtitle with gDSI info
        idx = iCell==ind;
        if any(idx)
            subtitle(['cell ' num2str(iCell) ', gDSI ' num2str(gDSI(iCell),'%.2f') ', pDSI ' num2str(gDSI_plaid(iCell),'%.2f')], 'fontweight', 'bold')
        else
            subtitle(['cell ' num2str(iCell) ', gDSI ' num2str(gDSI(iCell),'%.2f') ', pDSI ' num2str(gDSI_plaid(iCell),'%.2f')])
        end
    start = start+1;
    if start > 20
        sgtitle([mouse ' ' date ' Polar plots (gDSI pref dir): grating (black) + plaid (colored) + VS (green) + IOC (magenta)'], 'FontSize', 10)
        print(fullfile(outDir,[date '_' mouse '_' run_str '_PolarPlots_gDSI_' num2str(n) '.png']), '-dpng', '-r300')
        figure('Visible','off', 'Position', [0 0 1920 1080]);
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle([mouse ' ' date ' Polar plots (gDSI pref dir): grating (black) + plaid (colored) + VS (green) + IOC (magenta)'], 'FontSize', 10)
        print(fullfile(outDir,[date '_' mouse '_' run_str '_PolarPlots_gDSI_' num2str(n) '.png']), '-dpng', '-r300')
    end
end
%close all

% --- Old polar plots (all cells, weighted DSI, commented out — gDSI polar plots above cover this) ---
% figure('Visible','off', 'Position', [0 0 1920 1080]);
% start = 1;
% n = 1;
% for iCell = 1:nCells
%     subplot(5,4,start)
%         % Plaid tuning (colored, one per mask phase)
%         for im = 1:nMaskPhas
%             h = polarplot([x_rad x_rad(1)], [avg_resp_dir(iCell,:,im,2,1) avg_resp_dir(iCell,1,im,2,1)]);
%             if im == 1; h_plaid = h; end
%             hold on
%         end
%         % Grating tuning (black, thick)
%         h_grat = polarplot([x_rad x_rad(1)], [avg_resp_dir(iCell,:,1,1,1) avg_resp_dir(iCell,1,1,1,1)], 'k', 'LineWidth', 2);
%         % Weighted circular mean preferred direction lines
%         rmax = max(rlim);
%         vec_grat = sum(avg_resp_dir(iCell,:,1,1,1) .* exp(1i * x_rad));
%         polarplot([angle(vec_grat) angle(vec_grat)], [0 rmax], 'k-', 'LineWidth', 1.5)
%         colors = get(gca, 'ColorOrder');
%         for im = 1:nMaskPhas
%             vec_plaid = sum(avg_resp_dir(iCell,:,im,2,1) .* exp(1i * x_rad));
%             polarplot([angle(vec_plaid) angle(vec_plaid)], [0 rmax], '-', 'Color', colors(im,:), 'LineWidth', 1.5)
%         end
%         % VS and IOC predicted plaid preferred direction
%         theta_pref = angle(vec_grat);
%         h_vs = polarplot([theta_pref - alpha_VS, theta_pref - alpha_VS], [0 rmax], 'g--', 'LineWidth', 1.5);
%         h_ioc = polarplot([theta_pref - alpha_IOC, theta_pref - alpha_IOC], [0 rmax], 'm--', 'LineWidth', 1.5);
%         if iCell == 1; legend([h_plaid, h_grat, h_vs, h_ioc], {'plaid','grating','VS','IOC'}, 'FontSize', 5, 'Location', 'best'); end
%         % Subtitle with DSI info
%         idx = iCell==ind;
%         if any(idx)
%             subtitle(['cell ' num2str(iCell) ', DSI_w ' num2str(DSI_w(iCell),'%.2f') ', pDSI_w ' num2str(DSI_w_plaid(iCell),'%.2f')], 'fontweight', 'bold')
%         else
%             subtitle(['cell ' num2str(iCell) ', DSI_w ' num2str(DSI_w(iCell),'%.2f') ', pDSI_w ' num2str(DSI_w_plaid(iCell),'%.2f')])
%         end
%     start = start+1;
%     if start > 20
%         sgtitle([mouse ' ' date ' Polar plots (weighted pref dir): grating (black) + plaid (colored) + VS (green) + IOC (magenta)'], 'FontSize', 10)
%         print(fullfile(outDir,[date '_' mouse '_' run_str '_PolarPlots_weighted_' num2str(n) '.png']), '-dpng', '-r300')
%         figure('Visible','off', 'Position', [0 0 1920 1080]);
%         start = 1;
%         n = n+1;
%     end
%     if iCell == nCells
%         sgtitle([mouse ' ' date ' Polar plots (weighted pref dir): grating (black) + plaid (colored) + VS (green) + IOC (magenta)'], 'FontSize', 10)
%         print(fullfile(outDir,[date '_' mouse '_' run_str '_PolarPlots_weighted_' num2str(n) '.png']), '-dpng', '-r300')
%     end
% end

%% Polar plots — direction-selective subset
fprintf('Current stage: Plotting polar plots (direction-selective subset)\n')

% --- Old polar plots (DS subset, peak DSI argmax, commented out) ---
% figure('Visible','off', 'Position', [0 0 1920 1080]);
% start = 1;
% n = 1;
% for ii = 1:length(ind_dsi_both)
%     iCell = ind_dsi_both(ii);
%     subplot(5,4,start)
%         % Plaid tuning (colored, one per mask phase)
%         for im = 1:nMaskPhas
%             h = polarplot([x_rad x_rad(1)], [avg_resp_dir(iCell,:,im,2,1) avg_resp_dir(iCell,1,im,2,1)]);
%             if im == 1; h_plaid = h; end
%             hold on
%         end
%         % Grating tuning (black, thick)
%         h_grat = polarplot([x_rad x_rad(1)], [avg_resp_dir(iCell,:,1,1,1) avg_resp_dir(iCell,1,1,1,1)], 'k', 'LineWidth', 2);
%         % Peak preferred direction lines (argmax bin)
%         rmax = max(rlim);
%         theta_pref = deg2rad(stimDirs(DSI_maxInd(iCell)));
%         polarplot([theta_pref theta_pref], [0 rmax], 'k-', 'LineWidth', 1.5)
%         colors = get(gca, 'ColorOrder');
%         for im = 1:nMaskPhas
%             [~, plaid_max_ind] = max(avg_resp_dir(iCell,:,im,2,1));
%             theta_plaid_peak = deg2rad(stimDirs(plaid_max_ind));
%             polarplot([theta_plaid_peak theta_plaid_peak], [0 rmax], '-', 'Color', colors(im,:), 'LineWidth', 1.5)
%         end
%         % VS and IOC predicted plaid preferred direction
%         h_vs = polarplot([theta_pref - alpha_VS, theta_pref - alpha_VS], [0 rmax], 'g--', 'LineWidth', 1.5);
%         h_ioc = polarplot([theta_pref - alpha_IOC, theta_pref - alpha_IOC], [0 rmax], 'm--', 'LineWidth', 1.5);
%         if ii == 1; legend([h_plaid, h_grat, h_vs, h_ioc], {'plaid','grating','VS','IOC'}, 'FontSize', 5, 'Location', 'best'); end
%         subtitle(['cell ' num2str(iCell) ', DSI ' num2str(DSI(iCell),'%.2f') ', pDSI ' num2str(DSI_plaid(iCell),'%.2f')])
%     start = start+1;
%     if start > 20
%         sgtitle([mouse ' ' date ' Polar plots (DS cells, peak DSI): grating (black) + plaid (colored) + VS (green) + IOC (magenta)'], 'FontSize', 10)
%         print(fullfile(outDir,[date '_' mouse '_' run_str '_PolarPlots_DS_peak_' num2str(n) '.png']), '-dpng', '-r300')
%         figure('Visible','off', 'Position', [0 0 1920 1080]);
%         start = 1;
%         n = n+1;
%     end
%     if ii == length(ind_dsi_both)
%         sgtitle([mouse ' ' date ' Polar plots (DS cells, peak DSI): grating (black) + plaid (colored) + VS (green) + IOC (magenta)'], 'FontSize', 10)
%         print(fullfile(outDir,[date '_' mouse '_' run_str '_PolarPlots_DS_peak_' num2str(n) '.png']), '-dpng', '-r300')
%     end
% end

% --- New polar plots (DS subset, gDSI > 0.2) ---
figure('Visible','off', 'Position', [0 0 1920 1080]);
start = 1;
n = 1;
for ii = 1:length(ind_dsi_both)
    iCell = ind_dsi_both(ii);
    subplot(5,4,start)
        % Plaid tuning (colored, one per mask phase)
        for im = 1:nMaskPhas
            h = polarplot([x_rad x_rad(1)], [avg_resp_dir(iCell,:,im,2,1) avg_resp_dir(iCell,1,im,2,1)]);
            if im == 1; h_plaid = h; end
            hold on
        end
        % Grating tuning (black, thick)
        h_grat = polarplot([x_rad x_rad(1)], [avg_resp_dir(iCell,:,1,1,1) avg_resp_dir(iCell,1,1,1,1)], 'k', 'LineWidth', 2);
        % gDSI preferred direction lines (circular mean)
        rmax = max(rlim);
        theta_pref = deg2rad(gDSI_prefDir(iCell));
        polarplot([theta_pref theta_pref], [0 rmax], 'k-', 'LineWidth', 1.5)
        colors = get(gca, 'ColorOrder');
        for im = 1:nMaskPhas
            if im == 1
                theta_plaid_pref = deg2rad(gDSI_plaid_prefDir(iCell));
            else
                [~, plaid_max_ind] = max(avg_resp_dir(iCell,:,im,2,1));
                theta_plaid_pref = deg2rad(stimDirs(plaid_max_ind));
            end
            polarplot([theta_plaid_pref theta_plaid_pref], [0 rmax], '-', 'Color', colors(im,:), 'LineWidth', 1.5)
        end
        % VS and IOC predicted plaid preferred direction
        h_vs = polarplot([theta_pref - alpha_VS, theta_pref - alpha_VS], [0 rmax], 'g--', 'LineWidth', 1.5);
        h_ioc = polarplot([theta_pref - alpha_IOC, theta_pref - alpha_IOC], [0 rmax], 'm--', 'LineWidth', 1.5);
        if ii == 1; legend([h_plaid, h_grat, h_vs, h_ioc], {'plaid','grating','VS','IOC'}, 'FontSize', 5, 'Location', 'best'); end
        subtitle(['cell ' num2str(iCell) ', gDSI ' num2str(gDSI(iCell),'%.2f') ', pDSI ' num2str(gDSI_plaid(iCell),'%.2f')])
    start = start+1;
    if start > 20
        sgtitle([mouse ' ' date ' Polar plots (DS cells, gDSI > 0.2): grating (black) + plaid (colored) + VS (green) + IOC (magenta)'], 'FontSize', 10)
        print(fullfile(outDir,[date '_' mouse '_' run_str '_PolarPlots_DS_gDSI_' num2str(n) '.png']), '-dpng', '-r300')
        figure('Visible','off', 'Position', [0 0 1920 1080]);
        start = 1;
        n = n+1;
    end
    if ii == length(ind_dsi_both)
        sgtitle([mouse ' ' date ' Polar plots (DS cells, gDSI > 0.2): grating (black) + plaid (colored) + VS (green) + IOC (magenta)'], 'FontSize', 10)
        print(fullfile(outDir,[date '_' mouse '_' run_str '_PolarPlots_DS_gDSI_' num2str(n) '.png']), '-dpng', '-r300')
    end
end
%close all

% --- Old polar plots (DS subset, weighted DSI, commented out — gDSI DS polar plots above cover this) ---
% figure('Visible','off', 'Position', [0 0 1920 1080]);
% start = 1;
% n = 1;
% for ii = 1:length(ind_dsi_both_w)
%     iCell = ind_dsi_both_w(ii);
%     subplot(5,4,start)
%         % Plaid tuning (colored, one per mask phase)
%         for im = 1:nMaskPhas
%             h = polarplot([x_rad x_rad(1)], [avg_resp_dir(iCell,:,im,2,1) avg_resp_dir(iCell,1,im,2,1)]);
%             if im == 1; h_plaid = h; end
%             hold on
%         end
%         % Grating tuning (black, thick)
%         h_grat = polarplot([x_rad x_rad(1)], [avg_resp_dir(iCell,:,1,1,1) avg_resp_dir(iCell,1,1,1,1)], 'k', 'LineWidth', 2);
%         % Weighted circular mean preferred direction lines
%         rmax = max(rlim);
%         vec_grat = sum(avg_resp_dir(iCell,:,1,1,1) .* exp(1i * x_rad));
%         polarplot([angle(vec_grat) angle(vec_grat)], [0 rmax], 'k-', 'LineWidth', 1.5)
%         colors = get(gca, 'ColorOrder');
%         for im = 1:nMaskPhas
%             vec_plaid = sum(avg_resp_dir(iCell,:,im,2,1) .* exp(1i * x_rad));
%             polarplot([angle(vec_plaid) angle(vec_plaid)], [0 rmax], '-', 'Color', colors(im,:), 'LineWidth', 1.5)
%         end
%         % VS and IOC predicted plaid preferred direction
%         theta_pref = angle(vec_grat);
%         h_vs = polarplot([theta_pref - alpha_VS, theta_pref - alpha_VS], [0 rmax], 'g--', 'LineWidth', 1.5);
%         h_ioc = polarplot([theta_pref - alpha_IOC, theta_pref - alpha_IOC], [0 rmax], 'm--', 'LineWidth', 1.5);
%         if ii == 1; legend([h_plaid, h_grat, h_vs, h_ioc], {'plaid','grating','VS','IOC'}, 'FontSize', 5, 'Location', 'best'); end
%         subtitle(['cell ' num2str(iCell) ', DSI_w ' num2str(DSI_w(iCell),'%.2f') ', pDSI_w ' num2str(DSI_w_plaid(iCell),'%.2f')])
%     start = start+1;
%     if start > 20
%         sgtitle([mouse ' ' date ' Polar plots (DS cells, weighted DSI): grating (black) + plaid (colored) + VS (green) + IOC (magenta)'], 'FontSize', 10)
%         print(fullfile(outDir,[date '_' mouse '_' run_str '_PolarPlots_DS_weighted_' num2str(n) '.png']), '-dpng', '-r300')
%         figure('Visible','off', 'Position', [0 0 1920 1080]);
%         start = 1;
%         n = n+1;
%     end
%     if ii == length(ind_dsi_both_w)
%         sgtitle([mouse ' ' date ' Polar plots (DS cells, weighted DSI): grating (black) + plaid (colored) + VS (green) + IOC (magenta)'], 'FontSize', 10)
%         print(fullfile(outDir,[date '_' mouse '_' run_str '_PolarPlots_DS_weighted_' num2str(n) '.png']), '-dpng', '-r300')
%     end
% end

%% Grating vs Plaid preferred direction scatter
fprintf('Current stage: Plotting grating vs plaid preferred direction scatter\n')

% --- Old weighted circular-mean computation (commented out) ---
% theta_grat_all = zeros(1, nCells);
% theta_plaid_all = zeros(1, nCells);
% for iCell = 1:nCells
%     vec_g = sum(avg_resp_dir(iCell,:,1,1,1) .* exp(1i * x_rad));
%     vec_p = sum(avg_resp_dir(iCell,:,1,2,1) .* exp(1i * x_rad));
%     theta_grat_all(iCell) = rad2deg(angle(vec_g));
%     theta_plaid_all(iCell) = rad2deg(angle(vec_p));
% end
% theta_grat_all = mod(theta_grat_all, 360);
% theta_plaid_all = mod(theta_plaid_all, 360);

% --- Old peak preferred directions (argmax bin, commented out) ---
% theta_grat_all_peak = zeros(1, nCells);
% theta_plaid_all_peak = zeros(1, nCells);
% for iCell = 1:nCells
%     theta_grat_all_peak(iCell) = stimDirs(DSI_maxInd(iCell));
%     [~, plaid_max_ind] = max(avg_resp_dir(iCell,:,1,2,1));
%     theta_plaid_all_peak(iCell) = stimDirs(plaid_max_ind);
% end

% --- Old peak scatter plot (commented out) ---
% figure('Visible','off', 'Position', [0 0 1080 1080]);
% hold on
% plot([0 360], [0 360], 'k--', 'LineWidth', 1)
% plot([0 360], [0+shift_VS_deg 360+shift_VS_deg], 'g--', 'LineWidth', 1)
% plot([0 360], [0+shift_IOC_deg 360+shift_IOC_deg], 'm--', 'LineWidth', 1)
% scatter(theta_grat_all_peak, theta_plaid_all_peak, 30, 'k', 'LineWidth', 0.8)
% xlabel('Grating preferred direction (deg)')
% ylabel('Plaid preferred direction (deg)')
% xlim([0 360]); ylim([0 360])
% axis square
% legend({'component (y=x)', ['VS (y=x' sprintf('%+.1f', shift_VS_deg) ')'], ...
%         ['IOC (y=x' sprintf('%+.1f', shift_IOC_deg) ')'], ...
%         'all cells'}, 'Location', 'best', 'FontSize', 7)
% sgtitle([mouse ' ' date ' Grating vs Plaid preferred direction (peak pref dir)'], 'FontSize', 10)
% print(fullfile(outDir,...
%     [date '_' mouse '_' run_str '_GratVsPlaidPrefDir_peak.png']), '-dpng', '-r300')

% --- Old weighted scatter plot (commented out — gDSI scatter below covers this) ---
% figure('Visible','off', 'Position', [0 0 1080 1080]);
% hold on
% plot([0 360], [0 360], 'k--', 'LineWidth', 1)
% plot([0 360], [0+shift_VS_deg 360+shift_VS_deg], 'g--', 'LineWidth', 1)
% plot([0 360], [0+shift_IOC_deg 360+shift_IOC_deg], 'm--', 'LineWidth', 1)
% scatter(theta_grat_all, theta_plaid_all, 30, 'k', 'LineWidth', 0.8)
% xlabel('Grating preferred direction (deg)')
% ylabel('Plaid preferred direction (deg)')
% xlim([0 360]); ylim([0 360])
% axis square
% legend({'component (y=x)', ['VS (y=x' sprintf('%+.1f', shift_VS_deg) ')'], ...
%         ['IOC (y=x' sprintf('%+.1f', shift_IOC_deg) ')'], ...
%         'all cells'}, 'Location', 'best', 'FontSize', 7)
% sgtitle([mouse ' ' date ' Grating vs Plaid preferred direction (weighted pref dir)'], 'FontSize', 10)
% print(fullfile(outDir,...
%     [date '_' mouse '_' run_str '_GratVsPlaidPrefDir_weighted.png']), '-dpng', '-r300')

% --- New gDSI scatter plot ---
theta_grat_gDSI = gDSI_prefDir;          % already 1 x nCells, 0-360
theta_plaid_gDSI = gDSI_plaid_prefDir;   % already 1 x nCells, 0-360

figure('Visible','off', 'Position', [0 0 1080 1080]);
hold on
plot([0 360], [0+shift_VS_deg 360+shift_VS_deg], 'g--', 'LineWidth', 1)  % VS
plot([0 360], [0+shift_VS_deg-360 360+shift_VS_deg-360], 'g--', 'LineWidth', 1, 'HandleVisibility', 'off')  % VS wrap
plot([0 360], [0+shift_VS_deg+360 360+shift_VS_deg+360], 'g--', 'LineWidth', 1, 'HandleVisibility', 'off')  % VS wrap
plot([0 360], [0+shift_IOC_deg 360+shift_IOC_deg], 'm--', 'LineWidth', 1)  % IOC
plot([0 360], [0+shift_IOC_deg-360 360+shift_IOC_deg-360], 'm--', 'LineWidth', 1, 'HandleVisibility', 'off')  % IOC wrap
plot([0 360], [0+shift_IOC_deg+360 360+shift_IOC_deg+360], 'm--', 'LineWidth', 1, 'HandleVisibility', 'off')  % IOC wrap
ind_nonDS = setdiff(1:nCells, ind_dsi_both);
scatter(theta_grat_gDSI(ind_nonDS), theta_plaid_gDSI(ind_nonDS), 30, 'k', 'LineWidth', 0.8)
scatter(theta_grat_gDSI(ind_dsi_both), theta_plaid_gDSI(ind_dsi_both), 30, 'b', 'filled')

xlabel('Grating preferred direction (deg)')
ylabel('Plaid preferred direction (deg)')
xlim([0 360]); ylim([0 360])
axis square
legend({['VS (y=x' sprintf('%+.1f', shift_VS_deg) ')'], ...
        ['IOC (y=x' sprintf('%+.1f', shift_IOC_deg) ')'], ...
        ['non-DS (n=' num2str(length(ind_nonDS)) ')'], ...
        ['DS gDSI>0.2 (n=' num2str(length(ind_dsi_both)) ')']}, ...
       'Location', 'best', 'FontSize', 7)
sgtitle([mouse ' ' date ' Grat vs Plaid pref dir (gDSI) — DS: ' ...
    num2str(length(ind_dsi_both)) '/' num2str(nCells)], 'FontSize', 10)

print(fullfile(outDir,...
    [date '_' mouse '_' run_str '_GratVsPlaidPrefDir_gDSI.png']), '-dpng', '-r300')

% --- gDSI scatter plot (DS cells only, gDSI > 0.2 both) ---
figure('Visible','off', 'Position', [0 0 1080 1080]);
hold on
plot([0 360], [0+shift_VS_deg 360+shift_VS_deg], 'g--', 'LineWidth', 1)
plot([0 360], [0+shift_VS_deg-360 360+shift_VS_deg-360], 'g--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot([0 360], [0+shift_VS_deg+360 360+shift_VS_deg+360], 'g--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot([0 360], [0+shift_IOC_deg 360+shift_IOC_deg], 'm--', 'LineWidth', 1)
plot([0 360], [0+shift_IOC_deg-360 360+shift_IOC_deg-360], 'm--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot([0 360], [0+shift_IOC_deg+360 360+shift_IOC_deg+360], 'm--', 'LineWidth', 1, 'HandleVisibility', 'off')

scatter(theta_grat_gDSI(ind_dsi_both), theta_plaid_gDSI(ind_dsi_both), 30, 'k', 'LineWidth', 0.8)

xlabel('Grating preferred direction (deg)')
ylabel('Plaid preferred direction (deg)')
xlim([0 360]); ylim([0 360])
axis square
legend({['VS (y=x' sprintf('%+.1f', shift_VS_deg) ')'], ...
        ['IOC (y=x' sprintf('%+.1f', shift_IOC_deg) ')'], ...
        ['DS cells (n=' num2str(length(ind_dsi_both)) ')']}, 'Location', 'best', 'FontSize', 7)
sgtitle([mouse ' ' date ' Grating vs Plaid pref dir (gDSI > 0.2 both, n=' num2str(length(ind_dsi_both)) '/' num2str(nCells) ')'], 'FontSize', 10)

print(fullfile(outDir,...
    [date '_' mouse '_' run_str '_GratVsPlaidPrefDir_gDSI_DS.png']), '-dpng', '-r300')

% --- gDSI scatter: colored by gDSI (grating) ---
figure('Visible','off', 'Position', [0 0 1080 1080]);
hold on
plot([0 360], [0+shift_VS_deg 360+shift_VS_deg], 'g--', 'LineWidth', 1)
plot([0 360], [0+shift_VS_deg-360 360+shift_VS_deg-360], 'g--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot([0 360], [0+shift_VS_deg+360 360+shift_VS_deg+360], 'g--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot([0 360], [0+shift_IOC_deg 360+shift_IOC_deg], 'm--', 'LineWidth', 1)
plot([0 360], [0+shift_IOC_deg-360 360+shift_IOC_deg-360], 'm--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot([0 360], [0+shift_IOC_deg+360 360+shift_IOC_deg+360], 'm--', 'LineWidth', 1, 'HandleVisibility', 'off')

scatter(theta_grat_gDSI, theta_plaid_gDSI, 30, gDSI, 'filled')
colormap(parula)
cb = colorbar; cb.Label.String = 'gDSI (grating)';
set(gca, 'CLim', [0 1])

xlabel('Grating preferred direction (deg)')
ylabel('Plaid preferred direction (deg)')
xlim([0 360]); ylim([0 360])
axis square
legend({['VS (y=x' sprintf('%+.1f', shift_VS_deg) ')'], ...
        ['IOC (y=x' sprintf('%+.1f', shift_IOC_deg) ')']}, 'Location', 'best', 'FontSize', 7)
sgtitle({[mouse ' ' date ' Grating vs Plaid pref dir — colored by gDSI (grating)'], ...
    ['n=' num2str(nCells) ' cells']}, 'FontSize', 10)

print(fullfile(outDir,...
    [date '_' mouse '_' run_str '_GratVsPlaidPrefDir_gDSI_colorGrat.png']), '-dpng', '-r300')

% --- gDSI scatter: colored by gDSI (plaid) ---
figure('Visible','off', 'Position', [0 0 1080 1080]);
hold on
plot([0 360], [0+shift_VS_deg 360+shift_VS_deg], 'g--', 'LineWidth', 1)
plot([0 360], [0+shift_VS_deg-360 360+shift_VS_deg-360], 'g--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot([0 360], [0+shift_VS_deg+360 360+shift_VS_deg+360], 'g--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot([0 360], [0+shift_IOC_deg 360+shift_IOC_deg], 'm--', 'LineWidth', 1)
plot([0 360], [0+shift_IOC_deg-360 360+shift_IOC_deg-360], 'm--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot([0 360], [0+shift_IOC_deg+360 360+shift_IOC_deg+360], 'm--', 'LineWidth', 1, 'HandleVisibility', 'off')

scatter(theta_grat_gDSI, theta_plaid_gDSI, 30, gDSI_plaid, 'filled')
colormap(parula)
cb = colorbar; cb.Label.String = 'gDSI (plaid)';
set(gca, 'CLim', [0 1])

xlabel('Grating preferred direction (deg)')
ylabel('Plaid preferred direction (deg)')
xlim([0 360]); ylim([0 360])
axis square
legend({['VS (y=x' sprintf('%+.1f', shift_VS_deg) ')'], ...
        ['IOC (y=x' sprintf('%+.1f', shift_IOC_deg) ')']}, 'Location', 'best', 'FontSize', 7)
sgtitle({[mouse ' ' date ' Grating vs Plaid pref dir — colored by gDSI (plaid)'], ...
    ['n=' num2str(nCells) ' cells']}, 'FontSize', 10)

print(fullfile(outDir,...
    [date '_' mouse '_' run_str '_GratVsPlaidPrefDir_gDSI_colorPlaid.png']), '-dpng', '-r300')

% --- gDSI scatter (DS cells only): colored by gDSI (grating) ---
figure('Visible','off', 'Position', [0 0 1080 1080]);
hold on
plot([0 360], [0+shift_VS_deg 360+shift_VS_deg], 'g--', 'LineWidth', 1)
plot([0 360], [0+shift_VS_deg-360 360+shift_VS_deg-360], 'g--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot([0 360], [0+shift_VS_deg+360 360+shift_VS_deg+360], 'g--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot([0 360], [0+shift_IOC_deg 360+shift_IOC_deg], 'm--', 'LineWidth', 1)
plot([0 360], [0+shift_IOC_deg-360 360+shift_IOC_deg-360], 'm--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot([0 360], [0+shift_IOC_deg+360 360+shift_IOC_deg+360], 'm--', 'LineWidth', 1, 'HandleVisibility', 'off')

scatter(theta_grat_gDSI(ind_dsi_both), theta_plaid_gDSI(ind_dsi_both), 30, gDSI(ind_dsi_both), 'filled')
colormap(parula)
cb = colorbar; cb.Label.String = 'gDSI (grating)';
set(gca, 'CLim', [0 1])

xlabel('Grating preferred direction (deg)')
ylabel('Plaid preferred direction (deg)')
xlim([0 360]); ylim([0 360])
axis square
legend({['VS (y=x' sprintf('%+.1f', shift_VS_deg) ')'], ...
        ['IOC (y=x' sprintf('%+.1f', shift_IOC_deg) ')'], ...
        ['DS cells (n=' num2str(length(ind_dsi_both)) ')']}, 'Location', 'best', 'FontSize', 7)
sgtitle({[mouse ' ' date ' Grating vs Plaid pref dir — colored by gDSI (grating)'], ...
    ['gDSI > 0.2 both, n=' num2str(length(ind_dsi_both)) '/' num2str(nCells)]}, 'FontSize', 10)

print(fullfile(outDir,...
    [date '_' mouse '_' run_str '_GratVsPlaidPrefDir_gDSI_DS_colorGrat.png']), '-dpng', '-r300')

% --- gDSI scatter (DS cells only): colored by gDSI (plaid) ---
figure('Visible','off', 'Position', [0 0 1080 1080]);
hold on
plot([0 360], [0+shift_VS_deg 360+shift_VS_deg], 'g--', 'LineWidth', 1)
plot([0 360], [0+shift_VS_deg-360 360+shift_VS_deg-360], 'g--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot([0 360], [0+shift_VS_deg+360 360+shift_VS_deg+360], 'g--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot([0 360], [0+shift_IOC_deg 360+shift_IOC_deg], 'm--', 'LineWidth', 1)
plot([0 360], [0+shift_IOC_deg-360 360+shift_IOC_deg-360], 'm--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot([0 360], [0+shift_IOC_deg+360 360+shift_IOC_deg+360], 'm--', 'LineWidth', 1, 'HandleVisibility', 'off')

scatter(theta_grat_gDSI(ind_dsi_both), theta_plaid_gDSI(ind_dsi_both), 30, gDSI_plaid(ind_dsi_both), 'filled')
colormap(parula)
cb = colorbar; cb.Label.String = 'gDSI (plaid)';
set(gca, 'CLim', [0 1])

xlabel('Grating preferred direction (deg)')
ylabel('Plaid preferred direction (deg)')
xlim([0 360]); ylim([0 360])
axis square
legend({['VS (y=x' sprintf('%+.1f', shift_VS_deg) ')'], ...
        ['IOC (y=x' sprintf('%+.1f', shift_IOC_deg) ')'], ...
        ['DS cells (n=' num2str(length(ind_dsi_both)) ')']}, 'Location', 'best', 'FontSize', 7)
sgtitle({[mouse ' ' date ' Grating vs Plaid pref dir — colored by gDSI (plaid)'], ...
    ['gDSI > 0.2 both, n=' num2str(length(ind_dsi_both)) '/' num2str(nCells)]}, 'FontSize', 10)

print(fullfile(outDir,...
    [date '_' mouse '_' run_str '_GratVsPlaidPrefDir_gDSI_DS_colorPlaid.png']), '-dpng', '-r300')

%% Delta preferred direction polar histogram
fprintf('Current stage: Plotting delta preferred direction polar histogram\n')

% --- Old delta pref dir polar histogram (peak DSI, commented out) ---
% delta_pref = theta_grat_all_peak(ind_dsi_both) - theta_plaid_all_peak(ind_dsi_both);
% delta_pref = mod(delta_pref + 180, 360) - 180;
% figure('Visible','off', 'Position', [0 0 720 720]);
% polarhistogram(deg2rad(delta_pref), deg2rad(-180:11.25:180), ...
%     'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.6, 'EdgeColor', 'k')
% hold on
% r_max = rlim;
% r_max = r_max(2);
% vs_angle_deg = -shift_VS_deg;
% ioc_angle_deg = -shift_IOC_deg;
% polarplot(deg2rad([vs_angle_deg vs_angle_deg]), [0 r_max], 'b--', 'LineWidth', 1.5)
% polarplot(deg2rad([ioc_angle_deg ioc_angle_deg]), [0 r_max], 'r--', 'LineWidth', 1.5)
% legend({'', 'VS', 'IOC'}, 'Location', 'best', 'FontSize', 8)
% sgtitle({[mouse ' ' date ' \DeltaPrefDir (grating - plaid) (peak DSI)'], ...
%     ['DSI > 0.5 both, n=' num2str(length(ind_dsi_both)) '/' num2str(nCells)]}, 'FontSize', 10)
% print(fullfile(outDir,...
%     [date '_' mouse '_' run_str '_DeltaPrefDirPolar_peak.png']), '-dpng', '-r300')

% --- Old delta pref dir polar histogram (weighted DSI, commented out — gDSI version below covers this) ---
% delta_pref_w = theta_grat_all(ind_dsi_both_w) - theta_plaid_all(ind_dsi_both_w);
% delta_pref_w = mod(delta_pref_w + 180, 360) - 180;
% figure('Visible','off', 'Position', [0 0 720 720]);
% polarhistogram(deg2rad(delta_pref_w), deg2rad(-180:11.25:180), ...
%     'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.6, 'EdgeColor', 'k')
% hold on
% r_max = rlim;
% r_max = r_max(2);
% vs_angle_deg = -shift_VS_deg;
% ioc_angle_deg = -shift_IOC_deg;
% polarplot(deg2rad([vs_angle_deg vs_angle_deg]), [0 r_max], 'b--', 'LineWidth', 1.5)
% polarplot(deg2rad([ioc_angle_deg ioc_angle_deg]), [0 r_max], 'r--', 'LineWidth', 1.5)
% legend({'', 'VS', 'IOC'}, 'Location', 'best', 'FontSize', 8)
% sgtitle({[mouse ' ' date ' \DeltaPrefDir (grating - plaid) (weighted DSI)'], ...
%     ['DSI_w > 0.5 both, n=' num2str(length(ind_dsi_both_w)) '/' num2str(nCells)]}, 'FontSize', 10)
% print(fullfile(outDir,...
%     [date '_' mouse '_' run_str '_DeltaPrefDirPolar_weighted.png']), '-dpng', '-r300')

% --- New delta pref dir polar histogram (gDSI) with 4-category classification ---
delta_pref = theta_grat_gDSI(ind_dsi_both) - theta_plaid_gDSI(ind_dsi_both);
delta_pref = mod(delta_pref + 180, 360) - 180;

vs_angle_deg = -shift_VS_deg;
ioc_angle_deg = -shift_IOC_deg;
half_window_deg = (360 / nStimDir) / 4;  % quarter stimulus step (±5.625° for 16 dirs)

% Classify cells: 4-category (VS / IOC / C1 / C2) with pairwise degeneracy
clf_result = classifyDeltaPrefDir(delta_pref, vs_angle_deg, ioc_angle_deg, c2_angle_deltapref, half_window_deg);

figure('Visible','off', 'Position', [0 0 720 720]);

polarhistogram(deg2rad(delta_pref), deg2rad(-180:11.25:180), ...
    'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.6, 'EdgeColor', 'k')
hold on

r_max = rlim;
r_max = r_max(2);

% Prediction radial lines
h_vs  = polarplot(deg2rad([vs_angle_deg vs_angle_deg]), [0 r_max], 'b--', 'LineWidth', 1.5);
h_ioc = polarplot(deg2rad([ioc_angle_deg ioc_angle_deg]), [0 r_max], 'r--', 'LineWidth', 1.5);
% Colored scatter ring at r_max * 1.08
r_ring = r_max * 1.08;
if ~isempty(clf_result.ind_VS)
    polarplot(deg2rad(delta_pref(clf_result.ind_VS)), r_ring*ones(1,clf_result.n_VS), ...
        'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerEdgeColor', 'none')
end
if ~isempty(clf_result.ind_IOC)
    polarplot(deg2rad(delta_pref(clf_result.ind_IOC)), r_ring*ones(1,clf_result.n_IOC), ...
        'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.8 0.2 0.2], 'MarkerEdgeColor', 'none')
end
if ~isempty(clf_result.ind_C1)
    polarplot(deg2rad(delta_pref(clf_result.ind_C1)), r_ring*ones(1,clf_result.n_C1), ...
        'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerEdgeColor', 'none')
end
if ~isempty(clf_result.ind_pattern)
    polarplot(deg2rad(delta_pref(clf_result.ind_pattern)), r_ring*ones(1,clf_result.n_pattern), ...
        'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.5 0.3 0.7], 'MarkerEdgeColor', 'none')
end
if ~isempty(clf_result.ind_unclass)
    polarplot(deg2rad(delta_pref(clf_result.ind_unclass)), r_ring*ones(1,clf_result.n_unclass), ...
        'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.75 0.75 0.75], 'MarkerEdgeColor', 'none')
end

legend([h_vs, h_ioc], {'VS', 'IOC'}, 'Location', 'best', 'FontSize', 8)

% Build classification string
clf_str = sprintf('VS:%d  IOC:%d  C1:%d', clf_result.n_VS, clf_result.n_IOC, clf_result.n_C1);
if clf_result.n_pattern > 0
    clf_str = [clf_str sprintf('  Pattern(merged):%d', clf_result.n_pattern)];
end
clf_str = [clf_str sprintf('  Unclass:%d', clf_result.n_unclass)];
if ~isempty(clf_result.degenerate_pairs)
    degen_strs = cellfun(@(p) [p{1} '/' p{2}], clf_result.degenerate_pairs, 'UniformOutput', false);
    clf_str = [clf_str '  [near-degen: ' strjoin(degen_strs, ', ') ']'];
end

sgtitle({[mouse ' ' date ' \DeltaPrefDir (grating - plaid) (gDSI)'], ...
    ['gDSI > 0.2 both, n=' num2str(length(ind_dsi_both)) '/' num2str(nCells) ', window=\pm' num2str(half_window_deg) '\circ'], ...
    clf_str}, 'FontSize', 9)

fprintf('  Delta pref dir classification: %s\n', clf_str)

print(fullfile(outDir,...
    [date '_' mouse '_' run_str '_DeltaPrefDirPolar_gDSI.png']), '-dpng', '-r300')

%% Cell Classification Pie Chart & Bar
fprintf('Current stage: Plotting cell classification pie chart\n')

% Counts for pie/bar
pie_labels = {'VS', 'IOC', 'C1', 'Unclassified'};
pie_counts = [clf_result.n_VS, clf_result.n_IOC, clf_result.n_C1, clf_result.n_unclass + clf_result.n_C2];
pie_colors = [0.2 0.4 0.8; 0.8 0.2 0.2; 0.4 0.4 0.4; 0.75 0.75 0.75];

% If pattern (merged VS+IOC) exists, add as separate slice
if clf_result.n_pattern > 0
    pie_labels = [pie_labels, {'VS/IOC (merged)'}];
    pie_counts = [pie_counts, clf_result.n_pattern];
    pie_colors = [pie_colors; 0.5 0.3 0.7];
end

% Remove zero-count categories for cleaner pie
nz = pie_counts > 0;
pie_labels_nz = pie_labels(nz);
pie_counts_nz = pie_counts(nz);
pie_colors_nz = pie_colors(nz, :);
n_total_clf = sum(pie_counts_nz);

figure('Visible', 'off', 'Position', [0 0 1200 500]);

subplot(1, 2, 1)
if ~isempty(pie_counts_nz) && any(pie_counts_nz > 0)
    h_pie = pie(pie_counts_nz);
    % Color the pie slices
    for pp = 1:length(pie_counts_nz)
        h_pie(2*pp-1).FaceColor = pie_colors_nz(pp, :);
    end
    legend(pie_labels_nz, 'Location', 'bestoutside', 'FontSize', 8)
end
title(sprintf('Classification (n=%d, window=\\pm%d\\circ)', n_total_clf, half_window_deg))

subplot(1, 2, 2)
bh = barh(pie_counts_nz, 'FaceColor', 'flat');
bh.CData = pie_colors_nz;
yticks(1:length(pie_counts_nz))
yticklabels(pie_labels_nz)
xlabel('Number of cells')
for pp = 1:length(pie_counts_nz)
    text(pie_counts_nz(pp) + max(pie_counts_nz)*0.03, pp, ...
        sprintf('%d (%.1f%%)', pie_counts_nz(pp), 100*pie_counts_nz(pp)/n_total_clf), ...
        'VerticalAlignment', 'middle', 'FontSize', 9)
end
xlim([0 max(pie_counts_nz)*1.3])
title('Cell counts by category')

sgtitle({[mouse ' ' date ' — Cell Classification Summary (gDSI > 0.2 both)'], ...
    sprintf('n=%d cells, half-window=%d\\circ', n_total_clf, half_window_deg)}, 'FontSize', 10)

print(fullfile(outDir, ...
    [date '_' mouse '_' run_str '_CellClassification_gDSI.png']), '-dpng', '-r300')

%% Circular V-test
fprintf('Current stage: Circular V-test\n')

% V-test: tests whether distribution clusters around a specified direction
% Uses Mardia-Jupp (1999) approximation for one-tailed Rayleigh test
alpha_rad = deg2rad(delta_pref);
n_vtest = length(alpha_rad);

vtest_predictions = struct();
vtest_names = {'VS', 'IOC', 'C1', 'C2'};
vtest_angles = [vs_angle_deg, ioc_angle_deg, 0, c2_angle_deltapref];

for vv = 1:4
    mu_rad = deg2rad(vtest_angles(vv));
    C_bar = mean(cos(alpha_rad - mu_rad));
    V_stat = n_vtest * C_bar;
    u_stat = V_stat * sqrt(2 / n_vtest);
    % Mardia-Jupp approximation for p-value
    p_val = exp(sqrt(1 + 4*n_vtest + 4*(n_vtest^2 - V_stat^2)) - (1 + 2*n_vtest));
    p_val = min(p_val, 1);  % clamp

    vtest_predictions.(vtest_names{vv}).V = V_stat;
    vtest_predictions.(vtest_names{vv}).u = u_stat;
    vtest_predictions.(vtest_names{vv}).p = p_val;
    vtest_predictions.(vtest_names{vv}).mu_deg = vtest_angles(vv);

    fprintf('  V-test toward %s (%+.1f deg): V=%.1f, u=%.2f, p=%.2e\n', ...
        vtest_names{vv}, vtest_angles(vv), V_stat, u_stat, p_val)
end

% Find winner (highest V)
V_vals = [vtest_predictions.VS.V, vtest_predictions.IOC.V, ...
          vtest_predictions.C1.V, vtest_predictions.C2.V];
[~, winner_idx] = max(V_vals);
vtest_winner = vtest_names{winner_idx};
fprintf('  V-test winner: %s (V=%.1f, p=%.2e)\n', ...
    vtest_winner, V_vals(winner_idx), vtest_predictions.(vtest_winner).p)

%% gDSI Threshold Sweep
fprintf('Current stage: gDSI threshold sweep\n')

sweep_thresholds = [0.2, 0.3, 0.4, 0.5];
nThresh = length(sweep_thresholds);
sweep_n_cells  = zeros(1, nThresh);
sweep_n_VS     = zeros(1, nThresh);
sweep_n_IOC    = zeros(1, nThresh);
sweep_n_C1     = zeros(1, nThresh);
sweep_n_C2     = zeros(1, nThresh);
sweep_n_unc    = zeros(1, nThresh);
sweep_V_VS     = zeros(1, nThresh);
sweep_V_IOC    = zeros(1, nThresh);

for tt = 1:nThresh
    thr = sweep_thresholds(tt);
    ind_sweep = find(gDSI > thr & gDSI_plaid > thr);
    sweep_n_cells(tt) = length(ind_sweep);

    if sweep_n_cells(tt) < 3
        continue
    end

    dp_sweep = theta_grat_gDSI(ind_sweep) - theta_plaid_gDSI(ind_sweep);
    dp_sweep = mod(dp_sweep + 180, 360) - 180;

    clf_sweep = classifyDeltaPrefDir(dp_sweep, vs_angle_deg, ioc_angle_deg, c2_angle_deltapref, half_window_deg);
    sweep_n_VS(tt)  = clf_sweep.n_VS;
    sweep_n_IOC(tt) = clf_sweep.n_IOC;
    sweep_n_C1(tt)  = clf_sweep.n_C1;
    sweep_n_C2(tt)  = clf_sweep.n_C2;
    sweep_n_unc(tt) = clf_sweep.n_unclass + clf_sweep.n_pattern;

    % V-test at each threshold
    a_sweep = deg2rad(dp_sweep);
    n_sw = length(a_sweep);
    sweep_V_VS(tt)  = n_sw * mean(cos(a_sweep - deg2rad(vs_angle_deg)));
    sweep_V_IOC(tt) = n_sw * mean(cos(a_sweep - deg2rad(ioc_angle_deg)));

    fprintf('  gDSI > %.1f: n=%d, VS:%d IOC:%d C1:%d C2:%d unc:%d, V_VS=%.1f V_IOC=%.1f\n', ...
        thr, sweep_n_cells(tt), sweep_n_VS(tt), sweep_n_IOC(tt), sweep_n_C1(tt), sweep_n_C2(tt), ...
        sweep_n_unc(tt), sweep_V_VS(tt), sweep_V_IOC(tt))
end

figure('Visible', 'off', 'Position', [0 0 1200 900]);

% Top-left: Stacked bar of proportions
subplot(2, 2, 1)
prop_data = [sweep_n_VS; sweep_n_IOC; sweep_n_C1; sweep_n_unc + sweep_n_C2]';
bar_h = bar(prop_data, 'stacked');
bar_h(1).FaceColor = [0.2 0.4 0.8];  % VS
bar_h(2).FaceColor = [0.8 0.2 0.2];  % IOC
bar_h(3).FaceColor = [0.4 0.4 0.4];  % C1
bar_h(4).FaceColor = [0.75 0.75 0.75]; % unclass
xticks(1:nThresh)
xticklabels(arrayfun(@(x) sprintf('%.1f', x), sweep_thresholds, 'UniformOutput', false))
xlabel('gDSI threshold')
ylabel('Number of cells')
legend({'VS', 'IOC', 'C1', 'Unclass'}, 'Location', 'bestoutside', 'FontSize', 7)
title('Classification by threshold')

% Top-right: V-statistic for VS and IOC
subplot(2, 2, 2)
plot(sweep_thresholds, sweep_V_VS, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b')
hold on
plot(sweep_thresholds, sweep_V_IOC, 'r-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r')
xlabel('gDSI threshold')
ylabel('V-statistic')
legend({'V toward VS', 'V toward IOC'}, 'Location', 'best', 'FontSize', 8)
title('V-test strength vs threshold')

% Bottom-left: n_cells remaining
subplot(2, 2, 3)
bar(sweep_thresholds, sweep_n_cells, 0.5, 'FaceColor', [0.3 0.3 0.6])
xlabel('gDSI threshold')
ylabel('n cells (gDSI both > thr)')
title('Cells surviving threshold')
for tt = 1:nThresh
    text(sweep_thresholds(tt), sweep_n_cells(tt) + max(sweep_n_cells)*0.03, ...
        num2str(sweep_n_cells(tt)), 'HorizontalAlignment', 'center', 'FontSize', 9)
end

% Bottom-right: Mini polar histograms at each threshold
subplot(2, 2, 4)
hold on
colors_sweep = [0.2 0.4 0.8; 0.4 0.7 0.4; 0.9 0.6 0.1; 0.8 0.2 0.2];
for tt = 1:nThresh
    thr = sweep_thresholds(tt);
    ind_sw = find(gDSI > thr & gDSI_plaid > thr);
    if length(ind_sw) < 3, continue; end
    dp_sw = theta_grat_gDSI(ind_sw) - theta_plaid_gDSI(ind_sw);
    dp_sw = mod(dp_sw + 180, 360) - 180;
    [counts, edges] = histcounts(dp_sw, -180:11.25:180);
    centers = (edges(1:end-1) + edges(2:end)) / 2;
    % Normalize to density for comparison
    counts_norm = counts / sum(counts);
    plot(centers, counts_norm, '-', 'Color', colors_sweep(tt,:), 'LineWidth', 1.5)
end
xline(vs_angle_deg, 'b--', 'LineWidth', 1)
xline(ioc_angle_deg, 'r--', 'LineWidth', 1)
xline(0, 'k--', 'LineWidth', 0.5)
xlabel('\DeltaPrefDir (deg)')
ylabel('Normalized density')
legend([arrayfun(@(x) sprintf('gDSI>%.1f', x), sweep_thresholds, 'UniformOutput', false), ...
    {'VS', 'IOC', 'C1'}], 'Location', 'best', 'FontSize', 6)
title('Distribution shape vs threshold')
xlim([-180 180])

sgtitle({[mouse ' ' date ' — gDSI Threshold Sweep'], ...
    sprintf('Offset=%+d\\circ, half-window=%d\\circ', signed_offset(1), half_window_deg)}, 'FontSize', 10)

print(fullfile(outDir, ...
    [date '_' mouse '_' run_str '_gDSI_ThresholdSweep.png']), '-dpng', '-r300')

%% DSI vs gDSI scatter — colored by angular distance to VS / IOC
fprintf('Current stage: Plotting DSI vs gDSI scatter\n')

% Cell selection: responsive + ANOVA (no gDSI threshold — show full range)
ind_scatter = intersect(resp_ind_dir, p_dir);

% Delta pref dir for each cell (grating minus plaid, circular mean dirs)
delta_all = gDSI_prefDir(ind_scatter) - gDSI_plaid_prefDir(ind_scatter);
delta_all = mod(delta_all + 180, 360) - 180;  % wrap to [-180, 180]

% Circular angular distance to VS and IOC predictions
dist_VS  = abs(delta_all - vs_angle_deg);
dist_VS  = min(dist_VS, 360 - dist_VS);   % range [0, 180]
dist_IOC = abs(delta_all - ioc_angle_deg);
dist_IOC = min(dist_IOC, 360 - dist_IOC);  % range [0, 180]

figure('Visible','off', 'Position', [0 0 1440 720]);

% Left: colored by VS distance
subplot(1,2,1)
scatter(gDSI(ind_scatter), DSI(ind_scatter), 30, dist_VS, 'filled')
hold on
plot([0 1], [0 1], 'k--', 'LineWidth', 0.5)
colormap(subplot(1,2,1), flipud(jet))
cb = colorbar; cb.Label.String = 'Angular distance to VS (deg)';
clim([0 180])
xlim([0 1]); ylim([0 1]); axis square
xlabel('gDSI (grating)'); ylabel('Peak DSI (grating)')
title(sprintf('VS distance  (n=%d)', length(ind_scatter)))

% Right: colored by IOC distance
subplot(1,2,2)
scatter(gDSI(ind_scatter), DSI(ind_scatter), 30, dist_IOC, 'filled')
hold on
plot([0 1], [0 1], 'k--', 'LineWidth', 0.5)
colormap(subplot(1,2,2), flipud(jet))
cb = colorbar; cb.Label.String = 'Angular distance to IOC (deg)';
clim([0 180])
xlim([0 1]); ylim([0 1]); axis square
xlabel('gDSI (grating)'); ylabel('Peak DSI (grating)')
title(sprintf('IOC distance  (n=%d)', length(ind_scatter)))

sgtitle({[mouse ' ' date ' DSI vs gDSI (grating)'], ...
    ['resp + ANOVA, n=' num2str(length(ind_scatter)) '/' num2str(nCells)]}, 'FontSize', 10)

print(fullfile(outDir,...
    [date '_' mouse '_' run_str '_DSI_vs_gDSI_gradient_gDSI.png']), '-dpng', '-r300')

% %% direction tuning curves (ALIGNED — commented out)
% fprintf('Current stage: Plotting tuning curves\n')
%
% [avg_resp_grat, avg_resp_plaid] = getAlignedGratPlaidTuning(avg_resp_dir, '', plaid_offset);
%
% figure('Visible','off', 'Position', [0 0 1920 1080]);
% start = 1;
% n = 1;
% step = 360/nStimDir; x = -(180-step/2):step:(180-step/2);
% for iCell = 1:nCells
%     subplot(5,4,start)
%         for im = 1:nMaskPhas
%             plot(x, avg_resp_plaid(iCell,:,im))
%             hold on
%         end
%         plot(x, avg_resp_grat(iCell,:),'k')
%         if iCell ==1; legend([arrayfun(@(x) [num2str(x) ' deg'], maskPhas, 'UniformOutput', false), {'gratings'}]); end;
%         xlabel('plaid direction')
%         ylabel('df/f')
%         idx = iCell==ind;
%         if any(idx)
%             subtitle(['cell ' num2str(iCell)],'fontweight','bold')
%         else
%             subtitle(['cell ' num2str(iCell)]);
%         end
%     start = start+1;
%     if start >20
%         sgtitle([mouse ' ' date ' Grating (black) + Plaid (colored) tuning curves, aligned to pref grating dir'], 'FontSize', 10)
%         print(fullfile(outDir,[date '_' mouse '_' run_str '_DirectionTuningGratingsAndPlaids_' num2str(n) '.png']),'-dpng', '-r300')
%         figure('Visible','off', 'Position', [0 0 1920 1080]);
%
%         start = 1;
%         n = n+1;
%     end
%     if iCell == nCells
%         sgtitle([mouse ' ' date ' Grating (black) + Plaid (colored) tuning curves, aligned to pref grating dir'], 'FontSize', 10)
%         print(fullfile(outDir,[date '_' mouse '_' run_str '_DirectionTuningGratingsAndPlaids_' num2str(n) '.png']), '-dpng', '-r300')
%     end
% end
% %close all
%
% %% Grating tuning curves (separate, ALIGNED — commented out)
%
% figure('Visible','off', 'Position', [0 0 1920 1080]);
% start = 1;
% n = 1;
% step = 360/nStimDir; x = -(180-step/2):step:(180-step/2);
% for iCell = 1:nCells
%     subplot(5,4,start)
%         plot(x, avg_resp_grat(iCell,:),'k')
%         xlabel('direction')
%         ylabel('df/f')
%         idx = iCell==ind;
%         if any(idx)
%             subtitle(['cell ' num2str(iCell)],'fontweight','bold')
%         else
%             subtitle(['cell ' num2str(iCell)])
%         end
%     start = start+1;
%     if start >20
%         sgtitle([mouse ' ' date ' Grating-only tuning curves (df/f vs direction, aligned to pref dir)'], 'FontSize', 10)
%         print(fullfile(outDir,[date '_' mouse '_' run_str '_GratingTuningAligned_' num2str(n) '.png']),'-dpng', '-r300')
%         figure('Visible','off', 'Position', [0 0 1920 1080]);
%
%         start = 1;
%         n = n+1;
%     end
%     if iCell == nCells
%         sgtitle([mouse ' ' date ' Grating-only tuning curves (df/f vs direction, aligned to pref dir)'], 'FontSize', 10)
%         print(fullfile(outDir,[date '_' mouse '_' run_str '_GratingTuningAligned_' num2str(n) '.png']), '-dpng', '-r300')
%     end
% end
% %close all
%
% %% Plaid tuning curves (separate, ALIGNED — commented out)
%
% figure('Visible','off', 'Position', [0 0 1920 1080]);
% start = 1;
% n = 1;
% step = 360/nStimDir; x = -(180-step/2):step:(180-step/2);
% for iCell = 1:nCells
%     subplot(5,4,start)
%         for im = 1:nMaskPhas
%             plot(x, avg_resp_plaid(iCell,:,im))
%             hold on
%         end
%         if iCell == 1; legend(arrayfun(@(x) [num2str(x) ' deg'], maskPhas, 'UniformOutput', false)); end;
%         xlabel('plaid direction')
%         ylabel('df/f')
%         idx = iCell==ind;
%         if any(idx)
%             subtitle(['cell ' num2str(iCell)],'fontweight','bold')
%         else
%             subtitle(['cell ' num2str(iCell)])
%         end
%     start = start+1;
%     if start >20
%         sgtitle([mouse ' ' date ' Plaid-only tuning curves by mask phase (df/f vs direction, aligned to pref grating dir)'], 'FontSize', 10)
%         print(fullfile(outDir,[date '_' mouse '_' run_str '_PlaidTuningAligned_' num2str(n) '.png']),'-dpng', '-r300')
%         figure('Visible','off', 'Position', [0 0 1920 1080]);
%
%         start = 1;
%         n = n+1;
%     end
%     if iCell == nCells
%         sgtitle([mouse ' ' date ' Plaid-only tuning curves by mask phase (df/f vs direction, aligned to pref grating dir)'], 'FontSize', 10)
%         print(fullfile(outDir,[date '_' mouse '_' run_str '_PlaidTuningAligned_' num2str(n) '.png']), '-dpng', '-r300')
%     end
% end
% %close all
%
% %% Polar plots by individual cell (ALIGNED — commented out)
% fprintf('Current stage: Plotting polar plots\n')
%
% figure('Visible','off', 'Position', [0 0 1920 1080]);
% start = 1;
% n = 1;
%
% step = 360/nStimDir; x = -(180-step/2):step:(180-step/2);
% x_rad = deg2rad(x);
% for iCell =1:nCells
%     subplot(5,4,start)
%         for im = 1:nMaskPhas
%             polarplot([x_rad x_rad(1)], [avg_resp_plaid(iCell,:,im) avg_resp_plaid(iCell,1,im)])
%             hold on
%         end
%         polarplot([x_rad x_rad(1)], [avg_resp_grat(iCell,:) avg_resp_grat(iCell,1)],'k', 'LineWidth',2)
%         % Weighted circular mean preferred direction (theta_pref = arg(sum(R_theta * e^(i*theta))))
%         rmax = max(rlim);
%         vec_grat = sum(avg_resp_grat(iCell,:) .* exp(1i * x_rad));
%         polarplot([angle(vec_grat) angle(vec_grat)], [0 rmax], 'k-', 'LineWidth', 1.5)
%         colors = get(gca, 'ColorOrder');
%         for im = 1:nMaskPhas
%             vec_plaid = sum(avg_resp_plaid(iCell,:,im) .* exp(1i * x_rad));
%             polarplot([angle(vec_plaid) angle(vec_plaid)], [0 rmax], '-', 'Color', colors(im,:), 'LineWidth', 1.5)
%         end
%         idx = iCell==ind;
%         if any(idx)
%             subtitle(['cell ' num2str(iCell) ', gDSI ' num2str(DSI(iCell),'%.2f') ', pDSI ' num2str(DSI_plaid(iCell),'%.2f')],'fontweight','bold')
%         else
%             subtitle(['cell ' num2str(iCell) ', gDSI ' num2str(DSI(iCell),'%.2f') ', pDSI ' num2str(DSI_plaid(iCell),'%.2f')])
%         end
%     start = start+1;
%     if start>20
%         sgtitle([mouse ' ' date ' Polar plots: grating (black) + plaid (colored), aligned to pref grating dir'], 'FontSize', 10)
%         print(fullfile(outDir,[date '_' mouse '_' run_str '_PolarPlots_' num2str(n) '.png']), '-dpng', '-r300')
%         figure('Visible','off', 'Position', [0 0 1920 1080]);
%
%         start = 1;
%         n = n+1;
%     end
%     if iCell == nCells
%         sgtitle([mouse ' ' date ' Polar plots: grating (black) + plaid (colored), aligned to pref grating dir'], 'FontSize', 10)
%         print(fullfile(outDir,[date '_' mouse '_' run_str '_PolarPlots_' num2str(n) '.png']), '-dpng', '-r300')
%     end
% end
% %close all

% %% ZcZp by individual cell
% figure('Visible','off', 'Position', [0 0 1920 1080]);
% start = 1;
% n = 1;
%
% for iCell = 1:nCells
%     subplot(5,4,start)
%         for im = 1:nMaskPhas
%             scatter(Zc(im,iCell), Zp(im,iCell))
%             hold on
%         end
%         ylabel('Zp'); ylim([-4 8]);
%         xlabel('Zc'); xlim([-4 8]);
%         if iCell ==1; legend(arrayfun(@(x) [num2str(x) ' deg'], maskPhas, 'UniformOutput', false)); end;
%         idx = iCell==ind;
%         if any(idx)
%             subtitle(['cell ' num2str(iCell)],'fontweight','bold')
%         else
%             subtitle(['cell ' num2str(iCell)])
%         end
%         plotZcZpBorders
%     start = start+1;
%     if start>20
%         sgtitle([mouse ' ' date ' - Zp Zc by cell'])
%         print(fullfile(outDir,[date '_' mouse '_' run_str '_ZcZpByCell_' num2str(n) '.png']), '-dpng', '-r300')
%         figure('Visible','off', 'Position', [0 0 1920 1080]);
% 
%         start = 1;
%         n = n+1;
%     end
%     if iCell == nCells
%         sgtitle([mouse ' ' date ' - Zp Zc by cell'])
%         print(fullfile(outDir,[date '_' mouse '_' run_str '_ZcZpByCell_' num2str(n) '.png']), '-dpng', '-r300')
%     end
% end
% %close all

%% Trials per stimulus condition
fprintf('Current stage: Plotting trials per stimulus\n')

figure('Visible','off', 'Position', [0 0 1920 1080]);
subplot(2,1,1)
    bar(stimDirs, squeeze(trialsperstim(:,:,1)))
    title('Number of trials per direction — Gratings')
    ylabel('# of trials')
    xlabel('direction (deg)')
subplot(2,1,2)
    bar(stimDirs, squeeze(trialsperstim(:,:,2)))
    title('Number of trials per direction — Plaids (by mask phase)')
    ylabel('# of trials')
    xlabel('direction (deg)')
    legend(arrayfun(@(x) num2str(x), maskPhas, 'UniformOutput', false))
    
print(fullfile(outDir,[date '_' mouse '_' run_str '_TrialsPerStim.png']),'-dpng', '-r300')       

%% look at df/f across individual trials with comparison to statistical tests (ttest for grating, anova for plaid)
fprintf('Current stage: Plotting single-trial responses\n')


figure('Visible','off', 'Position', [0 0 1920 1080]);
start=1;
n=1;
for iCell = 1:nCells
    subplot(5,4,start)
        for iDir = 1:nStimDir
            scatter(stimDirs(iDir)*ones(1,size(resp_cell{iDir,1,1},2)),resp_cell{iDir,1,1}(iCell,:));
            hold on
            max_df = max(resp_cell{iDir,1,1}(iCell,:),[],2); %find max df/f to use to print t-test result
            tt = h_resp(iCell,iDir,1,1);
            if tt == 1; text(stimDirs(iDir), max_df,'*'); end;
        end
        idx = iCell==ind;
        if any(idx)
            subtitle(['cell ' num2str(iCell) ', p=' num2str(p_anova_dir(iCell),'%.3f')],'fontweight','bold')
        else
            subtitle(['cell ' num2str(iCell) ', p=' num2str(p_anova_dir(iCell),'%.3f')])
        end
        ylabel('df/f'); xlabel('direction');
    start = start+1;    
    if start>20
        sgtitle([mouse ' ' date ' Single-trial grating responses by direction (* = sig t-test vs blank, subtitle = ANOVA p)'], 'FontSize', 10)
        print(fullfile(outDir,[date '_' mouse '_' run_str '_GratingResponses_' num2str(n) '.png']), '-dpng', '-r300')
        figure('Visible','off', 'Position', [0 0 1920 1080]);

        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle([mouse ' ' date ' Single-trial grating responses by direction (* = sig t-test vs blank, subtitle = ANOVA p)'], 'FontSize', 10)
        print(fullfile(outDir,[date '_' mouse '_' run_str '_GratingResponses_' num2str(n) '.png']), '-dpng', '-r300')
    end     
end
    %close all
   
   
figure('Visible','off', 'Position', [0 0 1920 1080]);
start=1;
n=1;
for iCell = 1:nCells
    subplot(5,4,start)
        [p_min, p_ind] = min(p_anova_plaid(:,iCell),[],1);
        for iDir = 1:nStimDir
            scatter(stimDirs(iDir)*ones(1,size(resp_cell{iDir,p_ind,2},2)),resp_cell{iDir,p_ind,2}(iCell,:));
            hold on
        end
        idx = iCell==ind;
        if any(idx)
            subtitle(['cell ' num2str(iCell) ', p=' num2str(p_min,'%.3f')],'fontweight','bold')
        else
            subtitle(['cell ' num2str(iCell) ', p=' num2str(p_min,'%.3f')])
        end
        ylabel('df/f'); xlabel('direction');
    start = start+1;    
    if start>20
        sgtitle([mouse ' ' date ' Single-trial plaid responses by direction (best mask phase by ANOVA)'], 'FontSize', 10)
        print(fullfile(outDir,[date '_' mouse '_' run_str '_PlaidResponses_' num2str(n) '.png']), '-dpng', '-r300')
        figure('Visible','off', 'Position', [0 0 1920 1080]);

        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle([mouse ' ' date ' Single-trial plaid responses by direction (best mask phase by ANOVA)'], 'FontSize', 10)
        print(fullfile(outDir,[date '_' mouse '_' run_str '_PlaidResponses_' num2str(n) '.png']), '-dpng', '-r300')
    end        
end
  %close all

%% DSI histograms
fprintf('Current stage: Plotting DSI histograms\n')

% --- Old DSI histograms (peak, commented out) ---
% figure('Visible','off', 'Position', [0 0 1920 1080]);
% subplot(2,1,1)
%     histogram(DSI)
%     hold on
%     xline(0.5,'r--','LineWidth',1.5)
%     xlabel('DSI')
%     ylabel('# cells')
%     title(['Grating DSI (peak) — ' num2str(length(DSI_ind)) '/' num2str(nCells) ' cells > 0.5'])
% subplot(2,1,2)
%     histogram(DSI_plaid)
%     hold on
%     xline(0.5,'r--','LineWidth',1.5)
%     xlabel('DSI')
%     ylabel('# cells')
%     title(['Plaid DSI (peak) — ' num2str(length(DSI_plaid_ind)) '/' num2str(nCells) ' cells > 0.5'])
% sgtitle([mouse ' ' date ' DSI distributions: grating vs plaid (peak DSI)'], 'FontSize', 10)
% print(fullfile(outDir,[date '_' mouse '_' run_str '_DSIhistograms_peak.png']),'-dpng', '-r300')

% --- Old DSI histograms (weighted, commented out — gDSI histogram below covers this) ---
% figure('Visible','off', 'Position', [0 0 1920 1080]);
% subplot(2,1,1)
%     histogram(DSI_w)
%     hold on
%     xline(0.5,'r--','LineWidth',1.5)
%     xlabel('DSI')
%     ylabel('# cells')
%     title(['Grating DSI (weighted) — ' num2str(length(DSI_w_ind)) '/' num2str(nCells) ' cells > 0.5'])
% subplot(2,1,2)
%     histogram(DSI_w_plaid)
%     hold on
%     xline(0.5,'r--','LineWidth',1.5)
%     xlabel('DSI')
%     ylabel('# cells')
%     title(['Plaid DSI (weighted) — ' num2str(length(DSI_w_plaid_ind)) '/' num2str(nCells) ' cells > 0.5'])
% sgtitle([mouse ' ' date ' DSI distributions: grating vs plaid (weighted DSI)'], 'FontSize', 10)
% print(fullfile(outDir,[date '_' mouse '_' run_str '_DSIhistograms_weighted.png']),'-dpng', '-r300')

% --- New gDSI histograms (Li et al. 2025) ---
figure('Visible','off', 'Position', [0 0 1920 1080]);
subplot(2,1,1)
    histogram(gDSI)
    hold on
    xline(0.2,'r--','LineWidth',1.5)
    xlabel('gDSI')
    ylabel('# cells')
    title(['Grating gDSI — ' num2str(length(gDSI_ind)) '/' num2str(nCells) ' cells > 0.2'])
subplot(2,1,2)
    histogram(gDSI_plaid)
    hold on
    xline(0.2,'r--','LineWidth',1.5)
    xlabel('gDSI')
    ylabel('# cells')
    title(['Plaid gDSI — ' num2str(length(gDSI_plaid_ind)) '/' num2str(nCells) ' cells > 0.2'])
sgtitle([mouse ' ' date ' gDSI distributions (Li et al. 2025)'], 'FontSize', 10)
print(fullfile(outDir,[date '_' mouse '_' run_str '_gDSI_histograms.png']),'-dpng', '-r300')

%% VS/IOC-resolved Zp/Zc scatter figures
fprintf('Current stage: Plotting Zp/Zc VSIOC figures\n')

if ~ZpZcVSIOC.is_degenerate
    % Figure: Zp_VS vs Zp_IOC scatter, colored by classification
    figure('Visible','off', 'Position', [0 0 720 720]);
    hold on
    plot([-4 8], [-4 8], 'k--', 'LineWidth', 0.5)
    plot([-4 8], [1.28 1.28], 'k:', 'LineWidth', 0.5)
    plot([1.28 1.28], [-4 8], 'k:', 'LineWidth', 0.5)

    % Unclassified (all cells first, then overlay classified)
    scatter(ZpZcVSIOC.Zp_VS(1,:), ZpZcVSIOC.Zp_IOC(1,:), 20, [0.75 0.75 0.75], 'filled', 'MarkerFaceAlpha', 0.5)
    if ~isempty(ZpZcVSIOC.ind_VS_PDS)
        scatter(ZpZcVSIOC.Zp_VS(1,ZpZcVSIOC.ind_VS_PDS), ZpZcVSIOC.Zp_IOC(1,ZpZcVSIOC.ind_VS_PDS), ...
            30, [0.2 0.4 0.8], 'filled')
    end
    if ~isempty(ZpZcVSIOC.ind_IOC_PDS)
        scatter(ZpZcVSIOC.Zp_VS(1,ZpZcVSIOC.ind_IOC_PDS), ZpZcVSIOC.Zp_IOC(1,ZpZcVSIOC.ind_IOC_PDS), ...
            30, [0.8 0.2 0.2], 'filled')
    end
    if ~isempty(ZpZcVSIOC.ind_CDS)
        scatter(ZpZcVSIOC.Zp_VS(1,ZpZcVSIOC.ind_CDS), ZpZcVSIOC.Zp_IOC(1,ZpZcVSIOC.ind_CDS), ...
            30, [0.4 0.4 0.4], 'filled')
    end

    xlabel('Zp_{VS}'); ylabel('Zp_{IOC}')
    xlim([-4 8]); ylim([-4 8]); axis square
    legend({'y=x', '', 'all cells', ...
        sprintf('VS-pattern (n=%d)', length(ZpZcVSIOC.ind_VS_PDS)), ...
        sprintf('IOC-pattern (n=%d)', length(ZpZcVSIOC.ind_IOC_PDS)), ...
        sprintf('Component (n=%d)', length(ZpZcVSIOC.ind_CDS))}, ...
        'Location', 'best', 'FontSize', 7)
    sgtitle({[mouse ' ' date ' Zp_{VS} vs Zp_{IOC} (phase 1)'], ...
        sprintf('VS-PDS:%d, IOC-PDS:%d, CDS:%d / %d cells', ...
        length(ZpZcVSIOC.ind_VS_PDS), length(ZpZcVSIOC.ind_IOC_PDS), ...
        length(ZpZcVSIOC.ind_CDS), nCells)}, 'FontSize', 10)

    print(fullfile(outDir,[date '_' mouse '_' run_str '_ZpVS_vs_ZpIOC_gDSI.png']), '-dpng', '-r300')

    % Figure: Side-by-side Zp_VS vs Zc and Zp_IOC vs Zc
    figure('Visible','off', 'Position', [0 0 1440 720]);
    subplot(1,2,1)
        hold on
        scatter(ZpZcVSIOC.Zc(1,:), ZpZcVSIOC.Zp_VS(1,:), 20, [0.75 0.75 0.75], 'filled', 'MarkerFaceAlpha', 0.5)
        if ~isempty(ZpZcVSIOC.ind_VS_PDS)
            scatter(ZpZcVSIOC.Zc(1,ZpZcVSIOC.ind_VS_PDS), ZpZcVSIOC.Zp_VS(1,ZpZcVSIOC.ind_VS_PDS), ...
                30, [0.2 0.4 0.8], 'filled')
        end
        if ~isempty(ZpZcVSIOC.ind_CDS)
            scatter(ZpZcVSIOC.Zc(1,ZpZcVSIOC.ind_CDS), ZpZcVSIOC.Zp_VS(1,ZpZcVSIOC.ind_CDS), ...
                30, [0.4 0.4 0.4], 'filled')
        end
        plotZcZpBorders
        xlabel('Zc'); ylabel('Zp_{VS}')
        xlim([-4 8]); ylim([-4 8]); axis square
        title(sprintf('VS-pattern: %d PDS, %d CDS', length(ZpZcVSIOC.ind_VS_PDS), length(ZpZcVSIOC.ind_CDS)))
    subplot(1,2,2)
        hold on
        scatter(ZpZcVSIOC.Zc(1,:), ZpZcVSIOC.Zp_IOC(1,:), 20, [0.75 0.75 0.75], 'filled', 'MarkerFaceAlpha', 0.5)
        if ~isempty(ZpZcVSIOC.ind_IOC_PDS)
            scatter(ZpZcVSIOC.Zc(1,ZpZcVSIOC.ind_IOC_PDS), ZpZcVSIOC.Zp_IOC(1,ZpZcVSIOC.ind_IOC_PDS), ...
                30, [0.8 0.2 0.2], 'filled')
        end
        if ~isempty(ZpZcVSIOC.ind_CDS)
            scatter(ZpZcVSIOC.Zc(1,ZpZcVSIOC.ind_CDS), ZpZcVSIOC.Zp_IOC(1,ZpZcVSIOC.ind_CDS), ...
                30, [0.4 0.4 0.4], 'filled')
        end
        plotZcZpBorders
        xlabel('Zc'); ylabel('Zp_{IOC}')
        xlim([-4 8]); ylim([-4 8]); axis square
        title(sprintf('IOC-pattern: %d PDS, %d CDS', length(ZpZcVSIOC.ind_IOC_PDS), length(ZpZcVSIOC.ind_CDS)))

    sgtitle([mouse ' ' date ' Zp vs Zc (VS/IOC resolved, phase 1)'], 'FontSize', 10)
    print(fullfile(outDir,[date '_' mouse '_' run_str '_ZpZc_VSIOC_gDSI.png']), '-dpng', '-r300')
else
    fprintf('  Skipping Zp/Zc VSIOC figures (VS ≈ IOC, degenerate case)\n')
end

fprintf('Done.\n')
