%% speedMixing_comparison_DC.m
% Compares grating preferred direction and grat-vs-plaid scatter using
% stim-only (TF_stim) vs mask-only (TF_mask) grating tuning curves.
%
% Figures:
%   1) Grat vs Plaid pref dir scatter — stim-only grating tuning
%   2) Grat vs Plaid pref dir scatter — mask-only grating tuning
%   3) Stim-only pref dir vs mask-only pref dir (cell-by-cell)
%
% Output: Analysis\2P\[date]_[mouse]\...\gDSI\speedMixing\
% Dependencies: CrossOriRandDir_ExptList_DC, getGDSI, preprocessed data

%% Setup & data loading (from exptAnalysis_DC_gDSI.m)
clc; clear all; close all; clear all global;
startup
doRedChannel = 0;
ds = 'CrossOriRandDir_ExptList_DC';
eval(ds)
nexp = length(expt);
iexp = 3;
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
fprintf('=== Speed-Mixing Comparison ===\n')

%% Load data
fprintf('Loading data...\n')
load(fullfile(dataBase, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
load(fullfile(dataBase, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
load(fullfile(dataBase, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))

%% Derive plaid offset from data
signed_offset = unique(maskDiffs(maskDiffs ~= 0));
if isempty(signed_offset)
    signed_offset = 45;
end
plaid_offset = abs(signed_offset(1));
offset_rad = deg2rad(signed_offset(1));
fprintf('Plaid offset: %+.0f deg (mask - stim direction)\n', signed_offset(1))

% Component speeds for VA/IOC prediction lines
SF_val = expt(iexp).SF;
if isfield(expt, 'TF_stim')
    speed_stim = expt(iexp).TF_stim / SF_val;
    speed_mask = expt(iexp).TF_mask / SF_val;
else
    speed_stim = expt(iexp).TF / SF_val;
    speed_mask = speed_stim;
    fprintf('WARNING: TF_stim/TF_mask not set — speeds are equal, speed mixing is not an issue.\n')
end
fprintf('Speeds: stim=%.1f deg/s, mask=%.1f deg/s\n', speed_stim, speed_mask)

% VS/IOC angular shifts (radians) — stimulus property, computed once
alpha_VS = atan2(speed_mask * sin(offset_rad), speed_stim + speed_mask * cos(offset_rad));
alpha_IOC = atan2(speed_mask - speed_stim * cos(offset_rad), speed_stim * sin(offset_rad));
shift_VS_deg = -rad2deg(alpha_VS);
shift_IOC_deg = -rad2deg(alpha_IOC);
fprintf('VS shift: %+.2f deg, IOC shift: %+.2f deg\n', shift_VS_deg, shift_IOC_deg)

%% Response extraction
fprintf('Extracting responses...\n')
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

resp_win = prewin_frames+5:prewin_frames+nFramesOn;

%% Build stim-only and mask-only tuning curves
fprintf('Building tuning curves (stim-only, mask-only)...\n')

tc_stim = zeros(nCells, nStimDir);   % stim-alone trials only (TF_stim)
tc_mask = zeros(nCells, nStimDir);   % mask-alone trials only (TF_mask)
tc_plaid = zeros(nCells, nStimDir);  % plaid trials (same for both)

for iDir = 1:nStimDir
    ind_stimdir = find(stimDir_all == stimDirs(iDir));
    ind_maskdir = find(maskDir_all == stimDirs(iDir));

    % Stim-only: stim direction matches this bin, stim-alone condition
    ind_stim_only = intersect(ind_stimdir, ind_stimAlone);
    % Mask-only: mask direction matches this bin, mask-alone condition
    ind_mask_only = intersect(ind_maskdir, ind_maskAlone);
    % Plaid: stim direction matches, plaid trials
    ind_plaid_dir = intersect(ind_stimdir, ind_plaid);

    if ~isempty(ind_stim_only)
        tc_stim(:,iDir) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind_stim_only),1),3));
    end
    if ~isempty(ind_mask_only)
        tc_mask(:,iDir) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind_mask_only),1),3));
    end
    if ~isempty(ind_plaid_dir)
        tc_plaid(:,iDir) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind_plaid_dir),1),3));
    end
end

%% Compute gDSI + preferred direction
fprintf('Computing gDSI...\n')

gDSI_stim_struct  = getGDSI(tc_stim, stimDirs);
gDSI_mask_struct  = getGDSI(tc_mask, stimDirs);
gDSI_plaid_struct = getGDSI(tc_plaid, stimDirs);

prefDir_stim  = gDSI_stim_struct.prefDir_deg;
prefDir_mask  = gDSI_mask_struct.prefDir_deg;
prefDir_plaid = gDSI_plaid_struct.prefDir_deg;

%% Output directory
outDir = fullfile(base, 'Analysis\2P', [date '_' mouse], ...
    [date '_' mouse '_' run_str], 'gDSI', 'speedMixing');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% Figure 1: Grat vs Plaid pref dir — STIM-ONLY grating tuning
fprintf('Plotting Figure 1: Grat vs Plaid (stim-only)...\n')

figure('Visible','off', 'Position', [0 0 1080 1080]);
hold on
plot([0 360], [0 360], 'k--', 'LineWidth', 1)
plot([0 360], [0+shift_VS_deg 360+shift_VS_deg], 'g--', 'LineWidth', 1)
plot([0 360], [0+shift_IOC_deg 360+shift_IOC_deg], 'm--', 'LineWidth', 1)
scatter(prefDir_stim, prefDir_plaid, 30, 'k', 'LineWidth', 0.8)
xlabel('Grating preferred direction (deg)')
ylabel('Plaid preferred direction (deg)')
xlim([0 360]); ylim([0 360])
axis square
legend({'component (y=x)', ['VS (y=x' sprintf('%+.1f', shift_VS_deg) ')'], ...
        ['IOC (y=x' sprintf('%+.1f', shift_IOC_deg) ')'], ...
        'all cells'}, 'Location', 'best', 'FontSize', 7)
sgtitle([mouse ' ' date ' Grating vs Plaid pref dir — STIM-ONLY (TF=' num2str(expt(iexp).TF_stim) ')'], 'FontSize', 12)
print(fullfile(outDir, [date '_' mouse '_' run_str '_GratVsPlaid_stimOnly.png']), '-dpng', '-r300')

%% Figure 2: Grat vs Plaid pref dir — MASK-ONLY grating tuning
fprintf('Plotting Figure 2: Grat vs Plaid (mask-only)...\n')

figure('Visible','off', 'Position', [0 0 1080 1080]);
hold on
plot([0 360], [0 360], 'k--', 'LineWidth', 1)
plot([0 360], [0+shift_VS_deg 360+shift_VS_deg], 'g--', 'LineWidth', 1)
plot([0 360], [0+shift_IOC_deg 360+shift_IOC_deg], 'm--', 'LineWidth', 1)
scatter(prefDir_mask, prefDir_plaid, 30, 'k', 'LineWidth', 0.8)
xlabel('Grating preferred direction (deg)')
ylabel('Plaid preferred direction (deg)')
xlim([0 360]); ylim([0 360])
axis square
legend({'component (y=x)', ['VS (y=x' sprintf('%+.1f', shift_VS_deg) ')'], ...
        ['IOC (y=x' sprintf('%+.1f', shift_IOC_deg) ')'], ...
        'all cells'}, 'Location', 'best', 'FontSize', 7)
sgtitle([mouse ' ' date ' Grating vs Plaid pref dir — MASK-ONLY (TF=' num2str(expt(iexp).TF_mask) ')'], 'FontSize', 12)
print(fullfile(outDir, [date '_' mouse '_' run_str '_GratVsPlaid_maskOnly.png']), '-dpng', '-r300')

%% Figure 3: Cell-by-cell pref dir change (stim-only vs mask-only)
fprintf('Plotting Figure 3: Stim vs Mask pref dir (cell-by-cell)...\n')

figure('Visible','off', 'Position', [0 0 1080 1080]);
hold on
plot([0 360], [0 360], 'k--', 'LineWidth', 1)
scatter(prefDir_stim, prefDir_mask, 30, 'k', 'filled', 'MarkerFaceAlpha', 0.5)
xlabel(['Preferred direction — stim-only (TF=' num2str(expt(iexp).TF_stim) ', deg)'])
ylabel(['Preferred direction — mask-only (TF=' num2str(expt(iexp).TF_mask) ', deg)'])
xlim([0 360]); ylim([0 360])
axis square
legend({'unity (no change)'}, 'Location', 'best', 'FontSize', 8)
sgtitle([mouse ' ' date ' — Preferred direction: stim-only vs mask-only (n=' num2str(nCells) ' cells)'], 'FontSize', 12)
print(fullfile(outDir, [date '_' mouse '_' run_str '_prefDir_stimVsMask.png']), '-dpng', '-r300')

%% Done
fprintf('\nDone. Figures saved to:\n  %s\n', outDir)
