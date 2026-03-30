%% df/f Trace Viewer — trial-averaged traces per cell, colored by condition
% Standalone script: loads raw data, computes trial-aligned df/f,
% plots grating (black), plaid (blue), blank (gray) traces per cell.

clc; clear all; close all;
startup
ds = 'CrossOriRandDir_ExptList_DC';
eval(ds)
iexp = 3;
frame_rate = 15;

%% Path construction & data loading
mouse = expt(iexp).mouse;
date = expt(iexp).date;
ImgFolder = expt(iexp).coFolder;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

dataBase = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\' expt(iexp).saveLoc];
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\David';

sessionDir = fullfile(dataBase, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]);

fprintf([mouse ' ' date '\n'])
fprintf('Loading data...\n')
load(fullfile(sessionDir, [date '_' mouse '_' run_str '_TCs.mat']))
load(fullfile(sessionDir, [date '_' mouse '_' run_str '_dataStim.mat']))
load(fullfile(sessionDir, [date '_' mouse '_' run_str '_input.mat']))

%% Compute trial-aligned df/f
prewin_frames = unique(celleqel2mat_padded(input.tItiWaitFrames))./3;
postwin_frames = unique(celleqel2mat_padded(input.nStimOneFramesOn));
tt = (1-prewin_frames:postwin_frames).*(1/frame_rate);

data_resp = nan(prewin_frames+postwin_frames, nCells, nTrials);
for itrial = 1:nTrials
    if cStimOn(itrial) + postwin_frames < sz(3)
        data_resp(:,:,itrial) = npSub_tc(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
    end
end

data_f = mean(data_resp(1:prewin_frames,:,:),1);
data_dfof_tc = (data_resp - data_f) ./ data_f;

%% Separate trial indices by condition
ind_stimAlone = intersect(find(stimCon_all), find(maskCon_all==0));
ind_plaid     = intersect(find(stimCon_all), find(maskCon_all));
ind_blank     = intersect(find(stimCon_all==0), find(maskCon_all==0));

%% Compute per-direction, per-condition trial-averaged df/f (per cell)
nDirs = length(stimDirs);
nFrames = prewin_frames + postwin_frames;

% dfof_*_byDir: (frames x nCells x nDirs) — trial-averaged per cell per dir
dfof_grat_byDir  = nan(nFrames, nCells, nDirs);
dfof_plaid_byDir = nan(nFrames, nCells, nDirs);

for iDir = 1:nDirs
    ind_dir = find(stimDir_all == stimDirs(iDir));
    ind_grat_dir  = intersect(ind_dir, ind_stimAlone);
    ind_plaid_dir = intersect(ind_dir, ind_plaid);

    if ~isempty(ind_grat_dir)
        dfof_grat_byDir(:,:,iDir) = nanmean(data_dfof_tc(:,:,ind_grat_dir), 3);
    end
    if ~isempty(ind_plaid_dir)
        dfof_plaid_byDir(:,:,iDir) = nanmean(data_dfof_tc(:,:,ind_plaid_dir), 3);
    end
end

dfof_blank = squeeze(nanmean(data_dfof_tc(:,:,ind_blank), 3));  % no direction

%% Plot: one subplot per direction + one for blank
nSubplots = nDirs + 1;
nCol = ceil(sqrt(nSubplots));
nRow = ceil(nSubplots / nCol);

figure('Visible','off', 'Position', [0 0 1920 1080])

for iDir = 1:nDirs
    subplot(nRow, nCol, iDir)
    hold on

    % Individual cell traces (thin, transparent)
    for iCell = 1:nCells
        plot(tt, dfof_grat_byDir(:,iCell,iDir), 'k', 'LineWidth', 0.3, 'Color', [0 0 0 0.15])
        plot(tt, dfof_plaid_byDir(:,iCell,iDir), 'b', 'LineWidth', 0.3, 'Color', [0 0 1 0.15])
    end

    % Population mean (thick)
    plot(tt, nanmean(dfof_grat_byDir(:,:,iDir), 2), 'k', 'LineWidth', 1.5)
    plot(tt, nanmean(dfof_plaid_byDir(:,:,iDir), 2), 'b', 'LineWidth', 1.5)

    xline(0, 'r--')
    title([num2str(stimDirs(iDir)) '\circ'], 'FontSize', 8)
    set(gca, 'FontSize', 5)
    hold off
end

% Blank subplot
subplot(nRow, nCol, nDirs + 1)
hold on
for iCell = 1:nCells
    plot(tt, dfof_blank(:,iCell), 'Color', [0.6 0.6 0.6 0.15], 'LineWidth', 0.3)
end
plot(tt, nanmean(dfof_blank, 2), 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5)
xline(0, 'r--')
title('blank', 'FontSize', 8)
set(gca, 'FontSize', 5)
hold off

% Legend on first subplot
subplot(nRow, nCol, 1)
hold on
h1 = plot(nan, nan, 'k', 'LineWidth', 1.5);
h2 = plot(nan, nan, 'b', 'LineWidth', 1.5);
h3 = plot(nan, nan, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
legend([h1 h2 h3], {'grating','plaid','blank'}, 'FontSize', 5, 'Location', 'best')

sgtitle({[mouse ' ' date ' — Trial-averaged df/f traces by direction (' num2str(nCells) ' cells)'], ...
    ['grating: ' num2str(length(ind_stimAlone)) ' trials, ' ...
     'plaid: ' num2str(length(ind_plaid)) ' trials, ' ...
     'blank: ' num2str(length(ind_blank)) ' trials']}, 'FontSize', 12)

%% Save figure
outDir = fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], 'dfof_traces');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

print(fullfile(outDir, [date '_' mouse '_' run_str '_dfofTraces.png']), '-dpng', '-r300')
fprintf('Saved to %s\n', outDir)
