%% dF/F response maps by direction — diagnostic viewer
% Shows 4x4 Gaussian-filtered dF/F maps for gratings and plaids separately.
% Uses same ExptList and path conventions as exptAnalysis.

clc; clear all; close all; clear global;
startup
ds = 'CrossOriRandDir_ExptList_DC';
eval(ds)

iexp = 4;

%% Load
mouse = expt(iexp).mouse;
date = expt(iexp).date;
ImgFolder = expt(iexp).coFolder;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

dataBase = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home', expt(iexp).saveLoc);
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\David';
sessionDir = fullfile(dataBase, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]);

load(fullfile(sessionDir, [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof', 'mask_cell')
load(fullfile(sessionDir, [date '_' mouse '_' run_str '_dataStim.mat']))

outDir = fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], 'dfofMaps');
if ~exist(outDir, 'dir'); mkdir(outDir); end
fprintf([mouse ' ' date ' — dF/F maps\n'])
fprintf('Saving to: %s\n', outDir)

%% Compute indexing
nStim1 = nStimCon * nMaskCon * nMaskPhas;
nPerDir = nMaskDiff + 1;
myfilter = fspecial('gaussian', [20 20], 0.5);

fprintf('nStim1=%d (contrast combos), nPerDir=%d, nStimDir=%d\n', nStim1, nPerDir, nStimDir)
fprintf('data_dfof has %d slices total\n', size(data_dfof, 3))

%% Grating alone maps (4x4)
nCol = ceil(sqrt(nStimDir));
nRow = ceil(nStimDir / nCol);

figure('Position', [100 100 900 900]);
for is = 1:nStimDir
    idx = nStim1 + (is-1)*nPerDir + 1;
    subplot(nRow, nCol, is)
    imagesc(imfilter(data_dfof(:,:,idx), myfilter)); axis image off
    title(sprintf('%.1f%s', stimDirs(is), char(176)), 'FontSize', 9)
end
sgtitle(sprintf('%s %s — Grating alone dF/F', mouse, date))
print(fullfile(outDir, [date '_' mouse '_' run_str '_dfofMaps_grating.png']), '-dpng', '-r300')

%% Plaid maps (4x4)
figure('Position', [1050 100 900 900]);
for is = 1:nStimDir
    idx = nStim1 + (is-1)*nPerDir + 2;
    subplot(nRow, nCol, is)
    imagesc(imfilter(data_dfof(:,:,idx), myfilter)); axis image off
    title(sprintf('%.1f%s', stimDirs(is), char(176)), 'FontSize', 9)
end
sgtitle(sprintf('%s %s — Plaid dF/F', mouse, date))
print(fullfile(outDir, [date '_' mouse '_' run_str '_dfofMaps_plaid.png']), '-dpng', '-r300')

%% Max projection + cell overlay
figure('Position', [100 50 600 600]);
data_dfof_filt = imfilter(data_dfof, myfilter);
data_dfof_max = max(data_dfof_filt(:,:,nStim1+1:end), [], 3);
imagesc(data_dfof_max); axis image off
hold on
boundaries = bwboundaries(mask_cell > 0);
for k = 1:length(boundaries)
    plot(boundaries{k}(:,2), boundaries{k}(:,1), 'w', 'LineWidth', 0.5)
end
hold off
sgtitle(sprintf('%s %s — Max dF/F + ROIs (%d cells)', mouse, date, max(mask_cell(:))))
print(fullfile(outDir, [date '_' mouse '_' run_str '_dfofMaps_maxROI.png']), '-dpng', '-r300')
