%%This script is following the OriAnalysisStepByStep tutorial on Google
%%Docs
clc; clear all; close all; clear all global;
startup
set(0, 'DefaultFigureVisible', 'off');
%% Step 0: Load the data

baseDir = 'Z:\All_Staff\home\David\Analysis\Tutorial\221129_i2080\221129_i2080_runs-003';
figDir = fullfile(baseDir, 'figures');
if ~exist(figDir, 'dir'); mkdir(figDir); end

expt = load(fullfile(baseDir, '221129_i2080_runs-003_input'));
data = load(fullfile(baseDir,'221129_i2080_runs-003_TCs.mat'));
nOn = 30;
nOff = 60;

nFramesPerTrial = nOff + nOn;
nTrials = size(data.npSub_tc, 1) / nFramesPerTrial;

allDirections =  cell2mat(expt.input.tGratingDirectionDeg);
allOrientations = mod(allDirections, 180);  % converts each trial's direction to orientation
uniqueDirections = unique(allDirections);
uniqueOrientations = uniqueDirections(1,(1:8)); 
nDirections = length(uniqueDirections);
nOrientations = nDirections/2; % each orientation has 2 directions



%% Step 1: transform timecourse matrix

npSub_trials = reshape(data.npSub_tc, nFramesPerTrial, nTrials, [] );
nCells = size(npSub_trials,3);


f = mean(npSub_trials(nOff/2:nOff,:,:),1);

df = npSub_trials - f;

dff = df ./ f;

%% Step 2: orientation tuning

base_win = nOff/2:nOff; 

% for simplicity, let's start out with hardcoded halfway through stim, 
% then we can find a tailored window that works by averaging the space
resp_win = (nOff+1):(nOff + nOn/2);


% take average of the cell responses during these time windows (frames)
% this means we are only looking at nTrials and nCells, so we transpose to
% get the right nCells x nTrials matrix
base = squeeze(mean(dff(base_win, : , :),1))';
resp = squeeze(mean(dff(resp_win,:,:), 1 ))';


alpha_corrected = 0.05 / (nOrientations - 1);  % Bonferroni correction

sig_mat  = zeros(nCells, nOrientations);
avg_resp = zeros(nCells, nOrientations);
sem_resp = zeros(nCells, nOrientations);

for i = 1:nOrientations
    trial_idx = allOrientations == uniqueOrientations(i);
    nTrialsOri = sum(trial_idx);

    base_ori = base(:, trial_idx);
    resp_ori = resp(:, trial_idx);
    
    avg_resp(:, i) = mean(resp_ori - base_ori, 2);
    sem_resp(:, i) = std(resp_ori - base_ori, 0, 2) ./ sqrt(nTrialsOri);
    
    for c = 1:nCells
        [h, ~] = ttest(resp_ori(c,:), base_ori(c,:), 'Alpha', alpha_corrected);
        sig_mat(c, i) = h;
    end
end



% check for responsiveness
responsive_cells = any(sig_mat, 2);
fprintf('Responsive (Bonferroni): %d / %d\n', sum(responsive_cells), nCells);

% Without Bonferroni for comparison
sig_mat_uncorrected = zeros(nCells, nOrientations);
for i = 1:nOrientations
    trial_idx = allOrientations == uniqueOrientations(i);
    base_ori = base(:, trial_idx);
    resp_ori = resp(:, trial_idx);
    for c = 1:nCells
        [h, ~] = ttest(resp_ori(c,:), base_ori(c,:), 'Alpha', 0.05);
        sig_mat_uncorrected(c, i) = h;
    end
end
responsive_uncorrected = any(sig_mat_uncorrected, 2);
fprintf('Responsive (uncorrected): %d / %d\n', sum(responsive_uncorrected), nCells);


%% Step 3: Plotting tuning curves


% avg_resp is nCells x nTrials

nPlotPerFig = 16;
nFigs = ceil(nCells / nPlotPerFig);

for tuning_figs = 1:nFigs
    figure('Visible','off', 'Position', [0 0 1920 1080]);
    for p = 1:nPlotPerFig
        cellIdx = (tuning_figs-1) * nPlotPerFig + p;
        if cellIdx > nCells
            break
        end
        subplot(4, 4, p);
        errorbar(uniqueOrientations, avg_resp(cellIdx,:), sem_resp(cellIdx,:), 'o-');
        title(sprintf('Cell %d | %d sig', cellIdx, sum(sig_mat(cellIdx,:))));
        xlabel('Orientation');
        ylabel('dF/F');
    end
    print(fullfile(figDir, ['TuningCurves_' num2str(tuning_figs) '.png']), '-dpng', '-r300');
    close;
end

%% Step 4: Fit data

resp_sub = resp - base; 


theta = deg2rad(uniqueOrientations);  % 1 x nOri, in radians

b_fit  = zeros(nCells, 1);
k1_fit = zeros(nCells, 1);
R_fit  = zeros(nCells, 1);
u1_fit = zeros(nCells, 1);
sse_fit = zeros(nCells, 1);
R2_fit  = zeros(nCells, 1);

% Compute average response per orientation from resp_sub
avg_resp_ori = zeros(nCells, nOrientations);
for i = 1:nOrientations
    trial_idx = allOrientations == uniqueOrientations(i);
    avg_resp_ori(:, i) = mean(resp_sub(:, trial_idx), 2);
end

% Fit each cell
for c = 1:nCells
    [b_fit(c), k1_fit(c), R_fit(c), u1_fit(c), sse_fit(c), R2_fit(c)] = ...
        miaovonmisesfit_ori(theta, avg_resp_ori(c,:));
end



theta_fine = deg2rad(0:1:179);
nPlotPerFig = 16;
nFigs = ceil(nCells / nPlotPerFig);

for fig = 1:nFigs
    figure('Visible','off', 'Position', [0 0 1920 1080]);
    for p = 1:nPlotPerFig
        cellIdx = (fig-1) * nPlotPerFig + p;
        if cellIdx > nCells
            break
        end
        subplot(4, 4, p);
        errorbar(uniqueOrientations, avg_resp_ori(cellIdx,:), sem_resp(cellIdx,:), 'o-');
        hold on;
        y_fit = b_fit(cellIdx) + R_fit(cellIdx) * exp(k1_fit(cellIdx) * (cos(2*(theta_fine - u1_fit(cellIdx))) - 1));
        plot(rad2deg(theta_fine), y_fit, 'r-');
        hold off;
        title(sprintf('Cell %d | %d sig', cellIdx, sum(sig_mat(cellIdx,:))));
        xlabel('Orientation');
        xticks(uniqueOrientations);
        xtickangle(45);  % angle the labels so they don't overlap
        ylabel('dF/F');
    end
    print(fullfile(figDir, ['VonMisesFit_' num2str(fig) '.png']), '-dpng', '-r300');
    close;
end
pref_ori = rad2deg(u1_fit);

%% Step 5: Model Verification

nBoots = 1000;
u1_boot = zeros(nCells, nBoots);  % store preferred ori from each bootstrap

for b = 1:nBoots
    % For each orientation, resample trials with replacement
    avg_resp_boot = zeros(nCells, nOrientations);
    
    for i = 1:nOrientations
        trial_idx = find(allOrientations == uniqueOrientations(i));  % indices of trials for this ori
        resampled_idx = randsample(trial_idx, length(trial_idx), true);  % resample with replacement
        avg_resp_boot(:, i) = mean(resp_sub(:, resampled_idx), 2);
    end
    
    % Fit each cell with the resampled tuning curve
    for c = 1:nCells
        [~, ~, ~, u1_tmp, ~, ~] = miaovonmisesfit_ori(theta, avg_resp_boot(c,:));
        u1_boot(c, b) = rad2deg(u1_tmp);
    end
end

% reliablity verification
% Difference between each bootstrap's pref ori and the original fit
pref_ori = rad2deg(u1_fit);  % original preferred orientations from Step 4
boot_diff = abs(u1_boot - pref_ori);  % nCells x nBoots

% Sort and find 90th percentile for each cell
boot_diff_sorted = sort(boot_diff, 2);  % sort along bootstrap dimension
pct90 = boot_diff_sorted(:, 900);  % 900th out of 1000 = 90th percentile

% Reliably fit cells: 90% of bootstraps within one orientation step
reliable_cells = pct90 < 22.5;
fprintf('Reliably fit: %d / %d\n', sum(reliable_cells), nCells);

% Intersection: reliably fit AND responsive
reliable_responsive = reliable_cells & responsive_cells;
fprintf('Reliably fit & responsive: %d / %d\n', sum(reliable_responsive), nCells);


%% Step 6: Summary Plots

figure('Visible','off', 'Position', [0 0 1920 1080]);
bar([nCells, sum(responsive_cells), sum(reliable_cells), sum(reliable_responsive)]);
xticklabels({'Segmented', 'Responsive', 'Reliably Fit', 'Reliable & Responsive'});
ylabel('Number of Cells');
title('Cell Summary');
print(fullfile(figDir, 'CellSummary.png'), '-dpng', '-r300');
close;

figure('Visible','off', 'Position', [0 0 1920 1080]);

subplot(1, 2, 1);
% Wrap pref_ori into [0, 180) to handle circular orientation space
pref_ori_wrapped = mod(pref_ori(reliable_responsive), 180);
% Bin edges centered on each orientation: half-step before first, half-step after last
ori_step = uniqueOrientations(2) - uniqueOrientations(1);  % 22.5
ori_edges = [uniqueOrientations - ori_step/2, uniqueOrientations(end) + ori_step/2];
histogram(pref_ori_wrapped, ori_edges);
xlabel('Preferred Orientation (deg)');
ylabel('Number of Cells');
title('Preferred Orientation Distribution');
xticks(uniqueOrientations);
xtickangle(45);

subplot(1, 2, 2);
k1_edges = [0:0.5:5, 6:2:30];  % fine bins 0-5, coarser 6-30
histogram(k1_fit(reliable_responsive), k1_edges);
xlabel('k1 (Sharpness)');
ylabel('Number of Cells');
title('Tuning Sharpness Distribution');
print(fullfile(figDir, 'SummaryHistograms.png'), '-dpng', '-r300');
close;

%% Step 7: Inspect k1 ≈ 0 cells

% Find reliably-fit & responsive cells with k1 in the lowest histogram bin
k1_thresh = 1;  % adjust to match the bin edge you want to inspect
k1_zero_idx = find(reliable_responsive & k1_fit < k1_thresh);
fprintf('Reliable+responsive with k1 < %g: %d\n', k1_thresh, length(k1_zero_idx));

nPlotPerFig = 16;
nFigs_k1 = ceil(length(k1_zero_idx) / nPlotPerFig);

for fig = 1:nFigs_k1
    figure('Visible','off', 'Position', [0 0 1920 1080]);
    for p = 1:nPlotPerFig
        idx = (fig-1) * nPlotPerFig + p;
        if idx > length(k1_zero_idx)
            break
        end
        cellIdx = k1_zero_idx(idx);
        subplot(4, 4, p);
        errorbar(uniqueOrientations, avg_resp_ori(cellIdx,:), sem_resp(cellIdx,:), 'o-');
        hold on;
        y_fit = b_fit(cellIdx) + R_fit(cellIdx) * exp(k1_fit(cellIdx) * (cos(2*(theta_fine - u1_fit(cellIdx))) - 1));
        plot(rad2deg(theta_fine), y_fit, 'r-');
        hold off;
        title(sprintf('Cell %d | k1=%.2f | R=%.3f | R2=%.2f', cellIdx, k1_fit(cellIdx), R_fit(cellIdx), R2_fit(cellIdx)));
        xlabel('Orientation');
        xticks(uniqueOrientations);
        xtickangle(45);
        ylabel('dF/F');
    end
    sgtitle(sprintf('k1 < %g cells (%d/%d)', k1_thresh, fig, nFigs_k1));
    print(fullfile(figDir, ['k1Zero_TuningCurves_' num2str(fig) '.png']), '-dpng', '-r300');
    close;
end