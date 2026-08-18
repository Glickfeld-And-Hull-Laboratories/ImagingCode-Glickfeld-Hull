%%This script is following the OriAnalysisStepByStep tutorial on Google
%%Docs

%% Step 0: Load the data

baseDir = 'Z:\All_Staff\home\David\Analysis\Tutorial\221129_i2080\221129_i2080_runs-003';

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
resp_win = (nOff+1):(nOn/2)+nOff; 


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
    figure;
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
    figure;
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


%% Step 6: Summary Plots

figure;
bar([nCells, sum(responsive_cells), sum(reliable_cells)]);
xticklabels({'Segmented', 'Responsive', 'Reliably Fit'});
ylabel('Number of Cells');
title('Cell Summary');

figure;

subplot(1, 2, 1);
histogram(pref_ori(reliable_cells), uniqueOrientations);
xlabel('Preferred Orientation (deg)');
ylabel('Number of Cells');
title('Preferred Orientation Distribution');
xticks(uniqueOrientations);
xtickangle(45);

subplot(1, 2, 2);
histogram(k1_fit(reliable_cells));
xlabel('k1 (Sharpness)');
ylabel('Number of Cells');
title('Tuning Sharpness Distribution');