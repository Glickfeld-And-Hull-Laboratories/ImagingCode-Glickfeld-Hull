% test_ret_position.m
% Retinotopy analysis using masks and FOV from a separate reference run.
% Registers ret data to the reference FOV, applies reference masks,
% extracts timecourses, fits 2D Gaussian RFs, and finds optimal stimulus position.

clear all; clc; close all;

%% Parameters
mouse = 'i2244'
date = '260604'
time = '1241'
RetImgFolder = '005'
frame_rate = 15

refRun  = '006';
refDate = date;

if computer == 'GLNXA64'
    isilonName = '/home/cc735@dhe.duke.edu/GlickfeldLabShare/All_staff';
else
    isilonName = 'Z:';
end
base       = fullfile(isilonName, '/home/ACh/Data/2p_data/');
fnOut_base = fullfile(isilonName, '/home/ACh/Analysis/2p_analysis/VIP_YM90K');

run_str = ['runs-' RetImgFolder];
ref_str = ['runs-' refRun];
fnOut   = fullfile(fnOut_base, mouse, date, RetImgFolder);
mkdir(fnOut)

%% Load retinotopy imaging data
CD = fullfile(base, mouse, date, RetImgFolder);
cd(CD);
imgMatFile = [RetImgFolder '_000_000.mat'];
load(imgMatFile);
nframes = info.config.frames;
fprintf('Reading %d frames...\n', nframes)
data_temp = squeeze(sbxread([RetImgFolder '_000_000'], 0, nframes));

%% Load behavior data
fName = fullfile(isilonName, '/Behavior/Data', ['data-' mouse '-' date '-' time '.mat']);
load(fName);
retino_input = input;
clear input
nOn  = retino_input.nScansOn;
nOff = retino_input.nScansOff;

%% Find stimulus onsets via photodiode
fprintf('Finding stimulus onsets via photodiode...\n')
[stimOns, ~] = photoFrameFinder_Sanworks(info.frame);
ntrials = length(stimOns);
fprintf('%d trials detected\n', ntrials)

%% Load reference FOV and masks
fprintf('Loading reference FOV and masks from run %s...\n', refRun)
load(fullfile(fnOut_base, mouse, refDate, refRun, 'regOuts&Img.mat'));
load(fullfile(fnOut_base, mouse, refDate, refRun, 'mask_cell.mat'));   % mask_cell, mask_np, mask_label
load(fullfile(fnOut_base, mouse, refDate, refRun, 'resp_inds.mat'));   % resp (logical)
load(fullfile(fnOut_base, mouse, refDate, refRun, 'input.mat'));
ref_input = input;
resp_ind = find(resp);

referenceFOV = data_avg;
nCells = max(mask_cell(:));
fprintf('%d cells loaded from reference run\n', nCells)

%% Register ret data to reference FOV
fprintf('Registering ret data to reference FOV...\n')
retAvrg = mean(data_temp, 3);
[~, retAvrg_reg] = stackRegGPU(retAvrg, referenceFOV);
[~, data_reg]    = stackRegGPU(data_temp, retAvrg_reg);
sz = size(data_reg);
clear data_temp retAvrg retAvrg_reg

figure;
imagesc(mean(data_reg, 3)); colormap gray; caxis([200 4000])
hold on
bounds = bwboundaries(mask_cell > 0);
for i = 1:length(bounds)
    plot(bounds{i}(:,2), bounds{i}(:,1), 'r', 'LineWidth', 0.5)
end
title(sprintf('Registered FOV  %d imported masks', nCells))
set(gca, 'TickDir', 'out', 'XTick', [], 'YTick', []); box off

%% Extract timecourses with neuropil subtraction
fprintf('Extracting timecourses...\n')
down          = 5;
data_reg_down = stackGroupProject(data_reg, down);
data_tc       = stackGetTimeCourses(data_reg, mask_cell);
data_tc_down  = stackGetTimeCourses(data_reg_down, mask_cell);

np_tc      = zeros(sz(3), nCells);
np_tc_down = zeros(floor(sz(3)/down), nCells);
for i = 1:nCells
    np_tc(:,i)      = stackGetTimeCourses(data_reg, mask_np(:,:,i));
    np_tc_down(:,i) = stackGetTimeCourses(data_reg_down, mask_np(:,:,i));
    fprintf('  Cell %d / %d\n', i, nCells)
end
clear data_reg_down

ii = 0.01:0.01:1;
x  = zeros(100, nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down - tcRemoveDC(np_tc_down * ii(i)));
end
[~, ind] = max(x, [], 1);
np_w = 0.01 * ind;

npSub_tc = data_tc - tcRemoveDC(np_tc) .* np_w;

%% Build trial matrix
fprintf('Building trial matrix from PD onsets...\n')
tc_mat = nan(nOn+nOff, nCells, ntrials);
for itrial = 1:ntrials
    on = stimOns(itrial);
    if ~isnan(on) && (on + nOn + nOff/2) <= size(npSub_tc, 1)
        tc_mat(:,:,itrial) = npSub_tc((on - nOff/2):(on - 1 + nOn + nOff/2), :);
    end
end

baselineFrames = (nOff/4+1):nOff/2;
stimFrames     = (nOff/2+1):(nOff/2+nOn);

tc_f    = mean(tc_mat(baselineFrames, :, :), 1);
tc_dfof = (tc_mat - tc_f) ./ tc_f;
clear tc_mat tc_f

%% Stimulus info
Az  = celleqel2mat_padded(retino_input.tGratingAzimuthDeg);
El  = celleqel2mat_padded(retino_input.tGratingElevationDeg);
Azs = unique(Az);
Els = unique(El);
if min(Els) < 0, Els = fliplr(Els); end

nStim = length(Azs) * length(Els);
Stims = zeros(nStim, 2);
idx = 1;
for iEl = 1:length(Els)
    for iAz = 1:length(Azs)
        Stims(idx,:) = [Els(iEl) Azs(iAz)];
        idx = idx + 1;
    end
end
fprintf('%d unique Az+El stimuli\n', nStim)

%% Response per Az/El combination
ret_mean_resp = squeeze(mean(tc_dfof(stimFrames, :, :), 1));

resp_by_stim = nan(length(Els), length(Azs), nCells);
for i_el = 1:length(Els)
    indE = find(El == Els(i_el));
    for i_az = 1:length(Azs)
        indA = find(Az == Azs(i_az));
        these_trials = intersect(indE, indA);
        resp_by_stim(i_el, i_az, :) = nanmean(ret_mean_resp(:, these_trials), 2);
    end
end

%% Build Ind_struct
Ind_struct = [];
for iStim = 1:nStim
    indE = find(El == Stims(iStim,1));
    indA = find(Az == Stims(iStim,2));
    Ind_struct(iStim).all_trials = intersect(indE, indA);
end

%% Fit retinotopy: 2D ellipse fits with shuffling
fprintf('Begin fitting retinotopy data...\n')

resp_dFoverF = squeeze(mean(tc_dfof(stimFrames, :, :), 1));

[AzAz, ElEl] = meshgrid(Azs, Els);
grid2.AzAz   = AzAz;
grid2.ElEl   = ElEl;
dAz = median(diff(Azs));
dEl = median(diff(Els));
Az_vec00 = Azs(1):(dAz/10):Azs(end);
El_vec00 = Els(1):(dEl/10):Els(end);
[AzAz00, ElEl00]  = meshgrid(Az_vec00, El_vec00);
grid2.AzAz00 = AzAz00;
grid2.ElEl00 = ElEl00;

Nshuf  = 100;
h_all  = zeros(1, nCells);
h_all(resp_ind) = 1;
Fit_struct  = [];
ret_run_str = run_str;

fprintf('Nshuf = %d\n', Nshuf)
for count_shuf = 0:Nshuf
    fprintf('count_shuf: %d/%d\n', count_shuf, Nshuf)
    Im_mat_USE = zeros(nCells, nStim);
    for iCond = 1:nStim
        ind_all = Ind_struct(iCond).all_trials;
        if count_shuf > 0
            ind_all_1 = ind_all(randsample(length(ind_all), length(ind_all), 1));
        else
            ind_all_1 = ind_all;
        end
        Im_mat_USE(:,iCond) = mean(resp_dFoverF(:,ind_all_1), 2);
    end

    for iCell = 1:nCells
        if count_shuf > 0 && h_all(1,iCell) == 0
            continue
        end
        a = Im_mat_USE(iCell,:);
        if max(a,[],2) > 0
            b    = reshape(a', length(Azs), length(Els));
            data = b';
            PLOTIT_FIT  = 0;
            SAVEALLDATA = (count_shuf == 0);
            Fit_2Dellipse_ret_lbub
            if count_shuf == 0
                Fit_struct(iCell).True.s_ = s;
            else
                Fit_struct(iCell).Shuf(count_shuf).s_ = s;
            end
        end
    end
end
fprintf('Shuffling done\n')
save(fullfile(fnOut, [date '_' mouse '_' run_str '_Fit_struct.mat']), 'Fit_struct', '-v7.3')

%% Goodness-of-fit assessment (lbub)
fprintf('Assessing goodness of fit\n')
fit_true_vec = NaN(nCells, 10);
for iCell = 1:nCells
    if ~isempty(Fit_struct(iCell).True)
        tmp = [Fit_struct(iCell).True.s_.x, ...
               Fit_struct(iCell).True.s_.Elhicut_50, ...
               Fit_struct(iCell).True.s_.Azhicut_50, ...
               Fit_struct(iCell).True.s_.Elhicut_10, ...
               Fit_struct(iCell).True.s_.Azhicut_10];
        fit_true_vec(iCell,:) = tmp;
    end
end

fit_shuf_vec = NaN(nCells, 10, Nshuf);
for count_shuf = 1:Nshuf
    for iCell = 1:nCells
        if ~isempty(Fit_struct(iCell).Shuf)
            tmp = [Fit_struct(iCell).Shuf(count_shuf).s_.x, ...
                   Fit_struct(iCell).Shuf(count_shuf).s_.Elhicut_50, ...
                   Fit_struct(iCell).Shuf(count_shuf).s_.Azhicut_50, ...
                   Fit_struct(iCell).Shuf(count_shuf).s_.Elhicut_10, ...
                   Fit_struct(iCell).Shuf(count_shuf).s_.Azhicut_10];
            fit_shuf_vec(iCell,:,count_shuf) = tmp;
        end
    end
end

Npars       = size(fit_shuf_vec, 2);
lbub_fits   = NaN(nCells, Npars, 5);
alpha_bound = 0.025;
ind_shuf_lb = ceil(Nshuf * alpha_bound);
ind_shuf_ub = ceil(Nshuf * (1-alpha_bound));
for iCell = 1:nCells
    for count2 = 1:Npars
        tmp = squeeze(fit_shuf_vec(iCell,count2,:));
        i_sorted = sort(tmp);
        lbub_fits(iCell,count2,1) = i_sorted(ind_shuf_lb);
        lbub_fits(iCell,count2,2) = i_sorted(ind_shuf_ub);
        lbub_fits(iCell,count2,3) = mean(i_sorted);
        lbub_fits(iCell,count2,5) = std(i_sorted);
    end
    lbub_fits(iCell,:,4) = fit_true_vec(iCell,:);
end

lbub_diff = lbub_fits(:,:,2) - lbub_fits(:,:,1);

goodfit_ind = [];
for iCell = 1:nCells
    if lbub_diff(iCell,4) < retino_input.gratingAzimuthStepDeg*2 && ...
       lbub_diff(iCell,5) < retino_input.gratingAzimuthStepDeg*2
        goodfit_ind = [goodfit_ind iCell];
    end
end

goodfit_ind2 = zeros(size(goodfit_ind));
for i = 1:length(goodfit_ind)
    if sum(round(lbub_fits(goodfit_ind(i),4,4)) == [min(Azs) max(Azs)]) || ...
       sum(round(lbub_fits(goodfit_ind(i),5,4)) == [min(Els) max(Els)])
        continue
    end
    goodfit_ind2(i) = goodfit_ind(i);
end
goodfit_ind = goodfit_ind2(goodfit_ind2 ~= 0);
clear goodfit_ind2
fprintf('%d good-fit cells (final, lbub)\n', length(goodfit_ind))

fitAzimDeg = lbub_fits(:, 4, 4)';
fitElevDeg = lbub_fits(:, 5, 4)';

%% Optimal stimulus position (good-fit cells: all, then labeled only)
% plotRFdistanceMap uses col 5 = Az, col 4 = El from fit_true_vec
% All good-fit cells
rfAz_all = fit_true_vec(goodfit_ind, 5);
rfEl_all = fit_true_vec(goodfit_ind, 4);
opt_Az_all = mean(rfAz_all);
opt_El_all = mean(rfEl_all);
fprintf('Optimal position (all good-fit, n=%d): Az=%.2f deg, El=%.2f deg\n', ...
    length(goodfit_ind), opt_Az_all, opt_El_all)

% Labeled (red) cells among good-fit
goodfit_labeled_ind = goodfit_ind(mask_label(goodfit_ind) == 1);
rfAz_lab = fit_true_vec(goodfit_labeled_ind, 5);
rfEl_lab = fit_true_vec(goodfit_labeled_ind, 4);
opt_Az_lab = mean(rfAz_lab);
opt_El_lab = mean(rfEl_lab);
fprintf('Optimal position (labeled good-fit, n=%d): Az=%.2f deg, El=%.2f deg\n', ...
    length(goodfit_labeled_ind), opt_Az_lab, opt_El_lab)

%% RF distance maps: actual stim position vs optimal position
actual_Az = double(ref_input.gratingAzimuthDeg);
actual_El = double(ref_input.gratingElevationDeg);

[distMap_actual, dist_vec_actual] = plotRFdistanceMap( ...
    lbub_fits(:,:,4), goodfit_ind, mask_cell, actual_Az, actual_El);
title(sprintf('Actual stim position (Az=%.1f, El=%.1f)\nmedian dist=%.2f deg (n=%d goodfit)', ...
    actual_Az, actual_El, median(dist_vec_actual, 'omitnan'), length(goodfit_ind)))

[distMap_opt, dist_vec_opt] = plotRFdistanceMap( ...
    lbub_fits(:,:,4), goodfit_ind, mask_cell, opt_Az_all, opt_El_all);
title(sprintf('Optimal position (Az=%.1f, El=%.1f)\nmedian dist=%.2f deg (n=%d goodfit)', ...
    opt_Az_all, opt_El_all, median(dist_vec_opt, 'omitnan'), length(goodfit_ind)))

%% Save
save(fullfile(fnOut, [date '_' mouse '_' run_str '_lbub_fits.mat']), ...
    'lbub_fits', 'lbub_diff', 'goodfit_ind', 'resp_ind', 'fitAzimDeg', 'fitElevDeg', ...
    'opt_Az_all', 'opt_El_all', 'opt_Az_lab', 'opt_El_lab', ...
    'goodfit_labeled_ind', '-v7.3')

save(fullfile(fnOut, [date '_' mouse '_' run_str '_ret_fromRef.mat']), ...
    'resp_by_stim', 'npSub_tc', 'tc_dfof', 'resp_ind', ...
    'goodfit_ind', 'Azs', 'Els', 'Stims', 'fitAzimDeg', 'fitElevDeg', '-v7.3')
fprintf('Saved to %s\n', fnOut)

%% Save all figures
figHandles = findall(0, 'Type', 'figure');
for i = 1:length(figHandles)
    fig = figHandles(i);
    figName = get(fig, 'Name');
    if isempty(figName)
        figName = sprintf('fig%d', fig.Number);
    else
        figName = sprintf('fig%d_%s', fig.Number, strrep(figName, ' ', '_'));
    end
    print(fig, fullfile(fnOut, [date '_' mouse '_' run_str '_' figName '.pdf']), '-dpdf', '-bestfit')
end
fprintf('Saved %d figures\n', length(figHandles))