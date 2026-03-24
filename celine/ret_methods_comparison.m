% retOnly_fromRef.m
% Retinotopy analysis using masks and FOV from a separate reference run.
% Registers ret data to the reference FOV, applies reference masks,
% extracts timecourses, and filters to cells responsive in the reference run.

clear all; clc; close all;

%% Parameters � ret run
mouse        = 'i1428';
date         = '260320';
time         = '1506';
RetImgFolder = '002';
frame_rate   = 15;

% Reference run (provides FOV, masks, and responsive cell list)
refRun  = '003';
refDate = date;   % change if reference is from a different date

if computer == 'GLNXA64'
    isilonName = '/home/cc735@dhe.duke.edu/GlickfeldLabShare/All_staff';
else
    isilonName = 'Z:';
end
base       = fullfile(isilonName, '/home/ACh/Data/2p_data/');
fnOut_base = fullfile(isilonName, '/home/ACh/Analysis/2p_analysis/retTesting');

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
% 'retino_input' is now in workspace
retino_input=input;
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
load(fullfile(fnOut_base, mouse, refDate, refRun,'regOuts&Img.mat'));   % loads data_avg
load(fullfile(fnOut_base, mouse, refDate, refRun, 'mask_cell.mat'));    % loads mask_cell, mask_np
load(fullfile(fnOut_base, mouse, refDate, refRun, 'resp_inds.mat'));  % loads resp (logical
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

% Validation: registered FOV with imported masks overlaid
figure;
imagesc(mean(data_reg, 3)); colormap gray; caxis([200 4000])
hold on
bounds = bwboundaries(mask_cell > 0);
for i = 1:length(bounds)
    plot(bounds{i}(:,2), bounds{i}(:,1), 'r', 'LineWidth', 0.5)
end
title(sprintf('Registered FOV � %d imported masks', nCells))
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

% Neuropil weights by maximizing skewness
ii = 0.01:0.01:1;
x  = zeros(100, nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down - tcRemoveDC(np_tc_down * ii(i)));
end
[~, ind] = max(x, [], 1);
np_w = 0.01 * ind;

npSub_tc = data_tc - tcRemoveDC(np_tc) .* np_w;
clear data_reg data_reg_down data_tc_down np_tc_down

%% Build trial matrix using PD onsets
% Window: nOff/2 frames before onset, nOn+nOff/2 frames after
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

% Build Stims matrix (El, Az) matching data_dfof_avg stacking order
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
ret_mean_resp = squeeze(mean(tc_dfof(stimFrames, :, :), 1));  % nCells x ntrials

resp_by_stim = nan(length(Els), length(Azs), nCells);
for i_el = 1:length(Els)
    indE = find(El == Els(i_el));
    for i_az = 1:length(Azs)
        indA = find(Az == Azs(i_az));
        these_trials = intersect(indE, indA);
        resp_by_stim(i_el, i_az, :) = nanmean(ret_mean_resp(:, these_trials), 2);
    end
end

%% Weighted average preferred position
[azGrid, elGrid] = meshgrid(Azs, Els);
prefAzimDeg = nan(1, nCells);
prefElevDeg = nan(1, nCells);
for iCell = 1:nCells
    w = max(resp_by_stim(:,:,iCell), 0);
    total_w = sum(w(:));
    if total_w > 0
        prefAzimDeg(iCell) = sum(w(:) .* azGrid(:)) / total_w;
        prefElevDeg(iCell) = sum(w(:) .* elGrid(:)) / total_w;
    end
end
fprintf('%d/%d cells have NaN preferred position (all responses <= 0)\n', ...
    sum(isnan(prefAzimDeg)), nCells)


%% Build Ind_struct (trial indices per stimulus, needed for fitting)
Ind_struct = [];
for iStim = 1:nStim
    indE = find(El == Stims(iStim,1));
    indA = find(Az == Stims(iStim,2));
    Ind_struct(iStim).all_trials = intersect(indE, indA);
end

%% Fit retinotopy � 2D ellipse fits with shuffling
fprintf('Begin fitting retinotopy data...\n')

resp_dFoverF = squeeze(mean(tc_dfof(stimFrames, :, :), 1));  % nCells x ntrials

[AzAz, ElEl]   = meshgrid(Azs, Els);
grid2.AzAz     = AzAz;
grid2.ElEl     = ElEl;
dAz = median(diff(Azs));
dEl = median(diff(Els));
Az_vec00 = Azs(1):(dAz/10):Azs(end);
El_vec00 = Els(1):(dEl/10):Els(end);
[AzAz00, ElEl00]  = meshgrid(Az_vec00, El_vec00);
grid2.AzAz00 = AzAz00;
grid2.ElEl00 = ElEl00;

Nshuf   = 100;
h_all   = zeros(1, nCells);
h_all(resp_ind) = 1;
Fit_struct = [];
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

    ifig  = 1;
    start = 1;
    for iCell = 1:nCells
        if count_shuf > 0 && h_all(1,iCell) == 0
            continue
        end
        a = Im_mat_USE(iCell,:);
        if max(a,[],2) > 0
            b    = reshape(a', length(Azs), length(Els));
            data = b';
            if count_shuf == 0
                PLOTIT_FIT  = 1;
                SAVEALLDATA = 1;
                Fit_2Dellipse_ret_lbub
                eval(['Fit_struct(iCell).True.s_ = s;']);
            else
                PLOTIT_FIT  = 0;
                SAVEALLDATA = 0;
                Fit_2Dellipse_ret_CC
                eval(['Fit_struct(iCell).Shuf(count_shuf).s_ = s;']);
            end
        end
    end
    if count_shuf == 0
        Im_mat_true = Im_mat_USE;  % save true (non-shuffled) responses for R2
        set(gcf, 'Position', [0 0 800 1000]);
        print(fullfile(fnOut, [date '_' mouse '_' run_str '_RFfits' num2str(ifig) '.pdf']), '-dpdf')
    end
end
fprintf('Shuffling done\n')
save(fullfile(fnOut, [date '_' mouse '_' run_str '_Fit_struct.mat']), 'Fit_struct', '-v7.3')

%% Goodness-of-fit assessment
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

Npars         = size(fit_shuf_vec, 2);
lbub_fits     = NaN(nCells, Npars, 5);
alpha_bound   = 0.025;
ind_shuf_lb   = ceil(Nshuf * alpha_bound);
ind_shuf_ub   = ceil(Nshuf * (1-alpha_bound));
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

% Remove RFs at retinotopy perimeter
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

%% R-squared goodness-of-fit (alternative to lbub)
r2_threshold = 0.65;
r2_vec = nan(1, nCells);
for iCell = 1:nCells
    if isempty(Fit_struct(iCell).True)
        continue
    end
    x_fit = Fit_struct(iCell).True.s_.x;
    % x_fit = [A, sigma_Az, sigma_El, Az0, El0, xi]
    A       = x_fit(1);
    sig_az  = x_fit(2);
    sig_el  = x_fit(3);
    Az0     = x_fit(4);
    El0     = x_fit(5);
    xi      = x_fit(6);
    az_rot  =  (azGrid - Az0) .* cos(xi) + (elGrid - El0) .* sin(xi);
    el_rot  = -(azGrid - Az0) .* sin(xi) + (elGrid - El0) .* cos(xi);
    fit_surf = A .* exp(-(az_rot.^2 ./ (2*sig_az^2) + el_rot.^2 ./ (2*sig_el^2)));
    actual   = resp_by_stim(:,:,iCell);  % nEls x nAzs, matches azGrid/elGrid
    ss_res   = sum((actual(:) - fit_surf(:)).^2);
    ss_tot   = sum((actual(:) - mean(actual(:))).^2);
    if ss_tot > 0
        r2_vec(iCell) = 1 - ss_res/ss_tot;
    end
end

goodfit_ind_r2 = find(r2_vec >= r2_threshold);
fprintf('%d good-fit cells (R2 >= %.2f)\n', length(goodfit_ind_r2), r2_threshold)

% Select example cells: ~6 from goodfit_ind, ~3 non-goodfit responsive cells
rng('shuffle');
n_good_ex    = min(6, length(goodfit_ind));
non_good_ind = resp_ind(~ismember(resp_ind, goodfit_ind));
n_nongood_ex = min(3, length(non_good_ind));
ex_cells = [goodfit_ind(randperm(length(goodfit_ind), n_good_ex)), ...
            non_good_ind(randperm(length(non_good_ind), n_nongood_ex))'];
n_ex = length(ex_cells);

% Weighted average preferred position
figure('Name', 'Weighted average preferred position');
for i = 1:n_ex
    subplot(3, 3, i)
    imagesc(Azs, Els, resp_by_stim(:,:,ex_cells(i)))
    set(gca, 'YDir', 'normal', 'TickDir', 'out'); box off; grid off
    hold on
    plot(prefAzimDeg(ex_cells(i)), prefElevDeg(ex_cells(i)), 'r+', 'MarkerSize', 10, 'LineWidth', 2)
    is_good = ismember(ex_cells(i), goodfit_ind);
    title(sprintf('Cell %d%s', ex_cells(i), repmat('*', 1, is_good)))
    if i > 6, xlabel('Az (deg)'); end
    if mod(i,3) == 1, ylabel('El (deg)'); end
end
colormap parula

% 2D Gaussian fit (same cells)
figure('Name', '2D Gaussian fit');
for i = 1:n_ex
    subplot(3, 3, i)
    imagesc(Azs, Els, resp_by_stim(:,:,ex_cells(i)))
    set(gca, 'YDir', 'normal', 'TickDir', 'out'); box off; grid off
    hold on
    if ~isempty(Fit_struct(ex_cells(i)).True)
        x_fit_i  = Fit_struct(ex_cells(i)).True.s_.x;
        A_i      = x_fit_i(1);  sig_az_i = x_fit_i(2);  sig_el_i = x_fit_i(3);
        Az0_i    = x_fit_i(4);  El0_i    = x_fit_i(5);  xi_i     = x_fit_i(6);
        az_rot_i =  (AzAz00 - Az0_i) .* cos(xi_i) + (ElEl00 - El0_i) .* sin(xi_i);
        el_rot_i = -(AzAz00 - Az0_i) .* sin(xi_i) + (ElEl00 - El0_i) .* cos(xi_i);
        fit_surf_i = A_i .* exp(-(az_rot_i.^2 ./ (2*sig_az_i^2) + el_rot_i.^2 ./ (2*sig_el_i^2)));
        contour(Az_vec00, El_vec00, fit_surf_i, 3, 'w')
        plot(Az0_i, El0_i, 'w+', 'MarkerSize', 10, 'LineWidth', 2)
        r2_str = sprintf(' R�=%.2f', r2_vec(ex_cells(i)));
    else
        r2_str = ' no fit';
    end
    is_good = ismember(ex_cells(i), goodfit_ind);
    title(sprintf('Cell %d%s%s', ex_cells(i), repmat('*', 1, is_good), r2_str))
    if i > 6, xlabel('Az (deg)'); end
    if mod(i,3) == 1, ylabel('El (deg)'); end
end
colormap parula

% Extract fit-based preferred location from lbub_fits (col 4=Az0, col 5=El0, dim3=4 is true fit)
fitAzimDeg = lbub_fits(:, 4, 4)';  % 1 x nCells
fitElevDeg = lbub_fits(:, 5, 4)';

% Distance between fit center and weighted average, for goodfit cells
rf_dist = sqrt((fitAzimDeg(goodfit_ind) - prefAzimDeg(goodfit_ind)).^2 + ...
               (fitElevDeg(goodfit_ind) - prefElevDeg(goodfit_ind)).^2);

figure;
histogram(rf_dist, 'BinWidth', 2, 'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'none')
xlabel('Distance between fit center and weighted average (deg)')
ylabel('# cells')
title(sprintf('RF center method comparison (n=%d goodfit cells)', length(goodfit_ind)))
set(gca, 'TickDir', 'out'); box off; grid off
xline(median(rf_dist, 'omitnan'), 'k--', sprintf('median=%.1f�', median(rf_dist, 'omitnan')), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom')


save(fullfile(fnOut, [date '_' mouse '_' run_str '_lbub_fits.mat']), ...
    'lbub_fits', 'lbub_diff', 'goodfit_ind', 'goodfit_ind_r2', 'r2_vec', 'resp_ind', 'fitAzimDeg', 'fitElevDeg', '-v7.3')

%% Filter outputs to responsive cells
resp_by_stim_resp = resp_by_stim(:, :, resp_ind);
prefAzimDeg_resp  = prefAzimDeg(resp_ind);
prefElevDeg_resp  = prefElevDeg(resp_ind);
npSub_tc_resp     = npSub_tc(:, resp_ind);
tc_dfof_resp      = tc_dfof(:, resp_ind, :);
fprintf('%d/%d cells retained (responsive)\n', length(resp_ind), nCells)

%% Save
save(fullfile(fnOut, [date '_' mouse '_' run_str '_ret_fromRef.mat']), ...
    'resp_by_stim', 'resp_by_stim_resp', ...
    'prefAzimDeg',  'prefElevDeg', ...
    'prefAzimDeg_resp', 'prefElevDeg_resp', ...
    'npSub_tc', 'npSub_tc_resp', ...
    'tc_dfof',  'tc_dfof_resp', ...
    'resp_ind', 'goodfit_ind', 'goodfit_ind_r2', 'r2_vec', 'Azs', 'Els', 'Stims', ...
    'fitAzimDeg', 'fitElevDeg', '-v7.3')
fprintf('Saved to %s\n', fnOut)

%% Venn diagram: lbub vs R2 good-fit cells
lbub_set = ismember(resp_ind, goodfit_ind);
r2_set   = ismember(resp_ind, goodfit_ind_r2);

n_lbub_only = sum( lbub_set & ~r2_set);
n_r2_only   = sum(~lbub_set &  r2_set);
n_both      = sum( lbub_set &  r2_set);
n_neither   = sum(~lbub_set & ~r2_set);

fprintf('lbub only: %d | R2 only: %d | both: %d | neither: %d\n', ...
    n_lbub_only, n_r2_only, n_both, n_neither)

figure;
pie([n_lbub_only, n_r2_only, n_both, n_neither], ...
    {sprintf('lbub only (n=%d)', n_lbub_only), ...
     sprintf('R� only (n=%d)',   n_r2_only), ...
     sprintf('both (n=%d)',      n_both), ...
     sprintf('neither (n=%d)',   n_neither)})
colormap([0.2 0.5 0.8; 0.9 0.4 0.3; 0.5 0.3 0.7; 0.8 0.8 0.8])
title(sprintf('Good-fit cells (of %d responsive)', length(resp_ind)))

%% Nested circle diagram: total -> responsive -> goodfit
figure;
theta = linspace(0, 2*pi, 300);

% Scale radii by sqrt(n) so area is proportional to count
r_total = 1;
r_resp  = sqrt(length(resp_ind) / nCells);
r_good  = sqrt(length(goodfit_ind) / nCells);

patch(r_total*cos(theta), r_total*sin(theta), [0.85 0.85 0.85], 'EdgeColor', 'none'); hold on
patch(r_resp *cos(theta), r_resp *sin(theta), [0.4  0.6  0.85], 'EdgeColor', 'none')
patch(r_good *cos(theta), r_good *sin(theta), [0.2  0.45 0.7 ], 'EdgeColor', 'none')

text(0,  r_total-0.05, sprintf('all cells\nn=%d', nCells),           'HorizontalAlignment', 'center', 'VerticalAlignment', 'top',    'FontSize', 10)
text(0,  r_resp-0.05,  sprintf('responsive\nn=%d', length(resp_ind)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top',    'FontSize', 10, 'Color', 'w')
text(0, -0.05,         sprintf('goodfit\nn=%d', length(goodfit_ind)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top',    'FontSize', 10, 'Color', 'w')

axis equal off
title('Cell counts by criterion')

%% Compare RF centers from fits to selected stim location
plotRFdistanceMap(lbub_fits(:,:,4), goodfit_ind, mask_cell, double(ref_input.gratingAzimuthDeg), double(ref_input.gratingElevationDeg));