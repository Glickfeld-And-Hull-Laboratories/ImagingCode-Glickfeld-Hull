% ret_optimal_position.m
% Load ret data, fit RFs, find optimal stimulus location, plot heatmap.

clear all; clc; close all;

%% Parameters
mouse        = 'i2238';
date         = '260415';
time         = '1055';
RetImgFolder = '004';
refRun       = '006';
refDate      = date;

if computer == 'GLNXA64'
    isilonName = '/home/cc735@dhe.duke.edu/GlickfeldLabShare/All_staff';
else
    isilonName = 'Z:';
end
base       = fullfile(isilonName, '/home/ACh/Data/2p_data/');
fnOut_base = fullfile(isilonName, '/home/ACh/Analysis/2p_analysis/VIP_YM90K');
fnOut      = fullfile(fnOut_base, mouse, date, RetImgFolder);
mkdir(fnOut)
run_str    = ['runs-' RetImgFolder];

%% Load imaging data
CD = fullfile(base, mouse, date, RetImgFolder);
cd(CD);
load([RetImgFolder '_000_000.mat']);
nframes   = info.config.frames;
data_temp = squeeze(sbxread([RetImgFolder '_000_000'], 0, nframes));

%% Load behavior data
load(fullfile(isilonName, '/Behavior/Data', ['data-' mouse '-' date '-' time '.mat']));
retino_input = input; clear input
nOn  = retino_input.nScansOn;
nOff = retino_input.nScansOff;

%% Photodiode onsets
[stimOns, ~] = photoFrameFinder_Sanworks(info.frame);
ntrials = length(stimOns);
fprintf('%d trials detected\n', ntrials)

%% Load reference FOV and masks
load(fullfile(fnOut_base, mouse, refDate, refRun, 'regOuts&Img.mat'));
load(fullfile(fnOut_base, mouse, refDate, refRun, 'mask_cell.mat'));
load(fullfile(fnOut_base, mouse, refDate, refRun, 'resp_inds.mat'));
resp_ind     = find(resp);
referenceFOV = data_avg;
nCells       = max(mask_cell(:));

%% Register
retAvrg          = mean(data_temp, 3);
[~, retAvrg_reg] = stackRegGPU(retAvrg, referenceFOV);
[~, data_reg]    = stackRegGPU(data_temp, retAvrg_reg);
sz = size(data_reg);
clear data_temp retAvrg retAvrg_reg

%% Extract timecourses with neuropil subtraction
down         = 5;
data_tc      = stackGetTimeCourses(data_reg, mask_cell);
data_tc_down = stackGetTimeCourses(stackGroupProject(data_reg, down), mask_cell);

np_tc      = zeros(sz(3), nCells);
np_tc_down = zeros(floor(sz(3)/down), nCells);
for i = 1:nCells
    np_tc(:,i)      = stackGetTimeCourses(data_reg, mask_np(:,:,i));
    np_tc_down(:,i) = stackGetTimeCourses(stackGroupProject(data_reg, down), mask_np(:,:,i));
end

ii = 0.01:0.01:1;
x  = zeros(100, nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down - tcRemoveDC(np_tc_down * ii(i)));
end
[~, ind] = max(x, [], 1);
npSub_tc = data_tc - tcRemoveDC(np_tc) .* (0.01 * ind);
clear data_reg data_tc data_tc_down np_tc np_tc_down

%% Trial matrix and dF/F
tc_mat = nan(nOn+nOff, nCells, ntrials);
for itrial = 1:ntrials
    on = stimOns(itrial);
    if ~isnan(on) && (on - nOff/2 >= 1) && (on + nOn + nOff/2 - 1) <= size(npSub_tc,1)
        tc_mat(:,:,itrial) = npSub_tc((on - nOff/2):(on - 1 + nOn + nOff/2), :);
    end
end
baselineFrames = (nOff/4+1):nOff/2;
stimFrames     = (nOff/2+1):(nOff/2+nOn);
tc_f           = mean(tc_mat(baselineFrames, :, :), 1);
tc_dfof        = (tc_mat - tc_f) ./ tc_f;
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

%% Response per Az/El
resp_dFoverF = squeeze(mean(tc_dfof(stimFrames, :, :), 1));  % nCells x ntrials

Ind_struct = [];
for iStim = 1:nStim
    Ind_struct(iStim).all_trials = intersect(find(El == Stims(iStim,1)), find(Az == Stims(iStim,2)));
end

%% Fit RFs (no shuffling)
fprintf('Fitting RFs...\n')
[AzAz, ElEl] = meshgrid(Azs, Els);
grid2.AzAz   = AzAz;
grid2.ElEl   = ElEl;
dAz = median(diff(Azs));
dEl = median(diff(Els));
[AzAz00, ElEl00] = meshgrid(Azs(1):(dAz/10):Azs(end), Els(1):(dEl/10):Els(end));
grid2.AzAz00 = AzAz00;
grid2.ElEl00 = ElEl00;

h_all = zeros(1, nCells);
h_all(resp_ind) = 1;
Fit_struct  = [];
ret_run_str = run_str;
count_shuf  = 0;

Im_mat_USE = zeros(nCells, nStim);
for iCond = 1:nStim
    Im_mat_USE(:,iCond) = mean(resp_dFoverF(:, Ind_struct(iCond).all_trials), 2);
end

for iCell = 1:nCells
    a = Im_mat_USE(iCell,:);
    if max(a) > 0
        data        = reshape(a', length(Azs), length(Els))';
        PLOTIT_FIT  = 0;
        SAVEALLDATA = 1;
        Fit_2Dellipse_ret_lbub
        Fit_struct(iCell).True.s_ = s;
    end
end

%% Collect fit centers
fit_true_vec = NaN(nCells, 10);
for iCell = 1:nCells
    if ~isempty(Fit_struct(iCell).True)
        fit_true_vec(iCell,:) = [Fit_struct(iCell).True.s_.x, ...
                                  Fit_struct(iCell).True.s_.Elhicut_50, ...
                                  Fit_struct(iCell).True.s_.Azhicut_50, ...
                                  Fit_struct(iCell).True.s_.Elhicut_10, ...
                                  Fit_struct(iCell).True.s_.Azhicut_10];
    end
end

%% Optimal stimulus position
all_fit_ind = find(~isnan(fit_true_vec(:,4)));
[opt_Az, opt_El] = findOptimalStimLocation(fit_true_vec, all_fit_ind);

%% Pixel heatmap with optimal position
pix_dfof_avg = zeros(sz(1), sz(2), nStim);
idx = 1;
for i_el = 1:length(Els)
    for i_az = 1:length(Azs)
        these = find(El == Els(i_el) & Az == Azs(i_az));
        imgs  = zeros(sz(1), sz(2), length(these));
        for k = 1:length(these)
            on = stimOns(these(k));
            if ~isnan(on) && (on - nOff/2 >= 1) && (on + nOn - 1 <= sz(3))
                baseF = mean(data_reg(:,:, on-nOff/2 : on-1), 3);
                stimF = mean(data_reg(:,:, on : on+nOn-1), 3);
                imgs(:,:,k) = (double(stimF) - double(baseF)) ./ double(baseF);
            end
        end
        pix_dfof_avg(:,:,idx) = mean(imgs, 3);
        idx = idx + 1;
    end
end

pixThresh        = 0.2 * max(pix_dfof_avg(:));
img_avg_resp_pix = zeros(1, nStim);
for i = 1:nStim
    img = pix_dfof_avg(:,:,i);
    img_avg_resp_pix(i) = mean(img(img > pixThresh));
end

respMatrix = fliplr(rot90(reshape(img_avg_resp_pix, length(Els), length(Azs)), 3));
figure;
imagesc(Azs, flip(Els), respMatrix);
set(gca, 'YDir', 'normal', 'TickDir', 'out'); box off; grid off
hold on
plot(opt_Az, opt_El, 'rx', 'MarkerSize', 12, 'LineWidth', 2)
xlabel('Azimuth (deg)'); ylabel('Elevation (deg)')
title(sprintf('Optimal stim: Az=%.1f, El=%.1f deg', opt_Az, opt_El))
colorbar
