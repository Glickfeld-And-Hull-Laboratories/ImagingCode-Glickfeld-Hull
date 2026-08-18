function [fov_avg, fov_norm, fov_red, corrmap, dfmax, masks, maskNP, cellTCs_all, input] = ...
    correct_framerate_mismatch(fov_avg, fov_norm, fov_red, corrmap, dfmax, masks, maskNP, ...
    cellTCs_all, input, hz15p5_day, n_crop_lines, fn_backup)
% CORRECT_FRAMERATE_MISMATCH  Correct a day accidentally recorded at 15.5 Hz instead of 15 Hz.
%
% The 15.5 Hz recording produces a smaller FOV. This function:
%   1. Crops n_crop_lines from the top of the larger (15 Hz) day's images/masks
%   2. Resamples the 15.5 Hz day's timecourses to 15 Hz via linear interpolation
%   3. Rescales stimulus onset/offset frame indices for the 15.5 Hz day
%   4. Optionally backs up originals before overwriting
%
% NOTE: Pupil data (if extracted) must be corrected separately after extraction.
%
% Usage:
%   [fov_avg, fov_norm, fov_red, corrmap, dfmax, masks, maskNP, cellTCs_all, input] = ...
%       correct_framerate_mismatch(fov_avg, fov_norm, fov_red, corrmap, dfmax, masks, maskNP, ...
%       cellTCs_all, input, hz15p5_day, n_crop_lines)
%
%   [...] = correct_framerate_mismatch(..., fn_backup)
%
% Inputs:
%   fov_avg, fov_norm, fov_red - cell arrays {1,2} of FOV images
%   corrmap, dfmax             - cell arrays {1,2} of pixel correlation / max dF/F maps
%   masks                      - cell arrays {1,2} of labeled cell masks (from bwlabel)
%   maskNP                     - cell arrays {1,2} of neuropil masks [y x nCells] or [y x]
%   cellTCs_all                - cell arrays {1,2} of timecourses [nFrames x nCells]
%   input                      - struct array(1,2) with stimOns/Offs_photodiode fields
%   hz15p5_day                 - day index recorded at 15.5 Hz (1 or 2)
%   n_crop_lines               - lines to crop from top of the larger (15 Hz) FOV
%   fn_backup                  - (optional) folder path to save originals before modification

hz_normal = 15;
hz_fast   = 15.5;
crop_day  = 3 - hz15p5_day;   % the 15 Hz day has the larger FOV and gets cropped

% --- Backup originals ---
if nargin >= 12 && ~isempty(fn_backup)
    mkdir(fn_backup);
    orig_data.fov_avg      = fov_avg;
    orig_data.fov_norm     = fov_norm;
    orig_data.fov_red      = fov_red;
    orig_data.corrmap      = corrmap;
    orig_data.dfmax        = dfmax;
    orig_data.masks        = masks;
    orig_data.maskNP       = maskNP;
    orig_data.cellTCs_all  = cellTCs_all;
    orig_data.input        = input;
    save(fullfile(fn_backup, 'originals_before_framerate_correction.mat'), '-struct', 'orig_data');
    fprintf('Originals backed up to: %s\n', fn_backup);
end

% --- 1. Crop the larger (15 Hz) day FOV ---

% Check whether any cells sit entirely within the rows being cropped
cells_in_crop = unique(masks{crop_day}(1:n_crop_lines, :));
cells_in_crop(cells_in_crop == 0) = [];
if ~isempty(cells_in_crop)
    warning('Cropping will remove pixels from %d cell(s): [%s]. Check mask integrity.', ...
        numel(cells_in_crop), num2str(cells_in_crop));
end

crop_rows = n_crop_lines+1 : size(fov_avg{crop_day}, 1);

fov_avg{crop_day}  = fov_avg{crop_day}(crop_rows, :);
fov_norm{crop_day} = fov_norm{crop_day}(crop_rows, :);
fov_red{crop_day}  = fov_red{crop_day}(crop_rows, :);
corrmap{crop_day}  = corrmap{crop_day}(crop_rows, :);
dfmax{crop_day}    = dfmax{crop_day}(crop_rows, :);
masks{crop_day}    = masks{crop_day}(crop_rows, :);

if ndims(maskNP{crop_day}) == 3
    maskNP{crop_day} = maskNP{crop_day}(crop_rows, :, :);
else
    maskNP{crop_day} = maskNP{crop_day}(crop_rows, :);
end

fprintf('Cropped %d lines from top of day %d FOV. New size: [%d x %d]\n', ...
    n_crop_lines, crop_day, size(fov_avg{crop_day},1), size(fov_avg{crop_day},2));

% --- 2. Resample 15.5 Hz timecourses to 15 Hz ---

tc      = cellTCs_all{hz15p5_day};
nFrames = size(tc, 1);
original_t = (0:nFrames-1) / hz_fast;
new_t      = (0:nFrames-1) / hz_normal;

% new_t extends slightly beyond the original recording; clamp to valid range
new_t = min(new_t, original_t(end));

cellTCs_all{hz15p5_day} = interp1(original_t, tc, new_t, 'linear');

fprintf('Resampled day %d timecourses: %.1f Hz -> %.1f Hz (%d frames)\n', ...
    hz15p5_day, hz_fast, hz_normal, nFrames);

% --- 3. Correct stimulus timing for 15.5 Hz day ---
% Frame indices recorded at 15.5 Hz are rescaled to equivalent 15 Hz frame indices

scale = hz_normal / hz_fast;  % 15/15.5 ≈ 0.9677

input(hz15p5_day).stimOns_photodiode  = round(input(hz15p5_day).stimOns_photodiode  * scale);
input(hz15p5_day).stimOffs_photodiode = round(input(hz15p5_day).stimOffs_photodiode * scale);

fprintf('Corrected stimulus timing for day %d (scale = %.4f)\n', hz15p5_day, scale);
fprintf('Framerate mismatch correction complete.\n');
