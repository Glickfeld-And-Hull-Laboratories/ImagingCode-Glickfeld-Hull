function validate_framerate_correction(orig, corrected, hz15p5_day, n_crop_lines)
% VALIDATE_FRAMERATE_CORRECTION  Sanity checks after correct_framerate_mismatch.
%
% Usage:
%   validate_framerate_correction(orig, corr, hz15p5_day, n_crop_lines)
%
% Inputs:
%   orig         - struct with original data fields (loaded from backup mat):
%                    .fov_avg, .masks, .cellTCs_all, .input
%   corr         - struct with corrected data (current workspace variables):
%                    .fov_avg, .masks, .cellTCs_all, .input
%   hz15p5_day   - which day was recorded at 15.5 Hz (1 or 2)
%   n_crop_lines - number of lines cropped from the other day
%
% Example - run right after correct_framerate_mismatch:
%   orig = load(fullfile(fn_backup, 'originals_before_framerate_correction.mat'));
%   corr.fov_avg = fov_avg; corr.masks = masks;
%   corr.cellTCs_all = cellTCs_all; corr.input = input;
%   validate_framerate_correction(orig, corr, 1, 17);

crop_day  = 3 - hz15p5_day;
hz_normal = 15;
hz_fast   = 15.5;
scale     = hz_normal / hz_fast;
pass      = true;

fprintf('\n=== Framerate correction validation ===\n');
fprintf('15.5 Hz day: %d  |  Cropped day: %d  |  Crop lines: %d\n\n', ...
    hz15p5_day, crop_day, n_crop_lines);

% 1. FOV sizes match after cropping
sz1 = size(corrected.fov_avg{1});
sz2 = size(corrected.fov_avg{2});
if isequal(sz1, sz2)
    fprintf('[PASS] FOV sizes match: [%d x %d]\n', sz1(1), sz1(2));
else
    fprintf('[FAIL] FOV size mismatch: day1=[%d x %d], day2=[%d x %d]\n', ...
        sz1(1), sz1(2), sz2(1), sz2(2));
    pass = false;
end

% Check crop removed the right number of lines
orig_sz = size(orig.fov_avg{crop_day});
new_sz  = size(corrected.fov_avg{crop_day});
expected_rows = orig_sz(1) - n_crop_lines;
if new_sz(1) == expected_rows
    fprintf('[PASS] Crop: %d rows removed from day %d (%d -> %d)\n', ...
        n_crop_lines, crop_day, orig_sz(1), new_sz(1));
else
    fprintf('[FAIL] Crop row count wrong: expected %d, got %d\n', expected_rows, new_sz(1));
    pass = false;
end

% 2. Cell count unchanged after crop
nc_orig = max(orig.masks{crop_day}(:));
nc_corr = max(corrected.masks{crop_day}(:));
if nc_orig == nc_corr
    fprintf('[PASS] Cell count unchanged after crop: %d cells\n', nc_corr);
else
    fprintf('[FAIL] Cell count changed: %d -> %d (cells severed by crop)\n', nc_orig, nc_corr);
    pass = false;
end

% 3. No NaNs in resampled timecourses
tc_corr = corrected.cellTCs_all{hz15p5_day};
n_nan   = sum(isnan(tc_corr(:)));
if n_nan == 0
    fprintf('[PASS] No NaNs in resampled timecourses\n');
else
    fprintf('[FAIL] %d NaN values in resampled timecourses\n', n_nan);
    pass = false;
end

% 4. Stim onsets/offsets within valid frame range
nFrames = size(tc_corr, 1);
ons  = corrected.input(hz15p5_day).stimOns_photodiode;
offs = corrected.input(hz15p5_day).stimOffs_photodiode;
if all(ons >= 1) && all(ons <= nFrames) && all(offs >= 1) && all(offs <= nFrames)
    fprintf('[PASS] Stim onsets/offsets within valid frame range [1 %d]\n', nFrames);
else
    fprintf('[FAIL] Stim times out of range: ons [%d %d], offs [%d %d], nFrames=%d\n', ...
        min(ons), max(ons), min(offs), max(offs), nFrames);
    pass = false;
end

% 5. Inter-stimulus intervals scale correctly
orig_ons  = orig.input(hz15p5_day).stimOns_photodiode;
corr_ons  = corrected.input(hz15p5_day).stimOns_photodiode;
if numel(orig_ons) > 1
    orig_isis = diff(double(orig_ons));
    corr_isis = diff(double(corr_ons));
    actual_scale   = mean(corr_isis) / mean(orig_isis);
    expected_scale = scale;
    if abs(actual_scale - expected_scale) < 0.01
        fprintf('[PASS] ISI scaling correct: %.4f (expected %.4f)\n', actual_scale, expected_scale);
    else
        fprintf('[FAIL] ISI scaling off: got %.4f, expected %.4f\n', actual_scale, expected_scale);
        pass = false;
    end
else
    fprintf('[SKIP] Only one stimulus onset - cannot check ISI scaling\n');
end

% 6. Timecourse correlation: original vs resampled (mean across cells)
tc_orig = orig.cellTCs_all{hz15p5_day};
% Both have nFrames rows; compare directly (same length by construction)
r_vals = diag(corrected(tc_orig, tc_corr));   % correlation per cell
med_r  = median(r_vals);
if med_r > 0.99
    fprintf('[PASS] Median TC correlation (orig vs resampled): %.4f\n', med_r);
elseif med_r > 0.95
    fprintf('[WARN] Median TC correlation lower than expected: %.4f\n', med_r);
else
    fprintf('[FAIL] Median TC correlation very low: %.4f\n', med_r);
    pass = false;
end

fprintf('\n');
if pass
    fprintf('All checks passed.\n');
else
    fprintf('One or more checks failed - review above.\n');
end
fprintf('=========================================\n\n');
