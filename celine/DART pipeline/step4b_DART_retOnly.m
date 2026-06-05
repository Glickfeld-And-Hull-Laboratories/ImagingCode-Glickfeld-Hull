% ret_alignment_summary.m
% Loads the multiday matched alignment, fits retinotopy for each day, and reports:
%   - mean/median RF distance from actual stimulus location
%   - optimal stimulus position (Az/El that minimizes mean population RF distance)
%   - mean/median RF distance from that optimal position

clearvars -except instr sessNum
if ~exist('instr', 'var')
    instr = input('Enter name of instructions file: ', 's');
end
run(instr);
if exist('sessNum', 'var'), instructions.session = num2str(sessNum); end
ds = instructions.ds;
run(ds);
rc = behavConstsDART;

day_id           = str2double(instructions.session);
match_day        = expt(day_id).multiday_matchdays;
nd               = 2;
mouse            = expt(day_id).mouse;
ExperimentFolder = expt(day_id).exptType;

if expt(day_id).multiday_timesincedrug_hours > 0
    dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end
fn_multi = fullfile(rc.analysis, ExperimentFolder, mouse, ['multiday_' dart_str]);

switch instructions.refDay
    case '1'; allDays = [match_day, day_id];
    case '2'; allDays = [day_id, match_day];
end
fprintf('Sessions: %s\n', num2str(allDays));
clear instr

load(fullfile(fn_multi, 'timecourses.mat'))          % match_ind lives here
load(fullfile(fn_multi, 'multiday_alignment.mat'))  % fov_avg, masks, fitGeoTAf, match_ind
load(fullfile(fn_multi, 'input.mat'))
inputStructure = input; clear input

[ret_npSub_tc_matched, ret_distance_matched, resp_by_stim_matched, ret_dfof_trial_matched, ...
 trialIndSourceUsed, lbub_fits_matched, goodfit_ind_matched, r2_vec_matched, ...
 fitAzimDeg_matched, fitElevDeg_matched, prefAzimDeg_matched, prefElevDeg_matched, ...
 Azs_matched, Els_matched, distMap_matched, dist_vec_matched] = ...
    retinotopy_for_matched_data(nd, allDays, expt, mouse, fov_avg, masks, fitGeoTAf, ...
        instructions, inputStructure, match_ind, false, fn_multi);

fprintf('\n=== RETINOTOPY ALIGNMENT SUMMARY ===\n')
for id = 1:nd
    goodfit = goodfit_ind_matched{id};

    % Extract RF centers for goodfit cells from true fit params
    % lbub_fits dim3=4 is the true (unshuffled) fit; col4=Az0, col5=El0
    fit_true = lbub_fits_matched{id}(:,:,4);
    rfAz = fit_true(goodfit, 4);
    rfEl = fit_true(goodfit, 5);
    valid = ~isnan(rfAz) & ~isnan(rfEl);
    rfAz = rfAz(valid);
    rfEl = rfEl(valid);

    stimAz = double(inputStructure(id).gratingAzimuthDeg);
    stimEl = double(inputStructure(id).gratingElevationDeg);
    dist_act = sqrt((rfAz - stimAz).^2 + (rfEl - stimEl).^2);

    % Geometric median via Weiszfeld algorithm: minimizes mean Euclidean distance,
    % guaranteeing mean(dist_opt) <= mean(dist_act)
    p = [mean(rfAz); mean(rfEl)];
    for iter = 1:200
        d = max(sqrt((rfAz - p(1)).^2 + (rfEl - p(2)).^2), 1e-10);
        w = 1 ./ d;
        p_new = [sum(w .* rfAz); sum(w .* rfEl)] / sum(w);
        if norm(p_new - p) < 1e-8, break; end
        p = p_new;
    end
    opt_Az = p(1);  opt_El = p(2);
    dist_opt = sqrt((rfAz - opt_Az).^2 + (rfEl - opt_El).^2);

    fprintf('\nDay %d (session %d):\n', id, allDays(id))
    fprintf('  Actual stim:               Az = %.1f deg,  El = %.1f deg\n', stimAz, stimEl)
    fprintf('  Good-fit cells:            n = %d\n', sum(valid))
    fprintf('  RF dist from actual stim:  mean = %.2f deg,  median = %.2f deg\n', mean(dist_act), median(dist_act))
    fprintf('  Optimal stim position:     Az = %.1f deg,  El = %.1f deg\n', opt_Az, opt_El)
    fprintf('  RF dist from optimal stim: mean = %.2f deg,  median = %.2f deg\n', mean(dist_opt), median(dist_opt))
end
fprintf('\n')

save(fullfile(fn_multi, 'retino_aligned.mat'), ...
    'ret_npSub_tc_matched', 'ret_distance_matched', 'resp_by_stim_matched', 'ret_dfof_trial_matched', ...
    'trialIndSourceUsed', 'goodfit_ind_matched', 'r2_vec_matched', 'lbub_fits_matched', ...
    'fitAzimDeg_matched', 'fitElevDeg_matched', 'prefAzimDeg_matched', 'prefElevDeg_matched', ...
    'Azs_matched', 'Els_matched', 'distMap_matched', 'dist_vec_matched');
fprintf('Saved retino_aligned.mat to %s\n', fn_multi);

load(fullfile(fn_multi, 'cell_analysis.mat'), 'keep_cells');

ret_distance_keep = cell(1, nd);
for id = 1:nd
    ret_distance_keep{id} = ret_distance_matched{id}(:, keep_cells);
end

goodfit_both     = intersect(goodfit_ind_matched{1}, goodfit_ind_matched{2});
goodfit_ind_keep = find(ismember(keep_cells, goodfit_both));

save(fullfile(fn_multi, 'retino_aligned.mat'), ...
    'ret_npSub_tc_matched', 'ret_distance_matched', 'resp_by_stim_matched', 'ret_dfof_trial_matched', ...
    'trialIndSourceUsed', 'goodfit_ind_matched', 'goodfit_ind_keep', 'r2_vec_matched', 'lbub_fits_matched', ...
    'fitAzimDeg_matched', 'fitElevDeg_matched', 'prefAzimDeg_matched', 'prefElevDeg_matched', ...
    'Azs_matched', 'Els_matched', 'distMap_matched', 'dist_vec_matched', ...
    'ret_distance_keep');
fprintf('Saved retino_aligned.mat to %s\n', fn_multi);