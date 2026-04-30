%standaloe script to get match-aligned retinotopy data without redoing the
%extraction
clear all; clear global; clc
prompt = 'Enter name of instructions file: ';
instr = input(prompt, 's');
clear prompt
run(instr);

ds=instructions.ds;
run(ds); 

dataStructLabels = {'contrastxori'};
rc = behavConstsDART; % directories

% Input validation
if ~exist('expt', 'var')
    error('Dataset %s not found', ds);
end

day_id = str2double(instructions.session);


if length(expt) < day_id
    error('Day_id %d not valid for this dataset', day_id);
else
    match_day = expt(day_id).multiday_matchdays;

end

nd = 2; % hardcoding the number of days for now
mouse = expt(day_id).mouse;
experimentFolder = expt(day_id).exptType;
fnout = fullfile(rc.analysis, mouse);

% Set up file naming based on drug condition
if expt(day_id).multiday_timesincedrug_hours > 0
    dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end
fn_multi = fullfile(rc.analysis, experimentFolder, mouse, ['multiday_' dart_str]);

% Determine which session was used as reference for cell matching
x = instructions.refDay;
switch x
    case '1'
        pre = 1;  % baseline session, used as reference, is in the 1st position
        post = 2;
        fprintf('Baseline used as reference\n');
         allDays = [match_day,day_id];
        fprintf('Analyzing sessions: %s\n', num2str(allDays));
    case '2'
        pre = 2;
        post = 1;  % post-DART session, used as reference, is in the 1st position
        fprintf('Post-DART used as reference\n');
        allDays = [day_id, match_day];
        fprintf('Analyzing sessions: %s\n', num2str(allDays));
end
clear x instr

% Load the matched cell data
cd(fn_multi)
load(fullfile(fn_multi, 'timecourses.mat'))
load(fullfile(fn_multi, 'multiday_alignment.mat'))
load(fullfile(fn_multi, 'input.mat'))
load(fullfile(fn_multi,'cell_analysis.mat'))

% Rename input to inputStructure to avoid conflict with MATLAB's input function
inputStructure = input;
clear input
%%
 [ret_npSub_tc_matched, ret_distance_matched, resp_by_stim_matched, ret_dfof_trial_matched, trialIndSourceUsed, ...
     lbub_fits_matched, goodfit_ind_matched, r2_vec_matched, fitAzimDeg_matched, fitElevDeg_matched, ...
     prefAzimDeg_matched, prefElevDeg_matched, Azs_matched, Els_matched, distMap_matched, dist_vec_matched] = ...
        retinotopy_for_matched_data(nd, allDays, expt, mouse, fov_avg, masks, fitGeoTAf, ...
            instructions, inputStructure, match_ind, false, fn_multi);

    ret_npSub_tc_keep   = cell(1, nd);
    ret_distance_keep   = cell(1, nd);
    resp_by_stim_keep   = cell(1, nd);
    ret_dfof_trial_keep = cell(1, nd);
    for id = 1:nd
        ret_npSub_tc_keep{id}   = ret_npSub_tc_matched{id}(:, keep_cells);
        ret_distance_keep{id}   = ret_distance_matched{id}(:, keep_cells);
        resp_by_stim_keep{id}   = resp_by_stim_matched{id}(:, :, keep_cells);
        ret_dfof_trial_keep{id} = ret_dfof_trial_matched{id}(:, keep_cells, :);
    end
    % Convert goodfit indices from matched set to keep_cells subset
    goodfit_both = intersect(goodfit_ind_matched{1},goodfit_ind_matched{2});
    goodfit_ind_keep = find(ismember(keep_cells, goodfit_both));

    % Subset per-cell fit results to keep_cells
    for id = 1:nd
        r2_vec_keep{id}       = r2_vec_matched{id}(keep_cells);
        fitAzimDeg_keep{id}   = fitAzimDeg_matched{id}(keep_cells);
        fitElevDeg_keep{id}   = fitElevDeg_matched{id}(keep_cells);
        prefAzimDeg_keep{id}  = prefAzimDeg_matched{id}(keep_cells);
        prefElevDeg_keep{id}  = prefElevDeg_matched{id}(keep_cells);
        lbub_fits_keep{id}    = lbub_fits_matched{id}(keep_cells, :, :);
    end

    save(fullfile(fn_multi, 'retino_aligned.mat'), ...
        'ret_npSub_tc_matched', 'ret_distance_matched', 'resp_by_stim_matched', 'ret_dfof_trial_matched', ...
        'ret_npSub_tc_keep', 'ret_distance_keep', 'resp_by_stim_keep', 'ret_dfof_trial_keep', 'trialIndSourceUsed', ...
        'goodfit_ind_matched', 'goodfit_ind_keep', ...
        'r2_vec_matched', 'r2_vec_keep', 'lbub_fits_matched', 'lbub_fits_keep', ...
        'fitAzimDeg_matched', 'fitAzimDeg_keep', 'fitElevDeg_matched', 'fitElevDeg_keep', ...
        'prefAzimDeg_matched', 'prefAzimDeg_keep', 'prefElevDeg_matched', 'prefElevDeg_keep', ...
        'distMap_matched', 'dist_vec_matched');

    % Visual space figure: 4 panels (day x HTP group)
    % lbub_fits_matched{id}: [nMatch x 10 x 5], dim 3 slice 4 = true fit values
    %   col 4 = Az center, col 5 = El center, col 7 = Elhicut_50, col 8 = Azhicut_50
    stimAz    = arrayfun(@(id) double(inputStructure(id).gratingAzimuthDeg), 1:nd);
    stimEl    = arrayfun(@(id) double(inputStructure(id).gratingElevationDeg), 1:nd);
    theta     = linspace(0, 2*pi, 100);
    htpGroups = {find(red_cells_keep), find(green_cells_keep)};
    htpLabels = {'HTP+', 'HTP-'};
    dayColors = {'b', 'k'};
    dayLabels = {'Day 1', 'Day 2'};

    figure;
    for id = 1:nd
        lbub = lbub_fits_matched{id};
        col  = dayColors{id};
        for ig = 1:2
            subplot(nd, 2, (id-1)*2 + ig); hold on;
            plot(stimAz(id), stimEl(id), 'x', 'Color', col, 'MarkerSize', 14, 'LineWidth', 2);
            cellInd = htpGroups{ig};
            for i = 1:length(cellInd)
                mc    = keep_cells(cellInd(i));
                Az0   = lbub(mc, 4, 4);
                El0   = lbub(mc, 5, 4);
                sigAz = lbub(mc, 8, 4);
                sigEl = lbub(mc, 7, 4);
                if any(isnan([Az0 El0 sigAz sigEl])), continue; end
                xe = Az0 + sigAz * cos(theta);
                ye = El0 + sigEl * sin(theta);
                if ismember(mc, goodfit_ind_matched{id})
                    plot(xe, ye, col, 'LineWidth', 0.5);
                else
                    plot(xe, ye, col, 'LineWidth', 0.5, 'LineStyle', '--');
                end
            end
            xlim([-50 50]); ylim([-50 50]);
            xlabel('Azimuth (deg)')
            ylabel('Elevation (deg)')
            title(sprintf('%s  %s', dayLabels{id}, htpLabels{ig}))
            axis square
            set(gca, 'TickDir', 'out', 'box', 'off')
        end
    end