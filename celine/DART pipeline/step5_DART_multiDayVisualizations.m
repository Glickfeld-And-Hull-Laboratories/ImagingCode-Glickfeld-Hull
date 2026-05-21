% This script has the basic steps to investigate matched 2-photon data such
% as from DART experiments. It focusses on visualization, rather than
% statistical analyis.

clear all; close all; clc

prompt = 'Enter name of instructions file: ';
instr = input(prompt, 's');
clear prompt
run(instr);


ds = instructions.ds;
run(ds);

dataStructLabels = {'contrastxori'};
experimentFolder = instructions.experimentFolder;
rc = behavConstsDART;

sess_list = instructions.sess_list;
nSess = length(sess_list);
nd = 2;
frame_rate = instructions.frame_rate;
targetCon = instructions.targetCon;
nCon = length(targetCon);
targetSize = instructions.targetSize;
nSize = length(targetSize);
retDistThresh = 7.5; % degrees - RF distance threshold; cells must be < this on both days

% Determine which session was used as reference for cell matching
x = instructions.refDay;
switch x
    case '1'
        pre = 1;
        post = 2;
        fprintf('Baseline used as reference\n');
        matchDrx = 'f';
    case '2'
        pre = 2;
        post = 1;
        fprintf('Post-DART used as reference\n');
        matchDrx = 'r';
end
clear x instr

% Create output directory
sess_title = num2str(sess_list(1));
for iSess = 2:nSess
    sess_title = [sess_title '_' num2str(sess_list(iSess))];
end
d = char(string(datetime('today')));

if nSess == 1
    if expt(sess_list(1)).multiday_timesincedrug_hours > 0
        dart_str = [expt(sess_list(1)).drug '_' num2str(expt(sess_list(1)).multiday_timesincedrug_hours) 'Hr'];
    else
        dart_str = 'control';
    end
    fnout = fullfile(rc.analysis, experimentFolder, expt(sess_list(1)).mouse, ['multiday_' dart_str], d);
else
    fnout = fullfile(rc.analysis, experimentFolder, ['concat' sess_title], d);
end
mkdir(fnout); cd(fnout)

% Initialize concatenation variables
mice = []; red_concat = []; green_concat = []; nKeep_concat = [];
dirs_concat = []; cons_concat = []; sizes_concat = [];

% Initialize cell arrays for concatenation
cell_vars = {'tc_trial_avrg_stat', 'tc_trial_avrg_loc', 'tc_trial_avrg_stat_largePupil', 'tc_trial_avrg_stat_smallPupil', ...
    'pref_responses_stat', 'pref_responses_loc', 'pref_responses_stat_smallPupil', 'pref_responses_stat_largePupil', 'h', 'data_resp', 'RIx', 'norm_dir_resp', 'pref_dir', ...
    'noiseCorr', 'sigCorr', 'nonPref_trial_avrg_stat', 'nonPref_trial_avrg_loc', 'data_dfof_runOnset'};
for i = 1:length(cell_vars)
    eval([cell_vars{i} '_concat = cell(1,nd);']);
end
pref_allTrials_stat_concat = cell(nCon, nSize, nd);
pref_allTrials_loc_concat  = cell(nCon, nSize, nd);

% Initialize retinotopy concatenation variables if needed
if isfield(instructions, 'load_retino') && instructions.load_retino
    ret_npSub_tc_keep_concat   = cell(1, nd);
    ret_distance_keep_concat   = cell(1, nd);
    resp_by_stim_keep_concat   = cell(1, nd);
    ret_dfof_trial_keep_concat = cell(1, nd);
end

% Retinotopy filter variables (goodfit_ind_keep is intersection of both days,
% already in keep-cell space  written by step 4 to retino_aligned.mat)
goodfit_concat             = cell(1, nd);
ret_distance_retino_concat = cell(1, nd);
for id = 1:nd
    goodfit_concat{id}             = [];
    ret_distance_retino_concat{id} = [];
end

% Cell count tracking across sessions
nMatched_concat     = [];
nMatched_red_concat = [];
nMatched_grn_concat = [];

drug = cell(1, nSess);

% Main concatenation loop
for iSess = 1:nSess
    thisSess = sess_list(iSess);
    mouse = expt(thisSess).mouse;
    mice = [mice; mouse];
    drug{iSess} = expt(thisSess).drug;

    % Determine data path
    if expt(thisSess).multiday_timesincedrug_hours > 0
        dart_str = [expt(thisSess).drug '_' num2str(expt(thisSess).multiday_timesincedrug_hours) 'Hr'];
    else
        dart_str = 'control';
    end
    fn_multi = fullfile(rc.analysis, experimentFolder, mouse, ['multiday_' dart_str]);

    % Load data
    load(fullfile(fn_multi, 'tc_keep.mat'));
    load(fullfile(fn_multi, 'resp_keep.mat'));
    load(fullfile(fn_multi, 'input.mat'));
    load(fullfile(fn_multi, 'behavioral_state.mat'));
    load(fullfile(fn_multi, 'cell_analysis.mat'));
    load(fullfile(fn_multi, 'HT_pyr_relationship.mat'));
    nKeep = size(tc_trial_avrg_stat{post}, 2);

    % Load retinotopy data subsetted to keep_cells (written by step 4)
    % goodfit_ind_keep is already the intersection of both days in keep-cell space
    retino_file = fullfile(fn_multi, 'retino_aligned.mat');
    if exist(retino_file, 'file')
        retData = load(retino_file, 'goodfit_ind_keep', 'ret_distance_keep');
        gf_this = false(1, nKeep);
        gf_this(retData.goodfit_ind_keep) = true;
        for id = 1:nd
            goodfit_concat{id}             = [goodfit_concat{id}, gf_this];
            ret_distance_retino_concat{id} = [ret_distance_retino_concat{id}, retData.ret_distance_keep{id}];
        end
    else
        warning('retino_aligned.mat not found for session %d: %s', iSess, retino_file);
        for id = 1:nd
            goodfit_concat{id}             = [goodfit_concat{id}, false(1, nKeep)];
            ret_distance_retino_concat{id} = [ret_distance_retino_concat{id}, nan(1, nKeep)];
        end
    end

    % Track matched cell counts from cell_summary (loaded with cell_analysis.mat)
    nMatched_concat     = [nMatched_concat,     cell_summary.Total_Cells(1)];
    nMatched_red_concat = [nMatched_red_concat, cell_summary.Total_Cells(2)];
    nMatched_grn_concat = [nMatched_grn_concat, cell_summary.Total_Cells(3)];

    inputStructure = input;
    clear input

    % Load optional per-cell retinotopy data if requested
    if isfield(instructions, 'load_retino') && instructions.load_retino
        if exist(retino_file, 'file')
            load(retino_file);
            fprintf('Loaded retinotopy data for session %d\n', iSess);
        else
            warning('retino_aligned.mat not found for session %d', iSess);
        end
    end

    % Process trial conditions
    tCon_match = cell(1, nd);
    nTrials = [];
    for id = 1:nd
        nTrials = [nTrials, length(inputStructure(id).tBaseGratingContrast)];
        tCon_match{id} = celleqel2mat_padded(inputStructure(id).tGratingContrast(1:nTrials(id)));
    end
    dirs = unique(celleqel2mat_padded(inputStructure(post).tGratingDirectionDeg(1:nTrials(post))));
    cons = unique(tCon_match{post});
    sharedCon = find(ismember(cons, targetCon));

    if length(unique(cellfun(@class, inputStructure(1).tGratingDiameterDeg, 'UniformOutput', false))) > 1
        sizes = unique(cellfun(@double, inputStructure(1).tGratingDiameterDeg));
    else
        sizes = unique(cell2mat(inputStructure(1).tGratingDiameterDeg));
    end
    sharedSize = find(ismember(sizes, targetSize));

    % Concatenate basic variables
    dirs_concat  = [dirs_concat,  dirs];
    cons_concat  = [cons_concat,  cons(sharedCon)];
    sizes_concat = [sizes_concat, sizes(sharedSize)];
    red_concat   = [red_concat,   red_cells_keep];
    green_concat = [green_concat, green_cells_keep];
    nKeep_concat = [nKeep_concat, nKeep];

    % Concatenate day-specific data
    for id = 1:nd
        tc_trial_avrg_stat_concat{id}             = cat(2, tc_trial_avrg_stat_concat{id},             tc_trial_avrg_stat{id}(:,:,sharedCon,sharedSize));
        tc_trial_avrg_loc_concat{id}              = cat(2, tc_trial_avrg_loc_concat{id},              tc_trial_avrg_loc{id}(:,:,sharedCon,sharedSize));
        tc_trial_avrg_stat_largePupil_concat{id}  = cat(2, tc_trial_avrg_stat_largePupil_concat{id},  tc_trial_avrg_stat_largePupil{id}(:,:,sharedCon,sharedSize));
        tc_trial_avrg_stat_smallPupil_concat{id}  = cat(2, tc_trial_avrg_stat_smallPupil_concat{id},  tc_trial_avrg_stat_smallPupil{id}(:,:,sharedCon,sharedSize));
        pref_responses_loc_concat{id}             = cat(1, pref_responses_loc_concat{id},             pref_responses_loc{id}(:,sharedCon,sharedSize));
        pref_responses_stat_concat{id}            = cat(1, pref_responses_stat_concat{id},            pref_responses_stat{id}(:,sharedCon,sharedSize));
        pref_responses_stat_smallPupil_concat{id} = cat(1, pref_responses_stat_smallPupil_concat{id}, pref_responses_stat_smallPupil{id}(:,sharedCon,sharedSize));
        pref_responses_stat_largePupil_concat{id} = cat(1, pref_responses_stat_largePupil_concat{id}, pref_responses_stat_largePupil{id}(:,sharedCon,sharedSize));
        RIx_concat{id}                            = cat(1, RIx_concat{id},                            sum(RIx{id}));
        norm_dir_resp_concat{id}                  = cat(1, norm_dir_resp_concat{id},                  norm_dir_resp{id}(:,:,sharedCon,sharedSize));
        pref_dir_concat{id}                       = cat(1, pref_dir_concat{id},                       prefDir_keep{id});
        noiseCorr_concat{id}                      = cat(2, noiseCorr_concat{id},                      noiseCorr{id});
        sigCorr_concat{id}                        = cat(2, sigCorr_concat{id},                        sigCorr{id});
        h_concat{id}                              = cat(1, h_concat{id},                              h_keep{id});
        data_resp_concat{id}                      = cat(1, data_resp_concat{id},                      data_resp_keep{id});
    end

    % Concatenate retinotopy data if loaded
    if isfield(instructions, 'load_retino') && instructions.load_retino && exist('ret_npSub_tc_keep', 'var')
        for id = 1:nd
            ret_npSub_tc_keep_concat{id}   = cat(2, ret_npSub_tc_keep_concat{id},   ret_npSub_tc_keep{id});
            ret_distance_keep_concat{id}   = cat(2, ret_distance_keep_concat{id},   ret_distance_keep{id});
            resp_by_stim_keep_concat{id}   = cat(3, resp_by_stim_keep_concat{id},   resp_by_stim_keep{id});
            ret_dfof_trial_keep_concat{id} = cat(2, ret_dfof_trial_keep_concat{id}, ret_dfof_trial_keep{id});
        end
    end

    if iSess == 1
        norm_diff_concat = norm_diff;
    else
        norm_diff_concat = cat(4, norm_diff_concat, norm_diff(:,sharedCon,sharedSize,:));
    end

    fprintf('Session %d of %d completed\n', iSess, nSess);
end

% Final variables
red_ind_concat   = find(red_concat);
green_ind_concat = find(green_concat);
cons             = targetCon;
nKeep_total      = sum(nKeep_concat);

% Optional: restrict cells to those responsive to a specific stimulus size
% Set instructions.sizeFilter to a size value (in deg) present in targetSize, e.g.:
%   instructions.sizeFilter = 50;  % only keep cells responsive to the 50 deg stimulus
% Leave unset or empty to use all responsive cells.
if isfield(instructions, 'sizeFilter') && ~isempty(instructions.sizeFilter)
    sizeFilterIdx = find(targetSize == instructions.sizeFilter);
    if isempty(sizeFilterIdx)
        warning('sizeFilter value %.1f not found in targetSize. Skipping size filter.', instructions.sizeFilter);
    else
        respToFilterSize = logical( ...
            sum(sum(h_concat{pre}(:,:,:,sizeFilterIdx),  2), 3) | ...
            sum(sum(h_concat{post}(:,:,:,sizeFilterIdx), 2), 3));
        sizeFilter_cells = find(respToFilterSize);
        red_ind_concat   = intersect(red_ind_concat,   sizeFilter_cells);
        green_ind_concat = intersect(green_ind_concat, sizeFilter_cells);
        fprintf('Size filter (%.1f deg): %d/%d cells retained (%d HTP+, %d HTP-)\n', ...
            instructions.sizeFilter, length(sizeFilter_cells), nKeep_total, ...
            length(red_ind_concat), length(green_ind_concat));
    end
end

% Retinotopy-filtered cell selection
% goodfit_concat is in keep-cell space; since goodfit_ind_keep is already
% the day-1 & day-2 intersection, both days are identical  AND for safety
goodfit_both = goodfit_concat{1} & goodfit_concat{2};
closeRF_both = ret_distance_retino_concat{1} < retDistThresh & ...
               ret_distance_retino_concat{2} < retDistThresh;
%UNstableRF = abs(ret_distance_retino_concat{1} - ret_distance_retino_concat{2})>2.5;
retino_cells  = find(goodfit_both & closeRF_both);
retino_red    = intersect(retino_cells, red_ind_concat);
retino_green  = intersect(retino_cells, green_ind_concat);
fprintf('Retinotopy filter (thresh=%.0f deg): %d goodfit both days, %d within thresh both days\n', ...
    retDistThresh, sum(goodfit_both), length(retino_cells));
fprintf('  HTP+: %d, HTP-: %d\n', length(retino_red), length(retino_green));

% Cell count summary table and chart
goodfit_red   = intersect(find(goodfit_both), red_ind_concat);
goodfit_green = intersect(find(goodfit_both), green_ind_concat);

countLabels = {'Matched cells', 'Keep cells', 'Keep Goodfit (both days)', 'Retino'};
counts_all   = [sum(nMatched_concat),     nKeep_total,              sum(goodfit_both),      length(retino_cells)];
counts_red   = [sum(nMatched_red_concat), length(red_ind_concat),   length(goodfit_red),    length(retino_red)];
counts_green = [sum(nMatched_grn_concat), length(green_ind_concat), length(goodfit_green),  length(retino_green)];

cellCountSummary = table(counts_all', counts_red', counts_green', ...
    'VariableNames', {'All', 'HTP_pos', 'HTP_neg'}, ...
    'RowNames', countLabels);
disp(cellCountSummary);
writetable(cellCountSummary, fullfile(fnout, 'cellCountSummary.csv'), 'WriteRowNames', true);

figure;
b = bar([counts_all; counts_red; counts_green]');
b(1).FaceColor = [0.5 0.5 0.5];
b(2).FaceColor = [0.85 0.2 0.2];
b(3).FaceColor = [0.2 0.7 0.3];
set(gca, 'XTickLabel', countLabels, 'TickDir', 'out');
box off;
legend({'All', 'HTP+', 'HTP-'}, 'Location', 'northeast');
ylabel('Cell count');
title('Cell counts through analysis pipeline');
xtickangle(20);
print(fullfile(fnout, 'cellCountSummary.pdf'), '-dpdf', '-bestfit');

% Clear individual session data, keep only concatenated
clear tc_trial_avrg_stat tc_trial_avrg_loc tc_trial_avrg_stat_largePupil tc_trial_avrg_stat_smallPupil
clear pref_responses_stat pref_responses_loc h data_resp RIx norm_dir_resp pref_dir noiseCorr sigCorr
clear nonPref_trial_avrg_stat nonPref_trial_avrg_loc data_dfof_runOnset norm_diff
clear red_cells_keep green_cells_keep prefDir_keep h_keep data_resp_keep pref_responses_stat_largePupil pref_responses_stat_smallPupil
if exist('ret_npSub_tc_matched', 'var')
    clear ret_npSub_tc_matched resp_by_stim_matched ret_dfof_trial_matched
end

exp_idx = [];
for iSess = 1:nSess
    exp_idx = [exp_idx, iSess * ones(1, nKeep_concat(iSess))];
end

%  Cell selection - find cells with running and stationary data for both days
haveRunning = cell(1, nd);
haveStat    = cell(1, nd);
for id = 1:nd
    haveRunning{id} = sum(squeeze(sum(~isnan(pref_responses_loc_concat{id}), 2)), 2) == nSize * nCon;
    haveStat{id}    = sum(squeeze(sum(~isnan(pref_responses_stat_concat{id}), 2)), 2) == nSize * nCon;
end

haveRunning_both = find(haveRunning{pre} .* haveRunning{post});
haveStat_both    = find(haveStat{pre} .* haveStat{post});
runningCells     = intersect(haveStat_both, haveRunning_both);

% Find cells responsive to small and large sizes
respToSmall = logical(sum(squeeze(sum(h_concat{pre}(:,:,:,1:2), 2)), 2) + ...
                      sum(squeeze(sum(h_concat{post}(:,:,:,1:2), 2)), 2));
respToLarge = logical(sum(squeeze(sum(h_concat{pre}(:,:,:,nSize), 2)), 2) + ...
                      sum(squeeze(sum(h_concat{post}(:,:,:,nSize), 2)), 2));

% Cell type assignments
runningGreen = intersect(runningCells, green_ind_concat);
runningRed   = intersect(runningCells, red_ind_concat);
statGreen    = green_ind_concat;
statRed      = red_ind_concat;

% Mouse indices
sessInds = cell(1, nSess);
start = 1;
for iSess = 1:nSess
    sessInds{iSess} = start:(start - 1) + nKeep_concat(iSess);
    start = start + nKeep_concat(iSess);
end

% Cell counts for each mouse
cellCounts      = nan(nSess, 2);
cellCountsGreen = nan(nSess, 2);
mouseNames      = [];

for iMouse = 1:nSess
    cellCounts(iMouse, 1)      = length(intersect(runningRed',   sessInds{iMouse}));
    cellCounts(iMouse, 2)      = length(intersect(statRed',      sessInds{iMouse}));
    cellCountsGreen(iMouse, 1) = length(intersect(runningGreen', sessInds{iMouse}));
    cellCountsGreen(iMouse, 2) = length(intersect(statGreen',    sessInds{iMouse}));
    mouseNames = [mouseNames, string(mice(iMouse,:))];
end

cellCountTableRed   = table(cellCounts,      'RowNames', mouseNames)
cellCountTableGreen = table(cellCountsGreen, 'RowNames', mouseNames)
writetable(cellCountTableRed, fullfile(fnout, 'cellCounts.csv'), 'WriteRowNames', true);

% Running by condition matrix
runningByCondition = nan(nKeep_total, nCon, nSize);
for iCon = 1:nCon
    for iSize = 1:nSize
        runningPre  = ~isnan(pref_responses_loc_concat{pre}(:, iCon, iSize));
        runningPost = ~isnan(pref_responses_loc_concat{post}(:, iCon, iSize));
        runningByCondition(:, iCon, iSize) = runningPre .* runningPost;
    end
end

%% Visualizations of stationary responses for all cells
close all
plotNeuralTimecourse(tc_trial_avrg_stat_concat, tc_trial_avrg_stat_concat, ...
    red_ind_concat, green_ind_concat, 'DayOrder', matchDrx, ...
    'UseDashedLines', [false, true], ...
    'Colors1', {'k', 'b'}, ...
    'Colors2', {'k', 'b'}, ...
    'Titles', {'HTP+', 'HTP-'}, ...
    'StimStart', 31);

figs = findobj('Type', 'figure');
sizeTitles = length(figs):-1:1;
for i = 1:length(figs)
    figure(figs(i));
    saveas(gcf, sprintf('stationary_neural_timecourse_size_%d.pdf', sizeTitles(i)));
end

%% size and contrast tuning plots
plotContrastResponse(pref_responses_stat_concat, pref_responses_stat_concat, ...
    red_ind_concat, green_ind_concat, targetCon, targetSize, 'DayOrder', matchDrx, ...
    'UseDashedLines', [false, true], ...
    'Titles', {'HTP+', 'HTP-'}, ...
    'YLabel', 'dF/F');
sgtitle('Stationary')
saveas(gcf, sprintf('stationary_contrast_response.pdf'));

plotSizeResponse(pref_responses_stat_concat, pref_responses_stat_concat, ...
    red_ind_concat, green_ind_concat, targetCon, targetSize, 'DayOrder', matchDrx, ...
    'UseDashedLines', [false, true], ...
    'Titles', {'HTP+', 'HTP-'}, ...
    'YLabel', 'dF/F');
sgtitle('Stationary')
saveas(gcf, sprintf('stationary_size_response.pdf'));

%% Size plot for each experiment separately
plotSizeResponseByExperiment(pref_responses_stat_concat, pref_responses_stat_concat, ...
    red_ind_concat, green_ind_concat, targetCon, targetSize, exp_idx, nKeep_concat, 'DayOrder', matchDrx, ...
    'UseDashedLines', [false, true], ...
    'Titles', {'HTP+', 'HTP-'}, ...
    'YLabel', 'dF/F')

%% Scatterplot and ttest for stationary trials
cell_indices = {red_ind_concat, green_ind_concat};
cell_names   = {'HTP+', 'HTP-'};
text_pos     = {[0.1, 0.1], [0.15, 0.15]};

for iSize = 1:nSize
    for iCon = 1:nCon
        figure;
        for iCellType = 1:2
            subplot(1, 2, iCellType);
            these_cells = cell_indices{iCellType};
            pre_data    = pref_responses_stat_concat{pre}(these_cells, iCon, iSize);
            post_data   = pref_responses_stat_concat{post}(these_cells, iCon, iSize);
            scatter(pre_data, post_data);
            set(gca, 'TickDir', 'out');
            grid off; box off; axis square;
            refline(1);
            title([cell_names{iCellType}, ' con ', num2str(cons(iCon)), ' size ', num2str(sizes(iSize))]);
            [~, p] = ttest(pre_data, post_data);
            if p < 0.05
                text(text_pos{iCellType}(1), text_pos{iCellType}(2), ...
                    sprintf('p = %.3f\nn = %d', p, length(these_cells)));
            end
        end
    end
end

%% Visualizations of running responses for cells that have running data
close all
plotNeuralTimecourse(tc_trial_avrg_loc_concat, tc_trial_avrg_loc_concat, ...
    runningRed, runningGreen, 'DayOrder', matchDrx, ...
    'UseDashedLines', [false, true], ...
    'Colors1', {'k', 'b'}, ...
    'Colors2', {'k', 'b'}, ...
    'Titles', {'HTP+', 'HTP-'}, ...
    'StimStart', 31);

figs = findobj('Type', 'figure');
sizeTitles = length(figs):-1:1;
for i = 1:length(figs)
    figure(figs(i));
    saveas(gcf, sprintf('running_neural_timecourse_size_%d.pdf', sizeTitles(i)));
end

plotContrastResponse(pref_responses_loc_concat, pref_responses_loc_concat, ...
    runningRed, runningGreen, cons, sizes, 'DayOrder', matchDrx, ...
    'UseDashedLines', [false, true], ...
    'Titles', {'HTP+', 'HTP-'}, ...
    'YLabel', 'dF/F');
sgtitle('Running')
saveas(gcf, sprintf('running_contrast_response.pdf'));

plotSizeResponse(pref_responses_loc_concat, pref_responses_loc_concat, ...
    runningRed, runningGreen, cons, sizes, 'DayOrder', matchDrx, ...
    'UseDashedLines', [false, true], ...
    'Titles', {'HTP+', 'HTP-'}, ...
    'YLabel', 'dF/F');
sgtitle('Running')
saveas(gcf, sprintf('running_size_response.pdf'));

plotSizeResponse_byCondition(pref_responses_loc_concat, pref_responses_loc_concat, ...
    red_ind_concat, green_ind_concat, cons, sizes, runningByCondition)

%% Normalized direction tuning
dirs_for_plotting = dirs - (length(dirs) == 8) * 180;

green_dir_avrg_stat = cell(1, nd);
red_dir_avrg_stat   = cell(1, nd);
green_dir_se_stat   = cell(1, nd);
red_dir_se_stat     = cell(1, nd);

figure('Position', [100, 100, 900, 300 * nSize]);
for iSize = 1:nSize
    for id = 1:nd
        green_data = nanmean(norm_dir_resp_concat{id}(statGreen, :, :, iSize), [3, 4]);
        green_dir_avrg_stat{id} = circshift(nanmean(green_data, 1), 4);
        green_dir_se_stat{id}   = circshift(nanstd(green_data, [], 1) / sqrt(length(statGreen)), 4);

        red_data = nanmean(norm_dir_resp_concat{id}(statRed, :, :, iSize), [3, 4]);
        red_dir_avrg_stat{id} = circshift(nanmean(red_data, 1), 4);
        red_dir_se_stat{id}   = circshift(nanstd(red_data, [], 1) / sqrt(length(statRed)), 4);
    end

    all_data   = [];
    all_errors = [];
    for id = 1:nd
        all_data   = [all_data,   green_dir_avrg_stat{id}, red_dir_avrg_stat{id}];
        all_errors = [all_errors, green_dir_se_stat{id},   red_dir_se_stat{id}];
    end
    yMin = min(all_data - all_errors);
    yMax = max(all_data + all_errors);
    padding = 0.1 * (yMax - yMin);
    yMin = yMin - padding;
    yMax = yMax + padding;

    subplot(nSize, 2, (iSize - 1) * 2 + 1);
    errorbar(dirs_for_plotting, red_dir_avrg_stat{pre},  red_dir_se_stat{pre},  'k'); hold on;
    errorbar(dirs_for_plotting, red_dir_avrg_stat{post}, red_dir_se_stat{post}, 'b');
    title(['Stationary, HTP+, size ', num2str(sizes(iSize))]);
    if iSize == nSize; ylabel('dF/F'); end

    subplot(nSize, 2, (iSize - 1) * 2 + 2);
    errorbar(dirs_for_plotting, green_dir_avrg_stat{pre},  green_dir_se_stat{pre},  '--k'); hold on;
    errorbar(dirs_for_plotting, green_dir_avrg_stat{post}, green_dir_se_stat{post}, '--b');
    title(['Stationary, HTP-, size ', num2str(sizes(iSize))]);

    for i = 1:2
        subplot(nSize, 2, (iSize - 1) * 2 + i);
        set(gca, 'TickDir', 'out'); grid off; box off; axis square;
        xticks(dirs); ylim([yMin yMax]); xlim([-10 140]);
    end
end
sgtitle('Normalized direction tuning (averaged over contrast)');
print(fullfile(fnout, 'dirTuning_allSizes.pdf'), '-dpdf', '-bestfit');

%% Plot change in pref direction between the two days
pref_dir_change = abs(pref_dir_concat{pre} - pref_dir_concat{post});
figure;
subplot(1, 2, 1); polarhistogram(pref_dir_change(green_ind_concat)); title('HTP-')
subplot(1, 2, 2); polarhistogram(pref_dir_change(red_ind_concat));   title('HTP+')
print(fullfile(fnout, 'prefDirChange.pdf'), '-dpdf', '-bestfit')

%% Plot fraction HTP+ cells suppressed and facilitated
norm_diff_red = norm_diff_concat(:,:,:,red_ind_concat);
facil_red = norm_diff_red >= 1;
supp_red  = norm_diff_red <= -1;
N = length(red_ind_concat);

facil_table_stat = squeeze(sum(facil_red(1,:,:,:), 4) / N);
supp_table_stat  = squeeze(sum(supp_red(1,:,:,:),  4) / N);

nCon   = size(facil_table_stat, 1);
nSizes = size(facil_table_stat, 2);
colors = {'k', 'r', 'b', 'g', 'm', 'c', 'y'};
if nSizes > length(colors)
    colors = [colors, repmat({'k'}, 1, nSizes - length(colors))];
end

figure;
subplot(1, 2, 1);
b = bar(1:nCon, supp_table_stat, 'grouped', 'FaceColor', "#00AFEF", 'EdgeColor', [1 1 1]);
for i = 1:nSizes; b(i).FaceColor = colors{i}; end
xticklabels(arrayfun(@num2str, cons, 'UniformOutput', false));
title('Suppressed'); ylim([0 0.6]); ylabel('Fraction HTP+ cells'); xlabel('Contrast');
set(gca, 'TickDir', 'out'); grid off; box off;

subplot(1, 2, 2);
b = bar(1:nCon, facil_table_stat, 'grouped', 'FaceColor', "#00AFEF", 'EdgeColor', [1 1 1]);
for i = 1:nSizes; b(i).FaceColor = colors{i}; end
xticklabels(arrayfun(@num2str, cons, 'UniformOutput', false));
title('Facilitated'); ylim([0 0.6]); xlabel('Contrast');
set(gca, 'TickDir', 'out'); grid off; box off;
sgtitle('Stationary');
set(gcf, 'units', 'inches', 'position', [5, 5, 3, 1.75]);
print(fullfile(fnout, 'Facil_supp_stat.pdf'), '-dpdf');

%% Plot fraction HTP- cells suppressed and facilitated
norm_diff_green = norm_diff_concat(:,:,:,green_ind_concat);
facil_green = norm_diff_green >= 1;
supp_green  = norm_diff_green <= -1;
N = length(green_ind_concat);

facil_table_stat = squeeze(sum(facil_green(1,:,:,:), 4) / N);
supp_table_stat  = squeeze(sum(supp_green(1,:,:,:),  4) / N);

nCon   = size(facil_table_stat, 1);
nSizes = size(facil_table_stat, 2);

figure;
subplot(1, 2, 1);
b = bar(1:nCon, supp_table_stat, 'grouped', 'FaceColor', "#00AFEF", 'EdgeColor', [1 1 1]);
for i = 1:nSizes; b(i).FaceColor = colors{i}; end
xticklabels(arrayfun(@num2str, cons, 'UniformOutput', false));
title('Suppressed'); ylim([0 0.6]); ylabel('Fraction HTP- cells'); xlabel('Contrast');
set(gca, 'TickDir', 'out'); grid off; box off;

subplot(1, 2, 2);
b = bar(1:nCon, facil_table_stat, 'grouped', 'FaceColor', "#00AFEF", 'EdgeColor', [1 1 1]);
for i = 1:nSizes; b(i).FaceColor = colors{i}; end
xticklabels(arrayfun(@num2str, cons, 'UniformOutput', false));
title('Facilitated'); ylim([0 0.6]); xlabel('Contrast');
set(gca, 'TickDir', 'out'); grid off; box off;
sgtitle('Stationary');
set(gcf, 'units', 'inches', 'position', [5, 5, 3, 1.75]);
print(fullfile(fnout, 'Facil_supp_stat_green.pdf'), '-dpdf');

%% Running trials - fraction HTP+ suppressed and facilitated
norm_diff_red = norm_diff_concat(:,:,:,runningRed);
facil_red = norm_diff_red >= 1;
supp_red  = norm_diff_red <= -1;
N = length(runningRed);

facil_table_loc = squeeze(sum(facil_red(2,:,:,:), 4) / N);
supp_table_loc  = squeeze(sum(supp_red(2,:,:,:),  4) / N);

figure;
subplot(1, 2, 1);
b = bar(1:nCon, supp_table_loc, 'grouped', 'FaceColor', "#00AFEF", 'EdgeColor', [1 1 1]);
b(1).FaceColor = 'k'; b(2).FaceColor = 'r';
xticklabels({'25', '50', '100'});
title('Suppressed'); ylim([0 .6]); ylabel('Fraction HTP+ cells'); xlabel('Contrast');
set(gca, 'TickDir', 'out'); box off;

subplot(1, 2, 2);
b = bar(1:nCon, facil_table_loc, 'grouped', 'FaceColor', "#00AFEF", 'EdgeColor', [1 1 1]);
b(1).FaceColor = 'k'; b(2).FaceColor = 'r';
xticklabels({'25', '50', '100'});
title('Facilitated'); ylim([0 .6]); xlabel('Contrast');
set(gca, 'TickDir', 'out'); box off;
sgtitle('Running');
set(gcf, 'units', 'inches', 'position', [5, 5, 3, 1.75]);
print(fullfile(fnout, 'Facil_supp_loc.pdf'), '-dpdf');

%% Large and small pupil stationary timecourses
close all
plotNeuralTimecourse(tc_trial_avrg_stat_smallPupil_concat, tc_trial_avrg_stat_largePupil_concat, ...
    red_ind_concat, red_ind_concat, 'DayOrder', matchDrx, ...
    'UseDashedLines', [false, false], ...
    'Colors1', {'k', 'b'}, ...
    'Colors2', {'k', 'b'}, ...
    'Titles', {'HTP+ small ', 'HTP+ large '}, ...
    'StimStart', 31);

plotNeuralTimecourse(tc_trial_avrg_stat_smallPupil_concat, tc_trial_avrg_stat_largePupil_concat, ...
    green_ind_concat, green_ind_concat, 'DayOrder', matchDrx, ...
    'UseDashedLines', [true, true], ...
    'Colors1', {'k', 'b'}, ...
    'Colors2', {'k', 'b'}, ...
    'Titles', {'HTP- small ', 'HTP- large '}, ...
    'StimStart', 31);

plotSizeResponse(pref_responses_stat_smallPupil_concat, pref_responses_stat_largePupil_concat, ...
    red_ind_concat, red_ind_concat, targetCon, targetSize, 'DayOrder', matchDrx, ...
    'UseDashedLines', [false, false], ...
    'Titles', {'small', 'large'},...
    'YLabel', 'dF/F');
sgtitle(['Stationary by pupil HTP+'])
saveas(gcf, 'pupil_stationary_size_response_HTP.pdf');

plotSizeResponse(pref_responses_stat_smallPupil_concat, pref_responses_stat_largePupil_concat, ...
    green_ind_concat, green_ind_concat, targetCon, targetSize, 'DayOrder', matchDrx, ...
    'UseDashedLines', [true, true], ...
    'Titles', {'small', 'large'},...
    'YLabel', 'dF/F');
sgtitle(['Stationary by pupil HTP-'])
saveas(gcf, 'pupil_stationary_size_response_Pyr.pdf');
%% Split cells by noise correlation on the control day and examine DART effects
figure; histogram(noiseCorr_concat{pre}(1, red_ind_concat));   xlim([-.1  1])
figure; histogram(noiseCorr_concat{pre}(1, green_ind_concat)); xlim([-.2  1.2])
figure; cdfplot(noiseCorr_concat{pre}(1, red_ind_concat))

median(noiseCorr_concat{pre}(1, red_ind_concat))

highNoiseCorr_red = find(noiseCorr_concat{pre}(1, red_ind_concat) >  .5);
lowNoiseCorr_red  = find(noiseCorr_concat{pre}(1, red_ind_concat) <= .5);

plotNeuralTimecourse(tc_trial_avrg_stat_concat, tc_trial_avrg_stat_concat, ...
    lowNoiseCorr_red, highNoiseCorr_red, 'DayOrder', matchDrx, ...
    'UseDashedLines', [false, false], ...
    'Colors1', {'k', 'b'}, ...
    'Colors2', {'k', 'b'}, ...
    'Titles', {'lowCorr', 'highCorr'}, ...
    'StimStart', 31);

%% Contrast response plots for each cell type per mouse
for iMouse = 1:nSess
    figure;
    mouseInds_this = sessInds{iMouse};
    red_this   = intersect(red_ind_concat,   mouseInds_this);
    green_this = intersect(green_ind_concat, mouseInds_this);

    subplot(1, 2, 1);
    if ~isempty(green_this)
        green_resp_pre  = squeeze(mean(pref_responses_stat_concat{pre}(green_this,  :, 2), 1, 'omitnan'));
        green_resp_post = squeeze(mean(pref_responses_stat_concat{post}(green_this, :, 2), 1, 'omitnan'));
        green_se_pre    = squeeze(std(pref_responses_stat_concat{pre}(green_this,   :, 2), 0, 1, 'omitnan')) / sqrt(length(green_this));
        green_se_post   = squeeze(std(pref_responses_stat_concat{post}(green_this,  :, 2), 0, 1, 'omitnan')) / sqrt(length(green_this));
        errorbar(cons, green_resp_pre,  green_se_pre,  '--k', 'LineWidth', 1.5); hold on;
        errorbar(cons, green_resp_post, green_se_post, '--b', 'LineWidth', 1.5);
        title(['HTP- (n=' num2str(length(green_this)) ')']);
    else
        title('HTP- (n=0)');
    end
    ylabel('dF/F'); xlabel('Contrast (%)');
    set(gca, 'TickDir', 'out'); grid off; box off;

    subplot(1, 2, 2);
    if ~isempty(red_this)
        red_resp_pre  = squeeze(mean(pref_responses_stat_concat{pre}(red_this,  :, 2), 1, 'omitnan'));
        red_resp_post = squeeze(mean(pref_responses_stat_concat{post}(red_this, :, 2), 1, 'omitnan'));
        red_se_pre    = squeeze(std(pref_responses_stat_concat{pre}(red_this,   :, 2), 0, 1, 'omitnan')) / sqrt(length(red_this));
        red_se_post   = squeeze(std(pref_responses_stat_concat{post}(red_this,  :, 2), 0, 1, 'omitnan')) / sqrt(length(red_this));
        errorbar(cons, red_resp_pre,  red_se_pre,  'k', 'LineWidth', 1.5); hold on;
        errorbar(cons, red_resp_post, red_se_post, 'b', 'LineWidth', 1.5);
        title(['HTP+ (n=' num2str(length(red_this)) ')']);
    else
        title('HTP+ (n=0)');
    end
    xlabel('Contrast (%)');
    set(gca, 'TickDir', 'out'); grid off; box off;

    sgtitle(['Contrast Response - ' mouseNames{iMouse}]);
    print(fullfile(fnout, [char(mouseNames{iMouse}) '_contrastResponse.pdf']), '-dpdf', '-bestfit');
end

%% Retinotopically-matched cells - stationary timecourses and size tuning
close all

plotNeuralTimecourse(tc_trial_avrg_stat_concat, tc_trial_avrg_stat_concat, ...
    retino_red, retino_green, 'DayOrder', matchDrx, ...
    'UseDashedLines', [false, true], ...
    'Colors1', {'k', 'b'}, ...
    'Colors2', {'k', 'b'},  ...
    'StimStart', 31);

figs = findobj('Type', 'figure');
sizeTitles = length(figs):-1:1;
for i = 1:length(figs)
    figure(figs(i));
    saveas(gcf, sprintf('retino_stationary_timecourse_size_%d.pdf', sizeTitles(i)));
end

plotSizeResponse(pref_responses_stat_concat, pref_responses_stat_concat, ...
    retino_red, retino_green, targetCon, targetSize, 'DayOrder', matchDrx, ...
    'UseDashedLines', [false, true], ...
    'Titles', {'HTP+', 'HTP-'}, ...
    'YLabel', 'dF/F');
sgtitle(['Stationary - ret distance < ' num2str(retDistThresh)])
saveas(gcf, 'retino_stationary_size_response.pdf');

retDistTable = array2table(nan(nSess, 3), 'VariableNames', {'Pre', 'Post','delta'}, 'RowNames', mouseNames);
for iMouse = 1:nSess
    idx = sessInds{iMouse};
    retDistTable.Pre(iMouse)  = mean(ret_distance_retino_concat{pre}(idx),  'omitnan');
    retDistTable.Post(iMouse) = mean(ret_distance_retino_concat{post}(idx), 'omitnan');
    retDistTable.delta(iMouse) = mean(ret_distance_retino_concat{post}(idx)-ret_distance_retino_concat{pre}(idx), 'omitnan');
end
disp(retDistTable)
writetable(retDistTable, fullfile(fnout, 'retDistance_byMouse.csv'), 'WriteRowNames', true);

retinoCellTable = array2table(nan(nSess, 2), 'VariableNames', {'HTP_pos', 'HTP_neg'}, 'RowNames', mouseNames);
for iMouse = 1:nSess
    retinoCellTable.HTP_pos(iMouse) = length(intersect(retino_red,   sessInds{iMouse}));
    retinoCellTable.HTP_neg(iMouse) = length(intersect(retino_green, sessInds{iMouse}));
end
disp(retinoCellTable)
writetable(retinoCellTable, fullfile(fnout, 'retinoCells_byMouse.csv'), 'WriteRowNames', true);

%%
retino_running       = intersect(runningCells, retino_cells);
retino_running_red   = intersect(retino_running, red_ind_concat);
retino_running_green = intersect(retino_running, green_ind_concat);

fprintf('Retino + running: %d total, %d HTP+, %d HTP-\n', ...
    length(retino_running), length(retino_running_red), length(retino_running_green));

plotSizeResponse(pref_responses_loc_concat, pref_responses_loc_concat, ...
    retino_running_red, retino_running_green, targetCon, targetSize, 'DayOrder', matchDrx, ...
    'UseDashedLines', [false, true], ...
    'Titles', {'HTP+', 'HTP-'}, ...
    'YLabel', 'dF/F');
sgtitle(['Running - ret distance < ' num2str(retDistThresh)])
saveas(gcf, 'retino_running_size_response.pdf');

close all

plotNeuralTimecourse(tc_trial_avrg_loc_concat, tc_trial_avrg_loc_concat, ...
    retino_running_red, retino_running_green, 'DayOrder', matchDrx, ...
    'UseDashedLines', [false, true], ...
    'Colors1', {'k', 'b'}, ...
    'Colors2', {'k', 'b'},  ...
    'StimStart', 31);

figs = findobj('Type', 'figure');
sizeTitles = length(figs):-1:1;
for i = 1:length(figs)
    figure(figs(i));
    saveas(gcf, sprintf('retino_loc_timecourse_size_%d.pdf', sizeTitles(i)));
end