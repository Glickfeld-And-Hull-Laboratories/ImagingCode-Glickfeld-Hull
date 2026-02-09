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

sess_list=instructions.sess_list;
nSess = length(sess_list);
nd = 2;
frame_rate = instructions.frame_rate; 
targetCon = instructions.targetCon;
nCon = length(targetCon);
targetSize=instructions.targetSize;
nSize = length(targetSize); 

% Determine which session was used as reference for cell matching
x = instructions.refDay;
switch x
    case '1'
        pre = 1;
        post = 2;
        fprintf('Baseline used as reference\n');
        matchDrx='f';
    case '2'
        pre = 2;
        post = 1;
        fprintf('Post-DART used as reference\n');
        matchDrx='r';
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
dirs_concat = []; cons_concat = []; sizes_concat=[];

% Initialize cell arrays for concatenation
cell_vars = {'tc_trial_avrg_stat', 'tc_trial_avrg_loc', 'tc_trial_avrg_stat_largePupil', 'tc_trial_avrg_stat_smallPupil', ...
    'pref_responses_stat', 'pref_responses_loc', 'pref_responses_stat_smallPupil', 'pref_responses_stat_largePupil','h', 'data_resp', 'RIx', 'norm_dir_resp', 'pref_dir', ...
    'noiseCorr', 'sigCorr', 'nonPref_trial_avrg_stat', 'nonPref_trial_avrg_loc', 'data_dfof_runOnset'};
for i = 1:length(cell_vars)
    eval([cell_vars{i} '_concat = cell(1,nd);']);
end
pref_allTrials_stat_concat = cell(nCon, nSize, nd);
pref_allTrials_loc_concat = cell(nCon, nSize, nd);

% Initialize retinotopy concatenation variables if needed
if isfield(instructions, 'load_retino') && instructions.load_retino
    ret_npSub_tc_keep_concat = cell(1, nd);
    ret_distance_keep_concat = cell(1, nd);
    resp_by_stim_keep_concat = cell(1, nd);
    ret_dfof_trial_keep_concat = cell(1, nd);
end

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
    
    inputStructure = input;
    clear input
    % Load retinotopy data if requested
    if isfield(instructions, 'load_retino') && instructions.load_retino
        retino_file = fullfile(fn_multi, 'retino_aligned.mat');
        if exist(retino_file, 'file')
            load(fullfile(fn_multi, 'retino_aligned.mat'));
            fprintf('Loaded retinotopy data for session %d\n', iSess);
            fprintf('  ret_distance_keep structure check:\n');
        else
            warning('Retinotopy file not found for session %d: %s', iSess, retino_file);
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

    if length(unique(cellfun(@class,inputStructure(1).tGratingDiameterDeg,'UniformOutput',false))) > 1
        sizes = unique(cellfun(@double,inputStructure(1).tGratingDiameterDeg));
    else
        sizes = unique(cell2mat(inputStructure(1).tGratingDiameterDeg));
    end
    sharedSize=find(ismember(sizes, targetSize));
     
    % Concatenate basic variables
    dirs_concat = [dirs_concat, dirs];
    cons_concat = [cons_concat, cons(sharedCon)];
    sizes_concat = [sizes_concat, sizes(sharedSize)];
    red_concat = [red_concat, red_cells_keep];
    green_concat = [green_concat, green_cells_keep];
    nKeep_concat = [nKeep_concat, nKeep];
    
    % Concatenate day-specific data
    for id = 1:nd
        tc_trial_avrg_stat_concat{id} = cat(2, tc_trial_avrg_stat_concat{id}, tc_trial_avrg_stat{id}(:,:,sharedCon,sharedSize));
        tc_trial_avrg_loc_concat{id} = cat(2, tc_trial_avrg_loc_concat{id}, tc_trial_avrg_loc{id}(:,:,sharedCon,sharedSize));
        tc_trial_avrg_stat_largePupil_concat{id} = cat(2, tc_trial_avrg_stat_largePupil_concat{id}, tc_trial_avrg_stat_largePupil{id}(:,:,sharedCon,sharedSize));
        tc_trial_avrg_stat_smallPupil_concat{id} = cat(2, tc_trial_avrg_stat_smallPupil_concat{id}, tc_trial_avrg_stat_smallPupil{id}(:,:,sharedCon,sharedSize));
        pref_responses_loc_concat{id} = cat(1, pref_responses_loc_concat{id}, pref_responses_loc{id}(:,sharedCon,sharedSize));
        pref_responses_stat_concat{id} = cat(1, pref_responses_stat_concat{id}, pref_responses_stat{id}(:,sharedCon,sharedSize));
        pref_responses_stat_smallPupil_concat{id} = cat(1, pref_responses_stat_smallPupil_concat{id}, pref_responses_stat_smallPupil{id}(:,sharedCon,sharedSize));
        pref_responses_stat_largePupil_concat{id} = cat(1, pref_responses_stat_largePupil_concat{id}, pref_responses_stat_largePupil{id}(:,sharedCon,sharedSize));
        RIx_concat{id} = cat(1, RIx_concat{id}, sum(RIx{id}));
        norm_dir_resp_concat{id} = cat(1, norm_dir_resp_concat{id}, norm_dir_resp{id}(:,:,sharedCon,sharedSize));
        pref_dir_concat{id} = cat(2, pref_dir_concat{id}, prefDir_keep{id});
        noiseCorr_concat{id} = cat(2, noiseCorr_concat{id}, noiseCorr{id});
        sigCorr_concat{id} = cat(2, sigCorr_concat{id}, sigCorr{id});
        h_concat{id} = cat(1, h_concat{id}, h_keep{id});
        data_resp_concat{id} = cat(1, data_resp_concat{id}, data_resp_keep{id});
    end
    
    % Concatenate retinotopy data if loaded
    if isfield(instructions, 'load_retino') && instructions.load_retino && exist('ret_npSub_tc_keep', 'var')
        fprintf('  Concatenating retinotopy data for session %d\n', iSess);
        for id = 1:nd
            ret_npSub_tc_keep_concat{id} = cat(2, ret_npSub_tc_keep_concat{id}, ret_npSub_tc_keep{id});
            
            fprintf('    Day %d ret_distance_keep: before concat size = [%s], adding size = [%s]\n', ...
                id, num2str(size(ret_distance_keep_concat{id})), num2str(size(ret_distance_keep{id})));
            ret_distance_keep_concat{id} = cat(2, ret_distance_keep_concat{id}, ret_distance_keep{id});
            fprintf('    Day %d ret_distance_keep: after concat size = [%s]\n', id, num2str(size(ret_distance_keep_concat{id})));
            
            resp_by_stim_keep_concat{id} = cat(3, resp_by_stim_keep_concat{id}, resp_by_stim_keep{id});
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
red_ind_concat = find(red_concat);
green_ind_concat = find(green_concat);
cons = targetCon;
nKeep_total = sum(nKeep_concat);

% Summary of concatenated retinotopy data
if isfield(instructions, 'load_retino') && instructions.load_retino
    fprintf('\nFinal concatenated retinotopy data:\n');
    for id = 1:nd
        fprintf('  Day %d:\n', id);
        if iscell(ret_distance_keep_concat) && length(ret_distance_keep_concat) >= id
            fprintf('    ret_distance_keep_concat: size = [%s], numel = %d\n', ...
                num2str(size(ret_distance_keep_concat{id})), numel(ret_distance_keep_concat{id}));
        end
        if iscell(ret_npSub_tc_keep_concat) && length(ret_npSub_tc_keep_concat) >= id
            fprintf('    ret_npSub_tc_keep_concat: size = [%s]\n', num2str(size(ret_npSub_tc_keep_concat{id})));
        end
    end
end

% Clear individual session data, keep only concatenated
clear tc_trial_avrg_stat tc_trial_avrg_loc tc_trial_avrg_stat_largePupil tc_trial_avrg_stat_smallPupil
clear pref_responses_stat pref_responses_loc h data_resp RIx norm_dir_resp pref_dir noiseCorr sigCorr
clear nonPref_trial_avrg_stat nonPref_trial_avrg_loc data_dfof_runOnset norm_diff
clear red_cells_keep green_cells_keep prefDir_keep h_keep data_resp_keep pref_responses_stat_largePupil  pref_responses_stat_smallPupil
if exist('ret_npSub_tc_matched', 'var')
    clear ret_npSub_tc_matched  resp_by_stim_matched ret_dfof_trial_matched
end

exp_idx = [];
for iSess = 1:nSess
    exp_idx = [exp_idx, iSess * ones(1, nKeep_concat(iSess))];
end


%%  Cell selection - find cells with running and stationary data for both days
haveRunning = cell(1,nd);
haveStat = cell(1,nd);
for id = 1:nd
    haveRunning{id} = sum(squeeze(sum(~isnan(pref_responses_loc_concat{id}),2)),2) == nSize*nCon;
    haveStat{id} = sum(squeeze(sum(~isnan(pref_responses_stat_concat{id}),2)),2) == nSize*nCon;
end

haveRunning_both = find(haveRunning{pre} .* haveRunning{post});
haveStat_both = find(haveStat{pre} .* haveStat{post});
runningCells = intersect(haveStat_both, haveRunning_both);

% Find cells responsive to small and large sizes
respToSmall = logical(sum(squeeze(sum(h_concat{pre}(:,:,:,1:2),2)),2) + ...
                     sum(squeeze(sum(h_concat{post}(:,:,:,1:2),2)),2));
respToLarge = logical(sum(squeeze(sum(h_concat{pre}(:,:,:,nSize),2)),2) + ...
                     sum(squeeze(sum(h_concat{post}(:,:,:,nSize),2)),2));


% Cell type assignments
runningGreen = intersect(runningCells, green_ind_concat);
runningRed = intersect(runningCells, red_ind_concat);
statGreen = green_ind_concat;
statRed = red_ind_concat;

% Mouse indices
sessInds = cell(1,nSess);
start = 1;
for iSess = 1:nSess
    sessInds{iSess} = start:(start-1)+nKeep_concat(iSess);
    start = start + nKeep_concat(iSess);
end

% Cell counts for each mouse
cellCounts = nan(nSess,2);
cellCountsGreen = nan(nSess,2);
mouseNames = [];

for iMouse = 1:nSess
    cellCounts(iMouse,1) = length(intersect(runningRed', sessInds{iMouse}));
    cellCounts(iMouse,2) = length(intersect(statRed', sessInds{iMouse}));
    cellCountsGreen(iMouse,1) = length(intersect(runningGreen', sessInds{iMouse}));
    cellCountsGreen(iMouse,2) = length(intersect(statGreen', sessInds{iMouse}));
    mouseNames = [mouseNames, string(mice(iMouse,:))];
end

% Save cell count tables
cellCountTableRed = table(cellCounts, RowNames=mouseNames)
cellCountTableGreen = table(cellCountsGreen, RowNames=mouseNames)
writetable(cellCountTableRed, fullfile(fnout,'cellCounts.csv'), 'WriteRowNames', true);

% Running by condition matrix    - this indicates which cells are matched
% for behavioral state within a single stimulus condition, rather than the
% more strict cirterion of being matched across all stimulus conditions


runningByCondition = nan(nKeep_total, nCon, nSize);
for iCon = 1:nCon
    for iSize = 1:nSize
        runningPre = ~isnan(pref_responses_loc_concat{pre}(:,iCon,iSize));
        runningPost = ~isnan(pref_responses_loc_concat{post}(:,iCon,iSize));
        runningByCondition(:,iCon,iSize) = runningPre .* runningPost;
    end
end




%% Visualizations of stationary responses for all cells
close all
% Plot stationary timecourses for all cells
plotNeuralTimecourse(tc_trial_avrg_stat_concat, tc_trial_avrg_stat_concat, ...
    red_ind_concat, green_ind_concat,'DayOrder',matchDrx, ...
    'UseDashedLines', [false, true], ...
    'Colors1', {'k', 'b'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'b'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'HTP+', 'HTP-'}, ...
    'StimStart', 31);
% Saves all the figures

figs = findobj('Type', 'figure');
sizeTitles = length(figs):-1:1;
for i = 1:length(figs)
    figure(figs(i));
    saveas(gcf, sprintf('stationary_neural_timecourse_size_%d.pdf', sizeTitles(i)));
end
%%

% Plot stationary contrast response function for all cells
plotContrastResponse(pref_responses_stat_concat, pref_responses_stat_concat, ...
    red_ind_concat, green_ind_concat, targetCon,targetSize,'DayOrder',matchDrx, ...
    'UseDashedLines', [false, true], ...  % Dashed lines for the right plot
    'Titles', {'HTP+', 'HTP-'}, ...
    'YLabel', 'dF/F');
sgtitle('Stationary')
saveas(gcf, sprintf('stationary_contrast_response.pdf'));

% Plot stationary size response function for all cells
plotSizeResponse(pref_responses_stat_concat, pref_responses_stat_concat, ...
    red_ind_concat, green_ind_concat, targetCon,targetSize,'DayOrder',matchDrx, ...
    'UseDashedLines', [false, true], ...  % Dashed lines for the right plot
    'Titles', {'HTP+', 'HTP-'}, ...
    'YLabel', 'dF/F');
sgtitle('Stationary')
saveas(gcf, sprintf('stationary_size_response.pdf'));

%% Scatterplot and ttest for stationary trials

cell_indices = {red_ind_concat, green_ind_concat};
cell_names = {'HTP+', 'HTP-'};
text_pos = {[0.1, 0.1], [0.15, 0.15]};

for iSize = 1:nSize
    for iCon = 1:nCon
        figure;
        
        for iCellType = 1:2
            subplot(1, 2, iCellType);
            
            % Get data for this cell type
            these_cells = cell_indices{iCellType};
            pre_data = pref_responses_stat_concat{pre}(these_cells, iCon, iSize);
            post_data = pref_responses_stat_concat{post}(these_cells, iCon, iSize);
            
            % Calculate axis limits
            all_data = [pre_data; post_data];
            data_max = max(all_data) * 1.1;
            axis_lim = [0, max(0.2, data_max)];
            
            % Plot
            scatter(pre_data, post_data);
            
            
            % Format
            set(gca, 'TickDir', 'out');
            grid off;
            box off;
            axis square;
           % xlim(axis_lim);
            %ylim(axis_lim);
            refline(1);
            title([cell_names{iCellType}, ' con ', num2str(cons(iCon)), ' size ', num2str(sizes(iSize))]);
            
            % Statistical test
            [~, p] = ttest(pre_data, post_data);
            if p < 0.05
                p_text = sprintf('p = %.3f\nn = %d', p, length(these_cells));
                text(text_pos{iCellType}(1), text_pos{iCellType}(2), p_text);
            end
        end
    end
end
%% Visualizations of running responses for cells that have running data
close all
% Plot running timecourses for cells that have running data
plotNeuralTimecourse(tc_trial_avrg_loc_concat, tc_trial_avrg_loc_concat, ...
    runningRed, runningGreen,'DayOrder',matchDrx, ...
    'UseDashedLines', [false, true], ...
    'Colors1', {'k', 'b'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'b'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'HTP+', 'HTP-'}, ...
    'StimStart', 31);

figs = findobj('Type', 'figure');
sizeTitles = length(figs):-1:1;
for i = 1:length(figs)
    figure(figs(i));
    saveas(gcf, sprintf('running_neural_timecourse_size_%d.pdf', sizeTitles(i)));
end

% Plot running contrast response function for cells that have running data
plotContrastResponse(pref_responses_loc_concat, pref_responses_loc_concat, ...
    runningRed, runningGreen, cons,sizes,'DayOrder',matchDrx, ...
    'UseDashedLines', [false, true], ...  % Dashed lines for the right plot
    'Titles', {'HTP+', 'HTP-'}, ...
    'YLabel', 'dF/F');
sgtitle('Running')
saveas(gcf, sprintf('running_contrast_response.pdf'));

% Plot running size response function for for cells that have running data
plotSizeResponse(pref_responses_loc_concat, pref_responses_loc_concat, ...
    runningRed, runningGreen, cons,sizes,'DayOrder',matchDrx, ...
    'UseDashedLines', [false, true], ...  % Dashed lines for the right plot
    'Titles', {'HTP+', 'HTP-'}, ...
    'YLabel', 'dF/F');
sgtitle('Running')
saveas(gcf, sprintf('running_size_response.pdf'));

%% Normalized direction tuning 
dirs_for_plotting = dirs - (length(dirs) == 8) * 180;

green_dir_avrg_stat = cell(1, nd);
red_dir_avrg_stat = cell(1, nd);
green_dir_se_stat = cell(1, nd);
red_dir_se_stat = cell(1, nd);

figure('Position', [100, 100, 900, 300*nSize]);
for iSize = 1:nSize
    for id = 1:nd
        green_data = nanmean(norm_dir_resp_concat{id}(statGreen, :, :, iSize), [3, 4]);
        green_dir_avrg_stat{id} = circshift(nanmean(green_data, 1), 4);
        green_dir_se_stat{id} = circshift(nanstd(green_data, [], 1) / sqrt(length(statGreen)), 4);
        
        red_data = nanmean(norm_dir_resp_concat{id}(statRed, :, :, iSize), [3, 4]);
        red_dir_avrg_stat{id} = circshift(nanmean(red_data, 1), 4);
        red_dir_se_stat{id} = circshift(nanstd(red_data, [], 1) / sqrt(length(statRed)), 4);
    end
    
    all_data = [];
    all_errors = [];
    for id = 1:nd
        all_data = [all_data, green_dir_avrg_stat{id}, red_dir_avrg_stat{id}];
        all_errors = [all_errors, green_dir_se_stat{id}, red_dir_se_stat{id}];
    end
    yMin = min(all_data - all_errors);
    yMax = max(all_data + all_errors);
    padding = 0.1 * (yMax - yMin);
    yMin = yMin - padding;
    yMax = yMax + padding;
    
    subplot(nSize, 2, (iSize-1)*2 + 1);
    errorbar(dirs_for_plotting, red_dir_avrg_stat{pre}, red_dir_se_stat{pre}, 'k');
    hold on;
    errorbar(dirs_for_plotting, red_dir_avrg_stat{post}, red_dir_se_stat{post}, 'b');
    title(['Stationary, HTP+, size ', num2str(sizes(iSize))]);
    if iSize == nSize
        ylabel('dF/F');
    end
    
    subplot(nSize, 2, (iSize-1)*2 + 2);
    errorbar(dirs_for_plotting, green_dir_avrg_stat{pre}, green_dir_se_stat{pre}, '--k');
    hold on;
    errorbar(dirs_for_plotting, green_dir_avrg_stat{post}, green_dir_se_stat{post}, '--b');
    title(['Stationary, HTP-, size ', num2str(sizes(iSize))]);
    
    for i = 1:2
        subplot(nSize, 2, (iSize-1)*2 + i);
        set(gca, 'TickDir', 'out');
        grid off;
        box off;
        axis square;
        xticks(dirs);
        ylim([yMin yMax]);
        xlim([-10 140]);
    end
end
sgtitle('Normalized direction tuning (averaged over contrast)');
print(fullfile(fnout, 'dirTuning_allSizes.pdf'), '-dpdf', '-bestfit');


%% Plot change in pref direction between the two days
pref_dir_change=abs(pref_dir_concat{pre}-pref_dir_concat{post});
figure;subplot(1,2,1)
polarhistogram(pref_dir_change(green_ind_concat))
title('HTP-')
subplot(1,2,2)
polarhistogram(pref_dir_change(red_ind_concat))
title('HTP+')
print(fullfile(fnout,'prefDirChange.pdf'),'-dpdf','-bestfit')

%% Plot fraction HTP+ cells suppressed and facilitated 
% As a function of size and contrast, for stationary data
norm_diff_red = norm_diff_concat(:,:,:,red_ind_concat);
facil_red = norm_diff_red >= 1;
supp_red = norm_diff_red <= -1;
N = length(red_ind_concat);

facil_table_stat = squeeze(sum(facil_red(1,:,:,:), 4) / N);
supp_table_stat = squeeze(sum(supp_red(1,:,:,:), 4) / N);

% Get dimensions
nCon = size(facil_table_stat, 1);
nSizes = size(facil_table_stat, 2);

% Define colors for different sizes
colors = {'k', 'r', 'b', 'g', 'm', 'c', 'y'};
if nSizes > length(colors)
    colors = [colors, repmat({'k'}, 1, nSizes - length(colors))];
end

figure;

% Suppressed
subplot(1, 2, 1);
b = bar(1:nCon, supp_table_stat, 'grouped', 'FaceColor', "#00AFEF", 'EdgeColor', [1 1 1]);
for i = 1:nSizes
    b(i).FaceColor = colors{i};
end
xticklabels(arrayfun(@num2str, cons, 'UniformOutput', false));
title('Suppressed');
ylim([0 0.6]);
ylabel('Fraction HTP+ cells');
xlabel('Contrast');
set(gca, 'TickDir', 'out');
grid off;
box off;

% Facilitated
subplot(1, 2, 2);
b = bar(1:nCon, facil_table_stat, 'grouped', 'FaceColor', "#00AFEF", 'EdgeColor', [1 1 1]);
for i = 1:nSizes
    b(i).FaceColor = colors{i};
end
xticklabels(arrayfun(@num2str, cons, 'UniformOutput', false));
title('Facilitated');
ylim([0 0.6]);
xlabel('Contrast');
set(gca, 'TickDir', 'out');
grid off;
box off;

sgtitle('Stationary');

% Set figure size
x0 = 5;
y0 = 5;
width = 3;
height = 1.75;
set(gcf, 'units', 'inches', 'position', [x0, y0, width, height]);

print(fullfile(fnout, 'Facil_supp_stat.pdf'), '-dpdf');

%% Plot fraction HTP- cells suppressed and facilitated 
% As a function of size and contrast, for stationary data
norm_diff_green = norm_diff_concat(:,:,:,green_ind_concat);
facil_green = norm_diff_green >= 1;
supp_green = norm_diff_green <= -1;
N = length(green_ind_concat);

facil_table_stat = squeeze(sum(facil_green(1,:,:,:), 4) / N);
supp_table_stat = squeeze(sum(supp_green(1,:,:,:), 4) / N);

% Get dimensions
nCon = size(facil_table_stat, 1);
nSizes = size(facil_table_stat, 2);

% Define colors for different sizes
colors = {'k', 'r', 'b', 'g', 'm', 'c', 'y'};
if nSizes > length(colors)
    colors = [colors, repmat({'k'}, 1, nSizes - length(colors))];
end

figure;

% Suppressed
subplot(1, 2, 1);
b = bar(1:nCon, supp_table_stat, 'grouped', 'FaceColor', "#00AFEF", 'EdgeColor', [1 1 1]);
for i = 1:nSizes
    b(i).FaceColor = colors{i};
end
xticklabels(arrayfun(@num2str, cons, 'UniformOutput', false));
title('Suppressed');
ylim([0 0.6]);
ylabel('Fraction HTP- cells');
xlabel('Contrast');
set(gca, 'TickDir', 'out');
grid off;
box off;

% Facilitated
subplot(1, 2, 2);
b = bar(1:nCon, facil_table_stat, 'grouped', 'FaceColor', "#00AFEF", 'EdgeColor', [1 1 1]);
for i = 1:nSizes
    b(i).FaceColor = colors{i};
end
xticklabels(arrayfun(@num2str, cons, 'UniformOutput', false));
title('Facilitated');
ylim([0 0.6]);
xlabel('Contrast');
set(gca, 'TickDir', 'out');
grid off;
box off;

sgtitle('Stationary');

% Set figure size
x0 = 5;
y0 = 5;
width = 3;
height = 1.75;
set(gcf, 'units', 'inches', 'position', [x0, y0, width, height]);

print(fullfile(fnout, 'Facil_supp_stat_green.pdf'), '-dpdf');
%% Same as above for running trials, using only HTP+ cells that have running data
norm_diff_red = norm_diff_concat(:,:,:,runningRed);
facil_red=norm_diff_red(:,:,:,:)>=1;
supp_red=norm_diff_red(:,:,:,:)<=-1;

N=length(runningRed);

 facil_table_loc = squeeze(sum(facil_red(2,:,:,:),4)/N);
 supp_table_loc = squeeze(sum(supp_red(2,:,:,:),4)/N);
    
    figure;
    subplot(1,2,1)
    % b=bar([1,2,3],[supp_table_loc(:,1),supp_table_loc(:,2),supp_table_loc(:,3)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
    b = bar(1:nCon, supp_table_loc, 'grouped', 'FaceColor', "#00AFEF", 'EdgeColor', [1 1 1]);
    b(1).FaceColor='k'
    b(2).FaceColor='r'
    % b(3).FaceColor='b'
    xticklabels({'25','50','100'})
    hold on
    title('Suppressed')
    ylim([0 .6])
    ylabel(["Fraction HTP+ cells"]) 
    xlabel(["Contrast"])
    set(gca,'TickDir','out')
    box off
    
    subplot(1,2,2)
    % b=bar([1,2,3],[facil_table_loc(:,1),facil_table_loc(:,2),facil_table_loc(:,3)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
    b = bar(1:nCon, facil_table_loc, 'grouped', 'FaceColor', "#00AFEF", 'EdgeColor', [1 1 1]);
    b(1).FaceColor='k'
    b(2).FaceColor='r'
    % b(3).FaceColor='b'
    xticklabels({'25','50','100'})
    hold on
    title('Facilitated')
    ylim([0 .6])
    %ylabel(["Fraction HTP+ cells"]) 
    xlabel(["Contrast"])
    set(gca,'TickDir','out')
    box off
    sgtitle('Running')
    x0=5;
    y0=5;
    width=3;
    height=1.75;
    set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'Facil_supp_loc.pdf'),'-dpdf');



%% plot large and small pupil stationary timecourses
close all
plotNeuralTimecourse(tc_trial_avrg_stat_smallPupil_concat, tc_trial_avrg_stat_largePupil_concat, ...
    red_ind_concat, red_ind_concat,'DayOrder',matchDrx, ...
    'UseDashedLines', [false, false], ...
    'Colors1', {'k', 'b'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'b'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'HTP+ small ', 'HTP+ large '}, ...
    'StimStart', 31);

plotNeuralTimecourse(tc_trial_avrg_stat_smallPupil_concat, tc_trial_avrg_stat_largePupil_concat, ...
    green_ind_concat, green_ind_concat, 'DayOrder',matchDrx,...
    'UseDashedLines', [true, true], ...
    'Colors1', {'k', 'b'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'b'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'HTP- small ', 'HTP- large '}, ...
    'StimStart', 31);


%% Split cells by noise correlation on the control day and examine DART effects

%within each cell of noiseCorr_concat is a 2 X nNeurons matrix, where on the
%first dimension 1=R and 2=p values.
figure;histogram(noiseCorr_concat{pre}(1,red_ind_concat));xlim([-.1 1])
figure;histogram(noiseCorr_concat{pre}(1,green_ind_concat));xlim([-.2 1.2])
figure;cdfplot(noiseCorr_concat{pre}(1,red_ind_concat))

median(noiseCorr_concat{pre}(1,red_ind_concat))

highNoiseCorr_red=find(noiseCorr_concat{pre}(1,red_ind_concat)>.3);
lowNoiseCorr_red=find(noiseCorr_concat{pre}(1,red_ind_concat)<=.3);


plotNeuralTimecourse(tc_trial_avrg_stat_concat, tc_trial_avrg_stat_concat, ...
    lowNoiseCorr_red, highNoiseCorr_red,'DayOrder',matchDrx, ...
    'UseDashedLines', [false, false], ...
    'Colors1', {'k', 'b'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'b'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'lowCorr', 'highCorr'}, ...
    'StimStart', 31);

%% Split cells by noise correlation on the control day and examine DART effects

%within each cell of noiseCorr_concat is a 2 X nNeurons matrix, where on the
%first dimension 1=R and 2=p values.
figure;histogram(noiseCorr_concat{pre}(1,red_ind_concat));xlim([-.1 1])
figure;histogram(noiseCorr_concat{pre}(1,green_ind_concat));xlim([-.2 1.2])
figure;cdfplot(noiseCorr_concat{pre}(1,red_ind_concat))

median(noiseCorr_concat{pre}(1,red_ind_concat))

highNoiseCorr_red=find(noiseCorr_concat{pre}(1,red_ind_concat)>.3);
lowNoiseCorr_red=find(noiseCorr_concat{pre}(1,red_ind_concat)<=.3);


plotNeuralTimecourse(tc_trial_avrg_stat_concat, tc_trial_avrg_stat_concat, ...
    highNoiseCorr_red, highNoiseCorr_red, 'DayOrder',matchDrx,...
    'UseDashedLines', [false, false], ...
    'Colors1', {'k', 'b'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'b'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'lowCorr', 'highCorr'}, ...
    'StimStart', 31);

%% Contrast response plots for each cell type per mouse
for iMouse = 1:nSess
    figure;
    mouseInds_this = sessInds{iMouse};
    
    red_this = intersect(red_ind_concat, mouseInds_this);
    green_this = intersect(green_ind_concat, mouseInds_this);
    
    subplot(1, 2, 1);
    if ~isempty(green_this)
        green_resp_pre = squeeze(mean(pref_responses_stat_concat{pre}(green_this, :, 2), 1, 'omitnan'));
        green_resp_post = squeeze(mean(pref_responses_stat_concat{post}(green_this, :, 2), 1, 'omitnan'));
        green_se_pre = squeeze(std(pref_responses_stat_concat{pre}(green_this, :, 2), 0, 1, 'omitnan')) / sqrt(length(green_this));
        green_se_post = squeeze(std(pref_responses_stat_concat{post}(green_this, :, 2), 0, 1, 'omitnan')) / sqrt(length(green_this));
        
        errorbar(cons, green_resp_pre, green_se_pre, '--k', 'LineWidth', 1.5);
        hold on;
        errorbar(cons, green_resp_post, green_se_post, '--b', 'LineWidth', 1.5);
        title(['HTP- (n=' num2str(length(green_this)) ')']);
    else
        title('HTP- (n=0)');
    end
    ylabel('dF/F');
    xlabel('Contrast (%)');
    set(gca, 'TickDir', 'out');
    grid off;
    box off;
    
    subplot(1, 2, 2);
    if ~isempty(red_this)
        red_resp_pre = squeeze(mean(pref_responses_stat_concat{pre}(red_this, :, 2), 1, 'omitnan'));
        red_resp_post = squeeze(mean(pref_responses_stat_concat{post}(red_this, :, 2), 1, 'omitnan'));
        red_se_pre = squeeze(std(pref_responses_stat_concat{pre}(red_this, :, 2), 0, 1, 'omitnan')) / sqrt(length(red_this));
        red_se_post = squeeze(std(pref_responses_stat_concat{post}(red_this, :, 2), 0, 1, 'omitnan')) / sqrt(length(red_this));
        
        errorbar(cons, red_resp_pre, red_se_pre, 'k', 'LineWidth', 1.5);
        hold on;
        errorbar(cons, red_resp_post, red_se_post, 'b', 'LineWidth', 1.5);
        title(['HTP+ (n=' num2str(length(red_this)) ')']);
    else
        title('HTP+ (n=0)');
    end
    xlabel('Contrast (%)');
    set(gca, 'TickDir', 'out');
    grid off;
    box off;
    
    sgtitle(['Contrast Response - ' mouseNames{iMouse}]);
    print(fullfile(fnout, [char(mouseNames{iMouse}) '_contrastResponse.pdf']), '-dpdf', '-bestfit');
end

%% single cell scatter SSIx and DartIx

