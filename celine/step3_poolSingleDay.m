% Pool singleday_extraction.mat outputs across multiple sessions and
% visualize size tuning curves and timecourses for red (interneuron) vs
% green (pyramidal) cell populations.
%
% Session list file: plain text, one integer session index per line.

clear all; close all; clc

sess_file = input('Enter session list file: ', 's');

fid = fopen(sess_file, 'r');
if fid == -1
    error('Could not open session list file: %s', sess_file);
end

ds = '';
sess_list = [];
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line) && ~isempty(line) && ~startsWith(strtrim(line), '#')
        if isempty(ds)
            ds = strtrim(line);
        else
            num = str2double(strtrim(line));
            if ~isnan(num) && num > 0
                sess_list = [sess_list; num];
            end
        end
    end
end
fclose(fid);

fprintf('Using datasheet: %s\n', ds);
fprintf('Sessions: %s\n', mat2str(sess_list'));

run(ds);
rc = behavConstsDART;
nSess = length(sess_list);

frame_rate = 15;

% Optional: restrict analysis to cells responsive at specific stimulus sizes.
% Can specify different sizes for red (HTP+) and green (HTP-) cells.
% Uses the resp_by_size mask saved in step2.
filterRedSizeStr = input('Filter RED cells by size responsiveness? Enter size in deg or leave blank for none: ', 's');
filterGreenSizeStr = input('Filter GREEN cells by size responsiveness? Enter size in deg or leave blank for none: ', 's');

if isempty(strtrim(filterRedSizeStr))
    filterRedSize = [];
    fprintf('No size filter for red cells\n');
else
    filterRedSize = str2double(strtrim(filterRedSizeStr));
    fprintf('Filtering red cells to those responsive at %g deg\n', filterRedSize);
end

if isempty(strtrim(filterGreenSizeStr))
    filterGreenSize = [];
    fprintf('No size filter for green cells\n');
else
    filterGreenSize = str2double(strtrim(filterGreenSizeStr));
    fprintf('Filtering green cells to those responsive at %g deg\n', filterGreenSize);
end

% First pass: collect all unique contrasts and sizes across sessions
allCons = [];
allSizes = [];
sessionData = cell(nSess, 1);

fprintf('First pass: collecting all contrasts and sizes...\n');
for iSess = 1:nSess
    day_id = sess_list(iSess);
    mouse    = expt(day_id).mouse;
    expDate  = expt(day_id).date;
    runs     = expt(day_id).contrastxori_runs;
    fnin     = fullfile(rc.analysis, expt(day_id).exptType, mouse, expDate, runs{1});

    d = load(fullfile(fnin, 'singleday_extraction.mat'));
    sessionData{iSess} = d;

    allCons  = union(allCons,  d.cons);
    allSizes = union(allSizes, d.sizes);

    fprintf('  Session %d (%s %s): cons=%s, sizes=%s\n', ...
        day_id, mouse, expDate, mat2str(d.cons), mat2str(d.sizes));
end

targetCons  = allCons;
targetSizes = allSizes;
nCon  = length(targetCons);
nSize = length(targetSizes);

fprintf('\nUsing union - Contrasts: %s, Sizes: %s\n', mat2str(targetCons), mat2str(targetSizes));

% Per-session cell count tracking
sessTable = struct('mouse', {}, 'totalRed', {}, 'totalGreen', {}, 'filtRed', {}, 'filtGreen', {});

% Initialize concatenation arrays
tc_stat_concat        = [];
tc_loc_concat         = [];
tc_lp_concat          = [];
tc_sp_concat          = [];
conBySize_stat_concat = [];
conBySize_loc_concat  = [];
conBySize_lp_concat   = [];
conBySize_sp_concat   = [];
red_concat   = [];
nKeep_concat = [];
ret_distance_concat = [];

fprintf('\nSecond pass: concatenating data...\n');
for iSess = 1:nSess
    day_id = sess_list(iSess);
    d = sessionData{iSess};

    % Apply separate size filters for red and green cells if requested
    if ~isempty(filterRedSize) || ~isempty(filterGreenSize)
        if ~isfield(d, 'resp_by_size')
            error(['Session %d is missing resp_by_size. Re-run step2 to generate this field.'], day_id);
        end
        
        sizeMask = true(length(d.keep_cells), 1);  % start with all cells included
        
        % Apply red cell filter
        if ~isempty(filterRedSize)
            [~, redSizeIdx] = min(abs(d.sizes - filterRedSize));
            if abs(d.sizes(redSizeIdx) - filterRedSize) > 0.1
                warning('Session %d: requested red filter size %g not found (closest: %g) - skipping red filter', ...
                    day_id, filterRedSize, d.sizes(redSizeIdx));
            else
                redMask = logical(d.resp_by_size(:, redSizeIdx));
                sizeMask(d.red_cells & ~redMask) = false;  % exclude red cells that don't respond at target size
            end
        end
        
        % Apply green cell filter
        if ~isempty(filterGreenSize)
            [~, greenSizeIdx] = min(abs(d.sizes - filterGreenSize));
            if abs(d.sizes(greenSizeIdx) - filterGreenSize) > 0.1
                warning('Session %d: requested green filter size %g not found (closest: %g) - skipping green filter', ...
                    day_id, filterGreenSize, d.sizes(greenSizeIdx));
            else
                greenMask = logical(d.resp_by_size(:, greenSizeIdx));
                sizeMask(~d.red_cells & ~greenMask) = false;  % exclude green cells that don't respond at target size
            end
        end
    else
        sizeMask = true(length(d.keep_cells), 1);  % no filter - include all cells
    end

    nKeep   = sum(sizeMask);
    nFrames = size(d.tc_trial_avrg_stat, 1);

    fnin_i   = fullfile(rc.analysis, expt(day_id).exptType, expt(day_id).mouse, expt(day_id).date, expt(day_id).contrastxori_runs{1});
    ret_file = fullfile(fnin_i, 'retino_aligned.mat');
    if exist(ret_file, 'file')
        r = load(ret_file, 'ret_distance_keep');
        ret_distance_concat = [ret_distance_concat, r.ret_distance_keep(sizeMask')];
    else
        warning('Session %d: retino_aligned.mat not found - filling ret_distance with NaN', day_id);
        ret_distance_concat = [ret_distance_concat; nan(nKeep, 1)];
    end

    fprintf('Debug: sizeMask size = [%d %d], red_cells size = [%d %d]\n', ...
    size(sizeMask,1), size(sizeMask,2), size(d.red_cells,1), size(d.red_cells,2));
    temp_red = d.red_cells(sizeMask);
    red_concat = [red_concat; temp_red(:)];
    nKeep_concat = [nKeep_concat, nKeep];

    tc_this_stat = nan(nFrames, nKeep, nCon, nSize);
    tc_this_loc  = nan(nFrames, nKeep, nCon, nSize);
    tc_this_lp   = nan(nFrames, nKeep, nCon, nSize);
    tc_this_sp   = nan(nFrames, nKeep, nCon, nSize);

    resp_this_stat = nan(nKeep, nCon, nSize);
    resp_this_loc  = nan(nKeep, nCon, nSize);
    resp_this_lp   = nan(nKeep, nCon, nSize);
    resp_this_sp   = nan(nKeep, nCon, nSize);

    [~, conIdx]  = ismember(targetCons,  d.cons);
    [~, sizeIdx] = ismember(targetSizes, d.sizes);

    for iCon = 1:nCon
        for iSize = 1:nSize
            if conIdx(iCon) > 0 && sizeIdx(iSize) > 0
                cIdx = conIdx(iCon);
                sIdx = sizeIdx(iSize);

                tc_this_stat(:, :, iCon, iSize) = d.tc_trial_avrg_stat(:, sizeMask, cIdx, sIdx);
                tc_this_loc(:, :, iCon, iSize)  = d.tc_trial_avrg_loc(:, sizeMask, cIdx, sIdx);
                tc_this_lp(:, :, iCon, iSize)   = d.tc_trial_avrg_largePupil(:, sizeMask, cIdx, sIdx);
                tc_this_sp(:, :, iCon, iSize)   = d.tc_trial_avrg_smallPupil(:, sizeMask, cIdx, sIdx);

                resp_this_stat(:, iCon, iSize) = d.conBySize_resp_stat(sizeMask, cIdx, sIdx);
                resp_this_loc(:, iCon, iSize)  = d.conBySize_resp_loc(sizeMask, cIdx, sIdx);
                resp_this_lp(:, iCon, iSize)   = d.conBySize_resp_largePupil(sizeMask, cIdx, sIdx);
                resp_this_sp(:, iCon, iSize)   = d.conBySize_resp_smallPupil(sizeMask, cIdx, sIdx);
            end
        end
    end

    tc_stat_concat = cat(2, tc_stat_concat, tc_this_stat);
    tc_loc_concat  = cat(2, tc_loc_concat,  tc_this_loc);
    tc_lp_concat   = cat(2, tc_lp_concat,   tc_this_lp);
    tc_sp_concat   = cat(2, tc_sp_concat,   tc_this_sp);

    conBySize_stat_concat = cat(1, conBySize_stat_concat, resp_this_stat);
    conBySize_loc_concat  = cat(1, conBySize_loc_concat,  resp_this_loc);
    conBySize_lp_concat   = cat(1, conBySize_lp_concat,   resp_this_lp);
    conBySize_sp_concat   = cat(1, conBySize_sp_concat,   resp_this_sp);

    sessTable(iSess).mouse      = expt(day_id).mouse;
    sessTable(iSess).totalRed   = sum(d.red_cells);
    sessTable(iSess).totalGreen = sum(~d.red_cells);
    sessTable(iSess).filtRed    = sum(d.red_cells(sizeMask));
    sessTable(iSess).filtGreen  = sum(~d.red_cells(sizeMask));

    fprintf('  Session %d: %d/%d cells kept after size filter (%d red, %d green)\n', ...
        day_id, nKeep, length(d.keep_cells), sum(d.red_cells(sizeMask)), sum(~d.red_cells(sizeMask)));
end


red_ind   = find(red_concat);
green_ind = find(~red_concat);

% nearRet = find(ret_distance_concat<6);
% red_ind = intersect(red_ind, nearRet);
% green_ind = intersect(green_ind, nearRet);

nKeep_total = sum(nKeep_concat);
cons  = targetCons;
sizes = targetSizes;

fprintf('\nTotal: %d cells (%d red, %d green)\n', nKeep_total, length(red_ind), length(green_ind));

% Print cell count table
fprintf('\n%-20s %8s %8s %8s %10s %10s %10s %12s\n', ...
    'Mouse', 'TotRed', 'TotGreen', 'TotAll', 'FiltRed', 'FiltGreen', 'FiltAll', 'FiltFrac(%)');
fprintf('%s\n', repmat('-', 1, 90));
for iSess = 1:nSess
    totAll  = sessTable(iSess).totalRed  + sessTable(iSess).totalGreen;
    filtAll = sessTable(iSess).filtRed   + sessTable(iSess).filtGreen;
    fprintf('%-20s %8d %8d %8d %10d %10d %10d %11.1f\n', ...
        sessTable(iSess).mouse, ...
        sessTable(iSess).totalRed, sessTable(iSess).totalGreen, totAll, ...
        sessTable(iSess).filtRed,  sessTable(iSess).filtGreen,  filtAll, ...
        100 * filtAll / totAll);
end
fprintf('%s\n', repmat('-', 1, 90));
grandTotRed   = sum([sessTable.totalRed]);
grandTotGreen = sum([sessTable.totalGreen]);
grandTotAll   = grandTotRed + grandTotGreen;
grandFiltRed  = sum([sessTable.filtRed]);
grandFiltGreen= sum([sessTable.filtGreen]);
grandFiltAll  = grandFiltRed + grandFiltGreen;
fprintf('%-20s %8d %8d %8d %10d %10d %10d %11.1f\n', ...
    'TOTAL', grandTotRed, grandTotGreen, grandTotAll, ...
    grandFiltRed, grandFiltGreen, grandFiltAll, ...
    100 * grandFiltAll / grandTotAll);
fprintf('\n');
if ~isempty(filterRedSize) && ~isempty(filterGreenSize)
    fprintf('Applied size filters: RED cells responsive at %g deg, GREEN cells responsive at %g deg\n', ...
        filterRedSize, filterGreenSize);
elseif ~isempty(filterRedSize)
    fprintf('Applied size filter: RED cells responsive at %g deg (GREEN cells unfiltered)\n', filterRedSize);
elseif ~isempty(filterGreenSize)
    fprintf('Applied size filter: GREEN cells responsive at %g deg (RED cells unfiltered)\n', filterGreenSize);
else
    fprintf('No size filters applied - all responsive cells included\n');
end
fprintf('Final conditions - Contrasts: %s, Sizes: %s\n', mat2str(cons), mat2str(sizes));

% Output directory
d_str = char(string(datetime('today')));
sess_title = strjoin(arrayfun(@num2str, sess_list, 'UniformOutput', false), '_');
folder_suffix = '';
if ~isempty(filterRedSize) && ~isempty(filterGreenSize)
    folder_suffix = ['_red' num2str(filterRedSize) 'deg_green' num2str(filterGreenSize) 'deg'];
elseif ~isempty(filterRedSize)
    folder_suffix = ['_red' num2str(filterRedSize) 'deg'];
elseif ~isempty(filterGreenSize)
    folder_suffix = ['_green' num2str(filterGreenSize) 'deg'];
end
fnout = fullfile(rc.analysis, 'pooledSingleDay', ['sess' sess_title], [d_str folder_suffix]);
mkdir(fnout); cd(fnout);

% Save cell count table to text file
fidTable = fopen(fullfile(fnout, 'cellCounts.txt'), 'w');
fprintf(fidTable, '%-20s %8s %8s %8s %10s %10s %10s %12s\n', ...
    'Mouse', 'TotRed', 'TotGreen', 'TotAll', 'FiltRed', 'FiltGreen', 'FiltAll', 'FiltFrac(%)');
fprintf(fidTable, '%s\n', repmat('-', 1, 90));
for iSess = 1:nSess
    totAll  = sessTable(iSess).totalRed  + sessTable(iSess).totalGreen;
    filtAll = sessTable(iSess).filtRed   + sessTable(iSess).filtGreen;
    fprintf(fidTable, '%-20s %8d %8d %8d %10d %10d %10d %11.1f\n', ...
        sessTable(iSess).mouse, ...
        sessTable(iSess).totalRed, sessTable(iSess).totalGreen, totAll, ...
        sessTable(iSess).filtRed,  sessTable(iSess).filtGreen,  filtAll, ...
        100 * filtAll / totAll);
end
fprintf(fidTable, '%s\n', repmat('-', 1, 90));
fprintf(fidTable, '%-20s %8d %8d %8d %10d %10d %10d %11.1f\n', ...
    'TOTAL', grandTotRed, grandTotGreen, grandTotAll, ...
    grandFiltRed, grandFiltGreen, grandFiltAll, ...
    100 * grandFiltAll / grandTotAll);
fclose(fidTable);

stimStart = double(d.nOff) / 2;

%% Running data availability

runningByCondition = false(size(conBySize_loc_concat));
nRunningCells = zeros(nCon, nSize);
for iCon = 1:nCon
    for iSize = 1:nSize
        loc_vals = conBySize_loc_concat(:, iCon, iSize);
        valid_running = ~isnan(loc_vals) & abs(loc_vals) > 0.001;
        runningByCondition(:, iCon, iSize) = valid_running;
        nRunningCells(iCon, iSize) = sum(valid_running);
    end
end

% %% 1. Size tuning curves - stationary (all cells)
% fprintf('\nGenerating Plot 1: Size tuning - stationary (all cells)\n');
% plotSizeResponse_singleDay(conBySize_stat_concat, red_ind, green_ind, cons, sizes, ...
%     'Colors', {'r', 'g'}, 'Titles', {'HTP+ (red)', 'HTP- (green)'}, ...
%     'YLabel', 'dF/F', 'XLabel', 'Size (deg)');
% sgtitle('Size tuning - stationary (all cells)');
% saveas(gcf, fullfile(fnout, 'sizeResponse_stat.pdf'));
% 
% %% 2. Neural timecourses - stationary (all cells)
% fprintf('Generating Plot 2: Neural timecourses - stationary (all cells)\n');
% plotNeuralTimecourse_singleDay(tc_stat_concat, red_ind, green_ind, ...
%     'Colors', {'r', 'g'}, 'Titles', {'HTP+', 'HTP-'}, ...
%     'FrameRate', frame_rate, 'StimStart', stimStart);
% figs = findobj('Type', 'figure');
% for i = 1:length(figs)
%     saveas(figs(i), fullfile(fnout, sprintf('timecourse_stat_size%d.pdf', i)));
% end
% 
% %% 3. Size tuning curves - stationary vs running (condition-matched cells)
% fprintf('Generating Plot 3: Size tuning - stationary vs running (condition-matched cells)\n');
% plotSizeResponse_byCondition_singleDay(conBySize_stat_concat, conBySize_loc_concat, ...
%     red_ind, green_ind, cons, sizes, runningByCondition, ...
%     'Colors1', {'k', 'r'}, 'Colors2', {'k', 'r'}, ...
%     'Titles', {'HTP+ (red)', 'HTP- (green)'}, ...
%     'YLabel', 'dF/F', 'XLabel', 'Size (deg)');
% sgtitle('Size tuning - stationary vs running (condition-matched cells)');
% saveas(gcf, fullfile(fnout, 'sizeResponse_byCondition.pdf'));
% 
% %% 4. Neural timecourses - stationary vs running (condition-matched cells)
% fprintf('Generating Plot 4: Neural timecourses - stationary vs running (condition-matched cells)\n');
% plotNeuralTimecourse_byCondition_singleDay(tc_stat_concat, tc_loc_concat, ...
%     red_ind, green_ind, runningByCondition, ...
%     'Colors1', {'k', 'r'}, 'Colors2', {'k', 'r'}, ...
%     'Titles', {'HTP+ (red)', 'HTP- (green)'}, ...
%     'FrameRate', frame_rate, 'StimStart', stimStart);
% figs = findobj('Type', 'figure');
% for i = 1:length(figs)
%     saveas(figs(i), fullfile(fnout, sprintf('timecourse_byCondition_size%d.pdf', i)));
% end

%% 5. Size tuning - stationary with all contrasts overlaid (all cells)
fprintf('Generating Bonus Plot 5: Size tuning - stationary contrast overlay (all cells)\n');
plotSizeResponse_contrastOverlay_singleDay(conBySize_stat_concat, red_ind, green_ind, cons, sizes, ...
    'Titles', {'HTP+ (red)', 'HTP- (green)'}, ...
    'YLabel', 'dF/F', 'XLabel', 'Size (deg)');
sgtitle('Size tuning - stationary (all contrasts overlaid)');
saveas(gcf, fullfile(fnout, 'sizeResponse_stat_contrastOverlay.pdf'));

%% 6. Size tuning - stationary vs running overlay (condition-matched cells)
fprintf('Generating Bonus Plot 6: Size tuning - stationary vs running contrast overlay (condition-matched)\n');
plotSizeResponse_contrastOverlay_byCondition_singleDay(conBySize_stat_concat, conBySize_loc_concat, ...
    red_ind, green_ind, cons, sizes, runningByCondition, ...
    'Titles', {'HTP+ (red)', 'HTP- (green)'}, ...
    'YLabel', 'dF/F', 'XLabel', 'Size (deg)');
sgtitle('Size tuning - stationary vs running contrast overlay (condition-matched)');
saveas(gcf, fullfile(fnout, 'sizeResponse_contrastOverlay_byCondition.pdf'));

fprintf('\nAll plots completed (4 main + 2 bonus)!\n');
fprintf('Saved figures to %s\n', fnout);
%% Surround Suppression Index
% Surround Suppression Index (SSI) per cell per contrast
% SSI = (peak response - response at largest size) / peak response
% Computed from stationary trials. NaN if peak <= 0.

ssi = nan(nKeep_total, nCon);
largestSizeIdx = nSize;  % sizes are sorted ascending from step2

for iCon = 1:nCon
    resp_mat = squeeze(conBySize_stat_concat(:, iCon, :));  % nCells x nSize
    peak_resp = max(resp_mat, [], 2);
    large_resp = resp_mat(:, largestSizeIdx);
    valid = peak_resp > 0;
    ssi(valid, iCon) = (peak_resp(valid) - large_resp(valid)) ./ peak_resp(valid);
end

% Plot: bar + overlaid jittered datapoints per cell type per contrast
ssi_red   = ssi(red_ind,   :);
ssi_green = ssi(green_ind, :);

cellTypes  = {ssi_red, ssi_green};
typeLabels = {'HTP+ (red)', 'HTP- (green)'};
typeColors = [0.85 0.2 0.2; 0.2 0.65 0.2];

nGroups = nCon;
nTypes  = 2;
barWidth = 0.35;
groupX   = 1:nGroups;
offsets  = [-0.5 0.5] * barWidth;

figure; hold on;
for iType = 1:nTypes
    dat = cellTypes{iType};  % nCells x nCon
    xPos = groupX + offsets(iType);

    % Bar of mean SSI
    means = nanmean(dat, 1);
    sems  = nanstd(dat, 0, 1) ./ sqrt(sum(~isnan(dat), 1));
    bar(xPos, means, barWidth, 'FaceColor', typeColors(iType,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.6);
    errorbar(xPos, means, sems, 'k.', 'LineWidth', 1.2, 'CapSize', 5);

    % Jittered individual points (outliers removed)
    for iCon = 1:nCon
        vals = dat(:, iCon);
        vals = vals(~isnan(vals));
        outlierMask = isoutlier(vals, 'median', 'ThresholdFactor', 5);
        fprintf('  %s, contrast %g%%: removing %d/%d outliers\n', ...
            typeLabels{iType}, cons(iCon)*100, sum(outlierMask), length(vals));
        %vals = vals(~outlierMask);
        jitter = (rand(length(vals), 1) - 0.5) * barWidth * 0.6;
        scatter(xPos(iCon) + jitter, vals, 18, typeColors(iType,:), ...
            'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0);
    end
end

set(gca, 'XTick', groupX, 'XTickLabel', arrayfun(@(c) sprintf('%g%%', c*100), cons, 'UniformOutput', false), ...
    'TickDir', 'out', 'Box', 'off');
xlim([0.5 nGroups + 0.5]);
xlabel('Contrast');
ylabel('SSI');
ylim([0 3.5])
title('Surround suppression index by contrast and cell type');
hold off;

saveas(gcf, fullfile(fnout, 'SSI_byContrast.pdf'));