% Pool df/F and MAD variability outputs across sessions and plot together.
% Loads singleday_extraction.mat and singleday_extraction_MAD.mat per session.
%
% Session list txt format:
%   Line 1 (uncommented): datasheet .m file name
%   Remaining uncommented lines: session indices (day_id) to pool

if ~exist('sess_file', 'var')
    sess_file = input('Enter session list file: ', 's');
end

fid = fopen(sess_file, 'r');
lines = {};
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if ischar(line) && ~isempty(line) && ~startsWith(line, '#')
        lines{end+1} = line;
    end
end
fclose(fid);

ds      = lines{1};
day_ids = cellfun(@str2num, lines(2:end));

run(ds);
rc = behavConstsDART;
dataStructLabels = {'contrastxori'};

% First pass: load all sessions and collect union of cons/sizes
nSess       = length(day_ids);
sessionData = cell(nSess, 1);
allCons     = [];
allSizes    = [];
ExperimentFolder = '';

for i = 1:nSess
    iDay             = day_ids(i);
    mouse            = expt(iDay).mouse;
    expDate          = expt(iDay).date;
    ExperimentFolder = expt(iDay).exptType;
    runs             = eval(['expt(iDay).' cell2mat(dataStructLabels) '_runs']);
    fnday            = fullfile(rc.analysis, ExperimentFolder, mouse, expDate, runs{1});

    matDF  = fullfile(fnday, 'singleday_extraction.mat');
    matMAD = fullfile(fnday, 'singleday_extraction_MAD.mat');

    if ~exist(matDF, 'file') || ~exist(matMAD, 'file')
        fprintf('Session %d: missing mat file(s), skipping\n', iDay);
        continue
    end

    d = load(matDF,  'tc_trial_avrg_stat', 'conBySize_resp_stat', ...
                     'red_cells', 'cons', 'sizes', 'nOn', 'nOff');
    m = load(matMAD, 'tc_trial_MAD_stat', 'conBySize_MAD_stat');
    d.tc_MAD  = m.tc_trial_MAD_stat;
    d.cbsMAD  = m.conBySize_MAD_stat;
    sessionData{i} = d;

    allCons  = union(allCons,  d.cons);
    allSizes = union(allSizes, d.sizes);
    fprintf('Session %d (%s %s): %d cells\n', iDay, mouse, expDate, size(d.tc_trial_avrg_stat, 2));
end

sessionData = sessionData(~cellfun(@isempty, sessionData));
nSess   = length(sessionData);
nCon    = length(allCons);
nSize   = length(allSizes);
nFrames = size(sessionData{1}.tc_trial_avrg_stat, 1);

% Second pass: concatenate into (nFrames x nCells x nCon x nSize)
tc_dfof_all   = [];
tc_MAD_all    = [];
cbsDFOF_all   = [];
cbsMAD_all    = [];
red_cells_all = [];

for i = 1:nSess
    d         = sessionData{i};
    nKeep_day = size(d.tc_trial_avrg_stat, 2);

    [~, conIdx]  = ismember(allCons,  d.cons);
    [~, sizeIdx] = ismember(allSizes, d.sizes);

    tc_dfof_this  = nan(nFrames, nKeep_day, nCon, nSize);
    tc_MAD_this   = nan(nFrames, nKeep_day, nCon, nSize);
    cbs_dfof_this = nan(nKeep_day, nCon, nSize);
    cbs_MAD_this  = nan(nKeep_day, nCon, nSize);

    for iCon = 1:nCon
        for iSize = 1:nSize
            if conIdx(iCon) > 0 && sizeIdx(iSize) > 0
                cIdx = conIdx(iCon);
                sIdx = sizeIdx(iSize);
                tc_dfof_this(:, :, iCon, iSize)  = d.tc_trial_avrg_stat(:, :, cIdx, sIdx);
                tc_MAD_this(:, :, iCon, iSize)   = d.tc_MAD(:, :, cIdx, sIdx);
                cbs_dfof_this(:, iCon, iSize)    = d.conBySize_resp_stat(:, cIdx, sIdx);
                cbs_MAD_this(:, iCon, iSize)     = d.cbsMAD(:, cIdx, sIdx);
            end
        end
    end

    tc_dfof_all   = cat(2, tc_dfof_all,   tc_dfof_this);
    tc_MAD_all    = cat(2, tc_MAD_all,    tc_MAD_this);
    cbsDFOF_all   = cat(1, cbsDFOF_all,   cbs_dfof_this);
    cbsMAD_all    = cat(1, cbsMAD_all,    cbs_MAD_this);
    red_cells_all = [red_cells_all; d.red_cells(:)];

    fprintf('  %d cells (cumulative: %d)\n', nKeep_day, length(red_cells_all));
end

HTP_cells   = logical(red_cells_all);
nKeep_total = length(red_cells_all);
nOn   = sessionData{1}.nOn;
nOff  = sessionData{1}.nOff;
cons  = allCons;
sizes = allSizes;

fprintf('\nTotal pooled: %d cells (%d HTP+, %d HTP-)\n', nKeep_total, sum(HTP_cells), sum(~HTP_cells));

% Save
fnpool = fullfile(rc.analysis, ExperimentFolder);
save(fullfile(fnpool, 'pooled_dfof_and_MAD.mat'), ...
    'tc_dfof_all', 'tc_MAD_all', 'cbsDFOF_all', 'cbsMAD_all', ...
    'red_cells_all', 'HTP_cells', 'cons', 'sizes', 'nOn', 'nOff', '-v7.3');
fprintf('Saved to %s\n', fullfile(fnpool, 'pooled_dfof_and_MAD.mat'));

% Output folder
d_str      = char(string(datetime('today')));
sess_title = strjoin(arrayfun(@num2str, day_ids, 'UniformOutput', false), '_');
fnout      = fullfile(rc.analysis, 'pooled_dfof_and_MAD', ['sess' sess_title], d_str);
mkdir(fnout);

%% Plot: one figure per contrast
% Rows: [df/F ; MAD], Columns: sizes
% HTP+ (green), HTP- (blue), mean +/- SEM shading

t_frames = (1:(nOn+nOff)) - nOff/2;
clrs     = [0.2 0.65 0.2; 0.15 0.45 0.85];  % green = HTP+, blue = HTP-
masks    = {HTP_cells, ~HTP_cells};
lbls     = {'HTP+', 'HTP-'};
metricTC = {tc_dfof_all, tc_MAD_all};
rowYLbl  = {'\DeltaF/F', 'MAD (\DeltaF/F)'};

for iCon = 1:nCon
    figure('Position', [50 50 200*nSize 420]);
    axH = gobjects(2, nSize);
    for iRow = 1:2   % row 1 = df/F, row 2 = MAD
        for iSize = 1:nSize
            axH(iRow, iSize) = subplot(2, nSize, (iRow-1)*nSize + iSize);
            hold on
            for iType = 1:2
                tc = metricTC{iRow}(:, masks{iType}, iCon, iSize);
                m  = mean(tc, 2, 'omitnan');
                se = std(tc, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(tc(1,:))));
                fill([t_frames fliplr(t_frames)], [(m+se)' fliplr((m-se)')], ...
                    clrs(iType,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
                plot(t_frames, m, 'Color', clrs(iType,:), 'LineWidth', 1.2, ...
                    'DisplayName', lbls{iType});
            end
            xline(0,   'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
            xline(nOn, 'k:',  'LineWidth', 0.5, 'HandleVisibility', 'off');
            set(gca, 'TickDir', 'out', 'Box', 'off');
            xlabel('Frame');
            if iSize == 1
                ylabel(rowYLbl{iRow});
            end
            if iRow == 1
                title(sprintf('%g deg', sizes(iSize)));
            end

        end
    end
    % Shared y limits per row
    for iRow = 1:1
        yl = cell2mat(arrayfun(@(ax) ax.YLim, axH(iRow,:), 'UniformOutput', false)');
        set(axH(iRow,:), 'YLim', [min(yl(:,1)) max(yl(:,2))]);
    end
    sgtitle(sprintf('Contrast %g%% � stationary | green: HTP+, blue: HTP-', cons(iCon)*100));
    saveas(gcf, fullfile(fnout, sprintf('timecourse_dfof_MAD_con%g.pdf', cons(iCon)*100)));
end
fprintf('Plots saved to %s\n', fnout);