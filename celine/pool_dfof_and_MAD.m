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

%% Plot: contrasts overlaid, separate figures for HTP+ and HTP-
% Rows: [df/F ; MAD], Columns: sizes, contrasts as light-to-dark blue

t_frames  = (1:(nOn+nOff)) - nOff/2;
con_clrs  = [0.75 0.88 1; 0.2 0.5 0.85; 0.05 0.15 0.5];
rowYLbl   = {'\DeltaF/F', 'MAD (\DeltaF/F)'};
metricTC  = {tc_dfof_all, tc_MAD_all};
type_masks = {HTP_cells, ~HTP_cells};
type_lbls  = {'HTP+', 'HTP-'};
type_fnames = {'HTPplus', 'HTPminus'};

allAxH = gobjects(2, nSize, 2);  % (row, size, type)

for iType = 1:2
    mask = type_masks{iType};
    figure('Position', [50 50 200*nSize 420]);
    axH = gobjects(2, nSize);

    for iRow = 1:2
        for iSize = 1:nSize
            axH(iRow, iSize) = subplot(2, nSize, (iRow-1)*nSize + iSize);
            hold on
            for iCon = 1:nCon
                tc = metricTC{iRow}(:, mask, iCon, iSize);
                m  = mean(tc, 2, 'omitnan');
                se = std(tc, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(tc(1,:))));
                fill([t_frames fliplr(t_frames)], [(m+se)' fliplr((m-se)')], ...
                    con_clrs(iCon,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
                plot(t_frames, m, 'Color', con_clrs(iCon,:), 'LineWidth', 1.2, ...
                    'DisplayName', sprintf('%g%%', cons(iCon)*100));
            end
            yl = ylim;
            plot([0 nOn], [yl(1) yl(1)], 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
            set(gca, 'TickDir', 'out', 'Box', 'off');
            xlabel('Frame');
            if iSize == 1, ylabel(rowYLbl{iRow}); end
            if iRow == 1,  title(sprintf('%g deg', sizes(iSize))); end
        end
    end

    sgtitle(sprintf('%s cells, stationary | contrasts overlaid', type_lbls{iType}));
    allAxH(:, :, iType) = axH;
    figh(iType) = gcf; %#ok<AGROW>
end

% Shared y limits per row across both figures
for iRow = 1:2
    axRow = squeeze(allAxH(iRow, :, :));  % nSize x 2
    yl = cell2mat(arrayfun(@(ax) ax.YLim, axRow(:), 'UniformOutput', false));
    ylShared = [min(yl(:,1)) max(yl(:,2))];
    set(axRow(:), 'YLim', ylShared);
end

for iType = 1:2
    saveas(figh(iType), fullfile(fnout, sprintf('timecourse_dfof_MAD_%s_contrastOverlay.pdf', type_fnames{iType})));
end
fprintf('Plots saved to %s\n', fnout);