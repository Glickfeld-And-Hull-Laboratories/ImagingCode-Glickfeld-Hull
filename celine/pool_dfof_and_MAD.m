% Pool population-level dF/F and MAD outputs across sessions and plot.
% Loads singleday_extraction_MAD.mat per session.
% Stacks sessions and plots mean +/- SE across sessions, direction-averaged.
%
% Session list txt format:
%   Line 1 (uncommented): datasheet .m file name
%   Remaining uncommented lines: session indices (day_id) to pool
clear all
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

% First pass: load all sessions and collect union of cons/sizes/dirs
nSess       = length(day_ids);
sessionData = cell(nSess, 1);
allCons     = [];
allSizes    = [];
allDirs     = [];
ExperimentFolder = '';

for i = 1:nSess
    iDay             = day_ids(i);
    mouse            = expt(iDay).mouse;
    expDate          = expt(iDay).date;
    ExperimentFolder = expt(iDay).exptType;
    runs             = eval(['expt(iDay).' cell2mat(dataStructLabels) '_runs']);
    fnday            = fullfile(rc.analysis, ExperimentFolder, mouse, expDate, runs{1});

    matMAD = fullfile(fnday, 'singleday_extraction_MAD.mat');
    if ~exist(matMAD, 'file')
        fprintf('Session %d: missing singleday_extraction_MAD.mat, skipping\n', iDay);
        continue
    end

    d = load(matMAD, 'pop_MAD_tc_HTPminus',   'pop_MAD_resp_HTPminus', ...
                     'pop_dfof_tc_HTPminus',   'pop_dfof_resp_HTPminus', ...
                     'pop_MAD_tc_HTPplus',     'pop_MAD_resp_HTPplus', ...
                     'pop_dfof_tc_HTPplus',    'pop_dfof_resp_HTPplus', ...
                     'cons', 'sizes', 'dirs', 'nOn', 'nOff');
    sessionData{i} = d;

    allCons  = union(allCons,  d.cons);
    allSizes = union(allSizes, d.sizes);
    allDirs  = union(allDirs,  d.dirs);
    fprintf('Session %d (%s %s): loaded\n', iDay, mouse, expDate);
end

sessionData = sessionData(~cellfun(@isempty, sessionData));
nSess   = length(sessionData);
nCon    = length(allCons);
nSize   = length(allSizes);
nDir    = length(allDirs);
nFrames = size(sessionData{1}.pop_MAD_tc_HTPminus, 1);
nOn     = sessionData{1}.nOn;
nOff    = sessionData{1}.nOff;
cons    = allCons;
sizes   = allSizes;
dirs    = allDirs;

% Second pass: align to common cons/sizes/dirs and stack into nFrames x nSess x nCon x nSize x nDir
pop_MAD_tc_HTPminus_sess    = nan(nFrames, nSess, nCon, nSize, nDir);
pop_MAD_resp_HTPminus_sess  = nan(nSess, nCon, nSize, nDir);
pop_dfof_tc_HTPminus_sess   = nan(nFrames, nSess, nCon, nSize, nDir);
pop_dfof_resp_HTPminus_sess = nan(nSess, nCon, nSize, nDir);
pop_MAD_tc_HTPplus_sess     = nan(nFrames, nSess, nCon, nSize, nDir);
pop_MAD_resp_HTPplus_sess   = nan(nSess, nCon, nSize, nDir);
pop_dfof_tc_HTPplus_sess    = nan(nFrames, nSess, nCon, nSize, nDir);
pop_dfof_resp_HTPplus_sess  = nan(nSess, nCon, nSize, nDir);

for i = 1:nSess
    d = sessionData{i};
    [~, conIdx]  = ismember(allCons,  d.cons);
    [~, sizeIdx] = ismember(allSizes, d.sizes);
    [~, dirIdx]  = ismember(allDirs,  d.dirs);
    for iCon = 1:nCon
        for iSize = 1:nSize
            for iDir = 1:nDir
                if conIdx(iCon) > 0 && sizeIdx(iSize) > 0 && dirIdx(iDir) > 0
                    cIdx = conIdx(iCon);
                    sIdx = sizeIdx(iSize);
                    dIdx = dirIdx(iDir);
                    pop_MAD_tc_HTPminus_sess(:, i, iCon, iSize, iDir)   = d.pop_MAD_tc_HTPminus(:, cIdx, sIdx, dIdx);
                    pop_MAD_resp_HTPminus_sess(i, iCon, iSize, iDir)    = d.pop_MAD_resp_HTPminus(cIdx, sIdx, dIdx);
                    pop_dfof_tc_HTPminus_sess(:, i, iCon, iSize, iDir)  = d.pop_dfof_tc_HTPminus(:, cIdx, sIdx, dIdx);
                    pop_dfof_resp_HTPminus_sess(i, iCon, iSize, iDir)   = d.pop_dfof_resp_HTPminus(cIdx, sIdx, dIdx);
                    pop_MAD_tc_HTPplus_sess(:, i, iCon, iSize, iDir)    = d.pop_MAD_tc_HTPplus(:, cIdx, sIdx, dIdx);
                    pop_MAD_resp_HTPplus_sess(i, iCon, iSize, iDir)     = d.pop_MAD_resp_HTPplus(cIdx, sIdx, dIdx);
                    pop_dfof_tc_HTPplus_sess(:, i, iCon, iSize, iDir)   = d.pop_dfof_tc_HTPplus(:, cIdx, sIdx, dIdx);
                    pop_dfof_resp_HTPplus_sess(i, iCon, iSize, iDir)    = d.pop_dfof_resp_HTPplus(cIdx, sIdx, dIdx);
                end
            end
        end
    end
end

clear sessionData
fprintf('\nTotal sessions pooled: %d\n', nSess);

% Save
fnpool = fullfile(rc.analysis, ExperimentFolder);
save(fullfile(fnpool, 'pooled_dfof_and_MAD.mat'), ...
    'pop_MAD_tc_HTPminus_sess',   'pop_MAD_resp_HTPminus_sess', ...
    'pop_dfof_tc_HTPminus_sess',  'pop_dfof_resp_HTPminus_sess', ...
    'pop_MAD_tc_HTPplus_sess',    'pop_MAD_resp_HTPplus_sess', ...
    'pop_dfof_tc_HTPplus_sess',   'pop_dfof_resp_HTPplus_sess', ...
    'cons', 'sizes', 'dirs', 'nOn', 'nOff', '-v7.3');
fprintf('Saved to %s\n', fullfile(fnpool, 'pooled_dfof_and_MAD.mat'));

% Output folder
d_str      = char(string(datetime('today')));
sess_title = strjoin(arrayfun(@num2str, day_ids, 'UniformOutput', false), '_');
fnout      = fullfile(rc.analysis, 'pooled_dfof_and_MAD', ['sess' sess_title], d_str);
mkdir(fnout);

%% Plot: mean +/- SE across sessions, direction-averaged, contrasts overlaid, rows = [dF/F; MAD], cols = sizes

t_frames    = (1:(nOn+nOff)) - nOff/2;
con_clrs    = [0.75 0.88 1; 0.2 0.5 0.85; 0.05 0.15 0.5];
rowYLbl     = {'\DeltaF/F', 'MAD (frac. change)'};
type_sess   = {{pop_dfof_tc_HTPminus_sess, pop_MAD_tc_HTPminus_sess}, ...
               {pop_dfof_tc_HTPplus_sess,  pop_MAD_tc_HTPplus_sess}};
type_lbls   = {'HTP-', 'HTP+'};
type_fnames = {'HTPminus', 'HTPplus'};

for iType = 1:2
    metricSess = type_sess{iType};
    figure('Position', [50 50 200*nSize 420]);
    for iRow = 1:2
        for iSize = 1:nSize
            subplot(2, nSize, (iRow-1)*nSize + iSize);
            hold on
            for iCon = 1:nCon
                tc = squeeze(nanmean(metricSess{iRow}(:, :, iCon, iSize, :), 5));
                m  = mean(tc, 2, 'omitnan');
                se = std(tc, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(tc(1,:))));
                fill([t_frames fliplr(t_frames)], [(m+se)' fliplr((m-se)')], ...
                    con_clrs(iCon,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
                plot(t_frames, m, 'Color', con_clrs(iCon,:), 'LineWidth', 1.2, ...
                    'DisplayName', sprintf('%g%%', cons(iCon)*100));
            end
            set(gca, 'TickDir', 'out', 'Box', 'off');
            xlim([t_frames(1) - nOff/4, t_frames(end) + nOff/4]);
            xlabel('Frame');
            if iSize == 1, ylabel(rowYLbl{iRow}); end
            if iRow == 1,  title(sprintf('%g deg', sizes(iSize))); end
        end
    end
    for iRow = 1:2
        allMin = inf;
        allMax = -inf;
        for iSize = 1:nSize
            for iCon = 1:nCon
                tc = squeeze(nanmean(metricSess{iRow}(:, :, iCon, iSize, :), 5));
                m  = mean(tc, 2, 'omitnan');
                se = std(tc, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(tc(1,:))));
                allMin = min(allMin, min(m - se, [], 'omitnan'));
                allMax = max(allMax, max(m + se, [], 'omitnan'));
            end
        end
        pad = 0.05 * (allMax - allMin);
        ylShared = [allMin - pad, allMax + pad];
        for iSize = 1:nSize
            ax = subplot(2, nSize, (iRow-1)*nSize + iSize);
            set(ax, 'YLim', ylShared);
            plot(ax, [0 nOn], [ylShared(1) ylShared(1)], ...
                'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
        end
    end
    sgtitle(sprintf('%s population, stationary | n = %d sessions', type_lbls{iType}, nSess));
    saveas(gcf, fullfile(fnout, sprintf('timecourse_dfof_MAD_%s_contrastOverlay.pdf', type_fnames{iType})));
end
fprintf('Plots saved to %s\n', fnout);

%% Plot: contrast and size tuning for dF/F and MAD, rows = metric, cols = contrast/size tuning

type_resp  = {{pop_MAD_resp_HTPminus_sess, pop_dfof_resp_HTPminus_sess}, ...
              {pop_MAD_resp_HTPplus_sess,  pop_dfof_resp_HTPplus_sess}};

for iType = 1:2
    mad_resp  = squeeze(nanmean(type_resp{iType}{1}, 4));  % nSess x nCon x nSize, dir-averaged
    dfof_resp = squeeze(nanmean(type_resp{iType}{2}, 4));

    mad_con  = squeeze(nanmean(mad_resp,  3));  % nSess x nCon, size-averaged
    dfof_con = squeeze(nanmean(dfof_resp, 3));
    mad_sz   = squeeze(nanmean(mad_resp,  2));  % nSess x nSize, contrast-averaged
    dfof_sz  = squeeze(nanmean(dfof_resp, 2));

    % Tuning: rows = metric, col1 = contrast tuning (line per size), col2 = size tuning (line per contrast)
    sz_clrs  = gray(nSize + 2);
    sz_clrs  = sz_clrs(2:end-1, :);  % avoid pure white/black ends

    mad_resp_dir  = squeeze(nanmean(type_resp{iType}{1}, 4));   % nSess x nCon x nSize
    dfof_resp_dir = squeeze(nanmean(type_resp{iType}{2}, 4));
    metrics2      = {dfof_resp_dir, mad_resp_dir};
    ylbls_tune    = {'\DeltaF/F', 'MAD (frac. change)'};

    figure('Position', [50 50 500 600]);
    slopes = nan(2, nSize);  % rows = [dfof, mad], cols = sizes
    for iRow = 1:2
        dat3 = metrics2{iRow};  % nSess x nCon x nSize

        % contrast tuning, one line per size
        subplot(3, 2, (iRow-1)*2 + 1); hold on
        for iSize = 1:nSize
            m  = nanmean(squeeze(dat3(:, :, iSize)), 1);   % 1 x nCon
            se = nanstd(squeeze(dat3(:, :, iSize)), 0, 1) ./ sqrt(sum(~isnan(squeeze(dat3(:, :, iSize))), 1));
            xv = cons(:)';
            fill([xv fliplr(xv)], [(m+se) fliplr(m-se)], sz_clrs(iSize,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            plot(xv, m, '-o', 'Color', sz_clrs(iSize,:), 'MarkerFaceColor', sz_clrs(iSize,:), ...
                'MarkerSize', 5, 'LineWidth', 1.2, 'DisplayName', sprintf('%g deg', sizes(iSize)));
            p = polyfit(xv, m, 1);
            slopes(iRow, iSize) = p(1) / nanmean(m);
        end
        set(gca, 'TickDir', 'out', 'Box', 'off');
        xlabel('Contrast'); ylabel(ylbls_tune{iRow});
        if iRow == 1, legend('Location', 'northwest', 'Box', 'off'); end

        % size tuning, one line per contrast
        subplot(3, 2, (iRow-1)*2 + 2); hold on
        for iCon = 1:nCon
            m  = nanmean(squeeze(dat3(:, iCon, :)), 1);
            se = nanstd(squeeze(dat3(:, iCon, :)), 0, 1) ./ sqrt(sum(~isnan(squeeze(dat3(:, iCon, :))), 1));
            xv = sizes(:)';
            fill([xv fliplr(xv)], [(m+se) fliplr(m-se)], con_clrs(iCon,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            plot(xv, m, '-o', 'Color', con_clrs(iCon,:), 'MarkerFaceColor', con_clrs(iCon,:), ...
                'MarkerSize', 5, 'LineWidth', 1.2, 'DisplayName', sprintf('%g%%', cons(iCon)*100));
        end
        set(gca, 'TickDir', 'out', 'Box', 'off');
        xlabel('Size (deg)');
        if iRow == 1, legend('Location', 'northwest', 'Box', 'off'); end
    end

    % contrast slope vs size
    subplot(3, 2, 5); hold on
    plot(sizes(:)', slopes(1,:), 'k-o', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'LineWidth', 1.2, 'DisplayName', '\DeltaF/F');
    plot(sizes(:)', slopes(2,:), '-o', 'Color', [0.2 0.5 0.85], 'MarkerFaceColor', [0.2 0.5 0.85], 'MarkerSize', 5, 'LineWidth', 1.2, 'DisplayName', 'MAD');
    yline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    set(gca, 'TickDir', 'out', 'Box', 'off');
    xlabel('Size (deg)'); ylabel('Contrast slope / mean response');
    legend('Location', 'best', 'Box', 'off');
    sgtitle(sprintf('%s population, stationary | n = %d sessions', type_lbls{iType}, nSess));
    saveas(gcf, fullfile(fnout, sprintf('tuning_dfof_MAD_%s.pdf', type_fnames{iType})));
end