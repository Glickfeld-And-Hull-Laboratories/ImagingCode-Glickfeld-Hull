function plotNormDiffBySize(norm_diff_concat, red_ind_concat, green_ind_concat, exp_idx, nSess, cons, targetSize, fnout)
% plotNormDiffBySize  Plot mean +/- SE of norm_diff vs stimulus size
%   Generates two figures (stationary and locomotion), each with nCon rows
%   x 2 columns. Left column = HTP+ (red), right column = HTP- (green).
%   Overall mean +/- SE in black; per-experiment mean +/- SE in color.
%
%   INPUTS
%   norm_diff_concat  [2 x nCon x nSize x nNeurons]  state 1=stat, 2=loc
%   red_ind_concat    indices of HTP+ cells
%   green_ind_concat  indices of HTP- cells
%   exp_idx           [1 x nNeurons] session index per neuron
%   nSess             number of sessions
%   cons              [1 x nCon] contrast values
%   targetSize        [1 x nSize] stimulus sizes
%   fnout             output directory for saving PDFs

nCon  = size(norm_diff_concat, 2);
nSize = size(norm_diff_concat, 3);
colors = lines(nSess);

cell_indices = {red_ind_concat, green_ind_concat};
cell_names   = {'HTP+', 'HTP-'};
state_names  = {'Stationary', 'Locomotion'};
state_idx    = {1, 2};

for iState = 1:2
    figure('Position', [100, 100, 700, 300*nCon]);

    for iCon = 1:nCon
        data = squeeze(norm_diff_concat(state_idx{iState}, iCon, :, :)); % [nSize x nNeurons]

        for iCellType = 1:2
            these_cells = cell_indices{iCellType};
            data_cells  = data(:, these_cells); % [nSize x nCellsOfType]
            exp_cells   = exp_idx(these_cells);

            sp = subplot(nCon, 2, (iCon-1)*2 + iCellType);
            hold on;

            % Per-experiment mean +/- SE (drawn first, underneath overall)
            jitter_range = 0.03 * range(targetSize);
            jitter_step  = jitter_range / max(nSess - 1, 1);
            jitter_offsets = linspace(-jitter_range/2, jitter_range/2, nSess);
            for iSess = 1:nSess
                sess_mask = exp_cells == iSess;
                data_sess = data_cells(:, sess_mask);
                sess_mean = mean(data_sess, 2, 'omitnan');
                sess_se   = std(data_sess, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(data_sess), 2));
                errorbar(targetSize + jitter_offsets(iSess), sess_mean, sess_se, 'o', ...
                    'Color', colors(iSess,:), 'LineWidth', 0.75, ...
                    'MarkerFaceColor', colors(iSess,:), 'MarkerSize', 4, 'CapSize', 0);
            end

            % Overall mean +/- SE on top
            overall_mean = mean(data_cells, 2, 'omitnan');
            overall_se   = std(data_cells, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(data_cells), 2));
            errorbar(targetSize, overall_mean, overall_se, 'k-o', ...
                'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'CapSize', 4);

            yline(0, '--', 'Color', [.6 .6 .6], 'LineWidth', 1);

            if iCellType == 1
                ylabel(sprintf('Con %g\nNorm Diff', cons(iCon)));
            end
            if iCon == 1
                title(cell_names{iCellType});
            end
            if iCon == nCon
                xlabel('Stimulus size');
            end
            set(gca, 'TickDir', 'out', 'XTick', targetSize);
            box off;
        end
    end

    sgtitle(state_names{iState});
    if nargin >= 8 && ~isempty(fnout)
        print(fullfile(fnout, sprintf('normDiff_bySize_%s.pdf', state_names{iState})), '-dpdf', '-bestfit');
    end
end
end