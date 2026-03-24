function plotNormDiffBySize_scatter(norm_diff_concat, red_ind_concat, green_ind_concat, exp_idx, nSess, cons, targetSize, fnout)
% plotNormDiffBySize_scatter  Like plotNormDiffBySize but with individual
%   data points colored by experiment, and no line connecting grand mean.
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
            nCells      = size(data_cells, 2);

            subplot(nCon, 2, (iCon-1)*2 + iCellType);
            hold on;

            xpos = 1:nSize;

            % Individual points colored by experiment
            for iSess = 1:nSess
                sess_mask  = exp_cells == iSess;
                data_sess  = data_cells(:, sess_mask);
                nCellsSess = sum(sess_mask);
                for iSize = 1:nSize
                    jitter = 0.15 * (rand(1, nCellsSess) - 0.5);
                    scatter(xpos(iSize) + jitter, data_sess(iSize, :), 8, ...
                        colors(iSess,:), 'filled', 'MarkerFaceAlpha', 0.4);
                end
            end

            % Grand mean +/- SE on top (no connecting line)
            overall_mean = mean(data_cells, 2, 'omitnan');
            overall_se   = std(data_cells, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(data_cells), 2));
            errorbar(xpos, overall_mean, overall_se, 'ko', ...
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
            xlim([0.5, nSize+0.5]);
            set(gca, 'TickDir', 'out', 'XTick', xpos, 'XTickLabel', arrayfun(@num2str, targetSize, 'UniformOutput', false));
            box off;
        end
    end

    sgtitle(state_names{iState});
    if nargin >= 8 && ~isempty(fnout)
        print(fullfile(fnout, sprintf('normDiff_bySize_scatter_%s.pdf', state_names{iState})), '-dpdf', '-bestfit');
    end
end
end