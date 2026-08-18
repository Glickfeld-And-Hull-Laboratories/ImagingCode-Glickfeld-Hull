% Create experiment index for each cell
exp_idx = [];
for iSess = 1:nSess
    exp_idx = [exp_idx, iSess * ones(1, nKeep_concat(iSess))];
end

% Define colors for experiments
colors = lines(nSess);

cell_indices = {red_ind_concat, green_ind_concat};
cell_names = {'HTP+', 'HTP-'};

% Version 1: Delta retinotopic distance
delta_ret_distance = ret_distance_keep_concat{post}-ret_distance_keep_concat{pre};

for iCon = 1:nCon
    data = squeeze(norm_diff_concat(1,iCon,:,:));
    num_sizes = size(data,1);
    keep = all(data <= 5, 1);
    data_clean = data(:, keep);
    delta_ret_distance_clean = delta_ret_distance(keep);
    exp_idx_clean = exp_idx(keep);
    
    figure('Position', [100, 100, 1400, 600]);
    for iCellType = 1:2
        these_cells = cell_indices{iCellType};
        cell_mask = ismember(find(keep), these_cells);
        
        for i = 1:num_sizes
            subplot(2, num_sizes, (iCellType-1)*num_sizes + i)
            hold on
            
            % Plot points color coded by experiment
            for iSess = 1:nSess
                sess_mask = exp_idx_clean == iSess & cell_mask;
                scatter(delta_ret_distance_clean(sess_mask), data_clean(i,sess_mask), 25, colors(iSess,:), 'filled', 'MarkerFaceAlpha', 0.6)
            end
            
            mdl = fitlm(delta_ret_distance_clean(cell_mask), data_clean(i,cell_mask));
            r_squared = mdl.Rsquared.Ordinary;
            p_value = mdl.Coefficients.pValue(2);
            x_fit = linspace(min(delta_ret_distance_clean), max(delta_ret_distance_clean), 100)';
            [y_fit, y_ci] = predict(mdl, x_fit);
            fill([x_fit; flipud(x_fit)], [y_ci(:,1); flipud(y_ci(:,2))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
            plot(x_fit, y_fit, 'k-', 'LineWidth', 1.5)
            
            if i == 1
                ylabel([cell_names{iCellType} ' Norm Diff'])
            end
            if iCellType == 2
                xlabel('Delta Ret Distance')
            end
            title(sprintf('Size %d, R² = %.2f, p = %.2f', i, r_squared, p_value))
            ylim([min(data_clean(:)), max(data_clean(:))])
            xlim([-15, 15])
            set(gca, 'TickDir', 'out')
            grid off
            box off
        end
    end
    sgtitle(sprintf('Contrast %.2f', cons(iCon)));
    print(sprintf('normDiff_vsRet_byCellType_con%.2f', cons(iCon)), '-dpdf', '-bestfit')
end

% Version 2: Noise correlation
noise_corr = noiseCorr_concat{pre}(1,:);

for iCon = 1:nCon
    data = squeeze(norm_diff_concat(1,iCon,:,:));
    num_sizes = size(data,1);
    keep = all(data <= 5, 1);
    data_clean = data(:, keep);
    noise_corr_clean = noise_corr(keep);
    exp_idx_clean = exp_idx(keep);
    
    figure('Position', [100, 100, 1400, 600]);
    for iCellType = 1:2
        these_cells = cell_indices{iCellType};
        cell_mask = ismember(find(keep), these_cells);
        
        for i = 1:num_sizes
            subplot(2, num_sizes, (iCellType-1)*num_sizes + i)
            hold on
            
            % Plot points color coded by experiment
            for iSess = 1:nSess
                sess_mask = exp_idx_clean == iSess & cell_mask;
                scatter(noise_corr_clean(sess_mask), data_clean(i,sess_mask), 25, colors(iSess,:), 'filled', 'MarkerFaceAlpha', 0.6)
            end
            
            mdl = fitlm(noise_corr_clean(cell_mask), data_clean(i,cell_mask));
            r_squared = mdl.Rsquared.Ordinary;
            p_value = mdl.Coefficients.pValue(2);
            x_fit = linspace(min(noise_corr_clean), max(noise_corr_clean), 100)';
            [y_fit, y_ci] = predict(mdl, x_fit);
            fill([x_fit; flipud(x_fit)], [y_ci(:,1); flipud(y_ci(:,2))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
            plot(x_fit, y_fit, 'k-', 'LineWidth', 1.5)
            
            if i == 1
                ylabel([cell_names{iCellType} ' Norm Diff'])
            end
            if iCellType == 2
                xlabel('Noise Correlation')
            end
            title(sprintf('Size %d, R² = %.2f, p = %.2f', i, r_squared, p_value))
            ylim([min(data_clean(:)), max(data_clean(:))])
            set(gca, 'TickDir', 'out')
            grid off
            box off
        end
    end
    sgtitle(sprintf('Contrast %.2f', cons(iCon)));
    print(sprintf('normDiff_vsNoiseCorr_byCellType_con%.2f', cons(iCon)), '-dpdf', '-bestfit')
end
%% running
% Create experiment index for each cell
exp_idx = [];
for iSess = 1:nSess
    exp_idx = [exp_idx, iSess * ones(1, nKeep_concat(iSess))];
end

% Define colors for experiments
colors = lines(nSess);

cell_indices = {red_ind_concat, green_ind_concat};
cell_names = {'HTP+', 'HTP-'};

% Version 1: Locomotion - Delta retinotopic distance
delta_ret_distance = ret_distance_keep_concat{post}-ret_distance_keep_concat{pre};

for iCon = 1:nCon
    data = squeeze(norm_diff_concat(2,iCon,:,:));
    num_sizes = size(data,1);
    keep = all(data <= 5, 1);
    data_clean = data(:, keep);
    delta_ret_distance_clean = delta_ret_distance(keep);
    exp_idx_clean = exp_idx(keep);
    
    figure('Position', [100, 100, 1400, 600]);
    for iCellType = 1:2
        these_cells = cell_indices{iCellType};
        cell_mask = ismember(find(keep), these_cells);
        
        for i = 1:num_sizes
            subplot(2, num_sizes, (iCellType-1)*num_sizes + i)
            hold on
            
            % Plot points color coded by experiment
            for iSess = 1:nSess
                sess_mask = exp_idx_clean == iSess & cell_mask;
                scatter(delta_ret_distance_clean(sess_mask), data_clean(i,sess_mask), 25, colors(iSess,:), 'filled', 'MarkerFaceAlpha', 0.6)
            end
            
            mdl = fitlm(delta_ret_distance_clean(cell_mask), data_clean(i,cell_mask));
            r_squared = mdl.Rsquared.Ordinary;
            p_value = mdl.Coefficients.pValue(2);
            x_fit = linspace(min(delta_ret_distance_clean), max(delta_ret_distance_clean), 100)';
            [y_fit, y_ci] = predict(mdl, x_fit);
            fill([x_fit; flipud(x_fit)], [y_ci(:,1); flipud(y_ci(:,2))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
            plot(x_fit, y_fit, 'k-', 'LineWidth', 1.5)
            
            if i == 1
                ylabel([cell_names{iCellType} ' Norm Diff'])
            end
            if iCellType == 2
                xlabel('Delta Ret Distance')
            end
            title(sprintf('Size %d, R² = %.2f, p = %.2f', i, r_squared, p_value))
            ylim([min(data_clean(:)), max(data_clean(:))])
            xlim([-15, 15])
            set(gca, 'TickDir', 'out')
            grid off
            box off
        end
    end
    sgtitle(sprintf('Locomotion - Contrast %.2f', cons(iCon)));
    print(sprintf('normDiff_vsRet_byCellType_loc_con%.2f', cons(iCon)), '-dpdf', '-bestfit')
end

% Version 2: Locomotion - Noise correlation
noise_corr = noiseCorr_concat{pre}(1,:);

for iCon = 1:nCon
    data = squeeze(norm_diff_concat(2,iCon,:,:));
    num_sizes = size(data,1);
    keep = all(data <= 5, 1);
    data_clean = data(:, keep);
    noise_corr_clean = noise_corr(keep);
    exp_idx_clean = exp_idx(keep);
    
    figure('Position', [100, 100, 1400, 600]);
    for iCellType = 1:2
        these_cells = cell_indices{iCellType};
        cell_mask = ismember(find(keep), these_cells);
        
        for i = 1:num_sizes
            subplot(2, num_sizes, (iCellType-1)*num_sizes + i)
            hold on
            
            % Plot points color coded by experiment
            for iSess = 1:nSess
                sess_mask = exp_idx_clean == iSess & cell_mask;
                scatter(noise_corr_clean(sess_mask), data_clean(i,sess_mask), 25, colors(iSess,:), 'filled', 'MarkerFaceAlpha', 0.6)
            end
            
            mdl = fitlm(noise_corr_clean(cell_mask), data_clean(i,cell_mask));
            r_squared = mdl.Rsquared.Ordinary;
            p_value = mdl.Coefficients.pValue(2);
            x_fit = linspace(min(noise_corr_clean), max(noise_corr_clean), 100)';
            [y_fit, y_ci] = predict(mdl, x_fit);
            fill([x_fit; flipud(x_fit)], [y_ci(:,1); flipud(y_ci(:,2))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
            plot(x_fit, y_fit, 'k-', 'LineWidth', 1.5)
            
            if i == 1
                ylabel([cell_names{iCellType} ' Norm Diff'])
            end
            if iCellType == 2
                xlabel('Noise Correlation')
            end
            title(sprintf('Size %d, R² = %.2f, p = %.2f', i, r_squared, p_value))
            ylim([min(data_clean(:)), max(data_clean(:))])
            set(gca, 'TickDir', 'out')
            grid off
            box off
        end
    end
    sgtitle(sprintf('Locomotion - Contrast %.2f', cons(iCon)));
    print(sprintf('normDiff_vsNoiseCorr_byCellType_loc_con%.2f', cons(iCon)), '-dpdf', '-bestfit')
end