% Compute SSI per cell per day per contrast
SSI = cell(1, nd);
for id = 1:nd
    SSI{id} = nan(nKeep_total, nCon);
    for iCon = 1:nCon
        resp             = squeeze(pref_responses_stat_concat{id}(:, iCon, :)); % [nCells x nSize]
        peak_resp        = max(resp, [], 2);
        largest_resp     = resp(:, end);
        SSI{id}(:, iCon) = (peak_resp - largest_resp) ./ peak_resp;
    end
end

% Mask outliers
for id = 1:nd
    SSI{id}(abs(SSI{id}) > 2) = nan;
end

% Plot 1: SSI vs ret distance, per contrast, day, and cell type
day_labels  = {'Pre', 'Post'};
cell_groups = {red_goodfit, green_goodfit};
cell_names  = {'HTP+', 'HTP-'};

for iCon = 1:nCon
    figure;
    for id = 1:nd
        for iCT = 1:2
            subplot(nd, 2, (id-1)*2 + iCT);
            idx = cell_groups{iCT};
            scatter(ret_distance_retino_concat{id}(idx), SSI{id}(idx, iCon), 20, 'filled');
            xlabel('RF distance (deg)'); ylabel('SSI');
            title([day_labels{id} ' - ' cell_names{iCT} ' (n=' num2str(length(idx)) ')']);
            set(gca, 'TickDir', 'out'); box off;
        end
    end
    sgtitle(['SSI vs RF distance - contrast ' num2str(targetCon(iCon))]);
    saveas(gcf, fullfile(fnout, ['SSI_vs_retDistance_con' num2str(targetCon(iCon)) '.pdf']));
end

% Plot 2: change in SSI vs change in ret distance, per contrast and cell type
delta_retDist = ret_distance_retino_concat{post} - ret_distance_retino_concat{pre};

figure;
for iCon = 1:nCon
    delta_SSI = SSI{post}(:, iCon) - SSI{pre}(:, iCon);
    for iCT = 1:2
        subplot(nCon, 2, (iCon-1)*2 + iCT);
        idx = cell_groups{iCT};
        scatter(delta_retDist(idx), delta_SSI(idx), 20, 'filled');
        yline(0, 'k--'); xline(0, 'k--');
        xlabel('\Delta RF distance (deg)'); ylabel('\Delta SSI');
        title([cell_names{iCT} ' - con ' num2str(targetCon(iCon)) ' (n=' num2str(length(idx)) ')']);
        set(gca, 'TickDir', 'out'); box off;
    end
end
sgtitle('\Delta SSI vs \Delta RF distance');
saveas(gcf, fullfile(fnout, 'deltaSSI_vs_deltaRetDistance_allCon.pdf'));