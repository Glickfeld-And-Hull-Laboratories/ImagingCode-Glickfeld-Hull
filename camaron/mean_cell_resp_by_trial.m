function pooled_data = mean_cell_resp_by_trial(dFoF)
% Finds mean cell activity across trials for all experiments
    pooled_data = [];
    for i = 1:length(dFoF)
        data = dFoF{i};
        mean_data = mean(data, 3, 'omitnan');
        pooled_data = cat(2, pooled_data, mean_data);
    end
end
