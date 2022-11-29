function pooled_tuned_data = mean_tuned_cell_resp_by_trial(tuned_data)
    pooled_tuned_data = {};
    for bin = 1:4
        data_set = tuned_data(:,bin);
        pooled_tuned_data{bin} = mean_cell_resp_by_trial(data_set);
    end
end
