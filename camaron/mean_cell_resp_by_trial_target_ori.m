function pooled_target_data = mean_cell_resp_by_trial_target_ori(target_ori_data, bp_ind_sig_resp_stim_intersect)
% Finds mean cell activity across trials (for all targets) for all experiments
    pooled_target_data = cell(1,size(target_ori_data,2));
    for i = 1:size(target_ori_data,2)
        for k = 1:size(target_ori_data,1)
            cells = bp_ind_sig_resp_stim_intersect{k};
            data = target_ori_data{k,i};
            mean_data = mean(data(:,cells,:), 3, 'omitnan');
            pooled_target_data{i} = cat(2, pooled_target_data{i}, mean_data);
        end
    end
end