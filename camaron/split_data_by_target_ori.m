function [target_ori_data, target_ori_data_adapt, target_ori_data_contol] = split_data_by_target_ori(dFoF, ori_trial_ind, adapt_trial_ind, control_trial_ind)
% Finds mean cell activity across trials (for all targets) for all experiments
   for i = 1:size(ori_trial_ind, 2)
        trials_this_ori = ori_trial_ind{i};
        adapt_trials_this_ori = intersect(adapt_trial_ind, trials_this_ori);
        control_trials_this_ori = intersect(control_trial_ind, trials_this_ori);
        target_ori_data{i} = dFoF(:,:,trials_this_ori);
        target_ori_data_adapt{i} = dFoF(:,:,adapt_trials_this_ori);
        target_ori_data_contol{i} = dFoF(:,:, control_trials_this_ori);
    end
end