%% Load data...
% Save/load appropriate vars - pool_2AFC_2P_data.m



%% Percent responsive

prcnt_adapt_sig_b = nSigAdapt_b/nAllCells * 100;
prcnt_stim_sig_b = nSigStim_control_b/nAllCells * 100;

prcnt_adapt_sig_p = nSigAdapt_p/nAllCells * 100;
prcnt_stim_sig_p = nSigStim_control_p/nAllCells * 100;


for i = 1:length(ori_bins)
    prcnt_adapt_sig_tuned_b(i) = nSigTunedAdapt_b(i)/nAllCells * 100;
    prcnt_stim_sig_tuned_b(i) = nSigTunedStim_control_b(i)/nAllCells * 100;

    prcnt_adapt_sig_tuned_p(i) = nSigTunedAdapt_p(i)/nAllCells * 100;
    prcnt_stim_sig_tuned_p(i) = nSigTunedStim_control_p(i)/nAllCells * 100;
end

paired_adapt_sig_tuned = [prcnt_adapt_sig_tuned_b; prcnt_adapt_sig_tuned_p]';
paired_stim_sig_tuned = [prcnt_stim_sig_tuned_b; prcnt_stim_sig_tuned_p]';


figure()
bar(paired_adapt_sig_tuned)
ylabel('% Responsive')
xticks([1 2 3 4])
xticklabels({'0', '45', '90', '135'})
title("% cells responding to adaptor")

figure()
bar(paired_stim_sig_tuned)
ylabel('% Responsive')
xticks([1 2 3 4])
xticklabels({'0', '45', '90', '135'})
title("% cells repsonding to target")