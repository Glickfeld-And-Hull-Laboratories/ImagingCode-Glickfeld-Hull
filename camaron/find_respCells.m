function [resp_ind, p] = find_respCells(data, base_win, resp_win)
% Given a data structure of shape (Time X Cells X Trials) perform a ttest
% comparing activity for each cell across the base_win and resp_win time
% windows.

% data, double: dF/F data, (Time X Cells X Trials)
% base_win, double: range of timepoints to use as baseline activity
% reps_win, double: range of timepoints to use as response activity

data_base = squeeze(mean(data(base_win,:,:),1, "omitnan"));
data_resp = squeeze(mean(data(resp_win,:,:),1, "omitnan"));
[h, p] = ttest(data_resp, data_base, 'tail', 'right','dim',2);
resp_ind = find(h);

end
