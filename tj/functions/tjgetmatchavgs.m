function[a, b, c, d, e, f] = tjgetmatchavgs(tc, matches, nOn, nOff, nTrials)

%this function will take in a timecourse from phase reverse tj experiment
%animal, matching cell indices, number of frames, and number of trials ->
%will return the avg, std, and stderr for cell population avg across trials
%and the average of each cell across all trials

tc_match = tc.cellTCs_match{1,2}(:,matches);
nCells = size(tc_match,2);
data_trial = permute(reshape(tc_match,[nOff + nOn nTrials nCells]),[1 3 2]);
data_f = mean(data_trial(nOff-15:nOff,:,:),1);
data_dfof = (data_trial-data_f)./data_f;


pop_avg_per_trial = mean(data_dfof(nOff+1:nOn,:,:),1, 'omitnan');
pop_std_per_trial = std(pop_avg_per_trial, 'omitnan');
pop_se_per_trial = pop_std_per_trial./sqrt(size(pop_avg_per_trial,2));

pop_avg_per_trial = mean(pop_avg_per_trial,2, 'omitnan');

pop_avg_per_trial = squeeze(pop_avg_per_trial);
pop_std_per_trial = squeeze(pop_std_per_trial);
pop_se_per_trial = squeeze(pop_se_per_trial);


cell_avg_trial = mean(data_dfof(nOff+1:nOn,:,:),1, 'omitnan');
cell_std_trial = std(cell_avg_trial,[],3, 'omitnan');
cell_se_trial = cell_std_trial./sqrt(size(cell_avg_trial,3));

cell_avg_trial = mean(cell_avg_trial,3, 'omitnan');

cell_avg_trial = squeeze(cell_avg_trial);
cell_std_trial = squeeze(cell_std_trial);
cell_se_trial = squeeze(cell_se_trial);


a = pop_avg_per_trial;
b = pop_std_per_trial;
c = pop_se_per_trial;
d = cell_avg_trial;
e = cell_std_trial;
f = cell_se_trial;


end

