% prepare data for ANOVA
all_cells = union(red_all,green_all);

dfof_stat = horzcat(pref_responses_stat_concat{pre},pref_responses_stat_concat{post});
dfof_loc = horzcat(pref_responses_loc_concat{pre},pref_responses_loc_concat{post});

cellID=(1:nKeep_total)';
cell_type_col=categorical(red_concat)';

dfof_full = horzcat(dfof_stat,dfof_loc);
dfof_table = array2table(dfof_full,'VariableNames',{'Sd1c1','Sd1c2','Sd1c3', ...
    'Sd2c1','Sd2c2','Sd2c3','Rd1c1','Rd1c2','Rd1c3','Rd2c1','Rd2c2','Rd2c3'});

labels_table =table(mouseID,cellID,cell_type_col,'VariableNames',{'mouseID' 'cellID' 'cellType'});
stat_dfof_summary_full = [labels_table,dfof_table];
matched_dfof_summary = stat_dfof_summary_full(ismember(stat_dfof_summary_full.cellID,all_cells),:);


SST_matched_dfof = matched_dfof_summary((matched_dfof_summary.cellType=='1'),:);
Pyr_matched_dfof = matched_dfof_summary((matched_dfof_summary.cellType=='0'),:);
clear dfof_stat_table cell_type_col cellID dfof_stat_matched

%%
% run an ANOVA on the full model for each cell type
DART = categorical([1 1 1 2 2 2 1 1 1 2 2 2].');
contrast = categorical([1 2 3 1 2 3 1 2 3 1 2 3].');
behState = categorical([1 1 1 1 1 1 2 2 2 2 2 2].');

% Create a table with your categorical variables
w = table(DART, contrast, behState, ...
    'VariableNames', {'DART', 'contrast', 'behState'});% within-design

rm_SST_matched = fitrm(SST_matched_dfof, 'Sd1c1-Rd2c3 ~ 1', 'WithinDesign', w);
ranova(rm_SST_matched, 'withinmodel', 'behState*DART*contrast')

rm_Pyr_matched = fitrm(Pyr_matched_dfof, 'Sd1c1-Rd2c3 ~ 1', 'WithinDesign', w);
ranova(rm_Pyr_matched, 'withinmodel', 'behState*DART*contrast')

%%
w = table(categorical([1 1 1 2 2 2 ].'), categorical([1 2 3 1 2 3].'), 'VariableNames', {'DART', 'contrast'}); % within-design

rm_SST_stat = fitrm(SST_matched_dfof, 'Sd1c1-Sd2c3 ~ 1', 'WithinDesign', w)
ranova(rm_SST_stat, 'withinmodel', 'DART*contrast')

rm_SST_loc = fitrm(SST_matched_dfof, 'Rd1c1-Rd2c3 ~ 1', 'WithinDesign', w)
ranova(rm_SST_loc, 'withinmodel', 'DART*contrast')

rm_Pyr_stat = fitrm(Pyr_matched_dfof, 'Sd1c1-Sd2c3 ~ 1', 'WithinDesign', w)
ranova(rm_Pyr_stat, 'withinmodel', 'DART*contrast')

rm_Pyr_loc = fitrm(Pyr_matched_dfof, 'Rd1c1-Rd2c3 ~ 1', 'WithinDesign', w)
ranova(rm_Pyr_loc, 'withinmodel', 'DART*contrast')




