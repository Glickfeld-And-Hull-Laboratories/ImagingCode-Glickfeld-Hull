dfof_stat_pre = reshape(pref_responses_stat_concat{pre},[3*nKeep_total,1]);
dfof_stat_post = reshape(pref_responses_stat_concat{post},[3*nKeep_total,1]);
dfof_loc_pre = reshape(pref_responses_loc_concat{pre},[3*nKeep_total,1]);
dfof_loc_post = reshape(pref_responses_loc_concat{post},[3*nKeep_total,1]);

dfof_col = vertcat(dfof_stat_pre,dfof_stat_post,dfof_loc_pre,dfof_loc_post);

cellID=1:nKeep_total;
cellID_col=repmat(cellID',12,1);

mouseIDcol = repmat(mouseID,12,1);

day1 = repelem(["pre" "post"],1,(3*nKeep_total))';
day=repmat(day1,2,1);


contrast1 = repelem(cons,nKeep_total)';
contrast=repmat(contrast1,4,1);

cell_type_col=repmat(red_concat,1,12)';

behStateCol = repelem(["stat" "loc"],1,(6*nKeep_total))';

clear dfof_stat_pre dfof_stat_post dfof_loc_pre dfof_loc_post cellID1 day1 contrast1



dfof_summary = table(mouseIDcol,cellID_col,cell_type_col,contrast,behStateCol,day,dfof_col, ...
    'VariableNames',{'mouseID' 'cellID' 'cellType' 'contrast' 'behState' 'day' 'dfof'});

%% example mixed model for all cells, stationary only
stat_dfof_summary_stat=dfof_summary((dfof_summary.behState=='stat'),:);
SST_stat_dfof = stat_dfof_summary_stat((stat_dfof_summary_stat.cellType==1),:);
Pyr_stat_dfof = stat_dfof_summary_stat((stat_dfof_summary_stat.cellType==0),:);

%for each cell type, construct a mixed model, with random effects to
%account for having multiple cells per mouse and mutiple measurements per
%cell, and fixed effects for contrast, pre/post day, and their interaction
lme_sst_stat= fitlme(SST_stat_dfof,'dfof~contrast*day+(1|mouseID)+(1|cellID:mouseID)');
lme_pyr_stat= fitlme(Pyr_stat_dfof,'dfof~contrast*day+(1|mouseID)+(1|cellID:mouseID)');

%% example mixed model for cells matched across behavioral states
%from the full dataframe, extract rows for running and stationary trials for SST and
%Pyr cells seperately. Use all SST and all Pyr cells.
all_cells = union(red_all,green_all);
matched_dfof_summary = dfof_summary(ismember(dfof_summary.cellID,all_cells),:);

SST_dfof = matched_dfof_summary((matched_dfof_summary.cellType==1),:);
%for each cell type, construct a mixed model, with random effects to
%account for having multiple cells per mouse and mutiple measurements per
%cell, and fixed effects for contrast, pre/post day, and their interaction
lme_sst= fitlme(SST_dfof,'dfof~behState+contrast+day+(behState:day)+(1|mouseID)+(1|cellID:mouseID)');
anova(lme_sst)

%simple main effects, i.e. ANOVA within each level of the behavioral state
%variable
SST_dfof_stat=SST_dfof(SST_dfof.behState=='stat',:);
SST_dfof_loc=SST_dfof(SST_dfof.behState=='loc',:);

lme_sst_stat= fitlme(SST_dfof_stat,'dfof~contrast*day+(1|mouseID)+(1|cellID:mouseID)');
anova(lme_sst_stat)

lme_sst_loc= fitlme(SST_dfof_loc,'dfof~contrast*day+(1|mouseID)+(1|cellID:mouseID)');
anova(lme_sst_loc)


%% regression for pupil
dfof_small_pre = reshape(pref_responses_stat_lowPupil_concat{pre},[3*nKeep_total,1]);
dfof_small_post = reshape(pref_responses_stat_lowPupil_concat{post},[3*nKeep_total,1]);
dfof_large_pre = reshape(pref_responses_stat_hiPupil_concat{pre},[3*nKeep_total,1]);
dfof_large_post = reshape(pref_responses_stat_hiPupil_concat{post},[3*nKeep_total,1]);

dfof_col = vertcat(dfof_small_pre,dfof_small_post,dfof_large_pre,dfof_large_post);

cellID=1:nKeep_total;
cellID_col=repmat(cellID',12,1);

mouseIDcol = repmat(mouseID,12,1);

DART1 = repelem(["pre" "post"],1,(3*nKeep_total))';
DART=repmat(DART1,2,1);


contrast1 = repelem(cons,nKeep_total)';
contrast=repmat(contrast1,4,1);

cell_type_col=repmat(red_concat,1,12)';

pupilSize = repelem(["small" "large"],1,(6*nKeep_total))';

clear dfof_small_pre dfof_small_post dfof_large_pre dfof_large_post cellID1 DART1 contrast1



dfof_summary = table(mouseIDcol,cellID_col,cell_type_col,contrast,pupilSize,DART,dfof_col, ...
    'VariableNames',{'mouseID' 'cellID' 'cellType' 'contrast' 'pupilSize' 'DART' 'dfof'});

%%

matched_dfof_summary = dfof_summary(ismember(dfof_summary.cellID,have_allPupil),:);

SST_dfof = matched_dfof_summary((matched_dfof_summary.cellType==1),:);

lme_sst= fitlme(SST_dfof,'dfof~pupilSize+contrast+DART+(pupilSize:DART)+(1|mouseID)+(1|cellID:mouseID)');
anova(lme_sst)

SST_dfof_small=SST_dfof(SST_dfof.pupilSize=='small',:);
SST_dfof_large=SST_dfof(SST_dfof.pupilSize=='large',:);

lme_sst_small= fitlme(SST_dfof_small,'dfof~contrast*DART+(1|mouseID)+(1|cellID:mouseID)');
anova(lme_sst_small)

lme_sst_large= fitlme(SST_dfof_large,'dfof~contrast*DART+(1|mouseID)+(1|cellID:mouseID)');
anova(lme_sst_large)


Pyr_dfof = matched_dfof_summary((matched_dfof_summary.cellType==0),:);

lme_pyr= fitlme(Pyr_dfof,'dfof~pupilSize+contrast+DART+(pupilSize:DART)+(1|mouseID)+(1|cellID:mouseID)');
anova(lme_pyr)

pyr_dfof_small=Pyr_dfof(Pyr_dfof.pupilSize=='small',:);
pyr_dfof_large=Pyr_dfof(Pyr_dfof.pupilSize=='large',:);

lme_pyr_small= fitlme(pyr_dfof_small,'dfof~contrast*DART+(1|mouseID)+(1|cellID:mouseID)');
anova(lme_pyr_small)

lme_pyr_large= fitlme(pyr_dfof_large,'dfof~contrast*DART+(1|mouseID)+(1|cellID:mouseID)');
anova(lme_pyr_large)


