%% read in the dfof data tables
%one table for PEG, one for DART. These have hardcoded cell ID numbers so
%PEG has to go first. I have hardcoded to folders where these files are,
%based on the data I'm including as of 11/12/2023
clear all
PEG_data=readtable(fullfile('Z:\home\ACh\Analysis\2p_analysis\concat303_311_319_329_355\14-Nov-2023','dfof_PEG_output.csv'));
DRT_data = readtable(fullfile('Z:\home\ACh\Analysis\2p_analysis\concat138_142_163_171_178_190_294_307_323_333\14-Nov-2023','dfof_DRT_output.csv'));

dataTable=[DRT_data;PEG_data];
clear PEG_data DRT_data;

dataTable.drug = categorical(dataTable.drug);
dataTable.behState = categorical(dataTable.behState);
dataTable.cellID = categorical(dataTable.cellID);
dataTable.mouseID   = categorical(dataTable.mouseID);
dataTable.cellType = categorical(dataTable.cellType);
%% mixed effects model for dfof
lme = fitlme(dataTable,'dfof~cellType+contrast+behState+drug+day+(cellType*day)+(cellType*day*drug)+(1|mouseID)+(1|cellID:mouseID)')
%%
vars = ["dfof"];
factors = ["drug","behState","day"];
meanScoresByFactor = varfun(@nanmean, ...
                            dataTable, ...
                            "InputVariables",vars, ...
                            "GroupingVariables",factors)
%% post-hoc tests
%break down by cell type
SOMonly = dataTable(find(dataTable.cellType=='1'),:);
PYRonly = dataTable(find(dataTable.cellType=='0'),:);

lmeSOM = fitlme(SOMonly,'dfof~contrast+behState+drug+day+(day*drug)+(contrast*day*drug)+(behState*day*drug)+(1|mouseID)+(1|cellID:mouseID)');
lmePYR = fitlme(PYRonly,'dfof~contrast+behState+drug+day+(day*drug)+(contrast*day*drug)+(behState*day*drug)+(1|mouseID)+(1|cellID:mouseID)');
%ANOVA_SOM=anova(lmeSOM)
%ANOVA_PYR=anova(lmePYR)
lmeSOM

%% looking at SOM within drug condition seperately

SOM_DRT = SOMonly(find(SOMonly.drug=='DRT'),:);
SOM_PEG = SOMonly(find(SOMonly.drug=='PEG'),:);

lme_SOM_DRT = fitlme(SOM_DRT, 'dfof~day+(day*contrast)+(day*behState)+(1|mouseID)+(1|cellID:mouseID)');
lme_SOM_PEG = fitlme(SOM_PEG, 'dfof~day+(day*contrast)+(day*behState)+(1|mouseID)+(1|cellID:mouseID)');


%% post-hoc comaparisons for mixed model
%what is the mean df/f for each cell type on each day
vars = ["dfof"];
factors = ["behState","drug","day"];
meanScoresByFactor = varfun(@nanmean, ...
                            SOMonly, ...
                            "InputVariables",vars, ...
                            "GroupingVariables",factors)
%% SOM_stat DART
SOM_stat = SOM_DRT((SOM_DRT.behState=='stat'),:);
lme_SOM_stat= fitlme(SOM_stat, 'dfof~day+contrast+(day*contrast)+(1|mouseID)+(1|cellID:mouseID)');

SOM_loc = SOM_DRT((SOM_DRT.behState=='loc'),:);
lme_SOM_loc= fitlme(SOM_loc, 'dfof~day+contrast+(day*contrast)+(1|mouseID)+(1|cellID:mouseID)');

anova(lme_SOM_stat)
anova(lme_SOM_loc)

%% SOM stat PEG

SOM_stat_PEG = SOM_PEG((SOM_PEG.behState=='stat'),:);
lme_SOM_stat_PEG= fitlme(SOM_stat_PEG, 'dfof~day+contrast+(day*contrast)+(1|mouseID)+(1|cellID:mouseID)');

SOM_loc_PEG = SOM_PEG((SOM_PEG.behState=='loc'),:);
lme_SOM_loc_PEG= fitlme(SOM_loc_PEG, 'dfof~day+contrast+(day*contrast)+(1|mouseID)+(1|cellID:mouseID)');

anova(lme_SOM_stat_PEG)
anova(lme_SOM_loc_PEG)

%% looking at PYR within drug condition seperately

PYR_DRT = PYRonly(find(PYRonly.drug=='DRT'),:);
PYR_PEG = PYRonly(find(PYRonly.drug=='PEG'),:);

lme_PYR_DRT = fitlme(PYR_DRT, 'dfof~day+(day*contrast)+(day*behState)+(1|mouseID)+(1|cellID:mouseID)');
lme_PYR_PEG = fitlme(PYR_PEG, 'dfof~day+(day*contrast)+(day*behState)+(1|mouseID)+(1|cellID:mouseID)');

anova(lme_PYR_DRT)
anova(lme_PYR_PEG)
%% PYR_stat DART
PYR_stat = PYR_DRT((PYR_DRT.behState=='stat'),:);
lme_PYR_stat= fitlme(PYR_stat, 'dfof~day+contrast+(day*contrast)+(1|mouseID)+(1|cellID:mouseID)');

PYR_loc = PYR_DRT((PYR_DRT.behState=='loc'),:);
lme_PYR_loc= fitlme(PYR_loc, 'dfof~day+contrast+(day*contrast)+(1|mouseID)+(1|cellID:mouseID)');

anova(lme_PYR_stat)
anova(lme_PYR_loc)

%% PYR stat PEG

PYR_stat_PEG = PYR_PEG((PYR_PEG.behState=='stat'),:);
lme_PYR_stat_PEG= fitlme(PYR_stat_PEG, 'dfof~day+contrast+(day*contrast)+(1|mouseID)+(1|cellID:mouseID)');

PYR_loc_PEG = PYR_PEG((PYR_PEG.behState=='loc'),:);
lme_PYR_loc_PEG= fitlme(PYR_loc_PEG, 'dfof~day+contrast+(day*contrast)+(1|mouseID)+(1|cellID:mouseID)');

anova(lme_PYR_stat_PEG)
anova(lme_PYR_loc_PEG)



