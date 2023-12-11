%% read in the normDiff data tables
%one table for PEG, one for DART. These have hardcoded cell ID numbers so
%PEG has to go first. I have hardcoded to folders where these files are,
%based on the data I'm including as of 11/12/2023
clear all
PEG_data=readtable(fullfile('Z:\home\ACh\Analysis\2p_analysis\concat303_311_319_329_355_359\04-Dec-2023','DART_effect_output.csv'));
DRT_data = readtable(fullfile('Z:\home\ACh\Analysis\2p_analysis\YM90K_DART\30-Nov-2023','DART_effect_output.csv'));

dataTable=[DRT_data;PEG_data];
clear PEG_data DRT_data;

dataTable.drug = categorical(dataTable.drug);
%% mixed effects model for norm diff
lme = fitlme(dataTable,'normDiff~cellType+contrast+behState+drug+(cellType*drug)+(1|mouseID)+(1|cellID:mouseID)');
lme
%%
vars = ["normDiff"];
factors = ["cellType","contrast","behState"];
meanScoresByFactor = varfun(@nanmean, ...
                            dataTable, ...
                            "InputVariables",vars, ...
                            "GroupingVariables",factors)
%% post-hoc tests
%break down by cell type
SOMonly = dataTable(find(dataTable.cellType==1),:);
PYRonly = dataTable(find(dataTable.cellType==0),:);

lmeSOM = fitlme(SOMonly,'normDiff~contrast+behState+drug+(behState*drug)+(contrast*drug)+(1|mouseID)+(1|cellID:mouseID)');
lmePYR = fitlme(PYRonly,'normDiff~contrast+behState+drug+(behState*drug)+(contrast*drug)+(1|mouseID)+(1|cellID:mouseID)');
%ANOVA_SOM=anova(lmeSOM)
%ANOVA_PYR=anova(lmePYR)
lmeSOM

%%

SOM_DRT = SOMonly(find(SOMonly.drug=='DRT'),:);
SOM_PEG = SOMonly(find(SOMonly.drug=='PEG'),:);

lme_SOM_DRT = fitlme(SOM_DRT, 'normDiff~contrast+behState+(1|mouseID)+(1|cellID:mouseID)')
lme_SOM_PEG = fitlme(SOM_PEG, 'normDiff~contrast+behState+(1|mouseID)+(1|cellID:mouseID)')


