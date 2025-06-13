%% Figure S2
%load in data and indices
%hard code specific information
cons=[.25 .5 1] %the visual stimulus contrasts used
pre=2; post=1; %indexing the two imaging days

load('FigureS2_data.mat')
nCells=size(pref_responses_stat_concat{1},1);
%% Figure S2A - timecourses
plotNeuralTimecourse(tc_trial_avrg_stat_concat, tc_trial_avrg_stat_concat, ...
    sst_ind, pyr_ind, ...
    'UseDashedLines', [false, true], ...
    'Colors1', {'k', 'c'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'c'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'SST', 'Pyr'}, ...
    'StimStart', 31);

%% Figure S2B - scatterplots and signed hypoetneuse 
scatter_signedHypDist(pref_responses_stat_concat, pre,post,sst_ind,pyr_ind)
%% Figure S2B - t-tests with Bonferroni correction
% pairwise ttests for dfof response at each contrast for SST cells
[sst_h1, sst_p1]= ttest(pref_responses_stat_concat{pre}(sst_ind,1),pref_responses_stat_concat{post}(sst_ind,1));
[sst_h2, sst_p2]= ttest(pref_responses_stat_concat{pre}(sst_ind,2),pref_responses_stat_concat{post}(sst_ind,2));
[sst_h3, sst_p3]= ttest(pref_responses_stat_concat{pre}(sst_ind,3),pref_responses_stat_concat{post}(sst_ind,3));

%corrected for three tests
sst_pvalues = [(sst_p1*3);(sst_p2*3);(sst_p3*3)];

% pairwise ttests for dfof response at each contrast for Pyr cells
[pyr_h1, pyr_p1]= ttest(pref_responses_stat_concat{pre}(pyr_ind,1),pref_responses_stat_concat{post}(pyr_ind,1));
[pyr_h2, pyr_p2]= ttest(pref_responses_stat_concat{pre}(pyr_ind,2),pref_responses_stat_concat{post}(pyr_ind,2));
[pyr_h3, pyr_p3]= ttest(pref_responses_stat_concat{pre}(pyr_ind,3),pref_responses_stat_concat{post}(pyr_ind,3));

%corrected for three tests
pyr_pvalues = [(pyr_p1*3);(pyr_p2*3);(pyr_p3*3)];

contrasts = cons';
table(contrasts,sst_pvalues,pyr_pvalues)
%% Figure S2C - contrast response plots
plotContrastResponse(pref_responses_stat_concat, pref_responses_stat_concat, ...
    sst_ind, pyr_ind, cons, ...
    'UseDashedLines', [false, true], ...  % Dashed lines for the right plot
    'Titles', {'SST', 'Pyr'}, ...
    'YLabel', 'dF/F','Colors', {'k', 'c'});

%% Figure S2C - ANOVA
w = table(categorical([1 1 1 2 2 2 ].'), categorical([1 2 3 1 2 3].'), 'VariableNames', {'DART', 'contrast'}); % within-design
rm_SST_stat = fitrm(SST_stat_dfof, 'd1c1-d2c3 ~ 1', 'WithinDesign', w);
ranova(rm_SST_stat, 'withinmodel', 'DART*contrast')

rm_Pyr_stat = fitrm(Pyr_stat_dfof, 'd1c1-d2c3 ~ 1', 'WithinDesign', w);
ranova(rm_Pyr_stat, 'withinmodel', 'DART*contrast')

%% Figure S2D - Contrast response for modulation index (DART vs. PEG)
load('FigureS2_D_H_data.mat')
DART=1;
PEG = 2;

%Compute modulation index
for id = 1:2
    %create rectified versions where negative values are set to 0
    combined_pref_responses_stat_rect{id} = combined_pref_responses_stat_concat{id};
    pref_responses_loc_rect{id} = combined_pref_responses_loc_concat{id};
    combined_pref_responses_stat_rect{id}(combined_pref_responses_stat_rect{id} < 0)=0;
    pref_responses_loc_rect{id}(pref_responses_loc_rect{id} < 0)=0;
end

norm_delta_stat=(combined_pref_responses_stat_rect{post}-combined_pref_responses_stat_rect{pre})./(combined_pref_responses_stat_rect{post}+combined_pref_responses_stat_rect{pre});
norm_delta_loc=(pref_responses_loc_rect{post}-pref_responses_loc_rect{pre})./(pref_responses_loc_rect{post}+pref_responses_loc_rect{pre});

%plot contrast response

conResp_green_avrg_stat = cell(1,2);
conResp_red_avrg_stat = cell(1,2); 
conResp_green_se_stat = cell(1,2); 
conResp_red_se_stat = cell(1,2); 

for iDrug = 1:2
    
    thisDrug = find(drug == iDrug);
    pyr_thisDrug = intersect(thisDrug,combined_pyr_inds);
    sst_thisDrug = intersect(thisDrug,combined_sst_inds);
        
    conResp_green_avrg_stat{iDrug}=mean(norm_delta_stat(pyr_thisDrug,:),1,'omitmissing');
    green_std=std(norm_delta_stat(pyr_thisDrug,:),1,'omitmissing');
    conResp_green_se_stat{iDrug}=green_std/sqrt(length(pyr_thisDrug));
    
    conResp_red_avrg_stat{iDrug}=mean(norm_delta_stat(sst_thisDrug,:),1,'omitmissing');
    red_std=std(norm_delta_stat(sst_thisDrug,:),1,'omitmissing');
    conResp_red_se_stat{iDrug}=red_std/sqrt(length(sst_thisDrug));
    
    clear green_std red_std
 
end

figure
subplot(1,2,1)
errorbar(cons,conResp_red_avrg_stat{PEG},conResp_red_se_stat{PEG},'color',"#00FFFF");
hold on
errorbar(cons,conResp_red_avrg_stat{DART},conResp_red_se_stat{DART},'b');
title('SST')
ylabel('Mod. index')
xlabel('contrast') 
xlim([0 1.1])
ylim([-.4 .2])
xticks([.25 .5 1])
set(gca, 'TickDir', 'out')
box off

subplot(1,2,2) 
errorbar(cons,conResp_green_avrg_stat{PEG},conResp_green_se_stat{PEG},'--','color',"#00FFFF");
hold on
errorbar(cons,conResp_green_avrg_stat{DART},conResp_green_se_stat{DART},'--b');
title('Pyr')
xlabel('contrast') 
xlim([0 1.1])
ylim([-.4 .2])
xticks([.25 .5 1])
box off

x0=5;
y0=5;
width=2.2;
height=1.6;
set(gcf,'units','inches','position',[x0,y0,width,height])

%% Figure S2D - ANOVA

%stationary model for cells not matched for running
WithinModel = table(categorical([1 2 3]'), 'VariableNames', {'contrast'}); % within-design
statDeltaModel_SST_unmatched = fitrm(norm_delta_stat_unMatchedSST, 'S25-S100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(statDeltaModel_SST_unmatched, 'withinmodel', 'contrast');
disp(ranovatbl);

statDeltaModel_PYR_unmatched = fitrm(norm_delta_stat_unMatchedPyr, 'S25-S100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(statDeltaModel_PYR_unmatched, 'withinmodel', 'contrast');
disp(ranovatbl);
%% Figure S2E - timecourses
plotNeuralTimecourse(tc_trial_avrg_loc_concat, tc_trial_avrg_loc_concat, ...
    sst_running_ind, pyr_running_ind, ...
    'UseDashedLines', [false, true], ...
    'Colors1', {'k', 'c'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'c'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'SST', 'Pyr'}, ...
    'StimStart', 31);

%% Figure S2F - scatterplots and signed hypoetneuse 
scatter_signedHypDist(pref_responses_loc_concat, pre,post,sst_running_ind,pyr_running_ind)
%% Figure S2F - t-tests with Bonferroni correction
% pairwise ttests for dfof response at each contrast for SST cells
[sst_h1, sst_p1]= ttest(pref_responses_loc_concat{pre}(sst_running_ind,1),pref_responses_loc_concat{post}(sst_running_ind,1));
[sst_h2, sst_p2]= ttest(pref_responses_loc_concat{pre}(sst_running_ind,2),pref_responses_loc_concat{post}(sst_running_ind,2));
[sst_h3, sst_p3]= ttest(pref_responses_loc_concat{pre}(sst_running_ind,3),pref_responses_loc_concat{post}(sst_running_ind,3));

%corrected for three tests
sst_pvalues = [(sst_p1*3);(sst_p2*3);(sst_p3*3)];

% pairwise ttests for dfof response at each contrast for Pyr cells
[pyr_h1, pyr_p1]= ttest(pref_responses_loc_concat{pre}(pyr_running_ind,1),pref_responses_loc_concat{post}(pyr_running_ind,1));
[pyr_h2, pyr_p2]= ttest(pref_responses_loc_concat{pre}(pyr_running_ind,2),pref_responses_loc_concat{post}(pyr_running_ind,2));
[pyr_h3, pyr_p3]= ttest(pref_responses_loc_concat{pre}(pyr_running_ind,3),pref_responses_loc_concat{post}(pyr_running_ind,3));

%corrected for three tests
pyr_pvalues = [(pyr_p1*3);(pyr_p2*3);(pyr_p3*3)];

contrasts = cons';
table(contrasts,sst_pvalues,pyr_pvalues)
%% Figure S2G - contrast response plots
plotContrastResponse(pref_responses_loc_concat, pref_responses_loc_concat, ...
    sst_running_ind, pyr_running_ind, cons, ...
    'UseDashedLines', [false, true], ...  % Dashed lines for the right plot
    'Titles', {'SST', 'Pyr'}, ...
    'YLabel', 'dF/F', 'Colors',{'k', 'c'});

%% Figure S2G - ANOVA
w = table(categorical([1 1 1 2 2 2 ].'), categorical([1 2 3 1 2 3].'), 'VariableNames', {'DART', 'contrast'}); % within-design

rm_SST_loc = fitrm(SST_matched_dfof, 'Rd1c1-Rd2c3 ~ 1', 'WithinDesign', w);
ranova(rm_SST_loc, 'withinmodel', 'DART*contrast')


rm_Pyr_loc = fitrm(Pyr_matched_dfof, 'Rd1c1-Rd2c3 ~ 1', 'WithinDesign', w);
ranova(rm_Pyr_loc, 'withinmodel', 'DART*contrast')

%% Figure S2H - Contrast response for modulation index (DART vs. PEG), running
conResp_green_avrg_loc = cell(1,2); 
conResp_red_avrg_loc = cell(1,2); 
conResp_green_se_loc = cell(1,2); 
conResp_red_se_loc = cell(1,2); 

for iDrug = 1:2
    
    thisDrug = find(drug == iDrug);
    pyr_thisDrug = intersect(thisDrug,combined_pyr_running_inds);
    sst_thisDrug = intersect(thisDrug,combined_sst_running_inds);
        
    conResp_green_avrg_loc{iDrug}=mean(norm_delta_loc(pyr_thisDrug,:),1,'omitmissing');
    green_std=std(norm_delta_loc(pyr_thisDrug,:),1,'omitmissing');
    conResp_green_se_loc{iDrug}=green_std/sqrt(length(pyr_thisDrug));
    
    conResp_red_avrg_loc{iDrug}=mean(norm_delta_loc(sst_thisDrug,:),1,'omitmissing');
    red_std=std(norm_delta_loc(sst_thisDrug,:),1,'omitmissing');
    conResp_red_se_loc{iDrug}=red_std/sqrt(length(sst_thisDrug));
    
    clear green_std red_std
 
end

figure
subplot(1,2,1) %for the second day
errorbar(cons,conResp_red_avrg_loc{PEG},conResp_red_se_loc{PEG},'color',"#00FFFF");
hold on
errorbar(cons,conResp_red_avrg_loc{DART},conResp_red_se_loc{DART},'b');
title('SST')
ylabel('Mod. index')
xlabel('contrast') 
xlim([0 1.2])
ylim([-.25 .3])
xticks([.25 .5 1])
% set(gca, 'TickDir', 'out')
box off

subplot(1,2,2) %for the first day
errorbar(cons,conResp_green_avrg_loc{PEG},conResp_green_se_loc{PEG},'--','color',"#00FFFF");
hold on
errorbar(cons,conResp_green_avrg_loc{DART},conResp_green_se_loc{DART},'--b');
title('Pyr')
xlabel('contrast') 
xlim([0 1.2])
ylim([-.25 .3])
xticks([.25 .5 1])
box off

x0=5;
y0=5;
width=2.2;
height=1.6;
set(gcf,'units','inches','position',[x0,y0,width,height])

%% Figure S2H - ANOVA

%running model for cells fully matched for running
WithinModel = table(categorical([1 2 3]'), 'VariableNames', {'contrast'}); % within-design
locDeltaModel_SST = fitrm(norm_delta_loc_SST, 'L25-L100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(locDeltaModel_SST, 'withinmodel', 'contrast');
disp(ranovatbl);

locDeltaModel_Pyr = fitrm(norm_delta_loc_Pyr, 'L25-L100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(locDeltaModel_Pyr, 'withinmodel', 'contrast');
disp(ranovatbl);