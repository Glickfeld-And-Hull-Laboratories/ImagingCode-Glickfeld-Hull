%% Figure 2
%load in data and indices
%hard code specific information
cons=[.25 .5 1] %the visual stimulus contrasts used
pre=2; post=1; %indexing the two imaging days

load('Figure2_data.mat')
nCells=size(pref_responses_stat_concat{1},1);
%% Figure 2A - timecourses
plotNeuralTimecourse(tc_trial_avrg_stat_concat, tc_trial_avrg_stat_concat, ...
    sst_ind, pyr_ind, ...
    'UseDashedLines', [false, true], ...
    'Colors1', {'k', 'b'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'b'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'SST', 'Pyr'}, ...
    'StimStart', 31);

%% Figure 2B - scatterplots and signed hypoetneuse 
scatter_signedHypDist(pref_responses_stat_concat, pre,post,sst_ind,pyr_ind)
%% Figure 2B - t-tests with Bonferroni correction
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
%% Figure 2C - contrast response plots
plotContrastResponse(pref_responses_stat_concat, pref_responses_stat_concat, ...
    sst_ind, pyr_ind, cons, ...
    'UseDashedLines', [false, true], ...  % Dashed lines for the right plot
    'Titles', {'SST', 'Pyr'}, ...
    'YLabel', 'dF/F');

%% Figure 2C - ANOVA
w = table(categorical([1 1 1 2 2 2 ].'), categorical([1 2 3 1 2 3].'), 'VariableNames', {'DART', 'contrast'}); % within-design
rm_SST_stat = fitrm(SST_stat_dfof, 'd1c1-d2c3 ~ 1', 'WithinDesign', w);
ranova(rm_SST_stat, 'withinmodel', 'DART*contrast')

rm_Pyr_stat = fitrm(Pyr_stat_dfof, 'd1c1-d2c3 ~ 1', 'WithinDesign', w);
ranova(rm_Pyr_stat, 'withinmodel', 'DART*contrast')
