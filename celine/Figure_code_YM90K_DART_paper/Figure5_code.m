%% Figure 5
%load in data and indices
%hard code specific information
cons=[.25 .5 1] %the visual stimulus contrasts used
pre=2; post=1; %indexing the two imaging days
load('Figure5_data.mat')

%% Figure 5A - Timecourses for SST cells
plotNeuralTimecourse(tc_trial_avrg_stat_concat, tc_trial_avrg_loc_concat, ...
    sst_withRuning_ind, sst_withRuning_ind, ...
    'UseDashedLines', [false, false], ...
    'Colors1', {'k', 'b'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'b'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'Stationary', 'Running'}, ...
    'StimStart', 31);

%% Figure 5B and 5E - scatterplots and signed hypoetneuse 
scatter_signedHypDist(pref_responses_loc_concat, pre,post,sst_withRuning_ind,pyr_withRuning_ind)
%% Figure 5B and 5E- t-tests with Bonferroni correction

% pairwise ttests for dfof response at each contrast for SST cells
[sst_h1, sst_p1]= ttest(pref_responses_loc_concat{pre}(sst_withRuning_ind,1),pref_responses_loc_concat{post}(sst_withRuning_ind,1));
[sst_h2, sst_p2]= ttest(pref_responses_loc_concat{pre}(sst_withRuning_ind,2),pref_responses_loc_concat{post}(sst_withRuning_ind,2));
[sst_h3, sst_p3]= ttest(pref_responses_loc_concat{pre}(sst_withRuning_ind,3),pref_responses_loc_concat{post}(sst_withRuning_ind,3));

%corrected for three tests
sst_pvalues = [(sst_p1*3);(sst_p2*3);(sst_p3*3)];

% pairwise ttests for dfof response at each contrast for Pyr cells
[pyr_h1, pyr_p1]= ttest(pref_responses_loc_concat{pre}(pyr_withRuning_ind,1),pref_responses_loc_concat{post}(pyr_withRuning_ind,1));
[pyr_h2, pyr_p2]= ttest(pref_responses_loc_concat{pre}(pyr_withRuning_ind,2),pref_responses_loc_concat{post}(pyr_withRuning_ind,2));
[pyr_h3, pyr_p3]= ttest(pref_responses_loc_concat{pre}(pyr_withRuning_ind,3),pref_responses_loc_concat{post}(pyr_withRuning_ind,3));

%corrected for three tests
pyr_pvalues = [(pyr_p1*3);(pyr_p2*3);(pyr_p3*3)];

contrasts = cons';
table(contrasts,sst_pvalues,pyr_pvalues)

%% Figure 5C - contrast response plots for SST cells
plotContrastResponse(pref_responses_stat_concat, pref_responses_loc_concat, ...
    sst_withRuning_ind, sst_withRuning_ind, cons, ...
    'UseDashedLines', [false, false], ...  % Dashed lines for the right plot
    'Titles', {'Stationary', 'Running'}, ...
    'YLabel', 'dF/F');
%% Figure 5C - ANOVA
w = table(categorical([1 1 1 2 2 2 ].'), categorical([1 2 3 1 2 3].'), 'VariableNames', {'DART', 'contrast'}); % within-design

rm_SST_stat = fitrm(SST_matched_dfof, 'Sd1c1-Sd2c3 ~ 1', 'WithinDesign', w)
ranova(rm_SST_stat, 'withinmodel', 'DART*contrast')

rm_SST_loc = fitrm(SST_matched_dfof, 'Rd1c1-Rd2c3 ~ 1', 'WithinDesign', w)
ranova(rm_SST_loc, 'withinmodel', 'DART*contrast')
%% Figure 5D - Timecourses for Pyr cells
plotNeuralTimecourse(tc_trial_avrg_stat_concat, tc_trial_avrg_loc_concat, ...
    pyr_withRuning_ind, pyr_withRuning_ind, ...
    'UseDashedLines', [true, true], ...
    'Colors1', {'k', 'b'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'b'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'Stationary', 'Running'}, ...
    'StimStart', 31);

%% Figure 5F - contrast response plots for Pyr cells
plotContrastResponse(pref_responses_stat_concat, pref_responses_loc_concat, ...
    pyr_withRuning_ind, pyr_withRuning_ind, cons, ...
    'UseDashedLines', [true, true], ...  % Dashed lines for the right plot
    'Titles', {'Stationary', 'Running'}, ...
    'YLabel', 'dF/F');
%% Figure 5F - ANOVA
w = table(categorical([1 1 1 2 2 2 ].'), categorical([1 2 3 1 2 3].'), 'VariableNames', {'DART', 'contrast'}); % within-design
rm_Pyr_stat = fitrm(Pyr_matched_dfof, 'Sd1c1-Sd2c3 ~ 1', 'WithinDesign', w)
ranova(rm_Pyr_stat, 'withinmodel', 'DART*contrast')

rm_Pyr_loc = fitrm(Pyr_matched_dfof, 'Rd1c1-Rd2c3 ~ 1', 'WithinDesign', w)
ranova(rm_Pyr_loc, 'withinmodel', 'DART*contrast')

%% Figure 5G -  Bar charts for the fraction of cells suppressed and facilitated SST cells

%make a subset of normalized difference for the SST cells with running data only, then make
% find how many are facilitated or suppressed by more than 1 std from
% baseline
norm_diff_red = norm_diff(:,:,sst_withRuning_ind);

facil_sst=norm_diff_red(:,:,:)>=1;
supp_sst=norm_diff_red(:,:,:)<=-1;

N=length(sst_withRuning_ind);
facil_table_stat = sum(facil_sst(1,:,:),3)/N;
supp_table_stat = sum(supp_sst(1,:,:),3)/N;
facil_table_loc = sum(facil_sst(2,:,:),3)/N;
supp_table_loc = sum(supp_sst(2,:,:),3)/N;

figure;
subplot(1,2,1)
b=bar([1,2,3],[supp_table_stat; supp_table_loc],'grouped','FaceColor',"#00ffff",'EdgeColor', [1 1 1]);
b(1).FaceColor="#70D0F6"
b(2).FaceColor="#0C8ABB"
ylim([0 .6])
title('Suppressed')
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

subplot(1,2,2)
b=bar([1,2,3],[facil_table_stat; facil_table_loc],'FaceColor',"#a329cc",'EdgeColor', [1 1 1]);
b(1).FaceColor="#C983B1"
b(2).FaceColor="#883367"
xticklabels({'25','50','100'})
ylim([0 .6])
title('Facilitated')
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
%% Figure 5G - Chi-square for the fraction of cells suppressed and facilitated SST cells
%compute chi squares for suppression, stationary vs running
%25% contrast
% N previously set to be the number of SST+ cells. The same cells are
% included in stationary and running.
n1=supp_table_stat(1)*N;
n2=supp_table_loc(1)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%50
n1=supp_table_stat(2)*N;
n2=supp_table_loc(2)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%100
n1=supp_table_stat(3)*N;
n2=supp_table_loc(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);


%[chi2stat1, chi2stat2, chi2stat3; p1*3, p2*3,p3*3]
[chi2stat1, chi2stat2, chi2stat3; p1, p2,p3]

clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3

%compute chi squares for facilitation
%25
n1=facil_table_stat(1)*N;
n2=facil_table_loc(1)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%50
n1=facil_table_stat(2)*N;
n2=facil_table_loc(2)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%100
n1=facil_table_stat(3)*N;
n2=facil_table_loc(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);

[chi2stat1, chi2stat2, chi2stat3; p1, p2,p3]
clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3 n1 n2 x1 x2
%% Figure 5H - Bar charts for the fraction of cells suppressed and facilitated Pyr cells
norm_diff_green = norm_diff(:,:,pyr_withRuning_ind);

facil_pyr=norm_diff_green(:,:,:)>=1;
supp_pyr=norm_diff_green(:,:,:)<=-1;

N=length(pyr_withRuning_ind);
facil_table_stat = sum(facil_pyr(1,:,:),3)/N;
supp_table_stat = sum(supp_pyr(1,:,:),3)/N;
facil_table_loc = sum(facil_pyr(2,:,:),3)/N;
supp_table_loc = sum(supp_pyr(2,:,:),3)/N;

figure;
subplot(1,2,1)
b=bar([1,2,3],[supp_table_stat; supp_table_loc],'grouped','FaceColor',"#00ffff",'EdgeColor', [1 1 1]);
b(1).FaceColor="#70D0F6"
b(2).FaceColor="#0C8ABB"
ylim([0 .6])
title('Suppressed')
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

subplot(1,2,2)
b=bar([1,2,3],[facil_table_stat; facil_table_loc],'FaceColor',"#a329cc",'EdgeColor', [1 1 1]);
b(1).FaceColor="#C983B1"
b(2).FaceColor="#883367"
xticklabels({'25','50','100'})
ylim([0 .6])
title('Facilitated')
ylabel(["Fraction HTP+ cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])

%% Figure 5H - Chi-square for the fraction of cells suppressed and facilitated Pyr cells

%25% contrast

n1=int32(supp_table_stat(1)*N);
n2=int32(supp_table_loc(1)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%50
n1=int32(supp_table_stat(2)*N);
n2=int32(supp_table_loc(2)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%100
n1=int32(supp_table_stat(3)*N);
n2=int32(supp_table_loc(3)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);



%[chi2stat1, chi2stat2, chi2stat3; p1*3, p2*3,p3*3]
[chi2stat1, chi2stat2, chi2stat3; p1, p2,p3]

clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3

%compute chi squares for facilitation
%25
n1=int32(supp_table_stat(1)*N);
n2=int32(facil_table_loc(1)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%50
n1=int32(facil_table_stat(2)*N);
n2=int32(facil_table_loc(2)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%100
n1=int32(facil_table_stat(3)*N);
n2=int32(facil_table_loc(3)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);

[chi2stat1, chi2stat2, chi2stat3; p1, p2,p3]
clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3 n1 n2 x1 x2
