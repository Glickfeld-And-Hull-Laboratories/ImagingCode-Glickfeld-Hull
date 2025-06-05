%% Figure 3
%load in data and indices
%hard code specific information
cons=[.25 .5 1] %the visual stimulus contrasts used
nCon = length(cons);
pre=2; post=1; %indexing the two imaging days
load('Figure3_data.mat')
nCells=size(pref_allTrials_stat_concat{1},2);
load('Figure3_data_PEG.mat')
nCells_PEG = size(pref_allTrials_stat_concat_PEG{1},2);
%% Figure 3 A - Normalized difference scatters
% Calculate normalized differences
[norm_diff, bsln_std] = calculateNormalizedDifference(pref_allTrials_stat_concat, ...
    pref_allTrials_loc_concat, pre, post, nCon, nCells);
%plot the scatters for DART 
figure;
subplot(1,2,1)
boxchart(squeeze(norm_diff(1,:,sst_ind))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75],BoxEdgeColor=[0 0 0]);
hold on
scatter([1, 2, 3],squeeze(norm_diff(1,:,sst_ind))',20,[.26 .29 .33], 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)

xticklabels({'25','50','100'})
xlabel('Contrast(%)')
ylabel('Normalized difference')
ylim([-10 10])
title('SST')
hold off
set(gca,'TickDir','out')
box off

subplot(1,2,2)
boxchart(squeeze(norm_diff(1,:,pyr_ind))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75],BoxEdgeColor=[0 0 0]);
hold on
scatter([1, 2, 3],squeeze(norm_diff(1,:,pyr_ind))',20,[.26 .29 .33], 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)
boxchart(squeeze(norm_diff(1,:,pyr_ind))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75],BoxEdgeColor=[0 0 0]);
xticklabels({'25','50','100'})
xlabel('Contrast(%)')
ylim([-20 20])
title('Pyr')
hold off
set(gca,'TickDir','out')
box off
x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])

%% Figure 3 A -  Levene's test
% Extract data from the matrix for SST cells
data = squeeze(norm_diff(1, :, sst_ind));  % Extract data from the specified dimension

% Perform Levene's test
[p, stats] = vartestn(data','TestType','LeveneAbsolute');

% Display the p-value
disp(['Levene''s test for SST cells, p-value: ', num2str(p)]);

% Extract data from the matrix for SST cells
data = squeeze(norm_diff(1, :, pyr_ind));  % Extract data from the specified dimension

% Perform Levene's test
[p, stats] = vartestn(data','TestType','LeveneAbsolute');

% Display the p-value
disp(['Levene''s test for pyramidal cells, p-value: ', num2str(p)]);

%F-test for equality of variances, asking whether the higher contrast has
%greater variance than the lower contrast, for SST cells
[h,p1,ci,stats1] = vartest2(norm_diff(1,1,sst_ind),norm_diff(1,2,sst_ind),'Tail','left'); %25 vs 50
[h,p2,ci,stats2] = vartest2(norm_diff(1,1,sst_ind),norm_diff(1,3,sst_ind),'Tail','left'); %25 vs 100
[h,p3,ci,stats3] = vartest2(norm_diff(1,2,sst_ind),norm_diff(1,3,sst_ind),'Tail','left'); %50 vs 100


%F-test for equality of variances, asking whether the higher contrast has
%greater variance than the lower contrast, for Pyr cells
[h,p4,ci,stats4] = vartest2(norm_diff(1,1,pyr_ind),norm_diff(1,2,pyr_ind),'Tail','left'); %25 vs 50
[h,p5,ci,stats5] = vartest2(norm_diff(1,1,pyr_ind),norm_diff(1,3,pyr_ind),'Tail','left'); %25 vs 100
[h,p6,ci,stats6] = vartest2(norm_diff(1,2,pyr_ind),norm_diff(1,3,pyr_ind),'Tail','left'); %50 vs 100

format long 
[stats1.fstat, stats2.fstat, stats3.fstat; p1*3, p2*3,p3*3]
[stats4.fstat, stats5.fstat, stats6.fstat; p4*3, p5*3,p6*3]
format short

clear h p1 p2 p3 ci stats1 stats2 stats3 h p4 p5 p6 ci stats4 stats5 stats6

%% Figure 3B - Bar charts for the fraction of cells suppressed and facilitated SST cells

%make a subset of normalized difference for the SST cells only, then make
% find how many are facilitated or suppressed by more than 1 std from
% baseline
norm_diff_sst = norm_diff(:,:,sst_ind);
facil_sst=norm_diff_sst(:,:,:)>=1;
supp_sst=norm_diff_sst(:,:,:)<=-1;

N=length(sst_ind);
facil_table_stat = sum(facil_sst(1,:,:),3)/N;
supp_table_stat = sum(supp_sst(1,:,:),3)/N;

figure;
subplot(1,2,1)
bar([1,2,3],[supp_table_stat],'FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
xticklabels({'25','50','100'})
title('Suppressed')
ylim([0 .5])
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

subplot(1,2,2)
bar([1,2,3],[facil_table_stat],'FaceColor',"#A8518A",'EdgeColor', [1 1 1])
xticklabels({'25','50','100'})
title('Facilitated')
ylim([0 .5])
%ylabel(["Fraction HTP+ cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])

%% Figure 3B - Chi-square for the fraction of cells suppressed and facilitated SST cells
%compute chi squares for suppression
%25 vs 50
n1=supp_table_stat(1)*N;
n2=supp_table_stat(2)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%25 vs 10%
n1=supp_table_stat(1)*N;
n2=supp_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%50 vs 100
n1=supp_table_stat(2)*N;
n2=supp_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);

[chi2stat1, chi2stat2, chi2stat3; p1*3, p2*3,p3*3]

clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3

%compute chi squares for facilitation
%25 vs 50
n1=facil_table_stat(1)*N;
n2=facil_table_stat(2)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%25 vs 10%
n1=facil_table_stat(1)*N;
n2=facil_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%50 vs 100
n1=facil_table_stat(2)*N;
n2=facil_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);

[chi2stat1, chi2stat2, chi2stat3; p1*3, p2*3,p3*3]
clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3 n1 n2 x1 x2
%% Figure 3B - Bar charts for the fraction of cells suppressed and facilitated Pyramidal cells
norm_diff_pyr = norm_diff(:,:,pyr_ind);
facil_pyr=norm_diff_pyr(:,:,:)>=1;
supp_pyr=norm_diff_pyr(:,:,:)<=-1;

N=length(pyr_ind);
facil_table_stat = sum(facil_pyr(1,:,:),3)/N;
supp_table_stat = sum(supp_pyr(1,:,:),3)/N;

figure;
subplot(1,2,1)
bar([1,2,3],[supp_table_stat],'FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
xticklabels({'25','50','100'})
title('Suppressed')
ylim([0 .4])
ylabel(["Fraction Pyr cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

subplot(1,2,2)
bar([1,2,3],[facil_table_stat],'FaceColor',"#A8518A",'EdgeColor', [1 1 1])
xticklabels({'25','50','100'})
title('Facilitated')
ylim([0 .4])
%ylabel(["Fraction HTP+ cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])


%% Figure 3B - Chi-square for the fraction of cells suppressed and facilitated Pyramidal cells

%compute chi squares for suppression
%25 vs 50
n1=supp_table_stat(1)*N;
n2=supp_table_stat(2)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%25 vs 10%
n1=supp_table_stat(1)*N;
n2=supp_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%50 vs 100
n1=supp_table_stat(2)*N;
n2=supp_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);

[chi2stat1, chi2stat2, chi2stat3; p1*3, p2*3,p3*3]

clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3

%compute chi squares for facilitation
%25 vs 50
n1=facil_table_stat(1)*N;
n2=facil_table_stat(2)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%25 vs 100
n1=facil_table_stat(1)*N;
n2=facil_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%50 vs 100
n1=facil_table_stat(2)*N;
n2=facil_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);

[chi2stat1, chi2stat2, chi2stat3; p1*3, p2*3,p3*3]
clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3 n1 n2 x1 x2

%% Figure 3 C - Normalized difference scatters
% Calculate normalized differences
[norm_diff_PEG, bsln_std_PEG] = calculateNormalizedDifference(pref_allTrials_stat_concat_PEG, ...
    pref_allTrials_loc_concat_PEG, pre, post, nCon, nCells_PEG);
%plot the scatters for PEG 
figure;
subplot(1,2,1)
boxchart(squeeze(norm_diff_PEG(1,:,sst_ind_PEG))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75],BoxEdgeColor=[0 0 0]);
hold on
scatter([1, 2, 3],squeeze(norm_diff_PEG(1,:,sst_ind_PEG))',20,[.26 .29 .33], 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)

xticklabels({'25','50','100'})
xlabel('Contrast(%)')
ylabel('Normalized difference')
ylim([-10 10])
title('SST')
hold off
set(gca,'TickDir','out')
box off

subplot(1,2,2)
boxchart(squeeze(norm_diff_PEG(1,:,pyr_ind_PEG))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75],BoxEdgeColor=[0 0 0]);
hold on
scatter([1, 2, 3],squeeze(norm_diff_PEG(1,:,pyr_ind_PEG))',20,[.26 .29 .33], 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)
boxchart(squeeze(norm_diff_PEG(1,:,pyr_ind_PEG))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75],BoxEdgeColor=[0 0 0]);
xticklabels({'25','50','100'})
xlabel('Contrast(%)')
ylim([-20 20])
title('Pyr')
hold off
set(gca,'TickDir','out')
box off
x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])

%% Figure 3 C -  Levene's test
% Extract data from the matrix for SST cells
data = squeeze(norm_diff_PEG(1, :, sst_ind_PEG));  % Extract data from the specified dimension

% Perform Levene's test
[p, stats] = vartestn(data','TestType','LeveneAbsolute');

% Display the p-value
disp(['Levene''s test for SST cells, p-value: ', num2str(p)]);

% Extract data from the matrix for SST cells
data = squeeze(norm_diff_PEG(1, :, pyr_ind_PEG));  % Extract data from the specified dimension

% Perform Levene's test
[p, stats] = vartestn(data','TestType','LeveneAbsolute');

% Display the p-value
disp(['Levene''s test for pyramidal cells, p-value: ', num2str(p)]);

%F-test for equality of variances, asking whether the higher contrast has
%greater variance than the lower contrast, for SST cells
[h,p1,ci,stats1] = vartest2(norm_diff_PEG(1,1,sst_ind_PEG),norm_diff_PEG(1,2,sst_ind_PEG),'Tail','left'); %25 vs 50
[h,p2,ci,stats2] = vartest2(norm_diff_PEG(1,1,sst_ind_PEG),norm_diff_PEG(1,3,sst_ind_PEG),'Tail','left'); %25 vs 100
[h,p3,ci,stats3] = vartest2(norm_diff_PEG(1,2,sst_ind_PEG),norm_diff_PEG(1,3,sst_ind_PEG),'Tail','left'); %50 vs 100


%F-test for equality of variances, asking whether the higher contrast has
%greater variance than the lower contrast, for Pyr cells
[h,p4,ci,stats4] = vartest2(norm_diff_PEG(1,1,pyr_ind_PEG),norm_diff_PEG(1,2,pyr_ind_PEG),'Tail','left'); %25 vs 50
[h,p5,ci,stats5] = vartest2(norm_diff_PEG(1,1,pyr_ind_PEG),norm_diff_PEG(1,3,pyr_ind_PEG),'Tail','left'); %25 vs 100
[h,p6,ci,stats6] = vartest2(norm_diff_PEG(1,2,pyr_ind_PEG),norm_diff_PEG(1,3,pyr_ind_PEG),'Tail','left'); %50 vs 100

format long 
[stats1.fstat, stats2.fstat, stats3.fstat; p1*3, p2*3,p3*3]
[stats4.fstat, stats5.fstat, stats6.fstat; p4*3, p5*3,p6*3]
format short

clear h p1 p2 p3 ci stats1 stats2 stats3 h p4 p5 p6 ci stats4 stats5 stats6

%% Figure 3B - Bar charts for the fraction of cells suppressed and facilitated SST cells

%make a subset of normalized difference for the SST cells only, then make
% find how many are facilitated or suppressed by more than 1 std from
% baseline
norm_diff_sst_PEG = norm_diff_PEG(:,:,sst_ind_PEG);
facil_sst_PEG=norm_diff_sst_PEG(:,:,:)>=1;
supp_sst_PEG=norm_diff_sst_PEG(:,:,:)<=-1;

N=length(sst_ind_PEG);
facil_table_stat = sum(facil_sst_PEG(1,:,:),3)/N;
supp_table_stat = sum(supp_sst(1,:,:),3)/N;

figure;
subplot(1,2,1)
bar([1,2,3],[supp_table_stat],'FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
xticklabels({'25','50','100'})
title('Suppressed')
ylim([0 .5])
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

subplot(1,2,2)
bar([1,2,3],[facil_table_stat],'FaceColor',"#A8518A",'EdgeColor', [1 1 1])
xticklabels({'25','50','100'})
title('Facilitated')
ylim([0 .5])
%ylabel(["Fraction HTP+ cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])

%% Figure 3B - Chi-square for the fraction of cells suppressed and facilitated SST cells
%compute chi squares for suppression
%25 vs 50
n1=supp_table_stat(1)*N;
n2=supp_table_stat(2)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%25 vs 10%
n1=supp_table_stat(1)*N;
n2=supp_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%50 vs 100
n1=supp_table_stat(2)*N;
n2=supp_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);

[chi2stat1, chi2stat2, chi2stat3; p1*3, p2*3,p3*3]

clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3

%compute chi squares for facilitation
%25 vs 50
n1=facil_table_stat(1)*N;
n2=facil_table_stat(2)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%25 vs 10%
n1=facil_table_stat(1)*N;
n2=facil_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%50 vs 100
n1=facil_table_stat(2)*N;
n2=facil_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);

[chi2stat1, chi2stat2, chi2stat3; p1*3, p2*3,p3*3]
clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3 n1 n2 x1 x2
%% Figure 3D - Bar charts for the fraction of cells suppressed and facilitated Pyramidal cells
norm_diff_pyr = norm_diff_PEG(:,:,pyr_ind_PEG);
facil_pyr=norm_diff_pyr(:,:,:)>=1;
supp_pyr=norm_diff_pyr(:,:,:)<=-1;

N=length(pyr_ind_PEG);
facil_table_stat = sum(facil_pyr(1,:,:),3)/N;
supp_table_stat = sum(supp_pyr(1,:,:),3)/N;

figure;
subplot(1,2,1)
bar([1,2,3],[supp_table_stat],'FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
xticklabels({'25','50','100'})
title('Suppressed')
ylim([0 .4])
ylabel(["Fraction Pyr cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

subplot(1,2,2)
bar([1,2,3],[facil_table_stat],'FaceColor',"#A8518A",'EdgeColor', [1 1 1])
xticklabels({'25','50','100'})
title('Facilitated')
ylim([0 .4])
%ylabel(["Fraction HTP+ cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])


%% Figure 3D - Chi-square for the fraction of cells suppressed and facilitated Pyramidal cells

%compute chi squares for suppression
%25 vs 50
n1=supp_table_stat(1)*N;
n2=supp_table_stat(2)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%25 vs 10%
n1=supp_table_stat(1)*N;
n2=supp_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%50 vs 100
n1=supp_table_stat(2)*N;
n2=supp_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);

[chi2stat1, chi2stat2, chi2stat3; p1*3, p2*3,p3*3]

clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3

%compute chi squares for facilitation
%25 vs 50
n1=facil_table_stat(1)*N;
n2=facil_table_stat(2)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%25 vs 100
n1=facil_table_stat(1)*N;
n2=facil_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%50 vs 100
n1=facil_table_stat(2)*N;
n2=facil_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);

[chi2stat1, chi2stat2, chi2stat3; p1*3, p2*3,p3*3]
clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3 n1 n2 x1 x2