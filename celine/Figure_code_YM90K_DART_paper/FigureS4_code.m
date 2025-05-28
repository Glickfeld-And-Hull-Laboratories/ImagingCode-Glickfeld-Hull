%% Figure S4
%load in data and indices
%hard code specific information
cons=[.25 .5 1] %the visual stimulus contrasts used
pre=2; post=1; %indexing the two imaging days
nMice=10;
load('FigureS4_data.mat')
%% Figure S4A - Example pupil trace
load('FigureS4_pupilExample.mat');%pupil diamter for both imaging sessions concatenated
statPupilThreshold=prctile(i2107_pupil,50);
%plot only the first session
plot(i2107_pupil(1:400)*2); hline(statPupilThreshold*2); xlabel('trials');ylabel('pupil diam'); 
vline(165,'-m')
vline(350,'-g')

%% Figure S4B - pupil diameter

pupilMean_overall =mean(pupilMeans_clean,2,'omitmissing');
pupilSTD_overall =std(pupilMeans_clean,[],2,'omitmissing');
pupilSE_overall = pupilSTD_overall./sqrt(nMice);   

figure
plot(pupilMeans_clean(:,:),'-','Color',[.5 .5 .5])
xlim([.75 3.25])
box off
set(gca,'TickDir','out')
xticks([1 2 3])
xticklabels({'Small pupil','Large pupil','Running'})
ylabel('Pupil diameter (mm)')
hold on
ylim([0.2 .8])
scatter([1 2 3],pupilMeans_clean(:,:),10, "MarkerEdgeColor","none","MarkerFaceColor",[0 0 0],"MarkerFaceAlpha",.25)
errorbar(pupilMean_overall,pupilSE_overall,'.-k','MarkerSize',15)


x0=5;
y0=5;
width=2;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])

% confirm that there's no difference in size across days within the "small"
% pupil trials
[h p] = ttest(pupilMeans_concat(pre,2,:),pupilMeans_concat(post,2,:))
% large pupil trials are significantly larger
[h p] = ttest(pupilMeans_concat(pre,1,:),pupilMeans_concat(post,1,:))


%confirm that pupil is significantly larger in large trials, averaging over
%days
[h p] = ttest(pupilMeans_clean(1,:),pupilMeans_clean(2,:))
%test whether pupuil is different between large trials and running
[h p] = ttest(pupilMeans_clean(2,:),pupilMeans_clean(3,:))

%% Figure S4C - Locomotion by pupil size

pupilMotorMean_overall =mean(motorByPupil_clean,2,'omitmissing');
pupilMotorSTD_overall =std(motorByPupil_clean,[],2,'omitmissing');
pupilMotorSE_overall = pupilMotorSTD_overall./sqrt(nMice);   

figure
plot(motorByPupil_clean(:,:),'-','Color',[.5 .5 .5])
xlim([.75 2.25])
box off
set(gca,'TickDir','out')
xticks([1 2])
xticklabels({'Small pupil','Large pupil'})
ylabel('Wheel speed (cm/s)')
hold on
ylim([.1 .4])
scatter([1 2],motorByPupil_clean(:,:),10, "MarkerEdgeColor","none","MarkerFaceColor",[0 0 0],"MarkerFaceAlpha",.25)
errorbar(pupilMotorMean_overall,pupilMotorSE_overall,'.-k','MarkerSize',15)


x0=5;
y0=5;
width=1;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])


%look at how motor activity differs with pupil size control day
[h p] = ttest(motorByPupil_clean(1,:),motorByPupil_clean(2,:))

%% Figure S4D - timecourses for SST cells
plotNeuralTimecourse(tc_trial_avrg_stat_smallPupil_concat, tc_trial_avrg_stat_largePupil_concat, ...
    sst_inds, sst_inds, ...
    'UseDashedLines', [false, false], ...
    'Colors1', {'k', 'b'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'b'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'Small pupil', 'Large pupil'}, ...
    'StimStart', 31);

%% Figure S4E - Contrast response for SST cells
plotContrastResponse(pref_responses_stat_smallPupil_concat, pref_responses_stat_largePupil_concat, ...
    sst_inds, sst_inds, cons, ...
    'UseDashedLines', [false, false], ...  % Dashed lines for the right plot
    'Titles', {'Small pupil', 'Large pupil'}, ...
    'YLabel', 'dF/F');
%% Figure S4E - ANOVA for SST cells
w = table(categorical([1 1 1 2 2 2 ].'), categorical([1 2 3 1 2 3].'), 'VariableNames', {'DART', 'contrast'}); % within-design

rm_SST_small = fitrm(SST_matched_dfof_pupil, 'Sd1c1-Sd2c3 ~ 1', 'WithinDesign', w);
ranova(rm_SST_small, 'withinmodel', 'DART*contrast')

rm_SST_large = fitrm(SST_matched_dfof_pupil, 'Ld1c1-Ld2c3 ~ 1', 'WithinDesign', w);
ranova(rm_SST_large, 'withinmodel', 'DART*contrast')

%% Figure S4F - timecourses for Pyr cells
plotNeuralTimecourse(tc_trial_avrg_stat_smallPupil_concat, tc_trial_avrg_stat_largePupil_concat, ...
    pyr_inds, pyr_inds, ...
    'UseDashedLines', [true, true], ...
    'Colors1', {'k', 'b'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'b'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'Small pupil', 'Large pupil'}, ...
    'StimStart', 31);
%% Figure S4G - Contrast response for Pyr cells
plotContrastResponse(pref_responses_stat_smallPupil_concat, pref_responses_stat_largePupil_concat, ...
    pyr_inds, pyr_inds, cons, ...
    'UseDashedLines', [true, true], ...  % Dashed lines for the right plot
    'Titles', {'Small pupil', 'Large pupil'}, ...
    'YLabel', 'dF/F');
%% Figure S4G - ANOVA for Pyr  cells
w = table(categorical([1 1 1 2 2 2 ].'), categorical([1 2 3 1 2 3].'), 'VariableNames', {'DART', 'contrast'}); % within-design
rm_Pyr_small = fitrm(Pyr_matched_dfof_pupil, 'Sd1c1-Sd2c3 ~ 1', 'WithinDesign', w);
ranova(rm_Pyr_small, 'withinmodel', 'DART*contrast')

rm_Pyr_large = fitrm(Pyr_matched_dfof_pupil, 'Ld1c1-Ld2c3 ~ 1', 'WithinDesign', w);
ranova(rm_Pyr_large, 'withinmodel', 'DART*contrast')

%% Figure S4H - Normalized difference for SST cells
norm_diff_red = norm_diff_pupil(:,:,sst_inds);

facil_red=norm_diff_red(:,:,:)>=1;
supp_red=norm_diff_red(:,:,:)<=-1;

N=length(sst_inds);
facil_table_smallPupil = sum(facil_red(1,:,:),3)/N;
supp_table_smallPupil = sum(supp_red(1,:,:),3)/N;
facil_table_largePupil = sum(facil_red(2,:,:),3)/N;
supp_table_largePupil = sum(supp_red(2,:,:),3)/N;

figure;
subplot(1,2,1)
b=bar([1,2,3],[supp_table_smallPupil; supp_table_largePupil],'grouped','FaceColor',"#00ffff",'EdgeColor', [1 1 1]);
b(1).FaceColor="#70D0F6"
b(2).FaceColor="#0C8ABB"
xticklabels({'25','50','100'})
ylim([0 .4])
title('Suppressed')
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

subplot(1,2,2)
b=bar([1,2,3],[facil_table_smallPupil; facil_table_largePupil],'FaceColor',"#a329cc",'EdgeColor', [1 1 1]);
b(1).FaceColor="#C983B1"
b(2).FaceColor="#883367"
xticklabels({'25','50','100'})
ylim([0 .4])
title('Facilitated')
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
%% Figure S4H - Chi Squared for SST cells

%compute chi squares for suppression, small pupil vs large pupil
%25% contrast

n1=supp_table_smallPupil(1)*N;
n2=supp_table_largePupil(1)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil1,p1] = crosstab(x1,x2);

%50
n1=fix(supp_table_smallPupil(2)*N);
n2=fix(supp_table_largePupil(2)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil2,p2] = crosstab(x1,x2);

%100
n1=supp_table_smallPupil(3)*N;
n2=supp_table_largePupil(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil3,p3] = crosstab(x1,x2);



[chi2smallPupil1, chi2smallPupil2, chi2smallPupil3; p1*3, p2*3,p3*3]

clear h p1 p2 p3 chi2smallPupil1 chi2smallPupil2 chi2smallPupil3

%compute chi squares for facilitation
%25
n1=fix(facil_table_smallPupil(1)*N);
n2=fix(facil_table_largePupil(1)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil1,p1] = crosstab(x1,x2);

%50
n1=fix(facil_table_smallPupil(2)*N);
n2=fix(facil_table_largePupil(2)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil2,p2] = crosstab(x1,x2);

%100
n1=fix(facil_table_smallPupil(3)*N);
n2=fix(facil_table_largePupil(3)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil3,p3] = crosstab(x1,x2);

[chi2smallPupil1, chi2smallPupil2, chi2smallPupil3; p1, p2,p3]
clear h p1 p2 p3 chi2smallPupil1 chi2smallPupil2 chi2smallPupil3 n1 n2 x1 x2

%% Figure S4I - Normalized difference for Pyramidal Cells

norm_diff_green = norm_diff_pupil(:,:,pyr_inds);

facil_green=norm_diff_green(:,:,:)>=1;
supp_green=norm_diff_green(:,:,:)<=-1;

N=length(pyr_inds);
facil_table_smallPupil_green = sum(facil_green(1,:,:),3)/N;
supp_table_smallPupil_green = sum(supp_green(1,:,:),3)/N;
facil_table_largePupil_green = sum(facil_green(2,:,:),3)/N;
supp_table_largePupil_green = sum(supp_green(2,:,:),3)/N;

figure;
subplot(1,2,1)
b=bar([1,2,3],[supp_table_smallPupil_green; supp_table_largePupil_green],'grouped','FaceColor',"#00ffff",'EdgeColor', [1 1 1]);
b(1).FaceColor="#70D0F6"
b(2).FaceColor="#0C8ABB"
xticklabels({'25','50','100'})
ylim([0 .4])
title('Suppressed')
ylabel(["Fraction Pyr cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

subplot(1,2,2)
b=bar([1,2,3],[facil_table_smallPupil_green; facil_table_largePupil_green],'FaceColor',"#a329cc",'EdgeColor', [1 1 1]);
b(1).FaceColor="#C983B1"
b(2).FaceColor="#883367"
xticklabels({'25','50','100'})
ylim([0 .4])
title('Facilitated')
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
%% Figure S4I - Chi Squared for Pyramidal Cells
%25% contrast

n1=fix(supp_table_smallPupil_green(1)*N);
n2=fix(supp_table_largePupil_green(1)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil1,p1] = crosstab(x1,x2);

%50
n1=fix(supp_table_smallPupil_green(2)*N);
n2=fix(supp_table_largePupil_green(2)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil2,p2] = crosstab(x1,x2);

%100
n1=fix(supp_table_smallPupil_green(3)*N);
n2=fix(supp_table_largePupil_green(3)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil3,p3] = crosstab(x1,x2);


[chi2smallPupil1, chi2smallPupil2, chi2smallPupil3; p1*3, p2*3,p3*3]

clear h p1 p2 p3 chi2smallPupil1 chi2smallPupil2 chi2smallPupil3


%compute chi squares for facilitation
%25
n1=fix(facil_table_smallPupil_green(1)*N);
n2=fix(facil_table_largePupil_green(1)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil1,p1] = crosstab(x1,x2);

%50
n1=fix(facil_table_smallPupil_green(2)*N);
n2=fix(facil_table_largePupil_green(2)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil2,p2] = crosstab(x1,x2);

%100
n1=fix(facil_table_smallPupil_green(3)*N);
n2=fix(facil_table_largePupil_green(3)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil3,p3] = crosstab(x1,x2);

[chi2smallPupil1, chi2smallPupil2, chi2smallPupil3; p1, p2,p3]
clear h p1 p2 p3 chi2smallPupil1 chi2smallPupil2 chi2smallPupil3 n1 n2 x1 x2
