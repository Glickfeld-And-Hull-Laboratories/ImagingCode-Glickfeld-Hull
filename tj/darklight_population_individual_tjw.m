%%
%clear everything
clear all
clear all global
clc
close all

%%
%find folders to load and experiment info

fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight\multi_day_drift'; %folder to save files to
dataset = 'exp_list_darklight_actual_tjw'; %experiment list to pick files from
eval(dataset); %load dataset


%%

%we start with loading the proper information based on the mouse's info in
%the experiment list

d1 = 34; 
d2 = d1+1; 
d3 = d1+2;
d4 = d1+3;
mouse = expt(d1).mouse; 
ref_str_d1 = ['runs-',expt(d1).runs];
ref_str_d2 = ['runs-',expt(d2).runs]; 
ref_str_d3 = ['runs-',expt(d3).runs]; 
ref_str_d4 = ['runs-',expt(d4).runs];
img_area = expt(d1).img_loc{1};
img_layer = expt(d1).img_loc{2};
date_d1 = expt(d1).date; 
date_d2 = expt(d2).date; 
date_d3 = expt(d3).date; 
date_d4 = expt(d4).date; 
img_folder_d1 = expt(d1).runs; 
img_folder_d2 = expt(d2).runs; 
img_folder_d3 = expt(d3).runs; 
img_folder_d4 = expt(d4).runs; 

%load ori info for each day
d1_ori = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'oriTuningInfo.mat']));
d2_ori = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'oriTuningInfo.mat']));
d3_ori = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'oriTuningInfo.mat']));
d4_ori = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'oriTuningInfo.mat']));

d1_fits = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'oriTuningAndFits.mat']));
d2_fits = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'oriTuningAndFits.mat']));
d3_fits = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'oriTuningAndFits.mat']));
d4_fits = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'oriTuningAndFits.mat']));

%load multiday data for day 2
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));
d3_matches = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'multiday_alignment.mat']));
d4_matches = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'multiday_alignment.mat']));


%% matching indices

tuned_d1 = d1_ori.ind_theta90;
tuned_d2 = d2_ori.ind_theta90;
tuned_d3 = d3_ori.ind_theta90;
tuned_d4 = d4_ori.ind_theta90;

tuned_all = unique([tuned_d1 tuned_d2 tuned_d3 tuned_d4]);

match_d2 = find([d2_matches.cellImageAlign.pass]); 
match_d3 = find([d3_matches.cellImageAlign.pass]); 
match_d4 = find([d4_matches.cellImageAlign.pass]); 

match_all = intersect(intersect(intersect(match_d2,match_d3),match_d4),tuned_all);


%import avg response to each ori for each cell on each day
a = d1_fits.avgResponseEaOri(match_all,:);
b = d2_fits.avgResponseEaOri(match_all,:);
c = d3_fits.avgResponseEaOri(match_all,:);
d = d4_fits.avgResponseEaOri(match_all,:);

%%
%make correlation matrices of base1 to other sessions for avg ori
basecorr = corr(a,b);
figure;
subplot(2,2,1)
heatmap(basecorr)
xlabel('Orientation Base1')
ylabel('Orientation Base2')
colormap default
subplot(2,2,2)
darkcorr = corr(a,c);
heatmap(darkcorr)
xlabel('Orientation Base1')
ylabel('Orientation Post-dark')
colormap default
subplot(2,2,3)
lightcorr = corr(a,d);
heatmap(lightcorr)
xlabel('Orientation Base1')
ylabel('Orientation Post-dark+7d')
colormap default

%%
%make and plot diagonals (same ori across sessions)
diag_mat = [];
diag_mat(:,1) = diag(basecorr);
diag_mat(:,2) = diag(darkcorr);
diag_mat(:,3) = diag(lightcorr);

figure;
plot(diag_mat)
xlim([0,9])
ylim([0,1])
xlabel('Orientation')
ylabel('Population Vector Correlation')
legend(['Base1 - Base2'], ['Base1 - Post-dark'], ['Base1 - Post-dark+7d'], 'Location', 'Best')
newcolors = [0 0 0
             .5 .5 .5
             .8 .8 .8];
colororder(newcolors)


%%
%base1 correlated with itself
base1corr_base1corr = corr(a,a);
off_diag_base1_base1 = [];

for i = 1:(length(base1corr_base1corr)-1)
    off_diag_base1_base1(i) = base1corr_base1corr(i,i+1);
end
off_diag_base1_base1 = off_diag_base1_base1';

%base2 correlated w/ itself
base2corr_base2corr = corr(b,b);
off_diag_base2_base2 = [];

for i = 1:(length(base2corr_base2corr)-1)
    off_diag_base2_base2(i) = base2corr_base2corr(i,i+1);
end
off_diag_base2_base2 = off_diag_base2_base2';

%post-dark correlated w/ itself
darkcorr_darkcorr = corr(c,c);
off_diag_dark_dark = [];

for i = 1:(length(darkcorr_darkcorr)-1)
    off_diag_dark_dark(i) = darkcorr_darkcorr(i,i+1);
end
off_diag_dark_dark = off_diag_dark_dark';

%postdark +7d light correlated w/ itself
lightcorr_lightcorr = corr(d,d);
off_diag_light_light = [];

for i = 1:(length(lightcorr_lightcorr)-1)
    off_diag_light_light(i) = lightcorr_lightcorr(i,i+1);
end
off_diag_light_light = off_diag_light_light';

%%
off_diag_mat = [];
off_diag_mat(:,1) = off_diag_base1_base1;
off_diag_mat(:,2) = off_diag_base2_base2;
off_diag_mat(:,3) = off_diag_dark_dark;
off_diag_mat(:,4) = off_diag_light_light;

%%
figure;
errorbar(mean(off_diag_mat), (std(off_diag_mat)/sqrt(length(off_diag_mat))))
xlim([0 5])
ylim([0 1])

figure;
errorbar([2,4,6], mean(diag_mat), (std(diag_mat)/sqrt(length(diag_mat))), '-o')
hold on;
errorbar([1,3,5,7], mean(off_diag_mat), (std(off_diag_mat)/sqrt(length(off_diag_mat))), '-o')
xlim([0 8])
xticks([0:8])
xticklabels({'','Base1 w/self','Base1 to Base2','Base2 w/self','Base1 to Post-dark','Post-dark w/self','Base1 to Post-dark+7d','Post-dark+7d w/self',''})
legend(['Between-session (Same ori)'], ['Within-session (n-1 ori)'], 'Location', 'Best')
ylabel('Population Correlation')

%%

%make and plot diagonals (same ori across sessions)
diag_mat_prev_ses = [];
diag_mat_prev_ses(:,1) = diag(corr(a,b));
diag_mat_prev_ses(:,2) = diag(corr(b,c));
diag_mat_prev_ses(:,3) = diag(corr(c,d));

figure;
errorbar(mean(diag_mat), (std(diag_mat)/sqrt(length(diag_mat))))
xlim([0 4])
ylim([0 1])
hold on
errorbar(mean(diag_mat_prev_ses), (std(diag_mat_prev_ses)/sqrt(length(diag_mat_prev_ses))))

%%
%farther off-diagonals
%real one!!!


%base1_base1
off_diag_base1_base1_22point5diff = [];

for i = 1:length(base1corr_base1corr)
    if i == 8
        off_diag_base1_base1_22point5diff(i) = base1corr_base1corr(i-7,i);
    else
        off_diag_base1_base1_22point5diff(i) = base1corr_base1corr(i+1,i);
    end
       
end
off_diag_base1_base1_22point5diff = off_diag_base1_base1_22point5diff';


off_diag_base1_base1_45diff = [];

for i = 1:length(base1corr_base1corr)
    if i == 7
        off_diag_base1_base1_45diff(i) = base1corr_base1corr(i-6,i);
    elseif i == 8
        off_diag_base1_base1_45diff(i) = base1corr_base1corr(i-6,i);
    else
        off_diag_base1_base1_45diff(i) = base1corr_base1corr(i+2,i);
    end
       
end
off_diag_base1_base1_45diff = off_diag_base1_base1_45diff';


off_diag_base1_base1_67point5diff = [];

for i = 1:length(base1corr_base1corr)
    if i == 6
        off_diag_base1_base1_67point5diff(i) = base1corr_base1corr(i-5,i);
    elseif i == 7
        off_diag_base1_base1_67point5diff(i) = base1corr_base1corr(i-5,i);
    elseif i == 8
        off_diag_base1_base1_67point5diff(i) = base1corr_base1corr(i-5,i);
    else
        off_diag_base1_base1_67point5diff(i) = base1corr_base1corr(i+3,i);
    end
       
end
off_diag_base1_base1_67point5diff = off_diag_base1_base1_67point5diff';


off_diag_base1_base1_90diff = [];

for i = 1:length(base1corr_base1corr)
    if i == 5
        off_diag_base1_base1_90diff(i) = base1corr_base1corr(i-4,i);
    elseif i == 6
        off_diag_base1_base1_90diff(i) = base1corr_base1corr(i-4,i);
    elseif i == 7
        off_diag_base1_base1_90diff(i) = base1corr_base1corr(i-4,i);
    elseif i == 8
        off_diag_base1_base1_90diff(i) = base1corr_base1corr(i-4,i);
    else
        off_diag_base1_base1_90diff(i) = base1corr_base1corr(i+4,i);
    end
       
end
off_diag_base1_base1_90diff = off_diag_base1_base1_90diff';

%base2_base2
off_diag_base2_base2_22point5diff = [];

for i = 1:length(base2corr_base2corr)
    if i == 8
        off_diag_base2_base2_22point5diff(i) = base2corr_base2corr(i-7,i);
    else
        off_diag_base2_base2_22point5diff(i) = base2corr_base2corr(i+1,i);
    end
       
end
off_diag_base2_base2_22point5diff = off_diag_base2_base2_22point5diff';


off_diag_base2_base2_45diff = [];

for i = 1:length(base2corr_base2corr)
    if i == 7
        off_diag_base2_base2_45diff(i) = base2corr_base2corr(i-6,i);
    elseif i == 8
        off_diag_base2_base2_45diff(i) = base2corr_base2corr(i-6,i);
    else
        off_diag_base2_base2_45diff(i) = base2corr_base2corr(i+2,i);
    end
       
end
off_diag_base2_base2_45diff = off_diag_base2_base2_45diff';


off_diag_base2_base2_67point5diff = [];

for i = 1:length(base2corr_base2corr)
    if i == 6
        off_diag_base2_base2_67point5diff(i) = base2corr_base2corr(i-5,i);
    elseif i == 7
        off_diag_base2_base2_67point5diff(i) = base2corr_base2corr(i-5,i);
    elseif i == 8
        off_diag_base2_base2_67point5diff(i) = base2corr_base2corr(i-5,i);
    else
        off_diag_base2_base2_67point5diff(i) = base2corr_base2corr(i+3,i);
    end
       
end
off_diag_base2_base2_67point5diff = off_diag_base2_base2_67point5diff';


off_diag_base2_base2_90diff = [];

for i = 1:length(base2corr_base2corr)
    if i == 5
        off_diag_base2_base2_90diff(i) = base2corr_base2corr(i-4,i);
    elseif i == 6
        off_diag_base2_base2_90diff(i) = base2corr_base2corr(i-4,i);
    elseif i == 7
        off_diag_base2_base2_90diff(i) = base2corr_base2corr(i-4,i);
    elseif i == 8
        off_diag_base2_base2_90diff(i) = base2corr_base2corr(i-4,i);
    else
        off_diag_base2_base2_90diff(i) = base2corr_base2corr(i+4,i);
    end
       
end
off_diag_base2_base2_90diff = off_diag_base2_base2_90diff';

%dark_dark
off_diag_dark_dark_22point5diff = [];

for i = 1:length(darkcorr_darkcorr)
    if i == 8
        off_diag_dark_dark_22point5diff(i) = darkcorr_darkcorr(i-7,i);
    else
        off_diag_dark_dark_22point5diff(i) = darkcorr_darkcorr(i+1,i);
    end
       
end
off_diag_dark_dark_22point5diff = off_diag_dark_dark_22point5diff';


off_diag_dark_dark_45diff = [];

for i = 1:length(darkcorr_darkcorr)
    if i == 7
        off_diag_dark_dark_45diff(i) = darkcorr_darkcorr(i-6,i);
    elseif i == 8
        off_diag_dark_dark_45diff(i) = darkcorr_darkcorr(i-6,i);
    else
        off_diag_dark_dark_45diff(i) = darkcorr_darkcorr(i+2,i);
    end
       
end
off_diag_dark_dark_45diff = off_diag_dark_dark_45diff';


off_diag_dark_dark_67point5diff = [];

for i = 1:length(darkcorr_darkcorr)
    if i == 6
        off_diag_dark_dark_67point5diff(i) = darkcorr_darkcorr(i-5,i);
    elseif i == 7
        off_diag_dark_dark_67point5diff(i) = darkcorr_darkcorr(i-5,i);
    elseif i == 8
        off_diag_dark_dark_67point5diff(i) = darkcorr_darkcorr(i-5,i);
    else
        off_diag_dark_dark_67point5diff(i) = darkcorr_darkcorr(i+3,i);
    end
       
end
off_diag_dark_dark_67point5diff = off_diag_dark_dark_67point5diff';


off_diag_dark_dark_90diff = [];

for i = 1:length(darkcorr_darkcorr)
    if i == 5
        off_diag_dark_dark_90diff(i) = darkcorr_darkcorr(i-4,i);
    elseif i == 6
        off_diag_dark_dark_90diff(i) = darkcorr_darkcorr(i-4,i);
    elseif i == 7
        off_diag_dark_dark_90diff(i) = darkcorr_darkcorr(i-4,i);
    elseif i == 8
        off_diag_dark_dark_90diff(i) = darkcorr_darkcorr(i-4,i);
    else
        off_diag_dark_dark_90diff(i) = darkcorr_darkcorr(i+4,i);
    end
       
end
off_diag_dark_dark_90diff = off_diag_dark_dark_90diff';

%light_light
off_diag_light_light_22point5diff = [];

for i = 1:length(lightcorr_lightcorr)
    if i == 8
        off_diag_light_light_22point5diff(i) = lightcorr_lightcorr(i-7,i);
    else
        off_diag_light_light_22point5diff(i) = lightcorr_lightcorr(i+1,i);
    end
       
end
off_diag_light_light_22point5diff = off_diag_light_light_22point5diff';


off_diag_light_light_45diff = [];

for i = 1:length(lightcorr_lightcorr)
    if i == 7
        off_diag_light_light_45diff(i) = lightcorr_lightcorr(i-6,i);
    elseif i == 8
        off_diag_light_light_45diff(i) = lightcorr_lightcorr(i-6,i);
    else
        off_diag_light_light_45diff(i) = lightcorr_lightcorr(i+2,i);
    end
       
end
off_diag_light_light_45diff = off_diag_light_light_45diff';


off_diag_light_light_67point5diff = [];

for i = 1:length(lightcorr_lightcorr)
    if i == 6
        off_diag_light_light_67point5diff(i) = lightcorr_lightcorr(i-5,i);
    elseif i == 7
        off_diag_light_light_67point5diff(i) = lightcorr_lightcorr(i-5,i);
    elseif i == 8
        off_diag_light_light_67point5diff(i) = lightcorr_lightcorr(i-5,i);
    else
        off_diag_light_light_67point5diff(i) = lightcorr_lightcorr(i+3,i);
    end
       
end
off_diag_light_light_67point5diff = off_diag_light_light_67point5diff';


off_diag_light_light_90diff = [];

for i = 1:length(lightcorr_lightcorr)
    if i == 5
        off_diag_light_light_90diff(i) = lightcorr_lightcorr(i-4,i);
    elseif i == 6
        off_diag_light_light_90diff(i) = lightcorr_lightcorr(i-4,i);
    elseif i == 7
        off_diag_light_light_90diff(i) = lightcorr_lightcorr(i-4,i);
    elseif i == 8
        off_diag_light_light_90diff(i) = lightcorr_lightcorr(i-4,i);
    else
        off_diag_light_light_90diff(i) = lightcorr_lightcorr(i+4,i);
    end
       
end
off_diag_light_light_90diff = off_diag_light_light_90diff';


%%

off_diag_22point5diff(:,1) = off_diag_base1_base1_22point5diff;
off_diag_22point5diff(:,2) = off_diag_base2_base2_22point5diff;
off_diag_22point5diff(:,3) = off_diag_dark_dark_22point5diff;
off_diag_22point5diff(:,4) = off_diag_light_light_22point5diff;

off_diag_45diff(:,1) = off_diag_base1_base1_45diff;
off_diag_45diff(:,2) = off_diag_base2_base2_45diff;
off_diag_45diff(:,3) = off_diag_dark_dark_45diff;
off_diag_45diff(:,4) = off_diag_light_light_45diff;

off_diag_67point5diff(:,1) = off_diag_base1_base1_67point5diff;
off_diag_67point5diff(:,2) = off_diag_base2_base2_67point5diff;
off_diag_67point5diff(:,3) = off_diag_dark_dark_67point5diff;
off_diag_67point5diff(:,4) = off_diag_light_light_67point5diff;

off_diag_90diff(:,1) = off_diag_base1_base1_90diff;
off_diag_90diff(:,2) = off_diag_base2_base2_90diff;
off_diag_90diff(:,3) = off_diag_dark_dark_90diff;
off_diag_90diff(:,4) = off_diag_light_light_90diff;

%%
figure;
errorbar([2,4,6,8], mean(off_diag_22point5diff), (std(off_diag_22point5diff)/sqrt(length(off_diag_22point5diff))), '-o', 'Color', [0 0 0], "MarkerSize", 9)
hold on
errorbar([2,4,6,8], mean(off_diag_45diff), (std(off_diag_45diff)/sqrt(length(off_diag_45diff))), '-diamond', 'Color', [.3 .3 .3], "MarkerSize", 9)
errorbar([2,4,6,8], mean(off_diag_67point5diff), (std(off_diag_67point5diff)/sqrt(length(off_diag_67point5diff))), '-^', 'Color', [.6 .6 .6], "MarkerSize", 9)
errorbar([2,4,6,8], mean(off_diag_90diff), (std(off_diag_90diff)/sqrt(length(off_diag_90diff))), '-square', 'Color', [.9 .9 .9], "MarkerSize", 10)
xlim([0 10])
xticks([0:10])
xticklabels({'','','Base1 w/self','','Base2 w/self','','Post-dark w/self','','Post-dark+7d w/self'})
legend(['22.5 degrees away'], ['45 degrees away'], ['67.5 degrees away'],  ['90 degrees away'], 'FontSize', 12, 'Location', 'Best')
ylabel('Population Correlation')


%%
figure;
errorbar([2,4,6], mean(diag_mat), (std(diag_mat)/sqrt(length(diag_mat))), '-o', 'Color', 'red', "MarkerSize", 9)
hold on;
errorbar([1,3,5,7], mean(off_diag_22point5diff), (std(off_diag_22point5diff)/sqrt(length(off_diag_22point5diff))), '-o', 'Color', [0 0 0], "MarkerSize", 9)
errorbar([1,3,5,7], mean(off_diag_45diff), (std(off_diag_45diff)/sqrt(length(off_diag_45diff))), '-diamond', 'Color', [.3 .3 .3], "MarkerSize", 9)
errorbar([1,3,5,7], mean(off_diag_67point5diff), (std(off_diag_67point5diff)/sqrt(length(off_diag_67point5diff))), '-^', 'Color', [.6 .6 .6], "MarkerSize", 9)
errorbar([1,3,5,7], mean(off_diag_90diff), (std(off_diag_90diff)/sqrt(length(off_diag_90diff))), '-square', 'Color', [.9 .9 .9], "MarkerSize", 10)
xlim([0 8])
xticks([0:8])
xticklabels({'','Base1 w/self','Base1 to Base2','Base2 w/self','Base1 to Post-dark','Post-dark w/self','Base1 to Post-dark+7d','Post-dark+7d w/self',''})
legend(['Same ori (b/w session)'],['22.5 degrees away (w/i session)'], ['45 degrees away (w/i session)'], ['67.5 degrees away (w/i session)'],  ['90 degrees away (w/i session)'], 'FontSize', 12, 'Location', 'Best')
ylabel('Population Correlation')

%%
%save vars for group analyses

save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_diag_matrices.mat']), 'diag_mat', 'off_diag_mat', 'diag_mat_prev_ses', ...
    'off_diag_22point5diff', 'off_diag_45diff', 'off_diag_67point5diff', 'off_diag_90diff')