%This code is used to analyze individual subject data from the ketamine 10mg/kg 2p imaging experiment, so that these data are prepared for
%pooled data afterwards.
%%
%clear everything
clear all
clear all global
clc
close all

%%
%find folders to load and experiment info

d = [8:10 12];
newD1Match = 0;

dataset = 'exp_list_ket_1wk_tjw'; %experiment list to pick files from
eval(dataset); %load dataset
data_owner = expt(d(1)).folder;
drug_list = {expt(d).drug};

if strcmp(data_owner,'lindsey')
    fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P'; %folder to load files from
elseif strcmp(data_owner,'tj')
    fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
end
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\ket\'; %folder to save files to
% newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\reverse_match\ket_4sess\'; %folder to save files to

if ~exist(newfnout)
    mkdir(newfnout)
end


mouse = expt(d(1)).mouse; 
img_area = expt(d(1)).img_loc{1};
img_layer = expt(d(1)).img_loc{2};
for i = 1:length(d);
    ref_str{i} = ['runs-',expt(d(i)).runs];
    date{i} = expt(d(i)).date; 
    img_folder{i} = expt(d(i)).runs; 
    %load ori info for each day
    if newD1Match & i == 1
        fnout_use = fullfile(fnout, [date{i} '_' mouse], [date{i} '_' mouse '_' ref_str{i}], 'newD1');
    elseif  newD1Match & i > 1
        fnout_use = fullfile(fnout, [date{i} '_' mouse], [date{i} '_' mouse '_' ref_str{i}], 'newD1Match');
    else 
        fnout_use = fullfile(fnout, [date{i} '_' mouse], [date{i} '_' mouse '_' ref_str{i}]);
    end

    alld_ori{i} = load(fullfile(fnout_use, [date{i} '_' mouse '_' ref_str{i} '_' 'oriTuningInfo.mat']));
    alld_fits{i} = load(fullfile(fnout_use, [date{i} '_' mouse '_' ref_str{i} '_' 'oriTuningAndFits.mat']));
    alld_new_pref{i} = load(fullfile(fnout_use, [date{i} '_' mouse '_' ref_str{i} '_' 'new_pref.mat']));
    %k and max vals
    alld_k_and_max{i} = load(fullfile(fnout_use, [date{i} '_' mouse '_' ref_str{i} '_' 'k_and_max_vals.mat']));

    if i > 1
        alld_matches{i} = load(fullfile(fnout_use, [date{i} '_' mouse '_' ref_str{i} '_' 'multiday_alignment.mat']));
    end
    alld_tuned{i} = alld_ori{i}.ind_theta90;

    legend_str{i} = ['Session ' num2str(i)];
end

% matching indices
for i = 1:length(d)-1
    for ii = 2:length(d)
        if ii > i
            if i == 1
                d_match{ii} =  find([alld_matches{ii}.cellImageAlign.pass]);
            end
            tuned_comp{i,ii} = intersect(alld_tuned{i}, alld_tuned{ii});
            match_comp{i,ii} =  intersect(tuned_comp{ii},d_match{ii});
        end
    end
end

temp_tuned = [];
for ii = 1:length(d)
    temp_tuned = [temp_tuned alld_tuned{i}];
end
%tuned on at least one day
tuned_all = unique(temp_tuned);

if length(d)>2
    temp_match = [];
    for ii = 2:length(d)-1
        temp_match = intersect(d_match{ii+1},d_match{ii});
    end
else
    temp_match = d_match{2};
end
match_all = intersect(temp_match,tuned_all);
for i = 1:length(d)
    alld_k{i} = alld_k_and_max{i}.k1(match_all);
    alld_max{i} = alld_k_and_max{i}.max_dfof(match_all);
end


%%
%Tuning width
fig = figure;
sgtitle('Tuning Width (k) by Session-')
subplot(2,1,1)
for i = 1:length(d)
    cdfplot(alld_k{i});
    hold on    
end
legend(legend_str, 'Location', 'southeast')
xlabel(['Tuning Width (k)'])
ylabel(['% of cells'])
title([])
hold off
subplot(2,1,2)
for i = 2
    temp = (alld_k{1}-alld_k{i})./(alld_k{1}+alld_k{i});
    histogram(temp,[-1:.1:1])
    hold on
end
xlim([-1 1])
xlabel('D1-DN/D1+DN')
ylabel('Number of cells')
print(fullfile(newfnout, [mouse,  '_k_by_sess.pdf']), '-dpdf', '-bestfit')

%Max dF/F
fig = figure;
subplot(2,1,1)
sgtitle('Max df/f by Session')
for i = 1:length(d)
    cdfplot(alld_max{i});
    hold on    
end
legend(legend_str, 'Location', 'southeast')
xlabel(['Max df/f'])
ylabel(['% of cells'])
title([])
hold off
subplot(2,1,2)
for i = 2
    temp = (alld_max{1}-alld_max{i})./(alld_max{1}+alld_max{i});
    histogram(temp,[-1:.1:1])
    hold on
end
xlim([-1 1])
xlabel('D1-DN/D1+DN')
ylabel('Number of cells')
print(fullfile(newfnout, [mouse,  '_max_by_sess.pdf']), '-dpdf', '-bestfit')

%within session pref ori change
for i = 1:length(d)
    alld_within_prefori{i,1} = changeprefto90TJ(alld_new_pref{i}.prefOri_1);
    alld_within_prefori{i,2} = changeprefto90TJ(alld_new_pref{i}.prefOri_2);
    alld_score_prefori_within{i} = abs(alld_within_prefori{i,1}(alld_tuned{i})-alld_within_prefori{i,2}(alld_tuned{i}));
end

fig = figure;
sgtitle('Pref Ori Changes within Session')
for i = 1:length(d)
    cdfplot(alld_score_prefori_within{i});
    hold on    
end
legend(legend_str, 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori w/i Session'])
ylabel(['% of cells'])
title([])
hold off
print(fullfile(newfnout, [mouse,  '_prefori_change_within_sess.pdf']), '-dpdf', '-bestfit')

%pref ori changes
for i = 1:length(d)
    temp_prefori = alld_ori{i}.prefOri(1,:);
    alld_dscore_prefori{i} = double(temp_prefori>90);
    alld_dscore_prefori{i}(alld_dscore_prefori{i}>0) = 180;
    alld_dscore_prefori{i} = abs(alld_dscore_prefori{i}-temp_prefori);
end
for i = 1:length(d)-1
    d1_comp_prefori{i} = abs(alld_dscore_prefori{i}(match_comp{1,i+1})-alld_dscore_prefori{i+1}(match_comp{1,i+1}));
    d1_comp_legend{i} = ['D1 vs D' num2str(i+1)];
end

for i = 1:length(d)-1
    for ii = 2:length(d)
        if ii>i
            alld_comp_prefori{i,ii} = abs(alld_dscore_prefori{i}(match_comp{i,ii})-alld_dscore_prefori{ii}(match_comp{i,ii}));
        end
    end
end

nDiff= length(d)-1;
sessDiff_comp_prefori = cell(1,nDiff);
for i = 1:length(d)-1
    for ii = 2:length(d)
        if ii>i
            sessDiff_comp_prefori{ii-i} = [sessDiff_comp_prefori{ii-i} alld_comp_prefori{i,ii}];
        end
    end
end
for i = 1:nDiff
    sessDiff_comp_legend{i} = [num2str(i+1) ' session diff'];
end

fig = figure;
sgtitle('Pref Ori Changes from Day 1')
for i = 1:length(d)-1
    cdfplot(d1_comp_prefori{i});
    hold on    
end
xlim([0 90])
legend(d1_comp_legend, 'Location', 'southeast')
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title([])
hold off
print(fullfile(newfnout, [mouse,  '_prefori_change_fromd1.pdf']), '-dpdf', '-bestfit')

fig = figure;
sgtitle('Pref Ori Changes All Sessions')
for i = 1:nDiff
    cdfplot(sessDiff_comp_prefori{i});
    hold on    
end
xlim([0 90])
legend(sessDiff_comp_legend, 'Location', 'southeast')
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title([])
hold off
print(fullfile(newfnout, [mouse,  '_prefori_change_all.pdf']), '-dpdf', '-bestfit')
save(fullfile(newfnout, [mouse '_matchedData.mat']),'alld_k','alld_max','d1_comp_prefori','alld_score_prefori_within','drug_list')

%% fraction matched/unmatched
for i = 2:length(d)
    alld_nCells_tot{i} = size(alld_matches{i}.cellImageAlign,2);
    alld_ravg{i} = nan(1,alld_nCells_tot{i});
    alld_good_cells{i} = zeros(size(alld_ravg{i}));
    alld_matched_cells{i} = zeros(size(alld_ravg{i}));
    for iCell = 1:alld_nCells_tot{i}
        if alld_matches{i}.cellImageAlign(iCell).pass == 0
            if ~isempty(alld_matches{i}.cellImageAlign(iCell).d)
                alld_ravg{i}(iCell) = corr(alld_matches{i}.cellImageAlign(iCell).d(1).avg_img(:),alld_matches{i}.cellImageAlign(iCell).d(2).avg_img(:));
                alld_good_cells{i}(iCell) = 1;
            end
        else
            alld_good_cells{i}(iCell) = 1;
            alld_matched_cells{i}(iCell) = 1;
        end
    end
    alld_nCells_good{i} = sum(alld_good_cells{i});
    alld_nCells_matched{i} = sum(alld_matched_cells{i});
    alld_nCells_corr{i} = sum(alld_ravg{i}>0.85);

    alld_fract_matched{i} = alld_nCells_matched{i}./alld_nCells_good{i};
    alld_fract_missing{i} = alld_nCells_corr{i}./alld_nCells_good{i};
end

figure;
subplot(1,2,1)
for i = 2:length(d)
bar(i,alld_fract_matched{i})
hold on
end
xlabel('Session')
ylabel('Fraction cells')
ylim([0 1])
title('Matched')
%set(gca,'Xtick',2:length(d),'Xticklabel',num2str(2:length(d))')
subplot(1,2,2)
for i = 2:length(d)
bar(i,alld_fract_missing{i})
hold on
end
xlabel('Session')
ylabel('Fraction cells')
ylim([0 1])
title('Missing')
print(fullfile(newfnout, [mouse,  '_matched&missing.pdf']), '-dpdf', '-bestfit')



%%
%% NEW matching indices

new_tuned_d1 = d1_new_pref.new_tune;
new_tuned_d2 = d2_new_pref.new_tune;
new_tuned_d3 = d3_new_pref.new_tune;
new_tuned_d4 = d4_new_pref.new_tune;

new_tuned_d1_d2 = intersect(new_tuned_d1, new_tuned_d2);
new_tuned_d1_d3 = intersect(new_tuned_d1, new_tuned_d3);
new_tuned_d1_d4 = intersect(new_tuned_d1, new_tuned_d4);
new_tuned_d2_d3 = intersect(new_tuned_d2, new_tuned_d3);
new_tuned_d2_d4 = intersect(new_tuned_d2, new_tuned_d4);
new_tuned_d3_d4 = intersect(new_tuned_d3, new_tuned_d4);

new_tuned_all = unique([new_tuned_d1 new_tuned_d2 new_tuned_d3 new_tuned_d4]);

new_match_d2 = find([d2_matches.cellImageAlign.pass]); 
new_match_d3 = find([d3_matches.cellImageAlign.pass]); 
new_match_d4 = find([d4_matches.cellImageAlign.pass]); 

new_match_d1_d2 = intersect(new_tuned_d1_d2, new_match_d2);
new_match_d1_d3 = intersect(new_tuned_d1_d3, new_match_d3);
new_match_d1_d4 = intersect(new_tuned_d1_d4, new_match_d4);
new_match_d2_d3 = intersect(intersect(new_match_d2, new_match_d3), new_tuned_d2_d3);
new_match_d2_d4 = intersect(intersect(new_match_d2, new_match_d4), new_tuned_d2_d4);
new_match_d3_d4 = intersect(intersect(new_match_d3, new_match_d4), new_tuned_d3_d4);

new_match_all = intersect(intersect(intersect(new_match_d2,new_match_d3),new_match_d4),new_tuned_all);

%%
%within session theta change
d1_within_new_prefori_1 = changeprefto90TJ(d1_new_pref.theta_1);
d1_within_new_prefori_2 = changeprefto90TJ(d1_new_pref.theta_2);
d2_within_new_prefori_1 = changeprefto90TJ(d2_new_pref.theta_1);
d2_within_new_prefori_2 = changeprefto90TJ(d2_new_pref.theta_2);
d3_within_new_prefori_1 = changeprefto90TJ(d3_new_pref.theta_1);
d3_within_new_prefori_2 = changeprefto90TJ(d3_new_pref.theta_2);
d4_within_new_prefori_1 = changeprefto90TJ(d4_new_pref.theta_1);
d4_within_new_prefori_2 = changeprefto90TJ(d4_new_pref.theta_2);

d_score_new_prefori_within_d1 = abs(d1_within_new_prefori_1-d1_within_new_prefori_2);
d_score_new_prefori_within_d1 = d_score_new_prefori_within_d1(alld_tuned);
d_score_new_prefori_within_d2 = abs(d2_within_new_prefori_1-d2_within_new_prefori_2);
d_score_new_prefori_within_d2 = d_score_new_prefori_within_d2(tuned_d2);
d_score_new_prefori_within_d3 = abs(d3_within_new_prefori_1-d3_within_new_prefori_2);
d_score_new_prefori_within_d3 = d_score_new_prefori_within_d3(tuned_d3);
d_score_new_prefori_within_d4 = abs(d4_within_new_prefori_1-d4_within_new_prefori_2);
d_score_new_prefori_within_d4 = d_score_new_prefori_within_d4(tuned_d4);

fig = figure;
sgtitle('Theta Changes within Session')
% a=subplot(2,2,1)
h=cdfplot(d_score_new_prefori_within_d1);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d_score_new_prefori_within_d2);
set(j, 'LineStyle', '-', 'Color', [.2 .2 .2], 'LineWidth',1);
l=cdfplot(d_score_new_prefori_within_d3);
set(l, 'LineStyle', '-', 'Color', [.4 .4 .4], 'LineWidth',1);
m=cdfplot(d_score_new_prefori_within_d4);
set(m, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',1);
legend(['Session 1'], ['Session 2'], ['Session 3'], ['Session 4'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori w/i Session'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(newfnout, [mouse,  '_new_prefori_change_within_sess.pdf']), '-dpdf', '-bestfit')



%%
d1_new_prefori = d1_new_pref.theta;
dscore_new_prefori_1 = double(d1_new_prefori>90);
dscore_new_prefori_1(dscore_new_prefori_1>0) = 180;
dscore_new_prefori_1 = abs(dscore_new_prefori_1-d1_new_prefori);

d2_new_prefori = d2_new_pref.theta;
dscore_new_prefori_2 = double(d2_new_prefori>90);
dscore_new_prefori_2(dscore_new_prefori_2>0) = 180;
dscore_new_prefori_2 = abs(dscore_new_prefori_2-d2_new_prefori);

d3_new_prefori = d3_new_pref.theta;
dscore_new_prefori_3 = double(d3_new_prefori>90);
dscore_new_prefori_3(dscore_new_prefori_3>0) = 180;
dscore_new_prefori_3 = abs(dscore_new_prefori_3-d3_new_prefori);

d4_new_prefori = d4_new_pref.theta;
dscore_new_prefori_4 = double(d4_new_prefori>90);
dscore_new_prefori_4(dscore_new_prefori_4>0) = 180;
dscore_new_prefori_4 = abs(dscore_new_prefori_4-d4_new_prefori);


d_score_new_prefori_d1_d2 = abs(dscore_new_prefori_1-dscore_new_prefori_2);
d_score_new_prefori_d1_d2 = d_score_new_prefori_d1_d2(new_match_d1_d2);

d_score_new_prefori_d2_d3 = abs(dscore_new_prefori_2-dscore_new_prefori_3);
d_score_new_prefori_d2_d3 = d_score_new_prefori_d2_d3(new_match_d2_d3);

d_score_new_prefori_d3_d4 = abs(dscore_new_prefori_3-dscore_new_prefori_4);
d_score_new_prefori_d3_d4 = d_score_new_prefori_d3_d4(new_match_d3_d4);

d_score_new_prefori_d1_d3 = abs(dscore_new_prefori_1-dscore_new_prefori_3);
d_score_new_prefori_d1_d3 = d_score_new_prefori_d1_d3(new_match_d1_d3);

d_score_new_prefori_d2_d4 = abs(dscore_new_prefori_2-dscore_new_prefori_4);
d_score_new_prefori_d2_d4 = d_score_new_prefori_d2_d4(new_match_d2_d4);

d_score_new_prefori_d1_d4 = abs(dscore_new_prefori_1-dscore_new_prefori_4);
d_score_new_prefori_d1_d4 = d_score_new_prefori_d1_d4(new_match_d1_d4);

%%

fig = figure;
sgtitle('Vec Sum Pref Ori Changes from Day 1')
% a=subplot(2,2,1)
h=cdfplot(d_score_new_prefori_d1_d2);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d_score_new_prefori_d1_d3);
set(j, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',1);
l=cdfplot(d_score_new_prefori_d1_d4);
set(l, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',1);
legend(['1 sess'], ['2 sess'], ['3 sess'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(newfnout, [mouse,  '_vecsum_prefori_change_fromd1.pdf']), '-dpdf', '-bestfit')

%%

d_score_new_prefori_1sess = [d_score_new_prefori_d1_d2 d_score_new_prefori_d2_d3 d_score_new_prefori_d3_d4];
d_score_new_prefori_2sess = [d_score_new_prefori_d1_d3 d_score_new_prefori_d2_d4];
d_score_new_prefori_3sess = [d_score_new_prefori_d1_d4];

%%

fig = figure;
sgtitle('Vec Sum Pref Ori Changes All Sessions')
% a=subplot(2,2,1)
h=cdfplot(d_score_new_prefori_1sess);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',1);
hold on
j=cdfplot(d_score_new_prefori_2sess);
set(j, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',1);
l=cdfplot(d_score_new_prefori_3sess);
set(l, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',1);
legend(['1 sess'], ['2 sess'], ['3 sess'], 'Location', 'southeast')
xlim([0 90])
xlabel(['Change in Pref Ori'])
ylabel(['% of cells'])
title([])
hold off

print(fullfile(newfnout, [mouse,  '_vecsum_prefori_change_all.pdf']), '-dpdf', '-bestfit')


%%

wk1_overlap = length(intersect(alld_tuned, new_tuned_d1)) / ((length(intersect(alld_tuned, new_tuned_d1))) + (length(setdiff(new_tuned_d1, alld_tuned))) + (length(setdiff(alld_tuned, new_tuned_d1))));
wk1_vecsum =  length(setdiff(new_tuned_d1, alld_tuned)) / ((length(intersect(alld_tuned, new_tuned_d1))) + (length(setdiff(new_tuned_d1, alld_tuned))) + (length(setdiff(alld_tuned, new_tuned_d1))));
wk1_tradition =  length(setdiff(alld_tuned, new_tuned_d1)) / ((length(intersect(alld_tuned, new_tuned_d1))) + (length(setdiff(new_tuned_d1, alld_tuned))) + (length(setdiff(alld_tuned, new_tuned_d1))));

wk2_overlap = length(intersect(tuned_d2, new_tuned_d2)) / ((length(intersect(tuned_d2, new_tuned_d2))) + (length(setdiff(new_tuned_d2, tuned_d2))) + (length(setdiff(tuned_d2, new_tuned_d2))));
wk2_vecsum =  length(setdiff(new_tuned_d2, tuned_d2)) / ((length(intersect(tuned_d2, new_tuned_d2))) + (length(setdiff(new_tuned_d2, tuned_d2))) + (length(setdiff(tuned_d2, new_tuned_d2))));
wk2_tradition =  length(setdiff(tuned_d2, new_tuned_d2)) / ((length(intersect(tuned_d2, new_tuned_d2))) + (length(setdiff(new_tuned_d2, tuned_d2))) + (length(setdiff(tuned_d2, new_tuned_d2))));

wk3_overlap = length(intersect(tuned_d3, new_tuned_d3)) / ((length(intersect(tuned_d3, new_tuned_d3))) + (length(setdiff(new_tuned_d3, tuned_d3))) + (length(setdiff(tuned_d3, new_tuned_d3))));
wk3_vecsum =  length(setdiff(new_tuned_d3, tuned_d3)) / ((length(intersect(tuned_d3, new_tuned_d3))) + (length(setdiff(new_tuned_d3, tuned_d3))) + (length(setdiff(tuned_d3, new_tuned_d3))));
wk3_tradition =  length(setdiff(tuned_d3, new_tuned_d3)) / ((length(intersect(tuned_d3, new_tuned_d3))) + (length(setdiff(new_tuned_d3, tuned_d3))) + (length(setdiff(tuned_d3, new_tuned_d3))));

wk4_overlap = length(intersect(tuned_d4, new_tuned_d4)) / ((length(intersect(tuned_d4, new_tuned_d4))) + (length(setdiff(new_tuned_d4, tuned_d4))) + (length(setdiff(tuned_d4, new_tuned_d4))));
wk4_vecsum =  length(setdiff(new_tuned_d4, tuned_d4)) / ((length(intersect(tuned_d4, new_tuned_d4))) + (length(setdiff(new_tuned_d4, tuned_d4))) + (length(setdiff(tuned_d4, new_tuned_d4))));
wk4_tradition =  length(setdiff(tuned_d4, new_tuned_d4)) / ((length(intersect(tuned_d4, new_tuned_d4))) + (length(setdiff(new_tuned_d4, tuned_d4))) + (length(setdiff(tuned_d4, new_tuned_d4))));

%%

% figure;
% x_lab = categorical(["Week 1", "Week 2", "Week 3", "Week 4", "Week 5"]);
% bar(x_lab, [wk1_overlap*100, wk2_overlap*100, wk3_overlap*100, wk4_overlap*100, wk5_overlap*100]);
% ylabel("% of Cells");
% title(["Overlap of Tuning for Traditional and Vector Sum"])

figure;
x_lab = categorical(["Sess 1", "Sess 2", "Sess 3", "Sess 4"]);
y = [wk1_overlap*100, wk1_vecsum*100, wk1_tradition*100; wk2_overlap*100, wk2_vecsum*100, wk2_tradition*100; wk3_overlap*100, wk3_vecsum*100, ...
    wk3_tradition*100; wk4_overlap*100, wk4_vecsum*100, wk4_tradition*100;];
b = bar(x_lab, y, 'stacked','FaceColor','flat');
ylabel("% of Cells");
ylim([0 100])
legend(['Overlap'], ['Vector Sum Only'], ['Traditional Only'], 'Location', 'southeast')
title(["Overlap of Tuning for Traditional and Vector Sum"])

print(fullfile(newfnout, [mouse,  '_tuning_overlap.pdf']), '-dpdf', '-bestfit')

% for k = 1:size(y,2)
%     b(k).CData = k;
% end

%%

% figure;
% % subplot(3,2,1)
% h=scatter(dscore_prefori_1(intersect(tuned_d1,new_tuned_d1)), dscore_new_prefori_1(intersect(tuned_d1,new_tuned_d1)));
% set(h, "MarkerEdgeColor", [0 0 0]);
% hold on
% j=scatter(dscore_prefori_2(intersect(tuned_d2,new_tuned_d2)), dscore_new_prefori_2(intersect(tuned_d2,new_tuned_d2)));
% set(j, "MarkerEdgeColor", [.2 .2 .2]);
% l=scatter(dscore_prefori_3(intersect(tuned_d3,new_tuned_d3)), dscore_new_prefori_3(intersect(tuned_d3,new_tuned_d3)));
% set(l, "MarkerEdgeColor", [.4 .4 .4]);
% m=scatter(dscore_prefori_4(intersect(tuned_d4,new_tuned_d4)), dscore_new_prefori_4(intersect(tuned_d4,new_tuned_d4)));
% set(m, "MarkerEdgeColor", [.6 .6 .6]);
% o=scatter(dscore_prefori_5(intersect(tuned_d5,new_tuned_d5)), dscore_new_prefori_5(intersect(tuned_d5,new_tuned_d5)));
% set(o, "MarkerEdgeColor", [.8 .8 .8]);
% xlim([0 90])
% xlabel(['Change in Pref Ori'])
% ylabel(['% of cells'])
% title([])

fig = figure;
subplot(2,2,1)
h=scatter(dscore_prefori_1(intersect(alld_tuned,new_tuned_d1)), dscore_new_prefori_1(intersect(alld_tuned,new_tuned_d1)));
title('Sess 1')
hold on
subplot(2,2,2)
j=scatter(dscore_prefori_2(intersect(tuned_d2,new_tuned_d2)), dscore_new_prefori_2(intersect(tuned_d2,new_tuned_d2)));
title('Sess 2')
subplot(2,2,3)
l=scatter(dscore_prefori_3(intersect(tuned_d3,new_tuned_d3)), dscore_new_prefori_3(intersect(tuned_d3,new_tuned_d3)));
title('Sess 3')
subplot(2,2,4)
m=scatter(dscore_prefori_4(intersect(tuned_d4,new_tuned_d4)), dscore_new_prefori_4(intersect(tuned_d4,new_tuned_d4)));
title('Sess 4')
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Vec Sum Pref Ori');
xlabel(han,'Traditional Pref Ori');
title(han,'');

print(fullfile(newfnout, [mouse,  '_prefori_overlap.pdf']), '-dpdf', '-bestfit')


%%

%%

a = alld_fits.avgResponseEaOri(match_all,:);
b = d2_fits.avgResponseEaOri(match_all,:);
c = d3_fits.avgResponseEaOri(match_all,:);
d = d4_fits.avgResponseEaOri(match_all,:);


%%
a_cells_corr = corr(a.');
b_cells_corr = corr(b.');
c_cells_corr = corr(c.');
d_cells_corr = corr(d.');

%%
a_idx = eye(size(a_cells_corr));
a_Y = (1-a_idx).*a_cells_corr;
a_Z = a_cells_corr(~a_idx);

b_idx = eye(size(b_cells_corr));
b_Y = (1-b_idx).*b_cells_corr;
b_Z = b_cells_corr(~b_idx);

c_idx = eye(size(c_cells_corr));
c_Y = (1-c_idx).*c_cells_corr;
c_Z = c_cells_corr(~c_idx);

d_idx = eye(size(d_cells_corr));
d_Y = (1-d_idx).*d_cells_corr;
d_Z = d_cells_corr(~d_idx);


%%
abcd_mat = [a_Z b_Z c_Z d_Z];
abcd_corrs = corr(abcd_mat);

%%
figure;
scatter(1,[abcd_corrs(1,2), abcd_corrs(2,3), abcd_corrs(3,4)])
hold on
scatter(2, [abcd_corrs(1,3), abcd_corrs(2,4)])
scatter(3, [abcd_corrs(1,4)])
xlim([0 5])
ylim([0 1])
xlabel('Days Apart')
ylabel('Correlation')
xticklabels({'','','2','','4','','6',''})

print(fullfile(newfnout, [mouse,  '_tunecurvecorr_bysess.pdf']), '-dpdf', '-bestfit')

%%
figure;
errorbar(1,mean([abcd_corrs(1,2), abcd_corrs(2,3), abcd_corrs(3,4)]), std([abcd_corrs(1,2), abcd_corrs(2,3), abcd_corrs(3,4)]), 'o')
hold on
errorbar(2, mean([abcd_corrs(1,3), abcd_corrs(2,4)]), std([abcd_corrs(1,3), abcd_corrs(2,4)]), 'o')
scatter(3, [abcd_corrs(1,4)])
xlim([0 5])
ylim([0 1])
xlabel('Days Apart')
ylabel('Correlation')
xticklabels({'','','2','','4','','6',''})

print(fullfile(newfnout, [mouse,  '_tunecurvecorr_all.pdf']), '-dpdf', '-bestfit')

%%
%save all of the change scores from above for future use
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_d_scores.mat']), 'd_score_prefori_d1_d2', 'd_score_prefori_d2_d3', 'd_score_prefori_d3_d4', ...
    'd_score_prefori_d1_d3', 'd_score_prefori_d1_d4')

%save pref ori info
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_ori_changes.mat']), 'd1_prefori', 'd2_prefori', 'd3_prefori', 'd4_prefori')

%save tuned and matched cell IDs
save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_id_matches.mat']), 'match_all')


