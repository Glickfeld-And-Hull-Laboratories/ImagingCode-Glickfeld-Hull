%clear everything
clear all
clear all global
clc
%%
%find folders to load and experiment info

fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\2P\Arc_Multi_Day_Comparisons'; %folder to save files to
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\2P\Arc_Pooled_Multi_Day'; %folder to save files to
dataset = 'exp_list_arc_tjw'; %experiment list to pick files from
eval(dataset); %load dataset
arc_list = [7 10]; %mice from exp list to use
lacz_list = [1 4]; % mice from exp list to use

arc_scores = [];
lacz_scores = [];

%% HARD CODED
%arc mice

%mouse 1
mouse_arc_1 = expt(arc_list(1)).mouse;
img_area_arc_1 = expt(arc_list(1)).img_loc{1};
img_layer_arc_1 = expt(arc_list(1)).img_loc{2};
ses1_arc_1 = cell2mat(struct2cell(load(fullfile(newfnout, [mouse_arc_1 '_' img_area_arc_1 '_' img_layer_arc_1 '_' 'RW_d_scores']))));
ses2_arc_1 = cell2mat(struct2cell(load(fullfile(newfnout, [mouse_arc_1 '_' img_area_arc_1 '_' img_layer_arc_1 '_' 'RW_d_scores_2']))));
ses_all_arc_1 = [ses1_arc_1 ses2_arc_1];

%mouse 2
mouse_arc_2 = expt(arc_list(2)).mouse;
img_area_arc_2 = expt(arc_list(2)).img_loc{1};
img_layer_arc_2 = expt(arc_list(2)).img_loc{2};
ses1_arc_2 = cell2mat(struct2cell(load(fullfile(newfnout, [mouse_arc_2 '_' img_area_arc_2 '_' img_layer_arc_2 '_' 'RW_d_scores']))));
ses2_arc_2 = cell2mat(struct2cell(load(fullfile(newfnout, [mouse_arc_2 '_' img_area_arc_2 '_' img_layer_arc_2 '_' 'RW_d_scores_2']))));
ses_all_arc_2 = [ses1_arc_2 ses2_arc_2];

%combined arc
all_arc = [ses_all_arc_1 ses_all_arc_2];


%lacz mice

%mouse 1
mouse_lacz_1 = expt(lacz_list(1)).mouse;
img_area_lacz_1 = expt(lacz_list(1)).img_loc{1};
img_layer_lacz_1 = expt(lacz_list(1)).img_loc{2};
ses1_lacz_1 = cell2mat(struct2cell(load(fullfile(newfnout, [mouse_lacz_1 '_' img_area_lacz_1 '_' img_layer_lacz_1 '_' 'RW_d_scores']))));
ses2_lacz_1 = cell2mat(struct2cell(load(fullfile(newfnout, [mouse_lacz_1 '_' img_area_lacz_1 '_' img_layer_lacz_1 '_' 'RW_d_scores_2']))));
ses_all_lacz_1 = [ses1_lacz_1 ses2_lacz_1];

%mouse 2
mouse_lacz_2 = expt(lacz_list(2)).mouse;
img_area_lacz_2 = expt(lacz_list(2)).img_loc{1};
img_layer_lacz_2 = expt(lacz_list(2)).img_loc{2};
ses1_lacz_2 = cell2mat(struct2cell(load(fullfile(newfnout, [mouse_lacz_2 '_' img_area_lacz_2 '_' img_layer_lacz_2 '_' 'RW_d_scores']))));
ses2_lacz_2 = cell2mat(struct2cell(load(fullfile(newfnout, [mouse_lacz_2 '_' img_area_lacz_2 '_' img_layer_lacz_2 '_' 'RW_d_scores_2']))));
ses_all_lacz_2 = [ses1_lacz_2 ses2_lacz_2];

%combined mice
all_lacz = [ses_all_lacz_1 ses_all_lacz_2];

%%
%cdf plot
figure; cdfplot(arc_scores)
hold on
cdfplot(lacz_scores)
hold off
xlabel('Change in Pref Ori')
xlim([0 90]);
legend('Arc', 'LacZ')


%% LOOP
%arc loop
for iexp = arc_list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};
    ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_d_scores']));
    ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_d_scores_2']));
    ses_all = [ses1 ses2];
    ncells = length(ses_all);
    arc_scores = [arc_scores ses_all];
end

%lacz loop
for iexp = lacz_list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};
    ses1 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_d_scores']));
    ses2 = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_d_scores_2']));
    ses_all = [ses1 ses2];
    ncells = length(ses_all);
    lacz_scores = [lacz_scores ses_all];
end

%%
%cdf plot
figure; cdfplot(arc_scores)
hold on
cdfplot(lacz_scores)
hold off
xlabel('Change in Pref Ori')
xlim([0 90]);
legend(['Arc (n = ', num2str(length(arc_scores)), ')'], ['LacZ (n = ', num2str(length(lacz_scores)), ')'])
title('Pref Ori Changes for LacZ vs Arc (KRAB)')
print(fullfile(realfnout, ['all mice', '_change_ori.pdf']), '-dpdf', '-bestfit')

%% THIS NEEDS FIXED - I need to put d1-d2 and d1-d3 cdf plots separately and am having trouble
%with it
d1_d2_all = [];
d1_d3_all = [];
lacz_scores = struct2cell(lacz_scores);
for iexp = 1:length(lacz_scores)
    d1_d2 = lacz_scores{iexp,:,iexp}
    d1_d2_all = [d1_d2_all d1_d2]
end

%%