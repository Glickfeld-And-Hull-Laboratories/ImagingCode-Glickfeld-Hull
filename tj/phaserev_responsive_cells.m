%%
%clear everything
clear all
clear all global
clc
close all

%%
%find folders to load and experiment info

fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\phaserev\multi_day'; %folder to save files to
dataset = 'exp_list_phaserev_tjw'; %experiment list to pick files from
eval(dataset); %load dataset

%%

%we start with loading the proper information based on the mouse's info in
%the experiment list

%baseline1 to baseline1+4d comparison
d1 = 1; %base1 in expt list
d2 = 2; %base1+4d in expt list
d3 = 3;
d4 = 4;
d5 = 5;
d6 = 6;
d7 = 7;
mouse = expt(d1).mouse; 
ref_str_d1 = ['runs-',expt(d1).runs];
ref_str_d2 = ['runs-',expt(d2).runs]; 
ref_str_d3 = ['runs-',expt(d3).runs]; 
ref_str_d4 = ['runs-',expt(d4).runs]; 
ref_str_d5 = ['runs-',expt(d5).runs]; 
ref_str_d6 = ['runs-',expt(d6).runs]; 
ref_str_d7 = ['runs-',expt(d7).runs]; 
img_area = expt(d1).img_loc{1};
img_layer = expt(d1).img_loc{2};
date_d1 = expt(d1).date; 
date_d2 = expt(d2).date; 
date_d3 = expt(d3).date; 
date_d4 = expt(d4).date; 
date_d5 = expt(d5).date; 
date_d6 = expt(d6).date; 
date_d7 = expt(d7).date; 


%load multiday data for day 2
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));
d3_matches = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'multiday_alignment.mat']));
d4_matches = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'multiday_alignment.mat']));
d5_matches = load(fullfile(fnout, [date_d5 '_' mouse], [date_d5 '_' mouse '_' ref_str_d5], [date_d5 '_' mouse '_' ref_str_d5 '_' 'multiday_alignment.mat']));
d6_matches = load(fullfile(fnout, [date_d6 '_' mouse], [date_d6 '_' mouse '_' ref_str_d6], [date_d6 '_' mouse '_' ref_str_d6 '_' 'multiday_alignment.mat']));
d7_matches = load(fullfile(fnout, [date_d7 '_' mouse], [date_d7 '_' mouse '_' ref_str_d7], [date_d7 '_' mouse '_' ref_str_d7 '_' 'multiday_alignment.mat']));


%load TCs
d1_TCs = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'TCs.mat']));
d2_TCs = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'TCs.mat']));
d3_TCs = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'TCs.mat']));
d4_TCs = load(fullfile(fnout, [date_d4 '_' mouse], [date_d4 '_' mouse '_' ref_str_d4], [date_d4 '_' mouse '_' ref_str_d4 '_' 'TCs.mat']));
d5_TCs = load(fullfile(fnout, [date_d5 '_' mouse], [date_d5 '_' mouse '_' ref_str_d5], [date_d5 '_' mouse '_' ref_str_d5 '_' 'TCs.mat']));
d6_TCs = load(fullfile(fnout, [date_d6 '_' mouse], [date_d6 '_' mouse '_' ref_str_d6], [date_d6 '_' mouse '_' ref_str_d6 '_' 'TCs.mat']));
d7_TCs = load(fullfile(fnout, [date_d7 '_' mouse], [date_d7 '_' mouse '_' ref_str_d7], [date_d7 '_' mouse '_' ref_str_d7 '_' 'TCs.mat']));

load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'input.mat']));

%% matching indices
match_d2 = find([d2_matches.cellImageAlign.pass]); 
match_d3 = find([d3_matches.cellImageAlign.pass]); 
match_d4 = find([d4_matches.cellImageAlign.pass]); 
match_d5 = find([d5_matches.cellImageAlign.pass]); 
match_d6 = find([d6_matches.cellImageAlign.pass]); 
match_d7 = find([d7_matches.cellImageAlign.pass]); 

nOn = input.nScansOn;
nOff = input.nScansOff;
nTrials = input.stopAfterNTrials;

%%
%d2
tc_match = d2_TCs.cellTCs_match{1,2}(:,match_d2);
nCells = size(tc_match,2);
data_trial = permute(reshape(tc_match,[nOff + nOn nTrials nCells]),[1 3 2]);
data_f = mean(data_trial(nOff-15:nOff,:,:),1);
data_dfof = (data_trial-data_f)./data_f;

% h = zeros(nCells);
% p = zeros(nCells);

base_win = nOff-15:nOff;
resp_win = nOff+1:nOn;
base = squeeze(mean(data_dfof(base_win,:,:),1)); %avg across base
resp = squeeze(mean(data_dfof(resp_win,:,:),1));

ptrials = nTrials-1;
psig = double(0.05/ptrials);

for i = 1:nCells


                [h(i,:), p(i,:)] = ttest(resp(i,:), base(i,:),'tail','right');
                
end

sig_cell = zeros(1,nCells);

for i = 1:nCells
    sig_cell(i) = sum(h(i,:)); %sum of all direction significance levels (0 or 1) for each cell
end

n_sig_cell_d2 = sum(sig_cell>=1);
perc_sig_cell_d2 = n_sig_cell_d2/nCells;
%%
%d4
tc_match = d4_TCs.cellTCs_match{1,2}(:,match_d4);
nCells = size(tc_match,2);
data_trial = permute(reshape(tc_match,[nOff + nOn nTrials nCells]),[1 3 2]);
data_f = mean(data_trial(nOff-15:nOff,:,:),1);
data_dfof = (data_trial-data_f)./data_f;

% h = zeros(nCells);
% p = zeros(nCells);

base_win = nOff-15:nOff;
resp_win = nOff+1:nOn;
base = squeeze(mean(data_dfof(base_win,:,:),1)); %avg across base
resp = squeeze(mean(data_dfof(resp_win,:,:),1));

ptrials = nTrials-1;
psig = double(0.05/ptrials);

for i = 1:nCells


                [h(i,:), p(i,:)] = ttest(resp(i,:), base(i,:),'tail','right');
                
end

sig_cell = zeros(1,nCells);

for i = 1:nCells
    sig_cell(i) = sum(h(i,:)); %sum of all direction significance levels (0 or 1) for each cell
end

n_sig_cell_d4 = sum(sig_cell>=1);
perc_sig_cell_d4 = n_sig_cell_d4/nCells;

%%
%d5
tc_match = d5_TCs.cellTCs_match{1,2}(:,match_d5);
nCells = size(tc_match,2);
data_trial = permute(reshape(tc_match,[nOff + nOn nTrials nCells]),[1 3 2]);
data_f = mean(data_trial(nOff-15:nOff,:,:),1);
data_dfof = (data_trial-data_f)./data_f;

% h = zeros(nCells);
% p = zeros(nCells);

base_win = nOff-15:nOff;
resp_win = nOff+1:nOn;
base = squeeze(mean(data_dfof(base_win,:,:),1)); %avg across base
resp = squeeze(mean(data_dfof(resp_win,:,:),1));

ptrials = nTrials-1;
psig = double(0.05/ptrials);

for i = 1:nCells


                [h(i,:), p(i,:)] = ttest(resp(i,:), base(i,:),'tail','right');
                
end

sig_cell = zeros(1,nCells);

for i = 1:nCells
    sig_cell(i) = sum(h(i,:)); %sum of all direction significance levels (0 or 1) for each cell
end

n_sig_cell_d5 = sum(sig_cell>=1);
perc_sig_cell_d5 = n_sig_cell_d5/nCells;


%%
%d6
tc_match = d6_TCs.cellTCs_match{1,2}(:,match_d6);
nCells = size(tc_match,2);
data_trial = permute(reshape(tc_match,[nOff + nOn nTrials nCells]),[1 3 2]);
data_f = mean(data_trial(nOff-15:nOff,:,:),1);
data_dfof = (data_trial-data_f)./data_f;

% h = zeros(nCells);
% p = zeros(nCells);

base_win = nOff-15:nOff;
resp_win = nOff+1:nOn;
base = squeeze(mean(data_dfof(base_win,:,:),1)); %avg across base
resp = squeeze(mean(data_dfof(resp_win,:,:),1));

ptrials = nTrials-1;
psig = double(0.05/ptrials);

for i = 1:nCells


                [h(i,:), p(i,:)] = ttest(resp(i,:), base(i,:),'tail','right');
                
end

sig_cell = zeros(1,nCells);

for i = 1:nCells
    sig_cell(i) = sum(h(i,:)); %sum of all direction significance levels (0 or 1) for each cell
end

n_sig_cell_d6 = sum(sig_cell>=1);
perc_sig_cell_d6 = n_sig_cell_d6/nCells;

%%
%d1 and d7 need different info than above
input_d1 = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'input.mat']));
input_d7 = load(fullfile(fnout, [date_d7 '_' mouse], [date_d7 '_' mouse '_' ref_str_d7], [date_d7 '_' mouse '_' ref_str_d7 '_' 'input.mat']));

nOn = input_d1.input.nScansOn;
nOff = input_d1.input.nScansOff;
nTrials = input_d1.input.stopAfterNTrials;

ori_d1 = celleqel2mat_padded(input_d1.input.tGratingDirectionDeg);
ori_d7 = celleqel2mat_padded(input_d7.input.tGratingDirectionDeg);
oris = unique(ori_d1);

d1_45_id = find(ori_d1 == 45);
d7_45_id = find(ori_d7 == 45);

%%
%d1
tc_match = d1_TCs.npSub_tc;
nCells = size(tc_match,2);
data_trial = permute(reshape(tc_match,[nOff + nOn nTrials nCells]),[1 3 2]);
data_f = mean(data_trial(nOff-15:nOff,:,:),1);
data_dfof = (data_trial-data_f)./data_f;

% h = zeros(nCells);
% p = zeros(nCells);

base_win = nOff-15:nOff;
resp_win = nOff+1:nOn;
base = squeeze(mean(data_dfof(base_win,:,:),1)); %avg across base
resp = squeeze(mean(data_dfof(resp_win,:,:),1));

for i = 1:nCells


                [h(i,:), p(i,:)] = ttest(resp(i,d1_45_id), base(i,d1_45_id),'tail','right');
                
end

sig_cell = zeros(1,nCells);

for i = 1:nCells
    sig_cell(i) = sum(h(i,:)); %sum of all direction significance levels (0 or 1) for each cell
end

n_sig_cell_d1 = sum(sig_cell>=1);
perc_sig_cell_d1 = n_sig_cell_d1/nCells;

%%
%d7
tc_match = d7_TCs.cellTCs_match{1,2}(:,match_d7);
nCells = size(tc_match,2);
data_trial = permute(reshape(tc_match,[nOff + nOn nTrials nCells]),[1 3 2]);
data_f = mean(data_trial(nOff-15:nOff,:,:),1);
data_dfof = (data_trial-data_f)./data_f;

% h = zeros(nCells);
% p = zeros(nCells);

base_win = nOff-15:nOff;
resp_win = nOff+1:nOn;
base = squeeze(mean(data_dfof(base_win,:,:),1)); %avg across base
resp = squeeze(mean(data_dfof(resp_win,:,:),1));


for i = 1:nCells


                [h(i,:), p(i,:)] = ttest(resp(i,d7_45_id), base(i,d7_45_id),'tail','right');
                
end

sig_cell = zeros(1,nCells);

for i = 1:nCells
    sig_cell(i) = sum(h(i,:)); %sum of all direction significance levels (0 or 1) for each cell
end

n_sig_cell_d7 = sum(sig_cell>=1);
perc_sig_cell_d7 = n_sig_cell_d7/nCells;

%%
all_percs = [perc_sig_cell_d1 perc_sig_cell_d2 perc_sig_cell_d4 perc_sig_cell_d5 perc_sig_cell_d6 perc_sig_cell_d7];

figure;
bar(all_percs*100)
ylim([0,100])
ylabel('% of cells')
xticklabels({'Day 1 8ori', 'Day 1 SRP', 'Day 3 SRP', 'Day 4 SRP', 'Day 5 SRP', 'Day 5 8ori'})
title(['Cells Siginificantly Responsive to 45 Degrees (', mouse, ')'])

print(fullfile(newfnout, [mouse,  '_sig_cells.pdf']), '-dpdf', '-bestfit')

