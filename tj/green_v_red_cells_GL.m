clear all
clc

%%
dataset = 'exp_list_arc_tjw'; %experiment list to pick files from
eval(dataset); %load dataset
d1 = 79; %from expt list
mouse = expt(d1).mouse; %mouse
run_str = 'runs-003'; %string on file name to load
img_area = expt(d1).img_loc{1};
img_layer = expt(d1).img_loc{2};
date = expt(d1).date; %day 1 (ref day) date
img_folder = expt(d1).runs; %img folder of day 1
time = expt(d1).time_mat;
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];
fnout = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\2P', datemouse, datemouserun); %folder to load files from

%%
%load day 1 data
CD = fnout; %finds current directory
cd(CD); %sets current directory
files = dir('*.mat');
 for k = 1 : length(files)
    baseFileName = files(k).name;
    fullFileName = fullfile(files(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    load(fullFileName);
 end
 clear files, clear baseFileName, clear fullFileName, clear k
 
%%
% Extracting green cell indices from analyzed data
ncells = 1:length(avgResponseEaOri);
greenCells = ncells(~ismember(ncells,redCells));

save(fullfile(fnout, [datemouserun '_multiday_alignment.mat']), 'greenCells', 'goodCells', 'okayCells', 'redCells','redChImg');

%%
%load day 2 data
d2 = 80;
date = expt(d2).date; %day 1 (ref day) date
img_folder = expt(d2).runs; %img folder of day 1
time = expt(d2).time_mat;
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];
fnout = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\2P', datemouse, datemouserun); %folder to load files from

load(fullfile(fnout, [datemouserun '_TCs.mat']))
load(fullfile(fnout, [datemouserun '_multiday_alignment.mat']))

match = find([cellImageAlign.pass]); %finds matched cell indices
red_match_ind = intersect(match, redCells);
green_match_ind = intersect(match, greenCells);

save(fullfile(fnout, [datemouserun '_tj_matches.mat']), 'match', 'red_match_ind', 'green_match_ind');

%%
%load day 3 data
d3 = 81;
date = expt(d3).date; %day 1 (ref day) date
img_folder = expt(d3).runs; %img folder of day 1
time = expt(d3).time_mat;
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];
fnout = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\2P', datemouse, datemouserun); %folder to load files from

load(fullfile(fnout, [datemouserun '_TCs.mat']))
load(fullfile(fnout, [datemouserun '_multiday_alignment.mat']))

match = find([cellImageAlign.pass]); %finds matched cell indices
red_match_ind = intersect(match, redCells);
green_match_ind = intersect(match, greenCells);

save(fullfile(fnout, [datemouserun '_tj_matches.mat']), 'match', 'red_match_ind', 'green_match_ind');