% Extracting k value and max response from analyzed data
clear all
clc

%%
dataset = 'exp_list_arc_tjw'; %experiment list to pick files from
eval(dataset); %load dataset
d1 = 3; %day 1 in expt list
mouse = expt(d1).mouse; %mouse
ref_str = 'runs-003'; %string on file name to load
img_area = expt(d1).img_loc{1};
img_layer = expt(d1).img_loc{2};
date = expt(d1).date; %day 1 (ref day) date
img_folder = expt(d1).runs; %img folder of day 1
time = expt(d1).time_mat;
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' ref_str];
fnout = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\2P', datemouse, datemouserun); %folder to load files from


%%
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
Dir = celleqel2mat_padded(input.tGratingDirectionDeg); %directions of each trial
nDirs = unique(Dir); %how many directions
Ori = Dir;
Ori(Ori>=180) = Ori(Ori>=180) - 180; %turn Dirs into Oris by subtracting 180 from those over 180
Oris_list = unique(Ori); %what are the Oris
nOri = length(Oris_list); %how many Oris
nCells = size(avgResponseEaOri,1);

Ori_rads = deg2rad(Oris_list);
k1 = zeros(1,nCells);
u1 = zeros(1,nCells);
sse = zeros(1,nCells);
rsquare = zeros(1,nCells);
b = zeros(1,nCells);
r = zeros(1,nCells);

for i = 1:nCells
    [b(i) k1(i) r(i) u1(i) sse(i) rsquare(i)] = miaovonmisesfit_ori(Ori_rads,avgResponseEaOri(i,:));
end

max_dfof = squeeze(max(vonMisesFitAllCellsAllBoots,[],1));
max_dfof = max_dfof(1,:);

save(fullfile(fnout, [date '_' mouse '_' ref_str '_k_and_max_vals.mat']),'k1','max_dfof')


%%
