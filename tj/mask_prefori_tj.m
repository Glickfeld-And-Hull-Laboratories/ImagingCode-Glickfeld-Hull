%%
clear all
clc
close all

%%
dataset = 'exp_list_arc_tjw';
eval(dataset); %load dataset
d1 = 2; %from expt list
mouse = expt(d1).mouse; %mouse
ref_str = ['runs-' expt(d1).runs];
img_area = expt(d1).img_loc{1};
img_layer = expt(d1).img_loc{2};
date = expt(d1).date; %day 1 (ref day) date
img_folder = expt(d1).runs; %img folder of day 1
time = expt(d1).time_mat;
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' ref_str];
fnout = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P', datemouse, datemouserun); %folder to load files from

%%
%load files
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
prefori = prefOri(1,:);
cell_ids = nonzeros(unique(mask_cell))';
mask_prefs = mask_cell;

for i = 1:length(cell_ids)
    if cellImageAlign(i).pass == 1 && ismember(i, ind_theta90)
        mask_prefs(mask_prefs == i) = prefori(i);
    else
        mask_prefs(mask_prefs == i) = 0;
    end
end

%%
figure;
imagesc(mask_prefs);
colorbar