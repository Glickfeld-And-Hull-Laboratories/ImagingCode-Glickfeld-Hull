% Extracting k value and max response from analyzed data
clear all
clc
close all

%%
%dataset = 'exp_list_arc_tjw'; %experiment list to pick files from
% dataset = 'exp_list_tjw'; %experiment list to pick files from
% dataset = 'exp_list_darklight_tjw';
dataset = 'exp_list_darklight_actual_tjw';
eval(dataset); %load dataset
d1 = 14; %from expt list
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
%%
dir_mat = cell2mat(input.tGratingDirectionDeg);
ori_mat = dir_mat;
ori_mat(find(dir_mat>=180)) = ori_mat(dir_mat>=180)-180;
oris = unique(ori_mat);
nOri = length(oris);
dirs = unique(dir_mat);
nDir = length(dirs);
nCells = length(avgResponseEaOri);

%%

[maxResp prefOri_ind] = max(avgResponseEaOri,[],2);
newAvg = zeros(nCells,nOri);
for iCell = 1:nCells
    if prefOri_ind(iCell)<4
        index = 4 - prefOri_ind(iCell);
        newAvg(iCell,:) = avgResponseEaOri(iCell,[8-index+1:8 1:8-index]);
    elseif prefOri_ind(iCell)>4
        index = prefOri_ind(iCell) - 4;
        newAvg(iCell,:) = avgResponseEaOri(iCell,[index+1:8 1:index]);
    elseif prefOri_ind(iCell)==4
        newAvg(iCell,:) = avgResponseEaOri(iCell,:);
    end
end

%%
figure;
fast_errbar(1:8,newAvg,1,'color',[0.3 0 0.5]);
fix_axes(gcf,16,'Orientation','dF/F'); axis square

%%
save(fullfile(fnout, [date '_' mouse '_' ref_str '_newAvg.mat']),'newAvg')
