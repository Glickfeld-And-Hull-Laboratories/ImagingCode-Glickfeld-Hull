%clear everything
clear all
clear all global
clc
close all

%%
dataset = 'exp_list_darklight_tjw';
eval(dataset); %load dataset
d1 = 8; %from expt list
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
save(fullfile(fnout, [date '_' mouse '_' ref_str 'newAvg.mat']),'newAvg')
