%%
clear all
clc
close all

%%
% dataset = 'exp_list_slc_long_tjw';
dataset = 'exp_list_ket_1wk_tjw';
eval(dataset); %load dataset
d1 = 4; %from expt list
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
nOn = input.nScansOn;
nOff = input.nScansOff;
nCells = size(data_dfof,2);

ntrials = size(input.tGratingDirectionDeg,2);
Dir = cell2mat_padded(input.tGratingDirectionDeg);
Dir = Dir(1:ntrials);
Dirs = unique(Dir);
nDirs = length(Dirs);

nori = length(Dirs)/2;
Ori = Dir;
for iori = 1:nori
    ind = find(Dir == Dirs(iori+nori));
    Ori(ind) = Dirs(iori);
end
Oris = unique(Ori);

%%
resp_win = 70:90;
a = squeeze(mean(data_dfof(resp_win,:,:),1));

%%

x = Ori;
y = a;
for ii = 1:nCells
z= sum(y(ii,:).*cos(2.*deg2rad(x))+y(ii,:).*i.*sin(2.*deg2rad(x)));
theta(ii) = 0.5*rad2deg(angle(z));
end

theta(find(theta<0)) = 180+theta(find(theta<0));

%%

figure;
scatter(prefOri(1,:),theta)
xlabel('von Mises')
ylabel('Vector Sum')

%% within session variability

[aaa aaa_idx] = datasample(a, length(a), 2 ,Replace=false);
aaa_oris = Ori(aaa_idx);

aaa_1st_set = aaa(:,1:length(aaa)/2);
aaa_oris_1st_set = aaa_oris(:,1:length(aaa_oris)/2);
aaa_2nd_set = aaa(:,(length(aaa)/2)+1:length(aaa));
aaa_oris_2nd_set = aaa_oris(:,(length(aaa_oris)/2)+1:length(aaa_oris));

%%

for ii = 1:nCells
    z_1 = sum(aaa_1st_set(ii,:).*cos(2.*deg2rad(aaa_oris_1st_set))+aaa_1st_set(ii,:).*i.*sin(2.*deg2rad(aaa_oris_1st_set)));
    theta_1(ii) = 0.5*rad2deg(angle(z_1));
end

for ii = 1:nCells
    z_2 = sum(aaa_2nd_set(ii,:).*cos(2.*deg2rad(aaa_oris_2nd_set))+aaa_2nd_set(ii,:).*i.*sin(2.*deg2rad(aaa_oris_2nd_set)));
    theta_2(ii) = 0.5*rad2deg(angle(z_2));
end

theta_1(find(theta_1<0)) = 180+theta_1(find(theta_1<0));
theta_2(find(theta_2<0)) = 180+theta_2(find(theta_2<0));

%%
%for traditional pref ori w/i session reliability, I need to sample data_tc
%and run the vonmises part of the day1 code for each 'half' of the data

%remember to save both 'halves' of each pref ori method

if exist('cellTCs_match','var') == 1
    npSub_tc = cellTCs_match{2};
end

down = 10;
nframes = size(npSub_tc,1)./down;
nCells = size(npSub_tc,2);
data_tc_down = squeeze(mean(reshape(npSub_tc, [down,nframes,nCells]),1));
tuningDownSampFactor = down;

[avgResponseEaOri_1,semResponseEaOri_1,vonMisesFitAllCells_1,fitReliability_1,R_square_1,tuningTC_1, ...
    avgResponseEaOri_2,semResponseEaOri_2,vonMisesFitAllCells_2,fitReliability_2,R_square_2,tuningTC_2] = ...
    getOriTuningTJ_datasplit(data_tc_down,input,tuningDownSampFactor,aaa_idx); 
    % vonMisesFitAllCells_1 = squeeze(vonMisesFitAllCellsAllBoots_1(:,1,:));
    % vonMisesFitAllCells_2 = squeeze(vonMisesFitAllCellsAllBoots_2(:,1,:));

[max_resp_1 prefOri_1] = max(vonMisesFitAllCells_1,[],1);
[max_resp_2 prefOri_2] = max(vonMisesFitAllCells_2,[],1);

prefOri_1 = squeeze(prefOri_1(:,1,:))';
prefOri_2 = squeeze(prefOri_2(:,1,:))';

%%

for iboot = 1:1000;
    [aaa aaa_idx] = datasample(a, length(a), 2 ,Replace=true);
    aaa_oris = Ori(aaa_idx);
        for ii = 1:nCells
            z = sum(aaa(ii,:).*cos(2.*deg2rad(aaa_oris))+aaa(ii,:).*i.*sin(2.*deg2rad(aaa_oris)));
            theta_resample(iboot,ii) = 0.5*rad2deg(angle(z));
        end
end
theta_resample(find(theta_resample<0)) = 180+theta_resample(find(theta_resample<0));

[lowerErr, upperErr, lowerCI, upperCI] = ciFromBoot(theta_resample,95);
ci_diff = abs(upperCI-lowerCI);
ci_diff(find(ci_diff>90)) = 180-ci_diff(find(ci_diff>90));
new_tune = find(ci_diff <= 45);

figure; histogram(ci_diff);

%%
save(fullfile(fnout, [date '_' mouse '_' ref_str '_new_pref.mat']),'theta', 'new_tune', 'theta_1', 'theta_2', 'prefOri_1', 'prefOri_2');
