
%% V1 analysis for experiment inclusion

clc; clear all; close all;
doRedChannel = 1;
ds = 'CrossOriRandDirFourPhase_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = length(expt);
frame_rate = 15;
seed = rng;

max_dist = 5;

expt_good = [9 10 13 51 26 53 40];
expt_bad = [14 25 43 57];
expt_all = [expt_good, expt_bad];

for iexp = expt_all

mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc;
ImgFolder = expt(iexp).coFolder;
time = expt(iexp).coTime;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
fprintf([mouse ' ' date '\n'])

load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_respData.mat']))
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_stimData.mat']))

nCells = size(avg_resp_dir,1);
avg_mean=[];
avg_std=[];

for ic = 1:nCells
    [mc mi] = max(avg_resp_dir(ic,:,1,1,1));
    avg_mean = [avg_mean; mc];
    stdc = avg_resp_dir(ic,mi,1,1,2);
    avg_std = [avg_std; stdc];
end
avg_mean = mean(avg_mean);
avg_std = mean(avg_std);

avg_grat = mean(trialsperstim(:,1,1));
low_grat = min(trialsperstim(:,1,1));
avg_plaid = mean(reshape(trialsperstim(:,:,2),[],1));
low_plaid = min(reshape(trialsperstim(:,:,2),[],1));

%plot figure

if find(expt_good(:) == iexp) 
    w=2;
else; w=1; end

leg_str{iexp}=[num2str(iexp) ' n=' num2str(nCells)];

figure(1);
subplot(3,2,1)
    scatter(avg_mean,avg_std,'LineWidth',w)
    hold on
    text(avg_mean-0.004, avg_std+0.001, num2str(mouse),'FontSize',5)
    subtitle('grating responses')
    xlabel('mean df/f gratings')
    ylabel('mean std')
subplot(3,2,2)
    scatter(avg_grat,low_grat,'LineWidth',w)
    hold on
    subtitle('grating trials')
    xlabel('mean # of trials')
    ylabel('lowest # of trials')
    % legend(leg_str)
subplot(3,2,3)
    hold on
    scatter(avg_plaid,low_plaid,'LineWidth',w)
    subtitle('plaid trials')
    xlabel('mean # of trials')
    ylabel('lowest # of trials')
subplot(3,2,4)
    scatter(avg_grat, avg_std,'LineWidth',w)
    hold on
    subtitle('grating std by trial #')
    xlabel('mean # of trials')
    ylabel('mean std')
subplot(3,2,5)
    scatter(iexp, avg_mean,'LineWidth',w)
    hold on
    subtitle('mean df/f by expt number')
    xlabel('experiment number')
    ylabel('mean df/f gratings')    
subplot(3,2,6)
    scatter(iexp, nCells,'LineWidth',w)
    hold on
    subtitle('nCells by expt number')
    xlabel('experiment number')
    ylabel('nCells')    
end

% save figure
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';
figure(1); print(fullfile(outDir, [svName '_ExperimentInclusion_V1_prefdironly.pdf']),'-dpdf', '-fillpage') 



%% Higher visual area analysis for experiment inclusion

clc; clear all; close all;
doRedChannel = 1;
ds = 'CrossOriRandDirFourPhase_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = length(expt);
frame_rate = 15;
seed = rng;

max_dist = 5;

expt_good = [19 17 30];
expt_bad = [18 20 21 28 34 42 44];
expt_all = [expt_good, expt_bad];

for iexp = expt_all

mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc;
ImgFolder = expt(iexp).coFolder;
time = expt(iexp).coTime;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
fprintf([mouse ' ' date '\n'])

load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_respData.mat']))
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_stimData.mat']))

nCells = size(avg_resp_dir,1);
avg_mean=[];
avg_std=[];

for ic = 1:nCells
    avg_m = mean(reshape(avg_resp_dir(ic,:,1,1,1),1,[]));
    avg_mean = [avg_mean; avg_m];
    avg_s = mean(reshape(avg_resp_dir(ic,:,1,1,2),1,[]));
    avg_std = [avg_std; avg_s];
end
avg_mean = mean(avg_mean);
avg_std = mean(avg_std);

avg_grat = mean(trialsperstim(:,1,1));
low_grat = min(trialsperstim(:,1,1));
avg_plaid = mean(reshape(trialsperstim(:,:,2),[],1));
low_plaid = min(reshape(trialsperstim(:,:,2),[],1));

%plot figure

if find(expt_good(:) == iexp) 
    w=2;
else; w=1; end

leg_str{iexp}=[num2str(iexp) ' n=' num2str(nCells)];

figure(2);
subplot(3,2,1)
    scatter(avg_mean,avg_std,'LineWidth',w)
    hold on
    text(avg_mean-0.002, avg_std+0.001, num2str(mouse),'FontSize',5)
    subtitle('grating responses')
    xlabel('mean df/f gratings')
    ylabel('mean std')
subplot(3,2,2)
    scatter(avg_grat,low_grat,'LineWidth',w)
    hold on
    subtitle('grating trials')
    xlabel('mean # of trials')
    ylabel('lowest # of trials')
    % legend(leg_str)
subplot(3,2,3)
    hold on
    scatter(avg_plaid,low_plaid,'LineWidth',w)
    subtitle('plaid trials')
    xlabel('mean # of trials')
    ylabel('lowest # of trials')
subplot(3,2,4)
    scatter(avg_grat, avg_std,'LineWidth',w)
    hold on
    subtitle('grating std by trial #')
    xlabel('mean # of trials')
    ylabel('mean std')
subplot(3,2,5)
    scatter(iexp, avg_mean,'LineWidth',w)
    hold on
    subtitle('mean df/f by expt number')
    xlabel('experiment number')
    ylabel('mean df/f gratings')    
subplot(3,2,6)
    scatter(iexp, nCells,'LineWidth',w)
    hold on
    subtitle('nCells by expt number')
    xlabel('experiment number')
    ylabel('nCells')   
end

% save figure
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';
figure(2); print(fullfile(outDir, [svName '_ExperimentInclusion_HVAs.pdf']),'-dpdf', '-fillpage') 


