close all
clear all
%path names
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff';
th_fn = fullfile(fn_base, 'home\jerry');
data_fn = fullfile(th_fn, 'matlab_tut');
mworks_fn = fullfile(fn_base, 'Behavior\Data');
fnout = fullfile(th_fn, 'matlab_tut\Analysis');

%experiment info
date = '221129';
ImgFolder = '003';
time = '1526';
mouse = 'i2080';
frame_rate = 15.5;
run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];
%% 

% Load Time Course Data
load(fullfile(mworks_fn, ['data-' mouse '-' date '-' time '.mat'])) %from mworks to extract trial info
fname = fullfile(fnout, datemouse, datemouserun, [datemouserun '_TCs.mat']);
load(fname); %this loads in time course data np_tc npSub_tc and data_tc
ntrials = size(input.tGratingDirectionDeg,2) %each unique value is one trial

% Reshape time course data and calculate dF/F
% Reshape frames x cells to frames/trial x trials x cells 
twtc = reshape(npSub_tc,[],ntrials,size(npSub_tc,2)); % i.e. trialwise time course
% dF/F
nOn = input.nScansOn;
nOff = input.nScansOff;
twtc_f0 = mean(twtc(nOff/2:nOff,:,:),1);
twtc_df = twtc - twtc_f0;
twtc_dfof = twtc_df ./ twtc_f0;

%% 
% find out the specfics of the grating directions
Dir = celleqel2mat_padded(input.tGratingDirectionDeg); %transforms cell array into matrix (1 x ntrials)
Dirs = unique(Dir);
nDirs = length(Dirs);
% find cell response dynamics by qualitatively looking at avearge of all cells in all trials
tc_avg = mean(mean(twtc_dfof(:,:,:),3),2); %avg across cells then trials
tc_avg = squeeze(tc_avg');
plot(tc_avg)
% define response windows
base_win = 31:60; %WOULD BE BETTER NOT TO HARD CODE THESE, AND INSTEAD HAVE THEM BE RELATIVE TO NoN AND NoFF, EG resp_win = nOff+2:nOn+nOff;
resp_win = 61:90;
base = squeeze(mean(twtc_dfof(base_win,:,:),1))'; % calculate avg dfof in the base window...
resp = squeeze(mean(twtc_dfof(resp_win,:,:),1))'; % ...then transpose to get ncell x ntrials

% now convert the grating directions to orientations
ori = Dir;
ori(Dir>=180) = ori(Dir>=180) - 180;

% prep unique stim list
nStim = length(unique(ori));
Stims = unique(ori);
% Initiate empty arrays and for loops through unique orientations to compare average resp vs. base activity level (averaged over frames)
% n = # of trials with the specified orientation in each for loop
results = []; %IDEALLY WHEN INITIALIZING YOU SHOULD MAKE A MATRIX THAT IS ALREADY THE CORRECT DIMENTIONS; THIS WAY 1) IT RUNS FASTER AND 
% 2) YOU KNOW IF SOMETHING IS WRONG BECAUSE YOU GET AN ERROR THAT THE
% DIMENSIONS DON'T FIT. USUALLY YOU ONLY WANT TO INITIALIZE WITH [] IF YOU
% HAVE A TYPE OF DATA THAT WONT GO INTO A MATRIX WELL OR YOU REALLY DON'T
% KNOW WHAT THE FINAL DIMENSIONS WILL BE / YOU WANT IT TO BE VERY FELXIBLE,
% SUCH AS IF YOUR CONCATENATING OVER DIFFERENT SESSIONS.
pvals = [];
ses = [];
for i = 1:nStim
    ind = find(ori == Stims(i)); % find indices of trials with this ori
    [result,pval,ci,stats] = ttest(base(:,ind),resp(:,ind),'Alpha',0.05/(length(unique(ori))),'Dim',2); %"(length(unique(ori)))" CAN BE REPLACED WITH nStims; 
    results = [results result];
    pvals = [pvals pval];
    sd = stats.sd;
    se = sd/sqrt(length(ind));
    ses = [ses se];
end

% how many cells are significant for at least one orientation
orisig_cells = [];
for i = 1:length(results(:,1))
    if sum(results(i,:)) > 0
        orisig_cells = [orisig_cells i];
    else
        continue
    end
end

%YOU CAN REPLACE THE LOOP ABOVE WITH A SINGLE LINE:
%orisig_cells=logical(sum(results,2)). THIS DOES GIVE YOU A LOGICAL LIST OF
%CELL IDS FOR CELLS THAT ARE RESPONSIVE, RATHER THAN A LIST OF THE ID
%NUMBERS. OFTEN A LOGICAL IS EASIER TO DEAL WITH BUT IF YOU WANT A LIST OF
%THE ID NUBERS YOU CAN USE "FIND"

tot_nCells = length(base(:,1));
fprintf([num2str(length(orisig_cells)) ' out of ' num2str(tot_nCells) ' cells are responsive to at least one orientation\n'])

%% 
% Average tuning curve for all cells in the FOV
% (i.e. each cell's stimulus response to each orientation)
% find average dfof in each orientation for every cell (output should be nCell x nOri)
ori_dfof = []; 
for i = 1:nStim
    ind = find(ori == Stims(i));
    eachori_dfof = mean(resp(:,ind),2) - mean(base(:,ind),2);
    ori_dfof = [ori_dfof eachori_dfof];
end
clear eachori_dfof

% plot them
[n, n2] = subplotn(tot_nCells);
% initiate smoother index 

t_curve_idx = (0:1:179);
k1 = [];
u1 = [];

close all
figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:tot_nCells
    subplot(n,n2,i);
    hold on
    errorbar(Stims,ori_dfof(i,:),ses(i,:));
    xticks(Stims);
    title(num2str(sum(results(i,:))));
    ylabel('df/f');
    %von Mises functionin
    %"response - baseline dfof" = ori_dfof
    [b_tmp,k1_tmp,R1_tmp,u1_tmp,sse,r_sqr] = miaovonmisesfit_ori(deg2rad(Stims),ori_dfof(i,:));
    y_fit = b_tmp+R1_tmp.*exp(k1_tmp.*(cos(2.*(deg2rad(t_curve_idx)-u1_tmp))-1));
    k1(i) = k1_tmp;
    u1(i) = u1_tmp;
    plot(y_fit);
    [m,pref_ori(i)] = max(y_fit);
end
clear m
hold off

preserve_ori = pref_ori;
%%
% Bootstrapping
pref_ori = preserve_ori; %so this section can be re-run

nBoot = 1000;
boot_pref_ori = zeros(tot_nCells,nBoot); % preallocation for performance 
% boot_pref_ori is nCells x nIterations

for v = 1:nBoot
    ori_iter = [];
    for i = 1:nStim %DO YOU NEED TO DO THIS FOR EVERY LOOP? ISN'T IT THE SAME EACH TIME?
        ind = find(ori == Stims(i)); % find all trials with current ori
        ind_rd = randsample(ind,round(length(ind)/4*3),true);
        eachori = mean(resp(:,ind_rd),2) - mean(base(:,ind_rd),2);
        ori_iter = [ori_iter eachori]; %#ok<AGROW>
    end
    for j = 1:tot_nCells
        [b_tmp,k1_tmp,R1_tmp,u1_tmp,sse,r_sqr] = miaovonmisesfit_ori(deg2rad(Stims),ori_iter(j,:));
        y_fit = b_tmp+R1_tmp.*exp(k1_tmp.*(cos(2.*(deg2rad(t_curve_idx)-u1_tmp))-1));
        [val, tmp_idx] = max(y_fit); %#ok<SAGROW>
        boot_pref_ori(j,v) = tmp_idx - 1;
    end
    if mod(v,50) == 0 % indicate progress every 50 loops
        fprintf(['Bootstrapping is at the ' num2str(v) 'th iteration.\n'])
    end    
end

%transpose pref_ori to nCell x 1
pref_ori = pref_ori';

%duplicate the column nBoot times for later array-wide arithmetic 
%   calculation because every cell has the same
%   non-bootstrapped pref_ori
%pref_ori_big is nCell x nIterations
pref_ori_big = repmat(pref_ori,1,nBoot);

%boot_pref_ori is nCell x nIterations
%absoridiff is a custom function in my git folder
pref_ori_diff = arrayfun(@absoridiff,pref_ori_big,boot_pref_ori);

sorted = sort(pref_ori_diff,2);
keep_cells = find(sorted(:,nBoot*0.9)<22.5); %indices of cells with 90th percentile deviation less than 22.5

figure;
y = [tot_nCells length(orisig_cells) length(keep_cells)];
bar(y);
ylim([0 90]);
xticklabels({"segmented","responsive","reliable"});

figure;
subplot(1,2,1);
histogram(rad2deg(u1(keep_cells)),15);
title('preferred orientations') %this really shouldn't be an almost normal distribution but I can't tell where I did wrong..
subplot(1,2,2);
histogram(k1(keep_cells),15)
title('sharpness')
