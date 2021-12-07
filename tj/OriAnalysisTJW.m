%% 
% PART 1
% 
% -----------------------------------------------------------------------------------------------------------------------
% 
% Transform the neuropil-subtracted time course from the tutorial (npSub_tc, 
% which is nframes x ncells) into a matrix that is nframes/trial x ntrials x ncells. 
% In the experiment used to collect this dataset, all trials are the same length 
% (nOff + nOn) and so you can use the reshape function to create this new trial 
% matrix. 
% 
% 
% 
% Clear Workspace

clear
clear all
clear global
%% 
% Load TC data

load('Z:\All_staff\home\lindsey\Analysis\2P\test\200825_i1328\200825_i1328_runs-003\200825_i1328_runs-003_TCs.mat')
load('Z:\All_staff\home\lindsey\Analysis\2P\test\200825_i1328\200825_i1328_runs-003\200825_i1328_runs-003_input.mat')
load('Z:\All_staff\home\lindsey\Analysis\2P\test\200825_i1328\200825_i1328_runs-003\200825_i1328_runs-003_mask_cell.mat')
load('Z:\All_staff\home\lindsey\Analysis\2P\test\200825_i1328\200825_i1328_runs-003\200825_i1328_runs-003_reg_shifts.mat')
load('Z:\All_staff\home\lindsey\Analysis\2P\test\200825_i1328\200825_i1328_runs-003\200825_i1328_runs-003_stimActFOV.mat')
%% 
% Get nOn, nOff, nFramesPerTrial, nTrials, and nCells for reshaping

nOn = input.nScansOn;
nOff = input.nScansOff;
nFramesPerTrial = nOn+nOff;
ntrials = size(input.tGratingDirectionDeg,2);
nCells = size(data_tc,2)
%% 
% nFramesPerTrial x nTrial x nCells Matrix

F_matrix = reshape(data_tc,[nFramesPerTrial,ntrials,nCells]);
%% 
% Potentially another way to do F0?

%firstcell = data_matrix(:,:,1)
%last20off = firstcell(101:120,:)
%avg20off = mean(last20off)

%for i = 1:nCells
%    cell_ind = data_matrix(:,:,i)
%    cell_last20_off = cell_ind(101:120,:)
%    cell_avg20off = mean(cell_last20_off)
%end 
%% 
% Now create a version of this matrix that is in dF/F instead of F.  You can 
% create an F0 matrix by taking the average of a window at the end of the Off 
% period on each trial (remember: each trial is Off then On. Then you need to 
% subtract F0 from the F matrix, and then divide by F0.
% 
% 
% 
% Plot of avg F per frame

avg_F_per_frame = mean(F_matrix,2);
avg_F_per_frame = mean(avg_F_per_frame,3);
plot(avg_F_per_frame)
%% 
% F0 Matrix

F0_matrix = mean(F_matrix(101:120,:,:));
%% 
% dF/F Matrix

dfOverF_matrix = (F_matrix - F0_matrix) ./ F0_matrix;
%% 
% ------------------------------------------------------------------------------------------------------------------------
% 
% PART 2
% 
% -----------------------------------------------------------------------------------------------------------------------
%% 
% Define two analysis windows- base_win and resp_win that reference the frames 
% on each trial for the baseline and response windows (simplest is to choose the 
% duration of the stimulus, and the same window before the stimulus; could also 
% get fancier and tailor your response window according to the delay in visual 
% response by looking at the average. It is also potentially important to select 
% a baseline window where the fluorescence has recovered following the previous 
% trial).
% 
% 
% 
% Define BL and Resp Windows

BL_window = [101:120];
Resp_window = [121:140];
%% 
% Create two matrices: base (a ncells x ntrials matrix with the average response 
% during the baseline window) and resp (a ncells x ntrials matrix with the average 
% response during the response window)

base = squeeze(mean(dfOverF_matrix(BL_window,:,:)));
resp = squeeze(mean(dfOverF_matrix(Resp_window,:,:)));
%% 
% Create a “for” loop where you loop through each orientation and select the 
% set of trials for that orientation and use a “ttest” to compare the baseline 
% and stimulus response across those trials. Output of ttest function should give 
% you the thresholded significance (using the bonferonni corrected alpha value) 
% for ncells x nori (this will be 1s for significant values). You should create 
% three new variables here – one that is ncells x nori for the outcome of the 
% ttest, and two that are the same size that have the average response of each 
% cell to each orientation and the standard error for that response (you will 
% use these last two later).
% 
% 
% 
% Get the list of orientations

Dir = celleqel2mat_padded(input.tGratingDirectionDeg);
nDirs = unique(Dir);
Ori = Dir;
Ori(Ori>=180) = Ori(Ori>=180) - 180;
Oris_list = unique(Ori);
nOri = length(Oris_list);
bonferroni = (.05/(nOri-1));
%% 
% For Loop

h = zeros(nOri,nCells);
p = zeros(nOri,nCells);
avg_resp_ori = zeros(nOri,nCells);
std_err_ori = zeros(nOri,nCells);

for i = 1:nOri
    Ori_ID = find(Ori == Oris_list(i));
    b = base(Ori_ID,:);
    r = resp(Ori_ID,:);
    [h(i,:) p(i,:)] = ttest(b,r,'Alpha',bonferroni);
    avg_resp_ori(i,:) = mean(r,1);
    std_err_ori(i,:) = (std(r,1))/(sqrt(length(Ori_ID)));
end
%% 
% How many cells are significant for at least one orientation?

sig_ori = zeros(1,nCells);

for i = 1:nCells
    sig_ori(i) = sum(h(:,i));
end

n_sig_ori = sum(sig_ori>=1)
%% 
% How many are significant without Bonferroni?

%75 are sig without it
%%
% ----------------------------------------------------------------------------------
% 
% PART 3
% 
% ----------------------------------------------------------------------------------
%% 
% 
% 
% Plot average tuning curve for all cells in FoV (each cell's response to each 
% orientation)

[n n2] = subplotn(nCells);
for i = 1:nCells 
    subplot(n,n2,i)
    errorbar(Oris_list,avg_resp_ori(:,i),std_err_ori(:,i))
    title(['# Sig Oris = ', num2str(sig_ori(i))])
end

%% 
% ____________________________________________________________
% 
% PART 4
% 
% _________________________________________________________
%% 
% Fit data with von Mises function for preferred orientation and orientation 
% selectivity
%% 
% New matrix of baseline subtracted responses for each trial (BL-resp)

%BL_sub_resp = resp-base
%% 
% Loop through each cell and use Von Mises function
% 
% inputs are: orientation presented (theta, 1xnOri) and avg response of each 
% cell to each ori
% 
% outputs are: b (BL), k1 (sharpness of tuning), u1 (preferred ori), sse (sum 
% squ error), and R2 -> should be one value for each cell
% 
% change theta from degrees to radians using degr2rad

clear b r
%%
Ori_rads = deg2rad(Oris_list);
k1 = zeros(1,nCells);
u1 = zeros(1,nCells);
sse = zeros(1,nCells);
rsquare = zeros(1,nCells);
b = zeros(1,nCells);
r = zeros(1,nCells);
%%
%[b k1 r u1 sse rsquare] = miaovonmisesfit_ori(Ori_rads,avg_resp_ori(:,1))


for i = 1:nCells
    [b(i) k1(i) r(i) u1(i) sse(i) rsquare(i)] = miaovonmisesfit_ori(Ori_rads,avg_resp_ori(:,i));
end


%% 
y_fits = zeros(180,nCells);
new_oris = deg2rad(0:1:179);
for i = 1:nCells
    b_tmp = b(i);
    R1_tmp = r(i);
    k1_tmp = k1(i);
    u1_tmp = u1(i);
    y_fits(:,i) = b_tmp+R1_tmp.*exp(k1_tmp.*(cos(2.*(new_oris-u1_tmp))-1));
end    

%% preferred orientation (in degrees)
pref_ori = rad2deg(u1);
%% plotting the fits

[n n2] = subplotn(nCells)
for i = 1:nCells 
    subplot(n,n2,i)
    hold on
    errorbar(Oris_list,avg_resp_ori(:,i),std_err_ori(:,i))
    plot(y_fits(:,i))
    vline(pref_ori(i),'black')
    title(['# Sig Oris = ', num2str(sig_ori(i))])
end


%% 
%%bootstrap
%bootstrap = 1000
%for i = 1:bootstrap
%    
%end


for i = 1:nOri
    Ori_ID_resample = randsample(find(Ori == Oris_list(i)),10,1);
    b_resample = base(Ori_ID_resample,:);
    r_resample = resp(Ori_ID_resample,:);
    avg_resp_ori_resample(i,:) = mean(r_resample,1);
   % std_err_ori_resample(i,:) = (std(r_resample,1))/(sqrt(length(Ori_ID)));
end
%% 

[n n2] = subplotn(nCells)
for i = 1:nCells 
    subplot(n,n2,i)
    errorbar(Oris_list,avg_resp_ori_resample(:,i),std_err_ori(:,i))
end