%%  REMEMBER TO LABEL GRAPHS
% PART 1 - Transform TC into trial structure that is also in dF/F
% 
% -----------------------------------------------------------------------------------------------------------------------
% 
% Transform the neuropil-subtracted time course from the tutorial (npSub_tc, 
% which is nframes x ncells) into a matrix that is nframes/trial x ntrials x ncells. 
% In the experiment used to collect this dataset, all trials are the same length 
% (nOff + nOn) and so you can use the reshape function to create this new trial 
% matrix. 
% 
% We would like the data in dF/F rather than F (how much F changes with stimulus), and we would like them
% arranged in such a way that we can tell when each ori is presented
% 
% Clear Workspace

clear
clear all
clear global
%% 
% Load TC data

%path info
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff'; %base file name
tj_fn = fullfile(fn_base, 'home\tj'); %personal file path
data_fn = fullfile(tj_fn, '2P_Imaging\tutorial'); 
mworks_fn = fullfile(fn_base, 'Behavior\Data'); %behavior data
fnout = fullfile(tj_fn, 'Analysis\Analysis\2P\tutorial\211020_tj_092021'); %analysis folder

%experiment info
date = '211020';
ImgFolder = '002';
time = '0917';
mouse = 'tj_092021';
frame_rate = 15.5;
run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];

%load relevant files (.mat) based on CD
CD = fullfile(fnout, datemouserun); %finds current directory
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
% Get nOn, nOff, nFramesPerTrial, nTrials, and nCells for reshaping

nOn = input.nScansOn; %use input data
nOff = input.nScansOff; %use input data
nFramesPerTrial = nOn+nOff;
ntrials = size(input.tGratingDirectionDeg,2); %this is 1x160 (ntrials)
nCells = size(npSub_tc,2); %nFrames x nCells
%% 
% nFramesPerTrial x nTrial x nCells Matrix reshaping

F_matrix = reshape(npSub_tc,[nFramesPerTrial,ntrials,nCells]);

%% 
% Now create a version of this matrix that is in dF/F instead of F.  You can 
% create an F0 matrix by taking the average of a window at the end of the Off 
% period on each trial (remember: each trial is Off then On. Then you need to 
% subtract F0 from the F matrix, and then divide by F0.
% 
% 
% 
% Plot of avg F per frame to see optimal baseline/response windows

avg_F_per_frame = mean(F_matrix,2); %avg across trials
avg_F_per_frame = mean(avg_F_per_frame,3); %avg across cells
plot(avg_F_per_frame)
%% 
% F0 Matrix

F0_matrix = mean(F_matrix(30:59,:,:)); %avg across last part of OFF frames
%this part could be written in a NON-hard-coded manner where it finds the
%min and max values and takes frames from around those values
%% 
% dF/F Matrix

dfOverF_matrix = (F_matrix - F0_matrix) ./ F0_matrix; %subtract baseline from F and divide
plot(mean(dfOverF_matrix(:,:,1),2)) %an example of the type of data we can extract from this part
 %this is dF/F averaged across trials for each frame of cell 1
%% 
% ------------------------------------------------------------------------------------------------------------------------
% 
% PART 2 - Identify cells significantly responsive to one or more
% orientations
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
% We have visually-responsive cells in dF/F and now we want to see which
% ones are responsive to particular orientations
% 
% Define BL and Resp Windows

BL_window = [30:59];
Resp_window = [61:90];
%% 
% Create two matrices: base (a ntrials x ncells matrix with the average response 
% during the baseline window) and resp (a ntrials x nncells matrix with the average 
% response during the response window)

base = squeeze(mean(dfOverF_matrix(BL_window,:,:))); %avg across OFF frames selected
resp = squeeze(mean(dfOverF_matrix(Resp_window,:,:))); %avg across ON frames selected

%this plot just gives an example of one cell's base and resp periods over
%trials
subplot(1,2,1)
plot(base(:,1))
subplot(1,2,2)
plot(resp(:,1))

%% 
% Create a �for� loop where you loop through each orientation and select the 
% set of trials for that orientation and use a �ttest� to compare the baseline 
% and stimulus response across those trials. Output of ttest function should give 
% you the thresholded significance (using the bonferonni corrected alpha value) 
% for ncells x nori (this will be 1s for significant values). You should create 
% three new variables here � one that is ncells x nori for the outcome of the 
% ttest, and two that are the same size that have the average response of each 
% cell to each orientation and the standard error for that response (you will 
% use these last two later).
% 
% 
% 
% Get the list of orientations

Dir = celleqel2mat_padded(input.tGratingDirectionDeg); %directions of each trial
nDirs = unique(Dir); %how many directions
Ori = Dir;
Ori(Ori>=180) = Ori(Ori>=180) - 180; %turn Dirs into Oris by subtracting 180 from those over 180
Oris_list = unique(Ori); %what are the Oris
nOri = length(Oris_list); %how many Oris
bonferroni = (.05/(nOri-1)); %bonferroni correction for multiple comparisons to minimize type I error
%% 
% For Loop

h = zeros(nOri,nCells);
p = zeros(nOri,nCells);
avg_resp_ori = zeros(nOri,nCells);
std_err_ori = zeros(nOri,nCells);

for i = 1:nOri
    Ori_ID = find(Ori == Oris_list(i)); %array of trial numbers for specific ori
    b = base(Ori_ID,:); %df/f for BL window on trials of specific ori
    r = resp(Ori_ID,:); %df/f for Resp window on trials of specific ori
    [h(i,:) p(i,:)] = ttest(b,r,'Alpha',bonferroni,'tail','left'); %t test of BL vs Resp means for given ori
    avg_resp_ori(i,:) = mean(r,1); %avg response of each cell at each ori
    std_err_ori(i,:) = (std(r,1))/(sqrt(length(Ori_ID))); %std error of each cell at each ori
end
%% 
% How many cells are significant for at least one orientation?

sig_ori = zeros(1,nCells);

for i = 1:nCells
    sig_ori(i) = sum(h(:,i)); %sum of all direction significance levels (0 or 1) for each cell
end

n_sig_ori = sum(sig_ori>=1) %how many cells are >=1 (at least one significant ori)
%% 
% How many are significant without Bonferroni?

for i = 1:nOri
    Ori_ID = find(Ori == Oris_list(i));
    b = base(Ori_ID,:); 
    r = resp(Ori_ID,:);
    [h_2(i,:) p_2(i,:)] = ttest(b,r,'tail','left'); 
end

sig_ori_2 = zeros(1,nCells);

for i = 1:nCells
    sig_ori_2(i) = sum(h(:,i)); 
end

n_sig_ori_2 = sum(sig_ori_2>=1)
%%
% ----------------------------------------------------------------------------------
% 
% PART 3 - tuning curves to show responsivity of each cell to each
% orientation
% 
% ----------------------------------------------------------------------------------
%% 
% 
% 
% Plot average tuning curve for all cells in FoV (each cell's response to each 
% orientation)

[n n2] = subplotn(nCells); %finds optimal # of subplots to make
for i = 1:nCells 
    subplot(n,n2,i) %which subplot position
    errorbar(Oris_list,avg_resp_ori(:,i),std_err_ori(:,i)) %plot with error bars of all oris for particular cell
    title(['# Sig Oris = ', num2str(sig_ori(i))])
end

%% 
% ____________________________________________________________
% 
% PART 4 - Fit data with von Mises function for preferred orientation and
% orientation selectivity
% 
% _________________________________________________________
%% 
% The von Mises is a bell shaped curve that can fit the data and provide
% parameters in its output such as preferred ori for stimuli that were not
% necessarily shown -> make predictions
%% 
% New matrix of baseline subtracted responses for each trial (BL-resp)

BL_sub_resp = resp-base;
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

%this will loop through each cell with the von Mises fit function
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
[pref_ori_val pref_ori_deg] = max(y_fits,[],1);
%pref_ori_deg2 = rad2deg(u1); this one makes values close to 180 negative
%and is using u1 rather than y_fits
%pref_ori_deg3 = pref_ori_deg2; pref_ori_deg3(find(pref_ori_deg2<0)) = pref_ori_deg2(find(pref_ori_deg2<0))+180;
%% plotting the fits

[n n2] = subplotn(nCells);
for i = 1:nCells ;
    subplot(n,n2,i);
    hold on;
    errorbar(Oris_list,avg_resp_ori(:,i),std_err_ori(:,i)); %plot original data again
    plot(y_fits(:,i)); %plot the fits
    vline(pref_ori_deg(i),'black'); %make a line where pref ori is
    title(['# Sig Oris = ', num2str(sig_ori(i))]);
end;
hold off;
%% %
%Part 5 - Bootstrapping to determine reliability of fits

%% Just because we are able to fit a function to the orientation tuning data doesnt mean that its a good fit;
% not all cells are tuned and they will not be well-fit; we want to
% determine reliability of the fit for likelihood that the measured pref
% ori is the true pref ori -> how variable are the tuning curves when we
% sample from subsets of trials many times?

%%bootstrap
bootstrap = 10;

avg_resp_ori_resample = zeros(nOri,nCells);
std_err_ori_resample = zeros(nOri,nCells);
y_fits_resample = zeros(length(new_oris),nCells);
pref_ori_resample = zeros(bootstrap,nCells);




for boot = 1:bootstrap
    
    for i = 1:nOri %this is resampling random trials for each ori
        Ori_ID_resample = randsample(find(Ori == Oris_list(i)),length(find(Ori == Oris_list(i))),1); %changed 10 to ntrials
        %b_resample = base(Ori_ID_resample,:);
        r_resample = resp(Ori_ID_resample,:);
        avg_resp_ori_resample(i,:) = mean(r_resample,1);
        %std_err_ori_resample(i,:) = (std(r_resample,1))/(sqrt(length(Ori_ID_resample)));
    end
    
    for c = 1:nCells ;
        [b(c) k1(c) r(c) u1(c) sse(c) rsquare(c)] = miaovonmisesfit_ori(Ori_rads,avg_resp_ori_resample(:,c));
        b_tmp = b(c);
        R1_tmp = r(c);
        k1_tmp = k1(c);
        u1_tmp = u1(c);
        y_fits_resample(:,c) = b_tmp+R1_tmp.*exp(k1_tmp.*(cos(2.*(new_oris-u1_tmp))-1));
        pref_ori_resample(boot,c) = rad2deg(u1(c)); %could put this loop outside with (boot,:) instead of (boot:c)
        pref_ori_resample(pref_ori_resample<0) = pref_ori_resample(pref_ori_resample<0)+180;
    end;
    

    
end

%[pref_ori_val_resample(boot,c) pref_ori_deg_resample(boot,c)] = max(y_fits_resample(boot,c),[],1)
%[pref_ori_val pref_ori_deg] = max(y_fits,[],1)
%[pref_ori_val_resample pref_ori_deg_resample] = max(y_fits_resample,[],1)
%still has negative numbers come up
%% Determine reliability of pref ori by subtracting original pref oris from the bootstrapped
%pref oris on each bootstrap and take absolute value; sort in ascending
%order, and get the 90th %ile (value in 900 position)

pref_ori_90prct = prctile(sort(abs(pref_ori_resample - pref_ori_deg)),90);
well_fit_cells = pref_ori_90prct<22.5;
n_well_fit = sum(well_fit_cells);

%% Part 6 - Summarize responsivity and reliability of the cells

%bar plot of # of cells, number of responsive cells, and number of cells
%reliably fit by von Mises

x_lab = categorical(["Total", "Sig to at least one ori", "Well fit"])
x_lab = reordercats(x_lab,["Total", "Sig to at least one ori", "Well fit"])
bar(x_lab, [nCells, n_sig_ori, n_well_fit])
ylabel("Number of Cells")

%% 

%histogram of distribution of pref ori and sharpness of tuning for reliably
%fit cells

pref_ori_wellfit = pref_ori_deg(well_fit_cells);
sharpness_well_fit = k1(well_fit_cells);
subplot(1,2,1)
histogram(pref_ori_wellfit, 8)
subplot(1,2,2)
histogram(sharpness_well_fit, 8)
