%setting my paths
fn_base = 'Z:\home\celine\Analysis\';
fnout = 'Z:\home\celine\Analysis\';

%identifying animal and run
date = '200825';
ImgFolder = '003';
time = '1248';
mouse = 'i1328';
frame_rate = 15.5;
run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];

%paste file name needed and load files
mWorks_fn = fullfile([date '_' mouse '_runs-' ImgFolder '_' 'input.mat']);
load(mWorks_fn);

tc_fn = fullfile([date '_' mouse '_runs-' ImgFolder '_' 'TCs.mat']);
load(tc_fn);


%find frames per trial
nOn = input.nScansOn;
nOff = input.nScansOff;
framesPerTrial = nOn+nOff; 
nTrials = size(npSub_tc, 1)/framesPerTrial;
nCells = size(npSub_tc, 2);

%reshape the timecourse into a frames/trial X trials X cells matric
raw_F_mat = reshape(npSub_tc,[framesPerTrial,nTrials,nCells]); %in this matric each item in the Z dimension is a single cell, and for each cell
%each column is a trial and the rows represent frames within that trial
%the values within those amtrices represent flourscence. 

%Now create a version of this matrix that is in dF/F instead of F.
%You can create an F0 matrix by taking the average of a window at the end of the Off period on each trial 
%(remember: each trial is Off then On. Then you need to subtract F0 from
%the F matrix, and then divide by F0).

% --> I am going to cut down the number of cells so that I can see the
% plots more clearly. 

%raw_F_mat = raw_F_mat(:,:,1:50);
nCells = size(raw_F_mat,3);

%convert to dF/F
F = raw_F_mat(61:120,:,:); %pull out the baseline window - here using the last 60 frames
F_avrg = mean(F,1); %take the average over the window
dF = raw_F_mat - F_avrg;
dF_overF = dF./F_avrg;

%define base_win and resp_win
base_win = dF_overF(61:120,:,:);
resp_win = dF_overF(120:180,:,:);

%collapse over the frames within each windown, taking the average
base = transpose(squeeze(mean(base_win,1)));
resp = transpose(squeeze(mean(resp_win,1)));
%now arranged in cells X trials

drxns= cell2mat(input.tGratingDirectionDeg);
%%
oris = unique(drxns);
oris = oris(1:8);
n_ori = length(oris);


%8 orientations. Each time, the angle increases by 22.5
ttest_results = zeros(nCells, n_ori);
avrg_response = zeros(nCells, n_ori);
std_response = zeros(nCells, n_ori);
se_response = zeros(nCells, n_ori);

drxn1 = 0; %start at 0

for i  = 1:n_ori
    drxn2 = drxn1 + 180;
    inds = find(drxns == drxn1 | drxns == drxn2);
    this_base_win = base(:,inds);
    this_resp_win = resp(:,inds);
    [h,p] = ttest(this_base_win,this_resp_win, 'Dim',2,'alpha',0.05/7); 
       
    this_avrg = mean(this_resp_win,2);
    this_std = std(this_resp_win,0,2);

    ttest_results(:,i)=transpose(h);
    avrg_response(:,i)=transpose(this_avrg);
    std_response(:,i)=transpose(this_std);
    se_response(:,i)=transpose(this_std)/sqrt(length(inds));
    
    drxn1=drxn1+22.5;
end
clear this_base_win this_resp_win this_std this_avrg

nSig_cells = sum(max(ttest_results,[],2));

%% Von Mises fit
%make baseline-subtracted response per cell
resp_sub_base = resp - base;

%find average baseline-sbutracted response of each cell to each orientation
%I can put this loop into the loop above
avrg_resp_sub_base = zeros(nCells, n_ori);
drxn1=0;
for i  = 1:n_ori
    drxn2 = drxn1 + 180;
    inds = find(drxns == drxn1 | drxns == drxn2);
    this_resp_sub_base = resp_sub_base(:,inds);
    this_avrg = mean(this_resp_sub_base,2);
    avrg_resp_sub_base(:,i)=this_avrg;
    drxn1=drxn1+22.5;
end
clear inds this_avrg this_resp_sub_base this_avrg

%loop through each cell
%apply von mises fxn to it

vonMises_output = zeros(nCells, 6);
%The inputs should be the orientations that you presented (theta, 1 x nori) and the average response of each cell to each orientation. 
for i = 1:nCells
    thisCell = avrg_resp_sub_base(i,:);
    [b_hat, k1_hat, R1_hat,u1_hat,sse,R_square] = miaovonmisesfit_ori(oris,thisCell);
    vonMises_output(i,:)=[b_hat, k1_hat, R1_hat,u1_hat,sse,R_square];
end
clear thisCell b_hat k1_hat R1_hat u1_hat sse R_square %get rid of the last iteration of these looping variables

fit_oris = [0:1:180];
y_fits = zeros(nCells, size(fit_oris,2));
pref_oris = nan(nCells, 1);
for i = 1:nCells
    in = vonMises_output(i,:);
    b_tmp = in(1);
    k1_tmp = in(2);
    R1_tmp = in(3);
    u1_tmp = in(4);
    y_fits(i,:) = b_tmp+R1_tmp.*exp(k1_tmp.*(cos(2.*(deg2rad(fit_oris)-u1_tmp))-1)); 
    pref_oris(i)= fit_oris(find(y_fits(i,:)==max(y_fits(i,:))));

end
clear b_tmp k1_tmp R1_tmp u1_tmp
%% I split this into two figures since there were so many subplots
figure;
for i = 1:10
    subplot(5,2,i)
    %find how many oris each cell is significant for
    sig_oris = strcat(string(sum(ttest_results(i,:))), " orientations");
    %plotting - for one cell
    plot_cell_avrg = avrg_response(i,:);
    plot_se = se_response(i,:);
    plot_y_fit = y_fits(i,:);
    errorbar(oris,plot_cell_avrg,plot_se);
    hold;
    plot(fit_oris,plot_y_fit);
    vline(pref_oris(i));
    title({"Significant response to",sig_oris});
end

figure;
plot_i = 1;
for i = nCells-9:nCells
    subplot(5,2,plot_i)
    %find how many oris each cell is significant for
    sig_oris = strcat(string(sum(ttest_results(i,:))), " orientations");
    %plotting - for one cell
    plot_cell_avrg = avrg_response(i,:);
    plot_se = se_response(i,:);
    plot_y_fit = y_fits(i,:);
    errorbar(oris,plot_cell_avrg,plot_se);
    hold;
    plot(fit_oris,plot_y_fit);
    vline(pref_oris(i));
    title({"Significant response to",sig_oris});
    plot_i = plot_i+1;
end
clear sig_oris  plot_i   plot_cell_avrg plot_se plot_y_fit 
%% Von Mises bootstrap
reps = 1000;
shuff_pref_ori = nan(nCells,reps);



for r = 1:reps
    %make random sample
    shuff_avrg_resp_sub_base = zeros(nCells, n_ori);
    drxn1=0;
    for i  = 1:n_ori %for each orientation
        drxn2 = drxn1 + 180;
        OG_inds = find(drxns == drxn1 | drxns == drxn2);
        inds=randsample(OG_inds,16,1); %randomly select ~80% of the trials for this ori
        this_resp_sub_base = resp_sub_base(:,inds);
        this_avrg = mean(this_resp_sub_base,2);
        shuff_avrg_resp_sub_base(:,i)=this_avrg; 
        drxn1=drxn1+22.5;
    end
   
    
    %get VonM pref oris
    for i = 1:nCells
        thisCell = shuff_avrg_resp_sub_base(i,:);
        [b_hat, k1_hat, R1_hat,u1_hat,sse,R_square] = miaovonmisesfit_ori(oris,thisCell);
        this_fits = b_hat+R1_hat.*exp(k1_hat.*(cos(2.*(deg2rad(fit_oris)-u1_hat))-1));
        this_pref = fit_oris(min(find(this_fits==max(this_fits))));
        shuff_pref_ori(i,r)=this_pref;
    end
      
end 
clear drxn1 drxn2 x inds this_avrg this_resp_sub_base this_avrg OG_inds thisCell b_hat k1_hat R1_hat u1_hat sse R_square this_fits this_pref %get rid of the last iteration of these looping variables

%% subtract original preferred ori from shuffled preferred ori and take abs value

absDiff_pref_ori=abs(bsxfun(@minus, shuff_pref_ori, pref_oris));

% sort absDiff in acesnding order for each cell (row)
%collect a list of cell indeces that meet our criteria
%start with an empty list of cell
fit_cell_inds = [];
for i = 1:nCells
   thisCell = sort(absDiff_pref_ori(i,:));
   if thisCell(:,(.9*reps)) < 22.5
       %add this i to the fit_cells list
       fit_cell_inds = [fit_cell_inds, i];
   end
end


% Summarize the responsivity and reliability of the imaged cells.
% a. 	Make a bar plot that has:
%   i. 	The number of segmented cells.
%   ii. 	The number of cells that were responsive to at least one orientation.
%   iii. 	The number of cells that were reliably fit by a von Mises.

chart_mat = [nCells, nSig_cells, size(fit_cell_inds,2)];
labs = categorical({'Number of cells','Number significant','Number well fit'});
figure; bar(labs,chart_mat)
% b. 	Make a histogram of the distribution of preferred orientations (u1) and sharpness (k1) of orientation tuning for the reliably fit cells.

fit_cell_u1 = vonMises_output(fit_cell_inds,4);
hist(fit_cell_u1, 8)
xlim([0 180]);
title("Preferred orientations of well fit cells");

fit_cell_k1 = vonMises_output(fit_cell_inds,2);
hist(fit_cell_k1, 20)
title("Sharpness of well fit cells");

%% "advanced"


% 2.	Try these analyses from a new dataset where stimuli were presented
% at different directions and spatial frequencies (mouse- i1328; date-
% 201129; time- 1514)- you can use singleChannelTC_SFxDirection for time
% course extraction. The goal here will be to make a plot showing how
% orientation tuning changes as a function of spatial frequency. This can
% be done by:
    % a.	Aligning all cells to their preferred orientation and then
    % averaging tuning curves for all cells for each spatial frequency. 
    %b.     Plotting the average tuning width (k from the fit) as a function of
    % spatial frequency.

