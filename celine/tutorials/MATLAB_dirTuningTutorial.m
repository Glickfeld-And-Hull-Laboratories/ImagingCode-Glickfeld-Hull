%setting my paths
fn_data_base = 'Z:\home\celine\Data\2p_data\';
fnout_base = 'Z:\home\celine\Analysis\';

%identifying animal and run
date = '210501';
ImgFolder = '001';
time = '1201';
mouse = 'CC05';
frame_rate = 60.8;
run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];
%%
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

nCells = size(raw_F_mat,3);

%convert to dF/F
F = raw_F_mat(61:120,:,:); %pull out the baseline window - here using the last 60 frames
F_avrg = mean(F,1); %take the average over the window
dF = raw_F_mat - F_avrg;
dF_overF = dF./F_avrg;

%define base_win and resp_win
base_win = dF_overF(61:120,:,:);
resp_win = dF_overF(120:180,:,:);

%collapse over the frames within each window, taking the average
base = transpose(squeeze(mean(base_win,1)));
resp = transpose(squeeze(mean(resp_win,1)));
%now arranged in cells X trials

drxns= cell2mat(input.tGratingDirectionDeg);

Udrxns = unique(drxns); %there are 16 unique directions
n_drxns = size(Udrxns,2);
%%
ttest_results = zeros(nCells, n_drxns);
avrg_response = zeros(nCells, n_drxns);
std_response = zeros(nCells, n_drxns);
se_response = zeros(nCells, n_drxns);



for i  = 1:n_drxns
    drxn1 = drxns(i); 
    inds = find(drxns == drxn1);
    this_base_win = base(:,inds);
    this_resp_win = resp(:,inds);
    [h,p] = ttest(this_base_win,this_resp_win, 'Dim',2,'alpha',0.05/7);    
    this_avrg = mean(this_resp_win,2);
    this_std = std(this_resp_win,0,2);

    ttest_results(:,i)=transpose(h);
    avrg_response(:,i)=transpose(this_avrg);
    std_response(:,i)=transpose(this_std);
    se_response(:,i)=transpose(this_std)/sqrt(length(inds));
end
clear this_base_win this_resp_win this_std this_avrg

nSig_cells = sum(max(ttest_results,[],2));

%% Von Mises fit for direction
%make baseline-subtracted response per cell
resp_sub_base = resp - base;

%find average baseline-sbutracted response of each cell to each drxn
avrg_resp_sub_base = zeros(nCells, n_drxns);

for i  = 1:n_drxns
    drxn1 = Udrxns(i)
    inds = find(drxns == drxn1);
    this_resp_sub_base = resp_sub_base(:,inds);
    this_avrg = mean(this_resp_sub_base,2);
    avrg_resp_sub_base(:,i)=this_avrg;
end
clear inds this_avrg this_resp_sub_base this_avrg drxn1


%%
vonMises_output = nan(nCells, 8);
%The inputs should be the drxnentations that you presented (theta, 1 x ndrxn) and the average response of each cell to each direction. 
for i = 1:nCells
    thisCell = avrg_resp_sub_base(i,:);
    [b_hat,k1_hat,R1_hat,R2_hat,u1_hat,u2_hat,sse,R_square] = miaovonmisesfit_dir(deg2rad(Udrxns),thisCell);
    vonMises_output(i,:)=[b_hat,k1_hat,R1_hat,R2_hat,u1_hat,u2_hat,sse,R_square];
end
clear thisCell b_hat k1_hat R1_hat u1_hat sse R_square %get rid of the last iteration of these looping variables
% trying to find the preferred drxn and vonM curve for each cell 

fit_drxns = [0:1:359];
y_fits = nan(nCells, size(fit_drxns,2));
pref_drxns = nan(nCells, 1);
for i = 1:nCells
    in = vonMises_output(i,:);
    b_tmp = in(1);
    k_tmp = in(2);
    R1_tmp = in(3);
    R2_tmp = in(4);
    u1_tmp = in(5);
    y_fits(i,:) =  b_tmp+R1_tmp.*exp(k_tmp.*(cos(deg2rad(fit_drxns)-u1_tmp)-1))+R2_tmp.*exp(k_tmp.*(cos(deg2rad(fit_drxns)-u1_tmp-pi)-1));
    pref_drxns(i)= fit_drxns(find(y_fits(i,:)==max(y_fits(i,:))));
    % sometimes there are two preferred directions - here I am taking the
    % first one (i.e. smaller degree value)
end
clear b_tmp k1_tmp R1_tmp u1_tmp

%% plotting

figure;
for i = 1:10
    subplot(5,2,i)
    %find how many oris each cell is significant for
    sig_drxns = strcat(string(sum(ttest_results(i,:))), " directions");
    %plotting - for one cell
    plot_cell_avrg = avrg_response(i,:);
    plot_se = se_response(i,:);
    plot_y_fit = y_fits(i,:);
    errorbar(Udrxns,plot_cell_avrg,plot_se);
    hold;
    plot(fit_drxns,plot_y_fit);
    vline(pref_drxns(i));
    title({"Significant response to",sig_drxns});
end

figure;
plot_i = 1;
for i = nCells-9:nCells
    subplot(5,2,plot_i)
    %find how many oris each cell is significant for
    sig_drxns = strcat(string(sum(ttest_results(i,:))), " directions");
    %plotting - for one cell
    plot_cell_avrg = avrg_response(i,:);
    plot_se = se_response(i,:);
    plot_y_fit = y_fits(i,:);
    errorbar(Udrxns,plot_cell_avrg,plot_se);
    hold;
    plot(fit_drxns,plot_y_fit);
    vline(pref_drxns(i));
    title({"Significant response to",sig_drxns});
    plot_i = plot_i+1;
end
clear sig_drxns  plot_i   plot_cell_avrg plot_se plot_y_fit

%% Von Mises bootstrap
reps = 10;
shuff_pref_drxn = nan(nCells,reps);



for r = 1:reps
    %make random sample
 
   %make random sample
    shuff_avrg_resp_sub_base = zeros(nCells, n_drxns);
    for i  = 1:n_drxns %for each direction
        drxn1 = Udrxns(i);
        OG_inds = find(drxns == drxn1);
        sub_n = round(size(OG_inds,2)*.8);
        inds=randsample(OG_inds,sub_n,1); %randomly select ~80% of the trials for this ori
        this_resp_sub_base = resp_sub_base(:,inds);
        this_avrg = mean(this_resp_sub_base,2);
        shuff_avrg_resp_sub_base(:,i)=this_avrg; 
    end
   
    
    %get VonM pref oris
    for i = 1:nCells
        thisCell = shuff_avrg_resp_sub_base(i,:);
        [b_tmp,k_tmp,R1_tmp,R2_tmp,u1_tmp,u2_tmp,sse,R_square] = miaovonmisesfit_dir(deg2rad(Udrxns),thisCell);
        this_fits =  b_tmp+R1_tmp.*exp(k_tmp.*(cos(deg2rad(fit_drxns)-u1_tmp)-1))+R2_tmp.*exp(k_tmp.*(cos(deg2rad(fit_drxns)-u1_tmp-pi)-1));
        this_pref = fit_drxns(min(find(this_fits==max(this_fits))));
        shuff_pref_drxn(i,r)=this_pref;
    end
      
end 
clear drxn1 drxn2 x inds this_avrg this_resp_sub_base this_avrg OG_inds thisCell b_hat k1_hat R1_hat u1_hat sse R_square this_fits this_pref %get rid of the last iteration of these looping variables
%% subtract original preferred ori from shuffled preferred ori and take abs value

absDiff_pref_drxn=abs(bsxfun(@minus, shuff_pref_drxn, pref_drxns));

% sort absDiff in acesnding order for each cell (row)
%collect a list of cell indeces that meet our criteria
%start with an empty list of cell
fit_cell_inds = [];
for i = 1:nCells
   thisCell = sort(absDiff_pref_drxn(i,:));
   if thisCell(:,(.9*reps)) < 22.5
       %add this i to the fit_cells list
       fit_cell_inds = [fit_cell_inds, i];
   end
end


% Summarize the responsivity and reliability of the imaged cells.
% a. 	Make a bar plot that has:
%   i. 	The number of segmented cells.
%   ii. 	The number of cells that were responsive to at least one direction.
%   iii. 	The number of cells that were reliably fit by a von Mises.

chart_mat = [nCells, nSig_cells, size(fit_cell_inds,2)];
labs = categorical({'Number of cells','Number significant','Number well fit'});
figure; bar(labs,chart_mat)
% b. 	Make a histogram of the distribution of preferred directions (u1) and sharpness (k1) of direction tuning for the reliably fit cells.

fit_cell_u1 = rad2deg(vonMises_output(fit_cell_inds,5));
hist(fit_cell_u1)
% xlim([0 360]);
title("Preferred direction of well fit cells");


fit_cell_k1 = vonMises_output(fit_cell_inds,2);
hist(fit_cell_k1, 20)
title("Sharpness of well fit cells");






