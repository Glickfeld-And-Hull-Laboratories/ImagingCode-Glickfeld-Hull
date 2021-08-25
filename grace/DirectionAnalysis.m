% basic code for doing direction analysis. Good for just getting a tuning
% curve and fit with bootstrapping but nothing fancy.
clear all
clear global
%% variables
%get path names
date = '200114';
ImgFolder = strvcat('003');
run = strvcat('000');
mouse = 'i1314';
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);

%open TC data
CD = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' run_str];
cd(CD);
imgMatFile = [date '_' mouse '_' run_str '_TCs.mat'];
load(imgMatFile);
imgMatFile2 = [date '_' mouse '_' run_str '_input.mat'];
load(imgMatFile2);

%create necessary variables
tGratingDir = celleqel2mat_padded(input.tGratingDirectionDeg);
dirs = unique(tGratingDir);
nDir = length(dirs);
nOn = input.nScansOn;
nOff = input.nScansOff;
nFrames = nOn+nOff;
nCells = size(npSub_tc,2);
nTrials = size(tGratingDir,2);

%transform into dF/F for each trial
trial_tc = reshape(npSub_tc,[nFrames nTrials nCells]);
trial_f = mean(trial_tc(nOff/2:nOff,:,:),1);
trial_dfof = (trial_tc-trial_f)./trial_f;

%% Identify cells that are significantly responsive to at least 1 direction

% Define analysis windows
resp_wind = nOff+1:nOff+nOn;
base_wind = 1+nOff-nOn:nOff;

% Create two matrices
dfof_resp = squeeze(mean(trial_dfof(resp_wind,:,:),1));
dfof_base = squeeze(mean(trial_dfof(base_wind,:,:),1));
dfof_subtract = dfof_resp - dfof_base;

% ttest for response significance
dfof_dir = zeros(nDir, nCells, 2);
h = zeros(nDir, nCells);
p = zeros(nDir, nCells);
for idir = 1:nDir
    ind = find(tGratingDir==dirs(idir));
    x = dfof_base(ind,:);
    y = dfof_resp(ind,:);
    [h(idir,:),p(idir,:)] = ttest(x,y,'dim',1,'Alpha',0.05./(nDir-1),'tail','left');
    dfof_dir(idir,:,1) = mean(y-x,1);
    dfof_dir(idir,:,2) = std(y-x,[],1)./sqrt(length(ind));
end
figure; 
% subplot(1,2,1)
imagesc(h');
xlabel('Direction')
ylabel('Cells')
figAxForm
% subplot(1,2,2)
% imagesc(isTunedVariable) %make sure this is size cells x 1
% xlabel('')
% ylabel('Cells')
% figAxForm
% subplt
if exist(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' run_str '\directionAnalysis'])
    print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' run_str '\directionAnalysis\sig_heatmap'],'-dpdf') 

else mkdir(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' run_str '\directionAnalysis'])
    print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' run_str '\directionAnalysis\sig_heatmap'],'-dpdf') 
end
%% Average tuning curves for all cells
for iCell = 1:nCells
    figure;
    errorbar(dirs,dfof_dir(:,iCell,1),dfof_dir(:,iCell,2))
    xlabel('Direction (deg)')
    ylabel('dF/F')
    suptitle(['Cell #' num2str(iCell)]) 
    if exist(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' run_str '\directionAnalysis\directionTuning'])
        print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' run_str '\directionAnalysis\directionTuning\directionTuning_Cell' num2str(iCell)],'-dpdf','-fillpage') 

    else mkdir(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' run_str '\directionAnalysis\directionTuning'])
        print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' run_str '\directionAnalysis\directionTuning\directionTuning_Cell' num2str(iCell)],'-dpdf','-fillpage')
    end
end
%% Fit data with von mises function
theta = deg2rad(dirs); 
thetafine = deg2rad(1:360); 
dfof_mean = squeeze(dfof_dir(:,:,1));
for iCell = 1:nCells
    try
    [b_hat(iCell),k1_hat(iCell),R1_hat(iCell),R2_hat(iCell),u1_hat(iCell),u2_hat(iCell),sse(iCell),R_square(iCell)] = miaovonmisesfit_dir(theta,dfof_mean(:,iCell));
    y_fit(iCell,:) = b_hat(iCell)+R1_hat(iCell).*exp(k1_hat(iCell).*(cos(thetafine-u1_hat(iCell))-1))+R2_hat(iCell).*exp(k1_hat(iCell).*(cos(thetafine-u2_hat(iCell))-1));
    catch
        b_hat(iCell) = NaN;
        k1_hat(iCell) = NaN;
        R1_hat(iCell) = NaN;
        R2_hat(iCell) = NaN;
        u1_hat(iCell) = NaN;
        u2_hat(iCell) = NaN;
        sse_hat(iCell) = NaN;
        R_square_hat(iCell)= NaN; 
        y_fit(iCell,:) = NaN(1,length(thetafine));
    end
end
%% plot to see the tuning
idx = find((R_square>0.5)'==1);
dfof_polar = dfof_dir;
for i_cell = 1:length(idx)
    figure;
    iCell=idx(i_cell);
    errorbar(dirs,dfof_dir(:,iCell,1),dfof_dir(:,iCell,2), 'LineStyle','none', 'Marker','o')
    hold on
    plot(1:360, y_fit(iCell,:),'r')
    xlabel('Direction (deg)')
    ylabel('dF/F')
%     suptitle(['Cell #' num2str(iCell)])
    saveas(
    print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' run_str '\directionAnalysis\vonmises' num2str(iCell)],'-dpdf')
%     figure;
%     if min(dfof_dir(:,iCell,1)) < 0
%         dfof_polar(:,iCell,1) = squeeze(dfof_dir(:,iCell,1))-squeeze(min(dfof_dir(:,iCell,1)));
%         polarplot([theta 2*pi],[(dfof_polar(:,iCell)); (dfof_polar(1,iCell))])
%         suptitle(['Cell #' num2str(iCell)])
%     else
%         polarplot([theta 2*pi],[squeeze(dfof_dir(:,iCell,1)); squeeze(dfof_dir(1,iCell,1))])
%         suptitle(['Cell #' num2str(iCell)])
%     end
%     print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' run_str '\directionAnalysis\polarplot' num2str(iCell)],'-dpdf')
end
%% saving structure for future polar plots
DT.y_fit = y_fit; 
DT.dfof_dir = dfof_dir; 
DT.dirs = dirs;
DT.idx = idx;
DT.dfof_polar = dfof_polar;
save(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' run_str '\' date '_' mouse '_' run_str '_DirTuning.mat'],'DT')
%% Reliability of tuning curves
% theta = deg2rad(dirs); 
% thetafine = deg2rad(1:360); 
% dfof_mean = squeeze(dfof_dir(:,:,1));
% nBootstraps = 1000;
% b_boot = nan(nCells, nBootstraps);
% k1_boot = nan(nCells, nBootstraps);
% R1_boot = nan(nCells, nBootstraps);
% R2_boot = nan(nCells, nBootstraps);
% u1_boot = nan(nCells, nBootstraps);
% u2_boot = nan(nCells, nBootstraps);
% sse_boot = nan(nCells, nBootstraps);
% R_square_boot = nan(nCells, nBootstraps);
% y_fit_boot = nan(nCells,length(thetafine),nBootstraps);
% for i = 1:1000
%     fprintf([num2str(i) '\n'])
%     ind_rand = randsample(nTrials, nTrials, 1);
%     dfof_rand = zeros(nDir, nCells);
%     for idir = 1:nDir
%         ind = find(tGratingDir==dirs(idir));
%         rand_trials = ismember(ind_rand,ind);
%         dfof_rand(idir,:) = mean(dfof_subtract(ind_rand(find(rand_trials)),:),1);
%     end
% tic    
% for iCell = 1:nCells
%         try
%         [b_boot(iCell,i),k1_boot(iCell,i),R1_boot(iCell,i),R2_boot(iCell,i),u1_boot(iCell,i),u2_boot(iCell,i),sse_boot(iCell,i),R_square_boot(iCell,i)] = miaovonmisesfit_dir(theta,dfof_rand(:,iCell));
%         y_fit_boot(iCell,:,i) = b_boot(iCell,i)+R1_boot(iCell,i).*exp(k1_boot(iCell,i).*(cos(thetafine-u1_boot(iCell,i))-1))+R2_boot(iCell,i).*exp(k1_boot(iCell,i).*(cos(thetafine-u2_boot(iCell,i))-1));
%         catch 
%         b_boot(iCell,i) = NaN;
%         k1_boot(iCell,i) = NaN;
%         R1_boot(iCell,i) = NaN;
%         R2_boot(iCell,i) = NaN;
%         u1_boot(iCell,i) = NaN;
%         u2_boot(iCell,i) = NaN;
%         sse_boot(iCell,i) = NaN;
%         R_square_boot(iCell,i)= NaN; 
%         y_fit_boot(iCell,:,i) = repmat(NaN,1,length(thetafine));
%         end
%         
% end 
% toc
% end

%% Reliability of tuning curves within direction
theta = deg2rad(dirs); 
thetafine = deg2rad(1:360); 
dfof_mean = squeeze(dfof_dir(:,:,1));
nBootstraps = 1000;
b_boot = nan(nCells, nBootstraps);
k1_boot = nan(nCells, nBootstraps);
R1_boot = nan(nCells, nBootstraps);
R2_boot = nan(nCells, nBootstraps);
u1_boot = nan(nCells, nBootstraps);
u2_boot = nan(nCells, nBootstraps);
sse_boot = nan(nCells, nBootstraps);
R_square_boot = nan(nCells, nBootstraps);
y_fit_boot = nan(nCells,length(thetafine),nBootstraps);
for i = 1:1000
    fprintf([num2str(i) '\n'])
    dfof_rand = zeros(nDir, nCells);
    for idir = 1:nDir
        ind = find(tGratingDir==dirs(idir));
        ind_rand = randsample(length(ind), length(ind), 1);
        dfof_rand(idir,:) = mean(dfof_subtract(ind(ind_rand),:),1);
    end   
for iCell = 1:nCells
        try
        [b_boot(iCell,i),k1_boot(iCell,i),R1_boot(iCell,i),R2_boot(iCell,i),u1_boot(iCell,i),u2_boot(iCell,i),sse_boot(iCell,i),R_square_boot(iCell,i)] = miaovonmisesfit_dir(theta,dfof_rand(:,iCell));
        y_fit_boot(iCell,:,i) = b_boot(iCell,i)+R1_boot(iCell,i).*exp(k1_boot(iCell,i).*(cos(thetafine-u1_boot(iCell,i))-1))+R2_boot(iCell,i).*exp(k1_boot(iCell,i).*(cos(thetafine-u2_boot(iCell,i))-1));
        catch 
        b_boot(iCell,i) = NaN;
        k1_boot(iCell,i) = NaN;
        R1_boot(iCell,i) = NaN;
        R2_boot(iCell,i) = NaN;
        u1_boot(iCell,i) = NaN;
        u2_boot(iCell,i) = NaN;
        sse_boot(iCell,i) = NaN;
        R_square_boot(iCell,i)= NaN; 
        y_fit_boot(iCell,:,i) = repmat(NaN,1,length(thetafine));
        end
        
end 
end
%% reliability of tuning for each cell
idx = find((R_square>0.5)'==1);
pref_ori = zeros(nCells,nBootstraps);
spot_900 = nan(1,nCells);
reliably_fit_cells = nan(1,nCells);
for i_cell = 1:length(idx)
    iCell=idx(i_cell);
    pref_ori(iCell,:) = abs(u1_boot(iCell,:) - u1_hat(iCell));
    pref_ori_sort = sort(pref_ori(iCell,:));
    no_nan = pref_ori_sort(:,~isnan(pref_ori_sort));
    spot_900(1,iCell) = prctile(no_nan,90);
    reliably_fit_cells(1,iCell) = spot_900(1,iCell) <= deg2rad(22.5);
    figure;
%     subplot(2,1,1)
    iCell=idx(i_cell);
    histogram(rad2deg(u1_boot(iCell,:)),'FaceColor','r','BinWidth',22.5);
    figXAxis([ ],'Pref Direction (Degrees)',[0 360],0:45:360,0:45:360);
    ylabel('nBootstraps')
    suptitle(['Cell ' num2str(iCell) ' - 90% btsps within ' num2str(rad2deg(spot_900(1,iCell))) ' degrees'])
%     subplot(2,1,2)
%     histogram(k1_boot(iCell,:));
%     axis square
    print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' run_str '\directionAnalysis\bootstraps' num2str(iCell)],'-dpdf')
end

num_of_resp_cells = length(idx);
num_of_reliable_cells = nansum(reliably_fit_cells(:,:));
figure;
x2 = [1 2 3];
y2 = [nCells num_of_resp_cells num_of_reliable_cells];
cellnames = {'Segmented Cells'; 'Responsive Cells'; 'Reliably Fit Cells';};
bar(x2,y2,'r');
set(gca,'xticklabel',cellnames)
ylabel('Number of cells')
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' run_str '\directionAnalysis\barplot'],'-dpdf')

%% Multiple Comparison Test

