% %% Make expected trial distribution
% clc; clear all; close all;
% doRedChannel = 0;
% ds = 'CrossOriRandDirFourPhase_ExptList_SG';
% rc = behavConstsAV;
% eval(ds)
% nexp = length(expt);
% 
% iexp = 38;
% max_dist = 10;
% 
% frame_rate = 15;
% seed = rng;
% 
% mouse = expt(iexp).mouse;
% date = expt(iexp).date;
% area = expt(iexp).img_loc;
% ImgFolder = expt(iexp).coFolder;
% time = expt(iexp).coTime;
% nrun = length(ImgFolder);
% run_str = catRunName(cell2mat(ImgFolder), nrun);
% 
% base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
% 
% fprintf([mouse ' ' date '\n'])
% 
% load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_respData.mat']))
% load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_stimData.mat']))
% 
% ntrials = reshape(trialsperstim(:,:,2),[],1);
% 
% mean = mean(ntrials);
% std = std(ntrials);
% 
% pd = makedist('Normal','mu',mean,'sigma',std);
% % mean =  13.9167
% % std = 3.4073
% 
% % trialsperstimFOUR = trialsperstim;

%% Load one phase data
clearvars -except trialsperstimFOUR
clc;close all; 
doRedChannel = 0;
ds = 'CrossOriRandDirFourPhase_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = length(expt);

iexp = 112; 
max_dist = 10;

frame_rate = 15;
seed = rng;

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc;
ImgFolder = expt(iexp).coFolder;
time = expt(iexp).coTime;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';

fprintf([mouse ' ' date '\n'])

%% Pref direction analysis
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']))

%%
if doRedChannel == 0
    red_cells = [];
end

prewin_frames = unique(celleqel2mat_padded(input.tItiWaitFrames))./3;
nFramesOn = unique(celleqel2mat_padded(input.nStimOneFramesOn));
postwin_frames = unique(celleqel2mat_padded(input.nStimOneFramesOn));
tt = (1-prewin_frames:postwin_frames).*(1/frame_rate);
data_resp = nan(prewin_frames+postwin_frames,nCells,nTrials-1);
data_f = nan(1,nCells,nTrials-1);

for itrial = 1:nTrials
    if cStimOn(itrial) + postwin_frames < sz(3)
        data_resp(:,:,itrial) = npSub_tc(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
    end
end

data_f = mean(data_resp(1:prewin_frames,:,:),1);
data_dfof_tc = (data_resp-data_f)./data_f;

ind_stimAlone = intersect(intersect(find(stimCon_all),find(maskCon_all==0)),find(centroid_dist<max_dist));
ind_maskAlone = intersect(intersect(find(stimCon_all==0),find(maskCon_all)),find(centroid_dist<max_dist));
ind_plaid = intersect(intersect(find(stimCon_all),find(maskCon_all)),find(centroid_dist<max_dist));
ind_blank = intersect(find(stimCon_all==0),find(maskCon_all==0));
ind_p = cell(1,nMaskPhas);
% maskPhas_shuf = maskPhas_all(randperm(length(maskPhas_all)));
% ind_p_shuf = cell(1,nMaskPhas);
for ip = 1:nMaskPhas
    ind_p{ip} = intersect(find(maskPhas_all == maskPhas(ip)),find(centroid_dist<max_dist));
%     ind_p_shuf{ip} = intersect(find(maskPhas_shuf == maskPhas(ip)),find(centroid_dist<max_dist));
end

resp_cell = cell(nStimDir,nMaskPhas,2);
trialsperstim = zeros(nStimDir,nMaskPhas,2);
trialInd = cell(nStimDir,nMaskPhas,2);
h_resp =zeros(nCells,nStimDir,nMaskPhas,2);
p_resp =zeros(nCells,nStimDir,nMaskPhas,2);
resp_win = prewin_frames+5:prewin_frames+nFramesOn;
resp_blank = squeeze(mean(data_dfof_tc(prewin_frames+frame_rate:prewin_frames+nFramesOn,:,ind_blank),1));
avg_resp_dir = zeros(nCells, nStimDir,nMaskPhas, 2, 2);
% avg_resp_dir_shuf = zeros(nCells, nStimDir,nMaskPhas, 2, 2);
all_resp_dir = [];
all_resp_plaid = cell(1,nMaskPhas);
all_dir = [];
all_plaid = cell(1,nMaskPhas);
nStim = nStimDir;
for iDir = 1:nStimDir
    ind_stimdir = intersect(find(stimDir_all == stimDirs(iDir)),find(centroid_dist<max_dist));
    ind_maskdir = intersect(find(maskDir_all == stimDirs(iDir)),find(centroid_dist<max_dist));
    ind_diralone = [intersect(ind_stimdir, ind_stimAlone), intersect(ind_maskdir, ind_maskAlone)];
    ind_dirplaid = [intersect(ind_stimdir, ind_plaid)];
    trialsperstim(iDir,1,1) = length(ind_diralone);
    trialInd{iDir,1,1} = ind_diralone;
    resp_cell{iDir,1,1} = squeeze(mean(data_dfof_tc(resp_win,:,ind_diralone),1));
    [h_resp(:,iDir,1,1), p_resp(:,iDir,1,1)] = ttest2(resp_cell{iDir,1,1},resp_blank,'dim',2,'tail','right','alpha', 0.05./nStim);
    avg_resp_dir(:,iDir,1,1,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind_diralone),1),3));
    avg_resp_dir(:,iDir,1,1,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,ind_diralone),1),[],3)./sqrt(length(ind_diralone)));
    all_resp_dir = [all_resp_dir squeeze(mean(data_dfof_tc(resp_win,:,ind_diralone),1))];
    all_dir = [all_dir iDir.*ones(size(ind_diralone))];
    for ip = 1:nMaskPhas
        ind_dpplaid = intersect(ind_dirplaid,ind_p{ip});
%         ind_dpplaid_shuf = intersect(ind_dirplaid,ind_p_shuf{ip});
        trialsperstim(iDir,ip,2) = length(ind_dpplaid);
        trialInd{iDir,ip,2} = ind_dpplaid;
        resp_cell{iDir,ip,2} = squeeze(mean(data_dfof_tc(resp_win,:,ind_dpplaid),1));
        [h_resp(:,iDir,ip,2), p_resp(:,iDir,ip,2)] = ttest2(resp_cell{iDir,ip,2},resp_blank,'dim',2,'tail','right','alpha', 0.05./nStim);
        avg_resp_dir(:,iDir,ip,2,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind_dpplaid),1),3));
        avg_resp_dir(:,iDir,ip,2,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,ind_dpplaid),1),[],3)./sqrt(length(ind_dpplaid)));
%         avg_resp_dir_shuf(:,iDir,ip,2,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind_dpplaid_shuf),1),3));
        if iDir == 1
            all_resp_plaid{ip} = [];
            all_plaid{ip} = [];
        end
        all_resp_plaid{ip} = [all_resp_plaid{ip} squeeze(mean(data_dfof_tc(resp_win,:,ind_dpplaid),1))];
        all_plaid{ip} = [all_plaid{ip} iDir.*ones(size(ind_dpplaid))];
    end
end

%find number of trials per stim
ntrial_p = [];
ntrial_g = [];
resps = double.empty(nCells,0);
avg_resp_dir_shuf = zeros(nCells, nStimDir, nMaskPhas);
for id = 1:nStimDir
    ntrial_p = [ntrial_p; size(resp_cell{id,1,2},2)];
    ntrial_g = [ntrial_g; size(resp_cell{id,1,1},2)];
end



resp_ind = find(sum(sum(sum(h_resp,2),3),4));
resp_ind_dir = find(sum(h_resp(:,:,1,1),2)); %sig responsive to gratings
resp_ind_plaid = find(sum(sum(h_resp(:,:,:,2),2),3));
p_anova_dir = zeros(1,nCells);
p_anova_plaid = zeros(2,nCells);
for iCell = 1:nCells
    p_anova_dir(iCell) = anova1(all_resp_dir(iCell,:), all_dir, 'off'); %direction selective to gratings
    for ip = 1:nMaskPhas
        p_anova_plaid(ip,iCell) = anova1(all_resp_plaid{ip}(iCell,:), all_plaid{ip}, 'off'); %direction selective to plaids
    end
end

p_dir = find(p_anova_dir<0.05);
p_plaid1 = find(p_anova_plaid(1,:)<0.05);

for iCell = 1:nCells
    [max_val max_ind] = max(avg_resp_dir(iCell,:,1,1,1));
    null_ind = max_ind+(nStimDir./2);
    null_ind(find(null_ind>nStimDir)) = null_ind(find(null_ind>nStimDir))-nStimDir;
    min_val = avg_resp_dir(iCell,null_ind,1,1,1);
    if min_val < 0; min_val = 0; end
    DSI(iCell) = (max_val-min_val)./(max_val+min_val);
end

DSI_ind = find(DSI>0.5); %direction selective to gratings
OSI_ind = find(DSI<0.5);

if ~exist(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' max_dist]),'dir')
    mkdir(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' num2str(max_dist)]))
end
save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_respData.mat']), 'resp_cell', 'data_dfof_tc', 'resp_ind', 'frame_rate', 'h_resp', 'avg_resp_dir','p_anova_dir','p_anova_plaid');
save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_stimData.mat']), 'prewin_frames', 'postwin_frames', 'resp_win', 'ind_stimAlone', 'ind_maskAlone', 'ind_plaid', 'ind_blank', 'trialsperstim');


%%

avg_resp_dir_rand = zeros(nCells,nStimDir,2);
for i = 1:nStimDir
    n = size(resp_cell{i,1,2},2);
    ind1 = randsample(1:n,ceil(n/2));
    ind2 = setdiff(1:n,ind1);
    avg_resp_dir_rand(:,i,1) = mean(resp_cell{i,1,2}(:,ind1),2);
    avg_resp_dir_rand(:,i,2) = mean(resp_cell{i,1,2}(:,ind2),2);
end


int = unique(diff(stimDirs));
component = avg_resp_dir(:,:,1,1,1)+circshift(avg_resp_dir(:,:,1,1,1),-120./int,2);
pattern = circshift(avg_resp_dir(:,:,1,1,1),-60./int,2);

%Calculate pattern and component prediction
comp_corr = zeros(nMaskPhas,nCells);
patt_corr = zeros(nMaskPhas,nCells);
comp_patt_corr = zeros(nMaskPhas,nCells);
plaid_corr = zeros(1,nCells);
plaid_corr_rand = zeros(1,nCells);

for iCell = 1:nCells
    for ip = 1:nMaskPhas
        comp_corr(ip,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,ip,2,1),component(iCell,:)));
        patt_corr(ip,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,ip,2,1),pattern(iCell,:)));
        comp_patt_corr(ip,iCell) = triu2vec(corrcoef(component(iCell,:),pattern(iCell,:)));
    end
end
Rp = ((patt_corr)-(comp_corr.*comp_patt_corr))./sqrt((1-comp_corr.^2).*(1-comp_patt_corr.^2));
Rc = ((comp_corr)-(patt_corr.*comp_patt_corr))./sqrt((1-patt_corr.^2).*(1-comp_patt_corr.^2));
Zp = (0.5.*log((1+Rp)./(1-Rp)))./sqrt(1./(nStimDir-3));
Zc = (0.5.*log((1+Rc)./(1-Rc)))./sqrt(1./(nStimDir-3));

ZcZp_diff = Zc-Zp;
ind1 = intersect(find(Zp(1,:)>1.28),find(Zp(1,:)-Zc(1,:)>1.28));


%% Trials per stimulus condition

figure;
subplot(2,1,1)
    bar(squeeze(trialsperstim(:,:,1))) 
    title('Trials per stim - GRATINGS')
    ylabel('# of trials')
    xlabel('direction')
subplot(2,1,2)
    bar(squeeze(trialsperstim(:,:,2))) 
    title('Trials per stim - PLAIDS')
    ylabel('# of trials')
    xlabel('direction')
    
   
print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_TrialsPerStim.pdf']),'-dpdf', '-fillpage')       

    

%% ZcZp of population (resp_ind)

figure; 
movegui('center')
scatter(Zc(resp_ind), Zp(resp_ind))
    hold on
    scatter(Zc(ind1),Zp(ind1));
    xlabel('Zc')
    ylabel('Zp')
    ylim([-4 8])
    xlim([-4 8])
    hold on
    if i==1; title('pattern cells at 0'); end
    plotZcZpBorders


sgtitle('Pattern direction selective cells')
print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_ZcZp.pdf']),'-dpdf', '-fillpage')       

%% Set responsive index

ind = intersect(intersect(resp_ind_dir,DSI_ind),p_dir); %resp to 1 grating, DSI>0.5, resp to one direction of gratings
% ind = ind1;

    
%% Direction tuning curves
respgrat = avg_resp_dir(:,:,1,1,1);
[prefresp, prefdir] = max(respgrat,[],2);
vecshift = zeros(size(respgrat,1),1)+6; 
vecshift = vecshift - prefdir;
respplaid = avg_resp_dir(:,:,:,2,1);

avg_resp_grat =[];
avg_resp_plaid = [];

for iCell = 1:size(respgrat,1)
    avg_resp_grat(iCell,:) = circshift(respgrat(iCell,:),vecshift(iCell),2);
    for i = 1:nMaskPhas
        avg_resp_plaid(iCell,:,i) = circshift(respplaid(iCell,:,i),vecshift(iCell)+(60./int),2);
    end
end


figure;
start = 1;
n = 1;
x=[-150:30:180];
for iCell = 1:nCells
    subplot(5,4,start)
        for im = 1:nMaskPhas
            plot(x, avg_resp_plaid(iCell,:,im))
            hold on
        end
        plot(x, avg_resp_grat(iCell,:),'k') 
        if iCell ==1; legend('0 deg', 'gratings'); end;
        xlabel('plaid direction')
        ylabel('df/f')
        idx = iCell==ind;
        if any(idx)
            subtitle(['cell ' num2str(iCell)],'fontweight','bold')
        else
            subtitle(['cell ' num2str(iCell)]); 
        end
    start = start+1;
    if start >20
        sgtitle([mouse ' ' date ' Direction tuning curves aligned to preferred grating direction'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_DirectionTuningGratingsAndPlaids_' num2str(n) '.pdf']),'-dpdf', '-fillpage')       
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle([mouse ' ' date ' Direction tuning curves aligned to preferred grating direction'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_DirectionTuningGratingsAndPlaids_' num2str(n) '.pdf']), '-dpdf','-fillpage')
    end
end
close all

%% Polar plots by individual cell

figure;
start = 1;
n = 1;
x=[-150:30:180];
x_rad = deg2rad(x);

for iCell =1:nCells
    subplot(5,4,start)
        for im = 1:nMaskPhas
            polarplot([x_rad x_rad(1)], [avg_resp_plaid(iCell,:,im) avg_resp_plaid(iCell,1,im)])
            hold on
        end
        polarplot([x_rad x_rad(1)], [avg_resp_grat(iCell,:) avg_resp_grat(iCell,1)],'k', 'LineWidth',2) 
        idx = iCell==ind;
        if any(idx)
            subtitle(['cell ' num2str(iCell)],'fontweight','bold')
        else
            subtitle(['cell ' num2str(iCell)])
        end
    start = start+1;    
    if start>20
        sgtitle([mouse ' ' date])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_PolarPlots_' num2str(n) '.pdf']), '-dpdf','-fillpage')
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle([mouse ' ' date])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_PolarPlots_' num2str(n) '.pdf']), '-dpdf','-fillpage')
    end
end     
close all

%% look at df/f across individual trials with comparison to statistical tests (ttest for grating, anova for plaid)

% c = linspace(0,round(max(centroid_dist)),round(max(centroid_dist)));


figure;
start=1;
n=1;
for iCell = 1:nCells
    subplot(5,4,start)
        for iDir = 1:nStimDir
%             cent = round(centroid_dist(trialInd{iDir,1,1}));
            scatter(iDir*ones(1,size(resp_cell{iDir,1,1},2)),resp_cell{iDir,1,1}(iCell,:));
            hold on
            max_df = max(resp_cell{iDir,1,1}(iCell,:),[],2); %find max df/f to use to print t-test result
            tt = h_resp(iCell,iDir,1,1);
            if tt == 1; text(iDir, max_df,'*'); end;
        end
        idx = iCell==ind;
        if any(idx)
            subtitle(['cell ' num2str(iCell) ', DSI ' num2str(DSI(iCell),'%.2f') ', p=' num2str(p_anova_dir(iCell),'%.3f')],'fontweight','bold')
        else
            subtitle(['cell ' num2str(iCell) ', DSI ' num2str(DSI(iCell),'%.2f') ', p=' num2str(p_anova_dir(iCell),'%.3f')])
        end
        ylabel('df/f'); xlabel('direction');
    start = start+1;    
    if start>20
        sgtitle([mouse ' ' date ' - Df/f response to gratings, t-test result (0/1), and anova'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_GratingResponses_' num2str(n) '.pdf']), '-dpdf','-fillpage')
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle([mouse ' ' date ' - Df/f response to gratings, t-test result (0/1), and anova'])
        print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_GratingResponses_' num2str(n) '.pdf']), '-dpdf','-fillpage')
    end     
end
    close all
   

%% Bootstrap one phase --> four phase

pd = makedist('Normal','mu',22.6042,'sigma',4.1192);
nboots = 1;

[boot_base, boot_amp, boot_rsq, boot_sse, boot_Zp, boot_Zc] = bootstrap_onephase(resp_cell, 4, pd, nboots);

figure;
for ib = 1:nboots
    cdfplot(boot_amp(:,ib));
    hold on
end
    xlabel('boot amp')
    
save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_bootstrapFits.mat']), 'boot_base', 'boot_amp', 'boot_rsq', 'boot_sse', 'boot_Zp', 'boot_Zc');




%% Pair wise euclidean distance 

pd = makedist('Normal','mu',22.6042,'sigma',4.1192);

data = resp_cell; 
nPhases = 4;

nStimDir = size(data,1);
nCells = size(data{1,1,1},1);

resps = double.empty(nCells,0);
avg_resp_new = zeros(nCells, nStimDir, nPhases);
avg_resp_grat = zeros(nCells, nStimDir);

for id = 1:nStimDir
    resps = data{id,1,2};
    for ip = 1:nPhases
        ntrials = round(random(pd));
        if ntrials > size(resps,2)
            ntrials = size(resps,2);
        end
        resp_cell_new{id,ip} = resps(:,randperm(size(resps,2),ntrials));
        avg_resp_new(:,id,ip) =  mean(resp_cell_new{id,ip}(:,:),2); 
    end
    avg_resp_grat(:,id) =  mean(data{id,1,1}(:,:),2);
end

component = avg_resp_grat(:,:)+circshift(avg_resp_grat(:,:),-120./nStimDir,2);
pattern = circshift(avg_resp_grat(:,:),-60./nStimDir,2);

%Calculate pattern and component prediction
comp_corr = zeros(nPhases,nCells);
patt_corr = zeros(nPhases,nCells);
comp_patt_corr = zeros(nPhases,nCells);
plaid_corr = zeros(1,nCells);
plaid_corr1 = zeros(1,nCells);
plaid_corr2 = zeros(1,nCells);
plaid_corr3 = zeros(1,nCells);
plaid_corr4 = zeros(1,nCells);
plaid_corr5 = zeros(1,nCells);
plaid_corr6 = zeros(1,nCells);
plaid_corr_rand = zeros(1,nCells);

for iCell = 1:nCells
    for ip = 1:nPhases
        comp_corr(ip,iCell) = triu2vec(corrcoef(avg_resp_new(iCell,:,ip),component(iCell,:)));
        patt_corr(ip,iCell) = triu2vec(corrcoef(avg_resp_new(iCell,:,ip),pattern(iCell,:)));
        comp_patt_corr(ip,iCell) = triu2vec(corrcoef(component(iCell,:),pattern(iCell,:)));
    end
    plaid_corr1(1,iCell) = triu2vec(corrcoef(avg_resp_new(iCell,:,1),avg_resp_new(iCell,:,2)));
    plaid_corr2(1,iCell) = triu2vec(corrcoef(avg_resp_new(iCell,:,1),avg_resp_new(iCell,:,3)));
    plaid_corr3(1,iCell) = triu2vec(corrcoef(avg_resp_new(iCell,:,1),avg_resp_new(iCell,:,4)));
    plaid_corr4(1,iCell) = triu2vec(corrcoef(avg_resp_new(iCell,:,2),avg_resp_new(iCell,:,3)));
    plaid_corr5(1,iCell) = triu2vec(corrcoef(avg_resp_new(iCell,:,2),avg_resp_new(iCell,:,4)));
    plaid_corr6(1,iCell) = triu2vec(corrcoef(avg_resp_new(iCell,:,3),avg_resp_new(iCell,:,4)));
    plaid_corr(1,iCell) = (plaid_corr1(1,iCell)+plaid_corr2(1,iCell)+plaid_corr3(1,iCell)+plaid_corr4(1,iCell)+plaid_corr5(1,iCell)+plaid_corr6(1,iCell))/6;
end
Rp = ((patt_corr)-(comp_corr.*comp_patt_corr))./sqrt((1-comp_corr.^2).*(1-comp_patt_corr.^2));
Rc = ((comp_corr)-(patt_corr.*comp_patt_corr))./sqrt((1-patt_corr.^2).*(1-comp_patt_corr.^2));
Zp = (0.5.*log((1+Rp)./(1-Rp)))./sqrt(1./(nStimDir-3));
Zc = (0.5.*log((1+Rc)./(1-Rc)))./sqrt(1./(nStimDir-3));


ZpZc = [];
for iCell = 1:nCells
    ZpZc(:,1,iCell) = Zp(:,iCell);
    ZpZc(:,2,iCell) = Zc(:,iCell);
end

ZpZcPWdist = double.empty(6,0);
for iCell = 1:nCells
    ZpZcPWdist(:,iCell) = pdist(ZpZc(:,:,iCell));
end

save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_ZpZc_pairwiseDist.mat']), 'ind', 'Zp', 'Zc', 'ZpZcPWdist', 'plaid_corr');


stop
%% Bootstrap trials

num_trials = [11 22 33 44];
nboots = 500;

[std_avg] = bootstrap_trials(resp_cell,num_trials,nboots);

figure;
    plot(num_trials,std_avg);
    hold on
    xlabel('nTrials')
    xticks([11 22 33 44])
    xlim([0 50])
    ylim([0 2])
    ylabel('std avg')
    title(['nboots = ' num2str(nboots)])

print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_bootstrapTrials.pdf']), '-dpdf','-fillpage')
save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_bootstrapTrials.mat']), 'std_avg', 'nboots', 'num_trials');

 