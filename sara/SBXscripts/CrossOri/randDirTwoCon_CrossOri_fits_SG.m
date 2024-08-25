clc; clear all; close all;
doRedChannel = 1;
ds = 'CrossOriRandDirTwoCon_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = length(expt);
frame_rate = 15;
seed = rng;

max_dist = 5;

iexp = [3];

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

if ~exist('sz','var')
    nCells = size(data_tc,2);
    nFrames = size(data_tc,1);
else
end


SF_all = celleqel2mat_padded(input.tStimOneGratingSpatialFreqCPD);
SFs = unique(SF_all);
nSF = length(SFs);
if nSF>1
    doSF = 1;
    nF = nSF;
else 
    doSF = 0;
    nF = 1;
end
TF_all = celleqel2mat_padded(input.tMaskOneGratingTemporalFreqCPS);
TFs = unique(TF_all);
nTF = length(TFs);
if nTF>1
    doTF = 1;
    nF = nTF;
else 
    doTF = 0;
    nF = 1;
end

prewin_frames = frame_rate;
nFramesOn = unique(celleqel2mat_padded(input.nStimOneFramesOn));
postwin_frames = unique(celleqel2mat_padded(input.nStimOneFramesOn));
tt = (1-prewin_frames:postwin_frames).*(1/frame_rate);
data_resp = nan(prewin_frames+postwin_frames,nCells,nTrials);
data_f = nan(1,nCells,nTrials);

if exist('npSub_tc','var') == 0
    npSub_tc = data_tc;
    clear data_tc
end

if exist('sz', 'var')
    for itrial = 1:nTrials
        if cStimOn(itrial) + postwin_frames < sz(3)  %sz(3) should be = to nFrames
            data_resp(:,:,itrial) = npSub_tc(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
        end
    end
else
    for itrial = 1:nTrials
        if cStimOn(itrial) + postwin_frames < nFrames
            data_resp(:,:,itrial) = npSub_tc(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
        end
    end
end


data_f = mean(data_resp(1:prewin_frames,:,:),1);
data_dfof_tc = (data_resp-data_f)./data_f;

maskPhas_all = celleqel2mat_padded(input.tMaskOneGratingPhaseDeg);
maskPhas = unique(maskPhas_all);
nMaskPhas = length(maskPhas);

nMaskCon = length(maskCons);
nStimCon = length(stimCons);

ind_stimAlone = intersect(intersect(find(stimCon_all),find(maskCon_all==0)),find(centroid_dist<max_dist));
ind_maskAlone = intersect(intersect(find(stimCon_all==0),find(maskCon_all)),find(centroid_dist<max_dist));
ind_plaid = intersect(intersect(find(stimCon_all>0),find(maskCon_all>0)),find(centroid_dist<max_dist));
ind_blank = intersect(intersect(find(stimCon_all==0),find(maskCon_all==0)),find(centroid_dist<max_dist));

base_cell = cell(nStimDir,nMaskCon,nStimCon);
resp_cell = cell(nStimDir,nMaskCon,nStimCon);
trialInd = cell(nStimDir,nMaskCon,nStimCon);
trialsperstim = zeros(nStimDir,nMaskCon,nStimCon);
h_resp =zeros(nCells,nStimDir,nMaskCon,nStimCon);
p_resp =zeros(nCells,nStimDir,nMaskCon,nStimCon);
base_win = prewin_frames-ceil(frame_rate/2):prewin_frames;
resp_win = prewin_frames+5:prewin_frames+nFramesOn;
resp_blank = squeeze(mean(data_dfof_tc(prewin_frames+frame_rate:prewin_frames+nFramesOn,:,ind_blank),1));
data_dfof_con_ph_tc_avg = nan(prewin_frames+postwin_frames, nCells, nStimDir, nMaskCon, nStimCon);

avg_resp_dir = zeros(nCells,nStimDir,nMaskCon,nStimCon);
all_resp_dir = [];
all_resp_plaid = cell(1,nMaskCon,nStimCon);
all_dir = [];
all_plaid = cell(1,nMaskCon,nStimCon);
nStim = nStimDir;


for iDir = 1:nStimDir
    ind_stimdir = find(stimDir_all == stimDirs(iDir));
    ind_maskdir = find(maskDir_all == stimDirs(iDir));
    ind_diralone = [intersect(ind_stimdir, ind_stimAlone), intersect(ind_maskdir, ind_maskAlone)];
    ind_dirplaid = [intersect(ind_stimdir, ind_plaid)];
    all_resp_dir = [all_resp_dir squeeze(mean(data_dfof_tc(resp_win,:,ind_diralone),1))];
    all_dir = [all_dir iDir.*ones(size(ind_diralone))];
    for im = 1:nMaskCon
        ind_mask = find(maskCon_all == maskCons(im));
        for it = 1:nStimCon
            if maskCons(im) == stimCons(it)
                ind_stim = find(stimCon_all == stimCons(it));
                ind = intersect(ind_dirplaid, intersect(ind_stim,ind_mask));
                trialsperstim(iDir,im,it) = length(ind);
                resp_cell{iDir,im,it} = squeeze(mean(data_dfof_tc(resp_win,:,ind),1));
                base_cell{iDir,im,it} = squeeze(mean(data_dfof_tc(base_win,:,ind),1));
                data_dfof_con_ph_tc_avg(:,:,iDir,im,it,1) = squeeze(nanmean(data_dfof_tc(:,:,ind),3));
                data_dfof_con_ph_tc_avg(:,:,iDir,im,it,2) = squeeze(nanstd(data_dfof_tc(:,:,ind),[],3)./sqrt(length(ind)));
                [h_resp(:,iDir,im,it), p_resp(:,iDir,im,it)] = ttest2(resp_cell{iDir,im,it},base_cell{iDir,im,it},'dim',2,'tail','right','alpha', 0.05./nStim);
                avg_resp_dir(:,iDir,im,it,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind),1),3));
                avg_resp_dir(:,iDir,im,it,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,ind),1),[],3)./sqrt(length(ind)));
                trialInd{iDir,im,it} = ind;
            elseif maskCons(im)==0;
                ind_stim = find(stimCon_all == stimCons(it));
                ind = intersect(ind_diralone, intersect(ind_stim,ind_mask));
                trialsperstim(iDir,im,it) = length(ind);
                resp_cell{iDir,im,it} = squeeze(mean(data_dfof_tc(resp_win,:,ind),1));
                base_cell{iDir,im,it} = squeeze(mean(data_dfof_tc(base_win,:,ind),1));
                data_dfof_con_ph_tc_avg(:,:,iDir,im,it,1) = squeeze(nanmean(data_dfof_tc(:,:,ind),3));
                data_dfof_con_ph_tc_avg(:,:,iDir,im,it,2) = squeeze(nanstd(data_dfof_tc(:,:,ind),[],3)./sqrt(length(ind)));
                avg_resp_dir(:,iDir,im,it,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind),1),3));
                avg_resp_dir(:,iDir,im,it,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,ind),1),[],3)./sqrt(length(ind)));
                trialInd{iDir,im,it} = ind;
            end
            if iDir == 1
                all_resp_plaid{1,im,it} = [];
                all_plaid{1,im,it} = [];
            end
            all_resp_plaid{1,im,it} = [all_resp_plaid{1,im,it} squeeze(mean(data_dfof_tc(resp_win,:,ind),1))];
            all_plaid{1,im,it} = [all_plaid{1,im,it} iDir.*ones(size(ind))];
        end
    end
end
 
all_resp_plaid_re = reshape(all_resp_plaid,1,[]);
all_plaid_re = reshape(all_plaid,1,[]);
nCond = size(all_resp_plaid_re,2);

respgratlow_ind = find(sum(h_resp(:,:,1,1),2));
respgrathigh_ind = find(sum(h_resp(:,:,1,2),2));
respplaidlow_ind = find(sum(h_resp(:,:,2,1),2));
respplaidhigh_ind = find(sum(h_resp(:,:,3,2),2));

resp_ind = unique([respgratlow_ind; respgrathigh_ind; respplaidlow_ind; respplaidhigh_ind]);

p_anova_dir = zeros(1,nCells);
p_anova_plaid = zeros(2,nCells);
for iCell = 1:nCells
    p_anova_dir(iCell) = anova1(all_resp_dir(iCell,:), all_dir, 'off'); %direction selective to gratings
    for in = 1:nCond
        p_anova_plaid(in,iCell) = anova1(all_resp_plaid_re{in}(iCell,:), all_plaid_re{in}, 'off'); %direction selective to plaids
    end
end

p_dir = find(p_anova_dir<0.05);
p_plaid1 = find(p_anova_plaid(1,:)<0.05);
p_plaid2 = find(p_anova_plaid(2,:)<0.05);
p_plaid3 = find(p_anova_plaid(3,:)<0.05);
p_plaid4 = find(p_anova_plaid(4,:)<0.05);
p_plaid5 = find(p_anova_plaid(5,:)<0.05);
p_plaid6 = find(p_anova_plaid(6,:)<0.05);
p_all = unique([p_dir,p_plaid1,p_plaid2,p_plaid3,p_plaid4,p_plaid5,p_plaid6]); %significantly responsive to a direction (anova) for gratings or any plaid set


for iCell = 1:nCells
    [max_valL max_indL] = max(avg_resp_dir(iCell,:,1,1,1));
    null_indL = max_indL+(nStimDir./2);
    null_indL(find(null_indL>nStimDir)) = null_indL(find(null_indL>nStimDir))-nStimDir;
    min_valL = avg_resp_dir(iCell,null_indL,1,1,1);
    if min_valL < 0; min_valL = 0; end
    DSIL(iCell) = (max_valL-min_valL)./(max_valL+min_valL);
    DSI_maxIndL(iCell) = max_indL; 

    [max_valH max_indH] = max(avg_resp_dir(iCell,:,1,2,1));
    null_indH = max_indH+(nStimDir./2);
    null_indH(find(null_indH>nStimDir)) = null_indH(find(null_indH>nStimDir))-nStimDir;
    min_valH = avg_resp_dir(iCell,null_indH,1,1,1);
    if min_valH < 0; min_valH = 0; end
    DSIH(iCell) = (max_valH-min_valH)./(max_valH+min_valH);
    DSI_maxIndH(iCell) = max_indH; 
end

if ~exist(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' max_dist]),'dir')
    mkdir(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' num2str(max_dist)]))
end
save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_respData.mat']), 'resp_cell', 'data_dfof_tc', 'p_all', 'resp_ind', 'respgratlow_ind', 'respgrathigh_ind', 'respplaidlow_ind', 'respplaidhigh_ind', 'frame_rate', 'h_resp', 'avg_resp_dir');
save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_stimData.mat']), 'prewin_frames', 'postwin_frames', 'resp_win', 'ind_diralone', 'ind_plaid','trialsperstim', 'DSIL', 'DSIH', 'DSI_maxIndL', 'DSI_maxIndH');


%%

int = unique(diff(stimDirs));
componentL = avg_resp_dir(:,:,1,1,1)+circshift(avg_resp_dir(:,:,1,1,1),-120./int,2);
patternL = circshift(avg_resp_dir(:,:,1,1,1),-60./int,2);
componentH = avg_resp_dir(:,:,1,2,1)+circshift(avg_resp_dir(:,:,1,1,1),-120./int,2);
patternH = circshift(avg_resp_dir(:,:,1,2,1),-60./int,2);


%Calculate pattern and component prediction
comp_corr = zeros(nStimCon,nCells);
patt_corr = zeros(nStimCon,nCells);
comp_patt_corr = zeros(nStimCon,nCells);
plaid_corr = zeros(1,nCells);


for iCell = 1:nCells
    for it = 1:nStimCon
        if it == 1
            comp_corr(it,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,it+1,it,1),componentL(iCell,:)));
            patt_corr(it,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,it+1,it,1),patternL(iCell,:)));
            comp_patt_corr(it,iCell) = triu2vec(corrcoef(componentL(iCell,:),patternL(iCell,:)));
        else
            comp_corr(it,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,it+1,it,1),componentH(iCell,:)));
            patt_corr(it,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,it+1,it,1),patternH(iCell,:)));
            comp_patt_corr(it,iCell) = triu2vec(corrcoef(componentH(iCell,:),patternH(iCell,:)));
        end
    end
    plaid_corr(1,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,2,1,1),avg_resp_dir(iCell,:,3,2,1)));
end


Rp = ((patt_corr)-(comp_corr.*comp_patt_corr))./sqrt((1-comp_corr.^2).*(1-comp_patt_corr.^2));
Rc = ((comp_corr)-(patt_corr.*comp_patt_corr))./sqrt((1-patt_corr.^2).*(1-comp_patt_corr.^2));
Zp = (0.5.*log((1+Rp)./(1-Rp)))./sqrt(1./(nStimDir-3));
Zc = (0.5.*log((1+Rc)./(1-Rc)))./sqrt(1./(nStimDir-3));

ZcZp_diff = Zc-Zp;
indL = intersect(find(Zp(1,:)>1.28),find(Zp(1,:)-Zc(1,:)>1.28));
indH = intersect(find(Zp(2,:)>1.28),find(Zp(2,:)-Zc(2,:)>1.28));
ind_p = {indL,indH};
ind_arr = [indL,indH];
ind_pds = unique(ind_arr);

save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_plaidCorr.mat']), 'plaid_corr');

%%

DSI_ind = find(DSIL>0.5);
ind = intersect(intersect(DSI_ind,p_all),resp_ind);

%% Trials per stim

figure;
subplot(4,1,1)
    bar(squeeze(trialsperstim(:,1,1))) 
    title('Trials per stim - 0.25 GRATINGS')
    ylabel('# of trials'); ylim([0 50])
    xlabel('direction')
subplot(4,1,2)
    bar(squeeze(trialsperstim(:,1,2)))
    title('Trials per stim - 0.5 GRATINGS')
    ylabel('# of trials'); ylim([0 50])
    xlabel('direction')
subplot(4,1,3)
    bar(squeeze(trialsperstim(:,2,1))) 
    title('Trials per stim - 0.5 PLAIDS')
    ylabel('# of trials'); ylim([0 50])
    xlabel('direction')
subplot(4,1,4)
    bar(squeeze(trialsperstim(:,3,2))) 
    title('Trials per stim - 1.0 PLAIDS')
    ylabel('# of trials'); ylim([0 50])
    xlabel('direction')

print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_TrialsPerStim.pdf']),'-dpdf', '-fillpage')       


%%

for i = 1:2
    ind_p{i} = intersect(ind_p{i},ind);
end

figure; 
movegui('center')
for i = 1:2
    subplot(2,2,i)
    scatter(Zc(i,ind), Zp(i,ind))
    hold on
    scatter(Zc(i,ind_p{1}),Zp(i,ind_p{1}));
    xlabel('Zc')
    ylabel('Zp')
    ylim([-4 8])
    xlim([-4 8])
    hold on
    if i==1; title('pattern cells at 0.5 plaid'); end
    plotZcZpBorders
end
for i = 1:2
    subplot(2,2,i+2)
    scatter(Zc(i,ind), Zp(i,ind))
    hold on
    scatter(Zc(i,ind_p{2}),Zp(i,ind_p{2}));
    xlabel('Zc')
    ylabel('Zp')
    ylim([-4 8])
    xlim([-4 8])
    hold on
    if i==1; title('pattern cells at 1.0 plaid'); end
    plotZcZpBorders
end

sgtitle('Pattern direction selective cells at 2 contrasts')
print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_ZcZp.pdf']),'-dpdf', '-fillpage')       

%%
respgrat = avg_resp_dir(:,:,1,1,1);
[prefresp, prefdir] = max(respgrat,[],2);
vecshift = zeros(size(respgrat,1),1)+6; 
vecshift = vecshift - prefdir;
respplaid = avg_resp_dir(:,:,2,1,1);
respplaid(:,:,2) = avg_resp_dir(:,:,3,2,1);

avg_resp_grat =[];
avg_resp_plaid = [];

for iCell = 1:size(respgrat,1)
    avg_resp_grat(iCell,:) = circshift(respgrat(iCell,:),vecshift(iCell),2);
    for it = 1:nStimCon
        avg_resp_plaid(iCell,:,it) = circshift(respplaid(iCell,:,it),vecshift(iCell)+(60./int),2);
    end
end


figure;
start = 1;
n = 1;

x=[-150:30:180];
x_rad = deg2rad(x);
for iCell =1:nCells
    subplot(5,4,start)
        for it = 1:nStimCon
            if it == 1
                polarplot([x_rad x_rad(1)], [avg_resp_plaid(iCell,:,it) avg_resp_plaid(iCell,1,it)],'Color',[0 0.4470 0.7410],'LineWidth',1)
                hold on
                polarplot([x_rad x_rad(1)], [avg_resp_grat(iCell,:,it) avg_resp_grat(iCell,1,it)],'Color',[0 0.4470 0.7410]) 
            else
                polarplot([x_rad x_rad(1)], [avg_resp_plaid(iCell,:,it) avg_resp_plaid(iCell,1,it)],'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
                hold on
            end
        end
        idx = iCell==ind;
        if any(idx)
            subtitle(['cell ' num2str(iCell)],'fontweight','bold')
        else
            subtitle(['cell ' num2str(iCell)]); 
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

%% Population averaging (comparing across plaids)

avg_resp_plaid_all = squeeze(mean(avg_resp_plaid(ind,:,:),1));

PCI = Zp-Zc;


figure;
subplot(3,2,1)
    polarplot([x_rad x_rad(1)], [avg_resp_plaid_all(:,1)' avg_resp_plaid_all(1,1)],'Color',[0 0.4470 0.7410],'LineWidth',2)
    hold on
    polarplot([x_rad x_rad(1)], [avg_resp_plaid_all(:,2)' avg_resp_plaid_all(1,2)],'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
    subtitle('Population average')
subplot(3,2,2)
    cdfplot(Zp(1,:))
    hold on
    cdfplot(Zp(2,:))
    xlabel('Zp'); xlim([-5 5])
    ylabel('Fraction of cells')
    set(gca,'TickDir','out'); box off; axis square; grid off
subplot(3,2,3)
    cdfplot(Zc(1,:))
    hold on
    cdfplot(Zc(2,:))
    xlabel('Zc'); xlim([-5 5])
    ylabel('Fraction of cells')   
    set(gca,'TickDir','out'); box off; axis square; grid off
 subplot(3,2,4)
    cdfplot(PCI(1,:))
    hold on
    cdfplot(PCI(2,:))
    xlabel('PCI (Zp - Zc)'); xlim([-5 5])
    ylabel('Fraction of cells')   
    set(gca,'TickDir','out'); box off; axis square; grid off
  subplot(3,2,5)
    cdfplot(DSIL)
    hold on
    cdfplot(DSIH)
    xlabel('DSI'); 
    ylabel('Fraction of cells')   
    set(gca,'TickDir','out'); box off; axis square; grid off

print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_PopulationAverage.pdf']), '-dpdf','-fillpage')


