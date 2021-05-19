clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandDirRandPhase_ExptList';
rc = behavConstsAV;
eval(ds)
nexp = length(expt);

iexp = 32

frame_rate = 15;

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc;
ImgFolder = expt(iexp).coFolder;
time = expt(iexp).coTime;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
%LG_base = '\\CRASH.dhe.duke.edu\data\home\lindsey';

fprintf([mouse ' ' date '\n'])

%% Pref direction analysis
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))

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

ind_stimAlone = intersect(find(stimCon_all),find(maskCon_all==0));
ind_maskAlone = intersect(find(stimCon_all==0),find(maskCon_all));
ind_plaid = intersect(find(stimCon_all),find(maskCon_all));
ind_blank = intersect(find(stimCon_all==0),find(maskCon_all==0));

resp_cell = cell(nStimDir,2,nSF);
trialsperstim = zeros(nStimDir,2,nSF);
h_resp =zeros(nCells,nStimDir,2,nSF);
p_resp =zeros(nCells,nStimDir,2,nSF);
resp_win = prewin_frames+5:prewin_frames+nFramesOn;
resp_blank = squeeze(mean(data_dfof_tc(prewin_frames+frame_rate:prewin_frames+nFramesOn,:,ind_blank),1));
avg_resp_dir = zeros(nCells, nStimDir, 2, nSF, 2);

all_resp_dir = cell(1,nSF);
all_resp_plaid = cell(1,nSF);
all_dir = cell(1,nSF);
all_plaid = cell(1,nSF);

p_anova_dir = zeros(nSF,nCells);
p_anova_plaid = zeros(nSF,nCells);

resp_ind{iSF} = cell(1,nSF);
resp_ind_dir{iSF} = cell(1,nSF);
resp_ind_plaid{iSF} = cell(1,nSF);

nStim = nStimDir;
for iSF = 1:nSF
    ind_SF = find(SF_all == SFs(iSF));
    for iDir = 1:nStimDir
        ind_stimdir = find(stimDir_all == stimDirs(iDir));
        ind_maskdir = find(maskDir_all == maskDirs(iDir));
        ind_diralone = intersect(ind_SF,[intersect(ind_stimdir, ind_stimAlone), intersect(ind_maskdir, ind_maskAlone)]);
        ind_dirplaid = intersect(ind_SF,[intersect(ind_stimdir, ind_plaid)]);
        trialsperstim(iDir,1,iSF) = length(ind_diralone);
        trialsperstim(iDir,2,iSF) = length(ind_dirplaid);
        resp_cell{iDir,1,iSF} = squeeze(mean(data_dfof_tc(resp_win,:,ind_diralone),1));
        resp_cell{iDir,2,iSF} = squeeze(mean(data_dfof_tc(resp_win,:,ind_dirplaid),1));
        [h_resp(:,iDir,1,iSF), p_resp(:,iDir,1,iSF)] = ttest2(resp_cell{iDir,1},resp_blank,'dim',2,'tail','right','alpha', 0.05./nStim);
        [h_resp(:,iDir,2,iSF), p_resp(:,iDir,2,iSF)] = ttest2(resp_cell{iDir,2},resp_blank,'dim',2,'tail','right','alpha', 0.05./nStim);
        avg_resp_dir(:,iDir,1,iSF,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind_diralone),1),3));
        avg_resp_dir(:,iDir,1,iSF,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,ind_diralone),1),[],3)./sqrt(length(ind_diralone)));
        avg_resp_dir(:,iDir,2,iSF,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind_dirplaid),1),3));
        avg_resp_dir(:,iDir,2,iSF,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,ind_dirplaid),1),[],3)./sqrt(length(ind_dirplaid)));
        all_resp_dir{iSF} = [all_resp_dir{iSF} squeeze(mean(data_dfof_tc(resp_win,:,ind_diralone),1))];
        all_resp_plaid{iSF} = [all_resp_plaid{iSF} squeeze(mean(data_dfof_tc(resp_win,:,ind_dirplaid),1))];
        all_dir{iSF} = [all_dir{iSF} iDir.*ones(size(ind_diralone))];
        all_plaid{iSF} = [all_plaid{iSF} iDir.*ones(size(ind_dirplaid))];
    end

    resp_ind{iSF} = find(sum(sum(h_resp(:,:,:,iSF),2),3));
    resp_ind_dir{iSF} = find(sum(h_resp(:,:,1,iSF),2));
    resp_ind_plaid{iSF} = find(sum(h_resp(:,:,2,iSF),2));
    
    for iCell = 1:nCells
        p_anova_dir(iSF,iCell) = anova1(all_resp_dir{iSF}(iCell,:), all_dir{iSF}, 'off');
        p_anova_plaid(iSF,iCell) = anova1(all_resp_plaid{iSF}(iCell,:), all_plaid{iSF}, 'off');
    end
end
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']), 'resp_cell', 'data_dfof_tc', 'resp_ind', 'frame_rate', 'h_resp', 'avg_resp_dir');
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']), 'prewin_frames', 'postwin_frames', 'resp_win', 'ind_stimAlone', 'ind_maskAlone', 'ind_plaid', 'ind_blank');


%%
b_ori = cell(1,nSF);
k1_ori = cell(1,nSF);
R1_ori = cell(1,nSF);
u1_ori = cell(1,nSF);
sse_ori = cell(1,nSF);
R_square_ori = cell(1,nSF);
y_ori_fit = cell(1,nSF);

comp_corr = zeros(iSF,nCells);
patt_corr = zeros(iSF,nCells);
comp_patt_corr = zeros(iSF,nCells);

Rc = zeros(nSF,nCells);
Rp = zeros(nSF,nCells);
Zc = zeros(nSF,nCells);
Zp = zeros(nSF,nCells);

stim_OSI = zeros(iSF,nCells);

for iSF = 1:nSF
    b_ori{iSF} = nan(1,nCells);
    k1_ori{iSF} = nan(1,nCells);
    R1_ori{iSF} = nan(1,nCells);
    u1_ori{iSF} = nan(1,nCells);
    sse_ori{iSF} = nan(1,nCells);
    R_square_ori{iSF} = nan(1,nCells);
    range = 1:1:180;
    y_ori_fit{iSF} = nan(nCells,length(range));

    figure; 
    movegui('center')
    start = 1;
    n = 1;
    for iC = 1:length(resp_ind{iSF})
        iCell = resp_ind{iSF}(iC);
        if start>25
            suptitle([date ' ' mouse ' Ori Tuning; SF = ' num2str(SFs(iSF))])
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_plaidResp_SF' num2str(SFs(iSF)) '_oriTuning_' num2str(n) '.pdf']),'-dpdf', '-fillpage')
            n = n+1;
            figure; 
            movegui('center')
            start = 1;
        end
        subplot(5,5,start) 
        errorbar(stimDirs, avg_resp_dir(iCell,:,1,iSF,1), avg_resp_dir(iCell,:,1,iSF,2))
        hold on
        data = [avg_resp_dir(iCell,:,1,iSF,1) avg_resp_dir(iCell,1,1,iSF,1)];
        theta = [deg2rad(stimDirs) pi];
        [b_ori{iSF}(:,iCell),k1_ori{iSF}(:,iCell),R1_ori{iSF}(:,iCell),u1_ori{iSF}(:,iCell),sse_ori{iSF}(:,iCell),R_square_ori{iSF}(:,iCell)] ...
            = miaovonmisesfit_ori(theta,data);
        y_ori_fit{iSF}(iCell,:) = b_ori{iSF}(:,iCell)+R1_ori{iSF}(:,iCell).*exp(k1_ori{iSF}(:,iCell).*(cos(2*(deg2rad(range)-u1_ori{iSF}(:,iCell)))-1));
        plot(range,y_ori_fit{iSF}(iCell,:))
        errorbar(stimDirs, avg_resp_dir(iCell,:,2,iSF,1), avg_resp_dir(iCell,:,2,iSF,2))
        start = start+1;
        if p_anova_plaid(iSF,iCell)<0.05
            title('*')
        end
    end
    suptitle([date ' ' mouse ' Ori Tuning; SF = ' num2str(SFs(iSF))])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_plaidResp_SF' num2str(SFs(iSF)) '_oriTuning_' num2str(n) '.pdf']),'-dpdf', '-fillpage')       

    int = unique(diff(stimDirs));
    avg_resp_alldir = cat(2,avg_resp_dir,avg_resp_dir);
    component = avg_resp_alldir(:,:,1,iSF,1)+circshift(avg_resp_alldir(:,:,1,iSF,1),-90./int,2);
    pattern = circshift(avg_resp_alldir(:,:,1,iSF,1),-45./int,2);

    for iCell = 1:nCells
        comp_corr(iSF,iCell) = triu2vec(corrcoef(avg_resp_alldir(iCell,:,2,iSF,1),component(iCell,:)));
        patt_corr(iSF,iCell) = triu2vec(corrcoef(avg_resp_alldir(iCell,:,2,iSF,1),pattern(iCell,:)));
        comp_patt_corr(iSF,iCell) = triu2vec(corrcoef(component(iCell,:),pattern(iCell,:)));
    end
    Rp(iSF,:) = ((patt_corr(iSF,:))-(comp_corr(iSF,:).*comp_patt_corr(iSF,:)))./sqrt((1-comp_corr(iSF,:).^2).*(1-comp_patt_corr(iSF,:).^2));
    Rc(iSF,:) = ((comp_corr(iSF,:))-(patt_corr(iSF,:).*comp_patt_corr(iSF,:)))./sqrt((1-patt_corr(iSF,:).^2).*(1-comp_patt_corr(iSF,:).^2));
    Zp(iSF,:) = (0.5.*log((1+Rp(iSF,:))./(1-Rp(iSF,:))))./sqrt(1./(nStimDir-3));
    Zc(iSF,:) = (0.5.*log((1+Rc(iSF,:))./(1-Rc(iSF,:))))./sqrt(1./(nStimDir-3));

    figure; 
    movegui('center')
    scatter(Zc(iSF,resp_ind{iSF}), Zp(iSF,resp_ind{iSF}))
    xlabel('Zc')
    ylabel('Zp')
    ylim([-4 8])
    xlim([-4 8])
    hold on
    plotZcZpBorders
    suptitle([date ' ' mouse '; SF = ' num2str(SFs(iSF))])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_ZcZp_SF' num2str(SFs(iSF)) '.pdf']),'-dpdf', '-fillpage')       

    ind = find(Zp(iSF,:)>1);
    figure; 
    movegui('center')
    start = 1;
    [n n2] = subplotn(length(ind));
    for iC = 1:length(ind)
        iCell = ind(iC);
        subplot(n,n2,start) 
        errorbar(stimDirs, avg_resp_dir(iCell,:,1,iSF,1), avg_resp_dir(iCell,:,1,iSF,2))
        hold on
        errorbar(stimDirs, avg_resp_dir(iCell,:,2,iSF,1), avg_resp_dir(iCell,:,2,iSF,2))
        start = start+1;
        title(['Zc = ' num2str(chop(Zc(iSF,iCell),2)) '; Zp= ' num2str(chop(Zp(iSF,iCell),2))])
    end
    suptitle([date ' ' mouse ' Cells with Zp>1; SF = ' num2str(SFs(iSF))])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuning_Zpabove1_SF' num2str(SFs(iSF)) '.pdf']),'-dpdf', '-fillpage')       


    ind = find(Zc(iSF,:)>2);
    figure; 
    movegui('center')
    start = 1;
    [n n2] = subplotn(length(ind));
    for iC = 1:length(ind)
        iCell = ind(iC);
        subplot(n,n2,start) 
        errorbar(stimDirs, avg_resp_dir(iCell,:,1,iSF,1), avg_resp_dir(iCell,:,1,iSF,2))
        hold on
        errorbar(stimDirs, avg_resp_dir(iCell,:,2,iSF,1), avg_resp_dir(iCell,:,2,iSF,2))
        start = start+1;
        title(['Zc = ' num2str(chop(Zc(iSF,iCell),2)) '; Zp= ' num2str(chop(Zp(iSF,iCell),2))])
    end
    suptitle([date ' ' mouse ' Cells with Zc>2; SF = ' num2str(SFs(iSF))])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuning_Zcabove2_SF' num2str(SFs(iSF)) '.pdf']),'-dpdf', '-fillpage')       

    avg_resp_dir_rect = avg_resp_dir;
    avg_resp_dir_rect(find(avg_resp_dir<0)) = 0;
    nStimOri = size(avg_resp_dir,2);
    for iCell = 1:nCells
        [max_val max_ind] = max(avg_resp_dir_rect(iCell,:,1,iSF,1));
        orth_ind = max_ind+(nStimDir./2);
        orth_ind(find(orth_ind>nStimDir)) = orth_ind(find(orth_ind>nStimDir))-nStimDir;
        min_val = avg_resp_dir_rect(iCell,orth_ind,1,iSF,1);
        stim_OSI(iSF,iCell) = (max_val-min_val)./(max_val+min_val);
    end
end

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriAnalysis.mat']), 'component','pattern','Rp', 'Rc', 'Zp', 'Zc', 'stim_OSI', 'nCells', 'b_ori', 'k1_ori', 'R1_ori', 'u1_ori', 'sse_ori', 'R_square_ori', 'y_ori_fit');

%%

figure; 
for iSF = 1:nSF
    cdfplot(k1_ori{iSF})
    hold on
end

figure; 
for iSF = 1:nSF
    cdfplot(Zp(iSF,:))
    hold on
end
movegui('center')

figure; 
for iSF = 1:nSF
    cdfplot(Zc(iSF,:))
    hold on
end
movegui('center')

figure; 
for iSF = 1:nSF
    cdfplot(stim_OSI(iSF,:))
    hold on
end
movegui('center')

figure;
plot(stim_OSI)
movegui('center')


