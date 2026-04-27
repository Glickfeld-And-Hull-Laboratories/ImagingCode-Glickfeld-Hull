clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandDirType2_ExptList';
eval(ds)
nexp = length(expt);

iexp = 2;

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
%LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\David';

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

tDir = celleqel2mat_padded(input.tStimOneGratingDirectionDeg);
tTF = celleqel2mat_padded(input.tStimOneGratingTemporalFreqCPS);
tMTF = celleqel2mat_padded(input.tMaskOneGratingTemporalFreqCPS);
tCon = celleqel2mat_padded(input.tMaskOneGratingContrast);

Dirs = unique(tDir);
nDir = length(Dirs);
TFs = unique(tTF);
nTF = length(TFs);
mTFs = unique(tMTF);
nMTF = length(mTFs);
Cons = unique(tCon);
nCon = length(Cons);

trial_ind = cell(nDir,nTF,nMTF,nCon);


for iDir = 1:nDir
    ind_dir = find(tDir == Dirs(iDir));
    for iTF = 1:nTF
        ind_tf = find(tTF == TFs(iTF));
        for iMTF = 1:nMTF
            ind_mtf = find(tMTF == mTFs(iMTF));
            for iCon = 1:nCon
                ind_con = find(tCon == Cons(iCon));
                trial_ind{iDir,iTF,iMTF,iCon} = intersect(intersect(intersect(ind_dir,ind_tf),ind_mtf),ind_con);
            end
        end
    end
end

grating_ind = cell(nDir,nTF);
plaid_ind = cell(nDir,nTF,nMTF);
mDirs = circshift(1:nDir,-4);
for iDir = 1:nDir
    for iTF = 1:nTF
        for iMTF = 1:nMTF
            grating_ind{iDir,iTF} = [grating_ind{iDir,iTF} trial_ind{iDir,iTF,iMTF,1}];
            plaid_ind{iDir,iTF,iMTF} = trial_ind{iDir,iTF,iMTF,2};
        end
    end
end

resp_cell = cell(nDir,nTF,nMTF,nCon);
base_cell = cell(nDir,nTF,nMTF,nCon);
trialsperstim = zeros(nDir,nTF,nMTF,nCon);
h_resp =zeros(nCells,nDir,nTF);
p_resp =zeros(nCells,nDir,nTF);
resp_win = prewin_frames+5:prewin_frames+nFramesOn;
base_win = prewin_frames:prewin_frames;
avg_resp_dir = zeros(nCells, nStimDir,nTF,nMTF,nCon, 2);
nStim = nStimDir;
for iDir = 1:nDir
    for iTF = 1:nTF
        for iMTF = 1:nMTF
            for iCon = 1:nCon
                if iCon == 1 && iMTF == 1
                    trialsperstim(iDir,iTF,iMTF,iCon) = length(grating_ind{iDir,iTF});
                    resp_cell{iDir,iTF,iMTF,iCon} = squeeze(mean(data_dfof_tc(resp_win,:,grating_ind{iDir,iTF}),1));
                    base_cell{iDir,iTF,iMTF,iCon} = squeeze(mean(data_dfof_tc(base_win,:,grating_ind{iDir,iTF}),1));
                    [h_resp(:,iDir,iTF), p_resp(:,iDir,iTF)] = ttest2(resp_cell{iDir,iTF,iMTF,iCon},base_cell{iDir,iTF,iMTF,iCon},'dim',2,'tail','right','alpha', 0.05./(nDir*nTF));
                elseif iCon == 1 && iMTF == 2
                    trialsperstim(iDir,iTF,iMTF,iCon) = 0;
                    resp_cell{iDir,iTF,iMTF,iCon} = nan(size(resp_cell{iDir,iTF,1,iCon}));
                    base_cell{iDir,iTF,iMTF,iCon} = nan(size(resp_cell{iDir,iTF,1,iCon}));
                elseif iCon == 2
                    trialsperstim(iDir,iTF,iMTF,iCon) = length(plaid_ind{iDir,iTF,iMTF});
                    resp_cell{iDir,iTF,iMTF,iCon} = squeeze(mean(data_dfof_tc(resp_win,:,plaid_ind{iDir,iTF,iMTF}),1));
                    base_cell{iDir,iTF,iMTF,iCon} = squeeze(mean(data_dfof_tc(base_win,:,plaid_ind{iDir,iTF,iMTF}),1));
                end
                avg_resp_dir(:,iDir,iTF,iMTF,iCon,1) = mean(resp_cell{iDir,iTF,iMTF,iCon},2);
                avg_resp_dir(:,iDir,iTF,iMTF,iCon,2) = std(resp_cell{iDir,iTF,iMTF,iCon},[],2)./sqrt(size(resp_cell{iDir,iTF,iMTF,iCon},2));
            end
        end
    end
end

resp_ind = find(sum(sum(h_resp,2),3));

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']), 'resp_cell', 'data_dfof_tc', 'resp_ind', 'frame_rate', 'h_resp', 'avg_resp_dir');
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']), 'prewin_frames', 'postwin_frames', 'resp_win', 'plaid_ind','grating_ind','trial_ind');

%%

b_dir = nan(nCells,nTF,nMTF,nCon);
k1_dir = nan(nCells,nTF,nMTF,nCon);
R1_dir = nan(nCells,nTF,nMTF,nCon);
R2_dir = nan(nCells,nTF,nMTF,nCon);
u1_dir = nan(nCells,nTF,nMTF,nCon);
sse_dir = nan(nCells,nTF,nMTF,nCon);
R_square_dir = nan(nCells,nTF,nMTF,nCon);
range = 1:1:360;
y_dir_fit = nan(nCells,length(range),nTF,nMTF,nCon);

for iCell = 1:nCells
    for iTF  = 1:nTF
        for iMTF = 1:nMTF
            for iCon = 1:nCon
                if ~isnan(avg_resp_dir(iCell,1,iTF,iMTF,iCon,1))
                    data = [avg_resp_dir(iCell,:,iTF,iMTF,iCon,1) avg_resp_dir(iCell,1,iTF,iMTF,iCon,1)];
                    theta = [deg2rad(stimDirs) 2.*pi];
                    [b_dir(iCell,iTF,iMTF,iCon),k1_dir(iCell,iTF,iMTF,iCon),R1_dir(iCell,iTF,iMTF,iCon), ...
                        R2_dir(iCell,iTF,iMTF,iCon),u1_dir(iCell,iTF,iMTF,iCon),u2_dir(iCell,iTF,iMTF,iCon),sse_dir(iCell,iTF,iMTF,iCon),R_square_dir(iCell,iTF,iMTF,iCon)] ...
                        = miaovonmisesfit_dir(theta,data);
                    y_dir_fit(iCell,:,iTF,iMTF,iCon) = b_dir(iCell,iTF,iMTF,iCon)+R1_dir(iCell,iTF,iMTF,iCon).*exp(k1_dir(iCell,iTF,iMTF,iCon).*(cos(deg2rad(range)-u1_dir(iCell,iTF,iMTF,iCon))-1))...
                        +R2_dir(iCell,iTF,iMTF,iCon).*exp(k1_dir(iCell,iTF,iMTF,iCon).*(cos(deg2rad(range)-u1_dir(iCell,iTF,iMTF,iCon)-pi)-1));
                end
            end
        end
    end
end

DSI = (R1_dir-R2_dir)./(R1_dir+R2_dir);
dsi_grating_ind =  cell(1,nTF);
dsi_plaid_ind =  cell(1,nTF);
for iTF = 1:nTF
    dsi_grating_ind{iTF} = find(DSI(:,iTF,1,1)>0.5);
    if iTF == 1
        dsi_plaid_ind{iTF} = find(DSI(:,iTF,2,2)>0.5);
    else
        dsi_plaid_ind{iTF} = find(DSI(:,iTF,1,2)>0.5);
    end
end


TFI = (R1_dir(:,1,1,1)-R1_dir(:,2,1,1))./(R1_dir(:,1,1,1)+R1_dir(:,2,1,1));
TF_pref = nan(1,nCells);
TF_pref(find(TFI>0)) = 1;
TF_pref(find(TFI<0)) = 2;

maskDiff = double(input.tMaskOneGratingDirectionDeg{1} - input.tStimOneGratingDirectionDeg{1});
int = unique(diff(stimDirs));
component = nan(nCells,nDir,nTF);
pattern = nan(nCells,nDir,nTF);
for iTF = 1:nTF
    component(:,:,iTF) = avg_resp_dir(:,:,iTF,1,1,1)+circshift(avg_resp_dir(:,:,iTF,1,1,1),-maskDiff./int,2);
    pattern(:,:,iTF) = circshift(avg_resp_dir(:,:,iTF,1,1,1),-(maskDiff/2)./int,2);
end

comp_corr = zeros(nCells,nTF);
patt_corr = zeros(nCells,nTF);
comp_patt_corr = zeros(nCells,nTF);

for iCell = 1:nCells
    for iTF = 1:nTF
        comp_corr(iCell,iTF) = triu2vec(corrcoef(avg_resp_dir(iCell,:,iTF,iTF,2,1),component(iCell,:,iTF)));
        patt_corr(iCell,iTF) = triu2vec(corrcoef(avg_resp_dir(iCell,:,iTF,iTF,2,1),pattern(iCell,:,iTF)));
        comp_patt_corr(iCell,iTF) = triu2vec(corrcoef(component(iCell,:,iTF),pattern(iCell,:,iTF)));
    end
end
Rp = ((patt_corr)-(comp_corr.*comp_patt_corr))./sqrt((1-comp_corr.^2).*(1-comp_patt_corr.^2));
Rc = ((comp_corr)-(patt_corr.*comp_patt_corr))./sqrt((1-patt_corr.^2).*(1-comp_patt_corr.^2));
Zp = (0.5.*log((1+Rp)./(1-Rp)))./sqrt(1./(nStimDir-3));
Zc = (0.5.*log((1+Rc)./(1-Rc)))./sqrt(1./(nStimDir-3));

figure; 
movegui('center')
scatter(Zc(intersect(dsi_grating_ind{1},resp_ind),1), Zp(intersect(dsi_grating_ind{1},resp_ind),1))
hold on
scatter(Zc(intersect(dsi_grating_ind{2},resp_ind),2), Zp(intersect(dsi_grating_ind{2},resp_ind),2))
xlabel('Zc')
ylabel('Zp')
ylim([-4 8])
xlim([-4 8])
hold on
plotZcZpBorders
legend({['TF = ' num2str(TFs(1))], ['TF = ' num2str(TFs(2))]})
sgtitle([date ' ' mouse])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_ZcZp.pdf']),'-dpdf', '-fillpage')       

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirAnalysis.mat']), 'b_dir','k1_dir','R1_dir','R2_dir','u1_dir','sse_dir','R_square_dir','y_dir_fit','DSI','dsi_grating_ind','dsi_plaid_ind','TF_pref','TFI','Zp','Zc');

%%
[ioc_ang ioc_sp av_ang av_sp] = iocVectorCalc(deg2rad(0),1,deg2rad(maskDiff),4);

u1_g = zeros(1,nCells);
u1_p = zeros(1,nCells);
grating_ind = [];
plaid_ind = [];
for iTF = 1:nTF
    u1_g(find(TF_pref==iTF)) = u1_dir(find(TF_pref==iTF),iTF,1,1);
    grating_ind = [grating_ind; intersect(resp_ind,dsi_grating_ind{iTF})];
    if iTF == 1
        u1_p = u1_dir(:,iTF,2,2);
        plaid_ind = intersect(resp_ind,dsi_plaid_ind{iTF});
    end
end

prefDiff = u1_g-u1_p;
indUse = intersect(grating_ind,plaid_ind);

figure(1); 
subplot(2,2,1)
scatter(rad2deg(u1_g(indUse)),rad2deg(u1_p(indUse)))
refline(1)
refline(1,-maskDiff)
refline(1,-rad2deg(ioc_ang))
refline(1,-rad2deg(av_ang))
xlabel('Pref dir- Grating')
ylabel('Pref dir- Plaid')
xlim([0 360])
ylim([0 360])
axis square
title(['Plaid angle: ' num2str(maskDiff)])
legend('data','slow comp','fast comp','IOC','VA','Location','northeastoutside')

scale = 20;
subplot(2,2,2)
polarhistogram(prefDiff(indUse),32);
hold on
polarplot([0;0]*pi/180,[0;1]*scale)
polarplot([maskDiff;maskDiff]*pi/180,[0;1]*scale)
polarplot([rad2deg(ioc_ang);rad2deg(ioc_ang)]*pi/180,[0;1]*scale)
polarplot([rad2deg(av_ang);rad2deg(av_ang)]*pi/180,[0;1]*scale)


diff_thresh = 7.5;
prefDiff_deg = rad2deg(prefDiff);
ind_fast = intersect(indUse, intersect(find(prefDiff_deg>-diff_thresh), find(prefDiff_deg<diff_thresh)));
ind_slow = intersect(indUse,  intersect(find(prefDiff_deg>(maskDiff-diff_thresh)), find(prefDiff_deg<(maskDiff+diff_thresh))));
ind_ioc = intersect(indUse,  intersect(find(prefDiff_deg>(rad2deg(ioc_ang)-diff_thresh)), find(prefDiff_deg<(rad2deg(ioc_ang)+diff_thresh))));
ind_av = intersect(indUse,  intersect(find(prefDiff_deg>(rad2deg(av_ang)-diff_thresh)), find(prefDiff_deg<(rad2deg(av_ang)+diff_thresh))));
ind_fast_av = intersect(ind_fast,ind_av);
ind_fast = setdiff(ind_fast, ind_fast_av);
ind_av = setdiff(ind_av, ind_fast_av);

[prefResp prefDir] = max(avg_resp_dir(:,:,1,1,1),[],2);
nCells = size(avg_resp_dir,1);
avg_resp_dir_all_shift = zeros(size(avg_resp_dir(:,:,:,:,1)));
for i = 1:nCells
    avg_resp_dir_all_shift(i,:,:,:,1) = circshift(avg_resp_dir(i,:,:,:,1),1-prefDir(i),2);
end
avg_resp_dir_all_shift_circ = cat(2,avg_resp_dir_all_shift, avg_resp_dir_all_shift(:,1,:,:,1));
figure(2);
subplot(2,2,1)
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_fast,:,1,1),1))
hold on
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_fast,:,2,1),1))
title(['Fast comp- n =' num2str(length(ind_fast))])
rlim([0 0.2])
subplot(2,2,2)
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_slow,:,1,1),1))
hold on
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_slow,:,2,1),1))
title(['Slow comp- n =' num2str(length(ind_slow))])
rlim([0 0.2])
subplot(2,2,3)
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_ioc,:,1,1),1))
hold on
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_ioc,:,2,1),1))
title(['IOC- n =' num2str(length(ind_ioc))])
rlim([0 0.2])
subplot(2,2,4)
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_av,:,1,1),1))
hold on
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_av,:,2,1),1))
title(['VA- n =' num2str(length(ind_av))])
rlim([0 0.2])
sgtitle([num2str(maskDiff) ' deg'])
print(fullfile(fnout,['crossOriRandDir_Type2_' num2str(maskDiff) 'deg_TuningCurves.pdf']), '-fillpage', '-dpdf')

figure(3)
subplot(1,2,1)
piechart([length(ind_fast) length(ind_slow) length(ind_ioc) length(ind_av) length(ind_fast_av)],{'Fast','Slow','IOC','VA','Fast+AV'})
title([num2str(maskDiff) ' deg'])