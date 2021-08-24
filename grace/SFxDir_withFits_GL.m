% load TCs
mouse = 'i1803';
date = '210622';
ImgFolder = strvcat('003');
nrun = size(ImgFolder,1);
frame_rate = 15.5;
% if Day 1 say 1; else say 0
ref_day = 1;
run_str = catRunName(ImgFolder, nrun);
load(fullfile(fnout, [day '_' mouse], [day '_' mouse '_' run_str], [day '_' mouse '_' run_str '_input.mat']));
tGratingDir = celleqel2mat_padded(input.tGratingDirectionDeg);
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P';
TCs = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']));
if ref_day
npSub_tc = TCs.npSub_tc;
else
np = TCs.cellTCs_match;
npSub_tc = np{2};
end
nCells = size(npSub_tc,2);
ntrials = size(tGratingDir,2);
prewin_frames = nOff-10:nOff;
postwin_frames = nOff+10:nOff+nOn;
tt = (1-nOff:nOn).*(1/frame_rate);
data_tr = reshape(npSub_tc,[nOn+nOff ntrials nCells]);

data_f = mean(data_tr(prewin_frames,:,:),1);
data_dfof_tc = (data_tr-data_f)./data_f;
data_dfof_resp = squeeze(mean(data_dfof_tc(postwin_frames,:,:),1));
data_dfof_base = squeeze(mean(data_dfof_tc(prewin_frames,:,:),1));

SF_mat = celleqel2mat_padded(input.tGratingSpatialFreqCPD);
SFs = unique(SF_mat);
nSF = length(SFs);

Dir = celleqel2mat_padded(input.tGratingDirectionDeg); 
Dirs = unique(Dir);
nDirs = length(Dirs);

data_dfof_stim = zeros(nDirs,nSF,nCells);
h = zeros(nDirs,nSF,nCells);
p = zeros(nDirs,nSF,nCells);
trialInd = cell(nDirs,nSF);

for iDir = 1:nDirs
    ind_dir = find(Dir(1:640) == Dirs(iDir));
    for iSF = 1:nSF
        ind_SF = find(SF_mat(1:640) == SFs(iSF));
        trialInd{iDir,iSF} = intersect(ind_dir,ind_SF);
        data_dfof_stim(iDir, iSF, :) = mean(data_dfof_resp(trialInd{iDir,iSF},:),1);
        [h(iDir,iSF,:), p(iDir,iSF,:)] = ttest(data_dfof_resp(trialInd{iDir,iSF},:), data_dfof_base(trialInd{iDir,iSF},:),'dim',1, 'tail', 'right', 'alpha', 0.05./((nDirs.*nSF)-1));
    end
end

h_all = find(sum(sum(h,1),2));

figure;
movegui('center')
start = 1;
n = 1;
for iCell = 1:nCells
    if start>49
        figure;
        movegui('center')
        start = 1;
    end
    subplot(7,7,start)
    imagesc(data_dfof_stim(:,:,iCell))
    colormap gray
    if find(h_all == iCell)
        title('Sig')
    end
    start = start+1;
end

[max_dir_val max_dir_ind] = max(squeeze(mean(data_dfof_stim,2)),[],1);
max_sf_val = zeros(1,nCells);
max_sf_ind = zeros(1,nCells);
fit_out = cell(1,nCells);
g_fit = cell(1,nCells);
prefSF = zeros(1,nCells);
figure; movegui('center')
start = 1;
for iCell = 1:nCells
    if start >49
        figure;
        start = 1;
    end
    [max_sf_val(1,iCell) max_sf_ind(1,iCell)] = max(data_dfof_stim(max_dir_ind(iCell),:,iCell),[],2);
    [fit_out{iCell} g_fit{iCell}] =fit(log2(SFs)',data_dfof_stim(max_dir_ind(iCell),:,iCell)','gauss1');
    subplot(7,7,start)
    plot(log2(SFs),data_dfof_stim(max_dir_ind(iCell),:,iCell),'o')
    hold on
    plot(fit_out{iCell})
    prefSF(1,iCell)= fit_out{iCell}.b1;
    RsqSF(1,iCell)= g_fit{iCell}.rsquare;
    if find(h_all == iCell)
        title('Sig')
    end
    start = start+1;
    legend('off')
    if rem(iCell, 10) == 0
        fprintf([num2str(iCell) '\n'])
    end
end
prefSF_cut = prefSF;
prefSF_cut(find(prefSF>max(log2(SFs)))) = max(log2(SFs));
prefSF_cut(find(prefSF<min(log2(SFs)))) = min(log2(SFs));
figure; movegui('center');
RsqSF(find(RsqSF<0)) = 0;
subplot(2,1,1); histogram(RsqSF(h_all),'BinWidth',0.1); vline(0.8)
ind = intersect(find(RsqSF<0.8),h_all); xlabel('Rsq')
subplot(2,1,2); histogram(prefSF_cut(ind))
set(gca, 'XTick', log2(SFs), 'XTickLabels', SFs)
vline(mean(prefSF_cut(ind),2))
xlabel('SF (cpd)')

suptitle([date ' ' mouse])
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_SFfitDist.pdf']), '-dpdf')


save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_SFfits.mat']), 'fit_out','g_fit', 'prefSF', 'RsqSF','data_dfof_resp','data_dfof_stim','h_all', 'h', 'trialInd')

