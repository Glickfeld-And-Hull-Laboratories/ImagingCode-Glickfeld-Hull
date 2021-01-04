%Path names
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
lg_fn = fullfile(fn_base, 'home\lindsey');
data_fn = fullfile(lg_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data');
fnout = fullfile(lg_fn, 'Analysis\2P');

%Specific experiment information
date = '200731';
ImgFolder = '002';
time = '1126';
mouse = 'i1323';
frame_rate = 15;
run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];

nOn = input.nScansOn;
nOff = input.nScansOff;

Dir = celleqel2mat_padded(input.tGratingDirectionDeg); %transforms cell array into matrix (1 x ntrials)
Dirs = unique(Dir);
nDirs = length(Dirs);

SF_mat = celleqel2mat_padded(input.tGratingSpatialFreqCPD); %transforms cell array into matrix (1 x ntrials)
SFs = unique(SF_mat);
nSF = length(SFs);

Ori = Dir;
Ori(find(Dir>=180)) = Dir(find(Dir>=180))-180;
Oris = unique(Ori);
nOris = length(Oris);

ntrials = length(Dir);
nCells = size(npSub_tc,2);

frame_rate = 15;
prewin_frames = nOff-10:nOff;
postwin_frames = nOff+10:nOff+nOn;
tt = (1-nOff:nOn).*(1/frame_rate);
data_tr = reshape(npSub_tc,[nOn+nOff ntrials nCells]);

data_f = mean(data_tr(prewin_frames,:,:),1);
data_dfof_tc = (data_tr-data_f)./data_f;
data_dfof_resp = squeeze(mean(data_dfof_tc(postwin_frames,:,:),1));
data_dfof_base = squeeze(mean(data_dfof_tc(prewin_frames,:,:),1));

data_dfof_stim_dir = zeros(nDirs,nSF,nCells);
h_dir = zeros(nDirs,nSF,nCells);
p_dir = zeros(nDirs,nSF,nCells);
trialInd = cell(nDirs,nSF);
for iDir = 1:nDirs
    ind_dir = find(Dir == Dirs(iDir));
    for iSF = 1:nSF
        ind_SF = find(SF_mat == SFs(iSF));
        trialInd{iDir,iSF} = intersect(ind_dir,ind_SF);
        data_dfof_stim(iDir, iSF, :) = mean(data_dfof_resp(trialInd{iDir,iSF},:),1);
        [h_dir(iDir,iSF,:) p_dir(iDir,iSF,:)] = ttest(data_dfof_resp(trialInd{iDir,iSF},:), data_dfof_base(trialInd{iDir,iSF},:),'dim',1, 'tail', 'right', 'alpha', 0.05./((nDirs.*nSF)-1));
    end
end

h_dir_all = find(sum(sum(h_dir,1),2));

data_dfof_stim_ori = zeros(nOris,nSF,nCells);
h_ori = zeros(nOris,nSF,nCells);
p_ori = zeros(nOris,nSF,nCells);
trialInd = cell(nOris,nSF);
for iOri= 1:nOris
    ind_ori = find(Ori == Oris(iOri));
    for iSF = 1:nSF
        ind_SF = find(SF_mat == SFs(iSF));
        trialInd{iOri,iSF} = intersect(ind_ori,ind_SF);
        data_dfof_stim_ori(iOri, iSF, :) = mean(data_dfof_resp(trialInd{iOri,iSF},:),1);
        [h_ori(iOri,iSF,:) p_ori(iOri,iSF,:)] = ttest(data_dfof_resp(trialInd{iOri,iSF},:), data_dfof_base(trialInd{iOri,iSF},:),'dim',1, 'tail', 'right', 'alpha', 0.05./((nOris.*nSF)-1));
    end
end

h_ori_all = find(sum(sum(h_ori,1),2));

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
    imagesc(data_dfof_stim_dir(:,:,iCell))
    colormap gray
    if find(h_ori_all == iCell)
        title('Sig')
    end
    start = start+1;
end

[max_dir_val max_dir_ind] = max(squeeze(mean(data_dfof_stim_dir,2)),[],1);
max_sf_val = zeros(1,nCells);
max_sf_ind = zeros(1,nCells);
fit_out = cell(1,nCells);
g_fit = cell(1,nCells);
prefSF = zeros(1,nCells);
figure; movegui('center')
start = 1;
for iCell = 1:nCells
%     if start >49
%         figure;
%         start = 1;
%     end
    [max_sf_val(1,iCell) max_sf_ind(1,iCell)] = max(data_dfof_stim_dir(max_dir_ind(iCell),:,iCell),[],2);
    [fit_out{iCell} g_fit{iCell}] =fit(log2(SFs)',data_dfof_stim_dir(max_dir_ind(iCell),:,iCell)','gauss1');
%     subplot(7,7,start)
%     plot(log2(SFs),data_dfof_stim(max_dir_ind(iCell),:,iCell),'o')
%     hold on
%     plot(fit_out{iCell})
    prefSF(1,iCell)= fit_out{iCell}.b1;
    RsqSF(1,iCell)= g_fit{iCell}.rsquare;
%     if find(h_all == iCell)
%         title('Sig')
%     end
%     start = start+1;
%     legend('off')
    if rem(iCell, 10) == 0
        fprintf([num2str(iCell) '\n'])
    end
end
prefSF_cut = prefSF;
prefSF_cut(find(prefSF>max(log2(SFs)))) = max(log2(SFs));
prefSF_cut(find(prefSF<min(log2(SFs)))) = min(log2(SFs));
figure; movegui('center');
RsqSF(find(RsqSF<0)) = 0;
subplot(2,1,1); histogram(RsqSF(h_all)); vline(0.8)
ind = intersect(find(RsqSF<0.8),h_all); xlabel('Rsq')
subplot(2,1,2); histogram(prefSF_cut(ind))
set(gca, 'XTick', log2(SFs), 'XTickLabels', SFs)
vline(mean(prefSF_cut(ind),2))
xlabel('SF (cpd)')

suptitle([date ' ' mouse])
% print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_SFfitDist.pdf']), '-dpdf')


[max_val max_ori] = max(mean(data_dfof_stim_ori,2),[],1);
max_ori = squeeze(max_ori);
data_dfof_stim_ori_align = zeros(size(data_dfof_stim_ori));
for iCell = 1:nCells
    data_dfof_stim_ori_align(:,:,iCell) = circshift(data_dfof_stim_ori(:,:,iCell),1-max_ori(iCell),1);
end

figure; movegui('center')
for iSF = 1:nSF
    subplot(3,3,iSF)
    errorbar(Oris, mean(data_dfof_stim_ori_align(:,iSF,h_ori_all),3), std(data_dfof_stim_ori_align(:,iSF,h_ori_all),[],3)./sqrt(length(h_ori_all)),'-o')
    hold on
    title(num2str(SFs(iSF)))
    ylim([0 0.5])
    xlabel('Ori')
    ylabel('dF/F')
end
suptitle([date ' ' mouse '- n = ' num2str(length(h_ori_all)) ' cells'])
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_OriTuningByStimSF.pdf']), '-dpdf','-fillpage')


figure; movegui('center')
for iSF = 1:nSF
    ind_sf = intersect(h_ori_all,find(max_sf_ind==iSF));
    subplot(3,3,iSF)
    errorbar(Oris, mean(data_dfof_stim_ori_align(:,iSF,ind_sf),3), std(data_dfof_stim_ori_align(:,iSF,ind_sf),[],3)./sqrt(length(ind_sf)),'-o')
    title([num2str(SFs(iSF)) ' CPD; n = ' num2str(length(ind_sf))])
    ylim([0 1.2])
    xlabel('Ori')
    ylabel('dF/F')
end
suptitle([date ' ' mouse])
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_OriTuningByStimSF_atPrefSF.pdf']), '-dpdf','-fillpage')

data_dfof_stim_ori_norm = data_dfof_stim_ori_align./data_dfof_stim_ori_align(1,:,:);
figure; movegui('center')
for iSF = 1:nSF
    ind_sf = intersect(h_ori_all,find(max_sf_ind==iSF));
    subplot(3,3,iSF)
    errorbar(Oris, mean(data_dfof_stim_ori_norm(:,2,ind_sf),3), std(data_dfof_stim_ori_norm(:,2,ind_sf),[],3)./sqrt(length(ind_sf)),'-o')
    hold on
    errorbar(Oris, mean(data_dfof_stim_ori_norm(:,5,ind_sf),3), std(data_dfof_stim_ori_norm(:,5,ind_sf),[],3)./sqrt(length(ind_sf)),'-o')    
    title([num2str(SFs(iSF)) ' CPD; n = ' num2str(length(ind_sf))])
    ylim([0 2.5])
    xlabel('Ori')
    ylabel('Norm dF/F')
    if iSF == 1
        legend({num2str(SFs(2)), num2str(SFs(5))})
    end
end
suptitle([date ' ' mouse])
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_OriTuningByStimSF_atHighLowSF.pdf']), '-dpdf','-fillpage')



% save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_SFfits.mat']), 'fit_out','g_fit', 'prefSF', 'RsqSF','data_dfof_resp','data_dfof_stim','h_all', 'h', 'trialInd')

