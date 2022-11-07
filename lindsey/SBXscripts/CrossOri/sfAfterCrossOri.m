%% get path names
close all;clear all;clc;

ds = 'CrossOriRandDir_16x16_ExptList';
eval(ds)
nexp = length(expt);
iexp = 4;
rc = behavConstsAV;

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};
ImgFolder = expt(iexp).sfFolder;
coFolder = expt(iexp).coFolder;
time = expt(iexp).sfTime;
nrun = length(ImgFolder);
frameRateHz = params.frameRate;

run_str = catRunName(cell2mat(ImgFolder), nrun);
co_run_str = catRunName(cell2mat(coFolder), nrun);

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
%LG_base = '\\CRASH.dhe.duke.edu\data\home\lindsey';

fprintf(['2P imaging Dir analysis \nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([ImgFolder{irun} ' - ' time{irun} '\n'])
end

%% load
tic
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = [LG_base '\Data\2P_images\' date '_' mouse '\' ImgFolder{irun}];
    CD = [LG_base '\Data\2P_images\' mouse '\' date '\' ImgFolder{irun}];
    cd(CD);
    imgMatFile = [ImgFolder{irun} '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time{irun} '.mat'];
    load(fName);
    
    temp(irun) = input;
    nOn = temp(irun).nScansOn;
    nOff = temp(irun).nScansOff;
    ntrials = size(temp(irun).tGratingDirectionDeg,2);
    nframes = ntrials*(nOn+nOff);
    
    
    fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n'])
    data_temp = sbxread(imgMatFile(1,1:11),0,nframes);
    if size(data_temp,1)== 2
        data_temp = data_temp(1,:,:,:);
    end
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    fprintf('Complete')
end
input = concatenateDataBlocks(temp);
fprintf('\nAll runs read\n')
fprintf([num2str(size(data,3)) ' total frames\n'])
clear data_temp
clear temp

toc

% register to cross-ori experiment

load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_reg_shifts.mat']))
[out, data_reg] = stackRegister(data,data_avg);
data_reg_avg = mean(data_reg,3);
mkdir(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg', 'data_reg_avg')
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
clear data
% test stability
figure; 
subplot(2,2,1);
imagesc(data_reg_avg);
title('SF run avg')
subplot(2,2,2);
imagesc(data_avg)
title('Cross-ori run avg')
sz = size(data_avg);
rgb = zeros(sz(1),sz(2),3);
rgb(:,:,1) = data_reg_avg./max(data_reg_avg(:));
rgb(:,:,2) = data_avg./max(data_avg(:));
subplot(2,2,3);
image(rgb)
title('Overlay')

print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

% use cross-ori mask to get TCs

fprintf(['Loading masks from cross-ori runs: ' cell2mat(coFolder) '\n'])

load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_mask_cell.mat']))
fprintf('Cell and neuropil masks loaded\n')

nCells = max(mask_cell(:)); % take max label of mask_cell, should circumvent bwlabel
fprintf([num2str(nCells) ' total cells selected\n'])
fprintf('Cell segmentation complete\n')

%neuropil subtraction
down = 5;
sz = size(data_reg);

data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_reg_down  = stackGroupProject(data_reg,down);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
nCells = size(data_tc,2);
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '%s/n']) 
end
%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
[max_skew ind] =  max(x,[],1);
np_w = 0.01*ind;
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
clear data_reg data_reg_down

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

fprintf('\nNeuropil subtraction complete\n')

clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell

%% Direction x SF analysis
%load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
%load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
nOn = input.nScansOn;
nOff = input.nScansOff;
ntrials = length(input.tGratingDirectionDeg);
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

dir_mat = celleqel2mat_padded(input.tGratingDirectionDeg);
dirs = unique(dir_mat);
nDir = length(dirs);

SF_mat = celleqel2mat_padded(input.tGratingSpatialFreqCPD); %transforms cell array into matrix (1 x ntrials)
SFs = unique(SF_mat);
nSF = length(SFs);

base_win = [nOff-5:nOff];
resp_win = [nOff+4:nOff+nOn];
tt_dir = (1-nOff:nOn).*(1000/frameRateHz);

data_dfof_stim_dir = zeros(nDir,nSF,nCells);
h_dir_sf = zeros(nDir,nSF,nCells);
p_dir_sf = zeros(nDir,nSF,nCells);
h_dir = zeros(nDir,nCells);
p_dir = zeros(nDir,nCells);
trialInd = cell(nDir,nSF);
trial_n = zeros(nDir,nSF);
for iDir = 1:nDir
    ind_dir = find(dir_mat == dirs(iDir));
    [h_dir(iDir,:) p_dir(iDir,:)] = ttest(data_dfof_resp(ind_dir,:), data_dfof_base(ind_dir,:),'dim',1, 'tail', 'right', 'alpha', 0.05./(nDir-1));
    for iSF = 1:nSF
        ind_SF = find(SF_mat == SFs(iSF));
        trialInd{iDir,iSF} = intersect(ind_dir,ind_SF);
        trial_n(iDir,iSF) = length(trialInd{iDir,iSF});
        data_dfof_stim_dir(iDir, iSF, :) = mean(data_dfof_resp(trialInd{iDir,iSF},:),1);
        [h_dir_sf(iDir,iSF,:) p_dir_sf(iDir,iSF,:)] = ttest(data_dfof_resp(trialInd{iDir,iSF},:), data_dfof_base(trialInd{iDir,iSF},:),'dim',1, 'tail', 'right', 'alpha', 0.05./(nDir-1));
    end
end

h_dir_sf_all = find(sum(sum(h_dir_sf,1),2));
h_dir_all = find(sum(h_dir,1));

[max_dir_val max_dir_ind] = max(squeeze(mean(data_dfof_stim_dir,2)),[],1);
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
    [max_sf_val(1,iCell) max_sf_ind(1,iCell)] = max(data_dfof_stim_dir(max_dir_ind(iCell),:,iCell),[],2);
    [s i] = sort(data_dfof_stim_dir(max_dir_ind(iCell),:,iCell));
    if abs(i(nSF)-i(nSF-1))<3  || s(nSF-1)./s(nSF)<0.4
        x = log2([0.01 SFs 0.64]);
        options = fitoptions('gauss1', 'Lower', [s(end) x(i(nSF)) 0], 'Upper', [s(end).*1.5 x(i(nSF)+2) Inf]);
        [fit_out{iCell} g_fit{iCell}] =fit(log2(SFs)',data_dfof_stim_dir(max_dir_ind(iCell),:,iCell)','gauss1',options);
    else
        fit_out{iCell}.b1 = NaN;
        g_fit{iCell}.rsquare = NaN;
    end
    subplot(7,7,start)
    plot(log2(SFs),data_dfof_stim_dir(max_dir_ind(iCell),:,iCell),'o')
    hold on
    if ~isnan(fit_out{iCell}.b1)
        plot(fit_out{iCell})
    end
    prefSF(1,iCell)= fit_out{iCell}.b1;
    RsqSF(1,iCell)= g_fit{iCell}.rsquare;
    if find(h_dir_all == iCell)
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
subplot(2,1,1); histogram(RsqSF(h_dir_all)); vline(0.8)
ind_sf = intersect(find(RsqSF>0.8),h_dir_all); xlabel('Rsq')
subplot(2,1,2); histogram(prefSF_cut(ind_sf))
set(gca, 'XTick', log2(SFs), 'XTickLabels', SFs)
vline(mean(prefSF_cut(ind_sf)))
xlabel('SF (cpd)')

suptitle([date ' ' mouse])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_SFfitDist.pdf']), '-dpdf')

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']), 'SF_mat', 'SFs','nSF','dir_mat','dirs','nDir')
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_SFfits.mat']), 'fit_out','g_fit', 'prefSF', 'RsqSF','data_dfof_resp','data_dfof_resp','h_dir_all', 'h_dir', 'h_dir_sf_all','h_dir_sf', 'trialInd','max_sf_ind')

%%
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_respData.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_input.mat']))
resp_ind_anova = intersect(h_dir_all,find(p_anova_dir<0.05));
stimDirs = 0:30:330;
stimDirs_temp = [stimDirs 360];
nStimDir = length(stimDirs);
shift = nStimDir/2;
f(1) = figure;
i = 1;
n=1;
for iCell = 1:length(resp_ind_anova)
    iC = resp_ind_anova(iCell);
    if i == 37
        movegui('center')
        orient(f(n),'landscape')
        sgtitle([mouse ' ' date ' -SF = ' num2str(input.stimOneGratingSpatialFreqCPD)])
        print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_CrossOri&SF_tuningCurves_' num2str(n) '.pdf']),'-dpdf','-bestfit')
        n=1+n;
        i = 1;
        f(n) = figure;
    end
    [max_val max_ind] = max(diag(squeeze(avg_resp_dir(iC,:,:,1))));
    temp_sq = squeeze(circshift(circshift(avg_resp_dir(iC,:,:,1),shift-max_ind,2),shift-max_ind,3));
    subplot(6,6,i)
    imagesc(smooth2(flipud(temp_sq),'gauss',[3 3],3./sqrt(12)))
    axis square
    set(gca, 'XTick',1:nStimDir/4:nStimDir+1,'XTickLabels',stimDirs_temp(1:nStimDir/4:end)-180,'YTick',1:nStimDir/4:nStimDir+1,'YTickLabels',fliplr(stimDirs_temp(1:nStimDir/4:end)-180))
    title(num2str(iC))
    i = 1+i;
    subplot(6,6,i)
    imagesc(smooth2(flipud(circshift(data_dfof_stim_dir(:,:,iCell)',shift-max_ind,2)),'gauss',[3 3],3./sqrt(12)))
    set(gca, 'XTick',1:nStimDir/4:nStimDir+1,'XTickLabels',stimDirs_temp(1:nStimDir/4:end)-180,'YTick',1:nSF,'YTickLabels',fliplr(SFs))
    i = 1+i;
end
movegui('center')
orient(f(n),'landscape') 
sgtitle([mouse ' ' date ' -SF = ' num2str(input.stimOneGratingSpatialFreqCPD)])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_CrossOri&SF_tuningCurves_' num2str(n) '.pdf']),'-dpdf','-bestfit')
        

%%
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_dirAnalysis.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_respData.mat']))

zp_sf = zeros(nSF,2);
zc_sf = zeros(nSF,2);
for i = 1:nSF
    ind = intersect(intersect(ind_sf,resp_ind),find(max_sf_ind==i));
    zp_sf(i,1) = mean(Zp(ind));
    zp_sf(i,2) = std(Zp(ind))./sqrt(length(ind));
    zc_sf(i,1) = mean(Zc(ind));
    zc_sf(i,2) = std(Zc(ind))./sqrt(length(ind));
end
figure; movegui('center')
subplot(1,3,1)
scatter(log2(SFs(max_sf_ind(intersect(ind_sf,resp_ind)))),Zc(intersect(ind_sf,resp_ind)))
xlabel('Preferred SF')
ylabel('Zc')
ylim([-4 4])
xlim([-6 -1])
subplot(1,3,2)
scatter(log2(SFs(max_sf_ind(intersect(ind_sf,resp_ind)))),Zp(intersect(ind_sf,resp_ind)))
xlabel('Preferred SF')
ylabel('Zp')
ylim([-4 4])
xlim([-6 -1])
subplot(1,3,3)
errorbar(log2(SFs),zc_sf(:,1),zc_sf(:,2),'o')
hold on
errorbar(log2(SFs),zp_sf(:,1),zp_sf(:,2),'o')
xlabel('Preferred SF')
ylabel('Zc/Zp')
ylim([-4 4])
xlim([-6 -1])

figure; movegui('center')
subplot(1,3,1)
scatter(prefSF_cut(intersect(ind_sf,resp_ind)),Zc(intersect(ind_sf,resp_ind)))
xlabel('Peak SF')
ylabel('Zc')
ylim([-4 4])
xlim([-6 -1])
subplot(1,3,2)
scatter(prefSF_cut(intersect(ind_sf,resp_ind)),Zp(intersect(ind_sf,resp_ind)))
xlabel('Peak SF')
ylabel('Zp')
ylim([-4 4])
xlim([-6 -1])
subplot(1,3,3)
errorbar(log2(SFs),zc_sf(:,1),zc_sf(:,2),'o')
hold on
errorbar(log2(SFs),zp_sf(:,1),zp_sf(:,2),'o')
xlabel('Preferred SF')
ylabel('Zc/Zp')
ylim([-4 4])
xlim([-6 -1])
legend({'Zc','Zp'})
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_ZpZc_SF.pdf']),'-dpdf', '-bestfit')

figure;
scatter(prefSF_cut(resp_ind),Zc(resp_ind))
hold on
scatter(prefSF_cut(resp_ind),Zp(resp_ind))
xlabel('Preferred SF')
ylabel('Zc/Zp')
ylim([-4 4])

%%
Ori = dir_mat;
Ori(find(dir_mat>=180)) = dir_mat(find(dir_mat>=180))-180;
Oris = unique(Ori);
nOris = length(Oris);

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
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_OriTuningByStimSF.pdf']), '-dpdf','-fillpage')

figure; movegui('center')
for iSF = 1:nSF
    errorbar(Oris, mean(data_dfof_stim_ori_align(:,iSF,h_ori_all),3)./mean(data_dfof_stim_ori_align(1,iSF,h_ori_all),3), std(data_dfof_stim_ori_align(:,iSF,h_ori_all),[],3)./sqrt(length(h_ori_all))./std(data_dfof_stim_ori_align(1,iSF,h_ori_all),[],3)./sqrt(length(h_ori_all)),'-o')
    hold on
    
end
ylim([0 1.25])
xlim([-10 180])    
xlabel('Ori')
ylabel('dF/F')
legend(num2str(SFs'))
suptitle([date ' ' mouse '- n = ' num2str(length(h_ori_all)) ' cells'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_OriTuningByStimSF_Norm.pdf']), '-dpdf','-fillpage')

figure; movegui('center')
for iSF = 1:nSF
    ind_sf = intersect(h_ori_all,find(max_sf_ind==iSF));
    %subplot(3,3,iSF)
    errorbar(Oris, mean(data_dfof_stim_ori_align(:,iSF,ind_sf),3), std(data_dfof_stim_ori_align(:,iSF,ind_sf),[],3)./sqrt(length(ind_sf)),'-o')
    hold on
    %title([num2str(SFs(iSF)) ' CPD; n = ' num2str(length(ind_sf))])
    legend(num2str(SFs))
    ylim([0 1.2])
    xlabel('Ori')
    ylabel('dF/F')
end
suptitle([date ' ' mouse])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_OriTuningByStimSF_atPrefSF.pdf']), '-dpdf','-fillpage')

data_dfof_stim_ori_norm = data_dfof_stim_ori_align./data_dfof_stim_ori_align(1,:,:);
figure; movegui('center')
for iSF = 1:nSF
    ind_sf = intersect(h_ori_all,find(max_sf_ind==iSF));
    subplot(3,3,iSF)
    errorbar(Oris, mean(data_dfof_stim_ori_norm(:,1,ind_sf),3), std(data_dfof_stim_ori_norm(:,1,ind_sf),[],3)./sqrt(length(ind_sf)),'-o')
    hold on
    errorbar(Oris, mean(data_dfof_stim_ori_norm(:,3,ind_sf),3), std(data_dfof_stim_ori_norm(:,3,ind_sf),[],3)./sqrt(length(ind_sf)),'-o')    
    title([num2str(SFs(iSF)) ' CPD; n = ' num2str(length(ind_sf))])
    ylim([0 2.5])
    xlabel('Ori')
    ylabel('Norm dF/F')
    if iSF == 1
        legend({num2str(SFs(1)), num2str(SFs(3))})
    end
end
suptitle([date ' ' mouse])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_OriTuningByStimSF_atHighLowSF.pdf']), '-dpdf','-fillpage')
