close all
clear all
clc

ds = 'ChC_ToneXDir_ExptList';
eval(ds)
nexp = length(expt);

driver_list = unique([expt.driver]);
ndriver = length(driver_list);
driver_n = zeros(1,ndriver);
driver_ind = cell(1,ndriver);
driver_depth = zeros(2,ndriver);
z = [expt.z];
for i = 1:ndriver
    driver_ind{i} = strcmp([expt.driver],driver_list{i});
    driver_n(i) = sum(driver_ind{i});
    driver_depth(1,i) = mean(z(find(driver_ind{i})));
    driver_depth(2,i) = std(z(find(driver_ind{i})))./sqrt(driver_n(i));
end

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey\Analysis\2P';
data_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey\Data\2P_images';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\ChC';
sv_name = 'ToneXDir';

OSI_all = [];
DSI_all = [];
h_dir_use = [];
cell_n = [];
driver_all = [];
data_dfof_dir_all = [];
pupil_resp_all = [];
respTC_all = [];
for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc;
    driver = expt(iexp).driver;
    ImgFolder = expt(iexp).folder;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);
    
    datemouse = [date '_' mouse];
    datemouserun = [date '_' mouse '_' run_str];

    fprintf([datemouse '\n'])
    load(fullfile(LG_base,datemouse,datemouserun,[datemouserun '_oriResp.mat']))

    nCells = size(stim_OSI,2);
    data_dfof_dir_all = cat(1,data_dfof_dir_all, data_dfof_dir);
    OSI_all = [OSI_all; stim_OSI'];
    DSI_all = [DSI_all; stim_DSI'];
    h_dir_use = cat(1,h_dir_use,permute(h_dir,[3 1 2]));
    driver_all = [driver_all; repmat(driver,[nCells 1])];
    cell_n = [cell_n; size(h_dir,3) length(find(sum(sum(h_dir,1),2)))];
    
    if exist(fullfile(LG_base,datemouse,datemouserun,[datemouserun '_pupil.mat']))
        load(fullfile(LG_base,datemouse,datemouserun,[datemouserun '_pupil.mat']))
    else
        pupil_resp = nan(90,2);
    end
    pupil_resp_all = cat(3,pupil_resp_all,pupil_resp);

    load(fullfile(LG_base, datemouse, datemouserun, [datemouserun '_TCs.mat']))
    load(fullfile(LG_base, datemouse, datemouserun, [datemouserun '_input.mat']))
    nOn = unique(celleqel2mat_padded(input.nStimOneFramesOn));
    nOff = unique(celleqel2mat_padded(input.tItiWaitFrames));
    ntrials = size(input.tStimOneGratingDirectionDeg,2); %this is a cell array with one value per trial, so length = ntrials
    cStimOneOn = celleqel2mat_padded(input.cStimOneOn);
    tDir = celleqel2mat_padded(input.tStimOneGratingDirectionDeg); %transforms cell array into matrix (1 x ntrials)
    Dirs = unique(tDir);
    nDirs = length(Dirs);
    tCon = celleqel2mat_padded(input.tStimOneGratingContrast);
    tAud = ~celleqel2mat_padded(input.tBlock2TrialNumber);
    base_win = nOff/2:nOff;
    nCells = size(data_tc,2);
    data_trial = nan(nOn+nOff+nOff,nCells,ntrials);
    for itrial = 1:ntrials-1
        data_trial(:,:,itrial) = npSub_tc(cStimOneOn(itrial)-nOff+1:cStimOneOn(itrial)+nOn+nOff,:);
    end
    data_trial = permute(data_trial, [1 3 2]);
    data_f = mean(data_trial(base_win,:,:),1);
    data_dfof = bsxfun(@rdivide,bsxfun(@minus,data_trial,data_f),data_f);
    respTC = zeros(nOff+nOn+nOff,nCells,nDirs+1,2);
    for iA = 1:2
        for iDir = 1:nDirs
            ind{iDir,iA} = intersect(find(tAud == iA-1),intersect(find(tCon),find(tDir==Dirs(iDir))));
            respTC(:,:,iDir,iA) = squeeze(mean(data_dfof(:,ind{iDir,iA},:),2,'omitnan'));
        end
        ind{iDir+1,iA} = intersect(find(tAud == iA-1),find(tCon==0));
        respTC(:,:,iDir+1,iA) = squeeze(mean(data_dfof(:,ind{iDir+1,iA},:),2,'omitnan'));
    end
    save(fullfile(LG_base, datemouse, datemouserun, [datemouserun '_respTC.mat']),'respTC')
    respTC_all = cat(2,respTC_all,respTC);
end 
%%
totCells = zeros(2,ndriver);
respCells = zeros(2,ndriver);
for i = 1:ndriver
    totCells(1,i) = mean(cell_n(find(driver_ind{i}),1));
    totCells(2,i) = std(cell_n(find(driver_ind{i}),1))./sqrt(driver_n(i));
    respCells(1,i) = mean(cell_n(find(driver_ind{i}),2));
    respCells(2,i) = std(cell_n(find(driver_ind{i}),2))./sqrt(driver_n(i));
end
varNames = {'Total Expts', 'Cells/Expt','RespCells/Expt','Depth'};
cellTab = table(driver_n',chop(totCells(1,:),3)',chop(respCells(1,:),3)',chop(driver_depth(1,:),3)','RowNames',driver_list','VariableNames',varNames);
writetable(cellTab,fullfile(fnout,[sv_name '_cellDataTable.txt']))

nCells = size(OSI_all,1);

load(fullfile(LG_base,datemouse,datemouserun,[datemouserun '_input.mat']))
nOn = unique(celleqel2mat_padded(input.nStimOneFramesOn));
nOff = unique(celleqel2mat_padded(input.tItiWaitFrames));
Dir = celleqel2mat_padded(input.tStimOneGratingDirectionDeg); 
if find(Dir>+360)
    Dir(find(Dir>=360)) = Dir(find(Dir>=360))-360;
end
Dirs = unique(Dir);
nDirs = length(Dirs);
Oris = Dirs(1:nDirs/2);

cellTypes = unique(driver_all);
cellType_ind = cell(size(cellTypes));
cellType_vis_ind = cell(size(cellTypes));
cellType_novis_ind = cell(size(cellTypes));
cellN = zeros(size(cellTypes));
figure; movegui('center') 
for i = 1:length(cellTypes)
    cellType_ind{i} = find(strcmp(driver_all,cellTypes{i}));
    cellType_vis_ind{i} = intersect(find(sum(sum(h_dir_use,2),3)),find(strcmp(driver_all,cellTypes{i})));
    cellType_novis_ind{i} = setdiff(find(strcmp(driver_all,cellTypes{i})),find(sum(sum(h_dir_use,2),3)));
    cellN(i) = length(cellType_vis_ind{i});
    subplot(2,2,1)
    cdfplot(OSI_all(cellType_vis_ind{i}))
    hold on
    subplot(2,2,2)
    cdfplot(DSI_all(cellType_vis_ind{i}))
    hold on
end
subplot(2,2,1)
title('OSI')
xlabel('OSI')
legend(cellTypes,'location','southeast')
subplot(2,2,2)
title('DSI')
xlabel('DSI')
print(fullfile(fnout,[sv_name '_TuningSummary.pdf']),'-dpdf','-fillpage')


data_dfof_dir_shift = zeros(nCells,nDirs,2);
data_dfof_ori_shift = zeros(nCells,nDirs/2,2);
data_dfof_ori_all = squeeze(mean(reshape(data_dfof_dir_all(:,1:nDirs,:,1),[nCells nDirs/2 2 2]),3));
[max_dir_val max_dir_ind] = max(mean(data_dfof_dir_all(:,1:16,:,1),3),[],2);
[max_ori_val max_ori_ind] = max(mean(data_dfof_ori_all(:,:,:),3),[],2);
for iCell = 1:size(OSI_all,1)
    shift = nDirs/2 - max_dir_ind(iCell);
    data_dfof_dir_shift(iCell,:,:) = circshift(data_dfof_dir_all(iCell,1:16,:,1),shift,2);
    ori_shift = nDirs/4 - max_ori_ind(iCell);
    data_dfof_ori_shift(iCell,:,:) = circshift(data_dfof_ori_all(iCell,:,:),ori_shift,2);
end


figure; movegui('center')
norm_val =squeeze(max(max(data_dfof_dir_shift,[],2),[],3));
for ii = 1:length(cellTypes)
    subplot(3,length(cellTypes),ii)
    
    errorbar(Dirs, squeeze(mean(data_dfof_dir_shift(cellType_vis_ind{ii},:,1)./norm_val(cellType_vis_ind{ii}),1)), ...
        squeeze(std(data_dfof_dir_shift(cellType_vis_ind{ii},:,1)./norm_val(cellType_vis_ind{ii}),[],1)./sqrt(length(cellType_vis_ind{ii}))))
    hold on
    errorbar(Dirs, squeeze(mean(data_dfof_dir_shift(cellType_vis_ind{ii},:,2)./norm_val(cellType_vis_ind{ii}),1)), ...
        squeeze(std(data_dfof_dir_shift(cellType_vis_ind{ii},:,2)./norm_val(cellType_vis_ind{ii}),[],1)./sqrt(length(cellType_vis_ind{ii}))))
    xlabel('Direction')
    ylabel('Norm. dF/F')
    ylim([-0.5 1.5])
    title(cellTypes{ii})
end


norm_val =squeeze(max(max(data_dfof_ori_shift,[],2),[],3));
for ii = 1:length(cellTypes)
    subplot(3,length(cellTypes),ii+length(cellTypes))
    errorbar(Oris, squeeze(mean(data_dfof_ori_shift(cellType_vis_ind{ii},:,1)./norm_val(cellType_vis_ind{ii}),1)), ...
        squeeze(std(data_dfof_ori_shift(cellType_vis_ind{ii},:,1)./norm_val(cellType_vis_ind{ii}),[],1)./sqrt(length(cellType_vis_ind{ii}))))
    hold on
    errorbar(Oris, squeeze(mean(data_dfof_ori_shift(cellType_vis_ind{ii},:,2)./norm_val(cellType_vis_ind{ii}),1)), ...
        squeeze(std(data_dfof_ori_shift(cellType_vis_ind{ii},:,2)./norm_val(cellType_vis_ind{ii}),[],1)./sqrt(length(cellType_vis_ind{ii}))))
    xlabel('Orientation')
    ylabel('Norm. dF/F')
    ylim([-0.5 1.5])
end

for ii = 1:length(cellTypes)
    subplot(3,length(cellTypes),ii+(2*length(cellTypes)))
    errorbar(ii, squeeze(mean(data_dfof_dir_all(cellType_vis_ind{ii},end,1)./norm_val(cellType_vis_ind{ii}),1)), ...
        squeeze(std(data_dfof_dir_all(cellType_vis_ind{ii},end,1)./norm_val(cellType_vis_ind{ii}),[],1)./sqrt(length(cellType_vis_ind{ii}))),'o')
    hold on
    errorbar(ii, squeeze(mean(data_dfof_dir_all(cellType_vis_ind{ii},end,2)./norm_val(cellType_vis_ind{ii}),1)), ...
        squeeze(std(data_dfof_dir_all(cellType_vis_ind{ii},end,2)./norm_val(cellType_vis_ind{ii}),[],1)./sqrt(length(cellType_vis_ind{ii}))),'o')
    ylim([-0.1 1])
    ylabel('Norm. dF/F')
    legend({'Blank', 'Tone only'})
end

print(fullfile(fnout,[sv_name '_Dir&OriTuning.pdf']),'-dpdf','-bestfit')

figure;
tt = (1-nOff:nOn).*(1000/params.frameRate);
pupil_resp_all_sub = pupil_resp_all./mean(pupil_resp_all(nOff/2:nOff,:,:),1);
[n n2] = subplotn(nexp+1);
for i = 1:nexp
    subplot(n,n2,i)
    plot(tt,pupil_resp_all_sub(:,1,i))
    hold on
    plot(tt,pupil_resp_all_sub(:,2,i))
    xlabel('Time from stim')
    xlim([-2000 2000])
end
subplot(n,n2,i+1)
errorbar(tt,mean(pupil_resp_all_sub(:,1,:),3,'omitnan'),std(pupil_resp_all_sub(:,1,:),[],3,'omitnan')./sqrt(nexp))
hold on
errorbar(tt,mean(pupil_resp_all_sub(:,2,:),3,'omitnan'),std(pupil_resp_all_sub(:,2,:),[],3,'omitnan')./sqrt(nexp))
xlabel('Time from stim')
xlim([-2000 2000])

print(fullfile(fnout,[sv_name '_pupilDiameterSummary.pdf']),'-dpdf','-fillpage')

figure;
for ii = 1:length(cellTypes)
    subplot(2,length(cellTypes),ii)
    tt = (-nOff:nOn+nOff-1).*(1000./(15));
    errorbar(tt, mean(respTC_all(:,cellType_vis_ind{ii},end,1),2),std(respTC_all(:,cellType_vis_ind{ii},end,1),[],2)./sqrt(length(cellType_vis_ind{1})))
    hold on
    errorbar(tt, mean(respTC_all(:,cellType_vis_ind{ii},end,2),2),std(respTC_all(:,cellType_vis_ind{ii},end,2),[],2)./sqrt(length(cellType_vis_ind{1})))
    ylabel('dF/F')
    xlabel('Time (ms)')
    title(cellTypes{ii})
end

for ii = 1:length(cellTypes)
    subplot(2,length(cellTypes),ii+length(cellTypes))
    tt = (-nOff:nOn+nOff-1).*(1000./(15));
    errorbar(tt, mean(respTC_all(:,cellType_novis_ind{ii},end,1),2),std(respTC_all(:,cellType_novis_ind{ii},end,1),[],2)./sqrt(length(cellType_novis_ind{1})))
    hold on
    errorbar(tt, mean(respTC_all(:,cellType_novis_ind{ii},end,2),2),std(respTC_all(:,cellType_novis_ind{ii},end,2),[],2)./sqrt(length(cellType_novis_ind{1})))
    ylabel('dF/F')
    xlabel('Time (ms)')
end
suptitle("No visual stimulus")

figure;
for ii = 1:length(cellTypes)
    subplot(2,length(cellTypes),ii)
    tt = (-nOff:nOn+nOff-1).*(1000./(15));
    errorbar(tt, mean(mean(respTC_all(:,cellType_vis_ind{ii},1:nDirs,1),3),2),std(mean(respTC_all(:,cellType_vis_ind{ii},1:nDirs,1),3),[],2)./sqrt(length(cellType_vis_ind{1})))
    hold on
    errorbar(tt, mean(mean(respTC_all(:,cellType_vis_ind{ii},1:nDirs,2),3),2),std(mean(respTC_all(:,cellType_vis_ind{ii},1:nDirs,2),3),[],2)./sqrt(length(cellType_vis_ind{1})))
    ylabel('dF/F')
    xlabel('Time (ms)')
    title(cellTypes{ii})
end

for ii = 1:length(cellTypes)
    subplot(2,length(cellTypes),ii+length(cellTypes))
    tt = (-nOff:nOn+nOff-1).*(1000./(15));
     errorbar(tt, mean(mean(respTC_all(:,cellType_novis_ind{ii},1:nDirs,1),3),2),std(mean(respTC_all(:,cellType_novis_ind{ii},1:nDirs,1),3),[],2)./sqrt(length(cellType_novis_ind{1})))
    hold on
    errorbar(tt, mean(mean(respTC_all(:,cellType_novis_ind{ii},1:nDirs,2),3),2),std(mean(respTC_all(:,cellType_novis_ind{ii},1:nDirs,2),3),[],2)./sqrt(length(cellType_novis_ind{1})))
    ylabel('dF/F')
    xlabel('Time (ms)')
end
suptitle("Visual stimulus")