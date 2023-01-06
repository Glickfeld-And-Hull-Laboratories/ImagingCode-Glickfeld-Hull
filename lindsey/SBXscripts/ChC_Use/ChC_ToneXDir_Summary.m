close all
clear all
clc

ds = 'ToneXDir_ExptList';
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

data_dfof_dir_all = [];
data_dfof_ori_all = [];
OSI_all = [];
DSI_all = [];
h_dir_use = [];
cell_n = [];
driver_all = [];
pupil_resp_all = [];
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

    data_dfof_dir_all = cat(1,data_dfof_dir_all, data_dfof_dir);
    data_dfof_ori_all = cat(1,data_dfof_ori_all, data_dfof_ori);
    OSI_all = [OSI_all; stim_OSI'];
    DSI_all = [DSI_all; stim_DSI'];
    h_dir_use = cat(1,h_dir_use,permute(h_dir,[3 1 2]));
    driver_all = [driver_all; repmat(driver,size(stim_OSI'))];
    cell_n = [cell_n; size(h_dir,3) length(find(sum(sum(h_dir,1),2)))];
    
    load(fullfile(LG_base,datemouse,datemouserun,[datemouserun '_pupil.mat']))
    pupil_resp_all = cat(3,pupil_resp_all,pupil_resp);
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
cellN = zeros(size(cellTypes));
figure; movegui('center') 
for i = 1:length(cellTypes)
    cellType_ind{i} = find(strcmp(driver_all,cellTypes{i}));
    cellType_vis_ind{i} = intersect(find(sum(sum(h_dir_use,2),3)),find(strcmp(driver_all,cellTypes{i})));
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
[max_dir_val max_dir_ind] = max(mean(data_dfof_dir_all(:,1:16,:,1),3),[],2);
[max_ori_val max_ori_ind] = max(mean(data_dfof_ori_all(:,:,:,1),3),[],2);
for iCell = 1:size(OSI_all,1)
    shift = nDirs/2 - max_dir_ind(iCell);
    data_dfof_dir_shift(iCell,:,:) = circshift(data_dfof_dir_all(iCell,1:16,:,1),shift,2);
    ori_shift = nDirs/4 - max_ori_ind(iCell);
    data_dfof_ori_shift(iCell,:,:) = circshift(data_dfof_ori_all(iCell,:,:),ori_shift,2);
end


figure; movegui('center')
norm_val =squeeze(max(max(data_dfof_dir_shift,[],2),[],3));
for ii = 1:length(cellTypes)
    subplot(2,length(cellTypes),ii)
    
    errorbar(Dirs, squeeze(mean(data_dfof_dir_shift(cellType_vis_ind{ii},:,1)./norm_val(cellType_vis_ind{ii}),1)), ...
        squeeze(std(data_dfof_dir_shift(cellType_vis_ind{ii},:,1)./norm_val(cellType_vis_ind{ii}),[],1)./sqrt(length(cellType_vis_ind{ii}))))
    hold on
    errorbar(Dirs, squeeze(mean(data_dfof_dir_shift(cellType_vis_ind{ii},:,2)./norm_val(cellType_vis_ind{ii}),1)), ...
        squeeze(std(data_dfof_dir_shift(cellType_vis_ind{ii},:,2)./norm_val(cellType_vis_ind{ii}),[],1)./sqrt(length(cellType_vis_ind{ii}))))
    xlabel('Direction')
    ylabel('Norm. dF/F')
    ylim([-0.5 1.5])
end


norm_val =squeeze(max(max(data_dfof_ori_shift,[],2),[],3));
for ii = 1:length(cellTypes)
    subplot(2,length(cellTypes),ii+length(cellTypes))
    errorbar(Oris, squeeze(mean(data_dfof_ori_shift(cellType_vis_ind{ii},:,1)./norm_val(cellType_vis_ind{ii}),1)), ...
        squeeze(std(data_dfof_ori_shift(cellType_vis_ind{ii},:,1)./norm_val(cellType_vis_ind{ii}),[],1)./sqrt(length(cellType_vis_ind{ii}))))
    hold on
    errorbar(Oris, squeeze(mean(data_dfof_ori_shift(cellType_vis_ind{ii},:,2)./norm_val(cellType_vis_ind{ii}),1)), ...
        squeeze(std(data_dfof_ori_shift(cellType_vis_ind{ii},:,2)./norm_val(cellType_vis_ind{ii}),[],1)./sqrt(length(cellType_vis_ind{ii}))))
    xlabel('Orientation')
    ylabel('Norm. dF/F')
    ylim([-0.5 1.5])
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
errorbar(tt,mean(pupil_resp_all_sub(:,1,:),3),std(pupil_resp_all_sub(:,1,:),[],3)./sqrt(nexp))
hold on
errorbar(tt,mean(pupil_resp_all_sub(:,2,:),3),std(pupil_resp_all_sub(:,2,:),[],3)./sqrt(nexp))
xlabel('Time from stim')
xlim([-2000 2000])
print(fullfile(fnout,[sv_name '_pupilDiameterSummary.pdf']),'-dpdf','-fillpage')
