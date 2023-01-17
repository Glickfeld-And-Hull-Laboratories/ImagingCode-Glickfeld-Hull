close all
clear all
clc

ds = 'SFxDir_ExptList';
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
data_dfof_dir_sf_all = [];
data_dfof_ori_all = [];
prefSF_all = [];
OSI_all = [];
DSI_all = [];
h_dir_use = [];
driver_all = [];
r_wheel_all = [];
r_pupil_all = [];
pupil_mod_all = [];
cell_n = [];
data_dfof_tc_all = [];
r_all = [];
r_driver_all = [];
[n n2] = subplotn(nexp);
for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc;
    driver = expt(iexp).driver;
    ImgFolder = expt(iexp).sfFolder;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);
    
    datemouse = [date '_' mouse];
    datemouserun = [date '_' mouse '_' run_str];

    fprintf([datemouse '\n'])
    load(fullfile(LG_base,datemouse,datemouserun,[datemouserun '_TCs.mat']))

    r = corrcoef(smoothdata(npSub_tc,1));
    r_all = [r_all; triu2vec(r)];
    r_driver_all = [r_driver_all; repmat(driver,size(triu2vec(r)))];

    load(fullfile(LG_base,datemouse,datemouserun,[datemouserun '_stimData.mat']))
    data_dfof_tc_all = [data_dfof_tc_all squeeze(mean(data_dfof_tc,2))];

    load(fullfile(LG_base,datemouse,datemouserun,[datemouserun '_respData.mat']))
    load(fullfile(LG_base,datemouse,datemouserun,[datemouserun '_SFfits.mat']))

    data_dfof_dir_sf_all = cat(1,data_dfof_dir_sf_all,data_dfof_dir_sf);
    data_dfof_ori_all = cat(1,data_dfof_ori_all, data_dfof_ori);
    prefSF_all = [prefSF_all; prefSF'];
    OSI_all = [OSI_all; OSI];
    DSI_all = [DSI_all; DSI];
    h_dir_use = cat(1,h_dir_use,permute(h_dir,[3 1 2]));
    driver_all = [driver_all; repmat(driver,size(OSI))];
    cell_n = [cell_n; size(h_dir,3) length(h_dir_all)];

    load(fullfile(LG_base,datemouse,datemouserun,[datemouserun '_wheelSpeed.mat']))
    
    r_wheel_all = [r_wheel_all; r_wheel];
    
    load(fullfile(LG_base,datemouse,datemouserun,[datemouserun '_pupil.mat']))
    r_pupil_all = [r_pupil_all; r_pupil];
    
    if size(s_pos_vis_celltc,2) ~= size(OSI,1)
        break
    end
    %vis+ vis- blank+ blank-
    pupil_mod_all = cat(2,pupil_mod_all, cat(3,mean(s_pos_vis_celltc,3,'omitnan'), mean(s_neg_vis_celltc,3,'omitnan'), mean(s_pos_celltc,3,'omitnan'), mean(s_neg_celltc,3,'omitnan')));
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
writetable(cellTab,fullfile(fnout,'cellDataTable.txt'))

nCells = size(OSI_all,1);

load(fullfile(LG_base,datemouse,datemouserun,[datemouserun '_input.mat']))
Dir = celleqel2mat_padded(input.tGratingDirectionDeg); 
if find(Dir>+360)
    Dir(find(Dir>=360)) = Dir(find(Dir>=360))-360;
end
Dirs = unique(Dir);
nDirs = length(Dirs);
SF_mat = celleqel2mat_padded(input.tGratingSpatialFreqCPD);
SFs = unique(SF_mat);
nSF = length(SFs);

cellTypes = unique(driver_all);
cellType_ind = cell(size(cellTypes));
cellType_vis_ind = cell(size(cellTypes));
cellN = zeros(size(cellTypes));
figure; movegui('center') 
for i = 1:length(cellTypes)
    cellType_ind{i} = find(strcmp(driver_all,cellTypes{i}));
    cellType_vis_ind{i} = intersect(find(sum(sum(h_dir_use,2),3)),find(strcmp(driver_all,cellTypes{i})));
    cellN(i) = length(cellType_vis_ind{i});
    subplot(2,3,1)
    cdfplot(OSI_all(cellType_vis_ind{i}))
    hold on
    subplot(2,3,2)
    cdfplot(DSI_all(cellType_vis_ind{i}))
    hold on
    subplot(2,3,3)
    cdfplot(prefSF_all(cellType_vis_ind{i}))
    hold on
    subplot(2,3,4)
    cdfplot(r_wheel_all(cellType_vis_ind{i}))
    hold on
    subplot(2,3,5)
    cdfplot(r_pupil_all(cellType_vis_ind{i}))
    hold on
    subplot(2,3,6)
    ind = find(strcmp(r_driver_all, cellTypes{i}));
    cdfplot(r_all(ind))
    hold on
end
subplot(2,3,1)
title('OSI')
xlabel('OSI')
legend(cellTypes,'location','southeast')
subplot(2,3,2)
title('DSI')
xlabel('DSI')
subplot(2,3,3)
title('Pref SF')
xlabel('log2(Pref SF)')
subplot(2,3,4)
title('Running Corr')
xlabel('R')
subplot(2,3,5)
title('Pupil Corr')
xlabel('R')
subplot(2,3,6)
title('Cell-cell Corr')
xlabel('R')
print(fullfile(fnout,'TuningSummary.pdf'),'-dpdf','-fillpage')


data_dfof_dir_sf_shift = zeros(size(data_dfof_dir_sf_all));
data_dfof_dir_shift = zeros(nCells,nDirs);
data_dfof_ori_sf_shift = zeros(nCells,nDirs/2,nSF);
data_dfof_ori_shift = zeros(nCells,nDirs/2);
data_dfof_sf = zeros(nCells,nSF);
[max_dir_val max_dir_ind] = max(mean(data_dfof_dir_sf_all(:,:,:,1),3),[],2);
[max_ori_val max_ori_ind] = max(data_dfof_ori_all(:,:,1),[],2);
[max_sf_val max_sf_ind] = max(mean(data_dfof_dir_sf_all(:,:,:,1),2),[],3);
for iCell = 1:size(OSI_all,1)
    shift = nDirs/2 - max_dir_ind(iCell);
    data_dfof_dir_sf_shift(iCell,:,:,:) = circshift(data_dfof_dir_sf_all(iCell,:,:,:),shift,2);
    data_dfof_dir_shift(iCell,:) = data_dfof_dir_sf_shift(iCell,:,max_sf_ind(iCell),1);
    data_dfof_sf(iCell,:) = squeeze(data_dfof_dir_sf_shift(iCell,max_dir_ind(iCell),:,1))';
    ori_shift = nDirs/4 - max_ori_ind(iCell);
    data_dfof_ori_shift(iCell,:) = circshift(data_dfof_ori_all(iCell,:,1),ori_shift,2);
    data_dfof_ori_sf_shift(iCell,:,:) = circshift(squeeze(mean(reshape(squeeze(data_dfof_dir_sf_shift(iCell,:,:,1)),[nDirs/2 2 nSF]),2)),-2,1);
end

data_dfof_dir_sf_shift_rect = data_dfof_dir_sf_shift;
data_dfof_dir_sf_shift_rect(find(data_dfof_dir_sf_shift<0)) = 0;
data_dfof_ori_sf_shift_rect = data_dfof_ori_sf_shift;
data_dfof_ori_sf_shift_rect(find(data_dfof_ori_sf_shift<0)) = 0;
DSI_SF = squeeze((data_dfof_dir_sf_shift_rect(:,nDirs/2,:,1)-data_dfof_dir_sf_shift_rect(:,nDirs,:,1))./(data_dfof_dir_sf_shift_rect(:,nDirs/2,:,1)+data_dfof_dir_sf_shift_rect(:,nDirs,:,1)));
OSI_SF = squeeze((data_dfof_ori_sf_shift_rect(:,nDirs/4,:,1)-data_dfof_ori_sf_shift_rect(:,nDirs/2,:,1))./(data_dfof_ori_sf_shift_rect(:,nDirs/4,:,1)+data_dfof_ori_sf_shift_rect(:,nDirs/2,:,1)));
figure; movegui('center')
for i = 1:nSF
    for ii = 1:length(cellTypes)
        subplot(2,length(cellTypes),ii)
        errorbar(SFs,mean(DSI_SF(cellType_vis_ind{ii},:),1,'omitnan'),std(DSI_SF(cellType_vis_ind{ii},:),[],1,'omitnan')./sqrt(length(cellType_vis_ind{ii})),'ok')
        ylabel('DSI')
        xlabel('SF')
        title(num2str(cellTypes{ii}))
        ylim([-0.25 0.75])
        xlim([0 0.4])
        hold on
        subplot(2,length(cellTypes),ii+length(cellTypes))
        errorbar(SFs,mean(OSI_SF(cellType_vis_ind{ii},:),1,'omitnan'),std(OSI_SF(cellType_vis_ind{ii},:),[],1,'omitnan')./sqrt(length(cellType_vis_ind{ii})),'ok')
        ylabel('OSI')
        xlabel('SF')
        ylim([-0.25 0.75])
        xlim([0 0.4])
        hold on
    end
end
print(fullfile(fnout,'OSI&DSIbySF.pdf'),'-dpdf','-bestfit')

figure; movegui('center')
for i = 1:nSF
    for ii = 1:length(cellTypes)
        subplot(2,length(cellTypes),ii)
        cdfplot(DSI_SF(cellType_vis_ind{ii},i))
        xlabel('DSI')
        title(num2str(cellTypes{ii}))
        hold on
        subplot(2,length(cellTypes),ii+length(cellTypes))
        cdfplot(OSI_SF(cellType_vis_ind{ii},i))
        xlabel('OSI')
        title(num2str(cellTypes{ii}))
        hold on
    end
end



norm_val =squeeze(max(max(data_dfof_dir_sf_shift(:,:,:,1),[],2),[],3));
figure; movegui('center')
for i = 1:nSF
    for ii = 1:length(cellTypes)
        subplot(2,length(cellTypes),i)
        errorbar(Dirs, squeeze(mean(data_dfof_dir_sf_shift(cellType_vis_ind{ii},:,i,1)./norm_val(cellType_vis_ind{ii}),1)), ...
            squeeze(std(data_dfof_dir_sf_shift(cellType_vis_ind{ii},:,i,1)./norm_val(cellType_vis_ind{ii}),[],1)./sqrt(length(cellType_vis_ind{ii}))))
        hold on
    end
    ylim([0 1])
    title(num2str(SFs(i)))
    xlabel('Direction')
    ylabel('dF/F')
end

figure; movegui('center')
for ii = 1:length(cellTypes)
    for i = 1:nSF
        subplot(1,length(cellTypes),ii)
        errorbar(Dirs, squeeze(mean(data_dfof_dir_sf_shift(cellType_vis_ind{ii},:,i,1)./norm_val(cellType_vis_ind{ii}),1)), ...
            squeeze(std(data_dfof_dir_sf_shift(cellType_vis_ind{ii},:,i,1)./norm_val(cellType_vis_ind{ii}),[],1)./sqrt(length(cellType_vis_ind{ii}))))
        hold on
    end
    if ii == 1
        legend(num2str(SFs'))
    end
    ylim([0 1])
    title([num2str(cellTypes{ii}) ' n = ' num2str(length(cellType_vis_ind{ii}))])
    xlabel('Direction')
    ylabel('dF/F')
end
print(fullfile(fnout,'DirectionBySFTuning.pdf'),'-dpdf','-bestfit')

[n n2] = subplotn(length(cellTypes)+1);
norm_ori_val = max(data_dfof_ori_shift,[],2);
figure; movegui('center')
for ii = 1:length(cellTypes)
    for i = 1:nSF
        subplot(n,n2,ii)
        errorbar(Dirs(1:nDirs/2), squeeze(mean(data_dfof_ori_sf_shift(cellType_vis_ind{ii},:,i)./norm_val(cellType_vis_ind{ii}),1)), ...
            squeeze(std(data_dfof_ori_sf_shift(cellType_vis_ind{ii},:,i)./norm_val(cellType_vis_ind{ii}),[],1)./sqrt(length(cellType_vis_ind{ii}))))
        hold on
    end
    ylim([0 1])
    title(num2str(cellTypes{ii}))
    xlabel('Orientation')
    ylabel('dF/F')
    subplot(n,n2,length(cellTypes)+1)
    errorbar(Dirs(1:nDirs/2), squeeze(mean(data_dfof_ori_shift(cellType_vis_ind{ii},:)./norm_ori_val(cellType_vis_ind{ii}),1)), ...
            squeeze(std(data_dfof_ori_shift(cellType_vis_ind{ii},:)./norm_ori_val(cellType_vis_ind{ii}),[],1)./sqrt(length(cellType_vis_ind{ii}))))
    hold on
    xlabel('Orientation')
    ylabel('dF/F')
    ylim([0 1])
end



figure; movegui('center') 
for i = 1:length(cellTypes)
    subplot(2,length(cellTypes),i) 
    imagesc(squeeze(mean(data_dfof_dir_sf_shift(cellType_vis_ind{i},:,:,1),1)))
    set(gca,'XtickLabels',SFs,'YTickLabels',Dirs(2:2:end))
    title([cellTypes{i} ' ' num2str(cellN(i))])
    subplot(2,length(cellTypes),i+length(cellTypes)) 
    imagesc(squeeze(mean(data_dfof_dir_sf_shift(cellType_vis_ind{i},:,:,1)./max(max(data_dfof_dir_sf_shift(cellType_vis_ind{i},:,:,1),[],2),[],3),1)))
    set(gca,'XtickLabels',SFs,'YTickLabels',Dirs(2:2:end))
end
print(fullfile(fnout,'DirectionBySFHeatmap.pdf'),'-dpdf','-bestfit')



figure; movegui('center') 
for i = 1:length(cellTypes)
    subplot(2,length(cellTypes),i) 
    imagesc(squeeze(mean(data_dfof_ori_sf_shift(cellType_vis_ind{i},:,:),1)))
    set(gca,'XtickLabels',SFs,'YTickLabels',Dirs(1:nDirs/2))
    title([cellTypes{i} ' ' num2str(cellN(i))])
    subplot(2,length(cellTypes),i+length(cellTypes)) 
    imagesc(squeeze(mean(data_dfof_ori_sf_shift(cellType_vis_ind{i},:,:)./max(max(data_dfof_ori_sf_shift(cellType_vis_ind{i},:,:),[],2),[],3),1)))
    set(gca,'XtickLabels',SFs,'YTickLabels',Dirs(1:nDirs/2))
end

figure; movegui('center') 
data_dfof_sf_avg = squeeze(mean(data_dfof_dir_sf_shift(:,:,:,1),2));
data_dfof_dir_avg = squeeze(mean(data_dfof_dir_sf_shift(:,:,:,1),3));
for i = 1:length(cellTypes)
    subplot(2,2,1)
    errorbar(Dirs, mean(data_dfof_dir_shift(cellType_vis_ind{i},:)./max(data_dfof_dir_shift(cellType_vis_ind{i},:),[],2),1), std(data_dfof_dir_shift(cellType_vis_ind{i},:)./max(data_dfof_dir_shift(cellType_vis_ind{i},:),[],2),[],1)./sqrt(length(cellType_vis_ind{i})))
    hold on
    subplot(2,2,2)
    errorbar(SFs, mean(data_dfof_sf(cellType_vis_ind{i},:)./max(data_dfof_sf(cellType_vis_ind{i},:),[],2),1), std(data_dfof_sf(cellType_vis_ind{i},:)./max(data_dfof_sf(cellType_vis_ind{i},:),[],2),[],1)./sqrt(length(cellType_vis_ind{i})))
    hold on
    subplot(2,2,3)
    errorbar(Dirs, mean(data_dfof_dir_avg(cellType_vis_ind{i},:)./max(data_dfof_dir_avg(cellType_vis_ind{i},:),[],2),1), std(data_dfof_dir_avg(cellType_vis_ind{i},:)./max(data_dfof_dir_avg(cellType_vis_ind{i},:),[],2),[],1)./sqrt(length(cellType_vis_ind{i})))
    ylim([0 1])
    hold on
    subplot(2,2,4)
    errorbar(SFs, mean(data_dfof_sf_avg(cellType_vis_ind{i},:)./max(data_dfof_sf_avg(cellType_vis_ind{i},:),[],2),1), std(data_dfof_sf_avg(cellType_vis_ind{i},:)./max(data_dfof_sf_avg(cellType_vis_ind{i},:),[],2),[],1)./sqrt(length(cellType_vis_ind{i})))
    hold on
end
subplot(2,2,1)
ylim([0 1])
title('Pref SF')
xlabel('Direction')
subplot(2,2,2)
title('Pref Dir')
xlabel('SFs')
subplot(2,2,3)
ylim([0 1])
title('All SF')
xlabel('Direction')
subplot(2,2,4)
title('All Dir')
xlabel('SFs')

pupil_mod_all_avg = squeeze(mean(pupil_mod_all(40:50,:,:),1)-mean(pupil_mod_all(1:10,:,:),1));
mod_str = {'vis+', 'vis-', 'blank+', 'blank-'};
figure; movegui('center') 
for i = 1:length(cellTypes)
    for ii = 1:4
        subplot(2,2,ii)
        cdfplot(pupil_mod_all_avg(cellType_ind{i},ii))
        hold on
        title(mod_str{ii})
    end
end
print(fullfile(fnout,'EyeMovementSummary.pdf'),'-dpdf','-fillpage')

figure;
tt = (1-nOff/2:nOn+nOff/2).*1000/params.frameRate;
for i = 1:length(cellTypes)
    subplot(2,2,i)
    plot(tt,data_dfof_tc_all(:,cellType_ind{i}),'c')
    hold on
    plot(tt,mean(data_dfof_tc_all(:,cellType_ind{i}),2),'k')
    title(cellTypes{i})
    xlabel('Time (s)')
    ylabel('dF/F')
    if i>2
        ylim([-0.5 1.5])
    else
       ylim([-0.1 0.4]) 
    end
end
print(fullfile(fnout,'AllTimeCourses.pdf'),'-dpdf','-fillpage')

figure;
for i = 1:length(cellTypes)
    subplot(2,2,i)
    histogram(mean(data_dfof_tc_all(nOff/2+1:nOff/2+1+nOn,cellType_ind{i}),1),[-0.2:0.1:1.5])
    xlabel('dF/F')
    ylabel('Number of cells')
    title(cellTypes{i})
end
print(fullfile(fnout,'AllResponses.pdf'),'-dpdf','-fillpage')

figure;
for i = 1:length(cellTypes)
    subplot(2,2,1)
    plot(tt,mean(data_dfof_tc_all(:,cellType_ind{i}),2))
    hold on
    if i == length(cellTypes)
        xlabel('Time (s)')
        ylabel('dF/F')
        ylim([-0.1 0.4])
        legend([cellTypes],'location','northwest')
        title('All cells')
    end
    subplot(2,2,2)
    plot(tt,mean(data_dfof_tc_all(:,cellType_ind{i}),2)./max(mean(data_dfof_tc_all(:,cellType_ind{i}),2),[],1))
    hold on
    if i == length(cellTypes)
        xlabel('Time (s)')
        ylabel('Normalized dF/F')
        ylim([-0.1 1.1])
        title('All cells')
    end
    subplot(2,2,3)
    plot(tt,mean(data_dfof_tc_all(:,cellType_vis_ind{i}),2))
    hold on
    if i == length(cellTypes)
        xlabel('Time (s)')
        ylabel('dF/F')
        ylim([-0.1 0.4])
        legend([cellTypes],'location','northwest')
        title('Responsive')
    end
    subplot(2,2,4)
    plot(tt,mean(data_dfof_tc_all(:,cellType_vis_ind{i}),2)./max(mean(data_dfof_tc_all(:,cellType_vis_ind{i}),2),[],1))
    hold on
    if i == length(cellTypes)
        xlabel('Time (s)')
        ylabel('Normalized dF/F')
        ylim([-0.1 1.1])
        title('Responsive')
    end
end
print(fullfile(fnout,'TimeCourseSummary.pdf'),'-dpdf','-fillpage')
