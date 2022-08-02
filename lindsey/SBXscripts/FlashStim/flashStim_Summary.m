clear all; close all; clc;
ds = 'FlashStim_ExptList';
eval(ds)
frame_rate = params.frameRate;
nexp = size(expt,2);
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'FlashStim');
if ~exist(summaryDir)
    mkdir(summaryDir)
end
arcIx = [];
respOriIx = [];
respDirIx = [];
respPairOriIx = [];
respPairDirIx = [];
totCells = 0;
data_dfof_ori_all = [];
data_dfof_dir_all = [];
y_fit_all = [];
prefOri_all = [];
Rsq_all = [];
theta_90_all = [];
mouseIx = [];
mouse_list = cell(1,nexp);

for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    ImgFolder = expt(iexp).preFolder;
    nrun = length(ImgFolder);
    pre_run_str = catRunName(ImgFolder, nrun);
    ImgFolder = expt(iexp).postFolder;
    nrun = length(ImgFolder);
    post_run_str = catRunName(ImgFolder, nrun);
    mouse_list{iexp} = mouse;

    fprintf([mouse ' ' date '\n'])

    preData = load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' pre_run_str], [date '_' mouse '_' pre_run_str '_oriResp.mat']));
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' pre_run_str], [date '_' mouse '_' pre_run_str '_OriBootstrap.mat']));
    postData = load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' post_run_str], [date '_' mouse '_' post_run_str '_oriResp.mat']));
    
    Oris = preData.Oris;
    Dirs = preData.Dirs;
    indOri = find(Oris == expt(iexp).pairDeg);
    indDir = find(Dirs == expt(iexp).pairDeg);
    targetOri = ceil(length(Oris)/2);
    targetDir = ceil(length(Dirs)/2);
    oriShift = targetOri-indOri;
    dirShift = targetDir-indDir;

    nCells = size(preData.data_dfof_ori,1);
    arcIx = [arcIx; strcmp(expt(iexp).gRNA,'Arc').*ones(nCells,1)];
    respPairOriIx = [respPairOriIx; preData.h_ori(indOri,:)'];
    respOriIx = [respOriIx; (sum(preData.h_ori,1)>0)'];
    data_dfof_ori_all = cat(1, data_dfof_ori_all, cat(3, circshift(preData.data_dfof_ori(:,:,1),oriShift,2), circshift(postData.data_dfof_ori(:,:,1),oriShift,2)));
    y_fit_all = cat(1, y_fit_all, cat(3, preData.y_fit', postData.y_fit'));
    prefOri_all = [prefOri_all; [preData.u1_ori; postData.u1_ori]'];
    Rsq_all = [Rsq_all; [preData.R_square_ori; postData.R_square_ori]'];
    respPairDirIx = [respPairDirIx; preData.h_dir(indDir,:)'];
    respDirIx = [respDirIx; (sum(preData.h_dir,1)>0)'];
    data_dfof_dir_all = cat(1, data_dfof_dir_all, cat(3, circshift(preData.data_dfof_dir(:,:,1),dirShift,2), circshift(postData.data_dfof_dir(:,:,1),dirShift,2)));
    mouseIx = [mouseIx; iexp.*ones(nCells,1)];
    theta_90_all = [theta_90_all; theta_90];
end

%%

resp_ind_ori = intersect(find(std(y_fit_all(:,:,1),[],2)),intersect(find(theta_90_all<22.5),find(respOriIx)));
resp_ind_dir = intersect(find(std(y_fit_all(:,:,1),[],2)),intersect(find(theta_90_all<22.5),find(respDirIx)));
% resp_ind_ori = find(respPairOriIx);
% resp_ind_dir = find(respPairDirIx);
% figure; 
% errorbar(Oris', mean(data_dfof_ori_all(resp_ind,:,1),1), std(data_dfof_ori_all(resp_ind,:,1),[],1)./sqrt(length(resp_ind)),'-o')
% hold on
% errorbar(Oris', mean(data_dfof_ori_all(resp_ind,:,2),1), std(data_dfof_ori_all(resp_ind,:,2),[],1)./sqrt(length(resp_ind)),'-o')
% 
% resp_ind_use = intersect(resp_ind, find(sum(Rsq_all>0.5,2)==2));
% figure;
% histogram(rad2deg(prefOri_all(resp_ind_use,1)),Oris)
% hold on
% histogram(rad2deg(prefOri_all(resp_ind_use,2)),Oris)

data_dfof_ori_all_rect = data_dfof_ori_all;
data_dfof_ori_all_rect(find(data_dfof_ori_all_rect<0)) = 0;
tuningOri = data_dfof_ori_all_rect./sum(data_dfof_ori_all_rect,2);
deltaTuningOri = tuningOri(:,:,2)-tuningOri(:,:,1);

figure;
subplot(2,2,1)
for i = 1:nexp
    resp_ind_use = intersect(resp_ind_ori, find(mouseIx == i));
    errorbar(Oris-Oris(targetOri), mean(deltaTuningOri(resp_ind_use,:),1), std(deltaTuningOri(resp_ind_use,:),[],1)./sqrt(size(deltaTuningOri(resp_ind_use,:),1)),'-o')
    hold on
end
ylim([-0.1 0.1])
xlim([-90 110])
ylabel('Change in tuning')
xlabel('Orientation (deg)')
legend(mouse_list)

subplot(2,2,2)
n = zeros(1,2);
for i = 1:2
    resp_ind_use = intersect(resp_ind_ori, find(arcIx == i-1));
    errorbar(Oris-Oris(targetOri), mean(deltaTuningOri(resp_ind_use,:),1), std(deltaTuningOri(resp_ind_use,:),[],1)./sqrt(size(deltaTuningOri(resp_ind_use,:),1)),'-o')
    hold on
    n(i) = length(resp_ind_use);
end
ylim([-0.1 0.1])
xlim([-90 110])
ylabel('Change in tuning')
xlabel('Orientation (deg)')
legend({['LacZ- ' num2str(n(1))],['Arc- ' num2str(n(2))]})

data_dfof_dir_all_rect = data_dfof_dir_all;
data_dfof_dir_all_rect(find(data_dfof_dir_all_rect<0)) = 0;
tuningDir = data_dfof_dir_all_rect./sum(data_dfof_dir_all_rect,2);
deltaTuningDir = tuningDir(:,:,2)-tuningDir(:,:,1);

subplot(2,2,3)
for i = 1:nexp
    resp_ind_use = intersect(resp_ind_dir, find(mouseIx == i));
    errorbar(Dirs-Dirs(targetDir), mean(deltaTuningDir(resp_ind_use,:),1), std(deltaTuningDir(resp_ind_use,:),[],1)./sqrt(size(deltaTuningDir(resp_ind_use,:),1)),'-o')
    hold on
end
ylim([-0.1 0.1])
xlim([-180 200])
ylabel('Change in tuning')
xlabel('Direction (deg)')
legend(mouse_list)

subplot(2,2,4)
n = zeros(1,2);
for i = 1:2
    resp_ind_use = intersect(resp_ind_dir, find(arcIx == i-1));
    errorbar(Dirs-Dirs(targetDir), mean(deltaTuningDir(resp_ind_use,:),1), std(deltaTuningDir(resp_ind_use,:),[],1)./sqrt(size(deltaTuningDir(resp_ind_use,:),1)),'-o')
    hold on
    n(i) = length(resp_ind_use);
end
ylim([-0.1 0.1])
xlim([-180 200])
ylabel('Change in tuning')
xlabel('Direction (deg)')
legend({['LacZ- ' num2str(n(1))],['Arc- ' num2str(n(2))]})

sgtitle('Pre Pair Responsive Cells')
print(fullfile(summaryDir,'flashStimSummary_PrePairRespIx.pdf'),'-dpdf','-bestfit')

% resp_ind_use = intersect(resp_ind, find(sum(Rsq_all>0.5,2)==2));
% y_fit_all_norm = y_fit_all./max(y_fit_all,[],2);
% figure; 
% plot(y_fit_all_norm(resp_ind_use,:,1)','k')
% hold on
% plot(y_fit_all_norm(resp_ind_use,:,2)','c')
% 
% figure;
% prefOri_all_rect = rad2deg(prefOri_all);
% prefOri_all_rect(find(prefOri_all_rect<0)) = prefOri_all_rect(find(prefOri_all_rect<0))+180;
% prefOri_all_rect(find(prefOri_all_rect>90)) = 180-prefOri_all_rect(find(prefOri_all_rect>90));
% 
% scatter(prefOri_all_rect(resp_ind_ori,1),prefOri_all_rect(resp_ind_ori,2))
