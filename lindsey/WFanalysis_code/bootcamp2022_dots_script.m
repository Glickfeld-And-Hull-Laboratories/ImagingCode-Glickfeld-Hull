%% creates ROIs from retinotopy experiment- only needs to be done once for each mouse
%run iXXX_paths.m (for your mouse) before running code. 
%will also need to manually choose number of ROIs (is hard coded below) and
%assign them to areas 'area_list'
close all
clear all
clc

behav_folder = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
data_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Data\Widefield\';
date = '220830';
mouse = 'i1373';
run = 'Run_';
suffix = '_MMStack_Pos0.ome';
time = {'1839','1901'};
rate = 10;
%% Read files and load image stack
roi_data = [];
for i = 1:length(time)
    cd(fullfile(data_pn, [date '_' mouse]));
    imageStack = [mouse '_' run num2str(i) suffix];
    roi_data_temp = double(readtiff([imageStack '.tif']));
    % roi_data_avg = mean(roi_data,3);
    roi_data = cat(3,roi_data,roi_data_temp);
    load(fullfile(behav_folder,['data-' mouse '-' date '-' time{i}]));
    if i > 1
        ntrials = length(input.counterValues);
        for ii = 1:ntrials
            input.counterValues{ii} = input.counterValues{ii}+count;
            input.counterTimesUs{ii} = input.counterTimesUs{ii}+count_base;
            input.wheelSpeedTimesUs{ii} = input.wheelSpeedTimesUs{ii}+wheel_base;
        end
    end
    count = input.counterValues{end}(end);
    count_base = input.counterTimesUs{end}(end);
    wheel_base = input.wheelSpeedTimesUs{end}(end);
    temp(i) = input;
end
clear roi_data_temp
input = concatenateDataBlocks(temp);

figure; 
imagesc(mean(roi_data,3))
colormap gray
movegui('center')
%% Calculate DF/F for each direction in one cycle based on mean of the offs for that direction

nOns = input.nScansOn; %identify how many frames the stimulus was presented for
nOffs = input.nScansOff; % %identify how many frames the stimulus was off
nStim = 3; %input.gratingContrastStepN; %identify how many stimuli were shown
nRep = size(input.counterValues,2)./nStim; %identify the number of times the stimuli were repeated
nTrials = size(input.counterValues,2); %identify the total number of trials

nTotalFrames = (nOns + nOffs)*nStim*nRep; %confirm total number of frames collected
siz = size(roi_data); %dertermine the XY size of the image frames
 
%make a single stim average movie
%pre-allocate data structures in the correct dimentions, defined by the
%size of the XY size, number of frames in each trial, and number of trials
trial_data = zeros(siz(1), siz(2), nOns+nOffs, nTrials);
trial_dF = zeros(siz(1), siz(2), nOns+nOffs, nTrials);
trial_dFoverF = zeros(siz(1), siz(2), nOns+nOffs, nTrials);
trial_F = zeros(siz(1), siz(2), nTrials);
%loop through trials
for i = 1:nTrials
    trial_data(:,:,:,i) = roi_data(:,:,1+((i-1)*(nOns+nOffs)):i*(nOns+nOffs));
    trial_F(:,:,i) = mean(trial_data(:,:,1+nOffs/2:nOffs,i),3);
    trial_dF(:,:,:,i) = bsxfun(@minus, trial_data(:,:,:,i), trial_F(:,:,i));
    trial_dFoverF(:,:,:,i) = bsxfun(@rdivide, trial_dF(:,:,:,i), trial_F(:,:,i));
end
avg_dFoverF = mean(trial_dFoverF,4);
avg_resp = mean(avg_dFoverF(:,:,nOffs+2:nOffs+10),3);

figure;
x = nOffs/2:nOffs+nOns/2;
[n n2] = subplotn(length(x));
for i = 1:length(x)
    subplot(n,n2,i)
    imagesc(avg_dFoverF(:,:,x(i)))
    colormap gray
    title(num2str(x(i)))
end
movegui('center')

%% Subtract vasculature from DF image, select and store ROI
figure; imagesc(avg_resp); colormap(gray);
movegui('center')
clim([0 .08])
%adjust number of ROIs to choose
%%
nROI = 6;
for i = 1:nROI
    roi(i) = impoly;
end
 
siz = size(avg_resp);
mask = zeros(siz(1),siz(2),nROI);
for i = 1:nROI
    mask(:,:,i) = createMask(roi(i));
end
 
roi_cluster = sum(mask,3);
mask_cell = bwlabel(roi_cluster);
figure; imagesc(mask_cell);
 movegui('center')

%%
%adjust list of areas to track
area_list = strvcat('V1', 'LM', 'A1', 'A2','PM', 'RL');
% ind = find(vasc_mask); 
mask_cell_V = mask_cell;
% mask_cell_V(ind) = 0;
% figure; imagesc(mask_cell_V);
 
%neuropil mask
np = imCellNeuropil(mask_cell, 2, 5);
np_V = zeros(size(np));
for i = 1:nROI
    np_V_temp = np(:,:,i);
%     np_V_temp(ind) = 0;
    np_V(:,:,i) = np_V_temp;
end
figure; imagesc(max(np_V,[],3))
movegui('center')

% save(fullfile(anal_pn, mouse, [mouse '_' roi_date '_roi_masks.mat']), 'roi_cluster', 'mask_cell', 'area_list', 'mask_cell_V', 'np_V');
 
%% get timecourses from all ROIS with retinotopy
roiTC = stackGetTimeCourses(roi_data, mask_cell_V);
nROI = size(roiTC,2);
roiTC_np = zeros(size(roiTC));
for i = 1:nROI
    roiTC_np(:,i) = stackGetTimeCourses(roi_data, np_V(:,:,i));
end

dir_mat = celleqel2mat_padded(input.tDotDirectionDeg);
dirs = unique(dir_mat);
nDir = length(dirs);
coh_mat = round(celleqel2mat_padded(input.tDotCoherence));
cohs = unique(coh_mat);
nCoh = length(cohs);
nTrials = length(dir_mat);

roiTC_trial = reshape(roiTC, [nOns+nOffs nTrials nROI]);
roiTC_np_trial = reshape(roiTC_np, [nOns+nOffs nTrials nROI]);

ind_use = cell(1,4);
roiTC_stim = cell(1,4);
roiTC_F_stim = cell(1,4);
roiTC_DFoverF_stim = cell(1,4);
roiTC_np_stim = cell(1,4);
roiTC_np_F_stim = cell(1,4);
roiTC_np_DFoverF_stim = cell(1,4);
roiTC_npSub_DFoverF_stim = cell(1,4);
roiTC_npSub_DFoverF_avg_sem = zeros(nOns+nOffs,nROI,4,2);
start = 1;
for i = 1:nCoh
    ind_coh = find(coh_mat==cohs(i));
    if i == 1
        ind_use{start} = ind_coh;
        roiTC_stim{start} = roiTC_trial(:,ind_use{start},:);
        roiTC_F_stim{start} = mean(roiTC_trial(nOffs/2:nOffs,ind_use{start},:),1);
        roiTC_DFoverF_stim{start} = (roiTC_stim{start}-roiTC_F_stim{start})./roiTC_F_stim{start};
        roiTC_np_stim{start} = roiTC_np_trial(:,ind_use{start},:);
        roiTC_np_F_stim{start} = mean(roiTC_np_trial(nOffs/2:nOffs,ind_use{start},:),1);
        roiTC_np_DFoverF_stim{start} = (roiTC_np_stim{start}-roiTC_np_F_stim{start})./roiTC_np_F_stim{start};
        roiTC_npSub_DFoverF_stim{start} = roiTC_DFoverF_stim{start}-roiTC_np_DFoverF_stim{start};
        roiTC_npSub_DFoverF_avg_sem(:,:,start,1) = squeeze(mean(roiTC_npSub_DFoverF_stim{start},2));
        roiTC_npSub_DFoverF_avg_sem(:,:,start,2) = squeeze(std(roiTC_npSub_DFoverF_stim{start},[],2)./sqrt(length(ind_use{start})));
        start = start+1;
    else
        ind_use{start} = ind_coh;
        roiTC_stim{start} = roiTC_trial(:,ind_use{start},:);
        roiTC_F_stim{start} = mean(roiTC_trial(nOffs/2:nOffs,ind_use{start},:),1);
        roiTC_DFoverF_stim{start} = (roiTC_stim{start}-roiTC_F_stim{start})./roiTC_F_stim{start};
        roiTC_np_stim{start} = roiTC_np_trial(:,ind_use{start},:);
        roiTC_np_F_stim{start} = mean(roiTC_np_trial(nOffs/2:nOffs,ind_use{start},:),1);
        roiTC_np_DFoverF_stim{start} = (roiTC_np_stim{start}-roiTC_np_F_stim{start})./roiTC_np_F_stim{start};
        roiTC_npSub_DFoverF_stim{start} = roiTC_DFoverF_stim{start}-roiTC_np_DFoverF_stim{start};
        roiTC_npSub_DFoverF_avg_sem(:,:,start,1) = squeeze(mean(roiTC_npSub_DFoverF_stim{start},2));
        roiTC_npSub_DFoverF_avg_sem(:,:,start,2) = squeeze(std(roiTC_npSub_DFoverF_stim{start},[],2)./sqrt(length(ind_use{start})));
        start = start+1;
        for ii = 1:nDir
            ind_dir = find(dir_mat == dirs(ii));
            ind_use{start} = intersect(ind_coh,ind_dir);
            roiTC_stim{start} = roiTC_trial(:,ind_use{start},:);
            roiTC_F_stim{start} = mean(roiTC_trial(nOffs/2:nOffs,ind_use{start},:),1);
            roiTC_DFoverF_stim{start} = (roiTC_stim{start}-roiTC_F_stim{start})./roiTC_F_stim{start};
            roiTC_np_stim{start} = roiTC_np_trial(:,ind_use{start},:);
            roiTC_np_F_stim{start} = mean(roiTC_np_trial(nOffs/2:nOffs,ind_use{start},:),1);
            roiTC_np_DFoverF_stim{start} = (roiTC_np_stim{start}-roiTC_np_F_stim{start})./roiTC_np_F_stim{start};
            roiTC_npSub_DFoverF_stim{start} = roiTC_DFoverF_stim{start}-roiTC_np_DFoverF_stim{start};
            roiTC_npSub_DFoverF_avg_sem(:,:,start,1) = squeeze(mean(roiTC_npSub_DFoverF_stim{start},2));
            roiTC_npSub_DFoverF_avg_sem(:,:,start,2) = squeeze(std(roiTC_npSub_DFoverF_stim{start},[],2)./sqrt(length(ind_use{start})));
            start = start+1;
        end
    end
end

%%
% coh 0; coh 1 both dir; coh 1 dir 0; coh 1 dir 180 
tt= (1-nOffs:nOns)./(rate./1000);

col_mat = {'b','g','m','c'};
leg_str = {'coh 0'; 'coh 1'; 'dir 0'; 'dir 180'};
figure;
[n n2] = subplotn(nROI);
for i = 1:nROI
    subplot(n,n2,i)
    for ii = 1:2
        shadedErrorBar(tt, roiTC_npSub_DFoverF_avg_sem(:,i,ii,1), roiTC_npSub_DFoverF_avg_sem(:,i,ii,2),'lineprops',col_mat{ii});
        hold on
    end
    title(area_list(i,:))
    ylabel('dF/F')
    ylim([-0.01 0.05])
    vline(0)
end
movegui('center')
sgtitle('Coh 0 - blue; Coh 1- green')
figure
for i = 1:nROI
    subplot(n,n2,i)
    for ii = 3:4
        shadedErrorBar(tt, roiTC_npSub_DFoverF_avg_sem(:,i,ii,1), roiTC_npSub_DFoverF_avg_sem(:,i,ii,2),'lineprops',col_mat{ii});
        hold on
    end
    ylabel('dF/F')
    ylim([-0.01 0.05])
    vline(0)
end
sgtitle('Dir 0 - magenta; Dir 180- cyan')
movegui('center')
 
%% wheel speed
 
[wheel_speed] = wheelSpeedCalc(input,32,'orange'); %add purple wheel!
figure; plot(wheel_speed)
wheel_tc = zeros(nOns+nOffs, nTrials);
for iTrial = 1:nTrials
    wheel_tc(:,iTrial) = wheel_speed(1+((iTrial-1).*(nOns+nOffs)):iTrial.*(nOns+nOffs));
end
wheel_trial_avg = mean(wheel_tc(nOffs:nOns+nOffs,:),1);
figure; movegui('center')
plot(wheel_trial_avg)

RIx = abs(wheel_trial_avg)>0.5;

%%
roiTC_npSub_DFoverF_run_avg_sem = zeros(nOns+nOffs,nROI,4,2);
roiTC_npSub_DFoverF_norun_avg_sem = zeros(nOns+nOffs,nROI,4,2);
ind_run = cell(1,4);
ind_norun = cell(1,4);
for i = 1:4
    ind_run{i} = find(ismember(ind_use{i},find(RIx)));
    roiTC_npSub_DFoverF_run_avg_sem(:,:,i,1) = squeeze(mean(roiTC_npSub_DFoverF_stim{i}(:,ind_run{i},:),2));
    roiTC_npSub_DFoverF_run_avg_sem(:,:,i,2) = squeeze(std(roiTC_npSub_DFoverF_stim{i}(:,ind_run{i},:),[],2)./sqrt(length(ind_run{i})));
    ind_norun{i} = find(ismember(ind_use{i},find(~RIx)));
    roiTC_npSub_DFoverF_norun_avg_sem(:,:,i,1) = squeeze(mean(roiTC_npSub_DFoverF_stim{i}(:,ind_norun{i},:),2));
    roiTC_npSub_DFoverF_norun_avg_sem(:,:,i,2) = squeeze(std(roiTC_npSub_DFoverF_stim{i}(:,ind_norun{i},:),[],2)./sqrt(length(ind_norun{i})));
end

leg_str = {'Coh = 0', 'Coh = 1'};
figure;
for i = 1:nROI
    for ii = 1:2
        subplot(2,nROI,i+((ii-1).*nROI))
        shadedErrorBar(tt,roiTC_npSub_DFoverF_run_avg_sem(:,i,ii,1),roiTC_npSub_DFoverF_run_avg_sem(:,i,ii,2),'lineprops',{'r','markerfacecolor','r'});
        hold on
        shadedErrorBar(tt,roiTC_npSub_DFoverF_norun_avg_sem(:,i,ii,1),roiTC_npSub_DFoverF_norun_avg_sem(:,i,ii,2),'lineprops',{'k','markerfacecolor','k'});
        if ii == 1;
            title(area_list(i,:))
        end
        if i == 1;
            legend(leg_str(ii))
        end
        ylim([-0.01 0.1])
    end
end    
sgtitle('Red- running; black- stationary')
movegui('center')

leg_str = {'Dir = 0', 'Dir = 180'};
figure;
for i = 1:nROI
    for ii = 3:4
        subplot(2,nROI,i+((ii-3).*nROI))
        shadedErrorBar(tt,roiTC_npSub_DFoverF_run_avg_sem(:,i,ii,1),roiTC_npSub_DFoverF_run_avg_sem(:,i,ii,2),'lineprops',{'r','markerfacecolor','r'});
        hold on
        shadedErrorBar(tt,roiTC_npSub_DFoverF_norun_avg_sem(:,i,ii,1),roiTC_npSub_DFoverF_norun_avg_sem(:,i,ii,2),'lineprops',{'k','markerfacecolor','k'});
        if ii-2 == 1;
            title(area_list(i,:))
        end
        if i == 1;
            legend(leg_str(ii-2))
        end
        ylim([-0.01 0.1])
    end
end    
sgtitle('Red- running; black- stationary')
movegui('center')