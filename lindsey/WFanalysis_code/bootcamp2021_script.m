%% creates ROIs from retinotopy experiment- only needs to be done once for each mouse
%run iXXX_paths.m (for your mouse) before running code. 
%will also need to manually choose number of ROIs (is hard coded below) and
%assign them to areas 'area_list'
 
%% Read files and load image stack
irun = 1;
cd(fullfile(data_pn, mouse, [date '_' mouse], [mouse '_' run{irun}]));
imageStack = [mouse '_' run{irun} suffix];
roi_data = double(readtiff([imageStack '.tif']));
% roi_data_avg = mean(roi_data,3);

load(fullfile(behav_folder,['data-' mouse '-' date '-' time{irun}]));

figure; 
imagesc(mean(roi_data,3))
colormap gray
movegui('center')
%% Calculate DF/F for each direction in one cycle based on mean of the offs for that direction

nOns = input.nScansOn; %identify how many frames the stimulus was presented for
nOffs = input.nScansOff; % %identify how many frames the stimulus was off
nStim = 1; %input.gratingDirectionStepN; %identify how many stimuli were shown
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
% avg_ves_resp = mean(avg_dFoverF(:,:,nOffs+10:nOffs+30),3);
% writetiff(avg_ves_resp, fullfile(anal_pn, mouse, [mouse '_' roi_date '_roi_avg_ves_resp.tif']));
% %% Convert vessel traces on imageJ to vessel map and subtract from ROI mask
% 
% vasc_map = readtiff(fullfile(anal_pn,mouse, [mouse '_' roi_date '_vessel_trace.tif']));
% figure; imagesc(vasc_map);
% vasc_mask = zeros(size(vasc_map)); vasc_mask(find(vasc_map ==0)) = 1;
% figure; imagesc(vasc_mask); colormap(gray)
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
clim([0 .05])
%adjust number of ROIs to choose
%%
nROI = 5;
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
area_list = strvcat('LM','AL','V1','PM','AM');
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
%roiTC_npsub = roiTC-roiTC_np;
 
 
roiTC_stim = zeros(nOns+nOffs, nROI, nRep, nStim);
roiTC_F_stim = zeros(1, nROI, nRep, nStim);
roiTC_DF_stim = zeros(nOns+nOffs, nROI, nRep, nStim);
roiTC_DFoverF_stim = zeros(nOns+nOffs, nROI, nRep, nStim);
roiTCnp_stim = zeros(nOns+nOffs, nROI, nRep, nStim);
roiTCnp_F_stim = zeros(1, nROI, nRep, nStim);
roiTCnp_DF_stim = zeros(nOns+nOffs, nROI, nRep, nStim);
roiTCnp_DFoverF_stim = zeros(nOns+nOffs, nROI, nRep, nStim);
for i = 1:nStim
    for ii = 1:nRep
        roiTC_stim(:,:,ii,i) = roiTC(1+((i-1)*(nOns+nOffs))+((ii-1)*(nStim*(nOns+nOffs))):(i*(nOns+nOffs))+((ii-1)*(nStim*(nOns+nOffs))),:);
        roiTC_F_stim(:,:,ii,i) = mean(roiTC_stim((nOffs/2):nOffs,:,ii,i));
        roiTC_DF_stim(:,:,ii,i) = bsxfun(@minus, roiTC_stim(:,:,ii,i), roiTC_F_stim(:,:,ii,i));
        roiTC_DFoverF_stim(:,:,ii,i) = bsxfun(@rdivide, roiTC_DF_stim(:,:,ii,i), roiTC_F_stim(:,:,ii,i));
        roiTCnp_stim(:,:,ii,i) = roiTC_np(1+((i-1)*(nOns+nOffs))+((ii-1)*(nStim*(nOns+nOffs))):(i*(nOns+nOffs))+((ii-1)*(nStim*(nOns+nOffs))),:);
        roiTCnp_F_stim(:,:,ii,i) = mean(roiTCnp_stim((nOffs/2):nOffs,:,ii,i));
        roiTCnp_DF_stim(:,:,ii,i) = bsxfun(@minus, roiTCnp_stim(:,:,ii,i), roiTCnp_F_stim(:,:,ii,i));
        roiTCnp_DFoverF_stim(:,:,ii,i) = bsxfun(@rdivide, roiTCnp_DF_stim(:,:,ii,i), roiTCnp_F_stim(:,:,ii,i));
    end
end
roiTC_DFoverF_avg = squeeze(mean(reshape(roiTC_DFoverF_stim, [nOns+nOffs, nROI, nRep*nStim]),3));
roiTC_DFoverF_sem = squeeze(std(reshape(roiTC_DFoverF_stim, [nOns+nOffs, nROI, nRep*nStim]),[],3)./sqrt(double(nRep*nStim)));
 
roiTC_DF_avg = squeeze(mean(reshape(roiTC_DF_stim, [nOns+nOffs, nROI, nRep*nStim]),3));
roiTC_DF_sem = squeeze(std(reshape(roiTC_DF_stim, [nOns+nOffs, nROI, nRep*nStim]),[],3)./sqrt(double(nRep*nStim)));
 
tt= (1-nOffs:nOns)./(rate./1000);
 
roiTCnpsub_DFoverF_stim = roiTC_DFoverF_stim-roiTCnp_DFoverF_stim;
roiTCnpsub_DFoverF_avg = squeeze(mean(reshape(roiTCnpsub_DFoverF_stim, [nOns+nOffs, nROI, nRep*nStim]),3));
roiTCnpsub_DFoverF_sem = squeeze(std(reshape(roiTCnpsub_DFoverF_stim, [nOns+nOffs, nROI, nRep*nStim]),[],3)./sqrt(double(nRep*nStim)));
figure;
[n n2] = subplotn(nROI);
for i = 1:nROI
    subplot(n,n2,i)
    shadedErrorBar(tt, roiTCnpsub_DFoverF_avg(:,i), roiTCnpsub_DFoverF_sem(:,i));
    hold on
    vline(0)
    title(area_list(i,:))
    ylabel('dF/F')
end
movegui('center')

figure;
plot(tt, roiTCnpsub_DFoverF_avg)
ylim([-.01 0.05])
ylabel('dF/F')
hold on
vline(0)
legend(area_list, 'Location', 'BestOutside')
movegui('center')

%print(fullfile(anal_pn, mouse, [mouse '_' roi_date '_allAreaRetResp.pdf']),'-dpdf')
 
%% wheel speed
 
[wheel_speed] = wheelSpeedCalc(input,32,'red'); %add purple wheel!
figure; plot(wheel_speed)
wheel_tc = zeros(nOns+nOffs, nTrials);
for iTrial = 1:nTrials
    wheel_tc(:,iTrial) = wheel_speed(1+((iTrial-1).*(nOns+nOffs)):iTrial.*(nOns+nOffs));
end
wheel_trial_avg = mean(wheel_tc(nOffs:nOns+nOffs,:),1);
figure; movegui('center')
plot(wheel_trial_avg)

RIx = wheel_trial_avg>0.5;
ind1 = find(wheel_trial_avg<=0);
ind2 = find(wheel_trial_avg>0 & wheel_trial_avg<0.2);
ind3 = find(wheel_trial_avg>=0.2);

%%
roiTCnpsub_DFoverF_run_avg = mean(roiTCnpsub_DFoverF_stim(:,:,find(RIx),:),3);
roiTCnpsub_DFoverF_norun_avg = mean(roiTCnpsub_DFoverF_stim(:,:,find(RIx==0),:),3);
roiTCnpsub_DFoverF_run_sem = std(roiTCnpsub_DFoverF_stim(:,:,find(RIx),:),[],3)./sqrt(length(find(RIx)));
roiTCnpsub_DFoverF_norun_sem = std(roiTCnpsub_DFoverF_stim(:,:,find(RIx==0),:),[],3)./sqrt(length(find(RIx==0)));

figure;
for i = 1:nROI
    subplot(n,n2,i)
    shadedErrorBar(tt,roiTCnpsub_DFoverF_run_avg(:,i),roiTCnpsub_DFoverF_run_sem(:,i),{'r','markerfacecolor','r'});
    hold on
    shadedErrorBar(tt,roiTCnpsub_DFoverF_norun_avg(:,i),roiTCnpsub_DFoverF_norun_sem(:,i),{'k','markerfacecolor','k'});
    title(area_list(i,:))
end    
movegui('center')

%%
roiTCnpsub_DFoverF_runsome_avg = mean(roiTCnpsub_DFoverF_stim(:,:,ind2,:),3);
roiTCnpsub_DFoverF_norun_avg = mean(roiTCnpsub_DFoverF_stim(:,:,ind1,:),3);
roiTCnpsub_DFoverF_runsome_sem = std(roiTCnpsub_DFoverF_stim(:,:,ind2,:),[],3)./sqrt(length(ind2));
roiTCnpsub_DFoverF_norun_sem = std(roiTCnpsub_DFoverF_stim(:,:,ind1,:),[],3)./sqrt(length(ind1));
roiTCnpsub_DFoverF_runmore_avg = mean(roiTCnpsub_DFoverF_stim(:,:,ind3,:),3);
roiTCnpsub_DFoverF_runmore_sem = std(roiTCnpsub_DFoverF_stim(:,:,ind3,:),[],3)./sqrt(length(ind3));

figure;
for i = 1:nROI
    subplot(n,n2,i)
    shadedErrorBar(tt,roiTCnpsub_DFoverF_norun_avg(:,i),roiTCnpsub_DFoverF_norun_sem(:,i),{'k','markerfacecolor','k'});
    hold on
    shadedErrorBar(tt,roiTCnpsub_DFoverF_runsome_avg(:,i),roiTCnpsub_DFoverF_runsome_sem(:,i),{'r','markerfacecolor','r'});
    shadedErrorBar(tt,roiTCnpsub_DFoverF_runmore_avg(:,i),roiTCnpsub_DFoverF_runmore_sem(:,i),{'m','markerfacecolor','m'});
    title(area_list(i,:))
end    
movegui('center')

%%
roiTCnpsub_DFoverF_avg = squeeze(mean(roiTCnpsub_DFoverF_stim(nOffs+1:nOffs+nOns,:,:),1));
roiTCnpsub_DFoverF_avg2 = squeeze(mean(roiTCnpsub_DFoverF_stim(nOffs+1:nOffs+nOns/2,:,:),1));
[n n2] = subplotn(nROI);
r = zeros(1,nROI);
figure;
for i = 1:nROI
    subplot(n,n2,i)
    scatter(roiTCnpsub_DFoverF_avg(i,:),wheel_trial_avg);
    hold on
    scatter(roiTCnpsub_DFoverF_avg2(i,:),wheel_trial_avg);
    r(i) = triu2vec(corrcoef(roiTCnpsub_DFoverF_avg(i,:),wheel_trial_avg));
    xlabel('dF/F')
    ylabel('Speed')
    title([area_list(i,:) '- r = ' num2str(chop(r(i),2))])
    xlim([-0.01 0.03])
end
movegui('center')
roiTC_all = reshape(roiTCnpsub_DFoverF_avg, [nROI*nTrials 1]);
x = 1:nROI;
area_num = repmat(x', [1 nTrials]);
area_num = reshape(area_num,[nTrials*nROI 1]);
run_num = repmat(RIx,[nROI 1]);
run_num = reshape(run_num,[nTrials*nROI 1]);
[p t s] = anovan(roiTC_all,{area_num,run_num});

save([anal_pn '\' mouse '_exptData.mat'], 'roiTCnpsub_DFoverF_stim','RIx','area_list','nROI','tt','input');