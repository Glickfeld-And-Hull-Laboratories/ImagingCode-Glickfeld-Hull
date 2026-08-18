original_t = (0:81000-1) / 15.5;
new_t = (0:81000-1) / 15;
cellTCs_15hz = interp1(original_t, npSub_tc, new_t);
%%
np_tc = interp1(original_t, np_tc, new_t);
data_tc = interp1(original_t, data_tc, new_t);
%%
input.stimOffs_photodiode = round(input.stimOffs_photodiode * 15 / 15.5);
input.stimOns_photodiode = round(input.stimOns_photodiode * 15 / 15.5);
save('input.mat','input')
%%

data_tc_interp = nan(nOn+nOff,nCells,nTrials);
for itrial = 1:nTrials
  if ~isnan(stimOns(itrial)) & (stimOns_15hz(itrial)+nOn+nOff/2)<nFrames
    data_tc_interp(:,:,itrial) = cellTCs_15hz(stimOns_15hz(itrial)-nOff/2:stimOns_15hz(itrial)-1+nOn+nOff/2,:);
  end
end

data_f_trial_interp = nanmean(data_tc_interp(1:nOff/2,:,:),1);
data_dfof_trial_interp = bsxfun(@rdivide, bsxfun(@minus,data_tc_interp, data_f_trial_interp), data_f_trial_interp);
data_dfof_trial_interp = permute(data_dfof_trial_interp,[1 3 2]); %nFrames x nTrials x nCells

tc_cell_avrg_interp = nanmean(data_dfof_trial_interp(:,:,:),3);%average pver cells, one row per trial
tc_trial_avrg_interp = squeeze(nanmean(data_dfof_trial_interp(:,:,:),2));%average over trials, one row per cell
tc_cell_trial_avrg = nanmean(tc_cell_avrg_interp,2);%average over trials and cells

figure;
plot(tc_trial_avrg_interp, 'LineWidth',.005,'color',[.25 .25 .25]);
hold on;
plot(tc_cell_trial_avrg, 'LineWidth',2, 'color','k');
hold on;
vline(nOff/2,'g')
title('');
hold off
ylim([-.02 .18])
%% 

data_dfof = data_dfof(18:529,:);
mask_cell=mask_cell(18:529,:);
mask_cell_red = mask_cell_red(18:529,:);
mask_np=mask_np(18:529,:);
save(fullfile( 'mask_cell.mat'), 'data_dfof', 'mask_cell', 'mask_cell_red', 'mask_np','mask_label')
%%
data_avg=data_avg(18:529,:);
regImg=regImg(18:529,:);
save(fullfile('regOuts&Img.mat'),'outs','regImg','data_avg')
%%
redChImg=redChImg(18:529,:);
save(fullfile('redImage'),'redChImg')

