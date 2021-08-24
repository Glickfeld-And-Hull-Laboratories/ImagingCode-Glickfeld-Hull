clear;
analysis_out = 'A:\home\carlo\mikeAnalysis\2P\';
bdata_source = 'A:\home\mike\Data\Behavior\';

removing_bad_trials = false; %bad_trials = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]; %added to remove shitty trials, turn off if not needed

RCExptListMike_Inter; s=2;

img_fn = [expt.date(1,:) '_img' expt.mouse(1,:) '\getTC_' num2str(expt.run(1,:)) '\'];
load([analysis_out, img_fn, expt.region{1,1} 'figDataUntrimmed_' num2str(expt.date(1,:)) '_img' num2str(expt.mouse(1,:)) '.mat'], 'nIC','tt','targetAlign_tc','targetAlign_events','targetAligndFoverF','ind_block2','ind_rew','frameRateHz');
n = 1;
subplot(4,2,n)
title([expt.mouse(1,:) ' | ' num2str(nIC) ' neurons'])
hold on;
shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_rew),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_rew),3),[],2)./sqrt(nIC), 'k');
xlabel('Time from cue')
ylabel('Avg dF/F')
n = n+2;
subplot(4,2,n)
shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_block2),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_block2),3),[],2)./sqrt(nIC), 'r');
xlabel('Time from cue')
ylabel('Avg dF/F')

img_fn = [expt.date(2,:) '_img' expt.mouse(2,:) '\getTC_' num2str(expt.run(2,:)) '\'];
load([analysis_out, img_fn, (expt.region{1,2}) 'figDataUntrimmed_' num2str(expt.date(2,:)) '_img' num2str(expt.mouse(2,:)) '.mat'], 'nIC','tt','targetAlign_tc','targetAlign_events','targetAligndFoverF','ind_block2','ind_rew','frameRateHz');
n = n-1;
subplot(4,2,n)
title([expt.mouse(2,:) ' | ' num2str(nIC) 'neurons'])
hold on;
shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_rew),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_rew),3),[],2)./sqrt(nIC), 'k');
xlabel('Time from cue')
ylabel('Avg dF/F')
n = n+2;
subplot(4,2,n)
shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_block2),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_block2),3),[],2)./sqrt(nIC), 'r');
xlabel('Time from cue')
ylabel('Avg dF/F')

img_fn = [expt.date(3,:) '_img' expt.mouse(3,:) '\getTC_' num2str(expt.run(3,:)) '\'];
load([analysis_out, img_fn, (expt.region{1,3}) 'figDataUntrimmed_' num2str(expt.date(3,:)) '_img' num2str(expt.mouse(3,:)) '.mat'], 'nIC','tt','targetAlign_tc','targetAlign_events','targetAligndFoverF','ind_block2','ind_rew','frameRateHz');
n = n+1;
subplot(4,2,n)
title([expt.mouse(3,:) ' | ' num2str(nIC) 'neurons'])
hold on;
shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_rew),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_rew),3),[],2)./sqrt(nIC), 'k');
xlabel('Time from cue')
ylabel('Avg dF/F')
n = n+2;
subplot(4,2,n)
shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_block2),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_block2),3),[],2)./sqrt(nIC), 'r');
xlabel('Time from cue')
ylabel('Avg dF/F')

img_fn = [expt.date(4,:) '_img' expt.mouse(4,:) '\getTC_' num2str(expt.run(4,:)) '\'];
load([analysis_out, img_fn, (expt.region{1,4}) 'figDataUntrimmed_' num2str(expt.date(4,:)) '_img' num2str(expt.mouse(4,:)) '.mat'], 'nIC','tt','targetAlign_tc','targetAlign_events','targetAligndFoverF','ind_block2','ind_rew','frameRateHz');
n = n-1;
subplot(4,2,n)
title([expt.mouse(4,:) ' | ' num2str(nIC) 'neurons'])
hold on;
shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_rew),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_rew),3),[],2)./sqrt(nIC), 'k');
xlabel('Time from cue')
ylabel('Avg dF/F')
n = n+2;
subplot(4,2,n)
shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_block2),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_block2),3),[],2)./sqrt(nIC), 'r');
xlabel('Time from cue')
ylabel('Avg dF/F')

sgtitle(['dF over F for ' expt.name ' dendrites'])
savefig(fullfile(analysis_out,expt.name, 'indivMouseDFOF.fig'))

close all

img_fn = [expt.date(1,:) '_img' expt.mouse(1,:) '\getTC_' num2str(expt.run(1,:)) '\'];
load([analysis_out, img_fn, expt.region{1,1} 'figDataUntrimmed_' num2str(expt.date(1,:)) '_img' num2str(expt.mouse(1,:)) '.mat'], 'nIC','tt','targetAlign_tc','targetAlign_events','targetAligndFoverF','ind_block2','ind_rew','frameRateHz');
n = 1;
subplot(4,2,n)
title([expt.mouse(1,:) ' | ' num2str(nIC) ' neurons'])
hold on;
shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew),3),2).*(1000./frameRateHz), nanstd(nanmean(targetAlign_events(:,:,ind_rew),3),[],2)./sqrt(nIC).*(1000./frameRateHz), 'k');
xlabel('Time from cue')
ylabel('Spike rate (Hz)')
%ylim([0 8])
n = n+2;
subplot(4,2,n)
shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_block2),3),2).*(1000./frameRateHz), nanstd(nanmean(targetAlign_events(:,:,ind_block2),3),[],2)./sqrt(nIC).*(1000./frameRateHz), 'r');
xlabel('Time from cue')
ylabel('Spike rate (Hz)')
%ylim([0 8])

img_fn = [expt.date(2,:) '_img' expt.mouse(2,:) '\getTC_' num2str(expt.run(2,:)) '\'];
load([analysis_out, img_fn, (expt.region{1,2}) 'figDataUntrimmed_' num2str(expt.date(2,:)) '_img' num2str(expt.mouse(2,:)) '.mat'], 'nIC','tt','targetAlign_tc','targetAlign_events','targetAlign_events','ind_block2','ind_rew','frameRateHz');
n = n-1;
subplot(4,2,n)
title([expt.mouse(2,:) ' | ' num2str(nIC) 'neurons'])
hold on;
shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew),3),2).*(1000./frameRateHz), nanstd(nanmean(targetAlign_events(:,:,ind_rew),3),[],2)./sqrt(nIC).*(1000./frameRateHz), 'k');
xlabel('Time from cue')
ylabel('Spike rate (Hz)')
%ylim([0 8])
n = n+2;
subplot(4,2,n)
shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_block2),3),2).*(1000./frameRateHz), nanstd(nanmean(targetAlign_events(:,:,ind_block2),3),[],2)./sqrt(nIC).*(1000./frameRateHz), 'r');
xlabel('Time from cue')
ylabel('Spike rate (Hz)')
%ylim([0 8])

img_fn = [expt.date(3,:) '_img' expt.mouse(3,:) '\getTC_' num2str(expt.run(3,:)) '\'];
load([analysis_out, img_fn, (expt.region{1,3}) 'figDataUntrimmed_' num2str(expt.date(3,:)) '_img' num2str(expt.mouse(3,:)) '.mat'], 'nIC','tt','targetAlign_tc','targetAlign_events','targetAlign_events','ind_block2','ind_rew','frameRateHz');
n = n+1;
subplot(4,2,n)
title([expt.mouse(3,:) ' | ' num2str(nIC) 'neurons'])
hold on;
shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew),3),2).*(1000./frameRateHz), nanstd(nanmean(targetAlign_events(:,:,ind_rew),3),[],2)./sqrt(nIC).*(1000./frameRateHz), 'k');
xlabel('Time from cue')
ylabel('Spike rate (Hz)')
%ylim([0 8])
n = n+2;
subplot(4,2,n)
shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_block2),3),2).*(1000./frameRateHz), nanstd(nanmean(targetAlign_events(:,:,ind_block2),3),[],2)./sqrt(nIC).*(1000./frameRateHz), 'r');
xlabel('Time from cue')
ylabel('Spike rate (Hz)')
%ylim([0 8])

img_fn = [expt.date(4,:) '_img' expt.mouse(4,:) '\getTC_' num2str(expt.run(4,:)) '\'];
load([analysis_out, img_fn, (expt.region{1,4}) 'figDataUntrimmed_' num2str(expt.date(4,:)) '_img' num2str(expt.mouse(4,:)) '.mat'], 'nIC','tt','targetAlign_tc','targetAlign_events','targetAlign_events','ind_block2','ind_rew','frameRateHz');
n = n-1;
subplot(4,2,n)
title([expt.mouse(4,:) ' | ' num2str(nIC) 'neurons'])
hold on;
shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew),3),2).*(1000./frameRateHz), nanstd(nanmean(targetAlign_events(:,:,ind_rew),3),[],2)./sqrt(nIC).*(1000./frameRateHz), 'k');
xlabel('Time from cue')
ylabel('Spike rate (Hz)')
%ylim([0 8])
n = n+2;
subplot(4,2,n)
shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_block2),3),2).*(1000./frameRateHz), nanstd(nanmean(targetAlign_events(:,:,ind_block2),3),[],2)./sqrt(nIC).*(1000./frameRateHz), 'r');
xlabel('Time from cue')
ylabel('Spike rate (Hz)')
%ylim([0 8])

sgtitle(['dF over F for ' expt.name ' dendrites'])
savefig(fullfile(analysis_out,expt.name, 'indivMouseSpike.fig'))

close all