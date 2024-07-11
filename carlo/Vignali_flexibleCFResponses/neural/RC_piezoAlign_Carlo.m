clear; close all;
%% load data and check if piezo is working that day
analysis_out = 'A:\home\ned\analysis\2P\';
bxSourceBase = 'A:\home\ned\rawData\behavioral\';

%RCExptListCarlo_Block
%RCExptListCarlo_Inter
expType=input('monomodal(1) or dimodal(2) session: ');
seshType=input('interleaved (1), blocked (2), or generalization (3) session: ');
if seshType==1 && expType == 1
RCExptListCarlo_Inter;
elseif seshType==2 && expType == 1
RCExptListCarlo_Block;
elseif seshType==1 && expType == 2
RCExptListCarlo_VA_Inter;
elseif seshType==2 && expType == 2
RCExptListCarlo_VA_Block;
elseif seshType==3 && expType == 2
RCExptListCarlo_VA_Inter;
elseif seshType==1 && expType == 3
RCExptListCarlo_CSPlus_Inter;
elseif seshType==2 && expType == 3
RCExptListCarlo_CSPlus_Block;
end
pairings = loadRCList_Ned(expt,seshType);

exp_subset = pairings{1,1}; 
iexp = exp_subset; isesh = pairings{3,1};
mouse = strtrim(expt{1,iexp}.mouse{1,isesh});
date = expt{1,iexp}.date{1,isesh};
run = expt{1,iexp}.run{1,isesh};

crpTrainDayList;
sessions = [date '_img' mouse];% for behavior data

% prewin_frames = 30;
% postwin_frames = 120;

cue = input('which session is this?: '); % 1: CS+ | 2: CS- | 3: interleaved
if cue == 1
    figTitle = 'CS+';
elseif cue == 2
    figTitle = 'CS-';
elseif cue == 3
    figTitle = 'Interleaved';
end

fprintf([date ' ' mouse '\n']);
threshold = input('Which threshold do you wish to plot? ("-#" for decon | "#.#" for first der.): ');
if isempty(threshold)
img_fn = [date '_img' mouse '\getTC_' run '\'];
elseif threshold<0 %deconvolve data
img_fn = [date '_img' mouse '\getTC_' run '_' num2str(threshold) '\'];
elseif threshold>0 %first derivative data
img_fn = [date '_img' mouse '\getTC_' run '_FD' num2str(threshold) '\'];        
end

img_fn2 = [mouse '_' date '_' run];
mworks = getBxData_Ned(bxSourceBase, sessions);
nf = mworks.counterValues{end}(end);

cd(fullfile('A:\home\ned\rawData\2P\', [date '_img' mouse]));
if ~exist(fullfile([img_fn2 '.ephys']))
cd(fullfile([analysis_out date '_img' mouse '\getTC_' run '\']));
load('_cueAlignPiezo.mat');
else
fn_piezo = fopen([img_fn2 '.ephys']);
piezo_data = fread(fn_piezo,'single');
piezo_data_volts = piezo_data(2:2:end);
% verify if data looks good that day
tempFig=setFigure; hold on;plot(piezo_data_volts);
savefig(fullfile(analysis_out,img_fn, '_piezoTrace.fig'));

%% downsample piezo data to each frame and plot piezo and licking
% for data collected before 2021, the frame numbers were not recorded. but it seems there is a certain number of extra piezo reads before the first frame.
% so use this way the determine the beginning and then every 3 piezo read is 1 frame
% piezo data is in 90Hz and imaging data is in 30Hz. 
if str2double(date) - 202100 < 0 % if data was collected before year 2021
    first_read = 227; % the first 113 reads are always frame 0, piezo_data is a vector that has both voltage and read numbers
    last_read = first_read + 6*nf - 1;
else
    first_read = find(piezo_data==1,1,'first');
    last_read = find(piezo_data==nf,1,'first')+5;
end

piezo_data_temp = piezo_data(first_read:last_read);
piezo_data_temp = reshape(piezo_data_temp,[2 size(piezo_data_temp,1)./2]);
if str2double(date) - 202100 < 0
    a = 1:1:nf;
    bxSourceBase = [a;a;a];
    b2 = bxSourceBase(:);
    piezo_data_temp (1,:) = b2';
end
piezo_frames = zeros(nf,1); % average piezo reads for each frame
for iframe = 1:nf
    ind = find(piezo_data_temp(1,:)==iframe);
    piezo_frames(iframe,:) = mean(piezo_data_temp(2,ind),2,'omitnan');
end

end
load(fullfile([analysis_out,img_fn,  '_targetAlign.mat']));

rewDelay_frames =  round(0.6.*frameRateHz);
cTargetOn = mworks.cTargetOn;
if iscell(cTargetOn) % if it is a cell, it means cTargetOn wasn't being over written in extract TC. If it's not a cell, it's already over written in extract TC. can be used directly
    cTargetOn = celleqel2mat_padded(mworks.cTargetOn);
    cTargetOn(1) = nan; % get rid of first trial
end
    load([analysis_out date '_img' mouse '\getTC_' run '\' date '_img' mouse '_' run 'saveOutputs.mat']);
    load([analysis_out,date '_img' mouse '\getTC_' run '\' date '_img' mouse '_' run, '_nPCA', ...
    num2str(nPCA), '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', ...
    num2str(cluster_threshold), '_TCave.mat'], 'cTargetOn_cutted');
    nIC = size(targetAlign_events,2); %you get targetalign_events from neural align, this is the neural data (in spikes) aligned to cue
nTrials = length(cTargetOn_cutted)+1;
if ~exist('targetAlign_piezo')
    loaded=0;
targetAlign_piezo = nan(prewin_frames+postwin_frames,nTrials);
for itrial = 1:nTrials
    if ~isnan(cTargetOn(itrial))
        if cTargetOn(itrial)+postwin_frames-1 <= nf
            targetAlign_piezo(:,itrial) = piezo_frames(cTargetOn(itrial)-prewin_frames:cTargetOn(itrial)+postwin_frames-1,:);
        end
    end
end
else
    loaded=1;
end
ind_nan = find(isnan(targetAlign_piezo(1,:)));

g= [0.4660 0.6740 0.1880];
r= [0.6350 0.0780 0.1840];


targetAlign_piezo = abs(targetAlign_piezo);
load(fullfile(analysis_out,img_fn, [ '_cueAlignLick.mat']));
tempFig=setFigure; hold on;
subplot(2,1,1); hold on;
shadedErrorBar_CV(tt,nanmean(targetAlign_piezo,2),nanstd(targetAlign_piezo,[],2)./sqrt(nTrials),'k',1);
ylabel('Piezo voltage');
xlabel('Time from cue (ms)');
if cue == 1
    xline(0,'g');
elseif cue == 3
    xline(0,'k');
elseif cue == 2
    xline(0,'r');
end
vline(700,'k');
xlim([-1000 4000])
subplot(2,1,2); hold on;
shadedErrorBar_CV(tt,nanmean(lickCueAlign,2),nanstd(lickCueAlign,[],2)./sqrt(nTrials),'k',1);
ylabel('Lick rate');
xlabel('Time from cue (ms)');
if cue == 1
    xline(0,'g');
elseif cue == 3
    xline(0,'k');
elseif cue == 2
    xline(0,'r');
end
vline(700,'k');
sgtitle([mouse ' ' figTitle ' lickAlign and piezoAlign | whole trial']);
savefig(fullfile(analysis_out,img_fn, [ '_avgTrialPiezo_abs.fig']));
saveas(tempFig, [analysis_out img_fn  '_avgTrialPiezo_abs.pdf']);

preRew_lickSearchRange = prewin_frames+lickDelay_frames:prewin_frames+lickDelay_frames+rewDelay_frames;
postRew_lickSearchRange = prewin_frames+rewDelay_frames+lickDelay_frames:prewin_frames+rewDelay_frames+rewDelay_frames;
preRew_piezoAmp = mean(targetAlign_piezo(preRew_lickSearchRange,:),1);% average across frames for each trial
postRew_piezoAmp = mean(targetAlign_piezo(postRew_lickSearchRange,:),1);

tempFig=setFigure; hold on;
subplot(1,2,1);
scatter(preRew_lickBurstHz, preRew_piezoAmp,'ok'); 
xlabel('Lick rate');
ylabel('Piezo voltage');
title('Pre reward');
xlim([0 10]);
ylim([-0.2 0.2]);
axis square;
subplot(1,2,2);
scatter(postRew_lickBurstHz, postRew_piezoAmp,'ok');
xlabel('Lick rate');
ylabel('Piezo voltage');
title('Post reward');
xlim([0 10]);
ylim([-0.2 0.2]);
axis square;
supertitle([mouse figTitle ' Lick vs. Piezo | pre- and post-reward']);
savefig(fullfile(analysis_out,img_fn, [ '_LickvsPiezo_abs.fig']));


%% plot neural activity during trials of top and bottom 10/25% of movement amplitude

[sortPiezoAmp sortPiezoAmp_ind] = sort(preRew_piezoAmp,'ascend');
nnan = sum(isnan(preRew_piezoAmp));
ind_low25piezo_prerew = sortPiezoAmp_ind(1:floor(nTrials/4));
ind_high25piezo_prerew = sortPiezoAmp_ind(nTrials-floor(nTrials/4)+1:end-nnan);
HL_piezo.low25_prerew = nanmean(preRew_piezoAmp(:,ind_low25piezo_prerew),2);
HL_piezo.high25_prerew = nanmean(preRew_piezoAmp(:,ind_high25piezo_prerew),2);
ind_low10piezo_prerew = sortPiezoAmp_ind(1:floor(nTrials/10));
ind_high10piezo_prerew = sortPiezoAmp_ind(nTrials-floor(nTrials/10)+1:end-nnan);
HL_piezo.low10_prerew = nanmean(preRew_piezoAmp(:,ind_low10piezo_prerew),2);
HL_piezo.high10_prerew = nanmean(preRew_piezoAmp(:,ind_high10piezo_prerew),2);

[sortPiezoAmp sortPiezoAmp_ind] = sort(postRew_piezoAmp,'ascend');
nnan = sum(isnan(postRew_piezoAmp));
ind_low25piezo_postrew = sortPiezoAmp_ind(1:floor(nTrials/4));
ind_high25piezo_postrew = sortPiezoAmp_ind(nTrials-floor(nTrials/4)+1:end-nnan);
HL_piezo.low25_postrew = nanmean(postRew_piezoAmp(:,ind_low25piezo_postrew),2);
HL_piezo.high25_postrew = nanmean(postRew_piezoAmp(:,ind_high25piezo_postrew),2);
ind_low10piezo_postrew = sortPiezoAmp_ind(1:floor(nTrials/10));
ind_high10piezo_postrew = sortPiezoAmp_ind(nTrials -floor(nTrials/10)+1:end-nnan);
HL_piezo.low10_postrew = nanmean(postRew_piezoAmp(:,ind_low10piezo_postrew),2);
HL_piezo.high10_postrew = nanmean(postRew_piezoAmp(:,ind_high10piezo_postrew),2);

nIC = size(targetAlign_events,2);
tempFig=setFigure; hold on;
subplot(2,1,1);
shadedErrorBar_CV(tt,nanmean(nanmean(targetAlign_events(:,:,ind_low25piezo_prerew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_low25piezo_prerew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'b',1);
hold on;
shadedErrorBar_CV(tt,nanmean(nanmean(targetAlign_events(:,:,ind_high25piezo_prerew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_high25piezo_prerew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k',1);
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
if cue == 1
    vline(0,'g');
elseif cue == 3
    vline(0,'k');
elseif cue == 2
    vline(0,'r');
end
vline(700,'k');
ylim([0 inf]);
title(['Pre-reward: ' num2str(chop(HL_piezo.low25_prerew,2)) ' vs ' num2str(chop(HL_piezo.high25_prerew,2)) ' V']);
subplot(2,1,2);
shadedErrorBar_CV(tt,nanmean(nanmean(targetAlign_events(:,:,ind_low25piezo_postrew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_low25piezo_postrew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'b',1);
hold on;
shadedErrorBar_CV(tt,nanmean(nanmean(targetAlign_events(:,:,ind_high25piezo_postrew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_high25piezo_postrew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k',1);
xlabel('Time from cue (ms)');
ylabel('Spike rate (Hz)');
if cue == 1
    vline(0,'g');
elseif cue == 3
    vline(0,'k');
elseif cue == 2
    vline(0,'r');
end
vline(700,'k');
ylim([0 inf]);
title(['Post-reward: ' num2str(chop(HL_piezo.low25_postrew,2)) ' vs ' num2str(chop(HL_piezo.high25_postrew,2)) ' V']);
hold on;
sgtitle([mouse ' ' figTitle ' | cue aligned motion: low25(blue) high25(black)']);
hold off;
savefig(fullfile(analysis_out,img_fn, [ '_cueAlignSpiking_byPiezoAmp25_abs.fig']))

tempFig=setFigure; hold on;
subplot(2,1,1);
shadedErrorBar_CV(tt,nanmean(nanmean(targetAlign_events(:,:,ind_low10piezo_prerew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_low10piezo_prerew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'b',1);
hold on
shadedErrorBar_CV(tt,nanmean(nanmean(targetAlign_events(:,:,ind_high10piezo_prerew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_high10piezo_prerew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k',1);
xlabel('Time from cue (ms)');
ylabel('Spike rate (Hz)');
if cue == 1
    vline(0,'g');
elseif cue == 3
    vline(0,'k');
elseif cue == 2
    vline(0,'r');
end
vline(700,'k');
ylim([0 inf]);
title(['Pre-reward: ' num2str(chop(HL_piezo.low10_prerew,2)) ' vs ' num2str(chop(HL_piezo.high10_prerew,2)) ' V']);
subplot(2,1,2);
shadedErrorBar_CV(tt,nanmean(nanmean(targetAlign_events(:,:,ind_low10piezo_postrew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_low10piezo_postrew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'b',1);
hold on;
shadedErrorBar_CV(tt,nanmean(nanmean(targetAlign_events(:,:,ind_high10piezo_postrew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_high10piezo_postrew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k',1);
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
if cue == 1
    vline(0,'g');
elseif cue == 3
    vline(0,'k');
elseif cue == 2
    vline(0,'r');
end
vline(700,'k');
ylim([0 inf]);
title(['Post-reward: ' num2str(chop(HL_piezo.low10_postrew,2)) ' vs ' num2str(chop(HL_piezo.high10_postrew,2)) ' V']);
hold off;
sgtitle([mouse ' ' figTitle ' | piezo bursts by rate: low10 (blue) & high10 (black)']);
savefig(fullfile(analysis_out,img_fn, [ '_cueAlignSpiking_byPiezoAmp10_abs.fig']));

%% sort out trials by looking at time of big movements: 25% and 10%
if loaded==0
    piezo_base_std = nanstd(reshape(targetAlign_piezo(1:prewin_frames,:),[prewin_frames.*nTrials 1]),[],1);
    piezo_base_avg = nanmean(reshape(targetAlign_piezo(1:prewin_frames,:),[prewin_frames.*nTrials 1]),1);
elseif loaded==1
    piezo_base_std = nanstd(targetAlign_piezo(1:prewin_frames,:),[],1);
    piezo_base_avg = nanmean(targetAlign_piezo(1:prewin_frames,:),1);
end
targetAlign_piezo_thresh1 = zeros(size(targetAlign_piezo));%>1std <2std
targetAlign_piezo_thresh1ALL = zeros(size(targetAlign_piezo));%>1std
targetAlign_piezo_thresh2ALL = zeros(size(targetAlign_piezo)); %>2std
targetAlign_piezo_thresh2 = zeros(size(targetAlign_piezo));%>2std <3std
targetAlign_piezo_thresh3 = zeros(size(targetAlign_piezo));%>3std
piezoSearchRange = prewin_frames+lickDelay_frames:size(lickCueAlign,1)-lickSearch_frames;
piezoStart = nan(1,nTrials);
% for each trial, frames that has a voltage>n*std = 1, others = 0
for itrial = 1:nTrials
    targetAlign_piezo_thresh1(find(targetAlign_piezo(:,itrial)>=(piezo_base_avg + piezo_base_std.*1) & targetAlign_piezo(:,itrial)<(piezo_base_avg + piezo_base_std.*2)),itrial) = 1;
    targetAlign_piezo_thresh2(find(targetAlign_piezo(:,itrial)>=(piezo_base_avg + piezo_base_std.*2) & targetAlign_piezo(:,itrial)<(piezo_base_avg + piezo_base_std.*3)),itrial) = 1;
    targetAlign_piezo_thresh3(find(targetAlign_piezo(:,itrial)>=(piezo_base_avg + piezo_base_std.*3)),itrial) = 1;
    targetAlign_piezo_thresh1ALL(find(targetAlign_piezo(:,itrial)>=(piezo_base_avg + piezo_base_std.*1)),itrial) = 1;
    targetAlign_piezo_thresh2ALL(find(targetAlign_piezo(:,itrial)>=(piezo_base_avg + piezo_base_std.*2)),itrial) = 1;
    if find(targetAlign_piezo_thresh3(piezoSearchRange,itrial),1,'first')
        piezoStart(:,itrial) = find(targetAlign_piezo_thresh3(piezoSearchRange,itrial),1,'first');% the start of very big movements in each trial (3std>baseline)
    end
end

[sortPiezoStart sortPiezoStart_ind] = sort(piezoStart,'ascend');
nnan = sum(isnan(piezoStart));
ind_earlypiezo_rew = sortPiezoStart_ind(1:floor((length(sortPiezoStart_ind)-nnan)/4));% top 25% of trials with early big movement
ind_latepiezo_rew = sortPiezoStart_ind((length(sortPiezoStart_ind)-nnan)-floor((length(sortPiezoStart_ind)-nnan)/4)+1:end-nnan);
ind_allearlypiezo25_rew = sortPiezoStart_ind(1:floor(length(sortPiezoStart_ind)/4));
ind_alllatepiezo25_rew = sortPiezoStart_ind(length(sortPiezoStart_ind)-floor(length(sortPiezoStart_ind)/4)+1:end);
ind_allearlypiezo10_rew = sortPiezoStart_ind(1:floor(length(sortPiezoStart_ind)/10)); % top 10% of all trials, regardless of nans. does it make sense? 
ind_alllatepiezo10_rew = sortPiezoStart_ind(length(sortPiezoStart_ind)-floor(length(sortPiezoStart_ind)/10)+1:end);
HL_piezo.early_rew = nanmean(piezoStart(:,ind_earlypiezo_rew),2);
HL_piezo.late_rew = nanmean(piezoStart(:,ind_latepiezo_rew),2);

tempFig=setFigure; hold on;
subplot(2,1,1);
shadedErrorBar_CV(tt,nanmean(nanmean(targetAlign_events(:,:,ind_earlypiezo_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_earlypiezo_rew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'b',1);
hold on;
shadedErrorBar_CV(tt,nanmean(nanmean(targetAlign_events(:,:,ind_latepiezo_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_latepiezo_rew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k',1);
scatter(piezoStart(:,ind_earlypiezo_rew).*(1000./frameRateHz),zeros(1,length(ind_earlypiezo_rew))+(.5+max(nanmean(nanmean(targetAlign_events(:,:,ind_latepiezo_rew),3),2).*(1000./frameRateHz))),'ob');
scatter(piezoStart(:,ind_latepiezo_rew).*(1000./frameRateHz),zeros(1,length(ind_latepiezo_rew))+(.5+max(nanmean(nanmean(targetAlign_events(:,:,ind_latepiezo_rew),3),2).*(1000./frameRateHz))),'ok');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 inf]);
if cue == 1
    vline(0,'g');
elseif cue == 3
    vline(0,'k');
elseif cue == 2
    vline(0,'r');
end
vline(700,'k');
title(['earliest 25% vs. latest 25% of motion: ' num2str(chop(HL_piezo.early_rew.*(1000./frameRateHz),2)) ' vs ' num2str(chop(HL_piezo.late_rew.*(1000./frameRateHz),2)) ' ms']);

subplot(2,1,2);
shadedErrorBar_CV(tt,nanmean(nanmean(targetAlign_events(:,:,ind_allearlypiezo10_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_allearlypiezo10_rew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'b',1);
hold on;
shadedErrorBar_CV(tt,nanmean(nanmean(targetAlign_events(:,:,ind_alllatepiezo10_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_alllatepiezo10_rew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k',1);
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 inf]);
if cue == 1
    vline(0,'g');
elseif cue == 3
    vline(0,'k');
elseif cue == 2
    vline(0,'r');
end
vline(700,'k');
title('all earliest 10% vs latest 10% of motion');
sgtitle([mouse ' ' figTitle '- motion bursts by latency: early (blue) & late (black)']);
savefig(fullfile(analysis_out,img_fn, [ '_cueAlignSpiking_byPiezoLatency_abs.fig']));


[TrialMovement, TrialMovement_ind] = sort(nanmean(targetAlign_piezo,1));
mostTrialMovement_ind = TrialMovement_ind(:,end-(floor(size(targetAlign_piezo,2)./10)):end);
leastTrialMovement_ind = TrialMovement_ind(:,1:floor(size(targetAlign_piezo,2)./10));

mostMovementAlignedEvents = nanmean(targetAlign_events(:,:,mostTrialMovement_ind),3);
leastMovementAlignedEvents = nanmean(targetAlign_events(:,:,leastTrialMovement_ind),3);
mostTrialMovement = targetAlign_piezo(:,mostTrialMovement_ind);
leastTrialMovement = targetAlign_piezo(:,leastTrialMovement_ind);


tempFig = figure;
subplot(2,1,1);hold on;
plot(tt,nanmean(mostMovementAlignedEvents,2).*(1000./frameRateHz),'color',[1 0.5 0]);
hold on;
plot(tt,nanmean(leastMovementAlignedEvents,2).*(1000./frameRateHz),'color',[0 0.5 1]);
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
vline(0,'k');
vline(767,'k');
legend('top 10%','bottom 10%','Location','NorthEast');
title('Neural data aligned to trials with top 10% (org) and bottom 10% (blu) of movement');
subplot(2,1,2);hold on;
plot(tt,nanmean(mostTrialMovement,2),'color',[1 0.5 0]);
hold on;
plot(tt,nanmean(leastTrialMovement,2),'color',[0 0.5 1]);
xlabel('Time from cue');
ylabel('Piezo voltage');
vline(0,'k');
vline(767,'k');
title('Piezo data for trials with top 10% (org) and bottom 10% (blu) of movement');
hold off;
sgtitle([num2str(size(targetAlign_events,2)) ' neurons from ' num2str(size(expt,2)) ' animals | neural aligned to trials with most and least movement']);
savefig(fullfile(analysis_out,img_fn, [ '_mostLeastMotionAlign.fig']));
saveas(tempFig, [analysis_out img_fn  '_mostLeastMotionAlign.pdf']);


    piezoWholeSession = []; piezoWholeSessionPostReward = []; piezoWholeSessionPreCue = [];
for itrial = 1:size(targetAlign_piezo,2)
piezoWholeSession = [piezoWholeSession; targetAlign_piezo(:,itrial)];
piezoWholeSessionPreCue = [piezoWholeSessionPreCue; targetAlign_piezo(1:50,itrial)];
piezoWholeSessionPostReward = [piezoWholeSessionPostReward; targetAlign_piezo(73:150,itrial)];
end

uWholeSession = nanmean(piezoWholeSession);
oWholeSession = nanstd(piezoWholeSession,[],1);
uWholeSessionPreCue = nanmean(piezoWholeSessionPreCue);
oWholeSessionPreCue = nanstd(piezoWholeSessionPreCue,[],1);
uWholeSessionPostReward = nanmean(piezoWholeSessionPostReward);
oWholeSessionPostReward = nanstd(piezoWholeSessionPostReward,[],1);


tempFig = setFigure; hold on;
plot(1:length(piezoWholeSession),piezoWholeSession,'k');
hold on;
yline(uWholeSession,'r--');
yline((oWholeSession.*1)+uWholeSession,'Color',[0.8500, 0.3250, 0.0980]);
yline((oWholeSession.*2)+uWholeSession,'Color',[0.9290, 0.6940, 0.1250]);
yline((oWholeSession.*3)+uWholeSession,'g');
yline((oWholeSession.*4)+uWholeSession,'b');
yline((oWholeSession.*5)+uWholeSession,'c');
yline((oWholeSession.*6)+uWholeSession,'m');
legend({'piezo','mean','1SD','2SD','3SD','4SD','5SD','6SD'},'Location','northeast');
ylim([0 1])
title('piezo transients across imaging session');
supertitle('mean and s.d. found using whole dataset');


tempFig = setFigure; hold on;
plot(1:length(piezoWholeSession),piezoWholeSession,'k');
hold on;
yline(uWholeSessionPreCue,'r--');
yline((oWholeSessionPreCue.*1)+uWholeSessionPreCue,'Color',[0.8500, 0.3250, 0.0980]);
yline((oWholeSessionPreCue.*2)+uWholeSessionPreCue,'Color',[0.9290, 0.6940, 0.1250]);
yline((oWholeSessionPreCue.*3)+uWholeSessionPreCue,'g');
yline((oWholeSessionPreCue.*4)+uWholeSessionPreCue,'b');
yline((oWholeSessionPreCue.*5)+uWholeSessionPreCue,'c');
yline((oWholeSessionPreCue.*6)+uWholeSessionPreCue,'m');
legend({'piezo','mean','1SD','2SD','3SD','4SD','5SD','6SD'},'Location','northeast');
ylim([0 1])
title('piezo transients across of imaging session');
supertitle('mean and s.d. found using only pre-cue frames');


tempFig = setFigure; hold on;
plot(1:length(piezoWholeSessionPostReward),piezoWholeSessionPostReward,'k');
hold on;
yline(uWholeSessionPreCue,'r--');
yline((oWholeSessionPreCue.*1)+uWholeSessionPreCue,'Color',[0.8500, 0.3250, 0.0980]);
yline((oWholeSessionPreCue.*2)+uWholeSessionPreCue,'Color',[0.9290, 0.6940, 0.1250]);
yline((oWholeSessionPreCue.*3)+uWholeSessionPreCue,'g');
yline((oWholeSessionPreCue.*4)+uWholeSessionPreCue,'b');
yline((oWholeSessionPreCue.*5)+uWholeSessionPreCue,'c');
yline((oWholeSessionPreCue.*6)+uWholeSessionPreCue,'m');
legend({'piezo','mean','1SD','2SD','3SD','4SD','5SD','6SD'},'Location','northeast');
ylim([0 1])
title('piezo transients across post-reward frames of imaging session');
supertitle('mean and s.d. found using pre-cue frames');

%create a logical (150XnTrial) indicating frames where value resides 2sd or
%1sd above the average for that session

    motionPeaks = nan(150,size(targetAlign_piezo,2));
    cutoffPiezo = uWholeSessionPreCue+(2.*oWholeSessionPreCue);
for itrial = 1:size(targetAlign_piezo,2)
    for iframe = 1:150
        if targetAlign_piezo(iframe,itrial) > cutoffPiezo
        motionPeaks(iframe,itrial) = 1;
        elseif targetAlign_piezo(iframe,itrial) <= cutoffPiezo
        motionPeaks(iframe,itrial) = 0;
        end
    end
end


%% lick align formatting

    nFrames = mworks.counterValues{nTrials}(end);
    rewDelay_frames = round((767/1000).*frameRateHz);% there's 700ms between the cue and the reward delivery !!!!! if you change tooFastMs in the varibles in MWorks, this needs to be changed
    nIC = size(targetAlign_events,2); %you get targetalign_events from neural align, this is the neural data (in spikes) aligned to cue

    
    %window blocks
    lickDelay_frames =  round(0.1.*frameRateHz);
    lickSearch_frames =  round(0.3.*frameRateHz);
    postLick_frames = round(0.5.*frameRateHz);

    %search windows
    preRew_lickSearchRange_700ms = prewin_frames+lickDelay_frames:prewin_frames+lickDelay_frames+rewDelay_frames;
    postRew_lickSearchRange = prewin_frames+rewDelay_frames+lickDelay_frames:size(motionPeaks,1)-lickSearch_frames-postLick_frames;% need to '-lickSearch_frames-postLick_frames' because later the inds needs to + lickSearch_frames or +postLick_frames, this is just for the following inds to be within matrix dimentions

    % piezo align blanks
    firstPostRew_piezoAlignEvents = nan(3.*postLick_frames,nIC,nTrials); 
    firstPostRew_piezoAlign = nan(3.*postLick_frames,nTrials);    
    lastPreRew_piezoAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    lastPreRew_piezoAlign = nan(3.*postLick_frames,nTrials);    
    preRewPiezoTrials = zeros(1,nTrials);
    postRewPiezoTrials = zeros(1,nTrials);
    lastPreRewPiezoFrame = nan(1,nTrials);
    firstPostRewPiezoFrame = nan(1,nTrials);
   
    tTInx = 1; ntrials=0;
      
    for itrial = 1:nTrials %nTrials%
     
        ind_postP = [];
        ind_preP = [];
        
        if ~isnan(cTargetOn(itrial))
            if cTargetOn(itrial)+postwin_frames-1 < nFrames
                counterTimes = mworks.counterTimesUs{itrial};
                counterVals = mworks.counterValues{itrial};
            
                ind_preP = intersect(preRew_lickSearchRange_700ms,find(motionPeaks(:,itrial))); %finds every instance of a lick that occurs within the search window
                ind_postP = intersect(postRew_lickSearchRange,find(motionPeaks(:,itrial))); %finds every instance of a lick that occurs [after the reward delivery]
                
              
                ind_preP = find(motionPeaks(prewin_frames:prewin_frames+rewDelay_frames,itrial),1,'last');
         
                if ~isempty(ind_preP)
                    lastPreRewPiezoFrame(1,itrial) = ind_preP;
                    preRewPiezoTrials(1,itrial) = 1;
                    lastPreRew_piezoAlignEvents(:,:,itrial) = targetAlign_events(ind_preP-postLick_frames+prewin_frames:ind_preP+postLick_frames+postLick_frames-1+prewin_frames,:,itrial);
                    lastPreRew_piezoAlign(:,itrial) = motionPeaks(ind_preP-postLick_frames+prewin_frames:ind_preP+postLick_frames+postLick_frames-1+prewin_frames,itrial);
                end
                
                ind_postP = find(motionPeaks(prewin_frames+rewDelay_frames:prewin_frames+postwin_frames-postLick_frames-postLick_frames,itrial),1,'first');
          
                if ~isempty(ind_postP)
                    firstPostRewPiezoFrame(1,itrial) = ind_postP;
                    postRewPiezoTrials(1,itrial) = 1;
                    firstPostRew_piezoAlignEvents(:,:,itrial) = targetAlign_events(ind_postP-postLick_frames+prewin_frames+rewDelay_frames:ind_postP+postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,:,itrial);
                    firstPostRew_piezoAlign(:,itrial) = motionPeaks(ind_postP-postLick_frames+prewin_frames+rewDelay_frames:ind_postP+postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,itrial);
                end

            tTInx = tTInx+1;
            end
        end
    end
     
    
%% piezo align formatting

ind_nan = find(isnan(targetAlign_piezo(1,:)));

g= [0.4660 0.6740 0.1880];
r= [0.6350 0.0780 0.1840];

gBlock2 = logical(double(cell2mat(mworks.tBlock2TrialNumber(1,1:nTrials))));

ndends=0; ntrials=0; lastPrePiezoEvent=[]; firstPostPiezoEvent=[];  
csPluslastPrePiezoEvent=[]; csPlusfirstPostPiezoEvent=[]; csMinuslastPrePiezoEvent=[]; csMinusfirstPostPiezoEvent=[];

lastPrePiezoEvent = [lastPrePiezoEvent (nanmean(nanmean(lastPreRew_piezoAlignEvents,3),2))]; 
firstPostPiezoEvent = [firstPostPiezoEvent (nanmean(nanmean(firstPostRew_piezoAlignEvents,3),2))];
csPluslastPrePiezoEvent = [csPluslastPrePiezoEvent (nanmean(nanmean(lastPreRew_piezoAlignEvents(:,:,~gBlock2),3),2))];
csPlusfirstPostPiezoEvent = [csPlusfirstPostPiezoEvent (nanmean(nanmean(firstPostRew_piezoAlignEvents(:,:,~gBlock2),3),2))];
csMinuslastPrePiezoEvent = [csMinuslastPrePiezoEvent (nanmean(nanmean(lastPreRew_piezoAlignEvents(:,:,gBlock2),3),2))];
csMinusfirstPostPiezoEvent = [csMinusfirstPostPiezoEvent (nanmean(nanmean(firstPostRew_piezoAlignEvents(:,:,gBlock2),3),2))];

    tl_rew = (1-postLick_frames:(postLick_frames*2)).*(1000./frameRateHz); %orig :: -500ms to 1000 ms around reward
    %tl_rew = (1-postPiezo_frames*2:(postPiezo_frames)).*(1000./frameRateHz); %novel :: -1000ms to 500ms around reward
    tempFig=setFigure; hold on; % align neural and Piezoing data to the first and last Piezo, respectively
    subplot(2,2,1); hold on;
    shadedErrorBar_CV(tl_rew, lastPrePiezoEvent.*(1000./frameRateHz), (lastPrePiezoEvent.*(1000./frameRateHz))./sqrt(nIC),'k',1);
    hold on;
    title('Last movement before reward');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 3]);
    subplot(2,2,3);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(lastPreRew_piezoAlign,2).*(1000./frameRateHz), (nanstd(lastPreRew_piezoAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(lastPreRew_piezoAlign)))),'k',1);
    xlabel('Time from movement (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 35]);
    title([num2str(sum(preRewPiezoTrials)) 'trials with pre-reward movement']);
    subplot(2,2,2);hold on;
    shadedErrorBar_CV(tl_rew, firstPostPiezoEvent.*(1000./frameRateHz), (firstPostPiezoEvent.*(1000./frameRateHz))./sqrt(nIC),'k',1);
    title('First movement after reward');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 3]);
    subplot(2,2,4);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(firstPostRew_piezoAlign,2).*(1000./frameRateHz), (nanstd(firstPostRew_piezoAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(firstPostRew_piezoAlign)))),'k',1);
    xlabel('Time from movement (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 35]);
    title([num2str(sum(postRewPiezoTrials)) 'trials with post-reward movement']);
    sgtitle([num2str(max(nIC)) ' neurons from ' num2str(size(expt,2)) ' animals | Movement relative to reward']);
    hold off;
    if seshType ~= 2
    savefig(fullfile(analysis_out,img_fn, '_lastVsFirstPiezo.fig'));
    saveas(tempFig, [analysis_out img_fn  '_lastVsFirstPiezo.pdf']);
    else
    savefig(fullfile(analysis_out,img_fn, [ '_lastVsFirstPiezo.fig']));
    saveas(tempFig, [analysis_out img_fn  '_lastVsFirstPiezo.pdf']);
    end

 
    %% CS+ and CS- for spike rate aligned to first and last movement (2o above u) in reference to reward delivery
 b2Idx = double(cell2mat(mworks.tBlock2TrialNumber(1,1:nTrials))); rewIdx = logical(b2Idx~=1);
    b2Idx = logical(b2Idx);
    
    
csPlusfirstPostRew_piezoAlign = firstPostRew_piezoAlign(:,~gBlock2);
csMinusfirstPostRew_piezoAlign = firstPostRew_piezoAlign(:,gBlock2);
csPluslastPreRew_piezoAlign = lastPreRew_piezoAlign(:,~gBlock2);
csMinuslastPreRew_piezoAlign = lastPreRew_piezoAlign(:,gBlock2);

csPluspostRewPiezoTrials = postRewPiezoTrials(:,~gBlock2);
csMinuspostRewPiezoTrials = postRewPiezoTrials(:,gBlock2);
csPluspreRewPiezoTrials = preRewPiezoTrials(:,~gBlock2);
csMinuspreRewPiezoTrials = preRewPiezoTrials(:,gBlock2);

    %CS+
    if seshType == 2 && sum(rewIdx)>0
     tl_rew = (1-postLick_frames:(postLick_frames*2)).*(1000./frameRateHz); %orig :: -500ms to 1000 ms around reward
    %tl_rew = (1-postPiezo_frames*2:(postPiezo_frames)).*(1000./frameRateHz); %novel :: -1000ms to 500ms around reward
    tempFig=setFigure; hold on; % align neural and Piezoing data to the first and last Piezo, respectively
    subplot(2,2,1); hold on;
    shadedErrorBar_CV(tl_rew, csPluslastPrePiezoEvent.*(1000./frameRateHz), (csPluslastPrePiezoEvent.*(1000./frameRateHz))./sqrt(nIC),'k',1);
    hold on;
    title('Last movement before reward');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 4]);
    subplot(2,2,3);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(lastPreRew_piezoAlign(:,rewIdx),2).*(1000./frameRateHz), (nanstd(lastPreRew_piezoAlign(:,rewIdx),[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(lastPreRew_piezoAlign(:,rewIdx))))),'k',1);
    xlabel('Time from movement (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 35]);
    title([num2str(sum(csPluspreRewPiezoTrials)) 'CS+ trials with pre-reward movement']);
    subplot(2,2,2);hold on;
    shadedErrorBar_CV(tl_rew, csPlusfirstPostPiezoEvent.*(1000./frameRateHz), (csPlusfirstPostPiezoEvent.*(1000./frameRateHz))./sqrt(nIC),'k',1);
    title('First movement after reward');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 4]);
    subplot(2,2,4);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(firstPostRew_piezoAlign(:,rewIdx),2).*(1000./frameRateHz), (nanstd(firstPostRew_piezoAlign(:,rewIdx),[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(firstPostRew_piezoAlign(:,rewIdx))))),'k',1);
    xlabel('Time from movement (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 35]);
    title([num2str(sum(csPluspostRewPiezoTrials)) 'CS+ trials with post-reward movement']);
    sgtitle([num2str(max(nIC)) ' neurons from ' num2str(size(expt,2)) ' animals | Movement relative to reward']);
    hold off;
    if seshType ~= 2
    savefig(fullfile(analysis_out,img_fn, '_csPluslastVsFirstPiezo.fig'));
    saveas(tempFig, [analysis_out img_fn  '_csPluslastVsFirstPiezo.pdf']);
    else
    savefig(fullfile(analysis_out,img_fn, ['_csPluslastVsFirstPiezo.fig']));
    saveas(tempFig, [analysis_out img_fn  '_csPluslastVsFirstPiezo.pdf']);
    end
    
    end
   
    %CS-
    if seshType == 2 && sum(b2Idx)>0
    tl_rew = (1-postLick_frames:(postLick_frames*2)).*(1000./frameRateHz); %orig :: -500ms to 1000 ms around reward
    %tl_rew = (1-postPiezo_frames*2:(postPiezo_frames)).*(1000./frameRateHz); %novel :: -1000ms to 500ms around reward
    tempFig=setFigure; hold on; % align neural and Piezoing data to the first and last Piezo, respectively
    subplot(2,2,1); hold on;
    shadedErrorBar_CV(tl_rew, csMinuslastPrePiezoEvent.*(1000./frameRateHz), (csMinuslastPrePiezoEvent.*(1000./frameRateHz))./sqrt(nIC),'r');
    hold on;
    title('Last movement before reward');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 4]);
    subplot(2,2,3);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(lastPreRew_piezoAlign(:,b2Idx),2).*(1000./frameRateHz), (nanstd(lastPreRew_piezoAlign(:,b2Idx),[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(lastPreRew_piezoAlign(:,b2Idx))))),'r');
    xlabel('Time from movement (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 35]);
    title([num2str(sum(csMinuspreRewPiezoTrials)) 'CS- trials with pre-reward movement']);
    subplot(2,2,2);hold on;
    shadedErrorBar_CV(tl_rew, csMinusfirstPostPiezoEvent.*(1000./frameRateHz), (csMinusfirstPostPiezoEvent.*(1000./frameRateHz))./sqrt(nIC),'r');
    title('First movement after reward');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 4]);
    subplot(2,2,4);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(firstPostRew_piezoAlign(:,b2Idx),2).*(1000./frameRateHz), (nanstd(firstPostRew_piezoAlign(:,b2Idx),[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(firstPostRew_piezoAlign(:,b2Idx))))),'r');
    xlabel('Time from movement (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 35]);
    title([num2str(sum(csMinuspostRewPiezoTrials)) 'CS- trials with post-reward movement']);
    sgtitle([num2str(max(nIC)) ' neurons from ' num2str(size(expt,2)) ' animals | Movement relative to reward']);
    hold off;
    if seshType ~= 2
    savefig(fullfile(analysis_out,img_fn, '_csMinuslastVsFirstPiezo.fig'));
    saveas(tempFig, [analysis_out img_fn  '_csMinuslastVsFirstPiezo.pdf']);
    else
    savefig(fullfile(analysis_out,img_fn, [ '_csMinuslastVsFirstPiezo.fig']));
    saveas(tempFig, [analysis_out img_fn  '_csMinuslastVsFirstPiezo.pdf']);
    end
    
    end


    
    
    
save(fullfile(analysis_out,img_fn,['_cueAlignPiezo.mat']), ...
    'targetAlign_piezo','targetAlign_piezo_thresh1','targetAlign_piezo_thresh1ALL',...
    'targetAlign_piezo_thresh2','targetAlign_piezo_thresh2ALL','targetAlign_piezo_thresh3', ...
    'preRew_piezoAmp', 'postRew_piezoAmp', 'ind_low25piezo_prerew', 'ind_high25piezo_prerew',...
    'ind_low25piezo_postrew', 'ind_high25piezo_postrew', 'ind_low10piezo_prerew', ...
    'ind_high10piezo_prerew','ind_low10piezo_postrew', 'ind_high10piezo_postrew',...
    'ind_earlypiezo_rew', 'ind_latepiezo_rew','ind_allearlypiezo25_rew', ...
    'ind_alllatepiezo25_rew', 'ind_allearlypiezo10_rew', 'ind_alllatepiezo10_rew',...
    'HL_piezo','piezo_frames');

