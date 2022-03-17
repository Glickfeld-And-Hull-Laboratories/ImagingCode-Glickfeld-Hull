close all
clear all
clear all global
clc
dataset = 'plaidDiscrim_exptList';
%dataset = 'i484_passive_ExptList';
%dataset = 'plaidDiscrim_Passive_exptList';
eval(dataset);
iexp = 10;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
if strcmp(expt(iexp).folder,'lindsey')
    data_base = LG_base;
elseif strcmp(expt(iexp).folder,'camaron')
    data_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron';
end
mouse = expt(iexp).mouse;
date = expt(iexp).date;
nrun = size(expt(iexp).runs,1);
run_str = ['runs']; 
for irun = 1:nrun
    run_str = [run_str '-' expt(iexp).runs(irun,:)];
end

%%
CD = [data_base '\Data\2P_images\' expt(iexp).mouse '\' expt(iexp).date '\' expt(iexp).runs];
cd(CD);

fn = [expt(iexp).runs '_000_000_eye.mat'];

%load data
data_temp = load(fn);
data_temp = squeeze(data_temp.data);

%crop frames to match mworks data
load(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
nFrames = input.counterValues{end}(end);
data = data_temp(:,:,1:nFrames);      % the raw images...
clear data_temp
%% Crop image to isolate pupil 
%(bright spots can be mistaken for pupil)
[data_crop rect] = cropEyeData(data);

%% measure pupil position/diameter
rad_range = [5 20]; %adjust to expected range of pupil size (if low end is too small then may find noisy bright stuff)
Eye_data = extractEyeData(data_crop,rad_range);
%if pupil not found reliably, adjust the image cropping or the rad_range

%% align to stimulus presentation
[rad centroid] = alignEyeData(Eye_data,input);

%% average by stimulus
%wheel trajectory
[qVals qThresh] = wheelTrajectory(input, 15);
tt_wheel = -8000:10000;

load(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
dirs = unique(tDir);
dirs(find(isnan(dirs))) = [];
nDir = length(dirs);
centroid.dir = zeros(size(centroid.tc,1), 2, nDir,2);
rad.dir = zeros(size(centroid.tc,1), nDir,2);
wheel.dir = zeros(size(qVals,1), nDir,2);
for iDir = 1:nDir
    for iMask = 1:2
        ind = intersect(find(tMaskTrial == iMask-1),find(tDir==dirs(iDir)));
        centroid.dir(:,:,iDir,iMask) = nanmean(centroid.tc(:,:,ind),3);
        rad.dir(:,iDir,iMask) = nanmean(rad.tc(:,ind),2);
        wheel.dir(:,iDir,iMask) = nanmean(qVals(:,ind),2);
    end
end
centroid.dircat = zeros(size(centroid.tc,1), 2, 2,2);
rad.dircat = zeros(size(centroid.tc,1), 2,2);
wheel.dircat = zeros(size(qVals,1),2,2);
for iMask = 1:2
    ind = intersect(find(tMaskTrial == iMask-1),find(tDir<45));
    centroid.dircat(:,:,1,iMask) = nanmean(centroid.tc(:,:,ind),3);
    rad.dircat(:,1,iMask) = nanmean(rad.tc(:,ind),2);
    wheel.dircat(:,1,iMask) = nanmean(qVals(:,ind),2);
    ind = intersect(find(tMaskTrial == iMask-1),find(tDir>45));
    centroid.dircat(:,:,2,iMask) = nanmean(centroid.tc(:,:,ind),3);
    rad.dircat(:,2,iMask) = nanmean(rad.tc(:,ind),2);
    wheel.dircat(:,2,iMask) = nanmean(qVals(:,ind),2);
end
wheel.dir = wheel.dir./32;
wheel.dircat = wheel.dircat./32;

tt = (1-input.frameRateHz:input.frameRateHz*3).*(1000./frameRateHz);
newcolors = parula(nDir);
figure;
subplot(4,2,1)
plot(tt,squeeze(centroid.dir(:,1,:,1)-mean(centroid.dir(1:input.frameRateHz,1,:,1),1)))
colororder(newcolors)
title('Horizontal- grating')
xlim([-500 1000])
ylim([-1 3])
xlabel('Time from stimulus (ms)')
legend(num2str(dirs'),'location','northwest')
subplot(4,2,3)
plot(tt,squeeze(centroid.dir(:,2,:,1)-mean(centroid.dir(1:input.frameRateHz,2,:,1),1)))
title('Vertical- grating')
xlim([-500 1000])
ylim([-1 3])
xlabel('Time from stimulus (ms)')
subplot(4,2,2)
plot(tt,squeeze(centroid.dir(:,1,:,2)-mean(centroid.dir(1:input.frameRateHz,1,:,2),1)))
title('Horizontal- plaid')
xlim([-500 1000])
ylim([-1 3])
xlabel('Time from stimulus (ms)')
subplot(4,2,4)
plot(tt,squeeze(centroid.dir(:,2,:,2)-mean(centroid.dir(1:input.frameRateHz,2,:,2),1)))
title('Vertical- plaid')
xlim([-500 1000])
ylim([-1 3])
xlabel('Time from stimulus (ms)')
subplot(4,2,5)
plot(tt,squeeze(rad.dir(:,:,1)-mean(rad.dir(1:input.frameRateHz,:,1),1)))
title('Radius- grating')
xlim([-500 1000])
ylim([-1 3])
xlabel('Time from stimulus (ms)')
subplot(4,2,6)
plot(tt,squeeze(rad.dir(:,:,2)-mean(rad.dir(1:input.frameRateHz,:,2),1)))
title('Radius- plaid')
xlim([-500 1000])
ylim([-1 3])
xlabel('Time from stimulus (ms)')
subplot(4,2,7)
plot(tt_wheel,squeeze(wheel.dir(:,:,1)-mean(wheel.dir(7000:8000,:,1),1)))
title('Wheel- grating')
xlim([-500 1000])
ylim([-5 5])
xlabel('Time from stimulus (ms)')
subplot(4,2,8)
plot(tt_wheel,squeeze(wheel.dir(:,:,2)-mean(wheel.dir(7000:8000,:,2),1)))
title('Wheel- plaid')
xlim([-500 1000])
xlabel('Time from stimulus (ms)')
ylim([-5 5])
sgtitle([date ' ' mouse])
print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_eyePlots.pdf']),'-dpdf','-fillpage')

figure;
newcolors = defaultPlotColors();
colororder(newcolors)
subplot(4,2,1)
plot(tt,squeeze(centroid.dircat(:,1,:,1)-mean(centroid.dircat(1:input.frameRateHz,1,:,1),1)))
title('Horizontal- grating')
xlim([-500 1000])
ylim([-1 3])
xlabel('Time from stimulus (ms)')
legend({'Right','Up'},'location','northwest')
subplot(4,2,3)
plot(tt,squeeze(centroid.dircat(:,2,:,1)-mean(centroid.dircat(1:input.frameRateHz,2,:,1),1)))
title('Vertical- grating')
xlim([-500 1000])
ylim([-1 3])
xlabel('Time from stimulus (ms)')
subplot(4,2,2)
plot(tt,squeeze(centroid.dircat(:,1,:,2)-mean(centroid.dircat(1:input.frameRateHz,1,:,2),1)))
title('Horizontal- plaid')
xlim([-500 1000])
ylim([-1 3])
xlabel('Time from stimulus (ms)')
subplot(4,2,4)
plot(tt,squeeze(centroid.dircat(:,2,:,2)-mean(centroid.dircat(1:input.frameRateHz,2,:,2),1)))
title('Vertical- plaid')
xlim([-500 1000])
ylim([-1 3])
xlabel('Time from stimulus (ms)')
subplot(4,2,5)
plot(tt,squeeze(rad.dircat (:,:,1)-mean(rad.dircat (1:input.frameRateHz,:,1),1)))
colororder(newcolors)
title('Radius- grating')
xlim([-500 1000])
ylim([-1 3])
xlabel('Time from stimulus (ms)')
subplot(4,2,6)
plot(tt,squeeze(rad.dircat (:,:,2)-mean(rad.dircat (1:input.frameRateHz,:,2),1)))
title('Radius- plaid')
xlim([-500 1000])
xlabel('Time from stimulus (ms)')
ylim([-1 3])
subplot(4,2,7)
plot(tt_wheel,squeeze(wheel.dircat(:,:,1)-mean(wheel.dircat(7000:8000,:,1),1)))
title('Wheel- grating')
xlim([-500 1000])
ylim([-5 5])
subplot(4,2,8)
plot(tt_wheel,squeeze(wheel.dircat(:,:,2)-mean(wheel.dircat(7000:8000,:,2),1)))
title('Wheel- plaid')
xlim([-500 1000])
ylim([-5 5])
sgtitle([date ' ' mouse])
print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_eyePlots_dircat.pdf']),'-dpdf','-fillpage')


wheel.max_val = max(abs(qVals(8000:9000,:)),[],1);
ind_stat = find(wheel.max_val<30);
centroid.stat = zeros(size(centroid.tc,1), 2, 2,2);
rad.stat = zeros(size(centroid.tc,1), 2,2);
wheel.stat = zeros(size(qVals,1),2,2);

for iMask = 1:2
    ind = intersect(ind_stat,intersect(find(tMaskTrial == iMask-1),find(tDir<45)));
    length(ind)
    centroid.stat(:,:,1,iMask) = nanmean(centroid.tc(:,:,ind),3);
    rad.stat(:,1,iMask) = nanmean(rad.tc(:,ind),2);
    wheel.stat(:,1,iMask) = nanmean(qVals(:,ind),2);
    ind = intersect(ind_stat,intersect(find(tMaskTrial == iMask-1),find(tDir>45)));
    length(ind)
    centroid.stat(:,:,2,iMask) = nanmean(centroid.tc(:,:,ind),3);
    rad.stat(:,2,iMask) = nanmean(rad.tc(:,ind),2);
    wheel.stat(:,2,iMask) = nanmean(qVals(:,ind),2);
end
wheel.stat = wheel.stat./32;
figure;
newcolors = defaultPlotColors();
colororder(newcolors)
subplot(4,2,1)
plot(tt,squeeze(centroid.stat(:,1,:,1)-mean(centroid.stat(1:input.frameRateHz,1,:,1),1)))
title('Horizontal- grating')
xlim([-500 1000])
ylim([-1 3])
xlabel('Time from stimulus (ms)')
legend({'Right','Up'},'location','northwest')
subplot(4,2,3)
plot(tt,squeeze(centroid.stat(:,2,:,1)-mean(centroid.stat(1:input.frameRateHz,2,:,1),1)))
title('Vertical- grating')
xlim([-500 1000])
ylim([-1 3])
xlabel('Time from stimulus (ms)')
subplot(4,2,2)
plot(tt,squeeze(centroid.stat(:,1,:,2)-mean(centroid.stat(1:input.frameRateHz,1,:,2),1)))
title('Horizontal- plaid')
xlim([-500 1000])
ylim([-1 3])
xlabel('Time from stimulus (ms)')
subplot(4,2,4)
plot(tt,squeeze(centroid.stat(:,2,:,2)-mean(centroid.stat(1:input.frameRateHz,2,:,2),1)))
title('Vertical- plaid')
xlim([-500 1000])
ylim([-1 3])
xlabel('Time from stimulus (ms)')
subplot(4,2,5)
plot(tt,squeeze(rad.stat (:,:,1)-mean(rad.stat (1:input.frameRateHz,:,1),1)))
colororder(newcolors)
title('Radius- grating')
xlim([-500 1000])
ylim([-1 3])
xlabel('Time from stimulus (ms)')
subplot(4,2,6)
plot(tt,squeeze(rad.stat (:,:,2)-mean(rad.stat (1:input.frameRateHz,:,2),1)))
title('Radius- plaid')
xlim([-500 1000])
xlabel('Time from stimulus (ms)')
ylim([-1 3])
subplot(4,2,7)
plot(tt_wheel,squeeze(wheel.stat(:,:,1)-mean(wheel.stat(7000:8000,:,1),1)))
title('Wheel- grating')
xlim([-500 1000])
ylim([-5 5])
subplot(4,2,8)
plot(tt_wheel,squeeze(wheel.stat(:,:,2)-mean(wheel.stat(7000:8000,:,2),1)))
title('Wheel- plaid')
xlim([-500 1000])
ylim([-5 5])
sgtitle([date ' ' mouse])
print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_eyePlots_stat.pdf']),'-dpdf','-fillpage')


save(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_eyeData.mat']),'rect','rad','centroid','wheel')
