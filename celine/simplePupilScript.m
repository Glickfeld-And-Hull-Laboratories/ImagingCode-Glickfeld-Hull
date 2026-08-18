
% navigate to the folder in question
ImgFolder = '000';

fn = [ImgFolder '_000_000_eye.mat'];
data_temp = load(fn);
data = squeeze(data_temp.data);
clear data_temp


data_avg = mean(data,3);
imagesc(data_avg); movegui('center') % plot mean image

ax = gca;
rect = round(getrect(ax));
datat = data_avg(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));
figure; imagesc(datat)

data_crop = data(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3),:);

rad_range = [3 35]; %adjust to expected range of pupil size (if low end is too small then may find noisy bright stuff)
Eye_data = extractEyeData(data_crop,rad_range);

Rad_temp = sqrt(Eye_data.Area./pi);
Centroid_temp = Eye_data.Centroid;
Rad_temp(Eye_data.badFrames,:) =nan(length(Eye_data.badFrames),1);
Centroid_temp(Eye_data.badFrames,:) = nan(length(Eye_data.badFrames),2);

rad.tc = zeros((nOff/2)+nOn, ntrials);
centroid.tc = zeros((nOff/2)+nOn,2, ntrials);
calib = 1/26.6; %mm per pixel
for it = 1:ntrials
    crange = stimOn(it)-(nOff/2):stimOn(it)+nOn-1;
    rad.tc(:,it) = Rad_temp(crange,:).*calib;
    centroid.tc(:,:,it) = Centroid_temp(crange,:).*calib;
end
rad.stim = nanmean(rad.tc(1+nOff/2:end,:),1); %average over stim window
centroid.stim = squeeze(nanmean(centroid.tc(1+nOff/2:end,:,:),1))./0.025;

save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_eye.mat']), 'centroid', 'rad', 'rect', 'rad_range')

%plot eye position across trials 
figure; subplot(2,1,1)
scatter(centroid.stim(1,:),centroid.stim(2,:), [], rad.stim); colorbar
ind = find(~isnan(centroid.stim(1,:)));
centroid.med = findMaxNeighbors(centroid.stim(:,ind),2); %function to f
hold on;
scatter(centroid.med(1),centroid.med(2),'or')
centroid.dist = sqrt((centroid.stim(1,:)-centroid.med(1)).^2 + (centroid.stim(2,:)-centroid.med(2)).^2);
title('Color- radius')
xlabel('x-pos')
ylabel('y-pos')
subplot(2,1,2)
hist(centroid.dist,0:0.5:60)
title([num2str(sum(centroid.dist<4)) ' trials w/in 4 deg'])
sgtitle([num2str(sum(~isnan(centroid.dist))) '/' num2str(nTrials) ' measurable trials'])
xlabel('Centroid distance from median')
