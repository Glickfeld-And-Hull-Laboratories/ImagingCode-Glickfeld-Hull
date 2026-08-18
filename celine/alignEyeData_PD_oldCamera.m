function [rad centroid] = alignEyeData_PD(Eye_data,input,stimOns);
% Eye_data is output of extractEyeData
% input is from mworks

% Output: rad and centroid are structures with timecourses and base/stim 
% for pupil radius and centroid. 
% centroid also has a field (dists) that measures the distance from the 
% median pupil position on each trial during the stim.

calib = 1/26.6; %mm per pixel
sesh_id = input.saveTime;
msg1 = ['Session ' sesh_id 'used photodiode stimOns as cStimOn'];
msgbox(msg1);

if nargin == 3 
    cStimOn = stimOns;
    if isfield(input,'nScansOn') % this is true normally for visstimret
        prewin_frames = input.nScansOn;
        postwin_frames = input.nScansOn;
    else
        if isfield(input,'trialOutcomeCell') % need to double check for behavior
            prewin_frames = input.frameRateHz;
            postwin_frames = input.frameRateHz.*3;
        else
            prewin_frames = input.nFramesOn;
            postwin_frames = input.nFramesOn;
        end
    end
else
    if isfield(input,'nScansOn')
        cStimOn = input.nScansOff+1:input.nScansOff+input.nScansOn:input.counterValues{end}(end);
        prewin_frames = input.nScansOn;
        postwin_frames = input.nScansOn;
    elseif isfield(input,'cStimOn')
        cStimOn = celleqel2mat_padded(input.cStimOn);
        if isfield(input,'trialOutcomeCell')
            prewin_frames = input.frameRateHz;
            postwin_frames = input.frameRateHz.*3;
        else
            prewin_frames = input.nFramesOn;
            postwin_frames = input.nFramesOn;
        end
    elseif isfield(input,'cStimOneOn')
        cStimOn = celleqel2mat_padded(input.cStimOneOn);
        if isfield(input,'nStimOneFramesOn')
            prewin_frames = unique(celleqel2mat_padded(input.nStimOneFramesOn));
            postwin_frames = unique(celleqel2mat_padded(input.nStimOneFramesOn));
        else
            if length(input.frameRateHz)>1
                input.frameRateHz = unique(input.frameRateHz);
            end
            prewin_frames = input.frameRateHz-1;
            postwin_frames = input.frameRateHz;
        end
    end
end
nTrials = size(cStimOn,2);
Rad_temp = sqrt(Eye_data.Area./pi);
Centroid_temp = Eye_data.Centroid;
Rad_temp(Eye_data.badFrames,:) =nan(length(Eye_data.badFrames),1);
Centroid_temp(Eye_data.badFrames,:) = nan(length(Eye_data.badFrames),2);
rad.tc = zeros(prewin_frames+postwin_frames, nTrials);
centroid.tc = zeros(prewin_frames+postwin_frames,2, nTrials);
nframes = size(Rad_temp,1);

for itrial = 1:nTrials
    if itrial == nTrials
        crange = [double(cStimOn(itrial))-prewin_frames:nframes];
    else
        crange = [double(cStimOn(itrial))-prewin_frames: double(cStimOn(itrial+1)-prewin_frames-1)];
    end
    if sum(isnan(Rad_temp(crange)),1)>0
        if sum(isnan(Rad_temp(crange)),1)./length(crange)> 0.25
            Rad_temp(crange,1) = NaN(length(crange),1);
            Centroid_temp(crange,:) = NaN(length(crange),2);
        else
            nanind = intersect(crange,find(isnan(Rad_temp)));
            dataind = intersect(crange,find(~isnan(Rad_temp)));
            for inan = 1:length(nanind)
                gap = min(abs(nanind(inan)-dataind),[],1);
                good_ind_stim = find(abs(nanind(inan)-dataind) == gap);
                Rad_temp(nanind(inan),1) = mean(Rad_temp(dataind(good_ind_stim),1),1);
                Centroid_temp(nanind(inan),:) = mean(Centroid_temp(dataind(good_ind_stim),:),1);
            end
        end
    end
    if itrial < nTrials
        rad.tc(:,itrial) = Rad_temp(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
        centroid.tc(:,:,itrial) = Centroid_temp(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
    else
        if (cStimOn(itrial)+postwin_frames)<nframes
            rad.tc(:,itrial) = Rad_temp(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
            centroid.tc(:,:,itrial) = Centroid_temp(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
        else
            rad.tc(:,itrial) = nan(prewin_frames+postwin_frames,1);
            centroid.tc(:,:,itrial) = nan(prewin_frames+postwin_frames,2,1);
        end
    end
        
end

rad_mat_calib = bsxfun(@times, rad.tc, calib);
centroid_mat_calib = bsxfun(@times,centroid.tc,calib);
rad.base = mean(rad_mat_calib(1:prewin_frames,:),1);
rad.stim = mean(rad_mat_calib(prewin_frames+1:end,:),1);
centroid.base = squeeze(mean(centroid_mat_calib(1:prewin_frames,:,:),1))./0.025;
centroid.stim = squeeze(mean(centroid_mat_calib(prewin_frames+1:end,:,:),1))./0.025;

figure; subplot(2,1,1)
scatter(centroid.stim(1,:),centroid.stim(2,:), [], rad.stim); colorbar
ind = find(~isnan(centroid.stim(1,:)));
%centroid_med = geometric_median(centroid_stim(:,ind));
centroid.med = findMaxNeighbors(centroid.stim(:,ind),2);
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
%print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupilPosDist.pdf']),'-dpdf','-fillpage');
movegui('center')
