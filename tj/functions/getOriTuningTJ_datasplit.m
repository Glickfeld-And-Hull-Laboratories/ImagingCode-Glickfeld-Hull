function [avgResponseEaOri_1,semResponseEaOri_1,vonMisesFitAllCells_1,fitReliability_1,R_square_1,tuningTC_1, ...
    avgResponseEaOri_2,semResponseEaOri_2,vonMisesFitAllCells_2,fitReliability_2,R_square_2,tuningTC_2] = ...
    getOriTuningTJ_datasplit(tc,mworks,downSampleFactor,split_idx)


nBoot = 1;


[nFrames,nCells] = size(tc);

%params
nOn = (mworks.nScansOn)./downSampleFactor;
nOff = (mworks.nScansOff)./downSampleFactor;
basewin = nOff/2:nOff;
respwin = nOff+1:nOn+nOff;
nTrials = mworks.trialsSinceReset;
if mod(nFrames,nTrials) > 0
    nframesAllTrials = nTrials*(nOn+nOff);
    if nframesAllTrials > nFrames
        nTrials = floor(nFrames/(nOn+nOff));
    end 
end
tc = tc(1:(nOn+nOff)*nTrials,:);
tDirection = cell2mat_padded(mworks.tGratingDirectionDeg);
tDirection = tDirection(1:nTrials);
tOrientation = tDirection;
tOrientation(tOrientation > 179) = tOrientation(tOrientation > 179) - 180;
[orientationInd, orientations] = findgroups(tOrientation);
nStim = length(orientations);
% theta = [orientations 180];

trialTC = reshape(tc,nOff+nOn,nTrials,nCells);

%*********
% [trialTC_1 trialTC_1_idx] = datasample(trialTC, length(trialTC)/2, 2 ,Replace=false);
% trialTC_1_oris = orientationInd(trialTC_1_idx);

trialTC_1 = trialTC(:, split_idx(1:round(length(split_idx)/2)), :);
trialTC_1_oris = orientationInd(:, split_idx(1:round(length(split_idx)/2)), :);

trialTC_2 = trialTC(:, split_idx(round(length(split_idx)/2)+1:end), :);
trialTC_2_oris = orientationInd(:, split_idx(round(length(split_idx)/2)+1:end), :);

% aaa_oris_1st_set = aaa_oris(:,1:length(aaa_oris)/2);
% aaa_2nd_set = aaa(:,(length(aaa)/2)+1:length(aaa));
% aaa_oris_2nd_set = aaa_oris(:,(length(aaa_oris)/2)+1:length(aaa_oris));
% 
% trialTC_2_idx = setxor(trialTC_1_idx, double(cell2mat(mworks.counter)/90));
% trialTC_2 = trialTC(:,trialTC_2_idx,:);
% trialTC_2_oris = orientationInd(trialTC_2_idx);
%*********
%'half' 1
F_1 = mean(trialTC_1(basewin,:,:),1);
dFF_1 = bsxfun(@rdivide, bsxfun(@minus,trialTC_1,F_1),F_1);
dFF_1(find(isnan(dFF_1)))=0;

tuningTC_1 = nan(nOn+nOff,nCells,nStim);
avgResponseEaOri_1 = nan(nCells,nStim);
semResponseEaOri_1 = nan(nCells,nStim);
tuningResamp_1 = nan(nCells,nStim,nBoot);

for istim = 1:nStim
    ind = find(trialTC_1_oris == istim);
    tuningTC_1(:,:,istim) = squeeze(mean(dFF_1(:,ind,:),2));
    avgResponseEaOri_1(:,istim) = squeeze(mean(mean(dFF_1(respwin,ind,:),1),2));
    semResponseEaOri_1(:,istim) = squeeze(ste(mean(dFF_1(respwin,ind,:),1),2));
    for iboot = 1:nBoot
        n = length(ind);
        randTrials = randsample(ind,n,1);
        tuningResamp_1(:,istim,iboot) = squeeze(mean(mean(...
            dFF_1(respwin,randTrials,:),1),2));
    end
end

[vonMisesFitAllCells_1,~,fitReliability_1,R_square_1] = vonmisesReliableFit(avgResponseEaOri_1,...
    tuningResamp_1,orientations,nBoot);

%'half' 2
F_2 = mean(trialTC_2(basewin,:,:),1);
dFF_2 = bsxfun(@rdivide, bsxfun(@minus,trialTC_2,F_2),F_2);
dFF_2(find(isnan(dFF_2)))=0;

tuningTC_2 = nan(nOn+nOff,nCells,nStim);
avgResponseEaOri_2 = nan(nCells,nStim);
semResponseEaOri_2 = nan(nCells,nStim);
tuningResamp_2 = nan(nCells,nStim,nBoot);

for istim = 1:nStim
    ind = find(trialTC_2_oris == istim);
    tuningTC_2(:,:,istim) = squeeze(mean(dFF_2(:,ind,:),2));
    avgResponseEaOri_2(:,istim) = squeeze(mean(mean(dFF_2(respwin,ind,:),1),2));
    semResponseEaOri_2(:,istim) = squeeze(ste(mean(dFF_2(respwin,ind,:),1),2));
    for iboot = 1:nBoot
        n = length(ind);
        randTrials = randsample(ind,n,1);
        tuningResamp_2(:,istim,iboot) = squeeze(mean(mean(...
            dFF_2(respwin,randTrials,:),1),2));
    end
end

[vonMisesFitAllCells_2,~,fitReliability_2,R_square_2] = vonmisesReliableFit(avgResponseEaOri_2,...
    tuningResamp_2,orientations,nBoot);


end