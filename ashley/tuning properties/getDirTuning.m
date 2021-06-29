function [avgResponseEaDir,semResponseEaDir,vonMisesFitAllCells,fitReliability,R_square,tuningTC] = ...
    getDirTuning(tc,mworks,downSampleFactor)

nBoot = 1000;


[nFrames,nCells] = size(tc);

%params
nOn = (mworks.nScansOn)./downSampleFactor;
nOff = (mworks.nScansOff)./downSampleFactor;
basewin = nOff/2:nOff;
respwin = nOff+1:nOn+nOff;
nTrials = mworks.trialSinceReset;
if mod(nFrames,nTrials) > 0
    nframesAllTrials = nTrials*(nOn+nOff);
    if nframesAllTrials > nFrames
        nTrials = floor(nFrames/(nOn+nOff));
    end 
end
tc = tc(1:(nOn+nOff)*nTrials,:);
tDirection = cell2mat_padded(mworks.tGratingDirectionDeg);
tDirection = tDirection(1:nTrials);
[directionInd, directions] = findgroups(tDirection);
nStim = length(directions);

trialTC = reshape(tc,nOff+nOn,nTrials,nCells);
F = mean(trialTC(basewin,:,:),1);
dFF = bsxfun(@rdivide, bsxfun(@minus,trialTC,F),F);

tuningTC = nan(nOn+nOff,nCells,nStim);
avgResponseEaDir = nan(nCells,nStim);
semResponseEaDir = nan(nCells,nStim);
tuningResamp = nan(nCells,nStim,nBoot);
for istim = 1:nStim
    ind = find(directionInd == istim);
    tuningTC(:,:,istim) = squeeze(mean(dFF(:,ind,:),2));
    avgResponseEaDir(:,istim) = squeeze(mean(mean(dFF(respwin,ind,:),1),2));
    semResponseEaDir(:,istim) = squeeze(ste(mean(dFF(respwin,ind,:),1),2));
    for iboot = 1:nBoot
        n = length(ind);
        randTrials = randsample(ind,n,1);
        tuningResamp(:,istim,iboot) = squeeze(mean(mean(...
            dFF(respwin,randTrials,:),1),2));
    end
end

[vonMisesFitAllCells,~,fitReliability,R_square] = vonmisesReliableFit_direction(avgResponseEaDir,...
    tuningResamp,directions,nBoot);

end