nTrials = size(input.counterTimesUs,2);
counterTimes = [];
for i = 1:nTrials
    counterTimes = [counterTimes input.counterTimesUs{i}];
end
plot(diff(counterTimes))

%timecourse alignment
[cStimOn stimOffs] = photoFrameFinder_Sanworks(info.frame);
nTrials = length(cStimOn);
[nFrames nCells] = size(npSub_tc);
nOn = input.nScansOn;
nOff = input.nScansOff;
data_tc = nan(nOn+nOff,nTrials,nCells);
for itrial = 1:nTrials
  if ~isnan(cStimOn(itrial)) & (cStimOn(itrial)+nOn+nOff/2)<nFrames
    data_tc(:,itrial,:) = npSub_tc(cStimOn(itrial)-nOff/2:cStimOn(itrial)-1+nOn+nOff/2,:);
  end
end


%troubleshoot
tavg = squeeze(mean(data_tc_trial, 2));
for i = 1:107
    plot(tavg(i,:));
end
    