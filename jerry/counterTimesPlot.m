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
data_tc = nan(nOn+nOff,nCells,nTrials);
for itrial = 1:nTrials
  if ~isnan(cStimOn(itrial)) & (cStimOn(itrial)+nOn+nOff/2)<nFrames
    data_tc(:,:,itrial) = npSub_tc(cStimOn(itrial)-nOff/2:cStimOn(itrial)-1+nOn+nOff/2,:);
  end
end


%troubleshoot
tavg = squeeze(mean(data_tc_trial, 2));
for i = 1:107
    plot(tavg(i,:));
end
    



[cStimOn stimOffs] = photoFrameFinder_Sanworks(info.frame);
nTrials = length(cStimOn);
[nFrames nCells] = size(npSub_tc);
nOn = input.nScansOn;
nOff = input.nScansOff;
data_tc = nan(nOn+nOff,nCells,nTrials);

for itrial = 1:nTrials
  if ~isnan(cStimOn(itrial)) & (cStimOn(itrial)+nOn+nOff/2)<nFrames
    data_tc(:,:,itrial) = npSub_tc(cStimOn(itrial)-nOff/2:cStimOn(itrial)-1+nOn+nOff/2,:);
  end
end
%data_tc = data_tc(:,:,1:end-1);
data_f_trial = mean(data_tc(1:nOff/2,:,:),1);
data_dfof_trial = bsxfun(@rdivide, bsxfun(@minus,data_tc, data_f_trial), data_f_trial);

tavg = squeeze(nanmean(data_dfof_trial,3));
figure;
hold on
for icell = 1:nCells
    plot(tavg(:,icell), 'LineWidth',.005,'color',[.25 .25 .25]);
end
tc_cellavg = mean(tavg,2);
plot(tc_cellavg);