function [isMovingT,Trials,wholemovingtrial,NeuronMoving,wholestattrial,NeuronStat] = GetMovingStructs(expt_data,Trials,trigger_timesSD,trigger_timesSS,blwindow,sdwindow,sswindow,spikes)

beforeStimOnWheel = 400; % time window for before stim on
afterStimOnWheel = 100; % time window for after stim on

wheelSI = cellfun(@(x) floor(mean(diff(x))/1000),expt_data.wheelSpeedTimesUs,'un',0);
temp_LHS = cellfun(@(y,z,s,start) floor(find((y/1000-z)>start,1,'first') - beforeStimOnWheel/s),expt_data.wheelSpeedTimesUs,num2cell(expt_data.tThisTrialStartTimeMs),wheelSI,num2cell(expt_data.tItiTimeMs),'un',0);
stimWindow_wheelPts = cellfun(@(x,s) x:(x + (beforeStimOnWheel+afterStimOnWheel)/s),temp_LHS,wheelSI,'un',0);
try
    isMovingT = cell2mat(cellfun(@(x,y,z) mean(abs(x(y)))>=(500/z),expt_data.wheelSpeedValues,stimWindow_wheelPts,wheelSI,'un',0));
catch
    isMovingT = cell2mat(cellfun(@(x,y,z) mean(abs(x(y(1):end)))>=(500/z),expt_data.wheelSpeedValues,stimWindow_wheelPts,wheelSI,'un',0));
end
disp([num2str(sum(isMovingT)) ' moving trials']);
Trials(:,5)=isMovingT';


MoveTrial=find(isMovingT)';
for i=1:length(MoveTrial)
    MoveTrialType(i,1)=Trials(MoveTrial(i),4);
    trigger_timesMVSD(i)=trigger_timesSD(MoveTrial(i));
    trigger_timesMVSS(i)=trigger_timesSS(MoveTrial(i));
end
trigger_timesMVBL=trigger_timesMVSD-blwindow;

[NeuronMoving]=get_raster_NB(spikes.spiketime_S,trigger_timesMVBL,trigger_timesMVSD,trigger_timesMVSS,0,blwindow,0,sdwindow,0,sswindow);
[wholemovingtrial]=get_raster_NB(spikes.spiketime_S,trigger_timesMVBL,trigger_timesMVSD,trigger_timesMVSS,0,blwindow,0,1000,0,0);
%%

StatTriTimes1=trigger_timesSD;
StatTriTimes2=trigger_timesSS;
for i=1:length(MoveTrial);
      StatTriTimes1(MoveTrial(i))=NaN;
      StatTriTimes2(MoveTrial(i))=NaN;
end
StatTimeSD=StatTriTimes1(~isnan(StatTriTimes1));
StatTimeSS=StatTriTimes2(~isnan(StatTriTimes2));
%%

trigger_timesSTSD=StatTimeSD;
trigger_timesSTSS=StatTimeSS;
trigger_timesSTBL=StatTimeSD-blwindow;

[NeuronStat]=get_raster_NB(spikes.spiketime_S,trigger_timesSTBL,trigger_timesSTSD,trigger_timesSTSS,0,blwindow,0,sdwindow,0,sswindow);
[wholestattrial]=get_raster_NB(spikes.spiketime_S,trigger_timesSTBL,trigger_timesSTSD,trigger_timesSTSS,0,blwindow,0,1000,0,0);
end