% datename=input('Date?')
% load Z:\home\Naomi\Electrophysiology\Kilosort_Analysis\,datename,'\',[datename,'_datafile',spikeFiles,'_spikes.mat']

MUITrial=find(spikes.MUA)';
SUTrial=find(spikes.MUA==0);
%%
for i=1:length(MUITrial)
    MUNeuron.Neuron(i).Stim=Neuron.Neuron(MUITrial(i)).Stim
end
for i=1:length(SUTrial)
    SUNeuron.Neuron(i).Stim=Neuron.Neuron(SUTrial(i)).Stim
end
%%
% for i=1:length(MUITrial)
% MUSpikes.Neuron(i,1)=spikes.spiketime_S(MUITrial(i))
% end
% 
% for i=1:length(SUTrial)
% SUSpikes.Neuron(i,1)=spikes.spiketime_S(SUTrial(i))
% end
% %%
% %% get broad spont FR over time 
% blwindow=input('Baseline window? (in seconds)');
% sdwindow=input('Peristim window? (in seconds)');
% sswindow=input('Peri-second-stim window? (in seconds)');
% % % for spontaneous activity
% % stimOn=stimOn(1:10,1);
% trigger_timesBL=(stimOn(:,1)/data.info.sampleFreq)-blwindow;
% trigger_timesSD=(stimOn(:,1)/data.info.sampleFreq);
% trigger_timesSS=(stimOn(:,2)/data.info.sampleFreq);
% %%use whole trial for rasters
% [MUwholetrial]=get_raster_NB(MUSpikes.Neuron,trigger_timesBL,trigger_timesSD,trigger_timesSS,0,blwindow,0,1000,0,0);
% [MUNeuronSDR]=get_raster_NB(MUSpikes.Neuron,trigger_timesBL,trigger_timesSD,trigger_timesSS,0,blwindow,0,sdwindow,0,sswindow);
% [SIwholetrial]=get_raster_NB(MUSpikes.Neuron,trigger_timesBL,trigger_timesSD,trigger_timesSS,0,blwindow,0,1000,0,0);
% [SINeuronSDR]=get_raster_NB(MUSpikes.Neuron,trigger_timesBL,trigger_timesSD,trigger_timesSS,0,blwindow,0,sdwindow,0,sswindow);
% %%
% 
% dimensions=input('SF(1) or SF and DIR(2)?')
% [StimTrialTimes,Combos,Trials,StimIdent] = StimIdentities(dimensions,inputs,trigger_timesSD,trigger_timesBL);
% Trials=Trials(1:length(trigger_timesSD),:);
% [MUNeuron,Neuron4Raster] = TrialTypeStruct(StimIdent,Trials,NeuronSDR,wholetrial);
% for i=1:length(MUITrial)
%     MUITrialType(i,1)=Trials(MoveTrial(i),4);
%     trigger_timesMVSD(i)=trigger_timesSD(MoveTrial(i));
%     trigger_timesMVSS(i)=trigger_timesSS(MoveTrial(i));
% end
% % trigger_timesMVBL=trigger_timesMVSD-blwindow;
% % 
% % [NeuronMoving]=get_raster_NB(spikes.spiketime_S,trigger_timesMVBL,trigger_timesMVSD,trigger_timesMVSS,0,blwindow,0,sdwindow,0,sswindow);
% % [wholemovingtrial]=get_raster_NB(spikes.spiketime_S,trigger_timesMVBL,trigger_timesMVSD,trigger_timesMVSS,0,blwindow,0,1000,0,0);
% % %%
% % 
% % StatTriTimes1=trigger_timesSD;
% % StatTriTimes2=trigger_timesSS;
% % for i=1:length(MoveTrial);
% %       StatTriTimes1(MoveTrial(i))=NaN;
% %       StatTriTimes2(MoveTrial(i))=NaN;
% % end
% % StatTimeSD=StatTriTimes1(~isnan(StatTriTimes1));
% % StatTimeSS=StatTriTimes2(~isnan(StatTriTimes2));
% % %%
% % 
% % trigger_timesSTSD=StatTimeSD;
% % trigger_timesSTSS=StatTimeSS;
% % trigger_timesSTBL=StatTimeSD-blwindow;
% % 
% % [NeuronStat]=get_raster_NB(spikes.spiketime_S,trigger_timesSTBL,trigger_timesSTSD,trigger_timesSTSS,0,blwindow,0,sdwindow,0,sswindow);
% % [wholestattrial]=get_raster_NB(spikes.spiketime_S,trigger_timesSTBL,trigger_timesSTSD,trigger_timesSTSS,0,blwindow,0,1000,0,0);
% % end