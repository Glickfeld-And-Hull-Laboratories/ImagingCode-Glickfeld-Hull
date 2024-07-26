
function [Neuron] = get_raster_NB(spike_times,trigger_timesBL,trigger_timesSD,trigger_timesSS,BLstart,BLstop,SDstart,SDstop,SSstart,SSstop)
p = inputParser;

p.addParamValue('stop',[],@isnumeric);
p.addParamValue('start',[],@isnumeric);
p.addParamValue('plot',false,@islogical);
p.addParamValue('count',false,@islogical);
p.addParamValue('color',[0 0 0],@isnumeric);
p.addParamValue('trialOffset',0,@isnumeric);
% parse inputs
% % p.parse(varargin{:});
params = p.Results;

params.BLstart=BLstart;
params.BLstop=BLstop;
params.SDstart=SDstart;
params.SDstop=SDstop;
params.SSstart=SSstart;
params.SSstop=SSstop;

% if isempty(params.stop)
%     params.stop = mean(diff(trigger_times));
% end
% 
% if isempty(params.start)
%     params.start = 0;
% end

HzCalc=[];
%%
for j=1:size(spike_times,1)
    spike_times_temp=cell2mat(spike_times(j,1));
for i=1:numel(trigger_timesBL)       
    tempSpikes=spike_times_temp(spike_times_temp>params.BLstart+trigger_timesBL(i) & spike_times_temp<params.BLstop+trigger_timesBL(i))-trigger_timesBL(i);
    TF=isempty(tempSpikes);
    if TF==1
        HzCount=0;
    elseif TF==0
        HzCount=length(tempSpikes);
    end
        HzCalc=HzCount/params.BLstop;
        Neuron(j).BLResp(i).Stim=tempSpikes;
        % Neuron(j).BL4Raster(i).Stim=tempSpikes-params.BLstop;
        Neuron(j).BLHzCount(i).Stim=HzCount;
        Neuron(j).BLHz(i).Stim=HzCalc;
        % Neuron(j).Hz(i).Stim=HzCount;
        % spikesByTrial(i) = tempSpikes;
    end
end
%%
for j=1:size(spike_times,1)
    spike_times_temp=cell2mat(spike_times(j,1));
for i=1:numel(trigger_timesSD)       
    tempSpikes=spike_times_temp(spike_times_temp>params.SDstart+trigger_timesSD(i) & spike_times_temp<params.SDstop+trigger_timesSD(i))-trigger_timesSD(i);
    TF=isempty(tempSpikes);
    if TF==1
        HzCount=0;
    elseif TF==0
        HzCount=length(tempSpikes);
    end
        HzCalc=HzCount/params.SDstop;
        Neuron(j).SDResp(i).Stim=tempSpikes;
        Neuron(j).SDHzCount(i).Stim=HzCount;
        Neuron(j).SDHz(i).Stim=HzCalc;
        % Neuron(j).Hz(i).Stim=HzCount;
        % spikesByTrial(i) = tempSpikes;
    end
end
%%
for j=1:size(spike_times,1)
    spike_times_temp=cell2mat(spike_times(j,1));
for i=1:numel(trigger_timesSS)       
    tempSpikes=spike_times_temp(spike_times_temp>params.SSstart+trigger_timesSS(i) & spike_times_temp<params.SSstop+trigger_timesSS(i))-trigger_timesSS(i);
    TF=isempty(tempSpikes);
    if TF==1
        HzCount=0;
    elseif TF==0
        HzCount=length(tempSpikes);
    end
        HzCalc=HzCount/params.SSstop;
        Neuron(j).SSResp(i).Stim=tempSpikes;
        Neuron(j).SSHzCount(i).Stim=HzCount;
        Neuron(j).SSHz(i).Stim=HzCalc;
        % Neuron(j).Hz(i).Stim=HzCount;
        % spikesByTrial(i) = tempSpikes;
    end
end

%%

% % for i=1:size(spikes.spiketime_S,1)
% spike_times=cell2mat(spikes.spiketime_S(i,1));
% for rep = 1:numel(trigger_times)
%      tempSpikes = spike_times(spike_times>params.start+trigger_times(rep) & spike_times < params.stop + trigger_times(rep)) - (trigger_times(rep));
%      Neuron(j).RSP(rep,:) = tempSpikes;
%      spikesByTrial{rep} = tempSpikes;
% end
% end
     % spiketime = tempSpikes;

% if params.plot
%     % figure;
%     for rep = 1:numel(trigger_times)
%     hold on
%     if ~isempty(spikesByTrial{rep})
%         try
%         plot([spikesByTrial{rep}; spikesByTrial{rep}],[rep+params.trialOffset-0.35; rep+params.trialOffset+0.35],'Color',params.color);
%         catch
%             plot([spikesByTrial{rep}'; spikesByTrial{rep}'],[rep+params.trialOffset-0.35; rep+params.trialOffset+0.35],'Color',params.color);
%         end
%     end
%     end
% end
% 
% % if params.count
% %     spikesByTrial = cell2mat(cellfun(@numel,spikesByTrial,'un',0));
% % end
% 
% end
% 