function spikesByTrial = get_raster_JL(spike_times,trigger_times,varargin)

p = inputParser;

p.addParamValue('stop',[],@isnumeric);
p.addParamValue('start',[],@isnumeric);
p.addParamValue('plot',false,@islogical);
p.addParamValue('count',false,@islogical);
p.addParamValue('color',[0 0 0],@isnumeric);
p.addParamValue('trialOffset',0,@isnumeric);
% parse inputs
p.parse(varargin{:});
params = p.Results;

if isempty(params.stop)
    params.stop = mean(diff(trigger_times));
end

if isempty(params.start)
    params.start = 0;
end

for rep = 1:numel(trigger_times)
     tempSpikes = spike_times(spike_times>params.start+trigger_times(rep) & spike_times < params.stop + trigger_times(rep)) - (trigger_times(rep));
      spikesByTrial{rep} = tempSpikes;
end

if params.plot
%     figure;
    for rep = 1:numel(trigger_times)
    hold on
    if ~isempty(spikesByTrial{rep})
        try
        plot([spikesByTrial{rep}; spikesByTrial{rep}],[rep+params.trialOffset-0.35; rep+params.trialOffset+0.35],'Color',params.color);
        catch
            plot([spikesByTrial{rep}'; spikesByTrial{rep}'],[rep+params.trialOffset-0.35; rep+params.trialOffset+0.35],'Color',params.color);
        end
    end
    end
end

if params.count
    spikesByTrial = cell2mat(cellfun(@numel,spikesByTrial,'un',0));
end
            
end
            