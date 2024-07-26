 function [position,template,chID] = WFpositioner(avgwaveform,channelposlist)
        % position [x,y] is the average electrode position weighted by the wave
        % size
       
        wavesize = peak2peak(avgwaveform); % Measure the peak to peak size of the wave on each channel
        wavesize = wavesize.^10; % Bias the wavesizes so that big waves are weighted way stronger than small
        wavesize = bsxfun(@rdivide,wavesize,sum(wavesize)); % Normalize the wavesizes to add up to 1
        [~,chID] = max(wavesize);
        template = avgwaveform(:,chID);
        weightedx = nansum(wavesize'.*channelposlist(:,1)); % Weighted average position.
        weightedy = nansum(wavesize'.*channelposlist(:,2));
        position = [weightedx, weightedy];
 end