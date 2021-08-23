function [stimOnFrames stimOffFrames] = photoFrameFinder_Sanworks(events);
%new code for Sanworks photodiode
%events is info.frames from scanbox .mat file
    %fix discontinuity
    if find(diff(events)<0)
        events(find(diff(events)<0)+1:end) = events(find(diff(events)<0)+1:end)+2^16;
    end
        %remove any 0s from frame list
    if find(events == 0)
        events(find(events==0)) = [];
    end
    framesOn = unique(events);
    nframes = max(framesOn(:));
    frameMat = zeros(1,nframes);
    frameMat(framesOn) = 1;
    frameTrig = diff(frameMat);
    stimOnFrames = find(frameTrig == 1);
    stimOffFrames = find(frameTrig == -1);
end
    