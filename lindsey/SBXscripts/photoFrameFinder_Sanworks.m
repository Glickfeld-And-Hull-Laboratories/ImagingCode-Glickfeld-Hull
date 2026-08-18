function [stimOnFrames stimOffFrames] = photoFrameFinder_Sanworks(events);
%new code for Sanworks photodiode
%events is info.frame from scanbox .mat file
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
    %correct for single frame skip (on-off-on or off-on-off blip)
    dFrameTrig = diff(frameTrig);
    shortEvents = find(abs(dFrameTrig)==2);
    msg = [num2str(length(shortEvents)) ' short frame on/off events detected and corrected.'];
    if ~isempty(shortEvents)
        msgbox(msg,'photoFrameFinder Message','warn');
        for f = 1:length(shortEvents)
            replaceFrames = [shortEvents(f) shortEvents(f)+1];
            frameTrig(replaceFrames) = 0;
        end
    end
    stimOnFrames = find(frameTrig == 1) + 1; % +1 accounts for derivative
    stimOffFrames = find(frameTrig == -1) + 1;
end
    