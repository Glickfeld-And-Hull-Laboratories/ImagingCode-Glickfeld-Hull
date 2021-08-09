function [stimOnFrames] = photoFrameFinder_Sanworks(events);
%new code for Sanworks photodiode
%events is info.frames from scanbox .mat file
    framesOn = unique(events);
    nframes = max(framesOn(:));
    frameMat = zeros(1,nframes);
    frameMat(framesOn) = 1;
    frameTrig = diff(frameMat);
    stimOnFrames = find(frameTrig == 1);
end
    