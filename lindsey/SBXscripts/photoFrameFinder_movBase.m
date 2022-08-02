function [photoLoc, photoFrames] = photoFrameFinder_movBase(photoData,nframes)

    % uses ephys output from scanbox to find frames stim onsets

    %Clock
    clock = photoData(1:2:end);
    clockSize = size(clock);
    clockAvg = mean(clock);
    clockStd = std(clock);
    clockAvgStd = clockAvg-4*clockStd;

    %Photodiode
    photo = photoData(2:2:end);

    clockLoc = [];
    for i = 2:clockSize(1)
        if(clock(i) < clockAvgStd && clock(i-1) > clockAvgStd)
            clockLoc = [clockLoc i];
        end
    end

    photo_smooth = smoothdata(photo,'gaussian',500);
    photo_high = movmax(photo_smooth, 100000);
    photo_smooth = photo_smooth-photo_high;
    photo_diff = diff(photo_smooth);
    photo_diff_rect = photo_diff;
    photo_diff_rect(find(photo_diff>0)) = 0;
    photoSize = size(photo_diff);

    photoMin = min(photo_diff_rect(clockLoc(1):photoSize(1)),[],1);
    photoThresh = 0.5.*photoMin;
    photoRecover = 0.05.*photoMin;
    photoLoc = [];
    photoFrames = [];
    n = 0;
    for i = 50:nframes
        photoAmp_min = min(photo_diff_rect(clockLoc(i-1):clockLoc(i)),[],1);
        photoAmp_prev = photo_diff_rect(clockLoc(i-1));
        if n == 0 && photoAmp_min<photoThresh/2 && photoAmp_prev>photoThresh/2
            photoLoc = [photoLoc clockLoc(i)];
            photoFrames = [photoFrames i];
            n=1+n;
        elseif n>0 && max(photo_diff_rect(photoLoc(n):clockLoc(i-1)),[],1) > photoRecover && i-photoFrames(n) > 4
            if photoAmp_min<photoThresh && photoAmp_prev>photoThresh 
                photoLoc = [photoLoc clockLoc(i)];
                photoFrames = [photoFrames i];
                n = n+1;
            end
        end
        if rem(i,1000) == 0
            fprintf([num2str(i) '\n'])
        end
    end
    diff_photo = diff(photoFrames);
    ind = find(diff_photo>35);
    ind_double = ind(find(diff(ind)==1));
    photoFrames(ind_double+1) = [];
    photoLoc(ind_double+1) = [];
end