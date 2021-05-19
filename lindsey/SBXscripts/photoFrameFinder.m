function [photoLoc photoFrames] = photoFrameFinder(photoData,nframes);

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

    photo_smooth = smooth(photo,500);
    photo_high = movmax(photo_smooth, 100000);
    photo_smooth = photo_smooth-photo_high;
    photoSize = size(photo_smooth);
    photoMax = max(photo_smooth(clockLoc(1):photoSize(1)),[],1);
    photoMin = min(photo_smooth(clockLoc(1):photoSize(1)),[],1);
    photoAvg50= photoMax-(0.5*(photoMax-photoMin));
    photoAvg25 = photoMax-(0.25*(photoMax-photoMin));
    photoAvg75 = photoMax-(0.75*(photoMax-photoMin));
    photoAvg5 = photoMax-(0.05*(photoMax-photoMin));

    photoLoc = [];
    photoFrames = [];
    n = 0;
    for i = 2:nframes
        photoAmp_min = min(photo_smooth(clockLoc(i-1):clockLoc(i)),[],1);
        photoAmp_prev = photo_smooth(clockLoc(i-1));
        
        if n == 0 && photoAmp_min<photoAvg50 && photoAmp_prev>photoAvg50
            ind = find(photo_smooth(clockLoc(i):clockLoc(i+1))<photoAvg50);
            if length(ind)>100
                photoLoc = [photoLoc clockLoc(i)];
                photoFrames = [photoFrames i];
                n=1+n;
            end
        elseif n>0 && max(photo_smooth(photoLoc(n):clockLoc(i-1)),[],1) > photoAvg5 && i-photoFrames(n) > 4
            if photoAmp_min<photoAvg25 && photoAmp_prev>photoAvg25 
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
    ind_double_temp = ind_double;
    for i = 1:length(ind_double)
        if diff_photo(ind_double_temp(i))>25
            photoFrames(ind_double_temp(i)+1) = [];
            photoLoc(ind_double_temp(i)+1) = [];
            ind_double_temp(i+1:end) = ind_double(i+1:end)-1;
        elseif diff_photo(ind_double_temp(i))<25
            photoFrames = [photoFrames(1:ind_double_temp(i)) photoFrames(ind_double_temp(i))+11 photoFrames(ind_double_temp(i)+1:end)];
            photoLoc = [photoLoc(1:ind_double_temp(i)) photoLoc(ind_double_temp(i))+3665 photoLoc(ind_double_temp(i)+1:end)];
            ind_double_temp(i+1:end) = ind_double(i+1:end)+1;
        end
    end
end