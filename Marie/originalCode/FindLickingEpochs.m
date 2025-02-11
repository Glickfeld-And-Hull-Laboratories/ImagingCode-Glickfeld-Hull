function FirstLicksEpochs = FindLickingEpochs(licks)
FirstLicksEpochs = [];
k = 1;
for i = 2:(length(licks)-2)
    if (((licks(i) - licks(i-1)) > 1) || isnan(licks(i) - licks(i-1)))     % sets 1 as "reset" after which I would consider a potential new onset of licking
        if (((licks(i+1)-licks(i)) < .21) && ((licks(i+2)-licks(i+1))<0.21)) % sets boundary between licks to consider as onset of epoch
            FirstLicksEpochs(k) = licks(i);
            k = k+1;
        end
    end
end
FirstLicksEpochs = FirstLicksEpochs.';
end