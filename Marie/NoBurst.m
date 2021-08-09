function noburstSpikes = NoBurst(unit, time)
k=1;
for i = 2:(length(unit)-1)
    if (((unit(i) - unit(i-1)) > time) && ((unit(i+1) - unit(i)) > time))
        noburstSpikes(k,:)= unit(i);
        k = k+1;
    end
end
end
    