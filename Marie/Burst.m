function [burstSpikes, k] = Burst(unit, time)
k=1;
for i = 1:length(unit)-3
    if (((unit(i+1) - unit(i)) <= time)) %&& ((unit(i+2) - unit(i+1)) < time))
        unit(i);
        unit(i+1);
        unit(i+2);
        burstSpikes(k,:)= unit(i);
        burstSpikes(k+1, :) = unit(i+1);
        burstSpikes(k+2, :) = unit(k+2);
        k = k + 2;
        %i = i+2;
    end
end
end