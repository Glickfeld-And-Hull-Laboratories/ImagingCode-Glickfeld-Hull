function  unitNoBurst= RemoveBurst(unit, time)
k=1;
j = 2;
unitNoBurst(1) = unit(1);
for i = 2:length(unit)
    if (((unit(i) - unit(i-1)) <= time)) %&& ((unit(i+2) - unit(i+1)) < time)
        burstSpikes(k,:)= unit(i-1);
        burstSpikes(k+1, :) = unit(i);
        k = k + 2;
        %i = i+2;
    else
        unitNoBurst(j) = unit(i);
        j = j+1;
    end
end
end