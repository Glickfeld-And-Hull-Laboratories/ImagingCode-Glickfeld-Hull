figure
hold on
%tiledlayout('flow');
AlllHz = 0;
TotalN = 0;
k = 1;
title_ = 'SS response to JuiceDelivery';
for n=1:length(SS)%length(AboutOneHz)%length([GoodUnitStruct.unitID])
    %if find(divisors(n) == 2)
    %    figure
    %    tiledlayout('flow');
    %end
    %nexttile
    unit = SS(n,1);
    unitIN = find([AllUnitStruct.unitID] == unit); %AboutOneHz(n,1)); % finds index in GoodUnitStruct that corresponds with first unit on list
    [Hz,N,edges,HzSd] = GeneralLineForStruct(JuiceTimesAdj, unitIN, AllUnitStruct, -5, 5, .1, 'k');
    AllHz(k,:) = Hz;
    AllN(k,:) = N;
    AllSd(k,:) = HzSd;
    AllZ(k,:) = (N - Hz)/HzSd;
    hold on
    stunit = num2str(unit);                    % extract name/number of unit at index n (string that will be used in title of histogram)                         
    %stHz = num2str(Hz, '%.0f');
    %Depth = shallowToDeepGood(n,3);
    %stDepth = num2str(Depth);
    %title_ = strcat(stunit,', ', stDepth, ' deep, ', ',', stHz, ' Hz');   % make the title (can't get a space to show up after the comma)
    %title_ = ['Unit ', stunit, ', ', stDepth, ' deep, ', stHz, ' Hz'];   % make the title (can't get a space to show up after the comma)
    %title_ = ['Unit ', stunit, ' Resp to Juice Delivery'];
    title(title_);
    %GeneralHistForStruct(FirstJuiceAdj, unit, GoodUnitStruct, -5, 5, .1, 'b');
    %GeneralHistForStruct(NoJuiceAdj, unit, GoodUnitStruct, -5, 5, .1, 'r');
    k = k+1;
end
AvgHz = mean(AllHz);
AvgN = mean(AllN);
AvgZ = mean(AllZ);
stDev = std(AllN, 0, 1);
sterr = (stDev/sqrt(k));
ci95 = 1.96 * (stDev/sqrt(k));
edges = edges(1:(length(edges)-1)); % remove last traiiling edge so sizes of N and edges match)
%errorbar(edges, AvgN, sterr, 'r');
hold off
figure
hold on
plot(edges, AllZ, 'color', [.85 .85 .85]);
plot(edges, AvgZ, 'r');
xline(0,'b');
title(title_)
box off
%ax.TickDir = 'out'
ax = gca; 
ax.TickDir = 'out';
ax.FontName = 'Times New Roman';
ax.FontSize = 10;

