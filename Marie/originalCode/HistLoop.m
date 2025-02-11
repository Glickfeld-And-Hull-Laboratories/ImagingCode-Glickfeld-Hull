% simple code to generate all possible histograms using one stimulus and
% all the units in in a structure using GeneralHistFroStruct.

figure
tiledlayout('flow');
for n=1:length(shallowToDeepGood)%length([GoodUnitStruct.unitID])%length(AboutOneHz)%length([GoodUnitStruct.unitID])
    if find(divisors(n) == 50)
       figure
        tiledlayout('flow');
    end
    nexttile
    unit = shallowToDeepGood(n,1);
    unitIN = find([AllUnitStruct.unitID] == unit); %AboutOneHz(n,1)); % finds index in GoodUnitStruct that corresponds with first unit on list
    Hz = GeneralHistForStruct(LaserStimAdj, unitIN, AllUnitStruct, -.1, .1, .001, 'k');
    hold on
    stunit = num2str(unit);                    % extract name/number of unit at index n (string that will be used in title of histogram)                         
    stHz = num2str(Hz, '%.0f');
    Depth = shallowToDeepGood(n,4);
    stDepth = num2str(Depth);
    %title_ = strcat(stunit,', ', stDepth, ' deep, ', ',', stHz, ' Hz');   % make the title (can't get a space to show up after the comma)
    title_ = ['Unit ', stunit, ', ', stDepth, ' deep, ', stHz, ' Hz'];   % make the title (can't get a space to show up after the comma)
    title(title_);
    %GeneralHistForStruct(FirstJuiceAdj, unit, GoodUnitStruct, -5, 5, .1, 'b');
    %GeneralHistForStruct(NoJuiceAdj, unit, GoodUnitStruct, -5, 5, .1, 'r');
end


tiledlayout('flow');
for n=1:86%length([GoodUnitStruct.unitID])%length(AboutOneHz)%length([GoodUnitStruct.unitID])
    nexttile
    unit = shallowToDeepMua(n,1);
    unitIN = find([MultiUnitStruct.unitID] == unit); %AboutOneHz(n,1)); % finds index in GoodUnitStruct that corresponds with first unit on list
    Hz = GeneralHistForStruct(LaserStimAdj, unitIN, MultiUnitStruct, -.1, .3, .001, 'k');
    hold on
    stunit = num2str(unit);                    % extract name/number of unit at index n (string that will be used in title of histogram)                         
    stHz = num2str(Hz, '%.0f');
    Depth = shallowToDeepMua(n,3);
    stDepth = num2str(Depth);
    %title_ = strcat(stunit,', ', stDepth, ' deep, ', ',', stHz, ' Hz');   % make the title (can't get a space to show up after the comma)
    title_ = ['Unit ', stunit, ', ', stDepth, ' deep, ', stHz, ' Hz'];   % make the title (can't get a space to show up after the comma)
    title(title_);
    %GeneralHistForStruct(FirstJuiceAdj, unit, GoodUnitStruct, -5, 5, .1, 'b');
    %GeneralHistForStruct(NoJuiceAdj, unit, GoodUnitStruct, -5, 5, .1, 'r');
end


tiledlayout('flow');
for n=1:length([MultiUnitStruct.unitID])
    nexttile
    unit = find([MultiUnitStruct.unitID] == muaShallowToDeep(n,1));
    GeneralHistForStruct(LaserStimAdj, unit, MultiUnitStruct, -.1, .1, .005, 'k');
    %hold on
    %GeneralHistForStruct(LaserStim5, unit, MultiUnitStruct, -.1, .2, .005, 'r');
end

for i = 1:length(CSlist)
    for j = 1:length(SSlist)
    xCorrStructNewLimits(AllUnitStruct, -.02, .02, .001, CSlist(i), SSlist(j), 0, 300);
    end
end

for n=1:length([AllUnitStruct.unitID])
    FRstructRunNorun(ifrun, AllUnitStruct, AllUnitStruct(n).unitID);
end

figure
tiledlayout('flow');
for n=1:length([GoodUnitStruct.unitID])%length(AboutOneHz)%length([GoodUnitStruct.unitID])
    if find(divisors(n) == 50)
       figure
        tiledlayout('flow');
    end
    nexttile
    unitIN = n;
    unit = GoodUnitStruct(n,1);
    unit
    Hz = GeneralHistForStruct(LaserStimAdj, unitIN, AllUnitStruct, -.1, .3, .001, 'k');
    hold on
    stunit = num2str(unit);                    % extract name/number of unit at index n (string that will be used in title of histogram)                         
    stHz = num2str(Hz, '%.0f');
    %Depth = shallowToDeepMua(n,4);
    %stDepth = num2str(Depth);
    %title_ = strcat(stunit,', ', stDepth, ' deep, ', ',', stHz, ' Hz');   % make the title (can't get a space to show up after the comma)
    title_ = ['Unit ', stunit, ', ', stHz, ' Hz'];   % make the title (can't get a space to show up after the comma)
    title(title_);
    %GeneralHistForStruct(FirstJuiceAdj, unit, GoodUnitStruct, -5, 5, .1, 'b');
    %GeneralHistForStruct(NoJuiceAdj, unit, GoodUnitStruct, -5, 5, .1, 'r');
end
