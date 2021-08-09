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
    Hz = GeneralHistForStruct(LaserStimAdj, unitIN, AllUnitStruct, -.02, .02, .001, 'k');
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