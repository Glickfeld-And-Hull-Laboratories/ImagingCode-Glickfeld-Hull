figure
hold on
%tiledlayout('flow');
for n=1:length(CS)%length(AboutOneHz)%length([GoodUnitStruct.unitID])
    %if find(divisors(n) == 2)
    %    figure
    %    tiledlayout('flow');
    %end
    %nexttile
    unit = CS(n,1);
    unitIN = find([AllUnitStruct.unitID] == unit); %AboutOneHz(n,1)); % finds index in GoodUnitStruct that corresponds with first unit on list
    Hz = GeneralLineForStruct(JuiceTimesAdj, unitIN, AllUnitStruct, -5, 5, .1, 'k');
    hold on
    stunit = num2str(unit);                    % extract name/number of unit at index n (string that will be used in title of histogram)                         
    %stHz = num2str(Hz, '%.0f');
    %Depth = shallowToDeepGood(n,3);
    %stDepth = num2str(Depth);
    %title_ = strcat(stunit,', ', stDepth, ' deep, ', ',', stHz, ' Hz');   % make the title (can't get a space to show up after the comma)
    %title_ = ['Unit ', stunit, ', ', stDepth, ' deep, ', stHz, ' Hz'];   % make the title (can't get a space to show up after the comma)
    title_ = ['Unit ', stunit, ' Resp to Juice Delivery'];
    title(title_);
    %GeneralHistForStruct(FirstJuiceAdj, unit, GoodUnitStruct, -5, 5, .1, 'b');
    %GeneralHistForStruct(NoJuiceAdj, unit, GoodUnitStruct, -5, 5, .1, 'r');
end
