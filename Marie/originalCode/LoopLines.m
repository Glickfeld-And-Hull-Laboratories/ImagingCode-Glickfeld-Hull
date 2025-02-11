figure
tiledlayout('flow');
for n=1:10%length(chanlist)%length(AboutOneHz)%length([GoodUnitStruct.unitID])
    %if find(divisors(n) == 50)
     %  figure
      %  tiledlayout('flow');
    %end
    nexttile
   [time, Waveforms, SampleTS] = SampleWaveformsTS( .15, 100, LaserStim50ms, .05, chanlist(n)-5);
    hold on
    stchan = num2str(chanlist(n));                    % extract name/number of unit at index n (string that will be used in title of histogram)                         
    %Depth = shallowToDeepMua(n,4);
    %stDepth = num2str(Depth);
    %title_ = strcat(stunit,', ', stDepth, ' deep, ', ',', stHz, ' Hz');   % make the title (can't get a space to show up after the comma)
    title_ = ['chan ', stchan ];   % make the title (can't get a space to show up after the comma)
    title(title_);
    %GeneralHistForStruct(FirstJuiceAdj, unit, GoodUnitStruct, -5, 5, .1, 'b');
    %GeneralHistForStruct(NoJuiceAdj, unit, GoodUnitStruct, -5, 5, .1, 'r');
end
