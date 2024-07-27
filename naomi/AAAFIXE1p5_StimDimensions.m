

clearvars -except inputs

SF=cell2mat(inputs.tStimOneGratingSpatialFreqCPD);
Dir=inputs.tStimOneGratingDirectionDeg;
Dir=cell2mat(Dir);
Dir=double(Dir);
Trials=[SF' Dir'];
Trials(:,3)=Trials(:,1)+Trials(:,2);

c=unique(SF)
d=unique(Dir)
%%
k=[0:length(c):(length(c)*length(d))]
for j=1:length(d)
for i=1:length(c)
    Combos1(i,1)=c(i)
    Combos2(i,1)=d(j)
end
TempCombos=[Combos2+Combos1]
Combos(k(j)+1:length(TempCombos)*j,:)=TempCombos
end
Combos(:,2)=1:length(Combos)
clearvars -except inputs Trials Combos

%%
for i=1:length(Trials)
    Trials(i,4)=find(Combos==Trials(i,3))
end

