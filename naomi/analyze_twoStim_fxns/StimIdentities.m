function [StimTrialTimes,Combos,Trials,StimIdent] = StimIdentities(dimensions,inputs,trigger_timesSD,trigger_timesBL)
% dimensions=input('SF(1) or SF and DIR(2)?')
if dimensions==1
TI=cell2mat(inputs.tStimTwoGratingSpatialFreqCPD);
TI=TI(1:length(trigger_timesBL));
StimIdent=unique(TI);
elseif dimensions==2
SF=cell2mat(inputs.tStimOneGratingSpatialFreqCPD);
Dir=inputs.tStimOneGratingDirectionDeg;
Dir=cell2mat(Dir);
Dir=double(Dir);
Trials=[SF' Dir'];
Trials(:,3)=Trials(:,1)+Trials(:,2);
%%

c=unique(SF);
d=unique(Dir);
k=[0:length(c):(length(c)*length(d))];
for j=1:length(d)
for i=1:length(c)
    % if length(c==1)
    %     Combos2=d
    % else
    Combos1(i,1)=c(i);
    Combos2(i,1)=d(j);
    % end
end
TempCombos=[Combos2+Combos1 Combos2 Combos1];
% TempCombos=[Combos2+Combos1 ]

Combos(k(j)+1:length(TempCombos)*j,:)=TempCombos;

end
%%

Combos(:,4)=1:length(Combos);
for i=1:length(Trials)
    Trials(i,4)=find(Combos(:,1)==Trials(i,3));
end
TI=Trials(:,4)';
StimIdent=unique(TI);
%%

end
StimIdent=unique(TI);
for k=1:length(StimIdent)
    for h=1:length(trigger_timesSD);
    if TI(h)==StimIdent(k)
        TrialStimID(h,k)=trigger_timesSD(h);
    else TrialStimID(h,k)= NaN;
    end
    end
end



for i=1:length(StimIdent)
    A=rmmissing(TrialStimID(:,i));
    StimTrialTimes(i).ID=A;
end

end