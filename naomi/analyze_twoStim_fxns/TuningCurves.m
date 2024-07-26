function [vis] = TuningCurves(Combos,Ori,BLvsPS1,mSDPS1,k)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
for i=1:length(Combos)
    for j=1:length(Ori)
        if Combos(i,2)==Ori(j)
            TempMat(i,j)=BLvsPS1(i,k)
            TempTC(i,j)=mSDPS1(i,k)
        else TempMat(i,j)=NaN
            TempTC(i,j)=NaN
        end
    end
end
for i=1:length(Ori)
    TempCk=TempMat(:,i)
    TempCk=TempCk(~isnan(TempCk))
    TempTC2=TempTC(:,i)
    TempTC2=TempTC2(~isnan(TempTC2))
    if sum(TempCk)>0
        figure
        vis=scatter(1:length(TempTC2),TempTC2,'MarkerEdgeColor','k')
        hold on
        for m=1:length(TempCk)
            if TempCk(m)==1
                scatter(m,TempTC2(m),'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor','k')
                title(['Tuning Curve; Orientation ',num2str(Ori(i)),' degrees'])
                subtitle(['Neuron',num2str(k)])
            end
        end
        end
    end
end