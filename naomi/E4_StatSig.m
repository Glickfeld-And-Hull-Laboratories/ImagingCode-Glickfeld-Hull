-+Neuron1=Neuron
Neuron=Neuron1
Type=input('Normal (1), MU (2), or SU (3)?')
if Type==2
    Neuron=MUNeuron
elseif Type==3
    Neuron=SUNeuron
elseif Type==1
    Neuron=Neuron
end

for h=1:size(Neuron.Neuron,2)
for i=1:size(Neuron.Neuron(1).Stim,2)

mPrestim(i,h)=mean(struct2array(Neuron.Neuron(h).Stim(i).BLHz));
Prestim=struct2array(Neuron.Neuron(h).Stim(i).BLHz);

mPeristim1(i,h)=mean(struct2array(Neuron.Neuron(h).Stim(i).SDHz));
Peristim1=struct2array(Neuron.Neuron(h).Stim(i).SDHz);

mPeristim2(i,h)=mean(struct2array(Neuron.Neuron(h).Stim(i).SSHz));
Peristim2=struct2array(Neuron.Neuron(h).Stim(i).SSHz);

BLvsPS1(i,h)=ttest(Prestim',Peristim1','tail','left','Alpha',(0.5/size(Neuron.Neuron(1).Stim,2)));
PPS1vsPS2(i,h)=ttest(Peristim1',Peristim2','tail','left','Alpha',0.5/size(Neuron.Neuron(1).Stim,2));
APS1vsPS2(i,h)=ttest(Peristim1',Peristim2','tail','right','Alpha',0.5/size(Neuron.Neuron(1).Stim,2));

AI=(mPeristim2-mPeristim1)./mPeristim1;

end
end
%% Look at just the stimulus driven part (take out the baseline)
mSDPS1=mPeristim1-mPrestim
mSDPS2=mPeristim2-mPrestim
%%
clearvars -except mPrestim mSDPS1 mSDPS2 BLvsPS1 PPS1vsPS2 APS1vsPS2 AI Combos Trials Neuron MUNeuron SUNeuron
%% Find only the stimulus driven neurons and look at them in particular
%% Tuning Curves
Ori=unique(Combos(:,2));
SF=unique(Combos(:,3));
%%
TC=input('Would you like tuning curves? (yes=1, no=2)')
if TC==1
    addpath 'Z:\home\Naomi\Electrophysiology\New Analysis\analyze_twoStim_fxns\'
    k=input('Which Neuron? (You have ',num2str(sizemSDPS1,2),' neurons)')
[vis] = TuningCurves(Combos,Ori,BLvsPS1,mSDPS1,k)
end
%% first let's get the preferred stimulus for each neuron
for k=1:size(BLvsPS1,2)
    for i=1:size(BLvsPS1,1)
        if BLvsPS1(i,k)==1
            TempAI(i,k)=AI(i,k);
        else TempAI(i,k)=NaN;
        end
    end
end
%%
IDs=unique(Trials(:,4));
IDs=table2array(IDs);
for m=1:size(TempAI,2);
for i=length(TempAI)
    TempTempAI(:,1)=TempAI(:,m);
    TempTempAI(:,2)=APS1vsPS2(:,m);
    TempTempAI(:,3)=PPS1vsPS2(:,m);
    TempTempAI(1:size(TempAI,1),4)=zeros;
    TempPref=max(mSDPS1(:,m));
    PrefLoc=find(mSDPS1(:,m)==TempPref);
    TempTempAI(PrefLoc,4)=1;
end

    scatter(1:length(IDs),TempTempAI(:,1),'MarkerEdgeColor','k')
    xlim([0 length(IDs)])
    hold on
for k=1:length(TempTempAI)
    if TempTempAI(k,2)==1
        scatter(k,TempTempAI(k,1),'MarkerFaceColor','k','MarkerEdgeColor','k')
    end
    hold on
    if TempTempAI(k,3)==1
        scatter(k,TempTempAI(k,1),'MarkerFaceColor','k','MarkerEdgeColor','k')
    end
    hold on
    if TempTempAI(k,4)==1
        scatter(k,TempTempAI(k,1),'MarkerEdgeColor','b')
    end
    
end
title('Early Postnatal: P0-23')
% xline(length(SF),'--')
% xline(length(SF)*2,'--')
yline([0 0])
xlabel('Spatial Frequency')
ylabel('Adaptation Index')
end