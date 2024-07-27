Neuron1=Neuron
% Neuron=Neuron1
% Type=input('Normal (1), MU (2), or SU (3)?')
% if Type==2
%     Neuron=MUNeuron
% elseif Type==3
%     Neuron=SUNeuron
% elseif Type==1
%     Neuron=Neuron
% end

for z=1:size(Neuron.Neuron,2)
%%% find stim driven combos
a=find(BLvsPS1(:,z)==1)
if isempty(a)
    b=nan
else b=mSDPS1(a,z)
end
%%% get the firing rates for these SDR
% b=mSDPS1(a,z)

%%%find the largest SDR
c=max(b)
%%% find the location of that largest one
d=find(mSDPS1(:,z)==c)
if length(d)>1
    d=d(1)
end
%%% find the SF of that combo
e=Combos(d,3)
%%% find the other trials with that SF
if isempty(e)
    f=nan
else
f=find(Combos(:,3)==e)
end
%%% see if those combos are SD
if isnan(f)==1
    mAISFP(z)=NaN
else
    for cc=1:length(f)
     g(cc)=BLvsPS1(f(cc),z)
    end
%%% get all relevant firing rates
for dd=1:length(g)
    if g(dd)==1
        h(dd)=mSDPS1(f(dd),z)
        i(dd)=mSDPS2(f(dd),z)
    else h(dd)=nan
        i(dd)=nan
    end
end
%%% get rid of nans
j=h(~isnan(h))
k=i(~isnan(i))
%%% get the mean firing rate
l=mean(j)  %% SD1
m=mean(k)  %% SD2

mAISFP(z,1)=(m-l)/l
mAISFP(z,2)=e
mAISFP(z,3)=find(Combos(1:5,3)==e)
end
end
for y=1:length(mAISFP)
    scatter(mAISFP(y,3),mAISFP(y,1),'MarkerEdgeColor','k','LineWidth',1.5)
    hold on
    xlim([0 6])
    ylim([-2 1])
    box off
    yline([0 0])
    title('Adulthood: P61+')
end
