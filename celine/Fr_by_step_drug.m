function [counts,bsln,drug] = Fr_by_step_drug(dataIN,drugOn)

% a script to get firing rate vs. current injection, averaged over sweeps
% within a condition and plotted seperately for sweeps of different [drug]
% conditions. Inputs: data file name, sweep when drug condition starts,
% number of sweeps, pA at each sweep. Assumes that the steps happen at
% specific times.

steps = [10 15 20 25];


fullData =readtable(dataIN);
data=fullData(:,[1 5]); %pare down to just sweep and time
data=table2array(data);
nSweeps=max(data(:,1)); %find how many sweeps there are in this full dataset
stepInds = cell(1,4);
stepInds{1}=find(data(:,2)>5500 & data(:,2)<5750);
stepInds{2}=find(data(:,2)>8250 & data(:,2)<8500);
stepInds{3}=find(data(:,2)>11000 & data(:,2)<11250);
stepInds{4}=find(data(:,2)>13750 & data(:,2)<14000);

counts=zeros(nSweeps,4);

for iSweep = 1:nSweeps
    for iStep = 1:4
        if ~isempty(find(data(:,1)==iSweep))
            sweepInds = find(data==iSweep);
            tempStepInds = stepInds{iStep};
            counts(iSweep,iStep)=length(intersect(sweepInds,tempStepInds));
        end
    end
end
Fr = counts.*25;
%split into baseline and drug datasets

bsln = Fr((drugOn-10):(drugOn),:);

drug = Fr((drugOn+20):(drugOn+30),:);


figure;
plot(steps,bsln,'color', [.5 .5 .5], 'linewidth', 1)
hold on
plot(steps,mean(bsln,1),'k','linewidth', 1.5)
hold on
plot(steps,drug,'color', [ 0.5843    0.8157    0.9882], 'linewidth', 1)
hold on
plot(steps,mean(drug,1),'b','linewidth', 1.5)
hold off
title('Firing rate vs. step')
box off
xlabel("mV step")
ylabel("Firing rate")
set(gca,'TickDir','out')

%% get Vm and IR measurements 




