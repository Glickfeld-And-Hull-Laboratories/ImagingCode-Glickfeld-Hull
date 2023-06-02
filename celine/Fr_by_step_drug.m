function [counts] = Fr_by_step_drug(dataIN,nSweeps)

% a script to get firing rate vs. current injection, averaged over sweeps
% within a condition and plotted seperately for sweeps of different [drug]
% conditions. Inputs: data file name, sweep when drug condition starts,
% number of sweeps, pA at each sweep. Assumes that the steps happen at
% specific times.

drugOn=12;
steps = [10 15 20 25];

inputFile = [dataIN,'.csv']
fullData =readtable(inputFile);
data=fullData(:,[1 5]); %pare down to just sweep and time
data=table2array(data);

stepInds = cell(1,4);
stepInds{1}=find(data(:,2)>11000 & data(:,2)<11500);
stepInds{2}=find(data(:,2)>16500 & data(:,2)<17000);
stepInds{3}=find(data(:,2)>22000 & data(:,2)<22500);
stepInds{4}=find(data(:,2)>27500 & data(:,2)<28000);

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
Fr = counts.*2;
%split into baseline and drug datasets
bsln = Fr(1:(drugOn-1),:);
drug = Fr(drugOn:nSweeps,:);
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