%analysis script to get mean events/min seperated by drug condition for
%each DART type (can also seperate by ddHTP vs. +HPT)
close all
clear all

%%
expTable = readtable('expCellsList_231002.csv');

expCounts = zeros(size(expTable,1),2);
expAmps  = zeros(size(expTable,1),2);

for iCell = 1:size(expTable,1)
    fileName = expTable.filenames{iCell};
    drugStart = expTable.drugSweep(iCell);
    totalSweeps = expTable.totalSweep(iCell);

    [baselineCounts, drugCounts, baselineAmps, drugAmps] = eventCountAndAmpDRUG(fileName,drugStart,totalSweeps);
    
    expCounts(iCell,:)=[baselineCounts, drugCounts];
    expAmps(iCell,:)=[baselineAmps, drugAmps];

end


%manually adding two cell that had zero events, 23724001_cell2 and
%230615_cell4
expCounts(8:9,:)=[0,0;0,0];
expAmps(8:9,:)=[0,0;0,0];

%replace nans with 0
expCounts(isnan(expCounts))=0;
expAmps(isnan(expAmps))=0;
%%
cntrlTable = readtable('conCellsList_231002.csv');

cntrlCounts = zeros(size(cntrlTable,1),2);
cntrlAmps  = zeros(size(cntrlTable,1),2);

for iCell = 1:size(cntrlTable,1)
    fileName = cntrlTable.filenames{iCell};
    drugStart = cntrlTable.drugSweep(iCell);
    totalSweeps = cntrlTable.totalSweep(iCell);

    [baselineCounts, drugCounts, baselineAmps, drugAmps] = eventCountAndAmpDRUG(fileName,drugStart,totalSweeps);
    
    cntrlCounts(iCell,:)=[baselineCounts, drugCounts];
    cntrlAmps(iCell,:)=[baselineAmps, drugAmps];

end

cntrlCounts(isnan(cntrlCounts))=0;
cntrlAmps(isnan(cntrlAmps))=0;

%% plotting
figure;
scatter([1 2],cntrlCounts,'black')
hold on
scatter([1 2],expCounts,'blue')
plot(cntrlCounts','black')
plot(expCounts','blue')
plot(nanmean(cntrlCounts),'*k','MarkerSize',12)
plot(nanmean(expCounts),'*b','MarkerSize',12)
xticks([1,2])
xticklabels({'baseline','NBQX'})
ylabel("Spontaneous ePSCs / min")
set(gca,'TickDir','out')
box off
%ylim([-5 400])
xlim([0.5 2.5])

print('epscCountByDrug.pdf','-dpdf','-bestfit')


figure;
scatter([1 2],cntrlAmps,'black')
hold on
scatter([1 2],expAmps,'blue')
plot(cntrlAmps','black')
plot(expAmps','blue')
plot(nanmean(cntrlAmps),'*k','MarkerSize',12)
plot(nanmean(expAmps),'*b','MarkerSize',12)
xticks([1,2])
xticklabels({'baseline','NBQX'})
ylabel("ePSC amplitude (nA)")
set(gca,'TickDir','out')
box off
%ylim([-5 400])
xlim([0.5 2.5])
print('epscAmpByDrug.pdf','-dpdf','-bestfit')
