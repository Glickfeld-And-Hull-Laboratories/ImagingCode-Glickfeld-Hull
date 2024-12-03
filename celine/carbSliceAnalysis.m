wd = 'Z:\home\celine\Analysis\ePhys\VIP_cells_carbachol';
cd(wd);

cellList = [5,6,7];
currentResponse = nan(20,length(cellList));
for iCell=1:length(cellList)
    fileName = ['Cell',num2str(cellList(iCell)),'.csv'];
    tempData=readtable( fileName);
    currentResponse(:,iCell)=table2array(tempData(1:20,5));
    clear tempData
end
figure
plot(currentResponse,'b');
box off
set(gca,'TickDir','out')
ylabel('delta pA')
xlabel('Sweep')
hold on

cellList = [4];
currentResponse = nan(20,length(cellList));
for iCell=1:length(cellList)
    fileName = ['Cell',num2str(cellList(iCell)),'.csv'];
    tempData=readtable( fileName);
    currentResponse(:,iCell)=table2array(tempData(1:20,5));
    clear tempData
end

plot(currentResponse,'--b');
box off
set(gca,'TickDir','out')
ylabel('delta pA')
xlabel('Sweep')
ylim([-60 10])