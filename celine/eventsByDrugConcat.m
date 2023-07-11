%analysis script to get mean events/min seperated by drug condition for
%each DART type (can also seperate by ddHTP vs. +HPT)

%enter a list of cells for each DART condition
expList={'230623_cell2' '230623_cell3' '230623_cell4' '230615_cell3' '230615_cell6'}; %last two are NES

cntrlList={'230621_cell2' '230621_cell3'  '230623_cell1' '230615_cell2' '230616_cell1' '230616_cell2' '230616_cell3' '230616_cell4'}; %last 5 are NES

exp = zeros(size(expList,2),2);
for iCell = 1:3
    [exp(iCell,1),exp(iCell,2)]=eventsPerMinDRUG(expList{iCell},50,10,20,30,20);
end

%here, I know that the last two cells have 0 epsc's and so the files won't
%ready properly - to take care of this I made the exp dataframe all zeros
%to start and I leave those two rows unchanged.

cntrl = nan(size(cntrlList,2),2);

for iCell = 1:3
    [cntrl(iCell,1),cntrl(iCell,2)]=eventsPerMinDRUG(cntrlList{iCell},50,10,20,30,20);
end
%cells 5-8 have different parameters 
for iCell = 4:8
    [cntrl(iCell,1),cntrl(iCell,2)]=eventsPerMinDRUG(cntrlList{iCell},70,10,20,50,20);
end

%% looking at access resistance

expQC = nan(size(expList,2),2);
cntrlQC = nan(size(cntrlList,2),2);

for iCell = 1:size(expList,2)
    inputFile = [expList{iCell},'_statistics.csv']
    fullData =readtable(inputFile);
    data=fullData(:,[5;12]);
    data=table2array(data);
    expQC(iCell,1) = mean(data(:,2));
    expQC(iCell,2) = std(data(:,1));
    clear fullData data inputFile
end

for iCell = 1:size(cntrlList,2)
    inputFile = [cntrlList{iCell},'_statistics.csv']
    fullData =readtable(inputFile);
    data=fullData(:,[5;12]);
    data=table2array(data);
    cntrlQC(iCell,1) = mean(data(:,2));
    cntrlQC(iCell,2) = std(data(:,1));
    clear fullData data inputFile
end
%% removing cells that don't pass R_access QC 
% right now only requirement is series resistance less than 30MOhm

exp_pass = find(expQC < 30);
cntrl_pass = find(cntrlQC < 30);

%% plotting
figure;
scatter([1 2],cntrl(cntrl_pass,:),'black')
hold on
scatter([1 2],exp(exp_pass,:),'blue')
plot(cntrl(cntrl_pass,:)','black')
plot(exp(exp_pass,:)','blue')
plot(mean(cntrl(cntrl_pass,:)),'*k','MarkerSize',12)
plot(mean(exp(exp_pass,:)),'*b','MarkerSize',12)
xticks([1,2])
xticklabels({'baseline','NBQX'})
ylabel("Spontaneous ePSCs / min")
set(gca,'TickDir','out')
box off
ylim([-5 400])
xlim([0.5 2.5])
