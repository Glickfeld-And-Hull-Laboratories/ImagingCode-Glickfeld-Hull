clear all
fullData=readtable('experimentList_012224.csv');

output=nan(size(fullData,2),2);

for i=1:size(fullData,1)
    [count, amp]=eventCountAndAmp(cell2mat(fullData.eventFileName(i)), ...
        fullData.startSweep(i), fullData.stopSweep(i),fullData.ConversionToPA(i));
   output(i,:)=[count amp];
end
clear i

fullData.eventPerMin=output(:,1);
fullData.eventAmp=output(:,2);

cells_WO_events=readtable('cells_WO_events.csv');

fullData = vertcat(fullData, cells_WO_events);

reducedData = fullData(:,[1,4,7,9,11,12]);
reducedData.drugCondition = categorical(reducedData.drugCondition);
reducedData.group=categorical(reducedData.group);
reducedData.eventHz = reducedData.eventPerMin/60;


drugOnly=reducedData(reducedData.drugCondition=='NBQX',:);
bslnOnly=reducedData(reducedData.drugCondition=='none',:);
%only has the cells that have both baseline ad drug conditions
matched=join(drugOnly,bslnOnly,'Keys','cellID');
matched.drugCondition_bslnOnly=[];
matched.group_bslnOnly=[];
%% plotting

conMeans = [mean(matched.eventHz_bslnOnly(matched.group_drugOnly=='con'),'omitmissing'),
    mean(matched.eventHz_drugOnly(matched.group_drugOnly=='con'),'omitmissing')];

expMeans = [mean(matched.eventHz_bslnOnly(matched.group_drugOnly=='exp'),'omitmissing'),
    mean(matched.eventHz_drugOnly(matched.group_drugOnly=='exp'),'omitmissing')];

conSTD = [std(matched.eventHz_bslnOnly(matched.group_drugOnly=='con'),[],'omitmissing'),
    std(matched.eventHz_drugOnly(matched.group_drugOnly=='con'),[],'omitmissing')];

nCon = length(matched.eventHz_bslnOnly(matched.group_drugOnly=='con'));
conSE=conSTD / nCon;

expSTD = [std(matched.eventHz_bslnOnly(matched.group_drugOnly=='exp'),[],'omitmissing'),
    std(matched.eventHz_drugOnly(matched.group_drugOnly=='exp'),[],'omitmissing')];

nExp = length(matched.eventHz_drugOnly(matched.group_drugOnly=='exp'));
expSE=expSTD / nExp;


figure
subplot(1,2,1)
plot([matched.eventHz_bslnOnly(matched.group_drugOnly=='con'),matched.eventHz_drugOnly(matched.group_drugOnly=='con')]')
hold on
scatter([1,2],[matched.eventHz_bslnOnly(matched.group_drugOnly=='con'),matched.eventHz_drugOnly(matched.group_drugOnly=='con')]',10, ...
    "MarkerEdgeColor","none","MarkerFaceColor",[0 0 0],"MarkerFaceAlpha",.25)
errorbar([1,2],conMeans, conSE,'.-k','MarkerSize',15)
xticks([1,2])
xticklabels({'Baseline','NBQX'})
set(gca,'TickDir','out')
box off
ylim([0 10])
xlim([0.5 2.5])

subplot(1,2,2)
plot([matched.eventHz_bslnOnly(matched.group_drugOnly=='exp'),matched.eventHz_drugOnly(matched.group_drugOnly=='exp')]',"Color",[0 0 1 .25])
hold on
scatter([1,2],[matched.eventHz_bslnOnly(matched.group_drugOnly=='exp'),matched.eventHz_drugOnly(matched.group_drugOnly=='exp')]',10, ...
    "MarkerEdgeColor","none","MarkerFaceColor",[0 0 1],"MarkerFaceAlpha",.25)
errorbar([1,2],expMeans, expSE,'.-b','MarkerSize',15)
xticks([1,2])
xticklabels({'Baseline','NBQX'})
set(gca,'TickDir','out')
box off
ylim([0 10])
xlim([0.5 2.5])
x0=5;
y0=5;
width=1.5;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(fullfile('ePSCS_Hz.pdf'),'-dpdf')
%% two-way ANOVA for epcs per minute
matchedIDs=categorical(unique(matched.cellID));

testTable=vertcat(bslnOnly,drugOnly);
testTable.cellID = categorical(testTable.cellID);
testTable=testTable(ismember(testTable.cellID,matchedIDs),:);

lme1 = fitlme(testTable, 'eventHz~drugCondition*group+(1|cellID)');
anova(lme1)

anova(testTable, 'eventHz ~ group * drugCondition');


[h,p1,~,stats1] =ttest2(matched.eventHz_bslnOnly(matched.group_drugOnly=='con'),matched.eventHz_bslnOnly(matched.group_drugOnly=='exp'));
[h,p2,~,stats2] =ttest2(matched.eventHz_drugOnly(matched.group_drugOnly=='con'),matched.eventHz_drugOnly(matched.group_drugOnly=='exp'));
[p1*2, p2*2]


[h,p3,~,stats3] =ttest(matched.eventHz_bslnOnly(matched.group_drugOnly=='con'),matched.eventHz_drugOnly(matched.group_drugOnly=='con'));
[h,p4,~,stats4] =ttest(matched.eventHz_bslnOnly(matched.group_drugOnly=='exp'),matched.eventHz_drugOnly(matched.group_drugOnly=='exp'));
[p3*2, p4*2]
%%
bslnMeans = [mean(bslnOnly.eventAmp(bslnOnly.group=='con'),'omitmissing'),
    mean(bslnOnly.eventAmp(bslnOnly.group=='exp'),'omitmissing')];

bslnSTD = [std(bslnOnly.eventAmp(bslnOnly.group=='con'),[],'omitmissing'),
    std(bslnOnly.eventAmp(bslnOnly.group=='exp'),[],'omitmissing')];

%nByGroup = [length(bslnOnly.eventAmp(bslnOnly.group=='con')), length(bslnOnly.eventAmp(bslnOnly.group=='exp'))]';
nByGroup = [9,14]';
bslnSE=bslnSTD ./ nByGroup;

figure;
swarmchart(bslnOnly.group,bslnOnly.eventAmp,10,'filled')
hold on
errorbar([1,2],bslnMeans, bslnSE,'.-k','MarkerSize',15)
xticklabels({'Control','DART'})
title('ePSCs amplitude(pA)')
set(gca,'TickDir','out')
box off
ylim([0 50])
x0=5;
y0=5;
width=1;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(fullfile('ePSC_amp.pdf'),'-dpdf')
%% t-test for amplitude
[h,p,ci,stats] = ttest2(bslnOnly.eventAmp(bslnOnly.group=='con'),bslnOnly.eventAmp(bslnOnly.group=='exp'))
%%
figure
boxchart(matched.controlType_bslnOnly,matched.eventHz_bslnOnly)
figure
boxchart(matched.controlType_drugOnly,matched.eventHz_drugOnly)