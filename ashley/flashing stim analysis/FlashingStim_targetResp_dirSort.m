edit FlashingStim_dataSortedByCycle_combineDatasets.m
%%
for icyc = 1:length(cycles)
    ind = find(tCyclesOn == cycles(icyc));   
    cycTrialOutcome{icyc} = trialOutcome(ind);
    cycDirectionDeg{icyc} = DirectionDeg(ind);
    cycCatchDirectionDeg{icyc} = catchDirectionDeg(ind);
    cycCatchTrialOutcome{icyc} = catchTrialOutcome(ind);
    cycCatchCycle{icyc} = catchCycle(ind);
end

% Success_ind = find(strcmp('success',trialOutcome));

% DirectionDeg_success_ind = DirectionDeg(Success_ind);
% 
% load('cellSelectivityIndices.mat')
% cellsPrefZero = find(dirPref_ind == 1 | dirPref_ind == 5);
% cellsPrefNinety = find(dirPref_ind == 3 | dirPref_ind == 7);
% cellsSlctvZero = intersect(cellsPrefZero,dirSlctvCells);

%% color-code cells by orientation preference
DirFolder = '006';
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, DirFolder);
cd(fileSave);
% load('oriTuningPreferences.mat')
load('TuningPreferences.mat')

dataTrialStart = cycDataDFoverF_cmlvNoTarget{4};
v_ind = cycV_ind{4};
% a_ind = cycA_ind{1};
a_ind = cycAV_ind{4};


preStimResp_V = zeros(size(v_ind,2),size(dataTrialStart,2));
for itrial =1:size(v_ind,1);
    for icell = 1:size(dataTrialStart,2)
        preStimResp_V(itrial,icell) = mean(dataTrialStart(1:30,icell,v_ind(itrial)),1);
    end
end

baselineStimResp_V = zeros(size(v_ind,2),size(dataTrialStart,2));
for itrial = 1:size(v_ind,1);
    for icell = 1:size(dataTrialStart,2)
        baselineStimResp_V(itrial,icell) = mean(dataTrialStart(36:end,icell,v_ind(itrial)),1);
    end
end

baselineStimRespTtest_V= ttest(preStimResp_V,baselineStimResp_V,'alpha', 0.01);
baselineStimRespIndex_V = find(baselineStimRespTtest_V == 1);


cellsPrefZero = find(dirPref_ind == 1 | dirPref_ind == 5);
cellsSelectZero_dir = intersect(dirSlctvCells,cellsPrefZero);
cellsSelectZero_ori = intersect(oriSlctvCells,cellsPrefZero);
cellsSelectZero = union(cellsSelectZero_dir,cellsSelectZero_ori);
cellsPrefRespZero = intersect(baselineStimRespIndex_V,cellsPrefZero);

cellsPrefNinety = find(dirPref_ind == 3 | dirPref_ind == 7);
cellsSelectNinety_dir = intersect(dirSlctvCells,cellsPrefNinety);
cellsSelectNinety_ori = intersect(oriSlctvCells,cellsPrefNinety);
cellsSelectNinety = union(cellsSelectNinety_dir,cellsSelectNinety_ori);
cellsPrefRespNinety = intersect(baselineStimRespIndex_V,cellsPrefNinety);

nCells = size(cycDataDFoverF_cmlvNoTarget{7},2);
oriSlctvCellsAll = union(oriSlctvCells,dirSlctvCells);
notSlctvCells = setdiff([1:nCells],oriSlctvCellsAll);
notRespCells = setdiff([1:nCells],baselineStimRespIndex_V);

%% sort cell trials by target grating orientation (average each) (successes only)
% for icyc = 1:length(cycles)
%     tempdata = cycDataDFoverF{icyc};
%     cycdirs = cycDirectionDeg{icyc};
%     cycoutcome = cycTrialOutcome{icyc};
%     cycsuccess = find(strcmp(cycoutcome,'success'));
%     cycdirmean = zeros(size(tempdata,1),size(tempdata,2),length(Dirs));
%     for i = 1: length(Dirs)
%         ind = intersect(cycsuccess, find(cycdirs == Dirs(i)));
%         if isempty(ind)
%             cycdirmean = [];
%         else
%             cycdirmean(:,:,i) = mean(tempdata(:,:,ind),3);
%         end
%     end
%     cycTargetDirMean{icyc} = cycdirmean;
% end

%% plot target response, ea cycle
% figure;
% for icyc = 1:length(cycles)
%     tempdata = cycTargetDirMean{icyc};
%     tempdatamean = squeeze(mean(tempdata,2));
%     
%     subplot(2,5,icyc)
%     plot(tempdatamean(:,end))
%     hold on
%     vline(30,'k');
%     hold on
%     for i = 1:cycles(icyc)-1
%         L = (i*cycTime)+30;
%         vline(L,'k:');
%         hold on
%     end
%     vline((cycles(icyc)*cycTime+30),'c');
%     hold on
% end
%     
% hold on
% for i = 1: length(Dirs)
%     legendInfo{i} = [num2str(Dirs(i)) ' degrees'];
% end
% legend(legendInfo,'Location','SouthEast')

%% find each target resp

for i = 1:length(Dirs)
    ind = intersect(find(DirectionDeg == Dirs(i)),find(strcmp(trialOutcome,'success')));
    targetData = zeros(60,size(dataTC,2),length(ind));
    targetDataDF = zeros(60,size(dataTC,2),length(ind));
    targetDataDFoverF = zeros(60,size(dataTC,2),length(ind));
    targetBL = zeros(length(Dirs),size(dataTC,2),length(ind));
    targetNorm = zeros(60,size(dataTC,2),length(ind));
    
    if cTargetOn(ind(end),1)+40 > size(dataTC,1) 
        for itrial = 1:length(ind)-1
            targetData(:,:,itrial) = dataTC(cTargetOn(ind(itrial))-19:cTargetOn(ind(itrial))+40,:);
            targetBaseline = dataTC(cLeverDown(ind(itrial))-29:cLeverDown(ind(itrial)),:);
            targetDataDF(:,:,itrial) = bsxfun(@minus, targetData(:,:,itrial), mean(targetBaseline,1));
            targetDataDFoverF(:,:,itrial) = bsxfun(@rdivide, targetDataDF(:,:,itrial), mean(targetBaseline,1));
            targetBL(i,:,itrial) = mean(targetDataDFoverF(17:21,:,itrial),1);
            targetNorm(:,:,itrial) = bsxfun(@minus,targetDataDFoverF(:,:,itrial),targetBL(i,:,itrial));
        end
    else
        for itrial = 1:length(ind)
            targetData(:,:,itrial) = dataTC(cTargetOn(ind(itrial))-19:cTargetOn(ind(itrial))+40,:);
            targetBaseline = dataTC(cLeverDown(ind(itrial))-29:cLeverDown(ind(itrial)),:);
            targetDataDF(:,:,itrial) = bsxfun(@minus, targetData(:,:,itrial), mean(targetBaseline,1));
            targetDataDFoverF(:,:,itrial) = bsxfun(@rdivide, targetDataDF(:,:,itrial), mean(targetBaseline,1));
            targetBL(i,:,itrial) = mean(targetDataDFoverF(17:21,:,itrial),1);
            targetNorm(:,:,itrial) = bsxfun(@minus,targetDataDFoverF(:,:,itrial),targetBL(i,:,itrial));
        end
    end
    
    dirTargetData{i} = targetData;
    dirTargetDataDF{i} = targetDataDF;
    dirTargetDataDFoverF{i} = targetDataDFoverF;
    dirTargetInd{i} = ind;
    dirTargetNorm{i} = targetNorm;
end

%% Normalize target response
% for icyc = 1:length(cycles)
%     tempdata = cycTargetDirMean{icyc};
%     targetindex = cycles(icyc)*cycTime+30;
%     if isempty(tempdata) == 1
%         targetBaselineF(icyc,:,:) = NaN(1,nCells,length(Dirs));
%     else
%         targetBaselineF(icyc,:,:) = mean(tempdata(targetindex-3:targetindex+2,:,:),1);
%     end
% end
% 
% for icyc = 1:length(cycles)
%     tempdata = cycTargetDirMean{icyc};
%     targetindex = cycles(icyc)*cycTime+30;
%     if isempty(tempdata) == 1
%         dirtraces(:,:,:,icyc) = NaN(40,nCells,length(Dirs));
%     else
%         tempdata = tempdata(targetindex-10:targetindex+29,:,:);
%         tempdatasub = bsxfun(@minus,tempdata,targetBaselineF(icyc,:,:));
%         for idir = 1:length(Dirs)
%             dirtraces(:,:,idir,icyc) = tempdatasub(:,:,idir);
%         end
%     end
% end

%%
figure;
% colors = brewermap(length(Dirs),'*YlGn');
colors = brewermap(length(Dirs)+15,'*YlGn');
colorind = [3:2:length(Dirs)+12];
colors = colors(colorind,:);

for i = 1:length(Dirs)
    tempdata = dirTargetNorm{i};
    plot(nanmean(nanmean(tempdata,3),2),'Color',colors(i,:))
    hold on
end
hold on
for i = 1: length(Dirs)
    legendInfo{i} = [num2str(Dirs(i)) ' degrees'];
end
legend(legendInfo,'Location','NorthWest')
hold on
vline(20,'k')
vline(20+tooFastTime,'k:')
vline(20+maxReactTime,'k:')
title([mouse ' - ' date '; target resp'])
% set(gca,'Color','k')

%%
cellGroup = oriSlctvCellsAll;
% cellsSubgroup1 = find(ismember(cellGroup,intersect(cellGroup,cellsPrefZero)));
% cellsSubgroup2 = find(ismember(cellGroup,intersect(cellGroup,cellsPrefNinety)));
% cellsSubgroup3 = find(ismember(cellGroup,intersect(cellGroup,setdiff([1:nCells],cat(1,cellsPrefZero,cellsPrefNinety)))));
cellsSubgroup1 = cellsSelectZero;
cellsSubgroup2 = setdiff(cellGroup,cat(1,cellsPrefZero,cellsPrefNinety));
cellsSubgroup3 = cellsSelectNinety;

% mean response to target
cellTargetDirMean = zeros(60,nCells,length(Dirs));
cellTargetDirMean_point = zeros(nCells,length(Dirs));
cellTargetdirMean_diff = zeros(nCells,length(Dirs));
for i = 1:length(Dirs)
    cellTargetDirMean(:,:,i) = nanmean(dirTargetNorm{i},3);
    cellTargetDirMean_point(:,i) = squeeze(nanmean(cellTargetDirMean(21:25,:,i),1));
    cellTargetDirMean_diff(:,i) = squeeze(cellTargetDirMean(21,:,i)-cellTargetDirMean(25,:,i));
end

% cellTargetDirMean = nanmean(dirtraces,4);
% cellTargetDirMean_point = squeeze(nanmean(cellTargetDirMean(21:25,:,:),1));

%differential of response to target

%% plot cell response to target by tuning preference (mean resp to target)

celltargetrespfig = figure;
% colors = brewermap(length(Dirs),'*YlGn');
colors = brewermap(length(Dirs)+15,'*YlGn');
colorind = [3:2:length(Dirs)+12];
colors = colors(colorind,:);

for idir = 1:length(Dirs)
    if idir == 1
        plot([1:length(cellsSubgroup1)],cellTargetDirMean_point(cellsSubgroup1,idir),'ko','MarkerFaceColor', 'k')
    else
        plot([1:length(cellsSubgroup1)],cellTargetDirMean_point(cellsSubgroup1,idir),'Color',colors(idir,:),'LineStyle','none','Marker','o','MarkerFaceColor', colors(idir,:))
    end
    hold on
end

hold on
vline(length(cellsSubgroup1),'k')

for idir = 1:length(Dirs)
    if idir == 1
        plot([length(cellsSubgroup1)+1:length(cellsSubgroup1)+length(cellsSubgroup2)],cellTargetDirMean_point(cellsSubgroup2,idir),'ko','MarkerFaceColor', 'k')
    else
        plot([length(cellsSubgroup1)+1:length(cellsSubgroup1)+length(cellsSubgroup2)],cellTargetDirMean_point(cellsSubgroup2,idir),'Color',colors(idir,:),'LineStyle','none','Marker','o','MarkerFaceColor', colors(idir,:))
    end
    hold on
end

hold on
vline(length(cellsSubgroup1)+length(cellsSubgroup2),'k')

for idir = 1:length(Dirs)
    if idir == 1
        plot([length(cellsSubgroup1)+length(cellsSubgroup2)+1:length(cellsSubgroup1)+length(cellsSubgroup2)+length(cellsSubgroup3)],cellTargetDirMean_point(cellsSubgroup3,idir),'ko','MarkerFaceColor', 'k')
    else
        plot([length(cellsSubgroup1)+length(cellsSubgroup2)+1:length(cellsSubgroup1)+length(cellsSubgroup2)+length(cellsSubgroup3)],cellTargetDirMean_point(cellsSubgroup3,idir),'Color',colors(idir,:),'LineStyle','none','Marker','o','MarkerFaceColor', colors(idir,:))
    end
    hold on
end

legend(legendInfo)
title([mouse ' - ' date ';mean resp to target'])

%% plot cell response to target by peak resp 

celltargetrespfig = figure;
% colors = brewermap(length(Dirs),'*YlGn');
colors = brewermap(length(Dirs)+15,'*YlGn');
colorind = [3:2:length(Dirs)+12];
colors = colors(colorind,:);

for idir = 1:length(Dirs)
    if idir == 1
        plot([1:length(cellsSubgroup1)],cellTargetDirMean_diff(cellsSubgroup1,idir),'ko','MarkerFaceColor', 'k')
    else
        plot([1:length(cellsSubgroup1)],cellTargetDirMean_diff(cellsSubgroup1,idir),'Color',colors(idir,:),'LineStyle','none','Marker','o','MarkerFaceColor', colors(idir,:))
    end
    hold on
end

hold on
vline(length(cellsSubgroup1),'k')

for idir = 1:length(Dirs)
    if idir == 1
        plot([length(cellsSubgroup1)+1:length(cellsSubgroup1)+length(cellsSubgroup2)],cellTargetDirMean_diff(cellsSubgroup2,idir),'ko','MarkerFaceColor', 'k')
    else
        plot([length(cellsSubgroup1)+1:length(cellsSubgroup1)+length(cellsSubgroup2)],cellTargetDirMean_diff(cellsSubgroup2,idir),'Color',colors(idir,:),'LineStyle','none','Marker','o','MarkerFaceColor', colors(idir,:))
    end
    hold on
end

hold on
vline(length(cellsSubgroup1)+length(cellsSubgroup2),'k')

for idir = 1:length(Dirs)
    if idir == 1
        plot([length(cellsSubgroup1)+length(cellsSubgroup2)+1:length(cellsSubgroup1)+length(cellsSubgroup2)+length(cellsSubgroup3)],cellTargetDirMean_diff(cellsSubgroup3,idir),'ko','MarkerFaceColor', 'k')
    else
        plot([length(cellsSubgroup1)+length(cellsSubgroup2)+1:length(cellsSubgroup1)+length(cellsSubgroup2)+length(cellsSubgroup3)],cellTargetDirMean_diff(cellsSubgroup3,idir),'Color',colors(idir,:),'LineStyle','none','Marker','o','MarkerFaceColor', colors(idir,:))
    end
    hold on
end

legend(legendInfo)
title([mouse ' - ' date ';target differential'])



%%
% L = 60;
% data_aroundTarget = zeros(L,size(dataTC,2),length(Success_ind));
% dF_aroundTarget = zeros(L,size(dataTC,2),length(Success_ind)); 
% dFoverF_aroundTarget = zeros(L,size(dataTC,2),length(Success_ind));
% for i = 1:length(Success_ind)
%     trial = Success_ind(i);
%     data_aroundTarget(:,:,i) = dataTC(cTargetOn(trial)-(L/2):cTargetOn(trial)+((L/2)-1),:);
%     dF_aroundTarget(:,:,i) = bsxfun(@minus, data_aroundTarget(:,:,i), mean(data_aroundTarget((L/2)-5:(L/2),:,i),1));
%     dFoverF_aroundTarget(:,:,i) = bsxfun(@rdivide, dF_aroundTarget(:,:,i),mean(data_aroundTarget((L/2)-5:(L/2),:,i),1));
% end
% 
% 
% V_Success_ind = find(ismember(Success_ind,intersect(Success_ind,V_ind)));
% 
% V_targetDirAvg = zeros(L,length(Dirs));
% %% cells pref 90
% for i = 1:length(Dirs)
%     ind = find(DirectionDeg_success_ind == Dirs(i));
%     ind = intersect(V_Success_ind,ind);
%     V_targetDirAvg90(:,i) = mean(mean(dFoverF_aroundTarget(:,cellsPrefNinety,ind),3),2);
% end
% 
% for i = 1:length(Dirs)
%     direction = Dirs(i);
%     legendinfo{i} = num2str(direction);
% end
% 
% cMap = colormap(winter(length(Dirs)));
% figure;
% for i = 1:length(Dirs)
%     plot(V_targetDirAvg90((20:end),i),'Color',cMap(i,:))
%     hold on
% end
% vline(10,'c')
% hold on
% legend(legendinfo)
% ylabel('dF/F')
% xlabel('frames')
% title('average target response, cells pref 90')
% 
% %% cells pref 0
% for i = 1:length(Dirs)
%     ind = find(DirectionDeg_success_ind == Dirs(i));
%     ind = intersect(V_Success_ind,ind);
%     V_targetDirAvg0(:,i) = mean(mean(dFoverF_aroundTarget(:,cellsPrefZero,ind),3),2);
% end
% 
% for i = 1:length(Dirs)
%     direction = Dirs(i);
%     legendinfo{i} = num2str(direction);
% end
% 
% cMap = colormap(winter(length(Dirs)));
% figure;
% for i = 1:length(Dirs)
%     plot(V_targetDirAvg0((20:end),i),'Color',cMap(i,:))
%     hold on
% end
% vline(10,'c')
% hold on
% legend(legendinfo)
% ylabel('dF/F')
% xlabel('frames')
% title('average target response, cells pref 0')
%     
% 
% 
