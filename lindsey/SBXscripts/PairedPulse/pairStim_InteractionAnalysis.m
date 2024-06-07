close all 
clear all
date = '240530';
ImgFolder = {'002'};
mouse = 'i1395';
nrun = size(ImgFolder,2);
run_str = catRunName(ImgFolder, nrun);

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))

%%
nCells = size(dfof_resp_one,1);
tStimOne = celleqel2mat_padded(input.tstimOne);
stimOne = unique(tStimOne);
nStim = length(stimOne);
tStimTwo = celleqel2mat_padded(input.tstimTwo);

resp_cell = cell(1,nStim);
h_stim = zeros(nStim,nCells);
p_stim = zeros(nStim,nCells);
dfof_resp_stim = zeros(nStim,nCells);
for iStim = 1:nStim
    ind = find(tStimOne==stimOne(iStim));
    resp_cell{iStim} = cat(3,dfof_resp_one(:,ind),dfof_resp_two(:,ind));
    [h_stim(iStim,:), p_stim(iStim,:)] = ttest(dfof_resp_one(:,ind)',dfof_base_one(:,ind)','tail','right','alpha',0.05./(nStim-1));
    dfof_resp_stim(iStim,:) = nanmean(dfof_resp_one(:,ind),2)';
end

dfof_all_stim = zeros(nCells,nStim,nStim,2);
for is1 = 1:nStim
    for is2 = 1:nStim
        ind = intersect(find(tStimOne == stimOne(is1)),find(tStimTwo == stimOne(is2)));
        dfof_all_stim(:,is1,is2,1) = nanmean(dfof_resp_one(:,ind),2);
        dfof_all_stim(:,is1,is2,2) = nanmean(dfof_resp_two(:,ind),2);
    end
end
adapt_all_stim = (dfof_all_stim(:,:,:,2)-dfof_all_stim(:,:,:,1))./dfof_all_stim(:,:,:,1);

ind_combo = cell(nStim,nStim);
ind_combo_n = zeros(nStim,nStim);
h_resp = dfof_resp_stim>0.05;
for i = 1:nStim
    temp_resp = dfof_all_stim(:,i,i,1)>0.05;
    ind_combo{i,i} = find(h_stim(i,:)+temp_resp'==2);
    ind_combo_n(i,i) = length(find(h_stim(i,:)+temp_resp'==2));
    for j = 1:nStim
        if i~=j
            ind_combo{i,j} = find(h_stim(i,:)+h_stim(j,:)+h_resp(i,:)+h_resp(j,:)==4);
            if length(ind_combo{i,j})>0
                ind_combo_n(i,j) = length(find(h_stim(i,:)+h_stim(j,:)+h_resp(i,:)+h_resp(j,:)==4));
            end
        end
    end
end
figure; imagesc(ind_combo_n)
stimstr = {'g1','g2','n1','n2'};
set(gca,'Xtick',[1:4],'Ytick',[1:4],'XtickLabel',stimstr,'YTickLabel',stimstr)
colorbar
axis square
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_cellNCombo.pdf']),'-dpdf','-bestfit')


figure;
for is = 1:nStim
    ind_cell = ind_combo{is,is};
    scatter(is.*ones(size(ind_cell)), adapt_all_stim(ind_cell,is,is),'oc')
    hold on
    errorbar(is, nanmean(adapt_all_stim(ind_cell,is,is),1), nanstd(adapt_all_stim(ind_cell,is,is),[],1)./sqrt(length(ind_cell)),'ok')
end
set(gca, 'XTick', [1:4], 'XTickLabel', {'0.08', '0.32', 'nat1', 'nat2'})
xlim([0 5])

allStim_comp = [];
stimOne_comp = [];
start = 1;
for is1 = 1:nStim
    ind1 = intersect(find(tStimOne == stimOne(is1)),find(tStimTwo == stimOne(is1)));
    for is2= 1:nStim
        if is2 ~= is1
            ind2 = intersect(find(tStimTwo == stimOne(is1)),find(tStimOne == stimOne(is2)));
            ind_cell = ind_combo{is1,is2};
            figure(1)
            subplot(4,3,start)
            errorbar([1 2],[nanmean(nanmean(dfof_resp_two(ind_cell,ind1))) nanmean(nanmean(dfof_resp_two(ind_cell,ind2)))],[nanstd(nanmean(dfof_resp_two(ind_cell,ind1),2),1)./sqrt(length(ind_cell)) nanstd(nanmean(dfof_resp_two(ind_cell,ind2),2),1)./sqrt(length(ind_cell))])
            hold on
            errorbar(1,nanmean(nanmean(dfof_resp_one(ind_cell,ind1))),nanstd(nanmean(dfof_resp_one(ind_cell,ind1),2),1)./sqrt(length(ind_cell)))
            ylim([0 0.5])
            title(['S2: ' stimstr{is1}])
            xlabel(['Stim 1'])
            set(gca, 'XTick', [1 2], 'XTickLabel', [is1 is2])
            xlim([0 3])
            figure(2)
            subplot(4,3,start)
            plot(nanmean(nanmean(tc_one_dfof(:,ind_cell,ind1),2),3))
            hold on
            plot(nanmean(nanmean(tc_two_dfof(:,ind_cell,ind1),2),3)-nanmean(nanmean(nanmean(tc_two_dfof(base_win,ind_cell,ind1),2),3)))
            plot(nanmean(nanmean(tc_two_dfof(:,ind_cell,ind2),2),3)-nanmean(nanmean(nanmean(tc_two_dfof(base_win,ind_cell,ind2),2),3)))
            title(['S2: ' stimstr{is1} '- n = ' num2str(length(ind_cell))])
            legend({stimstr{is1},stimstr{is1},stimstr{is2}},'location','northwest')
            start = start +1;
            ylim([-0.05 0.6])
            axis square
            allStim_comp = [allStim_comp; nanmean(dfof_resp_two(ind_cell,ind1),2) nanmean(dfof_resp_two(ind_cell,ind2),2) is2.*ones(length(ind_cell),1) is1.*ones(length(ind_cell),1)];
            stimOne_comp = [stimOne_comp; nanmean(dfof_resp_one(ind_cell,ind1),2) nanmean(dfof_resp_two(ind_cell,ind2),2) is2.*ones(length(ind_cell),1) is1.*ones(length(ind_cell),1)];
        end
    end
end
figure(1)
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgRespComp.pdf']),'-dpdf','-fillpage')
figure(2)
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgTCComp.pdf']),'-dpdf','-fillpage')

figure;
ind_grating = intersect(find(allStim_comp(:,3)<3),find(allStim_comp(:,4)<3));
subplot(2,3,1)
plot([allStim_comp(ind_grating,1:2)'],'-ok')
hold on
errorbar(1:2, nanmean(allStim_comp(ind_grating,1:2),1), nanstd(allStim_comp(ind_grating,1:2),[],1)./sqrt(length(ind_grating)),'-or')
[h p_grating] = ttest(allStim_comp(ind_grating,1),allStim_comp(ind_grating,2));
title(['Grating comp- p = ' num2str(chop(p_grating,3))])
xlim([0 3])
ylim([0 1])
ind_nat = intersect(find(allStim_comp(:,3)>2),find(allStim_comp(:,4)>2));
subplot(2,3,2)
plot([allStim_comp(ind_nat,1:2)'],'-ok')
hold on
errorbar(1:2, nanmean(allStim_comp(ind_nat,1:2),1), nanstd(allStim_comp(ind_nat,1:2),[],1)./sqrt(length(ind_nat)),'-or')
[h p_image] = ttest(allStim_comp(ind_nat,1),allStim_comp(ind_nat,2));
title(['Image comp- p = ' num2str(chop(p_image,3))])
xlim([0 3])
ylim([0 1])
ind_mix = [intersect(find(allStim_comp(:,3)<3),find(allStim_comp(:,4)>2)); intersect(find(allStim_comp(:,3)>2),find(allStim_comp(:,4)<3))];
subplot(2,3,3)
plot([allStim_comp(ind_mix,1:2)'],'-ok')
hold on
errorbar(1:2, nanmean(allStim_comp(ind_mix,1:2),1), nanstd(allStim_comp(ind_mix,1:2),[],1)./sqrt(length(ind_mix)),'-or')
[h p_mix] = ttest(allStim_comp(ind_mix,1),allStim_comp(ind_mix,2));
title(['Cross comp- p = ' num2str(chop(p_mix,3))])
xlim([0 3])
ylim([0 1])
ind_grating = intersect(find(allStim_comp(:,3)<3),find(allStim_comp(:,4)<3));
subplot(2,3,4)
plot([stimOne_comp(ind_grating,1:2)'],'-ok')
hold on
errorbar(1:2, nanmean(stimOne_comp(ind_grating,1:2),1), nanstd(stimOne_comp(ind_grating,1:2),[],1)./sqrt(length(ind_grating)),'-or')
[h p_one_grating] = ttest(stimOne_comp(ind_grating,1),stimOne_comp(ind_grating,2));
title(['Grating comp- p = ' num2str(chop(p_one_grating,3))])
xlim([0 3])
ylim([0 1])
ind_nat = intersect(find(stimOne_comp(:,3)>2),find(stimOne_comp(:,4)>2));
subplot(2,3,5)
plot([stimOne_comp(ind_nat,1:2)'],'-ok')
hold on
errorbar(1:2, nanmean(stimOne_comp(ind_nat,1:2),1), nanstd(stimOne_comp(ind_nat,1:2),[],1)./sqrt(length(ind_nat)),'-or')
[h p_one_image] = ttest(stimOne_comp(ind_nat,1),stimOne_comp(ind_nat,2));
title(['Image comp- p = ' num2str(chop(p_one_image,3))])
xlim([0 3])
ylim([0 1])
ind_mix = [intersect(find(stimOne_comp(:,3)<3),find(stimOne_comp(:,4)>2)); intersect(find(stimOne_comp(:,3)>2),find(stimOne_comp(:,4)<3))];
subplot(2,3,6)
plot([stimOne_comp(ind_mix,1:2)'],'-ok')
hold on
errorbar(1:2, nanmean(stimOne_comp(ind_mix,1:2),1), nanstd(stimOne_comp(ind_mix,1:2),[],1)./sqrt(length(ind_mix)),'-or')
[h p_one_mix] = ttest(stimOne_comp(ind_mix,1),stimOne_comp(ind_mix,2));
title(['Cross comp- p = ' num2str(chop(p_one_mix,3))])
xlim([0 3])
ylim([0 1])

print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_allCellsAllComps.pdf']),'-dpdf','-fillpage')

combo_resp = cell(nStim,nStim);
for is1 = 1:nStim
    for is2 = 1:nStim
        if is2>=is1
            combo_resp{is1,is2} = squeeze(dfof_all_stim(ind_combo{is1,is2},is1,is2,:));
        end
    end
end

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respComp.mat']),'allStim_comp','stimOne_comp','combo_resp','ind_combo','ind_combo_n','h_stim','dfof_all_stim','adapt_all_stim')

