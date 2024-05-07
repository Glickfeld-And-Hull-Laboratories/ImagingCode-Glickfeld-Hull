date = '240428';
ImgFolder = {'002'};
mouse = 'i1392';
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
        if i<j
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
    subplot(2,2,is)
    plot(squeeze(dfof_all_stim(ind_combo{is,is},is,:,1))')
end

figure;
for is = 1:nStim
    ind_cell = ind_combo{is,is};
    scatter(is.*ones(size(ind_cell)), adapt_all_stim(ind_cell,is,is),'oc')
    hold on
    errorbar(is, nanmean(adapt_all_stim(ind_cell,is,is),1), nanstd(adapt_all_stim(ind_cell,is,is),[],1)./sqrt(length(ind_cell)),'ok')
end
set(gca, 'XTick', [1:4], 'XTickLabel', {'0.08', '0.32', 'nat1', 'nat2'})
xlim([0 5])

start = 1;
for is1 = 1:nStim
    ind1 = intersect(find(tStimOne == stimOne(is1)),find(tStimTwo == stimOne(is1)));
    for is2= 1:nStim
        if is2>is1
            ind2 = intersect(find(tStimTwo == stimOne(is1)),find(tStimOne == stimOne(is2)));
            ind_cell = ind_combo{is1,is2};
            figure(1)
            subplot(2,3,start)
            errorbar([1 2],[nanmean(nanmean(dfof_resp_two(ind_cell,ind1))) nanmean(nanmean(dfof_resp_two(ind_cell,ind2)))],[nanstd(nanmean(dfof_resp_two(ind_cell,ind1),2),1)./sqrt(length(ind_cell)) nanstd(nanmean(dfof_resp_two(ind_cell,ind2),2),1)./sqrt(length(ind_cell))])
            hold on
            errorbar(1,nanmean(nanmean(dfof_resp_one(ind_cell,ind1))),nanstd(nanmean(dfof_resp_one(ind_cell,ind1),2),1)./sqrt(length(ind_cell)))
            ylim([0 0.5])
            title(['S2: ' stimstr{is1}])
            xlabel(['Stim 1'])
            set(gca, 'XTick', [1 2], 'XTickLabel', [is1 is2])
            xlim([0 3])
            figure(2)
            subplot(2,3,start)
            plot(nanmean(nanmean(tc_one_dfof(:,ind_cell,ind1),2),3))
            hold on
            plot(nanmean(nanmean(tc_two_dfof(:,ind_cell,ind1),2),3)-nanmean(nanmean(nanmean(tc_two_dfof(base_win,ind_cell,ind1),2),3)))
            plot(nanmean(nanmean(tc_two_dfof(:,ind_cell,ind2),2),3)-nanmean(nanmean(nanmean(tc_two_dfof(base_win,ind_cell,ind2),2),3)))
            title(['S2: ' stimstr{is1} '- n = ' num2str(length(ind_cell))])
            legend({stimstr{is1},stimstr{is1},stimstr{is2}},'location','northwest')
            start = start +1;
            ylim([-0.05 0.6])
            axis square

        end
    end
end
figure(1)
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgRespComp.pdf']),'-dpdf','-fillpage')
figure(2)
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgTCComp.pdf']),'-dpdf','-fillpage')

combo_resp = cell(nStim,nStim);
for is1 = 1:nStim
    for is2 = 1:nStim
        if is2>=is1
            combo_resp{is1,is2} = squeeze(dfof_all_stim(ind_combo{is1,is2},is1,is2,:));
        end
    end
end

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respComp.mat']),'combo_resp','ind_combo','ind_combo_n','h_stim','dfof_all_stim','adapt_all_stim')

figure;
for i = 1:nStim
    ind1 = find(tStimOne== stimOne(i));
    for ii = 1:nStim
        ind2 = find(tStimTwo==stimOne(ii));
        ind = intersect(ind1,ind2);
        subplot(2,2,i)
        plot(nanmean(nanmean(tc_one_dfof(:,find(h_resp(i,:)),ind),2),3))
        hold on
    end
end
