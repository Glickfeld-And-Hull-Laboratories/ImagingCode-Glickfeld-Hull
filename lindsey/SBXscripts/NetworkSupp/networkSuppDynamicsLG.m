
load('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\celine\Analysis\2p_analysis\suppressionDataForLG_033122\supressionDataforLG.mat')

%% plot average data
szs = unique(conditions(:,1));
cons= unique(conditions(:,2));
nSz = length(szs);
nCon =length(cons);
cond = zscore(conditions);
cond_int = [cond cond(:,1).*cond(:,2)];
nCond = nSz*nCon;

nFrames = size(meanTC_byConditionLG,3);
conditions_s = reshape(permute(reshape(conditions,[nSz nCon 2]),[2 1 3]),[nSz*nCon 2]);
figure;
for i = 1:nCond
    meanTC_byConditionLG_s = reshape(permute(reshape(squeeze(mean(meanTC_byConditionLG,1)),[nSz nCon nFrames]),[2 1 3]),[nSz*nCon nFrames]);
    subplot(5,4,i)
    plot(squeeze(meanTC_byConditionLG_s(i,:)))
    ylim([-0.015 0.07])
    %title(num2str(conditions_s(i,:)))
    
end
print('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeByCon\SizeByCon.pdf','-dpdf','-fillpage')

%% response window analysis
resp_wins(1) = {62:63};
resp_wins(2) = {67:70};
resp_wins(3) = {76:106};
resp_win_str = {'early','middle','late'};
avgResponseWinsLG = squeeze(mean(responseWinsLG,1));

%plot dF/F by response amplitude and stats
figure;
for i = 1:3
    subplot(4,3,i)
    imagesc(reshape(avgResponseWinsLG(:,i),[5 4]))
    colormap redblue
    set(gca, 'Xtick', 1:nCon, 'XtickLabel',cons,'Ytick', 1:nSz, 'YtickLabel',szs )
    axis square
    clim([-0.05 0.05])
    title(resp_win_str{i})
    mdl = fitlm(cond_int,avgResponseWinsLG(:,i));
    subplot(4,3,i+3)
    text(0,1,['slope = ' num2str(chop(mdl.Coefficients.Estimate(2),2))])
    text(0,0,['p = '  num2str(chop(mdl.Coefficients.pValue(2),2))])
    title('Size')
    axis off
    ylim([-1 2])
    subplot(4,3,i+6)
    text(0,1,['slope = ' num2str(chop(mdl.Coefficients.Estimate(3),2))])
    text(0,0,['p = '  num2str(chop(mdl.Coefficients.pValue(3),2))])
    title('Contrast')
    axis off
    ylim([-1 2])
    subplot(4,3,i+9)
    text(0,1,['slope = ' num2str(chop(mdl.Coefficients.Estimate(4),2))])
    text(0,0,['p = '  num2str(chop(mdl.Coefficients.pValue(4),2))])
    title('Interaction')
    axis off
    ylim([-1 2])
end
print('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeByCon\SizeByConByInt_RespWindow.pdf','-dpdf','-fillpage')

%% PCA on 20 conditions 
avgMeanTC_byConditionLG = squeeze(mean(meanTC_byConditionLG,1));
[u s v] = pca(avgMeanTC_byConditionLG,'Centered',false);
s_s = reshape(s,[5 4 20]);
%PCxPC plots
figure;
start = 1;
for i = 1:5
    for j = 1:5
        subplot(5,5,start)
        if i>j
            for ii = 1:4
                scatter(s_s(:,ii,i),s_s(:,ii,j))
                hold on
            end
            xlabel(['PC' num2str(i)])
            ylabel(['PC' num2str(j)])
        elseif i==j
            plot(u(:,i))
            title(['PC' num2str(i)])
        elseif i<j
            for ii = 1:5
                scatter(s_s(ii,:,i),s_s(ii,:,j))
                hold on
            end
            xlabel(['PC' num2str(i)])
            ylabel(['PC' num2str(j)])
        end
        start = start+1;
    end
end

%PC weights for size x contrast and stats
figure;
for i=1:5
    subplot(5,5,i)
    plot(u(:,i))
    ylim([-0.2 0.5])
    title(['PC' num2str(i)])
    subplot(5,5,i+5)
    imagesc(s_s(:,:,i))
    clim([-0.1 0.1])
    colormap redblue
    set(gca, 'Xtick', 1:nCon, 'XtickLabel',cons,'Ytick', 1:nSz, 'YtickLabel',szs )
    subplot(5,5,i+10)
    mdl = fitlm(cond_int,s(:,i));
    text(0,1,['slope = ' num2str(chop(mdl.Coefficients.Estimate(2),2))])
    text(0,0,['p = '  num2str(chop(mdl.Coefficients.pValue(2),2))])
    title('Size')
    axis off
    ylim([-1 2])
    subplot(5,5,i+15)
    text(0,1,['slope = ' num2str(chop(mdl.Coefficients.Estimate(3),2))])
    text(0,0,['p = '  num2str(chop(mdl.Coefficients.pValue(3),2))])
    title('Contrast')
    axis off
    ylim([-1 2])
    subplot(5,5,i+20)
    text(0,1,['slope = ' num2str(chop(mdl.Coefficients.Estimate(4),2))])
    text(0,0,['p = '  num2str(chop(mdl.Coefficients.pValue(4),2))])
    title('Interaction')
    axis off
    ylim([-1 2])
end
print('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeByCon\SizeByConByInt_PCA.pdf','-dpdf','-fillpage')

%% Locomotion analysis

avgByCond_Stat = squeeze(nanmean(TC_byConditionStatLG,1))';
avgByCond_Loc = squeeze(nanmean(TC_byConditionLocLG,1))';

statPCA = (u'*avgByCond_Stat)';
locPCA = (u'*avgByCond_Loc)';
statPCA_s = reshape(statPCA,[5 4 20]);
locPCA_s = reshape(locPCA,[5 4 20]);

figure;
for i=1:5
    subplot(5,5,i)
    plot(u(:,i))
    ylim([-0.2 0.5])
    title(['PC' num2str(i)])
    subplot(5,5,i+5)
    imagesc(statPCA_s(:,:,i))
    clim([-0.1 0.1])
    colormap redblue
    set(gca, 'Xtick', 1:nCon, 'XtickLabel',cons,'Ytick', 1:nSz, 'YtickLabel',szs )
    subplot(5,5,i+10)
    mdl = fitlm(cond_int,statPCA(:,i));
    text(0,1,['slope = ' num2str(chop(mdl.Coefficients.Estimate(2),2))])
    text(0,0,['p = '  num2str(chop(mdl.Coefficients.pValue(2),2))])
    title('Size')
    axis off
    ylim([-1 2])
    subplot(5,5,i+15)
    text(0,1,['slope = ' num2str(chop(mdl.Coefficients.Estimate(3),2))])
    text(0,0,['p = '  num2str(chop(mdl.Coefficients.pValue(3),2))])
    title('Contrast')
    axis off
    ylim([-1 2])
    subplot(5,5,i+20)
    text(0,1,['slope = ' num2str(chop(mdl.Coefficients.Estimate(4),2))])
    text(0,0,['p = '  num2str(chop(mdl.Coefficients.pValue(4),2))])
    title('Interaction')
    axis off
    ylim([-1 2])
end
sgtitle('Stationary')
print('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeByCon\SizeByConByInt_PCA_Stat.pdf','-dpdf','-fillpage')


figure;
for i=1:5
    subplot(5,5,i)
    plot(u(:,i))
    ylim([-0.2 0.5])
    title(['PC' num2str(i)])
    subplot(5,5,i+5)
    imagesc(locPCA_s(:,:,i))
    clim([-0.1 0.1])
    colormap redblue
    set(gca, 'Xtick', 1:nCon, 'XtickLabel',cons,'Ytick', 1:nSz, 'YtickLabel',szs )
    subplot(5,5,i+10)
    mdl = fitlm(cond_int,locPCA(:,i));
    text(0,1,['slope = ' num2str(chop(mdl.Coefficients.Estimate(2),2))])
    text(0,0,['p = '  num2str(chop(mdl.Coefficients.pValue(2),2))])
    title('Size')
    axis off
    ylim([-1 2])
    subplot(5,5,i+15)
    text(0,1,['slope = ' num2str(chop(mdl.Coefficients.Estimate(3),2))])
    text(0,0,['p = '  num2str(chop(mdl.Coefficients.pValue(3),2))])
    title('Contrast')
    axis off
    ylim([-1 2])
    subplot(5,5,i+20)
    text(0,1,['slope = ' num2str(chop(mdl.Coefficients.Estimate(4),2))])
    text(0,0,['p = '  num2str(chop(mdl.Coefficients.pValue(4),2))])
    title('Interaction')
    axis off
    ylim([-1 2])
end
sgtitle('Locomotion')
print('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeByCon\SizeByConByInt_PCA_Loc.pdf','-dpdf','-fillpage')

%% Single cell analysis
nCells = size(meanTC_byConditionLG,1);
cellPCA = zeros(nCond,nCond,nCells);
cellPCA_stat = zeros(nCond,nCond,nCells);
cellPCA_loc = zeros(nCond,nCond,nCells);
for i = 1:nCells
    cellPCA(:,:,i) = (u'*squeeze(meanTC_byConditionLG(i,:,:))')';
    cellPCA_stat(:,:,i) = (u'*squeeze(TC_byConditionStatLG(i,:,:))')';
    cellPCA_loc(:,:,i) = (u'*squeeze(TC_byConditionLocLG(i,:,:))')';
end
cellPCA_s = reshape(cellPCA,[nSz nCon nCond nCells]);
cellSlopes = zeros(3,5,nCells);
cellPs = zeros(3,5,nCells);
cellPCA_stat_s = reshape(cellPCA_stat,[nSz nCon nCond nCells]);
cellSlopes_stat = zeros(3,5,nCells);
cellPs_stat = zeros(3,5,nCells);
cellPCA_loc_s = reshape(cellPCA_loc,[nSz nCon nCond nCells]);
cellSlopes_loc = zeros(3,5,nCells);
cellPs_loc = zeros(3,5,nCells);
doPlot = 0;
for i = 1:nCells
    for ii = 1:5
        mdl = fitlm(cond_int,cellPCA(:,ii,i));
        cellSlopes(:,ii,i) = mdl.Coefficients.Estimate(2:4);
        cellPs(:,ii,i) = mdl.Coefficients.pValue(2:4);
        mdl = fitlm(cond_int,cellPCA_stat(:,ii,i));
        cellSlopes_stat(:,ii,i) = mdl.Coefficients.Estimate(2:4);
        cellPs_stat(:,ii,i) = mdl.Coefficients.pValue(2:4);
        mdl = fitlm(cond_int,cellPCA_loc(:,ii,i));
        cellSlopes_loc(:,ii,i) = mdl.Coefficients.Estimate(2:4);
        cellPs_loc(:,ii,i) = mdl.Coefficients.pValue(2:4);
    end
end

pIX = cellPs<0.05;
anova_str = {'Size','Contrast','Int'};
leg_str_pc = {'Neither','PC1','PC2', 'Both'};
leg_str_st = {'Neither','Size','Con', 'Both'};
figure;
for i = 1:3
    subplot(3,3,i)
    ind1=find(pIX(i,1,:));
    ind2=find(pIX(i,2,:));
    ind_u1 = setdiff(ind1,ind2);
    ind_u2 = setdiff(ind2,ind1);
    ind_12 = intersect(ind2,ind1);
    ind_n = setdiff(1:nCells,[ind1; ind2]);
    scatter(squeeze(cellSlopes(i,1,ind_n)),squeeze(cellSlopes(i,2,ind_n)))
    hold on
    scatter(squeeze(cellSlopes(i,1,ind_u1)),squeeze(cellSlopes(i,2,ind_u1)))
    scatter(squeeze(cellSlopes(i,1,ind_u2)),squeeze(cellSlopes(i,2,ind_u2)))
    scatter(squeeze(cellSlopes(i,1,ind_12)),squeeze(cellSlopes(i,2,ind_12)))
    axis square
    xlabel('PC1 slope')
    ylabel('PC2 slope')
    title(anova_str{i})
    xlim([-0.2 0.2])
    ylim([-0.2 0.2])
    if i == 3
        legend(leg_str_pc,'location','southwest')
    end
    subplot(3,3,i+3)
    pie([length(ind_n)./nCells length(ind_u1)./nCells length(ind_u2)./nCells length(ind_12)./nCells],leg_str)
    subplot(3,3,i+6)
    ind_s=find(sum(pIX(1,i,:),1));
    ind_c=find(sum(pIX(2,i,:),1));
    ind_us = setdiff(ind_s,ind_c);
    ind_uc = setdiff(ind_c,ind_s);
    ind_sc = intersect(ind_s,ind_c);
    ind_n = setdiff(1:nCells,[ind_s; ind_c]);
    scatter(squeeze(cellSlopes(1,i,ind_n)),squeeze(cellSlopes(2,i,ind_n)))
    hold on
    scatter(squeeze(cellSlopes(1,i,ind_us)),squeeze(cellSlopes(2,i,ind_us)))
    scatter(squeeze(cellSlopes(1,i,ind_uc)),squeeze(cellSlopes(2,i,ind_uc)))
    scatter(squeeze(cellSlopes(1,i,ind_sc)),squeeze(cellSlopes(2,i,ind_sc)))
    xlabel('Size slope')
    ylabel('Contrast slope')
    xlim([-0.2 0.2])
    ylim([-0.2 0.2])
    axis square
    if i == 3
        legend(leg_str_st,'location','southwest')
    end
    title(['PC' num2str(i)])
end
print('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeByCon\SingleCellSlopes_SizeConInt.pdf','-dpdf','-fillpage')

for i = 1:3
    ind1=find(pIX(i,1,:));
    ind2=find(pIX(i,2,:));
    ind_u1 = setdiff(ind1,ind2);
    ind_u2 = setdiff(ind2,ind1);
    ind_12 = intersect(ind2,ind1);
    n_str{1} = [leg_str_pc{2} '- n = ' num2str(length(ind_u1))];
    n_str{2} = [leg_str_pc{3} '- n = ' num2str(length(ind_u2))];
    n_str{3} = [leg_str_pc{4} '- n = ' num2str(length(ind_12))];
    figure;
    for ii = 1:nCond
        meanTC_byConditionLG_s =reshape(permute(reshape(squeeze(mean(meanTC_byConditionLG(ind_u1,:,:),1)),[nSz nCon nFrames]),[2 1 3]),[nSz*nCon nFrames]);
        subplot(5,4,ii)
        plot(squeeze(meanTC_byConditionLG_s(ii,:)))
        hold on
        ylim([-0.05 0.15])
        title(num2str(conditions_s(ii,:)))
    end
    for ii = 1:nCond
        meanTC_byConditionLG_s =reshape(permute(reshape(squeeze(mean(meanTC_byConditionLG(ind_u2,:,:),1)),[nSz nCon nFrames]),[2 1 3]),[nSz*nCon nFrames]);
        subplot(5,4,ii)
        plot(squeeze(meanTC_byConditionLG_s(ii,:)))
    end
    for ii = 1:nCond
        meanTC_byConditionLG_s =reshape(permute(reshape(squeeze(mean(meanTC_byConditionLG(ind_12,:,:),1)),[nSz nCon nFrames]),[2 1 3]),[nSz*nCon nFrames]);
        subplot(5,4,ii)
        plot(squeeze(meanTC_byConditionLG_s(ii,:)))
    end
    subplot(5,4,1)
    legend(n_str)
    sgtitle(anova_str{i})
    print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeByCon\SizeByConResp_sigPC1PC2_' anova_str{i} '.pdf'],'-dpdf','-fillpage')
end

figure;
for i = 1:5
    ind_s = find(pIX(1,i,:));
    ind_c = find(pIX(2,i,:));
    ind_n = setdiff(1:nCells,unique([ind_s; ind_c]));
    ind_so = setdiff(ind_s,ind_c);
    ind_co = setdiff(ind_c,ind_s);
    subplot(4,5,i)
    imagesc(mean(cellPCA_s(:,:,i,ind_s),4))
    clim([-0.15 0.15])
    axis square
    title(['PC' num2str(i)])
    if i == 1
        ylabel('Size')
    end
    subplot(4,5,i+5)
    imagesc(mean(cellPCA_s(:,:,i,ind_c),4))
    clim([-0.15 0.15])
    axis square
    if i == 1
        ylabel('Con')
    end
    subplot(4,5,i+10)
    imagesc(mean(cellPCA_s(:,:,i,ind_so),4))
    clim([-0.15 0.15])
    axis square
    if i == 1
        ylabel('Size only')
    end
    subplot(4,5,i+15)
    imagesc(mean(cellPCA_s(:,:,i,ind_co),4))
    clim([-0.15 0.15])
    axis square
    if i == 1
        ylabel('Con only')
    end
end
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeByCon\SizeByConWeights_byPC_byModGroup.pdf'],'-dpdf','-fillpage')

figure
for i = 1:5
    ind_s = find(pIX(1,i,:));
    ind_c = find(pIX(2,i,:));
    ind_n = setdiff(1:nCells,unique([ind_s; ind_c]));
    subplot(3,5,i)
    imagesc(mean(cellPCA_s(:,:,i,ind_n),4))
    clim([-0.15 0.15])
    axis square
    if i == 1
        ylabel('None')
    end
    w = reshape(cellPCA_s(:,:,i,ind_n),[nCon*nSz length(ind_n)]);
    h = ttest(w);
    title([num2str(sum(h)) '/' num2str(length(ind_n))])
    subplot(3,5,i+5)
    imagesc(mean(cellPCA_s(:,:,i,ind_n(find(h))),4))
    clim([-0.15 0.15])
    axis square
    if i == 1
        ylabel('W>0')
    end
    subplot(3,5,i+10)
    imagesc(mean(cellPCA_s(:,:,i,ind_n(find(~h))),4))
    clim([-0.15 0.15])
    axis square
    if i == 1
        ylabel('W=0')
    end
end
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeByCon\SizeByConWeights_byPC_nonModGroup.pdf'],'-dpdf','-fillpage')

% clear n_str
% for i = 1:3
%     ind1=find(pIX(i,1,:));
%     ind2=find(pIX(i,2,:));
%     n_str{1} = [leg_str_pc{2} '- n = ' num2str(length(ind1))];
%     n_str{2} = [leg_str_pc{3} '- n = ' num2str(length(ind2))];
%     figure;
%     for ii = 1:nCond
%         meanTC_byConditionLG_s =reshape(permute(reshape(squeeze(mean(meanTC_byConditionLG(ind1,:,:),1)),[nSz nCon nFrames]),[2 1 3]),[nSz*nCon nFrames]);
%         subplot(5,4,ii)
%         plot(squeeze(meanTC_byConditionLG_s(ii,:)))
%         hold on
%         ylim([-0.05 0.15])
%         title(num2str(conditions_s(ii,:)))
%     end
%     for ii = 1:nCond
%         meanTC_byConditionLG_s =reshape(permute(reshape(squeeze(mean(meanTC_byConditionLG(ind2,:,:),1)),[nSz nCon nFrames]),[2 1 3]),[nSz*nCon nFrames]);
%         subplot(5,4,ii)
%         plot(squeeze(meanTC_byConditionLG_s(ii,:)))
%     end
%     subplot(5,4,1)
%     legend(n_str)
%     sgtitle(anova_str{i})
% end

% %compare running- too noisy
% figure;
% for i = 1:3
%     subplot(3,3,i)
%     pIX = cellPs<0.05;
%     ind1=find(pIX(i,1,:));
%     ind2=find(pIX(i,2,:));
%     ind_u1 = setdiff(ind1,ind2);
%     ind_u2 = setdiff(ind2,ind1);
%     ind_12 = intersect(ind2,ind1);
%     ind_n = setdiff(1:nCells,[ind1; ind2]);
%     scatter(squeeze(cellSlopes(i,1,ind_n)),squeeze(cellSlopes(i,2,ind_n)))
%     hold on
%     scatter(squeeze(cellSlopes(i,1,ind_u1)),squeeze(cellSlopes(i,2,ind_u1)))
%     scatter(squeeze(cellSlopes(i,1,ind_u2)),squeeze(cellSlopes(i,2,ind_u2)))
%     scatter(squeeze(cellSlopes(i,1,ind_12)),squeeze(cellSlopes(i,2,ind_12)))
%     axis square
%     xlabel('PC1 slope')
%     ylabel('PC2 slope')
%     title(anova_str{i})
%     xlim([-1 1])
%     ylim([-1 1])
%     if i == 1
%         legend(leg_str_pc)
%     end
%     subplot(3,3,i+3)
%     pIX = cellPs_stat<0.05;
%     ind1=find(pIX(i,1,:));
%     ind2=find(pIX(i,2,:));
%     ind_u1 = setdiff(ind1,ind2);
%     ind_u2 = setdiff(ind2,ind1);
%     ind_12 = intersect(ind2,ind1);
%     ind_n = setdiff(1:nCells,[ind1; ind2]);
%     scatter(squeeze(cellSlopes_stat(i,1,ind_n)),squeeze(cellSlopes_stat(i,2,ind_n)))
%     hold on
%     scatter(squeeze(cellSlopes_stat(i,1,ind_u1)),squeeze(cellSlopes_stat(i,2,ind_u1)))
%     scatter(squeeze(cellSlopes_stat(i,1,ind_u2)),squeeze(cellSlopes_stat(i,2,ind_u2)))
%     scatter(squeeze(cellSlopes_stat(i,1,ind_12)),squeeze(cellSlopes_stat(i,2,ind_12)))
%     axis square
%     xlabel('PC1 slope (Stat)')
%     ylabel('PC2 slope (Stat)')
%     title(anova_str{i})
%     xlim([-1 1])
%     ylim([-1 1])
%     subplot(3,3,i+6)
%     pIX = cellPs_loc<0.05;
%     ind1=find(pIX(i,1,:));
%     ind2=find(pIX(i,2,:));
%     ind_u1 = setdiff(ind1,ind2);
%     ind_u2 = setdiff(ind2,ind1);
%     ind_12 = intersect(ind2,ind1);
%     ind_n = setdiff(1:nCells,[ind1; ind2]);
%     scatter(squeeze(cellSlopes(i,1,ind_n)),squeeze(cellSlopes_loc(i,2,ind_n)))
%     hold on
%     scatter(squeeze(cellSlopes_loc(i,1,ind_u1)),squeeze(cellSlopes_loc(i,2,ind_u1)))
%     scatter(squeeze(cellSlopes_loc(i,1,ind_u2)),squeeze(cellSlopes_loc(i,2,ind_u2)))
%     scatter(squeeze(cellSlopes_loc(i,1,ind_12)),squeeze(cellSlopes_loc(i,2,ind_12)))
%     axis square
%     xlabel('PC1 slope (Loc)')
%     ylabel('PC2 slope (Loc)')
%     title(anova_str{i})
%     xlim([-1 1])
%     ylim([-1 1])
% end
    

cellmat = zeros(5,5,3);
for ix = 1:3
    for i = 1:5
        for j = 1:5
            indi = find(pIX(ix,i,:));
            indj = find(pIX(ix,j,:));
            if i == j
                cellmat(i,j,ix) = length(indi)./nCells;
            else
                cellmat(i,j,ix) = length(intersect(indi,indj))./nCells;
            end
        end
    end
end
figure; 
for i = 1:3
subplot(1,3,i)
imagesc(cellmat(:,:,i))
axis square
clim([0 0.8])
title(anova_str{i})
xlabel('PC')
ylabel('PC')
end
print('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeByCon\FractSigCells_SizeConInt.pdf','-dpdf','-fillpage')


ind_s = find(abs(cellSlopes(1,1,:))>0.08);
ind_1 = find(pIX(1,1,:));
ind_use = intersect(ind_s,ind_1);
for i = ind_use'
    conditions_s = reshape(permute(reshape(conditions,[nSz nCon 2]),[2 1 3]),[nSz*nCon 2]);
    figure;
    for ii = 1:nCond
        meanTC_byConditionLG_s =reshape(permute(reshape(squeeze(meanTC_byConditionLG(i,:,:)),[nSz nCon nFrames]),[2 1 3]),[nSz*nCon nFrames]);
        subplot(5,4,ii)
        plot(squeeze(meanTC_byConditionLG_s(ii,:)))
        ylim([-0.1 0.4])
        title(num2str(conditions_s(ii,:)))
    end
    sgtitle(num2str(i))

    figure;
    for ii = 1:5
        subplot(5,5,ii)
        plot(u(:,ii))
        ylim([-0.2 0.5])
        title(['PC' num2str(ii)])
        subplot(5,5,ii+5)
        imagesc(cellPCA_s(:,:,ii,i))
        clim([-0.1 0.1])
        colormap redblue
        set(gca, 'Xtick', 1:nCon, 'XtickLabel',cons,'Ytick', 1:nSz, 'YtickLabel',szs )
        subplot(5,5,ii+10)
        text(0,1,['slope = ' num2str(chop(cellSlopes(1,ii,i),2))])
        text(0,0,['p = '  num2str(chop(cellPs(1,ii,i),2))])
        title('Size')
        axis off
        ylim([-1 2])
        subplot(5,5,ii+15)
        text(0,1,['slope = ' num2str(chop(cellSlopes(2,ii,i),2))])
        text(0,0,['p = '  num2str(chop(cellPs(2,ii,i),2))])
        title('Contrast')
        axis off
        ylim([-1 2])
        subplot(5,5,ii+20)
        text(0,1,['slope = ' num2str(chop(cellSlopes(3,ii,i),2))])
        text(0,0,['p = '  num2str(chop(cellPs(3,ii,i),2))])
        title('Interaction')
        axis off
        ylim([-1 2])    
    end
sgtitle(num2str(i))
end

cellPCA_corr = zeros(5,nCells);
cellPCA_allcorr = zeros(1,nCells);
cellPCA_all = reshape(cellPCA(:,1:5,:),[20*5 nCells]);
s_all = reshape(s(:,1:5),[20*5 1]);
for iC = 1:nCells
    cellPCA_allcorr(1,iC) = triu2vec(corrcoef(cellPCA_all(:,iC),s_all));
    for i = 1:5
        cellPCA_corr(i,iC) = triu2vec(corrcoef(cellPCA(:,i,iC),s(:,i)));
    end
end
ind = find(sum(sum(pIX,1),2));
figure; plot(cellPCA_corr(:,ind))
ylabel('Correlation with average')
xlabel('PC')
set(gca,'Xtick',1:5);
print('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeByCon\SingleCell_CorrWithAvg.pdf','-dpdf','-fillpage')

figure
start = 1;
for i = 1:5
    for j = 1:5
        subplot(5,5,start)
        if i == j   
            scatter(cellPCA_allcorr(1,ind), cellPCA_corr(i,ind))
        else
            scatter(cellPCA_corr(i,ind), cellPCA_corr(j,ind))
        end
        xlim([-1 1])
        ylim([-1 1])
        axis square
        start = start+1;
    end
end

%% PCA on 20 conditions for all cells
sz = size(meanTC_byConditionLG);
meanTC_byConditionLG_allCells = reshape(meanTC_byConditionLG,[sz(1).*sz(2) sz(3)]);
meanTC_byConditionLG_allCells_z = meanTC_byConditionLG_allCells./std(meanTC_byConditionLG_allCells(:,20:60),[],2);
ind = [];
for i = 1:sz(1).*sz(2)
    if find(abs(meanTC_byConditionLG_allCells_z(i,60:end))>3)
        ind = [ind i];
    end
end
[u s v] = pca(meanTC_byConditionLG_allCells_z(ind,:),'Centered',false);
figure;
for i = 1:16
subplot(4,4,i)
plot(u(:,i));
end
avgPCA = (u'*avgMeanTC_byConditionLG')';
avgPCA_s = reshape(avgPCA(:,1:20),[5 4 20]);

figure;
for i=1:5
    subplot(5,5,i)
    plot(u(:,i))
    ylim([-0.2 0.5])
    title(['PC' num2str(i)])
    subplot(5,5,i+5)
    imagesc(avgPCA_s(:,:,i))%./max(max(avgPCA_s(:,:,i),[],1),[],2))
    clim([-0.1 0.1])
    colormap redblue
    set(gca, 'Xtick', 1:nCon, 'XtickLabel',cons,'Ytick', 1:nSz, 'YtickLabel',szs )
    subplot(5,5,i+10)
    mdl = fitlm(cond_int,avgPCA(:,i));
    text(0,1,['slope = ' num2str(chop(mdl.Coefficients.Estimate(2),2))])
    text(0,0,['p = '  num2str(chop(mdl.Coefficients.pValue(2),2))])
    title('Size')
    axis off
    ylim([-1 2])
    subplot(5,5,i+15)
    text(0,1,['slope = ' num2str(chop(mdl.Coefficients.Estimate(3),2))])
    text(0,0,['p = '  num2str(chop(mdl.Coefficients.pValue(3),2))])
    title('Contrast')
    axis off
    ylim([-1 2])
    subplot(5,5,i+20)
    text(0,1,['slope = ' num2str(chop(mdl.Coefficients.Estimate(4),2))])
    text(0,0,['p = '  num2str(chop(mdl.Coefficients.pValue(4),2))])
    title('Interaction')
    axis off
    ylim([-1 2])
end

%% ICA
PCuse = 1:10;
nIC = 5;
termtol = 0.00001; % termination tolerance
maxrounds = 1000; % #of iteration
sclass = class(u);

mixedsig = u(:,PCuse);
    
[icasig, A, W] = fastica(mixedsig','numOfIC',nIC);

icskew = skewness(icasig');
[icskew, ICord] = sort(abs(icskew), 'descend');
ica_A = A(:,ICord);
icasig = icasig(ICord,:);
ica_W = W(ICord,:);

figure;
[n n2] = subplotn(nIC);
for iC = 1:nIC
    subplot(n,n2,iC)
    plot(icasig(iC,:))
end
% subplot(n,n2,iC+1)
% plot(ica_sig')
figure;
ICs_w = reshape(ica_A,[5 4 nIC]);
for iC = 1:nIC
    subplot(n,n2,iC)
    imagesc(ICs_w(:,:,iC))
    colormap redblue
end

figure;
ICs_f = reshape(ica_filters,[5 4 nIC]);
for iC = 1:nIC
    subplot(n,n2,iC)
    imagesc(ICs_f(:,:,iC))
    colormap redblue
end

%% NNMF
[w h] = nnmf(mixedsig,5);
[icasig, A, W] = fastica(data','numOfIC',nIC);


%% Basis functions
t = 1:(123-62);
tdfaste = 7;
tdfasti = 7;
trfast = 4;
tdslow = 30;
trslow = 10;
allBuf = zeros(1,62);
Ibuf = zeros(1,2);
fastI = [allBuf -exp(-t./tdfasti)+exp(-t./trfast)];
fastI = [Ibuf fastI(1:end-2)];
fastE = [allBuf exp(-t./tdfaste)-exp(-t./trfast)];
Ibuf = zeros(1,4);
slowI = [allBuf -exp(-t./tdslow)+exp(-t./trslow)];
slowI = [Ibuf slowI(1:end-4)];
Ibuf = zeros(1,3);
slowE = [allBuf exp(-t./tdslow)-exp(-t./trslow)];
slowE = [Ibuf slowE(1:end-3)];

% figure;
% subplot(2,2,1)
% plot(t,fastE)
% subplot(2,2,2)
% plot(t,fastI)
% subplot(2,2,3)
% plot(t,slowE)
% subplot(2,2,4)
% plot(t,slowI)

figure;
A = .5;
B = .12;
C = .1;
D = .4;
% tot = A.*fastE;
% plot(t,tot)
% hold on
% tot = A.*fastE+B.*fastI;
% plot(t,tot)
% tot = A.*fastE+B.*fastI+C.*slowE+D.*slowI;
% plot(t,tot)

te = -62:123-62-1;
tot = A.*fastE;
plot(te,tot) 
hold on
tot = A.*fastE+B.*fastI;
plot(te,tot) 
plot(te,meanTC_byConditionLG_s(4,:))

A.*fastE+B.*fastI+C.*slowE+D.*slowI;


AmpA_guess = 0.5.*ones(1,20);
AmpB_guess = 0.5.*ones(1,20);
AmpC_guess = 0.5.*ones(1,20);
AmpD_guess = 0.5.*ones(1,20);
fastRiseTau_guess = 3;
slowRiseTau_guess = 10;
fastDecayTau_guess = 7;
slowDecayTau_guess = 30;

AmpA_lb = zeros(1,20);
AmpB_lb = zeros(1,20);
AmpC_lb = zeros(1,20);
AmpD_lb = zeros(1,20);
fastRiseTau_lb = 1;
slowRiseTau_lb = 3;
fastDecayTau_lb = 4;
slowDecayTau_lb = 10;

AmpA_ub = 5.*ones(1,20);
AmpB_ub = 5.*ones(1,20);
AmpC_ub = 5.*ones(1,20);
AmpD_ub = 5.*ones(1,20);
fastRiseTau_ub = 10;
slowRiseTau_ub = 30;
fastDecayTau_ub = 30;
slowDecayTau_ub = 100;

options = optimset('MaxFunEvals',inf,'MaxIter',100000);
[out,fval,success] = fminsearchbnd(@ExPlusIn,...
    [AmpA_guess, AmpB_guess,AmpC_guess,AmpD_guess, fastRiseTau_guess, slowRiseTau_guess,fastDecayTau_guess,slowDecayTau_guess],...
    [AmpA_lb, AmpB_lb,AmpC_lb,AmpD_lb, fastRiseTau_lb, slowRiseTau_lb,fastDecayTau_lb,slowDecayTau_lb],...
    [AmpA_ub, AmpB_ub,AmpC_ub,AmpD_ub, fastRiseTau_ub, slowRiseTau_ub,fastDecayTau_ub,slowDecayTau_ub], options);

function miaosse = ExPlusIn(in)
        % pull out the slope and intercept
        AmpA = in(1);
        AmpB = in(2);
        AmpC = in(3);
        AmpD = in(4);
        fastRiseTau = in(5);
        slowRiseTau = in(6);
        fastDecayTau = in(7);
        slowDecayTau = in(8);
       
        allBuf = zeros(1,62);
        Ibuf = zeros(1,2);
        fastI = [allBuf -exp(-t./fastDecayTau)+exp(-t./fastRiseTau)];
        fastI = [Ibuf fastI(1:end-2)];
        fastE = [allBuf exp(-t./fastDecayTau)-exp(-t./fastRiseTau)];
        Ibuf = zeros(1,4);
        slowI = [allBuf -exp(-t./slowDecayTau)+exp(-t./slowRiseTau)];
        slowI = [Ibuf slowI(1:end-4)];
        Ibuf = zeros(1,3);
        slowE = [allBuf exp(-t./slowDecayTau)-exp(-t./slowRiseTau)];
        slowE = [Ibuf slowE(1:end-3)];
        y_fit = AmpA.*fastE+AmpB.*fastI+AmpC.*slowE+AmpD.*slowI;
        
        residuals = data - y_fit;
        miaosse = sum(residuals.^2);
    end