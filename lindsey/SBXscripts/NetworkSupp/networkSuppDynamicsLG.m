
load('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\celine\Analysis\2p_analysis\suppressionDataForLG_033122\suppressionDataForLG.m')

%% plot average data
szs = unique(conditions(:,1));
cons= unique(conditions(:,2));
nSz = length(szs);
nCon =length(cons);
cond = zscore(conditions);
cond_int = [cond cond(:,1).*cond(:,2)];

figure; 
for i = 1:20
    subplot(4,5,i)
    plot(squeeze(mean(meanTC_byConditionLG(:,i,:),1)))
    ylim([-0.01 0.06])
    title(num2str(conditions(i,:)))
end
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
