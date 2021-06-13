clear all
close all
clc
doRedChannel = 0;
ds = 'CrossOriRandDir_ExptList';
driver = 'SLC';
area_list = strvcat('V1','LM','AL','RL','PM');
narea = length(area_list);
eval(ds)
rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'RandDirSummary');
leg_str = cell(4,narea);
Zc_all_all = [];
Zp_all_all = [];
ZcZp_all_all = [];
Zc_pca_all = [];
Zp_pca_all = [];
f2overf1_all_all = [];
resp_ind_all_all = [];
f1resp_ind_all_all = [];
totCells = 0;
area_ind = [];
pca_area = [];
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir,['randDir_Summary_' area_list(iA,:) '_' driver '.mat']))
    areaSummary(iA).name = area_list(iA,:);
    areaSummary(iA).mice = unique(mouse_list,'rows');
    areaSummary(iA).nmice = size(unique(mouse_list,'rows'),1);
    areaSummary(iA).nexp = size(mouse_list,1);
    areaSummary(iA).totCells = length(stim_OSI_all);
    areaSummary(iA).respCells = length(resp_ind_all);
    figure(1)
    ind = resp_ind_all;
    leg_str{1,iA} = [area_list(iA,:) '- ' num2str(length(ind)) ' cells'];
    subplot(3,3,1)
    cdfplot(stim_OSI_all(ind))
    hold on
    subplot(3,3,2)
    cdfplot(stim_DSI_all(ind))
    hold on
    subplot(3,3,3)
    cdfplot(k_all(ind))
    hold on
    subplot(3,3,4)
    cdfplot(Zc_all(ind))
    hold on
    subplot(3,3,5)
    cdfplot(Zp_all(ind))
    hold on
    subplot(3,3,6)
    errorbar(mean(Zc_all(ind),2),mean(Zp_all(ind),2),std(Zc_all(ind),[],2)./sqrt(length(ind)),std(Zc_all(ind),[],2)./sqrt(length(ind)),std(Zp_all(ind),[],2)./sqrt(length(ind)),std(Zp_all(ind),[],2)./sqrt(length(ind)),'o')
    hold on
    ind_dsi = intersect(resp_ind_all,find(stim_DSI_all>0.5));
    leg_str{3,iA} = [area_list(iA,:) '- ' num2str(length(ind_dsi)) ' cells'];
    subplot(3,3,7)
    cdfplot(Zc_all(ind_dsi))
    hold on
    subplot(3,3,8)
    cdfplot(Zp_all(ind_dsi))
    hold on
    subplot(3,3,9)
    errorbar(mean(Zc_all(ind_dsi),2),mean(Zp_all(ind_dsi),2),std(Zc_all(ind_dsi),[],2)./sqrt(length(ind_dsi)),std(Zc_all(ind_dsi),[],2)./sqrt(length(ind_dsi)),std(Zp_all(ind_dsi),[],2)./sqrt(length(ind_dsi)),std(Zp_all(ind_dsi),[],2)./sqrt(length(ind_dsi)),'o')
    hold on
    
    figure(2)
    ind = intersect(resp_ind_all,find(f1_all>0.02));
    f1_ind = ind;
    leg_str{2,iA} = [area_list(iA,:) '- ' num2str(length(ind)) ' cells'];
    if ~isempty(ind)
        subplot(2,2,1)
        cdfplot(f2overf1_all(ind))
        hold on
        subplot(2,2,2)
        errorbar(mean(f1_all(ind),2),mean(f2_all(ind),2),std(f1_all(ind),[],2)./sqrt(length(ind)),std(f1_all(ind),[],2)./sqrt(length(ind)),std(f2_all(ind),[],2)./sqrt(length(ind)),std(f2_all(ind),[],2)./sqrt(length(ind)),'o')
        hold on
        subplot(2,2,3)
        errorbar(mean(Zc_all(ind),2),mean(f2overf1_all(ind),2),std(Zc_all(ind),[],2)./sqrt(length(ind)),std(Zc_all(ind),[],2)./sqrt(length(ind)),std(f2overf1_all(ind),[],2)./sqrt(length(ind)),std(f2overf1_all(ind),[],2)./sqrt(length(ind)),'o')
        hold on
        subplot(2,2,4)
        errorbar(mean(Zp_all(ind),2),mean(f2overf1_all(ind),2),std(Zp_all(ind),[],2)./sqrt(length(ind)),std(Zp_all(ind),[],2)./sqrt(length(ind)),std(f2overf1_all(ind),[],2)./sqrt(length(ind)),std(f2overf1_all(ind),[],2)./sqrt(length(ind)),'o')
        hold on
    else
        subplot(2,2,1)
        x = 0:3;
        plot(x,nan(size(x)))
        subplot(2,2,2)
        scatter(nan,nan)
        subplot(2,2,3)
        scatter(nan,nan)
        subplot(2,2,4)
        scatter(nan,nan)
    end
   
    figure(3)
    subplot(3,narea,iA)
    ind_h = intersect(resp_ind_all,find(stim_OSI_all>0.5));
    ind_l = intersect(resp_ind_all,find(stim_OSI_all<0.5));
    cdfplot(plaid_SI_all(ind_l))
    hold on
    cdfplot(plaid_SI_all(ind_h))
    xlim([-1 1])
    legend(['OSI<0.5- n =' num2str(length(ind_l))],['OSI>0.5- n =' num2str(length(ind_h))],'location','southeast')
    xlabel('Selectivity index')
    ylabel('Fraction of cells')
    title(area_list(iA,:))
    subplot(3,narea,iA+narea)
    cdfplot(f2overf1_all(ind_l))
    hold on
    cdfplot(f2overf1_all(ind_h))
    xlabel('F2/F1')
    ylabel('Fraction of cells')
    legend(['OSI<0.5- n =' num2str(length(ind_l))],['OSI>0.5- n =' num2str(length(ind_h))],'location','southeast')
    xlim([0 5])
    ind_h = intersect(resp_ind_all,find(plaid_SI_all>0));
    ind_l = intersect(resp_ind_all,find(plaid_SI_all<0));
    subplot(3,narea,iA+(narea.*2))
    cdfplot(f2overf1_all(ind_l))
    hold on
    cdfplot(f2overf1_all(ind_h))
    xlabel('F2/F1')
    ylabel('Fraction of cells')
    xlim([0 5])
    legend(['SI<0- n =' num2str(length(ind_l))],['SI>0- n =' num2str(length(ind_h))],'location','southeast')
    
    Zc_all_all = [Zc_all_all Zc_all];
    Zp_all_all = [Zp_all_all Zp_all];
    ZcZp_all = Zc_all-Zp_all;
    ZcZp_all_all = [ZcZp_all_all ZcZp_all];
    f2overf1_all_all = [f2overf1_all_all f2overf1_all];
    resp_ind_all_all = [resp_ind_all_all resp_ind_all+totCells];
    f1resp_ind_all_all = [f1resp_ind_all_all f1_ind+totCells];
    area_ind = [area_ind iA.*ones(size(Zc_all))];
    totCells = totCells+size(Zc_all,2);
    
    figure(4)
    subplot(2,2,1)
    errorbar(iA, mean(Zc_all(resp_ind_all),2),std(Zc_all(resp_ind_all),[],2)./sqrt(length(resp_ind_all)),'ok')
    hold on
    subplot(2,2,2)
    errorbar(iA, mean(Zp_all(resp_ind_all),2),std(Zp_all(resp_ind_all),[],2)./sqrt(length(resp_ind_all)),'ok')
    hold on
    subplot(2,2,3)
    errorbar(iA, mean(ZcZp_all(resp_ind_all),2),std(ZcZp_all(resp_ind_all),[],2)./sqrt(length(resp_ind_all)),'ok')
    hold on
    subplot(2,2,4)
    errorbar(iA, mean(f2overf1_all(f1_ind),2),std(f2overf1_all(f1_ind),[],2)./sqrt(length(f1_ind)),'ok')
    hold on
    
    figure(5)
    subplot(2,narea,iA)
    scatter(Zc_all(resp_ind_all),Zp_all(resp_ind_all),'ok')
    hold on
    Zc_use = intersect(resp_ind_all,intersect(find(Zc_all>1.28),find(Zc_all-Zp_all>1.28)));
    scatter(Zc_all(Zc_use),Zp_all(Zc_use),'ob')
    Zp_use = intersect(resp_ind_all,intersect(find(Zp_all>1.28),find(Zp_all-Zc_all>1.28)));
    scatter(Zc_all(Zp_use),Zp_all(Zp_use),'or')
    xlim([-4 8])
    ylim([-4 8])
    xlabel('Zc')
    ylabel('Zp')
    title(area_list(iA,:))
    axis square
    plotZcZpBorders
    subplot(2,narea,iA+narea)
    scatter(Zc_all(ind_dsi),Zp_all(ind_dsi),'ok')
    hold on
    Zc_use = intersect(ind_dsi,intersect(find(Zc_all>1.28),find(Zc_all-Zp_all>1.28)));
    scatter(Zc_all(Zc_use),Zp_all(Zc_use),'ob')
    Zp_use = intersect(ind_dsi,intersect(find(Zp_all>1.28),find(Zp_all-Zc_all>1.28)));
    scatter(Zc_all(Zp_use),Zp_all(Zp_use),'or')
    xlim([-4 8])
    ylim([-4 8])
    xlabel('Zc')
    ylabel('Zp')
    title([area_list(iA,:) '- DSI>0.5'])
    axis square
    plotZcZpBorders
    
%     load(fullfile(summaryDir,['randDir_PCA_Summary_' area_list(iA,:) '.mat']))
%     Zc_pca_all = [Zc_pca_all Zc_all];
%     Zp_pca_all = [Zp_pca_all Zp_all];
%     pca_area = [pca_area iA.*ones(size(Zc_all))];
%     figure(6)
%     subplot(2,1,1)
%     errorbar(iA,nanmean(Zc_all,2),nanstd(Zc_all,[],2)./sqrt(sum(~isnan(Zc_all))),'ok')
%     hold on
%     subplot(2,1,2)
%     errorbar(iA,nanmean(Zp_all,2),nanstd(Zp_all,[],2)./sqrt(sum(~isnan(Zp_all))),'ok')
%     hold on
end
figure(1)
subplot(3,3,1)
legend(leg_str{1,:})
xlabel('OSI')
title('')
subplot(3,3,2)
xlabel('DSI')
title('')
subplot(3,3,3)
xlabel('Kappa')
title('')
subplot(3,3,4)
xlabel('Zc')
xlim([-2 6])
title('All cells')
subplot(3,3,5)
xlabel('Zp')
xlim([-2 6])
title('')
subplot(3,3,6)
xlabel('Zc')
ylabel('Zp')
xlim([-0.5 2])
ylim([-0.4 0.4])
subplot(3,3,7)
xlabel('Zc')
xlim([-2 6])
title('DSI>0.5')
subplot(3,3,8)
xlabel('Zp')
xlim([-2 6])
title('')
subplot(3,3,9)
xlabel('Zc')
ylabel('Zp')
xlim([-0.5 2])
ylim([-0.4 0.4])
print(fullfile(summaryDir, ['randDir_allArea_summary_' driver '.pdf']),'-dpdf', '-fillpage') 

figure(2)
subplot(2,2,1)
legend(leg_str{2,:},'location','southeast')
xlabel('f2overf1')
title('')
subplot(2,2,2)
xlabel('F1')
ylabel('F2')
xlim([0 0.2])
ylim([0 0.2])
subplot(2,2,3)
xlabel('Zc')
ylabel('F2overF1')
xlim([-0.5 2])
ylim([0 1.25])
subplot(2,2,4)
xlabel('Zp')
ylabel('F2overF1')
xlim([-0.4 0.4])
ylim([0 1.25])
print(fullfile(summaryDir, ['randDir_allArea_summary_F2F1_' driver '.pdf']),'-dpdf', '-fillpage')

figure(3)
print(fullfile(summaryDir, ['randDir_allArea_summary_SI&F2F1byOSI_' driver '.pdf']),'-dpdf', '-fillpage')

[p_Zc table_Zc stats_Zc] = anova1(Zc_all_all(resp_ind_all_all),area_ind(resp_ind_all_all),'off');
post_Zc = multcompare(stats_Zc,'display','off');
[p_Zp table_Zp stats_Zp] = anova1(Zp_all_all(resp_ind_all_all),area_ind(resp_ind_all_all),'off');
post_Zp = multcompare(stats_Zp,'display','off');
[p_ZcZp table_ZcZp stats_ZcZp] = anova1(ZcZp_all_all(resp_ind_all_all),area_ind(resp_ind_all_all),'off');
post_ZcZp = multcompare(stats_ZcZp,'display','off');
[p_f2f1 table_f2f1 stats_f2f1] = anova1(f2overf1_all_all(f1resp_ind_all_all),area_ind(f1resp_ind_all_all),'off');
post_f2f1 = multcompare(stats_f2f1,'display','off');
figure(4) 
groups = mat2cell(post_Zc(:,1:2),[ones(1,size(post_Zc,1))],[2]);
subplot(2,2,1)
ind = find(post_Zc(:,end)<0.05);
xlim([0 narea+1])
ylim([-0.5 3])
sigstar(groups(ind),post_Zc(ind,end),1)
set(gca,'XTick',1:narea,'XTickLabel',area_list)
ylabel('Zc')
subplot(2,2,2)
ind = find(post_Zp(:,end)<0.05);
xlim([0 narea+1])
ylim([-0.5 0.5])
sigstar(groups(ind),post_Zp(ind,end),1)
set(gca,'XTick',1:narea,'XTickLabel',area_list)
ylabel('Zp')
subplot(2,2,3)
ind = find(post_ZcZp(:,end)<0.05);
xlim([0 narea+1])
ylim([-0.5 3])
sigstar(groups(ind),post_ZcZp(ind,end),1)
set(gca,'XTick',1:narea,'XTickLabel',area_list)
ylabel('Zc-Zp')
subplot(2,2,4)
ind = find(post_f2f1(:,end)<0.05);
xlim([0 narea+1])
ylim([0 1])
sigstar(groups(ind),post_f2f1(ind,end),1)
set(gca,'XTick',1:narea,'XTickLabel',area_list)
ylabel('F2/F1')
print(fullfile(summaryDir, ['randDir_allArea_summary_ZcZpF2F1_withStats_' driver '.pdf']),'-dpdf', '-fillpage')

figure(5)
print(fullfile(summaryDir, ['randDir_allArea_summary_ZcZp_scatters_' driver '.pdf']),'-dpdf', '-fillpage')

% [p_Zc_pca table_Zc_pca stats_Zc_pca] = anova1(Zc_pca_all,pca_area,'off');
% post_Zc_pca = multcompare(stats_Zc_pca,'display','off');
% [p_Zp_pca table_Zp_pca stats_Zp_pca] = anova1(Zp_pca_all,pca_area,'off');
% post_Zp_pca = multcompare(stats_Zp_pca,'display','off');

% figure(6)
% subplot(2,1,1)
% ind = find(post_Zc_pca(:,end)<0.05);
% xlim([0 narea+1])
% ylim([0 3.5])
% sigstar(groups(ind),post_Zc_pca(ind,end),1)
% set(gca,'XTick',1:narea,'XTickLabel',area_list)
% ylabel('PCA- Zc')
% subplot(2,1,2)
% ind = find(post_Zp_pca(:,end)<0.05);
% xlim([0 narea+1])
% ylim([0 3.5])
% sigstar(groups(ind),post_Zp_pca(ind,end),1)
% set(gca,'XTick',1:narea,'XTickLabel',area_list)
% ylabel('PCA- Zp')
% print(fullfile(summaryDir, ['randDir_allArea_summary_PCA_ZcZp.pdf']),'-dpdf', '-fillpage') 





