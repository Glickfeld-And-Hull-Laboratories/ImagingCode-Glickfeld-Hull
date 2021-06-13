clc; clear all; close all; close all hidden
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'RandDirSummary');
area_list = ['V1'; 'LM'; 'AL'; 'PM'; 'RL'];

narea =length(area_list);
stimplaid_corrs_all_all = [];
stimplaid_dists_all_all = [];
stimplaid_angs_all_all = [];
Zc_all_all = [];
Zp_all_all = [];
area_ind = [];
leg_str = cell(1,narea);
    
for iarea = 1:narea
    area = area_list(iarea,:);
    load(fullfile(summaryDir,['randDir_PCA_Summary_' area '.mat']))
    stimplaid_corrs_all_all = [stimplaid_corrs_all_all; stimplaid_corrs_all'];
    stimplaid_dists_all_all = [stimplaid_dists_all_all; stimplaid_dists_all'];
    stimplaid_angs_all_all = [stimplaid_angs_all_all; stimplaid_angs_all'];
    Zc_all_all = [Zc_all_all; Zc_all(find(~isnan(Zc_all)))'];
    Zp_all_all = [Zp_all_all; Zp_all(find(~isnan(Zp_all)))'];
    area_ind = [area_ind; iarea*ones(size(Zc_all(find(~isnan(Zc_all))),2),1)];
    
    leg_str{iarea} = [area '- n = ' num2str(size(Zc_all(find(~isnan(Zc_all))),2))];
    
    figure(1)
    subplot(2,2,1) 
    errorbar([1:3],nanmean(stimplaid_corrs_all,2),nanstd(stimplaid_corrs_all,[],2)./sqrt(sum(~isnan(stimplaid_corrs_all(1,:)))),'-o')
    hold on
    subplot(2,2,2) 
    errorbar([1:3],nanmean(stimplaid_dists_all,2),nanstd(stimplaid_dists_all,[],2)./sqrt(sum(~isnan(stimplaid_dists_all(1,:)))),'-o')    
    hold on 
    subplot(2,2,3) 
    errorbar([1:3],nanmean(stimplaid_angs_all,2),nanstd(stimplaid_angs_all,[],2)./sqrt(sum(~isnan(stimplaid_angs_all(1,:)))),'-o')    
    hold on
    subplot(2,2,4)
    errorbar([1 2], [nanmean(Zc_all,2)' nanmean(Zp_all,2)'], [nanstd(Zc_all,[],2)'./sqrt(size(Zc_all,2)) nanstd(Zp_all,[],2)'./sqrt(size(Zc_all,2))],'-o')
    hold on
    
    figure(2)
    subplot(2,2,1)
    scatter(stimplaid_angs_all(1,:),Zc_all(find(~isnan(Zc_all))))
    hold on
    subplot(2,2,2)
    scatter(stimplaid_corrs_all(1,:),Zc_all(find(~isnan(Zc_all))))
    hold on
    subplot(2,2,3)
    scatter(stimplaid_angs_all(2,:),Zp_all(find(~isnan(Zp_all))))
    hold on
    subplot(2,2,4)
    scatter(stimplaid_corrs_all(2,:),Zp_all(find(~isnan(Zp_all))))
    hold on
end
figure(1)
subplot(2,2,1)
set(gca,'Xtick', 1:3, 'XTickLabels',{'Component','Vector','Opposite'})
xlim([0 4])
ylim([0 1])
ylabel('Plaid-Stim Correlation')
[p_corr table_corr] = anovan(stimplaid_corrs_all_all(:),{repmat(area_ind, [3 1]),[ones(size(area_ind)); 2.*ones(size(area_ind)); 3.*ones(size(area_ind))]});
title(['Area- p = ' num2str(chop(p_corr(1),2)) '; Cond- p = ' num2str(chop(p_corr(2),2))])
legend(leg_str)
subplot(2,2,2)
set(gca,'Xtick', 1:3, 'XTickLabels',{'Component','Vector','Opposite'})
xlim([0 4])
ylim([0 8])
ylabel('Plaid-Stim Neural distance')
[p_dist table_dist] = anovan(stimplaid_dists_all_all(:),{repmat(area_ind, [3 1]),[ones(size(area_ind)); 2.*ones(size(area_ind)); 3.*ones(size(area_ind))]});
title(['Area- p = ' num2str(chop(p_dist(1),2)) '; Cond- p = ' num2str(chop(p_dist(2),2))])
subplot(2,2,3)
set(gca,'Xtick', 1:3, 'XTickLabels',{'Component','Vector','Opposite'})
xlim([0 4])
ylim([0 90])
ylabel('Plaid-Stim Neural angle')
[p_ang table_ang] = anovan(stimplaid_angs_all_all(:),{repmat(area_ind, [3 1]),[ones(size(area_ind)); 2.*ones(size(area_ind)); 3.*ones(size(area_ind))]});
title(['Area- p = ' num2str(chop(p_ang(1),2)) '; Cond- p = ' num2str(chop(p_ang(2),2))])
subplot(2,2,4)
set(gca,'Xtick',1:2,'XTickLabel',{'Component','Pattern'})
xlim([0 3])
ylim([0 3])
ylabel('Correlation with plaid')
[p_Z table_Z] = anovan([Zc_all_all; Zp_all_all],{repmat(area_ind, [2 1]),[ones(size(area_ind)); 2.*ones(size(area_ind))]});
title(['Area- p = ' num2str(chop(p_Z(1),2)) '; Cond- p = ' num2str(chop(p_Z(2),2))])
suptitle('PCA summary- 16 Direction- SF = 0.1')    
%print(fullfile(summaryDir,['randDir_PCA_allArea_Summary.pdf']),'-dpdf', '-bestfit')


figure(2)
subplot(2,2,1)
xlabel('Component angle')
ylabel('Component Z')
xlim([0 90])
ylim([0 3])
axis square
subplot(2,2,2)
xlabel('Component correlation')
ylabel('Component Z')
xlim([0 1])
ylim([0 3])
axis square
subplot(2,2,3)
xlabel('Vector angle')
ylabel('Pattern Z')
xlim([0 90])
ylim([0 2])
axis square
subplot(2,2,4)
xlabel('Vector correlation')
ylabel('Pattern Z')
xlim([0 1])
ylim([0 2])
axis square
suptitle('PCA summary- 16 Direction- SF = 0.1')

