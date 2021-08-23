clear all
close all
clc
ds_list = {'randDirRandPhase'; 'randDirRandPhase'};
nCon = 2;
Cons = [0.5 0.125];
driver = {'SLC'};
area_list = {'V1'};
narea = length(area_list);

rc = behavConstsAV;
frame_rate = 15;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

leg_str = cell(4,nCon);
Zc_all_all = [];
Zp_all_all = [];
ZcZp_all_all = [];
resp_ind_all_all = [];
totCells = 0;
Con_ind = [];
for iCon = 1:nCon
    fprintf(['Con = ' num2str(Cons(iCon)) '\n'])
    ds = ds_list{iCon};
    con_str = num2str(Cons(iCon));
    con_str = ['pt' con_str(3:end)];
    summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', ['RandDirSummary']);
    load(fullfile(summaryDir,[ds '_Summary_' cell2mat(area_list) '_' cell2mat(driver) '_Con' con_str '.mat']))
    ConSummary(iCon).name = Cons(iCon);
    ConSummary(iCon).mice = unique(mouse_list,'rows');
    ConSummary(iCon).nmice = size(unique(mouse_list,'rows'),1);
    ConSummary(iCon).nexp = size(mouse_list,1);
    ConSummary(iCon).totCells = length(stim_OSI_all);
    ConSummary(iCon).respCells = length(resp_ind_all);
    figure(1)
    ind = resp_ind_all;
    leg_str{1,iCon} = ['Con = ' num2str(Cons(iCon)) '- ' num2str(length(ind)) ' cells'];
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
    leg_str{3,iCon} = ['Con = ' num2str(Cons(iCon)) '- ' num2str(length(ind_dsi)) ' cells'];
    subplot(3,3,7)
    cdfplot(Zc_all(ind_dsi))
    hold on
    subplot(3,3,8)
    cdfplot(Zp_all(ind_dsi))
    hold on
    subplot(3,3,9)
    errorbar(mean(Zc_all(ind_dsi),2),mean(Zp_all(ind_dsi),2),std(Zc_all(ind_dsi),[],2)./sqrt(length(ind_dsi)),std(Zc_all(ind_dsi),[],2)./sqrt(length(ind_dsi)),std(Zp_all(ind_dsi),[],2)./sqrt(length(ind_dsi)),std(Zp_all(ind_dsi),[],2)./sqrt(length(ind_dsi)),'o')
    hold on
    
    figure(3)
    subplot(2,2,1)
    cdfplot(stim_SI_all(ind))
    hold on
    subplot(2,2,2)
    cdfplot(plaid_SI_all(ind))
    hold on 
    subplot(2,2,3)
    cdfplot(Zc_all(ind))
    hold on
    subplot(2,2,4)
    cdfplot(Zp_all(ind))
    hold on
    
    figure(4)
    subplot(3,nCon,iCon)
    ind_h = intersect(resp_ind_all,find(stim_SI_all>0.5));
    ind_l = intersect(resp_ind_all,find(stim_SI_all<0.5));
    cdfplot(plaid_SI_all(ind_l))
    hold on
    cdfplot(plaid_SI_all(ind_h))
    xlim([-1 1])
    legend(['SI<0.5- n =' num2str(length(ind_l))],['SI>0.5- n =' num2str(length(ind_h))],'location','southeast')
    xlabel('Masking index')
    ylabel('Fraction of cells')
    title(['Con = ' num2str(Cons(iCon))])
    subplot(3,nCon,iCon+nCon)
    ind_h = intersect(resp_ind_all,find(stim_OSI_all>0.5));
    ind_l = intersect(resp_ind_all,find(stim_OSI_all<0.5));
    cdfplot(plaid_SI_all(ind_l))
    hold on
    cdfplot(plaid_SI_all(ind_h))
    xlim([-1 1])
    title('')
    legend(['OSI<0.5- n =' num2str(length(ind_l))],['OSI>0.5- n =' num2str(length(ind_h))],'location','southeast')
    xlabel('Masking index')
    ylabel('Fraction of cells')
    subplot(3,1,3)
    [n edges bin] = histcounts(stim_SI_all,[0:0.2:1]);
    temp = zeros(length(n),2,2);
    for i = 1:length(n)
        ind = intersect(resp_ind_all, find(bin == i));
        temp(i,1,1) = mean(stim_SI_all(ind));
        temp(i,1,2) = std(stim_SI_all(ind))./sqrt(length(ind));
        temp(i,2,1) = mean(plaid_SI_all(ind));
        temp(i,2,2) = std(plaid_SI_all(ind))./sqrt(length(ind));
    end
    errorbar(temp(:,1,1), temp(:,2,1), temp(:,2,2), temp(:,2,2), temp(:,1,2), temp(:,1,2), '-o')
    hold on
    xlabel('Stimulus selectivity')
    ylabel('Masking Index')
    
    
    Zc_all_all = [Zc_all_all Zc_all];
    Zp_all_all = [Zp_all_all Zp_all];
    ZcZp_all = Zc_all-Zp_all;
    ZcZp_all_all = [ZcZp_all_all ZcZp_all];
    resp_ind_all_all = [resp_ind_all_all resp_ind_all+totCells];
    Con_ind = [Con_ind iCon.*ones(size(Zc_all))];
    totCells = totCells+size(Zc_all,2);
    
    figure(6)
    subplot(2,2,1)
    errorbar(iCon, mean(Zc_all(resp_ind_all),2),std(Zc_all(resp_ind_all),[],2)./sqrt(length(resp_ind_all)),'ok')
    hold on
    subplot(2,2,2)
    errorbar(iCon, mean(Zp_all(resp_ind_all),2),std(Zp_all(resp_ind_all),[],2)./sqrt(length(resp_ind_all)),'ok')
    hold on
    subplot(2,2,3)
    errorbar(iCon, mean(ZcZp_all(resp_ind_all),2),std(ZcZp_all(resp_ind_all),[],2)./sqrt(length(resp_ind_all)),'ok')
    hold on
    
    figure(7)
    subplot(2,nCon,iCon)
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
    title(['Con = ' num2str(Cons(iCon))])
    axis square
    plotZcZpBorders
    subplot(2,nCon,iCon+nCon)
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
    title([Cons(iCon) '- DSI>0.5'])
    axis square
    plotZcZpBorders
    
    figure(8)
    subplot(2,2,iCon)
    mdl = fitlm(avg_resp_plaid_align_all(resp_ind_all,9),mean(avg_resp_dir_align_all(resp_ind_all,[9 13]),2));
    scatter(avg_resp_plaid_align_all(resp_ind_all,9),mean(avg_resp_dir_align_all(resp_ind_all,[9 13]),2));
    hold on
    [n edges bin] = histcounts(avg_resp_plaid_align_all(resp_ind_all,9),[0:0.1:1.2]);
    for i = 1:length(n)
        ind = find(bin==i);
        x = avg_resp_plaid_align_all(resp_ind_all(ind),9);
        y = mean(avg_resp_dir_align_all(resp_ind_all(ind),[9 13]),2);
        errorbar(mean(x),mean(y),std(y,[],1)./sqrt(length(ind)),std(y,[],1)./sqrt(length(ind)),std(x,[],1)./sqrt(length(ind)),std(x,[],1)./sqrt(length(ind)),'or')
    end
    xlabel('Plaid (dF/F)')
    ylabel('Avg Test & Mask')
    title([num2str(Cons(iCon)) '- ' num2str(chop(mdl.Coefficients.Estimate(2),2))])
    xlim([0 1.2])
    ylim([0 1.2])
    refline(1)
    subplot(2,2,iCon+2)
    ind = find(sum(h_resp_all(:,1,:),3));
    mdl = fitlm(avg_resp_dir_all(ind,1,2,1),mean(avg_resp_dir_all(ind,[1 5],1,1),2));
    scatter(avg_resp_dir_all(ind,1,2,1),mean(avg_resp_dir_all(ind,[1 5],1,1),2));
    hold on
    [n edges bin] = histcounts(avg_resp_dir_all(ind,1,2,1),[0:0.1:1.2]);
    for i = 1:length(n)
        temp_ind = find(bin==i);
        x = avg_resp_dir_all(ind(temp_ind),1,2,1);
        y = mean(avg_resp_dir_all(ind(temp_ind),[1 5],1,1),2);
        errorbar(mean(x),mean(y),std(y,[],1)./sqrt(length(ind)),std(y,[],1)./sqrt(length(ind)),std(x,[],1)./sqrt(length(ind)),std(x,[],1)./sqrt(length(ind)),'or')
    end
    xlabel('Plaid (dF/F)')
    ylabel('Avg Test & Mask')
    title([num2str(Cons(iCon)) '- ' num2str(chop(mdl.Coefficients.Estimate(2),2))])
    xlim([0 1.2])
    ylim([0 1.2])
    refline(1)
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
xlim([-1 3])
ylim([-1 3])
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
xlim([-1 3])
ylim([-1 3])
sgtitle([cell2mat(area_list) ' ' cell2mat(driver)])
print(fullfile(summaryDir, ['randDir_Con_summary_' cell2mat(area_list) '_' cell2mat(driver) '.pdf']),'-dpdf', '-fillpage') 

figure(3)
subplot(2,2,1)
legend(leg_str{1,:})
xlabel('Stim Selectivity Index')
title('')
subplot(2,2,2)
xlabel('Masking Index')
title('')
subplot(2,2,3)
xlabel('Zc')
title('')
subplot(2,2,4)
xlabel('Zp')
title('')
sgtitle([cell2mat(area_list) ' ' cell2mat(driver)])
print(fullfile(summaryDir, ['randDir_Con_summary_SelectivityComp' cell2mat(area_list) '_' cell2mat(driver) '.pdf']),'-dpdf', '-fillpage') 

figure(4)
legend(leg_str{1,:})
sgtitle([cell2mat(area_list) ' ' cell2mat(driver)])
print(fullfile(summaryDir, ['randDir_Con_summary_MIbySI_' cell2mat(area_list) '_' cell2mat(driver) '.pdf']),'-dpdf', '-fillpage')

[p_Zc table_Zc stats_Zc] = anova1(Zc_all_all(resp_ind_all_all),Con_ind(resp_ind_all_all),'off');
post_Zc = multcompare(stats_Zc,'display','off');
[p_Zp table_Zp stats_Zp] = anova1(Zp_all_all(resp_ind_all_all),Con_ind(resp_ind_all_all),'off');
post_Zp = multcompare(stats_Zp,'display','off');
[p_ZcZp table_ZcZp stats_ZcZp] = anova1(ZcZp_all_all(resp_ind_all_all),Con_ind(resp_ind_all_all),'off');
post_ZcZp = multcompare(stats_ZcZp,'display','off');

figure(6) 
groups = mat2cell(post_Zc(:,1:2),[ones(1,size(post_Zc,1))],[2]);
subplot(2,2,1)
ind = find(post_Zc(:,end)<0.05);
xlim([0 nCon+1])
ylim([-1 3])
sigstar(groups(ind),post_Zc(ind,end),1)
set(gca,'XTick',1:nCon,'XTickLabel',Cons)
ylabel('Zc')
xlabel('Con')
subplot(2,2,2)
ind = find(post_Zp(:,end)<0.05);
xlim([0 nCon+1])
ylim([-1 3])
sigstar(groups(ind),post_Zp(ind,end),1)
set(gca,'XTick',1:nCon,'XTickLabel',Cons)
ylabel('Zp')
xlabel('Con')
subplot(2,2,3)
ind = find(post_ZcZp(:,end)<0.05);
xlim([0 nCon+1])
ylim([-0.5 3])
sigstar(groups(ind),post_ZcZp(ind,end),1)
set(gca,'XTick',1:nCon,'XTickLabel',Cons)
ylabel('Zc-Zp')
xlabel('Con')
sgtitle([cell2mat(area_list) ' ' cell2mat(driver)])
print(fullfile(summaryDir, ['randDir_Con_summary_ZcZp_withStats_' cell2mat(area_list) '_' cell2mat(driver) '.pdf']),'-dpdf', '-fillpage')

figure(7)
sgtitle([cell2mat(area_list) ' ' cell2mat(driver)])
print(fullfile(summaryDir, ['randDir_Con_summary_ZcZp_scatters_' cell2mat(area_list) '_' cell2mat(driver) '.pdf']),'-dpdf', '-fillpage')

figure(8)
sgtitle([cell2mat(area_list) ' ' cell2mat(driver)])
print(fullfile(summaryDir, ['randDir_Con_summary_DarioFig_' cell2mat(area_list) '_' cell2mat(driver) '.pdf']),'-dpdf', '-fillpage')



