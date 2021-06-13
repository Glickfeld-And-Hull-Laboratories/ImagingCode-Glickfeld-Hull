clear all
close all
clc
doRedChannel = 0;
ds = 'CrossOriRandDir_ExptList';
driver_list = {'SLC';'SOM';'PV'};
ndriver = length(driver_list);
area_list = {'V1'};
sf = 'pt05';
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
f2overf1_all_all = [];
resp_ind_all_all = [];
f1resp_ind_all_all = [];
totCells = 0;
driver_ind = [];
for iD = 1:ndriver
    fprintf([driver_list{iD} '\n'])
    driver = driver_list{iD};
    load(fullfile(summaryDir,['randDir_Summary_' cell2mat(area_list) '_' driver '_SF' sf '.mat']))
    driverSummary(iD).name = driver_list{iD};
    driverSummary(iD).mice = unique(mouse_list,'rows');
    driverSummary(iD).nmice = size(unique(mouse_list,'rows'),1);
    driverSummary(iD).nexp = size(mouse_list,1);
    driverSummary(iD).totCells = length(stim_OSI_all);
    driverSummary(iD).respCells = length(resp_ind_all);
    figure(1)
    ind = resp_ind_all;
    leg_str{1,iD} = [driver_list{iD} '- ' num2str(length(ind)) ' cells'];
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
    leg_str{3,iD} = [driver_list{iD} '- ' num2str(length(ind_dsi)) ' cells'];
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
    leg_str{2,iD} = [driver_list{iD} '- ' num2str(length(ind)) ' cells'];
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
    subplot(3,ndriver,iD)
    ind_h = intersect(resp_ind_all,find(stim_SI_all>0.5));
    ind_l = intersect(resp_ind_all,find(stim_SI_all<0.5));
    cdfplot(plaid_SI_all(ind_l))
    hold on
    cdfplot(plaid_SI_all(ind_h))
    xlim([-1 1])
    legend(['SI<0.5- n =' num2str(length(ind_l))],['SI>0.5- n =' num2str(length(ind_h))],'location','southeast')
    xlabel('Masking index')
    ylabel('Fraction of cells')
    title(driver_list{iD})
    subplot(3,ndriver,iD+ndriver)
    ind_h = intersect(resp_ind_all,find(stim_OSI_all>0.5));
    ind_l = intersect(resp_ind_all,find(stim_OSI_all<0.5));
    cdfplot(plaid_SI_all(ind_l))
    hold on
    cdfplot(plaid_SI_all(ind_h))
    xlim([-1 1])
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
    
    figure(4)
    if ~isempty(f1_ind)
    subplot(2,ndriver,iD)
    cdfplot(f2overf1_all(ind_l))
    hold on
    cdfplot(f2overf1_all(ind_h))
    xlabel('F2/F1')
    ylabel('Fraction of cells')
    legend(['OSI<0.5- n =' num2str(length(ind_l))],['OSI>0.5- n =' num2str(length(ind_h))],'location','southeast')
    xlim([0 5])
    ind_h = intersect(resp_ind_all,find(plaid_SI_all>0));
    ind_l = intersect(resp_ind_all,find(plaid_SI_all<0));
    subplot(2,ndriver,iD+(ndriver))
    cdfplot(f2overf1_all(ind_l))
    hold on
    cdfplot(f2overf1_all(ind_h))
    xlabel('F2/F1')
    ylabel('Fraction of cells')
    xlim([0 5])
    legend(['MI<0- n =' num2str(length(ind_l))],['MI>0- n =' num2str(length(ind_h))],'location','southeast')
    end
    
    Zc_all_all = [Zc_all_all Zc_all];
    Zp_all_all = [Zp_all_all Zp_all];
    ZcZp_all = Zc_all-Zp_all;
    ZcZp_all_all = [ZcZp_all_all ZcZp_all];
    resp_ind_all_all = [resp_ind_all_all resp_ind_all+totCells];
    driver_ind = [driver_ind iD.*ones(size(Zc_all))];
    if ~isempty(f1_ind)
        f2overf1_all_all = [f2overf1_all_all f2overf1_all];
        f1resp_ind_all_all = [f1resp_ind_all_all f1_ind+totCells];
    else
        f2overf1_all_all = [f2overf1_all_all nan(size(Zc_all))];
        f1resp_ind_all_all = [f1resp_ind_all_all f1_ind+totCells];
    end
    totCells = totCells+size(Zc_all,2);
    
    figure(5)
    subplot(2,2,1)
    errorbar(iD, mean(Zc_all(resp_ind_all),2),std(Zc_all(resp_ind_all),[],2)./sqrt(length(resp_ind_all)),'ok')
    hold on
    subplot(2,2,2)
    errorbar(iD, mean(Zp_all(resp_ind_all),2),std(Zp_all(resp_ind_all),[],2)./sqrt(length(resp_ind_all)),'ok')
    hold on
    subplot(2,2,3)
    errorbar(iD, mean(ZcZp_all(resp_ind_all),2),std(ZcZp_all(resp_ind_all),[],2)./sqrt(length(resp_ind_all)),'ok')
    hold on
    if ~isempty(f1_ind)
        subplot(2,2,4)
        errorbar(iD, mean(f2overf1_all(f1_ind),2),std(f2overf1_all(f1_ind),[],2)./sqrt(length(f1_ind)),'ok')
        hold on
    end
    
    figure(6)
    subplot(2,ndriver,iD)
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
    title(driver_list{iD})
    axis square
    plotZcZpBorders
    subplot(2,ndriver,iD+ndriver)
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
    title([driver_list{iD} '- DSI>0.5'])
    axis square
    plotZcZpBorders
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
ylim([-0.5 2])
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
ylim([-0.5 2])
print(fullfile(summaryDir, ['randDir_allDriver_summary_' cell2mat(area_list) '_SF' sf '.pdf']),'-dpdf', '-fillpage') 

if ~isempty(f1resp_ind_all_all)
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
print(fullfile(summaryDir, ['randDir_allDriver_summary_F2F1_' cell2mat(area_list) '_SF' sf '.pdf']),'-dpdf', '-fillpage')
end

figure(3)
print(fullfile(summaryDir, ['randDir_allDriver_summary_MIbySI_' cell2mat(area_list) '_SF' sf '.pdf']),'-dpdf', '-fillpage')

figure(4)
print(fullfile(summaryDir, ['randDir_allDriver_summary_F2F1byOSI_' cell2mat(area_list) '_SF' sf '.pdf']),'-dpdf', '-fillpage')

[p_Zc table_Zc stats_Zc] = anova1(Zc_all_all(resp_ind_all_all),driver_ind(resp_ind_all_all),'off');
post_Zc = multcompare(stats_Zc,'display','off');
[p_Zp table_Zp stats_Zp] = anova1(Zp_all_all(resp_ind_all_all),driver_ind(resp_ind_all_all),'off');
post_Zp = multcompare(stats_Zp,'display','off');
[p_ZcZp table_ZcZp stats_ZcZp] = anova1(ZcZp_all_all(resp_ind_all_all),driver_ind(resp_ind_all_all),'off');
post_ZcZp = multcompare(stats_ZcZp,'display','off');
if ~isempty(f1resp_ind_all_all)
    [p_f2f1 table_f2f1 stats_f2f1] = anova1(f2overf1_all_all(f1resp_ind_all_all),driver_ind(f1resp_ind_all_all),'off');
    post_f2f1 = multcompare(stats_f2f1,'display','off');
end
figure(5) 
groups = mat2cell(post_Zc(:,1:2),[ones(1,size(post_Zc,1))],[2]);
subplot(2,2,1)
ind = find(post_Zc(:,end)<0.05);
xlim([0 ndriver+1])
ylim([-0 2])
sigstar(groups(ind),post_Zc(ind,end),1)
set(gca,'XTick',1:ndriver,'XTickLabel',driver_list)
ylabel('Zc')
subplot(2,2,2)
ind = find(post_Zp(:,end)<0.05);
xlim([0 ndriver+1])
ylim([0 2])
sigstar(groups(ind),post_Zp(ind,end),1)
set(gca,'XTick',1:ndriver,'XTickLabel',driver_list)
ylabel('Zp')
subplot(2,2,3)
ind = find(post_ZcZp(:,end)<0.05);
xlim([0 ndriver+1])
ylim([0 2])
sigstar(groups(ind),post_ZcZp(ind,end),1)
set(gca,'XTick',1:ndriver,'XTickLabel',driver_list)
ylabel('Zc-Zp')
if ~isempty(f1resp_ind_all_all)
    subplot(2,2,4)
    ind = find(post_f2f1(:,end)<0.05);
    xlim([0 ndriver+1])
    ylim([0 1])
    sigstar(groups(ind),post_f2f1(ind,end),1)
    set(gca,'XTick',1:ndriver,'XTickLabel',driver_list)
    ylabel('F2/F1')
end
print(fullfile(summaryDir, ['randDir_allDriver_summary_ZcZpF2F1_withStats_' cell2mat(area_list) '_SF' sf '.pdf']),'-dpdf', '-fillpage')

figure(6)
print(fullfile(summaryDir, ['randDir_allDriver_summary_ZcZp_scatters_' cell2mat(area_list) '_SF' sf '.pdf']),'-dpdf', '-fillpage')




