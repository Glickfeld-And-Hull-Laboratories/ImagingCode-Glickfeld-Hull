clear all
close all
clc
doRedChannel = 0;
ds = 'CrossOriRandDirRandPhase_ExptList';
driver_list = {'SLC';'SOM';'PV';'SCN'};
ndriver = length(driver_list);
area_list = {'V1'};
narea = length(area_list);
eval(ds)
frame_rate = 15;
nexp = size(expt,2);
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
leg_str = cell(3,ndriver);
Zc_all_all = [];
Zp_all_all = [];
ZcZp_all_all = [];
resp_ind_dirall = [];
resp_ind_phaseall = [];
amp_all_all = [];
amp_shuf_all_all = [];
b_all_all = [];
stim_SI_all_all = [];
plaid_SI_all_all = [];
phase_SI_all_all = [];
phase_MI_all_all = [];
mouse_ind_all = [];
totCells = 0;
driver_ind_all = [];
totm = 1;
for iD = 1:ndriver
    fprintf([driver_list{iD} '\n'])
    driver = driver_list{iD};
    load(fullfile(summaryDir,['randDirRandPhase_Summary_' cell2mat(area_list) '_' driver '.mat']))
    driverSummary(iD).name = driver_list{iD};
    driverSummary(iD).mice = unique(mouse_list,'rows');
    driverSummary(iD).nmice = size(unique(mouse_list,'rows'),1);
    driverSummary(iD).nexp = size(mouse_list,1);
    driverSummary(iD).totCells = length(stim_OSI_all);
    driverSummary(iD).respCells_dir = length(resp_ind_all_dir);
    driverSummary(iD).respCells_phase = length(resp_ind_all_phase);
    figure(1)
    ind_dir = resp_ind_all_dir(find(~isnan(plaid_SI_all(resp_ind_all_dir))));
    ind_phase = resp_ind_all_phase(find(~isnan(b_all(resp_ind_all_phase))));
    fprintf([num2str(length(ind_dir)) ' cells\n'])
    leg_str{1,iD} = [driver_list{iD} '- ' num2str(length(ind_dir)) ' cells'];
    leg_str{2,iD} = [driver_list{iD} '- ' num2str(length(ind_phase)) ' cells'];
    subplot(3,3,1)
    cdfplot(stim_OSI_all(ind_dir))
    hold on
    subplot(3,3,2)
    cdfplot(stim_DSI_all(ind_dir))
    hold on
    subplot(3,3,3)
    cdfplot(k_all(ind_dir))
    hold on
    subplot(3,3,4)
    cdfplot(Zc_all(ind_dir))
    hold on
    subplot(3,3,5)
    cdfplot(Zp_all(ind_dir))
    hold on
    subplot(3,3,6)
    errorbar(mean(Zc_all(ind_dir),2),mean(Zp_all(ind_dir),2),std(Zc_all(ind_dir),[],2)./sqrt(length(ind_dir)),std(Zc_all(ind_dir),[],2)./sqrt(length(ind_dir)),std(Zp_all(ind_dir),[],2)./sqrt(length(ind_dir)),std(Zp_all(ind_dir),[],2)./sqrt(length(ind_dir)),'o')
    hold on
    ind_dsi = intersect(ind_dir,find(stim_DSI_all>0.5));
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
    subplot(3,2,1)
    cdfplot(stim_SI_all(ind_dir))
    hold on
    subplot(3,2,2)
    cdfplot(plaid_SI_all(ind_dir))
    hold on
    subplot(3,2,3)
    cdfplot(phase_SI_all(ind_phase))
    hold on
    subplot(3,2,4)
    cdfplot(phase_MI_all(ind_phase))
    hold on
    subplot(3,2,5)
    cdfplot(b_all(ind_phase))
    hold on
    subplot(3,2,6)
    cdfplot(amp_all(ind_phase)-amp_shuf_all(ind_phase))
    hold on

    figure(3)
    subplot(2,ndriver,iD)
    ind_h = intersect(ind_dir,find(stim_SI_all>0.5));
    ind_l = intersect(ind_dir,find(stim_SI_all<0.5));
    cdfplot(plaid_SI_all(ind_l))
    hold on
    cdfplot(plaid_SI_all(ind_h))
    xlim([-1 1])
    legend(['SI<0.5- n =' num2str(length(ind_l))],['SI>0.5- n =' num2str(length(ind_h))],'location','southeast')
    xlabel('Masking index')
    ylabel('Fraction of cells')
    title(driver_list{iD})
    subplot(2,ndriver,iD+ndriver)
    ind_h = intersect(ind_dir,find(stim_OSI_all>0.5));
    ind_l = intersect(ind_dir,find(stim_OSI_all<0.5));
    cdfplot(plaid_SI_all(ind_l))
    hold on
    cdfplot(plaid_SI_all(ind_h))
    xlim([-1 1])
    legend(['OSI<0.5- n =' num2str(length(ind_l))],['OSI>0.5- n =' num2str(length(ind_h))],'location','southeast')
    xlabel('Masking index')
    ylabel('Fraction of cells')
    
    figure(4)
    subplot(3,1,1)
    [n edges bin] = histcounts(stim_SI_all,[0:0.2:1]);
    temp = zeros(length(n),2,2);
    for i = 1:length(n)
        ind = intersect(ind_dir, find(bin == i));
        temp(i,1,1) = nanmean(stim_SI_all(ind));
        temp(i,1,2) = nanstd(stim_SI_all(ind))./sqrt(length(ind));
        temp(i,2,1) = nanmean(plaid_SI_all(ind));
        temp(i,2,2) = nanstd(plaid_SI_all(ind))./sqrt(length(ind));
    end
    errorbar(temp(:,1,1), temp(:,2,1), temp(:,2,2), temp(:,2,2), temp(:,1,2), temp(:,1,2), '-o')
    hold on
    xlabel('Stimulus selectivity')
    ylabel('Masking Index')
    
    subplot(3,1,2)
    [n edges bin] = histcounts(phase_SI_all,[0:0.2:1]);
    temp = zeros(length(n),2,2);
    for i = 1:length(n)
        ind = intersect(ind_phase, find(bin == i));
        temp(i,1,1) = nanmean(phase_SI_all(ind));
        temp(i,1,2) = nanstd(phase_SI_all(ind))./sqrt(length(ind));
        temp(i,2,1) = nanmean(b_all(ind));
        temp(i,2,2) = nanstd(b_all(ind))./sqrt(length(ind));
    end
    errorbar(temp(:,1,1), temp(:,2,1), temp(:,2,2), temp(:,2,2), temp(:,1,2), temp(:,1,2), '-o')
    hold on
    xlabel('Stimulus selectivity')
    ylabel('Baseline')
    
    subplot(3,1,3)
    [n edges bin] = histcounts(phase_SI_all,[0:0.2:1]);
    temp = zeros(length(n),2,2);
    for i = 1:length(n)
        ind = intersect(ind_phase, find(bin == i));
        temp(i,1,1) = nanmean(phase_SI_all(ind));
        temp(i,1,2) = nanstd(phase_SI_all(ind))./sqrt(length(ind));
        temp(i,2,1) = nanmean(amp_all(ind));
        temp(i,2,2) = nanstd(amp_all(ind))./sqrt(length(ind));
    end
    errorbar(temp(:,1,1), temp(:,2,1), temp(:,2,2), temp(:,2,2), temp(:,1,2), temp(:,1,2), '-o')
    hold on
    xlabel('Stimulus selectivity')
    ylabel('Sine Amplitude')
    
    
    stim_SI_all_all = [stim_SI_all_all stim_SI_all];
    plaid_SI_all_all = [plaid_SI_all_all plaid_SI_all];
    Zc_all_all = [Zc_all_all Zc_all];
    Zp_all_all = [Zp_all_all Zp_all];
    ZcZp_all = Zc_all-Zp_all;
    ZcZp_all_all = [ZcZp_all_all ZcZp_all];
    b_all_all = [b_all_all b_all];
    amp_all_all = [amp_all_all amp_all];
    amp_shuf_all_all = [amp_shuf_all_all amp_shuf_all];
    phase_SI_all_all = [phase_SI_all_all phase_SI_all];
    resp_ind_dirall = [resp_ind_dirall ind_dir+totCells];
    resp_ind_phaseall = [resp_ind_phaseall ind_phase+totCells];
    
    driver_ind_temp = mat2cell(repmat(driver,size(Zc_all')),ones(size(Zc_all,2),1));
    driver_ind_all = [driver_ind_all;  driver_ind_temp];
    mouse_ind_temp = mat2cell(mouse_ind,ones(size(Zc_all,2),1));
    mouse_ind_all = [mouse_ind_all;  mouse_ind_temp];
    totCells = totCells+size(Zc_all,2);
    
    figure(5)
    subplot(2,2,1)
    errorbar(iD, mean(Zc_all(ind_dir),2),std(Zc_all(ind_dir),[],2)./sqrt(length(ind_dir)),'ok')
    hold on
    subplot(2,2,2)
    errorbar(iD, mean(Zp_all(ind_dir),2),std(Zp_all(ind_dir),[],2)./sqrt(length(ind_dir)),'ok')
    hold on
    subplot(2,2,3)
    errorbar(iD, mean(ZcZp_all(ind_dir),2),std(ZcZp_all(ind_dir),[],2)./sqrt(length(ind_dir)),'ok')
    hold on
    
    figure(9)
    mice = unique(mouse_ind_temp);
    nm = length(mice);
    for i = 1:nm
        ind_use = intersect(ind_dir,find(strcmp(mouse_ind_temp,mice(i))));
        subplot(2,1,1)
        swarmchart(totm.*ones(size(Zc_all(ind_use))),Zc_all(ind_use)',[],defaultPlotColors(iD))
        hold on
        subplot(2,1,2)
        swarmchart(totm.*ones(size(Zp_all(ind_use))),Zp_all(ind_use)',[],defaultPlotColors(iD))
        hold on
        totm = totm+1;
    end
    
    
    figure(6)
    subplot(2,ndriver,iD)
    scatter(Zc_all(ind_dir),Zp_all(ind_dir),'ok')
    hold on
    Zc_use = intersect(ind_dir,intersect(find(Zc_all>1.28),find(Zc_all-Zp_all>1.28)));
    scatter(Zc_all(Zc_use),Zp_all(Zc_use),'ob')
    Zp_use = intersect(ind_dir,intersect(find(Zp_all>1.28),find(Zp_all-Zc_all>1.28)));
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
    
    figure(7)
    driver_col = defaultPlotColors();
    subplot(ndriver,ndriver, 1+((iD-1)*ndriver))
    scatter(b_all(ind_phase),plaid_SI_all(ind_phase),[],'MarkerEdgeColor',driver_col(iD,:))
    title(['r = ' num2str(chop(triu2vec(corrcoef(b_all(ind_phase),plaid_SI_all(ind_phase))),2))])
    xlabel('Baseline')
    ylabel('Masking Index')
    xlim([-1 1])
    ylim([-1 1])
    subplot(ndriver,ndriver, 2+((iD-1)*ndriver))
    scatter(amp_all(ind_phase),plaid_SI_all(ind_phase),[],'MarkerEdgeColor',driver_col(iD,:))
    title(['r = ' num2str(chop(triu2vec(corrcoef(amp_all(ind_phase),plaid_SI_all(ind_phase))),2))])
    xlabel('Amplitude')
    ylabel('Masking Index')
    xlim([0 1])
    ylim([-1 1])
    subplot(ndriver,ndriver, 3+((iD-1)*ndriver))
    scatter(b_all(ind_phase),amp_all(ind_phase),[],'MarkerEdgeColor',driver_col(iD,:))
    title(['r = ' num2str(chop(triu2vec(corrcoef(b_all(ind_phase),amp_all(ind_phase))),2))])
    xlabel('Baseline')
    ylabel('Amplitude')
    xlim([-1 1])
    ylim([0 1])
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
print(fullfile(summaryDir, ['randDir_allDriver_summary_' cell2mat(area_list) '.pdf']),'-dpdf', '-fillpage') 


figure(2)
subplot(3,2,1)
legend(leg_str{1,:},'location','northwest')
xlabel('Stim Selectivity Index')
title('')
subplot(3,2,2)
xlabel('Masking Index')
title('')
subplot(3,2,3)
legend(leg_str{2,:},'location','southeast')
xlabel('Stim Selectivity Index')
subplot(3,2,4)
xlabel('Masking Index')
title('')
subplot(3,2,5)
xlabel('Baseline')
title('')
subplot(3,2,6)
xlabel('Amp-Shuf')
title('')
print(fullfile(summaryDir, ['randDir_allDriver_summary_SelectivityComp' cell2mat(area_list) '.pdf']),'-dpdf', '-fillpage') 

figure(3)
print(fullfile(summaryDir, ['randDir_allDriver_summary_MIbySI_cdfs_' cell2mat(area_list) '.pdf']),'-dpdf', '-fillpage')

[drivers,~,driverIx] = unique(driver_ind_all);
figure(4)
subplot(3,1,1)
legend(driver_list)
ylim([-0.6 0.1])
[h_MI,atab_MI,ctab_MI,stats_MI] = aoctool(stim_SI_all_all(resp_ind_dirall),plaid_SI_all_all(resp_ind_dirall),driverIx(resp_ind_dirall),[],'','','','off');
title(['Main effect SI: p = ' num2str(chop(atab_MI{3,end},2)) '; Cell type: p = ' num2str(chop(atab_MI{2,end},2)) '; Int: p = ' num2str(chop(atab_MI{4,end},2))])
subplot(3,1,2)
ylim([-0.6 0.1])
[h_B,atab_B,ctab_B,stats_B] = aoctool(phase_SI_all_all(resp_ind_phaseall),b_all_all(resp_ind_phaseall),driverIx(resp_ind_phaseall),[],'','','','off');
title(['Main effect SI: p = ' num2str(chop(atab_B{3,end},2)) '; Cell type: p = ' num2str(chop(atab_B{2,end},2)) '; Int: p = ' num2str(chop(atab_B{4,end},2))])
subplot(3,1,3)
ylim([0 0.3])
[h_A,atab_A,ctab_A,stats_A] = aoctool(phase_SI_all_all(resp_ind_phaseall),amp_all_all(resp_ind_phaseall),driverIx(resp_ind_phaseall),[],'','','','off');
title(['Main effect SI: p = ' num2str(chop(atab_A{3,end},2)) '; Cell type: p = ' num2str(chop(atab_A{2,end},2)) '; Int: p = ' num2str(chop(atab_A{4,end},2))])
print(fullfile(summaryDir, ['randDir_allDriver_summary_MIbySI_bins_' cell2mat(area_list) '.pdf']),'-dpdf', '-fillpage')


[p_Zc table_Zc stats_Zc] = anova1(Zc_all_all(resp_ind_dirall),driverIx(resp_ind_dirall),'off');
post_Zc = multcompare(stats_Zc,'display','off');
[p_Zp table_Zp stats_Zp] = anova1(Zp_all_all(resp_ind_dirall),driverIx(resp_ind_dirall),'off');
post_Zp = multcompare(stats_Zp,'display','off');
[p_ZcZp table_ZcZp stats_ZcZp] = anova1(ZcZp_all_all(resp_ind_dirall),driverIx(resp_ind_dirall),'off');
post_ZcZp = multcompare(stats_ZcZp,'display','off');

[mice,~,miceIx] = unique(mouse_ind_all);
tbl_Zc = table(Zc_all_all(resp_ind_dirall)',driver_ind_all(resp_ind_dirall),mouse_ind_all(resp_ind_dirall),'VariableNames',{'Zc','Driver','Mouse'});
lme_Zc = fitlme(tbl_Zc,'Zc~Driver+(1|Mouse)');
lme_Zc_nofe = fitlme(tbl_Zc,'Zc~Driver');
emm_Zc = emmeans(glme_Zc,{'Driver'});
outZc = contrasts_wald(glme_Zc,emm_Zc,[1 0 0 -1]);
outZc = contrasts_wald(glme_Zc,emm_Zc,[1 0 -1]);
outZc = contrasts_wald(glme_Zc,emm_Zc,[0 1 -1]);

tbl_Zp = table(Zp_all_all(resp_ind_dirall)',driver_ind_all(resp_ind_dirall)',miceIx(resp_ind_dirall),'VariableNames',{'Zp','Driver','Mouse'});
lme_Zp = fitlme(tbl_Zp,'Zp~Driver+(1|Mouse)');
tbl_ZcZp = table(Zc_all_all(resp_ind_dirall)'-Zp_all_all(resp_ind_dirall)',driver_ind_all(resp_ind_dirall)',miceIx(resp_ind_dirall),'VariableNames',{'ZcZp','Driver','Mouse'});
lme_ZcZp = fitlme(tbl_ZcZp,'ZcZp~Driver+(1|Mouse)');

figure(5) 
groups = mat2cell(post_Zc(:,1:2),[ones(1,size(post_Zc,1))],[2]);
subplot(2,2,1)
ind = find(post_Zc(:,end)<0.05);
xlim([0 ndriver+1])
ylim([-1 3])
sigstar(groups(ind),post_Zc(ind,end),1)
set(gca,'XTick',1:ndriver,'XTickLabel',driver_list)
ylabel('Zc')
subplot(2,2,2)
ind = find(post_Zp(:,end)<0.05);
xlim([0 ndriver+1])
ylim([-1 3])
sigstar(groups(ind),post_Zp(ind,end),1)
set(gca,'XTick',1:ndriver,'XTickLabel',driver_list)
ylabel('Zp')
subplot(2,2,3)
ind = find(post_ZcZp(:,end)<0.05);
xlim([0 ndriver+1])
ylim([-1 3])
sigstar(groups(ind),post_ZcZp(ind,end),1)
set(gca,'XTick',1:ndriver,'XTickLabel',driver_list)
ylabel('Zc-Zp')
print(fullfile(summaryDir, ['randDir_allDriver_summary_ZcZp_withStats_' cell2mat(area_list) '.pdf']),'-dpdf', '-fillpage')

figure(6)
print(fullfile(summaryDir, ['randDir_allDriver_summary_ZcZp_scatters_' cell2mat(area_list) '.pdf']),'-dpdf', '-fillpage')

figure(7)

figure(9)
subplot(2,1,1)
ylim([-8 8])
set(gca,'Xtick',1:length(unique(mouse_ind_all)),'XtickLabel',unique(mouse_ind_all))
subplot(2,1,2)
ylim([-8 8])
set(gca,'Xtick',1:length(unique(mouse_ind_all)),'XtickLabel',unique(mouse_ind_all))


