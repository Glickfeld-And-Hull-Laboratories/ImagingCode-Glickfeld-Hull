clear all
close all
clc
ds_list = {'randDirRandPhase'; 'randDirRandPhase'};
svName = 'randPhase';
nCon = 2;
Cons = [0.5 0.125];
driver = {'SLC'};
area_list = {'V1'};
narea = length(area_list);

rc = behavConstsAV;
frame_rate = 15;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

leg_str = cell(4,nCon);
b_all = [];
amp_all = [];
resp_ind_all_all = [];
totCells = 0;
Con_ind = [];
for iCon = 1:nCon
    fprintf(['Con = ' num2str(Cons(iCon)) '\n'])
    ds = ds_list{iCon};
    summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', ['RandPhaseSummary']);
    con_str = num2str(Cons(iCon));
    con_str = ['pt' con_str(3:end)];
    load(fullfile(summaryDir,[ds '_Summary_' cell2mat(area_list) '_' cell2mat(driver) '_Con' con_str '.mat']))
    ConSummary(iCon).name = Cons(iCon);
    ConSummary(iCon).mice = unique(mouse_list,'rows');
    ConSummary(iCon).nmice = size(unique(mouse_list,'rows'),1);
    ConSummary(iCon).nexp = size(mouse_list,1);
    ConSummary(iCon).totCells = length(b_all_all);
    ConSummary(iCon).respCells = length(resp_ind_all);
    figure(1)
    ind = resp_ind_all;
    leg_str{1,iCon} = ['Con = ' num2str(Cons(iCon)) '- ' num2str(length(ind)) ' cells'];
    subplot(2,2,1)
    cdfplot(testPI_all(ind))
    xlabel('Selectivity index')
    ylabel('Fraction of cells')
    title('')
    hold on
    subplot(2,2,2)
    cdfplot(plaidSI_all(ind))
    xlabel('Masking index')
    ylabel('Fraction of cells')
    title('')
    hold on
    subplot(2,2,3)
    cdfplot(b_all_all(ind))
    xlabel('Sine baseline')
    ylabel('Fraction of cells')
    title('')
    hold on
    subplot(2,2,4)
    cdfplot(amp_all_all(ind)-amp_shuf_all(ind))
    xlabel('Sine amplitude (-shuf)')
    ylabel('Fraction of cells')
    title('')
    hold on
    
    
    figure(2)
    subplot(2,nCon,iCon)
    ind_h = intersect(resp_ind_all,find(testPI_all>0.5));
    ind_l = intersect(resp_ind_all,find(testPI_all<0.5));
    cdfplot(b_all_all(ind_l))
    hold on
    cdfplot(b_all_all(ind_h))
    xlim([-1 1])
    legend(['SI<0.5- n =' num2str(length(ind_l))],['SI>0.5- n =' num2str(length(ind_h))],'location','southeast')
    xlabel('Sine baseline')
    ylabel('Fraction of cells')
    title(['Con = ' num2str(Cons(iCon))])
    subplot(2,nCon,iCon+nCon)
    cdfplot(amp_all_all(ind_l)-amp_shuf_all(ind_l))
    hold on
    cdfplot(amp_all_all(ind_h)-amp_shuf_all(ind_h))
    xlim([-1 1])
    title('')
    xlabel('Sine amplitude')
    ylabel('Fraction of cells')
    
    figure(3)
    subplot(2,nCon,iCon)
    ind_h = intersect(resp_ind_all,find(plaidSI_all>0));
    ind_l = intersect(resp_ind_all,find(plaidSI_all<0));
    cdfplot(b_all_all(ind_l))
    hold on
    cdfplot(b_all_all(ind_h))
    xlim([-1 1])
    legend(['MI<0- n =' num2str(length(ind_l))],['MI>0- n =' num2str(length(ind_h))],'location','southeast')
    xlabel('Sine baseline')
    ylabel('Fraction of cells')
    title(['Con = ' num2str(Cons(iCon))])
    subplot(2,nCon,iCon+nCon)
    cdfplot(amp_all_all(ind_l)-amp_shuf_all(ind_l))
    hold on
    cdfplot(amp_all_all(ind_h)-amp_shuf_all(ind_h))
    xlim([-1 1])
    title('')
    xlabel('Sine ampliutude')
    ylabel('Fraction of cells')
end
figure(1)
subplot(2,2,1)
legend(leg_str{1,:})
suptitle([cell2mat(area_list) ' ' cell2mat(driver)])
print(fullfile(summaryDir, [svName '_Con_summary_' cell2mat(area_list) '_' cell2mat(driver) '.pdf']),'-dpdf', '-fillpage') 

figure(2)
suptitle([cell2mat(area_list) ' ' cell2mat(driver)])
print(fullfile(summaryDir, [svName '_Con_summary_SelectivityComp' cell2mat(area_list) '_' cell2mat(driver) '.pdf']),'-dpdf', '-fillpage') 

figure(3)
suptitle([cell2mat(area_list) ' ' cell2mat(driver)])
print(fullfile(summaryDir, [svName '_Con_summary_MaskingComp' cell2mat(area_list) '_' cell2mat(driver) '.pdf']),'-dpdf', '-fillpage') 

% [p_Zc table_Zc stats_Zc] = anova1(Zc_all_all(resp_ind_all_all),Con_ind(resp_ind_all_all),'off');
% post_Zc = multcompare(stats_Zc,'display','off');
% [p_Zp table_Zp stats_Zp] = anova1(Zp_all_all(resp_ind_all_all),Con_ind(resp_ind_all_all),'off');
% post_Zp = multcompare(stats_Zp,'display','off');
% [p_ZcZp table_ZcZp stats_ZcZp] = anova1(ZcZp_all_all(resp_ind_all_all),Con_ind(resp_ind_all_all),'off');
% post_ZcZp = multcompare(stats_ZcZp,'display','off');



