%% Figures comparing areas with SI fit (re-load .mat files generated in first section)
close all; clear all; clc;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
svName = 'randPhase';
driver = strvcat('SLC_all','SLC_all','SLC_all'); 
area = 'all_areas';
area_list = strvcat('V1','LM','AL');
narea = length(area_list);
nCells = [];


figure;
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    ind = resptest_ind_all;
    leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))];

    subplot(3,2,1)
        Rsq_temp = R_square_all_all;
        Rsq_temp(find(Rsq_temp<0)) = 0;
        cdfplot(Rsq_temp(ind,:))
        n = sum(~isnan(Rsq_temp(ind,:)));
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
%         subtitle('Rsquared')
        xlabel('Rsquared')
        legend(leg_str, 'location', 'southeast')

    subplot(3,2,2)
        cdfplot(amp_all_all(ind,:))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
%         subtitle('Phase modulation amplitude')
        xlabel('Phase modulation amplitude')
        legend(leg_str, 'location', 'southeast')

    subplot(3,2,3)
        cdfplot(b_all_all(ind))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        ylabel('# of cells')
%         subtitle('Phase modulation baseline')
        xlabel('Phase modulation baseline')
        legend(area_list, 'location', 'southeast')  
        
    subplot(3,2,4)
        cdfplot(amp_shuf_all(ind))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        ylabel('# of cells')
        subtitle('SHUFFLED')
        xlabel('Phase modulation amplitude (shuf)')
        xlim([0 0.8])
        legend(area_list, 'location', 'southeast')  
end
    sgtitle('RandPhase spatial invariance measure of V1, AL, LM (ind=responsive to test)')
    print(fullfile(outDir, [svName '_' area '_RandPhaseFigure_ForLindsey.pdf']),'-dpdf', '-fillpage') 



%     
% driver = strvcat('SLC_008','SLC_010','SLC_058','SLC_014','SLC_042','SLC_043','SLC_052','SLC_053','SLC_054','SLC_044','SLC_045','SLC_046','SLC_047','SLC_048','SLC_007','SLC_029','SLC_031','SLC_049','SLC_050','SLC_055','SLC_056','SLC_057');
% area_list = strvcat('V1','V1','V1','V1','V1','V1','V1','V1','V1','LM','LM','LM','LM','LM','AL','AL','AL','AL','AL','AL','AL','AL');
% narea = length(area_list);
% 
% for iA = 1:narea
%     load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
%     ind = resptest_ind_all;
%     leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(ind))];  
% 
%     if area_list(iA,:) == 'V1'
%         ll = '-';
%         c = [0 0.4470 0.7410];
%     elseif area_list(iA,:) == 'LM'
%         ll = '--';
%         c = [0.8500 0.3250 0.0980];
%     else 
%         ll = ':';
%         c = [0.9290 0.6940 0.1250];
%     end
% 
%     subplot(3,2,5)
%         [f,x_cdf] = ecdf(amp_all_all(ind));
%         plot(x_cdf,f,'Color',c,'LineWidth',1)
%         hold on
%         ylabel('# of cells')
%         xlabel('Phase modulation amplitude')
% %         legend(area_list, 'location', 'southeast')          
% end
% 
% %     print(fullfile(outDir, [svName '_' area '_Compare.pdf']),'-dpdf', '-fillpage') 
