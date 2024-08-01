%% Analyzing PV and SOM cells
close all; clear all; clc;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';
driver = strvcat('syn'); 
area_list = strvcat('V1');
narea = length(area_list);


figure;
for iA = 1 %1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(ind_all))];
    
    ind_red = red_cells_all;
    ind_nonRed = setdiff(ind_all, red_cells_all);
    
    subplot(3,2,iA)
        cdfplot(DSI_all)
        hold on
        cdfplot(DSI_all(ind_nonRed))
        cdfplot(DSI_all(ind_red))
        subtitle('DSI')
        set(gca,'TickDir','out'); box off; grid off
        
    subplot(3,2,iA+2)
        cdfplot(k1_all)
        hold on
        cdfplot(k1_all(ind_nonRed))
        cdfplot(k1_all(ind_red))
        subtitle('k1')
        set(gca,'TickDir','out'); box off; grid off
        
    subplot(3,2,iA+4)
        cdfplot(b_all)
        hold on
        cdfplot(b_all(ind_nonRed))
        cdfplot(b_all(ind_red))
        subtitle('b')
        set(gca,'TickDir','out'); box off; grid off
end



% print(fullfile(outDir, [svName '_SummaryByExperiment.pdf']),'-dpdf', '-fillpage') 

