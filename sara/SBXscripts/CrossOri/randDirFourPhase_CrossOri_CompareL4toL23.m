%%% compare L4 to L2/3

close all; clear all; clc;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';
driver = strvcat('SLC', 'Scn'); 
area = 'all_areas';
area_list = strvcat('V1', 'V1');
narea = length(area_list);
nCells = [];


figure;
for iA = 1:2
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))];
    

    subplot(3,3,1)
        cdfplot(amp_all(resp_ind))
        hold on
        xlabel('fit amplitude')
        set(gca,'TickDir','out'); box off; axis square; grid off
    subplot(3,3,2)
        cdfplot(b_all(resp_ind))
        hold on
        xlabel('fit baseline')
        set(gca,'TickDir','out'); box off; axis square; grid off
    subplot(3,3,3)
        scatter(-amp_all(resp_ind),b_all(resp_ind),3,"filled")
        hold on
        ylabel('fit baseline')
        xlabel('-fit amp')
        set(gca,'TickDir','out'); box off; axis square; grid off 
    subplot(3,3,4)
        meanZc = mean(Zc_all,1);
        cdfplot(meanZc(resp_ind))
        hold on
        xlabel('mean Zc per cell')
        set(gca,'TickDir','out'); box off; axis square; grid off
    subplot(3,3,5)
        meanZp = mean(Zp_all,1);
        cdfplot(meanZp(resp_ind))
        hold on
        xlabel('mean Zp')
        set(gca,'TickDir','out'); box off; axis square; grid off
    subplot(3,3,6)
        maxZc = max(Zc_all,[],1);
        cdfplot(maxZc(resp_ind))
        hold on
        xlabel('max Zc')
        set(gca,'TickDir','out'); box off; axis square; grid off
    subplot(3,3,7)
        maxZp = max(Zp_all,[],1);
        cdfplot(maxZp(resp_ind))
        hold on
        xlabel('max Zp')
        set(gca,'TickDir','out'); box off; axis square; grid off

        movegui('center')
end
%print(fullfile(outDir, [svName '_VarianceByFitAmp.pdf']),'-dpdf', '-fillpage') 
