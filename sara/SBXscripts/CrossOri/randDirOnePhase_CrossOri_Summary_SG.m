close all; clear all; clc;
doRedChannel = 1;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary', 'summaries');

doPlot = 1;
ds = ['CrossOriRandDirFourPhase_ExptList_SG'];
svName = 'randPhase';
eval(ds)
driver = 'SCN';
img_area = {'LM';'L2/3'}; %LM
inj_area = 'LM';
img_layer = 'L2/3';

max_dist = 5;

rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);


b_all = [];
amp_all = [];
rsq_all = [];
sse_all = [];

ind_all = [];
ZpZcPWdist_all = [];
plaid_corr_all = [];

mouse_list = [];
totCells = zeros(nexp,1);

%V1 L2/3 - 5 24 46 47 79
%V1 L4 - 62 66
%LM L2/3 - 38 49 50 86
% AL L2/3 - 76 82 83 84

start=1;
for iexp = [38 49 50 86]
    mouse = expt(iexp).mouse;
    mouse_list = strvcat(mouse_list, mouse);
    date = expt(iexp).date;
    if isfield(expt,'copFolder') 
        ImgFolder = expt(iexp).copFolder;
    else
        ImgFolder = expt(iexp).coFolder;
    end
    nrun = length(ImgFolder);
    run_str = 'runs-002';
        
    fprintf([mouse ' ' date '\n'])
        
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_bootstrapFits.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_ZpZc_pairwiseDist.mat']))
    
    nCells = size(Zp,2);
    totCells(iexp,:) = nCells;
    
    b_all = [b_all; mean(boot_base,2)];
    amp_all = [amp_all; mean(boot_amp,2)];
    rsq_all = [rsq_all; mean(boot_rsq,2)];
    sse_all = [sse_all; mean(boot_sse,2)];
    

    ind_all = [ind_all; ind+sum(totCells(1:iexp-1,:),1)];

    plaid_corr_all = [plaid_corr_all; plaid_corr'];
    ZpZcPWdist_all = [ZpZcPWdist_all; mean(ZpZcPWdist(:,ind),1)'];
    
    start=start+1;
end
    save(fullfile(summaryDir,[svName '_OnePhaseSummary_' inj_area '_' driver '.mat']), 'ind_all', 'ZpZcPWdist_all', 'plaid_corr_all','b_all', 'amp_all', 'rsq_all', 'sse_all', 'mouse_list')



%% plotting to compare spatial invariance across areas

close all; clear all; clc;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';
driver = strvcat('SLC','SLC','SLC'); 
area = 'all_areas';
area_list = strvcat('V1', 'LM', 'AL');
narea = length(area_list);
nCells = [];

figure;
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_OnePhaseSummary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))

    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(ind_all))];

    subplot(3,3,1)
        histogram(b_all(ind_all,:))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('fit baseline')
        legend(leg_str)
    subplot(3,3,2)
        histogram(amp_all(ind_all,:))
        hold on
        xlim([0 5])
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('fit amplitude')
    subplot(3,3,3)
        cdfplot(amp_all(ind_all,:))
        hold on
        xlim([0 5])
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('fit amplitude')
    subplot(3,3,4)
        histogram(mean(ZpZcPWdist_all,1))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('PW dist of Zp Zc, avg per cell')
    subplot(3,3,5)
        cdfplot(mean(ZpZcPWdist_all,1))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('PW dist of Zp Zc, avg per cell')
    subplot(3,3,6)
        histogram(plaid_corr_all(ind_all))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('plaid-plaid correlation')
     subplot(3,3,7)
        cdfplot(plaid_corr_all(ind_all))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('plaid-plaid correlation')

    
    if iA==1
         V1amp = amp_all(ind_all);
    elseif iA==2
        LMamp = amp_all(ind_all);
        [h,p] = ttest2(LMamp,V1amp,'Tail','right','Alpha',0.01);
        fprintf(['V1&LM h=' num2str(h) ' p=' num2str(p) '\n'])
    elseif iA ==3
        ALamp = amp_all(ind_all);
        [h,p] = ttest2(ALamp,V1amp,'Tail','right','Alpha',0.01);
        fprintf(['V1&AL h=' num2str(h) ' p=' num2str(p) '\n'])
        [h,p] = ttest2(ALamp,LMamp,'Tail','right','Alpha',0.01);
        fprintf(['LM&AL h=' num2str(h) ' p=' num2str(p) '\n'])
    end

end

print(fullfile(outDir, [svName '_OnePhaseComparison_AllAreas.pdf']),'-dpdf', '-fillpage') 



%% 
close all; clear all; clc;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';
driver = strvcat('SLC'); %Scn
area = 'all_areas';
area_list = strvcat('AL');
img_layer = 'L23';
narea = length(area_list);
nCells = [];


figure(1);
for iA = 1
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:)  '_' driver(iA,:) '.mat'])))
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))];


    subplot(3,3,1)
        histogram(b_all(resp_ind,:),-5:0.25:5)
        hold on
        xline(median(b_all(resp_ind,:)),'--')
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('fit baseline')
        legend(leg_str{iA})
    subplot(3,3,2)
        histogram(amp_all(resp_ind,:),0:0.2:6)
        hold on
        xline(median(amp_all(resp_ind,:)),'--')
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('fit amplitude')
    subplot(3,3,3)
        cdfplot(amp_all(resp_ind,:))
        hold on
        xlim([0 5])
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('fit amplitude')
    subplot(3,3,4)
        histogram(mean(ZpZcPWdist_all(:,resp_ind),1),0:0.2:5)
        hold on
        xline(median(mean(ZpZcPWdist_all(:,resp_ind,:))),'--')
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('PW dist of Zp Zc, avg per cell')
    subplot(3,3,5)
        cdfplot(mean(ZpZcPWdist_all(:,resp_ind),1))
        hold on
        xlim([0 5])
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('PW dist of Zp Zc, avg per cell')
    subplot(3,3,6)
        histogram(plaid_corr_all(resp_ind),-1:0.1:1.5)
        hold on
        xline(median(plaid_corr_all(resp_ind,:)),'--')
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('plaid-plaid correlation')
     subplot(3,3,7)
        cdfplot(plaid_corr_all(resp_ind))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('plaid-plaid correlation')

        fouramp = amp_all(resp_ind);
end


figure(1);
for iA = 1
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_OnePhaseSummary_' area_list(iA,:)  '_' driver(iA,:) '.mat'])))
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(ind_all))];
    

    fprintf(['mean ' num2str(mean(amp_all(ind_all))) '\n sem ' num2str(std(amp_all(ind_all))/length(amp_all(ind_all))) '\n'])


    subplot(3,3,1)
        histogram(b_all(ind_all,:),-5:0.25:5)
        hold on
        xline(median(b_all(ind_all,:)),'--')
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('fit baseline')
    subplot(3,3,2)
        histogram(amp_all(ind_all,:),0:0.2:6)
        hold on
        xline(median(amp_all(ind_all,:)),'--')
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('fit amplitude')
        legend(leg_str{iA})
    subplot(3,3,3)
        cdfplot(amp_all(ind_all,:))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('fit amplitude')
    subplot(3,3,4)
        histogram(ZpZcPWdist_all,0:0.2:5)
        hold on
        xline(median(ZpZcPWdist_all),'--')
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('PW dist of Zp Zc, avg per cell')
    subplot(3,3,5)
        cdfplot(ZpZcPWdist_all)
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('PW dist of Zp Zc, avg per cell')    
    subplot(3,3,6)
        histogram(plaid_corr_all(ind_all),-1:0.1:1.5)
        hold on
        xline(median(plaid_corr_all(ind_all,:)),'--')
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('plaid-plaid correlation')
    subplot(3,3,7)
        cdfplot(plaid_corr_all(ind_all))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('plaid-plaid correlation')

        oneamp = amp_all(ind_all);

end


[h,p] = ttest2(fouramp,oneamp,'Tail','right','Alpha',0.01);
fprintf([' h=' num2str(h) ' p=' num2str(p) '\n'])


figure(1);
print(fullfile(outDir, [svName '_OneVsFourPhaseSummary_' area_list(iA,:)  '_' img_layer '.pdf']),'-dpdf', '-fillpage') 

stop


figure(2);
for iA = 1
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' img_layer '_' driver(iA,:) '.mat'])))
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))];
    
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end

    subplot(2,2,1)
        scatter(amp_all(resp_ind),Rsq_all(resp_ind),'.')
        hold on
        subtitle('Four Phase')
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('Amplitude')
        ylabel('R-squared')
        ylim([-0.5 1])

    subplot(2,2,2)
        scatter(amp_all(resp_ind),sse_all_all(resp_ind),'.')
        hold on
        subtitle('Four Phase')
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('Amplitude')
        ylabel('SSE')


        fouramp = amp_all(resp_ind);
end


figure(2);
for iA = 1
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_OnePhaseSummary_' area_list(iA,:) '_' img_layer '_' driver(iA,:) '.mat'])))
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(ind_all))];
    
    if exist('red_cells_all','var')
        ind_all = setdiff(ind_all, red_cells_all);
    end

    subplot(2,2,3)
        scatter(amp_all(ind_all),rsq_all(ind_all),'.','Color', [0.7 0.7 0.7])
        hold on
        subtitle('One Phase Control')
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('Amplitude')
        ylabel('R-squared')
        ylim([-0.5 1])

    subplot(2,2,4)
        scatter(amp_all(ind_all),sse_all(ind_all),'.','Color', [0.7 0.7 0.7])
        hold on
        subtitle('One Phase Control')
        set(gca,'TickDir','out'); box off; axis square; grid off
        ylim([0 25])
        xlim([0 5])
        xlabel('Amplitude')
        ylabel('SSE')
end

figure(2);
print(fullfile(outDir, [svName '_OneVsFourPhaseSummary_GoodnessOfFits_' area_list(iA,:)  '_' img_layer '.pdf']),'-dpdf', '-fillpage') 
