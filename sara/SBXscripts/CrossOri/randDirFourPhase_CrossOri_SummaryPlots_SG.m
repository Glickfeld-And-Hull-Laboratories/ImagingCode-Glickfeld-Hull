% 

%% summary CDFs, compare across areas
close all; clear all; clc;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';
driver = strvcat('SLC', 'SLC', 'SLC', 'Scn', 'SLC'); 
area = 'all_areas';
area_list = strvcat('V1', 'LM', 'AL', 'V1', 'PM');
narea = length(area_list);
nCells = [];


figure(1);
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))];
    
    fprintf(['mean ' num2str(mean(amp_all(resp_ind))) '\n sem ' num2str(std(amp_all(resp_ind))/length(amp_all(resp_ind))) '\n'])
    
    subplot(5,5,1)
        cdfplot(Rsq_all(resp_ind,:))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('fit Rsq')
        legend(leg_str, 'location', 'southeast'); set(gca,'TickDir','out'); box off; grid off
        
    subplot(5,5,2)
        cdfplot(sse_all_all(resp_ind,:))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('fit SSE')
    
    subplot(5,5,3)
        cdfplot(b_all(resp_ind,:))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('fit baseline')
    
    subplot(5,5,4)
        cdfplot(amp_all(resp_ind,:))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('fit amplitude')
        
    subplot(5,5,5)
        cdfplot(mean(Zp_all(:,resp_ind),1))
        hold on
        % xlim([-2 5])
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('Zp all')
        
    subplot(5,5,6)
        cdfplot(mean(PCI_all(:,resp_ind),1))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('PCI all')
    
    subplot(5,5,7)
        cdfplot(mean(Zc_all(:,resp_ind),1))
        hold on
        % xlim([-2 5])
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('Zc all')
    
    subplot(5,5,8)
        cdfplot(mean(ZpZcPWdist_all(:,resp_ind),1))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('PW dist of Zp Zc, avg per cell')
    
    subplot(5,5,9)
        cdfplot(plaid_corr_all(resp_ind))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('plaid-plaid correlation')

    subplot(5,5,10)
        cdfplot(max(Zc_all(:,resp_ind),[],1))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('Zc - max per cell')

    subplot(5,5,11)
        cdfplot(max(Zp_all(:,resp_ind),[],1))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('Zp - max per cell')

    subplot(5,5,12)
        cdfplot(k1_all(resp_ind))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('k1')

    subplot(5,5,13)
        cdfplot(max(PCI_all(:,resp_ind),[],1))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('PCI max')

    if iA==1
         V1amp = amp_all(resp_ind);
    elseif iA==2
        LMamp = amp_all(resp_ind);
        [h,p] = ttest2(LMamp,V1amp,'Tail','right','Alpha',0.01);
        fprintf(['V1&LM h=' num2str(h) ' p=' num2str(p) '\n'])
    elseif iA ==3
        ALamp = amp_all(resp_ind);
        [h,p] = ttest2(ALamp,V1amp,'Tail','right','Alpha',0.01);
        fprintf(['V1&AL h=' num2str(h) ' p=' num2str(p) '\n'])
        [h,p] = ttest2(ALamp,LMamp,'Tail','right','Alpha',0.01);
        fprintf(['LM&AL h=' num2str(h) ' p=' num2str(p) '\n'])
    end
end
print(fullfile(outDir, [svName '_Summary.pdf']),'-dpdf', '-fillpage') 


figure(2);
iAcount=narea+1;
iAcount2 = narea*2+1;
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    totCells = size(DSI_all,1);
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))];

    [popZc_std popZc_avg] = std(mean(Zc_all(:,resp_ind),1));
    [popZp_std popZp_avg] = std(mean(Zp_all(:,resp_ind),1));
    popZc_sem = popZc_std/(sqrt(length(resp_ind)));
    popZp_sem = popZp_std/(sqrt(length(resp_ind)));

    subplot(4,6,1)
        scatter(popZc_avg,popZp_avg,10,"filled")
        hold on
        ylim([0 0.5]); ylabel('avg Zp')
        xlim([0 2]); xlabel('avg Zc')
        errorbar(popZc_avg,popZp_avg,popZp_sem,'Color',[.7 .7 .7],"LineStyle","none")
        errorbar(popZc_avg,popZp_avg,popZc_sem,"horizontal",'Color',[.7 .7 .7],"LineStyle","none")
        set(gca,'TickDir','out'); box off; axis square; grid off
        subtitle([area_list(iA,:) ' ' driver(iA,:)])


    inclusion_crit(1) = size(sig_stim,1)/totCells;
    inclusion_crit(2) = size(sig_dir,2)/totCells;
    inclusion_crit(3) = size(find(DSI_all>0.5),1)/totCells;
    inclusion_crit(4) = length(resp_ind)/totCells;

    subplot(4,6,iA+2)
        bar(inclusion_crit);
        hold on
        ylim([0 1])
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('Inclusion Criteria (stim, dir, DSI)')
        subtitle([area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(totCells)])

    ind_patt(1) = length(intersect(intersect(find(Zp_all(1,:)>1.28),find(Zp_all(1,:)-Zc_all(1,:)>1.28)),resp_ind));
    ind_patt(2) = length(intersect(intersect(find(Zp_all(2,:)>1.28),find(Zp_all(2,:)-Zc_all(2,:)>1.28)),resp_ind));
    ind_patt(3) = length(intersect(intersect(find(Zp_all(3,:)>1.28),find(Zp_all(2,:)-Zc_all(2,:)>1.28)),resp_ind));
    ind_patt(4) = length(intersect(intersect(find(Zp_all(4,:)>1.28),find(Zp_all(4,:)-Zc_all(4,:)>1.28)),resp_ind));
    ind_comp(1) = length(intersect(intersect(find(Zc_all(1,:)>1.28),find(Zc_all(1,:)-Zp_all(1,:)>1.28)),resp_ind));
    ind_comp(2) = length(intersect(intersect(find(Zc_all(2,:)>1.28),find(Zc_all(2,:)-Zp_all(2,:)>1.28)),resp_ind));
    ind_comp(3) = length(intersect(intersect(find(Zc_all(3,:)>1.28),find(Zc_all(2,:)-Zp_all(2,:)>1.28)),resp_ind));
    ind_comp(4) = length(intersect(intersect(find(Zc_all(4,:)>1.28),find(Zc_all(4,:)-Zp_all(4,:)>1.28)),resp_ind));

    avg_classification(1) = mean(ind_patt(:))/length(resp_ind);
    avg_classification(2) = mean(ind_comp(:))/length(resp_ind);

    ind_patt1 = intersect(intersect(find(Zp_all(1,:)>1.28),find(Zp_all(1,:)-Zc_all(1,:)>1.28)),resp_ind);
    ind_patt2 = intersect(intersect(find(Zp_all(2,:)>1.28),find(Zp_all(2,:)-Zc_all(2,:)>1.28)),resp_ind);
    ind_patt3 = intersect(intersect(find(Zp_all(3,:)>1.28),find(Zp_all(2,:)-Zc_all(2,:)>1.28)),resp_ind);
    ind_patt4 = intersect(intersect(find(Zp_all(4,:)>1.28),find(Zp_all(4,:)-Zc_all(4,:)>1.28)),resp_ind);
    ind_comp1 = intersect(intersect(find(Zc_all(1,:)>1.28),find(Zc_all(1,:)-Zp_all(1,:)>1.28)),resp_ind);
    ind_comp2 = intersect(intersect(find(Zc_all(2,:)>1.28),find(Zc_all(2,:)-Zp_all(2,:)>1.28)),resp_ind);
    ind_comp3 = intersect(intersect(find(Zc_all(3,:)>1.28),find(Zc_all(2,:)-Zp_all(2,:)>1.28)),resp_ind);
    ind_comp4 = intersect(intersect(find(Zc_all(4,:)>1.28),find(Zc_all(4,:)-Zp_all(4,:)>1.28)),resp_ind);
    
    max_classification(1) = length(unique([ind_patt1;ind_patt2;ind_patt3;ind_patt4]))/length(resp_ind);
    max_classification(2) = length(unique([ind_comp1;ind_comp2;ind_comp3;ind_comp4]))/length(resp_ind);

    subplot(4,6,iAcount+2)
        bar(avg_classification);
        hold on
        ylim([0 0.6])
        yticks([0 0.2 0.4 0.6])
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('Average classification (PDS, CDS)')
        subtitle([area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))])    

    subplot(4,6,iAcount2+2)
        bar(max_classification);
        hold on
        ylim([0 1.])
        yticks([0 0.5 1])
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('Max classification (PDS, CDS)')
        subtitle([area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))])    



iAcount=iAcount+1;
iAcount2=iAcount2+1;
end

print(fullfile(outDir, [svName '_SummaryContinued.pdf']),'-dpdf', '-fillpage') 




figure(3);
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))];

    
    subplot(2,2,1)
        histogram(b_all(resp_ind,:))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('fit baseline')
    if iA==4    
        subplot(2,2,2)
            scatter(b_all(resp_ind,:),amp_all(resp_ind,:))
            hold on
            set(gca,'TickDir','out'); box off; axis square; grid off
            xlabel('fit baseline')
            ylabel('fit amp')
    end
end

figure(4);
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))];
    
    
    ind1 = intersect(find(Zp_all(1,:)>1.28),find(Zp_all(1,:)-Zc_all(1,:)>1.28));
    ind2 = intersect(find(Zp_all(2,:)>1.28),find(Zp_all(2,:)-Zc_all(2,:)>1.28));
    ind4 = intersect(find(Zp_all(4,:)>1.28),find(Zp_all(4,:)-Zc_all(4,:)>1.28));
    pattern1 = intersect(ind1,resp_ind);
    pattern2 = intersect(ind2,resp_ind);
    pattern4 = intersect(ind4,resp_ind);
    
    subplot(4,5,iA)
        scatter(Zc_all(1,resp_ind),Zp_all(1,resp_ind),15,[0.7 0.7 0.7])
        hold on
        scatter(Zc_all(1,pattern1),Zp_all(1,pattern1),15)
        xlabel('Zc'); ylabel('Zp'); xlim([-4 8]); ylim([-4 8]); plotZcZpBorders; set(gca,'TickDir','out'); axis square
        subtitle(leg_str(:,iA))
        
   subplot(4,5,iA+5)
        scatter(Zc_all(2,resp_ind),Zp_all(2,resp_ind),15,[0.7 0.7 0.7])
        hold on
        scatter(Zc_all(2,pattern1),Zp_all(2,pattern1),15)
        xlabel('Zc'); ylabel('Zp'); xlim([-4 8]); ylim([-4 8]); plotZcZpBorders; set(gca,'TickDir','out'); axis square
        
   subplot(4,5,iA+10)
        scatter(Zc_all(3,resp_ind),Zp_all(3,resp_ind),15,[0.7 0.7 0.7])
        hold on
        scatter(Zc_all(3,pattern1),Zp_all(3,pattern1),15)
        xlabel('Zc'); ylabel('Zp'); xlim([-4 8]); ylim([-4 8]); plotZcZpBorders; set(gca,'TickDir','out'); axis square
        
   subplot(4,5,iA+15)
        scatter(Zc_all(4,resp_ind),Zp_all(4,resp_ind),15,[0.7 0.7 0.7])
        hold on
        scatter(Zc_all(4,pattern1),Zp_all(4,pattern1),15)
        xlabel('Zc'); ylabel('Zp'); xlim([-4 8]); ylim([-4 8]); plotZcZpBorders; set(gca,'TickDir','out'); axis square
end
print(fullfile(outDir, [svName '_SummaryZpZc.pdf']),'-dpdf', '-fillpage') 


figure(5);
iAcount=narea+1;
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))];
    
%     Zp_avg = mean(Zp_all,1);
%     Zc_avg = mean(Zc_all,1);
%     pattind = (Zp_avg-Zc_avg);
    pattern_ind = (Zp_all-Zc_all);    
    pattind = mean(pattern_ind,1);
    pattpeak = max(pattern_ind,[],1);

    ind1 = intersect(find(Zp_all(1,:)>1.28),find(Zp_all(1,:)-Zc_all(1,:)>1.28));
    ind2 = intersect(find(Zp_all(2,:)>1.28),find(Zp_all(2,:)-Zc_all(2,:)>1.28));
    ind3 = intersect(find(Zp_all(3,:)>1.28),find(Zp_all(3,:)-Zc_all(3,:)>1.28));
    ind4 = intersect(find(Zp_all(4,:)>1.28),find(Zp_all(4,:)-Zc_all(4,:)>1.28));
    pattern1 = intersect(ind1,resp_ind);
    pattern2 = intersect(ind2,resp_ind);
    pattern3 = intersect(ind3,resp_ind);
    pattern4 = intersect(ind4,resp_ind);
    p = [pattern1; pattern2; pattern3; pattern4];
    [C,ia,ic] = unique(p);    
    a_counts = accumarray(ic,1);

    p1 = C(find(a_counts==1),:);
    p2 = C(find(a_counts==2),:);
    p3 = C(find(a_counts==3),:);
    p4 = C(find(a_counts==4),:);

    c1 = [0.9375    0.7813    0.7813];
    c2 = [0.9023    0.5742    0.5625];
    c3 = [0.8320    0.3672    0.3398];
    c4 = [0.7266    0.1094    0.1094];

    % Linear regression fit
        pCoef = polyfit(-amp_all(resp_ind),pattind(resp_ind),1);
        pFit = polyval(pCoef,-amp_all(resp_ind));
        fprintf(['   slope of fit = ' num2str(pCoef(1)) ' \n'])
    % Correlation
        [rho, pval] = corr(-amp_all(resp_ind),pattind(resp_ind)');
        fprintf(['   pval of correlation = ' num2str(pval) ' \n'])

    subplot(4,4,iA)
        scatter(-amp_all(resp_ind),pattind(resp_ind),10,'filled')
        hold on
        scatter(-amp_all(p1),pattind(p1),10,c1,'filled')
        scatter(-amp_all(p2),pattind(p2),10,c2,'filled')
        scatter(-amp_all(p3),pattind(p3),10,c3,'filled')
        scatter(-amp_all(p4),pattind(p4),10,c4,'filled')
        plot(-amp_all(resp_ind),pFit)
        ylabel('Mean pattern index (Zp-Zc)')
        xlabel('Spatial invariance (-amp)')
        ylim([-4 8])
        xlim([-6 0])
        set(gca,'TickDir','out'); axis square
        subtitle(leg_str(:,iA))
    
    subplot(4,4,iAcount)
        scatter(pattpeak(resp_ind),amp_all(resp_ind),20,'filled')
        hold on
        scatter(pattpeak(p1),amp_all(p1),20,c1,'filled')
        scatter(pattpeak(p2),amp_all(p2),20,c2,'filled')
        scatter(pattpeak(p3),amp_all(p3),20,c3,'filled')
        scatter(pattpeak(p4),amp_all(p4),20,c4,'filled')
        xlabel('Peak pattern index (Zp-Zc)')
        ylabel('Spatial variance (amp)')
        xlim([-4 8])
        set(gca,'TickDir','out'); axis square
        subtitle(leg_str(:,iA))
iAcount=iAcount+1;        
end
print(fullfile(outDir, [svName '_Summary_likeNicholas.pdf']),'-dpdf', '-fillpage') 




figure;
for iA = 1
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))];
    
%     Zp_avg = mean(Zp_all,1);
%     Zc_avg = mean(Zc_all,1);
%     pattind = (Zp_avg-Zc_avg);
    pattern_ind = (Zp_all-Zc_all);    
    pattind = mean(pattern_ind,1);
    pattpeak = max(pattern_ind,[],1);

    ind1 = intersect(find(Zp_all(1,:)>1.28),find(Zp_all(1,:)-Zc_all(1,:)>1.28));
    ind2 = intersect(find(Zp_all(2,:)>1.28),find(Zp_all(2,:)-Zc_all(2,:)>1.28));
    ind3 = intersect(find(Zp_all(3,:)>1.28),find(Zp_all(3,:)-Zc_all(3,:)>1.28));
    ind4 = intersect(find(Zp_all(4,:)>1.28),find(Zp_all(4,:)-Zc_all(4,:)>1.28));
    pattern1 = intersect(ind1,resp_ind);
    pattern2 = intersect(ind2,resp_ind);
    pattern3 = intersect(ind3,resp_ind);
    pattern4 = intersect(ind4,resp_ind);
    p = [pattern1; pattern2; pattern3; pattern4];
    [C,ia,ic] = unique(p);    
    a_counts = accumarray(ic,1);

    p1 = C(find(a_counts==1),:);
    p2 = C(find(a_counts==2),:);
    p3 = C(find(a_counts==3),:);
    p4 = C(find(a_counts==4),:);

    c1 = [0.9375    0.7813    0.7813];
    c2 = [0.9023    0.5742    0.5625];
    c3 = [0.8320    0.3672    0.3398];
    c4 = [0.7266    0.1094    0.1094];
    
    subplot(2,2,1)
        scatter(-amp_all(resp_ind),pattind(resp_ind),[],'filled')
        hold on
        scatter(-amp_all(p1),pattind(p1),[],c1,'filled')
        scatter(-amp_all(p2),pattind(p2),[],c2,'filled')
        scatter(-amp_all(p3),pattind(p3),[],c3,'filled')
        scatter(-amp_all(p4),pattind(p4),[],c4,'filled')
        ylabel('Mean pattern index (Zp-Zc)')
        xlabel('Spatial invariance (-amp)')
        ylim([-4 6])
        
        set(gca,'TickDir','out'); axis square
        subtitle(leg_str(:,iA))

    subplot(2,2,2)
        scatter(pattpeak(resp_ind),amp_all(resp_ind),[],'filled')
        hold on
        scatter(pattpeak(p1),amp_all(p1),[],c1,'filled')
        scatter(pattpeak(p2),amp_all(p2),[],c2,'filled')
        scatter(pattpeak(p3),amp_all(p3),[],c3,'filled')
        scatter(pattpeak(p4),amp_all(p4),[],c4,'filled')
        xlabel('Peak pattern index (Zp-Zc)')
        ylabel('Spatial variance (amp)')
        xlim([-4 6])
        set(gca,'TickDir','out'); axis square
        subtitle(leg_str(:,iA))
        
end
print(fullfile(outDir, [svName '_Summary_likeNicholas_justV1.pdf']),'-dpdf', '-fillpage') 



stop


figure;
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))];

    subplot(3,2,1)
        bm = mean(b_all(ind_all));
        bsem = std(b_all(ind_all))/length(b_all(ind_all));
        bar(iA,bm)
        errorbar(iA,bm,bsem,-bsem,'k')
        hold on
        ylim([-1.5 1.5])
        set(gca,'TickDir','out'); box off; axis square; grid off
        ylabel('avg PCI fit baseline')
    
    subplot(3,2,2)
        PCIm = mean(mean(PCI_all(:,ind_all),1));
        PCIsem = std(mean(PCI_all(:,ind_all),1))/size(PCI_all(:,ind_all),2);
        bar(iA,PCIm)
        hold on
        errorbar(iA,PCIm,PCIsem,-PCIsem,'k')
        ylim([-1.5 1.5])
        set(gca,'TickDir','out'); box off; axis square; grid off
        ylabel('avg PCI')

end
print(fullfile(outDir, [svName '_Summary2.pdf']),'-dpdf', '-fillpage') 


figure;
for iA = 1
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))];


    FWHM = rad2deg(2.*acos(1./k1_all .* log((exp(k1_all)./2)+(exp(-k1_all)./2))));

    subplot(2,2,1)
        scatter(amp_all(ind_all),FWHM(ind_all))
        ylabel('FWHM')
        xlabel('amp_all')
        set(gca,'TickDir','out'); box off; axis square; grid off
    subplot(2,2,2)
        scatter(plaid_corr_all(ind_all),FWHM(ind_all))
        ylabel('FWHM')
        xlabel('plaid_corr')
        set(gca,'TickDir','out'); box off; axis square; grid off
    subplot(2,2,3)
        scatter(mean(ZpZcPWdist_all(:,ind_all),1),FWHM(ind_all))
        ylabel('FWHM')
        xlabel('ZpZcPWdist')
        set(gca,'TickDir','out'); box off; axis square; grid off
        
    legend([leg_str(iA)])
end
print(fullfile(outDir, [svName '_tuningbroadness.pdf']),'-dpdf', '-fillpage') 

figure;
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))];

    x_arr = ones(1,length(ind_all))*iA;

    subplot(2,2,1)
        scatter(x_arr,mean(Zp_all(:,ind_all),1))
        hold on
        scatter(iA,mean(mean(Zp_all(:,ind_all),1)),"black")
        ylabel('mean Zp')
        xlabel('area')
        set(gca,'TickDir','out'); box off; axis square; grid off

    subplot(2,2,2)
        scatter(x_arr,max(Zp_all(:,ind_all),[],1))
        hold on
        scatter(iA,mean(max(Zp_all(:,ind_all),[],1)),"black")
        ylabel('max Zp')
        xlabel('area')
        set(gca,'TickDir','out'); box off; axis square; grid off

    subplot(2,2,3)
        scatter(x_arr,mean(PCI_all(:,ind_all),1))
        hold on
        scatter(iA,mean(mean(PCI_all(:,ind_all),1)),"black")
        ylabel('mean PCI')
        xlabel('area')
        set(gca,'TickDir','out'); box off; axis square; grid off

    subplot(2,2,4)
        scatter(x_arr,max(PCI_all(:,ind_all),[],1))
        hold on
        scatter(iA,mean(max(PCI_all(:,ind_all),[],1)),"black")
        ylabel('max PCI')
        xlabel('area')
        set(gca,'TickDir','out'); box off; axis square; grid off

    legend([leg_str(iA)])
end
print(fullfile(outDir, [svName '_motioninvarianceByArea.pdf']),'-dpdf', '-fillpage') 



%% summary CDFs per experiment
close all; clear all; clc;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';
driver = strvcat('009','010','013','026','040','051','053'); 
area_list = strvcat('V1','V1','V1','V1','V1','V1','V1');
narea = length(area_list);
nCells = [];


figure;
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))

    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))];
    
    subplot(3,2,1)
        cdfplot(Rsq_all(resp_ind,:))
        hold on
        xlabel('fit Rsq')
        legend(leg_str, 'location', 'southeast'); set(gca,'TickDir','out'); box off; grid off
    subplot(3,2,2)
        cdfplot(sse_all_all(resp_ind,:))
        hold on
        xlabel('fit SSE')
        set(gca,'TickDir','out'); box off; grid off
    subplot(3,2,3)
        cdfplot(b_all(resp_ind,:))
        hold on
        xlabel('fit baseline')
        set(gca,'TickDir','out'); box off; grid off
    subplot(3,2,4)
        cdfplot(amp_all(resp_ind,:))
        hold on
        xlabel('fit amplitude')
        set(gca,'TickDir','out'); box off; grid off
    subplot(3,2,5)
        cdfplot(reshape(Zp_all(:,resp_ind),[],1))
        hold on
        xlabel('Zp all')
        set(gca,'TickDir','out'); box off; grid off
    subplot(3,2,6)
        cdfplot(reshape(PCI_all(:,resp_ind),[],1))
        hold on
        xlabel('PCI all')
        set(gca,'TickDir','out'); box off; grid off
end

print(fullfile(outDir, [svName '_SummaryByExperiment.pdf']),'-dpdf', '-fillpage') 


%% Comparing LM 30 deg to 50 deg

close all; clear all; clc;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';
driver = strvcat('019','022','028','030'); 
area_list = strvcat('LM','LM','LM','LM');
narea = length(area_list);
nCells = [];

figure(2);
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))

    subplot(3,2,1)
        scatter(30,length(sig_stim)/nCells)
        hold on
        subtitle('sig stim')
        xlabel('stim degree')
        ylabel('% cells')
    subplot(3,2,2)
        scatter(30,length(sig_dir)/nCells)
        hold on
        subtitle('sig dir')
        xlabel('stim degree')
        ylabel('% cells')
    subplot(3,2,3)
        scatter(30,length(find(DSI_all>0.5))/nCells)
        hold on
        subtitle('DSI>0.5')
        xlabel('stim degree')
        ylabel('% cells')
    subplot(3,2,5)
        h = cdfplot(prefDir_resamp_all);
        hold on
        xlabel('stim degree')
        ylabel('% cells')
        subtitle('30 deg - r, 50 - b')
    
    ind = intersect(intersect(sig_dir,find(DSI_all>0.5)),sig_stim);

    subplot(3,2,4)
        scatter(30,length(ind)/nCells)
        hold on
        subtitle('DSI>0.5')
        xlabel('stim degree')
        ylabel('% cells')

    subplot(3,2,6)
        xax = ones(length(ind),1)+30;
        scatter(xax,prefDir_resamp_all(ind))
        hold on
        xlim([10 70])
        ylim([0 100])
        subtitle('prefDir for intersect ind')
        xlabel('stim degree')
        ylabel('% cells')
end

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';
driver = strvcat('035','039'); 
area_list = strvcat('LM','LM');
narea = length(area_list);
nCells = [];

figure(2);
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    
    subplot(3,2,1)
        scatter(50,length(sig_stim)/nCells)
        hold on
        xlim([10 70])
        ylim([0 1])
        subtitle('sig stim')
        xlabel('stim degree')
        ylabel('% cells')
    subplot(3,2,2)
        scatter(50,length(sig_dir)/nCells)
        hold on
        xlim([10 70])
        ylim([0 1])
        subtitle('sig dir')
        xlabel('stim degree')
        ylabel('% cells')
    subplot(3,2,3)
        scatter(50,length(find(DSI_all>0.5))/nCells)
        hold on
        subtitle('DSI>0.5')
        xlim([10 70])
        ylim([0 1])
        xlabel('stim degree')
        ylabel('% cells')
    subplot(3,2,5)
        h = cdfplot(prefDir_resamp_all);
        hold on
        set(h,'color','r')
        ylim([0 1])
        xlabel('num of resamp bootstraps that match pref dir')
        ylabel('fraction of cells')

    ind = intersect(intersect(sig_dir,find(DSI_all>0.5)),sig_stim);

    subplot(3,2,4)
        scatter(50,length(ind)/nCells)
        hold on
        xlim([10 70])
        ylim([0 1])
        subtitle('meet all 3 criteria')
        xlabel('stim degree')
        ylabel('% cells')

    subplot(3,2,6)
        xax = ones(length(ind),1)+50;
        scatter(xax,prefDir_resamp_all(ind))
        hold on
        xlim([10 70])
        ylim([0 100])
        subtitle('prefDir for intersect ind')
        xlabel('stim degree')
        ylabel('% cells')

end

figure(2); print(fullfile(outDir, [svName '_StimSizeComparisonLM.pdf']),'-dpdf', '-fillpage') 

%% V1 looking at cell inclusion criteria

close all; clear all; clc;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';
driver = strvcat('010','013','014','023','025','026','029'); 
area_list = strvcat('V1','V1','V1','V1','V1','V1','V1');
narea = length(area_list);
nCells = [];

figure(3);
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    leg_str{iA}=[driver(iA,:)];

    subplot(3,2,1)
        scatter(30,length(sig_stim)/nCells)
        hold on
        legend(leg_str)
        subtitle('sig stim')
        xlabel('stim degree')
        ylabel('% cells')
    subplot(3,2,2)
        scatter(30,length(sig_dir)/nCells)
        hold on
        subtitle('sig dir')
        xlabel('stim degree')
        ylabel('% cells')
    subplot(3,2,3)
        scatter(30,length(find(DSI_all>0.5))/nCells)
        hold on
        subtitle('DSI>0.5')
        xlabel('stim degree')
        ylabel('% cells')
    subplot(3,2,5)
        h = cdfplot(prefDir_resamp_all);
        hold on
        xlabel('stim degree')
        ylabel('% cells')
        subtitle('30 deg - r, 50 - b')
    
    ind = intersect(intersect(sig_dir,find(DSI_all>0.5)),sig_stim);

    subplot(3,2,4)
        scatter(30,length(ind)/nCells)
        hold on
        subtitle('ind')
        xlabel('stim degree')
        ylabel('% cells')

    subplot(3,2,6)
        xax = ones(length(ind),1)+iA;
        scatter(xax,prefDir_resamp_all(ind))
        hold on
        ylim([0 100])
        subtitle('prefDir for intersect ind')
        xlabel('stim degree')
        ylabel('% cells')
end

figure(3); print(fullfile(outDir, [svName '_InclusionCriteriaComparisonV1.pdf']),'-dpdf', '-fillpage') 



    
%%  
clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandDirFourPhase_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = length(expt);
frame_rate = 15;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randDir_randPhase';

max_dist = 5;

figure;
for iexp = [4]
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc;
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);

    fprintf([mouse ' ' date '\n'])
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_PatternCorrelationIndexFits.mat']))
    
    subplot(2,3,1)
        cdfplot(sse_all)
        hold on
    subplot(2,3,2)
        cdfplot(b_hat_all)
        hold on
    subplot(2,3,3)
        cdfplot(amp_hat_all)
        hold on
    subplot(2,3,4)
        cdfplot(sse_all(ind))
        hold on
    subplot(2,3,5)
        cdfplot(amp_hat_all(ind))
        hold on
        
    sse_og = sse_all;

    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_PatternCorrelationIndexFits_test2.mat']))

    subplot(2,3,1)
        cdfplot(sse_all)
        hold on
    subplot(2,3,2)
        cdfplot(b_hat_all)
        hold on
    subplot(2,3,3)
        cdfplot(amp_hat_all)
        hold on
    subplot(2,3,4)
        cdfplot(sse_all(ind))
        hold on
    subplot(2,3,5)
        cdfplot(amp_hat_all(ind))
        hold on
        
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_PatternCorrelationIndexFits_test3.mat']))

    subplot(2,3,1)
        cdfplot(sse_all)
        hold on
    subplot(2,3,2)
        cdfplot(b_hat_all)
        hold on
    subplot(2,3,3)
        cdfplot(amp_hat_all)
        hold on
    subplot(2,3,4)
        cdfplot(sse_all(ind))
        hold on
    subplot(2,3,5)
        cdfplot(amp_hat_all(ind))
        hold on
        
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_PatternCorrelationIndexFits_test.mat']))
    
    subplot(2,3,1)
        cdfplot(sse_all)
        hold on
        xlabel('summed square of residuals (SSE)')
        set(gca,'TickDir','out'); box off; axis square; grid off
        legend('sinefit','sinefit b guess','sinefit b mean','sinefit ub/lb')   
    subplot(2,3,2)
        cdfplot(b_hat_all)
        hold on
        xlabel('PCI (Zp-Zc) baseline')
        set(gca,'TickDir','out'); box off; axis square; grid off
    subplot(2,3,3)
        cdfplot(amp_hat_all)
        hold on
        xlabel('PCI (Zp-Zc) amplitude')
        set(gca,'TickDir','out'); box off; axis square; grid off
    subplot(2,3,4)
        cdfplot(sse_all(ind))
        hold on
        xlabel('summed square of residuals (SSE)')
        subtitle('ind only')
        set(gca,'TickDir','out'); box off; axis square; grid off          
   subplot(2,3,5)
        cdfplot(amp_hat_all(ind))
        hold on
        xlabel('PCI (Zp-Zc) amplitude')
        subtitle('ind only')
        set(gca,'TickDir','out'); box off; axis square; grid off        
   subplot(2,3,6)
        scatter(sse_all,sse_og)
        hold on
        xlabel('sse test (b guess)')
        ylabel('sse og')        
        set(gca,'TickDir','out'); box off; axis square; grid off
end
print(fullfile(outDir, [svName '_PCIfitsCompareSinefit.pdf']),'-dpdf', '-fillpage') 


%% shuffled analysis


clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandDirFourPhase_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = length(expt);
frame_rate = 15;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randDir_randPhase';

max_dist = 5;

figure;
for iexp = [4]
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc;
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);
    
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_bootstrapShuffled.mat']))
    shuf_base = boot_base;
    shuf_amp = boot_amp;
    shuf_Zp = boot_Zp;
    shuf_Zc = boot_Zc;

    
    iiexp = 5;
    mouse = expt(iiexp).mouse;
    date = expt(iiexp).date;
    area = expt(iiexp).img_loc;
    ImgFolder = expt(iiexp).coFolder;
    time = expt(iiexp).coTime;
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_bootstrapFits.mat']))
    
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc;
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;

    fprintf([mouse ' ' date '\n'])
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_PatternCorrelationIndexFits.mat']))
    
    
    subplot(3,3,1)
        cdfplot(sse_all)
        hold on
        xlabel('summed square of residuals (SSE)')
        set(gca,'TickDir','out'); box off; grid off
        legend('data','shuffled')

    subplot(3,3,2)
        cdfplot(b_hat_all)
        hold on
        for ib = 1:100
            h(1,1) = cdfplot(boot_base(:,ib));
            h(2,2) = cdfplot(shuf_base(:,ib));
            set(h(1,1),'Color',[0.7 0.7 0.7]);
            set(h(2,2),'Color',[0.4940 0.1840 0.5560]);
        end
        xlabel('PCI (Zp-Zc) baseline')
        set(gca,'TickDir','out'); box off; grid off
        
    subplot(3,3,3)
        cdfplot(amp_hat_all)
        hold on
        for ib = 1:100
            h(1,2) = cdfplot(boot_amp(:,ib));
            h(2,2) = cdfplot(shuf_amp(:,ib));
            set(h(1,2),'Color',[0.7 0.7 0.7]);
            set(h(2,2),'Color',[0.4940 0.1840 0.5560]);
        end
        xlabel('PCI (Zp-Zc) amplitude')
        set(gca,'TickDir','out'); box off; grid off
  
        boot_PCI = boot_Zp-boot_Zc;
        shuf_PCI = shuf_Zp-shuf_Zc;

    subplot(3,3,4)
        cdfplot(reshape(PCI,[],1))
        hold on
        h(1,2) = cdfplot(reshape(boot_PCI,[],1));
        h(2,2) = cdfplot(reshape(shuf_PCI,[],1));
        set(h(1,2),'Color',[0.7 0.7 0.7]);
        set(h(2,2),'Color',[0.4940 0.1840 0.5560]);
        xlabel('PCI at each phase (n=4*nCells)')
        set(gca,'TickDir','out'); box off; grid off
        
    subplot(3,3,5)
        cdfplot(reshape(Zp,[],1))
        hold on
        h(1,2) = cdfplot(reshape(boot_Zp,[],1));
        h(2,2) = cdfplot(reshape(shuf_Zp,[],1));
        set(h(1,2),'Color',[0.7 0.7 0.7]);
        set(h(2,2),'Color',[0.4940 0.1840 0.5560]);
        xlabel('Zp at each phase (n=4*nCells)')
        set(gca,'TickDir','out'); box off; grid off
        
    subplot(3,3,6)
        cdfplot(reshape(Zc,[],1))
        hold on
        h(1,2) = cdfplot(reshape(boot_Zc,[],1));
        h(2,2) = cdfplot(reshape(shuf_Zc,[],1));
        set(h(1,2),'Color',[0.7 0.7 0.7]);
        set(h(2,2),'Color',[0.4940 0.1840 0.5560]);
        xlabel('Zc at each phase (n=4*nCells)')
        set(gca,'TickDir','out'); box off; grid off
        hold on
        
    subplot(3,3,7)
        scatter(reshape(Zc,[],1),reshape(Zp,[],1))
        hold on
        xlabel('Zc')
        ylabel('Zp')
        xlim([-8 8])
        ylim([-8 8])
        plotZcZpBorders
        
    subplot(3,3,8)
        scatter(b_hat_all, amp_hat_all)
        hold on
        xlabel('PCI fit baseline')
        ylabel('PCI fit amplitude')
        subtitle('Data - Exp 4')
        
    subplot(3,3,9)
        histogram(mean(boot_amp,1),'FaceColor',[0.7 0.7 0.7])
        hold on
        histogram(mean(shuf_amp,1),'FaceColor',[0.4940 0.1840 0.5560])
        xline(mean(amp_hat_all))
        xlabel('mean PCI amplitude')
        subtitle('Data - Exp 4')

end
print(fullfile(outDir, [svName '_PCIfitsCompareShuffled.pdf']),'-dpdf', '-fillpage') 




%% plot k tuning value

close all; clear all; clc;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';
driver = strvcat('SLC'); 
area = 'all_areas';
area_list = strvcat('V1');
narea =1;%narea = length(area_list);
nCells = [];


figure;
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(ind_all))];
       
    subplot(3,3,1)
        scatter(k1_all(ind_all),Zp_all(1,ind_all))
        xlabel('k')
        ylabel('Zp, 1 phase'); set(gca,'TickDir','out')
%         xlim([0 20])
        
    subplot(3,3,2)
        scatter(k1_all(ind_all),mean(Zp_all(:,ind_all)))
        xlabel('k')
        ylabel('mean Zp'); set(gca,'TickDir','out')  
%         xlim([0 20])
        
    subplot(3,3,3)
        scatter(k1_all(ind_all),mean(PCI_all(:,ind_all)))
        xlabel('k')
        ylabel('mean PCI'); set(gca,'TickDir','out') 
%         xlim([0 20])
        
   subplot(3,3,4)
        scatter(k1_all(ind_all),Zc_all(1,ind_all))
        xlabel('k')
        ylabel('Zc, 1 phase'); set(gca,'TickDir','out')
%         xlim([0 20])
   subplot(3,3,5)
        scatter(k1_all(ind_all),amp_all(ind_all))
        xlabel('k')
        ylabel('PCI fit amp'); set(gca,'TickDir','out')
%         xlim([0 20])

   subplot(3,3,6)
        scatter(k1_all(ind_all),b_all(ind_all))
        xlabel('k')
        ylabel('PCI fit baseline'); set(gca,'TickDir','out')
%         xlim([0 20])
end
print(fullfile(outDir, [svName '_DirectionFits.pdf']),'-dpdf', '-fillpage') 



%% Pairwise euclidean distance of ZpZc compared to shuffled
clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandDirFourPhase_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = length(expt);
frame_rate = 15;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randDir_randPhase';

max_dist = 5;

iexp = [5];
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc;
ImgFolder = expt(iexp).coFolder;
time = expt(iexp).coTime;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

h_all = [];
p_all = [];

load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_ZpZc_pairwiseDist.mat']))
ZpZcPWdist_sh = ZpZcPWdist;
Zp_sh = Zp;
Zc_sh = Zc;

figure;
    subplot(1,2,1)
    	cdfplot(reshape(ZpZcPWdist_sh,[],1))
        hold on
    subplot(1,2,2)
    	cdfplot(mean(ZpZcPWdist_sh,1))
        hold on 
        
start=1;
for iexp = [19];
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc;
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);
    
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_ZpZc_pairwiseDist.mat']))
    
    [h, p] =  ttest2(mean(ZpZcPWdist_sh),mean(ZpZcPWdist));
        h_all = [h_all; h];
        p_all = [p_all; p];
    leg_str{start} = num2str(p_all(start,:));
     
%         leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(ind_all))];

    subplot(1,2,1)
        cdfplot(reshape(ZpZcPWdist,[],1))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('PW dist of Zp Zc, all 6 for all cells')
    subplot(1,2,2)
        cdfplot(mean(ZpZcPWdist,1))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('PW dist of Zp Zc, avg per cell')
%         legend(leg_str, 'location','southeast');
start=start+1;
end  
print(fullfile(outDir, [svName '_PWdistance.pdf']),'-dpdf', '-fillpage') 


figure;
start=1;
for iCell = 1:20
    subplot(5,4,start)
        for im = 1:4
            scatter(Zc_sh(im,iCell), Zp_sh(im,iCell))
            hold on
        end
        ylabel('Zp'); ylim([-4 8]);
        xlabel('Zc'); xlim([-4 8]);
        if iCell ==1; legend('0 deg','90 deg','180 deg', '270 deg'); end;
        subtitle(['cell ' num2str(iCell)])
        plotZcZpBorders
    start = start+1;    
end
print(fullfile(outDir, [svName '_ZpZcPlots_OneToFourPhase.pdf']),'-dpdf', '-fillpage') 

%% Comparing 2 contrast to 4 phase

clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandDirFourPhase_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = length(expt);
frame_rate = 15;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randDir_randPhase';
max_dist = 5;

iexp = [40];
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc;
ImgFolder = expt(iexp).coFolder;
time = expt(iexp).coTime;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_stimData.mat']))
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_plaidCorr.mat']))

ind = intersect(intersect(DSI_ind,p_dir),resp_ind_dir);

%plot plaid corr for ind

figure;
    subplot(2,2,1)
        cdfplot(plaid_corr(ind))
        hold on
        

clear all
ds = 'CrossOriRandDirTwoCon_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = length(expt);
frame_rate = 15;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randDir_randPhase';
max_dist = 5;

iexp = [1];
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc;
ImgFolder = expt(iexp).coFolder;
time = expt(iexp).coTime;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_respData.mat']))
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_stimData.mat']))
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_plaidCorr.mat']))

DSI_ind = find(DSIL>0.5);
ind = intersect(intersect(DSI_ind,p_all),resp_ind);

    subplot(2,2,1)
        cdfplot(plaid_corr(ind))
        hold on
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('Avg plaid correlation')

print(fullfile(outDir, [svName '_TwoConComparedtoFourPhase.pdf']),'-dpdf', '-fillpage') 

%create ind, plot plaidcorr



%% compare variance and fit amplitude (scatter of all cells), for all areas
close all; clear all; clc;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';
driver = strvcat('SLC', 'Scn', 'SLC', 'SLC'); 
area = 'all_areas';
area_list = strvcat('V1', 'V1', 'LM', 'AL');
narea = length(area_list);
nCells = [];


figure;
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))];
    
    mean_avg_all(mean_avg_all<0) = 0;

    subplot(4,4,iA)
        scatter(amp_all(resp_ind),(std_avg_all(resp_ind)./mean_avg_all(resp_ind)))
        hold on
        xlabel('Fit amplitude')
        ylabel('variance/mean')
        set(gca,'TickDir','out'); box off; axis square; grid off
        legend(leg_str{iA}, 'location', 'northeast'); set(gca,'TickDir','out'); box off; grid off
    subplot(4,4,iA+4)
        scatter(mean_avg_all(resp_ind),std_avg_all(resp_ind))
        hold on
        refline(1,0)
        xlabel('mean')
        ylabel('variance')
        set(gca,'TickDir','out'); box off; axis square; grid off
        legend(leg_str{iA}, 'location', 'northeast'); set(gca,'TickDir','out'); box off; grid off
         
end
print(fullfile(outDir, [svName '_VarianceByFitAmp.pdf']),'-dpdf', '-fillpage') 


%% layer 4 by experiment

close all; clear all; clc;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';
exptn = strvcat('063','064','107','109','113','114','115','116'); 
area = 'all_areas';
area_list = strvcat('V1','V1','V1','V1','V1','V1','V1','V1');
narea = length(area_list);
nCells = [];



figure;
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' exptn(iA,:) '.mat'])))
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end
    leg_str{iA}=[area_list(iA,:) ' ' exptn(iA,:) ' n=' num2str(length(resp_ind))];
    

    c(1,:)=[0.6350 0.0780 0.1840]; c(2,:)=[0.8500 0.3250 0.0980];
    c(3,:)=[0 0.4470 0.7410]; c(4,:)=[0.3010 0.7450 0.9330];

    mean_avg_all(mean_avg_all<0) = 0;

if iA<5 % Only plot first 4 experiments
    subplot(4,4,1)
        h(1,iA)=cdfplot(mean_avg_all(resp_ind));
        hold on
        set(h(1,iA),'Color',c(iA,:))
        xlabel('mean df/f')
        set(gca,'TickDir','out'); box off; axis square; grid off
     subplot(4,4,2)
        h(2,iA)=cdfplot(std_avg_all(resp_ind));
        hold on
        set(h(2,iA),'Color',c(iA,:))
        xlabel('variance (std)')
        set(gca,'TickDir','out'); box off; axis square; grid off 
    subplot(4,4,3)
        h(3,iA)=cdfplot(std_avg_all(resp_ind)./mean_avg_all(resp_ind));
        hold on
        set(h(3,iA),'Color',c(iA,:))
        xlabel('variance/mean')
        set(gca,'TickDir','out'); box off; axis square; grid off
    subplot(4,4,4)
        scatter(amp_all(resp_ind),std_avg_all(resp_ind)./mean_avg_all(resp_ind))
        hold on
        ylabel('variance/mean')
        xlabel('amp')
        set(gca,'TickDir','out'); box off; axis square; grid off
     
 
    subplot(4,4,5)
        h(5,iA)=cdfplot(b_all(resp_ind));
        hold on
        set(h(5,iA),'Color',c(iA,:))
        xlabel('phase mod baseline')
        set(gca,'TickDir','out'); box off; axis square; grid off
    subplot(4,4,6)
        h(6,iA)=cdfplot(amp_all(resp_ind));
        hold on
        set(h(6,iA),'Color',c(iA,:))
        xlabel('phase mod amplitude')
        set(gca,'TickDir','out'); box off; axis square; grid off
    subplot(4,4,7)
        h(7,iA)=cdfplot(plaid_corr_all(resp_ind));
        hold on
        set(h(7,iA),'Color',c(iA,:))
        xlabel('avg corr across plaid phase')
        set(gca,'TickDir','out'); box off; axis square; grid off
    subplot(4,4,8)
        h(8,iA)=cdfplot(mean(ZpZcPWdist_all(:,resp_ind),1));
        hold on
        set(h(8,iA),'Color',c(iA,:))
        xlabel('avg ZpZc dist')
        set(gca,'TickDir','out'); box off; axis square; grid off
end

% Now plot all amplitude modulations
    subplot(4,4,9)
        cdfplot(amp_all(resp_ind));
        hold on
        xlabel('phase mod amplitude')
        set(gca,'TickDir','out'); box off; axis square; grid off
end
sgtitle('L4 single cell comparison - local injection FLEX-6s (red), transgenic 8m (blue)')
print(fullfile(outDir, [svName '_Layer4_InjectionTransgenicComparison.pdf']),'-dpdf', '-fillpage') 

