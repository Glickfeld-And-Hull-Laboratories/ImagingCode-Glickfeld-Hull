close all; clear all; clc;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = ([base '\Analysis\Neuropixel\CrossOri\randDirFourPhase\summaries']);
outDir = ([base '\Analysis\Neuropixel\CrossOri\randDirFourPhase']);
svName = 'randDirFourPhase_CrossOri';


expts = strvcat('g01', 'g06', 'g12', 'g17'); %Four acute penetrations in V1, 1 marmoset
nexp = length(expts);

sig_stim = [];
sig_dir = [];
totCells = [];

DSI_all = [];
Zp_all = [];
Zc_all = [];
PCI_all = [];
ZpZcPWdist_all = [];

yfit_all_all = [];
b_all = [];
amp_all = [];
pha_all = [];
per_all = [];
sse_all_all = [];
Rsq_all  = [];

k1_all = [];
dir_Rsq_all = [];
dir_sse_all = [];
dir_yfits = [];

plaid_corr_all = [];


start=1;
for iexp=1:nexp

    load(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_fitsSG.mat']))

    plaid_corr_all = [plaid_corr_all; plaid_corr'];
    
    b_all = [b_all; b_hat_all];
    amp_all = [amp_all; amp_hat_all];
    Rsq_all = [Rsq_all; R_square_all];
    sse_all_all = [sse_all_all; sse_all];
    
    yfit_all_all = cat(1,yfit_all_all, yfit_all);
    
    if start > 1
        sig_stim = [sig_stim; resp_ind_dir + totCells]; 
        sig_dir = [sig_dir, p_dir + totCells];   
    else
        sig_stim = resp_ind_dir; 
        sig_dir = p_dir;   
    end
    
    if start > 1
        totCells = totCells + nCells;    
    else
        totCells = nCells;
    end
    
    Zp_all = cat(2,Zp_all, Zp);
    Zc_all = cat(2,Zc_all, Zc);
    ZpZcPWdist_all = cat(2,ZpZcPWdist_all, ZpZcPWdist);
    
    DSI_all = [DSI_all; DSI'];


    k1_all = [k1_all; k1_hat_all];
    dir_Rsq_all = [dir_Rsq_all; dir_R_square_all];
    dir_sse_all = [dir_sse_all; dir_sse_all];

start=start+1;
end

    save(fullfile(summaryDir,[svName '_Summary_V1_MAR.mat']),  'sig_dir','sig_stim','plaid_corr_all','ZpZcPWdist_all','k1_all','dir_Rsq_all','dir_sse_all','totCells','Zp_all','Zc_all','b_all','amp_all','Rsq_all','yfit_all_all','sse_all_all', 'DSI_all','expts')

%%


resp_ind = intersect(sig_stim,find(DSI_all>0.5));

    
figure(1);
    inclusion_crit(1) = size(sig_stim,1)/totCells;
    % inclusion_crit(2) = size(sig_dir,2)/totCells;
    inclusion_crit(2) = size(find(DSI_all>0.5),1)/totCells;
    inclusion_crit(3) = length(resp_ind)/totCells;
    subplot(4,4,1)
        bar(inclusion_crit);
        hold on
        ylim([0 1])
        set(gca,'TickDir','out'); box off; axis square; grid off
        xlabel('Inclusion Criteria (ttest, DSI>.5, int)')
    subtitle('marmoset V1')
print(fullfile(outDir, [svName '_marmoset_Summary.pdf']),'-dpdf', '-fillpage') 


    
figure(2); 
    ind1 = intersect(find(Zp_all(1,:)>1.28),find(Zp_all(1,:)-Zc_all(1,:)>1.28));
    ind2 = intersect(find(Zp_all(2,:)>1.28),find(Zp_all(2,:)-Zc_all(2,:)>1.28));
    ind4 = intersect(find(Zp_all(4,:)>1.28),find(Zp_all(4,:)-Zc_all(4,:)>1.28));
    pattern1 = intersect(ind1,resp_ind);
    pattern2 = intersect(ind2,resp_ind);
    pattern4 = intersect(ind4,resp_ind);    
    subplot(5,5,1)
        scatter(Zc_all(1,resp_ind),Zp_all(1,resp_ind),5,[0.7 0.7 0.7],'filled')
        hold on
        scatter(Zc_all(1,pattern4),Zp_all(1,pattern4),5,'filled')
        xlabel('Zc'); ylabel('Zp'); xlim([-4 8]); ylim([-4 8]); plotZcZpBorders; set(gca,'TickDir','out'); axis square
   subplot(5,5,2)
        scatter(Zc_all(2,resp_ind),Zp_all(2,resp_ind),5,[0.7 0.7 0.7],'filled')
        hold on
        scatter(Zc_all(2,pattern4),Zp_all(2,pattern4),5,'filled')
        xlabel('Zc'); ylabel('Zp'); xlim([-4 8]); ylim([-4 8]); plotZcZpBorders; set(gca,'TickDir','out'); axis square
   subplot(5,5,3)
        scatter(Zc_all(3,resp_ind),Zp_all(3,resp_ind),5,[0.7 0.7 0.7],'filled')
        hold on
        scatter(Zc_all(3,pattern4),Zp_all(3,pattern4),5,'filled')
        xlabel('Zc'); ylabel('Zp'); xlim([-4 8]); ylim([-4 8]); plotZcZpBorders; set(gca,'TickDir','out'); axis square
   subplot(5,5,4)
        scatter(Zc_all(4,resp_ind),Zp_all(4,resp_ind),5,[0.7 0.7 0.7],'filled')
        hold on
        scatter(Zc_all(4,pattern4),Zp_all(4,pattern4),5,'filled')
        xlabel('Zc'); ylabel('Zp'); xlim([-4 8]); ylim([-4 8]); plotZcZpBorders; set(gca,'TickDir','out'); axis square
print(fullfile(outDir, [svName '_marmoset_SummaryZpZc.pdf']),'-dpdf', '-fillpage') 


figure(3);
    pattern_ind = (Zp_all-Zc_all);    
    pattind = mean(pattern_ind,1);
    pattpeak = max(pattern_ind,[],1);   
    subplot(3,3,1)
        scatter((-1*amp_all(resp_ind)),pattpeak(resp_ind),7,'filled')
        hold on
        xlabel('Spatial invariance (-Amplitude)')
        ylabel('Direction invariance (Fit peak)')
        set(gca,'TickDir','out'); axis square
        subtitle(['marmoset V1'])
print(fullfile(outDir, [svName '_marmoset_Summary_likeNicholas.pdf']),'-dpdf', '-fillpage') 



%% pie chart

figure;

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


    pie0 = [ind_patt(1)/length(resp_ind) ind_comp(1)/length(resp_ind) (1-(ind_patt(1)/length(resp_ind))-(ind_comp(1)/length(resp_ind)))];
    pie90 = [ind_patt(2)/length(resp_ind) ind_comp(2)/length(resp_ind) (1-(ind_patt(2)/length(resp_ind))-(ind_comp(2)/length(resp_ind)))];
    pie270 = [ind_patt(4)/length(resp_ind) ind_comp(4)/length(resp_ind) (1-(ind_patt(4)/length(resp_ind))-(ind_comp(4)/length(resp_ind)))];

    ind_patt_cells0 = intersect(intersect(find(Zp_all(1,:)>1.28),find(Zp_all(1,:)-Zc_all(1,:)>1.28)),resp_ind);
    ind_patt_cells0_pattAt90 = intersect(intersect(find(Zp_all(2,:)>1.28),find(Zp_all(2,:)-Zc_all(2,:)>1.28)),ind_patt_cells0);
    ind_patt_cells0_compAt90 = intersect(intersect(find(Zc_all(2,:)>1.28),find(Zc_all(2,:)-Zp_all(2,:)>1.28)),ind_patt_cells0);
    ind_patt_cells0_pattAt180 = intersect(intersect(find(Zp_all(3,:)>1.28),find(Zp_all(3,:)-Zc_all(3,:)>1.28)),ind_patt_cells0);
    ind_patt_cells0_compAt180 = intersect(intersect(find(Zc_all(3,:)>1.28),find(Zc_all(3,:)-Zp_all(3,:)>1.28)),ind_patt_cells0);
    ind_patt_cells0_pattAt270 = intersect(intersect(find(Zp_all(4,:)>1.28),find(Zp_all(4,:)-Zc_all(4,:)>1.28)),ind_patt_cells0);
    ind_patt_cells0_compAt270 = intersect(intersect(find(Zc_all(4,:)>1.28),find(Zc_all(4,:)-Zp_all(4,:)>1.28)),ind_patt_cells0);
    

    ind_patt_cells270 = intersect(intersect(find(Zp_all(4,:)>1.28),find(Zp_all(4,:)-Zc_all(4,:)>1.28)),resp_ind);
    ind_patt_cells270_pattAt90 = intersect(intersect(find(Zp_all(2,:)>1.28),find(Zp_all(2,:)-Zc_all(2,:)>1.28)),ind_patt_cells270);
    ind_patt_cells270_compAt90 = intersect(intersect(find(Zc_all(2,:)>1.28),find(Zc_all(2,:)-Zp_all(2,:)>1.28)),ind_patt_cells270);
    ind_patt_cells270_pattAt180 = intersect(intersect(find(Zp_all(3,:)>1.28),find(Zp_all(3,:)-Zc_all(3,:)>1.28)),ind_patt_cells270);
    ind_patt_cells270_compAt180 = intersect(intersect(find(Zc_all(3,:)>1.28),find(Zc_all(3,:)-Zp_all(3,:)>1.28)),ind_patt_cells270);
    ind_patt_cells270_pattAt0 = intersect(intersect(find(Zp_all(1,:)>1.28),find(Zp_all(4,:)-Zc_all(1,:)>1.28)),ind_patt_cells270);
    ind_patt_cells270_compAt0 = intersect(intersect(find(Zc_all(1,:)>1.28),find(Zc_all(4,:)-Zp_all(1,:)>1.28)),ind_patt_cells270);

    pie_patt270at90 = [length(ind_patt_cells270_pattAt90)/length(ind_patt_cells270) length(ind_patt_cells270_compAt90)/length(ind_patt_cells270) (1-(length(ind_patt_cells270_pattAt90)/length(ind_patt_cells270))-( length(ind_patt_cells270_compAt90)/length(ind_patt_cells270)))];
    pie_patt270at180 = [length(ind_patt_cells270_pattAt180)/length(ind_patt_cells270) length(ind_patt_cells270_compAt180)/length(ind_patt_cells270) (1-(length(ind_patt_cells270_pattAt180)/length(ind_patt_cells270))-( length(ind_patt_cells270_compAt180)/length(ind_patt_cells270)))];
    pie_patt270at0 = [length(ind_patt_cells270_pattAt0)/length(ind_patt_cells270) length(ind_patt_cells270_compAt0)/length(ind_patt_cells270) (1-(length(ind_patt_cells270_pattAt0)/length(ind_patt_cells270))-( length(ind_patt_cells270_compAt0)/length(ind_patt_cells270)))];

    pie_patt0atOthers(1,:) = [pie_patt270at90];
    pie_patt0atOthers(2,:) = [pie_patt270at180];
    pie_patt0atOthers(3,:) = [pie_patt270at0];
    pie_patt0atOtherAvg = mean(pie_patt0atOthers,1);
    

    subplot(5,5,1)
        piechart(pie270)
        % subtitle('pie 270')
    subplot(5,5,2)
        piechart(pie_patt270at90)
        % subtitle('patt 0 at 90')
    subplot(5,5,3)
        piechart(pie_patt270at180)
        % subtitle('patt 0 at 180')
    subplot(5,5,4)
        piechart(pie_patt270at0)
        % subtitle('patt 0 at 270')
    subplot(5,5,5)
        piechart(pie_patt0atOtherAvg)
        % subtitle('patt 0 at others')
    subplot(5,5,6)
        piechart([max_classification(1) 1-max_classification(1)])


    subplot(5,5,11)
        piechart(pie90)
        % subtitle('pie 90')
    
        
    print(fullfile(outDir, [svName '_piecharts.pdf']),'-dpdf', '-fillpage') 



  
