close all; clear all; clc;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'PlaidDiscrimSummary');
svName = 'plaidDiscrim';
area_list = {'V1','LM'};%,'RL','PM'};
col_mat = {'k','c'};
narea = size(area_list,2);
stimOn = 25000;
stimSz = 1000;
stimSetDur = 'LT';
stimSetSz = 'FF';
Zc_all_all = [];
Zp_all_all = [];
C_stim_all_all = [];
leg_str = cell(3,narea);

for i = 1:narea
    area = area_list{i};
    load(fullfile(summaryDir,[svName '_Summary_' area '_stimOn' stimSetDur num2str(stimOn) '_' stimSetSz '.mat']))
    
    Zc_all_all = [Zc_all_all Zc_all];
    Zp_all_all = [Zp_all_all Zp_all];
    C_stim_all_all = [C_stim_all_all C_stim_all];
    x = (sum(stim_resp_dir_all(:,:,:,1),2)./size(stim_resp_dir_all(:,:,:,1),2)).^2;
    y = sum(((stim_resp_dir_all(:,:,:,1).^2)./size(stim_resp_dir_all(:,:,:,1),2)),2);
    z = 1-(1/size(stim_resp_dir_all(:,:,:,1),2));
    S_all = squeeze((1-(x./y))./z);
    pyr_resp_use = setdiff(resp_ind_all,red_cells_all);
    in_resp_use = intersect(resp_ind_all,red_cells_all);
    figure(1)
    subplot(2,2,1)
    cdfplot(Zc_all(pyr_resp_use))
    hold on
    xlabel('Zc')
    title('')
    axis square
    
    subplot(2,2,2)
    cdfplot(Zp_all(pyr_resp_use))
    hold on
    xlabel('Zp')
    title('')
    axis square
    
%     subplot(2,3,3)
%     cdfplot(Zc_all(pyr_resp_use)-Zp_all(pyr_resp_use))
%     hold on
%     xlabel('Zc-Zp')
%     title('')
%     axis square
%     
%     subplot(2,3,4)
%     cdfplot(abs(C_stim_all(2,pyr_resp_use)))
%     hold on
%     xlabel('abs(Target weight)')
%     title('')
%     axis square
    
    subplot(2,2,4)
    cdfplot(C_stim_all(3,pyr_resp_use))
    hold on
    xlabel('Mask weight')
    title('')
    axis square
    
    subplot(2,2,3)
    cdfplot(abs(C_stim_all(5,pyr_resp_use)))
    hold on
    xlabel('abs(Choice weight)')
    title('')
    axis square
    
    figure(2)
    subplot(3,2,1)
    cdfplot(S_all(pyr_resp_use,1))
    title('Grating')
    xlim([0 1.25])
    xlabel('Sparseness')
    hold on
    subplot(3,2,2)
    cdfplot(S_all(pyr_resp_use,2))
    title('Plaid')
    xlim([0 1.25])
    xlabel('Sparseness')
    hold on
    subplot(3,2,3)
    errorbar(dirs, mean(stim_resp_dir_all(pyr_resp_use,:,1,1),1),std(stim_resp_dir_all(pyr_resp_use,:,1,1),[],1)./sqrt(size(stim_resp_dir_all,1)))
    hold on
    ylim([0 .2])
    ylabel('dF/F')
    xlim([-10 100])
    subplot(3,2,4)
    errorbar(dirs, mean(stim_resp_dir_all(pyr_resp_use,:,2,1),1),std(stim_resp_dir_all(pyr_resp_use,:,2,1),[],1)./sqrt(size(stim_resp_dir_all,1)))
    hold on
    ylim([0 .2])
    xlim([-10 100])
    ylabel('dF/F')
    
    subplot(3,2,5)
    cdfplot(sum(h_dir_all(pyr_resp_use,:,1),2))
    title('Grating')
    xlabel('Number of significant stimuli')
    hold on
    subplot(3,2,6)
    cdfplot(sum(h_dir_all(pyr_resp_use,:,2),2))
    title('Plaid')
    xlabel('Number of significant stimuli')
    hold on
    
    figure(5)
    subplot(2,2,1)
    [values, edges] = histcounts(sum(h_dir_all(pyr_resp_use,:,1),2), [0:1:9],'Normalization', 'probability');
    plot(0:8, values)
    hold on
    subplot(2,2,2)
    [values, edges] = histcounts(sum(h_dir_all(pyr_resp_use,:,2),2), [0:1:9],'Normalization', 'probability');
    plot(0:8, values)
    hold on
    subplot(2,2,3)
    [values, edges] = histcounts(sum(h_dir_all(pyr_resp_use,:,1),2), [0:1:9]);
    plot(0:8, values)
    hold on
    subplot(2,2,4)
    [values, edges] = histcounts(sum(h_dir_all(pyr_resp_use,:,2),2), [0:1:9]);
    plot(0:8, values)
    hold on
    
    
    figure(3)
    train_all_mat = [mean(alltrain_traintest_pctCorr_all,2) mean(alltrain_alltest_pctCorr_all,2)];
    train_grating_mat = [mean(gratingtrain_traintest_pctCorr_all,2) mean(gratingtrain_gratingtest_pctCorr_all,2) mean(gratingtrain_plaidtest_pctCorr_all,2)];
    train_plaid_mat = [mean(plaidtrain_traintest_pctCorr_all,2) mean(plaidtrain_plaidtest_pctCorr_all,2) mean(plaidtrain_gratingtest_pctCorr_all,2)];
    subplot(3,1,1)
    plot(1:2, train_all_mat',col_mat{i})
    hold on
    errorbar(1:2, mean(train_all_mat,1)', std(train_all_mat,[],1)./sqrt(size(train_all_mat,1)),['o' col_mat{i}])
    subplot(3,1,2)
    plot(1:3, train_grating_mat',col_mat{i})
    hold on
    errorbar(1:3, mean(train_grating_mat,1)', std(train_grating_mat,[],1)./sqrt(size(train_grating_mat,1)),['o' col_mat{i}])
    subplot(3,1,3)
    plot(1:3, train_plaid_mat',col_mat{i})
    hold on
    errorbar(1:3, mean(train_plaid_mat,1)', std(train_plaid_mat,[],1)./sqrt(size(train_plaid_mat,1)),['o' col_mat{i}])
    
    figure(4)
    scatter(Zc_all(pyr_resp_use),Zp_all(pyr_resp_use))
    hold on
    
    %     figure(2)
%     subplot(2,3,1)
%     cdfplot(Zc_all(in_resp_use))
%     hold on
%     xlabel('Zc')
%     title('')
%     axis square
%     
%     subplot(2,3,2)
%     cdfplot(Zp_all(in_resp_use))
%     hold on
%     xlabel('Zp')
%     title('')
%     axis square
%     
%     subplot(2,3,3)
%     cdfplot(Zc_all(in_resp_use)-Zp_all(in_resp_use))
%     hold on
%     xlabel('Zc-Zp')
%     title('')
%     axis square
%     
%     subplot(2,3,4)
%     cdfplot(abs(C_stim_all(2,in_resp_use)))
%     hold on
%     xlabel('abs(Target weight)')
%     title('')
%     axis square
%     
%     subplot(2,3,5)
%     cdfplot(C_stim_all(3,in_resp_use))
%     hold on
%     xlabel('Mask weight')
%     title('')
%     axis square
%     
%     subplot(2,3,6)
%     cdfplot(abs(C_stim_all(5,in_resp_use)))
%     hold on
%     xlabel('abs(Choice weight)')
%     title('')
%     axis square
    
    nexp = sum(totCells>0);
    leg_str{1,i} = [area '- ' num2str(length(pyr_resp_use)) ', ' num2str(nexp)];
    leg_str{2,i} = [area '- ' num2str(length(in_resp_use)) ', ' num2str(nexp)];
    
%     resp_use = intersect(resp_ind_all, red_cells_all);
%     subplot(3,3,4)
%     cdfplot(Zc_all(resp_use))
%     hold on
%     xlabel('Zc')
%     title('')
%     
%     subplot(3,3,5)
%     cdfplot(Zp_all(resp_use))
%     hold on
%     xlabel('Zp')
%     title('')
%     
%     subplot(3,3,6)
%     cdfplot(C_stim_all(3,resp_use))
%     hold on
%     xlabel('Mask weight')
%     title('')
%     
%     nexp = sum(totCells>0);
%     leg_str{2,i} = [area '- ' num2str(length(resp_use)) ', ' num2str(nexp)];
%     
%     resp_use = setdiff(resp_ind_all,red_cells_all);
%     subplot(3,3,7)
%     cdfplot(Zc_all(resp_use))
%     hold on
%     xlabel('Zc')
%     title('')
%     
%     subplot(3,3,8)
%     cdfplot(Zp_all(resp_use))
%     hold on
%     xlabel('Zp')
%     title('')
%     
%     subplot(3,3,9)
%     cdfplot(C_stim_all(3,resp_use))
%     hold on
%     xlabel('Mask weight')
%     title('')
%     
%     nexp = sum(totCells>0);
%     leg_str{3,i} = [area '- ' num2str(length(resp_use)) ', ' num2str(nexp)];
    
%     figure(3)
%     movegui('center')
%     subplot(2,2,i)
%     scatter(Zc_all(pyr_resp_use),Zp_all(pyr_resp_use))
%     hold on
%     scatter(Zc_all(in_resp_use),Zp_all(in_resp_use))
%     xlabel('Zc')
%     ylabel('Zp')
%     title('')
%     ylim([-8 8])
%     xlim([-8 8])
%     plotZcZpBorders
%     axis square
%     title(area_list{i})

    
end

figure(1)
subplot(2,2,1)
legend(area_list,'location','southeast')
sgtitle('Pyr cells')
print(fullfile(summaryDir,[svName '_' cell2mat(area_list) '_Summary_stimOn' stimSetDur num2str(stimOn) '_' stimSetSz '_pyr.pdf']),'-dpdf','-fillpage')
figure(2)
subplot(3,2,1)
legend(area_list,'location','southeast')
sgtitle('Pyr cells')
print(fullfile(summaryDir,[svName '_' cell2mat(area_list) '_Summary_stimOn' stimSetDur num2str(stimOn) '_' stimSetSz '_pyrSelectivity.pdf']),'-dpdf','-fillpage')

figure(3)
subplot(3,1,1)
ylim([0 1])
set(gca, 'XTick',1:2,'XTickLabel',{'Train', 'Test'})
xlim([0 3])
ylabel('Percent Correct- Stimulus direction')
title('Train all')
subplot(3,1,2)
ylim([0 1])
set(gca, 'XTick',1:3,'XTickLabel',{'Train', 'Test grating', 'Test plaid'}) 
xlim([0 4])
ylabel('Percent Correct- Stimulus direction')
title('Train grating')
subplot(3,1,3)
ylim([0 1])
set(gca, 'XTick',1:3,'XTickLabel',{'Train', 'Test plaid', 'Test grating'}) 
xlim([0 4])
ylabel('Percent Correct- Stimulus direction')
title('Train plaid')
% print(fullfile(summaryDir,[svName '_' cell2mat(area_list) '_Summary_stimOn' stimSetDur num2str(stimOn) '_' stimSetSz '_popDecoder_8pcs.pdf']),'-dpdf','-fillpage')

figure(4)
xlabel('Zc')
ylabel('Zp')
axis square
title('')
ylim([-8 8])
xlim([-8 8])
plotZcZpBorders
% print(fullfile(summaryDir,[svName '_' cell2mat(area_list) '_Summary_stimOn' stimSetDur num2str(stimOn) '_' stimSetSz '_scatter.pdf']),'-dpdf','-fillpage')

figure(5)
subplot(2,2,1)
ylim([0 0.6])
title('Gratings')
xlabel('Number of stimuli')
ylabel('Fraction of cells')
subplot(2,2,2)
ylim([0 0.6])
title('Plaids')
xlabel('Number of stimuli')
ylabel('Fraction of cells')
subplot(2,2,3)
title('Gratings')
xlabel('Number of stimuli')
ylabel('Number of cells')
ylim([0 80])
subplot(2,2,4)
title('Plaids')
xlabel('Number of stimuli')
ylabel('Number of cells')
ylim([0 80])
print(fullfile(summaryDir,[svName '_' cell2mat(area_list) '_Summary_stimOn' stimSetDur num2str(stimOn) '_' stimSetSz '_pdfs.pdf']),'-dpdf','-fillpage')

% subplot(2,3,1)
% legend(leg_str(2,:),'location','southeast')
% sgtitle('SOM cells')
% print(fullfile(summaryDir,[svName '_AllAreaSummary_stimOn' stimSetDur num2str(stimOn) '_' stimSetSz '_SOM.pdf']),'-dpdf','-fillpage')

% subplot(3,3,4)
% title('SOM+')
% legend(leg_str(2,:),'location','southeast')
% subplot(3,3,7)
% title('Pyr')
% legend(leg_str(3,:),'location','southeast')
% figure(3)
% print(fullfile(summaryDir,[svName '_AllAreaSummary_stimOn' stimSetDur num2str(stimOn) '_' stimSetSz '_ZcZpScatter.pdf']),'-dpdf','-fillpage')
%     


    
