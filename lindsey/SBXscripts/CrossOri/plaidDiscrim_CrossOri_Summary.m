close all; clear all; clc;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'PlaidDiscrimSummary');

ds = ['plaidDiscrim_exptList'];
svName = 'plaidDiscrim';
eval(ds)
rc = behavConstsAV;
frame_rate = 30;
stimOn = 25000;
stimSz = 1000;
stimSetDur = 'LT';
stimSetSz = 'FF';
area = 'RL';
nexp = size(expt,2);

resp_ind_all = [];

red_cells_all = [];
z_all = [];

stim_resp_dir_all = [];
resp_stim_short_all = [];
resp_stim_long_all = [];
resp_dec_short_all = [];
resp_dec_long_all = [];
Zc_all = [];
Zp_all = [];
corr_stim_all = [];
corr_dec_all = [];
driver_all = [];

C_stim_all = [];
P_stim_all = [];
C_dec_all = [];
P_dec_all = [];

mouse_list = [];

totTrials = zeros(nexp,2);
quant = zeros(nexp,2);

totCells = zeros(nexp,1);
for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    layer = expt(iexp).img_loc{2};
    driver = expt(iexp).driver;
    ImgFolder = expt(iexp).runs;

    nrun = 1;
    run_str = catRunName(ImgFolder, nrun);

    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
    tDoFB = celleqel2mat_padded(input.tDoFeedbackMotion);
    nTrials = length(tDoFB);
    if strcmp(expt(iexp).img_loc(1), area) & (strcmp(stimSetDur,'GT') & input.stimOnTimeMs > stimOn & sum(tDoFB)<nTrials || strcmp(stimSetDur,'LT')...
            & input.stimOnTimeMs <= stimOn) & sum(tDoFB)<nTrials & input.gratingMaxDiameterDeg == stimSz
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirAnalysis.mat']))
        if input.gratingSpatialFreqCPD == 0.1
        fprintf([mouse ' ' date '\n'])
        mouse_list = strvcat(mouse_list, mouse);

        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirAnalysis.mat']))
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_GLM.mat']))
        
        red_cells = find(mask_label);

        nCells = size(stim_resp_dir,1);
        totCells(iexp,:) = nCells;

        resp_ind_all = [resp_ind_all; resp_ind+sum(totCells(1:iexp-1,:),1)];
        red_cells_all = [red_cells_all red_cells+sum(totCells(1:iexp-1,:),1)];
        z_all = [z_all expt(iexp).z.*ones(1,nCells)];
        driver_all = [driver_all repmat(mat2cell(expt(iexp).driver,1,3),ones(1,nCells))];
        stim_resp_dir_all = cat(1, stim_resp_dir_all, stim_resp_dir(:,1:8,:,1:2));

        Zc_all = [Zc_all Zc];
        Zp_all = [Zp_all Zp];
        
        C_stim_all = [C_stim_all C_stim];
        P_stim_all = [P_stim_all P_stim];
        C_dec_all = [C_dec_all C_dec];
        P_dec_all = [P_dec_all P_dec];
        
        tDecisionTime = celleqel2mat_padded(input.tDecisionTimeMs);
        IgnoreIx = strcmp(input.trialOutcomeCell,'ignore');
        quant(iexp,:) = quantile(tDecisionTime,[0.2 0.8]);
        ind_short = find(tDecisionTime<quant(iexp,1) & ~tDoFB & ~IgnoreIx & ~isnan(squeeze(data_dec_dfof(1,1,:))'));
        ind_long = find(tDecisionTime>quant(iexp,2) & ~tDoFB & ~IgnoreIx & ~isnan(squeeze(data_dec_dfof(1,1,:))'));
        totTrials(iexp,1) = length(ind_short);
        totTrials(iexp,2) = length(ind_long);
        resp_stim_short = mean(data_stim_dfof(:,:,ind_short),3);
        resp_stim_long = mean(data_stim_dfof(:,:,ind_long),3);
        resp_stim_short_all = [resp_stim_short_all resp_stim_short];
        resp_stim_long_all = [resp_stim_long_all resp_stim_long];
        resp_dec_short = mean(data_dec_dfof(:,:,ind_short),3);
        resp_dec_long = mean(data_dec_dfof(:,:,ind_long),3);
        resp_dec_short_all = [resp_dec_short_all resp_dec_short];
        resp_dec_long_all = [resp_dec_long_all resp_dec_long];
        corr_stim = zeros(1,nCells);
        corr_dec = zeros(1,nCells);
        for iC = 1:nCells
            corr_stim(1,iC) = triu2vec(corrcoef(resp_stim_short(:,iC),resp_stim_long(:,iC)));
            corr_dec(1,iC) = triu2vec(corrcoef(resp_dec_short(:,iC),resp_dec_long(:,iC)));
        end
        corr_stim_all = [corr_stim_all corr_stim];
        corr_dec_all = [corr_dec_all corr_dec];
        end
    end
end

save(fullfile(summaryDir,[svName '_Summary_' area '_stimOn' stimSetDur num2str(stimOn) '_' stimSetSz '.mat']),'totCells','driver_all','z_all','red_cells_all','corr_stim_all','corr_dec_all','resp_dec_short_all','resp_dec_long_all','resp_stim_short_all','resp_stim_long_all','resp_ind_all','Zp_all','Zc_all','stim_resp_dir_all','C_stim_all','P_stim_all','C_dec_all','P_dec_all', 'mouse_list')

%%

resp_ind_all_pyr = setdiff(resp_ind_all, red_cells_all);
resp_ind_all_red = intersect(resp_ind_all,red_cells_all);
figure;
subplot(2,2,1)
scatter(Zc_all(resp_ind_all_pyr),Zp_all(resp_ind_all_pyr))
hold on
scatter(Zc_all(resp_ind_all_red),Zp_all(resp_ind_all_red))
xlabel('Zc')
ylabel('Zp')
title('')
ylim([-8 8])
xlim([-8 8])
plotZcZpBorders
subplot(2,2,2)
histogram(Zc_all(resp_ind_all_pyr)-Zp_all(resp_ind_all_pyr),[-10:1:10])
hold on
histogram(Zc_all(resp_ind_all_red)-Zp_all(resp_ind_all_red),[-10:1:10])
xlabel('Zc-Zp')
subplot(2,2,3)
cdfplot(Zc_all(resp_ind_all_pyr))
hold on
cdfplot(Zc_all(resp_ind_all_red))
legend({'Pyr','SOM'},'location','northwest')
xlabel('Zc')
xlim([-10 10])
title('')
subplot(2,2,4)
cdfplot(Zp_all(resp_ind_all_pyr))
hold on
cdfplot(Zp_all(resp_ind_all_red))
xlabel('Zp')
xlim([-10 10])
title('')
suptitle({[area '- n = ' num2str(size(mouse_list,1)) ' expts; ' num2str(size(unique(mouse_list,'rows'),1)) ' mice'], ['All responsive cells- n = ' num2str(length(resp_ind_all))]})
print(fullfile(summaryDir,[svName '_Summary_' area '_stimOn' stimSetDur num2str(stimOn) '_' stimSetSz '_ZpZc.pdf']),'-dpdf')

figure;
subplot(2,2,1)
scatter(corr_stim_all(resp_ind_all_pyr),Zc_all(resp_ind_all_pyr)-Zp_all(resp_ind_all_pyr))
c_pyr = triu2vec(corrcoef(corr_stim_all(resp_ind_all_pyr),Zc_all(resp_ind_all_pyr)-Zp_all(resp_ind_all_pyr)));
hold on
scatter(corr_stim_all(resp_ind_all_red),Zc_all(resp_ind_all_red)-Zp_all(resp_ind_all_red))
c_som = triu2vec(corrcoef(corr_stim_all(resp_ind_all_red),Zc_all(resp_ind_all_red)-Zp_all(resp_ind_all_red)));
xlabel('Correlation (short v long RTs)')
ylabel('Zc-Zp')
ylim([-10 10])
xlim([-1 1])
legend({['Pyr- ' num2str(chop(c_pyr,2))],['SOM- ' num2str(chop(c_som,2))]},'location','northwest')
title('Stim response')
subplot(2,2,2)
scatter(corr_dec_all(resp_ind_all_pyr),Zc_all(resp_ind_all_pyr)-Zp_all(resp_ind_all_pyr));
c_pyr = triu2vec(corrcoef(corr_dec_all(resp_ind_all_pyr),Zc_all(resp_ind_all_pyr)-Zp_all(resp_ind_all_pyr)));
hold on
scatter(corr_dec_all(resp_ind_all_red),Zc_all(resp_ind_all_red)-Zp_all(resp_ind_all_red))
c_som = triu2vec(corrcoef(corr_dec_all(resp_ind_all_red),Zc_all(resp_ind_all_red)-Zp_all(resp_ind_all_red)));
xlabel('Correlation (short v long RTs)')
ylabel('Zc-Zp')
ylim([-10 10])
xlim([-1 1])
legend({['Pyr- ' num2str(chop(c_pyr,2))],['SOM- ' num2str(chop(c_som,2))]},'location','northwest')
title('Decision response')
subplot(2,2,3)
high_corr_ind = find(corr_stim_all>0.8);
scatter(Zc_all(intersect(high_corr_ind,resp_ind_all_pyr)),Zp_all(intersect(high_corr_ind,resp_ind_all_pyr)))
hold on
scatter(Zc_all(intersect(high_corr_ind,resp_ind_all_red)),Zp_all(intersect(high_corr_ind,resp_ind_all_red)))
xlabel('Zc')
xlabel('Zp')
title('')
ylim([-8 8])
xlim([-8 8])
plotZcZpBorders
subplot(2,2,4)
histogram(Zc_all(intersect(high_corr_ind,resp_ind_all_pyr))-Zp_all(intersect(high_corr_ind,resp_ind_all_pyr)),[-10:1:10])
hold on
histogram(Zc_all(intersect(high_corr_ind,resp_ind_all_red))-Zp_all(intersect(high_corr_ind,resp_ind_all_red)),[-10:1:10])
xlabel('Zc-Zp')
suptitle({[area '- n = ' num2str(size(mouse_list,1)) ' expts; ' num2str(size(unique(mouse_list,'rows'),1)) ' mice'], ['All responsive cells- n = ' num2str(length(resp_ind_all))]})
print(fullfile(summaryDir,[svName '_Summary_' area '_stimOn' stimSetDur num2str(stimOn) '_' stimSetSz '_CorrByRT.pdf']),'-dpdf')


figure;
Zc_use = intersect(find(Zc_all>1.28),find(Zc_all-Zp_all>1.28));
Zp_use = intersect(find(Zp_all>1.28),find(Zp_all-Zc_all>1.28));
subplot(2,2,1)
shadedErrorBar(tt, mean(resp_stim_short_all(:,Zc_use),2),std(resp_stim_short_all(:,Zc_use),[],2)./sqrt(length(Zc_use)))
hold on
shadedErrorBar(tt, mean(resp_stim_long_all(:,Zc_use),2),std(resp_stim_long_all(:,Zc_use),[],2)./sqrt(length(Zc_use)),'lineprops','b')
title(['Zc cells- ' num2str(length(Zc_use))])
legend({'short','long'},'location','northwest')
xlabel('Time from stim (s)')
subplot(2,2,2)
shadedErrorBar(tt, mean(resp_stim_short_all(:,Zp_use),2),std(resp_stim_short_all(:,Zp_use),[],2)./sqrt(length(Zp_use)))
hold on
shadedErrorBar(tt, mean(resp_stim_long_all(:,Zp_use),2),std(resp_stim_long_all(:,Zp_use),[],2)./sqrt(length(Zp_use)),'lineprops','b')
title(['Zp cells- ' num2str(length(Zp_use))])
xlabel('Time from stim (s)')
subplot(2,2,3)
shadedErrorBar(tt, mean(resp_dec_short_all(:,Zc_use),2),std(resp_dec_short_all(:,Zc_use),[],2)./sqrt(length(Zc_use)))
hold on
shadedErrorBar(tt, mean(resp_dec_long_all(:,Zc_use),2),std(resp_dec_long_all(:,Zc_use),[],2)./sqrt(length(Zc_use)),'lineprops','b')
title(['Zc cells- ' num2str(length(Zc_use))])
xlabel('Time from decision (s)')
subplot(2,2,4)
shadedErrorBar(tt, mean(resp_dec_short_all(:,Zp_use),2),std(resp_dec_short_all(:,Zp_use),[],2)./sqrt(length(Zp_use)))
hold on
shadedErrorBar(tt, mean(resp_dec_long_all(:,Zp_use),2),std(resp_dec_long_all(:,Zp_use),[],2)./sqrt(length(Zp_use)),'lineprops','b')
title(['Zp cells- ' num2str(length(Zp_use))])
xlabel('Time from decision (s)')
suptitle({[area '- n = ' num2str(size(mouse_list,1)) ' expts; ' num2str(size(unique(mouse_list,'rows'),1)) ' mice'], ['All responsive cells- n = ' num2str(length(resp_ind_all))]})
print(fullfile(summaryDir,[svName '_Summary_' area '_stimOn' stimSetDur num2str(stimOn) '_' stimSetSz '_RespByRT.pdf']),'-dpdf')

figure;
% subplot(2,3,1)
% plot(C_stim_all(:,resp_ind_all_pyr),'k')
% hold on
% plot(C_stim_all(:,resp_ind_all_red),'r')
% predictors = ({'Int','Dir','Mask','Rew','Choice','RT'});
% set(gca,'XTick',1:6,'XTickLabel',predictors);
% xlim([0 7])
% ylim([-0.75 0.75])
% ylabel('Weight')
% title('Stimulus')
subplot(2,3,2)
errorbar(mean(C_stim_all(:,resp_ind_all_pyr),2),std(C_stim_all(:,resp_ind_all_pyr),[],2)./sqrt(length(resp_ind_all_pyr)),'k')
hold on
errorbar(mean(C_stim_all(:,resp_ind_all_red),2),std(C_stim_all(:,resp_ind_all_red),[],2)./sqrt(length(resp_ind_all_red)),'r')
set(gca,'XTick',1:6,'XTickLabel',predictors);
xlim([0 7])
ylim([-0.3 0.3])
ylabel('Weight')
cell_mat = nan(1,size(C_stim_all,2));
cell_mat(resp_ind_all_pyr) = 0;
cell_mat(resp_ind_all_red) = 1;
cell_mat = reshape(repmat(cell_mat,[size(C_stim_all,1) 1]), [size(C_stim_all,1)*size(C_stim_all,2) 1]);
pred_mat = reshape(repmat([1:size(C_stim_all,1)]',[1 size(C_stim_all,2)]), [size(C_stim_all,1)*size(C_stim_all,2) 1]);
pred_mat(isnan(cell_mat)) = [];
resp_stim_mat = reshape(C_stim_all,[size(C_stim_all,1)*size(C_stim_all,2) 1]);
resp_stim_mat(isnan(cell_mat)) = [];
resp_dec_mat = reshape(C_dec_all,[size(C_dec_all,1)*size(C_dec_all,2) 1]);
resp_dec_mat(isnan(cell_mat)) = [];
cell_mat(isnan(cell_mat)) = [];

[p_stim tbl_stim stats_stim] = anovan(resp_stim_mat,[pred_mat cell_mat]);
title(['Cell type- p=' num2str(chop(p_stim(2),2)) '; Pred- p=' num2str(chop(p_stim(1),2)) ])

subplot(2,3,3)
P_stim_dig = P_stim_all<0.05;
plot(sum(P_stim_dig(:,resp_ind_all_pyr),2)./length(resp_ind_all_pyr),'k')
hold on
plot(sum(P_stim_dig(:,resp_ind_all_red),2)./length(resp_ind_all_red),'r')
ylabel('Fraction significant')
set(gca,'XTick',1:6,'XTickLabel',predictors);
ylim([0 1])
xlim([0 7])
ylabel('Weight')
% subplot(2,3,4)
% plot(C_dec_all(:,resp_ind_all_pyr),'k')
% hold on
% plot(C_dec_all(:,resp_ind_all_red),'r')
% set(gca,'XTick',1:6,'XTickLabel',predictors);
% xlim([0 7])
% ylim([-0.75 0.75])
% ylabel('Weight')
% title('Decision')
subplot(2,3,5)
errorbar(mean(C_dec_all(:,resp_ind_all_pyr),2),std(C_dec_all(:,resp_ind_all_pyr),[],2)./sqrt(length(resp_ind_all_pyr)),'k')
hold on
errorbar(mean(C_dec_all(:,resp_ind_all_red),2),std(C_dec_all(:,resp_ind_all_red),[],2)./sqrt(length(resp_ind_all_red)),'r')
set(gca,'XTick',1:6,'XTickLabel',predictors);
xlim([0 7])
ylim([-0.3 0.3])
ylabel('Weight')
[p_dec tbl_dec stats_dec] = anovan(resp_dec_mat,[pred_mat cell_mat]);
title(['Cell type- p=' num2str(chop(p_dec(2),2)) '; Pred- p=' num2str(chop(p_dec(1),2)) ])
subplot(2,3,6)
P_dec_dig = P_dec_all<0.05;
plot(sum(P_dec_dig(:,resp_ind_all_pyr),2)./length(resp_ind_all_pyr),'k')
hold on
plot(sum(P_dec_dig(:,resp_ind_all_red),2)./length(resp_ind_all_red),'r')
ylabel('Fraction significant')
set(gca,'XTick',1:6,'XTickLabel',predictors);
ylim([0 1])
xlim([0 7])

print(fullfile(summaryDir,[svName '_Summary_' area '_stimOn' stimSetDur num2str(stimOn) '_' stimSetSz '_GLM.pdf']), '-dpdf','-bestfit')

figure;
start = 1;
for i = 1:5
    for ii = 1:5
        if i ~= ii
            subplot(5,5,start)
            scatter(C_stim_all(i+1,resp_ind_all_pyr), C_stim_all(ii+1,resp_ind_all_pyr))
            hold on
            scatter(C_stim_all(i+1,resp_ind_all_red), C_stim_all(ii+1,resp_ind_all_red))
            xlabel(predictors(i+1))
            ylabel(predictors(ii+1))
            xlim([-1 1])
            ylim([-1 1])
            axis square
        end
        start = start+1;
    end
end

figure;
ind_p = cell(1,5);
for i = 1:5
    subplot(3,2,i)
    ind_p{i} = find(P_stim_all(i+1,:)<0.05);
    ind_pyr = intersect(ind_p{i},resp_ind_all_pyr);
    ind_red = intersect(ind_p{i},resp_ind_all_red);
    scatter(Zc_all(ind_pyr),Zp_all(ind_pyr))
    hold on
    scatter(Zc_all(ind_red),Zp_all(ind_red))
    xlabel('Zc')
    ylabel('Zp')
    title('')
    ylim([-8 8])
    xlim([-8 8])
    plotZcZpBorders
    title(predictors(i+1))
    axis square
end

ind_dir = setdiff(ind_p{1},[ind_p{3} ind_p{4} ind_p{5}]);
subplot(3,2,6)
scatter(Zc_all(intersect(ind_dir,resp_ind_all_pyr)),Zp_all(intersect(ind_dir,resp_ind_all_pyr)))
hold on
scatter(Zc_all(intersect(ind_dir,resp_ind_all_red)),Zp_all(intersect(ind_dir,resp_ind_all_red)))
xlabel('Zc')
ylabel('Zp')
title('')
ylim([-8 8])
xlim([-8 8])
plotZcZpBorders
axis square
title('Dir only')
print(fullfile(summaryDir,[svName '_Summary_' area '_stimOn' stimSetDur num2str(stimOn) '_' stimSetSz '_GLM_ZcZp.pdf']), '-dpdf','-bestfit')

figure;
for i = 1:5
    subplot(3,2,i)
    scatter(C_stim_all(i+1,resp_ind_all_pyr),Zc_all(resp_ind_all_pyr)-Zp_all(resp_ind_all_pyr))
    hold on
    scatter(C_stim_all(i+1,resp_ind_all_red),Zc_all(resp_ind_all_red)-Zp_all(resp_ind_all_red))
    xlabel(predictors{i+1})
    ylabel('Zc-Zp')
    ylim([-10 10])
    xlim([-1 1])
    axis square
    title(num2str(chop(triu2vec(corrcoef(C_stim_all(i+1,resp_ind_all),Zc_all(resp_ind_all)-Zp_all(resp_ind_all))),2)))
end
print(fullfile(summaryDir,[svName '_Summary_' area '_stimOn' stimSetDur num2str(stimOn) '_' stimSetSz '_GLM_ZcminusZp.pdf']), '-dpdf','-bestfit')

Zp_use_pyr = intersect(resp_ind_all_pyr, intersect(find(Zp_all>1.28), find(Zp_all-Zc_all>1.28)));
Zc_use_pyr = intersect(resp_ind_all_pyr, intersect(find(Zc_all>1.28), find(Zc_all-Zp_all>1.28)));
Zp_use_red = intersect(resp_ind_all_red, intersect(find(Zp_all>1.28), find(Zp_all-Zc_all>1.28)));
Zc_use_red = intersect(resp_ind_all_red, intersect(find(Zc_all>1.28), find(Zc_all-Zp_all>1.28)));
figure;
subplot(2,2,1)
errorbar(mean(C_stim_all(:,Zc_use_pyr),2),std(C_stim_all(:,Zc_use_pyr),[],2)./sqrt(length(Zc_use_pyr)))
hold on
errorbar(mean(C_stim_all(:,Zc_use_red),2),std(C_stim_all(:,Zc_use_red),[],2)./sqrt(length(Zc_use_red)))
set(gca,'XTick',1:6,'XTickLabel',predictors);
xlim([0 7])
ylim([-0.3 0.3])
ylabel('Weight')
title(['Zc cells- ' num2str(length(Zc_use))])
subplot(2,2,3)
P_stim_dig = P_stim_all<0.05;
plot(sum(P_stim_dig(:,Zc_use_pyr),2)./length(Zc_use_pyr))
hold on
plot(sum(P_stim_dig(:,Zc_use_red),2)./length(Zc_use_red))
ylabel('Fraction significant')
set(gca,'XTick',1:6,'XTickLabel',predictors);
ylim([0 1])
xlim([0 7])
subplot(2,2,2)
errorbar(mean(C_stim_all(:,Zp_use_pyr),2),std(C_stim_all(:,Zp_use_pyr),[],2)./sqrt(length(Zp_use_pyr)))
hold on
errorbar(mean(C_stim_all(:,Zp_use_red),2),std(C_stim_all(:,Zp_use_red),[],2)./sqrt(length(Zp_use_red)))
set(gca,'XTick',1:6,'XTickLabel',predictors);
xlim([0 7])
ylim([-0.3 0.3])
ylabel('Weight')
title(['Zp cells- ' num2str(length(Zp_use))])
subplot(2,2,4)
P_stim_dig = P_stim_all<0.05;
plot(sum(P_stim_dig(:,Zp_use_pyr),2)./length(Zp_use_pyr))
hold on
plot(sum(P_stim_dig(:,Zp_use_red),2)./length(Zp_use_red))
ylabel('Fraction significant')
set(gca,'XTick',1:6,'XTickLabel',predictors);
ylim([0 1])
xlim([0 7])
print(fullfile(summaryDir,[svName '_Summary_' area '_stimOn' stimSetDur num2str(stimOn) '_' stimSetSz '_GLM_ZcVsZp.pdf']), '-dpdf','-bestfit')


[idx c] = kmeans(C_stim_all',6);
figure;
for i = 1:6
    ind = find(idx==i);
    subplot(2,3,i)
    errorbar(1:6, mean(C_stim_all(:,ind),2),std(C_stim_all(:,ind),[],2)./sqrt(length(ind)))
    ylim([-0.5 0.5])
    set(gca,'XTick',1:6,'XTickLabel',predictors);
    xlim([0 7])
    title(num2str(length(ind)))
end

C_stim_abs = C_stim_all;
C_stim_abs([2 4 5],:) = abs(C_stim_all([2 4 5],:));
figure;
[idx c] = kmeans(C_stim_abs',4);
for i = 1:4
    ind = find(idx==i);
    subplot(2,3,i)
    errorbar(1:6, mean(C_stim_abs(:,ind),2),std(C_stim_abs(:,ind),[],2)./sqrt(length(ind)))
    ylim([-0.5 0.5])
    set(gca,'XTick',1:6,'XTickLabel',predictors);
    xlim([0 7])
    title(num2str(length(ind)))
end