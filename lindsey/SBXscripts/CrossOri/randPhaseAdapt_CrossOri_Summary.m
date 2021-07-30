close all; clear all; clc;
doRedChannel = 0;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');

ds = ['CrossOriRandPhaseAdapt_ExptList'];
svName = 'randPhase4Adapt';
eval(ds)
driver = 'SLC';
area = 'V1';
max_dist = 5;
rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);

sigCells = [];
fractSigCells = [];
respCells = [];
resp_ind_all = [];

p_anova_noadapt_all = [];
b_noadapt_all = [];
amp_noadapt_all = [];
pha_noadapt_all = [];
per_noadapt_all = [];
sse_noadapt_all = [];
R_square_noadapt_all = [];
yfit_noadapt_all = [];
resp_avg_noadapt_all = [];

p_anova_singadapt_all = [];
b_singadapt_all = [];
amp_singadapt_all = [];
pha_singadapt_all = [];
per_singadapt_all = [];
sse_singadapt_all = [];
R_square_singadapt_all = [];
yfit_singadapt_all = [];
resp_avg_singadapt_all = [];

p_anova_shuf_noadapt_all = [];
b_shuf_noadapt_all = [];
amp_shuf_noadapt_all = [];
pha_shuf_noadapt_all = [];
per_shuf_noadapt_all = [];
sse_shuf_noadapt_all = [];
R_square_shuf_noadapt_all = [];
yfit_shuf_noadapt_all = [];

p_anova_shuf_singadapt_all = [];
b_shuf_singadapt_all = [];
amp_shuf_singadapt_all = [];
pha_shuf_singadapt_all = [];
per_shuf_singadapt_all = [];
sse_shuf_singadapt_all = [];
R_square_shuf_singadapt_all = [];
yfit_shuf_singadapt_all = [];

plaidSI_noadapt_all = [];
testPI_noadapt_all = [];
plaidSI_singadapt_all = [];
testPI_singadapt_all = [];
preftest_ind_all = [];
prefmask_ind_all =[];
test_resp_noadapt_all = [];
mask_resp_noadapt_all = [];
test_resp_singadapt_all = [];
mask_resp_singadapt_all = [];

trN_noadapt_all = [];
trN_singadapt_all = [];

mouse_list = [];

totCells = zeros(nexp,1);
for iexp = 1:nexp
    if strcmp(expt(iexp).driver,driver)
        mouse = expt(iexp).mouse;
        mouse_list = strvcat(mouse_list, mouse);
        date = expt(iexp).date;
        area = expt(iexp).img_loc{1};
        if isfield(expt,'copFolder')
            ImgFolder = expt(iexp).copFolder;
        else
            ImgFolder = expt(iexp).coFolder;
        end
        nrun = length(ImgFolder);
        run_str = catRunName(cell2mat(ImgFolder), nrun);

        fprintf([mouse ' ' date '\n'])

        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_' num2str(max_dist) '_phaseFits.mat']))

        nCells = size(resp_cell_noadapt{end,end,end},1);
        totCells(iexp,:) = nCells;

        p_anova_noadapt(find(p_anova_noadapt==0)) = NaN;
        p_anova_noadapt_all = [p_anova_noadapt_all;  p_anova_noadapt];
        b_noadapt_all = [b_noadapt_all; b_hat_noadapt];
        amp_noadapt_all = [amp_noadapt_all; amp_hat_noadapt];
        pha_noadapt_all = [pha_noadapt_all; pha_hat_noadapt];
        per_noadapt_all = [per_noadapt_all; per_hat_noadapt];
        R_square_noadapt_all = [R_square_noadapt_all; R_square_noadapt];
        sse_noadapt_all = [sse_noadapt_all; sse_noadapt];
        resp_avg_noadapt_all = cat(1,resp_avg_noadapt_all, resp_avg_noadapt);
        
        p_anova_singadapt(find(p_anova_singadapt==0)) = NaN;
        p_anova_singadapt_all = [p_anova_singadapt_all;  p_anova_singadapt];
        b_singadapt_all = [b_singadapt_all; b_hat_singadapt];
        amp_singadapt_all = [amp_singadapt_all; amp_hat_singadapt];
        pha_singadapt_all = [pha_singadapt_all; pha_hat_singadapt];
        per_singadapt_all = [per_singadapt_all; per_hat_singadapt];
        R_square_singadapt_all = [R_square_singadapt_all; R_square_singadapt];
        sse_singadapt_all = [sse_singadapt_all; sse_singadapt];
        resp_avg_singadapt_all = cat(1,resp_avg_singadapt_all, resp_avg_singadapt);
        
        p_anova_shuf_noadapt(find(p_anova_shuf_noadapt==0)) = NaN;
        p_anova_shuf_noadapt_all = [p_anova_shuf_noadapt_all;  p_anova_shuf_noadapt];
        b_shuf_noadapt_all = [b_shuf_noadapt_all; b_hat_shuf_noadapt];
        amp_shuf_noadapt_all = [amp_shuf_noadapt_all; amp_hat_shuf_noadapt];
        pha_shuf_noadapt_all = [pha_shuf_noadapt_all; pha_hat_shuf_noadapt];
        per_shuf_noadapt_all = [per_shuf_noadapt_all; per_hat_shuf_noadapt];
        R_square_shuf_noadapt_all = [R_square_shuf_noadapt_all; R_square_shuf_noadapt];
        sse_shuf_noadapt_all = [sse_shuf_noadapt_all; sse_shuf_noadapt];
        
        p_anova_shuf_singadapt(find(p_anova_shuf_singadapt==0)) = NaN;
        p_anova_shuf_singadapt_all = [p_anova_shuf_singadapt_all;  p_anova_shuf_singadapt];
        b_shuf_singadapt_all = [b_shuf_singadapt_all; b_hat_shuf_singadapt];
        amp_shuf_singadapt_all = [amp_shuf_singadapt_all; amp_hat_shuf_singadapt];
        pha_shuf_singadapt_all = [pha_shuf_singadapt_all; pha_hat_shuf_singadapt];
        per_shuf_singadapt_all = [per_shuf_singadapt_all; per_hat_shuf_singadapt];
        R_square_shuf_singadapt_all = [R_square_shuf_singadapt_all; R_square_shuf_singadapt];
        sse_shuf_singadapt_all = [sse_shuf_singadapt_all; sse_shuf_singadapt];

        yfit_noadapt_all = cat(1,yfit_noadapt_all, yfit_noadapt);
        yfit_shuf_noadapt_all = cat(1,yfit_shuf_noadapt_all, yfit_shuf_noadapt);
        yfit_singadapt_all = cat(1,yfit_singadapt_all, yfit_singadapt);
        yfit_shuf_singadapt_all = cat(1,yfit_shuf_singadapt_all, yfit_shuf_singadapt);

        resp_ind_all = [resp_ind_all; resp_ind+sum(totCells(1:iexp-1,:),1)];
        preftest_ind_all = [preftest_ind_all; preftest_ind+sum(totCells(1:iexp-1,:),1)];
        prefmask_ind_all = [prefmask_ind_all; prefmask_ind+sum(totCells(1:iexp-1,:),1)];
        
        trN_noadapt_all = [trN_noadapt_all; trN_noadapt];
        trN_singadapt_all = [trN_singadapt_all; trN_singadapt];
        
        mask_resp_noadapt = nanmean(resp_cell_noadapt{end,1,1},2);
        test_resp_noadapt = nanmean(resp_cell_noadapt{1,end,1},2);
        mask_resp_noadapt(find(mask_resp_noadapt<0)) = 0;
        test_resp_noadapt(find(test_resp_noadapt<0)) = 0;
        mask_resp_noadapt_all = [mask_resp_noadapt_all; mask_resp_noadapt];
        test_resp_noadapt_all = [test_resp_noadapt_all; test_resp_noadapt];

        testPI_noadapt_all = [testPI_noadapt_all; ((test_resp_noadapt-mask_resp_noadapt) ./ (mask_resp_noadapt+test_resp_noadapt))];
        
        mask_resp_singadapt = nanmean(resp_cell_singadapt{end,1,1},2);
        test_resp_singadapt = nanmean(resp_cell_singadapt{1,end,1},2);
        mask_resp_singadapt(find(mask_resp_singadapt<0)) = 0;
        test_resp_singadapt(find(test_resp_singadapt<0)) = 0;
        mask_resp_singadapt_all = [mask_resp_singadapt_all; mask_resp_singadapt];
        test_resp_singadapt_all = [test_resp_singadapt_all; test_resp_singadapt];
        
        testPI_singadapt_all = [testPI_singadapt_all; ((test_resp_singadapt-mask_resp_singadapt) ./ (mask_resp_singadapt+test_resp_singadapt))];
        
        plaid_resp_noadapt = zeros(nCells,nMaskPhas);
        plaid_resp_singadapt = zeros(nCells,nMaskPhas);
        plaidSI_noadapt = zeros(nCells,nMaskPhas);
        plaidSI_singadapt = zeros(nCells,nMaskPhas);
        
        for ip = 1:nMaskPhas
            plaid_resp_noadapt(:,ip) = nanmean(resp_cell_noadapt{end,end,ip},2);
            plaid_resp_noadapt(find(plaid_resp_noadapt(:,ip)<0),ip) = 0;
            plaid_resp_singadapt(:,ip) = nanmean(resp_cell_singadapt{end,end,ip},2);
            plaid_resp_singadapt(find(plaid_resp_singadapt(:,ip)<0),ip) = 0;
            plaidSI_noadapt(:,ip) = (plaid_resp_noadapt(:,ip)-(mask_resp_noadapt+test_resp_noadapt)) ./ (plaid_resp_noadapt(:,ip) + mask_resp_noadapt + test_resp_noadapt);
            plaidSI_singadapt(:,ip) = (plaid_resp_singadapt(:,ip)-(mask_resp_singadapt+test_resp_singadapt)) ./ (plaid_resp_singadapt(:,ip) + mask_resp_singadapt + test_resp_singadapt);
        end
        plaidSI_noadapt_all = [plaidSI_noadapt_all; plaidSI_noadapt];
        plaidSI_singadapt_all = [plaidSI_singadapt_all; plaidSI_singadapt]; 
    end
end
save(fullfile(summaryDir,[svName '_Summary_' area '_' driver '.mat']),'p_anova_noadapt_all','b_noadapt_all','amp_noadapt_all','R_square_noadapt_all','p_anova_singadapt_all','b_singadapt_all','amp_singadapt_all','R_square_singadapt_all','p_anova_shuf_noadapt_all','b_shuf_noadapt_all','amp_shuf_noadapt_all','R_square_shuf_noadapt_all','p_anova_shuf_singadapt_all','b_shuf_singadapt_all','amp_shuf_singadapt_all','R_square_shuf_singadapt_all','yfit_noadapt_all','yfit_singadapt_all','yfit_shuf_noadapt_all','yfit_shuf_singadapt_all','plaidSI_noadapt_all','testPI_noadapt_all','plaidSI_singadapt_all','testPI_singadapt_all','resp_ind_all','preftest_ind_all','prefmask_ind_all','trN_noadapt_all', 'trN_singadapt_all','mouse_list')
%%
nopref_ind_all = setdiff(resp_ind_all,[preftest_ind_all; prefmask_ind_all]);

nCells = size(resp_avg_noadapt_all,1);
resp_avg_noadapt_align = zeros(nCells,nMaskPhas);
resp_avg_singadapt_align = zeros(nCells,nMaskPhas);
npos = zeros(nCells,1);
for iC = 1:nCells
    [max_val max_ind] = max(resp_avg_noadapt_all(iC,:),[],2);
    resp_avg_noadapt_align(iC,:) = circshift(resp_avg_noadapt_all(iC,:),-max_ind+1);
    resp_avg_singadapt_align(iC,:) = circshift(resp_avg_singadapt_all(iC,:),-max_ind+1);
    npos(iC,1) = sum(plaidSI_noadapt_all(iC,:)>0);
end
figure; 
subplot(3,2,1)
errorbar(maskPhas, mean(resp_avg_noadapt_align(intersect(find(npos>0),preftest_ind_all),:),1),std(resp_avg_noadapt_align(intersect(find(npos>0),preftest_ind_all),:),1)./sqrt(length(intersect(find(npos>0),preftest_ind_all))))
hold on 
errorbar(maskPhas, mean(resp_avg_singadapt_align(intersect(find(npos>0),preftest_ind_all),:),1),std(resp_avg_singadapt_align(intersect(find(npos>0),preftest_ind_all),:),1)./sqrt(length(intersect(find(npos>0),preftest_ind_all))))
xlabel('Phase')
ylim([0 0.6])
p_preftestF_noadapt = anova1(resp_avg_noadapt_align(intersect(find(npos>0),preftest_ind_all),2:end),[],'off');
p_preftestF_singadapt = anova1(resp_avg_singadapt_align(intersect(find(npos>0),preftest_ind_all),2:end),[],'off');
title(['Test pref F- p=' num2str(chop(p_preftestF_noadapt,2)) ' ' num2str(chop(p_preftestF_singadapt,2))])
subplot(3,2,2)
errorbar(maskPhas, mean(resp_avg_noadapt_align(intersect(find(npos==0),preftest_ind_all),:),1),std(resp_avg_noadapt_align(intersect(find(npos==0),preftest_ind_all),:),1)./sqrt(length(intersect(find(npos==0),preftest_ind_all))))
hold on 
errorbar(maskPhas, mean(resp_avg_singadapt_align(intersect(find(npos==0),preftest_ind_all),:),1),std(resp_avg_singadapt_align(intersect(find(npos==0),preftest_ind_all),:),1)./sqrt(length(intersect(find(npos==0),preftest_ind_all))))
xlabel('Phase')
ylim([0 0.6])
p_preftestS_noadapt = anova1(resp_avg_noadapt_align(intersect(find(npos==0),preftest_ind_all),2:end),[],'off');
p_preftestS_singadapt = anova1(resp_avg_singadapt_align(intersect(find(npos==0),preftest_ind_all),2:end),[],'off');
title(['Test pref S- p=' num2str(chop(p_preftestS_noadapt,2)) ' ' num2str(chop(p_preftestS_singadapt,2))])
subplot(3,2,3)
errorbar(maskPhas, mean(resp_avg_noadapt_align(intersect(find(npos>0),prefmask_ind_all),:),1),std(resp_avg_noadapt_align(intersect(find(npos>0),prefmask_ind_all),:),1)./sqrt(length(intersect(find(npos>0),prefmask_ind_all))))
hold on 
errorbar(maskPhas, mean(resp_avg_singadapt_align(intersect(find(npos>0),prefmask_ind_all),:),1),std(resp_avg_singadapt_align(intersect(find(npos>0),prefmask_ind_all),:),1)./sqrt(length(intersect(find(npos>0),prefmask_ind_all))))
xlabel('Phase')
ylim([0 0.6])
p_prefmaskF_noadapt = anova1(resp_avg_noadapt_align(intersect(find(npos>0),prefmask_ind_all),2:end),[],'off');
p_prefmaskF_singadapt = anova1(resp_avg_singadapt_align(intersect(find(npos>0),prefmask_ind_all),2:end),[],'off');
title(['Mask pref F- p=' num2str(chop(p_prefmaskF_noadapt,2)) ' ' num2str(chop(p_prefmaskF_singadapt,2))])
subplot(3,2,4)
errorbar(maskPhas, mean(resp_avg_noadapt_align(intersect(find(npos==0),prefmask_ind_all),:),1),std(resp_avg_noadapt_align(intersect(find(npos==0),prefmask_ind_all),:),1)./sqrt(length(intersect(find(npos==0),prefmask_ind_all))))
hold on 
errorbar(maskPhas, mean(resp_avg_singadapt_align(intersect(find(npos==0),prefmask_ind_all),:),1),std(resp_avg_singadapt_align(intersect(find(npos==0),prefmask_ind_all),:),1)./sqrt(length(intersect(find(npos==0),prefmask_ind_all))))
xlabel('Phase')
ylim([0 0.6])
p_prefmaskS_noadapt = anova1(resp_avg_noadapt_align(intersect(find(npos==0),prefmask_ind_all),2:end),[],'off');
p_prefmaskS_singadapt = anova1(resp_avg_singadapt_align(intersect(find(npos==0),prefmask_ind_all),2:end),[],'off');
title(['Mask pref S- p=' num2str(chop(p_prefmaskS_noadapt,2)) ' ' num2str(chop(p_prefmaskS_singadapt,2))])
subplot(3,2,5)
errorbar(maskPhas, mean(resp_avg_noadapt_align(intersect(find(npos>0),nopref_ind_all),:),1),std(resp_avg_noadapt_align(intersect(find(npos>0),nopref_ind_all),:),1)./sqrt(length(intersect(find(npos>0),nopref_ind_all))))
hold on 
errorbar(maskPhas, mean(resp_avg_singadapt_align(intersect(find(npos>0),nopref_ind_all),:),1),std(resp_avg_singadapt_align(intersect(find(npos>0),nopref_ind_all),:),1)./sqrt(length(intersect(find(npos>0),nopref_ind_all))))
xlabel('Phase')
ylim([0 0.6])
p_noprefF_noadapt = anova1(resp_avg_noadapt_align(intersect(find(npos>0),nopref_ind_all),2:end),[],'off');
p_noprefF_singadapt = anova1(resp_avg_singadapt_align(intersect(find(npos>0),nopref_ind_all),2:end),[],'off');
title(['No pref F- p=' num2str(chop(p_noprefF_noadapt,2)) ' ' num2str(chop(p_noprefF_singadapt,2))])
subplot(3,2,6)
errorbar(maskPhas, mean(resp_avg_noadapt_align(intersect(find(npos==0),nopref_ind_all),:),1),std(resp_avg_noadapt_align(intersect(find(npos==0),nopref_ind_all),:),1)./sqrt(length(intersect(find(npos==0),nopref_ind_all))))
hold on 
errorbar(maskPhas, mean(resp_avg_singadapt_align(intersect(find(npos==0),nopref_ind_all),:),1),std(resp_avg_singadapt_align(intersect(find(npos==0),nopref_ind_all),:),1)./sqrt(length(intersect(find(npos==0),nopref_ind_all))))
xlabel('Phase')
ylim([0 0.6])
p_noprefS_noadapt = anova1(resp_avg_noadapt_align(intersect(find(npos==0),nopref_ind_all),2:end),[],'off');
p_noprefS_singadapt = anova1(resp_avg_singadapt_align(intersect(find(npos==0),nopref_ind_all),2:end),[],'off');
title(['No pref S- p=' num2str(chop(p_noprefS_noadapt,2)) ' ' num2str(chop(p_noprefS_singadapt,2))])

figure;
subplot(2,2,1)
errorbar(mean(test_resp_noadapt_all(resp_ind_all,:),1),mean(mask_resp_noadapt_all(resp_ind_all,:),1), std(mask_resp_noadapt_all(resp_ind_all,:),[],1)./sqrt(length(resp_ind_all)),std(mask_resp_noadapt_all(resp_ind_all,:),[],1)./sqrt(length(resp_ind_all)),std(test_resp_noadapt_all(resp_ind_all,:),[],1)./sqrt(length(resp_ind_all)),std(test_resp_noadapt_all(resp_ind_all,:),[],1)./sqrt(length(resp_ind_all)),'o')
hold on
errorbar(mean(test_resp_singadapt_all(resp_ind_all,:),1),mean(mask_resp_singadapt_all(resp_ind_all,:),1), std(mask_resp_singadapt_all(resp_ind_all,:),[],1)./sqrt(length(resp_ind_all)),std(mask_resp_singadapt_all(resp_ind_all,:),[],1)./sqrt(length(resp_ind_all)),std(test_resp_singadapt_all(resp_ind_all,:),[],1)./sqrt(length(resp_ind_all)),std(test_resp_singadapt_all(resp_ind_all,:),[],1)./sqrt(length(resp_ind_all)),'o')
xlabel('Test response')
ylabel('Mask response')
xlim([0 .5])
ylim([0 .5])
title('All cells')
subplot(2,2,3)
errorbar(mean(test_resp_noadapt_all(preftest_ind_all,:),1),mean(mask_resp_noadapt_all(preftest_ind_all,:),1), std(mask_resp_noadapt_all(preftest_ind_all,:),[],1)./sqrt(length(preftest_ind_all)),std(mask_resp_noadapt_all(preftest_ind_all,:),[],1)./sqrt(length(preftest_ind_all)),std(test_resp_noadapt_all(preftest_ind_all,:),[],1)./sqrt(length(preftest_ind_all)),std(test_resp_noadapt_all(preftest_ind_all,:),[],1)./sqrt(length(preftest_ind_all)),'o')
hold on
errorbar(mean(test_resp_singadapt_all(preftest_ind_all,:),1),mean(mask_resp_singadapt_all(preftest_ind_all,:),1), std(mask_resp_singadapt_all(preftest_ind_all,:),[],1)./sqrt(length(preftest_ind_all)),std(mask_resp_singadapt_all(preftest_ind_all,:),[],1)./sqrt(length(preftest_ind_all)),std(test_resp_singadapt_all(preftest_ind_all,:),[],1)./sqrt(length(preftest_ind_all)),std(test_resp_singadapt_all(preftest_ind_all,:),[],1)./sqrt(length(preftest_ind_all)),'o')
xlabel('Test response')
ylabel('Mask response')
xlim([0 .5])
ylim([0 .5])
title('Test Pref cells')
subplot(2,2,4)
errorbar(mean(test_resp_noadapt_all(prefmask_ind_all,:),1),mean(mask_resp_noadapt_all(prefmask_ind_all,:),1), std(mask_resp_noadapt_all(prefmask_ind_all,:),[],1)./sqrt(length(prefmask_ind_all)),std(mask_resp_noadapt_all(prefmask_ind_all,:),[],1)./sqrt(length(prefmask_ind_all)),std(test_resp_noadapt_all(prefmask_ind_all,:),[],1)./sqrt(length(prefmask_ind_all)),std(test_resp_noadapt_all(prefmask_ind_all,:),[],1)./sqrt(length(prefmask_ind_all)),'o')
hold on
errorbar(mean(test_resp_singadapt_all(prefmask_ind_all,:),1),mean(mask_resp_singadapt_all(prefmask_ind_all,:),1), std(mask_resp_singadapt_all(prefmask_ind_all,:),[],1)./sqrt(length(prefmask_ind_all)),std(mask_resp_singadapt_all(prefmask_ind_all,:),[],1)./sqrt(length(prefmask_ind_all)),std(test_resp_singadapt_all(prefmask_ind_all,:),[],1)./sqrt(length(prefmask_ind_all)),std(test_resp_singadapt_all(prefmask_ind_all,:),[],1)./sqrt(length(prefmask_ind_all)),'o')
xlabel('Test response')
ylabel('Mask response')
xlim([0 .5])
ylim([0 .5])
title('Mask Pref cells')
print(fullfile(summaryDir, [svName '_TestMaskResp_AllTrialFits.pdf']),'-dpdf','-fillpage')

figure;
subplot(3,1,1)
cdfplot(testPI_noadapt_all(resp_ind_all,:))
n = sum(~isnan(testPI_noadapt_all(resp_ind_all,:)));
hold on
cdfplot(testPI_singadapt_all(resp_ind_all,:))
xlabel('Test/Mask preference')
xlim([-1 1])
legend({'Control','Adapt'},'location','southeast')
title(['All cells- n = ' num2str(n)])
subplot(3,1,2)
cdfplot(testPI_noadapt_all(preftest_ind_all,:))
n = sum(~isnan(testPI_noadapt_all(preftest_ind_all,:)));
hold on
cdfplot(testPI_singadapt_all(preftest_ind_all,:))
xlabel('Test/Mask preference')
xlim([-1 1])
legend({'Control','Adapt'},'location','southeast')
title(['Pref test- n = ' num2str(n)])
subplot(3,1,3)
cdfplot(testPI_noadapt_all(prefmask_ind_all,:))
n = sum(~isnan(testPI_noadapt_all(prefmask_ind_all,:)));
hold on
cdfplot(testPI_singadapt_all(prefmask_ind_all,:))
xlabel('Test/Mask preference')
xlim([-1 1])
legend({'Control','Adapt'},'location','southeast')
title(['Pref mask- n = ' num2str(n)])
print(fullfile(summaryDir, [svName '_TestPI_AllTrialFits.pdf']),'-dpdf','-fillpage')

figure;
subplot(3,1,1)
cdfplot(plaidSI_noadapt_all(resp_ind_all,1))
n = sum(~isnan(plaidSI_noadapt_all(resp_ind_all,1)));
hold on
cdfplot(plaidSI_singadapt_all(resp_ind_all,1))
xlabel('Suppression index')
xlim([-1 1])
legend({'Control','Adapt'},'location','southeast')
title(['All cells- n = ' num2str(n)])
subplot(3,1,2)
cdfplot(plaidSI_noadapt_all(preftest_ind_all,1))
n = sum(~isnan(plaidSI_noadapt_all(preftest_ind_all,1)));
hold on
cdfplot(plaidSI_singadapt_all(preftest_ind_all,1))
xlabel('Suppression index')
xlim([-1 1])
legend({'Control','Adapt'},'location','southeast')
title(['Pref test- n = ' num2str(n)])
subplot(3,1,3)
cdfplot(plaidSI_noadapt_all(prefmask_ind_all,1))
n = sum(~isnan(plaidSI_noadapt_all(prefmask_ind_all,1)));
hold on
cdfplot(plaidSI_singadapt_all(prefmask_ind_all,1))
xlabel('Suppression index')
xlim([-1 1])
legend({'Control','Adapt'},'location','southeast')
title(['Pref mask- n = ' num2str(n)])
print(fullfile(summaryDir, [svName '_PlaidSI_AllTrialFits.pdf']),'-dpdf','-fillpage')


figure;
subplot(2,2,1)
scatter(p_anova_noadapt_all(resp_ind_all,:), amp_noadapt_all(resp_ind_all,:))
hold on
scatter(p_anova_noadapt_all(preftest_ind_all,:), amp_noadapt_all(preftest_ind_all,:))
xlabel('anova p-value')
ylabel('Sine amplitude')
title('No adapt')
ylim([0 2])
subplot(2,2,2)
scatter(p_anova_singadapt_all(resp_ind_all,:), amp_singadapt_all(resp_ind_all,:))
hold on
scatter(p_anova_singadapt_all(preftest_ind_all,:), amp_singadapt_all(preftest_ind_all,:))
xlabel('anova p-value')
ylabel('Sine amplitude')
ylim([0 2])
title('Sing adapt')
subplot(2,2,3)
scatter(R_square_noadapt_all(resp_ind_all,:), amp_noadapt_all(resp_ind_all,:))
hold on
scatter(R_square_noadapt_all(preftest_ind_all,:), amp_noadapt_all(preftest_ind_all,:))
xlabel('R-squared')
ylabel('Sine amplitude')
title('No adapt')
ylim([0 2])
subplot(2,2,4)
scatter(R_square_singadapt_all(resp_ind_all,:), amp_singadapt_all(resp_ind_all,:))
hold on
scatter(R_square_singadapt_all(preftest_ind_all,:), amp_singadapt_all(preftest_ind_all,:))
xlabel('R-squared')
ylabel('Sine amplitude')
ylim([0 2])
title('Sing adapt')
movegui('center')


%Rsquare
figure;
subplot(3,2,1)
Rsq_temp = R_square_noadapt_all;
Rsq_temp(find(Rsq_temp<0)) = 0;
cdfplot(Rsq_temp(resp_ind_all,:))
n = sum(~isnan(Rsq_temp(resp_ind_all,:)));
hold on
Rsq_temp = R_square_shuf_noadapt_all;
Rsq_temp(find(Rsq_temp<0)) = 0;
cdfplot(Rsq_temp(resp_ind_all,:))
xlabel('Rsquared')
xlim([0 1])
legend({['Data- n = ' num2str(n)],'Shuffled'},'location','southeast')
title('No Adapt')
subplot(3,2,2)
Rsq_temp = R_square_singadapt_all;
Rsq_temp(find(Rsq_temp<0)) = 0;
cdfplot(Rsq_temp(resp_ind_all,:))
n = sum(~isnan(Rsq_temp(resp_ind_all,:)));
hold on
Rsq_temp = R_square_shuf_singadapt_all;
Rsq_temp(find(Rsq_temp<0)) = 0;
cdfplot(Rsq_temp(resp_ind_all,:))
xlabel('Rsquared')
xlim([0 1])
legend({['Data- n = ' num2str(n)],'Shuffled'},'location','southeast')
title('Sing Adapt')
subplot(3,2,3)
Rsq_temp = R_square_noadapt_all;
Rsq_temp(find(Rsq_temp<0)) = 0;
cdfplot(Rsq_temp(preftest_ind_all,:))
n = sum(~isnan(Rsq_temp(preftest_ind_all,:)));
hold on
Rsq_temp = R_square_shuf_noadapt_all;
Rsq_temp(find(Rsq_temp<0)) = 0;
cdfplot(Rsq_temp(preftest_ind_all,:))
xlabel('Rsquared')
xlim([0 1])
legend({['Data- n = ' num2str(n)],'Shuffled'},'location','southeast')
title('No Adapt- Pref test')
subplot(3,2,4)
Rsq_temp = R_square_singadapt_all;
Rsq_temp(find(Rsq_temp<0)) = 0;
cdfplot(Rsq_temp(preftest_ind_all,:))
n = sum(~isnan(Rsq_temp(preftest_ind_all,:)));
hold on
Rsq_temp = R_square_shuf_singadapt_all;
Rsq_temp(find(Rsq_temp<0)) = 0;
cdfplot(Rsq_temp(preftest_ind_all,:))
xlabel('Rsquared')
xlim([0 1])
legend({['Data- n = ' num2str(n)],'Shuffled'},'location','southeast')
title('Sing Adapt- Pref test')
subplot(3,2,5)
Rsq_temp = R_square_noadapt_all;
Rsq_temp(find(Rsq_temp<0)) = 0;
cdfplot(Rsq_temp(prefmask_ind_all,:))
n = sum(~isnan(Rsq_temp(prefmask_ind_all,:)));
hold on
Rsq_temp = R_square_shuf_noadapt_all;
Rsq_temp(find(Rsq_temp<0)) = 0;
cdfplot(Rsq_temp(prefmask_ind_all,:))
xlabel('Rsquared')
xlim([0 1])
legend({['Data- n = ' num2str(n)],'Shuffled'},'location','southeast')
title('No Adapt- Pref mask')
subplot(3,2,6)
Rsq_temp = R_square_singadapt_all;
Rsq_temp(find(Rsq_temp<0)) = 0;
cdfplot(Rsq_temp(prefmask_ind_all,:))
n = sum(~isnan(Rsq_temp(prefmask_ind_all,:)));
hold on
Rsq_temp = R_square_shuf_singadapt_all;
Rsq_temp(find(Rsq_temp<0)) = 0;
cdfplot(Rsq_temp(prefmask_ind_all,:))
xlabel('Rsquared')
xlim([0 1])
legend({['Data- n = ' num2str(n)],'Shuffled'},'location','southeast')
title('Sing Adapt- Pref mask')
print(fullfile(summaryDir, [svName '_Rsq_AllTrialFits.pdf']),'-dpdf','-fillpage')

%base
figure;
subplot(3,2,1)
b_temp = b_noadapt_all;
cdfplot(b_temp(resp_ind_all,:))
n = sum(~isnan(b_temp(resp_ind_all,:)));
xlabel('Baseline')
xlim([-1 1])
legend(['Data- n = ' num2str(n)],'location','southeast')
title('No Adapt')
subplot(3,2,2)
b_temp = b_singadapt_all;
cdfplot(b_temp(resp_ind_all,:))
n = sum(~isnan(b_temp(resp_ind_all,:)));
xlabel('Baseline')
xlim([-1 1])
legend(['Data- n = ' num2str(n)],'location','southeast')
title('Sing Adapt')
subplot(3,2,3)
b_temp = b_noadapt_all;
cdfplot(b_temp(preftest_ind_all,:))
n = sum(~isnan(b_temp(preftest_ind_all,:)));
xlabel('Baseline')
xlim([-1 1])
legend(['Data- n = ' num2str(n)],'location','southeast')
title('No Adapt- Pref test')
subplot(3,2,4)
b_temp = b_singadapt_all;
cdfplot(b_temp(preftest_ind_all,:))
n = sum(~isnan(b_temp(preftest_ind_all,:)));
xlabel('Baseline')
xlim([-1 1])
legend(['Data- n = ' num2str(n)],'location','southeast')
title('Sing Adapt- Pref test')
subplot(3,2,5)
b_temp = b_noadapt_all;
cdfplot(b_temp(prefmask_ind_all,:))
n = sum(~isnan(b_temp(prefmask_ind_all,:)));
xlabel('Baseline')
xlim([-1 1])
legend(['Data- n = ' num2str(n)],'location','southeast')
title('No Adapt- Pref mask')
subplot(3,2,6)
b_temp = b_singadapt_all;
cdfplot(b_temp(prefmask_ind_all,:))
n = sum(~isnan(b_temp(prefmask_ind_all,:)));
xlabel('Baseline')
xlim([-1 1])
legend(['Data- n = ' num2str(n)],'location','southeast')
title('Sing Adapt- Pref mask')
print(fullfile(summaryDir, [svName '_Baseline_AllTrialFits.pdf']),'-dpdf','-fillpage')

%amp
figure;
subplot(3,2,1)
amp_temp = amp_noadapt_all;
amp_temp(find(amp_temp<0)) = 0;
cdfplot(amp_temp(resp_ind_all,:))
n = sum(~isnan(amp_temp(resp_ind_all,:)));
hold on
amp_temp = amp_shuf_noadapt_all;
amp_temp(find(amp_temp<0)) = 0;
cdfplot(amp_temp(resp_ind_all,:))
xlabel('Amplitude')
xlim([0 1])
legend({['Data- n = ' num2str(n)],'Shuffled'},'location','southeast')
title('No Adapt')
subplot(3,2,2)
amp_temp = amp_singadapt_all;
amp_temp(find(amp_temp<0)) = 0;
cdfplot(amp_temp(resp_ind_all,:))
n = sum(~isnan(amp_temp(resp_ind_all,:)));
hold on
amp_temp = amp_shuf_singadapt_all;
amp_temp(find(amp_temp<0)) = 0;
cdfplot(amp_temp(resp_ind_all,:))
xlabel('Amplitude')
xlim([0 1])
legend({['Data- n = ' num2str(n)],'Shuffled'},'location','southeast')
title('Sing Adapt')
subplot(3,2,3)
amp_temp = amp_noadapt_all;
amp_temp(find(amp_temp<0)) = 0;
cdfplot(amp_temp(preftest_ind_all,:))
n = sum(~isnan(amp_temp(preftest_ind_all,:)));
hold on
amp_temp = amp_shuf_noadapt_all;
amp_temp(find(amp_temp<0)) = 0;
cdfplot(amp_temp(preftest_ind_all,:))
xlabel('Amplitude')
xlim([0 1])
legend({['Data- n = ' num2str(n)],'Shuffled'},'location','southeast')
title('No Adapt- Pref test')
subplot(3,2,4)
amp_temp = amp_singadapt_all;
amp_temp(find(amp_temp<0)) = 0;
cdfplot(amp_temp(preftest_ind_all,:))
n = sum(~isnan(amp_temp(preftest_ind_all,:)));
hold on
amp_temp = amp_shuf_singadapt_all;
amp_temp(find(amp_temp<0)) = 0;
cdfplot(amp_temp(preftest_ind_all,:))
xlabel('Amplitude')
xlim([0 1])
legend({['Data- n = ' num2str(n)],'Shuffled'},'location','southeast')
title('Sing Adapt- Pref test')
subplot(3,2,5)
amp_temp = amp_noadapt_all;
amp_temp(find(amp_temp<0)) = 0;
cdfplot(amp_temp(prefmask_ind_all,:))
n = sum(~isnan(amp_temp(prefmask_ind_all,:)));
hold on
amp_temp = amp_shuf_noadapt_all;
amp_temp(find(amp_temp<0)) = 0;
cdfplot(amp_temp(prefmask_ind_all,:))
xlabel('Amplitude')
xlim([0 1])
legend({['Data- n = ' num2str(n)],'Shuffled'},'location','southeast')
title('No Adapt- Pref mask')
subplot(3,2,6)
amp_temp = amp_singadapt_all;
amp_temp(find(amp_temp<0)) = 0;
cdfplot(amp_temp(prefmask_ind_all,:))
n = sum(~isnan(amp_temp(prefmask_ind_all,:)));
hold on
amp_temp = amp_shuf_singadapt_all;
amp_temp(find(amp_temp<0)) = 0;
cdfplot(amp_temp(prefmask_ind_all,:))
xlabel('Amplitude')
xlim([0 1])
legend({['Data- n = ' num2str(n)],'Shuffled'},'location','southeast')
title('Sing Adapt- Pref mask')
print(fullfile(summaryDir, [svName '_Amplitude_AllTrialFits.pdf']),'-dpdf','-fillpage')

figure;
subplot(3,2,1)
cdfplot(b_noadapt_all(preftest_ind_all,:))
hold on
cdfplot(b_singadapt_all(preftest_ind_all,:))
xlabel('Baseline')
xlim([-1 1])
legend({'Control','Adapt'},'location','southeast')
[h p] = ttest(b_noadapt_all(preftest_ind_all,:), b_singadapt_all(preftest_ind_all,:));
title(['Pref test- n = ' num2str(length(preftest_ind_all)) '; p = ' num2str(chop(p,2))])
subplot(3,2,2)
cdfplot(b_noadapt_all(prefmask_ind_all,:))
hold on
cdfplot(b_singadapt_all(prefmask_ind_all,:))
xlabel('Baseline')
xlim([-1 1])
legend({'Control','Adapt'},'location','southeast')
[h p] = ttest(b_noadapt_all(prefmask_ind_all,:), b_singadapt_all(prefmask_ind_all,:));
title(['Pref mask- n = ' num2str(length(prefmask_ind_all)) '; p = ' num2str(chop(p,2))])
subplot(3,2,3)
amp_temp = amp_noadapt_all-amp_shuf_noadapt_all;
cdfplot(amp_temp(preftest_ind_all,:))
n = sum(~isnan(amp_temp(preftest_ind_all,:)));
hold on
amp_temp_a = amp_singadapt_all-amp_shuf_singadapt_all;
cdfplot(amp_temp_a(preftest_ind_all,:))
xlabel('Amplitude-Shuf')
xlim([-1 1])
legend({'Control','Adapt'},'location','southeast')
[h p] = ttest(amp_temp(preftest_ind_all,:), amp_temp_a(preftest_ind_all,:));
title(['p = ' num2str(chop(p,2))])
subplot(3,2,4)
amp_temp = amp_noadapt_all-amp_shuf_noadapt_all;
cdfplot(amp_temp(prefmask_ind_all,:))
n = sum(~isnan(amp_temp(prefmask_ind_all,:)));
hold on
amp_temp_a = amp_singadapt_all-amp_shuf_singadapt_all;
cdfplot(amp_temp_a(prefmask_ind_all,:))
xlabel('Amplitude-Shuf')
xlim([-1 1])
legend({'Control','Adapt'},'location','southeast')
[h p] = ttest(amp_temp(prefmask_ind_all,:), amp_temp_a(prefmask_ind_all,:));
title(['p = ' num2str(chop(p,2))])
subplot(3,2,5)
amp_temp = amp_noadapt_all-amp_shuf_noadapt_all;
cdfplot(amp_temp(intersect(find(npos==0),preftest_ind_all),:))
n = sum(~isnan(amp_temp(intersect(find(npos==0),preftest_ind_all),:)));
hold on
amp_temp_a = amp_singadapt_all-amp_shuf_singadapt_all;
cdfplot(amp_temp_a(intersect(find(npos==0),preftest_ind_all),:))
xlabel('Amplitude-Shuf')
xlim([-1 1])
legend({'Control','Adapt'},'location','southeast')
[h p] = ttest(amp_temp(intersect(find(npos==0),preftest_ind_all),:), amp_temp_a(intersect(find(npos==0),preftest_ind_all),:));
title(['Supp- p = ' num2str(chop(p,2))])
subplot(3,2,6)
amp_temp = amp_noadapt_all-amp_shuf_noadapt_all;
cdfplot(amp_temp(intersect(find(npos>0),preftest_ind_all),:))
n = sum(~isnan(amp_temp(intersect(find(npos>0),preftest_ind_all),:)));
hold on
amp_temp_a = amp_singadapt_all-amp_shuf_singadapt_all;
cdfplot(amp_temp_a(intersect(find(npos>0),preftest_ind_all),:))
xlabel('Amplitude-Shuf')
xlim([-1 1])
legend({'Control','Adapt'},'location','southeast')
[h p] = ttest(amp_temp(intersect(find(npos>0),preftest_ind_all),:), amp_temp_a(intersect(find(npos>0),preftest_ind_all),:));
title(['Facil- p = ' num2str(chop(p,2))])

print(fullfile(summaryDir, [svName '_Summary_AllTrialFits.pdf']),'-dpdf','-fillpage')


[h p_no_pt] = ttest(amp_noadapt_all(preftest_ind_all),amp_shuf_noadapt_all(preftest_ind_all),'tail','right');
[h p_sing_pt] = ttest(amp_singadapt_all(preftest_ind_all),amp_shuf_singadapt_all(preftest_ind_all),'tail','right');
[h p_no_ptf] = ttest(amp_noadapt_all(intersect(find(npos>0),preftest_ind_all)),amp_shuf_noadapt_all(intersect(find(npos>0),preftest_ind_all)),'tail','right');
[h p_sing_ptf] = ttest(amp_singadapt_all(intersect(find(npos>0),preftest_ind_all)),amp_shuf_singadapt_all(intersect(find(npos>0),preftest_ind_all)),'tail','right');
[h p_no_pts] = ttest(amp_noadapt_all(intersect(find(npos==0),preftest_ind_all)),amp_shuf_noadapt_all(intersect(find(npos==0),preftest_ind_all)),'tail','right');
[h p_sing_pts] = ttest(amp_singadapt_all(intersect(find(npos==0),preftest_ind_all)),amp_shuf_singadapt_all(intersect(find(npos==0),preftest_ind_all)),'tail','right');

figure;
subplot(2,2,1)
scatter(plaidSI_noadapt_all(preftest_ind_all,1),plaidSI_singadapt_all(preftest_ind_all,1))
xlabel('Control MI')
ylabel('Adapt MI')
refline(1)
axis square
[h, p_SI] = ttest(plaidSI_noadapt_all(preftest_ind_all,1),plaidSI_singadapt_all(preftest_ind_all,1));
subplot(2,2,2)
scatter(amp_noadapt_all(preftest_ind_all),amp_singadapt_all(preftest_ind_all))
xlabel('Control Amp')
ylabel('Adapt Amp')
xlim([0 1])
ylim([0 1])
refline(1)
axis square
[h p_amp] = ttest(amp_noadapt_all(preftest_ind_all),amp_singadapt_all(preftest_ind_all));
subplot(2,2,3)
scatter(amp_noadapt_all(intersect(find(R_square_noadapt_all>0.1),preftest_ind_all)),amp_singadapt_all(intersect(find(R_square_noadapt_all>0.1),preftest_ind_all)))
xlabel('Control Amp')
ylabel('Adapt Amp')
xlim([0 1])
ylim([0 1])
refline(1)
axis square
subplot(2,2,4)
scatter(pha_noadapt_all(intersect(find(R_square_noadapt_all>0.1),preftest_ind_all)),pha_singadapt_all(intersect(find(R_square_noadapt_all>0.1),preftest_ind_all)))
xlabel('Control Phase')
ylabel('Adapt Phase')
xlim([0 2.*pi])
ylim([0 2.*pi])
refline(1)
axis square
movegui('center')


figure;
subplot(2,2,1)
cdfplot(b_noadapt_all(preftest_ind_all,:))
hold on
cdfplot(b_singadapt_all(preftest_ind_all,:))
xlabel('Sine Baseline')
xlim([-1 1])
legend({'Control','Adapt'},'location','southeast')
[h p] = ttest(b_noadapt_all(preftest_ind_all,:), b_singadapt_all(preftest_ind_all,:));
title(['Pref test- n = ' num2str(length(preftest_ind_all)) '; p = ' num2str(chop(p,2))])
subplot(2,2,2)
cdfplot(amp_noadapt_all(preftest_ind_all,:)-amp_shuf_noadapt_all(preftest_ind_all,:))
hold on
cdfplot(amp_singadapt_all(preftest_ind_all,:)-amp_shuf_singadapt_all(preftest_ind_all,:))
xlabel('Sine Amplitude (-shuf)')
xlim([-1 1])
[h p] = ttest(amp_noadapt_all(preftest_ind_all,:)-amp_shuf_noadapt_all(preftest_ind_all,:), amp_singadapt_all(preftest_ind_all,:)-amp_shuf_singadapt_all(preftest_ind_all,:));
title(['Pref test- n = ' num2str(length(preftest_ind_all)) '; p = ' num2str(chop(p,2))])
print(fullfile(summaryDir, [svName '_SineBaseAmpSummary.pdf']),'-dpdf','-fillpage')