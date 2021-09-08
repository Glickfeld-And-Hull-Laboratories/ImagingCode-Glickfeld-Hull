close all; clear all; clc;
doRedChannel = 0;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');

ds = ['CrossOriRandPhase2TF_ExptList'];
svName = 'randPhase2TF';
eval(ds)
driver = 'SLC';
area = 'V1';
doPlot = 1;
nTF = 2;
rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);

sigCells = [];
fractSigCells = [];
respCells = [];
resp_ind_all = [];

p_anova_all_all = [];
b_all_all = [];
amp_all_all = [];
pha_all_all = [];
per_all_all = [];
sse_all_all = [];
R_square_all_all = [];
yfit_all_all = [];

p_anova_shuf_all = [];
b_shuf_all = [];
amp_shuf_all = [];
pha_shuf_all = [];
per_shuf_all = [];
sse_shuf_all = [];
R_square_shuf_all = [];
yfit_shuf_all = [];

plaidSI_all = [];
testPI_all = [];

mouse_list = [];

data_dfof_con_ph_tc_avg_all = [];
SI_avg_all = [];
testPI_sign_all = [];

totCells = zeros(nexp,1);
for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    mouse_list = strvcat(mouse_list, mouse);
    date = expt(iexp).date;
    area = expt(iexp).img_loc{1};
    ImgFolder = expt(iexp).coFolder;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);

    fprintf([mouse ' ' date '\n'])

    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))

    TFs = unique(celleqel2mat_padded(input.tMaskOneGratingTemporalFreqCPS));
    
    nCells = size(resp_cell{end,end,end},1);
    totCells(iexp,:) = nCells;

    p_anova(find(p_anova==0)) = NaN;
    p_anova_all_all = [p_anova_all_all;  p_anova];
    b_all_all = [b_all_all; b_hat];
    amp_all_all = [amp_all_all; amp_hat];
    pha_all_all = [pha_all_all; pha_hat];
    per_all_all = [per_all_all; per_hat];
    R_square_all_all = [R_square_all_all; R_square];
    sse_all_all = [sse_all_all; sse];

    p_anova_shuf(find(p_anova_shuf==0)) = NaN;
    p_anova_shuf_all = [p_anova_shuf_all;  p_anova_shuf];
    b_shuf_all = [b_shuf_all; b_hat_shuf];
    amp_shuf_all = [amp_shuf_all; amp_hat_shuf];
    pha_shuf_all = [pha_shuf_all; pha_hat_shuf];
    per_shuf_all = [per_shuf_all; per_hat_shuf];
    R_square_shuf_all = [R_square_shuf_all; R_square_shuf];
    sse_shuf_all = [sse_shuf_all; sse_shuf];

    yfit_all_all = cat(1,yfit_all_all, yfit);
    yfit_shuf_all = cat(1,yfit_shuf_all, yfit_shuf);

    resp_ind_all = [resp_ind_all; resp_ind+sum(totCells(1:iexp-1,:),1)];

    plaid_resp = zeros(nCells,nTF);
    mask_resp = zeros(nCells,nTF);
    test_resp = zeros(nCells,nTF);
    for iTF = 1:nTF
        plaid_resp(:,iTF) = mean(resp_cell{end,end,1,iTF},2);
        mask_resp(:,iTF) = mean(resp_cell{end,1,1,iTF},2);
    end
   test_resp = mean(resp_cell{1,end,1,1},2);
    plaid_resp(find(plaid_resp<0)) = 0;
    mask_resp(find(mask_resp<0)) = 0;
    test_resp(find(test_resp<0)) = 0;

    plaidSI_all = [plaidSI_all; (plaid_resp-(mask_resp+test_resp)) ./ (plaid_resp + mask_resp + test_resp)];
    testPI_all = [testPI_all; abs((test_resp-mask_resp) ./ (mask_resp+test_resp))];
    testPI_sign_all = [testPI_sign_all; ((test_resp-mask_resp) ./ (mask_resp+test_resp))];
    
    data_dfof_con_ph_tc_avg_all = cat(2,data_dfof_con_ph_tc_avg_all,data_dfof_con_ph_tc_avg);
    SI_avg_all = [SI_avg_all; SI_avg];
end   
save(fullfile(summaryDir,[svName '_Summary_' area '_' driver '.mat']),'p_anova_all_all','b_all_all','amp_all_all','R_square_shuf_all','p_anova_shuf_all','b_shuf_all','amp_shuf_all','R_square_shuf_all','yfit_all_all','yfit_shuf_all','plaidSI_all','testPI_all','resp_ind_all','data_dfof_con_ph_tc_avg_all','mouse_list')


%%
if doPlot
figure; 
movegui('center')
subplot(2,2,1)
cdfplot(testPI_all(resp_ind_all,1))
hold on
cdfplot(testPI_all(resp_ind_all,2))
xlabel('Stimulus selectivity')
title('')
legend({[num2str(TFs(1)) ' Hz']; [num2str(TFs(2)) ' Hz']},'location','northwest')
subplot(2,2,2)
cdfplot(plaidSI_all(resp_ind_all,1))
hold on
cdfplot(plaidSI_all(resp_ind_all,2))
xlabel('Modulation index')
title('')
subplot(2,2,3)
cdfplot(b_all_all(resp_ind_all,1))
hold on
cdfplot(b_all_all(resp_ind_all,2))
xlabel('Sine baseline')
title('')
subplot(2,2,4)
cdfplot(amp_all_all(resp_ind_all,1)-amp_shuf_all(resp_ind_all,1))
[h_tf1 p_tf1] = ttest(amp_all_all(resp_ind_all,1)-amp_shuf_all(resp_ind_all,1));
hold on
cdfplot(amp_all_all(resp_ind_all,2)-amp_shuf_all(resp_ind_all,2))
[h_tf2 p_tf2] = ttest(amp_all_all(resp_ind_all,2)-amp_shuf_all(resp_ind_all,2));
xlabel('Sine amplitude (-shuf)')
title('')
legend({['p = ' num2str(chop(p_tf1,2))]; ['p = ' num2str(chop(p_tf2,2))]},'location','southeast')
sgtitle(['n = ' num2str(size(mouse_list,1)) ' mice; n = ' num2str(length(resp_ind_all)) ' cells'])
print(fullfile(summaryDir,[svName '_Summary_' area '_' driver '.pdf']),'-dpdf')
savefig(fullfile(summaryDir,[svName '_Summary_' area '_' driver '.fig']))

figure;
subplot(1,2,1)
scatter(amp_all_all(resp_ind_all,1),amp_all_all(resp_ind_all,2))
hold on
scatter(amp_all_all(intersect(find(testPI_sign_all>0),resp_ind_all),1),amp_all_all(intersect(find(testPI_sign_all>0),resp_ind_all),2))
xlabel('Sine amplitude- TF = 1 Hz')
ylabel('Sine amplitude- TF = 5 Hz')
axis square
xlim([0 1])
ylim([0 1])
refline(1)
[h p] = ttest(amp_all_all(resp_ind_all,1),amp_shuf_all(resp_ind_all,2));
subplot(1,2,2)
cdfplot(amp_all_all(resp_ind_all,1))
hold on
cdfplot(amp_all_all(resp_ind_all,2))
legend({[num2str(TFs(1)) ' Hz']; [num2str(TFs(2)) ' Hz']},'location','southeast')
axis square
xlabel('Sine amplitude')
title('')
sgtitle(['n = ' num2str(size(mouse_list,1)) ' mice; n = ' num2str(length(resp_ind_all)) ' cells; p = ' num2str(chop(p,2))])
print(fullfile(summaryDir,[svName '_Summary_' area '_' driver '_TFcomp.pdf']),'-dpdf')
savefig(fullfile(summaryDir,[svName '_Summary_' area '_' driver '_TFcomp.fig']))

ind_lowmod = intersect(find(testPI_sign_all>0),intersect(resp_ind_all,intersect(find(amp_all_all(:,1)>0.45),find(amp_all_all(:,2)<0.2))));
for i = 1:length(ind_lowmod)
    figure;
    x = ind_lowmod(i);
    temp = data_dfof_con_ph_tc_avg_all(:,x,:,:,:,:);
    ymax = 1.25.*max(temp(:));
    subplot(6,2,1)
    plot(tt,data_dfof_con_ph_tc_avg_all(:,x,1,2,1,1))
    hold on
    plot(tt,data_dfof_con_ph_tc_avg_all(:,x,2,1,1,1))
    plot(tt,mean(data_dfof_con_ph_tc_avg_all(:,x,2,2,:,1),5))
    title('1 Hz')
    ylim([-0.1 ymax])
    xlim([-1 4])
    subplot(6,2,2)
    plot(tt,data_dfof_con_ph_tc_avg_all(:,x,1,2,1,1))
    hold on
    plot(tt,data_dfof_con_ph_tc_avg_all(:,x,2,1,1,2))
    plot(tt,mean(data_dfof_con_ph_tc_avg_all(:,x,2,2,:,2),5))
    title('5 Hz')
    ylim([-0.1 ymax])
    xlim([-1 4])
    for ip = 1:4
        t1 = [3:2:11];
        t2 = [4:2:12];
        subplot(6,2,t1(ip))
        plot(tt,data_dfof_con_ph_tc_avg_all(:,x,2,2,ip,1))
        ylim([-0.1 ymax])
        xlim([-1 4])
        subplot(6,2,t2(ip))
        plot(tt,data_dfof_con_ph_tc_avg_all(:,x,2,2,ip,2))
        ylim([-0.1 ymax])
        xlim([-1 4])
    end
    subplot(6,2,t1(ip+1))
    scatter(maskPhas, SI_avg_all(x,:,1))
    hold on
    plot(0:359, yfit_all_all(x,:,1))
    xlabel('Phase (deg)')
    ylim([-1 1])
    subplot(6,2,t2(ip+1))
    scatter(maskPhas, SI_avg_all(x,:,2))
    hold on
    plot(0:359, yfit_all_all(x,:,2))
    xlabel('Phase (deg)')
    ylim([-1 1])
    sgtitle(num2str(x))
    print(fullfile(summaryDir,[svName '_Summary_' area '_' driver '_exCell_' num2str(x) '.pdf']),'-dpdf')
    savefig(fullfile(summaryDir,[svName '_Summary_' area '_' driver '_exCell_' num2str(x) '.fig']))
end

ind_highmod = intersect(resp_ind_all,intersect(find(amp_all_all(:,1)<0.1),find(amp_all_all(:,2)>0.3)));
for i = 1:length(ind_highmod)
    figure;
    x = ind_highmod(i);
    temp = data_dfof_con_ph_tc_avg_all(:,x,:,:,:,:);
    ymax = 1.25.*max(temp(:));
    subplot(6,2,1)
    plot(tt,data_dfof_con_ph_tc_avg_all(:,x,1,2,1,1))
    hold on
    plot(tt,data_dfof_con_ph_tc_avg_all(:,x,2,1,1,1))
    plot(tt,mean(data_dfof_con_ph_tc_avg_all(:,x,2,2,:,1),5))
    title('1 Hz')
    ylim([-0.1 ymax])
    xlim([-1 4])
    subplot(6,2,2)
    plot(tt,data_dfof_con_ph_tc_avg_all(:,x,1,2,1,1))
    hold on
    plot(tt,data_dfof_con_ph_tc_avg_all(:,x,2,1,1,2))
    plot(tt,mean(data_dfof_con_ph_tc_avg_all(:,x,2,2,:,2),5))
    title('5 Hz')
    ylim([-0.1 ymax])
    xlim([-1 4])
    for ip = 1:4
        t1 = [3:2:11];
        t2 = [4:2:12];
        subplot(6,2,t1(ip))
        plot(tt,data_dfof_con_ph_tc_avg_all(:,x,2,2,ip,1))
        ylim([-0.1 ymax])
        xlim([-1 4])
        subplot(6,2,t2(ip))
        plot(tt,data_dfof_con_ph_tc_avg_all(:,x,2,2,ip,2))
        ylim([-0.1 ymax])
        xlim([-1 4])
    end
    subplot(6,2,t1(ip+1))
    scatter(maskPhas, SI_avg_all(x,:,1))
    hold on
    plot(0:359, yfit_all_all(x,:,1))
    xlabel('Phase (deg)')
    ylim([-1 1])
    subplot(6,2,t2(ip+1))
    scatter(maskPhas, SI_avg_all(x,:,2))
    hold on
    plot(0:359, yfit_all_all(x,:,2))
    xlabel('Phase (deg)')
    ylim([-1 1])
    sgtitle(num2str(x))
    print(fullfile(summaryDir,[svName '_Summary_' area '_' driver '_exCell_' num2str(x) '.pdf']),'-dpdf')
    savefig(fullfile(summaryDir,[svName '_Summary_' area '_' driver '_exCell_' num2str(x) '.fig']))
end
end
