close all; clear all; clc;
doRedChannel = 0;
sfstr = {'high','low'};
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
nSF = size(sfstr,2);
eye_n_all =cell(1,nSF);
sigCells = cell(1,nSF);
totCells = cell(1,nSF);
fractSigCells = cell(1,nSF);
respCells = cell(1,nSF);
resp_ind_all = cell(1,nSF);

p_anova_each_all = cell(1,nSF);
b_each_all = cell(1,nSF);
amp_each_all = cell(1,nSF);
per_each_all = cell(1,nSF);
pha_each_all = cell(1,nSF);
sse_each_all = cell(1,nSF);
R_square_each_all = cell(1,nSF);
yfit_each_all = cell(1,nSF);

p_anova_all_all = cell(1,nSF);
b_all_all = cell(1,nSF);
amp_all_all = cell(1,nSF);
pha_all_all = cell(1,nSF);
per_all_all = cell(1,nSF);
sse_all_all = cell(1,nSF);
R_square_all_all = cell(1,nSF);
yfit_all_all = cell(1,nSF);

p_anova_shuf_all = cell(1,nSF);
b_shuf_all = cell(1,nSF);
amp_shuf_all = cell(1,nSF);
pha_shuf_all = cell(1,nSF);
per_shuf_all = cell(1,nSF);
sse_shuf_all = cell(1,nSF);
R_square_shuf_all = cell(1,nSF);
yfit_shuf_all = cell(1,nSF);
% phaseMod_all = cell(1,nSF);
% phaseMod_avg_all = cell(1,nSF);
plaidSI_all = cell(1,nSF);
testPI_all = cell(1,nSF);
f1f2_all = cell(1,nSF);
for iSF = 1:nSF
    ds = ['CrossOriRandPhase_' sfstr{iSF} 'SF_ExptList'];
    eval(ds)
    rc = behavConstsAV;
    frame_rate = 30;
    nexp = size(expt,2);
    eye_n_all{iSF} = zeros(nexp,4);
    sigCells{iSF} = zeros(nexp,3);
    totCells{iSF} = zeros(nexp,1);
    fractSigCells{iSF} = zeros(nexp,3,2);
    respCells{iSF} = zeros(nexp,1);
    
    p_anova_each_all{iSF} = [];
    b_each_all{iSF} = [];
    amp_each_all{iSF} = [];
    per_each_all{iSF} = [];
    pha_each_all{iSF} = [];
    sse_each_all{iSF} = [];
    R_square_each_all{iSF} = [];
    yfit_each_all{iSF} = [];
    
    p_anova_all_all{iSF} = [];
    yfit_all_all{iSF} = [];
    pha_all_all{iSF} = [];
    b_all_all{iSF} = [];
    amp_all_all{iSF} = [];
    per_all_all{iSF} = [];
    R_square_all_all{iSF} = [];
    sse_all_all{iSF} = [];
    
    p_anova_shuf_all{iSF} = [];
    yfit_shuf_all{iSF} = [];
    pha_shuf_all{iSF} = [];
    b_shuf_all{iSF} = [];
    amp_shuf_all{iSF} = [];
    per_shuf_all{iSF} = [];
    R_square_shuf_all{iSF} = [];
    sse_shuf_all{iSF} = [];
    
    resp_ind_all{iSF} = [];
    
    plaidSI_all{iSF} = [];
    testPI_all{iSF} = [];
    f1f2_all{iSF} = [];
%     phaseMod_all{iSF} = [];
%     phaseMod_avg_all{iSF} = [];
    fprintf([sfstr{iSF} ' SF\n'])
    for iexp = 1:nexp
        mouse = expt(iexp).mouse;
        date = expt(iexp).date;
        area = expt(iexp).img_loc{1};
        ImgFolder = expt(iexp).coFolder;
        nrun = length(ImgFolder);
        run_str = catRunName(cell2mat(ImgFolder), nrun);
        

        fprintf([mouse ' ' date '\n'])


        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits.mat']))
        
        
        nCells = size(resp_cell{end,end,end},1);
        totCells{iSF}(iexp,:) = nCells;
        eye_n_all{iSF}(iexp,:) = eye_n;
        
        p_anova(find(p_anova==0)) = NaN;
        p_anova_each_all{iSF} = [p_anova_each_all{iSF};  p_anova];
        b_each_all{iSF} = [b_each_all{iSF}; b_hat];
        amp_each_all{iSF} = [amp_each_all{iSF}; amp_hat];
        per_each_all{iSF} = [per_each_all{iSF}; per_hat];
        pha_each_all{iSF} = [pha_each_all{iSF}; pha_hat];
        sse_each_all{iSF} = [sse_each_all{iSF}; sse];
        R_square_each_all{iSF} = [R_square_each_all{iSF}; R_square];
        
        p_anova_all(find(p_anova_all==0)) = NaN;
        p_anova_all_all{iSF} = [p_anova_all_all{iSF};  p_anova_all];
        b_all_all{iSF} = [b_all_all{iSF}; b_hat_all];
        amp_all_all{iSF} = [amp_all_all{iSF}; amp_hat_all];
        pha_all_all{iSF} = [pha_all_all{iSF}; pha_hat_all];
        per_all_all{iSF} = [per_all_all{iSF}; per_hat_all];
        R_square_all_all{iSF} = [R_square_all_all{iSF}; R_square_all];
        sse_all_all{iSF} = [sse_all_all{iSF}; sse_all];
        
        p_anova_shuf(find(p_anova_shuf==0)) = NaN;
        p_anova_shuf_all{iSF} = [p_anova_shuf_all{iSF};  p_anova_shuf];
        b_shuf_all{iSF} = [b_shuf_all{iSF}; b_hat_shuf];
        amp_shuf_all{iSF} = [amp_shuf_all{iSF}; amp_hat_shuf];
        pha_shuf_all{iSF} = [pha_shuf_all{iSF}; pha_hat_shuf];
        per_shuf_all{iSF} = [per_shuf_all{iSF}; per_hat_shuf];
        R_square_shuf_all{iSF} = [R_square_shuf_all{iSF}; R_square_shuf];
        sse_shuf_all{iSF} = [sse_shuf_all{iSF}; sse_shuf];
        
        yfit_each_all{iSF} = cat(1,yfit_each_all{iSF}, yfit);
        yfit_all_all{iSF} = cat(1,yfit_all_all{iSF}, yfit_all);
        yfit_shuf_all{iSF} = cat(1,yfit_shuf_all{iSF}, yfit_shuf);
        
        resp_ind_all{iSF} = [resp_ind_all{iSF}; resp_ind+sum(totCells{iSF}(1:iexp-1,:),1)];
        
        plaid_resp = mean(resp_cell{end,end,1},2);
        mask_resp = mean(resp_cell{end,1,1},2);
        test_resp = mean(resp_cell{1,end,1},2);
        plaid_resp(find(plaid_resp<0)) = 0;
        mask_resp(find(mask_resp<0)) = 0;
        test_resp(find(test_resp<0)) = 0;
        
        plaidSI_all{iSF} = [plaidSI_all{iSF}; (plaid_resp-(mask_resp+test_resp)) ./ (plaid_resp + mask_resp + test_resp)];
        testPI_all{iSF} = [testPI_all{iSF}; abs((test_resp-mask_resp) ./ (mask_resp+test_resp))];
        
        if isfield(expt,'prFolder') & ~isempty(expt(iexp).prFolder)
            ImgFolder = expt(iexp).prFolder;
            nrun = length(ImgFolder);
            pr_run_str = catRunName(cell2mat(ImgFolder), nrun);
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' pr_run_str], [date '_' mouse '_' pr_run_str '_f1f2.mat']))
            f1f2_all{iSF} = [f1f2_all{iSF} f2overf1];
        else
            f1f2_all{iSF} = [f1f2_all{iSF} nan(1,nCells)];
        end
%         max_val = max(yfit,[],1);
%         min_val = min(yfit,[],1);
%         max_val(find(max_val<0)) = 0;
%         min_val(find(min_val<0)) = 0;
%         phaseMod_all{iSF} = [phaseMod_all{iSF}; squeeze((max_val-min_val)./(max_val+min_val))];
%         max_val = max(yfit_avg,[],1);
%         min_val = min(yfit_avg,[],1);
%         max_val(find(max_val<0)) = 0;
%         min_val(find(min_val<0)) = 0;
%         phaseMod_avg_all{iSF} = [phaseMod_avg_all{iSF}; squeeze((max_val-min_val)./(max_val+min_val))];
    end
end


%%

for iSF = 2%1:nSF
    figure;
    for im = 2:3
        for it = 2:3
        subplot(2,2,1)
        scatter(p_anova_each_all{iSF}(resp_ind_all{iSF},im,it), amp_each_all{iSF}(resp_ind_all{iSF},im,it))
        hold on
        xlabel('anova p-value')
        ylabel('Sine amplitude')
        ylim([0 2])

        subplot(2,2,2)
        Rsq_temp = R_square_each_all{iSF};
        Rsq_temp(find(Rsq_temp<0)) = 0;
        scatter(p_anova_each_all{iSF}(resp_ind_all{iSF},im,it), Rsq_temp(resp_ind_all{iSF},im,it))
        hold on
        xlabel('anova p-value')
        ylabel('Rsquare')
        ylim([0 1])

        subplot(2,2,3)
        scatter(p_anova_each_all{iSF}(resp_ind_all{iSF},im,it), amp_each_all{iSF}(resp_ind_all{iSF},im,it))
        hold on
        xlabel('anova p-value')
        ylabel('Sine amplitude')
        xlim([0 0.05])
        ylim([0 2])

        subplot(2,2,4)
        Rsq_temp = R_square_each_all{iSF};
        Rsq_temp(find(Rsq_temp<0)) = 0;
        scatter(p_anova_each_all{iSF}(resp_ind_all{iSF},im,it), Rsq_temp(resp_ind_all{iSF},im,it))
        hold on
        xlabel('anova p-value')
        ylabel('Rsquare')
        xlim([0 0.05])
        ylim([0 1])      
        end
    end
    suptitle(sfstr{iSF})
    %print(fullfile(summaryDir, ['randPhase_sineVsP_' sfstr{iSF} 'SF.pdf']),'-dpdf','-fillpage')
end

n = zeros(3,2,nSF);
colmat = strvcat('k','b');
sineAmp_avg = zeros(3,2,nSF);
figure;
for iSF = 1:nSF
    Rsq_temp = R_square_each_all{iSF};
    Rsq_temp(find(Rsq_temp<0)) = 0;
    for i = 1:3
        subplot(2,2,iSF)
        scatter(amp_each_all{iSF}(resp_ind_all{iSF},i),Rsq_temp(resp_ind_all{iSF},i))
        hold on
        ylabel('Rsq')
        xlabel('Sine amp')
        ylim([0 1])
        xlim([0 1])
        title(sfstr{iSF})
        ind = find(Rsq_temp(resp_ind_all{iSF},i)>0.1);
        ind_all = find(~isnan(Rsq_temp(resp_ind_all{iSF},i)));
        n(i,1,iSF) = length(ind);
        n(i,2,iSF) = length(ind_all);
        sineAmp_avg(i,1,iSF) = nanmean(amp_each_all{iSF}(resp_ind_all{iSF}(ind),i),1);
        sineAmp_avg(i,2,iSF) = nanstd(amp_each_all{iSF}(resp_ind_all{iSF}(ind),i),1)./sqrt(length(ind));
        CO = get(gca,'ColorOrder');
        text(.2,.8-((i-1).*0.1),num2str(n(i,2,iSF)),'Color',CO(i,:))
    end
    subplot(2,2,1)
    legend({'All', '<5deg', '<2deg'})
    subplot(2,2,3)
    scatter(1:3, n(:,1,iSF)./n(:,2,iSF),['o' colmat(iSF)])
    set(gca,'XTick', 1:3, 'XTickLabel',{'All', '<5deg', '<2deg'})
    ylabel('Fraction Rsq>0.1')
    legend(sfstr,'location', 'northwest')
    hold on
    xlim([0 4])
    ylim([0 1])
    subplot(2,2,4)
    errorbar(1:3, sineAmp_avg(:,1,iSF), sineAmp_avg(:,2,iSF),['o' colmat(iSF)])
    hold on
    set(gca,'XTick', 1:3, 'XTickLabel',{'All', '<5deg', '<2 deg'})
    ylabel('Sine amp')
    xlim([0 4])
    ylim([0 1])
end
print(fullfile(summaryDir, 'randPhase_compSineAmp_byRsq_SI.pdf'),'-dpdf','-fillpage')

% 
% figure;
% subplot(2,2,1)
% iSF = 2;
% Rsq_temp = R_square_all{iSF};
% Rsq_temp(find(Rsq_temp<0)) = 0;
% cdfplot(plaidSI_all{iSF}(intersect(resp_ind_all{iSF}, find(Rsq_temp<0.1))))
% hold on
% cdfplot(plaidSI_all{iSF}(intersect(resp_ind_all{iSF}, find(Rsq_temp>0.1))))
% xlabel('Suppression index')
% title('')
% legend({'Rsq<0.1','Rsq>0.1'},'location','northwest')
% subplot(2,2,2)
% cdfplot(testPI_all{iSF}(intersect(resp_ind_all{iSF}, find(Rsq_temp<0.1))))
% hold on
% cdfplot(testPI_all{iSF}(intersect(resp_ind_all{iSF}, find(Rsq_temp>0.1))))
% hold on
% xlabel('Preference index')
% title('')
% legend({'Rsq<0.1','Rsq>0.1'},'location','northwest')
% print(fullfile(summaryDir, 'randPhase_compPhaseMod_byRsq.pdf'),'-dpdf','-fillpage')

figure;
leg_str = cell(1,nSF);
for iSF = 1:nSF
    subplot(2,2,1)
    cdfplot(plaidSI_all{iSF}(resp_ind_all{iSF}))
    hold on
    subplot(2,2,2)
    cdfplot(testPI_all{iSF}(resp_ind_all{iSF}))
    hold on
    subplot(2,2,3)
    cdfplot(amp_all_all{iSF}(resp_ind_all{iSF}))
    hold on
    subplot(2,2,4)
    cdfplot(b_all_all{iSF}(resp_ind_all{iSF}))
    hold on
    leg_str{iSF} = [sfstr{iSF} '- n = ' num2str(length(resp_ind_all{iSF}))];
end
subplot(2,2,1)
xlabel('Modulation index')
title('')
legend(leg_str,'location', 'southeast')
subplot(2,2,2)
xlabel('Stim preference index')
title('')
subplot(2,2,3)
xlabel('Sine fit amplitudue')
title('')
subplot(2,2,4)    
xlabel('Sine fit baseline')
title('')
print(fullfile(summaryDir, 'randPhase_highlowSF_comp.pdf'),'-dpdf','-fillpage')
%%
nCon = 4;
stim = [];
for it = 2:3
    for im = 2:3
        stim = [stim; maskCons(im) stimCons(it)];
    end
end
for iSF = 1:2
    figure;
    pha_temp = reshape(pha_each_all{iSF}(resp_ind_all{iSF},2:3,2:3), [length(resp_ind_all{iSF}) 4]);
    for ic = 1:nCon
        subplot(nCon,2,1+(ic-1).*2)
        rose(pha_temp(:,ic) - pha_all_all{iSF}(resp_ind_all{iSF},:))
        title([num2str(chop(stim(ic,1),2)) ' Mask; ' num2str(chop(stim(ic,2),2)) ' Test'])
        subplot(nCon,2,2+(ic-1).*2)
        rose(pha_temp(:,ic)-pha_shuf_all{iSF}(resp_ind_all{iSF},:))
        title(['Shuffled'])
    end
    suptitle([sfstr{iSF} ' SF'])
    print(fullfile(summaryDir, ['randPhase_phaseDiff_SingVsAllTrialFits' [sfstr{iSF} ' SF'] '.pdf']),'-dpdf','-fillpage')
end


figure;
for iSF = 1:nSF
subplot(2,2,1+(iSF-1).*2)
Rsq_temp = R_square_all_all{iSF};
Rsq_temp(find(Rsq_temp<0)) = 0;
cdfplot(Rsq_temp(resp_ind_all{iSF},:))
n = sum(~isnan(Rsq_temp(resp_ind_all{iSF},:)));
hold on
Rsq_temp = R_square_shuf_all{iSF};
Rsq_temp(find(Rsq_temp<0)) = 0;
cdfplot(Rsq_temp(resp_ind_all{iSF},:))
xlabel('Rsquared')
title([sfstr{iSF} ' SF'])
xlim([0 1])
legend({['Data- n = ' num2str(n)],'Shuffled'},'location','southeast')
subplot(2,2,2+(iSF-1).*2)
cdfplot(amp_all_all{iSF}(resp_ind_all{iSF},:))
hold on
cdfplot(amp_shuf_all{iSF}(resp_ind_all{iSF},:))
xlabel('Sine Amplitude')
title([sfstr{iSF} ' SF'])
xlim([0 1])
legend({['Data- n = ' num2str(n)],'Shuffled'},'location','southeast')
end
print(fullfile(summaryDir, 'randPhase_RsqAmp_AllTrialFits.pdf'),'-dpdf','-fillpage')

figure;
for iSF = 1:nSF
subplot(2,2,1+(iSF-1).*2)
pha_temp = wrapTo2Pi(pha_all_all{iSF}(resp_ind_all{iSF},:));
hist(rad2deg(pha_temp))
xlabel('Phase (deg)')
title([sfstr{iSF} ' SF'])
subplot(2,2,2+(iSF-1).*2)
pha_temp2 = wrapTo2Pi(pha_shuf_all{iSF}(resp_ind_all{iSF},:));
hist(rad2deg(pha_temp2))
xlabel('Phase (deg)')
title('Shuffle')
end
print(fullfile(summaryDir, 'randPhase_phaseDist_AllTrialFits.pdf'),'-dpdf','-fillpage')


figure;
iSF = 2;
subplot(2,2,1)
ind = resp_ind_all{iSF};
cdfplot(amp_all_all{iSF}(ind(find(testPI_all{iSF}(ind,:)<0.3),:),1))
hold on
cdfplot(amp_all_all{iSF}(ind(find(testPI_all{iSF}(ind,:)>0.9),:),1))
xlabel('Sine amplitude')
title('')
legend({['PI<0.3- n=' num2str(length(find(testPI_all{iSF}(ind,:)<0.3)))] ,['PI>0.9- n=' num2str(length(find(testPI_all{iSF}(ind,:)>0.9)))]},'location','southeast')
subplot(2,2,3)
cdfplot(plaidSI_all{iSF}(ind(find(testPI_all{iSF}(ind,:)<0.3),:),1))
hold on
cdfplot(plaidSI_all{iSF}(ind(find(testPI_all{iSF}(ind,:)>0.9),:),1))
xlabel('Suppression index')
title('')
legend({'PI<0.3','PI>0.9'},'location','northwest')

subplot(2,2,2)
ind = resp_ind_all{iSF};
cdfplot(R_square_all_all{iSF}(ind(find(testPI_all{iSF}(ind,:)<0.3),:),1))
hold on
cdfplot(R_square_all_all{iSF}(ind(find(testPI_all{iSF}(ind,:)>0.9),:),1))
xlabel('Rsq')
title('')
legend({'PI<0.3', 'PI>0.9'},'location','southeast')

subplot(2,2,4)
ind = resp_ind_all{iSF};
cdfplot(amp_all_all{iSF}(ind(find(plaidSI_all{iSF}(ind,:)<0),:),1))
hold on
cdfplot(amp_all_all{iSF}(ind(find(plaidSI_all{iSF}(ind,:)>0),:),1))
xlabel('Sine Amp')
title('')
legend({['SI<0- n=' num2str(length(find(plaidSI_all{iSF}(ind,:)<0)))] ,['SI>0- n=' num2str(length(find(plaidSI_all{iSF}(ind,:)>0)))]},'location','southeast')
print(fullfile(summaryDir, 'randPhase_sineAmpvsPrefIndex.pdf'),'-dpdf','-fillpage')



figure;
[n edges bin] = histcounts(testPI_all{iSF});
for i = 1:length(n)
    ind = intersect(find(R_square_all_all{iSF}>0.1),intersect(resp_ind_all{iSF}, find(bin == i)));
    subplot(2,1,1)
    errorbar(nanmean(testPI_all{iSF}(ind,:),1), nanmean(amp_all_all{iSF}(ind,:),1),nanstd(amp_all_all{iSF}(ind,:),[],1)./sqrt(length(ind)),nanstd(amp_all_all{iSF}(ind,:),[],1)./sqrt(length(ind)),nanstd(testPI_all{iSF}(ind,:),[],1)./sqrt(length(ind)),nanstd(testPI_all{iSF}(ind,:),[],1)./sqrt(length(ind)),'ok')
    hold on
    subplot(2,1,2)
    errorbar(nanmean(testPI_all{iSF}(ind,:),1), nanmean(plaidSI_all{iSF}(ind,:),1),nanstd(plaidSI_all{iSF}(ind,:),[],1)./sqrt(length(ind)),nanstd(plaidSI_all{iSF}(ind,:),[],1)./sqrt(length(ind)),nanstd(testPI_all{iSF}(ind,:),[],1)./sqrt(length(ind)),nanstd(testPI_all{iSF}(ind,:),[],1)./sqrt(length(ind)),'ok')
    hold on
end
subplot(2,1,1)
ylim([0 1])
ylabel('Sine Amplitude')
xlabel('Ori preference index')
subplot(2,1,2)
ylim([-0.5 .5])
ylabel('Suppression index')
xlabel('Ori preference index')

figure;
iSF = 2;
ind = resp_ind_all{iSF};
subplot(2,2,1)
scatter(testPI_all{iSF}(ind,:), amp_all_all{iSF}(ind,:),'ok')
ylim([0 1])
xlim([0 1])
ylabel('Sine Amplitude')
xlabel('Ori preference index')
subplot(2,2,2)
scatter(testPI_all{iSF}(ind,:), plaidSI_all{iSF}(ind,:),'ok')
ylim([-1 1])
xlim([0 1])
ylabel('Suppression index')
xlabel('Ori preference index')
subplot(2,2,3)
scatter(testPI_all{iSF}(ind,:), f1f2_all{iSF}(:,ind),'ok')
ylim([0 5])
xlim([0 1])
ylabel('F2/F1')
xlabel('Ori preference index')
subplot(2,2,4)
scatter(amp_all_all{iSF}(ind,:), f1f2_all{iSF}(:,ind),'ok')
ylim([0 5])
xlim([0 1])
ylabel('F2/F1')
xlabel('Sine amplitude')
print(fullfile(summaryDir, 'randPhase_PrefIndexScatters.pdf'),'-dpdf','-fillpage')

figure;
subplot(2,2,1)
ind = resp_ind_all{iSF};
cdfplot(amp_all_all{iSF}(ind(find(f1f2_all{iSF}(:,ind)<1),:),1))
hold on
cdfplot(amp_all_all{iSF}(ind(find(f1f2_all{iSF}(:,ind)>1),:),1))
xlabel('Sine Amp')
title('')
legend({['F2/F1<1- n=' num2str(length(find(f1f2_all{iSF}(:,ind)<1)))] ,['F2/F1>1- n=' num2str(length(find(f1f2_all{iSF}(:,ind)>1)))]},'location','southeast')
subplot(2,2,2)
ind = resp_ind_all{iSF};
cdfplot(testPI_all{iSF}(ind(find(f1f2_all{iSF}(:,ind)<1),:),1))
hold on
cdfplot(testPI_all{iSF}(ind(find(f1f2_all{iSF}(:,ind)>1),:),1))
xlabel('Ori Pref')
title('')
subplot(2,2,3)
ind = resp_ind_all{iSF};
cdfplot(plaidSI_all{iSF}(ind(find(f1f2_all{iSF}(:,ind)<1),:),1))
hold on
cdfplot(plaidSI_all{iSF}(ind(find(f1f2_all{iSF}(:,ind)>1),:),1))
xlabel('Suppression index')
title('')
print(fullfile(summaryDir, 'randPhase_F2F1.pdf'),'-dpdf','-fillpage')


%%
figure;
amp_avg = zeros(3,2,nSF);
n = zeros(3,2,nSF);
for iSF = 1:nSF
    for i = 1:3
        n(i,1,iSF) = sum(~isnan(amp_all_all{iSF}(:,i)));
        amp_avg(i,1,iSF) = nanmean(amp_all_all{iSF}(:,i),1);
        amp_avg(i,2,iSF) = nanstd(amp_all_all{iSF}(:,i),[],1)./sqrt(n(i,1,iSF));
        n(i,2,iSF) = sum(~isnan(R_square_each_all{iSF}(resp_ind_all{iSF},i)));
    end
    subplot(2,2,1)
    plot(1:3, n(:,1,iSF)./n(:,2,iSF),['o' colmat(iSF)])
    hold on
    text(iSF,.8, num2str(n(:,2,iSF)),'Color', colmat(iSF))
    subplot(2,2,2)
    errorbar(1:3, amp_avg(:,1,iSF), amp_avg(:,2,iSF),['o' colmat(iSF)])
    hold on
    text(iSF,.3, num2str(n(:,1,iSF)),'Color', colmat(iSF))
end
subplot(2,2,1)
set(gca,'XTick', 1:3, 'XTickLabel',{'All', '<5deg', '<2deg'})
ylabel('Fract anova < 0.05')
xlim([0 4])
ylim([0 1])
subplot(2,2,2)
set(gca,'XTick', 1:3, 'XTickLabel',{'All', '<5deg', '<2deg'})
ylabel('Phase Mod')
ylabel('Sine Amp')
xlim([0 4])
ylim([0 1])
legend(sfstr)

amp_avg = zeros(3,2,nSF);
n = zeros(3,2,nSF);
for iSF = 1:nSF
    Rsq_temp = R_square_each_all{iSF};
    Rsq_temp(find(Rsq_temp<0)) = 0;
    for i = 1:3
        ind = find(Rsq_temp(resp_ind_all{iSF},i)>0.1);
        ind_all = find(~isnan(Rsq_temp(resp_ind_all{iSF},i)));
        n(i,1,iSF) = length(ind);
        n(i,2,iSF) = length(ind_all);
        amp_avg(i,1,iSF) = nanmean(amp_each_all{iSF}(resp_ind_all{iSF}(ind),i),1);
        amp_avg(i,2,iSF) = nanstd(amp_each_all{iSF}(resp_ind_all{iSF}(ind),i),1)./sqrt(length(ind));
    end
    subplot(2,2,3)
    scatter(1:3, n(:,1,iSF)./n(:,2,iSF),['o' colmat(iSF)])
    hold on
    text(iSF,.8, num2str(n(:,2,iSF)),'Color', colmat(iSF))
    subplot(2,2,4)
    errorbar(1:3, amp_avg(:,1,iSF), amp_avg(:,2,iSF),['o' colmat(iSF)])
    hold on
    text(iSF,.3, num2str(n(:,1,iSF)),'Color', colmat(iSF))
end
subplot(2,2,3)
set(gca,'XTick', 1:3, 'XTickLabel',{'All', '<5deg', '<2deg'})
ylabel('Fraction Rsq>0.1')
xlim([0 4])
ylim([0 1])
subplot(2,2,4)
set(gca,'XTick', 1:3, 'XTickLabel',{'All', '<5deg', '<2deg'})
ylabel('Sine Amp')
xlim([0 4])
ylim([0 1])
print(fullfile(summaryDir, 'randPhase_phaseFitSummary_SI.pdf'),'-dpdf','-fillpage')

%%
figure;
iSF = 2;
for i = 1:3
    subplot(1,3,i)
%     ind = find(Rsq_temp(resp_ind_all{iSF},i)>0.1);
    amp = amp_all_all{iSF}(resp_ind_all{iSF},i);
    PI = testPI_all{iSF}(resp_ind_all{iSF},:);
    [n edges bin] = histcounts(amp);
    for ii = 1:length(n)
        ind2 = find(bin == ii);
        errorbar(nanmean(amp(ind2,1),1), nanmean(PI(ind2,1),1), nanstd(PI(ind2,1),[],1)./sqrt(length(ind2)), nanstd(PI(ind2,1),[],1)./sqrt(length(ind2)), nanstd(amp(ind2,1),[],1)./sqrt(length(ind2)),  nanstd(amp(ind2,1),[],1)./sqrt(length(ind2)), 'ok');
        hold on
    end
    ylim([0 1])
    xlim([0 1])
    xlabel('Sine amp')
    ylabel('Preference index')
end

figure;
iSF = 2;
SI = plaidSI_all{iSF}(resp_ind_all{iSF},:);
PI = testPI_all{iSF}(resp_ind_all{iSF},:);
[n edges bin] = histcounts(PI);
subplot(2,2,1)
for ii = 1:length(n)
    ind2 = find(bin == ii);
    errorbar(nanmean(PI(ind2,1),1), nanmean(SI(ind2,1),1), nanstd(SI(ind2,1),[],1)./sqrt(length(ind2)), nanstd(SI(ind2,1),[],1)./sqrt(length(ind2)), nanstd(PI(ind2,1),[],1)./sqrt(length(ind2)),  nanstd(PI(ind2,1),[],1)./sqrt(length(ind2)), 'ok');
    hold on
end
ylim([-.5 0.5])
xlim([0 1])
ylabel('SI')
xlabel('Preference index')

amp = amp_each_all{iSF}(resp_ind_all{iSF},3);
subplot(2,2,2)
for ii = 1:length(n)
    ind2 = find(bin == ii);
    errorbar(nanmean(PI(ind2,1),1), nanmean(amp(ind2,1),1), nanstd(amp(ind2,1),[],1)./sqrt(length(ind2)), nanstd(amp(ind2,1),[],1)./sqrt(length(ind2)), nanstd(PI(ind2,1),[],1)./sqrt(length(ind2)),  nanstd(PI(ind2,1),[],1)./sqrt(length(ind2)), 'ok');
    hold on
end
ylim([0 1])
xlim([0 1])
ylabel('Sine amplitude')
xlabel('Preference index')

for  i = 3
    amp = amp_all_all{iSF}(resp_ind_all{iSF},i);
    ind3 = find(isnan(amp)+isnan(PI) == 0);
    lm = fitlm(amp(ind3),PI(ind3));
end


