close all; clear all; clc;
doRedChannel = 0;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');

doPlot = 1;
ds = ['CrossOriRandPhaseFF_ExptList'];
svName = 'randPhaseFF';
eval(ds)
driver = 'SCN';
area = 'V1';
doRedCells = 0;
SF = 0.05;
con = 0.5;
sz = 1000;
doSFSave = 1;
doConSave = 0;
doSzSave = 1;
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

p_anova_downsamp_all = [];
b_downsamp_all = [];
amp_downsamp_all = [];
pha_downsamp_all = [];
per_downsamp_all = [];
sse_downsamp_all = [];
R_square_downsamp_all = [];
yfit_downsamp_all = [];

p_anova_thresh_all = [];
b_thresh_all = [];
amp_thresh_all = [];
pha_thresh_all = [];
per_thresh_all = [];
R_square_thresh_all = [];
sse_thresh_all = [];

p_anova_shuf_thresh_all = [];
b_shuf_thresh_all = [];
amp_shuf_thresh_all = [];
pha_shuf_thresh_all = [];
per_shuf_thresh_all = [];
R_square_shuf_thresh_all = [];
sse_shuf_thresh_all = [];

test_thresh_resp_all = [];
mask_thresh_resp_all = [];

plaidSI_all = [];
testPI_all = [];
plaidSI_thresh_all = [];
testPI_thresh_all = [];
OSI_all = [];
max_dir_all = [];

red_cells_all = [];
z_all = [];

mouse_list = [];

totCells = zeros(nexp,1);
for iexp = 1:nexp
    if strcmp(expt(iexp).driver,driver) & (~isfield(expt,'con') || (isfield(expt,'con') & expt(iexp).con == con))...
            & (~isfield(expt,'SF') || (isfield(expt,'SF') & expt(iexp).SF == SF)) ...
            & (~isfield(expt,'size') || (isfield(expt,'size') & expt(iexp).size == sz))
        
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
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits.mat']))
        
        if doRedCells
            if isfield(expt,'coFolder')
                ImgFolder = expt(iexp).coFolder;
                nrun = length(ImgFolder);
                run_str = catRunName(cell2mat(ImgFolder), nrun);
            end
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']),'red_cells')
        else
            red_cells = [];
        end

        nCells = size(resp_cell{end,end,end},1);
        totCells(iexp,:) = nCells;

        p_anova_all(find(p_anova_all==0)) = NaN;
        p_anova_all_all = [p_anova_all_all;  p_anova_all];
        b_all_all = [b_all_all; b_hat_all];
        amp_all_all = [amp_all_all; amp_hat_all];
        pha_all_all = [pha_all_all; pha_hat_all];
        per_all_all = [per_all_all; per_hat_all];
        R_square_all_all = [R_square_all_all; R_square_all];
        sse_all_all = [sse_all_all; sse_all];

        p_anova_shuf(find(p_anova_shuf==0)) = NaN;
        p_anova_shuf_all = [p_anova_shuf_all;  p_anova_shuf];
        b_shuf_all = [b_shuf_all; b_hat_shuf];
        amp_shuf_all = [amp_shuf_all; amp_hat_shuf];
        pha_shuf_all = [pha_shuf_all; pha_hat_shuf];
        per_shuf_all = [per_shuf_all; per_hat_shuf];
        R_square_shuf_all = [R_square_shuf_all; R_square_shuf];
        sse_shuf_all = [sse_shuf_all; sse_shuf];
        
        if exist('p_anova_downsamp')
            p_anova_downsamp(find(p_anova_downsamp==0)) = NaN;
            p_anova_downsamp_all = [p_anova_downsamp_all;  p_anova_downsamp];
            b_downsamp_all = [b_downsamp_all; b_hat_downsamp];
            amp_downsamp_all = [amp_downsamp_all; amp_hat_downsamp];
            pha_downsamp_all = [pha_downsamp_all; pha_hat_downsamp];
            per_downsamp_all = [per_downsamp_all; per_hat_downsamp];
            R_square_downsamp_all = [R_square_downsamp_all; R_square_downsamp];
            sse_downsamp_all = [sse_downsamp_all; sse_downsamp];
            yfit_downsamp_all = cat(1,yfit_downsamp_all, yfit_downsamp);
        end
        
        yfit_all_all = cat(1,yfit_all_all, yfit_all);
        yfit_shuf_all = cat(1,yfit_shuf_all, yfit_shuf);

        resp_ind_all = [resp_ind_all; resp_ind+sum(totCells(1:iexp-1,:),1)];
        red_cells_all = [red_cells_all red_cells'+sum(totCells(1:iexp-1,:),1)];
        z_all = [z_all expt(iexp).z.*ones(1,sum(totCells(1:iexp-1,:),1))];

        plaid_resp = mean(resp_cell{end,end,1},2);
        mask_resp = mean(resp_cell{end,1,1},2);
        test_resp = mean(resp_cell{1,end,1},2);
        plaid_resp(find(plaid_resp<0)) = 0;
        mask_resp(find(mask_resp<0)) = 0;
        test_resp(find(test_resp<0)) = 0;

        plaidSI_all = [plaidSI_all; (plaid_resp-(mask_resp+test_resp)) ./ (plaid_resp + mask_resp + test_resp)];
        testPI_all = [testPI_all; abs((test_resp-mask_resp) ./ (mask_resp+test_resp))];
        
        plaid_thresh_resp = mean(resp_cell{end,end,1},2).*0.75;
        mask_thresh_resp = mean(resp_cell{end,1,1},2).*0.75;
        test_thresh_resp = mean(resp_cell{1,end,1},2).*0.75;
        plaid_thresh_resp(find(plaid_thresh_resp<0.04)) = 0;
        mask_thresh_resp(find(mask_thresh_resp<0.04)) = 0;
        test_thresh_resp(find(test_thresh_resp<0.04)) = 0;
        
        test_thresh_resp_all = [test_thresh_resp_all; test_thresh_resp];
        mask_thresh_resp_all = [mask_thresh_resp_all; mask_thresh_resp];
        
        plaidSI_thresh_all = [plaidSI_thresh_all; (plaid_thresh_resp-(mask_thresh_resp+test_thresh_resp)) ./ (plaid_thresh_resp + mask_thresh_resp + test_thresh_resp)];
        testPI_thresh_all = [testPI_thresh_all; abs((test_thresh_resp-mask_thresh_resp) ./ (mask_thresh_resp+test_thresh_resp))];
        
        if isfield(expt,'dirFolder')
            if ~isempty(expt(iexp).dirFolder)
                ImgFolder = expt(iexp).dirFolder;
                nrun = length(ImgFolder);
                run_str = catRunName(cell2mat(ImgFolder), nrun);
                load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']))
                OSI_all = [OSI_all; OSI_resp];
                max_dir_all = [max_dir_all; max_ind_dir]; 
            else
                OSI_all = [OSI_all; nan(size(b_all_all))];
                max_dir_all = [max_dir_all; nan(size(b_all_all))]; 
            end
        elseif ~isempty(expt(iexp).coFolder)
            ImgFolder = expt(iexp).coFolder;
            nrun = length(ImgFolder);
            run_str = catRunName(cell2mat(ImgFolder), nrun);
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
            OSI_all = [OSI_all; nan(size(b_all_all))];
            [max_val_dir, max_ind_dir] = max(avg_resp_dir(:,:,1,1),[],2);
            max_dir_all = [max_dir_all; max_ind_dir]; 
        else
            OSI_all = [OSI_all; nan(size(b_all_all))];
            max_dir_all = [max_dir_all; nan(size(b_all_all))]; 
        end
        if exist(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_thresh.mat']))
        	load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_thresh.mat']))
            p_anova_all(find(p_anova_all==0)) = NaN;
            p_anova_thresh_all = [p_anova_thresh_all;  p_anova_all];
            b_thresh_all = [b_thresh_all; b_hat_all];
            amp_thresh_all = [amp_thresh_all; amp_hat_all];
            pha_thresh_all = [pha_thresh_all; pha_hat_all];
            per_thresh_all = [per_thresh_all; per_hat_all];
            R_square_thresh_all = [R_square_thresh_all; R_square_all];
            sse_thresh_all = [sse_thresh_all; sse_all];

            p_anova_shuf(find(p_anova_shuf==0)) = NaN;
            p_anova_shuf_thresh_all = [p_anova_shuf_thresh_all;  p_anova_shuf];
            b_shuf_thresh_all = [b_shuf_thresh_all; b_hat_shuf];
            amp_shuf_thresh_all = [amp_shuf_thresh_all; amp_hat_shuf];
            pha_shuf_thresh_all = [pha_shuf_thresh_all; pha_hat_shuf];
            per_shuf_thresh_all = [per_shuf_thresh_all; per_hat_shuf];
            R_square_shuf_thresh_all = [R_square_shuf_thresh_all; R_square_shuf];
            sse_shuf_thresh_all = [sse_shuf_thresh_all; sse_shuf];
        end

        
    end
end
if doRedCells
    driver_str= 'noDriver';
else
    driver_str = driver;
end
sf_str = num2str(SF);
sf_str = ['pt' sf_str(3:end)];
con_str = num2str(con);
con_str = ['pt' con_str(3:end)];
sz_str = num2str(sz);

if doSFSave & doSzSave
    save(fullfile(summaryDir,[svName '_Summary_' area '_' driver_str '_SF' sf_str '_Sz' sz_str '.mat']),'z_all','red_cells_all','p_anova_all_all','b_all_all','amp_all_all','R_square_shuf_all','p_anova_shuf_all','b_shuf_all','amp_shuf_all','R_square_shuf_all','yfit_all_all','yfit_shuf_all','plaidSI_all','testPI_all','resp_ind_all','OSI_all','max_dir_all','mouse_list')
elseif doSzSave
    save(fullfile(summaryDir,[svName '_Summary_' area '_' driver_str '_Sz' sz_str '.mat']),'z_all','red_cells_all','p_anova_all_all','b_all_all','amp_all_all','R_square_shuf_all','p_anova_shuf_all','b_shuf_all','amp_shuf_all','R_square_shuf_all','yfit_all_all','yfit_shuf_all','plaidSI_all','testPI_all','resp_ind_all','OSI_all','max_dir_all','mouse_list')    
elseif doConSave
    save(fullfile(summaryDir,[svName '_Summary_' area '_' driver_str '_Con' con_str '.mat']),'z_all','red_cells_all','p_anova_all_all','b_all_all','amp_all_all','R_square_shuf_all','p_anova_shuf_all','b_shuf_all','amp_shuf_all','R_square_shuf_all','yfit_all_all','yfit_shuf_all','plaidSI_all','testPI_all','resp_ind_all','OSI_all','max_dir_all','mouse_list')
else    
    save(fullfile(summaryDir,[svName '_Summary_' area '_' driver_str '.mat']),'z_all','red_cells_all','p_anova_all_all','b_all_all','amp_all_all','R_square_shuf_all','p_anova_shuf_all','b_shuf_all','amp_shuf_all','R_square_shuf_all','yfit_all_all','yfit_shuf_all','plaidSI_all','testPI_all','resp_ind_all','OSI_all','max_dir_all','mouse_list')
end

%%
if doPlot
    resp_ind_all = setdiff(resp_ind_all, red_cells_all);
    figure; 
    subplot(2,2,1)
    cdfplot(testPI_all(resp_ind_all))
    xlabel('Stimulus selectivity')
    title('')
    subplot(2,2,2)
    cdfplot(plaidSI_all(resp_ind_all))
    xlabel('Modulation index')
    title('')
    subplot(2,2,3)
    cdfplot(b_all_all(resp_ind_all))
    xlabel('Sine baseline')
    xlim([-1 1])
    title('')
    subplot(2,2,4)
    cdfplot(amp_all_all(resp_ind_all)-amp_shuf_all(resp_ind_all))
    xlabel('Sine amplitude')
    xlim([-0.2 1])
    title('')
    suptitle({[area '- n = ' num2str(size(mouse_list,1)) ' expts; ' num2str(size(unique(mouse_list,'rows'),1)) ' mice'], ['All responsive cells- n = ' num2str(length(resp_ind_all))]})
    print(fullfile(summaryDir,[svName '_Summary_' area '_' driver_str '_SF' sf_str '_Sz' sz_str '_Con' con_str '.pdf']),'-dpdf')
    
    if length(red_cells_all)
        figure;
        subplot(2,2,1)
        cdfplot(testPI_all(red_cells_all))
        xlabel('Stimulus selectivity')
        title('')
        subplot(2,2,2)
        cdfplot(plaidSI_all(red_cells_all))
        xlabel('Modulation index')
        title('')
        subplot(2,2,3)
        cdfplot(b_all_all(red_cells_all))
        xlabel('Sine baseline')
        xlim([-1 1])
        title('')
        subplot(2,2,4)
        cdfplot(amp_all_all(red_cells_all)-amp_shuf_all(red_cells_all))
        xlabel('Sine amplitude')
        xlim([-0.2 1])
        title('')
        suptitle({[area '- n = ' num2str(size(mouse_list,1)) ' expts; ' num2str(size(unique(mouse_list,'rows'),1)) ' mice'], ['All ' driver ' cells- n = ' num2str(length(red_cells_all))]})
        print(fullfile(summaryDir,[svName '_Summary_' area '_' driver '_SF' sf_str '_Sz' sz_str '_Con' con_str '.pdf']),'-dpdf')
    end

% figure;
% subplot(2,2,1)
% Rsq_temp = R_square_all_all;
% Rsq_temp(find(Rsq_temp<0)) = 0;
% cdfplot(Rsq_temp(resp_ind_all,:))
% n = sum(~isnan(Rsq_temp(resp_ind_all,:)));
% hold on
% Rsq_temp = R_square_shuf_all;
% Rsq_temp(find(Rsq_temp<0)) = 0;
% cdfplot(Rsq_temp(resp_ind_all,:))
% if ~isempty(R_square_downsamp_all)
% Rsq_temp = R_square_downsamp_all;
% Rsq_temp(find(Rsq_temp<0)) = 0;
% cdfplot(Rsq_temp(resp_ind_all,:))
% end
% xlabel('Rsquared')
% xlim([0 1])
% legend({['Data- n = ' num2str(n)],'Shuffled','Downsampled'},'location','southeast')
% subplot(2,2,2)
% cdfplot(amp_all_all(resp_ind_all,:))
% hold on
% cdfplot(amp_shuf_all(resp_ind_all,:))
% if ~isempty(amp_downsamp_all)
% cdfplot(amp_downsamp_all(resp_ind_all,:))
% end
% xlabel('Sine Amplitude')
% xlim([0 1])
% legend({['Data- n = ' num2str(n)],'Shuffled','Downsampled'},'location','southeast')
% print(fullfile(summaryDir, [svName '_RsqAmp_AllTrialFits_8phase.pdf']),'-dpdf','-fillpage')
% 
% 
% figure;
% subplot(2,2,1)
% ind = resp_ind_all;
% cdfplot(amp_all_all(ind(find(testPI_all(ind,:)<0.3),:),1))
% hold on
% cdfplot(amp_all_all(ind(find(testPI_all(ind,:)>0.9),:),1))
% xlabel('Sine amplitude')
% title('')
% legend({['PI<0.3- n=' num2str(length(find(testPI_all(ind,:)<0.3)))] ,['PI>0.9- n=' num2str(length(find(testPI_all(ind,:)>0.9)))]},'location','southeast')
% subplot(2,2,3)
% cdfplot(plaidSI_all(ind(find(testPI_all(ind,:)<0.3),:),1))
% hold on
% cdfplot(plaidSI_all(ind(find(testPI_all(ind,:)>0.9),:),1))
% xlabel('Suppression index')
% title('')
% legend({'PI<0.3','PI>0.9'},'location','northwest')
% 
% subplot(2,2,2)
% ind = resp_ind_all;
% cdfplot(R_square_all_all(ind(find(testPI_all(ind,:)<0.3),:),1))
% hold on
% cdfplot(R_square_all_all(ind(find(testPI_all(ind,:)>0.9),:),1))
% xlabel('Rsq')
% title('')
% legend({'PI<0.3', 'PI>0.9'},'location','southeast')
% 
% subplot(2,2,4)
% ind = resp_ind_all;
% cdfplot(amp_all_all(ind(find(plaidSI_all(ind,:)<0),:),1))
% hold on
% cdfplot(amp_all_all(ind(find(plaidSI_all(ind,:)>0),:),1))
% xlabel('Sine Amp')
% title('')
% legend({['SI<0- n=' num2str(length(find(plaidSI_all(ind,:)<0)))] ,['SI>0- n=' num2str(length(find(plaidSI_all(ind,:)>0)))]},'location','southeast')
% print(fullfile(summaryDir, [svName '_sineAmpvsPrefIndex_8phase.pdf']),'-dpdf','-fillpage')
% 
% if ~isnan(OSI_all)
% figure;
% subplot(2,2,1)
% ind = resp_ind_all;
% cdfplot(amp_all_all(ind(find(OSI_all(ind,:)<0.5),:),1))
% hold on
% cdfplot(amp_all_all(ind(find(OSI_all(ind,:)>0.5),:),1))
% xlabel('Sine amplitude')
% title('')
% legend({['OSI<0.5- n=' num2str(length(find(OSI_all(ind,:)<0.5)))] ,['OSI>0.5- n=' num2str(length(find(OSI_all(ind,:)>0.5)))]},'location','southeast')
% subplot(2,2,3)
% cdfplot(plaidSI_all(ind(find(OSI_all(ind,:)<0.5),:),1))
% hold on
% cdfplot(plaidSI_all(ind(find(OSI_all(ind,:)>0.5),:),1))
% xlabel('Suppression index')
% title('')
% legend({'OSI<0.5','OSI>0.5'},'location','northwest')
% 
% subplot(2,2,2)
% ind = resp_ind_all;
% cdfplot(R_square_all_all(ind(find(OSI_all(ind,:)<0.5),:),1))
% hold on
% cdfplot(R_square_all_all(ind(find(OSI_all(ind,:)>0.5),:),1))
% xlabel('Rsq')
% title('')
% legend({'OSI<0.5', 'OSI>0.5'},'location','southeast')
% print(fullfile(summaryDir, [svName '_sineAmpvsOSI_8phase.pdf']),'-dpdf','-fillpage')
% 
% figure;
% [n edges bin] = histcounts(testPI_all);
% for i = 1:length(n)
%     ind = intersect(find(R_square_all_all>0.1),intersect(resp_ind_all, find(bin == i)));
%     subplot(2,2,1)
%     errorbar(nanmean(testPI_all(ind,:),1), nanmean(amp_all_all(ind,:),1),nanstd(amp_all_all(ind,:),[],1)./sqrt(length(ind)),nanstd(amp_all_all(ind,:),[],1)./sqrt(length(ind)),nanstd(testPI_all(ind,:),[],1)./sqrt(length(ind)),nanstd(testPI_all(ind,:),[],1)./sqrt(length(ind)),'ok')
%     hold on
%     subplot(2,2,2)
%     errorbar(nanmean(testPI_all(ind,:),1), nanmean(plaidSI_all(ind,:),1),nanstd(plaidSI_all(ind,:),[],1)./sqrt(length(ind)),nanstd(plaidSI_all(ind,:),[],1)./sqrt(length(ind)),nanstd(testPI_all(ind,:),[],1)./sqrt(length(ind)),nanstd(testPI_all(ind,:),[],1)./sqrt(length(ind)),'ok')
%     hold on
% end
% subplot(2,2,1)
% ylim([0 1])
% ylabel('Sine Amplitude')
% xlabel('Ori preference index')
% subplot(2,2,2)
% ylim([-0.5 .5])
% ylabel('Suppression index')
% xlabel('Ori preference index')
% [n edges bin] = histcounts(OSI_all,5);
% for i = 1:length(n)
%     ind = intersect(find(R_square_all_all>0.1),intersect(resp_ind_all, find(bin == i)));
%     subplot(2,2,3)
%     errorbar(nanmean(OSI_all(ind,:),1), nanmean(amp_all_all(ind,:),1),nanstd(amp_all_all(ind,:),[],1)./sqrt(length(ind)),nanstd(amp_all_all(ind,:),[],1)./sqrt(length(ind)),nanstd(OSI_all(ind,:),[],1)./sqrt(length(ind)),nanstd(OSI_all(ind,:),[],1)./sqrt(length(ind)),'ok')
%     hold on
%     subplot(2,2,4)
%     errorbar(nanmean(OSI_all(ind,:),1), nanmean(plaidSI_all(ind,:),1),nanstd(plaidSI_all(ind,:),[],1)./sqrt(length(ind)),nanstd(plaidSI_all(ind,:),[],1)./sqrt(length(ind)),nanstd(OSI_all(ind,:),[],1)./sqrt(length(ind)),nanstd(OSI_all(ind,:),[],1)./sqrt(length(ind)),'ok')
%     hold on
% end
% subplot(2,2,3)
% ylim([0 1])
% ylabel('Sine Amplitude')
% xlabel('OSI')
% subplot(2,2,4)
% ylim([-0.5 .5])
% ylabel('Suppression index')
% xlabel('OSI')
% print(fullfile(summaryDir, [svName '_OSIorPrefIndexBins_8phase.pdf']),'-dpdf','-fillpage')
% 
% 
% 
% figure;
% ind = resp_ind_all;
% max_dir_ind_reset = max_dir_all;
% max_dir_ind_reset(find(max_dir_all>180)) = max_dir_ind_reset(find(max_dir_all>180))-360;
% pref_dir_ind = intersect(ind,intersect(find(max_dir_ind_reset>-22.5),find(max_dir_ind_reset<22.5)));
% subplot(2,2,1)
% scatter(testPI_all(ind,:), amp_all_all(ind,:),'ok')
% hold on
% scatter(testPI_all(pref_dir_ind,:), amp_all_all(pref_dir_ind,:),'or')
% ylim([0 1])
% xlim([0 1])
% ylabel('Sine Amplitude')
% xlabel('Stim preference index')
% % subplot(2,3,2)
% % scatter(testPI_all(ind,:), amp_all_all(ind,:)-amp_shuf_all(ind,:),'ok')
% % hold on
% % scatter(testPI_all(pref_dir_ind,:), amp_all_all(pref_dir_ind,:)-amp_shuf_all(pref_dir_ind,:),'or')
% % ylim([0 1])
% % xlim([0 1])
% % ylabel('Sine Amplitude (Data-Shuf)')
% % xlabel('Stim preference index')
% subplot(2,2,2)
% scatter(testPI_all(ind,:), plaidSI_all(ind,:),'ok')
% hold on
% scatter(testPI_all(pref_dir_ind,:), plaidSI_all(pref_dir_ind,:),'or')
% ylim([-1 1])
% xlim([0 1])
% ylabel('Suppression index')
% xlabel('Stim preference index')
% subplot(2,2,3)
% scatter(OSI_all(ind,:), amp_all_all(ind,:),'ok')
% hold on
% scatter(OSI_all(pref_dir_ind,:), amp_all_all(pref_dir_ind,:),'or')
% ylim([0 1])
% xlim([0 1])
% ylabel('Sine Amplitude')
% xlabel('OSI')
% % subplot(2,3,5)
% % scatter(OSI_all(ind,:), amp_all_all(ind,:)-amp_shuf_all(ind,:),'ok')
% % hold on
% % scatter(OSI_all(pref_dir_ind,:), amp_all_all(pref_dir_ind,:)-amp_shuf_all(pref_dir_ind,:),'or')
% % ylim([0 1])
% % xlim([0 1])
% % ylabel('Sine Amplitude')
% % xlabel('OSI')
% subplot(2,2,4)
% scatter(OSI_all(ind,:), plaidSI_all(ind,:),'ok')
% hold on
% scatter(OSI_all(pref_dir_ind,:), plaidSI_all(pref_dir_ind,:),'or')
% ylim([-1 1])
% xlim([0 1])
% ylabel('Suppression index')
% xlabel('OSI')
% print(fullfile(summaryDir, [svName '_OSIorPrefIndexScatters_8phase.pdf']),'-dpdf','-fillpage')
end