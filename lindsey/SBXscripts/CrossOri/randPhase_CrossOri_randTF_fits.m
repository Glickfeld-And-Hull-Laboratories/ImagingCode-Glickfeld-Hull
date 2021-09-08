clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandPhase2TF_ExptList';
eval(ds)
rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);
seed = rng;
doTFTuning = 0;
for iexp = 1:nexp

    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc{1};
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);

    LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

    fprintf(['2P imaging sine fitting analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
    for irun=1:nrun
        fprintf([ImgFolder{irun} ' - ' time{irun} '\n'])
    end

    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']))
    
    nCells = size(resp_cell{end,end,end,end},1);
    nTrials = size(stimCon_all,2);
    trial_n = zeros(nMaskCon,nStimCon,nMaskPhas,nTF);
    trialInd = cell(nMaskCon,nStimCon,nMaskPhas,nTF);
    for itf = 1:nTF
        ind_tf = find(TF_all == TFs(itf));
        for im = 1:nMaskCon
            ind_mask = find(maskCon_all == maskCons(im));
            for it = 1:nStimCon
                ind_stim = find(stimCon_all == stimCons(it));
                ind_sm = intersect(ind_mask,ind_stim);
                if it>1 & im>1
                    for ip = 1:nMaskPhas
                        ind_phase = find(maskPhas_all == maskPhas(ip));
                        ind = intersect(ind_phase,intersect(ind_tf,ind_sm));
                        trialInd{im,it,ip,itf} = ind;
                        trial_n(im,it,ip,itf) = length(ind);
                    end
                else
                    trialInd{im,it,1,itf} = intersect(ind_tf,ind_sm);
                    trial_n(im,it,itf) = length(trialInd{im,it,1,itf});
                end
            end
        end
    end

    p_anova = nan(nCells,nTF);
    b_hat = nan(nCells,nTF); 
    amp_hat = nan(nCells,nTF); 
    per_hat = nan(nCells,nTF); 
    pha_hat = nan(nCells,nTF); 
    sse = nan(nCells,nTF); 
    R_square = nan(nCells,nTF);
    yfit = nan(nCells,length(0:1:359),nTF);
    p_anova_shuf = nan(nCells,nTF);
    b_hat_shuf = nan(nCells,nTF); 
    amp_hat_shuf = nan(nCells,nTF); 
    per_hat_shuf = nan(nCells,nTF); 
    pha_hat_shuf = nan(nCells,nTF); 
    sse_shuf = nan(nCells,nTF); 
    R_square_shuf = nan(nCells,nTF);
    yfit_shuf = nan(nCells,length(0:1:359),nTF);
    SI_avg = nan(nCells,nMaskPhas,nTF);
            
    eye_n = nan(nTF,nMaskPhas);
    phase_range = 0:1:359;
%     if ~exist('centroid_dist_tf')
%         centroid_dist_tf = cell(1,nTF);
%         centroid_med_tf = cell(1,nTF);
%         ind = find(~isnan(centroid_stim(1,:)));
%         centroid_med = findMaxNeighbors(centroid_stim(:,ind),2);
%         for itf = 1:nTF
%             ind_tf = intersect(find(TF_all == TFs(itf)),find(~isnan(centroid_stim(1,:))));
%             centroid_med_tf{itf} = findMaxNeighbors(centroid_stim(:,ind_tf),2);
%             centroid_dist_tf{itf} = sqrt((centroid_stim(1,:)-centroid_med(1)).^2 + (centroid_stim(2,:)-centroid_med(2)).^2);
%         end
%         save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']), 'rect', 'Area', 'Centroid', 'SNR', 'Val', 'frame_rate' , 'rad_mat_start','centroid_mat_start', 'cStimOn', 'rad_base','rad_stim','centroid_base', 'centroid_stim', 'centroid_dist','centroid_med', 'centroid_dist_tf', 'centroid_med_tf');
%     end
    
    %less than 2 deg
    for itf = 1:nTF
        im = 2;
        it = 2;
        [memb ind_test] = ismember(trialInd{1,it,1,itf},find(centroid_dist<2));
        [memb ind_mask] = ismember(trialInd{im,1,1,itf},find(centroid_dist<2));
        test_avg = mean(resp_cell{1,it,1,itf}(:,find(ind_test)),2);
        test_avg_rect = test_avg;
        test_avg_rect(find(test_avg<0)) = 0;
        mask_avg = mean(resp_cell{im,1,1,itf}(:,find(ind_mask)),2);
        mask_avg_rect = mask_avg;
        mask_avg_rect(find(mask_avg<0)) = 0;
        resp_avg = nan(nCells,nMaskPhas);
        resp_all = [];
        stim_all = [];
        for ip = 1:nMaskPhas
            [memb ind] = ismember(trialInd{im,it,ip,itf},find(centroid_dist<2));
            resp_all = [resp_all resp_cell{im,it,ip,itf}(:,find(ind))];
            stim_all = [stim_all ip.*ones(size(resp_cell{im,it,ip,itf}(1,find(ind))))];
            resp_avg(:,ip) = mean(resp_cell{im,it,ip,itf}(:,find(ind)),2);
        end

        resp_all_rect = resp_all;
        resp_all_rect(find(resp_all<0)) = 0;
        SI_all = (resp_all_rect-(test_avg_rect+mask_avg_rect))./(resp_all_rect+(test_avg_rect+mask_avg_rect));
        resp_avg_rect = resp_avg;
        resp_avg_rect(find(resp_avg<0)) = 0;
        SI_avg(:,:,itf) = (resp_avg_rect-(test_avg_rect+mask_avg_rect))./(resp_avg_rect+(test_avg_rect+mask_avg_rect));
        [eye_n(itf,:) edges bin] = histcounts(stim_all,[1:5]);
        fprintf([num2str(TFs(itf)) ' TF Eye-n: ' num2str(eye_n(itf,:)) '\n'])
        if sum(squeeze(eye_n(itf,:))<4)==0
            figure;
            start = 1;
            n = 1;
            for iCell = 1:nCells
                if start>25
                    suptitle([mouse ' ' date '- Mask ' num2str(im) ' Test ' num2str(it) ' TF ' num2str(TFs(itf))])
                    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_SI' num2str(n) '_M' num2str(im) 'T' num2str(it) 'TF' num2str(itf) '_LT2degEyeMvmt.pdf']), '-dpdf','-fillpage')
                    figure;
                    start = 1;
                    n = n+1;
                end
                p_anova(iCell,itf) = anova1(resp_all(iCell,:), stim_all,'off');
                if max(SI_avg(iCell,:,itf),[],2)>min(SI_avg(iCell,:,itf),[],2)
                    [b_hat(iCell,itf), amp_hat(iCell,itf), per_hat(iCell,itf),pha_hat(iCell,itf),sse(iCell,itf),R_square(iCell,itf)] = sinefit(deg2rad(maskPhas(stim_all)),SI_all(iCell,:));
                    subplot(5,5,start)
                    scatter(maskPhas(stim_all),SI_all(iCell,:));
                    hold on
                    scatter(maskPhas,SI_avg(iCell,:,itf))
                    yfit(iCell,:,itf) = b_hat(iCell,itf)+amp_hat(iCell,itf).*(sin(2*pi*deg2rad(phase_range)./per_hat(iCell,itf) + 2.*pi/pha_hat(iCell,itf)));
                    plot(phase_range, yfit(iCell,:,itf));
                    title(['Rsq = ' num2str(chop(R_square(iCell,itf),2)) '; p = ' num2str(chop(p_anova(iCell,itf),2))])
                else
                    b_hat(iCell,itf) = max(SI_avg(iCell,:,itf),[],2);
                    amp_hat(iCell,itf) = 0;
                    per_hat(iCell,itf) = NaN;
                    pha_hat(iCell,itf) = NaN;
                    sse(iCell,itf) = 0;
                    R_square(iCell,itf)= 0;
                    yfit(iCell,:,itf) = nan(length(phase_range),1);
                end
                start = start+1;
            end
            suptitle([mouse ' ' date '- Mask ' num2str(im) ' Test ' num2str(it) ' TF ' num2str(TFs(itf))])
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_SI' num2str(n) '_M' num2str(im) 'T' num2str(it) 'TF' num2str(itf) '_LT2degEyeMvmt.pdf']), '-dpdf','-fillpage')
            
            stim_all_shuf = stim_all(randperm(length(stim_all)));
            figure;
            start = 1;
            n = 1;
            for iCell = 1:nCells
                if start>25
                    suptitle([mouse ' ' date '- Shuffled- TF ' num2str(TFs(itf)) ])
                    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_SI_TF' num2str(itf) '_shuffled' num2str(n) '_LT2degEyeMvmt.pdf']), '-dpdf','-fillpage')
                    figure;
                    start = 1;
                    n = n+1;
                end
                p_anova_shuf(iCell,itf) = anova1(resp_all(iCell,:), stim_all_shuf,'off');
                if max(SI_avg(iCell,:,itf),[],2)>min(SI_avg(iCell,:,itf),[],2)
                    [b_hat_shuf(iCell,itf), amp_hat_shuf(iCell,itf), per_hat_shuf(iCell,itf),pha_hat_shuf(iCell,itf),sse_shuf(iCell,itf),R_square_shuf(iCell,itf)] = sinefit(deg2rad(maskPhas(stim_all_shuf)),SI_all(iCell,:));
                    subplot(5,5,start)
                    scatter(maskPhas(stim_all_shuf),SI_all(iCell,:));
                    hold on
                    scatter(maskPhas,SI_avg(iCell,:,itf))
                    yfit_shuf(iCell,:,itf) = b_hat_shuf(iCell,itf)+amp_hat_shuf(iCell,itf).*(sin(2*pi*deg2rad(phase_range)./per_hat_shuf(iCell,itf) + 2.*pi/pha_hat_shuf(iCell,itf)));
                    plot(phase_range, yfit_shuf(iCell,:,itf));
                    title(['Rsq = ' num2str(chop(R_square_shuf(iCell,itf),2)) '; p = ' num2str(chop(p_anova_shuf(iCell,itf),2))])
                else
                    b_hat_shuf(iCell,itf) = max(SI_avg(iCell,:,itf),[],2);
                    amp_hat_shuf(iCell,itf) = 0;
                    per_hat_shuf(iCell,itf) = NaN;
                    pha_hat_shuf(iCell,itf) = NaN;
                    sse_shuf(iCell,itf) = 0;
                    R_square_shuf(iCell,itf)= 0;
                    yfit_shuf(iCell,:,itf) = nan(length(phase_range),1);
                end
                start = start+1;
            end
            suptitle([mouse ' ' date '- Shuffled- TF ' num2str(TFs(itf)) ])
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_SI_TF' num2str(itf) '_shuffled' num2str(n) '_LT2degEyeMvmt.pdf']), '-dpdf','-fillpage')
        end
    end
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits.mat']),'SI_avg','yfit', 'b_hat', 'amp_hat', 'per_hat', 'pha_hat', 'sse', 'R_square',  'p_anova', 'yfit_shuf', 'b_hat_shuf', 'amp_hat_shuf', 'per_hat_shuf', 'pha_hat_shuf', 'sse_shuf', 'R_square_shuf', 'p_anova_shuf','trial_n', 'trialInd', 'eye_n')



%%

    if nTF >1 & doTFTuning

        sum(p_anova(iCell,:)<0.05,1)

        plaid_resp = zeros(nCells,nTF);
        mask_resp = zeros(nCells,nTF);
        test_resp = zeros(nCells,nTF);
        plaidSI = zeros(nCells,nTF);
        testPI = zeros(nCells,nTF);

        for itf = 1:nTF
            plaid_resp(:,itf) = mean(resp_cell{end,end,1,itf},2);
            mask_resp(:,itf) =  mean(resp_cell{end,1,1,itf},2);
            test_resp(:,itf) =  mean(resp_cell{1,end,1,itf},2);
            plaid_resp(find(plaid_resp(:,itf)<0),itf) = 0;
            mask_resp(find(mask_resp(:,itf)<0),itf) = 0;
            test_resp(find(test_resp(:,itf)<0),itf) = 0;
            plaidSI(:,itf) = (plaid_resp(:,itf)-(mask_resp(:,itf)+test_resp(:,itf))) ./ (plaid_resp(:,itf) + mask_resp(:,itf) + test_resp(:,itf));
            testPI(:,itf) = abs((test_resp(:,itf)-mask_resp(:,itf)) ./ (mask_resp(:,itf)+test_resp(:,itf)));
        end

        figure;
        for itf = 1:nTF
            subplot(1,2,1)
            cdfplot(plaidSI(resp_ind_tf{itf},itf))
            xlim([-1 1])
            hold on
            subplot(1,2,2)
            cdfplot(testPI(resp_ind_tf{itf},itf))
            xlim([0 1])
            hold on
        end
        subplot(1,2,1)
        xlabel('Suppression index')
        subplot(1,2,2)
        legend(num2str(TFs'),'location', 'northwest')
        xlabel('Preference index')

        data_dfof_stim = zeros(nCells,3,nTF);
        for itf = 1:nTF
            data_dfof_stim(:,1,itf) = mean(resp_cell{end,1,1,itf},2);
            data_dfof_stim(:,2,itf) = mean(resp_cell{1,end,1,itf},2);
            data_dfof_stim(:,3,itf) = mean(resp_cell{end,end,1,itf},2);
        end

        [max_dir_val max_dir_ind] = max(squeeze(mean(data_dfof_stim,3)),[],2);
        max_tf_val = zeros(1,nCells);
        max_tf_ind = zeros(1,nCells);
        fit_out = cell(1,nCells);
        g_fit = cell(1,nCells);
        prefTF = zeros(1,nCells);
        % figure; movegui('center')
        % start = 1;
        for iCell = 1:nCells
        %     if start >49
        %         figure;
        %         start = 1;
        %     end
            [max_tf_val(1,iCell) max_tf_ind(1,iCell)] = max(data_dfof_stim(iCell,max_dir_ind(iCell),:),[],3);
            [fit_out{iCell} g_fit{iCell}] =fit(log2(TFs)',squeeze(data_dfof_stim(iCell,max_dir_ind(iCell),:)),'gauss1','Lower',[0 log2(TFs(1)) 0],'Upper',[Inf log2(TFs(end)) Inf]);
        %     subplot(7,7,start)
        %     plot(log2(TFs),squeeze(data_dfof_stim(iCell,max_dir_ind(iCell),:)),'o')
        %     hold on
        %     plot(fit_out{iCell})
            prefTF(1,iCell)= fit_out{iCell}.b1;
            RsqTF(1,iCell)= g_fit{iCell}.rsquare;
        %     if find(h_all == iCell)
        %         title('Sig')
        %     end
        %     start = start+1;
        %     legend('off')
            if rem(iCell, 10) == 0
                fprintf([num2str(iCell) '\n'])
            end
        end
        save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TFfits.mat']),'fit_out','g_fit','prefTF','RsqTF','max_tf_ind','max_dir_ind','plaidSI','testPI')

        figure; movegui('center');
        RsqTF(find(RsqTF<0)) = 0;
        subplot(2,1,1); histogram(RsqTF(resp_ind)); vline(0.8)
        ind = intersect(find(RsqTF>0.8),resp_ind); xlabel('Rsq')
        subplot(2,1,2); histogram(prefTF(ind),nTF)
        set(gca, 'XTick', log2(TFs), 'XTickLabels', TFs)
        vline(mean(prefTF(ind),2))
        xlabel('TF (cps)')

        figure;
        amp_diff = amp_hat-amp_hat_shuf;
        legstr = [];
        for itf = 1:nTF
            if length(resp_ind_tf{itf})>1
                if sum(~isnan(amp_diff(resp_ind_tf{itf},itf)),1)>1
                cdfplot(amp_diff(resp_ind_tf{itf},itf))
                hold on
                legstr = [legstr; [num2str(TFs(itf)) '- n = ' num2str(length(resp_ind_tf{itf}))]];
                end
            end
        end
        xlabel('Sine Amp: Data-Shuf')
        legend(legstr,'location', 'southeast')
        title('')
        print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_diffAmpByStimTF_cdfs.pdf']),'-dpdf','-bestfit');


        figure;
        movegui('center')
        for itf = 1:nTF
            subplot(3,2,itf)
            for iCell = 1:nCells
                if find(resp_ind_tf{itf} == iCell)
                    if p_anova(iCell,itf)<0.05
                        scatter(2.^prefTF(iCell), amp_diff(iCell,itf),'or')
                    else
                        scatter(2.^prefTF(iCell), amp_diff(iCell,itf),'ok')
                    end
                    hold on
                end
            end
            ylim([-0.2 0.5])
            xlim([0.02 0.32])
            ylabel('Sine Amp: Data-Shuf')
            xlabel('Preferred TF')
            title(['Stim TF: ' num2str(TFs(itf))])
        end
        subplot(3,2,6)
        for iCell = 1:nCells
            if find(resp_ind == iCell)
                if p_anova(iCell,max_tf_ind(1,iCell))<0.05
                    scatter(2.^prefTF(iCell), amp_diff(iCell,max_tf_ind(1,iCell)),'or')
                else
                    scatter(2.^prefTF(iCell), amp_diff(iCell,max_tf_ind(1,iCell)),'ok')
                end
                hold on
            end
        end
        ylabel('Sine Amp: Data-Shuf')
        xlabel('Preferred TF')
        title('At Peak TF')
        ylim([-0.2 0.5])
        xlim([0.02 0.32])
        print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_diffAmpByPrefTF_scatter.pdf']),'-dpdf','-fillpage');
    end
    close all
end