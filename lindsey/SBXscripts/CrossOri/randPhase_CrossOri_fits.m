clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandDirRandPhase_ExptList';
eval(ds)
rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);
nanframes = zeros(1,nexp);
max_dist = 2;

for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc{1};
    ImgFolder = expt(iexp).copFolder;
    time = expt(iexp).copTime;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);
    

    LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

    fprintf(['2P imaging sine fitting analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
    for irun=1:nrun
        fprintf([ImgFolder{irun} ' - ' time{irun} '\n'])
    end


    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']))
    
    max_dist = 2;
    seed = rng;
    nCells = size(resp_cell{end,end,end},1);
    nTrials = size(stimCon_all,2);
    trial_n = zeros(nMaskCon,nStimCon,nMaskPhas);
    trialInd = cell(nMaskCon,nStimCon,nMaskPhas);
    for im = 1:nMaskCon
        ind_mask = find(maskCon_all == maskCons(im));
        for it = 1:nStimCon
            ind_stim = find(stimCon_all == stimCons(it));
            ind_sm = intersect(ind_mask,ind_stim);
            if it>1 & im>1
                for ip = 1:nMaskPhas
                    ind_phase = find(maskPhas_all == maskPhas(ip));
                    ind = intersect(ind_phase,ind_sm);
                    trialInd{im,it,ip} = ind;
                    trial_n(im,it,ip) = length(ind);
                end
            else
                trialInd{im,it,1} = ind_sm;
                trial_n(im,it,1) = length(ind_sm);
            end
        end
    end
    
    p_anova_all = nan(nCells,1);
    b_hat_all = nan(nCells,1); 
    amp_hat_all = nan(nCells,1); 
    per_hat_all = nan(nCells,1); 
    pha_hat_all = nan(nCells,1); 
    sse_all = nan(nCells,1); 
    R_square_all = nan(nCells,1);
    yfit_all = nan(nCells,length(0:1:359),1);
    
    eye_n_all = nan(1,nMaskPhas);
    phase_range = 0:1:359;
    %less than 2 deg
    resp_all = [];
    stim_all = [];
    resp_avg = cell(1,nMaskPhas);
    trN = zeros(1,nMaskPhas);
    for im = 2:nMaskCon
        for it = 2:nStimCon
            [memb ind_test] = ismember(trialInd{1,it,1},find(centroid_dist<max_dist));
            [memb ind_mask] = ismember(trialInd{im,1,1},find(centroid_dist<max_dist));
            test_avg = mean(resp_cell{1,it,1}(:,find(ind_test)),2);
            test_avg_rect = test_avg;
            test_avg_rect(find(test_avg<0)) = 0;
            mask_avg = mean(resp_cell{im,1,1}(:,find(ind_mask)),2);
            mask_avg_rect = mask_avg;
            mask_avg_rect(find(mask_avg<0)) = 0;
            for ip = 1:nMaskPhas
                [memb ind] = ismember(trialInd{im,it,ip},find(centroid_dist<max_dist));
                resp_all = [resp_all resp_cell{im,it,ip}(:,find(ind))];
                stim_all = [stim_all ip.*ones(size(resp_cell{im,it,ip}(1,find(ind))))];
                resp_avg{1,ip} = [resp_avg{1,ip} resp_cell{im,it,ip}(:,find(ind))];
                trN(ip) = length(find(ind));
            end
        end
    end

    resp_downsamp_rect = resp_all;
    resp_downsamp_rect(find(resp_all<0)) = 0;
    SI_all = (resp_downsamp_rect-(test_avg_rect+mask_avg_rect))./(resp_downsamp_rect+(test_avg_rect+mask_avg_rect));
    resp_avg_downsamp = nan(nCells,nMaskPhas);
    SI_all_avg = nan(nCells,nMaskPhas,2);
    for ip = 1:nMaskPhas
        resp_avg_downsamp(:,ip) = mean(resp_avg{1,ip},2);
        SI_all_avg(:,ip,1) = nanmean(SI_all(:,find(stim_all==ip)),2);
        SI_all_avg(:,ip,2) = nanstd(SI_all(:,find(stim_all==ip)),[],2)./sqrt(length(find(stim_all==ip)));
    end
    resp_avg_rect = resp_avg_downsamp;
    resp_avg_rect(find(resp_avg_downsamp<0)) = 0;
    SI_avg = (resp_avg_rect-(test_avg_rect+mask_avg_rect))./(resp_avg_rect+(test_avg_rect+mask_avg_rect));
    [eye_n edges bin] = histcounts(stim_all,[1:5]);
   
    figure;
    start = 1;
    n = 1;
    for iCell =1:nCells
        if start>25
            sgtitle([mouse ' ' date '- Thresh- Trials < ' num2str(max_dist) '  deg'])
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_SI_maxDist' num2str(max_dist) '_' num2str(n) '.pdf']), '-dpdf','-fillpage')
            figure;
            start = 1;
            n = n+1;
        end
        p_anova_all(iCell,1) = anova1(resp_all(iCell,:), stim_all,'off');
        if max(SI_avg(iCell,:),[],2)>min(SI_avg(iCell,:),[],2)
            [b_hat_all(iCell,1), amp_hat_all(iCell,1), per_hat_all(iCell,1),pha_hat_all(iCell,1),sse_all(iCell,1),R_square_all(iCell,1)] = sinefit(deg2rad(maskPhas(stim_all)),SI_all(iCell,:));
            subplot(5,5,start)
            scatter(maskPhas(stim_all),SI_all(iCell,:));
            hold on
            scatter(maskPhas,SI_avg(iCell,:))
            yfit_all(iCell,:,1) = b_hat_all(iCell,1)+amp_hat_all(iCell,1).*(sin(2*pi*deg2rad(phase_range)./per_hat_all(iCell,1) + 2.*pi/pha_hat_all(iCell,1)));
            plot(phase_range, yfit_all(iCell,:,1));
            title(['Rsq = ' num2str(chop(R_square_all(iCell,1),2)) '; p = ' num2str(chop(p_anova_all(iCell,1),2))])
        else
            b_hat_all(iCell,1) = max(SI_avg(iCell,:),[],2);
            amp_hat_all(iCell,1) = 0;
            per_hat_all(iCell,1) = NaN;
            pha_hat_all(iCell,1) = NaN;
            sse_all(iCell,1) = 0;
            R_square_all(iCell,1)= 0;
            yfit_all(iCell,:,1) = nan(length(phase_range),1);
        end
        start = start+1;
    end
    sgtitle([mouse ' ' date '- Thresh- Trials < ' num2str(max_dist) '  deg'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_SI_maxDist' num2str(max_dist) '_' num2str(n) '.pdf']), '-dpdf','-fillpage')

    p_anova_shuf = nan(nCells,1);
    b_hat_shuf = nan(nCells,1); 
    amp_hat_shuf = nan(nCells,1); 
    per_hat_shuf = nan(nCells,1); 
    pha_hat_shuf = nan(nCells,1); 
    sse_shuf = nan(nCells,1); 
    R_square_shuf = nan(nCells,1);
    yfit_shuf = nan(nCells,length(0:1:359),1);
    
    rng(seed);
    stim_all_shuf = stim_all(randperm(length(stim_all)));
    figure;
    start = 1;
    n = 1;
    for iCell = 1:nCells
        if start>25
            sgtitle([mouse ' ' date '- Thresh Shuffled- Trials < ' num2str(max_dist) '  deg'])
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_SI_maxDist' num2str(max_dist) '_' num2str(n) '_shuffled.pdf']), '-dpdf','-fillpage')
            figure;
            start = 1;
            n = n+1;
        end
        p_anova_shuf(iCell,1) = anova1(resp_all(iCell,:), stim_all_shuf,'off');
        if max(SI_avg(iCell,:),[],2)>min(SI_avg(iCell,:),[],2)
            [b_hat_shuf(iCell,1), amp_hat_shuf(iCell,1), per_hat_shuf(iCell,1),pha_hat_shuf(iCell,1),sse_shuf(iCell,1),R_square_shuf(iCell,1)] = sinefit(deg2rad(maskPhas(stim_all_shuf)),SI_all(iCell,:));
            subplot(5,5,start)
            scatter(maskPhas(stim_all_shuf),SI_all(iCell,:));
            hold on
            scatter(maskPhas,SI_avg(iCell,:))
            yfit_shuf(iCell,:,1) = b_hat_shuf(iCell,1)+amp_hat_shuf(iCell,1).*(sin(2*pi*deg2rad(phase_range)./per_hat_shuf(iCell,1) + 2.*pi/pha_hat_shuf(iCell,1)));
            plot(phase_range, yfit_shuf(iCell,:,1));
            title(['Rsq = ' num2str(chop(R_square_shuf(iCell,1),2)) '; p = ' num2str(chop(p_anova_shuf(iCell,1),2))])
        else
            b_hat_shuf(iCell,1) = max(SI_avg(iCell,:),[],2);
            amp_hat_shuf(iCell,1) = 0;
            per_hat_shuf(iCell,1) = NaN;
            pha_hat_shuf(iCell,1) = NaN;
            sse_shuf(iCell,1) = 0;
            R_square_shuf(iCell,1)= 0;
            yfit_shuf(iCell,:,1) = nan(length(phase_range),1);
        end
        start = start+1;
    end
    sgtitle([mouse ' ' date '- Thresh- Shuffled Trials < ' num2str(max_dist) '  deg'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_SI_maxDist' num2str(max_dist) '_' num2str(n) '_shuffled.pdf']), '-dpdf','-fillpage')

    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits.mat']),'trN', 'seed', 'yfit_all', 'b_hat_all', 'amp_hat_all', 'per_hat_all', 'pha_hat_all', 'sse_all', 'R_square_all', 'p_anova_all', 'yfit_shuf', 'b_hat_shuf', 'amp_hat_shuf', 'per_hat_shuf', 'pha_hat_shuf', 'sse_shuf', 'R_square_shuf', 'p_anova_shuf','trial_n', 'trialInd','SI_all_avg', 'max_dist')
    close all
end