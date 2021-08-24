clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandPhaseAdapt_ExptList';
eval(ds)
rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);
nanframes = zeros(1,nexp);
max_dist = 5;

for iexp = 5:6
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


    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']))
    
    max_dist = 5;
    seed = rng;
    nCells = size(resp_cell_noadapt{end,end,end},1);
    nTrials = size(stimCon_all,2);

    p_anova_noadapt = nan(nCells,1);
    b_hat_noadapt = nan(nCells,1); 
    amp_hat_noadapt = nan(nCells,1); 
    per_hat_noadapt = nan(nCells,1); 
    pha_hat_noadapt = nan(nCells,1); 
    sse_noadapt = nan(nCells,1); 
    R_square_noadapt = nan(nCells,1);
    yfit_noadapt = nan(nCells,length(0:1:359),1);
    
    eye_n_all = nan(1,nMaskPhas);
    phase_range = 0:1:359;
    %less than 4 deg
    resp_all = [];
    stim_all = [];
    resp_avg = cell(1,nMaskPhas);
    trN_noadapt = zeros(1,nMaskPhas);
    for im = 2:nMaskCon
        for it = 2:nStimCon
            [memb ind_test] = ismember(trialInd_noadapt{1,it,1},find(centroid_dist<max_dist));
            [memb ind_mask] = ismember(trialInd_noadapt{im,1,1},find(centroid_dist<max_dist));
            test_avg = mean(resp_cell_noadapt{1,it,1}(:,find(ind_test)),2);
            test_avg_rect = test_avg;
            test_avg_rect(find(test_avg<0)) = 0;
            mask_avg = mean(resp_cell_noadapt{im,1,1}(:,find(ind_mask)),2);
            mask_avg_rect = mask_avg;
            mask_avg_rect(find(mask_avg<0)) = 0;
            for ip = 1:nMaskPhas
                [memb ind] = ismember(trialInd_noadapt{im,it,ip},find(centroid_dist<max_dist));
                resp_all = [resp_all resp_cell_noadapt{im,it,ip}(:,find(ind))];
                stim_all = [stim_all ip.*ones(size(resp_cell_noadapt{im,it,ip}(1,find(ind))))];
                resp_avg{1,ip} = [resp_avg{1,ip} resp_cell_noadapt{im,it,ip}(:,find(ind))];
                trN_noadapt(ip) = length(find(ind));
            end
        end
    end
    trN_noadapt = [length(find(ind_test)) length(find(ind_mask)) trN_noadapt];
    
    resp_downsamp_rect = resp_all;
    resp_downsamp_rect(find(resp_all<0)) = 0;
    SI_all = (resp_downsamp_rect-(test_avg_rect+mask_avg_rect))./(resp_downsamp_rect+(test_avg_rect+mask_avg_rect));
    resp_avg_noadapt = nan(nCells,nMaskPhas);
    SI_all_avg = nan(nCells,nMaskPhas,2);
    for ip = 1:nMaskPhas
        resp_avg_noadapt(:,ip) = mean(resp_avg{1,ip},2);
        SI_all_avg(:,ip,1) = nanmean(SI_all(:,find(stim_all==ip)),2);
        SI_all_avg(:,ip,2) = nanstd(SI_all(:,find(stim_all==ip)),[],2)./sqrt(length(find(stim_all==ip)));
    end
    resp_avg_rect = resp_avg_noadapt;
    resp_avg_rect(find(resp_avg_noadapt<0)) = 0;
    SI_avg = (resp_avg_rect-(test_avg_rect+mask_avg_rect))./(resp_avg_rect+(test_avg_rect+mask_avg_rect));
    [eye_n edges bin] = histcounts(stim_all,[1:5]);
   
    figure;
    start = 1;
    n = 1;
    for iCell =1:nCells
        if start>25
            sgtitle([mouse ' ' date '- No Adapt- Trials < ' num2str(max_dist) '  deg'])
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_SI_lessThan' num2str(max_dist) 'degEyeMvmt' num2str(n) '_noAdapt.pdf']), '-dpdf','-fillpage')
            figure;
            start = 1;
            n = n+1;
        end
        p_anova_noadapt(iCell,1) = anova1(resp_all(iCell,:), stim_all,'off');
        if max(SI_avg(iCell,:),[],2)>min(SI_avg(iCell,:),[],2)
            [b_hat_noadapt(iCell,1), amp_hat_noadapt(iCell,1), per_hat_noadapt(iCell,1),pha_hat_noadapt(iCell,1),sse_noadapt(iCell,1),R_square_noadapt(iCell,1)] = sinefit(deg2rad(maskPhas(stim_all)),SI_all(iCell,:));
            subplot(5,5,start)
            scatter(maskPhas(stim_all),SI_all(iCell,:));
            hold on
            scatter(maskPhas,SI_avg(iCell,:))
            yfit_noadapt(iCell,:,1) = b_hat_noadapt(iCell,1)+amp_hat_noadapt(iCell,1).*(sin(2*pi*deg2rad(phase_range)./per_hat_noadapt(iCell,1) + 2.*pi/pha_hat_noadapt(iCell,1)));
            plot(phase_range, yfit_noadapt(iCell,:,1));
            title(['Rsq = ' num2str(chop(R_square_noadapt(iCell,1),2)) '; p = ' num2str(chop(p_anova_noadapt(iCell,1),2))])
        else
            b_hat_noadapt(iCell,1) = max(SI_avg(iCell,:),[],2);
            amp_hat_noadapt(iCell,1) = 0;
            per_hat_noadapt(iCell,1) = NaN;
            pha_hat_noadapt(iCell,1) = NaN;
            sse_noadapt(iCell,1) = 0;
            R_square_noadapt(iCell,1)= 0;
            yfit_noadapt(iCell,:,1) = nan(length(phase_range),1);
        end
        start = start+1;
    end
    sgtitle([mouse ' ' date '- No Adapt- Trials < ' num2str(max_dist) '  deg'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_SI_lessThan' num2str(max_dist) 'degEyeMvmt' num2str(n) '_noAdapt.pdf']), '-dpdf','-fillpage')

    p_anova_singadapt = nan(nCells,1);
    b_hat_singadapt = nan(nCells,1); 
    amp_hat_singadapt = nan(nCells,1); 
    per_hat_singadapt = nan(nCells,1); 
    pha_hat_singadapt = nan(nCells,1); 
    sse_singadapt = nan(nCells,1); 
    R_square_singadapt = nan(nCells,1);
    yfit_singadapt = nan(nCells,length(0:1:359),1);
    
    eye_n_all = nan(1,nMaskPhas);
    phase_range = 0:1:359;
    %less than 4 deg
    resp_all = [];
    stim_all = [];
    resp_avg = cell(1,nMaskPhas);
    trN_singadapt = zeros(1,nMaskPhas);
    for im = 2:nMaskCon
        for it = 2:nStimCon
            [memb ind_test] = ismember(trialInd_singadapt{1,it,1},find(centroid_dist<max_dist));
            [memb ind_mask] = ismember(trialInd_singadapt{im,1,1},find(centroid_dist<max_dist));
            test_avg = mean(resp_cell_singadapt{1,it,1}(:,find(ind_test)),2);
            test_avg_rect = test_avg;
            test_avg_rect(find(test_avg<0)) = 0;
            mask_avg = mean(resp_cell_singadapt{im,1,1}(:,find(ind_mask)),2);
            mask_avg_rect = mask_avg;
            mask_avg_rect(find(mask_avg<0)) = 0;
            for ip = 1:nMaskPhas
                [memb ind] = ismember(trialInd_singadapt{im,it,ip},find(centroid_dist<max_dist));
                resp_all = [resp_all resp_cell_singadapt{im,it,ip}(:,find(ind))];
                stim_all = [stim_all ip.*ones(size(resp_cell_singadapt{im,it,ip}(1,find(ind))))];
                resp_avg{1,ip} = [resp_avg{1,ip} resp_cell_singadapt{im,it,ip}(:,find(ind))];
                trN_singadapt(ip) = length(find(ind));
            end
        end
    end
    trN_singadapt = [length(find(ind_test)) length(find(ind_mask)) trN_singadapt];
    resp_downsamp_rect = resp_all;
    resp_downsamp_rect(find(resp_all<0)) = 0;
    SI_all = (resp_downsamp_rect-(test_avg_rect+mask_avg_rect))./(resp_downsamp_rect+(test_avg_rect+mask_avg_rect));
    resp_avg_singadapt = nan(nCells,nMaskPhas);
    SI_all_avg = nan(nCells,nMaskPhas,2);
    for ip = 1:nMaskPhas
        resp_avg_singadapt(:,ip) = mean(resp_avg{1,ip},2);
        SI_all_avg(:,ip,1) = nanmean(SI_all(:,find(stim_all==ip)),2);
        SI_all_avg(:,ip,2) = nanstd(SI_all(:,find(stim_all==ip)),[],2)./sqrt(length(find(stim_all==ip)));
    end
    resp_avg_rect = resp_avg_singadapt;
    resp_avg_rect(find(resp_avg_singadapt<0)) = 0;
    SI_avg = (resp_avg_rect-(test_avg_rect+mask_avg_rect))./(resp_avg_rect+(test_avg_rect+mask_avg_rect));
    [eye_n edges bin] = histcounts(stim_all,[1:5]);
   
    figure;
    start = 1;
    n = 1;
    for iCell =1:nCells
        if start>25
            sgtitle([mouse ' ' date '- Sing Adapt - Trials < ' num2str(max_dist) '  deg'])
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_SI_lessThan' num2str(max_dist) 'degEyeMvmt' num2str(n) '_singAdapt.pdf']), '-dpdf','-fillpage')
            figure;
            start = 1;
            n = n+1;
        end
        p_anova_singadapt(iCell,1) = anova1(resp_all(iCell,:), stim_all,'off');
        if max(SI_avg(iCell,:),[],2)>min(SI_avg(iCell,:),[],2)
            [b_hat_singadapt(iCell,1), amp_hat_singadapt(iCell,1), per_hat_singadapt(iCell,1),pha_hat_singadapt(iCell,1),sse_singadapt(iCell,1),R_square_singadapt(iCell,1)] = sinefit(deg2rad(maskPhas(stim_all)),SI_all(iCell,:));
            subplot(5,5,start)
            scatter(maskPhas(stim_all),SI_all(iCell,:));
            hold on
            scatter(maskPhas,SI_avg(iCell,:))
            yfit_singadapt(iCell,:,1) = b_hat_singadapt(iCell,1)+amp_hat_singadapt(iCell,1).*(sin(2*pi*deg2rad(phase_range)./per_hat_singadapt(iCell,1) + 2.*pi/pha_hat_singadapt(iCell,1)));
            plot(phase_range, yfit_singadapt(iCell,:,1));
            title(['Rsq = ' num2str(chop(R_square_singadapt(iCell,1),2)) '; p = ' num2str(chop(p_anova_singadapt(iCell,1),2))])
        else
            b_hat_singadapt(iCell,1) = max(SI_avg(iCell,:),[],2);
            amp_hat_singadapt(iCell,1) = 0;
            per_hat_singadapt(iCell,1) = NaN;
            pha_hat_singadapt(iCell,1) = NaN;
            sse_singadapt(iCell,1) = 0;
            R_square_singadapt(iCell,1)= 0;
            yfit_singadapt(iCell,:,1) = nan(length(phase_range),1);
        end
        start = start+1;
    end
    sgtitle([mouse ' ' date '- Sing Adapt- Trials < ' num2str(max_dist) '  deg'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_SI_lessThan' num2str(max_dist) 'degEyeMvmt' num2str(n) '_singAdapt.pdf']), '-dpdf','-fillpage')

    p_anova_shuf_noadapt = nan(nCells,1);
    b_hat_shuf_noadapt = nan(nCells,1); 
    amp_hat_shuf_noadapt = nan(nCells,1); 
    per_hat_shuf_noadapt = nan(nCells,1); 
    pha_hat_shuf_noadapt = nan(nCells,1); 
    sse_shuf_noadapt = nan(nCells,1); 
    R_square_shuf_noadapt = nan(nCells,1);
    yfit_shuf_noadapt = nan(nCells,length(0:1:359),1);
    
    rng(seed);
    stim_all_shuf = stim_all(randperm(length(stim_all)));
    figure;
    start = 1;
    n = 1;
    for iCell = 1:nCells
        if start>25
            sgtitle([mouse ' ' date '- No Adapt Shuffled- Trials < ' num2str(max_dist) '  deg'])
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_SI_lessThan' num2str(max_dist) 'degEyeMvmt' num2str(n) '_NoAdapt_shuffled.pdf']), '-dpdf','-fillpage')
            figure;
            start = 1;
            n = n+1;
        end
        p_anova_shuf_noadapt(iCell,1) = anova1(resp_all(iCell,:), stim_all_shuf,'off');
        if max(SI_avg(iCell,:),[],2)>min(SI_avg(iCell,:),[],2)
            [b_hat_shuf_noadapt(iCell,1), amp_hat_shuf_noadapt(iCell,1), per_hat_shuf_noadapt(iCell,1),pha_hat_shuf_noadapt(iCell,1),sse_shuf_noadapt(iCell,1),R_square_shuf_noadapt(iCell,1)] = sinefit(deg2rad(maskPhas(stim_all_shuf)),SI_all(iCell,:));
            subplot(5,5,start)
            scatter(maskPhas(stim_all_shuf),SI_all(iCell,:));
            hold on
            scatter(maskPhas,SI_avg(iCell,:))
            yfit_shuf_noadapt(iCell,:,1) = b_hat_shuf_noadapt(iCell,1)+amp_hat_shuf_noadapt(iCell,1).*(sin(2*pi*deg2rad(phase_range)./per_hat_shuf_noadapt(iCell,1) + 2.*pi/pha_hat_shuf_noadapt(iCell,1)));
            plot(phase_range, yfit_shuf_noadapt(iCell,:,1));
            title(['Rsq = ' num2str(chop(R_square_shuf_noadapt(iCell,1),2)) '; p = ' num2str(chop(p_anova_shuf_noadapt(iCell,1),2))])
        else
            b_hat_shuf_noadapt(iCell,1) = max(SI_avg(iCell,:),[],2);
            amp_hat_shuf_noadapt(iCell,1) = 0;
            per_hat_shuf_noadapt(iCell,1) = NaN;
            pha_hat_shuf_noadapt(iCell,1) = NaN;
            sse_shuf_noadapt(iCell,1) = 0;
            R_square_shuf_noadapt(iCell,1)= 0;
            yfit_shuf_noadapt(iCell,:,1) = nan(length(phase_range),1);
        end
        start = start+1;
    end
    sgtitle([mouse ' ' date '- No Adapt- Shuffled Trials < ' num2str(max_dist) '  deg'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_SI_lessThan' num2str(max_dist) 'degEyeMvmt' num2str(n) '_NoAdapt_shuffled.pdf']), '-dpdf','-fillpage')

    p_anova_shuf_singadapt = nan(nCells,1);
    b_hat_shuf_singadapt = nan(nCells,1); 
    amp_hat_shuf_singadapt = nan(nCells,1); 
    per_hat_shuf_singadapt = nan(nCells,1); 
    pha_hat_shuf_singadapt = nan(nCells,1); 
    sse_shuf_singadapt = nan(nCells,1); 
    R_square_shuf_singadapt = nan(nCells,1);
    yfit_shuf_singadapt = nan(nCells,length(0:1:359),1);
    
    rng(seed);
    stim_all_shuf = stim_all(randperm(length(stim_all)));
    figure;
    start = 1;
    n = 1;
    for iCell = 1:nCells
        if start>25
            sgtitle([mouse ' ' date '- Sing Adapt Shuffled- Trials < ' num2str(max_dist) '  deg'])
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_SI_lessThan' num2str(max_dist) 'degEyeMvmt' num2str(n) '_SingAdapt_shuffled.pdf']), '-dpdf','-fillpage')
            figure;
            start = 1;
            n = n+1;
        end
        p_anova_shuf_singadapt(iCell,1) = anova1(resp_all(iCell,:), stim_all_shuf,'off');
        if max(SI_avg(iCell,:),[],2)>min(SI_avg(iCell,:),[],2)
            [b_hat_shuf_singadapt(iCell,1), amp_hat_shuf_singadapt(iCell,1), per_hat_shuf_singadapt(iCell,1),pha_hat_shuf_singadapt(iCell,1),sse_shuf_singadapt(iCell,1),R_square_shuf_singadapt(iCell,1)] = sinefit(deg2rad(maskPhas(stim_all_shuf)),SI_all(iCell,:));
            subplot(5,5,start)
            scatter(maskPhas(stim_all_shuf),SI_all(iCell,:));
            hold on
            scatter(maskPhas,SI_avg(iCell,:))
            yfit_shuf_singadapt(iCell,:,1) = b_hat_shuf_singadapt(iCell,1)+amp_hat_shuf_singadapt(iCell,1).*(sin(2*pi*deg2rad(phase_range)./per_hat_shuf_singadapt(iCell,1) + 2.*pi/pha_hat_shuf_singadapt(iCell,1)));
            plot(phase_range, yfit_shuf_singadapt(iCell,:,1));
            title(['Rsq = ' num2str(chop(R_square_shuf_singadapt(iCell,1),2)) '; p = ' num2str(chop(p_anova_shuf_singadapt(iCell,1),2))])
        else
            b_hat_shuf_singadapt(iCell,1) = max(SI_avg(iCell,:),[],2);
            amp_hat_shuf_singadapt(iCell,1) = 0;
            per_hat_shuf_singadapt(iCell,1) = NaN;
            pha_hat_shuf_singadapt(iCell,1) = NaN;
            sse_shuf_singadapt(iCell,1) = 0;
            R_square_shuf_singadapt(iCell,1)= 0;
            yfit_shuf_singadapt(iCell,:,1) = nan(length(phase_range),1);
        end
        start = start+1;
    end
    sgtitle([mouse ' ' date '- Sing Adapt- Shuffled Trials < ' num2str(max_dist) '  deg'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_SI_lessThan' num2str(max_dist) 'degEyeMvmt' num2str(n) '_SingAdapt_shuffled.pdf']), '-dpdf','-fillpage')
   
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_' num2str(max_dist) '_phaseFits.mat']),'trN_noadapt','trN_singadapt', 'seed', 'yfit_noadapt', 'b_hat_noadapt', 'amp_hat_noadapt', 'per_hat_noadapt', 'pha_hat_noadapt', 'sse_noadapt', 'R_square_noadapt',  'p_anova_noadapt', 'yfit_singadapt', 'b_hat_singadapt', 'amp_hat_singadapt', 'per_hat_singadapt', 'pha_hat_singadapt', 'sse_singadapt', 'R_square_singadapt', 'p_anova_singadapt', 'yfit_shuf_noadapt', 'b_hat_shuf_noadapt', 'amp_hat_shuf_noadapt', 'per_hat_shuf_noadapt', 'pha_hat_shuf_noadapt', 'sse_shuf_noadapt', 'R_square_shuf_noadapt', 'p_anova_shuf_noadapt','yfit_shuf_singadapt', 'b_hat_shuf_singadapt', 'amp_hat_shuf_singadapt', 'per_hat_shuf_singadapt', 'pha_hat_shuf_singadapt', 'sse_shuf_singadapt', 'R_square_shuf_singadapt', 'p_anova_shuf_singadapt','trialsperstim_noadapt','trialsperstim_singadapt', 'trialInd_noadapt','trialInd_singadapt','SI_all_avg', 'max_dist','resp_avg_noadapt','resp_avg_singadapt')
    close all
end