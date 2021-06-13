clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandPhase_ExptList';
eval(ds)
rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);
nanframes = zeros(1,nexp);

ntrials = nan(nexp,3,2);
SI_avg_all = [];
resp_avg_all = [];
resp_ind_all = [];
max_dir_all = [];
dir_resp_ind_all = [];
mouse_list = [];
pupil_avg = nan(nexp,2);
totCells = 0;
for iexp = 1:nexp
    if strcmp(expt(iexp).driver,'SLC')
        mouse = expt(iexp).mouse;
        date = expt(iexp).date;
        area = expt(iexp).img_loc{1};
        ImgFolder = expt(iexp).coFolder;
        time = expt(iexp).coTime;
        nrun = length(ImgFolder);
        run_str = catRunName(cell2mat(ImgFolder), nrun);
        mouse_list = [mouse_list; mouse];
        LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

        fprintf([mouse ' ' date '\n'])

        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']))
        
        ImgFolder = expt(iexp).dirFolder;
        nrun = length(ImgFolder);
        run_str = catRunName(cell2mat(ImgFolder), nrun);
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']))
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
        
        [max_val max_dir] = max(mean(data_dfof_dir(140:end,:,:),1),[],3);
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

        p = prctile(rad_stim,[25 75]);
        SI_avg = zeros(nCells,2);
        resp_avg  = zeros(nCells,2);
        for i = 1:2
            if i==1
                ind_use = find(rad_stim<p(i));
            else
                ind_use = find(rad_stim>p(i));
            end
            pupil_avg(iexp,i) = mean(rad_stim(ind_use));
            [memb ind_test] = ismember(trialInd{1,nStimCon,1},ind_use);
            [memb ind_mask] = ismember(trialInd{nMaskCon,1,1},ind_use);
            ntrials(iexp,1,i) = length(find(ind_test));
            ntrials(iexp,2,i) = length(find(ind_mask));
            test_avg = mean(resp_cell{1,nStimCon,1}(:,find(ind_test)),2);
            test_avg_rect = test_avg;
            test_avg_rect(find(test_avg<0)) = 0;
            mask_avg = mean(resp_cell{nMaskCon,1,1}(:,find(ind_mask)),2);
            mask_avg_rect = mask_avg;
            mask_avg_rect(find(mask_avg<0)) = 0;
            [memb ind] = ismember(trialInd{nMaskCon,nStimCon,1},ind_use);
            ntrials(iexp,3,i) = length(find(ind));
            plaid_avg = mean(resp_cell{nMaskCon,nStimCon,1}(:,find(ind)),2);
            plaid_avg_rect = plaid_avg;
            plaid_avg_rect(find(plaid_avg<0)) = 0;
            resp_avg(:,i) = max([test_avg mask_avg],[],2);
            SI_avg(:,i) = (plaid_avg_rect-(test_avg_rect+mask_avg_rect))./(plaid_avg_rect+(test_avg_rect+mask_avg_rect));
        end
        SI_avg_all = [SI_avg_all; SI_avg];   
        resp_avg_all = [resp_avg_all; resp_avg];   
        resp_ind_all = [resp_ind_all; resp_ind+totCells];
        
        max_dir_all = [max_dir_all max_dir];
        dir_resp_ind_all = [dir_resp_ind_all dirresp_ind+totCells];
        totCells = totCells+nCells;
    end
end

dir_pref_ind = intersect(dir_resp_ind_all, [find(max_dir_all==1) find(max_dir_all==5)]);
figure;
subplot(2,2,3)
plot([1 2], pupil_avg,'-ok')
pupil_avg_avg = [nanmean(pupil_avg',2) nanstd(pupil_avg',[],2)./sqrt(sum(~isnan(pupil_avg'),2))];
ylabel('Pupil radius (mm)')
ylim([0 0.5])
set(gca,'XTick',1:2,'XTickLabel',[25 75])
xlabel('Percentile')
[h_SI p_SI] = ttest(SI_avg_all(resp_ind_all,1),SI_avg_all(resp_ind_all,2));
[h_resp p_resp] = ttest(resp_avg_all(resp_ind_all,1),resp_avg_all(resp_ind_all,2));
resp_avg = zeros(2,2);
SI_avg = zeros(2,2);
for i = 1:2
    subplot(2,2,1)
    cdfplot(SI_avg_all(resp_ind_all,i))
    SI_avg(i,:) = [nanmean(SI_avg_all(resp_ind_all,i),1)  nanstd(SI_avg_all(resp_ind_all,i),[],1)./sqrt(length(resp_ind_all))];
    hold on 
    subplot(2,2,2)
    cdfplot(SI_avg_all(intersect(dir_pref_ind,resp_ind_all),i))
    hold on 
    subplot(2,2,4)
    cdfplot(resp_avg_all(resp_ind_all,i))
    resp_avg(i,:) = [nanmean(resp_avg_all(resp_ind_all,i),1)  nanstd(resp_avg_all(resp_ind_all,i),[],1)./sqrt(length(resp_ind_all))];
    hold on 
end
subplot(2,2,1)
legend({'Small','Large'},'location','southeast')
title('All responsive cells')
xlabel('Masking index')
ylabel('Fraction of cells')
subplot(2,2,2)
title('Test/Mask Pref')
xlabel('Masking index')
ylabel('Fraction of cells')
subplot(2,2,4)
xlabel('Test/Mask response (dF/F)')
ylabel('Fraction of cells')

fnout = fullfile(LG_base, 'Analysis\2P\CrossOri\CrossOri_Figures');
print(fullfile(fnout,'CrossOriRandPhase_arousal.pdf'),'-dpdf')