close all; clear all; clc;
doRedChannel = 0;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
manuscriptDir = fullfile(LG_base, 'Manuscripts','2021','Barbera2021','FigureDataAndCode');
ds = ['CrossOriRandDirRandPhase_ExptList'];
eval(ds);

nexp = size(expt,2);
expt_ind = [];
resp_ind_all = [];
resp_ind_all_phase = [];
resp_tc_all = [];
resp_avg_all = [];
SI_avg_all = [];
amp_all = [];
amp_shuf_all = [];
yfit_allcells = [];
yfit_shuf_allcells = [];
anova_all = [];
pha_all = [];
b_all = [];
test_resp_all = [];
mask_resp_all = [];
plaid_resp_all = [];
phase_MI_all = [];
phase_MI_max_all = [];
phase_SI_all = [];
totCells = 0;
expUse = [];
mouse_list = [];

for iexp = 4:nexp
    if strcmp(expt(iexp).driver,'PV') & expt(iexp).sf == 0.05
        expUse = [expUse iexp];
        mouse = expt(iexp).mouse;
        mouse_list = strvcat(mouse_list, mouse);
        date = expt(iexp).date;
        fprintf([mouse ' ' date '\n'])
        ImgFolder = expt(iexp).copFolder;
        nrun = length(ImgFolder);
        run_str = catRunName(cell2mat(ImgFolder), nrun);

        clear resp_ind
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']));
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']));
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']));
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']),'npSub_tc');
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']),'centroid_dist');
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits.mat']));
        if strcmp(date,'210226') & strcmp(mouse,'i1340')
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
        end
        
        resp_ind_all = [resp_ind_all; resp_ind+totCells];
        resp_ind_all_phase = [resp_ind_all_phase; unique([resptest_ind; respmask_ind])+totCells];
        nCells = size(data_dfof_tc,2);
        totCells = totCells+nCells;
        expt_ind = [expt_ind; iexp.*ones(nCells,1)];
        expt(iexp).nP = nMaskPhas;
        postwin_frames = frame_rate.*6;
        resp_tc = nan(prewin_frames+postwin_frames,nCells,nMaskCon,nStimCon,nMaskPhas,2);
        resp_avg = nan(nCells,nMaskCon,nStimCon,nMaskPhas,2);
        SI_avg = nan(nCells,nMaskCon,nStimCon,nMaskPhas,2);
            
        if isempty(find(trN<4)) & length(resp_ind_all_phase)>4
            nOn = frame_rate.*4;
            data_trial = nan(prewin_frames+postwin_frames,nCells,nTrials);
            for iTrial = 1:nTrials-1
                data_trial(:,:,iTrial) = npSub_tc(cStimOn(iTrial)-prewin_frames:cStimOn(iTrial)+postwin_frames-1,:);
            end
            data_f = mean(data_trial(1:prewin_frames,:,:),1);
            data_df = data_trial-data_f;
            data_dfof = data_df./data_f;
            resp_win = prewin_frames+5:prewin_frames+nOn;
            resp_tc = zeros(prewin_frames+postwin_frames,nCells,nMaskCon,nStimCon,nMaskPhas,2);
            resp_avg = zeros(nCells,nMaskCon,nStimCon,nMaskPhas,2);
            resp_avg_rect = zeros(nCells,nMaskCon,nStimCon,nMaskPhas,2);
            resp_cell = cell(nMaskCon,nStimCon,nMaskPhas);
            resp_cell_allphase = cell(nMaskCon,nStimCon);
            SI_avg = zeros(nCells,nMaskCon,nStimCon,nMaskPhas,2);
            ind_n = zeros(nMaskCon,nStimCon,nMaskPhas);
            for im = 1:nMaskCon
                ind_m = find(maskCon_all == maskCons(im));
                for is = 1:nStimCon
                    ind_s = find(stimCon_all == stimCons(is));
                    resp_cell_allphase{im,is} = [];
                    if im>1 & is>1
                        for ip = 1:nMaskPhas
                            ind_p = find(maskPhas_all == maskPhas(ip));
                            ind_use = intersect(find(centroid_dist<2), intersect(ind_p,intersect(ind_m,ind_s)));
                            ind_n(im,is,ip) = length(ind_use);
                            resp_tc(:,:,im,is,ip,1) = nanmean(data_dfof(:,:,ind_use),3);
                            resp_tc(:,:,im,is,ip,2) = nanstd(data_dfof(:,:,ind_use),[],3)./sqrt(length(ind_use));
                            resp_cell{im,is,ip} = squeeze(nanmean(data_dfof(resp_win,:,ind_use),1));
                            if length(ind_use)>1
                                resp_cell_allphase{im,is} = [resp_cell_allphase{im,is} resp_cell{im,is,ip}];
                            else
                                resp_cell_allphase{im,is} = [resp_cell_allphase{im,is} resp_cell{im,is,ip}'];
                            end
                            resp_cell{im,is,ip}(find(resp_cell{im,is,ip}<0)) = 0;
                            resp_avg(:,im,is,ip,1) = nanmean(resp_cell{im,is,ip},2)./sqrt(length(ind_use));
                            resp_avg(:,im,is,ip,2) = nanstd(resp_cell{im,is,ip},[],2)./sqrt(length(ind_use));
                            SI_avg(:,im,is,ip,1) = squeeze(nanmean((resp_cell{im,is,ip}-(resp_avg_rect(:,1,is,1,1)+resp_avg_rect(:,im,1,1,1)))./...
                                (resp_cell{im,is,ip}+(resp_avg_rect(:,1,is,1,1)+resp_avg_rect(:,im,1,1,1))),2));
                            SI_avg(:,im,is,ip,2) = squeeze(nanstd((resp_cell{im,is,ip}-(resp_avg_rect(:,1,is,1,1)+resp_avg_rect(:,im,1,1,1)))./...
                                (resp_cell{im,is,ip}+(resp_avg_rect(:,1,is,1,1)+resp_avg_rect(:,im,1,1,1))),[],2))./sqrt(length(ind_use));
                        end
                    else
                        ind_use = intersect(ind_m,ind_s);
                        ind_n(im,is,1) = length(ind_use);
                        resp_tc(:,:,im,is,1,1) = nanmean(data_dfof(:,:,ind_use),3);
                        resp_tc(:,:,im,is,1,2) = nanstd(data_dfof(:,:,ind_use),[],3)./sqrt(length(ind_use));
                        resp_cell{im,is,1} = squeeze(nanmean(data_dfof(resp_win,:,ind_use),1));
                        resp_cell_allphase{im,is} = resp_cell{im,is,1};
                        resp_avg(:,im,is,1,1) = squeeze(nanmean(resp_cell{im,is,1},2));
                        resp_avg_rect(:,im,is,1,1) = resp_avg(:,im,is,1,1);
                        resp_avg_rect(find(resp_avg_rect(:,im,is,1,1)<0),im,is,1,1) = 0;
                        resp_avg(:,im,is,1,2) = squeeze(nanstd(resp_cell{im,is,1},[],2))./sqrt(size(resp_cell{im,is,1},2));
                    end
                end
            end

            test_resp = resp_avg_rect(:,1,end,1,1);
            mask_resp = resp_avg_rect(:,end,1,1,1);
            plaid_resp = nanmean(resp_cell_allphase{end,end},2);
            plaid_resp(find(plaid_resp<0)) = 0;

            stimSI = (abs(test_resp-mask_resp))./(test_resp+mask_resp);
            phase_SI_all = [phase_SI_all; stimSI];
            test_resp_all = [test_resp_all; test_resp];
            mask_resp_all = [mask_resp_all; mask_resp];
            plaid_resp_all = [plaid_resp_all; plaid_resp];
            plaidSI = (plaid_resp-(test_resp+mask_resp))./(plaid_resp+(test_resp+mask_resp));
            phase_MI_all = [phase_MI_all; plaidSI];
            phase_MI_max = (plaid_resp-(max([test_resp mask_resp],[],2)))./(plaid_resp+(max([test_resp mask_resp],[],2)));
            phase_MI_max_all = [phase_MI_max_all; phase_MI_max];

            resp_tc_all = cat(2,resp_tc_all,resp_tc);
            resp_avg_all = cat(1,resp_avg_all, resp_avg);
            SI_avg_all = cat(1,SI_avg_all, SI_avg);

            amp_all = [amp_all; amp_hat_all];
            amp_shuf_all = [amp_shuf_all; amp_hat_shuf];
            yfit_allcells = [yfit_allcells; yfit_all];
            yfit_shuf_allcells = [yfit_shuf_allcells; yfit_shuf];
            anova_all = [anova_all; p_anova_all];
            pha_all = [pha_all; pha_hat_all];
            b_all = [b_all; b_hat_all];
        else
            phase_SI_all = [phase_SI_all; nan(size(amp_hat_all))];
            test_resp_all = [test_resp_all; nan(size(amp_hat_all))];
            mask_resp_all = [mask_resp_all; nan(size(amp_hat_all))];
            plaid_resp_all = [plaid_resp_all; nan(size(amp_hat_all))];
            phase_MI_all = [phase_MI_all; nan(size(amp_hat_all))];
            phase_MI_max_all = [phase_MI_max_all; nan(size(amp_hat_all))];

            resp_tc_all = cat(2,resp_tc_all,resp_tc);
            resp_avg_all = cat(1,resp_avg_all, resp_avg);
            SI_avg_all = cat(1,SI_avg_all, SI_avg);

            amp_all = [amp_all; nan(size(amp_hat_all))];
            amp_shuf_all = [amp_shuf_all; nan(size(amp_hat_all))];
            yfit_allcells = [yfit_allcells; nan(size(amp_hat_all,1),360)];
            yfit_shuf_allcells = [yfit_shuf_allcells; nan(size(amp_hat_all,1),360)];
            anova_all = [anova_all; nan(size(amp_hat_all))];
            pha_all = [pha_all; nan(size(amp_hat_all))];
            b_all = [b_all; nan(size(amp_hat_all))];
        end

       
    end
end
tt= (-prewin_frames:postwin_frames-1).*(1000./frame_rate);
save(fullfile(manuscriptDir, 'Figure8_CalciumImaging_data.mat'),'anova_all','resp_ind_all_phase','amp_all','pha_all',...
'amp_shuf_all', 'b_all', 'phase_SI_all', 'phase_MI_all', 'phase_MI_max_all', 'mouse_list','expt','expUse','resp_tc_all','tt','expt_ind',...
'SI_avg_all','yfit_allcells','data_avg','maskPhas','nMaskPhas')