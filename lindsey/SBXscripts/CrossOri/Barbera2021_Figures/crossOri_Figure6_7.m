close all; clear all; clc;
doRedChannel = 0;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir_F6 = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'CrossOri_Figures', 'CrossOri_Figure6');
summaryDir_F7 = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'CrossOri_Figures', 'CrossOri_Figure7');
summaryDir_GA = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'CrossOri_Figures', 'CrossOri_GA');

%all experiments SF = 0.05
ds = ['CrossOriRandPhase_ExptList'];
eval(ds);
expt1 = expt;
clear expt
ds = ['CrossOriRandPhase_lowSF_ExptList'];
eval(ds);
expt2 = expt;
clear expt
expt = combineStructures(expt1,expt2);

%load all experiments
nexp = size(expt,2);
expt_ind = [];
nP_ind = [];
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
OSI_all = [];
DSI_all = [];
peak_dir_all = [];
max_dir_all = [];
f1_all = [];
f2_all = [];
f2overf1_all = [];
totCells = 0;
expUse = [];
OSIexp = [];
mouse_list = [];
for iexp = 1:nexp
    if strcmp(expt(iexp).driver,'SLC') || strcmp(expt(iexp).driver,'EMX')
        expUse = [expUse iexp];
        mouse = expt(iexp).mouse;
        mouse_list = strvcat(mouse_list, mouse);
        date = expt(iexp).date;
        fprintf([mouse ' ' date '\n'])
        ImgFolder = expt(iexp).coFolder;
        nrun = length(ImgFolder);
        run_str = catRunName(cell2mat(ImgFolder), nrun);

        clear resp_ind
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']));
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']));
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']));
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']),'npSub_tc');
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']),'centroid_dist');
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits.mat']));

        resp_ind_all = [resp_ind_all; resp_ind+totCells];
        resp_ind_all_phase = [resp_ind_all_phase; unique([resptest_ind; respmask_ind])+totCells];

        nOn = frame_rate.*4;
        postwin_frames = frame_rate.*6;
        nCells = size(data_dfof_tc,2);
        totCells = totCells+nCells;
        expt_ind = [expt_ind; iexp.*ones(nCells,1)];
        expt(iexp).nP = nMaskPhas;
        nP_ind = [nP_ind; nMaskPhas.*ones(nCells,1)];

        data_trial = nan(prewin_frames+postwin_frames,nCells,nTrials);
        for iTrial = 1:nTrials-1
            data_trial(:,:,iTrial) = npSub_tc(cStimOn(iTrial)-prewin_frames:cStimOn(iTrial)+postwin_frames-1,:);
        end
        data_f = mean(data_trial(1:prewin_frames,:,:),1);
        data_df = data_trial-data_f;
        data_dfof = data_df./data_f;
        resp_win = prewin_frames+5:prewin_frames+nOn;
        resp_tc = zeros(prewin_frames+postwin_frames,nCells,nMaskCon,nStimCon,nMaskPhas,1);
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

        if nMaskPhas == 8
            resp_tc_all = cat(2,resp_tc_all,resp_tc);
            resp_avg_all = cat(1,resp_avg_all, resp_avg);
            SI_avg_all = cat(1,SI_avg_all, SI_avg);
        end
        amp_all = [amp_all; amp_hat_all];
        amp_shuf_all = [amp_shuf_all; amp_hat_shuf];
        yfit_allcells = [yfit_allcells; yfit_all];
        yfit_shuf_allcells = [yfit_shuf_allcells; yfit_shuf];
        anova_all = [anova_all; p_anova_all];
        pha_all = [pha_all; pha_hat_all];
        b_all = [b_all; b_hat_all];

        if nMaskPhas == 8
            ImgFolder = expt1(iexp).dirFolder;
            nrun = length(ImgFolder);
            run_str = catRunName(cell2mat(ImgFolder), nrun);
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']))
            [max_val max_ind] = max(mean(data_dfof_dir(resp_win,:,:),1),[],3);
            max_dir_all = [max_dir_all; max_ind'];
            OSI_all = [OSI_all; OSI_resp];
            DSI_all = [DSI_all; DSI_resp];
            peak_dir_all = [peak_dir_all; max_ind_dir];
            OSIexp = [OSIexp iexp];
        else
            OSI_all = [OSI_all; nan(size(amp_hat_all))];
            DSI_all = [DSI_all; nan(size(amp_hat_all))];
            peak_dir_all = [peak_dir_all; nan(size(amp_hat_all))];
            max_dir_all = [max_dir_all; nan(size(amp_hat_all))];
        end

        if nMaskPhas == 4
            if ~isempty(expt2(iexp-size(expt1,2)).prFolder)
                ImgFolder = expt2(iexp-size(expt1,2)).prFolder;
                nrun = length(ImgFolder);
                run_str = catRunName(cell2mat(ImgFolder), nrun);
                load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_f1f2.mat']))
                f1_all = [f1_all; f1'];
                f2_all = [f2_all; f2'];
                f2overf1_all = [f2overf1_all; f2overf1'];
            else
                f1_all = [f1_all; nan(size(amp_hat_all))];
                f2_all = [f2_all; nan(size(amp_hat_all))];
                f2overf1_all = [f2overf1_all; nan(size(amp_hat_all))];
            end
        else
            f1_all = [f1_all; nan(size(amp_hat_all))];
            f2_all = [f2_all; nan(size(amp_hat_all))];
            f2overf1_all = [f2overf1_all; nan(size(amp_hat_all))];
        end
    else
        expt(iexp).nP = NaN;
    end
end
stim_OSI_all = OSI_all;
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
save(fullfile(summaryDir,['randPhase_Summary_V1_SLC.mat']),'anova_all','stim_OSI_all','resp_ind_all_phase','amp_all', 'amp_shuf_all', 'b_all', 'phase_SI_all', 'phase_MI_all', 'phase_MI_max_all', 'mouse_list')
%%
respN = length(resp_ind_all_phase);
modN = sum(anova_all(resp_ind_all_phase)<0.05);
medAmp = median(amp_all(resp_ind_all_phase));
medShufAmp = median(amp_shuf_all(resp_ind_all_phase));
avgB = mean(b_all(resp_ind_all_phase));
semB = std(b_all(resp_ind_all_phase))./sqrt(respN);
[p_amp h_amp] = signrank(amp_all(resp_ind_all_phase),amp_shuf_all(resp_ind_all_phase),'tail','right');
figure;
subplot(2,2,1)
histogram(amp_all(resp_ind_all_phase),[0:0.05:1])
hold on
histogram(amp_shuf_all(resp_ind_all_phase),[0:0.05:1])
xlabel('Sine amplitude')
ylabel('Number of cells')
xlim([0 1])
legend({'data','shuffle'},'location','northeast')
mouse_cell = {expt.mouse};
mice  = unique(mouse_cell(expUse));
nmice = length(mice);
title([num2str(nexp) ' expts; ' num2str(nmice) ' mice; ' num2str(length(resp_ind_all_phase)) ' cells'])
subplot(2,2,2)
histogram(b_all(resp_ind_all_phase),[-1:0.1:1])
xlabel('Sine baseline')
ylabel('Number of cells')
xlim([-1 1])
subplot(2,2,3)
pha_all_circ = wrapTo2Pi(pha_all);
ind = intersect(resp_ind_all_phase, find(anova_all<0.05));
histogram(rad2deg(pha_all_circ(ind)),[0:45:360])
xlabel('Preferred phase (deg)')
ylabel('Number of cells')
ylim([0 100])
title(['n = ' num2str(length(ind)) ' cells'])
print(fullfile(summaryDir_F6, 'Figure6_ModulationSummary.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F6, 'Figure6_ModulationSummary.fig'))

figure
subplot(2,2,1)
histogram(amp_all(intersect(find(nP_ind==8),resp_ind_all_phase)),[0:0.05:1])
hold on
histogram(amp_shuf_all(intersect(find(nP_ind==8),resp_ind_all_phase)),[0:0.05:1])
xlabel('Sine amplitude')
ylabel('Number of cells')
xlim([0 1])
allmice_8phase = mouse_cell(find(cell2mat({expt.nP})==8));
nexp_8phase = length(allmice_8phase);
mice_8phase  = unique(allmice_8phase);
nmice_8phase = length(mice_8phase);
title(['8 Phase- ' num2str(nexp_8phase) ' expts; ' num2str(nmice_8phase) ' mice; ' num2str(length(intersect(find(nP_ind==8),resp_ind_all_phase))) ' cells'])
subplot(2,2,2)
histogram(amp_all(intersect(find(nP_ind==4),resp_ind_all_phase)),[0:0.05:1])
hold on
histogram(amp_shuf_all(intersect(find(nP_ind==4),resp_ind_all_phase)),[0:0.05:1])
xlabel('Sine amplitude')
ylabel('Number of cells')
xlim([0 1])
allmice_4phase = mouse_cell(find(cell2mat({expt.nP})==4));
nexp_4phase = length(allmice_4phase);
mice_4phase  = unique(allmice_4phase);
nmice_4phase = length(mice_4phase);
title(['4 Phase- ' num2str(nexp_4phase) ' expts; ' num2str(nmice_4phase) ' mice; ' num2str(length(intersect(find(nP_ind==4),resp_ind_all_phase))) ' cells'])
subplot(2,2,3)
ind = intersect(find(nP_ind==8),intersect(resp_ind_all_phase, find(anova_all<0.05)));
histogram(rad2deg(pha_all_circ(ind)),[0:45:360])
ylim([0 40])
xlabel('Preferred phase (deg)')
ylabel('Number of cells')
title(['n = ' num2str(length(ind)) ' cells'])
subplot(2,2,4)
ind = intersect(find(nP_ind==4),intersect(resp_ind_all_phase, find(anova_all<0.05)));
histogram(rad2deg(pha_all_circ(ind)),[0:45:360])
ylim([0 60])
xlabel('Preferred phase (deg)')
ylabel('Number of cells')
title(['n = ' num2str(length(ind)) ' cells'])
print(fullfile(summaryDir_F6, 'Figure6_ModulationSummary_byDataset.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F6, 'Figure6_ModulationSummary_byDataset.fig'))


ind_h = find(amp_all-amp_shuf_all > 0);
ind_l = find(amp_all-amp_shuf_all < 0);
tt= (-prewin_frames:postwin_frames-1).*(1000./frame_rate);
nStimCon = 2;
nMaskCon = 2;
nMaskPhas = 8;
maskPhas = [0:45:315];
figure;
i = 5;
subplot(3,3,1)
shadedErrorBar(tt', resp_tc_all(:,ind_h(i),1,nStimCon,1,1), resp_tc_all(:,ind_h(i),1,nStimCon,1,2),'lineprops',{['b-'], 'markerfacecolor', 'b'});
hold on
shadedErrorBar(tt, resp_tc_all(:,ind_h(i),nMaskCon,1,1,1), resp_tc_all(:,ind_h(i),nMaskCon,1,1,2),'lineprops',{['r-'], 'markerfacecolor','r'});
xlabel('Time (ms)')
ylabel('dF/F')
title('Stim/Mask')
ylim([-0.1 1.2])
for ip = 1:nMaskPhas
    subplot(3,3,ip+1)
    shadedErrorBar(tt, resp_tc_all(:,ind_h(i),nMaskCon,nStimCon,ip,1), resp_tc_all(:,ind_h(i),nMaskCon,nStimCon,ip,2));
    xlabel('Time (ms)')
    ylabel('dF/F')
    title(['Phase ' num2str(maskPhas(ip)) ' deg'])
    ylim([-0.1 1.2])
end
iexp = expt_ind(ind_h(i));
suptitle([expt(iexp).mouse ' ' expt(iexp).date ' - cell #' num2str(ind_h(i))])
print(fullfile(summaryDir_F6, 'Figure6_ExampleCell_posMod.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F6, 'Figure6_ExampleCell_posMod.fig'))

figure;
subplot(3,3,1)
shadedErrorBar(tt', resp_tc_all(:,ind_l(i),1,nStimCon,1,1), resp_tc_all(:,ind_l(i),1,nStimCon,1,2),'lineprops',{['b-'], 'markerfacecolor', 'b'});
hold on
shadedErrorBar(tt, resp_tc_all(:,ind_l(i),nMaskCon,1,1,1), resp_tc_all(:,ind_l(i),nMaskCon,1,1,2),'lineprops',{['r-'], 'markerfacecolor','r'});
xlabel('Time (ms)')
ylabel('dF/F')
title('Stim/Mask')
ylim([-0.1 1.2])
for ip = 1:nMaskPhas
    subplot(3,3,ip+1)
    shadedErrorBar(tt, resp_tc_all(:,ind_l(i),nMaskCon,nStimCon,ip,1), resp_tc_all(:,ind_l(i),nMaskCon,nStimCon,ip,2));
    xlabel('Time (ms)')
    ylabel('dF/F')
    title(['Phase ' num2str(maskPhas(ip)) ' deg'])
    ylim([-0.1 1.2])
end
iexp = expt_ind(ind_l(i));
suptitle([expt(iexp).mouse ' ' expt(iexp).date ' - cell #' num2str(ind_l(i))])
print(fullfile(summaryDir_F6, 'Figure6_ExampleCell_noMod.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F6, 'Figure6_ExampleCell_noMod.fig'))

figure;
subplot(1,2,1)
errorbar(maskPhas, squeeze(SI_avg_all(ind_h(i),nMaskCon,nStimCon,:,1)),squeeze(SI_avg_all(ind_h(i),nMaskCon,nStimCon,:,2)),'ok')
hold on
plot(1:360,yfit_allcells(ind_h(i),:),'-r')
xlim([-22.5 360])
ylim([-1.2 1.2])
xlabel('Phase (deg)')
ylabel('Selectivity index')
title(['Cell #' num2str(ind_h(i)) '; Amp: ' num2str(chop(amp_all(ind_h(i)),2)) '; Base: ' num2str(chop(b_all(ind_h(i)),2))])
subplot(1,2,2)
errorbar(maskPhas, squeeze(SI_avg_all(ind_l(i),nMaskCon,nStimCon,:,1)),squeeze(SI_avg_all(ind_l(i),nMaskCon,nStimCon,:,2)),'ok')
hold on
plot(1:360,yfit_allcells(ind_l(i),:),'-r')
xlim([-22.5 360])
ylim([-1.2 1.2])
xlabel('Phase (deg)')
ylabel('Selectivity index')
title(['Cell #' num2str(ind_l(i)) '; Amp: ' num2str(chop(amp_all(ind_l(i)),2)) '; Base: ' num2str(chop(b_all(ind_l(i)),2))])
suptitle([expt(iexp).mouse ' ' expt(iexp).date])
print(fullfile(summaryDir_F6, 'Figure6_ExampleCells_Tuning&Fit.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F6, 'Figure6_ExampleCells_Tuning&Fit.fig'))

figure;
subplot(3,2,1)
scatter(OSI_all(resp_ind_all_phase), amp_all(resp_ind_all_phase))
hold on
[n edges bin] = histcounts(OSI_all,[0:0.1:1]);
sin_avg = zeros(length(n),2);
OSI_avg = zeros(length(n),2);
ind_use = intersect(find(~isnan(OSI_all)),resp_ind_all_phase);
for i = 1:length(n)
    ind = intersect(resp_ind_all_phase,find(bin == i));
    sin_avg(i,1) = nanmean(amp_all(ind),1);
    sin_avg(i,2) = nanstd(amp_all(ind),[],1)./sqrt(length(ind));
    OSI_avg(i,1) = nanmean(OSI_all(ind),1);
    OSI_avg(i,2) = nanstd(OSI_all(ind),[],1)./sqrt(length(ind));
end
errorbar(OSI_avg(:,1),sin_avg(:,1),sin_avg(:,2),sin_avg(:,2),OSI_avg(:,2),OSI_avg(:,2),'-o')
r = [ones(size(OSI_all(ind_use))) OSI_all(ind_use)]\amp_all(ind_use);
text(0.05,0.7,['r = ' num2str(chop(r(2),2))]);
xlabel('OSI')
ylabel('Sine amplitude')
ylim([-1 1])
title([num2str(sum(~isnan(OSI_all(resp_ind_all_phase)))) ' cells'])
subplot(3,2,2)
ind_use = intersect(find(~isnan(OSI_all)),intersect(resp_ind_all_phase,unique([find(max_dir_all==1); find(max_dir_all==5); find(max_dir_all==9); find(max_dir_all==13)])));
scatter(OSI_all(ind_use), amp_all(ind_use))
hold on
sin_avg = zeros(length(n),2);
OSI_avg = zeros(length(n),2);
for i = 1:length(n)
    ind = intersect(ind_use,find(bin == i));
    sin_avg(i,1) = nanmean(amp_all(ind),1);
    sin_avg(i,2) = nanstd(amp_all(ind),[],1)./sqrt(length(ind));
    OSI_avg(i,1) = nanmean(OSI_all(ind),1);
    OSI_avg(i,2) = nanstd(OSI_all(ind),[],1)./sqrt(length(ind));
end
errorbar(OSI_avg(:,1),sin_avg(:,1),sin_avg(:,2),sin_avg(:,2),OSI_avg(:,2),OSI_avg(:,2),'-o')
r = [ones(size(OSI_all(ind_use))) OSI_all(ind_use)]\amp_all(ind_use);
text(0.05,0.7,['r = ' num2str(chop(r(2),2))]);
xlabel('OSI')
ylabel('Sine amplitude')
ylim([-1 1])
title([num2str(sum(~isnan(OSI_all(ind_use)))) ' cells- test/mask ori is pref'])
subplot(3,2,3)
scatter(OSI_all(resp_ind_all_phase), phase_MI_all(resp_ind_all_phase))
hold on
ind_use = intersect(find(~isnan(OSI_all)),resp_ind_all_phase);
[n edges bin] = histcounts(OSI_all,[0:0.1:1]);
SI_avg = zeros(length(n),2);
OSI_avg = zeros(length(n),2);
for i = 1:length(n)
    ind = intersect(resp_ind_all_phase,find(bin == i));
    SI_avg(i,1) = nanmean(phase_MI_all(ind),1);
    SI_avg(i,2) = nanstd(phase_MI_all(ind),[],1)./sqrt(length(ind));
    OSI_avg(i,1) = nanmean(OSI_all(ind),1);
    OSI_avg(i,2) = nanstd(OSI_all(ind),[],1)./sqrt(length(ind));
end
errorbar(OSI_avg(:,1),SI_avg(:,1),SI_avg(:,2),SI_avg(:,2),OSI_avg(:,2),OSI_avg(:,2),'-o')
r = [ones(size(OSI_all(ind_use))) OSI_all(ind_use)]\phase_MI_all(ind_use);
text(0.05,0.7,['r = ' num2str(chop(r(2),2))]);
xlabel('OSI')
ylabel('Masking Index')
ylim([-1 1])
title([num2str(sum(~isnan(OSI_all(resp_ind_all_phase)))) ' cells'])
subplot(3,2,4)
ind_use = intersect(find(~isnan(OSI_all)),intersect(resp_ind_all_phase,unique([find(max_dir_all==1); find(max_dir_all==5); find(max_dir_all==9); find(max_dir_all==13)])));
scatter(OSI_all(ind_use), phase_MI_all(ind_use))
hold on
SI_avg = zeros(length(n),2);
OSI_avg = zeros(length(n),2);
for i = 1:length(n)
    ind = intersect(ind_use,find(bin == i));
    SI_avg(i,1) = nanmean(phase_MI_all(ind),1);
    SI_avg(i,2) = nanstd(phase_MI_all(ind),[],1)./sqrt(length(ind));
    OSI_avg(i,1) = nanmean(OSI_all(ind),1);
    OSI_avg(i,2) = nanstd(OSI_all(ind),[],1)./sqrt(length(ind));
end
errorbar(OSI_avg(:,1),SI_avg(:,1),SI_avg(:,2),SI_avg(:,2),OSI_avg(:,2),OSI_avg(:,2),'-o')
r = [ones(size(OSI_all(ind_use))) OSI_all(ind_use)]\phase_MI_all(ind_use);
text(0.05,0.7,['r = ' num2str(chop(r(2),2))]);
xlabel('OSI')
ylabel('Masking Index')
ylim([-1 1])
title([num2str(sum(~isnan(OSI_all(ind_use)))) ' cells- test/mask ori is pref'])
subplot(3,2,5)
scatter(OSI_all(resp_ind_all_phase), b_all(resp_ind_all_phase))
hold on
ind_use = intersect(find(~isnan(OSI_all)),resp_ind_all_phase);
[n edges bin] = histcounts(OSI_all,[0:0.1:1]);
b_avg = zeros(length(n),2);
OSI_avg = zeros(length(n),2);
for i = 1:length(n)
    ind = intersect(resp_ind_all_phase,find(bin == i));
    b_avg(i,1) = nanmean(b_all(ind),1);
    b_avg(i,2) = nanstd(b_all(ind),[],1)./sqrt(length(ind));
    OSI_avg(i,1) = nanmean(OSI_all(ind),1);
    OSI_avg(i,2) = nanstd(OSI_all(ind),[],1)./sqrt(length(ind));
end
errorbar(OSI_avg(:,1),b_avg(:,1),b_avg(:,2),b_avg(:,2),OSI_avg(:,2),OSI_avg(:,2),'-o')
r = [ones(size(OSI_all(ind_use))) OSI_all(ind_use)]\b_all(ind_use);
text(0.05,0.7,['r = ' num2str(chop(r(2),2))]);
xlabel('OSI')
ylabel('Sine baseline')
ylim([-1 1])
title([num2str(sum(~isnan(OSI_all(resp_ind_all_phase)))) ' cells'])
subplot(3,2,6)
ind_use = intersect(find(~isnan(OSI_all)),intersect(resp_ind_all_phase,unique([find(max_dir_all==1); find(max_dir_all==5); find(max_dir_all==9); find(max_dir_all==13)])));
scatter(OSI_all(ind_use), b_all(ind_use))
hold on
b_avg = zeros(length(n),2);
OSI_avg = zeros(length(n),2);
for i = 1:length(n)
    ind = intersect(ind_use,find(bin == i));
    b_avg(i,1) = nanmean(b_all(ind),1);
    b_avg(i,2) = nanstd(b_all(ind),[],1)./sqrt(length(ind));
    OSI_avg(i,1) = nanmean(OSI_all(ind),1);
    OSI_avg(i,2) = nanstd(OSI_all(ind),[],1)./sqrt(length(ind));
end
errorbar(OSI_avg(:,1),b_avg(:,1),b_avg(:,2),b_avg(:,2),OSI_avg(:,2),OSI_avg(:,2),'-o')
r = [ones(size(OSI_all(ind_use))) OSI_all(ind_use)]\b_all(ind_use);
text(0.05,0.7,['r = ' num2str(chop(r(2),2))]);
xlabel('OSI')
ylabel('Sine baseline')
ylim([-1 1])
title([num2str(sum(~isnan(phase_SI_all(ind_use)))) ' cells- test/mask ori is pref'])
print(fullfile(summaryDir_F7, 'Figure7_SineAmplitudeOSI_Summary.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F7, 'Figure7_SineAmplitudeOSI_Summary.fig'))


figure;
subplot(2,2,1)
scatter(b_all(resp_ind_all_phase), amp_all(resp_ind_all_phase))
hold on
[n edges bin] = histcounts(b_all,[-1:0.2:1]);
sin_avg = zeros(length(n),2);
b_avg = zeros(length(n),2);
ind_use = intersect(find(~isnan(b_all)),resp_ind_all_phase);
for i = 1:length(n)
    ind = intersect(resp_ind_all_phase,find(bin == i));
    sin_avg(i,1) = nanmean(amp_all(ind),1);
    sin_avg(i,2) = nanstd(amp_all(ind),[],1)./sqrt(length(ind));
    b_avg(i,1) = nanmean(b_all(ind),1);
    b_avg(i,2) = nanstd(b_all(ind),[],1)./sqrt(length(ind));
end
errorbar(b_avg(:,1),sin_avg(:,1),sin_avg(:,2),sin_avg(:,2),b_avg(:,2),b_avg(:,2),'-o')
r = [ones(size(b_all(ind_use))) b_all(ind_use)]\amp_all(ind_use);
text(0.05,0.7,['r = ' num2str(chop(r(2),2))]);
xlabel('Sine baseline')
ylabel('Sine amplitude')
xlim([-1 0.5])
ylim([-1 1])
title([num2str(sum(~isnan(b_all(resp_ind_all_phase)))) ' cells'])
subplot(2,2,2)
ind_use = intersect(find(~isnan(b_all)),intersect(resp_ind_all_phase,unique([find(max_dir_all==1); find(max_dir_all==5); find(max_dir_all==9); find(max_dir_all==13)])));
scatter(b_all(ind_use), amp_all(ind_use))
hold on
sin_avg = zeros(length(n),2);
b_avg = zeros(length(n),2);
for i = 1:length(n)
    ind = intersect(ind_use,find(bin == i));
    sin_avg(i,1) = nanmean(amp_all(ind),1);
    sin_avg(i,2) = nanstd(amp_all(ind),[],1)./sqrt(length(ind));
    b_avg(i,1) = nanmean(b_all(ind),1);
    b_avg(i,2) = nanstd(b_all(ind),[],1)./sqrt(length(ind));
end
errorbar(b_avg(:,1),sin_avg(:,1),sin_avg(:,2),sin_avg(:,2),b_avg(:,2),b_avg(:,2),'-o')
r = [ones(size(b_all(ind_use))) b_all(ind_use)]\amp_all(ind_use);
text(0.05,0.7,['r = ' num2str(chop(r(2),2))]);
xlabel('Sine baseline')
ylabel('Sine amplitude')
ylim([-1 1])
xlim([-1 0.5])
title([num2str(sum(~isnan(b_all(ind_use)))) ' cells- test/mask ori is pref'])
print(fullfile(summaryDir_F7, 'Figure7_SineAmplitudeBaseline_Summary.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F7, 'Figure7_SineAmplitudeBaseline_Summary.fig'))

figure;
scatter(phase_SI_all(resp_ind_all_phase), amp_all(resp_ind_all_phase))
hold on
scatter(phase_SI_all(intersect(find(nP_ind==8),resp_ind_all_phase)), amp_all(intersect(find(nP_ind==8),resp_ind_all_phase)))

figure;
subplot(3,2,1)
scatter(phase_SI_all(resp_ind_all_phase), amp_all(resp_ind_all_phase))
hold on
[n edges bin] = histcounts(phase_SI_all,[0:0.1:1]);
sin_avg = zeros(length(n),2);
stimSI_avg = zeros(length(n),2);
ind_use = resp_ind_all_phase;
for i = 1:length(n)
    ind = intersect(resp_ind_all_phase,find(bin == i));
    sin_avg(i,1) = nanmean(amp_all(ind),1);
    sin_avg(i,2) = nanstd(amp_all(ind),[],1)./sqrt(length(ind));
    stimSI_avg(i,1) = nanmean(phase_SI_all(ind),1);
    stimSI_avg(i,2) = nanstd(phase_SI_all(ind),[],1)./sqrt(length(ind));
end
errorbar(stimSI_avg(:,1),sin_avg(:,1),sin_avg(:,2),sin_avg(:,2),stimSI_avg(:,2),stimSI_avg(:,2),'-o')
r = [ones(size(phase_SI_all(ind_use))) phase_SI_all(ind_use)]\amp_all(ind_use);
text(0.05,0.7,['r = ' num2str(chop(r(2),2))]);
xlabel('SI- [abs(M-T)/(M+T)]')
ylabel('Sine amplitude')
ylim([-1 1])
title([num2str(sum(~isnan(phase_SI_all(resp_ind_all_phase)))) ' cells'])
subplot(3,2,2)
ind_use = intersect(resp_ind_all_phase,unique([find(max_dir_all==1); find(max_dir_all==5); find(max_dir_all==9); find(max_dir_all==13)]));
scatter(phase_SI_all(ind_use), amp_all(ind_use))
hold on
sin_avg = zeros(length(n),2);
stimSI_avg = zeros(length(n),2);
for i = 1:length(n)
    ind = intersect(ind_use,find(bin == i));
    sin_avg(i,1) = nanmean(amp_all(ind),1);
    sin_avg(i,2) = nanstd(amp_all(ind),[],1)./sqrt(length(ind));
    stimSI_avg(i,1) = nanmean(phase_SI_all(ind),1);
    stimSI_avg(i,2) = nanstd(phase_SI_all(ind),[],1)./sqrt(length(ind));
end
errorbar(stimSI_avg(:,1),sin_avg(:,1),sin_avg(:,2),sin_avg(:,2),stimSI_avg(:,2),stimSI_avg(:,2),'-o')
r = [ones(size(phase_SI_all(ind_use))) phase_SI_all(ind_use)]\amp_all(ind_use);
text(0.05,0.7,['r = ' num2str(chop(r(2),2))]);
xlabel('SI')
ylabel('Sine amplitude')
ylim([-1 1])
title([num2str(sum(~isnan(phase_SI_all(ind_use)))) ' cells- test/mask ori is pref'])
subplot(3,2,3)
scatter(phase_SI_all(resp_ind_all_phase), phase_MI_all(resp_ind_all_phase))
hold on
ind_use = resp_ind_all_phase;
[n edges bin] = histcounts(phase_SI_all,[0:0.1:1]);
SI_avg = zeros(length(n),2);
stimSI_avg = zeros(length(n),2);
for i = 1:length(n)
    ind = intersect(resp_ind_all_phase,find(bin == i));
    SI_avg(i,1) = nanmean(phase_MI_all(ind),1);
    SI_avg(i,2) = nanstd(phase_MI_all(ind),[],1)./sqrt(length(ind));
    stimSI_avg(i,1) = nanmean(phase_SI_all(ind),1);
    stimSI_avg(i,2) = nanstd(phase_SI_all(ind),[],1)./sqrt(length(ind));
end
errorbar(stimSI_avg(:,1),SI_avg(:,1),SI_avg(:,2),SI_avg(:,2),stimSI_avg(:,2),stimSI_avg(:,2),'-o')
r = [ones(size(phase_SI_all(ind_use))) phase_SI_all(ind_use)]\phase_MI_all(ind_use);
text(0.05,0.7,['r = ' num2str(chop(r(2),2))]);
xlabel('SI- [abs(M-T)/(M+T)]')
ylabel('Masking Index')
ylim([-1 1])
title([num2str(sum(~isnan(phase_SI_all(resp_ind_all_phase)))) ' cells'])
subplot(3,2,4)
ind_use = intersect(resp_ind_all_phase,unique([find(max_dir_all==1); find(max_dir_all==5); find(max_dir_all==9); find(max_dir_all==13)]));
scatter(phase_SI_all(ind_use), phase_MI_all(ind_use))
hold on
SI_avg = zeros(length(n),2);
stimSI_avg = zeros(length(n),2);
for i = 1:length(n)
    ind = intersect(ind_use,find(bin == i));
    SI_avg(i,1) = nanmean(phase_MI_all(ind),1);
    SI_avg(i,2) = nanstd(phase_MI_all(ind),[],1)./sqrt(length(ind));
    stimSI_avg(i,1) = nanmean(phase_SI_all(ind),1);
    stimSI_avg(i,2) = nanstd(phase_SI_all(ind),[],1)./sqrt(length(ind));
end
errorbar(stimSI_avg(:,1),SI_avg(:,1),SI_avg(:,2),SI_avg(:,2),stimSI_avg(:,2),stimSI_avg(:,2),'-o')
r = [ones(size(phase_SI_all(ind_use))) phase_SI_all(ind_use)]\phase_MI_all(ind_use);
text(0.05,0.7,['r = ' num2str(chop(r(2),2))]);
xlabel('SI')
ylabel('Masking Index')
ylim([-1 1])
title([num2str(sum(~isnan(phase_SI_all(ind_use)))) ' cells- test/mask ori is pref'])
subplot(3,2,5)
scatter(phase_SI_all(resp_ind_all_phase), b_all(resp_ind_all_phase))
hold on
ind_use = resp_ind_all_phase;
[n edges bin] = histcounts(phase_SI_all,[0:0.1:1]);
b_avg = zeros(length(n),2);
stimSI_avg = zeros(length(n),2);
for i = 1:length(n)
    ind = intersect(resp_ind_all_phase,find(bin == i));
    b_avg(i,1) = nanmean(b_all(ind),1);
    b_avg(i,2) = nanstd(b_all(ind),[],1)./sqrt(length(ind));
    stimSI_avg(i,1) = nanmean(phase_SI_all(ind),1);
    stimSI_avg(i,2) = nanstd(phase_SI_all(ind),[],1)./sqrt(length(ind));
end
errorbar(stimSI_avg(:,1),b_avg(:,1),b_avg(:,2),b_avg(:,2),stimSI_avg(:,2),stimSI_avg(:,2),'-o')
r = [ones(size(phase_SI_all(ind_use))) phase_SI_all(ind_use)]\b_all(ind_use);
text(0.05,0.7,['r = ' num2str(chop(r(2),2))]);
xlabel('SI- [abs(M-T)/(M+T)]')
ylabel('Sine baseline')
ylim([-1 1])
title([num2str(sum(~isnan(phase_SI_all(resp_ind_all_phase)))) ' cells'])
subplot(3,2,6)
ind_use = intersect(resp_ind_all_phase,unique([find(max_dir_all==1); find(max_dir_all==5); find(max_dir_all==9); find(max_dir_all==13)]));
scatter(phase_SI_all(ind_use), b_all(ind_use))
hold on
b_avg = zeros(length(n),2);
stimSI_avg = zeros(length(n),2);
for i = 1:length(n)
    ind = intersect(ind_use,find(bin == i));
    b_avg(i,1) = nanmean(b_all(ind),1);
    b_avg(i,2) = nanstd(b_all(ind),[],1)./sqrt(length(ind));
    stimSI_avg(i,1) = nanmean(phase_SI_all(ind),1);
    stimSI_avg(i,2) = nanstd(phase_SI_all(ind),[],1)./sqrt(length(ind));
end
errorbar(stimSI_avg(:,1),b_avg(:,1),b_avg(:,2),b_avg(:,2),stimSI_avg(:,2),stimSI_avg(:,2),'-o')
r = [ones(size(phase_SI_all(ind_use))) phase_SI_all(ind_use)]\b_all(ind_use);
text(0.05,0.7,['r = ' num2str(chop(r(2),2))]);
xlabel('SI')
ylabel('Sine baseline')
ylim([-1 1])
title([num2str(sum(~isnan(phase_SI_all(ind_use)))) ' cells- test/mask ori is pref'])
print(fullfile(summaryDir_F7, 'Figure7_SineAmplitudeSI_Summary.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F7, 'Figure7_SineAmplitudeSI_Summary.fig'))


figure;
subplot(2,2,1)
scatter(f2overf1_all(resp_ind_all_phase),amp_all(resp_ind_all_phase))
hold on
[n edges bin] = histcounts(f2overf1_all,[0:0.2:2]);
f2f1_avg = zeros(length(n),2);
amp_avg = zeros(length(n),2);
ind_use = intersect(find(~isnan(f2overf1_all)),resp_ind_all_phase);
for i = 1:length(n)
    ind = intersect(resp_ind_all_phase,find(bin == i));
    f2f1_avg(i,1) = nanmean(f2overf1_all(ind),1);
    f2f1_avg(i,2) = nanstd(f2overf1_all(ind),[],1)./sqrt(length(ind));
    amp_avg(i,1) = nanmean(amp_all(ind),1);
    amp_avg(i,2) = nanstd(amp_all(ind),[],1)./sqrt(length(ind));
end
errorbar(f2f1_avg(:,1),amp_avg(:,1),amp_avg(:,2),amp_avg(:,2),f2f1_avg(:,2),f2f1_avg(:,2),'-o')
r = [ones(size(amp_all(ind_use))) amp_all(ind_use)]\f2overf1_all(ind_use);
text(3,0.7,['r = ' num2str(chop(r(2),2))]);
xlabel('F2/F1')
ylabel('Sine amplitude')
title([num2str(sum(~isnan(f2overf1_all(resp_ind_all_phase)))) ' cells'])
subplot(2,2,2)
scatter(f2overf1_all(resp_ind_all_phase),phase_MI_all(resp_ind_all_phase))
hold on
[n edges bin] = histcounts(f2overf1_all,[0:0.2:2]);
f2f1_avg = zeros(length(n),2);
MI_avg = zeros(length(n),2);
ind_use = intersect(find(~isnan(f2overf1_all)),resp_ind_all_phase);
for i = 1:length(n)
    ind = intersect(resp_ind_all_phase,find(bin == i));
    f2f1_avg(i,1) = nanmean(f2overf1_all(ind),1);
    f2f1_avg(i,2) = nanstd(f2overf1_all(ind),[],1)./sqrt(length(ind));
    MI_avg(i,1) = nanmean(phase_MI_all(ind),1);
    MI_avg(i,2) = nanstd(phase_MI_all(ind),[],1)./sqrt(length(ind));
end
errorbar(f2f1_avg(:,1),MI_avg(:,1),MI_avg(:,2),MI_avg(:,2),f2f1_avg(:,2),f2f1_avg(:,2),'-o')
r = [ones(size(phase_MI_all(ind_use))) phase_MI_all(ind_use)]\f2overf1_all(ind_use);
text(3,0.7,['r = ' num2str(chop(r(2),2))]);
xlabel('F2/F1')
ylabel('Masking Index')
title([num2str(sum(~isnan(f2overf1_all(resp_ind_all_phase)))) ' cells'])
subplot(2,2,3)
scatter(f2overf1_all(resp_ind_all_phase),phase_MI_all(resp_ind_all_phase))
hold on
[n edges bin] = histcounts(f2overf1_all,[0:0.2:2]);
f2f1_avg = zeros(length(n),2);
b_avg = zeros(length(n),2);
ind_use = intersect(find(~isnan(f2overf1_all)),resp_ind_all_phase);
for i = 1:length(n)
    ind = intersect(resp_ind_all_phase,find(bin == i));
    f2f1_avg(i,1) = nanmean(f2overf1_all(ind),1);
    f2f1_avg(i,2) = nanstd(f2overf1_all(ind),[],1)./sqrt(length(ind));
    b_avg(i,1) = nanmean(b_all(ind),1);
    b_avg(i,2) = nanstd(b_all(ind),[],1)./sqrt(length(ind));
end
errorbar(f2f1_avg(:,1),b_avg(:,1),b_avg(:,2),b_avg(:,2),f2f1_avg(:,2),f2f1_avg(:,2),'-o')
r = [ones(size(b_all(ind_use))) b_all(ind_use)]\f2overf1_all(ind_use);
text(3,0.7,['r = ' num2str(chop(r(2),2))]);
xlabel('F2/F1')
ylabel('Sine baseline')
title([num2str(sum(~isnan(f2overf1_all(resp_ind_all_phase)))) ' cells'])
print(fullfile(summaryDir_F7, 'Figure7_SineAmplitudeF2F1_Summary.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F7, 'Figure7_SineAmplitudeF2F1_Summary.fig'))



figure;
subplot(2,3,1)
scatter(phase_SI_all(resp_ind_all_phase), phase_MI_all(resp_ind_all_phase))
hold on
[n edges bin] = histcounts(phase_SI_all,[0:0.1:1]);
plaidSI_avg = zeros(length(n),2);
stimSI_avg = zeros(length(n),2);
for i = 1:length(n)
    ind = intersect(resp_ind_all_phase,find(bin == i));
    plaidSI_avg(i,1) = nanmean(phase_MI_all(ind),1);
    plaidSI_avg(i,2) = nanstd(phase_MI_all(ind),[],1)./sqrt(length(ind));
    stimSI_avg(i,1) = nanmean(phase_SI_all(ind),1);
    stimSI_avg(i,2) = nanstd(phase_SI_all(ind),[],1)./sqrt(length(ind));
end
errorbar(stimSI_avg(:,1),plaidSI_avg(:,1),plaidSI_avg(:,2),plaidSI_avg(:,2),stimSI_avg(:,2),stimSI_avg(:,2),'-o')
xlabel('Stim Preference- [abs(M-T)/(M+T)]')
ylabel('Selectivity Index')
title([num2str(sum(~isnan(phase_SI_all(resp_ind_all_phase)))) ' cells'])
subplot(2,3,2)
scatter(OSI_all(resp_ind_all_phase), phase_MI_all(resp_ind_all_phase))
hold on
[n, edges bin] = histcounts(OSI_all,[0:0.1:1]);
plaidSI_avg = zeros(length(n),2);
OSI_avg = zeros(length(n),2);
for i = 1:length(n)
    ind = intersect(resp_ind_all_phase,find(bin == i));
    plaidSI_avg(i,1) = nanmean(phase_MI_all(ind),1);
    plaidSI_avg(i,2) = nanstd(phase_MI_all(ind),[],1)./sqrt(length(ind));
    OSI_avg(i,1) = nanmean(OSI_all(ind),1);
    OSI_avg(i,2) = nanstd(OSI_all(ind),[],1)./sqrt(length(ind));
end
errorbar(OSI_avg(:,1),plaidSI_avg(:,1),plaidSI_avg(:,2),plaidSI_avg(:,2),OSI_avg(:,2),OSI_avg(:,2),'-o')
xlabel('OSI')
ylabel('Selectivity Index')
title([num2str(sum(~isnan(OSI_all(resp_ind_all_phase)))) ' cells'])
subplot(2,3,3)
ind_use = intersect(resp_ind_all_phase,unique([find(peak_dir_all<11.25); find(peak_dir_all>348.75);...
    intersect(find(peak_dir_all>78.75),find(peak_dir_all<101.25));...
    intersect(find(DSI_all<0.5), intersect(find(peak_dir_all>168.75), find(peak_dir_all<191.25)));...
    intersect(find(DSI_all<0.5), intersect(find(peak_dir_all>258.75), find(peak_dir_all<281.25)))]));
scatter(OSI_all(ind_use), phase_MI_all(ind_use))
hold on
plaidSI_avg = zeros(length(n),2);
OSI_avg = zeros(length(n),2);
for i = 1:length(n)
    ind = intersect(ind_use,find(bin == i));
    plaidSI_avg(i,1) = nanmean(phase_MI_all(ind),1);
    plaidSI_avg(i,2) = nanstd(phase_MI_all(ind),[],1)./sqrt(length(ind));
    OSI_avg(i,1) = nanmean(OSI_all(ind),1);
    OSI_avg(i,2) = nanstd(OSI_all(ind),[],1)./sqrt(length(ind));
end
errorbar(OSI_avg(:,1),plaidSI_avg(:,1),plaidSI_avg(:,2),plaidSI_avg(:,2),OSI_avg(:,2),OSI_avg(:,2),'-o')
xlabel('OSI')
ylabel('Selectivity Index')
title([num2str(sum(~isnan(OSI_all(ind_use)))) ' cells- test/mask is pref'])
subplot(2,3,4)
ind1 = intersect(resp_ind_all_phase,find(phase_SI_all>0.5));
cdfplot(phase_MI_all(ind1));
hold on
ind2 = intersect(resp_ind_all_phase,find(phase_SI_all<0.5));
cdfplot(phase_MI_all(ind2));
legend({['stimSI>0.5 - n=' num2str(length(ind1))], ['stimSI<0.5 - n=' num2str(length(ind2))]},'location','southeast')
xlabel('Selectivity Index')
title('')
subplot(2,3,5)
ind1 = intersect(resp_ind_all_phase,find(OSI_all>0.5));
cdfplot(phase_MI_all(ind1));
hold on
ind2 = intersect(resp_ind_all_phase,find(OSI_all<0.5));
cdfplot(phase_MI_all(ind2));
legend({['OSI>0.5 - n=' num2str(length(ind1))], ['OSI<0.5 - n=' num2str(length(ind2))]},'location','southeast')
xlabel('Selectivity Index')
title('')
subplot(2,3,6)
ind1 = intersect(ind_use,find(OSI_all>0.5));
cdfplot(phase_MI_all(ind1));
hold on
ind2 = intersect(ind_use,find(OSI_all<0.5));
cdfplot(phase_MI_all(ind2));
legend({['OSI>0.5 - n=' num2str(length(ind1))], ['OSI<0.5 - n=' num2str(length(ind2))]},'location','southeast')
xlabel('Selectivity Index')
title('')
suptitle('Selectivity Index (matched phase)')
print(fullfile(summaryDir_F7, 'Figure7_plaidSIOSI_Summary.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F7, 'Figure7_plaidSIOSI_Summary.fig'))


figure;
subplot(3,2,1)
scatter(mask_resp_all(resp_ind_all_phase)+test_resp_all(resp_ind_all_phase), plaid_resp_all(resp_ind_all_phase))
xlim([0 3])
ylim([0 3])
refline(1)
xlabel('Mask + Test dF/F')
ylabel('Plaid dF/F')
title({'Test or Mask responsive',['n = ' num2str(length(resp_ind_all_phase))]})
subplot(3,2,2)
scatter(mask_resp_all(ind_use)+test_resp_all(ind_use), plaid_resp_all(ind_use))
xlim([0 3])
ylim([0 3])
refline(1)
xlabel('Mask + Test dF/F')
ylabel('Plaid dF/F')
title({'Test or Mask preferring',['n = ' num2str(length(ind_use))]})
subplot(3,2,3)
histogram(phase_MI_all(resp_ind_all_phase))
xlabel('Selectivity index')
ylabel('Number of cells')
xlim([-1 1])
subplot(3,2,4)
histogram(phase_MI_all(ind_use))
xlabel('Selectivity index')
ylabel('Number of cells')
xlim([-1 1])
subplot(3,2,5)
cdfplot(phase_MI_all(resp_ind_all_phase))
xlabel('Selectivity index')
ylabel('Fraction of cells')
xlim([-1 1])
subplot(3,2,6)
cdfplot(phase_MI_all(ind_use))
xlabel('Selectivity index')
ylabel('Fraction of cells')
xlim([-1 1])
print(fullfile(summaryDir_F7, 'Figure7_SI_allVpref_Summary.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F7, 'Figure7_allVpref_Summary.fig'))

figure;
subplot(2,3,1)
scatter(phase_SI_all(resp_ind_all_phase), b_all(resp_ind_all_phase))
hold on
[n edges bin] = histcounts(phase_SI_all,[0:0.1:1]);
b_all_avg = zeros(length(n),2);
stimSI_avg = zeros(length(n),2);
for i = 1:length(n)
    ind = intersect(resp_ind_all_phase,find(bin == i));
    b_all_avg(i,1) = nanmean(b_all(ind),1);
    b_all_avg(i,2) = nanstd(b_all(ind),[],1)./sqrt(length(ind));
    stimSI_avg(i,1) = nanmean(phase_SI_all(ind),1);
    stimSI_avg(i,2) = nanstd(phase_SI_all(ind),[],1)./sqrt(length(ind));
end
errorbar(stimSI_avg(:,1),b_all_avg(:,1),b_all_avg(:,2),b_all_avg(:,2),stimSI_avg(:,2),stimSI_avg(:,2),'-o')
xlabel('Stim Preference- [abs(M-T)/(M+T)]')
ylabel('Selectivity Index')
title([num2str(sum(~isnan(phase_SI_all(resp_ind_all_phase)))) ' cells'])
subplot(2,3,2)
scatter(OSI_all(resp_ind_all_phase), b_all(resp_ind_all_phase))
hold on
[n edges bin] = histcounts(OSI_all,[0:0.1:1]);
b_all_avg = zeros(length(n),2);
OSI_avg = zeros(length(n),2);
for i = 1:length(n)
    ind = intersect(resp_ind_all_phase,find(bin == i));
    b_all_avg(i,1) = nanmean(b_all(ind),1);
    b_all_avg(i,2) = nanstd(b_all(ind),[],1)./sqrt(length(ind));
    OSI_avg(i,1) = nanmean(OSI_all(ind),1);
    OSI_avg(i,2) = nanstd(OSI_all(ind),[],1)./sqrt(length(ind));
end
errorbar(OSI_avg(:,1),b_all_avg(:,1),b_all_avg(:,2),b_all_avg(:,2),OSI_avg(:,2),OSI_avg(:,2),'-o')
xlabel('OSI')
ylabel('Selectivity Index')
title([num2str(sum(~isnan(OSI_all(resp_ind_all_phase)))) ' cells'])
subplot(2,3,3)
ind_use = intersect(resp_ind_all_phase,unique([find(peak_dir_all<22.5); find(peak_dir_all>337.5);...
    intersect(find(DSI_all<0.5), [find(peak_dir_all>157.5); find(peak_dir_all<202.5)])]));
scatter(OSI_all(ind_use), b_all(ind_use))
hold on
b_all_avg = zeros(length(n),2);
OSI_avg = zeros(length(n),2);
for i = 1:length(n)
    ind = intersect(ind_use,find(bin == i));
    b_all_avg(i,1) = nanmean(b_all(ind),1);
    b_all_avg(i,2) = nanstd(b_all(ind),[],1)./sqrt(length(ind));
    OSI_avg(i,1) = nanmean(OSI_all(ind),1);
    OSI_avg(i,2) = nanstd(OSI_all(ind),[],1)./sqrt(length(ind));
end
errorbar(OSI_avg(:,1),b_all_avg(:,1),b_all_avg(:,2),b_all_avg(:,2),OSI_avg(:,2),OSI_avg(:,2),'-o')
xlabel('OSI')
ylabel('Selectivity Index')
title([num2str(sum(~isnan(OSI_all(ind_use)))) ' cells- test is pref'])
subplot(2,3,4)
ind1 = intersect(resp_ind_all_phase,find(phase_SI_all>0.5));
cdfplot(b_all(ind1));
hold on
ind2 = intersect(resp_ind_all_phase,find(phase_SI_all<0.5));
cdfplot(b_all(ind2));
legend({['stimSI>0.5 - n=' num2str(length(ind1))], ['stimSI<0.5 - n=' num2str(length(ind2))]},'location','southeast')
xlabel('Selectivity Index')
title('')
subplot(2,3,5)
ind1 = intersect(resp_ind_all_phase,find(OSI_all>0.5));
cdfplot(b_all(ind1));
hold on
ind2 = intersect(resp_ind_all_phase,find(OSI_all<0.5));
cdfplot(b_all(ind2));
legend({['OSI>0.5 - n=' num2str(length(ind1))], ['OSI<0.5 - n=' num2str(length(ind2))]},'location','southeast')
xlabel('Selectivity Index')
title('')
subplot(2,3,6)
ind1 = intersect(ind_use,find(OSI_all>0.5));
cdfplot(b_all(ind1));
hold on
ind2 = intersect(ind_use,find(OSI_all<0.5));
cdfplot(b_all(ind2));
legend({['OSI>0.5 - n=' num2str(length(ind1))], ['OSI<0.5 - n=' num2str(length(ind2))]},'location','southeast')
xlabel('Selectivity Index')
title('')
suptitle('Selectivity Index (avg all phase)')
print(fullfile(summaryDir_F7, 'Figure7_allPhaseSIOSI_Summary.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F7, 'Figure7_allPhaseSIOSI_Summary.fig'))

pr_run = [{expt2.prFolder}];
pr_mice = [{expt2.mouse}];
expt_pr = find(~cellfun(@isempty,pr_run));
nexpt_pr = length(expt_pr);
nmice_pr = length(unique(pr_mice(expt_pr)));

figure;
ind_use = intersect(resp_ind_all_phase,find(f1_all>0.02));
subplot(2,2,1)
scatter(log(f2overf1_all(ind_use)), amp_all(ind_use))
hold on
[n edges bin] = histcounts(log(f2overf1_all),[-4:0.5:1]);
amp_all_avg = zeros(length(n),2);
f2overf1_avg = zeros(length(n),2);
for i = 1:length(n)
    ind = intersect(ind_use,find(bin == i));
    amp_all_avg(i,1) = nanmean(amp_all(ind),1);
    amp_all_avg(i,2) = nanstd(amp_all(ind),[],1)./sqrt(length(ind));
    f2overf1_avg(i,1) = nanmean(log(f2overf1_all(ind)),1);
    f2overf1_avg(i,2) = nanstd(log(f2overf1_all(ind)),[],1)./sqrt(length(ind));
end
errorbar(f2overf1_avg(:,1),amp_all_avg(:,1),amp_all_avg(:,2),amp_all_avg(:,2),f2overf1_avg(:,2),f2overf1_avg(:,2),'-o')
xlabel('log(F2/F1)')
ylabel('Sine Amplitude')
xlim([-4 1])
[b_amp dev_amp stats_amp] = glmfit(f2overf1_all, amp_all);
title(['Slope = ' num2str(chop(b_amp(2),3)) '; p = ' num2str(chop(stats_amp.p(2),3))])

subplot(2,2,2)
scatter(log(f2overf1_all(ind_use)), b_all(ind_use))
hold on
b_all_avg = zeros(length(n),2);
f2overf1_avg = zeros(length(n),2);
for i = 1:length(n)
    ind = intersect(ind_use,find(bin == i));
    b_all_avg(i,1) = nanmean(b_all(ind),1);
    b_all_avg(i,2) = nanstd(b_all(ind),[],1)./sqrt(length(ind));
    f2overf1_avg(i,1) = nanmean(log(f2overf1_all(ind)),1);
    f2overf1_avg(i,2) = nanstd(log(f2overf1_all(ind)),[],1)./sqrt(length(ind));
end
errorbar(f2overf1_avg(:,1),b_all_avg(:,1),b_all_avg(:,2),b_all_avg(:,2),f2overf1_avg(:,2),f2overf1_avg(:,2),'-o')
xlabel('log(F2/F1)')
ylabel('Sine Baseline')
xlim([-4 1])
[b_base dev_base stats_base] = glmfit(f2overf1_all, b_all);
title(['Slope = ' num2str(chop(b_base(2),3)) '; p = ' num2str(chop(stats_base.p(2),3))])

subplot(2,2,3)
ind1 = intersect(ind_use,find(f2overf1_all>1));
cdfplot(amp_all(ind1));
hold on
ind2 = intersect(ind_use,find(f2overf1_all<0.5));
cdfplot(amp_all(ind2));
legend({['F2/F1>1 - n=' num2str(length(ind1))], ['F2/F1<0.5 - n=' num2str(length(ind2))]},'location','southeast')
xlabel('Sine Amplitude')
title('')

subplot(2,2,4)
ind1 = intersect(ind_use,find(f2overf1_all>1));
cdfplot(b_all(ind1));
hold on
ind2 = intersect(ind_use,find(f2overf1_all<0.5));
cdfplot(b_all(ind2));
xlabel('Sine Baseline')
title('')
suptitle(['n = ' num2str(nexpt_pr) ' expts; ' num2str(nmice_pr) ' mice; ' num2str(length(ind_use)) ' cells']) 


print(fullfile(summaryDir_F7, 'Figure7_f2f1_Summary.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F7, 'Figure7_f2f1_Summary.fig'))

nStimCon = 2;
nMaskCon = 2;
tt= (-prewin_frames:postwin_frames-1).*(1000./frame_rate);
resp_all = squeeze(mean(resp_tc_all(prewin_frames:prewin_frames+(4.*15),:,:,:,:,1),1));
ind1 = find(resp_all(:,1,2,1)>resp_all(:,2,2,1));
ind2 = find(resp_all(:,1,2,1)<resp_all(:,2,2,3));
ind = intersect(resp_ind_all_phase,intersect(ind1,ind2));
ind = [100 235 253 335 336];
figure;
movegui('center')
for i = 1:length(ind)
    max_val = resp_all(ind(i),2,2,3).*4.5;
    subplot(5,3,((i-1).*3)+1)
    shadedErrorBar(tt', resp_tc_all(:,ind(i),1,nStimCon,1,1), resp_tc_all(:,ind(i),1,nStimCon,1,2));
    ylim([-0.2 max_val])
    title(num2str(ind(i)))
    subplot(5,3,((i-1).*3)+2)
    shadedErrorBar(tt', resp_tc_all(:,ind(i),nMaskCon,nStimCon,1,1), resp_tc_all(:,ind(i),nMaskCon,nStimCon,1,2));
    ylim([-0.2 max_val])
    subplot(5,3,((i-1).*3)+3)
    shadedErrorBar(tt', resp_tc_all(:,ind(i),nMaskCon,nStimCon,3,1), resp_tc_all(:,ind(i),nMaskCon,nStimCon,3,2));
    ylim([-0.2 max_val])
end
print(fullfile(summaryDir_GA, 'GraphicAbstract_ExCells.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_GA, 'GraphicAbstract_ExCells.fig'))    
    
