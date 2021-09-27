close all
clear all
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
CrossOri_base = fullfile('Analysis', '2P', 'CrossOri');
PV_summaryDir = fullfile(LG_base, CrossOri_base, 'RandDirRandPhaseSummary','randDirRandPhase');
SLC_summaryDir = fullfile(LG_base, CrossOri_base, 'RandPhaseSummary','randPhase');
summaryDir_F8 = fullfile(LG_base, CrossOri_base, 'CrossOri_Figures', 'CrossOri_Figure8');
driver_list = {'SLC'; 'PV'};
ndriver = length(driver_list);
leg_str = cell(1,ndriver);
area ='V1';

for iD = 1:length(driver_list)
    driver = driver_list{iD};
    summaryDir = eval([driver_list{iD} '_summaryDir']);
    load([summaryDir '_Summary_' area '_' driver '.mat'])
    driverSummary(iD).name = driver_list{iD};
    driverSummary(iD).mice = unique(mouse_list,'rows');
    driverSummary(iD).nmice = size(unique(mouse_list,'rows'),1);
    driverSummary(iD).nexp = size(mouse_list,1);
    driverSummary(iD).totCells = size(phase_SI_all,1);
    driverSummary(iD).respCells_phase = length(resp_ind_all_phase);
    
    ind_phase = resp_ind_all_phase(find(~isnan(b_all(resp_ind_all_phase))));
    leg_str{1,iD} = [driver_list{iD} '- ' num2str(length(ind_phase)) ' cells; ' num2str(driverSummary(iD).nmice) ' mice'];
    figure(1)
    subplot(2,3,1)
    cdfplot(phase_SI_all(ind_phase))
    ylabel('Fraction of cells')
    xlabel('Stimulus selectivity index')
    title('')
    phase_SI_driver{iD} = phase_SI_all(ind_phase);
    [val dim] = max(size(phase_SI_all));
    phase_SI_all_avg(iD,1) = mean(phase_SI_all(ind_phase));
    phase_SI_all_avg(iD,2) = std(phase_SI_all(ind_phase),[],dim)./sqrt(length(ind_phase));
    hold on
    subplot(2,3,2)
    cdfplot(b_all(ind_phase))
    ylabel('Fraction of cells')
    xlabel('Sinusoid baseline')
    title('')
    b_driver{iD} = b_all(ind_phase);
    b_all_avg(iD,1) = mean(b_all(ind_phase));
    b_all_avg(iD,2) = std(b_all(ind_phase),[],dim)./sqrt(length(ind_phase));
    hold on
    subplot(2,3,3)
    cdfplot(amp_all(ind_phase))
    ylabel('Fraction of cells')
    xlabel('Sinusoid amplitude')
    title('')
    amp_driver{iD} = amp_all(ind_phase);
    amp_all_avg(iD,1) = mean(amp_all(ind_phase));
    amp_all_avg(iD,2) = std(amp_all(ind_phase),[],dim)./sqrt(length(ind_phase));
    amp_all_avg(iD,3) = median(amp_all(ind_phase));
    hold on
    subplot(2,3,4)
    cdfplot(amp_all(ind_phase)-amp_shuf_all(ind_phase))
    ylabel('Fraction of cells')
    xlabel({'Sinusoid amplitude','(-shuffle)'})
    title('')
    ampshuf_driver{iD} = amp_all(ind_phase)-amp_shuf_all(ind_phase);
    ampshuf_all_avg(iD,1) = mean(amp_all(ind_phase)-amp_shuf_all(ind_phase));
    ampshuf_all_avg(iD,2) = std(amp_all(ind_phase)-amp_shuf_all(ind_phase),[],dim)./sqrt(length(ind_phase));
    hold on
    subplot(2,3,5)
    cdfplot(phase_MI_all(ind_phase))
    ylabel('Fraction of cells')
    xlabel({'Masking index','test+mask vs plaid'})
    title('')
    phase_MI_driver{iD} = phase_MI_all(ind_phase);
    phase_MI_all_avg(iD,1) = mean(phase_MI_all(ind_phase));
    phase_MI_all_avg(iD,2) = std(phase_MI_all(ind_phase),[],dim)./sqrt(length(ind_phase));
    hold on
    subplot(2,3,6)
    cdfplot(phase_MI_max_all(ind_phase))
    ylabel('Fraction of cells')
    xlabel({'Masking index','test vs plaid'})
    title('')
    phase_MI_max_driver{iD} = phase_MI_max_all(ind_phase);
    phase_MI_max_all_avg(iD,1) = mean(phase_MI_max_all(ind_phase));
    phase_MI_max_all_avg(iD,2) = std(phase_MI_max_all(ind_phase),[],dim)./sqrt(length(ind_phase));
    hold on
    
    figure(2)
    subplot(2,2,1+iD-1)
    histogram(b_all(ind_phase),[-1:0.1:1])
    hold on
    ylabel('Number of cells')
    xlabel('Sinusoid baseline')
    title(driver)
    subplot(2,2,3+iD-1)
    histogram(amp_all(ind_phase),[0:0.05:1])
    hold on
    ylabel('Number of cells')
    xlabel('Sinusoid amplitude')
    title(driver)
    
    figure(3)
    cdfplot(stim_OSI_all(ind_phase))
    hold on
    stim_OSI_avg(iD,1) = nanmean(stim_OSI_all(ind_phase));
    stim_OSI_avg(iD,2) = nanstd(stim_OSI_all(ind_phase),[],dim)./sqrt(length(ind_phase));
    stim_OSI_avg(iD,3) = sum(~isnan(stim_OSI_all(ind_phase)));
    leg_str2{1,iD} = [driver_list{iD} '- ' num2str(sum(~isnan(stim_OSI_all(ind_phase)))) ' cells'];
    stim_OSI_driver{iD} = stim_OSI_all(ind_phase);
    anova_tot(iD,1) = sum(anova_all(ind_phase)<0.05);
    anova_tot(iD,2) = length(anova_all(ind_phase));
    if iD == 2
        z_avg(1) = mean(z_all);
        z_avg(2) = std(z_all)./sqrt(length(z_all));
    end
end
figure(1)
subplot(2,3,1)
legend(leg_str)
movegui('center')
print(fullfile(summaryDir_F8, 'Figure8_PV&SLC_PhaseSummaryCdfs.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_PV&SLC_PhaseSummaryCdfs.fig'))

figure(2)
movegui('center')
print(fullfile(summaryDir_F8, 'Figure8_PV&SLC_AmpBaseSummaryHists.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_PV&SLC_AmpBaseSummaryHists.fig'))

figure(3)
title('')
ylabel('Fraction of cells')
xlabel('OSI')
legend(leg_str2)
movegui('center')
print(fullfile(summaryDir_F8, 'Figure8_PV_OSI.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_PV_OSI.fig'))

p_stim_OSI = ranksum(stim_OSI_driver{1},stim_OSI_driver{2});
p_phase_SI = ranksum(phase_SI_driver{1},phase_SI_driver{2});
p_b = ranksum(b_driver{1},b_driver{2});
p_amp = ranksum(amp_driver{1},amp_driver{2});
p_ampshuf = ranksum(ampshuf_driver{1},ampshuf_driver{2});
p_MI = ranksum(phase_MI_driver{1},phase_MI_driver{2});
[h p_anova chi] = prop_test([22 468],[108 814],true);
fract_mod = anova_tot(2)./length(ind_phase);

%% example cells
mouse = 'i1340';
date = '210226';
run_str  = 'runs-003';
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']),'npSub_tc');
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']));
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']))

nCells = size(npSub_tc,2);
nTrials = length(maskPhas_all);
resp_ind = unique([resptest_ind; respmask_ind]);
nOn = frame_rate.*4;
postwin_frames = frame_rate.*6;
data_trial = nan(prewin_frames+postwin_frames,nCells,nTrials);
for iTrial = 1:nTrials-1
    data_trial(:,:,iTrial) = npSub_tc(cStimOn(iTrial)-prewin_frames:cStimOn(iTrial)+postwin_frames-1,:);
end
data_f = mean(data_trial(1:prewin_frames,:,:),1);
data_df = data_trial-data_f;
data_dfof = data_df./data_f;
tt= (-prewin_frames:postwin_frames-1).*(1000./frame_rate);
resp_tc = zeros(prewin_frames+postwin_frames,nCells,nMaskCon,nStimCon,nMaskPhas,2);
for im = 1:nMaskCon
    ind_m = find(maskCon_all == maskCons(im));
    for is = 1:nStimCon
        ind_s = find(stimCon_all == stimCons(is));
        if im>1 & is>1
            for ip = 1:nMaskPhas
                ind_p = find(maskPhas_all == maskPhas(ip));
                ind_use = intersect(find(centroid_dist<2), intersect(ind_p,intersect(ind_m,ind_s)));
                resp_tc(:,:,im,is,ip,1) = nanmean(data_dfof(:,:,ind_use),3);
                resp_tc(:,:,im,is,ip,2) = nanstd(data_dfof(:,:,ind_use),[],3)./sqrt(length(ind_use));
           end
        else
            ind_use = intersect(find(centroid_dist<2), intersect(ind_m,ind_s));
            resp_tc(:,:,im,is,1,1) = nanmean(data_dfof(:,:,ind_use),3);
            resp_tc(:,:,im,is,1,2) = nanstd(data_dfof(:,:,ind_use),[],3)./sqrt(length(ind_use));
        end
    end
end

figure;
imagesc(data_avg);
colormap gray
title([mouse ' ' date])
movegui('center')
clim([100 1200])
print(fullfile(summaryDir_F8, 'Figure8_PVexampleFOV.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_PVexampleFOV.fig'))

ex_ind = [14 20 21 24];
%ex_ind = [19 22 36 39];
for iC = 1:length(ex_ind)
    figure;
    iCell = ex_ind(iC);
    subplot(4,3,1)
    shadedErrorBar(tt, resp_tc(:,iCell,1,2,1,1),resp_tc(:,iCell,1,2,1,2));
    ylim([-0.1 max(max(max(max(resp_tc(:,iCell,:,:,:,1),[],1),[],3),[],4),[],5).*1.2])
    title('test')
    subplot(4,3,2)
    shadedErrorBar(tt, resp_tc(:,iCell,2,1,1,1),resp_tc(:,iCell,2,1,1,2));
    ylim([-0.1 max(max(max(max(resp_tc(:,iCell,:,:,:,1),[],1),[],3),[],4),[],5).*1.2])
    title('mask')
    for ip = 1:nMaskPhas
        subplot(4,3,2+ip)
        shadedErrorBar(tt, resp_tc(:,iCell,2,2,ip,1,1),resp_tc(:,iCell,2,2,ip,2));
        ylim([-0.1 max(max(max(max(resp_tc(:,iCell,:,:,:,1),[],1),[],3),[],4),[],5).*1.2])
        title(num2str(maskPhas(ip)))
    end
    errorbar(maskPhas, SI_all_avg(iCell,:,1),SI_all_avg(iCell,:,2),'o');
    hold on
    plot(1:360,yfit_all(iCell,:))
    ylim([-1 1])
    title(['Amp- ' num2str(chop(amp_hat_all(iCell),2)) '; Base- ' num2str(chop(b_hat_all(iCell),2))])
    suptitle([mouse ' ' date ' cell #' num2str(iCell)])
    movegui('center')
    print(fullfile(summaryDir_F8, ['Figure8_PVexamplePhaseResp_Cell' num2str(iCell) '.pdf']),'-dpdf','-bestfit')
    savefig(fullfile(summaryDir_F8, ['Figure8_PVexamplePhaseResp_Cell' num2str(iCell) '.fig']))
end

