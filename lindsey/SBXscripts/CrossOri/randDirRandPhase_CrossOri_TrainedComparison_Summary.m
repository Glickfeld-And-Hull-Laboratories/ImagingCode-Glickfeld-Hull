%% Rand Dir
%30 deg trained vs naive
randDirSummary = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\CrossOri\RandDirSummary';
load([randDirSummary '\randDirPassive_Summary_V1_noDriver_SFpt05_Sz30.mat'])
resp_ind_all= setdiff(resp_ind_all, red_cells_all);
figure;
movegui('center')
subplot(2,2,1)
cdfplot(stim_OSI_all(resp_ind_all))
hold on
xlabel('OSI')
title('')
subplot(2,2,3)
cdfplot(Zc_all(resp_ind_all))
xlabel('Zc')
xlim([-2 10])
hold on
subplot(2,2,4)
cdfplot(Zp_all(resp_ind_all))
hold on
xlabel('Zp')
xlim([-2 10])
leg_str{1} = ['Trained- ' num2str(length(resp_ind_all)) ' cells'];
load([randDirSummary '\randDirRandPhase_Summary_V1_SLC_SFpt05_Sz30.mat'])
subplot(2,2,1)
cdfplot(stim_OSI_all(resp_ind_all))
hold on
xlabel('OSI')
title('')
subplot(2,2,3)
cdfplot(Zc_all(resp_ind_all))
xlabel('Zc')
xlim([-2 10])
title('')
hold on
subplot(2,2,4)
cdfplot(Zp_all(resp_ind_all))
hold on
xlabel('Zp')
xlim([-2 10])
title('')
leg_str{2} = ['Naive- ' num2str(length(resp_ind_all)) ' cells'];
subplot(2,2,1)
legend(leg_str)
sgtitle('V1 - SF 0.05; Sz 30')
print(fullfile(randDirSummary, 'randDir_TrainedVsNaive_SFpt05_Sz30.pdf'),'-dpdf')

%1000 deg trained vs naive
load([randDirSummary '\randDirPassive_Summary_V1_noDriver_SFpt05_Sz1000.mat'])
resp_ind_all= setdiff(resp_ind_all, red_cells_all);
figure;
movegui('center')
subplot(2,2,1)
cdfplot(stim_OSI_all(resp_ind_all))
hold on
xlabel('OSI')
title('')
subplot(2,2,3)
cdfplot(Zc_all(resp_ind_all))
xlabel('Zc')
xlim([-2 10])
hold on
subplot(2,2,4)
cdfplot(Zp_all(resp_ind_all))
hold on
xlabel('Zp')
xlim([-2 10])
leg_str{1} = ['Trained- ' num2str(length(resp_ind_all)) ' cells'];
load([randDirSummary '\randDirFF_Summary_V1_SLC_SFpt05_Sz1000.mat'])
subplot(2,2,1)
cdfplot(stim_OSI_all(resp_ind_all))
hold on
xlabel('OSI')
title('')
subplot(2,2,3)
cdfplot(Zc_all(resp_ind_all))
xlabel('Zc')
xlim([-2 10])
title('')
hold on
subplot(2,2,4)
cdfplot(Zp_all(resp_ind_all))
hold on
xlabel('Zp')
xlim([-2 10])
title('')
leg_str{2} = ['Naive- ' num2str(length(resp_ind_all)) ' cells'];
subplot(2,2,1)
legend(leg_str)
sgtitle('V1 - SF 0.05; Sz 1000')
print(fullfile(randDirSummary, 'randDir_TrainedVsNaive_SFpt05_Sz1000.pdf'),'-dpdf')

%% Rand Phase
%30 deg trained vs naive
randPhaseSummary = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\CrossOri\RandPhaseSummary';
load([randPhaseSummary '\randPhasePassive_Summary_V1_noDriver_SFpt05_Sz30.mat'])
resp_ind_all= setdiff(resp_ind_all, red_cells_all);
figure;
movegui('center')
subplot(2,2,1)
cdfplot(testPI_all(resp_ind_all))
hold on
subplot(2,2,2)
cdfplot(plaidSI_all(resp_ind_all))
hold on
subplot(2,2,3)
cdfplot(b_all_all(resp_ind_all))
xlim([-1 1])
hold on
subplot(2,2,4)
cdfplot(amp_all_all(resp_ind_all)-amp_shuf_all(resp_ind_all))
xlim([-0.2 1])
hold on
leg_str{1} = ['Trained- ' num2str(length(resp_ind_all)) ' cells'];
load([randPhaseSummary '\randDirRandPhase_Summary_V1_SLC_SFpt05_Sz30.mat'])
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
leg_str{2} = ['Naive- ' num2str(length(resp_ind_all)) ' cells'];
subplot(2,2,1)
legend(leg_str)
sgtitle('V1 - SF 0.05; Sz 30')
print(fullfile(randDirSummary, 'randDir_TrainedVsNaive_SFpt05_Sz30.pdf'),'-dpdf')

%1000 deg trained vs naive
load([randPhaseSummary '\randPhasePassive_Summary_V1_noDriver_SFpt05_Sz1000.mat'])
resp_ind_all= setdiff(resp_ind_all, red_cells_all);
figure;
movegui('center')
subplot(2,2,1)
cdfplot(testPI_all(resp_ind_all))
hold on
subplot(2,2,2)
cdfplot(plaidSI_all(resp_ind_all))
hold on
subplot(2,2,3)
cdfplot(b_all_all(resp_ind_all))
xlim([-1 1])
hold on
subplot(2,2,4)
cdfplot(amp_all_all(resp_ind_all)-amp_shuf_all(resp_ind_all))
xlim([-0.2 1])
hold on
leg_str{1} = ['Trained- ' num2str(length(resp_ind_all)) ' cells'];
load([randPhaseSummary '\randPhaseFF_Summary_V1_SLC_SFpt05_Sz1000.mat'])
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
leg_str{2} = ['Naive- ' num2str(length(resp_ind_all)) ' cells'];
subplot(2,2,1)
legend(leg_str)
sgtitle('V1 - SF 0.05; Sz 1000')
print(fullfile(randPhaseSummary, 'randPhase_TrainedVsNaive_SFpt05_Sz1000.pdf'),'-dpdf')