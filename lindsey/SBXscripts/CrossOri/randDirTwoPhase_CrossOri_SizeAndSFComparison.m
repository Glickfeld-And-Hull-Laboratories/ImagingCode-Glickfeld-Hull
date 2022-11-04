clear all; close all; clc;
doRedChannel = 0;
ds = 'CrossOriRandDirTwoPhaseFF_ExptList';
eval(ds)
rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'RandDirSummary');
expt_group = {'randDirTwoPhase_Summary_lowSF', 'randDirTwoPhase_Summary_hiSF', 'randDirTwoPhaseFF_Summary_lowSF'};
area = 'V1';
driver = 'SLC';
sz_str = {'RF','FF'};
sf_str = {'lowSF','hiSF'};
sz_mat = [1 1 2];
sf_mat = [1 2 1];
ngroup = size(expt_group,2);
leg_str = cell(1,ngroup);
cell_n = zeros(1,ngroup);
colors = defaultPlotColors();
for i = 1:ngroup
    load(fullfile(summaryDir,[expt_group{i} '_' area '_' driver '.mat']))
    cell_n(i) = length(resp_ind_all);
    figure(1)
    subplot(2,3,1)
    cdfplot(Zc_all(1,resp_ind_all))
    hold on
    subplot(2,3,2)
    cdfplot(Zp_all(1,resp_ind_all))
    hold on
    subplot(2,3,3)
    cdfplot(Zp_all(1,resp_ind_all)-Zc_all(1,resp_ind_all))
    hold on
    subplot(2,3,4)
    cdfplot(stimOSI_all(resp_ind_all))
    hold on
    subplot(2,3,5)
    cdfplot(stimDSI_all(resp_ind_all))
    hold on
    subplot(2,3,6)
    p1 = cdfplot(plaid_corr_all(resp_ind_all));
    p1.LineStyle = '-';
    p1.Color = colors(i,:);
    hold on
    p2 = cdfplot(plaid_corr_rand_all(resp_ind_all));
    p2.LineStyle = ':';
    p2.Color = colors(i,:);
    hold on
    leg_str{i} = [sz_str{sz_mat(i)} ' ' sf_str{sf_mat(i)} '- n=' num2str(cell_n(i))];
end
figure(1)
subplot(2,3,1)
xlabel('Zc')
title('')
legend(leg_str,'location','southeast')
subplot(2,3,2)
xlabel('Zp')
title('')
subplot(2,3,3)
xlabel('Zp-Zc')
title('')
subplot(2,3,4)
xlabel('OSI')
title('')
subplot(2,3,5)
xlabel('DSI')
title('')
subplot(2,3,6)
xlabel('Phase Corr')
legend({'Across phase', 'Within phase'})
title('')
sgtitle(['RandDir Two Phase - Size and SF'])
print(fullfile(summaryDir, ['RandDirTwoPhase_SizeAndSFSummary_' area '_' driver '.pdf']),'-dpdf', '-fillpage')