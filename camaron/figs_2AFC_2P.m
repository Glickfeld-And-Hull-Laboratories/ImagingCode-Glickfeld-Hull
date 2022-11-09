% Plot figures and return quantifications

%% Load data...
% Save/load appropriate vars - pool_2AFC_2P_data.m

mouse = 'i475'
cd(['Z:\All_Staff\home\camaron\Analysis\2P\' mouse])
filename = [mouse '_pooled_2AFC_2P_data.mat'];
load(filename)

%% Percent responsive * Modify seciton to use 'stats' struct

prcnt_adapt_sigDuringAdapt_b = b_adapt_sigDuringAdapt_nCells/nAllCells * 100;
prcnt_stim_sigDuringControl_b = b_stim_sigDuringControl_nCells/nAllCells * 100;
prcnt_stim_sigDuringAdapt_b = b_stim_sigDuringAdapt_nCells/nAllCells * 100;


prcnt_adapt_sigDuringAdapt_p = p_adapt_sigDuringAdapt_nCells/nAllCells * 100;
prcnt_stim_sigDuringControl_p = p_stim_sigDuringControl_nCells/nAllCells * 100;
prcnt_stim_sigDuringAdapt_p = p_stim_sigDuringAdapt_nCells/nAllCells * 100;


for i = 1:length(ori_bins)
    prcnt_adapt_tuned_sigDuringAdapt_b(i) = b_tuned_adapt_sigDuringAdapt_nCells(i)/nAllCells * 100;
    prcnt_stim_tuned_sigDuringControl_b(i) = b_tuned_stim_sigDuringControl_nCells(i)/nAllCells * 100;
    prcnt_stim_tuned_sigDuringAdapt_b(i) = b_tuned_stim_sigDuringAdapt_nCells(i)/nAllCells * 100;


    prcnt_adapt_tuned_sigDuringAdapt_p(i) = p_tuned_adapt_sigDuringAdapt_nCells(i)/nAllCells * 100;
    prcnt_stim_tuned_sigDuringControl_p(i) = p_tuned_stim_sigDuringControl_nCells(i)/nAllCells * 100;
    prcnt_stim_tuned_sigDuringAdapt_p(i) = p_tuned_stim_sigDuringAdapt_nCells(i)/nAllCells * 100;

end

%% plotting percent responsive

f = figure()
f.Position(3:4) = [840 630];
tiledlayout(2,3)

% figure();
% subplot(2,3,1)
ax1 = nexttile;
b = bar([prcnt_adapt_sigDuringAdapt_b; prcnt_adapt_sigDuringAdapt_p]);
ylabel('% Responsive')
title("% resp to adaptor")
xticklabels({'Active', 'Passive'})
b.FaceColor = "flat";
b.CData = [0 0.4470 0.7410; 0.6350 0.0780 0.1840];
b.FaceAlpha = 0.9;
set(gca, 'box', 'off')


% subplot(2,3,2)
ax2 = nexttile;
b = bar([prcnt_stim_sigDuringControl_b; prcnt_stim_sigDuringControl_p]);
ylabel('% Responsive')
title("% resp to target during control trials")
xticklabels({'Active', 'Passive'})
b.FaceColor = "flat";
b.CData = [0 0.4470 0.7410; 0.6350 0.0780 0.1840];
b.FaceAlpha = 0.9;
set(gca, 'box', 'off')


% subplot(2,3,3)
ax3 = nexttile;
b = bar([prcnt_stim_sigDuringAdapt_b; prcnt_stim_sigDuringAdapt_p]);
ylabel('% Responsive')
title("% resp to target during adapt trials")
xticklabels({'Active', 'Passive'})
b.FaceColor = "flat";
b.CData = [0 0.4470 0.7410; 0.6350 0.0780 0.1840];
b.FaceAlpha = 0.9;
set(gca, 'box', 'off')

linkaxes([ax1 ax2 ax3],'y')

% subplot(2,3,4)
ax4 = nexttile;
paired_adapt_tuned_sigDuringAdapt = [prcnt_adapt_tuned_sigDuringAdapt_b; prcnt_adapt_tuned_sigDuringAdapt_p]';
b = bar(paired_adapt_tuned_sigDuringAdapt);
ylabel('% Responsive')
xticks([1 2 3 4])
xticklabels({'0', '45', '90', '135'})
set(gca, 'box', 'off')


% subplot(2,3,5)
ax5 = nexttile;
paired_stim_tuned_sigDuringControl = [prcnt_stim_tuned_sigDuringControl_b; prcnt_stim_tuned_sigDuringControl_p]';
b = bar(paired_stim_tuned_sigDuringControl);
ylabel('% Responsive')
xticks([1 2 3 4])
xticklabels({'0', '45', '90', '135'})
set(gca, 'box', 'off')


% subplot(2,3,6)
ax6 = nexttile;
paired_stim_tuned_sigDuringAdapt = [prcnt_stim_tuned_sigDuringAdapt_b; prcnt_stim_tuned_sigDuringAdapt_p]';
b = bar(paired_stim_tuned_sigDuringAdapt);
ylabel('% Responsive')
xticks([1 2 3 4])
xticklabels({'0', '45', '90', '135'})
set(gca, 'box', 'off')


linkaxes([ax4 ax5 ax6],'y')

sgtitle([mouse ': ' num2str(nAllCells) ' total cells'])


% SAVE

fn = [mouse '_cell_responsiveness.pdf'];
exportgraphics(gcf,fn)

%% All cells dF/F

f = figure()
f.Position(3:4) = [840 630];
tiledlayout(1,3)

ax1 = nexttile;
b = bar([b_data_adapt_stats.mean_df; p_data_adapt_stats.mean_df])
ylabel('mean response')
title("All cells, adaptor responses")
xticklabels({'Active', 'Passive'})
b.FaceColor = "flat";
b.CData = [0 0.4470 0.7410; 0.6350 0.0780 0.1840];
b.FaceAlpha = 0.9;
set(gca, 'box', 'off')

ax2 = nexttile;
b = bar([b_data_stim_control_stats.mean_df; p_data_stim_control_stats.mean_df])
ylabel('mean response')
title("All cells, target responses; control trials")
xticklabels({'Active', 'Passive'})
b.FaceColor = "flat";
b.CData = [0 0.4470 0.7410; 0.6350 0.0780 0.1840];
b.FaceAlpha = 0.9;
set(gca, 'box', 'off')


ax3 = nexttile;
b = bar([b_data_stim_adapt_stats.mean_df; p_data_stim_adapt_stats.mean_df])
ylabel('mean response')
title("All cells, target responses, adapt trials")
xticklabels({'Active', 'Passive'})
b.FaceColor = "flat";
b.CData = [0 0.4470 0.7410; 0.6350 0.0780 0.1840];
b.FaceAlpha = 0.9;
set(gca, 'box', 'off')

linkaxes([ax1 ax2 ax3],'y')

fn = [mouse '_cell_responsiveness.pdf'];
exportgraphics(gcf,fn, 'Append', true)

%% plot traces with target, adapter vline, and respp window shaded




%% Responsive by condition

f = figure()
f.Position(3:4) = [840 630];
tiledlayout(2,2)

ax1 = nexttile;
b = bar([b_adapt_sigDuringAdapt_adapt_stats.mean_df; p_adapt_sigDuringAdapt_adapt_stats.mean_df])
ylabel('mean response')
title("Adapt win, SigDuringAdapt")
xticklabels({'Active', 'Passive'})
b.FaceColor = "flat";
b.CData = [0 0.4470 0.7410; 0.6350 0.0780 0.1840];
b.FaceAlpha = 0.9;
set(gca, 'box', 'off')

ax2 = nexttile;
b = bar([b_stim_sigDuringControl_control_stats.mean_df; p_stim_sigDuringControl_control_stats.mean_df])
ylabel('mean response')
title("Target win, SigDuringControl, control trials")
xticklabels({'Active', 'Passive'})
b.FaceColor = "flat";
b.CData = [0 0.4470 0.7410; 0.6350 0.0780 0.1840];
b.FaceAlpha = 0.9;
set(gca, 'box', 'off')


ax3 = nexttile;
b = bar([b_stim_sigDuringControl_adapt_stats.mean_df; p_stim_sigDuringControl_adapt_stats.mean_df])
ylabel('mean response')
title("Target win, SigDuringControl, adapt trials")
xticklabels({'Active', 'Passive'})
b.FaceColor = "flat";
b.CData = [0 0.4470 0.7410; 0.6350 0.0780 0.1840];
b.FaceAlpha = 0.9;
set(gca, 'box', 'off')

ax4 = nexttile;
b = bar([b_stim_sigDuringAdapt_adapt_stats.mean_df; p_stim_sigDuringAdapt_adapt_stats.mean_df])
ylabel('mean response')
title("Target win, SigDuringAdapt, adapt trials")
xticklabels({'Active', 'Passive'})
b.FaceColor = "flat";
b.CData = [0 0.4470 0.7410; 0.6350 0.0780 0.1840];
b.FaceAlpha = 0.9;
set(gca, 'box', 'off')

linkaxes([ax1 ax2 ax3 ax4],'y')


fn = [mouse '_cell_responsiveness.pdf'];
exportgraphics(gcf,fn, 'Append', true)

%% Responsive + tuned by condition 

f = figure()
f.Position(3:4) = [840 630];
tiledlayout(2,2)

for i = 1:4
b_means(i) = b_tuned_adapt_sigDuringAdapt_adapt_stats(i).mean_df;
p_means(i) = p_tuned_adapt_sigDuringAdapt_adapt_stats(i).mean_df;
end

ax1 = nexttile;
b = bar([b_means; p_means]');
ylabel('mean response')
title("Adapt win, SigDuringAdapt")
xticks([1 2 3 4])
xticklabels({'0', '45', '90', '135'})
set(gca, 'box', 'off')

%--

for i = 1:4
b_means(i) = b_tuned_stim_sigDuringControl_control_stats(i).mean_df;
p_means(i) = p_tuned_stim_sigDuringControl_control_stats(i).mean_df;
end

ax2 = nexttile;
b = bar([b_means; p_means]');
ylabel('mean response')
title("Target win, SigDuringControl, control trials")
xticks([1 2 3 4])
xticklabels({'0', '45', '90', '135'})
set(gca, 'box', 'off')

%--

for i = 1:4
b_means(i) = b_tuned_stim_sigDuringControl_adapt_stats(i).mean_df;
p_means(i) = p_tuned_stim_sigDuringControl_adapt_stats(i).mean_df;
end

ax3 = nexttile;
b = bar([b_means; p_means]');
ylabel('mean response')
title("Target win, SigDuringControl, adapt trials")
xticks([1 2 3 4])
xticklabels({'0', '45', '90', '135'})
set(gca, 'box', 'off')

%--

for i = 1:4
b_means(i) = b_tuned_stim_sigDuringAdapt_adapt_stats(i).mean_df;
p_means(i) = p_tuned_stim_sigDuringAdapt_adapt_stats(i).mean_df;
end

ax4 = nexttile;
b = bar([b_means; p_means]');
ylabel('mean response')
title("Target win, SigDuringAdapt, adapt trials")
xticks([1 2 3 4])
xticklabels({'0', '45', '90', '135'})
set(gca, 'box', 'off')

linkaxes([ax1 ax2 ax3 ax4],'y')


fn = [mouse '_cell_responsiveness.pdf'];
exportgraphics(gcf,fn, 'Append', true)




