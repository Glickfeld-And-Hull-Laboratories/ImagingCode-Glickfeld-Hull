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

%% 1. plotting percent responsive

f = figure();
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

%% 2. plot traces with target, adapter vline, and resp window shaded

% All cells
f = figure();
f.Position(3:4) = [1050 525];
tiledlayout(1,3)

ax1 = nexttile;
plot_adapt_TC_bp(b_data_adapt_stats, p_data_adapt_stats, resp_win)
title("All cells, adaptor responses")

ax2 = nexttile;
plot_targ_TC_bp(b_data_stim_adapt_stats, p_data_stim_adapt_stats, resp_win)
title("All cells, target responses, adapt trials")

ax3 = nexttile;
plot_targ_TC_bp(b_data_stim_control_stats, p_data_stim_control_stats, resp_win)
title("All cells, target responses; control trials")


linkaxes([ax1 ax2 ax3],'y')

fn = [mouse '_cell_responsiveness.pdf'];
exportgraphics(gcf,fn, 'Append', true)


%% 3. All cells dF/F

f = figure();
f.Position(3:4) = [840 630];
tiledlayout(1,3)

ax1 = nexttile;
data = [b_data_adapt_stats.mean_delta_df; p_data_adapt_stats.mean_delta_df];
sem = [b_data_adapt_stats.sem_delta_df; p_data_adapt_stats.sem_delta_df];
barplot_w_sem(data, sem)
ylabel('avg change in amplitude')
title("All cells, adaptor responses")
xticklabels({'Active', 'Passive'})

ax2 = nexttile;
data = [b_data_stim_adapt_stats.mean_delta_df; p_data_stim_adapt_stats.mean_delta_df];
sem = [b_data_stim_adapt_stats.sem_delta_df; p_data_stim_adapt_stats.sem_delta_df];
barplot_w_sem(data, sem)
ylabel('avg change in amplitude')
title("All cells, target responses, adapt trials")
xticklabels({'Active', 'Passive'})


ax3 = nexttile;
data = [b_data_stim_control_stats.mean_delta_df; p_data_stim_control_stats.mean_delta_df];
sem = [b_data_stim_control_stats.sem_delta_df; p_data_stim_control_stats.sem_delta_df];
barplot_w_sem(data, sem)
ylabel('avg change in amplitude')
title("All cells, target responses; control trials")
xticklabels({'Active', 'Passive'})


linkaxes([ax1 ax2 ax3],'y')

fn = [mouse '_cell_responsiveness.pdf'];
exportgraphics(gcf,fn, 'Append', true)

%% 4. Adapt resp traces


% Responsive + tuned by condition, Cells sig responsive to adaptor
f = figure();
% f.Position(3:4) = [840 840];
tiledlayout(4,2)

ax1 = nexttile(1,[2,2]);
plot_adapt_TC_bp(b_adapt_sigDuringAdapt_adapt_stats, p_adapt_sigDuringAdapt_adapt_stats, resp_win)
title("Adapt win, SigDuringAdapt")

% 
for i = 1:4
% subplot(2,2,i)
nexttile
plot_adapt_TC_bp(b_tuned_adapt_sigDuringAdapt_adapt_stats(i), p_tuned_adapt_sigDuringAdapt_adapt_stats(i), resp_win)
title(num2str(ori_bins(i)))
end

% h = findobj('type', 'Axes');
% linkaxes([h(1), h(2), h(3), h(4), h(5)], 'y')

fn = [mouse '_cell_responsiveness.pdf'];
exportgraphics(gcf,fn, 'Append', true)


%% 5. Adaptor index (B / (A+B))

f = figure();
f.Position(3:4) = [840 630];
tiledlayout(2,1)

ax1 = nexttile;
data = [b_adapt_sigDuringAdapt_adapt_stats.Aix_mean; p_adapt_sigDuringAdapt_adapt_stats.Aix_mean];
sem = [b_adapt_sigDuringAdapt_adapt_stats.Aix_sem; p_adapt_sigDuringAdapt_adapt_stats.Aix_sem];
barplot_w_sem(data, sem)
ylabel('mean Aix (B/(A+B))')
title("Adapt win, SigDuringAdapt")
xticklabels({'Active', 'Passive'})


ax2 = nexttile;
for i = 1:4
b_means(i) = b_tuned_adapt_sigDuringAdapt_adapt_stats(i).Aix_mean;
p_means(i) = p_tuned_adapt_sigDuringAdapt_adapt_stats(i).Aix_mean;

b_sems(i) = b_tuned_adapt_sigDuringAdapt_adapt_stats(i).Aix_sem;
p_sems(i) = p_tuned_adapt_sigDuringAdapt_adapt_stats(i).Aix_sem;
end

data = [b_means; p_means]';
sem = [b_sems; p_sems]';

barplot_w_sem(data, sem)
ylabel('mean Aix (B/(A+B))')
title("Adapt win, SigDuringAdapt")
xticks([1 2 3 4])
xticklabels({'0', '45', '90', '135'})
set(gca, 'box', 'off')

fn = [mouse '_cell_responsiveness.pdf'];
exportgraphics(gcf,fn, 'Append', true)

%% 6. targ resp sig during adapt, adapt trials

% Responsive + tuned by condition, Cells sig responsive to adaptor
f = figure();
f.Position(3:4) = [840 840];
tiledlayout(4,2)

ax1 = nexttile(1,[2,2]);
plot_targ_TC_bp(b_stim_sigDuringAdapt_adapt_stats, p_stim_sigDuringAdapt_adapt_stats, resp_win)
title("Targ win, SigDuringAdapt, adapt trials")

% 
for i = 1:4
% subplot(2,2,i)
nexttile
plot_targ_TC_bp(b_tuned_stim_sigDuringAdapt_adapt_stats(i), p_tuned_stim_sigDuringAdapt_adapt_stats(i), resp_win)
title(num2str(ori_bins(i)))
end

% h = findobj('type', 'Axes');
% linkaxes([h(1), h(2), h(3), h(4), h(5)], 'y')

fn = [mouse '_cell_responsiveness.pdf'];
exportgraphics(gcf,fn, 'Append', true)

%% 7. targ resp sig during adapt, control trials

% Responsive + tuned by condition, Cells sig responsive to adaptor
f = figure();
f.Position(3:4) = [840 840];
tiledlayout(4,2)

ax1 = nexttile(1,[2,2]);
plot_targ_TC_bp(b_stim_sigDuringAdapt_control_stats, p_stim_sigDuringAdapt_control_stats, resp_win)
title("Targ win, SigDuringAdapt, control trials")

% 
for i = 1:4
% subplot(2,2,i)
nexttile
plot_targ_TC_bp(b_tuned_stim_sigDuringAdapt_control_stats(i), p_tuned_stim_sigDuringAdapt_control_stats(i), resp_win)
title(num2str(ori_bins(i)))
end

% h = findobj('type', 'Axes');
% linkaxes([h(1), h(2), h(3), h(4), h(5)], 'y')

fn = [mouse '_cell_responsiveness.pdf'];
exportgraphics(gcf,fn, 'Append', true)

%% 8. Responsive + tuned by condition, Cells sig responsive to adaptor


f = figure();
f.Position(3:4) = [840 630];
tiledlayout(2,3)

ax1 = nexttile;
data = [b_adapt_sigDuringAdapt_adapt_stats.mean_delta_df; p_adapt_sigDuringAdapt_adapt_stats.mean_delta_df];
sem = [b_adapt_sigDuringAdapt_adapt_stats.sem_delta_df; p_adapt_sigDuringAdapt_adapt_stats.sem_delta_df];
barplot_w_sem(data, sem)
ylabel('avg change in amplitude')
title("Adapt win, SigDuringAdapt")
xticklabels({'Active', 'Passive'})


ax2 = nexttile;
data = [b_stim_sigDuringAdapt_adapt_stats.mean_delta_df; p_stim_sigDuringAdapt_adapt_stats.mean_delta_df];
sem = [b_stim_sigDuringAdapt_adapt_stats.sem_delta_df; p_stim_sigDuringAdapt_adapt_stats.sem_delta_df];
barplot_w_sem(data, sem)
ylabel('avg change in amplitude')
title("Target win, SigDuringAdapt, adapt trials")
xticklabels({'Active', 'Passive'})


ax3 = nexttile;
data = [b_stim_sigDuringAdapt_control_stats.mean_delta_df; p_stim_sigDuringAdapt_control_stats.mean_delta_df];
sem = [b_stim_sigDuringAdapt_control_stats.sem_delta_df; p_stim_sigDuringAdapt_control_stats.sem_delta_df];
barplot_w_sem(data, sem)
ylabel('avg change in amplitude')
title("Target win, SigDuringAdapt, control trials")
xticklabels({'Active', 'Passive'})


ax4 = nexttile;
for i = 1:4
b_means(i) = b_tuned_adapt_sigDuringAdapt_adapt_stats(i).mean_delta_df;
p_means(i) = p_tuned_adapt_sigDuringAdapt_adapt_stats(i).mean_delta_df;

b_sems(i) = b_tuned_adapt_sigDuringAdapt_adapt_stats(i).sem_delta_df;
p_sems(i) = p_tuned_adapt_sigDuringAdapt_adapt_stats(i).sem_delta_df;
end

data = [b_means; p_means]';
sem = [b_sems; p_sems]';
barplot_w_sem(data, sem)
ylabel('avg change in amplitude')
title("Adapt win, SigDuringAdapt")
xticks([1 2 3 4])
xticklabels({'0', '45', '90', '135'})


ax5 = nexttile;
for i = 1:4
b_means(i) = b_tuned_stim_sigDuringAdapt_adapt_stats(i).mean_delta_df;
p_means(i) = p_tuned_stim_sigDuringAdapt_adapt_stats(i).mean_delta_df;

b_sems(i) = b_tuned_stim_sigDuringAdapt_adapt_stats(i).sem_delta_df;
p_sems(i) = p_tuned_stim_sigDuringAdapt_adapt_stats(i).sem_delta_df;
end

data = [b_means; p_means]';
sem = [b_sems; p_sems]';
barplot_w_sem(data, sem)
ylabel('avg change in amplitude')
title("Target win, SigDuringAdapt, adapt trials")
xticks([1 2 3 4])
xticklabels({'0', '45', '90', '135'})

ax6 = nexttile;
for i = 1:4
b_means(i) = b_tuned_stim_sigDuringAdapt_control_stats(i).mean_delta_df;
p_means(i) = p_tuned_stim_sigDuringAdapt_control_stats(i).mean_delta_df;

b_sems(i) = b_tuned_stim_sigDuringAdapt_control_stats(i).sem_delta_df;
p_sems(i) = p_tuned_stim_sigDuringAdapt_control_stats(i).sem_delta_df;
end

data = [b_means; p_means]';
sem = [b_sems; p_sems]';
barplot_w_sem(data, sem)
ylabel('avg change in amplitude')
title("Target win, SigDuringAdapt, control trials")
xticks([1 2 3 4])
xticklabels({'0', '45', '90', '135'})



linkaxes([ax1 ax2 ax3],'y')
linkaxes([ax4 ax5 ax6],'y')



fn = [mouse '_cell_responsiveness.pdf'];
exportgraphics(gcf,fn, 'Append', true)


%% 9. targ resp sig during control, adapt trials

% Responsive + tuned by condition, Cells sig responsive to adaptor
f = figure();
f.Position(3:4) = [840 840];
tiledlayout(4,2)

ax1 = nexttile(1,[2,2]);
plot_targ_TC_bp(b_stim_sigDuringControl_adapt_stats, p_stim_sigDuringControl_adapt_stats, resp_win)
title("Targ win, SigDuringControl, Adapt Trials")

% 
for i = 1:4
% subplot(2,2,i)
nexttile
plot_targ_TC_bp(b_tuned_stim_sigDuringControl_adapt_stats(i), p_tuned_stim_sigDuringControl_adapt_stats(i), resp_win)
title(num2str(ori_bins(i)))
end

% h = findobj('type', 'Axes');
% linkaxes([h(1), h(2), h(3), h(4), h(5)], 'y')

fn = [mouse '_cell_responsiveness.pdf'];
exportgraphics(gcf,fn, 'Append', true)

%% 10. targ resp sig during control, control trials

% Responsive + tuned by condition, Cells sig responsive to adaptor
f = figure();
f.Position(3:4) = [840 840];
tiledlayout(4,2)

ax1 = nexttile(1,[2,2]);
plot_targ_TC_bp(b_stim_sigDuringControl_control_stats, p_stim_sigDuringControl_control_stats, resp_win)
title("Targ win, SigDuringControl, Control Trials")

% 
for i = 1:4
% subplot(2,2,i)
nexttile
plot_targ_TC_bp(b_tuned_stim_sigDuringControl_control_stats(i), p_tuned_stim_sigDuringControl_control_stats(i), resp_win)
title(num2str(ori_bins(i)))
end

% h = findobj('type', 'Axes');
% linkaxes([h(1), h(2), h(3), h(4), h(5)], 'y')

fn = [mouse '_cell_responsiveness.pdf'];
exportgraphics(gcf,fn, 'Append', true)







%% 11. Responsive + tuned by condition  Cells sig responsive target

f = figure();
f.Position(3:4) = [840 630];
tiledlayout(2,2)



ax1 = nexttile;
data = [b_stim_sigDuringControl_adapt_stats.mean_delta_df; p_stim_sigDuringControl_adapt_stats.mean_delta_df];
sem = [b_stim_sigDuringControl_adapt_stats.sem_delta_df; p_stim_sigDuringControl_adapt_stats.sem_delta_df];
barplot_w_sem(data, sem)
ylabel('avg change in amplitude')
title("Target win, SigDuringControl, adapt trials")
xticklabels({'Active', 'Passive'})


%--
ax2 = nexttile;
data = [b_stim_sigDuringControl_control_stats.mean_delta_df; p_stim_sigDuringControl_control_stats.mean_delta_df];
sem = [b_stim_sigDuringControl_control_stats.sem_delta_df; p_stim_sigDuringControl_control_stats.sem_delta_df];
barplot_w_sem(data, sem)
ylabel('avg change in amplitude')
title("Target win, SigDuringControl, control trials")
xticklabels({'Active', 'Passive'})


%--


ax3 = nexttile;
for i = 1:4
b_means(i) = b_tuned_stim_sigDuringControl_adapt_stats(i).mean_delta_df;
p_means(i) = p_tuned_stim_sigDuringControl_adapt_stats(i).mean_delta_df;

b_sems(i) = b_tuned_stim_sigDuringControl_adapt_stats(i).sem_delta_df;
p_sems(i) = p_tuned_stim_sigDuringControl_adapt_stats(i).sem_delta_df;
end

data = [b_means; p_means]';
sem = [b_sems; p_sems]';
barplot_w_sem(data, sem)
ylabel('avg change in amplitude')
title("Target win, SigDuringControl, adapt trials")
xticks([1 2 3 4])
xticklabels({'0', '45', '90', '135'})




%--


ax4 = nexttile;
for i = 1:4
b_means(i) = b_tuned_stim_sigDuringControl_control_stats(i).mean_delta_df;
p_means(i) = p_tuned_stim_sigDuringControl_control_stats(i).mean_delta_df;

b_sems(i) = b_tuned_stim_sigDuringControl_control_stats(i).sem_delta_df;
p_sems(i) = p_tuned_stim_sigDuringControl_control_stats(i).sem_delta_df;
end

data = [b_means; p_means]';
sem = [b_sems; p_sems]';
barplot_w_sem(data, sem)
ylabel('avg change in amplitude')
title("Target win, SigDuringControl, control trials")
xticks([1 2 3 4])
xticklabels({'0', '45', '90', '135'})





linkaxes([ax1 ax2],'y')
linkaxes([ax3 ax4],'y')


fn = [mouse '_cell_responsiveness.pdf'];
exportgraphics(gcf,fn, 'Append', true)

%% 12. Responsive interneruons

f = figure();
f.Position(3:4) = [840 840];
tiledlayout(3,2)


ax1 = nexttile(1,[1,2]);
plot_adapt_TC_bp(b_adapt_sigDuringAdapt_adapt_interneurons_stats, p_adapt_sigDuringAdapt_adapt_interneurons_stats, resp_win)
title("Adapt win, INT: SigDuringAdapt")

ax2 = nexttile;
plot_targ_TC_bp(b_stim_sigDuringControl_control_interneurons_stats, p_stim_sigDuringControl_control_interneurons_stats, resp_win)
title("Stim win, INT: SigDuringControl, control trials")

ax3 = nexttile;
plot_targ_TC_bp(b_stim_sigDuringControl_adapt_interneurons_stats, p_stim_sigDuringControl_adapt_interneurons_stats, resp_win)
title("Stim win, INT: SigDuringControl, adapt trials")

ax4 = nexttile;
plot_targ_TC_bp(b_stim_sigDuringAdapt_control_interneurons_stats, p_stim_sigDuringAdapt_control_interneurons_stats, resp_win)
title("Stim win, INT: SigDuringAdapt, control trials,")

ax5 = nexttile;
plot_targ_TC_bp(b_stim_sigDuringAdapt_adapt_interneurons_stats, p_stim_sigDuringAdapt_adapt_interneurons_stats, resp_win)
title("Stim win, INT: SigDuringAdapt, adapt trials")

fn = [mouse '_cell_responsiveness.pdf'];
exportgraphics(gcf,fn, 'Append', true)


%% 13. Responsive pyramidal

f = figure();
f.Position(3:4) = [840 840];
tiledlayout(3,2)


ax1 = nexttile(1,[1,2]);
plot_adapt_TC_bp(b_adapt_sigDuringAdapt_adapt_pyramidal_stats, p_adapt_sigDuringAdapt_adapt_pyramidal_stats, resp_win)
title("Adapt win, Pyr: SigDuringAdapt")

ax2 = nexttile;
plot_targ_TC_bp(b_stim_sigDuringControl_control_pyramidal_stats, p_stim_sigDuringControl_control_pyramidal_stats, resp_win)
title("Stim win, Pyr: SigDuringControl, control trials")

ax3 = nexttile;
plot_targ_TC_bp(b_stim_sigDuringControl_adapt_pyramidal_stats, p_stim_sigDuringControl_adapt_pyramidal_stats, resp_win)
title("Stim win, Pyr: SigDuringControl, adapt trials")

ax4 = nexttile;
plot_targ_TC_bp(b_stim_sigDuringAdapt_control_pyramidal_stats, p_stim_sigDuringAdapt_control_pyramidal_stats, resp_win)
title("Stim win, Pyr: SigDuringAdapt, control trials,")

ax5 = nexttile;
plot_targ_TC_bp(b_stim_sigDuringAdapt_adapt_pyramidal_stats, p_stim_sigDuringAdapt_adapt_pyramidal_stats, resp_win)
title("Stim win, Pyr: SigDuringAdapt, adapt trials")

linkaxes([ax2 ax3],'y')
linkaxes([ax4 ax5],'y')
 
fn = [mouse '_cell_responsiveness.pdf'];
exportgraphics(gcf,fn, 'Append', true)

%% 14. Adaptor responses across cell type

f = figure();
f.Position(3:4) = [840 630];
tiledlayout(2,2)

ax1 = nexttile;
data = [b_adapt_sigDuringAdapt_adapt_interneurons_stats.mean_delta_df; p_adapt_sigDuringAdapt_adapt_interneurons_stats.mean_delta_df];
sem = [b_adapt_sigDuringAdapt_adapt_interneurons_stats.sem_delta_df; p_adapt_sigDuringAdapt_adapt_interneurons_stats.sem_delta_df];
barplot_w_sem(data, sem)
ylabel('avg change in amplitude')
title("INT: Adapt win, SigDuringAdapt")
xticklabels({'Active', 'Passive'})

ax2 = nexttile;
data = [b_adapt_sigDuringAdapt_adapt_pyramidal_stats.mean_delta_df; p_adapt_sigDuringAdapt_adapt_pyramidal_stats.mean_delta_df];
sem = [b_adapt_sigDuringAdapt_adapt_pyramidal_stats.sem_delta_df; p_adapt_sigDuringAdapt_adapt_pyramidal_stats.sem_delta_df];
barplot_w_sem(data, sem)
ylabel('avg change in amplitude')
title("Pyr: Adapt win, SigDuringAdapt")
xticklabels({'Active', 'Passive'})

ax3 = nexttile;
data = [b_adapt_sigDuringAdapt_adapt_interneurons_stats.Aix_mean; p_adapt_sigDuringAdapt_adapt_interneurons_stats.Aix_mean];
sem = [b_adapt_sigDuringAdapt_adapt_interneurons_stats.Aix_sem; p_adapt_sigDuringAdapt_adapt_interneurons_stats.Aix_sem];
barplot_w_sem(data, sem)
ylabel('mean Aix (B/(A+B))')
title("INT: Adapt win, SigDuringAdapt")
xticklabels({'Active', 'Passive'})

ax4 = nexttile;
data = [b_adapt_sigDuringAdapt_adapt_pyramidal_stats.Aix_mean; p_adapt_sigDuringAdapt_adapt_pyramidal_stats.Aix_mean];
sem = [b_adapt_sigDuringAdapt_adapt_pyramidal_stats.Aix_sem; p_adapt_sigDuringAdapt_adapt_pyramidal_stats.Aix_sem];
barplot_w_sem(data, sem)
ylabel('mean Aix (B/(A+B))')
title("Pyr: Adapt win, SigDuringAdapt")
xticklabels({'Active', 'Passive'})

linkaxes([ax1 ax2],'y')
linkaxes([ax3 ax4],'y')

fn = [mouse '_cell_responsiveness.pdf'];
exportgraphics(gcf,fn, 'Append', true)



%% Local Functions

function plot_adapt_TC_bp(b_TC_stats, p_TC_stats, resp_win)

shadedErrorBar([], b_TC_stats.mean_TC , b_TC_stats.sem_TC, 'lineProps', 'b')
hold on
shadedErrorBar([], p_TC_stats.mean_TC , p_TC_stats.sem_TC, 'lineProps', 'r')

adaptor_vline = [20,31,42,53];
target_vline = 64;

xline(adaptor_vline, '--k')
xline(target_vline, '--g')
set(gca, 'box', 'off')
ylabel('amplitude')
xlabel('Time')

base_win = [17 18 19 20 21];

t = 1:length(b_TC_stats.mean_TC);
ptchidx = (t >= resp_win(1)) & (t <= resp_win(end));     % Area To Shade
ptchidx2 = (t >= base_win(1)) & (t <= base_win(end));     % Area To Shade

y = gca().YLim;
y1 = ones(size(find(ptchidx == 1))) .* y(1);
y2 = ones(size(find(ptchidx == 1))) .* y(2);
patch([t(ptchidx) fliplr(t(ptchidx))], [y1, fliplr(y2)], 'k', 'FaceAlpha', 0.15, 'Edgecolor', 'none') 
patch([t(ptchidx2) fliplr(t(ptchidx2))], [y1, fliplr(y2)], 'k', 'FaceAlpha', 0.15, 'Edgecolor', 'none') 
hold off

end

function plot_targ_TC_bp(b_TC_stats, p_TC_stats, resp_win)

shadedErrorBar([], b_TC_stats.mean_TC , b_TC_stats.sem_TC, 'lineProps', 'b')
hold on
shadedErrorBar([], p_TC_stats.mean_TC , p_TC_stats.sem_TC, 'lineProps', 'r')

target_vline = 20;

xline(target_vline, '--g')
set(gca, 'box', 'off')
ylabel('amplitude')
xlabel('Time')

base_win = [17 18 19 20 21];


t = 1:length(b_TC_stats.mean_TC);
ptchidx = (t >= resp_win(1)) & (t <= resp_win(end));     % Area To Shade
ptchidx2 = (t >= base_win(1)) & (t <= base_win(end));     % Area To Shade
y = gca().YLim;
y1 = ones(size(find(ptchidx == 1))) .* y(1);
y2 = ones(size(find(ptchidx == 1))) .* y(2);
patch([t(ptchidx) fliplr(t(ptchidx))], [y1, fliplr(y2)], 'k', 'FaceAlpha', 0.15, 'Edgecolor', 'none') 
patch([t(ptchidx2) fliplr(t(ptchidx2))], [y1, fliplr(y2)], 'k', 'FaceAlpha', 0.15, 'Edgecolor', 'none') 
hold off

end

function barplot_w_sem(data, sem)

b = bar(data, 'grouped');
hold on 

[ngroups,nbars] = size(data);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',data,sem,'k','linestyle','none');
hold off
set(gca, 'box', 'off')


if nbars == 1 % bandaid for plotting issue
    b.FaceColor = "flat";
    b.CData = [0 0.4470 0.7410; 0.6350 0.0780 0.1840];
    b.FaceAlpha = 0.9;
end

end
