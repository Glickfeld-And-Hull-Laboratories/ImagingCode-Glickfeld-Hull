%% Load data...
% Save/load appropriate vars - pool_2AFC_2P_data.m

mouse = 'i475'
cd(['Z:\All_Staff\home\camaron\Analysis\2P\' mouse])
filename = [mouse '_pooled_2AFC_2P_data.mat'];
load(filename)

%% NEW Segment Plots

%% 1. plot traces with target, adapter vline, and resp window shaded

% All cells
f = figure(1);
f.Position(3:4) = [1050 525];
tiledlayout(1,3)

ax1 = nexttile;
plot_adapt_TC_bp(b_adapt_stats, p_adapt_stats, resp_win)
title("Adaptor responses")

ax2 = nexttile;
plot_targ_TC_bp(b_targ_adapt_stats, p_targ_adapt_stats, resp_win)
title("Target responses, adapt trials")

ax3 = nexttile;
plot_targ_TC_bp(b_targ_control_stats, p_targ_control_stats, resp_win)
title("Target responses; control trials")


linkaxes([ax1 ax2 ax3],'y')

fn = [mouse '_cell_responsiveness_new_segment.pdf'];
exportgraphics(gcf,fn)
close(f)

%% 2. dF/F

f = figure(2);
f.Position(3:4) = [840 630];
tiledlayout(1,3)

ax1 = nexttile;
data = [b_adapt_stats.mean_delta_df; p_adapt_stats.mean_delta_df];
sem = [b_adapt_stats.sem_delta_df; p_adapt_stats.sem_delta_df];
barplot_w_sem(data, sem)
ylabel('avg change in amplitude')
title("Adaptor responses")
xticklabels({'Active', 'Passive'})

ax2 = nexttile;
data = [b_targ_adapt_stats.mean_delta_df; p_targ_adapt_stats.mean_delta_df];
sem = [b_targ_adapt_stats.sem_delta_df; p_targ_adapt_stats.sem_delta_df];
barplot_w_sem(data, sem)
ylabel('avg change in amplitude')
title("Target responses, adapt trials")
xticklabels({'Active', 'Passive'})


ax3 = nexttile;
data = [b_targ_control_stats.mean_delta_df; p_targ_control_stats.mean_delta_df];
sem = [b_targ_control_stats.sem_delta_df; p_targ_control_stats.sem_delta_df];
barplot_w_sem(data, sem)
ylabel('avg change in amplitude')
title("Target responses; control trials")
xticklabels({'Active', 'Passive'})


linkaxes([ax1 ax2 ax3],'y')

fn = [mouse '_cell_responsiveness_new_segment.pdf'];
exportgraphics(gcf,fn, 'Append', true)
close(f)


%% 3. plot traces with target, adapter vline, and resp window shaded

% % All cells
% f = figure();
% f.Position(3:4) = [1050 525];
% % tiledlayout(1,3)
% 
% ax1 = nexttile;
% plot_adapt_TC_bp(b_adapt_AIX_stats, p_adapt_AIX_stats, resp_win)
% title("Adaptor responses - adapt resp cells only")
% 
% fn = [mouse '_cell_responsiveness_new_segment.pdf'];
% exportgraphics(gcf,fn, 'Append', true)

%% 3. Adaptor index (B / (A+B))

f = figure(3);
f.Position(3:4) = [840 630];

data = [b_adapt_stats.Aix_mean; p_adapt_stats.Aix_mean];
sem = [b_adapt_stats.Aix_sem; p_adapt_stats.Aix_sem];
barplot_w_sem(data, sem)
ylabel('mean Aix (B/(A+B))')
title("AIX(B/(A+B)")
xticklabels({'Active', 'Passive'})


fn = [mouse '_cell_responsiveness_new_segment.pdf'];
exportgraphics(gcf,fn, 'Append', true)

close(f)

%% 4. Distribition of Adaptor index (B / (A+B))

f = figure(4);
f.Position(3:4) = [840 630];
% tiledlayout(1,2)

% nexttile
data = [b_adapt_stats.Aix_by_cell; p_adapt_stats.Aix_by_cell];
scatter(data(1,:), data(2,:))
ylabel('passive')
title("Scatter AIX (B/(A+B))")
xlabel('behaving')
refline(1)

fn = [mouse '_cell_responsiveness_new_segment.pdf'];
exportgraphics(gcf,fn, 'Append', true)
close(f)


%% 5. Adaptor index (B /A))

f = figure(5);
f.Position(3:4) = [840 630];

data = [b_adapt_stats.Aix_alt_mean; p_adapt_stats.Aix_alt_mean];
sem = [b_adapt_stats.Aix_alt_sem; p_adapt_stats.Aix_alt_sem];
barplot_w_sem(data, sem)
ylabel('mean Aix (B/A)')
title("AIX(B/A)")
xticklabels({'Active', 'Passive'})


fn = [mouse '_cell_responsiveness_new_segment.pdf'];
exportgraphics(gcf,fn, 'Append', true)
close(f)


%% 6. Distribition of Adaptor index (B /A)

f = figure(6);
f.Position(3:4) = [840 630];
% tiledlayout(1,2)

% nexttile
data = [b_adapt_stats.Aix_alt; p_adapt_stats.Aix_alt];
scatter(data(1,:), data(2,:))
ylabel('passive')
title("Scatter AIX (B/A)")
xlabel('behaving')
refline(1)

fn = [mouse '_cell_responsiveness_new_segment.pdf'];
exportgraphics(gcf,fn, 'Append', true)
close(f)


%% 7. Target Time courses (adapt, control)

f = figure(7);
f.Position(3:4) = [840 630];
tiledlayout(2, 6)

for i = 1:length(b_target_ori_data_adapt_stats)
    nexttile
    plot_targ_TC_bp(b_target_ori_data_adapt_stats(i), p_target_ori_data_adapt_stats(i), resp_win)
    title(tOris(i))
end


for i = 1:length(b_target_ori_data_adapt_stats)
    nexttile
    plot_targ_TC_bp(b_target_ori_data_control_stats(i), p_target_ori_data_control_stats(i), resp_win)
end


h = findobj('type', 'Axes');
linkaxes([h(1), h(2), h(3), h(4), h(5), h(6),...
    h(7), h(8), h(9), h(10), h(11), h(12)], 'y')

sgtitle(['Target response, adapt v control, ' num2str(b_target_ori_data_adapt_stats(1).count_cells) ' cells'])



fn = [mouse '_cell_responsiveness_new_segment.pdf'];
exportgraphics(gcf,fn, 'Append', true)
close(f)


%% 8. 

f = figure(8);
f.Position(3:4) = [840 630];
tiledlayout(2, 6)

for i = 1:length(b_target_ori_data_adapt_stats)
    nexttile;
    data = [b_target_ori_data_adapt_stats(i).mean_delta_df; p_target_ori_data_adapt_stats(i).mean_delta_df];
    sem = [b_target_ori_data_adapt_stats(i).sem_delta_df; p_target_ori_data_adapt_stats(i).sem_delta_df];
    barplot_w_sem(data, sem)
    if i == 1
        ylabel('avg change in amplitude')
    %     title("Target responses, adapt trials")
    end
    xticklabels({'Active', 'Passive'})

    title(tOris(i))

end

for i = 1:length(b_target_ori_data_adapt_stats)
    nexttile;
    data = [b_target_ori_data_control_stats(i).mean_delta_df; p_target_ori_data_control_stats(i).mean_delta_df];
    sem = [b_target_ori_data_control_stats(i).sem_delta_df; p_target_ori_data_control_stats(i).sem_delta_df];
    barplot_w_sem(data, sem)
    if i == 1
        ylabel('avg change in amplitude')
    %     title("Target responses, control trials")
    end
    xticklabels({'Active', 'Passive'})
end

h = findobj('type', 'Axes');
linkaxes([h(1), h(2), h(3), h(4), h(5), h(6),...
    h(7), h(8), h(9), h(10), h(11), h(12)], 'y')

sgtitle(['Target response, adapt v control, ' num2str(b_target_ori_data_adapt_stats(1).count_cells) ' cells'])

fn = [mouse '_cell_responsiveness_new_segment.pdf'];
exportgraphics(gcf,fn, 'Append', true)
close(f)


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
