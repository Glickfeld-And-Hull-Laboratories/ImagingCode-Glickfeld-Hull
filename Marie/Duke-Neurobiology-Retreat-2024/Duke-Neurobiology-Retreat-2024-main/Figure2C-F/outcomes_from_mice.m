N_p = [];
N_r = [];
N_o = [];
N_b = [];
N_all = [];
%C = colororder;
These_Rlist = Rlist([Rlist.day] == 7 | [Rlist.day] == 8 | [Rlist.day] == 9);
for R = 1:length(These_Rlist)
    DayTrials = These_Rlist(R).TrialStruct;
    DayTrials = DayTrials(strcmp({DayTrials.TrialType}, 'b'));
    N_all = [N_all; DayTrials];
    Trials = DayTrials(strcmp({DayTrials.Outcome}, 'p'));
    if ~isempty(Trials)
        N_p = [N_p; Trials];
    end
    Trials = DayTrials(strcmp({DayTrials.Outcome}, 'r'));
    if ~isempty(Trials)
        N_r = [N_r; Trials];
    end
    Trials = DayTrials(strcmp({DayTrials.Outcome}, 'o'));
    if ~isempty(Trials)
        N_o = [N_o; Trials];
    end
    Trials = DayTrials(strcmp({DayTrials.Outcome}, 'b'));
    if ~isempty(Trials)
        N_b = [N_b; Trials];
    end
end
        [h_o, edges] = histcounts([N_o.RTj], [-2.5:.05:2.5]);
        [h_r, edges] = histcounts([N_r.RTj] , [-2.5:.05:2.5]);
        [h_p, edges] = histcounts([N_p.RTj] , [-2.5:.05:2.5]);
        [h_b, edges] = histcounts([N_b.RTj] , [-2.5:.05:2.5]);
        figure
        hold on
bar(edges(1:end-1), h_b/length([N_all]), 'FaceColor', C(3,:), 'EdgeColor', 'none', 'BarWidth', 1);
bar(edges(1:end-1), h_o/length([N_all]), 'FaceColor', C(3,:), 'EdgeColor', 'none', 'BarWidth', 1);
bar(edges(1:end-1), h_r/length([N_all]), 'FaceColor', C(2,:), 'EdgeColor', 'none', 'BarWidth', 1);
bar(edges(1:end-1), h_p/length([N_all]), 'FaceColor', C(1,:), 'EdgeColor', 'none', 'BarWidth', 1);
xline(-.682, 'g', 'LineWidth', 1);
xline(0, 'c', 'LineWidth', 1);
xlim([-1 0.5])
title(['mouse:' num2str(These_Rlist(R).mouse) ' day:'  num2str(These_Rlist(R).day)]);
xlabel('time from reward (s)');
ylabel('liklihood of lick onset');
FormatFigure(NaN, NaN);
% legend({'Predict'; 'React'; 'Outside Trial'; 'cue'; 'reward'})
% legend('boxoff')
ylim([0 .15]);
title('naive');
% FigureWrap(NaN, ['outcomes_from_naive_mice'], NaN, NaN, NaN, NaN, 2.5, 3.5);

N_p = [];
N_r = [];
N_o = [];
N_b = [];
N_all = [];
%C = colororder;
These_Rlist = Rlist([Rlist.day] >= 10);
for R = 1:length(These_Rlist)
    DayTrials = These_Rlist(R).TrialStruct;
    DayTrials = DayTrials(strcmp({DayTrials.TrialType}, 'b'));
    N_all = [N_all; DayTrials];
    Trials = DayTrials(strcmp({DayTrials.Outcome}, 'p'));
    if ~isempty(Trials)
        N_p = [N_p; Trials];
    end
    Trials = DayTrials(strcmp({DayTrials.Outcome}, 'r'));
    if ~isempty(Trials)
        N_r = [N_r; Trials];
    end
    Trials = DayTrials(strcmp({DayTrials.Outcome}, 'o'));
    if ~isempty(Trials)
        N_o = [N_o; Trials];
    end
    Trials = DayTrials(strcmp({DayTrials.Outcome}, 'b'));
    if ~isempty(Trials)
        N_b = [N_b; Trials];
    end
end
        [h_o, edges] = histcounts([N_o.RTj], [-2.5:.05:2.5]);
        [h_r, edges] = histcounts([N_r.RTj] , [-2.5:.05:2.5]);
        [h_p, edges] = histcounts([N_p.RTj] , [-2.5:.05:2.5]);
        [h_b, edges] = histcounts([N_b.RTj] , [-2.5:.05:2.5]);
        figure
        hold on
bar(edges(1:end-1), h_b/length([N_all]), 'FaceColor', C(3,:), 'EdgeColor', 'none', 'BarWidth', 1);
bar(edges(1:end-1), h_o/length([N_all]), 'FaceColor', C(3,:), 'EdgeColor', 'none', 'BarWidth', 1);
bar(edges(1:end-1), h_r/length([N_all]), 'FaceColor', C(2,:), 'EdgeColor', 'none', 'BarWidth', 1);
bar(edges(1:end-1), h_p/length([N_all]), 'FaceColor', C(1,:), 'EdgeColor', 'none', 'BarWidth', 1);
xline(-.682, 'g', 'LineWidth', 1);
xline(0, 'c', 'LineWidth', 1);
xlim([-1 0.5])
title(['mouse:' num2str(These_Rlist(R).mouse) ' day:'  num2str(These_Rlist(R).day)]);
xlabel('time from reward (s)');
ylabel('liklihood of lick onset');
FormatFigure(NaN, NaN);
% legend({'Predict'; 'React'; 'Outside Trial'; 'cue'; 'reward'})
% legend('boxoff')
ylim([0 .15]);
title('training');
% FigureWrap(NaN, ['outcomes from training mice'], NaN, NaN, NaN, NaN, 2.5, 3.5);

N_p = [];
N_r = [];
N_o = [];
N_b = [];
N_all = [];
%C = colororder;
These_Rlist = Rlist([Rlist.day] <= 3);
for R = 1:length(These_Rlist)
    DayTrials = These_Rlist(R).TrialStruct;
    DayTrials = DayTrials(strcmp({DayTrials.TrialType}, 'j'));
    N_all = [N_all; DayTrials];
    Trials = DayTrials(strcmp({DayTrials.Outcome}, 'p'));
    if ~isempty(Trials)
        N_p = [N_p; Trials];
    end
    Trials = DayTrials(strcmp({DayTrials.Outcome}, 'r'));
    if ~isempty(Trials)
        N_r = [N_r; Trials];
    end
    Trials = DayTrials(strcmp({DayTrials.Outcome}, 'o'));
    if ~isempty(Trials)
        N_o = [N_o; Trials];
    end
    Trials = DayTrials(strcmp({DayTrials.Outcome}, 'b'));
    if ~isempty(Trials)
        N_b = [N_b; Trials];
    end
end
        [h_o, edges] = histcounts([N_o.RTj], [-2.5:.05:2.5]);
        [h_r, edges] = histcounts([N_r.RTj] , [-2.5:.05:2.5]);
        [h_p, edges] = histcounts([N_p.RTj] , [-2.5:.05:2.5]);
        [h_b, edges] = histcounts([N_b.RTj] , [-2.5:.05:2.5]);
        figure
        hold on
bar(edges(1:end-1), h_b/length([N_all]), 'FaceColor', C(3,:), 'EdgeColor', 'none', 'BarWidth', 1);
bar(edges(1:end-1), h_o/length([N_all]), 'FaceColor', C(3,:), 'EdgeColor', 'none', 'BarWidth', 1);
bar(edges(1:end-1), h_r/length([N_all]), 'FaceColor', C(2,:), 'EdgeColor', 'none', 'BarWidth', 1);
bar(edges(1:end-1), h_p/length([N_all]), 'FaceColor', C(1,:), 'EdgeColor', 'none', 'BarWidth', 1);
%xline(-.682, 'g', 'LineWidth', 1);
xline(0, 'c', 'LineWidth', 1);
xlim([-1 0.5])
title(['mouse:' num2str(These_Rlist(R).mouse) ' day:'  num2str(These_Rlist(R).day)]);
xlabel('time from reward (s)');
ylabel('liklihood of lick onset');
FormatFigure(NaN, NaN);
% legend({'Predict'; 'React'; 'Outside Trial'; 'cue'; 'reward'})
% legend('boxoff')
title('naive 1-3');
ylim([0 .15]);
% FigureWrap(NaN, ['outcomes_from_naive1_3mice'], NaN, NaN, NaN, NaN, 2.5, 3.5);

N_p = [];
N_r = [];
N_o = [];
N_b = [];
N_all = [];
%C = colororder;
These_Rlist = Rlist([Rlist.TrainBoo] == 1);
for R = 1:length(These_Rlist)
    DayTrials = These_Rlist(R).TrialStruct;
    DayTrials = DayTrials(strcmp({DayTrials.TrialType}, 'b'));
    N_all = [N_all; DayTrials];
    Trials = DayTrials(strcmp({DayTrials.Outcome}, 'p'));
    if ~isempty(Trials)
        N_p = [N_p; Trials];
    end
    Trials = DayTrials(strcmp({DayTrials.Outcome}, 'r'));
    if ~isempty(Trials)
        N_r = [N_r; Trials];
    end
    Trials = DayTrials(strcmp({DayTrials.Outcome}, 'o'));
    if ~isempty(Trials)
        N_o = [N_o; Trials];
    end
    Trials = DayTrials(strcmp({DayTrials.Outcome}, 'b'));
    if ~isempty(Trials)
        N_b = [N_b; Trials];
    end
end
        [h_o, edges] = histcounts([N_o.RTj], [-2.5:.05:2.5]);
        [h_r, edges] = histcounts([N_r.RTj] , [-2.5:.05:2.5]);
        [h_p, edges] = histcounts([N_p.RTj] , [-2.5:.05:2.5]);
        [h_b, edges] = histcounts([N_b.RTj] , [-2.5:.05:2.5]);
        figure
        hold on
bar(edges(1:end-1), h_b/length([N_all]), 'FaceColor', C(3,:), 'EdgeColor', 'none', 'BarWidth', 1);
bar(edges(1:end-1), h_o/length([N_all]), 'FaceColor', C(3,:), 'EdgeColor', 'none', 'BarWidth', 1);
bar(edges(1:end-1), h_r/length([N_all]), 'FaceColor', C(2,:), 'EdgeColor', 'none', 'BarWidth', 1);
bar(edges(1:end-1), h_p/length([N_all]), 'FaceColor', C(1,:), 'EdgeColor', 'none', 'BarWidth', 1);
xline(-.682, 'g', 'LineWidth', 1);
xline(0, 'c', 'LineWidth', 1);
xlim([-1 0.5])
title(['mouse:' num2str(These_Rlist(R).mouse) ' day:'  num2str(These_Rlist(R).day)]);
xlabel('time from reward (s)');
ylabel('liklihood of lick onset');
FormatFigure(NaN, NaN);
% legend({'Predict'; 'React'; 'Outside Trial'; 'cue'; 'reward'})
% legend('boxoff')
title('TrainBoo');
ylim([0 .15]);
% FigureWrap(NaN, ['outcomes_from_TrainBoo_mice'], NaN, NaN, NaN, NaN, 2.5, 3.5);
