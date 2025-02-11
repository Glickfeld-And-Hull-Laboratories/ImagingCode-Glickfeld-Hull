figure
hold on
scatter([Rlist([Rlist.day] < 7).day].', [Rlist([Rlist.day] < 7).meanRT].', 'MarkerFaceColor', C(4,:), 'MarkerEdgeColor', 'none');
eb1 = errorbar([Rlist([Rlist.day] < 7).day].', [Rlist([Rlist.day] < 7).meanRT].',[Rlist([Rlist.day] < 7).sterrRT].', 'vertical', 'LineStyle', 'none');
set(eb1, 'color', C(4,:))
scatter([Rlist([Rlist.day] == 7 | [Rlist.day] == 8 | [Rlist.day] == 9).day].', [Rlist([Rlist.day] == 7 | [Rlist.day] == 8 | [Rlist.day] == 9).meanRT].', 'MarkerFaceColor', C(1,:), 'MarkerEdgeColor', 'none');
eb2 = errorbar([Rlist([Rlist.day] == 7 | [Rlist.day] == 8 | [Rlist.day] == 9).day].', [Rlist([Rlist.day] == 7 | [Rlist.day] == 8 | [Rlist.day] == 9).meanRT].', [Rlist([Rlist.day] == 7 | [Rlist.day] == 8 | [Rlist.day] == 9).sterrRT].', 'vertical', 'LineStyle', 'none');
set(eb2, 'color', C(1,:))
scatter([Rlist([Rlist.day] >= 10).day].', [Rlist([Rlist.day] >= 10).meanRT].', 'MarkerFaceColor', C(2,:), 'MarkerEdgeColor', 'none');
eb3 = errorbar([Rlist([Rlist.day] >= 10).day].', [Rlist([Rlist.day] >= 10).meanRT].',[Rlist([Rlist.day] >= 10).sterrRT].', 'vertical', 'LineStyle', 'none');
set(eb3, 'color', C(2,:))
yline(0, 'c', 'LineWidth', 1);
yline(-.68, 'g', 'LineWidth', 1);
xlabel('day');
ylabel('reaction time from rewarded solenoid')
title('title')
ylim([-1 2]);
% FigureWrap(NaN, 'overall_Lick_RTj', NaN, NaN, NaN, NaN, 2.5, 3.5);
