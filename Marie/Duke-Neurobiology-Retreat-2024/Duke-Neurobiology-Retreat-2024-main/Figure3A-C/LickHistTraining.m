bmin = -3;
bmax = 2;
binsize = .1;
Izero = ((abs(bmin)-1)/binsize) - 1;
smooth = [1, 5];
clear N_naive
clear N_hab
clear N_expert

counterL = 1;
list_trialsC = Rlist([Rlist.day] < 4);
for n = 1:length(list_trialsC)
    Trials = [list_trialsC(n).TrialStruct.JuiceTime];
    trigger = Trials(~isnan(Trials));
    if ~isempty(trigger)
        [N_naive(counterL).LickHist_N, N_naive(counterL).LickHist_edges] = LickHist(trigger.', list_trialsC(n).AllLicks, [-2 2], binsize, 'k', 0);
        [N_naive(counterL).DeliveryHist_N, N_naive(counterL).DeliveryHist_edges] = LickHist(trigger.', list_trialsC(n).JuiceTimes, [-2 2], binsize, 'k', 0);
        [N_naive(counterL).RunMean, N_naive(counterL).RunEdges, tester] = RunSpeedHistLines(trigger.', list_trialsC(n).RunningStruct.SpeedTimesAdj, list_trialsC(n).RunningStruct.SpeedValues, -2, 2);
        N_naive(counterL).trials_c = length(trigger);
        counterL = counterL + 1;
    end
end

counterL = 1;
list_trialsC = Rlist([Rlist.TrainBoo] == 1);
for n = 1:length(list_trialsC)
   Trials = [list_trialsC(n).TrialStruct.JuiceTime];
    trigger = Trials(~isnan(Trials));
    if ~isempty(trigger)
        [N_expert(counterL).LickHist_N, N_expert(counterL).LickHist_edges] = LickHist(trigger.', list_trialsC(n).AllLicks, [-2 2], binsize, 'k', 0);
        [N_expert(counterL).DeliveryHist_N, N_expert(counterL).DeliveryHist_edges] = LickHist(trigger.', list_trialsC(n).JuiceTimes, [-2 2], binsize, 'k', 0);
                [N_expert(counterL).RunMean, N_expert(counterL).RunEdges, ~] = RunSpeedHistLines(trigger.', list_trialsC(n).RunningStruct.SpeedTimesAdj, list_trialsC(n).RunningStruct.SpeedValues, -2, 2);
        N_expert(counterL).trials_c = length(trigger);
        counterL = counterL + 1;
    end
end

counterL = 1;
list_trialsC = Rlist([Rlist.day] == 7);
for n = 1:length(list_trialsC)
    Trials = [list_trialsC(n).TrialStruct.JuiceTime];
    trigger = Trials(~isnan(Trials));
    if ~isempty(trigger)
        [N_hab(counterL).LickHist_N, N_hab(counterL).LickHist_edges] = LickHist(trigger.', list_trialsC(n).AllLicks, [-2 2], binsize, 'k', 0);
        [N_hab(counterL).DeliveryHist_N, N_hab(counterL).DeliveryHist_edges] = LickHist(trigger.', list_trialsC(n).JuiceTimes, [-2 2], binsize, 'k', 0);
        [N_hab(counterL).RunMean, N_hab(counterL).RunEdges, tester] = RunSpeedHistLines(trigger.', list_trialsC(n).RunningStruct.SpeedTimesAdj, list_trialsC(n).RunningStruct.SpeedValues, -2, 2);
        N_hab(counterL).trials_c = length(trigger);
        counterL = counterL + 1;
    end
end

figure
hold on
shadedErrorBar2(N_naive(1).LickHist_edges(1:end-1), nanmean(cell2mat({N_naive.LickHist_N}.'), 1), nanstd(cell2mat({N_naive.LickHist_N}.'), 0, 1)/sqrt(size(cell2mat({N_naive.LickHist_N}.'), 1)), 'LineProp', {C(4,:)});
shadedErrorBar2(N_hab(1).LickHist_edges(1:end-1), nanmean(cell2mat({N_hab.LickHist_N}.'), 1), nanstd(cell2mat({N_hab.LickHist_N}.'), 0, 1)/sqrt(size(cell2mat({N_hab.LickHist_N}.'), 1)), 'LineProp', {C(1,:)});
shadedErrorBar2(N_expert(1).LickHist_edges(1:end-1), nanmean(cell2mat({N_expert.LickHist_N}.'), 1), nanstd(cell2mat({N_expert.LickHist_N}.'), 0, 1)/sqrt(size(cell2mat({N_expert.LickHist_N}.'), 1)), 'LineProp', {C(2,:)});
xline(0, 'Color', C(6,:), 'LineStyle', '-', 'LineWidth', 1);
xline(-.68, 'Color', C(5,:), 'LineStyle', '--', 'LineWidth', 1);
xlabel('time from rewarded solenoid');
ylabel('licks/s');
xlim([-1 .5]);
ylim([0 7]);
% legend({'naive'; 'trials with early predictive licking'; 'trials with late reactive licking'; 'trials with late predictive licking';}, 'Location', 'northoutside', 'FontSize', 12);
% legend('boxoff')
xlabel('time from rewarded solenoid');
ylabel('licks/s');

 % FigureWrap(NaN, 'LickHistTraining', NaN, NaN, NaN, NaN, 2.5, 3.5);

