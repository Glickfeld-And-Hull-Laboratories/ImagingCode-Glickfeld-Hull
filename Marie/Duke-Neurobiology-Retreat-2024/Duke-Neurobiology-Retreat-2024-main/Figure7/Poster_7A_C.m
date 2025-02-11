bmin = -3;
bmax = 2;
binsize = .01;
Izero = ((abs(bmin)-1)/binsize) - 1;
smooth = [1, 5];
clear N
clear N_pE
clear N_pL
%C = colororder;

counter = 1;
for n = 1:length(SS_DecTone)
    R = SS_DecTone(n).RecorNum;
    if [Rlist(R).PredPerc_potP] > .3
    Trials = [Rlist(R).TrialStruct];
    trigger = Trials(strcmp({Trials.TrialType}, 'b'));
    trigger = trigger(strcmp({trigger.Outcome}, 'p'));
    trigger = trigger([trigger.RTj] < -.2);
    if ~isempty(trigger)
        trigger = [trigger.JuiceTime];
        [N(counter,:), edges] = OneUnitHistStructTimeLimLineINDEX(trigger.', n, SS_DecTone, bmin, bmax, binsize, [0 inf], 4, 'k', NaN, 0, 0);
        % zscore if we want to switch back N(counter, :) = (N(counter, :) - mean(N(counter,1:Izero)))/std(N(counter,1:Izero));
      % N(counter, :) = (N(counter, :) - mean(N(counter,1:Izero)));
        counter = counter + 1;
    end
    end
end
N = N(~isnan(sum(N,2)),:);
if smooth(1) == 1
    N_pE.N = smoothdata(N, 2, 'sgolay', smooth(2));
else
    N_pE.N = N;
end
N_pE.edges = edges;
counterL = 1;
trials_c = 0;
list_trialsC = Rlist([Rlist.PredPerc_potP] > .3);
for n = 1:length(list_trialsC)
    Trials = [list_trialsC(n).TrialStruct];
        trigger = Trials(strcmp({Trials.TrialType}, 'b'));
    trigger = trigger(strcmp({trigger.Outcome}, 'p'));
    trigger = trigger([trigger.RTj] < -.2);
    if ~isempty(trigger)
        trigger = [trigger.JuiceTime];
        [N_pE(counterL).LickHist_N, N_pE(counterL).LickHist_edges] = LickHist(trigger.', list_trialsC(n).AllLicks, [-2 2], .01, 'k', 0);
        [N_pE(counterL).DeliveryHist_N, N_pE(counterL).DeliveryHist_edges] = LickHist(trigger.', list_trialsC(n).JuiceTimes, [-2 2], .02, 'k', 0);
                [N_pE(counterL).RunMean, N_pE(counterL).RunEdges, ~] = RunSpeedHistLines(trigger.', list_trialsC(n).RunningStruct.SpeedTimesAdj, list_trialsC(n).RunningStruct.SpeedValues, -2, 2);
        N_pE(counterL).trials_c = length(trigger);
        counterL = counterL + 1;
    end
end



clear N
counter = 1;
for n = 1:length(SS_DecTone)
    R = SS_DecTone(n).RecorNum;
    if [Rlist(R).PredPerc_potP] >= .3
    Trials = [Rlist(R).TrialStruct];
    trigger = Trials(strcmp({Trials.TrialType}, 'b'));
    trigger = trigger(strcmp({trigger.Outcome}, 'p'));
    trigger = trigger([trigger.RTj] > -.2);
    if ~isempty(trigger)
        trigger = [trigger.JuiceTime];
        [N(counter,:), edges] = OneUnitHistStructTimeLimLineINDEX(trigger.', n, SS_DecTone, bmin, bmax, binsize, [0 inf], 4, 'k', NaN, 0, 0);
      % N(counter, :) = (N(counter, :) - mean(N(counter,1:Izero)));
        counter = counter + 1;
    end
    end
end
if smooth(1) == 1
    N_pL.N = smoothdata(N, 2, 'sgolay', smooth(2));
else
    N_pL.N = N;
end
N_pL.edges = edges;
counterL = 1;
trials_c = 0;
list_trialsC = Rlist([Rlist.PredPerc_potP] > .3);
for n = 1:length(list_trialsC)
    Trials = [list_trialsC(n).TrialStruct];
        trigger = Trials(strcmp({Trials.TrialType}, 'b'));
    trigger = trigger(strcmp({trigger.Outcome}, 'p'));
    trigger = trigger([trigger.RTj] >= -.2);
    if ~isempty(trigger)
        trigger = [trigger.JuiceTime];
        [N_pL(counterL).LickHist_N, N_pL(counterL).LickHist_edges] = LickHist(trigger.', list_trialsC(n).AllLicks, [-2 2], .01, 'k', 0);
        [N_pL(counterL).DeliveryHist_N, N_pL(counterL).DeliveryHist_edges] = LickHist(trigger.', list_trialsC(n).JuiceTimes, [-2 2], .02, 'k', 0);
        [N_pL(counterL).RunMean, N_pL(counterL).RunEdges, tester] = RunSpeedHistLines(trigger.', list_trialsC(n).RunningStruct.SpeedTimesAdj, list_trialsC(n).RunningStruct.SpeedValues, -2, 2);
        N_pL(counterL).trials_c = length(trigger);
        counterL = counterL + 1;
    end
end

figure
hold on
shadedErrorBar2(edges(1:end-1), nanmean(cell2mat({N_pL.N}.'), 1), nanstd(cell2mat({N_pL.N}.'), 0, 1)/sqrt(size(cell2mat({N_pL.N}.'), 1)), 'LineProp', {C(9,:)});
shadedErrorBar2(edges(1:end-1), nanmean(cell2mat({N_pE.N}.'), 1), nanstd(cell2mat({N_pE.N}.'), 0, 1)/sqrt(size(cell2mat({N_pE.N}.'), 1)), 'LineProp', {C(8,:)});
% centerY = 1.4;
% plot([-.5; -.5], [centerY; centerY + .25], '-k',  [-.5; -.25], [centerY; centerY], '-k', 'LineWidth', 1)
xline(0, 'c', 'LineWidth', 1);
xline(-.68, 'g', 'LineWidth', 1);
% legend({'trials with reactive licking'; 'trials with predictive licking'; }, 'Location', 'northoutside', 'FontSize', 12);
% legend('boxoff')
xlabel('time from rewarded solenoid');
ylabel('Sspk/s');
xlim([-1 .5]);
ylim([70 125]);
title('trained click resp w_wo tone_inclPartTrain');
% FigureWrap(NaN, 'Trained_Sspk_pred_incPartTrain_Spk_DecTone_EarlyLate', NaN, NaN, NaN, NaN, 2.5, 3.5);


figure
hold on
shadedErrorBar2(N_pL(1).LickHist_edges(1:end-1), nanmean(cell2mat({N_pL.LickHist_N}.'), 1), nanstd(cell2mat({N_pL.LickHist_N}.'), 0, 1)/sqrt(size(cell2mat({N_pL.LickHist_N}.'), 1)), 'LineProp', {C(9,:)});
shadedErrorBar2(N_pE(1).LickHist_edges(1:end-1), nanmean(cell2mat({N_pE.LickHist_N}.'), 1), nanstd(cell2mat({N_pE.LickHist_N}.'), 0, 1)/sqrt(size(cell2mat({N_pE.LickHist_N}.'), 1)), 'LineProp', {C(8,:)});
% centerY = 1.4;
% plot([-.5; -.5], [centerY; centerY + .25], '-k',  [-.5; -.25], [centerY; centerY], '-k', 'LineWidth', 1)
xline(0, 'c', 'LineWidth', 1);
xline(-.68, 'g', 'LineWidth', 1);
% legend({'trials with reactive licking'; 'trials with predictive licking'; }, 'Location', 'northoutside', 'FontSize', 12);
% legend('boxoff')
xlabel('time from rewarded solenoid');
ylabel('licks/s');
xlim([-1 .5]);
%ylim([-1 7]);
title('trained click resp w/wo tone inclPartTrain_');
% FigureWrap(NaN, 'Trained_Sspk_pred_incPartTrain_deltaSpk_lickHist_DecTone_EarlyLate', NaN, NaN, NaN, NaN, 2.5, 3.5);


figure
hold on
shadedErrorBar2(N_pL(1).RunEdges.', nanmean(cell2mat({N_pL.RunMean}), 2), nanstd(cell2mat({N_pL.RunMean}), 0, 2)/sqrt(size(cell2mat({N_pL.LickHist_N}), 2)), 'LineProp', {C(9,:)});
shadedErrorBar2(N_pE(1).RunEdges.', nanmean(cell2mat({N_pE.RunMean}), 2), nanstd(cell2mat({N_pE.RunMean}), 0, 2)/sqrt(size(cell2mat({N_pE.LickHist_N}), 2)), 'LineProp', {C(8,:)});
% centerY = 1.4;
% plot([-.5; -.5], [centerY; centerY + .25], '-k',  [-.5; -.25], [centerY; centerY], '-k', 'LineWidth', 1)
xline(0, 'c', 'LineWidth', 1);
xline(-.68, 'g', 'LineWidth', 1);
% legend({'trials with reactive licking'; 'trials with predictive licking'; }, 'Location', 'northoutside', 'FontSize', 12);
% legend('boxoff')
xlabel('time from rewarded solenoid');
ylabel('m/s');
xlim([-1 .5]);
ylim([5 20]);
title('trained click resp w/wo tone inclPartTrain_');
% FigureWrap(NaN, 'Trained_Sspk_pred_incPartTrain_deltaSpk_locoSpeed_DecTone_EarlyLate', NaN, NaN, NaN, NaN, 2.5, 3.5);
