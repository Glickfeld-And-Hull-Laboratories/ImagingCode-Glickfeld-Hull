bmin = -3;
bmax = 2;
binsize = .01;
Izero = ((abs(bmin)-1)/binsize) - 1;
smooth = [1, 5];
clear N
clear N_r_Dec
clear N_o_Dec
clear N_r_Inc
clear N_o_Inc
clear N_p_Inc
clear N_p_Dec
%C = colororder;

clear N
counter = 1;
for n = 1:length(SS_DecTone)
    R = SS_DecTone(n).RecorNum;
    if [Rlist(R).PredPerc_potP] > .3
    Trials = [Rlist(R).TrialStruct];
    trigger = Trials(strcmp({Trials.TrialType}, 'b'));
    trigger = trigger(strcmp({trigger.Outcome}, 'r'));
    if ~isempty(trigger)
        trigger = [trigger.RTj_realTime];
        [N(counter,:), edges] = OneUnitHistStructTimeLimLineINDEX(trigger.', n, SS_DecTone, bmin, bmax, binsize, [0 inf], 4, 'k', NaN, 0, 0);
        % zscore if we want to switch back N(counter, :) = (N(counter, :) - mean(N(counter,1:Izero)))/std(N(counter,1:Izero));
      % N(counter, :) = (N(counter, :) - mean(N(counter,1:Izero)));
        counter = counter + 1;
    end
    end
end
N = N(~isnan(sum(N,2)),:);
if smooth(1) == 1
    N_r_Dec.N = smoothdata(N, 2, 'sgolay', smooth(2));
else
    N_r_Dec.N = N;
end
N_r_Dec.edges = edges;
counterL = 1;
trials_c = 0;
list_trialsC = Rlist([Rlist.PredPerc_potP] > .3);
for n = 1:length(list_trialsC)
    Trials = [list_trialsC(n).TrialStruct];
        trigger = Trials(strcmp({Trials.TrialType}, 'b'));
    trigger = trigger(strcmp({trigger.Outcome}, 'r'));
    if ~isempty(trigger)
        trigger = [trigger.RTj_realTime];
        [N_r_Dec(counterL).LickHist_N, N_r_Dec(counterL).LickHist_edges] = LickHist(trigger.', list_trialsC(n).AllLicks, [-2 2], .01, 'k', 0);
        [N_r_Dec(counterL).DeliveryHist_N, N_r_Dec(counterL).DeliveryHist_edges] = LickHist(trigger.', list_trialsC(n).JuiceTimes, [-2 2], .02, 'k', 0);
                [N_r_Dec(counterL).RunMean, N_r_Dec(counterL).RunEdges, ~] = RunSpeedHistLines(trigger.', list_trialsC(n).RunningStruct.SpeedTimesAdj, list_trialsC(n).RunningStruct.SpeedValues, -2, 2);
        N_r_Dec(counterL).trials_c = length(trigger);
        counterL = counterL + 1;
    end
end



clear N
counter = 1;
for n = 1:length(SS_DecTone)
    R = SS_DecTone(n).RecorNum;
    if [Rlist(R).PredPerc_potP] > .3
    Trials = [Rlist(R).LickOnsets];
    trigger = Trials(strcmp({Trials.Outcome}, 'o'));
    trigger = trigger([trigger.RewardBoo] == 0);
    if ~isempty(trigger)
        trigger = [trigger.time];
        [N(counter,:), edges] = OneUnitHistStructTimeLimLineINDEX(trigger.', n, SS_DecTone, bmin, bmax, binsize, [0 inf], 4, 'k', NaN, 0, 0);
      % N(counter, :) = (N(counter, :) - mean(N(counter,1:Izero)));
        counter = counter + 1;
    end
    end
end
if smooth(1) == 1
    N_o_Dec.N = smoothdata(N, 2, 'sgolay', smooth(2));
else
    N_o_Dec.N = N;
end
N_o_Dec.edges = edges;
counterL = 1;
trials_c = 0;
list_trialsC = Rlist([Rlist.PredPerc_potP] > .3);
for n = 1:length(list_trialsC)
    Trials = [list_trialsC(n).LickOnsets];
        trigger = Trials(strcmp({Trials.Outcome}, 'o'));
    trigger = trigger([trigger.RewardBoo] == 0);
    if ~isempty(trigger)
        trigger = [trigger.time];
        [N_o_Dec(counterL).LickHist_N, N_o_Dec(counterL).LickHist_edges] = LickHist(trigger.', list_trialsC(n).AllLicks, [-2 2], .01, 'k', 0);
        [N_o_Dec(counterL).DeliveryHist_N, N_o_Dec(counterL).DeliveryHist_edges] = LickHist(trigger.', list_trialsC(n).JuiceTimes, [-2 2], .02, 'k', 0);
        [N_o_Dec(counterL).RunMean, N_o_Dec(counterL).RunEdges, tester] = RunSpeedHistLines(trigger.', list_trialsC(n).RunningStruct.SpeedTimesAdj, list_trialsC(n).RunningStruct.SpeedValues, -2, 2);
        N_o_Dec(counterL).trials_c = length(trigger);
        counterL = counterL + 1;
    end
end

clear N
counter = 1;
for n = 1:length(SS_IncTone)
    R = SS_IncTone(n).RecorNum;
    if [Rlist(R).PredPerc_potP] > .3
    Trials = [Rlist(R).TrialStruct];
    trigger = Trials(strcmp({Trials.TrialType}, 'b'));
    trigger = trigger(strcmp({trigger.Outcome}, 'r'));
    if ~isempty(trigger)
        trigger = [trigger.RTj_realTime];
        [N(counter,:), edges] = OneUnitHistStructTimeLimLineINDEX(trigger.', n, SS_IncTone, bmin, bmax, binsize, [0 inf], 4, 'k', NaN, 0, 0);
        % zscore if we want to switch back N(counter, :) = (N(counter, :) - mean(N(counter,1:Izero)))/std(N(counter,1:Izero));
      % N(counter, :) = (N(counter, :) - mean(N(counter,1:Izero)));
        counter = counter + 1;
    end
    end
end
N = N(~isnan(sum(N,2)),:);
if smooth(1) == 1
    N_r_Inc.N = smoothdata(N, 2, 'sgolay', smooth(2));
else
    N_r_Inc.N = N;
end
N_r_Inc.edges = edges;
counterL = 1;
trials_c = 0;
list_trialsC = Rlist([Rlist.PredPerc_potP] > .3);
for n = 1:length(list_trialsC)
    Trials = [list_trialsC(n).TrialStruct];
        trigger = Trials(strcmp({Trials.TrialType}, 'b'));
    trigger = trigger(strcmp({trigger.Outcome}, 'r'));
    if ~isempty(trigger)
        trigger = [trigger.RTj_realTime];
        [N_r_Inc(counterL).LickHist_N, N_r_Inc(counterL).LickHist_edges] = LickHist(trigger.', list_trialsC(n).AllLicks, [-2 2], .01, 'k', 0);
        [N_r_Inc(counterL).DeliveryHist_N, N_r_Inc(counterL).DeliveryHist_edges] = LickHist(trigger.', list_trialsC(n).JuiceTimes, [-2 2], .02, 'k', 0);
                [N_r_Inc(counterL).RunMean, N_r_Inc(counterL).RunEdges, ~] = RunSpeedHistLines(trigger.', list_trialsC(n).RunningStruct.SpeedTimesAdj, list_trialsC(n).RunningStruct.SpeedValues, -2, 2);
        N_r_Inc(counterL).trials_c = length(trigger);
        counterL = counterL + 1;
    end
end

clear N
counter = 1;
for n = 1:length(SS_DecTone)
    R = SS_DecTone(n).RecorNum;
    if [Rlist(R).PredPerc_potP] > .3
    Trials = [Rlist(R).TrialStruct];
    trigger = Trials(strcmp({Trials.TrialType}, 'b'));
    trigger = trigger(strcmp({trigger.Outcome}, 'p'));
    if ~isempty(trigger)
        trigger = [trigger.RTj_realTime];
        [N(counter,:), edges] = OneUnitHistStructTimeLimLineINDEX(trigger.', n, SS_DecTone, bmin, bmax, binsize, [0 inf], 4, 'k', NaN, 0, 0);
        % zscore if we want to switch back N(counter, :) = (N(counter, :) - mean(N(counter,1:Izero)))/std(N(counter,1:Izero));
      % N(counter, :) = (N(counter, :) - mean(N(counter,1:Izero)));
        counter = counter + 1;
    end
    end
end
N = N(~isnan(sum(N,2)),:);
if smooth(1) == 1
    N_p_Dec.N = smoothdata(N, 2, 'sgolay', smooth(2));
else
    N_p_Dec.N = N;
end
N_p_Dec.edges = edges;
counterL = 1;
trials_c = 0;
list_trialsC = Rlist([Rlist.PredPerc_potP] > .3);
for n = 1:length(list_trialsC)
    Trials = [list_trialsC(n).TrialStruct];
        trigger = Trials(strcmp({Trials.TrialType}, 'b'));
    trigger = trigger(strcmp({trigger.Outcome}, 'p'));
    if ~isempty(trigger)
        trigger = [trigger.RTj_realTime];
        [N_p_Dec(counterL).LickHist_N, N_p_Dec(counterL).LickHist_edges] = LickHist(trigger.', list_trialsC(n).AllLicks, [-2 2], .01, 'k', 0);
        [N_p_Dec(counterL).DeliveryHist_N, N_p_Dec(counterL).DeliveryHist_edges] = LickHist(trigger.', list_trialsC(n).JuiceTimes, [-2 2], .02, 'k', 0);
                [N_p_Dec(counterL).RunMean, N_p_Dec(counterL).RunEdges, ~] = RunSpeedHistLines(trigger.', list_trialsC(n).RunningStruct.SpeedTimesAdj, list_trialsC(n).RunningStruct.SpeedValues, -2, 2);
        N_p_Dec(counterL).trials_c = length(trigger);
        counterL = counterL + 1;
    end
end


clear N
counter = 1;
for n = 1:length(SS_IncTone)
    R = SS_IncTone(n).RecorNum;
    if [Rlist(R).PredPerc_potP] > .3
    Trials = [Rlist(R).TrialStruct];
    trigger = Trials(strcmp({Trials.TrialType}, 'b'));
    trigger = trigger(strcmp({trigger.Outcome}, 'p'));
    if ~isempty(trigger)
        trigger = [trigger.RTj_realTime];
        [N(counter,:), edges] = OneUnitHistStructTimeLimLineINDEX(trigger.', n, SS_IncTone, bmin, bmax, binsize, [0 inf], 4, 'k', NaN, 0, 0);
        % zscore if we want to switch back N(counter, :) = (N(counter, :) - mean(N(counter,1:Izero)))/std(N(counter,1:Izero));
      % N(counter, :) = (N(counter, :) - mean(N(counter,1:Izero)));
        counter = counter + 1;
    end
    end
end
N = N(~isnan(sum(N,2)),:);
if smooth(1) == 1
    N_p_Inc.N = smoothdata(N, 2, 'sgolay', smooth(2));
else
    N_p_Inc.N = N;
end
N_p_Inc.edges = edges;
counterL = 1;
trials_c = 0;
list_trialsC = Rlist([Rlist.PredPerc_potP] > .3);
for n = 1:length(list_trialsC)
    Trials = [list_trialsC(n).TrialStruct];
        trigger = Trials(strcmp({Trials.TrialType}, 'b'));
    trigger = trigger(strcmp({trigger.Outcome}, 'p'));
    if ~isempty(trigger)
        trigger = [trigger.RTj_realTime];
        [N_p_Inc(counterL).LickHist_N, N_p_Inc(counterL).LickHist_edges] = LickHist(trigger.', list_trialsC(n).AllLicks, [-2 2], .01, 'k', 0);
        [N_p_Inc(counterL).DeliveryHist_N, N_p_Inc(counterL).DeliveryHist_edges] = LickHist(trigger.', list_trialsC(n).JuiceTimes, [-2 2], .02, 'k', 0);
                [N_p_Inc(counterL).RunMean, N_p_Inc(counterL).RunEdges, ~] = RunSpeedHistLines(trigger.', list_trialsC(n).RunningStruct.SpeedTimesAdj, list_trialsC(n).RunningStruct.SpeedValues, -2, 2);
        N_p_Inc(counterL).trials_c = length(trigger);
        counterL = counterL + 1;
    end
end


clear N
counter = 1;
for n = 1:length(SS_IncTone)
    R = SS_IncTone(n).RecorNum;
    if [Rlist(R).PredPerc_potP] > .3
    Trials = [Rlist(R).LickOnsets];
    trigger = Trials(strcmp({Trials.Outcome}, 'o'));
    trigger = trigger([trigger.RewardBoo] == 0);
    if ~isempty(trigger)
        trigger = [trigger.time];
        [N(counter,:), edges] = OneUnitHistStructTimeLimLineINDEX(trigger.', n, SS_IncTone, bmin, bmax, binsize, [0 inf], 4, 'k', NaN, 0, 0);
      % N(counter, :) = (N(counter, :) - mean(N(counter,1:Izero)));
        counter = counter + 1;
    end
    end
end
if smooth(1) == 1
    N_o_Inc.N = smoothdata(N, 2, 'sgolay', smooth(2));
else
    N_o_Inc.N = N;
end
N_o_Inc.edges = edges;
counterL = 1;
trials_c = 0;
list_trialsC = Rlist([Rlist.PredPerc_potP] > .3);
for n = 1:length(list_trialsC)
   Trials = [list_trialsC(n).LickOnsets];
        trigger = Trials(strcmp({Trials.Outcome}, 'o'));
    trigger = trigger([trigger.RewardBoo] == 0);
    if ~isempty(trigger)
        trigger = [trigger.time];
        [N_o_Inc(counterL).LickHist_N, N_o_Inc(counterL).LickHist_edges] = LickHist(trigger.', list_trialsC(n).AllLicks, [-2 2], .01, 'k', 0);
        [N_o_Inc(counterL).DeliveryHist_N, N_o_Inc(counterL).DeliveryHist_edges] = LickHist(trigger.', list_trialsC(n).JuiceTimes, [-2 2], .02, 'k', 0);
        [N_o_Inc(counterL).RunMean, N_o_Inc(counterL).RunEdges, tester] = RunSpeedHistLines(trigger.', list_trialsC(n).RunningStruct.SpeedTimesAdj, list_trialsC(n).RunningStruct.SpeedValues, -2, 2);
        N_o_Inc(counterL).trials_c = length(trigger);
        counterL = counterL + 1;
    end
end

figure
hold on
shadedErrorBar2(edges(1:end-1), nanmean(cell2mat({N_o_Dec.N}.'), 1), nanstd(cell2mat({N_o_Dec.N}.'), 0, 1)/sqrt(size(cell2mat({N_o_Dec.N}.'), 1)), 'LineProp', {C(7,:)});
shadedErrorBar2(edges(1:end-1), nanmean(cell2mat({N_o_Inc.N}.'), 1), nanstd(cell2mat({N_o_Inc.N}.'), 0, 1)/sqrt(size(cell2mat({N_o_Inc.N}.'), 1)), 'LineProp', {C(8,:)});
% centerY = 1.4;
% plot([-.5; -.5], [centerY; centerY + .25], '-k',  [-.5; -.25], [centerY; centerY], '-k', 'LineWidth', 1)
xline(0, 'k', 'LineWidth', 1);
%xline(-.68, 'g', 'LineWidth', 1);
% legend({'trials with reactive licking'; 'trials with predictive licking'; }, 'Location', 'northoutside', 'FontSize', 12);
% legend('boxoff')
xlabel('time from lick onset');
ylabel('Sspk/s');
xlim([-1 .5]);
ylim([70 125]);
title('misses');
% FigureWrap(NaN, 'Trained_Sspk_missTrials_incPartTrain_Spk_Inc_DecTone_LickAligned', NaN, NaN, NaN, NaN, 2.5, 3.5);

figure
hold on
shadedErrorBar2(edges(1:end-1), nanmean(cell2mat({N_r_Dec.N}.'), 1), nanstd(cell2mat({N_r_Dec.N}.'), 0, 1)/sqrt(size(cell2mat({N_r_Dec.N}.'), 1)), 'LineProp', {C(7,:)});
shadedErrorBar2(edges(1:end-1), nanmean(cell2mat({N_r_Inc.N}.'), 1), nanstd(cell2mat({N_r_Inc.N}.'), 0, 1)/sqrt(size(cell2mat({N_r_Inc.N}.'), 1)), 'LineProp', {C(8,:)});
% centerY = 1.4;
% plot([-.5; -.5], [centerY; centerY + .25], '-k',  [-.5; -.25], [centerY; centerY], '-k', 'LineWidth', 1)
xline(0, 'k', 'LineWidth', 1);
%xline(-.68, 'g', 'LineWidth', 1);
% legend({'trials with reactive licking'; 'trials with predictive licking'; }, 'Location', 'northoutside', 'FontSize', 12);
% legend('boxoff')
xlabel('time from lick onset');
ylabel('Sspk/s');
xlim([-1 .5]);
ylim([70 125]);
title('reactions');
% FigureWrap(NaN, 'Trained_Sspk_reactTrials_incPartTrain_Spk_Inc_DecTone_LickAligned', NaN, NaN, NaN, NaN, 2.5, 3.5);

figure
hold on
shadedErrorBar2(edges(1:end-1), nanmean(cell2mat({N_p_Dec.N}.'), 1), nanstd(cell2mat({N_p_Dec.N}.'), 0, 1)/sqrt(size(cell2mat({N_p_Dec.N}.'), 1)), 'LineProp', {C(7,:)});
shadedErrorBar2(edges(1:end-1), nanmean(cell2mat({N_p_Inc.N}.'), 1), nanstd(cell2mat({N_p_Inc.N}.'), 0, 1)/sqrt(size(cell2mat({N_p_Inc.N}.'), 1)), 'LineProp', {C(8,:)});
% centerY = 1.4;
% plot([-.5; -.5], [centerY; centerY + .25], '-k',  [-.5; -.25], [centerY; centerY], '-k', 'LineWidth', 1)
xline(0, 'k', 'LineWidth', 1);
%xline(-.68, 'g', 'LineWidth', 1);
% legend({'trials with reactive licking'; 'trials with predictive licking'; }, 'Location', 'northoutside', 'FontSize', 12);
% legend('boxoff')
xlabel('time from lick onset');
ylabel('Sspk/s');
xlim([-1 .5]);
ylim([70 125]);
title('predictions');
% FigureWrap(NaN, 'Trained_Sspk_predictTrials_incPartTrain_Spk_Inc_DecTone_LickAligned', NaN, NaN, NaN, NaN, 2.5, 3.5);


