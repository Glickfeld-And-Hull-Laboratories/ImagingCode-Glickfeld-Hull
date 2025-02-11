%Determine which CSpk are modulated by various things and compute the
%scatter plots in Figure 4 from Neurobiology Retreat 2024.

% Set parameters
bmin = -1;
bmax = 1;
binsize = .01;
Izero = abs(bmin/binsize) - 1;
SD = 3;

% Remove any old calculations
FieldNames = fields(CS);
if any(strcmp(FieldNames, 'JuiceTimes_both'))
    CS = rmfield(CS, 'JuiceTimes_both');
end
if any(strcmp(FieldNames, 'JuiceTimes_clk'))
    CS = rmfield(CS, 'JuiceTimes_clk');
end
if any(strcmp(FieldNames, 'JuiceTimes_clkMod'))
    CS = rmfield(CS, 'JuiceTimes_clkMod');
end
if any(strcmp(FieldNames, 'ToneTimes'))
    CS = rmfield(CS, 'ToneTimes');
end
if any(strcmp(FieldNames, 'ToneMod'))
    CS = rmfield(CS, 'ToneMod');
end


for n = 1:length(CS)
    R = CS(n).RecorNum;
    Trials = [Rlist(R).TrialStruct];
    if [Rlist(R).day] == 7 | [Rlist(R).day] == 8 | [Rlist(R).day] == 9
        trigger = Trials(strcmp({Trials.TrialType}, 'b'));
    else
        trigger = Trials(strcmp({Trials.TrialType}, 'j'));
    end
    if length(trigger) > 20
        trigger = [trigger.JuiceTime];
        [N, edges] = OneUnitHistStructTimeLimLineINDEX(trigger.', n, CS, bmin, bmax, binsize, [0 inf], 4, 'k', NaN, 0, 0);
        [struct.LatLow, struct.LatHigh, struct.modBoo, struct.Dir, struct.doubleBoo] = LatencyMod(0, N, edges, SD, [0 .4], 0);
        CS(n).JuiceTimes_clk.modLatStruct = struct;
        CS(n).JuiceTimes_clk.Dir = struct.Dir;
        if CS(n).JuiceTimes_clk.modLatStruct.modBoo == 1
            CS(n).JuiceTimes_clkMod = true;
        else
            if [Rlist(R).day] >= 7 %here I force it to be only false if Cspk is unresponsive to both juice alone and juice after cue
                trigger = Trials(strcmp({Trials.TrialType}, 'b'));
                trigger = [trigger.JuiceTime];
                [N, edges] = OneUnitHistStructTimeLimLineINDEX(trigger.', n, CS, bmin, bmax, binsize, [0 inf], 4, 'k', NaN, 0, 0);
                [struct.LatLow, struct.LatHigh, struct.modBoo, struct.Dir, struct.doubleBoo] = LatencyMod(0, N, edges, SD, [0 .4], 0);
                CS(n).JuiceTimes_both.modLatStruct = struct;
                CS(n).JuiceTimes_both.Dir = struct.Dir;
                if CS(n).JuiceTimes_both.modLatStruct.modBoo == 0
                    CS(n).JuiceTimes_clkMod = false;
                else
                    CS(n).JuiceTimes_clkMod = true;
                end
            else
                CS(n).JuiceTimes_clkMod = false;
            end
        end
    else
        CS(n).JuiceTimes_clkMod = NaN;
    end
end

for n = 1:length(CS)
    if ~isempty([Rlist(CS(n).RecorNum).ToneTimes])
        [N, edges] = OneUnitHistStructTimeLimLineINDEX(Rlist(CS(n).RecorNum).ToneTimes, n, CS, bmin, bmax, binsize, [0 inf], 4, 'k', NaN, 0, 0);
        [struct.LatLow, struct.LatHigh, struct.modBoo, struct.Dir, struct.doubleBoo] = LatencyMod(0, N, edges, SD, [0 .4], 0);
        CS(n).ToneTimes.modLatStruct = struct;
        CS(n).ToneTimes.Dir = struct.Dir;
        if CS(n).ToneTimes.modLatStruct.modBoo == 1
            CS(n).ToneMod = true;
        else
            CS(n).ToneMod = false;
        end
    else
        CS(n).ToneMod = NaN;
    end
end


CS_ModJuiceClk =  CS([CS.JuiceTimes_clkMod] == 1);
CS_NOmodJuiceClk =  CS([CS.JuiceTimes_clkMod] == 0);


CS_ModTone = CS([CS.ToneMod] == 1);
CS_NOtoneMod = CS([CS.ToneMod] == 0);

% Calculate percentage of Cspk modulated by Tone on recordings with more
% than 25 Cspks.
%
for n = 1:length(Rlist)
    CS_thisR = CS([CS.RecorNum] == n);
    if length(CS_thisR) > 25
        Rlist(n).ToneModCS_perc = length(CS_thisR([CS_thisR.ToneMod] == 1))/length(CS_thisR);
    else
        Rlist(n).ToneModCS_perc = NaN;
    end
    if isempty(Rlist(n).ToneTimes)
        Rlist(n).ToneModCS_perc = NaN;
    end
end

figure
hold on
%scatter([Rlist([Rlist.day] <=3).ToneModCS_perc]*100, [Rlist([Rlist.day] <= 3).PredPerc_potP]*100, 'MarkerFaceColor', C(4,:), 'MarkerEdgeColor', 'none')
xlabel('% tone modulated CS');
ylabel('% predicitive licking');
scatter([Rlist([Rlist.day] >= 7).ToneModCS_perc]*100, [Rlist([Rlist.day] >= 7).PredPerc_potP]*100, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none')
h = lsline;
set(h, 'LineWidth', 1.5, 'Color', [.5 .5 .5])
hold on
scatter([Rlist([Rlist.day] == 7 | [Rlist.day] == 8 | [Rlist.day] == 9 ).ToneModCS_perc]*100, [Rlist([Rlist.day] == 7 | [Rlist.day] == 8 | [Rlist.day] == 9 ).PredPerc_potP]*100, 'MarkerFaceColor', C(1,:), 'MarkerEdgeColor', 'none')
%scatter([Rlist([Rlist.day] <=3).ToneModCS_perc]*100, [Rlist([Rlist.day] <= 3).PredPerc_potP]*100, 'MarkerFaceColor', C(4,:), 'MarkerEdgeColor', 'none')
scatter([Rlist([Rlist.day] >= 10).ToneModCS_perc]*100, [Rlist([Rlist.day] >= 10).PredPerc_potP]*100, 'MarkerFaceColor', C(2,:), 'MarkerEdgeColor', 'none')
title('Prediction vs Cspk Modulation');
xlim([20 90]);
[rho, pVal] = corr([Rlist([Rlist.day] >= 7).ToneModCS_perc].', [Rlist([Rlist.day] >= 7).PredPerc_potP].', 'Rows', 'complete');
text(80,100,['rho = ' num2str(rho)],'Color','k','FontSize',8)
text(80, 95,['p = ' num2str(pVal)],'Color','k','FontSize',8)
% FigureWrap(NaN, 'PercCspkMod_PercPred', NaN, NaN, NaN, NaN, 2.5, 3.5);


figure
hold on
scatter([Rlist([Rlist.day] >= 7).day]-6, [Rlist([Rlist.day] >= 7).ToneModCS_perc]*100, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none')
h = lsline;
set(h, 'LineWidth', 1.5, 'Color', [.5 .5 .5])
scatter([Rlist([Rlist.day] == 7 | [Rlist.day] == 8 | [Rlist.day] == 9 ).day]-6, [Rlist([Rlist.day] == 7 | [Rlist.day] == 8 | [Rlist.day] == 9 ).ToneModCS_perc]*100, 'MarkerFaceColor', C(1,:), 'MarkerEdgeColor', 'none')
scatter([Rlist([Rlist.day] >= 10).day] -6, [Rlist([Rlist.day] >= 10).ToneModCS_perc]*100, 'MarkerFaceColor', C(2,:), 'MarkerEdgeColor', 'none')
xlabel('day of training')
ylabel('% Cspk modulated by tone')
title('Cspk Modulation vs. Day')
[rho, pVal] = corr([Rlist([Rlist.day] >= 7).day].'-6, [Rlist([Rlist.day] >= 7).ToneModCS_perc].'*100, 'Rows', 'complete')
text(17,85,['rho = ' num2str(rho)],'Color','k','FontSize',8)
text(17, 80,['p = ' num2str(pVal)],'Color','k','FontSize',8)
% FigureWrap(NaN, 'day_PercCspkMod', NaN, NaN, NaN, NaN, 2.5, 3.5);

figure
hold on
scatter([Rlist([Rlist.day] >= 7).day]-6, [Rlist([Rlist.day] >= 7).PredPerc_potP]*100, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none')
h = lsline;
set(h, 'LineWidth', 1.5, 'Color', [.5 .5 .5])
scatter([Rlist([Rlist.day] == 7 | [Rlist.day] == 8 | [Rlist.day] == 9 ).day]-6, [Rlist([Rlist.day] == 7 | [Rlist.day] == 8 | [Rlist.day] == 9 ).PredPerc_potP]*100, 'MarkerFaceColor', C(1,:), 'MarkerEdgeColor', 'none')
scatter([Rlist([Rlist.day] >= 10).day] -6, [Rlist([Rlist.day] >= 10).PredPerc_potP]*100, 'MarkerFaceColor', C(2,:), 'MarkerEdgeColor', 'none')
xlabel('day of training')
ylabel('% trials with predictive licking')
title('Prediction vs Day')
[rho, pVal] = corr([Rlist([Rlist.day] >= 7).day].'-6, [Rlist([Rlist.day] >= 7).PredPerc_potP].'*100, 'Rows', 'complete')
text(16, 20,['rho = ' num2str(rho)],'Color','k','FontSize',8)
text(16, 15,['p = ' num2str(pVal)],'Color','k','FontSize',8)
% FigureWrap(NaN, 'day_PercPredTrials', NaN, NaN, NaN, NaN, 2.5, 3.5);

%
% figure
% hold on
% scatter([Rlist([Rlist.day] >= 7).ToneModCS_perc]*100, [Rlist([Rlist.day] >= 7).meanRT], 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none')
% h = lsline;
% set(h, 'LineWidth', 1.5, 'Color', [.5 .5 .5])
% scatter([Rlist([Rlist.day] == 7 | [Rlist.day] == 8 | [Rlist.day] == 9 ).ToneModCS_perc]*100, [Rlist([Rlist.day] == 7 | [Rlist.day] == 8 | [Rlist.day] == 9).meanRT], 'MarkerFaceColor', C(1,:), 'MarkerEdgeColor', 'none')
% scatter([Rlist([Rlist.day] >= 10).ToneModCS_perc]*100, [Rlist([Rlist.day] >= 10).meanRT], 'MarkerFaceColor', C(2,:), 'MarkerEdgeColor', 'none')
% xlabel('Percent Cspk modulated by tone')
% ylabel('meanRT')
% % FigureWrap(NaN, 'PercCspkMod_meanRT', NaN, NaN, NaN, NaN, NaN, NaN);

