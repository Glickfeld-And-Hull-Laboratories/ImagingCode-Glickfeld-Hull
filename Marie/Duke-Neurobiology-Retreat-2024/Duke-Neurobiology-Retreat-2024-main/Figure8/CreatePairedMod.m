bmin = -1;
bmax = 1;
binsize = .01;
bsize = .01;
Izero = abs(bmin/bsize) - 1;

% use new CS_paired and SS_paired creation code from CreateDataset
% clear CS_paired
% clear SS_paired
% % create CS_paired & SS paired- only need to do once
% counter = 1;
% for n = 1:length(CS)
%     if ~isempty(CS(n).PCpair)
%         struct = SS([SS.RecorNum] == [CS(n).RecorNum]);
%         if ~isempty(find([struct.unitID] == CS(n).PCpair, 1))
%             CS_paired(counter) = CS(n);
%             SS_paired(counter) = struct(find([struct.unitID] == CS(n).PCpair, 1));
%             if length(find([struct.unitID] == CS(n).PCpair, 1)) > 1
%                 counter = counter + 1;
%                 CS_paired(counter) = CS(n);
%                 SS_paired(counter) = struct(find([struct.unitID] == CS(n).PCpair, 2));
%             end
%             counter = counter + 1;
%         end
%     end
% end
clear struct

CS_paired = rmfield(CS_paired, {'JuiceTimes_both', 'JuiceTimes_clk', 'JuiceTimes_clkMod', 'ToneTimes', 'ToneMod'});

SD = 4.5;
for n = 1:length(CS_paired)
    R = CS_paired(n).RecorNum;
    Trials = [Rlist(R).TrialStruct];
    if [Rlist(R).day] == 7 | [Rlist(R).day] == 8 | [Rlist(R).day] == 9
        trigger = Trials(strcmp({Trials.TrialType}, 'b'));
    else
        trigger = Trials(strcmp({Trials.TrialType}, 'j'));
    end
    if length(trigger) > 20
        trigger = [trigger.JuiceTime];
        [N, edges] = OneUnitHistStructTimeLimLineINDEX(trigger.', n, CS_paired, bmin, bmax, bsize, [0 inf], 4, 'k', NaN, 0, 0);
        [struct.LatLow, struct.LatHigh, struct.modBoo, struct.Dir, struct.doubleBoo] = LatencyMod(0, N, edges, SD, [0 .4], 0);
        CS_paired(n).JuiceTimes_clk.modLatStruct = struct;
        CS_paired(n).JuiceTimes_clk.Dir = struct.Dir;
        if CS_paired(n).JuiceTimes_clk.modLatStruct.modBoo == 1
            CS_paired(n).JuiceTimes_clkMod = true;
        else
            if [Rlist(R).day] >= 7 %here I force it to be only false if Cspk is unresponsive to both juice alone and juice after cue
                trigger = Trials(strcmp({Trials.TrialType}, 'b'));
                trigger = [trigger.JuiceTime];
                [N, edges] = OneUnitHistStructTimeLimLineINDEX(trigger.', n, CS_paired, bmin, bmax, bsize, [0 inf], 4, 'k', NaN, 0, 0);
                [struct.LatLow, struct.LatHigh, struct.modBoo, struct.Dir, struct.doubleBoo] = LatencyMod(0, N, edges, SD, [0 .4], 0);
                CS_paired(n).JuiceTimes_both.modLatStruct = struct;
                CS_paired(n).JuiceTimes_both.Dir = struct.Dir;
                if CS_paired(n).JuiceTimes_both.modLatStruct.modBoo == 0
                    CS_paired(n).JuiceTimes_clkMod = false;
                else 
                   CS_paired(n).JuiceTimes_clkMod = true; 
                end
            else
             CS_paired(n).JuiceTimes_clkMod = false;
            end
        end
    else
        CS_paired(n).JuiceTimes_clkMod = NaN;
    end
end

for n = 1:length(CS_paired)
    if ~isempty([Rlist(CS_paired(n).RecorNum).ToneTimes])
        [N, edges] = OneUnitHistStructTimeLimLineINDEX(Rlist(CS_paired(n).RecorNum).ToneTimes, n, CS_paired, bmin, bmax, bsize, [0 inf], 4, 'k', NaN, 0, 0);
        [struct.LatLow, struct.LatHigh, struct.modBoo, struct.Dir, struct.doubleBoo] = LatencyMod(0, N, edges, SD, [0 .4], 0);
        CS_paired(n).ToneTimes.modLatStruct = struct;
        CS_paired(n).ToneTimes.Dir = struct.Dir;
        if CS_paired(n).ToneTimes.modLatStruct.modBoo == 1
            CS_paired(n).ToneMod = true;
        else
            CS_paired(n).ToneMod = false;
        end
    else
        CS_paired(n).ToneMod = NaN;
    end
end


%CS_pairedpaired_mod = CS_paired([CS_paired.ToneMod] | [CS_paired.JuiceTimes_clkMod] | [CS_paired.NoJuiceClkMod]);
%SSpaired_mod = SS_paired([CS_paired.ToneMod] | [CS_paired.JuiceTimes_clkMod] | [CS_paired.NoJuiceClkMod]);

CS_paired_ModJuiceClk =  CS_paired([CS_paired.JuiceTimes_clkMod] == 1);
SS_paired_ModJuiceClk =  SS_paired([CS_paired.JuiceTimes_clkMod] == 1);

CS_paired_NOmodJuiceClk =  CS_paired([CS_paired.JuiceTimes_clkMod] == 0);
SS_paired_NOmodJuiceClk =  SS_paired([CS_paired.JuiceTimes_clkMod] == 0);


CS_paired_ModTone = CS_paired([CS_paired.ToneMod] == 1);
SS_paired_ModTone = SS_paired([CS_paired.ToneMod] == 1);

CS_paired_NOtoneMod = CS_paired([CS_paired.ToneMod] == 0);
SS_paired_NOtoneMod = SS_paired([CS_paired.ToneMod] == 0);
%
% CS_pairedpairedMod_train = CS_pairedpaired_mod([CS_pairedpaired_mod.TrainBoo]);
% CS_pairedpairedMod_naive = CS_paired([CS_paired.day] == 7 | [CS_paired.day] == 8 | [CS_paired.day] == 9);
% SSpairedMod_train = SSpaired_mod([SSpaired_mod.TrainBoo]);
% SSpairedMod_naive = SSpaired_mod([SSpaired_mod.day]  == 7 | [SSpaired_mod.day] == 8 | [SSpaired_mod.day] == 9);
%

% CS_pairedpaired_NOmod = CS_paired(~[CS_paired.ToneMod] & ~[CS_paired.JuiceTimes_clkMod] & ~ [CS_paired.NoJuiceClkMod]);
% SSpaired_NOmod = SS_paired(~[CS_paired.ToneMod] & ~[CS_paired.JuiceTimes_clkMod] & ~ [CS_paired.NoJuiceClkMod]);

% CS_pairedpairedNOMod_train = CS_pairedpaired_NOmod([CS_pairedpaired_NOmod.TrainBoo]);
% CS_pairedpairedNOMod_naive = CS_pairedpaired_NOmod([CS_pairedpaired_NOmod.day] <= 7);
% SSpairedNOMod_train = SSpaired_NOmod([SSpaired_NOmod.TrainBoo]);
% SSpairedNOMod_naive = SSpaired_NOmod([SSpaired_NOmod.day] <= 7);
