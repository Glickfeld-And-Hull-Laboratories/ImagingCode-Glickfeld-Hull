bmin = -1;
bmax = 1;
binsize = .01;
Izero = abs(bmin/binsize) - 1;

% create SS & SS paired- only need to do once
% counter = 1;
% for n = 1:length(SS)
%     if ~isempty(SS(n).PCpair)
%         struct = SS([SS.RecorNum] == [SS(n).RecorNum]);
%         if ~isempty(find([struct.unitID] == SS(n).PCpair, 1))
%             SS(counter) = SS(n);
%             SS_paired(counter) = SS(find([struct.unitID] == SS(n).PCpair, 1));
%             if length(find([struct.unitID] == SS(n).PCpair, 1)) > 1
%                 counter = counter + 1;
%                 SS(counter) = SS(n);
%                 SS_paired(counter) = SS(find([struct.unitID] == SS(n).PCpair, 2));
%             end
%             counter = counter + 1;
%         end
%     end
% end

% Remove any old calculations
FieldNames = fields(SS);
if any(strcmp(FieldNames, 'JuiceTimes_both'))
    SS = rmfield(SS, 'JuiceTimes_both');
end
if any(strcmp(FieldNames, 'JuiceTimes_clk'))
    SS = rmfield(CSSS, 'JuiceTimes_clk');
end
if any(strcmp(FieldNames, 'JuiceTimes_clkMod'))
    SS = rmfield(SS, 'JuiceTimes_clkMod');
end
if any(strcmp(FieldNames, 'ToneTimes'))
    SS = rmfield(SS, 'ToneTimes');
end
if any(strcmp(FieldNames, 'ToneMod'))
    SS = rmfield(SS, 'ToneMod');
end


SD = 4;
% for n = 1:length(SS)
%     R = SS(n).RecorNum;
%     Trials = [Rlist(R).TrialStruct];
%     if [Rlist(R).day] == 7 | [Rlist(R).day] == 8 | [Rlist(R).day] == 9
%         trigger = Trials(strcmp({Trials.TrialType}, 'b'));
%     else
%         trigger = Trials(strcmp({Trials.TrialType}, 'j'));
%     end
%     if length(trigger) > 20
%         trigger = [trigger.JuiceTime];
%         [N, edges] = OneUnitHistStructTimeLimLineINDEX(trigger.', n, SS, bmin, bmax, binsize, [0 inf], 4, 'k', NaN, 0, 0);
%         [struct.LatLow, struct.LatHigh, struct.modBoo, struct.Dir, struct.doubleBoo] = LatencyMod(0, N, edges, SD, [0 .4], 0);
%         SS(n).JuiceTimes_clk.modLatStruct = struct;
%         SS(n).JuiceTimes_clk.Dir = struct.Dir;
%         if SS(n).JuiceTimes_clk.modLatStruct.modBoo == 1
%             SS(n).JuiceTimes_clkMod = true;
%         else
%             if [Rlist(R).day] >= 7 %here I force it to be only false if Cspk is unresponsive to both juice alone and juice after cue
%                 trigger = Trials(strcmp({Trials.TrialType}, 'b'));
%                 trigger = [trigger.JuiceTime];
%                 [N, edges] = OneUnitHistStructTimeLimLineINDEX(trigger.', n, SS, bmin, bmax, binsize, [0 inf], 4, 'k', NaN, 0, 0);
%                 [struct.LatLow, struct.LatHigh, struct.modBoo, struct.Dir, struct.doubleBoo] = LatencyMod(0, N, edges, SD, [0 .4], 0);
%                 SS(n).JuiceTimes_both.modLatStruct = struct;
%                 SS(n).JuiceTimes_both.Dir = struct.Dir;
%                 if SS(n).JuiceTimes_both.modLatStruct.modBoo == 0
%                     SS(n).JuiceTimes_clkMod = false;
%                 else 
%                    SS(n).JuiceTimes_clkMod = true; 
%                 end
%             else
%              SS(n).JuiceTimes_clkMod = false;
%             end
%         end
%     else
%         SS(n).JuiceTimes_clkMod = NaN;
%     end
% end

for n = 1:length(SS)
    if ~isempty([Rlist(SS(n).RecorNum).ToneTimes])
        [N, edges] = OneUnitHistStructTimeLimLineINDEX(Rlist(SS(n).RecorNum).ToneTimes, n, SS, bmin, bmax, binsize, [0 inf], 4, 'k', NaN, 0, 0);
        [struct.LatLow, struct.LatHigh, struct.modBoo, struct.Dir, struct.doubleBoo] = LatencyMod(0, N, edges, SD, [0 .4], 0);
        SS(n).ToneTimes.modLatStruct = struct;
        SS(n).ToneTimes.Dir = struct.Dir;
        if SS(n).ToneTimes.modLatStruct.modBoo == 1
            SS(n).ToneMod = true;
            SS(n).ToneModDir = struct.Dir;
        else
            SS(n).ToneMod = false;
        end
    else
        SS(n).ToneMod = NaN;
    end
end

%
% SSpaired_train = SS([SS.TrainBoo]);
% SSpaired_naive = SS([SS.day] == 7 | [SS.day] == 8 | [SS.day] == 9);
% %
% % SSpaired_mod = SS([SS.LickOnsetEmptyMod] | [SS.ToneMod] | [SS.JuiceTimes_clkMod]);
% % SSpaired_mod = SS_paired([SS.LickOnsetEmptyMod] | [SS.ToneMod] | [SS.JuiceTimes_clkMod]);


%SSpaired_mod = SS([SS.ToneMod] | [SS.JuiceTimes_clkMod] | [SS.NoJuiceClkMod]);
%SSpaired_mod = SS_paired([SS.ToneMod] | [SS.JuiceTimes_clkMod] | [SS.NoJuiceClkMod]);

% SS_ModJuiceClk =  SS([SS.JuiceTimes_clkMod] == 1);
% %SS_paired_modJuiceClk =  SS_paired([SS.JuiceTimes_clkMod] == 1);
% 
% SS_NOmodJuiceClk =  SS([SS.JuiceTimes_clkMod] == 0);
% %SS_paired_NOmodJuiceClk =  SS_paired([SS.JuiceTimes_clkMod] == 0);


SS_ModTone = SS([SS.ToneMod] == 1);
SS_IncTone = SS_ModTone([SS_ModTone.ToneModDir] == 1);
SS_DecTone = SS_ModTone([SS_ModTone.ToneModDir] == -1);
SS_NOtoneMod = SS([SS.ToneMod] == 0);
%
% SSpairedMod_train = SSpaired_mod([SSpaired_mod.TrainBoo]);
% SSpairedMod_naive = SS([SS.day] == 7 | [SS.day] == 8 | [SS.day] == 9);
% SSpairedMod_train = SSpaired_mod([SSpaired_mod.TrainBoo]);
% SSpairedMod_naive = SSpaired_mod([SSpaired_mod.day]  == 7 | [SSpaired_mod.day] == 8 | [SSpaired_mod.day] == 9);
%

% SSpaired_NOmod = SS(~[SS.ToneMod] & ~[SS.JuiceTimes_clkMod] & ~ [SS.NoJuiceClkMod]);
% SSpaired_NOmod = SS_paired(~[SS.ToneMod] & ~[SS.JuiceTimes_clkMod] & ~ [SS.NoJuiceClkMod]);

% SSpairedNOMod_train = SSpaired_NOmod([SSpaired_NOmod.TrainBoo]);
% SSpairedNOMod_naive = SSpaired_NOmod([SSpaired_NOmod.day] <= 7);
% SSpairedNOMod_train = SSpaired_NOmod([SSpaired_NOmod.TrainBoo]);
% SSpairedNOMod_naive = SSpaired_NOmod([SSpaired_NOmod.day] <= 7);
