%Made by anthony so may be buggy
%Created:2/15/2023
%should be able to run everything up until Tprime generator, run those and
%then run the rest 

%open JuiceTimes and initialize string names

RecordingLabel = "2021_03_23_1615_g0";
FilePathnoRec = "D:\Anthony\CellSorting\2021_03_23_1615_g0\";
FilePath = FilePathnoRec + RecordingLabel;
Channel = "8";
fid = fopen(RecordingLabel + "_tcat.nidq.XD_" + ChannelNumber + "_1_0.txt");
JuiceTimes = fscanf(fid, '%f');
fclose(fid);
%if no tone times --> 
ToneTimes = JuiceTimes - 0.682;

%if tone times --> 
%fid = fopen('1673_230213_g0_tcat.nidq.XD_2_3_0.txt');
%ToneTimes = fscanf(fid, '%f');
%fclose(fid);

%remember to change the channel --> 2 if on 2, 5 if its on 8
%Look at Licks and make sure everything looks good, then generate trial
%struct) 
FindLickLevels(JuiceTimes, 3, 5, 20, .8);
[TrialStruct, JuiceAlone, ToneAlone, JuiceAfterTone, ToneBeforeJuice, FictiveJuice] = JuiceToneCreateTrialStAT(JuiceTimes, ToneTimes);
%Find All The Licks 
[AllLicks, LickDetectParams] = FindAllLicks(TrialStruct(1).FictiveJuice - 30, TrialStruct(end).FictiveJuice, .6, .8);
save AllLicks.txt AllLicks -ascii -double


%STOP AND DO TPRIME SYNC

JuiceTimesAdjtprimestring = "TPrime -syncperiod=1.000000 -tostream=" + FilePath + "_tcat.imec0.ap.SY_384_6_500.txt -fromstream=1," + FilePath + "_tcat.nidq.XD_" + Channel + "_0_0.txt -events=1," + FilePath + "_tcat.nidq.XD_" + Channel + "_1_0.txt," + FilePathnoRec + "JuiceTimesAdj.txt" ;
AllLicksAdjtprimestring =  "TPrime -syncperiod=1.000000 -tostream=" + FilePath + "_tcat.imec0.ap.SY_384_6_500.txt -fromstream=1," + FilePath + "_tcat.nidq.XD_" + Channel + "_0_0.txt -events=1," + FilePathnoRec + "AllLicks.txt," + FilePathnoRec + "AllLicksAdj.txt" ;


%RETURN AFTER TPRIME SYNC

%open adjusted file 
fid = fopen('JuiceTimesAdj.txt');
JuiceTimesAdj = fscanf(fid, '%f');
fclose(fid);
fid = fopen('AllLicksAdj.txt');
AllLicksAdj = fscanf(fid, '%f');
fclose(fid);
ToneTimesAdj = JuiceTimesAdj - 0.682;

%fid = fopen('ToneTimesAdj.txt');
%ToneTimesAdj = fscanf(fid, '%f');
%fclose(fid);
%Make our final trial structs

[TrialStructAdj, JuiceAloneAdj, ToneAloneAdj, JuiceAfterToneAdj, ToneBeforeJuiceAdj, FictiveJuiceAdj] = JuiceToneCreateTrialStAT(JuiceTimesAdj, ToneTimesAdj);

%Summary Structs
load MEH_chanMap.mat
SummaryStruct_SS = SummaryStructMakerClassical([ComplexStruct1.PCUnitID], 'SS_pause', 'expert', [], 'reward_ACWF_SS', GoodUnitStruct, 1, [0 inf],  MEH_chanMap, TrialStructRTtAdj, [ComplexStruct1.ComplexUnitID], AllLicksAdj);
SummaryStruct_MLI = SummaryStructMakerClassical([MLIstruct1.MLIUnitID].', 'MLI', 'expert', 'ccg', 'reward_ACWF_MLI', GoodUnitStruct, 1, [0 inf],  MEH_chanMap, TrialStructRTtAdj, [MLIstruct1.PCUnitID], AllLicksAdj);
SummaryStruct_CS = SummaryStructMakerClassical([ComplexStruct1.ComplexUnitID], 'CS_pause', 'expert', [], 'reward_ACWF_CS', GoodUnitStruct, 1, [0 inf],  MEH_chanMap, TrialStructRTtAdj, [ComplexStruct1.PCUnitID], AllLicksAdj);

