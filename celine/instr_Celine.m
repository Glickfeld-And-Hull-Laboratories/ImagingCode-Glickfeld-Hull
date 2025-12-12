%% Instruction sheet
instructions.sess_list=[55]; %the sessions to concatenate and analyze
instructions.runOnset_sess_list = []; %skipping sessions that didn't have any running onsets
instructions.ds='DART_V1_YM90K_Celine';
instructions.experimentFolder = 'SST_YM90K';
instructions.refDay='2';
instructions.stillTimeList = [5]; % seconds of stillness required before run onset
instructions.targetCon = [.25 .5 1]; %the contrasts to analyze, in case not all datasets have the same set of contrasts
instructions.targetSize = [7.5 15.0 30.0000 60.0 120.000]; %the sizes to analyze, in case not all datasets have the same set of sizes
instructions.frame_rate = 15;
%% Notes
%cooh=52 
% alx-DART = 47 44