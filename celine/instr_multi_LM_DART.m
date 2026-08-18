%% Instruction sheet
instructions.sess_list=[8 10 20 22 35]; %the sessions to concatenate and analyze
instructions.runOnset_sess_list = [10 22 35]; %skipping sessions that didn't have any running onsets
instructions.ds='DART_V1_YM90K_Celine';
instructions.experimentFolder = 'SST_YM90K';
instructions.refDay='2';
instructions.stillTimeList = [5]; % seconds of stillness required before run onset
instructions.targetCon = [.25 .5 1]; %the contrasts to analyze, in case not all datasets have the same set of contrasts
instructions.targetSize = [7.5000   30.0000  120.0000]; %the sizes to analyze, in case not all datasets have the same set of sizes
instructions.frame_rate = 15;
%% Notes
% This is a set of SST mice for which I collected data in area LM before
% and after YM90K-DART application.

%full session list 8 10 20 22 35