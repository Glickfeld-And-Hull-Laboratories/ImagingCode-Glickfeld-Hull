%% Instruction sheet
instructions.sess_list=[8 10 20 22]; %the sessions to concatenate and analyze
instructions.ds='DART_V1_YM90K_Celine';
instructions.experimentFolder = 'SST_YM90K';
instructions.refDay='2';
instructions.targetCon = [.25 .5 1]; %the contrasts to analyze, in case not all datasets have the same set of contrasts
instructions.targetSize = [7.5000   30.0000  120.0000]; %the sizes to analyze, in case not all datasets have the same set of sizes
instructions.frame_rate = 15;
%% Notes
% This is a set of SST mice for which I collected data in area LM before
% and after YM90K-DART application.