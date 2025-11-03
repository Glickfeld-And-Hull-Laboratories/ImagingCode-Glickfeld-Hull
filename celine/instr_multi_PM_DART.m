%% Instruction sheet
instructions.sess_list=[37]; %the sessions to concatenate and analyze
instructions.runOnset_sess_list = [37]; %skipping sessions that didn't have any running onsets
instructions.ds='DART_V1_YM90K_Celine';
instructions.experimentFolder = 'SST_YM90K';
instructions.refDay='2';
instructions.stillTimeList = [5]; % seconds of stillness required before run onset
instructions.targetCon = [.25 .5 1]; %the contrasts to analyze, in case not all datasets have the same set of contrasts
instructions.targetSize = [15.000   30.0000  60.0000]; %the sizes to analyze, in case not all datasets have the same set of sizes
instructions.frame_rate = 15;
%% Notes
% This is a set of SST mice for which I collected data in area PM before
% and after YM90K-DART application.

%note that session 37, mouse i2198, was epoched using the counter vals 