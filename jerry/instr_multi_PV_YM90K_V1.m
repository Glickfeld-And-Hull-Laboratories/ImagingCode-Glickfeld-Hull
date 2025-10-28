%% Instruction sheet
instructions.sess_list=[42 44]; %the sessions to concatenate and analyze
instructions.ds='DART_expt_info_jerry';
instructions.experimentFolder = 'PV_YM90K';
instructions.refDay='2';
instructions.targetCon = [.125 .25 .5 1]; %the contrasts to analyze, in case not all datasets have the same set of contrasts
instructions.targetSize = [20.0000  1000.0000]; %the sizes to analyze, in case not all datasets have the same set of sizes
instructions.frame_rate = 15;
%% Notes
% these are PV YM90K-DART experiments i3314 and i3315 (the two good ones)