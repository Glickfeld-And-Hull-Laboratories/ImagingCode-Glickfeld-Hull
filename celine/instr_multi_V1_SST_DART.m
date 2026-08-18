%% Instruction sheet
instructions.sess_list=[138 142 163 171 178 190 294 307 323 333]; %the sessions to concatenate and analyze
instructions.ds='DART_V1_contrast_ori_Celine';
instructions.experimentFolder = 'SST_YM90K';
instructions.refDay='2';
instructions.stillTimeList = [5]; % seconds of stillness required before run onset
instructions.targetCon = [.25 .5 1]; %the contrasts to analyze, in case not all datasets have the same set of contrasts
instructions.targetSize = []; %the sizes to analyze, in case not all datasets have the same set of sizes
instructions.frame_rate = 15;
%% Notes
% This is a set of SST mice for which I collected data in area V1 before
% and after YM90K-DART application. This is the data set used in Cammarata
% et al 2025

%full session list 138 142 163 171 178 190 294 307 323 333