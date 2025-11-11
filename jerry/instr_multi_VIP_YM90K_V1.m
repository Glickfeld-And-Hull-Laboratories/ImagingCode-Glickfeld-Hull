%% Instruction sheet
instructions.session= [39]; %the sessions to concatenate and analyze
instructions.ds='DART_V1_YM90K_Celine';
instructions.experimentFolder = 'VIP_YM90K';
instructions.refDay='2';
instructions.targetCon = [.125 .25 .5 1]; %the contrasts to analyze, in case not all datasets have the same set of contrasts
instructions.targetSize = [20.0000  1000.0000]; %the sizes to analyze, in case not all datasets have the same set of sizes
instructions.frame_rate = 15;
%% Notes
% this is for VIP PEG experiment i3333 (4mM YM-PEG, 4mM Blank-DART, 0.4mM
% Alx-DART)