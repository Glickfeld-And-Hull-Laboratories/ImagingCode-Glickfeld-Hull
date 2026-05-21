%% Instruction sheet
instructions.sess_list=[85 89 92]; %the sessions to concatenate and analyze
instructions.ds='DART_V1_YM90K_Celine'; %name of the datasheet file
instructions.experimentFolder = 'VIP_YM90K'; %name of the experiment folder
instructions.refDay='2'; %which of the two days of data was used as the reference for matching? 
instructions.load_retino = true; %whether to look for and load retinotopic alignment data
instructions.targetCon = [.25 .5 1]; %the contrasts to analyze, in case not all datasets have the same set of contrasts
instructions.targetSize = [7.5 15.0 30.0000 60.0 120.000]; %the sizes to analyze, in case not all datasets have the same set of sizes
instructions.sizeFilter = [] %set this to select for cells responsive to a partiular size
instructions.frame_rate = 15;

instructions.runOnset_sess_list = []; %the sessions to concatenate for running onset analyses - may be different from above if some sessions lack running trials
instructions.stillTimeList = [5]; % seconds of stillness required before run onset
%% Notes
%Use thuis space for a description of the dataset