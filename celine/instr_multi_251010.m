%% Instruction sheet for run onset concatenation
instructions.sess_list = [66, 57, 53, 44]; % the sessions to concatenate and analyze
instructions.ds = 'DART_V1_atropine_Celine';
instructions.experimentFolder = 'SST_atropine';
instructions.refDay = '2'; % '1' for baseline as reference, '2' for post-DART as reference
instructions.stillTimeList = [5]; % seconds of stillness required before run onset
instructions.frame_rate = 15; % imaging frame rate in Hz

%% Notes
% This instruction file is used for concatenating run onset data across
% multiple sessions. Make sure that runOnsetDataCollect_optionalITI has
% been run for each session in sess_list before running runOnsetConcat.