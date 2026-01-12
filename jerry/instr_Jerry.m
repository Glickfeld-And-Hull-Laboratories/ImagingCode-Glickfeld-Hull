%% Instruction sheet
instructions.session='2';
instructions.ds='DREADD_V1_datasheet';
instructions.tIdxSource = ["PD" "PD"]; % PD or MW or cS; ref day then match day
instructions.refDay='2'; % 1st day or 2nd day used as reference
instructions.tDropBool = false; % true or false
instructions.tDropRefDay = []; % ref day is always instructions.session
instructions.tDropMatchDay = [];
instructions.tDropAction = 'DEL'; % DEL or REP