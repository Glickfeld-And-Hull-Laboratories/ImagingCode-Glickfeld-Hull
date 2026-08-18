%% Instruction sheet
instructions.session='92' %the session you want to analyse
instructions.ds='DART_V1_YM90K_Celine' %the datasheet where information about this session is stored
instructions.refDay='2' %which of the two days should serve as the reference for matching. 
instructions.tIdxSource = "PD"; % PD, MW_pd, MW or cS
instructions.tDropBool = false; % true or false of whether there are any trials to drop from either day
instructions.tDropRefDay = []; % which trials to drop from the ref day, which in this contect always means instructions.session
instructions.tDropMatchDay = [];%which trials to drop from the other day
instructions.tDropAction = 'DEL'; % DEL or REP
%% notes