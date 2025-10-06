%% EXAMPLE i2205 control
% expt(1).mouse = 'i2205';
% expt(1).date = '250521';
% expt(1).img_loc  = {'V1';'L2/3'};
% expt(1).obj = '16x';
% expt(1).zoom = 2;
% expt(1).frame_rate = 15;
% expt(1).multiday_timesincedrug_hours = '0';
% expt(1).multiday_matchdays = [];
% expt(1).indicator = {'GCaMP8S'};
% expt(1).drug = 'YM90K-DART';
% expt(1).greenChannelLabel = 'any';
% expt(1).redChannelLabel = 'VIP';
% expt(1).redChannelTag = {'tdTomato','NES-HTP'};
% expt(1).redChannelRun = '001';
% expt(1).contrastxori_runs = {'004'};
% expt(1).contrastxori_time = {'1036'};
% expt(1).wheelColor = 'orange';
% expt(1).exptType = 'VIP_YM90K';
%% i2776_D1 Pre test <-- Enter descriptive title (optional but helps with clarity)
expt(1).mouse = 'i2776'; % <-- enter mouse ID number
expt(1).date = '250929'; % <-- enter date in the form yymmdd
expt(1).img_loc  = {'V1';'L2/3'}; %usually in the format {'region';'layer'}
expt(1).obj = '16x'; % <-- change if using the 25X objective
expt(1).zoom = 1.7;% <--  enter the magnification 
expt(1).frame_rate = 15; % <-- enter the 2p imaging framerate
expt(1).multiday_trainingday_days = '1'; % <-- enter the day since training starts
expt(1).multiday_matchdays = [4]; % <-- indicate other expt numbers (e.g., this is expt 1) that are matched with this one.
expt(1).indicator = {'GCaMP8f'}; % <-- enter GCaMP version used
%expt(1).drug = ''; % <-- enter [DART] drug used, or modify this to reflect other relevant manipulations
expt(1).greenChannelLabel = 'any'; % <-- if using a cell-type specific GCaMP, list the cell type here. Otherwise this is typically "any"
%expt(1).redChannelLabel = ''; % <-- if using a [cell-type specific] red indicator, list the cell type here. 
%expt(1).redChannelTag = {}; % <-- if using a red indicator, list the indicator, such as tdTomato.
%expt(1).redChannelRun = ''; % <-- if using a red indicator and there is a red snapshot to identify those cells, list the 2p data folder for that snapshot run here. 
expt(1).dir_runs = {'002'}; % <-- enter the 2p data folder for the run with visual stimuli to be analyzed.
expt(1).dir_time = {'1233'}; % <-- enter the timestamp for the mWorks file associated with the run lsited above.
expt(1).wheelColor = 'orange';  % <-- must be one of the colors used to differentiate running wheels in wheelSpeedCalc.m. Typically "orange".
expt(1).exptType = ''; % <-- if you use different folders for different projects, enter that folder name here
%% Copy the above and make a new entry for each new experiment, ie each time you collect data.
% Change the expt id number for each new entry (expt(2), expt(3), etc).
% Most DART experiments will have at least two entries, one for the control
% and one for the DART session. If you collect different datasets on the
% same day/mouse, such as collecting from two cortical regions, it is
% probably easiest to make a separate entry for each.
%% i2777_D1 Pre test_plane 175
expt(2).mouse = 'i2777'; % <-- enter mouse ID number
expt(2).date = '251002'; % <-- enter date in the form yymmdd
expt(2).img_loc  = {'V1';'L2/3'}; %usually in the format {'region';'layer'}
expt(2).obj = '16x'; % <-- change if using the 25X objective
expt(2).zoom = 2;% <--  enter the magnification 
expt(2).frame_rate = 15; % <-- enter the 2p imaging framerate
expt(2).multiday_trainingday_days = '1'; % <-- enter the day since training starts
expt(2).multiday_matchdays = [5]; % <-- indicate other expt numbers (e.g., this is expt 1) that are matched with this one.
expt(2).indicator = {'GCaMP8f'}; % <-- enter GCaMP version used
expt(2).greenChannelLabel = 'any'; % <-- if using a cell-type specific GCaMP, list the cell type here. Otherwise this is typically "any"
expt(2).dir_runs = {'002'}; % <-- enter the 2p data folder for the run with visual stimuli to be analyzed.
expt(2).dir_time = {'1208'}; % <-- enter the timestamp for the mWorks file associated with the run lsited above.
expt(2).wheelColor = 'orange';  % <-- must be one of the colors used to differentiate running wheels in wheelSpeedCalc.m. Typically "orange".
expt(2).exptType = ''; % <-- if you use different folders for different projects, enter that folder name here
%% i2777_D1 Pre test_plane 211.71
expt(3).mouse = 'i2777';
expt(3).date = '251002';
expt(3).img_loc  = {'V1';'L2/3'};
expt(3).obj = '16x';
expt(3).zoom = 2; 
expt(3).frame_rate = 15;
expt(3).multiday_trainingday_days = '1';
expt(3).multiday_matchdays = [6];
expt(3).indicator = {'GCaMP8f'};
expt(3).greenChannelLabel = 'any';
expt(3).dir_runs = {'005'};
expt(3).dir_time = {'1326'};
expt(3).wheelColor = 'orange';
expt(3).exptType = '';
%% i2776_D5 Post test
expt(4).mouse = 'i2776';
expt(4).date = '251003';
expt(4).img_loc  = {'V1';'L2/3'};
expt(4).obj = '16x';
expt(4).zoom = 1.7; 
expt(4).frame_rate = 15;
expt(4).multiday_trainingday_days = '5';
expt(4).multiday_matchdays = [1];
expt(4).indicator = {'GCaMP8f'};
expt(4).greenChannelLabel = 'any';
expt(4).dir_runs = {'001'};
expt(4).dir_time = {'1533'};
expt(4).wheelColor = 'orange';
expt(4).exptType = '';
%% i2777_D5 Post test_plane 175
expt(5).mouse = 'i2777';
expt(5).date = '251006';
expt(5).img_loc  = {'V1';'L2/3'};
expt(5).obj = '16x';
expt(5).zoom = 2; 
expt(5).frame_rate = 15;
expt(5).multiday_trainingday_days = '5';
expt(5).multiday_matchdays = [2];
expt(5).indicator = {'GCaMP8f'};
expt(5).greenChannelLabel = 'any';
%expt(5).dir_runs = {'001'};
%expt(5).dir_time = {'1533'};
expt(5).wheelColor = 'orange';
expt(5).exptType = '';
%% i2777_D5 Post test_plane 211.71
expt(6).mouse = 'i2777';
expt(6).date = '251006';
expt(6).img_loc  = {'V1';'L2/3'};
expt(6).obj = '16x';
expt(6).zoom = 2; 
expt(6).frame_rate = 15;
expt(6).multiday_trainingday_days = '5';
expt(6).multiday_matchdays = [3];
expt(6).indicator = {'GCaMP8f'};
expt(6).greenChannelLabel = 'any';
%expt(6).dir_runs = {'001'};
%expt(6).dir_time = {'1533'};
expt(6).wheelColor = 'orange';
expt(6).exptType = '';