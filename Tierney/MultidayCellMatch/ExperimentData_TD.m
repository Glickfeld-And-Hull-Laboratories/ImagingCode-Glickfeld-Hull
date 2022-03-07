%% i2509 is good, get corrImg and ori data for all other data sets

%% i2508 220204 - baseline 2
 
expt(1).mouse = 'i2508';
expt(1).date = '220204';
expt(1).experiment = 'baseline'
expt(1).folder = '2P_images'; %%%%%%%
expt(1).img_loc  = {'V1';'L2/3'};
expt(1).z = -200;
expt(1).obj = '16x';
expt(1).zoom = 2;
expt(1).frame_rate = 15;
expt(1).multiday_time_days = 0; % days since baseline
expt(1).matchday_baseline = []; % input should be the comparison day
expt(1).matchday_MD = [2]; % input should be the comparison day
expt(1).matchday_recovery = []; % input should be the comparison day
expt(1).regImgStartFrame = 39808; %%%%%%%
expt(1).motionThreshold = 0.05;
expt(1).areaBorders = 0;
expt(1).img_strct  = {'cells'};
expt(1).indicator = {'virus';'GCaMP8s'};
expt(1).drug = nan; %%%%%%%
expt(1).greenredsimultaneous = 0; %%%%%%%
expt(1).greenChannelLabel = 'ANY'; %%%%%%%
expt(1).rettuning = {'002';'0947'}; %%%%%%%
expt(1).stimruns = ['002';'003';'004'];
expt(1).time = {'0957','1016','1051'};
expt(1).eye_str = {'Ipsi','Contra','Ipsi'};
expt(1).data_loc = 'Tierney';

%% i2508 220211 - post-MD 2
 
expt(2).mouse = 'i2508';
expt(2).date = '220211';
expt(2).experiment = 'MD'
expt(2).folder = '2P_images'; %%%%%%%
expt(2).img_loc  = {'V1';'L2/3'};
expt(2).z = -203.90;
expt(2).obj = '16x';
expt(2).zoom = 2;
expt(2).frame_rate = 15;
expt(2).multiday_time_days = 7; % days since baseline
expt(2).matchday_baseline = [1]; % input should be the comparison day
expt(2).matchday_MD = []; % input should be the comparison day
expt(2).matchday_recovery = []; % input should be the comparison day
expt(2).regImgStartFrame = 39808; %%%%%%%
expt(2).motionThreshold = 0.05;
expt(2).areaBorders = 0;
expt(2).img_strct  = {'cells'};
expt(2).indicator = {'virus';'GCaMP8s'};
expt(2).drug = nan; %%%%%%%
expt(2).greenredsimultaneous = 0; %%%%%%%
expt(2).greenChannelLabel = 'ANY'; %%%%%%%
expt(2).rettuning = {'002';'0947'}; %%%%%%%
expt(2).stimruns = ['001';'002';'003'];
expt(2).time = {'1440','1502','1537'};
expt(2).eye_str = {'Ipsi','Contra','Ipsi'};
expt(2).data_loc = 'Tierney';

%% i2509 220215 - baseline 4
 
expt(3).mouse = 'i2509';
expt(3).date = '220215';
expt(3).experiment = 'baseline'
expt(3).folder = '2P_images'; %%%%%%%
expt(3).img_loc  = {'V1';'L2/3'};
expt(3).z = -164.06;
expt(3).obj = '16x';
expt(3).zoom = 2;
expt(3).frame_rate = 15;
expt(3).multiday_time_days = 0; % days since baseline
expt(3).matchday_baseline = []; % input should be the comparison day
expt(3).matchday_MD = [4]; % input should be the comparison day
expt(3).matchday_recovery = [5]; % input should be the comparison day
expt(3).regImgStartFrame = 39808; %%%%%%%
expt(3).motionThreshold = 0.05;
expt(3).areaBorders = 0;
expt(3).img_strct  = {'cells'};
expt(3).indicator = {'virus';'GCaMP8s'};
expt(3).drug = nan; %%%%%%%
expt(3).greenredsimultaneous = 0; %%%%%%%
expt(3).greenChannelLabel = 'ANY'; %%%%%%%
expt(3).rettuning = {'002';'0947'}; %%%%%%%
expt(3).stimruns = ['001';'002';'003'];
expt(3).time = {'1421','1445','1520'};
expt(3).eye_str = {'Ipsi','Contra','Ipsi'};
expt(3).data_loc = 'Tierney';

%% i2509 220222 - post-MD 1
 
expt(4).mouse = 'i2509';
expt(4).date = '220222';
expt(4).experiment = 'MD';
expt(4).folder = '2P_images'; %%%%%%%
expt(4).img_loc  = {'V1';'L2/3'};
expt(4).z = -168.75;
expt(4).obj = '16x';
expt(4).zoom = 2;
expt(4).frame_rate = 15;
expt(4).multiday_time_days = 7;
expt(4).matchday_baseline = [3]; % input should be the comparison day
expt(4).matchday_MD = []; % input should be the comparison day
expt(4).matchday_recovery = [5]; % input should be the comparison day
expt(4).regImgStartFrame = 39808; %%%%%%%
expt(4).motionThreshold = 0.05;
expt(4).areaBorders = 0;
expt(4).img_strct  = {'cells'};
expt(4).indicator = {'virus';'GCaMP8s'};
expt(4).drug = nan; %%%%%%%
expt(4).greenredsimultaneous = 0; %%%%%%%
expt(4).greenChannelLabel = 'ANY'; %%%%%%%
expt(4).rettuning = {'002';'0947'}; %%%%%%%
expt(4).stimruns = ['001';'002';'003'];
expt(4).time = {'1531','1552','1629'};
expt(4).eye_str = {'Ipsi','Contra','Ipsi'};
expt(4).data_loc = 'Tierney';

%% i2509 220225 - post-MD 3
 
expt(5).mouse = 'i2509';
expt(5).date = '220225';
expt(5).experiment = 'recovery'
expt(5).folder = '2P_images'; %%%%%%%
expt(5).img_loc  = {'V1';'L2/3'};
expt(5).z = -165.62;
expt(5).obj = '16x';
expt(5).zoom = 2;
expt(5).frame_rate = 15;
expt(5).multiday_time_days = 10;
expt(5).matchday_baseline = [3]; % input should be the comparison day
expt(5).matchday_MD = [4]; % input should be the comparison day
expt(5).matchday_recovery = []; % input should be the comparison day
expt(5).regImgStartFrame = 39808; %%%%%%%
expt(5).motionThreshold = 0.05;
expt(5).areaBorders = 0;
expt(5).img_strct  = {'cells'};
expt(5).indicator = {'virus';'GCaMP8s'};
expt(5).drug = nan; %%%%%%%
expt(5).greenredsimultaneous = 0; %%%%%%%
expt(5).greenChannelLabel = 'ANY'; %%%%%%%
expt(5).rettuning = {'002';'0947'}; %%%%%%%
expt(5).stimruns = ['001';'002';'003'];
expt(5).time = {'0825','0846','0920'};
expt(5).eye_str = {'Ipsi','Contra','Ipsi'};
expt(5).data_loc = 'Tierney';

%% i2509 220303 - post-MD 4
 
expt(6).mouse = 'i2509';
expt(6).date = '220303';
expt(6).experiment = 'recovery'
expt(6).folder = '2P_images'; %%%%%%%
expt(6).img_loc  = {'V1';'L2/3'};
expt(6).z = -175.00;
expt(6).obj = '16x';
expt(6).zoom = 2;
expt(6).frame_rate = 15;
expt(6).multiday_time_days = 16;
expt(6).matchday_baseline = [3]; % input should be the comparison day
expt(6).matchday_MD = [4]; % input should be the comparison day
expt(6).matchday_recovery = []; % input should be the comparison day
expt(6).regImgStartFrame = 39808; %%%%%%%
expt(6).motionThreshold = 0.05;
expt(6).areaBorders = 0;
expt(6).img_strct  = {'cells'};
expt(6).indicator = {'virus';'GCaMP8s'};
expt(6).drug = nan; %%%%%%%
expt(6).greenredsimultaneous = 0; %%%%%%%
expt(6).greenChannelLabel = 'ANY'; %%%%%%%
expt(6).rettuning = {'002';'0947'}; %%%%%%%
expt(6).stimruns = ['001';'002';'003'];
expt(6).time = {'1023','1041','1116'};
expt(6).eye_str = {'Ipsi','Contra','Ipsi'};
expt(6).data_loc = 'Tierney';


%% i2510 220223 - baseline 2
 
expt(7).mouse = 'i2510';
expt(7).date = '220223';
expt(7).experiment = 'baseline'
expt(7).folder = '2P_images'; %%%%%%%
expt(7).img_loc  = {'V1';'L2/3'};
expt(7).z = -141.4;
expt(7).obj = '16x';
expt(7).zoom = 2;
expt(7).frame_rate = 15;
expt(7).multiday_time_days = 0; % days since baseline
expt(7).matchday_baseline = []; % input should be the comparison day
expt(7).matchday_MD = [8]; % input should be the comparison day
expt(7).matchday_recovery = [9]; % input should be the comparison day
expt(7).regImgStartFrame = 39808; %%%%%%%
expt(7).motionThreshold = 0.05;
expt(7).areaBorders = 0;
expt(7).img_strct  = {'cells'};
expt(7).indicator = {'virus';'GCaMP8s'};
expt(7).drug = nan; %%%%%%%
expt(7).greenredsimultaneous = 0; %%%%%%%
expt(7).greenChannelLabel = 'ANY'; %%%%%%%
expt(7).rettuning = {'002';'0947'}; %%%%%%%
expt(7).stimruns = ['001';'002';'003'];
expt(7).time = {'0848','0907','0942'};
expt(7).eye_str = {'Ipsi','Contra','Ipsi'};
expt(7).data_loc = 'Tierney';

%% i2510 220302 - post-MD 1
 
expt(8).mouse = 'i2510';
expt(8).date = '220302';
expt(8).experiment = 'MD'
expt(8).folder = '2P_images'; %%%%%%%
expt(8).img_loc  = {'V1';'L2/3'};
expt(8).z = -153.90;
expt(8).obj = '16x';
expt(8).zoom = 2;
expt(8).frame_rate = 15;
expt(8).multiday_time_days = 7; % days since baseline
expt(8).matchday_baseline = [7]; % input should be the comparison day
expt(8).matchday_MD = []; % input should be the comparison day
expt(8).matchday_recovery = [9]; % input should be the comparison day
expt(8).regImgStartFrame = 39808; %%%%%%%
expt(8).motionThreshold = 0.05;
expt(8).areaBorders = 0;
expt(8).img_strct  = {'cells'};
expt(8).indicator = {'virus';'GCaMP8s'};
expt(8).drug = nan; %%%%%%%
expt(8).greenredsimultaneous = 0; %%%%%%%
expt(8).greenChannelLabel = 'ANY'; %%%%%%%
expt(8).rettuning = {'002';'0947'}; %%%%%%%
expt(8).stimruns = ['001';'002';'003'];
expt(8).time = {'1412','1433','1513'};
expt(8).eye_str = {'Ipsi','Contra','Ipsi'};
expt(8).data_loc = 'Tierney';

%% i2510 220304 - post-MD 2
 
expt(9).mouse = 'i2510';
expt(9).date = '220304';
expt(9).experiment = 'recovery'
expt(9).folder = '2P_images'; %%%%%%%
expt(9).img_loc  = {'V1';'L2/3'};
expt(9).z = -141.4; %CHANGE THISSSSSSSSSSS
expt(9).obj = '16x';
expt(9).zoom = 2;
expt(9).frame_rate = 15;
expt(9).multiday_time_days = 9; % days since baseline
expt(9).matchday_baseline = [7]; % input should be the comparison day
expt(9).matchday_MD = [8]; % input should be the comparison day
expt(9).matchday_recovery = []; % input should be the comparison day
expt(9).regImgStartFrame = 39808; %%%%%%%
expt(9).motionThreshold = 0.05;
expt(9).areaBorders = 0;
expt(9).img_strct  = {'cells'};
expt(9).indicator = {'virus';'GCaMP8s'};
expt(9).drug = nan; %%%%%%%
expt(9).greenredsimultaneous = 0; %%%%%%%
expt(9).greenChannelLabel = 'ANY'; %%%%%%%
expt(9).rettuning = {'002';'0947'}; %%%%%%%
expt(9).stimruns = ['001';'002';'003'];
expt(9).time = {'0848','0907','0942'}; %CHANGE THISSSSSSSSSSS
expt(9).eye_str = {'Ipsi','Contra','Ipsi'};
expt(9).data_loc = 'Tierney';

%% i2511 220224 - baseline 2
 
expt(10).mouse = 'i2511';
expt(10).date = '220224';
expt(10).experiment = 'baseline2'
expt(10).folder = '2P_images'; %%%%%%%
expt(10).img_loc  = {'V1';'L2/3'};
expt(10).z = -148.43;
expt(10).obj = '16x';
expt(10).zoom = 2;
expt(10).frame_rate = 15;
expt(10).multiday_time_days = 0; % days since baseline
expt(10).matchday_baseline = []; % input should be the comparison day
expt(10).matchday_MD = [11]; % input should be the comparison day
expt(10).matchday_recovery = [12]; % input should be the comparison day
expt(10).regImgStartFrame = 39808; %%%%%%%
expt(10).motionThreshold = 0.05;
expt(10).areaBorders = 0;
expt(10).img_strct  = {'cells'};
expt(10).indicator = {'virus';'GCaMP8s'};
expt(10).drug = nan; %%%%%%%
expt(10).greenredsimultaneous = 0; %%%%%%%
expt(10).greenChannelLabel = 'ANY'; %%%%%%%
expt(10).rettuning = {'002';'0947'}; %%%%%%%
expt(10).stimruns = ['001';'002';'003'];
expt(10).time = {'1240','1300','1335'};
expt(10).eye_str = {'Ipsi','Contra','Ipsi'};
expt(10).data_loc = 'Tierney';

%% i2511 220303 - post-MD 1
 
expt(11).mouse = 'i2511';
expt(11).date = '220303';
expt(11).experiment = 'MD'
expt(11).folder = '2P_images'; %%%%%%%
expt(11).img_loc  = {'V1';'L2/3'};
expt(11).z = -160.93;
expt(11).obj = '16x';
expt(11).zoom = 2;
expt(11).frame_rate = 15;
expt(11).multiday_time_days = 7; % days since baseline
expt(11).matchday_baseline = [10]; % input should be the comparison day
expt(11).matchday_MD = []; % input should be the comparison day
expt(11).matchday_recovery = [12]; % input should be the comparison day
expt(11).regImgStartFrame = 39808; %%%%%%%
expt(11).motionThreshold = 0.05;
expt(11).areaBorders = 0;
expt(11).img_strct  = {'cells'};
expt(11).indicator = {'virus';'GCaMP8s'};
expt(11).drug = nan; %%%%%%%
expt(11).greenredsimultaneous = 0; %%%%%%%
expt(11).greenChannelLabel = 'ANY'; %%%%%%%
expt(11).rettuning = {'002';'0947'}; %%%%%%%
expt(11).stimruns = ['001';'002';'003'];
expt(11).time = {'1156','1226','1303'};
expt(11).eye_str = {'Ipsi','Contra','Ipsi'};
expt(11).data_loc = 'Tierney';

%% i2511 220304 - post-MD 2
 
expt(12).mouse = 'i2511';
expt(12).date = '220304';
expt(12).experiment = 'recovery'
expt(12).folder = '2P_images'; %%%%%%%
expt(12).img_loc  = {'V1';'L2/3'};
expt(12).z = -148.43; %CHANGE THISSSSSSSSSSS
expt(12).obj = '16x';
expt(12).zoom = 2;
expt(12).frame_rate = 15;
expt(12).multiday_time_days = 8; % days since baseline
expt(12).matchday_baseline = [10]; % input should be the comparison day
expt(12).matchday_MD = [11]; % input should be the comparison day
expt(12).matchday_recovery = []; % input should be the comparison day
expt(12).regImgStartFrame = 39808; %%%%%%%
expt(12).motionThreshold = 0.05;
expt(12).areaBorders = 0;
expt(12).img_strct  = {'cells'};
expt(12).indicator = {'virus';'GCaMP8s'};
expt(12).drug = nan; %%%%%%%
expt(12).greenredsimultaneous = 0; %%%%%%%
expt(12).greenChannelLabel = 'ANY'; %%%%%%%
expt(12).rettuning = {'002';'0947'}; %%%%%%%
expt(12).stimruns = ['001';'002';'003'];
expt(12).time = {'1240','1300','1335'}; %CHANGE THISSSSSSSSSSS
expt(12).eye_str = {'Ipsi','Contra','Ipsi'};
expt(12).data_loc = 'Tierney';