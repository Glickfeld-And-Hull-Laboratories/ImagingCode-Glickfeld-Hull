%% INPUT MATCHDAYS FOR ALL EXPT EXCEPT 2509


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
expt(1).multiday_time_days = 0;
expt(1).multiday_matchdays = [1];
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
expt(2).multiday_time_days = 7;
expt(2).multiday_matchdays = [1]; % input should be the comparison day
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
expt(3).multiday_time_days = 0;
expt(4).matchday_baseline = []; % input should be the comparison day
expt(5).matchday_MD = [4]; % input should be the comparison day
expt(4).matchday_recovery = [5]; % input should be the comparison day
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
expt(5).matchday_MD = []; % input should be the comparison day
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

%% i2510 220223 - baseline 2
 
expt(6).mouse = 'i2510';
expt(6).date = '220223';
expt(6).experiment = 'baseline'
expt(6).folder = '2P_images'; %%%%%%%
expt(6).img_loc  = {'V1';'L2/3'};
expt(6).z = -141.4;
expt(6).obj = '16x';
expt(6).zoom = 2;
expt(6).frame_rate = 15;
expt(6).multiday_time_days = 0;
expt(6).multiday_matchdays = [];
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
expt(6).time = {'0848','0907','0942'};
expt(6).eye_str = {'Ipsi','Contra','Ipsi'};
expt(6).data_loc = 'Tierney';

%% i2511 220224 - baseline 2
 
expt(8).mouse = 'i2511';
expt(8).date = '220224';
expt(8).experiment = 'baseline2'
expt(8).folder = '2P_images'; %%%%%%%
expt(8).img_loc  = {'V1';'L2/3'};
expt(8).z = -148.43;
expt(8).obj = '16x';
expt(8).zoom = 2;
expt(8).frame_rate = 15;
expt(8).multiday_time_days = 0;
expt(8).multiday_matchdays = [];
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
expt(8).time = {'1240','1300','1335'};
expt(8).eye_str = {'Ipsi','Contra','Ipsi'};
expt(8).data_loc = 'Tierney';