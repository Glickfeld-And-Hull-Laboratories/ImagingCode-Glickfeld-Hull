%% i2508 220131 - baseline 1
 
expt(1).mouse = 'i2508';
expt(1).date = '220131';
expt(2).experiment = 'baseline'
expt(1).folder = '2P_images'; %%%%%%%
expt(1).img_loc  = {'V1';'L2/3'};
expt(1).z = -204.68;
expt(1).obj = '16x';
expt(1).zoom = 2;
expt(1).frame_rate = 15;
expt(1).multiday_time_days = 0;
expt(1).multiday_matchdays = [];
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
expt(1).time = {'1335','1359','1434'};
expt(1).eye_str = {'Ipsi','Contra','Ipsi'};
expt(1).data_loc = 'Tierney';

%% i2508 220204 - baseline 2
 
expt(2).mouse = 'i2508';
expt(2).date = '220204';
expt(2).experiment = 'baseline'
expt(2).folder = '2P_images'; %%%%%%%
expt(2).img_loc  = {'V1';'L2/3'};
expt(2).z = -200;
expt(2).obj = '16x';
expt(2).zoom = 2;
expt(2).frame_rate = 15;
expt(2).multiday_time_days = 4; %%%%%%%
expt(2).multiday_matchdays = [1];
expt(2).regImgStartFrame = 39808; %%%%%%%
expt(2).motionThreshold = 0.05;
expt(2).areaBorders = 0;
expt(2).img_strct  = {'cells'};
expt(2).indicator = {'virus';'GCaMP8s'};
expt(2).drug = nan; %%%%%%%
expt(2).greenredsimultaneous = 0; %%%%%%%
expt(2).greenChannelLabel = 'ANY'; %%%%%%%
expt(2).rettuning = {'002';'0947'}; %%%%%%%
expt(2).stimruns = ['002';'003';'004'];
expt(2).time = {'0957','1016','1051'};
expt(2).eye_str = {'Ipsi','Contra','Ipsi'};
expt(2).data_loc = 'Tierney';

%% i2508 220210 - post-MD
 
expt(3).mouse = 'i2508';
expt(3).date = '220210';
expt(3).experiment = 'MD'
expt(3).folder = '2P_images'; %%%%%%%
expt(3).img_loc  = {'V1';'L2/3'};
expt(3).z = -209.37;
expt(3).obj = '16x';
expt(3).zoom = 2;
expt(3).frame_rate = 15;
expt(3).multiday_time_days = 6;
expt(3).multiday_matchdays = [2]; % input should be the comparison day
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
expt(3).time = {'1540','1559','1634'};
expt(3).eye_str = {'Ipsi','Contra','Ipsi'};
expt(3).data_loc = 'Tierney';