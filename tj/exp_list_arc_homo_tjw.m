%this experiment list contains imaging done in the green channel for the
%Arc homozygous mouse line; single imaging sessions (not matched)

%%something is not right on here!!!***
%%
expt = [];
frameRateHz = 15;

%1 is 2566
%%
% i2566 - Arc Homo - Session 1; 15min blocks of 16dir, blank screen, 1dir, blank screen
expt(1).mouse = 'i2566';
expt(1).sex = 'male';
expt(1).date = '230914';
expt(1).ref_date = '230914';
expt(1).img_point = '1';
expt(1).img_loc  = {'V1';'L2_3'};
expt(1).img_strct  = {'cells'};
expt(1).green_indicator = {'AAV';'GFP-Arc'};
expt(1).time_mat = ['1112'];
expt(1).runs = ['003'];
expt(1).nrun = size(expt(1).runs,1);
expt(1).folder = 'tj';
expt(1).z = 241.4;
expt(1).obj = '25x';
expt(1).zoom = 5.7;

%%
% i2567 - Arc Homo - Session 1; 15min blocks of 16dir, blank screen, 1dir, blank screen
expt(2).mouse = 'i2567';
expt(2).sex = 'male';
expt(2).date = '230915';
expt(2).ref_date = '230915';
expt(2).img_point = '1';
expt(2).img_loc  = {'V1';'L2_3'};
expt(2).img_strct  = {'cells'};
expt(2).green_indicator = {'AAV';'GFP-Arc'};
expt(2).time_mat = ['1324'];
expt(2).runs = ['003'];
expt(2).nrun = size(expt(2).runs,1);
expt(2).folder = 'tj';
expt(2).z = 214.06;
expt(2).obj = '25x';
expt(2).zoom = 5.7;

%%
% i2566 - Arc Homo - Session 2; 15min blocks of 16dir, blank screen, 1dir,
% blank screen - only took 2000 frames at start of each block
expt(3).mouse = 'i2566';
expt(3).sex = 'male';
expt(3).date = '230924';
expt(3).ref_date = '230924';
expt(3).img_point = '2';
expt(3).img_loc  = {'V1';'L2_3'};
expt(3).img_strct  = {'cells'};
expt(3).green_indicator = {'AAV';'GFP-Arc'};
expt(3).time_mat = ['1023'];
expt(3).runs = ['003'];
expt(3).nrun = size(expt(3).runs,1);
expt(3).folder = 'tj';
expt(3).z = 203.9;
expt(3).obj = '25x';
expt(3).zoom = 6.7;

%%
% i2567 - Arc Homo - Session 2; 15min blocks of 16dir, blank screen, 1dir,
% blank screen - only took 2000 frames at start of each block
expt(4).mouse = 'i2567';
expt(4).sex = 'male';
expt(4).date = '230924';
expt(4).ref_date = '230924';
expt(4).img_point = '2';
expt(4).img_loc  = {'V1';'L2_3'};
expt(4).img_strct  = {'cells'};
expt(4).green_indicator = {'AAV';'GFP-Arc'};
expt(4).time_mat = ['1144'];
expt(4).runs = ['003'];
expt(4).nrun = size(expt(4).runs,1);
expt(4).folder = 'tj';
expt(4).z = 213.28;
expt(4).obj = '25x';
expt(4).zoom = 6.7;

%%
% i2566 - Arc Homo - Session 3; 60min of 1dir (0deg) - took 2000 frame
% stacks ~every 15min
expt(5).mouse = 'i2566';
expt(5).sex = 'male';
expt(5).date = '231004';
expt(5).ref_date = '231004';
expt(5).img_point = '3';
expt(5).img_loc  = {'V1';'L2_3'};
expt(5).img_strct  = {'cells'};
expt(5).green_indicator = {'AAV';'GFP-Arc'};
expt(5).time_mat = ['0800'];
expt(5).runs = ['003'];
expt(5).nrun = size(expt(5).runs,1);
expt(5).folder = 'tj';
expt(5).z = 206.25;
expt(5).obj = '25x';
expt(5).zoom = 5.7;

%%
% i2566 - Arc Homo - Session 4; 60min of 1dir (0deg) - took 2000 frame
% stacks ~every 15min --> ***mouse was in dark overnight
expt(6).mouse = 'i2566';
expt(6).sex = 'male';
expt(6).date = '231006';
expt(6).ref_date = '231006';
expt(6).img_point = '4';
expt(6).img_loc  = {'V1';'L2_3'};
expt(6).img_strct  = {'cells'};
expt(6).green_indicator = {'AAV';'GFP-Arc'};
expt(6).time_mat = ['0843'];
expt(6).runs = ['001'];
expt(6).nrun = size(expt(6).runs,1);
expt(6).folder = 'tj';
expt(6).z = 207.81;
expt(6).obj = '25x';
expt(6).zoom = 5.7;

%%
% i2567 - Arc Homo - Session 3; 60min of 1dir (0deg); all frames collected
% --> ***mouse was in dark for 4d
expt(7).mouse = 'i2567';
expt(7).sex = 'male';
expt(7).date = '231018';
expt(7).ref_date = '231018';
expt(7).img_point = '3';
expt(7).img_loc  = {'V1';'L2_3'};
expt(7).img_strct  = {'cells'};
expt(7).green_indicator = {'AAV';'GFP-Arc'};
expt(7).time_mat = ['0757'];
expt(7).runs = ['001'];
expt(7).nrun = size(expt(7).runs,1);
expt(7).folder = 'tj';
expt(7).z = 170.31;
expt(7).obj = '25x';
expt(7).zoom = 5.7;

%%
% i2570 - Arc Homo - Session 1; 60min of 1dir (0deg); all frames collected
%1st mouse injected w/ non-fluorescent Cre
expt(8).mouse = 'i2570';
expt(8).sex = 'male';
expt(8).date = '231103';
expt(8).ref_date = '231103';
expt(8).img_point = '1';
expt(8).img_loc  = {'V1';'L2_3'};
expt(8).img_strct  = {'cells'};
expt(8).green_indicator = {'AAV';'GFP-Arc'};
expt(8).time_mat = ['0846'];
expt(8).runs = ['001'];
expt(8).nrun = size(expt(8).runs,1);
expt(8).folder = 'tj';
expt(8).z = 200.78;
expt(8).obj = '25x';
expt(8).zoom = 5.7;

%%



%%
% i2571 - Arc Homo - Session 1; 60min of 1dir (0deg); all frames collected
%non-fluorescent Cre
expt(10).mouse = 'i2571';
expt(10).sex = 'male';
expt(10).date = '231105';
expt(10).ref_date = '231105';
expt(10).img_point = '1';
expt(10).img_loc  = {'V1';'L2_3'};
expt(10).img_strct  = {'cells'};
expt(10).green_indicator = {'AAV';'GFP-Arc'};
expt(10).time_mat = ['0931'];
expt(10).runs = ['001'];
expt(10).nrun = size(expt(10).runs,1);
expt(10).folder = 'tj';
expt(10).z = 200.78;
expt(10).obj = '25x';
expt(10).zoom = 5.7;


%%
% i2571 - Arc Homo - Session 2; 60min of 1dir (0deg); 2000 frames every 15
% min due to attempt at imaging L4
expt(11).mouse = 'i2571';
expt(11).sex = 'male';
expt(11).date = '231108';
expt(11).ref_date = '231108';
expt(11).img_point = '2';
expt(11).img_loc  = {'V1';'L2_3'};
expt(11).img_strct  = {'cells'};
expt(11).green_indicator = {'AAV';'GFP-Arc'};
expt(11).time_mat = ['0858'];
expt(11).runs = ['001'];
expt(11).nrun = size(expt(11).runs,1);
expt(11).folder = 'tj';
expt(11).z = 404.68;
expt(11).obj = '25x';
expt(11).zoom = 5.7;


%%
% i2571 - Arc Homo - Session 3; 60min of 1dir (0deg); mouse in dark
% overnight
expt(12).mouse = 'i2571';
expt(12).sex = 'male';
expt(12).date = '231109';
expt(12).ref_date = '231109';
expt(12).img_point = '3';
expt(12).img_loc  = {'V1';'L2_3'};
expt(12).img_strct  = {'cells'};
expt(12).green_indicator = {'AAV';'GFP-Arc'};
expt(12).time_mat = ['0834'];
expt(12).runs = ['001'];
expt(12).nrun = size(expt(12).runs,1);
expt(12).folder = 'tj';
expt(12).z = 185.93;
expt(12).obj = '25x';
expt(12).zoom = 5.7;

%%
% i2570 - Arc Homo - Session 3; 60min of 1dir (0deg); mouse in dark
% for 4d; non-fluorescent Cre
expt(13).mouse = 'i2570';
expt(13).sex = 'male';
expt(13).date = '231110';
expt(13).ref_date = '231110';
expt(13).img_point = '2';
expt(13).img_loc  = {'V1';'L2_3'};
expt(13).img_strct  = {'cells'};
expt(13).green_indicator = {'AAV';'GFP-Arc'};
expt(13).time_mat = ['0825'];
expt(13).runs = ['001'];
expt(13).nrun = size(expt(13).runs,1);
expt(13).folder = 'tj';
expt(13).z = 164.84;
expt(13).obj = '25x';
expt(13).zoom = 5.7;












