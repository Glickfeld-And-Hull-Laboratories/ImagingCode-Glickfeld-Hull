expt = [];
frameRateHz = 15.5;
%img_point: 1 = baseline1, 2 = baseline1 + 4d, 3 = baseline2, 4 = baseline2
%+ 4hrs, 5 = post 4d dark, 6 = post dark + 4hrs, 7 = post dark + 7d

%1-7 is 2537, 8-14 is 2538
%%
%% i2537 - C57 - Baseline1
expt(1).mouse = 'i2537';
expt(1).sex = 'male';
expt(1).date = '221209';
expt(1).ref_date = '221209';
expt(1).img_point = '1';
expt(1).img_loc  = {'V1';'L2_3'};
expt(1).img_strct  = {'cells'};
expt(1).green_indicator = {'AAV';'GCaMP7f'};
expt(1).time_mat = ['0754'];
expt(1).runs = ['001'];
expt(1).nrun = size(expt(1).runs,1);
expt(1).folder = 'tj';
expt(1).z = 220.31;
expt(1).obj = '16x';
expt(1).zoom = 2.0;

%% i2537 - C57 - Baseline1 + 4d
expt(2).mouse = 'i2537';
expt(2).sex = 'male';
expt(2).date = '221212';
expt(2).ref_date = '221209';
expt(2).img_point = '2';
expt(2).img_loc  = {'V1';'L2_3'};
expt(2).img_strct  = {'cells'};
expt(2).green_indicator = {'AAV';'GCaMP7f'};
expt(2).time_mat = ['1057'];
expt(2).runs = ['001'];
expt(2).nrun = size(expt(2).runs,1);
expt(2).folder = 'tj';
expt(2).z = 223.43;
expt(2).obj = '16x';
expt(2).zoom = 2.0;

%% i2537 - C57 - Baseline2
expt(3).mouse = 'i2537';
expt(3).sex = 'male';
expt(3).date = '221216';
expt(3).ref_date = '221209';
expt(3).img_point = '3';
expt(3).img_loc  = {'V1';'L2_3'};
expt(3).img_strct  = {'cells'};
expt(3).green_indicator = {'AAV';'GCaMP7f'};
expt(3).time_mat = ['0807'];
expt(3).runs = ['001'];
expt(3).nrun = size(expt(3).runs,1);
expt(3).folder = 'tj';
expt(3).z = 228.12;
expt(3).obj = '16x';
expt(3).zoom = 2.0;

%% i2537 - C57 - Baseline2 + 4hr
expt(4).mouse = 'i2537';
expt(4).sex = 'male';
expt(4).date = '221216';
expt(4).ref_date = '221216';
expt(4).img_point = '4';
expt(4).img_loc  = {'V1';'L2_3'};
expt(4).img_strct  = {'cells'};
expt(4).green_indicator = {'AAV';'GCaMP7f'};
expt(4).time_mat = ['1201'];
expt(4).runs = ['002'];
expt(4).nrun = size(expt(4).runs,1);
expt(4).folder = 'tj';
expt(4).z = 220.31;
expt(4).obj = '16x';
expt(4).zoom = 2.0;

%% i2537 - C57 - Post dark 4d
expt(5).mouse = 'i2537';
expt(5).sex = 'male';
expt(5).date = '221220';
expt(5).ref_date = '221216';
expt(5).img_point = '5';
expt(5).img_loc  = {'V1';'L2_3'};
expt(5).img_strct  = {'cells'};
expt(5).green_indicator = {'AAV';'GCaMP7f'};
expt(5).time_mat = ['1022'];
expt(5).runs = ['001'];
expt(5).nrun = size(expt(5).runs,1);
expt(5).folder = 'tj';
expt(5).z = 222.65;
expt(5).obj = '16x';
expt(5).zoom = 2.0;

%% i2537 - C57 - Post dark 4d + 4hr light
expt(6).mouse = 'i2537';
expt(6).sex = 'male';
expt(6).date = '221220';
expt(6).ref_date = '221220';
expt(6).img_point = '6';
expt(6).img_loc  = {'V1';'L2_3'};
expt(6).img_strct  = {'cells'};
expt(6).green_indicator = {'AAV';'GCaMP7f'};
expt(6).time_mat = ['1416'];
expt(6).runs = ['002'];
expt(6).nrun = size(expt(6).runs,1);
expt(6).folder = 'tj';
expt(6).z = 221.87;
expt(6).obj = '16x';
expt(6).zoom = 2.0;

%% i2537 - C57 - Post dark + 7d
expt(7).mouse = 'i2537';
expt(7).sex = 'male';
expt(7).date = '221227';
expt(7).ref_date = '221220';
expt(7).img_point = '7';
expt(7).img_loc  = {'V1';'L2_3'};
expt(7).img_strct  = {'cells'};
expt(7).green_indicator = {'AAV';'GCaMP7f'};
expt(7).time_mat = ['0938'];
expt(7).runs = ['001'];
expt(7).nrun = size(expt(7).runs,1);
expt(7).folder = 'tj';
expt(7).z = 221.09;
expt(7).obj = '16x';
expt(7).zoom = 2.0;

%% i2538 - C57 - Baseline1
expt(8).mouse = 'i2538';
expt(8).sex = 'male';
expt(8).date = '221209';
expt(8).ref_date = '221209';
expt(8).img_point = '1';
expt(8).img_loc  = {'V1';'L2_3'};
expt(8).img_strct  = {'cells'};
expt(8).green_indicator = {'AAV';'GCaMP7f'};
expt(8).time_mat = ['0917'];
expt(8).runs = ['001'];
expt(8).nrun = size(expt(8).runs,1);
expt(8).folder = 'tj';
expt(8).z = 314.84;
expt(8).obj = '16x';
expt(8).zoom = 2.0;

%% i2538 - C57 - Baseline1 + 4d
expt(9).mouse = 'i2538';
expt(9).sex = 'male';
expt(9).date = '221212';
expt(9).ref_date = '221209';
expt(9).img_point = '2';
expt(9).img_loc  = {'V1';'L2_3'};
expt(9).img_strct  = {'cells'};
expt(9).green_indicator = {'AAV';'GCaMP7f'};
expt(9).time_mat = ['1214'];
expt(9).runs = ['001'];
expt(9).nrun = size(expt(9).runs,1);
expt(9).folder = 'tj';
expt(9).z = 314.84;
expt(9).obj = '16x';
expt(9).zoom = 2.0;

%% i2538 - C57 - Baseline2
expt(10).mouse = 'i2538';
expt(10).sex = 'male';
expt(10).date = '221216';
expt(10).ref_date = '221209';
expt(10).img_point = '3';
expt(10).img_loc  = {'V1';'L2_3'};
expt(10).img_strct  = {'cells'};
expt(10).green_indicator = {'AAV';'GCaMP7f'};
expt(10).time_mat = ['0924'];
expt(10).runs = ['001'];
expt(10).nrun = size(expt(10).runs,1);
expt(10).folder = 'tj';
expt(10).z = 307.81;
expt(10).obj = '16x';
expt(10).zoom = 2.0;

%% i2538 - C57 - Baseline2 + 4hr
expt(11).mouse = 'i2538';
expt(11).sex = 'male';
expt(11).date = '221216';
expt(11).ref_date = '221216';
expt(11).img_point = '4';
expt(11).img_loc  = {'V1';'L2_3'};
expt(11).img_strct  = {'cells'};
expt(11).green_indicator = {'AAV';'GCaMP7f'};
expt(11).time_mat = ['1328'];
expt(11).runs = ['002'];
expt(11).nrun = size(expt(11).runs,1);
expt(11).folder = 'tj';
expt(11).z = 324.21;
expt(11).obj = '16x';
expt(11).zoom = 2.0;

%% i2538 - C57 - Post dark 4d
expt(12).mouse = 'i2538';
expt(12).sex = 'male';
expt(12).date = '221220';
expt(12).ref_date = '221216';
expt(12).img_point = '5';
expt(12).img_loc  = {'V1';'L2_3'};
expt(12).img_strct  = {'cells'};
expt(12).green_indicator = {'AAV';'GCaMP7f'};
expt(12).time_mat = ['0743'];
expt(12).runs = ['001'];
expt(12).nrun = size(expt(12).runs,1);
expt(12).folder = 'tj';
expt(12).z = 324.21;
expt(12).obj = '16x';
expt(12).zoom = 2.0;

%% i2538 - C57 - Post dark 4d + 4hr light
expt(13).mouse = 'i2538';
expt(13).sex = 'male';
expt(13).date = '221220';
expt(13).ref_date = '221220';
expt(13).img_point = '6';
expt(13).img_loc  = {'V1';'L2_3'};
expt(13).img_strct  = {'cells'};
expt(13).green_indicator = {'AAV';'GCaMP7f'};
expt(13).time_mat = ['1141'];
expt(13).runs = ['002'];
expt(13).nrun = size(expt(13).runs,1);
expt(13).folder = 'tj';
expt(13).z = 327.34;
expt(13).obj = '16x';
expt(13).zoom = 2.0;

%% i2538 - C57 - Post dark 4d + 7d
expt(14).mouse = 'i2538';
expt(14).sex = 'male';
expt(14).date = '221227';
expt(14).ref_date = '221220';
expt(14).img_point = '7';
expt(14).img_loc  = {'V1';'L2_3'};
expt(14).img_strct  = {'cells'};
expt(14).green_indicator = {'AAV';'GCaMP7f'};
expt(14).time_mat = ['1057'];
expt(14).runs = ['001'];
expt(14).nrun = size(expt(14).runs,1);
expt(14).folder = 'tj';
expt(14).z = 327.24;
expt(14).obj = '16x';
expt(14).zoom = 2.0;

%% i2543 - C57 - Baseline1
expt(15).mouse = 'i2543';
expt(15).sex = 'male';
expt(15).date = '230113';
expt(15).ref_date = '230113';
expt(15).img_point = '1';
expt(15).img_loc  = {'V1';'L2_3'};
expt(15).img_strct  = {'cells'};
expt(15).green_indicator = {'AAV';'GCaMP7f'};
expt(15).time_mat = ['0840'];
expt(15).runs = ['001'];
expt(15).nrun = size(expt(15).runs,1);
expt(15).folder = 'tj';
expt(15).z = 194.53;
expt(15).obj = '16x';
expt(15).zoom = 2.0;

%% i2543 - C57 - Baseline1+4d
expt(16).mouse = 'i2543';
expt(16).sex = 'male';
expt(16).date = '230117';
expt(16).ref_date = '230113';
expt(16).img_point = '2';
expt(16).img_loc  = {'V1';'L2_3'};
expt(16).img_strct  = {'cells'};
expt(16).green_indicator = {'AAV';'GCaMP7f'};
expt(16).time_mat = ['0737'];
expt(16).runs = ['001'];
expt(16).nrun = size(expt(16).runs,1);
expt(16).folder = 'tj';
expt(16).z = 211.71;
expt(16).obj = '16x';
expt(16).zoom = 2.0;

%% i2543 - C57 - Baseline2
expt(17).mouse = 'i2543';
expt(17).sex = 'male';
expt(17).date = '230120';
expt(17).ref_date = '230113';
expt(17).img_point = '3';
expt(17).img_loc  = {'V1';'L2_3'};
expt(17).img_strct  = {'cells'};
expt(17).green_indicator = {'AAV';'GCaMP7f'};
expt(17).time_mat = ['0659'];
expt(17).runs = ['001'];
expt(17).nrun = size(expt(17).runs,1);
expt(17).folder = 'tj';
expt(17).z = 193.75;
expt(17).obj = '16x';
expt(17).zoom = 2.0;

%% i2543 - C57 - Baseline2 + 4hr 
expt(18).mouse = 'i2543';
expt(18).sex = 'male';
expt(18).date = '230120';
expt(18).ref_date = '230120';
expt(18).img_point = '4';
expt(18).img_loc  = {'V1';'L2_3'};
expt(18).img_strct  = {'cells'};
expt(18).green_indicator = {'AAV';'GCaMP7f'};
expt(18).time_mat = ['1050'];
expt(18).runs = ['002'];
expt(18).nrun = size(expt(18).runs,1);
expt(18).folder = 'tj';
expt(18).z = 188.28;
expt(18).obj = '16x';
expt(18).zoom = 2.0;

%% i2543 - C57 - Post dark
expt(19).mouse = 'i2543';
expt(19).sex = 'male';
expt(19).date = '230124';
expt(19).ref_date = '230120';
expt(19).img_point = '5';
expt(19).img_loc  = {'V1';'L2_3'};
expt(19).img_strct  = {'cells'};
expt(19).green_indicator = {'AAV';'GCaMP7f'};
expt(19).time_mat = ['0743'];
expt(19).runs = ['001'];
expt(19).nrun = size(expt(19).runs,1);
expt(19).folder = 'tj';
expt(19).z = 189.06;
expt(19).obj = '16x';
expt(19).zoom = 2.0;

%% i2543 - C57 - Post dark + 4hr light
expt(20).mouse = 'i2543';
expt(20).sex = 'male';
expt(20).date = '230124';
expt(20).ref_date = '230124';
expt(20).img_point = '6';
expt(20).img_loc  = {'V1';'L2_3'};
expt(20).img_strct  = {'cells'};
expt(20).green_indicator = {'AAV';'GCaMP7f'};
expt(20).time_mat = ['1132'];
expt(20).runs = ['002'];
expt(20).nrun = size(expt(20).runs,1);
expt(20).folder = 'tj';
expt(20).z = 188.28;
expt(20).obj = '16x';
expt(20).zoom = 2.0;

%% i2543 - C57 - Post dark + 7d light
expt(21).mouse = 'i2543';
expt(21).sex = 'male';
expt(21).date = '230127';
expt(21).ref_date = '230124';
expt(21).img_point = '7';
expt(21).img_loc  = {'V1';'L2_3'};
expt(21).img_strct  = {'cells'};
expt(21).green_indicator = {'AAV';'GCaMP7f'};
expt(21).time_mat = ['0805'];
expt(21).runs = ['001'];
expt(21).nrun = size(expt(21).runs,1);
expt(21).folder = 'tj';
expt(21).z = 195.31;
expt(21).obj = '16x';
expt(21).zoom = 2.0;












