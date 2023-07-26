%this experiment list is a combination of the long (7 session) dark/light mice, as well as the 4
%session mice that is indicated in the middle of the list
%%
expt = [];
frameRateHz = 15.5;
%img_point: 1 = baseline1, 2 = baseline1 + 4d, 3 = baseline2, 4 = baseline2
%+ 4hrs, 5 = post 4d dark, 6 = post dark + 4hrs, 7 = post dark + 7d

%1-7 is 2537, 8-14 is 2538, 15-21 is 2543
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

%%
%% NEW IMAGING TIMELINE 1 = baseline1; 2 = baseline2; 3 = postdark 4 = postdark + 7d
%%
%% i2554 - SOM/CBA - Baseline 1
expt(22).mouse = 'i2554';
expt(22).sex = 'male';
expt(22).date = '230320';
expt(22).ref_date = '230320';
expt(22).img_point = '1';
expt(22).img_loc  = {'V1';'L2_3'};
expt(22).img_strct  = {'cells'};
expt(22).green_indicator = {'AAV';'GCaMP7f'};
expt(22).time_mat = ['1139'];
expt(22).runs = ['001'];
expt(22).nrun = size(expt(22).runs,1);
expt(22).folder = 'tj';
expt(22).z = 200.00;
expt(22).obj = '16x';
expt(22).zoom = 2.0;

%% i2554 - SOM/CBA - Baseline 2
expt(23).mouse = 'i2554';
expt(23).sex = 'male';
expt(23).date = '230327';
expt(23).ref_date = '230320';
expt(23).img_point = '2';
expt(23).img_loc  = {'V1';'L2_3'};
expt(23).img_strct  = {'cells'};
expt(23).green_indicator = {'AAV';'GCaMP7f'};
expt(23).time_mat = ['1141'];
expt(23).runs = ['001'];
expt(23).nrun = size(expt(23).runs,1);
expt(23).folder = 'tj';
expt(23).z = 201.56;
expt(23).obj = '16x';
expt(23).zoom = 2.0;

%% i2554 - SOM/CBA - Post-dark
expt(24).mouse = 'i2554';
expt(24).sex = 'male';
expt(24).date = '230403';
expt(24).ref_date = '230320';
expt(24).img_point = '3';
expt(24).img_loc  = {'V1';'L2_3'};
expt(24).img_strct  = {'cells'};
expt(24).green_indicator = {'AAV';'GCaMP7f'};
expt(24).time_mat = ['1217'];
expt(24).runs = ['001'];
expt(24).nrun = size(expt(24).runs,1);
expt(24).folder = 'tj';
expt(24).z = 197.65;
expt(24).obj = '16x';
expt(24).zoom = 2.0;

%% i2554 - SOM/CBA - Post-dark 1 week
expt(25).mouse = 'i2554';
expt(25).sex = 'male';
expt(25).date = '230410';
expt(25).ref_date = '230320';
expt(25).img_point = '4';
expt(25).img_loc  = {'V1';'L2_3'};
expt(25).img_strct  = {'cells'};
expt(25).green_indicator = {'AAV';'GCaMP7f'};
expt(25).time_mat = ['0847'];
expt(25).runs = ['001'];
expt(25).nrun = size(expt(25).runs,1);
expt(25).folder = 'tj';
expt(25).z = 200.78;
expt(25).obj = '16x';
expt(25).zoom = 2.0;


%% i2555 - SOM/CBA - Baseline 1
expt(26).mouse = 'i2555';
expt(26).sex = 'male';
expt(26).date = '230320';
expt(26).ref_date = '230320';
expt(26).img_point = '1';
expt(26).img_loc  = {'V1';'L2_3'};
expt(26).img_strct  = {'cells'};
expt(26).green_indicator = {'AAV';'GCaMP7f'};
expt(26).time_mat = ['1255'];
expt(26).runs = ['001'];
expt(26).nrun = size(expt(26).runs,1);
expt(26).folder = 'tj';
expt(26).z = 232.03;
expt(26).obj = '16x';
expt(26).zoom = 2.0;

%% i2555 - SOM/CBA - Baseline 2
expt(27).mouse = 'i2555';
expt(27).sex = 'male';
expt(27).date = '230327';
expt(27).ref_date = '230320';
expt(27).img_point = '2';
expt(27).img_loc  = {'V1';'L2_3'};
expt(27).img_strct  = {'cells'};
expt(27).green_indicator = {'AAV';'GCaMP7f'};
expt(27).time_mat = ['1304'];
expt(27).runs = ['001'];
expt(27).nrun = size(expt(27).runs,1);
expt(27).folder = 'tj';
expt(27).z = 249.21;
expt(27).obj = '16x';
expt(27).zoom = 2.0;

%% i2555 - SOM/CBA - Post-dark
expt(28).mouse = 'i2555';
expt(28).sex = 'male';
expt(28).date = '230403';
expt(28).ref_date = '230320';
expt(28).img_point = '3';
expt(28).img_loc  = {'V1';'L2_3'};
expt(28).img_strct  = {'cells'};
expt(28).green_indicator = {'AAV';'GCaMP7f'};
expt(28).time_mat = ['1344'];
expt(28).runs = ['001'];
expt(28).nrun = size(expt(28).runs,1);
expt(28).folder = 'tj';
expt(28).z = 233.59;
expt(28).obj = '16x';
expt(28).zoom = 2.0;

%% i2555 - SOM/CBA - Post-dark 1 week
expt(29).mouse = 'i2555';
expt(29).sex = 'male';
expt(29).date = '230410';
expt(29).ref_date = '230320';
expt(29).img_point = '4';
expt(29).img_loc  = {'V1';'L2_3'};
expt(29).img_strct  = {'cells'};
expt(29).green_indicator = {'AAV';'GCaMP7f'};
expt(29).time_mat = ['1002'];
expt(29).runs = ['001'];
expt(29).nrun = size(expt(29).runs,1);
expt(29).folder = 'tj';
expt(29).z = 242.96;
expt(29).obj = '16x';
expt(29).zoom = 2.0;

%% i2556 - SOM/CBA - Baseline 1
expt(30).mouse = 'i2556';
expt(30).sex = 'male';
expt(30).date = '230325';
expt(30).ref_date = '230325';
expt(30).img_point = '1';
expt(30).img_loc  = {'V1';'L2_3'};
expt(30).img_strct  = {'cells'};
expt(30).green_indicator = {'AAV';'GCaMP7f'};
expt(30).time_mat = ['0808'];
expt(30).runs = ['001'];
expt(30).nrun = size(expt(30).runs,1);
expt(30).folder = 'tj';
expt(30).z = 228.90;
expt(30).obj = '16x';
expt(30).zoom = 2.0;

%% i2556 - SOM/CBA - Baseline 2
expt(31).mouse = 'i2556';
expt(31).sex = 'male';
expt(31).date = '230329';
expt(31).ref_date = '230325';
expt(31).img_point = '2';
expt(31).img_loc  = {'V1';'L2_3'};
expt(31).img_strct  = {'cells'};
expt(31).green_indicator = {'AAV';'GCaMP7f'};
expt(31).time_mat = ['1218'];
expt(31).runs = ['001'];
expt(31).nrun = size(expt(31).runs,1);
expt(31).folder = 'tj';
expt(31).z = 232.81;
expt(31).obj = '16x';
expt(31).zoom = 2.0;

%% i2556 - SOM/CBA - Post-dark
expt(32).mouse = 'i2556';
expt(32).sex = 'male';
expt(32).date = '230405';
expt(32).ref_date = '230325';
expt(32).img_point = '3';
expt(32).img_loc  = {'V1';'L2_3'};
expt(32).img_strct  = {'cells'};
expt(32).green_indicator = {'AAV';'GCaMP7f'};
expt(32).time_mat = ['1044'];
expt(32).runs = ['001'];
expt(32).nrun = size(expt(32).runs,1);
expt(32).folder = 'tj';
expt(32).z = 228.90;
expt(32).obj = '16x';
expt(32).zoom = 2.0;

%% i2556 - SOM/CBA - Post-dark 1 week
expt(33).mouse = 'i2556';
expt(33).sex = 'male';
expt(33).date = '230412';
expt(33).ref_date = '230325';
expt(33).img_point = '4';
expt(33).img_loc  = {'V1';'L2_3'};
expt(33).img_strct  = {'cells'};
expt(33).green_indicator = {'AAV';'GCaMP7f'};
expt(33).time_mat = ['0924'];
expt(33).runs = ['001'];
expt(33).nrun = size(expt(33).runs,1);
expt(33).folder = 'tj';
expt(33).z = 224.21;
expt(33).obj = '16x';
expt(33).zoom = 2.0;

%%
%% COHORT 2
%%

%% i2557 - PV/CBA - Baseline 1
expt(34).mouse = 'i2557';
expt(34).sex = 'male';
expt(34).date = '230512';
expt(34).ref_date = '230512';
expt(34).img_point = '1';
expt(34).img_loc  = {'V1';'L2_3'};
expt(34).img_strct  = {'cells'};
expt(34).green_indicator = {'AAV';'GCaMP7f'};
expt(34).time_mat = ['1249'];
expt(34).runs = ['001'];
expt(34).nrun = size(expt(34).runs,1);
expt(34).folder = 'tj';
expt(34).z = 290.62;
expt(34).obj = '16x';
expt(34).zoom = 2.0;

%% i2557 - PV/CBA - Baseline 2/Pre dark
expt(35).mouse = 'i2557';
expt(35).sex = 'male';
expt(35).date = '230520';
expt(35).ref_date = '230512';
expt(35).img_point = '2';
expt(35).img_loc  = {'V1';'L2_3'};
expt(35).img_strct  = {'cells'};
expt(35).green_indicator = {'AAV';'GCaMP7f'};
expt(35).time_mat = ['0935'];
expt(35).runs = ['001'];
expt(35).nrun = size(expt(35).runs,1);
expt(35).folder = 'tj';
expt(35).z = 285.15;
expt(35).obj = '16x';
expt(35).zoom = 2.0;

%% i2557 - PV/CBA - Post Dark 15 mins
expt(36).mouse = 'i2557';
expt(36).sex = 'male';
expt(36).date = '230525';
expt(36).ref_date = '230512';
expt(36).img_point = '3';
expt(36).img_loc  = {'V1';'L2_3'};
expt(36).img_strct  = {'cells'};
expt(36).green_indicator = {'AAV';'GCaMP7f'};
expt(36).time_mat = ['0757'];
expt(36).runs = ['001'];
expt(36).nrun = size(expt(36).runs,1);
expt(36).folder = 'tj';
expt(36).z = 289.84;
expt(36).obj = '16x';
expt(36).zoom = 2.0;

%% i2557 - PV/CBA - Post Dark 1 week
expt(37).mouse = 'i2557';
expt(37).sex = 'male';
expt(37).date = '230601';
expt(37).ref_date = '230512';
expt(37).img_point = '4';
expt(37).img_loc  = {'V1';'L2_3'};
expt(37).img_strct  = {'cells'};
expt(37).green_indicator = {'AAV';'GCaMP7f'};
expt(37).time_mat = ['0829'];
expt(37).runs = ['001'];
expt(37).nrun = size(expt(37).runs,1);
expt(37).folder = 'tj';
expt(37).z = 283.59;
expt(37).obj = '16x';
expt(37).zoom = 2.0;


%%
%%LM STARTS HERE*****
%%
%% i2558 - PV/CBA - Baseline 1
expt(38).mouse = 'i2558';
expt(38).sex = 'male';
expt(38).date = '230512';
expt(38).ref_date = '230512';
expt(38).img_point = '1';
expt(38).img_loc  = {'LM';'L2_3'};
expt(38).img_strct  = {'cells'};
expt(38).green_indicator = {'AAV';'GCaMP7f'};
expt(38).time_mat = ['1406'];
expt(38).runs = ['001'];
expt(38).nrun = size(expt(38).runs,1);
expt(38).folder = 'tj';
expt(38).z = 325.78;
expt(38).obj = '16x';
expt(38).zoom = 2.0;

%% i2558 - PV/CBA - Baseline 2/Pre dark
expt(39).mouse = 'i2558';
expt(39).sex = 'male';
expt(39).date = '230520';
expt(39).ref_date = '230512';
expt(39).img_point = '2';
expt(39).img_loc  = {'LM';'L2_3'};
expt(39).img_strct  = {'cells'};
expt(39).green_indicator = {'AAV';'GCaMP7f'};
expt(39).time_mat = ['1056'];
expt(39).runs = ['001'];
expt(39).nrun = size(expt(39).runs,1);
expt(39).folder = 'tj';
expt(39).z = 323.43;
expt(39).obj = '16x';
expt(39).zoom = 2.0;

%% i2558 - PV/CBA - Post Dark 15 mins
expt(40).mouse = 'i2558';
expt(40).sex = 'male';
expt(40).date = '230525';
expt(40).ref_date = '230512';
expt(40).img_point = '3';
expt(40).img_loc  = {'LM';'L2_3'};
expt(40).img_strct  = {'cells'};
expt(40).green_indicator = {'AAV';'GCaMP7f'};
expt(40).time_mat = ['1054'];
expt(40).runs = ['001'];
expt(40).nrun = size(expt(40).runs,1);
expt(40).folder = 'tj';
expt(40).z = 318.75;
expt(40).obj = '16x';
expt(40).zoom = 2.0;

%% i2558 - PV/CBA - Post Dark 1 week
expt(41).mouse = 'i2558';
expt(41).sex = 'male';
expt(41).date = '230601';
expt(41).ref_date = '230512';
expt(41).img_point = '4';
expt(41).img_loc  = {'LM';'L2_3'};
expt(41).img_strct  = {'cells'};
expt(41).green_indicator = {'AAV';'GCaMP7f'};
expt(41).time_mat = ['0944'];
expt(41).runs = ['001'];
expt(41).nrun = size(expt(41).runs,1);
expt(41).folder = 'tj';
expt(41).z = 325.78;
expt(41).obj = '16x';
expt(41).zoom = 2.0;

%% i2559 - PV/CBA - Baseline 1
expt(42).mouse = 'i2559';
expt(42).sex = 'male';
expt(42).date = '230513';
expt(42).ref_date = '230513';
expt(42).img_point = '1';
expt(42).img_loc  = {'LM';'L2_3'};
expt(42).img_strct  = {'cells'};
expt(42).green_indicator = {'AAV';'GCaMP7f'};
expt(42).time_mat = ['0814'];
expt(42).runs = ['001'];
expt(42).nrun = size(expt(42).runs,1);
expt(42).folder = 'tj';
expt(42).z = 328.12;
expt(42).obj = '16x';
expt(42).zoom = 2.0;

%% i2559 - PV/CBA - Baseline 2/Pre dark
expt(43).mouse = 'i2559';
expt(43).sex = 'male';
expt(43).date = '230521';
expt(43).ref_date = '230513';
expt(43).img_point = '2';
expt(43).img_loc  = {'LM';'L2_3'};
expt(43).img_strct  = {'cells'};
expt(43).green_indicator = {'AAV';'GCaMP7f'};
expt(43).time_mat = ['0917'];
expt(43).runs = ['001'];
expt(43).nrun = size(expt(43).runs,1);
expt(43).folder = 'tj';
expt(43).z = 342.96;
expt(43).obj = '16x';
expt(43).zoom = 2.0;

%% i2559 - PV/CBA - Post Dark 15 mins
expt(44).mouse = 'i2559';
expt(44).sex = 'male';
expt(44).date = '230526';
expt(44).ref_date = '230513';
expt(44).img_point = '3';
expt(44).img_loc  = {'LM';'L2_3'};
expt(44).img_strct  = {'cells'};
expt(44).green_indicator = {'AAV';'GCaMP7f'};
expt(44).time_mat = ['0827'];
expt(44).runs = ['001'];
expt(44).nrun = size(expt(44).runs,1);
expt(44).folder = 'tj';
expt(44).z = 329.68;
expt(44).obj = '16x';
expt(44).zoom = 2.0;

%% i2559 - PV/CBA - Post Dark 1 week
expt(45).mouse = 'i2559';
expt(45).sex = 'male';
expt(45).date = '230602';
expt(45).ref_date = '230513';
expt(45).img_point = '4';
expt(45).img_loc  = {'LM';'L2_3'};
expt(45).img_strct  = {'cells'};
expt(45).green_indicator = {'AAV';'GCaMP7f'};
expt(45).time_mat = ['0825'];
expt(45).runs = ['001'];
expt(45).nrun = size(expt(45).runs,1);
expt(45).folder = 'tj';
expt(45).z = 336.71;
expt(45).obj = '16x';
expt(45).zoom = 2.0;

%% i2560 - PV/CBA - Baseline 1
expt(46).mouse = 'i2560';
expt(46).sex = 'male';
expt(46).date = '230513';
expt(46).ref_date = '230513';
expt(46).img_point = '1';
expt(46).img_loc  = {'LM';'L2_3'};
expt(46).img_strct  = {'cells'};
expt(46).green_indicator = {'AAV';'GCaMP7f'};
expt(46).time_mat = ['0931'];
expt(46).runs = ['001'];
expt(46).nrun = size(expt(46).runs,1);
expt(46).folder = 'tj';
expt(46).z = 287.50;
expt(46).obj = '16x';
expt(46).zoom = 2.0;

%% i2560 - PV/CBA - Baseline 2/Pre dark
expt(47).mouse = 'i2560';
expt(47).sex = 'male';
expt(47).date = '230521';
expt(47).ref_date = '230513';
expt(47).img_point = '2';
expt(47).img_loc  = {'LM';'L2_3'};
expt(47).img_strct  = {'cells'};
expt(47).green_indicator = {'AAV';'GCaMP7f'};
expt(47).time_mat = ['1035'];
expt(47).runs = ['001'];
expt(47).nrun = size(expt(47).runs,1);
expt(47).folder = 'tj';
expt(47).z = 282.81;
expt(47).obj = '16x';
expt(47).zoom = 2.0;

%% i2560 - PV/CBA - Post Dark 15 mins
expt(48).mouse = 'i2560';
expt(48).sex = 'male';
expt(48).date = '230526';
expt(48).ref_date = '230513';
expt(48).img_point = '3';
expt(48).img_loc  = {'LM';'L2_3'};
expt(48).img_strct  = {'cells'};
expt(48).green_indicator = {'AAV';'GCaMP7f'};
expt(48).time_mat = ['0950'];
expt(48).runs = ['001'];
expt(48).nrun = size(expt(48).runs,1);
expt(48).folder = 'tj';
expt(48).z = 278.12;
expt(48).obj = '16x';
expt(48).zoom = 2.0;

%% i2560 - PV/CBA - Post Dark 1 week
expt(49).mouse = 'i2560';
expt(49).sex = 'male';
expt(49).date = '230602';
expt(49).ref_date = '230513';
expt(49).img_point = '4';
expt(49).img_loc  = {'LM';'L2_3'};
expt(49).img_strct  = {'cells'};
expt(49).green_indicator = {'AAV';'GCaMP7f'};
expt(49).time_mat = ['1016'];
expt(49).runs = ['001'];
expt(49).nrun = size(expt(49).runs,1);
expt(49).folder = 'tj';
expt(49).z = 285.15;
expt(49).obj = '16x';
expt(49).zoom = 2.0;




























