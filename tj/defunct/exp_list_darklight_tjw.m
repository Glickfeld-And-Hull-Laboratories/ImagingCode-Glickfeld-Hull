%This experiment list is from the defunctportion of the dark light experiment where there are not
%useful comparisons to make due to the lack of standardization between control and experiment
%periods
%%
expt = [];
frameRateHz = 15.5;
%img_point: 1 = pre dark baseline, 2 = post dark, 3 = post 4hr light, 4 = 7
%day light
%% i2537 - C57 - Baseline in V1
expt(1).mouse = 'i2537';
expt(1).sex = 'male';
expt(1).date = '221128';
expt(1).ref_date = '221128';
expt(1).img_point = '1';
expt(1).img_loc  = {'V1';'L2_3'};
expt(1).img_strct  = {'cells'};
expt(1).green_indicator = {'AAV';'GCaMP7f'};
expt(1).time_mat = ['1407'];
expt(1).runs = ['001'];
expt(1).nrun = size(expt(1).runs,1);
expt(1).folder = 'tj';
expt(1).z = 230.46;
expt(1).obj = '16x';
expt(1).zoom = 2.0;
expt(1).wheel = 1;
%% i2538 - C57 - Baseline in V1
expt(2).mouse = 'i2538';
expt(2).sex = 'male';
expt(2).date = '221128';
expt(2).ref_date = '221128';
expt(2).img_point = '1';
expt(2).img_loc  = {'V1';'L2_3'};
expt(2).img_strct  = {'cells'};
expt(2).green_indicator = {'AAV';'GCaMP7f'};
expt(2).time_mat = ['1454'];
expt(2).runs = ['001'];
expt(2).nrun = size(expt(2).runs,1);
expt(2).folder = 'tj';
expt(2).z = 310.15;
expt(2).obj = '16x';
expt(2).zoom = 2.0;
expt(2).wheel = 1;
%% i2537 - C57 - Post-Dark in V1
expt(3).mouse = 'i2537';
expt(3).sex = 'male';
expt(3).date = '221202';
expt(3).ref_date = '221202';
expt(3).img_point = '2';
expt(3).img_loc  = {'V1';'L2_3'};
expt(3).img_strct  = {'cells'};
expt(3).green_indicator = {'AAV';'GCaMP7f'};
expt(3).time_mat = ['0804'];
expt(3).runs = ['001'];
expt(3).nrun = size(expt(3).runs,1);
expt(3).folder = 'tj';
expt(3).z = 220.31;
expt(3).obj = '16x';
expt(3).zoom = 2.0;
expt(3).wheel = 1;
%% i2538 - C57 - Post-Dark in V1
expt(4).mouse = 'i2538';
expt(4).sex = 'male';
expt(4).date = '221202';
expt(4).ref_date = '221202';
expt(4).img_point = '2';
expt(4).img_loc  = {'V1';'L2_3'};
expt(4).img_strct  = {'cells'};
expt(4).green_indicator = {'AAV';'GCaMP7f'};
expt(4).time_mat = ['0858'];
expt(4).runs = ['001'];
expt(4).nrun = size(expt(4).runs,1);
expt(4).folder = 'tj';
expt(4).z = 302.54;
expt(4).obj = '16x';
expt(4).zoom = 2.0;
expt(4).wheel = 1;
%% i2537 - C57 - Post-Light in V1
expt(5).mouse = 'i2537';
expt(5).sex = 'male';
expt(5).date = '221202';
expt(5).ref_date = '221202';
expt(5).img_point = '3';
expt(5).img_loc  = {'V1';'L2_3'};
expt(5).img_strct  = {'cells'};
expt(5).green_indicator = {'AAV';'GCaMP7f'};
expt(5).time_mat = ['1204'];
expt(5).runs = ['002'];
expt(5).nrun = size(expt(5).runs,1);
expt(5).folder = 'tj';
expt(5).z = 233.59;
expt(5).obj = '16x';
expt(5).zoom = 2.0;
expt(5).wheel = 1;
%% i2538 - C57 - Post-Light in V1
expt(6).mouse = 'i2538';
expt(6).sex = 'male';
expt(6).date = '221202';
expt(6).ref_date = '221202';
expt(6).img_point = '3';
expt(6).img_loc  = {'V1';'L2_3'};
expt(6).img_strct  = {'cells'};
expt(6).green_indicator = {'AAV';'GCaMP7f'};
expt(6).time_mat = ['1254'];
expt(6).runs = ['002'];
expt(6).nrun = size(expt(6).runs,1);
expt(6).folder = 'tj';
expt(6).z = 320.31;
expt(6).obj = '16x';
expt(6).zoom = 2.0;
expt(6).wheel = 1;
%% i2537 - C57 - 1 Week Post-Light in V1
expt(7).mouse = 'i2537';
expt(7).sex = 'male';
expt(7).date = '221209';
expt(7).ref_date = '221202';
expt(7).img_point = '4';
expt(7).img_loc  = {'V1';'L2_3'};
expt(7).img_strct  = {'cells'};
expt(7).green_indicator = {'AAV';'GCaMP7f'};
expt(7).time_mat = ['0754'];
expt(7).runs = ['001'];
expt(7).nrun = size(expt(7).runs,1);
expt(7).folder = 'tj';
expt(7).z = 220.31;
expt(7).obj = '16x';
expt(7).zoom = 2.0;
expt(7).wheel = 1;
%% i2538 - C57 - 1 Week Post-Light in V1
expt(8).mouse = 'i2538';
expt(8).sex = 'male';
expt(8).date = '221209';
expt(8).ref_date = '221202';
expt(8).img_point = '4';
expt(8).img_loc  = {'V1';'L2_3'};
expt(8).img_strct  = {'cells'};
expt(8).green_indicator = {'AAV';'GCaMP7f'};
expt(8).time_mat = ['0917'];
expt(8).runs = ['001'];
expt(8).nrun = size(expt(8).runs,1);
expt(8).folder = 'tj';
expt(8).z = 314.94;
expt(8).obj = '16x';
expt(8).zoom = 2.0;
expt(8).wheel = 1;













