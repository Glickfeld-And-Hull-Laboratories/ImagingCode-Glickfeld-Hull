exp_list = [];
exp_list.mouse_mat = {'Y13' 'X32' 'AC39' 'AC44' 'AC39' 'AC44' 'Y26' 'AC45' 'DR9' 'AC42' 'M13' 'M14' 'M22' 'M31' 'M31'};
exp_list.date_mat = {'110509' '110512' '110809' '110810' '110812' '110813' '110821'  '110822' '110922' '111002' '111004' '111127' '111201' '120107' '120107'};
exp_list.runs_mat = {[1:2] [3:4] [1:3] [1:4]  [8:10] [1:4] [7:9] [1:4] [1:4] [4:6] [4:5] [1:4] [1:3] [1:2] [3:4]};
exp_list.prot_mat = {1 2 1 1 2 1 1 1 1 2 2 1 1 1 2};
exp_list.run_mat = {0 0 1 1 1 1 1 0 1 1 1 1 1 1 1};
exp_list.blanks_mat = {1 1 1 1 1 1 1 1 1 1 1 1 1 1 1};
exp_list.dir_mat = {1 1 2 2 2 2 2 2 2 1 1 2 2 1 1};
exp_list.depth_mat = {75 85 75 85 60 75 60 80 95 100 100 90 80 75 150};
exp_list.zoom_mat = {1 1 1 1 1 1.5 1.5 1 1.5 1.5 1 1 1 1.5 1.5};

fn_out = 'G:\users\lindsey\analysisLG\experiments\SF5xTF5_2P_PM\SF5xTF5_2P_V1_PM_exp_list.mat';
save(fn_out, 'exp_list');
