firstStim_ind = 2735;
window = 20;
nstim = 10;

%% Cell 1
expt(1).date = '250218';
expt(1).abfdate = '2025_02_18';
expt(1).driver = {'SST'};
expt(1).marker = {'NES-HT';'tdTomato'};
expt(1).firstFile = 0;
expt(1).freqList = [10 50 100 nan 50 10 100];
expt(1).drugList = ['c' 'c' 'c' 'w' 'd' 'd' 'd'];

%% Cell 2
expt(2).date = '250218';
expt(2).abfdate = '2025_02_18';
expt(2).driver = {'SST'};
expt(2).marker = {'NES-HT';'tdTomato'};
expt(2).firstFile = 7;
expt(2).freqList = [10 50 100 nan nan 50 10 100];
expt(2).drugList = ['c' 'c' 'c' 'w' 'w' 'd' 'd' 'd'];

%% Cell 3
expt(3).date = '250228';
expt(3).abfdate = '2025_02_28';
expt(3).driver = {'SST'};
expt(3).marker = {'NES-HT';'tdTomato'};
expt(3).firstFile = 8;
expt(3).freqList = [100 50 10 nan nan nan 50 100 10];
expt(3).drugList = ['c' 'c' 'c' 'w' 'w' 'w' 'd' 'd' 'd'];

%% Cell 4
expt(4).date = '250319';
expt(4).abfdate = '2025_03_19';
expt(4).driver = {'SST'};
expt(4).marker = {'NES-HT';'tdTomato'};
expt(4).firstFile = 0;
expt(4).freqList = [100 50 10 nan 10 50 100];
expt(4).drugList = ['c' 'c' 'c' 'w' 'd' 'd' 'd'];
expt(4).omit_run_trace = [6 13; 7 1];

%% Cell 5
expt(5).date = '250319';
expt(5).abfdate = '2025_03_19';
expt(5).driver = {'SST'};
expt(5).marker = {'NES-HT';'tdTomato'};
expt(5).firstFile = 7;
expt(5).freqList = [100 50 10 nan 10 50 100];
expt(5).drugList = ['c' 'c' 'c' 'w' 'd' 'd' 'd'];
expt(5).omit_run_trace = [6 5;6 12;6 14;6 19;7 17];


