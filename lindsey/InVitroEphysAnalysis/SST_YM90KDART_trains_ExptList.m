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

