params.frameRate = 15;

%% i1369 220518 TF- 2Hz, 0.05cpd
expt(1).mouse = 'i1369';
expt(1).date = '220518';
expt(1).img_loc  = {'V1';'L2/3'};
expt(1).inj_loc = {'V1'};
expt(1).z = -200;
expt(1).img_strct  = {'cells'};
expt(1).driver = {'SLC'};
expt(1).indicator = {'tg';'GCaMP6s'};
expt(1).rfFolder = {'002'}; %map spatial receptive field (36 ret)
expt(1).rfTime = {'1616'};
expt(1).prFolder = {'003'}; %phase reversal (simple v. complex)
expt(1).prTime = {'1722'};
expt(1).saveLoc = 'sara';

%%