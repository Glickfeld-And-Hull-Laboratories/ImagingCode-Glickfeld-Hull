% Twelve directions, two contrasts of gratings (so really 4 types of
% plaids)
% 1s on, 3s off

params.frameRate = 15;


%% i1392 240104, V1
expt(1).mouse = 'i1392';
expt(1).date = '240104';
expt(1).img_loc  = {'V1'};
expt(1).z = -200;
expt(1).img_strct  = {'cells'};
expt(1).driver = {'SLC'};
expt(1).indicator = {'tg';'GCaMP6s'};
expt(1).coFolder = {'002'};
expt(1).coTime = {'1340'};
expt(1).SF = 0.05;
expt(1).TF = 2;
expt(1).saveLoc = 'sara';

%% i2714 240315, V1
expt(2).mouse = 'i2714';
expt(2).date = '240315';
expt(2).img_loc  = {'V1'};
expt(2).z = -342;
expt(2).img_strct  = {'cells'};
expt(2).driver = {'Scnn1a'};
expt(2).indicator = {'FLEX';'GCaMP8s'};
expt(2).coFolder = {'002'};
expt(2).coTime = {'0945'};
expt(2).SF = 0.05;
expt(2).TF = 2;
expt(2).saveLoc = 'sara';


%% i1392 240329, AL
expt(3).mouse = 'i1392';
expt(3).date = '240329';
expt(3).img_loc  = {'AL'};
expt(3).z = -200;
expt(3).img_strct  = {'cells'};
expt(3).driver = {'SLC'};
expt(3).indicator = {'tg';'GCaMP6s'};
expt(3).coFolder = {'002'};
expt(3).coTime = {'1034'};
expt(3).SF = 0.05;
expt(3).TF = 2;
expt(3).saveLoc = 'sara';

%% i1395 240412, V1
expt(4).mouse = 'i1395';
expt(4).date = '240412';
expt(4).img_loc  = {'V1'};
expt(4).z = -200;
expt(4).img_strct  = {'cells'};
expt(4).driver = {'SLC'};
expt(4).indicator = {'tg';'GCaMP6s'};
expt(4).coFolder = {'002'};
expt(4).coTime = {'1205'};
expt(4).SF = 0.05;
expt(4).TF = 2;
expt(4).saveLoc = 'sara';
