% Twelve directions, four phases
% 1s on, 3s off

params.frameRate = 15;

%% i1380 230606, V1
expt(1).mouse = 'i1380';
expt(1).date = '230608';
expt(1).img_loc  = {'V1'};
expt(1).z = -200;
expt(1).img_strct  = {'cells'};
expt(1).driver = {'SLC'};
expt(1).indicator = {'tg';'GCaMP6s'};
expt(1).coFolder = {'002'};
expt(1).coTime = {'1357'};
expt(1).SF = 0.05;
expt(1).TF = 2;
expt(1).saveLoc = 'sara';

%% i1387 230609, V1
expt(2).mouse = 'i1387';
expt(2).date = '230609';
expt(2).img_loc  = {'V1'};
expt(2).z = -200;
expt(2).img_strct  = {'cells'};
expt(2).driver = {'SLC'};
expt(2).indicator = {'tg';'GCaMP6s'};
expt(2).coFolder = {'002'};
expt(2).coTime = {'1042'};
expt(2).SF = 0.05;
expt(2).TF = 2;
expt(2).saveLoc = 'sara';

%% i1386 230609, V1
expt(3).mouse = 'i1386';
expt(3).date = '230609';
expt(3).img_loc  = {'V1'};
expt(3).z = -200;
expt(3).img_strct  = {'cells'};
expt(3).driver = {'SLC'};
expt(3).indicator = {'tg';'GCaMP6s'};
expt(3).coFolder = {'002'};
expt(3).coTime = {'1413'};
expt(3).SF = 0.05;
expt(3).TF = 2;
expt(3).saveLoc = 'sara';

%% i1375 230612, V1
expt(4).mouse = 'i1375';
expt(4).date = '230612';
expt(4).img_loc  = {'V1'};
expt(4).z = -200;
expt(4).img_strct  = {'cells'};
expt(4).driver = {'SLC'};
expt(4).indicator = {'tg';'GCaMP6s'};
expt(4).coFolder = {'002'};
expt(4).coTime = {'1424'};
expt(4).SF = 0.05;
expt(4).TF = 2;
expt(4).saveLoc = 'sara';

%% i1375 230621, V1  -- REALLY ONE PHASE
expt(5).mouse = 'i1375';
expt(5).date = '230621';
expt(5).img_loc  = {'V1'};
expt(5).z = -200;
expt(5).img_strct  = {'cells'};
expt(5).driver = {'SLC'};
expt(5).indicator = {'tg';'GCaMP6s'};
expt(5).coFolder = {'002'};
expt(5).coTime = {'1007'};
expt(5).SF = 0.05;
expt(5).TF = 2;
expt(5).saveLoc = 'sara';

%% i1373 230622, V1
expt(6).mouse = 'i1373';
expt(6).date = '230622';
expt(6).img_loc  = {'V1'};
expt(6).z = -200;
expt(6).img_strct  = {'cells'};
expt(6).driver = {'SLC'};
expt(6).indicator = {'tg';'GCaMP6s'};
expt(6).coFolder = {'002'};
expt(6).coTime = {'1008'};
expt(6).SF = 0.05;
expt(6).TF = 2;
expt(6).saveLoc = 'sara';

%% i2706 230626, V1
expt(7).mouse = 'i2706';
expt(7).date = '230626';
expt(7).img_loc  = {'V1'};
expt(7).z = -180;
expt(7).img_strct  = {'cells'};
expt(7).driver = {'SLC'};
expt(7).indicator = {'tg';'GCaMP8s'};
expt(7).coFolder = {'002'};
expt(7).coTime = {'1330'};
expt(7).SF = 0.05;
expt(7).TF = 2;
expt(7).saveLoc = 'sara';

%% i2706 230710, V1
expt(8).mouse = 'i2706';
expt(8).date = '230710';
expt(8).img_loc  = {'V1'};
expt(8).z = -180;
expt(8).img_strct  = {'cells'};
expt(8).driver = {'SLC'};
expt(8).indicator = {'tg';'GCaMP8s'};
expt(8).coFolder = {'002'};
expt(8).coTime = {'1442'};
expt(8).SF = 0.05;
expt(8).TF = 2;
expt(8).saveLoc = 'sara';

%% i1375 230718, V1 (162000 frames)
expt(9).mouse = 'i1375';
expt(9).date = '230718';
expt(9).img_loc  = {'V1'};
expt(9).z = -200;
expt(9).img_strct  = {'cells'};
expt(9).driver = {'SLC'};
expt(9).indicator = {'tg';'GCaMP6s'};
expt(9).coFolder = {'002'};
expt(9).coTime = {'1327'};
expt(9).SF = 0.05;
expt(9).TF = 2;
expt(9).saveLoc = 'sara';

%% i1373 230721, V1 (162000 frames)
expt(10).mouse = 'i1373';
expt(10).date = '230721';
expt(10).img_loc  = {'V1'};
expt(10).z = -200;
expt(10).img_strct  = {'cells'};
expt(10).driver = {'SLC'};
expt(10).indicator = {'tg';'GCaMP6s'};
expt(10).coFolder = {'002'};
expt(10).coTime = {'1342'};
expt(10).SF = 0.05;
expt(10).TF = 2;
expt(10).saveLoc = 'sara';

%% i2707 230727, V1 (108000 frames)
expt(11).mouse = 'i2707';
expt(11).date = '230727';
expt(11).img_loc  = {'V1'};
expt(11).z = -200;
expt(11).img_strct  = {'cells'};
expt(11).driver = {'syn'};
expt(11).indicator = {'tg';'GCaMP8s'};
expt(11).coFolder = {'002'};
expt(11).coTime = {'1427'};
expt(11).redImg = {'005'};
expt(11).SF = 0.05;
expt(11).TF = 2;
expt(11).saveLoc = 'sara';

%% i2708 230728, V1 (108000 frames)
expt(12).mouse = 'i2708';
expt(12).date = '230728';
expt(12).img_loc  = {'V1'};
expt(12).z = -200;
expt(12).img_strct  = {'cells'};
expt(12).driver = {'syn'};
expt(12).indicator = {'tg';'GCaMP8s'};
expt(12).coFolder = {'002'};
expt(12).coTime = {'0959'};
expt(12).redImg = {'005'};
expt(12).SF = 0.05;
expt(12).TF = 2;
expt(12).saveLoc = 'sara';

%% i1380 230728, V1 (162000 frames)
expt(13).mouse = 'i1380';
expt(13).date = '230728';
expt(13).img_loc  = {'V1'};
expt(13).z = -200;
expt(13).img_strct  = {'cells'};
expt(13).driver = {'SLC'};
expt(13).indicator = {'tg';'GCaMP6s'};
expt(13).coFolder = {'002'};
expt(13).coTime = {'1239'};
expt(13).SF = 0.05;
expt(13).TF = 2;
expt(13).saveLoc = 'sara';

%% i1386 230802, V1 (162000 frames)
expt(14).mouse = 'i1386';
expt(14).date = '230802';
expt(14).img_loc  = {'V1'};
expt(14).z = -200;
expt(14).img_strct  = {'cells'};
expt(14).driver = {'SLC'};
expt(14).indicator = {'tg';'GCaMP6s'};
expt(14).coFolder = {'002'};
expt(14).coTime = {'1038'};
expt(14).SF = 0.05;
expt(14).TF = 2;
expt(14).saveLoc = 'sara';

%% i1387 230803, V1 (162000 frames)
expt(15).mouse = 'i1387';
expt(15).date = '230803';
expt(15).img_loc  = {'V1'};
expt(15).z = -200;
expt(15).img_strct  = {'cells'};
expt(15).driver = {'SLC'};
expt(15).indicator = {'tg';'GCaMP6s'};
expt(15).coFolder = {'002'};
expt(15).coTime = {'1017'};
expt(15).SF = 0.05;
expt(15).TF = 2;
expt(15).saveLoc = 'sara';

%% i1375 230805, PM (162000 frames)
expt(16).mouse = 'i1375';
expt(16).date = '230805';
expt(16).img_loc  = {'PM'};
expt(16).z = -200;
expt(16).img_strct  = {'cells'};
expt(16).driver = {'SLC'};
expt(16).indicator = {'tg';'GCaMP6s'};
expt(16).coFolder = {'002'};
expt(16).coTime = {'1204'};
expt(16).SF = 0.05;
expt(16).TF = 2;
expt(16).saveLoc = 'sara';

%% i1373 230808, AL (162000 frames)
expt(17).mouse = 'i1373';
expt(17).date = '230808';
expt(17).img_loc  = {'AL'};
expt(17).z = -200;
expt(17).img_strct  = {'cells'};
expt(17).driver = {'SLC'};
expt(17).indicator = {'tg';'GCaMP6s'};
expt(17).coFolder = {'002'};
expt(17).coTime = {'1326'};
expt(17).SF = 0.05;
expt(17).TF = 2;
expt(17).saveLoc = 'sara';

%% i1380 230809, AL (162000 frames)
expt(18).mouse = 'i1380';
expt(18).date = '230809';
expt(18).img_loc  = {'AL'};
expt(18).z = -200;
expt(18).img_strct  = {'cells'};
expt(18).driver = {'SLC'};
expt(18).indicator = {'tg';'GCaMP6s'};
expt(18).coFolder = {'002'};
expt(18).coTime = {'1053'};
expt(18).SF = 0.05;
expt(18).TF = 2;
expt(18).saveLoc = 'sara';

%% i1373 230810, LM (162000 frames)
expt(19).mouse = 'i1373';
expt(19).date = '230810';
expt(19).img_loc  = {'LM'};
expt(19).z = -200;
expt(19).img_strct  = {'cells'};
expt(19).driver = {'SLC'};
expt(19).indicator = {'tg';'GCaMP6s'};
expt(19).coFolder = {'002'};
expt(19).coTime = {'1000'};
expt(19).SF = 0.05;
expt(19).TF = 2;
expt(19).saveLoc = 'sara';

%% i1386 230811, AL (162000 frames)
expt(20).mouse = 'i1386';
expt(20).date = '230811';
expt(20).img_loc  = {'AL'};
expt(20).z = -200;
expt(20).img_strct  = {'cells'};
expt(20).driver = {'SLC'};
expt(20).indicator = {'tg';'GCaMP6s'};
expt(20).coFolder = {'002'};
expt(20).coTime = {'1034'};
expt(20).SF = 0.05;
expt(20).TF = 2;
expt(20).saveLoc = 'sara';

%% i1387 230812, AL (162000 frames)
expt(21).mouse = 'i1387';
expt(21).date = '230812';
expt(21).img_loc  = {'AL'};
expt(21).z = -200;
expt(21).img_strct  = {'cells'};
expt(21).driver = {'SLC'};
expt(21).indicator = {'tg';'GCaMP6s'};
expt(21).coFolder = {'002'};
expt(21).coTime = {'1148'};
expt(21).SF = 0.05;
expt(21).TF = 2;
expt(21).saveLoc = 'sara';

%% i1380 230817, LM (162000 frames)
expt(22).mouse = 'i1380';
expt(22).date = '230817';
expt(22).img_loc  = {'AL'};
expt(22).z = -200;
expt(22).img_strct  = {'cells'};
expt(22).driver = {'SLC'};
expt(22).indicator = {'tg';'GCaMP6s'};
expt(22).coFolder = {'002'};
expt(22).coTime = {'1328'};
expt(22).SF = 0.05;
expt(22).TF = 2;
expt(22).saveLoc = 'sara';

%% i2115 230821, V1 (162000 frames)
expt(23).mouse = 'i2115';
expt(23).date = '230821';
expt(23).img_loc  = {'V1'};
expt(23).z = -200;
expt(23).img_strct  = {'cells'};
expt(23).driver = {'syn'};
expt(23).indicator = {'tg';'GCaMP8f'};
expt(23).coFolder = {'002'};
expt(23).coTime = {'1025'};
expt(23).redImg = {'005'};
expt(23).SF = 0.05;
expt(23).TF = 2;
expt(23).saveLoc = 'sara';

%% i1373 230822, V1 (162000 frames) -- ONE PHASE
expt(24).mouse = 'i1373';
expt(24).date = '230822';
expt(24).img_loc  = {'V1'};
expt(24).z = -200;
expt(24).img_strct  = {'cells'};
expt(24).driver = {'SLC'};
expt(24).indicator = {'tg';'GCaMP6s'};
expt(24).coFolder = {'002'};
expt(24).coTime = {'1216'};
expt(24).SF = 0.05;
expt(24).TF = 2;
expt(24).saveLoc = 'sara';

%% i1381 230823, V1 (162000 frames)
expt(25).mouse = 'i1381';
expt(25).date = '230823';
expt(25).img_loc  = {'V1'};
expt(25).z = -200;
expt(25).img_strct  = {'cells'};
expt(25).driver = {'SLC'};
expt(25).indicator = {'tg';'GCaMP6s'};
expt(25).coFolder = {'002'};
expt(25).coTime = {'1021'};
expt(25).SF = 0.05;
expt(25).TF = 2;
expt(25).saveLoc = 'sara';

%% i1389 230825, V1 (162000 frames)
expt(26).mouse = 'i1389';
expt(26).date = '230825';
expt(26).img_loc  = {'V1'};
expt(26).z = -200;
expt(26).img_strct  = {'cells'};
expt(26).driver = {'SLC'};
expt(26).indicator = {'tg';'GCaMP6s'};
expt(26).coFolder = {'002'};
expt(26).coTime = {'1129'};
expt(26).SF = 0.05;
expt(26).TF = 2;
expt(26).saveLoc = 'sara';


%% i1389 230829, AL (162000 frames)
expt(27).mouse = 'i1389';
expt(27).date = '230829';
expt(27).img_loc  = {'AL'};
expt(27).z = -200;
expt(27).img_strct  = {'cells'};
expt(27).driver = {'SLC'};
expt(27).indicator = {'tg';'GCaMP6s'};
expt(27).coFolder = {'002'};
expt(27).coTime = {'1204'};
expt(27).SF = 0.05;
expt(27).TF = 2;
expt(27).saveLoc = 'sara';

%% i1389 230830, LM (162000 frames)
expt(28).mouse = 'i1389';
expt(28).date = '230830';
expt(28).img_loc  = {'LM'};
expt(28).z = -200;
expt(28).img_strct  = {'cells'};
expt(28).driver = {'SLC'};
expt(28).indicator = {'tg';'GCaMP6s'};
expt(28).coFolder = {'002'};
expt(28).coTime = {'1008'};
expt(28).SF = 0.05;
expt(28).TF = 2;
expt(28).saveLoc = 'sara';

%% i2119 230912, V1 (162000 frames)
expt(29).mouse = 'i2119';
expt(29).date = '230912';
expt(29).img_loc  = {'V1'};
expt(29).z = -200;
expt(29).img_strct  = {'cells'};
expt(29).driver = {'syn'};
expt(29).indicator = {'tg';'GCaMP8f'};
expt(29).coFolder = {'002'};
expt(29).coTime = {'1455'};
expt(29).redImg = {'005'};
expt(29).SF = 0.05;
expt(29).TF = 2;
expt(29).saveLoc = 'sara';

%% i1380 230913, LM (162000 frames)
expt(30).mouse = 'i1380';
expt(30).date = '230913';
expt(30).img_loc  = {'LM'};
expt(30).z = -200;
expt(30).img_strct  = {'cells'};
expt(30).driver = {'SLC'};
expt(30).indicator = {'tg';'GCaMP6s'};
expt(30).coFolder = {'002'};
expt(30).coTime = {'1422'};
expt(30).SF = 0.05;
expt(30).TF = 2;
expt(30).saveLoc = 'sara';

%% i1381 230914, V1 (162000 frames) -- ONE PHASE
expt(31).mouse = 'i1381';
expt(31).date = '230914';
expt(31).img_loc  = {'V1'};
expt(31).z = -200;
expt(31).img_strct  = {'cells'};
expt(31).driver = {'SLC'};
expt(31).indicator = {'tg';'GCaMP6s'};
expt(31).coFolder = {'002'};
expt(31).coTime = {'1339'};
expt(31).SF = 0.05;
expt(31).TF = 2;
expt(31).saveLoc = 'sara';

%% i1387 230918, V1 (162000 frames) -- ONE PHASE
expt(32).mouse = 'i1387';
expt(32).date = '230918';
expt(32).img_loc  = {'V1'};
expt(32).z = -200;
expt(32).img_strct  = {'cells'};
expt(32).driver = {'SLC'};
expt(32).indicator = {'tg';'GCaMP6s'};
expt(32).coFolder = {'002'};
expt(32).coTime = {'1033'};
expt(32).SF = 0.05;
expt(32).TF = 2;
expt(32).saveLoc = 'sara';

%% i1386 230918, V1 (162000 frames) -- ONE PHASE
expt(33).mouse = 'i1386';
expt(33).date = '230918';
expt(33).img_loc  = {'V1'};
expt(33).z = -200;
expt(33).img_strct  = {'cells'};
expt(33).driver = {'SLC'};
expt(33).indicator = {'tg';'GCaMP6s'};
expt(33).coFolder = {'002'};
expt(33).coTime = {'1438'};
expt(33).SF = 0.05;
expt(33).TF = 2;
expt(33).saveLoc = 'sara';

%% i1387 230919, V1 (162000 frames)
expt(34).mouse = 'i1387';
expt(34).date = '230919';
expt(34).img_loc  = {'AL'};
expt(34).z = -200;
expt(34).img_strct  = {'cells'};
expt(34).driver = {'SLC'};
expt(34).indicator = {'tg';'GCaMP6s'};
expt(34).coFolder = {'002'};
expt(34).coTime = {'0934'};
expt(34).SF = 0.05;
expt(34).TF = 2;
expt(34).saveLoc = 'sara';

%% i1380 230929, LM (162000 frames)
expt(35).mouse = 'i1380';
expt(35).date = '230929';
expt(35).img_loc  = {'LM'};
expt(35).z = -200;
expt(35).img_strct  = {'cells'};
expt(35).driver = {'SLC'};
expt(35).indicator = {'tg';'GCaMP6s'};
expt(35).coFolder = {'002'};
expt(35).coTime = {'1000'};
expt(35).SF = 0.05;
expt(35).TF = 2;
expt(35).saveLoc = 'sara';

%% i1380 231002, LM (162000 frames)
expt(36).mouse = 'i1380';
expt(36).date = '231002';
expt(36).img_loc  = {'LM'};
expt(36).z = -200;
expt(36).img_strct  = {'cells'};
expt(36).driver = {'SLC'};
expt(36).indicator = {'tg';'GCaMP6s'};
expt(36).coFolder = {'002'};
expt(36).coTime = {'1023'};
expt(36).SF = 0.05;
expt(36).TF = 2;
expt(36).saveLoc = 'sara';

%% i2709 231004, V1 (162000 frames)
expt(37).mouse = 'i2709';
expt(37).date = '231004';
expt(37).img_loc  = {'V1'};
expt(37).z = -200;
expt(37).img_strct  = {'cells'};
expt(37).driver = {'syn'};
expt(37).indicator = {'tg';'GCaMP8s'};
expt(37).coFolder = {'002'};
expt(37).coTime = {'1012'};
expt(37).SF = 0.05;
expt(37).TF = 2;
expt(37).saveLoc = 'sara';

%% i1389 231009, LM (162000 frames)
expt(38).mouse = 'i1389';
expt(38).date = '231009';
expt(38).img_loc  = {'LM'};
expt(38).z = -200;
expt(38).img_strct  = {'cells'};
expt(38).driver = {'SLC'};
expt(38).indicator = {'tg';'GCaMP6s'};
expt(38).coFolder = {'002'};
expt(38).coTime = {'1426'};
expt(38).SF = 0.05;
expt(38).TF = 2;
expt(38).saveLoc = 'sara';

%% i1389 231015, LM (162000 frames)
expt(39).mouse = 'i1389';
expt(39).date = '231015';
expt(39).img_loc  = {'LM'};
expt(39).z = -200;
expt(38).img_strct  = {'cells'};
expt(39).driver = {'SLC'};
expt(39).indicator = {'tg';'GCaMP6s'};
expt(39).coFolder = {'002'};
expt(39).coTime = {'1405'};
expt(39).SF = 0.05;
expt(39).TF = 2;
expt(39).saveLoc = 'sara';

%% i1392 231201, V1 (162000 frames)
expt(40).mouse = 'i1392';
expt(40).date = '231201';
expt(40).img_loc  = {'V1'};
expt(40).z = -200;
expt(40).img_strct  = {'cells'};
expt(40).driver = {'SLC'};
expt(40).indicator = {'tg';'GCaMP6s'};
expt(40).coFolder = {'002'};
expt(40).coTime = {'1450'};
expt(40).SF = 0.05;
expt(40).TF = 2;
expt(40).saveLoc = 'sara';

%% i1391 231204, V1 (162000 frames)
expt(41).mouse = 'i1391';
expt(41).date = '231204';
expt(41).img_loc  = {'V1'};
expt(41).z = -200;
expt(41).img_strct  = {'cells'};
expt(41).driver = {'SLC'};
expt(41).indicator = {'tg';'GCaMP6s'};
expt(41).coFolder = {'002'};
expt(41).coTime = {'1429'};
expt(41).SF = 0.05;
expt(41).TF = 2;
expt(41).saveLoc = 'sara';

%% i1392 231208, LM (162000 frames)
expt(42).mouse = 'i1392';
expt(42).date = '231208';
expt(42).img_loc  = {'LM'};
expt(42).z = -200;
expt(42).img_strct  = {'cells'};
expt(42).driver = {'SLC'};
expt(42).indicator = {'tg';'GCaMP6s'};
expt(42).coFolder = {'002'};
expt(42).coTime = {'1407'};
expt(42).SF = 0.05;
expt(42).TF = 2;
expt(42).saveLoc = 'sara';

%% i1391 231219, V1 (162000 frames)
expt(43).mouse = 'i1391';
expt(43).date = '231219';
expt(43).img_loc  = {'V1'};
expt(43).z = -200;
expt(43).img_strct  = {'cells'};
expt(43).driver = {'SLC'};
expt(43).indicator = {'tg';'GCaMP6s'};
expt(43).coFolder = {'002'};
expt(43).coTime = {'1300'};
expt(43).SF = 0.05;
expt(43).TF = 2;
expt(43).saveLoc = 'sara';

%% i1392 240105, AL (162000 frames)
expt(44).mouse = 'i1392';
expt(44).date = '240105';
expt(44).img_loc  = {'AL'};
expt(44).z = -200;
expt(44).img_strct  = {'cells'};
expt(44).driver = {'SLC'};
expt(44).indicator = {'tg';'GCaMP6s'};
expt(44).coFolder = {'002'};
expt(44).coTime = {'1207'};
expt(44).SF = 0.05;
expt(44).TF = 2;
expt(44).saveLoc = 'sara';

%% i1392 240108, V1 (162000 frames) -- ONE PHASE
expt(45).mouse = 'i1392';
expt(45).date = '240108';
expt(45).img_loc  = {'V1'};
expt(45).z = -200;
expt(45).img_strct  = {'cells'};
expt(45).driver = {'SLC'};
expt(45).indicator = {'tg';'GCaMP6s'};
expt(45).coFolder = {'002'};
expt(45).coTime = {'1326'};
expt(45).SF = 0.05;
expt(45).TF = 2;
expt(45).saveLoc = 'sara';

%% i1380 240109, V1 (162000 frames) -- ONE PHASE
expt(46).mouse = 'i1380';
expt(46).date = '240109';
expt(46).img_loc  = {'V1'};
expt(46).z = -200;
expt(46).img_strct  = {'cells'};
expt(46).driver = {'SLC'};
expt(46).indicator = {'tg';'GCaMP6s'};
expt(46).coFolder = {'002'};
expt(46).coTime = {'1424'};
expt(46).SF = 0.05;
expt(46).TF = 2;
expt(46).saveLoc = 'sara';

%% i1391 240112, V1 (162000 frames) -- ONE PHASE
expt(47).mouse = 'i1391';
expt(47).date = '240112';
expt(47).img_loc  = {'V1'};
expt(47).z = -200;
expt(47).img_strct  = {'cells'};
expt(47).driver = {'SLC'};
expt(47).indicator = {'tg';'GCaMP6s'};
expt(47).coFolder = {'002'};
expt(47).coTime = {'1326'};
expt(47).SF = 0.05;
expt(47).TF = 2;
expt(47).saveLoc = 'sara';

%% i1392 240115, AL (162000 frames) -- ONE PHASE
expt(48).mouse = 'i1392';
expt(48).date = '240115';
expt(48).img_loc  = {'AL'};
expt(48).z = -200;
expt(48).img_strct  = {'cells'};
expt(48).driver = {'SLC'};
expt(48).indicator = {'tg';'GCaMP6s'};
expt(48).coFolder = {'002'};
expt(48).coTime = {'1310'};
expt(48).SF = 0.05;
expt(48).TF = 2;
expt(48).saveLoc = 'sara';

%% i1380 240116, LM (162000 frames) -- ONE PHASE
expt(49).mouse = 'i1380';
expt(49).date = '240116';
expt(49).img_loc  = {'LM'};
expt(49).z = -200;
expt(49).img_strct  = {'cells'};
expt(49).driver = {'SLC'};
expt(49).indicator = {'tg';'GCaMP6s'};
expt(49).coFolder = {'002'};
expt(49).coTime = {'1140'};
expt(49).SF = 0.05;
expt(49).TF = 2;
expt(49).saveLoc = 'sara';

%% i1392 240117, LM (162000 frames) -- ONE PHASE
expt(50).mouse = 'i1392';
expt(50).date = '240117';
expt(50).img_loc  = {'LM'};
expt(50).z = -200;
expt(50).img_strct  = {'cells'};
expt(50).driver = {'SLC'};
expt(50).indicator = {'tg';'GCaMP6s'};
expt(50).coFolder = {'002'};
expt(50).coTime = {'1300'};
expt(50).SF = 0.05;
expt(50).TF = 2;
expt(50).saveLoc = 'sara';

%% i1387 240122, V1 (162000 frames) 
expt(51).mouse = 'i1387';
expt(51).date = '240122';
expt(51).img_loc  = {'V1'};
expt(51).z = -200;
expt(51).img_strct  = {'cells'};
expt(51).driver = {'SLC'};
expt(51).indicator = {'tg';'GCaMP6s'};
expt(51).coFolder = {'002'};
expt(51).coTime = {'1143'};
expt(51).SF = 0.05;
expt(51).TF = 2;
expt(51).saveLoc = 'sara';

%% i1392 240123, V1 (162000 frames) -- ONE PHASE
expt(52).mouse = 'i1392';
expt(52).date = '240123';
expt(52).img_loc  = {'V1'};
expt(52).z = -200;
expt(52).img_strct  = {'cells'};
expt(52).driver = {'SLC'};
expt(52).indicator = {'tg';'GCaMP6s'};
expt(52).coFolder = {'002'};
expt(52).coTime = {'1158'};
expt(52).SF = 0.05;
expt(52).TF = 2;
expt(52).saveLoc = 'sara';

%% i1391 240124, V1 (162000 frames)
expt(53).mouse = 'i1391';
expt(53).date = '240124';
expt(53).img_loc  = {'V1'};
expt(53).z = -200;
expt(53).img_strct  = {'cells'};
expt(53).driver = {'SLC'};
expt(53).indicator = {'tg';'GCaMP6s'};
expt(53).coFolder = {'002'};
expt(53).coTime = {'1443'};
expt(53).SF = 0.05;
expt(53).TF = 2;
expt(53).saveLoc = 'sara';

%% i1381 240125, V1 (162000 frames)
expt(54).mouse = 'i1381';
expt(54).date = '240125';
expt(54).img_loc  = {'V1'};
expt(54).z = -200;
expt(54).img_strct  = {'cells'};
expt(54).driver = {'SLC'};
expt(54).indicator = {'tg';'GCaMP6s'};
expt(54).coFolder = {'002'};
expt(54).coTime = {'1445'};
expt(54).SF = 0.05;
expt(54).TF = 2;
expt(54).saveLoc = 'sara';

%% i1387 240130, V1 (162000 frames) -- ONE PHASE
expt(55).mouse = 'i1387';
expt(55).date = '240130';
expt(55).img_loc  = {'V1'};
expt(55).z = -200;
expt(55).img_strct  = {'cells'};
expt(55).driver = {'SLC'};
expt(55).indicator = {'tg';'GCaMP6s'};
expt(55).coFolder = {'002'};
expt(55).coTime = {'1431'};
expt(55).SF = 0.05;
expt(55).TF = 2;
expt(55).saveLoc = 'sara';


%% i1381 240208, V1 (162000 frames) -- ONE PHASE
expt(56).mouse = 'i1381';
expt(56).date = '240208';
expt(56).img_loc  = {'V1'};
expt(56).z = -200;
expt(56).img_strct  = {'cells'};
expt(56).driver = {'SLC'};
expt(56).indicator = {'tg';'GCaMP6s'};
expt(56).coFolder = {'002'};
expt(56).coTime = {'1400'};
expt(56).SF = 0.05;
expt(56).TF = 2;
expt(56).saveLoc = 'sara';

%% i1387 240209, V1 (162000 frames)
expt(57).mouse = 'i1387';
expt(57).date = '240209';
expt(57).img_loc  = {'AL'};
expt(57).z = -200;
expt(57).img_strct  = {'cells'};
expt(57).driver = {'SLC'};
expt(57).indicator = {'tg';'GCaMP6s'};
expt(57).coFolder = {'002'};
expt(57).coTime = {'1401'};
expt(57).SF = 0.05;
expt(57).TF = 2;
expt(57).saveLoc = 'sara';

%% i1391 240212, AL (162000 frames)
expt(58).mouse = 'i1391';
expt(58).date = '240212';
expt(58).img_loc  = {'V1'};
expt(58).z = -200;
expt(58).img_strct  = {'cells'};
expt(58).driver = {'SLC'};
expt(58).indicator = {'tg';'GCaMP6s'};
expt(58).coFolder = {'002'};
expt(58).coTime = {'1105'};
expt(58).SF = 0.05;
expt(58).TF = 2;
expt(58).saveLoc = 'sara';

%% i1392 240213, AL (162000 frames)
expt(59).mouse = 'i1392';
expt(59).date = '240213';
expt(59).img_loc  = {'AL'};
expt(59).z = -200;
expt(59).img_strct  = {'cells'};
expt(59).driver = {'SLC'};
expt(59).indicator = {'tg';'GCaMP6s'};
expt(59).coFolder = {'002'};
expt(59).coTime = {'1440'};
expt(59).SF = 0.05;
expt(59).TF = 2;
expt(59).saveLoc = 'sara';


%% i1392 240215, LM (162000 frames)
expt(60).mouse = 'i1392';
expt(60).date = '240215';
expt(60).img_loc  = {'LM'};
expt(60).z = -200;
expt(60).img_strct  = {'cells'};
expt(60).driver = {'SLC'};
expt(60).indicator = {'tg';'GCaMP6s'};
expt(60).coFolder = {'002'};
expt(60).coTime = {'1016'};
expt(60).SF = 0.05;
expt(60).TF = 2;
expt(60).saveLoc = 'sara';

%% i2714 240301, V1 L4 (162000 frames)
expt(61).mouse = 'i2714';
expt(61).date = '240301';
expt(61).img_loc  = {'V1'};
expt(61).z = -380;
expt(61).img_strct  = {'cells'};
expt(61).driver = {'Scnn1a'};
expt(61).indicator = {'FLEX';'GCaMP6s'};
expt(61).coFolder = {'002'};
expt(61).coTime = {'1258'};
expt(61).SF = 0.05;
expt(61).TF = 2;
expt(61).saveLoc = 'sara';

%% i2714 240305, V1 L4 (162000 frames) -- ONE PHASE
expt(62).mouse = 'i2714';
expt(62).date = '240305';
expt(62).img_loc  = {'V1'};
expt(62).z = -380;
expt(62).img_strct  = {'cells'};
expt(62).driver = {'Scnn1a'};
expt(62).indicator = {'FLEX';'GCaMP6s'};
expt(62).coFolder = {'002'};
expt(62).coTime = {'1429'};
expt(62).SF = 0.05;
expt(62).TF = 2;
expt(62).saveLoc = 'sara';

%% i2714 240308, V1 L4 (162000 frames)
expt(63).mouse = 'i2714';
expt(63).date = '240308';
expt(63).img_loc  = {'V1'};
expt(63).z = -368;
expt(63).img_strct  = {'cells'};
expt(63).driver = {'Scnn1a'};
expt(63).indicator = {'FLEX';'GCaMP6s'};
expt(63).coFolder = {'002'};
expt(63).coTime = {'1410'};
expt(63).SF = 0.05;
expt(63).TF = 2;
expt(63).saveLoc = 'sara';

%% i2715 240307, V1 L4 (162000 frames)
expt(64).mouse = 'i2715';
expt(64).date = '240307';
expt(64).img_loc  = {'V1'};
expt(64).z = -354;
expt(64).img_strct  = {'cells'};
expt(64).driver = {'Scnn1a'};
expt(64).indicator = {'FLEX';'GCaMP6s'};
expt(64).coFolder = {'002'};
expt(64).coTime = {'1345'};
expt(64).SF = 0.05;
expt(64).TF = 2;
expt(64).saveLoc = 'sara';

%% i2715 240311, V1 L4 (162000 frames)
expt(65).mouse = 'i2715';
expt(65).date = '240311';
expt(65).img_loc  = {'V1'};
expt(65).z = -343;
expt(65).img_strct  = {'cells'};
expt(65).driver = {'Scnn1a'};
expt(65).indicator = {'FLEX';'GCaMP6s'};
expt(65).coFolder = {'002'};
expt(65).coTime = {'1253'};
expt(65).SF = 0.05;
expt(65).TF = 2;
expt(65).saveLoc = 'sara';

%% i2715 240313, V1 L4 (162000 frames) -- ONE PHASE
expt(66).mouse = 'i2715';
expt(66).date = '240313';
expt(66).img_loc  = {'V1'};
expt(66).z = -347;
expt(66).img_strct  = {'cells'};
expt(66).driver = {'Scnn1a'};
expt(66).indicator = {'FLEX';'GCaMP6s'};
expt(66).coFolder = {'002'};
expt(66).coTime = {'1007'};
expt(66).SF = 0.05;
expt(66).TF = 2;
expt(66).saveLoc = 'sara';


%% i1393 240322, V1 (162000 frames)
expt(67).mouse = 'i1393';
expt(67).date = '240322';
expt(67).img_loc  = {'V1'};
expt(67).z = -180;
expt(67).img_strct  = {'cells'};
expt(67).driver = {'SLC'};
expt(67).indicator = {'tg';'GCaMP6s'};
expt(67).coFolder = {'002'};
expt(67).coTime = {'0921'};
expt(67).SF = 0.05;
expt(67).TF = 2;
expt(67).saveLoc = 'sara';

%% i1394 240405, V1 (162000 frames)
expt(68).mouse = 'i1394';
expt(68).date = '240405';
expt(68).img_loc  = {'V1'};
expt(68).z = -200;
expt(68).img_strct  = {'cells'};
expt(68).driver = {'SLC'};
expt(68).indicator = {'tg';'GCaMP6s'};
expt(68).coFolder = {'002'};
expt(68).coTime = {'1228'};
expt(68).SF = 0.05;
expt(68).TF = 2;
expt(68).saveLoc = 'sara';

%% i1395 240409, V1 (162000 frames)
expt(69).mouse = 'i1395';
expt(69).date = '240409';
expt(69).img_loc  = {'V1'};
expt(69).z = -200;
expt(69).img_strct  = {'cells'};
expt(69).driver = {'SLC'};
expt(69).indicator = {'tg';'GCaMP6s'};
expt(69).coFolder = {'002'};
expt(69).coTime = {'0959'};
expt(69).SF = 0.05;
expt(69).TF = 2;
expt(69).saveLoc = 'sara';

%% i1394 240409, AL (162000 frames)
expt(70).mouse = 'i1394';
expt(70).date = '240409';
expt(70).img_loc  = {'AL'};
expt(70).z = -200;
expt(70).img_strct  = {'cells'};
expt(70).driver = {'SLC'};
expt(70).indicator = {'tg';'GCaMP6s'};
expt(70).coFolder = {'002'};
expt(70).coTime = {'1417'};
expt(70).SF = 0.05;
expt(70).TF = 2;
expt(70).saveLoc = 'sara';

%% i1394 240416, V1 (162000 frames)
expt(71).mouse = 'i1394';
expt(71).date = '240416';
expt(71).img_loc  = {'V1'};
expt(71).z = -200;
expt(71).img_strct  = {'cells'};
expt(71).driver = {'SLC'};
expt(71).indicator = {'tg';'GCaMP6s'};
expt(71).coFolder = {'002'};
expt(71).coTime = {'1443'};
expt(71).SF = 0.05;
expt(71).TF = 2;
expt(71).saveLoc = 'sara';

%% i1396 240423, V1 (162000 frames)
expt(72).mouse = 'i1396';
expt(72).date = '240423';
expt(72).img_loc  = {'V1'};
expt(72).z = -200;
expt(72).img_strct  = {'cells'};
expt(72).driver = {'SLC'};
expt(72).indicator = {'tg';'GCaMP8m'};
expt(72).coFolder = {'002'};
expt(72).coTime = {'1435'};
expt(72).SF = 0.05;
expt(72).TF = 2;
expt(72).saveLoc = 'sara';

%% i1396 240501, LM (162000 frames)
expt(73).mouse = 'i1396';
expt(73).date = '240501';
expt(73).img_loc  = {'LM'};
expt(73).z = -200;
expt(73).img_strct  = {'cells'};
expt(73).driver = {'SLC'};
expt(73).indicator = {'tg';'GCaMP8m'};
expt(73).coFolder = {'002'};
expt(73).coTime = {'1331'};
expt(73).SF = 0.05;
expt(73).TF = 2;
expt(73).saveLoc = 'sara';

%% i1396 240502, AL (162000 frames)
expt(74).mouse = 'i1396';
expt(74).date = '240502';
expt(74).img_loc  = {'AL'};
expt(74).z = -200;
expt(74).img_strct  = {'cells'};
expt(74).driver = {'SLC'};
expt(74).indicator = {'tg';'GCaMP8m'};
expt(74).coFolder = {'002'};
expt(74).coTime = {'1348'};
expt(74).SF = 0.05;
expt(74).TF = 2;
expt(74).saveLoc = 'sara';

%% i1396 240507, V1 (162000 frames)
expt(75).mouse = 'i1396';
expt(75).date = '240507';
expt(75).img_loc  = {'V1'};
expt(75).z = -200;
expt(75).img_strct  = {'cells'};
expt(75).driver = {'SLC'};
expt(75).indicator = {'tg';'GCaMP8m'};
expt(75).coFolder = {'002'};
expt(75).coTime = {'1452'};
expt(75).SF = 0.05;
expt(75).TF = 2;
expt(75).saveLoc = 'sara';

%% i1396 240508, AL (162000 frames) -- ONE PHASE
expt(76).mouse = 'i1396';
expt(76).date = '240508';
expt(76).img_loc  = {'AL'};
expt(76).z = -200;
expt(76).img_strct  = {'cells'};
expt(76).driver = {'SLC'};
expt(76).indicator = {'tg';'GCaMP8m'};
expt(76).coFolder = {'002'};
expt(76).coTime = {'0930'};
expt(76).SF = 0.05;
expt(76).TF = 2;
expt(76).saveLoc = 'sara';

%% i1396 240510, LM (162000 frames) -- ONE PHASE
expt(77).mouse = 'i1396';
expt(77).date = '240510';
expt(77).img_loc  = {'LM'};
expt(77).z = -200;
expt(77).img_strct  = {'cells'};
expt(77).driver = {'SLC'};
expt(77).indicator = {'tg';'GCaMP8m'};
expt(77).coFolder = {'002'};
expt(77).coTime = {'0959'};
expt(77).SF = 0.05;
expt(77).TF = 2;
expt(77).saveLoc = 'sara';

%% i1393 240510, LM (162000 frames)
expt(78).mouse = 'i1393';
expt(78).date = '240510';
expt(78).img_loc  = {'LM'};
expt(78).z = -200;
expt(78).img_strct  = {'cells'};
expt(78).driver = {'SLC'};
expt(78).indicator = {'tg';'GCaMP8s'};
expt(78).coFolder = {'002'};
expt(78).coTime = {'1349'};
expt(78).SF = 0.05;
expt(78).TF = 2;
expt(78).saveLoc = 'sara';

%% i1393 240513, V1 (162000 frames)
expt(79).mouse = 'i1393';
expt(79).date = '240513';
expt(79).img_loc  = {'V1'};
expt(79).z = -200;
expt(79).img_strct  = {'cells'};
expt(79).driver = {'SLC'};
expt(79).indicator = {'tg';'GCaMP6s'};
expt(79).coFolder = {'002'};
expt(79).coTime = {'1031'};
expt(79).SF = 0.05;
expt(79).TF = 2;
expt(79).saveLoc = 'sara';

%% i1393 240515, AL (162000 frames) -- ONE PHASE
expt(80).mouse = 'i1393';
expt(80).date = '240515';
expt(80).img_loc  = {'AL'};
expt(80).z = -200;
expt(80).img_strct  = {'cells'};
expt(80).driver = {'SLC'};
expt(80).indicator = {'tg';'GCaMP6s'};
expt(80).coFolder = {'002'};
expt(80).coTime = {'0926'};
expt(80).SF = 0.05;
expt(80).TF = 2;
expt(80).saveLoc = 'sara';

%% i1397 240515, AL (162000 frames) -- ONE PHASE
expt(81).mouse = 'i1397';
expt(81).date = '240515';
expt(81).img_loc  = {'AL'};
expt(81).z = -200;
expt(81).img_strct  = {'cells'};
expt(81).driver = {'SLC'};
expt(81).indicator = {'tg';'GCaMP6s'};
expt(81).coFolder = {'002'};
expt(81).coTime = {'1337'};
expt(81).SF = 0.05;
expt(81).TF = 2;
expt(81).saveLoc = 'sara';

%% i1392 240516, AL (162000 frames) -- ONE PHASE
expt(82).mouse = 'i1392';
expt(82).date = '240516';
expt(82).img_loc  = {'AL'};
expt(82).z = -200;
expt(82).img_strct  = {'cells'};
expt(82).driver = {'SLC'};
expt(82).indicator = {'tg';'GCaMP6s'};
expt(82).coFolder = {'002'};
expt(82).coTime = {'1155'};
expt(82).SF = 0.05;
expt(82).TF = 2;
expt(82).saveLoc = 'sara';

%% i1394 240520, AL (162000 frames) -- ONE PHASE
expt(83).mouse = 'i1394';
expt(83).date = '240520';
expt(83).img_loc  = {'AL'};
expt(83).z = -200;
expt(83).img_strct  = {'cells'};
expt(83).driver = {'SLC'};
expt(83).indicator = {'tg';'GCaMP6s'};
expt(83).coFolder = {'002'};
expt(83).coTime = {'0928'};
expt(83).SF = 0.05;
expt(83).TF = 2;
expt(83).saveLoc = 'sara';

%% i1397 240524, AL (162000 frames) -- ONE PHASE
expt(84).mouse = 'i1397';
expt(84).date = '240524';
expt(84).img_loc  = {'AL'};
expt(84).z = -200;
expt(84).img_strct  = {'cells'};
expt(84).driver = {'SLC'};
expt(84).indicator = {'tg';'GCaMP6s'};
expt(84).coFolder = {'002'};
expt(84).coTime = {'1225'};
expt(84).SF = 0.05;
expt(84).TF = 2;
expt(84).saveLoc = 'sara';

%% i1397 240528, LM (162000 frames) -- ONE PHASE
expt(85).mouse = 'i1397';
expt(85).date = '240528';
expt(85).img_loc  = {'LM'};
expt(85).z = -200;
expt(85).img_strct  = {'cells'};
expt(85).driver = {'SLC'};
expt(85).indicator = {'tg';'GCaMP6s'};
expt(85).coFolder = {'002'};
expt(85).coTime = {'1450'};
expt(85).SF = 0.05;
expt(85).TF = 2;
expt(85).saveLoc = 'sara';

%% i1397 240530, LM (162000 frames) -- ONE PHASE
expt(86).mouse = 'i1397';
expt(86).date = '240530';
expt(86).img_loc  = {'LM'};
expt(86).z = -200;
expt(86).img_strct  = {'cells'};
expt(86).driver = {'SLC'};
expt(86).indicator = {'tg';'GCaMP6s'};
expt(86).coFolder = {'002'};
expt(86).coTime = {'1337'};
expt(86).SF = 0.05;
expt(86).TF = 2;
expt(86).saveLoc = 'sara';

%% i1397 240606, LM (162000 frames)
expt(87).mouse = 'i1397';
expt(87).date = '240606';
expt(87).img_loc  = {'LM'};
expt(87).z = -200;
expt(87).img_strct  = {'cells'};
expt(87).driver = {'SLC'};
expt(87).indicator = {'tg';'GCaMP6s'};
expt(87).coFolder = {'002'};
expt(87).coTime = {'1012'};
expt(87).SF = 0.05;
expt(87).TF = 2;
expt(87).saveLoc = 'sara';

%% i1395 240606, LM (162000 frames)
expt(88).mouse = 'i1395';
expt(88).date = '240606';
expt(88).img_loc  = {'LM'};
expt(88).z = -200;
expt(88).img_strct  = {'cells'};
expt(88).driver = {'SLC'};
expt(88).indicator = {'tg';'GCaMP6s'};
expt(88).coFolder = {'002'};
expt(88).coTime = {'1402'};
expt(88).SF = 0.05;
expt(88).TF = 2;
expt(88).saveLoc = 'sara';

%% i1397 240607, AL (162000 frames)
expt(89).mouse = 'i1397';
expt(89).date = '240607';
expt(89).img_loc  = {'AL'};
expt(89).z = -200;
expt(89).img_strct  = {'cells'};
expt(89).driver = {'SLC'};
expt(89).indicator = {'tg';'GCaMP6s'};
expt(89).coFolder = {'002'};
expt(89).coTime = {'1033'};
expt(89).SF = 0.05;
expt(89).TF = 2;
expt(89).saveLoc = 'sara';

%% i1395 240607, AL (162000 frames)
expt(90).mouse = 'i1395';
expt(90).date = '240607';
expt(90).img_loc  = {'AL'};
expt(90).z = -200;
expt(90).img_strct  = {'cells'};
expt(90).driver = {'SLC'};
expt(90).indicator = {'tg';'GCaMP6s'};
expt(90).coFolder = {'002'};
expt(90).coTime = {'1407'};
expt(90).SF = 0.05;
expt(90).TF = 2;
expt(90).saveLoc = 'sara';

%% i1394 240610, LM (162000 frames)
expt(91).mouse = 'i1394';
expt(91).date = '240610';
expt(91).img_loc  = {'LM'};
expt(91).z = -200;
expt(91).img_strct  = {'cells'};
expt(91).driver = {'SLC'};
expt(91).indicator = {'tg';'GCaMP6s'};
expt(91).coFolder = {'002'};
expt(91).coTime = {'1243'};
expt(91).SF = 0.05;
expt(91).TF = 2;
expt(91).saveLoc = 'sara';

%% i1394 240611, AL (162000 frames)
expt(92).mouse = 'i1394';
expt(92).date = '240611';
expt(92).img_loc  = {'AL'};
expt(92).z = -200;
expt(92).img_strct  = {'cells'};
expt(92).driver = {'SLC'};
expt(92).indicator = {'tg';'GCaMP6s'};
expt(92).coFolder = {'002'};
expt(92).coTime = {'1106'};
expt(92).SF = 0.05;
expt(92).TF = 2;
expt(92).saveLoc = 'sara';

%% i1397 240611, AL (162000 frames)
expt(93).mouse = 'i1397';
expt(93).date = '240611';
expt(93).img_loc  = {'AL'};
expt(93).z = -200;
expt(93).img_strct  = {'cells'};
expt(93).driver = {'SLC'};
expt(93).indicator = {'tg';'GCaMP6s'};
expt(93).coFolder = {'002'};
expt(93).coTime = {'1445'};
expt(93).SF = 0.05;
expt(93).TF = 2;
expt(93).saveLoc = 'sara';

%% i1393 240612, AL (162000 frames)
expt(94).mouse = 'i1393';
expt(94).date = '240612';
expt(94).img_loc  = {'AL'};
expt(94).z = -200;
expt(94).img_strct  = {'cells'};
expt(94).driver = {'SLC'};
expt(94).indicator = {'tg';'GCaMP6s'};
expt(94).coFolder = {'002'};
expt(94).coTime = {'1055'};
expt(94).SF = 0.05;
expt(94).TF = 2;
expt(94).saveLoc = 'sara';

%% i1397 240612, PM (162000 frames)
expt(95).mouse = 'i1397';
expt(95).date = '240612';
expt(95).img_loc  = {'PM'};
expt(95).z = -200;
expt(95).img_strct  = {'cells'};
expt(95).driver = {'SLC'};
expt(95).indicator = {'tg';'GCaMP6s'};
expt(95).coFolder = {'002'};
expt(95).coTime = {'1438'};
expt(95).SF = 0.05;
expt(95).TF = 2;
expt(95).saveLoc = 'sara';

%% i1394 240613, PM (162000 frames)
expt(96).mouse = 'i1394';
expt(96).date = '240613';
expt(96).img_loc  = {'PM'};
expt(96).z = -200;
expt(96).img_strct  = {'cells'};
expt(96).driver = {'SLC'};
expt(96).indicator = {'tg';'GCaMP6s'};
expt(96).coFolder = {'002'};
expt(96).coTime = {'1422'};
expt(96).SF = 0.05;
expt(96).TF = 2;
expt(96).saveLoc = 'sara';

%% i1395 240614, PM (162000 frames)
expt(97).mouse = 'i1395';
expt(97).date = '240614';
expt(97).img_loc  = {'PM'};
expt(97).z = -200;
expt(97).img_strct  = {'cells'};
expt(97).driver = {'SLC'};
expt(97).indicator = {'tg';'GCaMP6s'};
expt(97).coFolder = {'002'};
expt(97).coTime = {'1008'};
expt(97).SF = 0.05;
expt(97).TF = 2;
expt(97).saveLoc = 'sara';

%% i1396 240614, PM (162000 frames)
expt(98).mouse = 'i1396';
expt(98).date = '240614';
expt(98).img_loc  = {'PM'};
expt(98).z = -200;
expt(98).img_strct  = {'cells'};
expt(98).driver = {'SLC'};
expt(98).indicator = {'tg';'GCaMP8m'};
expt(98).coFolder = {'002'};
expt(98).coTime = {'1408'};
expt(98).SF = 0.05;
expt(98).TF = 2;
expt(98).saveLoc = 'sara';

%% i1397 240617, AL (162000 frames)
expt(99).mouse = 'i1397';
expt(99).date = '240617';
expt(99).img_loc  = {'AL'};
expt(99).z = -200;
expt(99).img_strct  = {'cells'};
expt(99).driver = {'SLC'};
expt(99).indicator = {'tg';'GCaMP6s'};
expt(99).coFolder = {'002'};
expt(99).coTime = {'1021'};
expt(99).SF = 0.05;
expt(99).TF = 2;
expt(99).saveLoc = 'sara';

%% i1394 240617, PM (162000 frames)
expt(100).mouse = 'i1394';
expt(100).date = '240617';
expt(100).img_loc  = {'PM'};
expt(100).z = -200;
expt(100).img_strct  = {'cells'};
expt(100).driver = {'SLC'};
expt(100).indicator = {'tg';'GCaMP6s'};
expt(100).coFolder = {'002'};
expt(100).coTime = {'1427'};
expt(100).SF = 0.05;
expt(100).TF = 2;
expt(100).saveLoc = 'sara';

%% i1397 240618, PM (162000 frames)
expt(101).mouse = 'i1397';
expt(101).date = '240618';
expt(101).img_loc  = {'PM'};
expt(101).z = -200;
expt(101).img_strct  = {'cells'};
expt(101).driver = {'SLC'};
expt(101).indicator = {'tg';'GCaMP6s'};
expt(101).coFolder = {'002'};
expt(101).coTime = {'1421'};
expt(101).SF = 0.05;
expt(101).TF = 2;
expt(101).saveLoc = 'sara';

%% i1396 240624, PM (162000 frames)
expt(102).mouse = 'i1396';
expt(102).date = '240624';
expt(102).img_loc  = {'PM'};
expt(102).z = -200;
expt(102).img_strct  = {'cells'};
expt(102).driver = {'SLC'};
expt(102).indicator = {'tg';'GCaMP6s'};
expt(102).coFolder = {'002'};
expt(102).coTime = {'0948'};
expt(102).SF = 0.05;
expt(102).TF = 2;
expt(102).saveLoc = 'sara';

