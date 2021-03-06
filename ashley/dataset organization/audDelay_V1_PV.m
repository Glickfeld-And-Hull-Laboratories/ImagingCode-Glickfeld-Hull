params.frameRate = 30;
params.stimOnTime = 100;
params.nBaselineMs = 1000;
params.nStimMs_visAlign = 3000;
params.nStimMs_audAlign = 9000;
params.nFramesVisDelay = 4;
params.motionCutoff = 0.2;

%% 905 180321
expt(1).SubNum = '905';
expt(1).mouse = '905';
expt(1).date = '180321';
expt(1).img_loc  = {'V1';'L2/3'};
expt(1).z = -284;
expt(1).img_strct  = {'cells'};
expt(1).indicator = {'virus';'PV-GCaMP6s'};
expt(1).time_mat = ['1019';'1051'];
expt(1).runs = ['002';'003'];
expt(1).nrun = size(expt(1).runs,1);
expt(1).frame_rate = 30;
expt(1).rettuning = {'001';'1007'};
expt(1).stimOnMs = 100;
expt(1).motionThreshold = 0.05;
expt(1).areaBorders = 0;
expt(1).regImgStartFrame = 68284;