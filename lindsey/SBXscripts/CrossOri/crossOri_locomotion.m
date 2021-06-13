clc; clear all; close all;
doRedChannel = 0;
ds{1} = 'CrossOriRandPhase_ExptList';
ds{2} = 'CrossOriRandPhase_lowSF_ExptList';
%ds = 'CrossOriRandDirRandPhase_ExptList';
ds{3} = 'CrossOriRandDirTwoPhase_ExptList';
ds{4} = 'CrossOriSingleStimAdapt_ExptList';
rc = behavConstsAV;
doWheel = [];
fractLoc = [];
for i = 1:length(ds)
    clear expt
eval(ds{i})
nexp = length(expt);

start = 1;
for iexp = 1:nexp
    x = 0;
    if isfield(expt,'SF')
        if expt(iexp).SF == 0.05 & strcmp(expt(iexp).img_loc(1),'V1') & strcmp(expt(iexp).driver,'SLC')
            x = 1;
        end
    elseif strcmp(expt(iexp).img_loc(1),'V1') & strcmp(expt(iexp).driver,'SLC')
        x = 1;
    end
    
    if x == 1

    frame_rate = 15;

    %%
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc;
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);

    LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
    %LG_base = '\\CRASH.dhe.duke.edu\data\home\lindsey';

    fprintf([mouse ' ' date '\n'])

    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
%     load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
%     
%     wheel_speed = wheelSpeedCalc(input,32,'red');
%     nTrials = length(input.tStimOneGratingDirectionDeg);
%     cStimOn = celleqel2mat_padded(input.cStimOneOn);
%     cStimOff = celleqel2mat_padded(input.cStimOneOff);
% %     nTrials = length(input.tTestStimGratingDirectionDeg);
% %     cStimOn = celleqel2mat_padded(input.cTestOn);
% %     cStimOff = celleqel2mat_padded(input.cTestOff);
%     nFrames = length(wheel_speed);
%     wheel_trial = nan(1,nTrials);
%     for itrial = 1:nTrials
%     	wheel_trial(1,itrial) = mean(wheel_speed(1,cStimOn(itrial):cStimOff(itrial)),2);
%     end
%     
    
    %save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_wheelData.mat']),'wheel_speed','wheel_trial')
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_wheelData.mat']))
    doWheel = [doWheel length(find(wheel_trial~=0))];
    fractLoc = [fractLoc length(find(wheel_trial>2))./length(wheel_trial)];
%     subplot(5,4,start)
%     plot(wheel_trial)
%     ylim([-10 30])
%     title([mouse ' ' date '- ' num2str(length(find(wheel_trial>2))) ' trials'])
%     start = start+1;   
    end
end
end
nSession = length(find(doWheel));
fractLoc_avg = [mean(fractLoc(find(doWheel)),2) std(fractLoc(find(doWheel)),[],2)./sqrt(length(find(doWheel)))];
fractLoc_med = median(fractLoc(find(doWheel)),2);
fractLoc_min_max = [min(fractLoc(find(doWheel)),[],2) max(fractLoc(find(doWheel)),[],2)];

% fnout = fullfile(LG_base, 'Analysis\2P\CrossOri\CrossOri_Figures');
% print(fullfile(fnout,'CrossOriRandPhase2_locomotion.pdf'),'-dpdf')