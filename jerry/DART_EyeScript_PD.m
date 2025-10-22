clear all; clear global; close all;
%paths
ds = 'DART_expt_info_jerry'; %dataset info
dataStructLabels = {'contrastxori'};
rc =  behavConstsDART; %directories
eval(ds);

experimentFolder = 'PV_YM90K';

day_id = input('Enter day id ');% alternative to run from command line.
pre_day = expt(day_id).multiday_matchdays;
nd=2; %hardcoding the number of days for now

mouse = expt(day_id).mouse;

allDays = [day_id,pre_day];
%%

for day = 1:2
    iexp = allDays(day);
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    run = expt(iexp).contrastxori_runs{1};
    time = expt(iexp).contrastxori_time{1};
        
    CD = fullfile(rc.data, mouse, date, run);
    data_out = fullfile(rc.analysis,experimentFolder, mouse, date,run);
    cd(CD);
    fn = [run '_000_000_eye.mat'];
    %load data
    data_temp = load(fn);
    data_temp = squeeze(data_temp.data);
    
    %crop frames to match mworks data
    infofName = fullfile(rc.achData,expt(iexp).mouse,expt(iexp).date,run,[run '_000_000.mat']);
    inputfName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time '.mat'];
    load(inputfName);
    if isfile(infofName)
        load(infofName);
        [stimOns,stimOffs] = photoFrameFinder_Sanworks(info.frame);
    else
        warningMessage = sprintf('Warning: info.frame does not exist:\n%s\n Using mWks StimOn', infofName);
        uiwait(msgbox(warningMessage));
    end
    
    %load('G:\home\ACh\Analysis\2p_analysis\i3321\241125\i3321_241125_runs-004_correctedTiming.mat');
    nFrames = input.counterValues{end}(end);
    %nFrames = 86400;
    nTrials = length(stimOns);
    %nTrials = 960;
    data = data_temp(:,:,1:nFrames);      % the raw images...
    
    % Crop image to isolate pupil 
    %(bright spots can be mistaken for pupil)
    if day == 1
        [data_crop rect] = cropEyeData(data);
    else
        [data_crop rect] = cropEyeData(data,rect);
    end
    
    % measure pupil position/diameter
    rad_range = [3 40]; %adjust to expected range of pupil size (if low end is too small then may find noisy bright stuff)
    Eye_data = extractEyeData_jerry(data_crop,rad_range);
    %if pupil not found reliably, adjust the image cropping or the rad_range
    
    % align to stimulus presentation
    [rad centroid] = alignEyeData_PD(Eye_data,input,stimOns);
    
    % wheel data
    wheel_data = wheelSpeedCalc(input,32,'orange');
    wheel_trial = mean(reshape(wheel_data,[input.nScansOff+input.nScansOn nTrials]),1);
    
    save(fullfile(data_out,'pupil.mat'),'rect','rad_range','Eye_data','rad','centroid','wheel_trial')

end

