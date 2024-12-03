clear all;
%paths
ds = 'DART_expt_info'; %dataset info
dataStructLabels = {'contrastxori'};
rc =  behavConstsDART; %directories178
eval(ds);

experimentFolder = 'PV_CMPDA';

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
        
    CD = fullfile(rc.achData, mouse, date, run);
    data_out = fullfile(rc.achAnalysis,experimentFolder, mouse, date,run);
    cd(CD);
    fn = [run '_000_000_eye.mat'];
    
    %load data
    data_temp = load(fn);
    data_temp = squeeze(data_temp.data);
    
    %crop frames to match mworks data
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time '.mat'];
    load(fName);
    %nFrames = input.counterValues{end}(end);
    nFrames = 86400;
    %nTrials = size(input.counterValues,2);
    nTrials = 960;
    data = data_temp(:,:,1:nFrames);      % the raw images...
    
    % Crop image to isolate pupil 
    %(bright spots can be mistaken for pupil)
    if day == 1
        [data_crop rect] = cropEyeData(data);
    else
        [data_crop rect] = cropEyeData(data,rect);
    end
    
    % measure pupil position/diameter
    rad_range = [3 30]; %adjust to expected range of pupil size (if low end is too small then may find noisy bright stuff)
    Eye_data = extractEyeData(data_crop,rad_range);
    %if pupil not found reliably, adjust the image cropping or the rad_range
    
    % align to stimulus presentation
    [rad centroid] = alignEyeData(Eye_data,input);
    
    % wheel data
    wheel_data = wheelSpeedCalc(input,32,'orange');
    wheel_trial = mean(reshape(wheel_data,[input.nScansOff+input.nScansOn nTrials]),1);
    
    save(fullfile(data_out,'pupil.mat'),'rect','rad_range','Eye_data','rad','centroid','wheel_trial')

end

