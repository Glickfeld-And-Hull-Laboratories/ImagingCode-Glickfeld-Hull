clear all; clear global; close all;

% Get dataset and experiment parameters
prompt = 'Enter ds (e.g., DART_V1_YM90K_Celine): ';
ds = input(prompt, 's');
clear prompt
rc = behavConstsDART; % directories
eval(ds);
day_id = input('Enter day id '); 
pre_day = expt(day_id).multiday_matchdays;
nd = 2; % Number of days hardcoded to = 2
mouse = expt(day_id).mouse;
experimentFolder = expt(day_id).exptType;
allDays = [day_id, pre_day];

% Process eye tracking data for each day
for day = 1:2
    iexp = allDays(day);
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    run = expt(iexp).contrastxori_runs{1};
    time = expt(iexp).contrastxori_time{1};
    
    % Set up paths
    CD = fullfile(rc.data, mouse, date, run);
    data_out = fullfile(rc.analysis, experimentFolder, mouse, date, run);
    cd(CD);
    fn = [run '_000_000_eye.mat'];
    
    % Load raw eye data
    data_temp = load(fn);
    data_temp = squeeze(data_temp.data);
    
    % Load behavioral data and determine frame timing
    infofName = fullfile(rc.data, expt(iexp).mouse, expt(iexp).date, run, [run '_000_000.mat']);
    inputfName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time '.mat'];
    load(inputfName);
    
    % Get stimulus timing information
    if isfile(infofName)
        load(infofName);
        [stimOns, stimOffs] = photoFrameFinder_Sanworks(info.frame);
    else
        warningMessage = sprintf('Warning: info.frame does not exist:\n%s\n Using mWks StimOn', infofName);
        uiwait(msgbox(warningMessage));
    end
    
    % Crop data to match behavioral recording length
    nFrames = input.counterValues{end}(end);
    nTrials = length(stimOns);
    data = data_temp(:, :, 1:nFrames); % the raw images
    
    % Crop image to isolate pupil (bright spots can be mistaken for pupil)
    if day == 1
        [data_crop, rect] = cropEyeData(data);
    else
        [data_crop, rect] = cropEyeData(data, rect);
    end
    
    % Extract pupil position and diameter
    rad_range = [3 15]; % adjust to expected range of pupil size (if low end is too small then may find noisy bright stuff)
    Eye_data = extractEyeData(data_crop, rad_range);
    % if pupil not found reliably, adjust the image cropping or the rad_range
    
    % Align eye data to stimulus presentation
    [rad, centroid] = alignEyeData_PD(Eye_data, input, stimOns);
    
    % Process wheel data
    wheel_data = wheelSpeedCalc(input, 32, 'orange');
    wheel_trial = mean(reshape(wheel_data, [input.nScansOff + input.nScansOn, nTrials]), 1);
    
    % Save processed data
    save(fullfile(data_out, 'pupil.mat'), 'rect', 'rad_range', 'Eye_data', 'rad', 'centroid', 'wheel_trial')
end

clear all;
