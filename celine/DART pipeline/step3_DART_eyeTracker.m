clear all; clear global; close all;

% Get dataset and experiment parameters
prompt = 'Enter name of instructions file: ';
instr = input(prompt, 's');
clear prompt
run(instr);
rc = behavConstsDART; % directories

ds=instructions.ds;
run(ds);

day_id = str2double(instructions.session);

if length(expt) < day_id
    error('Day_id %d not valid for this dataset', day_id);
else
    pre_day = expt(day_id).multiday_matchdays;
    allDays = [day_id, pre_day];
    fprintf('Analyzing sessions: %s\n', num2str(allDays));
end


nd = 2; % Number of days hardcoded to = 2
mouse = expt(day_id).mouse;
experimentFolder = expt(day_id).exptType;



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
    inputfName = fullfile(data_out,'input.mat');
    load(inputfName);

    nFrames = input.counterValues{end}(end);
switch instructions.tIdxSource
    case 'PD'
        stimOns=input.stimOns_photodiode;
    case 'MW'
        stimOns=cell2mat(input.stimOns_mwCounter);
        clear input_correct
    case 'cS'
        stimOns = cell2mat(input.cStimOn);
    otherwise
        error('No valid trial indexing source specificed in instr file. Use "PD" or "MW".');
end
    

    
    % Crop data to match behavioral recording length
    nTrials = length(stimOns);
    data = data_temp(:, :, 1:nFrames); % the raw images
    
    % Crop image to isolate pupil (bright spots can be mistaken for pupil)
    if day == 1
        [data_crop, rect] = cropEyeData(data);
    else
        [data_crop, rect] = cropEyeData(data, rect);
    end
    
    % Extract pupil position and diameter
    rad_range = [3 50]; % adjust to expected range of pupil size (if low end is too small then may find noisy bright stuff)
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
