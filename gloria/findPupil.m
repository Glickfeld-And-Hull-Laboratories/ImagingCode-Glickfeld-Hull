clear all
clear all global
close all

%Set paths
data_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\ACh\Aging';
fnout = [data_pn, '\Gloria\'];

%Set variables
mouse = 'i2207';
date = '230306';
run = '002';
time = '1522';
DOB = '2021-02-25';
record_date = '2023-03-06'; 
datemouse = [date '_' mouse ];
datemouserun = [date '_' mouse '_' run];

%% 
%load data
CD = [data_pn '\data\2p\' mouse '\' date '\' run];
cd(CD);
fn = [run '_000_000_eye.mat'];
data_temp = load(fn);
data_temp = squeeze(data_temp.data);

%crop frames to match mworks data
fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data\data-' mouse '-' date '-' time '.mat'];
load(fName);
nFrames = input.counterValues{end}(end);
data = data_temp(:,:,1:nFrames); % the raw images...

%% 
CD = [data_pn '\Gloria\eyeParameters'];
cd(CD);
eyeParameters = [datemouserun '_eyeParameters.mat'];

if isfile(eyeParameters)
    load(eyeParameters)
    fprintf(['Loaded ' eyeParameters '\r\n'])
else
    mkdir(fullfile(fnout, mouse, datemouse, datemouserun)) %creates new folder
    % Crop image to isolate pupil (bright spots can be mistaken for pupil)
    [data_crop rect] = cropEyeData(data);
    % measure pupil position/diameter
    rad_range = [0 0]; %adjust to expected range of pupil size (if low end is too small then may find noisy bright stuff)
end
%% 
Eye_data = extractEyeData_Gloria(data_crop,rad_range,fnout, mouse, datemouse, datemouserun);
%if pupil not found reliably, adjust the image cropping or the rad_range

% align to stimulus presentation
[rad centroid] = alignEyeData_Gloria(Eye_data,input);
print(fullfile(fnout, mouse, datemouse, datemouserun, [datemouserun '_eyeCentroid.pdf']),'-dpdf','-fillpage')

%% Cleaning data...

nan_ind_TC = Eye_data.badFrames;
Eye_data.Area(Eye_data.badFrames) = NaN; %makes bad frames into NaN values
pupil_rad = sqrt(Eye_data.Area./pi); %calculate r from Area

nan_ind_trials = find(isnan(rad.stim)); %finds the trial numbers with NaN

%% Incorporating wheel data

nOn = input.nScansOn; % # frames on
nOff = input.nScansOff; % # frames off
ntrials = size(input.tGratingDirectionDeg,2);

[wheel_speed] = wheelSpeedCalc(input,32,'orange'); % arguments: the name of your mWorks input structure, the number of clicks in a rotation of the encoderm (32), and which wheel used (orange, purple, or red)

wheel_tc = zeros(nOn+nOff, ntrials); % segregates wheel speed time course into trials
for iTrial = 1:ntrials
    wheel_tc(:,iTrial) = wheel_speed(1+((iTrial-1).*(nOn+nOff)):iTrial.*(nOn+nOff));
end
wheel_trial_avg = mean(wheel_tc(nOff:nOn+nOff,:),1); % average wheel speed over each trial
RIx = wheel_trial_avg>2.0; % threshold to determine which trials are running (0 = false, 1 = true)
run_ind = find(RIx); % indices of all running trials
notrun_ind = find(~RIx); % indices of all not running trials

%% Export data

save(fullfile(fnout, "eyeTCs", [datemouserun '_eyeTCs.mat']), 'mouse', 'date', 'DOB', 'record_date', 'datemouserun', 'run', 'time', 'input','nan_ind_TC', 'pupil_rad', 'wheel_speed') 
save(fullfile(fnout, "eyeTrials", [datemouserun '_eyeTrials.mat']), 'mouse', 'DOB', 'record_date', 'datemouserun', 'wheel_trial_avg', 'rad', 'run_ind', 'notrun_ind', 'nan_ind_trials') % saves a separate file of the analyzed pupil and wheel data
save(fullfile(fnout, "eyeParameters", [datemouserun '_eyeParameters.mat']), 'mouse', 'date', 'DOB', 'record_date', 'run', 'time', 'rad_range', 'data_crop', 'rect', 'Eye_data','input') % saves all the analysis parameters in case it needs to be replicated

