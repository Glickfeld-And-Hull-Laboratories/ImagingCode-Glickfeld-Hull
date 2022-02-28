clear all; clear global; close all
clc

ds = 'ExperimentData_TD'; %dataset info
dataStructLabels = {'stimruns'};
eval(ds)

% day_id_list = 

day_id(2) = 4;
day_id(1) = expt(day_id(2)).multiday_matchdays;

nd = length(day_id);

%Specific experiment information
mouse = expt(day_id(1)).mouse;
date = expt(day_id(1)).date;
ImgFolder = expt(day_id(1)).stimruns;
eye_str = expt(day_id).eye_str;

frame_rate = 15;
contra = strcmp(eye_str,'Contra'); %blocks for eye stimulation: 1 is when contra is open; 0 is contra closed
nrun = size(ImgFolder,1);
run_str = catRunName(ImgFolder, nrun);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];

%Path names
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
td_fn = fullfile(fn_base, 'home\Tierney');
data_fn = fullfile(td_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data');
fnIn = fullfile(fn_base, 'home\Tierney\Analysis\2P');
fnPop = fullfile(fn_base, 'home\Tierney\Analysis\2P\PopulationAnalysis');

if expt(day_id(2)).multiday_time_days>0
    time_str = [num2str(expt(day_id(2)).multiday_time_days) 'Days'];
else
    time_str = 'baseline';
end

fn_pop = fullfile(fnIn,mouse,['population_' time_str,'_',expt(day_id(2)).experiment]);

%% Load previous analyses

 fName = fullfile(fnIn, datemouse, datemouserun, [datemouserun, '_respData' '.mat']);
 load(fName);
 
 %% Get population averages for ODI, contra/ipsi inputs, tuning, etc.