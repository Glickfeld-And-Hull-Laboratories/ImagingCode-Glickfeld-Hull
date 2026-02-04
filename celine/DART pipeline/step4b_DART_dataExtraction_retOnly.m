%standaloe script to get match-aligned retinotopy data without redoing the
%extraction
clear all; clear global; clc
prompt = 'Enter name of instructions file: ';
instr = input(prompt, 's');
clear prompt
run(instr);

ds=instructions.ds;
run(ds); 

dataStructLabels = {'contrastxori'};
rc = behavConstsDART; % directories

% Input validation
if ~exist('expt', 'var')
    error('Dataset %s not found', ds);
end

day_id = str2double(instructions.session);


if length(expt) < day_id
    error('Day_id %d not valid for this dataset', day_id);
else
    match_day = expt(day_id).multiday_matchdays;

end

nd = 2; % hardcoding the number of days for now
mouse = expt(day_id).mouse;
experimentFolder = expt(day_id).exptType;
fnout = fullfile(rc.analysis, mouse);

% Set up file naming based on drug condition
if expt(day_id).multiday_timesincedrug_hours > 0
    dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end
fn_multi = fullfile(rc.analysis, experimentFolder, mouse, ['multiday_' dart_str]);

% Determine which session was used as reference for cell matching
x = instructions.refDay;
switch x
    case '1'
        pre = 1;  % baseline session, used as reference, is in the 1st position
        post = 2;
        fprintf('Baseline used as reference\n');
         allDays = [match_day,day_id];
        fprintf('Analyzing sessions: %s\n', num2str(allDays));
    case '2'
        pre = 2;
        post = 1;  % post-DART session, used as reference, is in the 1st position
        fprintf('Post-DART used as reference\n');
        allDays = [day_id, match_day];
        fprintf('Analyzing sessions: %s\n', num2str(allDays));
end
clear x instr

% Load the matched cell data
cd(fn_multi)
load(fullfile(fn_multi, 'timecourses.mat'))
load(fullfile(fn_multi, 'multiday_alignment.mat'))
load(fullfile(fn_multi, 'input.mat'))
load(fullfile(fn_multi,'cell_analysis.mat'))

% Rename input to inputStructure to avoid conflict with MATLAB's input function
inputStructure = input;
clear input
%%

[ret_npSub_tc_matched, ret_distance_matched,resp_by_stim_matched,ret_dfof_trial_matched] = retinotopy_for_matched_data(nd, ...
    allDays, expt, mouse, fov_avg, masks, fitGeoTAf, ...
    instructions, inputStructure,match_ind,false);

trialIndSourceUsed=instructions.tIdxSource;

%make "keep" subsets for the matched retino data 
ret_npSub_tc_keep = cell(1,nd);
ret_distance_keep = cell(1,nd);
resp_by_stim_keep=cell(1,nd);
ret_dfof_trial_keep=cell(1,nd);

for id = 1:nd
    ret_npSub_tc_keep{id}=ret_npSub_tc_matched{id}(:,keep_cells);
    ret_distance_keep{id}=ret_distance_matched{id}(:,keep_cells);
    resp_by_stim_keep{id}=resp_by_stim_matched{id}(:,:,keep_cells);
    ret_dfof_trial_keep{id}=ret_dfof_trial_matched{id}(:,keep_cells,:);
end

save(fullfile(fn_multi, 'retino_aligned.mat'), 'ret_npSub_tc_matched', ...
    'ret_distance_matched','resp_by_stim_matched','ret_dfof_trial_matched',...
    'ret_npSub_tc_keep','ret_distance_keep','resp_by_stim_keep','ret_dfof_trial_keep','trialIndSourceUsed');

