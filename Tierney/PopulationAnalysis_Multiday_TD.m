%% Experiment info
clear all; clear global; close all
clc

%Path names
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
fn_analysis = fullfile(fn_base, 'home\Tierney\Analysis\2P'); % CHANGE THIS TO MULTIDAY IMAGING FOLDER
fn_multi = fullfile(fn_base, 'home\Tierney\Analysis\2P\MultidayAnalysis');
fn_pop = fullfile(fn_base, 'home\Tierney\Analysis\2P\PopulationAnalysis');
% fn_out = fullfile(fn_pop,['population_' time_str,'_',expt(day_id(2)).experiment]);

ds = 'ExperimentData_TD'; % dataset info 
eval(ds)

sess_list = [2, 4, 8, 11];% enter all the sessions you want to pool
nSess = length(sess_list);
% 
% nd=2%hard coding for two days per experimental session
% frame_rate = 15;
% 
% sess_title = string(sess_list(1));
% for iSess = 2:nSess
%     sess_title = strcat(sess_title,'_',string(sess_list(iSess)));
% end
% 
% fnout = fullfile(rc.achAnalysis,strcat('pooled_', sess_title));
% mkdir(fnout);

% day_id_list = 

% Baseline (1), post-MD (2), and recovery (3) days
day_id(2) = 4;
day_id(1) = expt(day_id(2)).matchday_baseline;
% day_id(3) = expt(day_id(2)).matchday_recovery;

nd = length(day_id);

%Specific experiment information
mouse = expt(day_id(1)).mouse;

for id = 1:nd
    date = expt(day_id(id)).date;
    runs = expt(day_id(id)).stimruns;
    nrun = length(runs);
    contra = strcmp(expt(day_id(id)).eye_str,'Contra'); % 1 is contra eye open; 0 is ipsi eye open
    
    run_str = catRunName(runs, nrun);
    datemouse = [date '_' mouse];
    datemouserun = [date '_' mouse '_' run_str];
   
    % Load respData and ODI
    fn_respData = fullfile(fn_multi, mouse, fn_match, [mouse, '_respData_multiday.mat']);
    fn_ODI = fullfile(fn_multi, mouse, fn_match, [mouse, '_ODI_multiday.mat']);
    fn_oriResp = fullfile(fn_analysis, datemouse, datemouserun, [datemouserun, '_oriResp' '.mat']);
    temp_respData = load(fn_respData);
    temp_ODI = load(fn_ODI);
    temp_oriResp = load(fn_oriResp);
      
    respData{id} = temp_respData;
    ODI{id} = temp_ODI;
    oriResp{id} = temp_oriResp;
    
    clear temp_respData
    clear temp_ODI
    clear temp_oriResp
end

% if expt(day_id(2)).multiday_time_days>0
%     time_str = [num2str(expt(day_id(2)).multiday_time_days) 'Days'];
% else
%     time_str = 'baseline';
% end


 %% Population averages for ODI
 
 % Get ODIs
 ODI_list = [];
 
 for id = 1:nd
     ODI_list = ODI{id}.ODI;
     mean_ODI(id) = nanmean(ODI_list);
 end
 
swarmchart(day_id(1:3),mean_ODI)
axis square
xlabel('Timepoint')
set(gca,'XTick',3:5,'XTickLabel',{'Pre','Post','Rec'})
ylabel('Mean ODI')
ylim([-1 1])
title('Mean ODI across days')

% MAKE EACH MOUSE A DIFFERENT COLOR

print(fullfile(fn_pop, [mouse '_meanODI.pdf']),'-dpdf','-bestfit')
