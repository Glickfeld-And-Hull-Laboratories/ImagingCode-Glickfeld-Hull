clear all; clear global; close all;
clc

ds = 'DART_expt_info'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories
eval(ds);

% day_id = [24 29 32 35 38 40 42 44]; % concat days, enter one post-DART day id for each mouse

fnroot = fullfile(rc.achAnalysis,'PV_YM90K','summary_analyses');
fn_epi = fullfile(fnroot,'epileptiform');

nOn = 30;
nOff = 60;
nTot = nOn+nOff;
ntrials = 960;

%% PV_DART off_np_tc epileptiform quantification

% off_np_tc: each row is a mouse. Each mouse cell contains a n x 2 cell array
% where n is the number of imaging sessions for the mouse, the 1st column
% is an averaged-across-cells timecourse made from a concatenation of the last 30 off-frames of
% every trial, and the 2nd column is the mean-subtracted timecourse of the
% 1st column. 

cd(fn_epi);
load('off_np_tc.mat');
nMice = size(off_np_tc,1);
perc_mat = NaN(nMice,3);
std_mat = NaN(nMice,3);

for mouse = 1:nMice %iterate through mouse
    mouse_data = off_np_tc{mouse,1};
    mouse_id = off_np_tc{mouse,2};
    this_sessions = query_expt(str2double(mouse_id(2:end)));

    all_raw_f = cat(1,mouse_data{:,1});
    ylim_max = 20 + max(all_raw_f);
    ylim_min = min(all_raw_f)-20;

    figure;
    sgtitle(mouse_id);
    for sess_idx = 1:size(mouse_data,1)
        curr_sesh = mouse_data{sess_idx,1}; % data for the current imaging session of the current mouse
        this_session_num = this_sessions(sess_idx);
        t_since_drug_hrs = expt(this_session_num).multiday_timesincedrug_hours;
        nplots = length(this_sessions);

        if t_since_drug_hrs == '0'
            this_tp_title = 'baseline';
        else
            this_tp_title = [t_since_drug_hrs ' hrs post-DART'];
        end

        mean_curr_sesh = mean(curr_sesh);
        std_curr_sesh = std(curr_sesh);
        frac_past = sum(curr_sesh > (mean_curr_sesh+3*std_curr_sesh))/length(curr_sesh);
        perc_mat(mouse, sess_idx) = frac_past;
        std_mat(mouse, sess_idx) = std_curr_sesh;

        subplot(nplots,1,sess_idx);
        plot(curr_sesh);
        ylim([ylim_min ylim_max]);
        title([this_tp_title ' std= ' num2str(std_curr_sesh) ' frac=' num2str(frac_past)]);
    end
    saveas(gcf,fullfile(fn_epi,'plots',[mouse_id '.pdf']));
end

