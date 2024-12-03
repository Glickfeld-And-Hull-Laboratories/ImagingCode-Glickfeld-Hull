clear all; clear global; close all;
clc

ds = 'DART_expt_info'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories
eval(ds);

% day_id = [24 29 32 35 38 40 42 44]; % concat days, enter one post-DART day id for each mouse

fnroot = fullfile(rc.achAnalysis,'PV_CMPDA','summary_analyses');
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
t_since_drug_mat = strings(nMice,3);

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
        t_since_drug_mat(mouse,sess_idx) = t_since_drug_hrs;
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
    % saveas(gcf,fullfile(fn_epi,'plots',[mouse_id '.pdf']));
end

save(fullfile(fn_epi,'epi_data.mat'),"std_mat","perc_mat","t_since_drug_mat");
clear all_raw_f ylim_max ylim_min std_mat perc_mat curr_sesh off_np_tc
%% plot epilepsy data as modulation index

cd(fn_epi);
load('epi_data.mat');
load('DART_expt_info.mat');
% for dart_expt_info, 2nd column is dart dosage, 3rd column is virus titer
% DART:
% 1 = YM90K-PEG 10:1 alx647-COOH
% 2 = 1:10:1 YM90K-DART:blank-DART:alx-DART
% 3 = 4:6:1 YM90K-DART:blank-DART:alx-DART
% 4 = 5:5:1 YM90K-DART:blank-DART:alx-DART
% 5 = YM90K-DART 10:1 alx647-DART
% virus:
% 1 = GCaMP6s 1:1 MB HTP
% 2 = GCaMP8f 1:1 NES HTP
% 3 = 2:1:1 GCaMP8f:NES HTP:aCSF

if exist('off_np_tc','var') == 1
    mouse_list = convertCharsToStrings(off_np_tc(:,2)); 
    clear off_np_tc
else
    load('off_np_tc.mat');
    mouse_list = convertCharsToStrings(off_np_tc(:,2));
    clear off_np_tc
end

nMice = size(mouse_list,1);
TPs_mat = t_since_drug_mat(:,2:3);
std_index_mat = nan(nMice,3);
frac_index_mat = nan(nMice,3);

% Fig 1 frac norm
figure;
set(gcf, 'Name', 'EpiPlot_frac_norm');
title('Modulation of Epileptic Activity, frac, normalized');
xlim([0 26]);
ylim([0 3.7]);
xlabel('Time Since DART Infusion (hrs)');
ylabel('Fraction Modulation Index (AU)');
hold on

for mouse = 1:nMice
    if sum(~isnan(std_mat(mouse,:))) == 2
        nTP = 2;
    elseif sum(~isnan(std_mat(mouse,:))) == 3
        nTP = 3;
    else
        error('Incorrect number of timepoints, check data matrix dimensions')
    end
    % calculate normalized fraction value
    for tp = 1:nTP
        std_index_mat(mouse,tp) = std_mat(mouse,tp)/std_mat(mouse,1);
        frac_index_mat(mouse,tp) = perc_mat(mouse,tp)/perc_mat(mouse,1);
    end

    % color code by dart dosage 
    if expt_info{mouse,2} == 1
        this_col = "#D95319"; % mute orange-ish
    elseif expt_info{mouse,2} == 2
        this_col = "#40E0D0"; % turquoise
    elseif expt_info{mouse,2} == 3
        this_col = "#89CFF0"; % baby blue
    elseif expt_info{mouse,2} == 4
        this_col = "#7393B3"; % blue gray
    elseif expt_info{mouse,2} == 5
        this_col = "#0047AB"; % cobalt blue
    else
        error('Unrecognized DART dosage value. Check metadata values.')
    end

    % line style code by virus titer/type
    if expt_info{mouse,3} == 1
        this_marker = "o"; 
    elseif expt_info{mouse,3} == 2
        this_marker = "^"; 
    elseif expt_info{mouse,3} == 3
        this_marker = "square"; 
    else
        error('Unrecognized viral titer. Check metadata values.');
    end
    plot(str2double(t_since_drug_mat(mouse,:)),frac_index_mat(mouse,:),"Color",this_col,"Marker",this_marker,"LineWidth",1.5);
end

legend
hold off

% fig 2 frac un-norm

figure;
set(gcf, 'Name', 'EpiPlot_frac_raw');
title('Modulation of Epileptic Activity, frac, raw value');
xlim([0 26]);
ylim([0 0.027]);
xlabel('Time Since DART Infusion (hrs)');
ylabel('Fraction of Frames');
hold on

for mouse = 1:nMice
    if sum(~isnan(std_mat(mouse,:))) == 2
        nTP = 2;
    elseif sum(~isnan(std_mat(mouse,:))) == 3
        nTP = 3;
    else
        error('Incorrect number of timepoints, check data matrix dimensions')
    end

    % color code by dart dosage 
    if expt_info{mouse,2} == 1
        this_col = "#D95319"; % mute orange-ish
    elseif expt_info{mouse,2} == 2
        this_col = "#40E0D0"; % turquoise
    elseif expt_info{mouse,2} == 3
        this_col = "#89CFF0"; % baby blue
    elseif expt_info{mouse,2} == 4
        this_col = "#7393B3"; % blue gray
    elseif expt_info{mouse,2} == 5
        this_col = "#0047AB"; % cobalt blue
    else
        error('Unrecognized DART dosage value. Check metadata values.')
    end

    % line style code by virus titer/type
    if expt_info{mouse,3} == 1
        this_marker = "o"; 
    elseif expt_info{mouse,3} == 2
        this_marker = "^"; 
    elseif expt_info{mouse,3} == 3
        this_marker = "square"; 
    else
        error('Unrecognized viral titer. Check metadata values.');
    end
    plot(str2double(t_since_drug_mat(mouse,:)),perc_mat(mouse,:),"Color",this_col,"Marker",this_marker,"LineWidth",1.5);
end
legend
hold off

% Fig 3 std norm
figure;
set(gcf, 'Name', 'EpiPlot_std_norm');
title('Modulation of Epileptic Activity, std, normalized');
xlim([0 26]);
% ylim([0 3.7]);
xlabel('Time Since DART Infusion (hrs)');
ylabel('Std Modulation Index (AU)');
hold on

for mouse = 1:nMice
    if sum(~isnan(std_mat(mouse,:))) == 2
        nTP = 2;
    elseif sum(~isnan(std_mat(mouse,:))) == 3
        nTP = 3;
    else
        error('Incorrect number of timepoints, check data matrix dimensions')
    end
    % calculate normalized fraction value
    for tp = 1:nTP
        std_index_mat(mouse,tp) = std_mat(mouse,tp)/std_mat(mouse,1);
        frac_index_mat(mouse,tp) = perc_mat(mouse,tp)/perc_mat(mouse,1);
    end

    % color code by dart dosage 
    if expt_info{mouse,2} == 1
        this_col = "#D95319"; % mute orange-ish
    elseif expt_info{mouse,2} == 2
        this_col = "#40E0D0"; % turquoise
    elseif expt_info{mouse,2} == 3
        this_col = "#89CFF0"; % baby blue
    elseif expt_info{mouse,2} == 4
        this_col = "#7393B3"; % blue gray
    elseif expt_info{mouse,2} == 5
        this_col = "#0047AB"; % cobalt blue
    else
        error('Unrecognized DART dosage value. Check metadata values.')
    end

    % line style code by virus titer/type
    if expt_info{mouse,3} == 1
        this_marker = "o"; 
    elseif expt_info{mouse,3} == 2
        this_marker = "^"; 
    elseif expt_info{mouse,3} == 3
        this_marker = "square"; 
    else
        error('Unrecognized viral titer. Check metadata values.');
    end
    plot(str2double(t_since_drug_mat(mouse,:)),std_index_mat(mouse,:),"Color",this_col,"Marker",this_marker,"LineWidth",1.5);
end

legend
hold off

% fig 4 std un-norm

figure;
set(gcf, 'Name', 'EpiPlot_std_raw');
title('Modulation of Epileptic Activity, std, raw value');
xlim([0 26]);
% ylim([0 0.027]);
xlabel('Time Since DART Infusion (hrs)');
ylabel('Standard Deviation of F');
hold on

for mouse = 1:nMice
    if sum(~isnan(std_mat(mouse,:))) == 2
        nTP = 2;
    elseif sum(~isnan(std_mat(mouse,:))) == 3
        nTP = 3;
    else
        error('Incorrect number of timepoints, check data matrix dimensions')
    end

    % color code by dart dosage 
    if expt_info{mouse,2} == 1
        this_col = "#D95319"; % mute orange-ish
    elseif expt_info{mouse,2} == 2
        this_col = "#40E0D0"; % turquoise
    elseif expt_info{mouse,2} == 3
        this_col = "#89CFF0"; % baby blue
    elseif expt_info{mouse,2} == 4
        this_col = "#7393B3"; % blue gray
    elseif expt_info{mouse,2} == 5
        this_col = "#0047AB"; % cobalt blue
    else
        error('Unrecognized DART dosage value. Check metadata values.')
    end

    % line style code by virus titer/type
    if expt_info{mouse,3} == 1
        this_marker = "o"; 
    elseif expt_info{mouse,3} == 2
        this_marker = "^"; 
    elseif expt_info{mouse,3} == 3
        this_marker = "square"; 
    else
        error('Unrecognized viral titer. Check metadata values.');
    end
    plot(str2double(t_since_drug_mat(mouse,:)),std_mat(mouse,:),"Color",this_col,"Marker",this_marker,"LineWidth",1.5);
end
legend
hold off

% figHandles = findall(0,'Type','figure');
% this_figname = get(figHandles,'Name');
% for i = 1:numel(figHandles)
%     t = this_figname{i,1};
%     saveas(figHandles(i),fullfile(fn_epi,'plots',[t '.pdf']));
% end
%
% or use "save_all_open_figs(directory_to_be_saved_to)"


%% 


