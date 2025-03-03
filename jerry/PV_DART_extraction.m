clear all; clear global; close all;
clc

ds = 'DART_expt_info'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories
eval(ds);

day_id = [24 28 31 34 37 40 42 44 46]; % concat days, enter one post-DART day id for each mouse

fnroot = fullfile(rc.achAnalysis,'PV_YM90K','summary_analyses');
fn_epi = fullfile(fnroot,'epileptiform');

mkdir(fn_epi);

nOn = 30;
nOff = 60;
nTot = nOn+nOff;
ntrials = 960;
%% extract off_np_tc data

% read data for every experiment session and put into cells (1st column
% data, 2nd column mouse ID)
prepped_np_TCs = cell(size(day_id,2),2);


for idx = 1:length(day_id)
    %find all of the sessions this current mouse had
    this_day = day_id(idx);
    this_mouse = expt(this_day).mouse;
    this_mouse_num = str2double(this_mouse(2:end));
    these_sesh = query_expt(this_mouse_num); %query_expt is a custom function by TH, see "help query_expt"
    mouse_data_cell = cell(size(these_sesh,1),5); %first column stores raw F TC, second stores mean-subtracted, third stores nan-padded (full trial length), fourth is full np tc; fifth is mean-subtracted full np tc
    for j = 1:length(these_sesh) % extract data from every session this mouse had
        this_sesh = these_sesh(j);
        data_path_cell = fullfile(rc.achAnalysis,'PV_YM90K',this_mouse, expt(this_sesh).date, expt(this_sesh).contrastxori_runs, 'TCs.mat');
        data_path = data_path_cell{1};
        load(data_path);
        avg_np_tc = mean(np_tc,2);
        avg_np_tc_trial = reshape(avg_np_tc,nTot,[]);
        avg_np_tc_f0 = avg_np_tc_trial(nOff/2+1:nOff,:);
        avg_np_tc_nanpad_trial = nan(nTot,ntrials);
        for itrial = 1:ntrials
            avg_np_tc_nanpad_trial(nOff/2+1:nOff,itrial) = avg_np_tc_f0(:,itrial);
        end
        avg_np_tc_nanpad = reshape(avg_np_tc_nanpad_trial,[],1);
        avg_np_tc_off = reshape(avg_np_tc_f0,[],1);
        mean_sub_tc = avg_np_tc_off - mean(avg_np_tc_off);
        full_avg_tc = mean(np_tc,2);
        meanSub_full_tc = full_avg_tc - mean(full_avg_tc);

        % save np tc during the last 30 off frames of every trial
        mouse_data_cell{j,1} = avg_np_tc_off;
        mouse_data_cell{j,2} = mean_sub_tc;
        mouse_data_cell{j,3} = avg_np_tc_nanpad;
        mouse_data_cell{j,4} = np_tc;
        mouse_data_cell{j,5} = meanSub_full_tc;
    end
    prepped_np_TCs{idx,1} = mouse_data_cell; % save off frames np tc for all mice into the big cell array
    prepped_np_TCs{idx,2} = this_mouse;
end

% d=string(datetime('today'));
save(fullfile(fn_epi,'prepped_np_TCs'),'prepped_np_TCs','day_id','-v7.3');




