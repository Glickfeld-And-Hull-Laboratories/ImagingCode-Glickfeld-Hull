clear all; clear global; close all;
clc

ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories
eval(ds);

day_id = [169 177 183]; % concat days, enter one post-DART day id for each mouse

fnroot = fullfile(rc.achAnalysis,'SST_YM90K','summary_analyses');
fn_epi = fullfile(fnroot,'epileptiform');

mkdir(fn_epi);

nOn = 30;
nOff = 60;
nTot = nOn+nOff;
ntrials = 960;
%% extract off_np_tc data

% read data for every experiment session and put into cells (1st column
% data, 2nd column mouse ID)
prepped_np_TCs_SST = cell(size(day_id,2),2);


for idx = 1:length(day_id)
    %find all of the sessions this current mouse had
    this_day = day_id(idx);
    this_mouse = expt(this_day).mouse;
    this_mouse_num = str2double(this_mouse(2:end));
    these_sesh = query_expt_celine(this_mouse_num); %query_expt is a custom function by TH, see "help query_expt"
    mouse_data_cell = cell(size(these_sesh,1),4); %first column stores raw F TC, second stores mean-subtracted, third stores np_tc (everything), fourth is mean-subtracted full np tc
    for j = 1:length(these_sesh) % extract data from every session this mouse had
        this_sesh = these_sesh(j);
        data_path_cell = fullfile(rc.achAnalysis,'SST_YM90K',this_mouse, expt(this_sesh).date, expt(this_sesh).contrastxori_runs, 'TCs.mat');
        data_path = data_path_cell{1};
        load(data_path,'np_tc');
        avg_np_tc = mean(np_tc,2);
        avg_np_tc_trial = reshape(avg_np_tc,nTot,[]);
        avg_np_tc_f0 = avg_np_tc_trial(nOff/2+1:nOff,:);
        avg_np_tc_off = reshape(avg_np_tc_f0,[],1);
        mean_sub_tc = avg_np_tc_off - mean(avg_np_tc_off);
        full_avg_tc = mean(np_tc,2);
        meanSub_full_tc = full_avg_tc - mean(full_avg_tc);

        % save np tc during the last 30 off frames of every trial
        mouse_data_cell{j,1} = avg_np_tc_off;
        mouse_data_cell{j,2} = mean_sub_tc;
        mouse_data_cell{j,3} = np_tc;
        mouse_data_cell{j,4} = meanSub_full_tc;
    end
    prepped_np_TCs_SST{idx,1} = mouse_data_cell; % save off frames np tc for all mice into the big cell array
    prepped_np_TCs_SST{idx,2} = this_mouse;
end

% d=string(datetime('today'));wheelspeedcalc
save(fullfile(fn_epi,'prepped_np_TCs_SST'),'prepped_np_TCs_SST','day_id','-v7.3');


%% behavior file

runTrials_SST = cell(length(day_id),2);

for idx = 1:length(day_id)
    mouse = expt(day_id(idx)).mouse;
    runTrials_SST{idx,2} = mouse;
    this_mouse_num = str2double(mouse(2:end));
    % these_sesh = query_expt_celine(this_mouse_num);
    
    if mouse == "i2062"
        these_sesh = [167;169];
    elseif mouse == "i2066"
        these_sesh = [182;183];
    elseif mouse == "i2067"
        these_sesh = [176;177];
    else
        errow('Wrong Mouse');
    end

    nSesh = length(these_sesh);
    run_by_sesh = cell(nSesh,1);
    for sesh = 1:nSesh
        this_day = these_sesh(sesh);
        expDate = expt(this_day).date;
        expTime = expt(this_day).contrastxori_time{1};
        bRoot = 'data-';
        behFName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\' bRoot mouse '-' expDate '-' expTime '.mat'];
        load(behFName);
        ws = wheelSpeedCalc(input,32,expt(this_day).wheelColor); 
        ws_trial = reshape(ws,90,[]);
        ws_stimoff = ws_trial(31:60,:);
        haveRunning = mean(ws_stimoff,1) > 2;
        run_by_sesh{sesh,1} = haveRunning;
    end
    runTrials_SST{idx,1} = run_by_sesh;
end

% save(fullfile(fn_epi,'runTrials_SST'),'runTrials_SST','-v7.3');
save(fullfile('G:\home\ACh\Analysis\2p_analysis\epileptiform_analysis','runTrials_SST'),'runTrials_SST','-v7.3');
