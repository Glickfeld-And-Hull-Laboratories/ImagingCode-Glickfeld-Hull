clear all; clear global; close all;
clc

dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories

session_id_TH = [28 31 34 69 71]; % enter post-DART session IDs
session_id_CC = [169 177 183]; % enter post-DART session IDs


%% Jerry's
ds = 'DART_expt_info';
eval(ds);

fnroot_TH = 'G:\home\ACh\Analysis\2p_analysis\PV_YM90K';
segment_TCs_PV = cell(length(session_id_TH),2);

for sesh = 1:length(session_id_TH)
    TCmat = cell(2,1);
    pre_day = expt(session_id_TH(sesh)).multiday_matchdays;
    post_day = session_id_TH(sesh);
    this_days = [pre_day post_day];
    for id = 1:length(this_days)
        this_session = this_days(id);
        tc_fn = fullfile(fnroot_TH,expt(this_session).mouse,expt(this_session).date,expt(this_session).contrastxori_runs,'TCs.mat');
        load(cell2mat(tc_fn));
        TCmat{id,1} = data_tc;
    end
    segment_TCs_PV{sesh,1} = TCmat;
    segment_TCs_PV{sesh,2} = expt(post_day).mouse;
end

%% Celine's 
ds = 'DART_V1_contrast_ori_Celine';
eval(ds);

fnroot_CC = 'G:\home\ACh\Analysis\2p_analysis\SST_YM90K';
segment_TCs_SST = cell(length(session_id_CC),2);

for sesh = 1:length(session_id_CC)
    TCmat = cell(2,1);
    pre_day = expt(session_id_CC(sesh)).multiday_matchdays;
    post_day = session_id_CC(sesh);
    this_days = [pre_day post_day];
    disp(this_days);
    for id = 1:length(this_days)
        this_session = this_days(id);
        tc_fn = fullfile(fnroot_CC,expt(this_session).mouse,expt(this_session).date,expt(this_session).contrastxori_runs,'TCs.mat');
        load(cell2mat(tc_fn));
        TCmat{id,1} = data_tc;
    end
    segment_TCs_SST{sesh,1} = TCmat;
    segment_TCs_SST{sesh,2} = expt(post_day).mouse;
end

