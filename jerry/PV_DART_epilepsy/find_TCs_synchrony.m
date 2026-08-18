clear all; clear global; close all;
clc

dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories

session_id_TH = [28 31 34 69 71]; % enter post-DART session IDs
session_id_CC = [169 177 183]; % enter post-DART session IDs

nOn = 30;
nOff = 60;
nTot = nOn+nOff;
%% Jerry's
ds = 'DART_expt_info';
eval(ds);

fnroot_TH = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\ACh\Analysis\2p_analysis\PV_YM90K';
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
        TCmat{id,1} = npSub_tc;
    end
    segment_TCs_PV{sesh,1} = TCmat;
    segment_TCs_PV{sesh,2} = expt(post_day).mouse;
end

nMice_PV = size(segment_TCs_PV,1);
nDay = 2;

%% number of cells active
iti_act_PV = cell(nMice_PV,nDay);
iti_synch_PV = cell(nMice_PV,nDay);
for iM = 1:nMice_PV
    for iD = 1:nDay
        nCells = size(segment_TCs_PV{iM,1}{iD},2);
        trial_tc = reshape(segment_TCs_PV{iM,1}{iD},nTot,[],nCells);
        iti_resp = reshape(trial_tc(1+nOff/2:nOff,:,:),nOff/2,[],10,nCells);
        sz = size(iti_resp);
        iti_resp = reshape(iti_resp,sz(1)*sz(2),10,nCells);
        iti_act_PV{iM,iD} = zeros(size(iti_resp));
        nFrames = size(iti_resp,1);
        for iT = 1:10
            for iC = 1:nCells
                iti_win = iti_resp(:,iT,iC);
                pct = prctile(iti_win,75);
                iti_pct = iti_win(find(iti_win<pct));
                iti_pct_mean = mean(iti_pct);
                iti_pct_std = std(iti_pct);
                temp_ind = find(iti_win>iti_pct_mean+(iti_pct_std*5));
                iti_act_PV{iM,iD}(temp_ind,iT,iC) = 1;
            end
        end
        iti_act_PV{iM,iD} = reshape(iti_act_PV{iM,iD},nFrames*10, nCells);
        iti_synch_PV{iM,iD} = sum(iti_act_PV{iM,iD},2)./nCells;
    end
end

start = 1;
figure;
for iM = 1:nMice_PV
    for iD = 1:nDay
        subplot(nMice_PV,nDay,start)
        histogram(iti_synch_PV{iM,iD},[0:0.01:1])
        xlim([0 1])
        ylim([0 100])
        start = start+1;
        med_synch(iM,iD) = mean(iti_synch_PV{iM,iD});
        if iD == 1
            title(segment_TCs_PV{iM,2})
        end
        if iD == 2
            title(num2str(chop((med_synch(iM,2)./med_synch(iM,1)),3)))
        end
    end
end







%% Celine's 
ds = 'DART_V1_contrast_ori_Celine';
eval(ds);

fnroot_CC = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\ACh\Analysis\2p_analysis\SST_YM90K';
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

nMice_SST = size(segment_TCs_SST,1);
nDay = 2;

%%
iti_act_SST = cell(nMice_SST,nDay);
iti_synch_SST = cell(nMice_SST,nDay);
for iM = 1:nMice_SST
    for iD = 1:nDay
        nCells = size(segment_TCs_SST{iM,1}{iD},2);
        trial_tc = reshape(segment_TCs_SST{iM,1}{iD},nTot,[],nCells);
        iti_resp = reshape(trial_tc(1+nOff/2:nOff,:,:),nOff/2,[],10,nCells);
        sz = size(iti_resp);
        iti_resp = reshape(iti_resp,sz(1)*sz(2),10,nCells);
        iti_act_SST{iM,iD} = zeros(size(iti_resp));
        nFrames = size(iti_resp,1);
        for iT = 1:10
            for iC = 1:nCells
                iti_win = iti_resp(:,iT,iC);
                pct = prctile(iti_win,75);
                iti_pct = iti_win(find(iti_win<pct));
                iti_pct_mean = mean(iti_pct);
                iti_pct_std = std(iti_pct);
                temp_ind = find(iti_win>iti_pct_mean+(iti_pct_std*5));
                iti_act_SST{iM,iD}(temp_ind,iT,iC) = 1;
            end
        end
        iti_act_SST{iM,iD} = reshape(iti_act_SST{iM,iD},nFrames*10, nCells);
        iti_synch_SST{iM,iD} = sum(iti_act_SST{iM,iD},2)./nCells;
    end
end

start = 1;
figure;
for iM = 1:nMice_SST
    for iD = 1:nDay
        subplot(nMice_SST,nDay,start)
        histogram(iti_synch_SST{iM,iD},[0:0.01:1])
        xlim([0 1])
        ylim([0 100])
        start = start+1;
        med_synch(iM,iD) = mean(iti_synch_SST{iM,iD});
        if iD == 1
            title(segment_TCs_SST{iM,2})
        end
        if iD == 2
            title(num2str(chop((med_synch(iM,2)./med_synch(iM,1)),3)))
        end
    end
end