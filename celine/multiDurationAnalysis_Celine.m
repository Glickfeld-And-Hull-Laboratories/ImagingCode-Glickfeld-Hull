% to analyze individual day's data
clear all; clear global; close all

%identifying animal and run
date = '211025';
imgFolder = '003';
time = '1543';
mouse = 'i499';

%setting my paths
fn_base = 'Z:\home\Celine\Analysis\2p_analysis\';


fn = fullfile(fn_base,mouse,date,imgFolder);
mkdir(fn);
cd(fn);

run_str = catRunName(imgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];

%% load data

load(fullfile([datemouserun '_trial_TCs.mat']));
load(fullfile([datemouserun '_TCs.mat']));
beh_prefix = strcat('Z:\Behavior\Data\data-');
beh_file = [beh_prefix mouse '-' date '-' time '.mat'];
load(beh_file); %load the mworks behavioral file


stimOneTime = celleqel2mat_padded(input.tStimOneGratingOnTimeMs); %duration for each trial
stimOneTimes = unique(stimOneTime); %list of durations used
nTime = length(stimOneTimes); %how many different durations were used
ntrials = size(input.tThisTrialStartTimeMs,2);
frame_rate = double(input.frameRateHz);
nCells=size(data_tc_trial,3);


%right now I'm only doing one contrast
stimOneContrast = celleqel2mat_padded(input.tStimOneGratingContrast);
stimOneCons = unique(stimOneContrast);

%right now I'm only doing one orientation
tOri = celleqel2mat_padded(input.tStimOneGratingDirectionDeg);
oris = unique(tOri);
nOri = length(oris);

%% re-split into trials
nFrames = size(npSub_tc,1);
nTrials = size(input.tThisTrialStartTimeMs,2); 
cStimOne = cell2mat(input.cStimOneOn);   
%split into trials using maximum length
nOnMax = max(stimOneTime)*(frame_rate/1000); %find the maximum number of frames for the on period
nFrameTrial = nOnMax+3*frame_rate; %to add one second on either side

data_tc_trial = nan(nFrameTrial,nTrials,nCells); %make empty datafrmae
for iTrial = 1:nTrials
    if cStimOne(iTrial)+nOnMax+frame_rate < nFrames %to remove trials too close to the end
        tempF =  mean(npSub_tc(cStimOne(iTrial)-20 :cStimOne(iTrial)-1,:),1);
        data_tc_trial_temp = npSub_tc(cStimOne(iTrial)-frame_rate : cStimOne(iTrial)+nOnMax+(2*frame_rate)-1,:);
        data_tc_trial(:,iTrial,:) = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial_temp, tempF), tempF);
    end 
end

%% identifying visually responsive cells
% for each duration

% data_resp = zeros(nCells,nTime,2);
% h = zeros(nCells, nTime);
% p = zeros(nCells, nTime);
% % to use a response window that is half the stimulus window
% for  it = 1:nTime
%    x = round(stimOneTimes(it)*frame_rate/1000)
%    resp_win = frame_rate+1:frame_rate+x;
%    base_win = frame_rate/2 : frame_rate;%half second before the stim
%    inds = find(stimOneTime == stimOneTimes(it));
%    data_resp(:,it,1) =squeeze(nanmean(nanmean(data_tc_trial(resp_win,inds,:)),2));
%    data_resp(:,it,2) =squeeze(std(nanmean(data_tc_trial(resp_win,inds,:)),[],2))./sqrt(length(inds));
%    
%    [h(:,it), p(:,it)] = ttest(squeeze(nanmean(data_tc_trial(resp_win,inds,:),1)), squeeze(nanmean(data_tc_trial(base_win,inds,:),1)),'dim',1,'tail','right','alpha',0.05./(nOri.*3-1));
%     
% end

data_resp = zeros(nCells,nOri,nTime,2);
h = zeros(nCells, nOri,nTime);
p = zeros(nCells, nOri,nTime);
pref_ori = zeros(1,nCells);

for  iOri = 1:nOri
    for iTime = 1:nTime
       %resp_win = frame_rate+1:frame_rate+(max(stimOneTimes)/2000)*frame_rate;
       x = round(stimOneTimes(iTime)*frame_rate/1000);
       resp_win = frame_rate+1:frame_rate+x;
       base_win = frame_rate/2 : frame_rate;%half second before the stim
       %inds = find(tOri == oris(iOri)); 
       inds = intersect(find(tOri == oris(iOri)),find(stimOneTime == stimOneTimes(iTime)));
       data_resp(:,iOri,iTime,1) =squeeze(nanmean(nanmean(data_tc_trial(resp_win,inds,:)),2));
       data_resp(:,iOri,iTime,2) =squeeze(std(nanmean(data_tc_trial(resp_win,inds,:)),[],2))./sqrt(length(inds));
       [h(:,iOri,iTime), p(:,iOri,iTime)] = ttest(squeeze(nanmean(data_tc_trial(resp_win,inds,:),1)), squeeze(nanmean(data_tc_trial(base_win,inds,:),1)),'dim',1,'tail','right','alpha',0.05./(nOri.*3-1));
       [max_val, pref_ori] = max(mean(data_resp(:,:,:,1),3),[],2);
    end
   
end

h_all = sum(h,2:3);
resp=logical(h_all);
nResp = sum(resp)
respInds = find(resp);
%% put in a tuning curve averaging over stimulus duration
figure;
movegui('center')
start = 1;


if nResp<36
    [n n2] = subplotn(nResp);
    tot = n.*n2;
else
    n = 6;
    n2 = 6;
    tot = 36;
end

for iCell = 1:nResp
    if start>tot
        figure; movegui('center')
        start = 1;
    end
    subplot(n,n2,start)
    for iTime = 1:nTime
            errorbar(oris, data_resp(iCell,:,iTime,1), data_resp(iCell,:,iTime,2))
            hold on
            title (['Cell ' num2str(iCell)]);
    end
        
    start= start+1;
end
    
    %ylim([-0.1 inf])

   
%% looking at time courses 

% at each cell at it's preferred ori

t = 1:(size(data_tc_trial,1));
t=(t-frame_rate)/frame_rate;

[n n2] = subplotn(nTime);
figure;
    for it = 1:nTime %loop through the stim durations
        inds1 = find(stimOneTime == stimOneTimes(it)); %find trials with that duration
        for iCell = 1:nResp
        inds2 = find(tOri == oris(pref_ori(respInds(iCell))));
        inds = intersect(inds1,inds2);
        cell_temp_trials = squeeze(nanmean(data_tc_trial(:,inds,respInds(iCell)),2));
        temp_trials(:,iCell) = cell_temp_trials;
        end
        subplot(n,n2,it)
        plot(t, temp_trials)
        hold on
        plot(t, mean(squeeze(nanmean(data_tc_trial(:,inds,resp),2)),2),'k');
        hold on
        ylim([-.2 .25])
        x = [0 stimOneTimes(it)/1000 stimOneTimes(it)/1000 0];
        y = [.18 .18 .2 .2 ];
        patch(x,y,'b')
        hold off
        title(num2str(stimOneTimes(it)))
    end
    

print(fullfile(fn, [datemouserun '_prefOri_timecourses']),'-dpdf');

figure;
    for it = 1:nTime %loop through the stim durations
        inds1 = find(stimOneTime == stimOneTimes(it)); %find trials with that duration
        for iCell = 1:nResp
        inds2 = find(tOri == oris(pref_ori(respInds(iCell))));
        inds = intersect(inds1,inds2);
        cell_temp_trials = squeeze(nanmean(data_tc_trial(:,inds,respInds(iCell)),2));
        temp_trials(:,iCell) = cell_temp_trials;
        end
        temp_mean = nanmean(temp_trials,2);
        temp_se = std(temp_trials,[],2)/sqrt(sum(resp));

        subplot(n,n2,it)
        shadedErrorBar(t,temp_mean,temp_se);
        hold on
        ylim([-.05 .1])
        x = [0 stimOneTimes(it)/1000 stimOneTimes(it)/1000 0];
        y = [.08 .08 .082 .082 ];
        patch(x,y,'b')
        hold off
        title(num2str(stimOneTimes(it)))
    end
print(fullfile(fn, [datemouserun '_prefOri_shadedEB_timecourses']),'-dpdf');
clear temp_trials temp_mean temp_se inds inds1 inds2 iCell it
%% sliding average 
data_tc_trial_smooth = movmean(data_tc_trial,2,1); 


[n n2] = subplotn(nTime);
figure;
    for it = 1:nTime %loop through the stim durations
        inds = find(stimOneTime == stimOneTimes(it)); %find trials with that duration
        temp_trials = squeeze(nanmean(data_tc_trial_smooth(:,inds,resp),2));
        subplot(n,n2,it)
        plot(t, temp_trials)
        hold on
        plot(t, mean(squeeze(nanmean(data_tc_trial_smooth(:,inds,resp),2)),2),'k');
        hold on
        ylim([-.2 .25])
        x = [0 stimOneTimes(it)/1000 stimOneTimes(it)/1000 0];
        y = [.18 .18 .19 .19 ];
        patch(x,y,'b')
        hline(0)
        hold off
        title(num2str(stimOneTimes(it)))
    end
    
print(fullfile(fn, [datemouserun 'smoothed_timecourses']),'-dpdf'); 

figure;
    for it = 1:nTime %loop through the stim durations
        inds = find(stimOneTime == stimOneTimes(it)); %find trials with that duration
        temp_mean = mean(squeeze(nanmean(data_tc_trial_smooth(:,inds,resp),2)),2);
        temp_se= std(squeeze(nanmean(data_tc_trial_smooth(:,inds,resp),2)),[],2)/sqrt(sum(resp));
        subplot(n,n2,it)
        shadedErrorBar(t,temp_mean,temp_se);
        hold on
        ylim([-.03 .05])
        x = [0 stimOneTimes(it)/1000 stimOneTimes(it)/1000 0];
        y = [.046 .046 .048 .048 ];
        patch(x,y,'b')
        hline(0)
        hold off
        title(num2str(stimOneTimes(it)))
    end
    
print(fullfile(fn, [datemouserun 'smoothed_shadedEB_timecourses']),'-dpdf');
%% proceed with smoothed data
%width at half max of response
%time of 0 crossing
%max
%min
%% what does the activity of individual cells look like?
close all
start = 1;
if nResp<36
    [n n2] = subplotn(nResp);
    tot = n.*n2;
else
    n = 6;
    n2 = 6;
    tot = 36;
end
respInds = find(resp);


for iCell = 1:nResp
    if start>tot
        figure; movegui('center')
        start = 1;
    end
    subplot(n,n2,start)

    inds1 = find(stimOneTime == stimOneTimes(2)); %find trials with that duration
    inds2 = find(tOri == oris(pref_ori(respInds(iCell))));
    inds = intersect(inds1,inds2);
    temp_mean = mean(nanmean(data_tc_trial(:,inds,respInds(iCell)),2),2);
    temp_se= nanstd(data_tc_trial(:,inds,respInds(iCell)),[],2)/sqrt(length(inds));
    shadedErrorBar(t,temp_mean,temp_se);    
    hold on
    title (['Cell ' num2str(respInds(iCell))]);


    start= start+1;
end


    