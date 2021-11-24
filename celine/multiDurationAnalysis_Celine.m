% to analyze individual day's data
clear all; clear global; close all

%identifying animal and run
date = '211120';
imgFolder = '003';
time = '1753';
mouse = 'WK08';

%setting my paths
fn_base = 'Z:\home\Celine\Analysis\2p_analysis\';


fn = fullfile(fn_base,mouse,date,imgFolder);
mkdir(fn);
cd(fn);

run_str = catRunName(imgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];

%% load data

%load(fullfile([datemouserun '_trial_TCs.mat']));
load(fullfile([datemouserun '_TCs.mat']));
load(fullfile([datemouserun '_mask_cell.mat']));

beh_prefix = strcat('Z:\Behavior\Data\data-');
beh_file = [beh_prefix mouse '-' date '-' time '.mat'];
load(beh_file); %load the mworks behavioral file


stimOneTime = celleqel2mat_padded(input.tStimOneGratingOnTimeMs); %duration for each trial
stimOneTimes = unique(stimOneTime); %list of durations used
nTime = length(stimOneTimes); %how many different durations were used
ntrials = size(input.tThisTrialStartTimeMs,2);
frame_rate = double(input.frameRateHz);
nCells=size(npSub_tc,2);


%right now I'm only doing one contrast
stimOneContrast = celleqel2mat_padded(input.tStimOneGratingContrast);
stimOneCons = unique(stimOneContrast);

%right now I'm only doing one orientation
tOri = celleqel2mat_padded(input.tStimOneGratingDirectionDeg);
oris = unique(tOri);
nOri = length(oris);

%% looking at masks
inputfile_temp = input;
clear input
cell_stats=regionprops(mask_cell);
area=[cell_stats.Area];
hist(area)
cutoff=input('enter soma/dendrite cuttoff value:');
dend=find(area < cutoff-5);
soma=find(area > cutoff+5);
input = inputfile_temp;
length(soma)
length(dend)
clear inputfile_temp
%% looking at wheel speed
[wheel_speed] = wheelSpeedCalc(input,32,'purple'); 
figure; plot(wheel_speed)
print(fullfile(fn, [datemouserun 'running']),'-dpdf');
mean(wheel_speed)


%% re-split into trials
nFrames = size(npSub_tc,1);
nTrials = size(input.tThisTrialStartTimeMs,2); 
cStimOne = cell2mat(input.cStimOneOn);   
%split into trials using maximum length
nOnMax = max(stimOneTime)*(frame_rate/1000); %find the maximum number of frames for the on period
nFrameTrial = nOnMax+3*frame_rate; %to add one second on either side

data_tc_trial = nan(nFrameTrial,nTrials,nCells); %make empty dataframe
wheel_tc = nan(nFrameTrial,nTrials);
for iTrial = 1:nTrials
    if cStimOne(iTrial)+nOnMax+(2*frame_rate) < nFrames %to remove trials too close to the end
        tempF =  mean(npSub_tc(cStimOne(iTrial)-20 :cStimOne(iTrial)-1,:),1);
        data_tc_trial_temp = npSub_tc(cStimOne(iTrial)-frame_rate : cStimOne(iTrial)+nOnMax+(2*frame_rate)-1,:);
        data_tc_trial(:,iTrial,:) = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial_temp, tempF), tempF);
        wheel_tc(:,iTrial) = wheel_speed(1,cStimOne(iTrial)-frame_rate : cStimOne(iTrial)+nOnMax+(2*frame_rate)-1);
    end 
end

wheel_trial_avg = mean(wheel_tc(frame_rate-10:nOnMax+10,:),1);
RIx = wheel_trial_avg>1.5;
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
min_ori = zeros(1,nCells);
threshPass=zeros(nCells, nOri,nTime);

for  iOri = 1:nOri
    for iTime = 1:nTime
       %resp_win = frame_rate+1:frame_rate+(max(stimOneTimes)/2000)*frame_rate;
       x = round(stimOneTimes(iTime)*frame_rate/1000);
       resp_win = frame_rate+5 : frame_rate+9;
%        resp_win = frame_rate+1:frame_rate+x;
       base_win = frame_rate/2 : frame_rate;%half second before the stim
       %inds = find(tOri == oris(iOri)); 
       inds = intersect(find(tOri == oris(iOri)),find(stimOneTime == stimOneTimes(iTime)));
       data_resp(:,iOri,iTime,1) =squeeze(nanmean(nanmean(data_tc_trial(resp_win,inds,:)),2));
       data_resp(:,iOri,iTime,2) =squeeze(std(nanmean(data_tc_trial(resp_win,inds,:)),[],2))./sqrt(length(inds));
       [h(:,iOri,iTime), p(:,iOri,iTime)] = ttest(squeeze(nanmean(data_tc_trial(resp_win,inds,:),1)), squeeze(nanmean(data_tc_trial(base_win,inds,:),1)),'dim',1,'tail','right','alpha',0.05./(nOri.*3-1));
       [max_val, pref_ori] = max(mean(data_resp(:,:,:,1),3),[],2);
       [min_val, min_ori] = min(mean(data_resp(:,:,:,1),3),[],2);
       for iCell = 1:nCells
          thresh =  1.5*std(nanmean(data_tc_trial(base_win,inds,iCell),2));
          threshPass(iCell,iOri,iTime)=logical(data_resp(iCell,iOri,iTime,1)>thresh);
       end
    end
   
end

h_all = sum(h,2:3);
h_resp=logical(h_all);

thresh_all=sum(threshPass,2:3);
thresh_resp=logical(thresh_all);


respInds = find(h_resp);
nonRespInds = find(~h_resp);
nResp=length(respInds)

countsTable = table([nCells], [sum(h_resp)],[sum(thresh_resp)],'VariableNames',{'total cells' 'responsive cells ttest' 'responsive cells threshold'})
writetable(countsTable,fullfile(fn,'counts.csv'))

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
            errorbar(oris, data_resp(respInds(iCell),:,iTime,1), data_resp(respInds(iCell),:,iTime,2))
            hold on
            if ismember(respInds(iCell),dend)
                titleExtra = 'dend';
            else
                titleExtra = 'soma';
            end
            title (['Cell ' num2str(respInds(iCell)) titleExtra]);
    end
        
    start= start+1;
end
    
    %ylim([-0.1 inf])

   
%% looking at time courses - somas

% at each cell at it's preferred ori

t = 1:(size(data_tc_trial,1));
t=(t-frame_rate)/frame_rate;

respSoma = intersect(respInds,soma);

[n n2] = subplotn(nTime);
figure;
    for it = 1:nTime %loop through the stim durations
        inds1 = find(stimOneTime == stimOneTimes(it)); %find trials with that duration
        for iCell = 1:length(respSoma)
        inds2 = find(tOri == oris(pref_ori(respSoma(iCell))));
        inds = intersect(inds1,inds2);
        cell_temp_trials = squeeze(nanmean(data_tc_trial(:,inds,respSoma(iCell)),2));
        temp_trials(:,iCell) = cell_temp_trials;
        end
        subplot(n,n2,it)
        plot(t, temp_trials)
        hold on
        plot(t, mean(squeeze(nanmean(temp_trials,2)),2),'k');
        hold on
        ylim([-.2 .45])
        x = [0 stimOneTimes(it)/1000 stimOneTimes(it)/1000 0];
        y = [.43 .43 .432 .432 ];
        patch(x,y,'b')
        hline(0)
        hold off
        title(num2str(stimOneTimes(it)))
    end
    

print(fullfile(fn, [datemouserun '_prefOri_timecourses']),'-dpdf');

figure;
    for it = 1:nTime %loop through the stim durations
        inds1 = find(stimOneTime == stimOneTimes(it)); %find trials with that duration
        for iCell = 1:length(respSoma)
        inds2 = find(tOri == oris(pref_ori(respSoma(iCell))));
        inds = intersect(inds1,inds2);
        cell_temp_trials = squeeze(nanmean(data_tc_trial(:,inds,respSoma(iCell)),2));
        temp_trials(:,iCell) = cell_temp_trials;
        end
        temp_mean = nanmean(temp_trials,2);
        temp_se = std(temp_trials,[],2)/sqrt(length(respSoma));

        subplot(n,n2,it)
        shadedErrorBar(t,temp_mean,temp_se);
        hold on
        ylim([-.05 .15])
        x = [0 stimOneTimes(it)/1000 stimOneTimes(it)/1000 0];
        y = [.13 .13 .132 .132 ];
        patch(x,y,'b')
        hline(0)
        hold off
        title(num2str(stimOneTimes(it)))
    end
print(fullfile(fn, [datemouserun '_prefOri_shadedEB_timecourses']),'-dpdf');
clear temp_trials temp_mean temp_se inds inds1 inds2 iCell it
%% dendrites

t = 1:(size(data_tc_trial,1));
t=(t-frame_rate)/frame_rate;

[n n2] = subplotn(nTime);
figure;
    for it = 1:nTime %loop through the stim durations
        inds1 = find(stimOneTime == stimOneTimes(it)); %find trials with that duration
        for iCell = 1:length(dend)
        inds2 = find(tOri == oris(pref_ori(dend(iCell))));
        inds = intersect(inds1,inds2);
        cell_temp_trials = squeeze(nanmean(data_tc_trial(:,inds,dend(iCell)),2));
        temp_trials(:,iCell) = cell_temp_trials;
        end
        subplot(n,n2,it)
        plot(t, temp_trials)
        hold on
        plot(t, mean(squeeze(nanmean(temp_trials,2)),2),'k');
        hold on
        ylim([-.2 .45])
        x = [0 stimOneTimes(it)/1000 stimOneTimes(it)/1000 0];
        y = [.43 .43 .432 .432 ];
        patch(x,y,'b')
        hline(0)
        hold off
        title(num2str(stimOneTimes(it)))
    end
    

print(fullfile(fn, [datemouserun 'dend_prefOri_timecourses']),'-dpdf');

figure;
    for it = 1:nTime %loop through the stim durations
        inds1 = find(stimOneTime == stimOneTimes(it)); %find trials with that duration
        for iCell = 1:length(dend)
        inds2 = find(tOri == oris(pref_ori(dend(iCell))));
        inds = intersect(inds1,inds2);
        cell_temp_trials = squeeze(nanmean(data_tc_trial(:,inds,dend(iCell)),2));
        temp_trials(:,iCell) = cell_temp_trials;
        end
        temp_mean = nanmean(temp_trials,2);
        temp_se = std(temp_trials,[],2)/sqrt(length(dend));

        subplot(n,n2,it)
        shadedErrorBar(t,temp_mean,temp_se);
        hold on
        ylim([-.05 .15])
        x = [0 stimOneTimes(it)/1000 stimOneTimes(it)/1000 0];
        y = [.13 .13 .132 .132 ];
        patch(x,y,'b')
        hline(0)
        hold off
        title(num2str(stimOneTimes(it)))
    end
print(fullfile(fn, [datemouserun 'dend_prefOri_shadedEB_timecourses']),'-dpdf');
clear temp_trials temp_mean temp_se inds inds1 inds2 iCell it


%% at off-preferred ori


figure;
    for it = 1:nTime %loop through the stim durations
        inds1 = find(stimOneTime == stimOneTimes(it)); %find trials with that duration
        for iCell = 1:length(respSoma)
            pref = pref_ori(respSoma(iCell));
            if pref < 4
                off_ori = pref+1;
            else
                off_ori = 1;
            end
            inds2 = find(tOri == oris(off_ori));
            inds = intersect(inds1,inds2);
            cell_temp_trials = squeeze(nanmean(data_tc_trial(:,inds,respSoma(iCell)),2));
            temp_trials(:,iCell) = cell_temp_trials;
        end
        temp_mean = nanmean(temp_trials,2);
        temp_se = std(temp_trials,[],2)/sqrt(length(respSoma));

        subplot(n,n2,it)
        shadedErrorBar(t,temp_mean,temp_se);
        hold on
        ylim([-.05 .05])
        x = [0 stimOneTimes(it)/1000 stimOneTimes(it)/1000 0];
        y = [.045 .045 .0453 .0453 ];
        patch(x,y,'b')
        hline(0)
        hold off
        title(num2str(stimOneTimes(it)))
    end
    
print(fullfile(fn, [datemouserun '_offOri_shadedEB_timecourses']),'-dpdf');
clear temp_trials temp_mean temp_se inds inds1 inds2 iCell it
%% running vs. not running, all oris
figure;
    for it = 1:nTime %loop through the stim durations
        inds1 = find(stimOneTime == stimOneTimes(it)); %find trials with that duration
        for iCell = 1:length(respSoma)
        inds2 = find(RIx);
        inds = intersect(inds1,inds2);
        cell_temp_trials = squeeze(nanmean(data_tc_trial(:,inds,respSoma(iCell)),2));
        temp_trials(:,iCell) = cell_temp_trials;
        inds3 = find(~RIx);
        inds4 = intersect(inds1,inds3);
        cell_temp_trials2 = squeeze(nanmean(data_tc_trial(:,inds4,respSoma(iCell)),2));
        temp_trials2(:,iCell) = cell_temp_trials2;
        end
        temp_mean = nanmean(temp_trials,2);
        temp_se = std(temp_trials,[],2)/sqrt(length(respSoma));
        
        temp_mean2 = nanmean(temp_trials2,2);
        temp_se2 = std(temp_trials2,[],2)/sqrt(length(respSoma));

        subplot(n,n2,it)
        shadedErrorBar(t,temp_mean,temp_se);
        hold on
        shadedErrorBar(t,temp_mean2,temp_se2,'b');
        ylim([-.05 .15])
        x = [0 stimOneTimes(it)/1000 stimOneTimes(it)/1000 0];
        y = [.13 .13 .132 .132 ];
        patch(x,y,'b')
        hline(0)
        
        hold off
        title(num2str(stimOneTimes(it)))
    end
print(fullfile(fn, [datemouserun 'running_shadedEB_timecourses']),'-dpdf');


%% what does the activity of individual cells look like?



for iTime = 1:nTime
    figure
    start = 1;
    if nResp<36
        [n n2] = subplotn(length(respSoma));
        tot = n.*n2;
    else
        n = 6;
        n2 = 6;
        tot = 36;
    end
    


    for iCell = 1:(length(respSoma))
        if start>tot
            figure; movegui('center')
            start = 1;
        end
        subplot(n,n2,start)

        inds1 = find(stimOneTime == stimOneTimes(iTime)); %find trials with that duration
        inds2 = find(tOri == oris(pref_ori((respSoma(iCell)))));
        inds = intersect(inds1,inds2);
        temp_mean = mean(nanmean(data_tc_trial(:,inds,respInds(iCell)),2),2);
        temp_se= nanstd(data_tc_trial(:,inds,respSoma(iCell)),[],2)/sqrt(length(inds));
        shadedErrorBar(t,temp_mean,temp_se);    
        hold on
        title (['Cell ' num2str(respSoma(iCell))]);


        start= start+1;
    end
end

%% define oscillation 

%% sliding average 
data_tc_trial_smooth = movmean(data_tc_trial,2,1); 

t = 1:(size(data_tc_trial,1));
t=(t-frame_rate)/frame_rate;

respSoma = intersect(respInds,soma);

[n n2] = subplotn(nTime);

figure;
    for it = 1:nTime %loop through the stim durations
        inds = find(stimOneTime == stimOneTimes(it)); %find trials with that duration
        temp_trials = squeeze(nanmean(data_tc_trial_smooth(:,inds,respInds),2));
        subplot(n,n2,it)
        plot(t, temp_trials)
        hold on
        plot(t, mean(squeeze(nanmean(data_tc_trial_smooth(:,inds,respInds),2)),2),'k');
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
        inds1 = find(stimOneTime == stimOneTimes(it)); %find trials with that duration
        for iCell = 1:length(respSoma)
        inds2 = find(tOri == oris(pref_ori(respSoma(iCell))));
        inds = intersect(inds1,inds2);
        cell_temp_trials = squeeze(nanmean(data_tc_trial_smooth(:,inds,respSoma(iCell)),2));
        temp_trials(:,iCell) = cell_temp_trials;
        end
        temp_mean = nanmean(temp_trials,2);
        temp_se = std(temp_trials,[],2)/sqrt(length(respSoma));
        
        
        subplot(n,n2,it)
        shadedErrorBar(t,temp_mean,temp_se);
        hold on
        %ylim([-.03 .05])
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