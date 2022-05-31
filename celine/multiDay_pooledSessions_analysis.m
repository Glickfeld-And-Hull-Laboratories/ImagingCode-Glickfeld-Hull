%multiDay_pooledSessions_analysis
%purpose: to anlayze pre- vs. post-DART data for matched cells, pooling
%data from multiple mice/imaging sessions

clear all; clear global; close all
clc
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc =  behavConstsDART; %directories
eval(ds);

sess_list = [131 133];%enter all the sessions you want to pool
nSess=length(sess_list);

nd=2;%hard coding for two days per experimental session
frame_rate = 15;

sess_title = string(sess_list(1));
for iSess = 2:nSess
    sess_title = strcat(sess_title,'_',string(sess_list(iSess)));
end

fnout = fullfile(rc.achAnalysis,strcat('pooled_', sess_title));
mkdir(fnout);

mouse_list = cell(1,nSess);
tCon_pooled = cell(1,nSess);
tOri_pooled = cell(1,nSess);
prefOri_pooled = cell(1,nSess);
prefCon_pooled = cell(1,nSess);
stimStart_pooled = []; %stimStart should be the same across sessions but I want to make sure
trial_avrg_tc = cell(1,nSess);
trial_tc = cell(1,nSess);
green_ind = cell(1,nSess);
red_ind = cell(1,nSess);

%to load each day and combine the data
for iSess = 1:nSess
    day_id = sess_list(iSess);
    
    mouse = expt(day_id).mouse;
    mouse_list{iSess} = mouse;
    if expt(day_id).multiday_timesincedrug_hours>0
        dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
    else
        dart_str = 'control';
    end

    fn_multi = fullfile(rc.achAnalysis,mouse,['multiday_' dart_str]);
    cd(fn_multi)
    load(fullfile(fn_multi,'tc_keep.mat'))
    tCon_pooled{iSess} = tCon_match;
    tOri_pooled{iSess} = tOri_match;
    stimStart_pooled = [stimStart_pooled, stimStart];
    trial_avrg_tc{iSess} = tc_trial_avrg_keep;
    trial_tc{iSess} = data_trial_keep;
    green_ind{iSess} = green_ind_keep;
    red_ind{iSess} = red_ind_keep;
    prefOri_pooled{iSess} = pref_ori_keep;
    prefCon_pooled{iSess} = pref_con_keep;
    clear 'sess_title' 'explanation1' 'pref_ori_keep' 'pref_con_keep' 'tOri_match' 'tCon_match' 'data_trial_keep' 'nTrials' 'tc_trial_avrg_keep'  'green_keep_logical' 'red_keep_logical' 'green_ind_keep' 'red_ind_keep' 'stimStart';
end
%% to find all the contrasts represented in the data across sessions -
%assumes two days per session
cons = [];
all_cons = cell(nSess,1);
for iSess = 1:nSess
   cons1 = unique(tCon_pooled{iSess}{1}); 
   cons2 = unique(tCon_pooled{iSess}{2}); 
   all_cons{iSess} =  intersect(cons1,cons2); %find contrasts that were used in both days of this experiment
   
   clear cons1 cons2
end


%figure out which contrasts were covered in all experiments

%cons = unique(cell2mat_padded(all_cons)); %use this if different mice
%experienced different contrasts
cons = unique(cell2mat(all_cons));
nCon = length(cons);

%check whether stimStart is the same for all sessions

if length(unique(stimStart_pooled))>1
   error('Error. \n stim start time do not match');
else
    stimStart=(unique(stimStart_pooled));
    clear stimStart_pooled
end
%should do the same to check frame rate

nGreen = 0;
nRed=0;
for iSess = 1:nSess
     nGreen=nGreen+length(green_ind{iSess});
     nRed=nRed+length(red_ind{iSess});
end

%% plotting timecourses of response to preferred ori
tc_green_avrg = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg = cell(1,nd); %same for red
tc_green_se = cell(1,nd); %this will be the se across all green cells
tc_red_se = cell(1,nd); %same for red
green_resp_full = cell(1,nd);
red_resp_full = cell(1,nd);
green_resp_firstSec = cell(1,nd);
red_resp_firstSec = cell(1,nd);
cellCounts=nan(2,nCon);
tc_green_all=cell(nd,nCon);
tc_red_all=cell(nd,nCon);

for id = 1:nd
    for iCon=1:nCon
        green_trials=[]; %start with empty array for each day X contrast
        red_trials=[]; 
        for iSess = 1:nSess
            
            conInd = find(all_cons{iSess}==cons(iCon)); %figure out which index will
            %correspond to the desired contrast in this session
            %grab all those trials for this session/day
            green_trials=[green_trials, trial_avrg_tc{iSess}{id}(:,green_ind{iSess},conInd)];
            red_trials=[red_trials, trial_avrg_tc{iSess}{id}(:,red_ind{iSess},conInd)];
            
        end
        tc_green_all{id,conInd}=green_trials;
        tc_green_avrg{id}(:,conInd)=nanmean(green_trials,2);
        green_resp_full{id}(:,conInd)=nanmean(green_trials(stimStart:(stimStart+30),:),1);
        green_resp_firstSec{id}(:,conInd)=nanmean(green_trials(stimStart:(stimStart+15),:),1);
        green_std=std(green_trials,[],2);
        tc_green_se{id}(:,conInd)=green_std/sqrt(size(green_trials,2));
        cellCounts(1,iCon)=size(green_trials,2);
        
        tc_red_all{id,conInd}=red_trials;
        tc_red_avrg{id}(:,conInd)=nanmean(red_trials,2);
        red_resp_full{id}(:,conInd)=nanmean(red_trials(stimStart:(stimStart+30),:),1);
        red_resp_firstSec{id}(:,conInd)=nanmean(red_trials(stimStart:(stimStart+15),:),1);
        red_std=std(red_trials,[],2);
        tc_red_se{id}(:,conInd)=red_std/sqrt(size(red_trials,2));
        cellCounts(2,iCon)=size(red_trials,2);
        clear green_std  red_std
        end
end
    %%


%create a time axis in seconds
t=1:(size(tc_green_avrg{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(1,2,1) %for the first day

shadedErrorBar(t,tc_red_avrg{2}(:,iCon),tc_red_se{2}(:,iCon),'r');
ylim([-.02 .3]);
hold on
shadedErrorBar(t,tc_green_avrg{2}(:,iCon),tc_green_se{2}(:,iCon));
title(['Pre-DART contrast = ' num2str(cons(iCon))])
txt1 = ['HT- ' num2str(cellCounts(1,iCon))];
text(-1.5,0.25,txt1);
txt2 = ['HT+ ' num2str(cellCounts(2,iCon))];
text(-1.5,0.23,txt2,'Color','r');
ylabel('dF/F') 
xlabel('s') 
line([0,2],[-.01,-.01])
axis square


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg{1}(:,iCon),tc_red_se{1}(:,iCon),'r');
ylim([-.02 .3]);
hold on
shadedErrorBar(t,tc_green_avrg{1}(:,iCon),tc_green_se{1}(:,iCon));
ylabel('dF/F') 
xlabel('s') 
title(['Post-DART contrast = ' num2str(cons(iCon))])
line([0,2],[-.01,-.01])
axis square

print(fullfile(fnout,[num2str(cons(iCon)) '_timecourses.pdf']),'-dpdf');
end 
clear txt1 txt2 

%% make a plot of individual timecourses 

setYmin = -.2; %indicate y axes you want
setYmax = 0.8;


for iCon = 1:nCon

figure
subplot(2,2,1)
plot(t, tc_green_all{2,iCon},'color',[0 0 0 .2])
hold on
plot(t, tc_green_avrg{2}(:,iCon),'color',[0 0 0],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 1')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
hold off


subplot(2,2,2)
plot(t, tc_red_all{2,iCon},'color',[.7 .05 .05 .2])
hold on
plot(t, tc_red_avrg{2}(:,iCon),'color',[.7 .05 .05],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 1')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
hold off


subplot(2,2,3)
plot(t, tc_green_all{1,iCon},'color',[0 0 0 .2])
hold on
plot(t, tc_green_avrg{1}(:,iCon),'color',[0 0 0],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 2')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
hold off


subplot(2,2,4)
plot(t, tc_red_all{1,iCon},'color',[.7 .05 .05 .2])
hold on
plot(t, tc_red_avrg{1}(:,iCon),'color',[.7 .05 .05],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 2')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
sgtitle(['contrast = '  num2str(cons(iCon))])
hold off
print(fullfile(fnout,[num2str(cons(iCon)) '_indiv_timecourses.pdf']),'-dpdf');
end

%% scatterplot of responses

for iCon = 1:nCon
figure; movegui('center') 
subplot(1,2,1)
scatter(green_resp_full{2}(:,iCon),green_resp_full{1}(:,iCon),'k')
% hold on
% scatter(green_resp_firstSec{2}(:,iCon),green_resp_firstSec{1}(:,iCon),'k','filled')
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 .5])
xlim([-.1 .5])
refline(1)
title('HT- ')
axis square
hold off


subplot(1,2,2)
scatter(red_resp_full{2}(:,iCon),red_resp_full{1}(:,iCon),'MarkerEdgeColor',[.7 .05 .05])
% hold on
% scatter(red_resp_firstSec{2}(:,iCon),red_resp_firstSec{1}(:,iCon),'filled','MarkerFaceColor',[.7 .05 .05])

ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 .5])
xlim([-.1 .5])

refline(1)
title('HT+')
axis square
hold off

sgtitle(num2str(cons(iCon)))
print(fullfile(fnout,[num2str(cons(iCon)) 'maxResp_crossDay.pdf']),'-dpdf','-bestfit')

end


%% time to peak
tHalfMax = cell(1,nd);

for id = 1:nd
    tHalfMaxTemp=nan(nKeep,1);

        for iCell=1:nKeep
            %pull the data for a given cell at a given contrast (pref ori)
            tempData=tc_trial_avrg_keep_allCon{id}(stimStart:stimStart+nOn,iCell);
            if resp_keep{id}(iCell)
                smoothData=smoothdata(tempData,'movmean',5) ;
                halfMax = max(smoothData(3:length(smoothData)))/2;
                tHalfMaxCell =double(min(find(smoothData>halfMax)))/double(frame_rate);
                if rem(iCell, 10) == 0 
                figure;plot(tempData)
                hold on
                plot(smoothData)
                hold on
                vline(min(find(smoothData>halfMax)))
                end
                if length(tHalfMaxCell)>0
                    tHalfMaxTemp(iCell)=tHalfMaxCell;
                end
            end
 
       end
    tHalfMax{id}=tHalfMaxTemp;
end
clear tHalfMaxCell tHalfMaxTemp tempData smoothData halfMax
%% scatter for tMax



figure; movegui('center') 
subplot(1,2,1)
scatter((tHalfMax{2}(green_ind_keep)),(tHalfMax{1}(green_ind_keep)),'k')
ylabel('post-DART half-max(s)')
xlabel('pre-DART half-max(s)')
ylim([0 2.5])
xlim([0 2.5])
refline(1)
title('HT- ')
axis square


subplot(1,2,2)
scatter((tHalfMax{2}(red_ind_keep)),(tHalfMax{1}(red_ind_keep)),'MarkerEdgeColor',[.7 .05 .05])

ylabel('post-DART half-max(s)')
xlabel('pre-DART half-max(s)')
ylim([0 2.5])
xlim([0 2.5])


refline(1)
title('HT+')
axis square

sgtitle('time to half max (s), averaged over contrast')
print(fullfile(fnout, ['tHalfMax_crossDay.pdf']),'-dpdf','-bestfit')