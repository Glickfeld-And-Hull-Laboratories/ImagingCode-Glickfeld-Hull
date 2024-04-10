
clear all; clear global; close all
clc
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc =  behavConstsDART; %directories
eval(ds);
%285 300 295 308 324 334 DART YM90K 
% 299 289 304 312 320 330
sess_list = [366 368];%enter1 all the sessions you want to concatenate
nSess=length(sess_list);

nd=2;%hard coding for two days per experimental session

% INDICATE THE PRE VS. POST DAYS DEPENDING ON THE ORDER OF MATCHING
prompt = 'Which sesson was used as reference for matching: 0- baseline, 1- post-DART';
            x = input(prompt);
            switch x
                case 0
                    pre=1; %baeline session, used as reference, is in the 1st position
                    post=2;
                    "baseline used as reference"
                case 1
                  pre=2;
                  post=1; %post-DART session, used as reference, is in the 1st position  
                  "post-DART used as reference"
            end
clear x prompt

targetCon = [.125, .25, .5, 1]%what contrast to extract for all data - must be one that all datasets had

frame_rate = 15;

sess_title = string(sess_list(1));
for iSess = 2:nSess
    sess_title = strcat(sess_title,'_',string(sess_list(iSess)));
end
d=string(datetime('today'));
if nSess == 1
         if expt(sess_list(1)).multiday_timesincedrug_hours>0
            dart_str = [expt(sess_list(1)).drug '_' num2str(expt(sess_list(1)).multiday_timesincedrug_hours) 'Hr'];
        else
            dart_str = 'control';
        end
        
        fnout = fullfile(rc.achAnalysis,expt(sess_list(1)).mouse,['multiday_' dart_str],d);
else
    fnout= fullfile(rc.achAnalysis,strcat('concat', sess_title),d);
end
mkdir(fnout);
cd(fnout)
clear d sess_title

%%
mice=[];
red_concat=[];
green_concat=[];
nKeep_concat=[];

tc_trial_avrg_stat_concat=cell(1,nd);
tc_trial_avrg_loc_concat=cell(1,nd);
conBySize_resp_stat_concat=cell(1,nd);
conBySize_resp_loc_concat=cell(1,nd);
h_concat=cell(1,nd);

for iSess = 1:nSess
    day_id = sess_list(iSess)
    mouse = expt(day_id).mouse;
    mice=[mice;mouse];

    if expt(day_id).multiday_timesincedrug_hours>0
        dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
    else
        dart_str = 'control';
    end
    fn_multi = fullfile(rc.achAnalysis,mouse,['multiday_' dart_str]);

    load(fullfile(fn_multi,'tc_keep.mat'))
    load(fullfile(fn_multi,'resp_keep.mat'))
    load(fullfile(fn_multi,'input.mat'))
   
    nKeep = size(tc_trial_avrg_stat{post},2);

    %tells the contrast, direction and orientation for each trial each day
    tCon_match = cell(1,nd);
    tSize_match = cell(1,nd);

    %find the contrasts, directions and orientations for each day
    for id = 1:nd
        tCon_match{id} = celleqel2mat_padded(input(id).tGratingContrast(1:nTrials(id)));
        tSize_match{id} = celleqel2mat_padded(input(id).tGratingDiameterDeg(1:nTrials(id)));
         
    end

    cons = unique(tCon_match{post});
    sharedCon=find(ismember(cons, targetCon));
    sizes = unique(tSize_match{1});


    nOn = input(1).nScansOn;
    nOff = input(1).nScansOff;

    red_concat = [red_concat, red_keep_logical];
    green_concat = [green_concat, green_keep_logical];
    nKeep_concat = [nKeep_concat,nKeep];

    
    for id = 1:nd
        tc_trial_avrg_stat_concat{id}=cat(2,tc_trial_avrg_stat_concat{id},tc_trial_avrg_stat{id});
        tc_trial_avrg_loc_concat{id}=cat(2,tc_trial_avrg_loc_concat{id},tc_trial_avrg_loc{id});
        conBySize_resp_stat_concat{id}=cat(1,conBySize_resp_stat_concat{id},conBySize_resp_stat_keep{id});
        conBySize_resp_loc_concat{id}=cat(1,conBySize_resp_loc_concat{id},conBySize_resp_loc_keep{id});
        h_concat{id}=cat(1,h_concat{id},h_keep{id});
    end

end
%

red_ind_concat = find(red_concat);
green_ind_concat = find(green_concat);
%

nCon = length(cons)
nSize = length(sizes)
nKeep_total = sum(nKeep_concat);

clear data_resp_keep data_trial_keep green_ind_keep red_ind_keep conBySize_resp_loc_keep conBySize_resp_stat_keep tc_trial_avrg_stat tc_trial_avrg_loc conBySize_resp_loc_match conBySize_resp_stat_match red_keep_logical green_keep_logical pref_responses_stat pref_responses_loc resp_keep mouse pref_con_keep pref_dir_keep pref_size_keep explanation1 explanation2
%% cell selection
% find cells that I have running and stationary data for on both days - these are called pass cells

haveRunning_pre=(logical(sum(sum(~isnan(conBySize_resp_loc_concat{pre}),3),2)==nSize)); %find cells that have a no NAN values for any sise (the total of non-NAN should == the number of size)
haveRunning_post=(logical(sum(sum(~isnan(conBySize_resp_loc_concat{post}),3),2)==nSize));
haveRunning_both= find(haveRunning_pre.* haveRunning_post); %find cells that meet this criteria for both days - now in indices, not logical

haveStat_pre=(logical(sum(sum(~isnan(conBySize_resp_stat_concat{pre}),3),2)==nSize)); %find cells that have a no NAN values for any sise (the total of non-NAN should == the number of size)
haveStat_post=(logical(sum(sum(~isnan(conBySize_resp_stat_concat{post}),3),2)==nSize));
haveStat_both= find(haveStat_pre.* haveStat_post); %find cells that meet this criteria for both days - now in indices, not logical

runningCells = intersect(haveStat_both, haveRunning_both);

%find cells responsive to a given size
mySize = 1; %1 is the smallest size
respToSizeBothDays = cell(1,nd);
for id = 1:nd
    respToSizeBothDays{id}=sum(h_concat{id}(:,:,:,mySize),2);
end
respToSmall = logical(respToSizeBothDays{pre}+respToSizeBothDays{post}); %to find cells that were responsive to this size on either day


runningGreen = intersect(runningCells, green_ind_concat);
runningRed= intersect(runningCells, red_ind_concat);

statGreen = green_ind_concat;
statRed = red_ind_concat;

clear haveRunning_pre haveRunning_post haveRunning_both haveStat_both haveStat_pre haveStat_post

%%
%get an array with indices for each mouse
mouseInds=cell(1,nSess);
start=1;
for iMouse = 1:nSess
    mouseInds{iMouse}=start:(start-1)+nKeep_concat(iMouse);
    start = start+nKeep_concat(iMouse)
end
clear start iMouse
% %%alternate for times when I know I don't have enough running data
% pass_green = green_ind_concat;
% runningRed = red_ind_concat;

%find how many haveRunning red cells I got for each mouse
cellCounts = nan(nSess,2);
mouseNames=[];
for iMouse = 1:nSess
   cellCounts(iMouse,1)=length(intersect(runningRed',(mouseInds{iMouse})));
   cellCounts(iMouse,2)=length(intersect(statRed',(mouseInds{iMouse})));
   mouseNames=[mouseNames, string(mice(iMouse,:))];
end


cellCountsGreen = nan(nSess,nCon);
mouseNames=[];
for iMouse = 1:nSess
   cellCountsGreen(iMouse,1)=length(intersect(runningGreen',(mouseInds{iMouse})));
   cellCountsGreen(iMouse,2)=length(intersect(statGreen',(mouseInds{iMouse})));

   mouseNames=[mouseNames, string(mice(iMouse,:))];
end


cellCountTable = table(cellCounts, RowNames=mouseNames)
cellCountTableGreen = table(cellCountsGreen, RowNames=mouseNames)
%writetable(cellCountTable,fullfile(fnout,['cellCounts.csv']),'WriteRowNames',true)
clear cellCounts cellCountsGreen

%% plotting some data


% make figure with se shaded, one figure per contrast - stationary

tc_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_stat = cell(1,nd); %same for red
tc_green_se_stat = cell(1,nd); %this will be the se across all green cells
tc_red_se_stat = cell(1,nd); %same for red




for id = 1:nd
  for iCon = 1:nCon
      for iSize = 1:nSize
        
        tc_green_avrg_stat{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_stat_concat{id}(:,statGreen,iCon,iSize),2);
        green_std=nanstd(tc_trial_avrg_stat_concat{id}(:,statGreen,iCon,iSize),[],2);
        tc_green_se_stat{id}(:,iCon,iSize)=green_std/sqrt(length(statGreen));
        
        tc_red_avrg_stat{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_stat_concat{id}(:,statRed,iCon,iSize),2);
        red_std=nanstd(tc_trial_avrg_stat_concat{id}(:,statRed,iCon,iSize),[],2);
        tc_red_se_stat{id}(:,iCon,iSize)=red_std/sqrt(length(statRed));
        
        clear green_std red_std
      end 
    end
end

z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
    for iSize = 1:nSize
    figure
    subplot(1,2,1) %for the first day
    
    
    
    ylim([-.05 .25]);
    hold on
    shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon,iSize),tc_green_se_stat{pre}(:,iCon,iSize),'k');
    hold on
    shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon,iSize),tc_green_se_stat{post}(:,iCon,iSize),'b','transparent');
    hold on
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    title(['-HTP',' n = ', num2str(length(statGreen))])
    
    ylabel('dF/F') 
    xlabel('s') 
    set(gca,'XColor', 'none','YColor','none')
    
    
    subplot(1,2,2) %+HTP
    shadedErrorBar(t,tc_red_avrg_stat{pre}(:,iCon,iSize),tc_red_se_stat{pre}(:,iCon,iSize),'k');
    hold on
    shadedErrorBar(t,tc_red_avrg_stat{post}(:,iCon,iSize),tc_red_se_stat{post}(:,iCon,iSize),'b');
    ylim([-.05 .25]);
    hold on
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    ylabel('dF/F') 
    xlabel('s') 
    title(['+HTP',' n = ', num2str(length(statRed))])
    
    x0=5;
    y0=5;
    width=4;
    height=3;
    set(gcf,'units','inches','position',[x0,y0,width,height])
    set(gca,'XColor', 'none','YColor','none')
    
    sgtitle(['stationary, con ' num2str(cons(iCon)) ' size ' num2str(sizes(iSize))])
    
    print(fullfile(fnout,[num2str(cons(iCon)) '_' num2str(sizes(iSize)) '_stat_cellType_timecourses.pdf']),'-dpdf');
    end
end 



%% plots for running trials

tc_green_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_loc = cell(1,nd); %same for red
tc_green_se_loc = cell(1,nd); %this will be the se across all green cells
tc_red_se_loc = cell(1,nd); %same for red



for id = 1:nd
  for iCon = 1:nCon
      for iSize = 1:nSize
        
        tc_green_avrg_loc{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_loc_concat{id}(:,runningGreen,iCon,iSize),2);
        green_std=nanstd(tc_trial_avrg_loc_concat{id}(:,runningGreen,iCon,iSize),[],2);
        tc_green_se_loc{id}(:,iCon,iSize)=green_std/sqrt(length(runningGreen));
        
        tc_red_avrg_loc{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_loc_concat{id}(:,runningRed,iCon,iSize),2);
        red_std=nanstd(tc_trial_avrg_loc_concat{id}(:,runningRed,iCon,iSize),[],2);
        tc_red_se_loc{id}(:,iCon,iSize)=red_std/sqrt(length(runningRed));
        
        clear green_std red_std
      end 
    end
end

z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_loc{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
    for iSize = 1:nSize
    figure
    subplot(1,2,1) %for the first day
    
    
    
    ylim([-.05 .35]);
    hold on
    shadedErrorBar(t,tc_green_avrg_loc{pre}(:,iCon,iSize),tc_green_se_loc{pre}(:,iCon,iSize),'k');
    hold on
    shadedErrorBar(t,tc_green_avrg_loc{post}(:,iCon,iSize),tc_green_se_loc{post}(:,iCon,iSize),'b','transparent');
    hold on
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    title(['-HTP',' n = ', num2str(length(runningGreen))])
    
    ylabel('dF/F') 
    xlabel('s') 
    set(gca,'XColor', 'none','YColor','none')
    
    
    subplot(1,2,2) %+HTP
    shadedErrorBar(t,tc_red_avrg_loc{pre}(:,iCon,iSize),tc_red_se_loc{pre}(:,iCon,iSize),'k');
    hold on
    shadedErrorBar(t,tc_red_avrg_loc{post}(:,iCon,iSize),tc_red_se_loc{post}(:,iCon,iSize),'b');
    ylim([-.05 .35]);
    hold on
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    ylabel('dF/F') 
    xlabel('s') 
    title(['+HTP',' n = ', num2str(length(runningRed))])
    
    x0=5;
    y0=5;
    width=4;
    height=3;
    set(gcf,'units','inches','position',[x0,y0,width,height])
    set(gca,'XColor', 'none','YColor','none')
    
    sgtitle(['running, con ' num2str(cons(iCon)) ' size ' num2str(sizes(iSize))])
    
    print(fullfile(fnout,[num2str(cons(iCon)) '_' num2str(sizes(iSize)) '_loc_cellType_timecourses.pdf']),'-dpdf');
    end
end 

%% population size tuning - averaged over contrast

%errorbar for stat resp and loc resp vs size, where error is across mice

%get average for each day

sizeResp_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
sizeResp_red_avrg_stat = cell(1,nd); %same for red
sizeResp_green_se_stat = cell(1,nd); %this will be the se across all green cells
sizeResp_red_se_stat = cell(1,nd); %same for red



for id = 1:nd
    green_data=squeeze(mean(conBySize_resp_stat_concat{id}(statGreen,:,:),2));%pulling the green cells and averaging over contrast
    sizeResp_green_avrg_stat{id}=nanmean(green_data,1);
    green_std=nanstd(green_data,1);
    sizeResp_green_se_stat{id}=green_std/sqrt(length(statGreen));
    
    
    red_data=squeeze(mean(conBySize_resp_stat_concat{id}(statRed,:,:),2));%pulling the red cells and averaging over contrast
    sizeResp_red_avrg_stat{id}=nanmean(red_data,1);
    red_std=nanstd(red_data,1);
    sizeResp_red_se_stat{id}=red_std/sqrt(length(statRed));
    
    clear green_std red_std green_data red_data
 
end



sizeResp_green_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
sizeResp_red_avrg_loc = cell(1,nd); %same for red
sizeResp_green_se_loc = cell(1,nd); %this will be the se across all green cells
sizeResp_red_se_loc = cell(1,nd); %same for red



for id = 1:nd
    green_data=squeeze(mean(conBySize_resp_loc_concat{id}(runningGreen,:,:),2));%pulling the green cells and averaging over contrast
    sizeResp_green_avrg_loc{id}=nanmean(green_data,1);
    green_std=nanstd(green_data,1);
    sizeResp_green_se_loc{id}=green_std/sqrt(length(runningGreen));
    
    red_data=squeeze(mean(conBySize_resp_loc_concat{id}(runningRed,:,:),2));%pulling the red cells and averaging over contrast
    sizeResp_red_avrg_loc{id}=nanmean(red_data,1);
    red_std=nanstd(red_data,1);
    sizeResp_red_se_loc{id}=red_std/sqrt(length(runningRed));
    
    clear green_std red_std green_data red_data
 
end



figure
subplot(2,2,1) %for the first day
errorbar(sizes,sizeResp_green_avrg_stat{pre},sizeResp_green_se_stat{pre},'k');
hold on
errorbar(sizes,sizeResp_green_avrg_stat{post},sizeResp_green_se_stat{post},'b');
title(['Stationary -HTP, DART',' n = ', num2str(length(statGreen))])
ylabel('dF/F, pref dir') 
xlabel('size (deg)') 
set(gca, 'TickDir', 'out')
box off
ylim([-.05 .2])

subplot(2,2,2) %for the second day
errorbar(sizes,sizeResp_red_avrg_stat{pre},sizeResp_red_se_stat{pre},'k');
hold on
errorbar(sizes,sizeResp_red_avrg_stat{post},sizeResp_red_se_stat{post},'b');
title(['Stationary +HTP',' n = ', num2str(length(statRed))])
ylabel('dF/F, pref dir') 
xlabel('size (deg)') 
set(gca, 'TickDir', 'out')
box off
ylim([-.05 .2])

subplot(2,2,3) %for the first day
errorbar(sizes,sizeResp_green_avrg_loc{pre},sizeResp_green_se_loc{pre},'k');
hold on
errorbar(sizes,sizeResp_green_avrg_loc{post},sizeResp_green_se_loc{post},'b');
title(['Running -HTP',' n = ', num2str(length(runningGreen))])
ylabel('dF/F, pref dir') 
xlabel('size (deg)') 
set(gca, 'TickDir', 'out')
box off
ylim([-.05 .3])

subplot(2,2,4) %for the second day
errorbar(sizes,sizeResp_red_avrg_loc{pre},sizeResp_red_se_loc{pre},'k');
hold on
errorbar(sizes,sizeResp_red_avrg_loc{post},sizeResp_red_se_loc{post},'b');
title(['Running +HTP',' n = ', num2str(length(runningRed))])
ylabel('dF/F, pref dir') 
xlabel('size (deg)') 
set(gca, 'TickDir', 'out')
box off
ylim([-.05 .3])

x0=5;
y0=5;
width=6;
height=6;
set(gcf,'units','inches','position',[x0,y0,width,height])



sgtitle(['population size tuning' ])

print(fullfile(fnout,['sizeTuning.pdf']),'-dpdf');

%% plotting size tuning for running vs. stationary, seperated by drug condition

figure
subplot(2,2,1) %for the first day
errorbar(sizes,sizeResp_green_avrg_stat{pre},sizeResp_green_se_stat{pre},'k');
hold on
errorbar(sizes,sizeResp_green_avrg_loc{pre},sizeResp_green_se_loc{pre},'m');
title(['Pre-DART -HTP',' n = ', num2str(length(runningGreen))])
ylabel('dF/F, pref dir') 
xlabel('size (deg)') 
set(gca, 'TickDir', 'out')
box off

subplot(2,2,2) %for the second day
errorbar(sizes,sizeResp_red_avrg_stat{pre},sizeResp_red_se_stat{pre},'k');
hold on
errorbar(sizes,sizeResp_red_avrg_loc{pre},sizeResp_red_se_loc{pre},'m');
title(['Pre-DART +HTP',' n = ', num2str(length(runningRed))])
ylabel('dF/F, pref dir') 
xlabel('size (deg)') 
set(gca, 'TickDir', 'out')
box off

subplot(2,2,3) %for the first day

errorbar(sizes,sizeResp_green_avrg_stat{post},sizeResp_green_se_stat{post},'k');
hold on
errorbar(sizes,sizeResp_green_avrg_loc{post},sizeResp_green_se_loc{post},'m');
title(['Post-DART -HTP',' n = ', num2str(length(runningGreen))])
ylabel('dF/F, pref dir') 
xlabel('size (deg)') 
set(gca, 'TickDir', 'out')
box off

subplot(2,2,4) %for the second day
errorbar(sizes,sizeResp_red_avrg_stat{post},sizeResp_red_se_stat{post},'k');
hold on
errorbar(sizes,sizeResp_red_avrg_loc{post},sizeResp_red_se_loc{post},'m');
title(['Post-DART +HTP',' n = ', num2str(length(runningRed))])
ylabel('dF/F, pref dir') 
xlabel('size (deg)') 
set(gca, 'TickDir', 'out')
box off

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])



sgtitle(['population size tuning' ])

print(fullfile(fnout,['sizeTuningVsBehState.pdf']),'-dpdf');

%% contrast response
%errorbar for stat resp and loc resp vs size, where error is across mice
conResp_green_avrg_stat = cell(nSize,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_stat = cell(nSize,nd); %same for red
conResp_green_se_stat = cell(nSize,nd); %this will be the se across all green cells
conResp_red_se_stat = cell(nSize,nd); %same for red

consForPlotting = [12.5 25 50 100];
for id = 1:nd
    for iSize = 1:nSize
        green_data=conBySize_resp_stat_concat{id}(green_ind_concat,:,iSize);%pulling the green cells at this size
        conResp_green_avrg_stat{id}(iSize,:)=nanmean(green_data,1);
        green_std=nanstd(green_data,1);
        conResp_green_se_stat{id}(iSize,:)=green_std/sqrt(length(statGreen));
        
        red_data=conBySize_resp_stat_concat{id}(red_ind_concat,:,iSize);%pulling the red cells at this size
        conResp_red_avrg_stat{id}(iSize,:)=nanmean(red_data,1);
        red_std=nanstd(red_data,1);
        conResp_red_se_stat{id}(iSize,:)=red_std/sqrt(length(statRed));
        
        clear green_std red_std green_data red_data
    end
end


figure
subplot(1,2,1) %for the first day
errorbar(consForPlotting,conResp_green_avrg_stat{pre}(1,:),conResp_green_se_stat{pre}(1,:),'k');
hold on
errorbar(consForPlotting,conResp_green_avrg_stat{post}(1,:),conResp_green_se_stat{post}(1,:),'b');
title('20 deg')
ylabel('dF/F, pref dir') 
set(gca, 'TickDir', 'out')
box off
ylim([-.005 .2])

subplot(1,2,2) %for the first day
errorbar(consForPlotting,conResp_green_avrg_stat{pre}(2,:),conResp_green_se_stat{pre}(2,:),'k');
hold on
errorbar(consForPlotting,conResp_green_avrg_stat{post}(2,:),conResp_green_se_stat{post}(2,:),'b');
title('fullfield')
ylabel('dF/F, pref dir') 
set(gca, 'TickDir', 'out')
box off
ylim([-.005 .2])

han=axes('visible','off'); 
han.XLabel.Visible='on';
xlabel(han,'Contrast (%)');


% 
% subplot(2,2,3) %for the first day
% errorbar(sizes,sizeResp_green_avrg_loc{pre},sizeResp_green_se_loc{pre},'k');
% hold on
% errorbar(sizes,sizeResp_green_avrg_loc{post},sizeResp_green_se_loc{post},'b');
% title(['Running -HTP',' n = ', num2str(length(runningGreen))])
% ylabel('dF/F, pref dir') 
% xlabel('size (deg)') 
% set(gca, 'TickDir', 'out')
% box off
% ylim([-.05 .3])
% 
% subplot(2,2,4) %for the second day
% errorbar(sizes,sizeResp_red_avrg_loc{pre},sizeResp_red_se_loc{pre},'k');
% hold on
% errorbar(sizes,sizeResp_red_avrg_loc{post},sizeResp_red_se_loc{post},'b');
% title(['Running +HTP',' n = ', num2str(length(runningRed))])
% ylabel('dF/F, pref dir') 
% xlabel('size (deg)') 
% set(gca, 'TickDir', 'out')
% box off
% ylim([-.05 .3])
% 
x0=5;
y0=5;
width=6;
height=3.5;
set(gcf,'units','inches','position',[x0,y0,width,height])

sgtitle(['Stationary',' n = ', num2str(length(statGreen))])

print(fullfile(fnout,['contrastTuning.pdf']),'-dpdf');

