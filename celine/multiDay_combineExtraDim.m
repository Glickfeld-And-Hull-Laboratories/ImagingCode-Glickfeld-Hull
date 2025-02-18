
clear all; clear global; close all
clc
%ds = 'DART_V1_contrast_ori_Celine'; %dataset info
ds = 'DART_V1_atropine_Celine';
dataStructLabels = {'contrastxori'};

eval(ds);

experimentFolder = 'SST_atropine';
sess_list = [4,8,12,44,46,55];%enter all the sessions you want to concatenate
nSess=length(sess_list);

oldNewCutoff = 4; %list "old" session (3 con, 1 size, 8 dir) first, and indicate here where the first "new" session (4 con, 2 sze, 4 dir) is

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

% INDICATE THE PRE VS. POST DAYS DEPENDING ON THE ORDER OF MATCHING
prompt = 'EMX mice? 0- no, 1- yes';
            isEMX = input(prompt);
            switch isEMX
                case 0
                    pre=2; %baeline session, used as reference, is in the 1st position
                    post=1;
                    "baseline used as reference"
                case 1
                  pre=1;
                  post=2; %post-DART session, used as reference, is in the 1st position  
                  "post-DART used as reference"
            end
clear x prompt

targetCon = [.25 .5 1]%what contrast to extract for all data - must be one that all datasets had


frame_rate = 15;

sess_title = string(sess_list(1));
for iSess = 2:nSess
    sess_title = strcat(sess_title,'_',string(sess_list(iSess)));
end
d=string(datetime('today'));

if computer == 'GLNXA64'
    isilonName =  '/home/cc735@dhe.duke.edu/GlickfeldLabShare';
    base = fullfile(isilonName, '\All_Staff\home\ACh\Analysis\2p_analysis',experimentFolder);
    
else
    isilonName = 'Z:';
    base = fullfile(isilonName, '\home\ACh\Analysis\2p_analysis',experimentFolder);
      
end


if nSess == 1
         if expt(sess_list(1)).multiday_timesincedrug_hours>0
            dart_str = [expt(sess_list(1)).drug '_' num2str(expt(sess_list(1)).multiday_timesincedrug_hours) 'Hr'];
        else
            dart_str = 'control';
        end
        
        fnout = fullfile(base,expt(sess_list(1)).mouse,['multiday_' dart_str],d);
else
    fnout= fullfile(base,strcat('concat', sess_title),d);
end
mkdir(fnout);
cd(fnout)
clear d sess_title

% concatenating data
nCon = length(targetCon)

mice={};
tc_trial_avrg_stat_concat=cell(1,nd);
tc_trial_avrg_stat_largePupil_concat=cell(1,nd);
tc_trial_avrg_stat_smallPupil_concat=cell(1,nd);
tc_trial_avrg_loc_concat=cell(1,nd);
pref_responses_loc_concat=cell(1,nd);
pref_responses_stat_concat=cell(1,nd);
pref_responses_stat_largePupil_concat=cell(1,nd);
pref_responses_stat_smallPupil_concat=cell(1,nd);
dirs_concat=[];
cons_concat=[];
green_concat=[];
red_concat=[];
nKeep_concat=[];
norm_dir_resp_stat_concat = cell(1,nd);
norm_dir_resp_loc_concat = cell(1,nd);
pref_dir_concat=cell(1,nd);
noiseCorr_concat = cell(4,nd);

pref_allTrials_stat_concat =cell(nCon,nd);
pref_allTrials_loc_concat =cell(nCon,nd);

%for the "old" sessions

for iSess = 1:oldNewCutoff-1
    targetSize=1;
    day_id = sess_list(iSess)
    mouse = expt(day_id).mouse;
    mice=[mice;mouse];




    if expt(day_id).multiday_timesincedrug_hours>0
        dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
    else
        dart_str = 'control';
    end
    fn_multi = fullfile(base,mouse,['multiday_' dart_str]);

    load(fullfile(fn_multi,'tc_keep.mat'));
    load(fullfile(fn_multi,'resp_keep.mat'));
    load(fullfile(fn_multi,'input.mat'));
    load(fullfile(fn_multi,'locomotion.mat'));
    load(fullfile(fn_multi,'fluor_intensity.mat'));
    load(fullfile(fn_multi,'HT_pyr_relationship.mat'));


    nKeep = size(tc_trial_avrg_stat{post},2);

    
%    tells the contrast, direction and orientation for each trial each day
    tCon_match = cell(1,nd);
    tDir_match = cell(1,nd);
    tOri_match = cell(1,nd);

 %   find the contrasts, directions and orientations for each day
    for id = 1:nd
        tCon_match{id} = celleqel2mat_padded(input(id).tGratingContrast(1:nTrials(id)));
        tDir_match{id} = celleqel2mat_padded(input(id).tGratingDirectionDeg(1:nTrials(id)));
        tOri_match{id} = tDir_match{id};
        tOri_match{id}(find(tDir_match{id}>=180)) = tDir_match{id}(find(tDir_match{id}>=180))-180;
    end
    dirs = unique(tDir_match{post});
    cons = unique(tCon_match{post});
    sharedCon=find(ismember(cons, targetCon));

    nOn = input(1).nScansOn;
    nOff = input(1).nScansOff;
  %  start conatenating
    dirs_concat = [dirs_concat,dirs]; 
    cons_concat = [cons_concat,cons(sharedCon)];
    red_concat = [red_concat, red_keep_logical];
    green_concat = [green_concat, green_keep_logical];
    nKeep_concat = [nKeep_concat,nKeep];

    clear cons
        
    for id = 1:nd
        
        tc_trial_avrg_stat_concat{id} =cat(2,tc_trial_avrg_stat_concat{id},tc_trial_avrg_stat{id}(:,:,sharedCon,targetSize));
        tc_trial_avrg_stat_largePupil_concat{id} = cat(2,tc_trial_avrg_stat_largePupil_concat{id},tc_trial_avrg_stat_largePupil{id}(:,:,sharedCon,targetSize));
        tc_trial_avrg_stat_smallPupil_concat{id} = cat(2,tc_trial_avrg_stat_smallPupil_concat{id},tc_trial_avrg_stat_smallPupil{id}(:,:,sharedCon,targetSize));
        tc_trial_avrg_loc_concat{id} =cat(2,tc_trial_avrg_loc_concat{id},tc_trial_avrg_loc{id}(:,:,sharedCon,targetSize));

        pref_responses_loc_concat{id}=cat(1,pref_responses_loc_concat{id},pref_responses_loc{id}(:,sharedCon,targetSize));
        pref_responses_stat_concat{id}=cat(1,pref_responses_stat_concat{id},pref_responses_stat{id}(:,sharedCon,targetSize));
        pref_responses_stat_largePupil_concat{id}=cat(1,pref_responses_stat_largePupil_concat{id},pref_responses_stat_largePupil{id}(:,sharedCon,targetSize));
        pref_responses_stat_smallPupil_concat{id}=cat(1,pref_responses_stat_smallPupil_concat{id},pref_responses_stat_smallPupil{id}(:,sharedCon,targetSize));
        
        norm_dir_resp_stat_concat{id}=cat(1,norm_dir_resp_stat_concat{id},norm_dir_resp_stat{id}(:,1:4,:,targetSize));
        norm_dir_resp_loc_concat{id}=cat(1,norm_dir_resp_loc_concat{id},norm_dir_resp_loc{id}(:,1:4,:,targetSize));
 
        pref_dir_concat{id}=cat(2,pref_dir_concat{id},pref_dir_keep{id});
        noiseCorr_concat{id}=cat(2,noiseCorr_concat{id},noiseCorr_OG{id});


        for i = 1:length(sharedCon)
            iCon=sharedCon(i);
            pref_allTrials_stat_concat{i,id}=[pref_allTrials_stat_concat{i,id},pref_allTrials_stat{iCon,id}];
            pref_allTrials_loc_concat{i,id}=[pref_allTrials_loc_concat{i,id},pref_allTrials_loc{iCon,id}];
        end
        clear meanF i
    end

   
end

%for the "new" sessions

for iSess = oldNewCutoff:nSess
    targetSize = 2; %this is the DIMENSION - ie for matrices with two sizes, enter 1 to get the smaller size or 2 to get the larger size

    day_id = sess_list(iSess)
    mouse = expt(day_id).mouse;
    mice=[mice;mouse];




    if expt(day_id).multiday_timesincedrug_hours>0
        dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
    else
        dart_str = 'control';
    end
    fn_multi = fullfile(base,mouse,['multiday_' dart_str]);

   load(fullfile(fn_multi,'tc_keep.mat'));
    load(fullfile(fn_multi,'resp_keep.mat'));
    load(fullfile(fn_multi,'input.mat'));
    load(fullfile(fn_multi,'locomotion.mat'));
    load(fullfile(fn_multi,'fluor_intensity.mat'));
   load(fullfile(fn_multi,'HT_pyr_relationship.mat'));


    nKeep = size(tc_trial_avrg_stat{post},2);

    
%    tells the contrast, direction and orientation for each trial each day
    tCon_match = cell(1,nd);
    tDir_match = cell(1,nd);
    tOri_match = cell(1,nd);

 %   find the contrasts, directions and orientations for each day
    for id = 1:nd
        tCon_match{id} = celleqel2mat_padded(input(id).tGratingContrast(1:nTrials(id)));
        tDir_match{id} = celleqel2mat_padded(input(id).tGratingDirectionDeg(1:nTrials(id)));
        tOri_match{id} = tDir_match{id};
        tOri_match{id}(find(tDir_match{id}>=180)) = tDir_match{id}(find(tDir_match{id}>=180))-180;
    end
    dirs = unique(tDir_match{post});
    cons = unique(tCon_match{post});
    sharedCon=find(ismember(cons, targetCon));

    nOn = input(1).nScansOn;
    nOff = input(1).nScansOff;
  %  start conatenating
    dirs_concat = [dirs_concat,dirs]; 
    cons_concat = [cons_concat,cons(sharedCon)];
    red_concat = [red_concat, red_keep_logical];
    green_concat = [green_concat, green_keep_logical];
    nKeep_concat = [nKeep_concat,nKeep];

    clear cons
        
    for id = 1:nd
        
        tc_trial_avrg_stat_concat{id} =cat(2,tc_trial_avrg_stat_concat{id},tc_trial_avrg_stat{id}(:,:,sharedCon,targetSize));
        % tc_trial_avrg_stat_largePupil_concat{id} = cat(2,tc_trial_avrg_stat_largePupil_concat{id},tc_trial_avrg_stat_largePupil{id}(:,:,sharedCon,targetSize));
        % tc_trial_avrg_stat_smallPupil_concat{id} = cat(2,tc_trial_avrg_stat_smallPupil_concat{id},tc_trial_avrg_stat_smallPupil{id}(:,:,sharedCon,targetSize));
        % tc_trial_avrg_loc_concat{id} =cat(2,tc_trial_avrg_loc_concat{id},tc_trial_avrg_loc{id}(:,:,sharedCon,targetSize));

        pref_responses_loc_concat{id}=cat(1,pref_responses_loc_concat{id},pref_responses_loc{id}(:,sharedCon,targetSize));
        pref_responses_stat_concat{id}=cat(1,pref_responses_stat_concat{id},pref_responses_stat{id}(:,sharedCon,targetSize));
        % pref_responses_stat_largePupil_concat{id}=cat(1,pref_responses_stat_largePupil_concat{id},pref_responses_stat_largePupil{id}(:,sharedCon,targetSize));
        % pref_responses_stat_smallPupil_concat{id}=cat(1,pref_responses_stat_smallPupil_concat{id},pref_responses_stat_smallPupil{id}(:,sharedCon,targetSize));
        % 
        norm_dir_resp_stat_concat{id}=cat(1,norm_dir_resp_stat_concat{id},norm_dir_resp_stat{id}(:,:,sharedCon,targetSize));
        norm_dir_resp_loc_concat{id}=cat(1,norm_dir_resp_loc_concat{id},norm_dir_resp_loc{id}(:,:,sharedCon,targetSize));
 
        pref_dir_concat{id}=cat(2,pref_dir_concat{id},pref_dir_keep{id});
        noiseCorr_concat{id}=cat(2,noiseCorr_concat{id},noiseCorr{id});

        for i = 1:length(sharedCon)
            iCon=sharedCon(i);
            pref_allTrials_stat_concat{i,id}=[pref_allTrials_stat_concat{i,id},pref_allTrials_stat{iCon,targetSize,id}];
            pref_allTrials_loc_concat{i,id}=[pref_allTrials_loc_concat{i,id},pref_allTrials_loc{iCon,targetSize,id}];
        end

        clear meanF i
    end

   
end
%
clear mouse day_id nKeep iSess fn_multi cons oris
clear explanation1 resp_keep tc_trial_avrg_keep_allCond pref_responses_allCond sig_diff pref_con_keep pref_ori_keep tOri_match tCon_match data_trial_keep nTrials tc_trial_avrg_keep green_keep_logical red_keep_logical green_ind_keep red_ind_keep
clear LMI RIx locCounts locResp locTCs statResp statTCs wheel_tc ttest_results_stat ttest_results_loc ttest_results_allCon_stat ttest_results_allCon_loc
clear data_con_resp_keep data_ori_resp_keep data_rep_keep dfof_max_diff dfof_max_diff_raw explanation2 resp_max_keep data_resp_keep pref_responses_stat pref_responses_loc
clear tc_trial_avrg_stat tc_trial_avrg_loc fullTC_keep norm_dir_resp sigCorr noiseCorr
clear red_fluor_all red_fluor_match green_fluor_match green_fluor_match red_fluor_keep green_fluor_keep R_p_values
red_ind_concat = find(red_concat);
green_ind_concat = find(green_concat);
%
cons=unique(cons_concat);

dirs=unique(dirs_concat);
nDir=length(dirs);
nKeep_total = sum(nKeep_concat);


%find cells that I have running data for on both days
haveRunning_pre = ~isnan(pref_responses_loc_concat{pre});

haveRunning_post = ~isnan(pref_responses_loc_concat{post});
for iCon =1:nCon
    haveRunning_both{iCon}= find(haveRunning_pre(:,iCon).* haveRunning_post(:,iCon));
    haveRunning_green{iCon} = intersect(haveRunning_both{iCon}, green_ind_concat);
    haveRunning_red{iCon} = intersect(haveRunning_both{iCon}, red_ind_concat);
end
%

clear haveRunning_pre haveRunning_post haveRunning_both

% %to swap green and red label for EMX experiments
if isEMX==1
    for iCon = 1:nCon
        haveRunning_green{iCon}=haveRunning_red{iCon};
        haveRunning_red{iCon}=[];
    end
    green_ind_concat = red_ind_concat;
    red_ind_concat=[];
end

% %for instances when I know I don't have enough running trials, juse use
% all keep cells.
% 
% for iCon =1:nCon
% 
% haveRunning_green{iCon} = green_ind_concat;
% haveRunning_red{iCon} = red_ind_concat;
% end

%get an array with indices for each mouse
mouseInds=cell(1,nSess);
start=1;
for iMouse = 1:nSess
    mouseInds{iMouse}=start:(start-1)+nKeep_concat(iMouse);
    start = start+nKeep_concat(iMouse)
end
clear start iMouse


%find how many haveRunning red cells I got for each mouse
cellCountsRed = nan(nSess,nCon);
mouseNames=[];
for iMouse = 1:nSess
    for iCon = 1:nCon
        cellCountsRed(iMouse, iCon,1)=length(intersect(haveRunning_red{iCon},(mouseInds{iMouse})));
        
    end
    mouseNames=[mouseNames, string(mice(iMouse,:))]
end


cellCountsGreen = nan(nSess,nCon);
mouseNames=[];
for iMouse = 1:nSess
    for iCon = 1:nCon
        cellCountsGreen(iMouse, iCon,1)=length(intersect(haveRunning_green{iCon},(mouseInds{iMouse})));
        
    end
    mouseNames=[mouseNames, string(mice(iMouse,:))]
end


cellCountTableRed = table(cellCountsRed, RowNames=mouseNames)
cellCountTableGreen = table(cellCountsGreen, RowNames=mouseNames)
writetable(cellCountTableRed,fullfile(fnout,['cellCounts.csv']),'WriteRowNames',true)
writetable(cellCountTableGreen,fullfile(fnout,['cellCounts_Green.csv']),'WriteRowNames',true)
clear cellCountsRed cellCountsGreen

green_all = intersect(haveRunning_green{1},haveRunning_green{2});
green_all = intersect(green_all, haveRunning_green{3});

red_all = intersect(haveRunning_red{1},haveRunning_red{2});
red_all = intersect(red_all, haveRunning_red{3});



%find cells that have running data on both days
have_HI_pre = ~isnan(pref_responses_stat_largePupil_concat{pre});
have_HI_post = ~isnan(pref_responses_stat_largePupil_concat{post});

have_LOW_pre = ~isnan(pref_responses_stat_smallPupil_concat{pre});
have_LOW_post = ~isnan(pref_responses_stat_smallPupil_concat{post});

have_bothPupil=cell(1,3);
for iCon =1:nCon
    have_HI_both= find(have_HI_pre(:,iCon).* have_HI_post(:,iCon));
    have_LOW_both=find(have_LOW_pre(:,iCon).* have_LOW_post(:,iCon));
    have_bothPupil{iCon}=intersect(have_HI_both,have_LOW_both);
end

clear have_HI_pre have_HI_post have_LOW_pre have_LOW_post have_HI_both have_LOW_both

have_allPupil = intersect(have_bothPupil{1},have_bothPupil{2});
have_allPupil = intersect(have_allPupil, have_bothPupil{3});

% make dfof summary table for statistics
mouseID=[];
for imouse =1:nSess
   ID_string=mouseNames(imouse);
   thisID = repelem(ID_string,nKeep_concat(imouse));
   mouseID=[mouseID,thisID];
end
mouseID=mouseID';

%% ALL KEEP CELLS

tc_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_stat = cell(1,nd); %same for red
tc_green_se_stat = cell(1,nd); %this will be the se across all green cells
tc_red_se_stat = cell(1,nd); %same for red



for id = 1:nd
    for iCon=1:nCon
        
    tc_green_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,green_ind_concat,iCon),2);
    green_std=nanstd(tc_trial_avrg_stat_concat{id}(:,green_ind_concat,iCon),[],2);
    tc_green_se_stat{id}(:,iCon)=green_std/sqrt(length(green_ind_concat));
    
    tc_red_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,red_ind_concat,iCon),2);
    red_std=nanstd(tc_trial_avrg_stat_concat{id}(:,red_ind_concat,iCon),[],2);
    tc_red_se_stat{id}(:,iCon)=red_std/sqrt(length(red_ind_concat));
    
    clear green_std red_std
    end
end

%create a 2-second stimulus marker
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);
positions=[1,2;3,4;5,6];
figure
for iCon = 1:nCon
p1=positions(iCon,1);
p2=positions(iCon,2);

subplot(3,2,p1) 
shadedErrorBar(t,tc_red_avrg_stat{pre}(:,iCon),tc_red_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_stat{post}(:,iCon),tc_red_se_stat{post}(:,iCon),'b');
ylim([-.02 .17]);
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);

if iCon==1
    title(['SST',' n = ', num2str(length(red_ind_concat))])
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
end
set(gca,'XColor', 'none','YColor','none')

subplot(3,2,p2) 
ylim([-.02 .25]);
hold on
shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon),tc_green_se_stat{pre}(:,iCon),'--k');
hold on
shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon),tc_green_se_stat{post}(:,iCon),'--b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2)

if iCon==1
    title(['Pyr',' n = ', num2str(length(green_ind_concat))])
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
end
set(gca,'XColor', 'none','YColor','none')

x0=5;
y0=0;
width=4;
height=9;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')


end  

print(fullfile(fnout,'_stat_allKeep_timecourses.pdf'),'-dpdf');

%% plot large pupil TCs

tc_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_stat = cell(1,nd); %same for red
tc_green_se_stat = cell(1,nd); %this will be the se across all green cells
tc_red_se_stat = cell(1,nd); %same for red



for id = 1:nd
    for iCon=1:nCon
    temp_green = intersect(green_ind_concat, have_bothPupil{iCon});
    temp_red = intersect(red_ind_concat, have_bothPupil{iCon});
    tc_green_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_largePupil_concat{id}(:,temp_green,iCon),2);
    green_std=nanstd(tc_trial_avrg_stat_largePupil_concat{id}(:,temp_green,iCon),[],2);
    tc_green_se_stat{id}(:,iCon)=green_std/sqrt(length(temp_green));
    
    tc_red_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_largePupil_concat{id}(:,temp_red,iCon),2);
    red_std=nanstd(tc_trial_avrg_stat_largePupil_concat{id}(:,temp_red,iCon),[],2);
    tc_red_se_stat{id}(:,iCon)=red_std/sqrt(length(temp_red));
    
    clear green_std red_std
    end
end
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(1,2,2) %for the first day

ylim([-.02 .25]);;
hold on
shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon),tc_green_se_stat{pre}(:,iCon),'--k');
hold on
shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon),tc_green_se_stat{post}(:,iCon),'--b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['Pyr',' n = ', num2str(length(temp_green))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,1) %for the second day
shadedErrorBar(t,tc_red_avrg_stat{pre}(:,iCon),tc_red_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_stat{post}(:,iCon),tc_red_se_stat{post}(:,iCon),'b');
ylim([-.02 .25]);;
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['SST',' n = ', num2str(length(temp_red))])

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['Large pupil stationary, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) '_stat_largePupil_timecourses.pdf']),'-dpdf');
end  

%% plot small pupil TCs

tc_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_stat = cell(1,nd); %same for red
tc_green_se_stat = cell(1,nd); %this will be the se across all green cells
tc_red_se_stat = cell(1,nd); %same for red



for id = 1:nd
    for iCon=1:nCon
    temp_green = intersect(green_ind_concat, have_bothPupil{iCon});
    temp_red = intersect(red_ind_concat, have_bothPupil{iCon});
    tc_green_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_smallPupil_concat{id}(:,temp_green,iCon),2);
    green_std=nanstd(tc_trial_avrg_stat_smallPupil_concat{id}(:,temp_green,iCon),[],2);
    tc_green_se_stat{id}(:,iCon)=green_std/sqrt(length(temp_green));
    
    tc_red_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_smallPupil_concat{id}(:,temp_red,iCon),2);
    red_std=nanstd(tc_trial_avrg_stat_smallPupil_concat{id}(:,temp_red,iCon),[],2);
    tc_red_se_stat{id}(:,iCon)=red_std/sqrt(length(temp_red));
    
    clear green_std red_std
    end
end
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(1,2,2) %for the first day

ylim([-.02 .25]);;
hold on
shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon),tc_green_se_stat{pre}(:,iCon),'--k');
hold on
shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon),tc_green_se_stat{post}(:,iCon),'--b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['Pyr',' n = ', num2str(length(temp_green))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,1) %for the second day
shadedErrorBar(t,tc_red_avrg_stat{pre}(:,iCon),tc_red_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_stat{post}(:,iCon),tc_red_se_stat{post}(:,iCon),'b');
ylim([-.02 .25]);;
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['SST',' n = ', num2str(length(temp_red))])

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['Small pupil stationary, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) '_stat_smallPupil_timecourses.pdf']),'-dpdf');
end

% for the cells that have stationary and running

% make figure with se shaded, one figure per contrast - stationary

%% Stationary timecourses for matched cells

tc_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_stat = cell(1,nd); %same for red
tc_green_se_stat = cell(1,nd); %this will be the se across all green cells
tc_red_se_stat = cell(1,nd); %same for red



for id = 1:nd
    for iCon=1:nCon
        
    tc_green_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,haveRunning_green{iCon},iCon),2);
    green_std=nanstd(tc_trial_avrg_stat_concat{id}(:,haveRunning_green{iCon},iCon),[],2);
    tc_green_se_stat{id}(:,iCon)=green_std/sqrt(length(haveRunning_green{iCon}));
    
    tc_red_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,haveRunning_red{iCon},iCon),2);
    red_std=nanstd(tc_trial_avrg_stat_concat{id}(:,haveRunning_red{iCon},iCon),[],2);
    tc_red_se_stat{id}(:,iCon)=red_std/sqrt(length(haveRunning_red{iCon}));
    
    clear green_std red_std
    end
end


%create a 2-second stimulus marker
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);
positions=[1,2;3,4;5,6];
figure
for iCon = 1:nCon
p1=positions(iCon,1);
p2=positions(iCon,2);

subplot(3,2,p1) 
shadedErrorBar(t,tc_red_avrg_stat{pre}(:,iCon),tc_red_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_stat{post}(:,iCon),tc_red_se_stat{post}(:,iCon),'b');
ylim([-.02 .25]);
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);

title(['SST',' n = ', num2str(length(haveRunning_red{iCon}))])
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);

set(gca,'XColor', 'none','YColor','none')

subplot(3,2,p2) 
ylim([-.02 .25]);
hold on
shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon),tc_green_se_stat{pre}(:,iCon),'--k');
hold on
shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon),tc_green_se_stat{post}(:,iCon),'--b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2)

title(['Pyr',' n = ', num2str(length(haveRunning_green{iCon}))])
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);

set(gca,'XColor', 'none','YColor','none')

x0=5;
y0=0;
width=4;
height=9;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')
sgtitle('stationary')

end  
print(fullfile(fnout,['stat_timecourses.pdf']),'-dpdf');




%% shaded error tc for loc

tc_green_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_loc = cell(1,nd); %same for red
tc_green_se_loc = cell(1,nd); %this will be the se across all green cells
tc_red_se_loc = cell(1,nd); %same for red

for id = 1:nd
    for iCon=1:nCon
    tc_green_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,haveRunning_green{iCon},iCon),2);
    green_std=nanstd(tc_trial_avrg_loc_concat{id}(:,haveRunning_green{iCon},iCon),[],2);
    tc_green_se_loc{id}(:,iCon)=green_std/sqrt(length(haveRunning_green{iCon}));
    
    tc_red_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,haveRunning_red{iCon},iCon),2);
    red_std=nanstd(tc_trial_avrg_loc_concat{id}(:,haveRunning_red{iCon},iCon),[],2);
    tc_red_se_loc{id}(:,iCon)=red_std/sqrt(length(haveRunning_red{iCon}));
    
    clear green_std red_std
    end
end


%create a 2-second stimulus marker
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);
positions=[1,2;3,4;5,6];
figure
for iCon = 1:nCon
p1=positions(iCon,1);
p2=positions(iCon,2);

subplot(3,2,p1) 
shadedErrorBar(t,tc_red_avrg_loc{pre}(:,iCon),tc_red_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_loc{post}(:,iCon),tc_red_se_loc{post}(:,iCon),'b');
ylim([-.02 .25]);
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);

title(['SST',' n = ', num2str(length(haveRunning_red{iCon}))])
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);

set(gca,'XColor', 'none','YColor','none')

subplot(3,2,p2) 
ylim([-.02 .35]);
hold on
shadedErrorBar(t,tc_green_avrg_loc{pre}(:,iCon),tc_green_se_loc{pre}(:,iCon),'--k');
hold on
shadedErrorBar(t,tc_green_avrg_loc{post}(:,iCon),tc_green_se_loc{post}(:,iCon),'--b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2)

title(['Pyr',' n = ', num2str(length(haveRunning_green{iCon}))])
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);

set(gca,'XColor', 'none','YColor','none')

x0=5;
y0=0;
width=4;
height=9;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle('running')
end
print(fullfile(fnout,['loc_timecourses.pdf']),'-dpdf');

clear txt1 txt2


%% scatterplot of max df/f for day 1 vs day 2
green_ex_list=[]; %to highlight particular cells
red_ex_list=[];

for iCon = 1:nCon
figure; movegui('center') 
subplot(2,2,2)
scatter((pref_responses_stat_concat{pre}(haveRunning_green{iCon},iCon)),(pref_responses_stat_concat{post}(haveRunning_green{iCon},iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
mean_pre_stat = nanmean(pref_responses_stat_concat{pre}(haveRunning_green{iCon},iCon));
mean_post_stat = nanmean(pref_responses_stat_concat{post}(haveRunning_green{iCon},iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_responses_stat_concat{pre}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
stderror_post= std(pref_responses_stat_concat{post}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
limMin=min(min(pref_responses_stat_concat{pre}(haveRunning_green{iCon},iCon)),min(pref_responses_stat_concat{post}(haveRunning_green{iCon},iCon)));
limMax=max(max(pref_responses_stat_concat{pre}(haveRunning_green{iCon},iCon)),max(pref_responses_stat_concat{post}(haveRunning_green{iCon},iCon)));
ylim([limMin limMax])
xlim([limMin limMax])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('-HTP stationary')
axis square
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
hold off


subplot(2,2,1)
scatter((pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon)),(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
mean_pre_stat = nanmean(pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon));
mean_post_stat = nanmean(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
stderror_post= std(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);% ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
limMin=min(min(pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon)),min(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon)));
limMax=max(max(pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon)),max(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon)));
ylim([limMin limMax])
xlim([limMin limMax])

hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
title('+HTP stationary')
axis square
hold off

subplot(2,2,4)
scatter((pref_responses_loc_concat{pre}(haveRunning_green{iCon},iCon)),(pref_responses_loc_concat{post}(haveRunning_green{iCon},iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
mean_pre_stat = nanmean(pref_responses_loc_concat{pre}(haveRunning_green{iCon},iCon));
mean_post_stat = nanmean(pref_responses_loc_concat{post}(haveRunning_green{iCon},iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_responses_loc_concat{pre}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
stderror_post= std(pref_responses_loc_concat{post}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
limMin=min(min(pref_responses_loc_concat{pre}(haveRunning_green{iCon},iCon)),min(pref_responses_loc_concat{post}(haveRunning_green{iCon},iCon)));
limMax=max(max(pref_responses_loc_concat{pre}(haveRunning_green{iCon},iCon)),max(pref_responses_loc_concat{post}(haveRunning_green{iCon},iCon)));
ylim([limMin limMax])
xlim([limMin limMax])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('-HTP running')
axis square
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
hold off

subplot(2,2,3)
scatter((pref_responses_loc_concat{pre}(haveRunning_red{iCon},iCon)),(pref_responses_loc_concat{post}(haveRunning_red{iCon},iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
mean_pre_stat = nanmean(pref_responses_loc_concat{pre}(haveRunning_red{iCon},iCon));
mean_post_stat = nanmean(pref_responses_loc_concat{post}(haveRunning_red{iCon},iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
stderror_post= std(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);
xlabel('pre-DART  dF/F')
limMin=min(min(pref_responses_loc_concat{pre}(haveRunning_red{iCon},iCon)),min(pref_responses_loc_concat{post}(haveRunning_red{iCon},iCon)));
limMax=max(max(pref_responses_loc_concat{pre}(haveRunning_red{iCon},iCon)),max(pref_responses_loc_concat{post}(haveRunning_red{iCon},iCon)));
ylim([limMin limMax])
xlim([limMin limMax])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('+HTP running')
uistack(hline,'bottom');
axis square
hold off
set(gca, 'TickDir', 'out')



sgtitle(num2str(cons(iCon)))
print(fullfile(fnout,[num2str(cons(iCon)) 'maxResp_crossDay.pdf']),'-dpdf','-bestfit')
clear mean_pre_stat mean_post_stat stderror_post stderror_pre
end


%% t-test for day within each contrast

[h1, p1]= ttest2(pref_responses_stat_concat{pre}(haveRunning_red{1},1),pref_responses_stat_concat{post}(haveRunning_red{1},1));
[h2,p2]= ttest2(pref_responses_stat_concat{pre}(haveRunning_red{2},2),pref_responses_stat_concat{post}(haveRunning_red{2},2));
[h3,p3]= ttest2(pref_responses_stat_concat{pre}(haveRunning_red{3},3),pref_responses_stat_concat{post}(haveRunning_red{3},3));
[h4,p4]= ttest2(pref_responses_loc_concat{pre}(haveRunning_red{1},1),pref_responses_loc_concat{post}(haveRunning_red{1},1));
[h5,p5]= ttest2(pref_responses_loc_concat{pre}(haveRunning_red{2},2),pref_responses_loc_concat{post}(haveRunning_red{2},2));
[h3,p3]= ttest2(pref_responses_loc_concat{pre}(haveRunning_red{3},3),pref_responses_loc_concat{post}(haveRunning_red{3},3));


%correct for six tests
table([p1*3 p2*3 p3*3],[p4*3 p5*3 p3*3])

clear h1 p1 h2 p2 h3 p3 h4 p4

[h1, p1]= ttest2(pref_responses_stat_concat{pre}(haveRunning_green{1},1),pref_responses_stat_concat{post}(haveRunning_green{1},1));
[h2,p2]= ttest2(pref_responses_stat_concat{pre}(haveRunning_green{2},2),pref_responses_stat_concat{post}(haveRunning_green{2},2));
[h3,p3]= ttest2(pref_responses_stat_concat{pre}(haveRunning_green{3},3),pref_responses_stat_concat{post}(haveRunning_green{3},3));
[h4,p4]= ttest2(pref_responses_loc_concat{pre}(haveRunning_green{1},1),pref_responses_loc_concat{post}(haveRunning_green{1},1));
[h5,p5]= ttest2(pref_responses_loc_concat{pre}(haveRunning_green{2},2),pref_responses_loc_concat{post}(haveRunning_green{2},2));
[h3,p3]= ttest2(pref_responses_loc_concat{pre}(haveRunning_green{3},3),pref_responses_loc_concat{post}(haveRunning_green{3},3));


%correct for three tests
table([p1*3 p2*3 p3*3],[p4*3 p5*3 p3*3])

clear h1 p1 h2 p2 h3 p3 h4 p4
%% split timecourses by correlation
% Identify high and low correlation cells

% cells with high correlation in the baseline day
highRInds = find(noiseCorr_concat{pre}(1,:)>0.5);
lowRInds = find(noiseCorr_concat{pre}(1,:)<=0.5);

redHigh=intersect(highRInds, red_ind_concat);
redLow=intersect(lowRInds, red_ind_concat);

% finding how the high and low R cells are distributed over epxeriments

RbyExp = zeros(2,nSess);
for iSess = 1:nSess
    mouseIndsTemp = mouseInds{iSess};
    RbyExp(1,iSess) = length(intersect(mouseIndsTemp,redHigh));
    RbyExp(2,iSess) = length(intersect(mouseIndsTemp,redLow));
end
array2table(RbyExp,RowNames={'high R'  'low R'})
%% split timecourses by correlation
hi_avrg_stat = cell(1,nd);
low_avrg_stat = cell(1,nd); 
hi_se_stat = cell(1,nd); 
low_se_stat = cell(1,nd);

%find cells with high correlation in the baseline day
highRInds = find(noiseCorr_concat{pre}(1,:)>0.5);
lowRInds = find(noiseCorr_concat{pre}(1,:)<=0.5);

redHigh=cell(nCon,1);
redLow=cell(nCon,1);

for id = 1:nd

   for iCon=1:nCon
       redHigh{iCon} =  intersect(haveRunning_red{iCon},highRInds); %find intersection of red cells with high R 
       % and red cells that I have running data for on both days, for this contrast
        redLow{iCon} =  intersect(haveRunning_red{iCon},lowRInds); %find intersection of red cells with low R 
       % and red cells that I have running data for on both days, for this contrast


        hi_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:, redHigh{iCon},iCon),2);
        high_std=nanstd(tc_trial_avrg_stat_concat{id}(:, redHigh{iCon},iCon),[],2);
        hi_se_stat{id}(:,iCon)=high_std/sqrt(length(redHigh{iCon}));
        
        low_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,redLow{iCon},iCon),2);
        low_std=nanstd(tc_trial_avrg_stat_concat{id}(:, redLow{iCon},iCon),[],2);
        low_se_stat{id}(:,iCon)=low_std/sqrt(length(redLow{iCon}));
        
        clear low_std high_std
    end
end
z=double(nOn)/double(frame_rate);

%creat a time axis in seconds
t=1:(size(hi_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure



subplot(1,2,1) 
shadedErrorBar(t,low_avrg_stat{pre}(:,iCon),low_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,low_avrg_stat{post}(:,iCon),low_se_stat{post}(:,iCon),'b');
ylim([-.02 .5]);
hold on
% line([0,.2],[-.01,-.01],'Color','black','LineWidth',2);
% hold on
line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['Weakly correlated',' n = ', num2str(length(redLow{iCon}))])
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) 
ylim([-.02 .5]);
hold on
shadedErrorBar(t,hi_avrg_stat{pre}(:,iCon),hi_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,hi_avrg_stat{post}(:,iCon),hi_se_stat{post}(:,iCon),'b');
hold on
% line([0,.2],[-.01,-.01],'Color','black','LineWidth',2);
% hold on
line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title([' Strongly correlated',' n = ', num2str(length(redHigh{iCon}))])

xlabel('s') 
set(gca,'XColor', 'none','YColor','none')

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['stationary, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) 'stat_R_timecourses.pdf']),'-dpdf');
clear txt1 highRed lowRed
end 


%% scatterplot of max df/f for day 1 vs day 2 for high and low R red cells

for iCon = 1:nCon
   redHigh =  intersect(haveRunning_red{iCon},highRInds); %find intersection of red cells with high R 
   % and red cells that I have running data for on both days, for this contrast
    redLow =  intersect(haveRunning_red{iCon},lowRInds); %find intersection of red cells with low R 
   % and red cells that I have running data for on both days, for this contrast


figure; movegui('center') 
subplot(1,2,2)
scatter((pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon)),(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
scatter((pref_responses_stat_concat{pre}(redHigh,iCon)),(pref_responses_stat_concat{post}(redHigh,iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on

mean_pre_stat = nanmean(pref_responses_stat_concat{pre}(redHigh,iCon));
mean_post_stat = nanmean(pref_responses_stat_concat{post}(redHigh,iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_responses_stat_concat{pre}(redHigh,iCon)) / sqrt( length(redHigh));
stderror_post= std(pref_responses_stat_concat{post}(redHigh,iCon)) / sqrt( length(redHigh));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);


xlabel('pre-DART  dF/F')
ylim([-.2 0.2])
xlim([-.2 0.2])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('+HTP high r, stationary')
axis square
set(gca, 'TickDir', 'out')



subplot(1,2,1)
scatter((pref_responses_stat_concat{pre}(redLow,iCon)),(pref_responses_stat_concat{post}(redLow,iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
mean_pre_stat = nanmean(pref_responses_stat_concat{pre}(redLow,iCon));
mean_post_stat = nanmean(pref_responses_stat_concat{post}(redLow,iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_responses_stat_concat{pre}(redLow,iCon)) / sqrt( length(redLow));
stderror_post= std(pref_responses_stat_concat{post}(redLow,iCon)) / sqrt( length(redLow));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);% ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylabel('post-DART dF/F')
ylim([-.2 0.2])
xlim([-.2 0.2])
set(gca, 'TickDir', 'out')
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('+HTP low r, stationary')
axis square


sgtitle(num2str(cons(iCon)))
print(fullfile(fnout,[num2str(cons(iCon)) 'stat_maxResp_R.pdf']),'-dpdf','-bestfit')
clear txt1 highRed lowRed
end

[h1, p1]= ttest(pref_responses_stat_concat{pre}(redLow,1),pref_responses_stat_concat{post}(redLow,1));
[h2,p2]= ttest(pref_responses_stat_concat{pre}(redHigh,1),pref_responses_stat_concat{post}(redHigh,1));

%correct for four tests
p1*2
p2*2
clear h1 p1 h2 p2

%% 
hi_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
low_avrg_loc = cell(1,nd); %same for red
hi_se_loc = cell(1,nd); %this will be the se across all green cells
low_se_loc = cell(1,nd); %same for red

redHigh=cell(nCon);
redLow=cell(nCon);
for id = 1:nd
    for iCon=1:nCon
   redHigh{iCon} =  intersect(haveRunning_red{iCon},highRInds); %find intersection of red cells with high R 
   % and red cells that I have running data for on both days, for this contrast
    redLow{iCon} =  intersect(haveRunning_red{iCon},lowRInds); %find intersection of red cells with low R 
   % and red cells that I have running data for on both days, for this contrast


    hi_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:, redHigh{iCon},iCon),2);
    high_std=nanstd(tc_trial_avrg_loc_concat{id}(:, redHigh{iCon},iCon),[],2);
    hi_se_loc{id}(:,iCon)=high_std/sqrt(length(redHigh{iCon}));
    
    low_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,redLow{iCon},iCon),2);
    low_std=nanstd(tc_trial_avrg_loc_concat{id}(:, redLow{iCon},iCon),[],2);
    low_se_loc{id}(:,iCon)=low_std/sqrt(length(redLow{iCon}));
    
    clear low_std high_std
    end
end
z=double(nOn)/double(frame_rate);

%creat a time axis in seconds
t=1:(size(hi_avrg_loc{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure

subplot(1,2,1) %for +HTP
shadedErrorBar(t,low_avrg_loc{pre}(:,iCon),low_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,low_avrg_loc{post}(:,iCon),low_se_loc{post}(:,iCon),'b');
ylim([-.02 .5]);
hold on
line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
%ylabel('dF/F') 
xlabel('s') 
title(['Weakly correlated',' n = ', num2str(length(redLow{iCon}))])
set(gca,'XColor', 'none','YColor','none')

subplot(1,2,2) 
ylim([-.02 .5]);
hold on
shadedErrorBar(t,hi_avrg_loc{pre}(:,iCon),hi_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,hi_avrg_loc{post}(:,iCon),hi_se_loc{post}(:,iCon),'b');
line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['Strongly correlated',' n = ', num2str(length(redHigh{iCon}))])
ylabel('dF/F') 
xlabel('s') 




x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['running, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) 'loc_R_timecourses.pdf']),'-dpdf');
clear txt1 highRed lowRed
end 


%% scatterplot of max df/f for day 1 vs day 2 for high and low R red cells, running

for iCon = 1:nCon
   redHigh =  intersect(haveRunning_red{iCon},highRInds); %find intersection of red cells with high R 
   % and red cells that I have running data for on both days, for this contrast
    redLow =  intersect(haveRunning_red{iCon},lowRInds); %find intersection of red cells with low R 
   % and red cells that I have running data for on both days, for this contrast


figure; movegui('center') 
subplot(1,2,2)
scatter((pref_responses_loc_concat{pre}(haveRunning_red{iCon},iCon)),(pref_responses_loc_concat{post}(haveRunning_red{iCon},iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
scatter((pref_responses_loc_concat{pre}(redHigh,iCon)),(pref_responses_loc_concat{post}(redHigh,iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on

mean_pre_stat = nanmean(pref_responses_loc_concat{pre}(redHigh,iCon));
mean_post_stat = nanmean(pref_responses_loc_concat{post}(redHigh,iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_responses_loc_concat{pre}(redHigh,iCon)) / sqrt( length(redHigh));
stderror_post= std(pref_responses_loc_concat{post}(redHigh,iCon)) / sqrt( length(redHigh));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);


xlabel('pre-DART  dF/F')
ylim([-.1 0.4])
xlim([-.1 0.4])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('+HTP high r, running')
axis square
set(gca, 'TickDir', 'out')



subplot(1,2,1)
scatter((pref_responses_loc_concat{pre}(redLow,iCon)),(pref_responses_loc_concat{post}(redLow,iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
mean_pre_stat = nanmean(pref_responses_loc_concat{pre}(redLow,iCon));
mean_post_stat = nanmean(pref_responses_loc_concat{post}(redLow,iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_responses_loc_concat{pre}(redLow,iCon)) / sqrt( length(redLow));
stderror_post= std(pref_responses_loc_concat{post}(redLow,iCon)) / sqrt( length(redLow));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);% ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylabel('post-DART dF/F')
ylim([-.1 0.4])
xlim([-.1 0.4])
set(gca, 'TickDir', 'out')
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('+HTP low r, running')
axis square


sgtitle(num2str(cons(iCon)))
print(fullfile(fnout,[num2str(cons(iCon)) 'loc_maxResp_R.pdf']),'-dpdf','-bestfit')
clear txt1 highRed lowRed
end

[h1, p1]= ttest(pref_responses_loc_concat{pre}(redLow,3),pref_responses_loc_concat{post}(redLow,3));
[h2,p2]= ttest(pref_responses_loc_concat{pre}(redHigh,3),pref_responses_loc_concat{post}(redHigh,3));

%correct for four tests
p1*2
p2*2
clear h1 p1 h2 p2
%% timecourses for locomotion - stationary on a cell-by-cell average basis
tc_green_avrg_subtracted = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_subtracted = cell(1,nd); %same for red
tc_green_se_subtracted = cell(1,nd); %this will be the se across all green cells
tc_red_se_subtracted = cell(1,nd); %same for red

tc_subtracted = cell(1,nd);

for id = 1:nd
    for iCon=1:nCon
    tc_subtracted{id}(:,:,iCon)  = tc_trial_avrg_loc_concat{id}(:,:,iCon) - tc_trial_avrg_stat_concat{id}(:,:,iCon);
    tc_green_avrg_subtracted{id}(:,iCon)=nanmean(tc_subtracted{id}(:,haveRunning_green{iCon},iCon),2);
    green_std=nanstd(tc_subtracted{id}(:,haveRunning_green{iCon},iCon),[],2);
    tc_green_se_subtracted{id}(:,iCon)=green_std/sqrt(length(haveRunning_green{iCon}));
    
    tc_red_avrg_subtracted{id}(:,iCon)=nanmean(tc_subtracted{id}(:,haveRunning_red{iCon},iCon),2);
    red_std=nanstd(tc_subtracted{id}(:,haveRunning_red{iCon},iCon),[],2);
    tc_red_se_subtracted{id}(:,iCon)=red_std/sqrt(length(haveRunning_red{iCon}));
    
    clear green_std red_std
    end
end
 
%plotting
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_subtracted{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(1,2,1) %for the first day



ylim([-.05 .3]);
hold on
shadedErrorBar(t,tc_green_avrg_subtracted{pre}(:,iCon),tc_green_se_subtracted{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_subtracted{post}(:,iCon),tc_green_se_subtracted{post}(:,iCon),'b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['-HTP',' n = ', num2str(length(haveRunning_green{iCon}))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_subtracted{pre}(:,iCon),tc_red_se_subtracted{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_subtracted{post}(:,iCon),tc_red_se_subtracted{post}(:,iCon),'b');
ylim([-.05 .3]);
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['+HTP',' n = ', num2str(length(haveRunning_red{iCon}))])

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['loc - stat, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) '_subtracted_cellType_timecourses.pdf']),'-dpdf');
end 

%% subtracted scatters

pref_resp_subtracted = cell(1,nd);

for id = 1:nd
    for iCon=1:nCon
        pref_resp_subtracted{id}(:,iCon)  = pref_responses_loc_concat{id}(:,iCon) - pref_responses_stat_concat{id}(:,iCon);
    end
end


for iCon = 1:nCon
figure; movegui('center') 
subplot(1,2,2)
scatter((pref_resp_subtracted{pre}(haveRunning_green{iCon},iCon)),(pref_resp_subtracted{post}(haveRunning_green{iCon},iCon)),10,'MarkerFaceColor',[.2 .2 .2],'MarkerEdgeColor',[.2 .2 .2],'jitter', 'on', 'jitterAmount',.01)
hold on
%hline=lsline;
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
mean_pre_stat = nanmean(pref_resp_subtracted{pre}(haveRunning_green{iCon},iCon));
mean_post_stat = nanmean(pref_resp_subtracted{post}(haveRunning_green{iCon},iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_resp_subtracted{pre}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
stderror_post= std(pref_resp_subtracted{post}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
limMin=min(min(pref_resp_subtracted{pre}(haveRunning_green{iCon},iCon)),min(pref_resp_subtracted{post}(haveRunning_green{iCon},iCon)));
limMax=max(max(pref_resp_subtracted{pre}(haveRunning_green{iCon},iCon)),max(pref_resp_subtracted{post}(haveRunning_green{iCon},iCon)));
% ylim([limMin limMax])
% xlim([limMin limMax])
ylim([0 .5])
xlim([0 .5])


title('-HTP loc - stat')
axis square
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
hold off


subplot(1,2,1)
scatter((pref_resp_subtracted{pre}(haveRunning_red{iCon},iCon)),(pref_resp_subtracted{post}(haveRunning_red{iCon},iCon)),10,'MarkerFaceColor',[.2 .2 .2],'MarkerEdgeColor',[.2 .2 .2],'jitter', 'on', 'jitterAmount',.01)
hold on
%hline=lsline;
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
mean_pre_stat = nanmean(pref_resp_subtracted{pre}(haveRunning_red{iCon},iCon));
mean_post_stat = nanmean(pref_resp_subtracted{post}(haveRunning_red{iCon},iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_resp_subtracted{pre}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
stderror_post= std(pref_resp_subtracted{post}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);% ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
limMin=min(min(pref_resp_subtracted{pre}(haveRunning_red{iCon},iCon)),min(pref_resp_subtracted{post}(haveRunning_red{iCon},iCon)));
limMax=max(max(pref_resp_subtracted{pre}(haveRunning_red{iCon},iCon)),max(pref_resp_subtracted{post}(haveRunning_red{iCon},iCon)));
ylim([0 .5])
xlim([0 .5])
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
title('+HTP loc-stat')
axis square
hold off




sgtitle(num2str(cons(iCon)))
print(fullfile(fnout,[num2str(cons(iCon)) 'subtracted_resp_crossDay.pdf']),'-dpdf','-bestfit')

% edges = [limMin:((limMax-limMin)/20):limMax];
% figure;
% subplot(1,2,1)
% hist(pref_resp_subtracted{pre}(haveRunning_green{iCon},iCon),edges)
% xlim([limMin limMax])
% subplot(1,2,2)
% hist(pref_resp_subtracted{post}(haveRunning_green{iCon},iCon),edges)
% xlim([limMin limMax])
% sgtitle(num2str(cons(iCon)))
% print(fullfile(fnout,[num2str(cons(iCon)) 'subtracted_distribution.pdf']),'-dpdf','-bestfit')


clear mean_pre_stat mean_post_stat stderror_post stderror_pre
end

%%
iCon=1
[h1, p1]= ttest(pref_resp_subtracted{pre}(haveRunning_green{iCon},iCon),pref_resp_subtracted{post}(haveRunning_green{iCon},iCon));
[h2,p2]= ttest(pref_resp_subtracted{pre}(haveRunning_red{iCon},iCon),pref_resp_subtracted{post}(haveRunning_red{iCon},iCon));

%correct for two tests
p1*2
p2*2

clear h1 p1 h2 p2 

meanEffectSize(pref_resp_subtracted{pre}(haveRunning_green{iCon},iCon),pref_resp_subtracted{post}(haveRunning_green{iCon},iCon))
meanEffectSize(pref_resp_subtracted{pre}(haveRunning_red{iCon},iCon),pref_resp_subtracted{post}(haveRunning_red{iCon},iCon))

%% calculate norm_diff
norm_diff = nan(2,nCon,nKeep_total);
for i = 1:nKeep_total
    for iCon = 1:nCon
        %for stationary trials
        mean_pre_stat = mean(pref_allTrials_stat_concat{iCon,pre}{i});
        mean_post_stat=mean(pref_allTrials_stat_concat{iCon,post}{i});
        std_pre_stat = std(pref_allTrials_stat_concat{iCon,pre}{i});
        norm_diff_stat = (mean_post_stat-mean_pre_stat) / std_pre_stat;

        %for running trials
        mean_pre_loc = mean(pref_allTrials_loc_concat{iCon,pre}{i});
        mean_post_loc=mean(pref_allTrials_loc_concat{iCon,post}{i});
        std_pre_loc = std(pref_allTrials_loc_concat{iCon,pre}{i});
        norm_diff_loc = (mean_post_loc-mean_pre_loc)/ std_pre_loc;

        %putting data into matrix
        norm_diff(1,iCon,i)=norm_diff_stat; %first is stationary
        norm_diff(2,iCon,i)=norm_diff_loc; %second is running
clear mean_pre_stat mean_post_stat std_pre_stat mean_pre_loc mean_post_loc std_pre_loc norn_diff_stat norm_diff_loc
    end 
end
%remove any infiinty values resulting from divisions by zero, and turn
%those into NANs instead
norm_diff(find(norm_diff == -Inf))=NaN;
norm_diff(find(norm_diff == Inf))=NaN;

%make a table of suppresses and facilitated cells for the cells that have
%both stationary and running within each condition
facil=norm_diff(:,:,:)>=1;
supp=norm_diff(:,:,:)<=-1;

supp_table_stat_red=nan(nCon,1);
facil_table_stat_red=nan(nCon,1);
supp_table_loc_red=nan(nCon,1);
facil_table_loc_red=nan(nCon,1);
supp_table_stat_green=nan(nCon,1);
facil_table_stat_green=nan(nCon,1);
supp_table_loc_green=nan(nCon,1);
facil_table_loc_green=nan(nCon,1);

for iCon = 1:nCon
        supp_table_stat_green(iCon)=sum(supp(1,iCon,haveRunning_green{iCon}),3)/length(haveRunning_green{iCon});
        supp_table_loc_green(iCon)=sum(supp(2,iCon,haveRunning_green{iCon}),3)/length(haveRunning_green{iCon});

        facil_table_stat_green(iCon)=sum(facil(1,iCon,haveRunning_green{iCon}),3)./length(haveRunning_green{iCon});
        facil_table_loc_green(iCon)=sum(facil(2,iCon,haveRunning_green{iCon}),3)/length(haveRunning_green{iCon});

        supp_table_stat_red(iCon)=sum(supp(1,iCon,haveRunning_red{iCon}),3)/length(haveRunning_red{iCon});
        supp_table_loc_red(iCon)=sum(supp(2,iCon,haveRunning_red{iCon}),3)/length(haveRunning_red{iCon});

        facil_table_stat_red(iCon)=sum(facil(1,iCon,haveRunning_red{iCon}),3)./length(haveRunning_red{iCon});
        facil_table_loc_red(iCon)=sum(facil(2,iCon,haveRunning_red{iCon}),3)/length(haveRunning_red{iCon});

end


    %%
figure;
subplot(2,2,1)
b=bar([1,2,3],[supp_table_stat_red(:,1),supp_table_loc_red(:,1)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
b(1).FaceColor="#70D0F6";
b(2).FaceColor="#0C8ABB";
xticklabels({'25','50','100'})
hold on
title('Suppressed')
ylim([0 .4])
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast"])
set(gca,'TickDir','out')
box off

subplot(2,2,2)
b=bar([1,2,3],[facil_table_stat_red(:,1),facil_table_loc_red(:,1)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
b(1).FaceColor="#C983B1";
b(2).FaceColor="#883367";
xticklabels({'25','50','100'})
hold on
title('Facilitated')
ylim([0 .4])
%ylabel(["Fraction SST cells"]) 
xlabel(["Contrast"])
set(gca,'TickDir','out')
box off

subplot(2,2,3)
b=bar([1,2,3],[supp_table_stat_green(:,1),facil_table_stat_green(:,1)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
b(1).FaceColor="#70D0F6";
b(2).FaceColor="#0C8ABB";
xticklabels({'25','50','100'})
hold on
%title('Suppressed')
ylim([0 .4])
ylabel(["Fraction Pyr cells"]) 
xlabel(["Contrast"])
set(gca,'TickDir','out')
box off

subplot(2,2,4)
b=bar([1,2,3],[facil_table_stat_green(:,1),facil_table_loc_green(:,1)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
b(1).FaceColor="#C983B1";
b(2).FaceColor="#883367";
xticklabels({'25','50','100'})
hold on
%title('Facilitated')
ylim([0 .4])
xlabel(["Contrast"])
set(gca,'TickDir','out')
box off


print(fullfile(fnout,'mached_facil_supp_byState.pdf'),'-dpdf');
%% alternative layout of same plot as above
figure;
subplot(2,2,1)
b=bar([1,2,3],[supp_table_stat_red(:,1),facil_table_stat_red(:,1)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
b(1).FaceColor="#70D0F6";
b(2).FaceColor="#C983B1";
xticklabels({'25','50','100'})
hold on
title('Stationary')
ylim([0 .4])
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast"])
set(gca,'TickDir','out')
box off

subplot(2,2,2)
b=bar([1,2,3],[supp_table_loc_red(:,1),facil_table_loc_red(:,1)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
b(1).FaceColor="#0C8ABB";
b(2).FaceColor="#883367";
xticklabels({'25','50','100'})
hold on
title('Running')
ylim([0 .4])
%ylabel(["Fraction SST cells"]) 
xlabel(["Contrast"])
set(gca,'TickDir','out')
box off

subplot(2,2,3)
b=bar([1,2,3],[supp_table_stat_green(:,1),facil_table_stat_green(:,1)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
b(1).FaceColor="#70D0F6";
b(2).FaceColor="#C983B1";
xticklabels({'25','50','100'})
hold on
%title('Suppressed')
ylim([0 .4])
ylabel(["Fraction Pyr cells"]) 
xlabel(["Contrast"])
set(gca,'TickDir','out')
box off

subplot(2,2,4)
b=bar([1,2,3],[supp_table_loc_green(:,1),facil_table_loc_green(:,1)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
b(1).FaceColor="#0C8ABB";
b(2).FaceColor="#883367";
xticklabels({'25','50','100'})
hold on
%title('Facilitated')
%ylim([0 .4])
xlabel(["Contrast"])
set(gca,'TickDir','out')
box off
%% compute chi squares 
%compute chi squares for stat red
for iCon = 1:nCon
    N=length(haveRunning_red{iCon});
    n1=round(supp_table_stat_red(iCon,1)*N);
    n2=round(facil_table_stat_red(2)*N);
    x1 = [repmat('a',N,1); repmat('b',N,1)];
    x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
    [tbl,chi2stat1,p1] = crosstab(x1,x2);
    p1*3
end

%compute chi squares for running red
for iCon = 1:nCon
    N=length(haveRunning_red{iCon});
    n1=round(supp_table_loc_red(iCon,1)*N);
    n2=round(facil_table_loc_red(2)*N);
    x1 = [repmat('a',N,1); repmat('b',N,1)];
    x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
    [tbl,chi2stat1,p1] = crosstab(x1,x2);
    p1*3
end


%compute chi squares for stat green
for iCon = 1:nCon
    N=length(haveRunning_green{iCon});
    n1=round(supp_table_stat_green(iCon,1)*N);
    n2=round(facil_table_stat_green(2)*N);
    x1 = [repmat('a',N,1); repmat('b',N,1)];
    x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
    [tbl,chi2stat1,p1] = crosstab(x1,x2);
    p1*3
end

%compute chi squares for running green
for iCon = 1:nCon
    N=length(haveRunning_green{iCon});
    n1=round(supp_table_loc_green(iCon,1)*N);
    n2=round(facil_table_loc_green(2)*N);
    x1 = [repmat('a',N,1); repmat('b',N,1)];
    x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
    [tbl,chi2stat1,p1] = crosstab(x1,x2);
    p1*3
end

%% plot contrast response

conResp_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_stat = cell(1,nd); %same for red
conResp_green_se_stat = cell(1,nd); %this will be the se across all green cells
conResp_red_se_stat = cell(1,nd); %same for red


conResp_green_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_loc = cell(1,nd); %same for red
conResp_green_se_loc = cell(1,nd); %this will be the se across all green cells
conResp_red_se_loc = cell(1,nd); %same for red


for id = 1:nd
    conResp_green_avrg_stat{id} = nan(1,nCon);
    conResp_red_avrg_stat{id} = nan(1,nCon);
    conResp_green_avrg_loc{id} = nan(1,nCon);
    conResp_red_avrg_loc{id} = nan(1,nCon);
        
    conResp_green_avrg_loc{id}=nanmean(pref_responses_loc_concat{id}(haveRunning_green{iCon} ,:),1);
    green_std=nanstd(pref_responses_loc_concat{id}(haveRunning_green{iCon},:),1);
    conResp_green_se_loc{id}=green_std/sqrt(length(haveRunning_green{iCon}));
    
    conResp_red_avrg_loc{id}=nanmean(pref_responses_loc_concat{id}(haveRunning_red{iCon},:),1);
    red_std=nanstd(pref_responses_loc_concat{id}(haveRunning_red{iCon},:),1);
    conResp_red_se_loc{id}=red_std/sqrt(length(haveRunning_red{iCon}));
    
    clear green_std red_std
 
        
    conResp_green_avrg_stat{id}=mean(pref_responses_stat_concat{id}(haveRunning_green{iCon},:),1,'omitnan');
    green_std=std(pref_responses_stat_concat{id}(haveRunning_green{iCon},:),1,'omitnan');
    conResp_green_se_stat{id}=green_std/sqrt(length(haveRunning_green{iCon}));
    
    conResp_red_avrg_stat{id}=mean(pref_responses_stat_concat{id}(haveRunning_red{iCon},:),1,'omitnan');
    red_std=std(pref_responses_stat_concat{id}(haveRunning_red{iCon},:),1,'omitnan');
    conResp_red_se_stat{id}=red_std/sqrt(length(haveRunning_red{iCon}));
    
    clear green_std red_std
 
end
%%

figure
subplot(1,2,2) %for the first day
errorbar(cons,conResp_green_avrg_stat{pre},conResp_green_se_stat{pre},'Ok');
hold on
errorbar(cons,conResp_green_avrg_stat{post},conResp_green_se_stat{post},'Ob');
title('-HTP')
xlabel('contrast') 
set(gca, 'TickDir', 'out')
xlim([0 1.1])
ylim([0 .1])
xticks([.25 .5 1])
box off

subplot(1,2,1) %for the second day
errorbar(cons,conResp_red_avrg_stat{pre},conResp_red_se_stat{pre},'Ok');
hold on
errorbar(cons,conResp_red_avrg_stat{post},conResp_red_se_stat{post},'Ob');
title('+HTP')
xlabel('contrast') 
xlim([0 1.1])
ylim([0 .1])
xticks([.25 .5 1])
set(gca, 'TickDir', 'out')
ylabel('dF/F, pref ori') 
box off

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])



sgtitle(['stationary' ])

print(fullfile(fnout,['contrast_response.pdf']),'-dpdf');



figure
subplot(1,2,2) %for the first day
errorbar(cons,conResp_green_avrg_loc{pre},conResp_green_se_loc{pre},'Ok');
hold on
errorbar(cons,conResp_green_avrg_loc{post},conResp_green_se_loc{post},'Ob');
title('-HTP')
xlabel('contrast') 
xlim([0 1.1])
ylim([0 .2])
xticks([.25 .5 1])
set(gca, 'TickDir', 'out')
box off

subplot(1,2,1) 
errorbar(cons,conResp_red_avrg_loc{pre},conResp_red_se_loc{pre},'Ok');
hold on
errorbar(cons,conResp_red_avrg_loc{post},conResp_red_se_loc{post},'Ob');
title('+HTP')
 ylim([0 .2])
xlim([0 1.1])
xticks([.25 .5 1])
xlabel('contrast') 
ylabel('dF/F, pref ori') 
set(gca, 'TickDir', 'out')
box off

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])



sgtitle(['Running' ])

print(fullfile(fnout,['loc_contrast_resposnse.pdf']),'-dpdf');

%% mouse by mouse means scatterplot of max df/f for day 1 vs day 2, and each subplot is one cell type
% this is for all keep cells, not only those matched across running
green_means_stat=cell(1,nd);
red_means_stat=cell(1,nd);
green_means_loc=cell(1,nd);
red_means_loc=cell(1,nd);

for id=1:nd
    temp_means1 = nan(nSess,nCon);
    temp_means2 = nan(nSess,nCon);
    temp_means3 = nan(nSess,nCon);
    temp_means4 = nan(nSess,nCon);
    for iMouse = 1:nSess
        for iCon = 1:nCon
            inds1 = intersect(haveRunning_green{iCon}, mouseInds{iMouse});
            temp_means1(iMouse, iCon)=mean(pref_responses_stat_concat{id}(inds1,iCon),'omitnan');

            inds2 = intersect(haveRunning_red{iCon}, mouseInds{iMouse});
            temp_means2(iMouse, iCon)=mean(pref_responses_stat_concat{id}(inds2,iCon),'omitnan');

            temp_means3(iMouse, iCon)=mean(pref_responses_loc_concat{id}(inds1,iCon),'omitnan');

            temp_means4(iMouse, iCon)=mean(pref_responses_loc_concat{id}(inds2,iCon),'omitnan');
        end
    end
    green_means_stat{id}=temp_means1;
    red_means_stat{id}=temp_means2;
    green_means_loc{id}=temp_means3;
    red_means_loc{id}=temp_means4;
end
z=12;
y=z-1;

cmap = jet(length(mice)); % Make nMice colors.


for iCon = 1:nCon
figure; movegui('center') 
subplot(2,2,2)
scatter((pref_responses_stat_concat{pre}(haveRunning_green{iCon},iCon)),(pref_responses_stat_concat{post}(haveRunning_green{iCon},iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter((green_means_stat{pre}(:,iCon)),(green_means_stat{post}(:,iCon)),10,cmap,'filled')
mean_pre_stat = nanmean(pref_responses_stat_concat{pre}(haveRunning_green{iCon},iCon));
mean_post_stat = nanmean(pref_responses_stat_concat{post}(haveRunning_green{iCon},iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_responses_stat_concat{pre}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
stderror_post= std(pref_responses_stat_concat{post}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
limMin=min(min(pref_responses_stat_concat{pre}(haveRunning_green{iCon},iCon)),min(pref_responses_stat_concat{post}(haveRunning_green{iCon},iCon)));
limMax=max(max(pref_responses_stat_concat{pre}(haveRunning_green{iCon},iCon)),max(pref_responses_stat_concat{post}(haveRunning_green{iCon},iCon)));
ylim([limMin limMax])
xlim([limMin limMax])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('-HTP stationary')
axis square
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
hold off


subplot(2,2,1)
scatter((pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon)),(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter(red_means_stat{pre}(:,iCon),red_means_stat{post}(:,iCon),10,cmap, 'filled')
mean_pre_stat = nanmean(pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon));
mean_post_stat = nanmean(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
stderror_post= std(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);% ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
limMin=min(min(pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon)),min(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon)));
limMax=max(max(pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon)),max(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon)));
ylim([limMin limMax])
xlim([limMin limMax])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
title('+HTP stationary')
axis square
hold off

subplot(2,2,4)
scatter((pref_responses_loc_concat{pre}(haveRunning_green{iCon},iCon)),(pref_responses_loc_concat{post}(haveRunning_green{iCon},iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter((green_means_loc{pre}(:,iCon)),(green_means_loc{post}(:,iCon)),10,cmap, 'filled')
mean_pre_stat = nanmean(pref_responses_loc_concat{pre}(haveRunning_green{iCon},iCon));
mean_post_stat = nanmean(pref_responses_loc_concat{post}(haveRunning_green{iCon},iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_responses_loc_concat{pre}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
stderror_post= std(pref_responses_loc_concat{post}(haveRunning_green{iCon},iCon)) / sqrt( length(haveRunning_green{iCon}));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
limMin=min(min(pref_responses_loc_concat{pre}(haveRunning_green{iCon},iCon)),min(pref_responses_loc_concat{post}(haveRunning_green{iCon},iCon)));
limMax=max(max(pref_responses_loc_concat{pre}(haveRunning_green{iCon},iCon)),max(pref_responses_loc_concat{post}(haveRunning_green{iCon},iCon)));
ylim([limMin limMax])
xlim([limMin limMax])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('-HTP running')
axis square
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
hold off

subplot(2,2,3)
scatter((pref_responses_loc_concat{pre}(haveRunning_red{iCon},iCon)),(pref_responses_loc_concat{post}(haveRunning_red{iCon},iCon)),10,'MarkerEdgeColor',[.5 .5 .5],'jitter', 'on', 'jitterAmount',.01)
hold on
scatter((red_means_loc{pre}(:,iCon)),(red_means_loc{post}(:,iCon)),10, cmap,'filled')
mean_pre_stat = nanmean(pref_responses_loc_concat{pre}(haveRunning_red{iCon},iCon));
mean_post_stat = nanmean(pref_responses_loc_concat{post}(haveRunning_red{iCon},iCon));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= std(pref_responses_stat_concat{pre}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
stderror_post= std(pref_responses_stat_concat{post}(haveRunning_red{iCon},iCon)) / sqrt( length(haveRunning_red{iCon}));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);
xlabel('pre-DART  dF/F')
limMin=min(min(pref_responses_loc_concat{pre}(haveRunning_red{iCon},iCon)),min(pref_responses_loc_concat{post}(haveRunning_red{iCon},iCon)));
limMax=max(max(pref_responses_loc_concat{pre}(haveRunning_red{iCon},iCon)),max(pref_responses_loc_concat{post}(haveRunning_red{iCon},iCon)));
ylim([limMin limMax])
xlim([limMin limMax])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('+HTP running')
uistack(hline,'bottom');
axis square
hold off
set(gca, 'TickDir', 'out')



sgtitle(num2str(cons(iCon)))
print(fullfile(fnout,[num2str(cons(iCon)) 'maxResp_byMouse.pdf']),'-dpdf','-bestfit')
clear mean_pre_stat mean_post_stat stderror_post stderror_pre
end

%% direction tuning
%% linear direction plot 

green_dir_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
red_dir_avrg_stat = cell(1,nd); %same for red
green_dir_se_stat = cell(1,nd); %this will be the se across all green cells
red_dir_se_stat = cell(1,nd); %same for red

for iCon = 1:nCon
    for id = 1:nd
       
        green_dir_avrg_stat{id}=nanmean(norm_dir_resp_stat_concat{id}(haveRunning_green{iCon},:,iCon),1);
        green_std=nanstd(norm_dir_resp_stat_concat{id}(haveRunning_green{iCon},:,iCon),[],1);
        green_dir_se_stat{id}=green_std/sqrt(length(haveRunning_green{iCon}));
        green_dir_avrg_stat{id}=circshift(green_dir_avrg_stat{id},4);
        green_dir_se_stat{id}=circshift(green_dir_se_stat{id},4);
        
        red_dir_avrg_stat{id}=nanmean(norm_dir_resp_stat_concat{id}(haveRunning_red{iCon},:,iCon),1);
        red_std=nanstd(norm_dir_resp_stat_concat{id}(haveRunning_red{iCon},:,iCon),[],1);
        red_dir_se_stat{id}=red_std/sqrt(length(haveRunning_red{iCon}));
        red_dir_avrg_stat{id}=circshift(red_dir_avrg_stat{id},4);
        red_dir_se_stat{id}=circshift(red_dir_se_stat{id},4);
        clear green_std red_std
        
    end
    
    
    
    green_dir_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
    red_dir_avrg_loc = cell(1,nd); %same for red
    green_dir_se_loc = cell(1,nd); %this will be the se across all green cells
    red_dir_se_loc = cell(1,nd); %same for red
    
    for id = 1:nd
       
        green_dir_avrg_loc{id}=nanmean(norm_dir_resp_loc_concat{id}(haveRunning_green{iCon},:,iCon),1);
        green_std=nanstd(norm_dir_resp_loc_concat{id}(haveRunning_green{iCon},:,iCon),[],1);
        green_dir_se_loc{id}=green_std/sqrt(length(haveRunning_green{iCon}));
        green_dir_avrg_loc{id}=circshift(green_dir_avrg_loc{id},4);
        green_dir_se_loc{id}=circshift(green_dir_se_loc{id},4);
        
        red_dir_avrg_loc{id}=nanmean(norm_dir_resp_loc_concat{id}(haveRunning_red{iCon},:,iCon),1);
        red_std=nanstd(norm_dir_resp_loc_concat{id}(haveRunning_red{iCon},:,iCon),[],1);
        red_dir_se_loc{id}=red_std/sqrt(length(haveRunning_red{iCon}));
        red_dir_avrg_loc{id}=circshift(red_dir_avrg_loc{id},4);
        red_dir_se_loc{id}=circshift(red_dir_se_loc{id},4);
        clear green_std red_std
        
    end
    
    
    
    figure
    subplot(2,2,1)
    errorbar(dirs,green_dir_avrg_stat{pre},green_dir_se_stat{pre},'k')
    hold on
    errorbar(dirs,green_dir_avrg_stat{post},green_dir_se_stat{post},'b')
    title(['Stationary, ', num2str(length(haveRunning_green{iCon})),' Pyr'])
    set(gca, 'TickDir', 'out')
    axis square
    box off
    ylabel('dF/F')
    % ylim([-0.05 .1])
    
    subplot(2,2,2)
    errorbar(dirs,red_dir_avrg_stat{pre},red_dir_se_stat{pre},'k')
    hold on
    errorbar(dirs,red_dir_avrg_stat{post},red_dir_se_stat{post},'b')
    title(['Stationary, ', num2str(length(haveRunning_red{iCon})),' SST'])
    set(gca, 'TickDir', 'out')
    axis square
    box off
    % ylim([-0.05 .1])
    
    subplot(2,2,3)
    errorbar(dirs,green_dir_avrg_loc{pre},green_dir_se_loc{pre},'k')
    hold on
    errorbar(dirs,green_dir_avrg_loc{post},green_dir_se_loc{post},'b')
    title('Running, Pyr')
    set(gca, 'TickDir', 'out')
    axis square
    box off
    xlabel('normalized direction')
    ylabel('dF/F')
    % ylim([-0.05 .2])
    
    
    subplot(2,2,4)
    errorbar(dirs,red_dir_avrg_loc{pre},red_dir_se_loc{pre},'k')
    hold on
    errorbar(dirs,red_dir_avrg_loc{post},red_dir_se_loc{post},'b')
    title('Running, SST')
    set(gca, 'TickDir', 'out')
    axis square
    box off
    xlabel('normalized direction')
    % ylim([-0.05 .2])
    
    sgtitle(['Normalized direction tuning ',num2str(cons(iCon))])
    
    
    print(fullfile(fnout,[num2str(cons(iCon)),'dirTuning.pdf']),'-dpdf','-bestfit')
end

%% compared effect during stationary to effect during running
diffStat = pref_responses_stat_concat{post}-pref_responses_stat_concat{pre};
diffLoc = pref_responses_loc_concat{post}-pref_responses_loc_concat{pre};

for iCon = 1:nCon
    figure
    scatter(diffStat(haveRunning_red{iCon},iCon),diffLoc(haveRunning_red{iCon},iCon));
    ylabel('diff run')
    xlabel('diff stat')
    xlim([-.8 .8])
    ylim([-.8 .8])
end 