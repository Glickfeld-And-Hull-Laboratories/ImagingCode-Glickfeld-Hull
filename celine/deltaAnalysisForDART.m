 
clear all; clear global; close all
clc
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc =  behavConstsDART; %directories
eval(ds);


sess_list = [138 142 163 171 178 190, 294,307,333,323,303,311,319,329,355,359];%enter all the sessions you want to concatenate
nSess=length(sess_list);

nd=2;%hard coding for two days per experimental session

% INDICATE THE PRE VS. POST DAYS DEPENDING ON THE ORDER OF MATCHING
pre=2;
post=1;

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
        
        fnout = fullfile(rc.achAnalysis,'SST_YM90K',expt(sess_list(1)).mouse,['multiday_' dart_str],d);
else
    fnout= fullfile(rc.achAnalysis,strcat('concat', sess_title),d);
end
mkdir(fnout);
cd(fnout)
clear d sess_title



%% concatenating data
mice=[];
tc_trial_avrg_stat_concat=cell(1,nd);
tc_trial_avrg_loc_concat=cell(1,nd);
tc_trial_avrg_keep_allCon_stat_concat=cell(1,nd);
tc_trial_avrg_keep_allCon_loc_concat=cell(1,nd);
resp_keep_concat=cell(1,nd);
resp_max_keep_concat=cell(1,nd);
pref_responses_loc_concat=cell(1,nd);
pref_responses_stat_concat=cell(1,nd);
pref_responses_allCon__stat_concat=cell(2,nd);
pref_responses_allCon__loc_concat=cell(2,nd);
RIx_concat=cell(1,nd);
dirs_concat=[];
cons_concat=[];
green_concat=[];
red_concat=[];
redRsq_concat=[];
dfof_max_diff_concat=[];
nKeep_concat=[];
%LMI_concat = cell(1,nd);
data_resp_concat = cell(1,nd);
red_fluor_concat=[];
green_fluor_concat=[];


for iSess = 1:nSess
    day_id = sess_list(iSess)
    mouse = expt(day_id).mouse;
    mice=[mice;mouse];

    if expt(day_id).multiday_timesincedrug_hours>0
        dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
    else
        dart_str = 'control';
    end
    fn_multi = fullfile(rc.achAnalysis,'SST_YM90K',mouse,['multiday_' dart_str]);

    load(fullfile(fn_multi,'tc_keep.mat'))
    load(fullfile(fn_multi,'resp_keep.mat'))
    load(fullfile(fn_multi,'input.mat'))
    % load(fullfile(fn_multi,'locomotion.mat'))
    % load(fullfile(fn_multi,'fluor_intensity.mat'))
    % load(fullfile(fn_multi,'HT_pyr_relationship.mat'))

    nKeep = size(tc_trial_avrg_stat{post},2);


    %tells the contrast, direction and orientation for each trial each day
    tCon_match = cell(1,nd);
    tDir_match = cell(1,nd);
    tOri_match = cell(1,nd);

    %find the contrasts, directions and orientations for each day
    for id = 1:nd
        tCon_match{id} = celleqel2mat_padded(input(id).tGratingContrast(1:nTrials(id)));
        tDir_match{id} = celleqel2mat_padded(input(id).tGratingDirectionDeg(1:nTrials(id)));
        tOri_match{id} = tDir_match{id};
        tOri_match{id}(find(tDir_match{id}>=180)) = tDir_match{id}(find(tDir_match{id}>=180))-180;
    end
    dirs = unique(tDir_match{post});
    cons = unique(tCon_match{post});

    nOn = input(1).nScansOn;
    nOff = input(1).nScansOff;
    %start conatenating
    dirs_concat = [dirs_concat;dirs]; %I will organize thisas rows so I can subsequently make sure everything matches
    cons_concat = [cons_concat;cons];
    red_concat = [red_concat, red_keep_logical];
    green_concat = [green_concat, green_keep_logical];
    nKeep_concat = [nKeep_concat,nKeep];
%    redRsq_concat = [redRsq_concat,Rsq_red];
    clear cons
for id = 1:nd
        tc_trial_avrg_keep_allCon_stat_concat{id} =cat(2,tc_trial_avrg_keep_allCon_stat_concat{id},tc_trial_avrg_keep_allCon_stat{id}(:,:));
        tc_trial_avrg_keep_allCon_loc_concat{id} =cat(2,tc_trial_avrg_keep_allCon_loc_concat{id},tc_trial_avrg_keep_allCon_loc{id}(:,:));
        tc_trial_avrg_stat_concat{id} =cat(2,tc_trial_avrg_stat_concat{id},tc_trial_avrg_stat{id}(:,:,:));
        tc_trial_avrg_loc_concat{id} =cat(2,tc_trial_avrg_loc_concat{id},tc_trial_avrg_loc{id}(:,:,:));
        resp_keep_concat{id}=cat(1,resp_keep_concat{id},resp_keep{id});
        resp_max_keep_concat{id}=cat(1,resp_max_keep_concat{id},resp_max_keep{id}(:,:));
%        LMI_concat{id}=cat(1,LMI_concat{id},LMI{id}(:,:));
        pref_responses_loc_concat{id}=cat(1,pref_responses_loc_concat{id},pref_responses_loc{id}(:,:));
        pref_responses_stat_concat{id}=cat(1,pref_responses_stat_concat{id},pref_responses_stat{id}(:,:));
        pref_responses_allCon__stat_concat{1,id}=cat(2,pref_responses_allCon__stat_concat{1,id},pref_responses_allCon_stat{1,id}(:,:));
        pref_responses_allCon__stat_concat{2,id}=cat(2,pref_responses_allCon__stat_concat{2,id},pref_responses_allCon_stat{2,id}(:,:));
        pref_responses_allCon__loc_concat{1,id}=cat(2,pref_responses_allCon__loc_concat{1,id},pref_responses_allCon_loc{1,id}(:,:));
        pref_responses_allCon__loc_concat{2,id}=cat(2,pref_responses_allCon__loc_concat{2,id},pref_responses_allCon_loc{2,id}(:,:));
%        RIx_concat{id}=cat(1,RIx_concat{id},sum(RIx{id}));
    end
    dfof_max_diff_concat=cat(1,dfof_max_diff_concat,dfof_max_diff);
%    green_fluor_concat=cat(2,green_fluor_concat,green_fluor_keep);
 %   red_fluor_concat=cat(2,red_fluor_concat,red_fluor_keep);
end

clear mouse day_id nKeep iSess fn_multi cons dirs
clear explanation1 resp_keep tc_trial_avrg_keep_allCon pref_responses_allCon sig_diff pref_con_keep pref_dir_keep tDir_match tOri_match tCon_match data_trial_keep nTrials tc_trial_avrg_keep green_keep_logical red_keep_logical green_ind_keep red_ind_keep
clear LMI RIx locCounts locResp locTCs statResp statTCs wheel_tc
clear data_con_resp_keep data_dir_resp_keep data_rep_keep dfof_max_diff dfof_max_diff_raw explanation2 resp_max_keep data_resp_keep
clear red_fluor_all red_fluor_match green_fluor_match green_fluor_match red_fluor_keep green_fluor_keep
clear tc_trial_avrg_loc tc_trial_avrg_stat tc_trial_avrg_keep_allCond pref_responses_allCond
clear linCellProps responseByCond responseByCondProps sig_corr_red Rsq_red

red_ind_concat = find(red_concat);
green_ind_concat = find(green_concat);
%% check that contrasts and dirs are matched
cons=unique(cons_concat);
nCon = length(cons)
dirs=unique(dirs_concat);
nDir=length(dirs);
nKeep_total = sum(nKeep_concat);

mouseNames=[];
for iMouse = 1:nSess
    mouseNames=[mouseNames, string(mice(iMouse,:))]
end
clear  iMouse

% make dfof summary table for statistics
mouseID=[];
for imouse =1:nSess
   ID_string=mouseNames(imouse);
   thisID = repelem(ID_string,nKeep_concat(imouse));
   mouseID=[mouseID,thisID];
end
mouseID=mouseID';
%% find cells that I have running data for on both days
haveRunning_pre = ~isnan(pref_responses_loc_concat{pre});
haveRunning_post = ~isnan(pref_responses_loc_concat{post});
haveRunning_both = find(haveRunning_pre.* haveRunning_post);
haveRunning_green = intersect(haveRunning_both, green_ind_concat);
haveRunning_red = intersect(haveRunning_both, red_ind_concat);

%% see how many cells are in the DART sessions (first 10 sessions)
nKeep_total=sum(nKeep_concat);
nDART=sum(nKeep_concat(1:10));
%% get delta values
for id = 1:nd
    pref_responses_stat_rect{id} = pref_responses_stat_concat{id};
    pref_responses_loc_rect{id} = pref_responses_loc_concat{id};
    pref_responses_stat_rect{id}(pref_responses_stat_rect{id} < 0)=0;
    pref_responses_loc_rect{id}(pref_responses_loc_rect{id} < 0)=0;
end

delta_stat=pref_responses_stat_concat{post}-pref_responses_stat_concat{pre};
norm_delta_stat=(pref_responses_stat_rect{post}-pref_responses_stat_rect{pre})./(pref_responses_stat_rect{post}+pref_responses_stat_rect{pre});
delta_loc=pref_responses_loc_concat{post}-pref_responses_loc_concat{pre};
norm_delta_loc=(pref_responses_loc_rect{post}-pref_responses_loc_rect{pre})./(pref_responses_loc_rect{post}+pref_responses_loc_rect{pre});

drug = ones(nKeep_total,1);
drug(nDART+1:nKeep_total) = 2; %change the drug label for PEG sessions from 1 to 2
%% plot contrast response for stationary cells, not matched for locomotion

conResp_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_stat = cell(1,nd); %same for red
conResp_green_se_stat = cell(1,nd); %this will be the se across all green cells
conResp_red_se_stat = cell(1,nd); %same for red

for iDrug = 1:2
    
    thisDrug = find(drug == iDrug);
    green_thisDrug = intersect(thisDrug,green_ind_concat);
    red_thisDrug = intersect(thisDrug,red_ind_concat);
        
    conResp_green_avrg_stat{iDrug}=mean(norm_delta_stat(green_thisDrug,:),1,'omitmissing');
    green_std=std(norm_delta_stat(green_thisDrug,:),1,'omitmissing');
    conResp_green_se_stat{iDrug}=green_std/sqrt(length(green_thisDrug));
    
    conResp_red_avrg_stat{iDrug}=mean(norm_delta_stat(red_thisDrug,:),1,'omitmissing');
    red_std=std(norm_delta_stat(red_thisDrug,:),1,'omitmissing');
    conResp_red_se_stat{iDrug}=red_std/sqrt(length(red_thisDrug));
    
    clear green_std red_std
 
end


PEG = 2;
DART=1;

figure
subplot(1,2,1) %for the second day
errorbar(cons,conResp_red_avrg_stat{PEG},conResp_red_se_stat{PEG},'color',"#EA9010");
hold on
errorbar(cons,conResp_red_avrg_stat{DART},conResp_red_se_stat{DART},'b');
title('SST')
% ylabel('Post - pre/post+pre')
% xlabel('contrast') 
xlim([0 1.1])
ylim([-.4 .2])
xticks([.25 .5 1])
set(gca, 'TickDir', 'out')
box off

subplot(1,2,2) %for the first day
errorbar(cons,conResp_green_avrg_stat{PEG},conResp_green_se_stat{PEG},'--','color',"#EA9010");
hold on
errorbar(cons,conResp_green_avrg_stat{DART},conResp_green_se_stat{DART},'--b');
title('Pyr')
% xlabel('contrast') 
% set(gca, 'TickDir', 'out')
xlim([0 1.1])
ylim([-.4 .2])
xticks([.25 .5 1])
box off

x0=5;
y0=5;
width=2.2;
height=1.6;
set(gcf,'units','inches','position',[x0,y0,width,height])

%sgtitle(['Stationary, not fully matched' ])

% green_PEG = intersect(find(drug==0),green_ind_concat);
% green_DART = intersect(find(drug==1),green_ind_concat);
% red_PEG = intersect(find(drug==0),red_ind_concat);
% red_DART = intersect(find(drug==1),red_ind_concat);
% 
% %two-sample t-tests for red, stationary
% [h1, p1]= ttest2(norm_delta_stat(red_PEG,1),norm_delta_stat(red_DART,1));
% [h2, p2]= ttest2(norm_delta_stat(red_PEG,2),norm_delta_stat(red_DART,2));
% [h3, p3]= ttest2(norm_delta_stat(red_PEG,3),norm_delta_stat(red_DART,3));
% %correct for 3 tests
% p1*3
% p2*3
% p3*3
% clear h1 p1 h2 p2 h3 p3
% 
% %two-sample t-tests for green, stationary
% [h1, p1]= ttest2(norm_delta_stat(green_PEG,1),norm_delta_stat(green_DART,1));
% [h2, p2]= ttest2(norm_delta_stat(green_PEG,2),norm_delta_stat(green_DART,2));
% [h3, p3]= ttest2(norm_delta_stat(green_PEG,3),norm_delta_stat(green_DART,3));
% %correct for 3 tests
% p1*3
% p2*3
% p3*3
% clear h1 p1 h2 p2 h3 p3
%% plot contrast response for fully matched cells

conResp_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_stat = cell(1,nd); %same for red
conResp_green_se_stat = cell(1,nd); %this will be the se across all green cells
conResp_red_se_stat = cell(1,nd); %same for red

for iDrug = 1:2
    
    thisDrug = find(drug == iDrug);
    green_thisDrug = intersect(thisDrug,haveRunning_green);
    red_thisDrug = intersect(thisDrug,haveRunning_red);
        
    conResp_green_avrg_stat{iDrug}=mean(norm_delta_stat(green_thisDrug,:),1,'omitmissing');
    green_std=std(norm_delta_stat(green_thisDrug,:),1,'omitmissing');
    conResp_green_se_stat{iDrug}=green_std/sqrt(length(green_thisDrug));
    
    conResp_red_avrg_stat{iDrug}=mean(norm_delta_stat(red_thisDrug,:),1,'omitmissing');
    red_std=std(norm_delta_stat(red_thisDrug,:),1,'omitmissing');
    conResp_red_se_stat{iDrug}=red_std/sqrt(length(red_thisDrug));
    
    clear green_std red_std
 
end



figure
subplot(1,2,1) %for the second day
errorbar(cons,conResp_red_avrg_stat{PEG},conResp_red_se_stat{PEG},'color',"#EA9010");
hold on
errorbar(cons,conResp_red_avrg_stat{DART},conResp_red_se_stat{DART},'b');
% title('SST')
% ylabel('Post - pre/post+pre')
% xlabel('contrast') 
xlim([0 1.2])
ylim([-.4 .2])
xticks([.25 .5 1])
% set(gca, 'TickDir', 'out')
box off

subplot(1,2,2) %for the first day
errorbar(cons,conResp_green_avrg_stat{PEG},conResp_green_se_stat{PEG},'--','color',"#EA9010");
hold on
errorbar(cons,conResp_green_avrg_stat{DART},conResp_green_se_stat{DART},'--b');
% title('Pyr')
% xlabel('contrast') 
% set(gca, 'TickDir', 'out')
xlim([0 1.2])
ylim([-.4 .2])
xticks([.25 .5 1])
box off
x0=5;
y0=5;
width=2.2;
height=1.6;
set(gcf,'units','inches','position',[x0,y0,width,height])
sgtitle('stationary')

green_PEG = intersect(find(drug==2),haveRunning_green);
green_DART = intersect(find(drug==1),haveRunning_green);
red_PEG = intersect(find(drug==2),haveRunning_red);
red_DART = intersect(find(drug==1),haveRunning_red);

%two-sample t-tests for red, stationary
[h1, p1]= ttest2(norm_delta_stat(red_PEG,1),norm_delta_stat(red_DART,1));
[h2, p2]= ttest2(norm_delta_stat(red_PEG,2),norm_delta_stat(red_DART,2));
[h3, p3]= ttest2(norm_delta_stat(red_PEG,3),norm_delta_stat(red_DART,3));
%correct for 3 tests
p1*3
p2*3
p3*3
clear h1 p1 h2 p2 h3 p3

%two-sample t-tests for green, stationary
[h1, p1]= ttest2(norm_delta_stat(green_PEG,1),norm_delta_stat(green_DART,1));
[h2, p2]= ttest2(norm_delta_stat(green_PEG,2),norm_delta_stat(green_DART,2));
[h3, p3]= ttest2(norm_delta_stat(green_PEG,3),norm_delta_stat(green_DART,3));
%correct for 3 tests
p1*3
p2*3
p3*3
clear h1 p1 h2 p2 h3 p3
%%
% for running
conResp_green_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_loc = cell(1,nd); %same for red
conResp_green_se_loc = cell(1,nd); %this will be the se across all green cells
conResp_red_se_loc = cell(1,nd); %same for red

for iDrug = 1:2
    
    thisDrug = find(drug == iDrug);
    green_thisDrug = intersect(thisDrug,haveRunning_green);
    red_thisDrug = intersect(thisDrug,haveRunning_red);
        
    conResp_green_avrg_loc{iDrug}=mean(norm_delta_loc(green_thisDrug,:),1,'omitmissing');
    green_std=std(norm_delta_loc(green_thisDrug,:),1,'omitmissing');
    conResp_green_se_loc{iDrug}=green_std/sqrt(length(green_thisDrug));
    
    conResp_red_avrg_loc{iDrug}=mean(norm_delta_loc(red_thisDrug,:),1,'omitmissing');
    red_std=std(norm_delta_loc(red_thisDrug,:),1,'omitmissing');
    conResp_red_se_loc{iDrug}=red_std/sqrt(length(red_thisDrug));
    
    clear green_std red_std
 
end

figure
subplot(1,2,1) %for the second day
errorbar(cons,conResp_red_avrg_loc{PEG},conResp_red_se_loc{PEG},'color',"#EA9010");
hold on
errorbar(cons,conResp_red_avrg_loc{DART},conResp_red_se_loc{DART},'b');
% title('SST')
% ylabel('Post - pre/post+pre')
% xlabel('contrast') 
xlim([0 1.2])
ylim([-.25 .25])
xticks([.25 .5 1])
% set(gca, 'TickDir', 'out')
box off

subplot(1,2,2) %for the first day
errorbar(cons,conResp_green_avrg_loc{PEG},conResp_green_se_loc{PEG},'--','color',"#EA9010");
hold on
errorbar(cons,conResp_green_avrg_loc{DART},conResp_green_se_loc{DART},'--b');
% title('Pyr')
% xlabel('contrast') 
% set(gca, 'TickDir', 'out')
xlim([0 1.2])
ylim([-.25 .25])
xticks([.25 .5 1])
box off

x0=5;
y0=5;
width=2.2;
height=1.6;
set(gcf,'units','inches','position',[x0,y0,width,height])
sgtitle('running')

%two-sample t-tests for red, running 
[h1, p1]= ttest2(norm_delta_loc(red_PEG,1),norm_delta_loc(red_DART,1));
[h2, p2]= ttest2(norm_delta_loc(red_PEG,2),norm_delta_loc(red_DART,2));
[h3, p3]= ttest2(norm_delta_loc(red_PEG,3),norm_delta_loc(red_DART,3));
%correct for 3 tests
p1*3
p2*3
p3*3
clear h1 p1 h2 p2 h3 p3

%two-sample t-tests for green, running
[h1, p1]= ttest2(norm_delta_loc(green_PEG,1),norm_delta_loc(green_DART,1));
[h2, p2]= ttest2(norm_delta_loc(green_PEG,2),norm_delta_loc(green_DART,2));
[h3, p3]= ttest2(norm_delta_loc(green_PEG,3),norm_delta_loc(green_DART,3));
%correct for 3 tests
p1*3
p2*3
p3*3
clear h1 p1 h2 p2 h3 p3



%% make data table to do ANOVA for norm_delta metric
cellID=(1:nKeep_total)';
drugCol=categorical(drug);
norm_delta_stat_table = array2table(norm_delta_stat,'VariableNames',{'S25','S50','S100'});
labels_table =table(mouseID,cellID,drugCol,'VariableNames',{'mouseID' 'cellID' 'drugCond'});
norm_delta_summary_stat = [labels_table,norm_delta_stat_table];

%beh_state_col = repmat('loc',nKeep_total,1);
norm_delta_loc_table = array2table(norm_delta_loc,'VariableNames',{'L25','L50','L100'});
labels_table =table(mouseID,cellID,drugCol,'VariableNames',{'mouseID' 'cellID' 'drugCond'});
norm_delta_summary_loc = [labels_table,norm_delta_loc_table];

norm_delta_full =horzcat(norm_delta_summary_stat,norm_delta_loc_table);

norm_delta_full_SST = norm_delta_full(haveRunning_red,:);
norm_delta_stat_SST = norm_delta_summary_stat(haveRunning_red,:);
norm_delta_stat_pyr = norm_delta_summary_stat(haveRunning_green,:);
norm_delta_loc_SST = norm_delta_summary_loc(haveRunning_red,:);
norm_delta_loc_Pyr = norm_delta_summary_loc(haveRunning_green,:);
norm_delta_stat_unMatchedSST = norm_delta_summary_stat(red_ind_concat,:);
norm_delta_stat_unMatchedPyr = norm_delta_summary_stat(green_ind_concat,:);



%% run ANOVA for norm_delta metric
%full model
WithinModel = table(categorical([1 1 1 2 2 2 ]'), categorical([1 2 3 1 2 3]'), 'VariableNames', {'behState', 'contrast'}); % within-design
fullDeltaModel_SST = fitrm(norm_delta_full_SST, 'S25-L100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(fullDeltaModel_SST, 'withinmodel', 'behState*contrast');
% Display the results
disp(ranovatbl);

%stationary model for matched cells
WithinModel = table(categorical([1 2 3]'), 'VariableNames', {'contrast'}); % within-design
statDeltaModel_SST = fitrm(norm_delta_stat_SST, 'S25-S100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(statDeltaModel_SST, 'withinmodel', 'contrast');
% Display the results
disp(ranovatbl);

%running model for matched cells
WithinModel = table(categorical([1 2 3]'), 'VariableNames', {'contrast'}); % within-design
locDeltaModel_SST = fitrm(norm_delta_loc_SST, 'L25-L100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(locDeltaModel_SST, 'withinmodel', 'contrast');
% Display the results
disp(ranovatbl);

%stationary model for un-matched cells
WithinModel = table(categorical([1 2 3]'), 'VariableNames', {'contrast'}); % within-design
statDeltaModel_SST_unmatched = fitrm(norm_delta_stat_unMatchedSST, 'S25-S100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(statDeltaModel_SST_unmatched, 'withinmodel', 'contrast');
% Display the results
disp(ranovatbl);

%stationary model for un-matched pyr
WithinModel = table(categorical([1 2 3]'), 'VariableNames', {'contrast'}); % within-design
statDeltaModel_PYR_unmatched = fitrm(norm_delta_stat_unMatchedPyr, 'S25-S100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(statDeltaModel_PYR_unmatched, 'withinmodel', 'contrast');
% Display the results
disp(ranovatbl);


%stationary model for matched Pyr cells
WithinModel = table(categorical([1 2 3]'), 'VariableNames', {'contrast'}); % within-design
statDeltaModel_Pyr = fitrm(norm_delta_stat_pyr, 'S25-S100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(statDeltaModel_Pyr, 'withinmodel', 'contrast');
% Display the results
disp(ranovatbl);


%running model for matched Pyr cells
WithinModel = table(categorical([1 2 3]'), 'VariableNames', {'contrast'}); % within-design
locDeltaModel_Pyr = fitrm(norm_delta_loc_Pyr, 'L25-L100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(locDeltaModel_Pyr, 'withinmodel', 'contrast');
% Display the results
disp(ranovatbl);
%% norm_delta contrast value for weakly and strongly correlated cells

red_DART = intersect(red_ind_concat,find(drug));
red_PEG= intersect(red_ind_concat,find(~drug));

highRInds = find(noiseCorr_OG_concat{pre}(1,:)>0.5);
lowRInds = find(noiseCorr_OG_concat{pre}(1,:)<=0.5);
redHigh=intersect(highRInds, red_ind_concat);
redLow=intersect(lowRInds, red_ind_concat);

conResp_redHigh_avrg_stat = cell(1,nd); %this will be the average across all redHigh cells - a single line
conResp_redLow_avrg_stat = cell(1,nd); %same for redLow
conResp_redHigh_se_stat = cell(1,nd); %this will be the se across all redHigh cells
conResp_redLow_se_stat = cell(1,nd); %same for redLow

for iDrug = 1:2
    
    thisDrug = find(drug == iDrug-1);
    lowThisDrug = intersect(thisDrug,redLow);
    highThisDrug = intersect(thisDrug,redHigh);
    length(lowThisDrug)
    length(highThisDrug)

    conResp_redHigh_avrg_stat{iDrug}=mean(norm_delta_stat(highThisDrug,:),1,'omitnan');
    redHigh_std=std(norm_delta_stat(highThisDrug,:),1,'omitnan');
    conResp_redHigh_se_stat{iDrug}=redHigh_std/sqrt(length(highThisDrug));
    
    conResp_redLow_avrg_stat{iDrug}=mean(norm_delta_stat(lowThisDrug,:),1,'omitnan');
    redLow_std=std(norm_delta_stat(lowThisDrug,:),1,'omitnan');
    conResp_redLow_se_stat{iDrug}=redLow_std/sqrt(length(lowThisDrug));
    
    clear redHigh_std redLow_std
 
end

PEG = 1;
DART = 2;
figure
subplot(1,2,1) %
errorbar(cons,conResp_redLow_avrg_stat{PEG},conResp_redLow_se_stat{PEG},'color',"#EA9010");
hold on
errorbar(cons,conResp_redLow_avrg_stat{DART},conResp_redLow_se_stat{DART},'b');
% title(['Weak Corr'])
% xlabel('contrast') 
% ylabel('dF/F, pref ori') 
xlim([0 1.2])
ylim([-.6 .2])
xticks([.25 .5 1])
box off

subplot(1,2,2) %
errorbar(cons,conResp_redHigh_avrg_stat{PEG},conResp_redHigh_se_stat{PEG},'color',"#EA9010");
hold on
errorbar(cons,conResp_redHigh_avrg_stat{DART},conResp_redHigh_se_stat{DART},'b');
% title(['Strong Corr'])
% 
% xlabel('contrast') 
% set(gca, 'TickDir', 'out')
xlim([0 1.2])
ylim([-.6 .2])
xticks([.25 .5 1])
box off

x0=5;
y0=5;
width=2.2;
height=1.6;
set(gcf,'units','inches','position',[x0,y0,width,height])
%% anova for weakly and strongly correlated cells


%stationary model for matched cells
WithinModel = table(categorical([1 2 3]'), 'VariableNames', {'contrast'}); % within-design
corrDeltaModel_SST = fitrm(norm_delta_stat_unMatchedSST, 'S25-S100 ~drugCond*highR', 'WithinDesign', WithinModel);
ranovatbl = ranova(corrDeltaModel_SST, 'withinmodel', 'contrast');
% Display the results
disp(ranovatbl);

%unpaired t-test for norm_delta vs drug in each corr group at 100% contrast

[h1, p1]= ttest2(norm_delta_stat(intersect(red_PEG,highRInds),3),norm_delta_stat(intersect(red_DART,highRInds),3));


%stationary for weak corr
WithinModel = table(categorical([1 2 3]'), 'VariableNames', {'contrast'}); % within-design
weakCorrDeltaModel_SST = fitrm(norm_deltat_stat_weakCorr, 'S25-S100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(weakCorrDeltaModel_SST, 'withinmodel', 'contrast');
% Display the results
disp(ranovatbl);


%stationary for strong corr
WithinModel = table(categorical([1 2 3]'), 'VariableNames', {'contrast'}); % within-design
strongCorrDeltaModel_SST = fitrm(norm_deltat_stat_strongCorr, 'S25-S100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(strongCorrDeltaModel_SST, 'withinmodel', 'contrast');
% Display the results
disp(ranovatbl);

