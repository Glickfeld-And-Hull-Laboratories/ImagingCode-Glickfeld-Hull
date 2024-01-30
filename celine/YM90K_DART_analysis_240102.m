%% read in and concatenate data

clear all; clear global; close all
clc

% ENTER DATASET NAME
load('ds_YM90K_DART.mat');
ds_name='YM90K_DART';
nSess=length(ds_YM90K_DART);

% setting key variables and output folder
%hard-coded variables specific to YM90K-DART/PEG SST project
nd=2;
pre=2;
post=1;
targetCon = [.25 .5 1];
frame_rate = 15;

if computer == 'GLNXA64'
    isilonName =  '/home/cc735@dhe.duke.edu/GlickfeldLabShare';
    base = fullfile(isilonName, '/All_Staff/home/ACh/Data/2p_data');
    
else
    isilonName = 'Z:';
    base = fullfile(isilonName, '\home\ACh\Analysis\2p_analysis');
      
end

d=string(datetime('today'));
fnout= fullfile(base,ds_name,d);

mkdir(fnout);
cd(fnout)
clear d sess_title
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
% concatenating data
nCon = length(targetCon)

mice={};
tc_trial_avrg_stat_concat=cell(1,nd);
tc_trial_avrg_stat_hiPupil_concat=cell(1,nd);
tc_trial_avrg_stat_lowPupil_concat=cell(1,nd);
tc_trial_avrg_loc_concat=cell(1,nd);
resp_keep_concat=cell(1,nd);
resp_max_keep_concat=cell(1,nd);
pref_responses_loc_concat=cell(1,nd);
pref_responses_stat_concat=cell(1,nd);
pref_responses_stat_hiPupil_concat=cell(1,nd);
pref_responses_stat_lowPupil_concat=cell(1,nd);
RIx_concat=cell(1,nd);
dirs_concat=[];
cons_concat=[];
green_concat=[];
red_concat=[];
dfof_max_diff_concat=[];
nKeep_concat=[];
LMI_concat = cell(1,nd);
data_resp_concat = cell(1,nd);
red_fluor_concat=[];
green_fluor_concat=[];
wheel_corr_concat=cell(1,nd);
meanF_concat=cell(1,nd);
norm_dir_resp_stat_concat = cell(1,nd);
norm_dir_resp_loc_concat = cell(1,nd);
pref_nonPref_stat_concat=cell(1,nd);
pref_nonPref_loc_concat=cell(1,nd);
pref_dir_concat=cell(1,nd);
noiseCorr_concat = cell(1,nd);
sigCorr_concat = cell(1,nd);
pref_allTrials_stat_concat =cell(nCon,nd);
pref_allTrials_loc_concat =cell(nCon,nd);
dataTableConat=[];
drug=cell(1,nSess);
pupilMeans_concat=nan(nd,2,nSess);
pupilCounts_concat=nan(nd,2,nSess);

nonPref_trial_avrg_stat_concat=cell(1,nd);
nonPref_trial_avrg_loc_concat=cell(1,nd);

responseByCondProps_concat=nan(6,2,nSess);

cellID_adjustment=0;
for iSess = 1:nSess
    mouse = ds_YM90K_DART(iSess).mouse;
    mice=[mice;mouse];
    thisDrug = ds_YM90K_DART(iSess).drug;
    drug{iSess}=thisDrug;

    if iSess > 1
        cellID_adjustment=max(temp_table.cell_ID_unique); %this should get saved until the next loop;
    end

    if ds_YM90K_DART(iSess).multiday_timesincedrug_hours>0
        dart_str = [ds_YM90K_DART(iSess).drug '_' num2str(ds_YM90K_DART(iSess).multiday_timesincedrug_hours) 'Hr'];
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
    load(fullfile(fn_multi,'pupilMeans.mat'));

    temp_table =readtable(fullfile(fn_multi,'dataTable.csv'));

    temp_table.z_speed=zscor_xnan(temp_table.speed);
    temp_table.z_pupil=zscor_xnan(temp_table.pupil);
    temp_table.cell_ID_unique=temp_table.cellID + cellID_adjustment;

    dataTableConat=[dataTableConat; temp_table];

    nKeep = size(tc_trial_avrg_stat{post},2);

    pupilMeans_concat(:,:,iSess)=pupilMeans;
    pupilCounts_concat(:,:,iSess)=pupilCounts;


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
    sharedCon=find(ismember(cons, targetCon));

    nOn = input(1).nScansOn;
    nOff = input(1).nScansOff;
    %start conatenating
    dirs_concat = [dirs_concat,dirs]; 
    cons_concat = [cons_concat,cons(sharedCon)];
    red_concat = [red_concat, red_keep_logical];
    green_concat = [green_concat, green_keep_logical];
    nKeep_concat = [nKeep_concat,nKeep];
    responseByCondProps_concat(:,:,iSess)=responseByCondProps;

    clear cons
    
    
    for id = 1:nd
        
        tc_trial_avrg_stat_concat{id} =cat(2,tc_trial_avrg_stat_concat{id},tc_trial_avrg_stat{id}(:,:,sharedCon));
        tc_trial_avrg_stat_hiPupil_concat{id} = cat(2,tc_trial_avrg_stat_hiPupil_concat{id},tc_trial_avrg_stat_hiPupil{id}(:,:,sharedCon));
        tc_trial_avrg_stat_lowPupil_concat{id} = cat(2,tc_trial_avrg_stat_lowPupil_concat{id},tc_trial_avrg_stat_lowPupil{id}(:,:,sharedCon));
        tc_trial_avrg_loc_concat{id} =cat(2,tc_trial_avrg_loc_concat{id},tc_trial_avrg_loc{id}(:,:,sharedCon));

        nonPref_trial_avrg_stat_concat{id} =cat(2,nonPref_trial_avrg_stat_concat{id},nonPref_trial_avrg_stat{id}(:,:,sharedCon));
        nonPref_trial_avrg_loc_concat{id} =cat(2,nonPref_trial_avrg_loc_concat{id},nonPref_trial_avrg_loc{id}(:,:,sharedCon));

        resp_keep_concat{id}=cat(1,resp_keep_concat{id},resp_keep{id});
        resp_max_keep_concat{id}=cat(1,resp_max_keep_concat{id},resp_max_keep{id}(:,sharedCon));
        LMI_concat{id}=cat(1,LMI_concat{id},LMI{id}(:,sharedCon));
        pref_responses_loc_concat{id}=cat(1,pref_responses_loc_concat{id},pref_responses_loc{id}(:,sharedCon));
        pref_responses_stat_concat{id}=cat(1,pref_responses_stat_concat{id},pref_responses_stat{id}(:,sharedCon));
        pref_responses_stat_hiPupil_concat{id}=cat(1,pref_responses_stat_hiPupil_concat{id},pref_responses_stat_hiPupil{id}(:,sharedCon));
        pref_responses_stat_lowPupil_concat{id}=cat(1,pref_responses_stat_lowPupil_concat{id},pref_responses_stat_lowPupil{id}(:,sharedCon));
        RIx_concat{id}=cat(1,RIx_concat{id},sum(RIx{id}));
        wheel_corr_concat{id}=cat(2,wheel_corr_concat{id},wheel_corr{id});
        meanF=mean(fullTC_keep{id},1);
        meanF_concat{id}=cat(2,meanF_concat{id}, meanF);
        norm_dir_resp_stat_concat{id}=cat(1,norm_dir_resp_stat_concat{id},norm_dir_resp_stat{id});
        norm_dir_resp_loc_concat{id}=cat(1,norm_dir_resp_loc_concat{id},norm_dir_resp_loc{id});
        pref_nonPref_stat_concat{id}=cat(1,pref_nonPref_stat_concat{id},pref_nonPref_stat{id});
        pref_nonPref_loc_concat{id}=cat(1,pref_nonPref_loc_concat{id},pref_nonPref_loc{id});
        pref_dir_concat{id}=cat(2,pref_dir_concat{id},pref_dir_keep{id});
        noiseCorr_concat{id}=cat(2,noiseCorr_concat{id},noiseCorr{id});
        sigCorr_concat{id}=cat(2,sigCorr_concat{id},sigCorr{id});
        for i = 1:length(sharedCon)
            iCon=sharedCon(i);
            pref_allTrials_stat_concat{i,id}=[pref_allTrials_stat_concat{i,id},pref_allTrials_stat{iCon,id}];
            pref_allTrials_loc_concat{i,id}=[pref_allTrials_loc_concat{i,id},pref_allTrials_loc{iCon,id}];
        end
        clear meanF i
    end
    dfof_max_diff_concat=cat(1,dfof_max_diff_concat,dfof_max_diff(:,sharedCon));
   green_fluor_concat=cat(2,green_fluor_concat,green_fluor_keep);
   red_fluor_concat=cat(2,red_fluor_concat,red_fluor_keep);
    
iSess
end
%
clear mouse iSess nKeep iSess fn_multi cons oris pupilMeans norm_dir_resp_loc norm_dir_resp_stat
clear explanation1 resp_keep sig_diff pref_con_keep pref_ori_keep tOri_match tCon_match data_trial_keep nTrials tc_trial_avrg_keep green_keep_logical red_keep_logical green_ind_keep red_ind_keep
clear LMI RIx locCounts locResp locTCs statResp statTCs wheel_tc ttest_results_stat ttest_results_loc ttest_results_allCon_stat ttest_results_allCon_loc
clear data_con_resp_keep data_ori_resp_keep data_rep_keep dfof_max_diff dfof_max_diff_raw explanation2 resp_max_keep data_resp_keep pref_responses_stat pref_responses_loc
clear tc_trial_avrg_stat tc_trial_avrg_loc fullTC_keep norm_dir_resp sigCorr noiseCorr responseByCondProps
clear red_fluor_all red_fluor_match green_fluor_match green_fluor_match red_fluor_keep green_fluor_keep R_p_values pre_nonPref_stat pre_nonPref_loc
red_ind_concat = find(red_concat);
green_ind_concat = find(green_concat);
%
cons=unique(cons_concat);

dirs=unique(dirs_concat);
oris=dirs-180; oris=oris(oris>=0);
nDir=length(dirs);
nKeep_total = sum(nKeep_concat);
mean(RIx_concat{pre})
mean(RIx_concat{post})

%find cells that have running data on both days
haveRunning_pre = ~isnan(pref_responses_loc_concat{pre});

haveRunning_post = ~isnan(pref_responses_loc_concat{post});
for iCon =1:nCon
haveRunning_both{iCon}= find(haveRunning_pre(:,iCon).* haveRunning_post(:,iCon));
haveRunning_green{iCon} = intersect(haveRunning_both{iCon}, green_ind_concat);
haveRunning_red{iCon} = intersect(haveRunning_both{iCon}, red_ind_concat);
end

clear haveRunning_pre haveRunning_post haveRunning_both


%get an array with indices for each mouse
mouseInds=cell(1,nSess);
start=1;
for iMouse = 1:nSess
    mouseInds{iMouse}=start:(start-1)+nKeep_concat(iMouse);
    start = start+nKeep_concat(iMouse)
end
clear start iMouse


%find how many haveRunning red cells exist for each mouse
cellCountsRed = nan(nSess,nCon);
mouseNames=[];
for iMouse = 1:nSess
    for iCon = 1:nCon
        cellCountsRed(iMouse, iCon,1)=length(intersect(haveRunning_red{iCon},(mouseInds{iMouse})));
        
    end
    mouseNames=[mouseNames, string(mice(iMouse,:))]
end
clear  iMouse


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
have_HI_pre = ~isnan(pref_responses_stat_hiPupil_concat{pre});
have_HI_post = ~isnan(pref_responses_stat_hiPupil_concat{post});

have_LOW_pre = ~isnan(pref_responses_stat_lowPupil_concat{pre});
have_LOW_post = ~isnan(pref_responses_stat_lowPupil_concat{post});

have_bothPupil=cell(1,3);
for iCon =1:nCon
    have_HI_both= find(have_HI_pre(:,iCon).* have_HI_post(:,iCon));
    have_LOW_both=find(have_LOW_pre(:,iCon).* have_LOW_post(:,iCon));
    have_bothPupil{iCon}=intersect(have_HI_both,have_LOW_both);
end
 
clear have_HI_pre have_HI_post have_LOW_pre have_LOW_post have_HI_both have_LOW_both

have_allPupil = intersect(have_bothPupil{1},have_bothPupil{2});
have_allPupil = intersect(have_allPupil, have_bothPupil{3});


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
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(1,2,1) %for the first day



ylim([-.02 .25]);;
hold on
shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon),tc_green_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon),tc_green_se_stat{post}(:,iCon),'b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['Pyr',' n = ', num2str(length(green_ind_concat))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
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
title(['SST',' n = ', num2str(length(red_ind_concat))])

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['stationary, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) '_stat_allKeep_timecourses.pdf']),'-dpdf');
end  

%% ALL KEEP CELLS at non-preferred directions

nonPref_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
nonPref_red_avrg_stat = cell(1,nd); %same for red
nonPref_green_se_stat = cell(1,nd); %this will be the se across all green cells
nonPref_red_se_stat = cell(1,nd); %same for red



for id = 1:nd
    for iCon=1:nCon
        
    nonPref_green_avrg_stat{id}(:,iCon)=nanmean(nonPref_trial_avrg_stat_concat{id}(:,green_ind_concat,iCon),2);
    green_std=nanstd(nonPref_trial_avrg_stat_concat{id}(:,green_ind_concat,iCon),[],2);
    nonPref_green_se_stat{id}(:,iCon)=green_std/sqrt(length(green_ind_concat));
    
    nonPref_red_avrg_stat{id}(:,iCon)=nanmean(nonPref_trial_avrg_stat_concat{id}(:,red_ind_concat,iCon),2);
    red_std=nanstd(nonPref_trial_avrg_stat_concat{id}(:,red_ind_concat,iCon),[],2);
    nonPref_red_se_stat{id}(:,iCon)=red_std/sqrt(length(red_ind_concat));
    
    clear green_std red_std
    end
end
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(nonPref_green_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(1,2,1) %for the first day



ylim([-.02 .25]);;
hold on
shadedErrorBar(t,nonPref_green_avrg_stat{pre}(:,iCon),nonPref_green_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,nonPref_green_avrg_stat{post}(:,iCon),nonPref_green_se_stat{post}(:,iCon),'b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['Pyr',' n = ', num2str(length(green_ind_concat))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
shadedErrorBar(t,nonPref_red_avrg_stat{pre}(:,iCon),nonPref_red_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,nonPref_red_avrg_stat{post}(:,iCon),nonPref_red_se_stat{post}(:,iCon),'b');
ylim([-.02 .25]);;
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['SST',' n = ', num2str(length(red_ind_concat))])

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['stationary, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) 'nonPref_stat_allKeep_timecourses.pdf']),'-dpdf');
end 


%% pupil stats
figure

for iSess = 1:nSess
    subplot(2,2,1)
    plot(pupilMeans_concat(pre,:,iSess),'o-')
    xlim([.5 2.5])
    box off
    set(gca,'TickDir','out')
    xticks([1 2])
    xticklabels({'Large pupil','Small pupil'})
    ylabel('Pupil diameter')
    hold on
    title('Pre-DART')
    ylim([0.15 0.45])

    subplot(2,2,2)
    plot(pupilMeans_concat(post,:,iSess),'o-')
    xlim([.5 2.5])
    box off
    set(gca,'TickDir','out')
    xticks([1 2])
    xticklabels({'Large pupil','Small pupil'})
    hold on
    title('Post-DART')
    ylim([0.15 0.45])

    subplot(2,2,3)
    plot(pupilMeans_concat([pre post],1,iSess),'o-')
    xlim([.5 2.5])
    box off
    set(gca,'TickDir','out')
    xticks([1 2])
    xticklabels({'Pre-DART','Post-DART'})
    ylabel('Pupil diameter')
    hold on
    title('Large pupil')
    ylim([0.15 0.45])

    subplot(2,2,4)
    plot(pupilMeans_concat([pre post],2,iSess),'o-')
    xlim([.5 2.5])
    box off
    set(gca,'TickDir','out')
    xticks([1 2])
    xticklabels({'Pre-DART','Post-DART'})
    hold on
    title('Small pupil')
    ylim([0.15 0.45])
end
print(fullfile(fnout,['pupilStats.pdf']),'-dpdf');



%% plot large pupil TCs

tc_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_stat = cell(1,nd); %same for red
tc_green_se_stat = cell(1,nd); %this will be the se across all green cells
tc_red_se_stat = cell(1,nd); %same for red

temp_green = intersect(green_ind_concat,have_bothPupil{iCon});

for id = 1:nd
    for iCon=1:nCon
    temp_green = intersect(green_ind_concat,have_bothPupil{iCon});
    temp_red = intersect(red_ind_concat,have_bothPupil{iCon});

    tc_green_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_hiPupil_concat{id}(:,temp_green,iCon),2);
    green_std=nanstd(tc_trial_avrg_stat_hiPupil_concat{id}(:,temp_green,iCon),[],2);
    tc_green_se_stat{id}(:,iCon)=green_std/sqrt(length(temp_green));
    
    tc_red_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_hiPupil_concat{id}(:,temp_red,iCon),2);
    red_std=nanstd(tc_trial_avrg_stat_hiPupil_concat{id}(:,temp_red,iCon),[],2);
    tc_red_se_stat{id}(:,iCon)=red_std/sqrt(length(temp_red));
    
    clear green_std red_std
    end
end
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon

temp_green = intersect(green_ind_concat,have_bothPupil{iCon});
temp_red = intersect(red_ind_concat,have_bothPupil{iCon});

figure
subplot(1,2,1) %for the first day



ylim([-.02 .25]);;
hold on
shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon),tc_green_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon),tc_green_se_stat{post}(:,iCon),'b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['Pyr',' n = ', num2str(length(temp_green))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
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

print(fullfile(fnout,[num2str(cons(iCon)) '_stat_hiPupil_timecourses.pdf']),'-dpdf');
end  

%% plot small pupil TCs

tc_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_stat = cell(1,nd); %same for red
tc_green_se_stat = cell(1,nd); %this will be the se across all green cells
tc_red_se_stat = cell(1,nd); %same for red

temp_green = intersect(green_ind_concat,have_bothPupil{iCon});

for id = 1:nd
    for iCon=1:nCon
    temp_green = intersect(green_ind_concat,have_bothPupil{iCon});
    temp_red = intersect(red_ind_concat,have_bothPupil{iCon});

    tc_green_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_lowPupil_concat{id}(:,temp_green,iCon),2);
    green_std=nanstd(tc_trial_avrg_stat_lowPupil_concat{id}(:,temp_green,iCon),[],2);
    tc_green_se_stat{id}(:,iCon)=green_std/sqrt(length(temp_green));
    
    tc_red_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_lowPupil_concat{id}(:,temp_red,iCon),2);
    red_std=nanstd(tc_trial_avrg_stat_lowPupil_concat{id}(:,temp_red,iCon),[],2);
    tc_red_se_stat{id}(:,iCon)=red_std/sqrt(length(temp_red));
    
    clear green_std red_std
    end
end
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon

temp_green = intersect(green_ind_concat,have_bothPupil{iCon});
temp_red = intersect(red_ind_concat,have_bothPupil{iCon});

figure
subplot(1,2,1) %for the first day



ylim([-.02 .25]);;
hold on
shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon),tc_green_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon),tc_green_se_stat{post}(:,iCon),'b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['Pyr',' n = ', num2str(length(temp_green))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
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

print(fullfile(fnout,[num2str(cons(iCon)) '_stat_lowPupil_timecourses.pdf']),'-dpdf');
end  
%% for the cells that have stationary and running

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
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(1,2,1) %for the first day



ylim([-.02 .3]);;
hold on
shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon),tc_green_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon),tc_green_se_stat{post}(:,iCon),'b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['Pyr',' n = ', num2str(length(haveRunning_green{iCon}))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_stat{pre}(:,iCon),tc_red_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_stat{post}(:,iCon),tc_red_se_stat{post}(:,iCon),'b');
ylim([-.02 .3]);;
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['SST',' n = ', num2str(length(haveRunning_red{iCon}))])

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['stationary, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) '_stat_timecourses.pdf']),'-dpdf');
end 



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


%creat a time axis in seconds
t=1:(size(tc_green_avrg_loc{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);


for iCon = 1:nCon
figure
subplot(1,2,1) %for the first day
shadedErrorBar(t,tc_green_avrg_loc{pre}(:,iCon),tc_green_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_loc{post}(:,iCon),tc_green_se_loc{post}(:,iCon),'b');
ylim([-.02 .5]);;
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
hold on
title(['Pyr',' n = ', num2str(length(haveRunning_green{iCon}))])
ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_loc{pre}(:,iCon),tc_red_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_loc{post}(:,iCon),tc_red_se_loc{post}(:,iCon),'b');
ylim([-.02 .5]);;
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
hold on
ylabel('dF/F') 
xlabel('s') 
title(['SST',' n = ', num2str(length(haveRunning_red{iCon}))])
x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['running, contrast = ' num2str(cons(iCon))])
print(fullfile(fnout,[num2str(cons(iCon)) '_loc_timecourses.pdf']),'-dpdf');
end 
clear txt1 txt2

%% shaded error tc for loc non-preferred

nonPref_green_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
nonPref_red_avrg_loc = cell(1,nd); %same for red
nonPref_green_se_loc = cell(1,nd); %this will be the se across all green cells
nonPref_red_se_loc = cell(1,nd); %same for red

for id = 1:nd
    for iCon=1:nCon
    nonPref_green_avrg_loc{id}(:,iCon)=nanmean(nonPref_trial_avrg_loc_concat{id}(:,haveRunning_green{iCon},iCon),2);
    green_std=nanstd(nonPref_trial_avrg_loc_concat{id}(:,haveRunning_green{iCon},iCon),[],2);
    nonPref_green_se_loc{id}(:,iCon)=green_std/sqrt(length(haveRunning_green{iCon}));
    
    nonPref_red_avrg_loc{id}(:,iCon)=nanmean(nonPref_trial_avrg_loc_concat{id}(:,haveRunning_red{iCon},iCon),2);
    red_std=nanstd(nonPref_trial_avrg_loc_concat{id}(:,haveRunning_red{iCon},iCon),[],2);
    nonPref_red_se_loc{id}(:,iCon)=red_std/sqrt(length(haveRunning_red{iCon}));
    
    clear green_std red_std
    end
end


%creat a time axis in seconds
t=1:(size(nonPref_green_avrg_loc{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);


for iCon = 1:nCon
figure
subplot(1,2,1) %for the first day
shadedErrorBar(t,nonPref_green_avrg_loc{pre}(:,iCon),nonPref_green_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,nonPref_green_avrg_loc{post}(:,iCon),nonPref_green_se_loc{post}(:,iCon),'b');
ylim([-.02 .5]);;
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
hold on
title(['Pyr',' n = ', num2str(length(haveRunning_green{iCon}))])
ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) %for the second day
shadedErrorBar(t,nonPref_red_avrg_loc{pre}(:,iCon),nonPref_red_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,nonPref_red_avrg_loc{post}(:,iCon),nonPref_red_se_loc{post}(:,iCon),'b');
ylim([-.02 .5]);;
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
hold on
ylabel('dF/F') 
xlabel('s') 
title(['SST',' n = ', num2str(length(haveRunning_red{iCon}))])
x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['running, contrast = ' num2str(cons(iCon))])
print(fullfile(fnout,[num2str(cons(iCon)) 'nonPref_loc_timecourses.pdf']),'-dpdf');
end 
clear txt1 txt2

%% scatterplot of mean response to pref-ori stimuli pre-DART vs post-DART, and each subplot is one cell type X behavioral state
%will probably drop these, byt maybe keep around for performing stats
green_ex_list=[]; %to highlight particular cells
red_ex_list=[];

for iCon = 1:nCon
figure; movegui('center') 
subplot(2,2,1)
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
title('Pyr stationary')
axis square
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
hold off


subplot(2,2,2)
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
title('SST stationary')
axis square
hold off

subplot(2,2,3)
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
title('Pyr running')
axis square
set(gca, 'TickDir', 'out')
uistack(hline,'bottom');
hold off

subplot(2,2,4)
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
title('SST running')
uistack(hline,'bottom');
axis square
hold off
set(gca, 'TickDir', 'out')



sgtitle(num2str(cons(iCon)))
print(fullfile(fnout,[num2str(cons(iCon)) 'scatter_crossDay.pdf']),'-dpdf','-bestfit')
clear mean_pre_stat mean_post_stat stderror_post stderror_pre
end




%% linear direction plot 


dirs_for_plotting=dirs-180;

green_dir_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
red_dir_avrg_stat = cell(1,nd); %same for red
green_dir_se_stat = cell(1,nd); %this will be the se across all green cells
red_dir_se_stat = cell(1,nd); %same for red


    for id = 1:nd
       
        green_dir_avrg_stat{id}=nanmean(nanmean(norm_dir_resp_stat_concat{id}(green_all,:,:),3),1);
        green_std=nanstd(nanmean(norm_dir_resp_stat_concat{id}(green_all,:,:),3),[],1);
        green_dir_se_stat{id}=green_std/sqrt(length(green_all));
        green_dir_avrg_stat{id}=circshift(green_dir_avrg_stat{id},4);
        green_dir_se_stat{id}=circshift(green_dir_se_stat{id},4);
        
        red_dir_avrg_stat{id}=nanmean(nanmean(norm_dir_resp_stat_concat{id}(red_all,:,:),3),1);
        red_std=nanstd(nanmean(norm_dir_resp_stat_concat{id}(red_all,:,:),3),[],1);
        red_dir_se_stat{id}=red_std/sqrt(length(red_all));
        red_dir_avrg_stat{id}=circshift(red_dir_avrg_stat{id},4);
        red_dir_se_stat{id}=circshift(red_dir_se_stat{id},4);
        clear green_std red_std
        
    end
    
    
    
    green_dir_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
    red_dir_avrg_loc = cell(1,nd); %same for red
    green_dir_se_loc = cell(1,nd); %this will be the se across all green cells
    red_dir_se_loc = cell(1,nd); %same for red
    
     for id = 1:nd
       
        green_dir_avrg_loc{id}=nanmean(nanmean(norm_dir_resp_loc_concat{id}(green_all,:,:),3),1);
        green_std=nanstd(nanmean(norm_dir_resp_loc_concat{id}(green_all,:,:),3),[],1);
        green_dir_se_loc{id}=green_std/sqrt(length(green_all));
        green_dir_avrg_loc{id}=circshift(green_dir_avrg_loc{id},4);
        green_dir_se_loc{id}=circshift(green_dir_se_loc{id},4);
        
        red_dir_avrg_loc{id}=nanmean(nanmean(norm_dir_resp_loc_concat{id}(red_all,:,:),3),1);
        red_std=nanstd(nanmean(norm_dir_resp_loc_concat{id}(red_all,:,:),3),[],1);
        red_dir_se_loc{id}=red_std/sqrt(length(red_all));
        red_dir_avrg_loc{id}=circshift(red_dir_avrg_loc{id},4);
        red_dir_se_loc{id}=circshift(red_dir_se_loc{id},4);
        clear green_std red_std
        
    end
    
    
    
    figure
    subplot(2,2,1)
    errorbar(dirs_for_plotting,green_dir_avrg_stat{pre},green_dir_se_stat{pre},'k')
    hold on
    errorbar(dirs_for_plotting,green_dir_avrg_stat{post},green_dir_se_stat{post},'b')
    title(['Stationary, ', num2str(length(green_all)),' fully matched Pyr'])
    set(gca, 'TickDir', 'out')
    axis square
    box off
    ylabel('dF/F')
    %ylim([-0.01 .2])
    
    subplot(2,2,2)
    errorbar(dirs_for_plotting,red_dir_avrg_stat{pre},red_dir_se_stat{pre},'k')
    hold on
    errorbar(dirs_for_plotting,red_dir_avrg_stat{post},red_dir_se_stat{post},'b')
    title(['Stationary, ', num2str(length(red_all)),' fully matched SST'])
    set(gca, 'TickDir', 'out')
    axis square
    box off
    %ylim([-0.01 .2])
    
    subplot(2,2,3)
    errorbar(dirs_for_plotting,green_dir_avrg_loc{pre},green_dir_se_loc{pre},'k')
    hold on
    errorbar(dirs_for_plotting,green_dir_avrg_loc{post},green_dir_se_loc{post},'b')
    title('Running, Pyr')
    set(gca, 'TickDir', 'out')
    axis square
    box off
    xlabel('normalized direction')
    ylabel('dF/F')
   % ylim([-0.01 .3])
    
    
    subplot(2,2,4)
    errorbar(dirs_for_plotting,red_dir_avrg_loc{pre},red_dir_se_loc{pre},'k')
    hold on
    errorbar(dirs_for_plotting,red_dir_avrg_loc{post},red_dir_se_loc{post},'b')
    title('Running, SST')
    set(gca, 'TickDir', 'out')
    axis square
    box off
    xlabel('normalized direction')
    %ylim([-0.01 .3])
    
    sgtitle(['Normalized direction tuning, no rectification'])
    
    
    print(fullfile(fnout,['dirTuning.pdf']),'-dpdf','-bestfit')


%% comparing the preferred direction across days

change_pref = pref_dir_concat{pre}-pref_dir_concat{post};
change_pref=deg2rad(change_pref);
figure
subplot(1,2,1)
polarhistogram(change_pref(green_ind_concat),'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5])
title('Pyr')
set(gca, 'TickDir', 'out')

subplot(1,2,2)
polarhistogram(change_pref(red_ind_concat),'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5])
title('SST')
set(gca, 'TickDir', 'out')
print(fullfile(fnout,['prefDir_change.pdf']),'-dpdf','-bestfit')

%% looking at DART effect on preferred and non-preferred orientations
%preferred = preferred direction and preferred - 180, ie same orientation;
%non-preferred = all other directions
green_means_stat = cell(1,nd);
green_se_stat = cell(1,nd);
green_means_loc = cell(1,nd);
green_se_loc = cell(1,nd);

red_means_stat = cell(1,nd);
red_se_stat = cell(1,nd);
red_means_loc = cell(1,nd);
red_se_loc = cell(1,nd);

for id = 1:nd
    red_means_stat{id}=mean(pref_nonPref_stat_concat{id}(red_ind_concat,:,:),1,'omitmissing');
    red_std_stat=std(pref_nonPref_stat_concat{id}(red_ind_concat,:,:),[],1,'omitmissing');
    red_se_stat{id}=red_std_stat/sqrt(length(red_ind_concat));
    
    red_means_loc{id}=mean(pref_nonPref_loc_concat{id}(red_ind_concat,:,:),1,'omitmissing');
    red_std_loc=std(pref_nonPref_loc_concat{id}(red_ind_concat,:,:),[],1,'omitmissing');
    red_se_loc{id}=red_std_loc/sqrt(length(red_ind_concat));

    green_means_stat{id}=mean(pref_nonPref_stat_concat{id}(green_ind_concat,:,:),1,'omitmissing');
    green_std_stat=std(pref_nonPref_stat_concat{id}(green_ind_concat,:,:),[],1,'omitmissing');
    green_se_stat{id}=green_std_stat/sqrt(length(green_ind_concat));
    
    green_means_loc{id}=mean(pref_nonPref_loc_concat{id}(green_ind_concat,:,:),1,'omitmissing');
    green_std_loc=std(pref_nonPref_loc_concat{id}(green_ind_concat,:,:),[],1,'omitmissing');
    green_se_loc{id}=green_std_loc/sqrt(length(green_ind_concat));
end

figure;
subplot(2,2,1)
errorbar(cons,squeeze(green_means_stat{pre}(1,1,:)),squeeze(green_se_stat{pre}(1,1,:)),'k')
hold on
errorbar(cons,squeeze(green_means_stat{pre}(1,2,:)),squeeze(green_se_stat{pre}(1,2,:)),'--k')
errorbar(cons,squeeze(green_means_stat{post}(1,1,:)),squeeze(green_se_stat{post}(1,1,:)),'b')
errorbar(cons,squeeze(green_means_stat{post}(1,2,:)),squeeze(green_se_stat{post}(1,2,:)),'--b')
title('Stationary, Pyr')
set(gca, 'TickDir', 'out')
axis square
box off
ylabel('dF/F')
xlim([0 1.25])
ylim([-.005 .28])

subplot(2,2,3)
errorbar(cons,squeeze(green_means_loc{pre}(1,1,:)),squeeze(green_se_loc{pre}(1,1,:)),'k')
hold on
errorbar(cons,squeeze(green_means_loc{pre}(1,2,:)),squeeze(green_se_loc{pre}(1,2,:)),'--k')
errorbar(cons,squeeze(green_means_loc{post}(1,1,:)),squeeze(green_se_loc{post}(1,1,:)),'b')
errorbar(cons,squeeze(green_means_loc{post}(1,2,:)),squeeze(green_se_loc{post}(1,2,:)),'--b')
title('Running, Pyr')
set(gca, 'TickDir', 'out')
axis square
box off
ylabel('dF/F')
xlabel('Contrast')
xlim([0 1.25])
ylim([-.005 .28])

subplot(2,2,2)
errorbar(cons,squeeze(red_means_stat{pre}(1,1,:)),squeeze(red_se_stat{pre}(1,1,:)),'k')
hold on
errorbar(cons,squeeze(red_means_stat{pre}(1,2,:)),squeeze(red_se_stat{pre}(1,2,:)),'--k')
errorbar(cons,squeeze(red_means_stat{post}(1,1,:)),squeeze(red_se_stat{post}(1,1,:)),'b')
errorbar(cons,squeeze(red_means_stat{post}(1,2,:)),squeeze(red_se_stat{post}(1,2,:)),'--b')
title('Stationary, SST')
set(gca, 'TickDir', 'out')
axis square
box off
xlim([0 1.25])
ylim([-.005 .28])


subplot(2,2,4)
errorbar(cons,squeeze(red_means_loc{pre}(1,1,:)),squeeze(red_se_loc{pre}(1,1,:)),'k')
hold on
errorbar(cons,squeeze(red_means_loc{pre}(1,2,:)),squeeze(red_se_loc{pre}(1,2,:)),'--k')
errorbar(cons,squeeze(red_means_loc{post}(1,1,:)),squeeze(red_se_loc{post}(1,1,:)),'b')
errorbar(cons,squeeze(red_means_loc{post}(1,2,:)),squeeze(red_se_loc{post}(1,2,:)),'--b')
title('Running, SST')
set(gca, 'TickDir', 'out')
axis square
box off
ylabel('dF/F')
xlabel('Contrast')
xlim([0 1.25])
ylim([-.005 .28])


clear green_means_stat green_means_loc red_mean_stat red_means_loc green_se_stat green_se_loc red_se_stat red_se_loc green_std red_std
% now with post-pre y axis


red_means_stat=mean(pref_nonPref_stat_concat{post}(red_ind_concat,:,:)-pref_nonPref_stat_concat{pre}(red_ind_concat,:,:),1,'omitmissing');
red_std_stat=std(pref_nonPref_stat_concat{post}(red_ind_concat,:,:)-pref_nonPref_stat_concat{pre}(red_ind_concat,:,:),[],1,'omitmissing');
red_se_stat=red_std_stat/sqrt(length(red_ind_concat));

red_means_loc=mean(pref_nonPref_loc_concat{post}(red_ind_concat,:,:)-pref_nonPref_loc_concat{pre}(red_ind_concat,:,:),1,'omitmissing');
red_std_loc=std(pref_nonPref_loc_concat{post}(red_ind_concat,:,:)-pref_nonPref_loc_concat{pre}(red_ind_concat,:,:),[],1,'omitmissing');
red_se_loc=red_std_loc/sqrt(length(red_ind_concat));

green_means_stat=mean(pref_nonPref_stat_concat{post}(green_ind_concat,:,:)-pref_nonPref_stat_concat{pre}(green_ind_concat,:,:),1,'omitmissing');
green_std_stat=std(pref_nonPref_stat_concat{post}(green_ind_concat,:,:)-pref_nonPref_stat_concat{pre}(green_ind_concat,:,:),[],1,'omitmissing');
green_se_stat=green_std_stat/sqrt(length(green_ind_concat));

green_means_loc=mean(pref_nonPref_loc_concat{post}(green_ind_concat,:,:)-pref_nonPref_loc_concat{pre}(green_ind_concat,:,:),1,'omitmissing');
green_std_loc=std(pref_nonPref_loc_concat{pre}(green_ind_concat,:,:)-pref_nonPref_loc_concat{post}(green_ind_concat,:,:),[],1,'omitmissing');
green_se_loc=green_std_loc/sqrt(length(green_ind_concat));


figure;
subplot(2,2,1)
errorbar(cons,squeeze(green_means_stat(1,1,:)),squeeze(green_se_stat(1,1,:)),'k')
hold on
errorbar(cons,squeeze(green_means_stat(1,2,:)),squeeze(green_se_stat(1,2,:)),'--k')
title('Stationary, Pyr')
set(gca, 'TickDir', 'out')
axis square
box off
ylabel('dF/F post-pre')
xlim([0 1.25])
ylim([-.04 .04])

subplot(2,2,3)
errorbar(cons,squeeze(green_means_loc(1,1,:)),squeeze(green_se_loc(1,1,:)),'k')
hold on
errorbar(cons,squeeze(green_means_loc(1,2,:)),squeeze(green_se_loc(1,2,:)),'--k')
title('Running, Pyr')
set(gca, 'TickDir', 'out')
axis square
box off
ylabel('dF/F post-pre')
xlabel('Contrast')
xlim([0 1.25])
ylim([-.04 .04])

subplot(2,2,2)
errorbar(cons,squeeze(red_means_stat(1,1,:)),squeeze(red_se_stat(1,1,:)),'k')
hold on
errorbar(cons,squeeze(red_means_stat(1,2,:)),squeeze(red_se_stat(1,2,:)),'--k')
title('Stationary, SST')
set(gca, 'TickDir', 'out')
axis square
box off
xlim([0 1.25])
ylim([-.04 .04])


subplot(2,2,4)
errorbar(cons,squeeze(red_means_loc(1,1,:)),squeeze(red_se_loc(1,1,:)),'k')
hold on
errorbar(cons,squeeze(red_means_loc(1,2,:)),squeeze(red_se_loc(1,2,:)),'--k')
title('Running, SST')
set(gca, 'TickDir', 'out')
axis square
box off
ylabel('dF/F')
xlabel('Contrast')
xlim([0 1.25])
ylim([-.04 .04])

clear green_means_stat green_means_loc red_mean_stat red_means_loc green_se_stat green_se_loc red_se_stat red_se_loc green_std red_std
print(fullfile(fnout,['conXpref.pdf']),'-dpdf','-bestfit')
%% Define OSI 

osi_stat= cell(1,nd);
osi_loc= cell(1,nd);

norm_dir_resp_stat_concat_rect=cell(1,nd);
norm_dir_resp_loc_concat_rect=cell(1,nd);
for id = 1:nd
    norm_dir_resp_stat_concat_rect{id}=norm_dir_resp_stat_concat{id};
    norm_dir_resp_stat_concat_rect{id}(norm_dir_resp_stat_concat_rect{id}<0)=0;

    norm_dir_resp_loc_concat_rect{id}=norm_dir_resp_loc_concat{id};
    norm_dir_resp_loc_concat_rect{id}(norm_dir_resp_loc_concat_rect{id}<0)=0;
end

problems=[];
for id = 1:nd
    osi_stat{id}=nan(nKeep_total,1);
    osi_loc{id}=nan(nKeep_total,1);

   prefs = mean(norm_dir_resp_stat_concat_rect{id}(:,1,:),3); %find the peak value, averaging over contrast
   orth= mean(norm_dir_resp_stat_concat_rect{id}(:,3,:),3); %find the +90-degree value, averaging over contrast
   problems_stat=find(prefs < orth);
   prefs(prefs < orth)=nan;
   osi_stat{id}(:,1)=(prefs-orth)./(prefs+orth);


   prefs_loc = mean(norm_dir_resp_loc_concat_rect{id}(:,1,:),3); %find the peak value, averaging over contrast
   orth_loc= mean(norm_dir_resp_loc_concat_rect{id}(:,3,:),3); %find the +90-degree value, averaging over contrast
   problesm_loc=find(prefs_loc < orth_loc);
   prefs_loc(prefs_loc < orth_loc)=nan;
   osi_loc{id}(:,1)=(prefs_loc-orth_loc)./(prefs_loc+orth_loc);
end



%% delta OSI 
%identify cells with OSI > 0.5 on the baseline day
OSI_stat_include = find(osi_stat{pre}>0.5);
OSI_loc_include = find(osi_loc{pre}>0.5);

%% plot distribution of OSI pre and post


figure
subplot(2,2,1)
h1=cdfplot(mean(osi_stat{pre}(green_ind_concat,:),2,'omitmissing'));
hold on
h2 = cdfplot(mean(osi_stat{post}(green_ind_concat,:),2,'omitmissing'));
set(h1, 'Color', 'k');
set(h2, 'Color', 'b');
title('Pyr stationary')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('OSI')
ylabel('Cumulative distribution')

subplot(2,2,2)
h1=cdfplot(mean(osi_stat{pre}(red_ind_concat,:),2,'omitmissing'));
hold on
h2=cdfplot(mean(osi_stat{post}(red_ind_concat,:),2,'omitmissing'));
set(h1, 'Color', 'k');
set(h2, 'Color', 'b');
title('SST stationary')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('OSI')
ylabel('Cumulative distribution')

subplot(2,2,3)
h1=cdfplot(mean(osi_loc{pre}(green_ind_concat,:),2,'omitmissing'));
hold on
h2=cdfplot(mean(osi_loc{post}(green_ind_concat,:),2,'omitmissing'));
set(h1, 'Color', 'k');
set(h2, 'Color', 'b');
title('Pyr running')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('OSI')
ylabel('Cumulative distribution')

subplot(2,2,4)
h1=cdfplot(mean(osi_loc{pre}(red_ind_concat,:),2,'omitmissing'));
hold on
h2=cdfplot(mean(osi_loc{post}(red_ind_concat,:),2,'omitmissing'));
set(h1, 'Color', 'k');
set(h2, 'Color', 'b');
title('SST running')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('OSI')
ylabel('Cumulative distribution')

sgtitle('Distribution of OSI')
print(fullfile(fnout,['OSI_distribution.pdf']),'-dpdf','-bestfit')

%% plot contrast response

conResp_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_stat = cell(1,nd); %same for red
conResp_green_se_stat = cell(1,nd); %this will be the se across all green cells
conResp_red_se_stat = cell(1,nd); %same for red



for id = 1:nd
   
        
    conResp_green_avrg_stat{id}=mean(pref_responses_stat_concat{id}(green_all,:),1,'omitnan');
    green_std=std(pref_responses_stat_concat{id}(green_all,:),1,'omitnan');
    conResp_green_se_stat{id}=green_std/sqrt(length(green_all));
    
    conResp_red_avrg_stat{id}=mean(pref_responses_stat_concat{id}(red_all,:),1,'omitnan');
    red_std=std(pref_responses_stat_concat{id}(red_all,:),1,'omitnan');
    conResp_red_se_stat{id}=red_std/sqrt(length(red_all));
    
    clear green_std red_std
 
end


figure
subplot(1,2,1) %for the first day
errorbar(cons,conResp_green_avrg_stat{pre},conResp_green_se_stat{pre},'k');
hold on
errorbar(cons,conResp_green_avrg_stat{post},conResp_green_se_stat{post},'b');
title(['Pyr',' n = ', num2str(length(green_all))])
ylabel('dF/F, pref ori') 
xlabel('contrast') 
set(gca, 'TickDir', 'out')
xlim([0 1])
ylim([0 .3])
xticks([.25 .5 1])
box off

subplot(1,2,2) %for the second day
errorbar(cons,conResp_red_avrg_stat{pre},conResp_red_se_stat{pre},'k');
hold on
errorbar(cons,conResp_red_avrg_stat{post},conResp_red_se_stat{post},'b');
title(['SST',' n = ', num2str(length(red_all))])
xlabel('contrast') 
xlim([0 1])
ylim([0 .3])
xticks([.25 .5 1])
set(gca, 'TickDir', 'out')
box off

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])



sgtitle(['stationary' ])

print(fullfile(fnout,['contrast_response.pdf']),'-dpdf');



conResp_green_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_loc = cell(1,nd); %same for red
conResp_green_se_loc = cell(1,nd); %this will be the se across all green cells
conResp_red_se_loc = cell(1,nd); %same for red



for id = 1:nd
   
        
    conResp_green_avrg_loc{id}=nanmean(pref_responses_loc_concat{id}(green_all ,:),1);
    green_std=nanstd(pref_responses_loc_concat{id}(green_all,:),1);
    conResp_green_se_loc{id}=green_std/sqrt(length(green_all));
    
    conResp_red_avrg_loc{id}=nanmean(pref_responses_loc_concat{id}(red_all,:),1);
    red_std=nanstd(pref_responses_loc_concat{id}(red_all,:),1);
    conResp_red_se_loc{id}=red_std/sqrt(length(red_all));
    
    clear green_std red_std
 
end


figure
subplot(1,2,1) %for the first day
errorbar(cons,conResp_green_avrg_loc{pre},conResp_green_se_loc{pre},'k');
hold on
errorbar(cons,conResp_green_avrg_loc{post},conResp_green_se_loc{post},'b');
title(['Pyr',' n = ', num2str(length(green_all))])
ylabel('dF/F, pref ori') 
xlabel('contrast') 
xlim([0 1])
ylim([0 .4])
xticks([.25 .5 1])
set(gca, 'TickDir', 'out')
box off

subplot(1,2,2) 
errorbar(cons,conResp_red_avrg_loc{pre},conResp_red_se_loc{pre},'k');
hold on
errorbar(cons,conResp_red_avrg_loc{post},conResp_red_se_loc{post},'b');
title(['SST',' n = ', num2str(length(red_all))])
 ylim([0 .4])
xlim([0 1])
xticks([.25 .5 1])
xlabel('contrast') 
set(gca, 'TickDir', 'out')
box off

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])



sgtitle(['Running' ])

print(fullfile(fnout,['loc_contrast_resposnse.pdf']),'-dpdf');
%% contrast response pupil

conResp_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_stat = cell(1,nd); %same for red
conResp_green_se_stat = cell(1,nd); %this will be the se across all green cells
conResp_red_se_stat = cell(1,nd); %same for red

temp_green = intersect(green_ind_concat,have_allPupil);
temp_red = intersect(red_ind_concat,have_allPupil);

for id = 1:nd
    conResp_green_avrg_stat{id}=mean(pref_responses_stat_lowPupil_concat{id}(temp_green,:),1,'omitnan');
    green_std=std(pref_responses_stat_lowPupil_concat{id}(temp_green,:),1,'omitnan');
    conResp_green_se_stat{id}=green_std/sqrt(length(temp_green));
    
    conResp_red_avrg_stat{id}=mean(pref_responses_stat_lowPupil_concat{id}(temp_red,:),1,'omitnan');
    red_std=std(pref_responses_stat_lowPupil_concat{id}(temp_red,:),1,'omitnan');
    conResp_red_se_stat{id}=red_std/sqrt(length(temp_red));
    
    clear green_std red_std
 
end


figure
subplot(2,2,1) %for the first day
errorbar(cons,conResp_green_avrg_stat{pre},conResp_green_se_stat{pre},'k');
hold on
errorbar(cons,conResp_green_avrg_stat{post},conResp_green_se_stat{post},'b');
title(['Small pupil Pyr',' n = ', num2str(length(temp_green))])
ylabel('dF/F, pref ori') 
%xlabel('contrast') 
set(gca, 'TickDir', 'out')
xlim([0 1])
ylim([0 .3])
xticks([.25 .5 1])
box off

subplot(2,2,2) %for the second day
errorbar(cons,conResp_red_avrg_stat{pre},conResp_red_se_stat{pre},'k');
hold on
errorbar(cons,conResp_red_avrg_stat{post},conResp_red_se_stat{post},'b');
title(['Small pupil SST',' n = ', num2str(length(temp_red))])
%xlabel('contrast') 
xlim([0 1])
ylim([0 .3])
xticks([.25 .5 1])
set(gca, 'TickDir', 'out')
box off

conResp_green_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_loc = cell(1,nd); %same for red
conResp_green_se_loc = cell(1,nd); %this will be the se across all green cells
conResp_red_se_loc = cell(1,nd); %same for red



for id = 1:nd
   
        
    conResp_green_avrg_loc{id}=nanmean(pref_responses_stat_hiPupil_concat{id}(green_all ,:),1);
    green_std=nanstd(pref_responses_stat_hiPupil_concat{id}(green_all,:),1);
    conResp_green_se_loc{id}=green_std/sqrt(length(green_all));
    
    conResp_red_avrg_loc{id}=nanmean(pref_responses_stat_hiPupil_concat{id}(red_all,:),1);
    red_std=nanstd(pref_responses_stat_hiPupil_concat{id}(red_all,:),1);
    conResp_red_se_loc{id}=red_std/sqrt(length(red_all));
    
    clear green_std red_std
 
end


subplot(2,2,3) %for the first day
errorbar(cons,conResp_green_avrg_loc{pre},conResp_green_se_loc{pre},'k');
hold on
errorbar(cons,conResp_green_avrg_loc{post},conResp_green_se_loc{post},'b');
title(['Large pupil Pyr',' n = ', num2str(length(green_all))])
ylabel('dF/F, pref ori') 
xlabel('contrast') 
xlim([0 1])
ylim([0 .4])
xticks([.25 .5 1])
set(gca, 'TickDir', 'out')
box off

subplot(2,2,4) 
errorbar(cons,conResp_red_avrg_loc{pre},conResp_red_se_loc{pre},'k');
hold on
errorbar(cons,conResp_red_avrg_loc{post},conResp_red_se_loc{post},'b');
title(['Large pupil SST',' n = ', num2str(length(red_all))])
 ylim([0 .4])
xlim([0 1])
xticks([.25 .5 1])
xlabel('contrast') 
set(gca, 'TickDir', 'out')
box off

x0=5;
y0=5;
width=5;
height=6;
set(gcf,'units','inches','position',[x0,y0,width,height])



print(fullfile(fnout,['pupil_contrast_resposnse.pdf']),'-dpdf');

%% response by condition

a=mean(pref_responses_stat_concat{pre}(green_all,:), "omitnan");
b=mean(pref_responses_loc_concat{pre}(green_all,:), "omitnan");

c=mean(pref_responses_stat_concat{pre}(red_all,:), "omitnan");
d=mean(pref_responses_loc_concat{pre}(red_all,:), "omitnan");

e=mean(pref_responses_stat_concat{post}(green_all,:), "omitnan");
f=mean(pref_responses_loc_concat{post}(green_all,:), "omitnan");

g=mean(pref_responses_stat_concat{post}(red_all,:), "omitnan");
h=mean(pref_responses_loc_concat{post}(red_all,:), "omitnan");

responseByCond = horzcat([a';b'],[c';d'],[e';f'],[g';h']);
clear a b c d e f g h

figure;
scatter(responseByCond(:,1),responseByCond(:,2),'k')
hold on
scatter(responseByCond(:,3),responseByCond(:,4),'b')
hold on

h2 = lsline;
set(h2(2),'color','k','LineWidth',2)
set(h2(1),'color','b','LineWidth',2)
hold on
scatter(responseByCond(4:6,1),responseByCond(4:6,2),'MarkerEdgeColor', 'k','MarkerFaceColor', 'k')
hold on
scatter(responseByCond(4:6,3),responseByCond(4:6,4),'MarkerEdgeColor','b','MarkerFaceColor','b')
set(gca,'TickDir','out')
xlabel('Mean Pyr dF/F')
ylabel('Mean SOM dF/F')
title('Mean response to different conditions')




for iMouse = 1:nSess 
    x = responseByCondProps_concat(5,1,iMouse):.005:responseByCondProps_concat(6,1,iMouse);
    m  = responseByCondProps_concat(1,1,iMouse);
    y = (m*x)+responseByCondProps_concat(2,1,iMouse);
    plot(x,y,'Color',"#808080") 


    x = responseByCondProps_concat(5,2,iMouse):.005:responseByCondProps_concat(6,2,iMouse);
    m  = responseByCondProps_concat(1,2,iMouse);
    y = (m*x)+responseByCondProps_concat(2,2,iMouse);
    plot(x,y,'Color',"#0080ff")  

    % x_population=min(responseByCondProps_concat(5,1,:)):.005:max(responseByCondProps_concat(6,1,:));
    % m_mean=mean(responseByCondProps_concat(1,1,:),'omitmissing');
    % b_mean=mean(responseByCondProps_concat(2,1,:),'omitmissing');
    % y_mean=(m_mean*x_population)+b_mean;
    % plot(x_population,y_mean,'k','LineWidth',2)
    % 
    % x_population=min(responseByCondProps_concat(5,2,:)):.005:max(responseByCondProps_concat(6,2,:));
    % m_mean=mean(responseByCondProps_concat(1,2,:),'omitmissing');
    % b_mean=mean(responseByCondProps_concat(2,2,:),'omitmissing');
    % y_mean=(m_mean*x_population)+b_mean;
    % plot(x_population,y_mean,'b','LineWidth',2)

end

ylim([-.01 .4])
xlim([-.01 .4])

print(fullfile(fnout, ['responseByCondition.pdf']),'-dpdf','-bestfit')



%% response by condidiont properties
% property X pre/post (1=pre, 2=post) X mouse
%properties: range = 5-6, slope = 1, intercept = 2
figure
subplot(1,2,1)
slopes=squeeze(responseByCondProps_concat(1,:,:));
plot(slopes, "-o",'Color',"#7b7b7b")
xlim([.5 2.5])
box off
set(gca,'TickDir','out')
xticks([1 2])
xticklabels({'Pre-DART','Post-DART'})
ylabel('Slope off Pyr vs SST activity across conditions')
hold on
plot(nanmean(slopes,2),"-*K",'LineWidth',2)
title("Slope change")

subplot(1,2,2)
intercepts=squeeze(responseByCondProps_concat(2,:,:));
plot(intercepts, "-o",'Color',"#7b7b7b")
xlim([.5 2.5])
box off
set(gca,'TickDir','out')
xticks([1 2])
xticklabels({'Pre-DART','Post-DART'})
ylabel('Intercept off Pyr vs SST activity across conditions')
hold on
plot(nanmean(intercepts,2),"-*K",'LineWidth',2)
title("Intercept change")

print(fullfile(fnout, ['SlopeAndIntChange.pdf']),'-dpdf','-bestfit')
%% calculate fract chage
raw_diff_stat = nan(nCon,nKeep_total);
raw_diff_loc = nan(nCon,nKeep_total);
fract_diff_stat = nan(nCon,nKeep_total);
fract_diff_loc = nan(nCon,nKeep_total);
for iCon = 1:nCon
    raw_diff_stat(iCon,:)=((pref_responses_stat_concat{post}(:,iCon))-(pref_responses_stat_concat{pre}(:,iCon)))';
    raw_diff_loc(iCon,:)=((pref_responses_loc_concat{post}(:,iCon))-(pref_responses_loc_concat{pre}(:,iCon)));
    
    fract_diff_stat(iCon,:)=((pref_responses_stat_concat{post}(:,iCon))-(pref_responses_stat_concat{pre}(:,iCon)))./(pref_responses_stat_concat{pre}(:,iCon));
    fract_diff_loc(iCon,:)=((pref_responses_loc_concat{post}(:,iCon))-(pref_responses_loc_concat{pre}(:,iCon)))./(pref_responses_loc_concat{pre}(:,iCon));

end


 %% bar chart of fraction red cells significantly suppressed or faciltiated
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
norm_diff(find(norm_diff == -Inf))=nan;
norm_diff(find(norm_diff == Inf))=nan;

%find how many cells are suppressed ( normalized diff < -1) or facilitated
%(normalized diff >1) at each contrast and behavioral state

suppressed = logical(norm_diff < -1);
facilitated = logical(norm_diff > 1);


%pull out red cells, get fractions that are suppressed vs. favilitated, for
%the red cells that have running trials at all contrasts
fractSupp_red = nan(2,3);
fractFacil_red = nan(2,3);
for iCon = 1:nCon
    fractSupp_red(:,iCon) = sum(((suppressed(:,iCon,red_all))),3)/length(red_all);
    fractFacil_red(:,iCon) = sum(((facilitated(:,iCon,red_all))),3)/length(red_all);
end


%averaging over contrast
supp_red_mean = mean(fractSupp_red,2)
facil_red_mean = mean(fractFacil_red,2)



figure;
for iCon = 1:nCon
    subplot(1,3,iCon)
    bar([1,2],[fractSupp_red(1,iCon) fractFacil_red(1,iCon);fractSupp_red(2,iCon) fractFacil_red(2,iCon)],'stacked')
    xticklabels({'Stationary','Running'})
    ylabel(["Fraction HTP+ cells"]) 
    set(gca,'TickDir','out')
    box off
    %ylim([0 .7])

end
sgtitle('fraction SST suppressed/facilitated by > 1std')
x0=5;
y0=5;
width=9;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,[ 'factionSuppFacilBar.pdf']),'-dpdf','-bestfit')

norm_diff_red_allKeep = nanmedian(norm_diff(1,:,red_ind_concat),3)



%% plot distributions of normalized difference 
figure;
subplot(1,2,1)
boxchart((squeeze(mean(norm_diff(:,:,red_all),2,"omitnan")))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75]);hold on;
scatter([1 2],(squeeze(mean(norm_diff(:,:,red_all),2,"omitnan")))',"red",'jitter', 'on', 'jitterAmount',.1)
%scatter([1 2],(squeeze(mean(norm_diff(:,:,green_all),2,"omitnan")))',"black",'jitter', 'on', 'jitterAmount',.1)
%plot(squeeze(mean(norm_diff(:,:,red_all),2,"omitnan")),'--k')
ylim([-6 10])
xticklabels({'Stationary','Running'})
ylabel('(Post-Pre) / std dev Pre')
title('SST')
hold off
set(gca,'TickDir','out')
box off

subplot(1,2,2)
boxchart((squeeze(mean(norm_diff(:,:,green_all),2,"omitnan")))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75]);hold on;
scatter([1 2],(squeeze(mean(norm_diff(:,:,green_all),2,"omitnan")))',"black",'jitter', 'on', 'jitterAmount',.1)
%plot(squeeze(mean(norm_diff(:,:,red_all),2,"omitnan")),'--k')
ylim([-6 10])
xticklabels({'Stationary','Running'})
ylabel('(Post-Pre) / std dev Pre')
title('Pyr')
hold off
set(gca,'TickDir','out')
box off

print(fullfile(fnout,[ 'normDiff_by_behState.pdf']),'-dpdf')

figure;
subplot(1,2,1)
boxchart(squeeze(norm_diff(1,:,red_all))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75]);
hold on
scatter([1, 2, 3],squeeze(norm_diff(1,:,red_all))',"red",'jitter', 'on', 'jitterAmount',.1)
%scatter([1, 2, 3],squeeze(norm_diff(1,:,green_ind_concat))',"black",'jitter', 'on', 'jitterAmount',.1)
%plot(squeeze(norm_diff(1,[1,3],red_ind_concat)))
xticklabels({'25%','50%','100%'})
ylim([-6 10])
xlabel('Contrast')
ylabel('(Post-Pre) / std dev Pre')
title('SST')
hold off
set(gca,'TickDir','out')
box off

subplot(1,2,2)
boxchart(squeeze(norm_diff(1,:,green_all))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75]);
hold on
scatter([1, 2, 3],squeeze(norm_diff(1,:,green_all))',"black",'jitter', 'on', 'jitterAmount',.1)
%plot(squeeze(norm_diff(1,[1,3],red_ind_concat)))
xticklabels({'25%','50%','100%'})
ylim([-6 10])
xlabel('Contrast')
ylabel('(Post-Pre) / std dev Pre')
title('Pyr')
hold off
set(gca,'TickDir','out')
box off

print(fullfile(fnout,[ 'normDiff_by_contrast.pdf']),'-dpdf')



%% time to peak - this is combining running, stationary, and all stimulus conditions

tHalfMax = cell(1,nd);

for id = 1:nd
    tHalfMaxTemp=nan(nKeep_total,1);

        for iCell=1:nKeep_total
            %pull the data for a given cell at a given contrast (pref ori)
            tempData=tc_trial_avrg_stat_concat{id}(stimStart:stimStart+nOn,iCell);
            if resp_keep_concat{id}(iCell)
                if mean(tempData >0)
                smoothData=smoothdata(tempData,'movmean',5) ;
                halfMax = max(smoothData(3:length(smoothData)))/2;
                tHalfMaxCell =double(min(find(smoothData>halfMax)))/double(frame_rate);
                if length(tHalfMaxCell)>0
                    tHalfMaxTemp(iCell)=tHalfMaxCell;
                end
                end
            end
 
       end
    tHalfMax{id}=tHalfMaxTemp;
end
clear tHalfMaxCell tHalfMaxTemp tempData smoothData halfMax
% scatter for tMax



figure; movegui('center') 
subplot(1,2,1)
scatter((tHalfMax{pre}(green_ind_concat)),(tHalfMax{post}(green_ind_concat)),10,'MarkerEdgeColor',[.4 .4 .4],'jitter', 'on', 'jitterAmount',.1)
hold on
mean_pre_stat = mean(tHalfMax{pre}(green_ind_concat),"omitnan");
mean_post_stat = mean(tHalfMax{post}(green_ind_concat),"omitnan");
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= nanstd(tHalfMax{pre}(green_ind_concat)) / sqrt( length(green_ind_concat));
stderror_post= nanstd(tHalfMax{post}(green_ind_concat)) / sqrt( length(green_ind_concat));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);
ylabel('post-DART half-max(s)')
xlabel('pre-DART half-max(s)')
ylim([0 2.15])
xlim([0 2.15])
hline=refline(1);
hline.Color = 'k';
hline.LineStyle = ':';
title('Pyr ')
axis square
hold off

subplot(1,2,2)
scatter((tHalfMax{pre}(red_ind_concat)),(tHalfMax{post}(red_ind_concat)),10,'MarkerEdgeColor',[.4 .4 .4],'jitter', 'on', 'jitterAmount',.1)
hold on
mean_pre_stat = nanmean(tHalfMax{pre}(red_ind_concat));
mean_post_stat = nanmean(tHalfMax{post}(red_ind_concat));
scatter(mean_pre_stat,mean_post_stat,20,'r','filled');
stderror_pre= nanstd(tHalfMax{pre}(red_ind_concat)) / sqrt( length(red_ind_concat));
stderror_post= nanstd(tHalfMax{post}(red_ind_concat)) / sqrt( length(red_ind_concat));
line([(mean_pre_stat-stderror_pre),(mean_pre_stat+stderror_pre)],[mean_post_stat, mean_post_stat],'Color','black','LineWidth',2);
line([mean_pre_stat, mean_pre_stat],[(mean_post_stat-stderror_post),(mean_post_stat+stderror_post)],'Color','blue','LineWidth',2);

%ylabel('post-DART half-max(s)')
xlabel('pre-DART half-max(s)')
ylim([0 2.15])
xlim([0 2.15])
hline=refline(1)
hline.Color = 'k';
hline.LineStyle = ':';
title('SST')
axis square
hold off

clear txt1 NrPoints
set(gca,'TickDir','out')

sgtitle('time to half max (s), stat')
print(fullfile(fnout,[ 'stat_tHalfMax_crossDay.pdf']),'-dpdf','-bestfit')

[h1, p1]= ttest(tHalfMax{pre}(green_all),tHalfMax{post}(green_all));
[h2,p2]= ttest(tHalfMax{pre}(red_all),tHalfMax{post}(red_all));

%correct for two tests
p1*2
p2*2
clear h1 p1 h2 p2

%% locomotion modulation index


for iCon = 1:nCon

figure; movegui('center') 
subplot(1,2,1)
scatter((LMI_concat{pre}(green_ind_concat,iCon)),(LMI_concat{post}(green_ind_concat,iCon)),10,'MarkerEdgeColor',[.4 .4 .4],'jitter', 'on', 'jitterAmount',.01)
ylabel('post-DART LMI')
xlabel('pre-DART  LMI')
ylim([-1 1])
xlim([-1 1])
% hline(0)
% vline(0)
refline(1)
hold on
scatter(nanmean(LMI_concat{pre}(green_ind_concat,iCon)),nanmean(LMI_concat{post}(green_ind_concat,iCon)),15,'r*')
axis square
title('Pyr ')



subplot(1,2,2)
scatter((LMI_concat{pre}(red_ind_concat,iCon)),(LMI_concat{post}(red_ind_concat,iCon)),10,'MarkerEdgeColor',[.4 .4 .4],'jitter', 'on', 'jitterAmount',.01)
hold on
ylabel('post-DART LMI')
xlabel('pre-DART  LMI')
ylim([-1 1])
xlim([-1 1])
% hline(0)
% vline(0)
refline(1)
hold on
axis square
scatter(nanmean(LMI_concat{pre}(red_ind_concat,iCon)),nanmean(LMI_concat{post}(red_ind_concat,iCon)),15,'r*')
title('SST')


clear txt1 NrPoints
sgtitle([num2str(cons(iCon))])
set(gca,'TickDir','out')
print(fullfile(fnout,[num2str(cons(iCon)) '_LMI.pdf']),'-dpdf');

end
%% Identify high and low correlation cells

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
RbyExp
%% stat high and low R
hi_avrg_stat = cell(1,nd);
low_avrg_stat = cell(1,nd); 
hi_se_stat = cell(1,nd); 
low_se_stat = cell(1,nd);

for id = 1:nd

   for iCon=1:nCon
        tempInds_high = intersect(redHigh, haveRunning_red{iCon});
        hi_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,tempInds_high,iCon),2);
        high_std=nanstd(tc_trial_avrg_stat_concat{id}(:,tempInds_high,iCon),[],2);
        hi_se_stat{id}(:,iCon)=high_std/sqrt(length(tempInds_high));
        
        tempInds_low = intersect(redLow, haveRunning_red{iCon});
        low_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,tempInds_low,iCon),2);
        low_std=nanstd(tc_trial_avrg_stat_concat{id}(:,tempInds_low,iCon),[],2);
        low_se_stat{id}(:,iCon)=low_std/sqrt(length(tempInds_low));
        
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
num=length(intersect(redLow, haveRunning_red{iCon}));
title([' Weakly correlated',' n = ', num2str(num)])
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
num=length(intersect(redHigh, haveRunning_red{iCon}));
title([' Strongly correlated',' n = ', num2str(num)])

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



%% loc high and low R
hi_avrg_loc = cell(1,nd);
low_avrg_loc = cell(1,nd); 
hi_se_loc = cell(1,nd); 
low_se_loc = cell(1,nd);

for id = 1:nd

   for iCon=1:nCon
        tempInds_high = intersect(redHigh, haveRunning_red{iCon});
        hi_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,tempInds_high,iCon),2);
        high_std=nanstd(tc_trial_avrg_loc_concat{id}(:,tempInds_high,iCon),[],2);
        hi_se_loc{id}(:,iCon)=high_std/sqrt(length(tempInds_high));
        
        tempInds_low = intersect(redLow, haveRunning_red{iCon});
        low_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,tempInds_low,iCon),2);
        low_std=nanstd(tc_trial_avrg_loc_concat{id}(:,tempInds_low,iCon),[],2);
        low_se_loc{id}(:,iCon)=low_std/sqrt(length(tempInds_low));
        
        clear low_std high_std
    end
end
z=double(nOn)/double(frame_rate);

%creat a time axis in seconds
t=1:(size(hi_avrg_loc{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure



subplot(1,2,1) 
shadedErrorBar(t,low_avrg_loc{pre}(:,iCon),low_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,low_avrg_loc{post}(:,iCon),low_se_loc{post}(:,iCon),'b');
ylim([-.02 .5]);
hold on
% line([0,.2],[-.01,-.01],'Color','black','LineWidth',2);
% hold on
line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
num=length(intersect(redLow, haveRunning_red{iCon}));
title([' Weakly correlated',' n = ', num2str(num)])
set(gca,'XColor', 'none','YColor','none')


subplot(1,2,2) 
ylim([-.02 .5]);
hold on
shadedErrorBar(t,hi_avrg_loc{pre}(:,iCon),hi_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,hi_avrg_loc{post}(:,iCon),hi_se_loc{post}(:,iCon),'b');
hold on
% line([0,.2],[-.01,-.01],'Color','black','LineWidth',2);
% hold on
line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
num=length(intersect(redHigh, haveRunning_red{iCon}));
title([' Strongly correlated',' n = ', num2str(num)])

xlabel('s') 
set(gca,'XColor', 'none','YColor','none')

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['Running, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) 'loc_R_timecourses.pdf']),'-dpdf');
clear txt1 highRed lowRed
end 





%% making dfof output table
mouseID=[];
for imouse =1:nSess
   ID_string=mouseNames(imouse);
   thisID = repelem(ID_string,nKeep_concat(imouse));
   mouseID=[mouseID,thisID];
end
mouseID=mouseID';

dfof_stat_pre = reshape(pref_responses_stat_concat{pre},[3*nKeep_total,1]);
dfof_stat_post = reshape(pref_responses_stat_concat{post},[3*nKeep_total,1]);
dfof_loc_pre = reshape(pref_responses_loc_concat{pre},[3*nKeep_total,1]);
dfof_loc_post = reshape(pref_responses_loc_concat{post},[3*nKeep_total,1]);

dfof_col = vertcat(dfof_stat_pre,dfof_stat_post,dfof_loc_pre,dfof_loc_post);

cellID=1:nKeep_total;
cellID_col=repmat(cellID',12,1);
cellID_col = cellID_col+622; %to account fro DRT cells

mouseIDcol = repmat(mouseID,12,1);

day1 = repelem(["pre" "post"],1,(3*nKeep_total))';
day=repmat(day1,2,1);

drug="PEG"; %%need to get drug coded correctly
drugCol=repmat(drug,size(mouseIDcol));

contrast1 = repelem(cons,nKeep_total)';
contrast=repmat(contrast1,4,1);

cell_type_col=repmat(red_concat,1,12)';

behStateCol = repelem(["stat" "loc"],1,(6*nKeep_total))';





dfof_output = table(mouseIDcol,cellID_col,cell_type_col,contrast,behStateCol,drugCol,day,dfof_col, ...
    'VariableNames',{'mouseID' 'cellID' 'cellType' 'contrast' 'behState' 'drug' 'day' 'dfof'});

writetable(dfof_output,fullfile(fnout,'dfof_DART_output.csv'))



