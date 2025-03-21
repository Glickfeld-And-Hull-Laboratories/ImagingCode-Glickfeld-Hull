% %% read in and concatenate data

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
    base = fullfile(isilonName, '/All_Staff/home/ACh\Analysis\2p_analysis\SST_YM90K');
    
else
    isilonName = 'Z:';
    base = fullfile(isilonName, '\home\ACh\Analysis\2p_analysis\SST_YM90K');
      
end

d=string(datetime('today'));
fnout= fullfile(base,ds_name,d);

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
resp_keep_concat=cell(1,nd);
resp_max_keep_concat=cell(1,nd);
pref_responses_loc_concat=cell(1,nd);
pref_responses_stat_concat=cell(1,nd);
pref_peak_stat_concat=cell(1,nd);
pref_peak_loc_concat=cell(1,nd);
pref_responses_stat_largePupil_concat=cell(1,nd);
pref_responses_stat_smallPupil_concat=cell(1,nd);
RIx_concat=cell(1,nd);
dirs_concat=[];
cons_concat=[];
green_concat=[];
red_concat=[];
dfof_max_diff_concat=[];
nKeep_concat=[];
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
noiseCorr_concat = cell(4,nd);
noiseCorr_OG_concat = cell(1,nd);
noiseCorrContrast_concat = cell(4,nCon,nd); 
sigCorr_concat = cell(1,nd);
pref_allTrials_stat_concat =cell(nCon,nd);
pref_allTrials_loc_concat =cell(nCon,nd);
pref_allTrials_largePupil_concat =cell(nCon,nd);
pref_allTrials_smallPupil_concat =cell(nCon,nd);
dataTableConat=[];
drug=cell(1,nSess);
pupilMeans_concat=nan(nd,3,nSess);
motorByPupil_concat=nan(nd,2,nSess);
pupilCounts_concat=nan(nd,2,nSess);
nonPref_trial_avrg_stat_concat=cell(1,nd);
nonPref_trial_avrg_loc_concat=cell(1,nd);

responseByCondProps_concat=nan(6,2,nSess);

cellID_adjustment=0;
for iSess = 1:nSess
    mouse = ds_YM90K_DART(iSess).mouse;
    mice=[mice;mouse];
    thisDrug = ds_YM90K_DART(iSess).drug;

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


    nKeep = size(tc_trial_avrg_stat{post},2);

    pupilMeans_concat(:,:,iSess)=pupilMeans;
    motorByPupil_concat(:,:,iSess)=motorByPupil;
%   pupilCounts_concat(:,:,iSess)=pupilCounts;


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
        tc_trial_avrg_stat_largePupil_concat{id} = cat(2,tc_trial_avrg_stat_largePupil_concat{id},tc_trial_avrg_stat_largePupil{id}(:,:,sharedCon));
        tc_trial_avrg_stat_smallPupil_concat{id} = cat(2,tc_trial_avrg_stat_smallPupil_concat{id},tc_trial_avrg_stat_smallPupil{id}(:,:,sharedCon));
        tc_trial_avrg_loc_concat{id} =cat(2,tc_trial_avrg_loc_concat{id},tc_trial_avrg_loc{id}(:,:,sharedCon));
        nonPref_trial_avrg_stat_concat{id} =cat(2,nonPref_trial_avrg_stat_concat{id},nonPref_trial_avrg_stat{id}(:,:,sharedCon));
        nonPref_trial_avrg_loc_concat{id} =cat(2,nonPref_trial_avrg_loc_concat{id},nonPref_trial_avrg_loc{id}(:,:,sharedCon));
        resp_keep_concat{id}=cat(1,resp_keep_concat{id},resp_keep{id});
        resp_max_keep_concat{id}=cat(1,resp_max_keep_concat{id},resp_max_keep{id}(:,sharedCon));
        pref_responses_loc_concat{id}=cat(1,pref_responses_loc_concat{id},pref_responses_loc{id}(:,sharedCon));
        pref_responses_stat_concat{id}=cat(1,pref_responses_stat_concat{id},pref_responses_stat{id}(:,sharedCon));
        pref_peak_stat_concat{id}=cat(1,pref_peak_stat_concat{id},pref_peak_stat{id}(:,sharedCon));
        pref_peak_loc_concat{id}=cat(1,pref_peak_loc_concat{id},pref_peak_loc{id}(:,sharedCon));
        pref_responses_stat_largePupil_concat{id}=cat(1,pref_responses_stat_largePupil_concat{id},pref_responses_stat_largePupil{id}(:,sharedCon));
        pref_responses_stat_smallPupil_concat{id}=cat(1,pref_responses_stat_smallPupil_concat{id},pref_responses_stat_smallPupil{id}(:,sharedCon));
        RIx_concat{id}=cat(1,RIx_concat{id},sum(RIx{id}));
        wheel_corr_concat{id}=cat(2,wheel_corr_concat{id},wheel_corr{id});
        meanF=mean(fullTC_keep{id},1);
        meanF_concat{id}=cat(2,meanF_concat{id}, meanF);
        norm_dir_resp_stat_concat{id}=cat(1,norm_dir_resp_stat_concat{id},norm_dir_resp_stat{id});
        norm_dir_resp_loc_concat{id}=cat(1,norm_dir_resp_loc_concat{id},norm_dir_resp_loc{id});
        pref_dir_concat{id}=cat(2,pref_dir_concat{id},pref_dir_keep{id});
        for i=1:4
            noiseCorr_concat{i,id}=cat(2,noiseCorr_concat{i,id},noiseCorr{i,id});
            for iCon = 1:nCon
                noiseCorrContrast_concat{i,iCon,id}=cat(2,noiseCorrContrast_concat{i,iCon,id},noiseCorrContrast{i,iCon,id});
            end
        end
        noiseCorr_OG_concat{id}=cat(2,noiseCorr_OG_concat{id},noiseCorr_OG{id});
        sigCorr_concat{id}=cat(2,sigCorr_concat{id},sigCorr{id});
        for i = 1:length(sharedCon)
            iCon=sharedCon(i);
            pref_allTrials_stat_concat{i,id}=[pref_allTrials_stat_concat{i,id},pref_allTrials_stat{iCon,id}];
            pref_allTrials_loc_concat{i,id}=[pref_allTrials_loc_concat{i,id},pref_allTrials_loc{iCon,id}];
            pref_allTrials_largePupil_concat{i,id}=[pref_allTrials_largePupil_concat{i,id},pref_allTrials_largePupil{iCon,id}];
            pref_allTrials_smallPupil_concat{i,id}=[pref_allTrials_smallPupil_concat{i,id},pref_allTrials_smallPupil{iCon,id}];
        end
        clear meanF i
    end
    dfof_max_diff_concat=cat(1,dfof_max_diff_concat,dfof_max_diff(:,sharedCon));
   green_fluor_concat=cat(2,green_fluor_concat,green_fluor_keep);
   red_fluor_concat=cat(2,red_fluor_concat,red_fluor_keep);
   if size(pref_responses_stat_concat{id},1) ~= size(pref_responses_stat_smallPupil_concat{id},1)
       break
   end


iSess
end

clear mouse  nKeep  fn_multi cons oris pupilMeans norm_dir_resp_loc norm_dir_resp_stat
clear explanation1 resp_keep sig_diff pref_con_keep pref_ori_keep tOri_match tCon_match data_trial_keep nTrials tc_trial_avrg_keep green_keep_logical red_keep_logical green_ind_keep red_ind_keep
clear LMI RIx locCounts locResp locTCs statResp statTCs wheel_tc ttest_results_stat ttest_results_loc ttest_results_allCon_stat ttest_results_allCon_loc
clear data_con_resp_keep data_ori_resp_keep data_rep_keep dfof_max_diff dfof_max_diff_raw explanation2 resp_max_keep data_resp_keep pref_responses_stat pref_responses_loc
clear tc_trial_avrg_stat tc_trial_avrg_loc fullTC_keep norm_dir_resp sigCorr noiseCorr responseByCondProps
clear red_fluor_all red_fluor_match green_fluor_match green_fluor_match red_fluor_keep green_fluor_keep R_p_values pre_nonPref_stat pre_nonPref_loc
clear pref_responses_stat_hiPupil pref_responses_stat_largePupil pref_responses_stat_lowPupil pref_responses_stat_smallPupil noiseCorr_OG noiseCorrContrast
red_ind_concat = find(red_concat);
green_ind_concat = find(green_concat);

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
haveRunning_both = cell(1,nd);
haveRunning_green= cell(1,nd);
haveRunning_red= cell(1,nd);
for iCon =1:nCon
haveRunning_both{iCon}= find(haveRunning_pre(:,iCon).* haveRunning_post(:,iCon));
haveRunning_green{iCon} = intersect(haveRunning_both{iCon}, green_ind_concat);
haveRunning_red{iCon} = intersect(haveRunning_both{iCon}, red_ind_concat);
end

clear haveRunning_pre haveRunning_post haveRunning_both


%create an array with indices for each mouse
mouseInds=cell(1,nSess);
start=1;
for iMouse = 1:nSess
    mouseInds{iMouse}=start:(start-1)+nKeep_concat(iMouse);
    start = start+nKeep_concat(iMouse)
end
clear start iMouse

green_all = intersect(haveRunning_green{1},haveRunning_green{2});
green_all = intersect(green_all, haveRunning_green{3});

red_all = intersect(haveRunning_red{1},haveRunning_red{2});
red_all = intersect(red_all, haveRunning_red{3});

%find cells that have pupil data on both days
have_LARGE_pre = ~isnan(pref_responses_stat_largePupil_concat{pre});
have_LARGE_post = ~isnan(pref_responses_stat_largePupil_concat{post});

have_SMALL_pre = ~isnan(pref_responses_stat_smallPupil_concat{pre});
have_SMALL_post = ~isnan(pref_responses_stat_smallPupil_concat{post});

have_bothPupil=cell(1,3);
for iCon =1:nCon
    have_HI_both= find(have_LARGE_pre(:,iCon).* have_LARGE_post(:,iCon));
    have_LOW_both=find(have_SMALL_pre(:,iCon).* have_SMALL_post(:,iCon));
    have_bothPupil{iCon}=intersect(have_HI_both,have_LOW_both);
end
 
clear have_LARGE_pre have_LARGE_post have_SMALL_pre have_SMALL_post have_HI_both have_LOW_both

have_allPupil = intersect(have_bothPupil{1},have_bothPupil{2});
have_allPupil = intersect(have_allPupil, have_bothPupil{3});

have_allPupil_green = intersect(have_allPupil, green_ind_concat);
have_allPupil_red = intersect(have_allPupil, red_ind_concat);

%find how many total, running, and pupil cells exist for each mouse
cellCountsRed = nan(nSess,nCon);
mouseNames=[];
for iMouse = 1:nSess
    
    cellCountsRed(iMouse, 1,1)=length(intersect(red_ind_concat,(mouseInds{iMouse})));
    cellCountsRed(iMouse, 2,1)=length(intersect(red_all,(mouseInds{iMouse})));
    cellCountsRed(iMouse, 3,1)=length(intersect(have_allPupil_red,(mouseInds{iMouse})));
        
    mouseNames=[mouseNames, string(mice(iMouse,:))]
end
clear  iMouse


cellCountsGreen = nan(nSess,nCon);
mouseNames=[];
for iMouse = 1:nSess
    
    cellCountsGreen(iMouse, 1,1)=length(intersect(green_ind_concat,(mouseInds{iMouse})));
    cellCountsGreen(iMouse, 2,1)=length(intersect(green_all,(mouseInds{iMouse})));
    cellCountsGreen(iMouse, 3,1)=length(intersect(have_allPupil_green,(mouseInds{iMouse})));
        
    mouseNames=[mouseNames, string(mice(iMouse,:))]
end


cellCountTableRed = table(cellCountsRed, RowNames=mouseNames)
cellCountTableGreen = table(cellCountsGreen, RowNames=mouseNames)
writetable(cellCountTableRed,fullfile(fnout,'cellCounts.csv'),'WriteRowNames',true)
writetable(cellCountTableGreen,fullfile(fnout,'cellCounts_Green.csv'),'WriteRowNames',true)


%find how many haveRunning red cells exist for each mouse
cellCountsRedAll = nan(nSess,1);
for iMouse = 1:nSess
    cellCountsRedAll(iMouse)=length(intersect(red_all,(mouseInds{iMouse})));
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


%calculate norm_diff
norm_diff = nan(2,nCon,nKeep_total);
bsln_std =  nan(2,nCon,nKeep_total);
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
        bsln_std(1,iCon,i)=std_pre_stat; %std for stationary trials on basline day
        bsln_std(2,iCon,i)=std_pre_loc; %std for running trials on basline day
clear mean_pre_stat mean_post_stat std_pre_stat mean_pre_loc mean_post_loc std_pre_loc norn_diff_stat norm_diff_loc norm_diff_stat norm_diff_loc
    end 
end
%remove any infiinty values resulting from divisions by zero, and turn
%those into NANs instead
norm_diff(find(norm_diff == -Inf))=NaN;
norm_diff(find(norm_diff == Inf))=NaN;

bsln_std(find(bsln_std == -Inf))=NaN;
bsln_std(find(bsln_std == Inf))=NaN;

% cells with high correlation in the baseline day
highRInds = find(noiseCorr_OG_concat{pre}(1,:)>0.5);
lowRInds = find(noiseCorr_OG_concat{pre}(1,:)<=0.5);

redHigh=intersect(highRInds, red_ind_concat);
redLow=intersect(lowRInds, red_ind_concat);

greenHigh=intersect(highRInds, green_ind_concat);
greenLow=intersect(lowRInds, green_ind_concat);

% finding how the high and low R cells are distributed over epxeriments

RbyExp = zeros(2,nSess);
for iSess = 1:nSess
    mouseIndsTemp = mouseInds{iSess};
    RbyExp(1,iSess) = length(intersect(mouseIndsTemp,redHigh));
    RbyExp(2,iSess) = length(intersect(mouseIndsTemp,redLow));
end
array2table(RbyExp,RowNames={'high R'  'low R'})
%% 
fractHigh = RbyExp(1,:)./sum(RbyExp,1);
mean(fractHigh)
std(fractHigh)
%% Figure 3C - timecourses

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
print(fullfile(fnout,'Fig_3C.pdf'),'-dpdf');

%% Fig 3D - contrast response for SST and Pyr cells

conResp_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_stat = cell(1,nd); %same for red
conResp_green_se_stat = cell(1,nd); %this will be the se across all green cells
conResp_red_se_stat = cell(1,nd); %same for red

for id = 1:nd
   
        
    conResp_green_avrg_stat{id}=mean(pref_responses_stat_concat{id}(green_ind_concat,:),1,'omitnan');
    green_std=std(pref_responses_stat_concat{id}(green_ind_concat,:),1,'omitnan');
    conResp_green_se_stat{id}=green_std/sqrt(length(green_ind_concat));
    
    conResp_red_avrg_stat{id}=mean(pref_responses_stat_concat{id}(red_ind_concat,:),1,'omitnan');
    red_std=std(pref_responses_stat_concat{id}(red_ind_concat,:),1,'omitnan');
    conResp_red_se_stat{id}=red_std/sqrt(length(red_ind_concat));
    
    clear green_std red_std
 
end


figure

subplot(1,2,1) %
errorbar(cons,conResp_red_avrg_stat{pre},conResp_red_se_stat{pre},'k');
hold on
errorbar(cons,conResp_red_avrg_stat{post},conResp_red_se_stat{post},'b');
%title(['SST',' n = ', num2str(length(red_ind_concat))])
%xlabel('contrast') 
ylabel('dF/F') 
xlim([0 1.2])
ylim([.0 .18])
xticks([.25 .5 1])
set(gca, 'TickDir', 'out')
box off

subplot(1,2,2) %
errorbar(cons,conResp_green_avrg_stat{pre},conResp_green_se_stat{pre},'k');
hold on
errorbar(cons,conResp_green_avrg_stat{post},conResp_green_se_stat{post},'b');
%title(['Pyr',' n = ', num2str(length(green_ind_concat))])
%xlabel('contrast') 
set(gca, 'TickDir', 'out')
xlim([0 1.2])
ylim([.0 .18])
xticks([.25 .5 1])
box off


x0=5;
y0=5;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(fullfile(fnout,'Fig_3D.pdf'),'-dpdf');
%% scatterplot with histogram
scatter_signedHypDist(pref_responses_stat_concat, pre,post,red_ind_concat,green_ind_concat)

%% Figure 3D statistics
% prepare data for ANOVA
mouseID=[];
for imouse =1:nSess
   ID_string=mouseNames(imouse);
   thisID = repelem(ID_string,nKeep_concat(imouse));
   mouseID=[mouseID,thisID];
end
mouseID=mouseID';
dfof_stat = horzcat(pref_responses_stat_concat{pre},pref_responses_stat_concat{post});
cellID=(1:nKeep_total)';
cell_type_col=categorical(red_concat)';
dfof_stat_table = array2table(dfof_stat,'VariableNames',{'d1c1','d1c2','d1c3','d2c1','d2c2','d2c3'});
labels_table =table(mouseID,cellID,cell_type_col,'VariableNames',{'mouseID' 'cellID' 'cellType'});
stat_dfof_summary_stat = [labels_table,dfof_stat_table];
SST_stat_dfof = stat_dfof_summary_stat((stat_dfof_summary_stat.cellType=='1'),:);
Pyr_stat_dfof = stat_dfof_summary_stat((stat_dfof_summary_stat.cellType=='0'),:);
clear dfof_stat_table cell_type_col cellID dfof_stat

% run the ANOVA
w = table(categorical([1 1 1 2 2 2 ].'), categorical([1 2 3 1 2 3].'), 'VariableNames', {'DART', 'contrast'}); % within-design
rm_SST_stat = fitrm(SST_stat_dfof, 'd1c1-d2c3 ~ 1', 'WithinDesign', w);
ranova(rm_SST_stat, 'withinmodel', 'DART*contrast')

rm_Pyr_stat = fitrm(Pyr_stat_dfof, 'd1c1-d2c3 ~ 1', 'WithinDesign', w);
ranova(rm_Pyr_stat, 'withinmodel', 'DART*contrast')


% pairwise ttests for dfof response at each contrast for SST cells
[sst_h1, sst_p1]= ttest(pref_responses_stat_concat{pre}(red_ind_concat,1),pref_responses_stat_concat{post}(red_ind_concat,1));
[sst_h2, sst_p2]= ttest(pref_responses_stat_concat{pre}(red_ind_concat,2),pref_responses_stat_concat{post}(red_ind_concat,2));
[sst_h3, sst_p3]= ttest(pref_responses_stat_concat{pre}(red_ind_concat,3),pref_responses_stat_concat{post}(red_ind_concat,3));

%corrected for three tests
sst_pvalues = [(sst_p1*3);(sst_p2*3);(sst_p3*3)];

% pairwise ttests for dfof response at each contrast for Pyr cells
[pyr_h1, pyr_p1]= ttest(pref_responses_stat_concat{pre}(green_ind_concat,1),pref_responses_stat_concat{post}(green_ind_concat,1));
[pyr_h2, pyr_p2]= ttest(pref_responses_stat_concat{pre}(green_ind_concat,2),pref_responses_stat_concat{post}(green_ind_concat,2));
[pyr_h3, pyr_p3]= ttest(pref_responses_stat_concat{pre}(green_ind_concat,3),pref_responses_stat_concat{post}(green_ind_concat,3));

%corrected for three tests
pyr_pvalues = [(pyr_p1*3);(pyr_p2*3);(pyr_p3*3)];

contrasts = cons';
table(contrasts,sst_pvalues,pyr_pvalues)




%% Figure 3F

%make a subset of normalized difference for the SST cells only, then make
% find how many are facilitated or suppressed by more than 1 std from
% baseline
norm_diff_red = norm_diff(:,:,red_ind_concat);
facil_red=norm_diff_red(:,:,:)>=1;
supp_red=norm_diff_red(:,:,:)<=-1;

N=length(red_ind_concat);
facil_table_stat = sum(facil_red(1,:,:),3)/N;
supp_table_stat = sum(supp_red(1,:,:),3)/N;

figure;
subplot(1,2,1)
bar([1,2,3],[supp_table_stat],'FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
xticklabels({'25','50','100'})
title('Suppressed')
ylim([0 .5])
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

subplot(1,2,2)
bar([1,2,3],[facil_table_stat],'FaceColor',"#A8518A",'EdgeColor', [1 1 1])
xticklabels({'25','50','100'})
title('Facilitated')
ylim([0 .5])
%ylabel(["Fraction HTP+ cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'Fig_3F.pdf'),'-dpdf')

%% Figure 3F statistics

%compute chi squares for suppression
%25 vs 50
n1=supp_table_stat(1)*N;
n2=supp_table_stat(2)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%25 vs 10%
n1=supp_table_stat(1)*N;
n2=supp_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%50 vs 100
n1=supp_table_stat(2)*N;
n2=supp_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);

[chi2stat1, chi2stat2, chi2stat3; p1*3, p2*3,p3*3]

clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3

%compute chi squares for facilitation
%25 vs 50
n1=facil_table_stat(1)*N;
n2=facil_table_stat(2)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%25 vs 10%
n1=facil_table_stat(1)*N;
n2=facil_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%50 vs 100
n1=facil_table_stat(2)*N;
n2=facil_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);

[chi2stat1, chi2stat2, chi2stat3; p1*3, p2*3,p3*3]
clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3 n1 n2 x1 x2
%% 3F pyramidal
norm_diff_green = norm_diff(:,:,green_ind_concat);
facil_green=norm_diff_green(:,:,:)>=1;
supp_green=norm_diff_green(:,:,:)<=-1;

N=length(green_ind_concat);
facil_table_stat = sum(facil_green(1,:,:),3)/N;
supp_table_stat = sum(supp_green(1,:,:),3)/N;

figure;
subplot(1,2,1)
bar([1,2,3],[supp_table_stat],'FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
xticklabels({'25','50','100'})
title('Suppressed')
ylim([0 .4])
ylabel(["Fraction Pyr cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

subplot(1,2,2)
bar([1,2,3],[facil_table_stat],'FaceColor',"#A8518A",'EdgeColor', [1 1 1])
xticklabels({'25','50','100'})
title('Facilitated')
ylim([0 .4])
%ylabel(["Fraction HTP+ cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'Fig_3F_pyramidal.pdf'),'-dpdf')


%% Figure 3E green statistics


%compute chi squares for suppression
%25 vs 50
n1=supp_table_stat(1)*N;
n2=supp_table_stat(2)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%25 vs 10%
n1=supp_table_stat(1)*N;
n2=supp_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%50 vs 100
n1=supp_table_stat(2)*N;
n2=supp_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);

[chi2stat1, chi2stat2, chi2stat3; p1*3, p2*3,p3*3]

clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3

%compute chi squares for facilitation
%25 vs 50
n1=facil_table_stat(1)*N;
n2=facil_table_stat(2)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%25 vs 100
n1=facil_table_stat(1)*N;
n2=facil_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%50 vs 100
n1=facil_table_stat(2)*N;
n2=facil_table_stat(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);

[chi2stat1, chi2stat2, chi2stat3; p1*3, p2*3,p3*3]
clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3 n1 n2 x1 x2

%% Figure 3E

figure;
subplot(1,2,1)
boxchart(squeeze(norm_diff(1,:,red_ind_concat))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75],BoxEdgeColor=[0 0 0]);
hold on
scatter([1, 2, 3],squeeze(norm_diff(1,:,red_ind_concat))',20,[.26 .29 .33], 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)

xticklabels({'25','50','100'})
xlabel('Contrast(%)')
ylabel('Normalized difference')
ylim([-10 10])
title('SST')
hold off
set(gca,'TickDir','out')
box off

subplot(1,2,2)
boxchart(squeeze(norm_diff(1,:,green_ind_concat))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75],BoxEdgeColor=[0 0 0]);
hold on
scatter([1, 2, 3],squeeze(norm_diff(1,:,green_ind_concat))',20,[.26 .29 .33], 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)
boxchart(squeeze(norm_diff(1,:,green_ind_concat))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75],BoxEdgeColor=[0 0 0]);
xticklabels({'25','50','100'})
xlabel('Contrast(%)')
ylim([-20 20])
title('Pyr')
hold off
set(gca,'TickDir','out')
box off
x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
%must manually export this figure in order to have it vectorized because of
%the large amount of data 
%%
% Extract data from the matrix for SST cells
data = squeeze(norm_diff(1, :, red_ind_concat));  % Extract data from the specified dimension

% Perform Levene's test
[p, stats] = vartestn(data','TestType','LeveneAbsolute');

% Display the p-value
disp(['Levene''s test for SST cells, p-value: ', num2str(p)]);

% Extract data from the matrix for SST cells
data = squeeze(norm_diff(1, :, green_ind_concat));  % Extract data from the specified dimension

% Perform Levene's test
[p, stats] = vartestn(data','TestType','LeveneAbsolute');

% Display the p-value
disp(['Levene''s test for pyramidal cells, p-value: ', num2str(p)]);

%F-test for equality of variances, asking whether the higher contrast has
%greater variance than the lower contrast, for SST cells
[h,p1,ci,stats1] = vartest2(norm_diff(1,1,red_ind_concat),norm_diff(1,2,red_ind_concat),'Tail','left'); %25 vs 50
[h,p2,ci,stats2] = vartest2(norm_diff(1,1,red_ind_concat),norm_diff(1,3,red_ind_concat),'Tail','left'); %25 vs 100
[h,p3,ci,stats3] = vartest2(norm_diff(1,2,red_ind_concat),norm_diff(1,3,red_ind_concat),'Tail','left'); %50 vs 100


%F-test for equality of variances, asking whether the higher contrast has
%greater variance than the lower contrast, for Pyr cells
[h,p4,ci,stats4] = vartest2(norm_diff(1,1,green_ind_concat),norm_diff(1,2,green_ind_concat),'Tail','left'); %25 vs 50
[h,p5,ci,stats5] = vartest2(norm_diff(1,1,green_ind_concat),norm_diff(1,3,green_ind_concat),'Tail','left'); %25 vs 100
[h,p6,ci,stats6] = vartest2(norm_diff(1,2,green_ind_concat),norm_diff(1,3,green_ind_concat),'Tail','left'); %50 vs 100

format long 
[stats1.fstat, stats2.fstat, stats3.fstat; p1*3, p2*3,p3*3]
[stats4.fstat, stats5.fstat, stats6.fstat; p4*3, p5*3,p6*3]
format short

clear h p1 p2 p3 ci stats1 stats2 stats3 h p4 p5 p6 ci stats4 stats5 stats6
%% Fig S4 histogram
figure;
histogram(noiseCorr_OG_concat{pre}(1,red_ind_concat),'FaceColor','black','BinWidth',.1)
xlabel('Correlation to Pyr activity')
ylabel('# Cells')
set(gca,'TickDir','out')
box off
x0=3;
y0=3;
width=2.5;
height=2;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'Fig_S4A.pdf'),'-dpdf');
%% Fig 4B

% stat high and low R
hi_avrg_stat = cell(1,nd);
low_avrg_stat = cell(1,nd); 
hi_se_stat = cell(1,nd); 
low_se_stat = cell(1,nd);

for id = 1:nd

   for iCon=1:nCon
        %here we are using all SST cells
        hi_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,redHigh,iCon),2);
        high_std=nanstd(tc_trial_avrg_stat_concat{id}(:,redHigh,iCon),[],2);
        hi_se_stat{id}(:,iCon)=high_std/sqrt(length(redHigh));
        
        low_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,redLow,iCon),2);
        low_std=nanstd(tc_trial_avrg_stat_concat{id}(:,redLow,iCon),[],2);
        low_se_stat{id}(:,iCon)=low_std/sqrt(length(redLow));
        
        clear low_std high_std
    end
end
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(hi_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

%Makes a plot for each contrast - 25 contrast used in paper
for iCon = 1:nCon
    figure
    subplot(1,2,1) 
    shadedErrorBar(t,low_avrg_stat{pre}(:,iCon),low_se_stat{pre}(:,iCon),'k');
    hold on
    shadedErrorBar(t,low_avrg_stat{post}(:,iCon),low_se_stat{post}(:,iCon),'b');
    ylim([-.02 .17]);
    hold on
    % line([0,.2],[-.01,-.01],'Color','black','LineWidth',2);
    % hold on
    line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    ylabel('dF/F') 
    xlabel('s') 
    num=length(redLow);
    title([' Weakly correlated',' n = ', num2str(num)])
    set(gca,'XColor', 'none','YColor','none')
    
    
    subplot(1,2,2) 
    ylim([-.02 .17]);
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
    num=length(redHigh);
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
     if iCon==1
        print(fullfile(fnout,'Fig_4B.pdf'),'-dpdf');
     % else %option to print the plots for 50 and 100 contrast as well,
     % rather than just displaying them
     %     print(fullfile(fnout,[num2str(cons(iCon)) 'stat_R_timecourses.pdf']),'-dpdf');
     end
    clear txt1 highRed lowRed
end 



%% Fig S5 A, % As a control, randomly pull SST cells in equal amounts to the number of 
% % weakly and strongly correlated cells.
% nShuff=100;
% 
% hi_shuff_stat = cell(1,nd);
% low_shuff_stat = cell(1,nd);
% shuff_normDiff=nan(2,nCon,nShuff);
% pValues_low=nan(1,nShuff);
% pValues_high=nan(1,nShuff);
% 
% 
% for iShuff = 1:nShuff
%     redLowShuff=randsample(red_ind_concat,length(redLow)); %set of cell IDs randomly selected from among the red cells
%     redHighShuff=setdiff(red_ind_concat,redLowShuff);
% 
%     for id = 1:nd
% 
%        for iCon=1:nCon
%             hi_shuff_stat{id}(:,iCon,iShuff)=nanmean(tc_trial_avrg_stat_concat{id}(:,redHighShuff,iCon),2);
%             low_shuff_stat{id}(:,iCon,iShuff)=nanmean(tc_trial_avrg_stat_concat{id}(:,redLowShuff,iCon),2);
% 
%             shuff_normDiff(1,iCon,iShuff)=nanmean(norm_diff(1,iCon,redLowShuff));
%             shuff_normDiff(2,iCon,iShuff)=nanmean(norm_diff(1,iCon,redHighShuff));
%            % clear redHighShuff redLowShuff
%         end
%     end
% 
% 
% end
% 
% og_means= [mean(norm_diff(1,2,redLow),'omitmissing'),mean(norm_diff(1,2,redHigh),'omitmissing')];
% 
% %plot normalized difference distributions for each shuffle
% figure;
% boxchart(squeeze(shuff_normDiff(:,2,:))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75],BoxEdgeColor=[0 0 0]);
% hold on
% scatter([1, 2],squeeze(shuff_normDiff(:,2,:))',20, [.26 .29 .33],'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)
% scatter([1,2],og_means,'r')
% xticklabels({'Weak corr.', 'Strong corr.'})
% set(gca,'TickDir','out')
% ylim([-1 0.5])
% box off
% x0=5;
% y0=5;
% width=2;
% height=1.5;
% title('Random sample')
% set(gcf,'units','inches','position',[x0,y0,width,height])
% print(fullfile(fnout,'Fig_S5_A.pdf'),'-dpdf');

%% Fig S5 A, % As a control, randomly pull SST cells in equal amounts to the number of 
% weakly and strongly correlated cells.
nShuff=100;

hi_shuff_stat = cell(1,nd);
low_shuff_stat = cell(1,nd);
shuff_normDiff=nan(2,nCon,nShuff);
pValues_low=nan(1,nShuff);
pValues_high=nan(1,nShuff);


for iShuff = 1:nShuff
    redLowShuff=randsample(red_ind_concat,length(redLow)); %set of cell IDs randomly selected from among the red cells
    redHighShuff=setdiff(red_ind_concat,redLowShuff);
    
    for id = 1:nd
    
       for iCon=1:nCon
            hi_shuff_stat{id}(:,iCon,iShuff)=nanmean(tc_trial_avrg_stat_concat{id}(:,redHighShuff,iCon),2);
            low_shuff_stat{id}(:,iCon,iShuff)=nanmean(tc_trial_avrg_stat_concat{id}(:,redLowShuff,iCon),2);
            
            pre_stat_low = NaN(1,length(redLowShuff));
            post_stat_low = NaN(1,length(redLowShuff));
            for iCell=1:length(redLowShuff)
                pre_stat_low(iCell) = mean(pref_allTrials_stat_concat{iCon,pre}{redLowShuff(iCell)});
                post_stat_low(iCell)=mean(pref_allTrials_stat_concat{iCon,post}{redLowShuff(iCell)});
            end
            
            std_pre_stat_low = std(pre_stat_low,[],2);
            shuff_normDiff(1,iCon,iShuff) = (mean(post_stat_low)-mean(pre_stat_low))/std_pre_stat_low;

            pre_stat_high = NaN(1,length(redHighShuff));
            post_stat_high = NaN(1,length(redHighShuff));
            for iCell=1:length(redHighShuff)
                pre_stat_high(iCell) = mean(pref_allTrials_stat_concat{iCon,pre}{redHighShuff(iCell)});
                post_stat_high(iCell)=mean(pref_allTrials_stat_concat{iCon,post}{redHighShuff(iCell)});
            end
            
            std_pre_stat_high = std(pre_stat_high,[],2);
            shuff_normDiff(2,iCon,iShuff) = (mean(post_stat_high)-mean(pre_stat_high))/std_pre_stat_high;

           
        end
    end
clear redHighShuff redLowShuff

end

og_means= [mean(norm_diff(1,2,redLow),'omitmissing'),mean(norm_diff(1,2,redHigh),'omitmissing')];

%plot normalized difference distributions for each shuffle
figure;
boxchart(squeeze(shuff_normDiff(:,2,:))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75],BoxEdgeColor=[0 0 0]);
hold on
scatter([1, 2],squeeze(shuff_normDiff(:,2,:))',20, [.26 .29 .33],'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)
scatter([1,2],og_means,'r')
xticklabels({'Weak corr.', 'Strong corr.'})
set(gca,'TickDir','out')
ylim([-1 0.5])
box off
x0=5;
y0=5;
width=2;
height=1.5;
title('Random sample')
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'Fig_S4B.pdf'),'-dpdf');
%% confidence interval comparison for bottstrapped distributions
CI_level = 95;
ci_low = ((100-CI_level)/2);
ci_high = 100-ci_low;

boot_data1 = squeeze(shuff_normDiff(1,2,:));
boot_data2 = squeeze(shuff_normDiff(2,2,:));

% Calculate confidence intervals
ci1 = prctile(boot_data1, [ci_low, ci_high]);
ci2 = prctile(boot_data2, [ci_low, ci_high]);

fprintf('CI for Data 1: [%f, %f]\n', ci1(1), ci1(2));
fprintf('CI for Data 2: [%f, %f]\n', ci2(1), ci2(2));

% Check for overlap
overlap = ~(ci1(2) < ci2(1) || ci2(2) < ci1(1));
fprintf('Do the confidence intervals overlap? %s\n', string(overlap));
ci1(2)-ci2(1)

%%

%plot timecourses averaged over shuffles
hi_avrg_stat = cell(1,nd);
low_avrg_stat = cell(1,nd); 
hi_std_stat = cell(1,nd);
low_std_stat = cell(1,nd);

for id = 1:nd
   for iCon=1:nCon
        hi_avrg_stat{id}(:,iCon)=mean(hi_shuff_stat{id}(:,iCon,:),3,"omitmissing");
        hi_std_stat{id}(:,iCon)=std(hi_shuff_stat{id}(:,iCon,:),[],3,'omitmissing');
        
        low_avrg_stat{id}(:,iCon)=mean(low_shuff_stat{id}(:,iCon,:),3,"omitmissing");
        low_std_stat{id}(:,iCon)=std(low_shuff_stat{id}(:,iCon,:),[],3,'omitmissing');
        
        clear low_std high_std
    end
end

z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(hi_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

%Makes a plot for each contrast - 50 contrast used in paper
for iCon = 1:nCon
    figure
    subplot(1,2,1) 
    shadedErrorBar(t,low_avrg_stat{pre}(:,iCon),low_std_stat{pre}(:,iCon),'k');
    hold on
    shadedErrorBar(t,low_avrg_stat{post}(:,iCon),low_std_stat{post}(:,iCon),'b');
    ylim([-.02 .17]);
    hold on
    % line([0,.2],[-.01,-.01],'Color','black','LineWidth',2);
    % hold on
    line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    ylabel('dF/F') 
    xlabel('s') 
    num=length(redLowShuff);
    title(['Weakly correlated',' n = ', num2str(num)])
    set(gca,'XColor', 'none','YColor','none')
    
    
    subplot(1,2,2) 
    ylim([-.02 .17]);
    hold on
    shadedErrorBar(t,hi_avrg_stat{pre}(:,iCon),hi_std_stat{pre}(:,iCon),'k');
    hold on
    shadedErrorBar(t,hi_avrg_stat{post}(:,iCon),hi_std_stat{post}(:,iCon),'b');
    hold on
    % line([0,.2],[-.01,-.01],'Color','black','LineWidth',2);
    % hold on
    line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    num=length(redHighShuff);
    title(['Strongly correlated',' n = ', num2str(num)])
    
    xlabel('s') 
    set(gca,'XColor', 'none','YColor','none')
    
    x0=5;
    y0=5;
    width=4;
    height=3;
    set(gcf,'units','inches','position',[x0,y0,width,height])
    set(gca,'XColor', 'none','YColor','none')
    
    sgtitle(['stationary, contrast = ' num2str(cons(iCon))])

    clear txt1 highRed lowRed
end 
    
% %% Fig. S5 B subsampling with replacement 
% sampleFract = 1; %what fraction of the total number of low/high cells we
% %will subsample each time
% 
% nShuff = 100; %how many bootstrap samples to draw
% 
% %how many cells we need to draw each time
% subNlow = round(sampleFract*length(redLow));
% subNhigh = round(sampleFract*length(redHigh));
% 
% low_shuff_stat = cell(1,nd);
% hi_shuff_stat = cell(1,nd);
% Resamp_normDiff=nan(2,nCon,nShuff);
% 
% 
% for iShuff = 1:nShuff
%     redLowResamp = randsample(redLow,subNlow,true);
%     redHighResamp = randsample(redHigh,subNhigh,true);
% 
% 
%     for id = 1:nd
% 
%        for iCon=1:nCon
%             low_shuff_stat{id}(:,iCon,iShuff)=nanmean(tc_trial_avrg_stat_concat{id}(:,redLowResamp,iCon),2);
%             hi_shuff_stat{id}(:,iCon,iShuff)=nanmean(tc_trial_avrg_stat_concat{id}(:,redHighResamp,iCon),2);
% 
%             Resamp_normDiff(1,iCon,iShuff)=nanmean(norm_diff(1,iCon,redLowResamp));
%             Resamp_normDiff(2,iCon,iShuff)=nanmean(norm_diff(1,iCon,redHighResamp));
% 
%         end
%     end
%     clear redHighResamp redLowResamp
% end
% 
% %plot normalized difference distributions for each shuffle
% figure;
% boxchart(squeeze(Resamp_normDiff(:,2,:))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75],BoxEdgeColor=[0 0 0]);
% hold on
% scatter([1, 2],squeeze(Resamp_normDiff(:,2,:))',20,[.26 .29 .33], 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)
% scatter([1,2],og_means,'r')
% xticklabels({'Weak corr.', 'Strong corr.'})
% set(gca,'TickDir','out')
% %ylim([-1 0.25])
% box off
% x0=5;
% y0=5;
% width=2;
% height=1.5;
% title('Subsample with replacement')
% set(gcf,'units','inches','position',[x0,y0,width,height])
% [h, p] = ttest(squeeze(Resamp_normDiff(1,2,:))',squeeze(Resamp_normDiff(2,2,:))')
% print(fullfile(fnout,'Fig_S5_B.pdf'),'-dpdf');
%% Fig. S5 B subsampling with replacement 
sampleFract = 1; %what fraction of the total number of low/high cells we
%will subsample each time

nResamp = 100; %how many bootstrap samples to draw

%how many cells we need to draw each time
subNlow = round(sampleFract*length(redLow));
subNhigh = round(sampleFract*length(redHigh));

low_resamp_stat = cell(1,nd);
hi_resamp_stat = cell(1,nd);
Resamp_normDiff=nan(2,nCon,nResamp);



for iResamp = 1:nResamp
    redLowResamp = randsample(redLow,subNlow,true);
    redHighResamp = randsample(redHigh,subNhigh,true);
    
    for id = 1:nd
    
       for iCon=1:nCon
            hi_resamp_stat{id}(:,iCon,iResamp)=nanmean(tc_trial_avrg_stat_concat{id}(:,redHighResamp,iCon),2);
            low_resamp_stat{id}(:,iCon,iResamp)=nanmean(tc_trial_avrg_stat_concat{id}(:,redLowResamp,iCon),2);
            
            pre_stat_low = NaN(1,length(redLowResamp));
            post_stat_low = NaN(1,length(redLowResamp));
            for iCell=1:length(redLowResamp)
                pre_stat_low(iCell) = mean(pref_allTrials_stat_concat{iCon,pre}{redLowResamp(iCell)});
                post_stat_low(iCell)=mean(pref_allTrials_stat_concat{iCon,post}{redLowResamp(iCell)});
            end
            
            std_pre_stat_low = std(pre_stat_low,[],2);
            Resamp_normDiff(1,iCon,iResamp) = (mean(post_stat_low)-mean(pre_stat_low))/std_pre_stat_low;

            pre_stat_high = NaN(1,length(redHighResamp));
            post_stat_high = NaN(1,length(redHighResamp));
            for iCell=1:length(redHighResamp)
                pre_stat_high(iCell) = mean(pref_allTrials_stat_concat{iCon,pre}{redHighResamp(iCell)});
                post_stat_high(iCell)=mean(pref_allTrials_stat_concat{iCon,post}{redHighResamp(iCell)});
            end
            
            std_pre_stat_high = std(pre_stat_high,[],2);
            Resamp_normDiff(2,iCon,iResamp) = (mean(post_stat_high)-mean(pre_stat_high))/std_pre_stat_high;

           
        end
    end
clear redHighResamp redLowResamp

end


%plot normalized difference distributions for each resample
figure;
boxchart(squeeze(Resamp_normDiff(:,2,:))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75],BoxEdgeColor=[0 0 0]);
hold on
scatter([1, 2],squeeze(Resamp_normDiff(:,2,:))',20,[.26 .29 .33], 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)
scatter([1,2],og_means,'r')
xticklabels({'Weak corr.', 'Strong corr.'})
set(gca,'TickDir','out')
%ylim([-1 0.25])
box off
x0=5;
y0=5;
width=2;
height=1.5;
title('Subsample with replacement')
set(gcf,'units','inches','position',[x0,y0,width,height])
[h, p] = ttest(squeeze(Resamp_normDiff(1,2,:))',squeeze(Resamp_normDiff(2,2,:))')
print(fullfile(fnout,'Fig_S4_c.pdf'),'-dpdf');
%% confidence interval comparison for resampling with replacement
CI_level = 65;
ci_low = ((100-CI_level)/2);
ci_high = 100-ci_low;

resamp_data1 = squeeze(Resamp_normDiff(1,2,:));
resamp_data2 = squeeze(Resamp_normDiff(2,2,:));

% Calculate confidence intervals
ci1 = prctile(resamp_data1, [ci_low, ci_high]);
ci2 = prctile(resamp_data2, [ci_low, ci_high]);

fprintf('CI for Data 1: [%f, %f]\n', ci1(1), ci1(2));
fprintf('CI for Data 2: [%f, %f]\n', ci2(1), ci2(2));
% Check for overlap
overlap = ~(ci1(2) < ci2(1) || ci2(2) < ci1(1));
fprintf('Do the confidence intervals overlap? %s\n', string(overlap));
ci1(2)-ci2(1)

%% plot timecourses of each shuffle
hi_avrg_stat = cell(1,nd);
low_avrg_stat = cell(1,nd); 
hi_std_stat = cell(1,nd);
low_std_stat = cell(1,nd);

for id = 1:nd
   for iCon=1:nCon
        hi_avrg_stat{id}(:,iCon)=mean(hi_shuff_stat{id}(:,iCon,:),3,"omitmissing");
        hi_std_stat{id}(:,iCon)=std(hi_shuff_stat{id}(:,iCon,:),[],3,'omitmissing');
        
        low_avrg_stat{id}(:,iCon)=mean(low_shuff_stat{id}(:,iCon,:),3,"omitmissing");
        low_std_stat{id}(:,iCon)=std(low_shuff_stat{id}(:,iCon,:),[],3,'omitmissing');
        
        clear low_std high_std
    end
end

z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(hi_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

%Makes a plot for each contrast - 50 contrast used in paper
for iCon = 1:nCon
    figure
    subplot(1,2,1) 
    plot(squeeze(low_shuff_stat{pre}(:,iCon,:)),'k');
    hold on
    plot(squeeze(low_shuff_stat{post}(:,iCon,:)),'b');
    ylim([-.02 .17]);
    hold on
    % line([0,.2],[-.01,-.01],'Color','black','LineWidth',2);
    % hold on
    line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    ylabel('dF/F') 
    xlabel('s') 
    title([' Weakly correlated',' n = ', num2str(subNlow)])
    set(gca,'XColor', 'none','YColor','none')
    
    
    subplot(1,2,2) 
    ylim([-.02 .17]);
    hold on
    plot(squeeze(hi_shuff_stat{pre}(:,iCon,:)),'k');
    hold on
    plot(squeeze(hi_shuff_stat{post}(:,iCon,:)),'b');
    hold on
    % line([0,.2],[-.01,-.01],'Color','black','LineWidth',2);
    % hold on
    line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    title([' Strongly correlated',' n = ', num2str(subNhigh)])
    
    xlabel('s') 
    set(gca,'XColor', 'none','YColor','none')
    
    x0=5;
    y0=5;
    width=4;
    height=3;
    set(gcf,'units','inches','position',[x0,y0,width,height])
    set(gca,'XColor', 'none','YColor','none')
    
    sgtitle(['stationary, contrast = ' num2str(cons(iCon))])
     if iCon==2
         print(fullfile(fnout,'Fig_S3_B.pdf'),'-dpdf');
     else
    print(fullfile(fnout,[num2str(cons(iCon)) 'shuf_stat_R_timecourses.pdf']),'-dpdf');
     end
    clear txt1 highRed lowRed
end 
    


hi_avrg_stat = cell(1,nd);
low_avrg_stat = cell(1,nd); 
hi_std_stat = cell(1,nd);
low_std_stat = cell(1,nd);

for id = 1:nd
   for iCon=1:nCon
        hi_avrg_stat{id}(:,iCon)=mean(hi_shuff_stat{id}(:,iCon,:),3,"omitmissing");
        hi_std_stat{id}(:,iCon)=std(hi_shuff_stat{id}(:,iCon,:),[],3,'omitmissing');
        
        low_avrg_stat{id}(:,iCon)=mean(low_shuff_stat{id}(:,iCon,:),3,"omitmissing");
        low_std_stat{id}(:,iCon)=std(low_shuff_stat{id}(:,iCon,:),[],3,'omitmissing');
        
        clear low_std high_std
    end
end

z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(hi_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

%Makes a plot for each contrast - 50 contrast used in paper
for iCon = 1:nCon
    figure
    subplot(1,2,1) 
    shadedErrorBar(t,low_avrg_stat{pre}(:,iCon),low_std_stat{pre}(:,iCon),'k');
    hold on
    shadedErrorBar(t,low_avrg_stat{post}(:,iCon),low_std_stat{post}(:,iCon),'b');
    ylim([-.02 .17]);
    hold on
    % line([0,.2],[-.01,-.01],'Color','black','LineWidth',2);
    % hold on
    line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    ylabel('dF/F') 
    xlabel('s') 
    title([' Weakly correlated',' n = ', num2str(subNlow)])
    set(gca,'XColor', 'none','YColor','none')
    
    
    subplot(1,2,2) 
    ylim([-.02 .17]);
    hold on
    shadedErrorBar(t,hi_avrg_stat{pre}(:,iCon),hi_std_stat{pre}(:,iCon),'k');
    hold on
    shadedErrorBar(t,hi_avrg_stat{post}(:,iCon),hi_std_stat{post}(:,iCon),'b');
    hold on
    % line([0,.2],[-.01,-.01],'Color','black','LineWidth',2);
    % hold on
    line([0,z],[-.015,-.015],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    title([' Strongly correlated',' n = ', num2str(subNhigh)])
    
    xlabel('s') 
    set(gca,'XColor', 'none','YColor','none')
    
    x0=5;
    y0=5;
    width=4;
    height=3;
    set(gcf,'units','inches','position',[x0,y0,width,height])
    set(gca,'XColor', 'none','YColor','none')
    
    sgtitle(['stationary, contrast = ' num2str(cons(iCon))])
    clear txt1 highRed lowRed
end 

%% Find CI for difference (ie permutation test) for shuffled control
shuff_normDiff_subtracted = squeeze(shuff_normDiff(2,2,:))-squeeze(shuff_normDiff(1,2,:));

CI_level = 50;
ci_low = ((100-CI_level)/2);
ci_high = 100-ci_low;

ci_subtracted_shuff = prctile(shuff_normDiff_subtracted, [ci_low, ci_high])
%% Find CI for difference (ie permutation test) for resampled contorl
Resamp_normDiff_subtracted = squeeze(Resamp_normDiff(2,2,:))-squeeze(Resamp_normDiff(1,2,:));

CI_level = 82;
ci_low = ((100-CI_level)/2);
ci_high = 100-ci_low;

ci_subtracted_resamp = prctile(Resamp_normDiff_subtracted, [ci_low, ci_high])

%% Find cohen's D for shuffled control
effect_shuffled = meanEffectSize(squeeze(shuff_normDiff(2,2,:)),squeeze(shuff_normDiff(1,2,:)),Effect="cohen")
%% Find cohen's D for resampled control
effect_resamp = meanEffectSize(squeeze(Resamp_normDiff(2,2,:)),squeeze(Resamp_normDiff(1,2,:)),Effect="cohen")
    
%% Fig 4C - contrast response split by correlation value

conResp_redHigh_avrg_stat = cell(1,nd); %this will be the average across all redHigh cells - a single line
conResp_redLow_avrg_stat = cell(1,nd); %same for redLow
conResp_redHigh_se_stat = cell(1,nd); %this will be the se across all redHigh cells
conResp_redLow_se_stat = cell(1,nd); %same for redLow

for id = 1:nd
    conResp_redHigh_avrg_stat{id}=mean(pref_responses_stat_concat{id}(redHigh,:),1,'omitnan');
    redHigh_std=std(pref_responses_stat_concat{id}(redHigh,:),1,'omitnan');
    conResp_redHigh_se_stat{id}=redHigh_std/sqrt(length(redHigh));
    
    conResp_redLow_avrg_stat{id}=mean(pref_responses_stat_concat{id}(redLow,:),1,'omitnan');
    redLow_std=std(pref_responses_stat_concat{id}(redLow,:),1,'omitnan');
    conResp_redLow_se_stat{id}=redLow_std/sqrt(length(redLow));
    
    clear redHigh_std redLow_std
 
end

figure
subplot(1,2,1) %
errorbar(cons,conResp_redLow_avrg_stat{pre},conResp_redLow_se_stat{pre},'k');
hold on
errorbar(cons,conResp_redLow_avrg_stat{post},conResp_redLow_se_stat{post},'b');
title(['Weak Corr',' n = ', num2str(length(redLow))])
xlabel('contrast') 
ylabel('dF/F, pref ori') 
xlim([0 1.2])
ylim([0 .15])
xticks([.25 .5 1])
set(gca, 'TickDir', 'out')
box off

subplot(1,2,2) %
errorbar(cons,conResp_redHigh_avrg_stat{pre},conResp_redHigh_se_stat{pre},'k');
hold on
errorbar(cons,conResp_redHigh_avrg_stat{post},conResp_redHigh_se_stat{post},'b');
title(['Strong Corr,' 'n = ', num2str(length(redHigh))])

xlabel('contrast') 
set(gca, 'TickDir', 'out')
xlim([0 1.2])
ylim([0 .15])
xticks([.25 .5 1])
box off

x0=5;
y0=5;
width=2.5;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(fullfile(fnout,'Fig_3C.pdf'),'-dpdf');

%% Figure 4C statistics
%from the full dataframe, extract rows for stationary trials for weakly and
%strongly correlated SST cells
SST_low_stat_dfof = SST_stat_dfof(ismember(SST_stat_dfof.cellID,redLow),:);
SST_high_stat_dfof = SST_stat_dfof(ismember(SST_stat_dfof.cellID,redHigh),:);

fullCorrTable = vertcat(SST_low_stat_dfof,SST_high_stat_dfof);
fullCorrTable.corrType = categorical([repmat(0,length(redLow),1);repmat(1,length(redHigh),1)]);

% run an ANOVA for corr type X contrast, baseline day only
w = table(categorical([1 2 3].'), 'VariableNames', {'contrast'}); % within-design
rm_fullCorr = fitrm(fullCorrTable, 'd1c1-d1c3 ~ corrType', 'WithinDesign', w)
ranova(rm_fullCorr, 'withinmodel', 'contrast')

% run the ANOVAs on low and high corr seperately
w = table(categorical([1 1 1 2 2 2 ].'), categorical([1 2 3 1 2 3].'), 'VariableNames', {'DART', 'contrast'}); % within-design
rm_SST_low = fitrm(SST_low_stat_dfof, 'd1c1-d2c3 ~ 1', 'WithinDesign', w)
ranova(rm_SST_low, 'withinmodel', 'DART*contrast')

rm_SST_high = fitrm(SST_high_stat_dfof, 'd1c1-d2c3 ~ 1', 'WithinDesign', w)
ranova(rm_SST_high, 'withinmodel', 'DART*contrast')


% pairwise ttests for dfof response at 50 contrast for each correlation
% category
[h, p1]= ttest(pref_responses_stat_concat{pre}(redLow,2),pref_responses_stat_concat{post}(redLow,2));
[h, p2]= ttest(pref_responses_stat_concat{pre}(redHigh,2),pref_responses_stat_concat{post}(redHigh,2));

% pairwise ttests for dfof response at each contrast for weakly correlated cells
[low_h1, low_p1]= ttest(pref_responses_stat_concat{pre}(redLow,1),pref_responses_stat_concat{post}(redLow,1));
[low_h2, low_p2]= ttest(pref_responses_stat_concat{pre}(redLow,2),pref_responses_stat_concat{post}(redLow,2));
[low_h3, low_p3]= ttest(pref_responses_stat_concat{pre}(redLow,3),pref_responses_stat_concat{post}(redLow,3));

%corrected for three tests
low_pvalues = [(low_p1*3);(low_p2*3);(low_p3*3)];

% pairwise ttests for dfof response at each contrast for strongly correlated cells
[high_h1, high_p1]= ttest(pref_responses_stat_concat{pre}(redHigh,1),pref_responses_stat_concat{post}(redHigh,1));
[high_h2, high_p2]= ttest(pref_responses_stat_concat{pre}(redHigh,2),pref_responses_stat_concat{post}(redHigh,2));
[high_h3, high_p3]= ttest(pref_responses_stat_concat{pre}(redHigh,3),pref_responses_stat_concat{post}(redHigh,3));

%corrected for three tests
high_pvalues = [(high_p1*3);(high_p2*3);(high_p3*3)]
contrasts = cons';

table(contrasts,low_pvalues,high_pvalues)

%% suppressed and facilitated for weak and strong correlation
%make a subset of normalized difference for the SST cells only, then make
% find how many are facilitated or suppressed by more than 1 std from
% baseline
%reselect norm_diff)red to be only the "red all" subset
norm_diff_redHigh = norm_diff(:,:,redHigh);
facil_red_high=norm_diff_redHigh(:,:,:)>=1;
supp_red_high=norm_diff_redHigh(:,:,:)<=-1;

N1=length(redHigh);
facil_table_high = sum(facil_red_high(1,:,:),3)/N1;
supp_table_high = sum(supp_red_high(1,:,:),3)/N1;

norm_diff_redLow = norm_diff(:,:,redLow);
facil_red_low=norm_diff_redLow(:,:,:)>=1;
supp_red_low=norm_diff_redLow(:,:,:)<=-1;
N2=length(redLow);
facil_table_low = sum(facil_red_low(1,:,:),3)/N2;
supp_table_low = sum(supp_red_low(1,:,:),3)/N2;

figure;
subplot(1,2,1)
b=bar([1,2,3],[supp_table_low; supp_table_high],'grouped','FaceColor',"#00ffff",'EdgeColor', [1 1 1]);
b(1).FaceColor="#70D0F6"
b(2).FaceColor="#0C8ABB"
ylim([0 .6])
title('Suppressed')
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

subplot(1,2,2)
b=bar([1,2,3],[facil_table_low; facil_table_high],'FaceColor',"#a329cc",'EdgeColor', [1 1 1]);
b(1).FaceColor="#C983B1"
b(2).FaceColor="#883367"
xticklabels({'25','50','100'})
ylim([0 .6])
title('Facilitated')
ylabel(["Fraction HTP+ cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'Facil_supp_byR.pdf'),'-dpdf')

%% Fig 5A - timecourses for running trials

tc_red_avrg_stat = cell(1,nd);
tc_red_se_stat = cell(1,nd); 
tc_red_avrg_loc = cell(1,nd);
tc_red_se_loc = cell(1,nd); 

for id = 1:nd
    for iCon=1:nCon
        
    tc_red_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,red_all,iCon),2);
    red_std_stat=nanstd(tc_trial_avrg_stat_concat{id}(:,red_all,iCon),[],2);
    tc_red_se_stat{id}(:,iCon)=red_std_stat/sqrt(length(red_all));
    
    tc_red_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,red_all,iCon),2);
    red_std_loc=nanstd(tc_trial_avrg_loc_concat{id}(:,red_all,iCon),[],2);
    tc_red_se_loc{id}(:,iCon)=red_std_loc/sqrt(length(red_all));
    
    clear red_std_stat red_std_loc
    end
end

%create a 2-second stimulus marker
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_red_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);
positions=[1,2;3,4;5,6];
figure
for iCon = 1:nCon
p1=positions(iCon,1);
p2=positions(iCon,2);

subplot(3,2,p1) 
ylim([-.02 .35]);
hold on
shadedErrorBar(t,tc_red_avrg_stat{pre}(:,iCon),tc_red_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_stat{post}(:,iCon),tc_red_se_stat{post}(:,iCon),'b','transparent');
hold on

if iCon==1
    title("Stationary")
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
elseif iCon==3
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
end
set(gca,'XColor', 'none','YColor','none')

subplot(3,2,p2) 
shadedErrorBar(t,tc_red_avrg_loc{pre}(:,iCon),tc_red_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_loc{post}(:,iCon),tc_red_se_loc{post}(:,iCon),'b');
ylim([-.02 .35]);
hold on

if iCon==1
    title("Running")
elseif iCon==3
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
end
set(gca,'XColor', 'none','YColor','none')

sgtitle(['SST',' n = ', num2str(length(red_all))])

x0=5;
y0=0;
width=4;
height=9;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')


end  
print(fullfile(fnout,'Fig_5A.pdf'),'-dpdf');
%% scatterplots for matched cells
scatter_signedHypDist(pref_responses_stat_concat, pre,post,red_all,green_all,'stationaryMatched')

scatter_signedHypDist(pref_responses_loc_concat, pre,post,red_all,green_all,'runningMatched')

%% Fig 5C- timecourses for Pyr cells running trials
tc_green_avrg_stat = cell(1,nd);
tc_green_se_stat = cell(1,nd); 
tc_green_avrg_loc = cell(1,nd);
tc_green_se_loc = cell(1,nd); 

for id = 1:nd
    for iCon=1:nCon
        
    tc_green_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,green_all,iCon),2);
    green_std_stat=nanstd(tc_trial_avrg_stat_concat{id}(:,green_all,iCon),[],2);
    tc_green_se_stat{id}(:,iCon)=green_std_stat/sqrt(length(green_all));
    
    tc_green_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,green_all,iCon),2);
    green_std_loc=nanstd(tc_trial_avrg_loc_concat{id}(:,green_all,iCon),[],2);
    tc_green_se_loc{id}(:,iCon)=green_std_loc/sqrt(length(green_all));
    
    clear green_std_stat green_std_loc
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
ylim([-.02 .33]);
hold on
shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon),tc_green_se_stat{pre}(:,iCon),'--k');
hold on
shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon),tc_green_se_stat{post}(:,iCon),'--b','transparent');
hold on
if iCon==1
    title("Stationary")
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
elseif iCon==3
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
end
set(gca,'XColor', 'none','YColor','none')

subplot(3,2,p2) 
shadedErrorBar(t,tc_green_avrg_loc{pre}(:,iCon),tc_green_se_loc{pre}(:,iCon),'--k');
hold on
shadedErrorBar(t,tc_green_avrg_loc{post}(:,iCon),tc_green_se_loc{post}(:,iCon),'--b');
ylim([-.02 .33]);
hold on
if iCon==1
    title("Running")
elseif iCon==3
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
end
set(gca,'XColor', 'none','YColor','none')

sgtitle(['Pyr',' n = ', num2str(length(green_all))])

x0=5;
y0=0;
width=4;
height=9;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')


end  
print(fullfile(fnout,'Fig_5C.pdf'),'-dpdf');

%% Figure 5B
conResp_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_stat = cell(1,nd); %same for red
conResp_green_se_stat = cell(1,nd); %this will be the se across all green cells
conResp_red_se_stat = cell(1,nd); %same for red

conResp_green_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_loc = cell(1,nd); %same for red
conResp_green_se_loc = cell(1,nd); %this will be the se across all green cells
conResp_red_se_loc = cell(1,nd); %same for red


for id = 1:nd
   
        
    conResp_green_avrg_stat{id}=mean(pref_responses_stat_concat{id}(green_all,:),1,'omitnan');
    green_std=std(pref_responses_stat_concat{id}(green_all,:),1,'omitnan');
    conResp_green_se_stat{id}=green_std/sqrt(length(green_all));
    
    conResp_red_avrg_stat{id}=mean(pref_responses_stat_concat{id}(red_all,:),1,'omitnan');
    red_std=std(pref_responses_stat_concat{id}(red_all,:),1,'omitnan');
    conResp_red_se_stat{id}=red_std/sqrt(length(red_all));
    
    clear green_std red_std
    
    conResp_green_avrg_loc{id}=nanmean(pref_responses_loc_concat{id}(green_all ,:),1);
    green_std=nanstd(pref_responses_loc_concat{id}(green_all,:),1);
    conResp_green_se_loc{id}=green_std/sqrt(length(green_all));
    
    conResp_red_avrg_loc{id}=nanmean(pref_responses_loc_concat{id}(red_all,:),1);
    red_std=nanstd(pref_responses_loc_concat{id}(red_all,:),1);
    conResp_red_se_loc{id}=red_std/sqrt(length(red_all));
    
    clear green_std red_std
 
end

figure
subplot(1,2,1)
errorbar(cons,conResp_red_avrg_stat{pre},conResp_red_se_stat{pre},'k');
hold on
errorbar(cons,conResp_red_avrg_stat{post},conResp_red_se_stat{post},'b');
title(['Stationary'])
%ylabel('dF/F') 
%xlabel('contrast') 
xlim([0 1.2])
ylim([0 .3])
xticks([.25 .50 1.00])
xticklabels([25 50 100])
set(gca, 'TickDir', 'out')
box off


subplot(1,2,2) 
errorbar(cons,conResp_red_avrg_loc{pre},conResp_red_se_loc{pre},'k');
hold on
errorbar(cons,conResp_red_avrg_loc{post},conResp_red_se_loc{post},'b');
title(['Running'])
ylim([0 .3])
xlim([0 1.2])
xticks([.25 .50 1.00])
xticklabels([25 50 100])
%xlabel('contrast') 
set(gca, 'TickDir', 'out')
box off

x0=5;
y0=5;
width=2.5;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(fullfile(fnout,'Fig_4B.pdf'),'-dpdf');
%% Fig 5B statistics

% prepare data for ANOVA
all_cells = union(red_all,green_all);

dfof_stat = horzcat(pref_responses_stat_concat{pre},pref_responses_stat_concat{post});
dfof_loc = horzcat(pref_responses_loc_concat{pre},pref_responses_loc_concat{post});

cellID=(1:nKeep_total)';
cell_type_col=categorical(red_concat)';

dfof_full = horzcat(dfof_stat,dfof_loc);
dfof_table = array2table(dfof_full,'VariableNames',{'Sd1c1','Sd1c2','Sd1c3', ...
    'Sd2c1','Sd2c2','Sd2c3','Rd1c1','Rd1c2','Rd1c3','Rd2c1','Rd2c2','Rd2c3'});

labels_table =table(mouseID,cellID,cell_type_col,'VariableNames',{'mouseID' 'cellID' 'cellType'});
dfof_summary_full = [labels_table,dfof_table];
matched_dfof_summary = dfof_summary_full(ismember(dfof_summary_full.cellID,all_cells),:);


SST_matched_dfof = matched_dfof_summary((matched_dfof_summary.cellType=='1'),:);
Pyr_matched_dfof = matched_dfof_summary((matched_dfof_summary.cellType=='0'),:);
clear dfof_stat_table cell_type_col cellID dfof_stat_matched

% run an ANOVA on the full model for each cell type
DART = categorical([1 1 1 2 2 2 1 1 1 2 2 2].');
contrast = categorical([1 2 3 1 2 3 1 2 3 1 2 3].');
behState = categorical([1 1 1 1 1 1 2 2 2 2 2 2].');

% prepare new categorical variables
w = table(DART, contrast, behState, ...
    'VariableNames', {'DART', 'contrast', 'behState'});% within-design

rm_SST_matched = fitrm(SST_matched_dfof, 'Sd1c1-Rd2c3 ~ 1', 'WithinDesign', w);
ranova(rm_SST_matched, 'withinmodel', 'behState*DART*contrast')

rm_Pyr_matched = fitrm(Pyr_matched_dfof, 'Sd1c1-Rd2c3 ~ 1', 'WithinDesign', w);
ranova(rm_Pyr_matched, 'withinmodel', 'behState*DART*contrast')

% models for stationary and running seperately 
w = table(categorical([1 1 1 2 2 2 ].'), categorical([1 2 3 1 2 3].'), 'VariableNames', {'DART', 'contrast'}); % within-design

rm_SST_stat = fitrm(SST_matched_dfof, 'Sd1c1-Sd2c3 ~ 1', 'WithinDesign', w)
ranova(rm_SST_stat, 'withinmodel', 'DART*contrast')

rm_SST_loc = fitrm(SST_matched_dfof, 'Rd1c1-Rd2c3 ~ 1', 'WithinDesign', w)
ranova(rm_SST_loc, 'withinmodel', 'DART*contrast')

rm_Pyr_stat = fitrm(Pyr_matched_dfof, 'Sd1c1-Sd2c3 ~ 1', 'WithinDesign', w)
ranova(rm_Pyr_stat, 'withinmodel', 'DART*contrast')

rm_Pyr_loc = fitrm(Pyr_matched_dfof, 'Rd1c1-Rd2c3 ~ 1', 'WithinDesign', w)
ranova(rm_Pyr_loc, 'withinmodel', 'DART*contrast')


% pairwise ttests for dfof response at each contrast for SST cells
[sst_h1, sst_p1]= ttest(pref_responses_stat_concat{pre}(red_all,1),pref_responses_stat_concat{post}(red_all,1));
[sst_h2, sst_p2]= ttest(pref_responses_stat_concat{pre}(red_all,2),pref_responses_stat_concat{post}(red_all,2));
[sst_h3, sst_p3]= ttest(pref_responses_stat_concat{pre}(red_all,3),pref_responses_stat_concat{post}(red_all,3));

%corrected for three tests
sst_pvalues_stat = [(sst_p1*3);(sst_p2*3);(sst_p3*3)];
contrasts = cons';
table(contrasts,sst_pvalues_stat)

% pairwise ttests for dfof response at each contrast for SST cells
[sst_h1, sst_p1]= ttest(pref_responses_loc_concat{pre}(red_all,1),pref_responses_loc_concat{post}(red_all,1));
[sst_h2, sst_p2]= ttest(pref_responses_loc_concat{pre}(red_all,2),pref_responses_loc_concat{post}(red_all,2));
[sst_h3, sst_p3]= ttest(pref_responses_loc_concat{pre}(red_all,3),pref_responses_loc_concat{post}(red_all,3));

%corrected for three tests
sst_pvalues_loc = [(sst_p1*3);(sst_p2*3);(sst_p3*3)];
contrasts = cons';
table(contrasts,sst_pvalues_loc)

%% Supplemental figure related to 4B
%stationary and running contrast response for Pyr cells
figure
subplot(1,2,1) %for the first day
errorbar(cons,conResp_green_avrg_stat{pre},conResp_green_se_stat{pre},'--k');
hold on
errorbar(cons,conResp_green_avrg_stat{post},conResp_green_se_stat{post},'--b');
%title(['-HTP',' n = ', num2str(length(green_all))])
ylabel('dF/F') 
xlabel('contrast') 
set(gca, 'TickDir', 'out')
ylim([0 .3])
xlim([0 1.2])
xticks([.25 .50 1.00])
xticklabels([25 50 100])
%xlabel('contrast') 
set(gca, 'TickDir', 'out')
box off


subplot(1,2,2) %for the first day
errorbar(cons,conResp_green_avrg_loc{pre},conResp_green_se_loc{pre},'--k');
hold on
errorbar(cons,conResp_green_avrg_loc{post},conResp_green_se_loc{post},'--b');
%title(['-HTP',' n = ', num2str(length(green_all))])
xlabel('contrast') 
ylim([0 .3])
xlim([0 1.2])
xticks([.25 .50 1.00])
xticklabels([25 50 100])
%xlabel('contrast') 
set(gca, 'TickDir', 'out')
box off

x0=5;
y0=5;
width=2.8;
height=1.45;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(fullfile(fnout,'Fig_S6B.pdf'),'-dpdf');
%%
% pairwise ttests for dfof response at each contrast for SST cells
[pyr_h1, pyr_p1]= ttest(pref_responses_stat_concat{pre}(green_all,1),pref_responses_stat_concat{post}(green_all,1));
[pyr_h2, pyr_p2]= ttest(pref_responses_stat_concat{pre}(green_all,2),pref_responses_stat_concat{post}(green_all,2));
[pyr_h3, pyr_p3]= ttest(pref_responses_stat_concat{pre}(green_all,3),pref_responses_stat_concat{post}(green_all,3));

%corrected for three tests
pyr_pvalues_stat = [(pyr_p1*3);(pyr_p2*3);(pyr_p3*3)];
contrasts = cons';
table(contrasts,pyr_pvalues_stat)

% pairwise ttests for dfof response at each contrast for SST cells
[pyr_h1, pyr_p1]= ttest(pref_responses_loc_concat{pre}(green_all,1),pref_responses_loc_concat{post}(green_all,1));
[pyr_h2, pyr_p2]= ttest(pref_responses_loc_concat{pre}(green_all,2),pref_responses_loc_concat{post}(green_all,2));
[pyr_h3, pyr_p3]= ttest(pref_responses_loc_concat{pre}(green_all,3),pref_responses_loc_concat{post}(green_all,3));

%corrected for three tests
pyr_pvalues_loc = [(pyr_p1*3);(pyr_p2*3);(pyr_p3*3)];
contrasts = cons';
table(contrasts,pyr_pvalues_loc)
x0=5;
y0=5;
width=2.5;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(fullfile(fnout,'Fig_S5B.pdf'),'-dpdf');
%% Fig 4C - Supp/facil stationary vs. running

%make a subset of normalized difference for the SST cells only, then make
% find how many are facilitated or suppressed by more than 1 std from
% baseline
%reselect norm_diff)red to be only the "red all" subset
norm_diff_red = norm_diff(:,:,red_all);

facil_red=norm_diff_red(:,:,:)>=1;
supp_red=norm_diff_red(:,:,:)<=-1;

N=length(red_all);
facil_table_stat = sum(facil_red(1,:,:),3)/N;
supp_table_stat = sum(supp_red(1,:,:),3)/N;
facil_table_loc = sum(facil_red(2,:,:),3)/N;
supp_table_loc = sum(supp_red(2,:,:),3)/N;

figure;
subplot(1,2,1)
b=bar([1,2,3],[supp_table_stat; supp_table_loc],'grouped','FaceColor',"#00ffff",'EdgeColor', [1 1 1]);
b(1).FaceColor="#70D0F6"
b(2).FaceColor="#0C8ABB"
ylim([0 .6])
title('Suppressed')
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

subplot(1,2,2)
b=bar([1,2,3],[facil_table_stat; facil_table_loc],'FaceColor',"#a329cc",'EdgeColor', [1 1 1]);
b(1).FaceColor="#C983B1"
b(2).FaceColor="#883367"
xticklabels({'25','50','100'})
ylim([0 .6])
title('Facilitated')
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'Fig_4C.pdf'),'-dpdf')



%% Fig 4C stats 

%compute chi squares for suppression, stationary vs running
%25% contrast
% N previously set to be the number of SST+ cells. The same cells are
% included in stationary and running.
n1=supp_table_stat(1)*N;
n2=supp_table_loc(1)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%50
n1=supp_table_stat(2)*N;
n2=supp_table_loc(2)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%100
n1=supp_table_stat(3)*N;
n2=supp_table_loc(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);


%[chi2stat1, chi2stat2, chi2stat3; p1*3, p2*3,p3*3]
[chi2stat1, chi2stat2, chi2stat3; p1, p2,p3]

clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3

%compute chi squares for facilitation
%25
n1=facil_table_stat(1)*N;
n2=facil_table_loc(1)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%50
n1=facil_table_stat(2)*N;
n2=facil_table_loc(2)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%100
n1=facil_table_stat(3)*N;
n2=facil_table_loc(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);

[chi2stat1, chi2stat2, chi2stat3; p1, p2,p3]
clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3 n1 n2 x1 x2

%% facilitated and suppressed for pyr cells
%make a subset of normalized difference for the SST cells only, then make
% find how many are facilitated or suppressed by more than 1 std from
% baseline
%reselect norm_diff)green to be only the "green all" subset
norm_diff_green = norm_diff(:,:,green_all);

facil_green=norm_diff_green(:,:,:)>=1;
supp_green=norm_diff_green(:,:,:)<=-1;

N=length(green_all);
facil_table_stat = sum(facil_green(1,:,:),3)/N;
supp_table_stat = sum(supp_green(1,:,:),3)/N;
facil_table_loc = sum(facil_green(2,:,:),3)/N;
supp_table_loc = sum(supp_green(2,:,:),3)/N;

figure;
subplot(1,2,1)
b=bar([1,2,3],[supp_table_stat; supp_table_loc],'grouped','FaceColor',"#00ffff",'EdgeColor', [1 1 1]);
b(1).FaceColor="#70D0F6"
b(2).FaceColor="#0C8ABB"
ylim([0 .6])
title('Suppressed')
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

subplot(1,2,2)
b=bar([1,2,3],[facil_table_stat; facil_table_loc],'FaceColor',"#a329cc",'EdgeColor', [1 1 1]);
b(1).FaceColor="#C983B1"
b(2).FaceColor="#883367"
xticklabels({'25','50','100'})
ylim([0 .6])
title('Facilitated')
ylabel(["Fraction HTP+ cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
%compute chi squares for suppression, stationary vs running
%25% contrast

n1=int32(supp_table_stat(1)*N);
n2=int32(supp_table_loc(1)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%50
n1=int32(supp_table_stat(2)*N);
n2=int32(supp_table_loc(2)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%100
n1=int32(supp_table_stat(3)*N);
n2=int32(supp_table_loc(3)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);



%[chi2stat1, chi2stat2, chi2stat3; p1*3, p2*3,p3*3]
[chi2stat1, chi2stat2, chi2stat3; p1, p2,p3]

clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3

%compute chi squares for facilitation
%25
n1=int32(supp_table_stat(1)*N);
n2=int32(facil_table_loc(1)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat1,p1] = crosstab(x1,x2);

%50
n1=int32(facil_table_stat(2)*N);
n2=int32(facil_table_loc(2)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat2,p2] = crosstab(x1,x2);

%100
n1=int32(facil_table_stat(3)*N);
n2=int32(facil_table_loc(3)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2stat3,p3] = crosstab(x1,x2);

[chi2stat1, chi2stat2, chi2stat3; p1, p2,p3]
clear h p1 p2 p3 chi2stat1 chi2stat2 chi2stat3 n1 n2 x1 x2



%% Fig Sx A, comparing the preferred direction across days

change_pref = pref_dir_concat{pre}-pref_dir_concat{post};
change_pref=deg2rad(change_pref);
figure

subplot(1,2,1)
polarhistogram(change_pref(red_ind_concat),'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5])
title('SST')
set(gca, 'TickDir', 'out')

subplot(1,2,2)
polarhistogram(change_pref(green_ind_concat),'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5])
title('Pyr')
set(gca, 'TickDir', 'out')
x0=5;
y0=5;
width=3;
height=2;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,['Fig_Sx_B.pdf']),'-dpdf','-bestfit')
%% Figure Sx B, direction tuning and DART effect
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
subplot(1,2,1)
errorbar(dirs_for_plotting,red_dir_avrg_stat{pre},red_dir_se_stat{pre},'k')
hold on
errorbar(dirs_for_plotting,red_dir_avrg_stat{post},red_dir_se_stat{post},'b')
title(['SST'])
ylabel('dF/F')
set(gca, 'TickDir', 'out')
axis square
box off

subplot(1,2,2)
errorbar(dirs_for_plotting,green_dir_avrg_stat{pre},green_dir_se_stat{pre},'k')
hold on
errorbar(dirs_for_plotting,green_dir_avrg_stat{post},green_dir_se_stat{post},'b')
title(['Pyr'])
set(gca, 'TickDir', 'out')
axis square
box off
%ylim([-0.01 .2])

x0=5;
y0=5;
width=3;
height=2;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,['Fig_Sx_B.pdf']),'-dpdf','-bestfit')

%% Fig Sx C all keep cells at non-preferred directions

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
subplot(1,2,1) 

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


subplot(1,2,2) %
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
     if iCon==2
        print(fullfile(fnout,'Fig_Sx_C.pdf'),'-dpdf');
     % else %option to print the plots for 25 and 100 contrast as well,
     % rather than just displaying them
     %     print(fullfile(fnout,[num2str(cons(iCon)) 'stat_R_timecourses.pdf']),'-dpdf');
     end
end 






%% Fig 5 B, SST cell timecourses split by pupil size, and SX B, Pyr responses split by pupil size

tc_green_avrg_stat_large = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_stat_large = cell(1,nd); %same for red
tc_green_se_stat_large = cell(1,nd); %this will be the se across all green cells
tc_red_se_stat_large = cell(1,nd); %same for red
tc_green_avrg_stat_small = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_stat_small = cell(1,nd); %same for red
tc_green_se_stat_small = cell(1,nd); %this will be the se across all green cells
tc_red_se_stat_small = cell(1,nd); %same for red



    for id = 1:nd
        for iCon=1:nCon
    
    
        tc_green_avrg_stat_large{id}(:,iCon)=nanmean(tc_trial_avrg_stat_largePupil_concat{id}(:,have_allPupil_green,iCon),2);
        green_std=nanstd(tc_trial_avrg_stat_largePupil_concat{id}(:,have_allPupil_green,iCon),[],2);
        tc_green_se_stat_large{id}(:,iCon)=green_std/sqrt(length(have_allPupil_green));
        
        tc_red_avrg_stat_large{id}(:,iCon)=nanmean(tc_trial_avrg_stat_largePupil_concat{id}(:,have_allPupil_red,iCon),2);
        red_std=nanstd(tc_trial_avrg_stat_largePupil_concat{id}(:,have_allPupil_red,iCon),[],2);
        tc_red_se_stat_large{id}(:,iCon)=red_std/sqrt(length(have_allPupil_red));
        
        clear green_std red_std
    
        tc_green_avrg_stat_small{id}(:,iCon)=nanmean(tc_trial_avrg_stat_smallPupil_concat{id}(:,have_allPupil_green,iCon),2);
        green_std=nanstd(tc_trial_avrg_stat_smallPupil_concat{id}(:,have_allPupil_green,iCon),[],2);
        tc_green_se_stat_small{id}(:,iCon)=green_std/sqrt(length(have_allPupil_green));
        
        tc_red_avrg_stat_small{id}(:,iCon)=nanmean(tc_trial_avrg_stat_smallPupil_concat{id}(:,have_allPupil_red,iCon),2);
        red_std=nanstd(tc_trial_avrg_stat_smallPupil_concat{id}(:,have_allPupil_red,iCon),[],2);
        tc_red_se_stat_small{id}(:,iCon)=red_std/sqrt(length(have_allPupil_red));
        
        clear green_std red_std
        end
    end
    
    
    z=double(nOn)/double(frame_rate);
    
    %create a time axis in seconds
    t=1:(size(tc_green_avrg_stat_large{1,1,1},1));
    t=(t-(double(stimStart)-1))/double(frame_rate);
    
    positions=[1,2;3,4;5,6];
    
    figure
    for iCon = 1:nCon
    p1=positions(iCon,1);
    p2=positions(iCon,2);
    
    subplot(3,2,p1) 
    ylim([-.02 .2]);
    hold on
    shadedErrorBar(t,tc_red_avrg_stat_small{pre}(:,iCon),tc_red_se_stat_small{pre}(:,iCon),'k');
    hold on
    shadedErrorBar(t,tc_red_avrg_stat_small{post}(:,iCon),tc_red_se_stat_small{post}(:,iCon),'b','transparent');
    hold on
    if iCon==1
        title("Small pupil")
        line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    elseif iCon==3
        line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    end
    set(gca,'XColor', 'none','YColor','none')
    
    subplot(3,2,p2) 
    
    hold on
    shadedErrorBar(t,tc_red_avrg_stat_large{pre}(:,iCon),tc_red_se_stat_large{pre}(:,iCon),'k');
    hold on
    shadedErrorBar(t,tc_red_avrg_stat_large{post}(:,iCon),tc_red_se_stat_large{post}(:,iCon),'b','transparent');
    hold on
    if iCon==1
        title("Large pupil")
    elseif iCon==3
        line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    end
    set(gca,'XColor', 'none','YColor','none')
    ylim([-.02 .2]);
    sgtitle(['SST',' n = ', num2str(length(have_allPupil_red))])
    
    x0=5;
    y0=0;
    width=4;
    height=9;
    set(gcf,'units','inches','position',[x0,y0,width,height])
    set(gca,'XColor', 'none','YColor','none')
    
    
    end  

print(fullfile(fnout,'Fig_5B.pdf'),'-dpdf');

figure
for iCon = 1:nCon
p1=positions(iCon,1);
p2=positions(iCon,2);

subplot(3,2,p1) 
ylim([-.02 .25]);
hold on
shadedErrorBar(t,tc_green_avrg_stat_small{pre}(:,iCon),tc_green_se_stat_small{pre}(:,iCon),'--k');
hold on
shadedErrorBar(t,tc_green_avrg_stat_small{post}(:,iCon),tc_green_se_stat_small{post}(:,iCon),'--b','transparent');
hold on
if iCon==1
    title("Small pupil")
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
elseif iCon==3
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
end
set(gca,'XColor', 'none','YColor','none')

subplot(3,2,p2) 

hold on
shadedErrorBar(t,tc_green_avrg_stat_large{pre}(:,iCon),tc_green_se_stat_large{pre}(:,iCon),'--k');
hold on
shadedErrorBar(t,tc_green_avrg_stat_large{post}(:,iCon),tc_green_se_stat_large{post}(:,iCon),'--b','transparent');
hold on
if iCon==1
    title("Large pupil")
elseif iCon==3
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
end
set(gca,'XColor', 'none','YColor','none')
ylim([-.02 .25]);
sgtitle(['Pyr',' n = ', num2str(length(have_allPupil_green))])

x0=5;
y0=0;
width=4;
height=9;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')


end  


print(fullfile(fnout,'Fig_S6C.pdf'),'-dpdf');

%% pupil stats
pupilMeans_clean = squeeze(mean(pupilMeans_concat(:,:,:),1,'omitmissing')); %average the two days
pupilMeans_clean = pupilMeans_clean([2,1,3],:); %rearrange rows
pupilMeans_clean = pupilMeans_clean .*2;

%pupilMeans_clean = pupilMeans_clean./(pupilMeans_clean(1,:)); %normalize to small pupil for each mouse
pupilMean_overall =mean(pupilMeans_clean,2,'omitmissing');
pupilSTD_overall =std(pupilMeans_clean,[],2,'omitmissing');
pupilSE_overall = pupilSTD_overall./sqrt(nSess);   

figure
plot(pupilMeans_clean(:,:),'-','Color',[.5 .5 .5])
xlim([.75 3.25])
box off
set(gca,'TickDir','out')
xticks([1 2 3])
xticklabels({'Small pupil','Large pupil','Running'})
ylabel('Pupil diameter (mm)')
hold on
ylim([0.2 .8])
scatter([1 2 3],pupilMeans_clean(:,:),10, "MarkerEdgeColor","none","MarkerFaceColor",[0 0 0],"MarkerFaceAlpha",.25)
errorbar(pupilMean_overall,pupilSE_overall,'.-k','MarkerSize',15)


x0=5;
y0=5;
width=2;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(fullfile(fnout,'Fig_S5b.pdf'),'-dpdf');

% confirm that there's no difference in size across days within the "small"
% pupil trials
[h p] = ttest(pupilMeans_concat(pre,2,:),pupilMeans_concat(post,2,:))
% large pupil trials are significantly
[h p] = ttest(pupilMeans_concat(pre,1,:),pupilMeans_concat(post,1,:))



%confirm that pupil is significantly large in large trials, averaging over
%days
[h p] = ttest(pupilMeans_clean(1,:),pupilMeans_clean(2,:))
%test whether pupuil is different between large trials and running
[h p] = ttest(pupilMeans_clean(2,:),pupilMeans_clean(3,:))

%%
motorByPupil_clean = squeeze(mean(pupilMeans_concat(:,:,:),1,'omitmissing')); %average the two days
motorByPupil_clean = motorByPupil_clean([2,1],:); %rearrange rows


%pupilMeans_clean = pupilMeans_clean./(pupilMeans_clean(1,:)); %normalize to small pupil for each mouse
pupilMotorMean_overall =mean(motorByPupil_clean,2,'omitmissing');
pupilMotorSTD_overall =std(motorByPupil_clean,[],2,'omitmissing');
pupilMotorSE_overall = pupilMotorSTD_overall./sqrt(nSess);   

figure
plot(motorByPupil_clean(:,:),'-','Color',[.5 .5 .5])
xlim([.75 2.25])
box off
set(gca,'TickDir','out')
xticks([1 2])
xticklabels({'Small pupil','Large pupil'})
ylabel('Wheel speed (cm/s)')
hold on
ylim([.1 .4])
scatter([1 2],motorByPupil_clean(:,:),10, "MarkerEdgeColor","none","MarkerFaceColor",[0 0 0],"MarkerFaceAlpha",.25)
errorbar(pupilMotorMean_overall,pupilMotorSE_overall,'.-k','MarkerSize',15)



x0=5;
y0=5;
width=1;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(fullfile(fnout,'Fig_s5c.pdf'),'-dpdf');


%look at how motor activity differs with pupil size
%control day
[h p] = ttest(motorByPupil_clean(1,:),motorByPupil_clean(2,:))

%% contrast response pupil

conResp_green_avrg_small = cell(1,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_small = cell(1,nd); %same for red
conResp_green_se_small = cell(1,nd); %this will be the se across all green cells
conResp_red_se_small = cell(1,nd); %same for red

conResp_green_avrg_large = cell(1,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_large = cell(1,nd); %same for red
conResp_green_se_large = cell(1,nd); %this will be the se across all green cells
conResp_red_se_large = cell(1,nd); %same for red

have_allPupil_green = intersect(have_allPupil, green_ind_concat);
have_allPupil_red = intersect(have_allPupil, red_ind_concat);

for id = 1:nd
   
        
    conResp_green_avrg_small{id}=mean(pref_responses_stat_smallPupil_concat{id}(have_allPupil_green,:),1,'omitnan');
    green_std=std(pref_responses_stat_smallPupil_concat{id}(have_allPupil_green,:),1,'omitnan');
    conResp_green_se_small{id}=green_std/sqrt(length(have_allPupil_green));
    
    conResp_red_avrg_small{id}=mean(pref_responses_stat_smallPupil_concat{id}(have_allPupil_red,:),1,'omitnan');
    red_std=std(pref_responses_stat_smallPupil_concat{id}(have_allPupil_red,:),1,'omitnan');
    conResp_red_se_small{id}=red_std/sqrt(length(have_allPupil_red));
    
    clear green_std red_std
 
    conResp_green_avrg_large{id}=nanmean(pref_responses_stat_largePupil_concat{id}(have_allPupil_green ,:),1);
    green_std=nanstd(pref_responses_stat_largePupil_concat{id}(have_allPupil_green,:),1);
    conResp_green_se_large{id}=green_std/sqrt(length(have_allPupil_green));
    
    conResp_red_avrg_large{id}=nanmean(pref_responses_stat_largePupil_concat{id}(have_allPupil_red,:),1);
    red_std=nanstd(pref_responses_stat_largePupil_concat{id}(have_allPupil_red,:),1);
    conResp_red_se_large{id}=red_std/sqrt(length(have_allPupil_red));
    
    clear green_std red_std
end



figure
subplot(1,2,1) %for the first day
errorbar(cons,conResp_red_avrg_small{pre},conResp_red_se_small{pre},'k');
hold on
errorbar(cons,conResp_red_avrg_small{post},conResp_red_se_small{post},'b');
title('Small pupil')
ylabel('dF/F') 
xlabel('contrast') 
set(gca, 'TickDir', 'out')
ylim([0 .15])
xlim([0 1.2])
xticks([.25 .50 1.00])
xticklabels([25 50 100])
%xlabel('contrast') 
set(gca, 'TickDir', 'out')
box off


subplot(1,2,2) %for the first day
errorbar(cons,conResp_red_avrg_large{pre},conResp_red_se_large{pre},'k');
hold on
errorbar(cons,conResp_red_avrg_large{post},conResp_red_se_large{post},'b');
title('Large pupil')
xlabel('contrast') 
ylim([0 .15])
xlim([0 1.2])
xticks([.25 .50 1.00])
xticklabels([25 50 100])
%xlabel('contrast') 
set(gca, 'TickDir', 'out')
box off

x0=5;
y0=5;
width=2.8;
height=1.45;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(fullfile(fnout,'Fig_5C.pdf'),'-dpdf');

figure
subplot(1,2,1) %for the first day
errorbar(cons,conResp_green_avrg_small{pre},conResp_green_se_small{pre},'--k');
hold on
errorbar(cons,conResp_green_avrg_small{post},conResp_green_se_small{post},'--b');
%title('Small pupil')
ylabel('dF/F') 
xlabel('contrast') 
set(gca, 'TickDir', 'out')
ylim([0 .3])
xlim([0 1.2])
xticks([.25 .50 1.00])
xticklabels([25 50 100])
%xlabel('contrast') 
set(gca, 'TickDir', 'out')
box off


subplot(1,2,2) %for the first day
errorbar(cons,conResp_green_avrg_large{pre},conResp_green_se_large{pre},'--k');
hold on
errorbar(cons,conResp_green_avrg_large{post},conResp_green_se_large{post},'--b');
%title('Large pupil')
xlabel('contrast') 
ylim([0 .3])
xlim([0 1.2])
xticks([.25 .50 1.00])
xticklabels([25 50 100])
%xlabel('contrast') 
set(gca, 'TickDir', 'out')
box off

x0=5;
y0=5;
width=2.8;
height=1.45;
set(gcf,'units','inches','position',[x0,y0,width,height])


print(fullfile(fnout,'Fig_S6D.pdf'),'-dpdf');

%% stats for data split by pupil

% prepare data for ANOVA
dfof_small = horzcat(pref_responses_stat_smallPupil_concat{pre},pref_responses_stat_smallPupil_concat{post});
dfof_large = horzcat(pref_responses_stat_largePupil_concat{pre},pref_responses_stat_largePupil_concat{post});

cellID=(1:nKeep_total)';
cell_type_col=categorical(red_concat)';

dfof_full_pupil = horzcat(dfof_small,dfof_large);
dfof_table_pupil = array2table(dfof_full_pupil,'VariableNames',{'Sd1c1','Sd1c2','Sd1c3', ...
    'Sd2c1','Sd2c2','Sd2c3','Ld1c1','Ld1c2','Ld1c3','Ld2c1','Ld2c2','Ld2c3'});

labels_table =table(mouseID,cellID,cell_type_col,'VariableNames',{'mouseID' 'cellID' 'cellType'});
dfof_summary_full_pupil = [labels_table,dfof_table_pupil];
matched_dfof_summary_pupil = dfof_summary_full_pupil(ismember(dfof_summary_full_pupil.cellID,have_allPupil),:);


SST_matched_dfof_pupil = matched_dfof_summary_pupil((matched_dfof_summary_pupil.cellType=='1'),:);
Pyr_matched_dfof_pupil = matched_dfof_summary_pupil((matched_dfof_summary_pupil.cellType=='0'),:);
clear dfof_stat_table cell_type_col cellID dfof_stat_matched

% run an ANOVA on the full model for each cell type
DART = categorical([1 1 1 2 2 2 1 1 1 2 2 2].');
contrast = categorical([1 2 3 1 2 3 1 2 3 1 2 3].');
pupilSize = categorical([1 1 1 1 1 1 2 2 2 2 2 2].');

% prepare new categorical variables
w = table(DART, contrast, pupilSize, ...
    'VariableNames', {'DART', 'contrast', 'pupilSize'});% within-design

rm_SST_matched = fitrm(SST_matched_dfof_pupil, 'Sd1c1-Ld2c3 ~ 1', 'WithinDesign', w);
ranova(rm_SST_matched, 'withinmodel', 'pupilSize*DART*contrast')

rm_Pyr_matched = fitrm(Pyr_matched_dfof_pupil, 'Sd1c1-Ld2c3 ~ 1', 'WithinDesign', w);
ranova(rm_Pyr_matched, 'withinmodel', 'pupilSize*DART*contrast')

% models for small and large pupil seperately 
w = table(categorical([1 1 1 2 2 2 ].'), categorical([1 2 3 1 2 3].'), 'VariableNames', {'DART', 'contrast'}); % within-design

rm_SST_small = fitrm(SST_matched_dfof_pupil, 'Sd1c1-Sd2c3 ~ 1', 'WithinDesign', w);
ranova(rm_SST_small, 'withinmodel', 'DART*contrast')

rm_SST_large = fitrm(SST_matched_dfof_pupil, 'Ld1c1-Ld2c3 ~ 1', 'WithinDesign', w);
ranova(rm_SST_large, 'withinmodel', 'DART*contrast')


rm_Pyr_small = fitrm(Pyr_matched_dfof_pupil, 'Sd1c1-Sd2c3 ~ 1', 'WithinDesign', w);
ranova(rm_Pyr_small, 'withinmodel', 'DART*contrast')

rm_Pyr_large = fitrm(Pyr_matched_dfof_pupil, 'Ld1c1-Ld2c3 ~ 1', 'WithinDesign', w);
ranova(rm_Pyr_large, 'withinmodel', 'DART*contrast')


% pairwise ttests for dfof response at each contrast for SST cells
[sst_h1, sst_p1]= ttest(pref_responses_stat_smallPupil_concat{pre}(have_allPupil_red,1),pref_responses_stat_smallPupil_concat{post}(have_allPupil_red,1));
[sst_h2, sst_p2]= ttest(pref_responses_stat_smallPupil_concat{pre}(have_allPupil_red,2),pref_responses_stat_smallPupil_concat{post}(have_allPupil_red,2));
[sst_h3, sst_p3]= ttest(pref_responses_stat_smallPupil_concat{pre}(have_allPupil_red,3),pref_responses_stat_smallPupil_concat{post}(have_allPupil_red,3));

%corrected for three tests
sst_pvalues_small = [(sst_p1*3);(sst_p2*3);(sst_p3*3)];
contrasts = cons';
table(contrasts,sst_pvalues_small)

% pairwise ttests for dfof response at each contrast for SST cells
[sst_h1, sst_p1]= ttest(pref_responses_stat_largePupil_concat{pre}(have_allPupil_red,1),pref_responses_stat_largePupil_concat{post}(have_allPupil_red,1));
[sst_h2, sst_p2]= ttest(pref_responses_stat_largePupil_concat{pre}(have_allPupil_red,2),pref_responses_stat_largePupil_concat{post}(have_allPupil_red,2));
[sst_h3, sst_p3]= ttest(pref_responses_stat_largePupil_concat{pre}(have_allPupil_red,3),pref_responses_stat_largePupil_concat{post}(have_allPupil_red,3));

%corrected for three tests
sst_pvalues_large = [(sst_p1*3);(sst_p2*3);(sst_p3*3)];
contrasts = cons';
table(contrasts,sst_pvalues_large)

% pairwise ttests for dfof response at each contrast for SST cells
[pyr_h1, pyr_p1]= ttest(pref_responses_stat_smallPupil_concat{pre}(have_allPupil_green,1),pref_responses_stat_smallPupil_concat{post}(have_allPupil_green,1));
[pyr_h2, pyr_p2]= ttest(pref_responses_stat_smallPupil_concat{pre}(have_allPupil_green,2),pref_responses_stat_smallPupil_concat{post}(have_allPupil_green,2));
[pyr_h3, pyr_p3]= ttest(pref_responses_stat_smallPupil_concat{pre}(have_allPupil_green,3),pref_responses_stat_smallPupil_concat{post}(have_allPupil_green,3));

%corrected for three tests
pyr_pvalues_small = [(pyr_p1*3);(pyr_p2*3);(pyr_p3*3)];
contrasts = cons';
table(contrasts,pyr_pvalues_small)

% pairwise ttests for dfof response at each contrast for SST cells
[pyr_h1, pyr_p1]= ttest(pref_responses_stat_largePupil_concat{pre}(have_allPupil_green,1),pref_responses_stat_largePupil_concat{post}(have_allPupil_green,1));
[pyr_h2, pyr_p2]= ttest(pref_responses_stat_largePupil_concat{pre}(have_allPupil_green,2),pref_responses_stat_largePupil_concat{post}(have_allPupil_green,2));
[pyr_h3, pyr_p3]= ttest(pref_responses_stat_largePupil_concat{pre}(have_allPupil_green,3),pref_responses_stat_largePupil_concat{post}(have_allPupil_green,3));

%corrected for three tests
pyr_pvalues_large = [(pyr_p1*3);(pyr_p2*3);(pyr_p3*3)];
contrasts = cons';
table(contrasts,pyr_pvalues_large)

%calculate norm_diff
norm_diff_pupil = nan(2,nCon,nKeep_total);
for i = 1:nKeep_total
    for iCon = 1:nCon
        %for smallPupil trials
        mean_pre_smallPupil = mean(pref_allTrials_smallPupil_concat{iCon,pre}{i});
        mean_post_smallPupil=mean(pref_allTrials_smallPupil_concat{iCon,post}{i});
        std_pre_smallPupil = std(pref_allTrials_smallPupil_concat{iCon,pre}{i});
        norm_diff_smallPupil = (mean_post_smallPupil-mean_pre_smallPupil) / std_pre_smallPupil;

        %for running trials
        mean_pre_largePupil = mean(pref_allTrials_largePupil_concat{iCon,pre}{i});
        mean_post_largePupil=mean(pref_allTrials_largePupil_concat{iCon,post}{i});
        std_pre_largePupil = std(pref_allTrials_largePupil_concat{iCon,pre}{i});
        norm_diff_largePupil = (mean_post_largePupil-mean_pre_largePupil)/ std_pre_largePupil;

        %putting data into matrix
        norm_diff_pupil(1,iCon,i)=norm_diff_smallPupil; %first is smallPupilionary
        norm_diff_pupil(2,iCon,i)=norm_diff_largePupil; %second is running
clear mean_pre_smallPupil mean_post_smallPupil std_pre_smallPupil mean_pre_largePupil mean_post_largePupil std_pre_largePupil norn_diff_smallPupil norm_diff_largePupil
    end 
end
%remove any infiinty values resulting from divisions by zero, and turn
%those into NANs instead
norm_diff_pupil(find(norm_diff_pupil == -Inf))=NaN;
norm_diff_pupil(find(norm_diff_pupil == Inf))=NaN;

%%

%make a subset of normalized difference for the SST cells only, then make
% find how many are facilitated or suppressed by more than 1 std from
% baseline
%reselect norm_diff)red to be only the "red all" subset
norm_diff_red = norm_diff_pupil(:,:,red_all);

facil_red=norm_diff_red(:,:,:)>=1;
supp_red=norm_diff_red(:,:,:)<=-1;

N=length(red_all);
facil_table_smallPupil = sum(facil_red(1,:,:),3)/N;
supp_table_smallPupil = sum(supp_red(1,:,:),3)/N;
facil_table_largePupil = sum(facil_red(2,:,:),3)/N;
supp_table_largePupil = sum(supp_red(2,:,:),3)/N;

figure;
subplot(1,2,1)
b=bar([1,2,3],[supp_table_smallPupil; supp_table_largePupil],'grouped','FaceColor',"#00ffff",'EdgeColor', [1 1 1]);
b(1).FaceColor="#70D0F6"
b(2).FaceColor="#0C8ABB"
xticklabels({'25','50','100'})
ylim([0 .4])
title('Suppressed')
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

subplot(1,2,2)
b=bar([1,2,3],[facil_table_smallPupil; facil_table_largePupil],'FaceColor',"#a329cc",'EdgeColor', [1 1 1]);
b(1).FaceColor="#C983B1"
b(2).FaceColor="#883367"
xticklabels({'25','50','100'})
ylim([0 .4])
title('Facilitated')
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'pupil_normDiff.pdf'),'-dpdf')

%%

%compute chi squares for suppression, small pupil vs large pupil
%25% contrast

n1=supp_table_smallPupil(1)*N;
n2=supp_table_largePupil(1)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil1,p1] = crosstab(x1,x2);

%50
n1=supp_table_smallPupil(2)*N;
n2=supp_table_largePupil(2)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil2,p2] = crosstab(x1,x2);

%100
n1=supp_table_smallPupil(3)*N;
n2=supp_table_largePupil(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil3,p3] = crosstab(x1,x2);



[chi2smallPupil1, chi2smallPupil2, chi2smallPupil3; p1*3, p2*3,p3*3]
%[chi2smallPupil1, chi2smallPupil2, chi2smallPupil3; p1, p2,p3]

clear h p1 p2 p3 chi2smallPupil1 chi2smallPupil2 chi2smallPupil3

%compute chi squares for facilitation
%25
n1=facil_table_smallPupil(1)*N;
n2=facil_table_largePupil(1)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil1,p1] = crosstab(x1,x2);

%50
n1=facil_table_smallPupil(2)*N;
n2=facil_table_largePupil(2)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil2,p2] = crosstab(x1,x2);

%100
n1=facil_table_smallPupil(3)*N;
n2=facil_table_largePupil(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil3,p3] = crosstab(x1,x2);

[chi2smallPupil1, chi2smallPupil2, chi2smallPupil3; p1, p2,p3]
clear h p1 p2 p3 chi2smallPupil1 chi2smallPupil2 chi2smallPupil3 n1 n2 x1 x2

%% For pyramidal cells

norm_diff_green = norm_diff_pupil(:,:,green_all);

facil_green=norm_diff_green(:,:,:)>=1;
supp_green=norm_diff_green(:,:,:)<=-1;

N=length(green_all);
facil_table_smallPupil_green = sum(facil_green(1,:,:),3)/N;
supp_table_smallPupil_green = sum(supp_green(1,:,:),3)/N;
facil_table_largePupil_green = sum(facil_green(2,:,:),3)/N;
supp_table_largePupil_green = sum(supp_green(2,:,:),3)/N;

figure;
subplot(1,2,1)
b=bar([1,2,3],[supp_table_smallPupil_green; supp_table_largePupil_green],'grouped','FaceColor',"#00ffff",'EdgeColor', [1 1 1]);
b(1).FaceColor="#70D0F6"
b(2).FaceColor="#0C8ABB"
xticklabels({'25','50','100'})
ylim([0 .4])
title('Suppressed')
ylabel(["Fraction Pyr cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

subplot(1,2,2)
b=bar([1,2,3],[facil_table_smallPupil_green; facil_table_largePupil_green],'FaceColor',"#a329cc",'EdgeColor', [1 1 1]);
b(1).FaceColor="#C983B1"
b(2).FaceColor="#883367"
xticklabels({'25','50','100'})
ylim([0 .4])
title('Facilitated')
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'pupil_normDiff_green.pdf'),'-dpdf')
%% compute chi squares for suppression, small pupil vs large pupil
%25% contrast

n1=supp_table_smallPupil_green(1)*N;
n2=supp_table_largePupil_green(1)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil1,p1] = crosstab(x1,x2);

%50
n1=supp_table_smallPupil_green(2)*N;
n2=supp_table_largePupil_green(2)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil2,p2] = crosstab(x1,x2);

%100
n1=supp_table_smallPupil_green(3)*N;
n2=round(supp_table_largePupil_green(3)*N);
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil3,p3] = crosstab(x1,x2);



[chi2smallPupil1, chi2smallPupil2, chi2smallPupil3; p1*3, p2*3,p3*3]
%[chi2smallPupil1, chi2smallPupil2, chi2smallPupil3; p1, p2,p3]

clear h p1 p2 p3 chi2smallPupil1 chi2smallPupil2 chi2smallPupil3



%compute chi squares for facilitation
%25
n1=facil_table_smallPupil_green(1)*N;
n2=facil_table_largePupil_green(1)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil1,p1] = crosstab(x1,x2);

%50
n1=facil_table_smallPupil_green(2)*N;
n2=facil_table_largePupil_green(2)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil2,p2] = crosstab(x1,x2);

%100
n1=facil_table_smallPupil_green(3)*N;
n2=facil_table_largePupil_green(3)*N;
x1 = [repmat('a',N,1); repmat('b',N,1)];
x2 = [repmat(1,n1,1); repmat(2,N-n1,1); repmat(1,n2,1); repmat(2,N-n2,1)];
[tbl,chi2smallPupil3,p3] = crosstab(x1,x2);

[chi2smallPupil1, chi2smallPupil2, chi2smallPupil3; p1, p2,p3]
clear h p1 p2 p3 chi2smallPupil1 chi2smallPupil2 chi2smallPupil3 n1 n2 x1 x2


%% get a table of capture values
capture = getCaptureValues_annulus(mice);
table(mice,capture(3,:)')
edges = linspace(1, 2, 10); % Create 20 bins.
% Plot the histogram.
histogram(capture(3,:),'BinEdges',edges);
xlim([1 2])
box off
set(gca, 'TickDir', 'out')
x0=5;
y0=5;
width=1.1;
height=1.1;
set(gcf,'units','inches','position',[x0,y0,width,height])

%%
captureByCell = [];
for iSess = 1:nSess
    temp = repmat(capture(3,iSess),[nKeep_concat(iSess),1]);
    captureByCell=[captureByCell;temp];
end

%% timecourses seperated by capture


%make figure with se shaded, one figure per contrast - stationary

tc_green_avrg_stat = cell(1,nd); 
tc_red_avrg_stat = cell(1,nd); 
tc_green_se_stat = cell(1,nd); 
tc_red_se_stat = cell(1,nd); 

tc_green_avrg_stat_b = cell(1,nd); 
tc_red_avrg_stat_b = cell(1,nd); 
tc_green_se_stat_b = cell(1,nd); 
tc_red_se_stat_b = cell(1,nd);


for id = 1:nd
    for iCon=1:nCon
        
    tc_green_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})<cutOff),iCon),2);
    green_std=nanstd(tc_trial_avrg_stat_concat{id}(:,haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})<cutOff),iCon),[],2);
    tc_green_se_stat{id}(:,iCon)=green_std/sqrt(length(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})<cutOff)));
    
    tc_red_avrg_stat{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})<cutOff),iCon),2);
    red_std=nanstd(tc_trial_avrg_stat_concat{id}(:,haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})<cutOff),iCon),[],2);
    tc_red_se_stat{id}(:,iCon)=red_std/sqrt(length(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})<cutOff)));


    tc_green_avrg_stat_b{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})>=cutOff),iCon),2);
    green_std=nanstd(tc_trial_avrg_stat_concat{id}(:,haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})>=cutOff),iCon),[],2);
    tc_green_se_stat_b{id}(:,iCon)=green_std/sqrt(length(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})>=cutOff)));
    
    tc_red_avrg_stat_b{id}(:,iCon)=nanmean(tc_trial_avrg_stat_concat{id}(:,haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})>=cutOff),iCon),2);
    red_std=nanstd(tc_trial_avrg_stat_concat{id}(:,haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})>=cutOff),iCon),[],2);
    tc_red_se_stat_b{id}(:,iCon)=red_std/sqrt(length(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})>=cutOff)));

    clear green_std red_std
    end
end
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(2,2,1) %for the first day

ylim([-.05 .3]);
hold on
shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon),tc_green_se_stat{pre}(:,iCon),'--k');
hold on
shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon),tc_green_se_stat{post}(:,iCon),'--b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['low capture -HTP',' n = ', num2str(length(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})<cutOff)))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(2,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_stat{pre}(:,iCon),tc_red_se_stat{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_stat{post}(:,iCon),tc_red_se_stat{post}(:,iCon),'b');
ylim([-.05 .3]);
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['low capture +HTP',' n = ', num2str(length(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})<cutOff)))])
set(gca,'XColor', 'none','YColor','none')

subplot(2,2,3) %for the first day

ylim([-.05 .3]);
hold on
shadedErrorBar(t,tc_green_avrg_stat_b{pre}(:,iCon),tc_green_se_stat_b{pre}(:,iCon),'--k');
hold on
shadedErrorBar(t,tc_green_avrg_stat_b{post}(:,iCon),tc_green_se_stat_b{post}(:,iCon),'--b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['high capture -HTP',' n = ', num2str(length(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})>=cutOff)))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(2,2,4) %for the second day
shadedErrorBar(t,tc_red_avrg_stat_b{pre}(:,iCon),tc_red_se_stat_b{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_stat_b{post}(:,iCon),tc_red_se_stat_b{post}(:,iCon),'b');
ylim([-.05 .3]);
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['high capture +HTP',' n = ', num2str(length(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})>=cutOff)))])

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['stationary, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) '_stat_capt_timecourses.pdf']),'-dpdf');
end 

clear tc_green_avrg_stat tc_red_avrg_stat tc_green_se_stat tc_red_se_stat tc_green_avrg_stat_b tc_red_avrg_stat_b tc_green_se_stat_b tc_red_se_stat_b 

% make figure with se shaded, one figure per contrast - loc

tc_green_avrg_loc = cell(1,nd); 
tc_red_avrg_loc = cell(1,nd); 
tc_green_se_loc = cell(1,nd); 
tc_red_se_loc = cell(1,nd); 

tc_green_avrg_loc_b = cell(1,nd); 
tc_red_avrg_loc_b = cell(1,nd); 
tc_green_se_loc_b = cell(1,nd); 
tc_red_se_loc_b = cell(1,nd); 


for id = 1:nd
    for iCon=1:nCon
        
    tc_green_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})<cutOff),iCon),2);
    green_std=nanstd(tc_trial_avrg_loc_concat{id}(:,haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})<cutOff),iCon),[],2);
    tc_green_se_loc{id}(:,iCon)=green_std/sqrt(length(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})<cutOff)));
    
    tc_red_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})<cutOff),iCon),2);
    red_std=nanstd(tc_trial_avrg_loc_concat{id}(:,haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})<cutOff),iCon),[],2);
    tc_red_se_loc{id}(:,iCon)=red_std/sqrt(length(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})<cutOff)));


    tc_green_avrg_loc_b{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})>=cutOff),iCon),2);
    green_std=nanstd(tc_trial_avrg_loc_concat{id}(:,haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})>=cutOff),iCon),[],2);
    tc_green_se_loc_b{id}(:,iCon)=green_std/sqrt(length(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})>=cutOff)));
    
    tc_red_avrg_loc_b{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})>=cutOff),iCon),2);
    red_std=nanstd(tc_trial_avrg_loc_concat{id}(:,haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})>=cutOff),iCon),[],2);
    tc_red_se_loc_b{id}(:,iCon)=red_std/sqrt(length(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})>=cutOff)));
    



    clear green_std red_std
    end
end
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_loc{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);



for iCon = 1:nCon
figure
subplot(2,2,1) %for the first day

ylim([-.05 .35]);
hold on
shadedErrorBar(t,tc_green_avrg_loc{pre}(:,iCon),tc_green_se_loc{pre}(:,iCon),'--k');
hold on
shadedErrorBar(t,tc_green_avrg_loc{post}(:,iCon),tc_green_se_loc{post}(:,iCon),'--b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['low capture -HTP',' n = ', num2str(length(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})<cutOff)))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(2,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_loc{pre}(:,iCon),tc_red_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_loc{post}(:,iCon),tc_red_se_loc{post}(:,iCon),'b');
ylim([-.05 .35]);
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['low capture +HTP',' n = ', num2str(length(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})<cutOff)))])
set(gca,'XColor', 'none','YColor','none')

subplot(2,2,3) %for the first day

ylim([-.05 .35]);
hold on
shadedErrorBar(t,tc_green_avrg_loc_b{pre}(:,iCon),tc_green_se_loc_b{pre}(:,iCon),'--k');
hold on
shadedErrorBar(t,tc_green_avrg_loc_b{post}(:,iCon),tc_green_se_loc_b{post}(:,iCon),'--b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title(['high capture -HTP',' n = ', num2str(length(haveRunning_green{iCon}(captureByCell(haveRunning_green{iCon})>=cutOff)))])

ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')


subplot(2,2,4) %for the second day
shadedErrorBar(t,tc_red_avrg_loc_b{pre}(:,iCon),tc_red_se_loc_b{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_loc_b{post}(:,iCon),tc_red_se_loc_b{post}(:,iCon),'b');
ylim([-.05 .35]);
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
ylabel('dF/F') 
xlabel('s') 
title(['high capture +HTP',' n = ', num2str(length(haveRunning_red{iCon}(captureByCell(haveRunning_red{iCon})>=cutOff)))])

x0=5;
y0=5;
width=4;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')

sgtitle(['running, contrast = ' num2str(cons(iCon))])


print(fullfile(fnout,[num2str(cons(iCon)) '_loc_capt_timecourses.pdf']),'-dpdf');
end 


clear tc_green_avrg_loc tc_red_avrg_loc tc_green_se_loc tc_red_se_loc tc_green_avrg_loc_b tc_red_avrg_loc_b tc_green_se_loc_b tc_red_se_loc_b 

%% additional correlation analyses

%does the average R value differ between behavioral states?
R_by_state = NaN(4,nd,nKeep_total);

for id = 1:nd
    for i = 1:4
        R_by_state(i,id,:)=noiseCorr_concat{i,id}(1,:);
    end
end


figure;
subplot(1,2,1)
boxplot([squeeze(R_by_state(1,pre,red_all)),squeeze(R_by_state(4,pre,red_all))]);
hold on
scatter([1, 2],squeeze(R_by_state([1,4],pre,red_all)),20,'k', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)
xticklabels({'Stationary','Running'})
%ylim([-.3 1.3])
ylabel(["Mean R value"]) 
set(gca,'TickDir','out')        
box off
hold on
subplot(1,2,2)
boxplot([squeeze(R_by_state(1,post,red_all)),squeeze(R_by_state(4,post,red_all))]);
hold on
scatter([1, 2],squeeze(R_by_state([1,4],post,red_all)),20,'k', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)
xticklabels({'Stationary','Running'})
%ylim([-.3 1.3])
ylabel(["Mean R value"]) 
set(gca,'TickDir','out')        
box off

[h1,p1] =ttest(squeeze(R_by_state(1,pre,red_all)),squeeze(R_by_state(4,pre,red_all)))

%%
%does the average R value differ between contrasts within a behvaioral
%state?
R_by_contrast = NaN(nCon,nd,nKeep_total);

for id = 1:nd
    for iCon = 1:nCon
        R_by_contrast(iCon,id,:)=noiseCorrContrast_concat{1,iCon,id}(1,:);
    end
end


figure;
scatter([1, 2, 3],squeeze(R_by_contrast(:,pre,red_ind_concat)),20,'k', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)
hold on
boxplot([squeeze(R_by_contrast(1,pre,red_ind_concat)),squeeze(R_by_contrast(2,pre,red_ind_concat)),squeeze(R_by_contrast(3,pre,red_ind_concat))]);
xticklabels({'25%','50%','100%'})
ylim([-1 1])
ylabel(["Mean R value"]) 
set(gca,'TickDir','out')
box off


[h1,p1] =ttest(squeeze(R_by_contrast(1,pre,red_ind_concat)),squeeze(R_by_contrast(2,pre,red_ind_concat)));
[h2,p2] =ttest(squeeze(R_by_contrast(1,pre,red_ind_concat)),squeeze(R_by_contrast(3,pre,red_ind_concat)));
[h3,p3] =ttest(squeeze(R_by_contrast(2,pre,red_ind_concat)),squeeze(R_by_contrast(3,pre,red_ind_concat)));
[p1*3,p2*3,p3*3]
%%
% ANOVA for noise corr change across contrasts within stationary

data = squeeze(R_by_contrast(:,pre,red_ind_concat))';
dataTable = array2table(data,'VariableNames',{'C1','C2','C3'});
w = table(categorical([1 2 3 ].'), 'VariableNames', {'contrast'}); % within-design
rm_SST_stat = fitrm(dataTable, 'C1-C3 ~ 1', 'WithinDesign', w)
ranova(rm_SST_stat, 'withinmodel', 'contrast')
%% normDiff comparisons
%compare norm_diff across contrasts and states
        
figure;
subplot(2,2,1)
scatter(squeeze(norm_diff(1,1,green_ind_concat)),squeeze(norm_diff(1,3,green_ind_concat)))
xlabel('norm diff 25% contrast')
ylabel('norm diff 100% contrast')
xlim([-5 7.5])
ylim([-5 7.5])
title('Pyr')

subplot(2,2,2)
scatter(squeeze(norm_diff(1,1,green_ind_concat)),squeeze(norm_diff(2,1,green_ind_concat)))
xlabel('norm diff 25% contrast stationary')
ylabel('norm diff 25% contrast running')
xlim([-5 7.5])
ylim([-5 7.5])
title('Pyr')

subplot(2,2,3)
scatter(squeeze(norm_diff(1,1,red_ind_concat)),squeeze(norm_diff(1,3,red_ind_concat)))
xlabel('norm diff 25% contrast')
ylabel('norm diff 100% contrast')
xlim([-5 7.5])
ylim([-5 7.5])
title('SST')

subplot(2,2,4)
scatter(squeeze(norm_diff(1,1,red_ind_concat)),squeeze(norm_diff(2,1,red_ind_concat)))
xlabel('norm diff 25% contrast stationary')
ylabel('norm diff 25% contrast running')
xlim([-5 7.5])
ylim([-5 7.5])
title('SST')

%% comparing normalized difference to noise correlation

mySample = red_ind_concat;
% Extract data (stationary state only)
stationary_norm_diff = squeeze(norm_diff(1,:,mySample)); % [contrast, neuron]
stationary_bsln_std = squeeze(bsln_std(1,:,mySample));   % [contrast, neuron]
baseline_noise_corr = noiseCorr_OG_concat{2}(1,mySample); % [1, neuron]

% Number of contrasts and neurons
[num_contrasts, num_neurons] = size(stationary_norm_diff);

% Store results
slopes = zeros(num_contrasts, 1);
intercepts = zeros(num_contrasts, 1);
r_squared = zeros(num_contrasts, 1);
p_values = zeros(num_contrasts, 1);
std_errors = zeros(num_contrasts, 1);

% Create figure
figure('Position', [100, 100, 1200, 400]);

% Analyze each contrast
for contrast = 1:num_contrasts
    % Extract data for this contrast
    y = stationary_norm_diff(contrast, :)';  % neural response difference
    x = baseline_noise_corr';                 % noise correlations
    
    % Handle NaNs and missing data
    valid_idx = ~isnan(x) & ~isnan(y) & isfinite(x) & isfinite(y);
    
    % Check if we have enough valid data points
    if sum(valid_idx) < 3
        warning('Not enough valid data points for contrast %d. Skipping.', contrast);
        title(sprintf('Contrast %d: Insufficient Data', contrast));
        continue;
    end
    
    % Use only valid data
    x_valid = x(valid_idx);
    y_valid = y(valid_idx);
    
    % Perform unweighted linear regression
    [p, stats] = polyfit(x_valid, y_valid, 1);
    
    % Store results
    intercepts(contrast) = p(2);
    slopes(contrast) = p(1);
    
    % Calculate standard error of the slope
    yfit = polyval(p, x_valid);
    residuals = y_valid - yfit;
    SSE = sum(residuals.^2);
    n = length(x_valid);
    
    % Standard error of the regression
    SE_regression = sqrt(SSE / (n-2));
    
    % Sum of squares of x deviations
    SS_x = sum((x_valid - mean(x_valid)).^2);
    
    % Standard error of the slope
    std_errors(contrast) = SE_regression / sqrt(SS_x);
    
    % Calculate R-squared
    SS_total = sum((y_valid - mean(y_valid)).^2);
    r_squared(contrast) = 1 - SSE/SS_total;
    
    % Calculate p-value for slope
    t_stat = slopes(contrast) / std_errors(contrast);
    p_values(contrast) = 2 * (1 - tcdf(abs(t_stat), length(x_valid) - 2));
    
    % Plot regression
    subplot(1, 3, contrast);
    
    % Create scatter plot with uniform size
    scatter(x_valid, y_valid, 50, 'filled', 'MarkerFaceAlpha', 0.7);
    hold on;
    
    % Add regression line
    x_range = linspace(min(x_valid), max(x_valid), 100);
    y_line = p(1) * x_range + p(2);
    plot(x_range, y_line, 'r-', 'LineWidth', 2);
    
    % Labels and title
    contrasts = [25, 50, 100];
    title(sprintf('Contrast: %d%%', contrasts(contrast)), 'FontSize', 12);
    xlabel('Noise Correlation (Baseline)', 'FontSize', 11);
    ylabel('Normalized Difference', 'FontSize', 11);
    
    xlim([-.2 1])
    ylim([-8 8])

    % Add horizontal line at y=0
    hline(0)
    
    % Add regression equation and stats
    % Calculate text position dynamically
    x_range = max(x_valid) - min(x_valid);
    y_range = range(ylim);
    text_x = min(x_valid);
    text_y_start = max(ylim) - 0.2 * y_range;
    text_y_step = 0.08 * y_range;
    
    text(text_x, text_y_start, sprintf('y = %.3fx + %.3f', slopes(contrast), intercepts(contrast)), 'FontSize', 10);
    text(text_x, text_y_start - text_y_step, sprintf('R = %.3f', r_squared(contrast)), 'FontSize', 10);
    text(text_x, text_y_start - 2*text_y_step, sprintf('p = %.4f', p_values(contrast)), 'FontSize', 10);
    
    % Add sample size
    text(text_x, text_y_start - 3*text_y_step, sprintf('n = %d', sum(valid_idx)), 'FontSize', 10);
    
    % Set box style and tick direction
    box off;
    set(gca, 'TickDir', 'out');
end

% Adjust spacing
sgtitle('Unweighted Regression: Noise Correlation vs. Normalized Difference (Stationary State)', 'FontSize', 14);
set(gcf, 'Color', 'w');

% Create a results table
contrast_labels = {'25%', '50%', '100%'};
results_table = table(contrast_labels', slopes, std_errors, intercepts, r_squared, p_values, ...
    'VariableNames', {'Contrast', 'Slope', 'StdError', 'Intercept', 'RSquared', 'PValue'});

% Display results table
disp('Unweighted Regression Results:');
disp(results_table);

% Save figure
saveas(gcf, 'unweighted_regression_plot.png');
fprintf('Figure saved as unweighted_regression_plot.png\n');

% Optional: Save results to CSV
writetable(results_table, 'unweighted_regression_results.csv');
fprintf('Results saved to unweighted_regression_results.csv\n');

% Additional analysis: Test if slopes are significantly different across contrasts
fprintf('\nComparing slopes across contrast levels:\n');

% Pairwise comparison of slopes (z-test)
for i = 1:num_contrasts
    for j = i+1:num_contrasts
        % Calculate Z-statistic for difference between slopes
        z_diff = (slopes(i) - slopes(j)) / sqrt(std_errors(i)^2 + std_errors(j)^2);
        p_diff = 2 * (1 - normcdf(abs(z_diff)));
        fprintf('Contrast %s vs %s: Difference = %.3f, z = %.3f, p = %.4f\n', ...
            contrast_labels{i}, contrast_labels{j}, slopes(i) - slopes(j), z_diff, p_diff);
    end
end

% Additional figure: Compare slopes across contrasts
figure('Position', [100, 500, 500, 400]);
bar(1:num_contrasts, slopes);
ylim([-2,2])
hold on;

% Add error bars
errorbar(1:num_contrasts, slopes, std_errors, 'k.', 'LineWidth', 1.5);

% Add labels and title
xlabel('Contrast', 'FontSize', 12);
ylabel('Regression Slope', 'FontSize', 12);
title('Comparison of Regression Slopes Across Contrasts', 'FontSize', 14);
xticks(1:num_contrasts);
xticklabels(contrast_labels);
grid on;
box off;
set(gca, 'TickDir', 'out');

% Save comparison figure
saveas(gcf, 'slope_comparison.png');
fprintf('Slope comparison figure saved as slope_comparison.png\n');
%% Plot norm_diff by contrast for high and low noiseCorr cells
% Ask user if they want to remove outliers
remove_outliers = input('Do you want to remove outliers? (1 for yes, 0 for no): ');

% If user chooses to remove outliers, ask for standard deviation threshold
std_threshold = 2; % Default value
if remove_outliers
    std_threshold = input('Enter the number of standard deviations to use for outlier detection (default is 2): ');
    if isempty(std_threshold)
        std_threshold = 2;
    end
end

% Create filtered copies of the data
high_filtered = norm_diff(1, :, redHigh);
low_filtered = norm_diff(1, :, redLow);

% Initialize counters for outliers
high_outliers_count = 0;
low_outliers_count = 0;

% For each contrast level, optionally remove outliers
for contrast = 1:size(norm_diff, 2)
    % For high correlation neurons
    high_contrast_data = squeeze(high_filtered(1, contrast, :));
    
    if remove_outliers
        high_mean = nanmean(high_contrast_data);
        high_std = nanstd(high_contrast_data);

        % Find outliers and set them to NaN
        outlier_indices = abs(high_contrast_data - high_mean) > std_threshold * high_std;
        high_contrast_outliers = sum(outlier_indices);
        high_outliers_count = high_outliers_count + high_contrast_outliers;

        % Print number of outliers for this contrast level
        fprintf('Contrast %d, High correlation neurons: %d outliers removed\n', contrast, high_contrast_outliers);

        high_contrast_data(outlier_indices) = NaN;
    end

    % Update the filtered data
    high_filtered(1, contrast, :) = high_contrast_data;

    % For low correlation neurons
    low_contrast_data = squeeze(low_filtered(1, contrast, :));
    
    if remove_outliers
        low_mean = nanmean(low_contrast_data);
        low_std = nanstd(low_contrast_data);

        % Find outliers and set them to NaN
        outlier_indices = abs(low_contrast_data - low_mean) > std_threshold * low_std;
        low_contrast_outliers = sum(outlier_indices);
        low_outliers_count = low_outliers_count + low_contrast_outliers;

        % Print number of outliers for this contrast level
        fprintf('Contrast %d, Low correlation neurons: %d outliers removed\n', contrast, low_contrast_outliers);

        low_contrast_data(outlier_indices) = NaN;
    end

    % Update the filtered data
    low_filtered(1, contrast, :) = low_contrast_data;
end

% Print total outliers removed in first dataset (if applicable)
if remove_outliers
    fprintf('\nFirst dataset (Stationary):\n');
    fprintf('Total outliers removed from high correlation neurons: %d\n', high_outliers_count);
    fprintf('Total outliers removed from low correlation neurons: %d\n', low_outliers_count);
    fprintf('Combined total outliers removed: %d\n\n', high_outliers_count + low_outliers_count);
end

% Calculate means for high and low correlation neurons at each contrast, ignoring NaNs
high_means = squeeze(nanmean(high_filtered, 3));
low_means = squeeze(nanmean(low_filtered, 3));

% Calculate SEM for high and low correlation neurons at each contrast, ignoring NaNs
% First get the standard deviation ignoring NaNs
high_std = squeeze(nanstd(high_filtered, 0, 3));
low_std = squeeze(nanstd(low_filtered, 0, 3));

% Count the actual number of non-NaN values for each contrast
high_count = squeeze(sum(~isnan(high_filtered), 3));
low_count = squeeze(sum(~isnan(low_filtered), 3));

% Calculate SEM using the actual counts of non-NaN values
high_sem = high_std ./ sqrt(high_count);
low_sem = low_std ./ sqrt(low_count);

% Define x-axis (contrast levels)
contrasts = 1:size(norm_diff, 2);
contrast_labels = {'25%', '50%', '100%'};

% Create a figure optimized for PDF output (using standard paper size proportions)
figure('Position', [100, 100, 800, 600], 'PaperPositionMode', 'auto', 'PaperSize', [8.5 11], 'PaperUnits', 'inches');

% First subplot
subplot(1, 2, 1)
hold on;

% Plot high correlation neurons with error bars
errorbar(contrasts, high_means, high_sem, '-o', 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410], 'MarkerFaceColor', [0, 0.4470, 0.7410]);

% Plot low correlation neurons with error bars
errorbar(contrasts, low_means, low_sem, '-o', 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980], 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
xlim([.75 3.25])

% Add labels and legend
xlabel('Contrast');
ylabel('Mean Raw Difference');
title('Stationary (SST)');
legend('High Correlation Neurons', 'Low Correlation Neurons', 'Location', 'best');

% Set x-axis ticks and labels
set(gca, 'XTick', contrasts, 'XTickLabel', contrast_labels, 'TickDir', 'out');

% Customize plot appearance
grid off;
box off;
set(gca, 'FontSize', 12);

% Create filtered copies of the data for the second subplot
high_filtered2 = norm_diff(2, :, highRInds);
low_filtered2 = norm_diff(2, :, lowRInds);

% Initialize counters for outliers in second dataset
high_outliers_count2 = 0;
low_outliers_count2 = 0;

% For each contrast level, optionally remove outliers
for contrast = 1:size(norm_diff, 2)
    % For high correlation neurons
    high_contrast_data = squeeze(high_filtered2(1, contrast, :));
    
    if remove_outliers
        high_mean = nanmean(high_contrast_data);
        high_std = nanstd(high_contrast_data);

        % Find outliers and set them to NaN
        outlier_indices = abs(high_contrast_data - high_mean) > std_threshold * high_std;
        high_contrast_outliers = sum(outlier_indices);
        high_outliers_count2 = high_outliers_count2 + high_contrast_outliers;

        % Print number of outliers for this contrast level
        fprintf('Contrast %d, High correlation neurons (Running): %d outliers removed\n', contrast, high_contrast_outliers);

        high_contrast_data(outlier_indices) = NaN;
    end

    % Update the filtered data
    high_filtered2(1, contrast, :) = high_contrast_data;

    % For low correlation neurons
    low_contrast_data = squeeze(low_filtered2(1, contrast, :));
    
    if remove_outliers
        low_mean = nanmean(low_contrast_data);
        low_std = nanstd(low_contrast_data);

        % Find outliers and set them to NaN
        outlier_indices = abs(low_contrast_data - low_mean) > std_threshold * low_std;
        low_contrast_outliers = sum(outlier_indices);
        low_outliers_count2 = low_outliers_count2 + low_contrast_outliers;

        % Print number of outliers for this contrast level
        fprintf('Contrast %d, Low correlation neurons (Running): %d outliers removed\n', contrast, low_contrast_outliers);

        low_contrast_data(outlier_indices) = NaN;
    end

    % Update the filtered data
    low_filtered2(1, contrast, :) = low_contrast_data;
end

% Print total outliers removed in second dataset (if applicable)
if remove_outliers
    fprintf('\nSecond dataset (Running):\n');
    fprintf('Total outliers removed from high correlation neurons: %d\n', high_outliers_count2);
    fprintf('Total outliers removed from low correlation neurons: %d\n', low_outliers_count2);
    fprintf('Combined total outliers removed: %d\n\n', high_outliers_count2 + low_outliers_count2);

    % Print grand total
    fprintf('GRAND TOTAL outliers removed across both datasets: %d\n', high_outliers_count + low_outliers_count + high_outliers_count2 + low_outliers_count2);
end

% Calculate means for high and low correlation neurons at each contrast, ignoring NaNs
high_means2 = squeeze(nanmean(high_filtered2, 3));
low_means2 = squeeze(nanmean(low_filtered2, 3));

% Calculate SEM for high and low correlation neurons at each contrast, ignoring NaNs
% First get the standard deviation ignoring NaNs
high_std2 = squeeze(nanstd(high_filtered2, 0, 3));
low_std2 = squeeze(nanstd(low_filtered2, 0, 3));

% Count the actual number of non-NaN values for each contrast
high_count2 = squeeze(sum(~isnan(high_filtered2), 3));
low_count2 = squeeze(sum(~isnan(low_filtered2), 3));

% Calculate SEM using the actual counts of non-NaN values
high_sem2 = high_std2 ./ sqrt(high_count2);
low_sem2 = low_std2 ./ sqrt(low_count2);

% Second subplot
subplot(1, 2, 2)
hold on;

% Plot high correlation neurons with error bars
errorbar(contrasts, high_means2, high_sem2, '-o', 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410], 'MarkerFaceColor', [0, 0.4470, 0.7410]);

% Plot low correlation neurons with error bars
errorbar(contrasts, low_means2, low_sem2, '-o', 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980], 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
xlim([.75 3.25])

% Add labels and legend
xlabel('Contrast');
ylabel('Mean Raw Difference');
title('Running (SST)');

% Set x-axis ticks and labels
set(gca, 'XTick', contrasts, 'XTickLabel', contrast_labels, 'TickDir', 'out');

% Customize plot appearance
grid off;
box off;
set(gca, 'FontSize', 12);

% Add more space between subplots
set(gcf, 'Position', get(gcf, 'Position') .* [1 1 1.2 1]);

% Add title to indicate whether outliers were removed
if remove_outliers
    suptitle(sprintf('Analysis with outliers (>%g SD) removed', std_threshold));
else
    suptitle('Analysis without outlier removal');
end

% Save figure with appropriate filename
if remove_outliers
    print('-dpdf', '-bestfit', sprintf('normDiff_vs_contrast_outliers_removed_%gSD.pdf', std_threshold));
else
    print('-dpdf', '-bestfit', 'normDiff_vs_contrast_no_outlier_removal.pdf');
end


%% contrast modulation vs noise corr
% modulation index will be slope of the line between norm diff at 25% and
% at 100%
normDiff_slopes = calculateContrastSlopes(raw_diff);

% Create a new figure
figure;

% Extract x and y data for regression
x_data = noiseCorr_OG_concat{pre}(1,red_ind_concat);
y_data = normDiff_slopes(1,red_ind_concat);

% Remove NaN values for regression analysis
valid_idx = ~isnan(x_data) & ~isnan(y_data) & ~isinf(x_data) & ~isinf(y_data);
x_valid = x_data(valid_idx);
y_valid = y_data(valid_idx);

% Check if we have enough valid data points
if length(x_valid) < 2
    error('Not enough valid data points for regression analysis');
end

% Create scatter plot (using all data points, NaNs will be automatically excluded)
scatter(x_data, y_data, 50, 'filled', 'MarkerFaceAlpha', 0.7);

% Hold the plot to add the regression line
hold on;

% Compute linear regression on valid data
[p, S] = polyfit(x_valid, y_valid, 1);

% Create a sequence of x values for smoother line plotting
x_range = linspace(min(x_valid), max(x_valid), 100);
y_fit = polyval(p, x_range);

% Plot regression line
plot(x_range, y_fit, 'r-', 'LineWidth', 2);

% Calculate R-squared
yresid = y_valid - polyval(p, x_valid);
SSresid = sum(yresid.^2);
SStotal = (length(y_valid)-1) * var(y_valid);
Rsq = 1 - SSresid/SStotal;

% Calculate p-value for the correlation
[R, P] = corrcoef(x_valid, y_valid);
correlation_coef = R(1,2);
p_value = P(1,2);

% Calculate standard error of the slope
n = length(x_valid);
se_slope = sqrt(SSresid/(n-2)) / sqrt(sum((x_valid - mean(x_valid)).^2));

% Calculate t-statistic and p-value for the slope
t_stat = p(1) / se_slope;
slope_p_value = 2 * (1 - tcdf(abs(t_stat), n-2));

% Display only essential regression statistics on the plot
equation = sprintf(' = %.4f', p(1));
r_squared = sprintf('R^2 = %.4f', Rsq);
p_val = sprintf('p = %.4g', p_value);
n_points = sprintf('n = %d', n);

% Create simplified text box with statistics
dim = [.2 .02 .3 .3];
str = {equation, r_squared, p_val, n_points};
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'BackgroundColor', 'white');

% Print only the requested statistics to command window
fprintf('\nRegression Statistics:\n');
fprintf('Beta: %.4f\n', p(1));
fprintf('R-squared: %.4f\n', Rsq);
fprintf('p-value: %.6f\n', p_value);
fprintf('Sample size: %d\n', n);

% Add prediction intervals (95%)
[y_fit_pred, delta] = polyval(p, x_range, S);
plot(x_range, y_fit_pred + delta, 'r--', 'LineWidth', 1);
plot(x_range, y_fit_pred - delta, 'r--', 'LineWidth', 1);

% Customize appearance
set(gca, 'TickDir', 'out');
box off;
xlabel('R value');
ylabel('Raw Diff Slope');
title('R value vs Raw Difference Slope');


hold off;


