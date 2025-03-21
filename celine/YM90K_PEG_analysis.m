
clear all; clear global; close all
clc

% ENTER DATASET NAME
load('ds_YM90K_PEG.mat');
ds_name='YM90K_PEG';
nSess=length(ds_YM90K_PEG);

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
    base = fullfile(isilonName, '\home\ACh\Analysis\2p_analysis\SST_YM90K');
      
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
tc_trial_avrg_stat_largePupil_concat=cell(1,nd);
tc_trial_avrg_stat_smallPupil_concat=cell(1,nd);
tc_trial_avrg_loc_concat=cell(1,nd);
resp_keep_concat=cell(1,nd);
resp_max_keep_concat=cell(1,nd);
pref_responses_loc_concat=cell(1,nd);
pref_responses_stat_concat=cell(1,nd);
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
dataTableConat=[];
drug=cell(1,nSess);
pupilMeans_concat=nan(nd,3,nSess);
pupilCounts_concat=nan(nd,2,nSess);
nonPref_trial_avrg_stat_concat=cell(1,nd);
nonPref_trial_avrg_loc_concat=cell(1,nd);

responseByCondProps_concat=nan(6,2,nSess);

cellID_adjustment=0;
for iSess = 1:nSess
    mouse = ds_YM90K_PEG(iSess).mouse;
    mice=[mice;mouse];
    thisDrug = ds_YM90K_PEG(iSess).drug;
    drug{iSess}=thisDrug;

    if iSess > 1
        cellID_adjustment=max(temp_table.cell_ID_unique); %this should get saved until the next loop;
    end

    if ds_YM90K_PEG(iSess).multiday_timesincedrug_hours>0
        dart_str = [ds_YM90K_PEG(iSess).drug '_' num2str(ds_YM90K_PEG(iSess).multiday_timesincedrug_hours) 'Hr'];
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
%    pupilCounts_concat(:,:,iSess)=pupilCounts;


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
        end
        clear meanF i
    end
    dfof_max_diff_concat=cat(1,dfof_max_diff_concat,dfof_max_diff(:,sharedCon));
   green_fluor_concat=cat(2,green_fluor_concat,green_fluor_keep);
   red_fluor_concat=cat(2,red_fluor_concat,red_fluor_keep);
    
iSess
end

clear mouse iSess nKeep iSess fn_multi cons oris pupilMeans norm_dir_resp_loc norm_dir_resp_stat
clear explanation1 resp_keep sig_diff pref_con_keep pref_ori_keep tOri_match tCon_match data_trial_keep nTrials tc_trial_avrg_keep green_keep_logical red_keep_logical green_ind_keep red_ind_keep
clear LMI RIx locCounts locResp locTCs statResp statTCs wheel_tc ttest_results_stat ttest_results_loc ttest_results_allCon_stat ttest_results_allCon_loc
clear data_con_resp_keep data_ori_resp_keep data_rep_keep dfof_max_diff dfof_max_diff_raw explanation2 resp_max_keep data_resp_keep pref_responses_stat pref_responses_loc
clear tc_trial_avrg_stat tc_trial_avrg_loc fullTC_keep norm_dir_resp sigCorr noiseCorr responseByCondProps
clear red_fluor_all red_fluor_match green_fluor_match green_fluor_match red_fluor_keep green_fluor_keep R_p_values pre_nonPref_stat pre_nonPref_loc
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

%find how many haveRunning red cells exist for each mouse
cellCountsRed = nan(nSess,nCon);
mouseNames=[];
for iMouse = 1:nSess
    for iCon = 1:nCon
        cellCountsRed(iMouse, iCon,1)=length(intersect(haveRunning_red{iCon},(mouseInds{iMouse})));
        
    end
    mouseNames=[mouseNames, string(mice(iMouse,:))];
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
writetable(cellCountTableRed,fullfile(fnout,'cellCounts.csv'),'WriteRowNames',true)
writetable(cellCountTableGreen,fullfile(fnout,'cellCounts_Green.csv'),'WriteRowNames',true)


%find how many haveRunning red cells exist for each mouse
cellCountsRedAll = nan(nSess,1);
for iMouse = 1:nSess
    cellCountsRedAll(iMouse)=length(intersect(red_all,(mouseInds{iMouse})));
end
clear  iMouse


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

dfof_stat_pre = reshape(pref_responses_stat_concat{pre},[3*nKeep_total,1]);
dfof_stat_post = reshape(pref_responses_stat_concat{post},[3*nKeep_total,1]);
dfof_loc_pre = reshape(pref_responses_loc_concat{pre},[3*nKeep_total,1]);
dfof_loc_post = reshape(pref_responses_loc_concat{post},[3*nKeep_total,1]);

dfof_col = vertcat(dfof_stat_pre,dfof_stat_post,dfof_loc_pre,dfof_loc_post);

cellID=1:nKeep_total;
cellID_col=repmat(cellID',12,1);

mouseIDcol = repmat(mouseID,12,1);

day1 = repelem(["pre" "post"],1,(3*nKeep_total))';
day=repmat(day1,2,1);


contrast1 = repelem(cons,nKeep_total)';
contrast=repmat(contrast1,4,1);

cell_type_col=repmat(red_concat,1,12)';

behStateCol = repelem(["stat" "loc"],1,(6*nKeep_total))';

clear dfof_stat_pre dfof_stat_post dfof_loc_pre dfof_loc_post cellID1 day1 contrast1



dfof_summary = table(mouseIDcol,cellID_col,cell_type_col,contrast,behStateCol,day,dfof_col, ...
    'VariableNames',{'mouseID' 'cellID' 'cellType' 'contrast' 'behState' 'day' 'dfof'});



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

% cells with high correlation in the baseline day
highRInds = find(noiseCorr_OG_concat{pre}(1,:)>0.5);
lowRInds = find(noiseCorr_OG_concat{pre}(1,:)<=0.5);

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
%% Figure 2C

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
ylim([-.02 .17]);
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
print(fullfile(fnout,'Fig_S3A_vertical.pdf'),'-dpdf');


%% Fig 2D - contrast response for SST and Pyr cells

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
%ylabel('dF/F') 
xlim([0 1.1])
ylim([.0 .12])
xticks([.25 .5 1])
%set(gca, 'TickDir', 'out')
box off

subplot(1,2,2) %
errorbar(cons,conResp_green_avrg_stat{pre},conResp_green_se_stat{pre},'k');
hold on
errorbar(cons,conResp_green_avrg_stat{post},conResp_green_se_stat{post},'b');
%title(['Pyr',' n = ', num2str(length(green_ind_concat))])
%xlabel('contrast') 
%set(gca, 'TickDir', 'out')
xlim([0 1.1])
ylim([.0 .12])
xticks([.25 .5 1])
box off


x0=5;
y0=5;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(fullfile(fnout,'Fig_S3B.pdf'),'-dpdf');


%% Figure 2D statistics
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

%% Figure 2E

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
ylim([0 .4])
ylabel(["Fraction SST cells"]) 
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
print(fullfile(fnout,'Fig_2E.pdf'),'-dpdf')

%% Figure 2E statistics

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

%% Supplemental figure related to Fig 2 - normalized difference

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
%% Fig 3B histogram
figure;
histogram(noiseCorr_concat{pre}(1,red_ind_concat),'FaceColor','black','BinWidth',.1)
xlabel('Correlation to Pyr activity')
ylabel('# Cells')
set(gca,'TickDir','out')
box off
x0=3;
y0=3;
width=2.5;
height=2;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'Fig_3B.pdf'),'-dpdf');
%% Fig 3C

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

%Makes a plot for each contrast - 50 contrast used in paper
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
     if iCon==2
        print(fullfile(fnout,'Fig_3C.pdf'),'-dpdf');
     % else %option to print the plots for 25 and 100 contrast as well,
     % rather than just displaying them
     %     print(fullfile(fnout,[num2str(cons(iCon)) 'stat_R_timecourses.pdf']),'-dpdf');
     end
    clear txt1 highRed lowRed
end 

%% Fig 3D - contrast response split by correlation value

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
% xlabel('contrast') 
% ylabel('dF/F, pref ori') 
xlim([0 1.1])
ylim([0 .15])
xticks([.25 .5 1])
% set(gca, 'TickDir', 'out')
box off

subplot(1,2,2) %
errorbar(cons,conResp_redHigh_avrg_stat{pre},conResp_redHigh_se_stat{pre},'k');
hold on
errorbar(cons,conResp_redHigh_avrg_stat{post},conResp_redHigh_se_stat{post},'b');
title(['Strong Corr,' 'n = ', num2str(length(redHigh))])

xlabel('contrast') 
set(gca, 'TickDir', 'out')
xlim([0 1.1])
ylim([0 .15])
xticks([.25 .5 1])
box off

x0=5;
y0=5;
width=2.5;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])


print(fullfile(fnout,'Fig_S4D.pdf'),'-dpdf');

%% Figure 3D statistics
%from the full dataframe, extract rows for stationary trials for weakly and
%strongly correlated SST cells
SST_low_stat_dfof = SST_stat_dfof(ismember(SST_stat_dfof.cellID,redLow),:);
SST_high_stat_dfof = SST_stat_dfof(ismember(SST_stat_dfof.cellID,redHigh),:);

fullCorrTable = vertcat(SST_low_stat_dfof,SST_high_stat_dfof);
fullCorrTable.corrType = categorical([repmat(0,length(redLow),1);repmat(1,length(redHigh),1)]);

% run an ANOVA for corr type X contrast, baseline day only
w = table(categorical([1 2 3].'), 'VariableNames', {'contrast'}); % within-design
rm_fullCorr = fitrm(fullCorrTable, 'd1c1-d1c3 ~ corrType', 'WithinDesign', w)
ranova(rm_fullCorr, 'withinmodel', 'contrast');

% run the ANOVAs on low and high corr seperately
w = table(categorical([1 1 1 2 2 2 ].'), categorical([1 2 3 1 2 3].'), 'VariableNames', {'DART', 'contrast'}); % within-design
rm_SST_low = fitrm(SST_low_stat_dfof, 'd1c1-d2c3 ~ 1', 'WithinDesign', w);
ranova(rm_SST_low, 'withinmodel', 'DART*contrast')

rm_SST_high = fitrm(SST_high_stat_dfof, 'd1c1-d2c3 ~ 1', 'WithinDesign', w);
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



%% Fig 4 - timecourses for running trials

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
print(fullfile(fnout,'Fig_4A_vertical.pdf'),'-dpdf');

%% Fig S4 A - related to Fig 4A - timecourses for Pyr cells running trials
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
print(fullfile(fnout,'Fig_S4A_vertical.pdf'),'-dpdf');

%% Figure 4B
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
subplot(2,2,1)
errorbar(cons,conResp_red_avrg_stat{pre},conResp_red_se_stat{pre},'k');
hold on
errorbar(cons,conResp_red_avrg_stat{post},conResp_red_se_stat{post},'b');
title(['Stationary'])
ylabel('dF/F') 
xlabel('contrast') 
xlim([0 1])
ylim([0 .3])
xticks([.25 .5 1])
set(gca, 'TickDir', 'out')
box off


subplot(2,2,2) 
errorbar(cons,conResp_red_avrg_loc{pre},conResp_red_se_loc{pre},'k');
hold on
errorbar(cons,conResp_red_avrg_loc{post},conResp_red_se_loc{post},'b');
title(['Running'])
ylim([0 .3])
xlim([0 1])
xticks([.25 .5 1])
xlabel('contrast') 
set(gca, 'TickDir', 'out')
box off

subplot(2,2,3)
errorbar(cons,conResp_green_avrg_stat{pre},conResp_green_se_stat{pre},'--k');
hold on
errorbar(cons,conResp_green_avrg_stat{post},conResp_green_se_stat{post},'--b');
%title(['Stationary'])
ylabel('dF/F') 
xlabel('contrast') 
xlim([0 1])
ylim([0 .3])
xticks([.25 .5 1])
set(gca, 'TickDir', 'out')
box off


subplot(2,2,4) 
errorbar(cons,conResp_green_avrg_loc{pre},conResp_green_se_loc{pre},'--k');
hold on
errorbar(cons,conResp_green_avrg_loc{post},conResp_green_se_loc{post},'--b');
%title(['Running'])
ylim([0 .3])
xlim([0 1])
xticks([.25 .5 1])
xlabel('contrast') 
set(gca, 'TickDir', 'out')
box off



x0=5;
y0=5;
width=2.5;
height=2.5;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(fullfile(fnout,'Fig_4B.pdf'),'-dpdf');
%% Fig 4B statistics
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
stat_dfof_summary_full = [labels_table,dfof_table];
matched_dfof_summary = stat_dfof_summary_full(ismember(stat_dfof_summary_full.cellID,all_cells),:);


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

rm_SST_stat = fitrm(SST_matched_dfof, 'Sd1c1-Sd2c3 ~ 1', 'WithinDesign', w);
ranova(rm_SST_stat, 'withinmodel', 'DART*contrast')

rm_SST_loc = fitrm(SST_matched_dfof, 'Rd1c1-Rd2c3 ~ 1', 'WithinDesign', w);
ranova(rm_SST_loc, 'withinmodel', 'DART*contrast')

rm_Pyr_stat = fitrm(Pyr_matched_dfof, 'Sd1c1-Sd2c3 ~ 1', 'WithinDesign', w);
ranova(rm_Pyr_stat, 'withinmodel', 'DART*contrast')

rm_Pyr_loc = fitrm(Pyr_matched_dfof, 'Rd1c1-Rd2c3 ~ 1', 'WithinDesign', w);
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



% pairwise ttests for dfof response at each contrast for pyr cells
[pyr_h1, pyr_p1]= ttest(pref_responses_stat_concat{pre}(green_all,1),pref_responses_stat_concat{post}(green_all,1));
[pyr_h2, pyr_p2]= ttest(pref_responses_stat_concat{pre}(green_all,2),pref_responses_stat_concat{post}(green_all,2));
[pyr_h3, pyr_p3]= ttest(pref_responses_stat_concat{pre}(green_all,3),pref_responses_stat_concat{post}(green_all,3));

%corrected for three tests
pyr_pvalues_stat = [(pyr_p1*3);(pyr_p2*3);(pyr_p3*3)];
contrasts = cons';
table(contrasts,pyr_pvalues_stat)

% pairwise ttests for dfof response at each contrast for pyr cells
[pyr_h1, pyr_p1]= ttest(pref_responses_loc_concat{pre}(green_all,1),pref_responses_loc_concat{post}(green_all,1));
[pyr_h2, pyr_p2]= ttest(pref_responses_loc_concat{pre}(green_all,2),pref_responses_loc_concat{post}(green_all,2));
[pyr_h3, pyr_p3]= ttest(pref_responses_loc_concat{pre}(green_all,3),pref_responses_loc_concat{post}(green_all,3));

%corrected for three tests
pyr_pvalues_loc = [(pyr_p1*3);(pyr_p2*3);(pyr_p3*3)];
contrasts = cons';
table(contrasts,pyr_pvalues_loc)
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
xticklabels({'25','50','100'})
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

%% Fig 4D  response by condition

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
%title('Mean response to different conditions')


% 
% 
% for iMouse = 1:nSess 
%     x = responseByCondProps_concat(5,1,iMouse):.005:responseByCondProps_concat(6,1,iMouse);
%     m  = responseByCondProps_concat(1,1,iMouse);
%     y = (m*x)+responseByCondProps_concat(2,1,iMouse);
%     plot(x,y,'Color',"#808080") 
% 
% 
%     x = responseByCondProps_concat(5,2,iMouse):.005:responseByCondProps_concat(6,2,iMouse);
%     m  = responseByCondProps_concat(1,2,iMouse);
%     y = (m*x)+responseByCondProps_concat(2,2,iMouse);
%     plot(x,y,'Color',"#0080ff")  
% 
%     x_population=min(responseByCondProps_concat(5,1,:)):.005:max(responseByCondProps_concat(6,1,:));
%     m_mean=mean(responseByCondProps_concat(1,1,:),'omitmissing');
%     b_mean=mean(responseByCondProps_concat(2,1,:),'omitmissing');
%     y_mean=(m_mean*x_population)+b_mean;
%     plot(x_population,y_mean,'k','LineWidth',2)
% 
%     x_population=min(responseByCondProps_concat(5,2,:)):.005:max(responseByCondProps_concat(6,2,:));
%     m_mean=mean(responseByCondProps_concat(1,2,:),'omitmissing');
%     b_mean=mean(responseByCondProps_concat(2,2,:),'omitmissing');
%     y_mean=(m_mean*x_population)+b_mean;
%     plot(x_population,y_mean,'b','LineWidth',2)
% 
% end

ylim([-.01 .3])
xlim([-.01 .3])
x0=5;
y0=5;
width=2;
height=2;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'Fig_4D.pdf'),'-dpdf')

%% Fig 4E % F
% property X pre/post (1=pre, 2=post) X mouse
%properties: range = 5-6, slope = 1, intercept = 2

slopes=squeeze(responseByCondProps_concat(1,:,:));
slopes_mean = mean(slopes,2,'omitmissing');
slopes_std = std(slopes,[],2,"omitmissing");
slopes_se = slopes_std/sum(cellCountsRedAll>0);

intercepts=squeeze(responseByCondProps_concat(2,:,:));
intercepts_mean = mean(intercepts,2,'omitmissing');
intercepts_std = std(intercepts,[],2,"omitmissing");
intercepts_se = intercepts_std/sum(cellCountsRedAll>0);

figure
subplot(1,2,1)
plot(slopes, "-o",'Color',"#7b7b7b")
% xlim([.5 2.5])
box off
set(gca,'TickDir','out')
xticks([1 2])
xticklabels({'Pre-DART','Post-DART'})
ylabel('Slope')
hold on
errorbar(slopes_mean, slopes_se,"-*K",'LineWidth',2)
%plot(nanmean(slopes,2),"-*K",'LineWidth',2)
%title("Slope change")

subplot(1,2,2)
plot(intercepts, "-o",'Color',"#7b7b7b")
xlim([.5 2.5])
% ylim([-.3 .1])
box off
set(gca,'TickDir','out')
xticks([1 2])
xticklabels({'Pre-DART','Post-DART'})
ylabel('Intercept')
hold on
errorbar(intercepts_mean, intercepts_se,"-*K",'LineWidth',2)
%plot(nanmean(intercepts,2),"-*K",'LineWidth',2)
%title("Intercept change")



x0=5;
y0=5;
width=3;
height=2;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'Fig_4EF.pdf'),'-dpdf')

%% Fig 4E-F stats
[h1, p1]= ttest(slopes(1,:),slopes(2,:));
[h2, p2]= ttest(intercepts(1,:),intercepts(2,:))

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
     if iCon==2
        print(fullfile(fnout,'Fig_Sx_C.pdf'),'-dpdf');
     % else %option to print the plots for 25 and 100 contrast as well,
     % rather than just displaying them
     %     print(fullfile(fnout,[num2str(cons(iCon)) 'stat_R_timecourses.pdf']),'-dpdf');
     end
end 






%% Fig Sy A, SST cell timecourses split by pupil size, and Sy B, Pyr responses split by pupil size

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
    temp_green = intersect(green_ind_concat,have_bothPupil{iCon});
    temp_red = intersect(red_ind_concat,have_bothPupil{iCon});

    tc_green_avrg_stat_large{id}(:,iCon)=nanmean(tc_trial_avrg_stat_largePupil_concat{id}(:,temp_green,iCon),2);
    green_std=nanstd(tc_trial_avrg_stat_largePupil_concat{id}(:,temp_green,iCon),[],2);
    tc_green_se_stat_large{id}(:,iCon)=green_std/sqrt(length(temp_green));
    
    tc_red_avrg_stat_large{id}(:,iCon)=nanmean(tc_trial_avrg_stat_largePupil_concat{id}(:,temp_red,iCon),2);
    red_std=nanstd(tc_trial_avrg_stat_largePupil_concat{id}(:,temp_red,iCon),[],2);
    tc_red_se_stat_large{id}(:,iCon)=red_std/sqrt(length(temp_red));
    
    clear green_std red_std

    tc_green_avrg_stat_small{id}(:,iCon)=nanmean(tc_trial_avrg_stat_smallPupil_concat{id}(:,temp_green,iCon),2);
    green_std=nanstd(tc_trial_avrg_stat_smallPupil_concat{id}(:,temp_green,iCon),[],2);
    tc_green_se_stat_small{id}(:,iCon)=green_std/sqrt(length(temp_green));
    
    tc_red_avrg_stat_small{id}(:,iCon)=nanmean(tc_trial_avrg_stat_smallPupil_concat{id}(:,temp_red,iCon),2);
    red_std=nanstd(tc_trial_avrg_stat_smallPupil_concat{id}(:,temp_red,iCon),[],2);
    tc_red_se_stat_small{id}(:,iCon)=red_std/sqrt(length(temp_red));
    
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
sgtitle(['SST',' n = ', num2str(length(red_all))])

x0=5;
y0=0;
width=4;
height=9;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')


end  
print(fullfile(fnout,'Fig_Sy_A.pdf'),'-dpdf');

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
sgtitle(['Pyr',' n = ', num2str(length(green_all))])

x0=5;
y0=0;
width=4;
height=9;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')


end  


print(fullfile(fnout,'Fig_Sy_B.pdf'),'-dpdf');

%% get a table of capture values
capture = getCaptureValues_annulus_peg(mice(2:length(mice)));
table(mice(2:length(mice)),capture(3,:)')
edges = linspace(1, 2, 10); % Create 20 bins.
histogram(capture(3,:),'BinEdges',edges);
xlim([1 2])
box off
set(gca, 'TickDir', 'out')
x0=5;
y0=5;
width=1.1;
height=1.1;
set(gcf,'units','inches','position',[x0,y0,width,height])

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
    text(text_x, text_y_start - text_y_step, sprintf('R² = %.3f', r_squared(contrast)), 'FontSize', 10);
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