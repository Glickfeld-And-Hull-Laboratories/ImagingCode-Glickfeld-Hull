clear all; clear global; 
close all
clc
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc =  behavConstsDART; %directories
eval(ds);

%day_id = 169; %enter post-DART day
day_id = input('Enter day id ');% alternative to run from command line.
pre_day = expt(day_id).multiday_matchdays;

nd=2; %hardcoding the number of days for now

mouse = expt(day_id).mouse;

fnout = fullfile(rc.achAnalysis,mouse);
if expt(day_id).multiday_timesincedrug_hours>0
    dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end


fn_multi = fullfile(rc.achAnalysis,mouse,['multiday_' dart_str]);

cd(fn_multi)
load(fullfile(fn_multi,'timecourses.mat'))
%load(fullfile(fn_multi,'multiday_alignment.mat'))
load(fullfile(fn_multi,'input.mat'))
frame_rate = input.frameImagingRateMs;
%% finding red fluorescence level
 allDays = [day_id,pre_day];

for id = 1 %currently only doing this for the baseline day
mouse = expt(allDays(id)).mouse;
date = expt(allDays(id)).date;
imgFolder = expt(allDays(id)).contrastxori_runs{1};
fn = fullfile(rc.achAnalysis,mouse,date,imgFolder);
cd(fn);
load(fullfile(fn,'redImage.mat'));
load(fullfile(fn,'mask_cell.mat'));

%for each day individually, load the masks and red image for that day
    
    
% cell_stats=regionprops(mask_cell_red);
% figure; imagesc(redChImg), colormap gray; caxis([200 1000]);
% hold on
% bound = cell2mat(bwboundaries(mask_cell_red(:,:,1)));
% plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',.5); hold on;
% for iC = 1:max(max(mask_cell_red))
%     text(cell_stats(iC).Centroid(1), cell_stats(iC).Centroid(2), num2str(iC), 'Color', 'red',...
%             'Fontsize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
%     
% end
 
%use stackGetTimeCourses to extract the red fluorescence within each mask
red_fluor_mask = stackGetTimeCourses(redChImg, mask_cell);
nCells=max(max(mask_cell));
for i = 1:nCells
red_fluor_np(i) = stackGetTimeCourses(redChImg, mask_np(:,:,i));
end

red_fluor_all = red_fluor_mask-red_fluor_np;

% clear mask_cell mask_np nCells red_fluor_np red_fluor_mask
end
%using the reference day
red_fluor_match=red_fluor_all(:,match_ind);
z_red_fluor=zscore(red_fluor_match);
load(fullfile(fn_multi,'multiday_alignment.mat'))
clear red_fluor_all red_fluor_mask red_fluor_np
%% get green fluor level
%using the reference day
green_fluor_match=mean(cellTCs_match{1},1);   




%% stimulus props

nOn = input(1).nScansOn;
nOff = input(1).nScansOff;

%tells the contrast, direction and orientation for each trial each day
tCon_match = cell(1,nd);
tDir_match = cell(1,nd);
tOri_match = cell(1,nd);

%find the contrasts, directions and orientations for each day
%in case of instances where the number of trails actually collected was not
%consistent with the number mWorks thinks occured, only take 1:nTrials as
%dictated above based on the number of frames recorded
for id = 1:nd
    nTrials(id) = size(cellTCs_match{id},1)/(nOn+nOff); %to account for times 
%when there is a disruption before the full set of trials is collcted, I'm 
%determining the number of trials each day by how many frames of data I 
%have divided by the number of frames per trial
    tCon_match{id} = celleqel2mat_padded(input(id).tGratingContrast(1:nTrials(id)));
    tDir_match{id} = celleqel2mat_padded(input(id).tGratingDirectionDeg(1:nTrials(id)));
    tOri_match{id} = tDir_match{id};
    tOri_match{id}(find(tDir_match{id}>=180)) = tDir_match{id}(find(tDir_match{id}>=180))-180;
end
oris = unique(tOri_match{1});
dirs = unique(tDir_match{1});
cons = unique(tCon_match{1});
nOri = length(oris);
nCon = length(cons);
nDir = length(dirs);

%% convert to trials

%change this to use a padded array, where I add zeros at the end. test=padarray(cellTCs_match{1},30,0,'post');
stimStart = (nOff/2)+1; %this indicates both the perdiod to trim off the start and the stim on period after trimming
stimEnd=stimStart+nOn-1;

cellTCs_match_OG = cellTCs_match;
data_dfof_trial_match = cell(1,nd); %make an empty array that is 1 by however many days there are (1X2 usually)

fractTimeActive_match = cell(1,nd);
cellstd_match = cell(1,nd);
for id = 1:nd %cycle through days
    
  %nTrials(id) = length(tDir_match{id}); %use the list of direction by trial to figure out how many trials there are
   %currently the way I center the stim on period requires me to cut out
   %one trial, hence the -1
   
    %for matched cell
    nCells = size(cellTCs_match{id},2);
    fractTimeActive_match{id} = zeros(1,nCells);
    %I will trim 30 frames of the start and add 30 frames of padding to the
    %end (padding is nans)
    cellTCs_match{id} = cellTCs_match{id}(stimStart:size(cellTCs_match{id},1),:);
    %cellTCs_match{id} = padarray(cellTCs_match{id},30,0,'post');
    cellTCs_match{id} = padarray(cellTCs_match{id},30,999,'post'); %I had been 
    %padding with zeroas but that altered the TC so I changed this to pad with NANs 
    cellTCs_match{id}(cellTCs_match{id}==999)=nan;
    nFrames = size(cellTCs_match{id},1);
    
    
    data_trial_match = reshape(cellTCs_match{id},[nOn+nOff nTrials(id) nCells]);
    data_f_match = mean(data_trial_match(1: (nOff/2),:,:),1);
    data_dfof_trial_match{id} = bsxfun(@rdivide,bsxfun(@minus,data_trial_match,data_f_match),data_f_match);
    meansub_match = cellTCs_match{id}-nanmean(cellTCs_match{id},1);
    cellstd = nanstd(meansub_match,[],1);
    cellstd_match{id}=cellstd;
    for iCell = 1:nCells
        fractTimeActive_match{id}(:,iCell) = length(find(meansub_match(:,iCell)>3.*cellstd(1,iCell)))./nFrames;
    end
    
end
clear  data_f_match cellstd 

% MUST USE NANMEAN INSTEAD OF MEAN MOVING FORWARD SINCE I SET THE PADDING
% VALUES TO NAN

%% find significant responses and preferred stimuli

resp_win = stimStart:(stimEnd+5);
base_win = 1: stimStart-1;

%make cell arrays to keep everything in
data_resp_match = cell(1,nd);
h_match = cell(1,nd);
p_match = cell(1,nd);
h_all_match = cell(1,nd);
resp_match = cell(1,nd);
pref_dir_match = cell(1,nd);
pref_con_match = cell(1,nd);
data_dir_resp_match = cell(1,nd);
data_con_resp_match = cell(1,nd);

%rect_dfof_trial=data_dfof_trial_match;
%rect_dfof_trial(find(rect_dfof_trial<0))=0;
% for each day
for id = 1:nd
    data_resp = zeros(nCells, nDir, nCon,2);
    h = zeros(nCells, nDir, nCon);
    p = zeros(nCells, nDir, nCon);
    passThresh=zeros(nCells, nDir, nCon);
    tCon = tCon_match{id}(:,1:nTrials(id));
    tOri = tOri_match{id}(:,1:nTrials(id));
    tDir = tDir_match{id}(:,1:nTrials(id));
    data_dfof_trial = data_dfof_trial_match{id}; 

    for iDir = 1:nDir
        ind_dir = find(tDir == dirs(iDir));
        for iCon = 1:nCon
            ind_con = find(tCon == cons(iCon));
            ind = intersect(ind_dir,ind_con); %for every orientation and then every contrast, find trials with that con/ori combination
            data_resp(:,iDir,iCon,1) = squeeze(nanmean(nanmean(data_dfof_trial(resp_win,ind,:),1),2));
            data_resp(:,iDir,iCon,2) = squeeze(std(nanmean(data_dfof_trial(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
            [h(:,iDir,iCon), p(:,iDir,iCon)] = ttest(nanmean(data_dfof_trial(resp_win,ind,:),1), nanmean(data_dfof_trial(base_win,ind,:),1),'dim',2,'tail','right','alpha',0.01./(nDir*nCon-1));
%             baseStd=squeeze(std(nanmean(data_dfof_trial(base_win,ind,:),1),[],2));
%             baseMean=squeeze(nanmean(nanmean(data_dfof_trial(base_win,ind,:),1),2));
%             thresh=baseMean + (3.*baseStd);
%             passThresh(:,iOri,iCon) = logical(data_resp(:,iOri,iCon,1) > thresh);
%             passThresh(:,iOri,iCon) = logical(data_resp(:,iOri,iCon,1) > .05);
%             h(:,iOri,iCon) = h(:,iOri,iCon) .*passThresh(:,iOri,iCon);
        
        end
    end 
    
    resp_sig = data_resp(:,:,:,1).*h;
    h_pass = sum(sum((h),2),3); 
    resp=logical(h_pass);
    
    pref_dir = zeros(1,nCells);
    pref_con = zeros(1,nCells);
    data_dir_resp = zeros(nCells,nDir); %at pref con
    data_con_resp = zeros(nCells,nCon); %at pref ori
    data_orth_resp=zeros(nCells,1);
    %I want to pull out the responses for each cell at it's preferred orientations, for
    %all contrasts, and at it's preferred contrast, for all orientations
    for iCell = 1:nCells
          [max_val, pref_dir(1,iCell)] = max(max(resp_sig(iCell,:,:),[],3));
          [max_val_con, pref_con(1,iCell)] = max(max(resp_sig(iCell,:,:),[],2)); 
          data_dir_resp(iCell,:)=data_resp(iCell,:,pref_con(iCell),1);
          data_con_resp(iCell,:)=data_resp(iCell,pref_dir(iCell),:,1);
          
%       if pref_ori(iCell)<= nDir/2
%           orth_ori(iCell)=pref_ori(iCell)+2;
%       elseif pref_ori(iCell)>nOri/2
%           orth_ori(iCell)=pref_ori(iCell)-2;
%       end
%     data_orth_resp(iCell,:)=nanmean(data_resp(iCell,orth_ori(iCell),:,1));
    end

 %then put into a cell matrix for the two days
data_resp_match{id} = data_resp;
h_match{id} = h;
p_match{id} = p;
h_all_match{id} = h_pass;
resp_match{id} = resp;
pref_dir_match{id} = pref_dir;
pref_con_match{id} = pref_con;
data_dir_resp_match{id} = data_dir_resp;
data_con_resp_match{id} = data_con_resp;
% data_orth_resp_match{id}=data_orth_resp;
 
end
clear data_resp p h_all resp pref_dir pref_con data_dir_resp data_con_resp data_dfof_trial tCon tOri tDir data_orth_resp  baseStd baseMean thresh pass
 
%% get basic counts 
red_ind_match_list = find(red_ind_match==1);

%this is a list of indices of all the green cells
green_ind_match = ~(red_ind_match);
green_ind_match_list = find(green_ind_match);


%find cells that were matched and were active
green_match_respd1 = intersect(green_ind_match_list,find(resp_match{2}==1));
green_match_respd2 = intersect(green_ind_match_list,find(resp_match{1}==1));

red_match_respd1 = intersect(red_ind_match_list,find(resp_match{2}==1));%this is for reverse matching
red_match_respd2 = intersect(red_ind_match_list,find(resp_match{1}==1));

%%
%find cells that were active on at least one day
resp_green_either_temp = union(green_match_respd1,green_match_respd2); %these are the green cells I will include from now on
resp_red_either_temp = union(red_match_respd1,red_match_respd2); %these are the red cells I will include from now on

keep_cells_temp = union(resp_green_either_temp,resp_red_either_temp);

outliers=cell(1,nd);
for id = 1:nd
  data_resp_keep_temp=data_resp_match{id}(keep_cells_temp,:,:,:,1);
  resp_max_keep_temp = max(squeeze(max(data_resp_keep_temp(:,:,:,1),[],3)),[],2);
  mean_max = nanmean(resp_max_keep_temp);
  std_max=std(resp_max_keep_temp);
  thresh=mean_max+(3*std_max);
  outliers{id}=keep_cells_temp(find(resp_max_keep_temp>thresh));
end
outliers_either=union(outliers{1},outliers{2});
fprintf(['removing ' num2str(length(outliers_either)) ' outlier cells'])

resp_green_either=setdiff(resp_green_either_temp,outliers_either);
resp_red_either=setdiff(resp_red_either_temp,outliers_either);

clear resp_green_either_temp resp_red_either_temp keep_cells_temp data_resp_keep_temp resp_max_keep_temp mean_max std_max thresh outliers
%Impose a constraint where if the max response of a cell is > 2 std from
%the max resposne of other cells I'm keeping, I drop it? 
%%
%these are the cells I will include from now on
%find cells that were active on at least one day
% resp_green_either = union(green_match_respd1,green_match_respd2); %these are the green cells I will include from now on
% resp_red_either = union(red_match_respd1,red_match_respd2); %these are the red cells I will include from now on
keep_cells = union(resp_green_either,resp_red_either);
nGreen_keep = length(resp_green_either); %how many green cells to keep
nRed_keep = length(resp_red_either);%how many red cells to keep
nKeep = length(keep_cells)

%making a list of the day 1 indices for the cells I want to keep
%find the indices
[resp_green_either,green_ind_keep] = intersect(keep_cells,resp_green_either,'stable');
[resp_red_either,red_ind_keep] = intersect(keep_cells,resp_red_either,'stable');

nGreen_keep_respd1 = length(green_match_respd1); %how many of those responded on d1
nGreen_keep_respd2 = length(green_match_respd2);%how many of those responded on d2

nRed_keep_respd1 = length(red_match_respd1); %how many of those responded on d1
nRed_keep_respd2 = length(red_match_respd2);%how many of those responded on d2


% make table of values
countsTable = table([nGreen_keep;nRed_keep],[nGreen_keep_respd1;nRed_keep_respd1],[nGreen_keep_respd2;nRed_keep_respd2],'VariableNames',{'Keep' 'Responsive pre' 'Responsive post'}, 'RowNames',{'Pyramidal cells'  'HT+ cells'})
writetable(countsTable,fullfile(fn_multi,'match_counts.csv'),'WriteRowNames',true)
clear  nGreen_match_respd1 nGreen_match_respd2  nRed_match_respd1 nRed_match_respd2

%% make a data structure subsets for only the keep cells
data_trial_keep=cell(1,nd);
pref_dir_keep=cell(1,nd);
pref_con_keep=cell(1,nd);
resp_keep=cell(1,nd);


for id = 1:nd
    data_trial_keep{id} = data_dfof_trial_match{id}(:,:,keep_cells);
    pref_dir_keep{id} = pref_dir_match{id}(:,keep_cells);
    pref_dir_keep{id}=dirs(pref_dir_keep{id});
    pref_con_keep{id} = pref_con_match{id}(:,keep_cells);
    pref_con_keep{id}=cons(pref_con_keep{id});
    resp_keep{id} = resp_match{id}(keep_cells);

end

red_fluor_keep=red_fluor_match(keep_cells);
green_fluor_keep=green_fluor_match(keep_cells);

conTable = table([mean(pref_con_keep{2}(green_ind_keep));mean(pref_con_keep{2}(red_ind_keep))],[mean(pref_con_keep{1}(green_ind_keep));mean(pref_con_keep{1}(red_ind_keep))],'VariableNames',{'mean pref con pre' 'mean pref con post'}, 'RowNames',{'Pyramidal cells'  'HT+ cells'})
writetable(conTable,fullfile(fn_multi,'conPref.csv'),'WriteRowNames',true)
save(fullfile(fn_multi,'fluor_intensity.mat'),'red_fluor_match','green_fluor_match','green_fluor_match','red_fluor_keep','green_fluor_keep')

%% looking at wheel speed
wheel_speed = cell(1,nd);

for id = 1:nd
    wheel_speed{id} = wheelSpeedCalc(input(id),32,expt(allDays(1)).wheelColor); 
    nanmean(wheel_speed{id})
end



wheel_tc = cell(1,nd);
wheel_trial_avg= cell(1,nd);
RIx = cell(1,nd);

for id = 1:nd
    wheel_tc{id}=zeros(nOn+nOff, nTrials(id));
    for iTrial = 1:nTrials(id)
        wheel_tc{id}(:,iTrial) = wheel_speed{id}(1+((iTrial-1).*(nOn+nOff)):iTrial.*(nOn+nOff));
    end
    wheel_trial_avg{id} = mean(wheel_tc{id}(nOff:nOn+nOff,:),1);
    RIx{id} = wheel_trial_avg{id}>2; %.55 is the noise level in the wheel movement
    mean(RIx{id})
end

%% narrow down to the stimuli preferred for each cell each day
%we will get one tc per cell per contrast. This represents that cell's tc
%averaged over trials at the preferred orientation, at each contrast
%only include stationary trials

tc_trial_avrg_stat=cell(1,nd);
tc_trial_avrg_loc=cell(1,nd);
rect_tc_trial_avrg_keep=cell(1,nd);
tc_trial_avrg_keep_allCond=cell(1,nd);
pref_responses_stat = cell(1,nd);
pref_responses_loc = cell(1,nd);
pref_responses_allCond = cell(1,nd);
trialCounts=cell(2,nd);


for id = 1:nd
    trialCounts{1,id}=[];
    trialCounts{2,id}=[];
    temp_tc_stat=nan((nOn+nOff),nKeep,nCon);
    temp_tc_loc=nan((nOn+nOff),nKeep,nCon);
    temp_pref_responses_stat=zeros(nKeep,nCon,1);
    temp_pref_responses_loc=zeros(nKeep,nCon,1);
    for iCon = 1:nCon
        for i=1:nKeep
            
            temp_TCs=data_trial_keep{id}(:,:,i); %only pulling from dfof data of keep cells
            tCon=tCon_match{id}(1:nTrials(id));
            tDir=tDir_match{id}(1:nTrials(id));
            %identify the trials where ori = pref ori
            temp_dir= pref_dir_keep{id}(i); %find the preferred ori of this cell (already in degrees)
            dir_inds = find(tDir==temp_dir); %these are the trials at that ori
            con_inds=find(tCon==cons(iCon));
            stat_inds = find(~RIx{id});
            loc_inds = find(RIx{id});
            temp_trials1 = intersect(dir_inds, con_inds); %preferred ori for this cell, looping through all cons
            temp_trials_stat = intersect(temp_trials1,stat_inds);
            temp_trials_loc = intersect(temp_trials1,loc_inds);
            temp_pref_responses_stat(i,iCon)=nanmean(nanmean(temp_TCs(stimStart:stimEnd,temp_trials_stat),1),2);
            temp_pref_responses_loc(i,iCon)=nanmean(nanmean(temp_TCs(stimStart:stimEnd,temp_trials_loc),1),2);
            temp_pref_responses_allCond(i,iCon)=nanmean(nanmean(temp_TCs(stimStart:stimEnd,temp_trials1),1),2);
            
            temp_tc_stat(:,i,iCon)=nanmean(temp_TCs(:,temp_trials_stat),2);
            temp_tc_loc(:,i,iCon)=nanmean(temp_TCs(:,temp_trials_loc),2);
            temp_tc_allCond(:,i,iCon)=nanmean(temp_TCs(:,temp_trials1),2);%average over all trials regardless of stationary vs. locomotion
            rect_tc_trial_avrg=temp_tc_stat;
            rect_tc_trial_avrg(rect_tc_trial_avrg<0)=0;
            
            
            trialCounts{1,id}=[trialCounts{1,id},length(temp_trials_stat)];
            trialCounts{2,id}=[trialCounts{2,id},length(temp_trials_loc)];
            length(temp_trials1);
        end

    end
    
    
tc_trial_avrg_stat{id}=temp_tc_stat; %this is a cell array with one cell 
tc_trial_avrg_loc{id}=temp_tc_loc;
tc_trial_avrg_keep_allCond{id}=temp_tc_allCond;
%per day; each cell contains the average tc for each cell at that individual cell's preferred orientation and contrast
rect_tc_trial_avrg_keep{1,iCon,id}=rect_tc_trial_avrg; %rectified version of above
pref_responses_stat{id} = temp_pref_responses_stat;
pref_responses_loc{id} = temp_pref_responses_loc;
pref_responses_allCond{id} = temp_pref_responses_allCond;

end

trialCountTable = table([mean(trialCounts{1,2});std(trialCounts{1,2})],[mean(trialCounts{2,2});std(trialCounts{2,2})],[mean(trialCounts{1,1});std(trialCounts{1,1})],[mean(trialCounts{2,1});std(trialCounts{2,1})],'VariableNames',{'pre stat' 'pre loc' 'post stat' 'post loc'}, 'RowNames',{'Mean'  'std'})
writetable(trialCountTable,fullfile(fn_multi,'trialCounts.csv'),'WriteRowNames',true)

%%
% all contrasts, stationary
tc_trial_avrg_keep_allCon_stat=cell(1,nd);
pref_responses_allCon_stat = cell(2,nd);

for id = 1:nd
    pref_responses_temp=nan(1,nKeep);
    tc_trial_avrg_temp=nan((nOn+nOff),nKeep);
    mean_resp_temp=nan(nKeep,1);
    for i=1:nKeep
        
            temp_TCs=data_trial_keep{id}(:,:,i); %only pulling from dfof data of keep cells
            tDir=tDir_match{id}(1:nTrials(id));
            %identify the trials where ori = pref ori
            temp_dir= pref_dir_keep{id}(i); %find the preferred ori of this cell and convert to degrees
            dir_inds = find(tDir==temp_dir); %these are the trials at that ori
            inds=intersect(find(~RIx{id}),dir_inds);

            pref_responses_temp(i)=nanmean(nanmean(temp_TCs(stimStart:stimEnd,inds),1),2);
            pref_sd_temp(i)=std(nanmean(temp_TCs(stimStart:stimEnd,inds),1));
            tc_trial_avrg_temp(:,i)=nanmean(temp_TCs(:,inds),2);
        

    end
    
pref_responses_allCon_stat{1,id}=pref_responses_temp;
pref_responses_allCon_stat{2,id}=pref_sd_temp;
tc_trial_avrg_keep_allCon_stat{id}=tc_trial_avrg_temp; %this is a cell array with one cell 
%per day; each cell contains the average tc for each cell at that individual cell's preferred orientation averaged over all contrast
end

% all contrasts, running
% 
tc_trial_avrg_keep_allCon_loc=cell(1,nd);
pref_responses_allCon_loc = cell(2,nd);

for id = 1:nd
    pref_responses_temp=nan(1,nKeep);
    tc_trial_avrg_temp=nan((nOn+nOff),nKeep);
    mean_resp_temp=nan(nKeep,1);
    for i=1:nKeep
        
            temp_TCs=data_trial_keep{id}(:,:,i); %only pulling from dfof data of keep cells
            tDir=tDir_match{id}(1:nTrials(id));
            %identify the trials where ori = pref ori
            temp_dir= pref_dir_keep{id}(i); %find the preferred ori of this cell and convert to degrees
            dir_inds = find(tDir==temp_dir); %these are the trials at that ori
            inds=intersect(find(RIx{id}),dir_inds);

            pref_responses_temp(i)=nanmean(nanmean(temp_TCs(stimStart:stimEnd,inds),1),2);
            pref_sd_temp(i)=std(nanmean(temp_TCs(stimStart:stimEnd,inds),1));
            tc_trial_avrg_temp(:,i)=nanmean(temp_TCs(:,inds),2);
        

    end
    
pref_responses_allCon_loc{1,id}=pref_responses_temp;
pref_responses_allCon_loc{2,id}=pref_sd_temp;
tc_trial_avrg_keep_allCon_loc{id}=tc_trial_avrg_temp; %this is a cell array with one cell 
%per day; each cell contains the average tc for each cell at that individual cell's preferred orientation averaged over all contrast
end




clear tc_trial_avrg temp_trials con_inds temp_con dir_inds temp_dir mean_resp_temp temp_TCs
%%

red_keep_logical = zeros(1,nKeep);
for i = 1:length(red_ind_keep)
   red_keep_logical(red_ind_keep(i))=1;
end
green_keep_logical = ~red_keep_logical;

% msubset of full timecourses for keep cells only - this is not shifted and
% padded
fullTC_keep=cell(1,nd);

for id = 1:nd
    fullTC_keep{id} = cellTCs_match_OG{id}(:,keep_cells);
end


%sig_diff will be a logical vector indicating which cells were
%significantly modulated at each contrast, in terms of the call's response
%to its preferred orientation
% sig_diff = cell(1,nCon);
% for iCon = 1:nCon
%     for i = 1:nKeep
%         sig_diff{iCon}(i)=ttest2(pref_responses{iCon,i,1},pref_responses{iCon,i,2});
%     end
% end


explanation1 = 'tc_trial_keep contains the timecourses for all "keep" cells for each day. The tOri_match and tCon_match data structures can be used to find trials of particular stim conditions within this. tc_trial_avrg_keep only has the timecourses averaged over tirals for each cell at its preferred orientation and at each contrast.';
save(fullfile(fn_multi,'tc_keep.mat'),'explanation1','fullTC_keep','pref_responses_stat','pref_responses_loc','resp_keep','tc_trial_avrg_keep_allCond','pref_responses_allCond','tc_trial_avrg_keep_allCon_stat','pref_responses_allCon_stat','tc_trial_avrg_keep_allCon_loc','pref_responses_allCon_loc', 'pref_con_keep','pref_dir_keep','tDir_match','tOri_match','tCon_match','data_trial_keep','nTrials','tc_trial_avrg_stat','tc_trial_avrg_loc', 'green_keep_logical', 'red_keep_logical','green_ind_keep', 'red_ind_keep','stimStart')


%% make and save response matrix for keep cells
% this is averaging over stationry and running trials

data_resp_keep = cell(1,nd);
resp_max_keep = cell(1,nd);


for id = 1:nd
  data_resp_keep{id}=data_resp_match{id}(keep_cells,:,:,:);
  resp_max_keep{id} = squeeze(max(data_resp_keep{id}(:,:,:,1),[],2));
  
end



resp_max_keep_rect = resp_max_keep;
for id = 1:nd
    resp_max_keep_rect{id}(find(resp_max_keep_rect{id}<0))=0;
 
end

dfof_max_diff = (resp_max_keep_rect{1}-resp_max_keep_rect{2})./(resp_max_keep_rect{1}+resp_max_keep_rect{2}); % (post-pre)/(post+pre), nCell X nCon
dfof_max_diff_raw = (resp_max_keep{1}-resp_max_keep{2});

%make a data frame for the keep cells only
data_con_resp_keep = cell(1,nd);
data_dir_resp_keep = cell(1,nd);
for id = 1:nd
    data_con_resp_keep{id} = data_con_resp_match{id}(keep_cells,:);   
    data_dir_resp_keep{id} = data_dir_resp_match{id}(keep_cells,:);  
end

explanation2 = 'data_resp_keep gives the df/f averaged over the full response window for all conditions in the form nCells X nOris X nCons X mean vs. std. resp_max_keep gives the df/f averaged over the full stim period for each cell at the preferred ori only, with a dimension for each contrast';
save(fullfile(fn_multi,'resp_keep.mat'),'explanation2','data_resp_keep','resp_max_keep','dfof_max_diff','dfof_max_diff_raw','data_con_resp_keep','data_dir_resp_keep')
%% making mask maps for various measurements
%show masks
%get masks of matched cells
mask_match = cell(1,nd);
mask_match{1}= zeros(size(corrmap{1}));
mask_match{2}=masks{2}; %the second cell in the "masks" array already is only for matched cells
for i = 1:size(match_ind,2)
   ind = match_ind(i);
   temp_mask_inds = find(masks{1}==ind);
   mask_match{1}(temp_mask_inds)=i;
end

figure;
imagesc(corrmap{3});
colormap gray
%caxis([0.05 .3])
title('average FOV reference day');
hold on
bound = cell2mat(bwboundaries(mask_match{1}(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','b','MarkerSize',2);
bound = cell2mat(bwboundaries(mask_match{2}(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','g','MarkerSize',2);
hold off
print(fullfile(fn_multi,'matchCells.pdf'),'-dpdf');



keep_masks = zeros(size(corrmap{1}));
keep_green_masks = zeros(size(corrmap{1}));
keep_red_masks = zeros(size(corrmap{1}));
keep_masks_fract_change_red = zeros(size(corrmap{1}));
keep_masks_fract_change_green = zeros(size(corrmap{1}));
keep_masks_raw_change_red = zeros(size(corrmap{1}));
keep_masks_raw_change_green = zeros(size(corrmap{1}));
keep_masks_d1_red = zeros(size(corrmap{1}));

for i = 1:length(keep_cells)
   ind = keep_cells(i);
   temp_mask_inds = find(masks{2}==ind); %pulling from the masks of matched cells from the baseline day
   keep_masks(temp_mask_inds)=i;
   
end



%I am converting these to be labelled by their position in the keep cell
%index

for i = 1:length(keep_cells)
  temp_mask_inds = find(keep_masks==i);
   if ismember(i,red_ind_keep)
       keep_red_masks(temp_mask_inds)=i;
       keep_masks_fract_change_red(temp_mask_inds) = dfof_max_diff(i,(size(dfof_max_diff,2)));
       keep_masks_raw_change_red(temp_mask_inds) = dfof_max_diff_raw(i,(size(dfof_max_diff,2)));
       keep_masks_d1_red(temp_mask_inds)=resp_max_keep{2}(i,size(resp_max_keep{2},2));
   else
       keep_green_masks(temp_mask_inds)=i;
       keep_masks_fract_change_green(temp_mask_inds) = dfof_max_diff(i,(size(dfof_max_diff,2)));
       keep_masks_raw_change_green(temp_mask_inds) = dfof_max_diff_raw(i,(size(dfof_max_diff,2)));
   end
end

figure;
imagesc(fov_red{3});
colormap gray
caxis([10 100])
title('matched red cells');
hold on
bound = cell2mat(bwboundaries(keep_red_masks));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
hold off
print(fullfile(fn_multi,'matchRedCells.pdf'),'-dpdf');

save(fullfile(fn_multi,'mask_measuremens.mat'),'keep_masks','keep_red_masks','keep_masks_fract_change_red','keep_masks_raw_change_red','keep_masks_d1_red','keep_green_masks','keep_masks_fract_change_green','keep_masks_raw_change_green')
%% extract stationary and running trials seperately
%this is averaging over all dirs, I don't think I need this anymore
% 
% 
% locTCs = cell(1,nd);
% statTCs = cell(1,nd);
% 
% locResp = cell(1,nd);
% statResp = cell(1,nd);
% locCounts=cell(1,nd);
% 
% for id = 1:nd
%     %make this the average for running state for each contrast
%     locTrials=nan((nOn+nOff),nKeep,nCon);
%     statTrials=nan((nOn+nOff),nKeep,nCon);
%     for iCon = 1:nCon
%         tCon = tCon_match{id}(:,1:nTrials(id));
%         ind_con = find(tCon == cons(iCon));
%         loc=squeeze(nanmean(data_trial_keep{id}(:,intersect(find(RIx{id}), ind_con),:),2)); %average of locomotion at this contrast trials for each cell
%         stat=squeeze(nanmean(data_trial_keep{id}(:,intersect(find(~RIx{id}),ind_con),:),2)); %average of stationary  at this contrast trials for each cell
% 
%         locTCs{id}(:,:,iCon)=loc;
%         statTCs{id}(:,:,iCon)=stat;
%         locCounts{id}(:,iCon)=[length(intersect(find(RIx{id}), ind_con)),length(intersect(find(~RIx{id}),ind_con))];
% 
%         locResp{id}(:,iCon)=nanmean(loc(stimStart:stimStart+nOn,:),1);
%         statResp{id}(:,iCon)=nanmean(stat(stimStart:stimStart+nOn,:),1);
%         clear loc stat ind_con tCon
%     end
% end
% 
% 




%for each day, extract the LMI for each cell, calculated at the preferred
%directon

LMI = cell(1,nd);
 %for the loc-stat/loc+stat version

for id = 1:nd
    for iCon=1:nCon
            locRectified = pref_responses_loc{id}(:,iCon);
            locRectified(find(locRectified<0))=0;
            statRectified = pref_responses_stat{id}(:,iCon);
            statRectified(find(statRectified<0))=0;
            LMI{id}(:,iCon)=(locRectified-statRectified)./(locRectified+statRectified);
        
    end
end

%for loc/stat version
% for id = 1:nd
%     for iCon=1:nCon
%             locTemp = pref_responses_loc{id}(:,iCon);
%             statTemp = pref_responses_stat{id}(:,iCon);
%             
%             LMI{id}(:,iCon)=(locTemp)./(statTemp);
%         
%     end
% end



save(fullfile(fn_multi,'locomotion.mat'),'LMI','RIx','wheel_tc','wheel_speed')
% %% comparing F and df/f for HT+ and HT-
% figure;
% subplot(1,2,1)
% dfof_pre = squeeze( nanmean(nanmean(data_dfof_trial_match{2}(resp_win,:,:),1),2)); %gets the average dfof during all stim on perdiods
% f_pre = nanmean(cellTCs_match{2},1)'; %gets the average fluorescence for the entire recording time
% scatter(f_pre,dfof_pre,'k')%plots all cells
% hold on
% scatter(f_pre(red_ind_match),dfof_pre(red_ind_match),'r') %replots only HT+ cells, now colored red
% hold on
% plot(mean(f_pre(red_ind_match)),mean(dfof_pre(red_ind_match)),'r.','MarkerSize',25);
% hold on
% plot(mean(f_pre(~red_ind_match)),mean(dfof_pre(~red_ind_match)),'k.','MarkerSize',25);
% xlabel({'Mean F for the';'entire timecourse'})
% ylabel('Mean dF/F for all stim on periods');
% xlim([0 20000])
% ylim([-.2 .5])
% title('pre-DART');
% axis square
% hold off
% 
% subplot(1,2,2)
% dfof_post = squeeze( nanmean(nanmean(data_dfof_trial_match{1}(resp_win,:,:),1),2)); %gets the average dfof during all stim on perdiods
% f_post = nanmean(cellTCs_match{1},1)'; %gets the average fluorescence for the entire recording time
% scatter(f_post,dfof_post,'k')%plots all cells
% hold on
% scatter(f_post(red_ind_match),dfof_post(red_ind_match),'r') %replots only HT+ cells, now colored red
% hold on
% plot(mean(f_post(red_ind_match)),mean(dfof_post(red_ind_match)),'r.','MarkerSize',25);
% hold on
% plot(mean(f_post(~red_ind_match)),mean(dfof_post(~red_ind_match)),'k.','MarkerSize',25);
% xlabel({'Mean F for the';'entire timecourse'})
% ylabel('Mean dF/F for all stim on periods');
% xlim([0 20000])
% ylim([-.2 .5])
% title('post-DART');
% axis square
% hold off
% print(fullfile(fn_multi, 'F_vs_dFoF.pdf'),'-dpdf');