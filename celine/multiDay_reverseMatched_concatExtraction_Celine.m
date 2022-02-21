clear all; clear global; close all
clc
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc =  behavConstsDART; %directories
eval(ds);

day_id_list = [142; %enter post-DART day
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
load(fullfile(fn_multi,'multiday_alignment.mat'))
load(fullfile(fn_multi,'input.mat'))
frame_rate = input.frameImagingRateMs;
%% finding red fluorescence level
%  
% red_fluor_all = cell(1,nd);
% 
% for i = 1:nd
% mouse = expt(day_id(i)).mouse;
% date = expt(day_id(i)).date;
% 
% imgFolder = expt(day_id(i)).contrastxori_runs{1};
% fn = fullfile(rc.achAnalysis,mouse,date,imgFolder);
% cd(fn);
% load('redImage.mat');  
% load('mask_cell.mat');
%     
%     
% % cell_stats=regionprops(mask_cell_red);
% % figure; imagesc(redChImg), colormap gray; caxis([200 1000]);
% % hold on
% % bound = cell2mat(bwboundaries(mask_cell_red(:,:,1)));
% % plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',.5); hold on;
% % for iC = 1:max(max(mask_cell_red))
% %     text(cell_stats(iC).Centroid(1), cell_stats(iC).Centroid(2), num2str(iC), 'Color', 'red',...
% %             'Fontsize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
% %     
% % end
% 
% red_fluor_all{i} = stackGetTimeCourses(redChImg, mask_cell);
% end
% red_fluor_match_d1=red_fluor_all{1}(:,match_ind);
% z_red_fluor_d1=zscore(red_fluor_match_d1);


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
cons = unique(tCon_match{1});
nOri = length(oris);
nCon = length(cons);


%% convert to trials

%change this to use a padded array, where I add zeros at the end. test=padarray(cellTCs_match{1},30,0,'post');
stimStart = (nOff/2)+1; %this indicates both the perdiod to trim off the start and the stim on period after trimming
stimEnd=stimStart+nOn-1;

cellTCs_match_OG = cellTCs_match;
data_dfof_trial_match = cell(1,nd); %make an empty array that is 1 by however many days there are (1X2 usually)

fractTimeActive_match = cell(1,nd);
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
    fractTimeActive_match{id} = zeros(1,nCells);
    meansub_match = cellTCs_match{id}-mean(cellTCs_match{id},1);
    cellstd_match = std(meansub_match,[],1);
    for iCell = 1:nCells
        fractTimeActive_match{id}(:,iCell) = length(find(meansub_match(:,iCell)>3.*cellstd_match(1,iCell)))./nFrames;
    end
    
end
clear data_trial_match data_f_match

% MUST USE NANMEAN INSTEAD OF MEAN MOVING FORWARD SINCE I SET THE PADDING
% VALUES TO NAN

%% find significant responses and preferred stimuli

resp_win = stimStart:stimEnd;
base_win = 1: stimStart-1;

%make cell arrays to keep everything in
data_resp_match = cell(1,nd);
h_match = cell(1,nd);
p_match = cell(1,nd);
h_all_match = cell(1,nd);
resp_match = cell(1,nd);
pref_ori_match = cell(1,nd);
pref_con_match = cell(1,nd);
data_ori_resp_match = cell(1,nd);
data_con_resp_match = cell(1,nd);

%rect_dfof_trial=data_dfof_trial_match;
%rect_dfof_trial(find(rect_dfof_trial<0))=0;
% for each day
for id = 1:nd
    data_resp = zeros(nCells, nOri, nCon,2);
    h = zeros(nCells, nOri, nCon);
    p = zeros(nCells, nOri, nCon);
    tCon = tCon_match{id}(:,1:nTrials(id));
    tOri = tOri_match{id}(:,1:nTrials(id));
    tDir = tDir_match{id}(:,1:nTrials(id));
    data_dfof_trial = data_dfof_trial_match{id}; 

    for iOri = 1:nOri
        ind_ori = find(tOri == oris(iOri));
        for iCon = 1:nCon
            ind_con = find(tCon == cons(iCon));
            ind = intersect(ind_ori,ind_con); %for every orientation and then every contrast, find trials with that con/ori combination
            data_resp(:,iOri,iCon,1) = squeeze(nanmean(nanmean(data_dfof_trial(resp_win,ind,:),1),2));
            data_resp(:,iOri,iCon,2) = squeeze(std(nanmean(data_dfof_trial(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
            [h(:,iOri,iCon), p(:,iOri,iCon)] = ttest(nanmean(data_dfof_trial(resp_win,ind,:),1), nanmean(data_dfof_trial(base_win,ind,:),1),'dim',2,'tail','right','alpha',0.05./(nOri*nCon.*3-1));
        end
    end

    h_all = sum(sum(h,2),3);

    resp=logical(h_all);
    
    pref_ori = zeros(1,nCells);
    pref_con = zeros(1,nCells);
    data_ori_resp = zeros(nCells,nOri); %at pref con
    data_con_resp = zeros(nCells,nCon); %at pref ori
    data_orth_resp=zeros(nCells,1);
    %I want to pull out the responses for each cell at it's preferred orientations, for
    %all contrasts, and at it's preferred contrast, for all orientations
    for iCell = 1:nCells
          [max_val, pref_ori(1,iCell)] = max(nanmean(data_resp(iCell,:,:,1),3),[],2);
          [max_val_con, pref_con(1,iCell)] = max(squeeze(nanmean(data_resp(iCell,:,:,1),2))',[],2);
          data_ori_resp(iCell,:)=data_resp(iCell,:,pref_con(iCell),1);
          data_con_resp(iCell,:)=nanmean(data_resp(iCell,:,:,1),2);
          
      if pref_ori(iCell)<= nOri/2
          orth_ori(iCell)=pref_ori(iCell)+2;
      elseif pref_ori(iCell)>nOri/2
          orth_ori(iCell)=pref_ori(iCell)-2;
      end
    data_orth_resp(iCell,:)=nanmean(data_resp(iCell,orth_ori(iCell),:,1));
    end

 %then put into a cell matrix for the two days
data_resp_match{id} = data_resp;
h_match{id} = h;
p_match{id} = p;
h_all_match{id} = h_all;
resp_match{id} = resp;
pref_ori_match{id} = pref_ori;
pref_con_match{id} = pref_con;
data_ori_resp_match{id} = data_ori_resp;
data_con_resp_match{id} = data_con_resp;
data_orth_resp_match{id}=data_orth_resp;
 
end
clear data_resp p h_all resp pref_ori pref_con data_ori_resp data_con_resp data_dfof_trial tCon tOri tDir data_orth_resp
 
%% get basic counts 
red_ind_match_list = find(red_ind_match==1);

%this is a list of indices of all the green cells
green_ind_match = ~(red_ind_match);
green_ind_match_list = find(green_ind_match);


%creating the arrays for red cells
red_match_respd1 = intersect(red_ind_match_list,find(resp_match{2}==1));%this is for reverse matching
red_match_respd2 = intersect(red_ind_match_list,find(resp_match{1}==1));

green_match_respd1 = intersect(green_ind_match_list,find(resp_match{2}==1));
green_match_respd2 = intersect(green_ind_match_list,find(resp_match{1}==1));

%matched
nGreen_keep = length(green_ind_match_list); %how many matched green cells
nGreen_keep_respd1 = length(green_match_respd1); %how many of those responded on d1
nGreen_keep_respd2 = length(green_match_respd2);%how many of those responded on d2
nRed_match = length(red_ind_match_list);%how many red cells matched
nRed_match_respd1 = length(red_match_respd1); %how many of the matched red cells responded on d1
nRed_match_respd2 = length(red_match_respd2); %how many of the matched red cells responded on d1

%% extract data for cells I want to keep

%first narrow down to the cells in question - find the green cells that
%were responsive on at least one day

resp_green_either = union(green_match_respd1,green_match_respd2); %these are the green cells I will include from now on


%these are the cells I will include from now on
keep_cells = union(resp_green_either,red_ind_match_list);

nKeep = length(keep_cells)

%making a list of the day 1 indices for the cells I want to keep
% keep_d1 = match_ind(:,keep_cells);
% 
[red_ind_match_list,red_ind_keep] = intersect(keep_cells,red_ind_match_list,'stable'); %find the indices *within keep_cells* that are red
[resp_green_either,green_ind_keep] = intersect(keep_cells,resp_green_either,'stable'); %same for green

%for now I will keep both sets of indices - if the scripts runs slowly I
%can clear these


%% counts for keep cells

green_match_respd1 = intersect(resp_green_either,find(resp_match{2}==1)); %this is for reverse matching
green_match_respd2 = intersect(resp_green_either,find(resp_match{1}==1));

%matched
nGreen_keep = length(resp_green_either); %how many matched green cells
nGreen_keep_respd1 = length(green_match_respd1); %how many of those responded on d1
nGreen_keep_respd2 = length(green_match_respd2);%how many of those responded on d2

% make table of values
countsTable = table([nGreen_keep;nRed_match],[nGreen_keep_respd1;nRed_match_respd1],[nGreen_keep_respd2;nRed_match_respd2],'VariableNames',{'Keep' 'Responsive pre' 'Responsive post'}, 'RowNames',{'Pyramidal cells'  'HT+ cells'})
writetable(countsTable,fullfile(fn_multi,'match_counts.csv'),'WriteRowNames',true)
clear  nGreen_match_respd1 nGreen_match_respd2  nRed_match_respd1 nRed_match_respd2

%% make a data structure subsets for only the keep cells
data_trial_keep=cell(1,nd);
pref_ori_keep=cell(1,nd);
pref_con_keep=cell(1,nd);

for id = 1:nd
    data_trial_keep{id} = data_dfof_trial_match{id}(:,:,keep_cells);
    pref_ori_keep{id} = pref_ori_match{id}(:,keep_cells);
    pref_ori_keep{id}=oris(pref_ori_keep{id});
    pref_con_keep{id} = pref_con_match{id}(:,keep_cells);
    pref_con_keep{id}=cons(pref_con_keep{id});
end


conTable = table([mean(pref_con_keep{2}(green_ind_keep));mean(pref_con_keep{2}(red_ind_keep))],[mean(pref_con_keep{1}(green_ind_keep));mean(pref_con_keep{1}(red_ind_keep))],'VariableNames',{'mean pref con pre' 'mean pref con post'}, 'RowNames',{'Pyramidal cells'  'HT+ cells'})
writetable(conTable,fullfile(fn_multi,'conPref.csv'),'WriteRowNames',true)

%% narrow down to the stimuli preferred for each cell each day
%we will get one tc per cell per contrast. This represents that cell's tc
%averaged over trials at the preferred orientation, at each contrast

tc_trial_avrg_keep=cell(1,nd);
rect_tc_trial_avrg_keep=cell(1,nd);
resp_prefStim_keep=cell(1,nd);

for id = 1:nd
    
    tc_trial_avrg=nan((nOn+nOff),nKeep,nCon);
    mean_resp_temp=nan(nKeep,nCon,1);
    for i=1:nKeep
        for iCon = 1:nCon
            temp_TCs=data_trial_keep{id}(:,:,i); %only pulling from dfof data of keep cells
            tCon=tCon_match{id}(1:nTrials(id));
            tDir=tDir_match{id}(1:nTrials(id));
            %identify the trials where ori = pref ori
            temp_ori= pref_ori_keep{2}(i); %find the preferred ori of this cell and convert to degrees
            ori_inds = find(tDir==temp_ori); %these are the trials at that ori

            con_inds=find(tCon==cons(iCon));
            temp_trials = intersect(ori_inds, con_inds); %preferred ori for this cell, looping through all cons

            tc_trial_avrg(:,i,iCon)=nanmean(temp_TCs(:,temp_trials),2);
            rect_tc_trial_avrg=tc_trial_avrg;
            rect_tc_trial_avrg(rect_tc_trial_avrg<0)=0;
            mean_resp_temp(i,iCon) = nanmean(tc_trial_avrg(stimStart:stimEnd,i,iCon));
        end

    end
    
    
tc_trial_avrg_keep{id}=tc_trial_avrg; %this is a cell array with one cell 
%per day; each cell contains the average tc for each cell at that individual cell's preferred orientation and contrast
rect_tc_trial_avrg_keep{1,iCon,id}=rect_tc_trial_avrg; %rectified version of above
resp_prefStim_keep{id}=mean_resp_temp; %a single column array with a value for each cell that represents the mean df/f of the response window for preferred stimuli
end


clear tc_trial_avrg temp_trials con_inds temp_con ori_inds temp_ori mean_resp_temp temp_TCs

red_keep_logical = zeros(1,nKeep);
for i = 1:length(red_ind_keep)
   red_keep_logical(red_ind_keep(i))=1;
end
green_keep_logical = ~red_keep_logical;


save(fullfile(fn_multi,'tc_keep.mat'),'nTrials','tc_trial_avrg_keep', 'green_keep_logical', 'red_keep_logical','green_ind_keep', 'red_ind_keep','stimStart')

%% make and save response matrix for keep cells

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
data_ori_resp_keep = cell(1,nd);
for id = 1:nd
    data_con_resp_keep{id} = data_con_resp_match{id}(keep_cells,:);   
    data_ori_resp_keep{id} = data_ori_resp_match{id}(keep_cells,:);  
end
save(fullfile(fn_multi,'resp_keep.mat'),'data_resp_keep','resp_max_keep','dfof_max_diff','dfof_max_diff_raw','data_con_resp_keep','data_ori_resp_keep')
%% making mask maps for various measurements
%show masks
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
   temp_mask_inds = find(masks{2}==ind);
   keep_masks(temp_mask_inds)=i;
   
end
%I am converting these to be labelled by their position in the keep cell
%index

for i = 1:length(keep_cells)
  temp_mask_inds = find(keep_masks==i);
   if ismember(i,red_ind_keep)
       keep_red_masks(temp_mask_inds)=i;
       keep_masks_fract_change_red(temp_mask_inds) = dfof_max_diff(i,2);
       keep_masks_raw_change_red(temp_mask_inds) = dfof_max_diff_raw(i,2);
       keep_masks_d1_red(temp_mask_inds)=resp_max_keep{2}(i,2);
   else
       keep_green_masks(temp_mask_inds)=i;
       keep_masks_fract_change_green(temp_mask_inds) = dfof_max_diff(i,2);
       keep_masks_raw_change_green(temp_mask_inds) = dfof_max_diff_raw(i,2);
   end
end

save(fullfile(fn_multi,'mask_measuremens.mat'),'keep_red_masks','keep_masks_fract_change_red','keep_masks_raw_change_red','keep_masks_d1_red','keep_green_masks','keep_masks_fract_change_green','keep_masks_raw_change_green')