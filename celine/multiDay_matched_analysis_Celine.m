clear all; clear global; close all
clc
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsAV; %directories
eval(ds)

day_id(2) = 123;
day_id(1) = expt(day_id(2)).multiday_matchdays;
nd = size(day_id,2);
day_id

mouse = expt(day_id(1)).mouse;

fnout = fullfile(rc.celineAnalysis,mouse);
if expt(day_id(2)).multiday_timesincedrug_hours>0
    dart_str = [expt(day_id(2)).drug '_' num2str(expt(day_id(2)).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end


fn_multi = fullfile(rc.celineAnalysis,mouse,['multiday_' dart_str]);

cd(fn_multi)
load(fullfile(fn_multi,'timecourses.mat'))
load(fullfile(fn_multi,'multiday_alignment.mat'))
load(fullfile(fn_multi,'input.mat'))

%% finding red fluorescent level
 
red_fluor_all = cell(1,nd);

for i = 1:nd
mouse = expt(day_id(i)).mouse;
date = expt(day_id(i)).date;

imgFolder = expt(day_id(i)).contrastxori_runs{1};
fn = fullfile(rc.celineAnalysis,mouse,date,imgFolder);
cd(fn);
load('redImage.mat');  
load('mask_cell.mat');
    
    
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

red_fluor_all{i} = stackGetTimeCourses(redChImg, mask_cell);
end
red_fluor_match_d1=red_fluor_all{1}(:,match_ind);
z_red_fluor_d1=zscore(red_fluor_match_d1);


%% stimulus props

%tells the contrast, direction and orientation for each trial each day
tCon_match = cell(1,nd);
tDir_match = cell(1,nd);
tOri_match = cell(1,nd);

%find the contrasts, directions and orientations for each day
for id = 1:nd
    tCon_match{id} = celleqel2mat_padded(input(id).tGratingContrast);
    tDir_match{id} = celleqel2mat_padded(input(id).tGratingDirectionDeg);
    tOri_match{id} = tDir_match{id};
    tOri_match{id}(find(tDir_match{id}>=180)) = tDir_match{id}(find(tDir_match{id}>=180))-180;
end
oris = unique(tOri_match{1});
cons = unique(tCon_match{1});
nOri = length(oris);
nCon = length(cons);

nOn = input(1).nScansOn;
nOff = input(1).nScansOff;

%% convert to trials

%change this to use a padded array, where I add zeros at the end. test=padarray(cellTCs_match{1},30,0,'post');
stimStart = (nOff/2)+1; %this indicates both the perdiod to trim off the start and the stim on period after trimming
stimEnd=stimStart+nOn-1;

cellTCs_match_OG = cellTCs_match;
data_dfof_trial_match = cell(1,nd); %make an empty array that is 1 by however many days there are (1X2 usually)

fractTimeActive_match = cell(1,nd);
for id = 1:nd %cycle through days
    
   nTrials = length(tDir_match{id}); %use the list of direction by trial to figure out how many trials there are
   %currently the way I center the stim on period requires me to cut out
   %one trial, hence the -1
   
    %for matched cell
    nCells = size(cellTCs_match{id},2);
    fractTimeActive_match{id} = zeros(1,nCells);
    %I will trim 30 frames of the start and add 30 frames of padding to the
    %end (padding is zeros)
    cellTCs_match{id} = cellTCs_match{id}(stimStart:size(cellTCs_match{id},1),:);
    cellTCs_match{id} = padarray(cellTCs_match{id},30,0,'post');
    nFrames = size(cellTCs_match{id},1);
    data_trial_match = reshape(cellTCs_match{id},[nOn+nOff nTrials nCells]);
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

%% find significant responses and preferred stimuli

resp_win = stimStart:stimEnd;
base_win = 1: (nOff/2);

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

rect_dfof_trial=data_dfof_trial_match;
%rect_dfof_trial(find(rect_dfof_trial<0))=0;
% for each day
for id = 1:nd
    data_resp = zeros(nCells, nOri, nCon,2);
    h = zeros(nCells, nOri, nCon);
    p = zeros(nCells, nOri, nCon);
    tCon = tCon_match{id}(:,1:nTrials);
    tOri = tOri_match{id}(:,1:nTrials);
    tDir = tDir_match{id}(:,1:nTrials);
    data_dfof_trial = rect_dfof_trial{id}; %will use rectified data

    for iOri = 1:nOri
        ind_ori = find(tOri == oris(iOri));
        for iCon = 1:nCon
            ind_con = find(tCon == cons(iCon));
            ind = intersect(ind_ori,ind_con); %for every orientation and then every contrast, find trials with that con/ori combination
            data_resp(:,iOri,iCon,1) = squeeze(mean(mean(data_dfof_trial(resp_win,ind,:),1),2));
            data_resp(:,iOri,iCon,2) = squeeze(std(mean(data_dfof_trial(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
            [h(:,iOri,iCon), p(:,iOri,iCon)] = ttest(mean(data_dfof_trial(resp_win,ind,:),1), mean(data_dfof_trial(base_win,ind,:),1),'dim',2,'tail','right','alpha',0.05./(nOri.*3-1));
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
          [max_val, pref_ori(1,iCell)] = max(mean(data_resp(iCell,:,:,1),3),[],2);
          [max_val_con, pref_con(1,iCell)] = max(squeeze(mean(data_resp(iCell,:,:,1),2))',[],2);
          data_ori_resp(iCell,:)=data_resp(iCell,:,pref_con(iCell),1);
          data_con_resp(iCell,:)=data_resp(iCell,pref_ori(iCell),:,1);
          
      if pref_ori(iCell)<= 4
          orth_ori(iCell)=pref_ori(iCell)+4;
      elseif pref_ori(iCell)>4
          orth_ori(iCell)=pref_ori(iCell)-4;
      end
    data_orth_resp(iCell,:)=mean(data_resp(iCell,orth_ori(iCell),:,1));
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
clear data_resp h p h_all resp pref_ori pref_con data_ori_resp data_con_resp data_dfof_trial tCon tOri tDir data_orth_resp
 
%% get basic counts 
red_ind_match_list = find(red_ind_match==1);

%this is a list of indices of all the green cells
green_ind_match = ~(red_ind_match);
green_ind_match_list = find(green_ind_match);


%creating the arrays for red cells
red_match_respd1 = intersect(red_ind_match_list,find(resp_match{1}==1));
red_match_respd2 = intersect(red_ind_match_list,find(resp_match{2}==1));

green_match_respd1 = intersect(green_ind_match_list,find(resp_match{1}==1));
green_match_respd2 = intersect(green_ind_match_list,find(resp_match{2}==1));

%matched
nGreen_match = length(green_ind_match_list); %how many matched green cells
nGreen_match_respd1 = length(green_match_respd1); %how many of those responded on d1
nGreen_match_respd2 = length(green_match_respd2);%how many of those responded on d2
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

green_match_respd1 = intersect(resp_green_either,find(resp_match{1}==1));
green_match_respd2 = intersect(resp_green_either,find(resp_match{2}==1));


%matched
nGreen_match = length(resp_green_either); %how many matched green cells
nGreen_match_respd1 = length(green_match_respd1); %how many of those responded on d1
nGreen_match_respd2 = length(green_match_respd2);%how many of those responded on d2

% make table of values
countsTable = table([nGreen_match;nRed_match],[nGreen_match_respd1;nRed_match_respd1],[nGreen_match_respd2;nRed_match_respd2],'VariableNames',{'Keep' 'Responsive day 1' 'Responsive day 2'}, 'RowNames',{'Pyramidal cells'  'HT+ cells'})
writetable(countsTable,fullfile(fn_multi,'match_counts.csv'),'WriteRowNames',true)
clear nGreen_match nGreen_match_respd1 nGreen_match_respd2 nRed_match nRed_match_respd1 nRed_match_respd2

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

%narrow down to the stimuli preferred for each cell each day - note that
%the contrast and orientation preferred for a given cell may differ across
%days
tc_trial_avrg_keep=cell(1,nd);
rect_tc_trial_avrg_keep=cell(1,nd);
resp_prefStim_keep=cell(1,nd);

for id = 1:nd
    tc_trial_avrg=nan((nOn+nOff),nKeep);
    mean_resp_temp=nan(nKeep,1);
    for i=1:nKeep
        temp_TCs=data_trial_keep{id}(:,:,i); %only pulling from dfof data of keep cells
        tCon=tCon_match{id}(1:nTrials);
        tDir=tDir_match{id}(1:nTrials);
        %identify the trials where ori = pref ori
        temp_ori= pref_ori_keep{id}(i); %find the preferred ori of this cell and convert to degrees
        ori_inds = find(tDir==temp_ori); %these are the trials at that ori
        temp_con = pref_con_keep{id}(i);%find the preferred contrast of this cell and convert to contrast value
        con_inds=find(tCon==temp_con);
        temp_trials = intersect(ori_inds, con_inds);
        tc_trial_avrg(:,i)=mean(temp_TCs(:,temp_trials),2);
        rect_tc_trial_avrg=tc_trial_avrg;
        rect_tc_trial_avrg(rect_tc_trial_avrg<0)=0;
        mean_resp_temp(i) = mean(tc_trial_avrg(stimStart:stimEnd,i));

    end
tc_trial_avrg_keep{id}=tc_trial_avrg; %this is a cell array with one cell 
%per day; each cell contains the average tc for each cell at that individual cell's preferred orientation and contrast
rect_tc_trial_avrg_keep{id}=rect_tc_trial_avrg; %rectified version of above
resp_prefStim_keep{id}=mean_resp_temp; %a single column array with a value for each cell that represents the mean df/f of the response window for preferred stimuli
end

clear tc_trial_avrg temp_trials con_inds temp_con ori_inds temp_ori mean_resp_temp temp_TCs
red_keep_logical = zeros(1,nKeep);
for i = 1:length(red_ind_keep)
   red_keep_logical(red_ind_keep(i))=1;
end
green_keep_logical = ~red_keep_logical;

z_red_fluor_keep = z_red_fluor_d1(keep_cells);

save(fullfile(fn_multi,'tc_keep.mat'),'tc_trial_avrg_keep', 'green_keep_logical', 'red_keep_logical')
%% prepare to plot the timecourses 

tc_green_avrg_match = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_match = cell(1,nd); %same for red
tc_green_se_match = cell(1,nd); %this will be the se across all green cells
tc_red_se_match = cell(1,nd); %same for red

for id = 1:nd
    tc_green_avrg_match{id}=mean(tc_trial_avrg_keep{id}(:,green_ind_keep),2);
    green_std=std(tc_trial_avrg_keep{id}(:,green_ind_keep),[],2);
    tc_green_se_match{id}=green_std/sqrt(length(green_ind_keep));
    
    tc_red_avrg_match{id}=mean(tc_trial_avrg_keep{id}(:,red_ind_keep),2);
    red_std=std(tc_trial_avrg_keep{id}(:,red_ind_keep),[],2);
    tc_red_se_match{id}=red_std/sqrt(length(red_ind_keep));
    
    clear green_std red_std
end

%% make figure with se shaded
figure
subplot(1,2,1) %for the first day

x=1:(size(tc_green_avrg_match{1},1));
x=(x-30)/15;
shadedErrorBar(x,tc_red_avrg_match{1},tc_red_se_match{1},'r');
ylim([-.02 .35]);
hold on
shadedErrorBar(x,tc_green_avrg_match{1},tc_green_se_match{1});
title('day 1')



subplot(1,2,2) %for the second day
shadedErrorBar(x,tc_red_avrg_match{2},tc_red_se_match{2},'r');
ylim([-.02 .35]);
hold on
shadedErrorBar(x,tc_green_avrg_match{2},tc_green_se_match{2});
title('day 2')


print(fullfile(fn_multi,['timecourses']),'-dpdf');
saveas(gcf,fullfile(fn_multi,[mouse '-' date 'tcPlot.jpg']));
%% make a plot of individual timecourses 
figure
subplot(2,2,1)
plot(tc_trial_avrg_keep{1}(:,green_ind_keep),'k')
ylim([-.1 .5]);
title('day 1')

subplot(2,2,2)
plot(tc_trial_avrg_keep{1}(:,red_ind_keep),'color',[.7 .05 .05])
ylim([-.1 .5]);
title('day 1')

subplot(2,2,3)
plot(tc_trial_avrg_keep{2}(:,green_ind_keep),'k')
ylim([-.1 .5]);
title('day 2')

subplot(2,2,4)
plot(tc_trial_avrg_keep{2}(:,red_ind_keep),'color',[.7 .05 .05])
ylim([-.1 .5]);
title('day 2')
print(fullfile(fn_multi,['indiv_timecourses']),'-dpdf');
%% makes a scatterplot of max df/f for day 1 vs day 2, and each subplot is one day
%this is for all cells I'm keeping, red and green


data_resp_keep = cell(1,nd);
resp_max_match = cell(1,nd);
resp_max_keep = cell(1,nd);

for id = 1:nd
    data_resp_keep{id}=data_resp_match{id}(keep_cells,:,:,:);
    resp_max_match{id} = squeeze(max(max(data_resp_match{id}(:,:,:,1),[],2),[],3));
   resp_max_keep{id} = squeeze(max(max(data_resp_keep{id}(:,:,:,1),[],2),[],3));
end


figure; movegui('center') 
subplot(1,2,1)
scatter(resp_max_keep{1}(green_ind_keep),resp_max_keep{2}(green_ind_keep),'k')
% hold on
% scatter(resp_max_keep{1}(red_ind_keep),resp_max_keep{2}(red_ind_keep),'MarkerEdgeColor',[.7 .05 .05])
% hold off
xlabel('D1- max dF/F')
ylabel('D2- max dF/F')
xlim([-.10 .5])
ylim([-.10 .5])
refline(1)
title('Max df/f for responsive HT- ')


subplot(1,2,2)
scatter(resp_max_keep{1}(red_ind_keep),resp_max_keep{2}(red_ind_keep),'MarkerEdgeColor',[.7 .05 .05])
% hold on
% scatter(resp_max_match{1}(red_ind_match_list),resp_max_match{2}(red_ind_match_list),'MarkerEdgeColor',[.7 .05 .05])
% hold off
xlabel('D1- max dF/F')
ylabel('D2- max dF/F')
xlim([-.10 .5])
ylim([-.10 .5])
refline(1)
title('Max df/f HT+')

print(fullfile(fn_multi,'maxResp_crossDay.pdf'),'-dpdf','-bestfit')
save(fullfile(fn_multi,'resp_keep.mat'),'data_resp_keep','resp_max_keep')
% extract the max df/f values for analysis
%% looking at change in dfof
for i = 1:nKeep
    dfof_max_diff(i) = (resp_max_keep{1}(i)-resp_max_keep{2}(i))/resp_max_keep{1}(i);
end

figure
x = [mean(dfof_max_diff(green_ind_keep)), mean(dfof_max_diff(red_ind_keep))];
y = [(std(dfof_max_diff(green_ind_keep)))/sqrt(length(green_ind_keep)), (std(dfof_max_diff(red_ind_keep)))/sqrt(length(red_ind_keep))];
%y = [std(dfof_max_diff(green_ind_keep)), std(dfof_max_diff(red_ind_keep))]
labs =categorical({'HT-','HT+'})
bar(labs,x)                
hold on
er = errorbar(labs,x,-y,y);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
%ylim([0 .2])
hold off
title('change in max dfof')

print(fullfile(fn_multi,'change_max_resp.pdf'),'-dpdf','-bestfit')


%% plot fractional change in max dfof vs. z-score of red fluorescence
wIn_red_z = zscore(red_fluor_match_d1(red_ind_keep));

figure
subplot(1,2,1)
scatter(z_red_fluor_keep(red_ind_keep),dfof_max_diff(red_ind_keep),'MarkerEdgeColor',[.7 .05 .05])
hold on
scatter(z_red_fluor_keep(green_ind_keep),dfof_max_diff(green_ind_keep),'k')
hold on
x=z_red_fluor_keep;
y=dfof_max_diff;
p = polyfit(x, y, 1);
px = [min(x) max(x)];
py = polyval(p, px);
plot(px, py, 'LineWidth', 2);
clear x y p px py
xlabel('red fluorescence z-score')
ylabel('fractional change max dfof')

subplot(1,2,2)
scatter(wIn_red_z,dfof_max_diff(red_ind_keep),'MarkerEdgeColor',[.7 .05 .05])
hold on
x=wIn_red_z;
y=dfof_max_diff(red_ind_keep);
p = polyfit(x, y, 1);
px = [min(x) max(x)];
py = polyval(p, px);
plot(px, py, 'LineWidth', 2);
clear x y p px py
xlabel('red fluorescence z-score within red cells')
ylabel('fractional change max dfof')
hold off

print(fullfile(fn_multi,'frac_change_vs_red.pdf'),'-dpdf','-bestfit')
%% plotting  ori response 
if nKeep<36
    [n n2] = subplotn(nKeep);
    tot = n.*n2;
else
    n = 4;
    n2 = 4;
    tot = 16;
end

% plot ori tuning at each contrast
figure;
movegui('center')
start = 1;

for iCell = 1:nKeep
    for id = 1:nd
        if start>tot
            figure; movegui('center')
            start = 1;
        end
        subplot(n,n2,start)

        for iCon = 1:nCon
            errorbar(oris, data_resp_keep{id}(iCell,:,iCon,1), data_resp_keep{id}(iCell,:,iCon,2),'-o')
            hold on
        end
        if ismember(iCell,red_ind_keep)
            extra_title='HT+';
        else
            extra_title='HT-';
        end
            start= start+1;
        ylim([-0.1 inf])
        title(['Day ' num2str(id) extra_title])
    end

end

% plot ori tuning at preferred contrast only
data_ori_resp_keep = cell(1,nd);
for id = 1:nd
    data_ori_resp_keep{id} = data_ori_resp_match{id}(keep_cells,:)    
end

figure;
movegui('center')
start = 1;

for iCell = 1:nKeep
        if start>tot
            figure; movegui('center')
            start = 1;
        end
        subplot(n,n2,start)
        for  id = 1:nd
            plot(oris, data_ori_resp_keep{id}(iCell,:))
            hold on
        end
        if ismember(iCell,red_ind_keep)
            extra_title='HT+';
        else
            extra_title='HT-';
        end
            start= start+1;
        ylim([-0.1 inf])
        title([extra_title])
        
end

%% plot observed pref ori 
%the orientation that gave the strongest response in the observed data

fig=figure; movegui('center') 
scatter(pref_ori_keep{1}(green_ind_keep),pref_ori_keep{2}(green_ind_keep),'k','jitter', 'on', 'jitterAmount', 5)
hold on
scatter(pref_ori_keep{1}(red_ind_keep),pref_ori_keep{2}(red_ind_keep),'MarkerEdgeColor',[.7 .05 .05],'jitter', 'on', 'jitterAmount', 5)
hold off
xlabel('D1- pref ori')
ylabel('D2- pref ori')
% xlim([0 .4])
% ylim([0 .4])
refline(1)
title('Observed pref ori')
saveas(fig, 'obsvOri.png')


%% extract fit orientation preference

%1 initial von M fit on real data/real orientations collected
vonMises_output = cell(1,nd);

for id = 1:nd
    vonMises_output_temp=nan(nKeep, 6);
    for i = 1:nKeep
        thisCell = data_ori_resp_keep{id}(i,:);
        [b_hat, k1_hat, R1_hat,u1_hat,sse,R_square] = miaovonmisesfit_ori(deg2rad(oris),thisCell);
        vonMises_output_temp(i,:)=[b_hat, k1_hat, R1_hat,u1_hat,sse,R_square];
    end
    vonMises_output{id}=vonMises_output_temp;

end
clear thisCell b_hat k1_hat R1_hat u1_hat sse R_square vonMises_output_temp %get rid of the last iteration of these looping variables
save(fullfile(fn_multi,'vonM_output.mat'),'vonMises_output')

%I can extract sharpness values from this structure, it is the second
%column


%2 fit with higher resolution
fit_oris = (0:1:179);

y_fits_keep = cell(1,nd); %an n cell X 180 degree matrix for each day representing the fit tuning curve
fit_pref_oris_keep = cell(1,nd); %a matrix for each day with the fit preferred ori

for id = 1:nd
    y_fits = zeros(nKeep, size(fit_oris,2));
    pref_oris = nan(nKeep, 1);
    for i = 1:nKeep
        in = vonMises_output{id}(i,:);
        b_tmp = in(1);
        k1_tmp = in(2);
        R1_tmp = in(3);
        u1_tmp = in(4);
        y_fits(i,:) = b_tmp+R1_tmp.*exp(k1_tmp.*(cos(2.*(deg2rad(fit_oris)-u1_tmp))-1));
        if size(fit_oris(find(y_fits(i,:)==max(y_fits(i,:)))),2)>1
            fit_pref_oris(i)=nan;
        else
            fit_pref_oris(i)= fit_oris(find(y_fits(i,:)==max(y_fits(i,:))));
        end
    end
    y_fits_keep{id}=y_fits;
    fit_pref_oris_keep{id}=fit_pref_oris;
clear b_tmp k1_tmp R1_tmp u1_tmp y_fits fit_pref_oris   
end

%% plots observed tuning curve with von M fit curve


figure;
movegui('center')
start = 1;

for iCell = 1:nKeep
        if start>tot
            figure; movegui('center')
            start = 1;
        end
        subplot(n,n2,start)
        for  id = 1:nd
            plot(oris, data_ori_resp_keep{id}(iCell,:))
            hold on
            plot(fit_oris,y_fits_keep{id}(iCell,:),':');
            hold on
        end
        
        if ismember(iCell,red_ind_keep)
            extra_title='HT+';
        else
            extra_title='HT-';
        end
            start= start+1;
        ylim([-0.1 inf])
        title([extra_title])
        
end


%% compare Von M fit pref ori across days
fig=figure; movegui('center') 
scatter(fit_pref_oris_keep{1}(green_ind_keep),fit_pref_oris_keep{2}(green_ind_keep),'k','jitter', 'on', 'jitterAmount', 5)
hold on
scatter(fit_pref_oris_keep{1}(red_ind_keep),fit_pref_oris_keep{2}(red_ind_keep),'MarkerEdgeColor',[.7 .05 .05],'jitter', 'on', 'jitterAmount', 5)
hold off
xlabel('D1- pref ori')
ylabel('D2- pref ori')
% xlim([0 .4])
% ylim([0 .4])
refline(1)
title('fit pref ori')
%saveas(fig, 'fitOri.png')

%% Von Mises bootstrap
reps = 1000;

shuff_pref_ori= cell(1, nd);
for id = 1:nd
    shuff_pref_ori{id} = nan(reps,nKeep);
    for r = 1:reps
        %make random sample - this will be redone every loop
        tCon = tCon_match{id}(:,1:nTrials);
        tDir = tDir_match{id}(:,1:nTrials);
        for i = 1:nKeep
            shuff_resp=nan(1,nOri);
            for iOri = 1:nOri
                con_inds=find(tCon==(pref_con_keep{id}(i)));
                ori_inds=find(tDir == oris(iOri));
                intersect_inds=intersect(con_inds,ori_inds);
                x=round(.8*length(intersect_inds));
                trial_inds = randsample(intersect_inds, x,1); %randomly select 80% of the trials at this contrast and ori
                temp_TCs=data_trial_keep{id}(resp_win,trial_inds,i); %extract the response window of these trials for this cell
               shuff_resp(1,iOri) = mean(squeeze(mean(temp_TCs)));
            end
            %now I can use the shuffled response by ori data to calculate
            %the von M fit for this cell on this rep
            
            [b_hat, k1_hat, R1_hat,u1_hat,sse,R_square] = miaovonmisesfit_ori(deg2rad(oris),shuff_resp);
            this_fits = b_hat+R1_hat.*exp(k1_hat.*(cos(2.*(deg2rad(fit_oris)-u1_hat))-1));
            this_pref = fit_oris(min(find(this_fits==max(this_fits))));
            shuff_pref_ori{id}(r,i)=this_pref;
        end
    end
end
save(fullfile(fnout,'bootstrap_fit_ori.mat'),'shuff_pref_ori');
clear shuff_resp this_pref this_fits b_hat k1_hat R1_hat u1_hat sse R_square x temp_TCs intersect_inds ori_inds con_inds iOri i r id 
%% compare von M bootstraps to original fit 
diff_pref = cell(1,nd);
well_fit = cell(1,nd);
for id = 1:nd
    diff_pref{id}=sort(abs(shuff_pref_ori{id} - fit_pref_oris_keep{id}),1);
    well_fit{id} = find(diff_pref{id}((reps*.9),:)<22.5);
end
%well_fit gives a list of the well-fit cell
figure; movegui('center') 
scatter(fit_pref_oris_keep{1}(well_fit{2}),fit_pref_oris_keep{2}(well_fit{2}),'k','jitter', 'on', 'jitterAmount', 5)
xlabel('D1- pref ori')
ylabel('D2- pref ori')
refline(1)
title('fit pref ori well-fit cells')

figure; movegui('center') 
subplot(1,2,1)
scatter(fit_pref_oris_keep{1}(green_ind_keep),pref_ori_keep{1}(green_ind_keep),'k','jitter', 'on', 'jitterAmount', 5)
hold on
scatter(fit_pref_oris_keep{1}(red_ind_keep),pref_ori_keep{1}(red_ind_keep),'MarkerEdgeColor',[.7 .05 .05],'jitter', 'on', 'jitterAmount', 5)
hold off
xlabel('fit pref ori')
ylabel('raw pref ori')
title('day 1')
refline(1)
subplot(1,2,2)
scatter(fit_pref_oris_keep{2}(green_ind_keep),pref_ori_keep{2}(green_ind_keep),'k','jitter', 'on', 'jitterAmount', 5)
hold on
scatter(fit_pref_oris_keep{2}(red_ind_keep),pref_ori_keep{2}(red_ind_keep),'MarkerEdgeColor',[.7 .05 .05],'jitter', 'on', 'jitterAmount', 5)
hold off
xlabel('fit pref ori')
ylabel('raw pref ori')
title('day 2')
refline(1)
%% contrast response curves



