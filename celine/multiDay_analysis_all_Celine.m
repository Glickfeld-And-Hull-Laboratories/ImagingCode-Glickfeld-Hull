%%
cd(fn_multi)
PC_ind_all = cell(1,2);

PC_ind_all{1} = find(red_ind_all{1} ==0);
PC_ind_all{2} = find(red_ind_all{2} ==0);
%% first plot timecourses

PCinds = PC_ind_all{1};
redinds= find(red_ind_all{1});

data_tc_trial{1} = reshape(data_dfof_all{1}, [nOn+nOff,nTrials,nCells_d1]);
data_tc_trial{2} = reshape(data_dfof_all{2}, [nOn+nOff,nTrials,nCells_d2]);
%looking at day 1
tc_cell_avrg1 = mean(data_tc_trial{1},3);%average pver cells, one row per trial
tc_trial_avrg1 = squeeze(mean(data_tc_trial{1},2));%average over trials, one row per cell
%tc_trial_std = squeeze(std(data_tc_trial{1},2));%std over trials for each cell
%std over cells, after averaging over trials
tc_cell_trial_avrg1 = mean(tc_cell_avrg1,2);%average over trials and cells
tc_trial_avrg_PC1 = tc_trial_avrg1(:,PCinds);
tc_trial_avrg_IN1 = tc_trial_avrg1(:,redinds);


figure;
plot(tc_trial_avrg_PC1, 'LineWidth',.005,'color',[.25 .25 .25]);
hold on;
plot(tc_trial_avrg_IN1, 'LineWidth',.005,'color',[.7 0 0]);
hold on;
plot(tc_cell_trial_avrg1, 'LineWidth',2, 'color','k');
hold on;
vline(60,'g')
title('Timecourses after np subtraction');
dim = [.2 .5 .3 .3];
str = strcat(num2str(nCells),' cells averaged over', num2str(nTrials), ' trials') ;
title('Day 1 timecourses');
hold off

% looking at day 2
PCinds = PC_ind_all{2};
redinds= find(red_ind_all{2});
tc_cell_avrg2 = mean(data_tc_trial{2},3);%average pver cells, one row per trial
tc_trial_avrg2 = squeeze(mean(data_tc_trial{2},2));%average over trials, one row per cell
tc_cell_trial_avrg2 = mean(tc_cell_avrg2,2);%average over trials and cells
tc_trial_avrg_PC2 = tc_trial_avrg2(:,PCinds);
tc_trial_avrg_IN2 = tc_trial_avrg2(:,redinds);

figure;
plot(tc_trial_avrg_PC2, 'LineWidth',.005,'color',[.25 .25 .25]);
hold on;
plot(tc_trial_avrg_IN2, 'LineWidth',.005,'color',[.7 0 0]);
hold on;
plot(tc_cell_trial_avrg2, 'LineWidth',2, 'color','k');
hold on;
vline(60,'g')
title('Timecourses after np subtraction');
dim = [.2 .5 .3 .3];
str = strcat(num2str(nCells),' cells averaged over', num2str(nTrials), ' trials') ;
title('Day 2 timecourses');
hold off

%% 1 evoked activity

%was there a significant difference in visually-evoked responsiveness


%using the previously determined responsiveness, resp_ind_all
%isolate the responsiveness of PV vs pyr cells
PC_h_all = cell(1,2);
red_h_all = cell(1,2);
for i=1:2
  PCinds = PC_ind_all{i};
  redinds= find(red_ind_all{i});
  PC_h_all{i} = h_all{i}(:,:,PCinds);
  red_h_all{i} = h_all{i}(:,:,redinds);
end
clear PCinds redinds

PC_resp_all = cell(1,2);
red_resp_all = cell(1,2);

for i=1:2
  PC_resp_all{i} = squeeze(logical(sum(sum(PC_h_all{i},1),2)));
  red_resp_all{i} = squeeze(logical(sum(sum(red_h_all{i},1),2)));
end

%this prodiveds logical values indicating whether a cell was responsive to
%at least one stim, which I can then use to determing proportions of
%responsive cells

% now I want to compare the percentage of cells that were significantly
%responsive on day1 vs day 2

% comparing dfof in the response window - I want to determine whether the
% magnitude of the response was different
% not sure that comparing dfof values directly makes sense - maybe I should
% look at the magnitude of change from the baseline condition
%I will look at the average response (averaging over trials of each
%stimulus condition)

%for putative pyr cells and for PV cells, find the mean data
PC_avgResp_all = cell(1,2);
red_avgResp_all = cell(1,2);


for i=1:2
  PCinds = PC_ind_all{i};
  redinds= find(red_ind_all{i});
  PC_avgResp_all{i} = squeeze(mean(mean(resp_avg_all{i}(:,:,PCinds,1))));
  red_avgResp_all{i} = squeeze(mean(mean(resp_avg_all{i}(:,:,redinds,1))));
end
clear PCinds redinds




%% 2 spontaneous activity
%was there a significant different in spontaneous (not visually evoked)


PC_all_fractActive = cell(1,2);
red_all_fractActive = cell(1,2);
for i=1:2
   PCinds = PC_ind_all{i};
   redinds= find(red_ind_all{i});
   PC_all_fractActive{i} = fractTimeActive_all{i}(PCinds);
   mean(PC_all_fractActive{i});
   red_all_fractActive{i} = fractTimeActive_all{i}(redinds);
   mean(red_all_fractActive{i});
end



%comparing spontaneous activity between pyr and PV cells

%% 3 tuning
% was there a significant difference in how many orientations a cells
% responded to?

%comparing the number of stimuli that each cell was responsive to
PC_nResp_all = cell(1,2);
red_nResp_all = cell(1,2);
for i=1:2

  PC_nResp_all{i} = squeeze(sum(sum(PC_h_all{i},1),2));
  red_nResp_all{i} = squeeze(sum(sum(red_h_all{i},1),2));
end

% the number of stimulus conditions a cell responds to could be influenced
% by the number of orientations or by the number of contrasts, so I want to
% differentiate those

%comparing the number of stimuli that each cell was responsive to
%how many different orientations did each cell respond to, regardless of
%contrast?
PC_nOris_all = cell(1,2);
red_nOris_all = cell(1,2);
for i=1:2
  PC_nOris_all{i} = sum(logical(squeeze(sum(PC_h_all{i},2))),1);
  red_nOris_all{i} = sum(logical(squeeze(sum(red_h_all{i},2))),1);
end

% find preferred orientation



% insert von mises sharpess, preferred ori comparisons



%% 4 contrast response function
% looking only at cells that were responsive - at their preferred ori (the
% ori that gave the highest response) what did the relationship of dfof to
% contrast look like?
fprintf(fid,'\n\n\nLooking at contrast sensitivity');

PC_pref_resp_all = cell(1,2);
red_pref_resp_all = cell(1,2);
for i=1:2
  PC_pref_resp_all{i} = resp_avg_all_pref{i}(:,find(red_ind_all==0));
  %PC_pref_resp_all{i} = PC_pref_resp_all{i}(:,all(~isnan(PC_pref_resp_all{i})));   % remove NaN columns
  red_pref_resp_all{i} = resp_avg_all_pref{i}(:,find(red_ind_all==1));
  %red_pref_resp_all{i} = red_pref_resp_all{i}(:,all(~isnan(red_pref_resp_all{i})));   % remove NaN columns
end

%plotting to make sure that I'm looking at the same thing as in the
%multi-panel figure
figure;
errorbar(cons, nanmean(PC_pref_resp_all{1},2), nanstd(PC_pref_resp_all{1},[],2)./sqrt(length(find(red_ind_all==0))),'-o')
hold on
title('putative pyr cells')
ylim([-0.05 0.3])
errorbar(cons, nanmean(PC_pref_resp_all{2},2), nanstd(PC_pref_resp_all{2},[],2)./sqrt(length(find(red_ind_all==0))),'-o')
hold off

figure;
errorbar(cons, nanmean(red_pref_resp_all{1},2), nanstd(red_pref_resp_all{1},[],2)./sqrt(length(find(red_ind_all==0))),'-o')
hold on
title('PV/SOM cells')
ylim([-0.05 0.3])
errorbar(cons, nanmean(red_pref_resp_all{2},2), nanstd(red_pref_resp_all{2},[],2)./sqrt(length(find(red_ind_all==0))),'-o')
hold off

%convert to slopes
PC_all_slopes =cell(1,2);
red_all_slopes =cell(1,2);
for i = 1:2
    for iCell = 1:size(PC_pref_resp_all{i},2)
        PC_all_slopes{i}(iCell)=(PC_pref_resp_all{i}(5,iCell)-PC_pref_resp_all{i}(1,iCell))/5;
    end
    for iCell = 1:size(red_pref_resp_all{i},2)
        red_all_slopes{i}(iCell)=(red_pref_resp_all{i}(5,iCell)-red_pref_resp_all{i}(1,iCell))/5;
    end
end

%% 
PC1_mat = [PC_resp_all{1}, PC_nResp_all{1},PC_avgResp_all{1},PC_all_fractActive{1}',PC_nOris_all{1}'];
cell_type=zeros(size(PC1_mat,1),1);
PC1_mat = [cell_type, PC1_mat];

PC2_mat = [PC_resp_all{2}, PC_nResp_all{2},PC_avgResp_all{2},PC_all_fractActive{2}',PC_nOris_all{2}'];
cell_type=zeros(size(PC2_mat,1),1);
PC2_mat = [cell_type, PC2_mat];

IN1_mat = [red_resp_all{1}, red_nResp_all{1},red_avgResp_all{1},red_all_fractActive{1}',red_nOris_all{1}'];
cell_type=ones(size(IN1_mat,1),1);
IN1_mat = [cell_type, IN1_mat];

IN2_mat = [red_resp_all{2}, red_nResp_all{2},red_avgResp_all{2},red_all_fractActive{2}',red_nOris_all{2}'];
cell_type=ones(size(IN2_mat,1),1);
IN2_mat = [cell_type, IN2_mat];

day1 = [PC1_mat;IN1_mat];
day=zeros(size(day1,1),1);
day1 = [day, day1];

day2 = [PC2_mat;IN2_mat];
day=zeros(size(day2,1),1);
day2 = [day, day2];

finalData = [day1;day2];

clear day1 day2 PC1_mat IN1_mat PC2_mat IN2_mat cell_type day

colNames = {'Day_0_1','CellType_0_1','Responsive_Y_N', 'nResp_conds','avrg_resp_mag','fracActive','nOris_resp'};
output = array2table(finalData, 'VariableNames',colNames);