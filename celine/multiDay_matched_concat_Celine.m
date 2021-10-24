clear all; clear global; close all
clc
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsAV; %directories
eval(ds)

%% depth 1
day_id1(2) = 113;
day_id1(1) = expt(day_id1(2)).multiday_matchdays;
nd = size(day_id1,2);
day_id1

mouse = expt(day_id1(1)).mouse;

fnout = fullfile(rc.celineAnalysis,mouse);
if expt(day_id1(2)).multiday_timesincedrug_hours>0
    dart_str = [expt(day_id1(2)).drug '_' num2str(expt(day_id1(2)).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end


fn_multi_in = fullfile(rc.celineAnalysis,mouse,['multiday_' dart_str]);

fn_multi = fullfile(rc.celineAnalysis,mouse,['multiday_concat_' dart_str]);
mkdir(fn_multi);
load(fullfile(fn_multi_in,'tc_keep.mat'));
load(fullfile(fn_multi_in,'resp_keep.mat'));
tc_trial_avrg_1 = tc_trial_avrg_keep;
red_keep_logical_1=red_keep_logical;
green_keep_logical_1=green_keep_logical;
data_resp_keep_1=data_resp_keep;
resp_max_concat_1=resp_max_keep;

clear tc_trial_avrg_keep red_keep_logical green_keep_logical resp_max_concat data_resp_keep
%% depth 2
day_id2(2) = 115;
day_id2(1) = expt(day_id2(2)).multiday_matchdays;
nd = size(day_id2,2);
day_id2

if expt(day_id1(2)).multiday_timesincedrug_hours>0
    dart_str = [expt(day_id2(2)).drug '_' num2str(expt(day_id2(2)).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end


fn_multi_in = fullfile(rc.celineAnalysis,mouse,['multiday_' dart_str]);


load(fullfile(fn_multi_in,'tc_keep.mat'));
load(fullfile(fn_multi_in,'resp_keep.mat'));
tc_trial_avrg_2 = tc_trial_avrg_keep;
red_keep_logical_2=red_keep_logical;
green_keep_logical_2=green_keep_logical;
data_resp_keep_2=data_resp_keep;
resp_max_concat_2=resp_max_keep;
clear tc_trial_avrg_keep  red_keep_logical green_keep_logical resp_max_concat data_resp_keep

%% combining
tc_concat = cell(1,2);
tc_concat{1}=[tc_trial_avrg_1{1}, tc_trial_avrg_2{1}];
tc_concat{2}=[tc_trial_avrg_1{2}, tc_trial_avrg_2{2}];

red_concat = [red_keep_logical_1, red_keep_logical_2];
green_concat = [green_keep_logical_1, green_keep_logical_2];

resp_concat = cell(1,2);
resp_max_concat = cell(1,2);
%resp_concat{1}=[data_resp_keep_1{1}, data_resp_keep_2{1}]
resp_max_concat{1}=[resp_max_concat_1{1}; resp_max_concat_2{1}]
resp_max_concat{2}=[resp_max_concat_1{2}; resp_max_concat_2{2}]



%% print counts to console
fprintf(['Green cells ' num2str(sum(green_concat)) '; Red cell ' num2str(sum(red_concat))])


%% prepare to plot the timecourses 

tc_green_avrg_match = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_match = cell(1,nd); %same for red
tc_green_se_match = cell(1,nd); %this will be the se across all green cells
tc_red_se_match = cell(1,nd); %same for red

for id = 1:nd
    tc_green_avrg_match{id}=mean(tc_concat{id}(:,find(green_concat)),2);
    green_std=std(tc_concat{id}(:,find(green_concat)),[],2);
    tc_green_se_match{id}=green_std/sqrt(sum(green_concat));
    
    tc_red_avrg_match{id}=mean(tc_concat{id}(:,find(red_concat)),2);
    red_std=std(tc_concat{id}(:,find(red_concat)),[],2);
    tc_red_se_match{id}=red_std/sqrt(sum(red_concat));
    
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
%% make a plot of individual timecourses 
figure
subplot(2,2,1)
plot(tc_concat{1}(:,find(green_concat)),'k')
ylim([-.1 .5]);
title('day 1')

subplot(2,2,2)
plot(tc_concat{1}(:,find(red_concat)),'color',[.7 .05 .05])
ylim([-.1 .5]);
title('day 1')

subplot(2,2,3)
plot(tc_concat{2}(:,find(green_concat)),'k')
ylim([-.1 .5]);
title('day 2')

subplot(2,2,4)
plot(tc_concat{2}(:,find(red_concat)),'color',[.7 .05 .05])
ylim([-.1 .5]);
title('day 2')
print(fullfile(fn_multi,['indiv_timecourses']),'-dpdf');
%% makes a scatterplot of max df/f for day 1 vs day 2, and each subplot is one day
%this is for all cells I'm keeping, red and green




figure; movegui('center') 
subplot(1,2,1)
scatter(resp_max_concat{1}(find(green_concat)),resp_max_concat{2}(find(green_concat)),'k')
% hold on
% scatter(resp_max_concat{1}(find(red_concat)),resp_max_concat{2}(find(red_concat)),'MarkerEdgeColor',[.7 .05 .05])
% hold off
xlabel('D1- max dF/F')
ylabel('D2- max dF/F')
xlim([0 .5])
ylim([0 .5])
refline(1)
title('Max df/f for responsive HT- ')


subplot(1,2,2)
scatter(resp_max_concat{1}(find(red_concat)),resp_max_concat{2}(find(red_concat)),'MarkerEdgeColor',[.7 .05 .05])
% hold on
% scatter(resp_max_match{1}(red_ind_match_list),resp_max_match{2}(red_ind_match_list),'MarkerEdgeColor',[.7 .05 .05])
% hold off
xlabel('D1- max dF/F')
ylabel('D2- max dF/F')
xlim([0 .5])
ylim([0 .5])
refline(1)
title('Max df/f HT+')

print(fullfile(fn_multi,'maxResp_crossDay.pdf'),'-dpdf','-bestfit')

% extract the max df/f values for analysis
%% looking at change in dfof
dfof_max_diff = (resp_max_concat{1}-resp_max_concat{2})./resp_max_concat{1};

figure
x = [mean(dfof_max_diff(find(green_concat))), mean(dfof_max_diff(find(red_concat)))];
y = [(std(dfof_max_diff(find(green_concat))))/sqrt(length(find(green_concat))), (std(dfof_max_diff(find(red_concat))))/sqrt(length(find(red_concat)))];
%y = [std(dfof_max_diff(find(green_concat))), std(dfof_max_diff(find(red_concat)))]
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

fprintf(['Green cells ' num2str(mean(dfof_max_diff(find(green_concat)))) 's.e.m.' num2str(std(dfof_max_diff(find(green_concat)))/sqrt(sum(green_concat))  )])
fprintf(['Red cells ' num2str(mean(dfof_max_diff(find(red_concat)))) 's.e.m.' num2str(std(dfof_max_diff(find(red_concat)))/sqrt(sum(red_concat))  )])