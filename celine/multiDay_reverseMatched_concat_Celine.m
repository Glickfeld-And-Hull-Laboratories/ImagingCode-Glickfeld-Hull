
clear all; clear global; close all
clc
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc =  behavConstsDART; %directories
eval(ds);

sess_list = [131 133];%enter all the sessions you want to concatenate
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
fnout= fullfile(rc.achAnalysis,strcat('concat', sess_title),d);
mkdir(fnout);
cd(fnout)
clear d sess_title
%% concatenating data
tc_trial_avrg_keep_concat=cell(1,nd);
tc_trial_avrg_keep_allCon_concat=cell(1,nd);
resp_keep_concat=cell(1,nd);
resp_max_keep_concat=cell(1,nd);
oris_concat=[];
cons_concat=[];
green_concat=[];
red_concat=[];
dfof_max_diff_concat=[];
nKeep_concat=[];
LMI_concat = cell(1,nd);
data_resp_concat = cell(1,nd);
locResp_concat = cell(1,nd);
statResp_concat = cell(1,nd);
pref_con_concat = cell(1,nd);
%add locomotion and and stationary TCs?
for iSess = 1:nSess
    day_id = sess_list(iSess)
    mouse = expt(day_id).mouse;

    if expt(day_id).multiday_timesincedrug_hours>0
        dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
    else
        dart_str = 'control';
    end
    fn_multi = fullfile(rc.achAnalysis,mouse,['multiday_' dart_str]);

    load(fullfile(fn_multi,'tc_keep.mat'))
    load(fullfile(fn_multi,'resp_keep.mat'))
    load(fullfile(fn_multi,'input.mat'))
    load(fullfile(fn_multi,'locomotion.mat'))

    nKeep = size(tc_trial_avrg_keep{post},2);


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
    oris = unique(tOri_match{post});
    cons = unique(tCon_match{post});

    nOn = input(1).nScansOn;
    nOff = input(1).nScansOff;
    %start conatenating
    oris_concat = [oris_concat;oris]; %I will organize thisas rows so I can subsequently make sure everything matches
    cons_concat = [cons_concat;cons];
    red_concat = [red_concat, red_keep_logical];
    green_concat = [green_concat, green_keep_logical];
    nKeep_concat = [nKeep_concat,nKeep];
    
    for id = 1:nd
        tc_trial_avrg_keep_concat{id} =cat(2,tc_trial_avrg_keep_concat{id},tc_trial_avrg_keep{id});
        tc_trial_avrg_keep_allCon_concat{id} =cat(2,tc_trial_avrg_keep_allCon_concat{id},tc_trial_avrg_keep_allCon{id});
        resp_keep_concat{id}=cat(1,resp_keep_concat{id},resp_keep{id});
        resp_max_keep_concat{id}=cat(1,resp_max_keep_concat{id},resp_max_keep{id});
        LMI_concat{id}=cat(1,LMI_concat{id},LMI{id});
        data_resp_concat{id}=cat(1,data_resp_concat{id},data_resp_keep{id});
        locResp_concat{id}=cat(1,locResp_concat{id},locResp{id});
        statResp_concat{id}=cat(1,statResp_concat{id},statResp{id});
        pref_con_keep{id} = pref_con_keep{id}'; %because I output it in the wrong dimensions 
        pref_con_concat{id}=cat(1,pref_con_concat{id},pref_con_keep{id});
    end
    dfof_max_diff_concat=cat(1,dfof_max_diff_concat,dfof_max_diff);
end

clear mouse day_id nKeep iSess fn_multi cons oris
clear explanation1 resp_keep tc_trial_avrg_keep_allCon pref_responses_allCon sig_diff pref_con_keep pref_ori_keep tOri_match tCon_match data_trial_keep nTrials tc_trial_avrg_keep green_keep_logical red_keep_logical green_ind_keep red_ind_keep
clear LMI RIx locCounts locResp locTCs statResp statTCs wheel_tc
clear data_con_resp_keep data_ori_resp_keep data_rep_keep dfof_max_diff dfof_max_diff_raw explanation2 resp_max_keep data_resp_keep
red_ind_concat = find(red_concat);
green_ind_concat = find(green_concat);
%% check that contrasts and oris are matched
cons=unique(cons_concat);
nCon = length(cons)
oris=unique(oris_concat);
nOri=length(oris);
nKeep_total = sum(nKeep_concat);


%% make figure with se shaded, one figure per contrast
%ADD stim on
tc_green_avrg_keep = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_keep = cell(1,nd); %same for red
tc_green_se_keep = cell(1,nd); %this will be the se across all green cells
tc_red_se_match = cell(1,nd); %same for red

for id = 1:nd
    for iCon=1:nCon
    tc_green_avrg_keep{id}(:,iCon)=nanmean(tc_trial_avrg_keep_concat{id}(:,green_ind_concat,iCon),2);
    green_std=std(tc_trial_avrg_keep_concat{id}(:,green_ind_concat,iCon),[],2);
    tc_green_se_keep{id}(:,iCon)=green_std/sqrt(length(green_ind_concat));
    
    tc_red_avrg_keep{id}(:,iCon)=nanmean(tc_trial_avrg_keep_concat{id}(:,red_ind_concat,iCon),2);
    red_std=std(tc_trial_avrg_keep_concat{id}(:,red_ind_concat,iCon),[],2);
    tc_red_se_match{id}(:,iCon)=red_std/sqrt(length(red_ind_concat));
    
    clear green_std red_std
    end
end


%create a time axis in seconds
t=1:(size(tc_green_avrg_keep{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
figure
subplot(1,2,1) %for the first day

shadedErrorBar(t,tc_red_avrg_keep{pre}(:,iCon),tc_red_se_match{pre}(:,iCon),'r');
ylim([-.02 .3]);
hold on
shadedErrorBar(t,tc_green_avrg_keep{pre}(:,iCon),tc_green_se_keep{pre}(:,iCon));
title(['Pre-DART contrast = ' num2str(cons(iCon))])
txt1 = ['HT- ' num2str(length(green_ind_concat))];
text(-1.5,0.25,txt1);
txt2 = ['HT+ ' num2str(length(red_ind_concat))];
text(-1.5,0.23,txt2,'Color','r');
ylabel('dF/F') 
xlabel('s') 

axis square


subplot(1,2,2) %for the second day
shadedErrorBar(t,tc_red_avrg_keep{post}(:,iCon),tc_red_se_match{post}(:,iCon),'r');
ylim([-.02 .3]);
hold on
shadedErrorBar(t,tc_green_avrg_keep{post}(:,iCon),tc_green_se_keep{post}(:,iCon));
ylabel('dF/F') 
xlabel('s') 
title(['Post-DART contrast = ' num2str(cons(iCon))])
axis square
x0=5;
y0=5;
width=6;
height=3
set(gcf,'units','inches','position',[x0,y0,width,height])

print(fullfile(fnout,[num2str(cons(iCon)) '_timecourses.pdf']),'-dpdf');
end 
clear txt1 txt2
%% make a plot of individual timecourses 

setYmin = -.2; %indicate y axes you want
setYmax = 0.85;


for iCon = 1:nCon

figure
subplot(2,2,1)
plot(t, tc_trial_avrg_keep_concat{pre}(:,green_ind_concat,iCon),'color',[0 0 0 .2])
hold on
plot(t, tc_green_avrg_keep{pre}(:,iCon),'color',[0 0 0],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 1')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
hold off


subplot(2,2,3)
plot(t, tc_trial_avrg_keep_concat{pre}(:,red_ind_concat,iCon),'color',[.7 .05 .05 .2])
hold on
plot(t, tc_red_avrg_keep{pre}(:,iCon),'color',[.7 .05 .05],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 1')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
hold off


subplot(2,2,2)
plot(t, tc_trial_avrg_keep_concat{post}(:,green_ind_concat,iCon),'color',[0 0 0 .2])
hold on
plot(t, tc_green_avrg_keep{post}(:,iCon),'color',[0 0 0],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 2')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
hold off


subplot(2,2,4)
plot(t, tc_trial_avrg_keep_concat{post}(:,red_ind_concat,iCon),'color',[.7 .05 .05 .2])
hold on
plot(t, tc_red_avrg_keep{post}(:,iCon),'color',[.7 .05 .05],'LineWidth',2)
ylim([setYmin setYmax]);
title('day 2')
xlim([-2 4])
ylabel('dF/F') 
xlabel('s') 
sgtitle(['contrast = '  num2str(cons(iCon))])
hold off

x0=5;
y0=5;
width=6;
height=6;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,[num2str(cons(iCon)) '_indiv_timecourses.pdf']),'-dpdf');
end
%% time to peak

tHalfMax = cell(1,nd);

for id = 1:nd
    tHalfMaxTemp=nan(nKeep_total,1);

        for iCell=1:nKeep_total
            %pull the data for a given cell at a given contrast (pref ori)
            tempData=tc_trial_avrg_keep_allCon_concat{id}(stimStart:stimStart+nOn,iCell);
            if resp_keep_concat{id}(iCell)
                smoothData=smoothdata(tempData,'movmean',5) ;
                halfMax = max(smoothData(3:length(smoothData)))/2;
                tHalfMaxCell =double(min(find(smoothData>halfMax)))/double(frame_rate);
                if rem(iCell, 10) == 0 
%                 figure;plot(tempData)
%                 hold on
%                 plot(smoothData)
%                 hold on
%                 vline(min(find(smoothData>halfMax)))
                end
                if length(tHalfMaxCell)>0
                    tHalfMaxTemp(iCell)=tHalfMaxCell;
                end
            end
 
       end
    tHalfMax{id}=tHalfMaxTemp;
end
clear tHalfMaxCell tHalfMaxTemp tempData smoothData halfMax
%% scatter for tMax




figure; movegui('center') 
subplot(1,2,1)
scatter((tHalfMax{pre}(green_ind_concat)),(tHalfMax{post}(green_ind_concat)),'k','jitter', 'on', 'jitterAmount', 0.3)
hold on
ylabel('post-DART half-max(s)')
xlabel('pre-DART half-max(s)')
ylim([0 2.5])
xlim([0 2.5])
refline(1)
title('HT- ')
axis square
hold on
plot(nanmean(tHalfMax{pre}(green_ind_concat)),nanmean(tHalfMax{post}(green_ind_concat)),'b*')
hs = findobj(gca, 'Type','scatter')
NrPoints = numel(hs.XData);
txt1 = ['n = ', num2str(NrPoints)];
text(1,2.2,txt1);
hold off



subplot(1,2,2)
scatter((tHalfMax{pre}(red_ind_concat)),(tHalfMax{post}(red_ind_concat)),'MarkerEdgeColor',[.7 .05 .05],'jitter', 'on', 'jitterAmount', 0.3)
hold on
ylabel('post-DART half-max(s)')
xlabel('pre-DART half-max(s)')
ylim([0 2.5])
xlim([0 2.5])
refline(1)
title('HT+')
axis square
hold on
plot(nanmean(tHalfMax{pre}(red_ind_concat)),nanmean(tHalfMax{post}(red_ind_concat)),'b*')
hs = findobj(gca, 'Type','scatter')
NrPoints = numel(hs.XData);
txt1 = ['n = ' num2str(NrPoints)];
text(1,2.2,txt1,'Color','r');
hold off

clear txt1 NrPoints

sgtitle('time to half max (s), averaged over contrast')
print(fullfile(fnout, ['tHalfMax_crossDay.pdf']),'-dpdf','-bestfit')
%% scatterplot of max df/f for day 1 vs day 2, and each subplot is one cell type


for iCon = 1:nCon
figure; movegui('center') 
subplot(1,2,1)
swarmchart((resp_max_keep_concat{pre}(green_ind_concat,iCon)),(resp_max_keep_concat{post}(green_ind_concat,iCon)),'k')
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 .5])
xlim([-.1 .5])
refline(1)
title('HT- ')
axis square


subplot(1,2,2)
swarmchart((resp_max_keep_concat{pre}(red_ind_concat,iCon)),(resp_max_keep_concat{post}(red_ind_concat,iCon)),'MarkerEdgeColor',[.7 .05 .05])
hold on
ylabel('post-DART dF/F')
xlabel('pre-DART  dF/F')
ylim([-.1 .5])
xlim([-.1 .5])

refline(1)
title('HT+')
axis square

sgtitle(num2str(cons(iCon)))
print(fullfile(fnout,[num2str(cons(iCon)) 'maxResp_crossDay.pdf']),'-dpdf','-bestfit')
end
%% locomotion modulation index

for iCon = 1:nCon

figure; movegui('center') 
subplot(1,2,1)
scatter((LMI_concat{pre}(green_ind_concat,iCon)),(LMI_concat{post}(green_ind_concat,iCon)),'k','jitter', 'on', 'jitterAmount', 0.3)
ylabel('post-DART LMI')
xlabel('pre-DART  LMI')
ylim([-1 1])
xlim([-1 1])
hline(0)
vline(0)
refline(1)
title('HT- ')
axis square


subplot(1,2,2)
scatter((LMI_concat{pre}(red_ind_concat,iCon)),(LMI_concat{post}(red_ind_concat,iCon)),'MarkerEdgeColor',[.7 .05 .05],'jitter', 'on', 'jitterAmount', 0.3)
ylabel('post-DART LMI')
xlabel('pre-DART  LMI')
ylim([-1 1])
xlim([-1 1])
hline(0)
vline(0)
refline(1)
title('HT+')
axis square

sgtitle(num2str(cons(iCon)))
print(fullfile(fnout,[num2str(cons(iCon)) '_LMI.pdf']),'-dpdf');
end
%% scatter of stat vs. loc 

for iCon = 1:nCon

figure
subplot(2,2,1)
swarmchart(statResp_concat{pre}(green_ind_concat,iCon),locResp_concat{pre}(green_ind_concat,iCon),'k')
title('Pre-DART')
ylim([-.1 .5]);
xlim([-.1 .5]);
ylabel('Running') 
xlabel('Stationary')
axis square
refline(1)


subplot(2,2,3)
swarmchart(statResp_concat{pre}(red_ind_concat,iCon),locResp_concat{pre}(red_ind_concat,iCon),'MarkerEdgeColor',[.7 .05 .05])
title('Pre-DART')
ylim([-.1 .5]);
xlim([-.1 .5]);
ylabel('Running') 
xlabel('Stationary')
axis square
refline(1)


subplot(2,2,2)
swarmchart(statResp_concat{post}(green_ind_concat,iCon),locResp_concat{post}(green_ind_concat,iCon),'k')
title('Post-DART')
ylim([-.1 .5]);
xlim([-.1 .5]);
ylabel('Running') 
xlabel('Stationary')  
axis square
refline(1)



subplot(2,2,4)
swarmchart(statResp_concat{post}(red_ind_concat,iCon),locResp_concat{post}(red_ind_concat,iCon),'MarkerEdgeColor',[.7 .05 .05])
title('Post-DART')
ylim([-.1 .5]);
xlim([-.1 .5]);
ylabel('Running') 
xlabel('Stationary') 
axis square
refline(1)


sgtitle(['contrast = '  num2str(cons(iCon))])

print(fullfile(fnout,[num2str(cons(iCon)) '_locVSstat.pdf']),'-dpdf');
end


%% plot all responses 
if nKeep_total<36
    [n n2] = subplotn(nKeep_total);
    tot = n.*n2;
else
    n = 5;
    n2 = 5;
    tot = 25;
end

figure;
movegui('center')
start = 1;

for iCell = 1:nKeep_total
     
        if start>tot

            figure; movegui('center')
           
            start = 1;
            
        end
        subplot(n,n2,start)
            for iCon = 1:nCon
                scatter(data_resp_concat{pre}(iCell,:,iCon,1),data_resp_concat{post}(iCell,:,iCon,1))
                %the contrasts are plotted indiviudally so that the points
                %corresponding to different contrasts are different colors,
                hold on
                lsline
            end
                if ismember(iCell,red_ind_concat)
                    title('\color{red}HT+')
                else
                    title('HT-')
                end
            %max_list and upperLim are to find a single axis upper limit
            %that will work for botht eh x and y axes
            max_list = [max(max(data_resp_concat{pre}(iCell,:,:,1))),max(max(data_resp_concat{post}(iCell,:,:,1)))];
            upperLim =max(max_list)+0.1;
            ylim([-.05 upperLim]);
            xlim([-.05 upperLim]);
            axis square
            hold off
            start=start+1;
                
end
%%
linCell=zeros(1,nKeep_total);%will provide a boolean index of which cells have a linear relationship
linCellProps = nan(4,nKeep_total); %this will contain the slopes and intercepts of cells that are deemed linear
data_pre_all=nan(nOri*nCon, nKeep_total);
count=1;
figure;
for iCell = 1:nKeep_total
    if ismember(iCell,red_ind_concat)
 data_pre=[];
 
 data_post=[];
  for iCon = 1:nCon
                data_pre=[data_pre,data_resp_concat{pre}(iCell,:,iCon,1)]; %collects the values for day 1 over all contrasts
                data_post=[data_post,data_resp_concat{post}(iCell,:,iCon,1)]; %collects the values for day 2 over all contrasts
  end
  data_pre_all(:,iCell)=data_pre;
   [R,p]=corrcoef(data_pre,data_post); %gets the correlation data for
   %the correlation between days, pooling across contrasts
    if R > 0.5
        linCell(iCell)=1;
        linFit = polyfit(data_pre, data_post, 1);
        linCellProps(1,iCell)=linFit(1); %slope
        linCellProps(2,iCell)=linFit(2); %intercept
        linCellProps(3,iCell)=R(2);
        linCellProps(4,iCell)=p(2);
        f=polyval(linFit,data_pre);
        plot(data_pre,f,'color',[0 0 0 .2],'LineWidth',.5)
        hold on
        scatter(data_pre,data_post,'ko')
        alpha(.3)
        count=count+1;

    end
    hold on 
    
    end

end
% collecting info for the linear cells
        
meanFit= [mean(linCellProps(1,find(linCell))),  mean(linCellProps(2,find(linCell)))];
xRange=[min(min(data_pre_all)), max(max(data_pre_all))]
f2=polyval(meanFit,xRange);
plot(xRange,f2,'b-','LineWidth',2)
xlim([-.1 .3])
txt1 = ['y=', num2str(round(meanFit(1),2)),'x',num2str(round(meanFit(2),4))];
text(.15,0.01,txt1);
hold off

myRef = refline(1)
myRef.LineStyle = ':'

title(['\color{red}HT+ cells, all stimulus conditions, n= ' num2str(length(intersect(red_ind_concat,find(linCell))))])
xlabel('pre-DART')
ylabel('post-DART')
myRef = refline(1)
myRef.LineStyle = ':'

axis square
print(fullfile(fnout,'HT+_allResp.pdf'),'-dpdf','-bestfit')


%% same as above for HT- cells


linCell=zeros(1,nKeep_total);%will provide a boolean index of which cells have a linear relationship
linCellProps = nan(4,nKeep_total); %this will contain the slopes and intercepts of cells that are deemed linear
data_pre_all=nan(nOri*nCon, nKeep_total);
count=1;
figure;
for iCell = 1:nKeep_total
    if ~ismember(iCell,red_ind_concat)
 data_pre=[];
 
 data_post=[];
  for iCon = 1:nCon
                data_pre=[data_pre,data_resp_concat{pre}(iCell,:,iCon,1)]; %collects the values for day 1 over all contrasts
                data_post=[data_post,data_resp_concat{post}(iCell,:,iCon,1)]; %collects the values for day 2 over all contrasts
  end
  data_pre_all(:,iCell)=data_pre;
   [R,p]=corrcoef(data_pre,data_post); %gets the correlation data for
   %the correlation between days, pooling across contrasts
    if R > 0.5
        linCell(iCell)=1;
        linFit = polyfit(data_pre, data_post, 1);
        linCellProps(1,iCell)=linFit(1); %slope
        linCellProps(2,iCell)=linFit(2); %intercept
        linCellProps(3,iCell)=R(2);
        linCellProps(4,iCell)=p(2);
        f=polyval(linFit,data_pre);
        plot(data_pre,f,'color',[0 0 0 .2],'LineWidth',.5)
        hold on
        scatter(data_pre,data_post,'ko')
        alpha(.3)
        %FIND A WAY TO CHANGE ALPHA FOR SCATTERS
        
        count=count+1;

    end
    hold on    
    end

end
% collecting info for the linear cells
        
meanFit= [mean(linCellProps(1,find(linCell))),  mean(linCellProps(2,find(linCell)))];
xRange=[min(min(data_pre_all)), max(max(data_pre_all))]
f2=polyval(meanFit,xRange);
plot(xRange,f2,'b-','LineWidth',2)
txt1 = ['y=', num2str(round(meanFit(1),2)),'x',num2str(round(meanFit(2),4))];
text(.2,-0.01,txt1);
hold off
axis square
myRef = refline(1)
myRef.LineStyle = ':'

title(['HT- cells, all stimulus conditions, n= ' num2str(length(intersect(green_ind_concat,find(linCell))))])
xlabel('pre-DART')
ylabel('post-DART')
myRef = refline(1)
myRef.LineStyle = ':'

print(fullfile(fnout,'HT-_allResp.pdf'),'-dpdf','-bestfit')
%% comparing contrast preference
%need more contrasts

conTable = table([mean(pref_con_concat{pre}(green_ind_concat));mean(pref_con_concat{pre}(red_ind_concat))],[mean(pref_con_concat{post}(green_ind_concat));mean(pref_con_concat{post}(red_ind_concat))],'VariableNames',{'mean pref con pre' 'mean pref con post'}, 'RowNames',{'HT- cells cells'  'HT+ cells'})
writetable(conTable,fullfile(fnout,'conPref.csv'),'WriteRowNames',true)

boxplot(pref_con_concat{pre},red_concat)
boxplot(pref_con_concat{post},red_concat)