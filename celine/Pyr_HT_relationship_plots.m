%pyr som relationship plots
clc
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc =  behavConstsDART; %directories
eval(ds);
%131 133 138 142 163 171
sess_list = [131 133 138 142 163 171];%enter all the sessions you want to concatenate
nSess=length(sess_list);

%% loop through sessions

linPropsAll = cell(1,nSess);
respByCondPropsAll=cell(1,nSess);
for iSess = 1:nSess
    day_id = sess_list(iSess)



    ds = 'DART_V1_contrast_ori_Celine'; %dataset info
    dataStructLabels = {'contrastxori'};
    
    rc = behavConstsDART; %directories
    eval(ds);
        
    
    pre_day = expt(day_id).multiday_matchdays;
    
    nd=2; %hardcoding the number of days for now
    
    % INDICATE THE PRE VS. POST DAYS DEPENDING ON THE ORDER OF MATCHING
    pre=2;
    post=1;
    
    mouse = expt(day_id).mouse;

    if expt(day_id).multiday_timesincedrug_hours>0
        dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
     else
        dart_str = 'control';
    end
    fn_multi = fullfile(rc.achAnalysis,mouse,['multiday_' dart_str]);
    
    
    load(fullfile(fn_multi,'tc_keep.mat'))
    load(fullfile(fn_multi,'multiday_alignment.mat'))
    load(fullfile(fn_multi,'resp_keep.mat'))
    load(fullfile(fn_multi,'input.mat'))
    load(fullfile(fn_multi,'locomotion.mat'))
    behInput = input;
    clear input;
    %tells the contrast, direction and orientation for each trial each day
    tCon_match = cell(1,nd);
    tDir_match = cell(1,nd);
    tOri_match = cell(1,nd);
    
    %find the contrasts, directions and orientations for each day
    for id = 1:nd
        tCon_match{id} = celleqel2mat_padded(behInput(id).tGratingContrast);
        tDir_match{id} = celleqel2mat_padded(behInput(id).tGratingDirectionDeg);
        tOri_match{id} = tDir_match{id};
        tOri_match{id}(find(tDir_match{id}>=180)) = tDir_match{id}(find(tDir_match{id}>=180))-180;
    end
    oris = unique(tOri_match{post});
    cons = unique(tCon_match{post});
    nOri = length(oris);
    nCon = length(cons);
    
    nOn = behInput(1).nScansOn;
    nOff = behInput(1).nScansOff;
    
   
    trialResp=cell(1,2);
    green_trialResp=cell(2,2);
    red_trialResp=cell(2,2);
    linCellProps = nan(6,4);
    
    for id = 1:nd
    trialResp{id} = mean(data_trial_keep{id}(stimStart:(stimStart+nOn),:,:),1);
    green_trialResp{1,id}=mean(trialResp{id}(:,~RIx{id},green_ind_keep),3);
    green_trialResp{2,id}=mean(trialResp{id}(:,RIx{id},green_ind_keep),3);
    red_trialResp{1,id}=mean(trialResp{id}(:,~RIx{id},red_ind_keep),3);
    red_trialResp{2,id}=mean(trialResp{id}(:,RIx{id},red_ind_keep),3);
    end
    
    
    idx = isnan(red_trialResp{1,pre});
    linfit = polyfit(green_trialResp{1,pre}(~idx),red_trialResp{1,pre}(~idx),1);
    y1 = polyval(linfit,green_trialResp{1,pre});
    [R,p]=corrcoef(green_trialResp{1,pre}(~idx),red_trialResp{1,pre}(~idx)); 
    linCellProps(1,1)=linfit(1); %slope
    linCellProps(2,1)=linfit(2); %intercept
    linCellProps(3,1)=R(2);
    linCellProps(4,1)=p(2);
    linCellProps(5,1)=min(green_trialResp{1,pre});
    linCellProps(6,1)=max(green_trialResp{1,pre});
    
    idx2 = isnan(red_trialResp{1,post});
    linfit = polyfit(green_trialResp{1,post}(~idx2),red_trialResp{1,post}(~idx2),1);
    y2 = polyval(linfit,green_trialResp{1,post});
    
    [R,p]=corrcoef(green_trialResp{1,post}(~idx2),red_trialResp{1,post}(~idx2)); 
    linCellProps(1,2)=linfit(1); %slope
    linCellProps(2,2)=linfit(2); %intercept
    linCellProps(3,2)=R(2);
    linCellProps(4,2)=p(2);
    linCellProps(5,2)=min(green_trialResp{1,post});
    linCellProps(6,2)=max(green_trialResp{1,post});

    
    idx = isnan(red_trialResp{2,pre});
    linfit = polyfit(green_trialResp{2,pre}(~idx),red_trialResp{2,pre}(~idx),1);
    y1 = polyval(linfit,green_trialResp{2,pre});
    [R,p]=corrcoef(green_trialResp{2,pre}(~idx),red_trialResp{2,pre}(~idx)); 
    linCellProps(1,3)=linfit(1); %slope
    linCellProps(2,3)=linfit(2); %intercept
    linCellProps(3,3)=R(2);
    linCellProps(4,3)=p(2);
    linCellProps(5,3)=min(green_trialResp{2,pre});
    linCellProps(6,3)=max(green_trialResp{2,pre});
    

    idx2 = isnan(red_trialResp{2,post});
    linfit = polyfit(green_trialResp{2,post}(~idx2),red_trialResp{2,post}(~idx2),1);
    y2 = polyval(linfit,green_trialResp{2,post});
    [R,p]=corrcoef(green_trialResp{2,post}(~idx2),red_trialResp{2,post}(~idx2)); 
    linCellProps(1,4)=linfit(1); %slope
    linCellProps(2,4)=linfit(2); %intercept
    linCellProps(3,4)=R(2);
    linCellProps(4,4)=p(2);
    linCellProps(5,4)=min(green_trialResp{2,post});
    linCellProps(6,4)=max(green_trialResp{2,post});

linPropsAll{iSess} = linCellProps;
 

% response by condition for cells matched across all conditions
% find cells that I ahve running data for on both days
haveRunning_pre = ~isnan(pref_responses_loc{pre});
haveRunning_post = ~isnan(pref_responses_loc{post});
haveRunning_both = find(haveRunning_pre.* haveRunning_post);
haveRunning_green = intersect(haveRunning_both, green_ind_keep);
haveRunning_red = intersect(haveRunning_both, red_ind_keep);


responseByCond = nan((nCon*2),4);

for iCon = 1:nCon
    if iCon == 1
        counter=1
    else
        counter=counter+2
    end
    
    responseByCond(counter,:)=[mean(pref_responses_stat{pre}(haveRunning_green,iCon), "omitnan") mean(pref_responses_stat{pre}(haveRunning_red,iCon), "omitnan") mean(pref_responses_stat{post}(haveRunning_green,iCon), "omitnan") mean(pref_responses_stat{post}(haveRunning_red,iCon), "omitnan")];
    responseByCond((counter+1),:)=[mean(pref_responses_loc{pre}(haveRunning_green,iCon), "omitnan") mean(pref_responses_loc{pre}(haveRunning_red,iCon), "omitnan") mean(pref_responses_loc{post}(haveRunning_green,iCon), "omitnan") mean(pref_responses_loc{post}(haveRunning_red,iCon), "omitnan")];

end

responseByCondProps = nan(6,2);
linfit = polyfit(responseByCond(:,1),responseByCond(:,2),1);
y1 = polyval(linfit,responseByCond(:,1));
[R,p]=corrcoef(responseByCond(:,1),responseByCond(:,2)); 
responseByCondProps(1,1)=linfit(1); %slope
responseByCondProps(2,1)=linfit(2); %intercept
responseByCondProps(3,1)=R(2);
responseByCondProps(4,1)=p(2);
responseByCondProps(5,1)=min(responseByCond(:,1));
responseByCondProps(6,1)=max(responseByCond(:,1));

linfit = polyfit(responseByCond(:,3),responseByCond(:,4),1);
y2 = polyval(linfit,responseByCond(:,3));
[R,p]=corrcoef(responseByCond(:,3),responseByCond(:,4)); 
responseByCondProps(1,2)=linfit(1); %slope
responseByCondProps(2,2)=linfit(2); %intercept
responseByCondProps(3,2)=R(2);
responseByCondProps(4,2)=p(2);
responseByCondProps(5,2)=min(responseByCond(:,3));
responseByCondProps(6,2)=max(responseByCond(:,3));

respByCondPropsAll{iSess}=responseByCondProps;




end


clearvars -except respByCondPropsAll linPropsAll ds dataStructLabels rc respByCondPropsAll sess_list nSess
%% plotting individual lines and finding mean
figure;
subplot(1,2,1)
meanFit1 = [];
meanFit2 = [];
meanFit3 = [];
meanFit4=[];
meanXRange1=[];
meanXRange2=[];

for iSess = 1:nSess
    fit = [linPropsAll{iSess}(1,1),  linPropsAll{iSess}(2,1)];
    xRange=[linPropsAll{iSess}(5,1), linPropsAll{iSess}(6,1)];
    yValues=polyval(fit,xRange);
    p=plot(xRange,yValues,'LineWidth',.5)
    p.Color=[0,0,0,0.25]
    meanFit1 = [meanFit1;fit];
    meanXRange1=[meanXRange1;xRange];
    hold on
    fit = [linPropsAll{iSess}(1,3),  linPropsAll{iSess}(2,3)];
    xRange=[linPropsAll{iSess}(5,3), linPropsAll{iSess}(6,3)];
    yValues=polyval(fit,xRange);
    p=plot(xRange,yValues,'LineWidth',.5)
    p.Color=[1,0,1,0.25]
    hold on
    meanFit2 = [meanFit2;fit];
    meanXRange2=[meanXRange2;xRange];
end
ylim([-.1,.25])
xlim([-.1,.25])
axis square
xRange1=mean(meanXRange1,1);
xRange2=mean(meanXRange2,1);
fit1=mean(meanFit1,1);
fit2=mean(meanFit2,1);
yValues = polyval(fit1,xRange1);
plot(xRange1,yValues,'k','LineWidth',1.5)

yValues = polyval(fit2,xRange2);
plot(xRange2,yValues,'m','LineWidth',1.5)
set(gca, 'TickDir', 'out','Box', 'off');
title("Pre-DART")
xlabel("Putative Pyr")
ylabel("HT+ SOM")


meanXRange1=[];
meanXRange2=[];

subplot(1,2,2)
for iSess = 1:nSess
    fit = [linPropsAll{iSess}(1,2),  linPropsAll{iSess}(2,2)];
    xRange=[linPropsAll{iSess}(5,2), linPropsAll{iSess}(6,2)];
    yValues=polyval(fit,xRange);
    p=plot(xRange,yValues,'LineWidth',.5)
    p.Color=[0,0,0,0.25]
    meanFit3 = [meanFit3;fit];
    meanXRange1=[meanXRange1;xRange];
    hold on
    fit = [linPropsAll{iSess}(1,4),  linPropsAll{iSess}(2,4)];
    xRange=[linPropsAll{iSess}(5,4), linPropsAll{iSess}(6,4)];
    yValues=polyval(fit,xRange);
    p=plot(xRange,yValues,'LineWidth',.5)
    p.Color=[1,0,1,0.25]
    hold on
    meanFit4 = [meanFit4;fit];
    meanXRange2=[meanXRange2;xRange];
end
ylim([-.1,.25])
xlim([-.1,.25])
axis square
xRange1=mean(meanXRange1,1);
xRange2=mean(meanXRange2,1);
fit3=mean(meanFit3,1);
fit4=mean(meanFit4,1);
yValues = polyval(fit3,xRange1);
plot(xRange1,yValues,'k','LineWidth',1.5)

yValues = polyval(fit4,xRange2);
plot(xRange2,yValues,'m','LineWidth',1.5)
set(gca, 'TickDir', 'out','Box', 'off');
title("Post-DART")
xlabel("Putative Pyr")
ylabel("HT+ SOM")


%% plotting individual lines and finding mean
figure;
subplot(1,2,1)
meanFit1 = [];
meanFit2 = [];
meanFit3 = [];
meanFit4=[];
meanXRange1=[];
meanXRange2=[];

for iSess = 1:nSess
    fit = [linPropsAll{iSess}(1,1),  linPropsAll{iSess}(2,1)];
    xRange=[linPropsAll{iSess}(5,1), linPropsAll{iSess}(6,1)];
    yValues=polyval(fit,xRange);
    p=plot(xRange,yValues,'LineWidth',.25)
    p.Color=[0,0,0,0.25]
    meanFit1 = [meanFit1;fit];
    meanXRange1=[meanXRange1;xRange];
    hold on
    fit = [linPropsAll{iSess}(1,2),  linPropsAll{iSess}(2,2)];
    xRange=[linPropsAll{iSess}(5,2), linPropsAll{iSess}(6,2)];
    yValues=polyval(fit,xRange);
    p=plot(xRange,yValues,'LineWidth',.25)
    p.Color=[0,0,1,0.25]
    hold on
    meanFit2 = [meanFit2;fit];
    meanXRange2=[meanXRange2;xRange];
end
ylim([-.1,.25])
xlim([-.1,.25])
axis square
xRange1=mean(meanXRange1,1);
xRange2=mean(meanXRange2,1);
fit1=mean(meanFit1,1);
fit2=mean(meanFit2,1);
yValues = polyval(fit1,xRange1);
plot(xRange1,yValues,'k','LineWidth',1.5)

yValues = polyval(fit2,xRange2);
plot(xRange2,yValues,'b','LineWidth',1.5)
set(gca, 'TickDir', 'out','Box', 'off');
title("Pre-DART")
xlabel("Putative Pyr")
ylabel("HT+ SOM")


meanXRange1=[];
meanXRange2=[];
subplot(1,2,2)
for iSess = 1:nSess
    fit = [linPropsAll{iSess}(1,3),  linPropsAll{iSess}(2,3)];
    xRange=[linPropsAll{iSess}(5,3), linPropsAll{iSess}(6,3)];
    yValues=polyval(fit,xRange);
    p=plot(xRange,yValues,'LineWidth',.25)
    p.Color=[0,0,0,0.25]
    meanFit3 = [meanFit3;fit];
    meanXRange1=[meanXRange1;xRange];
    hold on
    fit = [linPropsAll{iSess}(1,4),  linPropsAll{iSess}(2,4)];
    xRange=[linPropsAll{iSess}(5,4), linPropsAll{iSess}(6,4)];
    yValues=polyval(fit,xRange);
    p=plot(xRange,yValues,'LineWidth',.25)
    p.Color=[0,0,1,0.25]
    hold on
    meanFit4 = [meanFit4;fit];
    meanXRange2=[meanXRange2;xRange];
end
ylim([-.1,.25])
xlim([-.1,.25])
axis square
xRange1=mean(meanXRange1,1);
xRange2=mean(meanXRange2,1);
fit3=mean(meanFit3,1);
fit4=mean(meanFit4,1);
yValues = polyval(fit3,xRange1);
plot(xRange1,yValues,'k','LineWidth',1.5)

yValues = polyval(fit4,xRange2);
plot(xRange2,yValues,'b','LineWidth',1.5)
set(gca, 'TickDir', 'out','Box', 'off');
title("Running")
xlabel("Putative Pyr")
ylabel("HT+ SOM")
%% alternate version without running and stationary split
sess_list = [131 133 138 142 163 171];%enter all the sessions you want to concatenate
nSess=length(sess_list);
%% loop through sessions

linPropsAll = cell(1,nSess);
respByCondPropsAll=cell(1,nSess);
for iSess = 1:nSess
    day_id = sess_list(iSess)



    ds = 'DART_V1_contrast_ori_Celine'; %dataset info
    dataStructLabels = {'contrastxori'};
    
    rc = behavConstsDART; %directories
    eval(ds);
        
    
    pre_day = expt(day_id).multiday_matchdays;
    
    nd=2; %hardcoding the number of days for now
    
    % INDICATE THE PRE VS. POST DAYS DEPENDING ON THE ORDER OF MATCHING
    pre=2;
    post=1;
    
    mouse = expt(day_id).mouse;

    if expt(day_id).multiday_timesincedrug_hours>0
        dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
     else
        dart_str = 'control';
    end
    fn_multi = fullfile(rc.achAnalysis,mouse,['multiday_' dart_str]);
    
    
    load(fullfile(fn_multi,'tc_keep.mat'))
    load(fullfile(fn_multi,'multiday_alignment.mat'))
    load(fullfile(fn_multi,'resp_keep.mat'))
    load(fullfile(fn_multi,'input.mat'))
    load(fullfile(fn_multi,'locomotion.mat'))
    behInput = input;
    clear input;
    %tells the contrast, direction and orientation for each trial each day
    tCon_match = cell(1,nd);
    tDir_match = cell(1,nd);
    tOri_match = cell(1,nd);
    
    %find the contrasts, directions and orientations for each day
    for id = 1:nd
        tCon_match{id} = celleqel2mat_padded(behInput(id).tGratingContrast);
        tDir_match{id} = celleqel2mat_padded(behInput(id).tGratingDirectionDeg);
        tOri_match{id} = tDir_match{id};
        tOri_match{id}(find(tDir_match{id}>=180)) = tDir_match{id}(find(tDir_match{id}>=180))-180;
    end
    oris = unique(tOri_match{post});
    cons = unique(tCon_match{post});
    nOri = length(oris);
    nCon = length(cons);
    
    nOn = behInput(1).nScansOn;
    nOff = behInput(1).nScansOff;
    
   
    trialResp=cell(1,2);
    green_trialResp=cell(1,2);
    red_trialResp=cell(1,2);
    linCellProps = nan(6,2);
    
    for id = 1:nd
    trialResp{id} = mean(data_trial_keep{id}(stimStart:(stimStart+nOn),:,:),1);
    green_trialResp{1,id}=mean(trialResp{id}(:,:,green_ind_keep),3);
    red_trialResp{1,id}=mean(trialResp{id}(:,:,red_ind_keep),3);
    end
    
    
    idx = isnan(red_trialResp{1,pre});
    linfit = polyfit(green_trialResp{1,pre}(~idx),red_trialResp{1,pre}(~idx),1);
    y1 = polyval(linfit,green_trialResp{1,pre});
    [R,p]=corrcoef(green_trialResp{1,pre}(~idx),red_trialResp{1,pre}(~idx)); 
    linCellProps(1,1)=linfit(1); %slope
    linCellProps(2,1)=linfit(2); %intercept
    linCellProps(3,1)=R(2);
    linCellProps(4,1)=p(2);
    linCellProps(5,1)=min(green_trialResp{1,pre});
    linCellProps(6,1)=max(green_trialResp{1,pre});
    
    idx2 = isnan(red_trialResp{1,post});
    linfit = polyfit(green_trialResp{1,post}(~idx2),red_trialResp{1,post}(~idx2),1);
    y2 = polyval(linfit,green_trialResp{1,post});
    
    [R,p]=corrcoef(green_trialResp{1,post}(~idx2),red_trialResp{1,post}(~idx2)); 
    linCellProps(1,2)=linfit(1); %slope
    linCellProps(2,2)=linfit(2); %intercept
    linCellProps(3,2)=R(2);
    linCellProps(4,2)=p(2);
    linCellProps(5,2)=min(green_trialResp{1,post});
    linCellProps(6,2)=max(green_trialResp{1,post});

    

linPropsAll{iSess} = linCellProps;
 


end


clearvars -except respByCondPropsAll linPropsAll ds dataStructLabels rc respByCondPropsAll sess_list nSess
%% plotting individual lines and finding mean
figure;
meanFit1 = [];
meanFit2 = [];

meanXRange1=[];
meanXRange2=[];

for iSess = 2:nSess
    fit = [linPropsAll{iSess}(1,1),  linPropsAll{iSess}(2,1)];
    xRange=[linPropsAll{iSess}(5,1), linPropsAll{iSess}(6,1)];
    yValues=polyval(fit,xRange);
    p=plot(xRange,yValues,'LineWidth',.5)
    p.Color=[0,0,0,0.25]
    meanFit1 = [meanFit1;fit];
    meanXRange1=[meanXRange1;xRange];
    hold on
    fit = [linPropsAll{iSess}(1,2),  linPropsAll{iSess}(2,2)];
    xRange=[linPropsAll{iSess}(5,2), linPropsAll{iSess}(6,2)];
    yValues=polyval(fit,xRange);
    p=plot(xRange,yValues,'LineWidth',.5)
    p.Color=[0,0,1,0.25]
    hold on
    meanFit2 = [meanFit2;fit];
    meanXRange2=[meanXRange2;xRange];
end
ylim([-.1,.35])
xlim([-.1,.35])
axis square
xRange1=mean(meanXRange1,1);
xRange2=mean(meanXRange2,1);
fit1=mean(meanFit1,1);
fit2=mean(meanFit2,1);
yValues = polyval(fit1,xRange1);
plot(xRange1,yValues,'k','LineWidth',1.5)

yValues = polyval(fit2,xRange2);
plot(xRange2,yValues,'b','LineWidth',1.5)
set(gca, 'TickDir', 'out','Box', 'off');
title("SOM:Pyr relationship")
xlabel("Putative Pyr")
ylabel("HT+ SOM")



