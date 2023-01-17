function [linProps] = trialRelationship(day_id)
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
    linfit1 = polyfit(green_trialResp{2,post}(~idx2),red_trialResp{2,post}(~idx2),1);
    y2 = polyval(linfit,green_trialResp{2,post});
    [R,p]=corrcoef(green_trialResp{2,post}(~idx2),red_trialResp{2,post}(~idx2)); 
    linCellProps(1,4)=linfit(1); %slope
    linCellProps(2,4)=linfit(2); %intercept
    linCellProps(3,4)=R(2);
    linCellProps(4,4)=p(2);
    linCellProps(5,4)=min(green_trialResp{2,post});
    linCellProps(6,4)=max(green_trialResp{2,post});

linProps = linCellProps;

end