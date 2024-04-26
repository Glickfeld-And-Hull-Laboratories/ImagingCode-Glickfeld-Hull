
nCon = length(targetCon)

mice=[];
red_concat=[];
green_concat=[];
nKeep_concat=[];

tc_trial_avrg_stat_concat=cell(1,nd);
tc_trial_avrg_loc_concat=cell(1,nd);
conBySize_resp_stat_concat=cell(1,nd);
conBySize_resp_loc_concat=cell(1,nd);
h_concat=cell(1,nd);
data_resp_concat=cell(1,nd);

tc_trial_avrg_stat_largePupil_concat=cell(1,nd);
tc_trial_avrg_stat_smallPupil_concat=cell(1,nd);
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
dfof_max_diff_concat=[];
red_fluor_concat=[];
green_fluor_concat=[];
wheel_corr_concat=cell(1,nd);
meanF_concat=cell(1,nd);
norm_dir_resp_stat_concat = cell(1,nd);
norm_dir_resp_loc_concat = cell(1,nd);
pref_nonPref_stat_concat=cell(1,nd);
pref_nonPref_loc_concat=cell(1,nd);
pref_dir_concat=cell(1,nd);
noiseCorr_concat = cell(1,nd);
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
    thisSess = sess_list(iSess);
    mouse = expt(thisSess).mouse;
    mice=[mice;mouse];
    thisDrug = expt(thisSess).drug;
    drug{iSess}=thisDrug;
    

    if iSess > 1
        cellID_adjustment=max(temp_table.cell_ID_unique); %this should get saved until the next loop;
    end

    if expt(thisSess).multiday_timesincedrug_hours>0
        dart_str = [expt(thisSess).drug '_' num2str(expt(thisSess).multiday_timesincedrug_hours) 'Hr'];
    else
        dart_str = 'control';
    end
    fn_multi = fullfile(rc.achAnalysis,mouse,['multiday_' dart_str]);

    load(fullfile(fn_multi,'tc_keep.mat'));
    load(fullfile(fn_multi,'resp_keep.mat'));
    load(fullfile(fn_multi,'input.mat'));
    %load(fullfile(fn_multi,'locomotion.mat'));
    load(fullfile(fn_multi,'fluor_intensity.mat'));
    % load(fullfile(fn_multi,'HT_pyr_relationship.mat'));
    load(fullfile(fn_multi,'pupilMeans.mat'));

    nKeep = size(tc_trial_avrg_stat{post},2);

    pupilMeans_concat(:,:,iSess)=pupilMeans;
    motorByPupil_concat(:,:,iSess)=motorByPupil;
%    pupilCounts_concat(:,:,iSess)=pupilCounts;


    %tells the contrast, direction and orientation for each trial each day
    tCon_match = cell(1,nd);
    tDir_match = cell(1,nd);
    tOri_match = cell(1,nd);
    tSize_match = cell(1,nd);

    %find the contrasts, directions and orientations for each day
    for id = 1:nd
        tCon_match{id} = celleqel2mat_padded(input(id).tGratingContrast(1:nTrials(id)));
        tDir_match{id} = celleqel2mat_padded(input(id).tGratingDirectionDeg(1:nTrials(id)));
        tOri_match{id} = tDir_match{id};
        tOri_match{id}(find(tDir_match{id}>=180)) = tDir_match{id}(find(tDir_match{id}>=180))-180;
        tSize_match{id} = celleqel2mat_padded(input(id).tGratingDiameterDeg(1:nTrials(id)));
    end
    dirs = unique(tDir_match{post});
    cons = unique(tCon_match{post});
    sharedCon=find(ismember(cons, targetCon));
    sizes = unique(tSize_match{1});

    nOn = input(1).nScansOn;
    nOff = input(1).nScansOff;
    %start conatenating
    dirs_concat = [dirs_concat,dirs]; 
    cons_concat = [cons_concat,cons(sharedCon)];
    red_concat = [red_concat, red_keep_logical];
    green_concat = [green_concat, green_keep_logical];
    nKeep_concat = [nKeep_concat,nKeep];
    % responseByCondProps_concat(:,:,iSess)=responseByCondProps;

    clear cons
    
    
    for id = 1:nd
        
        tc_trial_avrg_stat_concat{id} =cat(2,tc_trial_avrg_stat_concat{id},tc_trial_avrg_stat{id}(:,:,sharedCon,:));
        tc_trial_avrg_stat_largePupil_concat{id} = cat(2,tc_trial_avrg_stat_largePupil_concat{id},tc_trial_avrg_stat_largePupil{id}(:,:,sharedCon,:));
        tc_trial_avrg_stat_smallPupil_concat{id} = cat(2,tc_trial_avrg_stat_smallPupil_concat{id},tc_trial_avrg_stat_smallPupil{id}(:,:,sharedCon,:));
        tc_trial_avrg_loc_concat{id} =cat(2,tc_trial_avrg_loc_concat{id},tc_trial_avrg_loc{id}(:,:,sharedCon,:));
        nonPref_trial_avrg_stat_concat{id} =cat(2,nonPref_trial_avrg_stat_concat{id},nonPref_trial_avrg_stat{id}(:,:,sharedCon,:));
        nonPref_trial_avrg_loc_concat{id} =cat(2,nonPref_trial_avrg_loc_concat{id},nonPref_trial_avrg_loc{id}(:,:,sharedCon,:));
        resp_keep_concat{id}=cat(1,resp_keep_concat{id},resp_keep{id});
        % resp_max_keep_concat{id}=cat(1,resp_max_keep_concat{id},resp_max_keep{id}(:,sharedCon,:));
        pref_responses_loc_concat{id}=cat(1,pref_responses_loc_concat{id},pref_responses_loc{id}(:,sharedCon,:));
        pref_responses_stat_concat{id}=cat(1,pref_responses_stat_concat{id},pref_responses_stat{id}(:,sharedCon,:));
        pref_peak_stat_concat{id}=cat(1,pref_peak_stat_concat{id},pref_peak_stat{id}(:,sharedCon,:));
        pref_peak_loc_concat{id}=cat(1,pref_peak_loc_concat{id},pref_peak_loc{id}(:,sharedCon,:));
        pref_responses_stat_largePupil_concat{id}=cat(1,pref_responses_stat_largePupil_concat{id},pref_responses_stat_largePupil{id}(:,sharedCon,:));
        pref_responses_stat_smallPupil_concat{id}=cat(1,pref_responses_stat_smallPupil_concat{id},pref_responses_stat_smallPupil{id}(:,sharedCon,:));
        % RIx_concat{id}=cat(1,RIx_concat{id},sum(RIx{id}));
        % wheel_corr_concat{id}=cat(2,wheel_corr_concat{id},wheel_corr{id});
        % meanF=mean(fullTC_keep{id},1);
        % meanF_concat{id}=cat(2,meanF_concat{id}, meanF);
        % norm_dir_resp_stat_concat{id}=cat(1,norm_dir_resp_stat_concat{id},norm_dir_resp_stat{id});
        % norm_dir_resp_loc_concat{id}=cat(1,norm_dir_resp_loc_concat{id},norm_dir_resp_loc{id});
        % pref_dir_concat{id}=cat(2,pref_dir_concat{id},pref_dir_keep{id});
        % noiseCorr_concat{id}=cat(2,noiseCorr_concat{id},noiseCorr{id});
        % sigCorr_concat{id}=cat(2,sigCorr_concat{id},sigCorr{id});
        for i = 1:length(sharedCon)
            iCon=sharedCon(i);
            for iSize = 1:length(sizes)
                pref_allTrials_stat_concat{i,id}=[pref_allTrials_stat_concat{i,id},pref_allTrials_stat{iCon,iSize,id}];
                pref_allTrials_loc_concat{i,id}=[pref_allTrials_loc_concat{i,id},pref_allTrials_loc{iCon,iSize,id}];
                pref_allTrials_largePupil_concat{i,id}=[pref_allTrials_largePupil_concat{i,id},pref_allTrials_largePupil{iCon,iSize,id}];
                pref_allTrials_smallPupil_concat{i,id}=[pref_allTrials_smallPupil_concat{i,id},pref_allTrials_smallPupil{iCon,iSize,id}];
            end
        end
        clear meanF i
    end
    % dfof_max_diff_concat=cat(1,dfof_max_diff_concat,dfof_max_diff(:,sharedCon,:));
   green_fluor_concat=cat(2,green_fluor_concat,green_fluor_keep);
   red_fluor_concat=cat(2,red_fluor_concat,red_fluor_keep);
    
iSess
end