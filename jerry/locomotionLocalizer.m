function [wheelspeeds,nRunningTrialsInConditions] = locomotionLocalizer(person,session)
%LOCOMOTIONLOCALIZER returns locomotion status in the given experiment sessions. 
% [wheelspeeds,nRunningTrialsInConditions] = locomotionLocalizer(person,session)
% PERSON is the intial of the experimenter given as string (e.g. "CC" or
% "TH"). Valid experimenters are defined in the function, edit function to add.
% SESSION is the sessions to be checked given as a row vector of doubles,
% defined by each experimenter's DART dataset info file.
% WHEELSPEEDS is an nSessions x 2 cell array that has mouse-day in the
% first column and nTrials x 1 double arrays with trial average wheelspeeds
% in the second column.
% NRUNNINGTRIALSINCONDITIONS is nSessions x 2 cell array that has mouse-day in the
% first column and nested cell arrays in the second column. The nested arrays
% have dimensions nStimConditions x 2 and have all stimulus conditions in
% the first column and number of running trials in that condition in the
% second column.

if isa(session,'double') == 0
    error('Invalid session class, supply session as a row vector of doubles.')
end

experimenters_array = ["TH" "CC"];

if ismember(person,experimenters_array) == false
    error('Invalid experimenter input, use either TH or CC as string ("CC" or "TH").')
end

isilon_home = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\';
dat_root = [isilon_home 'ACh\Data\2p_data'];
analysis_root = [isilon_home 'ACh\Analysis\2p_analysis'];
rc = behavConstsDART;
if person == "TH"
    ds = 'DART_expt_info';
    code_cd = [isilon_home 'jerry\repositories\ImagingCode-Glickfeld-Hull\jerry'];
elseif person == "CC"
    ds = 'DART_V1_atropine_Celine';
    code_cd = [isilon_home 'celine\Repositories\ImagingCode-Glickfeld-Hull\celine'];
end
cd(code_cd);
eval(ds);
disp(['Dataset name is ' ds ' and is hardcoded in function. Edit function if file name changed.']);


nSesh = length(session);

nRunningTrialsInConditions = cell(nSesh,2);
wheelspeeds = cell(nSesh,2);

for sesh = 1:length(session)
    mouse = expt(session(sesh)).mouse;
    day = expt(session(sesh)).date;
    run = expt(session(sesh)).contrastxori_runs;
    run = run{1,1};
    time = expt(session(sesh)).contrastxori_time;
    time = time{1,1};
    ExperimentFolder = expt(session(sesh)).exptType;
    %load([dat_root '\' mouse '\' day '\' run '\' run '_000_000.mat']);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\' 'data-' mouse '-' day '-' time '.mat'];
    expt(53).exptType = 'PV_CMPDA';
    fn = fullfile(rc.achAnalysis,ExperimentFolder,mouse,day);
    fnout = fullfile(fn,run);

    load(fName);
    %plot trial running speed
    loc_dat = input.wheelSpeedValues;
    wheelspd = zeros(length(loc_dat),1);
    for iTrial = 1:length(wheelspd)
        wheelspd(iTrial) = mean(loc_dat{iTrial});
    end
    figure
    plot(wheelspd);
    sgtitle('Trial Average Wheel Speed')
    print(fullfile(fnout,'wheelspd.pdf'),'-dpdf','-bestfit');

    %plot # of running trials per condition
    runidx = find(wheelspd>2);
    runTrialBoolean = wheelspd > 2;
    tCon = cell2mat(input.tGratingContrast);
    contrasts = unique(tCon);
    nCon = length(contrasts);
    nTrials = length(wheelspd);

    tDir = celleqel2mat_padded(input.tGratingDirectionDeg(1:nTrials));
    dirs = unique(tDir);
    nDir = length(dirs);

    tSize = celleqel2mat_padded(input.tGratingDiameterDeg(1:nTrials));
    sizes = unique(tSize);
    nSize = length(sizes);

    tot_conds = nSize * nDir * nCon;
    all_conds = cell(tot_conds,1);
    RunTrialsN = nan(tot_conds,1);
    cond_counter = 1;
    %all_idx = [];
    
    for iSize = 1:nSize
        ind_size = find(tSize == sizes(iSize));
        for iDir = 1:nDir
            ind_dir = find(tDir == dirs(iDir));
            for iCon = 1:nCon
                ind_con = find(tCon == contrasts(iCon));
                this_ind = intersect(intersect(ind_size,ind_dir,'stable'),ind_con,'stable'); 
                this_condition = ['Con-' num2str(contrasts(iCon)) '-Dir-' num2str(dirs(iDir)) '-Size-' num2str(sizes(iSize))];
                all_conds{cond_counter} = this_condition;
                haveRun = intersect(this_ind,runidx);
                RunTrialsN(cond_counter) = length(haveRun);
                cond_counter = cond_counter + 1;
                %all_idx = [all_idx length(this_ind)];
            end
        end
    end
    
    figure
    bar(RunTrialsN)
    sgtitle([mouse '-' day '-Number of Running Trials in Each Stimulus Condition'])
    ylim([0 max(RunTrialsN+1)])
    yticks([0:1:max(RunTrialsN+1)])
    xticks([1:1:length(all_conds)])
    xticklabels(all_conds)
    ylabel('# of Running Trials')
    
    ax = gca;
    ax.FontSize = 8; 
    print(fullfile(fnout,'nRunTrialsInCond.pdf'),'-dpdf','-bestfit');

    thisRunningTrialsInConditions = cell(tot_conds,2);
    thisRunningTrialsInConditions(:,1) = all_conds;
    thisRunningTrialsInConditions(:,2) = num2cell(RunTrialsN);
    nRunningTrialsInConditions{sesh,1} = [mouse '-' day];
    nRunningTrialsInConditions{sesh,2} = thisRunningTrialsInConditions;
    wheelspeeds{sesh,1} = [mouse '-' day];
    wheelspeeds{sesh,2} = wheelspd;
end



end
