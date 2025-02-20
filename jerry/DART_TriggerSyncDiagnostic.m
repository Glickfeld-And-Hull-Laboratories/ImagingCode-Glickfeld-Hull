function results = DART_TriggerSyncDiagnostic(person,session,saveplots)
%DART_TriggerSyncDiagnostic returns locomotion status in the given experiment sessions. 
% results = DART_TriggerSyncDiagnostic(person,session)
% PERSON is the intial of the experimenter given as string (e.g. "CC" or
% "TH"). Valid experimenters are defined in the function, edit function to add.
% SESSION is the sessions to be checked given as a row vector of doubles,
% defined by each experimenter's DART dataset info file.
% SAVEPLOTS takes value 1 or 0. Determines whether plots are saved or just
% displayed.


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
disp(['Dataset name is ' ds ' and is hardcoded in function. Edit function if file name changed.'])

nSesh = length(session);

% nRunningTrialsInConditions = cell(nSesh,2);
% wheelspeeds = cell(nSesh,2);

for sesh = 1:length(session)
    mouse = expt(session(sesh)).mouse;
    day = expt(session(sesh)).date;
    run = expt(session(sesh)).contrastxori_runs;
    run = run{1,1};
    time = expt(session(sesh)).contrastxori_time;
    time = time{1,1};
    ExperimentFolder = expt(session(sesh)).exptType;
    load([dat_root '\' mouse '\' day '\' run '\' run '_000_000.mat']);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\' 'data-' mouse '-' day '-' time '.mat'];

    fn = fullfile(rc.achAnalysis,ExperimentFolder,mouse,day);
    fnout = fullfile(fn,run);

    load(fName);

    % thisRunningTrialsInConditions = cell(tot_conds,2);
    % thisRunningTrialsInConditions(:,1) = all_conds;
    % thisRunningTrialsInConditions(:,2) = num2cell(RunTrialsN);
    % nRunningTrialsInConditions{sesh,1} = [mouse '-' day];
    % nRunningTrialsInConditions{sesh,2} = thisRunningTrialsInConditions;
    % wheelspeeds{sesh,1} = [mouse '-' day];
    % wheelspeeds{sesh,2} = wheelspd;
end


end
