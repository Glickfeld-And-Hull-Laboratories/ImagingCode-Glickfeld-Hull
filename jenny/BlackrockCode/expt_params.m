
% for establishing the following info from recordings during experiments
    % probe/headstage combination
    % genotype
    % experiment list
    % spike sorted file name
    
%%
clear params
% where to save the info for this experiment
save_path = 'Z:\All_staff\home\jen\Analysis\Spyking_circus\';

probeLog_path = 'Z:\All_staff\home\jen\Notebook files\probeLog.xlsx';

params.info.recDate = '230420';
params.info.probeType = 2;
% types =
    % (1) H4/optrode
    % (2) A32/optrode
    % (3) A32
params.info.headstage = 2;
% types = 
    % (1) old blackrock (pre 12/2021)
    % (2) new blackrock same side omnetics  
    % (3) new blackrock opp side omnetics (uncommon)
    
params.info.Genotype = 'PVxAi32';
% types =
    % 'Scnn1axAi32'
    % 'SOMxAi32'
% params.info.expt_list = {'oriTuning';'twoStim_multiContrast';'eOPN_ctrl'};
% params.info.expt_list = {'flashingBar';'sizeTuning_1secOn';'sizeCon_100msOn';'visStimRet';'optoTag'};
% params.info.expt_list = {'flashingBar';'sizeTuning_100msOn';'sizeTuning_1secOn'};

params.info.expt_list = {'oriTuning';'twoStim_PVChR2'};

% params.info.expt_list = {'oriTuning';'laserOnly'};
% params.info.expt_list = {'oriTuning';'twoStim_multiContrast';'twoStim_L4stim_pre';'driftingGratings_pre';'driftingGratings_postDrug';'optoTag'};
% params.info.expt_list = {'oriTuning';'twoStim_multiContrast';'twoStim_L4stim_pre'};
% params.info.expt_list = {'oriTuning';'twoStim_L23stim_100';'twoStim_L23stim_250';'twoStim_L23stim_1000'};
% params.info.expt_list = {'dirTuning_25Contrast';'sizeTuning_100msOn_80Con';'sizeTuning_100msOn_drugOn';'dirTuning_25Contrast_postDrug';'optoTag'};

% 
% types = 
    % 'flashingBar'
    % 'sizeTuning_1secOn';
    % 'oriTuning'
    % 'optoTag'
    % 'sizeCon_100msOn'
    % 'twoStim_L4stim_pre'
    % 'twoStim_L4stim_pre_shortISI'
    % 'twoStim_L4stim_peri'
    % 'twoStim_phaseShift'
    % 'eOpnStim'
    
%% wrap and save experiment info, write date to probeLog 

wrapExptParams(params,save_path,probeLog_path);

%%
function wrapExptParams(params,save_path,probeLog_path)
probeList = {'H4-optrode';'A32-optrode';'A32';};
hsList = {'old Blackrock';'new Blackrock (same side)';'old Blackrock (opposite side)'};

params.info.probeType = probeList{params.info.probeType};
params.info.headstage = hsList{params.info.headstage};

params.info.spikeFiles = ['datafile001'];
if numel(params.info.expt_list)>1
for file_i = 2:numel(params.info.expt_list)
    params.info.spikeFiles = [params.info.spikeFiles '-',num2str(file_i)];
end
end

temp_table = readtable(probeLog_path,'Sheet',1);
rowID = find(cellfun(@(x) strcmp(x,params.info.probeType),temp_table.Probe).*cellfun(@(x) strcmp(x,params.info.headstage),temp_table.HS));
if ~contains(temp_table.Dates{rowID},params.info.recDate)
try
    xlswrite(probeLog_path,...
    {['''' temp_table.Dates{rowID},',',params.info.recDate]},1,...
    ['F' num2str(rowID + 1)]);
    
catch
    input('close the excel file. press enter when done.');
    xlswrite(probeLog_path,...
    {['''' temp_table.Dates{rowID},',',params.info.recDate]},1,...
    ['F' num2str(rowID + 1)]);    
end
end

if ~exist([save_path,params.info.recDate])
    mkdir([save_path,params.info.recDate])
end
save(fullfile(save_path,params.info.recDate,[params.info.recDate,'_',params.info.spikeFiles,'_exptInfo.mat']))

disp(['wrapped experiment info for ',params.info.recDate, ' :)']);
end
