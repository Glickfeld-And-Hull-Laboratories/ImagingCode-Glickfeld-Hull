function make_params_JL(rec_date,varargin)

%% set expt variables
p = inputParser;
p.addParameter('filename', [], @ischar);

% parse inputs
p.parse(varargin{:});
inputParam = p.Results;

analysis_path = ['Z:\All_Staff\home\jen\Analysis\Spyking_Circus\', rec_date,'\'];
rawData_path = ['Z:\All_Staff\home\jen\Raw Data\Multiunit\',rec_date,'\'];

%% load experiment info 

fs = dir(fullfile(analysis_path, '*exptInfo.mat'));

if isempty(inputParam.filename) && ~isempty(fs) % if experiment record exists 
    load(fullfile(analysis_path, fs(1).name));
    spikeSortFile = params.info.spikeFiles;
    autoFill = true;
    
elseif isempty(inputParam.filename) && isempty(fs) % if record doesnt exist and no filename input, request filename 
    spikeSortFile = input('spike sorted filename?','s');
    disp('need to manually edit probeLog to include headstage info');
    autoFill = false;
else
    spikeSortFile = inputParam.filename; % if record doesnt exist and filename was provided, use this 
    autoFill = false;
end

% get filenames for spike sorted data 
if ~contains(spikeSortFile,'-')
    disp('found 1 file')
    layersFile = spikeSortFile;
    temp_files = {spikeSortFile};
else
    temp_files = cellfun(@(x) ['datafile00',x(end)],strsplit(spikeSortFile,'-'),'un',0);
    disp(['found ',num2str(numel(temp_files)),' files']);
    disp(temp_files);
end

%% fill out recording experiment data

if ~exist([analysis_path,rec_date,'_',spikeSortFile,'_metaData.mat']) % if a params file wasn't already created prior 

% if metadata doesn't exist yet, fill it out either manually or through the
% info in experiment info 
if autoFill
    exptType = params.info.expt_list'; 
else 

    for file_i = 1:numel(temp_files)

    exptType{file_i} = input(['which experiment type is ',temp_files{file_i},'? previous ones include: \n'...
        'sizeCon \n'...
        'flashingBar \n' ...
        'optoTag \n' ...
        'twoStim_L4stim \n' ...
        'twoStim_phaseShift \n' ...
        'twoStim_multiContrast \n' ...
        'oriTuning \n '],'s');
    end
end

params.info.spikesFile = [analysis_path,rec_date,'_',spikeSortFile,'_spikes.mat'];
params.info.triggersFile = [analysis_path,rec_date,'_',spikeSortFile,'_triggers.mat'];
params.info.layersFile = [analysis_path,rec_date,'_',spikeSortFile,'_layerBound.mat'];
params.info.recDate = rec_date;
params.info.experimentTypes = unique(exptType);

for exp_i = 1:numel(params.info.experimentTypes)
    sub_idx = cell2mat(cellfun(@(x) strcmp(x,params.info.experimentTypes{exp_i}),exptType,'un',0));
    params.info.([params.info.experimentTypes{exp_i},'Idx']) = sub_idx;
    params.info.([params.info.experimentTypes{exp_i},'InputFile']) = cellfun(@(x) [analysis_path,rec_date,'_',x,'_input.mat'],temp_files(sub_idx),'un',0);
end

else
    load([analysis_path,rec_date,'_',spikeSortFile,'_metaData.mat']);
end

% make spikes if not done yet
if ~exist(params.info.spikesFile)
    if ~exist([analysis_path,rec_date,'_',spikeSortFile,'_metaData.mat'])
    disp('running get units...must use same datafile names as indicated earlier');
    SS_type = input('Kilosort (1) or spyking-circus (2)?');
    if SS_type == 1
        loadPhyOutput
    elseif SS_type ==2 
        GetUnits_JL
    end
    else
        params.info.spikesFile = [analysis_path,rec_date,'_',spikeSortFile,'_spikes.mat'];
    end
end
    
if ~exist(params.info.layersFile)
    if ~exist([analysis_path,rec_date,'_',spikeSortFile,'_layerBound.mat']) % run layer analysis if layerBound doesn't exist yet 
        if numel(temp_files)>1
            disp(temp_files);
            layersFile = temp_files{input('which index for layers? \n')};
        end
        try
            LFPanalysis([rawData_path,layersFile,'.ns4'],params.info.triggersFile); 
        catch
            disp('could not complete LFP analysis initially...is this one of KF recordings?')
            waitFlag = input('y/n','s');
            if strcmp(waitFlag,'y')
            LFPanalysis(['Z:\All_staff\home\jen\Raw Data\Multiunit\KF\211011\',layersFile,'.ns4'],params.info.triggersFile);
            end
            
        end
    else
        params.info.layersFile = [analysis_path,rec_date,'_',spikeSortFile,'_layerBound.mat'];
    end
end
    
    load(params.info.layersFile); load(params.info.spikesFile);

% extract unit waveforms 
if size(spikes.template_big,1) == 82 % phy output
    params.waveform = spikes.template_big(22:end-11,:)';
    params.pt_ms = spikes.wave_stats.pt_ms/30;
else % spyking circus output 
    params.waveform = spikes.template_big(12:end,:)';
    params.pt_ms = spikes.wave_stats.pt_ms;
end

try
    temp_depth = spikes.template_pos(:,2); % two column probe 
catch
    temp_depth = spikes.template_pos; % one column probe
end

params.cellDepths = temp_depth;
params.L23_cells = arrayfun(@(x) x>layerBound.layer23_depth(1) && x<layerBound.layer23_depth(2),params.cellDepths);
params.L4_cells = arrayfun(@(x) x>layerBound.layer4_depth(1) && x<layerBound.layer4_depth(2),params.cellDepths);
params.L56_cells = params.cellDepths >= layerBound.layer4_depth(2);

disp([num2str(sum(params.L23_cells)) ' L23 cells']);
disp([num2str(sum(params.L4_cells)) ' L4 cells']);
disp([num2str(sum(params.L56_cells)) ' L56 cells']);
    
%% separate by narrow and broad spiking 

params.Pyr = spikeWaveform_separator(spikes);
disp('separated by waveform');

%% save to analysis path 

save([analysis_path,rec_date,'_',spikeSortFile,'_metaData.mat'],'params');