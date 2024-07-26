% function E2_make_params(rec_date,varargin)
%%% This code creates a couple of parameters potentially important for analysis.
%%% This includes reading and outputting the spike data from Kilosort
    %%% A similar thing is done for neuropixel data - this could (should?)
    %%% be added to this file?
%%% Also as output is a layer analysis. This is not set up for neuropixels
%%% yet, but should be. The output is layer_bound and give the location
%%% boundaries of layer 2/3 and 4
clearvars

rec_date=input('Date?')


%% set expt variables
p = inputParser;
p.addParameter('filename', [], @ischar);

% parse inputs
% p.parse(varargin{:});
inputParam = p.Results;

addpath('Z:\home\Naomi\Electrophysiology\New Analysis\npy-matlab-master\npy-matlab\')
addpath('Z:\home\Naomi\Electrophysiology\New Analysis\make_params_fxns\')

analysis_path=['Z:\home\Naomi\Electrophysiology\Kilosort_Analysis\',rec_date,'\'];
% if SS_type==1
% rawData_path = ['Z:\home\Naomi\Electrophysiology\Neuropixel\',rec_date,'\Raw Data\'];
% elseif SS_type==2
rawData_path =['Z:\home\Naomi\Electrophysiology\BlackRock\', rec_date,'\Raw Data\'];
% en%% load experiment info 

fs = dir(fullfile(analysis_path, '*exptInfo.mat'));

%% Loads params - experiment information from E1
    load(fullfile(analysis_path, fs(1).name));
%% Names the spike sort file
    spikeSortFile = params.info.spikeFiles;
%% Not sure below
    % autoFill = true;

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

exptType = params.info.expt_list'; 
params.info.spikesFile = [analysis_path,rec_date,'_',spikeSortFile,'_spikes.mat'];
params.info.triggersFile = [analysis_path,rec_date,'_',spikeSortFile,'_triggers.mat'];
params.info.layersFile = [analysis_path,rec_date,'_',spikeSortFile,'_layerBound.mat'];
params.info.inputsfile = [analysis_path,rec_date,'_',spikeSortFile,'_input.mat'];
params.info.recDate = rec_date;
params.info.experimentTypes = unique(exptType);

for exp_i = 1:numel(params.info.experimentTypes)
    sub_idx = cell2mat(cellfun(@(x) strcmp(x,params.info.experimentTypes{exp_i}),exptType,'un',0));
    params.info.([params.info.experimentTypes{exp_i},'Idx']) = sub_idx;
    params.info.([params.info.experimentTypes{exp_i},'InputFile']) = cellfun(@(x) [analysis_path,rec_date,'_',x,'_input.mat'],temp_files(sub_idx),'un',0);
end

% make spikes 
loadPhyOutput_NB


%% LFP / Layers Analysis
            LFPanalysis_NB([rawData_path,layersFile,'.ns4'],params.info.triggersFile); 
%%

load(params.info.layersFile); load(params.info.spikesFile);

% extract unit waveforms 
% if size(spikes.template_big,1) == 82 % phy output
    params.waveform = spikes.template_big(22:end-11,:)';
    params.pt_ms = spikes.wave_stats.pt_ms/30;
% else % spyking circus output 
%     params.waveform = spikes.template_big(12:end,:)';
%     params.pt_ms = spikes.wave_stats.pt_ms;
% end

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
% end
