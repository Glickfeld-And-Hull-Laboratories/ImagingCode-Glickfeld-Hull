% function analyze_twoStim(datename) 

addpath 'Z:\home\Naomi\Electrophysiology\New Analysis\analyze_twoStim_fxns\'
%% load data
Rig=input('Neuropixel(1) or Blackrock(2)?')

if Rig==1
save_path = ['Z:\home\Naomi\Electrophysiology\Neuropixel\',datename,'\'];
tempPath = ['Z:\home\Naomi\Electrophysiology\Neuropixel\',datename,'\'];
else 
save_path = ['Z:\home\Naomi\Electrophysiology\BlackRock\',datename,'\'];
tempPath = ['Z:\home\Naomi\Electrophysiology\BlackRock\',datename,'\'];
end
% tempPath = ['Z:\All_staff\home\jen\Analysis\Spyking_Circus\',datename,'\'];
% tempPath = ['Z:\home\jen\Analysis\Spyking_Circus\',datename,'\'];
paramsFile = dir(fullfile(tempPath,'*metaData.mat'));
data.params_path = [paramsFile.folder, '\', paramsFile.name];

load(data.params_path);
disp(params.info.expt_list');
exptIdx = input('which expt?');

exptname = params.info.expt_list{exptIdx};
data.info.expt = exptname;

load(params.info.spikesFile);
load(params.info.triggersFile);
load(params.info.layersFile);
 
disp('loaded');

%% set analysis params

data.info.thresholdFactor = 0.3;
data.info.sampleFreq = 30000;
data.info.bin_sizeS = 0.01;
data.info.preStim_timeS = 0.1;
data.info.postStim_timeS = 0.5;

data.info.maxWindowLHS = 0.05;
data.info.maxWindowRHS = 0.18;
data.info.nCells = numel(spikes.spiketime_S);
data.info.maxWindowBins= 2;
data.info.gaussSmoothWin = 0.06;
data.info.nBaselineTrials = 30;
data.info.nLaserOnTrials = 30;

%% find stim on and extract info from input file 

keep_fields = {'tISITimeMs','tBlock2TrialNumber','wheelSpeedValues','tStimOneGratingContrast','wheelSpeedTimesUs','tThisTrialStartTimeMs','tItiTimeMs'};

expt_idx = find(params.info.([exptname,'Idx']));

%%
count=1
for file_i = expt_idx
load(params.info.([exptname,'InputFile']){count});

if iscell(Trigger.data)
    % trigger data is from multi-experiment day
    shift_idx = sum(cell2mat(cellfun(@(x) size(x,2),Trigger.data(1:file_i-1),'un',0)));
else 
    % trigger data is from the only experiment that day
    Trigger.data = {Trigger.data};
    shift_idx = sum(cell2mat(cellfun(@(x) size(x,2),Trigger.data(1:file_i-1),'un',0)));
end


% find stim on times for stimulus and laser
[temp_stimOn,~] = photodioOn_JL(Trigger.data{file_i},data.info.sampleFreq,data.info.thresholdFactor,Trigger.photoID);

% reshape + shift
try
temp_stimOn = [temp_stimOn(1:2:end)' temp_stimOn(2:2:end)'];
catch
temp_stimOn = [temp_stimOn(1:2:end-1)' temp_stimOn(2:2:end)'];
end

temp_stimOn = temp_stimOn + shift_idx;

figure; plot(Trigger.data{file_i}(Trigger.photoID,1:2000000)); vline(temp_stimOn(1,:)-shift_idx); 


if inputs.trialSinceReset > size(temp_stimOn,1)
    disp('found more stimulus info than stimulus on (typical) - cropping input file...')
    [crop_dio{count},crop_input(count)] = match_dio_stim(temp_stimOn,inputs,keep_fields,'input_trials',1:size(temp_stimOn,1));
elseif inputs.trialSinceReset < size(temp_stimOn,1)
    disp('found more stimulus on than stimulus info - cropping dio on...')
    [crop_dio{count},crop_input(count)] = match_dio_stim(temp_stimOn,inputs,keep_fields,'dio_trials',1:input.trialSinceReset);    
elseif size(temp_stimOn,1) == size(temp_stimOn,1)
    disp('found same number of stimulus on and stimulus info - no cropping')
    [crop_dio{count},crop_input(count)] = match_dio_stim(temp_stimOn,inputs,keep_fields);  
end

count = count+1;
end

for field = keep_fields
    fname = field{1};
expt_data.(fname) = [crop_input.(fname)];
end

stimOn = cell2mat(crop_dio');
data.info.nTrials = size(stimOn,1);


disp([num2str(data.info.nTrials) ' stim on total']);
[~,all_contrasts] = get_idx_from_input(expt_data.tStimOneGratingContrast,'min_trials',10);
data.info.all_contrasts = all_contrasts;


%% find moving trials 
beforeStimOnWheel = 400; % time window for before stim on
afterStimOnWheel = 100; % time window for after stim on

wheelSI = cellfun(@(x) floor(mean(diff(x))/1000),expt_data.wheelSpeedTimesUs,'un',0);
temp_LHS = cellfun(@(y,z,s,start) floor(find((y/1000-z)>start,1,'first') - beforeStimOnWheel/s),expt_data.wheelSpeedTimesUs,num2cell(expt_data.tThisTrialStartTimeMs),wheelSI,num2cell(expt_data.tItiTimeMs),'un',0);
stimWindow_wheelPts = cellfun(@(x,s) x:(x + (beforeStimOnWheel+afterStimOnWheel)/s),temp_LHS,wheelSI,'un',0);
try
    isMovingT = cell2mat(cellfun(@(x,y,z) mean(abs(x(y)))>=(500/z),expt_data.wheelSpeedValues,stimWindow_wheelPts,wheelSI,'un',0));
catch
    isMovingT = cell2mat(cellfun(@(x,y,z) mean(abs(x(y(1):end)))>=(500/z),expt_data.wheelSpeedValues,stimWindow_wheelPts,wheelSI,'un',0));
end
disp([num2str(sum(isMovingT)) ' moving trials']);


%% get broad spont FR over time 

% % for spontaneous activity
[raster_spontHz,prespiketime] = cellfun(@(x) get_raster_NB(x,(stimOn(:,1)/data.info.sampleFreq)-0.5,'stop',0.5,'count',true),spikes.spiketime_S,'un',0);
% raster_spontHz is the count of spike
% preSpikeTimes is the actual timestamp
%below convert to Hz
data.spontHz_all = cell2mat(raster_spontHz)/0.5;
% 

% % for stim evoked activity 
[raster_visHz,SDspike] = cellfun(@(x) get_raster_NB(x,(stimOn(:,1)/data.info.sampleFreq)+0.03,'stop',0.3,'count',true),spikes.spiketime_S,'un',0);
% 
%this (these) need to match the time window sampled above.
data.visHz_all = cell2mat(raster_visHz)/0.3;



%% find visually responsive cells 

clear temp_ttest
raster_visStim = cellfun(@(x) get_raster_NB(x,(stimOn(:,1)/data.info.sampleFreq),'start',-0.1,'stop',0.1),spikes.spiketime_S,'un',0);
preStimMean = cellfun(@(y) cell2mat(cellfun(@(x) sum(x<0),y,'un',0)), raster_visStim,'un',0);
periStimMean = cellfun(@(y) cell2mat(cellfun(@(x) sum(x>0),y,'un',0)), raster_visStim,'un',0); 
temp_ttest = cell2mat(cellfun(@(x,y) ttest2(x,y,'tail','left'),preStimMean,periStimMean,'un',0));
temp_ttest(isnan(temp_ttest)) = 0;
params.visResp = logical(temp_ttest);
disp([num2str(sum(params.visResp)) ' vis resp']);


%% get psth

% for stim_i = [1 2]
%     [temp_psth,~] = cellfun(@(x) get_psth(x,(stimOn(:,stim_i)/data.info.sampleFreq),'start',.preStim_timeS,'bin_size',data.info.bin_sizeS,'stop',data.info.postStim_timeS,'smoothWin',gausswin),spikes.spiketime_S,'un',0);
%     data.stimOn_psth{stim_i} = cell2mat(temp_psth);
% end

