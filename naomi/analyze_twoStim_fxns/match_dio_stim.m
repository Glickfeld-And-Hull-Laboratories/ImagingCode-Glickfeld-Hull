function [crop_photodio,crop_input] = match_dio_stim(photodio,input,keep_fields,varargin)

p = inputParser;
% p.addParamValue('crop_dio', false, @islogical);
% p.addParamValue('crop_input', true, @islogical);
p.addParamValue('dio_trials', 0, @isnumeric);
p.addParamValue('input_trials',0,@isnumeric);
p.addParamValue('verbose',true,@islogical);

% parse inputs
p.parse(varargin{:});
params = p.Results;

if size(photodio,2) > size(photodio,1)
        photodio = photodio';
end
if params.dio_trials == 0 
    params.dio_trials = 1:size(photodio,1);
end
if params.input_trials == 0
    params.input_trials = 1:input.trialSinceReset;
end

crop_photodio = photodio(params.dio_trials,:);

for field_i = 1:numel(keep_fields)
    if params.verbose
    disp(['cropping ' keep_fields{field_i}]);
    end
    tempCrop = input.(keep_fields{field_i})(params.input_trials);
    if iscell(tempCrop)
        if numel(tempCrop{1})==1
            tempCrop = cell2mat(cellfun(@(x) double(x),tempCrop,'un',0));
        end
    end
    crop_input.(keep_fields{field_i}) = tempCrop;
end

% check to make sure same # of trials
if numel(params.dio_trials)~=numel(params.input_trials)
    disp('something wrong? different # of trials for stim on/photodio');
end
crop_input.nTrials = numel(params.dio_trials);

