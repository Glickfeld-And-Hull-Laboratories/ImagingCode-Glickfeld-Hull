function [outStruct] = trialDropper(inputStruct,trials,action)
    % trialDropper: modify specified trials in an experiment info struct
    % 
    % Inputs:
    %   inputStruct - mWorks .mat file
    %   trials      - trials to be removed from the struct. e.g. [20 40 91 92]
    %   action      - provided as an integer:
    %                 1: deletes the trials (changes total length)
    %                 2: fills trials with NaN (preserves total length)
    %                 DEFAULT is 1 if there is no given value.
    %
    % Output:
    %   outStruct   - struct with all specified trials removed from all fields


% load('G:\Behavior\Data\data-i2197-250715-1622.mat') %input
% load('G:\home\ACh\Data\2p_data\i2197\250715\003\003_000_000.mat') %info
% inputStruct = input(1);

[sz1 sz2] = size(inputStruct);
if sz1 ~= 1 | sz2 ~= 1
    error('Input struct has unexpected size.')
end

if nargin < 3
    action = 1;
end

assignedNTrials = inputStruct.stopAfterNTrials; % get manually assigned ntrials for this experiment
nTrials = length(inputStruct.tAnnulusGratingDiameterDeg);

% error checking

if nTrials ~= assignedNTrials % compare collected with assigned nTrials
    warning('Stored nTrials is inconsistent with assigned nTrials. Check for abberant events and modification history.\n');
    disp('CAUTION: Trial deletion indices should be relative to STORED, NOT ASSIGNED nTrials.\n')
    fprintf('Stored %d trials, assigned %d trials.\n',nTrials,assignedNTrials);
end

if isfield(inputStruct,'modifiedStruct')
    disp('Past modification of the struct detected.')
    if inputStruct.modifiedAction == 1
        warning('Current version had some original elements deleted. See "modifiedTrials".');
    elseif inputStruct.modifiedAction == 2
        warning('Current version had some original elements replaced with NaN. See "modifiedTrials".');
    end
end

% modify struct

fdnames = fieldnames(inputStruct)';
changedFields = cell(1,1);

for i=1:length(fdnames)
    this_fd = inputStruct.(fdnames{i});
    if length(this_fd)==nTrials && ~strcmp(fdnames(i), 'savedEvents')
        changedFields{i} = fdnames{i};
        if action == 1
            this_fd(trials) = [];
        elseif action == 2
            n = length(trials);
            nanVector = repmat({NaN},1,n);
            this_fd(trials) = nanVector;
        else
            error('Unspecified action input. Read function documentation.')
        end
        outStruct.(fdnames{i}) = this_fd;
    else
        outStruct.(fdnames{i}) = inputStruct.(fdnames{i});
    end    
end

% add modification history
changedFields = changedFields(~cellfun("isempty",changedFields));
outStruct.modifiedAction = action;
outStruct.modifiedFields = changedFields;
outStruct.modifiedStruct = 1;
outStruct.modifiedTime = char(datetime('now'));
outStruct.modifiedTrials = trials;

end
