function [outStruct] = trialDropper(inputStruct,trials,action)
    % trialDropper: modify specified trials in an experiment info
    % struct
    % 
    % Inputs:
    %   inputStruct - scanbox .mat file
    %   trials      - trials to be removed from the struct. e.g. [20 40 91
    %   92]
    %   action      - provided as an integer:
    %                 1: deletes the trials (changes total length)
    %                 2: fills trials with NaN (preserves total length)
    %                 DEFAULT is 1 if there is no given value.
    %
    % Output:
    %   outStruct   - struct with all specified trials removed from all
    %   fields


% load('G:\Behavior\Data\data-i2197-250715-1622.mat')
% load('G:\home\ACh\Data\2p_data\i2197\250715\003\003_000_000.mat')
% inputStruct = input(1);

[sz1 sz2] = size(inputStruct);
if sz1 ~= 1 | sz2 ~= 1
    error('Input struct has unexpected size.')
end

if nargin < 3
    action = 1;
end

nTrials = length(inputStruct.tGratingContrast);
fdnames = fieldnames(inputStruct)';

for i=1:length(fdnames)
    this_fd = inputStruct.(fdnames{i});
    if length(this_fd)==nTrials && ~strcmp(fdnames(i), 'savedEvents')
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

end
