function [idxOut,uniqueOut] = get_idx_from_input(trialList,varargin)
p = inputParser;
p.addParamValue('min_trials', 10, @isnumeric);

% parse inputs
p.parse(varargin{:});
params = p.Results;

uniqueOut = unique(trialList);
idxOut = arrayfun(@(x) trialList==x,uniqueOut,'un',0);
keep_i = cell2mat(cellfun(@(x) sum(x) >= params.min_trials ,idxOut,'un',0));
uniqueOut = uniqueOut(keep_i);
idxOut = idxOut(keep_i);