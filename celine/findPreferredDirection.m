function [pref_dir, pref_response] = findPreferredDirection(response_data, cell_indices)
%FINDPREFERREDDIRECTION Find the preferred direction for a set of cells
%
% Inputs:
%   response_data - 4D array [nCells x nDir x nCon x nSize] of response values
%                   (should be significant responses only, i.e., resp_sig)
%   cell_indices - (optional) vector of cell indices to analyze
%                  if empty or not provided, analyzes all cells
%
% Outputs:
%   pref_dir - vector of preferred direction indices for each analyzed cell
%   pref_response - vector of maximum response values at preferred direction
%
% Method: For each cell, averages responses across contrasts and sizes,
%         then finds the direction with maximum response

[nCells, nDir, nCon, nSize] = size(response_data);

% Handle optional cell_indices input
if nargin < 2 || isempty(cell_indices)
    cell_indices = 1:nCells;
end

nAnalyzedCells = length(cell_indices);
pref_dir = zeros(1, nAnalyzedCells);
pref_response = zeros(1, nAnalyzedCells);

% Special case: only one direction tested
if nDir == 1
    pref_dir(:) = 1;
    % Calculate mean response across contrasts and sizes for the single direction
    for i = 1:nAnalyzedCells
        iCell = cell_indices(i);
        pref_response(i) = mean(mean(squeeze(response_data(iCell, 1, :, :))));
    end
    return;
end

% Find preferred direction for each cell
for i = 1:nAnalyzedCells
    iCell = cell_indices(i);
    
    % Get responses for this cell: [nDir x nCon x nSize]
    cell_responses = squeeze(response_data(iCell, :, :, :));
    
    % Average across contrasts (dim 2) and sizes (dim 3)
    dir_means = mean(mean(cell_responses, 2), 3);
    
    % Find direction with maximum response
    [max_val, pref_dir(i)] = max(dir_means);
    pref_response(i) = max_val;
end

end