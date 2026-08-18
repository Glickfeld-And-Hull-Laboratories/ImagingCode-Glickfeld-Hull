function mouse_name = findMouseFromIndex(cell_index, nKeep_concat, mice)
% FINDMOUSEFROMINDEX Find which mouse a given cell index came from
%   mouse_name = findMouseFromIndex(cell_index, nKeep_concat, mice)
%   
%   Inputs:
%   cell_index - the index of the cell (1 to sum(nKeep_concat))
%   nKeep_concat - vector containing number of cells from each mouse
%   mice - cell array or string array of mouse ID names
%   
%   Output:
%   mouse_name - name/ID of the mouse the cell came from

% Calculate cumulative sum to get boundaries
cumulative_cells = cumsum(nKeep_concat);

% Find which mouse this index belongs to
mouse_idx = find(cell_index <= cumulative_cells, 1, 'first');

% Return the mouse name
mouse_name = mice(mouse_idx,:);

end

% Example usage:
% mouse_name = findMouseFromIndex(50, nKeep_concat, mice);
% fprintf('Cell index 50 came from mouse %s\n', mouse_name);

% For multiple indices at once:
function mouse_names = findMouseFromIndices(cell_indices, nKeep_concat, mice)
% FINDMOUSEFROMINDICES Find which mouse multiple cell indices came from
%   mouse_names = findMouseFromIndices(cell_indices, nKeep_concat, mice)

% Calculate cumulative sum to get boundaries
cumulative_cells = cumsum(nKeep_concat);

% Initialize output
mouse_names = cell(size(cell_indices));

% Find which mouse each index belongs to
for i = 1:length(cell_indices)
    mouse_idx = find(cell_indices(i) <= cumulative_cells, 1, 'first');
    mouse_names{i} = mice{mouse_idx};
end

end