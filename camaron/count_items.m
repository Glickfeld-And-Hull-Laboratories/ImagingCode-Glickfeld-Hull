function count_items(cell_array)
% count items in each field of cell array and display total
    
    for i = 1:length(cell_array)
        count(i) = length(cell_array{i});
    end
    
    total = sum(count);
    disp(total)
end

