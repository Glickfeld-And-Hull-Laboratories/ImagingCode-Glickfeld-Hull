function result_cell = isAdjacent(nums);
% results = isAdjacent(nums);
% nums should be a numerical vector
% check if entries are numerically adjacent to each other
% outputs a cell where each row is an adjacent chunk
% column 1: index of adjacent elements
% column 2: length of the chunk

nNums = length(nums);
if nNums < 2
    error('Vector has only 1 element.');
end

result_cell = cell(1,2);
count = 1;
renew = 0;
temp_nums = [];
for i= 1:nNums
    if renew == 1
        temp_nums = [];
    end
    if i == nNums
        continue
    end
    if abs(nums(i)-nums(i+1)) == 1
        temp_nums = [temp_nums nums(i) nums(i+1)];
        renew = 0;
    else
        result_cell{count,1} = unique(temp_nums);
        result_cell{count,2} = length(temp_nums);
        count = count+1;
        renew = 1;
    end
end


end