function [avg_resp_dir_shuf order] = doPlaidRespShuffle(avg_resp_dir,order)
    % creates new matrix with shuffled plaid responses across cells
    % order can be used to redo same shuffle
    
    nCells = size(avg_resp_dir,1);
    if nargin<2
        order = zeros(1,nCells);
        used_list = [];
        for i = 1:nCells
            use_list = setdiff(1:nCells,[i used_list]);
            use = randi(length(use_list),1);
            order(i) = use_list(use);
            used_list = [used_list use_list(use)];
        end
    end
    
    avg_resp_dir_shuf = avg_resp_dir;
    
    for i = 1:nCells
        avg_resp_dir_shuf(i,:,:,2,:) = avg_resp_dir_shuf(order(i),:,:,2,:);
    end
end