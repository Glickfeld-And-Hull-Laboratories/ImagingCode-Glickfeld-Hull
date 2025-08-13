function [h_match, p_match, resp_match] = findResponsiveCells(data_dfof_trial_match, ...
    tCon_match, tDir_match, tSize_match, stimStart, stimEnd, nTrials, ...
    nCells, nDir, nCon, nSize, dirs, cons, sizes)
%FINDRESPONSIVECELLS Identify significantly responsive cells across days
%
% Inputs:
%   data_dfof_trial_match - cell array of dfof data for each day 
%   tCon_match - cell array of contrast trial indices for each day
%   tDir_match - cell array of direction trial indices for each day  
%   tSize_match - cell array of size trial indices for each day
%   stimStart - frame number when stimulus starts
%   stimEnd - frame number when stimulus ends
%   nTrials - array of trial counts for each day
%   nCells - number of cells
%   nDir, nCon, nSize - number of directions, contrasts, sizes
%   dirs, cons, sizes - arrays of stimulus parameter values
%
% Outputs:
%   h_match - cell array of significance test results (1=significant, 0=not)
%   p_match - cell array of p-values from t-tests
%   resp_match - cell array of logical arrays indicating responsive cells

nd = length(data_dfof_trial_match); % number of days

resp_win = (stimStart+3):(stimEnd+3); %at 15 hz, 3 frames = ~200 ms.
base_win = 1: stimStart-1;

% Initialize output cell arrays
h_match = cell(1,nd);
p_match = cell(1,nd);
resp_match = cell(1,nd);

% for each day
for id = 1:nd
    h = zeros(nCells, nDir, nCon, nSize);
    p = zeros(nCells, nDir, nCon, nSize);
    
    tCon = tCon_match{id}(:,1:nTrials(id));
    tSize = tSize_match{id}(:,1:nTrials(id));  
    tDir = tDir_match{id}(:,1:nTrials(id));
    data_dfof_trial = data_dfof_trial_match{id}; 

    for iDir = 1:nDir
        ind_dir = find(tDir == dirs(iDir));
        for iCon = 1:nCon
            ind_con = find(tCon == cons(iCon));
            for iSize = 1:nSize                
                ind_size = find(tSize == sizes(iSize));
                ind = intersect(intersect(ind_dir, ind_con), ind_size);
                
                % Only calculate the t-test for responsiveness
                [h(:,iDir,iCon,iSize), p(:,iDir,iCon,iSize)] = ttest(...
                    nanmean(data_dfof_trial(resp_win,ind,:),1), ...
                    nanmean(data_dfof_trial(base_win,ind,:),1), ...
                    'dim',2,'tail','right','alpha',0.01./(nDir*nCon*nSize-1));
            end
        end
    end 
    
    % Determine which cells are responsive to at least one condition
    h_pass = sum(sum(sum(h(:,:,:,:),2),3),4);
    resp = logical(h_pass);
    
    % Store results
    h_match{id} = h;
    p_match{id} = p;
    resp_match{id} = resp;
end

end