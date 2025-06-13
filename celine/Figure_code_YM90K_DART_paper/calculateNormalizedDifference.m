function [norm_diff, bsln_std] = calculateNormalizedDifference(stat_data, loc_data, pre, post, nCon, nCells, varargin)
% CALCULATENORMALIZEDDIFFERENCE Calculate normalized differences between pre and post responses
%
% Usage:
%   [norm_diff, bsln_std] = calculateNormalizedDifference(stat_data, 
%                          loc_data, pre, post, nCon, nCells)
%   [norm_diff, bsln_std] = calculateNormalizedDifference(stat_data, 
%                          loc_data, pre, post, nCon, nCells, nSizes)
%
% Inputs:
%   stat_data - Cell array of stationary trial responses {contrast, day}{cell}
%               or {contrast, size, day}{cell} if size dimension present
%   loc_data - Cell array of locomotion trial responses {contrast, day}{cell}
%              or {contrast, size, day}{cell} if size dimension present
%   pre - Index for pre-treatment day
%   post - Index for post-treatment day
%   nCon - Number of contrast conditions
%   nCells - Total number of cells
%   nSizes - (Optional) Number of size conditions. If not provided, assumes no size dimension.
%
% Outputs:
%   norm_diff - Normalized difference matrix (2, nCon, nCells) or (2, nCon, nSizes, nCells)
%               First index: 1=stationary, 2=locomotion
%   bsln_std - Baseline standard deviation matrix with same dimensions as norm_diff
%
% This function calculates the normalized change in response for each cell
% between pre and post conditions, normalized by the baseline variability.
% It can handle data with or without a size dimension.

% Check for optional nSizes parameter
if nargin > 6
    nSizes = varargin{1};
    has_size_dim = true;
else
    nSizes = 1;
    has_size_dim = false;
end

% Initialize output arrays with appropriate dimensions
if has_size_dim
    norm_diff = nan(2, nCon, nSizes, nCells);
    bsln_std = nan(2, nCon, nSizes, nCells);
else
    norm_diff = nan(2, nCon, nCells);
    bsln_std = nan(2, nCon, nCells);
end

% Calculate for each cell and stimulus condition
for i = 1:nCells
    for iCon = 1:nCon
        for iSize = 1:nSizes
            % Handle indexing based on whether there's a size dimension
            if has_size_dim
                stat_idx = {iCon, iSize, pre};
                stat_idx_post = {iCon, iSize, post};
                loc_idx = {iCon, iSize, pre};
                loc_idx_post = {iCon, iSize, post};
            else
                stat_idx = {iCon, pre};
                stat_idx_post = {iCon, post};
                loc_idx = {iCon, pre};
                loc_idx_post = {iCon, post};
            end
            
            % For stationary trials
            mean_pre_stat = mean(stat_data{stat_idx{:}}{i}, 'omitnan');
            mean_post_stat = mean(stat_data{stat_idx_post{:}}{i}, 'omitnan');
            std_pre_stat = std(stat_data{stat_idx{:}}{i}, 'omitnan');
            
            % Calculate normalized difference for stationary condition
            if std_pre_stat > 0
                norm_diff_stat = (mean_post_stat - mean_pre_stat) / std_pre_stat;
            else
                norm_diff_stat = NaN;
            end

            % For running trials
            mean_pre_loc = mean(loc_data{loc_idx{:}}{i}, 'omitnan');
            mean_post_loc = mean(loc_data{loc_idx_post{:}}{i}, 'omitnan');
            std_pre_loc = std(loc_data{loc_idx{:}}{i}, 'omitnan');
            
            % Calculate normalized difference for locomotion condition
            if std_pre_loc > 0
                norm_diff_loc = (mean_post_loc - mean_pre_loc) / std_pre_loc;
            else
                norm_diff_loc = NaN;
            end

            % Store results in output matrices with appropriate indexing
            if has_size_dim
                norm_diff(1, iCon, iSize, i) = norm_diff_stat; % 1 = stationary
                norm_diff(2, iCon, iSize, i) = norm_diff_loc;  % 2 = running
                bsln_std(1, iCon, iSize, i) = std_pre_stat;    % Std for stationary trials
                bsln_std(2, iCon, iSize, i) = std_pre_loc;     % Std for running trials
            else
                norm_diff(1, iCon, i) = norm_diff_stat; % 1 = stationary
                norm_diff(2, iCon, i) = norm_diff_loc;  % 2 = running
                bsln_std(1, iCon, i) = std_pre_stat;    % Std for stationary trials
                bsln_std(2, iCon, i) = std_pre_loc;     % Std for running trials
            end
        end
    end 
end

% Replace infinity values with NaN
norm_diff(norm_diff == -Inf) = NaN;
norm_diff(norm_diff == Inf) = NaN;

bsln_std(bsln_std == -Inf) = NaN;
bsln_std(bsln_std == Inf) = NaN;

end