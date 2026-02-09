function plotSizeResponseByExperiment(data1, data2, cell_indices1, cell_indices2, contrasts, sizes, exp_idx, nKeep_concat, varargin)
% PLOTSIZERESPONSEBYEXPERIMENT Plot size response curves for each experiment separately
% Creates one figure per experiment
% 
% Usage:
%   plotSizeResponseByExperiment(data1, data2, cell_indices1, cell_indices2, contrasts, sizes, exp_idx, nKeep_concat) 
%   plotSizeResponseByExperiment(..., 'Name', Value, ...)
%
% All name-value pairs are passed directly to plotSizeResponse

nSess = length(nKeep_concat);

for iSess = 1:nSess
    sess_cells1 = cell_indices1(ismember(cell_indices1, find(exp_idx == iSess)));
    sess_cells2 = cell_indices2(ismember(cell_indices2, find(exp_idx == iSess)));
    
    if isempty(sess_cells1) && isempty(sess_cells2)
        continue;
    end
    
    plotSizeResponse(data1, data2, sess_cells1, sess_cells2, contrasts, sizes, varargin{:});
    sgtitle(sprintf('Experiment %d', iSess));
end

end