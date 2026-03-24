function [distMap, dist_vec] = plotRFdistanceMap(fit_true_vec, goodfit_ind, mask_cell, ref_Az, ref_El, varargin)
% plotRFdistanceMap - compute and visualize linear distance of each cell's
% fitted RF center from a reference Az/El point
%
% Usage:
%   distMap = plotRFdistanceMap(fit_true_vec, goodfit_ind, mask_cell, ref_Az, ref_El)
%   distMap = plotRFdistanceMap(..., 'savePath', fullfile(isilonName, ...))
%
% Inputs:
%   fit_true_vec  - [nCells x 10] matrix of fit parameters (col 4=El0, col 5=Az0)
%   goodfit_ind   - indices of cells passing quality filter
%   mask_cell     - [y x] image where pixel value = cell index
%   ref_Az        - reference azimuth (deg)
%   ref_El        - reference elevation (deg)
%
% Optional:
%   'savePath'    - full path string to save PDF
%
% Output:
%   distMap       - [y x] image of distances (NaN for non-goodfit cells)

savePath = [];
if nargin > 5
    for i = 1:2:length(varargin)
        if strcmp(varargin{i}, 'savePath')
            savePath = varargin{i+1};
        end
    end
end

distMap = NaN(size(mask_cell));
dist_vec = NaN(max(goodfit_ind), 1);

for i = 1:length(goodfit_ind)
    c = goodfit_ind(i);
    dAz = fit_true_vec(c, 5) - ref_Az;
    dEl = fit_true_vec(c, 4) - ref_El;
    d   = sqrt(dAz^2 + dEl^2);
    dist_vec(c) = d;
    ind = find(mask_cell == c);
    distMap(ind) = d;
end

imAlpha = ones(size(distMap));
imAlpha(isnan(distMap)) = 0;

figure; clf;
subplot(1,2,1)
imagesc(distMap, 'AlphaData', imAlpha)
axis image off
colormap(jet)
h = colorbar;
ylabel(h, 'Distance (deg)', 'Rotation', 270.0, 'VerticalAlignment', 'bottom')
title(sprintf('RF distance from Az=%.1f, El=%.1f', ref_Az, ref_El))
set(gca, 'color', [0 0 0], 'TickDir', 'out')

subplot(1,2,2)
hist(dist_vec(~isnan(dist_vec)))
xlabel('Distance (deg)')
ylabel('# cells')
title(sprintf('n=%d', sum(~isnan(dist_vec))))
set(gca, 'TickDir', 'out', 'box', 'off')

set(gcf, 'Position', [100 300 1000 400])

if ~isempty(savePath)
    print(savePath, '-dpdf')
end