function [opt_Az, opt_El] = findOptimalStimLocation(fit_true_vec, goodfit_ind)
% findOptimalStimLocation - find the Az/El point that minimizes the sum of
% squared distances to all fitted RF centers
%
% Usage:
%   [opt_Az, opt_El] = findOptimalStimLocation(fit_true_vec, goodfit_ind)
%
% Inputs:
%   fit_true_vec  - [nCells x 10] matrix of fit parameters (col 4=El0, col 5=Az0)
%   goodfit_ind   - indices of cells passing quality filter
%
% Outputs:
%   opt_Az        - optimal azimuth (deg)
%   opt_El        - optimal elevation (deg)

Az_centers = fit_true_vec(goodfit_ind, 5);
El_centers = fit_true_vec(goodfit_ind, 4);

% minimize sum of absolute distances - solution is the median
opt_Az = median(Az_centers);
opt_El = median(El_centers);

fprintf('Optimal stimulus location: Az=%.2f deg, El=%.2f deg\n', opt_Az, opt_El)
fprintf('Mean distance from optimal: %.2f deg\n', ...
    mean(sqrt((Az_centers-opt_Az).^2 + (El_centers-opt_El).^2)))