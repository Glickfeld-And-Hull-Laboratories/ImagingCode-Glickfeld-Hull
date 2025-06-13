function [t_stat, crit_vals, result] = adftest_TH(x, p, includeTrend)
% ADFTEST_TH Performs a manual Augmented Dickey-Fuller test
%
% Inputs:
%   x            - Time series vector
%   p            - Number of lagged Δx terms
%   includeTrend - Boolean: include deterministic time trend if true
%
% Outputs:
%   t_stat    - t-statistic for γ (lagged level)
%   crit_vals - Struct of critical values (crit_1, crit_5, crit_10)
%   result    - String stating if null hypothesis is rejected

    x = x(:); % ensure column vector
    n = length(x);

    dx = diff(x);                       % Δx
    y = dx((p+1):end);                  % dependent variable

    x_lag1 = x((p+1):(end-1));          % lagged level x_{t-1}
    X = x_lag1;                         % initialize regressor matrix

    % Add lagged differences
    for i = 1:p
        lagged_dx = dx((p+1-i):(end-i));
        X = [X lagged_dx];
    end

    % Add intercept
    X = [ones(size(X,1),1) X];

    % Add deterministic trend if requested
    if includeTrend
        trend = (p+1):(n-1);  % Time trend
        X = [X trend'];
    end

    % OLS regression
    b = (X' * X) \ (X' * y);
    resid = y - X * b;
    s2 = sum(resid.^2) / (length(y) - size(X,2));
    covb = s2 * inv(X' * X);

    % Determine position of γ in regressor matrix:
    % b = [intercept, γ, δ1, δ2, ..., trend?]
    gamma_pos = 2;

    t_stat = b(gamma_pos) / sqrt(covb(gamma_pos, gamma_pos));

    % Valid MATLAB struct field names
    if includeTrend
        crit_vals = struct('crit_1', -3.96, 'crit_5', -3.41, 'crit_10', -3.13);
    else
        crit_vals = struct('crit_1', -3.43, 'crit_5', -2.86, 'crit_10', -2.57);
    end

    % Test decision
    if t_stat < crit_vals.crit_5
        result = 'Reject null: time series is stationary.';
    else
        result = 'Fail to reject null: time series is non-stationary.';
    end

    % Display result
    fprintf('ADF test statistic: %.4f\n', t_stat);
    fprintf('Critical values: 1%% = %.2f, 5%% = %.2f, 10%% = %.2f\n', ...
        crit_vals.crit_1, crit_vals.crit_5, crit_vals.crit_10);
    fprintf('%s\n', result);
end