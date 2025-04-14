function [p_value] = ll_test(full_model,reduced_model)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ll_full = full_model.LogLikelihood;
ll_reduced = reduced_model.LogLikelihood;

% Calculate test statistic
lr_stat = -2 * (ll_reduced - ll_full);

% Degrees of freedom = difference in number of parameters
df = 1; % One variance component for random effect

% Calculate p-value (note: testing on boundary, so divide by 2)
p_value = 0.5 * (1 - chi2cdf(lr_stat, df));
% Display results
fprintf('Likelihood Ratio Test for Random Effect of mouse ID:\n');
fprintf('LR statistic: %.4f\n', lr_stat);
fprintf('Degrees of freedom: %d\n', df);
fprintf('p-value: %.4f\n', p_value);

end