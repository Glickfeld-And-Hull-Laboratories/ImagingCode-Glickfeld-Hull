function result = SampleEn_TH(data, r, dim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes Sample Entropy of a time series                %
% The Sample Entropy is the negative average natural logarithm of the   %
% conditional probability that two sequences that are similar for m     %
% points remain similar within a tolerance r at the next point, with    %
% self-matches not included in the probability calculation. It is a     % 
% time measure of complexity.                                           %
%                                                                       %
% Richman, J. S., & Moorman, J. R. (2000). Physiological time-series    %
%         analysis using approximate entropy and sample entropy. Am J   %
%         Physiol Heart Circ Physiol, 278, H2039-H2049.                 %
%         doi:10.1152/ajpheart.2000.278.6.H2039                         %
% Keller, K., Unakafov, A., & Unakafova, V. (2014). Ordinal Patterns,   %
%         Entropy, and EEG. Entropy, 16(12), 6212-6239.                 %
%         doi:10.3390/e16126212                                         %
%                                                                       %
%                                                                       %
% INPUTS                                                                %
%  Parameters:                                                          %
%   data = time series                                                  %                                               %
%   r = tolerance => percent of std to accept as a match                %
%   dim = embedded dimension => length of comparison window             %
%                                                                       %
%   Standard values are r=0.2, dim=3                                    %
%                                                                       %
% OUTPUTS                                                               %
%   result = sample entropy value                                       %
%                                                                       %
%                       Created by Brian Lord                           %
%                       University of Arizona                           %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sample Entropy calculation
% create matrix of data by embedding dimension
data_matrix = NaN(dim+1,length(data));
for a = 1:dim+1
    data_matrix(a,1:length(data)+1-a)=data(a:end);
end
data_matrix = data_matrix';
% take pairwise distances across each window
matrix_B = pdist(data_matrix(:,1:dim), 'chebychev');
matrix_A = pdist(data_matrix(:,1:dim+1), 'chebychev');
% take sum of counts of distances that fall within tolerance range
B = sum(matrix_B <= r*std(data));
A = sum(matrix_A <= r*std(data));
% calculate ratio of natural logarithm, with normalization correction
result = -log((A/B))*((length(data)-dim+1)/(length(data)-dim-1));
% correct inf to a maximum value
if isinf(result)
    result = -log(2 / ((length(data)-dim-1)*(length(data)-dim)));
end
