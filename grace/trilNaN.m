function [X,Y] = trilNaN(inputArg1,inputArg2)
%Tril but NaNs instead of 0s
%   Detailed explanation goes here
Y = ones(size(X));
Y = tril(Y);
ind = find(Y);
X = X(ind);
end

