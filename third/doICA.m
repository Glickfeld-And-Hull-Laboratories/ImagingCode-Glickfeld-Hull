function [B, iternum] = doICA(X, nIC, ica_A_guess, termtol, maxrounds, sclass)
%from mukamel PCA-ICA code- was fpica_standardica         
numSamples = size(X,2);

B = ica_A_guess;
BOld = zeros(size(B));

iternum = 0;
minAbsCos = 0;

errvec = zeros(maxrounds,1, sclass);
while (iternum < maxrounds) && ((1 - minAbsCos)>termtol)
    iternum = iternum + 1;

    % Symmetric orthogonalization.
    B = (X * ((X' * B) .^ 2)) / numSamples;
    B = B * real(inv(B' * B)^(1/2));

    % Test for termination condition.
    minAbsCos = min(abs(diag(B' * BOld)));

    BOld = B;
    errvec(iternum) = (1 - minAbsCos);
end

if iternum<maxrounds
    fprintf('Convergence in %d rounds.\n', iternum)
else
    fprintf('Failed to converge; terminating after %d rounds, current change in estimate %3.3g.\n', ...
        iternum, 1-minAbsCos)
end
