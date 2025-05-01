function acf = autocorr_TH(y, numLags, plotFlag)
    % autocorr_TH: Computes the autocorrelation function
    % 
    % Inputs:
    %   y        - time series vector
    %   numLags  - number of lags to compute
    %   plotFlag - (optional) if true, plots the ACF
    %
    % Output:
    %   acf - autocorrelation values from lag 0 to numLags

    if nargin < 3
        plotFlag = true;
    end

    y = y(:);  % ensure column vector
    y = y - mean(y);  % demean the series
    N = length(y);
    
    % Full autocorrelation using unbiased estimator
    [c, lags] = xcorr(y, numLags, 'unbiased');

    % Normalize by variance to get correlation
    acf = c(lags >= 0) / c(numLags + 1);  % c(numLags+1) is lag 0

    if plotFlag
        stem(0:numLags, acf, 'filled');
        xlabel('Lag');
        ylabel('Autocorrelation');
        title('Autocorrelation Function (ACF)');
        grid on;
        hold on;
        % 95% confidence bounds (approximate)
        confBound = 1.96 / sqrt(N);
        yline(confBound, 'r--');
        yline(-confBound, 'r--');
        hold off;
    end
end