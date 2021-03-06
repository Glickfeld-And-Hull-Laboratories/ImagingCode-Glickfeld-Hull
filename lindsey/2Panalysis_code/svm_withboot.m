function [pct_corr,svmStruct] = svm_withboot(data, groups, nboot)
% svm applied to binary data- train and test on holdout data (50%)
% data: n x p matrix where n is trials/observation and p are predictors
% groups: n x 1 matrix to identify trial category
% nboot: number of repeats of train/test outcome

pct_corr = nan(1,nboot);
for i = 1:nboot
    [train, test] = crossvalind('holdOut',groups);
%     cp = classperf(groups);
    svmStruct = fitcsvm(data(train,:),groups(train));
    out = predict(svmStruct, data(test,:));
    pct_corr(1,i) = 1-(sum(abs(groups(test,:)-out))./sum(test));
end

