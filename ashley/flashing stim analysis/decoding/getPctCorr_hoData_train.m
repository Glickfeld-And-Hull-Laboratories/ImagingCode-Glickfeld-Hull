function [pctCorr yhat_temp, targetGLM] = getPctCorr_hoData_train(X,Y,ind_sub,ind_train,dv)
%X- responses for all trials all cells
%Y- category for all trials
%ind_sub- set of trials to test
%ind_train- set of trials to train;
nt = length(ind_sub);
yhat = nan(nt,1);
C = eye(size(X,2));
[~,~,targetGLM] = glmfit(X(ind_train,:)*C,Y(ind_train),'binomial');
for i = 1:nt
    ind_test = ind_sub(i);
    X_holdout = X(ind_test,:);
    if find(ind_train==ind_test)
        othersInd = setdiff(ind_train,ind_test);

        [~,~,othersGLM] = glmfit(X(othersInd,:),Y(othersInd),'binomial');
        yhat_temp(i,:) = glmval(othersGLM.beta,X_holdout,'logit');
        yhat(i,:) = glmval(othersGLM.beta,X_holdout,'logit') > dv;
    else
        yhat_temp(i,:) = glmval(targetGLM.beta,X_holdout,'logit'); 
        yhat(i,:) = glmval(targetGLM.beta,X_holdout,'logit') > dv;
    end
end
pctCorr = mean((Y(ind_sub)-yhat) == 0);

end