function [mdl pct_corr] = fitrsvm_withcrossval(X,Y);
%X- trials x cells 
%Y- trials x 1
mdl = fitrsvm(X,Y);
Y_hat = predict(mdl,X);
Y_dig = zeros(size(Y));
Y_dig(find(Y>0.5)) = 1;
all_crit = zeros(1,100);
crit_list = min(Y_hat):(max(Y_hat)-min(Y_hat))./99:max(Y_hat);
for i = 1:100;
    all_crit(i) = mean((Y_hat>crit_list(i)) == Y_dig);
end
[max_val crit] = max(all_crit(:));
out = zeros(size(Y,1),1);
for i = 1:size(Y,1)
    X_temp = X;
    Y_temp = Y;
    X_temp(i,:) = [];
    Y_temp(i) = [];
    mdl_temp = fitrsvm(X_temp,Y_temp);
    pred_temp = predict(mdl_temp,X(i,:));
    out(i,1) = pred_temp>crit_list(crit);
end
pct_corr = mean(out==Y_dig);
    
