%create response matrix for all trials
trial_dfof = squeeze(mean(data_dfof(resp_win,:,:),1));

%transform into pca space for dimensionality reduction
[u s v] = pca(trial_dfof);
nPCs = 6; %number of PCs to use

%find number of trials per comparison
nTrialPerStim = zeros(1,nDirs);
for iDir = 1:nDirs
    nTrialPerStim(iDir) = length(find(tDir == Dirs(iDir)));
end
for i = 1:nDirs
    for j = 1:nDirs
        trial_mat(i,j) = nTrialPerStim(i)+nTrialPerStim(j);
    end
end
max_trialn = max(trial_mat(:));

%set of comparisons
nDiff = nDirs/2;
Diffs = Dirs(2:nDirs/2 + 1)-Dirs(1);

%preallocate outputs
y_all = nan(max_trialn,nDirs,nDiff,2); %tracks predicted and actual dir  
pct_corr = nan(nDiff,1); %percent correct from binomial
err = nan(nDiff,2); %error from binomial

for iDiff = 1:nDirs/2 %loop through direction comparisons
    for iDir = 1:nDirs %loop through all pairs of directions
        ind_1 = find(tDir == Dirs(iDir)); %find trials for dir 1
        if iDir+iDiff<=nDirs %find trials for dir 2
            ind_2 = find(tDir == Dirs(iDir+iDiff));
        else
            ind_2 = find(tDir == Dirs(iDir-(nDirs-iDiff)));
        end
        for i = 1:length([ind_1 ind_2]) %loop through all trials of type with hold-one-out
            if i<=length(ind_1) %if holding out trial from dir 1
                ind_1_use = setdiff(ind_1,ind_1(i)); %hold out
                ind_2_use = ind_2; 
                ind_test = ind_1(i); %held out
                t = 0; %tracks held out trial type
            else %if holding out trial from dir 2
                ind_2_use = setdiff(ind_2,ind_2(i-length(ind_1))); %hold out
                ind_1_use = ind_1;
                ind_test = ind_2(i-length(ind_1)); %held out
                t = 1; %tracks held out trial type
            end
            x = [s(ind_1_use,1:nPCs); s(ind_2_use,1:nPCs)]; %responses for training trials (in pc space)
            y = [zeros(length(ind_1_use),1); ones(length(ind_2_use),1)]; %directions for training trials
            C=eye(size(x,2)); %transform for GLM
            x_trans = x*C;
            [~,~,dirGLM] = glmfit(x_trans,y,'binomial'); %train
            yhat = glmval(dirGLM.beta,s(ind_test,1:nPCs),'logit'); %test
            y_all(i,iDir,iDiff,1) = yhat>0.5; %threshold for prediction of test trial outcome
            y_all(i,iDir,iDiff,2) = t; %true test trial outcome
        end
    end
    y_test = y_all(:,:,iDiff,1) == y_all(:,:,iDiff,2); %comparison of predicted vs actual
    y_test = double(y_test); %conversion to incorporate nans on next line
    y_test(find(isnan(y_all(:,:,iDiff,1)))) = y_all(find(isnan(y_all(:,:,iDiff,1)))); %dealing with nans for different trial numbers
    [pct_corr(iDiff,1) err(iDiff,:)] = binofit(sum(y_test(:),'omitnan'),sum(~isnan(y_test(:)))); %binomial computation of pct correct and error
end
figure; errorbar(Diffs,pct_corr,pct_corr-err(:,1),err(:,2)-pct_corr);
