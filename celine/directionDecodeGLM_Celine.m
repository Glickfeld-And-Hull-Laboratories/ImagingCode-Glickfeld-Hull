%will use trialResp, which has a cell for each day
for id = 1:nd
    %extract putative pyr cells from response matrix for all trials
    trial_dfof = trialResp{id}(:,green_ind_keep);
    trial_dfof(find(RIx{id}),:) = NaN;
    
    %transform into pca space for dimensionality reduction
    [u s v] = pca(trial_dfof);
    nPCs =6; %number of PCs to use
    for iCon = 1:nCon    
        ind_con1 = find(tCon_match{id}==cons(iCon));
        ind_con = intersect(ind_con1,find(~RIx{id}));

        %find number of trials per comparison
        nTrialPerStim = zeros(1,nDir);
        for iDir = 1:nDir
            nTrialPerStim(iDir) = length(find(tDir == dirs(iDir)));
        end
        
        for i = 1:nDir
            for j = 1:nDir
                trial_mat(i,j) = nTrialPerStim(i)+nTrialPerStim(j);
            end
        end
        max_trialn = max(trial_mat(:));
        
        %set of comparisons
        nDiff = nDir/2;
        Diffs = dirs(2:nDir/2 + 1)-dirs(1);
        
        %preallocate outputs
        y_all = nan(max_trialn,nDir,nDiff,2); %tracks predicted and actual dir  
        pct_corr = nan(nDiff,1); %percent correct from binomial
        err = nan(nDiff,2); %error from binomial
        figure;
        for iDiff = 1:nDir/2 %loop through direction comparisons
            for iDir = 1:nDir %loop through all pairs of directions
                ind_1_temp = find(tDir_match{id} == dirs(iDir)); %find trials for dir 1
                ind_1 = intersect(ind_1_temp, ind_con);
                if iDir+iDiff<=nDir %find trials for dir 2
                    ind_2_temp = find(tDir_match{id} == dirs(iDir+iDiff));
                else
                    ind_2_temp = find(tDir_match{id} == dirs(iDir-(nDir-iDiff)));
                end
                ind_2 = intersect(ind_2_temp, ind_con);
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
                    x = [s(ind_1_use,2:nPCs+1); s(ind_2_use,2:nPCs+1)]; %responses for training trials (in pc space)
                    y = [zeros(length(ind_1_use),1); ones(length(ind_2_use),1)]; %directions for training trials
                    C=eye(size(x,2)); %transform for GLM
                    x_trans = x*C;
                    [~,~,dirGLM] = glmfit(x_trans,y,'binomial'); %train
                    yhat = glmval(dirGLM.beta,s(ind_test,2:nPCs+1),'logit'); %test
                    y_all(i,iDir,iDiff,1) = yhat>0.5; %threshold for prediction of test trial outcome
                    y_all(i,iDir,iDiff,2) = t; %true test trial outcome
                    
                end
                scatter(x(:,1),x(:,2)); hold on;
            end
            if id == 1
                title(['Post-DART, contrast ' num2str(cons(iCon))])
            else
                title(['Baseline, contrast ' num2str(cons(iCon))])
            end
            xlabel('PC 2')
            ylabel('PC 3')
            print(['sess_',num2str(id),'_',num2str(cons(iCon)),'_PCs.pdf'],'-dpdf');
            y_test = y_all(:,:,iDiff,1) == y_all(:,:,iDiff,2); %comparison of predicted vs actual
            y_test = double(y_test); %conversion to incorporate nans on next line
            y_test(find(isnan(y_all(:,:,iDiff,1)))) = y_all(find(isnan(y_all(:,:,iDiff,1)))); %dealing with nans for different trial numbers
            [pct_corr(iDiff,1) err(iDiff,:)] = binofit(sum(y_test(:),'omitnan'),sum(~isnan(y_test(:)))); %binomial computation of pct correct and error
        end


        figure; errorbar(Diffs,pct_corr,pct_corr-err(:,1),err(:,2)-pct_corr);
        ylabel('percent correct')
        xlabel('direction difference (deg)')
        ylim([0 1.1])
        xticks([45, 90, 135, 180])
        set(gca, 'TickDir', 'out')
        box off
        if id == 1
            title(['Post-DART, contrast ' num2str(cons(iCon))])
        else
            title(['Baseline, contrast ' num2str(cons(iCon))])
        end
        print(['sess_',num2str(id),'_',num2str(cons(iCon)),'_decoding.pdf'],'-dpdf');
    end
end
