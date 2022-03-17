close all
clear all
clear all global
clc
dataset = 'plaidDiscrim_exptList';
eval(dataset);
nexp = size(expt,2);
used = nan(1,nexp);
for iexp = 1:nexp
    LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
    if strcmp(expt(iexp).folder,'lindsey')
        data_base = LG_base;
    elseif strcmp(expt(iexp).folder,'camaron')
        data_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron';
    end
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    run_str = ['runs'];
    nrun = size(expt(iexp).runs,1);
    for irun = 1:nrun
        run_str = [run_str '-' expt(iexp).runs(irun,:)];
    end
    load(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
    if sum(~celleqel2mat_padded(input.tDoFeedbackMotion))
        fprintf([mouse ' ' date '\n'])
        load(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirAnalysis.mat']))
        load(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
        load(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
        tDirCat = zeros(size(tDir));
        tDirCat(find(tDir>45)) = 1;
        tDirCat(find(tDir==45)) = NaN;
        trial_all = find(~isnan(tDir));
       
        z_stim_resp = zscore(data_stim_resp');
        ind_pyr = intersect(find(~mask_label),resp_ind);
        [u pcs lat t exp] = pca(z_stim_resp(:,ind_pyr));
        %ind_use = 1:find(cumsum(exp)>70,1,'first');
        nboot = 100;
        if length(ind_pyr)>8
            ind_use = 1:8;
        else
            ind_use = 1:length(ind_pyr);
        end
        glm_out = zeros(length(ind_use)+1,nboot);
        glm_grating_out = zeros(length(ind_use)+1,nboot);
        glm_plaid_out = zeros(length(ind_use)+1,nboot);
        alltrain_traintest_pctCorr = zeros(1,nboot);
        alltrain_alltest_pctCorr = zeros(1,nboot);
        
        gratingtrain_traintest_pctCorr = zeros(1,nboot);
        gratingtrain_gratingtest_pctCorr = zeros(1,nboot);
        gratingtrain_plaidtest_pctCorr = zeros(1,nboot);
        plaidtrain_traintest_pctCorr = zeros(1,nboot);
        plaidtrain_plaidtest_pctCorr = zeros(1,nboot);
        plaidtrain_gratingtest_pctCorr = zeros(1,nboot);
        
        trial_train = intersect(find(~isnan(tDirCat)),trial_all(randsample(length(trial_all),ceil(length(trial_all)/2))));
        trial_test = intersect(find(~isnan(tDirCat)),setdiff(trial_all,trial_train));
        trial_grating = intersect(find(~tMaskTrial),find(~isnan(tDirCat)));
        trial_grating_train = trial_grating(randsample(length(trial_grating),ceil(length(trial_grating)/2)));
        trial_gratingtrain_grating_test = setdiff(trial_grating,trial_grating_train);
        trial_gratingtrain_plaid_test = intersect(find(tMaskTrial),find(~isnan(tDirCat)));
        trial_plaid = intersect(find(tMaskTrial),find(~isnan(tDirCat)));
        trial_plaid_train = trial_plaid(randsample(length(trial_plaid),ceil(length(trial_plaid)/2)));
        trial_plaidtrain_plaidtest = setdiff(trial_plaid,trial_plaid_train);
        trial_plaidtrain_grating_test = intersect(find(~tMaskTrial),find(~isnan(tDirCat)));
        
        if length(ind_pyr)>5
            for iboot = 1:nboot
%                 if length(ind_pyr)>12
%                     ind_use = ind_pyr(randsample(length(ind_pyr),12));
%                 else
%                     ind_use = ind_pyr;
%                 end
                
                trial_train = intersect(find(~isnan(tDirCat)),trial_all(randsample(length(trial_all),ceil(length(trial_all)/2))));
                trial_test = intersect(find(~isnan(tDirCat)),setdiff(trial_all,trial_train));
                trial_grating = intersect(find(~tMaskTrial),find(~isnan(tDirCat)));
                trial_grating_train = trial_grating(randsample(length(trial_grating),ceil(length(trial_grating)/2)));
                trial_gratingtrain_grating_test = setdiff(trial_grating,trial_grating_train);
                trial_gratingtrain_plaid_test = intersect(find(tMaskTrial),find(~isnan(tDirCat)));
                trial_plaid = intersect(find(tMaskTrial),find(~isnan(tDirCat)));
                trial_plaid_train = trial_plaid(randsample(length(trial_plaid),ceil(length(trial_plaid)/2)));
                trial_plaidtrain_plaidtest = setdiff(trial_plaid,trial_plaid_train);
                trial_plaidtrain_grating_test = intersect(find(~tMaskTrial),find(~isnan(tDirCat)));
                
%                 glm_out(:,iboot) = glmfit(z_stim_resp(trial_train,ind_use),tDirCat(:,trial_train)','binomial');
%                 y_hat_train = glmval(glm_out(:,iboot),z_stim_resp(trial_train,ind_use),'logit');
%                 y_hat_test = glmval(glm_out(:,iboot),z_stim_resp(trial_test,ind_use),'logit');
                
                glm_out(:,iboot) = glmfit(pcs(trial_train,ind_use),tDirCat(:,trial_train)','binomial');
                y_hat_train = glmval(glm_out(:,iboot),pcs(trial_train,ind_use),'logit');
                y_hat_test = glmval(glm_out(:,iboot),pcs(trial_test,ind_use),'logit');
                alltrain_traintest_pctCorr(:,iboot) = sum(((y_hat_train>0.5)-tDirCat(trial_train)')==0)./length(trial_train);
                alltrain_alltest_pctCorr(:,iboot) = sum(((y_hat_test>0.5)-tDirCat(trial_test)')==0)./length(trial_test);

%                 glm_grating_out(:,iboot) = glmfit(z_stim_resp(trial_grating_train,ind_use),tDirCat(:,trial_grating_train)','binomial');
%                 y_hat_grating_train = glmval(glm_grating_out(:,iboot),z_stim_resp(trial_grating_train,ind_use),'logit');
%                 y_hat_grating_test = glmval(glm_grating_out(:,iboot),z_stim_resp(trial_gratingtrain_grating_test,ind_use),'logit');
%                 y_hat_plaid_test = glmval(glm_grating_out(:,iboot),z_stim_resp(trial_gratingtrain_plaid_test,ind_use),'logit');

                glm_grating_out(:,iboot) = glmfit(pcs(trial_grating_train,ind_use),tDirCat(:,trial_grating_train)','binomial');
                y_hat_grating_train = glmval(glm_grating_out(:,iboot),pcs(trial_grating_train,ind_use),'logit');
                y_hat_grating_test = glmval(glm_grating_out(:,iboot),pcs(trial_gratingtrain_grating_test,ind_use),'logit');
                y_hat_plaid_test = glmval(glm_grating_out(:,iboot),pcs(trial_gratingtrain_plaid_test,ind_use),'logit');
                
                gratingtrain_traintest_pctCorr(:,iboot) = sum(((y_hat_grating_train>0.5)-tDirCat(trial_grating_train)')==0)./length(trial_grating_train);
                gratingtrain_gratingtest_pctCorr(:,iboot) = sum(((y_hat_grating_test>0.5)-tDirCat(trial_gratingtrain_grating_test)')==0)./length(trial_gratingtrain_grating_test);
                gratingtrain_plaidtest_pctCorr(:,iboot) = sum(((y_hat_plaid_test>0.5)-tDirCat(trial_gratingtrain_plaid_test)')==0)./length(trial_gratingtrain_plaid_test);
                                
%                 glm_plaid_out(:,iboot) = glmfit(z_stim_resp(trial_plaid_train,ind_use),tDirCat(:,trial_plaid_train)','binomial');
%                 y_hat_plaid_train = glmval(glm_plaid_out(:,iboot),z_stim_resp(trial_plaid_train,ind_use),'logit');
%                 y_hat_plaid_test = glmval(glm_plaid_out(:,iboot),z_stim_resp(trial_plaidtrain_plaidtest,ind_use),'logit');
%                 y_hat_grating_test = glmval(glm_plaid_out(:,iboot),z_stim_resp(trial_plaidtrain_grating_test,ind_use),'logit');

                glm_plaid_out(:,iboot) = glmfit(pcs(trial_plaid_train,ind_use),tDirCat(:,trial_plaid_train)','binomial');
                y_hat_plaid_train = glmval(glm_plaid_out(:,iboot),pcs(trial_plaid_train,ind_use),'logit');
                y_hat_plaid_test = glmval(glm_plaid_out(:,iboot),pcs(trial_plaidtrain_plaidtest,ind_use),'logit');
                y_hat_grating_test = glmval(glm_plaid_out(:,iboot),pcs(trial_plaidtrain_grating_test,ind_use),'logit');
                
                plaidtrain_traintest_pctCorr(:,iboot) = sum(((y_hat_plaid_train>0.5)-tDirCat(trial_plaid_train)')==0)./length(trial_plaid_train);
                plaidtrain_plaidtest_pctCorr(:,iboot) = sum(((y_hat_plaid_test>0.5)-tDirCat(trial_plaidtrain_plaidtest)')==0)./length(trial_plaidtrain_plaidtest);
                plaidtrain_gratingtest_pctCorr(:,iboot) = sum(((y_hat_grating_test>0.5)-tDirCat(trial_plaidtrain_grating_test)')==0)./length(trial_plaidtrain_grating_test);
            end
            save(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_popDecode.mat']),'glm_out','alltrain_traintest_pctCorr','alltrain_alltest_pctCorr','glm_grating_out','gratingtrain_traintest_pctCorr','gratingtrain_gratingtest_pctCorr','gratingtrain_plaidtest_pctCorr','plaidtrain_traintest_pctCorr','plaidtrain_plaidtest_pctCorr','plaidtrain_gratingtest_pctCorr')
        end
    end
end