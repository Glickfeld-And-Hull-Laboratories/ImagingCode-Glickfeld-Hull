%create templates for fitting
if (length(Cond_USE))== nCond
    Cond = 'all_cond';
    SF_vec0 = flipud(uvar(:,2));
    TF_vec0 = uvar(:,1);
else
    Cond = 'subset_cond';
    SF_vec0 = [.32 .08 .02]; %flipped to have low to high SF in square  %flipud
    TF_vec0 = [1 4 15];
end

[tftf,sfsf] = meshgrid(TF_vec0,SF_vec0); 
grid2.sfsf = sfsf;
grid2.tftf = tftf;

dSF = median(diff(log2(SF_vec0)));
dTF = median(diff(log2(TF_vec0)));
SF_vec00 = log2(SF_vec0(1)):(dSF/10):log2(SF_vec0(end));
TF_vec00 = log2(TF_vec0(1)):(dTF/10):log2(TF_vec0(end));
[sfsf00,tftf00]=meshgrid(SF_vec00,TF_vec00);
grid2.sfsf00 = sfsf00;
grid2.tftf00 = tftf00;


%test for significant dF/F
    alphaB = .05./(nCond);
    Info_ttest_mat = zeros(size(roi_avg,2),nCond);
    for iCell = 1:size(roi_avg,2);
        start = 1;
        p_ttestB = zeros(1,nCond);
        for iCond = 1:nCond
            nRep = stim_reps(1,iCond);
            resp_off_sq = squeeze(resp_off);
            if size(roi_avg,2)==1
                resp_off_sq =resp_off_sq';
            end
            [h_ttestB1,p_ttestB1] = ttest(resp_off_sq(iCell,start:start-1+nRep),resp_on(iCell,start:start-1+nRep),alphaB,'left');
            p_ttestB(iCond) = p_ttestB1;
            start = start+nRep;
        end   
        Info_ttest_mat(iCell,:) = p_ttestB;
    end
    H_ttest = zeros(size(roi_avg,2),3);
    H_ttest(:,1) = min(Info_ttest_mat,[],2) < alphaB;
    if run ==1;
        Info_ttest_mat_run = zeros(size(roi_avg,2),nCond);
        Info_ttest_mat_norun = zeros(size(roi_avg,2),nCond);        
        for iCell = 1:size(roi_avg,2);
            start_run = 1;
            start_norun = 1;
            p_ttest_run = zeros(1,nCond);
            p_ttest_norun = zeros(1,nCond);
            resp_off_run_sq = squeeze(resp_off_run);
            resp_off_norun_sq = squeeze(resp_off_norun);
            if size(roi_avg,2)==1
                resp_off_run_sq =resp_off_run_sq';
                resp_off_norun_sq =resp_off_norun_sq';
                resp_on_run_sq =resp_on_run';
                resp_on_norun_sq =resp_on_norun';
            else
                resp_on_run_sq =resp_on_run;
                resp_on_norun_sq =resp_on_norun;
            end
            for iCond = 1:nCond
                nRep_run = stim_reps_run(1,iCond);
                nRep_norun = stim_reps_norun(1,iCond);
                [h_ttest_B1,p_ttest_B1] = ttest(resp_off_run_sq(iCell,start_run:start_run-1+nRep_run),resp_on_run_sq(iCell,start_run:start_run-1+nRep_run),alphaB,'left');
                p_ttest_run(1,iCond) = p_ttest_B1;
                [h_ttest_B1,p_ttest_B1] = ttest(resp_off_norun_sq(iCell,start_norun:start_norun-1+nRep_norun),resp_on_norun_sq(iCell,start_norun:start_norun-1+nRep_norun),alphaB,'left');
                p_ttest_norun(1,iCond) = p_ttest_B1;
                start_run = start_run+nRep_run;
                start_norun = start_norun+nRep_norun;
            end   
        Info_ttest_mat_run(iCell,:) = p_ttest_run;
        Info_ttest_mat_norun(iCell,:) = p_ttest_norun;
        end
        H_ttest(:,2) = min(Info_ttest_mat_run,[],2) < alphaB;
        H_ttest(:,3) = min(Info_ttest_mat_norun,[],2) < alphaB;
    end 
    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_POST' num2str(post_win) '_ttest.mat']);   
    save(fn_out, 'H_ttest')
%alpha2 = .05./(nCond); %Bonferroni correction, doesn't require indep; don't include blank in deg freedom

Mask_cells_USE = ones(1,size(roi_avg,2));
cells_USE = sum(Mask_cells_USE);

%choose data- add running sorting etc in here later
str_run = strvcat('allrun', 'running', 'norun');
Ind_struct = [];
start = 1;
for iCond = 1:nCond
    nRep = stim_reps(iCond);
    Ind_struct(iCond).all_trials = [start:start-1+nRep];
    start = start+nRep;
end
if run ==1;
    fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_running_ind.mat']);
    load(fn);
    start_run = 1;
    start_norun = 1;
    for iCond = 1:nCond
        nRep_run = stim_reps_run(iCond);
        nRep_norun = stim_reps_norun(iCond);
        Ind_struct(iCond).run_trials = run_ind(start_run:start_run-1+nRep_run);    
        Ind_struct(iCond).norun_trials = norun_ind(start_norun:start_norun-1+nRep_norun);
        start_run = start_run+nRep_run;
        start_norun = start_norun+nRep_norun;
    end
end
%fit data for real data and n shuffles
Fit_struct = [];
CM_mat = [];
for count_shuf = 0:Nshuf
    fprintf('.')
    if run == 1
        nRunning = 3;
    else
        nRunning = 1;
    end
    Im_mat_USE = zeros(size(roi_avg,2), length(Cond_USE), nRunning);
    Im_mat_std = zeros(size(roi_avg,2), length(Cond_USE), nRunning);
    for iCond = 1:length(Cond_USE)        
        ind_all = Ind_struct(Cond_USE(iCond)).all_trials;
        if run == 1;
            ind_run = Ind_struct(Cond_USE(iCond)).run_trials;
            ind_norun = Ind_struct(Cond_USE(iCond)).norun_trials;
        end
        if count_shuf > 0 %resample with replacement, don't resample by trial for now because running-rejection may be uneven for various trials..
            ind_all_1 = ind_all(randsample(length(ind_all),length(ind_all),1));
            if run == 1;
                ind_run_1 = ind_run(randsample(length(ind_run),length(ind_run),1));
                ind_norun_1 = ind_norun(randsample(length(ind_norun),length(ind_norun),1));
            end
        else
            ind_all_1 = ind_all;        
            if run == 1
                ind_run_1 = ind_run;
                ind_norun_1 = ind_norun;
            end
        end
        if P == 1;
            Im_mat_USE(:,iCond,1) = mean(resp_dF(:,ind_all_1),2);
            Im_mat_std(:,iCond,1) = std(resp_dF(:,ind_all_1),[],2);
            if run == 1;
                Im_mat_USE(:,iCond,2) = mean(resp_dF(:,ind_run_1),2);
                Im_mat_std(:,iCond,2) = std(resp_dF(:,ind_run_1),[],2);
                Im_mat_USE(:,iCond,3) = mean(resp_dF(:,ind_norun_1),2);
                Im_mat_std(:,iCond,3) = std(resp_dF(:,ind_norun_1),[],2);
            end
        elseif P==2;
            Im_mat_USE(:,iCond,1) = mean(resp_norm(:,ind_all_1),2);
            Im_mat_std(:,iCond,1) = std(resp_norm(:,ind_all_1),[],2);
            if run == 1;
                Im_mat_USE(:,iCond,2) = mean(resp_norm(:,ind_run_1),2);
                Im_mat_std(:,iCond,2) = std(resp_norm(:,ind_run_1),[],2);
                Im_mat_USE(:,iCond,3) = mean(resp_norm(:,ind_norun_1),2);
                Im_mat_std(:,iCond,3) = std(resp_norm(:,ind_norun_1),[],2);
            end
        end
    end
    for iRun = 1:nRunning;
        figure;
        start = 1;
        for iCell = 1:nCells;
            if H_ttest(iCell,iRun) == 1
                a = Im_mat_USE(iCell,:,iRun);
                b = reshape(a',length(SF_vec0),length(TF_vec0));
                %b2 = b( ind_SFuse(:,1),ind_TFuse(:,1));
                data = b';
                ind0 = find(data<0);
                data(ind0) = NaN;
                if count_shuf == 0
                    PLOTIT_FIT = 1;
                    SAVEALLDATA = 1;
                    Fit_2Dellipse_LG
                    eval(['Fit_struct(iCell).True.s_',str_run(iRun,:),' = s;']);
                else
                    SAVEALLDATA = 0;
                    PLOTIT_FIT = 0;
                    Fit_2Dellipse_LG
                    eval(['Fit_struct(iCell).Shuf(count_shuf).s_',str_run(iRun,:),' = s;']);
                end
            else
                s.CM_data = zeros(1,2);
                s.x = zeros(1,6);
                if nCond == 9;
                    oversamp = 33;
                elseif nCond == 25;
                    oversamp = 65;
                elseif nCond == 50;
                    oversamp = 65;
                end
                s.k2_plot_oversamp= zeros(oversamp);
                eval(['Fit_struct(iCell).True.s_',str_run(iRun,:),' = s;']);
            end
        end
    end
end

fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_POST' num2str(post_win) '_Fit_struct.mat']);   
save(fn_out, 'Fit_struct')


%test for goodness of fit

if Nshuf>1;
    for iCell = 1:cells_USE
        if Mask_cells_USE(iCell) == 1
            if ~isempty(Fit_struct(iCell).True)                
                eval(['tmp = Fit_struct(iCell).True.s_',str_run,'.x;']);
                eval(['tmp = [tmp Fit_struct(iCell).True.s_',str_run,'.SFhicut_50];']);
                eval(['tmp = [tmp Fit_struct(iCell).True.s_',str_run,'.TFhicut_50];']);
                eval(['tmp = [tmp Fit_struct(iCell).True.s_',str_run,'.SFhicut_10];']);
                eval(['tmp = [tmp Fit_struct(iCell).True.s_',str_run,'.TFhicut_10];']);
                fit_true_vec(iCell,:) = tmp;
            end
        end
    end

    for count_shuf = 1:Nshuf
        for iCell = 1:cells_USE
            if Mask_cells_USE(iCell) == 1
                if ~isempty(Fit_struct(iCell).Shuf)
                    eval(['tmp = Fit_struct(iCell).Shuf(count_shuf).s_',str_run,'.x;']);
                    eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_',str_run,'.SFhicut_50];']);
                    eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_',str_run,'.TFhicut_50];']);
                    eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_',str_run,'.SFhicut_10];']);
                    eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_',str_run,'.TFhicut_10];']);
                    %fit is: %A sigma_SF sigma_TF sf0 tf0 xi
                    fit_shuf_vec(iCell,:,count_shuf) = tmp;
                end
            end
        end
    end

    Npars = size(fit_shuf_vec,2);
    lbub_fits = zeros(cells_USE,Npars,5);
    alpha_bound = .025;
    for iCell = 1:cells_USE
        if Mask_cells_USE(iCell) == 1
            for count2 = 1:Npars
                tmp = squeeze(fit_shuf_vec(iCell,count2,:));
                [i,j] = sort(tmp);
                ind_shuf_lb = ceil(Nshuf*alpha_bound);
                ind_shuf_ub = ceil(Nshuf*(1-alpha_bound));
                lbub_fits(iCell,count2,1) = i(ind_shuf_lb);
                lbub_fits(iCell,count2,2) = i(ind_shuf_ub);
                lbub_fits(iCell,count2,3) = mean(i); 
                lbub_fits(iCell,count2,5) = std(i);
            end
            %now take means from truedata fit:
            lbub_fits(iCell,:,4) = fit_true_vec(iCell,:);
        end 
    end
end