clear all
areas = ['PM'; 'LM'; 'AL'];
inj = 'V1';
P = 2;
nPlanes = 1;
begin = 1;

for iArea = 2;
    matrix = 'SFTF_DIR_CON';
    image = areas(iArea,:);
    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);
    
    for iexp = 1:nexp
        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        run = exp_list.run_mat{iexp};
        blanks = exp_list.blanks_mat{iexp};
        SFs = exp_list.SF_mat{iexp};
        TFs = exp_list.TF_mat{iexp};
        dirs = exp_list.dirs_mat{iexp};
        nframes = exp_list.nframes_mat{iexp};
        
        nON = nframes;
        nOFF = nframes;
        if nframes == 12;
            pre_win = [7 12];
            post_win = [13 24];
        elseif nframes == 15;
            pre_win = [8 15];
            post_win = [16 30];
        end
        
        nSFTF = length(SFs).*length(TFs);
        nCond = dirs.*nSFTF;
        
        step = 360./dirs;
        ORI_vec0 = [0:step:(15.*step)]*pi/180;
        
        base = 'G:\users\lindsey\analysisLG\active mice';    
        outDir = fullfile(base, mouse,date);
        
        userun = exp_list.sftf_runs_mat{iexp};
        fn_local = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_local_max.mat']);
        load(fn_local);
        
        %direction runs
        userun = exp_list.dir_runs_mat{iexp};
        count_protocol = exp_list.dir_prot_mat{iexp};
        
        fn_reps= fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_reps.mat']);
        load(fn_reps);
        
        fn_stack = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_sorted.tif']);
        stack_sorted = readtiff(fn_stack);
        
        nReps = sum(stim_reps(1,:));
        siz = size(stack_sorted);
        roi_stim = zeros(siz(1), siz(2), nOFF+nON,nReps);
        start = 1;
        rep = 1;
        for iRep = 1:nReps;
            roi_stim(:,:,:, rep) = stack_sorted(:,:,start:start-1+nON+nOFF);
            start = start+nON+nOFF;
            rep = rep+1;
        end
        stim_off = squeeze(mean(roi_stim(:,:,pre_win(1):pre_win(2),:),3));
        stim_on = squeeze(mean(roi_stim(:,:,post_win(1):post_win(2),:),3));
        
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_allstim.mat']);
        save(fn_out, 'stim_on', 'stim_off');
        
        nRep = stim_reps(1,1);
        stack_var=var(stim_on,[],3);
        stack_mean=mean(stim_on,3);
        stack_mean_long = reshape(stack_mean, [siz(1)*siz(2) 1]);
        stack_var_long = reshape(stack_var, [siz(1)*siz(2) 1]);
        b=robustfit(stack_mean_long,stack_var_long);

        nPhoton=round((stack_mean_long./b(2)).*nRep.*(nON/2.667));
        fn_stim = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_photon.mat']);
        save(fn_stim, 'b', 'nPhoton');
          
        resp_dF = zeros(n_pix, size(stim_on,3));
        resp_dFoverF = zeros(n_pix, size(stim_on,3));
        for ipix = 1:n_pix
            sub_y = [i(ipix)-1:i(ipix)+1]; 
            sub_x = [j(ipix)-1:j(ipix)+1];
            roi_on = squeeze(mean(mean(stim_on(sub_y,sub_x,:),2),1));
            roi_off = squeeze(mean(mean(stim_off(sub_y,sub_x,:),2),1));
            resp_dF(ipix,:) = (roi_on-roi_off)';
            resp_dFoverF(ipix,:) = ((roi_on-roi_off)./roi_off)';
        end
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp.mat']);
        save(fn_out, 'resp_dF', 'resp_dFoverF');
        clear stim_off
        clear stim_on
        clear roi_stim
        
        %contrast runs
        userun = exp_list.con_runs_mat{iexp};
        count_protocol = exp_list.con_prot_mat{iexp};
        
        fn_reps= fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_reps.mat']);
        load(fn_reps);
        
        fn_stack = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_sorted.tif']);
        stack_sorted = readtiff(fn_stack);
        
        nReps = sum(stim_reps(1,:));
        siz = size(stack_sorted);
        roi_stim = zeros(siz(1), siz(2), nOFF+nON,nReps);
        start = 1;
        rep = 1;
        for iRep = 1:nReps;
            roi_stim(:,:,:, rep) = stack_sorted(:,:,start:start-1+nON+nOFF);
            start = start+nON+nOFF;
            rep = rep+1;
        end
        stim_off = squeeze(mean(roi_stim(:,:,pre_win(1):pre_win(2),:),3));
        stim_on = squeeze(mean(roi_stim(:,:,post_win(1):post_win(2),:),3));
        
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_allstim.mat']);
        save(fn_out, 'stim_on', 'stim_off');
        
        nRep = stim_reps(1,1);
        stack_var=var(stim_on,[],3);
        stack_mean=mean(stim_on,3);
        stack_mean_long = reshape(stack_mean, [siz(1)*siz(2) 1]);
        stack_var_long = reshape(stack_var, [siz(1)*siz(2) 1]);
        b=robustfit(stack_mean_long,stack_var_long);

        nPhoton=round((stack_mean_long./b(2)).*nRep.*(nON/2.667));
        fn_stim = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_photon.mat']);
        save(fn_stim, 'b', 'nPhoton');
          
        resp_dF = zeros(n_pix, size(stim_on,3));
        resp_dFoverF = zeros(n_pix, size(stim_on,3));
        for ipix = 1:n_pix
            sub_y = [i(ipix)-1:i(ipix)+1]; 
            sub_x = [j(ipix)-1:j(ipix)+1];
            roi_on = squeeze(mean(mean(stim_on(sub_y,sub_x,:),2),1));
            roi_off = squeeze(mean(mean(stim_off(sub_y,sub_x,:),2),1));
            resp_dF(ipix,:) = (roi_on-roi_off)';
            resp_dFoverF(ipix,:) = ((roi_on-roi_off)./roi_off)';
        end
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp.mat']);
        save(fn_out, 'resp_dF', 'resp_dFoverF');
        clear stim_off
        clear stim_on
        clear roi_stim
    end
end
