clear all
areas = ['PM'; 'LM'; 'AL'; 'RL'; 'AM'];
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

for iArea = 1:5;
    image = areas(iArea,:);
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);
    mouse_list = [];
    mice = exp_list.mouse_mat;
    all_fits(iArea).nexp = nexp;
    all_fits(iArea).name = image;
    all_boutons = [0 0];
    speed_list = [];
    for iexp = 1:nexp;
        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        userun = exp_list.runs_mat{iexp};
        count_prot = exp_list.prot_mat{iexp};
        run = exp_list.run_mat{iexp};
        blanks = exp_list.blanks_mat{iexp};
        dirs = exp_list.dir_mat{iexp};
        zoom = exp_list.zoom_mat{iexp};
        
        base = 'G:\users\lindsey\analysisLG\active mice';    
        outDir = fullfile(base, mouse,date);
        
        all_fits(iArea).expt(iexp).mouse = mouse;
        all_fits(iArea).expt(iexp).date = date;
        all_fits(iArea).expt(iexp).userun = userun;
        all_fits(iArea).expt(iexp).count_prot = count_prot;
        all_fits(iArea).expt(iexp).dirs = dirs;
        all_fits(iArea).expt(iexp).zoom = zoom;

        fn_lbub = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']);
        load(fn_lbub);
        fn_fit  = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_Fit_struct.mat']);
        load(fn_fit); 
        fn_local  = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_local_max.mat']);
        load(fn_local); 
        fn_resp = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp.mat']);
        load(fn_resp);
        fn_reps= fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_reps.mat']);
        load(fn_reps);
        
        nReps = sum(stim_reps(1,:));
        if dirs == 2
            nCond = 50;
            add = 1;
            stim_reps_dir = zeros(1,26);
            for iCond = 1:nCond/2
            stim_reps_dir(1,iCond) = sum(stim_reps(1,add:add+1));
            add = add+2;
            end
            stim_reps_dir(1,end) = stim_reps(1,end);
            stim_reps = stim_reps_dir;
        end
        
        dF_mat = zeros(size(resp_dFoverF,1), 25);
        start = 1;
        for iCond = 1:25
            nRep = stim_reps(iCond);
            dF_mat(:,iCond) = mean(resp_dF(:,start:start-1+nRep),2);
            start = start+nRep;
        end
        
        all_fits(iArea).expt(iexp).n = [size(lbub_fits,1) length(goodfit_ind)];
        
        all_boutons = all_boutons+all_fits(iArea).expt(iexp).n;
        
        for iCell = 1:size(lbub_fits,1)
            all_fits(iArea).expt(iexp).bouton(iCell).dF_fit = Fit_struct(iCell).True.s_.x(:,1);
            all_fits(iArea).expt(iexp).bouton(iCell).sigma_SF = Fit_struct(iCell).True.s_.x(:,2);
            all_fits(iArea).expt(iexp).bouton(iCell).sigma_TF = Fit_struct(iCell).True.s_.x(:,3);
            all_fits(iArea).expt(iexp).bouton(iCell).SF_fit = 2.^Fit_struct(iCell).True.s_.x(:,4);
            all_fits(iArea).expt(iexp).bouton(iCell).TF_fit = 2.^Fit_struct(iCell).True.s_.x(:,5);
            all_fits(iArea).expt(iexp).bouton(iCell).xi_fit = Fit_struct(iCell).True.s_.x(:,6);
            all_fits(iArea).expt(iexp).bouton(iCell).speed = 2.^lbub_fits(iCell,5,4)./2.^lbub_fits(iCell,4,4);
            all_fits(iArea).expt(iexp).bouton(iCell).pos = [i(iCell, :) j(iCell, :)];
            all_fits(iArea).expt(iexp).bouton(iCell).plotfit = Fit_struct(iCell).True.s_.k2b_plot;
            all_fits(iArea).expt(iexp).bouton(iCell).dFoverF = Fit_struct(iCell).True.s_.orig;
            all_fits(iArea).expt(iexp).bouton(iCell).dF = dF_mat';
            if find(goodfit_ind == iCell)>0
                all_fits(iArea).expt(iexp).bouton(iCell).goodfit = 1;
                speed_list = [speed_list all_fits(iArea).expt(iexp).bouton(iCell).speed];
            else
                all_fits(iArea).expt(iexp).bouton(iCell).goodfit = 0;
            end
        end
        all_fits(iArea).expt(iexp).median_speed = median(speed_list);
    end
    all_fits(iArea).n = all_boutons;
end
    fn_out = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
    save(fn_out, 'all_fits');
        
        