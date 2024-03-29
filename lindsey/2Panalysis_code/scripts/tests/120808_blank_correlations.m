areas = ['PM'; 'AL'; 'LM'];
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
for iArea = 3
    image = areas(iArea,:);
    
    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);
    for iexp = 1:nexp;
    
        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        userun = exp_list.runs_mat{iexp};
        count_protocol = exp_list.prot_mat{iexp};
        run = exp_list.run_mat{iexp};
        blanks = exp_list.blanks_mat{iexp};
        dirs = exp_list.dir_mat{iexp};

        base = 'G:\users\lindsey\analysisLG\active mice';    
        outDir = fullfile(base, mouse,date);
        fn_resp = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_resp.mat']);
        load(fn_resp);
        fn_reps = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_reps.mat']);
        load(fn_reps);
        fn_lbub = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']);
        load(fn_lbub);

        if dirs == 1
            nCond = 25;
        elseif dirs == 2
            nCond = 50;
        end

        n2 = length(goodfit_ind);
        if n2>1;
            r_blank = zeros(n2);
            for ipix1 = 1:n2
                for ipix2 = 1:n2
                    r_blank(ipix1,ipix2)= triu2vec(corrcoef(resp_dFoverF(goodfit_ind(:,ipix1),:),resp_dFoverF(goodfit_ind(:,ipix2),:)));
                end
            end 
        else
            r_blank_avg = NaN;
        end
        fn_out = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_bouton_corr_blank.mat']);
        save(fn_out, 'r_blank_avg');
    end
end
