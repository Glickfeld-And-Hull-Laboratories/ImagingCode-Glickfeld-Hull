clear all
areas = ['PM'; 'LM'; 'AL'];
for iArea = 1:3
    P = 2;
    matrix = 'SF5xTF5';
    inj = 'V1';
    image = areas(iArea,:);
    nON = 12;
    nOFF = 12;
    nPlanes = 1;

    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);
    
    npmask_dF_all = zeros(25, nexp);
    nphalo_dF_all = zeros(25, nexp);
    for iexp = 1:nexp;

        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        userun = exp_list.runs_mat{iexp};
        count_protocol = exp_list.prot_mat{iexp};
        run = exp_list.run_mat{iexp};
        blanks = exp_list.blanks_mat{iexp};
        dirs = exp_list.dir_mat{iexp};

        if dirs ==1
            nCond = 25;
        elseif dirs ==2
            nCond = 50;
        end
        
        base = 'G:\users\lindsey\analysisLG\active mice';
        outDir = fullfile(base, mouse,date);

        fn_stack= fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_stack_dF_all.tif']);
        stack_dF = readtiff(fn_stack);

         if dirs == 2
            add = 1;
            stack_dF_temp = zeros(size(stack_dF));
            for iCond = 1:nCond/2
            stack_dF_temp(:,:,iCond) = mean(stack_dF(:,:,add:add+1),3);
            add = add+2;
            end
            stack_dF_temp(:,:,end) = stack_dF(:,:,end);
            stack_dF = stack_dF_temp;
         end
        
        i = [];
        j=[];
        local_max = [];
        fn_local = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_local_max.mat']);
        load(fn_local);
        
        bouton_map = zeros(size(local_max));
        for ipix = 1:n_pix
            sub_y = i(ipix)-1:i(ipix)+1; 
            sub_x = j(ipix)-1:j(ipix)+1;
            bouton_map(sub_y, sub_x) = ipix;
        end

        siz = size(bouton_map); 
        npmasks = zeros(siz(1),siz(2),n_pix);
        mask_dil = zeros(siz(1),siz(2),n_pix);
        for ipix = 1:n_pix
            mask_dil(i(ipix)-2:i(ipix)+2, j(ipix)-2:j(ipix)+2,ipix) = 1;
        end
        all_mask_dil = sum(mask_dil,3);

        npmask = ones(size(all_mask_dil));
        npmask(find(all_mask_dil>0))=0;
        npmask(1:5,:)=0;
        npmask(236:240,:)=0;
        npmask(:,1:5)=0;
        npmask(:,252:256)=0;
        
        npmask_temp_all = zeros(siz(1), siz(2), n_pix);
        for ipix = 1:n_pix
            np_plus = zeros(siz(1),siz(2));
            np_plus(i(ipix)-4:i(ipix)+4, j(ipix)-4:j(ipix)+4) = 1;
            npmask_temp = np_plus-mask_dil(:,:,ipix);
            npmask_temp(find(all_mask_dil>0)) = 0;
            npmask_temp(1:5,:)=0;
            npmask_temp(236:240,:)=0;
            npmask_temp(:,1:5)=0;
            npmask_temp(:,252:256)=0;
            npmask_temp_all(:,:,ipix) = npmask_temp;
        end
        
        nphalo = sum(npmask_temp_all,3);
        nphalo(find(nphalo>0))=1;
        
        npmask_dF = zeros(1,25);
        nphalo_dF = zeros(1,25);
        for iCond = 1:25
            stack_dF_temp = stack_dF(:,:,iCond);
            npmask_dF(:,iCond) = mean(stack_dF_temp(find(npmask>0)));
            nphalo_dF(:,iCond) = mean(stack_dF_temp(find(nphalo>0)));
        end
        
        npmask_dF_all(:,iexp) = npmask_dF';
        nphalo_dF_all(:,iexp) = nphalo_dF';
    
        fn_neuropil = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_neuropil.mat']);
        save(fn_neuropil, 'npmask', 'nphalo', 'npmask_dF', 'nphalo_dF');
    end
    fn_out = fullfile('G:\users\lindsey\analysisLG\experiments', [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_neuropil.mat']);
    save(fn_out, 'npmask_dF_all', 'nphalo_dF_all');
end


