clear all;
close all;
clc
ds = 'CrossOriRandDirFF_ExptList';
eval(ds)
frame_rate = 15;
nexp = size(expt,2);
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'RandDirSummary');
svName = 'randDirFF';
if ~exist(summaryDir)
    mkdir(summaryDir)
end

doPlot = 0;
area_list = ['V1']; % 'LM'; 'AL'; 'PM'; 'RL'];
driver = 'SLC';
doRedCells =0;
SF = 0.05;
con = 0.5;
sz = 1000;
doSFSave = 1;
doConSave = 0;
doSzSave = 1;
if min(size(area_list)) == 1
    narea = 1;
else
    narea =length(area_list);
end

for iarea = narea
    area = area_list(iarea,:);
    fprintf([area ' ' driver '\n'])
    stim_OSI_all = [];
    plaid_OSI_all = [];
    stim_DSI_all = [];
    plaid_DSI_all = [];
    Zc_all = [];
    Zp_all = [];
    k_all = [];
    avg_resp_dir_align_all = [];
    avg_resp_plaid_align_all = [];
    norm_resp_plaid_align_all = [];
%     R1_all = [];
%     R2_all = [];
    stim_SI_all = [];
    plaid_SI_all = [];
    totCells = 0;
    resp_ind_all = [];
    resp_ind_dir_all = [];
    resp_ind_plaid_all = [];
    f1_all = [];
    f2_all = [];
    f2overf1_all = [];
    mouse_list = [];
    component_all = [];
    pattern_all = [];
    avg_resp_dir_all = [];
    h_resp_all = [];
    red_cells_all = [];
    z_all = [];
    for iexp = 1:nexp
        if sum(strcmp(expt(iexp).img_loc,area)) & strcmp(expt(iexp).driver,driver) ... 
                & (~isfield(expt,'SF') || (isfield(expt,'SF') & expt(iexp).SF == SF)) ...
                & (~isfield(expt,'con') || (isfield(expt,'con') & expt(iexp).con == con)) ...
                & (~isfield(expt,'size') || (isfield(expt,'size') & expt(iexp).size == sz))
            mouse = expt(iexp).mouse;
            mouse_list = strvcat(mouse_list, mouse);
            date = expt(iexp).date;
            ImgFolder = expt(iexp).coFolder;
            time = expt(iexp).coTime;
            nrun = length(ImgFolder);
            run_str = catRunName(cell2mat(ImgFolder), nrun);

            fprintf([mouse ' ' date '\n'])

            % load data
            if doRedCells
                load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']),'red_cells')
            else
                red_cells = [];
            end
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirAnalysis.mat']), 'Zc', 'Zp','k1_dir', 'stim_SI', 'stim_OSI', 'stim_DSI', 'plaid_OSI', 'plaid_DSI', 'plaid_SI', 'nCells');
            int = unique(diff(stimDirs));
            component = avg_resp_dir(:,:,1,1,1)+circshift(avg_resp_dir(:,:,1,1,1),-90./int,2);
            pattern = circshift(avg_resp_dir(:,:,1,1,1),-45./int,2);
            if length(expt(iexp).img_loc)>1
                i = find(strcmp(expt(iexp).img_loc,area));
                load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_splitImage.mat']))
                ind = find(maskCat==i);
                stim_OSI = stim_OSI(ind);
                plaid_OSI = plaid_OSI(ind);
                stim_DSI = stim_DSI(ind);
                plaid_DSI = plaid_DSI(ind);
                Zc = Zc(ind);
                Zp = Zp(ind);
                avg_resp_dir = avg_resp_dir(ind,:,:,:);
                k1_dir = k1_dir(ind);
                plaid_SI = plaid_SI(ind);
                stim_SI = stim_SI(ind);
                h_resp = h_resp(ind,:,:);
                nCells = length(ind);
                component = component(ind,:,:,:);
                pattern = pattern(ind,:,:,:);
            end
            fprintf(['n = ' num2str(nCells) '\n'])
            stim_SI_all = [stim_SI_all stim_SI];
            stim_OSI_all = [stim_OSI_all stim_OSI];
            plaid_OSI_all = [plaid_OSI_all plaid_OSI];
            stim_DSI_all = [stim_DSI_all stim_DSI];
            plaid_DSI_all = [plaid_DSI_all plaid_DSI];
            Zc_all = [Zc_all Zc];
            Zp_all = [Zp_all Zp];
            k_all = [k_all k1_dir];
            plaid_SI_all = [plaid_SI_all plaid_SI];
            avg_resp_dir_all = cat(1,avg_resp_dir_all,avg_resp_dir);
            avg_resp_dir_align = zeros(nCells,size(avg_resp_dir,2));
            avg_resp_plaid_align = zeros(nCells,size(avg_resp_dir,2));
            norm_resp_plaid_align = zeros(nCells,size(avg_resp_dir,2));
            component_all = [component_all; component];
            pattern_all = [pattern_all; pattern];
            h_resp_all = cat(1,h_resp_all, h_resp);
            for iC = 1:size(avg_resp_dir,1)
                [max_val max_ind] = max(avg_resp_dir(iC,:,1,1),[],2);
                avg_resp_dir_align(iC,:) = circshift(avg_resp_dir(iC,:,1,1),9-max_ind);
                avg_resp_plaid_align(iC,:) = circshift(avg_resp_dir(iC,:,2,1),9-max_ind);
                norm_resp_plaid_align(iC,:) = circshift(avg_resp_dir(iC,:,2,1),9-max_ind)./max_val;
            end
            avg_resp_dir_align_all = cat(1,avg_resp_dir_align_all,avg_resp_dir_align);
            avg_resp_plaid_align_all = cat(1,avg_resp_plaid_align_all,avg_resp_plaid_align);
            norm_resp_plaid_align_all = cat(1,norm_resp_plaid_align_all,norm_resp_plaid_align);
            resp_ind = find(sum(sum(h_resp,2),3));
            resp_ind_dir = find(sum(h_resp(:,:,1),2));
            resp_ind_plaid = find(sum(h_resp(:,:,2),2));

            resp_ind_all = [resp_ind_all resp_ind'+totCells];
            resp_ind_dir_all = [resp_ind_dir_all resp_ind_dir'+totCells];
            resp_ind_plaid_all = [resp_ind_plaid_all resp_ind_plaid'+totCells];
            
            red_cells_all = [red_cells_all red_cells'+totCells];
            z_all = [z_all expt(iexp).z.*ones(1,totCells)];

            if isfield(expt,'prFolder')
                if ~isempty(expt(iexp).prFolder)
                    ImgFolder = expt(iexp).prFolder;
                    nrun = length(ImgFolder);
                    run_str = catRunName(cell2mat(ImgFolder), nrun);
                    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_f1f2.mat']))
                    if length(expt(iexp).img_loc)>1
                        f1 = f1(ind);
                        f2 = f2(ind);
                        f2overf1 = f2overf1(ind);
                    end
                    f1_all = [f1_all f1];
                    f2_all = [f2_all f2];
                    f2overf1_all = [f2overf1_all f2overf1];
                else
                f1_all = [f1_all nan(size(stim_OSI))];
                f2_all = [f2_all nan(size(stim_OSI))];
                f2overf1_all = [f2overf1_all nan(size(stim_OSI))];
                end
            else
                f1_all = [f1_all nan(size(stim_OSI))];
                f2_all = [f2_all nan(size(stim_OSI))];
                f2overf1_all = [f2overf1_all nan(size(stim_OSI))];
            end

            totCells = totCells+nCells;

        end
    end
    if doRedCells
        driver_str= 'noDriver';
    else
        driver_str = driver;
    end
    sf_str = num2str(SF);
    sf_str = ['pt' sf_str(3:end)];
    con_str = num2str(con);
    con_str = ['pt' con_str(3:end)];
    sz_str = num2str(sz);
    if doSFSave & doSzSave
        save(fullfile(summaryDir,[svName '_Summary_' area '_' driver_str '_SF' sf_str '_Sz' sz_str '.mat']),'z_all','red_cells_all','stim_SI_all','stim_OSI_all','plaid_OSI_all','stim_DSI_all','plaid_DSI_all','Zc_all','Zp_all','plaid_SI_all','resp_ind_all','resp_ind_dir_all','resp_ind_plaid_all', 'f1_all','f2_all','f2overf1_all','k_all','avg_resp_dir_all','avg_resp_dir_align_all','avg_resp_plaid_align_all','norm_resp_plaid_align_all','component_all','pattern_all','stimDirs','h_resp_all','mouse_list')
    elseif doSFSave
        save(fullfile(summaryDir,[svName '_Summary_' area '_' driver_str '_SF' sf_str '.mat']),'z_all','red_cells_all','stim_SI_all','stim_OSI_all','plaid_OSI_all','stim_DSI_all','plaid_DSI_all','Zc_all','Zp_all','plaid_SI_all','resp_ind_all','resp_ind_dir_all','resp_ind_plaid_all', 'f1_all','f2_all','f2overf1_all','k_all','avg_resp_dir_all','avg_resp_dir_align_all','avg_resp_plaid_align_all','norm_resp_plaid_align_all','component_all','pattern_all','stimDirs','h_resp_all','mouse_list')
    elseif doConSave
        save(fullfile(summaryDir,[svName '_Summary_' area '_' driver_str '_Con' con_str '.mat']),'z_all','red_cells_all','stim_SI_all','stim_OSI_all','plaid_OSI_all','stim_DSI_all','plaid_DSI_all','Zc_all','Zp_all','plaid_SI_all','resp_ind_all','resp_ind_dir_all','resp_ind_plaid_all', 'f1_all','f2_all','f2overf1_all','k_all','avg_resp_dir_all','avg_resp_dir_align_all','avg_resp_plaid_align_all','norm_resp_plaid_align_all','component_all','pattern_all','stimDirs','h_resp_all','mouse_list')
    else
        save(fullfile(summaryDir,[svName '_Summary_' area '_' driver_str '.mat']),'z_all','red_cells_all','stim_SI_all','stim_OSI_all','plaid_OSI_all','stim_DSI_all','plaid_DSI_all','Zc_all','Zp_all','plaid_SI_all','resp_ind_all','resp_ind_dir_all','resp_ind_plaid_all', 'f1_all','f2_all','f2overf1_all','k_all','avg_resp_dir_all','avg_resp_dir_align_all','avg_resp_plaid_align_all','norm_resp_plaid_align_all','component_all','pattern_all','stimDirs','h_resp_all','mouse_list')        
    end
%%
if doPlot
    resp_ind_all= setdiff(resp_ind_all, red_cells_all);
    figure;
    subplot(2,2,1)
    cdfplot(stim_OSI_all(resp_ind_all))
    hold on
    cdfplot(plaid_OSI_all(resp_ind_all))
    xlabel('OSI')
    legend({'Stim','Plaid'},'Location','southeast')
    title('')
    subplot(2,2,2)
    cdfplot(stim_DSI_all(resp_ind_all))
    hold on
    cdfplot(plaid_DSI_all(resp_ind_all))
    xlabel('DSI')
    title('')
    legend({'Stim','Plaid'},'Location','southeast')
    subplot(2,2,3)
    cdfplot(Zc_all(resp_ind_all))
    hold on
    cdfplot(Zp_all(resp_ind_all))
    xlabel('Zc/Zp')
    xlim([-2 10])
    title('')
    legend({'Zc','Zp'},'Location','southeast')
    subplot(2,2,4)
    cdfplot(plaid_SI_all(resp_ind_all))
    hold on
    cdfplot(plaid_SI_all(intersect(resp_ind_all,find(stim_OSI_all<0.5))))
    cdfplot(plaid_SI_all(intersect(resp_ind_all,find(stim_DSI_all<0.5))))
    xlabel('Suppression Index')
    title('')
    legend({'All','stim OSI<0.5', 'stim DSI<0.5'},'Location','southeast')
    suptitle({[area '- n = ' num2str(size(mouse_list,1)) ' expts; ' num2str(size(unique(mouse_list,'rows'),1)) ' mice'], ['All responsive cells- n = ' num2str(length(resp_ind_all))]})
    print(fullfile(summaryDir, [svName '_OSI-DSI-Zc-Zp-SI_Summary_' area '_' driver_str '_SF' sf_str '_Sz' sz_str '_Con' con_str '.pdf']),'-dpdf', '-fillpage')       

    figure;
    subplot(2,2,1)
    scatter(Zc_all(resp_ind_all),Zp_all(resp_ind_all))
    xlabel('Zc')
    ylabel('Zp')
    xlim([-5 10])
    ylim([-5 10])
    hold on
    plotZcZpBorders
    axis square
    subplot(2,2,2)
    cdfplot(Zc_all(resp_ind_all))
    hold on
    cdfplot(Zp_all(resp_ind_all))
    xlabel('Zc/Zp')
    xlim([-5 10])
    title('')
    sgtitle(['All responsive cells- n = ' num2str(length(resp_ind_all))])
    print(fullfile(summaryDir, [svName '_Zc-Zp_Scatter_' area '_' driver_str '_SF' sf_str '_Sz' sz_str '_Con' con_str '.pdf']),'-dpdf', '-fillpage') 
    
    if length(red_cells_all)
        figure;
        subplot(2,2,1)
        cdfplot(stim_OSI_all(red_cells_all))
        hold on
        cdfplot(plaid_OSI_all(red_cells_all))
        xlabel('OSI')
        legend({'Stim','Plaid'},'Location','southeast')
        title('')
        subplot(2,2,2)
        cdfplot(stim_DSI_all(red_cells_all))
        hold on
        cdfplot(plaid_DSI_all(red_cells_all))
        xlabel('DSI')
        title('')
        legend({'Stim','Plaid'},'Location','southeast')
        subplot(2,2,3)
        cdfplot(Zc_all(red_cells_all))
        hold on
        cdfplot(Zp_all(red_cells_all))
        xlabel('Zc/Zp')
        xlim([-2 10])
        title('')
        legend({'Zc','Zp'},'Location','southeast')
        subplot(2,2,4)
        cdfplot(plaid_SI_all(red_cells_all))
        hold on
        cdfplot(plaid_SI_all(intersect(red_cells_all,find(stim_OSI_all<0.5))))
        cdfplot(plaid_SI_all(intersect(red_cells_all,find(stim_DSI_all<0.5))))
        xlabel('Suppression Index')
        title('')
        legend({'All','stim OSI<0.5', 'stim DSI<0.5'},'Location','southeast')
        suptitle({[area '- n = ' num2str(size(mouse_list,1)) ' expts; ' num2str(size(unique(mouse_list,'rows'),1)) ' mice'], ['All ' driver ' cells- n = ' num2str(length(red_cells_all))]})
        print(fullfile(summaryDir, [svName '_OSI-DSI-Zc-Zp-SI_Summary_' area '_' driver '_SF' sf_str '_Sz' sz_str '_Con' con_str '.pdf']),'-dpdf', '-fillpage')       

        figure;
        subplot(2,2,1)
        scatter(Zc_all(red_cells_all),Zp_all(red_cells_all))
        xlabel('Zc')
        ylabel('Zp')
        xlim([-5 10])
        ylim([-5 10])
        hold on
        plotZcZpBorders
        axis square
        subplot(2,2,2)
        cdfplot(Zc_all(red_cells_all))
        hold on
        cdfplot(Zp_all(red_cells_all))
        xlabel('Zc/Zp')
        xlim([-5 10])
        title('')
        sgtitle(['All ' driver ' cells- n = ' num2str(length(red_cells_all))])
        print(fullfile(summaryDir, [svName '_Zc-Zp_Scatter_' area '_' driver '_SF' sf_str '_Sz' sz_str '_Con' con_str '.pdf']),'-dpdf', '-fillpage') 

    end
end
%     figure;
%     subplot(3,2,1)
%     cdfplot(stim_OSI_all(intersect(resp_ind_all,find(plaid_SI_all<0))))
%     hold on
%     cdfplot(stim_OSI_all(intersect(resp_ind_all,find(plaid_SI_all>0))))
%     xlabel('stim OSI')
%     legend({'SI<0', 'SI>0'},'Location','northwest')
%     title('')
%     subplot(3,2,2)
%     cdfplot(plaid_OSI_all(intersect(resp_ind_all,find(plaid_SI_all<0))))
%     hold on
%     cdfplot(plaid_OSI_all(intersect(resp_ind_all,find(plaid_SI_all>0))))
%     xlabel('plaid OSI')
%     title('')
%     subplot(3,2,3)
%     cdfplot(stim_DSI_all(intersect(resp_ind_all,find(plaid_SI_all<0))))
%     hold on
%     cdfplot(stim_DSI_all(intersect(resp_ind_all,find(plaid_SI_all>0))))
%     xlabel('stim DSI')
%     title('')
%     subplot(3,2,4)
%     cdfplot(plaid_DSI_all(intersect(resp_ind_all,find(plaid_SI_all<0))))
%     hold on
%     cdfplot(plaid_DSI_all(intersect(resp_ind_all,find(plaid_SI_all>0))))
%     xlabel('plaid DSI')
%     title('')
%     subplot(3,2,5)
%     cdfplot(Zc_all(intersect(resp_ind_all,find(plaid_SI_all<0))))
%     hold on
%     cdfplot(Zc_all(intersect(resp_ind_all,find(plaid_SI_all>0))))
%     xlabel('Zc')
%     xlim([-2 10])
%     title('')
%     subplot(3,2,6)
%     cdfplot(Zp_all(intersect(resp_ind_all,find(plaid_SI_all<0))))
%     hold on
%     cdfplot(Zp_all(intersect(resp_ind_all,find(plaid_SI_all>0))))
%     xlabel('Zp')
%     xlim([-2 10])
%     title('')
%     sgtitle({'High vs low Suppression index',['All responsive cells- n = ' num2str(length(resp_ind_all))]})
%     print(fullfile(summaryDir, [svName '_highVlowSI' area '_' driver_str '_SF' sf_str '_Sz' sz_str '_Con' con_str '.pdf']),'-dpdf', '-fillpage') 
% 
%     figure;
%     subplot(2,2,1)
%     cdfplot(plaid_SI_all(intersect(resp_ind_all,find(stim_OSI_all<0.5))))
%     hold on
%     cdfplot(plaid_SI_all(intersect(resp_ind_all,find(stim_OSI_all>0.5))))
%     xlabel('Suppression index')
%     legend({'OSI<0.5', 'OSI>0.5'},'Location','northwest')
%     title('')
%     subplot(2,2,2)
%     cdfplot(stim_DSI_all(intersect(resp_ind_all,find(stim_OSI_all<0.5))))
%     hold on
%     cdfplot(stim_DSI_all(intersect(resp_ind_all,find(stim_OSI_all>0.5))))
%     xlabel('stim DSI')
%     title('')
%     subplot(2,2,3)
%     cdfplot(Zc_all(intersect(resp_ind_all,find(stim_OSI_all<0.5))))
%     hold on
%     cdfplot(Zc_all(intersect(resp_ind_all,find(stim_OSI_all>0.5))))
%     xlabel('Zc')
%     xlim([-2 10])
%     title('')
%     subplot(2,2,4)
%     cdfplot(Zp_all(intersect(resp_ind_all,find(stim_OSI_all<0.5))))
%     hold on
%     cdfplot(Zp_all(intersect(resp_ind_all,find(stim_OSI_all>0.5))))
%     xlabel('Zp')
%     xlim([-2 10])
%     title('')
%     sgtitle({'High vs low OSI', ['All responsive cells- n = ' num2str(length(resp_ind_all))]})
%     print(fullfile(summaryDir, [svName '_highVlowOSI_' area '_' driver_str '_SF' sf_str '_Sz' sz_str '_Con' con_str '.pdf']),'-dpdf', '-fillpage') 
% 
%     figure;
%     subplot(2,2,1)
%     cdfplot(plaid_SI_all(intersect(resp_ind_all,find(stim_DSI_all<0.5))))
%     hold on
%     cdfplot(plaid_SI_all(intersect(resp_ind_all,find(stim_DSI_all>0.5))))
%     xlabel('Suppression index')
%     legend({'DSI<0.5', 'DSI>0.5'},'Location','northwest')
%     title('')
%     subplot(2,2,2)
%     cdfplot(plaid_DSI_all(intersect(resp_ind_all,find(stim_DSI_all<0.5))))
%     hold on
%     cdfplot(plaid_DSI_all(intersect(resp_ind_all,find(stim_DSI_all>0.5))))
%     xlabel('plaid DSI')
%     title('')
%     subplot(2,2,3)
%     cdfplot(Zc_all(intersect(resp_ind_all,find(stim_DSI_all<0.5))))
%     hold on
%     cdfplot(Zc_all(intersect(resp_ind_all,find(stim_DSI_all>0.5))))
%     xlabel('Zc')
%     xlim([-2 10])
%     title('')
%     subplot(2,2,4)
%     cdfplot(Zp_all(intersect(resp_ind_all,find(stim_DSI_all<0.5))))
%     hold on
%     cdfplot(Zp_all(intersect(resp_ind_all,find(stim_DSI_all>0.5))))
%     xlabel('Zp')
%     xlim([-2 10])
%     title('')
%     sgtitle({'High vs low DSI', ['All responsive cells- n = ' num2str(length(resp_ind_all))]})
%     print(fullfile(summaryDir, [svName '_highVlowDSI_' area '_' driver_str '_SF' sf_str '_Sz' sz_str '_Con' con_str '.pdf']),'-dpdf', '-fillpage') 
%     
%     Zp_use = intersect(resp_ind_all, intersect(find(Zp_all>1.28), find(Zp_all-Zc_all>1.28)));
%     Zc_use = intersect(resp_ind_all, intersect(find(Zc_all>1.28), find(Zc_all-Zp_all>1.28)));
%     figure;
%     subplot(2,2,1)
%     cdfplot(stim_OSI_all(Zc_use))
%     hold on
%     cdfplot(stim_OSI_all(Zp_use))
%     xlabel('OSI')
%     xlim([0 1])
%     title('')
%     subplot(2,2,2)
%     cdfplot(stim_DSI_all(Zc_use))
%     hold on
%     cdfplot(stim_DSI_all(Zp_use))
%     xlabel('DSI')
%     xlim([0 1])
%     title('')
%     subplot(2,2,3)
%     cdfplot(k_all(Zc_use))
%     hold on
%     cdfplot(k_all(Zp_use))
%     xlabel('Kappa')
%     xlim([0 30])
%     title('')
%     sgtitle(['Tuning of Zc (n= ' num2str(length(Zc_use)) '); Zp (n = ' num2str(length(Zp_use)) ')'])
%     print(fullfile(summaryDir, [svName '_ZcZp_Tuning_' area '_' driver_str '_SF' sf_str '_Sz' sz_str '_Con' con_str '.pdf']),'-dpdf', '-fillpage')
end