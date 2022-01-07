clear all;
close all;
plotFits = 0;
ds = 'CrossOriRandDirRandPhase_ExptList';
eval(ds)
rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'RandDirRandPhase';
area_list = ['V1']; % 'LM'; 'AL'; 'PM'; 'RL'];
driver = 'SCN';
if min(size(area_list)) == 1
    narea = 1;
else
    narea =length(area_list);
end

for iarea = 1:narea
    area = area_list(iarea,:);
    fprintf([area ' ' driver '\n'])
    stim_OSI_all = [];
    plaid_OSI_all = [];
    stim_DSI_all = [];
    plaid_DSI_all = [];
    Zc_all = [];
    Zp_all = [];
    k_all = [];
%     R1_all = [];
%     R2_all = [];
    stim_SI_all = [];
    plaid_SI_all = [];
    totCells = 0;
    totExp = 0;
    expUse = [];
    resp_ind_all_dir = [];
    resp_ind_all_phase = [];
    resp_ind_dir_all = [];
    resp_ind_plaid_all = [];
    amp_all = [];
    amp_shuf_all = [];
    phase_SI_all = [];
    phase_MI_all = [];
    phase_MI_max_all = [];
    b_all = [];
    anova_all = [];
    mouse_list = [];
    respCellN = [];
    z_all = [];
    mouse_ind = [];
    for iexp = 1:nexp
        if sum(strcmp(expt(iexp).img_loc,area)) & strcmp(expt(iexp).driver,driver) & expt(iexp).con == 0.5
            mouse = expt(iexp).mouse;
            mouse_list = strvcat(mouse_list, mouse);
            date = expt(iexp).date;
            ImgFolder = expt(iexp).coFolder;
            z_all = [z_all; expt(iexp).z];
            nrun = length(ImgFolder);
            run_str = catRunName(cell2mat(ImgFolder), nrun);
            totExp = totExp+1;
            fprintf([mouse ' ' date '\n'])

            % load data

            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirAnalysis.mat']), 'Zc', 'Zp','k1_dir', 'stim_SI', 'stim_OSI', 'stim_DSI', 'plaid_OSI', 'plaid_DSI', 'plaid_SI', 'nCells');
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
                k1_dir = k1_dir(ind);
                plaid_SI = plaid_SI(ind);
                stim_SI = stim_SI(ind);
                h_resp = h_resp(ind,:,:);
                nCells = length(ind);
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

            resp_ind = find(sum(sum(h_resp,2),3));
            resp_ind_dir = find(sum(h_resp(:,:,1),2));
            resp_ind_plaid = find(sum(h_resp(:,:,2),2));

            resp_ind_all_dir = [resp_ind_all_dir resp_ind'+totCells];
            resp_ind_dir_all = [resp_ind_dir_all resp_ind_dir'+totCells];
            resp_ind_plaid_all = [resp_ind_plaid_all resp_ind_plaid'+totCells];

            mouse_ind = [mouse_ind; repmat(mouse, size(Zc'))];
            ImgFolder = expt(iexp).copFolder;
            nrun = length(ImgFolder);
            if ~isempty(ImgFolder)
                run_str = catRunName(cell2mat(ImgFolder), nrun);
                load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
                load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
                load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits.mat']))
                load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']))
                resp_ind_temp = unique([resptest_ind; respmask_ind]);
                resp_ind_all_phase = [resp_ind_all_phase resp_ind_temp'+totCells];
                if isempty(find(trN<4)) & length(resp_ind_all_phase)>4
                    test_resp = mean(resp_cell{1,end,1,1},2);
                    test_resp_rect = test_resp;
                    test_resp_rect(find(test_resp<0)) = 0;
                    mask_resp = mean(resp_cell{end,1,1,1},2);
                    mask_resp_rect = mask_resp;
                    mask_resp_rect(find(mask_resp<0)) = 0;
                    plaid_resp = mean(resp_cell{end,end,1,1},2);
                    plaid_resp_rect = plaid_resp;
                    plaid_resp_rect(find(plaid_resp<0)) = 0;
                    if length(expt(iexp).img_loc)>1
                        amp_hat_all = amp_hat_all(ind);
                        amp_shuf_all = amp_shuf_all(ind);
                        b_all = b_hat_all(ind);
                        test_resp_rect = test_resp_rect(ind);
                        mask_resp_rect = mask_resp_rect(ind);
                        plaid_resp_rect = plaid_resp_rect(ind);
                        p_anova_all = p_anova_all(ind);
                    end
                    anova_all = [anova_all p_anova_all'];
                    amp_all = [amp_all amp_hat_all'];
                    amp_shuf_all = [amp_shuf_all amp_hat_shuf'];
                    b_all = [b_all b_hat_all'];
                    phase_SI = abs(test_resp_rect-mask_resp_rect)./(test_resp_rect+mask_resp_rect);
                    phase_SI_all = [phase_SI_all phase_SI'];
                    phase_MI = (plaid_resp_rect-(test_resp_rect+mask_resp_rect))./(plaid_resp_rect+test_resp_rect+mask_resp_rect);
                    phase_MI_all = [phase_MI_all phase_MI'];
                    phase_MI_max = (plaid_resp_rect-max([test_resp_rect mask_resp_rect],[],2))./(plaid_resp_rect+max([test_resp_rect mask_resp_rect],[],2));
                    phase_MI_max_all = [phase_MI_max_all phase_MI_max'];
                    expUse(totExp) = 1;
                    respCellN = [respCellN length(resp_ind_temp')];
                else
                    anova_all = [anova_all nan(nCells,1)'];
                    amp_all = [amp_all nan(nCells,1)'];
                    amp_shuf_all = [amp_shuf_all nan(nCells,1)'];
                    b_all = [b_all nan(nCells,1)'];
                    phase_SI_all = [phase_SI_all nan(nCells,1)'];
                    phase_MI_all = [phase_MI_all nan(nCells,1)'];
                    phase_MI_max_all = [phase_MI_max_all nan(nCells,1)'];
                    expUse(totExp) = 0;
                    respCellN = [respCellN 0];
                end
            else
                anova_all = [anova_all nan(nCells,1)'];
                amp_all = [amp_all nan(nCells,1)'];
                amp_shuf_all = [amp_shuf_all nan(nCells,1)'];
                b_all = [b_all nan(nCells,1)'];
                phase_SI_all = [phase_SI_all nan(nCells,1)'];
                phase_MI_all = [phase_MI_all nan(nCells,1)'];
                phase_MI_max_all = [phase_MI_max_all nan(nCells,1)'];
                expUse(totExp) = 0;
                respCellN = [respCellN 0];
            end
            
            totCells = totCells+nCells;
            
            

            if plotFits & expUse(totExp)
                [n n2] = subplotn(length(resp_ind_temp));
                figure(totExp+1)
                for i = 1:length(resp_ind_temp)
                   subplot(n,n2,i)
                   errorbar(maskPhas, SI_all_avg(i,:,1),SI_all_avg(i,:,2),'o')
                   hold on
                   plot(1:360,yfit_all(i,:))
                   ylim([-1 1])
                   title(num2str(i))
                end
                suptitle([mouse ' ' date])
                print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_allCellsFitsPlusAvg.pdf']),'-dpdf','-bestfit')

                load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
                figure(1)
                subplot(3,4,totExp)
                imagesc(data_avg)
                title([mouse ' ' date])
            end
        end
    end
    mouse_list(find(expUse==0),:) = [];
    z_all(find(expUse==0),:) = [];
    save(fullfile(summaryDir,[svName '_Summary_' area '_' driver '.mat']),'mouse_ind','stim_SI_all','stim_OSI_all','plaid_OSI_all','stim_DSI_all','plaid_DSI_all','Zc_all','Zp_all','plaid_SI_all','resp_ind_all_dir','resp_ind_all_phase','resp_ind_dir_all','resp_ind_plaid_all', 'k_all','amp_all', 'amp_shuf_all', 'b_all', 'phase_SI_all', 'phase_MI_all', 'phase_MI_max_all','anova_all', 'mouse_list', 'z_all')
end
%%

    figure;
    subplot(2,2,1)
    cdfplot(stim_OSI_all(resp_ind_all_dir))
    hold on
    cdfplot(plaid_OSI_all(resp_ind_all_dir))
    xlabel('OSI')
    legend({'Stim','Plaid'},'Location','southeast')
    title('')
    subplot(2,2,2)
    cdfplot(stim_DSI_all(resp_ind_all_dir))
    hold on
    cdfplot(plaid_DSI_all(resp_ind_all_dir))
    xlabel('DSI')
    title('')
    legend({'Stim','Plaid'},'Location','southeast')
    subplot(2,2,3)
    cdfplot(Zc_all(resp_ind_all_dir))
    hold on
    cdfplot(Zp_all(resp_ind_all_dir))
    xlabel('Zc/Zp')
    xlim([-2 10])
    title('')
    legend({'Zc','Zp'},'Location','southeast')
    subplot(2,2,4)
    cdfplot(plaid_SI_all(resp_ind_all_dir))
    hold on
    cdfplot(plaid_SI_all(intersect(resp_ind_all_dir,find(stim_OSI_all<0.5))))
    cdfplot(plaid_SI_all(intersect(resp_ind_all_dir,find(stim_DSI_all<0.5))))
    xlabel('Suppression Index')
    title('')
    legend({'All','stim OSI<0.5', 'stim DSI<0.5'},'Location','southeast')
    suptitle({[area '- n = ' num2str(size(mouse_list,1)) ' expts; ' num2str(size(unique(mouse_list,'rows'),1)) ' mice'], ['All responsive cells- n = ' num2str(length(resp_ind_all_dir))]})
    print(fullfile(summaryDir, [svName '_OSI-DSI-Zc-Zp-SI_Summary_' area '_' driver '.pdf']),'-dpdf', '-fillpage')       

    figure;
    subplot(2,2,1)
    scatter(Zc_all(resp_ind_all_dir),Zp_all(resp_ind_all_dir))
    xlabel('Zc')
    ylabel('Zp')
    xlim([-5 10])
    ylim([-5 10])
    hold on
    plotZcZpBorders
    axis square
    subplot(2,2,2)
    cdfplot(Zc_all(resp_ind_all_dir))
    hold on
    cdfplot(Zp_all(resp_ind_all_dir))
    xlabel('Zc/Zp')
    xlim([-5 10])
    title('')
    suptitle(['All responsive cells- n = ' num2str(length(resp_ind_all_dir))])
    print(fullfile(summaryDir, [svName '_Zc-Zp_Scatter_' area '_' driver '.pdf']),'-dpdf', '-fillpage') 

    figure;
    subplot(3,2,1)
    cdfplot(stim_OSI_all(intersect(resp_ind_all_dir,find(plaid_SI_all<0))))
    hold on
    cdfplot(stim_OSI_all(intersect(resp_ind_all_dir,find(plaid_SI_all>0))))
    xlabel('stim OSI')
    legend({'SI<0', 'SI>0'},'Location','northwest')
    title('')
    subplot(3,2,2)
    cdfplot(plaid_OSI_all(intersect(resp_ind_all_dir,find(plaid_SI_all<0))))
    hold on
    cdfplot(plaid_OSI_all(intersect(resp_ind_all_dir,find(plaid_SI_all>0))))
    xlabel('plaid OSI')
    title('')
    subplot(3,2,3)
    cdfplot(stim_DSI_all(intersect(resp_ind_all_dir,find(plaid_SI_all<0))))
    hold on
    cdfplot(stim_DSI_all(intersect(resp_ind_all_dir,find(plaid_SI_all>0))))
    xlabel('stim DSI')
    title('')
    subplot(3,2,4)
    cdfplot(plaid_DSI_all(intersect(resp_ind_all_dir,find(plaid_SI_all<0))))
    hold on
    cdfplot(plaid_DSI_all(intersect(resp_ind_all_dir,find(plaid_SI_all>0))))
    xlabel('plaid DSI')
    title('')
    subplot(3,2,5)
    cdfplot(Zc_all(intersect(resp_ind_all_dir,find(plaid_SI_all<0))))
    hold on
    cdfplot(Zc_all(intersect(resp_ind_all_dir,find(plaid_SI_all>0))))
    xlabel('Zc')
    xlim([-2 10])
    title('')
    subplot(3,2,6)
    cdfplot(Zp_all(intersect(resp_ind_all_dir,find(plaid_SI_all<0))))
    hold on
    cdfplot(Zp_all(intersect(resp_ind_all_dir,find(plaid_SI_all>0))))
    xlabel('Zp')
    xlim([-2 10])
    title('')
    suptitle({'High vs low Suppression index',['All responsive cells- n = ' num2str(length(resp_ind_all_dir))]})
    print(fullfile(summaryDir, [svName '_highVlowSI' area '_' driver '.pdf']),'-dpdf', '-fillpage') 

    figure;
    subplot(2,2,1)
    cdfplot(plaid_SI_all(intersect(resp_ind_all_dir,find(stim_OSI_all<0.5))))
    hold on
    cdfplot(plaid_SI_all(intersect(resp_ind_all_dir,find(stim_OSI_all>0.5))))
    xlabel('Suppression index')
    legend({'OSI<0.5', 'OSI>0.5'},'Location','northwest')
    title('')
    subplot(2,2,2)
    cdfplot(stim_DSI_all(intersect(resp_ind_all_dir,find(stim_OSI_all<0.5))))
    hold on
    cdfplot(stim_DSI_all(intersect(resp_ind_all_dir,find(stim_OSI_all>0.5))))
    xlabel('stim DSI')
    title('')
    subplot(2,2,3)
    cdfplot(Zc_all(intersect(resp_ind_all_dir,find(stim_OSI_all<0.5))))
    hold on
    cdfplot(Zc_all(intersect(resp_ind_all_dir,find(stim_OSI_all>0.5))))
    xlabel('Zc')
    xlim([-2 10])
    title('')
    subplot(2,2,4)
    cdfplot(Zp_all(intersect(resp_ind_all_dir,find(stim_OSI_all<0.5))))
    hold on
    cdfplot(Zp_all(intersect(resp_ind_all_dir,find(stim_OSI_all>0.5))))
    xlabel('Zp')
    xlim([-2 10])
    title('')
    suptitle({'High vs low OSI', ['All responsive cells- n = ' num2str(length(resp_ind_all_dir))]})
    print(fullfile(summaryDir, ['randDir_highVlowOSI_' area '_' driver '.pdf']),'-dpdf', '-fillpage') 

    figure;
    subplot(2,2,1)
    cdfplot(plaid_SI_all(intersect(resp_ind_all_dir,find(stim_DSI_all<0.5))))
    hold on
    cdfplot(plaid_SI_all(intersect(resp_ind_all_dir,find(stim_DSI_all>0.5))))
    xlabel('Suppression index')
    legend({'DSI<0.5', 'DSI>0.5'},'Location','northwest')
    title('')
    subplot(2,2,2)
    cdfplot(plaid_DSI_all(intersect(resp_ind_all_dir,find(stim_DSI_all<0.5))))
    hold on
    cdfplot(plaid_DSI_all(intersect(resp_ind_all_dir,find(stim_DSI_all>0.5))))
    xlabel('plaid DSI')
    title('')
    subplot(2,2,3)
    cdfplot(Zc_all(intersect(resp_ind_all_dir,find(stim_DSI_all<0.5))))
    hold on
    cdfplot(Zc_all(intersect(resp_ind_all_dir,find(stim_DSI_all>0.5))))
    xlabel('Zc')
    xlim([-2 10])
    title('')
    subplot(2,2,4)
    cdfplot(Zp_all(intersect(resp_ind_all_dir,find(stim_DSI_all<0.5))))
    hold on
    cdfplot(Zp_all(intersect(resp_ind_all_dir,find(stim_DSI_all>0.5))))
    xlabel('Zp')
    xlim([-2 10])
    title('')
    suptitle({'High vs low DSI', ['All responsive cells- n = ' num2str(length(resp_ind_all_dir))]})
    print(fullfile(summaryDir, [svName '_highVlowDSI_' area '_' driver '.pdf']),'-dpdf', '-fillpage') 
    
    Zp_use = intersect(resp_ind_all_dir, intersect(find(Zp_all>1.28), find(Zp_all-Zc_all>1.28)));
    Zc_use = intersect(resp_ind_all_dir, intersect(find(Zc_all>1.28), find(Zc_all-Zp_all>1.28)));
    figure;
    subplot(2,2,1)
    cdfplot(stim_OSI_all(Zc_use))
    hold on
    cdfplot(stim_OSI_all(Zp_use))
    xlabel('OSI')
    xlim([0 1])
    title('')
    subplot(2,2,2)
    cdfplot(stim_DSI_all(Zc_use))
    hold on
    cdfplot(stim_DSI_all(Zp_use))
    xlabel('DSI')
    xlim([0 1])
    title('')
    subplot(2,2,3)
    cdfplot(k_all(Zc_use))
    hold on
    cdfplot(k_all(Zp_use))
    xlabel('Kappa')
    xlim([0 30])
    title('')
    suptitle(['Tuning of Zc (n= ' num2str(length(Zc_use)) '); Zp (n = ' num2str(length(Zp_use)) ')'])
    print(fullfile(summaryDir, [svName '_ZcZp_Tuning_' area '_' driver '.pdf']),'-dpdf', '-fillpage')
    end