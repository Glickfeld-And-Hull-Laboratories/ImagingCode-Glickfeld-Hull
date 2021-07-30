clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandDirTwoPhase_ExptList';
eval(ds)
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'RandDirSummary');
rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);
area_list = ['V1'];
driver = 'SLC';
narea =length(area_list);

for iarea = 1:narea
    area = area_list(iarea,:);
    fprintf([area '\n'])
    area_cell = [expt.img_loc];
    nexp_area = sum(strcmp(area_cell(1,:),area));
    dir_corr_all = nan(16,16,nexp_area);
    corr_dir_all = nan(16,16,nexp_area);
    stimplaid_dist_all = nan(16,16,nexp_area);
    stimstim_dist_all = nan(16,16,nexp_area);
    stimplaid_ang_all = nan(16,16,nexp_area);
    stimstim_ang_all = nan(16,16,nexp_area);
    Zc_all = nan(1,nexp_area);
    Zp_all = nan(1,nexp_area);
    u1_dir_all = [];
    avg_resp_dir_all = [];
    resp_norm_avg_all = zeros(16,16,nexp_area);
    pca_exp_stim_all = nan(100,nexp_area);
    pca_exp_plaid_all = nan(100,nexp_area);
    stim_pca_amp_all = nan(100,16,nexp_area);
    plaid_pca_amp_all = nan(100,16,nexp_area);
    stimplaid_dists = [];
    stimplaid_angs = [];
    stimplaid_corrs = [];
    exp_ind = [];
    start = 1;
    %%
    for iexp = 1:nexp
        if sum(strcmp(expt(iexp).img_loc,area)) & strcmp(expt(iexp).driver,driver)& expt(iexp).SF == 0.05
            exp_ind = [exp_ind iexp];
            mouse = expt(iexp).mouse;
            date = expt(iexp).date;
            ImgFolder = expt(iexp).coFolder;
            time = expt(iexp).coTime;
            nrun = length(ImgFolder);
            run_str = catRunName(cell2mat(ImgFolder), nrun);

            fprintf([mouse ' ' date '\n'])

            %% load data

            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
            
            if length(expt(iexp).img_loc)>1
                i = find(strcmp(expt(iexp).img_loc,area));
                load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_splitImage.mat']))
                ind = find(maskCat==i);
                npSub_tc = npSub_tc(:,ind);
            end
            nOn = unique(cell2mat(input.nStimOneFramesOn));
            nOff = unique(cell2mat(input.tItiWaitFrames));
            nCells = size(npSub_tc,2);
            data_trial = nan(nOn+nOff,nCells,nTrials);
            for iTrial = 1:nTrials-1
                data_trial(:,:,iTrial) = npSub_tc(cStimOn(iTrial)-frame_rate:cStimOn(iTrial)+nOn+nOff-frame_rate-1,:);
            end
            data_f = mean(data_trial(1:frame_rate,:,:),1);
            data_df = data_trial-data_f;
            data_dfof = data_df./data_f;

            sz = size(data_dfof);
            nCells = sz(2);
            nFrames = sz(1);
            nStimDir = length(stimDirs);
            tc_dir_dfof = zeros(nOn+nOff,nCells,nStimDir,2);
            for iDir = 1:nStimDir
                ind_stimdir = find(stimDir_all == stimDirs(iDir));
                ind_maskdir = find(maskDir_all == stimDirs(iDir));
                ind_diralone = [intersect(ind_stimdir, ind_stimAlone), intersect(ind_maskdir, ind_maskAlone)];
                ind_dirplaid = [intersect(ind_stimdir, ind_plaid)];
                tc_dir_dfof(:,:,iDir,1) = nanmean(data_dfof(:,:,ind_diralone),3);
                tc_dir_dfof(:,:,iDir,2) = nanmean(data_dfof(:,:,ind_dirplaid),3);
            end

            tc_stim_dfof_long = reshape(permute(tc_dir_dfof(:,:,:,1),[1 3 2]),[nStimDir.*(nFrames) nCells])';
            tc_plaid_dfof_long = reshape(permute(tc_dir_dfof(:,:,:,2),[1 3 2]),[nStimDir.*(nFrames) nCells])';
    %         [pca_coeff pca_score pca_latent pca_tsq pca_exp] = pca(tc_dir_dfof_long);
    %         stim_pca = permute(reshape(pca_score, [nFrames nStimDir nCells]),[1 3 2]);

            [X_norm, mu, sigma] = featureNormalize(tc_stim_dfof_long');
            covar = (X_norm'*X_norm);
            [U_stim,S_stim,V_stim] = svd(covar);
            stim_pca = (U_stim'*tc_stim_dfof_long);
            stim_pca = permute(reshape(stim_pca', [nFrames nStimDir nCells]),[1 3 2]);

            ev = diag(S_stim);
            %percent variance explained
            totvr = sum(sum(ev));
            pca_exp_stim = zeros(nCells,1);
            for i = 1:size(ev,1)
                pca_exp_stim(i,1) = sum(ev(1:i))/totvr;
            end


            [X_norm, mu, sigma] = featureNormalize(tc_plaid_dfof_long');
            covar = (X_norm'*X_norm);
            [U_plaid,S_plaid,V_plaid] = svd(covar);
            plaid_pca_orig = (U_plaid'*tc_plaid_dfof_long);
            plaid_pca_orig = permute(reshape(plaid_pca_orig', [nFrames nStimDir nCells]),[1 3 2]);

            ev = diag(S_plaid);
            %percent variance explained
            totvr = sum(sum(ev));
            pca_exp_plaid = zeros(nCells,1);
            for i = 1:size(ev,1)
                pca_exp_plaid(i,1) = sum(ev(1:i))/totvr;
            end

            x = get(gca, 'ColorOrder');
            x = [x; [1 0 1]; x; [1 0 1];];
            x = blues(70);
            x = x(70-15:2:70,:);
            ln_mat = ['k','r'];
            figure;
            for idir = 1:nStimDir
                plot3(stim_pca(:,1,idir), stim_pca(:,2,idir), stim_pca(:,3,idir),ln_mat(1+(idir>8)),'Marker','o','MarkerEdgeColor','k', 'MarkerFaceColor',x(idir,:))
                hold on
            end
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirResp_PCAspace_3D.pdf']),'-dpdf','-fillpage')
            figure;
            for idir = 1:nStimDir/2
                plot(stim_pca(:,1,idir), stim_pca(:,2,idir),ln_mat(1+(idir>8)),'Marker','o','MarkerEdgeColor','k', 'MarkerFaceColor',x(idir,:), 'MarkerSize', 10)
                hold on
            end
            ylim([-2.5 1.5])
            xlim([-2.5 1.5])
            axis square
            xlabel('PC1')
            ylabel('PC2')
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirResp_PCAspace_2D_PC1PC2.pdf']),'-dpdf','-fillpage')
            figure;
            for idir = 1:nStimDir/2
                plot(stim_pca(:,2,idir), stim_pca(:,3,idir),ln_mat(1+(idir>8)),'Marker','o','MarkerEdgeColor','k', 'MarkerFaceColor',x(idir,:), 'MarkerSize', 10)
                hold on
            end
            ylim([-2.5 1.5])
            xlim([-2.5 1.5])
            axis square
            xlabel('PC2')
            ylabel('PC3')
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirResp_PCAspace_2D_PC2PC3.pdf']),'-dpdf','-fillpage')
            
            pca_dir = reshape(stim_pca,[nFrames.*nCells nStimDir]);
            corr_dir = corrcoef(pca_dir);
    % 
%             figure; imagesc(corr_dir);
%             set(gca, 'XTick', 1:2:nStimDir, 'XTickLabel', num2str(stimDirs(1:2:nStimDir)'), 'YTick', 1:2:nStimDir, 'YTickLabel', num2str(stimDirs(1:2:nStimDir)')) 
%             xlabel('Grating Direction (deg)')
%             ylabel('Grating Direction (deg)')
    %         print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirResp_PCAspace_corr.pdf']),'-dpdf','-fillpage')
            corr_dir_all(:,:,start) = corr_dir;
    %         plaid_pca = zeros(nFrames,nCells,nMaskDir);
    %         for iDir = 1:nMaskDir
    %             plaid_pca(:,:,iDir) = tc_dir_dfof(:,:,iDir,2)*pca_coeff;
    %         end
            plaid_pca = (U_stim'*tc_plaid_dfof_long);
            plaid_pca = permute(reshape(plaid_pca', [nFrames nStimDir nCells]),[1 3 2]);

            stim_pca_amp = squeeze(mean(stim_pca(prewin_frames+postwin_frames:prewin_frames+postwin_frames+5,:,:),1));
            plaid_pca_amp = squeeze(mean(plaid_pca(prewin_frames+postwin_frames:prewin_frames+postwin_frames+5,:,:),1));

            stimplaid_dist = zeros(nStimDir,nStimDir);
            stimstim_dist = zeros(nStimDir,nStimDir);
            stimplaid_ang = zeros(nStimDir,nStimDir);
            stimstim_ang = zeros(nStimDir,nStimDir);
            for iplaid = 1:nStimDir
                for idir = 1:nStimDir
                    stimplaid_dist(iplaid,idir) = pdist2(plaid_pca_amp(:,iplaid)',stim_pca_amp(:,idir)');
                    stimstim_dist(iplaid,idir) = pdist2(stim_pca_amp(:,iplaid)',stim_pca_amp(:,idir)');
                    stimplaid_ang(iplaid,idir) = atan2(norm(cross(stim_pca_amp(1:3,idir),plaid_pca_amp(1:3,iplaid))),dot(stim_pca_amp(1:3,idir),plaid_pca_amp(1:3,iplaid)));
                    stimstim_ang(iplaid,idir) = atan2(norm(cross(stim_pca_amp(1:3,idir),stim_pca_amp(1:3,iplaid))),dot(stim_pca_amp(1:3,idir),stim_pca_amp(1:3,iplaid)));
                end
            end

            stimplaid_dist_all(:,:,start) = stimplaid_dist;
            stimstim_dist_all(:,:,start) = stimstim_dist;
            stimplaid_ang_all(:,:,start) = stimplaid_ang;
            stimstim_ang_all(:,:,start) = stimstim_ang;
%             figure; imagesc(stimplaid_dist)
%             colorbar
%             set(gca,'Xtick',1:2:16,'Xticklabel',num2str(stimDirs(1:2:16)'),'Ytick',1:2:16,'Yticklabel',num2str(stimDirs(1:2:16)'))
%             ylabel('Plaid Grating Direction')
%             xlabel('Single Grating Direction')
    %         

            if nCells > 99
                pca_exp_stim_all(:,start) = pca_exp_stim(1:100,:);
                pca_exp_plaid_all(:,start) = pca_exp_plaid(1:100,:);
                stim_pca_amp_all(:,:,start) = stim_pca_amp(1:100,:);
                plaid_pca_amp_all(:,:,start) = plaid_pca_amp(1:100,:);
            else
                pca_exp_stim_all(1:nCells,start) = pca_exp_stim(1:nCells,:);
                pca_exp_plaid_all(1:nCells,start) = pca_exp_plaid(1:nCells,:);
                stim_pca_amp_all(1:nCells,:,start) = stim_pca_amp(1:nCells,:);
                plaid_pca_amp_all(1:nCells,:,start) = plaid_pca_amp(1:nCells,:);
            end
%             figure;
%             for iDir = 1:nStimDir
%                 mDir = iDir+(nStimDir/4);
%                 if mDir>nStimDir
%                     mDir = mDir-nStimDir;
%                 end
%                 vDir = iDir+(nStimDir/8);
%                 if vDir>nStimDir
%                     vDir = vDir-nStimDir;
%                 end
%                 subplot(4,4,iDir)
%                 plot3(stim_pca(:,1,iDir),stim_pca(:,2,iDir), stim_pca(:,3,iDir),'-o')
%                 hold on
%                 plot3(stim_pca(:,1,mDir),stim_pca(:,2,mDir),stim_pca(:,3,mDir),'-o')
%                 plot3(stim_pca(:,1,vDir),stim_pca(:,2,vDir),stim_pca(:,3,vDir),'-o')
%                 plot3(plaid_pca(:,1,iDir),plaid_pca(:,2,iDir),plaid_pca(:,3,iDir),'-o')
%                 title(num2str(stimDirs(iDir)))
%             end
%             suptitle([mouse ' ' date])
%             print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_plaidResp_PCAspace_3D.pdf']),'-dpdf','-fillpage')
%             
%             
%             figure;
%             for iDir = 1%:nStimDir
%                 mDir = iDir+(nStimDir/4);
%                 if mDir>nStimDir
%                     mDir = mDir-nStimDir;
%                 end
%                 vDir = iDir+(nStimDir/8);
%                 if vDir>nStimDir
%                     vDir = vDir-nStimDir;
%                 end
%                 %subplot(4,4,iDir)
%                 plot(stim_pca(:,1,iDir),stim_pca(:,2,iDir),ln_mat(1),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',x(1,:), 'MarkerSize', 10)
%                 hold on
%                 plot(stim_pca(:,1,mDir),stim_pca(:,2,mDir),ln_mat(1),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',x(2,:), 'MarkerSize', 10)
%                 plot(stim_pca(:,1,vDir),stim_pca(:,2,vDir),ln_mat(1),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',x(3,:), 'MarkerSize', 10)
%                 plot(plaid_pca(:,1,iDir),plaid_pca(:,2,iDir),ln_mat(1),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',x(4,:), 'MarkerSize', 10)
%                 title(num2str(stimDirs(iDir)))
%                 ylim([-2.5 1.5])
%                 xlim([-2.5 1.5])
%                 axis square
%             end
%             suptitle([mouse ' ' date])
%             print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_plaidResp_PCAspace_2D_0deg.pdf']),'-dpdf','-fillpage')
%             print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_plaidResp_PCAspace_2D.pdf']),'-dpdf','-fillpage')
%             
%             figure;
%             hsv_map = [hsv(8); hsv(8)];
%             for iDir = 1:nStimDir
%                 plot3(stim_pca(:,1,iDir),stim_pca(:,2,iDir),stim_pca(:,3,iDir),'-o','Color', hsv_map(iDir,:))
%                 hold on
%             end
%             legend(num2str(stimDirs(1:nStimDir/2)'),'location','southwest')
%             xlabel('PCA 1')
%             ylabel('PCA 2')
%             title([mouse ' ' date '- Gratings'])
%             
%             figure;
%             hsv_map = [hsv(8); hsv(8)];
%             for iDir = 1:nStimDir
%                 plot3(plaid_pca(:,1,iDir),plaid_pca(:,2,iDir),plaid_pca(:,3,iDir),'-o','Color', hsv_map(iDir,:))
%                 hold on
%             end
%             legend(num2str(stimDirs(1:nStimDir/2)'),'location','southeast')
%             xlabel('PCA 1')
%             ylabel('PCA 2')
%             title([mouse ' ' date '- Plaids'])
%             
%             
%             
%             figure;
%             plot(squeeze(stim_pca(:,1,:)),squeeze(stim_pca(:,3,:)))
%         
%             figure;
%             for iDir = 1:nStimDir
%                 mDir = iDir+(nStimDir/4);
%                 if mDir>nStimDir
%                     mDir = mDir-nStimDir;
%                 end
%                 vDir = iDir+(nStimDir/8);
%                 if vDir>nStimDir
%                     vDir = vDir-nStimDir;
%                 end
%                 subplot(4,4,iDir)
%                 plot(stim_pca(:,2,iDir),stim_pca(:,3,iDir),'-o')
%                 hold on
%                 plot(stim_pca(:,2,mDir),stim_pca(:,3,mDir),'-o')
%                 plot(stim_pca(:,2,vDir),stim_pca(:,3,vDir),'-o')
%                 plot(plaid_pca(:,2,iDir),plaid_pca(:,3,iDir),'-o')
%                 title(num2str(stimDirs(iDir)))
%             end
%             suptitle([mouse ' ' date])
    %         print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseResp_PCAspace_2D.pdf']),'-dpdf','-fillpage')
    %         
    %         figure; 
    %         for iDir = 1:nStimDir/2
    %             subplot(1,2,1)
    %             plot(stim_pca(:,2,iDir),stim_pca(:,3,iDir),'-o')
    %             hold on
    %             subplot(1,2,2)
    %             plot(plaid_pca(:,2,iDir),plaid_pca(:,3,iDir),'-o')
    %             hold on
    %         end
    %         subplot(1,2,1)

            dir_corr = zeros(nStimDir);
            n = nCells;
            stim_pca_long = reshape(stim_pca(:,1:n,:),[nFrames.*n nMaskDir]);
            plaid_pca_long = reshape(plaid_pca(:,1:n,:),[nFrames.*n nMaskDir]);
            for iplaid = 1:nStimDir
                for istim = 1:nStimDir
                    dir_corr(iplaid,istim) = triu2vec(corrcoef(stim_pca_long(:,istim), plaid_pca_long(:,iplaid)));
                end
            end

            comp_d = [];
            vect_d = [];
            opp_d = [];
            comp_a = [];
            vect_a = [];
            opp_a = [];
            comp_c = [];
            vect_c = [];
            opp_c = [];
            for ip = 1:nStimDir
                op = ip+nStimDir/4;
                if op>nStimDir
                    op = op-nStimDir;
                end
                vp = ip+nStimDir/8;
                if vp>nStimDir
                    vp = vp-nStimDir;
                end
                cp1 = ip+nStimDir/2;
                if cp1>nStimDir
                    cp1 = cp1-nStimDir;
                end
                cp2 = ip+nStimDir/2;
                if cp2>nStimDir
                    cp2 = cp2-nStimDir;
                end
                comp_d = [comp_d stimplaid_dist(ip,ip) stimplaid_dist(ip,op)];
                vect_d = [vect_d stimplaid_dist(ip,vp)];
                opp_d = [opp_d stimplaid_dist(ip,cp1) stimplaid_dist(ip,cp2)];
                comp_a = [comp_a stimplaid_ang(ip,ip) stimplaid_ang(ip,op)];
                vect_a = [vect_a stimplaid_ang(ip,vp)];
                opp_a = [opp_a stimplaid_ang(ip,cp1) stimplaid_ang(ip,cp2)];
                comp_c = [comp_c dir_corr(ip,ip) dir_corr(ip,op)];
                vect_c = [vect_c dir_corr(ip,vp)];
                opp_c = [opp_c dir_corr(ip,cp1) dir_corr(ip,cp2)];
            end
            
            
            stimplaid_dists(start).comp = comp_d;
            stimplaid_dists(start).vect = vect_d;
            stimplaid_dists(start).opp = opp_d;
            stimplaid_angs(start).comp = comp_a;
            stimplaid_angs(start).vect = vect_a;
            stimplaid_angs(start).opp = opp_a;
            stimplaid_corrs(start).comp = comp_c;
            stimplaid_corrs(start).vect = vect_c;
            stimplaid_corrs(start).opp = opp_c;
            
    %         figure;
    %         imagesc(dir_corr)
    %         set(gca, 'XTick', 1:2:nStimDir, 'XTickLabel', num2str(stimDirs(1:2:nStimDir)'), 'YTick', 1:2:nStimDir, 'YTickLabel', num2str(stimDirs(1:2:nStimDir)')) 
    %         ylabel('Plaid Direction (deg)')
    %         xlabel('Grating Direction (deg)')
    %         title([mouse ' ' date])
    %         print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_grating&plaidCorr.pdf']),'-dpdf','-fillpage')

            dir_corr_all(:,:,start) = dir_corr;

            circ_pca_long = circshift(stim_pca_long,-nStimDir/4,2);
            comp_pca_long = stim_pca_long+circ_pca_long;
            patt_pca_long = circshift(stim_pca_long,-nStimDir/8,2);

            comp_corr = triu2vec(corrcoef(plaid_pca_long,comp_pca_long));
            patt_corr = triu2vec(corrcoef(plaid_pca_long,patt_pca_long));
            comp_patt_corr = triu2vec(corrcoef(patt_pca_long,comp_pca_long));

            Rp = ((patt_corr)-(comp_corr.*comp_patt_corr))./sqrt((1-comp_corr.^2).*(1-comp_patt_corr.^2));
            Rc = ((comp_corr)-(patt_corr.*comp_patt_corr))./sqrt((1-patt_corr.^2).*(1-comp_patt_corr.^2));
            Zp = (0.5.*log((1+Rp)./(1-Rp)))./sqrt(1./(nStimDir-3));
            Zc = (0.5.*log((1+Rc)./(1-Rc)))./sqrt(1./(nStimDir-3));
            
            Rp = ((patt_corr)-(comp_corr.*comp_patt_corr))./sqrt((1-comp_corr.^2).*(1-comp_patt_corr.^2));
            Rc = ((comp_corr)-(patt_corr.*comp_patt_corr))./sqrt((1-patt_corr.^2).*(1-comp_patt_corr.^2));
            Zp = (0.5.*log((1+Rp)./(1-Rp)))./sqrt(1./(nStimDir-3));
            Zc = (0.5.*log((1+Rc)./(1-Rc)))./sqrt(1./(nStimDir-3));

            Zc_all(:,start) = Zc;
            Zp_all(:,start) = Zp;
            save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pcaResp.mat']),'stim_pca','plaid_pca','U_stim','S_stim','V_stim','U_plaid','S_plaid','V_plaid','pca_exp_stim','pca_exp_plaid','dir_corr', 'Rc', 'Rp', 'Zc', 'Zp','stim_pca_amp','plaid_pca_amp','stimplaid_dist')
            start = start+1;
    %         b_dir = nan(1,nCells);
    %         k1_dir = nan(1,nCells);
    %         R1_dir = nan(1,nCells);
    %         R2_dir = nan(1,nCells);
    %         u1_dir = nan(1,nCells);
    %         sse_dir = nan(1,nCells);
    %         R_square_dir = nan(1,nCells);
    %         range = 1:1:360;
    %         y_dir_fit = nan(nCells,length(range));

    %         figure; 
    %         start = 1;
    %         n = 1;
    %         for iCell = 1:nCells
    %             data = [avg_resp_dir(iCell,:,1,1) avg_resp_dir(iCell,1,1,1)];
    %             theta = [deg2rad(stimDirs) 2.*pi];
    %             [b_dir(:,iCell),k1_dir(:,iCell),R1_dir(:,iCell),R2_dir(:,iCell),u1_dir(:,iCell),u2_dir(:,iCell),sse_dir(:,iCell),R_square_dir(:,iCell)] ...
    %                 = miaovonmisesfit_dir(theta,data);
    %             y_dir_fit(iCell,:) = b_dir(:,iCell)+R1_dir(:,iCell).*exp(k1_dir(:,iCell).*(cos(deg2rad(range)-u1_dir(:,iCell))-1))+R2_dir(:,iCell).*exp(k1_dir(:,iCell).*(cos(deg2rad(range)-u1_dir(:,iCell)-pi)-1));
    %             if start>25
    %                 suptitle([mouse ' ' date ' Dir'])
    %                 print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirTuningAndFitsCells' num2str(n) '.pdf']),'-dpdf', '-bestfit')
    %                 figure;
    %                 n = n+1;
    %                 start = 1;
    %             end
    %             subplot(5,5,start)
    %             errorbar(stimDirs, avg_resp_dir(iCell,:,1,1), avg_resp_dir(iCell,:,1,2), 'ok')
    %             hold on
    %             plot(range,y_dir_fit(iCell,:))
    %             start = start+1;
    %         end
    %         suptitle([mouse ' ' date ' Dir'])
    %          print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirTuningAndFitsCells' num2str(n) '.pdf']),'-dpdf', '-bestfit')
    % 
    %          figure;
    %          resp_norm_avg = zeros(nStimDir);
    %          [n edges bin]  = histcounts(rad2deg(u1_dir),[-11.25:22.5:371.25]);
    %         for idir1 = 1:nStimDir
    %             odir1 = idir1+nStimDir/2;
    %             if odir1>nStimDir
    %                 odir1 = odir1-nStimDir;
    %             end
    %             resp_norm = avg_resp_dir(:,idir1,2,1)./max(avg_resp_dir(:,:,1,1),[],2);
    %             for idir2 = 1:length(n)-1
    %                 subplot(4,4,idir1)
    %                 ind = find(bin == idir2);
    %                 resp_norm_avg(idir1,idir2) = nanmean(resp_norm(ind,:),1);
    %                 errorbar(mean(rad2deg(u1_dir(:,ind)),2), mean(resp_norm(ind,:),1),...
    %                     std(resp_norm(ind,:),[],1)./sqrt(length(ind)), std(resp_norm(ind,:),[],1)./sqrt(length(ind)),...
    %                     std(rad2deg(u1_dir(:,ind)),[],2)./sqrt(length(ind)),std(rad2deg(u1_dir(:,ind)),[],2)./sqrt(length(ind)),'ok');
    %                 hold on
    %             end
    %             resp_norm_avg(idir1,:) = circshift(resp_norm_avg(idir1,:),1-idir1,2);
    %             title(num2str(stimDirs(idir1)))
    %             xlabel('Pref dir')
    %             ylabel('Norm response')
    %             ylim([0 2])
    %             vline(stimDirs([idir1 odir1]),'k')
    %             x = stimDirs(idir1)+45;
    %             if x>360
    %                 x = x-360;
    %             end
    %             vline(x,'r')
    %             x = stimDirs([idir1 odir1])+90;
    %             if find(x>360)
    %                 x(find(x>360)) = x(find(x>360))-360;
    %             end
    %             vline(x, 'g')
    %         end
    %         suptitle([mouse ' ' date])
    %         print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_popTuning.pdf']),'-dpdf', '-bestfit')
    %         
    %         figure;
    %         errorbar(stimDirs, mean(resp_norm_avg,1), std(resp_norm_avg,[],1)./sqrt(nStimDir),'ok');
    %         ylim([0 1])
    %         idir1 = 1;
    %         odir1 = idir1+nStimDir/2;
    %         vline(stimDirs([idir1 odir1]),'k')
    %         x = stimDirs(idir1)+45;
    %         vline(x,'r')
    %         x = stimDirs([idir1 odir1])+90;
    %         vline(x, 'g')
    %         title([mouse ' ' date])
    %         xlabel('Direction (deg)')
    %         ylabel('Normalized response')
    %         print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_popTuningAlignAvg.pdf']),'-dpdf', '-bestfit')
    %         
    %         u1_dir_all = [u1_dir_all u1_dir(:,resp_ind)];
    %         avg_resp_dir_all = cat(1, avg_resp_dir_all, avg_resp_dir(resp_ind,:,:,:));
    %         resp_norm_avg_all(:,:,iexp) = resp_norm_avg;
        end
    end
    
    %%
    figure; 
    subplot(2,2,1)
    imagesc(nanmean(corr_dir_all,3))
    set(gca, 'XTick', 1:4:nStimDir, 'XTickLabel', num2str(stimDirs(1:4:nStimDir)'), 'YTick', 1:4:nStimDir, 'YTickLabel', num2str(stimDirs(1:4:nStimDir)')) 
    xlabel('Grating Direction (deg)')
    ylabel('Grating Direction (deg)')
    title('Grating-Grating')
    axis square
    caxis([0 1])
    subplot(2,2,2)
    imagesc(nanmean(dir_corr_all,3))
    set(gca, 'XTick', 1:4:nStimDir, 'XTickLabel', num2str(stimDirs(1:4:nStimDir)'), 'YTick', 1:4:nStimDir, 'YTickLabel', num2str(stimDirs(1:4:nStimDir)')) 
    xlabel('Grating Direction (deg)')
    ylabel('Plaid Direction (deg)')
    title('Grating-Plaid')
    axis square
    caxis([0 1])
    subplot(2,2,3)
    norm_corr = dir_corr_all./max(dir_corr_all,[],1);
    imagesc(nanmean(norm_corr,3));
    set(gca, 'XTick', 1:4:nStimDir, 'XTickLabel', num2str(stimDirs(1:4:nStimDir)'), 'YTick', 1:4:nStimDir, 'YTickLabel', num2str(stimDirs(1:4:nStimDir)')) 
    xlabel('Grating Direction (deg)')
    ylabel('Plaid Direction (deg)')
    title({'Grating-Plaid','Normalized'})
    axis square
    caxis([0 1])
    subplot(2,2,4)
    imagesc(nanmean(dir_corr_all,3)./max(nanmean(dir_corr_all,3),[],1))
    set(gca, 'XTick', 1:4:nStimDir, 'XTickLabel', num2str(stimDirs(1:4:nStimDir)'), 'YTick', 1:4:nStimDir, 'YTickLabel', num2str(stimDirs(1:4:nStimDir)')) 
    xlabel('Grating Direction (deg)')
    ylabel('Plaid Direction (deg)')
    title({'Grating-Plaid','Average Normalized'})
    axis square
    caxis([0 1])
    suptitle('Correlations')
    print(fullfile(summaryDir, ['grating&plaidCorr_PCA' area '.pdf']),'-dpdf','-fillpage')

    figure; 
    subplot(2,2,1)
    imagesc(nanmean(stimstim_dist_all,3))
    set(gca, 'XTick', 1:4:nStimDir, 'XTickLabel', num2str(stimDirs(1:4:nStimDir)'), 'YTick', 1:4:nStimDir, 'YTickLabel', num2str(stimDirs(1:4:nStimDir)')) 
    xlabel('Grating Direction (deg)')
    ylabel('Grating Direction (deg)')
    title('Grating-Grating')
    axis square
    caxis([0 5])
    subplot(2,2,2)
    imagesc(nanmean(stimplaid_dist_all,3))
    set(gca, 'XTick', 1:4:nStimDir, 'XTickLabel', num2str(stimDirs(1:4:nStimDir)'), 'YTick', 1:4:nStimDir, 'YTickLabel', num2str(stimDirs(1:4:nStimDir)')) 
    xlabel('Grating Direction (deg)')
    ylabel('Plaid Direction (deg)')
    title('Grating-Plaid')
    axis square
    caxis([0 5])
    subplot(2,2,3)
    norm_dist = stimplaid_dist_all./max(stimplaid_dist_all,[],1);
    imagesc(nanmean(norm_dist,3));
    set(gca, 'XTick', 1:4:nStimDir, 'XTickLabel', num2str(stimDirs(1:4:nStimDir)'), 'YTick', 1:4:nStimDir, 'YTickLabel', num2str(stimDirs(1:4:nStimDir)')) 
    xlabel('Grating Direction (deg)')
    ylabel('Plaid Direction (deg)')
    title({'Grating-Plaid','Normalized'})
    axis square
    caxis([0 1])
    subplot(2,2,4)
    imagesc(nanmean(stimplaid_dist_all,3)./max(nanmean(stimplaid_dist_all,3),[],1))
    set(gca, 'XTick', 1:4:nStimDir, 'XTickLabel', num2str(stimDirs(1:4:nStimDir)'), 'YTick', 1:4:nStimDir, 'YTickLabel', num2str(stimDirs(1:4:nStimDir)')) 
    xlabel('Grating Direction (deg)')
    ylabel('Plaid Direction (deg)')
    title({'Grating-Plaid','Average Normalized'})
    axis square
    caxis([0 1])
    suptitle('Distances')
    print(fullfile(summaryDir, ['grating&plaidDist_PCA' area '.pdf']),'-dpdf','-fillpage')

    figure; 
    subplot(2,2,1)
    imagesc(nanmean(stimstim_ang_all,3))
    set(gca, 'XTick', 1:4:nStimDir, 'XTickLabel', num2str(stimDirs(1:4:nStimDir)'), 'YTick', 1:4:nStimDir, 'YTickLabel', num2str(stimDirs(1:4:nStimDir)')) 
    xlabel('Grating Direction (deg)')
    ylabel('Grating Direction (deg)')
    title('Grating-Grating')
    axis square
    caxis([0 5])
    subplot(2,2,2)
    imagesc(nanmean(stimplaid_ang_all,3))
    set(gca, 'XTick', 1:4:nStimDir, 'XTickLabel', num2str(stimDirs(1:4:nStimDir)'), 'YTick', 1:4:nStimDir, 'YTickLabel', num2str(stimDirs(1:4:nStimDir)')) 
    xlabel('Grating Direction (deg)')
    ylabel('Plaid Direction (deg)')
    title('Grating-Plaid')
    axis square
    caxis([0 5])
    subplot(2,2,3)
    norm_ang = stimplaid_ang_all./max(stimplaid_ang_all,[],1);
    imagesc(nanmean(norm_ang,3));
    set(gca, 'XTick', 1:4:nStimDir, 'XTickLabel', num2str(stimDirs(1:4:nStimDir)'), 'YTick', 1:4:nStimDir, 'YTickLabel', num2str(stimDirs(1:4:nStimDir)')) 
    xlabel('Grating Direction (deg)')
    ylabel('Plaid Direction (deg)')
    title({'Grating-Plaid','Normalized'})
    axis square
    caxis([0 1])
    subplot(2,2,4)
    imagesc(nanmean(stimplaid_ang_all,3)./max(nanmean(stimplaid_ang_all,3),[],1))
    set(gca, 'XTick', 1:4:nStimDir, 'XTickLabel', num2str(stimDirs(1:4:nStimDir)'), 'YTick', 1:4:nStimDir, 'YTickLabel', num2str(stimDirs(1:4:nStimDir)')) 
    xlabel('Grating Direction (deg)')
    ylabel('Plaid Direction (deg)')
    title({'Grating-Plaid','Average Normalized'})
    axis square
    caxis([0 1])
    suptitle('Angles')
    print(fullfile(summaryDir, ['grating&plaidAng_PCA' area '.pdf']),'-dpdf','-fillpage')

    % [n edges bin]  = histcounts(rad2deg(u1_dir_all),[-11.25:22.5:371.25]);
    % figure;
    % for idir1 = 1:nStimDir
    %     odir1 = idir1+nStimDir/2;
    %     if odir1>nStimDir
    %         odir1 = odir1-nStimDir;
    %     end
    %     resp_norm = avg_resp_dir_all(:,idir1,2,1)./max(avg_resp_dir_all(:,:,1,1),[],2);
    %     for idir2 = 1:length(n)
    %         subplot(4,4,idir1)
    %         ind = find(bin == idir2);
    %         errorbar(mean(rad2deg(u1_dir_all(:,ind)),2), mean(resp_norm(ind,:),1),...
    %             std(resp_norm(ind,:),[],1)./sqrt(length(ind)), std(resp_norm(ind,:),[],1)./sqrt(length(ind)),...
    %             std(rad2deg(u1_dir_all(:,ind)),[],2)./sqrt(length(ind)),std(rad2deg(u1_dir_all(:,ind)),[],2)./sqrt(length(ind)),'ok');
    %         hold on
    %     end
    %     title(num2str(stimDirs(idir1)))
    %     xlabel('Pref dir')
    %     ylabel('Norm response')
    %     vline(stimDirs([idir1 odir1]),'k')
    %     x = stimDirs(idir1)+45;
    %     if x>360
    %         x = x-360;
    %     end
    %     vline(x,'r')
    %     x = stimDirs([idir1 odir1])+90;
    %     if find(x>360)
    %         x(find(x>360)) = x(find(x>360))-360;
    %     end
    %     vline(x, 'g')
    % end

    % figure;
    % errorbar(stimDirs, nanmean(nanmean(resp_norm_avg_all,1),3), nanstd(nanmean(resp_norm_avg_all,1),[],3)./sqrt(nexp),'ok');
    % ylim([0 1])
    % idir1 = 1;
    % odir1 = idir1+nStimDir/2;
    % vline(stimDirs([idir1 odir1]),'k')
    % x = stimDirs(idir1)+45;
    % vline(x,'r')
    % x = stimDirs([idir1 odir1])+90;
    % vline(x, 'g')
    % xlabel('Direction (deg)')
    % ylabel('Normalized response')
    % print(fullfile(summaryDir,['popTuningAlignAvg' area '.pdf']),'-dpdf', '-bestfit')


    figure;
    subplot(2,2,1)
    errorbar(1:100, nanmean(pca_exp_stim_all,2), nanstd(pca_exp_stim_all,[],2)./sqrt(nexp_area),'-o')
    hold on
    errorbar(1:100, nanmean(pca_exp_plaid_all,2), nanstd(pca_exp_plaid_all,[],2)./sqrt(nexp_area),'-o')
    ylabel('Explained variance')
    ylim([0 1])
    xlabel('Dimensions')
    legend({'stim','plaid'},'location','southeast')
%     [p t s] = anova2([pca_exp_stim_all(:,exp_ind) pca_exp_plaid_all(:,exp_ind)]',length(exp_ind),'off');
%     title(['p = ' num2str(chop(p,2))])
    subplot(2,2,2)
    errorbar(1:30, nanmean(pca_exp_stim_all(1:30,:,:),2), nanstd(pca_exp_stim_all(1:30,:,:),[],2)./sqrt(nexp_area),'-o')
    hold on
    errorbar(1:30, nanmean(pca_exp_plaid_all(1:30,:,:),2), nanstd(pca_exp_plaid_all(1:30,:,:),[],2)./sqrt(nexp_area),'-o')
    ylabel('Explained variance')
    ylim([0 1])
    xlabel('Dimensions')
    subplot(2,2,3)
    errorbar(1:100, nanmean(nanmean(abs(stim_pca_amp_all),2),3), nanstd(nanmean(abs(stim_pca_amp_all),2),[],3)./sqrt(nexp_area),'-o')
    hold on
    errorbar(1:100, nanmean(nanmean(abs(plaid_pca_amp_all),2),3), nanstd(nanmean(abs(plaid_pca_amp_all),2),[],3)./sqrt(nexp_area),'-o')
    ylabel('Amplitude')
    xlabel('Dimensions')
%     [p t s] = anova2([stim_pca_amp_all(:,exp_ind) plaid_pca_amp_all(:,exp_ind)]',length(exp_ind),'off');
%     title(['p = ' num2str(chop(p,2))])
    subplot(2,2,4)
    errorbar(1:30, nanmean(nanmean(abs(stim_pca_amp_all(1:30,:,:)),2),3), nanstd(nanmean(abs(stim_pca_amp_all(1:30,:,:)),2),[],3)./sqrt(nexp_area),'-o')
    hold on
    errorbar(1:30, nanmean(nanmean(abs(plaid_pca_amp_all(1:30,:,:)),2),3), nanstd(nanmean(abs(plaid_pca_amp_all(1:30,:,:)),2),[],3)./sqrt(nexp_area),'-o')
    ylabel('Amplitude')
    xlabel('Dimensions')
    print(fullfile(summaryDir,['explainedVariance&AmpPCA' area '.pdf']),'-dpdf', '-bestfit')

    stimplaid_corrs_comp = [{stimplaid_corrs.comp}];
    stimplaid_corrs_vect = [{stimplaid_corrs.vect}];
    stimplaid_corrs_opp = [{stimplaid_corrs.opp}];

    stimplaid_corrs_comp_avg = cellfun(@mean,stimplaid_corrs_comp);
    stimplaid_corrs_vect_avg = cellfun(@mean,stimplaid_corrs_vect);
    stimplaid_corrs_opp_avg = cellfun(@mean,stimplaid_corrs_opp);
    stimplaid_corrs_all = [stimplaid_corrs_comp_avg; stimplaid_corrs_vect_avg; stimplaid_corrs_opp_avg];

    stimplaid_dists_comp = [{stimplaid_dists.comp}];
    stimplaid_dists_vect = [{stimplaid_dists.vect}];
    stimplaid_dists_opp = [{stimplaid_dists.opp}];

    stimplaid_dists_comp_avg = cellfun(@mean,stimplaid_dists_comp);
    stimplaid_dists_vect_avg = cellfun(@mean,stimplaid_dists_vect);
    stimplaid_dists_opp_avg = cellfun(@mean,stimplaid_dists_opp);
    stimplaid_dists_all = [stimplaid_dists_comp_avg; stimplaid_dists_vect_avg; stimplaid_dists_opp_avg];

    stimplaid_angs_comp = [{stimplaid_angs.comp}];
    stimplaid_angs_vect = [{stimplaid_angs.vect}];
    stimplaid_angs_opp = [{stimplaid_angs.opp}];

    stimplaid_angs_comp_avg = rad2deg(cellfun(@mean,stimplaid_angs_comp));
    stimplaid_angs_vect_avg = rad2deg(cellfun(@mean,stimplaid_angs_vect));
    stimplaid_angs_opp_avg = rad2deg(cellfun(@mean,stimplaid_angs_opp));
    stimplaid_angs_all = [stimplaid_angs_comp_avg; stimplaid_angs_vect_avg; stimplaid_angs_opp_avg];

    figure;
    subplot(2,2,1) 
    errorbar(1,nanmean(stimplaid_corrs_comp_avg),nanstd(stimplaid_corrs_comp_avg)./sqrt(sum(~isnan(stimplaid_corrs_comp_avg))),'or')
    hold on
    errorbar(2,nanmean(stimplaid_corrs_vect_avg),nanstd(stimplaid_corrs_vect_avg)./sqrt(sum(~isnan(stimplaid_corrs_vect_avg))),'or')
    errorbar(3,nanmean(stimplaid_corrs_opp_avg),nanstd(stimplaid_corrs_opp_avg)./sqrt(sum(~isnan(stimplaid_corrs_opp_avg))),'or')
    hold on; plot(stimplaid_corrs_all,'k')
    set(gca,'Xtick', 1:3, 'XTickLabels',{'Component','Vector','Opposite'})
    xlim([0 4])
    ylim([0 1])
    [p t s] = anova1(stimplaid_corrs_all',[],'off');
    title(['p = ' num2str(chop(p,2))])
    ylabel('Plaid-Stim Correlation')

    subplot(2,2,2) 
    errorbar(1,nanmean(stimplaid_dists_comp_avg),nanstd(stimplaid_dists_comp_avg)./sqrt(sum(~isnan(stimplaid_dists_comp_avg))),'or')
    hold on
    errorbar(2,nanmean(stimplaid_dists_vect_avg),nanstd(stimplaid_dists_vect_avg)./sqrt(sum(~isnan(stimplaid_dists_vect_avg))),'or')
    errorbar(3,nanmean(stimplaid_dists_opp_avg),nanstd(stimplaid_dists_opp_avg)./sqrt(sum(~isnan(stimplaid_dists_opp_avg))),'or')
    hold on; plot(stimplaid_dists_all,'k')
    set(gca,'Xtick', 1:3, 'XTickLabels',{'Component','Vector','Opposite'})
    xlim([0 4])
    ylim([0 8])
    [p t s] = anova1(stimplaid_dists_all',[],'off');
    title(['p = ' num2str(chop(p,2))])
    ylabel('Plaid-Stim Neural distance')

    subplot(2,2,3) 
    errorbar(1,nanmean(stimplaid_angs_comp_avg),nanstd(stimplaid_angs_comp_avg)./sqrt(sum(~isnan(stimplaid_angs_comp_avg))),'or')
    hold on
    errorbar(2,nanmean(stimplaid_angs_vect_avg),nanstd(stimplaid_angs_vect_avg)./sqrt(sum(~isnan(stimplaid_angs_vect_avg))),'or')
    errorbar(3,nanmean(stimplaid_angs_opp_avg),nanstd(stimplaid_angs_opp_avg)./sqrt(sum(~isnan(stimplaid_angs_opp_avg))),'or')
    hold on; plot(stimplaid_angs_all,'k')
    set(gca,'Xtick', 1:3, 'XTickLabels',{'Component','Vector','Opposite'})
    xlim([0 4])
    ylim([0 90])
    [p_ang t s] = anova1(stimplaid_angs_all',[],'off');
    title(['p = ' num2str(chop(p_ang,2))])
    ylabel('Plaid-Stim Neural angle')

    subplot(2,2,4)
    plot([ones(nexp_area,1) 2.*ones(nexp_area,1)]', [Zc_all' Zp_all']', '-ok')
    hold on
    errorbar([1 2], [nanmean(Zc_all,2)' nanmean(Zp_all,2)'], [nanstd(Zc_all,[],2)'./sqrt(nexp_area) nanstd(Zp_all,[],2)'./sqrt(nexp_area)],'-or')
    set(gca,'Xtick',1:2,'XTickLabel',{'Component','Pattern'})
    xlim([0 3])
    ylim([0 5])
    ylabel('Correlation with plaid')
    [h, p_zcorr] = ttest(Zc_all,Zp_all);
    title(['p = ' num2str(chop(p_zcorr,2))])
    print(fullfile(summaryDir,['stimplaidDistance' area '.pdf']),'-dpdf', '-bestfit')
    suptitle(area)
    
    size(Zp_all)
    size(stimplaid_angs_all)
    
    save(fullfile(summaryDir,['randDir_PCA_Summary_' area '.mat']),'Zc_all','Zp_all','stimplaid_angs_all','stimplaid_corrs_all','stimplaid_dists_all')
end
