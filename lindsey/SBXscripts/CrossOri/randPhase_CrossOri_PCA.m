clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandPhase_ExptList';
eval(ds)
rc = behavConstsAV;
frame_rate = 30;
nexp = size(expt,2);
plaid_corr_all = zeros(8,8,2,nexp);
mask_corr_all = zeros(8,2,nexp);
test_corr_all = zeros(8,2,nexp);
plaid_dist_all = zeros(8,8,2,nexp);
mask_dist_all = zeros(8,2,nexp);
test_dist_all = zeros(8,2,nexp);
plaid_amp_all = zeros(8,2,nexp);
mask_amp_all = zeros(1,2,nexp);
test_amp_all = zeros(1,2,nexp);
plaid_ang_all = zeros(8,8,2,nexp);
mask_ang_all = zeros(8,2,nexp);
test_ang_all = zeros(8,2,nexp);
exp_ind = [];

%%
for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc{1};
    if strcmp(expt(iexp).driver,'SLC')
        exp_ind = [exp_ind iexp];
    ImgFolder = expt(iexp).dirFolder;
    time = expt(iexp).dirTime;
    nrun = length(ImgFolder);
    dir_run_str = catRunName(cell2mat(ImgFolder), nrun);

    LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

    fprintf([mouse ' ' date '\n'])

    %% load data

    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dir_run_str], [date '_' mouse '_' dir_run_str '_TCs.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dir_run_str], [date '_' mouse '_' dir_run_str '_input.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dir_run_str], [date '_' mouse '_' dir_run_str '_stimData.mat']))
    
    nOn = input.nScansOn;
    nOff = input.nScansOff;
    nCells = size(npSub_tc,2);
    nTrials = length(dir_mat)-1;
    dir_mat = dir_mat(:,1:end-1);
    tc_crop = npSub_tc(nOff-frame_rate+1:end-nOn-frame_rate,:);
    tc_trial = reshape(tc_crop,[nOn+nOff nTrials nCells]);
    tc_f = mean(tc_trial(1:frame_rate,:,:),1);
    tc_df = tc_trial-tc_f;
    tc_dfof = tc_df./tc_f;
    
    tc_dir_dfof = zeros(nOn+nOff,nCells,nDir);
    for idir = 1:nDir
        ind = find(dir_mat == dirs(idir));
        tc_dir_dfof(:,:,idir) = squeeze(mean(tc_dfof(:,ind,:),2));
    end
    
    tc_dir_dfof_long = reshape(permute(tc_dir_dfof,[1 3 2]),[nDir.*(nOn+nOff) nCells])';
    [X_norm, mu, sigma] = featureNormalize(tc_dir_dfof_long');
    covar = (X_norm'*X_norm);
    [U,S,V] = svd(covar);
    Z = (U'*tc_dir_dfof_long);
    
    ev = diag(S);
    %percent variance explained
    totvr = sum(sum(ev));
    for i = 1:size(ev,1)
        pca_exp(i) = sum(ev(1:i))/totvr;
    end
    
%     [pca_coeff pca_score pca_latent pca_tsq pca_exp] = pca(tc_dir_dfof_long);
%     pca_dir_dfof = reshape(pca_score, [nOn+nOff nDir nCells]);
    pca_dir_dfof = reshape(Z', [nOn+nOff nDir nCells]);
    pca_dir_amp = squeeze(mean(pca_dir_dfof(nOn:frame_rate+nOn,:,:),1));
%     figure;
%     x = get(gca, 'ColorOrder');
%     x = [x; [1 0 1]; x; [1 0 1];];
%     ln_mat = ['k','r'];
%     for idir = 1:nDir
%         plot3(pca_dir_dfof(:,idir,1), pca_dir_dfof(:,idir,2), pca_dir_dfof(:,idir,3),ln_mat(1+(idir>8)),'Marker','o','MarkerEdgeColor','k', 'MarkerFaceColor',x(idir,:))
%         hold on
%     end
%     print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dir_run_str], [date '_' mouse '_' dir_run_str '_dirResp_PCAspace_3D.pdf']),'-dpdf','-fillpage')
%     
%     figure;
%     start = 1;
%     for id1 = 1:6
%         for id2 = 1:6
%             if id1<id2
%                 subplot(4,4,start)
%                 for idir = 1:nDir
%                     plot(pca_dir_dfof(:,idir,id1), pca_dir_dfof(:,idir,id2),ln_mat(1+(idir>8)),'Color', x(idir,:))
%                     hold on
%                 end
%                 title(num2str([id1 id2]))
%                 start = start+1;
%             end
%         end
%     end
%     suptitle([mouse ' ' date])
%     print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dir_run_str], [date '_' mouse '_' dir_run_str '_dirResp_PCAspace_2D_allCombo.pdf']),'-dpdf','-fillpage')

    
    pca_dir = reshape(permute(pca_dir_dfof,[1 3 2]),[(nOn+nOff).*nCells nDir]);
    corr_dir = corrcoef(pca_dir);
%     figure; imagesc(corr_dir);
%     print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dir_run_str], [date '_' mouse '_' dir_run_str '_dirResp_PCAspace_corr.pdf']),'-dpdf','-fillpage')

    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
    nrun = length(ImgFolder);
    co_run_str = catRunName(cell2mat(ImgFolder), nrun);

    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_TCs.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_stimData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_input.mat']))
    
    prewin_frames = frame_rate;
    nFramesOn = unique(celleqel2mat_padded(input.nStimOneFramesOn));
    nFramesOff = unique(celleqel2mat_padded(input.nStimOneFramesOn))-frame_rate;
    cStimOn = celleqel2mat_padded(input.cStimOneOn);
    postwin_frames = nFramesOn+nFramesOff;
    nTrials = size(cStimOn,2);
    data_resp = nan(prewin_frames+postwin_frames,nCells,nTrials);
    for itrial = 1:nTrials-1
        if cStimOn(itrial) + postwin_frames < sz(3)
            data_resp(:,:,itrial) = npSub_tc(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
        end
    end

    data_f = mean(data_resp(1:prewin_frames,:,:),1);
    data_dfof_tc = (data_resp-data_f)./data_f;
    data_dfof_con_ph_tc_avg = nan(prewin_frames+postwin_frames, nCells, nMaskCon, nStimCon, nMaskPhas);

    for im = 1:nMaskCon
        ind_mask = find(maskCon_all == maskCons(im));
        for it = 1:nStimCon
            ind_stim = find(stimCon_all == stimCons(it));
            if it>1 & im>1
                for ip = 1:nMaskPhas
                    ind_phas = find(maskPhas_all == maskPhas(ip));
                    ind = intersect(ind_phas,intersect(ind_stim,ind_mask));
                    data_dfof_con_ph_tc_avg(:,:,im,it,ip) = squeeze(nanmean(data_dfof_tc(:,:,ind),3));
                end
            elseif it>1 || im>1
                ind = intersect(ind_mask,ind_stim);
                data_dfof_con_ph_tc_avg(:,:,im,it,1) = squeeze(nanmean(data_dfof_tc(:,:,ind),3));
            end
        end
    end
    
    test_alone = data_dfof_con_ph_tc_avg(:,:,1,2,1)';
    mask_alone = data_dfof_con_ph_tc_avg(:,:,2,1,1)';
    all_plaid = permute(squeeze(data_dfof_con_ph_tc_avg(:,:,2,2,:)), [2 1 3]);
    
    nF = size(test_alone,1);

%     test_pca = test_alone*pca_coeff;
%     mask_pca = mask_alone*pca_coeff;
%     plaid_pca = zeros(nF,nCells,nMaskPhas);
%     for ip = 1:nMaskPhas
%         plaid_pca(:,:,ip) = all_plaid(:,:,ip)*pca_coeff;
%     end
    test_pca = (U'*test_alone)';
    mask_pca = (U'*mask_alone)';
    
    plaid_pca = zeros([size(test_pca) nMaskPhas]);
    for ip = 1:nMaskPhas
        plaid_pca(:,:,ip) = permute(U'*all_plaid(:,:,ip), [2 1 3]);
    end
   
    
%     figure;
%     for ip = 1:nMaskPhas
%         subplot(3,3,ip)
%         plot3(test_pca(:,1),test_pca(:,2), test_pca(:,3),'-o')
%         hold on
%         plot3(mask_pca(:,1),mask_pca(:,2),mask_pca(:,3),'-o')
%         plot3(plaid_pca(:,1,ip),plaid_pca(:,2,ip),plaid_pca(:,3,ip),'-o')
%     end
%     subplot(3,3,9)
%     plot3(test_pca(:,1),test_pca(:,2), test_pca(:,3),'-o')
%     hold on
%     plot3(mask_pca(:,1),mask_pca(:,2),mask_pca(:,3),'-o')
%     for ip = 1:2:nMaskPhas
%         plot3(plaid_pca(:,1,ip),plaid_pca(:,2,ip),plaid_pca(:,3,ip),'-o')
%     end
%     suptitle([mouse ' ' date])
%     print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_phaseResp_PCAspace_3D.pdf']),'-dpdf','-fillpage')
% 
%     figure; 
%     for id = 1:9
%         subplot(3,3,id)
%         plot(test_pca(:,id))
%         hold on
%         plot(mask_pca(:,id))
%         if id == 1
%             leg_str = {'stim','mask'};
%         end
%         start = 1;
%         for ip = 1:2:nMaskPhas
%             plot(plaid_pca(:,id,ip))
%             leg_str{2+start} = [num2str(maskPhas(ip)) 'deg plaid'];
%             start=start+1;
%         end
%         xlabel('Frames')
%         ylabel('Amplitude')
%         if id == 1
%         legend(leg_str, 'location', 'northwest')
%         end
%         title(['Dim ' num2str(id) '- ' num2str(chop(pca_exp(id),2)) '%'])
%     end
%     suptitle([mouse ' ' date])
%     print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_pcaTC_eachDim.pdf']),'-dpdf','-fillpage')
% 
%     figure;
%     for ip = 1:nMaskPhas
%         subplot(3,3,ip)
%         plot(test_pca(:,2),test_pca(:,3),'-o')
%         hold on
%         plot(mask_pca(:,2),mask_pca(:,3),'-o')
%         plot(plaid_pca(:,1,ip),plaid_pca(:,2,ip),'-o')
%         xlabel('1st dim')
%     ylabel('2nd dim')
%     end
%     subplot(3,3,9)
%     plot(test_pca(:,2),test_pca(:,3),'-o')
%     hold on
%     plot(mask_pca(:,2),mask_pca(:,3),'-o')
%     for ip = 1:2:nMaskPhas
%         plot(plaid_pca(:,2,ip),plaid_pca(:,3,ip),'-o')
%     end
%     xlabel('2nd dim')
%     ylabel('3rd dim')
%     suptitle([mouse ' ' date])
%     print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_phaseResp_PCAspace_2D.pdf']),'-dpdf','-fillpage')
%     
     x = blues(70);
     x = x(70-15:2:70,:);
    figure;
    for ip = 1:nMaskPhas
        plot(plaid_pca(:,1,ip),plaid_pca(:,2,ip),'k','Marker','o','MarkerEdgeColor','k','MarkerFaceColor',x(ip,:), 'MarkerSize', 10)
        hold on
    end
    xlim([-4 0.5])
    ylim([-4 0.5])
    axis square
    suptitle([mouse ' ' date])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_phaseRespOnly_PCAspace_2D.pdf']),'-dpdf','-fillpage')
%     

    test_pca_avg = mean(test_pca(101:150,:),1);
    mask_pca_avg = mean(mask_pca(101:150,:),1);
    plaid_pca_avg = mean(plaid_pca(101:150,:,:),1);
    
    dims = [3 nCells];
    plaid_pca_amp = zeros(nMaskPhas,length(dims));
    test_pca_amp = zeros(1,length(dims));
    mask_pca_amp = zeros(1,length(dims));
    mask_corr = zeros(nMaskPhas,length(dims));
    test_corr = zeros(nMaskPhas,length(dims));
    plaid_corr = zeros(nMaskPhas,nMaskPhas,length(dims));
    test_pca_dist = zeros(nMaskPhas,length(dims));
    mask_pca_dist = zeros(nMaskPhas,length(dims));
    plaid_pca_dist = zeros(nMaskPhas,nMaskPhas,length(dims));
    test_pca_ang = zeros(nMaskPhas,length(dims));
    mask_pca_ang = zeros(nMaskPhas,length(dims));
    plaid_pca_ang = zeros(nMaskPhas,nMaskPhas,length(dims));

    for i = 1:length(dims)
        
        n = dims(i);
        test_pca_amp(1,i) = pdist2(test_pca_avg(:,1:n), zeros(size(test_pca_avg(:,1:n))));
        mask_pca_amp(1,i) = pdist2(mask_pca_avg(:,1:n), zeros(size(mask_pca_avg(:,1:n))));
        for ip = 1:nMaskPhas
            plaid_pca_amp(ip,i) = pdist2(plaid_pca_avg(:,1:n,ip), zeros(size(plaid_pca_avg(:,1:n,1))));
            test_pca_dist(ip,i) = pdist2(test_pca_avg(:,1:n), plaid_pca_avg(:,1:n,ip));
            mask_pca_dist(ip,i) = pdist2(mask_pca_avg(:,1:n), plaid_pca_avg(:,1:n,ip));
            if i==1
                test_pca_ang(ip,i) = atan2(norm(cross(test_pca_avg(:,1:n),plaid_pca_avg(:,1:n,ip))),dot(test_pca_avg(:,1:n),plaid_pca_avg(:,1:n,ip)));
                mask_pca_ang(ip,i) = atan2(norm(cross(mask_pca_avg(:,1:n),plaid_pca_avg(:,1:n,ip))),dot(mask_pca_avg(:,1:n),plaid_pca_avg(:,1:n,ip)));
            end
            for ip2 = 1:nMaskPhas
                plaid_pca_dist(ip,ip2,i) = pdist2(plaid_pca_avg(:,1:n,ip2), plaid_pca_avg(:,1:n,ip));
                if i == 1
                    plaid_pca_ang(ip,ip2,i) = atan2(norm(cross(plaid_pca_avg(:,1:n,ip2),plaid_pca_avg(:,1:n,ip))),dot(plaid_pca_avg(:,1:n,ip2),plaid_pca_avg(:,1:n,ip)));
                end
            end
        end

%         figure;
%         subplot(2,2,1)
%         plot(maskPhas, plaid_pca_amp(:,i)','-ok');
%         hold on
%         hline(test_pca_amp(:,i),'b')
%         hline(mask_pca_amp(:,i),'r')
%         hline(mask_pca_amp(:,i)+test_pca_amp(:,i),'m')
%         ylim([0 4])
%         xlabel('Phase (deg)')
%         ylabel('Amplitude')
%     
    
        sz = size(test_pca);
        test_pca_long = reshape(test_pca(:,1:n),[sz(1).*n 1]);
        mask_pca_long = reshape(mask_pca(:,1:n),[sz(1).*n 1]);
        plaid_pca_long = reshape(plaid_pca(:,1:n,:),[sz(1).*n nMaskPhas]);
        for ip = 1:nMaskPhas
            mask_corr(ip,i) = triu2vec(corrcoef(mask_pca_long, plaid_pca_long(:,ip)));
            test_corr(ip,i) = triu2vec(corrcoef(test_pca_long, plaid_pca_long(:,ip)));
        end
        plaid_corr(:,:,i) = corrcoef(plaid_pca_long);
        

%         subplot(2,2,3)
%         imagesc(plaid_corr(:,:,i))
%         caxis([0 1])
%         colorbar
%         set(gca,'Xtick', 1:2:nMaskPhas, 'XTickLabel', num2str(maskPhas(1:2:nMaskPhas)'), 'Ytick', 1:2:nMaskPhas, 'YTickLabel', num2str(maskPhas(1:2:nMaskPhas)'))
%         xlabel('Mask Phase (deg)')
%         ylabel('Mask Phase (deg)')
%     
%         subplot(2,2,4)
%         plot(maskPhas,test_corr(:,i))
%         hold on
%         plot(maskPhas,mask_corr(:,i))
%         plaid_corr_skip = plaid_corr(:,:,i);
%         plaid_corr_skip(find(plaid_corr(:,:,i)==1)) = nan;
%         plot(repmat(maskPhas, [nMaskPhas 1])',plaid_corr_skip,'k')
%         ylabel('Correlation')
%         xlabel('Phase (deg)')
%         %legend({'test','mask'},'location','southeast')
%         ylim([0 1])
%         
%         suptitle([mouse ' ' date ' - ' num2str(n) ' dimensions'])
%         print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_test&maskCorr_' num2str(n) 'dim.pdf']),'-dpdf','-fillpage')    
    end

    plaid_corr_all(:,:,:,iexp) = plaid_corr;
    mask_corr_all(:,:,iexp) = mask_corr;
    test_corr_all(:,:,iexp) = test_corr;
    plaid_amp_all(:,:,iexp) = plaid_pca_amp;
    mask_amp_all(:,:,iexp) = mask_pca_amp;
    test_amp_all(:,:,iexp) = test_pca_amp;
    plaid_dist_all(:,:,:,iexp) = plaid_pca_dist;
    mask_dist_all(:,:,iexp) = mask_pca_dist;
    test_dist_all(:,:,iexp) = test_pca_dist;
    plaid_ang_all(:,:,:,iexp) = plaid_pca_ang;
    mask_ang_all(:,:,iexp) = mask_pca_ang;
    test_ang_all(:,:,iexp) = test_pca_ang;

    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_pcaResp.mat']),'test_pca','mask_pca','plaid_pca','test_pca_amp','mask_pca_amp','plaid_pca_amp','U','S','V','pca_exp','mask_corr','test_corr','plaid_corr')
    end
end
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
% 
% figure;
% errorbar(maskPhas,mean(plaid_corr_all(:,1,1:5),3),std(plaid_corr_all(:,1,1:5),[],3)./sqrt(5),'-o')
% hold on
% errorbar(maskPhas,mean(plaid_corr_all(:,2,1:5),3),std(plaid_corr_all(:,2,1:5),[],3)./sqrt(5),'-o')
% ylabel('Correlation')
% xlabel('Phase (deg)')
% legend({'test','mask'},'location','southeast')
% ylim([0 1])
% 
% print(fullfile(summaryDir, 'test&maskCorr_PCA.pdf'),'-dpdf','-fillpage')
%    
plaid_corr_all_skip = plaid_corr_all;
plaid_corr_all_skip(find(plaid_corr_all == 1))= nan;
plaid_corr_all_diff = permute(permute(plaid_corr_all_skip, [1 3 4 2])-mask_corr_all,[1 4 2 3]);
dims = {'3', 'All'};
figure;
movegui('center')
for i = 1:length(dims)
    subplot(2,2,i)
    p = anova1(squeeze(plaid_amp_all(:,i,:))',[],'off');
    plot(repmat(maskPhas', [1 nexp]), squeeze(plaid_amp_all(:,i,:)./mask_amp_all(1,i,:)),'b')
    hold on
    errorbar(maskPhas, mean(squeeze(plaid_amp_all(:,i,:)./mask_amp_all(1,i,:)),2), std(squeeze(plaid_amp_all(:,i,:)./mask_amp_all(1,i,:)),[],2)./sqrt(nexp),'-ok')
    title({[dims{i} ' dimensions-'], ['p =' num2str(chop(p,2))]})
    ylim([0 3])
    ylabel('Plaid amp- norm to mask')
    xlabel('Mask phase (deg)')
    subplot(2,2,i+2)
    for ip = 1:nMaskPhas
        errorbar(maskPhas, mean(plaid_corr_all_diff(:,ip,i,:),4), std(plaid_corr_all_diff(:,ip,i,:),[],4)./sqrt(n))
        hold on
    end
    xlabel('Mask phase (deg)')
    ylabel({'Plaid/Plaid correlation -', 'Mask/Plaid Corr'})
    ylim([-0.1 0.2])
end

print(fullfile(summaryDir, 'randPhase_PCA_summary.pdf'),'-dpdf','-fillpage')

mask_dist_all_avg = squeeze(mean(mask_dist_all,1));
test_dist_all_avg = squeeze(mean(test_dist_all,1));
plaid_dist_all_avg = nan(size(test_dist_all_avg));
mask_corr_all_avg = squeeze(mean(mask_corr_all,1));
test_corr_all_avg = squeeze(mean(test_corr_all,1));
plaid_corr_all_avg = nan(size(test_corr_all_avg));
mask_ang_all_avg = rad2deg(squeeze(mean(mask_ang_all,1)));
test_ang_all_avg = rad2deg(squeeze(mean(test_ang_all,1)));
plaid_ang_all_avg = nan(size(test_ang_all_avg));
for ii = 1:2
    plaid_dist_all_tri = [];
    plaid_corr_all_tri = [];
    plaid_ang_all_tri = [];
    for i = 1:nexp
        plaid_dist_all_tri = [plaid_dist_all_tri triu2vec(plaid_dist_all(:,:,ii,i))];
        plaid_corr_all_tri = [plaid_corr_all_tri triu2vec(plaid_corr_all(:,:,ii,i))];
        plaid_ang_all_tri = [plaid_ang_all_tri triu2vec(plaid_ang_all(:,:,ii,i))];
    end
    plaid_dist_all_avg(ii,:) = mean(plaid_dist_all_tri,1);
    plaid_corr_all_avg(ii,:) = mean(plaid_corr_all_tri,1);
    plaid_ang_all_avg(ii,:) = rad2deg(mean(plaid_ang_all_tri,1));
end

stim_list = ['T','M','P'];
figure;
movegui('center')
for i = 1:2
subplot(3,2,i)
plot([test_dist_all_avg(i,:); mask_dist_all_avg(i,:); plaid_dist_all_avg(i,:)],'k')
hold on
errorbar(1,mean(test_dist_all_avg(i,:),2),std(test_dist_all_avg(i,:),[],2)./sqrt(nexp),'ok')
errorbar(2,mean(mask_dist_all_avg(i,:),2),std(mask_dist_all_avg(i,:),[],2)./sqrt(nexp),'ok')
errorbar(3,mean(plaid_dist_all_avg(i,:),2),std(plaid_dist_all_avg(i,:),[],2)./sqrt(nexp),'ok')
[p t s] = anova1([test_dist_all_avg(i,:)' mask_dist_all_avg(i,:)' plaid_dist_all_avg(i,:)'],[],'off');
t = multcompare(s,'display','off');
for ii = 1:3
    if t(ii,6)<0.05
        c= 'r';
    else
        c='k';
    end
    text(3.25,4-((ii-1).*1),[stim_list(t(ii,1)) 'v' stim_list(t(ii,2))],'Color',c)
end
set(gca,'Xtick',1:3,'XTickLabel',{'test','mask','plaid'})
xlim([0 4])
ylim([0 5])
title({['Distance (' dims{i} ' dimensions)'], ['p =' num2str(chop(p,2))]})
subplot(3,2,i+2)
plot([test_corr_all_avg(i,:); mask_corr_all_avg(i,:); plaid_corr_all_avg(i,:)],'k')
hold on
errorbar(1,mean(test_corr_all_avg(i,:),2),std(test_corr_all_avg(i,:),[],2)./sqrt(nexp),'ok')
errorbar(2,mean(mask_corr_all_avg(i,:),2),std(mask_corr_all_avg(i,:),[],2)./sqrt(nexp),'ok')
errorbar(3,mean(plaid_corr_all_avg(i,:),2),std(plaid_corr_all_avg(i,:),[],2)./sqrt(nexp),'ok')
[p t s] = anova1([test_corr_all_avg(i,:)' mask_corr_all_avg(i,:)' plaid_corr_all_avg(i,:)'],[],'off');
t = multcompare(s,'display','off');
for ii = 1:3
    if t(ii,6)<0.05
        c= 'r';
    else
        c='k';
    end
    text(3.25,.9-((ii-1).*.2),[stim_list(t(ii,1)) 'v' stim_list(t(ii,2))],'Color',c)
end
set(gca,'Xtick',1:3,'XTickLabel',{'test','mask','plaid'})
xlim([0 4])
ylim([0 1])
title({['Correlation'], ['p =' num2str(chop(p,2))]})
if i == 1
    subplot(3,2,i+4)
    plot([test_ang_all_avg(i,:); mask_ang_all_avg(i,:); plaid_ang_all_avg(i,:)],'k')
    hold on
    errorbar(1,mean(test_ang_all_avg(i,:),2),std(test_ang_all_avg(i,:),[],2)./sqrt(nexp),'ok')
    errorbar(2,mean(mask_ang_all_avg(i,:),2),std(mask_ang_all_avg(i,:),[],2)./sqrt(nexp),'ok')
    errorbar(3,mean(plaid_ang_all_avg(i,:),2),std(plaid_ang_all_avg(i,:),[],2)./sqrt(nexp),'ok')
    [p t s] = anova1([test_ang_all_avg(i,:)' mask_ang_all_avg(i,:)' plaid_ang_all_avg(i,:)'],[],'off');
    t = multcompare(s,'display','off');
    for ii = 1:3
        if t(ii,6)<0.05
            c= 'r';
        else
            c='k';
        end
        text(3.25,60-((ii-1).*15),[stim_list(t(ii,1)) 'v' stim_list(t(ii,2))],'Color',c)
    end
    set(gca,'Xtick',1:3,'XTickLabel',{'test','mask','plaid'})
    xlim([0 4])
    ylim([0 90])
    title({['Angle'], ['p =' num2str(chop(p,2))]})
end
end
print(fullfile(summaryDir, 'randPhase_PCA_DistCorrAng_summary.pdf'),'-dpdf','-fillpage')

figure;
plot([test_ang_all_avg(1,exp_ind); plaid_ang_all_avg(1,exp_ind)],'k')
hold on
errorbar(1,mean(test_ang_all_avg(1,exp_ind),2),std(test_ang_all_avg(1,exp_ind),[],2)./sqrt(nexp),'ok')
errorbar(2,mean(plaid_ang_all_avg(1,exp_ind),2),std(plaid_ang_all_avg(1,exp_ind),[],2)./sqrt(nexp),'ok')
set(gca,'Xtick',1:2,'XTickLabel',{'test','plaid'})
xlim([0 4])
ylim([0 90])
print(fullfile(summaryDir, 'randPhase_PCA_TestPlaidAng_summary.pdf'),'-dpdf','-fillpage')

