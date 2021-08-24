%% START HERE with frequency-dfof response accounting for frequency and direction
dfof_stim = zeros(nDir,nSF,nCells1,3);
trialInd = cell(nDir,nSF,3);
  for idir = 1:nDir
    ind_dir = find(dir_mat(1:640) == dirs(idir));
    for iSF = 1:nSF
        ind_SF = find(SF_mat(1:640) == SFs(iSF));
        trialInd{idir,iSF,3} = intersect(ind_dir,ind_SF);
        ntot = length(trialInd{idir,iSF,3});
        ntemp = ceil(ntot/2);
        ind_temp = randperm(ntot,ntemp);
        trialInd{idir,iSF,1} = trialInd{idir,iSF,3}(ind_temp);
        ind_temp2 = setdiff(1:ntot,ind_temp);
        trialInd{idir,iSF,2} = trialInd{idir,iSF,3}(ind_temp2);
        dfof_stim(idir,iSF,:,1) = mean(data_dfof_resp(trialInd{idir,iSF,1},:),1);
        dfof_stim(idir,iSF,:,2) = mean(data_dfof_resp(trialInd{idir,iSF,2},:),1);
        dfof_stim(idir,iSF,:,3) = mean(data_dfof_resp(trialInd{idir,iSF,3},:),1);
    end
  end
dfof_stim2 = zeros(nDir,nSF,length(cell_list),3);
trialInd2 = cell(nDir,nSF,3);
  for idir = 1:nDir2
    ind_dir2 = find(dir_mat2(1:640) == dirs2(idir));
    for iSF = 1:nSF
        ind_SF2 = find(SF_mat2(1:640) == SFs2(iSF));
        trialInd2{idir,iSF,3} = intersect(ind_dir2,ind_SF2);
        ntot = length(trialInd2{idir,iSF,3});
        ntemp = ceil(ntot/2);
        ind_temp = randperm(ntot,ntemp);
        trialInd2{idir,iSF,1} = trialInd2{idir,iSF,3}(ind_temp);
        ind_temp2 = setdiff(1:ntot,ind_temp);
        trialInd2{idir,iSF,2} = trialInd2{idir,iSF,3}(ind_temp2);
        dfof_stim2(idir,iSF,:,1) = mean(data_dfof_resp2(trialInd2{idir,iSF,1},:),1);
        dfof_stim2(idir,iSF,:,2) = mean(data_dfof_resp2(trialInd2{idir,iSF,2},:),1);
        dfof_stim2(idir,iSF,:,3) = mean(data_dfof_resp2(trialInd2{idir,iSF,3},:),1);
    end
  end
dfof_stim_rect = dfof_stim;
dfof_stim_rect(find(dfof_stim<0))=0;
dfof_stim2_rect = dfof_stim2;
dfof_stim2_rect(find(dfof_stim2<0))=0;
f_max = squeeze(max(max(dfof_stim_rect)));
f_max2 = squeeze(max(max(dfof_stim2_rect)));
c_pass = NaN(nCells1,100);
c_pass2 = NaN(nCells2,100);
S = NaN(nCells1,1);
S2 = NaN(nCells2,1);
t = NaN(nCells1,100);
t2 = NaN(nCells2,100);
[n n2] = subplotn(nCells1);
figure;
for iCell = 1:nCells1
    t(iCell,:) = linspace(0.0,f_max(iCell,3),100);
    for it = 1:100
        iT = t(iCell,it);
    c_pass(iCell,it) = sum(squeeze(max(dfof_stim_rect(:,:,iCell,3),[],2))>=iT)/size(dfof_stim_rect,1);
    end
    S(iCell) = 1 - (2*(sum(c_pass(iCell,:),2)/100));
    subplot(n,n2,iCell)
    bar(t(iCell,:),c_pass(iCell,:),'r')
    title(['cell ' num2str(iCell)])
    xlabel(['Threshold'])
    ylabel(['Responses'])
    print(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], ['threshold.pdf']),'-dpdf', '-bestfit')
end
figure;
ecdf(S(goodCells),'Bounds','on');hold on;ecdf(S(okayCells),'Bounds','on')
xlabel('selectivity');ylabel('fraction of cells')
legend({'goodCells' 'okayCells'})
figure;
for iCell = 1:nCells2
    t2(iCell,:) = linspace(0.0,f_max2(iCell,3),100);
    for it = 1:100
        iT = t2(iCell,it);
    c_pass2(iCell,it) = sum(squeeze(max(dfof_stim2_rect(:,:,iCell,3),[],2))>=iT)/size(dfof_stim2_rect,1);
    end
    S2(iCell) = 1 - (2*(sum(c_pass2(iCell,:),2)/100));
    subplot(n,n2,iCell)
    bar(t2(iCell,:),c_pass2(iCell,:),'r')
    title(['cell ' num2str(iCell)])
    xlabel(['Threshold'])
    ylabel(['Responses'])
    print(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], ['threshold2.pdf']),'-dpdf', '-bestfit')
end
save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_threshold.mat']),'t','c_pass','t2','c_pass2','S','S2');

%% ttest with freq
base_wind = 1+nOff1-nOn1:nOff1;
dfof_base = squeeze(mean(data_dfof1(base_wind,cell_list,:),1))';
resp_wind = nOff1+10:nOff1+nOn1;
dfof_resp = squeeze(mean(data_dfof1(resp_wind,cell_list,:),1))';
h1 = NaN(nDir,nCells2);
p1 = NaN(nDir,nCells2);
for iDir = 1:nDir
    ind = find(tGratingDir1(1:640)==dirs(iDir));
    x1 = dfof_base(ind,:);
    y1 = dfof_resp(ind,:);
    [h1(iDir,:),p1(iDir,:)] = ttest(x1,y1,'dim',1,'Alpha',0.05./(nDir-1),'tail','left');
end
dfof_base = squeeze(mean(data_dfof2(base_wind,:,:),1))';
dfof_resp = squeeze(mean(data_dfof2(resp_wind,:,:),1))';
h2 = NaN(nDir,nCells2);
p2 = NaN(nDir,nCells2);
for iDir = 1:nDir
    ind = find(tGratingDir1(1:640)==dirs(iDir));
    x2 = dfof_base(ind,:);
    y2 = dfof_resp(ind,:);
    [h2(iDir,:),p2(iDir,:)] = ttest(x2,y2,'dim',1,'Alpha',0.05./(nDir-1),'tail','left');
end
save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_ttest.mat']),'h1','h2','p1','p2');

%% signal corr with freq
max_dfof_stim = squeeze(max(dfof_stim(:,:,cell_list),[],2));
max_dfof_stim2 = squeeze(max(dfof_stim2,[],2));
dfof_stim_short = dfof_stim(:,:,cell_list);
signal_corr_half1 = NaN(nCells2,1);
signal_corr_half2 = NaN(nCells2,1);
signal_corr_day = NaN(nCells2,1);
signal_corr_day_half = NaN(nCells2,1);
signal_corr_maxSF = NaN(nCells2,1);
signal_corr_avgSF = NaN(nCells2,1);
for iCell = 1:nCells2
        signal_corr_half1(iCell,1) = triu2vec(corrcoef(dfof_stim(:,:,iCell,1),dfof_stim(:,:,iCell,2)));
        signal_corr_half2(iCell,1) = triu2vec(corrcoef(dfof_stim2(:,:,iCell,1),dfof_stim2(:,:,iCell,2)));
        signal_corr_day(iCell,1) = triu2vec(corrcoef(dfof_stim_short(:,:,iCell,3),dfof_stim2(:,:,iCell,3))); 
        signal_corr_day_half(iCell,1) = triu2vec(corrcoef(dfof_stim_short(:,:,iCell,1),dfof_stim2(:,:,iCell,1)));
        signal_corr_maxSF(iCell,1) = triu2vec(corrcoef(max_dfof_stim(:,iCell,1),max_dfof_stim2(:,iCell,2)));
        signal_corr_avgSF(iCell,1) = triu2vec(corrcoef(mean(dfof_stim(:,:,iCell,1),2),mean(dfof_stim2(:,:,iCell,1),2)));
end
% save(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_signalCorr.mat']),'signal_corr_half1','signal_corr_day');
save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_signalCorr.mat']),'signal_corr_half1','signal_corr_half2','signal_corr_day','signal_corr_day_half','signal_corr_avgSF','signal_corr_maxSF');

%% vonmises with frequency
theta = deg2rad(dirs); 
thetafine = deg2rad(1:360); 
for iCell = 1:nCells1
    try
    [b_hat(iCell),k1_hat(iCell),R1_hat(iCell),R2_hat(iCell),u1_hat(iCell),u2_hat(iCell),sse(iCell),R_square(iCell)] = miaovonmisesfit_dir(theta,squeeze(mean(dfof_stim(:,:,iCell,3),2))');
    y_fit(iCell,:) = b_hat(iCell)+R1_hat(iCell).*exp(k1_hat(iCell).*(cos(thetafine-u1_hat(iCell))-1))+R2_hat(iCell).*exp(k1_hat(iCell).*(cos(thetafine-u2_hat(iCell))-1));
    catch
        b_hat(iCell) = NaN;
        k1_hat(iCell) = NaN;
        R1_hat(iCell) = NaN;
        R2_hat(iCell) = NaN;
        u1_hat(iCell) = NaN;
        u2_hat(iCell) = NaN;
        sse_hat(iCell) = NaN;
        R_square_hat(iCell)= NaN; 
        y_fit(iCell,:) = NaN(1,length(thetafine));
    end
end
figure;
ecdf(k1_hat(goodCells),'Bounds','on');hold on;ecdf(k1_hat(okayCells),'Bounds','on')
xlabel('k value');ylabel('fraction of cells')
legend({'goodCells' 'okayCells'})
for iCell = 1:nCells2
    try
    [b_hat2(iCell),k1_hat2(iCell),R1_hat2(iCell),R2_hat2(iCell),u1_hat2(iCell),u2_hat2(iCell),sse2(iCell),R_square2(iCell)] = miaovonmisesfit_dir(theta,squeeze(mean(dfof_stim2(:,:,iCell,3),2))');
    y_fit2(iCell,:) = b_hat2(iCell)+R1_hat2(iCell).*exp(k1_hat2(iCell).*(cos(thetafine-u1_hat2(iCell))-1))+R2_hat2(iCell).*exp(k1_hat2(iCell).*(cos(thetafine-u2_hat2(iCell))-1));
    catch
        b_hat2(iCell) = NaN;
        k1_hat2(iCell) = NaN;
        R1_hat2(iCell) = NaN;
        R2_hat2(iCell) = NaN;
        u1_hat2(iCell) = NaN;
        u2_hat2(iCell) = NaN;
        sse_hat2(iCell) = NaN;
        R_square_hat2(iCell)= NaN; 
        y_fit2(iCell,:) = NaN(1,length(thetafine));
    end
end
save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_vonmises.mat']),'k1_hat','k1_hat2');
