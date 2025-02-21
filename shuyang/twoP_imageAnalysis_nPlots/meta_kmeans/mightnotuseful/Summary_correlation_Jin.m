
% sFrame1d = reshape(sFrame, x*y,size(sFrame,3));
% fFrame1d = reshape(fFrame, x*y,size(fFrame,3));
% inclusterSCorr = [];
% inclusterfCorr = [];
% SCorrindx = zeros(x*y, x*y);
% FCorrindx = zeros(x*y, x*y);
% Compindx = zeros(x*y, x*y);
% Compindx2 = zeros(x*y, x*y);

clear;
days = {'151009_img30', '151011_img30', '151212_img32', '160314_img38', '160315_img38', ...
    '160606_img46', '160722_img53', '160904_img55'};
mouseID = {'img30', 'img30', 'img32', 'img38', 'img38', 'img46', 'img53', 'img55'};
out_dir = 'Y:\home\jake\Analysis\WF Lever Analysis\Meta-kmeans_output_dir\';
out_dir2 = 'Y:\public\Motor Timing Paper\Ziye\WF_MetaKmeans\';
% for ii = 1:length(days)
%     load([out_dir, days{ii}, '\', 'cluster_wholeMovie_clean.mat']);
%     load([out_dir, days{ii}, '\', 'Correlation.mat']);
% 
%     [x,y,~] = size(sFrame);
%     sFrame1d = reshape(sFrame, x*y,size(sFrame,3));
%     fFrame1d = reshape(fFrame, x*y,size(fFrame,3));
% 
%     bwRP =  regionprops(clusters, 'Area');
%     cSize = cell2mat(struct2cell(bwRP));
%     [x,y] = size(clusters);
%     cSizeNorm = cSize/(x*y)*100;
% 
%     indx = reshape(clusters, x*y,1);
%     allclusters = reshape(clusters, x*y,1);
%     nC = max(allclusters);
%     indx(indx>0) = 2;
% %
%     k = 1;
%     col = [0.9,0,0];
%     outSCorr_avg = []; outFCorr_avg = []; outSCorrN = []; outSCorrNN = []; outFCorrN = []; outFCorrNN = [];
%     for i = 1:nC-1
%         insCluster = (allclusters == i);
%         if sum(insCluster)/(x*y) < 0.03
%             continue
%         else
%             indx(insCluster) = 9;
%             for j = i+1 : nC
%                 outCluster = (allclusters == j);
%                 if sum(outCluster)/(x*y) < 0.03
%                     continue
%                 else
%                     indx(outCluster) = 10;
%                     outSCorr{k} = corr(sFrame1d(insCluster,:)', sFrame1d(outCluster,:)');
%                     outFCorr{k} = corr(fFrame1d(insCluster,:)', fFrame1d(outCluster,:)');
%                     outSCorr_avg = [outSCorr_avg nanmean(outSCorr{k}(:))];
%                     outFCorr_avg = [outFCorr_avg nanmean(outFCorr{k}(:))];
% 
%                     outSCorr_temp = outSCorr{k}(:);
%                     outFCorr_temp = outFCorr{k}(:);
% 
%                     removeIdx = isnan(outSCorr_temp) | isnan(outFCorr_temp);
% 
%                     outSCorr_temp(removeIdx) = [];
% 
%                     outFCorr_temp(removeIdx) = [];
% 
%                     if abs(i-j) == 1
%                         outSCorrN = [outSCorrN nanmean(outSCorr{k}(:))];
%                         outFCorrN = [outFCorrN nanmean(outFCorr{k}(:))];
%                     else
%                         outSCorrNN = [outSCorrNN nanmean(outSCorr{k}(:))];
%                         outFCorrNN = [outFCorrNN nanmean(outFCorr{k}(:))];
%                     end
% 
%                     fig = figure;
%                     subplot(1,2,2); %scatter(outSCorr_temp, outFCorr_temp, 4, 'MarkerEdgeColor', [0.5,0.5,0.5]);
%                     hold on; plot([0:0.1:1], [0:0.1:1], 'k');xlabel('Success Corr'); ylabel('Fail Corr')
%                     h = errorbarxy(mean(outSCorr_temp), mean(outFCorr_temp),  std(outSCorr_temp)./sqrt(length(outSCorr_temp))...
%                         , std(outFCorr_temp)./sqrt(length(outFCorr_temp)),{'o', col, col, col});
%                     set(h.hMain,'LineWidth', 1);
% 
%                     indx = reshape(indx, x, y);
%                     subplot(1,2,1); imagesc(indx); xlabel(['Cluster #', num2str(i), ' and Cluster #', num2str(j)]);
%                     indx(outCluster) = 2;
%                     k = k +1;
% 
% %                     saveas(fig, [out_dir, days{ii}, '\', 'cluster', num2str(i), '_' num2str(j), '.tif']);
%                     saveas(fig, [out_dir2, days{ii}, '\', 'cluster', num2str(i), '_' num2str(j), '.fig']);
%                 end
%             end
%             indx(insCluster) = 2;
%         end
%     end
% 
%     close all;
% %
%     inSCorr_avg = []; inFCorr_avg = [];
%     fig = figure;
%     k =1;
%     if nC <= 6
%         subRow = 3;
%     else
%         subRow = 4;
%     end
%     for i = 1:nC
% 
%         clusterIdx = (allclusters == i);
%         if sum(clusterIdx) < 0.03
%             continue
%         else
%             indx(clusterIdx) = 10;
% 
%             df_f_Smean = mean(sFrame1d(clusterIdx,:));
% 
%             inSCorr{k} = corr(sFrame1d(clusterIdx,:)', sFrame1d(clusterIdx,:)');
%             inFCorr{k} = corr(fFrame1d(clusterIdx,:)', fFrame1d(clusterIdx,:)');
% 
%             inSCorr_temp = inSCorr{k}(:);
%             inFCorr_temp = inFCorr{k}(:);
%             removeIdx = (isnan(inSCorr_temp) | isnan(inFCorr_temp));
% 
%             inSCorr_temp(removeIdx) = [];
%             inFCorr_temp(removeIdx) = [];
% 
%             inSCorr_avg = [inSCorr_avg nanmean(inSCorr{k}(:))];
%             inFCorr_avg = [inFCorr_avg nanmean(inFCorr{k}(:))];
% 
% 
%             subplot(subRow,4,2*k); %scatter(inSCorr_temp, inFCorr_temp, 4, 'MarkerEdgeColor', [0.5,0.5,0.5]); %axis([0 1 0 1])
%             hold on; plot([0:0.1:1], [0:0.1:1], 'k'); xlabel('Success Corr'); ylabel('Fail Corr')
%             h = errorbarxy(mean(inSCorr_temp), mean(inFCorr_temp),  std(inSCorr_temp)./sqrt(length(inSCorr_temp))...
%                         , std(inFCorr_temp)./sqrt(length(inFCorr_temp)),{'o', col, col, col});
%             set(h.hMain,'LineWidth', 1);
% 
%             clusters = reshape(indx, x, y);
%             subplot(subRow,4,2*(k-1)+1); imagesc(clusters); xlabel(['Cluster #', num2str(i)]);
%             %
%             % %     SCorrTri = tril(inSCorr{k});
%             % %     SCorrVec = SCorrTri(:);
%             %     SCorrindx = zeros(x*y, x*y);
%             %     FCorrindx = zeros(x*y, x*y);
%             %
%             %     SCorrindx(clusterIdx, clusterIdx) = inSCorr{k};
%             %     FCorrindx(clusterIdx, clusterIdx) = inFCorr{k};
%             %
%             %     Compindx(SCorrindx > FCorrindx) = i;
%             %     Compindx2(SCorrindx < FCorrindx) = i;
% 
%             indx(clusterIdx) = 2;
%             %     inclusterSCorr = [];
%             %     inclusterfCorr = [];
%             %     for m = 1:length(clusterIdx)-1
%             %
%             %         for n = m+1 : length(clusterIdx)
%             %
%             %             inclusterSCorr = [inclusterSCorr corr(sFrame1d(clusterIdx(m),:)', sFrame1d(clusterIdx(n),:)')];
%             %             inclusterfCorr = [inclusterfCorr corr(fFrame1d(clusterIdx(m),:)', fFrame1d(clusterIdx(n),:)')];
%             %         end
%             %     end
%             %     SCorrCell{i} = inclusterSCorr;
%             %     FCorrCell{i} = inclusterfCorr;
%             k = k+1;
%         end
%     end
% %     saveas(fig, [out_dir, days{ii}, '\', 'withinCluster_Correlation.tif']);
%     saveas(fig, [out_dir2, days{ii}, '\', 'withinCluster_Correlation.fig']);
%
%     save([out_dir, days{ii}, '\', 'Correlation.mat'], 'outSCorr', 'outFCorr', 'outSCorr_avg', 'outFCorr_avg', 'inSCorr', 'inFCorr', ...
%         'inSCorr_avg', 'inFCorr_avg', 'sFrame', 'fFrame', 'cSizeNorm', 'outSCorrN', 'outSCorrNN', 'outFCorrN', 'outFCorrNN');
%
%     close all;
%     % fig = figure; s1 = scatter(outSCorr_avg, outFCorr_avg, 12, 'MarkerEdgeColor', [0.5,0.5,0.5]);
%     % hold on; plot(0:0.1:1, 0:0.1:1, 'k');
%     %
%     % s2 = scatter(inSCorr_avg, inFCorr_avg, 12, 'MarkerEdgeColor', [0.9,0,0]);
%     % xlabel('Average Correlation during Correct'); ylabel('Average Correlation during Early');
%     % legend([s1, s2], 'between clusters', 'within cluster'); title([subID(1:6), ' ', subID(8:end)]);
%     %
%     % saveas(fig, [outdir, 'Cluster_MeanCorrelation_Scatter.fig']);
%     % print([outdir 'Cluster_MeanCorrelation_Scatter.eps'], '-depsc');
%     % print([outdir 'Cluster_MeanCorrelation_Scatter.pdf'], '-dpdf');
%
% end
% fig = figure;
for ii = [8]
    load([out_dir, days{ii}, '\', 'cluster_wholeMovie_clean.mat']);
    load([out_dir, days{ii}, '\', 'Correlation.mat']);
    
    [x,y,failT] = size(fFrame);
    
    fFrame2d = reshape(fFrame,x*y,failT);
    
    SuccT = size(sFrame,3);
    
    sFrame2d = reshape(sFrame,x*y,SuccT);
    
    randS = randperm(SuccT/5-1);
    sFrame_rand = [];
    fFrame_Rmean = [];
    sFrame_Rmean = [];
    
    indx = reshape(clusters, x*y,1);
    allclusters = reshape(clusters, x*y,1);
    nC = max(allclusters);
    
    for i = 1:nC
        clusterIdx = (allclusters == i);
        sFrame_rand = [];
        for j = 1:failT/5
            sFrame_rand = cat(2, sFrame_rand, sFrame2d(clusterIdx,randS(j)*5+1:randS(j)*5+5));
        end
        fFrame_all = fFrame2d(clusterIdx,:);
%         sFrame_Rmean = [sFrame_Rmean nanmean(nanmean(sFrame_rand))];
%         fFrame_Rmean = [fFrame_Rmean nanmean(nanmean(fFrame_all))];
        sFrame_rand( ~any(sFrame_rand,2), : ) = []; 
        sFrame_rand( :, ~any(sFrame_rand,1) ) = [];
        [row, col] = find(sFrame_rand==0);
        sFrame_rand(row,:) = []; sFrame_rand(:,col) = 0;
        fFrame_all( ~any(fFrame_all,2), : ) = []; 
        fFrame_all( :, ~any(fFrame_all,1) ) = [];
        [row, col] = find(fFrame_all==0);
        fFrame_all(row,:) = []; fFrame_all(:,col) = 0;
        
%         [sC,sT] = size(sFrame_rand);
%         if mod(sT,2)==0
%             geoMeanS1 = nthroot(prod(sFrame_rand(:,1:sT-1),2),sT-1);
%             geoMeanF1 = nthroot(prod(fFrame_all(:,1:sT-1),2),sT-1);
%         else
%             geoMeanS1 = nthroot(prod(sFrame_rand(:,1:sT),2),sT);
%             geoMeanF1 = nthroot(prod(fFrame_all(:,1:sT),2),sT);
%             
%         end
%         geoMeanS1(geoMeanS1==0) = [];
%         geoMeanF1(geoMeanF1==0) = [];
%         sCS = length(geoMeanS1);
%         sCF = length(geoMeanF1);
%         if mod(sCS,2) == 0
%             geoMeanS2 = nthroot(prod(geoMeanS1(1:sCS - 1)*20), sCS - 1)/20;
%         else
%             geoMeanS2 = nthroot(prod(geoMeanS1*20), sCS)/20;
%         end
%         
%         if mod(sCF,2) == 0
%             geoMeanF2 = nthroot(prod(geoMeanF1(1:sCF - 1)*20), sCF - 1)/20;
%         else
%             geoMeanF2 = nthroot(prod(geoMeanF1*20), sCF)/20;
%         end
%         
%         if ii == 3 && i == 2
%             if mod(sCF,2) == 0
%                 geoMeanF2 = nthroot(prod(geoMeanF1(1:sCF - 1)*40), sCF - 1)/40;
%             else
%                 geoMeanF2 = nthroot(prod(geoMeanF1*40), sCF)/40;
%             end
%         end
%         sFrame_Rmean = [sFrame_Rmean geoMeanS2];
%         fFrame_Rmean = [fFrame_Rmean geoMeanF2];
        sFrame_Rmean = [sFrame_Rmean nanmean(nanmean(sFrame_rand))];
        fFrame_Rmean = [fFrame_Rmean nanmean(nanmean(fFrame_all))];
    end
    
    df_f_Smean{ii} = sFrame_Rmean;
    df_f_Fmean{ii} = fFrame_Rmean;
    inSCorr_avg_tot{ii} = inSCorr_avg;
    
    inFCorr_avg_tot{ii} = inFCorr_avg;
    
%     if ii ==1
%         s1 = scatter( inSCorr_avg_tot{ii}, df_f_Smean{ii}, 12, 'MarkerEdgeColor', [0.5, 0.5, 0.5]);
%         hold on;
%         s2 = scatter(inFCorr_avg_tot{ii}, df_f_Fmean{ii}, 12, 'MarkerEdgeColor', [0.9, 0, 0]);
%     else
%         scatter(inSCorr_avg_tot{ii}, df_f_Smean{ii}, 12, 'MarkerEdgeColor', [0.5, 0.5, 0.5]);
%         hold on;
%         scatter(inFCorr_avg_tot{ii}, df_f_Fmean{ii}, 12, 'MarkerEdgeColor', [0.9, 0, 0]);
%     end
        outSCorrN_avg_tot{ii} = outSCorrN;
    
        outSCorrNN_avg_tot{ii} = outSCorrNN;
    
        outFCorrN_avg_tot{ii} = outFCorrN;
    
        outFCorrNN_avg_tot{ii} = outFCorrNN;
end
totCorrS = [cell2mat(inSCorr_avg_tot)' cell2mat(df_f_Smean)'];
[~, ix] = sort(totCorrS(:,2));
sortCorrS = totCorrS(ix, :);

totCorrF = [cell2mat(inFCorr_avg_tot)' cell2mat(df_f_Fmean)'];
[~, ix] = sort(totCorrF(:,2));
sortCorrF = totCorrF(ix, :);

range = max(sortCorrS(end,2)-sortCorrS(1,2), sortCorrF(end,2)-sortCorrF(1,2));
grid = range / 5;
[~, idx] = histc(sortCorrS(:,2), sortCorrS(1,2):grid:range);
binCorrS = accumarray(idx(:),sortCorrS(:,1),[],@mean);
bindf_f_S = accumarray(idx(:),sortCorrS(:,2),[],@mean);
binCorrS(binCorrS==0) = [];
bindf_f_S(bindf_f_S==0)=[];

[~, idx] = histc(sortCorrF(:,2), sortCorrF(1,2):grid:range);
binCorrF = accumarray(idx(:),sortCorrF(:,1),[],@mean);
bindf_f_F = accumarray(idx(:),sortCorrF(:,2),[],@mean);
binCorrF(binCorrF==0)=[];
bindf_f_F(bindf_f_F==0)=[];

figure;
scatter(bindf_f_S, binCorrS, 12, 'MarkerEdgeColor', [0.5, 0.5, 0.5]);
hold on;
scatter(bindf_f_F, binCorrF, 12, 'MarkerEdgeColor', [0.9, 0, 0]);
xlabel('binned dF/F'); ylabel('binned correlation');

%% linear fit correlation and df/f
LinearCoeffS = polyfit(totCorrS(:,1), totCorrS(:,2), 1);
CorrSfit = polyval(LinearCoeffS, totCorrS(:,1));

LinearCoeffF = polyfit(totCorrF(:,1), totCorrF(:,2), 1);
CorrFfit = polyval(LinearCoeffF, totCorrF(:,1));
figure; s1 = scatter(totCorrS(:,1), totCorrS(:,2), 12, 'MarkerEdgeColor', [0.5, 0.5, 0.5]);
hold on;
plot(totCorrS(:,1), CorrSfit, 'color', [0.5, 0.5,0.5]);
s2 = scatter(totCorrF(:,1), totCorrF(:,2), 12, 'MarkerEdgeColor', [0.9, 0, 0]);
plot(totCorrF(:,1), CorrFfit, 'color', [0.9, 0, 0]);
ylabel('dF/F'); xlabel('correlation coefficient');
legend([s1, s2], 'Correct', 'Early');

%%
plot([-1:0.1:1], [-1:0.1:1], 'k');
hold on;
y2 = df_f_Fmean;
y1 = inFCorr_avg_tot;
y1 = cell2mat(y1);
y1_size = size(y1,2);
y1_mean = mean(y1,2);
y2 = cell2mat(y2);
y2_size = size(y2,2);
y2_mean = mean(y2,2);
col = [0.9 0 0];
h = errorbarxy(y1_mean, y2_mean,  std(y1)./sqrt(y1_size), std(y2)./sqrt(y2_size),{'o', col, col, col});
set(h.hMain,'LineWidth', 1);

hold on;
y2 = df_f_Smean;
y1 = inSCorr_avg_tot;
y1 = cell2mat(y1);
y1_size = size(y1,2);
y1_mean = mean(y1,2);
y2 = cell2mat(y2);
y2_size = size(y2,2);
y2_mean = mean(y2,2);
col = [0.5 0.5 0.5];
h = errorbarxy(y1_mean, y2_mean,  std(y1)./sqrt(y1_size), std(y2)./sqrt(y2_size),{'o', col, col, col});
set(h.hMain,'LineWidth', 1);
xlabel('Averaged Within Cluster Correlation'); ylabel('Averaged dF/F');
legend([s1, s2], 'Correct', 'Early');

col_mat = [ 0.9  0.9  0;
    1  0  1;
    0  1  1;
    0.5  0  0;
    0  1  0;
    0  0  1;
    1  0.6  1;
    0  0  0;
    1  0.8 0.4
    0  0.5 0.7
    0.5 0.4 0];
fig = figure;
scatter_plot(mouseID, inSCorr_avg_tot, inFCorr_avg_tot, col_mat);
hold on; plot([0:0.1:1], [0:0.1:1], 'k'); xlabel('Success Corr'); ylabel('Fail Corr')

fig = figure;
subplot(1,3,1);
scatter_plot(mouseID, inSCorr_avg_tot, inFCorr_avg_tot, col_mat)
hold on; plot([0:0.1:1], [0:0.1:1], 'k'); xlabel('Correct Corr'); ylabel('Early Corr')
axis square
title('Within Cluster');

subplot(1,3,2);
scatter_plot(mouseID, outSCorrN_avg_tot, outFCorrN_avg_tot, col_mat)
hold on; plot([0:0.1:1], [0:0.1:1], 'k'); xlabel('Correct Corr'); ylabel('Early Corr')
axis square
title('Neighboring Cluster');

subplot(1,3,3);
scatter_plot(mouseID, outSCorrNN_avg_tot, outFCorrNN_avg_tot, col_mat)
hold on; plot([0:0.1:1], [0:0.1:1], 'k'); xlabel('Correct Corr'); ylabel('Early Corr')
axis square
title('Not Neighboring Cluster');


% Compindx( ~any(Compindx,2), : ) = []; % rows
% Compindx( :, ~any(Compindx,1) ) = []; % columns
% figure;imagesc(Compindx)
% Compindx2( ~any(Compindx2,2), : ) = []; % rows
% Compindx2( :, ~any(Compindx2,1) ) = []; % columns
% figure;imagesc(Compindx2)

% SCorrindx( ~any(SCorrindx,2), : ) = []; % rows
% SCorrindx( :, ~any(SCorrindx,1) ) = []; % columns
% % SCorrindx(SCorrindx > 0.5) = 1;
% % SCorrindx(SCorrindx < 0 ) = 0;
% figure;imagesc(SCorrindx);
%
% FCorrindx( ~any(FCorrindx,2), : ) = []; % rows
% FCorrindx( :, ~any(FCorrindx,1) ) = []; % columns
% FCorrindx(FCorrindx > 0.5) = 1;
% FCorrindx(FCorrindx < 0 ) = 0;
% figure;imagesc(FCorrindx);
% k=1;
% figure;
% subplot(1,2,2); scatter(inSCorr{k}(:), inFCorr{k}(:), 4, 'MarkerEdgeColor', [0.5,0.5,0.5]); %axis([0 1 0 1])
%     hold on; plot([-1:0.1:1], [-1:0.1:1], 'k'); xlabel('Success Corr'); ylabel('Fail Corr')
%
%     clusters = reshape(indx, x, y);
%     subplot(1,2,1); imagesc(clusters); xlabel(['Cluster #', num2str(k)]);

% sc = [downsampled_movie(10,25,:);downsampled_movie(28,45,:);downsampled_movie(37,40,:);downsampled_movie(40,65,:);downsampled_movie(48,77,:);downsampled_movie(51,81,:)];