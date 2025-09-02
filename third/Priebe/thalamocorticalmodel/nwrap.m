rng(75001)
figure;
c50_list = [50];
sy_list = [2];
sx_list = [2];


nsamps = 1000;
thresh_list = [0.5];
PA_list = [120];

for ii = 1:length([phasen_list])
    subplot(1,2,ii)
    
    for i = 1:length(sy_list)
        sx = sx_list(1);
        PA = PA_list(1);
        thresh = thresh_list(1);
        sy = sy_list(i);
        c50 = c50_list(1);
        clear f1 f1v os osv ds Zp Zc Zpv Zcv ar
        nphases = 8;
        sv = 1;
        % sx = size_list(1)./(pi.*sy_list(i));
        % sx_list = [sx_list sx];
        %sx = sx_list(1);
        
        ar = zeros(1*nsamps,2);
        os = zeros(1,1*nsamps);
        osv = zeros(1,1*nsamps);
        f1 = zeros(1,1*nsamps);
        f1v = zeros(1,1*nsamps);
        ds = zeros(1,1*nsamps);
        Zp = zeros(1*nsamps,nphases);
        Zc = zeros(1*nsamps,nphases);
        Zpv = zeros(1*nsamps,nphases);
        Zpc = zeros(1*nsamps,nphases);
        sx = sx*sv; sy = sy *sv;
        for nc = [12]
            for sampnum = 1:nsamps
                [os(sampnum),osv(sampnum),ds(sampnum),Zp(sampnum,:),Zc(sampnum,:),Zpv(sampnum,:),Zcv(sampnum,:),ar(sampnum,:)] = lgnaggregatephgen2(sx,sy,nc,thresh,c50);
            end
        end
    
        for j=1:nsamps
          ff = fft(Zp(j,:)  - Zc(j,:));
          Zmod(j) = 2*abs(ff(2))/length(ff);
          Zmean(j) = ff(1)/length(ff);
          ff = fft(Zpv(j,:)  - Zcv(j,:));
          Zmodv(j) = 2*abs(ff(2))/length(ff);
          Zmeanv(j) = ff(1)/length(ff);
        end
    
        scatter(-Zmod,Zmean,'.')
        hold on
    end
    ylim([-4 2])
    xlim([-5 0])
    ylabel('Zmean')
    xlabel('Zmod')
    %AR_list = sy_list./sx_list(ii);
    leg = legend(num2str(thresh_list'),'Location','southwest');
    title(leg,'X-size')
    title(['Y - ' num2str(sy_list) '; X- ' num2str(sx_list) '; Plaid angle- ' num2str(PA_list)])

    
end
outpath = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\PriebeModel';
print(fullfile(outpath,'CompThresh_smrange.pdf'),'-dpdf','-fillpage')
save(fullfile(outpath,'Y2_Xpt5_PA120_threshpt3_conpt5'))

ars = sqrt(ar(:,2))./sqrt(ar(:,1));
area = pi.*ar(:,2).*ar(:,1);
figure;
subplot(2,2,1)
scatter(-Zmod,Zmean,[],ar(:,2))
colorbar
title('Y')
xlabel('Zmod')
ylabel('Zmean')
subplot(2,2,2)
scatter(-Zmod,Zmean,[],ar(:,1))
colorbar
title('X')
xlabel('Zmod')
ylabel('Zmean')
subplot(2,2,3)
scatter(-Zmod,Zmean,[],ars)
colorbar
title('AR')
xlabel('Zmod')
ylabel('Zmean')
subplot(2,2,4)
scatter(-Zmod,Zmean,[],area)
colorbar
title('Area')
xlabel('Zmod')
ylabel('Zmean')
sgtitle(['Y - ' num2str(sy_list) '; X- ' num2str(sx_list(ii)) '; Plaid angle- ' num2str(PA)])
print(fullfile(outpath,'CompXYARSize_PA120.pdf'),'-dpdf','-fillpage')

y_pred = fitlm([zscore(Zmean)' zscore(Zmod)'],zscore(ar(:,2)));
x_pred = fitlm([zscore(Zmean)' zscore(Zmod)'],zscore(ar(:,1)));
ar_pred = fitlm([zscore(Zmean)' zscore(Zmod)'],zscore(ars));
os_pred = fitlm([zscore(Zmean)' zscore(Zmod)'],zscore(os));
sz_pred = fitlm([zscore(Zmean)' zscore(Zmod)'],zscore(area));

figure;
subplot(2,2,1)
scatter(-Zmod,Zmean,[],ar(:,2))
colorbar
title('Y')
xlabel('Zmod')
ylabel('Zmean')
hold on
text(-15,0,['Mean: ' num2str(chop(y_pred.Coefficients.Estimate(2),2))])
text(-15,-1.25,['Mod: ' num2str(chop(y_pred.Coefficients.Estimate(3),2))])
subplot(2,2,2)
scatter(-Zmod,Zmean,[],area)
colorbar
title('Area')
xlabel('Zmod')
ylabel('Zmean')
hold on
text(-15,0,['Mean: ' num2str(chop(sz_pred.Coefficients.Estimate(2),2))])
text(-15,-1.25,['Mod: ' num2str(chop(sz_pred.Coefficients.Estimate(3),2))])
subplot(2,2,3)
scatter(-Zmod,Zmean,[],ars)
colorbar
title('AR')
xlabel('Zmod')
ylabel('Zmean')
hold on
text(-15,0,['Mean: ' num2str(chop(ar_pred.Coefficients.Estimate(2),2))])
text(-15,-1.25,['Mod: ' num2str(chop(ar_pred.Coefficients.Estimate(3),2))])
subplot(2,2,4)
scatter(-Zmod,Zmean,[],os)
colorbar
title('OS')
xlabel('Zmod')
ylabel('Zmean')
hold on
text(-15,0,['Mean: ' num2str(chop(os_pred_mean.Coefficients.Estimate(2),2))])
text(-15,-1.25,['Mod: ' num2str(chop(os_pred.Coefficients.Estimate(3),2))])
sgtitle(['Y - ' num2str(sy_list) '; X- ' num2str(sx_list(ii)) '; Plaid angle- ' num2str(PA)])
print(fullfile(outpath,'CompYARSzOS_wLM.pdf'),'-dpdf','-fillpage')


figure;
subplot(2,2,1)
scatter(-Zmod,Zmean,[],os)
colorbar
title('OS')
xlabel('Zmod')
ylabel('Zmean')

figure
subplot(2,2,1)
scatter(os,ar(:,2))
xlabel('OS')
ylabel('Y')
title(num2str(chop(triu2vec(corrcoef(os,ar(:,2))),2)))
subplot(2,2,2)
scatter(osv,ar(:,2))
xlabel('OSV')
ylabel('Y')
title(num2str(chop(triu2vec(corrcoef(osv,ar(:,2))),2)))
subplot(2,2,3)
scatter(os,ars)
xlabel('OS')
ylabel('AR')
title(num2str(chop(triu2vec(corrcoef(os,ars)),2)))
subplot(2,2,4)
scatter(os,area)
xlabel('OS')
ylabel('Area')
title(num2str(chop(triu2vec(corrcoef(os,area)),2)))
sgtitle(['Y - ' num2str(sy_list) '; X- ' num2str(sx_list(ii)) '; Plaid angle- ' num2str(PA)])
print(fullfile(outpath,'OSvY.pdf'),'-dpdf','-fillpage')

os_pred = fitlm([zscore(ar(:,2)) zscore(ar(:,1))],zscore(os));
figure; scatter(ar(:,1),ar(:,2),[],os); 


[nx edgesx binx] = histcounts(ar_all(:,1),10);
[ny edgesy biny] = histcounts(ar_all(:,2),10);
mean_mat = zeros(10,10);
mod_mat = zeros(10,10);
for i = 1:10
    indx = find(binx == i);
    for ii = 1:10
        indy = find(biny == ii);
        indu = intersect(indx,indy);
        mean_mat(ii,i) = mean(Zmean_all(indu));
        mod_mat(ii,i) = mean(Zmod_all(indu));
    end
end


figure; 
imagesc(flipud(mean_mat),'XData',.5, 'YData', .5);
suptitle('Zmean')

figure; 
for i = 1:9
    subplot(3,3,i)
    imagesc(flipud(mod_mat(:,:,i)),'XData',.5, 'YData', .5);
end
suptitle('Zmod')


[nx edgesmod binmod] = histcounts(Zmod_all,10);
[ny edgesmean binmean] = histcounts(Zmean_all,10);
for i = 1:10
    indx = find(binmod == i);
    for ii = 1:10
        indy = find(binmean == ii);
        indu = intersect(indx,indy);
        ary_mat(ii,i) = mean(ar_all(indu,2),1);
        ar_mat(ii,i) = mean(ars(indu));
    end
end
figure; 
imagesc(flipud(ary_mat),'XData',.5, 'YData', .5);
set(gca,'Xtick',0:2:10,'XTickLabels', edgesmod(1:2:end),'Ytick',0:2:10,'YTickLabels', fliplr(edgesmean(1:2:end)))
colorbar
suptitle('Y')

figure; 
imagesc(flipud(ar_mat),'XData',.5, 'YData', .5);
set(gca,'Xtick',0:2:10,'XTickLabels', edgesmod(1:2:end),'Ytick',0:2:10,'YTickLabels', fliplr(edgesmean(1:2:end)))
colorbar
suptitle('AR')