function [fig,m] = FigShotNoise(array,clim);
%[fig,m] = FigShotNoise(array,clim);

up = 99;

darray = double(array);
av = mean(darray,3);
avmax = prctile(av(:),up);
va = var(darray,0,3);
vamax = prctile(va(:),up);
f = inline('m*x','m','x');
m = lsqcurvefit(f,1,av(:),va(:));

avnorm = av(:)/m;
vanorm = va(:)/m^2;

avnormmax = prctile(avnorm,up);
vanormmax = prctile(vanorm,up);

% images
fig = figure
subplot(2,3,1);
imagesc(darray(:,:,1));axis square;colormap gray;title('Single frame');
subplot(2,3,4);
imagesc(mean(darray,3));axis square;colormap gray;title('Average frame');

% scatter plots
subplot(2,3,2);
plot(avnorm(:),vanorm(:),'k.','markersize',1);
hold on
plot([0 avmax]/m,f(m,[0 avmax])/m^2,'r','linewidth',2);
%xlim([0 avmax]/m);ylim([0 vamax]/m^2);

axis square;grid;
xlabel('Mean (counts)');
ylabel('Variance');
set(gca,'xtick',[00:10:100]);

snr = sqrt(av(:).^2./va(:));
snrmax = max(snr);

subplot(2,3,5);
plot(avnorm,snr,'k.','markersize',1);
axis square;grid;
ylabel('SNR');
xlabel('Mean (counts)');
xlim([0 avnormmax]);ylim([0 snrmax]);
set(fig,'paperpositionmode','auto');

set(gca,'xtick',[00:10:100]);

dx = 1;
au = [0:.5:avnormmax];
nn = hist(avnorm,au);
density = nn/sum(nn);
cumulative  = cumsum(density);

% cumulative distributions
subplot(2,3,3);
plot(au,cumulative,'k');
xlabel('FOV Photon counts/ pixel (.3\mus)')
ylabel('CDF');
axis square; box off;
xlim([0 1]);
ylim([0 1.1]);

dx = .1;
au = [0:dx:15];
nn = hist(sqrt(av(:).^2./va(:)),au);
density = nn/sum(nn);
cumulative  = cumsum(density);

set(gca,'xtick',[00:10:100]);
grid

subplot(2,3,6);
plot(au,cumulative,'k');
xlabel('Pixel SNR')
ylabel('CDF');
axis square; box off;
xlim([0 2]);
ylim([0 1.1]);
set(gca,'xtick',[0:2:20]);

grid
return;