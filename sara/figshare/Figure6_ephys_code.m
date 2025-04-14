%%% This will make the plots from figure 6              %%%

%%% contains the variable r which contains the following:%%
%%% r.vtest - membrane potential response amplitude to  %%%
%%%           the test stimulus at 48% contrast for each%%%
%%%           neuron                                    %%%
%%% r.vmask - same as r.vtest for the mask stimulus     %%%
%%% r.vplaid -Single cycle responses to plaids at 4     %%%
%%%           mask spatial phases: 0, 90, 180 and 270   %%%
%%% r.stest, r.smask,r.splaid and iplaid are the same for%%%
%%%           spike rates and current. itest and imask are single   %%%
%%%           cycle amplitudes                          %%%

clear all; close all

cd '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara\Data\fromNicholas\figshare'
load Figure6_ephys_data.mat 

%%% we need to calculate the phase modulation baseline and amplitude  %%%
%%% for our spike rates. To do this we fit a sine wave to the masking %%%
%%% index at each plaid phase using a bootstrap method                %%%

for k = 1:length(r.vtest)
    nboot = 5000;
    x = [0 90 180 270];
    x = deg2rad(x);
    vplaid = r.vplaid{k};
    vtest = r.vtest(k);
    vmask = r.vmask(k);
    vplaid = vplaid';
    si_v(k) = (vtest-vmask)/(vtest+vmask);
    si_v(k) = abs(si_v(k));
    ln = length(vplaid(:,1));

for j = 1:4
    bootsmp(:,j) = bootstrp(nboot,@mean,vplaid(:,j));
end

bootsmp = rect(bootsmp);

 for j = 1:length(bootsmp(:,1))
     for jj = 1:4
         mi(j,jj) = (bootsmp(j,jj)-(vmask+vtest))/(bootsmp(j,jj)+(vmask+vtest));
     end
 end
 mi_vm{k} = mi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zz = mean(mi);

xx = zeros(size(mi));
for j = 1:length(xx)
    xx(j,:) = x;
end

xx = xx(:);
y = mi(:);

amp_guess = (max(zz)-min(zz))/2; %%estimate the amplitude
off_guess = mean(zz);    
per_guess=2*pi;
tt = fft(zz);
pha_guess = angle(tt(2)); 

if off_guess>0
    off_lb = off_guess.*.5;
    off_ub = off_guess.*2;
else
    off_ub = off_guess.*.5;
    off_lb = off_guess.*2;
end

amp_lb = 0;
amp_ub = amp_guess.*2.5;
per_lb = (2*pi)-.0001; 
per_ub = 2*pi+.0001; 
pha_lb = -2*pi; 
pha_ub = 2*pi;  

[xData, yData] = prepareCurveData(xx,y);

% Set up fittype and options.
ft = fittype( 'off+amp.*(sin(2*pi*x./per + 2.*pi/pha))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Lower = [amp_lb off_lb per_lb pha_lb];  
opts.Upper = [amp_ub off_ub per_ub pha_ub];  
opts.StartPoint = [amp_guess off_guess per_guess pha_guess];

[fitresult,gof] = fit( xData, yData, ft, opts);
coeffs = coeffvalues(fitresult);
rsqr = gof.rsquare;
xp = linspace(min(x),max(x));
ftres{k} = feval(fitresult,xp);
clear fitresult gof opts ft xData yData

bb = fft(ftres{k});
sinf1(k) = 2*abs(bb(2)/100);
sindc(k) = bb(1)/100;
clearvars -except n3 n4 n1 n2 sinf1 sindc ftres mi_vm r si_v
end


%%%  now we have to do it all again for spike rate  %%%

for k = 1:length(r.stest)
    nboot = 5000;
    x = [0 90 180 270];
    x = deg2rad(x);
    splaid = r.splaid{k};
    stest = r.stest(k);
    smask = r.smask(k);
    splaid = splaid';
    si_spk(k) = (stest-smask)/(stest+smask);
    si_spk(k) = abs(si_spk(k));
    ln = length(splaid(:,1));

for j = 1:4
    bootsmp(:,j) = bootstrp(nboot,@mean,splaid(:,j));
end

bootsmp = rect(bootsmp);

 for j = 1:length(bootsmp(:,1))
     for jj = 1:4
         mi(j,jj) = (bootsmp(j,jj)-(smask+stest))/(bootsmp(j,jj)+(smask+stest));
     end
 end
 mi_spk{k} = mi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zz = mean(mi);

xx = zeros(size(mi));
for j = 1:length(xx)
    xx(j,:) = x;
end

xx = xx(:);
y = mi(:);

amp_guess = (max(zz)-min(zz))/2; %%estimate the amplitude
off_guess = mean(zz);   
per_guess=2*pi;
tt = fft(zz);
pha_guess = angle(tt(2)); 

if off_guess>0
    off_lb = off_guess.*.5;
    off_ub = off_guess.*2;
else
    off_ub = off_guess.*.5;
    off_lb = off_guess.*2;
end

amp_lb = 0;
amp_ub = amp_guess.*2.5;
per_lb = (2*pi)-.0001; 
per_ub = 2*pi+.0001; 
pha_lb = -2*pi; 
pha_ub = 2*pi; 

[xData, yData] = prepareCurveData(xx,y);

% Set up fittype and options.
ft = fittype( 'off+amp.*(sin(2*pi*x./per + 2.*pi/pha))', 'independent', 'x', 'dependent', 'y' );


  [b_hat_all(iCell,1), amp_hat_all(k,1), per_hat_all(k,1),pha_hat_all(k,1),sse_all(k,1),R_square_all(k,1)] = sinefit_PCI(deg2rad(phase),PCI(:,k));
   yfit_all(iCell,:,1) = b_hat_all(k,1)+amp_hat_all(k,1).*(sin(2*pi*deg2rad(phase_range)./per_hat_all(k,1) + 2.*pi/pha_hat_all(k,1)));

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Lower = [amp_lb off_lb per_lb pha_lb];  
opts.Upper = [amp_ub off_ub per_ub pha_ub];  
opts.StartPoint = [amp_guess off_guess per_guess pha_guess];

[fitresult,gof] = fit( xData, yData, ft, opts);
coeffs = coeffvalues(fitresult);
rsqr = gof.rsquare;
xp = linspace(min(x),max(x));
ftres_spk{k} = feval(fitresult,xp);
clear fitresult gof opts ft xData yData

bb = fft(ftres_spk{k});
sinf1_spk(k) = 2*abs(bb(2)/100);
sindc_spk(k) = bb(1)/100;
clearvars -except n1 n2 n4 n3 sinf1 sindc ftres mi_vm r si_v sinf1_spk sindc_spk mi_spk ftres_spk si_spk
end



%% now we have to do it all again again for current with opto  %%%

for k = 1:6
    nboot = 5000;
    x = [0 90 180 270];
    x = deg2rad(x);
    iplaid = r.iplaid{k};
    itest = mean(r.itest{k});
    imask = mean(r.imask{k});
    iplaid = iplaid';
    si_i(k) = (itest-imask)/(itest+imask);
    si_i(k) = abs(si_i(k));
    ln = length(iplaid(:,1));

for j = 1:4
    bootsmp(:,j) = bootstrp(nboot,@mean,iplaid(:,j));
end

bootsmp = rect(bootsmp);

 for j = 1:length(bootsmp(:,1))
     for jj = 1:4
         mi(j,jj) = (bootsmp(j,jj)-(imask+itest))/(bootsmp(j,jj)+(imask+itest));
     end
 end
 mi_i{k} = mi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zz = mean(mi);

xx = zeros(size(mi));
for j = 1:length(xx)
    xx(j,:) = x;
end

xx = xx(:);
y = mi(:);

amp_guess = (max(zz)-min(zz))/2; %%estimate the amplitude
off_guess = mean(zz);   
per_guess=2*pi;
tt = fft(zz);
pha_guess = angle(tt(2)); 

if off_guess>0
    off_lb = off_guess.*.5;
    off_ub = off_guess.*2;
else
    off_ub = off_guess.*.5;
    off_lb = off_guess.*2;
end

amp_lb = 0;
amp_ub = amp_guess.*2.5;
per_lb = (2*pi)-.0001; 
per_ub = 2*pi+.0001; 
pha_lb = -2*pi; 
pha_ub = 2*pi; 

[xData, yData] = prepareCurveData(xx,y);

% % Set up fittype and options.
% ft = fittype( 'off+amp.*(sin(2*pi*x./per + 2.*pi/pha))', 'independent', 'x', 'dependent', 'y' );
% 
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Lower = [amp_lb off_lb per_lb pha_lb];  
% opts.Upper = [amp_ub off_ub per_ub pha_ub];  
% opts.StartPoint = [amp_guess off_guess per_guess pha_guess];
% 
% [fitresult,gof] = fit( xData, yData, ft, opts);
% coeffs = coeffvalues(fitresult);
% rsqr = gof.rsquare;
% xp = linspace(min(x),max(x));
% ftres_i{k} = feval(fitresult,xp);
% clear fitresult gof opts ft xData yData
% 
% bb = fft(ftres_i{k});
% sinf1_i(k) = 2*abs(bb(2)/100);
% sindc_i(k) = bb(1)/100;
% clearvars -except n3 n4 sinf1_i sindc_i ftres_i mi_i n1 n2 sinf1 sindc ftres mi_vm r si_v sinf1_spk sindc_spk mi_spk ftres_spk si_spk
% end



%% now plot the example cells                                   %%%%

figure %example cell one
subplot(2,7,1)
plot(n1.vm(1,:),'-k'); hold on; plot(n1.vm(2,:),'-r');set(gca,'TickDir','out'); box off;  axis([0 100 -65 -40]); xticks([0 100]); xticklabels([0 500]); yticks([-65 -40]); xlabel('time (ms)'); ylabel('mV');
subplot(2,7,2)
plot(n1.vm(3,:),'-k'); set(gca,'TickDir','out'); box off; axis([0 100 -65 -40]); xticks([0 100]); xticklabels([0 500]); yticks([-65 -40]); xlabel('time (ms)'); ylabel('mV');
subplot(2,7,3)
plot(n1.vm(4,:),'-k'); set(gca,'TickDir','out'); box off; axis([0 100 -65 -40]); xticks([0 100]); xticklabels([0 500]); yticks([-65 -40]); xlabel('time (ms)'); ylabel('mV');
subplot(2,7,4)
plot(n1.vm(5,:),'-k'); set(gca,'TickDir','out'); box off; axis([0 100 -65 -40]); xticks([0 100]); xticklabels([0 500]); yticks([-65 -40]); xlabel('time (ms)'); ylabel('mV');
subplot(2,7,5)
plot(n1.vm(6,:),'-k'); set(gca,'TickDir','out'); box off; axis([0 100 -65 -40]); xticks([0 100]); xticklabels([0 500]); yticks([-65 -40]); xlabel('time (ms)'); ylabel('mV');
subplot(2,7,6:7)
x = [0,90,180,270];
% errorbar(x,mean(mi_vm{8}),std(mi_vm{8}),'o')
% xp = linspace(min(x),max(x));
% hold on; plot(xp,ftres{8},'-k')
axis([-15 285 -1 1]); xticks([0 90 180 270]); yticks([-1 0 1]); xlabel('mask phase (^o)'); ylabel('MI'); set(gca,'TickDir','out'); box off;
subplot(2,7,8)
bar(n1.spk(1,:),1,'FaceColor','k','EdgeAlpha',[0]); hold on; bar(n1.spk(2,:),1,'FaceColor','r','EdgeAlpha',[0]); set(gca,'TickDir','out'); box off; axis([0 20 0 16]); xticks([0 20]); xticklabels([0 16]); yticks([0 16]); xlabel('time (ms)'); ylabel('spks/s');
subplot(2,7,9)
bar(n1.spk(3,:),1,'FaceColor','k','EdgeAlpha',[0]); set(gca,'TickDir','out'); box off; axis([0 20 0 16]); xticks([0 20]); xticklabels([0 500]); yticks([0 16]);xlabel('time (ms)'); ylabel('spks/s');
subplot(2,7,10)
bar(n1.spk(4,:),1,'FaceColor','k','EdgeAlpha',[0]); set(gca,'TickDir','out'); box off; axis([0 20 0 16]); xticks([0 20]); xticklabels([0 500]); yticks([0 16]);xlabel('time (ms)'); ylabel('spks/s');
subplot(2,7,11)
bar(n1.spk(5,:),1,'FaceColor','k','EdgeAlpha',[0]); set(gca,'TickDir','out'); box off; axis([0 20 0 16]); xticks([0 20]); xticklabels([0 500]); yticks([0 16]);xlabel('time (ms)'); ylabel('spks/s');
subplot(2,7,12)
bar(n1.spk(6,:),1,'FaceColor','k','EdgeAlpha',[0]); set(gca,'TickDir','out'); box off; axis([0 20 0 16]); xticks([0 20]); xticklabels([0 500]); yticks([0 16]);xlabel('time (ms)'); ylabel('spks/s');
subplot(2,7,13:14)
% errorbar(x,mean(mi_spk{6}),std(mi_spk{6}),'o')
% xp = linspace(min(x),max(x));
% hold on; plot(xp,ftres_spk{6},'-k')
axis([-15 285 -1 1]); xticks([0 90 180 270]); yticks([-1 0 1]); xlabel('mask phase (^o)'); ylabel('MI'); set(gca,'TickDir','out'); box off;
set(gcf,'color','w');


%% now example cell two %%
figure
subplot(2,7,1)
plot(n2.vm(1,:),'-k'); hold on; plot(n2.vm(2,:),'-r');set(gca,'TickDir','out'); box off;  axis([0 100 -70 -40]); xticks([0 100]); xticklabels([0 500]); yticks([-70 -40]); xlabel('time (ms)'); ylabel('mV');
subplot(2,7,2)
plot(n2.vm(3,:),'-k'); set(gca,'TickDir','out'); box off; axis([0 100 -70 -40]); xticks([0 100]); xticklabels([0 500]); yticks([-70 -40]); xlabel('time (ms)'); ylabel('mV');
subplot(2,7,3)
plot(n2.vm(4,:),'-k'); set(gca,'TickDir','out'); box off; axis([0 100 -70 -40]); xticks([0 100]); xticklabels([0 500]); yticks([-70 -40]); xlabel('time (ms)'); ylabel('mV');
subplot(2,7,4)
plot(n2.vm(5,:),'-k'); set(gca,'TickDir','out'); box off; axis([0 100 -70 -40]); xticks([0 100]); xticklabels([0 500]); yticks([-70 -40]); xlabel('time (ms)'); ylabel('mV');
subplot(2,7,5)
plot(n2.vm(6,:),'-k'); set(gca,'TickDir','out'); box off; axis([0 100 -70 -40]); xticks([0 100]); xticklabels([0 500]); yticks([-70 -40]); xlabel('time (ms)'); ylabel('mV');
subplot(2,7,6:7)
x = [0,90,180,270];
errorbar(x,mean(mi_vm{17}),std(mi_vm{17}),'o')
xp = linspace(min(x),max(x));
hold on; plot(xp,ftres{17},'-k')
axis([-15 285 -1 1]); xticks([0 90 180 270]); yticks([-1 0 1]); xlabel('mask phase (^o)'); ylabel('MI'); set(gca,'TickDir','out'); box off;
subplot(2,7,8)
bar(n2.spk(1,:),1,'FaceColor','k','EdgeAlpha',[0]); hold on; bar(n2.spk(2,:),1,'FaceColor','r','EdgeAlpha',[0]); set(gca,'TickDir','out'); box off; axis([0 20 0 16]); xticks([0 20]); xticklabels([0 16]); yticks([0 16]); xlabel('time (ms)'); ylabel('spks/s');
subplot(2,7,9)
bar(n2.spk(3,:),1,'FaceColor','k','EdgeAlpha',[0]); set(gca,'TickDir','out'); box off; axis([0 20 0 16]); xticks([0 20]); xticklabels([0 500]); yticks([0 16]);xlabel('time (ms)'); ylabel('spks/s');
subplot(2,7,10)
bar(n2.spk(4,:),1,'FaceColor','k','EdgeAlpha',[0]); set(gca,'TickDir','out'); box off; axis([0 20 0 16]); xticks([0 20]); xticklabels([0 500]); yticks([0 16]);xlabel('time (ms)'); ylabel('spks/s');
subplot(2,7,11)
bar(n2.spk(5,:),1,'FaceColor','k','EdgeAlpha',[0]); set(gca,'TickDir','out'); box off; axis([0 20 0 16]); xticks([0 20]); xticklabels([0 500]); yticks([0 16]);xlabel('time (ms)'); ylabel('spks/s');
subplot(2,7,12)
bar(n2.spk(6,:),1,'FaceColor','k','EdgeAlpha',[0]); set(gca,'TickDir','out'); box off; axis([0 20 0 16]); xticks([0 20]); xticklabels([0 500]); yticks([0 16]);xlabel('time (ms)'); ylabel('spks/s');
subplot(2,7,13:14)
errorbar(x,mean(mi_spk{15}),std(mi_spk{15}),'o')
xp = linspace(min(x),max(x));
hold on; plot(xp,ftres_spk{15},'-k')
axis([-15 285 -1 1]); xticks([0 90 180 270]); yticks([-1 0 1]); xlabel('mask phase (^o)'); ylabel('MI'); set(gca,'TickDir','out'); box off;
set(gcf,'color','w');


%% now example cell three and four together (both are opto current) %%
figure
subplot(2,7,1)
plot(n3.curr(1,:)); hold on; plot(n3.curr(2,:)); set(gca,'TickDir','out'); box off; axis([0 100 40 70]); ylabel('current (pA)'); xticks([0 100]); xticklabels([0 500]); xlabel('time (ms)')
subplot(2,7,2)
plot(n3.curr(3,:)); set(gca,'TickDir','out'); box off; axis([0 100 10 40]); xticks([0 100]); xticklabels([0 500]); xlabel('time (ms)')
subplot(2,7,3)
plot(n3.curr(4,:)); set(gca,'TickDir','out'); box off; axis([0 100 10 40]);xticks([0 100]); xticklabels([0 500]); xlabel('time (ms)')
subplot(2,7,4)
plot(n3.curr(5,:)); set(gca,'TickDir','out'); box off; axis([0 100 10 40]);xticks([0 100]); xticklabels([0 500]); xlabel('time (ms)')
subplot(2,7,5)
plot(n3.curr(6,:)); set(gca,'TickDir','out'); box off; axis([0 100 10 40]);xticks([0 100]); xticklabels([0 500]); xlabel('time (ms)')
subplot(2,7,6:7)
% errorbar(x,mean(mi_i{6}),std(mi_i{6}),'o')
% xp = linspace(min(x),max(x));
% hold on; plot(xp,ftres_i{6},'-k')
axis([-15 285 -1 1]); xticks([0 90 180 270]); yticks([-1 0 1]); xlabel('mask phase (^o)'); ylabel('MI'); set(gca,'TickDir','out'); box off;
set(gcf,'color','w');

subplot(2,7,8)
plot(n4.curr(1,:)); hold on; plot(n3.curr(2,:)); axis([0 100 40 70]); set(gca,'TickDir','out'); box off;
subplot(2,7,9)
plot(n4.curr(3,:)); set(gca,'TickDir','out'); box off; axis([0 100 10 40]); axis([0 100 35 65]);  xticks([0 100]); xticklabels([0 500]); xlabel('time (ms)')
subplot(2,7,10)
plot(n4.curr(4,:));set(gca,'TickDir','out'); box off; axis([0 100 10 40]);axis([0 100 35 65]);  xticks([0 100]); xticklabels([0 500]); xlabel('time (ms)')
subplot(2,7,11)
plot(n4.curr(5,:));set(gca,'TickDir','out'); box off; axis([0 100 10 40]);axis([0 100 35 65]);  xticks([0 100]); xticklabels([0 500]); xlabel('time (ms)')
subplot(2,7,12)
plot(n4.curr(6,:));set(gca,'TickDir','out'); box off; axis([0 100 10 40]);axis([0 100 35 65]);  xticks([0 100]); xticklabels([0 500]); xlabel('time (ms)')
subplot(2,7,13:14)
% errorbar(x,mean(mi_i{3}),std(mi_i{3}),'o')
% xp = linspace(min(x),max(x));
% hold on; plot(xp,ftres_i{3},'-k')
axis([-15 285 -1 1]); xticks([0 90 180 270]); yticks([-1 0 1]); xlabel('mask phase (^o)'); ylabel('MI'); set(gca,'TickDir','out'); box off;
set(gcf,'color','w');

%%%%now we can make the histograms of sin amplitude and baseline %%%%
figure
subplot(221)
histogram(sinf1,'BinWidth',.16); axis([0 1 0 20]); xlabel('Phase modulation amplitude'); ylabel('number of neurons'); set(gca,'TickDir','out'); box off; xticks([0 .5 1]); yticks([0 10 20]);
title('membrane potential')
subplot(222)
histogram(sinf1_spk,'BinWidth',.16); axis([0 1 0 20]); xlabel('Phase modulation amplitude'); ylabel('number of neurons'); set(gca,'TickDir','out'); box off; xticks([0 .5 1]); yticks([0 10 20]);
title('spike rate')
subplot(223)
histogram(sindc,'BinWidth',.19); axis([-1 1 0 12]); xlabel('Phase modulation amplitude'); ylabel('number of neurons'); set(gca,'TickDir','out'); box off; xticks([0 .5 1]); yticks([0 10 20]);
title('membrane potential')
subplot(224)
histogram(sindc_spk,'BinWidth',.19); axis([-1 1 0 10]); xlabel('Phase modulation amplitude'); ylabel('number of neurons'); set(gca,'TickDir','out'); box off; xticks([0 .5 1]); yticks([0 10 20]);
title('spike rate')
set(gcf,'color','w');