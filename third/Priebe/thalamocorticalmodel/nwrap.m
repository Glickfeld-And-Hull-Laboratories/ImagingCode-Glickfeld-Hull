rng(75001)
figure;
c50_list = [20 50 80];
sy_list = [1:.5:2];
sx_list = 1./(pi.*sy_list);
nsamps = 1000;
thresh = 0.2;

for i = 1:length(sy_list)
    c50 = c50_list(i);
    clear f1 f1v os osv ds Zp Zc Zpv Zcv ar
    nphases = 8;
    sv = 1;
    sx = sy_list(2);
    sy = sy_list(2);
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


offset  = 200;
figure(504)
subplot(1,2,1)

tt= (copper(100));
for j=1:nsamps
    if ~isnan(os(j))
%        li = plot(-Zmodb(j),Zmodb(j)+Zmeanb(j),'.');
        li = plot(-Zmod(j),Zmean(j),'.');
        set(li,'Color',tt(1+floor(os(j)*0.99*100),:))
        hold on
    end
end
hold off



subplot(1,2,2)
ars = sqrt(ar(:,2))./sqrt(ar(:,1));
tt= (copper(100));
for j=1:nsamps
    if ~isnan(os(j))
%        li = plot(-Zmodb(j),Zmodb(j)+Zmeanb(j),'.');
        li = plot(-Zmod(j),Zmean(j),'.');
        set(li,'Color',tt(1+floor(ars(j)*(0.99*8)),:));
        hold on
    end
end
hold off




positions = [normrnd(5,0.01,8,1) normrnd(0,1,8,1)+(0:1.25:9)'-5 ];
clear ar

[os,osv,ds,Zp,Zc,Zpv,Zcv,ar(1:2)] = lgnaggregatephgen2(sx,sy,nc,thresh,positions);
% 
% 
% 
% positions = [normrnd(5,0.01,8,1) normrnd(0,1,8,1)+(0:1.25:9)'-5 ];
% clear ar
% 
% [os,osv,ds,Zp,Zc,Zpv,Zcv,ar(1:2)] = lgnaggregatephgen2(sx,sy,nc,thresh,positions);
