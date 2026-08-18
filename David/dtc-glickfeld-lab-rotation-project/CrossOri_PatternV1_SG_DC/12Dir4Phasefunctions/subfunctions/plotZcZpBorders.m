plot(1.28.*ones(1,length(-10:1:0)),-10:1:0,'--k')
hold on
plot(-10:1:0,1.28.*ones(1,length(-10:1:0)),'--k')
xp = [0:10];
xc = [1.28:10];
yp = xp+1.28;
yc = xc-1.28;
plot(xp,yp,'--k')
plot(xc,yc,'--k')
