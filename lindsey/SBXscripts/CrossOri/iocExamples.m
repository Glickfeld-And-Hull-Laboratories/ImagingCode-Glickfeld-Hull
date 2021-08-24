ang1 = 7*pi/8;
sp1 = 10;
ang2 = 0;
sp2 = 30;
[ioc_ang ioc_sp av_ang av_sp] = iocVectorCalc(ang1,sp1,ang2,sp2);
refline(1)
rad2deg(ioc_ang)
rad2deg(av_ang)

ang1 = -3*pi/8;
sp1 = 10;
ang2 = pi/2;
sp2 = 30;
[ioc_ang ioc_sp av_ang av_sp] = iocVectorCalc(ang1,sp1,ang2,sp2);
refline(1)
rad2deg(ioc_ang)
rad2deg(av_ang)