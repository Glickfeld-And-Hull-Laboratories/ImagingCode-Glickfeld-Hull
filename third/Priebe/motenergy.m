%symunit = g2d(1,30,1,30/8,30/4,10,d2r(0),[15 15],-r2d(atan(slope)));
%asymunit = g2d(1,30,1,30/8,30/4,10,d2r(90),[15 15],-r2d(atan(slope)));
%tfuncs = sin(0:2*pi/29:2*pi);
%tfuncl = sin(0:2*pi/50:2*pi); tfuncl = tfuncl(1:30);

si = 151;
width =0.3*si;
asym = diff(gengauss(width,si,si/2))*3000;
sym = diff(diff(gengauss(width,si,si/2)))*3000*3000;
asym = asym(1:149);
asym = asym / sum(sum(abs(asym)));
sym = sym / sum(sum(abs(sym)));
%asym = rot(sym,-2);
%asym = -sym;

%t = 0:0.5:15;
t = [0:(15/150):15];
n = 3;
tfuncs = t.^n .*exp(-t) .* (1/(factorial(n)) - (t.^2)/factorial(n+2)) * 7;
n = 5;
tfuncl = t.^n .*exp(-t) .* (1/(factorial(n)) - (t.^2)/factorial(n+2)) * 7;
symunits = tfuncs'*sym;
symunitl = tfuncl'*sym;
asymunits = tfuncs'*asym;
asymunitl = tfuncl'*asym;

subplot(3,4,9)
imagesc(symunits);caxis([-1 1]*0.02);
colormap gray
subplot(3,4,10)
imagesc(symunitl);caxis([-1 1]*0.02);
colormap gray
subplot(3,4,11)
imagesc(symunits);caxis([-1 1]*0.02);
colormap gray
subplot(3,4,12)
imagesc(symunitl);caxis([-1 1]*0.02);
colormap gray

subplot(3,4,5)
imagesc(asymunitl);caxis([-1 1]*0.03);
colormap gray
subplot(3,4,6)
imagesc(asymunits);caxis([-1 1]*0.03);
colormap gray
subplot(3,4,7)
imagesc(asymunitl);caxis([-1 1]*0.03);
colormap gray
subplot(3,4,8)
imagesc(asymunits);caxis([-1 1]*0.03);
colormap gray


subplot(3,4,1)
imagesc(symunits-asymunitl);caxis([-1 1]*0.03);
colormap gray
subplot(3,4,2)
imagesc(symunitl+asymunits);caxis([-1 1]*0.03);
colormap gray
subplot(3,4,3)
imagesc(symunits + asymunitl);caxis([-1 1]*0.03);
colormap gray
subplot(3,4,4)
imagesc(symunitl - asymunits);caxis([-1 1]*0.03);
colormap gray

print('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Grants\CrossOri_R01\BrainR01_2022\motenergy.pdf','-dpdf')

fts =  fft(tfuncs,32);
ftl =  fft(tfuncl,32);
Pfts = fts.*conj(fts)/ 32;
Pftl = ftl.*conj(ftl)/ 32;

leftneg= symunits-asymunitl;
leftpos = symunitl+asymunits;
rightpos = symunits + asymunitl;
rightneg = symunitl - asymunits;

