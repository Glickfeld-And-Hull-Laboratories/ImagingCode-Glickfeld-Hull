f1amp = 1;
f2amp = 1;
tau = 120;
t = 1:120;
g = 0.05*exp(-t/tau);
for offset = [0 1 2 4]
cnt1 = 1;
for f1amp = 0.05:0.05:5
f1ph= 0;
f2ph = 0;

step = pi/60;
stepf2 = 2*step;

cnt = 1;
vv = zeros(13,120);
ss = zeros(13,120);
cc = zeros(13,120);
for f1ph = 0:pi/6:(2*pi-pi/6)
v1 = f1amp*(sin( (step:step:(2*pi))+f1ph )) + f2amp* ( sin((stepf2:stepf2:(4*pi))+f2ph)) + offset;

s1 = rect(v1).^3;
ca = conv(s1,g);

ff = fft(v1);
vf1(cnt) = abs(2*ff(2)/length(ff));
vf1ang(cnt) = angle(ff(2));
vf2(cnt) = abs(2*ff(3)/length(ff));
vf2ang(cnt) = angle(ff(3));

ff = fft(s1);
sf1(cnt) = abs(2*ff(2)/length(ff));
sf1ang(cnt) = angle(ff(2));
sf2(cnt) = abs(2*ff(3)/length(ff));
sf2ang(cnt) = angle(ff(3));

ff = fft(ca);
cf1(cnt) = abs(2*ff(2)/length(ff));
cf1ang(cnt) = angle(ff(2));
cf2(cnt) = abs(2*ff(3)/length(ff));
cf2ang(cnt) = angle(ff(3));

cnt = cnt + 1;
vv(cnt,:) = v1;
ss(cnt,:) = s1;
cc(cnt,:) = ca(1:120);
end

[vf1comp(cnt1), vf2comp(cnt1) ] = compcontrastrevf1f2_12(vv);
[sf1comp(cnt1), sf2comp(cnt1) ] = compcontrastrevf1f2_12(ss);
[cf1comp(cnt1), cf2comp(cnt1) ] = compcontrastrevf1f2_12(cc);

% Now compute the predicted F1/F0 components.  Assume F2 amp is the same as
% F0
vresp = f2amp + f1amp*sin((step:step:(2*pi)));
sresp = rect(vresp).^3;
ff = fft(sresp);

vvf1(cnt1) = f1amp;
vvf0(cnt1) = f2amp;
ssf1(cnt1) = 2*abs(ff(2))/length(ff);
ssf0(cnt1) = ff(1)/length(ff);


caconv = conv(sresp,g);
caff = fft(caconv);
caf1(cnt1) = 2*abs(caff(2))/length(caff);
caf0(cnt1) = caff(1)/length(caff);


cnt1 = cnt1 + 1;
end
figure(100)
subplot(2,2,2)
plot(vf2comp./vf1comp,sf2comp./sf1comp,'x-')
hold on
subplot(2,2,4)
plot(sf2comp./sf1comp,cf2comp./cf1comp,'x-')
hold on

%hold on
%cnt1
cnt1 = 1;
end
%pause
subplot(2,2,1)
plot(vvf1./vvf0,ssf1./ssf0,'x-')
subplot(2,2,3)
plot(ssf1./ssf0,caf1./caf0,'x-')
subplot(2,2,2)
xlabel('VF2 / VF1')
ylabel('SF2 / SF1')
set(gca,'XScale','log','YScale','log')
legend('0','1','2','4')
hold off
subplot(2,2,4)
xlabel('SF2 / SF1')
ylabel('CF2 / CF1')
set(gca,'XScale','log','YScale','log')
legend('0','1','2','4')
hold off
subplot(2,2,1)
xlabel('VF1 / VF0')
ylabel('SF1 / SF0')
xlim([0.02 10])
ylim([0.02 10])
set(gca,'XScale','log','YScale','log')
subplot(2,2,3)
xlabel('SF1 / SF0')
ylabel('CF1 / CF0')
xlim([0.02 10])
ylim([0.02 10])
set(gca,'XScale','log','YScale','log')

