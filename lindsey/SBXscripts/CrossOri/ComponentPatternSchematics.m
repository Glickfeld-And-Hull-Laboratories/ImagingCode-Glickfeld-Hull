range = 1:360;
b = 0;
R1 = 1;
R2 = 0;
k1 = 10;
u1 = deg2rad(180);

y_fit = b+R1.*exp(k1.*(cos(deg2rad(range)-u1)-1))+R2.*exp(k1.*(cos(deg2rad(range)-u1-pi)-1));
figure; 
plot(range,y_fit)

component = circshift(y_fit,45)+ circshift(y_fit,-45);
hold on
plot(range,component)

pre_pattern = y_fit.*0.7 + component.*0.3;
plot(range,pre_pattern./max(pre_pattern(:)))
% k1 = 3;
% pre_pattern = b+R1.*exp(k1.*(cos(deg2rad(range)-u1)-1))+R2.*exp(k1.*(cos(deg2rad(range)-u1-pi)-1));
% plot(range,pre_pattern)
axis square
ylim([0 1])
xlim([0 360])
print('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Grants\CrossOri_R01\inhib_schem.pdf','-dpdf')