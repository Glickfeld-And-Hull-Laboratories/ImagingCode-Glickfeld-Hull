figure
semilogx(data(:,1),data(:,3),'-o')
hold on
% plot(data(:,1),data(:,5),'-o')
% hold on
% plot(data(:,1),data(:,6),'-o')
xlabel('Hours since infusion, log')
ylabel('Farred fluorescence intensity')
title('i2055 normalized fluorescence')
hold off

