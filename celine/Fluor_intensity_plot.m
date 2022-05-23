data_2055 = table2array(i2055intensity);
data_2056 = table2array(i2056intensity);

figure;
semilogx(data_2055(:,1),data_2055(:,4),'-o')
hold on
semilogx(data_2056(:,1),data_2056(:,4),'-o')
legend
ylim([0 1])
set(gca, 'TickDir', 'out')
xlabel('Hours since infusion')
ylabel('Normalized fluorescence intensity')
hold off
box off

%% 
data_2052 = table2array(i2052intensity);
data_2053 = table2array(i2053intensity);

data_2052(:,7)=data_2052(:,2)./max(data_2052(:,2));
data_2052(:,8)=data_2052(:,3)./max(data_2052(:,2));
data_2052(:,9)=data_2052(:,2)-data_2052(:,3);
data_2052(:,10)=data_2052(:,9)./max(data_2052(:,2));

figure;
plot(data_2052(:,1 ),data_2052(:,7),'-o')
hold on
plot(data_2052(:,1 ),data_2052(:,8),'-o')
hold on
plot(data_2052(:,1 ),data_2052(:,10),'-o')




data_2053(:,7)=data_2053(:,2)./max(data_2053(:,2));
data_2053(:,8)=data_2053(:,3)./max(data_2053(:,2));
data_2053(:,9)=data_2053(:,2)-data_2053(:,3);
data_2053(:,10)=data_2053(:,9)./max(data_2053(:,2));

figure;
plot(data_2053(:,1 ),data_2053(:,7),'-or')
hold on
plot(data_2053(:,1 ),data_2053(:,8),'-ok')
hold on
plot(data_2053(:,1 ),data_2053(:,10),'-o','Color',[0.5 0 0.8])
title({'i2053 fluorescence', ' normalized to max HT+ fluorescence'})
legend('HT+','control','delta fluor')
set(gca, 'TickDir', 'out')
xlabel('Hours since infusion')
ylabel('Normalized fluorescence intensity')
hold off
box off

figure;
plot(data_2052(:,1 ),(data_2052(:,10) ./max(data_2052(:,10))),'-o','Color',[0.25 0 1])
hold on
plot(data_2053(:,1 ),(data_2053(:,10) ./max(data_2053(:,10))),'-o','Color',[0.5 0 0.8])
legend('i2052','i2053')
title({'Delta fluorescence', ' normalized to max delta fluorescence'})
set(gca, 'TickDir', 'out')
xlabel('Hours since infusion')
ylabel('Normalized fluorescence intensity')
hold off
box off