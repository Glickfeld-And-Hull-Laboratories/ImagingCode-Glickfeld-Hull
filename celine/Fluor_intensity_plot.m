plot(data(:,1),data(:,2),'-o')
hold on
plot(data(:,1),data(:,3),'-o')
hold on
plot(data(:,1),(data(:,2)-data(:,3)),'-o')
xlabel('Hours since infusion')
ylabel('Farred fluorescence intensity')