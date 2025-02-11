function FirstTrain = FirstInTrain(LaserStim, preTrainTime);

k=2;
FirstTrain = LaserStim(1);
for i =2:length(LaserStim)
    if LaserStim(i)-LaserStim(i-1) > preTrainTime
        FirstTrain(k,1)=LaserStim(i);
        k = k + 1;
    end
end
