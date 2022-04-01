%% looking at wheel speed

[wheel_speed] = wheelSpeedCalc(input,32,'purple');
figure; plot(wheel_speed)
wheel_tc = zeros(nOns+nOffs, nTrials);
for iTrial = 1:nTrials
    wheel_tc(:,iTrial) = wheel_speed(1+((iTrial-1).*(nOns+nOffs)):iTrial.*(nOns+nOffs));
end
wheel_trial_avg = mean(wheel_tc(nOffs:nOns+nOffs,:),1);
figure; movegui('center')
plot(wheel_trial_avg)

RIx = wheel_trial_avg>0.5;
ind1 = find(wheel_trial_avg<=0);
ind2 = find(wheel_trial_avg>0 & wheel_trial_avg<0.2);
ind3 = find(wheel_trial_avg>=0.2);

