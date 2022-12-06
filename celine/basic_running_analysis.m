%% looking at wheel speed

[wheel_speed] = wheelSpeedCalc(input,32,'orange');
%arguments for wheelSpeedCals:  the name of your mWorks input
%structure (by default, it's just called "input"; the number of clicks in a
%rotation of the encoderm which is 32; which wheel you used, because thye
%are different sizes so a full rotation corresponds to a different distance
%of running depending on which wheel. %orange is teh hard plastic organe or
%blue ones, "purple" is the spongey purple one, "red" is the tilted red
%plastic one.
figure; plot(wheel_speed)
wheel_tc = zeros(nOns+nOffs, nTrials);
for iTrial = 1:nTrials
    wheel_tc(:,iTrial) = wheel_speed(1+((iTrial-1).*(nOns+nOffs)):iTrial.*(nOns+nOffs));
end
wheel_trial_avg = mean(wheel_tc(nOffs:nOns+nOffs,:),1);
figure; movegui('center')
plot(wheel_trial_avg)

RIx = wheel_trial_avg>0.5;


