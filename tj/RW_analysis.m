clear global
clear all
close all
%%    %Enter variables
date = '220318';
time = strvcat('0951');
folder = 'running';
alignToRef = 1;
irun = 1;
ImgFolder = strvcat('001');
mouse = 'i1369';
nrun = size(ImgFolder,1);
frame_rate = 15.5;
ref_str = 'runs-001';
run_str = catRunName(ImgFolder, nrun);
tj_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\2P_Imaging';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\2P';
behav_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';

fName = fullfile(behav_fn, ['data-i' '''' mouse '''' '-' date '-' time(irun,:) '.mat']);
load(fName);

out_path = fullfile(fnout, [date '_' mouse], folder);

mkdir(out_path);

nOns= input.nScansOn;
nOffs= input.nScansOff;
nTrials = input.stopAfterNTrials;


%%    %looking at wheel speed during on trials (plots on wheel trial average)
    
[wheel_speed] = wheelSpeedCalc(input,32,'purple');
figure(1); 
plot(wheel_speed);
wheel_tc = zeros(nOns+nOffs, nTrials);
for iTrial = 1:nTrials;
    wheel_tc(:,iTrial) = wheel_speed(1+((iTrial-1).*(nOns+nOffs)):iTrial.*(nOns+nOffs));
end
wheel_trial_avg = mean(wheel_tc(nOffs:nOns+nOffs,:),1);
tmp = (wheel_tc(nOffs:nOns+nOffs,:));

figure(2); movegui('center');
plot(wheel_trial_avg);
cd(out_path);
[t,s] = title([mouse ' running average across "on" period'], ['date: ' date '  ' 'time: ' time]);
xlabel('trials');
ylabel('cm/s');
print(fullfile(out_path, [mouse '_' date '_' time '_running.pdf']),'-dpdf','-bestfit')

RIx = wheel_trial_avg>5.5;
% ind1 = find(wheel_trial_avg<=0);
% ind2 = find(wheel_trial_avg>0 & wheel_trial_avg<0.2);
% ind3 = find(wheel_trial_avg>=0.2);


%%    %graph of wheel avg on vs off

sep_on = wheel_tc((nOffs+1):(nOns+nOffs), :);
sep_off = wheel_tc((1:nOffs), :);

wheel_trial_avg_on = mean(sep_on,1);
wheel_trial_avg_off = mean(sep_off,1);

wheel_avg_fat_on = mean(sep_on,2);
wheel_avg_fat_off = mean(sep_off,2);

total_avg_on = mean(wheel_trial_avg_on);
total_avg_off = mean(wheel_trial_avg_off);

%%    %for line plot avg on vs off
figure(3);
plot(wheel_trial_avg_on, 'color', [0.5 0.7 0.2]);
hold on;
plot(wheel_trial_avg_off, 'color', [0.6350 0.0780 0.1840]);
cd(out_path);
[t,s] = title([mouse ' running average "on" (green) and "off" (red) periods across trials'], ['date: ' date '  ' 'time: ' time]);
xlabel('trials');
ylabel('cm/s');


%%    %for on vs off unity line
figure(4);
h(1)=scatter(wheel_trial_avg_on,wheel_trial_avg_off) 
hold on;

hl = refline(1,0);
hl.Color='r';

[t,s] = title([mouse ' "off" vs "on" periods running avg'], ['date: ' date '  ' 'time: ' time]);
xlabel('on trials');
ylabel('off trials');

%%    %averaging frames across trials for off and on (frame across trials = fat)
figure(5);
x=[wheel_avg_fat_off;wheel_avg_fat_on];
plot(x);
hold on;
xline(60, 'r');
[t,s] = title([mouse ' running average of frame across trials for "off" (before line) and "on" (after line) periods'], ['date: ' date '  ' 'time: ' time]);
xlabel('frame');
ylabel('cm/s');

%%    %Average across all trials on vs off
figure(6)
bar(1,total_avg_on, 'g');
hold on;
bar(2, total_avg_off, 'r');
[t,s] = title([mouse ' running average across all trials "on"(green) and "off"(red) periods'], ['date: ' date '  ' 'time: ' time]);
xlabel('');
ylabel('cm/s');
