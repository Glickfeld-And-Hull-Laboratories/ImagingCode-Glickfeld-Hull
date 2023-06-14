clear all
clear all global
close all

% Import runOnsets
data_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\ACh\Aging';
CD = [data_pn '\Gloria\runOnset'];
cd(CD);
runOnsets = dir('*.mat'); 
nFiles = length(runOnsets);

for i = 1:nFiles 
    all_data{i} = load(runOnsets(i).name);
end

%% Calculate fractional pupil change and latency from running onset

nFrames = 60; %total frames for running onset window
befoRun = 30; %frames before run onset
aftRun = 30; %frames after run onset

%create vector for time axis
t = [1:nFrames];
t = (t-30)./15; %re-centers start of running (frame 30) as 0 and converts frames to seconds

avg_runTrigger_pupil = zeros(nFrames, nFiles);
avg_runTrigger_speed = zeros(nFrames, nFiles);

[row col] = subplotn(nFiles); 

for i = 1:nFiles
    avg_runTrigger_pupil(:,i) = all_data{i}.avg_runTrigger_pupil;
    avg_runTrigger_speed(:,i) = all_data{i}.avg_runTrigger_speed;

end

for i = 1:nFiles
    [r, lags] = xcorr(avg_runTrigger_speed(:,i), avg_runTrigger_pupil(:,i));
    [max_r(i), t_max(i)] = max(r); 
    subplot(row,col,i)
    stem(lags,r)
    hold on
    title([all_data{i}.mouse ' (' num2str(all_data{i}.age) 'wks' ')' ' - ' num2str(all_data{i}.nWindows) ' windows'], "FontSize", 10)
end

%% Excluding data based on minimum nWindows 
% exclude mouse if nWindows is <10

data = {};

for i = 1:nFiles
    nWindows = all_data{i}.nWindows;
    if nWindows<10
    continue 
    else
        data = cat(2, data, all_data{i});
    end
end 
%% 

nFiles = length(data);
avg_runTrigger_pupil = zeros(nFrames, nFiles);
avg_runTrigger_speed = zeros(nFrames, nFiles);
age = zeros(1, nFiles);
mouse = cell(1,nFiles);

for i = 1:nFiles
    runTrigger_speed = data{i}.runTrigger_speed;
    runTrigger_pupil_dfof = data{i}.runTrigger_pupil_dfof;
    avg_runTrigger_pupil(:,i) = data{i}.avg_runTrigger_pupil;
    avg_runTrigger_speed(:,i) = data{i}.avg_runTrigger_speed;
    age(i) = data{i}.age;
    mouse{i} = data{i}.mouse;
    nWindows = data{i}.nWindows;
end 
%% 

avg_baseline_pupil = mean(avg_runTrigger_pupil(1:befoRun/2,:));
[avg_max_pupil, max_ind] = max(avg_runTrigger_pupil(aftRun+1:nFrames,:));
max_ind = max_ind+befoRun;

frac_change = avg_max_pupil-avg_baseline_pupil; %finding the fractional pupil change for avg TC
avg_run_speed = mean(avg_runTrigger_speed(aftRun+1:nFrames,:)); %average overall run speed for each mouse


for i = 1:nFiles
   runTrigger_speed = data{i}.runTrigger_speed;
   runTrigger_pupil_dfof = data{i}.runTrigger_pupil_dfof;
   baseline_pupil = mean(runTrigger_pupil_dfof(1:befoRun/2,:),1);
   resp_pupil = mean(runTrigger_pupil_dfof(aftRun+1:nFrames,:),1);
    if frac_change(i)<0
        [h(i),p(i)] = ttest(baseline_pupil,resp_pupil,'tail','right');
    else 
        [h(i),p(i)] = ttest(baseline_pupil,resp_pupil,'tail','left');
    end
end 
%% 

latency = [];

for i = 1:nFiles
    if h(i) == 1
        latency(i) = t(max_ind(i)); %find the latency of hitting max pupil size from run onset
    end
end 

%% Plot all average TCs from each mice in one figure

[row col] = subplotn(nFiles); 
left_color = [0 0 0];
right_color = [0 0 0];

set(figure,'defaultAxesColorOrder',[left_color; right_color]);
for i = 1:nFiles
    subplot(row,col,i)
    shadedErrorBar(t,data{i}.avg_runTrigger_speed(:),data{i}.std_err_speed(:),'lineProps','-r'); % plot average time course with error bars of speed across all running windows
    axis([t(1) t(end) -2 20])
    vline(0,'black'); % indicate running onset
    hold on
    ylabel("Speed (pulses/s)")
    xlabel('Time (s)')
    title([data{i}.mouse ' (' num2str(data{i}.age) 'm' ')' ' - ' num2str(data{i}.nWindows) ' windows'], "FontSize", 10)
 
    yyaxis right
    axis([t(1) t(end) -0.2 0.8])
    vline(latency(i), '--') %indicate max pupil
    ylabel("Pupil Fractional Change", 'FontSize', 9)
    shadedErrorBar(t,data{i}.avg_runTrigger_pupil(:),data{i}.std_err_pupil(:), 'lineProps','-b'); % plot average time course with error bars of speed across all running windows
    
    sgtitle('Average time course for running onset');

    if i==1
        legend({'Running Speed', 'Pupil Size'}, 'Location', 'northwest', 'FontSize', 8)
    end
end 

%% Plot pupil fractional change and pupil latency vs. age

figure;
for i = 1:nFiles
    subplot(1,2,1)
    scatter(age(i), frac_change(i)) % plot age as x, fraction change as y
    hold on
    axis([0 26 -0.1 0.3])
    text(age(i)+0.75, frac_change(i), mouse(i)) %labels each point with mouse ID
    xlabel('Age')
    ylabel("Max Pupil Fractional Change")
    title(['Pupil Change vs. Age'], "FontSize", 12)
    
    subplot(1,2,2)
    if latency(i) == 0
        continue
    else
    scatter(age(i), latency(i)) % plot age as x, fraction running as y
    hold on
    axis([0 26 0 2.0])
    text(age(i)+0.75, latency(i), mouse(i)) %labels each point with mouse ID
    xlabel('Age')
    ylabel("Time (s)")
    title(['Max Pupil Latency vs. Age'], "FontSize", 12)
    end
end 

%% Plot Pupil Fractional Change vs. Avg Running Speed, and the ratio of these against age

figure;
for i = 1:nFiles
    subplot(1,2,1)
    scatter(avg_run_speed(i), frac_change(i))
    hold on
    axis([0 12 -0.1 0.3])
    xlabel('Average Running Speed')
    ylabel('Pupil Fractional Change')
    text(avg_run_speed(i)+0.25,frac_change(i), mouse(i)) %labels each point with mouse ID

    subplot(1,2,2)
    scatter(age(i), (frac_change(i)/avg_run_speed(i)))
    hold on
    axis([0 26 0 0.1])
    xlabel('Age (months)')
    ylabel('Pupil Fractional Change/Running Speed')
    text(age(i)+0.25,(frac_change(i)/avg_run_speed(i)), mouse(i)) %labels each point with mouse ID
end 

%% Plot pupil fractional change vs. latency time

figure;
for i = 1:nFiles
    if latency(i) == 0
        continue
    else 
    scatter(frac_change(i), latency(i))
    hold on
    axis([-0.05 0.3 0 2])
    xlabel('Pupil Fractional Change')
    ylabel('Latency Time (s)')
    text(frac_change(i)+0.005, latency(i), mouse(i)) %labels each point with mouse ID
    end 
end 

%% Cross correlogram between pupil TC and running TC

for i = 1:nFiles
    [r, lags] = xcorr(data{i}.avg_runTrigger_speed(:), data{i}.avg_runTrigger_pupil(:));
    subplot(row,col,i)
    stem(lags,r)
    hold on
    title([data{i}.mouse ' (' num2str(data{i}.age) 'm' ')' ' - ' num2str(data{i}.nWindows) ' windows'], "FontSize", 10)
end