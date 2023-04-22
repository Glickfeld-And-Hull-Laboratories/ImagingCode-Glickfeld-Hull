clear all
clear all global
close all

% Import eye summaries
data_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\ACh\Aging';
CD = [data_pn '\Gloria\eyeTrials'];
cd(CD);
eyeTrials = dir('*.mat'); 
nFiles = length(eyeTrials);

for i = 1:nFiles 
  data{i} = load(eyeTrials(i).name); 
end

%% Distribution of pupil radius (mean and STD)
%make a for loop which loops through every mouse, finds the mean pupil radius and SD and plots it on a graph

mouse = cell(1,nFiles);
age = zeros(1,nFiles);
DOB = cell(1,nFiles);
record_date = cell(1,nFiles);
rad = cell(1,nFiles);
mean_rad = zeros(1,nFiles);
std_rad = zeros(1,nFiles);

figure;
for i = 1:nFiles
    mouse{i} = data{i}.mouse;
    DOB{i} = data{i}.DOB;
    record_date{i} = data{i}.record_date;
    timepoints = datetime({DOB{i}; record_date{i}}, 'InputFormat','yyyy-MM-dd');
    time_diff = caldiff(timepoints, 'weeks'); % calculate age at time of recording in months
    [age(i)] = split(time_diff, 'weeks');

    rad{i} = data{i}.rad.stim;
    mean_rad(i) = mean(rad{i}, 'omitnan'); %mean pupil radius
    std_rad(i) = std(rad{i}, 'omitnan'); %STD of pupil radius

    %plot mean pupil radius across all trials vs. age
    errorbar(age(i), mean_rad(i), std_rad(i), 'o') % plot age as x, mean as y, std as error
    hold on
    axis([0 115 0 15])
    text(age(i)+2.0, mean_rad(i), mouse{i}) %labels each point with mouse(i) ID
    xlabel('Age (wks)')
    ylabel("Mean Pupil Radius + STD")
    title(['Pupil Radius vs. Age'], "FontSize", 12)

end 
    
%% Dependence of Running on Age
%plot fraction of run trials vs. age

run_ind = cell(1,nFiles);
notrun_ind = cell(1,nFiles);
fraction_run = zeros(1,nFiles);

figure;
for i = 1:nFiles
    run_ind{i} = data{i}.run_ind;
    notrun_ind{i} = data{i}.notrun_ind;
    nRun = length(run_ind{i});
    nNotRun = length(notrun_ind{i});
    nTotal = nRun + nNotRun;
    fraction_run(i) = nRun/nTotal; %proportion of run trials over total

    scatter(age(i), fraction_run(i)) % plot age as x, fraction running as y
    hold on
    axis([0 115 0 0.5])
    text(age(i)+2.0, fraction_run(i), mouse{i}) %labels each point with mouse ID
    xlabel('Age (wks)')
    ylabel("Fraction of Run Trials")
    title(['Proportion of Run Trials vs. Age'], "FontSize", 12)
end 

%% Average Run Speed vs. Age
%plot average speed of running trials vs. age

wheel_trial_avg = cell(1,nFiles);
run_trials = cell(1,nFiles);
mean_run_speed = zeros(1,nFiles);
std_run_speed = zeros(1,nFiles);

figure;
for i = 1:nFiles
    wheel_trial_avg{i} = data{i}.wheel_trial_avg;
    run_trials{i} = wheel_trial_avg{i}(run_ind{i}); %gives wheel speeds of running trials only
    mean_run_speed(i) = mean(run_trials{i}); %mean running trial speed
    std_run_speed(i) = std(run_trials{i}); %STD running trial speed

    errorbar(age(i), mean_run_speed(i), std_run_speed(i), 'o') % plot age as x, mean speed as y, std speed as error
    hold on
    axis([0 115 0 25])
    text(age(i)+2.0, mean_run_speed(i), mouse{i}) %labels each point with mouse ID
    xlabel('Age (wks)')
    ylabel("Mean Running Speed + STD")
    title(['Running Trial Speed vs. Age'], "FontSize", 12)

end

%% Relationship of pupil radius and running speed across age

pearson_coef = zeros(1, nFiles); %array of correlation coefficient between pupil and running for each mouse
clean_wheel_trial_avg = cell(1,nFiles);

figure;
for i = 1:nFiles
    clean_wheel_trial_avg{i} = data{i}.wheel_trial_avg;
    clean_wheel_trial_avg{i}(data{i}.nan_ind) = []; %removes the NaN trials from rad.stim for wheel data
    clean_rad = rmmissing(rad{i});
    R = corrcoef(clean_rad, clean_wheel_trial_avg{i}); %returns coef matrix between pupil diameter and wheel speed across trials
    pearson_coef(i) = triu2vec(R); %extracts the coef value from the matrix as a single vector

    %plot correlation coeffients vs. age
    scatter(age(i), pearson_coef(i)) % plot age as x, coefs as y
    hold on
    axis([0 115 0 1.0])
    text(age(i)+2.0, pearson_coef(i), mouse{i}) %labels each point with mouse ID
    xlabel('Age (wks)')
    ylabel("Correlation Coefficient")
    title(['Relationship of Pupil Radius and Running Speed vs. Age'], "FontSize", 12)

end
%% Variability in pupil size vs. variability in running speed

figure;
for i = 1:nFiles
    scatter(std_rad(i), std_run_speed(i)) % plot STD pupil as x, STD running speed as y
    hold on
    axis([0 2 0 7.0])
    text(std_rad(i)+0.04, std_run_speed(i), mouse{i}) %labels each point with mouse ID
    xlabel('STD of Pupil Radius')
    ylabel("STD of Running Trial Speed")
    title(['Variability of Pupil Radius and Running Speed'], "FontSize", 12)
end 

%% Mean Pupil Radius Stationary vs. Running
% plot mean pupil radius for running trials and for stationary trials 

mean_rad_run = zeros(1,nFiles);
std_rad_run = zeros(1,nFiles);

mean_rad_notrun = zeros(1,nFiles);
std_rad_notrun = zeros(1,nFiles);

figure;
for i = 1:nFiles
    mean_rad_run(i) = mean(rad{i}(run_ind{i}), 'omitnan'); %mean pupil diameter across run trials
    std_rad_run(i) = std(rad{i}(run_ind{i}), 'omitnan'); %STD pupil diameter across run trials

    mean_rad_notrun(i) = mean(rad{i}(notrun_ind{i}), 'omitnan'); % mean pupil diameter across stationary trials
    std_rad_notrun(i) = std(rad{i}(notrun_ind{i}), 'omitnan'); %STD pupil diameter across stationary trials

    subplot(1,2,1) %plot mean pupil radius during stationary trials
    errorbar(age(i), mean_rad_notrun(i), std_rad_notrun(i), 'o') % plot age as x, mean as y, std as error
    hold on
    axis([0 115 0 18])
    text(age(i)+3.0, mean_rad_notrun(i), mouse{i}) %labels each point with mouse ID
    xlabel('Age (wks)')
    ylabel("Mean Pupil Radius + STD")
    title(['Stationary Trials'], "FontSize", 12)

    subplot(1,2,2) %plot mean pupil radius during running trials
    errorbar(age(i), mean_rad_run(i), std_rad_run(i), 'o') % plot age as x, mean as y, std as error
    hold on
    axis([0 115 0 18])
    text(age(i)+3.0, mean_rad_run(i), mouse{i}) %labels each point with mouse ID
    xlabel('Age (wks)')
    ylabel("Mean Pupil Radius + STD")
    title(['Running Trials'], "FontSize", 12)

    sgtitle('Mean Pupil Radius vs. Age')

end

%% Sort trials into small vs. large pupil
% small pupils will be bottom 30th percentile, large pupils will be top 30th percentile (i.e. 70th percentile)

small_trials = cell(1,nFiles);
large_trials = cell(1,nFiles);

mean_rad_small = zeros(1,nFiles);
std_rad_small = zeros(1,nFiles);

mean_rad_large = zeros(1,nFiles);
std_rad_large = zeros(1,nFiles);

figure;
for i = 1:nFiles
    bottom_30 = prctile(rad{i}, 30); %bottom 30th percentile threshold
    top_30 = prctile(rad{i}, 70); %top 30th percentile threshold
    
    small_trials{i} = find(rad{i}<=bottom_30); % indices of trials that are bottom 30th percentile
    mean_rad_small(i) = mean(rad{i}(small_trials{i}), 'omitnan'); %mean pupil diameter for small pupil trials
    std_rad_small(i) = std(rad{i}(small_trials{i}), 'omitnan'); %STD of pupil diameter for small pupil trials

    large_trials{i} = find(rad{i}>=top_30); % indices of trials that are top 30th percentile
    mean_rad_large(i) = mean(rad{i}(large_trials{i}), 'omitnan'); %mean pupil diameter for large pupil trials
    std_rad_large(i) = std(rad{i}(large_trials{i}), 'omitnan'); %STD of pupil diameter for large pupil trials

    subplot(1,2,1) %plot mean pupil radius in small pupil trials
    errorbar(age(i), mean_rad_small(i), std_rad_small(i), 'o') % plot age as x, mean as y, std as error
    hold on
    axis([0 115 0 18])
    text(age(i)+3.0, mean_rad_small(i), mouse{i}) %labels each point with mouse ID
    xlabel('Age')
    ylabel("Mean Pupil Radius + STD")
    title(['Small Pupil Trials'], "FontSize", 12)

    subplot(1,2,2) %plot mean pupil radius in large pupil trials
    errorbar(age(i), mean_rad_large(i), std_rad_large(i), 'o') % plot age as x, mean as y, std as error
    hold on
    axis([0 115 0 18])
    text(age(i)+3.0, mean_rad_large(i), mouse{i}) %labels each point with mouse ID
    xlabel('Age')
    ylabel("Mean Pupil Radius + STD")
    title(['Large Pupil Trials'], "FontSize", 12)

    sgtitle('Pupil Radius vs. Age')

end
%% 
% plot mean running speed for small pupil vs large pupil trials across age

mean_run_small = zeros(1,nFiles);
std_run_small = zeros(1,nFiles);

mean_run_large = zeros(1,nFiles);
std_run_large = zeros(1,nFiles);

figure;
for i = 1:nFiles
    mean_run_small(i) = mean(wheel_trial_avg{i}(small_trials{i})); %mean running speed for small pupil trials
    std_run_small(i) = std(wheel_trial_avg{i}(small_trials{i})); %SD of running speed for small pupil trials

    mean_run_large(i) = mean(wheel_trial_avg{i}(large_trials{i})); %mean running speed for large pupil trials
    std_run_large(i) = std(wheel_trial_avg{i}(large_trials{i})); %SD of running speed for large pupil trials

    subplot(1,2,1) %plot mean running speed in small pupil trials
    errorbar(age(i), mean_run_small(i), std_run_small(i), 'o') % plot age as x, mean as y, std as error
    hold on
    axis([0 115 -2 25])
    text(age(i)+2.0, mean_run_small(i), mouse{i}) %labels each point with mouse ID
    xlabel('Age')
    ylabel("Mean Running Speed + STD")
    title(['Small Pupil Trials'], "FontSize", 12)

    subplot(1,2,2) %plot mean running speed in large pupil trials
    errorbar(age(i), mean_run_large(i), std_run_large(i), 'o') % plot age as x, mean as y, std as error
    hold on
    axis([0 115 -2 25])
    text(age(i)+2.0, mean_run_large(i), mouse{i}) %labels each point with mouse ID
    xlabel('Age')
    ylabel("Mean Running Speed + STD")
    title(['Large Pupil Trials'], "FontSize", 12)

    sgtitle('Running Speed vs. Age')

end

%% Export data 
fnout = [data_pn, '\Gloria\'];

varNames = {'Mouse ID','DOB','Record Date', 'Age (wks)'};
mouse_overview = table(mouse', DOB', record_date', age', 'VariableNames', varNames)
writetable(mouse_overview,'mouse_overview.xlsx')

for i = 1:nFiles
    mouse = mouse{i};
    age = age(i);
    datemouserun = data{i}.datemouserun;
    small_trials = small_trials{i};
    large_trials = large_trials{i};
    rad = rad{i};
    save(fullfile(fnout, "pupilSort", [data{i}.datemouserun '_pupilSort.mat']), 'mouse', 'age', 'datemouserun', 'small_trials', 'large_trials', 'rad') % saves a separate file of the pupil sorting
end
