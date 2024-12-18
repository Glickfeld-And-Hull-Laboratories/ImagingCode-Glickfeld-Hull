%get_peak_of_mean_lick
%this script find the peak of the average lick trace for each animal and
%then averages them together. 
%modified to be able to plot late df/f vs late lick rate

%set directories. Img30 has corrupted licking data. Exclude. 
clear;
lick_dir_base = 'Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\';
TC_dir = 'Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\';
days = {'151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36', '160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46', '160722_img53', '160904_img55'}; %'150718_img27', '150719_img27', '150716_img28', '150717_img28', 
ROI = {[2], [2], [1:3], [1,3], [1:4], [1:5], [1:2], [1:2], [1:2], [1:2], [3:4], [2,3,5], [1], [1:2], [3:4], [1,2,3], [3:5]}; 
max_rates_succ = [];
max_rates_fail = [];
late_lick_succ = [];
late_lick_fail = [];
late_lick_succ_sem = [];
late_lick_fail_sem = [];
late_f_succ = [];
late_f_fail = []; 
late_f_succ_sem = [];
late_f_fail_sem = [];
num_trials_succ = 0;
num_trials_fail = 0;

%main forloop for extracting data from crash
for ii = [1:2, 5:length(days)];
    %load lick traces
    lick_dir = [lick_dir_base, days{ii}, '_bx_outputs'];
    load(lick_dir, 'licking_data');
    load([TC_dir, days{ii}, '_success']);
    load([TC_dir, days{ii}, '_fail']);
    
    %find mean and sem for lick traces
    mean_lick_succ = mean(licking_data.lick_trace_succ,1);
    mean_lick_fail = mean(licking_data.lick_trace_fail,1);
    assert(length(mean_lick_succ)==16);
    assert(length(mean_lick_fail)==16);
    
    %find mean and sem for f on correct trials
    mean_late_f_succ = squeeze(mean(mean(success_roi(:,ROI{ii},[end-2:end]),3),2)); %ROIs and frames in a given trial
    this_late_f_succ_sem = std(mean_late_f_succ)/sqrt(length(mean_late_f_succ));
    mean_late_f_succ = mean(mean_late_f_succ); %avg together all the frames, then the ROIs
    
    %find mean and sem for f on early trials
    mean_late_f_fail = squeeze(mean(mean(fail_roi(:,ROI{ii},[end-2:end]),3),2)); %ROIs and frames in a given trial
    this_late_f_fail_sem = std(mean_late_f_fail)/sqrt(length(mean_late_f_fail));
    mean_late_f_fail = mean(mean_late_f_fail); %avg together all the frames, then the ROIs
    
    %find peak of the mean lick trace
    max_rates_succ = [max_rates_succ, max(mean_lick_succ(7:16))];
    max_rates_fail = [max_rates_fail, max(mean_lick_fail(7:16))];
    
    %find lick rate for late in the TC and sem
    mean_lick_succ = mean(licking_data.lick_trace_succ(:,[end-2:end]),2)*10;
    mean_lick_fail = mean(licking_data.lick_trace_fail(:,[end-2:end]),2)*10;
    late_lick_succ = [late_lick_succ, mean(mean_lick_succ)];
    late_lick_fail = [late_lick_fail, mean(mean_lick_fail)];
    late_lick_succ_sem = [late_lick_succ_sem, std(mean_lick_succ)/sqrt(length(mean_lick_succ))];
    late_lick_fail_sem = [late_lick_fail_sem, std(mean_lick_fail)/sqrt(length(mean_lick_fail))];
    
    %store the mean/sem df/f for late in the TC
    late_f_succ = [late_f_succ, mean_late_f_succ];
    late_f_fail = [late_f_fail, mean_late_f_fail];
    late_f_succ_sem = [late_f_succ_sem, this_late_f_succ_sem];
    late_f_fail_sem = [late_f_fail_sem, this_late_f_fail_sem];
    
    %record total trial nums
    num_trials_succ = num_trials_succ + size(success_roi, 1);
    num_trials_fail = num_trials_fail + size(fail_roi, 1);
end

%in the WF_lever_plotting_TCs the lick rates were divided by 10 in order to fit on the same plot as the df/f. correcting for that here. 
max_rates_succ= max_rates_succ*10;
max_rates_fail= max_rates_fail*10;

%make scatterplot to check that the data make sense 
figure; 
scatter(max_rates_succ, max_rates_fail);
xlabel('correct trials');
ylabel('incorrect trials');
title('peak lick rate of average lick traces');
ylim([0 20]); xlim([0 20]);
[h,p] = ttest(max_rates_succ, max_rates_fail)

%report mean and SEM. 
disp(['mean lick rate for correct trials = ' num2str(mean(max_rates_succ)) 'Hz  ', 'SEM=', num2str(std(max_rates_succ)/sqrt(length(max_rates_succ)))]);
disp(['mean lick rate for incorrect trials = ' num2str(mean(max_rates_fail)) 'Hz  ', 'SEM=', num2str(std(max_rates_fail)/sqrt(length(max_rates_fail)))]);

%% scatterplot of late df/f vs late lick rate
figure; scatter(late_f_succ, late_lick_succ, 'k'); hold on;
scatter(late_f_fail, late_lick_fail, 'r');
errorbar(late_f_succ, late_lick_succ, late_lick_succ_sem, 'Linestyle', 'none', 'Color', 'k');
errorbar(late_f_fail, late_lick_fail, late_lick_fail_sem, 'Linestyle', 'none', 'Color', 'r');
for ii = 1:length(late_f_succ)
    plot( [(late_f_succ(ii)-late_f_succ_sem(ii)),(late_f_succ(ii)+late_f_succ_sem(ii))] , [late_lick_succ(ii),late_lick_succ(ii)] , 'k'); %'Linestyle', 'none', 'Color', 'k'
    plot( [(late_f_fail(ii)-late_f_fail_sem(ii)),(late_f_fail(ii)+late_f_fail_sem(ii))] , [late_lick_fail(ii),late_lick_fail(ii)] , 'r');
end
xlabel('df/f');
ylabel('lick rate Hz');
title('df/f vs lick rate for frames 800 to 1000ms following lever release');

%linear regression of lick rate vs d/f for late trial frames
fitlm(late_f_succ, late_lick_succ)
fitlm(late_f_fail, late_lick_fail)

%% scatterplot of late df/f ratio vs late lick rate ratio
figure;
scatter((late_f_fail./late_f_succ), (late_lick_fail./late_lick_succ), 'k');
xlabel('df/f (early/correct)');
ylabel('lick rate Hz (early/correct)');
title('df/f ratio vs lick rate ratio 800 to 1000ms after lever release');




