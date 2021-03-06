%plot % of failed trials with a lick vs the difference: peak df/f on
%early trials with a lick - early trials without a lick.
clear;

%define some variables and pathnames
WF_plotting_lists_of_days;
early_TC_dir = 'Z:\Analysis\WF Lever Analysis\licking_investigation\early_trials_lick_v_no_lick\no_licks_-2_to_2\';
corr_TC_dir = 'Z:\Analysis\WF Lever Analysis\licking_investigation\correct_trials_lick_v_no_lick\no_licks_-2_to_2\';
struct_dir = 'Z:\Analysis\WF Lever Analysis\StructuresForPlotting\';
output_dir = 'Z:\Analysis\WF Lever Analysis\licking_investigation\percent licks on earlies scatterplots\';
colors = {'k','k', 'b','b', 'g','g', 'c','c', 'm','m', 'r','r', 'y','y', 'k', 'b', 'r'}; 
peak_window = [4:7];

%generate some empty cell arrays to be filled in
peak_diff_early = {};  peak_diff_early{1} = [];
peak_diff_corr = {};   peak_diff_corr{1} = [];
percent_early_with_lick = {};
percent_early_with_lick{1} = [];
percent_corr_with_lick = {}; percent_corr_with_lick{1} = [];
percent_correct= {}; percent_correct{1} = [];
corr_early_diff = {}; corr_early_diff{1}=[];
corr_early_ratio = {}; corr_early_ratio{1}=[];
peak_ratio_lick_corr_early = {}; peak_ratio_lick_corr_early{1} = [];

%load tc data and determine differences
for session_num = [[1,2], [5:length(days)]]
    %load indeces and struct
    load([corr_TC_dir, days{session_num}, '_trial_indeces']);
    load([early_TC_dir, days{session_num}, '_trial_indeces']);
    load([struct_dir, days{session_num}]);
    
    %get percent of early/corr trials with lick
    percent_corr_with_lick{session_num} = length(lick_corr_trials)/(length(no_lick_corr_trials) + length(lick_corr_trials));
    percent_early_with_lick{session_num} = length(lick_fail_trials)/(length(no_lick_fail_trials) + length(lick_fail_trials));
    
    %get percent correct
    percent_correct{session_num} = curr_struct.percent_correct;
    
    %get the diff of peak df/f for correct - early and corr/early
    corr_early_diff{session_num} = curr_struct.corr_magnitude_mean - curr_struct.early_magnitude_mean;
    corr_early_ratio{session_num} = curr_struct.corr_magnitude_mean./curr_struct.early_magnitude_mean;
    
    %if there are early TC data
    if exist([early_TC_dir, days{session_num}, '.mat']) ~= 0
        %load data
        load([early_TC_dir, days{session_num}, '.mat']);
        no_lick_TCs = early_TCs.no_lick_TCs;
        lick_TCs = early_TCs.lick_TCs;
        
        %make sure there are enough trials
        if size(lick_TCs,1) > 3 & size(no_lick_TCs,1) > 3
            %get rid of non_LS ROIs
            no_lick_TCs = no_lick_TCs(:,LS_ROIs{session_num},:);
            lick_TCs = lick_TCs(:,LS_ROIs{session_num},:);
            
            %calculate the the difference in peak df/f for lick vs no lick earlies
            peak_lick_fail = max(lick_TCs(:,:,peak_window),[],3);
            peak_no_lick = max(no_lick_TCs(:,:,peak_window),[],3);
            peak_lick_mean = mean(peak_lick_fail);
            peak_no_lick_mean = mean(peak_no_lick);
            if length(LS_ROIs{session_num}) > 1
                peak_lick_mean = mean(peak_lick_mean);
                peak_no_lick_mean = mean(peak_no_lick_mean);
            end
            peak_diff_early{session_num} = peak_lick_mean-peak_no_lick_mean;
        end
    end
    
    %if there are corr TC data  - get the diff of peak df/f for correct with licks - corrects without licks
    if exist([corr_TC_dir, days{session_num}, '.mat']) ~= 0
        %load data
        load([corr_TC_dir, days{session_num}, '.mat']);
        no_lick_TCs = corr_TCs.no_lick_TCs;
        lick_TCs = corr_TCs.lick_TCs;
        
        %make sure there are enough trials
        if size(lick_TCs,1) > 3 & size(no_lick_TCs,1) > 3
            %get rid of non_LS ROIs
            no_lick_TCs = no_lick_TCs(:,LS_ROIs{session_num},:);
            lick_TCs = lick_TCs(:,LS_ROIs{session_num},:);
            
            %calculate the the difference in peak df/f for lick vs no lick earlies
            peak_lick_corr = max(lick_TCs(:,:,peak_window),[],3);
            peak_no_lick = max(no_lick_TCs(:,:,peak_window),[],3);
            peak_lick_mean = mean(peak_lick_corr);
            peak_no_lick_mean = mean(peak_no_lick);
            if length(LS_ROIs{session_num}) > 1
                peak_lick_mean = mean(peak_lick_mean);
                peak_no_lick_mean = mean(peak_no_lick_mean);
            end
            peak_diff_corr{session_num} = peak_lick_mean-peak_no_lick_mean;
        end
    end
    
    %calculate the ratio of corr/early peak df/f for lick trials
    if exist([corr_TC_dir, days{session_num}, '.mat']) ~= 0 &  exist([early_TC_dir, days{session_num}, '.mat']) ~= 0
        if size(corr_TCs.lick_TCs,1) > 3 & size(early_TCs.lick_TCs,1) > 3
            lick_TCs = corr_TCs.lick_TCs;
            lick_TCs = lick_TCs(:,LS_ROIs{session_num},:);
            peak_lick_corr = max(lick_TCs(:,:,peak_window),[],3);
            lick_TCs = early_TCs.lick_TCs;
            lick_TCs = lick_TCs(:,LS_ROIs{session_num},:);
            peak_lick_fail = max(lick_TCs(:,:,peak_window),[],3);
            peak_ratio_lick_corr_early{session_num} = mean(peak_lick_corr)./mean(peak_lick_fail);
        end
    end
end

%% scatterplots
%plot diff of peak df/f for failed trials with licks - without licks  vs   percent of early trials with licks
peak_diff_early_mat = [];
percent_early_with_lick_mat = [];
figure; diag_line = [-1:10];
for session_num = find(~cellfun(@isempty, peak_diff_early));
    plot(peak_diff_early{session_num}, percent_early_with_lick{session_num}*100, strcat('o', colors{session_num}) ); hold on;
    peak_diff_early_mat = [peak_diff_early_mat, peak_diff_early{session_num}];
    percent_early_with_lick_mat = [percent_early_with_lick_mat, percent_early_with_lick{session_num}];
end
[Rho, P_value] = corr(peak_diff_early_mat', percent_early_with_lick_mat');
%plot(diag_line, diag_line); 
title(['Rho=', num2str(Rho), ' P-value=', num2str(P_value)]);
xlabel('peak df/f for earlies with licks - earlies without licks');
ylabel('% of early trials with a lick');
xlim([-0.02 0.15]);
ylim([0 100]);
%savefig([output_dir, 'percent_earlies_with_lick_vs_diff_in_peak_f_earlies_lick_minus_no_lick']);

%plot percent correct vs   percent of early trials with licks
percent_correct_mat = [];
percent_early_with_lick_mat = [];
figure;
for session_num = find(~cellfun(@isempty, percent_correct));
    plot(percent_correct{session_num}, percent_early_with_lick{session_num}*100, strcat('o', colors{session_num}) ); hold on;
    percent_correct_mat = [percent_correct_mat, percent_correct{session_num}];
    percent_early_with_lick_mat = [percent_early_with_lick_mat, percent_early_with_lick{session_num}];
end
[Rho, P_value] = corr(percent_correct_mat', percent_early_with_lick_mat');
ylim([0 100]);
xlim([0 100]);
xlabel('percent correct');
ylabel('% of early trials with a lick');
title(['Rho=', num2str(Rho), ' P-value=', num2str(P_value)]);
%savefig([output_dir, 'percent_earlies_with_lick_vs_percent_correct']);
    
%plot diff of peak df/f for correct with licks - corrects without licks vs   percent of early trials with licks
peak_diff_corr_mat = [];
percent_early_with_lick_mat = [];
figure; diag_line = [-1:10];
for session_num = find(~cellfun(@isempty, peak_diff_corr));
    plot(peak_diff_corr{session_num}, percent_early_with_lick{session_num}*100, strcat('o', colors{session_num}) ); hold on;
    peak_diff_corr_mat = [peak_diff_corr_mat, peak_diff_corr{session_num}];
    percent_early_with_lick_mat = [percent_early_with_lick_mat, percent_early_with_lick{session_num}];
end
[Rho, P_value] = corr(peak_diff_corr_mat', percent_early_with_lick_mat');
%plot(diag_line, diag_line); 
title(['Rho=', num2str(Rho), ' P-value=', num2str(P_value)]);
xlabel('peak df/f for corrects with licks - corrects without licks');
ylabel('% of early trials with a lick');
xlim([-0.02 0.15]);
ylim([0 100]);
%savefig([output_dir, 'percent_earlies_with_lick_vs_diff_in_peak_f_corrects_lick_minus_no_lick']);

%plot diff of peak df/f for correct - early  vs   percent of early trials with licks
%switched to using a ratio corr/early instead. Only using trials with licks
peak_diff_mat = [];
percent_early_with_lick_mat = [];
figure; diag_line = [-1:10];
for session_num = find(~cellfun(@isempty, peak_ratio_lick_corr_early));
    plot(peak_ratio_lick_corr_early{session_num}, percent_early_with_lick{session_num}*100, strcat('o', colors{session_num}) ); hold on;
    peak_diff_mat = [peak_diff_mat, mean(peak_ratio_lick_corr_early{session_num})];
    percent_early_with_lick_mat = [percent_early_with_lick_mat, percent_early_with_lick{session_num}];
end
[Rho, P_value] = corr(peak_diff_mat', percent_early_with_lick_mat');
%plot(diag_line, diag_line); 
title(['Rho=', num2str(Rho), ' P-value=', num2str(P_value)]);
xlabel('peak df/f for corrects/earlies');
ylabel('% of early trials with a lick');
xlim([0 5]);
ylim([0 100]);
savefig([output_dir, 'percent_earlies_with_lick_vs_corr_early_ratio_licks_only']);








