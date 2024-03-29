%plots session averages of timecourses for various conditions as well as
%early-earlies and late-earlies

clear
file_info
out_base = fullfile('Z:\home\ziye\2P_Analysis\2P_Analysis\');
mouseID = mouseID(1:30);
for id = 1:size(mouseID,2)
    for rID  = 1:2
        dest_sub  = fullfile('Z:\home\ziye\2P_Analysis\2P_Analysis',[date{id}, '_', runID{rID}, '_', mouseID{id}],'\');
        %         dest_sub = ['Z:\home\jake\Analysis\2P Analysis\Ziye_2P_figure\', date{id}, '_', runID{rID}, '_', mouseID{id}, '\'];
        if exist(dest_sub)
            load([dest_sub '_cell_TCs.mat']);
            load([dest_sub '_cell_resp.mat']);
            %             load([dest_sub '_cell_resp_amp.mat']);
            load([dest_sub '_cell_categories.mat']); % use only press and release cells
            load([dest_sub '_release_movies.mat']);
            load([dest_sub 'parse_behavior.mat']);
            load([dest_sub '_spike_variance.mat']);
            
            ncells{id} = size(press_resp,2);
            
            RS_cells{id} = unique([release_resp_cells success_resp_cells fail_resp_cells press_resp_cells tooFast_resp_cells]);
            
            if id <=20
                RS_cells_fail{id} = setdiff(RS_cells{id}, fake_cell);
            end
            
            tot_resp{id} = length(RS_cells{id});
            tot_success_cell{id} = length(success_resp_cells);
            tot_early_cell{id} = length(fail_resp_cells);
            tot_press_cell{id} = length(press_resp_cells);
            success_cell_only{id} = success_only_cells;
            success_fail_cell{id} = fail_and_success_cells;
            
            TC_length{id} = size(avg_release,2);
            TC_ifi{id} = ifi;
            pre_frames{id} = pre_release_frames;
            post_frames{id} = post_release_frames;
            pre_framesL{id} = pre_release_frames2;
            post_framesL{id} = post_release_frames2;
            
            RL_cells{id} = release_resp_cells;
            tot_RL_resp{id} = length(RL_cells{id});
            RL_trans_cells{id} = trans_cell_ind;
            RL_sustain_cells{id} = sustain_cell_ind;
            RL_trans_cellNum(id) = length(trans_cell_ind);
            RL_sustain_cellNum(id) = length(sustain_cell_ind);
            
            press_trans_cells(id) = length(intersect(press_resp_cells, trans_cell_ind'));
            press_sus_cells(id) = length(intersect(press_resp_cells, sustain_cell_ind'));
            
            release_resp_all{id} = mean((release_resp-release_base),1);
            release_resp_mean{id} = mean(release_resp_all{id},2);
            release_resp_sem{id} = std(release_resp_all{id},[],2)./sqrt(ncells{id});
            release_resp_RL{id} = mean((release_resp(:,release_resp_cells)-release_base(:,release_resp_cells)),1);
            release_resp_RL_mean{id} = mean(release_resp_RL{id},2);
            release_resp_RL_sem{id} = std(release_resp_RL{id},[],2)./sqrt(length(release_resp_cells));
            release_resp_RS{id} = mean((release_resp(:,RS_cells{id})-release_base(:,RS_cells{id})),1);
            release_resp_RS_mean{id} = mean(release_resp_RS{id},2);
            release_resp_RS_sem{id} = std(release_resp_RS{id},[],2)./sqrt(length(RS_cells{id}));
            %             release_resp_peak_mean{id} = mean(mean(release_peak));
            
            press_TC{id} = avg_press;
            press_TC_RS{id} = avg_press(RS_cells{id},:);
            press_TC_RS_L{id} = avg_press_long(RS_cells{id},:);
            press_TC_RL{id} = avg_press(RL_cells{id},:);
            press_TC_mean{id} = mean(avg_press,1);
            press_TC_sem{id} = std(avg_press,1)./sqrt(size(avg_press,1));
            press_TC_RS_mean{id} = mean(avg_press(RS_cells{id},:),1);
            press_TC_RS_sem{id} = std(avg_press(RS_cells{id},:),1)./sqrt(size(RS_cells{id},2));
            press_TC_RL_mean{id} = mean(avg_press(release_resp_cells,:),1);
            press_TC_RL_sem{id} = std(avg_press(release_resp_cells,:),1)./sqrt(size(release_resp_cells,2));
            press_resp_all{id} = mean((press_resp-press_base),1);
            press_resp_mean{id} = mean(press_resp_all{id},2);
            press_resp_sem{id} = std(press_resp_all{id},[],2)./sqrt(ncells{id});
            press_resp_RL{id} = mean((press_resp(:,release_resp_cells)-press_base(:,release_resp_cells)),1);
            press_resp_RL_mean{id} = mean(press_resp_RL{id},2);
            press_resp_RL_sem{id} = std(press_resp_RL{id},[],2)./sqrt(length(release_resp_cells));
            press_resp_RS{id} = mean((press_resp(:,RS_cells{id})-press_base(:,RS_cells{id})),1);
            press_resp_RS_mean{id} = mean(press_resp_RS{id},2);
            press_resp_RS_sem{id} = std(press_resp_RS{id},[],2)./sqrt(length(RS_cells{id}));
            %             press_resp_peak_mean{id} = mean(mean(press_peak));
            
            success_TC{id} = avg_success;
            success_TC_RS{id} = avg_success(RS_cells{id},:);
            success_TC_RS_L{id} = avg_success_long(RS_cells{id},:);
            success_TC_RL{id} = avg_success(RL_cells{id},:);
            success_TC_trans{id} = avg_success(trans_cell_ind,:);
            success_TC_sus{id} = avg_success(sustain_cell_ind,:);
            success_TC_mean{id} = mean(avg_success,1);
            success_TC_sem{id} = std(avg_success,1)./sqrt(size(avg_success,1));
            success_TC_RS_mean{id} = mean(avg_success(RS_cells{id},:),1);
            success_TC_RS_sem{id} = std(avg_success(RS_cells{id},:),1)./sqrt(size(RS_cells{id},2));
            success_TC_RL_mean{id} = mean(avg_success(release_resp_cells,:),1);
            success_TC_RL_sem{id} = std(avg_success(release_resp_cells,:),1)./sqrt(size(release_resp_cells,2));
            success_resp_all{id} = mean((success_resp-success_base),1);
            success_resp_mean{id} = mean(success_resp_all{id},2);
            success_resp_sem{id} = std(success_resp_all{id},[],2)./sqrt(ncells{id});
            success_resp_RS{id} = mean((success_resp(:,RS_cells{id})-success_base(:,RS_cells{id})),1);
            success_resp_RS_mean{id} = mean(success_resp_RS{id},2);
            success_resp_RS_sem{id} = std(success_resp_RS{id},[],2)./sqrt(length(RS_cells{id}));
            success_resp_RL{id} = mean((success_resp(:,release_resp_cells)-success_base(:,release_resp_cells)),1);
            success_resp_RL_mean{id} = mean(success_resp_RL{id},2);
            success_resp_RL_sem{id} = std(success_resp_RL{id},[],2)./sqrt(length(release_resp_cells));
            success_resp_RS_succonly{id} = mean((success_resp(:,success_only_cells)-success_base(:,success_only_cells)),1);
            success_resp_RS_succfail{id} = mean((success_resp(:,fail_and_success_cells)-success_base(:,fail_and_success_cells)),1);
            success_resp_SC{id} = mean((success_resp(:,success_resp_cells)-success_base(:,success_resp_cells)),1);
            success_resp_SC_mean{id} = mean(success_resp_SC{id},2);
            success_resp_SC_sem{id} = std(success_resp_SC{id},[],2)./sqrt(length(success_resp_cells));
            success_resp_RL_trans{id} = mean((success_resp(:,trans_cell_ind)-success_base(:,trans_cell_ind)),1);
            success_resp_RL_sus{id} = mean((success_resp(:,sustain_cell_ind)-success_base(:,sustain_cell_ind)),1);
            
            %             success_resp_peak_mean{id} = mean(mean(success_peak));
            
            tooFast_TC{id} = avg_tooFast;
            tooFast_TC_RS{id} = avg_tooFast(RS_cells{id},:);
            tooFast_TC_RS_L{id} = avg_tooFast_long(RS_cells{id},:);
            tooFast_TC_RL{id} = avg_tooFast(RL_cells{id},:);
            tooFast_resp_all{id} = mean((tooFast_resp-tooFast_base),1);
            tooFast_resp_RS{id} = mean((tooFast_resp(:,RS_cells{id})-tooFast_base(:,RS_cells{id})),1);
            tooFast_resp_RL{id} = mean((tooFast_resp(:,release_resp_cells)-tooFast_base(:,release_resp_cells)),1);
            
            fail_TC{id} = avg_fail;
            fail_TC_RS{id} = avg_fail(RS_cells{id},:);
            fail_TC_RS_L{id} = avg_fail_long(RS_cells{id},:);
            fail_TC_RL{id} = avg_fail(RL_cells{id},:);
            fail_TC_trans{id} = avg_fail(trans_cell_ind,:);
            fail_TC_sus{id} = avg_fail(sustain_cell_ind,:);
            fail_TC_mean{id} = mean(avg_fail,1);
            fail_TC_sem{id} = std(avg_fail,1)./sqrt(size(avg_fail,1));
            fail_TC_RS_mean{id} = mean(avg_fail(RS_cells{id},:),1);
            fail_TC_RS_sem{id} = std(avg_fail(RS_cells{id},:),1)./sqrt(size(RS_cells{id},2));
            fail_TC_RL_mean{id} = mean(avg_fail(release_resp_cells,:),1);
            fail_TC_RL_sem{id} = std(avg_fail(release_resp_cells,:),1)./sqrt(size(release_resp_cells,2));
            fail_resp_all{id} = mean((fail_resp-fail_base),1);
            fail_resp_mean{id} = mean(fail_resp_all{id},2);
            fail_resp_sem{id} = std(fail_resp_all{id},[],2)./sqrt(ncells{id});
            fail_resp_RS{id} = mean((fail_resp(:,RS_cells{id})-fail_base(:,RS_cells{id})),1);
            fail_resp_RS_mean{id} = mean(fail_resp_RS{id},2);
            fail_resp_RS_sem{id} = std(fail_resp_RS{id},[],2)./sqrt(length(RS_cells{id}));
            fail_resp_RL{id} = mean((fail_resp(:,release_resp_cells)-fail_base(:,release_resp_cells)),1);
            fail_resp_RL_mean{id} = mean(fail_resp_RL{id},2);
            fail_resp_RL_sem{id} = std(fail_resp_RL{id},[],2)./sqrt(length(release_resp_cells));
            fail_resp_SC{id} = mean((fail_resp(:,success_resp_cells)-fail_base(:,success_resp_cells)),1);
            fail_resp_SC_mean{id} = mean(fail_resp_SC{id},2);
            fail_resp_SC_sem{id} = std(fail_resp_SC{id},[],2)./sqrt(length(success_resp_cells));
            fail_resp_RL_trans{id} = mean((fail_resp(:,trans_cell_ind)-fail_base(:,trans_cell_ind)),1);
            fail_resp_RL_sus{id} = mean((fail_resp(:,sustain_cell_ind)-fail_base(:,sustain_cell_ind)),1);
            
            
            succ_hold_dur{id} = trial_outcome.succ_hold_dur'/1000;
            fail_hold_dur{id} = trial_outcome.fail_hold_dur'/1000;
            tf_hold_dur{id} = trial_outcome.tf_hold_dur'/1000;
            success_resp_RS_2{id} = mean((success_resp(:,RS_cells{id})-success_base(:,RS_cells{id})),2)';
            late_success_resp_RS_earlyWin{id} = mean((late_success_resp(:,RS_cells{id})-late_success_base(:,RS_cells{id})),1);
            fail_resp_RS_w1{id} = mean((fail_resp1(:,RS_cells{id})-fail_base(:,RS_cells{id})),2)';
            fail_resp_RS_w2{id} = mean((fail_resp2(:,RS_cells{id})-fail_base(:,RS_cells{id})),2)';
            fail_resp_RS_2{id} = mean((fail_resp(:,RS_cells{id})-fail_base(:,RS_cells{id})),2)';
            early_fail_resp_RS_earlyWin{id} = mean((early_fail_resp1(:,RS_cells{id})-early_fail_base(:,RS_cells{id})),1);
            late_fail_resp_RS_earlyWin{id} = mean((late_fail_resp1(:,RS_cells{id})-late_fail_base(:,RS_cells{id})),1);
            early_fail_resp_RS_lateWin{id} = mean((early_fail_resp2(:,RS_cells{id})-early_fail_base(:,RS_cells{id})),1);
            late_fail_resp_RS_lateWin{id} = mean((late_fail_resp2(:,RS_cells{id})-late_fail_base(:,RS_cells{id})),1);
            tooFast_resp_RS_2{id} = mean((tooFast_resp2(:,RS_cells{id})-tooFast_base(:,RS_cells{id})),2)';
            %                 late_tooFast_resp_RS_2{id} = mean((late_tooFast_resp(:,RS_cells{id})-late_tooFast_base(:,RS_cells{id})),2)';
            hold_dur_tot{id} = trial_outcome.req_holdTime(lever.cue2release>200);
            react_time_tot{id} = lever.cue2release(lever.cue2release>200);
            if id <= 20
                early_fail_TC_RS{id} = avg_early_fail(RS_cells_fail{id},:);
                late_fail_TC_RS{id} = avg_late_fail(RS_cells_fail{id},:);
            else
                early_fail_TC_RS{id} = avg_early_fail(RS_cells{id},:);
                late_fail_TC_RS{id} = avg_late_fail(RS_cells{id},:);
            end
            early_success_TC_RS{id} = avg_early_success(RS_cells{id},:);
            late_success_TC_RS{id} = avg_late_success(RS_cells{id},:);
            
            
            
            peak2release_1{id} = peak2release_std1;
            peak2release_2{id} = peak2release_std2;
            peak2cue_1{id} = peak2cue_std1;
            peak2cue_2{id} = peak2cue_std2;
            
        end
    end
end
col_mat_s = [0.5 0.5 0.5];

total_cells = sum(cell2mat(ncells));
total_resp = sum(cell2mat(tot_resp));
total_RL_resp = sum(cell2mat(tot_RL_resp));

% save([out_base, 'cell_count.mat'], 'total_cells', 'total_resp', 'total_RL_resp');

% stats for paper
% correct and correct&early cells resp to correct trials
a = cell2mat(success_resp_RS_succonly);
b = cell2mat(success_resp_RS_succfail);
disp(['mean of success_resp_RS_succonly is ', num2str(mean(a)), ' and SE is ', num2str(std(a)/sqrt(length(a))), ' N = ', num2str(length(a))]);
disp(['mean of success_resp_RS_succfail is ', num2str(mean(b)), ' and SE is ', num2str(std(b)/sqrt(length(b))), ' N = ', num2str(length(b))]);
[h,p,~,stats] = ttest2(a,b)

% col_mat = strvcat('r', 'b', 'r', 'b', 'g', 'm', 'c'); %hardcoded
col_mat = [ 0.9  0.9  0;
    1  0  1;
    0  1  1;
    0.5  0  0;
    0  1  0;
    0  0  1;
    1  0.6  1;
    0  0  0;
    1  0.8 0.4
    0  0.5 0.7
    0.5 0.4 0; 0.5 0.5 0.5; 0.3 0.5 1; 0.1 0.5 0.7; 0 0.6 0.2;0.8 0.8 0.4;0.1 0.1 0.1;0.3 0.7 0; 0 0 0; 0.1 0.5 0.5];

% scatter plot of dff for early early and late early using late window
fig = figure;
[h_elf2, p_elf2, stats_elf2] = scatter_plot(mouseID, early_fail_resp_RS_lateWin, late_fail_resp_RS_lateWin, col_mat_s);
hold on
xlim([-0.2 0.5]); ylim([-0.2 0.5])
x = -0.2:0.1:0.5; y=x;
plot(x,y,'-k')
vline(0,'--k')
hline(0,'--k')
xlabel('early < 1s dF/F')
ylabel('early > 3.5s dF/F')
title(['p value = ', num2str(p_elf2)]);
% title(['lick resp cells- n = ' num2str(total_lick_cell) ' out of total resp cells n = ' num2str(total_resp)])
saveas(fig, [out_base 'Summary_earlyLateEarlies_lateWindow_amp_scatter.fig']);
print([out_base 'Summary_earlyLateEarlies_lateWindow_amp_scatter.eps'], '-depsc');
print([out_base 'Summary_earlyLateEarlies_lateWindow_amp_scatter.pdf'], '-dpdf');

% scatter plot of dff for early early and late early using early window
fig = figure;
[h_elf1, p_elf1, stats_elf1] = scatter_plot(mouseID, early_fail_resp_RS_earlyWin, late_fail_resp_RS_earlyWin, col_mat_s);
hold on
xlim([-0.2 0.5]); ylim([-0.2 0.5])
x = -0.2:0.1:0.5; y=x;
plot(x,y,'-k')
vline(0,'--k')
hline(0,'--k')
xlabel('early < 1s dF/F')
ylabel('early > 3.5s dF/F')
title(['p value = ', num2str(p_elf1), ', early window' ]);
% title(['lick resp cells- n = ' num2str(total_lick_cell) ' out of total resp cells n = ' num2str(total_resp)])
saveas(fig, [out_base 'Summary_earlyLateEarlies_earlyWindow_amp_scatter_100165.fig']);
print([out_base 'Summary_earlyLateEarlies_earlyWindow_amp_scatter_100165.eps'], '-depsc');
print([out_base 'Summary_earlyLateEarlies_earlyWindow_amp_scatter_100165.pdf'], '-dpdf');

% scatter plot of dff for late correct and late early using early window
fig = figure;
[h_sfl, p_sfl, stats_sfl] = scatter_plot(mouseID, late_success_resp_RS_earlyWin, late_fail_resp_RS_earlyWin, col_mat_s);
hold on
xlim([-0.2 0.5]); ylim([-0.2 0.5])
x = -0.2:0.1:0.5; y=x;
plot(x,y,'-k')
vline(0,'--k')
hline(0,'--k')
xlabel('correct > 3.5s dF/F')
ylabel('early > 3.5s dF/F')
title(['p value = ', num2str(p_sfl), ', early window']);
% title(['lick resp cells- n = ' num2str(total_lick_cell) ' out of total resp cells n = ' num2str(total_resp)])
saveas(fig, [out_base 'Summary_LateEarliesCorrect_earlyWindow_amp_scatter_100165.fig']);
print([out_base 'Summary_LateEarliesCorrect_earlyWindow_amp_scatter_100165.eps'], '-depsc');
print([out_base 'Summary_LateEarliesCorrect_earlyWindow_amp_scatter_100165.pdf'], '-dpdf');


[tfHold_mean, tfResp_mean, tfResp_std, tIdx1, tIdx2] = avg_holdtime(tooFast_resp_RS_2, tf_hold_dur);
tfHold_mean = tfHold_mean(2:end);
tfResp_mean = tfResp_mean(2:end); 
tfResp_std = tfResp_std(2:end);
tfHold = (0.25*double(tIdx1)):0.25:0.25*double(tIdx2);
tfHold = tfHold(1:end)';
[LinearCoeffTF, fitTF] = polyfit(tfHold, tfResp_mean, 1);
CorrTFfit = polyval(LinearCoeffTF, tfHold);

fig = figure;
hold on
errorbar(tfHold, tfResp_mean, tfResp_std, '.', 'MarkerSize', 10, 'color', [0.5,0.5,0.5])
plot(tfHold, CorrTFfit, 'color', [0.5, 0.5,0.5], 'Linewidth', 1.5);

title('dF/F and hold time for tooFast');
xlabel('hold time(s)');
ylabel('dF/F');
saveas(fig, [out_base 'Summary_tooFast_avgdFF_holdtime_avgSession.fig']);
print([out_base 'Summary_tooFast_avgdFF_holdtime_avgSession.eps'], '-depsc');
print([out_base 'Summary_tooFast_avgdFF_holdtime_avgSession.pdf'], '-dpdf');

% col_mat_s = [0.5 0.5 0.5];%repmat([0.5 0.5 0.5], id, 1);
%
% sResp = cell2mat(success_resp_RS_2);
% sResp(isnan(sResp)) = [];
% sHold = cell2mat(succ_hold_dur);
% sHold(isnan(sHold)) = [];
% [LinearCoeffS1, fitS1] = polyfit(sHold, sResp, 1);
% CorrSfit1 = polyval(LinearCoeffS1, sHold);
% [BS,idx] = histc(sHold,0:0.25:round(max(sHold)));
% sHold_bins = accumarray(idx(:),sHold,[],@mean);
% sResp_std = accumarray(idx(:),sResp,[],@sem);
% sResp_bins = accumarray(idx(:),sResp,[],@mean);
% [LinearCoeffS2, fitS2] = polyfit(sHold_bins(BS>=10), sResp_bins(BS>=10), 1);
% CorrSfit2 = polyval(LinearCoeffS2, sHold_bins(BS>=10));
%
% fResp = cell2mat(fail_resp_RS_w1);
% fHold = cell2mat(fail_hold_dur);
% fHold(isnan(fResp)) = [];
% fResp(isnan(fResp)) = [];
% [LinearCoeffF1, fitF1] = polyfit(fHold, fResp, 1);
% CorrFfit1 = polyval(LinearCoeffF1, fHold);
% [BF,idx] = histc(fHold,0:0.25:round(max(fHold)));
% fHold_bins = accumarray(idx(:),fHold,[],@mean);
% fResp_std = accumarray(idx(:),fResp,[],@sem);
% fResp_bins = accumarray(idx(:),fResp,[],@mean);
% [LinearCoeffF2, fitF2] = polyfit(fHold_bins(BF>=10), fResp_bins(BF>=10), 1);
% CorrFfit2 = polyval(LinearCoeffF2, fHold_bins(BF>=10));
%
%
% fig = figure;
% subplot(1,2,1); axis square; ylim([-0.05 0.2])
% % scatter_plot(mouseID, succ_hold_dur, success_resp_RS_2, col_mat_s);
% hold on
% errorbar(sHold_bins(BS>=10), sResp_bins(BS>=10), sResp_std(BS>=10), '.', 'MarkerSize', 10, 'color', [0.5,0.5,0.5])
% plot(sHold_bins(BS>=10), CorrSfit2, 'color', [0.5, 0.5,0.5], 'Linewidth', 1.5);
%
% errorbar(fHold_bins(BF>=10), fResp_bins(BF>=10), fResp_std(BF>=10), '.', 'MarkerSize', 10, 'color', [0.9,0,0])
% plot(fHold_bins(BF>=10), CorrFfit2, 'color', [0.9, 0,0], 'Linewidth', 1.5);
% title('window 1, black- correct, red- early');
% xlabel('hold time(s)');
% ylabel('dF/F');
%
%
% fResp = cell2mat(fail_resp_RS_w2);
% fHold = cell2mat(fail_hold_dur);
% fHold(isnan(fResp)) = [];
% fResp(isnan(fResp)) = [];
% [LinearCoeffF1, fitF1] = polyfit(fHold, fResp, 1);
% CorrFfit1 = polyval(LinearCoeffF1, fHold);
% [BF,idx] = histc(fHold,0:0.25:round(max(fHold)));
% fHold_bins = accumarray(idx(:),fHold,[],@mean);
% fResp_std = accumarray(idx(:),fHold,[],@sem);
% fResp_bins = accumarray(idx(:),fResp,[],@mean);
% fResp_bins(13) = abs(fResp_bins(13));
% [LinearCoeffF2, fitF2] = polyfit(fHold_bins(BF>=10), fResp_bins(BF>=10), 1);
% CorrFfit2 = polyval(LinearCoeffF2, fHold_bins(BF>=10));
%
% subplot(1,2,2); axis square; ylim([-0.05 0.2])
% % scatter_plot(mouseID, succ_hold_dur, success_resp_RS_2, col_mat_s);
% hold on
% errorbar(sHold_bins(BS>=10), sResp_bins(BS>=10), sResp_std(BS>=10), '.', 'MarkerSize', 10, 'color', [0.5,0.5,0.5])
% plot(sHold_bins(BS>=10), CorrSfit2, 'color', [0.5, 0.5,0.5], 'Linewidth', 1.5);
%
% errorbar(fHold_bins(BF>=10), fResp_bins(BF>=10), fResp_std(BF>=10), '.', 'MarkerSize', 10, 'color', [0.9,0,0])
% plot(fHold_bins(BF>=10), CorrFfit2, 'color', [0.9, 0,0], 'Linewidth', 1.5);
%
% %plot(x,y,'-k')
% title('window 2, black-correct, red- early');
% xlabel('hold time(s)');
% ylabel('dF/F');
% saveas(fig, [out_base 'Summary_earlylateEarly_holdtime_all20.fig']);
% print([out_base 'Summary_earlylateEarly_holdtime_all20.eps'], '-depsc');
% print([out_base 'Summary_earlylateEarly_holdtime_all20.pdf'], '-dpdf');
%
%
%
%
% %summary of average response amplitudes
% fig=figure;
% %resp amplitudes- press/release
% subplot(3,2,1)
% errorbar(1:size(mouseID,2), cell2mat(release_resp_mean), cell2mat(release_resp_sem), 'ok')
% hold on
% errorbar(1:size(mouseID,2), cell2mat(press_resp_mean), cell2mat(press_resp_sem), 'oc')
% xlabel('Expt #')
% ylabel('dF/F')
% xlim([0 size(mouseID,2)+1])
% title('Press/release amplitude')
% %resp amplitudes- success/fail
% subplot(3,2,2)
% errorbar(1:size(mouseID,2), cell2mat(success_resp_mean), cell2mat(success_resp_sem), 'ok')
% hold on
% errorbar(1:size(mouseID,2), cell2mat(fail_resp_mean), cell2mat(fail_resp_sem), 'or')
% xlabel('Expt #')
% ylabel('dF/F')
% xlim([0 size(mouseID,2)+1])
% title('Success/fail amplitude')
%
% %resp amplitudes- press/release - driven by press or release
% subplot(3,2,3)
% errorbar(1:size(mouseID,2), cell2mat(release_resp_RS_mean), cell2mat(release_resp_RS_sem), 'ok')
% hold on
% errorbar(1:size(mouseID,2), cell2mat(press_resp_RS_mean), cell2mat(press_resp_RS_sem), 'oc')
% xlabel('Expt #')
% ylabel('dF/F')
% xlim([0 size(mouseID,2)+1])
% title('Press/release: driven')
% %resp amplitudes- success/fail - driven press or release
% subplot(3,2,4)
% errorbar(1:size(mouseID,2), cell2mat(success_resp_RS_mean), cell2mat(success_resp_RS_sem), 'ok')
% hold on
% errorbar(1:size(mouseID,2), cell2mat(fail_resp_RS_mean), cell2mat(fail_resp_RS_sem), 'or')
% xlabel('Expt #')
% ylabel('dF/F')
% xlim([0 size(mouseID,2)+1])
% title('Success/fail: driven ')
%
% %resp amplitudes- success/fail - driven by release
% subplot(3,2,5)
% errorbar(1:size(mouseID,2), cell2mat(success_resp_RL_mean), cell2mat(success_resp_RL_sem), 'ok')
% hold on
% errorbar(1:size(mouseID,2), cell2mat(fail_resp_RL_mean), cell2mat(fail_resp_RL_sem), 'or')
% xlabel('Expt #')
% ylabel('dF/F')
% xlim([0 size(mouseID,2)+1])
% title('Success/fail driven by release')
% %resp amplitudes- success/fail - driven by fail
% subplot(3,2,6)
% errorbar(1:size(mouseID,2), cell2mat(success_resp_SC_mean), cell2mat(success_resp_SC_sem), 'ok')
% hold on
% errorbar(1:size(mouseID,2), cell2mat(fail_resp_SC_mean), cell2mat(fail_resp_SC_sem), 'or')
% xlabel('Expt #')
% ylabel('dF/F')
% xlim([0 size(mouseID,2)+1])
% title('Success/fail driven by success')
%
% supertitle(['Summary of cell response amplitudes'])
% saveas(fig, [out_base 'Summary_cell_response_amp.fig']);
% print([out_base 'Summary_cell_response_amp.eps'], '-depsc');
% print([out_base 'Summary_cell_response_amp.pdf'], '-dpdf');
%
%scatter of response amplitudes
fig=figure;
x = [-0.27:.01:1];
y = x;
r = [];
p = [];

subplot(3,3,1)
[h_prall, p_prall] = scatter_plot(mouseID, release_resp_all, press_resp_all, col_mat_s);
hold on
plot(x,y,'-k')
xlim([-.27 0.91]);
ylim([-.27 0.91]);
xlabel('Release dF/F')
ylabel('Press dF/F')
hold on
vline(0,'--k')
hline(0,'--k')
title(['All cells- p = ' num2str(p_prall)])

subplot(3,3,2)
[h_prRS, p_prRS, stats_prRS]=scatter_plot(mouseID, release_resp_RS, press_resp_RS, col_mat_s);
hold on
plot(x,y,'-k')
xlim([-.27 0.91]);
ylim([-.27 0.91]);
xlabel('Release dF/F')
ylabel('Press dF/F')

hold on
vline(0,'--k')
hline(0,'--k')
title(['Resp cells- p = ' num2str(p_prRS)])
subplot(3,3,3)
[h_prRL, p_prRL] = scatter_plot(mouseID, release_resp_RL, press_resp_RL, col_mat_s);
hold on
plot(x,y,'-k')
xlim([-.27 0.91]);
ylim([-.27 0.91]);
xlabel('Release dF/F')
ylabel('Press dF/F')
hold on
vline(0,'--k')
hline(0,'--k')
title(['Release resp cells- p = ' num2str(p_prRL)])

subplot(3,3,4)
[h_sfall, p_sfall] = scatter_plot(mouseID, success_resp_all, fail_resp_all, col_mat_s);
hold on
xlim([-.27 0.91]);
ylim([-.27 0.91]);
plot(x,y,'-k')
xlabel('Correct dF/F')
ylabel('Early dF/F')

hold on
vline(0,'--k')
hline(0,'--k')
title(['All cells- p = ' num2str(chop(p_sfall,2))])
subplot(3,3,5)
[h_sfRS, p_sfRS, stats_sfRS] = scatter_plot(mouseID, success_resp_RS, fail_resp_RS, col_mat_s);
hold on
plot(x,y,'-k')
xlim([-.27 0.91]);
ylim([-.27 0.91]);
xlabel('Correct dF/F')
ylabel('Early dF/F')

hold on
vline(0,'--k')
hline(0,'--k')
title(['Resp cells- p = ' num2str(chop(p_sfRS,2))])

subplot(3,3,6)
[h_sfRL, p_sfRL] = scatter_plot(mouseID, success_resp_RL, fail_resp_RL, col_mat_s);
hold on
plot(x,y,'-k')
xlim([-.27 0.91]);
ylim([-.27 0.91]);
xlabel('Correct dF/F')
ylabel('Early dF/F')

hold on
vline(0,'--k')
hline(0,'--k')
title(['Release resp cells- p = ' num2str(chop(p_sfRL,2))])

subplot(3,3,7)
[h_stall, p_sfall, stats_sfall] = scatter_plot(mouseID, success_resp_all, tooFast_resp_all, col_mat_s);
hold on
plot(x,y,'-k')
xlim([-.27 0.91]);
ylim([-.27 0.91]);
xlabel('Correct dF/F')
ylabel('Too Fast dF/F')

hold on
vline(0,'--k')
hline(0,'--k')
title(['All cells- p = ' num2str(chop(p_sfall,2))])

subplot(3,3,8)
[h_sfRS, p_sfRS, stats_sfRS] = scatter_plot(mouseID, success_resp_RS, tooFast_resp_RS, col_mat_s);
hold on
plot(x,y,'-k')
xlim([-.27 0.91]);
ylim([-.27 0.91]);
xlabel('Correct dF/F')
ylabel('Too Fast dF/F')

hold on
vline(0,'--k')
hline(0,'--k')
title(['Resp cells- p = ' num2str(chop(p_sfRS,2))])

subplot(3,3,9)
[h_stRL, p_stRL, stats_stRL] = scatter_plot(mouseID, success_resp_RL, tooFast_resp_RL, [0.5 0.5 0.5]);
hold on
plot(x,y,'-k')
xlim([-.27 0.91]);
ylim([-.27 0.91]);
xlabel('Correct dF/F')
ylabel('Too Fast dF/F')

hold on
vline(0,'--k')
hline(0,'--k')
title(['Release resp cells- p = ' num2str(chop(p_stRL,2))])

supertitle(['Summary of cell response amplitudes'])

saveas(fig, [out_base 'Summary_cell_response_amp_scatter.fig']);
print([out_base 'Summary_cell_response_amp_scatter.eps'], '-depsc');
print([out_base 'Summary_cell_response_amp_scatter.pdf'], '-dpdf');
%
%
% %new figure for comparing transient and sustain response
% fig = figure;
% x = [-.05:.01:.3];
% y = x;
% subplot(1,2,1)
% scatter_plot(mouseID, success_resp_RL_trans, fail_resp_RL_trans, col_mat_s);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% axis square
% xlabel('Success dF/F')
% ylabel('Fail dF/F')
% title(['Transient response cells n=',num2str(sum(RL_trans_cellNum))])
% subplot(1,2,2)
% scatter_plot(mouseID, success_resp_RL_sus, fail_resp_RL_sus, col_mat_s);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% axis square
% xlabel('Success dF/F')
% ylabel('Fail dF/F')
% title(['Sustain response cells n=',num2str(sum(RL_sustain_cellNum))])
% saveas(fig, [out_base 'Summary_transient_sustain_resp_amp.fig']);
% print([out_base 'Summary_transient_sustain_resp_amp.eps'], '-depsc');
% print([out_base 'Summary_transient_sustain_resp_amp.pdf'], '-dpdf');
%
% %avg scatter of amplitudes by mouse
% fig=figure;
% x = [-.3:.01:1];
% y = x;
% % col_mat = strvcat('r', 'b', 'r', 'b', 'g', 'm', 'c');   %HARDCODED
% subplot(3,3,1)
% scatter_plot(mouseID, release_resp_all, press_resp_all, col_mat);
%
% hold on
% plot(x,y,'-k')
% xlim([-.05 .4]);
% ylim([-.05 .4]);
% xlabel('Release dF/F')
% ylabel('Press dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['All cells'])
% subplot(3,3,2)
% scatter_plot(mouseID, release_resp_RS, press_resp_RS, col_mat);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .4]);
% ylim([-.05 .4]);
% xlabel('Release dF/F')
% ylabel('Press dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Resp cells'])
% subplot(3,3,3)
% [h_prRS, p_prRS]=scatter_plot(mouseID, release_resp_RS, press_resp_RS, col_mat_s);
% hold on
% plot(x,y,'-k')
% xlim([-.27 .91]);
% ylim([-.27 .91]);
% vline(0,'--k')
% hline(0,'--k')
% xlabel('Release dF/F')
% ylabel('Press dF/F')
% title(['Resp cells'])
% % subplot(3,3,3)
% % scatter_plot(mouseID, release_resp_RL, press_resp_RL, col_mat);
% % hold on
% % plot(x,y,'-k')
% % xlim([-.05 .3]);
% % ylim([-.05 .3]);
% % xlabel('Release dF/F')
% % ylabel('Press dF/F')
% % hold on
% % vline(0,'--k')
% % hline(0,'--k')
% % title(['Release resp cells'])
% subplot(3,3,4)
% scatter_plot(mouseID, success_resp_all, fail_resp_all, col_mat);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .4]);
% ylim([-.05 .4]);
% xlabel('Correct dF/F')
% ylabel('Early dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['All cells'])
% % subplot(3,3,5)
% % scatter_plot(mouseID, success_resp_RS, fail_resp_RS, col_mat);
% % hold on
% % plot(x,y,'-k')
% % xlim([-.05 .3]);
% % ylim([-.05 .3]);
% % xlabel('Success dF/F')
% % ylabel('Fail dF/F')
% % hold on
% % vline(0,'--k')
% % hline(0,'--k')
% % title(['Resp cells'])
% subplot(3,3,5)
% scatter_plot(mouseID, success_resp_RL, fail_resp_RL, col_mat);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .4]);
% ylim([-.05 .4]);
% xlabel('Correct dF/F')
% ylabel('Early dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Release resp cells'])
%
% subplot(3,3,6)
% [h_sfRL, p_sfRL] = scatter_plot(mouseID, success_resp_RL, fail_resp_RL, col_mat_s);
% hold on
% plot(x,y,'-k')
% xlim([-.27 .91]);
% ylim([-.27 .91]);
% xlabel('Correct dF/F')
% ylabel('Early dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Release resp cells'])
% % title(['Release resp cells- p = ' num2str(chop(p_sfRL,2))])
%
% subplot(3,3,7)
% scatter_plot(mouseID, success_resp_all, tooFast_resp_all, col_mat);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .4]);
% ylim([-.05 .4]);
% xlabel('Correct dF/F')
% ylabel('TooFast Correct dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['All cells'])
%
% % subplot(3,3,8)
% % scatter_plot(mouseID, success_resp_RS, tooFast_resp_RS, col_mat);
% % hold on
% % plot(x,y,'-k')
% % xlim([-.05 .3]);
% % ylim([-.05 .3]);
% % xlabel('Success dF/F')
% % ylabel('Toofast Correct dF/F')
% % hold on
% % vline(0,'--k')
% % hline(0,'--k')
% % title(['Resp cells'])
%
% subplot(3,3,8)
% scatter_plot(mouseID, success_resp_RL, tooFast_resp_RL, col_mat);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .4]);
% ylim([-.05 .4]);
% xlabel('Correct dF/F')
% ylabel('Toofast Correct dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Release resp cells'])
%
% subplot(3,3,9)
% [h_stRL, p_stRL] = scatter_plot(mouseID, success_resp_RL, tooFast_resp_RL, col_mat_s);
% hold on
% plot(x,y,'-k')
% xlim([-.27 .91]);
% ylim([-.27 .91]);
% xlabel('Correct dF/F')
% ylabel('TooFast Correct dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Release resp cells'])
% % title(['Release resp cells- p = ' num2str(chop(p_stRL,2))])
%
% supertitle(['Summary of cell response amplitudes'])
% saveas(fig, [out_base 'Summary_avg_response_amp_scatter.fig']);
% print([out_base 'Summary_avg_response_amp_scatter.eps'], '-depsc');
% print([out_base 'Summary_avg_response_amp_scatter.pdf'], '-dpdf');
%
% % plot avg response amp scatter for only fail/success and tooFast/success
% fig = figure;
% subplot(1,2,1)
% scatter_plot(mouseID, success_resp_RS, fail_resp_RS, col_mat);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .35]);
% ylim([-.05 .35]);
% xlabel('Success dF/F')
% ylabel('Fail dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Resp cells'])
% subplot(1,2,2)
% scatter_plot(mouseID, success_resp_RL, tooFast_resp_RL, col_mat);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .35]);
% ylim([-.05 .35]);
% xlabel('Success dF/F')
% ylabel('Toofast Correct dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Release resp cells'])
%
% supertitle(['Summary of cell response amplitudes'])
% saveas(fig, [out_base 'Summary_avg_response_amp_scatter_2.fig']);
% print([out_base 'Summary_avg_response_amp_scatter_2.eps'], '-depsc');
% print([out_base 'Summary_avg_response_amp_scatter_2.pdf'], '-dpdf');
% %average timecourse for expts
%
% fig=figure;
% for id = 1:size(mouseID,2)
%     subplot(5,5,id)
%     tt =((-pre_frames{id}:post_frames{id}).*double(TC_ifi{id}))./1000;
%     shadedErrorBar(tt, success_TC_mean{id},success_TC_sem{id}, 'k');
%     hold on
%     shadedErrorBar(tt, fail_TC_mean{id},fail_TC_sem{id}, 'r');
%     hold on
%     shadedErrorBar(tt, press_TC_mean{id},press_TC_sem{id}, 'c');
%     xlim([-pre_frames{id}.*double(TC_ifi{id})./1000 post_frames{id}.*double(TC_ifi{id})./1000])
%     xlabel('Time (ms)')
%     ylabel('dF/F')
%     title([date{id} ' ' mouseID{id}])
% end
% supertitle(['Summary of all cell timecourses'])
% saveas(fig, [out_base 'Summary_allexptTCs_allcells.fig']);
% print([out_base 'Summary_allexptTCs_allcells.eps'], '-depsc');
% print([out_base 'Summary_allexptTCs_allcells.pdf'], '-dpdf');
%
% fig=figure;
% for id = 1:size(mouseID,2)
%     subplot(5,5,id)
%     tt =((-pre_frames{id}:post_frames{id}).*double(TC_ifi{id}))./1000;
%     shadedErrorBar(tt, success_TC_RS_mean{id},success_TC_RS_sem{id}, 'k');
%     hold on
%     shadedErrorBar(tt, fail_TC_RS_mean{id},fail_TC_RS_sem{id}, 'r');
%     hold on
%     shadedErrorBar(tt, press_TC_RS_mean{id},press_TC_RS_sem{id}, 'c');
%     xlim([-pre_frames{id}.*double(TC_ifi{id})./1000 post_frames{id}.*double(TC_ifi{id})./1000])
%     xlabel('Time (ms)')
%     ylabel('dF/F')
%     title([date{id} ' ' mouseID{id}])
% end
% supertitle(['Summary of responsive cell timecourses'])
% saveas(fig, [out_base 'Summary_allexptTCs_respcells.fig']);
% print([out_base 'Summary_allexptTCs_respcells.eps'], '-depsc');
% print([out_base 'Summary_allexptTCs_respcells.pdf'], '-dpdf');
%
% fig=figure;
% for id = 1:size(mouseID,2)
%     subplot(5,5,id)
%     tt =((-pre_frames{id}:post_frames{id}).*double(TC_ifi{id}))./1000;
%     shadedErrorBar(tt, success_TC_RL_mean{id},success_TC_RL_sem{id}, 'k');
%     hold on
%     shadedErrorBar(tt, fail_TC_RL_mean{id},fail_TC_RL_sem{id}, 'r');
%     hold on
%     shadedErrorBar(tt, press_TC_RL_mean{id},press_TC_RL_sem{id}, 'c');
%     xlim([-pre_frames{id}.*double(TC_ifi{id})./1000 post_frames{id}.*double(TC_ifi{id})./1000])
%     xlabel('Time (ms)')
%     ylabel('dF/F')
%     title([date{id} ' ' mouseID{id}])
% end
% supertitle(['Summary of release responsive cell timecourses'])
% saveas(fig, [out_base 'Summary_allexptTCs_relcells.fig']);
% print([out_base 'Summary_allexptTCs_relcells.eps'], '-depsc');
% print([out_base 'Summary_allexptTCs_relcells.pdf'], '-dpdf');

%% commented for now until a decision is made on how to average across experiments with different acquisition rates
%averaging across all cells- specific to different acquisition rates

frame_size = cell2mat(cellfun(@size, success_TC, 'UniformOutput', 0));
max_frame = max(frame_size(2:2:end));

success_TC_all = interp_frame(success_TC, max_frame);
fail_TC_all    = interp_frame(fail_TC, max_frame);
press_TC_all   = interp_frame(press_TC, max_frame);
tooFast_TC_all = interp_frame(tooFast_TC, max_frame);

success_TC_all_RS = interp_frame(success_TC_RS, max_frame);
early_success_TC_all_RS = interp_frame(early_success_TC_RS, max_frame);
late_success_TC_all_RS = interp_frame(late_success_TC_RS, max_frame);
fail_TC_all_RS    = interp_frame(fail_TC_RS, max_frame);
early_fail_TC_all_RS = interp_frame(early_fail_TC_RS, max_frame);
late_fail_TC_all_RS = interp_frame(late_fail_TC_RS, max_frame);
press_TC_all_RS   = interp_frame(press_TC_RS, max_frame);
tooFast_TC_all_RS = interp_frame(tooFast_TC_RS, max_frame);
early_tooFast_TC_all_RS = interp_frame(early_tooFast_TC_RS, max_frame);
late_tooFast_TC_all_RS = interp_frame(late_tooFast_TC_RS, max_frame);

success_TC_all_RL = interp_frame(success_TC_RL, max_frame);
fail_TC_all_RL    = interp_frame(fail_TC_RL, max_frame);
press_TC_all_RL   = interp_frame(press_TC_RL, max_frame);
tooFast_TC_all_RL = interp_frame(tooFast_TC_RL, max_frame);

frame_sizeL = cell2mat(cellfun(@size, success_TC_RS_L, 'UniformOutput', 0));
max_frameL = max(frame_sizeL(2:2:end));

success_TC_all_RS_L = interp_frame(success_TC_RS_L, max_frameL);
fail_TC_all_RS_L    = interp_frame(fail_TC_RS_L, max_frameL);
press_TC_all_RS_L   = interp_frame(press_TC_RS_L, max_frameL);
tooFast_TC_all_RS_L = interp_frame(tooFast_TC_RS_L, max_frameL);

pre_f = double(cell2mat(cellfun(@int64,pre_frames,'UniformOutput',0)));
post_f = double(cell2mat(cellfun(@int64,post_frames,'UniformOutput',0)));
ifi_f = double(cell2mat(cellfun(@int64,TC_ifi,'UniformOutput',0)));

tt =(-max(pre_f):max(post_f)).*double(min(ifi_f))./1000;
% ttL =(-max(cell2mat(pre_framesL)):max(cell2mat(post_framesL))).*double(min(cell2mat(TC_ifi)))./1000;

% plot early early and late early TCs
fig = figure;
hold on;
shadedErrorBar(tt, mean(early_fail_TC_all_RS,1)-0.005, std(early_fail_TC_all_RS,[],1)./sqrt(size(early_fail_TC_all_RS,1)), 'm');
shadedErrorBar(tt, mean(late_fail_TC_all_RS,1)+0.007, std(late_fail_TC_all_RS,[],1)./sqrt(size(late_fail_TC_all_RS,1)), 'r');
title(['Early early TC- magenta ', 'n=', num2str(size(early_fail_TC_all_RS,1)), ' and late early TC- red', ' n=', num2str(size(late_fail_TC_all_RS,1))]);
xlabel('Time (s)')
ylabel('dF/F')
xlim([-0.5 1])
saveas(fig, [out_base 'Summary_avgEarlyLateEarlies_TCs.fig']);
print([out_base 'Summary_avgEarlyLateEarlies_TCs.eps'], '-depsc');
print([out_base 'Summary_avgEarlyLateEarlies_TCs.pdf'], '-dpdf');


% plot early correct and late correct TCs
fig = figure;
shadedErrorBar(tt, mean(early_success_TC_all_RS,1), std(early_success_TC_all_RS,[],1)./sqrt(size(early_success_TC_all_RS,1)), 'b');
hold on;
shadedErrorBar(tt, mean(late_success_TC_all_RS,1), std(late_success_TC_all_RS,[],1)./sqrt(size(late_success_TC_all_RS,1)), 'k');
title('Early correct TC- blue and late correct TC- black');
xlabel('Time (s)')
ylabel('dF/F')
xlim([-0.5 1]); ylim([-0.05 0.3])
saveas(fig, [out_base 'Summary_avgEarlyLateCorrects_TCs_all20.fig']);
print([out_base 'Summary_avgEarlyLateCorrects_TCs_all20.eps'], '-depsc');
print([out_base 'Summary_avgEarlyLateCorrects_TCs_all20.pdf'], '-dpdf');

% plot early toofast and late toofast TCs
fig = figure;
shadedErrorBar(tt, mean(early_tooFast_TC_all_RS,1), std(early_tooFast_TC_all_RS,[],1)./sqrt(size(early_tooFast_TC_all_RS,1)), 'b');
hold on;
shadedErrorBar(tt, mean(late_tooFast_TC_all_RS,1), std(late_tooFast_TC_all_RS,[],1)./sqrt(size(late_tooFast_TC_all_RS,1)), 'k');
title('Early tooFast TC- blue and late tooFast TC- black');
xlabel('Time (s)')
ylabel('dF/F')
xlim([-0.5 1]);ylim([-0.05 0.3])
saveas(fig, [out_base 'Summary_avgEarlyLateTooFast_TCs_all20.fig']);
print([out_base 'Summary_avgEarlyLateTooFast_TCs_all20.eps'], '-depsc');
print([out_base 'Summary_avgEarlyLateTooFast_TCs_all20.pdf'], '-dpdf');


fig=figure;
subplot(3,1,1)
shadedErrorBar(tt, mean(success_TC_all,1), std(success_TC_all,[],1)./sqrt(size(success_TC_all,1)), 'k');
hold on;
shadedErrorBar(tt, mean(fail_TC_all,1), std(fail_TC_all,[],1)./sqrt(size(fail_TC_all,1)), 'r');
hold on
shadedErrorBar(tt, mean(press_TC_all,1), std(press_TC_all,[],1)./sqrt(size(press_TC_all,1)), 'c');
shadedErrorBar(tt, nanmean(tooFast_TC_all,1), nanstd(tooFast_TC_all,[],1)./sqrt(size(tooFast_TC_all,1)), 'g');
% title(['All cells- n = ' num2str(size(press_TC_all,1))])
title(['All cells- n = ' num2str(total_cells)])
xlabel('Time (ms)')
ylabel('dF/F')

subplot(3,1,2)
shadedErrorBar(tt, mean(success_TC_all_RS,1), std(success_TC_all_RS,[],1)./sqrt(size(success_TC_all_RS,1)), 'k');
hold on;
shadedErrorBar(tt, mean(fail_TC_all_RS,1), std(fail_TC_all_RS,[],1)./sqrt(size(fail_TC_all_RS,1)), 'r');
hold on
shadedErrorBar(tt, mean(press_TC_all_RS,1), std(press_TC_all_RS,[],1)./sqrt(size(press_TC_all_RS,1)), 'c');
shadedErrorBar(tt, nanmean(tooFast_TC_all_RS,1), nanstd(tooFast_TC_all_RS,[],1)./sqrt(size(tooFast_TC_all_RS,1)), 'g');
% title(['Responsive cells- n = ' num2str(size(press_TC_all_RS,1))])
title(['All cells- n = ' num2str(total_resp)])
xlabel('Time (ms)')
ylabel('dF/F')

subplot(3,1,3)
shadedErrorBar(tt, mean(success_TC_all_RL,1), std(success_TC_all_RL,[],1)./sqrt(size(success_TC_all_RL,1)), 'k');
hold on;
shadedErrorBar(tt, mean(fail_TC_all_RL,1), std(fail_TC_all_RL,[],1)./sqrt(size(fail_TC_all_RL,1)), 'r');
hold on
shadedErrorBar(tt, mean(press_TC_all_RL,1), std(press_TC_all_RL,[],1)./sqrt(size(press_TC_all_RL,1)), 'c');
shadedErrorBar(tt, nanmean(tooFast_TC_all_RL,1), nanstd(tooFast_TC_all_RL,[],1)./sqrt(size(tooFast_TC_all_RL,1)), 'g');
% title(['Release responsive cells- n = ' num2str(size(press_TC_all_RL,1))])
title(['All cells- n = ' num2str(total_RL_resp)])
xlabel('Time (ms)')
ylabel('dF/F')
% supertitle('Average all cells collected at 30 Hz')
supertitle(['Timecourses of average all cells across experiments Black- Correct; Red- Early; Cyan- Press; Green- TooFast']);
saveas(fig, [out_base 'Summary_avgexpt_TCs.fig']);
print([out_base 'Summary_avgexpt_TCs.eps'], '-depsc');
print([out_base 'Summary_avgexpt_TCs.pdf'], '-dpdf');


%% plot dF/F vs hold time for correct and early
% fig = figure;
% col_mat_s = [0.5,0.5,0.5];
% scatter_plot(mouseID, succ_hold_dur,  success_resp_RS_2, col_mat_s);
%
% hold on
% col_mat_s = [0.9,0,0];
% scatter_plot(mouseID, fail_hold_dur, fail_resp_RS_2, col_mat_s);
% %plot(x,y,'-k')
% title('dF/F and hold time for responsive cells, each point is one session, black-correct, red-early');
% xlabel('hold time(ms)');
% ylabel('dF/F');
% saveas(fig, [out_base 'Summary_avgdFF_holdtime_session.fig']);
% print([out_base 'Summary_avgdFF_holdtime_session.eps'], '-depsc');
% print([out_base 'Summary_avgdFF_holdtime_session.pdf'], '-dpdf');
[sHold_mean, sResp_mean, sResp_std, sIdx1, sIdx2] = avg_holdtime(success_resp_RS_2, succ_hold_dur);
sHold = (0.25*double(sIdx1)):0.25:0.25*double(sIdx2);
sHold = sHold';
[fHold_mean, fResp_mean, fResp_std, fIdx1, fIdx2] = avg_holdtime(fail_resp_RS_2, fail_hold_dur);
fHold = (0.25*double(fIdx1)):0.25:0.25*double(fIdx2); 
fHold = fHold(2:end-2)';


[LinearCoeffS, fitS] = polyfit(sHold, sResp_mean, 1);
CorrSfit = polyval(LinearCoeffS, sHold);

[LinearCoeffF, fitF] = polyfit(fHold, fResp_mean, 1);
CorrFfit = polyval(LinearCoeffF, fHold);

% calculate fitting and R^2
yfit = CorrSfit; y = sResp_mean;
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal

yfit = CorrFfit; y = fResp_mean;
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal

fig = figure;
hold on
errorbar(sHold, sResp_mean, sResp_std, '.', 'MarkerSize', 10, 'color', [0.5,0.5,0.5])
plot(sHold, CorrSfit, 'color', [0.5, 0.5,0.5], 'Linewidth', 1.5);

errorbar(fHold, fResp_mean, fResp_std, '.', 'MarkerSize', 10, 'color', [0.9,0,0])
plot(fHold, CorrFfit, 'color', [0.9, 0, 0], 'Linewidth', 1.5);

title('dF/F and hold time for responsive cells, each point is one trial, black-correct, red-early');
xlabel('hold time(s)');
ylabel('dF/F');
saveas(fig, [out_base 'Summary_avgdFF_holdtime_avgSession.fig']);
print([out_base 'Summary_avgdFF_holdtime_avgSession.eps'], '-depsc');
print([out_base 'Summary_avgdFF_holdtime_avgSession.pdf'], '-dpdf');

fig = figure;
col_mat_s = [0.5, 0.5, 0.5];
scatter_plot(mouseID, succ_hold_dur, success_resp_RS_2, col_mat_s);
hold on
plot(sHold, CorrSfit1, 'color', [0.5, 0.5,0.5], 'Linewidth', 1.5);

col_mat_s = [0.9 0 0];
scatter_plot(mouseID, fail_hold_dur, fail_resp_RS_2, col_mat_s);
plot(fHold, CorrFfit1, 'color', [0.9, 0,0], 'Linewidth', 1.5);

%plot(x,y,'-k')
title('event rate and hold time for responsive cells, each point is one trial, black-correct, red-early');
xlabel('hold time(s)');
ylabel('dF/F');
saveas(fig, [out_base 'Summary_avgdFF_holdtime_trial_scatter.fig']);
print([out_base 'Summary_avgdFF_holdtime_trial_scatter.eps'], '-depsc');
print([out_base 'Summary_avgdFF_holdtime_trial_scatter.pdf'], '-dpdf');

col_mat_s = repmat([0.5,0.5,0.5],11,1);


% plot hold time vs react time
col_mat_s = [0.5, 0.5, 0.5];

fig = figure;
rt_tot = cell2mat(react_time_tot);
rt_tot(isnan(rt_tot)) = [];
hold_tot = cell2mat(hold_dur_tot);
hold_tot(isnan(hold_tot)) = [];
[hold_tot,idx] = sort(hold_tot); % sort just the first column

rt_tot = rt_tot(idx);   % sort the whole matrix using the sort indices
[LinearCoefft, fitT] = polyfit(hold_tot, rt_tot, 1);
CorrTfit = polyval(LinearCoefft, hold_tot);

scatter_plot(mouseID, hold_dur_tot, react_time_tot, col_mat_s);
hold on
plot(hold_tot, CorrTfit, 'color', [0.5, 0.5,0.5], 'Linewidth', 1.5);
title('required hold time and lever react time');
xlabel('required hold time(ms)');
ylabel('react time(ms)');
saveas(fig, [out_base 'Summary_RT_holdtime_trial.fig']);
print([out_base 'Summary_RT_holdtime_trial.eps'], '-depsc');
print([out_base 'Summary_RT_holdtime_trial.pdf'], '-dpdf');

%% plot spike variability to cue and release
fig = figure;
% subplot(1,2,1)
scatter_plot(mouseID, peak2release_1, peak2cue_1, col_mat);
hold on
x = 0:1:250;
y = x;
plot(x,y,'k')
axis square
xlabel('std of spike relative to release')
ylabel('std of spike relative to cue')
title('cell variablitiy');

% subplot(1,2,2)
% scatter_plot(mouseID, peak2release_2, peak2cue_2, col_mat);
% hold on
% % x = 0:1:20;
% % y = x;
% plot(x,y,'k')
% axis square
% xlabel('std of spike relative to release')
% ylabel('std of spike relative to cue')
% title('event variability');

saveas(fig, [out_base 'Summary_spike_variablitiy.fig']);
print([out_base 'Summary_spike_variablitiy.eps'], '-depsc');
print([out_base 'Summary_spike_variablitiy.pdf'], '-dpdf');


%% plot TC for sustain and transient cells
success_TC_sus = success_TC_sus(~cellfun('isempty', success_TC_sus));
fail_TC_sus = fail_TC_sus(~cellfun('isempty', fail_TC_sus));

frame_size = cell2mat(cellfun(@size, success_TC_trans, 'UniformOutput', 0));
max_frame = max(frame_size(2:2:end));

success_TC_all_sus = interp_frame(success_TC_sus, max_frame);
fail_TC_all_sus    = interp_frame(fail_TC_sus, max_frame);

success_TC_all_trans = interp_frame(success_TC_trans, max_frame);
fail_TC_all_trans    = interp_frame(fail_TC_trans, max_frame);

fig = figure;
tt =(-max(cell2mat(pre_frames)):max(cell2mat(post_frames))).*double(min(cell2mat(TC_ifi)))./1000;
subplot(1,2,1)
shadedErrorBar(tt, mean(success_TC_all_trans,1), std(success_TC_all_trans,[],1)./sqrt(size(success_TC_all_trans,1)), 'k');
hold on;
shadedErrorBar(tt, mean(fail_TC_all_trans,1), std(fail_TC_all_trans,[],1)./sqrt(size(fail_TC_all_trans,1)), 'r');
ylim([-0.03 0.15])
xlim([-0.5,2.5])
xlabel('time (ms)')
ylabel('dF/F')
title(['transient cells n=',num2str(sum(RL_trans_cellNum))])

subplot(1,2,2)
shadedErrorBar(tt, mean(success_TC_all_sus,1), std(success_TC_all_sus,[],1)./sqrt(size(success_TC_all_sus,1)), 'k');
hold on;
shadedErrorBar(tt, mean(fail_TC_all_sus,1), std(fail_TC_all_sus,[],1)./sqrt(size(fail_TC_all_sus,1)), 'r');
ylim([-0.03 0.15])
xlim([-0.5,2.5])
xlabel('time (ms)')
ylabel('dF/F')
title(['sustain cells n=',num2str(sum(RL_sustain_cellNum))])


supertitle(['Timecourses of sustain and transient cells'])
saveas(fig, [out_base 'Summary_avgexpt_transient_sus_TCs.fig']);
print([out_base 'Summary_avgexpt_transient_sus_TCs.eps'], '-depsc');
print([out_base 'Summary_avgexpt_transient_sus_TCs.pdf'], '-dpdf');

% % scatter plot for early earlies and late earlies amplitude
% fig = figure;
% x = [-.05:.01:.3];
% y = x;
% subplot(1,2,1)
% scatter_plot(mouseID, early_fail_resp_RS_w1, late_fail_resp_RS_w1, col_mat_s);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% axis square
% xlabel('early early dF/F')
% ylabel('late early dF/F')
% title('first window')
% subplot(1,2,2)
% scatter_plot(mouseID, early_fail_resp_RS_w2, late_fail_resp_RS_w2, col_mat_s);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% axis square
% xlabel('early early dF/F')
% ylabel('late early dF/F')
% title('second window')
% saveas(fig, [out_base 'Summary_transient_sustain_resp_amp.fig']);
% print([out_base 'Summary_transient_sustain_resp_amp.eps'], '-depsc');
% print([out_base 'Summary_transient_sustain_resp_amp.pdf'], '-dpdf');