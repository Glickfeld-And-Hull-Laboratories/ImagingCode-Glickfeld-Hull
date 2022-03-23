% make a scatter plot to visualize amount of suppression during running before and after drug infusion
clear;
session_nodrug = '220126_img1911_1'; 
date = extractBetween(session_nodrug,1,14);
image_analysis_base_nodrug = 'Z:\home\shuyang\Analysis\motorizedWheel_Analysis\WF\imaging_analysis\'; 
image_analysis_dest_nodrug = [image_analysis_base_nodrug, session_nodrug, '\'];
mworks_nodrug = '1911-220126-1558_1'; 
behav_dest_nodrug = ['Z:\home\shuyang\Analysis\motorizedWheel_Analysis\WF\behavioral_analysis\' mworks_nodrug '\'];
dfOvF_nodrug = load([image_analysis_dest_nodrug session_nodrug '_dfOvF.mat']);
suppress_fast_ROIs_nodrug = dfOvF_nodrug.suppress_fastrun_ROIs;
suppress_slow_ROIs_nodrug = dfOvF_nodrug.suppress_slowrun_ROIs;

session_drug = '220126_img1911_1mMGabazine1uMAlexa1.5ul_1';
image_analysis_base_drug = 'Z:\home\shuyang\Analysis\motorizedWheel_Analysis\WF_drug\imaging_analysis\'; 
image_analysis_dest_drug = [image_analysis_base_drug, session_drug, '\'];
mworks_drug = '1911-220126-1629_1'; 
behav_dest_drug = ['Z:\home\shuyang\Analysis\motorizedWheel_Analysis\WF_drug\behavioral_analysis\' mworks_drug '\'];
dfOvF_drug = load([image_analysis_dest_drug session_drug '_dfOvF.mat']);
suppress_fast_ROIs_drug = dfOvF_drug.suppress_fastrun_ROIs;
suppress_slow_ROIs_drug = dfOvF_drug.suppress_slowrun_ROIs;
sessionnameL = strlength(session_drug);
drug_condition = extractBetween(session_drug,16,sessionnameL-2);

suppress_fast_ROIs_all = [suppress_fast_ROIs_nodrug,suppress_fast_ROIs_drug];
suppress_slow_ROIs_all = [suppress_slow_ROIs_nodrug,suppress_slow_ROIs_drug];

barplot = figure; 
subplot(1,2,1);
bar(suppress_fast_ROIs_all);
box off;
ylabel({'amount of suppression' 'during running in df/F'});
xlabel('individual ROIs');
ylim([-0.25 0.25]);
title(['fast running ',char(date)]);
%legend('before infusion',char(drug_condition)); legend boxoff;
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',12);
% fast_barplot.Units = 'centimeters';
% fast_barplot.Position = [1 3 6 4];
%savefig([image_analysis_dest_drug session_drug '_amountSuppression_fastRun.fig']);

subplot(1,2,2)
bar(suppress_slow_ROIs_all);
box off;
%ylabel('amount of suppression during running in df/F');
xlabel('individual ROIs');
title('slow running');
ylim([-0.25 0.25]);
set(gca,'XTickLabel',a,'FontSize',12);
legend('before infusion',char(drug_condition)); legend boxoff;
savefig([image_analysis_dest_drug session_drug '_amountSuppression.fig']);


%{
dfOvF_fig = figure;
x = [1,2,3];
dfOvF_plot = [dfOvF_fast_all,dfOvF_slow_all,dfOvF_stay_all];
dfOvF_ste_plot = [dfOvF_fast_ste,dfOvF_slow_ste,dfOvF_stay_ste];
% x_plot = repmat(x,size(dfOvF,1),1);
errorbar(x,dfOvF_plot,dfOvF_ste_plot,'.','LineStyle','none','MarkerSize',20);

xlim([0.5 3.5]);
x1= [1,2,3];
set(gca,'XTick',x1,'XTicklabel',{'fast running','slow running','stationary'});
ylabel('df/f');
title(sessions);
savefig([image_analysis_dest sessions '_dfOvF_scatter.fig']);
%}
