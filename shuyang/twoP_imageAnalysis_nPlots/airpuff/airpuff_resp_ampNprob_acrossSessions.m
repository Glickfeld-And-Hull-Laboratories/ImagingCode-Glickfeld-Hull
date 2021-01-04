clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = {'191114_img1040','191115_img1039','191115_img1041','191115_img1042'}; % ,'200316_img1064_airpuff_2'
image_analysis_base  = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
color_code = {'c','r','y','g'};

%% response probability of each cell across sessions
resp_prob_run_2std_all = [];
resp_prob_stay_2std_all = [];
resp_prob_run_3std_all = [];
resp_prob_stay_3std_all = [];

for ii = 1:4
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    dfOvF_strct = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    resp_prob_run_2std = dfOvF_strct.resp_prob_run_2std;
    resp_prob_stay_2std = dfOvF_strct.resp_prob_stay_2std;
    resp_prob_run_3std = dfOvF_strct.resp_prob_run_3std;
    resp_prob_stay_3std = dfOvF_strct.resp_prob_stay_3std;
    resp_prob_run_2std_all = cat(2,resp_prob_run_2std_all,resp_prob_run_2std);
    resp_prob_stay_2std_all = cat(2,resp_prob_stay_2std_all,resp_prob_stay_2std);
    resp_prob_run_3std_all = cat(2,resp_prob_run_3std_all,resp_prob_run_3std);
    resp_prob_stay_3std_all = cat(2,resp_prob_stay_3std_all,resp_prob_stay_3std);
end

mean_resp_prob_run_2std = mean(resp_prob_run_2std_all);
mean_resp_prob_stay_2std = mean(resp_prob_stay_2std_all);
mean_resp_prob_run_3std = mean(resp_prob_run_3std_all);
mean_resp_prob_stay_3std = mean(resp_prob_stay_3std_all);

% scatter_prob_2std = figure;
% scatter(resp_prob_stay_2std_all,resp_prob_run_2std_all,20,'filled','MarkerEdgeColor',...
%     [0.45 0.45 0.45],'MarkerFaceColor',[0.45 0.45 0.45]); hold on;
% scatter(mean_resp_prob_stay_2std,mean_resp_prob_run_2std,20,'filled','MarkerEdgeColor',...
%     [0.8431 0.0980 0.1098],'MarkerFaceColor',[0.8431 0.0980 0.1098]); hold on;
% axis square; xlabel('response probability (stationary)');ylabel('response probability (running)');
% xlim([0 1]);ylim([0 1]);
% line = refline(1,0);
% line.Color = 'r';
% line.LineWidth = 1;
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'FontSize',18);
% %scatter_prob_2std.Units = 'centimeters';
% %scatter_prob_2std.Position = [1 3 5 5];
% %fig_name = ['airpuff_resp_prob_scatter_' sessions{ii}];
% %path = 'Z:\Analysis\figures\figure5_airpuff_resp\';
% %print(scatter_prob,[path,fig_name],'-r600','-depsc');
% title('airpuff response probability, 2std across sessions');
% savefig(['Z:\Analysis\Airpuff_analysis\across_sessions\' 'airpuff_resp_prob_2std_scatter']);

scatter_prob_3std = figure;
scatter_prob_3std.Units = 'centimeters';
scatter_prob_3std.Position = [1 3 5 5];
scatter(resp_prob_stay_3std_all,resp_prob_run_3std_all,8,'filled','MarkerEdgeColor',...
    [0.45 0.45 0.45],'MarkerFaceColor',[0.45 0.45 0.45]); hold on;
scatter(mean_resp_prob_stay_3std,mean_resp_prob_run_3std,8,'filled','MarkerEdgeColor',...
    [0.8431 0.0980 0.1098],'MarkerFaceColor',[0.8431 0.0980 0.1098]); hold on;
axis square; xlabel('response probability (stationary)');ylabel('response probability (running)');
xlim([0 1]);ylim([0 1]);
line = refline(1,0);
line.Color = 'r';
line.LineWidth = 1;
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',8);
fig_name = 'airpuff_resp_prob_3std_scatter_across';
path = 'Z:\Analysis\figures\figure5_airpuff_resp\';
print(scatter_prob_3std,[path,fig_name],'-r600','-depsc');
% title('airpuff response probability, 3std across sessions');
% savefig(['Z:\Analysis\Airpuff_analysis\across_sessions\' 'airpuff_resp_prob_3std_scatter']);

%% response amplitude: pull out ONLY trials that has a response

nPC = 0;
resp_stay_ONLYh1_all = [];
resp_run_ONLYh1_all = [];
h_stay_all = [];
h_run_all = [];
for ii = 1:4
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    dfOvF_strct = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    resp_stay_ONLYh1 = dfOvF_strct.resp_stay_ONLYh1;
    resp_run_ONLYh1 = dfOvF_strct.resp_run_ONLYh1;
    h_stay = dfOvF_strct.h_stay;
    h_run = dfOvF_strct.h_run;
    resp_stay_ONLYh1_all = cat(2,resp_stay_ONLYh1_all,resp_stay_ONLYh1);
    resp_run_ONLYh1_all = cat(2,resp_run_ONLYh1_all,resp_run_ONLYh1);
    h_stay_all = cat(2,h_stay_all,h_stay);
    h_run_all = cat(2,h_run_all,h_run);
    nPC = nPC+length(h_stay);
end

%min(resp_stay_ONLYh1)
%min(resp_run_ONLYh1)
colorcode = {[0.8941 0.1020 0.1098],[0.2157 0.4941 0.7216],[0.3020 0.6863 0.2902],[0.4 0.4 0.4]};
amp_color_onlyresp1 = figure;hold on;
% scatter(100,100,'filled','MarkerEdgeColor',colorcode{1},'MarkerFaceColor',colorcode{1});
% scatter(100,100,'filled','MarkerEdgeColor',colorcode{2},'MarkerFaceColor',colorcode{2});
% scatter(100,100,'filled','MarkerEdgeColor',colorcode{3},'MarkerFaceColor',colorcode{3});
% scatter(100,100,'filled','MarkerEdgeColor',colorcode{4},'MarkerFaceColor',colorcode{4});
% having 2 for loops here because if there're more than 333 points in my
% scatter plot, it will become super blurry when importing into coreldraw, don't know why, fucking stupid. This has nothing to do with the colors or
% having too many points, somehow this is only a problem if there's if loops. So plotting 2 scatter plots seperatly and overlap them in coreldraw. 
for c = 1:333
    if h_stay_all(c) == 1&& h_run_all(c) == 0
        scatter(resp_stay_ONLYh1_all(c),resp_run_ONLYh1_all(c),8,'filled','MarkerEdgeColor',...
            colorcode{1},'MarkerFaceColor',colorcode{1});
    elseif h_stay_all(c) == 0 && h_run_all(c) == 1
        scatter(resp_stay_ONLYh1_all(c),resp_run_ONLYh1_all(c),8,'filled','MarkerEdgeColor',...
            colorcode{2},'MarkerFaceColor',colorcode{2});
    elseif h_run_all(c) == 1 && h_stay_all(c) == 1
        scatter(resp_stay_ONLYh1_all(c),resp_run_ONLYh1_all(c),8,'filled','MarkerEdgeColor',...
            colorcode{3},'MarkerFaceColor',colorcode{3});
    else%if h_run_all(c) == 0 && h_stay_all(c) == 0
        scatter(resp_stay_ONLYh1_all(c),resp_run_ONLYh1_all(c),8,'filled','MarkerEdgeColor',...
            colorcode{4},'MarkerFaceColor',colorcode{4});
    end
end
axis square; xlabel('resp. amp. in df/F (stationary)');ylabel('resp. amp. in df/F (running)');
xlim([0 2]); ylim([0 2]);
line = refline(1,0);
line.Color = [0 0 0];
line.LineWidth = 1;
%title(['response amplitude' sessions{ii}]);
amp_color_onlyresp1.Units = 'centimeters';
amp_color_onlyresp1.Position = [3 3 5 5];
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',8);
fig_name = 'resp_amp_colored_only_responsive_trials_across1';
%path = [image_analysis_dest 'resp\'];
path = 'Z:\Analysis\figures\figure5_airpuff_resp\';
%orient(amp_color_onlyh1,'landscape')
print(amp_color_onlyresp1,[path,fig_name],'-r600','-depsc');
hold off;
%figure;
amp_color_onlyresp2 = figure;hold on;
for c = 334:400
    if h_stay_all(c) == 1&& h_run_all(c) == 0
        scatter(resp_stay_ONLYh1_all(c),resp_run_ONLYh1_all(c),8,'filled','MarkerEdgeColor',...
            colorcode{1},'MarkerFaceColor',colorcode{1});
    elseif h_stay_all(c) == 0 && h_run_all(c) == 1
        scatter(resp_stay_ONLYh1_all(c),resp_run_ONLYh1_all(c),8,'filled','MarkerEdgeColor',...
            colorcode{2},'MarkerFaceColor',colorcode{2});
    elseif h_run_all(c) == 1 && h_stay_all(c) == 1
        scatter(resp_stay_ONLYh1_all(c),resp_run_ONLYh1_all(c),8,'filled','MarkerEdgeColor',...
            colorcode{3},'MarkerFaceColor',colorcode{3});
    elseif h_run_all(c) == 0 && h_stay_all(c) == 0
        scatter(resp_stay_ONLYh1_all(c),resp_run_ONLYh1_all(c),8,'filled','MarkerEdgeColor',...
            colorcode{4},'MarkerFaceColor',colorcode{4});
    end
end
%legend('responsive only stationary','responsive only running','responsive both','responsive neither'); legend('boxoff');
axis square; xlabel('resp. amp. in df/F (stationary)');ylabel('resp. amp. in df/F (running)');
xlim([0 2]); ylim([0 2]);
line = refline(1,0);
line.Color = [0 0 0];
line.LineWidth = 1;
%title(['response amplitude' sessions{ii}]);
amp_color_onlyresp2.Units = 'centimeters';
amp_color_onlyresp2.Position = [3 3 5 5];
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',8);
fig_name = 'resp_amp_colored_only_responsive_trials_across2';
%path = [image_analysis_dest 'resp\'];
path = 'Z:\Analysis\figures\figure5_airpuff_resp\';
%orient(amp_color_onlyh1,'landscape')
print(amp_color_onlyresp2,[path,fig_name],'-r600','-depsc');
hold off;
%savefig(['Z:\Analysis\Airpuff_analysis\across_sessions\' 'across_airpuff_resp_amp_onlyh1.fig']);

%% look at the percentage of cells that are responsive in either condition and both conditions
pcell_resp_stay_all = []; pcell_resp_run_all = [];
pcell_resp_both_all = [];

for ii = 1:4
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    dfOvF_strct = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    h_stay = dfOvF_strct.h_stay;
    h_run = dfOvF_strct.h_run;
    pcell_resp_stay = sum(h_stay==1)/length(h_stay);
    pcell_resp_run = sum(h_run==1)/length(h_run);
    pcell_resp_stay_all = cat(2,pcell_resp_stay_all,pcell_resp_stay);
    pcell_resp_run_all = cat(2,pcell_resp_run_all,pcell_resp_run);
    
    cell_resp_both = 0;
    for c = 1:length(h_stay)
        if h_stay(c) == 1 && h_run(c) == 1
            cell_resp_both = cell_resp_both + 1;
        end
    end
    pcell_resp_both = cell_resp_both/length(h_stay);
    pcell_resp_both_all = cat(2,pcell_resp_both_all,pcell_resp_both);

end

pcell_resp_stay_all_ave = mean(pcell_resp_stay_all);
pcell_resp_stay_all_ste = std(pcell_resp_stay_all)/length(pcell_resp_stay_all);

pcell_resp_run_all_ave = mean(pcell_resp_run_all);
pcell_resp_run_all_ste = std(pcell_resp_run_all)/length(pcell_resp_run_all);

pcell_resp_both_all_ave = mean(pcell_resp_both_all);
pcell_resp_both_all_ste = std(pcell_resp_both_all)/length(pcell_resp_both_all);


percent_resp = figure;hold on;
x1 = 1;
bar(x1,pcell_resp_stay_all_ave*100,'FaceColor',[0.8941 0.1020 0.1098]);
errorbar(x1,pcell_resp_stay_all_ave*100,pcell_resp_stay_all_ste*100,'.','Color','k','LineStyle','none','MarkerSize',4,'linewidth', 1);
x2 = 2;
bar(x2,pcell_resp_run_all_ave*100,'FaceColor',[0.2157 0.4941 0.7216]);
errorbar(x2,pcell_resp_run_all_ave*100,pcell_resp_run_all_ste*100,'.','Color','k','LineStyle','none','MarkerSize',4,'linewidth', 1);
x3 = 3;
bar(x3,pcell_resp_both_all_ave*100,'FaceColor',[0.3020 0.6863 0.2902]);
errorbar(x3,pcell_resp_both_all_ave*100,pcell_resp_both_all_ste*100,'.','Color','k','LineStyle','none','MarkerSize',4,'linewidth', 1);
x = [x1 x2 x3];
xlim([0.4 3.6]);
set(gca,'XTick',x,'XTicklabel',{'stationary','run','both'});
ylabel('% of cells responsive to airpuff');
xlabel('behavioral states');
ylim([0 100]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',8);
percent_resp.Units = 'centimeters';
percent_resp.Position = [3 3 5 5];
fig_name = 'across_session_resp_percent';
path = 'Z:\Analysis\figures\figure5_airpuff_resp\';
%fig_name = 'across_sessions_resp_perc';
%path = 'Z:\Analysis\Airpuff_analysis\across_sessions\';
print(percent_resp,[path,fig_name],'-r600','-dpdf');

