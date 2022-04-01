% plots for wide field imaging, self-paced experiments, acorss sessions
% staytionary vs. running across sessions
% df/F and speed aligned to running onset and offset across sessions
% df/f vs/ speed
% before running &. during running &. after running across sessions
% running window from 300ms before onset to 300ms after offset, each point is an average of 300ms

%% assign pathnames and datasets to be analyzed/written for moving dot experiments
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = {'190617_img1021_1','190617_img1023_1','190617_img1024_1',...
    '190618_img1025_1','200321_img1042_1','200321_img1049_1',...
    '200321_img1064_1'}; 
days = {'1021-190617_1','1023-190617_1','1024-190617_1',...
    '1025-190618_1','1042-200321_1','1049-200321_1',...
    '1064-200321_1'};
image_dest_base  = 'Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\'; %stores the data on crash in the movingDots analysis folder
% behavior analysis results 
color_code = {'b','r','k'};

%% df/f stay vs.run, line plot
ACS_dest = [image_dest_base 'acrossSessions\'];
ave_dfOvF_states_all = [];
for i = 1: length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    img_anal = load([image_dest '_imgAnalysis.mat' ]);
    ave_dfOvF_behav = img_anal.ave_dfOvF_behav;
    ave_dfOvF_states_all = cat(1,ave_dfOvF_states_all,ave_dfOvF_behav);
end
ave_dfOvF_states_allsessions = mean(ave_dfOvF_states_all);
ste_dfOvF_states_allsessions = std(ave_dfOvF_states_all)/sqrt(size(ave_dfOvF_states_all,1));

x = [1,2]; x_plot = repmat(x,size(ave_dfOvF_states_all,1),1);
dfOvF_behavStates = figure;
plot(x_plot',ave_dfOvF_states_all','.','LineStyle','-','linewidth', 0.8,'MarkerSize',6,'color',[0.9843 0.7059 0.6824]);hold on;
errorbar(x,ave_dfOvF_states_allsessions,ste_dfOvF_states_allsessions,'.','LineStyle','-',...
    'linewidth', 1,'Color',[0 0 0],'MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
xlabel ('behavioral state');xlim([0.5 2.5]); axis square;
set(gca,'XTick',x,'XTicklabel',{'stationary','run'});
ylabel('df/F');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
%title('df/f for each beahvioral state all sessions'); %legend('ROI1','ROI2','ROI3','ROI4','ROI5');

dfOvF_behavStates.Units = 'centimeters';
dfOvF_behavStates.Position = [1 1 5 5];
fig_name = 'across_session_stayVsrun_';
path = 'Z:\Analysis\figures\figure1_WF\';
%orient(dfOvF_behavStates,'landscape');
print(dfOvF_behavStates,[path,fig_name],'-r600','-depsc');
    
%saveas(dfOvF_behavStates, [ACS_dest 'dfOvF_behavStates_20200416']);
save([ACS_dest 'ACSanalysis_20200527.mat'],'ave_dfOvF_states_all',...
    'ave_dfOvF_states_allsessions','ste_dfOvF_states_allsessions');


%% run trig ave across sessions
ave_dfOvF_runTrigger_all = []; ave_speed_runTrigger_all = []; 
ACS_dest = [image_dest_base 'acrossSessions\'];
for i = 1: length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    img_anal = load([image_dest '_imgAnalysis.mat' ]);
    ave_dfOvF_runTrigger = img_anal.ave_dfOvF_runTrigger;
    if size(ave_dfOvF_runTrigger,2) == 1
        ave_dfOvF_runTrigger = ave_dfOvF_runTrigger';
    end
    ave_dfOvF_runTrigger_all = cat(1,ave_dfOvF_runTrigger_all,ave_dfOvF_runTrigger);
    ave_speed_runTrigger = img_anal.ave_speed_runTrigger;
    ave_speed_runTrigger_all = cat(1,ave_speed_runTrigger_all,ave_speed_runTrigger );
end
%acs: across sessions
ave_dfOvF_runTrigger_ACS = mean(ave_dfOvF_runTrigger_all);
ste_dfOvF_runTrigger_ACS = std(ave_dfOvF_runTrigger_all)/sqrt(length(ave_dfOvF_runTrigger_all));
ave_speed_runTrigger_ACS = mean(ave_speed_runTrigger_all);
ste_speed_runTrigger_ACS = std(ave_speed_runTrigger_all)/sqrt(length(ave_speed_runTrigger_all));

x = ((0:length(ave_dfOvF_runTrigger_ACS)-1)/10)-1;
dfOvF_runTrigger_ACS = figure;
subplot(2,1,1);hold on;
shadedErrorBar(x,ave_dfOvF_runTrigger_ACS,ste_dfOvF_runTrigger_ACS,{'color',[0.8431    0.0980    0.1098]},{'Linewidth',1}); 
% for errorbar: ,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
xlim([-1.1 1.6]);
ylim([0 0.35]);
vline(0,'k');
ylabel('df/F'); 
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7); %with the size of [1 1 5.5 5],font 8 won't plot the whole time scale 

subplot(2,1,2);hold on;
shadedErrorBar(x,ave_speed_runTrigger_ACS*2*3.1415926*7.5/128,ste_speed_runTrigger_ACS*2*3.1415926*7.5/128,...
    {'color',[0 0 0]},{'Linewidth',1});
%'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
xlabel('time from running onset(s)');
ylabel('speed(cm/s)');
ylim([0 15]);
vline(0,'k');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
xlim([-1.1 1.6]);
dfOvF_runTrigger_ACS.Units = 'centimeters';
dfOvF_runTrigger_ACS.Position = [1 1 5.5 5];
fig_name = 'across_session_runonset';
path = 'Z:\Analysis\figures\figure1_WF\';
%orient(dfOvF_behavStates,'landscape');
print(dfOvF_runTrigger_ACS,[path,fig_name],'-r600','-depsc');

%supertitle('run triggered average across sessions'); 
%saveas(dfOvF_runTrigger_ACS, [ ACS_dest 'dfOvF_runTrigger_ACS_20200416']);
save([ ACS_dest 'ACSanalysis_20200527.mat' ],'ave_dfOvF_runTrigger_ACS','ste_dfOvF_runTrigger_ACS',...
   'ave_speed_runTrigger_ACS','ste_speed_runTrigger_ACS', '-append' );


%% run offset ave across sessions
ave_dfOvF_runOff_all = []; ave_speed_runOff_all = []; 
ACS_dest = [image_dest_base 'acrossSessions\'];
for i = 1: length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    img_anal = load([image_dest '_imgAnalysis.mat' ]);
    ave_dfOvF_runOff = img_anal.ave_dfOvF_runOff;
    if size(ave_dfOvF_runOff,2) == 1
        ave_dfOvF_runOff = ave_dfOvF_runOff';
    end
    ave_dfOvF_runOff_all = cat(1,ave_dfOvF_runOff_all,ave_dfOvF_runOff);
    ave_speed_runOff = img_anal.ave_speed_runoff;
    ave_speed_runOff_all = cat(1,ave_speed_runOff_all,ave_speed_runOff );
end
%acs: across sessions
ave_dfOvF_runOff_ACS = mean(ave_dfOvF_runOff_all);
ste_dfOvF_runOff_ACS = std(ave_dfOvF_runOff_all)/sqrt(length(ave_dfOvF_runOff_all));
ave_speed_runOff_ACS = mean(ave_speed_runOff_all);
ste_speed_runOff_ACS = std(ave_speed_runOff_all)/sqrt(length(ave_speed_runOff_all));

x = ((0:length(ave_dfOvF_runOff_ACS)-1)/10)-1;
dfOvF_runOff_ACS = figure;
subplot(2,1,1);hold on;
shadedErrorBar(x,ave_dfOvF_runOff_ACS,ste_dfOvF_runOff_ACS,{'color',[0.8431    0.0980    0.1098]},{'Linewidth',1});
%'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
xlim([-1.1 1.6]);
ylim([0 0.35]);
vline(0, 'k');
ylabel('df/f'); 
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);

subplot(2,1,2);hold on;
shadedErrorBar(x,ave_speed_runOff_ACS*2*3.1415926*7.5/128,ste_speed_runOff_ACS*2*3.1415926*7.5/128,{'color',[0 0 0]},{'Linewidth',1});
%'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
xlabel('time from running offset(s)');
ylabel('speed(cm/s)');
ylim([0 15]);
vline(0, 'k');
xlim([-1.1 1.6]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);

dfOvF_runOff_ACS.Units = 'centimeters';
dfOvF_runOff_ACS.Position = [1 1 5.5 5];
fig_name = 'across_session_runoffset';
path = 'Z:\Analysis\figures\figure1_WF\';
%orient(dfOvF_behavStates,'landscape');
print(dfOvF_runOff_ACS,[path,fig_name],'-r600','-depsc');

%supertitle('run offset average across sessions'); 
%saveas(dfOvF_runOff_ACS, [ ACS_dest 'dfOvF_runOff_ACS_20200416']);
save([ ACS_dest 'ACSanalysis_20200527.mat' ],'ave_dfOvF_runOff_ACS','ste_dfOvF_runOff_ACS',...
   'ave_speed_runOff_ACS','ste_speed_runOff_ACS', '-append' );


%% SECTION I - draw df/f vs. speed across sessions 
ACS_dest = [image_dest_base 'acrossSessions\'];
% the coding for this section is super not straight forward because
% there're different number of RIOs and different number of unique speeds
% in each session, so have to used cells and NaNs and can't cell2mat. At
% the end, put all speeds and df/F into a vector. In which each df/f value
% of an ROI is correlated with the speed value. This is all what the first
% three hard-to-understand for loops are doing.
% e.g.: you get isp_all_mat, and it's -20,-20,-20,-10,-10,-10,0,0,0...
% because there's 3 ROIs in the first session. and then -20,-20,-10,-10...
% 2 ROIs in the second session. 
dfOvF_spd_all = {}; isp_all = {};
for i = 1: length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    img_anal = load([image_dest '_imgAnalysis.mat' ]);
    isp = img_anal.isp; isp = double(isp);
    dfOvF_spdmean = img_anal.dfOvF_spdmean; % ROI*number of unique speeds in each session
    dfOvF_spd_all = cat(2,dfOvF_spd_all, dfOvF_spdmean); %different sessions have different number of ROIs and number of unique speeds, so need to use cell
    isp_all = cat(2,isp_all, isp); % a cell variable, each element is one session
end
%generate matrix for speeds and df/F
isp_temp = cell2mat(cellfun(@size,isp_all, 'UniformOutput',0));
isp_length = max(isp_temp);
isp_mat = nan(size(isp_all,2), isp_length);
for n = 1: size(isp_mat,1) % for each session
    temp_speed = isp_all{n};
    isp_mat(n,1:size(temp_speed,2)) = temp_speed; % session * number of unique speeds with NaNs
end

dfOvF_spd_all_mat = []; isp_all_mat = [];
for f = 1: size(dfOvF_spd_all,2) % for each session
    temp_dfOvF = dfOvF_spd_all{f}; 
    dfOvF_spd_all_mat = [dfOvF_spd_all_mat; temp_dfOvF(:)]; % dfOvF values related to unique speeds in each session
    %make the matrix of speed exactly the same as df/f
    [isp_mat1,isp_mat2] = size(temp_dfOvF);
    temp_isp = isp_mat(f,:);
    temp_isp(isnan(temp_isp)) = []; %delate NaN
    temp_isp2 = repmat(temp_isp,isp_mat1,1); %repeat because there're more than 1 ROI per session
    isp_all_mat = [isp_all_mat; temp_isp2(:)]; % all uniques speeds of all sessions, speed values can be repeated if multiple session has this value
end

spd_all_plotx = sort(unique(isp_all_mat));
dfOvF_all_ploty = zeros(length(spd_all_plotx),1);
dfOvF_errbar = zeros(length(spd_all_plotx),1);
for j = 1: length(spd_all_plotx)
    c = find(isp_all_mat == spd_all_plotx(j));
    dfOvF_all_ploty(j) = mean(dfOvF_spd_all_mat(c)); % average across all ROIs from all sessions
    dfOvF_errbar(j) = std(dfOvF_spd_all_mat(c))/sqrt(length(c));
end

dfOvF_vs_spd_ACS = figure; 
shadedErrorBar(spd_all_plotx*2*3.1415926*7.5/128,dfOvF_all_ploty,dfOvF_errbar,{'color',[0.1922 0.6392 0.3294]});
%,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
xlabel ('speed(cm/s)');
ylabel('ave df/f');
xlim([0 35]);
title('df/f vs. speed across sessions'); 
saveas(dfOvF_vs_spd_ACS, [ ACS_dest 'dfOvF_vs_spd_ACS_20200527']);
save([ ACS_dest 'ACSanalysis_20200527.mat' ],'spd_all_plotx','dfOvF_all_ploty','dfOvF_errbar','-append');

% bin speeds and draw scatter plot---------------------------------------------------
% speed category (in quadrature units): <0, 0, 0-30,30-60,60-90
ncategory = 6;
dfOvF_catSpeed = zeros(1,ncategory);
stedfOvF_catSpeed = zeros(1,ncategory);
across_anal = load([ACS_dest 'ACSanalysis_20200527.mat']);
ave_dfOvF_states_allsessions = across_anal.ave_dfOvF_states_allsessions;
ste_dfOvF_states_allsessions = across_anal.ste_dfOvF_states_allsessions;
s1 = 0;
isp_all_mat_cm = isp_all_mat*2*3.1415926*7.5/128;
for c = 1:ncategory
    if c == 1 %negative speed
        idx = find(isp_all_mat_cm < 0);
        dfOvF_catSpeed(c) =  mean(dfOvF_spd_all_mat(idx)); % df/f related to this speed category
        stedfOvF_catSpeed(c) = std(dfOvF_spd_all_mat(idx))/sqrt(length(idx));
    elseif c == 2 %use ave df/f during stationary for speed = 0
        dfOvF_catSpeed(c) = ave_dfOvF_states_allsessions(1);
        stedfOvF_catSpeed(c) = ste_dfOvF_states_allsessions(1);
    else
        s1 = s1+10;
        idx = find(isp_all_mat_cm > s1-10 & isp_all_mat_cm <= s1);
        dfOvF_catSpeed(c) =  mean(dfOvF_spd_all_mat(idx));
        stedfOvF_catSpeed(c) = std(dfOvF_spd_all_mat(idx))/sqrt(length(idx));
    end
end

x = [1,2,3,4,5,6];
dfOvF_speed_category = figure;
errorbar(x,dfOvF_catSpeed,stedfOvF_catSpeed,'.','Color',[0.8431 0.0980 0.1098],...
    'MarkerSize',8,'MarkerEdgeColor',[0.8431 0.0980 0.1098],'MarkerFaceColor',[0.8431 0.0980 0.1098]);
xlabel ('speed(cm/s)');xlim([0.5 6.5]); axis square;
set(gca,'XTick',x,'XTicklabel',{'< 0','0','0-10','10-20','20-30','> 30'});
ylabel('df/F');ylim([0 0.3]);xlim([0.5 5.5]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
set(gca,'XTickLabelRotation',45);

dfOvF_speed_category.Units = 'centimeters';
dfOvF_speed_category.Position = [1 1 5 5];
fig_name = 'across_session_dfOvFvsSpd';
path = 'Z:\Analysis\figures\figure1_WF\';
%orient(dfOvF_behavStates,'landscape');
print(dfOvF_speed_category,[path,fig_name],'-r600','-depsc');

%title('df/f vs. speed all sessions'); 
%saveas(dfOvF_speed_category, [ACS_dest 'dfOvF_vs_spdCategory_20200416']);
save([ACS_dest 'ACSanalysis_20200527.mat'],'dfOvF_catSpeed','stedfOvF_catSpeed','-append');


%% SECTION III - df/f for beforeRunAfter across sessions 
ave_dfOvF_befoRunAft_all = [];
ACS_dest = [image_dest_base 'acrossSessions\'];
for i = 1: length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    img_anal = load([image_dest '_imgAnalysis.mat' ]);
    ave_befoRunaft = img_anal.ave_befoRunaft;
    ave_dfOvF_befoRunAft_all = cat(1,ave_dfOvF_befoRunAft_all,ave_befoRunaft ); 
end
%acs: across sessions
ave_dfOvF_befoRunAft_ACS = mean(ave_dfOvF_befoRunAft_all);
ste_dfOvF_befoRunAft_ACS = std(ave_dfOvF_befoRunAft_all)/sqrt(length(ave_dfOvF_befoRunAft_all));

x = [1,2,3];
dfOvF_befoRunAft_ACS = figure;
errorbar(x,ave_dfOvF_befoRunAft_ACS,ste_dfOvF_befoRunAft_ACS,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
plot(x,ave_dfOvF_befoRunAft_ACS,'.','Color','b','MarkerSize',20,'LineStyle','none');

xlabel ('behavioral state');xlim([0.5 3.5]);
set(gca,'XTick',x,'XTicklabel',{'before','run','after'});
ylabel('df/f');
title('df/f around running across sessions'); 
saveas(dfOvF_befoRunAft_ACS, [ ACS_dest 'dfOvF_befoRunAft_ACS_20200416']);
save([ ACS_dest 'ACSanalysis_20200416.mat' ],'ave_dfOvF_befoRunAft_ACS','ste_dfOvF_befoRunAft_ACS','-append' );



%{
%% SECTION IV - df/f and speed 300ms across sessions
dfOvF_bl300ms_all = []; speed_bl300ms_all = [];
ACS_dest = [image_dest_base 'acrossSessions\'];
period = 9;
for i = 1: length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    img_anal = load([image_dest '_imgAnalysis.mat' ]);
    speed_plot300ms = img_anal.speed_plot300ms;
    dfOvF_plot300ms = img_anal.dfOvF_plot300ms;
    dfOvF_bgPart300ms = dfOvF_plot300ms(:,(1: 1+period)); %first point is before running
    dfOvF_lstPart300ms = dfOvF_plot300ms(:,(end-period:end));%last point is after running
    dfOvF_blPart300ms = [dfOvF_bgPart300ms,dfOvF_lstPart300ms];%bl: begin and last
    speed_bgPart300ms = speed_plot300ms(:,(1: 1+period));
    speed_lstPart300ms = speed_plot300ms(:,(end-period:end));
    speed_blPart300ms = [speed_bgPart300ms,speed_lstPart300ms];%bl: begin and last
    
    dfOvF_bl300ms_all = cat(1,dfOvF_bl300ms_all,dfOvF_blPart300ms );
    speed_bl300ms_all = cat(1,speed_bl300ms_all,speed_blPart300ms );
    
end

ave_dfOvF_bl300ms = mean(dfOvF_bl300ms_all);
ave_speed_bl300ms = mean(speed_bl300ms_all);
ste_dfOvF_bl300ms = std(dfOvF_bl300ms_all)/sqrt(length(dfOvF_bl300ms_all));
ste_speed_bl300ms = std(speed_bl300ms_all)/sqrt(length(speed_bl300ms_all));

x = (1: 1: 20);
dfOvF_run300ms_ACS = figure;
subplot(2,1,1);hold on;
errorbar(x,ave_dfOvF_bl300ms,ste_dfOvF_bl300ms,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
%xlim([-5 10]);
ylim([-0.1 0]);
ylabel('df/f'); 

subplot(2,1,2);hold on;
errorbar(x,ave_speed_bl300ms,ste_speed_bl300ms,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
xlabel('frames');
ylabel('speed');
%xlim([-5 10]);

supertitle(['run 300ms average across sessions']); 
saveas(dfOvF_run300ms_ACS, [ ACS_dest 'dfOvF_run300ms_ACS']);
save([ ACS_dest 'ACSanalysis.mat' ],'ave_dfOvF_bl300ms','ave_speed_bl300ms',...
   'ste_dfOvF_bl300ms','ste_speed_bl300ms', '-append' );
%}

