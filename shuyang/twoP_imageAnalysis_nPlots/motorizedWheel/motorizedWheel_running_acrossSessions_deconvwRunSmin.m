%% assign document paths and experimental sessions
clear;
sessions = {'200116_img1041','200217_img1061','200225_img1049',...
    '200319_img1064','200319_img1064_2'}; %'200116_img1041'
image_analysis_base = 'Z:\Analysis\motorizedWheel_Analysis\running\imaging_analysis\';
% behavior analysis results
days = {'1041-200116_1','1061-200217_1','1049-200225_1',...
    '1064-200319_1','1064-200319_2'};%
color_code = {'b','r','k','c'};

%% df/f runtrig ave and runoffset ave and speed transit (df/f parts are the same as when doing deconvolution together just using 4std of the whole baseline as threshold)
dfOvF_runoff_fast_cells_allsession = []; %frame*cell
dfOvF_runoff_slow_cells_allsession = [];
dfOvF_runtrig_fast_cells_allsession = [];
dfOvF_runtrig_slow_cells_allsession = [];
dfOvF_spd_decrease_cells_allsession = [];
dfOvF_spd_increase_cells_allsession = [];

nPCs_total = 0;

for ii = 1:length(sessions)
    %load things
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    dfOvF_output = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    dfOvF_runoff_fast_cells = dfOvF_output.dfOvF_runoff_fast_cells; %frames*cells
    dfOvF_runoff_slow_cells = dfOvF_output.dfOvF_runoff_slow_cells; %frames*cells
    dfOvF_runtrig_fast_cells = dfOvF_output.dfOvF_runtrig_fast_cells;
    dfOvF_runtrig_slow_cells = dfOvF_output.dfOvF_runtrig_slow_cells;
    dfOvF_spd_decrease_cells = dfOvF_output.dfOvF_spd_decrease_cells;
    dfOvF_spd_increase_cells = dfOvF_output.dfOvF_spd_increase_cells;
    
    nPCs_total = nPCs_total + size(dfOvF_runoff_fast_cells,2);
    
    dfOvF_runoff_fast_cells_allsession = cat(2,dfOvF_runoff_fast_cells_allsession,dfOvF_runoff_fast_cells);%when cat matrix, the matrix with NaN will be cat. And you have matrix with NaN because for some sessions there's no running trials fullfill the runtrig/runoff criteria
    dfOvF_runoff_slow_cells_allsession = cat(2,dfOvF_runoff_slow_cells_allsession,dfOvF_runoff_slow_cells);
    dfOvF_runtrig_fast_cells_allsession = cat(2,dfOvF_runtrig_fast_cells_allsession,dfOvF_runtrig_fast_cells);
    dfOvF_runtrig_slow_cells_allsession = cat(2,dfOvF_runtrig_slow_cells_allsession,dfOvF_runtrig_slow_cells);
    dfOvF_spd_decrease_cells_allsession = cat(2,dfOvF_spd_decrease_cells_allsession,dfOvF_spd_decrease_cells);
    dfOvF_spd_increase_cells_allsession = cat(2,dfOvF_spd_increase_cells_allsession,dfOvF_spd_increase_cells); 
end
dfOvF_runoff_fast_cells_acrosssession = mean(dfOvF_runoff_fast_cells_allsession,2);%when cat matrix, the matrix with NaN will be cat. And you have matrix with NaN because for some sessions there's no running trials fullfill the runtrig/runoff criteria
ste_dfOvF_runoff_fast_cells = std(dfOvF_runoff_fast_cells_allsession,0,2)/sqrt(size(dfOvF_runoff_fast_cells_allsession,2));
dfOvF_runoff_slow_cells_acrossssession = mean(dfOvF_runoff_slow_cells_allsession,2);
ste_dfOvF_runoff_slow_cells = std(dfOvF_runoff_slow_cells_allsession,0,2)/sqrt(size(dfOvF_runoff_slow_cells_allsession,2));
dfOvF_runtrig_fast_cells_acrossssession = mean(dfOvF_runtrig_fast_cells_allsession,2);
ste_dfOvF_runtrig_fast_cells = std(dfOvF_runtrig_fast_cells_allsession,0,2)/sqrt(size(dfOvF_runtrig_fast_cells_allsession,2));
dfOvF_runtrig_slow_cells_acrossssession = mean(dfOvF_runtrig_slow_cells_allsession,2);
ste_dfOvF_runtrig_slow_cells = std(dfOvF_runtrig_slow_cells_allsession,0,2)/sqrt(size(dfOvF_runtrig_slow_cells_allsession,2));
dfOvF_spd_decrease_cells_acrossssession = mean(dfOvF_spd_decrease_cells_allsession,2);
ste_dfOvF_spd_decrease_cells = std(dfOvF_spd_decrease_cells_allsession,0,2)/sqrt(size(dfOvF_spd_decrease_cells_allsession,2));
dfOvF_spd_increase_cells_acrossssession = mean(dfOvF_spd_increase_cells_allsession,2);
ste_dfOvF_spd_increase_cells = std(dfOvF_spd_increase_cells_allsession,0,2)/sqrt(size(dfOvF_spd_increase_cells_allsession,2));
%%
x = (1:size(dfOvF_spd_increase_cells,1))/30 -1;
dfOvF_across = figure;
dfOvF_across.Units = 'centimeters';
dfOvF_across.Position = [3 3 11 12];
subplot(3,2,1);
shadedErrorBar(x,dfOvF_runtrig_slow_cells_acrossssession,ste_dfOvF_runtrig_slow_cells,{'color',[0.8431 0.0980 0.1098]});
ylabel('df/f'); ylim([0 0.85]); xlim([-1 1.5]);
vline(0,'k');
title('stay --> slow speed');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
subplot(3,2,2);
shadedErrorBar(x,dfOvF_runoff_slow_cells_acrossssession,ste_dfOvF_runoff_slow_cells,{'color',[0.8431 0.0980 0.1098]});
%fast_errbar(x,dfOvF_runoff_slow_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylim([0 0.85]);xlim([-1 1.5]);
vline(0,'k');
title('slow speed --> stay');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
subplot(3,2,3);
shadedErrorBar(x,dfOvF_runtrig_fast_cells_acrossssession,ste_dfOvF_runtrig_fast_cells,{'color',[0.8431 0.0980 0.1098]});
%fast_errbar(x,dfOvF_runtrig_fast_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylim([0 0.85]); xlim([-1 1.5]); vline(0,'k'); 
ylabel('df/f');
%xlabel('time from behavioral state change(s)');xlim([0 2.5]);
title('stay --> fast speed');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
subplot(3,2,4);
shadedErrorBar(x,dfOvF_runoff_fast_cells_acrosssession,ste_dfOvF_runoff_fast_cells,{'color',[0.8431 0.0980 0.1098]});
%fast_errbar(x,dfOvF_runoff_fast_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylim([0 0.85]);xlim([-1 1.5]); vline(0,'k'); 
title('fast speed --> stay');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
subplot(3,2,5);
shadedErrorBar(x,dfOvF_spd_decrease_cells_acrossssession,ste_dfOvF_spd_decrease_cells,{'color',[0.8431 0.0980 0.1098]})
%fast_errbar(x,dfOvF_spd_decrease_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylim([0 0.85]); xlim([-1 1.5]);vline(0,'k');
ylabel('df/f');
title('speed decrease');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
subplot(3,2,6);
shadedErrorBar(x,dfOvF_spd_increase_cells_acrossssession,ste_dfOvF_spd_increase_cells,{'color',[0.8431 0.0980 0.1098]});
%fast_errbar(x,dfOvF_spd_increase_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
xlim([-1 1.5]); ylim([0 0.85]); vline(0,'k');
xlabel('time from behavioral state change(s)');
title('speed increase');
%supertitle(['dfOvF across sessions nPCs = ' num2str(nPCs_total)]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);

fig_name = 'across_session_motorized_running_dfOvF';
path = 'Z:\Analysis\figures\figure6_motorized_running\';
print(dfOvF_across,[path,fig_name],'-r600','-depsc');

%savefig(['Z:\Analysis\motorizedWheel_Analysis\running\across_sessions\' 'across_sessions_dfOvF.fig']);

%% FR run trig ave and runoff ave (this part is the same as before because we're not putting this in the paper)

FR_runoff_fast_cells_allsession = []; %frame*cell
FR_runoff_slow_cells_allsession = [];
FR_runtrig_fast_cells_allsession = [];
FR_runtrig_slow_cells_allsession = [];
FR_spd_decrease_cells_allsession = [];
FR_spd_increase_cells_allsession = [];

nPCs_total = 0;

for ii = 1:length(sessions)
    %load things
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    spk_output = load([image_analysis_dest sessions{ii},'_spk_deconvolve_threshold-4.mat']);
    FR_runoff_fast_cells = spk_output.FR_runoff_fast_cells; %frames*cells
    FR_runoff_slow_cells = spk_output.FR_runoff_slow_cells; %frames*cells
    FR_runtrig_fast_cells = spk_output.FR_runtrig_fast_cells;
    FR_runtrig_slow_cells = spk_output.FR_runtrig_slow_cells;
    FR_spd_decrease_cells = spk_output.FR_spd_decrease_cells;
    FR_spd_increase_cells = spk_output.FR_spd_increase_cells;
    
    nPCs_total = nPCs_total + size(FR_runoff_fast_cells,2);
    
    FR_runoff_fast_cells_allsession = cat(2,FR_runoff_fast_cells_allsession,FR_runoff_fast_cells);%when cat matrix, the matrix with NaN will be cat. And you have matrix with NaN because for some sessions there's no running trials fullfill the runtrig/runoff criteria
    FR_runoff_slow_cells_allsession = cat(2,FR_runoff_slow_cells_allsession,FR_runoff_slow_cells);
    FR_runtrig_fast_cells_allsession = cat(2,FR_runtrig_fast_cells_allsession,FR_runtrig_fast_cells);
    FR_runtrig_slow_cells_allsession = cat(2,FR_runtrig_slow_cells_allsession,FR_runtrig_slow_cells);
    FR_spd_decrease_cells_allsession = cat(2,FR_spd_decrease_cells_allsession,FR_spd_decrease_cells);
    FR_spd_increase_cells_allsession = cat(2,FR_spd_increase_cells_allsession,FR_spd_increase_cells); 
end

FR_runoff_fast_cells_acrosssession = mean(FR_runoff_fast_cells_allsession,2);%when cat matrix, the matrix with NaN will be cat. And you have matrix with NaN because for some sessions there's no running trials fullfill the runtrig/runoff criteria
ste_FR_runoff_fast_cells = std(FR_runoff_fast_cells_allsession,0,2)/sqrt(size(FR_runoff_fast_cells_allsession,2));
FR_runoff_slow_cells_acrossssession = mean(FR_runoff_slow_cells_allsession,2);
ste_FR_runoff_slow_cells = std(FR_runoff_slow_cells_allsession,0,2)/sqrt(size(FR_runoff_slow_cells_allsession,2));
FR_runtrig_fast_cells_acrossssession = mean(FR_runtrig_fast_cells_allsession,2);
ste_FR_runtrig_fast_cells = std(FR_runtrig_fast_cells_allsession,0,2)/sqrt(size(FR_runtrig_fast_cells_allsession,2));
FR_runtrig_slow_cells_acrossssession = mean(FR_runtrig_slow_cells_allsession,2);
ste_FR_runtrig_slow_cells = std(FR_runtrig_slow_cells_allsession,0,2)/sqrt(size(FR_runtrig_slow_cells_allsession,2));
FR_spd_decrease_cells_acrossssession = mean(FR_spd_decrease_cells_allsession,2);
ste_FR_spd_decrease_cells = std(FR_spd_decrease_cells_allsession,0,2)/sqrt(size(FR_spd_decrease_cells_allsession,2));
FR_spd_increase_cells_acrossssession = mean(FR_spd_increase_cells_allsession,2);
ste_FR_spd_increase_cells = std(FR_spd_increase_cells_allsession,0,2)/sqrt(size(FR_spd_increase_cells_allsession,2));

%%
x = (1:size(FR_spd_increase_cells,1))/30 -1;
FR_across = figure;
FR_across.Units = 'centimeters';
FR_across.Position = [3 3 11 12];
subplot(3,2,1);
shadedErrorBar(x,FR_runtrig_slow_cells_acrossssession,ste_FR_runtrig_slow_cells,{'color',[0.1922 0.6392 0.3294]});
ylabel('firing rate(Hz)'); ylim([0 5.5]); xlim([-1 1.5]);
vline(0,'k');
title('stay --> slow speed');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
subplot(3,2,2);
shadedErrorBar(x,FR_runoff_slow_cells_acrossssession,ste_FR_runoff_slow_cells,{'color',[0.1922 0.6392 0.3294]});
%fast_errbar(x,dfOvF_runoff_slow_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylim([0 5.5]);xlim([-1 1.5]);
vline(0,'k');
title('slow speed --> stay');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
subplot(3,2,3);
shadedErrorBar(x,FR_runtrig_fast_cells_acrossssession,ste_FR_runtrig_fast_cells,{'color',[0.1922 0.6392 0.3294]});
%fast_errbar(x,dfOvF_runtrig_fast_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylim([0 5.5]); xlim([-1 1.5]); vline(0,'k'); 
ylabel('firing rate(Hz)');
%xlabel('time from behavioral state change(s)');xlim([0 2.5]);
title('stay --> fast speed');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
subplot(3,2,4);
shadedErrorBar(x,FR_runoff_fast_cells_acrosssession,ste_FR_runoff_fast_cells,{'color',[0.1922 0.6392 0.3294]});
%fast_errbar(x,dfOvF_runoff_fast_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylim([0 5.5]);xlim([-1 1.5]); vline(0,'k'); 
title('fast speed --> stay');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
subplot(3,2,5);
shadedErrorBar(x,FR_spd_decrease_cells_acrossssession,ste_FR_spd_decrease_cells,{'color',[0.1922 0.6392 0.3294]})
%fast_errbar(x,dfOvF_spd_decrease_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylim([0 5.5]); xlim([-1 1.5]);vline(0,'k');
ylabel('firing rate(Hz)');
title('speed decrease');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
subplot(3,2,6);
shadedErrorBar(x,FR_spd_increase_cells_acrossssession,ste_FR_spd_increase_cells,{'color',[0.1922 0.6392 0.3294]});
%fast_errbar(x,dfOvF_spd_increase_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
xlim([-1 1.5]); ylim([0 5.5]); vline(0,'k');
xlabel('time from behavioral state change(s)');
title('speed increase');
%supertitle(['dfOvF across sessions nPCs = ' num2str(nPCs_total)]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);

fig_name = 'across_session_motorized_running_FR';
path = 'Z:\Analysis\figures\figure6_motorized_running\';
%fig_name = 'across_sessions_resp_perc';
%path = 'Z:\Analysis\Airpuff_analysis\across_sessions\';
print(FR_across,[path,fig_name],'-r600','-depsc');

% savefig(['Z:\Analysis\motorizedWheel_Analysis\running\across_sessions\' 'across_sessions_FR.fig']);

%% dfOvF scatter plot

dfOvF_fast_cells_allsession = []; %frame*cell
dfOvF_slow_cells_allsession = [];
dfOvF_stay_cells_allsession = [];
nPCs_total = 0;

for ii = 1:length(sessions)
    %load things
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    dfOvF_output = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    dfOvF_fast_cells = dfOvF_output.dfOvF_fast_cells; %1*n cells
    dfOvF_slow_cells = dfOvF_output.dfOvF_slow_cells;
    dfOvF_stay_cells = dfOvF_output.dfOvF_stay_cells;
    
    nPCs_total = nPCs_total + length(dfOvF_fast_cells);
    
    dfOvF_fast_cells_allsession = cat(2,dfOvF_fast_cells_allsession,dfOvF_fast_cells);
    dfOvF_slow_cells_allsession = cat(2,dfOvF_slow_cells_allsession,dfOvF_slow_cells);
    dfOvF_stay_cells_allsession = cat(2,dfOvF_stay_cells_allsession,dfOvF_stay_cells);   
end

dfOvF_fast_cells_acrosssession = mean(dfOvF_fast_cells_allsession);
stedfOvF_fast_cells_allsession = std(dfOvF_fast_cells_allsession)/sqrt(length(dfOvF_fast_cells_allsession));

dfOvF_slow_cells_acrosssession = mean(dfOvF_slow_cells_allsession);
stedfOvF_slow_cells_allsession = std(dfOvF_slow_cells_allsession)/sqrt(length(dfOvF_slow_cells_allsession));

dfOvF_stay_cells_acrosssession = mean(dfOvF_stay_cells_allsession);
stedfOvF_stay_cells_allsession = std(dfOvF_stay_cells_allsession)/sqrt(length(dfOvF_stay_cells_allsession));

dfOvF_fig = figure;
x = [1,2,3];
dfOvF_plot = [dfOvF_fast_cells_acrosssession,dfOvF_slow_cells_acrosssession,dfOvF_stay_cells_acrosssession];
dfOvF_ste_plot = [stedfOvF_fast_cells_allsession,stedfOvF_slow_cells_allsession,stedfOvF_stay_cells_allsession];
errorbar(x,dfOvF_plot,dfOvF_ste_plot,'.','LineStyle','none','linewidth', 1,...
    'Color',[0.8431 0.0980 0.1098],'MarkerSize',8,'MarkerEdgeColor',[0.8431 0.0980 0.1098],'MarkerFaceColor',[0.8431 0.0980 0.1098]);
xlim([0.5 3.5]); ylim([0 0.6]);
x1= [1,2,3];
set(gca,'XTick',x1,'XTicklabel',{'fast','slow','stationary'});
ylabel('df/F');
axis square;
dfOvF_fig.Units = 'centimeters';
dfOvF_fig.Position = [1 1 5.5 5];
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
%set(gca,'XTickLabelRotation',45);
fig_name = 'across_session_fastvsSlowvsStay_dfOvF';
path = 'Z:\Analysis\figures\figure6_motorized_running\';
%orient(dfOvF_behavStates,'landscape');
print(dfOvF_fig,[path,fig_name],'-r600','-depsc');


%% FR scatter plot

FR_fast_cells_allsession = []; %frame*cell
FR_slow_cells_allsession = [];
FR_stay_cells_allsession = [];
nPCs_total = 0;

for ii = 1:length(sessions)
    %load things
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    spk_output = load([image_analysis_dest,'deconv_wfastRunSmin\' sessions{ii}, '_spk_deconvolve_staynrun_seperate.mat']);
    FR_fast_cells = spk_output.FR_fast_cells; %1*n cells
    FR_slow_cells = spk_output.FR_slow_cells;
    FR_stay_cells = spk_output.FR_stay_cells;
    
    nPCs_total = nPCs_total + length(FR_fast_cells);
    
    FR_fast_cells_allsession = cat(2,FR_fast_cells_allsession,FR_fast_cells);
    FR_slow_cells_allsession = cat(2,FR_slow_cells_allsession,FR_slow_cells);
    FR_stay_cells_allsession = cat(2,FR_stay_cells_allsession,FR_stay_cells);   
end

FR_fast_cells_acrosssession = mean(FR_fast_cells_allsession);
steFR_fast_cells_allsession = std(FR_fast_cells_allsession)/sqrt(length(FR_fast_cells_allsession));

FR_slow_cells_acrosssession = mean(FR_slow_cells_allsession);
steFR_slow_cells_allsession = std(FR_slow_cells_allsession)/sqrt(length(FR_slow_cells_allsession));

FR_stay_cells_acrosssession = mean(FR_stay_cells_allsession);
steFR_stay_cells_allsession = std(FR_stay_cells_allsession)/sqrt(length(FR_stay_cells_allsession));

FR_fig = figure;
x = [1,2,3];
FR_plot = [FR_fast_cells_acrosssession,FR_slow_cells_acrosssession,FR_stay_cells_acrosssession];
FR_ste_plot = [steFR_fast_cells_allsession,steFR_slow_cells_allsession,steFR_stay_cells_allsession];
errorbar(x,FR_plot,FR_ste_plot,'.','LineStyle','none','linewidth', 1,...
    'Color',[0.1922 0.6392 0.3294],'MarkerSize',8,'MarkerEdgeColor',[0.1922 0.6392 0.3294],'MarkerFaceColor',[0.1922 0.6392 0.3294]);
xlim([0.5 3.5]); %ylim([0 0.6]);
x1= [1,2,3];
set(gca,'XTick',x1,'XTicklabel',{'fast','slow','stationary'});
ylabel('firing rate(Hz)');
axis square;
FR_fig.Units = 'centimeters';
FR_fig.Position = [1 1 5.5 5];
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
%set(gca,'XTickLabelRotation',45);
fig_name = 'across_session_fastvsSlowvsStay_FR';
path = 'Z:\Analysis\motorizedWheel_Analysis\running\across_sessions\deconv_wRunSmin\';
%orient(dfOvF_behavStates,'landscape');
%print(FR_fig,[path,fig_name],'-r600','-depsc');



