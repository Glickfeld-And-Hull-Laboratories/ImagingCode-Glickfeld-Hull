% plot an example raw fluorescence trace and denoised trace for the final
% paper.
% raster plot of 2 running trails(Fig 3B)
% spike rate of example session
% spike rate histogram of example session
%Don't need this for the actual analysis
%% supplemental: raw F and denoised signal e.g.
clear;
sessions = '190603_img1025'; 
days = '1025-190603_1';
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
image_analysis_dest = [image_analysis_base, sessions, '\'];
% behavior analysis results
behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days '\'];
%load data
filename = dir([image_analysis_dest 'getTC\' '*' '_TCave.mat']);
TCave = load([image_analysis_dest 'getTC\' filename.name]);
TCave = TCave.tc_avg;

behav_output = load([behav_dest days '_behavAnalysis.mat']);
frm_stay_cell = behav_output.frames_stay_cell;
len_stayTrials = cellfun(@length,frm_stay_cell);
trials = 1:length(frm_stay_cell);
longest = trials(len_stayTrials == max(len_stayTrials));
frame_range = frm_stay_cell{longest}(70:end-70);% get rid of the first 2ish second at the beginning and and end statioanry trials

threshold = -4;% the threshold you used in deconvolveCa function
load([image_analysis_dest sessions,'_spk_deconvolve_threshold' num2str(threshold) '.mat']);

picked_cell = 43; % chose in between 16,31,42,43
x = (1:1:length(frame_range))/30;
plotred = TCave(frame_range,picked_cell).*spk_logic(frame_range,picked_cell);
plotred(plotred==0) = NaN;
decon_fig = figure;
subplot(2,1,1); hold on;
plot(x,TCave(frame_range,picked_cell),'linewidth', 1,'color',[0.8431 0.0980 0.1098]); hold on;
plot(x,plotred,'ko','MarkerSize',4); ylabel('raw fluoscence');
xlim([0,43]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);

subplot(2,1,2)
plot(x,kernel(frame_range,picked_cell),'linewidth', 1,'color',[0.1922 0.6392 0.3294]); hold on; 
ylabel('denoised signal');xlabel('time(s)');
xlim([0,43]); ylim([0 4100]);
hold off

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
decon_fig.Units = 'centimeters';
decon_fig.Position = [1 1 15 8];
fig_name = 'eg_deconvolution_190603_img1025_longest_stationary_trial_cell_43';
path = 'Z:\Analysis\figures\figure3_2Prun_spike\';
orient(decon_fig,'landscape');
print(decon_fig,[path,fig_name],'-r600','-dpdf');

%% figure 3B: 
% spike logic
% plot spike logic
allcells = 1:size(TCave,2);
goodcells = setdiff(allcells,badPCs);
x1 = (20250:21150);
spk_logic_plot = spk_logic_cl(x1,:);
spk_logic_plot = spk_logic_plot';
spklogic_fig = figure;
colormap gray;
imagesc(imcomplement(spk_logic_plot));
trial1_start = 172; trial1_length = 330-172;
trial2_start = 625; trial2_length = 774-625;
r1 = rectangle('Position',[trial1_start 1 trial1_length 71] , 'EdgeColor',[0.5882    0.5882    0.5882], 'FaceColor', [0.5882    0.5882    0.5882]);
%uistack(r1,'bottom');
r2 = rectangle('Position',[trial2_start 1 trial2_length 71] , 'EdgeColor',[0.5882    0.5882    0.5882], 'FaceColor', [0.5882    0.5882    0.5882]);
%uistack(r2,'bottom');
set(gca,'Xtick',[0 150 300 450 600 750 900],'XTicklabel',{'0','5','10','15','20','25','30'});

xlabel('time(s)');ylabel('cell #');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
spklogic_fig.Units = 'centimeters';
spklogic_fig.Position = [1 3 8 5];
fig_name = 'eg_spk_logic_190603_img1025_20300-21100';
path = 'Z:\Analysis\figures\figure3_2Prun_spike\';
orient(spklogic_fig,'landscape')
%print(spklogic_fig,[path,fig_name],'-r600','-depsc');
savefig([path,fig_name]);
% the spike logics look super grey and blur when import eps or pdf, so used copy figure in matlab figure

%% figure 3B speed and average FR
clear;
sessions = '190603_img1025'; 
days = '1025-190603_1';
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
image_analysis_dest = [image_analysis_base, sessions, '\'];
threshold = -4;% the threshold you used in deconvolveCa function
load([image_analysis_dest 'deconvolution_staywRunSmin\' sessions,'_spk_deconvolve_staynrun_seperate' '.mat']);
frames = 20250:21150;%range of frames you want to plot
spk_eg = spk_logic(frames,:);
FR_eg = mean(spk_eg,2)*30;
std_FR = std(spk_eg,0,2)*30/sqrt(size(spk_eg,2));

behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days '\'];
behav_output = load([behav_dest days '_behavAnalysis.mat']);
speed = behav_output.speed;
speed_plot = speed(frames);
%smooth speed into every 100ms
bin = length(speed_plot)/3;
mean_spd_every100ms = zeros(1,length(speed_plot));
for x1 = 1:bin
    y = 3*x1 -2;
    mean_spd_every100ms(y:y+2) = mean(speed_plot(y:y+2));
    if abs(mean_spd_every100ms(y))<2 %on average moves less than 2 units per 0.1s, consider as stationary
     mean_spd_every100ms(y:y+2) = 0;
    end
end
% convert frames to seconds
x_plot = frames/30 - frames(1)/30;

FR_speed_fig = figure;
yyaxis left; hold on;
plot(x_plot,mean_spd_every100ms*2*3.1415926*7.5/128,'linewidth', 1,'color','k');
trial1_start = 172/30; trial1_length = (330-172)/30;
trial2_start = 625/30; trial2_length = (774-625)/30;
r1 = rectangle('Position',[trial1_start -10.3333*2*3.1415926*7.5/128 trial1_length 35+10.3333*2*3.1415926*7.5/128] , ...
    'EdgeColor',[0.9922 0.7059 0.3843], 'FaceColor', [0.9922 0.7059 0.3843]);
uistack(r1,'bottom');
r2 = rectangle('Position',[trial2_start -10.3333*2*3.1415926*7.5/128 trial2_length 35+10.3333*2*3.1415926*7.5/128] , ...
    'EdgeColor',[0.9922 0.7059 0.3843], 'FaceColor', [0.9922 0.7059 0.3843]);
uistack(r2,'bottom');
ylim([-5 35]);xlim([min(x_plot) max(x_plot)]);
ylabel('speed(cm/s)')
set(gca,'YColor','k');
yyaxis right;
shadedErrorBar(x_plot,FR_eg,std_FR,{'color',[0.1922 0.6392 0.3294]},{'Linewidth',1}); 
%p = plot(x_plot,avedfOvF_btm_cl(x1),'Color',);
ylabel('firing rate (Hz)'); set(gca,'YColor',[0.1922 0.6392 0.3294]);
ylim([-0.2 10.5]);
xlim([min(x_plot) max(x_plot)]);

xlabel('time(s)');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
FR_speed_fig.Units = 'centimeters';
FR_speed_fig.Position = [1 3 8 4];
fig_name = 'eg_speed_and_FR_190603_img1025_20300-21100';
path = 'Z:\Analysis\figures\figure3_2Prun_spike\';
orient(FR_speed_fig,'landscape')
print(FR_speed_fig,[path,fig_name],'-r600','-depsc');

%% figure 3C: not putting in the main figure here
% plot spike rates for each cell during stationary
spkrate_fig = figure;
scatter(1:length(FRstay_cell_cl), FRstay_cell_cl,60,'.', 'k');
ylabel('Spike Rate (Hz)');
xlabel('Cell #'); axis square;
ylim([0 1.5]);
%title ([sessions ' deconvolution']);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
spkrate_fig.Units = 'centimeters';
spkrate_fig.Position = [1 3 9 9];
fig_name = 'eg_spk_ratescatter_190603_img1025';
path = 'Z:\Analysis\figures\figure3_2Prun_spike\';
print(spkrate_fig,[path,fig_name],'-r600','-depsc');

%% histogram
spk_ratehist = figure;
histogram(FRstay_cell_cl,'Facecolor',[0.5 0.5 0.5],'BinWidth',0.1);
xlim([0 1.5]);
xlabel('Spike Rate during stationary(Hz)');
ylabel('number of cells');
axis square;
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
spk_ratehist.Units = 'centimeters';
spk_ratehist.Position = [1 3 4.5 4];
fig_name = 'eg_spk_ratehist_190603_img1025';
path = 'Z:\Analysis\figures\figure3_2Prun_spike\';
print(spk_ratehist,[path,fig_name],'-r600','-dpdf'); % depsc looks blur when import in coreldraw
