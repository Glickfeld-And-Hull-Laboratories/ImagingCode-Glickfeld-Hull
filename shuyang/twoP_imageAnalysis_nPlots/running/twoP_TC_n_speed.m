% plot the example 2P session for paper. Don't need this script for the
% actual analysis.
% choose a period of time that contains running trials and stationary
% trials, plot TC of selected cells, and plot speed with running trials shaded

clear;
sessions = '190603_img1025'; 
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
image_analysis_dest = [image_analysis_base, sessions, '\'];

% behavior analysis results 
days = '1025-190603_1';
behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days '\'];
color_code = {'b','r','k','c'};

%%
load('Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\190603_img1025\190603_img1025_deconvolution_thresh-4_TCave_cl.mat');
image_analysis_dest = [image_analysis_base, sessions, '\'];
% behavior analysis results
%behav_dest = ['Z:\Analysis\Airpuff_analysis\behavioral_analysis\' days{i} '\'];
behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days '\'];
%load data
filename = dir([image_analysis_dest 'getTC\' '*' '_TCave.mat']);
TCave = load([image_analysis_dest 'getTC\' filename.name]);
tc_avg = TCave.tc_avg;
allcells = 1:size(tc_avg,2);
goodcells = setdiff(allcells,badPCs);
behav_output = load([behav_dest days '_behavAnalysis.mat']);
speed = double(behav_output.speed);

plotTC_Jin(tc_avg(14200:14900,:),1,1,goodcells,1); title('14200-14900');
plotTC_Jin(tc_avg(20300:21100,:),1,1,goodcells(1:24),1); title('20300-21100,1-24');
plotTC_Jin(tc_avg(20300:21100,:),1,1,goodcells(25:47),1); title('20300-21100,25-47');
plotTC_Jin(tc_avg(20300:21100,:),1,1,goodcells(48:71),1); title('20300-21100,48-71');
plotTC_Jin(tc_avg(43800:44400,:),1,1,goodcells,1); title('43800-44400');
plotTC_Jin(tc_avg(56500:58300,:),1,1,goodcells,1); title('56500-58300');

x1 = (20250:21150);x_plot = (x1/30) - min(x1)/30;
plotTC_Jin_finalFig(x_plot,tc_avg(x1,:),1,1,[8,13,15,16,17,18,20,21,31,38,42,43,44,50,51,62],1); %title('20300-21100,selected cells');

speed_plot = speed(x1);
bin = length(speed_plot)/3;
mean_spd_every100ms = zeros(1,length(speed_plot));
for x1 = 1:bin
    y = 3*x1 -2;
    mean_spd_every100ms(y:y+2) = mean(speed_plot(y:y+2));
    if abs(mean_spd_every100ms(y))<2 %on average moves less than 2 units per 0.1s, consider as stationary
     mean_spd_every100ms(y:y+2) = 0;
    end
end
eg_speed = figure;plot(x_plot,mean_spd_every100ms*2*3.1415926*7.5/128,'linewidth', 1,'color','k');hold on;
trial1_start = 172/30; trial1_length = (330-172)/30;
trial2_start = 625/30; trial2_length = (774-625)/30;
r1 = rectangle('Position',[trial1_start -10.3333*2*3.1415926*7.5/128 trial1_length 35+10.3333*2*3.1415926*7.5/128] , ...
    'EdgeColor',[0.8000 0.9216 0.7725], 'FaceColor', [0.8000 0.9216 0.7725]);
uistack(r1,'bottom');
r2 = rectangle('Position',[trial2_start -10.3333*2*3.1415926*7.5/128 trial2_length 35+10.3333*2*3.1415926*7.5/128] , ...
    'EdgeColor',[0.8000 0.9216 0.7725], 'FaceColor', [0.8000 0.9216 0.7725]);
uistack(r2,'bottom');
xlabel('time(s)');ylabel('speed(cm/s)');
xlim([min(x_plot) max(x_plot)]);
ylim([-5 35]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
eg_speed.Units = 'centimeters';
eg_speed.Position = [1 3 7 4];
fig_name = 'eg_speed_190603_img1025_20300-21100';
path = 'Z:\Analysis\figures\figure2_2Prun\';
orient(eg_speed,'landscape')
print(eg_speed,[path,fig_name],'-r600','-depsc');

%figure;rectangle('Position',[0 0 10 10] , ...
%    'EdgeColor',[0.8097    0.9255    0.7831], 'FaceColor', [0.8097    0.9255    0.7831]);


