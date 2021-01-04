% plot an example raw fluorescence trace and processed trace from different
% spike identification methods.
% (deconvolution ( together or seperate) or first derivatives)

%% deconvolution together
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

%len_stayTrials = cellfun(@length,frm_stay_cell);
%trials = 1:length(frm_stay_cell);
%longest = trials(len_stayTrials == max(len_stayTrials));
frame_range = (45600:1:46360);

threshold = -4;% the threshold you used in deconvolveCa function
load([image_analysis_dest sessions,'_spk_deconvolve_threshold' num2str(threshold) '.mat']);

picked_cell = 43; % chose in between 16,31,42,43. this index is including all of the cells (good and bad)
x = (1:1:length(frame_range))/30;
plotred = TCave(frame_range,picked_cell).*spk_logic(frame_range,picked_cell);
plotred(plotred==0) = NaN;
decon_fig = figure;
subplot(2,1,1); hold on;
plot(x,TCave(frame_range,picked_cell),'linewidth', 1,'color',[0.8431 0.0980 0.1098]); hold on;
plot(x,plotred,'ko','MarkerSize',4); ylabel('raw fluoscence');
xlim([0,26]); ylim([0 5500]);
trial1_start = (46002-45600)/30; trial1_length = (46152-46002)/30;
trial2_start = (45685-45600)/30; trial2_length = (45833-45685)/30;
r1 = rectangle('Position',[trial1_start 1 trial1_length 5500] , 'EdgeColor',[0.9922 0.7059 0.3843], 'FaceColor', [0.9922 0.7059 0.3843]);
uistack(r1,'bottom');
r2 = rectangle('Position',[trial2_start 1 trial2_length 5500] , 'EdgeColor',[0.9922 0.7059 0.3843], 'FaceColor', [0.9922 0.7059 0.3843]);
uistack(r2,'bottom');
fs = get(gca,'XTickLabel');
set(gca,'XTickLabel',fs,'FontSize',7);

subplot(2,1,2)
plot(x,kernel(frame_range,picked_cell),'linewidth', 1,'color',[0.1922 0.6392 0.3294]); hold on; 
plot(x,spk(frame_range,picked_cell),'linewidth', 1,'color','k'); hold on;
ylabel('denoised and deconvolved trace');xlabel('time(s)');
xlim([0,26]); ylim([0 4100]);

r1 = rectangle('Position',[trial1_start 1 trial1_length 4000] , 'EdgeColor',[0.9922 0.7059 0.3843], 'FaceColor', [0.9922 0.7059 0.3843]);
uistack(r1,'bottom');
r2 = rectangle('Position',[trial2_start 1 trial2_length 4000] , 'EdgeColor',[0.9922 0.7059 0.3843], 'FaceColor', [0.9922 0.7059 0.3843]);
uistack(r2,'bottom');
hold off

fs = get(gca,'XTickLabel');
set(gca,'XTickLabel',fs,'FontSize',7);
decon_fig.Units = 'centimeters';
decon_fig.Position = [1 1 15 8];
fig_name = 'eg_deconvolution_190603_img1025_cell43_frm45600-46500';
path = [image_analysis_dest '\eg_traces\'];
orient(decon_fig,'landscape');
print(decon_fig,[path,fig_name],'-r600','-dpdf');

%% first derivative
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

%len_stayTrials = cellfun(@length,frm_stay_cell);
%trials = 1:length(frm_stay_cell);
%longest = trials(len_stayTrials == max(len_stayTrials));
frame_range = (45600:1:46360);

load('Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\190603_img1025\derivative\190603_img1025_spikes.mat');

picked_cell = 43; 
picked_cell_cl = 43-9; %43th cell in all cells = 34th cells in good cells because there're 9 cells getting deleted before this cell
x = (1:1:length(frame_range))/30;
plotred = TCave(frame_range,picked_cell).*spk_bi_cellmat(frame_range,picked_cell_cl);
plotred(plotred==0) = NaN;
deriv_fig = figure;
subplot(2,1,1); hold on;
plot(x,TCave(frame_range,picked_cell),'linewidth', 1,'color',[0.8431 0.0980 0.1098]); hold on;
plot(x,plotred,'ko','MarkerSize',4); ylabel('raw fluoscence');
xlim([0,26]); ylim([0 5500]);
trial1_start = (46002-45600)/30; trial1_length = (46152-46002)/30;
trial2_start = (45685-45600)/30; trial2_length = (45833-45685)/30;
r1 = rectangle('Position',[trial1_start 1 trial1_length 5500] , 'EdgeColor',[0.9922 0.7059 0.3843], 'FaceColor', [0.9922 0.7059 0.3843]);
uistack(r1,'bottom');
r2 = rectangle('Position',[trial2_start 1 trial2_length 5500] , 'EdgeColor',[0.9922 0.7059 0.3843], 'FaceColor', [0.9922 0.7059 0.3843]);
uistack(r2,'bottom');
fs = get(gca,'XTickLabel');
set(gca,'XTickLabel',fs,'FontSize',7);

% plot first derivative and threshold (2std of first derivative)
first_deriv = diff(TCave(:,picked_cell));
subplot(2,1,2);
plot(x,first_deriv(frame_range),'linewidth', 1,'color',[0.1922 0.6392 0.3294]);
hline(std_best(picked_cell_cl),'k');
r1 = rectangle('Position',[trial1_start 1 trial1_length 1500] , 'EdgeColor',[0.9922 0.7059 0.3843], 'FaceColor', [0.9922 0.7059 0.3843]);
uistack(r1,'bottom');
r2 = rectangle('Position',[trial2_start 1 trial2_length 1500] , 'EdgeColor',[0.9922 0.7059 0.3843], 'FaceColor', [0.9922 0.7059 0.3843]);
uistack(r2,'bottom');
xlim([0,26]); ylim([-500 2000]);

fs = get(gca,'XTickLabel');
set(gca,'XTickLabel',fs,'FontSize',7);

deriv_fig.Units = 'centimeters';
deriv_fig.Position = [1 1 15 8];
fig_name = 'eg_derivative_190603_img1025_cell43_frm45600-46360';
path = [image_analysis_dest '\eg_traces\'];
orient(deriv_fig,'landscape');
print(deriv_fig,[path,fig_name],'-r600','-dpdf');

%% deconvolution seperate
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
load('Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\190603_img1025\deconvolution_sep\190603_img1025_spk_deconvolve_staynrun_seperate.mat');
picked_cell = 43; 
picked_cell_cl = 34; 
behav_output = load([behav_dest days '_behavAnalysis.mat']);
frame_range = (45600:1:46360);
frms_run = cell2mat(behav_output.frames_run_cell);
frms_stay = cell2mat(behav_output.frames_stay_cell);
run_in_frm_range = intersect(frms_run,frame_range); % this is the index in the whole experiment
% get the index of running frames in frms_run
run_inx = zeros(1,length(run_in_frm_range));
for i = 1:length(run_in_frm_range)
    run_inx(i) = find(frms_run == run_in_frm_range(i)); % this inx will always be contineous
end

% do the same thing for stationary
stay_in_frm_range = intersect(frms_stay,frame_range); % this is the index in the whole experiment
% get the index of running frames in frms_run
stay_inx = zeros(1,length(stay_in_frm_range));
for i = 1:length(stay_in_frm_range)
    stay_inx(i) = find(frms_stay == stay_in_frm_range(i));
end
break_points_stay = find(diff(stay_in_frm_range)>1);
break_points_run = find(diff(run_in_frm_range)>1);

run_in_frm_range(1) < stay_in_frm_range(1) % the first part is stay
% stay_in_frm_range(1:85) = 45600:45684
% run_in_frm_range(1:149) = 45685:45833
% stay_in_frm_range(86:253) = 45834:46001
% run_in_frm_range(150:end) = 46002:46152
% stay_in_frm_range(254:453) = 46153:46352
% 46352-46354; 46374-46376;46460-46475 : move back and forth or stationary
% too short

% put back kernel and spk together in time course
kernel_stay_cell = (kernel_stay(:,picked_cell_cl))';
kernel_run_cell = (kernel_run(:,picked_cell_cl))';
kernel_plot = [kernel_stay_cell(stay_inx(1:break_points_stay(1))),kernel_run_cell(run_inx(1:break_points_run(1))),...
    kernel_stay_cell(stay_inx(break_points_stay(1)+1:break_points_stay(2))),...
    kernel_run_cell(run_inx(break_points_run(1)+1:end)),...
    kernel_stay_cell(stay_inx(break_points_stay(2)+1:end))];

spk_stay_cell = (spk_stay(:,picked_cell_cl))';
spk_run_cell = (spk_run(:,picked_cell_cl))';
spk_plot = [spk_stay_cell(stay_inx(1:break_points_stay(1))),spk_run_cell(run_inx(1:break_points_run(1))),...
    spk_stay_cell(stay_inx(break_points_stay(1)+1:break_points_stay(2))),...
    spk_run_cell(run_inx(break_points_run(1)+1:end)),...
    spk_stay_cell(stay_inx(break_points_stay(2)+1:end))];

spk_inx_pickedcell = spk_inx_neurons{picked_cell_cl};
spk_inx_logic = double(ismember(frame_range,spk_inx_pickedcell));
x1 = (1:1:length(frame_range))/30;
plotred = TCave(frame_range,picked_cell).*spk_inx_logic';
plotred(plotred==0) = NaN;

deconv_sep_fig = figure;
subplot(2,1,1); hold on;
plot(x1,TCave(frame_range,picked_cell),'linewidth', 1,'color',[0.8431 0.0980 0.1098]); hold on;
plot(x1,plotred,'ko','MarkerSize',4); ylabel('raw fluoscence');
xlim([0,26]); ylim([0 5500]);
trial1_start = (46002-45600)/30; trial1_length = (46152-46002)/30;
trial2_start = (45685-45600)/30; trial2_length = (45833-45685)/30;
r1 = rectangle('Position',[trial1_start 1 trial1_length 5500] , 'EdgeColor',[0.9922 0.7059 0.3843], 'FaceColor', [0.9922 0.7059 0.3843]);
uistack(r1,'bottom');
r2 = rectangle('Position',[trial2_start 1 trial2_length 5500] , 'EdgeColor',[0.9922 0.7059 0.3843], 'FaceColor', [0.9922 0.7059 0.3843]);
uistack(r2,'bottom');
fs = get(gca,'XTickLabel');
set(gca,'XTickLabel',fs,'FontSize',7);

subplot(2,1,2);
x2 = (1:1:length(kernel_plot))/30;
plot(x2,kernel_plot,'linewidth', 1,'color',[0.1922 0.6392 0.3294]); hold on; 
plot(x2,spk_plot,'linewidth', 1,'color','k'); hold on;
ylabel('denoised and deconvolved trace');xlabel('time(s)');
xlim([0,26]); ylim([0 4100]);

r1 = rectangle('Position',[trial1_start 1 trial1_length 4000] , 'EdgeColor',[0.9922 0.7059 0.3843], 'FaceColor', [0.9922 0.7059 0.3843]);
uistack(r1,'bottom');
r2 = rectangle('Position',[trial2_start 1 trial2_length 4000] , 'EdgeColor',[0.9922 0.7059 0.3843], 'FaceColor', [0.9922 0.7059 0.3843]);
uistack(r2,'bottom');
hold off

fs = get(gca,'XTickLabel');
set(gca,'XTickLabel',fs,'FontSize',7);
deconv_sep_fig.Units = 'centimeters';
deconv_sep_fig.Position = [1 1 15 8];
fig_name = 'eg_deconv_sep_190603_img1025_cell43_frm45600-46360';
path = [image_analysis_dest '\eg_traces\'];
orient(deconv_sep_fig,'landscape');
print(deconv_sep_fig,[path,fig_name],'-r600','-dpdf');

%% deconvolution stay use run smin
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
load('Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\190603_img1025\deconvolution_staywRunSmin\190603_img1025_spk_deconvolve_staynrun_seperate.mat');
picked_cell = 43; 
picked_cell_cl = 43-9; %43th cell in all cells = 34th cells in good cells because there're 9 cells getting deleted before this cell

behav_output = load([behav_dest days '_behavAnalysis.mat']);
frame_range = (45600:1:46360);
frms_run = cell2mat(behav_output.frames_run_cell);
frms_stay = cell2mat(behav_output.frames_stay_cell);

x = (1:1:length(frame_range))/30;
plotred = TCave(frame_range,picked_cell).*spk_logic(frame_range,picked_cell_cl);
plotred(plotred==0) = NaN;
deconv_sep_fig = figure;
subplot(2,1,1); hold on;
plot(x,TCave(frame_range,picked_cell),'linewidth', 1,'color',[0.8431 0.0980 0.1098]); hold on;
plot(x,plotred,'ko','MarkerSize',4); ylabel('raw fluoscence');
xlim([0,26]); ylim([0 5500]);
trial1_start = (46002-45600)/30; trial1_length = (46152-46002)/30;
trial2_start = (45685-45600)/30; trial2_length = (45833-45685)/30;
r1 = rectangle('Position',[trial1_start 1 trial1_length 5500] , 'EdgeColor',[0.9922 0.7059 0.3843], 'FaceColor', [0.9922 0.7059 0.3843]);
uistack(r1,'bottom');
r2 = rectangle('Position',[trial2_start 1 trial2_length 5500] , 'EdgeColor',[0.9922 0.7059 0.3843], 'FaceColor', [0.9922 0.7059 0.3843]);
uistack(r2,'bottom');
fs = get(gca,'XTickLabel');
set(gca,'XTickLabel',fs,'FontSize',7);

subplot(2,1,2);
plot(x,kernel(frame_range,picked_cell_cl),'linewidth', 1,'color',[0.1922 0.6392 0.3294]); hold on; 
plot(x,spk(frame_range,picked_cell_cl),'linewidth', 1,'color','k'); hold on;
ylabel('denoised and deconvolved trace');xlabel('time(s)');
xlim([0,26]); ylim([0 4100]);
r1 = rectangle('Position',[trial1_start 1 trial1_length 4000] , 'EdgeColor',[0.9922 0.7059 0.3843], 'FaceColor', [0.9922 0.7059 0.3843]);
uistack(r1,'bottom');
r2 = rectangle('Position',[trial2_start 1 trial2_length 4000] , 'EdgeColor',[0.9922 0.7059 0.3843], 'FaceColor', [0.9922 0.7059 0.3843]);
uistack(r2,'bottom');
hold off

fs = get(gca,'XTickLabel');
set(gca,'XTickLabel',fs,'FontSize',7);
deconv_sep_fig.Units = 'centimeters';
deconv_sep_fig.Position = [1 1 15 8];
fig_name = 'eg_deconv_wRunSmin_190603_img1025_cell43_frm45600-46360';
path = [image_analysis_dest '\eg_traces\'];
orient(deconv_sep_fig,'landscape');
print(deconv_sep_fig,[path,fig_name],'-r600','-dpdf');
