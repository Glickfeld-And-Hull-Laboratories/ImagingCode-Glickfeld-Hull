% Analyze raw wide field imaging data. give you df/f and a ROI heatmap at
% the end.

%% SECTION ONE - assign pathnames and datasets to be analyzed/written. 
clear;
sessions = '220126_img1911_1mMGabazine1uMAlexa1.5ul_1';
image_source_base  = 'Z:\home\shuyang\Data\WF imaging\'; %location of permanently stored image files for retreiving meta data
%image_dest_base    = 'Z:\home\shuyang\Analysis\motorizedWheel_Analysis\WF\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
image_dest_base    = 'Z:\home\shuyang\Analysis\motorizedWheel_Analysis\WF_drug\imaging_analysis\'; 
image_dest = [image_dest_base sessions '\' sessions];
image_source = [image_source_base, sessions];

%% SECTION TWO - Uses a gui to allow user to draw ROIs 
if exist([image_dest_base sessions], 'file') ~= 7
    mkdir(image_dest_base, sessions);
end
WF_draw_ROIs_for_movingDots_Jin(sessions, image_source, image_dest); %automatically saves the ROIs

%% SECTION THREE - find/save frame times and meta data collected from the tiff files 

[all_files, meta_data, meta_data2] = obtain_tif_meta_data_Jin(image_source);
frame_times = get_frame_time_by_movie_info(meta_data);
dest =  [image_dest_base sessions '\' sessions];
if exist([image_dest_base sessions], 'file') ~= 7
    mkdir(image_dest_base, sessions);
end
save([dest '_frame_times.mat'],  'frame_times');
save([dest '_meta_data.mat'],  'meta_data', 'meta_data2', 'all_files');

%% SECTION FOUR - get raw F

meta_data_dir = [image_dest_base sessions '\' sessions '_meta_data.mat'];
%cluster_dir = [image_dest_base sessions '\' sessions '_cluster.mat'];
cluster_dir = 'Z:\home\shuyang\Analysis\motorizedWheel_Analysis\WF\imaging_analysis\220126_img1911_1\220126_img1911_1_cluster.mat';% use the ROI from the session before drug infusion
[avg_img, roi_sz,data_tc] = calculate_data_tc_jin(meta_data_dir, cluster_dir);
rawF_dir = [image_dest, '_raw_F.mat'];
save(rawF_dir, 'data_tc', 'roi_sz', 'avg_img');%automatically saves data_tc to bxOutputs

%% calculate df/F
% use bottom 10% as F
raw_F = load([image_dest, '_raw_F.mat']);
data_tc = raw_F.data_tc;
F_btm = zeros(1,size(data_tc,1));
dfOvF_btmbase = zeros(size(data_tc));
frames = 1:size(data_tc,2);
for n = 1:size(data_tc,1)
    data_tc_sort = sort(data_tc(n,:),'ascend');
    btm = data_tc_sort(1:floor(length(data_tc_sort)*0.1));
    F_btm(n) = mean(btm);
    dfOvF_btmbase(n,:) = (data_tc(n,:) - F_btm(n))/F_btm(n);
end
figure; plot(dfOvF_btmbase');
xlabel('frames (0.1s)');
ylabel('df/F');
savefig([image_dest, 'dfOvF.fig']);
save([image_dest,'_dfOvF.mat'],'dfOvF_btmbase','F_btm');


%% calculate df/F - using F from the control condition
%{
% use bottom 10% as F
raw_F = load([image_dest, '_raw_F.mat']);
data_tc = raw_F.data_tc;
dfOvF_ctrlbase = zeros(size(data_tc));
control = load('Z:\home\shuyang\Analysis\motorizedWheel_Analysis\WF\imaging_analysis\220107_img1907_1\220107_img1907_1_dfOvF.mat');
F = control.F_btm;
frames = 1:size(data_tc,2);
for n = 1:size(data_tc,1)
    dfOvF_ctrlbase(n,:) = (data_tc(n,:) - F(n))/F(n);
end
figure; plot(dfOvF_ctrlbase');
xlabel('frames (0.1s)');
ylabel('df/F');
savefig([image_dest, 'dfOvF_ctrlF.fig']);
save([image_dest,'_dfOvF_ctrlbase.mat'],'dfOvF_ctrlbase');
%}
%%  SECTION VI - draw heatmap for ROI

image_dest = [image_dest_base sessions '\' ];
%cluster = load([image_dest, '_cluster.mat']);
cluster = load(cluster_dir);
cluster = cluster.cluster;
raw_F = load([image_dest, sessions, '_raw_F.mat']);
avg_img = raw_F.avg_img;
ROI_fig = figure;
imagesc(avg_img);
shading flat; hold on;
for i=1:cluster.num_cluster
    line(cluster.roi_position{i}(:,1), cluster.roi_position{i}(:,2),'Color','k','LineWidth', 1.15);
    text(mean(cluster.roi_position{i}(:,1)),mean(cluster.roi_position{i}(:,2)), ...
        num2str(i), 'color', 'k', 'FontSize', 18, 'Color','k');
end
%title(['selected ROIs ',sessions{ii}]);
colormap jet; %colomap jet makes the image brighter
saveas(ROI_fig, [image_dest, sessions 'ROIplot.fig']);
%fig_name = ['eg_session_ROI_',sessions];
%path = 'Z:\Analysis\figures\figure1_WF\';
%print(ROI_fig,[path,fig_name],'-r600','-dpdf');




