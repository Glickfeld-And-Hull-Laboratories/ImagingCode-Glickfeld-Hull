% Analyze raw wide field imaging data. give you df/f and a ROI heatmap at
% the end.

%% SECTION ONE - assign pathnames and datasets to be analyzed/written. 
clear;
sessions = '220310_img1914_0.7mW_1'; 
days = '1914-220310-1358'; 
MWorks_source = 'Z:\home\shuyang\Data\behavior\WF_opto\';
image_source_base  = 'Z:\home\shuyang\Data\WF imaging\'; %location of permanently stored image files for retreiving meta data
image_dest_base    = 'Z:\home\shuyang\Analysis\WF_opto_analysis\'; %stores the data on crash in the movingDots analysis folder
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
cluster_dir = [image_dest_base sessions '\' sessions '_cluster.mat'];
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
save([image_dest,'_dfOvF.mat'],'dfOvF_btmbase');

% use the average of laseroff period as F ;; or the bottom 10% of laseroff period as F
MWorks_output = load([MWorks_source 'data-i' days '.mat' ]);
MWorks_output = MWorks_output.input;
laseroff = double(MWorks_output.nScansOff);
laseron = double(MWorks_output.nScansOn);
ntrials = double(MWorks_output.trialSinceReset);
laseroff_frame = zeros(1,ntrials);
laseroff_period = zeros(1,ntrials*(laseroff-20));
for i = 1:ntrials
    laseroff_frame(i) = 1+(laseroff+laseron)*(i-1); % the indexes of laser off, not the whole time of laseroff
    laseroff_period(1+(laseroff-20)*(i-1):((laseroff-20)*i)) = laseroff_frame(i)+20:laseroff_frame(i)+laseroff-1; % remove the rebound period
end
dfOvF_avglaseroff = zeros(size(data_tc));
for n = 1:size(data_tc,1)
    F_avglaseroff = data_tc(n,laseroff_period);
    dfOvF_avglaseroff(n,:) = (data_tc(n,:) - mean(F_avglaseroff))./mean(F_avglaseroff);
end
figure; plot(dfOvF_avglaseroff');
xlabel('frames (0.1s)');
ylabel('df/F');
save([image_dest,'_dfOvF.mat'],'dfOvF_avglaseroff','-append');

%% average dfOvF across trials, align to laser onset
MWorks_output = load([MWorks_source 'data-i' days '.mat' ]);
MWorks_output = MWorks_output.input;
dfOvF_temp = load([image_dest,'_dfOvF.mat']);
dfOvF = dfOvF_temp.dfOvF_btmbase;
%dfOvF = dfOvF_temp.dfOvF_avglaseroff;
% laseroff = 300;
% laseron = 50;
laseroff = double(MWorks_output.nScansOff);
laseron = double(MWorks_output.nScansOn);
%which frames are the first frame when laser is turned on each time: 151,311...
%ntrials = double(MWorks_output.trialSinceReset);
ntrials = 15;
laseron_frame = zeros(1,ntrials-1);
for i = 1:ntrials-1
    laseron_frame(i) = (laseroff+laseron)*(i-1)+1+laseroff; % the indexes of laser onset, not the whole time of laseron
end

laser_align_dfOvF = zeros(size(dfOvF,1),ntrials-1,30+laseron+15);% trial*frame*ROI, 1s before laser on and 1.5s after laser on 
for r = 1:size(dfOvF,1)
    for f = 1:length(laseron_frame)
        laser_align_dfOvF(r,f,:) = dfOvF(r,laseron_frame(f)-29:laseron_frame(f)+laseron+15);%-9
    end
end

align_fig = figure; 
align_fig.Units = 'centimeters';
align_fig.Position = [1 1 12 3.45];
x = double(1:(30+laseron+15));
x2 = (x./10)-0.9;
colorcode = {'r','b'};
for r = 1:size(dfOvF,1)
    shadedErrorBar(x2,squeeze(mean(laser_align_dfOvF(r,:,:),2)),squeeze(std(laser_align_dfOvF(r,:,:),0,2))/sqrt(size(laser_align_dfOvF(r,:,:),2)),colorcode{r});hold on; 
end
xlabel('time from laser on(s)');
ylabel('df/F');
title(sessions);
ylim([0 0.5]);
%ylim([-0.1 0.3]);
vline(0,'b'); vline(laseron/10,'b');
box off; hold off;
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',8,'FontName','Arial');
print(align_fig,[image_dest 'dfOvF_btmbase_align_laseron.pdf'],'-r600','-dpdf');
savefig([image_dest 'dfOvF_btmbase_align_laseron.fig']);
% print(align_fig,[image_dest 'dfOvFavglaseroff_align_laseron_frames.pdf'],'-r600','-dpdf');
% savefig([image_dest 'dfOvFavglaseroff_align_laseron_frames.fig']);


%%  SECTION VI - draw heatmap for ROI
%{
for ii = 1: length(sessions)
    image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    cluster = load([image_dest, '_cluster.mat']);
    cluster = cluster.cluster;
    raw_F = load([image_dest, '_raw_F.mat']);
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
    %saveas(ROI_fig, [image_dest, '_ROIplot'])
    fig_name = ['eg_session_ROI_',sessions{ii}];
    path = 'Z:\Analysis\figures\figure1_WF\';
    print(ROI_fig,[path,fig_name],'-r600','-dpdf');
    
end
%}

