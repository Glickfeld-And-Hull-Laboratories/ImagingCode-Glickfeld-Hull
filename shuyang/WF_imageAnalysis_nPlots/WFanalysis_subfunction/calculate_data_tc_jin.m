function [avg_img, roi_sz,data_tc] = calculate_data_tc_jin(meta_data_dir, cluster_dir)
%loads the raw data and uses that to generate a TC for each ROI. These TCs
%represent raw F and NOT a df/f. Saves data_tc and the frame size to the
%bx_outputs mat file. 

%load the meta_data and cluster info. meta_data used to load tiff files.
load(meta_data_dir, 'meta_data2');
load(cluster_dir);

%load tiff files
img = [];
for ii = 1:length(meta_data2)
    img_subset = readtiff(meta_data2{ii}(1).Filename);
    img = cat(3,img,img_subset);
    disp(['file #', num2str(ii), ' download complete']);
end
sz = size(img);
avg_img = mean(img,3);

%log the #of pixels in each ROI
roi_sz = sum(cluster.roi_mask,2);

% ----- cluster data by ROI and get raw TCs
disp('extracting ROIs separately');
data_tc = [];
for i = 1:cluster.num_cluster
    data_tc(i,:) = stackGetTimeCourses(img, reshape(cluster.roi_mask(i,:),[sz(1) sz(2)]));
end

%save sz and data_tc to bx_outputs
end