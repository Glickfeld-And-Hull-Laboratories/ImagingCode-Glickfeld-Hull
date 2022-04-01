function WF_draw_ROIs_for_movingDots_Jin(sessions, image_source, image_dest)
%Use this function to manually draw ROIs. 

% ---------------
image_source  = [image_source '\' sessions '_MMStack_Pos0.ome.tif'];
%----loads the a sample of the imaging data
img = [];
for ii = 1:50
    img_temp = double(imread(image_source,ii)); %imread: read the first image in the file, double: double precision
    img = cat(3,img,img_temp);%overlap the images
end

%alternative method for selecting image sample
% for i = 1:20
%     for ii = 51:80;
%         img_temp = double(imread(image_source,ii));
%         img = cat(3,img,img_temp);
%     end
%     if min(mean(mean(img,2),1))>0
%         break
%     end
% end

% calculate average and other states
avg_img = squeeze(mean(img,3));
std_img = std(img,[], 2);
all_sd = std(img(:));
clear img;

%plot the avg_img to draw the ROIs
f = figure; im = imagesc(avg_img); axis ij; colormap jet; colormap jet; shading flat;
num_cluster = input('Enter no of clusters:');

for i=1:num_cluster
    roi(i) =  impoly; %You can add vertices and adjust the size and position of the polygon by using the mouse.
    disp(['Finished ROI ' num2str(i)]);%disp: make the command window print something
end

%save size/location information for each ROI
roi_mask = [];
roi_position = {};
for i=1:length(roi)
    c_mask = roi(i).createMask(im);
    roi_mask(i,:) = c_mask(:);
    roi_position{i} = roi(i).getPosition;
    roi_position{i} = [roi_position{i}; roi_position{i}(1,:)];
end
cluster.roi_mask =roi_mask;
cluster.roi_mask_combined = sum(roi_mask,1);
cluster.roi_position = roi_position;
cluster.num_cluster = num_cluster;
close gcf;
save([image_dest '_cluster.mat'], 'cluster');
return
