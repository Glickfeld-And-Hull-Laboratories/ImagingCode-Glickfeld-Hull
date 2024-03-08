clear all
close all

% kernel AFTER MIP to further denoise
kernel = fspecial("gaussian",[5 5],1);
cellsizethresh = 500;

datafile = 'G:/home/ACh/Analysis/histology_quant/IHC03_confocal_TH/MIPs.mat';
outpath = 'G:/home/ACh/Analysis/histology_quant/IHC03_confocal_TH/masks';
mask_fig_path = 'G:/home/ACh/Analysis/histology_quant/IHC03_confocal_TH/masks/mask_pngs';

load(datafile);

% "img_data" is a nImageSeries x 3 cell. 1st column is red channel MIP
% for the corresponding image series, 2nd column is cyan channel, and 3rd
% column is series identifier

%%

% processing pipeline: in "confocal_zProj.m", each image plane was
% processed with a 5x5 gaussian filter with  = 1. MIP was done for each
% image SERIES. 

% This code chunk filters each MIP AGAIN with the same gaussian filter
% before segmentation.
size1 = size(img_data,1);
size2 = size(img_data,2);

img_processed = cell(size1,size2);

for i = 1:size(img_data,1)
    img_processed{i,1} = imfilter(img_data{i,1},kernel);
    img_processed{i,2} = img_data{i,2};
    img_processed{i,3} = img_data{i,3};
end

clear img_data

% all_masks = cell(size1,size2);

nCh = size2 - 1;
% everything below is included in the loop

% LOOP1: AROUND IMAGE SERIES "ser"
% LOOP2: LOOP THE SECOND CHANNEL "ch" 

% segmentation

% for ser = 1:size1
%     for ch = 1:nCh
%         img4segmenting = cat(3,img_processed{ser,ch},img_processed{ser,ch},img_processed{ser,ch});
%         mask_exp = zeros(2048,2048);
%         mask_all = zeros(2048,2048); % hard-coding image size
% 
%         for iter = 1:size(img4segmenting,3)
%             mask_temp = img4segmenting(:,:,iter);
%             mask_temp(find(mask_exp >= 1)) = 0;
%             bwOut = imCellEditInteractiveLG(mask_temp);
%             mask_all = mask_all+bwOut;
%             mask_exp = imCellBuffer(mask_all,3)+mask_all;
%             close all
%         end
% 
%         % gets rid of area smaller than 500 pixels
%         mask_cleaned = bwareaopen(mask_all,500);
%         mask_labeled = bwlabel(mask_cleaned);
%         nCells = max(unique(mask_labeled));
% 
%         % mask_path = fullfile(outpath,['mask_s' num2str(ser) '_ch' num2str(ch)]);
%         % save mask for reproducibility
% 
%         % save(fullfile(outpath,['mask_s_' num2str(ser) '_ch_' num2str(ch)]),'mask_labeled');
%         all_masks{ser,ch} = mask_labeled;
%         imshow(mask_labeled);
%         saveas(gcf,fullfile(mask_fig_path,['s_' num2str(ser) '_ch_' num2str(ch) '.png']));
%         close all;
%     end
%     all_masks{i,3} = img_processed{i,3};
% end
%below outside of loop
% save(fullfile(outpath,'masks'),'all_masks');
% clear all_masks;
load(fullfile(outpath,'masks.mat'));
%% processing
all_tcs = cell(23,3);
all_ttests = zeros(23,2);

%everything below in a new loop
for ser = 1:size1
    for ch = 1:nCh
        if ch == 1
            extract_img = 2;
        else
            extract_img = 1;
        end
        % gets intensity of each object in THE OTHER CHANNEL using the mask
        % obtained by segmenting the current channel 
        otherCh_intensity = stackGetTimeCourses(img_processed{ser,extract_img},all_masks{ser,ch});
        all_tcs{ser,extract_img} = otherCh_intensity;
        maskCh_intensity = stackGetTimeCourses(img_processed{ser,ch},all_masks{ser,ch});
        % rotates the mask counterclock wise 90 degrees 3 times
        mask_rot = rot90(all_masks{ser,ch});
        control_intensity = stackGetTimeCourses(img_processed{ser,extract_img},mask_rot);
        all_ttests(ser,extract_img) = ttest(otherCh_intensity,control_intensity);
    end
end



%below is outside of the loop

save(fullfile(outpath,'intensities'),'all_tcs');
save(fullfile(outpath,'ttests'),'all_ttests');


%% debugging & testing chunk

% segmentation pictures 
% kernel2 = fspecial("gaussian",[5 5],0.5);
% kernel3 = fspecial("gaussian",[3 3],0.5);
% kernel4 = fspecial("gaussian",[5 5],1); % this looked the best
% kernel5 = fspecial("gaussian",[3 3],1);
% kernel6 = fspecial("gaussian",[5 5],2);
% kernel7 = fspecial("gaussian",[3 3],2);
% 
% test_img = img_data{1,1};
% test1 = test_img;
% test2 = imfilter(test_img,kernel2);
% test3 = imfilter(test_img,kernel3);
% test4 = imfilter(test_img,kernel4);
% test5 = imfilter(test_img,kernel5);
% test6 = imfilter(test_img,kernel6);
% test7 = imfilter(test_img,kernel7);
% 
% figure;
% imshow(test1);
% figure;
% imshow(test2);
% figure;
% imshow(test3);
% figure;
% imshow(test4);
% figure;
% imshow(test5);
% figure;
% imshow(test6);
% figure;
% imshow(test7);
% 
% figure;
% imshow(img_data{1,1})
% figure;
% imshow(img_processed{1,1})

% rotating matrix
% A = {'a' 'b' 'c';'d' 'e' 'f';'g' 'h' 'i'};
% B = rot90(A,3);

% find 0s in an array
% a = [1 2 3 4 5 6 0 1 2 3 0 0 1 2];
% b = find(a == 0);
% a(a == 0) = 99;

% deprecated manual bwareaopen..
% mask_cell = bwlabel(mask_all);
% obj_idx = unique(mask_cell);
% obj_idx = obj_idx(2:end); % gets rid of 0
% % This is so the last cell size is also shown, because histcounts
% % processes arrays as histogram edges (so will return n-1 values if the
% % last endpoint is not added)
% temp_obj_idx = [obj_idx' max(obj_idx+1)]; 
% nPix = histcounts(mask_cell,temp_obj_idx); % creates a list of ROI sizes
% clear temp_obj_idx
% logi_idx = nPix > cellsizethresh; % logical array that access all cells larger than 500 pix
% nCells = sum(logi_idx);
% fprintf(['Segmented ' num2str(nCells) ' cells from series ' num2str(ser) ' channel ' num2str(ch) '.\n']);
% 
% mask_cleaned = mask_cell;
% nulls = find(logi_idx == 0);
% for k = 1:length(nulls)
%     mask_cleaned(mask_cleaned == nulls(k)) = 0;
% end

