% clear all
% close all

% kernel AFTER MIP to further denoise
kernel = fspecial("gaussian",[5 5],1);
cellsizethresh = 300;

datafile = 'G:/home/ACh/Analysis/histology_quant/IHC03_confocal_TH/MIPs.mat';
outpath = 'G:/home/ACh/Analysis/histology_quant/IHC03_confocal_TH/masks';
mask_fig_path = 'G:/home/ACh/Analysis/histology_quant/IHC03_confocal_TH/masks/mask_pngs';
mask_fig_path_vSize = 'G:/home/ACh/Analysis/histology_quant/IHC03_confocal_TH/masks/mask_pngs/px300';
load(datafile);

mkdir(mask_fig_path_vSize);
% "img_data" is a nImageSeries x 3 cell. 1st column is red channel MIP
% for the corresponding image series, 2nd column is cyan channel, and 3rd
% column is series identifier

%%

% processing pipeline: in "confocal_zProj.m", each image plane was
% processed with a 5x5 gaussian filter with  = 1. MIP was done for each
% image SERIES. 

% This code chunk filters each ch1 MIP AGAIN with the same gaussian filter
% before segmentation, BUT DOESN'T DO SO FOR CH2 MIP
size1 = size(img_data,1);
size2 = size(img_data,2);
% 
% img_processed = cell(size1,size2);
% 
% for i = 1:size(img_data,1)
%     img_processed{i,1} = imfilter(img_data{i,1},kernel);
%     img_processed{i,2} = img_data{i,2};
%     img_processed{i,3} = img_data{i,3};
% end
% 
clear img_data;
raw_masks = cell(size1,size2);
nCh = size2 - 1;
% % everything below is included in the loop
% 
% % LOOP1: AROUND IMAGE SERIES "ser"
% % LOOP2: LOOP THE SECOND CHANNEL "ch" 
% 
% % segmentation
% 
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
%         % gets rid of area smaller than 300 pixels
%         % mask_cleaned = bwareaopen(mask_all,300);
%         % mask_labeled = bwlabel(mask_cleaned);
%         % nCells = max(unique(mask_labeled));
% 
%         % mask_path = fullfile(outpath,['mask_s' num2str(ser) '_ch' num2str(ch)]);
%         % save mask for reproducibility
% 
%         % save(fullfile(outpath,['mask_s_' num2str(ser) '_ch_' num2str(ch)]),'mask_labeled');
%         raw_masks{ser,ch} = mask_all;
%         % imshow(mask_all);
%         % saveas(gcf,fullfile(mask_fig_path_vSize,['s_' num2str(ser) '_ch_' num2str(ch) '.png']));
%         % close all;
%     end
%     raw_masks{i,3} = img_processed{i,3};
% end
% % below outside of loop
% save(fullfile(outpath,'masks_unprocessed'),'raw_masks');
% clear raw_masks;

load(fullfile(outpath,'masks_unprocessed.mat'));

%% mask info display and pre-processing

% FLASHING IMAGES WARNING WHEN RUNNING THIS CODE CHUNK

final_masks = cell(23,3);
nCells_allCh = zeros(23,2);
% get rid of area smaller than 300 pixels, make labeld masks, get total
% number of cells, 
for ser = 1:size1
    for ch = 1:nCh
        mask_cleaned = bwareaopen(raw_masks{ser,ch},300);
        mask_labeled = bwlabel(mask_cleaned);
        nCells = max(unique(mask_labeled));
        nCells_allCh(ser,ch) = nCells;
        % save masks as png
        % save(fullfile(mask_fig_path_vSize,['mask_s_' num2str(ser) '_ch_' num2str(ch)]),'mask_labeled');
        final_masks{ser,ch} = mask_labeled;
        imshow(mask_labeled);
        saveas(gcf,fullfile(mask_fig_path_vSize,['s_' num2str(ser) '_ch_' num2str(ch) '.png']));
        close all;
    end
    final_masks{ser,3} = raw_masks{ser,3};
end

%save the final processed masks and cell counts after the loop
save(fullfile(mask_fig_path_vSize,['final_masks_' num2str(cellsizethresh)]),'final_masks');
save(fullfile(mask_fig_path_vSize,['cell_counts_' num2str(cellsizethresh)]),'nCells_allCh');

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

