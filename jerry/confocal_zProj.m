% Quantification of confocal images for DART IHC03
% using bioformats in matlab https://docs.openmicroscopy.org/bio-formats/5.8.0/developers/matlab-dev.html
clear all

datapath_lif = char('G:/home/ACh/Data/Histology/Confocal_IHC_DARTSOM/DART_IHC.lif');
outpath = 'G:/home/ACh/Analysis/histology_quant/IHC03_confocal_TH';
cd(outpath);
%biof = bfopen(datapath_lif);
%biof = biof{:,[1 2]};
%c1 = biof(:,1:2);
load(fullfile(outpath,'img_series'));
biof = c1;
clear c1

zStackSize = 10;
ch_indices = [3 4];
n_channels = 4;
windowWidth = 3; % kernel size if applying conv filter
kernel = ones(windowWidth) / windowWidth^2; % conv kernel

% biof is a n*4 cell array where n is the # of image series in the .lif
% {n,1} contains the actual image
% {n,2} contains the image label
% {n,3} contains the colormap (if applicable)
% {n,4} contains associated metadata

% for example: 
series1 = biof{1,1};
% series1_plane19 = series1{19, 1};% actual image
% series1_label19 = series1{19, 2};% image label

% IF WANT COLORMAP
% series1_colorMaps = biof{1, 3};
% if (isempty(series1_colorMaps{1}))
%   colormap(gray);
% else
%   colormap(series1_colorMaps{1}(1,:));
% end

%% archived single series operation
% SAVE THE Maximum intensity projections
% imlist = [3 7 11 15 19 23 27 31 35 39];
% img_stack = uint8.empty;
%windowWidth = 3; %define filter size
% kernel = ones(windowWidth) / windowWidth^2;
% 
% for i = 1:length(imlist)
%     this_plane = series1{imlist(i),1};
%     %apply a convolution filter
%     filtered_img = conv2(this_plane, kernel,'same');
%     this_image = uint8(filtered_img);
%     img_stack = cat(3,img_stack,this_image);
% end
% take max projection of FILTERED images
% example is from series1 red channel
% mip = max(img_stack,[],3);
% figure();
% imshow(mip);

% clear img_stack filtered_img this_img this_plane

%% looping to get MIP for all 

nSeries = size(biof,1);
fprintf(['There are ' num2str(nSeries) ' image series total.\n']);
mip_storage = cell(nSeries,length(ch_indices)+1);
%erase_str = ['G:/home/ACh/Data/Histology/Confocal_IHC_DARTSOM/DART_IHC.lif; ','; plane 1/40; Z=1/10; C=1/4'];
% identify indices of images to be analyzed based on zStackSize & ch_indices

for i = 1:nSeries
    this_series = biof{i,1};
    metdat = this_series{1,2};
    seriesID = string(metdat(63:71)); % hardcoded
    mip_storage{i,3} = seriesID;
    for j = 1:length(ch_indices)
        start_point = ch_indices(j);
        end_point = start_point + (zStackSize-1) * n_channels;
        img_indices = [start_point:n_channels:end_point];
        this_stack = uint8.empty;
        for k = 1:length(img_indices)
            this_image = this_series{img_indices(k),1};
            filtered_img = uint8(conv2(this_image, kernel,'same'));
            this_stack = cat(3,this_stack,filtered_img);
        end
        mip = max(this_stack,[],3);
        mip_storage{i,j} = mip;
    end
end

%% see if mip_storage was calculated correctly

for q = 1:length(mip_storage)-18
    figure();
    imshow(mip_storage{q,1});
    figure();
    imshow(mip_storage{q,2});
end

%% If all looks good, save

save(fullfile(outpath,'MIPs'),"mip_storage");

%% DEBUGGING CHUNK
% archived convolution filter test
% windowWidth = 3; %define filter size
% kernel = ones(windowWidth) / windowWidth^2;
% outputImage = conv2(series1_plane19, kernel,'same');
% outputImage = uint8(outputImage);
% imshow(outputImage);

% archived label ID 
% for i = 1:19
%     fprintf([series1{i,2} '\n']);
% end
% class(series1{1,2})

% save image as matrices in cells
% cell_test = cell(5,7);
% a = 5;
% b = 7;
% for i = 1:a
%     for j = 1:b
%         cell_test{i,j} = rand(10,10);
%     end
% end

% see if mip_storage was calculated correctly
% for q = 1:length(mip_storage)-15
%     figure();
%     imshow(mip_storage{q,1});
% end

