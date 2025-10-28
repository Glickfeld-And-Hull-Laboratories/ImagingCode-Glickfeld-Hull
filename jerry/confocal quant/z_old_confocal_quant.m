datapath = 'G:/home/ACh/Data/Histology/Confocal_IHC_DARTSOM/tiffs';
outpath = 'G:/home/ACh/Analysis/histology_quant';
cd(datapath);

%find all .tif files in the data folder that has "DART" in the filename
filePattern = fullfile(datapath,'*DART*.tif');
txtFilesList = dir(filePattern);
allfilenames = {txtFilesList.name}; %this gets a 1*n cell array

%convert the cell array into a string array of the same size
fns = strings(size(allfilenames));
[fns{:}] = allfilenames{:};


%% make image stacks from tif files
% All files are named DART_IHC_THxxxSx_x_zx_chxx
imlist = [fns(3) fns(7) fns(11) fns(15) fns(19) fns(23) fns(27) fns(31) fns(35) fns(39)];
im = imread(fns(19));

img_stack = uint8.empty;
for i = 1:length(imlist)
    this_image = imread(imlist(i));
    img_stack = cat(3,img_stack,this_image);
end

mip = max(img_stack,[],3);
clear img_stack this_image

%% testing some parameters

figure();
imshow(im);
figure();
imshow(mip);
figure();
imshow(mip,[0 510]);

%% segmentation
