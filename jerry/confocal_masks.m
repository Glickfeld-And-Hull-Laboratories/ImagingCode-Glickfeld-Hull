clear all
close all

datafile = 'G:/home/ACh/Analysis/histology_quant/IHC03_confocal_TH/MIPs.mat';
outpath = 'G:/home/ACh/Analysis/histology_quant/IHC03_confocal_TH/masks';

load(datafile);
img_data = mip_storage;
clear mip_storage;

%%

imgtest = cat(3,img_data{1,1},img_data{1,1},img_data{1,1});

for i = 1:length(imgtest,3)
    bwOut = imCellEditInteractiveLG(imgtest);

end