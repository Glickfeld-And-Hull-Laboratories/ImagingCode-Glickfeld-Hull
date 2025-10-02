clear all
date = '250929';%yymmdd
mouse = 'i2776' %ixxxx
greenFolder = '003'; %enter the first three digits


base = 'Z:\All_Staff\home\Xinyu\Data\2P_images'
data_path = fullfile(base, mouse,date,greenFolder);
cd(data_path)

nframes = 500;
load([greenFolder '_000_000.mat'])
data_920 = sbxread([greenFolder '_000_000'],0, nframes);

data_g_920 = squeeze(data_920(1,:,:,:));
data_g_920_avg = mean(data_g_920(:,:,:),3);

[out, data_reg] = stackRegister(data_g_920,data_g_920_avg);
reg_avg = mean(data_reg(:,:,:),3);

figure; imagesc(reg_avg);
colormap gray;
title([mouse,' avg FOV ',date])
%clim([200 500])