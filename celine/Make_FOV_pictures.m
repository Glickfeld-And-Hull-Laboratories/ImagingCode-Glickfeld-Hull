%% day 1- create target images
clear global
clear all
close all

mouse = 'TH19';
date = '210422';
run = '004';
nframes = 1000;
%base = 'Z:\home\celine\Data\2p_data\';
base = 'Z:\home\ACh\Data\2p_data\';
data_path = fullfile(base,mouse, date, run);


cd(data_path)
load([run '_000_000.mat'])

data = sbxread([run '_000_000'],0, nframes);

data_g = squeeze(data(1,:,:,:));
data_r = squeeze(data(2,:,:,:));


data_g_avg = mean(data_g(:,:,:),3);
data_r_avg = mean(data_r(:,:,:),3);

%plot both averages so that you can choose which one to register to
%(whichever has cleaner average)
figure
imagesc(data_g_avg);
title('green raw average');
colorbar
figure
imagesc(data_r_avg);
title('red raw average');
colorbar

out_base = 'Z:\home\ACh\Analysis\2p_analysis';
out_path = fullfile(out_base,mouse, date, run);
mkdir(out_path);
%% to get green channel from a seperate file

run = '003';
nframes = 200;

%base = 'Z:\home\celine\Data\2p_data\';
base = 'Z:\home\ACh\Data\2p_data\';
data_path = fullfile(base,mouse, date, run);


cd(data_path)
load([run '_000_002.mat'])

data_920 = sbxread([run '_000_002'],0, nframes);

data_g_920 = squeeze(data_920(1,:,:,:));
data_g_920_avg = mean(data_g_920(:,:,:),3);

nframes = 200;
[out data_g_reg_920] = stackRegister(data_g_920,data_g_920_avg);
regImg = mean(data_g_reg_920,3);

figure; imagesc(regImg);
colormap gray
caxis([200 1000])
title('green, registered at 920 nm');
print(fullfile(out_path, [date '_' mouse '_runs-' run '_green_reg_920.pdf']),'-dpdf','-bestfit')

%% to register to green
[out data_g_reg] = stackRegister(data_g,data_g_avg);
[outs data_r_reg] = stackRegister_MA(data_r,[],[],out);

data_g_reg_avg = mean(data_g_reg,3);
data_r_reg_avg = mean(data_r_reg,3);

%% to register to red
[out data_r_reg] = stackRegister(data_r,data_r_avg);
[outs data_g_reg] = stackRegister_MA(data_g,[],[],out);

data_g_reg_avg = mean(data_g_reg,3);
data_r_reg_avg = mean(data_r_reg,3);

%% to register to green from seperate file

nframes=1000;
[out,data_g_reg] = stackRegister(data_g,regImg);
[outs data_r_reg] = stackRegister_MA(data_r,[],[],out);
data_g_reg_avg = mean(data_g_reg,3);
data_r_reg_avg = mean(data_r_reg,3);

%% plotting



%to plot registered images of green and red
figure; imagesc(data_g_reg_avg); title('green, registered'); truesize;colorbar;
print(fullfile(out_path, [date '_' mouse '_runs-' run '_greenFOV.pdf']),'-dpdf','-bestfit')

figure; imagesc(data_r_reg_avg); title('red, registered'); truesize;colorbar;
print(fullfile(out_path, [date '_' mouse '_runs-' run '_redFOV.pdf']),'-dpdf','-bestfit')


% to print greyscale images of the registered activity 
figure; imagesc(data_g_reg_avg); title('green, registered'); truesize;
colormap gray
print(fullfile(out_path, [date '_' mouse '_runs-' run '_green_greyscale.pdf']),'-dpdf','-bestfit')

figure; imagesc(data_r_reg_avg); title('red, registered'); truesize;caxis([184 500]);
colormap gray
print(fullfile(out_path, [date '_' mouse '_runs-' run '_red_greyscale.pdf']),'-dpdf','-bestfit')

%% to plot red and green together

sz = size(data_g_reg);
rgb = zeros(sz(1),sz(2),3);
rgb(:,:,1) = data_r_reg_avg./max(data_r_reg_avg(:));
rgb(:,:,2)= data_g_reg_avg./max(data_g_reg_avg(:));
figure; image(rgb);



% data_r_reg_avg./data_r_reg_avg(:)

