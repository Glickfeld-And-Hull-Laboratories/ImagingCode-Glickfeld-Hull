% %% day 1- create target images
clear global
clear all
close all

mouse = 'ACH11';
date = '210622';
redrun = '007'; %enter the first three digits for the red run, assumes the
%last three at 000 - if not, change that on lines 23 and 24
greenrun = '005'; %enter the first three digits for the green run, 
%assumes the last three at 000 - if not change that on lines 52 and 53

base = 'Z:\home\ACh\Data\2p_data\';
data_path = fullfile(base,mouse, date, redrun);
cd(data_path)

out_base = 'Z:\home\ACh\Analysis\2p_analysis';
out_path = fullfile(out_base,mouse, date, redrun);
mkdir(out_path);

%load red run
nframes=1000;
load([redrun '_000_000.mat'])
data = sbxread([redrun '_000_000'],0, nframes);

data_g = squeeze(data(1,:,:,:));
data_r = squeeze(data(2,:,:,:));

data_g_avg = mean(data_g(:,:,:),3);
data_r_avg = mean(data_r(:,:,:),3);

% red snapshot, registered to red channel
[out data_r_reg] = stackRegister(data_r,data_r_avg);
[outs data_g_reg] = stackRegister_MA(data_g,[],[],out);

data_g_reg_avg = mean(data_g_reg,3);
data_r_reg_avg = mean(data_r_reg,3);


figure; imagesc(data_r_reg_avg); title('red at 1040 nm'); 
colormap gray;
print(fullfile(out_path, [date '_' mouse  '_red_FOV.pdf']),'-dpdf','-bestfit')


%% green snapshot
nframes = 200;

base = 'Z:\home\ACh\Data\2p_data\';
data_path = fullfile(base,mouse, date, greenrun);

cd(data_path)
load([ greenrun '_000_000.mat'])
data_920 = sbxread([greenrun '_000_000' ],0, nframes);

data_g_920 = squeeze(data_920(1,:,:,:));
data_g_920_avg = mean(data_g_920(:,:,:),3);

nframes = 200;
[out data_g_reg_920] = stackRegister(data_g_920,data_g_920_avg);
regImg = mean(data_g_reg_920,3);

figure; imagesc(regImg);
colormap gray;
title('green at 920 nm');
print(fullfile(out_path, [date '_' mouse '_green_FOV.pdf']),'-dpdf','-bestfit')


%%
sz=size(regImg);
rgb = zeros(sz(1),sz(2),3);
rgb(:,:,1) = data_r_reg_avg./max(data_r_reg_avg(:));
rgb(:,:,2) = regImg./max(regImg(:));
figure; image(rgb);  movegui('center')
title('Green-920 + Red-1040')
print(fullfile(out_path, [date '_' mouse 'red_and_green_FOV.pdf']),'-dpdf','-bestfit')


