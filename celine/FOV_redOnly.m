% %% day 1- create target images
clear global
clear all
close all

mouse = 'ACh11';
date = '210622';
redrun = '002';
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

data_r = squeeze(data(1,:,:,:));

data_r_avg = mean(data_r(:,:,:),3);

% red snapshot, registered to red channel
[out data_r_reg] = stackRegister(data_r,data_r_avg);


data_r_reg_avg = mean(data_r_reg,3);


figure; imagesc(data_r_reg_avg); title('red at 1040 nm'); 
colormap gray;
caxis([100 500])
print(fullfile(out_path, [date '_' mouse '_runs-' redrun '_red_FOV.pdf']),'-dpdf','-bestfit')

