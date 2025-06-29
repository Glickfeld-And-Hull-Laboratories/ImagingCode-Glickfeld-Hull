
mouse = 'i2207';
date = '250627';
redFolder = '000'; %enter the first three digits
redrun = '000'; %enter the last three digits for the red run
greenFolder = '001'; %enter the first three digits
greenrun = '000'; %enter the LAST three digits for the green run
depth='212.5';
experimentType = 'SST_YM90K';

%base= 'Z:/All_Staff/home/ACh/Aging/data/2p'
% base = 
% 'Z:/home/celine/Data/2p_data';
if computer == 'GLNXA64'
    base = '/home/cc735@dhe.duke.edu/GlickfeldLabShare/All_Staff/home/ACh/Data/2p_data';
    out_base = fullfile('/home/cc735@dhe.duke.edu/GlickfeldLabShare/All_Staff/home/ACh/Analysis/2p_analysis',experimentType);
else
    base = 'G:/home/ACh/Data/2p_data/';
    out_base = ['G:/home/ACh/Analysis/2p_analysis/' experimentType];
end

data_path = fullfile(base,mouse, date, redFolder);
out_path = fullfile(out_base,mouse, date, redFolder);

cd(data_path)
mkdir(out_path);

%load red run
nframes=1000;
load([redFolder '_000_' redrun '.mat'])
data = sbxread([redFolder '_000_' redrun],0, nframes);


data_g = squeeze(data(1,:,:,:));
data_r = squeeze(data(2,:,:,:));

data_g_avg = mean(data_g(:,:,:),3);
data_r_avg = mean(data_r(:,:,:),3);

% red snapshot, registered to red channel
[out, data_r_reg] = stackRegGPU(data_r,data_r_avg);
[outs, data_g_reg] = stackRegGPU(data_g,[],[],out);

data_g_reg_avg = mean(data_g_reg,3);
data_r_reg_avg = mean(data_r_reg,3);

% change for each figure
cd(out_path);
fig1=figure; imagesc(data_r_reg_avg); title([' ' depth ' red at 1040']); 
colormap gray;
caxis([50 800])
print(fullfile(out_path, [date '_' mouse  '_red_FOV.pdf']),'-dpdf','-bestfit')
saveas(fig1, 'redFOV.png')


%% green snapshot
nframes = 200;

data_path = fullfile(base,mouse, date, greenFolder);
cd(data_path)
load([greenFolder '_000_' greenrun '.mat'])
data_920 = sbxread([greenFolder '_000_' greenrun],0, nframes);


data_g_920 = squeeze(data_920(1,:,:,:));
data_g_920_avg = mean(data_g_920(:,:,:),3);

[out, data_g_reg_920] = stackRegGPU(data_g_920,data_g_920_avg);
regImg = mean(data_g_reg_920,3);

fig2=figure; imagesc(regImg);
colormap gray;
caxis([100 3500])
cd(out_path);
title([' ' depth ' green at 920']);
print(fullfile(out_path, [date '_' mouse '_FOV.pdf']),'-dpdf','-bestfit')
saveas(fig2, 'greenFOV.png');


%%
sz=size(regImg);
rgb = zeros(sz(1),sz(2),3);
rgb(:,:,1) = data_r_reg_avg./(max(data_r_reg_avg(:)));
rgb(:,:,2) = regImg./max(regImg(:));
fig3=figure; image(rgb);  movegui('center')
% change for each runc
title([' ' depth]);
print(fullfile(out_path, [date '_' mouse 'red_and_green_FOV.pdf']),'-dpdf','-bestfit')
saveas(fig3, 'redAndgreenFOV.png');


