%% day 1- create target images
clear all
mouse = 'i472';
date = '210222';
run = '003';
nframes = 1000;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Data\2P_images';
data_path = fullfile(base,mouse,date,run);
cd(data_path)
load([run '_000_000.mat'])
data = sbxread([run '_000_000'],0, nframes);
data_g = squeeze(data(1,:,:,:));
data_r = squeeze(data(2,:,:,:));

data_r_avg = mean(data_r(:,:,1:100),3);
[out data_r_reg] = stackRegister(data_r,data_r_avg);
[outs data_r_reg] = stackRegister_MA(data_r,[],[],out);

data_g_reg_avg = mean(data_g_reg,3);
data_r_reg_avg = mean(data_r_reg,3);

out_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\';
out_path = fullfile(out_base,[date '_' mouse], [date '_' mouse '_' run]);
mkdir(out_path)

figure; imagesc(data_g_reg_avg); title('green- 920nm'); truesize
print(fullfile(out_path, [date '_' mouse '_runs-' run '_greenFOV.pdf']),'-dpdf','-bestfit')
figure; imagesc(data_r_reg_avg); title('red- 920nm'); truesize
print(fullfile(out_path, [date '_' mouse '_runs-' run '_redFOV.pdf']),'-dpdf','-bestfit')
sz = size(data_g_reg_avg);
rgb = zeros(sz(1),sz(2),3);
rgb(:,:,1) = data_r_reg_avg./max(data_r_reg_avg(:));
rgb(:,:,2) = data_g_reg_avg./max(data_g_reg_avg(:));
figure; imagesc(rgb); title('red- 920nm'); truesize
print(fullfile(out_path, [date '_' mouse '_runs-' run '_rgbFOV.pdf']),'-dpdf','-bestfit')

%% day n create match images
clear all
mouse = 'i70220';
date = '201231';
run = '001';
nframes = 500;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\2P_Imaging';
data_path = fullfile(base,mouse,date,run);
cd(data_path)
for i = 1:3
load([run '_000_00' num2str(i) '.mat'])
data = sbxread([run '_000_00' num2str(i)],0, 500);
data_g = squeeze(data(1,:,:,:));
data_r = squeeze(data(2,:,:,:));

data_g_avg = mean(data_g(:,:,1:100),3);
[out data_g_reg] = stackRegister(data_g,data_g_avg);
[outs data_r_reg] = stackRegister_MA(data_r,[],[],out);

data_g_reg_avg = mean(data_g_reg,3);
data_r_reg_avg = mean(data_r_reg,3);


figure; imagesc(data_g_reg_avg); title(['green - 00' num2str(i)]); truesize
figure; imagesc(data_r_reg_avg); title(['red - 00' num2str(i)]); truesize

end
