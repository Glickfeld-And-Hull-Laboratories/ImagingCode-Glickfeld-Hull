%% day 1- create target images
clear all
mouse = 'i70220';
date = '201231';
run = '002';
nframes = 500;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\2P_Imaging';
data_path = fullfile(base,mouse, [date '_' mouse], run);
cd(data_path)
load([run '_000_000.mat'])
data = sbxread([run '_000_000'],0, nframes);
data_g = squeeze(data(1,:,:,:));
data_r = squeeze(data(2,:,:,:));

data_g_avg = mean(data_g(:,:,1:100),3);
[out data_g_reg] = stackRegister(data_g,data_g_avg);
[outs data_r_reg] = stackRegister_MA(data_r,[],[],out);

data_g_reg_avg = mean(data_g_reg,3);
data_r_reg_avg = mean(data_r_reg,3);

out_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\';
out_path = fullfule(out_base,[date '_' mouse], [date '_' mouse '_runs-' run]);

figure; imagesc(data_g_reg_avg); title('green - 000'); truesize
print(fullfile(out_path, [date '_' mouse '_runs-' run '_greenFOV.pdf']),'-dpdf','-bestfit')
figure; imagesc(data_r_reg_avg); title('red'); truesize
print(fullfile(out_path, [date '_' mouse '_runs-' run '_redFOV.pdf']),'-dpdf','-bestfit')

%% day n create match images
clear all
mouse = 'i70220';
date = '201231';
run = '001';
nframes = 500;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\2P_Imaging';
data_path = fullfile(base,mouse, [date '_' mouse], run);
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
