clear all;
clear global;
close all;
clc;

%% get path names
ref_date = '240908';
date = '240908';
time = strvcat('1610');
alignToRef = 1;
ImgFolder = strvcat('005');
mouse = 'i2581';
nrun = size(ImgFolder,1);
frame_rate = 15.5;
ref_str = 'runs-001';
run_str = catRunName(ImgFolder, nrun);
tj_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\2P_Imaging';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %MAKE SURE TO SET THIS
behav_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';

%%
data = [];
clear temp
trial_n = [];
irun = 1;
CD = fullfile(tj_fn, [mouse '\' date '_' mouse '\' ImgFolder(irun,:)]);
cd(CD);
imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
load(imgMatFile);
nframes = info.config.frames;
fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n'])
data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);

data_temp = squeeze(data_temp);
data = cat(3,data,data_temp);
trial_n = [trial_n nframes];

clear data_temp
%%
data_avg = squeeze(mean(reshape(data,[529 796 100 11]),3));
figure;
for i = 1:9
    subplot(3,3,i)
    imagesc(data_avg(:,:,i))
end

data_by_z = reshape(data,[529 796 100 11]);

%%
data_reg = [];
outs = [];

for i = 1:11;
    [out, reg] = stackRegister(data_by_z(:,:,:,i), data_avg(:,:,i));
    data_reg = cat(3,data_reg,reg);
end

%%
data_reg_avg = squeeze(mean(reshape(data_reg,[529 796 100 11]),3));
figure;
for i = 1:9
    subplot(3,3,i)
    imagesc(data_reg_avg(:,:,i))
end

%%
data_reg_reg = [];

for i = 1:11;
    if i == 1
        data_reg_reg = data_reg_avg(:,:,1);
    else
        [out, reg] = stackRegister(data_reg_avg(:,:,i), data_reg_avg(:,:,i-1));
        data_reg_reg = cat(3,data_reg_reg,reg);
    end
end

figure;
for i = 1:9
    subplot(3,3,i)
    imagesc(data_reg_reg(:,:,i))
end


%%
mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
CD = fullfile(fnout, [date '_' mouse '\' date '_' mouse '_' 'runs-' ImgFolder(irun,:)]);
cd(CD);
writetiff(data_reg_reg, '5_7x_postlight_32min.tiff', 'double');


