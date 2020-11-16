clear all; close all; clear global
mouse = '1243';
expDate = '201020';
imgFolder = '001';
fName = [imgFolder '_000_000'];
pmt = 2;
obj = '16x';
zm = 2;

cd(fullfile('Z:\home\ashley\data\',mouse,'two-photon imaging',expDate,imgFolder))
fnout = fullfile('Z:\home\ashley\Analysis\',mouse,'two-photon imaging',expDate,imgFolder);
% fnout = fullfile('Z:\home\ashley\data\',mouse,'two-photon imaging',expDate,imgFolder);
mkdir(fnout)
%%
load([fName '.mat'])

nframes = 100;

data = sbxread(fName,0,nframes);

data = squeeze(data(pmt,:,:,:));

%%

% [~,data_reg] = stackRegister(data,mean(data(:,:,1:100),3));

%%
img = mean(data,3);

figure;colormap gray; imagesc(img)
writetiff(img,fullfile(fnout,'meanImg'))

%% code for scale bar
%  [~,~,sb_img_50um] = scalebarCalib(expDate,obj,img,zm);
% 
% figure;colormap gray; imagesc(sb_img_50um)
% writetiff(sb_img_50um,fullfile(fnout,'scalebar_50um'))