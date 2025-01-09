mouse = 'i2178';
date = '241211';
run = '001';
fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Data\2P_images';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P';
CD = fullfile(fn, [mouse '\' date '\' run]); %change current dierectory;
cd(CD); %set CD
imgMatFile = [run '_000_001.mat']; %add the 0s to the imaging file
load(imgMatFile); %**load this file from CD
nframes = info.config.frames; %find nframes from info
data_temp = sbxread([run '_000_001'],0,nframes); %reads data from 0 to nframes from raw file
data = squeeze(data_temp);
figure; 
subplot(2,1,1)
imagesc(mean(data,3))
title('Green 920')
imgMatFile = [run '_000_002.mat']; %add the 0s to the imaging file
load(imgMatFile); %**load this file from CD
nframes = info.config.frames; %find nframes from info
data_temp = sbxread([run '_000_002'],0,nframes); %reads data from 0 to nframes from raw file
data = squeeze(data_temp(2,:,:,:));
subplot(2,1,2)
imagesc(mean(data,3))
title('Red 1040')
if ~exist(fullfile(fnout,[date '_' mouse],[date '_' mouse '_' run]))
    mkdir(fullfile(fnout,[date '_' mouse],[date '_' mouse '_' run]))
end
print(fullfile(fnout,[date '_' mouse],[date '_' mouse '_' run],[date '_' mouse '_' run '_FOVAvg.pdf']),'-dpdf','-fillpage')