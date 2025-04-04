%make low contrast noise grating at three SFs
targetSize = [918 1174];
SFs = [0.00001 0.05 0.15 0.5];
nSF = length(SFs);
noise_mat = [];

outPath = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Software\MWorks\Images\SFNoise';
for i = 1:nSF
    noise_mat(:,:,i) = circular_linear_grating(targetSize(1),200,SFs(i),0,.1,.5,0);
end
noise_mat = noise_mat(1:targetSize(1),1:targetSize(2),:).*255;
noise_mat = noise_mat-mean(noise_mat(:));
close all

start = 1;
for i  = 1:nSF
    figure;
    imagesc(noise_mat(:,:,i))
    axis off
    colormap gray
    saveas(gcf,fullfile(outPath,[num2str(start) '.png']))
    start = start+1;
end

%add noise gratings to images
pathName = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Software\MWorks\Images\BaseImages';
ls = dir(pathName);
outPath = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Software\MWorks\Images\NatImgPlusNoise';
files = find(~cellfun(@isempty,strfind({ls.name}, '.png')));
nfiles = length(files);
totfiles = size(ls,1);
start = 1;
for i = totfiles-nfiles+1:totfiles
    x = imread(fullfile(pathName,ls(i).name));
    xplusnoise = zeros([size(x) nSF]);
    for ii = 1:nSF
        xplusnoise(:,:,ii) = noise_mat(:,:,ii)+double(x);
        xplusnoise(find(xplusnoise<0)) = 0;
        xplusnoise(find(xplusnoise>255)) = 255;
        figure;
        imagesc(xplusnoise(:,:,ii))
        axis off
        colormap gray
        saveas(gcf,fullfile(outPath,[num2str(start) '.png']))
        start = start+1;
    end
end


%%
%make low contrast noise grating at three SFs and 2 oris
targetSize = [918 1174];
SFs = [0.12 0.48];
nSF = length(SFs);
Oris = [0 30];
nOri = length(Oris);
noise_mat = [];

outPath = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Software\MWorks\Images\2SF2OriNoise';
if ~exist(outPath)
    mkdir(outPath)
end

for i = 1:nSF
    for ii = 1:nOri
        noise_mat(:,:,i,ii) = circular_linear_grating(targetSize(1),200,SFs(i),Oris(ii),.1,.5,0);
    end
end
noise_mat = noise_mat(1:targetSize(1),1:targetSize(2),:,:).*255;
noise_mat = noise_mat-mean(noise_mat(:));
close all

start = 1;
for i  = 1:nSF
    for ii = 1:nOri
        figure;
        imagesc(noise_mat(:,:,i,ii))
        axis off
        colormap gray
        saveas(gcf,fullfile(outPath,[num2str(start) '.png']))
        start = start+1;
    end
end

%add noise gratings to images
pathName = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Software\MWorks\Images\BaseImages';

outPath = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Software\MWorks\Images\NatImgPlusNoise2Ori';
if ~exist(outPath)
    mkdir(outPath)
end

%choose deer image
x = imread(fullfile(pathName,'2.png'));
start = 1;
xplusnoise = zeros([size(x) nSF nOri]);
for i = 1:nSF
    for ii = 1:nOri
        xplusnoise(:,:,i,ii) = noise_mat(:,:,i,ii)+double(x);
        xplusnoise(find(xplusnoise<0)) = 0;
        xplusnoise(find(xplusnoise>255)) = 255;
        figure;
        imagesc(xplusnoise(:,:,i,ii))
        axis off
        colormap gray
        saveas(gcf,fullfile(outPath,[num2str(start) '.png']))
        start = start+1;
    end
end

%choose flower image
x = imread(fullfile(pathName,'5.png'));
xplusnoise = zeros([size(x) nSF nOri]);
for i = 1:nSF
    for ii = 1:nOri
        xplusnoise(:,:,i,ii) = noise_mat(:,:,i,ii)+double(x);
        xplusnoise(find(xplusnoise<0)) = 0;
        xplusnoise(find(xplusnoise>255)) = 255;
        figure;
        imagesc(xplusnoise(:,:,i,ii))
        axis off
        colormap gray
        saveas(gcf,fullfile(outPath,[num2str(start) '.png']))
        start = start+1;
    end
end

% all ori noise
Oris = [0:30:150];
nOri = length(Oris);
noise_mat = [];

outPath = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Software\MWorks\Images\SF2Ori6Noise';
if ~exist(outPath)
    mkdir(outPath)
end

for i = 1:nSF
    for ii = 1:nOri
        noise_mat(:,:,i,ii) = circular_linear_grating(targetSize(1),200,SFs(i),Oris(ii),.1,.5,0);
    end
end
noise_mat = noise_mat(1:targetSize(1),1:targetSize(2),:,:).*255;
noise_mat = noise_mat-mean(noise_mat(:));
close all

start = 1;
for i  = 1:nSF
    for ii = 1:nOri
        figure;
        imagesc(noise_mat(:,:,i,ii))
        axis off
        colormap gray
        saveas(gcf,fullfile(outPath,[num2str(start) '.png']))
        start = start+1;
    end
end