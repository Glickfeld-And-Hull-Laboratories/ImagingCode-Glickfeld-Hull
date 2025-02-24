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


