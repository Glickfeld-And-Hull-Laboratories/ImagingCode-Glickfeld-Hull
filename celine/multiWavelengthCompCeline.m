pn = 'Z:\home\ACh\Data\2p_data';
mouse = 'i3323';
date = '250319';

ImgType = 'flpMCh';
wl = [800 980 1040];
Runs = [0 1 2];


for i = 1:length(wl)
    ImgFolder = ['00' num2str(Runs(i))];
    CD = fullfile(pn, mouse, date, ImgFolder);
    cd(CD);
    imgMatFile = [ImgFolder '_000_000.mat'];
    load(imgMatFile);
    totframes = info.config.frames;
    data = sbxread([ImgFolder '_000_000'],0,totframes);
    data_r = squeeze(data(2,:,:,:));
    data_g = squeeze(data(1,:,:,:));
    if i == 1
        data_r_avg = zeros(size(data_r,1),size(data_r,2),length(wl));
        data_g_avg = zeros(size(data_r,1),size(data_r,2),length(wl));
        data_avg = mean(data_r(:,:,1:100),3);
        [out data_r_reg] = stackRegister(data_r,data_avg);
        data_r_orig_avg = mean(data_r_reg,3);
        data_r_avg(:,:,i) = data_r_orig_avg;
    else
        [out data_r_reg] = stackRegister(data_r,data_r_orig_avg);
        data_r_avg(:,:,i) = mean(data_r_reg,3);
    end
    [out2 data_g_reg] = stackRegister_MA(data_g,[],[],out);
    data_g_avg(:,:,i) = mean(data_g_reg,3);
end

fn_out = 'Z:\home\ACh\Analysis\2p_analysis';
mkdir(fullfile(fn_out,mouse,date));
% writetiff(data_r_avg, fullfile(fn_out,mouse,date,[ImgType '_Red.tif']))
% writetiff(data_g_avg, fullfile(fn_out,mouse,date,[ImgType '_Green.tif']))

sz = size(data_r_avg);
mask_data = data_r_avg;
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
for iStim = 1:size(data_r_avg,3)    
    mask_data_temp = mask_data(:,:,iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0; %blacks out old cells
    bwout = imCellEditInteractiveLG(mask_data_temp); %selection GUI
    mask_all = mask_all+bwout; %adds new cells to old cells
    mask_exp = imCellBuffer(mask_all,3)+mask_all; %creates buffer around cells to avoid fusing
    close all
end
mask_cell = bwlabel(mask_all); %turns logical into numbered cells
figure;
imagesc(mask_cell)

redCells = stackGetTimeCourses(data_r_avg, mask_cell);
redCells_norm = redCells./max(redCells,[],1);
nCells = size(redCells,2);
cormat = zeros(nCells,nCells);

for i = 1:nCells
    for j = 1:nCells
        cormat(i,j) = triu2vec(corrcoef(redCells(:,i),redCells(:,j)));
    end
end

figure; 
imagesc(redCells_norm')
set(gca, 'XTickLabels',wl)
xlabel('Wavelength')
ylabel('Cell #')
title([mouse ' ' date ' ' ImgType])
print(fullfile(fn_out,mouse,date,[ImgType '_RedCellHeatmap.pdf']),'-dpdf','-bestfit')
save(fullfile(fn_out,mouse,date,[ImgType '_RedCellHeatmap.pdf']),'mask_cell', 'data_r_avg','data_g_avg','redCells','redCells_norm')

figure
for i = 1:3
    subplot(2,3,i)
    imagesc(data_r_avg(:,:,i))
    colormap gray;
    title(num2str(wl(i)))
    if i==1
        ylabel('red')
    end

    subplot(2,3,i+3)
    imagesc(data_g_avg(:,:,i))
    colormap gray;
    if i==1
        ylabel('green')
    end
end
print(fullfile(fn_out,mouse,date,[ImgType '_comparison.pdf']),'-dpdf','-bestfit')