%pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Data\2P_images';
pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\ACh\Data\2p_data\';
mouse = 'i3320';
date = '241115';

ImgFolder ='003';

ImgType = 'GCaMP&flpMCh&flexTdt';
%ImgType = 'flxGCaMP&flpHT';
%wl = 800:40:1040;%
%wl = [920 800 1040];
wl = [800 900 1040];
%Run = 0:length(wl)-1;
Run = [3 2 1];
%RunType = mat2cell(wl, [1 length(wl)],[1]);

for i = 1:length(wl)
    ImgFolder = ['00' num2str(Run(i))];
    CD = fullfile(pn, mouse, date, ImgFolder);
    cd(CD);
    %imgMatFile = [ImgFolder '_000_00' num2str(Run(i)) '.mat'];
    imgMatFile = [ImgFolder '_000_000.mat'];
    load(imgMatFile);
    totframes = info.config.frames;
    %data = sbxread([ImgFolder '_000_00' num2str(Run(i))],0,totframes);
    data = sbxread([ImgFolder '_000_000'],0,totframes);
    data_r = squeeze(data(2,:,:,:));
    data_g = squeeze(data(1,:,:,:));
    if i == 1
        data_r_avg = zeros(size(data_r,1),size(data_r,2),length(wl));
        data_g_avg = zeros(size(data_r,1),size(data_r,2),length(wl));
        data_avg = mean(data_r(:,:,1:100),3);
        [out data_r_reg] = stackRegister(data_r,data_avg);
        data_r_orig_avg = mean(data_r_reg,3);
        data_r_avg(:,:,Run(i)) = data_r_orig_avg;
    else
        [out data_r_reg] = stackRegister(data_r,data_r_orig_avg);
        data_r_avg(:,:,Run(i)) = mean(data_r_reg,3);
    end
    if i == 3
        [out2 data_g_reg] = stackRegister(data_g, mean(data_g_reg,3));
        data_g_avg(:,:,Run(i)) = mean(data_g_reg,3);
        [out data_r_reg] = stackRegister_MA(data_r,[],[],out2);
        data_r_avg(:,:,Run(i)) = mean(data_r_reg,3);
    else
        [out2 data_g_reg] = stackRegister_MA(data_g,[],[],out);
        data_g_avg(:,:,Run(i)) = mean(data_g_reg,3);
    end
end

fn_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P';
mkdir(fullfile(fn_out,[date '_' mouse],[date '_' mouse '_' ImgFolder]));
writetiff(data_r_avg, fullfile(fn_out,[date '_' mouse],[date '_' mouse '_' ImgFolder],[date '_' mouse '_' ImgType '_Red.tif']))
writetiff(data_g_avg, fullfile(fn_out,[date '_' mouse],[date '_' mouse '_' ImgFolder],[date '_' mouse '_' ImgType '_Green.tif']))

for i = 1:3
subplot(2,3,i)
imagesc(data_g_avg(:,:,i))
title(num2str(wl(i)))
colormap gray
subplot(2,3,i+3)
imagesc(data_r_avg(:,:,i))
colormap gray
end
sgtitle([mouse ' ' date])
print(fullfile(fn_out,[date '_' mouse],[date '_' mouse '_' ImgFolder],[date '_' mouse '_' ImgType '_FOVbyWL.pdf']),'-dpdf','-bestfit')

%%
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
mask_np = imCellNeuropil(mask_cell, 3, 5);
nCells = size(redCells,2);
redCells_sub = zeros(size(redCells));
for i = 1:nCells
    for ii = 1:3
        redCells_sub(ii,i) = redCells(ii,i)-stackGetTimeCourses(data_r_avg(:,:,ii),mask_np(:,:,i));
    end
end

redCells_norm_sub = redCells_sub./max(redCells_sub,[],2);
figure;
scatter(redCells_norm_sub(1,:),redCells_norm_sub(3,:))
xlim([0 1]); ylim([0 1]); 
xlabel(num2str(wl(1)))
ylabel(num2str(wl(3)))
sgtitle([mouse ' ' date])
print(fullfile(fn_out,[date '_' mouse],[date '_' mouse '_' ImgFolder],[date '_' mouse '_' ImgType '_RedCellF.pdf']),'-dpdf','-bestfit')


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

figure; scatter(redCells(1,:),redCells(3,:))
xlim([0 600])
ylim([0 2000])
print(fullfile(fn_out,[date '_' mouse],[date '_' mouse '_' ImgFolder],[date '_' mouse '_' ImgType '_RedCellHeatmap.pdf']),'-dpdf','-bestfit')
save(fullfile(fn_out,[date '_' mouse],[date '_' mouse '_' ImgFolder],[date '_' mouse '_' ImgType '_RedCellF.mat']),'mask_cell', 'data_r_avg','data_g_avg','redCells','redCells_norm')