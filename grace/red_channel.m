%% 
clear all;clear global
%% 
date = '210513';
ImgFolder = strvcat('003');
mouse = 'i1800';
nrun = size(ImgFolder,1);
run_str = catRunName(ImgFolder, nrun);
gl_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\2P_Imaging';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P';

mask = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']));
shifts = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']));
dfof = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimActFOV.mat']));
fitGeo = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_multiday_alignment.mat']));

trans = fitGeo.fitGeoTAf;
fov_avg = shifts.data_reg_avg;
% fov_avg = imwarp(fov_avg1,trans, 'OutputView', imref2d(size(fov_avg1))); 
masks = mask.mask_cell;
cell_stats = regionprops(masks);
dfof_max1 = dfof.data_dfof_max;
dfof_max = imwarp(dfof_max1,trans, 'OutputView', imref2d(size(dfof_max1))); 
df_max = dfof.data_df_max;

%% red channel
day = 1;
irun = 1;
ImgFolder = strvcat('001');
imgMatFile = [ImgFolder '_000_000.mat'];
run = catRunName(ImgFolder, nrun);
CD = fullfile(gl_fn, [mouse '\' date '_' mouse '\' ImgFolder(irun,:)]);
cd(CD);
load(imgMatFile);
% fprintf(['Reading run ' num2str(irun) '- ' num2str(info.config.frames) ' frames \r\n'])
data_temp = sbxread(imgMatFile(1,1:11),0,info.config.frames);

data_rg = squeeze(data_temp(1,:,:,:));
data_rr = squeeze(data_temp(2,:,:,:));

[out data_g_reg] = stackRegister(data_rg,fov_avg);
[out2 data_r_reg] = stackRegister_MA(data_rr,[],[],out);
red = mean(data_r_reg,3);
redChImg = imwarp(red,trans, 'OutputView', imref2d(size(red))); 
green = mean(data_g_reg,3);
greenChImg = imwarp(green,trans, 'OutputView', imref2d(size(green))); 
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_red.mat']),'redChImg','greenChImg');
clear data_temp

mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run]));
sz = size(data_rg);
figure;
colormap gray; imagesc(redChImg); movegui('center')
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run], ['902_red.pdf']),'-dpdf','-bestfit')
figure;
colormap gray; imagesc(greenChImg); movegui('center')
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run], ['902_green.pdf']),'-dpdf','-bestfit')

sz = size(data_rr);
rgb = zeros(sz(1),sz(2),3);
rgb(:,:,1) = redChImg./max(redChImg(:));
rgb(:,:,2) = greenChImg./max(greenChImg(:));
figure; image(rgb); movegui('center')
hold on
bound = cell2mat(bwboundaries(masks(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
title('RGB from 1000 frames')
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run], [mouse ' red_green2.pdf']),'-dpdf','-bestfit')

figure;
subplot(1,3,1)
imagesc(redChImg);title('red channel');axis image
subplot(1,3,2)
imagesc(greenChImg);title('green channel');axis image
subplot(1,3,3)
image(rgb);title('rgb');axis image
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run], [mouse ' red_dfof.pdf']),'-dpdf','-bestfit')


%% figures
data_avg = mean(data_rg(:,:,400:450),3);
[out data_reg] = stackRegister(data_rg,data_avg);
data_reg_avg = mean(data_reg,3);

data_avg = mean(data_rr(:,:,400:450),3);
[out data_reg] = stackRegister(data_rr,data_avg);
data_reg_avg2 = mean(data_reg,3);

ImgFolder = strvcat('003');
imgMatFile = [ImgFolder '_000_000.mat'];
CD = fullfile(gl_fn, [mouse '\' date '_' mouse '\' ImgFolder(irun,:)]);
cd(CD);
load(imgMatFile);
% fprintf(['Reading run ' num2str(irun) '- ' num2str(info.config.frames) ' frames \r\n'])
data_temp = sbxread(imgMatFile(1,1:11),0,info.config.frames);
data_rg = squeeze(data_temp(1,:,:,:));
data_rr = squeeze(data_temp(2,:,:,:));
clear data_temp

data_avg = mean(data_rg(:,:,400:450),3);
[out data_reg] = stackRegister(data_rg,data_avg);
data_reg_avg3 = mean(data_reg,3);

data_avg = mean(data_rr(:,:,400:450),3);
[out data_reg] = stackRegister(data_rr,data_avg);
data_reg_avg4 = mean(data_reg,3);

figure;
suptitle([mouse])
subplot(3,2,1)
imagesc(data_reg_avg);title('green 920');axis image
subplot(3,2,2)
imagesc(data_reg_avg2);title('red 920');axis image
subplot(3,2,3)
imagesc(data_reg_avg3);title('green 1040');axis image
subplot(3,2,4)
imagesc(data_reg_avg4);title('red 1040');axis image
subplot(3,2,5)
imagesc(fov_avg);title('green 920 from expt');axis image
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['FOV_compare.pdf']),'-dpdf','-bestfit')


%% cell map
redChImg1 = imwarp(redChImg,trans, 'OutputView', imref2d(size(redChImg))); 
figure;
start = 1;
for iCell = 1:10
    width = 24; height = 24;
    xCenter = round(cell_stats(iCell).Centroid(2));
    yCenter = round(cell_stats(iCell).Centroid(1));
    xLeft(iCell) = (xCenter - width/2);
    yBottom(iCell) = (yCenter - height/2);
    cell_mask = masks(xLeft(iCell):(xLeft(iCell)+width),yBottom(iCell):(height+yBottom(iCell)));
    if xLeft(iCell) > 12 && xLeft(iCell) < 488 && yBottom(iCell) > 12 && yBottom(iCell) < 772
    subplot(10,6,start);
    x = redChImg1(xLeft(iCell):(xLeft(iCell)+width),yBottom(iCell):(height+yBottom(iCell)));
    imagesc(x)
    pos = get(gca, 'Position');
    pos(1) = 0.025;
    pos(3) = 0.05;
    set(gca, 'Position', pos)
    axis off
    axis square
    hold on
    bound = cell2mat(bwboundaries(cell_mask(:,:,1)));
    plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
    title('red')
    subplot(10,6,start+1);
    y = dfof_max(xLeft(iCell):(xLeft(iCell)+width),yBottom(iCell):(height+yBottom(iCell)));
    imagesc(y)
    pos = get(gca, 'Position');
    pos(1) = 0.075;
    pos(3) = 0.05;
    set(gca, 'Position', pos)
    axis off
    axis square
    hold on
    bound = cell2mat(bwboundaries(cell_mask(:,:,1)));
    plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
    title('dfof')
    subplot(10,6,start+2);
    z = rgb(xLeft(iCell):(xLeft(iCell)+width),yBottom(iCell):(height+yBottom(iCell)),:);
    imagesc(z)
    pos = get(gca, 'Position');
    pos(1) = 0.125;
    pos(3) = 0.05;
    set(gca, 'Position', pos)
    axis off
    axis square
    hold on
    bound = cell2mat(bwboundaries(cell_mask(:,:,1)));
    plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
    title('rgb')
    end
    
    xCenter = round(cell_stats(iCell+10).Centroid(2));
    yCenter = round(cell_stats(iCell+10).Centroid(1));
    xLeft(iCell+10) = (xCenter - width/2);
    yBottom(iCell+10) = (yCenter - height/2);
    cell_mask = masks(xLeft(iCell+10):(xLeft(iCell+10)+width),yBottom(iCell+10):(height+yBottom(iCell+10)));
    if xLeft(iCell+10) > 12 && xLeft(iCell+10) < 488 && yBottom(iCell+10) > 12 && yBottom(iCell+10) < 772
    subplot(10,6,start+3);
    x = redChImg1(xLeft(iCell+10):(xLeft(iCell+10)+width),yBottom(iCell+10):(height+yBottom(iCell+10)));
    imagesc(x)
    pos = get(gca, 'Position');
    pos(1) = .175;
    pos(3) = .05;
    set(gca, 'Position', pos)
    axis off
    axis square
    hold on
    bound = cell2mat(bwboundaries(cell_mask(:,:,1)));
    plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
    title('red')
    subplot(10,6,start+4);
    y = dfof_max(xLeft(iCell+10):(xLeft(iCell+10)+width),yBottom(iCell+10):(height+yBottom(iCell+10)));
    imagesc(y)
    pos = get(gca, 'Position');
    pos(1) = 0.225;
    pos(3) = 0.05;
    set(gca, 'Position', pos)
    axis off
    axis square
    hold on
    bound = cell2mat(bwboundaries(cell_mask(:,:,1)));
    plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
    title('dfof')
    subplot(10,6,start+5);
    z = rgb(xLeft(iCell+10):(xLeft(iCell+10)+width),yBottom(iCell+10):(height+yBottom(iCell+10)),:);
    imagesc(z)
    pos = get(gca, 'Position');
    pos(1) = 0.275;
    pos(3) = 0.05;
    set(gca, 'Position', pos)
    axis off
    axis square
    hold on
    bound = cell2mat(bwboundaries(cell_mask(:,:,1)));
    plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
    title('rgb')
    end
    
    start = start+6;
end
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [mouse '_rgb_stamps.pdf']),'-dpdf','-bestfit')

