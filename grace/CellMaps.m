clear all
clear global
%% 
mouse = 'i1313';
day2 = '200120';
day3 = '200201';
ImgFolder = strvcat('002');
ImgFolder2 = strvcat('002');
ref_date = '200118';
ref_run = strvcat('002');
nrun = size(ImgFolder,1);
nrun2 = size(ImgFolder2,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);
run_str2 = catRunName(ImgFolder2, nrun2);
run_str3 = catRunName(ImgFolder, nrun);
ref_str = catRunName(ref_run, size(ref_run,1));
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P';

% oriTuning_D1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_oriTuningAndFits.mat']));
% oriTuning_D2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_oriTuningAndFits.mat']));
% oriTuning_D3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_oriTuningAndFits.mat']));

TCs_D1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_TCs.mat']));
TCs_D2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_TCs.mat']));
TCs_D3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_TCs.mat']));

% [maxResp_D1 prefOri_D1] = max(squeeze(oriTuning_D1.vonMisesFitAllCellsAllBoots(:,1,:)),[],1);
% [maxResp_D2 prefOri_D2] = max(squeeze(oriTuning_D2.vonMisesFitAllCellsAllBoots(:,1,:)),[],1);
% [maxResp_D3 prefOri_D3] = max(squeeze(oriTuning_D3.vonMisesFitAllCellsAllBoots(:,1,:)),[],1);

% day 1
load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_input.mat']));
tGratingDir1 = celleqel2mat_padded(input.tGratingDirectionDeg);
dirs1 = unique(tGratingDir1);
nDir1 = length(dirs1);
nOn1 = input.nScansOn;
nOff1 = input.nScansOff;
nFrames1 = nOn1+nOff1;
npSub_tc1 = TCs_D1.npSub_tc;
nCells1 = size(npSub_tc1,2);
nTrials1 = size(tGratingDir1,2);

trial_tc1 = reshape(npSub_tc1,[nFrames1 nTrials1 nCells1]);
trial_f1 = mean(trial_tc1(nOff1/2:nOff1,:,:),1);
trial_dfof1 = (trial_tc1-trial_f1)./trial_f1;

% day 2
load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_input.mat']));
tGratingDir2 = celleqel2mat_padded(input.tGratingDirectionDeg);
dirs2 = unique(tGratingDir1);
nDir2 = length(dirs2);
nOn2 = input.nScansOn;
nOff2 = input.nScansOff;
nFrames2 = nOn2+nOff2;
npSub_tc2 = TCs_D2.npSub_tc;
nCells2 = size(npSub_tc2,2);
nTrials2 = size(tGratingDir2,2);

trial_tc2 = reshape(npSub_tc2,[nFrames2 nTrials2 nCells2]);
trial_f2 = mean(trial_tc2(nOff2/2:nOff2,:,:),1);
trial_dfof2 = (trial_tc2-trial_f2)./trial_f2;

% day 3
load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_input.mat']));
tGratingDir3 = celleqel2mat_padded(input.tGratingDirectionDeg);
dirs3 = unique(tGratingDir1);
nDir3 = length(dirs3);
nOn3 = input.nScansOn;
nOff3 = input.nScansOff;
nFrames3 = nOn3+nOff3;
npSub_tc3 = TCs_D3.npSub_tc;
nCells3 = size(npSub_tc3,2);
nTrials3 = size(tGratingDir3,2);

trial_tc3 = reshape(npSub_tc3,[nFrames3 nTrials3 nCells3]);
trial_f3 = mean(trial_tc3(nOff3/2:nOff3,:,:),1);
trial_dfof3 = (trial_tc3-trial_f3)./trial_f3;


%% Cell Maps
% loading data
maskD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_mask_cell.mat']));
mask_cell = maskD1.mask_cell;
dfof = maskD1.data_dfof_max;
reg_shifts = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_reg_shifts.mat']));
reg1 = reg_shifts.data_reg_avg;
reg1(find(reg1>7000)) = 0;
reg1 = (reg1./max(max(abs(reg1))));
cell_list = intersect(1:nCells1, unique(mask_cell));
cell_stats = regionprops(mask_cell);

maskD2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_mask_cell.mat']));
mask_cell2 = maskD2.mask_cell;
reg_shifts2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_transform.mat']));
dfof2 = reg_shifts2.r2rFGTA_dfof;
reg2 = reg_shifts2.r2rFGTA;
reg2(find(reg2>7000)) = 0;
reg2 = (reg2./max(max(abs(reg2))));
cell_list2 = intersect(1:nCells1, unique(mask_cell2));
cell_stats2 = regionprops(mask_cell2);
og_reg2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str], [day2 '_' mouse '_' run_str '_reg_shifts.mat']));
data_reg_avg2 = og_reg2.data_reg_avg;
stimAct2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str], [day2 '_' mouse '_' run_str '_stimActFOV.mat']));
dfof_max2 = stimAct2.data_dfof_max;

maskD3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_mask_cell.mat']));
mask_cell3 = maskD3.mask_cell;
reg_shifts3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_transform.mat']));
dfof3 = reg_shifts3.r2rFGTA_dfof;
reg3 = reg_shifts3.r2rFGTA;
reg3(find(reg3>7000)) = 0;
reg3 = (reg3./max(max(abs(reg3))));
cell_list3 = intersect(1:nCells1, unique(mask_cell3));
cell_stats3 = regionprops(mask_cell3);
og_reg3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_reg_shifts.mat']));
data_reg_avg3 = og_reg3.data_reg_avg;


%% Pixel Correlation

load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_pixel.mat']));
pix_3hz1 = pix;
pix_3hz1(isnan(pix_3hz1))=0;
load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_pixel.mat']));
pix2 = pix;
pix_fgta2 = pix_fgta;
load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_pixel.mat']));
pix3 = pix;
pix_fgta3 = pix_fgta;

r = triu2vec(corrcoef(pix_3hz1(:),pix_fgta2(:)));
til_str1 = num2str(chop(r,3));
r = triu2vec(corrcoef(pix_3hz1(:),pix_fgtn2(:)));
til_str2 = num2str(chop(r,3));
r = triu2vec(corrcoef(pix_3hz1(:),pix_imrt2(:)));
til_str3 = num2str(chop(r,3));

figure; 
subplot(2,2,1);imagesc(pix_3hz1);axis off;title(['ref date pixel corr'])
subplot(2,2,2);imagesc(pix_fgta2);axis off;title(['fitgeo affine pixel corr ' til_str1])
subplot(2,2,3);imagesc(pix_fgtn2);axis off;title(['fitgeo nonreflective pixel corr ' til_str2])
subplot(2,2,4);imagesc(pix_imrt2);axis off;title(['imregtform pixel corr ' til_str3])
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['CellMaps'], [ref_date '_' mouse '_registration_corr.pdf']),'-dpdf', '-bestfit')

% lineplot
figure;
x = [1 2 3];
xticks([1 2 3]);
y = [r2 r3 r1];
plot(x,y)
cellnames = {'fitgeo aff';'fitgeo nr';'imregt'};
ylabel('correlation coefficient')
set(gca,'xticklabel',cellnames)
hold on
err = [std(r2,[],1)./sqrt(nCells1) std(r3,[],1)./sqrt(nCells1) std(r1,[],1)./sqrt(nCells1)];
y_avg = mean(y,1);
errorbar(x,y_avg,err,'-ok','LineWidth',2);
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['CellMaps'], [ref_date '_' mouse '_reg_lineplot.pdf']),'-dpdf', '-bestfit')


% re-registering squares
height = 30; width = 30;
p2 = zeros(nCells1,1);
p3 = zeros(nCells1,1);
r2 = zeros(nCells1,1);
r3 = zeros(nCells1,1);
for iCell = 1:nCells1
    xCenter = round(cell_stats(iCell).Centroid(2));
    yCenter = round(cell_stats(iCell).Centroid(1));
    xLeft = (xCenter - width/2);
    yBottom = (yCenter - height/2);
 if xLeft > 1 && xLeft < 482 && yBottom > 1 && yBottom < 766
    a = pix_3hz1(xLeft:(xLeft+width),yBottom:(height+yBottom));
    a1 = reg1(xLeft:(xLeft+width),yBottom:(height+yBottom));
    x = pix_fgta2(xLeft:(xLeft+width),yBottom:(height+yBottom));
    [reg shift] = shift_opt(x,a,4);
    regpix = reg;
    y = reg2(xLeft:(xLeft+width),yBottom:(height+yBottom));
    [reg shift] = shift_opt(y,a1,4);
    regr2r = reg;
    z = pix_fgta3(xLeft:(xLeft+width),yBottom:(height+yBottom));
    [reg shift] = shift_opt(z,a,4);
    regpix3 = reg;
    h = reg3(xLeft:(xLeft+width),yBottom:(height+yBottom));
    [reg shift] = shift_opt(h,a1,4);
    regr2r3 = reg;
    p2(iCell) = triu2vec(corrcoef(a(:),regpix(:)));
    p3(iCell) = triu2vec(corrcoef(a(:),regpix3(:)));
    r2(iCell) = triu2vec(corrcoef(a1(:),regr2r(:)));
    r3(iCell) = triu2vec(corrcoef(a1(:),regr2r3(:)));
 end
end

% pixel corr scatter plot
figure;
scatter(r2,r3)
xlabel('reg correlation D1 and D2')
ylabel('reg correlation D1 and D3')
set(gca,'FontSize',16)
refline(1,0)
R = corrcoef(r2,r3);
disp(R(1,2));
str = ['    r = ',num2str(R(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
axis square
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['CellMaps'], [ref_date '_' mouse '_reg_scat.pdf']),'-dpdf', '-bestfit')

figure;
scatter(p2,p3)
xlabel('pixel correlation D1 and D2')
ylabel('pixel correlation D1 and D3')
set(gca,'FontSize',16)
refline(1,0)
R = corrcoef(p2,p3);
disp(R(1,2));
str = ['    r = ',num2str(R(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
axis square
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['CellMaps'], [ref_date '_' mouse '_pixel_scat.pdf']),'-dpdf', '-bestfit')

% comparing squares post fine registartion (day 1 vs fitgeo affine)

r1_2 = find(r2>0.1 & r2<0.2);
r2_2 = find(r2>0.2 & r2<0.3);
r3_2 = find(r2>0.3 & r2<0.4);
r4_2 = find(r2>0.4 & r2<0.5);
r5_2 = find(r2>0.5 & r2<0.6);
r6_2 = find(r2>0.6 & r2<0.7);
r7_2 = find(r2>0.7 & r2<0.8);
r8_2 = find(r2>0.8 & r2<0.9);
r9_2 = find(r2>0.9 & r2<1);

r1_3 = find(r3>0.1 & r3<0.2);
r2_3 = find(r3>0.2 & r3<0.3);
r3_3 = find(r3>0.3 & r3<0.4);
r4_3 = find(r3>0.4 & r3<0.5);
r5_3 = find(r3>0.5 & r3<0.6);
r6_3 = find(r3>0.6 & r3<0.7);
r7_3 = find(r3>0.7 & r3<0.8);
r8_3 = find(r3>0.8 & r3<0.9);
r9_3 = find(r3>0.9 & r3<1);

p1_2 = find(p2>0.1 & p2<0.2);
p2_2 = find(p2>0.2 & p2<0.3);
p3_2 = find(p2>0.3 & p2<0.4);
p4_2 = find(p2>0.4 & p2<0.5);
p5_2 = find(p2>0.5 & p2<0.6);
p6_2 = find(p2>0.6 & p2<0.7);
p7_2 = find(p2>0.7 & p2<0.8);
p8_2 = find(p2>0.8 & p2<0.9);
p9_2 = find(p2>0.9 & p2<1);

p1_3 = find(p3>0.1 & p3<0.2);
p2_3 = find(p3>0.2 & p3<0.3);
p3_3 = find(p3>0.3 & p3<0.4);
p4_3 = find(p3>0.4 & p3<0.5);
p5_3 = find(p3>0.5 & p3<0.6);
p6_3 = find(p3>0.6 & p3<0.7);
p7_3 = find(p3>0.7 & p3<0.8);
p8_3 = find(p3>0.8 & p3<0.9);
p9_3 = find(p3>0.9 & p3<1);

badP = [p1_2;p2_2;p3_2;p1_2;p2_3;p3_3];
badR = [r1_2;r2_2;r3_2;r1_3;r2_3;r3_3];
badPR = intersect(badR,badP);

decentP = [p4_2;p5_2;p4_3;p5_3];
decentR = [r4_2;r5_2;r6_2;r4_3;r5_3;r6_3];
decentPR = intersect(decentP,decentR);

goodP = [p6_2;p7_2;p8_2;p9_2;p6_3;p7_3;p8_3;p9_3];
goodR = [r7_2;r8_2;r9_2;r7_3;r8_3;r9_3];
goodPR = intersect(goodP,goodR);
%% 
start = 1;
height = 30; width = 30;
figure;
for iC = 2:5
iCell = goodPR(iC);
xCenter = round(cell_stats(iCell).Centroid(2));
yCenter = round(cell_stats(iCell).Centroid(1));
xLeft = (xCenter - width/2);
yBottom = (yCenter - height/2);
cell_mask = mask_cell(xLeft:(xLeft+width),yBottom:(height+yBottom));
x = reg1(xLeft:(xLeft+width),yBottom:(height+yBottom));
subplot(4,6,start)
imagesc(x)
pos = get(gca, 'Position');
pos(1) = 0.1;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis square
axis off
hold on
    bound = cell2mat(bwboundaries(cell_mask(:,:,1)));
    plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
y = reg2(xLeft:(xLeft+width),yBottom:(height+yBottom));
[reg shift] = shift_opt(y,x,4);
subplot(4,6,start+1)
imagesc(reg)
pos = get(gca, 'Position');
pos(1) = 0.16;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis off
axis square
r = triu2vec(corrcoef(x(:),reg(:)));
til_str = num2str(chop(r,2));
title(til_str,'FontSize',8);
hold on
    bound = cell2mat(bwboundaries(cell_mask(:,:,1)));
    plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
z = reg3(xLeft:(xLeft+width),yBottom:(height+yBottom));
[hi shift] = shift_opt(z,x,4);
subplot(4,6,start+2)
imagesc(hi)
pos = get(gca, 'Position');
pos(1) = 0.22;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis square
axis off
r = triu2vec(corrcoef(x(:),hi(:)));
til_str = num2str(chop(r,2));
title(til_str,'FontSize',8);
hold on
    bound = cell2mat(bwboundaries(cell_mask(:,:,1)));
    plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);

xCenter = round(cell_stats(iCell).Centroid(2));
yCenter = round(cell_stats(iCell).Centroid(1));
xLeft = (xCenter - width/2);
yBottom = (yCenter - height/2);
x = pix_3hz1(xLeft:(xLeft+width),yBottom:(height+yBottom));
subplot(4,6,start+3)
imagesc(x)
pos = get(gca, 'Position');
pos(1) = 0.3;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis square
axis off
hold on
    bound = cell2mat(bwboundaries(cell_mask(:,:,1)));
    plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
y = pix_fgta2(xLeft:(xLeft+width),yBottom:(height+yBottom));
[hey shift] = shift_opt(y,x,4);
subplot(4,6,start+4)
imagesc(hey)
pos = get(gca, 'Position');
pos(1) = 0.36;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis off
axis square
r = triu2vec(corrcoef(x(:),hey(:)));
til_str = num2str(chop(r,2));
title(til_str,'FontSize',8);
hold on
    bound = cell2mat(bwboundaries(cell_mask(:,:,1)));
    plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
z = pix_fgta3(xLeft:(xLeft+width),yBottom:(height+yBottom));
[mad shift] = shift_opt(z,x,4);
subplot(4,6,start+5)
imagesc(mad)
pos = get(gca, 'Position');
pos(1) = 0.42;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis square
axis off
r = triu2vec(corrcoef(x(:),mad(:)));
til_str = num2str(chop(r,2));
title(til_str,'FontSize',8);
hold on
    bound = cell2mat(bwboundaries(cell_mask(:,:,1)));
    plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);

start = start+6;
end 
% print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['CellMaps'], [ref_date '_' mouse '_bad_cells.pdf']),'-dpdf', '-bestfit')

%% 
start = 1;
height = 30; width = 30;
figure;
% bad cells
for iC = 1:4
iCell = badPR(iC);
xCenter = round(cell_stats(iCell).Centroid(2));
yCenter = round(cell_stats(iCell).Centroid(1));
xLeft = (xCenter - width/2);
yBottom = (yCenter - height/2);
x = reg1(xLeft:(xLeft+width),yBottom:(height+yBottom));
subplot(4,8,start)
imagesc(x)
pos = get(gca, 'Position');
pos(1) = 0.1;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis square
axis off
y = data_reg_avg2(xLeft:(xLeft+width),yBottom:(height+yBottom));
subplot(4,8,start+1)
imagesc(y)
pos = get(gca, 'Position');
pos(1) = 0.16;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis off
axis square
r = triu2vec(corrcoef(x(:),y(:)));
til_str = num2str(chop(r,2));
title(til_str,'FontSize',6);
z = reg2(xLeft:(xLeft+width),yBottom:(height+yBottom));
subplot(4,8,start+2)
imagesc(z)
pos = get(gca, 'Position');
pos(1) = 0.22;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis square
axis off
r = triu2vec(corrcoef(x(:),z(:)));
til_str = num2str(chop(r,2));
title(til_str,'FontSize',6);
[reg shift] = shift_opt(z,x,4);
subplot(4,8,start+3)
imagesc(reg)
pos = get(gca, 'Position');
pos(1) = 0.28;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis off
axis square
r = triu2vec(corrcoef(x(:),reg(:)));
til_str = num2str(chop(r,2));
title(til_str,'FontSize',6);

xCenter = round(cell_stats(iCell).Centroid(2));
yCenter = round(cell_stats(iCell).Centroid(1));
xLeft = (xCenter - width/2);
yBottom = (yCenter - height/2);
x = pix_3hz1(xLeft:(xLeft+width),yBottom:(height+yBottom));
subplot(4,8,start+4)
imagesc(x)
pos = get(gca, 'Position');
pos(1) = 0.36;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis square
axis off
y = pix2(xLeft:(xLeft+width),yBottom:(height+yBottom));
subplot(4,8,start+5)
imagesc(y)
pos = get(gca, 'Position');
pos(1) = 0.42;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis off
axis square
r = triu2vec(corrcoef(x(:),y(:)));
til_str = num2str(chop(r,2));
title(til_str,'FontSize',6);
z = pix_fgta2(xLeft:(xLeft+width),yBottom:(height+yBottom));
subplot(4,8,start+6)
imagesc(z)
pos = get(gca, 'Position');
pos(1) = 0.48;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis square
axis off
r = triu2vec(corrcoef(x(:),z(:)));
til_str = num2str(chop(r,2));
title(til_str,'FontSize',6);
[reg shift] = shift_opt(z,x,4);
subplot(4,8,start+7)
imagesc(reg)
pos = get(gca, 'Position');
pos(1) = 0.54;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis off
axis square
r = triu2vec(corrcoef(x(:),reg(:)));
til_str = num2str(chop(r,2));
title(til_str,'FontSize',6);

start = start+8;
end 
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['CellMaps'], [ref_date '_' mouse '_bad_cells.pdf']),'-dpdf', '-bestfit')

start = 1;
height = 24; width = 24;
mat = [r7_3(1:3) r7_5(1:3)];
mat = reshape(mat,[],1);
figure;
suptitle('re-registered pixel corr and data reg avg, 0.7 corr coefs')
for iC = 1:3
iCell = mat(iC);
xCenter = round(cell_stats(iCell).Centroid(2));
yCenter = round(cell_stats(iCell).Centroid(1));
xLeft = (xCenter - width/2);
yBottom = (yCenter - height/2);
x = pix_3hz1(xLeft:(xLeft+width),yBottom:(height+yBottom));
subplot(8,8,start)
imagesc(x)
pos = get(gca, 'Position');
pos(1) = 0.1;
pos(3) = 0.07;
set(gca, 'Position', pos)
axis square
axis off
xCenter2 = round(cell_stats2(iCell).Centroid(2));
yCenter2 = round(cell_stats2(iCell).Centroid(1));
xLeft2 = (xCenter2 - width/2);
yBottom2 = (yCenter2 - height/2);
y = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
[reg shift] = shift_opt(y,x,3);
subplot(8,8,start+1)
imagesc(reg)
pos = get(gca, 'Position');
pos(1) = .18;
pos(3) = 0.07;
set(gca, 'Position', pos)
axis off
axis square
r = triu2vec(corrcoef(x(:),reg(:)));
til_str = num2str(chop(r,2));
title(til_str,'FontSize',6);
x2 = reg1(xLeft:(xLeft+width),yBottom:(height+yBottom));
subplot(8,8,start+2)
imagesc(x2)
pos = get(gca, 'Position');
pos(1) = 0.26;
pos(3) = 0.07;
set(gca, 'Position', pos)
axis square
axis off
y2 = r2r_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
[reg shift] = shift_opt(y2,x2,3);
subplot(8,8,start+3)
imagesc(reg)
pos = get(gca, 'Position');
pos(1) = 0.34;
pos(3) = 0.07;
set(gca, 'Position', pos)
axis off
axis square
r = triu2vec(corrcoef(x2(:),reg(:)));
til_str = num2str(chop(r,2));
title(til_str,'FontSize',6);

iCell = mat(iC+3);
xCenter = round(cell_stats(iCell).Centroid(2));
yCenter = round(cell_stats(iCell).Centroid(1));
xLeft = (xCenter - width/2);
yBottom = (yCenter - height/2);
x = pix_3hz1(xLeft:(xLeft+width),yBottom:(height+yBottom));
subplot(8,8,start+4)
imagesc(x)
pos = get(gca, 'Position');
pos(1) = 0.5;
pos(3) = 0.07;
set(gca, 'Position', pos)
axis square
axis off
xCenter2 = round(cell_stats2(iCell).Centroid(2));
yCenter2 = round(cell_stats2(iCell).Centroid(1));
xLeft2 = (xCenter2 - width/2);
yBottom2 = (yCenter2 - height/2);
y = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
[reg shift] = shift_opt(y,x,3);
subplot(8,8,start+5)
imagesc(reg)
pos = get(gca, 'Position');
pos(1) = .58;
pos(3) = 0.07;
set(gca, 'Position', pos)
axis off
axis square
r = triu2vec(corrcoef(x(:),reg(:)));
til_str = num2str(chop(r,2));
title(til_str,'FontSize',6);
x2 = reg1(xLeft:(xLeft+width),yBottom:(height+yBottom));
subplot(8,8,start+6)
imagesc(x2)
pos = get(gca, 'Position');
pos(1) = 0.66;
pos(3) = 0.07;
set(gca, 'Position', pos)
axis square
axis off
y2 = r2r_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
[reg shift] = shift_opt(y2,x2,3);
subplot(8,8,start+7)
imagesc(reg)
pos = get(gca, 'Position');
pos(1) = 0.74;
pos(3) = 0.07;
set(gca, 'Position', pos)
axis off
axis square
r = triu2vec(corrcoef(x2(:),reg(:)));
til_str = num2str(chop(r,2));
title(til_str,'FontSize',6);

start = start+8;
end 
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['CellMaps'], [ref_date '_' mouse '_7UP.pdf']),'-dpdf', '-bestfit')

figure;
y = [length(r1_3) length(r1_5); length(r2_3) length(r2_5); length(r3_3) length(r3_5); length(r4_3) length(r4_5); length(r5_3) length(r5_5); length(r6_3) length(r6_5); length(r7_3) length(r7_5); length(r8_3) length(r8_5);length(r9_3) length(r9_5)];
b = bar(y);
b(1).FaceColor = [0 0.8047 0.8164];
b(2).FaceColor = [0.5977 0.3984 0.7969];
ylabel('nCells')
cellnames = {'0.1';'0.2';'0.3';'0.4';'0.5';'0.6';'0.7';'0.8';'0.9'};
set(gca,'xticklabel',cellnames);
legend('pixel img','data reg avg','location','northwest')
legend boxoff
title('CorrCoef Distribution of Re-Registered Cells')
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['CellMaps'], [ref_date '_' mouse '_coef_dist.pdf']),'-dpdf', '-bestfit')

% theshold distribution curve-finding reg distribution
height = 30; width = 30;
mat = [p1_2;p2_2;p3_2;p4_2;p5_2];
mat = reshape(mat,[],1);
rr2r2 = zeros(length(mat),1);
for iC = 1:length(mat)
    iCell = mat(iC);
    xCenter = round(cell_stats(iCell).Centroid(2));
    yCenter = round(cell_stats(iCell).Centroid(1));
    xCenter2 = round(cell_stats2(iCell).Centroid(2));
    yCenter2 = round(cell_stats2(iCell).Centroid(1));
    xLeft = (xCenter - width/2);
    yBottom = (yCenter - height/2);
    xLeft2 = (xCenter2 - width/2);
    yBottom2 = (yCenter2 - height/2);
 if xLeft > 1 && xLeft < 482 && yBottom > 1 && yBottom < 766
    a1 = reg1(xLeft:(xLeft+width),yBottom:(height+yBottom));
    y = reg2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    [reg shift] = shift_opt(y,a1,4);
    regr2r = reg;
    rr2r2(iCell) = triu2vec(corrcoef(a1(:),regr2r(:)));
 end
end
rr2r2(rr2r2==0)=[];

% finding pix distribution
height = 30; width = 30;
mat = [r1_2;r2_2;r3_2;r4_2;r5_2];
mat = reshape(mat,[],1);
rpix2 = zeros(length(mat),1);
for iC = 1:length(mat)
    iCell = mat(iC);
    xCenter = round(cell_stats(iCell).Centroid(2));
    yCenter = round(cell_stats(iCell).Centroid(1));
    xCenter2 = round(cell_stats2(iCell).Centroid(2));
    yCenter2 = round(cell_stats2(iCell).Centroid(1));
    xLeft = (xCenter - width/2);
    yBottom = (yCenter - height/2);
    xLeft2 = (xCenter2 - width/2);
    yBottom2 = (yCenter2 - height/2);
 if xLeft > 1 && xLeft < 482 && yBottom > 1 && yBottom < 766
    a = pix_3hz1(xLeft:(xLeft+width),yBottom:(height+yBottom));
    x = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    [reg shift] = shift_opt(x,a,4);
    regpix = reg;
    rpix2(iCell) = triu2vec(corrcoef(a(:),regpix(:)));
 end
end
rpix2(rpix2==0)=[];

figure;
hist = histogram(rr2r2,11,'BinLimits',[-.1,1],'facealpha',0.5);
hold on
hist2 = histogram(rpix2,11,'BinLimits',[-.1,1],'facealpha',0.5);
legend('data reg avg','pixel img','location','northwest')
legend boxoff
ylabel('nCells')
xlabel('correlation coefficient')
set(gca,'FontSize',16)
hist.FaceColor = [1 0 0];
hist2.FaceColor = [0.2539 0.4102 0.8789];
title('CorrCoef Dirstribution of Cells below 0.6 Thresh of Opposite Condition') 
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['CellMaps'], [ref_date '_' mouse '_below6thresh.pdf']),'-dpdf', '-bestfit')

r79 = find(p2>0.7 & r2>0.9);
r78 = find(p2>0.7 & r2>0.8);
r68 = find(p2>0.6 & r2>0.8);
r78rescue = [find(r3<0.7 & r2>0.8);find(r3>0.7 & r2<0.8)];

figure;
y = [length(r79) length(r78) length(r68)];
b = bar(y);
b.FaceColor = [0.3320 0.4180 0.1836];
ylabel('nCells')
cellnames = {'(p>0.7, r>0.9)      ';'(p>0.7, r>0.8)';'      (p>0.6, r>0.8)'};
set(gca,'xticklabel',cellnames,'FontSize',14)
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['CellMaps'], [ref_date '_' mouse '_nCells_per_&threshold.pdf']),'-dpdf', '-bestfit')

start = 1;
height = 24; width = 24;
figure;
suptitle('cells rescued by pixel corr>0.7 or data reg>0.8')
for iC = 1:10
iCell = r78rescue(iC);
xCenter = round(cell_stats(iCell).Centroid(2));
yCenter = round(cell_stats(iCell).Centroid(1));
xLeft = (xCenter - width/2);
yBottom = (yCenter - height/2);
x = pix_3hz1(xLeft:(xLeft+width),yBottom:(height+yBottom));
subplot(12,4,start)
imagesc(x)
pos = get(gca, 'Position');
pos(1) = 0.1;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis square
axis off
xCenter2 = round(cell_stats2(iCell).Centroid(2));
yCenter2 = round(cell_stats2(iCell).Centroid(1));
xLeft2 = (xCenter2 - width/2);
yBottom2 = (yCenter2 - height/2);
y = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
[reg shift] = shift_opt(y,x,4);
subplot(12,4,start+1)
imagesc(reg)
pos = get(gca, 'Position');
pos(1) = .15;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis off
axis square
r = triu2vec(corrcoef(x(:),reg(:)));
til_str = num2str(chop(r,2));
title(til_str,'FontSize',5);

xCenter = round(cell_stats(iCell).Centroid(2));
yCenter = round(cell_stats(iCell).Centroid(1));
xLeft = (xCenter - width/2);
yBottom = (yCenter - height/2);
x = reg1(xLeft:(xLeft+width),yBottom:(height+yBottom));
subplot(12,4,start+2)
imagesc(x)
pos = get(gca, 'Position');
pos(1) = 0.2;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis square
axis off
xCenter2 = round(cell_stats2(iCell).Centroid(2));
yCenter2 = round(cell_stats2(iCell).Centroid(1));
xLeft2 = (xCenter2 - width/2);
yBottom2 = (yCenter2 - height/2);
y = r2r_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
[reg shift] = shift_opt(y,x,4);
subplot(12,4,start+3)
imagesc(reg)
pos = get(gca, 'Position');
pos(1) = .25;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis off
axis square
r = triu2vec(corrcoef(x(:),reg(:)));
til_str = num2str(chop(r,2));
title(til_str,'FontSize',5);

start = start+4;
end
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['CellMaps'], [ref_date '_' mouse '_r78rescue.pdf']),'-dpdf', '-bestfit')
