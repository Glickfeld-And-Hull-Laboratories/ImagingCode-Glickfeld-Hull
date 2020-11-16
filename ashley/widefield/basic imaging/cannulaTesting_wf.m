mouse = '1225';
exptDate = '200702';
frameRateHz = 2;
rc = behavConstsAV;
fnout = fullfile(rc.ashleyAnalysis,mouse,'widefield imaging',exptDate);
mkdir(fnout)
%% load data, save dff movie
exptName = 'texasred_movie_duringInfusion';
dataPath = 'Z:\home\ashley\data\1225\widefield imaging\1225_200702\cannulatest_infusion_texasred_movie_during_1';
data = readtiff(dataPath);

dff = quickDFF_wf(data);
writetiff(dff,fullfile(fnout,[exptName '_dff']));
writetiff(data,fullfile(fnout,exptName));

img_after = data(:,:,120);
%% black out injection site
[ypix,xpix,nfr] = size(dff);
testFrame = 311;
figure; imagesc(data(:,:,testFrame)); colormap gray

injsite_mask = imCellPolyEditInteractive(data(:,:,testFrame));

dff_coverInj = dff;
dff_coverInj(repmat(injsite_mask,[1,nfr]) > 0) = min(dff(:));

figure; imagesc(dff_coverInj(:,:,testFrame)); colormap gray
writetiff(dff_coverInj,fullfile(fnout,[exptName '_coverInj']));

%% load before and after images

dataPath = 'Z:\home\ashley\data\1225\widefield imaging\1225_200702\cannulatest_infusion_texasred_movie_before_1';
img_before = readtiff(dataPath);

%register before to after images
[movingpoints,fixedpoints] = cpselect(img_before(:,:,1),img_after,'Wait',true);
img_before_reg_crop = transformRegImgAndCropToOrigSize(...
    img_before,img_after,movingpoints,fixedpoints);
figure;imagesc(img_before_reg_crop(:,:,1));colormap gray

img_before_sample = img_before_reg_crop(:,:,1);
% dff_before = quickDFF_wf(img_before_reg_crop);

%% get ROI
after_label = imCellPolyEditInteractive(img_after);
after_mask = bwlabel(after_label);
%% plot before and after images of same area of the window, and mask
img_fig = figure;colormap gray
subplot 331
imagesc(img_before_sample)
figAxForm
title('Before Infusion')
subplot 332
imagesc(img_after)
figAxForm
title('After Infusion')
figure(img_fig)
subplot 333
imagesc(after_mask);
figAxForm
colorbar

%% compare f
exROI = 1;
beforeFrames = [20,45,110];
afterFrames = [10,120,250];
baselineFrames = 1:(10*frameRateHz);
f_tc_before = stackGetTimeCourses(img_before_reg_crop,after_mask);
f_tc_after = stackGetTimeCourses(data,after_mask);
F = mean(f_tc_after(1:20,:),1);
dff_tc_after = (f_tc_after-F)./F;
F = mean(f_tc_before(1:20,:),1);
dff_tc_before = (f_tc_before-F)./F;
befores = f_tc_before(beforeFrames,exROI);
afters = f_tc_after(afterFrames,exROI);

figure(img_fig)
subplot 334
hold on
plot(ones(1,3),befores,'k.','MarkerSize',20)
plot(ones(1,3).*2,afters,'k.','MarkerSize',20)
figXAxis([],'infustion time',[0 3],1:2,{'before','after'})
figYAxis([],'F (a.u.)',[])
figAxForm

%% plot movie tc 
figure(img_fig)

tt_movie = (1:size(dff_tc_before,1))./frameRateHz;
subplot 335
plot(tt_movie,f_tc_before,'-','LineWidth',2)
figXAxis([],'time since infusion (s)',[-5 tt_movie(end)])
figYAxis([],'F (a.u.)',[])
figAxForm
vline(beforeFrames./frameRateHz,'k:')
title('Before Infusion')
subplot 338
plot(tt_movie,dff_tc_before,'-','LineWidth',2)
figXAxis([],'time since infusion (s)',[-5 tt_movie(end)])
figYAxis([],'dF/F',[-0.2 0.5])
figAxForm
vline(beforeFrames./frameRateHz,'k:')
title('Before Infusion')

tt_movie = (1:size(dff_tc_after,1))./frameRateHz;
subplot 336
plot(tt_movie,f_tc_after,'-','LineWidth',2)
figXAxis([],'time since infusion (s)',[-5 tt_movie(end)])
figYAxis([],'F (a.u.)',[])
figAxForm
vline(afterFrames./frameRateHz,'k:')
title('During Infusion')
subplot 339
plot(tt_movie,dff_tc_after,'-','LineWidth',2)
figXAxis([],'time since infusion (s)',[-5 tt_movie(end)])
figYAxis([],'dF/F',[-0.2 7])
figAxForm
vline(afterFrames./frameRateHz,'k:')
vline(30,'r-')
title('During Infusion')

%% get later movie tc
dataPath = 'Z:\home\ashley\data\1225\widefield imaging\1225_200702\cannulatest_infusion_texasred_movie_20minafter_1';
img_after_15min = readtiff(dataPath);
dff_15min = quickDFF_wf(img_after_15min);
exptName = 'texasred_movie_15min';
writetiff(dff_15min,fullfile(fnout,exptName));

%% plot later movie
F = mean(f_tc_after(1:20,:),1);
f_tc_15min = stackGetTimeCourses(img_after_15min,after_mask);
dff_tc_15min = (f_tc_15min-F)./F;
tt_movie = ((1:size(dff_tc_15min,1))./frameRateHz)+(20*60);
figure(img_fig)
subplot 337
plot(tt_movie,dff_tc_15min,'-','LineWidth',2)
figXAxis([],'time since infusion (s)',[tt_movie(1)-5 tt_movie(end)])
figYAxis([],'dF/F',[-0.5 7])
figAxForm
title('15 Min After Infusion')

setFigParams4Print('landscape')
print(fullfile(fnout,'infusion_expt_results'),'-dpdf','-fillpage')