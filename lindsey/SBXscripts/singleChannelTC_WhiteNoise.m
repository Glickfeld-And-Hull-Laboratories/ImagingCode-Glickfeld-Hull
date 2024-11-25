%% get path names
close all 
clear all global
clc
date = '241104';
ImgFolder = {'002'};
time = strvcat('1509');
mouse = 'i1396';
doFromRef = 0;
ref = strvcat('002');
nrun = size(ImgFolder,2);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);


LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

%% load
tic
data = [];
CD = [LG_base '\Data\2P_images\' mouse '\' date '\' ImgFolder{1}];
%CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\AW68\two-photon imaging\' date '\' ImgFolder(irun,:)];
%CD = [LG_base '\Data\2P_images\' mouse '-KM' '\' date '_' mouse '\' ImgFolder(irun,:)];
%CD = ['\\CRASH.dhe.duke.edu\data\home\kevin\Data\2P\' date '_' mouse '\' ImgFolder(irun,:)];
cd(CD);
imgMatFile = [ImgFolder{1} '_000_000.mat'];
load(imgMatFile);
fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time '.mat'];
load(fName);
nframes = [input.counterValues(end) info.config.frames];
fprintf(['Reading run ' ImgFolder{1} '- ' num2str(min(nframes)) ' frames \r\n'])
data = sbxread(imgMatFile(1,1:11),1,min(nframes));
if size(data,1)== 2
    data = data(1,:,:,:);
end
data = squeeze(data);
toc

%% Choose register interval
nep = floor(size(data,3)./5000);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*5000):500+((i-1)*5000)),3)); title([num2str(1+((i-1)*5000)) '-' num2str(500+((i-1)*5000))]); end

%% Register data
data_avg = mean(data(:,:,10001:10500),3);
if exist(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(data,[],[],out);
    clear out outs
elseif doFromRef
    ref_str = ['runs-' ref];
    if size(ref,1)>1
        ref_str = [ref_str '-' ref(size(ref,1),:)];
    end
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_reg_shifts.mat']))
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    %load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_mask_cell.mat']))
    %load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_trialData.mat']))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
else
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
clear data out

%% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*5000):500+((i-1)*5000)),3)); title([num2str(1+((i-1)*5000)) '-' num2str(500+((i-1)*5000))]); end
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_byFrame.pdf']),'-dpdf', '-bestfit')

figure; imagesq(mean(data_reg(:,:,1:5000),3)); truesize;
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

figure; n = [1 nep]; for ii = 1:2; subplot(2,1,ii); i = n(ii); imagesc(mean(data_reg(:,:,1+((i-1)*5000):500+((i-1)*5000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*5000))]); end

%% find activated cells

data_dfof = getPixelCorrelationImage(double(data_reg));
figure; imagesc(data_dfof)
%% cell segmentation 
sz = size(data_dfof);
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
mask_data = data_dfof;

for iStim = 1:size(data_dfof,3)
    mask_data_temp = mask_data(:,:,end+1-iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0;
    bwout = imCellEditInteractiveLG(mask_data_temp);
    mask_all = mask_all+bwout;
    mask_exp = imCellBuffer(mask_all,3)+mask_all;
    close all
end
mask_cell= bwlabel(mask_all);
figure; imagesc(mask_cell)
% bwout = imCellEditInteractive(data_dfof_max);
% mask_cell = bwlabel(bwout);


%% neuropil mask and subtraction
mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof', 'mask_cell', 'mask_np')

clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_data_temp mask_exp data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 

% neuropil subtraction
down = 5;
sz = size(data_reg);

data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_reg_down  = stackGroupProject(data_reg,down);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
nCells = size(data_tc,2);
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '%s/n']) 
end
%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
[max_skew ind] =  max(x,[],1);
np_w = 0.01*ind;
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
clear data_reg data_reg_down

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell
%% white noise analysis
nFrames = input.counterValues(end);
counter2imageID = nan(1,nFrames);
diffStimFrames = diff(input.stimUpdateFrames);
nStim = length(input.stimUpdateFrames);
breaks = [0 find(diffStimFrames>10)];
for i = 1:3
    start = breaks(i);
    for ii = 1:3000
        counter2imageID(input.stimUpdateFrames(start+ii):input.stimUpdateFrames(start+ii)+diffStimFrames(start+ii)-1) = ii+(i-1)*3000;
    end
end
counter2imageID_nonan = counter2imageID;
counter2imageID_nonan(find(isnan(counter2imageID))) = [];
nStim = length(unique(counter2imageID_nonan));

stimOff = find(isnan(counter2imageID));
stimOn = find(~isnan(counter2imageID));

%load images
imagefn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Software\MWorks\Images\NoiseStimulusCapture\5min_2deg_3rep';
if exist(fullfile(imagefn, 'fullmovie.mat'))
    load(fullfile(imagefn, 'fullmovie.mat'))
else
    files = dir(imagefn);
    nfiles = size(files,1)-2;
    imagestack = [];
    for itrial = 1:nfiles
        trialfn = [imagefn '\trial_' num2str(itrial)];
        images = dir(trialfn);
        nimage = size(images,1)-2;
        for ii = 1:nimage
            pngfn = [trialfn '\frame_' num2str(ii) '.png'];
            temp = imread(pngfn);
            imagestack = cat(3,imagestack, temp(:,:,1));
            if rem(ii,100) == 0
                fprintf([num2str(ii) '\n'])
            end
        end
    end
    save(fullfile(imagefn, 'fullmovie.mat'),'imagestack', '-v7.3');
end

% for i = 1:3
%     figure;
%     imagesc(imagestack(:,:,i))
%     colormap gray
%     axis off
%     print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_whitenoise' num2str(i) '.pdf']), '-dpdf')
% end

imagestack_down = imagestack(1:10:end,1:10:end,:);
imagesz = size(imagestack_down);

meanOn = mean(npSub_tc(stimOn,:),1);
meanOff = mean(npSub_tc(stimOff,:),1);
respDiff = (meanOn-meanOff)./meanOff;

tt = 0.033:0.033:0.033*size(np_tc,1);
figure; plot(tt,mean(np_tc,2)); vline([find(counter2imageID == 3000)./30 find(counter2imageID == 3001)./30 find(counter2imageID == 6000)./30 find(counter2imageID == 6001)./30])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_npTCAvg.pdf']), '-dpdf')
figure; histogram(respDiff); xlabel('(On-Off)/Off'); ylabel('Number of cells');
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respDiff.pdf']), '-dpdf')


nCells = size(npSub_tc,2);
meanF = mean(npSub_tc,1);
stdF = std(npSub_tc,[],1);
cellRespInd = cell(1,nCells);
eventMovieAvg = nan(imagesz(1),imagesz(2),10,nCells);
for iCell = 1:nCells
    fprintf(['Cell # ' num2str(iCell) '\n'])
    thresh_temp = find(npSub_tc(:,iCell)>meanF(iCell)+stdF(iCell)*2);
    cellRespInd{iCell} = thresh_temp(find(diff(thresh_temp)>10)+1);

    nEvent = length(cellRespInd{iCell});
    eventMovie = nan(imagesz(1),imagesz(2),10,nEvent);
    for i = 1:nEvent
        if cellRespInd{iCell}(i)>60
            frameList = counter2imageID(cellRespInd{iCell}(i)-9:cellRespInd{iCell}(i));
            for ii = 1:10
                if ~isnan(frameList(ii)) & frameList(ii)<imagesz(3)
                    eventMovie(:,:,ii,i) = imagestack_down(:,:,frameList(ii));
                end
            end
        end
        if rem(i,100) == 0
            fprintf([num2str(i) '\n'])
        end
    end
    eventMovieAvg(:,:,:,iCell) = mean(eventMovie,4,'omitnan');
end
tt = [-297:33:0];
respInd = find(respDiff>0.1);
cells = respInd([8 10 12 15 17 19]);

figure;
for iCell = 73:nCells

    %iCell = cells(i);
    subplot(6,6,iCell-72)
    imagesc(eventMovieAvg(:,:,9,iCell),[145 190])
    colormap gray
    axis off
    %sgtitle(['Cell #' num2str(iCell) '; deltaResp: ' num2str(respDiff(iCell))])
    %print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respFrameImg_Cell' num2str(iCell) '.pdf']), '-dpdf')
end
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respFrameImg_Cell73-84.pdf']), '-dpdf')

tt = 0.033:0.033:0.033*size(npSub_tc,1);
i =1;
figure;
iCell = cells(i);
plot(tt,npSub_tc(:,iCell))
hline(meanF(iCell))
hline(meanF(iCell)+stdF(iCell)*2)
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCthresh_Cell' num2str(iCell) '.pdf']), '-dpdf')


%full movie
sz = size(imagestack);
eventFrameAvg = zeros(sz(1),sz(2),6);
for i = 1:6
    iCell = cells(i);
    nEvent = length(cellRespInd{iCell});
    eventFrame = nan(sz(1),sz(2),nEvent);
    for i = 1:nEvent
        if cellRespInd{iCell}(i)>60
            frame = counter2imageID(cellRespInd{iCell}(i)-1);
            if ~isnan(frame) & frame<sz(3)
                eventFrame(:,:,i) = imagestack(:,:,frame);
            end
        end
    end
    eventFrameAvg(:,:,iCell) = mean(eventFrame,3,'omitnan');
end

for i = 1
    figure;
    iCell = cells(i);
    imagesc(eventFrameAvg(:,:,iCell),[145 190])
    truesize;
    colormap gray
    axis off
    sgtitle(['Cell #' num2str(iCell) '; deltaResp: ' num2str(respDiff(iCell))])
end