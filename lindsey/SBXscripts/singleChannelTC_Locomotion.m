%% Load, register, segment and neuropil correct 2P data
clc
close all
clear all global

%Path names

fn_base = findIsilon;
lg_fn = fullfile(fn_base, 'home', 'lindsey');
data_fn = fullfile(lg_fn, 'Data', '2P_images');
mworks_fn = fullfile(fn_base, 'Behavior', 'Data');
fnout = fullfile(lg_fn, 'Analysis', '2P');

%Specific experiment information
date = '250207';
ImgFolder = '002';
time = '1659';
mouse = 'i2188';
frame_rate = 15;
run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];

%Load 2P data
%Load mworks data- this has information about experiment (e.g. the visual stimuli presented and synchronization with the microscope)
fName = fullfile(mworks_fn, ['data-' mouse '-' date '-' time '.mat']);
load(fName);
%Load 2P metadata- this has information about the information about the imaging session (e.g. frame count, zoom)
CD = fullfile(data_fn, mouse, date, ImgFolder);
cd(CD);
imgMatFile = [ImgFolder '_000_000.mat'];
load(imgMatFile);
%Load 2P images
totframes = input.counterValues{end}(end); %this is from the mworks structure- finds the last value clocked for frame count
fprintf(['Reading ' num2str(totframes) ' frames \r\n'])
data = sbxread([ImgFolder '_000_000'],0,totframes); %loads the .sbx files with imaging data (path, nframes to skip, nframes to load)
%Data is nPMT x nYpix x nXpix x nframes. 
fprintf(['Data is ' num2str(size(data)) '\n'])
%When imaging a single channel, nPMT = 1, so squeeze:
data = squeeze(data);

%% Register 2P data
 
%Goal here is to remove X-Y movement artifacts
%1. Find a stable target
%    a. Plot average of 500 frames throughout stack
nframes = 500; %nframes to average for target
nskip = double(ceil(totframes/5)); %nframes to skip for each average

nep = floor(size(data,3)./nskip);
[n, n2] = subplotn(nep); 
figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:nep 
    subplot(n,n2,i); 
    imagesc(mean(data(:,:,1+((i-1)*nskip):nframes+((i-1)*nskip)),3)); 
    title([num2str(1+((i-1)*nskip)) '-' num2str(nframes+((i-1)*nskip))]); 
end

%    b. GUI to select target image- choose one that is sharp and close to center of stack
f=gcf;
w = waitforbuttonpress; %click on subplot
if w == 0
    axesClicked = gca;
    allAxes = flipud(findobj(f.Children,'Type','axes'));
    numClicked = find(axesClicked==allAxes);
    close all
end
fprintf(['Selected subplot ' num2str(numClicked) '\n'])
%    c. Create target image
data_avg = mean(data(:,:,1+((numClicked-1)*nskip):nframes+((numClicked-1)*nskip)),3); %average 500 frames to make target
%2. stackRegister minimizes the difference of each frame from the target
[out, data_reg] = stackRegister(data,data_avg);
%New average image after registration
data_reg_avg = mean(data_reg,3);
%Save registration shifts and target, and mworks data
mkdir(fullfile(fnout, datemouse, datemouserun))
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_reg_shifts.mat']), 'data_reg_avg', 'out', 'data_avg')
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_input.mat']), 'input')
%Test registration
%    a. Make sure first and last images are sharp and similar. Focus on vasculature and nuclei- not cells since F changes
ind = [1 nep];
for i = 1:length(ind) 
    subplot(2,1,i); 
    ix = ind(i);
    imagesc(mean(data_reg(:,:,1+((ix-1)*nskip):nframes+((ix-1)*nskip)),3)); 
    title([num2str(1+((ix-1)*nskip)) '-' num2str(nframes+((ix-1)*nskip))]); 
end
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_FOV_first&last.pdf']), '-dpdf')
%    b. Average of all frames should be sharp
figure;
imagesq(data_reg_avg); 
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_FOV_avg.pdf']), '-dpdf')

clear data
%% Align to locomotion
wheel_speed = wheelSpeedCalc(input,32,'purple');

figure; plot(wheel_speed)
hold on;
plot(squeeze(mean(mean(data_reg,1),2)))

delay = frame_rate*3;
ind = find(wheel_speed>10);
ind_use = [];
for i = 1:length(ind)
    if ~isempty(ind_use)
        if ind_use(end)<ind(i)-delay
            if max(wheel_speed(ind(i)-delay:ind(i)-1))<10
                ind_use = [ind_use ind(i)];
            end
        end
    elseif ind(i)>delay
        if max(wheel_speed(ind(i)-delay:ind(i)-1))<10
            ind_use = [ind_use ind(i)];
        end
    else
        ind_use = [ind_use ind(i)];
    end
end

vline(ind_use)

sz = size(data_reg);
data_align = nan(sz(1),sz(2),delay*2,length(ind_use));
for i = 1:length(ind_use)
    if ind_use(i)+delay-1<sz(3)
        data_align(:,:,:,i) = data_reg(:,:,ind_use(i)-delay:ind_use(i)+delay-1);
    else
        n = sz(3)-(ind_use(i)-delay-1);
        data_align(:,:,1:n,i) = data_reg(:,:,ind_use(i)-delay:end);
    end
end
data_align_p = permute(data_align,[1 2 4 3]);
data_align_f = mean(data_align_p(:,:,:,1:delay),4);
data_align_dfof = (data_align_p-data_align_f)./data_align_f;
data_align_dfof_avg = squeeze(nanmean(data_align_dfof,3));
data_align_avg = squeeze(nanmean(data_align_p,3));
writetiff(data_align_avg,fullfile(fnout, datemouse, datemouserun, [datemouserun '_avgMovie.tif']))
data_align_dfof_down = stackGroupProject(data_align_dfof_avg,frame_rate);
n = (delay*2)./frame_rate;
figure;
for i = 1:n
    subplot(3,2,i)
    imagesc(data_align_dfof_down(:,:,i))
    clim([0 .3])
end
tt = -delay:delay-1;
figure; plot(tt,squeeze(mean(mean(data_align_dfof_avg,1),2)))
figure; 
subplot(2,1,1)
imagesc(mean(data_align_dfof_avg(:,:,delay/2:delay),3))
title('pre')
clim([0 .3])
subplot(2,1,2)
imagesc(mean(data_align_dfof_avg(:,:,delay+1:end),3))
title('post')
clim([0 .3])

data_dfof_avg = mean(data_align_dfof_avg(:,:,delay+1:end),3);

clear data_dfof
%        Filtering data helps make cells more visible for selection
myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_avg_all = imfilter(data_dfof_avg,myfilter);
figure; movegui('center'); imagesc(data_dfof_avg_all);
title([mouse ' ' date])
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_runAlignFOV_filtered.pdf']), '-dpdf')
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_stimActFOV.mat']),'data_dfof_avg_all','ind_use')

thresh = 0.08;
mask_data = data_dfof_avg_all>thresh;
mask_data(1:10,:) = 0;
mask_data(:,1:10) = 0;
mask_data(sz(1)-9:sz(1),:) = 0;
mask_data(:,sz(2)-9:sz(2)) = 0;
figure; imagesc(mask_data)
mask_cell = bwlabel(mask_data);
figure; imagesc(mask_cell)
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_masks.pdf']), '-dpdf')

save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_mask_cell.mat']), 'mask_data', 'mask_cell', 'thresh')
clear data_align data_align_p data_align_f data_align_dfof data_align_dfof_avg data_align_dfof_down
%% visual response
if input.gratingContrast 
    [cStimOn cStimOff] = photoFrameFinder_Sanworks(info.frame);
    ntrials = length(cStimOn);
    prewin = input.nScansOff./2;
    postwin = input.nScansOn;
    data_resp = nan(sz(1),sz(2),prewin+postwin,ntrials);
    wheel_resp = nan(prewin+postwin,ntrials);
    for i = 1:ntrials
        data_resp(:,:,:,i) = data_reg(:,:,cStimOn(i)-prewin:cStimOn(i)+postwin-1);
        wheel_resp(:,i) = wheel_speed(1,cStimOn(i)-prewin:cStimOn(i)+postwin-1);
    end
    data_resp_f = mean(data_resp(:,:,1:prewin,:),3);
    data_resp_dfof = (data_resp-data_resp_f)./data_resp_f;
    
    data_tc = squeeze(mean(mean(mean(data_resp_dfof,1),2),4));
    wheel_tc = squeeze(mean(wheel_resp,2));
    figure; plot(data_tc./max(data_tc(:)));
    hold on
    plot(wheel_tc./max(wheel_tc(:)));
end


%% Extract cell timecourses
data_tc = stackGetTimeCourses(data_reg, mask_cell); %applies mask to stack (averages all pixels in each frame for each cell) to get timecourses
            %Timecourses are nFrames x nCells
fprintf(['data_tc is ' num2str(size(data_tc))]) 
[nFrames, nCells] = size(data_tc);

save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_TCs.mat']), 'data_tc')

%% make tuning curves
data_align = nan(delay*2,nCells,length(ind_use));
for i = 1:length(ind_use)
    if ind_use(i)+delay-1<sz(3)
        data_align(:,:,i) = data_tc(ind_use(i)-delay:ind_use(i)+delay-1,:);
    else
        n = sz(3)-(ind_use(i)-delay-1);
        data_align(1:n,:,i) = data_tc(ind_use(i)-delay:end,:);
    end
end

base_win = delay-frame_rate+1:delay;
resp_win = delay+1:delay+frame_rate;
data_f = mean(data_align(base_win,:,:),1);
data_dfof = bsxfun(@rdivide,bsxfun(@minus,data_align,data_f),data_f);

[h p] = ttest(nanmean(data_dfof(resp_win,:,:),1),nanmean(data_dfof(base_win,:,:),1),'Dim',3,'tail','right');
resp_ind = find(h);

tt_ms = tt.*(1000./frame_rate);
figure; plot(tt_ms,nanmean(nanmean(data_dfof(:,resp_ind,:),2),3))
xlabel('Time from running onset')
ylabel('dF/F')
title([mouse ' ' date])
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_avgTC_respCells.pdf']), '-dpdf')

save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_respData.mat']), 'data_align', 'h','resp_ind','resp_win','base_win','tt')

