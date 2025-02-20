%% Load, register, segment and neuropil correct 2P data
clc
close all
clear all global

%Path names

fn_base = findIsilon;
base_fn = fullfile(fn_base, 'home', 'LindseyW'); %edit
data_fn = fullfile(base_fn, 'Data', 'Widefield'); %edit
mworks_fn = fullfile(fn_base, 'Behavior', 'Data'); 
fnout = fullfile(base_fn, 'Analysis', 'Widefield'); %edit

%Specific experiment information
date = '250205'; %edit
ImgFolder = '001'; %edit
time = '1151'; %edit
mouse = 'i2194'; %edit
frame_rate = 10; %15; %edit
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' ImgFolder];

%Load data
%Load mworks data- this has information about experiment (e.g. the visual stimuli presented and synchronization with the microscope)
%fName = fullfile(mworks_fn, ['data-' mouse '-' date '-' time '.mat']);
fName = fullfile("Z:\All_Staff\Behavior\Data\data-i2194-250205-1511.mat") %temp
load(fName);

%Load tiff stack
%expt_fn = fullfile(data_fn,datemouse,datemouserun); %edit
expt_fn = fullfile('Z:\All_Staff\home\ACh\Data\WF_data\i2194\Pre_injection\i2194_250205_NoStim_1\i2194_250205_NoStim_1_MMStack_Pos0.ome.tif') %temp
data = readtiff(expt_fn); %loads the .sbx files with imaging data (path, nframes to skip, nframes to load)
%Data is nYpix x nXpix x nframes. 
fprintf(['Data is ' num2str(size(data)) '\n'])


%% Align to locomotion
wheel_speed = wheelSpeedCalc(input,32,'purple');

figure; plot(wheel_speed)
hold on;
plot(squeeze(mean(mean(data,1),2)))


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

sz = size(data);
data_align = nan(sz(1),sz(2),delay*2,length(ind_use));
for i = 1:length(ind_use)
    if ind_use(i)-delay>0
    if ind_use(i)+delay-1<sz(3)
        data_align(:,:,:,i) = data(:,:,ind_use(i)-delay:ind_use(i)+delay-1);
    else
        n = sz(3)-(ind_use(i)-delay-1);
        data_align(:,:,1:n,i) = data(:,:,ind_use(i)-delay:end);
    end
    end
end
data_align_p = permute(data_align,[1 2 4 3]);
data_align_f = mean(data_align_p(:,:,:,1:delay),4);
data_align_dfof = (data_align_p-data_align_f)./data_align_f;
data_align_dfof_avg = squeeze(nanmean(data_align_dfof,3));
data_align_avg = squeeze(nanmean(data_align_p,3));
writetiff(data_align_avg,fullfile(fnout, datemouse, datemouserun, [datemouserun '_avgMovie.tif'])) %edit
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
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_runAlignFOV_filtered.pdf']), '-dpdf') %edit
%save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_stimActFOV.mat']),'data_dfof_avg_all','ind_use') %edit

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

%save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_mask_cell.mat']), 'mask_data', 'mask_cell', 'thresh')
clear data_align data_align_p data_align_f data_align_dfof data_align_dfof_avg data_align_dfof_down

%% Extract cell timecourses
data_tc = stackGetTimeCourses(data, mask_cell); %applies mask to stack (averages all pixels in each frame for each cell) to get timecourses
            %Timecourses are nFrames x nCells
fprintf(['data_tc is ' num2str(size(data_tc))]) 
[nFrames, nCells] = size(data_tc);

%save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_TCs.mat']), 'data_tc')

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

%save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_respData.mat']), 'data_align', 'h','resp_ind','resp_win','base_win','tt')

