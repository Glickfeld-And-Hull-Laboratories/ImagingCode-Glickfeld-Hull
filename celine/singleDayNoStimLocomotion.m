%% Load, register, segment and neuropil correct 2P data
clc
close all
clear all global

expName = 'GRAB';

%Specific experiment information
date = '250512'
ImgFolder = '007';
time = '1506';
mouse = 'i2188';
frame_rate = 15;

fn_base = findIsilon;
data_fn = fullfile(fn_base, 'home', 'ACh');
data_fn = fullfile(data_fn, 'Data', '2p_data');
mworks_fn = fullfile(fn_base, 'Behavior', 'Data');

isilonName = findIsilon;
%Path names
if computer == 'GLNXA64'
    base = fullfile(isilonName, '/home/ACh\Analysis\2p_analysis');
    
else
    base = fullfile(isilonName, '\home\ACh\Analysis\2p_analysis');
      
end

d=string(datetime('today'));
fnout= fullfile(base,expName,mouse,d);

mkdir(fnout);
cd(fnout)
clear d sess_title



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
save(fullfile(fnout,'reg_shifts.mat'), 'data_reg_avg', 'out', 'data_avg')
save(fullfile(fnout,'input.mat'), 'input')
%Test registration
%    a. Make sure first and last images are sharp and similar. Focus on vasculature and nuclei- not cells since F changes
ind = [1 nep];
for i = 1:length(ind) 
    subplot(2,1,i); 
    ix = ind(i);
    imagesc(mean(data_reg(:,:,1+((ix-1)*nskip):nframes+((ix-1)*nskip)),3)); 
    title([num2str(1+((ix-1)*nskip)) '-' num2str(nframes+((ix-1)*nskip))]); 
end
print(fullfile(fnout,'FOV_first&last.pdf'), '-dpdf')
%    b. Average of all frames should be sharp
figure;
imagesc(data_reg_avg); 
print(fullfile(fnout,'FOV_avg.pdf'), '-dpdf')

clear data
%% Align to locomotion
% This section aligns imaging data with locomotion events to analyze neural activity changes
% during movement initiation
wheel_speed = wheelSpeedCalc(input,32,'orange');  % Calculate wheel speed from behavioral data

% Plot wheel speed and mean fluorescence to visualize relationship between movement and activity
figure; 
plot(wheel_speed)
hold on;
plot(squeeze(mean(mean(data_reg,1),2)))  % Plot average fluorescence over time
title('Wheel Speed and Mean Fluorescence Signal')
xlabel('Frame Number')
ylabel('Speed / Fluorescence')
legend('Wheel Speed', 'Mean Fluorescence')

% Find locomotion onsets with minimum delay of 3 seconds from previous locomotion
delay = frame_rate*3;  % 3 seconds delay in frames
ind = find(wheel_speed>10);  % Find frames where mouse is running (speed > 10)
ind_use = [];  % Initialize array for locomotion onset frames
for i = 1:length(ind)
    if ~isempty(ind_use)
        if ind_use(end)<ind(i)-delay  % If it's been more than 'delay' frames since last onset
            if max(wheel_speed(ind(i)-delay:ind(i)-1))<10  % And mouse was not running during the delay
                ind_use = [ind_use ind(i)];  % Add this frame as a locomotion onset
            end
        end
    elseif ind(i)>delay  % For the first onset, check if it's far enough from start
        if max(wheel_speed(ind(i)-delay:ind(i)-1))<10  % And mouse was not running
            ind_use = [ind_use ind(i)];  % Add as onset
        end
    end
    % No 'else' clause - we skip onsets that are too close to the start
end

vline(ind_use)  % Mark locomotion onsets on the plot
title('Wheel Speed with Identified Locomotion Onsets')

% Extract data segments around locomotion onsets with proper handling of boundaries
sz = size(data_reg);
data_align = nan(sz(1),sz(2),delay*2,length(ind_use));  % Initialize 4D array
for i = 1:length(ind_use)
    start_idx = ind_use(i)-delay;
    end_idx = ind_use(i)+delay-1;
    
    % Check if start index is valid (not too close to beginning)
    if start_idx < 1
        valid_start = 1;
        offset_start = 1 - start_idx;  % How many frames we're missing at the start
    else
        valid_start = start_idx;
        offset_start = 0;
    end
    
    % Check if end index is valid (not too close to end)
    if end_idx > sz(3)
        valid_end = sz(3);
    else
        valid_end = end_idx;
    end
    
    % Calculate how many frames we can actually get
    n_frames = valid_end - valid_start + 1;
    
    % Place the valid frames in the correct position in the output array
    data_align(:,:,offset_start+1:offset_start+n_frames,i) = data_reg(:,:,valid_start:valid_end);
end

% Process extracted data
data_align_p = permute(data_align,[1 2 4 3]);  % Rearrange dimensions
data_align_f = mean(data_align_p(:,:,:,1:delay),4);  % Calculate baseline (pre-locomotion)
data_align_dfof = (data_align_p-data_align_f)./data_align_f;  % Calculate dF/F
data_align_dfof_avg = squeeze(nanmean(data_align_dfof,3));  % Average across trials
data_align_avg = squeeze(nanmean(data_align_p,3));  % Average raw data across trials

% Downsample and visualize results
data_align_dfof_down = stackGroupProject(data_align_dfof_avg,frame_rate);  % Downsample
n = (delay*2)./frame_rate;
figure;
for i = 1:n  % Display downsampled data
    subplot(3,2,i)
    imagesc(data_align_dfof_down(:,:,i))
    clim([0 .3])  % Set color limits
    title(['Time bin ' num2str(i) ': ' num2str((i-1-n/2)) 's to ' num2str(i-n/2) 's'])
end
sgtitle('Neural Activity (dF/F) Around Locomotion Onset')  % Add super title

% Plot and visualize average response
tt = -delay:delay-1;  % Time axis
figure; 
plot(tt,squeeze(mean(mean(data_align_dfof_avg,1),2)))  % Plot average dF/F over time
title('Average Neural Response Around Locomotion Onset')
xlabel('Time from Running Onset (frames)')
ylabel('Mean dF/F')
grid on

figure; 
subplot(2,1,1)
imagesc(mean(data_align_dfof_avg(:,:,delay/2:delay),3))  % Show pre-running average
title('Pre-Running Neural Activity (averaged)')
clim([0 .3])
colorbar
xlabel('X Position (pixels)')
ylabel('Y Position (pixels)')

subplot(2,1,2)
imagesc(mean(data_align_dfof_avg(:,:,delay+1:end),3))  % Show during-running average
title('During-Running Neural Activity (averaged)')
clim([0 .3])
colorbar
xlabel('X Position (pixels)')
ylabel('Y Position (pixels)')
sgtitle(['Mouse ' mouse ': Activity Before vs. During Locomotion'])

data_dfof_avg = mean(data_align_dfof_avg(:,:,delay+1:end),3);  % Average of post-locomotion frames

% Filter data to make cells more visible for segmentation
myfilter = fspecial('gaussian',[20 20], 0.5);  % Create Gaussian filter
data_dfof_avg_all = imfilter(data_dfof_avg,myfilter);  % Apply filter
figure; 
movegui('center'); 
imagesc(data_dfof_avg_all);  % Display filtered image
title(['Filtered Running-Related Activity Map - Mouse ' mouse ' ' date])
colorbar
xlabel('X Position (pixels)')
ylabel('Y Position (pixels)')
print(fullfile(fnout,'runAlignFOV_filtered.pdf'), '-dpdf')  % Save as PDF
save(fullfile(fnout,'stimActFOV.mat'),'data_dfof_avg_all','ind_use')  % Save data

% Create cell masks based on activity threshold
thresh = 0.05;  % Threshold for activity
mask_data = data_dfof_avg_all>thresh;  % Create binary mask where activity > threshold
% Remove edges to avoid edge artifacts
mask_data(1:10,:) = 0;
mask_data(:,1:10) = 0;
mask_data(sz(1)-9:sz(1),:) = 0;
mask_data(:,sz(2)-9:sz(2)) = 0;
figure; 
imagesc(mask_data)  % Display mask
title(['Activity Threshold Mask (>' num2str(thresh) ' dF/F)'])
xlabel('X Position (pixels)')
ylabel('Y Position (pixels)')

mask_cell = bwlabel(mask_data);  % Label connected components (cells)
figure; 
imagesc(mask_cell)  % Display labeled cells
title(['Identified Cell Regions - Mouse ' mouse ' ' date])
xlabel('X Position (pixels)')
ylabel('Y Position (pixels)')
colorbar
colormap(parula)
print(fullfile(fnout,'masks.pdf'), '-dpdf')  % Save as PDF

save(fullfile(fnout,'mask_cell.mat'), 'mask_data', 'mask_cell', 'thresh')  % Save masks
%clear data_align data_align_p data_align_f data_align_dfof data_align_dfof_avg data_align_dfof_down  % Free memory

%% Extract cell timecourses
data_tc = stackGetTimeCourses(data_reg, mask_cell);  % Extract average timecourse for each cell
            % Timecourses are nFrames x nCells
fprintf(['data_tc is ' num2str(size(data_tc))])  % Report dimensions
[nFrames, nCells] = size(data_tc);  % Get dimensions

save(fullfile(fnout, 'TCs.mat'), 'data_tc')  % Save timecourses

%% Analyze cell activity around locomotion onset
% This section analyzes how individual cells respond to locomotion onset
data_align = nan(delay*2,nCells,length(ind_use));  % Initialize array for aligned timecourses

% Extract cell activity around locomotion onsets with proper boundary handling
for i = 1:length(ind_use)
    start_idx = ind_use(i)-delay;
    end_idx = ind_use(i)+delay-1;
    
    % Check if start index is valid (not too close to beginning)
    if start_idx < 1
        valid_start = 1;
        offset_start = 1 - start_idx;  % How many frames we're missing at the start
    else
        valid_start = start_idx;
        offset_start = 0;
    end
    
    % Check if end index is valid (not too close to end)
    if end_idx > nFrames
        valid_end = nFrames;
    else
        valid_end = end_idx;
    end
    
    % Calculate how many frames we can actually get
    n_frames = valid_end - valid_start + 1;
    
    % Place the valid frames in the correct position in the output array
    data_align(offset_start+1:offset_start+n_frames,:,i) = data_tc(valid_start:valid_end,:);
end

% Define analysis windows
base_win = delay-frame_rate+1:delay;  % 1-second window before locomotion onset
resp_win = delay+1:delay+frame_rate;  % 1-second window after locomotion onset
data_f = nanmean(data_align(base_win,:,:),1);  % Calculate baseline (using nanmean to handle NaNs)
data_dfof = bsxfun(@rdivide,bsxfun(@minus,data_align,data_f),data_f);  % Calculate dF/F

% Statistical test to find responsive cells
[h p] = ttest(nanmean(data_dfof(resp_win,:,:),1),nanmean(data_dfof(base_win,:,:),1),'Dim',3,'tail','right');  
resp_ind = find(h);  % Find cells with significant increase in activity

% Plot average response of responsive cells
tt_ms = tt.*(1000./frame_rate);  % Convert time to milliseconds
figure; 
plot(tt_ms,nanmean(nanmean(data_dfof(:,resp_ind,:),2),3))  % Plot average
xlabel('Time from Running Onset (ms)')
ylabel('dF/F')
title(['Average Response of Locomotion-Responsive Masks (n=' num2str(length(resp_ind)) ') - Mouse ' mouse ' ' date])
grid on
print(fullfile(fnout,'avgTC_respCells.pdf'), '-dpdf')  % Save as PDF

% Save analysis results
save(fullfile(fnout,'respData.mat'), 'data_align', 'h','resp_ind','resp_win','base_win','tt')

%% Create video of average activity around locomotion onset
% This section creates a video showing dF/F changes from 3 seconds before 
% to 3 seconds after running onset

% First, ensure we have the data we need
if ~exist('data_align_dfof_avg', 'var')
    disp('Loading required data...')
    load(fullfile(fnout,'stimActFOV.mat'))
    % If data isn't available, would need to recompute it
end

% Create a directory for the video
vid_dir = fullfile(fnout, 'videos');
if ~exist(vid_dir, 'dir')
    mkdir(vid_dir)
end

% Set up video parameters
vid_filename = fullfile(vid_dir, [mouse '_' date '_RunningOnset.mp4']);
frame_rate_output = 10; % Output video frame rate (fps)
v = VideoWriter(vid_filename, 'MPEG-4');
v.FrameRate = frame_rate_output;
v.Quality = 95;
open(v);

% Set up figure for video frames
fig = figure('Position', [100, 100, 800, 600], 'Color', 'w');

% Create colormap for better visualization
cmap_hot = hot(256);
cmap_hot = cmap_hot(end:-1:1,:); % Flip hot colormap

% Set color limits for consistency
clim_min = 0;
clim_max = 0.3; % Adjust based on your dF/F range

% Create time labels
time_points = -delay:1:delay-1;
time_seconds = time_points / frame_rate;

% Prepare average dF/F data
data_to_plot = data_align_dfof_avg;

% Create the video frames
disp('Creating video frames...')
for frame_idx = 1:size(data_to_plot, 3)
    % Get current time in seconds
    curr_time = time_seconds(frame_idx);
    
    % Clear the figure for this frame
    clf;
    
    % Plot the current frame
    imagesc(data_to_plot(:,:,frame_idx));
    colormap(cmap_hot);
    clim([clim_min clim_max]);
    
    % Add colorbar and labels
    cb = colorbar;
    ylabel(cb, 'dF/F');
    
    % Add title showing time relative to locomotion onset
    if curr_time < 0
        title(sprintf('%.1f seconds before locomotion onset', abs(curr_time)), 'FontSize', 14);
    elseif curr_time == 0
        title('Locomotion onset (t = 0)', 'FontSize', 14, 'FontWeight', 'bold');
    else
        title(sprintf('%.1f seconds after locomotion onset', curr_time), 'FontSize', 14);
    end
    
    % Add a timeline bar to visualize progress
    axes('Position', [0.2, 0.05, 0.6, 0.03], 'Color', 'none');
    xlim([-delay/frame_rate, delay/frame_rate]);
    ylim([0, 1]);
    hold on;
    % Draw timeline
    plot(time_seconds, zeros(size(time_seconds)), 'k-', 'LineWidth', 2);
    % Mark current position
    plot(curr_time, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    % Mark zero point
    plot(0, 0, 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
    
    % Add text labels
    text(-delay/frame_rate, 0.5, sprintf('%.1fs', -delay/frame_rate), 'HorizontalAlignment', 'left');
    text(0, 0.5, '0s', 'HorizontalAlignment', 'center');
    text(delay/frame_rate, 0.5, sprintf('%.1fs', delay/frame_rate), 'HorizontalAlignment', 'right');
    
    % No axis ticks needed for timeline
    set(gca, 'XTick', [], 'YTick', [], 'Box', 'off');
    
    % Add text in the corner showing mouse ID and date
    annotation('textbox', [0.02, 0.02, 0.2, 0.05], 'String', ...
        ['Mouse: ' mouse ', Date: ' date], 'EdgeColor', 'none', 'FontSize', 10);
    
    % Capture and write the frame
    drawnow;
    frame = getframe(fig);
    writeVideo(v, frame);
end

% Close the video writer
close(v);
disp(['Video saved to: ' vid_filename]);

% Create additional visualization showing key time points for quick reference
figure('Position', [100, 100, 1200, 400]);

% Time points to show (before, onset, after)
key_times = [-3, -2, -1, 0, 1, 2, 3]; % in seconds
key_frames = round((key_times * frame_rate) + delay + 1);
key_frames = max(1, min(size(data_to_plot, 3), key_frames)); % Ensure within bounds

for i = 1:length(key_times)
    subplot(1, length(key_times), i);
    frame_idx = key_frames(i);
    imagesc(data_to_plot(:,:,frame_idx));
    colormap(cmap_hot);
    clim([clim_min clim_max]);
    
    if key_times(i) < 0
        title(sprintf('%.1fs before', abs(key_times(i))));
    elseif key_times(i) == 0
        title('Onset (t=0)');
    else
        title(sprintf('%.1fs after', key_times(i)));
    end
    
    % Only add colorbar to the last subplot
    if i == length(key_times)
        colorbar;
    end
    
    % Remove tick labels for cleaner appearance
    set(gca, 'XTick', [], 'YTick', []);
end

sgtitle(['Neural Activity Around Locomotion Onset - Mouse ' mouse], 'FontSize', 14);
print(fullfile(fnout, 'RunningOnset_KeyFrames.pdf'), '-dpdf', '-bestfit');

disp('Analysis complete!');