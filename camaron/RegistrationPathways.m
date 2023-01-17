% 2P Registation Pathways
% Look at different registration routes, choose best, then continue with TC extraction.

%% Load 920 Data 

CD = [data_base '\Data\2P_images\' expt(iexp).mouse '\' expt(iexp).date '\' expt(iexp).redImg{1}];
cd(CD);
imgMatFile = [expt(iexp).redImg{1} expt(iexp).redImg_suffix{1} '.mat'];
load(imgMatFile);
nframes_red = info.config.frames;
fprintf(['Reading run ' expt(iexp).redImg{1} '- ' num2str(min(nframes_red)) ' frames \r\n'])
data = sbxread(imgMatFile(1,1:11),0,nframes_red);

%% Get green and red image averages for 920
% size(data,1) == 2

red_920 = squeeze(data(2,:,:,:));
green_920_behav = squeeze(data(1,:,:,:));
[green_920_shifts, green_920_reg_to_data_avg] = stackRegister(green_920_behav,data_avg); %green_data = Behavior green @ 920 (1000 frames); data_avg = mean image (1 avg frames) fom middle of stack (500Frames)
green_920_avg = mean(green_920_reg_to_data_avg,3);
%figure; imagesc(green_data_avg_920)
%title('Green at 920: Route 1')

[red_920_shifts, red_920_apply_green_920_shifts] = stackRegister_MA(red_920,[],[],green_920_shifts); % red_data = red @ 920 (1000 frames); out = green @ 920 shifts_xy
red_920_avg = mean(red_920_apply_green_920_shifts,3); % Average of registered red image @ 920
%figure; imagesc(red_data_avg_920) % Show image
%title('Red at 920: Route 2&3')

%% Load 1040 Data
% size(expt(iexp).redImg,2) == 2 % 1040 run exists

CD = [data_base '\Data\2P_images\' expt(iexp).mouse '\' expt(iexp).date '\' expt(iexp).redImg{2}];
cd(CD);
imgMatFile = [expt(iexp).redImg{2} expt(iexp).redImg_suffix{2} '.mat'];
load(imgMatFile);
nframes_red = info.config.frames;
fprintf(['Reading run ' expt(iexp).redImg{2} '- ' num2str(min(nframes_red)) ' frames \r\n'])
data = sbxread(imgMatFile(1,1:11),0,nframes_red);

green_1040 = squeeze(data(1,:,:,:)); % green @ 1040
red_1040 = squeeze(data(2,:,:,:)); % red data @ 1040 (1000 Frames)

%% Run each below and view images

%% Route 1 (apply green 1040 shifts)

[green_1040_shifts, green_1040_reg_to_green_920_avg] = stackRegister(green_1040, green_920_avg); % Green 1040 to green_data_avg 920
green_1040_avg = mean(green_1040_reg_to_green_920_avg,3);
%figure; imagesc(green_data_avg_1040)
%title('Green at 1040: Route 1')

[red_1040_shifts, red_1040_apply_green_1040_shifts] = stackRegister_MA(red_1040,[],[],green_1040_shifts); % red_data = red @ 1040 (1000 frames); out = green @ 1040 shifts_xy
red_1040_avg_apply_green_1040_shifts = mean(red_1040_apply_green_1040_shifts,3); % Average of registered red image @ 1040
figure; imagesc(red_1040_avg_apply_green_1040_shifts)
title('1: Red at 1040 applied green 1040')

%% Route 2 (original route in single channel TC; apply green 920 shifts)
[red_1040_shifts_from_red_920, red_1040_reg_to_red_920_avg] = stackRegister(red_1040,red_920_avg); % red_data = red @ 1040 (1000 frames); red_data_avg = red @ 920 (1 avg frame)
red_data_avg_1040_red_920 = mean(red_1040_reg_to_red_920_avg,3); % Average of registered red image @ 1040
figure; imagesc(red_data_avg_1040_red_920)
title('2: Red at 1040 w Red 920: Original')

%% Route 3 (Register avg of red 1040)

red_1040_mean = mean(red_1040, 3); % mean red 1040
[red_1040_to_red_920_shifts, red_1040_mean_reg_to_red_920_avg] = stackRegister(red_1040_mean,red_920_avg); % red_1040_mean = red @ 1040 (1 avg frames); red_data_avg = red @ 920 (1 avg frame)
red_1040_mean_mean = mean(red_1040_mean_reg_to_red_920_avg,3); % Average of registered red image @ 1040
figure; imagesc(red_1040_mean_mean)
title('3: Red at 1040 w Red 920Avg')
%% Select best image to keep (Index for var_cell{}; 1, 2, or 3) and save data
            
var_cell = {red_1040_avg_apply_green_1040_shifts, red_data_avg_1040_red_920, red_1040_mean_mean};
red_data_avg = var_cell{}; % Select best red_data avg (var_cell{1} ,var_cell{2} or var_cell{3})!!!
green_data_avg = green_920_avg; %

save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_redData.mat']), 'green_data_avg', 'red_data_avg')

%% Save pdf of best image
data_avg = mean(data_reg(:,:,size(data_reg,3)-10000:end),3);
figure; 
subplot(2,2,1)
sz = size(data_avg);
rgb = zeros(sz(1),sz(2),3);
rgb(:,:,1) = red_data_avg./max(red_data_avg(:));
imagesc(red_data_avg);
colormap gray
subplot(2,2,2)
rgb(:,:,2) = data_avg./max(data_avg(:));
imagesc(rgb);
if size(expt(iexp).redImg,2) == 2
    title('Red at 1040; Green at 920')
else
    title('Red at 920; Green at 920')
end

print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

clear data      

