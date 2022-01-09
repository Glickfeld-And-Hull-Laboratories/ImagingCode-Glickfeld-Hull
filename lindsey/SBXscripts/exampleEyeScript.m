%paths
mouse = 'i484';
date = '220107';
run = '001';
time = '1029';
data_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\camaron';
CD = [data_pn '\Data\2P_images\' mouse '\' date '\' run];
cd(CD);
fn = [run '_000_000_eye.mat'];

%load data
data_temp = load(fn);
data_temp = squeeze(data_temp.data);

%crop frames to match mworks data
fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time '.mat'];
load(fName);
nFrames = input.counterValues{end}(end);
data = data_temp(:,:,1:nFrames);      % the raw images...

%% Crop image to isolate pupil 
%(bright spots can be mistaken for pupil)
[data_crop rect] = cropEyeData(data);

%% measure pupil position/diameter
rad_range = [5 20]; %adjust to expected range of pupil size (if low end is too small then may find noisy bright stuff)
Eye_data = extractEyeData(data_crop,rad_range);
%if pupil not found reliably, adjust the image cropping or the rad_range

%% align to stimulus presentation
[rad centroid] = alignEyeData(Eye_data,input);



    
            
