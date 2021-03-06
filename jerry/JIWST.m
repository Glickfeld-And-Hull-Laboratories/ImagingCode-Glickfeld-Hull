%Just Imaging With Split Transform
clear;
%% get path names
% date is day2
date = '200108'; 
ImgFolder = strvcat('003');
time = strvcat('1158');
mouse = 'i1316';
% set align to 1 if aligning any day to day 1. set to 0 if you are just
% working with day 1
alignToRef = 1;
% ref_date is day1
ref_date = '200106';
ref_run = strvcat('003');
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);
ref_str = catRunName(ref_run, size(ref_run,1));
% Type in your own stuff here:
gl_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\jerry\2P_Imaging_Grace';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\jerry\Analysis\2P';
behav_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';

%% load and register
data = [];
clear temp
trial_n = [];
offset = 0;
%     loading data_dfof_max and data_reg_avg from day 1 (_reg)
    load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimActFOV.mat']))
    load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    data_dfof_max_reg = data_dfof_max;
    data_avg_reg = data_reg_avg;
%     loading data_dfof_max and data_reg_avg from day 1 (_ref)
    load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_stimActFOV.mat']))
    load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_reg_shifts.mat']))
    data_dfof_max_ref = data_dfof_max;
    data_avg_ref = data_reg_avg;

%% Thresholding    
%     defining day 1 and day 2 data_reg_avg and taking out the brightest
%     points of the image
    reg = data_avg_reg;
    ref = data_avg_ref;
    reg(find(reg>2000)) = 0;
    reg = (reg./max(max(abs(reg))));
    ref(find(ref>2000)) = 0;
    ref = (ref./max(max(abs(ref)))); 
%     figure(1);clf;imshowpair(reg, ref, 'montage');
%     BW1 = edge(ref,'Canny', [0.01 0.15]); BW2 = edge(reg,'Canny', [0.01 0.1]);
%     figure(2);clf;imshowpair(BW1, BW2, 'montage');
%     refTL = ref(1:512/2, 1:796/2);
%     refTR = ref(1:512/2, 796/2:end);
%     refBL = ref(512/2:end, 1:796/2);
%     refBR = ref(512/2:end, 796/2:end);
%% Transforming
%     selecting the matching landmarks on day 1 and day 2 data_reg_avg and
%     shifting the position of day 2 data_reg_avg and data_dfof_max accordingly
    [reg2ref, reg2ref_dfof] = splitTransform(ref, reg, data_dfof_max_reg);

%% Creating RGB Images
%     Red/green images created here. rgb_ref2ref is using data_reg_avg, and rgb_reg2ref_dfof is using
%     data_dfof_max
    sz_target = size(reg);
    rgb_reg2ref = zeros(sz_target(1), sz_target(2), 3);
    rgb_reg2ref_dfof = zeros(sz_target(1), sz_target(2), 3);
    rgb_reg2ref(:,:,1) = ref;
    rgb_reg2ref(:,:,2) = reg2ref;
    rgb_reg2ref_dfof(:,:,1) = data_dfof_max_ref;
    rgb_reg2ref_dfof(:,:,2) = reg2ref_dfof;
    figure; imagesc(rgb_reg2ref); title(['day 2 on day 1, data reg rgb overlay'])
%     mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['RGB Overlays']))
%     print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['RGB Overlays'], [date '_' mouse '_' run_str '_registered_reg2ref_overlay.pdf']), '-dpdf','-bestfit')
    figure; imagesc(rgb_reg2ref_dfof); title(['day 2 on day 1, dfof rgb overlay'])
%     print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['RGB Overlays'], [date '_' mouse '_' run_str '_registered_dfof_overlay.pdf']), '-dpdf','-bestfit')

%% Plotting
%     the first figure I showed
    figure; 
    subplot(2,3,1); imagesc(ref); axis off; axis equal; title('day 1 data avg')
    subplot(2,3,2); imagesc(reg); axis off; axis equal; title('day 2 data avg')
    subplot(2,3,3); imagesc(reg2ref); axis off; axis equal; title('registered day 2 data avg')
    subplot(2,3,4); imagesc(rgb_reg2ref_dfof(:,:,1)); axis off; axis equal; title('day 1 data dfof')
    subplot(2,3,5); imagesc(data_dfof_max_reg); axis off; axis equal; title('day 2 data dfof')
    subplot(2,3,6); imagesc(rgb_reg2ref_dfof(:,:,2)); axis off; axis equal; title('registered day 2 data dfof')
%     print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['RGB Overlays'], [date '_' mouse '_' run_str '_rgb_individuals.pdf']), '-dpdf','-bestfit')
