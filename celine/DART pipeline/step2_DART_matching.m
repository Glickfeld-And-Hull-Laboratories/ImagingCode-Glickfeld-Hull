%% Introduction
% This script matches two 2-photon datasets by registering one day to a reference day.
% Designed for experiments with "red" (interneurons with red fluorescent protein) and 
% "green" (all cells with GCaMP) cell populations. Assumes two runs per day:
% 1) Stimulus/behavior run 2) Red snapshot run (to identify red cells)

% User input sections are marked with: -->

%% Set-up

% Configure paths and load experiment datasheet
clear all; clear global; close all
clc

prompt = 'Enter name of instructions file: ';
instr = input(prompt, 's');
clear prompt
eval(instr);

ds=instructions.ds;
eval(ds);

rc = behavConstsDART; %this sets directory pathes based on the user's netID

% --> Edit if your stimulus run field has a different name
dataStructLabels = {'contrastxori'}; % identifies stimulus runs in datasheet


% Basic processing parameters
doCorrImg = true; 
brightnessScaleFactor = 0.3;

% Define days to match: day_id(1) = reference, day_id(2) = matched day
day_id(1) = str2double(instructions.session);

if length(expt) < day_id
    error('Day_id %d not valid for this dataset', day_id(1));
else
    day_id(2) = expt(day_id).multiday_matchdays;% automatically finds matched day
end

% Extract experiment info from datasheet
ExperimentFolder = expt(day_id(1)).exptType;
nd = length(day_id); % number of days (typically 2)
mouse = expt(day_id(1)).mouse;


% Create output folder with drug condition naming
fnout = fullfile(rc.analysis, ExperimentFolder,mouse);

if expt(day_id(2)).multiday_timesincedrug_hours>0
    dart_str = [expt(day_id(2)).drug '_' num2str(expt(day_id(1)).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end
fn_multi = fullfile(rc.analysis,mouse,['multiday_' dart_str]);
mkdir(fn_multi)

% Initialize storage structures for both days
data = cell(1,nd);           % raw registered data
fov_avg = cell(1,nd);        % average field of view
fov_norm = cell(1,nd);       % brightness-normalized FOV
fov_red = cell(1,nd);        % red channel images
align_fov = cell(1,nd);      % aligned FOV images
corrmap = cell(1,nd);        % pixel correlation maps
dfmax = cell(1,nd);          % max delta F/F images
masks = cell(1,nd);          % cell masks
maskNP = cell(1,nd);         % neuropil masks
red_ind = cell(1,nd);        % red cell indicators
cellTCs_all = cell(1,nd);    % all cell timecourses

% Log which sessions are being matched
days_text = strcat('Reference day: ', string(day_id(1)), ' matched day: ', string(day_id(2)))
cd(fn_multi)
fid = fopen('sessionsMatched.txt','wt');
fprintf(fid, days_text);
fclose(fid);

%% Load all data 

% Load and register data for each day
runFolder = [];
for id = 1:nd
    clear global
    
    % Get experiment info
    expDate = expt(day_id(id)).date;
    runs = eval(['expt(day_id(' num2str(id) ')).' cell2mat(dataStructLabels) '_runs']);
    nrun = length(runs);
    out_all = [];
    data_g = [];

    % Load data from all runs for this day
    for irun = 1:nrun
        imgFolder = runs{irun};
        fName = [imgFolder '_000_000'];
        cd(fullfile(rc.data, mouse, expDate, imgFolder))
        
        load(fName)
        if nrun == 1
            load(fullfile(fnout,expDate,imgFolder,'regOuts&Img.mat'))
            nframes = size(outs,1);
            runFolder = imgFolder;
        else
            nframes = info.config.nframes;
            runFolder = [runFolder '_' imgFolder];
        end
        
        fprintf(['Loading day ' num2str(id) ' data \n Mouse: ' expt(day_id(id)).mouse ' Date: ' expt(day_id(id)).date '\n'])
        
        % Read 2-photon data using scanbox function
        data_temp = sbxread(fName(1,1:11),0,nframes);
        data_g = cat(3, data_g, squeeze(data_temp(1,:,:,:)));
        clear data_temp
        out_all = [out_all; outs];
    end
    
    % Register data using GPU (corrects for X/Y motion using previously computed shifts)
    [~,data{id}] = stackRegGPU(data_g,[],[],double(out_all));
    data_avg = mean(data{id},3);
    figure; imagesc(data_avg); title(['Avg FOV day ' num2str(id)])
    clear data_g
    
    % Store processed images for this day
    fov_avg{id} = data_avg;
    fov_norm{id} = uint8((fov_avg{id}./max(fov_avg{id}(:))).*255);
    fov_norm{id}(fov_norm{id} > (brightnessScaleFactor*255)) = brightnessScaleFactor*255;
    
    % Load red channel snapshot if available
    if exist(fullfile(fnout,expDate,runFolder,'redImage.mat'))
        load(fullfile(fnout,expDate,runFolder,'redImage.mat'))
    	fov_red{id} = uint8((redChImg./max(redChImg(:))).*255);
    else
        fov_red{id} = zeros(size(data_avg));
    end
    
    % Load previously computed masks and timecourses
    load(fullfile(fnout,expDate,runFolder,'mask_cell.mat'))
    dfmax{id} = data_dfof(:,:,1);
    corrmap{id} = data_dfof(:,:,end);
    masks{id} = mask_cell;
    maskNP{id} = mask_np;
    red_ind{id} = mask_label;
    
    load(fullfile(fnout,expDate,runFolder,'TCs.mat'))
    cellTCs_all{id} = npSub_tc;
    
    load(fullfile(fnout,expDate,runFolder,'input.mat'))
    input_temp(id) = input;
end

input = input_temp;
save(fullfile(fn_multi,'input.mat'),'input')
clear input

%% Make correlation image and normalized FOV for manual alignment

% Create correlation maps: pixel brightness = correlation with neighboring pixels
% High correlation indicates pixels belonging to the same cell
corrmap_norm{1} = uint8((corrmap{1}./max(corrmap{1}(:))).*255);
corrmap_norm{2} = uint8((corrmap{2}./max(corrmap{2}(:))).*255);

% Create contrast-normalized FOV images
fov_norm{1} = uint8((fov_avg{1}./max(fov_avg{1}(:))).*255);
fov_norm{2} = uint8((fov_avg{2}./max(fov_avg{2}(:))).*255);

%% Manual alignment

% Manual coarse alignment to account for rotation and X/Y shifts between days

% --> INTERACTIVE STEP: A 4-panel window will appear showing low/high zoom views
% of both days. Select matching anatomical landmarks (blood vessel corners, 
% cell features, etc.) across days. Choose 3-7 points per image, prioritizing
% points toward the edges over center. Be precise - select specific parts of 
% features, not just general areas.
% 
% You'll select points on three image types:
% 1) Red snapshot 2) Correlation map 3) Normalized FOV
% Press Ctrl+W when done with each image to proceed to the next.

if exist(fullfile(fn_multi,'multiday_alignment.mat'))
    % Load existing alignment if previously computed
    load(fullfile(fn_multi,'multiday_alignment.mat'))
else
    % Perform new manual alignment
    [input_points_1, base_points_1] = cpselect(fov_red{2},fov_red{1},'Wait', true);
    [input_points_2, base_points_2] = cpselect(fov_norm{2},fov_norm{1},'Wait', true);
    [input_points_3, base_points_3] = cpselect(corrmap_norm{2},corrmap_norm{1},'Wait', true);
    input_points = [input_points_1; input_points_2; input_points_3];
    base_points = [base_points_1; base_points_2; base_points_3];
end

%% Apply manual alignment

% Transform matched day FOV to align with reference day using selected points
fitGeoTAf = fitgeotrans(input_points(:,:), base_points(:,:),'affine'); 
data{3} = imwarp(data{2},fitGeoTAf, 'OutputView', imref2d(size(data{2}))); 
fov_avg{3} = mean(data{3},3);
fov_norm{3} = uint8((fov_avg{3}./max(fov_avg{3}(:))).*255);
corrmap{3} = double(imwarp(corrmap{2},fitGeoTAf, 'OutputView', imref2d(size(corrmap{2}))));
corrmap_norm{3} = uint8((corrmap{3}./max(corrmap{3}(:))).*255);
red_trans = (imwarp(fov_red{2},fitGeoTAf, 'OutputView', imref2d(size(fov_red{2}))));
fov_red{3} = uint8(red_trans);
dfmax{3} = (imwarp(dfmax{2},fitGeoTAf, 'OutputView', imref2d(size(dfmax{2}))));

% Display alignment results
figure;colormap gray
movegui('center')
subplot 221
imshow(fov_norm{1}); title('Day 1 Data Avg')
subplot 222
imshow(fov_norm{3}); title('Transformed Day 2 Data Avg')
subplot 223
filler = zeros(size(fov_norm{1}));
imshow(cat(3,fov_norm{1},fov_norm{3},filler))
title('Overlay')
print(fullfile(fn_multi,'FOV manual alignment'),'-dpdf','-fillpage')

figure;colormap gray
movegui('center')
d1 = uint8((corrmap{1}./max(corrmap{1}(:))).*255);
d2 = uint8((corrmap{2}./max(corrmap{2}(:))).*255);
d21 = uint8((corrmap{3}./max(corrmap{3}(:))).*255);
subplot 221
imshow(d1); title('Day 1 Pixel Correlation Map')
subplot 222
imshow(d21); title('Transformed Day 2 Pixel Correlation Map')
subplot 223
imshow(cat(3,d1,d21,filler))
title('Overlay')
print(fullfile(fn_multi,'FOV correlation map manual alignment'),'-dpdf','-fillpage')

figure;colormap gray
movegui('center')
subplot 221
imshow(fov_red{1}); title('Day 1 Data Avg')
subplot 222
imshow(fov_red{3}); title('Transformed Day 2 Data Avg')
subplot 223
filler = zeros(size(fov_red{1}));
imshow(cat(3,fov_red{1},fov_red{3},filler))
title('Overlay')
print(fullfile(fn_multi,'red channel manual alignment'),'-dpdf','-fillpage')

%% Cell-by-cell matching

% --> INTERACTIVE STEP: Individual cell matching with user input
% For each "good cell" (not cut off by transformation), you'll see 6 images:
% Left column: Reference day cell in 3 views (correlation, avg/red, max df/f)
% Right column: Corresponding matched day views with correlation values
% 
% The target cell is always in the center of each image patch.
% Select which right-side view best shows the matching cell, then manually
% segment that cell. Only one cell per patch allowed.

% Parameters for cell patches
close all
w=30; h=30; % patch dimensions
buf = 3; np = 5;

% Get cell positions from reference day
cellPosition = regionprops(masks{1});
nc = length(cellPosition);
xCenter = cellfun(@(a) round(a(1)),{cellPosition.Centroid});
yCenter = cellfun(@(a) round(a(2)),{cellPosition.Centroid});

% Identify cells not too close to edges and not in transformation black regions
[ypix,xpix] = size(fov_avg{1});
goodCells = xCenter>(w/2) & xCenter<xpix-(w/2) & yCenter>(h/2) & yCenter<ypix-(h/2);
goodCells = goodCells & ...
    arrayfun(@(x) sum(sum(fov_norm{2}(masks{1}==x)))>0,1:nc);

% Initialize cell matching
mask_exp = zeros(size(fov_avg{1}));
mask_all = zeros(size(fov_avg{1}));

% Loop through each good cell for matching
start = 1;
for icell = 1:nc
    if goodCells(icell)
        % Extract cell patches and find optimal alignment for each view
        day1_cell_avg = fov_avg{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        day2_cell_avg = fov_avg{3}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        [reg_avg, shift_avg] = shift_opt(day2_cell_avg,day1_cell_avg,2);
        r_avg = corr(reg_avg(:),day1_cell_avg(:));
        
        day1_cell_max = dfmax{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        day2_cell_max = dfmax{3}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        [reg_max, shift_max] = shift_opt(day2_cell_max,day1_cell_max,2);
        r_max = corr(reg_max(:),day1_cell_max(:));

        day1_cell_corr = corrmap{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        day2_cell_corr = corrmap{3}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        [reg_corr, shift_corr] = shift_opt(day2_cell_corr,day1_cell_corr,2);
        r_corr = corr(reg_corr(:),day1_cell_corr(:));
        
        % Handle red channel for red cells
        if red_ind{1}(icell)
            day1_red_avg = fov_red{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
            day2_red_avg = fov_red{3}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
            [red_reg_avg, shift_red] = shift_opt(double(day2_red_avg),double(day1_red_avg),2);
            r_red = corr(red_reg_avg(:),double(day1_red_avg(:)));
        else
            day1_red_avg = nan; day2_red_avg = nan; red_reg_avg = nan; r_red = nan;
        end
        
        % Show images for manual matching if correlation is sufficient
        [max_val max_ind] = max([r_avg r_max r_corr r_red]);
        if max_val>0.1 & (r_corr>0.1 || r_red>0.1 || r_max>0.1)
            pass = true;
            figure;
            movegui('center')
            start = 1;
            
            % Display 6-panel comparison
            subplot(3,2,start)
            imagesc(day1_cell_corr); title('Corr')
            subplot(3,2,start+1)
            imagesc(reg_corr); title(num2str(r_corr))
            subplot(3,2,start+2)
            if red_ind{1}(icell)
                imagesc(day1_red_avg); title('Red')
            else
                imagesc(day1_cell_avg); title('Avg')
            end
            subplot(3,2,start+3)
            if red_ind{1}(icell)
                imagesc(red_reg_avg); title(num2str(r_red))
            else
                imagesc(reg_avg); title(num2str(r_avg))
            end
            subplot(3,2,start+4)
            imagesc(day1_cell_max); title('Max')
            subplot(3,2,start+5)
            imagesc(reg_max); title(num2str(r_max))
            drawnow
         
            % Get user input for image selection
            prompt = 'Choose image: 1- Corr, 2- Avg/Red, 3- Max, 0- skip: ';
            x = input(prompt);
            switch x
                case 0
                    pass = false; shifts = nan;
                case 1
                    img_select = corrmap{3}; shifts = shift_corr;
                case 2
                    if red_ind{1}(icell)
                        img_select = fov_red{3}; shifts = shift_red;
                    else
                        img_select = fov_avg{3}; shifts = shift_avg;
                    end
                case 3     
                    img_select = dfmax{3}; shifts = shift_max;
            end
        else
            pass = false; shifts = nan;
        end

        % Manual cell segmentation if cell passes correlation threshold
        if pass
            mask_data_temp = img_select;
            mask_data_temp(find(mask_exp >= 1)) = 0; % avoid double-segmenting
            mask_data_square = mask_data_temp(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
            movegui('center');
            bwout = imCellEditInteractive(mask_data_square);
            
            if sum(bwout(:))>1 % check if cell was actually segmented
                bwout_full = zeros(size(fov_avg{1}));
                bwout_full(...
                yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
                xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1) = bwout;
                mask_all = mask_all+bwout_full;
                mask_exp = imCellBuffer(mask_all,3)+mask_all;
            else
                pass = false;
            end
        end
        
        close all
        
        % Store matching results
        cellImageAlign(icell).center_yx = [yCenter(icell),xCenter(icell)];
        cellImageAlign(icell).d(1).avg_img = day1_cell_avg;
        cellImageAlign(icell).d(1).corr_img = day1_cell_corr;
        cellImageAlign(icell).d(1).red_img = day1_red_avg;
        cellImageAlign(icell).d(1).max_img = day1_cell_max;
        cellImageAlign(icell).d(2).avg_img = reg_avg;
        cellImageAlign(icell).d(2).corr_img = reg_corr;
        cellImageAlign(icell).d(2).red_img = red_reg_avg;
        cellImageAlign(icell).d(2).max_img = reg_max;
        cellImageAlign(icell).r_avg = r_avg;
        cellImageAlign(icell).r_corr = r_corr;
        cellImageAlign(icell).r_red = r_red;
        cellImageAlign(icell).shifts = shifts;
        cellImageAlign(icell).pass = pass;
    else
        cellImageAlign(icell).pass = false;
        cellImageAlign(icell).r_red = 0;
    end
    
    % Error checking for mask consistency
    if length(find([cellImageAlign.pass])) ~= max(max(bwlabel(mask_all)))
        mask_all = mask_all-bwout_full;
        mask_exp = imCellBuffer(mask_all,3)+mask_all;
        error(['Mismatch in cell numbers- redo cell ' num2str(icell)])
    end
end

% Finalize matched masks
mask_cell = bwlabel(mask_all);
masks{id} = mask_cell;
mask_np = imCellNeuropil(mask_cell, 3, 5);
maskNP{id} = mask_np;

% Display final mask comparison
figure; movegui('center')
subplot(2,2,1)
imagesc(masks{1}); title('Day 1 masks')
subplot(2,2,2)
imagesc(mask_cell); title('Day 2 masks after transform')
print(fullfile(fn_multi,'masksAfterTransform.pdf'),'-dpdf','-fillpage')

%% Extract timecourses from matched cells

% Extract fluorescence timecourses from successfully matched cells
% Day 1: use existing timecourses from previous analysis
% Day 2: extract new timecourses from newly segmented masks with neuropil subtraction

match_ind = find([cellImageAlign.pass]);
cellTCs_match{1} = cellTCs_all{1}(:,match_ind);

% Extract new timecourses from transformed Day 2 data
data_tc = stackGetTimeCourses(data{3}, mask_cell);
[nFrames nCells] = size(data_tc);

% Downsample for neuropil weight optimization
down = 5;
data_reg_down  = stackGroupProject(data{3},down);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);

% Calculate neuropil timecourses
np_tc = zeros(nFrames,nCells);
np_tc_down = zeros(floor(nFrames./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data{3},mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '%s/n']) 
end

% Optimize neuropil subtraction weights by maximizing skewness
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
[max_skew ind] =  max(x,[],1);
np_w = 0.01*ind;

% Apply neuropil subtraction
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
cellTCs_match{2} = npSub_tc;

% Track red cell indicators
red_ind_match = ismember(match_ind,find(~isnan([cellImageAlign.r_red])));
red_ind_all = red_ind;

% Save results
save(fullfile(fn_multi,'timecourses.mat'),'cellTCs_match', 'cellTCs_all', 'red_ind_all','red_ind_match','match_ind')
save(fullfile(fn_multi,'multiday_alignment.mat'),'cellImageAlign','fitGeoTAf', 'input_points','base_points', 'fov_avg', 'fov_norm','fov_red','dfmax','corrmap','masks','mask_np');

clear data_reg_down data