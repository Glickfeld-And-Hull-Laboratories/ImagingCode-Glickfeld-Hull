%% Introduction
% This script is designed to match two 2-photon datasets together. It uses
% one dataset as the reference and matches the other dataset to that. This
% was originally designed for experiments in which there are two
% classes of cells, "red" and "green", where red - typically some class of
% interneurons - are expressing a red fluorescence protein, and all cells
% are expressing GCaMP. This also assumes you have two runs on each day:
% one with your visual or behavioral stimuli, and one that is just used to
% tell which cells are red, refered to as a red snapshot.

%Throughout the script, parts that need user input are indiated with a
%comment that begins in an arrod like this: -->

%% Set-up

%This chunk is sets the paths where data will be found and stored.
% It assumes that the user has all experiment information organized in a
% datasheet ("ds"). DART_V1_YM90K_Celine is an exampleof a datasheet, ie a
% structure where there one entry for each experiment, identified by
% the experiment number.Each experiment has multiple fields for various
% relevant information.
clear all; clear global; close all
clc
% --> Edit this to the correct datasheet name
ds = 'DART_V1_YM90K_Celine'; 

% --> Edit this is your stimulus run field is called something different
dataStructLabels = {'contrastxori'};% identifies the stimulus runswithin your datasheet
eval(ds)

%these lines set up some basic parameters for the subsequent processing
doGreenOnly = false;
doCorrImg = true; 
brightnessScaleFactor = 0.3;

%This script assumes that you will match two days; one will be there
%reference, and the other will get matched to that. This is where you enter
%the reference day's experiment number from the datasheet. You only need to
%enter the experiment number for the reference day, in the following line
%the script automatically find that matched day based on what is entered in
%the datasheet.

% --> Edit this for the correct experiment ID
day_id(1) = 26; 
day_id(2) = expt(day_id(1)).multiday_matchdays;

%For the DART projects we have many different sub-projects that we have
%split into different experiment folders. If your datasheet has a similar
%strucutre, that will be identified here.
ExperimentFolder = expt(day_id(1)).exptType;

%this identifies how many days are going to be matched - typically 2.
nd = length(day_id);

%this identifies the mouse ID number
mouse = expt(day_id(1)).mouse;

%this is to find the correct path to the analysis folder, depending on
%whether the script is being run on Hubel or on Nuke and depending on the
%time since the previous experiment/infusion - that information is used to
%make a meaningulf folder name. 

% --> This should be edited to point to the desired analysis folder for your project.
if computer == 'GLNXA64'
    isilonName =  '/home/cc735@dhe.duke.edu/GlickfeldLabShare';
    database = fullfile('/All_Staff/home/ACh/Data/2p_data');
    base = fullfile(['/All_Staff/home/ACh/Analysis/2p_analysis/' ExperimentFolder]);
    beh_prefix = strcat(isilonName,'/All_Staff/Behavior/Data/data-');
else
    isilonName = 'duhs-user-nc1.dhe.duke.edu/';
    base = fullfile(['/home/ACh/Analysis/2p_analysis/' ExperimentFolder]);
    database = fullfile('/home/ACh/Data/2p_data');
   beh_prefix = strcat('Z:\Behavior\Data\data-');
end
fnout = fullfile(base,mouse);
if expt(day_id(2)).multiday_timesincedrug_hours>0
    dart_str = [expt(day_id(2)).drug '_' num2str(expt(day_id(1)).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end

% now an output folder with that name is created.
fn_multi = fullfile(base,mouse,['multiday_' dart_str]);
mkdir(fn_multi)

% these initialize the various output strucutres that will be filled over
% the course of the script
data = cell(1,nd);
fov_avg = cell(1,nd);
fov_norm = cell(1,nd);
fov_red = cell(1,nd);
align_fov = cell(1,nd);
corrmap = cell(1,nd);
dfmax = cell(1,nd);
masks = cell(1,nd);
maskNP = cell(1,nd);
red_ind = cell(1,nd);
cellTCs_all = cell(1,nd);

days_text = strcat('Reference day: ', string(day_id(1)), ' matched day: ', string(day_id(2)))
cd(fn_multi)
fid = fopen('sessionsMatched.txt','wt');
fprintf(fid, days_text);
fclose(fid);
%% load all data 

% This loops through the days and loads / registers each set of data.
runFolder = [];
for id = 1:nd % this 1:nd motif will come up many times; it's how we loop through the days
    clear global
    % find the informtion need to identify the corredt folder for the 2p
    % data
    expDate = expt(day_id(id)).date;
    runs = eval(['expt(day_id(' num2str(id) ')).' cell2mat(dataStructLabels) '_runs']);
    nrun = length(runs);
    out_all = [];
    data_g = [];

    for irun = 1:nrun
        imgFolder = runs{irun};
        fName = [imgFolder '_000_000'];
        cd(fullfile(database, mouse,expDate, imgFolder))
        %navigate to the folder that contains the 2p data for this day and
        %load the information for that run
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
        %sbxread is the scanbox function to read in the 2p data
        data_temp = sbxread(fName(1,1:11),0,nframes);
        %data_g will be the 2p data for the stimulus run for this
        %experiment
        data_g = cat(3, data_g, squeeze(data_temp(1,:,:,:)));
        clear data_temp
        out_all = [out_all; outs];
    end
    % this uses a GPU to register the data (that is, correct for X/Y
    % motion) by breaking the full data into smaller chunks. Assuming your
    % data has previously been reigstered, this will use the "outs"
    % generated in that prior registration and re-apply them to the data.
    [~,data{id}] = stackRegGPU(data_g,[],[],double(out_all));
    data_avg = mean(data{id},3);
    figure; imagesc(data_avg); title(['Avg FOV day ' num2str(id)])
    clear data_g
    % Not the {id} format here. This will be used repeatedly. We are using
    % a cell array with one cell per day, and then storing something - in
    % this case, the registered field of view averaged over frames - within
    % the cell for each day. 
    fov_avg{id} = data_avg;
    fov_norm{id} = uint8((fov_avg{id}./max(fov_avg{id}(:))).*255);
    fov_norm{id}(fov_norm{id} > (brightnessScaleFactor*255)) = brightnessScaleFactor*255;
    
    % Now we're checking whether there is a red snapshot for this day, and
    % loading it if there is.
    if exist(fullfile(fnout,expDate,runFolder,'redImage.mat'))
        load(fullfile(fnout,expDate,runFolder,'redImage.mat'))
    	fov_red{id} = uint8((redChImg./max(redChImg(:))).*255);
    else
        fov_red{id} = zeros(size(data_avg));
    end
    load(fullfile(fnout,expDate,runFolder,'mask_cell.mat'))

    % This next set of lines saves all the data for the individual days
    % into a cell array for the two days. 
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

% Here, we're making a correlation image for each day. This will be a
% version of the field of view image where the brighteness of r given pixel 
% corresponds to the magnitude of correlation between that pixel and those
% surrounding it. The correlation here is how the fluorescence of that
% pixel changes over time/frames compared to its neighbors. The concept is
% that pixels that together form a single cell should have highly
% correlated changes in fluorescence. The correlation image is one of three
% that will be used for the manual algiment.
corrmap_norm{1} = uint8((corrmap{1}./max(corrmap{1}(:))).*255);
corrmap_norm{2} = uint8((corrmap{2}./max(corrmap{2}(:))).*255);

% This creates a contrast-normalized version of the average FOV, another of
% the images that will be used for manual alignment.
fov_norm{1} = uint8((fov_avg{1}./max(fov_avg{1}(:))).*255);
fov_norm{2} = uint8((fov_avg{2}./max(fov_avg{2}(:))).*255);


%% manual align

% This chunk will allow you to make coarse adjustments to the field of
% view, accounting for rotation or variation in the X and Y dimensions.

% --> This will pop up a window for user input. The window will have four
% panels: a low zoom and a high zoom image of the FOV for each day. Your
% job is to find points that are the same across the two days; for
% example, you may see the same corner of a blood vessel or that same part
% of a given cell. You will flag that part of the image in the FOV for each
% day to tell the script "these points are the same across days." When
% you're done, the script will find the way to rotate the matched day that
% optimally lines up all those points to the reference day. 

% You will choose points on three images: the red snapshot, the correlation
% image and the normalized FOV. Picking 3-7 points on each image is usually
% enough, and points that are farther to the sides are more helpful than
% points in the center of the FOV. Be detailed - eg, you don't just want to
% find the same cell, you want to find the same part of an axon coming off
% that cell. When you're done with each image hit control + W to close it
% and the next image will come up.

if exist(fullfile(fn_multi,'multiday_alignment.mat'))
    % This checks whether there is already an aligmentment that was done
    % previously.If so, the script will load that rather than making the user
    % do it again.
    load(fullfile(fn_multi,'multiday_alignment.mat'))
else
    %If there is not a previously existing aligment, this is where you will
    %do it. 
    [input_points_1, base_points_1] = cpselect(fov_red{2},fov_red{1},'Wait', true);
    [input_points_2, base_points_2] = cpselect(fov_norm{2},fov_norm{1},'Wait', true);
    [input_points_3, base_points_3] = cpselect(corrmap_norm{2},corrmap_norm{1},'Wait', true);
    input_points = [input_points_1; input_points_2; input_points_3];
    base_points = [base_points_1; base_points_2; base_points_3];
end
%% Apply manual alignment
% This is where the script uses all the points you identified and actually
% rotates the matched FOV to align with the reference FOV. Then, it will
% pop out some images to show you the two FOVs and how the matched one was
% rotated. You'll most likely see some empty black space along the edges of
% the matched FOV where parts got cut off due to the rotation. 

fitGeoTAf = fitgeotrans(input_points(:,:), base_points(:,:),'affine'); 
data{3} = imwarp(data{2},fitGeoTAf, 'OutputView', imref2d(size(data{2}))); 
fov_avg{3} = mean(data{3},3);
fov_norm{3} = uint8((fov_avg{3}./max(fov_avg{3}(:))).*255);
corrmap{3} = double(imwarp(corrmap{2},fitGeoTAf, 'OutputView', imref2d(size(corrmap{2}))));
corrmap_norm{3} = uint8((corrmap{3}./max(corrmap{3}(:))).*255);
red_trans = (imwarp(fov_red{2},fitGeoTAf, 'OutputView', imref2d(size(fov_red{2}))));
fov_red{3} = uint8(red_trans);
dfmax{3} = (imwarp(dfmax{2},fitGeoTAf, 'OutputView', imref2d(size(dfmax{2}))));

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
%suptitle(mouse)
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

%% cell-by-cell correlation
% size of cell box
close all
w=30;
h=30;
buf = 3;
np = 5;
% green channel
% get cell centroids

cellImageAlign = struct;

cellPosition = regionprops(masks{1});
nc = length(cellPosition);

xCenter = cellfun(@(a) round(a(1)),{cellPosition.Centroid});
yCenter = cellfun(@(a) round(a(2)),{cellPosition.Centroid});

% index cells NOT too close to edge and NOT in black part of transformation
[ypix,xpix] = size(fov_avg{1});
goodCells = xCenter>(w/2) & xCenter<xpix-(w/2) & yCenter>(h/2) & yCenter<ypix-(h/2);

goodCells = goodCells & ...
    arrayfun(@(x) sum(sum(fov_norm{2}(masks{1}==x)))>0,1:nc);

    % fine register each cell    
mask_exp = zeros(size(fov_avg{1}));
mask_all = zeros(size(fov_avg{1}));

threshPercentile = 99;

thresholdedRed=cell(1,3);
for i=1:3
%     currentRed=fov_red{i};
%     highValues = find(currentRed>prctile(currentRed,threshPercentile,'all'));
%     redThresh = currentRed;
%     redThresh(highValues)=prctile(currentRed,threshPercentile,'all');
%     thresholdedRed{i}=redThresh;
    thresholdedRed{i}=fov_red{i};
end


            
start = 1;
for icell = 1:nc
    if goodCells(icell)
        % find best shift
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
        if red_ind{1}(icell)
            day1_red_avg = thresholdedRed{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
            day2_red_avg = thresholdedRed{3}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
            [red_reg_avg, shift_red] = shift_opt(double(day2_red_avg),double(day1_red_avg),2);
            r_red = corr(red_reg_avg(:),double(day1_red_avg(:)));
        else
            day1_red_avg = nan;
            day2_red_avg = nan;
            red_reg_avg = nan;
            r_red = nan;
        end
        
        [max_val max_ind] = max([r_avg r_max r_corr r_red]);
        if max_val>0.1 & (r_corr>0.1 || r_red>0.1 || r_max>0.1)
            pass = true;
            figure;
            movegui('center')
            start = 1;
            subplot(3,2,start)
            imagesc(day1_cell_corr)
            title('Corr')
            subplot(3,2,start+1)
            imagesc(reg_corr)
            title(num2str(r_corr))
            subplot(3,2,start+2)
            if red_ind{1}(icell)
                imagesc(day1_red_avg)
                 title('Red')
            else
                imagesc(day1_cell_avg)
                title('Avg')
            end
            subplot(3,2,start+3)
            if red_ind{1}(icell)
                imagesc(red_reg_avg)
                title(num2str(r_red))
            else
                imagesc(reg_avg)
                title(num2str(r_avg))
            end
            subplot(3,2,start+4)
            imagesc(day1_cell_max)
            title('Max')
            subplot(3,2,start+5)
            imagesc(reg_max)
            title(num2str(r_max))
            drawnow
         
            prompt = 'Choose image: 1- Corr, 2- Avg/Red, 3- Max, 0- skip: ';
            x = input(prompt);
            switch x
                case 0
                    pass = false;
                    shifts = nan;
                case 1
                    img_select = corrmap{3};
                    shifts = shift_corr;
                case 2
                    if red_ind{1}(icell)
                        img_select = thresholdedRed{3};
                        shifts = shift_red;
                    else
                        img_select = fov_avg{3};
                        shifts = shift_avg;
                    end
                case 3     
                    img_select = dfmax{3};
                    shifts = shift_max;
            end
        else
            pass = false;
            shifts = nan;
        end

        % shift data, get tc, all days
        if pass
            mask_data_temp = img_select;
            mask_data_temp(find(mask_exp >= 1)) = 0;
            mask_data_square = mask_data_temp(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
            movegui('center');
            bwout = imCellEditInteractive(mask_data_square);
            
            if sum(bwout(:))>1 %in case you chose not to select anything
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
    if length(find([cellImageAlign.pass])) ~= max(max(bwlabel(mask_all)))
        mask_all = mask_all-bwout_full;
        mask_exp = imCellBuffer(mask_all,3)+mask_all;
        error(['Mismatch in cell numbers- redo cell ' num2str(icell)])
    end
end
mask_cell = bwlabel(mask_all);
masks{id} = mask_cell;
mask_np = imCellNeuropil(mask_cell, 3, 5);
maskNP{id} = mask_np;
figure; movegui('center')
subplot(2,2,1)
imagesc(masks{1})
title('Day 1 masks')
subplot(2,2,2)
imagesc(mask_cell)
title('Day 2 masks after transform')
print(fullfile(fn_multi,'masksAfterTransform.pdf'),'-dpdf','-fillpage')

%%
%old TCs


match_ind = find([cellImageAlign.pass]);
cellTCs_match{1} = cellTCs_all{1}(:,match_ind);

%new TCs
data_tc = stackGetTimeCourses(data{3}, mask_cell);
[nFrames nCells] = size(data_tc);
down = 5;
data_reg_down  = stackGroupProject(data{3},down);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);

np_tc = zeros(nFrames,nCells);
np_tc_down = zeros(floor(nFrames./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data{3},mask_np(:,:,i));
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

cellTCs_match{2} = npSub_tc;

red_ind_match = ismember(match_ind,find(~isnan([cellImageAlign.r_red])));
red_ind_all = red_ind;

save(fullfile(fn_multi,'timecourses.mat'),'cellTCs_match', 'cellTCs_all', 'red_ind_all','red_ind_match','match_ind')
save(fullfile(fn_multi,'multiday_alignment.mat'),'cellImageAlign','fitGeoTAf', 'input_points','base_points', 'fov_avg', 'fov_norm','fov_red','dfmax','corrmap','masks','mask_np');

clear data_reg_down data