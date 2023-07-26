%clear everything
clear all
clear all global
clc
close all

%%
tj_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\2P_Imaging';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P';
behav_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
dataset = 'exp_list_phaserev_tjw'; %experiment list to pick files from
eval(dataset); %load dataset
%%
%session 2 info
expt_id = 14;
date = expt(expt_id).date;
ImgFolder = expt(expt_id).runs; %could we use char() instead here?
time = expt(expt_id).time_mat;
mouse = expt(expt_id).mouse;
run = ImgFolder; %multiple depths?***
nrun = size(ImgFolder,1); %what is this?***
run_str = catRunName(ImgFolder, nrun);
%%
%session 1 info
first_expt_id = 8;
ref_date = expt(first_expt_id).date;
ref_str = ['runs-', expt(first_expt_id).runs];
%%
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun 
    CD = fullfile(tj_fn, [mouse '\' date '_' mouse '\' ImgFolder(irun,:)]); %identify current dierectory;
    cd(CD); %set CD
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat']; %add the 0s to the imaging file
    load(imgMatFile); %**load this file from CD
    fName = fullfile(behav_fn, ['data-'  mouse  '-' date '-' time(irun,:) '.mat']); %find behavior data
    load(fName); %load behavior data
    
    %input is behavioral parameters
    %info is imaging parameters

    nframes = info.config.frames; %find nframes from info
    fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n']) %graphic display of frames reading
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes); %reads data from 0 to nframes from raw file
    %nPMT x nYpix x nXpix x nframes
    
    
    temp(irun) = input; 
    if isfield(input, 'nScansOn') %checks that nScansOn is in input var
        nOn = temp(irun).nScansOn; %find nOn in input
        nOff = temp(irun).nScansOff; %find nOff in input
        ntrials = size(temp(irun).tGratingDirectionDeg,2); %find ntrials based on input

        data_temp = squeeze(data_temp); %use only 1 PMT so squeeze data
        if nframes>ntrials*(nOn+nOff) %make sure nframes matches
            nframes = ntrials*(nOn+nOff); %if not make it match
            data_temp = data_temp(:,:,1:ntrials*(nOn+nOff)); %will restructure it to nYpix x nXpix x nframes
        elseif nframes<ntrials*(nOn+nOff) %similar to above
            temp(irun) = trialChopper(temp(irun),1:ceil(nframes./(nOn+nOff))); %rounds up to number of frames
        end
    end
    
    offset = offset+nframes;

    data_temp = squeeze(data_temp); 
    data = cat(3,data,data_temp); %concatenate data and data_temp along 3rd dimension;
    trial_n = [trial_n nframes]; 
end
input = concatenateDataBlocks(temp); %combines mat files;
clear data_temp
%clear temp

 %% Choose register interval
t = 2000; %nframes to skip for each average; could add nframes to not hard code number to average
nep = floor(size(data,3)./t); %divides frames by skips and rounds down
[n n2] = subplotn(nep); %finds ideal number of subplots to make
figure; %makes figure
for i = 1:nep; %for the number of plots
    subplot(n,n2,i); %this subplot
    imagesc(mean(data(:,:,1+((i-1)*t):500+((i-1)*t)),3)); %scaled color image of mean of frames for specified range
    title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]); %titled based on frame numbers
end
%these figures are taking averages of each pixel value across certain sets
%of frames; ex: what is the avg pixel value of pixel 1,1 for these 500
%frames; what about pixel 1,2 etc.

%% Register data

data_avg = mean(data(:,:,18001:18500),3); %mean of pixel values over selected range of frames

if exist(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str])) %if this folder exists)
    load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat'])) %load this mat file
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input') %save input?
    [outs, data_reg]=stackRegister_MA(data(:,:,:),[],[],out); %using shifts data to move all frames to target image
%    clear out outs
else
    [out, data_reg] = stackRegister(data,data_avg); %stacks 3d frames data to 2d avg target
    data_reg_avg = mean(data_reg,3); %mean of all registered frames
    reg = data_reg_avg; %sets reg = to above
    mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str])) %make new directory and save
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'data_reg_avg', 'out', 'data_avg')
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
%clear data


%% test stability
figure; imagesq(data_reg_avg); truesize; %why not imagesc?; avg pixel value of all frames registered***
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf','-bestfit') %save as pdf that fits the page

%%
nOn = double(input.nScansOn);
nOff = double(input.nScansOff);
phaseCyc = double(input.nScansPhaseCyc);
dir_mat = celleqel2mat_padded(input.tGratingDirectionDeg);
nTrials = length(dir_mat);
sz = size(data_reg);

data_resp = reshape(data_reg, [sz(1) sz(2) nOn+nOff nTrials]);
data_f = mean(data_resp(:,:,nOff-15:nOff,:),3);
data_resp_dfof = (double(data_resp)-data_f)./data_f; %Ypix x Xpix x frames x trials

clear data_resp data_f

dirs = unique(dir_mat);
nDir = length(dirs);

nStim = double(nDir); %number of stimulus combination selected randomly from list of 90 w/ replacement

% [n n2] = subplotn(nStim);
% data_dfof_stim = nan(sz(1),sz(2),nDir);
% start = 1;
% for iDir = 1:nDir
%     ind = find(dir_mat == dirs(iDir));
%     data_dfof_stim(:,:,iDir) = mean(mean(data_resp_dfof(:,:,nOff+1:nOff+nOn,ind),4),3);
%     subplot(n,n2,start)
%     imagesc(data_dfof_stim(:,:,iDir))
%     start = start+1;
% end

[n n2] = subplotn(nStim);
data_dfof_stim = nan(sz(1),sz(2),nDir);
start = 1;
for iDir = 1:nDir
    ind = find(dir_mat == dirs(iDir));
    data_dfof_stim(:,:,iDir,1) = mean(mean(data_resp_dfof(:,:,nOff+1:nOff+input.nScansPhaseCyc,ind),4),3);
    subplot(n,n2,start)
    imagesc(data_dfof_stim(:,:,iDir,1))
    start = start+1;
%     data_dfof_stim(:,:,iDir,2) = mean(mean(data_resp_dfof(:,:,nOff+1+input.nScansPhaseCyc:nOff+input.nScansPhaseCyc.*2,ind),4),3);
%     subplot(n,n2,start)
%     imagesc(data_dfof_stim(:,:,iDir,2))
%     start = start+1;
end

    
data_dfof_stim_all = reshape(data_dfof_stim, [sz(1) sz(2) nDir]);
data_dfof_max = max(data_dfof_stim_all,[],3);
clear data_dfof_stim data_resp_dfof

figure; imagesc(data_dfof_max); movegui('center')

data_dfof = cat(3,data_dfof_max,data_dfof_stim_all);
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_segmentData.mat']), 'data_dfof_max', 'data_dfof_stim_all', 'nStim')


%%
%pixel correlation
data_reg_3hz = stackGroupProject(data_reg,5); %averaging every 5 frames for less noisy stack
pix = getPixelCorrelationImage(data_reg_3hz); %how much each pixel is correlated with surrounding? -> could save this figure**
pix(isnan(pix))=0;
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pixel.mat']),'pix')

clear data_reg_3hz

%% rename day2 variables and load files from ref date 
fov_avg{2} = data_reg_avg; %avg fov for all reg frames in a cell array format from D2
dfmax{2} = data_dfof_max; %max fov for all reg frames from D2
corrmap{2} = pix; %pix corr map from D2
corrmap_norm{2} = uint8((corrmap{2}./max(corrmap{2}(:))).*255);
brightnessScaleFactor = 0.5; 
fov_norm{2} = uint8((fov_avg{2}./max(fov_avg{2}(:))).*255); 
fov_norm{2}(fov_norm{2} > (brightnessScaleFactor*255)) = brightnessScaleFactor*255; 


maskD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_mask_cell.mat']));
TCs_D1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_TCs.mat']));
dfofD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_segmentData.mat']));
shiftsD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_reg_shifts.mat']));
pixelD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_pixel.mat']));
inputD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_input.mat']));

%Need alternative load info if matching previously done!!
fov_avg{1} = shiftsD1.data_reg_avg; %d1 fov across all frames
dfmax{1} = dfofD1.data_dfof_max; %d1 max
cellTCs_all{1} = TCs_D1.npSub_tc; %np-subtracted TC from d1
%cellTCs_all{1} = TCs_D1.cellTCs_match{2}; %np-sub TC from D1 when it is registered as D2
input_temp(1) = inputD1; %input (behavioral setup) from d1
corrmap{1} = pixelD1.pix; %pix corr from d1
masks{1} = maskD1.mask_cell; %cell mask from d1
corrmap_norm{1} = uint8((corrmap{1}./max(corrmap{1}(:))).*255); %normalizing again?
brightnessScaleFactor = 0.5;
fov_norm{1} = uint8((fov_avg{1}./max(fov_avg{1}(:))).*255);
fov_norm{1}(fov_norm{1} > (brightnessScaleFactor*255)) = brightnessScaleFactor*255;

%% align Day 2 data stack to Day 1
[input_points_1, base_points_1] = cpselect(fov_norm{2},fov_norm{1},'Wait', true); %select control points from 2 images (fov images); 2 is matched to 1
[input_points_2, base_points_2] = cpselect(dfmax{2},dfmax{1},'Wait', true); %select control points from max projection
input_points = [input_points_1' input_points_2']; %x,y pairings of points on d2
base_points = [base_points_1' base_points_2']; %x,y pairings of points on d1
fitGeoTAf = fitgeotrans(input_points(:,:)', base_points(:,:)','affine'); %fits your input points to the base points? ***

% transform data_reg and make a new D2 FOV from the avg of that
data3 = imwarp(data_reg,fitGeoTAf, 'OutputView', imref2d(size(data_reg))); %displacing d2 data reg by the fits above 
fov_avg{3} = mean(data3,3); %new fov
fov_norm{3} = uint8((fov_avg{3}./max(fov_avg{3}(:))).*255); %normalized?***
corrmap{3} = double(imwarp(corrmap{2},fitGeoTAf, 'OutputView', imref2d(size(corrmap{2})))); %new pix corr map
corrmap_norm{3} = uint8((corrmap{3}./max(corrmap{3}(:))).*255);
dfmax{3} = imwarp(dfmax{2},fitGeoTAf, 'OutputView', imref2d(size(dfmax{2}))); %new dfof max
%redChImg = imwarp(fov_red{2},fitGeoTAf, 'OutputView', imref2d(size(fov_red{2})));

%SAVE all of the this!!

fov_avg_shift = fov_avg{3};
fov_norm_shift = fov_norm{3};
corrmap_shift = corrmap{3};
corrmap_norm_shift = corrmap_norm{3};
dfmax_shift = dfmax{3};
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_rematch_shifts.mat']),'input_points', 'base_points', 'fitGeoTAf', 'fov_avg_shift', 'fov_norm_shift', 'corrmap_shift', 'corrmap_norm_shift', 'dfmax_shift')


sz = size(fov_avg{1});
rgb = zeros(sz(1),sz(2),3);
%rgb(:,:,1) = redChImg./max(redChImg(:));
rgb(:,:,2) = fov_avg{3}./max(fov_avg{3}(:));
figure; image(rgb); movegui('center') %i think this would overlay your red and green?*** -> no, it overlays d1 and d2 images

figure;colormap gray %not sure what this part is doing - why the 3 and filler?*** also why is it colored?***
filler = zeros(size(dfmax{1}));
movegui('center')
subplot 221
imshow(cat(3,fov_norm{1},fov_norm{3},filler)) %makes rgb images
subplot 222
%imshow(cat(3,redChImg, fov_norm{1},filler))
subplot 223
imshow(cat(3,dfmax{1},dfmax{3},filler))
title('Overlay')
%print(fullfile(fn_multi,'overlays'),'-dpdf','-fillpage')

%the 3 concatenates that dimension; colormap gray is wrong here
%might want to save the fitGeoTrans data***
%fov_avg{1} is the d1 data
%fov_avg{2} is the NON-shifted d2/3 data
%fov_avg{3} is the shifted d2/3 data

%% make new cell masks for D2
clear input %why? -> i think it is used below***
% size of cell box
close all
w=30; %width?***
h=30; %height?***
%buf = 3; %buffer
buf = 4; %buffer
np = 5; %neuropil

% green channel
% get cell centroids

cellImageAlign = struct; %create structure array

cellPosition = regionprops(masks{1}); %properties of image regions - results in each cell centroid defined along with an 'area'?***
nc = length(cellPosition);

xCenter = cellfun(@(a) round(a(1)),{cellPosition.Centroid}); %apply function to each cell of cell array - rounds x and y centroids
yCenter = cellfun(@(a) round(a(2)),{cellPosition.Centroid}); %what is a here?***

% index cells NOT too close to edge and NOT in black part of transformation
% - ?***
[ypix,xpix] = size(fov_avg{1}); %npixels
goodCells = xCenter>(w/2) & xCenter<xpix-(w/2) & yCenter>(h/2) & yCenter<ypix-(h/2); %gets rid of cells close to edge

goodCells = goodCells & ...
    arrayfun(@(x) sum(sum(fov_norm{2}(masks{1}==x)))>0,1:nc); %?***

    % fine register each cell    
mask_exp = zeros(size(fov_avg{1}));
mask_all = zeros(size(fov_avg{1}));
            
start = 1; 
figure; movegui('center');
for icell = 1:nc %for each cell
    if goodCells(icell) %if cell is 'good' - remember that all below is going to be for each cell
        %if exist('multiDayD1') end is line 495
            %if multiDayD1.cellImageAlign(icell).pass end is line 494
                % find best shift - ?***
                day1_mask = masks{1}(... %%need original D1 mask trans to D2 - ?***
                    yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,... %is this identifying each cell as a cluster of pixels in mask?***
                    xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
                day1_cell_avg = fov_avg{1}(... %doing same as above but with fov on d1
                    yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
                    xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
                day2_cell_avg = fov_avg{3}(... %same as above but with d3
                    yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
                    xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
                [reg_avg, shift_avg] = shift_opt(day2_cell_avg,day1_cell_avg,2); % finds optimal shift to align data and target using correlation 
                r_avg = corr(reg_avg(:),day1_cell_avg(:)); %corr between reg_avg and cell_avg ?***
                
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
                
        %         only use cells with with one corr coef of 0.4 or higher
                [max_val max_ind] = max([r_avg r_max r_corr]); %find max corr values?*** - check values with Celine
                if max_val>0.55 & (r_corr>0.4 || r_avg>0.4 || r_max>0.4)
                    pass = true;
                    figure;
                    movegui('center')
                    start = 1;
                    subplot(3,2,start)
                    imagesc(day1_cell_corr) %this is the interactive part where d1 and d2/3 cells are plotted side by side
                    title('Corr')
                    subplot(3,2,start+1)
                    imagesc(reg_corr)
                    title(num2str(r_corr))
                    subplot(3,2,start+2)
                    imagesc(day1_cell_avg)
                    title('Avg')
                    subplot(3,2,start+3)
                    imagesc(reg_avg)
                    title(num2str(r_avg))
                    subplot(3,2,start+4)
                    imagesc(day1_cell_max)
                    title('Max')
                    subplot(3,2,start+5)
                    imagesc(reg_max)
                    title(num2str(r_max))
                    suptitle(['Cell #' num2str(icell)])
                    prompt = 'Choose image: 1- Corr, 2- Avg/Red, 3- Max, 0- skip: ';
                    drawnow
                    x = input(prompt); %this is why we cleared input above?***
                    switch x %switch between cases
                        case 0 %if 0 then no shift made and pass is = to false
                            pass = false;
                            shifts = nan;
                        case 1 %if 1,2,or 3 then use shifts from above
                            img_select = corrmap{3};
                            shifts = shift_corr;
                        case 2
                            img_select = fov_avg{3};
                            shifts = shift_avg;
                        case 3     
                            img_select = dfmax{3};
                            shifts = shift_max;
                    end
                else
                    pass = false;
                    shifts = nan;
                end
        
                %***THIS IS WHERE ALL THE SINGLE CELL DATA ARE KEPT, INCLUDING IF THE CELLS MATCH ACROSS DAYS***
                cellImageAlign(icell).center_yx = [yCenter(icell),xCenter(icell)]; %cell centers
                cellImageAlign(icell).d(1).avg_img = day1_cell_avg; %this and the ones below are avg,corr,max for d1 and d2
                cellImageAlign(icell).d(1).corr_img = day1_cell_corr;
                cellImageAlign(icell).d(1).max_img = day1_cell_max;
                cellImageAlign(icell).d(2).avg_img = reg_avg;
                cellImageAlign(icell).d(2).corr_img = reg_corr;
                cellImageAlign(icell).d(2).max_img = reg_max;
                cellImageAlign(icell).shifts = shifts; %shifts?***
                cellImageAlign(icell).pass = pass; %did not press 0
        
        
                % shift data, get tc, all days
                if pass %i think this is where you are selecting cell pixels***
                    mask_data_temp = img_select; %depends on which image stack you picked
                    mask_data_temp(find(mask_exp >= 1)) = 0; %darken out cells on above image
                    mask_data_square = mask_data_temp(... %make small square for selection
                    yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
                    xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
                    bwout = imCellEditInteractive(mask_data_square);
                    if sum(bwout(:))>1 %in case you chose not to select anything
                        bwout_full = zeros(size(fov_avg{1}));
                        bwout_full(...
                        yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
                        xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1) = bwout*icell;
                        mask_all = mask_all+bwout_full;
                        mask_temp = mask_all;
                        mask_temp(find(mask_all>=1)) = 1;
                        mask_exp = imCellBuffer(mask_temp,buf)+mask_temp;
                    else
                        cellImageAlign(icell).pass = false;
                        corner1 = 2;
                        corner2 = 2;
                        temp_mask = zeros(size(day1_mask));
%                         temp_mask(15,15) = 1;
                        temp_mask(corner1,corner2) = 1;
                        bwout_full = zeros(size(fov_avg{1}));
                        bwout_full(...
                        yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
                        xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1) = temp_mask*icell;
                        mask_all = mask_all+bwout_full;
                        mask_temp = mask_all;
                        mask_temp(find(mask_all>=1)) = 1;
                        mask_exp = imCellBuffer(mask_temp,buf)+mask_temp;
                        corner1+1;
                        corner2+1;
                    end
                close all
                else
                cellImageAlign(icell).pass = false;
                temp_mask = zeros(size(day1_mask));
                temp_mask(15,15) = 1;
                bwout_full = zeros(size(fov_avg{1}));
                bwout_full(...
                yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
                xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1) = temp_mask*icell;
                mask_all = mask_all+bwout_full;
                mask_temp = mask_all;
                mask_temp(find(mask_all>=1)) = 1;
                mask_exp = imCellBuffer(mask_temp,buf)+mask_temp;
                end
            %end
        %end
    else
        %cells are placeholders despite not being in new FOV
        a = yCenter(icell)-(h/2); %here to line 446 for cells on edge - ?***
        b = yCenter(icell)+(h/2)-1;
        c = xCenter(icell)-(w/2);
        d = xCenter(icell)+(w/2)-1;
        if a < 1
            a = 1;
        end
        if c < 1
            c = 1;
        end
        if b > size(fov_avg{3},1)
            b = size(fov_avg{3},1);
        end
        if d > size(fov_avg{3},2)
            d = size(fov_avg{3},2);
        end
        day1_mask = masks{1}(... %%need original D1 mask trans to D2
        a:b,...
        c:d);
        cellImageAlign(icell).pass = false;
        cellImageAlign(icell).r_red = 0;
        temp_mask = zeros(size(day1_mask));
        temp_mask(15,15) = 1;
        bwout_full = zeros(size(fov_avg{1}));
        bwout_full(...
        a:b,...
        c:d) = temp_mask*icell;
        cellImageAlign(icell).day2_mask = bwout_full;
        mask_all = mask_all+bwout_full;
        mask_temp = mask_all;
        mask_temp(find(mask_all>=1)) = 1;
        mask_exp = imCellBuffer(mask_temp,buf)+mask_temp;
    end
end
mask_cell = mask_all;
mask_np = imCellNeuropil(mask_cell, buf, np);
figure;
subplot(2,2,1)
imagesc(masks{1})
title('Day 1 masks')
subplot(2,2,2)
imagesc(mask_cell)
title('Day 2 masks after transform')
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['mask_comparison.pdf']),'-dpdf','-bestfit')
figure;
filler = zeros(size(masks{1}));
imshow(cat(3,masks{1},mask_cell,filler)) %what is filler and where are the red/green from ?***
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['mask_overlay.pdf']),'-dpdf','-bestfit')
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'mask_cell', 'mask_np')



%%
%old TCs
match = find([cellImageAlign.pass]); %finds matched cell indices
match_ind = match;
%[match_ind,red_ind] = intersect(match,redCells);
cellTCs_match{1} = cellTCs_all{1}(:,match_ind);

%new TCs
data_tc = stackGetTimeCourses(data3, mask_cell); %tc of d2/3 data
[nFrames nCells] = size(data_tc);
down = 5;
data_reg_down  = stackGroupProject(data3,down);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell); %downsampled data

np_tc = zeros(nFrames,nCells);
np_tc_down = zeros(floor(nFrames./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data3,mask_np(:,:,i)); 
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
npSub_tc = data_tc(:,:)-bsxfun(@times,tcRemoveDC(np_tc(:,:)),np_w); %np-subtracted data based on max skew
 
cellTCs_match{2} = npSub_tc(:,:);

save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']),'cellTCs_match', 'cellTCs_all')
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_multiday_alignment.mat']),'cellImageAlign','fitGeoTAf', 'input_points','base_points', 'fov_avg', 'fov_norm','dfmax','corrmap');
clear data_reg_down
%removed 'match_ind' from after 'cellTCs_all' in line 529
%removed 'redChImg' from before 'cellImageAlign' in line 530



%%
load(fName);

data_trial = permute(reshape(npSub_tc,[nOn+nOff ntrials nCells]),[1 3 2]);
ntrials = size(input.tGratingDirectionDeg,2);
Dir = cell2mat_padded(input.tGratingDirectionDeg);
Dir = Dir(1:ntrials);
    Dirs = unique(Dir);
    nDirs = length(Dirs);
    nOn = input.nScansOn;
    nOff = input.nScansOff;
    nCells = size(npSub_tc,2);
    data_f = mean(data_trial(nOff./2:nOff,:,:),1);
    data_dfof = (data_trial-data_f)./data_f;

 
    ndir = length(Dirs);
    [n, n2] = subplotn(nCells);
    h_dir = zeros(nCells, ndir);
    p_dir = zeros(nCells, ndir);
    base_win = nOff-15:nOff;
    resp_win = nOff+1:nOn;
    base = squeeze(mean(data_dfof(base_win,:,:),1)); %averaging across baseline window
    resp = squeeze(mean(data_dfof(resp_win,:,:),1)); %averaging across response window
    dir_resp = zeros(nCells,ndir);
    [x y] = ttest(resp', base', 'tail','right'); %ttest comparing base to resp for significance
    no_match = find(isnan(x)); %?***
    max_dir = zeros(nCells,1);
    figure;
    for i = 1:nCells
        if ~sum(no_match==i) %not getting this?***
        subplot(n, n2, i)
            for idir = 1:ndir %for each dir
                if nOn>29
                    ind = find(Dir == Dirs(idir)); %find trials of that dir
                else
                    ind = find(Dir(1:ntrials-1) == Dirs(idir));
                end
                [h_dir(i,idir), p_dir(i,idir)] = ttest(resp(i,ind), base(i,ind),'tail','right','alpha', 0.05/(ndir-1)); %ttest of each cell at each dir for base vs. resp
                if h_dir(i,idir) %if this cell/dir is sig make red, if not make black and plot
                    errorbar(Dirs(idir), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'or')
                else
                    errorbar(Dirs(idir), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'ok')
                end
                dir_resp(i,idir) = mean(resp(i,ind)-base(i,ind),2); %avg response at each dir for each cell
                hold on
            end
            if sum(h_dir(i,:),2)>0 %if cell has one sig dir
                temp_resp = dir_resp(i,:);
                temp_resp(find(h_dir(i,:)==0)) = NaN;
                [max_val max_ind] = max(temp_resp,[],2); %find max index and value
                max_dir(i,:) = max_ind; %which index/dir was highest?
            else
                [max_val max_ind] = max(dir_resp(i,:),[],2);
                max_dir(i,:) = max_ind;
            end
            title([num2str(Dirs(max_dir(i,:))) ' deg'])
        end
    end


%%
nOri = nDir;
Oris = dirs;
nori = nOri;
Ori = dir_mat;
Oris = unique(Ori);
    h_ori = zeros(nCells, nori);
    p_ori = zeros(nCells, nori);
    ori_resp = zeros(nCells,nori);

    max_ori = zeros(nCells,1);
    figure;
    for i = 1:nCells
        if ~sum(no_match==i)
            subplot(n, n2, i)
            for iori = 1:nori
                if nOn>29
                    ind = find(Ori == Oris(iori));
                else
                    ind = find(Ori(1:ntrials-1) == Oris(iori));
                end
                [h_ori(i,iori), p_ori(i,iori)] = ttest(resp(i,ind), base(i,ind),'tail','right','alpha', 0.05/(nori-1));
                if h_ori(i,iori)
                    errorbar(Oris(iori), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'or')
                else %should this be avg resp or avg resp-base?***
                    errorbar(Oris(iori), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'ok')
                end
                ori_resp(i,iori) = mean(resp(i,ind)-base(i,ind),2);
                ori_stderror(i,iori) = std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind));
                hold on
               
            end
            if sum(h_ori(i,:),2)>0
                temp_resp = ori_resp(i,:);
                temp_resp(find(h_ori(i,:)==0)) = NaN;
                [max_val max_ind] = max(temp_resp,[],2);
                max_ori(i,:) = max_ind;
            else
                [max_val, max_ind] = max(ori_resp(i,:),[],2);
                max_ori(i,:) = max_ind;
            end
            title([num2str(Oris(max_ori(i,:))) ' deg'])
        end
    end

avg_resp_ori = ori_resp.';
std_err_ori = ori_stderror.';

print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuning.pdf']),'-dpdf')
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_trialData.mat']),'data_dfof', 'h_ori', 'max_ori','avg_resp_ori',"std_err_ori")

%%
%replace stim ids with ori values

[Aval, ~, indAval] = unique(max_ori);

 Avalnew = [unique(dir_mat)]; 
 max_ori = Avalnew(indAval);
pref_ori = max_ori;

save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_prefori.mat']),'pref_ori')
