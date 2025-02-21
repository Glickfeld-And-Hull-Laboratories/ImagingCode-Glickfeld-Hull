%% This is essentially the same as the normal Day23 Image Processing, but this involves a special case where data that were previously 
%matched to day1 are used as 'day1' for the current comparison being made. 

%this script is used when only 1 channel imaging was performed (green)*	
clear all;
clear global;
close all;
clc;
%% get path names D2
ref_date = '230124';
date = '230127';
time = strvcat('0805');
alignToRef = 1;
ImgFolder = strvcat('001');
mouse = 'i2543';
nrun = size(ImgFolder,1);
frame_rate = 15.5;
ref_str = 'runs-001';
run_str = catRunName(ImgFolder, nrun);
tj_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\2P_Imaging';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P';
behav_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
%% load and register - same as d1 code
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = fullfile(tj_fn, [mouse '\' date '_' mouse '\' ImgFolder(irun,:)]);
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    fName = fullfile(behav_fn, ['data-'  mouse  '-' date '-' time(irun,:) '.mat']); %find behavior data
    load(fName);

    nframes = info.config.frames;
    fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n'])
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    
    
    temp(irun) = input;
    if isfield(input, 'nScansOn')
        nOn = temp(irun).nScansOn;
        nOff = temp(irun).nScansOff;
        ntrials = size(temp(irun).tGratingDirectionDeg,2);

        data_temp = squeeze(data_temp);
        if nframes>ntrials*(nOn+nOff)
            nframes = ntrials*(nOn+nOff);
            data_temp = data_temp(:,:,1:ntrials*(nOn+nOff));
        elseif nframes<ntrials*(nOn+nOff)
            temp(irun) = trialChopper(temp(irun),1:ceil(nframes./(nOn+nOff)));
        end
    end
    
    offset = offset+nframes;

    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    trial_n = [trial_n nframes];
end
input = concatenateDataBlocks(temp);
clear data_temp
clear temp

 %% Choose register interval - averaging 500 frames, also skipping 2000
t = 2000;
nep = floor(size(data,3)./t);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*t):500+((i-1)*t)),3)); title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]); end

%% Register data - identify clearest stack and align frames to that

data_avg = mean(data(:,:,34001:34500),3); 

if exist(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
    load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(data(:,:,:),[],[],out);
    clear out outs
else
    [out, data_reg] = stackRegister(data,data_avg);
    data_reg_avg = mean(data_reg,3);
    reg = data_reg_avg;
    mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'data_reg_avg', 'out', 'data_avg')
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
clear data

%data_avg = selected stack to register to
%data_reg = all frames registered
%data_reg_avg = mean of all registered frames

%% test stability - see how the image looks averaged across all frames now
figure; imagesq(data_reg_avg); truesize;
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf','-bestfit')

%% find activated cells - again this is same as d1 code
%max by trial type    
if isfield(input, 'nScansOn')
    nOn = input.nScansOn;
    nOff = input.nScansOff;
    if nOn>29
        sz = size(data_reg);
%         make data_tr smaller for dfof images to not use all of the memory
        data_tr = reshape(single(data_reg(:,:,1:28800)),[sz(1), sz(2), nOn+nOff, 320]);
        data_f = mean(double(data_tr(:,:,nOff/2:nOff,:)),3);
        data_df = bsxfun(@minus, double(data_tr), data_f); 
        data_dfof = bsxfun(@rdivide,data_df, data_f); 
        clear data_tr clear data_f clear data_df
    else
        sz = size(data_reg);
        data_tr = zeros(sz(1),sz(2), 100, ntrials-1);
        for itrial = 1:ntrials-1
            data_tr(:,:,:,itrial) = data_reg(:,:,((itrial-1)*(nOn+nOff))+71:170+((itrial-1)*(nOn+nOff)));
        end
        data_f = mean(data_tr(:,:,1:50,:),3);
        data_df = bsxfun(@minus, double(data_tr), data_f); 
        data_dfof = bsxfun(@rdivide,data_df, data_f); 
        clear data_tr clear data_f clear data_df
    end
end

if input.doDirStim
%     obtaining avg and max dfof images for each dir
    Dir = cell2mat_padded(input.tGratingDirectionDeg);
    Dir = Dir(1:ntrials);
    Dirs = unique(Dir);
    data_dfof_avg = zeros(sz(1),sz(2),length(Dirs));
    nDirs = length(Dirs);
    [n n2] = subplotn(nDirs);
    figure;
    for idir = 1:length(Dirs)
        if nOn>29
            ind = find(Dir(1:160) == Dirs(idir));
        else
            ind = find(Dir(1:ntrials-1) == Dirs(idir));
        end
        data_dfof_avg(:,:,idir) = mean(mean(data_dfof(:,:,nOff+1:nOn+nOff,ind),3),4);
        subplot(n,n2,idir)
        imagesc(data_dfof_avg(:,:,idir))
    end
    clear data_dfof
    myfilter = fspecial('gaussian',[20 20], 0.5);
    data_dfof_avg_all = imfilter(data_dfof_avg,myfilter);
    data_dfof_max = max(data_dfof_avg_all,[],3);
    
%     average dfof image for each direction - filtered
    figure; 
    Stims = Dirs;
    nStim = length(Dirs);
    [n n2] = subplotn(nDirs);
    data_dfof_avg_ori = zeros(sz(1), sz(2), nDirs/2);
    for i = 1:nStim 
        subplot(n,n2,i); 
        imagesc(data_dfof_avg_all(:,:,i));
        clim([0 max(data_dfof_avg_all(:))])
        title(num2str(Dirs(i)))
        colormap(gray)
        if i<(nDirs/2)+1
            data_dfof_avg_ori(:,:,i) = mean(data_dfof_avg_all(:,:,[i i+nDirs/2]),3); 
        end
    end

    print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_16Stim.pdf']), '-dpdf')


%     average dfof image for each orientation
    figure;
    [n n2] = subplotn(nDirs/2);
    for i = 1:nStim/2
        subplot(n,n2,i)
        imagesc(data_dfof_avg_ori(:,:,i));
        clim([0 max(data_dfof_avg_ori(:))])
        title(num2str(Dirs(i)))
        axis off
    end
    subplot(n,n2,i+1)
    imagesc(max(data_dfof_avg_ori,[],3))
    title('dfof Max')
    axis off
    data_dfof = cat(3,data_dfof_avg_ori,max(data_dfof_avg_ori,[],3));

    print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_activeCells.pdf']), '-dpdf')

    figure;
    imagesc(max(data_dfof_avg_ori,[],3))
    title('dfof Max')
    axis off

    print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofMax.pdf']), '-dpdf')

end

 save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg_all')

 %data_dfof = avg image for each ori + max
 %data_dfof_avg = avg image for each dir
 %data_dfof_avg_all = filtered version of dfof_avg
 %data_dfof_avg_ori = avg image for each ori
 %data_dfof_max = max of each pixel
%% create pixel correlation image - same as d1
data_reg_3hz = stackGroupProject(data_reg,5);
pix = getPixelCorrelationImage(data_reg_3hz);
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


%the main difference between this and the original day2/3 processing script is that this one uses
%the rematch_shifts as input
%load cell mask, time course, dfof, shifts, pix corr, and input from D1
maskD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_mask_cell.mat']));
TCs_D1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_TCs.mat']));
dfofD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_stimActFOV.mat']));
shiftsD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_reg_shifts.mat']));
pixelD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_pixel.mat']));
inputD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_input.mat']));
multiDayD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_multiday_alignment.mat']));
rematch_d1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_rematch_shifts.mat']));

%Need alternative load info if matching previously done!!
fov_avg{1} = rematch_d1.fov_avg_shift; %this uses shifted* data
dfmax{1} = rematch_d1.dfmax_shift; %d1 max
cellTCs_all{1} = TCs_D1.cellTCs_match{2}; %np-sub TC from D1 when it is registered as D2, which it is for this script*
input_temp(1) = inputD1; %input (behavioral setup) from d1
corrmap{1} = pixelD1.pix; %pix corr from d1
masks{1} = maskD1.mask_cell; %cell mask from d1
corrmap_norm{1} = uint8((corrmap{1}./max(corrmap{1}(:))).*255); %normalizing again
brightnessScaleFactor = 0.5;
fov_norm{1} = uint8((fov_avg{1}./max(fov_avg{1}(:))).*255);
fov_norm{1}(fov_norm{1} > (brightnessScaleFactor*255)) = brightnessScaleFactor*255;


%% align Day 2 data stack to Day 1
[input_points_1, base_points_1] = cpselect(fov_norm{2},fov_norm{1},'Wait', true); %select control points from 2 images (fov images); 2 is matched to 1
[input_points_2, base_points_2] = cpselect(dfmax{2},dfmax{1},'Wait', true); %select control points from max projection
input_points = [input_points_1' input_points_2']; %x,y pairings of points on d2
base_points = [base_points_1' base_points_2']; %x,y pairings of points on d1
fitGeoTAf = fitgeotrans(input_points(:,:)', base_points(:,:)','affine'); %fits your input points to the base points

% transform data_reg and make a new D2 FOV from the avg of that
data3 = imwarp(data_reg,fitGeoTAf, 'OutputView', imref2d(size(data_reg))); %displacing d2 data reg by the fits above 
fov_avg{3} = mean(data3,3); %new fov
fov_norm{3} = uint8((fov_avg{3}./max(fov_avg{3}(:))).*255); %normalized
corrmap{3} = double(imwarp(corrmap{2},fitGeoTAf, 'OutputView', imref2d(size(corrmap{2})))); %new pix corr map
corrmap_norm{3} = uint8((corrmap{3}./max(corrmap{3}(:))).*255);
dfmax{3} = imwarp(dfmax{2},fitGeoTAf, 'OutputView', imref2d(size(dfmax{2}))); %new dfof max
%redChImg = imwarp(fov_red{2},fitGeoTAf, 'OutputView', imref2d(size(fov_red{2})));

%SAVE all of the this!! - you may need these new shifts to now use this session as a 'day1' for
%another comparison

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
figure; image(rgb); movegui('center') % overlays d1 and d2 images

figure;colormap gray 
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


%fov_avg{1} is the d1 data
%fov_avg{2} is the NON-shifted d2/3 data
%fov_avg{3} is the shifted d2/3 data

%% make new cell masks for D2
clear input 
% size of cell box
close all
w=30; %width
h=30; %height
%buf = 3; %buffer
buf = 4; %buffer
np = 5; %neuropil

% green channel
% get cell centroids

cellImageAlign = struct; %create structure

cellPosition = regionprops(masks{1}); %properties of image regions - results in each cell centroid defined along with an 'area'
nc = length(cellPosition);

xCenter = cellfun(@(a) round(a(1)),{cellPosition.Centroid}); %apply function to each cell of cell array - rounds x and y centroids
yCenter = cellfun(@(a) round(a(2)),{cellPosition.Centroid}); 

% index cells NOT too close to edge and NOT in black part of transformation

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
                day1_mask = masks{1}(... %%need original D1 mask trans to D2
                    yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,... %identifying each cell as a cluster of pixels in mask
                    xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
                day1_cell_avg = fov_avg{1}(... %doing same as above but with fov on d1
                    yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
                    xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
                day2_cell_avg = fov_avg{3}(... %same as above but with d3
                    yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
                    xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
                [reg_avg, shift_avg] = shift_opt(day2_cell_avg,day1_cell_avg,2); % finds optimal shift to align data and target using correlation 
                r_avg = corr(reg_avg(:),day1_cell_avg(:)); %corr between reg_avg and cell_avg
                
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
                [max_val max_ind] = max([r_avg r_max r_corr]); %find max corr values
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
                    x = input(prompt); 
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
                cellImageAlign(icell).shifts = shifts; %shifts
                cellImageAlign(icell).pass = pass; %did not press 0
        
        
                % shift data, get tc, all days
                if pass %this is where you are selecting cell pixels
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
                        temp_mask = zeros(size(day1_mask));
                        temp_mask(15,15) = 1; %*this is different here because we were having trouble with cells overlapping in this version of matching (aka rematching)
                        bwout_full = zeros(size(fov_avg{1}));
                        bwout_full(...
                        yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
                        xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1) = temp_mask*icell;
                        mask_all = mask_all+bwout_full;
                        mask_temp = mask_all;
                        mask_temp(find(mask_all>=1)) = 1;
                        mask_exp = imCellBuffer(mask_temp,buf)+mask_temp;
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
                mask_temp(mask_all>=1) = 1;
                mask_exp = imCellBuffer(mask_temp,buf)+mask_temp;
                end
            %end
        %end
    else
        %cells are placeholders despite not being in new FOV
        a = yCenter(icell)-(h/2); %here to line 446 for cells on edge
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
imshow(cat(3,masks{1},mask_cell,filler)) 
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
load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']));
ntrials = size(input.tGratingDirectionDeg,2);
npSub_tc = cellTCs_match{2};
Dir = cell2mat_padded(input.tGratingDirectionDeg);
Dir = Dir(1:ntrials);
    Dirs = unique(Dir);
    nDirs = length(Dirs);
if isfield(input, 'nScansOn')
    nOn = input.nScansOn;
    nOff = input.nScansOff;
    nCells = size(npSub_tc,2);
    if nOn>29
        data_mat = zeros(nOn+nOff, nCells, ntrials); %this is trial-level data
        for itrial = 1:ntrials
            data_mat(:,:,itrial) = npSub_tc(1+((itrial-1).*(nOn+nOff)):(itrial.*(nOn+nOff)),:);
        end
        data_f = mean(data_mat(nOff/2:nOff,:,:),1); %2nd half of off frames for baseline
    else
        data_mat = zeros(100, nCells, ntrials-1);
        for itrial = 1:ntrials-1
            data_mat(:,:,itrial) = npSub_tc(((itrial-1)*(nOn+nOff))+71:170+((itrial-1)*(nOn+nOff)),:);
        end
        data_f = mean(data_mat(1:50,:,:),1);
    end
    data_df = bsxfun(@minus, data_mat, data_f);
    data_dfof = bsxfun(@rdivide, data_df, data_f); %df/f data matrix
    
    ndir = length(Dirs);
    [n, n2] = subplotn(nCells);
    h_dir = zeros(nCells, ndir);
    p_dir = zeros(nCells, ndir);
    base_win = 50:60;
    resp_win = 70:90;
    base = squeeze(mean(data_dfof(base_win,:,:),1)); %avg across base
    resp = squeeze(mean(data_dfof(resp_win,:,:),1)); %avg across resp
    resp_mat = squeeze(mean(data_mat(resp_win,:,:),1)); %non-standardized data
%     figure;scatter(mean(resp,2),mean(resp_mat,2));xlabel('dfof resp');ylabel('f resp')
    dir_resp = zeros(nCells,ndir);
    [x y] = ttest(resp', base', 'tail','right'); %what is this for?
    no_match = find(isnan(x));
    max_dir = zeros(nCells,1);
    figure;
    for i = 1:nCells
        if ~sum(no_match==i)
        subplot(n, n2, i)
            for idir = 1:ndir
                if nOn>29
                    ind = find(Dir == Dirs(idir));
                else
                    ind = find(Dir(1:ntrials-1) == Dirs(idir));
                end
                [h_dir(i,idir), p_dir(i,idir)] = ttest(resp(i,ind), base(i,ind),'tail','right','alpha', 0.05/(ndir-1));
                if h_dir(i,idir)
                    errorbar(Dirs(idir), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'or')
                else
                    errorbar(Dirs(idir), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'ok')
                end
                dir_resp(i,idir) = mean(resp(i,ind)-base(i,ind),2);
                hold on
            end
            if sum(h_dir(i,:),2)>0
                temp_resp = dir_resp(i,:);
                temp_resp(find(h_dir(i,:)==0)) = NaN;
                [max_val max_ind] = max(temp_resp,[],2);
                max_dir(i,:) = max_ind;
            else
                [max_val max_ind] = max(dir_resp(i,:),[],2);
                max_dir(i,:) = max_ind;
            end
            title([num2str(Dirs(max_dir(i,:))) ' deg'])
        end
    end
        print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirTuning.pdf']),'-dpdf')

    
    nori = length(Dirs)/2;
    Ori = Dir;
    for iori = 1:nori
        ind = find(Dir == Dirs(iori+nori));
        Ori(ind) = Dirs(iori);
    end
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
                else
                    errorbar(Oris(iori), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'ok')
                end
                ori_resp(i,iori) = mean(resp(i,ind)-base(i,ind),2);
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
    
    good_ind = unique([find(x)'; find(sum(h_dir,2)>0); find(sum(h_ori,2)>0)]); %cells with at least 1 sig dir/ori
    print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuning.pdf']),'-dpdf')
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_trialData.mat']),'data_dfof','resp_mat','max_dir','h_dir', 'h_ori', 'max_ori','good_ind','dir_resp')
end

%% ori fitting
nOn = input.nScansOn;
nOff = input.nScansOff;
dir_mat = celleqel2mat_padded(input.tGratingDirectionDeg);
nTrials = length(dir_mat);
input.trialSinceReset = nTrials;

down = 10;
nframes = size(npSub_tc,1)./down;
nCells = size(npSub_tc,2);
data_tc_down = squeeze(mean(reshape(npSub_tc, [down,nframes,nCells]),1));
tuningDownSampFactor = down;

[avgResponseEaOri,semResponseEaOri,vonMisesFitAllCellsAllBoots,fitReliability,R_square,tuningTC] = ...
    getOriTuningLG(data_tc_down,input,tuningDownSampFactor);
    vonMisesFitAllCells = squeeze(vonMisesFitAllCellsAllBoots(:,1,:));

save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuningAndFits.mat']),...
            'avgResponseEaOri','semResponseEaOri','vonMisesFitAllCellsAllBoots','fitReliability','R_square', 'tuningTC')

%%
dir_mat = celleqel2mat_padded(input.tGratingDirectionDeg);
ori_mat = dir_mat;
ori_mat(find(dir_mat>=180)) = ori_mat(dir_mat>=180)-180;
oris = unique(ori_mat);
figure; 
if nCells<49
    [n n2] = subplotn(nCells);
else
    [n, n2] = subplotn(49);
end
start = 1;
x = 0;
for ic = 1:nCells
    if start > 49 
        suptitle([mouse ' ' date ' n = ' num2str(length(find(fitReliability<22.5))) '/' num2str(nCells) '- well-fit'])
        print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuningFits_cells' num2str(start-49) '-' num2str(start-1) '.pdf']),'-dpdf','-fillpage')
        start = 1;
        x = x+1;
        figure;
    end
    subplot(n,n2,ic-(x.*49))
    errorbar(oris,avgResponseEaOri(ic,:), semResponseEaOri(ic,:),'-o')
    hold on
    plot(0:180,vonMisesFitAllCellsAllBoots(:,1,ic));
    tit_str = num2str(chop(R_square(1,ic),2));
    if fitReliability(ic)<22.5
        tit_str = [tit_str '- R'];
    end
    title(tit_str)
    start = start+1;
end
suptitle([mouse ' ' date ' n = ' num2str(length(find(fitReliability<22.5))) '/' num2str(nCells) '- well-fit'])
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuningFits' num2str(start-49) '-' num2str(start-1) '.pdf']),'-dpdf','-fillpage')

[max_resp prefOri] = max(vonMisesFitAllCellsAllBoots,[],1);
prefOri = squeeze(prefOri)-1;
prefOri_bootdiff = abs(prefOri(2:end,:)-prefOri(1,:)); %how different is each bootstrap from original?***
prefOri_bootdiff(find(prefOri_bootdiff>90)) = 180-prefOri_bootdiff(find(prefOri_bootdiff>90)); 
ind_theta90 = find(prctile(prefOri_bootdiff,90,1)<22.5); 
edges = [0 22.5:45:225]; 
[bin ind_bin] = histc(prefOri(1,:),edges); 
ind_bin(find(ind_bin==5)) = 1;
bin(1) = bin(1)+bin(5);
bin(5) = [];

tunedCells = cell(1,length(bin));
for i = 1:length(bin)
    tunedCells{i} = intersect(find(ind_bin==i),ind_theta90);
end

save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuningInfo.mat']),...
    'prefOri', 'prefOri_bootdiff', 'ind_theta90', 'tunedCells');
