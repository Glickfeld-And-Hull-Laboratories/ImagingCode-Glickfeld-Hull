clear all
clear global
close all
ds = 'DART_V1_PV_contrast'; %dataset info
dataStructLabels = 'contrastxori';
rc = behavConstsAV; %directories
eval(ds)
doAlignRedImage = true;
doRegRegImage = true;
%%
day1_id = 8;
day2_id = 9;
day3_id = nan;
%%
mouse = expt(day1_id).mouse;
fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging');
fnout = fullfile(fn,'multiday');
mkdir(fnout)
%% load up registered data
% day 1
expDate_day1 = expt(day1_id).date;
runs = eval(['expt(day1_id).' dataStructLabels '_runs']);
runFolder = runs{1};
fName = [runFolder '_000_000'];
data_temp = loadsbx_choosepmt(1,mouse,expDate_day1,runFolder,fName);
data_sub = data_temp - min(data_temp(:));
load(fullfile(fn,expDate_day1,runFolder,'regOuts&Img'),'outs')
[~,data_day1] = stackRegister_MA(data_sub,[],[],double(outs));

% day 2
expDate_day2 = expt(day2_id).date;
runs = eval(['expt(day2_id).' dataStructLabels '_runs']);
runFolder = runs{1};
fName = [runFolder '_000_000'];
data_temp = loadsbx_choosepmt(1,mouse,expDate_day2,runFolder,fName);
data_sub = data_temp - min(data_temp(:));
load(fullfile(fn,expDate_day2,runFolder,'regOuts&Img'),'outs')
[~,data_day2] = stackRegister_MA(data_sub,[],[],double(outs));
%% avg FOV each day
brightnessScaleFactor = 0.2;
day1_avg = double(mean(data_day1,3));
day1_norm = uint8((day1_avg./max(day1_avg(:))).*255);
day1_norm(day1_norm > (brightnessScaleFactor*255)) = brightnessScaleFactor*255;
figure;imagesc(day1_norm);colormap gray;title('Day 1 Data Avg')
day2_avg = double(mean(data_day2,3));
day2_norm = uint8((day2_avg./max(day2_avg(:))).*255);
day2_norm(day2_norm > (brightnessScaleFactor*255)) = brightnessScaleFactor*255;
figure;imagesc(day2_norm);colormap gray;title('Day 2 Data Avg')

%% choose alignment method
brightnessScaleFactor = 0.2;
if doAlignRedImage
    load(fullfile(fn,expDate_day1,'data processing','redImage.mat'))
    day1_red = uint8((redChImg./max(redChImg(:))).*255);
    load(fullfile(fn,expDate_day2,'data processing','redImage.mat'))
    day2_red = uint8((redChImg./max(redChImg(:))).*255);
    figure; colormap gray
    subplot 121; imagesc(day1_red); title('Day 1 Tagged')
    subplot 122; imagesc(day2_red); title('Day 2 Tagged')
    day1_align = day1_red;
    day2_align = day2_red;
else
    day1_align = day1_norm;
    day2_align = day2_norm;
end

%% corr map each day
data_day1_down = stackGroupProject(data_day1,floor(size(data_day1,3)./100));
day1_corrmap = getPixelCorrelationImage(double(data_day1_down));
figure;imagesc(day1_corrmap);colormap gray;title('Day 1 Pixel Correlation')
data_day2_down = stackGroupProject(data_day2,floor(size(data_day2,3)./100));
day2_corrmap = getPixelCorrelationImage(double(data_day2_down));
figure;imagesc(day2_corrmap);colormap gray;title('Day 2 Pixel Correlation')

%% manually align
[input_points, base_points] = cpselect(day2_align,day1_align,'Wait', true);
fitGeoTAf = fitgeotrans(input_points(:,:), base_points(:,:),'affine'); 
    
data_day2_trans = imwarp(data_day2,fitGeoTAf, 'OutputView', imref2d(size(data_day2))); 
day2_trans_avg = mean(data_day2_trans,3);
day2_trans_norm = uint8((day2_trans_avg./max(day2_trans_avg(:))).*255);
% day2_trans_corrmap = getPixelCorrelationImage(double(data_day2_trans));
day2_corrmap_norm = uint8((day2_corrmap./max(day2_corrmap(:))).*255);
day2_trans_corrmap = double(imwarp(day2_corrmap_norm,...
    fitGeoTAf, 'OutputView', imref2d(size(day2_corrmap_norm))));

figure;colormap gray
subplot 311
imshow(day1_norm); title('Day 1 Data Avg')
subplot 312
imshow(day2_trans_norm); title('Transformed Day 2 Data Avg')
subplot 313
filler = zeros(size(day1_avg));
imshow(cat(3,day1_norm,day2_trans_norm,filler))
title('Overlay')
print(fullfile(fnout,'FOV manual alignment'),'-dpdf','-fillpage')

figure;colormap gray
d1 = uint8((day1_corrmap./max(day1_corrmap(:))).*255);
d2 = uint8((day2_trans_corrmap./max(day2_trans_corrmap(:))).*255);
subplot 311
imshow(d1); title('Day 1 Pixel Correlation Map')
subplot 312
imshow(d2); title('Transformed Day 2 Pixel Correlation Map')
subplot 313
imshow(cat(3,d1,d2,filler))
title('Overlay')
print(fullfile(fnout,'FOV correlation map manual alignment'),'-dpdf','-fillpage')
%% cell-by-cell correlation
% size of cell box
w=30;
h=30;
buf = 4;
np = 6;
% green channel
% get cell centroids
load(fullfile(fn,expDate_day1,'data processing','final_mask_cells.mat'))
cellPosition = regionprops(green_mask_cell);
nc = length(cellPosition);

xCenter = cellfun(@(a) round(a(1)),{cellPosition.Centroid});
yCenter = cellfun(@(a) round(a(2)),{cellPosition.Centroid});

% index cells NOT too close to edge and NOT in black part of transformation
[ypix,xpix] = size(day1_avg);
goodCells = xCenter>(w/2) & xCenter<xpix-(w/2) & yCenter>(h/2) & yCenter<ypix-(h/2);

goodCells = goodCells & ...
    arrayfun(@(x) sum(sum(day2_trans_norm(green_mask_cell==x)))>0,1:nc);

% fine register each cell
exampleCells = [2,22,33];
green_cellImageAlign = struct;
green_cellTCs = struct;
for icell = 1:nc
    if goodCells(icell)
        % find best shift
        day1_cell_avg = day1_avg(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        day2_cell_avg = day2_trans_avg(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        [reg_avg, shift_avg] = shift_opt(day2_cell_avg,day1_cell_avg,4);
        r_avg = corr(reg_avg(:),day1_cell_avg(:));
        
        day1_cell_corr = day1_corrmap(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        day2_cell_corr = day2_trans_corrmap(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        [reg_corr, shift_corr] = shift_opt(day2_cell_corr,day1_cell_corr,4);
        r_corr = corr(reg_corr(:),day1_cell_corr(:));
        
        if r_avg > 0.8 && r_corr > 0.4
            pass = true;
            if r_avg > r_corr
                shifts = shift_avg;
            else
                shifts = shift_corr;
            end
        elseif r_corr > 0.6 && r_avg >0.4
            pass = true;
            shifts = shift_corr;
        else
            pass = false;
            shifts = nan;
        end
        
        
        mask = green_mask_cell(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        mask(mask > icell | mask < icell) = 0;
        mask(mask == icell) = 1;
        
        green_cellImageAlign(icell).center_yx = [yCenter(icell),xCenter(icell)];
        green_cellImageAlign(icell).mask = mask;
        green_cellImageAlign(icell).day1.avg_img = day1_cell_avg;
        green_cellImageAlign(icell).day1.corr_img = day1_cell_corr;
        green_cellImageAlign(icell).day2.avg_img = day2_cell_avg;
        green_cellImageAlign(icell).day2.corr_img = day2_cell_corr;
        green_cellImageAlign(icell).r_avg = r_avg;
        green_cellImageAlign(icell).r_corr = r_corr;
        green_cellImageAlign(icell).shifts = shifts;
        
        % shift data, get tc, all days
        if pass
            [~,data_day2_cell_shift] = stackRegister_MA(data_day2_trans(...
                yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
                xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1,:),[],[],...
                repmat(shifts,[length(data_day2_trans),1]));
            d1 = data_day1(yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
                xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1,:);
            tc_d1 = stackGetTimeCourses(d1,mask);
            tc_d2 = stackGetTimeCourses(data_day2_cell_shift,mask);
            tc_day1 = getWeightedNeuropilTimeCourse(d1,tc_d1,mask,buf,np);
            tc_day2 = getWeightedNeuropilTimeCourse(data_day2_cell_shift,tc_d2,mask,buf,np);

            
            green_cellTCs(icell).pass = true;
            green_cellTCs(icell).tc_day1 = tc_day1;
            green_cellTCs(icell).tc_day2 = tc_day2;
            if any(exampleCells==icell)
                figure
                suptitle(sprintf('Cell %s, pix map corr = %s, avg corr = %s',...
                    num2str(icell),sigfigString(r_corr),sigfigString(r_avg)))
                subplot 121
                imagesc(mean(d1,3))
                title('Day 1, Avg')
                subplot 122
                imagesc(mean(data_day2_cell_shift,3))
                title('Day 2, Avg')
            end
        else
            green_cellTCs(icell).pass = false;
            green_cellTCs(icell).tc_day1 = [];
            green_cellTCs(icell).tc_day2 = [];            
        end
    else
        green_cellTCs(icell).pass = false;
        green_cellTCs(icell).tc_day1 = [];
        green_cellTCs(icell).tc_day2 = []; 
    end
end

% red channel
cellPosition = regionprops(red_mask_cell);
nc = length(cellPosition);

xCenter = cellfun(@(a) round(a(1)),{cellPosition.Centroid});
yCenter = cellfun(@(a) round(a(2)),{cellPosition.Centroid});

% index cells NOT too close to edge and NOT in black part of transformation
[ypix,xpix] = size(day1_avg);
goodCells = xCenter>(w/2) & xCenter<xpix-(w/2) & yCenter>(h/2) & yCenter<ypix-(h/2);

goodCells = goodCells & ...
    arrayfun(@(x) sum(sum(day2_trans_norm(red_mask_cell==x)))>0,1:nc);

% fine register each cell
exampleCells = [2,4,6];
red_cellImageAlign = struct;
red_cellTCs = struct;
for icell = 1:nc
    if goodCells(icell)
        % find best shift
        day1_cell_avg = day1_avg(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        day2_cell_avg = day2_trans_avg(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        [reg_avg, shift_avg] = shift_opt(day2_cell_avg,day1_cell_avg,4);
        r_avg = corr(reg_avg(:),day1_cell_avg(:));
        
        day1_cell_corr = day1_corrmap(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        day2_cell_corr = day2_trans_corrmap(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        [reg_corr, shift_corr] = shift_opt(day2_cell_corr,day1_cell_corr,4);
        r_corr = corr(reg_corr(:),day1_cell_corr(:));
        
        if r_avg > 0.8 && r_corr > 0.4
            pass = true;
            if r_avg > r_corr
                shifts = shift_avg;
            else
                shifts = shift_corr;
            end
        elseif r_corr > 0.6 && r_avg >0.4
            pass = true;
            shifts = shift_corr;
        else
            pass = false;
            shifts = nan;
        end
        
        
        mask = red_mask_cell(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        mask(mask > icell | mask < icell) = 0;
        mask(mask == icell) = 1;
        
        red_cellImageAlign(icell).center_yx = [yCenter(icell),xCenter(icell)];
        red_cellImageAlign(icell).mask = mask;
        red_cellImageAlign(icell).day1.avg_img = day1_cell_avg;
        red_cellImageAlign(icell).day1.corr_img = day1_cell_corr;
        red_cellImageAlign(icell).day2.avg_img = day2_cell_avg;
        red_cellImageAlign(icell).day2.corr_img = day2_cell_corr;
        red_cellImageAlign(icell).r_avg = r_avg;
        red_cellImageAlign(icell).r_corr = r_corr;
        red_cellImageAlign(icell).shifts = shifts;
        
        % shift data, get tc, all days
        if pass
            [~,data_day2_cell_shift] = stackRegister_MA(data_day2_trans(...
                yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
                xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1,:),[],[],...
                repmat(shifts,[length(data_day2_trans),1]));
            d1 = data_day1(yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
                xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1,:);
            tc_d1 = stackGetTimeCourses(d1,mask);
            tc_d2 = stackGetTimeCourses(data_day2_cell_shift,mask);
            tc_day1 = getWeightedNeuropilTimeCourse(d1,tc_d1,mask,buf,np);
            tc_day2 = getWeightedNeuropilTimeCourse(data_day2_cell_shift,tc_d2,mask,buf,np);

            
            red_cellTCs(icell).pass = true;
            red_cellTCs(icell).tc_day1 = tc_day1;
            red_cellTCs(icell).tc_day2 = tc_day2;
            if any(exampleCells==icell)
                figure
                suptitle(sprintf('Cell %s, pix map corr = %s, avg corr = %s',...
                    num2str(icell),sigfigString(r_corr),sigfigString(r_avg)))
                subplot 121
                imagesc(mean(d1,3))
                title('Day 1, Avg')
                subplot 122
                imagesc(mean(data_day2_cell_shift,3))
                title('Day 2, Avg')
            end
        else
            red_cellTCs(icell).pass = false;
            red_cellTCs(icell).tc_day1 = [];
            red_cellTCs(icell).tc_day2 = [];            
        end
    else
        red_cellTCs(icell).pass = false;
        red_cellTCs(icell).tc_day1 = [];
        red_cellTCs(icell).tc_day2 = []; 
    end
end


save(fullfile(fnout,'timecourses'),'green_cellTCs','red_cellTCs')
save(fullfile(fnout,'multiday_alignment'),'green_cellImageAlign',...
    'red_cellImageAlign')