clear all; clear global; close all
clc
ds = 'ExperimentData_MultiDayMatch_TD'; %dataset info
dataStructLabels = {'stimruns'};
eval(ds)


day_id(2) = 3;
day_id(1) = expt(day_id(2)).multiday_matchdays;

nd = length(day_id);
brightnessScaleFactor = 0.3;
mouse = expt(day_id(1)).mouse;

%Path names
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
td_fn = fullfile(fn_base, 'home\Tierney');
data_fn = fullfile(td_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data');
fnIn = fullfile(fn_base, 'home\Tierney\Analysis\2P');

if expt(day_id(2)).multiday_time_days>0
    time_str = [num2str(expt(day_id(2)).multiday_time_days) 'Days'];
else
    time_str = 'control';
end

fn_multi = fullfile(fn_base, 'home\Tierney\Analysis\2P',mouse,['multiday_' time_str,'_',expt(day_id(2)).experiment]);
mkdir(fn_multi)
data = cell(1,nd);
fov_avg = cell(1,nd);
fov_norm = cell(1,nd);
align_fov = cell(1,nd);
corrmap = cell(1,nd);
dfmax = cell(1,nd);
masks = cell(1,nd);
maskNP = cell(1,nd);
cellTCs_all = cell(1,nd);
%% load all data 
for id = 1:nd 
    clear global
    ImgFolder = expt(day_id(id)).stimruns
    expDate = expt(day_id(id)).date;
    ImgFolder = eval(['expt(day_id(' num2str(id) ')).' cell2mat(dataStructLabels)]);
    nrun = length(ImgFolder);
    time = expt(day_id(id)).time
    contra = strcmp(expt(day_id(id)).eye_str,'Contra');

    %Load 2P data
data_day = [];
clear temp
offset_frames = 0;
offset_time = 0;

    for irun = 1:nrun
        %Load mworks data- this has information about experiment (e.g. the visual stimuli presented and synchronization with the microscope)
        fName = fullfile(mworks_fn, ['data-' mouse '-' expDate '-' time{irun} '.mat']);
        load(fName);
        ntrials = size(input.tGratingDirectionDeg,2);
        totframes = input.counterValues{end}(end); %this is from the mworks structure- finds the last value clocked for frame count
        input.contraTrialNumber = mat2cell(contra(irun).*ones(1,ntrials),1,ones(1,ntrials)); %create a field for tracking which eye stimulated
        temp(irun) = input;

        if irun>1
            for itrial = 1:ntrials
                temp(irun).counterValues{itrial} = bsxfun(@plus,temp(irun).counterValues{itrial},offset_frames);
                temp(irun).counterTimesUs{itrial} = bsxfun(@plus,temp(irun).counterTimesUs{itrial},offset_time-temp(irun).counterValues{1}(1)+1000);
                temp(irun).wheelSpeedTimesUs{itrial} = bsxfun(@plus,temp(irun).wheelSpeedTimesUs{itrial},offset_time-temp(irun).counterValues{1}(1)+1000);
            end
        end
        offset_frames = offset_frames+totframes;
        offset_time = offset_time+temp(irun).counterTimesUs{end}(end);

        %Load 2P metadata- this has information about the information about the imaging session (e.g. frame count, zoom)
        CD = fullfile(data_fn, mouse, expDate, ImgFolder(irun,:));
        cd(CD);
        imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
        load(imgMatFile);
       
        %Load 2P images

        fprintf(['Reading ' num2str(totframes) ' frames \r\n'])
        data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,totframes); %loads the .sbx files with imaging data (path, nframes to skip, nframes to load)
        data_temp = squeeze(data_temp);
        data_day = cat(3,data_day,data_temp);
    end

input = concatenateDataBlocks(temp);
clear data_temp
clear temp
    
% Registration
run_str = catRunName(ImgFolder, nrun);
datemouse = [expDate '_' mouse];
datemouserun = [expDate '_' mouse '_' run_str];

    if exist(fullfile(fnIn, datemouse, datemouserun)) 
        load(fullfile(fnIn, datemouse, datemouserun, [datemouserun '_reg_shifts.mat']))
        fprintf(['Starting registration for day ' num2str(id)])
        [~,data{id}] = stackRegister_MA(data_day,[],[],double(out));
        save(fullfile(fnIn, datemouse, datemouserun, [datemouserun '_input.mat']), 'input')
    else
        fprintf(['No registration found'])
    end 
    
    data_avg = mean(data{id},3);
    figure; imagesc(data_avg); title(['Avg FOV day ' num2str(id)])
    fov_avg{id} = data_avg;
    fov_norm{id} = uint8((fov_avg{id}./max(fov_avg{id}(:))).*255);
    fov_norm{id}(fov_norm{id} > (brightnessScaleFactor*255)) = brightnessScaleFactor*255;

    load(fullfile(fnIn,datemouse, datemouserun, [datemouserun '_mask_cell.mat']))
    dfmax{id} = data_dfof(:,:,1);
    corrmap{id} = data_dfof(:,:,end);
    masks{id} = mask_cell;
    maskNP{id} = mask_np;
    load(fullfile(fnIn,datemouse,datemouserun,[datemouserun '_TCs.mat']))
    cellTCs_all{id} = npSub_tc;
    load(fullfile(fnIn,datemouse,datemouserun,[datemouserun '_input.mat']))
    input_temp(id) = input;
end
input = input_temp;
save(fullfile(fn_multi,'input.mat'),'input')
clear input
%% manual align
corrmap_norm{1} = uint8((corrmap{1}./max(corrmap{1}(:))).*255);
corrmap_norm{2} = uint8((corrmap{2}./max(corrmap{2}(:))).*255);
if exist(fullfile(fn_multi,'multiday_alignment.mat'))
    load(fullfile(fn_multi,'multiday_alignment.mat'))
else
    [input_points_1, base_points_1] = cpselect(fov_norm{2},fov_norm{1},'Wait', true);
    [input_points_2, base_points_2] = cpselect(corrmap_norm{2},corrmap_norm{1},'Wait', true);
    input_points = [input_points_1; input_points_2];
    base_points = [base_points_1; base_points_2];
end
%%
fitGeoTAf = fitgeotrans(input_points(:,:), base_points(:,:),'affine'); 
data{3} = imwarp(data{2},fitGeoTAf, 'OutputView', imref2d(size(data{2}))); 
fov_avg{3} = mean(data{3},3);
fov_norm{3} = uint8((fov_avg{3}./max(fov_avg{3}(:))).*255);
corrmap{3} = double(imwarp(corrmap{2},fitGeoTAf, 'OutputView', imref2d(size(corrmap{2}))));
corrmap_norm{3} = uint8((corrmap{3}./max(corrmap{3}(:))).*255);
red_trans = (imwarp(fov_red{2},fitGeoTAf, 'OutputView', imref2d(size(fov_red{2})))); %%%%%%%%%%%% DELETE THIS AND NEXT LINE?
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

%%%%%%% DELETE THIS SECTION???
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
            day1_red_avg = fov_red{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
            day2_red_avg = fov_red{3}(...
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
        if max_val>0.55 & (r_corr>0.4 || r_red>0.4 || r_max>0.4)
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
                        img_select = fov_red{3};
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