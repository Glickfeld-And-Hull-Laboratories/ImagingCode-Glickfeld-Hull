clear all; clear global; close all
clc
ds = 'DART_V1_contrast_ori'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsAV; %directories
eval(ds)
doGreenOnly = false;
doCorrImg = true;


day_id(2) = 20;
day_id(1) = expt(day_id(2)).multiday_matchdays;

nd = length(day_id);
brightnessScaleFactor = 0.3;
mouse = expt(day_id(1)).mouse;

if strcmp(expt(day_id(1)).data_loc,'lindsey')
    root = fullfile(rc.data,mouse);
elseif strcmp(expt(day_id(1)).data_loc,'ashley')
    root = fullfile(rc.ashleyData,mouse,'two-photon imaging');
end
fnout = fullfile(rc.lindseyAnalysis,mouse);
if expt(day_id(2)).multiday_timesincedrug_hours>0
    dart_str = [expt(day_id(2)).drug '_' num2str(expt(day_id(2)).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end
fn_multi = fullfile(rc.lindseyAnalysis,mouse,['multiday_' dart_str]);
mkdir(fn_multi)
data = cell(1,nd);
fov_avg = cell(1,nd);
fov_norm = cell(1,nd);
fov_red = cell(1,nd);
align_fov = cell(1,nd);
corrmap = cell(1,nd);
dfmax = cell(1,nd);
masks = cell(1,nd);
red_ind = cell(1,nd);
cellTCs_all = cell(1,nd);
%% load all data
runFolder = [];
for id = 1:nd
    clear global
    expDate = expt(day_id(id)).date;
    runs = eval(['expt(day_id(' num2str(id) ')).' cell2mat(dataStructLabels) '_runs']);
    nrun = length(runs);
    out_all = [];
    data_g = [];
    for irun = 1:nrun
        imgFolder = runs{irun};
        fName = [imgFolder '_000_000'];
        cd(fullfile(root,expDate,imgFolder))
        load(fName)
        if nrun == 1
            load(fullfile(fnout,expDate,imgFolder,'regOuts&Img.mat'))
            nframes = size(outs,1);
            runFolder = imgFolder;
        else
            nframes = info.config.nframes;
            runFolder = [runFolder '_' imgFolder];
        end
        fprintf(['Loading day ' num2str(id) ' data \n Mouse: ' expt(day_id(id)).mouse ' Date: ' expt(day_id(id)).date])
        data_temp = sbxread(fName(1,1:11),0,nframes);
        data_g = cat(3, data_g, squeeze(data_temp(1,:,:,:)));
        clear data_temp
        out_all = [out_all; outs];
    end
    [~,data{id}] = stackRegister_MA(data_g,[],[],double(out_all));
    data_avg = mean(data{id},3);
    figure; imagesc(data_avg); title(['Avg FOV day ' num2str(id)])
    clear data_g
    fov_avg{id} = data_avg;
    fov_norm{id} = uint8((fov_avg{id}./max(fov_avg{id}(:))).*255);
    fov_norm{id}(fov_norm{id} > (brightnessScaleFactor*255)) = brightnessScaleFactor*255;

    if exist(fullfile(fnout,expDate,runFolder,'redImage.mat'))
        load(fullfile(fnout,expDate,runFolder,'redImage.mat'))
    	fov_red{id} = uint8((redChImg./max(redChImg(:))).*255);
    else
        fov_red{id} = zeros(size(data_avg));
    end
    load(fullfile(fnout,expDate,runFolder,'mask_cell.mat'))
    dfmax{id} = data_dfof(:,:,1);
    corrmap{id} = data_dfof(:,:,end);
    masks{id} = mask_cell;
    red_ind{id} = mask_label;
    load(fullfile(fnout,expDate,runFolder,'TCs.mat'))
    cellTCs_all{id} = npSub_tc;
    load(fullfile(fnout,expDate,runFolder,'input.mat'))
    input_temp(id) = input;
end
input = input_temp;
save(fullfile(fn_multi,'input.mat'),'input')
%% manual align
corrmap_norm{1} = uint8((corrmap{1}./max(corrmap{1}(:))).*255);
corrmap_norm{2} = uint8((corrmap{2}./max(corrmap{2}(:))).*255);
[input_points_1, base_points_1] = cpselect(fov_red{2},fov_red{1},'Wait', true);
[input_points_2, base_points_2] = cpselect(fov_norm{2},fov_norm{1},'Wait', true);
[input_points_3, base_points_3] = cpselect(corrmap_norm{2},corrmap_norm{1},'Wait', true);
input_points = [input_points_1; input_points_2; input_points_3];
base_points = [base_points_1; base_points_2; base_points_3];
fitGeoTAf = fitgeotrans(input_points(:,:), base_points(:,:),'affine'); 

data{3} = imwarp(data{2},fitGeoTAf, 'OutputView', imref2d(size(data{2}))); 
fov_avg{3} = mean(data{3},3);
fov_norm{3} = uint8((fov_avg{3}./max(fov_avg{3}(:))).*255);
corrmap{3} = double(imwarp(corrmap{2},fitGeoTAf, 'OutputView', imref2d(size(corrmap{2}))));
corrmap_norm{3} = uint8((corrmap{3}./max(corrmap{3}(:))).*255);
red_trans = (imwarp(fov_red{2},fitGeoTAf, 'OutputView', imref2d(size(fov_red{2}))));
fov_red{3} = uint8(red_trans);
df_max{3} = (imwarp(df_max{2},fitGeoTAf, 'OutputView', imref2d(size(df_max{2}))));

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
suptitle(mouse)
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
    
start = 1;
figure; movegui('center');
for icell = 1:nc
    if goodCells{id}(icell)
        % find best shift
        day1_cell_avg = fov_avg{id}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        day2_cell_avg = fov_avg{od(id)}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        [reg_avg, shift_avg] = shift_opt(day2_cell_avg,day1_cell_avg,2);
        r_avg = corr(reg_avg(:),day1_cell_avg(:));

        day1_cell_corr = corrmap{id}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        day2_cell_corr = corrmap{od(id)}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        [reg_corr, shift_corr] = shift_opt(day2_cell_corr,day1_cell_corr,2);
        r_corr = corr(reg_corr(:),day1_cell_corr(:));
        if red_ind{id}(icell)
            day1_red_avg = fov_red{id}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
            day2_red_avg = fov_red{od(id)}(...
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
        [max_val max_ind] = max([r_avg r_corr r_red]);
        if max_val>0.55 & (r_corr>0.4 || r_red>0.4)
            pass = true;
            switch max_ind
                case 1
                    shifts = shift_avg;
                case 2
                    shifts = shift_corr;
                case 3
                    shifts = shift_red;
            end
        else
            pass = false;
            shifts = nan;
        end


        mask = masks{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        mask(mask > icell | mask < icell) = 0;
        mask(mask == icell) = 1;
        mask_np = imCellNeuropil(mask, 3, 5);
        cellImageAlign(icell).center_yx = [yCenter(icell),xCenter(icell)];
        cellImageAlign(icell).d(1).mask = mask;
        cellImageAlign(icell).d(1).np_mask = mask_np;
        cellImageAlign(icell).d(1).avg_img = day1_cell_avg;
        cellImageAlign(icell).d(1).corr_img = day1_cell_corr;
        cellImageAlign(icell).d(1).red_img = day1_red_avg;
        cellImageAlign(icell).d(2).avg_img = reg_avg;
        cellImageAlign(icell).d(2).corr_img = reg_corr;
        cellImageAlign(icell).d(2).red_img = red_reg_avg;
        cellImageAlign(icell).r_avg = r_avg;
        cellImageAlign(icell).r_corr = r_corr;
        cellImageAlign(icell).r_red = r_red;
        cellImageAlign(icell).shifts = shifts;
        cellImageAlign(icell).pass = pass;

        if start>48
            figure;
            movegui('center')
            start = 1;
        end
        subplot(8,6,start)
        imagesc(day1_cell_corr)
        title(num2str(pass))
        subplot(8,6,start+1)
        imagesc(reg_corr)
        title(num2str(r_corr))
        subplot(8,6,start+2)
        imagesc(day1_cell_avg)
        title(num2str(pass))
        subplot(8,6,start+3)
        imagesc(reg_avg)
        title(num2str(r_avg))
        subplot(8,6,start+4)
        imagesc(day1_red_avg)
        title(num2str(pass))
        subplot(8,6,start+5)
        imagesc(red_reg_avg)
        title(num2str(r_red))

        start= start+6;

        % shift data, get tc, all days
        if pass
            data_day1_cell = data{id}(...
                yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
                xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1,:);
            tc_d1 = stackGetTimeCourses(data_day1_cell,mask);
            tc_d1_down = stackGetTimeCourses(stackGroupProject(data_day1_cell,5),mask);
            np_tc_d1 = stackGetTimeCourses(data_day1_cell,mask_np);
            np_tc_d1_down = stackGetTimeCourses(stackGroupProject(data_day1_cell,5),mask_np);
           
            
            data_day2_cell_shift = data{2}(...
                yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
                xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1,:);
            tc_d2 = stackGetTimeCourses(data_day2_cell_shift,mask);
            tc_d2_down = stackGetTimeCourses(stackGroupProject(data_day2_cell_shift,5),mask);
            np_tc_d2 = stackGetTimeCourses(data_day2_cell_shift,mask_np);
            np_tc_d2_down = stackGetTimeCourses(stackGroupProject(data_day2_cell_shift,5),mask_np);
            
            ii= 0.01:0.01:1;
            x_d1 = zeros(length(ii), 1);
            x_d2 = zeros(length(ii), 1);
            for i = 1:100
                x_d1(i,1) = skewness(tc_d1_down-tcRemoveDC(np_tc_d1_down*ii(i)));
                x_d2(i,1) = skewness(tc_d2_down-tcRemoveDC(np_tc_d2_down*ii(i)));
            end
            [max_skew ind_d1] =  max(x_d1,[],1);
            [max_skew ind_d2] =  max(x_d2,[],1);
            np_d1_w = 0.01*ind_d1;
            np_d2_w = 0.01*ind_d2;
            npSub_tc1 = tc_d1-(tcRemoveDC(np_tc_d1).*np_d1_w);
            npSub_tc2 = tc_d2-(tcRemoveDC(np_tc_d2).*np_d2_w);

            cellTCs_match{id} = [cellTCs_match{id} npSub_tc1];
            cellTCs_match{od(id)} = [cellTCs_match{od(id)} npSub_tc2];
        end
    end
end
match_ind = find([cellImageAlign.pass]);
red_ind_match{id} = ismember(match_ind,find(~isnan([cellImageAlign.r_red])));
red_ind_all = red_ind;

save(fullfile(fn_multi,'timecourses.mat'),'cellTCs_match', 'cellTCs_all', 'red_ind_all','red_ind_match')
save(fullfile(fn_multi,'multiday_alignment.mat'),'cellImageAlign','fitGeoTAf12','fitGeoTAf21', 'input_points','base_points', 'fov_avg', 'fov_norm','fov_red');
