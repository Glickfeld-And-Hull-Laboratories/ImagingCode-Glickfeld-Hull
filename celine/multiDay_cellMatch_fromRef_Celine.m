 clear all; clear global; close all
clc
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories
eval(ds)
doGreenOnly = false;
doCorrImg = true;


day_id(1) = 330; %enter the post-DART day ID here
day_id(2) = expt(day_id(1)).multiday_matchdays;

ref_day = 329;

nd = length(day_id);
brightnessScaleFactor = 0.3;
mouse = expt(day_id(1)).mouse;

if strcmp(expt(day_id(1)).data_loc,'lindsey')
    root = fullfile(rc.data,mouse);
elseif strcmp(expt(day_id(1)).data_loc,'ashley')
    root = fullfile(rc.ashleyData,mouse,'two-photon imaging');
elseif strcmp(expt(day_id(1)).data_loc,'tammy')
        root = rc.tammyData;
        expDate = expt(day_id(1)).date;
        runFolder = expt(day_id(1)).contrastxori_runs;
        dat = 'data-i';
elseif strcmp(expt(day_id(1)).data_loc,'ACh')
        root = rc.achData;
        expDate = expt(day_id(1)).date;
        runFolder = expt(day_id(1)).contrastxori_runs;
        dat = 'data-i';
end

fnout = fullfile(rc.achAnalysis,mouse);

if expt(day_id(2)).multiday_timesincedrug_hours>0
    dart_str = [expt(day_id(2)).drug '_' num2str(expt(day_id(1)).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end

fn_multi = fullfile(rc.achAnalysis,mouse,['multiday_' dart_str]);


dart_str_ref = [expt(ref_day).drug '_' num2str(expt(ref_day).multiday_timesincedrug_hours) 'Hr'];
fnref = fullfile(rc.achAnalysis,mouse,['multiday_' dart_str_ref]);


mkdir(fn_multi)


load(fullfile(fnref,'timecourses.mat')); %load the .mat that has the match inds and red inds 
%then clear the tc dataframes from that becuase we want to make new ones
%now
clear cellTCs_all cellTCs_match data_tc

data = cell(1,nd);
fov_avg = cell(1,nd);
fov_norm = cell(1,nd);
fov_red = cell(1,nd);
align_fov = cell(1,nd);
corrmap = cell(1,nd);
dfmax = cell(1,nd);
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
        cd(fullfile(root,mouse,expDate, imgFolder))
        load(fName)
        if nrun == 1
            load(fullfile(fnout,expDate,imgFolder,'regOuts&Img.mat'))
            nFrames = size(outs,1);
            runFolder = imgFolder;
        else
            nFrames = info.config.nFrames;
            runFolder = [runFolder '_' imgFolder];
        end
       
        fprintf(['Loading day ' num2str(id) ' data \n Mouse: ' expt(day_id(id)).mouse ' Date: ' expt(day_id(id)).date])
        data_temp = sbxread(fName(1,1:11),0,nFrames);
        data_g = cat(3, data_g, squeeze(data_temp(1,:,:,:)));
        clear data_temp
        out_all = [out_all; outs];
    end
    data_g = data_g(:,:,1:nFrames);%if one of the days is the wrong frame number
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
    clear mask_cell mask_cell_red_mask_label mask_np
    dfmax{id} = data_dfof(:,:,1);
%     corrmap{id} = data_dfof(:,:,end);
%     masks{id} = mask_cell;
%     maskNP{id} = mask_np;
%     red_ind{id} = mask_label;
    load(fullfile(fnout,expDate,runFolder,'TCs.mat'))
    cellTCs_all{id} = npSub_tc;
    load(fullfile(fnout,expDate,runFolder,'input.mat'))
    input_temp(id) = input;
end
input = input_temp;
save(fullfile(fn_multi,'input.mat'),'input')
clear input
%% manual align
corrmap_norm{1} = uint8((corrmap{1}./max(corrmap{1}(:))).*255);
corrmap_norm{2} = uint8((corrmap{2}./max(corrmap{2}(:))).*255);
if exist(fullfile(fn_multi,'multiday_alignment.mat')) %first check for an alignment file alrady existing for this data
    load(fullfile(fn_multi,'multiday_alignment.mat'))
elseif exist(fullfile(fnref,'multiday_alignment.mat')) %if one does not exist, pull from the reference folder
    load(fullfile(fnref,'multiday_alignment.mat'))
else
    [input_points_1, base_points_1] = cpselect(fov_red{2},fov_red{1},'Wait', true);
    [input_points_2, base_points_2] = cpselect(fov_norm{2},fov_norm{1},'Wait', true);
    [input_points_3, base_points_3] = cpselect(corrmap_norm{2},corrmap_norm{1},'Wait', true);
    input_points = [input_points_1; input_points_2; input_points_3];
    base_points = [base_points_1; base_points_2; base_points_3];
end
%%
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

%% copy masks from ref

 
figure; movegui('center')
subplot(2,2,1)
imagesc(masks{1})
title('Bsln day masks')
subplot(2,2,2)
imagesc(masks{2})
title('Matched day masks after transform')
print(fullfile(fn_multi,'masksAfterTransform.pdf'),'-dpdf','-fillpage')

%%

 %for the reference day, use the timecourses from when this dataset was
 %segmented individually, extracting the matched cells by index
cellTCs_match{1} = cellTCs_all{1}(:,match_ind);

%for the day being matched to the reference, apply the new masks
data_tc = stackGetTimeCourses(data{3}, masks{2});
[nFrames nCells] = size(data_tc);
down = 5;
data_reg_down  = stackGroupProject(data{3},down);
data_tc_down = stackGetTimeCourses(data_reg_down, masks{2});

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

%%
save(fullfile(fn_multi,'timecourses.mat'),'cellTCs_match', 'cellTCs_all', 'red_ind_all','red_ind_match','match_ind')
save(fullfile(fn_multi,'multiday_alignment.mat'),'cellImageAlign','fitGeoTAf', 'input_points','base_points', 'fov_avg', 'fov_norm','fov_red','dfmax','corrmap','masks','mask_np');

clear data_reg_down data