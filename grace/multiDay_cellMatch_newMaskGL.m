%% get path names D2
ref_date = '210621';
date = '210625';
alignToRef = 1;
ImgFolder = strvcat('003');
% time = strvcat('1204');
mouse = 'i1803';
nrun = size(ImgFolder,1);
frame_rate = 15.5;
ref_str = 'runs-003';
run_str = catRunName(ImgFolder, nrun);
gl_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\2P_Imaging';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P';
behav_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';

%% load files 
maskD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_mask_cell.mat']));
TCs_D1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_TCs.mat']));
dfofD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_stimActFOV.mat']));
shiftsD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_reg_shifts.mat']));
pixelD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_pixel.mat']));
inputD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_input.mat']));

fov_avg{1} = shiftsD1.data_reg_avg;
dfmax{1} = dfofD1.data_dfof_max;
cellTCs_all{1} = TCs_D1.npSub_tc;
input_temp(1) = inputD1;
corrmap{1} = pixelD1.pix;
masks{1} = maskD1.mask_cell;

dfofD2 = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimActFOV.mat']));
shiftsD2 = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']));
pixelD2 = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pixel.mat']));
inputD2 = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']));

fov_avg{2} = shiftsD2.data_reg_avg;
dfmax{2} = dfofD2.data_dfof_max;
input_temp(2) = inputD2;
corrmap{2} = pixelD2.pix;

corrmap_norm{1} = uint8((corrmap{1}./max(corrmap{1}(:))).*255);
corrmap_norm{2} = uint8((corrmap{2}./max(corrmap{2}(:))).*255);

brightnessScaleFactor = 0.5;
fov_norm{1} = uint8((fov_avg{1}./max(fov_avg{1}(:))).*255);
fov_norm{1}(fov_norm{1} > (brightnessScaleFactor*255)) = brightnessScaleFactor*255;
fov_norm{2} = uint8((fov_avg{2}./max(fov_avg{2}(:))).*255);
fov_norm{2}(fov_norm{2} > (brightnessScaleFactor*255)) = brightnessScaleFactor*255;

%% red channel
day = 2;
irun = 1;
WL = '920';
ImgFolder = strvcat('001');
run = catRunName(ImgFolder, nrun);
imgMatFile = [ImgFolder '_000_000.mat'];
CD = fullfile(gl_fn, [mouse '\' date '_' mouse '\' ImgFolder(irun,:)]);
cd(CD);
load(imgMatFile);
% fprintf(['Reading run ' num2str(irun) '- ' num2str(info.config.frames) ' frames \r\n'])
data_temp = sbxread(imgMatFile(1,1:11),0,info.config.frames);
mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run]));

if size(data_temp,1)>1
data_rg = squeeze(data_temp(1,:,:,:));
data_rr = squeeze(data_temp(2,:,:,:));
[out data_g_reg] = stackRegister(data_rg,fov_avg{day});
[out2 data_r_reg] = stackRegister_MA(data_rr,[],[],out);
red = mean(data_r_reg,3);
greenChImg = mean(data_g_reg,3);
clear data_temp 
fov_red{day} = uint8((red./max(red(:))).*255);
end

%% red/green cells
if alignToRef
load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_multiday_alignment.mat']));
else
% size of cell box
clear input
close all
w=30;
h=30;
buf = 3;
np = 5;
% green channel
% get cell centroids

redGreenCells = struct;

cellPosition = regionprops(masks{1});
nc = length(cellPosition);

xCenter = cellfun(@(a) round(a(1)),{cellPosition.Centroid});
yCenter = cellfun(@(a) round(a(2)),{cellPosition.Centroid});

% index cells NOT too close to edge and NOT in black part of transformation
[ypix,xpix] = size(fov_avg{1});
goodCells = xCenter>(w/2) & xCenter<xpix-(w/2) & yCenter>(h/2) & yCenter<ypix-(h/2);

goodCells = goodCells & ...
    arrayfun(@(x) sum(sum(fov_norm{1}(masks{1}==x)))>0,1:nc);

    % fine register each cell    
mask_exp = zeros(size(fov_avg{1}));
mask_all = zeros(size(fov_avg{1}));

for icell = 1:nc
    if goodCells(icell)
        % find best shift
        day1_cell_avg = fov_norm{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        
        corr = corrmap_norm{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        
        day1_cell_max = dfmax{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        
        day1_red_avg = fov_red{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        
        mask = masks{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);

            pass = true;
            figure;
            movegui('center')
            start = 1;
            subplot(1,4,start)
            imagesc(day1_cell_avg);axis image
            hold on
            bound = cell2mat(bwboundaries(mask(:,:,1)));
            plot(bound(:,2),bound(:,1),'-','color','r','MarkerSize',2);axis image
            title('avg')
            subplot(1,4,start+1)
            imagesc(corr);axis image
            hold on
            plot(bound(:,2),bound(:,1),'-','color','r','MarkerSize',2);axis image
            title('corr')
            subplot(1,4,start+2)
            imagesc(day1_cell_max);axis image
            hold on
            plot(bound(:,2),bound(:,1),'-','color','r','MarkerSize',2);
            title('dfof')
            subplot(1,4,start+3)
            imagesc(day1_red_avg);axis image
            hold on
            plot(bound(:,2),bound(:,1),'-','color','r','MarkerSize',2);
            title('Red')
            
            prompt = 'Choose one: 1- good, 2-okay, 3- none: ';
            x = input(prompt);
            switch x
                case 1
                    pass = true;
                    faint = false;
                case 2
                    pass = false;
                    faint = true;
                case 3
                    pass = false;
                    faint = false;
            end
    if pass
        redGreenCells(icell).pass = pass;
        redGreenCells(icell).faint = false;  
    elseif faint
        redGreenCells(icell).pass = false;
        redGreenCells(icell).faint = faint;  
    else 
        redGreenCells(icell).pass = false;
        redGreenCells(icell).faint = false;  
    end
    close all
    end  
end

goodCells = find([redGreenCells.pass]);
okayCells = find([redGreenCells.faint]);
redCells = sort([goodCells okayCells]);
redChImg = fov_red{1};
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_multiday_alignment.mat']),'goodCells','okayCells','redCells','redChImg');
end
%% manual align
[input_points_1, base_points_1] = cpselect(fov_avg{2},fov_avg{1},'Wait', true);
[input_points_2, base_points_2] = cpselect(dfmax{2},dfmax{1},'Wait', true);
% [input_points_3, base_points_3] = cpselect(corrmap_norm{2},corrmap_norm{1},'Wait', true);
input_points = [input_points_1 input_points_2];
base_points = [input_points_1 base_points_2];
fitGeoTAf = fitgeotrans(input_points(:,:), base_points(:,:),'affine');
% fitGeoTAf_FOV = fitgeotrans(input_points_2(:,:), base_points_2(:,:),'affine');
% fitGeoTAf_pix = fitgeotrans(input_points_3(:,:), base_points_3(:,:),'affine');

data3 = imwarp(data_reg,fitGeoTAf, 'OutputView', imref2d(size(data_reg))); 
fov_avg{3} = mean(data3,3);
fov_norm{3} = uint8((fov_avg{3}./max(fov_avg{3}(:))).*255);

corrmap{3} = double(imwarp(corrmap{2},fitGeoTAf, 'OutputView', imref2d(size(corrmap{2}))));
corrmap_norm{3} = uint8((corrmap{3}./max(corrmap{3}(:))).*255);
% red_trans = (imwarp(fov_red{2},fitGeoTAf, 'OutputView', imref2d(size(fov_red{2}))));
% fov_red{3} = uint8(red_trans);
dfmax{3} = imwarp(dfmax{2},fitGeoTAf, 'OutputView', imref2d(size(dfmax{2})));
redChImg = imwarp(fov_red{2},fitGeoTAf, 'OutputView', imref2d(size(fov_red{2})));

figure;
filler = zeros(size(fov_norm{1}));
imshow(cat(3,dfmax{1},dfmax{3},filler))
title('Overlay with all points')
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['FOV_alignment.pdf']),'-dpdf','-bestfit')

% figure;colormap gray
% movegui('center')
% subplot 221
% imshow(fov_red{1}); title('Day 1 Red Channel')
% subplot 222
% imshow(fov_red{3}); title('Transformed Day 2 Red Channel')
% subplot 223
% filler = zeros(size(fov_red{1}));
% imshow(cat(3,fov_red{1},fov_red{3},filler))
% title('Overlay')
% print(fullfile(fn_multi,'red channel manual alignment'),'-dpdf','-fillpage')

%% cell-by-cell correlation
clear input
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
figure; movegui('center');
for icell = 1:nc
    if goodCells(icell)
        % find best shift
        day1_mask = masks{1}(... %%need original D1 mask trans to D2
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
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
%         if red_ind{1}(icell)
%             day1_red_avg = fov_red{1}(...
%             yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
%             xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
%             day2_red_avg = fov_red{3}(...
%             yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
%             xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
%             [red_reg_avg, shift_red] = shift_opt(double(day2_red_avg),double(day1_red_avg),2);
%             r_red = corr(red_reg_avg(:),double(day1_red_avg(:)));
%         else
            day1_red_avg = nan;
            day2_red_avg = nan;
            red_reg_avg = nan;
            r_red = nan;
%         end
        
        [max_val max_ind] = max([r_avg r_max r_corr r_red]);
        if max_val>0.55 & (r_corr>0.4 || r_avg>0.4 || r_max>0.4)
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
%             if red_ind{1}(icell)
%                 imagesc(day1_red_avg)
%                  title('Red')
%             else
                imagesc(day1_cell_avg)
                title('Avg')
%             end
            subplot(3,2,start+3)
%             if red_ind{1}(icell)
%                 imagesc(red_reg_avg)
%                 title(num2str(r_red))
%             else
                imagesc(reg_avg)
                title(num2str(r_avg))
%             end
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
%                     if red_ind{1}(icell)
%                         img_select = fov_red{3};
%                     else
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

        cellImageAlign(icell).center_yx = [yCenter(icell),xCenter(icell)];
        cellImageAlign(icell).d(1).avg_img = day1_cell_avg;
        cellImageAlign(icell).d(1).corr_img = day1_cell_corr;
        cellImageAlign(icell).d(1).red_img = day1_red_avg;
        cellImageAlign(icell).d(1).max_img = day1_cell_max;
        cellImageAlign(icell).d(2).avg_img = reg_avg;
        cellImageAlign(icell).d(2).corr_img = reg_corr;
        cellImageAlign(icell).d(2).red_img = red_reg_avg;
        cellImageAlign(icell).d(2).max_img = reg_max;
%         cellImageAlign(icell).r_avg = r_avg;
%         cellImageAlign(icell).r_corr = r_corr;
%         cellImageAlign(icell).r_red = r_red;
        cellImageAlign(icell).shifts = shifts;
        cellImageAlign(icell).pass = pass;

        

        % shift data, get tc, all days
        if pass
            mask_data_temp = img_select;
            mask_data_temp(find(mask_exp >= 1)) = 0;
            mask_data_square = mask_data_temp(...
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
                mask_exp = imCellBuffer(mask_temp,3)+mask_temp;
            else
                cellImageAlign(icell).pass = false;
                temp_mask = zeros(size(day1_mask));
                temp_mask(find(day1_mask==icell)) = 1;
                bwout_full = zeros(size(fov_avg{1}));
                bwout_full(...
                yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
                xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1) = temp_mask*icell;
                mask_all = mask_all+bwout_full;
                mask_temp = mask_all;
                mask_temp(find(mask_all>=1)) = 1;
                mask_exp = imCellBuffer(mask_temp,3)+mask_temp;
            end
        close all
    else
        cellImageAlign(icell).pass = false;
        cellImageAlign(icell).r_red = 0;
        temp_mask = zeros(size(day1_mask));
        temp_mask(find(day1_mask==icell)) = 1;
        bwout_full = zeros(size(fov_avg{1}));
        bwout_full(...
        yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
        xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1) = temp_mask*icell;
        mask_all = mask_all+bwout_full;
        mask_temp = mask_all;
        mask_temp(find(mask_all>=1)) = 1;
        mask_exp = imCellBuffer(mask_temp,3)+mask_temp;
        end
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
        temp_mask(find(day1_mask==icell)) = 1;
        bwout_full = zeros(size(fov_avg{1}));
        bwout_full(...
        a:b,...
        c:d) = temp_mask*icell;
        cellImageAlign(icell).day2_mask = bwout_full;
        mask_all = mask_all+bwout_full;
        mask_temp = mask_all;
        mask_temp(find(mask_all>=1)) = 1;
        mask_exp = imCellBuffer(mask_temp,3)+mask_temp;
    end
end
%mask_cell = bwlabel(mask_all);
mask_cell = mask_all;
mask_np = imCellNeuropil(mask_cell, 3, 5);
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
match = find([cellImageAlign.pass]);
[match_ind,red_ind] = intersect(match,redCells);
cellTCs_match{1} = cellTCs_all{1}(:,match_ind);

% x = [size(cell2mat(cellTCs_all),2) length(redCells) length(match) length(match_ind)];
% figure;bar(x)
% cellnames = {'D1 cells','D1 red cells','D2 match cells','D2 red/match cells'};
% set(gca,'xticklabel',cellnames,'FontSize',12)

%new TCs
data_tc = stackGetTimeCourses(data3, mask_cell);
[nFrames nCells] = size(data_tc);
down = 5;
data_reg_down  = stackGroupProject(data3,down);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);

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
npSub_tc = data_tc(:,:)-bsxfun(@times,tcRemoveDC(np_tc(:,:)),np_w);
 
% cellTCs_match{2} = cellTCs_match{2}(:,red_ind);
cellTCs_match{2} = npSub_tc(:,:);

% red_ind_match = ismember(match_ind,find(~isnan([cellImageAlign.r_red])));
% red_ind_all = red_ind;

save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']),'cellTCs_match', 'cellTCs_all','match_ind')
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_multiday_alignment.mat']),'cellImageAlign','fitGeoTAf', 'input_points','base_points', 'fov_avg', 'fov_norm','dfmax','corrmap','redChImg');
clear data_reg_down

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
        data_mat = zeros(nOn+nOff, nCells, ntrials);
        for itrial = 1:ntrials
            data_mat(:,:,itrial) = npSub_tc(1+((itrial-1).*(nOn+nOff)):(itrial.*(nOn+nOff)),:);
        end
        data_f = mean(data_mat(nOff/2:nOff,:,:),1);
    else
        data_mat = zeros(100, nCells, ntrials-1);
        for itrial = 1:ntrials-1
            data_mat(:,:,itrial) = npSub_tc(((itrial-1)*(nOn+nOff))+71:170+((itrial-1)*(nOn+nOff)),:);
        end
        data_f = mean(data_mat(1:50,:,:),1);
    end
    data_df = bsxfun(@minus, data_mat, data_f);
    data_dfof = bsxfun(@rdivide, data_df, data_f);
    
    ndir = length(Dirs);
    [n, n2] = subplotn(nCells);
    h_dir = zeros(nCells, ndir);
    p_dir = zeros(nCells, ndir);
    base_win = 50:60;
    resp_win = 70:90;
    base = squeeze(mean(data_dfof(base_win,:,:),1));
    resp = squeeze(mean(data_dfof(resp_win,:,:),1));
    resp_mat = squeeze(mean(data_mat(resp_win,:,:),1));
%     figure;scatter(mean(resp,2),mean(resp_mat,2));xlabel('dfof resp');ylabel('f resp')
    dir_resp = zeros(nCells,ndir);
    [x y] = ttest(resp', base', 'tail','right');
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
    
    good_ind = unique([find(x)'; find(sum(h_dir,2)>0); find(sum(h_ori,2)>0)]);
    print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuning.pdf']),'-dpdf')
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_trialData.mat']),'data_dfof','resp_mat','max_dir','h_dir', 'h_ori', 'max_ori','good_ind')
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

% if mean(npSub_tc1,1) = 0
% [avgResponseEaOri,semResponseEaOri,vonMisesFitAllCellsAllBoots,fitReliability,R_square,tuningTC] = 0;
% else
[avgResponseEaOri,semResponseEaOri,vonMisesFitAllCellsAllBoots,fitReliability,R_square,tuningTC] = ...
    getOriTuningLG(data_tc_down,input,tuningDownSampFactor);
    vonMisesFitAllCells = squeeze(vonMisesFitAllCellsAllBoots(:,1,:));
% end

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
prefOri_bootdiff = abs(prefOri(2:end,:)-prefOri(1,:));
prefOri_bootdiff(find(prefOri_bootdiff>90)) = 180-prefOri_bootdiff(find(prefOri_bootdiff>90));
ind_theta90 = find(prctile(prefOri_bootdiff,90,1)<22.5);
edges = [0 22.5:45:180]; 
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
