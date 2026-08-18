clear all; clear global; close all
clc
ds = 'DART_behavior_ExptList'; %dataset info
dataStructLabels = {'behavior_runs'};
rc = behavConstsDART; %directories
eval(ds)
doGreenOnly = false;
doCorrImg = true;

%to use the post-DART timepoint as the template

day_id(1) = 19; %enter the refrence day ID here
day_id(2) = 17;%expt(day_id(1)).multiday_matchdays;

experimentFolder = 'SST_behavior';


nd = length(day_id);
brightnessScaleFactor = 0.3;
mouse = expt(day_id(1)).mouse;


if computer == 'GLNXA64'
    isilonName =  '/home/cc735@dhe.duke.edu/GlickfeldLabShare';
    database = fullfile('/All_Staff/home/ACh/Data/2p_data');
    base = fullfile('/All_Staff/home/ACh/Analysis/2p_analysis',experimentFolder);
    beh_prefix = strcat(isilonName,'/All_Staff/Behavior/Data/data-');
elseif string(hostname) == 'NB-NUKE'
    isilonName = 'Z:/All_Staff';
    base = fullfile(isilonName,'/home/ACh/Analysis/2p_analysis',experimentFolder);
    database = fullfile(isilonName,'/home/ACh/Data/2p_data');
else
    isilonName = 'duhs-user-nc1.dhe.duke.edu/dusom_glickfeldlab';
    base = fullfile('/All_Staff/home/ACh/Analysis/2p_analysis',experimentFolder);
    database = fullfile('/All_Staff/home/ACh/Data/2p_data/');
   
   beh_prefix = strcat('/All_Staff/Behavior/Data/data-');
end

fnout = fullfile(base,mouse);

if expt(day_id(1)).multiday_timesincedrug_hours>0
    dart_str = [expt(day_id(2)).drug '_' num2str(expt(day_id(1)).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end

fn_multi = fullfile(base,mouse,['multiday_' dart_str '_vsD2']);
mkdir(fn_multi)
data = cell(1,nd);
fov_avg = cell(1,nd);
fov_norm = cell(1,nd);
fov_red = cell(1,nd);
align_fov = cell(1,nd);
corrmap = cell(1,nd);
dfmax = cell(1,nd);
masks = cell(1,nd);
red_masks = cell(1,nd);
maskNP = cell(1,nd);
red_ind = cell(1,nd);
cellTCs_all = cell(1,nd);

days_text = strcat('Reference day: ', string(day_id(1)), ' matched day: ', string(day_id(2)))
cd(fn_multi)
fid = fopen('sessionsMatched.txt','wt');
fprintf(fid, days_text);
fclose(fid);
%% load all data 
runFolder = [];
for id = 1:nd 
    clear global
    expDate = expt(day_id(id)).date;
    runs = eval(['expt(day_id(' num2str(id) ')).' cell2mat(dataStructLabels) ]);
    nrun = length(runs);
    out_all = [];
    data_g = [];
    for irun = 1:nrun
        imgFolder = runs{irun};
        fName = [imgFolder '_000_000'];
        cd(fullfile(database, mouse,expDate, imgFolder))
        load(fName)
        if nrun == 1
            load(fullfile(fnout,expDate,imgFolder,'regOuts&Img.mat'))
            nframes = size(out,1);
            runFolder = imgFolder;
        else
            nframes = info.config.nframes;
            runFolder = [runFolder '_' imgFolder];
        end
        fprintf(['Loading day ' num2str(id) ' data \n Mouse: ' expt(day_id(id)).mouse ' Date: ' expt(day_id(id)).date])
        data_temp = sbxread(fName(1,1:11),0,nframes);
        data_g = cat(3, data_g, squeeze(data_temp(1,:,:,:)));
        clear data_temp
        out_all = [out_all; out];
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
    dfmax{id} = data_dfof(:,:,end);
    corrmap{id} = data_dfof(:,:,4);
    masks{id} = mask_cell;
    red_masks{id} = mask_cell_red;
    maskNP{id} = mask_np;
    red_ind{id} = mask_label;
    load(fullfile(fnout,expDate,runFolder,'TCs.mat'))
    cellTCs_all{id} = npSub_tc;
    load(fullfile(fnout,expDate,runFolder,'input.mat'))
    input_temp(id) = input;
    load(fullfile(fnout,expDate,runFolder,'stimuli.mat'))
    stimData{id}.cons = tContrast;
    stimData{id}.block = b1;
    stimData{id}.SIx = SIx;
    stimData{id}.MIx = MIx;
    stimData{id}.FIx = FIx;
    stimData{id}.cTarget = cTarget;
    stimData{id}.tquit = tquit;
    load(fullfile(fnout,expDate,runFolder, 'stimData.mat'))
    stimData{id}.base_win = base_win;
    stimData{id}.resp_win = resp_win;    
end
input = input_temp;
save(fullfile(fn_multi,'input.mat'),'input')
save(fullfile(fn_multi,'stimData.mat'),'stimData')
clear input
%% manual align
corrmap_norm{1} = uint8((corrmap{1}./max(corrmap{1}(:))).*255);
corrmap_norm{2} = uint8((corrmap{2}./max(corrmap{2}(:))).*255);
fov_norm{1} = uint8((fov_avg{1}./max(fov_avg{1}(:))).*255);
fov_norm{2} = uint8((fov_avg{2}./max(fov_avg{2}(:))).*255);

if exist(fullfile(fn_multi,'multiday_alignment.mat'))
    load(fullfile(fn_multi,'multiday_alignment.mat'),'fitGeoTAf','input_points','base_points')
else
    [input_points_1, base_points_1] = cpselect(fov_red{2},fov_red{1},'Wait', true);
    [input_points_2, base_points_2] = cpselect(fov_norm{2},fov_norm{1},'Wait', true);
    [input_points_3, base_points_3] = cpselect(corrmap_norm{2},corrmap_norm{1},'Wait', true);
    input_points = [input_points_1; input_points_2; input_points_3];
    base_points = [base_points_1; base_points_2; base_points_3];
end

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

%% target analysis
nCells= length(match_ind);
figure;
for id = 1:nd
    cTarget = stimData{id}.cTarget;
    nTrials = length(cTarget);
    data_targ = nan(50,nCells,nTrials);
    for itrial = 1:nTrials
        if ~isnan(cTarget(itrial))
            if cTarget(itrial)+50 <= size(cellTCs_match{id},1)
                data_targ(:,:,itrial) = cellTCs_match{id}(cTarget(itrial)-20:cTarget(itrial)+29,:);
            end
        end
    end
    data_f_targ = nanmean(data_targ(1:20,:,:),1);
    data_dfof_targ{id} = bsxfun(@rdivide,bsxfun(@minus,data_targ,data_f_targ),data_f_targ);

    frameRateHz = 30;

    tt = [-20:29].*(1000/frameRateHz);
    subplot(2,1,id)
    plot(tt,squeeze(nanmean(mean(data_dfof_targ{id},2),3)));
    vline((stimData{id}.base_win-20).*(1000/frameRateHz))
    vline((stimData{id}.resp_win-20).*(1000/frameRateHz))
end
print(fullfile(fn_multi, 'avgTC.pdf'),'-dpdf','-bestfit')

%% Contrast response
h = cell(1,nd);
p = cell(1,nd);
data_dfof_con = cell(1,nd);
con_resp_mat = cell(1,nd);
data_dfof_con_block = cell(1,nd);
con_block_resp_mat = cell(1,nd);
trialn = cell(1,nd);
base_win = stimData{1}.base_win;
for id = 1:nd
    b1 = stimData{id}.block;
    tContrast = stimData{id}.cons;
    cons = unique(tContrast);
    nCon = length(cons);
    SIx = stimData{id}.SIx;
    MIx = stimData{id}.MIx;
    tquit = stimData{id}.tquit;
    nblock = length(unique(b1));
    con_resp_mat{id} = zeros(nCon,nCells,2);
    con_block_resp_mat{id} = zeros(nCon,nCells,nblock,2,2); %ncon,ncells,nblock,hit/miss,mean/sem
    data_dfof_con{id} = zeros(size(data_dfof_targ{id},1), nCells, nCon);
    data_dfof_con_block{id} = zeros(size(data_dfof_targ{id},1), nCells, nCon, nblock,2); % nframes,ncells,ncon,nblock,hit/miss
    h{id} = zeros(nCon,nCells);
    p{id} = zeros(nCon,nCells);
    trialn{id} = zeros(nCon, nblock,2); %ncon,nblock,hit/miss
    [max_val max_ind] = max(mean(nanmean(data_dfof_targ{id},3),2),[],1);
    resp_win = max_ind-1:max_ind+1;
    for iCon = 1:nCon
        ind_con = find(tContrast == cons(iCon));
        data_dfof_con{id}(:,:,iCon) = nanmean(data_dfof_targ{id}(:,:,ind_con),3);
        con_resp_mat{id}(iCon,:,1) = mean(nanmean(data_dfof_targ{id}(resp_win,:,ind_con)-data_dfof_targ{id}(base_win,:,ind_con),3),1);
        con_resp_mat{id}(iCon,:,2) = nanstd(nanmean(data_dfof_targ{id}(resp_win,:,ind_con)-data_dfof_targ{id}(base_win,:,ind_con),1),[],3)./sqrt(length(ind_con));
        for iCell = 1:nCells
            [h{id}(iCon,iCell), p{id}(iCon,iCell)] = ttest(mean(permute(data_dfof_targ{id}(resp_win,iCell,ind_con),[1 3 2]),1),mean(permute(data_dfof_targ{id}(base_win,iCell,ind_con),[1 3 2]),1),'tail','right','alpha',0.05./(nCon));
        end
        for ib = 1:nblock
            ind_s = intersect(find(SIx),intersect(1:tquit,intersect(ind_con,find(b1==ib-1))));
            ind_m = intersect(find(MIx),intersect(1:tquit,intersect(ind_con,find(b1==ib-1))));
            data_dfof_con_block{id}(:,:,iCon,ib,1) = nanmean(data_dfof_targ{id}(:,:,ind_s),3);
            con_block_resp_mat{id}(iCon,:,ib,1,1) = mean(nanmean(data_dfof_targ{id}(resp_win,:,ind_s)-data_dfof_targ{id}(base_win,:,ind_s),3),1);
            con_block_resp_mat{id}(iCon,:,ib,1,2) = nanstd(nanmean(data_dfof_targ{id}(resp_win,:,ind_s)-data_dfof_targ{id}(base_win,:,ind_s),1),[],3)./sqrt(length(ind_s));
            data_dfof_con_block{id}(:,:,iCon,ib,2) = nanmean(data_dfof_targ{id}(:,:,ind_m),3);
            con_block_resp_mat{id}(iCon,:,ib,2,1) = mean(nanmean(data_dfof_targ{id}(resp_win,:,ind_m)-data_dfof_targ{id}(base_win,:,ind_m),3),1);
            con_block_resp_mat{id}(iCon,:,ib,2,2) = nanstd(nanmean(data_dfof_targ{id}(resp_win,:,ind_m)-data_dfof_targ{id}(base_win,:,ind_m),1),[],3)./sqrt(length(ind_m));
            trialn{id}(iCon,ib,1) = length(ind_s);
            trialn{id}(iCon,ib,2) = length(ind_m);
        end
    end
end

good_ind = unique([find(sum(h{1},1)) find(sum(h{2},1))]);
good_ind_red = intersect(good_ind,find(red_ind_match));
good_ind_green = intersect(good_ind,find(~red_ind_match));
%compare all trials all blocks
figure;
colors = {'b','k'};
line_style1 = '';
for id = 1:nd
    start = 1;
    for iCon = 1:nCon
        subplot(nCon,2,start)
        shadedErrorBar(tt, mean(data_dfof_con{id}(:,good_ind_green,iCon),2),std(data_dfof_con{id}(:,good_ind_green,iCon),[],2)./sqrt(length(good_ind_green)),[line_style1, colors{id}])
        ylim([-0.02 .15])
        ylabel('dF/F')
        xlabel('Time from target (ms)')
        hold on
        subplot(nCon,2,start+1)
        shadedErrorBar(tt, mean(data_dfof_con{id}(:,good_ind_red,iCon),2),std(data_dfof_con{id}(:,good_ind_red,iCon),[],2)./sqrt(length(good_ind_red)),[line_style1, colors{id}])
        hold on
        ylabel('dF/F')
        xlabel('Time from target (ms)')
        %title(num2str(cons(iCon)))
        ylim([-0.02 .15])
        start = start+2;
    end
end
suptitle(['All trials: D1- black; DN- blue; Pyr- n = ' num2str(length(good_ind_green)) ' left; SST- n = ' num2str(length(good_ind_red)) ' right'])
print(fullfile(fn_multi, 'contrastRespTC_allTrials.pdf'),'-dpdf', '-bestfit')

%compare success only within block across session

for ib = 1:2
    figure;
    for id = 1:nd
        start = 1;
        for iCon = 1:nCon
            subplot(nCon,2,start)
            shadedErrorBar(tt, mean(data_dfof_con_block{id}(:,good_ind_green,iCon,ib,1),2),std(data_dfof_con_block{id}(:,good_ind_green,iCon,ib,1),[],2)./sqrt(length(good_ind_green)),[line_style1, colors{id}])
            ylim([-0.02 .15])
            xlim([-800 800])
            ylabel('dF/F')
            xlabel('Time from target (ms)')
            hold on
            subplot(nCon,2,start+1)
            shadedErrorBar(tt, mean(data_dfof_con_block{id}(:,good_ind_red,iCon,ib,1),2),std(data_dfof_con_block{id}(:,good_ind_red,iCon,ib,1),[],2)./sqrt(length(good_ind_red)),[line_style1, colors{id}])
            hold on
            ylabel('dF/F')
            xlabel('Time from target (ms)')
            %title(num2str(cons(iCon)))
            ylim([-0.02 .15])
            xlim([-800 800])
            start = start+2;
        end
    end
    suptitle(['Block ' num2str(ib-1) ': D1- black; DN- blue; Pyr- n = ' num2str(length(good_ind_green)) ' left; SST- n = ' num2str(length(good_ind_red)) ' right'])
    print(fullfile(fn_multi, ['contrastRespTC_Block' num2str(ib-1) '.pdf']),'-dpdf', '-bestfit')
end

colors = {'k','r'};
days = {'N','1'};
for id = 1:2
    figure;
    for ib = 1:2
        start = 1;
        for iCon = 1:nCon
            subplot(nCon,2,start)
            shadedErrorBar(tt, mean(data_dfof_con_block{id}(:,good_ind_green,iCon,ib,1),2),std(data_dfof_con_block{id}(:,good_ind_green,iCon,ib,1),[],2)./sqrt(length(good_ind_green)),[line_style1, colors{ib}])
            ylim([-0.02 .15])
            xlim([-800 800])
            ylabel('dF/F')
            xlabel('Time from target (ms)')
            hold on
            subplot(nCon,2,start+1)
            shadedErrorBar(tt, mean(data_dfof_con_block{id}(:,good_ind_red,iCon,ib,1),2),std(data_dfof_con_block{id}(:,good_ind_red,iCon,ib,1),[],2)./sqrt(length(good_ind_red)),[line_style1, colors{ib}])
            hold on
            ylabel('dF/F')
            xlabel('Time from target (ms)')
            %title(num2str(cons(iCon)))
            ylim([-0.02 .15])
            xlim([-800 800])
            start = start+2;
        end
    end
    suptitle(['Day ' days{id} ': B0- black; B1- red; Pyr- n = ' num2str(length(good_ind_green)) ' left; SST- n = ' num2str(length(good_ind_red)) ' right'])
    print(fullfile(fn_multi, ['contrastRespTC_Day' days{id} '.pdf']),'-dpdf', '-bestfit')
end

colors = {'k','r'};
figure; 
for id = 1:nd
    for ib = 1:nblock
        subplot(2,2,id)
        errorbar(cons, mean(con_block_resp_mat{id}(:,good_ind_green,ib,1,1),2), std(con_block_resp_mat{id}(:,good_ind_green,ib,1,1),[],2)./sqrt(length(good_ind_green)),colors{ib})
        title(['D' days{id} ' HTP-'])
        xlabel('Contrast')
        ylim([-0.02 0.15])
        hold on
        subplot(2,2,2+id)
        errorbar(cons, mean(con_block_resp_mat{id}(:,good_ind_red,ib,1,1),2), std(con_block_resp_mat{id}(:,good_ind_red,ib,1,1),[],2)./sqrt(length(good_ind_red)),colors{ib})
        title(['D' days{id} ' HTP+'])
        ylim([-0.02 0.15])
        xlabel('Contrast')
        hold on
    end
end
suptitle(['B0- black; B1- red; Pyr- n = ' num2str(length(good_ind_green)) ' top; SST- n = ' num2str(length(good_ind_red)) ' bottom'])
print(fullfile(fn_multi, ['contrastResp_byBlock.pdf']),'-dpdf', '-bestfit')

figure; 
colors = {'b','k'};
for ib = 1:nblock
    for id = 1:nd
        subplot(2,2,ib)
        errorbar(cons, mean(con_block_resp_mat{id}(:,good_ind_green,ib,1,1),2), std(con_block_resp_mat{id}(:,good_ind_green,ib,1,1),[],2)./sqrt(length(good_ind_green)),colors{id})
        title(['B' num2str(ib-1) ' HTP-'])
        xlabel('Contrast')
        ylim([-0.02 0.15])
        hold on
        subplot(2,2,2+ib)
        errorbar(cons, mean(con_block_resp_mat{id}(:,good_ind_red,ib,1,1),2), std(con_block_resp_mat{id}(:,good_ind_red,ib,1,1),[],2)./sqrt(length(good_ind_red)),colors{id})
        title(['B' num2str(ib-1) ' HTP+'])
        ylim([-0.02 0.15])
        xlabel('Contrast')
        hold on
    end
end
suptitle(['D1- black; DN- blue; Pyr- n = ' num2str(length(good_ind_green)) ' top; SST- n = ' num2str(length(good_ind_red)) ' bottom'])
print(fullfile(fn_multi, ['contrastResp_byDay.pdf']),'-dpdf', '-bestfit')


figure;
for id = 1:2
    subplot(2,2,id)
    scatter(mean(con_block_resp_mat{id}(:,good_ind_green,1,1,1),1),mean(con_block_resp_mat{id}(:,good_ind_green,2,1,1),1))
    hold on
    y_sem = std(mean(con_block_resp_mat{id}(:,good_ind_green,2,1,1),1),[],2)./sqrt(length(good_ind_green));
    x_sem = std(mean(con_block_resp_mat{id}(:,good_ind_green,1,1,1),1),[],2)./sqrt(length(good_ind_green));
    errorbar(mean(mean(con_block_resp_mat{id}(:,good_ind_green,1,1,1),1),2),mean(mean(con_block_resp_mat{id}(:,good_ind_green,2,1,1),1),2),x_sem,x_sem,y_sem,y_sem,'ok');
    [h_green p_green] = ttest(mean(con_block_resp_mat{id}(:,good_ind_green,1,1,1),1),mean(con_block_resp_mat{id}(:,good_ind_green,2,1,1),1));
    xlabel('Mean resp B0')
    ylabel('Mean resp B1')
    xlim([-0.05 0.3])
    ylim([-0.05 0.3])
    axis square
    refline(1)
    title(['D' days{id} ' HTP-; p = ' num2str(chop(p_green,2))])
    subplot(2,2,2+id)
    scatter(mean(con_block_resp_mat{id}(:,good_ind_red,1,1,1),1),mean(con_block_resp_mat{id}(:,good_ind_red,2,1,1),1))
    hold on
    y_sem = std(mean(con_block_resp_mat{id}(:,good_ind_red,2,1,1),1),[],2)./sqrt(length(good_ind_red));
    x_sem = std(mean(con_block_resp_mat{id}(:,good_ind_red,1,1,1),1),[],2)./sqrt(length(good_ind_red));
    errorbar(mean(mean(con_block_resp_mat{id}(:,good_ind_red,1,1,1),1),2),mean(mean(con_block_resp_mat{id}(:,good_ind_red,2,1,1),1),2),x_sem,x_sem,y_sem,y_sem,'ok')
    [h_red p_red] = ttest(mean(con_block_resp_mat{id}(:,good_ind_red,1,1,1),1),mean(con_block_resp_mat{id}(:,good_ind_red,2,1,1),1));
    xlabel('Mean resp B0')
    ylabel('Mean resp B1')
    xlim([-0.02 0.15])
    ylim([-0.02 0.15])
    axis square
    refline(1)
    title(['D' days{id} ' HTP+; p = ' num2str(chop(p_red,2))])
end
print(fullfile(fn_multi, ['contrastResp_byBlock_scatter.pdf']),'-dpdf', '-bestfit')

figure;
for ib = 1:2
    subplot(2,2,ib)
    scatter(mean(con_block_resp_mat{2}(:,good_ind_green,ib,1,1),1),mean(con_block_resp_mat{1}(:,good_ind_green,ib,1,1),1))
    hold on
    y_sem = std(mean(con_block_resp_mat{1}(:,good_ind_green,ib,1,1),1),[],2)./sqrt(length(good_ind_green));
    x_sem = std(mean(con_block_resp_mat{2}(:,good_ind_green,ib,1,1),1),[],2)./sqrt(length(good_ind_green));
    errorbar(mean(mean(con_block_resp_mat{2}(:,good_ind_green,ib,1,1),1),2),mean(mean(con_block_resp_mat{1}(:,good_ind_green,ib,1,1),1),2),x_sem,x_sem,y_sem,y_sem,'ok');
    [h_green p_green] = ttest(mean(con_block_resp_mat{2}(:,good_ind_green,ib,1,1),1),mean(con_block_resp_mat{1}(:,good_ind_green,ib,1,1),1));
    xlabel('Mean resp D1')
    ylabel('Mean resp DN')
    xlim([-0.05 0.3])
    ylim([-0.05 0.3])
    axis square
    refline(1)
    title(['B' num2str(ib-1) ' HTP-; p = ' num2str(chop(p_green,2))])
    subplot(2,2,2+ib)
    scatter(mean(con_block_resp_mat{2}(:,good_ind_red,ib,1,1),1),mean(con_block_resp_mat{1}(:,good_ind_red,ib,1,1),1))
    hold on
    y_sem = std(mean(con_block_resp_mat{2}(:,good_ind_red,ib,1,1),1),[],2)./sqrt(length(good_ind_red));
    x_sem = std(mean(con_block_resp_mat{2}(:,good_ind_red,ib,1,1),1),[],2)./sqrt(length(good_ind_red));
    errorbar(mean(mean(con_block_resp_mat{2}(:,good_ind_red,ib,1,1),1),2),mean(mean(con_block_resp_mat{1}(:,good_ind_red,ib,1,1),1),2),x_sem,x_sem,y_sem,y_sem,'ok')
    [h_red p_red] = ttest(mean(con_block_resp_mat{2}(:,good_ind_red,ib,1,1),1),mean(con_block_resp_mat{1}(:,good_ind_red,ib,1,1),1));
    xlabel('Mean resp D1')
    ylabel('Mean resp DN')
    xlim([-0.02 0.15])
    ylim([-0.02 0.15])
    axis square
    refline(1)
    title(['B' num2str(ib-1) ' HTP+; p = ' num2str(chop(p_red,2))])
end
print(fullfile(fn_multi, ['contrastResp_byDay_scatter.pdf']),'-dpdf', '-bestfit')
save(fullfile(fn_multi, 'contrastResp.mat'),'h','p','data_dfof_con','con_resp_mat','data_dfof_con_block','con_block_resp_mat','trialn','good_ind','good_ind_red','good_ind_green','cons','tt','stimData')