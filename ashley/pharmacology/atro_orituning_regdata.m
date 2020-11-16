clear all
clear global
close all
ds = 'nbqx_oriTuning_V1'; %dataset info
rc = behavConstsAV; %directories
eval(ds)
slct_expt = 1; %which expt from ds to analyze
doPreviousReg = true;
doGreenOnly = false;
doDownsample = true;
%%
mouse = expt(slct_expt).mouse;
expDate = expt(slct_expt).date;
fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
fnout = fullfile(fn,'data processing');
mkdir(fnout)

%% downsample rate
if doDownsample
    dsFactor = 10;
else
    dsFactor = 1;
end

%% load and register gcamp data

% no stim runs
nruns = length(expt(slct_expt).nostim_runs);
data_nostim = cell(1,nruns);
for irun = 1:nruns
    
    runFolder = expt(slct_expt).nostim_runs{irun};
    mkdir(fullfile(fn,runFolder))
    fName = [runFolder '_000_000'];
        
    if expt(slct_expt).greenredsimultaneous == 1
        data_temp_g = loadsbx_choosepmt(3,mouse,expDate,runFolder,fName);
        data_temp_r = squeeze(data_temp_r(2,:,:,:));
        data_temp_g = squeeze(data_temp_g(1,:,:,:));
    else
        data_temp_g = loadsbx_choosepmt(1,mouse,expDate,runFolder,fName);
    end
    % remove negative data by subtraction
    data_sub = data_temp_g-min(min(min(data_temp_g,[],1),[],2),[],3);
    data_temp_g = data_sub;
    if expt(slct_expt).greenredsimultaneous == 1
        data_sub = data_temp_r-min(min(min(data_temp_r,[],1),[],2),[],3);
        data_temp_r = data_sub;
    end
    clear data_sub
    
    % downsample
    data_down_g = stackGroupProject(data_temp_g,dsFactor);
    clear data_temp_g
    
    % reg images
    nFr_regimg = round(99/dsFactor);
    if ~doPreviousReg && irun == 1
        regImgStartFrame = compareRegImg_2color(rc,ds,mouse,expDate,data_down_g);
        regImg = mean(data_down_g(:,:,regImgStartFrame:(regImgStartFrame+nFr_regimg)),3);
        save(fullfile(fnout,'regImg'),'regImg')
    elseif irun == 1
        load(fullfile(fnout,'regImg'))
    end
    
    if doPreviousReg
        load(fullfile(fn,runFolder,'regOuts&Img'))
        [~,data_g_reg] = stackRegister_MA(data_down_g,[],[],double(outs));
        figure;imagesc(regImg);colormap gray; title('reg image')
    else
        [outs,data_g_reg] = stackRegister(data_down_g,regImg);
        figure;imagesc(mean(data_g_reg,3));colormap gray; title(['no vis stim - ' expt(slct_expt).nostim_condition{irun}])
        clear data_temp_g
        save(fullfile(fn,runFolder,'regOuts&Img'),'outs','regImg','dsFactor')
    end
    
    data_nostim{irun} = data_g_reg;
    
    nfr = size(data_g_reg,3);
    nimg = 100;
    frInterval = 1:floor((nfr-100)/nimg):nfr-100;
    F = nan(size(data_g_reg,1),size(data_g_reg,2),nimg);
    for iimg = 1:nimg
        F(:,:,iimg) = mean(data_g_reg(:,:,frInterval(iimg):frInterval(iimg)+99),3);
    end
    writetiff(F,fullfile(fn,runFolder,'FImages'))
end
clear data_g_reg

% tuning runs
nruns = length(expt(slct_expt).dirtuning_runs);
data_tuning = cell(1,nruns);
data_mw = cell(1,nruns);
for irun = 1:nruns
    
    runFolder = expt(slct_expt).dirtuning_runs{irun};
    mkdir(fullfile(fn,runFolder))
    fName = [runFolder '_000_000'];
    nframes = expt(slct_expt).dirtuning_nframes{irun};
        
    if expt(slct_expt).greenredsimultaneous == 1
        data_temp_g = loadsbx_choosepmt(3,mouse,expDate,runFolder,fName,nframes);
        data_temp_r = squeeze(data_temp_r(2,:,:,:));
        data_temp_g = squeeze(data_temp_g(1,:,:,:));
    else
        data_temp_g = loadsbx_choosepmt(1,mouse,expDate,runFolder,fName,nframes);
    end
    % remove negative data by subtraction
    data_sub = data_temp_g-min(min(min(data_temp_g,[],1),[],2),[],3);
    data_temp_g = data_sub;
    if expt(slct_expt).greenredsimultaneous == 1
        data_sub = data_temp_r-min(min(min(data_temp_r,[],1),[],2),[],3);
        data_temp_r = data_sub;
    end
    clear data_sub
        
    % downsample
    data_down_g = stackGroupProject(data_temp_g,dsFactor);
    clear data_temp_g
    
    % reg images
    load(fullfile(fnout,'regImg'))
    
    if doPreviousReg
        load(fullfile(fn,runFolder,'regOuts&Img'))
        [~,data_g_reg] = stackRegister_MA(data_down_g,[],[],double(outs));
        figure;imagesc(regImg);colormap gray; title('reg image')
    else
        [outs,data_g_reg] = stackRegister(data_down_g,regImg);
        figure;imagesc(mean(data_g_reg,3));colormap gray; title(['no vis stim - ' expt(slct_expt).nostim_condition{irun}])
        clear data_temp_g
        save(fullfile(fn,runFolder,'regOuts&Img'),'outs','regImg','dsFactor')
    end
    
    data_tuning{irun} = data_g_reg;
    
    % save downsampled movie
    nfr = size(data_g_reg,3);
    nimg = 100;
    frInterval = 1:floor((nfr-100)/nimg):nfr-100;
    F = nan(size(data_g_reg,1),size(data_g_reg,2),nimg);
    for iimg = 1:nimg
        F(:,:,iimg) = mean(data_g_reg(:,:,frInterval(iimg):frInterval(iimg)+99),3);
    end
    writetiff(F,fullfile(fn,runFolder,'FImages'))
    
    % load mworks data
    mw_temp = loadMworksFile(mouse,expDate,expt(slct_expt).dirtuning_time{irun});
    data_mw{irun} = mw_temp;
end
clear data_g_reg mw_temp
close all

%% get red image
if ~isnan(expt(slct_expt).redChannelRun)
    redFolder = expt(slct_expt).redChannelRun;
    fName = [redFolder '_000_000'];
    data_r = loadsbx_choosepmt(2,mouse,expDate,redFolder,fName);
    load(fullfile(fnout,'regImg'))
    data_r_down = stackGroupProject(data_r,10);
    [~,data_rch_reg] = stackRegister(data_r_down,regImg);
    redChImg = mean(data_rch_reg,3);
    figure; colormap gray; imagesc(redChImg);
    save(fullfile(fnout,'redImage'),'redChImg')
    writetiff(redChImg,fullfile(fnout,'redImage'))
    clear data_r data_r down
end
%% get max dF/F for each gcamp run
[ypix,xpix,~] = size(data_nostim{1});
nfr = length(data_nostim)+length(data_tuning);
maxDFF_all = reshape(cell2mat(cellfun(@(x) quickMaxDFF(double(x)),cat(2,data_nostim,data_tuning),'unif',0)),[ypix,xpix,nfr]);

figure;colormap gray
[nrows,ncols] = optimizeSubplotDim(size(maxDFF_all,3));
for i = 1:size(maxDFF_all,3)
    subplot(nrows,ncols,i)
    imagesc(maxDFF_all(:,:,i));
end

writetiff(maxDFF_all,fullfile(fnout,'maxDFF_images'))
close all
%% cells selection image
if ~isnan(expt(slct_expt).redChannelRun)
    selection_img_temp = cat(3,redChImg,maxDFF_all);
else
    selection_img_temp = maxDFF_all;
end
[ypix,xpix,nfr] = size(selection_img_temp);

% smooth it
selection_img_smooth = medfilt3(selection_img_temp, [3 3 1]);

% crop it
figure;imagesc(selection_img_smooth(:,:,2));colormap gray

%**enter vals here***
xcrop = [1:10 780:xpix];
ycrop = [1:10 520:ypix];

selection_img = selection_img_smooth;
selection_img(:,xcrop) = 0;
selection_img(ycrop,:) = 0;

figure;imagesc(selection_img(:,:,2));colormap gray
% save it
save(fullfile(fnout,'selection_images.mat'),'selection_img','xcrop','ycrop');
clear selection_img_temp selection_img_smooth
%% get cells mask
if ~isnan(expt(slct_expt).redChannelRun)
    [red_mask_cell, green_mask_cell] = maskFromMultiMaxDFFStack_2color(selection_img);
else
    [green_mask_cell] = maskFromMultiMaxDFFStack(selection_img);
end

close all
figure; setFigParams4Print('portrait')
subplot 211
imagesc(green_mask_cell);
title({sprintf('%s %s- cells with behavior',...
    num2str(length(unique(green_mask_cell(:)))-1), ...
    expt(slct_expt).redChannelLabel);[mouse '-' expDate]})
if ~isnan(expt(slct_expt).redChannelRun)
subplot 212
imagesc(red_mask_cell);
title({sprintf('%s %s+ cells with behavior',...
    num2str(length(unique(red_mask_cell(:)))-1), ...
    expt(slct_expt).redChannelLabel);[mouse '-' expDate]})
end
print(fullfile(fnout,'final_mask_cells'),'-dpdf')

save(fullfile(fnout,'final_mask_cells.mat'),'green_mask_cell','red_mask_cell');
%% get time-courses
buf = 4;
np = 6;

tc_tagneg_nostim = cellfun(@(x) stackGetTimeCourses(x,green_mask_cell), ...
    data_nostim,'unif',0);
tc_tagneg_nostim_subnp = cellfun(@(x,y) getWeightedNeuropilTimeCourse(...
    x,y,green_mask_cell,buf,np), ...
    data_nostim, tc_tagneg_nostim,'unif',0);
tc_tagneg_tuning = cellfun(@(x) stackGetTimeCourses(x,green_mask_cell), ...
    data_tuning,'unif',0);
tc_tagneg_tuning_subnp = cellfun(@(x,y) getWeightedNeuropilTimeCourse(...
    x,y,green_mask_cell,buf,np), ...
    data_tuning, tc_tagneg_tuning,'unif',0);

if ~isnan(expt(slct_expt).redChannelRun)
tc_tagpos_nostim = cellfun(@(x) stackGetTimeCourses(x,red_mask_cell), ...
    data_nostim,'unif',0);
tc_tagpos_nostim_subnp = cellfun(@(x,y) getWeightedNeuropilTimeCourse(...
    x,y,red_mask_cell,buf,np), ...
    data_nostim, tc_tagpos_nostim,'unif',0);
tc_tagpos_tuning = cellfun(@(x) stackGetTimeCourses(x,red_mask_cell), ...
    data_tuning,'unif',0);
tc_tagpos_tuning_subnp = cellfun(@(x,y) getWeightedNeuropilTimeCourse(...
    x,y,red_mask_cell,buf,np), ...
    data_tuning, tc_tagpos_tuning,'unif',0);
end

% save time-courses
if ~isnan(expt(slct_expt).redChannelRun)
save(fullfile(fnout,'timecourses_cells'),'tc_tagneg_nostim_subnp',...
    'tc_tagpos_nostim_subnp','tc_tagneg_tuning_subnp','tc_tagpos_tuning_subnp')
else
save(fullfile(fnout,'timecourses_cells'),'tc_tagneg_nostim_subnp',...
    'tc_tagneg_tuning_subnp')
end