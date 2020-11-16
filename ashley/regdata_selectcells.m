clear all
clear global
close all
ds = 'DART_V1_PV_contrast'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsAV; %directories
eval(ds)
slct_expt = 9; %which expt from ds to analyze
doPreviousReg = false;
doGreenOnly = false;
doCorrImg = true;
%%
mouse = expt(slct_expt).mouse;
expDate = expt(slct_expt).date;
fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
fnout = fullfile(fn,'data processing');
mkdir(fnout)
%% load and register gcamp data
data_label = [];
for idatalabel = 1:length(dataStructLabels)
    runs = eval(['expt(slct_expt).' dataStructLabels{idatalabel} '_runs']);
    times = eval(['expt(slct_expt).' dataStructLabels{idatalabel} '_time']);
    nruns = length(runs);
    data = cell(1,nruns);
    for irun = 1:nruns
        runFolder = runs{irun};
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

        % reg images
        if ~doPreviousReg && irun == 1 && idatalabel == 1
            regImgStartFrame = compareRegImg_2color(rc,ds,mouse,expDate,data_temp_g);
            regImg = mean(data_temp_g(:,:,regImgStartFrame:(regImgStartFrame+99)),3);
            save(fullfile(fnout,'regImg'),'regImg')
        elseif irun == 1
            load(fullfile(fnout,'regImg'))
        end

        if doPreviousReg
            load(fullfile(fn,runFolder,'regOuts&Img'))
            [~,data_g_reg] = stackRegister_MA(data_temp_g,[],[],double(outs));
            figure;imagesc(regImg);colormap gray; title('reg image')
        else
            [outs,data_g_reg] = stackRegister(data_temp_g,regImg);
            figure;imagesc(mean(data_g_reg,3));colormap gray
            clear data_temp_g
            save(fullfile(fn,runFolder,'regOuts&Img'),'outs','regImg')
        end

        data{irun} = data_g_reg;

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


    %%
    data_mw = cell(1,nruns);
    for irun = 1:nruns
        data_mw{irun} = loadMworksFile(mouse,expDate,times{irun});
    end
    %%
    if strcmp(dataStructLabels{idatalabel},'contrast')|strcmp(dataStructLabels{idatalabel},'contrastxsize')|strcmp(dataStructLabels{idatalabel},'contrastxori')
        nBaseFr = round(expt(slct_expt).frame_rate);
        nStimFr = round(expt(slct_expt).frame_rate);
        sz = cellfun(@size,data,'unif',0);
        off = cellfun(@(x) x.nScansOff,data_mw,'unif',0);
        on = cellfun(@(x) x.nScansOn,data_mw,'unif',0);
        nt = cellfun(@(x,y,z) floor(x(3)./(y+z)),sz,off,on,'unif',0);

        maxConInd = cellfun(@(x,y) celleqel2mat_padded(x.tGratingContrast) == ...
            max(celleqel2mat_padded(x.tGratingContrast)),data_mw,'unif',0);
        trialStartFr = cellfun(@(x,y,z) (x+1):(x+y):(z*(x+y)),...
            off,on,nt,'unif',0);
        trialStartFr = cellfun(@(x,y) x(y),trialStartFr,maxConInd,'unif',0);
        stimID = cellfun(@(x,y) cell2mat_padded(x.tGratingContrast(y)),...
            data_mw,maxConInd,'unif',0);
        ind = cellfun(@(fr,d) (fr+nStimFr) < size(d,3),trialStartFr,data,'unif',0);

        maxDFF_stim = reshape(cell2mat(cellfun(@(d,fr,id,n,i) maxDFF_fromStim(d,fr(i),id(i),nBaseFr,nStimFr),...
            data,trialStartFr,stimID,nt,ind,'unif',0)),[sz{1}(1),sz{1}(2),length(sz)]);
        
        
        figure;colormap gray
        [nrows,ncols] = optimizeSubplotDim(size(maxDFF_stim,3));
        for i = 1:size(maxDFF_stim,3)
            subplot(nrows,ncols,i)
            imagesc(maxDFF_stim(:,:,i));
        end

        writetiff(maxDFF_stim,fullfile(fnout,'maxDFF_images'))
        selection_images{idatalabel} = maxDFF_stim;
    else
        selection_images{idatalabel} = zeros(size(data{1},1),size(data{1},2));
    end
    data_label{idatalabel} = data;
    clear data;
end
%% get red image
if ~doGreenOnly
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

%% 
if doCorrImg
    corrImg = getPixelCorrelationImage(double(data_label{1}{1}(:,:,1:5000)));
    figure; imagesc(corrImg); colormap gray
    print(fullfile(fnout,'corrMap'),'-dpdf')
    selection_images{end+1} = corrImg;
end
%% cells selection image
[ypix,xpix,~] = size(selection_images{1});
nfr = sum(cellfun(@(x) size(x,3),selection_images));

if doGreenOnly
    selection_img_temp = reshape(cell2mat(selection_images),[ypix,xpix,nfr]);
else
    selection_img_temp = redChImg;
    for i = 1:length(selection_images)
        selection_img_temp = cat(3,selection_img_temp,...
        selection_images{i});
    end
end

% smooth it
selection_img_smooth = medfilt3(selection_img_temp, [3 3 1]);

% crop it
figure;imagesc(selection_img_smooth(:,:,1));colormap gray

%**enter vals here***
xcrop = [1:10 780:xpix];
ycrop = [1:15 520:ypix];

selection_img = selection_img_smooth;
selection_img(:,xcrop) = 0;
selection_img(ycrop,:) = 0;

figure;imagesc(selection_img(:,:,1));colormap gray
% save it
save(fullfile(fnout,'selection_images.mat'),'selection_img','xcrop','ycrop');
clear selection_img_temp selection_img_smooth

%% get cells mask
if doGreenOnly
    green_mask_cell = maskFromMultiMaxDFFStack(selection_img);
else
    [red_mask_cell, green_mask_cell] = maskFromMultiMaxDFFStack_2color(selection_img); % first image should be red
end
close all
figure; setFigParams4Print('portrait')
subplot 211
if ~doGreenOnly
    imagesc(red_mask_cell);
    title({sprintf('%s %s+ cells with behavior',...
        num2str(length(unique(red_mask_cell(:)))-1), ...
        expt(slct_expt).redChannelLabel);[mouse '-' expDate]})
end
subplot 212
imagesc(green_mask_cell);
title({sprintf('%s %s- cells with behavior',...
    num2str(length(unique(green_mask_cell(:)))-1), ...
    expt(slct_expt).redChannelLabel);[mouse '-' expDate]})
print(fullfile(fnout,'final_mask_cells'),'-dpdf')
if doGreenOnly
    save(fullfile(fnout,'final_mask_cells.mat'),'green_mask_cell');
else
    save(fullfile(fnout,'final_mask_cells.mat'),'green_mask_cell','red_mask_cell');
end
%% get time-courses
buf = 4;
np = 6;
    
for idatalabel = 1:length(dataStructLabels)
    runs = eval(['expt(slct_expt).' dataStructLabels{idatalabel} '_runs']);
    times = eval(['expt(slct_expt).' dataStructLabels{idatalabel} '_time']);
    nruns = length(runs);
    for irun = 1:nruns
        runFolder = runs{irun};
        d_temp = data_label{idatalabel}{irun};
        tc_green = stackGetTimeCourses(d_temp,green_mask_cell);
        tc_subnp_green = getWeightedNeuropilTimeCourse(d_temp,tc_green,...
            green_mask_cell,buf,np);
        if doGreenOnly
            save(fullfile(fn,runFolder,'timecourses_cells'),'tc_green','tc_subnp_green')   
        else
            tc_red = stackGetTimeCourses(d_temp,red_mask_cell);
            tc_subnp_red = getWeightedNeuropilTimeCourse(d_temp,tc_red,...
                red_mask_cell,buf,np);
            save(fullfile(fn,runFolder,'timecourses_cells'),...
                'tc_green','tc_subnp_green','tc_red','tc_subnp_red') 
        end
    end
end
