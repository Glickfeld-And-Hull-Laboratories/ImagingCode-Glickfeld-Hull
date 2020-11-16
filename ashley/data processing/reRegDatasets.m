clear all
close all
ds = 'FSAV_V1_naive_GCaMP6m';
rc = behavConstsAV;
eval(ds)
slct_expt = 9;
doPreviousReg = true;
%%
iexp = slct_expt

SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
dirFolder = expt(iexp).dirtuning;
dirTime = expt(iexp).dirtuning_time;
if strcmp(rc.name, 'ashle')
    fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
    if ~exist(fullfile(fn,'data processing'),'dir')
        mkdir(fn,'data processing')
    end
elseif strcmp(rc.name, 'carolyn')
    fn = fullfile(rc.carolynAnalysis,mouse,'Two-Photon Imaging',expDate);
    if ~exist(fullfile(fn,'data processing'),'dir')
        mkdir(fn,'data processing')
    end  
else
    error('you do not belong here')
end
fnout = fullfile(fn,'data processing');
%% load and register all behavior and tuning data
% concatenate data files
nFrPerRun = zeros(1,expt(iexp).nrun);
for irun = 1:expt(iexp).nrun

    runFolder = expt(iexp).runs(irun,:);
    expTime = expt(iexp).time_mat(irun,:);
    fName = [runFolder '_000_000'];
    
    if (strcmp(ds, '_V1') | strcmp(ds(1:4),'FSAV')) & irun == 3 & strcmp(expDate,'150508')
        nframes = 13577;
        [input, data_temp, t] = Load_SBXdataPlusMWorksData(SubNum,expDate,expTime,mouse,runFolder,fName,nframes);
%     elseif (strcmp(ds, '_V1') | strcmp(ds,'')) & irun == 2 & strcmp(expDate, '150508')
%         continue
    elseif isfield(expt,'nframesPerRun')
        if ~isempty(expt(iexp).nframesPerRun)
            nframes = expt(iexp).nframesPerRun(irun);
            [input, data_temp, t] = Load_SBXdataPlusMWorksData(...
                SubNum,expDate,expTime,mouse,runFolder,fName,nframes);
        else
            [input, data_temp, t] = Load_SBXdataPlusMWorksData(...
                SubNum,expDate,expTime,mouse,runFolder,fName);
        end
    else
        [input, data_temp, t] = Load_SBXdataPlusMWorksData(...
            SubNum,expDate,expTime,mouse,runFolder,fName);
    end
    
    disp(t)
    nFrPerRun(irun) = size(data_temp,3);
    if irun == 1
        data_bx = data_temp;
        input_bx = input;
    else
        data_bx = cat(3,data_bx,data_temp);
        try
            input_bx = [input_bx input];
        catch
            inpNames1 = fieldnames(input_bx);
            inpNames2 = fieldnames(input);
            inpLong = gt(length(inpNames1),length(inpNames2));
            if inpLong == 1
                inpPlusInd = ismember(inpNames1,inpNames2);
                inpPlus = inpNames1(~inpPlusInd);
                for i = 1:length(inpPlus)
                    input.(genvarname(inpPlus{i})) = cell(1,input.trialSinceReset);
                end
            else
                inpPlusInd = ismember(inpNames2,inpNames1);
                inpPlus = inpNames2(~inpPlusInd);
                for i = 1:length(inpPlus)
                    input_bx.(char(genvarname(inpPlus(i)))) = cell(1,length(inpNames1));
                end
            end
            input_bx = [input_bx input];
        end
    end
    clear data_temp input
end
input_bx = concatenateDataBlocks(input_bx);
nfr_bx = size(data_bx,3);

% remove negative data by subtraction
data_sub = data_bx-min(min(min(data_bx,[],1),[],2),[],3);
data_bx = data_sub;
clear data_sub

% load tuning data
down = 10;
if strcmp(rc.name, 'ashle')
    fntun = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate,dirFolder);
elseif strcmp(rc.name, 'ashle')
    fntun = fullfile(rc.carolynAnalysis,mouse,'two-photon imaging',expDate,dirFolder);
end

fName = [dirFolder '_000_000'];

if any(strcmp(fieldnames(expt),'nTunFrames'))
    [input_tun, data_tun] = Load_SBXdataPlusMWorksData(SubNum,expDate,dirTime,mouse,dirFolder,fName,expt(iexp).nTunFrames);  
else
    [input_tun, data_tun] = Load_SBXdataPlusMWorksData(SubNum,expDate,dirTime,mouse,dirFolder,fName);
end

% down-sample
data_tun_down = stackGroupProject(data_tun,down);
clear data_tun
nfr_tun = size(data_tun_down,3);

% remove negative data by subtraction
data_sub = data_tun_down-min(min(min(data_tun_down,[],1),[],2),[],3);
data_tun = data_sub;
clear data_sub

if ~doPreviousReg
    regImgStartFrame = compareRegImg_2color(rc,ds,SubNum,expDate,data_bx);
else
    regImgStartFrame = expt(iexp).regImgStartFrame;
end

% register 
if doPreviousReg
    load(fullfile(fnout,'regOuts&Img.mat'))
    out_bx = double(out_bx);
    out_tun = double(out_tun);
    [~,data_bx_reg] = stackRegister_MA(data_bx,[],[],out_bx);
    [~,data_tun_reg] = stackRegister_MA(data_tun,[],[],out_tun);
    clear data_bx data_tun
else
    regImg = mean(data_bx(:,:,regImgStartFrame:(regImgStartFrame+99)),3);
    figure;imagesc(regImg);colormap gray
    [out_bx,data_bx_reg] = stackRegister(data_bx,regImg);
    [out_tun,data_tun_reg] = stackRegister(data_tun,regImg);
    clear data_bx
    save(fullfile(fnout,'regOuts&Img.mat'),'out_bx','out_tun', 'regImg')
end

% save some F images as tifs
nImg = 200;
frInterval = 1:floor((nfr_bx-100)/nImg):nfr_bx-100;
F_bx = nan(size(data_bx_reg,1),size(data_bx_reg,2),nImg);
for iimg = 1:nImg
    F_bx(:,:,iimg) = mean(data_bx_reg(:,:,frInterval(iimg):frInterval(iimg)+99),3);
end
writetiff(F_bx,fullfile(fnout,'FImages_bx'))
F_tun = mean(data_tun_reg,3);
writetiff(F_tun,fullfile(fnout,'FImages_tun'))
%% get task max dF/F for 1st resp, mid resp, and target resp

getTaskMaxDFF

dirTuningMaxDFF

close all
%% crop cells image

%crop it
tun_img = max(dFF_dirmax,[],3);
bx_img = max(dFF_bxMax,[],3);
% tuning image
figure;colormap gray; imagesc(tun_img)

%**enter vals here***
xcrop = [1:31 790:xpix];
ycrop = [1:15 260:ypix];

tun_crop = tun_img;
tun_crop(:,xcrop) = 0;
tun_crop(ycrop,:) = 0;

imagesc(tun_crop)

% check if behavior image still needs cropping
figure;colormap gray; imagesc(bx_img)
bx_crop = bx_img;
bx_crop(:,xcrop) = 0;
bx_crop(ycrop,:) = 0;

imagesc(bx_crop)

% crop orig max images
dir_crop = dFF_dirmax;
dir_crop(:,xcrop,:) = 0;
dir_crop(ycrop,:,:) = 0;

bx_crop = cat(3, start_max, long_max, tar_max);
bx_crop(:,xcrop,:) = 0;
bx_crop(ycrop,:,:) = 0;

save(fullfile(fnout,'max_images_crop.mat'),'dir_crop', 'bx_crop','xcrop','ycrop');
close all

%% get cells mask
dFF_stack = cat(3,dir_crop,bx_crop);
mask_cell = maskFromMultiMaxDFFStack(dFF_stack);

figure; setFigParams4Print('portrait')
imagesc(mask_cell);
title({[num2str(length(unique(mask_cell(:)))-1) ' cells with behavior'];[mouse '-' expDate]})
print(fullfile(fnout,'final_mask_cells'),'-dpdf')

save(fullfile(fnout,'final_mask_cells.mat'),'mask_cell');
close all
%% get dendrites mask
dFF_2D = reshape(dFF_stack,xpix*ypix,[]);
dFF_2D(mask_cell(:) > 0,:) = 0;

dFF_stack_den = reshape(dFF_2D,ypix,xpix,[]);
mask_den = maskFromMultiMaxDFFStack(dFF_stack_den);
figure; setFigParams4Print('portrait')
imagesc(mask_den);
title({[num2str(length(unique(mask_den(:)))-1) ' dendrites with behavior'];[mouse '-' expDate]})
print(fullfile(fnout,'final_mask_dendrites'),'-dpdf')

save(fullfile(fnout,'final_mask_dendrites.mat'),'mask_den');
close all
%% get time-courses
buf = 4;
np = 6;

data_bx_tc = stackGetTimeCourses(data_bx_reg,mask_cell);
data_bx_tc_subnp = getWeightedNeuropilTimeCourse(data_bx_reg,data_bx_tc,mask_cell,buf,np);

data_tun_tc = stackGetTimeCourses(data_tun_reg,mask_cell);
data_tun_tc_subnp = getWeightedNeuropilTimeCourse(data_tun_reg,data_tun_tc,mask_cell,buf,np);

data_bx_den_tc = stackGetTimeCourses(data_bx_reg,mask_den);
data_bx_den_tc_subnp = getWeightedNeuropilTimeCourse(data_bx_reg,data_bx_den_tc,mask_den,buf,np);

data_tun_den_tc = stackGetTimeCourses(data_tun_reg,mask_den);
data_tun_den_tc_subnp = getWeightedNeuropilTimeCourse(data_tun_reg,data_tun_den_tc,mask_den,buf,np);

if ~exist(fullfile(fn,dirFolder),'dir')
    mkdir(fntun,dirFolder)
end
save(fullfile(fnout,'timecourses_bx_cells.mat'),'data_bx_tc_subnp','data_bx_tc','buf','np')
save(fullfile(fntun,'timecourses_tun_cells.mat'),'data_tun_tc_subnp','data_tun_tc')

save(fullfile(fnout,'timecourses_bx_dendrites.mat'),'data_bx_den_tc_subnp','data_bx_den_tc','buf','np')
save(fullfile(fntun,'timecourses_tun_dendrites.mat'),'data_tun_den_tc_subnp','data_tun_den_tc')
