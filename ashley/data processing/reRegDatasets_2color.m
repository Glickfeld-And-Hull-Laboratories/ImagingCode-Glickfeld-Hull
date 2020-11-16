clear 
clear global
close all
ds = 'FSV_V1'; %dataset info
rc = behavConstsAV; %directories
eval(ds)
slct_expt = 30; %which expt from ds to analyze
doPreviousReg = true;
doGreenOnly = false;
nRedFrames = 20000;
% regPath = 'Z:\home\former-lab-members\carolyn\Analysis\1207\two-photon imaging\181018\data processing\regOuts&Img.mat';
%%

SubNum = expt(slct_expt).SubNum;
mouse = expt(slct_expt).mouse;
expDate = expt(slct_expt).date;
dirFolder = expt(slct_expt).dirtuning;
dirTime = expt(slct_expt).dirtuning_time;
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
nFrPerRun = zeros(1,expt(slct_expt).nrun);
for irun = 1:expt(slct_expt).nrun

    runFolder = expt(slct_expt).runs(irun,:);
    expTime = expt(slct_expt).time_mat(irun,:);
    fName = [runFolder '_000_000'];
    
    if ~isempty(expt(slct_expt).nframesPerRun)
        nframes = expt(slct_expt).nframesPerRun(irun);
        data_temp_g = loadsbx_choosepmt(1,mouse,expDate,runFolder,fName,expt(slct_expt).nframesPerRun);
        if expt(slct_expt).greenredsimultaneous == 1
            data_temp_r = loadsbx_choosepmt(2,mouse,expDate,runFolder,fName,expt(slct_expt).nframesPerRun);
    %         [input, data_temp, t] = Load_SBXdataPlusMWorksData(...
    %             SubNum,expDate,expTime,mouse,runFolder,fName,nframes);
        end
    else
        if expt(slct_expt).greenredsimultaneous == 1
            data_temp_g = loadsbx_choosepmt(3,mouse,expDate,runFolder,fName);
            data_temp_r = squeeze(data_temp_r(2,:,:,:));
            data_temp_g = squeeze(data_temp_g(1,:,:,:));
        else
            data_temp_g = loadsbx_choosepmt(1,mouse,expDate,runFolder,fName);
        end
%         [input, data_temp, t] = Load_SBXdataPlusMWorksData(...
%             SubNum,expDate,expTime,mouse,runFolder,fName);
    end
    input = loadMworksFile(SubNum,expDate,expTime);
    
    nFrPerRun(irun) = size(data_temp_g,3);
    if irun == 1
        data_bx_g = data_temp_g;
        if expt(slct_expt).greenredsimultaneous == 1
            data_bx_r = data_temp_r;
        end
        input_bx = input;
    else
        data_bx_g = cat(3,data_bx_g,data_temp_g);
        if expt(slct_expt).greenredsimultaneous == 1
            data_bx_r = cat(3,data_bx_r,data_temp_r);
        end
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
    clear data_temp_g data_temp_r input
end
input_bx = concatenateDataBlocks(input_bx);
nfr_bx = size(data_bx_g,3);

% remove negative data by subtraction
data_sub = data_bx_g-min(min(min(data_bx_g,[],1),[],2),[],3);
data_bx_g = data_sub;
if expt(slct_expt).greenredsimultaneous == 1
    data_sub = data_bx_r-min(min(min(data_bx_r,[],1),[],2),[],3);
    data_bx_r = data_sub;
end
clear data_sub

% select registration image
if ~doPreviousReg
    if expt(slct_expt).greenredsimultaneous == 1
        regImgStartFrame = compareRegImg_2color(rc,ds,SubNum,expDate,data_bx_g,data_bx_r);
    else
        regImgStartFrame = compareRegImg_2color(rc,ds,SubNum,expDate,data_bx_g);
    end
end
%% register
% register 
if doPreviousReg
    load([fnout '\regOuts&Img'])
    [~,data_bx_g_reg] = stackRegister_MA(data_bx_g,[],[],double(out_bx));
    if ~doGreenOnly
        if expt(slct_expt).greenredsimultaneous == 1
            data_bx_r = data_bx_r(:,:,1:nRedFrames);
            [~,data_bx_r_reg] = stackRegister_MA(data_bx_r,[],[],double(out_bx(1:nRedFrames,:)));   
        end
    end
    clear data_bx_g
elseif expt(slct_expt).greenredsimultaneous == 0
    regImg = mean(data_bx_g(:,:,regImgStartFrame:(regImgStartFrame+99)),3);
    figure;imagesc(regImg);colormap gray
    [out_bx,data_bx_g_reg] = stackRegister(data_bx_g,regImg);
    clear data_bx_g
    save([fnout '\regOuts&Img'],'out_bx','regImg')
else
    data_bx_r = data_bx_r(:,:,1:nRedFrames);
    regImg = mean(data_bx_g(:,:,regImgStartFrame:(regImgStartFrame+99)),3);
    figure;imagesc(regImg);colormap gray
    [out_bx,data_bx_g_reg] = stackRegister(data_bx_g,regImg);
    [~,data_bx_r_reg] = stackRegister_MA(data_bx_r,[],[],out_bx(1:nRedFrames,:));
%     [out_tun,data_tun_reg] = stackRegister(data_tun,regImg);
    clear data_bx_g data_bx_r
    save([fnout '\regOuts&Img'],'out_bx','regImg')
end

% save some F images as tifs
nImg = 200;
frInterval = 1:floor((nfr_bx-100)/nImg):nfr_bx-100;
F_bx = nan(size(data_bx_g_reg,1),size(data_bx_g_reg,2),nImg);
for iimg = 1:nImg
    F_bx(:,:,iimg) = mean(data_bx_g_reg(:,:,frInterval(iimg):frInterval(iimg)+99),3);
end
writetiff(F_bx,fullfile(fnout,'FImages_bx'))
% F_tun = mean(data_tun_reg,3);

%% load red image
if doGreenOnly
    redImgSelect = zeros(size(regImg));
else
    if ~isnan(expt(slct_expt).redChannelRun)
        redFolder = expt(slct_expt).redChannelRun;
        fName = [redFolder '_000_000'];
        data_r = loadsbx_choosepmt(2,mouse,expDate,redFolder,fName);
        [~,data_rch_reg] = stackRegister(data_r,regImg);
        redChImg = mean(data_rch_reg,3);
        figure; colormap gray; imagesc(redChImg);
    elseif expt(slct_expt).greenredsimultaneous
        redImgSelect = mean(data_bx_r_reg,3);
        [redChImg,s_star,s_regress] = greenBleedIntoRed(redImgSelect,data_bx_r_reg,data_bx_g_reg(:,:,1:nRedFrames));
        save([fnout 'greenChBleedthrough'],'redChImg','s_star','s_regress')
        clear data_bx_r_reg
    end
    writetiff(redChImg,fullfile(fnout,'redChannelRun'))
end
%% get task max dF/F for 1st resp, mid resp, and target resp
data_bx_reg = data_bx_g_reg;
getTaskMaxDFF

close all
%% crop cells image

%crop it
% tun_img = max(dFF_dirmax,[],3);
bx_img = max(dFF_bxMax,[],3);
% tuning image
% figure;colormap gray; imagesc(tun_img)

figure;colormap gray; imagesc(bx_img)
%**enter vals here***
xcrop = [1:20 780:xpix];
ycrop = [1:10 250:ypix];

% tun_crop = tun_img;
% tun_crop(:,xcrop) = 0;
% tun_crop(ycrop,:) = 0;
% 
% imagesc(tun_crop)

% check if behavior image still needs cropping
bx_crop = bx_img;
bx_crop(:,xcrop) = 0;
bx_crop(ycrop,:) = 0;

imagesc(bx_crop)

% crop orig max images
% dir_crop = dFF_dirmax;
% dir_crop(:,xcrop,:) = 0;
% dir_crop(ycrop,:,:) = 0;

bx_crop = cat(3, start_max, long_max, tar_max);
bx_crop(:,xcrop,:) = 0;
bx_crop(ycrop,:,:) = 0;

save(fullfile(fnout,'max_images_crop.mat'),'bx_crop','xcrop','ycrop');
close all

%% get cells mask
if doGreenOnly
    green_mask_cell = maskFromMultiMaxDFFStack(bx_crop);
    close all
    
    figure; setFigParams4Print('portrait')
    imagesc(green_mask_cell);
    title({sprintf('%s %s cells with behavior',...
        num2str(length(unique(green_mask_cell(:)))-1), ...
        expt(slct_expt).greenChannelLabel);[mouse '-' expDate]})
    print(fullfile(fnout,'final_mask_cells'),'-dpdf')
    save(fullfile(fnout,'final_mask_cells.mat'),'green_mask_cell');
else
    
    dFF_stack = cat(3,redChImg,bx_crop);
    [red_mask_cell, green_mask_cell] = maskFromMultiMaxDFFStack_2color(dFF_stack);
    close all
    figure; setFigParams4Print('portrait')
    subplot 211
    imagesc(green_mask_cell);
    title({sprintf('%s %s+ cells with behavior',...
        num2str(length(unique(green_mask_cell(:)))-1), ...
        expt(slct_expt).redChannelLabel);[mouse '-' expDate]})
    subplot 212
    imagesc(red_mask_cell);
    title({sprintf('%s %s- cells with behavior',...
        num2str(length(unique(red_mask_cell(:)))-1), ...
        expt(slct_expt).redChannelLabel);[mouse '-' expDate]})
    print(fullfile(fnout,'final_mask_cells'),'-dpdf')

    save(fullfile(fnout,'final_mask_cells.mat'),'green_mask_cell','red_mask_cell');
end

%% get time-courses
buf = 4;
np = 6;
if ~doGreenOnly
    data_bx_r_tc = stackGetTimeCourses(data_bx_g_reg,red_mask_cell);
    data_bx_r_tc_subnp = getWeightedNeuropilTimeCourse(data_bx_g_reg,data_bx_r_tc,red_mask_cell,buf,np);
end

data_bx_g_tc = stackGetTimeCourses(data_bx_g_reg,green_mask_cell);
data_bx_g_tc_subnp = getWeightedNeuropilTimeCourse(data_bx_g_reg,data_bx_g_tc,green_mask_cell,buf,np);

% if ~exist(fullfile(fn,dirFolder),'dir')
%     mkdir(fntun,dirFolder)
% end
if doGreenOnly    
    save(fullfile(fnout,'timecourses_bx_cells.mat'),...
        'data_bx_g_tc_subnp','data_bx_g_tc','buf','np')
else
    save(fullfile(fnout,'timecourses_bx_cells.mat'),...
        'data_bx_g_tc_subnp','data_bx_g_tc',...
        'data_bx_r_tc_subnp','data_bx_r_tc','buf','np')
end
clear data_bx_g_reg data_bx_r_reg
%% register passive expt, get TCs
load(fullfile(fnout,'final_mask_cells.mat'),'green_mask_cell','red_mask_cell');
if isfield(expt,'passExpt') 
    if ~isnan(expt(slct_expt).passExpt)
        fName = [expt(slct_expt).passExpt '_000_000'];
        data_pass_g = loadsbx_choosepmt(1,mouse,expDate,expt(slct_expt).passExpt,fName);
%         data_pass_r = loadsbx_choosepmt(2,mouse,expDate,expt(slct_expt).passExpt,fName);
        nfr_pass = size(data_pass_g,3);

%         [out_pass,data_pass_r_reg] = stackRegister(data_pass_r,regImg);
        [out_pass,data_pass_g_reg] = stackRegister(data_pass_g,regImg);
        save(fullfile(fnout,'regOuts&Img_pass.mat'),'out_pass','regImg')
        clear data_pass_g data_pass_r

        frInterval_pass = 1:floor((nfr_pass-100)/nImg):nfr_pass-100;
        F_pass = nan(size(data_pass_g_reg,1),size(data_pass_g_reg,2),nImg);
        for iimg = 1:nImg
            F_pass(:,:,iimg) = mean(data_pass_g_reg(:,:,frInterval_pass(iimg):frInterval_pass(iimg)+99),3);
        end
        writetiff(F_pass,fullfile(fnout,'FImages_pass'))

        if doGreenOnly
            data_pass_g_tc = stackGetTimeCourses(data_pass_g_reg,green_mask_cell);
            data_pass_g_tc_subnp = getWeightedNeuropilTimeCourse(data_pass_g_reg,data_pass_g_tc,green_mask_cell,buf,np);
        else
            data_pass_r_tc = stackGetTimeCourses(data_pass_g_reg,red_mask_cell);
            data_pass_r_tc_subnp = getWeightedNeuropilTimeCourse(data_pass_g_reg,data_pass_r_tc,red_mask_cell,buf,np);

            data_pass_g_tc = stackGetTimeCourses(data_pass_g_reg,green_mask_cell);
            data_pass_g_tc_subnp = getWeightedNeuropilTimeCourse(data_pass_g_reg,data_pass_g_tc,green_mask_cell,buf,np);
        end

        save(fullfile(fnout,'timecourses_pass_cells.mat'),...
            'data_pass_g_tc_subnp','data_pass_g_tc',...
            'data_pass_r_tc_subnp','data_pass_r_tc','buf','np')
        
%         save(fullfile(fnout,'timecourses_pass_cells.mat'),...
%             'data_pass_g_tc_subnp','data_pass_g_tc','buf','np')
    end
end