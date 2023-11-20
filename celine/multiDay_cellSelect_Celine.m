clear all; clear global;  close all
clc
%ds = 'con_ori_nonDART'
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories
eval(ds);
doGreenOnly = true;
doCorrImg = true;

day_id = 355;
%% load data for day

mouse = expt(day_id).mouse;
expDate = expt(day_id).date;

fn = fullfile(rc.achAnalysis,mouse,expDate); %can make this flexible if folder structure is different
mkdir(fn)

runs = eval(['expt(day_id).' cell2mat(dataStructLabels) '_runs']);
times = eval(['expt(day_id).' cell2mat(dataStructLabels) '_time']);
nruns = length(runs);
runFolder = [];
for irun = 1:nruns 
    imgFolder = runs{irun};
    if nruns == 1
        runFolder = imgFolder;
        fnout = fullfile(fn,runFolder);
        mkdir(fullfile(fn,runFolder))
        fName = [imgFolder '_000_000'];
    elseif irun < nruns
        runFolder = [runFolder '_' imgFolder];
        fName = [imgFolder '_000_000'];
        else1
        runFolder = [runFolder '_' imgFolder];
        fnout = fullfile(fn,runFolder);
        mkdir(fullfile(fn,runFolder))
        fName = [imgFolder '_000_000'];
    end


    if strcmp(expt(day_id).data_loc,'lindsey')
        root = rc.data;
        CD = fullfile(root,mouse,expDate,runFolder);
        dat = 'data-';
    elseif strcmp(expt(day_id).data_loc,'ashley')
        root = rc.ashleyData;
        CD = fullfile(root,mouse,'two-photon imaging', expDate,runFolder);
        dat = 'data-i';
    elseif strcmp(expt(day_id).data_loc,'tammy')
        root = rc.tammyData;
        CD = fullfile(root, mouse, '2P',expDate, runFolder);
        dat = 'data-i';
    elseif strcmp(expt(day_id).data_loc,'celine')
        root = rc.Data;
        CD = fullfile(root, mouse, expDate, runFolder);
        dat = 'data-i';
    elseif strcmp(expt(day_id).data_loc,'ACh')
        root = rc.achData;
        CD = fullfile(root, mouse, expDate, runFolder);
        dat = 'data-';
    end
    cd(CD);

    imgMatFile = [imgFolder '_000_000.mat'];
    load(imgMatFile);
    tHostname = lower(hostname);
    [s,tUsername] = dos('ECHO %USERNAME%');
    switch tHostname
        case {'nuke'}
            if username == 'celine' 
                fName = ['Z:\Behavior\Data\' dat mouse '-' expDate '-' times{irun} '.mat'];

            else
                fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\' dat mouse '-' expDate '-' times{irun} '.mat'];
            end
        case{'nb-hubel'}
                if username == 'cc735'
                    fName = ['Z:\Behavior\Data\' dat mouse '-' expDate '-' times{irun} '.mat'];
                else
                    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\' dat mouse '-' expDate '-' times{irun} '.mat'];
                 end
    end
    
    load(fName); %load the mworks behavioral file

    temp(irun) = input; %load the data from the mworks file into temp
    clear input
    nframes = [temp(irun).counterValues{end}(end) info.config.frames];

    fprintf(['Reading run ' num2str(irun) '- ' num2str(min(nframes)) ' frames \r\n'])

    if info.config.pmt1_gain > 0.5
        data_temp_g = sbxread(imgMatFile(1,1:11),0,min(nframes));
        data_temp_r = squeeze(data_temp_g(2,:,:,:));
        data_temp_g = squeeze(data_temp_g(1,:,:,:));
    else
        data_temp_g = sbxread(imgMatFile(1,1:11),0,min(nframes));
        data_temp_g = squeeze(data_temp_g(1,:,:,:));
    end

    fprintf(['Loaded ' num2str(min(nframes)) ' frames \r\n'])

    if nruns == 1 || irun == 1
        data_g = data_temp_g;
        clear data_temp_g
        if info.config.pmt1_gain > 0.5
            data_r = data_temp_r;
            clear data_temp_r
        end
    else
        data_g = cat(3, data_g, data_temp_g);
        clear data_temp_g
        if info.config.pmt1_gain > 0.5
            data_r = cat(3, data_r, data_temp_r);
            clear data_temp_r
        end
    end
end
%%
% register data for each day
%reg green data

if exist(fullfile(fnout,'regOuts&Img.mat')) %check if there is already registration info for this data
    load(fullfile(fnout,'regOuts&Img.mat'))
    [~,data_g_reg] = stackRegister_MA(data_g,[],[],double(outs));
    data_avg = mean(data_g_reg,3);
    figure;imagesc(data_avg);colormap gray
    save(fullfile(fnout,'regOuts&Img.mat'),'outs','regImg','data_avg')
    input = concatenateStructuresLG(temp);    
    save(fullfile(fnout,'input.mat'),'input')
    clear data_g
else %if not, must register. Start by showing average for each of four 500-frame bins so user can choose one
    nep = floor(size(data_g,3)./3000);
    [n n2] = subplotn(nep);
    figure; 
    movegui('center')
    for i = 1:nep 
        subplot(n,n2,i); 
        imagesc(mean(data_g(:,:,1+((i-1)*3000):100+((i-1)*3000)),3)); 
        title([num2str(1+((i-1)*3000)) '-' num2str(50+((i-1)*30000))]); 
        colormap gray; 
        clim([100 3000]); 
    end
    drawnow;
    
    regImgStartFrame = input('Enter Registration Image Start Frame, ENTER INTO DS:');
    
    regImg = mean(data_g(:,:,regImgStartFrame:(regImgStartFrame+49)),3);
    [outs,data_g_reg] = stackRegister(data_g,regImg);
    data_avg = mean(data_g_reg,3);
    figure;imagesc(data_avg);colormap gray; truesize; clim([100 3000]);
    print(fullfile(fnout,'avgFOV.pdf'),'-dpdf','-bestfit')
    clear data_g
    save(fullfile(fnout,'regOuts&Img.mat'),'outs','regImg','data_avg')
    input = concatenateStructuresLG(temp);    
    save(fullfile(fnout,'input.mat'),'input')
end
%
%reg red data 
%register the red data from the 920 nm run (same run used for green
%above)to the output of the green registration
% if info.config.pmt1_gain > 0.5
%     [~,data_r_reg] = stackRegister_MA(data_r,[],[],double(outs));
%     redChImg = mean(data_r_reg,3);
%     clear data_r clear data_r_reg
% end
    
% find activated cells
%find number of frames per trial and temporarily reshape data into trials
%overal goal here is to get green data in terms of df/f
nOn = input.nScansOn;
nOff = input.nScansOff;
sz = size(data_g_reg);
ntrials = size(input.tGratingContrast,2);
%ntrials = 374;
data_g_trial = reshape(data_g_reg, [sz(1) sz(2) nOn+nOff ntrials]);
data_g_f = squeeze(mean(data_g_trial(:,:,nOff/2:nOff,:),3));
data_g_on = squeeze(mean(data_g_trial(:,:,nOff+2:nOff+nOn,:),3));
data_g_dfof = (data_g_on-data_g_f)./data_g_f;
clear data_g_trial data_g_on data_g_f

%find the different contrasts and orientations
tCon = celleqel2mat_padded(input.tGratingContrast(1:ntrials));
%tCon = tCon(1:160);
cons = unique(tCon);
nCon = length(cons);
ind_con = find(tCon == max(cons(:)));
tDir = celleqel2mat_padded(input.tGratingDirectionDeg(1:ntrials));
%tDir = tDir(1:160);
tOri = tDir;
tOri(find(tDir>=180)) = tDir(find(tDir>=180))-180;
oris = unique(tOri);
nOri = length(oris);
data_g_ori = zeros(sz(1),sz(2), nOri);
data_temp = zeros(sz(1),sz(2), nOri, nCon);
[n n2] = subplotn(nOri);
figure; movegui('center');
for iOri = 1:nOri %for every orientation
    data_g_con = zeros(sz(1),sz(2), nCon); %make a data frame is the size of one imaging frame X the number of contrasts
    for iCon = 1:nCon %for every contrast
        ind_con = find(tCon == cons(iCon)); %find the indices of trials with that contrast
        ind_ori = intersect(ind_con, find(tOri == oris(iOri)));%find ind of intersection between that contrast and the ori we're looking at
        data_g_con(:,:,iCon) = mean(data_g_dfof(:,:,ind_ori),3); %pull out all those trials and average over the trials
        data_temp(:,:,iOri,iCon) = mean(data_g_dfof(:,:,ind_ori),3);
    end
    data_g_ori(:,:,iOri) = max(data_g_con,[],3); %find the contrast with the max response for that orientation
    subplot(n,n2,iOri)
    imagesc(data_g_ori(:,:,iOri))
    title(num2str(oris(iOri)))
end

%make a colored image comparing cell activity at the low (1-2), medium (3)
%and high (4-5) contrasts. There is a non-linear relationship between
%contrast and cells activated due to surround suppression. 
% rgb(:,:,1) = squeeze(max(mean(data_temp(:,:,:,1:2),4),[],3));
% rgb(:,:,2) = squeeze(max(mean(data_temp(:,:,:,3),4),[],3));
% rgb(:,:,3) = squeeze(max(mean(data_temp(:,:,:,4:5),4),[],3));
% figure;  movegui('center'); imagesc(rgb./max(max(rgb(:,:,3)))); 
% title('Comparing contrasts');

data_ori_max = max(data_g_ori,[],3);
data_dfof = cat(3, data_ori_max,data_ori_max,data_ori_max,data_g_ori);
figure; imagesc(data_ori_max); movegui('center');title('data ori max');
clear data_g_dfof

%get pixel correlation image, another method to identify cells
data_g_down = stackGroupProject(data_g_reg,100);
corrImg = getPixelCorrelationImage(data_g_down);
figure; imagesc(corrImg); movegui('center');title('pixel correlation');
data_dfof = cat(3, data_dfof, data_avg,data_avg,corrImg,corrImg);
clear data_g_down



%% load red cells
%this is where we use the 1040, 1000-frame run
if exist(fullfile(fnout,'redImage.mat'))
    load(fullfile(fnout,'redImage'))
elseif ~isempty(expt(day_id).redChannelRun) %if there IS a red channel run, find and load it
    redRun = expt(day_id).redChannelRun;
    imgMatFile = [redRun '_000_000.mat'];
    if strcmp(expt(day_id).data_loc,'lindsey')
        cd(fullfile(root,mouse, expDate,redRun));
    elseif strcmp(expt(day_id).data_loc,'ashley')
        cd(fullfile(root,mouse,'two-photon imaging', expDate,redRun));
    elseif strcmp(expt(day_id).data_loc,'tammy')
        cd(fullfile(root, mouse, '2P',expDate, redRun));
    elseif strcmp(expt(day_id).data_loc,'ACh')
        root = rc.achData;
        cd(fullfile(root, mouse, expDate, redRun));
    end
    load(imgMatFile);

    fprintf(['Reading run ' num2str(irun) '- ' num2str(info.config.frames) ' frames \r\n'])
    data_temp = sbxread(imgMatFile(1,1:11),0,info.config.frames);
    if size(data_temp,1) == 2
        data_rg = squeeze(data_temp(1,:,:,:));
        data_rr = squeeze(data_temp(2,:,:,:));
    else
        data_rr = squeeze(data_temp(1,:,:,:));
    end
    clear data_temp

    if exist('redChImg')
        [out, data_rr_reg] = stackRegister(stackGroupProject(data_rr,100), redChImg);
        redChImg = mean(data_rr_reg,3);
    elseif info.config.pmt0_gain>0.5 %if there is a green channel in this run, it gets registered to the registration image from green channel from the 920 run
        %data_rr = padarray(data_rr,9,0,'pre');
        redAvg = mean(data_rr,3);%take the average of the red at 1040
        [out, data_rr_reg] = stackRegister(data_rr,redAvg); %register red at 1040 to average of itself
        [~, data_rg_reg] = stackRegister_MA(data_rg,[],[],out); %register green at 1040 using the shifts from the registration of the red at 1040
        redChImgTemp = mean(data_rr_reg,3);  %take an average of the registered red at 1040
        rg_avg = mean(data_rg_reg,3); %take an average of the registered green at 1040
        [out2, ~] = stackRegister(rg_avg,data_avg); %register that green average to the registered data from my experimental run
        [~,redChImg]=stackRegister_MA(redChImgTemp,[],[],out2); %apply those shifts to the mean registered red at 1040
        
%         [out, data_rg_reg] = stackRegister(data_rg,data_avg); %register the green channel from the 1040 run to the green channel from the 920 run
%         [~, data_rr_reg]=stackRegister_MA(data_rr,[],[],out); %use those shifts to register the red 1040 run
%         redChImg = mean(data_rr_reg,3);
        
        
    else %if there is no green channel in this run
        redAvg = mean(data_rr,3);
        [out, data_rr_reg] = stackRegister(data_rr,redAvg);
        redChImgTemp = mean(data_rr_reg,3);
        [~,redChImg] = stackRegister(redChImgTemp,data_avg);
    end
    
    figure; colormap gray; imagesc(redChImg);  movegui('center');title('registration image for red channel');
   
    rgb = zeros(sz(1),sz(2),3);
    rgb(:,:,1) = redChImg./max(redChImg(:));
    rgb(:,:,2) = regImg./max(regImg(:));
    figure; image(rgb);  movegui('center')
    title('Green-920 + Red-1040')
    print(fullfile(fnout,'red_green_FOV.pdf'),'-dpdf','-bestfit')

    
    save(fullfile(fnout,'redImage'),'redChImg')
elseif ~exist('redChImg')
    redChImg = zeros(size(regImg));
end


%create red image where any pixel value above a certain percentile of the max is set to 90%
%of the max - removing the highest 10% of pixel values to create a lower
%contrast image for segmenting
threshPercentile = 99;

highValues = find(redChImg>prctile(redChImg,threshPercentile,'all'));
redThresh = redChImg;
redThresh(highValues)=prctile(redChImg,threshPercentile,'all');
figure; imagesc(redChImg);colormap gray;
figure; imagesc(redThresh);colormap gray;

clear data_rr data_rg data_rg_reg data_rr_reg
%% segment cells
close all

redForSegmenting = cat(3, redChImg,redChImg,redChImg); %make a dataframe that repeats the red channel image twice
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
%find and label the red cells - this is the first segmentation figure that
%comes up
if ~isempty(expt(day_id).redChannelRun)
    for iStim=1:size(redForSegmenting,3)
        mask_data_temp=redForSegmenting(:,:,iStim);
        mask_data_temp(find(mask_exp >= 1)) = 0;
        bwout = imCellEditInteractiveLG(mask_data_temp);
        mask_all = mask_all+bwout; 
        mask_exp = imCellBuffer(mask_all,3)+mask_all;
        close all
    end
end

%this version does not pad the red cells


mask_cell_red = bwlabel(mask_all);
mask_data = data_dfof; %this is the registered data from the 920 run
%after making masks for all the red cells, go through different stimuli and
%identify cells that are visible for each

for iStim = 1:size(data_dfof,3)
    mask_data_temp = mask_data(:,:,iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0;
    bwout = imCellEditInteractiveLG(mask_data_temp);
    mask_all = mask_all+bwout;
    mask_exp = imCellBuffer(mask_all,3)+mask_all;
    close all
end
mask_cell = bwlabel(mask_all);
figure; imagesc(mask_cell) 

nCells = max(mask_cell(:));
%mask_label = ones(1,nCells); %this is for the EMX and similar lines only
mask_label = zeros(1,nCells);
for i = 1:nCells
    if mask_cell_red(find(mask_cell == i, 1))
        mask_label(1,i) = 1;
    end
end

mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile(fnout, 'mask_cell.mat'), 'data_dfof', 'mask_cell', 'mask_cell_red', 'mask_np','mask_label')


rgb = zeros(sz(1),sz(2),3);
    rgb(:,:,1) = redChImg./(max(redChImg(:))*.5);
    rgb(:,:,2) = regImg./max(regImg(:));
    figure; image(rgb);  movegui('center')

hold on
bound = cell2mat(bwboundaries(mask_cell_red));
%plot(bound(:,2),bound(:,1),'.','color','b','MarkerSize',2);
hold off


%% extract timecourses
%data_g_reg = registered full image stack
%mask_cells = set of masks you segmented

data_tc = stackGetTimeCourses(data_g_reg, mask_cell); %data_tc timecourses of registered data based on masks you made
nCells = size(data_tc,2);
data_tc_down = stackGetTimeCourses(stackGroupProject(data_g_reg,5), mask_cell);
%stackGetTimecourses = fxn to turn masks and image stack into TCs
clear np_tc np_tc_down
sz = size(data_g_reg);
down = 5;
data_reg_down  = stackGroupProject(data_g_reg,down);
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_g_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['     Cell #' num2str(i) '%s/n']) 
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

save(fullfile(fnout, 'TCs.mat'), 'data_tc','np_tc','npSub_tc')

%clear data_g_reg data_reg_down
%% reshape by trials
%getting df/f for each trial, using a baseline window
data_tc_trial = reshape(npSub_tc, [nOn+nOff,ntrials,nCells]);
data_f_trial = mean(data_tc_trial(nOff/2:nOff,:,:),1);
data_dfof_trial = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial, data_f_trial), data_f_trial);

%split into baseline and response windows, run paired t-test to see if
%cells have a significant response (elevation in df/f comapred to baseline)
%for any orientations/contrasts
resp_win = nOff+2:nOn+nOff;
base_win = nOff/2:nOff;
data_resp = zeros(nCells, nOri, nCon,2);
h = zeros(nCells, nOri, nCon);
p = zeros(nCells, nOri, nCon);
for iOri = 1:nOri
    ind_ori = find(tOri == oris(iOri));
    for iCon = 1:nCon
        ind_con = find(tCon == cons(iCon));
        ind = intersect(ind_ori,ind_con); %for every orientation and then every contrast, find trials with that con/ori combination
        data_resp(:,iOri,iCon,1) = squeeze(mean(mean(data_dfof_trial(resp_win,ind,:),1),2));
        data_resp(:,iOri,iCon,2) = squeeze(std(mean(data_dfof_trial(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
        [h(:,iOri,iCon), p(:,iOri,iCon)] = ttest(mean(data_dfof_trial(resp_win,ind,:),1), mean(data_dfof_trial(base_win,ind,:),1),'dim',2,'tail','right','alpha',0.05./(nOri.*3-1));
    end
end
%%
h_all = sum(sum(h,2),3);

resp=logical(h_all);
red=mask_label';
resp_red=logical(resp.*red);
sum(resp)
sum(resp_red)

if length(find(n))<36
    [n n2] = subplotn(length(find(h_all)));
    tot = n.*n2;
else
    n = 6;
    n2 = 6;
    tot = 36;
end

% plot ori tuning at each contrast
figure;
movegui('center')
start = 1;
pref_ori = zeros(1,nCells);
for iCell = 1:nCells
    if start>tot
        figure; movegui('center')
        start = 1;
    end
    subplot(n,n2,start)
    %if find(find(h_all)==iCell)
        for iCon = 1:nCon
            errorbar(oris, data_resp(iCell,:,iCon,1), data_resp(iCell,:,iCon,2),'-o')
            hold on
        end
        if find(find(mask_label)==iCell)
            title('R')
        end
        start= start+1;
        ylim([-0.1 inf])
    %end
    [max_val, pref_ori(1,iCell)] = max(mean(data_resp(iCell,:,:,1),3),[],2);
end

% plots contrast preference at preferred orientation
figure;
movegui('center')
start = 1;
for iCell = 1:nCells
    if start>tot
        figure; movegui('center')
        start = 1;
    end
    subplot(n,n2,start)
    if find(find(h_all)==iCell)
        errorbar(cons, squeeze(data_resp(iCell,pref_ori(iCell),:,1)), squeeze(data_resp(iCell,pref_ori(iCell),:,2)),'-o')
        if find(find(mask_label)==iCell)
            title('R')
        end
        ylim([-0.1 inf])
        start = start+1;
    end
end
% %% LG tuning fit
% [avgResponseEaOri,semResponseEaOri,vonMisesFitAllCellsAllBoots,fitReliability,R_square,tuningTC]...
%     = getOriTuningLG(npSub_tc,input,5);
% %% extracting output from vonM fits
% % the output is in a 181 orientations X 1001 bootstraps X nCells matrix,
% % and there is a fitReliability output that tells the 90th percentile of
% % the difference between the bootstrapped fits and the original fit.
% % Ultimately I probably only want to take cells that have a fitReliability
% % value below 22.5 but right now none of them do
% median(fitReliability,2)
% 
% prefOris = nan(size(vonMisesFitAllCellsAllBoots,2),nCells);
% ogFitPref = nan(1,nCells);
% for i = 1:nCells
%     for fit = 2:size(vonMisesFitAllCellsAllBoots,2)
%         prefOris(fit,i) = mean(find(vonMisesFitAllCellsAllBoots(:,fit,i)==max(vonMisesFitAllCellsAllBoots(:,fit,i))));
%         ogFitPref(i) = mean(find(vonMisesFitAllCellsAllBoots(:,1,i)==max(vonMisesFitAllCellsAllBoots(:,1,i))));
%     end
% end
% meanPref = mean(prefOris,1);
% %I now have meanPref which tells me the average orientation indicated by
% %the 1000 bootstraps, and the max from the original fit. I'm not sure how
% %to get the oringal U value, or the K value
% vonMisesFitAllCells = squeeze(vonMisesFitAllCellsAllBoots(:,1,:)); % takes just the first of the bootstraps, which is the original fit 
% [maxResp prefOri] = max(vonMisesFitAllCells,[],1);
% 
% %combine max Resp, prefOri, and fit reliability 
% vonM_output = [maxResp; prefOri;fitReliability];
% 
% save(fullfile(fn,'vonM_fits.mat'),'vonMisesFitAllCells');
% save(fullfile(fn,'vonM_output.mat'),'vonM_output');

%% looking at time courses
red_tcs = npSub_tc(:,find(mask_label));
green_inds = 1:nCells;
green_inds = setdiff(green_inds, find(mask_label));
green_tcs = npSub_tc(:,green_inds);

data_tc_trial = reshape(npSub_tc, [nOn+nOff,ntrials,nCells]);

data_f_trial = mean(data_tc_trial(nOff/2:nOff,:,:),1);
data_dfof_trial = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial, data_f_trial), data_f_trial);

%looking at data with np subtracted
tc_cell_avrg = mean(data_dfof_trial(:,:,resp),3);%average pver cells, one row per trial
tc_trial_avrg = squeeze(mean(data_dfof_trial(:,:,resp),2));%average over trials, one row per cell
tc_cell_trial_avrg = mean(tc_cell_avrg,2);%average over trials and cells

figure;
plot(tc_trial_avrg, 'LineWidth',.005,'color',[.25 .25 .25]);
hold on;
plot(tc_cell_trial_avrg, 'LineWidth',2, 'color','k');
hold on;
vline(60,'g')
title('');
hold off
ylim([-.04 .1])



%% 
% figure;
% imagesc(data_avg);
% colormap gray
% title('average FOV');
% hold on
% bound = cell2mat(bwboundaries(mask_cell(:,:,1)));
% plot(bound(:,2),bound(:,1),'.','color','g','MarkerSize',2);
% bound = cell2mat(bwboundaries(mask_cell_red(:,:,1)));
% plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
% hold off

%% plotting contrast and ori functions for red cells only
% if length(find(mask_label))<36
%     [n n2] = subplotn(length(find(mask_label)));
%     tot = n.*n2;
% else
%     n = 6;
%     n2 = 6;
%     tot = 36;
% end
% 
% % plot ori tuning at each contrast
% figure;
% movegui('center')
% start = 1;
% pref_ori = zeros(1,nCells);
% for iCell = 1:nCells
%     if start>tot
%         figure; movegui('center')
%         start = 1;
%     end
%     subplot(n,n2,start)
%     if find(find(mask_label)==iCell)
%         for iCon = 1:nCon
%             errorbar(oris, data_resp(iCell,:,iCon,1), data_resp(iCell,:,iCon,2),'-o')
%             hold on
%         end
%         start= start+1;
%         ylim([-0.1 inf])
%     end
%     [max_val, pref_ori(1,iCell)] = max(mean(data_resp(iCell,:,:,1),3),[],2);
% end
% 
% % plots contrast preference at preferred orientation
% figure;
% movegui('center')
% start = 1;
% for iCell = 1:nCells
%     if start>tot
%         figure; movegui('center')
%         start = 1;
%     end
%     subplot(n,n2,start)
%     if find(find(mask_label)==iCell)
%         errorbar(cons, squeeze(data_resp(iCell,pref_ori(iCell),:,1)), squeeze(data_resp(iCell,pref_ori(iCell),:,2)),'-o')
%         ylim([-0.1 inf])
%         start = start+1;
%     end
% end
% 
% %% plotting contrast and ori functions for green cells only
% if length(green_inds)<36
%     [n n2] = subplotn(length(green_inds));
%     tot = n.*n2;
% else
%     n = 6;
%     n2 = 6;
%     tot = 36;
% end
% 
% % plot ori tuning at each contrast
% figure;
% movegui('center')
% start = 1;
% pref_ori = zeros(1,length(green_inds));
% for indx = 1:length(green_inds)
%     iCell = green_inds(indx);
%     if start>tot
%         figure; movegui('center')
%         start = 1;
%     end
%     subplot(n,n2,start)
%         for iCon = 1:nCon
%             errorbar(oris, data_resp(iCell,:,iCon,1), data_resp(iCell,:,iCon,2),'-o')
%             hold on
%         end
%         start= start+1;
%         ylim([-0.1 inf])
%     [max_val, pref_ori(1,iCell)] = max(mean(data_resp(iCell,:,:,1),3),[],2);
% end
% 
% % plots contrast preference at preferred orientation
% figure;
% movegui('center')
% start = 1;
% for indx = 1:length(green_inds)
%     iCell = green_inds(indx);
%     if start>tot
%         figure; movegui('center')
%         start = 1;
%     end
%     subplot(n,n2,start)
%         errorbar(cons, squeeze(data_resp(iCell,pref_ori(iCell),:,1)), squeeze(data_resp(iCell,pref_ori(iCell),:,2)),'-o')
%         ylim([-0.1 inf])
%         start = start+1;
%     
% end
