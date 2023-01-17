clear all; clear global;  close all
clc
%ds = 'con_ori_nonDART'
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories
eval(ds);
doGreenOnly = false;
doCorrImg = true;

day_id = 186;

refFolder = '002'; %ENTER THE FOLDER NUMBER, FROM THE SAME RECORDING SESSION, 
%TO IMPORT MASKS FROM


%% load data for day and register to the reference

mouse = expt(day_id).mouse;
expDate = expt(day_id).date;

fn = fullfile(rc.achAnalysis,mouse,expDate); %can make this flexible if folder structure is different
mkdir(fn)
fnref = fullfile(fn,refFolder);
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
    else
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
        %data_temp_r = squeeze(data_temp_g(2,:,:,:));
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
            %data_r = data_temp_r;
            %clear data_temp_r
        end
    else
        data_g = cat(3, data_g, data_temp_g);
        clear data_temp_g
        if info.config.pmt1_gain > 0.5
            %data_r = cat(3, data_r, data_temp_r);
            %clear data_temp_r
        end
    end
end

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
else %if not, must register to the reference now
    %load in reference data
    
    load(fullfile(fnref,'regOuts&Img.mat'))
    regImg = data_avg; clear data_avg; %rename the data_avg from the reference run, and then clear than variable
    [outs,data_g_reg] = stackRegister(data_g,regImg); %use regImg (which was the average of 
    %the registered data from the reference) to register the current data 
    data_avg = mean(data_g_reg,3);
    data_avg = mean(data_g_reg,3);
    figure;imagesc(data_avg);colormap gray; truesize; clim([100 3000]);
    print(fullfile(fnout,'avgFOV.pdf'),'-dpdf','-bestfit')
    clear data_g
    save(fullfile(fnout,'regOuts&Img.mat'),'outs','regImg','data_avg')
    input = concatenateStructuresLG(temp);    
    save(fullfile(fnout,'input.mat'),'input')
end
%%
% find activated cells
%find number of frames per trial and temporarily reshape data into trials
%overal goal here is to get green data in terms of df/f
nOn = input.nScansOn;
nOff = input.nScansOff;

% data_g_reg = data_g_reg(:,:,1:10080);
% ntrials = 160;
ntrials = size(input.tGratingContrast,2);
sz = size(data_g_reg);
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
data_dfof = cat(3, data_ori_max,data_g_ori);
figure; imagesc(data_ori_max); movegui('center');title('data ori max');
clear data_g_dfof

%get pixel correlation image, another method to identify cells
data_g_down = stackGroupProject(data_g_reg,100);
corrImg = getPixelCorrelationImage(data_g_down);
figure; imagesc(corrImg); movegui('center');title('pixel correlation');
data_dfof = cat(3, data_dfof, corrImg);
clear data_g_down

data_dfof = cat(3, data_dfof,data_avg);

%% load red cells and register them to the green signal
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
       redAvg = mean(data_rr,3);
        [out, data_rr_reg] = stackRegister(data_rr,redAvg);
        [~, data_rg_reg] = stackRegister_MA(data_rg,[],[],out);
        redChImgTemp = mean(data_rr_reg,3); 
        rg_avg = mean(data_rg_reg,3);
        [out2, ~] = stackRegister(rg_avg,data_avg);
        [~,redChImg]=stackRegister_MA(redChImgTemp,[],[],out2);
        
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
clear data_rr data_rg data_rg_reg data_rr_reg
%% import masks from a reference run
load(fullfile(fnref,'mask_cell.mat'));
nCells = max(mask_cell(:)); % take max label of mask_cell, should circumvent bwlabel
fprintf([num2str(nCells) ' total cells selected\n'])

figure;
imagesc(data_dfof(:,:,1)); hold on;
bound = cell2mat(bwboundaries(mask_cell(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','b','MarkerSize',.5); 
colormap(gray)
figure;
imagesc(data_dfof(:,:,3)); hold on;
bound = cell2mat(bwboundaries(mask_cell(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','b','MarkerSize',.5); 
colormap(gray)

save(fullfile(fnout, 'mask_cell.mat'), 'redChImg', 'data_dfof', 'mask_cell', 'mask_cell_red', 'mask_np','mask_label')

%% extract timecourses


data_tc = stackGetTimeCourses(data_g_reg, mask_cell);
nCells = size(data_tc,2);
data_tc_down = stackGetTimeCourses(stackGroupProject(data_g_reg,5), mask_cell);
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



h_all = sum(sum(h,2),3);

resp=logical(h_all);
red=mask_label';
resp_red=resp.*red;
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

%% looking at time courses
red_tcs = npSub_tc(:,find(mask_label));
green_inds = 1:nCells;
green_inds = setdiff(green_inds, find(mask_label));
green_tcs = npSub_tc(:,green_inds);

data_tc_trial = reshape(npSub_tc, [nOn+nOff,ntrials,nCells]);

data_f_trial = mean(data_tc_trial(nOff/2:nOff,:,:),1);
data_dfof_trial = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial, data_f_trial), data_f_trial);

%looking at data with np subtracted
tc_cell_avrg = mean(data_dfof_trial,3);%average pver cells, one row per trial
tc_trial_avrg = squeeze(mean(data_dfof_trial,2));%average over trials, one row per cell
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
figure;
imagesc(data_avg);
colormap gray
title('average FOV');
hold on
bound = cell2mat(bwboundaries(mask_cell(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','g','MarkerSize',2);
bound = cell2mat(bwboundaries(mask_cell_red(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
hold off

