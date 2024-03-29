clear all; clear global; close all
clc
ds = 'DART_V1_contrast_ori'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories
eval(ds)
doGreenOnly = false;
doCorrImg = true;

day_id = 15;
%% load data for day

mouse = expt(day_id).mouse;
expDate = expt(day_id).date;

fn = fullfile(rc.analysis,mouse,expDate); %can make this flexible if folder structure is different
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
    end
    cd(CD);

    imgMatFile = [imgFolder '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\' dat mouse '-' expDate '-' times{irun} '.mat'];
    load(fName);

    temp(irun) = input;
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

%% register data for each day
%reg green data
if exist(fullfile(fnout,'regOuts&Img.mat'))
    load(fullfile(fnout,'regOuts&Img.mat'))
    [~,data_g_reg] = stackRegister_MA(data_g,[],[],double(outs));
    data_avg = mean(data_g_reg,3);
    figure;imagesc(data_avg);colormap gray
    save(fullfile(fnout,'regOuts&Img.mat'),'outs','regImg','data_avg')
    input = concatenateStructuresLG(temp);    
    save(fullfile(fnout,'input.mat'),'input')
    clear data_g
else
    nep = floor(size(data_g,3)./10000);
    [n n2] = subplotn(nep);
    figure; 
    movegui('center')
    for i = 1:nep 
        subplot(n,n2,i); 
        imagesc(mean(data_g(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); 
        title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); 
        colormap gray; 
        clim([0 3000]); 
    end
    regImgStartFrame = input('Enter Registration Image Start Frame, ENTER INTO DS:');
    regImg = mean(data_g(:,:,regImgStartFrame:(regImgStartFrame+499)),3);
    [outs,data_g_reg] = stackRegister(data_g,regImg);
    data_avg = mean(data_g_reg,3);
    figure;imagesc(data_avg);colormap gray; truesize; clim([0 3000]);
    print(fullfile(fnout,'avgFOV.pdf'),'-dpdf','-bestfit')
    clear data_g
    save(fullfile(fnout,'regOuts&Img.mat'),'outs','regImg','data_avg')
    input = concatenateStructuresLG(temp);    
    save(fullfile(fnout,'input.mat'),'input')
end

%reg red data 
if info.config.pmt1_gain > 0.5
    [~,data_r_reg] = stackRegister_MA(data_r,[],[],double(outs));
    redChImg = mean(data_r_reg,3);
    clear data_r clear data_r_reg
end
    
%% find activated cells
nOn = input.nScansOn;
nOff = input.nScansOff;
sz = size(data_g_reg);
ntrials = size(input.tGratingContrast,2);
data_g_trial = reshape(data_g_reg, [sz(1) sz(2) nOn+nOff ntrials]);
data_g_f = squeeze(mean(data_g_trial(:,:,nOff/2:nOff,:),3));
data_g_on = squeeze(mean(data_g_trial(:,:,nOff+5:nOff+nOn,:),3));
data_g_dfof = (data_g_on-data_g_f)./data_g_f;
clear data_g_trial data_g_on data_g_f

tCon = celleqel2mat_padded(input.tGratingContrast);
cons = unique(tCon);
nCon = length(cons);
ind_con = find(tCon == max(cons(:)));
tDir = celleqel2mat_padded(input.tGratingDirectionDeg);
tOri = tDir;
tOri(find(tDir>=180)) = tDir(find(tDir>=180))-180;
oris = unique(tOri);
nOri = length(oris);
data_g_ori = zeros(sz(1),sz(2), nOri);
data_temp = zeros(sz(1),sz(2), nOri, nCon);
[n n2] = subplotn(nOri);
figure; movegui('center')
for iOri = 1:nOri
    data_g_con = zeros(sz(1),sz(2), nCon);
    for iCon = 1:nCon
        ind_con = find(tCon == cons(iCon));
        ind_ori = intersect(ind_con, find(tOri == oris(iOri)));
        data_g_con(:,:,iCon) = mean(data_g_dfof(:,:,ind_ori),3);
        data_temp(:,:,iOri,iCon) = mean(data_g_dfof(:,:,ind_ori),3);
    end
    data_g_ori(:,:,iOri) = max(data_g_con,[],3);
    subplot(n,n2,iOri)
    imagesc(data_g_ori(:,:,iOri))
    title(num2str(oris(iOri)))
end

rgb(:,:,1) = squeeze(max(mean(data_temp(:,:,:,1:2),4),[],3));
rgb(:,:,2) = squeeze(max(mean(data_temp(:,:,:,3),4),[],3));
rgb(:,:,3) = squeeze(max(mean(data_temp(:,:,:,4:5),4),[],3));
figure;  movegui('center'); imagesc(rgb./max(max(rgb(:,:,3)))); 

data_ori_max = max(data_g_ori,[],3);
data_dfof = cat(3, data_ori_max,data_g_ori);
figure; imagesc(data_ori_max); movegui('center')
clear data_g_dfof

data_g_down = stackGroupProject(data_g_reg,100);
corrImg = getPixelCorrelationImage(data_g_down);
figure; imagesc(corrImg); movegui('center')
data_dfof = cat(3, data_dfof, corrImg);
clear data_g_down

%% load red cells
if exist(fullfile(fnout,'redImage.mat'))
    load(fullfile(fnout,'redImage'))
elseif ~isempty(expt(day_id).redChannelRun)
    redRun = expt(day_id).redChannelRun;
    imgMatFile = [redRun '_000_000.mat'];
    if strcmp(expt(day_id).data_loc,'lindsey')
        cd(fullfile(root,mouse, expDate,redRun));
    elseif strcmp(expt(day_id).data_loc,'ashley')
        cd(fullfile(root,mouse,'two-photon imaging', expDate,redRun));
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
    elseif info.config.pmt0_gain>0.5
        [out, data_rg_reg] = stackRegister(data_rg, regImg);
        [~, data_rr_reg] = stackRegister_MA(data_rr,[],[],out);
        redChImg = mean(data_rr_reg,3); 
    else
        [out, data_rr_reg] = stackRegister(stackGroupProject(data_rr,100), regImg);
        redChImg = mean(data_rr_reg,3);
    end
    figure; colormap gray; imagesc(redChImg);  movegui('center')
    rgb = zeros(sz(1),sz(2),3);
    rgb(:,:,1) = redChImg./max(redChImg(:));
    rgb(:,:,2) = regImg./max(regImg(:));
    figure; image(rgb);  movegui('center')
    title('Green-920 + Red-1020')
    rgb(:,:,1) = redChImg./max(redChImg(:));
    gImg = mean(data_rg_reg,3);
    rgb(:,:,2) = gImg./max(gImg(:));
    figure; image(rgb);  movegui('center')
    title('Green-1020 + Red-1020')
    save(fullfile(fnout,'redImage'),'redChImg')
elseif ~exist('redChImg')
    redChImg = zeros(size(regImg));
end
clear data_rr data_rg data_rg_reg data_rr_reg
%% segment cells
close all
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));

if ~isempty(expt(day_id).redChannelRun)
    bwout = imCellEditInteractiveLG(redChImg);
    mask_all = mask_all+bwout;
    mask_exp = imCellBuffer(mask_all,3)+mask_all;
    close all
end

mask_cell_red = bwlabel(mask_all);
mask_data = data_dfof;

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
mask_label = zeros(1,nCells);
for i = 1:nCells
    if mask_cell_red(find(mask_cell == i, 1))
        mask_label(1,i) = 1;
    end
end

mask_np = imCellNeuropil(mask_cell, 3, 5);
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
save(fullfile(fnout, 'TCs.mat'), 'data_tc','np_tc','npSub_tc')

clear data_g_reg data_reg_down
%% 
data_tc_trial = reshape(npSub_tc, [nOn+nOff,ntrials,nCells]);
data_f_trial = mean(data_tc_trial(nOff/2:nOff,:,:),1);
data_dfof_trial = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial, data_f_trial), data_f_trial);

resp_win = nOff+5:nOn+nOff;
base_win = nOff/2:nOff;
data_resp = zeros(nCells, nOri, nCon,2);
h = zeros(nCells, nOri, nCon);
p = zeros(nCells, nOri, nCon);
for iOri = 1:nOri
    ind_ori = find(tOri == oris(iOri));
    for iCon = 1:nCon
        ind_con = find(tCon == cons(iCon));
        ind = intersect(ind_ori,ind_con);
        data_resp(:,iOri,iCon,1) = squeeze(mean(mean(data_dfof_trial(resp_win,ind,:),1),2));
        data_resp(:,iOri,iCon,2) = squeeze(std(mean(data_dfof_trial(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
        [h(:,iOri,iCon), p(:,iOri,iCon)] = ttest(mean(data_dfof_trial(resp_win,ind,:),1), mean(data_dfof_trial(base_win,ind,:),1),'dim',2,'tail','right','alpha',0.05./(nOri.*3-1));
    end
end
h_all = sum(sum(h,2),3);
if length(find(h_all))<36
    [n n2] = subplotn(length(find(h_all)));
    tot = n.*n2;
else
    n = 6;
    n2 = 6;
    tot = 36;
end
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
    if find(find(h_all)==iCell)
        for iCon = 1:nCon
            errorbar(oris, data_resp(iCell,:,iCon,1), data_resp(iCell,:,iCon,2),'-o')
            hold on
        end
        if find(find(mask_label)==iCell)
            title('R')
        end
        start= start+1;
        ylim([-0.1 inf])
    end
    [max_val, pref_ori(1,iCell)] = max(mean(data_resp(iCell,:,:,1),3),[],2);
end

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
