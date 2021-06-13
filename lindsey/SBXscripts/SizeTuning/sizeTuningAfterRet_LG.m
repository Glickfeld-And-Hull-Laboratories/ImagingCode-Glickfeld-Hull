%% size tuning script by kevin
% modifying newScript.m, singleChannelTCScript.m
% reads in Size Tuning Curves from single channel imaging data

% here we extract size tuning, using the neuron mask from retOnly run

%% get path names
clear all;clc; close all;
ds = 'retAndSzTuning_ExptList';
iexp = 26; 
rc = behavConstsAV;
eval(ds)

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
ImgFolder = cell2mat(expt(iexp).szFolder);
time = cell2mat(expt(iexp).szTime);
nrun = length(ImgFolder);
frameRateHz = params.frameRate;
RetImgFolder = cell2mat(expt(iexp).retFolder);
doFromRef = 1;

nrun = size(ImgFolder,1);
frame_rate = params.frameRate;
run_str = catRunName(ImgFolder, nrun);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];

fprintf(['2P imaging size tuning analysis - by KM, Glickfeld Lab\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([ImgFolder(irun,:) ' - ' time(irun,:) '\n'])
end

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
lg_fn = fullfile(fn_base, 'home\lindsey');
data_fn = fullfile(lg_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data');
fnout = fullfile(lg_fn, 'Analysis\2P');


%% load data, read with sbxread, and concatenate selected runs
data = [];
clear temp
trial_n = zeros(1,nrun);

fprintf(['\nBegin reading ' num2str(nrun) ' runs...'])
for irun = 1:nrun
    % load 2p imaging data
    CD = fullfile(data_fn, mouse, date, ImgFolder(irun,:));
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    
    % load behavior/experimental data
    fName = fullfile(mworks_fn, ['data-' mouse '-' date '-' time(irun,:) '.mat']);
    load(fName);
    
    % read in frames with sbxread
    nframes = info.config.frames;
    fprintf(['\nReading run ' num2str(irun) ' - ' num2str(nframes) ' frames \n'])
    tic
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    toc
    
    fprintf('Reshaping data...\n')
    temp(irun) = input;
    
    % store values on nOn + nOff, and measure number of trials
    nOn = temp(irun).nScansOn;
    nOff = temp(irun).nScansOff;
    ntrials = size(temp(irun).tGratingDirectionDeg,2);
    
    % squeeze because only 1 pmt channel
    data_temp = squeeze(data_temp);
    
    % if nframes =/= ntrials*(frames in trial), then resize
    if nframes>ntrials*(nOn+nOff)
        fprintf('Too many frames, truncating...\n')
        data_temp = data_temp(:,:,1:ntrials*(nOn+nOff));
        nframes = (ntrials*(nOn+nOff))
    elseif nframes<ntrials*(nOn+nOff)
        fprintf('Too many trials, chop chop...\n')
        temp(irun) = trialChopper(temp(irun),1:ceil(nframes./(nOn+nOff)));
        ntrials = ceil(nframes./(nOn+nOff))
    end
    
    data = cat(3,data,data_temp);
    trial_n(irun) = ntrials;
    fprintf('Complete\n')
end
fprintf('All runs read\n')
fprintf([num2str(size(data,3)) ' total frames\n'])
input = concatenateDataBlocks(temp);
for i=1:length(input.tGratingContrast) % replace int64(con=1) with double
    if ~(class(input.tGratingContrast{i})=="double")
        input.tGratingContrast{i} = double(input.tGratingContrast{i});
    end
end
clear data_temp
clear temp

%% Register data
fprintf(['Reference specified: ' RetImgFolder '\n'])

ref_str = ['runs-' RetImgFolder];
datemouserefrun = [date '_' mouse '_' ref_str];
% load from folder specified by ref_str
fprintf(['Loading from folder: ' ref_str '\n'])
load(fullfile(fnout, datemouse, datemouserefrun, [datemouserefrun '_reg_shifts.mat']))
load(fullfile(fnout, datemouse, datemouserefrun, [datemouserefrun '_mask_cell.mat']))
%load(fullfile(fnout, datemouse, datemouserun, [datemouserun '_trialData.mat']))

% register
fprintf('stackRegister with reference image\n')
[out, data_reg] = stackRegister(data,data_avg);

% save
fprintf('Registration complete, now saving...\n')
mkdir(fullfile(fnout, datemouse, datemouserun))
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_reg_shifts.mat']), 'out', 'data_avg')
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_input.mat']), 'input')

clear data % depending on memory

%% test stability
% figure 2 shows the registered images to check the stability
fprintf('\nExamine registered images for stability\n')
figure(2);clf;
regIntv = 3000;
nep = floor(size(data_reg,3)./regIntv);
[n n2] = subplotn(nep);
for i = 1:nep
    subplot(n,n2,i);
    imagesc(mean(data_reg(:,:,(1:500)+((i-1)*regIntv)),3));
    title([num2str(i) ': ' num2str(1+((i-1)*regIntv)) '-' num2str(500+((i-1)*regIntv))]);
end

% figure 3 shows imagesq (scaled?) of the registered data 1-10000, then prints to FOV_avg.mat
fprintf('Examine FOV\n')
figure(3);
if nframes>10000
    imagesq(mean(data_reg(:,:,1:10000),3));
else
    imagesq(mean(data_reg,3));
end
truesize;
% print to file
set(gcf, 'Position', [0 0 800 1000]);
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_FOV_avg.pdf']))
sz = size(data_avg);
rgb = zeros(sz(1),sz(2),3);
rgb(:,:,1) = data_avg./max(data_avg(:));
rgb(:,:,2) = mean(data_reg,3)./max(max(mean(data_reg,3),[],2),[],1);
figure; imagesc(rgb)

%% find activated cells

% check useFilt (only for dFoF calculation and segmentation, timecourse extraction not filtered)
useFilt = 0;

% calculate dF/F
fprintf('\nBegin image analysis...\n')

% max by trial type
if isfield(input, 'nScansOn')
    % nScansOn -> passive vis stim ret
    fprintf('nScansOn method - get dF/F\n')
    
    % load defining variables
    nOn = input.nScansOn;
    nOff = input.nScansOff;
    ntrials = size(input.tGratingDirectionDeg,2);
    sz = size(data_reg);
    
    fprintf('Calculating dF/F...\n')
    
    % split into trials by reshape
    data_tr = reshape(data_reg,[sz(1), sz(2), nOn+nOff, ntrials]);
    data_f = mean(data_tr(:,:,nOff/2:nOff,:),3);
    data_fOn = mean(data_tr(:,:,nOff:(nOn+nOff),:),3);
    data_dfof = squeeze((data_fOn-data_f)./data_f);
    
    % previous code would do data_dfof for every frame
    % data_dfof = (double(data_tr)-data_f)./data_f;
    
    clear data_tr data_f data_fOn % data_fLate
    
    if useFilt
        fprintf('Filtering images and time averaging...\n')
        switch useFilt
            case 1
                % filter with a 20x20 gaussian sigma=0.7
                fprintf('useFilt=1: Gaussian smoothing filter\n')
                myfilter = fspecial('gaussian',[20 20], 0.7);
                data_dfof = imfilter(data_dfof,myfilter);
            case 2
                % filter with median filter (default 3x3)
                fprintf('useFilt=2: Median filter\n')
                data_dfof = medfilt2(data_dfof, [3 3]);
            case 3
                % filter with wiener adaptive filter (default 3x3)
                fprintf('useFilt=3: Wiener adaptive filter\n')
                data_dfof = wiener2(data_dfof,[5 5]);
        end
    end
    
    fprintf('done\n')
end

% with dF/F, average by each stimulus, depending on experiment
if input.doRetStim
    % doRetStim -> retinotopy
    % requires data_dfof from above nScansOn method
    fprintf('input.doRetStim method - varying Az+El position\n')
    fprintf('ERROR: this script only does size tuning not retinotopy')
    
elseif input.doSizeStim
    % doSizeStim -> size tuning
    % requires data_dfof from above nScansOn method
    fprintf('input.doSizeStim method - varying grating diameter\n')
    
    % store vectors for grating diameter
    % what is celleqel2mat_padded?
    sz_mat = celleqel2mat_padded(input.tGratingDiameterDeg);
    szs = unique(sz_mat);
    nSize = length(szs);
    con_mat = celleqel2mat_padded(input.tGratingContrast);
    cons = unique(con_mat);
    nCon = length(cons);
    
    nStim = length(szs);
    fprintf([num2str(nStim) ' unique size stimuli\n'])
    
    fprintf('Averaging dF/F for each size...\n')
    data_dfof_avg = zeros(sz(1), sz(2), nStim);
    for isz = 1:nStim
        ind = find(sz_mat == szs(isz));
        data_dfof_avg(:,:,isz) = mean(data_dfof(:,:,ind),3);
    end
    
    % plot dF/F for all stimuli
    figure(5);clf;
    movegui('center')
    [n, n2] = subplotn(nStim);
    for i = 1:nStim
        subplot(n,n2,i);
        imagesc(data_dfof_avg(:,:,i));
        clim([0 max(data_dfof_avg(:))])
        title(['Size: ' num2str(szs(i)) ' deg'])
    end
    % print to pdf
    print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_FOV_resp_Size.pdf']), '-dpdf')
    
elseif input.doTFStim && ~input.doMatrix
    % doTFStim + !doMatrix -> temporal frequency?
    % requires data_dfof from above nScansOn method
    fprintf('input.doTFStim method - varying grating TF+SF\n')
    fprintf('ERROR: this script only does size tuning not TF tuning')
    
elseif input.doDirStim
    % doDirStim -> directional stimulation tuning curve
    % requires data_dfof from above nScansOn method
    fprintf('input.doDirStim method - varying grating direction\n')
    fprintf('ERROR: this script only does size tuning not direction tuning')
    
end

% take max across stimuli
fprintf('Final step: take max across stimuli\n')
data_dfof_max = max(data_dfof_avg,[],3);
figure(9);clf;
movegui('center')
if input.doDirStim
    imagesc(max(data_dfof_avg_ori,[],3))
    clim([0 max(data_dfof_avg_ori(:))])
else
    imagesc(data_dfof_max)
    clim([0 max(data_dfof_max(:))])
end
title('Maximum dF/F across all stimuli')

% save stimActFOV.mat containing: data_dfof_max, data_dfof_avg, nStim
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg', 'nStim')
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_stimData.mat']), 'sz_mat', 'szs', 'nSize', 'con_mat', 'cons', 'nCon');

%% load cell masks from retinotopy runs

%RetImgFolder = char('001','002','003');
nret = size(RetImgFolder,1);
ret_str = catRunName(RetImgFolder, nret);
fprintf(['Loading masks from retinotopy runs: ' ret_str '\n'])

% loads 'mask_cell', 'mask_np'
load(fullfile(fnout, datemouse, datemouserefrun, [datemouserefrun '_mask_cell.mat']))

fprintf('Cell and neuropil masks loaded\n')

% translate if necessary (should not be necessary after register with ref)
% mask_cell = imtranslate(mask_cell, [0 2]); % [x(+right) y(+down)]
% mask_np = imtranslate(mask_np, [0 2]);

% load ret fit data, in order to select only goodfit_ind cells
% fprintf(['Loading fits from retinotopy runs: ' ret_str '\n'])
% fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_str], [date '_' mouse '_' ret_str '_lbub_fits.mat']);
% load(fn_out);

% [mask_cell, nCells] = bwlabel(mask_cell); % bwlabel labels all individual cells
nCells = max(mask_cell(:)); % take max label of mask_cell, should circumvent bwlabel
fprintf([num2str(nCells) ' total cells selected\n'])
fprintf('Cell segmentation complete\n')

figure(11);clf;
movegui('center')
[n, n2] = subplotn(nStim);
for i = 1:nStim
    subplot(n,n2,i);
    shade_img = imShade(data_dfof_avg(:,:,i), mask_cell);
    imagesc(shade_img)
    if input.doSizeStim
        title([num2str(szs(i)) ' deg'])
    elseif input.doRetStim
        title(num2str(Stims(i,:)))
    end
    clim([0 max(data_dfof_avg(:))])
    colormap(gray)
end


%% Get time courses, including neuropil subtraction

fprintf('\nBegin time course extraction...\n')
down = 5; % downsample rate
data_reg_down  = stackGroupProject(data_reg,down);
fprintf(['Downsampling at M=' num2str(down) '\n'])

fprintf('Extracting cell signal...\n')
data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
nCells = size(data_tc,2);
fprintf([num2str(nCells) ' total cells extracted\n'])

fprintf('Extracting neuropil signal...\n')
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
    np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
    np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
    fprintf(['Cell #' num2str(i) ' / ' num2str(nCells) '\n'])
end

fprintf('Subtract neuropil signal, maximizing skewness\n')
%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
[max_skew, ind] =  max(x,[],1);
fprintf(['Maximum skews: ' num2str(max_skew) '\n'])
fprintf(['at inds: ' num2str(ind) '\n'])

np_w = 0.01*ind;
fprintf(['np_w = ' num2str(np_w) '\n'])

npSub_tc = data_tc-(tcRemoveDC(np_tc).*np_w);
clear data_reg_down data_reg

fprintf('Neuropil subtraction complete, saving data...\n')

save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_input.mat']), 'input')

fprintf('Calculating dF/F...\n')
%get dF/F
nCells = size(npSub_tc,2);
tc_mat = zeros(nOn+nOff, nCells, ntrials);
% reshape into trials, could do this in one line with reshape?
for itrial = 1:ntrials
    tc_mat(:,:,itrial) = npSub_tc(1+((itrial-1)*(nOn+nOff)):(itrial*(nOn+nOff)),:);
end
tc_f = mean(tc_mat(nOff/2:nOff,:,:),1);
tc_dfof = (tc_mat - tc_f) ./ tc_f;
clear tc_mat tc_f

fprintf('Time course extraction complete.\n')


%% now select window to extract response
fprintf(['Loading fits from retinotopy runs: ' RetImgFolder '\n'])
load(fullfile(fnout, datemouse, datemouserefrun, [datemouserefrun '_lbub_fits.mat']))
fprintf('\nExamine average cell dF/F timecourse to select response window\n')
respWindow = int64(7:15); % this will be the window to define dF/F response (offset by nITI + nDelay)
baseWindow = int64(-5:0);

tt = (1-nOff:nOn)*(1000./frame_rate);
avTC = squeeze(mean(mean(tc_dfof(:,goodfit_ind,:),2),3)); % average all cells,trials
% plots all trials for one condition + one cell at a time
figure(2);clf;
for i = 1:length(goodfit_ind)
    plot(tt', squeeze(mean(tc_dfof(:,goodfit_ind(i),:),3)), 'LineWidth',1)
    hold on
end
plot(tt', avTC, 'LineWidth', 5)
title('Average cell timecourses')
ylim([-0.05 0.6])
vline(0)
vline(tt(nOff+respWindow(1)))
vline(tt(nOff+respWindow(end)))
h1=vfill([tt(nOff+respWindow(1)),tt(nOff+respWindow(end))],'gray','FaceAlpha',0.5);
uistack(h1,'bottom')

clear data_reg

%% find good cells
cellAz = lbub_fits(goodfit_ind,4,4);
cellEl = lbub_fits(goodfit_ind,5,4);
fprintf('Retinotopy fits loaded, found cell receptive field coordinates\n')
nCells = length(goodfit_ind);
fprintf(['# goodfit cells = ' num2str(nCells) '\n'])

% input stimulus location based on experimental choice
stimEl = double(input.gratingElevationDeg);
stimAz = double(input.gratingAzimuthDeg);
fprintf(['Stimulus at: El ' num2str(stimEl) ', Az ' num2str(stimAz) '\n'])

% calculate cell distances
fprintf('Calculating cell RF distances to stimulus...\n')
cellDists = sqrt((cellAz-stimAz).^2+(cellEl-stimEl).^2);

n = length(find(cellDists<=10));
fprintf(['# goodfit cells within 10 deg = ' num2str(n) '\n'])

%% plot tuning
nSize = length(szs);

tuning_mat = zeros(nSize,nCon, 2, nCells);

Ind_struct = [];
if nCells<36
    [n, n2] = subplotn(nCells);
else
    [n, n2] = subplotn(36);
end
tt = (1-nOff:nOn)*(1000./frame_rate);
figure;
start = 1;
f = 1;
sizeTune = cell(nSize,nCon,nCells);
baseResp = cell(nSize,nCon,nCells);
for iCell = 1:nCells
    iC = goodfit_ind(iCell);
    if start >36
        set(gcf, 'Position', [0 0 800 1000]);
        print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_TCs' num2str(f) '.pdf']), '-dpdf')
        start = 1;
        f= f+1;
        figure;
    end
    subplot(n, n2, start)
    for iSize = 1:nSize
        ind_temp = find(sz_mat == szs(iSize));
        plot(tt', squeeze(mean(tc_dfof(:,iC,ind_temp),3)))
        title(iC)
        hold on
        for iCon = 1:nCon
            ind = intersect(ind_temp, find(con_mat==cons(iCon)));
            stimOn = mean(tc_dfof(nOff+(respWindow),iCell,ind),1);
            stimOff = mean(tc_dfof(nOff+(baseWindow),iCell,ind),1);
            % this cell matrix stores all stimOn for all trials in ind_all
            % dims (szs, cons, cells)
            sizeTune{iSize, iCon, iCell} = squeeze(stimOn);
            baseResp{iSize, iCon, iCell} = squeeze(stimOff);
            tuning_mat(iSize,iCon,1,iCell) = mean(stimOn,3);
            tuning_mat(iSize,iCon,2,iCell) = std(stimOn,[],3)./sqrt(length(ind));
            if iCell == 1
                Ind_struct(iSize).all_trials = ind;
            end
        end
    end
    ylim([-0.05 0.25])
    vline(nOff)
    start = start + 1;
end
set(gcf, 'Position', [0 0 800 1000]);
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_TCs' num2str(f) '.pdf']), '-dpdf')

figure;
start = 1;
f = 1;
for iCell = 1:nCells
    if start >36
        set(gcf, 'Position', [0 0 800 1000]);
        print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_Tuning' num2str(f) '.pdf']), '-dpdf')
        start = 1;
        f= f+1;
        figure;
    end
    subplot(n, n2, start)
    for iCon = 1:nCon
        ret_mat = tuning_mat(:,iCon,1,iCell);
        ret_mat = ret_mat';
        plot(szs,ret_mat)
        hold on 
    end
    start = start+1;
end
set(gcf, 'Position', [0 0 800 1000]);
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_Tuning' num2str(f) '.pdf']), '-dpdf')
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_Tuning.mat']), 'tc_dfof', 'tuning_mat', 'sizeTune', 'baseResp', 'szs', 'sz_mat', 'nSize', 'Ind_struct', 'respWindow', 'baseWindow','cellDists')

%% Fit size tuning curves
%integrated from sizeCurveFitting_SPbootstrap
override = 1; %override_all; % over-ride for creating new data (1=override, 0=skip if exist)

% runs through all cells at each contrast condition
% calls Fit_SizeTuneSmooth_KM script which outputs fit structure

% First model: single sigmoid (3 fit params)
% use 90% cutoff as measure of preferred size, no suppression index
% Second model: sum of + and - sigmoids (6 fit params)
% use peak of fit as preferred size, and measure suppression index

szRng = linspace(0,max(szs));

% single sigmoid fits Ae, ke=k1, xe=x1
logfit1 = @(coefs,xdata) coefs(1)./(1+exp(-coefs(2)*(xdata-coefs(3))))
%logfit1 = @(coefs,xdata) 2*coefs(1)./(1+exp(-coefs(2)*(xdata-coefs(3)))) - coefs(1)
% double sigmoid fits Ae, ke=k1+k2, xe=x1, Ai, ki=k2, xi=x1+x2
% ke and xi defined so ex curve is steeper and in curve is centered higher
logfit2 = @(coefs,xdata) coefs(1)./(1+exp(-(coefs(2)+coefs(5))*(xdata-coefs(3)))) - coefs(4)./(1+exp(-coefs(5)*(xdata-(coefs(3)+coefs(6)))))

% store # trials at each size for this run (same for all cells)
nTr = zeros(nSize,nCon);
for iSz = 1:nSize
    for iCon = 1:nCon
        nTr(iSz,iCon) = length(sizeTune{iSz,iCon,1});
    end
end

%cd('K:\Code')
opts = optimset('Display','off');
    
filename = fullfile(fnout, datemouse, datemouserun, [datemouserun '_sizeFitResults.mat']);
% if exist(filename, 'file') && ~override
%     fprintf('Found sizeFitResults_SP.mat, loading previous results...\n')
%     load(filename, 'sizeFits')
%     fprintf('Loaded sizeFits struct\n')
% else
Fit_struct = [];
Nshuf = 500;
fprintf(['Nshuf = ' num2str(Nshuf) '\n'])
sizeMean = squeeze(tuning_mat(:,:,1,:));
sizeSEM = squeeze(tuning_mat(:,:,2,:));
if ~exist('cellDists')
    cellDists = nan(1,nCells);
end
% store # trials at each size and highest con (same for all cells)
fprintf(['Sizes: ' num2str(szs) '\n# trials: ' num2str(nTr(:,nCon)')])
shuf_ind = cell(nSize,1);

fprintf('\nBegin shuffling...\n')
figure;

fprintf('Creating new size-tuning curve fit data...\n')
fprintf('Begin fitting size-tuning curves at all cells, all runs...')
for count_shuf = 0:Nshuf
    fprintf(['count_shuf: ' num2str(count_shuf) '/' num2str(Nshuf) '\n'])
    for iSz = 1:nSize
        for iCon = 1:nCon
            if count_shuf > 0
                shuf_ind{iSz,iCon} = randsample(nTr(iSz,iCon),nTr(iSz,iCon),1); % resample with replacement
            elseif count_shuf == 0
                shuf_ind{iSz,iCon} = 1:nTr(iSz,iCon); % shuf_count==0 use all trials
            end
        end
    end
    ifig = 1;
    start = 1;
    for iCell = 1:nCells
        fprintf(['\nCell# ' num2str(iCell) '/' num2str(nCells)]) %, (RF-stim dist: ' num2str(cellDists(iCell)) ' deg)'])
        
        for iCon = 1:nCon
            nPts = floor(mean(nTr(:,iCon)));
            dumdum = zeros(1,nPts);%[0]; % define zero point for dF/F
            szs0 = zeros(1,nPts);%[0.1]; % and size
            for iSz = 1:nSize
                nPts = nPts + nTr(iSz,iCon);
                dum = sizeTune{iSz,iCon,iCell}';
                dumdum = [dumdum dum(shuf_ind{iSz,iCon})];
                szs0 = [szs0 szs(iSz)*ones(1,nTr(iSz,iCon))];
            end

            % max of each size mean for use in initial guesses
            [maxMean maxVal] = max(tuning_mat(:,iCon,1,iCell));
            x_max = szs(maxVal);   

            if count_shuf == 0
                if iCon == nCon
                    PLOTIT_FIT = 1;
                else
                    PLOTIT_FIT = 0;
                end
                SAVEALLDATA = 1;
                Fit_SizeTuneSmooth_KM % call fit script, returns fit structure s, saves plots to kevin analysis folder
                eval(['Fit_struct(iCell,iCon).True.s_',' = s;']);
            else
                SAVEALLDATA = 0;
                PLOTIT_FIT = 0;
                Fit_SizeTuneSmooth_KM
                eval(['Fit_struct(iCell,iCon).Shuf(count_shuf).s_',' = s;']);
            end
        end
        start = start+1;
    end
    if count_shuf == 0 
        set(gcf, 'Position', [0 0 800 1000]);
        fn_out = fullfile(fnout, datemouse, datemouserun, [datemouserun '_SizeTuneFits' num2str(ifig) '.pdf']);
        print(fn_out,'-dpdf')
    end
end
fprintf('\nShuffling done, saving fit results\n')
props = whos('Fit_struct');
fn_out = fullfile(fnout, datemouse, datemouserun, [datemouserun '_Fit_struct.mat']);
if props.bytes>2000000000
    save(fn_out, 'Fit_struct','-v7.3')
else
    save(fn_out, 'Fit_struct')
end
%% assess fits
% extract values for prefSize (use for decision), Ftest
% also prefSize(1/2), suppInd(1,2), fit1.c1/OF1, fit2.c2/OF2, Rsq12

fprintf('Assessing goodness of fit\n')
fprintf('Reading in variables of interest\n')
fit_true_vec = NaN(nCells,10,nCon);
    for iCon = 1:nCon
        for iCell = 1:nCells
            if ~isempty(Fit_struct(iCell,iCon).True)
                eval('tmp = Fit_struct(iCell,iCon).True.s_.prefSize;');
                eval('tmp = [tmp Fit_struct(iCell,iCon).True.s_.prefSize1];');
                eval('tmp = [tmp Fit_struct(iCell,iCon).True.s_.prefSize2];');
                eval('tmp = [tmp Fit_struct(iCell,iCon).True.s_.suppInd];');
                eval('tmp = [tmp Fit_struct(iCell,iCon).True.s_.suppInd1];');
                eval('tmp = [tmp Fit_struct(iCell,iCon).True.s_.suppInd2];');
                eval('tmp = [tmp Fit_struct(iCell,iCon).True.s_.Fscore];');
                eval('tmp = [tmp Fit_struct(iCell,iCon).True.s_.Ftest];');
                eval('tmp = [tmp Fit_struct(iCell,iCon).True.s_.maxResp1];');
                eval('tmp = [tmp Fit_struct(iCell,iCon).True.s_.maxResp2];');

                % prefSize PS1 PS2 suppInd SI1 SI2 Fscore Ftest
                fit_true_vec(iCell,:,iCon) = tmp;
            end
        end
    end
    
    fit_shuf_vec = NaN(nCells,10,Nshuf,nCon);
    for iCon = 1:nCon
        for count_shuf = 1:Nshuf
            for iCell = 1:nCells
                if ~isempty(Fit_struct(iCell,iCon).Shuf)
                    eval('tmp = Fit_struct(iCell,iCon).Shuf(count_shuf).s_.prefSize;');
                    eval('tmp = [tmp Fit_struct(iCell,iCon).Shuf(count_shuf).s_.prefSize1];');
                    eval('tmp = [tmp Fit_struct(iCell,iCon).Shuf(count_shuf).s_.prefSize2];');
                    eval('tmp = [tmp Fit_struct(iCell,iCon).Shuf(count_shuf).s_.suppInd];');
                    eval('tmp = [tmp Fit_struct(iCell,iCon).Shuf(count_shuf).s_.suppInd1];');
                    eval('tmp = [tmp Fit_struct(iCell,iCon).Shuf(count_shuf).s_.suppInd2];');
                    eval('tmp = [tmp Fit_struct(iCell,iCon).Shuf(count_shuf).s_.Fscore];');
                    eval('tmp = [tmp Fit_struct(iCell,iCon).Shuf(count_shuf).s_.Ftest];');
                    eval('tmp = [tmp Fit_struct(iCell,iCon).Shuf(count_shuf).s_.maxResp1];');
                    eval('tmp = [tmp Fit_struct(iCell,iCon).Shuf(count_shuf).s_.maxResp2];');

                    % prefSize PS1 PS2 suppInd SI1 SI2 Fscore Ftest
                    fit_shuf_vec(iCell,:,count_shuf,iCon) = tmp;
                end
            end
        end
    end
    
    %% plot cells with size tuning curve and shuffle results
    %chosen=[44 54]; %[31 41 45 46 52 64 67 71 72 73 75 77 79 83 89];
    %chosen = goodfit_ind_size;
    chosen = [1];
    Npars = size(fit_shuf_vec,2);
    lbub_fits = NaN(nCells,Npars,5);
    alpha_bound = .025;
    ind_shuf_lb = ceil(Nshuf*alpha_bound); % 0.025 percentile
    ind_shuf_ub = ceil(Nshuf*(1-alpha_bound)); % 0.975 percentile
    for iCell = 1:nCells
        if sum(iCell==chosen)
            s = Fit_struct(iCell,nCon).True.s_;
            figure(1);clf;
            subplot(3,3,1)
            errorbar([0 szs],[0 sizeMean(:,nCon,iCell)'],[0 sizeSEM(:,nCon,iCell)'])
            hold on
            plot(s.szs0,s.data,'.b')
            plot(szRng,s.fitout1,'-r')
            plot(szRng,s.fitout2,'-g')
            hold off
            ylim([min([-0.5*s.maxResp1 min(s.data)]) 1.2*max([s.maxResp2 max(s.data)])])
            title(['Cell #' num2str(iCell) ': Dist= ' num2str(cellDists(iCell,1)) '- Size Tuning for high Con (Ftest=' num2str(fit_true_vec(iCell,8)) ')']);
            xlabel('Stimulus size (deg)')
            ylabel('dF/F')
        end
        
        for count2 = 1:Npars
            tmp = squeeze(fit_shuf_vec(iCell,count2,:,nCon));
            [i,j] = sort(tmp); % sort in order
            lbub_fits(iCell,count2,1) = i(ind_shuf_lb); %lower 0.025
            lbub_fits(iCell,count2,2) = i(ind_shuf_ub); %upper 0.975
            lbub_fits(iCell,count2,3) = mean(i); %mean
            lbub_fits(iCell,count2,5) = std(i); %stdev
            
            if sum(iCell==chosen)
                switch count2
                    case 1
                        subplot(3,3,4)
                        histogram(i,0:max(szs)+1)
                        xlim([0 100]);
                        line([0.5*mean(i) 0.5*mean(i)], ylim, 'color','red')
                        line([2*mean(i) 2*mean(i)], ylim, 'color','red')
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        title(['PrefSize shuffles cell #' num2str(iCell)])
                        xlabel('Size (deg)')
                        ylabel('count')
                    case 2
                        subplot(3,3,5)
                        histogram(i,0:max(szs)+1)
                        xlim([0 100]);
                        line([0.5*mean(i) 0.5*mean(i)], ylim, 'color','red')
                        line([2*mean(i) 2*mean(i)], ylim, 'color','red')
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        title(['PrefSize1 shuffles cell #' num2str(iCell)])
                        xlabel('Size (deg)')
                    case 3
                        subplot(3,3,6)
                        histogram(i,0:max(szs)+1)
                        xlim([0 100]);
                        line([0.5*mean(i) 0.5*mean(i)], ylim, 'color','red')
                        line([2*mean(i) 2*mean(i)], ylim, 'color','red')
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        title(['PrefSize2 shuffles cell #' num2str(iCell)])
                        xlabel('Size (deg)')
                    case 4
                        subplot(3,3,7)
                        histogram(i,0:0.01:2)
                        xlim([0 1])
                        line([mean(i) mean(i)], ylim)
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        title(['SuppInd shuffles cell #' num2str(iCell)])
                        xlabel('SI')
                        ylabel('count')
                    case 5
                        subplot(3,3,8)
                        histogram(i,0:0.01:2)
                        xlim([0 1])
                        line([mean(i) mean(i)], ylim)
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        title(['SuppInd1 shuffles cell #' num2str(iCell)])
                        xlabel('SI1')
                    case 6
                        subplot(3,3,9)
                        histogram(i,0:0.01:2)
                        xlim([0 1])
                        line([mean(i) mean(i)], ylim)
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        title(['SuppInd2 shuffles cell #' num2str(iCell)])
                        xlabel('SI2')
                    case 7
                        Fcrit = finv(0.95,nPts-6,nPts-3);
                        subplot(3,3,2)
                        histogram(i,50)
                        xlim([0 max(i)]);
                        line([mean(i) mean(i)], ylim)
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        line([Fcrit Fcrit], ylim, 'LineWidth', 3, 'color', 'red')
                        title(['Fscore shuffles cell #' num2str(iCell)])
                        xlabel('Fscore')
                        ylabel('count')
                    case 8
                        subplot(3,3,3)
                        histogram(i,[-0.5 0.5 1.5])
                        xlim([-0.5 1.5])
                        ylim([0 Nshuf])
                        title(['Ftest shuffles cell #' num2str(iCell)])
                        xlabel('Ftest')
                        pause
                end
            end
        end
        lbub_fits(iCell,:,4) = fit_true_vec(iCell,:,nCon); % true (no shuffle)
    end

%% determine good fits
% first sort cells based on Ftest same as True fit in >50% of shuffles
% then use that model's prefSize confidence interval (e.g. prefSize1 or 2)
% check bounds of confidence interval are within 1 octave of "True" fit
% e.g. lower(2.5th)>0.5*prefSizeTrue and upper(97.5th)<2*prefSizeTrue
goodfit_ind_size = [];
for iCell = 1:nCells
    Ftest = lbub_fits(iCell,8,4); %8=Ftest, 4=True fit
    Ftestshuf = lbub_fits(iCell,8,3); %8=Ftest, 3=mean of shuffles
    lOct = 0.5*lbub_fits(iCell,1,4); %1=prefSize, 4=True fit, lower octave
    hOct = 2*lbub_fits(iCell,1,4); %1=prefSize, 4=True fit, upper octave
    switch Ftest
        case 0 % model1
            if Ftestshuf<0.5 %model1 >50% of shuffles
                if (lbub_fits(iCell,2,1)>lOct) && (lbub_fits(iCell,2,2)<hOct) %2=prefSize1, 1=lb/2=ub
                    goodfit_ind_size = [goodfit_ind_size iCell];
                end
            end
        case 1 % model2
            if Ftestshuf>0.5 % model2 >50% of shuffles
                if (lbub_fits(iCell,3,1)>lOct) && (lbub_fits(iCell,3,2)<hOct) %3=prefSize2, 1=lb/2=ub
                    goodfit_ind_size = [goodfit_ind_size iCell];
                end
            end
    end
end

% is model1 + is model2
ism1 = intersect(goodfit_ind_size, find(~lbub_fits(:,8,4)));
ism2 = intersect(goodfit_ind_size, find(lbub_fits(:,8,4)));

fprintf(['#Good cells = ' num2str(length(goodfit_ind_size)) '\nModel 1: ' num2str(length(ism1)) ', Model 2: ' num2str(length(ism2)) '\nSaving good fits\n'])
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_lbub_fits']), 'goodfit_ind_size', 'lbub_fits')

%% plot summary
tuning_norm = squeeze(tuning_mat(:,:,1,:)./max(max(tuning_mat(:,:,1,:),[],1),[],2));
figure;

subplot(2,2,1)
ind = intersect(ism1, find(cellDists<=10));
ret_avg = mean(tuning_norm(:,ind),2)';
ret_sem= std(tuning_norm(:,ind),[],2)./sqrt(length(ind))';
errorbar(szs,ret_avg,ret_sem)
hold on
ylabel('dF/F normalized')
xlabel('Size (deg)')
ylim([-0.3 1.2])
axis square
title(['M1 cells within 10 deg of stim- n = ' num2str(length(ind))])
subplot(2,2,2)
ind = intersect(ism2, find(cellDists<=10));
ret_avg = mean(tuning_norm(:,ind),2)';
ret_sem= std(tuning_norm(:,ind),[],2)./sqrt(length(ind))';
errorbar(szs,ret_avg,ret_sem)
hold on
ylabel('dF/F normalized')
xlabel('Size (deg)')
ylim([-0.3 1.2])
axis square
title(['M2 cells within 10 deg- n = ' num2str(length(ind))])
subplot(2,2,3)
ind = intersect(goodfit_ind_size, find(cellDists<=10));
ret_avg = mean(tuning_norm(:,ind),2)';
ret_sem= std(tuning_norm(:,ind),[],2)./sqrt(length(ind))';
errorbar(szs,ret_avg,ret_sem)
hold on
ylabel('dF/F normalized')
xlabel('Size (deg)')
ylim([-0.3 1.2])
axis square
title(['All well-fit cells within 10 deg- n = ' num2str(length(ind))])

subplot(2,2,4)
hist(lbub_fits(ind,1,4))
xlim([0 80])
axis square
xlabel('Pref size (deg)')
ylabel('# neurons')
expLoc = expt(iexp).img_loc{1};
expStrct = cell2mat(expt(iexp).img_strct);
suptitle([mouse ' ' date ' ' expLoc ' ' expStrct])

print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_avgTuning_goodfits.pdf']), '-dpdf','-bestfit')