%% Choose experiment
clc; clear all; close all;

doRedChannel = 0;

mouse = 'i2707';
date = '230808';
runs = {'002'};
time = {'1709'};
nrun = length(runs);
frame_rate = 15;
ImgFolder = runs;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);
fprintf([mouse ' ' date '\n'])
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
fnout = fullfile(base, '\Analysis\2P\', [date '_' mouse], [date '_' mouse '_' run_str]);


%% Load temp data

load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataTEMP.mat']))

% load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
% load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
%%

% data_f = OFF
% data_resp = ON

b = 5;
mask_cell = zeros(size(data_f,1), size(data_f,2));

for i = 2:(size(data_f,1)-1)
    for j = 2:(size(data_f,2)-1)
        for is = 1:nStimCon
            ind_stim = find(stimCon_all == stimCons(is));
            for im = 1:nMaskCon
                ind_mask = find(maskCon_all == maskCons(im));
                if im>1 & is>1
                    for ip = 1:nMaskPhas
                        ind_p = find(maskPhas_all == maskPhas(ip));
                        ind_all = intersect(ind_p,intersect(ind_stim,ind_mask));
                        data_on = mean(reshape(data_f(i-1:i+1,j-1:j+1,ind_all),9,[]));
                        data_off = mean(reshape(data_resp(i-1:i+1,j-1:j+1,ind_all),9,[]));
                        [h(is,im,ip), p] = ttest(data_on,data_off,'Alpha',0.01);
                    end
                else
                    ind_all = intersect(ind_stim,ind_mask);
                    data_on = mean(reshape(data_f(i-1:i+1,j-1:j+1,ind_all),9,[]));
                    data_off = mean(reshape(data_resp(i-1:i+1,j-1:j+1,ind_all),9,[]));
                    [h(is,im,1), p] = ttest(data_on,data_off,'Alpha',0.01);
                end
            end
        end
        idx = find(h==1);
        if idx>5
            mask_cell(i,j) = 1;
        end
    end
end


figure; imagesc(mask_cell)
print(fullfile(fnout, [date '_' mouse '_' run_str  '_bouton_mask_resp6.pdf']), '-dpdf')


min_df = 0.1;
sz = size(mask_cell);
temp_data = ones(sz(1),sz(2));
temp_max = squeeze(max(max(data_resp_dfof,[],1),[],2));
[a, max_sort] = sort(temp_max,'descend');
data_dfof_avg_sort = cat(3,data_dfof_max,data_dfof(:,:,max_sort));



for ia = 1:size(data_dfof_avg_sort,3)
    fprintf([ '\n img:' num2str(ia)])
    temp_data_log = ~isnan(temp_data);
    [x, y, v] = find(data_dfof_avg_sort(:,:,ia).*temp_data_log);
    [a, ind_sort] = sort(v,'descend');
    ind = find(a>min_df);
    fprintf(['- ' num2str(length(ind)) ' pix: '])
    for a = 1:length(ind)
        i = x(ind_sort(a));
        j = y(ind_sort(a));
        if i>5 & j>5 & i<sz(1)-b & j<sz(2)-b
            if ~isnan(temp_data(i-1:i+1,j-1:j+1))
                all_pix = data_dfof_avg_sort(i-1:i+1,j-1:j+1,ia);
                [max_val max_ind] = max(all_pix(:));
                if max_ind == 5
                    h = zeros(length(stimCons),length(maskCons),length(maskPhas));
                    for is = 1:nStimCon
                        ind_stim = find(stimCon_all == stimCons(is));
                        for im = 1:nMaskCon
                            ind_mask = find(maskCon_all == maskCons(im));
                            if im>1 & is>1
                                for ip = 1:nMaskPhas
                                    ind_p = find(maskPhas_all == maskPhas(ip));
                                    ind_all = intersect(ind_p,intersect(ind_stim,ind_mask));
                                    data_on = mean(reshape(data_f(i-1:i+1,j-1:j+1,ind_all),9,[]));
                                    data_off = mean(reshape(data_resp(i-1:i+1,j-1:j+1,ind_all),9,[]));
                                    [h(is,im,ip), p] = ttest(data_on,data_off,'Alpha',0.01);
%                                     [h(is,im,ip) p] = ttest(squeeze(mean(mean(data_resp_dfeduof(i-1:i+1,j-1:j+1,ind_all),1),2)),0,'tail','right','alpha',0.001./(size(data_dfof,3)));
                                end
                            else 
                                ind_all = intersect(ind_stim,ind_mask);
                                [h(is,im,1) p] = ttest(squeeze(mean(mean(data_resp_dfof(i-1:i+1,j-1:j+1,ind_all),1),2)),0,'tail','right','alpha',0.001./(size(data_dfof,3)));
                            end
                        end
                    end
                    if sum(h(:))==1
                        mask_cell_hOne(i,j) = 1;
                    end
                    if sum(h(:))>2
                        mask_cell(i,j) = 1;
                        mask_all(i-1:i+1,j-1:j+1) = ones(3,3);
                        temp_data(i-2:i+2,j-2:j+2) = NaN(5,5);
                        fprintf('.')
                    end
                    mask_cell_bright(i,j) = 1;
                end
            end
        end
    end
end      
   





%%
% COPYING AND PASTING ORIGINAL SCRIPT
% Changing it so that instead of choosing the brightest of 8 pixels, it
% chooses the most significant of the 8 pixels...
%
%
%
%
%
%
%
%
%
%
%%
% retOnly script by kevin 
%get path names
close all;clear all;clc;

mouse = 'i2707';
date = '230808';
runs = {'002'};
time = {'1709'};
nrun = length(runs);
frame_rate = 15;
ImgFolder = runs;

run_str = catRunName(runs, 1);

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\';
SG_base = fullfile(fn_base, '\home\sara');
data_fn = fullfile(SG_base, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data');
% fnout = fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '-' run_str]);  %   why is this turning into a cell array?????
base = ('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\sara\Analysis\2P\');
fnout = fullfile(SG_base, '\Analysis\2P\', [date '_' mouse], [date '_' mouse '_' run_str]);

if ~exist(fnout)
    mkdir(fnout)
end

fprintf(['2P imaging retinotopy analysis - by KM, Glickfeld Lab\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([runs{irun} ' - ' time{irun} '\n'])
end

%% load data, read with sbxread, and concatenate selected runs
data = [];
clear temp
trial_n = zeros(1,nrun);

fprintf(['\nBegin reading ' num2str(nrun) ' runs...'])
for irun = 1:nrun
    %load imaging data
    dataFolder = runs{irun};
    fName = [dataFolder '_000_000'];

    % load behavior data
    mwName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data\data-' mouse '-' date '-' time{irun} '.mat'];
    load(mwName);

    temp(irun) = input;
    nframes = input.counterValues{end}(end);

    cd(fullfile(data_fn,mouse,date,runs{irun}))
    imgMatFile = [fName '.mat'];
    load(imgMatFile);
    data_temp = sbxread(fName,0,nframes);
    data_temp = squeeze(data_temp);

    data = cat(3,data,data_temp);
    fprintf('Complete\n')
end
fprintf('All runs read\n')
fprintf([num2str(size(data,3)) ' total frames\n'])
input = concatenateDataBlocks(temp);
clear data_temp
clear temp
%% Choose register interval
regIntv = 3000; %3000 for axons
nep = floor(size(data,3)./regIntv);
    % plot 500 frame means at each register interval
[n, n2] = subplotn(nep);
figure(1);clf;
% colormap(gray)
for i = 1:nep
    subplot(n,n2,i);
    imagesc(mean(data(:,:,(1:500)+(i-1)*regIntv),3));
    title([num2str(i) ': ' num2str(1+((i-1)*regIntv)) '-' num2str(500+((i-1)*regIntv))]);
%     clim([600 1500])
    clim([100 600])
end
movegui('center')
%% Register data

chooseInt = 8; %nep/2 % interval chosen for data_avg =[epoch of choice]-1

fprintf('\nBegin registering...\n')
    
    meanrng = regIntv*(chooseInt)+(1:500);
    data_avg = mean(data(:,:,meanrng),3);
    fprintf(['\nRegister frame averaged from ' num2str(meanrng(1)) ' - ' num2str(meanrng(end)) '\n'])
    
    %check if data has already been registered
        if exist(fullfile(fnout, [date '_' mouse '_' run_str '_reg_shifts.mat']))
            fprintf('load registered data\n')
            load(fullfile(fnout, [date '_' mouse '_' run_str '_reg_shifts.mat']))
            [outs, data_reg] = stackRegister_MA(data,[],[],out);
        else %register
            fprintf('stackRegister\n')
            [out, data_reg] = stackRegister(data,data_avg);
            save(fullfile(fnout, [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
            save(fullfile(fnout, [date '_' mouse '_' run_str '_input.mat']), 'input')
        end

    % save
    fprintf('Registration complete, now saving...\n')
    save(fullfile(fnout, [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg','meanrng')
    save(fullfile(fnout, [date '_' mouse '_' run_str '_input.mat']), 'input')
%end
clear data % depending on memory

stop

%% subsample data_reg to be a smaller window for analysis

% data_reg = data_reg(200:300,550:650,:);

%end subsample
%% test stability

% figure 2 shows the registered images to check the stability
fprintf('\nExamine registered images for stability\n')
figure(2);clf;
[n,n2] = optimizeSubplotDim(nep);
for i = 1:nep
    subplot(n,n2,i);
    imagesc(mean(data_reg(:,:,(1:500)+((i-1)*regIntv)),3));
    title([num2str(i) ': ' num2str(1+((i-1)*regIntv)) '-' num2str(500+((i-1)*regIntv))]);
end

sz = size(data_reg);

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
print(fullfile(fnout, [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf','-fillpage')

%% find activated cells
% calculate dF/F
fprintf('\nBegin image analysis...\n')

ntrials = input.trialsSinceReset;

if isfield(input, 'cStimOneOn') 
    cStimOn = celleqel2mat_padded(input.cStimOneOn);
    nTrials = length(cStimOn);
    sz = size(data_reg);

    data_resp = nan(sz(1),sz(2),nTrials);
    data_f = nan(sz(1),sz(2),nTrials);

    for itrial = 1:nTrials
        if cStimOn(itrial) + 20 < sz(3)
            data_resp(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)+5:cStimOn(itrial)+20),3);
        end
        data_f(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)-15:cStimOn(itrial)),3);
        data_resp_dfof = (data_resp-data_f)./data_f;
    end
else
        if isfield(input, 'nScansOn') 
          nOn = input.nScansOn; 
          nOff = input.nScansOff; 
          if nOn>29
            sz = size(data_reg);
            data_tr = reshape(data_reg,[sz(1), sz(2), nOn+nOff, ntrials]); %trial level data - now a 4d structure
            data_tr = data_tr(:,:,nOff/2:end,:);
            data_f = mean(data_tr(:,:,1:nOff/2,:),3); %avg across last half of OFF frames (baseline)
            data_df = bsxfun(@minus, double(data_tr), data_f); 
            data_resp_dfof = bsxfun(@rdivide,data_df, squeeze(data_f)); 
            clear data_tr data_f data_df
          else 
          end
        else
            stop
        end
end
        

stimCon_all = celleqel2mat_padded(input.tStimOneGratingContrast);
maskCon_all = celleqel2mat_padded(input.tMaskOneGratingContrast);
stimCons = unique(stimCon_all);
maskCons = unique(maskCon_all);
nStimCon = length(stimCons);
nMaskCon = length(maskCons);
maskPhas_all = celleqel2mat_padded(input.tMaskOneGratingPhaseDeg);
maskPhas = unique(maskPhas_all);
nMaskPhas = length(maskPhas);
stimDir_all = celleqel2mat_padded(input.tStimOneGratingDirectionDeg);
maskDir_all = rad2deg(wrapTo2Pi(deg2rad(celleqel2mat_padded(input.tMaskOneGratingDirectionDeg))));
maskDir_all(find(maskDir_all==360)) = 0;
stimDirs = unique(stimDir_all);
nStimDir = length(stimDirs);
maskDirs = unique(maskDir_all);
nMaskDir = length(maskDirs);
maskDiff_all = celleqel2mat_padded(input.tMaskOneGratingDirectionDeg) - celleqel2mat_padded(input.tStimOneGratingDirectionDeg);
maskDiffs = unique(maskDiff_all);
nMaskDiff = length(maskDiffs);
SF_all = celleqel2mat_padded(input.tStimOneGratingSpatialFreqCPD);
SFs = unique(SF_all);
nSF = length(SFs);
TF_all = celleqel2mat_padded(input.tMaskOneGratingTemporalFreqCPS);
TFs = unique(TF_all);
nTF = length(TFs);

if nStimCon >= 2  || nMaskCon >=2 
    nStim1 = nStimCon*nMaskCon*nMaskPhas;
    data_dfof = nan(sz(1),sz(2),nStim1);
    Stims = zeros(nStim1,3);
    start = 1;
    for is = 1:nStimCon
        ind_stim = find(stimCon_all == stimCons(is));
        for im = 1:nMaskCon
            ind_mask = find(maskCon_all == maskCons(im));
            if im>1 & is>1
                for ip = 1:nMaskPhas
                    ind_p = find(maskPhas_all == maskPhas(ip));
                    ind_all = intersect(ind_p,intersect(ind_stim,ind_mask));
                    data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,ind_all),3);
                    Stims(start,:) = [stimCons(is) maskCons(im) maskPhas(ip)];
                    start = start+1;
                end
            elseif im>1 & is==1 
                ind_all = intersect(ind_stim,ind_mask);
                data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,ind_all),3);
                Stims(start,:) = [stimCons(is) maskCons(im) 0];
                start = start+1;
            elseif im==1 & is>1 
                ind_all = intersect(ind_stim,ind_mask);
                data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,ind_all),3);
                Stims(start,:) = [stimCons(is) maskCons(im) 0];
                start = start+1;
            end
        end
    end
else 
    nStim1 = 0;
    data_dfof = [];
end

data_dfof(:,:,isnan(mean(mean(data_dfof,1),2))) = [];
myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_max = max(imfilter(data_dfof,myfilter),[],3);
figure;
movegui('center')
imagesc(data_dfof_max)
print(fullfile(fnout, [date '_' mouse '_' run_str '_maxdfof.pdf']),'-dpdf')

save(fullfile(fnout, [date '_' mouse '_' run_str '_dataStim.mat']), 'cStimOn', 'maskCon_all', 'stimCon_all', 'stimCons', 'maskCons', 'nStimCon', 'nMaskCon', 'stimDir_all', 'stimDirs', 'nStimDir', 'maskDir_all', 'maskDirs', 'nMaskDir', 'maskPhas_all', 'maskPhas', 'nMaskPhas', 'maskDiff_all','maskDiffs','nMaskDiff','SF_all', 'SFs', 'nSF','TF_all', 'TFs', 'nTF', 'frame_rate', 'nTrials')

%% axon segmentation
% find hotspots

fprintf('\nBegin axon segmentation...')
min_df = 0.2; %CHANGED from 0.05 to 0
b = 5;
mask_cell_bright = zeros(sz(1), sz(2));
mask_cell_sig = zeros(sz(1), sz(2));
mask_cell = zeros(sz(1), sz(2));
mask_all = zeros(sz(1), sz(2));
temp_max = squeeze(max(max(data_dfof,[],1),[],2));
[a, max_sort] = sort(temp_max,'descend');
data_dfof_avg_sort = cat(3,data_dfof_max,data_dfof(:,:,max_sort));

for ia = 1:size(data_dfof_avg_sort,3)
    fprintf([ '\n img:' num2str(ia)])
    [x, y, v] = find(data_dfof_avg_sort(:,:,ia));
    [a, ind_sort] = sort(v,'descend');
    ind = find(a>min_df);
    fprintf(['- ' num2str(length(ind)) ' pix: '])
    for a = 1:length(ind)
        i = x(ind_sort(a));
        j = y(ind_sort(a));
        mask_cell_bright(i,j) = mask_cell_bright(i,j) + 1;
    end
end   

for i = 1:size(mask_cell_bright,1)
    for j = 1:size(mask_cell_bright,2)
        if mask_cell_bright(i,j) <2
            mask_cell_bright(i,j) = 0;
        end
    end
end

figure;
imagesc(mask_cell_bright)


for i = 2:(size(data_f,1)-1)
    for j = 2:(size(data_f,2)-1)
        pix = mask_cell_bright(i,j);
        if pix>0
            for is = 1:nStimCon
                ind_stim = find(stimCon_all == stimCons(is));
                for im = 1:nMaskCon
                    ind_mask = find(maskCon_all == maskCons(im));
                    if im>1 & is>1
                        for ip = 1:nMaskPhas
                            ind_p = find(maskPhas_all == maskPhas(ip));
                            ind_all = intersect(ind_p,intersect(ind_stim,ind_mask));
                            data_on = mean(reshape(data_f(i-1:i+1,j-1:j+1,ind_all),9,[]));
                            data_off = mean(reshape(data_resp(i-1:i+1,j-1:j+1,ind_all),9,[]));
                            [h(is,im,ip), p] = ttest(data_on,data_off,'Alpha',0.001./(size(data_dfof,3)));
                        end
                    else
                        ind_all = intersect(ind_stim,ind_mask);
                        data_on = mean(reshape(data_f(i-1:i+1,j-1:j+1,ind_all),9,[]));
                        data_off = mean(reshape(data_resp(i-1:i+1,j-1:j+1,ind_all),9,[]));
                        [h(is,im,1), p] = ttest(data_on,data_off,'Alpha',0.001./(size(data_dfof,3)));
                    end
                end
            end
            idx = find(h==1);
            if idx>1
                mask_cell_sig(i,j) = 1;
            else
            end
        else
        end
    end
end
figure;
imagesc(mask_cell_sig)


temp_data = mask_cell_sig;
temp_data(temp_data==0) = NaN;

for ia = 1:size(data_dfof_avg_sort,3)
    fprintf([ '\n img:' num2str(ia)])
    temp_data_log = ~isnan(temp_data);
    [x, y, v] = find(data_dfof_avg_sort(:,:,ia).*temp_data_log);
    [a, ind_sort] = sort(v,'descend');
    ind = find(a>min_df);
    fprintf(['- ' num2str(length(ind)) ' pix: '])
    for a = 1:length(ind)
        i = x(ind_sort(a));
        j = y(ind_sort(a));
        if i>5 & j>5 & i<sz(1)-b & j<sz(2)-b
            all_pix = data_dfof_avg_sort(i-1:i+1,j-1:j+1,ia);
            [max_val max_ind] = max(all_pix(:));
            if max_ind == 5
                mask_cell(i,j) = 1;
                mask_all(i-1:i+1,j-1:j+1) = ones(3,3);
                temp_data(i-2:i+2,j-2:j+2) = NaN(5,5);
                fprintf('.')
            end
        end
    end
end      

figure; imagesc(mask_cell)
print(fullfile(fnout, [date '_' mouse '_' run_str  '_bouton_mask.pdf']), '-dpdf')


shadeimg = imShade(data_dfof_max,mask_cell); figure; imagesc(shadeimg)

save(fullfile(fnout, [date '_' mouse '_' run_str '_mask_cell.mat']), 'mask_cell', 'mask_all')
save(fullfile(fnout, [date '_' mouse '_' run_str '_FOVs.mat']), 'data_dfof','data_dfof_max','Stims')

% look at masks overlaid on dfof average FOV for each stimulus
figure;
start = 1;
for i = 1:size(data_dfof_avg_sort,3)
    subplot(3,3,start);
    shade_img = imShade(data_dfof_avg_sort(:,:,i), mask_all);
    imagesc(shade_img)
    if i==1
        title('max')
    else
        title(num2str(Stims(i-1,:)))
    end
    clim([0 max(data_dfof(:))]);
    colormap(gray)
    if start<9
        start = start+1;
    else
        figure;
        start = 1;
    end
end
set(gcf, 'Position', [0 0 800 1000]);
%print(fullfile(fnout, dataFolder, [mouse '_' expDate '_FOVresp.pdf']), '-dpdf')

%% Get time courses
nCells = sum(mask_cell(:));
fprintf(['Found ' num2str(nCells) ' boutons\n'])

fprintf('Extracting cell signal...\n')
data_tc = zeros(sz(3),nCells);
iC = 1;
for i = 1:sz(1)
    ind = find(mask_cell(i,:));
    if length(ind)>0
        for ii = 1:length(ind)
            fprintf([num2str(iC) ' '])
            j = ind(ii);
            data_tc(:,iC) = squeeze(mean(mean(data_reg(i-1:i+1,j-1:j+1,:),1),2));
            iC = 1+iC;
        end
    end
end
fprintf([num2str(nCells) ' total cells extracted\n'])
s = whos('data_tc');
if s.bytes < 2000000000
    save(fullfile(fnout, [date '_' mouse '_' run_str  '_TCs.mat']), 'data_tc', 'sz', 'nCells')
else
    save(fullfile(fnout, [date '_' mouse '_' run_str  '_TCs.mat']), 'data_tc' ,'-v7.3', 'sz', 'nCells')
end

fprintf('Time course extraction complete.\n')
clear data_reg
