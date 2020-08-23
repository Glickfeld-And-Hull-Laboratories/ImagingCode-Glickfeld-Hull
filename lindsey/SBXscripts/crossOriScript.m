clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandPhase_ExptList';
iexp = 5; 
rc = behavConstsAV;
eval(ds)

frame_rate = 30;

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};
ImgFolder = expt(iexp).coFolder;
time = expt(iexp).coTime;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
%LG_base = '\\CRASH.dhe.duke.edu\data\home\lindsey';

fprintf(['2P imaging retinotopy analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([ImgFolder{irun} ' - ' time{irun} '\n'])
end

%% load and register
tic
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = [LG_base '\Data\2P_images\' mouse '\' date '\' ImgFolder{irun}];
    %CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\AW68\two-photon imaging\' date '\' ImgFolder(irun,:)];
    %CD = [LG_base '\Data\2P_images\' mouse '-KM' '\' date '_' mouse '\' ImgFolder(irun,:)];
    %CD = ['\\CRASH.dhe.duke.edu\data\home\kevin\Data\2P\' date '_' mouse '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder{irun} '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time{irun} '.mat'];
    load(fName);
    
    temp(irun) = input;
    nframes = [temp(irun).counterValues{end}(end) info.config.frames];
    
    fprintf(['Reading run ' num2str(irun) '- ' num2str(min(nframes)) ' frames \r\n'])
    data_temp = sbxread(imgMatFile(1,1:11),0,min(nframes));
    if size(data_temp,1)== 2
        data_temp = data_temp(1,:,:,:);
    end
    
    if isfield(input, 'cStimOneOn') 
        if irun>1
            ntrials = size(input.cStimOneOn,2);
            for itrial = 1:ntrials
                temp(irun).cStimOneOn{itrial} = temp(irun).cStimOneOn{itrial}+offset;
                temp(irun).cStimOneOff{itrial} = temp(irun).cStimOneOff{itrial}+offset;
            end
        end
    end
    offset = offset+min(nframes);
        
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    trial_n = [trial_n nframes];
end
input = concatenateStructuresLG(temp);
clear data_temp
clear temp
toc

%% Choose register interval
nep = floor(size(data,3)./10000);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); colormap gray; clim([0 3000]); end

data_avg = mean(data(:,:,50001:50500),3);
%% Register data
if exist(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    [outs, data_reg] = stackRegister_MA(data,[],[],out);
else
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
clear data out

%% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_byFrame.pdf']),'-dpdf', '-bestfit')

figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

i = 1;
sz = size(data_reg);
rg = zeros(sz(1),sz(2),3);
first = mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3);
rg(:,:,1) = first./max(first(:));
i = nep; 
last = mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3);
rg(:,:,2) = last./max(last(:));
figure; image(rg)
%% if red channel data
if doRedChannel
    CD = [LG_base '\Data\2P_images\' date '_' mouse '\' ImgFolderRed(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolderRed(irun,:) '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
    load(fName);
    nframes = info.config.frames;
    fprintf(['Reading red data- ' num2str(nframes) ' frames \r\n'])
    data_red = sbxread(imgMatFile(1,1:11),0,nframes);
    data_red_g = squeeze(data_red(1,:,:,:));
    data_red_r = squeeze(data_red(2,:,:,:));

    [rg_out rg_reg] = stackRegister(data_red_g,data_avg);
    [rr_out rr_reg] = stackRegister_MA(data_red_r,[],[],rg_out);

    data_red_avg = mean(rr_reg,3);
    figure; imagesc(data_red_avg);
end

%% find activated cells

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
end

data_resp_dfof = (data_resp-data_f)./data_f;

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
maskDir_all = celleqel2mat_padded(input.tMaskOneGratingDirectionDeg);
stimDirs = unique(stimDir_all);
nStimDir = length(stimDirs);
maskDirs = unique(maskDir_all);
nMaskDir = length(maskDirs);
SF_all = celleqel2mat_padded(input.tStimOneGratingSpatialFreqCPD);
SFs = unique(SF_all);
nSF = length(SFs);

if nStimCon >= 2  || nMaskCon >=2 
    nStim1 = nStimCon*nMaskCon*nMaskPhas;
    data_dfof = nan(sz(1),sz(2),nStim1);

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
                    start = start+1;
                end
            else
                ind_all = intersect(ind_stim,ind_mask);
                data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,ind_all),3);
                start = start+1;
            end
        end
    end
else 
    nStim1 = 0;
    data_dfof = [];
end


if nStimDir > 1
    nStim2 = nStimDir.*2;
    data_dfof = cat(3, data_dfof, zeros(sz(1),sz(2), nStim2));
    start = nStim1+1;
    for is = 1:nStimDir
        ind_stimalone = intersect(intersect(find(stimCon_all),find(maskCon_all==0)),find(stimDir_all == stimDirs(is)));
        ind_maskalone = intersect(intersect(find(stimCon_all==0),find(maskCon_all)),find(maskDir_all == stimDirs(is)));
        data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,[ind_stimalone ind_maskalone]),3);
        ind_plaidstim = intersect(intersect(find(stimCon_all),find(maskCon_all)),find(stimDir_all == stimDirs(is)));
        ind_plaidmask = intersect(intersect(find(stimCon_all),find(maskCon_all)),find(maskDir_all == stimDirs(is)));
        data_dfof(:,:,start+1) = nanmean(data_resp_dfof(:,:,[ind_plaidstim ind_plaidmask]),3);
        start = start+2;
    end
else
    nStim2 = 0;
end

if nSF>1
    nStim3 = nSF;
    data_dfof = cat(3, data_dfof, zeros(sz(1),sz(2), nStim3));
    start = nStim1+nStim2+1;
    for iSF = 1:nSF
        ind_SF = find(SF_all == SFs(iSF));
        data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,ind_SF),3);
        start = start+1;
    end
else
    nStim3 = 0;
end

data_dfof(:,:,isnan(mean(mean(data_dfof,1),2))) = [];
myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_max = max(imfilter(data_dfof,myfilter),[],3);
figure;
movegui('center')
imagesc(data_dfof_max)
data_dfof = cat(3, data_dfof, data_dfof_max);
if doRedChannel
    data_dfof = cat(3,data_dfof,data_red_avg);
end

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']), 'cStimOn', 'maskCon_all', 'stimCon_all', 'stimCons', 'maskCons', 'nStimCon', 'nMaskCon', 'stimDir_all', 'stimDirs', 'nStimDir', 'maskDir_all', 'maskDirs', 'nMaskDir', 'SF_all', 'SFs', 'nSF', 'frame_rate', 'nTrials')

%% cell segmentation 
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
mask_data = data_dfof;

for iStim = 1:size(data_dfof,3)   
    mask_data_temp = mask_data(:,:,end+1-iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0;
    bwout = imCellEditInteractiveLG(mask_data_temp);
    if doRedChannel & iStim==1
        red_mask = bwout;
    end
    mask_all = mask_all+bwout;
    mask_exp = imCellBuffer(mask_all,3)+mask_all;
    close all
end
mask_cell= bwlabel(mask_all);
if doRedChannel
    red_cells = unique(mask_cell(find(red_mask)));
else
    red_cells = [];
end
figure; movegui('center')
imagesc(mask_cell)

clear data_adapt data_adapt_dfof data_test data_test_dfof data_test_avg
%% neuropil subtraction
mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof', 'mask_cell', 'mask_np', 'red_cells')

clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_data_temp mask_exp data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 


down = 5;
sz = size(data_reg);

data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_reg_down  = stackGroupProject(data_reg,down);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
nCells = size(data_tc,2);
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
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
clear data_reg data_reg_down

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc', 'nCells', 'sz')

clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell

