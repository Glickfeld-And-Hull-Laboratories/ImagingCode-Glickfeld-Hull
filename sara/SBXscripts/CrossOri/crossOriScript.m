clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandDirFourPhase_ExptList_SG';
iexp = 122;  
doPhaseAfterDir = 0;
doDirAfterPass = 0;
eval(ds)

frame_rate = params.frameRate;

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};
if doPhaseAfterDir
    ImgFolder = expt(iexp).coFolder;
    nrun = length(ImgFolder);
    ref_str = catRunName(cell2mat(ImgFolder), nrun);
    ImgFolder = expt(iexp).copFolder;
    time = expt(iexp).copTime;
elseif doDirAfterPass
    ImgFolder = expt(iexp).passFolder;
    nrun = length(ImgFolder);
    ref_str = catRunName(cell2mat(ImgFolder), nrun);
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
else
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
end
nrun = length(ImgFolder);
% run_str = catRunName(cell2mat(ImgFolder), nrun);
run_str = 'runs-002';

base = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\' expt(iexp).saveLoc];

fprintf(['2P imaging cross ori analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
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
    CD = [base '\Data\2P_images\' mouse '\' date '\' ImgFolder{irun}];
    %CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\AW68\two-photon imaging\' date '\' ImgFolder(irun,:)];
    %CD = [base '\Data\2P_images\' mouse '-KM' '\' date '_' mouse '\' ImgFolder(irun,:)];
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

% Choose register interval
regIntv = 5000;
nep = floor(size(data,3)./regIntv);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*regIntv):500+((i-1)*regIntv)),3)); title([num2str(1+((i-1)*regIntv)) '-' num2str(500+((i-1)*regIntv))]); colormap gray; clim([0 3000]); end
movegui('center')
%% Register data
data_avg = mean(data(:,:,85001:85500),3);
if doPhaseAfterDir || doDirAfterPass
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_reg_shifts.mat']))
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
elseif exist(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    [outs, data_reg] = stackRegister_MA(data,[],[],out);
else
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end

% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*regIntv):500+((i-1)*regIntv)),3)); title([num2str(1+((i-1)*regIntv)) '-' num2str(500+((i-1)*regIntv))]); end
print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_byFrame.pdf']),'-dpdf', '-bestfit')
movegui('center')
figure; imagesq(mean(data_reg(:,:,1:regIntv),3)); truesize;
print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

i = 1;
sz = size(data_reg);
rg = zeros(sz(1),sz(2),3);
first = mean(data_reg(:,:,1+((i-1)*regIntv):500+((i-1)*regIntv)),3);
rg(:,:,1) = first./max(first(:));
i = nep; 
last = mean(data_reg(:,:,1+((i-1)*regIntv):500+((i-1)*regIntv)),3);
rg(:,:,2) = last./max(last(:));
figure; image(rg)
movegui('center')
%% if red channel data
if doRedChannel & (~doPhaseAfterDir & ~doDirAfterPass)
    ImgFolderRed = expt(iexp).redImg;
    CD = [base '\Data\2P_images\' mouse '\' date '\' ImgFolderRed{1}];
    cd(CD);
    imgMatFile = [ImgFolderRed{1} '_000_000.mat'];
    load(imgMatFile);

    nframes = info.config.frames;
    fprintf(['Reading red data- ' num2str(nframes) ' frames \r\n'])
    data_red = sbxread(imgMatFile(1,1:11),0,nframes);
    data_red_g = squeeze(data_red(1,:,:,:));
    data_red_r = squeeze(data_red(2,:,:,:));

    [rg_out rg_reg] = stackRegister(data_red_g,data_avg);
    [rr_out rr_reg] = stackRegister_MA(data_red_r,[],[],rg_out);

    
    data_red_avg = mean(rr_reg,3);
    figure; imagesc(data_red_avg);
    
    ImgFolderRed = expt(iexp).redImg;
    CD = [base '\Data\2P_images\' mouse '\' date '\' ImgFolderRed{1}];        
    cd(CD);
    imgMatFile = [expt(iexp).redImg{1} '_000_000.mat'];
    load(imgMatFile);
    nframes = info.config.frames;
    fprintf(['Reading run ' expt(iexp).redImg{1} '- ' num2str(min(nframes)) ' frames \r\n'])
    data = sbxread(imgMatFile(1,1:11),0,nframes);
    if size(data,1) == 2
        red_data = squeeze(data(2,:,:,:));
        green_data = squeeze(data(1,:,:,:));
        [out, green_data_reg] = stackRegister(green_data,data_avg);
        [out2, red_data_reg] = stackRegister_MA(red_data,[],[],out);
        red_data_avg = mean(red_data_reg,3);
        figure; imagesc(red_data_avg)
        title('Red')
        green_data_avg = mean(green_data_reg,3);
        figure; imagesc(green_data_avg)
        title('Green')
        if size(expt(iexp).redImg,2) == 2
            CD = [base '\Data\2P_images\' mouse '\' date '\' ImgFolderRed{2}];                  
            cd(CD);
            imgMatFile = [expt(iexp).redImg{2} '_000_000.mat'];
            load(imgMatFile);
            nframes = info.config.frames;
            fprintf(['Reading run ' expt(iexp).redImg{2} '- ' num2str(min(nframes)) ' frames \r\n'])
            data = sbxread(imgMatFile(1,1:11),0,nframes);
            red_data = squeeze(data(2,:,:,:));
            green_data = squeeze(data(1,:,:,:));
%             [out, green_data_reg] = stackRegister(mean(green_dat,3),green_data_avg);
%             [outs, red_data_reg] = stackRegister_MA(mean(red_data,3),[],[],out);
%             red_data_avg = red_data_reg;
            [out, green_data_reg] = stackRegister(green_data,green_data_avg);
            [outs, red_data_reg] = stackRegister_MA(red_data,[],[],out);
            red_data_avg = mean(red_data_reg,3);
            figure; imagesc(red_data_avg); movegui('center')
            title('Red')
        end
        save(fullfile([base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_redData.mat']), 'green_data_avg', 'red_data_avg')
    end

    data_avg = mean(data_reg(:,:,size(data_reg,3)-10000:end),3);
    figure; 
    subplot(2,2,1)
    sz = size(data_avg);
    rgb = zeros(sz(1),sz(2),3);
    rgb(:,:,1) = red_data_avg./max(red_data_avg(:));
    imagesc(red_data_avg);
    colormap gray
    subplot(2,2,2)
    rgb(:,:,2) = data_avg./max(data_avg(:));
    imagesc(rgb);
    if size(expt(iexp).redImg,2) == 2
        title('Red at 1040; Green at 920')
    else
        title('Red at 920; Green at 920')
    end
    print(fullfile([base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

end

%% find activated cells
clear data out
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

if input.doTwoStimTogether
    maskCon_all = celleqel2mat_padded(input.tStimTwoGratingContrast);
    maskCons = unique(maskCon_all);
    nMaskCon = length(maskCons);
    maskPhas_all = celleqel2mat_padded(input.tStimTwoGratingPhaseDeg);
    maskPhas = unique(maskPhas_all);
    nMaskPhas = length(maskPhas);
    maskDir_all = rad2deg(wrapTo2Pi(deg2rad(celleqel2mat_padded(input.tStimTwoGratingDirectionDeg))));
    maskDir_all(find(maskDir_all==360)) = 0;
    maskDirs = unique(maskDir_all);
    nMaskDir = length(maskDirs);
end

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

if nStimDir > 1 & ~input.doTwoStimTogether
    nStim2 = nStimDir.*(nMaskDiff+1);
    data_dfof = cat(3, data_dfof, zeros(sz(1),sz(2), nStim2));
    start = nStim1+1;
    for is = 1:nStimDir
        ind_stimalone = intersect(intersect(find(stimCon_all),find(maskCon_all==0)),find(stimDir_all == stimDirs(is)));
        ind_maskalone = intersect(intersect(find(stimCon_all==0),find(maskCon_all)),find(maskDir_all == maskDirs(is)));
        ind_alone(is) = length([ind_stimalone ind_maskalone]);
        data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,[ind_stimalone ind_maskalone]),3);
        start = start+1;
        for id = 1:nMaskDiff
            ind_plaidstim = intersect(intersect(intersect(find(stimCon_all),find(maskCon_all)),find(stimDir_all == stimDirs(is))),find(maskDiff_all == maskDiffs(id)));
            ind_n(is,id) = length(ind_plaidstim);
            data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,ind_plaidstim),3);
            start = start+1;
        end
    end
    figure; 
    subplot(2,1,1); 
    imagesc(mean(data_dfof(:,:,nStim1+1:1+nMaskDiff:end),3))
    title('Grating')
    colormap gray
    % clim([0 1])
    subplot(2,1,2); 
    imagesc(mean(data_dfof(:,:,nStim1+2:1+nMaskDiff:end),3))
    title('Plaid')
    colormap gray
    % clim([0 1])
    sgtitle([mouse ' ' date])
    print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_GratingVsPlaid.pdf']),'-dpdf')
elseif nStimDir > 1 & input.doTwoStimTogether
    nStim2 = nStimDir.*2;
    data_dfof = cat(3, data_dfof, zeros(sz(1),sz(2), nStim2));
    start = nStim1+1;
    for is = 1:nStimDir
        ind_same = intersect(find(stimDir_all == stimDirs(is)), find(maskDir_all == stimDirs(is)));
        data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,ind_same),3);
        ind_both = find(stimCon_all&maskCon_all);
        ind_diff = setdiff(find(stimDir_all == stimDirs(is)), find(maskDir_all == stimDirs(is)));
        data_dfof(:,:,start+1) = nanmean(data_resp_dfof(:,:,intersect(ind_both,ind_diff)),3);
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

if nTF>1
    nStim3 = nTF;
    data_dfof = cat(3, data_dfof, zeros(sz(1),sz(2), nStim3));
    start = nStim1+nStim2+1;
    for iTF = 1:nTF
        ind_TF = find(TF_all == TFs(iTF));
        data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,ind_TF),3);
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
print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_maxdfof.pdf']),'-dpdf')

data_dfof = cat(3, data_dfof, data_dfof_max);
if doRedChannel & (~doPhaseAfterDir & ~doDirAfterPass)
    data_dfof = cat(3,data_dfof,red_data_avg);
end
if (strcmp(expt(iexp).driver,'SOM') || strcmp(expt(iexp).driver,'PV')) & ~doRedChannel
    data_dfof = cat(3,data_dfof,data_avg);
end

save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']), 'cStimOn', 'maskCon_all', 'stimCon_all', 'stimCons', 'maskCons', 'nStimCon', 'nMaskCon', 'stimDir_all', 'stimDirs', 'nStimDir', 'maskDir_all', 'maskDirs', 'nMaskDir', 'maskPhas_all', 'maskPhas', 'nMaskPhas', 'maskDiff_all','maskDiffs','nMaskDiff','SF_all', 'SFs', 'nSF','TF_all', 'TFs', 'nTF', 'frame_rate', 'nTrials')

%% cell segmentation 
if ~doPhaseAfterDir & ~doDirAfterPass
    mask_exp = zeros(sz(1),sz(2));
    mask_all = zeros(sz(1), sz(2));
    mask_data = data_dfof;
    
    if strcmp(expt(iexp).driver{1},'SOM') || strcmp(expt(iexp).driver{1},'PV')
        if ~doRedChannel
            bwout = imCellEditInteractiveLG(mean(data_reg,3));
            mask_all = mask_all+bwout;
            mask_exp = imCellBuffer(mask_all,3)+mask_all;
            close all
        end
    end
    
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
end
    

clear data_adapt data_adapt_dfof data_test data_test_dfof data_test_avg data_resp data_resp_dfof bwout

%% neuropil subtraction
 
if ~doPhaseAfterDir & ~doDirAfterPass
    mask_np = imCellNeuropil(mask_cell, 3, 5);
    save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof', 'mask_cell', 'mask_np', 'red_cells')
else
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_mask_cell.mat']))
end
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

save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc', 'nCells', 'sz')

clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell


 