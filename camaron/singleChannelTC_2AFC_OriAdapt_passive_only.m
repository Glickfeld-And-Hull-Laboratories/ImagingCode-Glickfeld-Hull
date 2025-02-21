%%
% close all
clearvars
clc
dataset = 'oriAdapt_V1_cam';
eval(dataset);
iexp = 137;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
if strcmp(expt(iexp).folder,'lindsey')
    data_base = LG_base;
elseif strcmp(expt(iexp).folder,'camaron')
    data_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron';
end
mouse = expt(iexp).mouse;
date = expt(iexp).date;
%% Passive condition - Load runs
nrun = size(expt(iexp).pass_run,1);
for irun = 1:nrun
    pass_str = ['runs-' expt(iexp).pass_run];
    CD = [data_base '\Data\2P_images\' expt(iexp).mouse '\' expt(iexp).date '\' expt(iexp).pass_run];
    cd(CD);
    imgMatFile = [expt(iexp).pass_run expt(iexp).pass_run_suffix(irun,:) '.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' expt(iexp).mouse '-' expt(iexp).date '-' expt(iexp).pass_time '.mat'];
    load(fName);
    pass_input = input;
    clear input
    nframes = [pass_input.counterValues{end}(end) info.config.frames];
    
    if min(nframes)<pass_input.counterValues{end}(end)
        ntrials = size(pass_input.trialOutcomeCell,2);
        for itrial = ntrials:-1:1
            if pass_input.counterValues{itrial}(end) <= nframes
                break
            end
        end
        pass_input = trialChopper(pass_input,[1 itrial]);
    end
    
    fprintf(['Reading run ' expt(iexp).pass_run '- ' num2str(min(nframes)) ' frames \r\n'])
    data = sbxread(imgMatFile(1,1:11),0,min(nframes));
    data = squeeze(data);
end

fprintf('Runs loaded\n')

%% Choose register interval
% Should be defined via code based on the total number of frames
nep = floor(size(data,3)./5000); 
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*5000):500+((i-1)*5000)),3)); title([num2str(1+((i-1)*5000)) '-' num2str(500+((i-1)*5000))]); end

%% Register data
if exist(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_reg_shifts.mat']))
    load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_reg_shifts.mat']))
    [out, data_reg] = stackRegister_MA(data,[],[],out); 
else
    data_avg = mean(data(:,:,20001:20500),3); % User defined interval!!! (Why not hard code?)
    [out, data_reg] = stackRegister(data,data_avg);
    smooth_out_x = smooth(out(:,3),200);
    smooth_out_y = smooth(out(:,4),200);
    diff_out = max([abs(smooth_out_x-out(:,3)) abs(smooth_out_y-out(:,4))], [],2);
    move_ind = find(diff_out>10);
    mkdir(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str]))
    save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_reg_shifts.mat']), 'out', 'data_avg', 'diff_out', 'move_ind')
    save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_input.mat']), 'pass_input')
end
clear data

%% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*5000):500+((i-1)*5000)),3)); title([num2str(1+((i-1)*5000)) '-' num2str(500+((i-1)*5000))]); end
figure; imagesq(mean(data_reg(:,:,1:5000),3)); truesize;
print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

%% Code for red data....

% *****
%% Check 2AFC photodiode (info must exist)
tic
irun = 1;
nrun = 1;
% run_str = ['runs']; 
% run_str = [run_str '-' expt(iexp).runs(irun,:)];

if exist(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_input.mat']))
    load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
    load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_input.mat']), 'pass_input')
end

nframes = [pass_input.counterValues{end}(end) info.config.frames];
ntrials = size(pass_input.trialOutcomeCell,2);




if exist([data_base '\Data\2P_images\' expt(iexp).mouse '\' expt(iexp).date '\' expt(iexp).pass_run(irun,:) '\' [expt(iexp).pass_run(irun,:) expt(iexp).pass_run_suffix(irun,:) '.ephys']]);
    photoData = [];
    for irun = 1:nrun
        filename = [data_base '\Data\2P_images\' expt(iexp).mouse '\' expt(iexp).date '\' expt(iexp).pass_run(irun,:) '\' [expt(iexp).pass_run(irun,:) expt(iexp).pass_run_suffix(irun,:) '.ephys']];
        fileID = fopen(filename, 'r', 'ieee-le');
        if fileID == -1, error('Cannot open file: %s', filename); end
        format = 'uint32';
        photoData = [photoData; fread(fileID, Inf, format)];
        fclose(fileID);
    end

    [photoLoc stimOnFrames] = photoFrameFinder_movBase(photoData,min(nframes));
    frameDiff = diff(stimOnFrames);
    ind_long = find(frameDiff>20);
    ind_long_long = ind_long(find(frameDiff(ind_long-1)>20)); % ??
    photoLoc(ind_long_long) = [];
    stimOnFrames(ind_long_long) = [];

    nf = rem(size(stimOnFrames,2),5); % rem = remainder after division; Why divide by 5? Adaptors plus target?
    photoLoc_rs = reshape(photoLoc(1:end-nf),[5 length(photoLoc(1:end-nf))./5])'; 
    photoLoc_diff = diff(photoLoc_rs,1,2);
    figure; plot(photoLoc_diff'); ylim([0 6000])

    tDoFB = celleqel2mat_padded(pass_input.tDoFeedbackMotion);
    tFramesStimOn = celleqel2mat_padded(pass_input.cStimOff)-celleqel2mat_padded(pass_input.cStimOn);
    ind_fast = tFramesStimOn<pass_input.nFramesTooFast;
    FBfast = tDoFB & ind_fast;
    ntrials = size(pass_input.tGratingContrast,2);
    cAdaptOn = nan(1,ntrials);
    cStimOn = nan(1,ntrials);
    n1 = 1; % First adaptor (distractor)
    n2 = 5; % Stimulus (presentation)
    cAdaptOn(1) = stimOnFrames(n1);
    cStimOn(1) = stimOnFrames(n2);
    for itrial = 2:ntrials
        if FBfast(itrial-1)
            n1 = n1+6;
            n2 = n2+6;
        else
            n1 = n1+5;
            n2 = n2+5;
        end
        cAdaptOn(itrial) = stimOnFrames(n1);
        cStimOn(itrial) = stimOnFrames(n2);
    end

    unique(cStimOn-cAdaptOn)
    print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str 'photoLoc_diff.pdf']),'-dpdf', '-bestfit')
    save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_photoData.mat']), 'photoLoc', 'stimOnFrames', 'cStimOn','cAdaptOn')
elseif isfield(info, "frame")
    [stimOnFrames stimOffFrames] = photoFrameFinder_Sanworks(info.frame);

    tDoFB = celleqel2mat_padded(pass_input.tDoFeedbackMotion);
    tFramesStimOn = celleqel2mat_padded(pass_input.cStimOff)-celleqel2mat_padded(pass_input.cStimOn);
    ind_fast = tFramesStimOn<pass_input.nFramesTooFast;
    FBfast = tDoFB & ind_fast;
    cAdaptOn = nan(1,ntrials);
    cStimOn = nan(1,ntrials);
    n1 = 1; % First adaptor (distractor)
    n2 = 5; % Stimulus (presentation)
    cAdaptOn(1) = stimOnFrames(n1);
    cStimOn(1) = stimOnFrames(n2);
    for itrial = 2:ntrials 
        if FBfast(itrial-1)
            n1 = n1+6;
            n2 = n2+6;
        else
            n1 = n1+5;
            n2 = n2+5;
        end
        cAdaptOn(itrial) = stimOnFrames(n1);
        cStimOn(itrial) = stimOnFrames(n2);
    end

    unique(cStimOn-cAdaptOn)
    save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_photoData.mat']), 'stimOnFrames', 'cStimOn','cAdaptOn')
    

else
    error("No photodiode data!!!")
end

clearvars photoData
toc

%% find activated cells
% cAdaptOn = celleqel2mat_padded(pass_input.cAdaptOn);
% cStimOn = celleqel2mat_padded(pass_input.cStimOn);

close all
tGratingOri = celleqel2mat_padded(pass_input.tGratingDirectionStart);
Oris = unique(tGratingOri);
nOri = length(Oris);
nTrials = length(tGratingOri);
sz = size(data_reg);
data_f = nan(sz(1),sz(2),nTrials);
data_adapt = nan(sz(1),sz(2),nTrials);
data_targ = nan(sz(1),sz(2),nTrials);
rem_trial = zeros(1,nTrials);
for itrial = 1:nTrials
    ind = cAdaptOn(itrial)-20:cStimOn(itrial);
    if sum(find(ind == move_ind)) == 0 
        if cAdaptOn(itrial)+25<sz(3)
            data_f(:,:,itrial) = mean(data_reg(:,:,cAdaptOn(itrial)-20:cAdaptOn(itrial)-1),3);
            data_adapt(:,:,itrial) = mean(data_reg(:,:,cAdaptOn(itrial)+5:cAdaptOn(itrial)+25),3);
        end
        if cStimOn(itrial)+15<sz(3)
            data_targ(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)+5:cStimOn(itrial)+15),3);
        end
    else
        rem_trial(1,itrial) = 1;
    end
end
data_adapt_dfof = (data_adapt-data_f)./data_f;
data_targ_dfof = (data_targ-data_f)./data_f;

[n n2] = subplotn(nOri+3);
figure;
data_dfof = zeros(sz(1),sz(2),nOri+3);
for iori = 1:nOri
    subplot(n,n2,iori)
    ind = find(tGratingOri == Oris(iori));
    data_dfof(:,:,iori)= nanmean(data_targ_dfof(:,:,ind),3);
    imagesc(data_dfof(:,:,iori));
    title([num2str(Oris(iori)) ' deg'])
end
subplot(n,n2,iori+1)
data_dfof(:,:,iori+1)= nanmean(data_targ_dfof,3);
imagesc(data_dfof(:,:,iori+1));
title('All')
aContrast = celleqel2mat_padded(pass_input.aGratingContrast);
b2Ix = celleqel2mat_padded(pass_input.tBlock2TrialNumber);
dirs = [0 90];
for idir = 1:2
    ind = find(aContrast == 1 & b2Ix == idir-1);
    data_dfof(:,:,iori+idir+1)= nanmean(data_adapt_dfof(:,:,ind),3);
    subplot(n,n2,iori+idir+1)
    imagesc(data_dfof(:,:,iori+idir+1));
    title(['Adapt- ' num2str(dirs(idir)) ' deg'])
end

data_dfof_max = max(data_dfof,[],3);
figure; imagesc(data_dfof_max)
data_dfof = cat(3,data_dfof_max,data_dfof);
print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_activated_cells.pdf']),'-dpdf', '-bestfit')

%% direction tuning
clear data_targ data_adapt data_f
dir_str = ['runs-' expt(iexp).dirtuning];
CD = [data_base '\Data\2P_images\' expt(iexp).mouse '\' expt(iexp).date '\' expt(iexp).dirtuning];
cd(CD);
imgMatFile = [expt(iexp).dirtuning expt(iexp).dirtuning_suffix(irun,:) '.mat'];
load(imgMatFile);
fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' expt(iexp).mouse '-' expt(iexp).date '-' expt(iexp).dirtuning_time '.mat'];
load(fName);
dir_input = input;
clear input
nframes_dir = info.config.frames;

fprintf(['Reading run ' expt(iexp).dirtuning '- ' num2str(nframes_dir) ' frames \r\n'])
data = sbxread(imgMatFile(1,1:11),0,nframes_dir);
data = squeeze(data);

if exist(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' dir_str], [date '_' mouse '_' dir_str '_stimData.mat']))
    load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' dir_str], [date '_' mouse '_' dir_str '_stimData.mat']))
    load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' dir_str], [date '_' mouse '_' dir_str '_reg_shifts.mat']))
    [out, data_reg_dir] = stackRegister_MA(data,[],[],out);
    clear data
else  
    run_str = ['runs']; 
    for irun = 1:nrun
        run_str = [run_str '-' expt(iexp).pass_run(irun,:)];
    end
    load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))

    [out, data_reg_dir] = stackRegister(data,data_avg);
    mkdir(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' dir_str]))
    save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' dir_str], [date '_' mouse '_' dir_str '_reg_shifts.mat']), 'out', 'data_avg')
    clear data
    
    dir_mat = celleqel2mat_padded(dir_input.tGratingDirectionDeg);
    Dirs = unique(dir_mat);
    nDir = length(Dirs);
    nTrials = length(dir_mat);
    sz = size(data_reg_dir);
    nOn = dir_input.nScansOn;
    nOff = dir_input.nScansOff;
    data_tr = reshape(data_reg_dir, [sz(1) sz(2) nOn+nOff nTrials]);
    data_f = mean(data_tr(:,:,ceil(nOff/2):nOff, :),3);
    dir_dfof = (double(data_tr)-data_f)./data_f;
    dir_dfof_all = zeros(sz(1),sz(2),nDir);
    [n n2] = subplotn(nDir);
    figure;
    for idir = 1:nDir
        ind = find(dir_mat == Dirs(idir));
        subplot(n,n2,idir)
        dir_dfof_all(:,:,idir) = mean(mean(dir_dfof(:,:,nOff+1:nOn+nOff, ind),3),4);
        imagesc(dir_dfof_all(:,:,idir))
    end
    save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' dir_str], [date '_' mouse '_' dir_str '_input.mat']), 'dir_input')
    save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' dir_str], [date '_' mouse '_' dir_str '_stimData.mat']), 'dir_mat', 'Dirs', 'nDir', 'dir_dfof_all')
    clear data_tr data_f
end

data_dfof = cat(3,cat(3,data_dfof,max(dir_dfof_all,[],3)),dir_dfof_all);

%% cell segmentation  
% if exist(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_mask_cell.mat']))
%     load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_mask_cell.mat']))
%     disp("Loaded mask_cell")

if strcmp(cell2mat(expt(iexp).img_strct),'cells')
    mask_exp = zeros(sz(1),sz(2));
    mask_all = zeros(sz(1), sz(2));

    if ~isempty(expt(iexp).redImg)
        bwout = imCellEditInteractiveLG(red_data_avg);
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
    print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_bouton_mask.pdf']), '-dpdf')

    mask_np = imCellNeuropil(mask_cell, 3, 5);
    if ~isempty(expt(iexp).redImg)
        save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_cell_red', 'mask_np')
    else
        save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_np')  
    end
elseif strcmp(cell2mat(expt(iexp).img_strct),'axons')
    fprintf('\nBegin axon segmentation...')
    min_df = 0.05;
    b = 5;
    mask_cell = zeros(sz(1), sz(2));
    mask_all = zeros(sz(1), sz(2));
    temp_data = ones(sz(1),sz(2));
    data_dfof_avg = data_dfof;
    data_dfof_max = max(data_dfof_avg,[],3);
    temp_max = squeeze(max(max(data_dfof_avg,[],1),[],2));
    [a, max_sort] = sort(temp_max,'descend');
    data_dfof_avg_sort = cat(3,data_dfof_max,data_dfof_avg(:,:,max_sort));
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
            if i>b & j>b & i<sz(1)-b & j<sz(2)-b
                if ~isnan(temp_data(i-1:i+1,j-1:j+1))
                    all_pix = data_dfof_avg_sort(i-1:i+1,j-1:j+1,ia);
                    [max_val max_ind] = max(all_pix(:));
                    if max_ind == 5
                        h = zeros(1, 2+nOri+nDir);
                        dirs = [0 90];
                        start = 0;
                        for idir = 1:2
                            ind = find(aContrast == 1 & b2Ix == idir-1);
                            [h(1,idir) p] = ttest(squeeze(mean(mean(data_adapt_dfof(i-1:i+1,j-1:j+1,ind),1),2)),0,'tail','right','alpha',0.05./(1+nOri+nDir));
                        end
                        start = start+2;
                        for iOri = 1:nOri
                            ind = find(tGratingOri == Oris(iori));
                            [h(1,start+iOri) p] = ttest(squeeze(mean(mean(data_targ_dfof(i-1:i+1,j-1:j+1,ind),1),2)),0,'tail','right','alpha',0.05./(1+nOri+nDir));
                        end
                        start = start+nOri;
                        for iDir = 1:nDir
                            ind = find(dir_mat == Dirs(iori));
                            [h(1,start+iDir) p] = ttest(squeeze(mean(mean(dir_dfof(i-1:i+1,j-1:j+1,ind),1),2)),0,'tail','right','alpha',0.05./(1+nOri+nDir));
                        end
                        if sum(h(:))>2
                            mask_cell(i,j) = 1;
                            mask_all(i-1:i+1,j-1:j+1) = ones(3,3);
                            temp_data(i-2:i+2,j-2:j+2) = NaN(5,5);
                            fprintf('.')
                        end
                    end
                end
            end
        end
    end      
    mask_cell = bwlabel(mask_all);
    figure; imagesc(mask_cell)
    mask_cell_red = zeros(size(mask_cell)); 
    print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_bouton_mask.pdf']), '-dpdf')

    save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_max', 'mask_cell') 
end


%% neuropil subtraction
clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_2 data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 
data_tc = stackGetTimeCourses(data_reg, mask_cell);
nCells = size(data_tc,2);
if strcmp(cell2mat(expt(iexp).img_strct),'cells')
    data_tc_down = stackGetTimeCourses(stackGroupProject(data_reg,5), mask_cell);
    clear np_tc np_tc_down
    sz = size(data_reg);
    down = 5;
    data_reg_down  = stackGroupProject(data_reg,down);
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
    save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc','np_tc','npSub_tc')
elseif strcmp(cell2mat(expt(iexp).img_strct),'axons')
    npSub_tc = data_tc;
    save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'npSub_tc','-v7.3')
end
clear data_reg data_reg_down
clear data_tc np_tc
%% Passive Condition *Update TC's (dfof) with photodiode data - start here

run_str = ['runs']; 
for irun = 1:nrun
    run_str = [run_str '-' expt(iexp).pass_run(irun,:)];
end
load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'mask_cell', 'mask_np') 

if ~isempty(expt(iexp).pass_run)

    frameRateHz = pass_input.frameRateHz;    
    tGratingOri = celleqel2mat_padded(pass_input.tGratingDirectionStart);
%     cStimOn = celleqel2mat_padded(pass_input.cStimOn);
%     cAdaptOn = celleqel2mat_padded(pass_input.cAdaptOn);
    b2Ix = celleqel2mat_padded(pass_input.tBlock2TrialNumber);
    tOris = unique(tGratingOri);
    nOri = length(tOris);
    aGratingOri = celleqel2mat_padded(pass_input.aGratingDirectionDeg);
    aGratingContrast = celleqel2mat_padded(pass_input.aGratingContrast);
    aCons = unique(aGratingContrast);
    naCon = length(aCons);
    aOris = unique(aGratingOri);
    naOri = length(aOris);
    nCells = size(npSub_tc,2);
    nframes = size(npSub_tc,1);
    nTrials = size(aGratingOri,2);
    data_stim = nan(50,nCells,nTrials);
    data_adapt = nan(100,nCells,nTrials);
    for itrial = 1:nTrials
        if cStimOn(itrial)+29 < nframes
            data_stim(:,:,itrial) = npSub_tc(cStimOn(itrial)-20:cStimOn(itrial)+29,:);
        end
        if cAdaptOn(itrial)+79< nframes
            data_adapt(:,:,itrial) = npSub_tc(cAdaptOn(itrial)-20:cAdaptOn(itrial)+79,:);
        end
    end
    dataf = mean(data_adapt(1:20,:,:),1);
    data_stim_dfof = bsxfun(@rdivide, bsxfun(@minus, data_stim, dataf), dataf);
    data_adapt_dfof = bsxfun(@rdivide, bsxfun(@minus, data_adapt, dataf), dataf);
    tt = [-20:29].*(1000./frameRateHz);
    tt_adapt = [-20:79].*(1000./frameRateHz);
    figure;
    subplot(1,2,1)
    plot(nanmean(mean(data_adapt_dfof,2),3));
%     vline([16 20 25 29])
    vline([20 31 42 53])
    title('Adapt')
    subplot(1,2,2)
    plot(nanmean(mean(data_stim_dfof,2),3));
    vline([20 31 42 53])
    title('Target')

    nt = cell(3,nOri);
    x = 1;
    start = 1;
%     figure;
%     for icon = 1:naCon
%         if aCons(icon) == 1
%             for iaOri = 1:naOri
%                 ind_con = find(aGratingContrast == aCons(icon));
%                 ind_con = intersect(ind_con, find(aGratingOri == aOris(iaOri)));
%                 for iori = 1:nOri
%                     ind_ori = intersect(ind_con, find(tGratingOri == tOris(iori)));
%                     subplot(3,nOri,start)
%                     plot(tt, mean(nanmean(data_stim_dfof(:,:,ind_ori),3),2))
%                     hold on
%                     title(['Adapt: Ori = '  num2str(aOris(iaOri))])
%                     ylim([-0.01 0.1])
%                     start = start+1;
%                     nt{x,iori} = ind_ori;
%                 end
%                 x = 1+x;
%             end
%         else
%             ind_con = find(aGratingContrast == aCons(icon));
%             for iori = 1:nOri
%                 ind_ori = intersect(ind_con, find(tGratingOri == tOris(iori)));
%                 subplot(3,nOri,start)
%                 plot(tt, mean(nanmean(data_stim_dfof(:,:,ind_ori),3),2))
%                 hold on
%                 title(['No adapt; Ori = ' num2str(tOris(iori))])
%                 ylim([-0.01 0.1])
%                 start = start+1;
%                 nt{x,iori} = ind_ori;
%             end
%             x = 1+x;
%         end
%     end
%     print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_avgTargetResponse_byTarget.pdf']), '-dpdf','-bestfit')

    SIx = strcmp(pass_input.trialOutcomeCell,'success');

    ind_aCon0 = find(aGratingContrast == 0);
    ind_aCon1 = find(aGratingContrast);
    ind_cond{1} = intersect(find(SIx),ind_aCon0);
    ind_cond{2} = intersect(find(SIx),intersect(ind_aCon1,find(aGratingOri==0)));
    ind_cond{3} = intersect(find(SIx),intersect(ind_aCon1,find(aGratingOri==90)));

    aCon_p1 = [NaN aGratingContrast];
    aOri_p1 = [NaN aGratingOri];
    ind_aCon0_p1{1} = intersect(ind_aCon0,find(aCon_p1==0));
    ind_aCon0_p1{2} = intersect(ind_aCon0,intersect(find(aCon_p1==1),find(aOri_p1 ==0)));
    ind_aCon0_p1{3} = intersect(ind_aCon0,intersect(find(aCon_p1==1),find(aOri_p1 ==90)));


%     figure;
%     subplot(1,2,1)
%     for i = 1:3
%      plot(tt, mean(mean(data_stim_dfof(:,:,ind_cond{i}),3),2))
%      hold on
%     end
%     title('Current trial')
%     ylabel('dF/F')
%     xlabel('Time from target (ms)')
%     legend({'Con = 0','Ori = 0', 'Ori = 90'})
%     subplot(1,2,2)
%     for i = 1:3
%      plot(tt, mean(mean(data_stim_dfof(:,:,ind_aCon0_p1{i}),3),2))
%      hold on
%     end
%     title('Previous trial, for Current Con = 0')
%     ylabel('dF/F')
%     xlabel('Time from target (ms)')
%     print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_avgTargetResponse.pdf']), '-dpdf','-bestfit')

    mask_label = zeros(1,nCells);
   
%     for i = 1:nCells
%         if mask_cell_red(find(mask_cell == i, 1))
%             mask_label(1,i) = 1;
%         end
%     end

    interval = ceil(64/nOri)+1;
    x = 1:interval:64;
    x(end) = 63;
    y = bluered;

%     figure;
%     if nCells <49
%         [n n2] = subplotn(nCells);
%     else
%         [n n2] = subplotn(49);
%     end
%     start = 1;
%     for iC = 1:nCells
%         if start > 49
%             sgtitle('Adapt response')
%             print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_adaptResp_allCells_success_cells' num2str(iC-49) '-' num2str(iC-1) '.pdf']), '-dpdf','-bestfit')
%             figure;
%             start = 1;
%         end
%         subplot(n,n2,start)
%         ind_adapt = find(SIx & aGratingContrast & b2Ix==0);
%         plot(mean(data_adapt_dfof(:,iC,ind_adapt),3),'Color',y(x(1),:))
%         hold on
%         ind_adapt = find(SIx & aGratingContrast & b2Ix);
%         plot(tt_adapt, mean(data_adapt_dfof(:,iC,ind_adapt),3),'Color',y(x(end),:))
%         ylim([-0.05 0.3])
%         if mask_label(iC)
%             title([num2str(iC) '-' expt(iexp).driver])
%         else
%             title(num2str(iC))
%         end
%         start = start+1;
%     end
%     sgtitle('Adapt response')
%     print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_adaptResp_allCells_success_cells' num2str(iC-48) '-' num2str(iC) '.pdf']), '-dpdf','-bestfit')

%     figure;
%     start = 1;
%     for iC = 1:nCells
%         if start > 49
%             sgtitle('Target response- All trials- All correct')
%             print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_targetResp_allCells_success_cells' num2str(iC-49) '-' num2str(iC-1) '.pdf']), '-dpdf','-bestfit')
%             figure;
%             start = 1;
%         end
%         subplot(n,n2,start)
%         for iori = 1:nOri
%             ind_ori = find(SIx & tGratingOri == tOris(iori));
%             plot(tt, mean(data_stim_dfof(:,iC,ind_ori),3),'Color',y(x(iori),:))
%             hold on
%         end
%         if mask_label(iC)
%             title([num2str(iC) '-' expt(iexp).driver])
%         else
%             title(num2str(iC))
%         end
%         start = start+1;
%     end
%     sgtitle('Target response- All trials- All correct')
%     print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_targetResp_allCells_success_cells' num2str(iC-48) '-' num2str(iC) '.pdf']), '-dpdf','-bestfit')


%     figure;
%     subplot(3,2,1)
%     ind_adapt = find(SIx & aGratingContrast & b2Ix==0);
%     plot(tt_adapt, mean(mean(data_adapt_dfof(:,find(mask_label),ind_adapt),3),2),'Color',y(x(1),:))
%     hold on
%     ind_adapt = find(SIx & aGratingContrast & b2Ix);
%     plot(tt_adapt, mean(mean(data_adapt_dfof(:,find(mask_label),ind_adapt),3),2),'Color',y(x(end),:))
%     title([expt(iexp).driver '+ n = ' num2str(sum(mask_label))])
%     xlabel('Time from adapt on (ms)')
%     ylim([-0.02 0.15])
%     subplot(3,2,2)
%     ind_adapt = find(SIx & aGratingContrast & b2Ix==0);
%     plot(tt_adapt, mean(mean(data_adapt_dfof(:,find(~mask_label),ind_adapt),3),2),'Color',y(x(1),:))
%     hold on
%     ind_adapt = find(SIx & aGratingContrast & b2Ix);
%     plot(tt_adapt, mean(mean(data_adapt_dfof(:,find(~mask_label),ind_adapt),3),2),'Color',y(x(end),:))
%     title([expt(iexp).driver '+ n = ' num2str(sum(~mask_label))])
%     xlabel('Time from adapt on (ms)')
%     ylim([-0.02 0.15])
%     subplot(3,2,3) 
%     for iori = 1:nOri
%         ind_ori = find(SIx & tGratingOri == tOris(iori));
%         plot(tt, mean(mean(data_stim_dfof(:,find(mask_label),ind_ori),3),2),'Color',y(x(iori),:))
%         hold on
%     end
%     xlabel('Time from target (ms)')
%     ylim([-0.02 0.15])
%     subplot(3,2,4)
%     for iori = 1:nOri
%         ind_ori = find(SIx & tGratingOri == tOris(iori));
%         plot(tt, mean(mean(data_stim_dfof(:,find(~mask_label),ind_ori),3),2),'Color',y(x(iori),:))
%         hold on
%     end
%     xlabel('Time from target (ms)')
%     ylim([-0.02 0.15])
%     sgtitle(['All corrects - n = ' num2str(sum(SIx))])
%     print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_avgResp_allCells_success.pdf']), '-dpdf','-bestfit')
end
%% adapt analysis - passive *Update TC's with photodiode data - end here
if ~isempty(expt(iexp).pass_run)

    close all
    base_win = 19:21;
    resp_win = 25:27;
    data_adapt_base = squeeze(nanmean(data_adapt_dfof(base_win,:,:),1));
    data_adapt_resp = squeeze(nanmean(data_adapt_dfof(resp_win,:,:),1));
    ind_con = find(aGratingContrast == 1);
    h_adapt = zeros(nCells,naOri);
    p_adapt = zeros(nCells,naOri);
    adapt_resp_ind = cell(1,naOri);
    for iOri = 1:naOri
        ind_ori = find(aGratingOri==aOris(iOri));
        ind = intersect(ind_ori,ind_con);
        [h_adapt(:,iOri), p_adapt(:,iOri)] = ttest(data_adapt_resp(:,ind), data_adapt_base(:,ind), 'tail', 'right','dim',2);
        adapt_resp_ind{iOri} = find(h_adapt(:,iOri));
    end

    data_stim_base = squeeze(nanmean(data_stim_dfof(base_win,:,:),1));
    data_stim_resp = squeeze(nanmean(data_stim_dfof(resp_win,:,:),1));
    h_stim = zeros(nCells,nOri);
    p_stim = zeros(nCells,nOri);
    for iOri = 1:nOri
        ind_ori = find(tGratingOri==tOris(iOri));
        [h_stim(:,iOri), p_stim(:,iOri)] = ttest(data_stim_resp(:,ind_ori), data_stim_base(:,ind_ori), 'tail', 'right','dim',2);
    end
    stim_resp_ind = find(sum(h_stim,2));
    stim_resp = data_stim_resp-data_stim_base;

    % if ~isfield(pass_input,'cAdaptOn_all')
        nFramesPerAdapt = double(ceil(((double(pass_input.dynAdaptPeriodMs + pass_input.dynAdaptFlashOffMs))./1000).*frameRateHz));
        base_win_all = zeros(length(base_win),4,2);
        resp_win_all = zeros(length(resp_win),4,2);
        for i = 1:4
            base_win_all(:,i,1) = base_win+((i-1).*nFramesPerAdapt);
            resp_win_all(:,i,1) = resp_win+((i-1).*nFramesPerAdapt);
            base_win_all(:,i,2) = base_win+((i-1).*nFramesPerAdapt);
            resp_win_all(:,i,2) = resp_win-2+((i-1).*nFramesPerAdapt);
        end
    % else
    %     error('need code for cAdaptOn_all')
    % end

    figure;
    start = 1;
    for iOri = 1:naOri
        for IN = 1:2
            subplot(2,2,start)
            ind_ori = find(aGratingOri==aOris(iOri));
            trial_ind = intersect(ind_ori,ind_con);
            cell_ind = intersect(adapt_resp_ind{iOri}, find(mask_label==IN-1));
            plot(nanmean(nanmean(data_adapt_dfof(:,cell_ind,trial_ind),3),2))
            if IN == 1
                IN_str = [expt(iexp).driver '-'];
            else
                IN_str = [expt(iexp).driver '+'];
            end
            vline(reshape(base_win_all(:,:,IN),[size(base_win_all,1).*size(base_win_all,2) 1]), 'b')
            vline(reshape(resp_win_all(:,:,IN),[size(resp_win_all,1).*size(resp_win_all,2) 1]), 'r')
            title(['Ori = ' num2str(aOris(iOri)) '; ' num2str(length(cell_ind)) ' ' IN_str ' cells'])
            start = start+1;
        end
    end

    adapt_cyc_resp = zeros(nCells,nTrials,4);
    for IN = 1:2
        for i = 1:4
            adapt_cyc_resp(find(mask_label==IN-1),:,i) = squeeze(nanmean(data_adapt_dfof(resp_win_all(:,i,IN),find(mask_label==IN-1),:),1)...
                -nanmean(data_adapt_dfof(base_win_all(:,i,IN),find(mask_label==IN-1),:),1));
        end
    end


    save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_adaptResp.mat']), 'mask_label', 'data_adapt_dfof', 'adapt_cyc_resp', 'base_win_all', 'resp_win_all','adapt_resp_ind')
    save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_stimResp.mat']), 'stim_resp', 'data_stim_dfof', 'stim_resp_ind', 'base_win', 'resp_win');
    save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_stimData.mat']), 'tGratingOri', 'tOris', 'aGratingOri', 'aOris', 'aGratingContrast', 'ind_cond');

end

%% dir tuning neuropil subtraction
close all
data_tc = stackGetTimeCourses(data_reg_dir, mask_cell);
if strcmp(cell2mat(expt(iexp).img_strct),'cells')
    data_tc_down = stackGetTimeCourses(stackGroupProject(data_reg_dir,5), mask_cell);
    nCells = size(data_tc,2);
    %np_tc = stackGetTimeCourses(data_reg,mask_np);
    clear np_tc np_tc_down
    sz = size(data_reg_dir);
    down = 5;
    data_reg_down  = stackGroupProject(data_reg_dir,down);
    np_tc = zeros(sz(3),nCells);
    np_tc_down = zeros(floor(sz(3)./down), nCells);
    for i = 1:nCells
         np_tc(:,i) = stackGetTimeCourses(data_reg_dir,mask_np(:,:,i));
         np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
         fprintf(['Cell #' num2str(i) '\n']) 
    end
    %get weights by maximizing skew
    ii= 0.01:0.01:1;
    x = zeros(length(ii), nCells);
    for i = 1:100
        x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
    end
    [max_skew ind] =  max(x,[],1);
    np_w = 0.01*ind;
    npSub_tc_dir = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
    clear data_reg_dir data_reg_down
    save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' dir_str], [date '_' mouse '_' dir_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc_dir')
elseif strcmp(cell2mat(expt(iexp).img_strct),'axons')
    npSub_tc_dir = data_tc;
    save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' dir_str], [date '_' mouse '_' dir_str '_TCs.mat']), 'npSub_tc_dir', '-v7.3')
end
clear data_tc np_tc data_reg_dir dir_dfof

%% Direction analysis
frameRateHz = 30;
nOn = dir_input.nScansOn;
nOff = dir_input.nScansOff;
dir_mat = celleqel2mat_padded(dir_input.tGratingDirectionDeg);
nTrials = length(dir_mat);
dir_input.trialSinceReset = nTrials;

down = 10;
nframes = size(npSub_tc_dir,1)./down;
data_tc_down = squeeze(mean(reshape(npSub_tc_dir, [down,nframes,nCells]),1));

tuningDownSampFactor = down;
[avgResponseEaOri,semResponseEaOri,vonMisesFitAllCellsAllBoots,fitReliability,R_square,tuningTC] = ...
    getOriTuning(data_tc_down,dir_input,tuningDownSampFactor,frameRateHz);
    vonMisesFitAllCells = squeeze(vonMisesFitAllCellsAllBoots(:,1,:));

if nCells<500
save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' dir_str], [date '_' mouse '_' dir_str '_oriTuningAndFits.mat']),...
            'avgResponseEaOri','semResponseEaOri','vonMisesFitAllCellsAllBoots','fitReliability','R_square', 'tuningTC')
else        
save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' dir_str], [date '_' mouse '_' dir_str '_oriTuningAndFits.mat']),...
            'avgResponseEaOri','semResponseEaOri','vonMisesFitAllCellsAllBoots','fitReliability','R_square', 'tuningTC','-v7.3')        
end
%% Ori Analysis
% if strcmp(cell2mat(expt(iexp).img_strct),'axons')
%     mask_label = zeros(1,nCells);
% 
% elseif strcmp(cell2mat(expt(iexp).img_strct),'cells')
%     mask_label = zeros(1,nCells);
%     for i = 1:nCells
%         if mask_cell_red(find(mask_cell == i, 1))
%             mask_label(1,i) = 1;
%         end
%     end
% end


dir_mat = celleqel2mat_padded(dir_input.tGratingDirectionDeg);
ori_mat = dir_mat;
ori_mat(find(dir_mat>=180)) = ori_mat(find(dir_mat>=180))-180;
oris = unique(ori_mat);
figure; 
if nCells<49
    [n n2] = subplotn(nCells);
else
    [n, n2] = subplotn(49);
end
if nCells>250
    nC = 250;
else
    nC = nCells;
end
start = 1;
x = 0;
for ic = 1:nC
    if start > 49
        sgtitle([expt(iexp).mouse ' ' expt(iexp).date ' n = ' num2str(length(find(fitReliability(find(mask_label))<22.5))) '/' num2str(sum(mask_label)) expt(iexp).driver '+ and ' num2str(length(find(fitReliability(find(~mask_label))<22.5))) '/' num2str(sum(~mask_label)) expt(iexp).driver '- well-fit'])
        print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' dir_str], [date '_' mouse '_' dir_str '_oriTuningFits_cells' num2str(start-49) '-' num2str(start-1) '.pdf']),'-dpdf','-fillpage')
        start = 1;
        x = x+1;
        figure;
    end
    subplot(n,n2,ic-(x.*49))
    errorbar(oris,avgResponseEaOri(ic,:), semResponseEaOri(ic,:),'-o')
    hold on
    plot(0:180,vonMisesFitAllCellsAllBoots(:,1,ic));
    tit_str = num2str(chop(R_square(1,ic),2));
    if mask_label(ic)
        tit_str = [tit_str '-' expt(iexp).driver];
    end
    if fitReliability(ic)<22.5
        tit_str = [tit_str '- R'];
    end
    title(tit_str)
    start = start+1;
end
sgtitle([expt(iexp).mouse ' ' expt(iexp).date ' n = ' num2str(length(find(fitReliability(find(mask_label))<22.5))) '/' num2str(sum(mask_label)) expt(iexp).driver '+ and ' num2str(length(find(fitReliability(find(~mask_label))<22.5))) '/' num2str(sum(~mask_label)) expt(iexp).driver '- well-fit'])
print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' dir_str], [date '_' mouse '_' dir_str '_oriTuningFits' num2str(start-49) '-' num2str(start-1) '.pdf']),'-dpdf','-fillpage')

[max_resp prefOri] = max(vonMisesFitAllCellsAllBoots,[],1);
prefOri = squeeze(prefOri)-1;
prefOri_bootdiff = abs(prefOri(2:end,:)-prefOri(1,:));
prefOri_bootdiff(find(prefOri_bootdiff>90)) = 180-prefOri_bootdiff(find(prefOri_bootdiff>90));
ind_theta90 = find(prctile(prefOri_bootdiff,90,1)<22.5);
% edges = [0 22.5:45:180]; 
% edges = 0:45:180;
edges = [0 22.5:45:180 180];

[bin edges ind_bin] = histcounts(prefOri(1,:),edges);

% [bin ind_bin] = histc(prefOri(1,:),edges);
ind_bin(find(ind_bin==5)) = 1;
bin(1) = bin(1)+bin(5);
bin(5) = [];

tunedCells = cell(1,length(bin));
for i = 1:length(bin)
    tunedCells{i} = intersect(find(ind_bin==i),ind_theta90);
end

figure;
plot(edges(1:4),bin, ':o') % PrefOri
xticks(edges(1:4))
hold on

for i = 1:length(tunedCells)
    tunedCells_bin(i) = length(cell2mat(tunedCells(i)));
end

plot(edges(1:4), tunedCells_bin, ':ok') % TunedCells
legend([{'prefOri'}, {'tunedCells'}])
title([mouse ': ' date])
ylim([0 inf])
xlabel('Oris')
ylabel('nCells')
%Check oris and binning for 210506

print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' dir_str], [date '_' mouse '_' dir_str '_prefOri_v_tuned.pdf']),'-dpdf','-fillpage')



save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' dir_str], [date '_' mouse '_' dir_str '_oriTuningInfo.mat']),...
    'prefOri', 'prefOri_bootdiff', 'ind_theta90', 'tunedCells', 'edges');




