close all;clear all;clc;

%%
base = findIsilon;
LG_base = fullfile(base,'home','lindsey');
behav_base = fullfile(base,'Behavior','Data');

mouse = 'i1373';
date = '221221';
ImgFolder = {'002'};
time = {'1429'};
nrun = length(ImgFolder);
frameRateHz = 15;
run_str = catRunName(cell2mat(ImgFolder), nrun);
        
outDir = fullfile(LG_base, 'Analysis', '2P', [date '_' mouse], [date '_' mouse '_' run_str]);

%% load
tic
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    fName = fullfile(behav_base, ['data-' mouse '-' date '-' time{irun} '.mat']);
    load(fName);
    
    data_base = fullfile(LG_base, 'Data', '2P_images', mouse, date, ImgFolder{irun});
    imgMatFile = [ImgFolder{irun} '_000_000.mat'];
    load(fullfile(data_base,imgMatFile));
    
    temp(irun) = input;
    nOn = temp(irun).nScansOn;
    nOff = temp(irun).nScansOff;
    ntrials = size(temp(irun).tGratingDirectionDeg,2);
    nframes = ntrials*(nOn+nOff);

    fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n'])
    data_temp = sbxread(fullfile(data_base, imgMatFile(1,1:11)),0,nframes);
    if size(data_temp,1)== 2
        data_temp = data_temp(1,:,:,:);
    end
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    fprintf('Complete')
end
input = concatenateDataBlocks(temp);
fprintf('\nAll runs read\n')
fprintf([num2str(size(data,3)) ' total frames\n'])
clear data_temp
clear temp

toc
        
%% Choose register interval
nskip = 5000;
nep = floor(size(data,3)./nskip);
[n n2] = subplotn(nep);
f= figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*nskip):500+((i-1)*nskip)),3)); title([num2str(1+((i-1)*nskip)) '-' num2str(500+((i-1)*nskip))]); colormap gray; clim([0 3000]); end
movegui('center')
f.WindowState= 'maximize';
data_avg = mean(data(:,:,45001:45500),3);
%% Register data
[out, data_reg] = stackRegister(data,data_avg);
figure; imagesc(mean(data_reg,3))
mkdir(fullfile(outDir))
print(fullfile(outDir, [date '_' mouse '_' run_str '_FOVavg.pdf']),'-dpdf')

%%

save(fullfile(outDir, [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
save(fullfile(outDir, [date '_' mouse '_' run_str '_input.mat']), 'input')
clear data out

%%
nOn = double(input.nScansOn);
nOff = double(input.nScansOff);
phaseCyc = double(input.nScansPhaseCyc);
dir_mat = celleqel2mat_padded(input.tGratingDirectionDeg);
nTrials = length(dir_mat);
sz = size(data_reg);

data_resp = reshape(data_reg, [sz(1) sz(2) nOn+nOff nTrials]);
data_f = mean(data_resp(:,:,nOff-15:nOff,:),3);
data_resp_dfof = (double(data_resp)-data_f)./data_f;

clear data_resp data_f

dirs = unique(dir_mat);
nDir = length(dirs);

tStimNum = celleqel2mat_padded(input.tStimulusNumber);
nPhase = input.gratingStartingPhaseStepN;
phase_mat = zeros(1,ntrials);
for itrial = 1:ntrials
    if tStimNum(itrial) < (nDir.*nPhase)
        temp = tStimNum(itrial);
    else
        temp = mod(tStimNum(itrial),nDir.*nPhase);
    end
    if temp < nPhase
        phase_mat(itrial) = input.gratingStartingPhaseDeg + (input.gratingStartingPhaseStepDeg.*temp);
    else
        phase_mat(itrial) = input.gratingStartingPhaseDeg + (input.gratingStartingPhaseStepDeg.*(mod(temp,input.gratingStartingPhaseStepN)));
    end
end
phases = unique(phase_mat);
mask_mat = celleqel2mat_padded(input.tMaskContrast);
masks = unique(mask_mat);
nmask = length(masks);

nStim = double(nPhase.*nDir*nmask);
[n n2] = subplotn(nStim);
data_dfof_stim = nan(sz(1),sz(2),nPhase,nDir,2);
start = 1;
for iPhase = 1:nPhase
    ind_p = find(phase_mat == phases(iPhase));
    for iDir = 1:nDir
        ind_d = find(dir_mat == dirs(iDir));
        for iM = 1:nmask
            ind_m = find(mask_mat == masks(iM));
            ind = intersect(ind_m,intersect(ind_p,ind_d));
            data_dfof_stim(:,:,iPhase,iDir,iM) = mean(mean(data_resp_dfof(:,:,nOff+1:nOff+nOn,ind),4),3);
            subplot(n,n2,start)
            imagesc(data_dfof_stim(:,:,iPhase,iDir,iM))
            start = start+1;
        end
    end
end
    
data_dfof_stim_all = reshape(data_dfof_stim, [sz(1) sz(2) nPhase.*nDir.*nmask]);
data_dfof_max = max(data_dfof_stim_all,[],3);
clear data_dfof_stim data_resp_dfof

figure; imagesc(data_dfof_max); movegui('center')

data_dfof = cat(3,data_dfof_max,data_dfof_stim_all);
save(fullfile(outDir, [date '_' mouse '_' run_str '_segmentData.mat']), 'data_dfof_max', 'data_dfof_stim_all', 'nStim')

%%
sz = size(data_dfof);
mask_data = data_dfof;
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
for iStim = 1:size(data_dfof,3)    
    mask_data_temp = mask_data(:,:,iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0; %blacks out old cells
    bwout = imCellEditInteractiveLG(mask_data_temp); %selection GUI
    mask_all = mask_all+bwout; %adds new cells to old cells
    mask_exp = imCellBuffer(mask_all,3)+mask_all; %creates buffer around cells to avoid fusing
    close all
end
mask_cell = bwlabel(mask_all); %turns logical into numbered cells
figure;
imagesc(mask_cell)
nMaskPix = 5; %thickness of neuropil ring in pixels
nBuffPix = 3; %thickness of buffer between cell and ring
mask_np = imCellNeuropil(mask_cell,nBuffPix,nMaskPix);

save(fullfile(outDir, [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_np')
clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_2 data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 

%% Neuropil subtraction
%Goal is to remove contamination from out-of-focus fluorescence
%Extract cell timecourses
data_tc = stackGetTimeCourses(data_reg, mask_cell); %applies mask to stack (averages all pixels in each frame for each cell) to get timecourses
            %Timecourses are nFrames x nCells
fprintf(['data_tc is ' num2str(size(data_tc))]) 
[nFrames, nCells] = size(data_tc);
%        1.  Downsampled timecourses for neuropil subtraction- averageing decreases noise
down = 5; %number of frames to average
data_reg_down = stackGroupProject(data_reg,down); %averages every 5 frames in stack  
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);    
%        2.  Extract neuropil timecourses (full and downsampled)
np_tc = zeros(nFrames,nCells);
np_tc_down = zeros(floor(nFrames./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '%s\n']) 
end
%        3. Find best neuropil weights by maximizing skew on downsampled subtractions
%            Assumes calcium signals are 1) sparse and 2) positive.  Too little subtraction will decrease sparseness and too much will make signals negative. 
%            Skewness decribes the shape of a distribution- a gaussian has a skew of 0; long tail to the right is a positive skew; long tail to the left if a negative skew.  
%            Thus, the best neuropil subtraction should maximize sparseness and minimize negative values.  
%            The best subtraction will therefore yield the highest skew for the distribution of fluorescence values for each cell.
%        a. Measure skew for all weights:
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
%        b. Find index with highest skew for each cell
[max_skew, ind] =  max(x,[],1); 
%        c. convert to weight
np_w = 0.01*ind; 
%        4. Subtract weighted neuropil response from full timecourses
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
clear data_reg data_reg_down data_tc_down np_tc_down

save(fullfile(outDir, [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

%% Phase reversal analysis
nOn = double(input.nScansOn);
nOff = double(input.nScansOff);
phaseCyc = double(input.nScansPhaseCyc);
ntrials = length(input.tGratingDirectionDeg);
nCells = size(npSub_tc,2);

cycPerTrial = floor(nOn/(phaseCyc*2));
data_trial = permute(reshape(npSub_tc,[nOn+nOff ntrials nCells]),[1 3 2]);
data_f = mean(data_trial(nOff./2:nOff,:,:),1);
data_dfof = (data_trial-data_f)./data_f;

data_dfof_cyc = zeros((phaseCyc.*2)+phaseCyc/2, nCells, ntrials, cycPerTrial);
for icyc = 1:cycPerTrial
    data_dfof_cyc(:,:,:,icyc) = data_dfof(nOff+(phaseCyc/2)+((icyc-1).*(phaseCyc.*2)):nOff+phaseCyc+(icyc.*(phaseCyc.*2))-1,:,:);
end
data_dfof_cycavg = mean(data_dfof_cyc,4);

dir_mat = celleqel2mat_padded(input.tGratingDirectionDeg);
dirs = unique(dir_mat);
nDir = length(dirs);

tStimNum = celleqel2mat_padded(input.tStimulusNumber);
nPhase = input.gratingStartingPhaseStepN;
phase_mat = zeros(1,ntrials);
for itrial = 1:ntrials
    if tStimNum(itrial) < (nDir.*nPhase)
        temp = tStimNum(itrial);
    else
        temp = mod(tStimNum(itrial),nDir.*nPhase);
    end
    if temp < nPhase
        phase_mat(itrial) = input.gratingStartingPhaseDeg + (input.gratingStartingPhaseStepDeg.*temp);
    else
        phase_mat(itrial) = input.gratingStartingPhaseDeg + (input.gratingStartingPhaseStepDeg.*(mod(temp,input.gratingStartingPhaseStepN)));
    end
end
phases = unique(phase_mat);

base_win = cell(1,2);
resp_win = cell(1,2);
base_win{1} = [phaseCyc/2-2:phaseCyc/2+3];
resp_win{1} = [(phaseCyc/2)+10:(phaseCyc/2)+15];
base_win{2} = [1.5*phaseCyc-2:1.5*phaseCyc+3];
resp_win{2} = [1.5*phaseCyc+10:1.5*phaseCyc+15];
tt = (1-phaseCyc/2:phaseCyc*2).*(1000/frameRateHz);

data_dfof_phasedir = zeros(2.5.*phaseCyc, nCells, nPhase, nDir, nmask);
h_dir = zeros(nCells,nPhase,nDir,nmask,2);
p_dir = nan(nCells,nPhase,nDir,nmask,2);
phasedir_resp_avg = zeros(nCells,nPhase,nDir,nmask,2,2);
trial_n = zeros(nPhase,nDir,nmask);
for iPhase = 1:nPhase
    ind_phase = find(phase_mat == phases(iPhase));
    for iDir = 1:nDir
        ind_dir = find(dir_mat == dirs(iDir));
        for iM = 1:nmask
            ind_m = find(mask_mat == masks(iM));
            ind = intersect(ind_m,intersect(ind_phase,ind_dir));
            trial_n(iPhase,iDir,iM) = length(ind);
            data_dfof_phasedir(:,:,iPhase,iDir,iM) = mean(data_dfof_cycavg(:,:,ind)-mean(data_dfof_cycavg(base_win{1},:,ind),1),3);
            for i = 1:2
                phasedir_resp_avg(:,iPhase,iDir,iM,i,1) = squeeze(mean((mean(data_dfof_cycavg(resp_win{i},:,ind),1)-mean(data_dfof_cycavg(base_win{i},:,ind),1)),3));
                phasedir_resp_avg(:,iPhase,iDir,iM,i,2) = squeeze(std((mean(data_dfof_cycavg(resp_win{i},:,ind),1)-mean(data_dfof_cycavg(base_win{i},:,ind),1)),[],3))./sqrt(length(ind));
                if length(ind)>3
                    [h_dir(:,iPhase,iDir,iM,i), p_dir(:,iPhase,iDir,iM,i)] = ttest2(squeeze(mean(data_dfof_cycavg(resp_win{i},:,ind),1)), squeeze(mean(data_dfof_cycavg(base_win{i},:,ind),1)),'dim', 2, 'tail', 'right', 'alpha', 0.05./(nDir-1));
                end
            end
        end
    end
end

[max_val, max_ind] = max(max(max(phasedir_resp_avg(:,:,:,1,:,1),[],5),[],2),[],3);

resp_ind_phase = find(sum(sum(sum(h_dir(:,:,:,1,:),2),3),5));

if length(resp_ind_phase)<64
[n n2] = subplotn(length(resp_ind_phase));
else
    n = 8;
    n2 = 8;
end
figure;
start = 1;
t = 1;
for iCell = 1:nCells
    if start>64
        suptitle([date ' ' mouse '- Phase reversal'])
        %print(fullfile(outDir,  [date '_' mouse '_' run_str '_avgRespByPhase_prefDir_' num2str(t) '.pdf']),'-dpdf','-fillpage');
        t = t+1;
        figure;
        start = 1;
    end
    subplot(n,n2,start)
    plot(tt, squeeze(data_dfof_phasedir(:,iCell,:,max_ind(iCell),1)))
    hold on
    vline(tt(resp_win{1}))
    vline(tt(resp_win{2}))
    %title([num2str(resp_ind_phase(iCell))])
    start = start+1;
end
suptitle([date ' ' mouse '- Phase reversal'])
print(fullfile(outDir,  [date '_' mouse '_' run_str '_avgRespByPhase_prefDir_' num2str(t) '.pdf']),'-dpdf','-fillpage');

% figure;
% [n n2] = subplotn(nCells);
% for iCell = 1:nCells
%     subplot(n,n2,iCell)
%     plot(tt, squeeze(data_dfof_phasedir(:,iCell,:,max_ind(iCell),1)))
%     hold on
%     vline(tt(resp_win{1}))
%     vline(tt(resp_win{2}))
%     title([num2str(iCell)])
% end
% suptitle([date ' ' mouse '- Phase reversal'])
% print(fullfile(outDir,  [date '_' mouse '_' run_str '_avgRespByPhase_prefDir_allCells.pdf']),'-dpdf','-fillpage');

save(fullfile(outDir,  [date '_' mouse '_' run_str '_dfofData.mat']), 'data_dfof', 'resp_win', 'base_win', 'tt', 'h_dir', 'p_dir', 'resp_ind_phase', 'data_dfof_phasedir', 'phasedir_resp_avg', 'max_ind', 'nCells')
save(fullfile(outDir,  [date '_' mouse '_' run_str '_stimData.mat']), 'mask_mat','masks', 'nmask','phase_mat', 'phases', 'nPhase', 'dir_mat', 'dirs', 'nDir', 'nOn', 'nOff','frameRateHz', 'phaseCyc')

%% measure F2/F1
f1_dir = zeros(nDir,nCells);
f2_dir = zeros(nDir,nCells);
f1 = zeros(1,nCells);
f2 = zeros(1,nCells);
for iCell = 1:nCells
    for iDir = 1:nDir
        cyc = squeeze(data_dfof_phasedir(phaseCyc./2+1:end,iCell,:,iDir,1))';
        [f1_dir(iDir,iCell),f2_dir(iDir,iCell),f1ang,projectedf1,f1mat,f2mat] = compcontrastrevf1f2(cyc);
    end
end
[max_dir_val max_dir_ind] = max(f1_dir,[],1);
orth_dir_ind = max_dir_ind+2;
orth_dir_ind(find(orth_dir_ind>4)) = orth_dir_ind(find(orth_dir_ind>4))-4;
f1 = indOnly(f1_dir',max_dir_ind');
f1_orth = indOnly(f1_dir',orth_dir_ind');
f2 = indOnly(f2_dir',max_dir_ind');
f2_orth = indOnly(f2_dir',orth_dir_ind');
f2overf1 = f2./f1;
f1_OSI = (f1-f1_orth)./(f1+f1_orth);
figure; scatter(f2overf1(find(f1>0.02)),f1_OSI(find(f1>0.02)))

figure;
subplot(2,2,1)
histogram(f1)
xlabel('F1')
subplot(2,2,2)
histogram(f2)
xlabel('F2')
subplot(2,2,3)
scatter(f1,f2)
vline(0.02)
xlabel('F1')
ylabel('F2')
xlim([0 0.3])
ylim([0 0.3])
refline(1)
subplot(2,2,4)
histogram(f2overf1(find(f1>0.02)))
xlabel('F2/F1')
suptitle([date ' ' mouse '- Phase reversal'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_f1f2.pdf']),'-dpdf','-bestfit');
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_f1f2.mat']), 'f1_dir', 'f2_dir', 'f1', 'f2', 'f2overf1', 'max_dir_ind', 'f1_orth','f2_orth','f1_OSI', 'resp_ind_phase')

%%
if input.doMask
    f1_mask_dir = zeros(1,nCells);
    f2_mask_dir = zeros(1,nCells);
    f1_mask = zeros(1,nCells);
    f2_mask = zeros(1,nCells);
    for iCell = 1:nCells
        for iDir = 1:nDir
            cyc = squeeze(data_dfof_phasedir(phaseCyc./2+1:end,iCell,:,iDir,2))';
            [f1_mask_dir(iDir,iCell),f2_mask_dir(iDir,iCell),f1ang,projectedf1,f1mat,f2mat] = compcontrastrevf1f2(cyc);
        end
    end

    mask_dir_ind = max_dir_ind+2;
    mask_dir_ind(find(mask_dir_ind>4)) = mask_dir_ind(find(mask_dir_ind>4))-4;

    f1_mask = indOnly(f1_mask_dir',mask_dir_ind');
    f2_mask = indOnly(f2_mask_dir',mask_dir_ind');
    f2overf1_mask = f2_mask./f1_mask;

    f1_mask90 = indOnly(f1_mask_dir',max_dir_ind');
    f2_mask90 = indOnly(f2_mask_dir',max_dir_ind');
    f2overf1_mask90 = f2_mask90./f1_mask90;

    for iCell = x
    figure;
    start = 1;
    
    for iDir = 1:4
        subplot(4,2,start)
        plot(tt, squeeze(data_dfof_phasedir(:,iCell,:,iDir,1)))
        subplot(4,2,start+1)
        plot(tt, squeeze(data_dfof_phasedir(:,iCell,:,iDir,2)))
        start = start+2;
    end
    suptitle([num2str(iCell) ' ' num2str(chop(f1(iCell),3))])
    end

    figure;
    subplot(2,2,1)
    histogram(f1_mask(find(f1>0.02)))
    xlabel('F1')
    subplot(2,2,2)
    histogram(f2_mask(find(f1>0.02)))
    xlabel('F2')
    subplot(2,2,3)
    scatter(f1_mask(find(f1>0.02)),f2_mask(find(f1>0.02)))
    xlabel('F1')
    ylabel('F2')
    xlim([0 0.3])
    ylim([0 0.3])
    refline(1)
    subplot(2,2,4)
    histogram(f2overf1_mask(find(f1>0.02)))
    xlabel('F2/F1')
    suptitle([date ' ' mouse '- Mask Phase reversal'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_f1f2_mask.pdf']),'-dpdf','-bestfit');

    figure;
    subplot(2,2,1)
    histogram(f1_mask90(find(f1>0.02)))
    xlabel('F1')
    subplot(2,2,2)
    histogram(f2_mask90(find(f1>0.02)))
    xlabel('F2')
    subplot(2,2,3)
    scatter(f1_mask90(find(f1>0.02)),f2_mask90(find(f1>0.02)))
    xlabel('F1')
    ylabel('F2')
    xlim([0 0.3])
    ylim([0 0.3])
    refline(1)
    subplot(2,2,4)
    histogram(f2overf1_mask90(find(f1>0.02)))
    xlabel('F2/F1')
    suptitle([date ' ' mouse '- Mask Phase 90 reversal'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_f1f2_mask90.pdf']),'-dpdf','-bestfit');


    cell_use = intersect(find(f1>0.02),find(f1_mask>0.02));
    figure; 
    subplot(2,3,1)
    scatter(f1(cell_use),f1_mask(cell_use))
    xlim([0 0.4])
    ylim([0 0.4])
    refline(1)
    title('F1')
    xlabel('grating')
    ylabel('plaid')
    subplot(2,3,2)
    scatter(f2(cell_use),f2_mask(cell_use))
    xlim([0 0.4])
    ylim([0 0.4])
    refline(1)
    title('F2')
    xlabel('grating')
    ylabel('plaid')
    subplot(2,3,3)
    scatter(f2overf1(cell_use),f2overf1_mask(cell_use))
    xlim([0 3])
    ylim([0 3])
    refline(1)
    title('F2/F1')
    xlabel('grating')
    ylabel('plaid')
    subplot(2,3,4)
    scatter(f1(cell_use),f1_mask90(cell_use))
    xlim([0 0.4])
    ylim([0 0.4])
    refline(1)
    title('F1')
    xlabel('grating')
    ylabel('plaid 90')
    subplot(2,3,5)
    scatter(f2(cell_use),f2_mask90(cell_use))
    xlim([0 0.4])
    ylim([0 0.4])
    refline(1)
    title('F2')
    xlabel('grating')
    ylabel('plaid 90')
    subplot(2,3,6)
    scatter(f2overf1(cell_use),f2overf1_mask90(cell_use))
    xlim([0 3])
    ylim([0 3])
    refline(1)
    title('F2/F1')
    xlabel('grating')
    ylabel('plaid 90')
    suptitle([date ' ' mouse '- Phase reversal'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_f1f2_CompareTestMask.pdf']),'-dpdf','-bestfit');
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_f1f2mask.mat']), 'f1_mask_dir', 'f2_mask_dir', 'f1_mask', 'f2_mask', 'f2overf1_mask','f1_mask90', 'f2_mask90', 'f2overf1_mask90', 'mask_dir_ind')
end