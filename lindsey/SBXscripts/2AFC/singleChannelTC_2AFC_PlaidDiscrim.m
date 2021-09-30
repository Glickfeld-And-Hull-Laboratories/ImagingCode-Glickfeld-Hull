%%
close all
clear all
clc
dataset = 'plaidDiscrim_exptList';
eval(dataset);
iexp = 1;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
if strcmp(expt(iexp).folder,'lindsey')
    data_base = LG_base;
elseif strcmp(expt(iexp).folder,'camaron')
    data_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron';
end
mouse = expt(iexp).mouse;
date = expt(iexp).date;
%% get path names
tic
data = [];
clear temp
offset = 0;

%load data
nrun = size(expt(iexp).runs,1);
for irun = 1:nrun
    CD = [data_base '\Data\2P_images\' expt(iexp).mouse '\' expt(iexp).date '\' expt(iexp).runs(irun,:)];
    cd(CD);
    imgMatFile = [expt(iexp).runs(irun,:) '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' expt(iexp).mouse '-' expt(iexp).date '-' expt(iexp).time_mat(irun,:) '.mat'];
    load(fName);
    
    nframes = [input.counterValues{end}(end) info.config.frames];
    
    if min(nframes)<input.counterValues{end}(end)
        ntrials = size(input.trialOutcomeCell,2);
        for itrial = ntrials:-1:1
            if input.counterValues{itrial}(end) <= nframes
                break
            end
        end
        input = trialChopper(input,[1 itrial]);
    end
    temp(irun) = input;
    fprintf(['Reading run ' num2str(irun) '- ' num2str(min(nframes)) ' frames \r\n'])
    data_temp = sbxread(imgMatFile(1,1:11),0,min(nframes));
    if size(data_temp,1)== 2
        data_temp = data_temp(1,:,:,:);
    end
    
    if irun>1
        if isfield(input, 'tLeftTrial')
            ntrials = size(input.trialOutcomeCell,2);
            for itrial = 1:ntrials
                temp(irun).counterValues{itrial} = bsxfun(@plus,temp(irun).counterValues{itrial},offset);
                temp(irun).cTrialStart{itrial} = temp(irun).cTrialStart{itrial}+offset;
                temp(irun).cAdaptOn{itrial} = temp(irun).cAdaptOn{itrial}+offset;
                temp(irun).cStimOn{itrial} = temp(irun).cStimOn{itrial}+offset;
                temp(irun).cDecision{itrial} = temp(irun).cDecision{itrial}+offset;
            end
        end
    end
    
    offset = offset+min(nframes);
        
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
end 

input = concatenateDataBlocks(temp);
clear data_temp
clear temp

fprintf('Runs loaded\n')
%% Plot outcome by trial number
SIx = strcmp(input.trialOutcomeCell, 'success');
MIx = strcmp(input.trialOutcomeCell, 'ignore');

figure;
plot(smooth(SIx,10));
hold on
plot(smooth(MIx,10));
hold on
plot(celleqel2mat_padded(input.tDoFeedbackMotion))

%% Choose register interval
nep = floor(size(data,3)./10000);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end

%% Register data
run_str = ['runs']; 
for irun = 1:nrun
    run_str = [run_str '-' expt(iexp).runs(irun,:)];
end

if exist(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    load(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    [out, data_reg] = stackRegister_MA(data,[],[],out); 
else
    data_avg = mean(data(:,:,40001:40500),3);
    [out, data_reg] = stackRegister(data,data_avg);
    smooth_out_x = smooth(out(:,3),200);
    smooth_out_y = smooth(out(:,4),200);
    diff_out = max([abs(smooth_out_x-out(:,3)) abs(smooth_out_y-out(:,4))], [],2);
    move_ind = find(diff_out>10);
    mkdir(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg', 'diff_out', 'move_ind')
    save(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
clear data

%% test stability

figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end
figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')
%%
if ~isempty(expt(iexp).redImg)
    if exist(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_redData.mat']))
        load(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_redData.mat']))
    else
        CD = [data_base '\Data\2P_images\' expt(iexp).mouse '\' expt(iexp).date '\' expt(iexp).redImg{1}];
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
                CD = [data_base '\Data\2P_images\' expt(iexp).mouse '\' expt(iexp).date '\' expt(iexp).redImg{2}];
                cd(CD);
                imgMatFile = [expt(iexp).redImg{2} '_000_000.mat'];
                load(imgMatFile);
                nframes = info.config.frames;
                fprintf(['Reading run ' expt(iexp).redImg{2} '- ' num2str(min(nframes)) ' frames \r\n'])
                data = sbxread(imgMatFile(1,1:11),0,nframes);
                red_data = squeeze(data(2,:,:,:));
                green_data = squeeze(data(1,:,:,:));
                [out, green_data_reg] = stackRegister(mean(green_data,3),green_data_avg);
                [outs, red_data_avg] = stackRegister_MA(mean(red_data,3),[],[],out);
                figure; imagesc(red_data_avg)
                title('Red')
            end
            save(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_redData.mat']), 'green_data_avg', 'red_data_avg')
        end
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
%     rgb(:,:,1) = green_data_avg./max(green_data_avg(:));
%     rgb(:,:,2) = data_avg./max(data_avg(:));
%     subplot(2,2,3); imagesc(rgb);
%     title('Green at 1040 and 920')
%     rgb(:,:,1) = red_data_avg./max(red_data_avg(:));
%     rgb(:,:,2) = green_data_avg./max(green_data_avg(:));
%     subplot(2,2,4); imagesc(rgb);
%     title('Red and Green at 1040')
    print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

    clear red_data red_data_reg green_data green_data_reg data

end 
%% counter check
ntrials = size(input.tGratingContrast,2);
counterVals = [];
counterTimes = [];
cStart = [];
cEnd = [];
for itrial = 1:ntrials
    counterTimes = [counterTimes input.counterTimesUs{itrial}./1000];
    counterVals = [counterVals input.counterValues{itrial}];
    cStart = [cStart input.counterValues{itrial}(1)];
    cEnd = [cEnd input.counterValues{itrial}(end)];
end

SIx = strcmp(input.trialOutcomeCell,'success');
dCount = diff(counterTimes);
dVal = diff(counterVals);
figure; plot(dCount); ylim([0 70]); vline(cEnd(SIx),':r'); 
% hold on; plot(dVal.*30)
% 
% short_ind = find(dCount<20);
% short_after_long_ind = find(dCount(short_ind-1)>40);
% short_ind(short_after_long_ind) = [];
% counterVals_fixed = counterVals;
% cStimOn = celleqel2mat_padded(input.cStimOn);
% cAdaptOn = celleqel2mat_padded(input.cAdaptOn);
% cDecision = celleqel2mat_padded(input.cDecision);
% 
% for i = 1:length(short_ind)
%     val = short_ind(i)+1;
%     counterVals_fixed(val:end) = counterVals_fixed(val:end)-1;
%     cStimOn(find(cStimOn>=val)) = cStimOn(find(cStimOn>=val))-1;
%     cAdaptOn(find(cAdaptOn>=val)) = cAdaptOn(find(cAdaptOn>=val))-1;
%     cDecision(find(cDecision>=val)) = cDecision(find(cDecision>=val))-1;
% end
% counterTimes_fixed = counterTimes;
% counterTimes_fixed(short_ind) = [];
% counterVals_fixed(short_ind) = [];
cAdaptOn = celleqel2mat_padded(input.cAdaptOn);
cStimOn = celleqel2mat_padded(input.cStimOn);

%% 2AFC photodiode check
if isfield(info,'frame')
    [stimOnFrames stimOffFrames] = photoFrameFinder_Sanworks(info.frame);
    save(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_photoData.mat']), 'stimOnFrames', 'stimOffFrames')
end
%% find activated cells
close all
tGratingDir = chop(celleqel2mat_padded(input.tGratingDirectionStart),3);
tPlaidDir = chop(celleqel2mat_padded(input.tPlaidDirectionStart),3);
tMaskTrial = celleqel2mat_padded(input.tMaskContrast)>0;
tDir = tGratingDir;
tDir(find(tMaskTrial)) = tPlaidDir(find(tMaskTrial));
dirs = unique(tDir);
nDir = length(dirs);

nTrials = length(tGratingDir);
sz = size(data_reg);
data_f = nan(sz(1),sz(2),nTrials);
data_targ = nan(sz(1),sz(2),nTrials);
data_tc = nan(sz(1),sz(2),21,nTrials);
for itrial = 1:nTrials
    if cStimOn(itrial)+15<sz(3)
        data_f(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)-frameRateHz:cStimOn(itrial)-1),3);
        data_targ(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)+5:cStimOn(itrial)+15),3);
        data_tc(:,:,:,itrial) = data_reg(:,:,cStimOn(itrial)-5:cStimOn(itrial)+15);
    end
end
data_targ_dfof = (data_targ-data_f)./data_f;
data_tc_dfof = mean((data_tc-mean(data_f,3))./mean(data_f,3),4);
figure;
for i = 1:21
    subplot(5,5,i)
    imagesc(mean(data_tc(:,:,i,:),4));
end

[n n2] = subplotn(nDir*2);
figure;
data_dfof = zeros(sz(1),sz(2),nDir*2);
grating_type = {'Grating','Plaid'};
start = 1;
for iDir = 1:nDir
    for i = 1:2
        subplot(n,n2,start)
        ind = intersect(find(tMaskTrial==i-1),find(tDir == dirs(iDir)));
        data_dfof(:,:,iDir)= nanmean(data_targ_dfof(:,:,ind),3);
        imagesc(data_dfof(:,:,iDir));
        title([grating_type{i} '- ' num2str(dirs(iDir)) ' deg'])
        start = start+1;
    end
end

data_dfof_max = max(data_dfof,[],3);
figure; imagesc(data_dfof_max)
data_dfof = cat(3,data_dfof_max,data_dfof);


%% cell segmentation  
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
    print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_bouton_mask.pdf']), '-dpdf')

    mask_np = imCellNeuropil(mask_cell, 3, 5);
    if ~isempty(expt(iexp).redImg)
        save(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_cell_red', 'mask_np')
    else
        save(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_np')  
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
                        h = zeros(1, 2+nGratingDir+nGratingDir);
                        dirs = [0 90];
                        start = 0;
                        for idir = 1:2
                            ind = find(aContrast == 1 & b2Ix == idir-1);
                            [h(1,idir) p] = ttest(squeeze(mean(mean(data_adapt_dfof(i-1:i+1,j-1:j+1,ind),1),2)),0,'tail','right','alpha',0.05./(1+nGratingDir+nGratingDir));
                        end
                        start = start+2;
                        for iOri = 1:nGratingDir
                            ind = find(tGratingDir == dirs(iDir));
                            [h(1,start+iOri) p] = ttest(squeeze(mean(mean(data_targ_dfof(i-1:i+1,j-1:j+1,ind),1),2)),0,'tail','right','alpha',0.05./(1+nGratingDir+nGratingDir));
                        end
                        start = start+nGratingDir;
                        for iDir = 1:nGratingDir
                            ind = find(dir_mat == dirs(iDir));
                            [h(1,start+iDir) p] = ttest(squeeze(mean(mean(dir_dfof(i-1:i+1,j-1:j+1,ind),1),2)),0,'tail','right','alpha',0.05./(1+nGratingDir+nGratingDir));
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
    print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_bouton_mask.pdf']), '-dpdf')

    save(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_max', 'mask_cell') 
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
    save(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc','np_tc','npSub_tc')
elseif strcmp(cell2mat(expt(iexp).img_strct),'axons')
    npSub_tc = data_tc;
    save(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'npSub_tc','-v7.3')
end
clear data_reg data_reg_down
clear data_tc np_tc

%% 2AFC analysis
cDecision = cStimOn + (celleqel2mat_padded(input.cDecision)-celleqel2mat_padded(input.cStimOn));
nCells = size(npSub_tc,2);
nframes = size(npSub_tc,1);
nTrials = size(cDecision,2);
data_stim = nan(frameRateHz*2,nCells,nTrials);
data_dec = nan(frameRateHz*2,nCells,nTrials);
for itrial = 1:nTrials
    if cStimOn(itrial)+frameRateHz-1 < nframes
        data_stim(:,:,itrial) = npSub_tc(cStimOn(itrial)-frameRateHz:cStimOn(itrial)+frameRateHz-1,:);
    end
    if ~isnan(cDecision(itrial))
        if cDecision(itrial)+29 < nframes
            data_dec(:,:,itrial) = npSub_tc(cDecision(itrial)-frameRateHz:cDecision(itrial)+frameRateHz-1,:);
        end
    end
end
dataf = mean(data_stim(1:frameRateHz,:,:),1);
data_stim_dfof = bsxfun(@rdivide, bsxfun(@minus, data_stim, dataf), dataf);
data_dec_dfof = bsxfun(@rdivide, bsxfun(@minus, data_dec, dataf), dataf);
tt = [-frameRateHz:frameRateHz-1].*(1000./frameRateHz);
figure;
subplot(1,2,1)
plot(tt,nanmean(mean(data_stim_dfof,2),3));
title('Target')
subplot(1,2,2)
plot(tt,nanmean(mean(data_dec_dfof,2),3));
title('Decision')


figure;
for iDir = 1:nDir
    ind = intersect(find(tMaskTrial == 0), find(tDir == dirs(iDir)));
    subplot(2,nDir,iDir)
    plot(tt, mean(nanmean(data_stim_dfof(:,:,ind),3),2))
    title([num2str(dirs(iDir))])
    ylim([-0.05 0.1])
    ind = intersect(find(tMaskTrial == 1), find(tDir == dirs(iDir)));
    subplot(2,nDir,iDir+nDir)
    plot(tt, mean(nanmean(data_stim_dfof(:,:,ind),3),2))
    title([num2str(dirs(iDir))])
    ylim([-0.05 0.1])
end
print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_avgTargetResponse_byTarget.pdf']), '-dpdf','-bestfit')

SIx = strcmp(input.trialOutcomeCell,'success');
MIx = strcmp(input.trialOutcomeCell,'incorrect');
IIx = strcmp(input.trialOutcomeCell,'ignore');
tLeftTrial = celleqel2mat_padded(input.tLeftTrial);
tLeftResp = zeros(size(SIx));
tLeftResp(find(tLeftTrial&SIx)) = 1;
tLeftResp(find(~tLeftTrial&~SIx)) = 1;

figure;
for iDir = 1:nDir
    ind = intersect(find(tMaskTrial == 0 & SIx), find(tDir == dirs(iDir)));
    subplot(2,nDir,iDir)
    plot(tt, mean(nanmean(data_stim_dfof(:,:,ind),3),2))
    hold on
    ind = intersect(find(tMaskTrial == 0 & MIx), find(tDir == dirs(iDir)));
    plot(tt, mean(nanmean(data_stim_dfof(:,:,ind),3),2))
    title([num2str(dirs(iDir))])
    ylim([-0.05 0.1])
    ind = intersect(find(tMaskTrial == 1 & SIx), find(tDir == dirs(iDir)));
    subplot(2,nDir,iDir+nDir)
    plot(tt, mean(nanmean(data_stim_dfof(:,:,ind),3),2))
    hold on
    ind = intersect(find(tMaskTrial == 1 & MIx), find(tDir == dirs(iDir)));
    plot(tt, mean(nanmean(data_stim_dfof(:,:,ind),3),2))
    title([num2str(dirs(iDir))])
    ylim([-0.05 0.1])
end
legend({'Hit','Miss'})
print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_avgTargetResponse_byTarget_byOutcome.pdf']), '-dpdf','-bestfit')


mask_label = zeros(1,nCells);
for i = 1:nCells
    if mask_cell_red(find(mask_cell == i, 1))
        mask_label(1,i) = 1;
    end
end

interval = ceil(64/nDir)+1;
x = 1:interval:64;
x(end) = 63;
y = bluered;

figure;
if nCells <49
    [n n2] = subplotn(nCells);
else
    [n n2] = subplotn(49);
end
figure;
start = 1;
for iC = 1:nCells
    if start > 49
        suptitle('Target response- Grating- All trials')
        print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_targetRespGrating_cells' num2str(iC-49) '-' num2str(iC-1) '.pdf']), '-dpdf','-bestfit')
        figure;
        start = 1;
    end
    subplot(n,n2,start)
    for iDir = 1:nDir
        ind_dir = find(tDir == dirs(iDir) & tMaskTrial==0);
        plot(tt, mean(data_stim_dfof(:,iC,ind_dir),3),'Color',y(x(iDir),:))
        hold on
    end
    if mask_label(iC)
        title([num2str(iC) '-' expt(iexp).driver])
    else
        title(num2str(iC))
    end
    start = start+1;
end
suptitle('Target response- Grating- All trials')
print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_targetRespGrating_cells' num2str(iC-48) '-' num2str(iC) '.pdf']), '-dpdf','-bestfit')

figure;
start = 1;
for iC = 1:nCells
    if start > 49
        suptitle('Target response- Plaid- All trials')
        print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_targetRespPlaid_cells' num2str(iC-49) '-' num2str(iC-1) '.pdf']), '-dpdf','-bestfit')
        figure;
        start = 1;
    end
    subplot(n,n2,start)
    for iDir = 1:nDir
        ind_dir = find(tDir == dirs(iDir) & tMaskTrial==1);
        plot(tt, mean(data_stim_dfof(:,iC,ind_dir),3),'Color',y(x(iDir),:))
        hold on
    end
    if mask_label(iC)
        title([num2str(iC) '-' expt(iexp).driver])
    else
        title(num2str(iC))
    end
    start = start+1;
end
suptitle('Target response- Plaid- All trials')
print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_targetRespPlaid_cells' num2str(iC-48) '-' num2str(iC) '.pdf']), '-dpdf','-bestfit')


figure;
start = 1;
for iC = 1:nCells
    if start > 49
        suptitle('Target grating response by choice- Trials <20 deg')
        print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_targetRespGrating_bychoice_cells' num2str(iC-49) '-' num2str(iC-1) '.pdf']), '-dpdf','-bestfit')
        figure;
        start = 1;
    end
    subplot(n,n2,start)
    plot(tt, mean(data_stim_dfof(:,iC,find(~tMaskTrial&~tLeftResp&~IIx&abs(tDir-45)<20)),3))
    hold on
    plot(tt, mean(data_stim_dfof(:,iC,find(~tMaskTrial&tLeftResp&~IIx&abs(tDir-45)<20)),3))
    if mask_label(iC)
        title([num2str(iC) '-' expt(iexp).driver])
    else
        title(num2str(iC))
    end
    start = start+1;
end
suptitle('Target grating response by choice- Trials <20 deg')
print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_targetRespGrating_bychoice_cells' num2str(iC-48) '-' num2str(iC) '.pdf']), '-dpdf','-bestfit')

figure;
start = 1;
for iC = 1:nCells
    if start > 49
        suptitle('Target plaid response by choice- Trials <20 deg')
        print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_targetRespPlaid_bychoice_cells' num2str(iC-49) '-' num2str(iC-1) '.pdf']), '-dpdf','-bestfit')
        figure;
        start = 1;
    end
    subplot(n,n2,start)
    plot(tt, mean(data_stim_dfof(:,iC,find(tMaskTrial&~tLeftResp&~IIx&abs(tDir-45)<20)),3))
    hold on
    plot(tt, mean(data_stim_dfof(:,iC,find(tMaskTrial&tLeftResp&~IIx&abs(tDir-45)<20)),3))
    if mask_label(iC)
        title([num2str(iC) '-' expt(iexp).driver])
    else
        title(num2str(iC))
    end
    start = start+1;
end
suptitle('Target plaid response by choice- Trials <20 deg')
print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_targetRespPlaid_bychoice_cells' num2str(iC-48) '-' num2str(iC) '.pdf']), '-dpdf','-bestfit')


figure;
start = 1;
for iC = 1:nCells
    if start > 49
        suptitle('Decision response- All trials')
        print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_decResp_allCells_success_cells' num2str(iC-49) '-' num2str(iC-1) '.pdf']), '-dpdf','-bestfit')
        figure;
        start = 1;
    end
    subplot(n,n2,start)
    for iDir = 1:nDir
        ind_dir = find(SIx & tDir == dirs(iDir));
        plot(tt, mean(data_dec_dfof(:,iC,ind_dir),3),'Color',y(x(iDir),:))
        hold on
    end
    if mask_label(iC)
        title([num2str(iC) '-' expt(iexp).driver])
    else
        title(num2str(iC))
    end
    start = start+1;
end
suptitle('Decision response- All trials')
print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_decResp_allCells_success_cells' num2str(iC-48) '-' num2str(iC) '.pdf']), '-dpdf','-bestfit')

figure;
subplot(3,2,1) 
for iDir = 1:nDir
    ind_dir = find(~tMaskTrial & tDir == dirs(iDir));
    plot(tt, mean(mean(data_stim_dfof(:,find(mask_label),ind_dir),3),2),'Color',y(x(iDir),:))
    hold on
end
xlabel('Time from target grating (ms)')
ylim([-0.02 0.15])
title([expt(iexp).driver '+: n = ' num2str(sum(mask_label))])
subplot(3,2,2)
for iDir = 1:nDir
    ind_dir = find(~tMaskTrial & tDir == dirs(iDir));
    plot(tt, mean(mean(data_stim_dfof(:,find(~mask_label),ind_dir),3),2),'Color',y(x(iDir),:))
    hold on
end
xlabel('Time from target grating (ms)')
ylim([-0.02 0.15])
title(['Pyr: n = ' num2str(sum(mask_label))])
subplot(3,2,3) 
for iDir = 1:nDir
    ind_dir = find(tMaskTrial & tDir == dirs(iDir));
    plot(tt, mean(mean(data_stim_dfof(:,find(mask_label),ind_dir),3),2),'Color',y(x(iDir),:))
    hold on
end
xlabel('Time from target plaid (ms)')
ylim([-0.02 0.15])
subplot(3,2,4)
for iDir = 1:nDir
    ind_dir = find(tMaskTrial & tDir == dirs(iDir));
    plot(tt, mean(mean(data_stim_dfof(:,find(~mask_label),ind_dir),3),2),'Color',y(x(iDir),:))
    hold on
end
xlabel('Time from target plaid (ms)')
ylim([-0.02 0.15])
subplot(3,2,5)
for iDir = 1:nDir
    ind_dir = find(SIx & tDir == dirs(iDir));
    plot(tt, mean(mean(data_dec_dfof(:,find(mask_label),ind_dir),3),2),'Color',y(x(iDir),:))
    hold on
end
xlabel('Time from decision (ms)')
ylim([-0.02 0.15])
subplot(3,2,6)
for iDir = 1:nDir
    ind_dir = find(SIx & tDir == dirs(iDir));
    plot(tt, mean(mean(data_dec_dfof(:,find(~mask_label),ind_dir),3),2),'Color',y(x(iDir),:))
    hold on
end
xlabel('Time from decision (ms)')
ylim([-0.02 0.15])
print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_avgResp_allCells.pdf']), '-dpdf','-bestfit')

%%
resp_win = frameRateHz+5:frameRateHz*2;
data_stim_resp = squeeze(mean(data_stim_dfof(resp_win,:,:),1));
stim_resp_dir = zeros(nCells,nDir,2,2);
for iDir = 1:nDir
    ind_dir = find(tDir == dirs(iDir) & tMaskTrial==0);
    stim_resp_dir(:,iDir,1,1) = nanmean(data_stim_resp(:,ind_dir),2);
    stim_resp_dir(:,iDir,1,2) = nanstd(data_stim_resp(:,ind_dir),[],2)./sqrt(length(ind_dir));
    ind_dir = find(tDir == dirs(iDir) & tMaskTrial);
    stim_resp_dir(:,iDir,2,1) = nanmean(data_stim_resp(:,ind_dir),2);
    stim_resp_dir(:,iDir,2,2) = nanstd(data_stim_resp(:,ind_dir),[],2)./sqrt(length(ind_dir));
end

figure;
start = 1;
for iC = 1:nCells
    if start > 49
        suptitle('Target response- All trials')
        print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_targetRespTuning_cells' num2str(iC-49) '-' num2str(iC-1) '.pdf']), '-dpdf','-bestfit')
        figure;
        start = 1;
    end
    subplot(n,n2,start)
    errorbar(dirs, stim_resp_dir(iC,:,1,1),stim_resp_dir(iC,:,1,2),'-o')
    hold on
    errorbar(dirs, stim_resp_dir(iC,:,2,1),stim_resp_dir(iC,:,2,2),'-o')
    ylim([-0.05 0.1])
    if mask_label(iC)
        title([num2str(iC) '-' expt(iexp).driver])
    else
        title(num2str(iC))
    end
    start = start+1;
end
suptitle('Target response- All trials')