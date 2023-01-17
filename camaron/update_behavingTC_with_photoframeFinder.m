load('Z:\All_Staff\home\camaron\Analysis\2P\behavior_pFF_days_to_fix.mat')

%% Load experiment info
dataset = 'oriAdapt_V1_cam';
eval(dataset); % run file to load expt.structure

for iday = 12:length(days_to_fix)

%%
clearvars -except expt iday days_to_fix
clc
close all

disp([num2str(iday) ' of ' num2str(length(days_to_fix))])


% %% Load experiment info
% dataset = 'oriAdapt_V1_cam';
% eval(dataset); % run file to load expt.structure

iexp = days_to_fix(iday); % Enter experiment number from oriAdapt_V1

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
CM_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron';

if strcmp(expt(iexp).folder,'lindsey')
    data_base = LG_base;
elseif strcmp(expt(iexp).folder,'camaron')
    data_base = CM_base;
end

mouse = expt(iexp).mouse;
date = expt(iexp).date;

%% Get imaging data (data, nframes, adapt_input)
% tic
% data = [];
% clear temp
% offset = 0;
% 
% 
% 
% 
% % load data
% nrun = size(expt(iexp).runs,1);
% for irun = 1:nrun
%     CD = [data_base '\Data\2P_images\' expt(iexp).mouse '\' expt(iexp).date '\' expt(iexp).runs(irun,:)];
%     cd(CD);
%     imgMatFile = [expt(iexp).runs(irun,:) expt(iexp).runs_suffix(irun,:) '.mat']; % DONE; Make variable to pull for oriAdapt_V1 that points to imgMatFile of restarted runs (ex: 001_000_001)
%     load(imgMatFile);
%     fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' expt(iexp).mouse '-' expt(iexp).date '-' expt(iexp).time_mat(irun,:) '.mat'];
%     load(fName);
%     
%     nframes = [input.counterValues{end}(end) info.config.frames];
% 
%   
%     
%     if min(nframes)<input.counterValues{end}(end)
%         ntrials = size(input.trialOutcomeCell,2);
%         for itrial = ntrials:-1:1
%             if input.counterValues{itrial}(end) <= nframes
%                 break
%             end
%         end
%         input = trialChopper(input,[1 itrial]);
%     end
%     temp(irun) = input;
% 
%     fprintf(['Reading run ' num2str(irun) '- ' num2str(min(nframes)) ' frames \r\n'])
%     data_temp = sbxread(imgMatFile(1,1:11),0,min(nframes));
%     if size(data_temp,1)== 2
%         data_temp = data_temp(1,:,:,:);
%     end
%     
%     if irun>1
%         if isfield(input, 'tLeftTrial')
%             ntrials = size(input.trialOutcomeCell,2);
%             for itrial = 1:ntrials
%                 temp(irun).counterValues{itrial} = bsxfun(@plus,temp(irun).counterValues{itrial},offset);
%                 temp(irun).cTrialStart{itrial} = temp(irun).cTrialStart{itrial}+offset;
%                 temp(irun).cAdaptOn{itrial} = temp(irun).cAdaptOn{itrial}+offset;
%                 temp(irun).cStimOn{itrial} = temp(irun).cStimOn{itrial}+offset;
%                 temp(irun).cDecision{itrial} = temp(irun).cDecision{itrial}+offset;
%             end
%         end
%     end
%     
%     offset = offset+min(nframes);
%         
%     data_temp = squeeze(data_temp);
%     data = cat(3,data,data_temp);
% end 
% 
% adapt_input = concatenateDataBlocks(temp);
% clear data_temp
% clear temp
% fprintf('Runs loaded\n')
% toc
% 
% 
% %% Choose register interval
% stack_spacing = 10000; % Number of frames in each image stack
% nep = floor(size(data,3)./stack_spacing);
% [n n2] = subplotn(nep);
% figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*stack_spacing):500+((i-1)*stack_spacing)),3)); title([num2str(1+((i-1)*stack_spacing)) '-' num2str(500+((i-1)*stack_spacing))]); end
% 
% %% Register green data (data_avg, data_reg, move_ind)
% run_str = ['runs']; 
% nrun = size(expt(iexp).runs,1);
% for irun = 1:nrun
%     run_str = [run_str '-' expt(iexp).runs(irun,:)];
% end
% 
% if exist(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
%     load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
%     [out, data_reg] = stackRegister_MA(data,[],[],out); 
% else
%     registration_interval = 50001:50500; % as determined above
%     data_avg = mean(data(:,:,registration_interval),3); 
%     [out, data_reg] = stackRegister(data,data_avg);
%     smooth_out_x = smooth(out(:,3),200);
%     smooth_out_y = smooth(out(:,4),200);
%     diff_out = max([abs(smooth_out_x-out(:,3)) abs(smooth_out_y-out(:,4))], [],2);
%     move_ind = find(diff_out>10); % keep
%     mkdir(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str]))
%     save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg', 'diff_out', 'move_ind')
%     save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'adapt_input')
% end
% clearvars data smooth_out_x smooth_out_y;
% 
% %% Check image for stability
% 
% figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*stack_spacing):500+((i-1)*stack_spacing)),3)); title([num2str(1+((i-1)*stack_spacing)) '-' num2str(500+((i-1)*stack_spacing))]); end
% figure; imagesq(mean(data_reg(:,:,1:stack_spacing),3)); truesize;
% print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')
% 
% clearvars nep n n2
% 
% %% If unstable, use different registration pathway...
% 
% % CODE GOES HERE
% 
% %% Register Red to green (green_data_avg, red_data_avg)
% 
% if ~isempty(expt(iexp).redImg)
%     if exist(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_redData.mat']))
%         load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_redData.mat']))
%         else
%             CD = [data_base '\Data\2P_images\' expt(iexp).mouse '\' expt(iexp).date '\' expt(iexp).redImg{1}];
%             cd(CD);
%             imgMatFile = [expt(iexp).redImg{1} expt(iexp).redImg_suffix{1} '.mat'];
%             load(imgMatFile);
%             nframes_red = info.config.frames;
%             fprintf(['Reading run ' expt(iexp).redImg{1} '- ' num2str(min(nframes_red)) ' frames \r\n'])
%             data = sbxread(imgMatFile(1,1:11),0,nframes_red);
%             if size(data,1) == 2
%                 red_data = squeeze(data(2,:,:,:));
%                 green_data = squeeze(data(1,:,:,:));
%                 [out, green_data_reg] = stackRegister(green_data,data_avg); %green_data = Behavior green @ 920 (1000 frames); data_avg = mean image (1 avg frames) fom middle of stack (500Frames)
%                 [out2, red_data_reg] = stackRegister_MA(red_data,[],[],out); % red_data = red @ 920 (1000 frames); out = green @ 920 shifts_xy
%                 red_data_avg = mean(red_data_reg,3); % Average of registered red image @ 920
%                 figure; imagesc(red_data_avg) % Show image
%                 title('Red at 920')
%                 green_data_avg = mean(green_data_reg,3);
%                 figure; imagesc(green_data_avg)
%                 title('Green at 920')
%                 if size(expt(iexp).redImg,2) == 2 % 1040 run exists
%                     CD = [data_base '\Data\2P_images\' expt(iexp).mouse '\' expt(iexp).date '\' expt(iexp).redImg{2}];
%                     cd(CD);
%                     imgMatFile = [expt(iexp).redImg{2} expt(iexp).redImg_suffix{2} '.mat'];
%                     load(imgMatFile);
%                     nframes_red = info.config.frames;
%                     fprintf(['Reading run ' expt(iexp).redImg{2} '- ' num2str(min(nframes_red)) ' frames \r\n'])
%                     data = sbxread(imgMatFile(1,1:11),0,nframes_red);
%                     red_data = squeeze(data(2,:,:,:)); % red data @ 1040 (1000 Frames)
%                     %mean(red_data)
%                     [out, red_data_reg] = stackRegister(red_data,red_data_avg); % red_data = red @ 1040 (1000 frames); red_data_avg = red @ 920 (1 avg frame)
%                     red_data_avg = mean(red_data_reg,3); % Average of registered red image @ 1040
%                     figure; imagesc(red_data_avg)
%                     title('Red at 1040')
%                 end
%                 save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_redData.mat']), 'green_data_avg', 'red_data_avg')
%             end
%     end    
%     %% Check red for stability
% 
%     data_avg = mean(data_reg(:,:,size(data_reg,3)-stack_spacing:end),3);
%     figure; 
%     subplot(2,2,1)
%     sz = size(data_avg);
%     rgb = zeros(sz(1),sz(2),3);
%     rgb(:,:,1) = red_data_avg./max(red_data_avg(:));
%     imagesc(red_data_avg);
%     colormap gray
%     subplot(2,2,2)
%     rgb(:,:,2) = data_avg./max(data_avg(:));
%     imagesc(rgb);
%     if size(expt(iexp).redImg,2) == 2
%         title('Red at 1040; Green at 920')
%     else
%         title('Red at 920; Green at 920')
%     end
%     print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')
% 
%     clearvars sz rgb
% end 


%% Check 2AFC photodiode (info must exist)
tic
irun = 1;
nrun = 1;
run_str = ['runs']; 
run_str = [run_str '-' expt(iexp).runs(irun,:)];

CD = [data_base '\Data\2P_images\' expt(iexp).mouse '\' expt(iexp).date '\' expt(iexp).runs(irun,:)];
cd(CD);
imgMatFile = [expt(iexp).runs(irun,:) expt(iexp).runs_suffix(irun,:) '.mat']; % DONE; Make variable to pull for oriAdapt_V1 that points to imgMatFile of restarted runs (ex: 001_000_001)
load(imgMatFile);
fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' expt(iexp).mouse '-' expt(iexp).date '-' expt(iexp).time_mat(irun,:) '.mat'];
load(fName);

if exist(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
    load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
    load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'adapt_input')
end

nframes = [adapt_input.counterValues{end}(end) info.config.frames];



if exist([data_base '\Data\2P_images\' expt(iexp).mouse '\' expt(iexp).date '\' expt(iexp).runs(irun,:) '\' [expt(iexp).runs(irun,:) expt(iexp).runs_suffix(irun,:) '.ephys']]);
    photoData = [];
    for irun = 1:nrun
        filename = [data_base '\Data\2P_images\' expt(iexp).mouse '\' expt(iexp).date '\' expt(iexp).runs(irun,:) '\' [expt(iexp).runs(irun,:) expt(iexp).runs_suffix(irun,:) '.ephys']];
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

    tDoFB = celleqel2mat_padded(adapt_input.tDoFeedbackMotion);
    tFramesStimOn = celleqel2mat_padded(adapt_input.cStimOff)-celleqel2mat_padded(adapt_input.cStimOn);
    ind_fast = tFramesStimOn<adapt_input.nFramesTooFast;
    FBfast = tDoFB & ind_fast;
    ntrials = size(adapt_input.tGratingContrast,2);
    cAdaptOn = nan(1,ntrials);
    cStimOn = nan(1,ntrials);
    n1 = 1; % First adaptor (distractor)
    n2 = 5; % Stimulus (presentation)
    cAdaptOn(1) = stimOnFrames(n1);
    cStimOn(1) = stimOnFrames(n2);
    for itrial = 2:ntrials % Is dropping the last trial (ntrials-1) the incorrect way to fix this indexing issue? 3/16/22 - CLM
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
    print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str 'photoLoc_diff.pdf']),'-dpdf', '-bestfit')
    save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_photoData.mat']), 'photoLoc', 'stimOnFrames', 'cStimOn','cAdaptOn')
elseif isfield(info, "frame")
    [stimOnFrames stimOffFrames] = photoFrameFinder_Sanworks(info.frame);

    tDoFB = celleqel2mat_padded(adapt_input.tDoFeedbackMotion);
    tFramesStimOn = celleqel2mat_padded(adapt_input.cStimOff)-celleqel2mat_padded(adapt_input.cStimOn);
    ind_fast = tFramesStimOn<adapt_input.nFramesTooFast;
    FBfast = tDoFB & ind_fast;
    cAdaptOn = nan(1,ntrials);
    cStimOn = nan(1,ntrials);
    n1 = 1; % First adaptor (distractor)
    n2 = 5; % Stimulus (presentation)
    cAdaptOn(1) = stimOnFrames(n1);
    cStimOn(1) = stimOnFrames(n2);
    for itrial = 2:ntrials % Is dropping the last trial (ntrials-1) the incorrect way to fix this indexing issue? 3/16/22 - CLM
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
    save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_photoData.mat']), 'stimOnFrames', 'cStimOn','cAdaptOn')
    

else
    error("No photodiode data!!!")
end

clearvars photoData
toc

%%

%% 2AFC analysis *Quick look at time courses

load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat'])) 

frameRateHz = double(adapt_input.frameRateHz);
cDecision = cStimOn + (celleqel2mat_padded(adapt_input.cDecision)-celleqel2mat_padded(adapt_input.cStimOn));
tGratingOri = celleqel2mat_padded(adapt_input.tGratingDirectionStart);
b2Ix = celleqel2mat_padded(adapt_input.tBlock2TrialNumber);
tOris = unique(tGratingOri);
nOri = length(tOris);
aGratingOri = celleqel2mat_padded(adapt_input.aGratingDirectionDeg);
aGratingContrast = celleqel2mat_padded(adapt_input.aGratingContrast);
aCons = unique(aGratingContrast);
naCon = length(aCons);
aOris = unique(aGratingOri);
naOri = length(aOris);
nCells = size(npSub_tc,2);
nframes = size(npSub_tc,1);
nTrials = size(aGratingOri,2);
data_stim = nan(50,nCells,nTrials);
data_stim_z = nan(50,nCells,nTrials);
data_adapt = nan(100,nCells,nTrials);
data_dec = nan(50,nCells,nTrials);
tc_z = npSub_tc./std(npSub_tc,[],1);
for itrial = 1:nTrials
    if cStimOn(itrial)+29 < nframes
        data_stim(:,:,itrial) = npSub_tc(cStimOn(itrial)-20:cStimOn(itrial)+29,:);
        data_stim_z(:,:,itrial) = tc_z(cStimOn(itrial)-20:cStimOn(itrial)+29,:);
    end
    if cAdaptOn(itrial)+79< nframes
        data_adapt(:,:,itrial) = npSub_tc(cAdaptOn(itrial)-20:cAdaptOn(itrial)+79,:);
    end
    if ~isnan(cDecision(itrial))
        if cDecision(itrial)+29 < nframes
            data_dec(:,:,itrial) = npSub_tc(cDecision(itrial)-20:cDecision(itrial)+29,:);
        end
    end
end
dataf = mean(data_adapt(1:20,:,:),1);
data_stim_dfof = bsxfun(@rdivide, bsxfun(@minus, data_stim, dataf), dataf);
data_adapt_dfof = bsxfun(@rdivide, bsxfun(@minus, data_adapt, dataf), dataf);
data_dec_dfof = bsxfun(@rdivide, bsxfun(@minus, data_dec, dataf), dataf);
tt = [-20:29].*(1000./frameRateHz);
tt_adapt = [-20:79].*(1000./frameRateHz);
figure;
subplot(1,2,1)
plot(nanmean(mean(data_adapt_dfof,2),3));
vline([20 31 42 53])
title('Adapt')
subplot(1,2,2)
plot(nanmean(mean(data_stim_dfof,2),3));
vline([20 31 42 53])
title('Target')

% base_win = [20:23]; % Unused
% resp_win = [27:31];

nt = cell(naOri,nOri);
x = 1;
start = 1;
figure;
for icon = 1:naCon
    if aCons(icon) == 1
        for iaOri = 1:naOri
            ind_con = find(aGratingContrast == aCons(icon));
            ind_con = intersect(ind_con, find(aGratingOri == aOris(iaOri)));
            for iori = 1:nOri
                ind_ori = intersect(ind_con, find(tGratingOri == tOris(iori)));
                subplot(naOri+1,nOri,start)
                plot(tt, mean(nanmean(data_stim_dfof(:,:,ind_ori),3),2))
                hold on
                title(['Adapt: Ori = '  num2str(aOris(iaOri))])
                ylim([-0.01 0.1])
                start = start+1;
                nt{x,iori} = ind_ori;
            end
            x = 1+x;
        end
    else
        ind_con = find(aGratingContrast == aCons(icon));
        for iori = 1:nOri
            ind_ori = intersect(ind_con, find(tGratingOri == tOris(iori)));
            subplot(naOri+1,nOri,start)
            plot(tt, mean(nanmean(data_stim_dfof(:,:,ind_ori),3),2))
            hold on
            title(['No adapt; Ori = ' num2str(tOris(iori))])
            ylim([-0.01 0.1])
            start = start+1;
            nt{x,iori} = ind_ori;
        end
        x = 1+x;
    end
end
print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_avgTargetResponse_byTarget_update.pdf']), '-dpdf','-bestfit')

SIx = strcmp(adapt_input.trialOutcomeCell,'success');
MIx = strcmp(adapt_input.trialOutcomeCell,'incorrect');
IIx = strcmp(adapt_input.trialOutcomeCell,'ignore');
tLeftTrial = celleqel2mat_padded(adapt_input.tLeftTrial);
tLeftResp = zeros(size(SIx));
tLeftResp(find(tLeftTrial&SIx)) = 1;
tLeftResp(find(~tLeftTrial&~SIx)) = 1;

x = 1;
start = 1;
figure;
for icon = 1:naCon
    if aCons(icon) == 1
        for iaOri = 1:naOri
            ind_con = find(aGratingContrast == aCons(icon) & SIx);
            ind_con = intersect(ind_con, find(aGratingOri == aOris(iaOri)));
            for iori = 1:nOri
                ind_ori = intersect(ind_con, find(tGratingOri == tOris(iori)));
                subplot(naOri+1,nOri,start)
                plot(tt, mean(nanmean(data_stim_dfof(:,:,ind_ori),3),2))
                hold on
                title(['Adapt: Ori = '  num2str(aOris(iaOri))])
                ylim([-0.01 0.1])
                start = start+1;
                nt{x,iori} = ind_ori;
            end
            x = 1+x;
        end
    else
        ind_con = find(aGratingContrast == aCons(icon) & SIx);
        for iori = 1:nOri
            ind_ori = intersect(ind_con, find(tGratingOri == tOris(iori)));
            subplot(naOri+1,nOri,start)
            plot(tt, mean(nanmean(data_stim_dfof(:,:,ind_ori),3),2))
            hold on
            title(['No adapt; Ori = ' num2str(tOris(iori))])
            ylim([-0.01 0.1])
            start = start+1;
            nt{x,iori} = ind_ori;
        end
        x = 1+x;
    end
end
x = 1;
start = 1;
for icon = 1:naCon
    if aCons(icon) == 1
        for iaOri = 1:naOri
            ind_con = find(aGratingContrast == aCons(icon) & MIx);
            ind_con = intersect(ind_con, find(aGratingOri == aOris(iaOri)));
            for iori = 1:nOri
                ind_ori = intersect(ind_con, find(tGratingOri == tOris(iori)));
                subplot(naOri+1,nOri,start)
                plot(tt, mean(nanmean(data_stim_dfof(:,:,ind_ori),3),2))
                hold on
                title(['Adapt: Ori = '  num2str(aOris(iaOri))])
                ylim([-0.01 0.1])
                start = start+1;
                nt{x,iori} = ind_ori;
            end
            x = 1+x;
        end
    else
        ind_con = find(aGratingContrast == aCons(icon) & MIx);
        for iori = 1:nOri
            ind_ori = intersect(ind_con, find(tGratingOri == tOris(iori)));
            subplot(naOri+1,nOri,start)
            plot(tt, mean(nanmean(data_stim_dfof(:,:,ind_ori),3),2))
            hold on
            title(['No adapt; Ori = ' num2str(tOris(iori))])
            ylim([-0.01 0.1])
            start = start+1;
            nt{x,iori} = ind_ori;
        end
        x = 1+x;
    end
end
legend({'Hit','Miss'})
print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_avgTargetResponse_byTarget_byOutcome_update.pdf']), '-dpdf','-bestfit')


ind_aCon0 = find(aGratingContrast == 0);
ind_aCon1 = find(aGratingContrast);
ind_cond{1} = intersect(find(SIx),ind_aCon0);
leg_str = {'Con = 0'};
for i = 1:naOri
    ind_cond{i+1} = intersect(find(SIx),intersect(ind_aCon1,find(aGratingOri==aOris(i))));
    leg_str{i+1} = ['Ori = ' num2str(aOris(i))];
end

aCon_p1 = [NaN aGratingContrast];
aOri_p1 = [NaN aGratingOri];
ind_aCon0_p1{1} = intersect(ind_aCon0,find(aCon_p1==0));
for i = 1:naOri
    ind_aCon0_p1{i+1} = intersect(ind_aCon0,intersect(find(aCon_p1==1),find(aOri_p1 ==aOris(i))));
end

figure;
subplot(1,2,1)
for i = 1:naOri+1
 plot(tt, mean(mean(data_stim_dfof(:,:,ind_cond{i}),3),2))
 hold on
end
title('Current trial')
ylabel('dF/F')
xlabel('Time from target (ms)')
legend(leg_str)
subplot(1,2,2)
for i = 1:naOri+1
 plot(tt, mean(mean(data_stim_dfof(:,:,ind_aCon0_p1{i}),3),2))
 hold on
end
title('Previous trial, for Current Con = 0')
ylabel('dF/F')
xlabel('Time from target (ms)')
print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_avgTargetResponse_updated.pdf']), '-dpdf','-bestfit')

interval = ceil(64/nOri)+1;
x = 1:interval:64;
x(end) = 63;
y = bluered;

if strcmp(cell2mat(expt(iexp).img_strct),'axons')
    mask_label = zeros(1,nCells);

elseif strcmp(cell2mat(expt(iexp).img_strct),'cells')
    mask_label = zeros(1,nCells);
    if ~isempty(expt(iexp).redImg)

        for i = 1:nCells
            if mask_cell_red(find(mask_cell == i, 1))
                mask_label(1,i) = 1;
            end
        end
    end

    figure;
    if nCells <49
        [n n2] = subplotn(nCells);
    else
        [n n2] = subplotn(49);
    end

    start = 1;
    for iC = 1:nCells
        if start > 49
             sgtitle('Adapt response')
            print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_adaptResp_allCells_success_cells_updated' num2str(iC-49) '-' num2str(iC-1) '.pdf']), '-dpdf','-bestfit')
            figure;
            start = 1;
        end
        subplot(n,n2,start)
        col_mat = [x(1) x(end)];
        for i = 1:naOri
            ind_adapt = find(SIx & aGratingContrast & aGratingOri==aOris(i));
            plot(tt_adapt, mean(data_adapt_dfof(:,iC,ind_adapt),3),'Color',y(col_mat(i),:))
            hold on
        end
        ylim([-0.05 0.3])
        if mask_label(iC)
            title([num2str(iC) '-' expt(iexp).driver])
        else
            title(num2str(iC))
        end
        start = start+1;
    end
     sgtitle('Adapt response')
    print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_adaptResp_allCells_success_cells_updated' num2str(iC-48) '-' num2str(iC) '.pdf']), '-dpdf','-bestfit')

    figure;
    start = 1;
    for iC = 1:nCells
        if start > 49
             sgtitle('Target response- All trials- All correct')
            print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_targetResp_allCells_success_cells_updated' num2str(iC-49) '-' num2str(iC-1) '.pdf']), '-dpdf','-bestfit')
            figure;
            start = 1;
        end
        subplot(n,n2,start)
        for iori = 1:nOri
            ind_ori = find(SIx & tGratingOri == tOris(iori));
            plot(tt, mean(data_stim_dfof(:,iC,ind_ori),3),'Color',y(x(iori),:))
            hold on
        end
        if mask_label(iC)
            title([num2str(iC) '-' expt(iexp).driver])
        else
            title(num2str(iC))
        end
        start = start+1;
    end
     sgtitle('Target response- All trials- All correct')
    print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_targetResp_allCells_success_cells_updated' num2str(iC-48) '-' num2str(iC) '.pdf']), '-dpdf','-bestfit')


    figure;
    start = 1;
    for iC = 1:nCells
        if start > 49
             sgtitle('Target response by choice- Trials <20 deg')
            print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_targetResp_allCells_bychoice_cells_updated' num2str(iC-49) '-' num2str(iC-1) '.pdf']), '-dpdf','-bestfit')
            figure;
            start = 1;
        end
        subplot(n,n2,start)
        plot(tt, mean(data_stim_dfof(:,iC,find(~tLeftResp&~IIx&abs(tGratingOri)<20)),3))
        hold on
        plot(tt, mean(data_stim_dfof(:,iC,find(tLeftResp&~IIx&abs(tGratingOri)<20)),3))
        if mask_label(iC)
            title([num2str(iC) '-' expt(iexp).driver])
        else
            title(num2str(iC))
        end
        start = start+1;
    end
     sgtitle('Target response by choice- Trials <20 deg')
    print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_targetResp_allCells_bychoice_cells_updated' num2str(iC-48) '-' num2str(iC) '.pdf']), '-dpdf','-bestfit')

    
    figure;
    start = 1;
    for iC = 1:nCells
        if start > 49
             sgtitle('Target response by target- Trials <20 deg')
            print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_targetResp_allCells_bytarget_cells_updated' num2str(iC-49) '-' num2str(iC-1) '.pdf']), '-dpdf','-bestfit')
            figure;
            start = 1;
        end
        subplot(n,n2,start)
        plot(tt, mean(data_stim_dfof(:,iC,find(~tLeftTrial&~IIx&abs(tGratingOri)<20)),3))
        hold on
        plot(tt, mean(data_stim_dfof(:,iC,find(tLeftTrial&~IIx&abs(tGratingOri)<20)),3))
        if mask_label(iC)
            title([num2str(iC) '-' expt(iexp).driver])
        else
            title(num2str(iC))
        end
        start = start+1;
    end
    sgtitle('Target response by target- Trials <20 deg')
    print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_targetResp_allCells_bytarget_cells_updated' num2str(iC-48) '-' num2str(iC) '.pdf']), '-dpdf','-bestfit')

    
    figure;
    start = 1;
    for iC = 1:nCells
        if start > 49
            sgtitle('Decision response- All trials')
            print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_decResp_allCells_success_cells_updated' num2str(iC-49) '-' num2str(iC-1) '.pdf']), '-dpdf','-bestfit')
            figure;
            start = 1;
        end
        subplot(n,n2,start)
        for iori = 1:nOri
            ind_ori = find(SIx & tGratingOri == tOris(iori));
            plot(tt, mean(data_dec_dfof(:,iC,ind_ori),3),'Color',y(x(iori),:))
            hold on
        end
        if mask_label(iC)
            title([num2str(iC) '-' expt(iexp).driver])
        else
            title(num2str(iC))
        end
        start = start+1;
    end
    sgtitle('Decision response- All trials')
    print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_decResp_allCells_success_cells_updated' num2str(iC-48) '-' num2str(iC) '.pdf']), '-dpdf','-bestfit')
end

figure;
subplot(3,2,1)
for i = 1:naOri
    ind_adapt = find(SIx & aGratingContrast & aGratingOri==aOris(i));
    plot(tt_adapt, mean(mean(data_adapt_dfof(:,find(mask_label),ind_adapt),3),2),'Color',y(col_mat(i),:))
    hold on
end
title([expt(iexp).driver '+ n = ' num2str(sum(mask_label))])
xlabel('Time from adapt on (ms)')
ylim([-0.02 0.15])
subplot(3,2,2)
for i = 1:naOri
    ind_adapt = find(SIx & aGratingContrast & aGratingOri==aOris(i));
    plot(tt_adapt, mean(mean(data_adapt_dfof(:,find(~mask_label),ind_adapt),3),2),'Color',y(col_mat(i),:))
    hold on
end
title([expt(iexp).driver '- n = ' num2str(sum(~mask_label))])
xlabel('Time from adapt on (ms)')
ylim([-0.02 0.15])
subplot(3,2,3) 
for iori = 1:nOri
    ind_ori = find(SIx & tGratingOri == tOris(iori));
    plot(tt, mean(mean(data_stim_dfof(:,find(mask_label),ind_ori),3),2),'Color',y(x(iori),:))
    hold on
end
xlabel('Time from target (ms)')
ylim([-0.02 0.15])
subplot(3,2,4)
for iori = 1:nOri
    ind_ori = find(SIx & tGratingOri == tOris(iori));
    plot(tt, mean(mean(data_stim_dfof(:,find(~mask_label),ind_ori),3),2),'Color',y(x(iori),:))
    hold on
end
xlabel('Time from target (ms)')
ylim([-0.02 0.15])
subplot(3,2,5)
for iori = 1:nOri
    ind_ori = find(SIx & tGratingOri == tOris(iori));
    plot(tt, mean(mean(data_dec_dfof(:,find(mask_label),ind_ori),3),2),'Color',y(x(iori),:))
    hold on
end
xlabel('Time from decision (ms)')
ylim([-0.02 0.15])
subplot(3,2,6)
for iori = 1:nOri
    ind_ori = find(SIx & tGratingOri == tOris(iori));
    plot(tt, mean(mean(data_dec_dfof(:,find(~mask_label),ind_ori),3),2),'Color',y(x(iori),:))
    hold on
end
xlabel('Time from decision (ms)')
ylim([-0.02 0.15])
sgtitle(['All corrects - n = ' num2str(sum(SIx))])
print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_avgResp_allCells_success_updated.pdf']), '-dpdf','-bestfit')


figure;
subplot(3,2,1)
for i = 1:naOri
    ind_adapt = find(MIx & aGratingContrast & aGratingOri==aOris(i));
    plot(tt_adapt, mean(mean(data_adapt_dfof(:,find(mask_label),ind_adapt),3),2),'Color',y(col_mat(i),:))
    hold on
end
title([expt(iexp).driver '+ n = ' num2str(sum(mask_label))])
xlabel('Time from adapt on (ms)')
ylim([-0.02 0.15])
subplot(3,2,2)
for i = 1:naOri
    ind_adapt = find(MIx & aGratingContrast & aGratingOri==aOris(i));
    plot(tt_adapt, mean(mean(data_adapt_dfof(:,find(~mask_label),ind_adapt),3),2),'Color',y(col_mat(i),:))
    hold on
end
title([expt(iexp).driver '- n = ' num2str(sum(~mask_label))])
xlabel('Time from adapt on (ms)')
ylim([-0.02 0.15])
subplot(3,2,3) 
for iori = 1:nOri
    ind_ori = find(MIx & tGratingOri == tOris(iori));
    plot(tt, mean(mean(data_stim_dfof(:,find(mask_label),ind_ori),3),2),'Color',y(x(iori),:))
    hold on
end
xlabel('Time from target (ms)')
ylim([-0.02 0.15])
subplot(3,2,4)
for iori = 1:nOri
    ind_ori = find(MIx & tGratingOri == tOris(iori));
    plot(tt, mean(mean(data_stim_dfof(:,find(~mask_label),ind_ori),3),2),'Color',y(x(iori),:))
    hold on
end
xlabel('Time from target (ms)')
ylim([-0.02 0.15])
subplot(3,2,5)
for iori = 1:nOri
    ind_ori = find(MIx & tGratingOri == tOris(iori));
    plot(tt, mean(mean(data_dec_dfof(:,find(mask_label),ind_ori),3),2),'Color',y(x(iori),:))
    hold on
end
xlabel('Time from decision (ms)')
ylim([-0.02 0.15])
subplot(3,2,6)
for iori = 1:nOri
    ind_ori = find(MIx & tGratingOri == tOris(iori));
    plot(tt, mean(mean(data_dec_dfof(:,find(~mask_label),ind_ori),3),2),'Color',y(x(iori),:))
    hold on
end
xlabel('Time from decision (ms)')
ylim([-0.02 0.15])
sgtitle(['All incorrects - n = ' num2str(sum(MIx))])
print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_avgResp_allCells_incorrect_updated.pdf']), '-dpdf','-bestfit')

%% adapt analysis - 2AFC * Perform separately from TC extraction
close all
base_win = 19:21;
resp_win = 25:27;
% base_win = 15:21;
% resp_win = 25:31;
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

figure;
start = 1;
for IN = 1:2
    subplot(naOri+1,2,start)
    cell_ind = intersect(stim_resp_ind, find(mask_label==IN-1));
    plot(nanmean(nanmean(data_stim_dfof(:,cell_ind,ind_con),3),2))
    if IN == 1
        IN_str = [expt(iexp).driver '-'];
    else
        IN_str = [expt(iexp).driver '+'];
    end
    vline(base_win, 'b')
    vline(resp_win, 'r')
    title(['Target resp; ' num2str(length(cell_ind)) ' ' IN_str ' cells'])
    start = start+1;
end

% if ~isfield(adapt_input,'cAdaptOn_all')
    nFramesPerAdapt = double(ceil(((double(adapt_input.dynAdaptPeriodMs + adapt_input.dynAdaptFlashOffMs))./1000).*frameRateHz));
    base_win_all = zeros(length(base_win),4,2);
    resp_win_all = zeros(length(resp_win),4,2);
    for i = 1:4
        base_win_all(:,i,1) = base_win+((i-1).*nFramesPerAdapt);
        resp_win_all(:,i,1) = resp_win+((i-1).*nFramesPerAdapt);
        base_win_all(:,i,2) = base_win+((i-1).*nFramesPerAdapt);
        resp_win_all(:,i,2) = resp_win+((i-1).*nFramesPerAdapt);
    end
% else
%     error('need code for cAdaptOn_all')
% end

for iOri = 1:naOri
    for IN = 1:2
        subplot(naOri+1,2,start)
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
print(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_avgResp_Adapt&Target_updated.pdf']), '-dpdf','-bestfit')


adapt_cyc_resp = zeros(nCells,nTrials,4);
for IN = 1:2
    for i = 1:4
        adapt_cyc_resp(find(mask_label==IN-1),:,i) = squeeze(nanmean(data_adapt_dfof(resp_win_all(:,i,IN),find(mask_label==IN-1),:),1)...
            -nanmean(data_adapt_dfof(base_win_all(:,i,IN),find(mask_label==IN-1),:),1));
    end
end

save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_adaptResp.mat']), 'mask_label', 'data_adapt_dfof', 'adapt_cyc_resp', 'base_win_all', 'resp_win_all','adapt_resp_ind')
save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimResp.mat']), 'stim_resp', 'stim_resp_ind', 'base_win', 'resp_win');
save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']), 'tGratingOri', 'tOris', 'aGratingOri', 'aOris', 'aGratingContrast', 'ind_cond', 'SIx', 'MIx');

end