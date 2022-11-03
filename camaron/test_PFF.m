% all_days = [i475_good_list i472_good_list];
dataset = 'oriAdapt_V1_cam';
eval(dataset); % run file to load expt.structure

% days_to_fix = [];
fixed_glitch = [];
u_photoLoc_diff_cell_move_base = {};

load('Z:\All_Staff\home\camaron\Analysis\2P\behavior_pFF_days_to_fix.mat')

all_days = days_to_fix;

for iday = 1:length(all_days)

clearvars -except all_days iday expt days_to_fix u_photoLoc_diff_cell u_photoLoc_diff_cell_move_base fixed_glitch
clc
close all
%%

disp([num2str(iday) ' of ' num2str(length(all_days))])


iexp = all_days(iday); % Enter experiment number from oriAdapt_V1

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
CM_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron';

if strcmp(expt(iexp).folder,'lindsey')
    data_base = LG_base;
elseif strcmp(expt(iexp).folder,'camaron')
    data_base = CM_base;
end

mouse = expt(iexp).mouse;
date = expt(iexp).date;
%%

% load data
nrun = size(expt(iexp).runs,1);
for irun = 1:nrun
    CD = [data_base '\Data\2P_images\' expt(iexp).mouse '\' expt(iexp).date '\' expt(iexp).runs(irun,:)];
    cd(CD);
    imgMatFile = [expt(iexp).runs(irun,:) expt(iexp).runs_suffix(irun,:) '.mat']; % DONE; Make variable to pull for oriAdapt_V1 that points to imgMatFile of restarted runs (ex: 001_000_001)
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' expt(iexp).mouse '-' expt(iexp).date '-' expt(iexp).time_mat(irun,:) '.mat'];
    load(fName);
    
    nframes = [input.counterValues{end}(end) info.config.frames];
end

%%
irun = 1;
nrun = 1;
run_str = ['runs']; 
run_str = [run_str '-' expt(iexp).runs(irun,:)];

if exist(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
    load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
    load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'adapt_input')
end

nframes = [adapt_input.counterValues{end}(end) info.config.frames];



% with photodiode
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

% Original
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
    title("original")
    u_photoLoc_diff = unique(photoLoc_diff);

if ~all(u_photoLoc_diff < 4333)
    fprintf("Photodiode glitch found!")
%     days_to_fix = [days_to_fix iexp]
    fixed_glitch = [fixed_glitch iexp]
    u_photoLoc_diff_cell_move_base{iday} = u_photoLoc_diff;
end

u_photoLoc_diff_cell_move_base =  u_photoLoc_diff_cell_move_base(~cellfun('isempty',u_photoLoc_diff_cell_move_base));

% Move Baseline
%     [photoLoc stimOnFrames] = photoFrameFinder_movBase(photoData,min(nframes));
% 
%     frameDiff = diff(stimOnFrames);
%     ind_long = find(frameDiff>20);
%     ind_long_long = ind_long(find(frameDiff(ind_long-1)>20)); % ??
%     photoLoc(ind_long_long) = [];
%     stimOnFrames(ind_long_long) = [];
% 
%     nf = rem(size(stimOnFrames,2),5); % rem = remainder after division; Why divide by 5? Adaptors plus target?
%     photoLoc_rs = reshape(photoLoc(1:end-nf),[5 length(photoLoc(1:end-nf))./5])'; 
%     photoLoc_diff = diff(photoLoc_rs,1,2);
%     figure; plot(photoLoc_diff'); ylim([0 6000])
%     title("move base")
%     u_photoLoc_diff_moveBase = unique(photoLoc_diff);

    %%

%     tDoFB = celleqel2mat_padded(adapt_input.tDoFeedbackMotion);
%     tFramesStimOn = celleqel2mat_padded(adapt_input.cStimOff)-celleqel2mat_padded(adapt_input.cStimOn);
%     ind_fast = tFramesStimOn<adapt_input.nFramesTooFast;
%     FBfast = tDoFB & ind_fast;
%     ntrials = size(adapt_input.tGratingContrast,2);
%     cAdaptOn = nan(1,ntrials);
%     cStimOn = nan(1,ntrials);
%     n1 = 1; % First adaptor (distractor)
%     n2 = 5; % Stimulus (presentation)
%     cAdaptOn(1) = stimOnFrames(n1);
%     cStimOn(1) = stimOnFrames(n2);
%     for itrial = 2:ntrials % Is dropping the last trial (ntrials-1) the incorrect way to fix this indexing issue? 3/16/22 - CLM
%         if FBfast(itrial-1)
%             n1 = n1+6;
%             n2 = n2+6;
%         else
%             n1 = n1+5;
%             n2 = n2+5;
%         end
%         cAdaptOn(itrial) = stimOnFrames(n1);
%         cStimOn(itrial) = stimOnFrames(n2);
%     end
% 
%     unique(cStimOn-cAdaptOn)
end

% pFF_cAdaptOn = cAdaptOn;

end

    %%

% 
% frameRateHz = double(adapt_input.frameRateHz);
% cDecision = cStimOn + (celleqel2mat_padded(adapt_input.cDecision)-celleqel2mat_padded(adapt_input.cStimOn));
% tGratingOri = celleqel2mat_padded(adapt_input.tGratingDirectionStart);
% b2Ix = celleqel2mat_padded(adapt_input.tBlock2TrialNumber);
% tOris = unique(tGratingOri);
% nOri = length(tOris);
% aGratingOri = celleqel2mat_padded(adapt_input.aGratingDirectionDeg);
% aGratingContrast = celleqel2mat_padded(adapt_input.aGratingContrast);
% aCons = unique(aGratingContrast);
% naCon = length(aCons);
% aOris = unique(aGratingOri);
% naOri = length(aOris);
% nCells = size(npSub_tc,2);
% nframes = size(npSub_tc,1);
% nTrials = size(aGratingOri,2);
% data_stim = nan(50,nCells,nTrials);
% data_stim_z = nan(50,nCells,nTrials);
% data_adapt = nan(100,nCells,nTrials);
% data_dec = nan(50,nCells,nTrials);
% tc_z = npSub_tc./std(npSub_tc,[],1);
% for itrial = 1:nTrials
%     if cStimOn(itrial)+29 < nframes
%         data_stim(:,:,itrial) = npSub_tc(cStimOn(itrial)-20:cStimOn(itrial)+29,:);
%         data_stim_z(:,:,itrial) = tc_z(cStimOn(itrial)-20:cStimOn(itrial)+29,:);
%     end
%     if cAdaptOn(itrial)+79< nframes
%         data_adapt(:,:,itrial) = npSub_tc(cAdaptOn(itrial)-20:cAdaptOn(itrial)+79,:);
%     end
%     if ~isnan(cDecision(itrial))
%         if cDecision(itrial)+29 < nframes
%             data_dec(:,:,itrial) = npSub_tc(cDecision(itrial)-20:cDecision(itrial)+29,:);
%         end
%     end
% end
% dataf = mean(data_adapt(1:20,:,:),1);
% data_stim_dfof = bsxfun(@rdivide, bsxfun(@minus, data_stim, dataf), dataf);
% data_adapt_dfof = bsxfun(@rdivide, bsxfun(@minus, data_adapt, dataf), dataf);
% data_dec_dfof = bsxfun(@rdivide, bsxfun(@minus, data_dec, dataf), dataf);
% figure;
% plot(nanmean(mean(data_adapt_dfof,2),3));
% vline([20 31 42 53])
% title('Adapt')
% 
% 
% 
% % Without photodiode
% cStimOn = celleqel2mat_padded(adapt_input.cStimOn);
% cAdaptOn = celleqel2mat_padded(adapt_input.cAdaptOn);
% data_stim = nan(50,nCells,nTrials);
% data_stim_z = nan(50,nCells,nTrials);
% data_adapt = nan(100,nCells,nTrials);
% data_dec = nan(50,nCells,nTrials);
% tc_z = npSub_tc./std(npSub_tc,[],1);
% for itrial = 1:nTrials
%     if cStimOn(itrial)+29 < nframes
%         data_stim(:,:,itrial) = npSub_tc(cStimOn(itrial)-20:cStimOn(itrial)+29,:);
%         data_stim_z(:,:,itrial) = tc_z(cStimOn(itrial)-20:cStimOn(itrial)+29,:);
%     end
%     if cAdaptOn(itrial)+79< nframes
%         data_adapt(:,:,itrial) = npSub_tc(cAdaptOn(itrial)-20:cAdaptOn(itrial)+79,:);
%     end
%     if ~isnan(cDecision(itrial))
%         if cDecision(itrial)+29 < nframes
%             data_dec(:,:,itrial) = npSub_tc(cDecision(itrial)-20:cDecision(itrial)+29,:);
%         end
%     end
% end
% dataf = mean(data_adapt(1:20,:,:),1);
% data_stim_dfof = bsxfun(@rdivide, bsxfun(@minus, data_stim, dataf), dataf);
% data_adapt_dfof = bsxfun(@rdivide, bsxfun(@minus, data_adapt, dataf), dataf);
% data_dec_dfof = bsxfun(@rdivide, bsxfun(@minus, data_dec, dataf), dataf);
% hold on
% plot(nanmean(mean(data_adapt_dfof,2),3), 'r');
% 
% input_cAdaptOn = cAdaptOn;
% 
% diff_cAdaptOn = input_cAdaptOn - pFF_cAdaptOn;
% unique(diff_cAdaptOn)

