%%
close all
clearvars
clc
%% Load experiment info
dataset = 'oriAdapt_V1_cam';
eval(dataset); % run file to load expt.structure

load('Z:\All_Staff\home\camaron\Analysis\2P\good_expt_list.mat')


list = sort([i472_good_list, i475_good_list]);

for iexp = list
    clearvars -except expt list iexp
    
    LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
    CM_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron';
    
    if strcmp(expt(iexp).folder,'lindsey')
        data_base = LG_base;
    elseif strcmp(expt(iexp).folder,'camaron')
        data_base = CM_base;
    end
    
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    
    %% Load neccesary variables 
    irun = 1;
    nrun = 1;
    run_str = ['runs']; 
    run_str = [run_str '-' expt(iexp).runs(irun,:)];
    load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'adapt_input')
    load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
    load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_photoData.mat']), 'stimOnFrames', 'cStimOn','cAdaptOn')
    
    
    
    %% 2AFC analysis *Quick look at time courses
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
    
    

    save(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimResp.mat']), 'data_stim_dfof','data_dec_dfof', '-append');


end

%%
% load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimResp.mat']))
