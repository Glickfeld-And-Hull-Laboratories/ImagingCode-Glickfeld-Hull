function [ret_npSub_tc_matched, ret_distance,resp_by_stim_matched,ret_dfof_trial_matched] = retinitopy_for_matched_data(nd, allDays, expt, mouse, fov_avg, masks, fitGeoTAf, instructions, inputStructure, validation_choice)

ret_npSub_tc_matched = cell(1,nd);
ret_distance = cell(1,nd);
resp_by_stim_matched=cell(1,nd);
ret_dfof_trial_matched=cell(1,nd);

for id = 1:nd
    sess = allDays(id);
    date = expt(sess).date;
    subnum = expt(sess).mouse;
    ImgFolder = expt(sess).ret_run;
    imgMatFile = [ImgFolder '_000_000.mat'];
    time = expt(sess).ret_time;

    retDataPath = fullfile('Z:\home\ACh\Data\2p_data',mouse,date,ImgFolder);
    load(fullfile(retDataPath,imgMatFile));
    nframes = info.config.frames;
    ret_data_temp = squeeze(sbxread(fullfile(retDataPath,[ImgFolder '_000_000']),0,nframes));
    fprintf(['Loaded ' num2str(nframes) ' frames \r\n'])

    fName = ['Z:\Behavior\Data\data-' subnum '-' date '-' time '.mat'];
    load(fName);
    ret_inputStructure = input;
    clear input

    referenceFOV = fov_avg{id};
    retAvrg = mean(ret_data_temp,3);
    [~, retAvrg_registered] = stackRegGPU(retAvrg,referenceFOV);
    [~, ret_data_registered] = stackRegGPU(ret_data_temp,retAvrg_registered);

    if id == 2
        ret_data_registered = imwarp(ret_data_registered,fitGeoTAf, 'OutputView', imref2d(size(ret_data_registered)));
    end
    
    referenceMasks = masks{id};
    ret_data_tc = stackGetTimeCourses(ret_data_registered, referenceMasks);

    mask_np = imCellNeuropil(referenceMasks, 3, 5);
    
    nCells = size(ret_data_tc,2);
    data_tc_down = stackGetTimeCourses(stackGroupProject(ret_data_registered,5), referenceMasks);
    sz = size(ret_data_registered);
    down = 5;
    data_reg_down = stackGroupProject(ret_data_registered,down);
    np_tc = zeros(sz(3),nCells);
    np_tc_down = zeros(floor(sz(3)./down), nCells);
    for i = 1:nCells
        np_tc(:,i) = stackGetTimeCourses(ret_data_registered,mask_np(:,:,i));
        np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
        fprintf(['     Cell #' num2str(i) '%s/n'])
    end
    
    ii = 0.01:0.01:1;
    x = zeros(length(ii), nCells);
    for i = 1:100
        x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
    end
    [max_skew, ind] = max(x,[],1);
    np_w = 0.01*ind;

    ret_npSub_tc = ret_data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);

    nOn = ret_inputStructure(1).nScansOn;
    nOff = ret_inputStructure(1).nScansOff;
    
    switch instructions.tIdxSource
        case 'PD'
            [ret_stimOns_temp, ret_stimOffs] = photoFrameFinder_Sanworks(info.frame);
            ret_inputStructure.stimOns_photodiode = ret_stimOns_temp;
            ret_inputStructure.stimOffs_photodiode = ret_stimOffs;
            ret_inputStructure.stimOns_mwCounter = [];
            ret_inputStructure.stimTimingSource = 'PD';
        case 'MW'
            input_correct = counterValCorrect_noPhotodiode(input);
            ret_inputStructure.stimOns_mwCounter = cell2mat(input_correct.cStimOn);
            ret_inputStructure.counterValues = input_correct.counterValues;
            ret_inputStructure.counterTimesUs = input_correct.counterTimesUs;
            ret_inputStructure.stimOns_photodiode = [];
            ret_inputStructure.stimOffs_photodiode = [];
            ret_inputStructure.stimTimingSource = 'MW';
            ret_stimOns_temp = ret_inputStructure.stimOns_mwCounter;
            clear input_correct
        case 'cS'
            ret_inputStructure.stimTimingSource = 'cS';
            ret_stimOns_temp = cell2mat(ret_inputStructure.cStimOn);
        otherwise
            error('No valid trial indexing source specificed in instr file.');
    end
    
    nTrials = length(ret_stimOns_temp);
    ret_data_trial = nan(nOn+nOff,nCells,nTrials);
    for itrial = 1:nTrials
        if ~isnan(ret_stimOns_temp(itrial)) & (ret_stimOns_temp(itrial)+nOn+nOff/2)<size(ret_npSub_tc,1)
            ret_data_trial(:,:,itrial) = ret_npSub_tc((ret_stimOns_temp(itrial)-nOff/2:ret_stimOns_temp(itrial)-1+nOn+nOff/2),:);
        end
    end
    
    baselineFrames = (nOff/4+1):nOff/2;
    stimFrames = (nOff/2+1):(nOff/2+nOn);
    
    F0 = mean(ret_data_trial(baselineFrames,:,:), 1);
    ret_dfof_trial = (ret_data_trial - F0) ./ F0;
    
    ret_mean_resp = squeeze(mean(ret_dfof_trial(stimFrames,:,:), 1));
    
    trialAz = celleqel2mat_padded(ret_inputStructure.tGratingAzimuthDeg);
    azs = unique(trialAz);
    trialEl = celleqel2mat_padded(ret_inputStructure.tGratingElevationDeg);
    els = unique(trialEl);
    
    resp_by_stim = nan(length(els),length(azs),nCells);
    for i_el = 1:length(els)
        this_el_trials = find(trialEl == els(i_el));
        for i_az = 1:length(azs)
            this_az_trials = find(trialAz==azs(i_az));
            these_trials = intersect(this_el_trials, this_az_trials);
            resp_by_stim(i_el, i_az,:) = nanmean(ret_mean_resp(:,these_trials),2);
        end
    end

    [nElev, nAzim, nCells] = size(resp_by_stim);
    resp_reshaped = reshape(resp_by_stim, [], nCells);
    [~, maxIdx] = max(resp_reshaped, [], 1);
    [maxElev, maxAzim] = ind2sub([nElev, nAzim], maxIdx);
    
    if validation_choice && nCells >= 5
        rng('shuffle');
        selected_cells = randperm(nCells, min(5, nCells));
        
        for i_cell = 1:length(selected_cells)
            this_cell = selected_cells(i_cell);
            
            figure('Name', sprintf('Day %d - Cell %d', id, this_cell));
            for i_el = 1:length(els)
                for i_az = 1:length(azs)
                    subplot(length(els), length(azs), (i_el-1)*length(azs) + i_az);
                    
                    this_el_trials = find(trialEl == els(i_el));
                    this_az_trials = find(trialAz == azs(i_az));
                    these_trials = intersect(this_el_trials, this_az_trials);
                    
                    if ~isempty(these_trials)
                        mean_trace = squeeze(nanmean(ret_dfof_trial(:, this_cell, these_trials), 3));
                        plot(mean_trace, 'k', 'LineWidth', 1);
                        hold on;
                        xline(nOff/2, 'r--');
                        xline(nOff/2 + nOn, 'r--');
                    end
                    
                    title(sprintf('Az:%d El:%d', azs(i_az), els(i_el)));
                    set(gca, 'TickDir', 'out');
                    grid off;
                    box off;
                    
                    if i_el == length(els)
                        xlabel('Frame');
                    end
                    if i_az == 1
                        ylabel('dF/F');
                    end
                end
            end
        end
    end
    
    finalAzim = double(inputStructure(id).gratingAzimuthDeg);
    finalElev = double(inputStructure(id).gratingElevationDeg);
    
    maxAzimDeg = azs(maxAzim);
    maxElevDeg = els(maxElev);
    
    ret_distance{id} = sqrt((maxAzimDeg - finalAzim).^2 + (maxElevDeg - finalElev).^2);
    ret_npSub_tc_matched{id} = ret_npSub_tc;
    resp_by_stim_matched{id}=resp_by_stim;
    ret_dfof_trial_matched{id}=ret_dfof_trial;
end

end