
%% retinotopy analysis, work in progress
mask_np_matched=cell(1,nd);
ret_npSub_tc_matched=cell(1,nd);
ret_stimOns =cell(1,nd);
for id = 1:nd
    sess = allDays(id);
    date = expt(sess).date;
    subnum = expt(sess).mouse;
    ImgFolder =  expt(sess).ret_run;
    imgMatFile = [ImgFolder '_000_000.mat'];
    time = expt(sess).ret_time;

    %load neural data for retinotopy
    retDataPath = fullfile('Z:\home\ACh\Data\2p_data',mouse,date,ImgFolder);
    load(fullfile(retDataPath,imgMatFile));
    nframes = info.config.frames;
    ret_data_temp = squeeze(sbxread(fullfile(retDataPath,[ImgFolder '_000_000']),0,nframes));
    fprintf(['Loaded ' num2str(nframes) ' frames \r\n'])

    %load mWorks stimulus data for retinotopy
    fName = ['Z:\Behavior\Data\data-' subnum '-' date '-' time '.mat'];
    load(fName);
    ret_inputStructure = input;
    clear input


   %register the retinotopy data to the stimulus run avrg FOV, using the un-rotated FOV of the matched day 
   referenceFOV = fov_avg{id};
   retAvrg = mean(ret_data_temp,3);
   [~, retAvrg_registered]=stackRegGPU(retAvrg,referenceFOV);
   [~, ret_data_registered]=stackRegGPU(ret_data_temp,retAvrg_registered);

   if id ==2
       apply the rotation
        ret_data_registered = imwarp(ret_data_registered,fitGeoTAf, 'OutputView', imref2d(size(ret_data_registered)));
   end
    % apply the masks
    referenceMasks = masks{id};
    mask_np = imCellNeuropil(referenceMasks, 3, 5);
    mask_np_matched{id}=mask_np;
    ret_data_tc = stackGetTimeCourses(ret_data_registered, referenceMasks);
    %re-create np masks and get np-subtracted timecourses for the ret run

    nCells = size(ret_data_tc,2);
    data_tc_down = stackGetTimeCourses(stackGroupProject(ret_data_registered,5), referenceMasks);
    clear np_tc np_tc_down
    sz = size(ret_data_registered);
    down = 5;
    data_reg_down  = stackGroupProject(ret_data_registered,down);
    np_tc = zeros(sz(3),nCells);
    np_tc_down = zeros(floor(sz(3)./down), nCells);
    for i = 1:nCells
        np_tc(:,i) = stackGetTimeCourses(ret_data_registered,mask_np(:,:,i));
        np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
        fprintf(['     Cell #' num2str(i) '%s/n']) 
    end
    %get weights by maximizing skew
    ii= 0.01:0.01:1;
    x = zeros(length(ii), nCells);
    for i = 1:100
        x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
    end
    [max_skew, ind] =  max(x,[],1);
    np_w = 0.01*ind;
    ret_npSub_tc = ret_data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);

% break into trials

    nOn = ret_inputStructure(1).nScansOn;
    nOff = ret_inputStructure(1).nScansOff;
    
     %this is the set of trial start times that will be used throughout the rest of this scr
    switch instructions.tIdxSource
        case 'PD'
            [ret_stimOns_temp, ret_stimOffs] = photoFrameFinder_Sanworks(info.frame);
            input.stimOns_photodiode = ret_stimOns_temp;
            input.stimOffs_photodiode = ret_stimOffs;
            input.stimOns_mwCounter = [];
            input.stimTimingSource = 'PD';
        case 'MW'
            input_correct = counterValCorrect_noPhotodiode(input);
            input.stimOns_mwCounter = cell2mat(input_correct.cStimOn);
            input.counterValues = input_correct.counterValues;
            input.counterTimesUs = input_correct.counterTimesUs;
            input.stimOns_photodiode = [];
            input.stimOffs_photodiode = [];
            input.stimTimingSource = 'MW';
            ret_stimOns_temp = input.stimOns_mwCounter;
            clear input_correct
        case 'cS'
            input.stimTimingSource = 'cS';
            ret_stimOns_temp = cell2mat(input.cStimOn);
        otherwise
            error('No valid trial indexing source specificed in instr file.');
    end
    nTrials = length(ret_stimOns_temp);
    ret_data_trial = nan(nOn+nOff,nCells,nTrials);
    for itrial = 1:nTrials
      if ~isnan(stimOns(itrial)) & (stimOns(itrial)+nOn+nOff/2)<nFrames
        ret_data_trial(:,:,itrial) = ret_npSub_tc(stimOns(itrial)-nOff/2:stimOns(itrial)-1+nOn+nOff/2,:);
      end
    end
    
    % Get unique stimulus parameters
    azs = unique(ret_inputStructure.tGratingAzimuthDeg);
    els = unique(ret_inputStructure.tGratingElevationDeg);
    
   % make a heatmap for each cell?
   %find the most responsive location
   %find the euclidian distance from the max to the position that was
   %used


    % the matched day will need to be handled differently becuase it first
    % needs to be rotated.

    %for the reference day
    %first will register the ret data mean to the stimulus run FOV, then
    %register the ret data to that shifted mean
    ret_npSub_tc_matched{id}=ret_npSub_tc;
    inputStructure(id).ret_stimOns=ret_stimOns_temp;
    %SAVE THESE VARIABLES


end

    