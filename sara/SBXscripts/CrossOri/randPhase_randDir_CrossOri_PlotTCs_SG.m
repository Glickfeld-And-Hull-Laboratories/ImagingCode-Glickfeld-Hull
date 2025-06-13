%% Plotting timecourses

clc; clear all; close all;
doRedCShannel = 0;
ds = 'CrossOriRandDirFourPhase_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = size(expt,2);
max_dist = 10;
frame_rate = 15;

%%

for iexp = [64]; 
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc{1};
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
    nrun = length(ImgFolder);
    run_str = 'runs-002';

    base = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara'];

    fprintf([mouse ' ' date '\n'])

% Load timecourse data, stim data
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_respData.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']))

    clear data_tc data_dfof_con_ph_tc_avg base_cell np_tc npSub_tc Val

    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_stimData.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_PatternCorrelationIndexFits.mat']))
    
    resp_ind = intersect(intersect(resp_ind_dir,p_dir),find(DSI>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end

    fprintf(['\nMouse: ' mouse '\nDate: ' date '\n'])
    fprintf(['nCells = ' num2str(nCells) '\n'])
    fprintf(['nCells resp = ' num2str(length(resp_ind)) '\n'])


    nTrials_pha = zeros(nStimDir,nMaskPhas);
    for id = 1:nStimDir
        for im = 1:nMaskPhas
           [memb ind] = ismember(trialInd{id,im,2},find(centroid_dist<max_dist));
           nTrials_pha(id,im) = sum(ind~=0);
        end
         fprintf(['Trials per mask phase = ' num2str(nTrials_pha(id,:)) '\n']) 
    end


    cStimOn = celleqel2mat_padded(input.cStimOneOn);
    cStimOff = celleqel2mat_padded(input.cStimOneOff);
    nTrials = length(cStimOn);
    if ~exist('sz','var')
        nCells = size(data_tc,2);
        nFrames = size(data_tc,1);
    else
    end

    framesbytrial = nan(nTrials(1), cStimOn(1)+cStimOff(1));

    for i = 1:nTrials
        framesbytrial(i,:) = (cStimOn(i)-cStimOn(1))+1:(cStimOn+cStimOff(i));
    end

    ind_mask = [];
    ind_stim = [];
    ind_phas = [];
    ind = [];
    trialInd = [];

    ind_p = cell(1,nMaskPhas);

    for ip = 1:nMaskPhas
        ind_p{ip} = intersect(find(maskPhas_all == maskPhas(ip)),find(centroid_dist<max_dist));
    end

    for iDir = 1:nStimDir
        ind_stimdir = intersect(find(stimDir_all == stimDirs(iDir)),find(centroid_dist<max_dist));
        ind_maskdir = intersect(find(maskDir_all == stimDirs(iDir)),find(centroid_dist<max_dist));
        ind_diralone = [intersect(ind_stimdir, ind_stimAlone), intersect(ind_maskdir, ind_maskAlone)];
        ind_dirplaid = [intersect(ind_stimdir, ind_plaid)];
        trialsperstim(iDir,1,1) = length(ind_diralone);
        trialInd{iDir,1,1} = ind_diralone;
        for ip = 1:nMaskPhas
            ind_dpplaid = intersect(ind_dirplaid,ind_p{ip});
            trialsperstim(iDir,ip,2) = length(ind_dpplaid);
            trialInd{iDir,ip,2} = ind_dpplaid;
        end
    end



% Plot individual and average time courses for all phases for cells

    % vsize = numel(respplaid_ind);
    % idx = randperm(vsize);
    % rand_ind = respplaid_ind(idx(1:15));
    
%     figure;
%     s=1;
    for ic = [46] %cell number %[rand_ind']
        figure;
        s=1;
        for id = 1:nStimDir
            for ip = 1:nMaskPhas
                if ip==1
                    ptrials = trialInd{id,1,1}; %set which stim you want to look at
                    ptrials = setdiff(ptrials, find(centroid_dist>max_dist));
                        X = 1:30; %taken from how data_dfof_tc is calculated in crossori_expt script (should equal off + on frames)
                        dfof_trialavg = mean(data_dfof_tc(:,ic,ptrials),3);
                        numtrials = size(squeeze(data_dfof_tc(:,ic,ptrials)),2);
                        dfof_trialstd = std(squeeze(data_dfof_tc(:,ic,ptrials)),0,2);
                        dfof_trialsem = dfof_trialstd ./ sqrt(numtrials);
                   subplot(12,5,s)
                       plot(X,dfof_trialavg, 'b','LineWidth',1);
                       shadedErrorBar(X,dfof_trialavg,dfof_trialsem);
                       xline(15,'--b')
                       ylim([-0.1 0.4])
                       set(gca,'TickDir','out'); box off; axis square
                       movegui('center')
                       if id==1
                        subtitle('grating')
                       end
                       if s>1
                           axis off
                       end
                    s=s+1;
                end
                ptrials = trialInd{id,ip,2}; %set which stim you want to look at
                ptrials = setdiff(ptrials, find(centroid_dist>max_dist));
                    X = [1:30]; %taken from how data_dfof_tc is calculated in crossori_expt script (should equal off + on frames)
                    dfof_trialavg = mean(data_dfof_tc(:,ic,ptrials),3);
                    numtrials = size(squeeze(data_dfof_tc(:,ic,ptrials)),2);
                    dfof_trialstd = std(squeeze(data_dfof_tc(:,ic,ptrials)),0,2);
                    dfof_trialsem = dfof_trialstd ./ sqrt(numtrials);
               subplot(12,5,s)
                   plot(X,dfof_trialavg, 'b','LineWidth',1);
                   shadedErrorBar(X,dfof_trialavg,dfof_trialsem);
                   xline(15,'--b')
                   ylim([-0.1 0.4])
                   set(gca,'TickDir','out'); box off; axis square; axis off
                   movegui('center')
                   s=s+1;
            end
        end
        sgtitle(['cell #' num2str(ic)])
        print(fullfile(base, 'Analysis\2P\CrossOri\timecourses', [date '_' mouse '_' run_str '_TCs_maxdist' num2str(max_dist) '_Cell' num2str(ic) '.pdf']), '-dpdf','-fillpage')
    end
end

