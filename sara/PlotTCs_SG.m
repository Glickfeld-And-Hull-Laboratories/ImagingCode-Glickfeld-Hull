
%% plot timecourses of all trials per cell (one cell per pdf page)

clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandPhase_15Hz_ExptList_SG';

rc = behavConstsAV;
eval(ds)
nexp = size(expt,2);

frame_rate = 15;

for iexp = [47]; 
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
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']))

    clear data_tc data_dfof_con_ph_tc_avg base_cell np_tc npSub_tc Val

    fprintf(['\nMouse: ' mouse '\nDate: ' date '\n'])
    fprintf(['nCells =' num2str(nCells) '\n'])
    fprintf(['nCells resp =' num2str(length(resp_ind)) '\n'])
    tm_ind = resptest_ind;
    x = setdiff(respmask_ind, resptest_ind);
    tm_ind = [tm_ind; x];
    fprintf(['nCells t or m resp =' num2str(length(tm_ind)) '\n'])
    fprintf(['nCells plaid resp =' num2str(length(respplaid_ind)) '\n'])
    

%         
%     phases = [0 45 90 135 180 225 270 315];
%     s = 1;
%     for p = 1:1:8
%         fprintf(['mask phase ' num2str(phases(s)) ', ' num2str(size(resp_all{2,2,p},2)) ' trials \n'])
%         s=s+1;
%     end



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

    for im = 1:nMaskCon
        ind_mask = find(maskCon_all == maskCons(im));
        for it = 1:nStimCon
            ind_stim = find(stimCon_all == stimCons(it));
            if it>1 & im>1
                for ip = 1:nMaskPhas
                    ind_phas = find(maskPhas_all == maskPhas(ip));
                    ind = intersect(ind_phas, intersect(ind_stim,ind_mask));
                    trialInd{im,it,ip} = ind;
                end
            end
        end
    end


% pull trial numbers for plaid phases and exclude bad pupil trials


    max_dist = 2;
    % for ip = 1:length(maskPhas)
    %     ptrials = trialInd{2,2,ip}; %set which specific plaid (1-8) you want to look at
    %     ptrials = setdiff(ptrials, find(centroid_dist<max_dist));
    % end


% How many cells are significantly responsive to plaids?
    fprintf(['\nNumber of cells significantly responsive to plaids: ' num2str(length(respplaid_ind)) '/' num2str(nCells) '\nListed in var respplaid_ind\n\n']);


% Plot individual and average time courses for all phases for cells

    vsize = numel(respplaid_ind);
    idx = randperm(vsize);
    rand_ind = respplaid_ind(idx(1:20));


    for ic = [rand_ind']; %cell number % rand_ind(:)
        figure;
        s=1;
        for ip = 1:nMaskPhas
            ptrials = trialInd{2,2,ip}; %set which specific plaid (1-8) you want to look at
            ptrials = setdiff(ptrials, find(centroid_dist<max_dist));
            for it = 1:length(ptrials)
                X = [1:75]; %taken from how data_dfof_tc is calculated in crossori_expt script (should equal on + off frames)
                subplot(4,2,s)
                plot(X,squeeze(data_dfof_tc(:,ic,ptrials(it))));
                dfof_trialavg = mean(data_dfof_tc(:,ic,ptrials),3);
                plot(X,dfof_trialavg, 'b','LineWidth',2);
                xline(15,'--b')
                xlabel('frames')
                ylabel('df/f')
                title(['cell #' num2str(ic), ', mask phase ' num2str(ip)])
                hold on
            end
           s=s+1;
        end
        print(fullfile(base, 'Analysis\2P\CrossOri\timecourses', [date '_' mouse '_' run_str '_TCs_maxdist' num2str(max_dist) '_Cell' num2str(ic) '.pdf']), '-dpdf','-fillpage')
        close all
    end
end


%% WORKING ON THIS; plot timecourses for preferred plaid phase for each cell (one cell per pdf page)

clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandPhase_15Hz_ExptList_SG';

rc = behavConstsAV;
eval(ds)
nexp = size(expt,2);

frame_rate = 15;

for iexp = 9 
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc{1};
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
    nrun = length(ImgFolder);
    run_str = 'runs-002';

    base = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara'];

    fprintf([mouse ' ' date '\n']);

% Load timecourse data, stim data
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'nCells', 'sz')
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']), 'nTrials', 'stimCon_all', 'nMaskPhas')
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']), 'centroid_dist')


    fprintf(['Mouse: ' mouse '\nDate: ' date '\n'])
    fprintf(['nCells = ' num2str(nCells) '\n'])
    fprintf(['nCells resp = ' num2str(length(resp_ind)) '\n'])
    fprintf(['nCells plaid resp = ' num2str(length(respplaid_ind)) '\n'])


for iDir = 1:nDir
    ind = find(dir_mat == dirs(iDir));
    data_dfof_dir(:,:,iDir) = mean(data_dfof(:,:,ind),3);
    dir_resp_avg(:,iDir,1) = squeeze(mean(mean(data_dfof(resp_win,:,ind),1),3));
    dir_resp_avg(:,iDir,2) = squeeze(std(mean(data_dfof(resp_win,:,ind),1),[],3))./sqrt(length(ind));
    dir_resp_mat = [dir_resp_mat; squeeze(mean(data_dfof(resp_win,:,ind),1))'];
    dir_list = [dir_list; iDir.*ones(length(ind),1)];
    [h_dir(:,iDir), p_dir(:,iDir)] = ttest(squeeze(mean(data_dfof(resp_win,:,ind),1)), squeeze(mean(data_dfof(base_win,:,ind),1)),'dim', 2, 'tail', 'right', 'alpha', 0.05./(nDir-1));
end



for im = 1:nMaskCon
        ind_mask = find(maskCon_all == maskCons(im));
        for it = 1:nStimCon
            ind_stim = find(stimCon_all == stimCons(it));
            if it>1 & im>1
                for ip = 1:nMaskPhas
                    ind_phas = find(maskPhas_all == maskPhas(ip));
                    ind = intersect(ind_phas, intersect(ind_f,intersect(ind_stim,ind_mask)));
                    trialsperstim(im,it,ip,iff) = length(ind);
                    resp_cell{im,it,ip} = squeeze(mean(data_dfof_tc(resp_win,:,ind),1));
                    base_cell{im,it,ip} = squeeze(mean(data_dfof_tc(base_win,:,ind),1));
                    data_dfof_con_ph_tc_avg(:,:,im,it,ip,iff,1) = squeeze(nanmean(data_dfof_tc(:,:,ind),3));
                    data_dfof_con_ph_tc_avg(:,:,im,it,ip,iff,2) = squeeze(nanstd(data_dfof_tc(:,:,ind),[],3)./sqrt(length(ind)));
                    trialInd{im,it,ip,iff} = ind; 
                end
            end
        end 
end


p_plaidanova = zeros(1,nCells); 

for iDir = 1:nDir
    ind = find(dir_mat == dirs(iDir));
%     data_dfof_dir(:,:,iDir) = mean(data_dfof(:,:,ind),3);
%     dir_resp_avg(:,iDir,1) = squeeze(mean(mean(data_dfof(resp_win,:,ind),1),3));
%     dir_resp_avg(:,iDir,2) = squeeze(std(mean(data_dfof(resp_win,:,ind),1),[],3))./sqrt(length(ind));
    dir_resp_mat = [dir_resp_mat; squeeze(mean(data_dfof(resp_win,:,ind),1))'];
%     dir_list = [dir_list; iDir.*ones(length(ind),1)];
%     [h_dir(:,iDir), p_dir(:,iDir)] = ttest(squeeze(mean(data_dfof(resp_win,:,ind),1)), squeeze(mean(data_dfof(base_win,:,ind),1)),'dim', 2, 'tail', 'right', 'alpha', 0.05./(nDir-1));
end


p_diranova(iCell) = anova1(dir_resp_mat(:,iCell),dir_list,'off');

end