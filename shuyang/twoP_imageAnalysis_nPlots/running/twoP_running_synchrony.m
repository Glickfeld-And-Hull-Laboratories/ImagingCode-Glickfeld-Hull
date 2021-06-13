% the purpose of the script is to see if cells are more/less synchronized
% when behavior state changes. See paper "On-beam synchrony in the cerebellum as the mechanism for the timing and oordination of movement, 2009 PNAS"
% synchrony is basically calculated by doing pair-wise (of cells)
% correlation within each behavioral state, and if the correlation of these
% pairs increase/decrese, synchrony increases/decreases
%% SECTION - assign pathnames and datasets to be analyzed/written. 
clear;
sessions = {'190429_img1021','190430_img1023','190507_img1024','190603_img1025'};
days = {'1021-190429_1','1023-190430_1','1024-190507_1','1025-190603_1'};

%sessions = {'190429_img1021','190507_img1024','190603_img1025'};
%days = {'1021-190429_1','1024-190507_1','1025-190603_1'};
%there might be more than 1 sessions on a single subject on the same day
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder

% behavior analysis results 
color_code = {'r','g','m','y','b'};

%% for each session, 1s running and 1s in stay

for i = 1:length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{i}, '\'];
    dfOvF_output = load([image_analysis_dest sessions{i} '_dfOvF.mat']);
    dfOvF = dfOvF_output.dfOvF_btm_cl;
    threshold = -4;% the threshold you used in deconvolveCa function
    spk_deconv_output = load([image_analysis_dest sessions{i},'_spk_deconvolve_threshold' num2str(threshold) '.mat']);
    spk_logic_cl = spk_deconv_output.spk_logic_cl; %frames*cells
    
    behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days{i} '\'];
    behav_output = load([behav_dest days{i} '_behavAnalysis.mat']);
    frames_runoff_mat = behav_output.frms_runoff_mat;%frames*trials
    
    tlen = 15; % how long the time bin is
    run = zeros(tlen,size(frames_runoff_mat,2));
    stay = zeros(tlen,size(frames_runoff_mat,2));
    % generate frame matrices for last n=tlen frames of running and n = tlen frames starting from 1s after running ends
    for n = 1: size(frames_runoff_mat,2)                                        % frames*trials
        run(:,n) = frames_runoff_mat(30,n)-(tlen-1): frames_runoff_mat(30,n);         % the 30th frame is the last frame of each running trial
        stay(:,n) = frames_runoff_mat(30,n)+30:frames_runoff_mat(30,n)+30+tlen-1;
    end
    
    % generate spike logic matrices    
    spk_run = zeros(size(run,1),size(run,2),size(spk_logic_cl,2)); %frames*trials*cells
    spk_stay = zeros(size(stay,1),size(stay,2),size(spk_logic_cl,2)); 
    for c = 1:size(spk_logic_cl,2)
        spk_logic_c = spk_logic_cl(:,c);
        spk_run(:,:,c) = spk_logic_c(run);
        spk_stay(:,:,c) = spk_logic_c(stay);
    end
    
    % generate dfOvF logic matrices    
    dfOvF_run = zeros(size(run,1),size(run,2),size(dfOvF,2)); %frames*trials*cells
    dfOvF_stay = zeros(size(stay,1),size(stay,2),size(dfOvF,2)); 
    for c = 1:size(dfOvF,2)
        dfOvF_c = dfOvF(:,c);
        dfOvF_run(:,:,c) = dfOvF_c(run);
        dfOvF_stay(:,:,c) = dfOvF_c(stay);
    end
    
    %=================================================================================================================
    pairs = nchoosek((1:size(dfOvF,2)),2); %pairs of cells, pair*2
    %raw correlation: do correlation for each pair of cells in each trial. conditions:0.5s before running ends, 1-1.5s after running ends, dfOvF --> spk
    r_dfOvF_run      = zeros(size(pairs,1),size(dfOvF_run,2)); %pairs*trials
    pval_dfOvF_run   = zeros(size(pairs,1),size(dfOvF_run,2));
    r_dfOvF_stay    = zeros(size(pairs,1),size(dfOvF_run,2));
    pval_dfOvF_stay = zeros(size(pairs,1),size(dfOvF_run,2));
    
    r_spk_run        = zeros(size(pairs,1),size(spk_run,2)); %pairs*trials
    pval_spk_run     = zeros(size(pairs,1),size(spk_run,2));
    r_spk_stay    = zeros(size(pairs,1),size(dfOvF_run,2));
    pval_spk_stay = zeros(size(pairs,1),size(dfOvF_run,2));
    
    for p = 1:size(pairs,1) %for each pair of cells
        for t = 1:size(dfOvF_run,2) %for each trial
            [r_dfOvF_run(p,t),pval_dfOvF_run(p,t)] = corr(dfOvF_run(:,t,pairs(p,1)),dfOvF_run(:,t,pairs(p,2)));
            [r_dfOvF_stay(p,t),pval_dfOvF_stay(p,t)] = corr(dfOvF_stay(:,t,pairs(p,1)),dfOvF_stay(:,t,pairs(p,2)));
            
            [r_spk_run(p,t),pval_spk_run(p,t)] = corr(spk_run(:,t,pairs(p,1)),spk_run(:,t,pairs(p,2)));
            [r_spk_stay(p,t),pval_spk_stay(p,t)] = corr(spk_stay(:,t,pairs(p,1)),spk_stay(:,t,pairs(p,2)));
        end
    end
    %spike gives a lot of NaNs so use df/f
    
    r_dfOvF_run_pairs     = mean(r_dfOvF_run,2);
    r_dfOvF_stay_pairs    = mean(r_dfOvF_stay,2);
    
    rshiftdfOvF_run  = zeros(size(pairs,1),size(dfOvF_run,2)*(size(dfOvF_run,2)-1));
    pshift_dfOvF_run = zeros(size(pairs,1),size(dfOvF_run,2)*(size(dfOvF_run,2)-1));%pairs*trials combinations, combinations = ntrials*ntrials-1 (C ntrial,2)
    rshiftdfOvF_stay = zeros(size(pairs,1),size(dfOvF_run,2)*(size(dfOvF_run,2)-1));
    pshiftdfOvF_stay = zeros(size(pairs,1),size(dfOvF_run,2)*(size(dfOvF_run,2)-1));
    
    %shift predictor: this can take fucking long
    for p = 1:size(pairs,1)
        a = 0;
        
        for t1 = 1:size(dfOvF_run,2)
            for t2 = 1:size(dfOvF_run,2)
                if t1~=t2 %for cellA trial1, corr with cellB trial2-N
                    a = a+1;
              
                    [rshiftdfOvF_run(p,a),pshift_dfOvF_run(p,a)] = corr(dfOvF_run(:,t1,pairs(p,1)),dfOvF_run(:,t2,pairs(p,2)));
                    [rshiftdfOvF_stay(p,a),pshiftdfOvF_stay(p,a)] = corr(dfOvF_stay(:,t1,pairs(p,1)),dfOvF_stay(:,t2,pairs(p,2)));
                end
            end
        end
    end
   
    rshiftdfOvF_run_pairs    = mean(rshiftdfOvF_run,2);
    rshiftdfOvF_stay_pairs  = mean(rshiftdfOvF_stay,2);
    
    eCorr_dfOvF_run_pairs     = r_dfOvF_run_pairs - rshiftdfOvF_run_pairs;
    eCorr_dfOvF_stay_pairs    = r_dfOvF_stay_pairs - rshiftdfOvF_stay_pairs;
    eCorr = [eCorr_dfOvF_run_pairs,eCorr_dfOvF_stay_pairs];
    
    x = [1,2]; x_plot = repmat(x,size(eCorr,1),1);
    eCorr_fig = figure;
    plot(x_plot',eCorr','.','LineStyle','-','linewidth', 1.25,'MarkerSize',20,'color',[0.5294 0.5294 0.5294]);hold on;
    fast_errbar(x,eCorr,1,'color',[0.7922 0 0.1255]);
    xlabel ('behavioral state');xlim([0.5 2.5]); 
    axis square;
    ylim([-1 1]);
    set(gca,'XTick',x,'XTicklabel',{'run','stationary'});
    ylabel('correlation between pairs of cells');
    title(['correlation ' sessions{i}]);
    savefig([image_analysis_dest sessions{i}, '_synch_corr_tlen' num2str(tlen)]);
    
    save([image_analysis_dest sessions{i} '_dfOvF.mat'] ,'eCorr','tlen','-append');
end

%{
%% across sessions
eCorr_all = [];
for i = 1:length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{i}, '\'];
    dfOvF_output = load([image_analysis_dest sessions{i} '_dfOvF.mat']);
    eCorr = dfOvF_output.eCorr;
    eCorr_all = cat(1,eCorr_all,eCorr);
end

x = [1,2]; x_plot = repmat(x,size(eCorr_all,1),1);
eCorr_fig = figure;
plot(x_plot',eCorr_all','.','LineStyle','none','MarkerSize',10,'color',[0.5294 0.5294 0.5294]);hold on;
fast_errbar(x,eCorr_all,1,'color',[0.7922 0 0.1255]);
xlabel ('behavioral state');xlim([0.5 2.5]);
axis square;
ylim([-1 1]);
set(gca,'XTick',x,'XTicklabel',{'run','stationary'});
ylabel('correlation between pairs of cells');
title('correlation across sessions');
savefig(['Z:\Analysis\2P_MovingDots_Analysis\across_sessions\' 'across_sessions_synch_corr.fig']);
%}

