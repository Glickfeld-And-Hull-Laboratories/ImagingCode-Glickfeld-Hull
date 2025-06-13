%% Plotting timecourses, gettting preferred phase...

clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandPhase_15Hz_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = size(expt,2);
frame_rate = 15;

%% plot timecourses with shaded error bar; one cell per page

clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandPhase_15Hz_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = size(expt,2);
frame_rate = 15;

for iexp = [92]; 
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
    fprintf(['nCells = ' num2str(nCells) '\n'])
    fprintf(['nCells resp = ' num2str(length(resp_ind)) '\n'])
    tm_ind = resptest_ind;
    x = setdiff(respmask_ind, resptest_ind);
    tm_ind = [tm_ind; x];
    fprintf(['nCells t or m resp = ' num2str(length(tm_ind)) '\n'])
    fprintf(['nCells plaid resp = ' num2str(length(respplaid_ind)) '\n'])

    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits.mat']))

    max_dist = 4;
    
    nTrials_pha = zeros(1,nMaskPhas);
    for i = 1:nMaskPhas
       [memb ind] = ismember(trialInd{2,2,i},find(centroid_dist<max_dist));
       nTrials_pha(1,i) = sum(ind~=0);
    end
    
    fprintf(['Trials per mask phase = ' num2str(nTrials_pha) '\n']) 

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

% Plot individual and average time courses for all phases for cells

    vsize = numel(respplaid_ind);
    idx = randperm(vsize);
    rand_ind = respplaid_ind(idx(1:15));
%     figure;
%     s=1;
    for ic = [1:15]; %cell number %[rand_ind']
        figure;
        s=1;
        for ip = 1:nMaskPhas
            ptrials = trialInd{2,2,ip}; %set which specific plaid (1-8) you want to look at
            ptrials = setdiff(ptrials, find(centroid_dist>max_dist));
                X = [1:75]; %taken from how data_dfof_tc is calculated in crossori_expt script (should equal on + off frames)
                dfof_trialavg = mean(data_dfof_tc(:,ic,ptrials),3);
                numtrials = size(squeeze(data_dfof_tc(:,ic,ptrials)),2);
                dfof_trialstd = std(squeeze(data_dfof_tc(:,ic,ptrials)),0,2);
                dfof_trialsem = dfof_trialstd ./ sqrt(numtrials);
           subplot(5,8,ip)
               plot(X,dfof_trialavg, 'b','LineWidth',1);
               shadedErrorBar(X,dfof_trialavg,dfof_trialsem);
               xline(15,'--b')
                ylim([-0.2 1])
               xlabel('frames')
               ylabel('df/f')
               set(gca,'TickDir','out'); box off; 
               title(['cell #' num2str(ic), ' ' num2str(ip)])
               movegui('center')
               s=s+1;
        end
        print(fullfile(base, 'Analysis\2P\CrossOri\timecourses', [date '_' mouse '_' run_str '_TCs_maxdist' num2str(max_dist) '_Cell' num2str(ic) '.pdf']), '-dpdf','-fillpage')
        close all
    end
end


%% plot timecourses of all trials per cell (one cell per pdf page)

for iexp = [85]; 
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
    fprintf(['nCells = ' num2str(nCells) '\n'])
    fprintf(['nCells resp = ' num2str(length(resp_ind)) '\n'])
    tm_ind = resptest_ind;
    x = setdiff(respmask_ind, resptest_ind);
    tm_ind = [tm_ind; x];
    fprintf(['nCells t or m resp = ' num2str(length(tm_ind)) '\n'])
    fprintf(['nCells plaid resp = ' num2str(length(respplaid_ind)) '\n'])

    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits.mat']))

    max_dist = 2;
    
    nTrials_pha = zeros(1,nMaskPhas);
    for i = 1:nMaskPhas
       [memb ind] = ismember(trialInd{2,2,i},find(centroid_dist<max_dist));
       nTrials_pha(1,i) = sum(ind~=0);
    end
    
    fprintf(['Trials per mask phase = ' num2str(nTrials_pha) '\n']) 

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

% Plot individual and average time courses for all phases for cells

    vsize = numel(respplaid_ind);
    idx = randperm(vsize);
    rand_ind = respplaid_ind(idx(1:15));
    
    for ic = [16 17]; %cell number %[rand_ind']
        figure;
        s=1;
        for ip = 1:nMaskPhas
            ptrials = trialInd{2,2,ip}; %set which specific plaid (1-8) you want to look at
            ptrials = setdiff(ptrials, find(centroid_dist>max_dist));
            for it = 1:length(ptrials)
                X = [1:75]; %taken from how data_dfof_tc is calculated in crossori_expt script (should equal on + off frames)
                subplot(4,2,s)
                dfof_trialavg = mean(data_dfof_tc(:,ic,ptrials),3);
                plot(X,squeeze(data_dfof_tc(:,ic,ptrials(it))));
                hold on
            end
           plot(X,dfof_trialavg, 'b','LineWidth',2);
           xline(15,'--b')
%             ylim([-1 6])
           xlabel('frames')
           ylabel('df/f')
           title(['cell #' num2str(ic), ', mask phase ' num2str(ip)])
           s=s+1;
        end
        print(fullfile(base, 'Analysis\2P\CrossOri\timecourses', [date '_' mouse '_' run_str '_TCs_maxdist' num2str(max_dist) '_Cell' num2str(ic) '.pdf']), '-dpdf','-fillpage')
        close all
    end
end

figure;
plot([1:length(centroid_dist)],centroid_dist)
hold on
h = refline(0,2);
h.Color = 'r';
h.LineWidth = 2;
% 
% figure;
% for ip = [1]
%     avg_dfof_t_all = [];
%     I_all = [];
%     for iexp = [8]
%         mouse = expt(iexp).mouse; date = expt(iexp).date; area = expt(iexp).img_loc{1}; ImgFolder = expt(iexp).coFolder; time = expt(iexp).coTime;
%         nrun = length(ImgFolder); run_str = 'runs-002'; base = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara'];
% 
%         % Load timecourse data, stim data
%         load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
%         load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
%         load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
%         load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
%         load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']))
%         load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits.mat']))
%         clear data_tc data_dfof_con_ph_tc_avg base_cell np_tc npSub_tc Val
%         
%         for im = 1:nMaskCon
%             ind_mask = find(maskCon_all == maskCons(im));
%             for it = 1:nStimCon
%                 ind_stim = find(stimCon_all == stimCons(it));
%                 if it>1 & im>1
%                     for ipp = 1:nMaskPhas
%                         ind_phas = find(maskPhas_all == maskPhas(ipp));
%                         ind = intersect(ind_phas, intersect(ind_stim,ind_mask));
%                         trialInd{im,it,ipp} = ind;
%                     end
%                 end
%             end
%         end
%         totCells(iexp,:) = nCells;
%         
%         ptrials = trialInd{2,2,ip}; %set which specific plaid (1-8) you want to look at
%         ptrials = setdiff(ptrials, find(centroid_dist>max_dist));
%         avg_dfof = mean(data_dfof_tc(:,:,ptrials),3);
%         avg_dfof_t = avg_dfof';
%         [B,I] = sort(pp_ind);
%         
%         I_all = [I_all; I+sum(totCells(1:iexp-1,:),1)];
%         avg_dfof_t_all = [avg_dfof_t_all; avg_dfof_t+sum(totCells(1:iexp-1,:),1)];
%     end
%         avg_dfof_sorted = avg_dfof_t_all(I_all,:);
% %         avg_dfof_sorted_norm = avg_dfof_sorted./max(max(avg_dfof_sorted,[],2));
%         
%         subplot(3,4,ip)
%             imagesc(avg_dfof_sorted(:,10:30))
%             hold on
% end




%% Plot timetraces by neuron ID (sorted by preferred phase)
clc; clear all; close all;
ds = 'CrossOriRandPhase_15Hz_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = size(expt,2);
frame_rate = 15;
max_dist = 2;

avg_dfof_t_all = [];
pp_ind_all = [];
lat_all_all = [];

phaseresp_all = cell(8,1);
prefp = [];
prefp_all = [];
amp_all_all = [];
latency_all = cell(8,1);

for iexp = [8 10 14 42 43 52 53 54 58]
    mouse = expt(iexp).mouse; date = expt(iexp).date; area = expt(iexp).img_loc{1}; ImgFolder = expt(iexp).coFolder; time = expt(iexp).coTime;
    nrun = length(ImgFolder); run_str = 'runs-002'; base = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara'];

    % Load timecourse data, stim data
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits.mat']))
    clear data_tc data_dfof_con_ph_tc_avg base_cell np_tc npSub_tc Val

    for im = 1:nMaskCon
        ind_mask = find(maskCon_all == maskCons(im));
        for it = 1:nStimCon
            ind_stim = find(stimCon_all == stimCons(it));
            if it>1 & im>1
                for ipp = 1:nMaskPhas
                    ind_phas = find(maskPhas_all == maskPhas(ipp));
                    ind = intersect(ind_phas, intersect(ind_stim,ind_mask));
                    trialInd{im,it,ipp} = ind;
                end
            end
        end
    end
    
    for ip = 1:8
        ptrials = trialInd{2,2,ip}; %set which specific plaid (1-8) you want to look at
        ptrials = setdiff(ptrials, find(centroid_dist>max_dist));
        avg_dfof = mean(data_dfof_tc(:,:,ptrials),3);
        avg_dfof_t = avg_dfof';
        maxdfof = max(max(avg_dfof_t, [], 2));
    %     norm_avg_dfof_t = avg_dfof_t./maxdfof;  %avg_dfof_t./maxdfof
%         avg_dfof_t_all = [avg_dfof_t_all; avg_dfof_t(resp_ind,:)];
%         pp_ind_all = [pp_ind_all; pp_ind(resp_ind)];

        %code for finding 1/2 time to peak for first cycle
        lat_dfof = avg_dfof_t(resp_ind,16:30); %avg_dfof_t(resp_ind,31:45);
        lat_max = max(lat_dfof, [], 2);
        lat_I_all = [];
        for c = 1:length(lat_max)
            [minDist,lat_I] = min(abs(lat_dfof(c,:)-(lat_max(c)/2)));
            lat_I_all = [lat_I_all; lat_I];
        end
        
        phaseresp_all{ip} =[phaseresp_all{ip}; avg_dfof_t(resp_ind,:)];
        latency_all{ip} = [latency_all{ip}; lat_I_all];
    end
    prefp = pp_ind(resp_ind);
    prefp_all = [prefp_all; prefp];
    amp_all_all = [amp_all_all; amp_hat_all(resp_ind)];
end

phaseresp_mat = permute(cat(3,phaseresp_all{:}),[3 1 2]);
latency_mat = permute(squeeze(cat(3,latency_all{:})),[2 1]);

stop
for ip = 1:8
    [B,I] = sort(prefp_all);
    
    figure(1)
    subplot(2,4,ip)
        imagesc(squeeze(phaseresp_mat(ip,I,:)))
        title(strcat('Mask phase ', num2str(ip)))
        hold on
end 

for ip = 1:8
    
    figure(2)    
    subplot(4,4,ip)
        c = [0 0.4470 0.7410];        
        histogram(squeeze(latency_mat(ip,prefp_all==ip)),7, 'FaceColor',c)
        title(strcat('Mask phase ', num2str(ip)))
        hold on
        if ip <= 4
            c = [0.6350 0.0780 0.1840];
            histogram(squeeze(latency_mat(ip,prefp_all==ip+4)),7, 'FaceColor',c)
        else
            c = [0.6350 0.0780 0.1840];
            histogram(squeeze(latency_mat(ip,prefp_all==ip-4)),7, 'FaceColor',c)
        end
        xlabel('Time to half peak (frames)')
    subplot(4,4,ip+8)
        c = [0 0.4470 0.7410];        
        cdfplot(squeeze(latency_mat(ip,prefp_all==ip)))
        title(strcat('Mask phase ', num2str(ip)))
        hold on
        if ip <= 4
            c = [0.6350 0.0780 0.1840];
            cdfplot(squeeze(latency_mat(ip,prefp_all==ip+4)))
        else
            c = [0.6350 0.0780 0.1840];
            cdfplot(squeeze(latency_mat(ip,prefp_all==ip-4)))
        end
        xlabel('Time to half peak (frames)')
        legend({'pref','anti'})
end


for ip = 1:8
    figure(3)
    subplot(2,2,1)
        tcs = cell2mat(phaseresp_all(ip,:,:));
        lats = cell2mat(squeeze(latency_all(ip,:)));
        scatter(lats, mean(tcs,2))
        hold on
        title('Time to half peak')
        xlabel('Time (frames)')
        ylabel('avg df/f')
    
    subplot(2,2,2)
        tcs = cell2mat(phaseresp_all(ip,:,:));
        tcs = tcs(:,16:end);
        lats = cell2mat(squeeze(latency_all(ip,:)));
        scatter(lats, max(tcs,[],2))
        hold on
        title('Time to half peak')
        xlabel('Time (frames)')
        ylabel('max df/f')
   
    subplot(2,2,3)
        edges = [10:5:60];
        avgtcs = mean(tcs,2);
        [N, e, bin] = histcounts(lats, edges);
        avgamp = [];
        avglat = [];
        for i = 1:length(edges)  
            amp = mean(avgtcs(bin==i));
            avgamp = [avgamp; amp];
            lat = mean(lats(bin==i));
            avglat = [avglat; lat];
        end
        scatter(avglat,avgamp,25,c)
        hold on
        ylim([-0.5 2])
        text(avglat(1:(end-1)),avgamp(1:(end-1)),sprintfc('  %d',N))
        title('Time to half peak, binned by latency (5 frames)')
end


figure(4); %plot average latency half peak
C = {[0.877 0.877 0.430],[0.707 0.906 0.465],[0.465 0.816 0.512],[0.246 0.715 0.551],[0 0.609 0.559],[0 0.496 0.523],[0.109 0.387 0.449],[0.164 0.281 0.344]};
avglat_byprefp = mean(latency_mat);
for ipp = 1:length(unique(prefp_all))
    subplot(2,1,1)
        h(ipp,1) = cdfplot(avglat_byprefp(find(prefp_all==ipp)));
        set(h(ipp,1), 'Color',C{ipp})
        hold on
end
    legend('1','2','3','4','5','6','7','8')
    set(h(:,1), 'LineWidth',1.5)
    set(gca,'TickDir','out');
    title('Latency by preferred plaid phase')
    

figure(5);
    edges = [0:3:15];
    [N, e, bin] = histcounts(avglat_byprefp, edges);
    avgMI = [];
    std_MI = [];
    sem_MI = [];
    for i = 1:length(edges)
        MI = mean(amp_all_all(find(bin==i)));
        stdev = std(amp_all_all(find(bin==i)));
        sem = stdev ./ sqrt(length(amp_all_all(find(bin==i))));
        avgMI = [avgMI; MI];
        std_MI = [std_MI; stdev];
        sem_MI = [sem_MI; sem];
    end
    subplot(2,1,1)
    scatter(edges,avgMI',25)
    hold on
    errorbar(edges,avgMI,sem_MI,'Color',[.7 .7 .7],"LineStyle","none")
    text(edges(1:(end-1)),avgMI(1:(end-1)),sprintfc('  %d',N))
    ylabel('avg MI modulation amp')
    xlabel('latency to 1/2 max peak')
    ylim([0 0.5])
    xlim([1 15])
    xticks([0:3:15]) 

stop

figure(1); print(fullfile(base, 'Analysis\2P\CrossOri\RandPhaseSummary', ['randPhase_Raster_AL.pdf']), '-dpdf','-fillpage')
figure(2); print(fullfile(base, 'Analysis\2P\CrossOri\RandPhaseSummary', ['randPhase_HalfRiseLatency_allcycles_AL.pdf']), '-dpdf','-fillpage')
figure(3); print(fullfile(base, 'Analysis\2P\CrossOri\RandPhaseSummary', ['randPhase_HalfRiseLatency_ByAvgDfof_allcycles_AL.pdf']), '-dpdf','-fillpage')
figure(4); print(fullfile(base, 'Analysis\2P\CrossOri\RandPhaseSummary', ['randPhase_HalfRiseLatency_AvgByPreferredPhase_firstcycle_AL.pdf']), '-dpdf','-fillpage')
figure(5); print(fullfile(base, 'Analysis\2P\CrossOri\RandPhaseSummary', ['randPhase_CompareInvarianceByHalfRiseLatency_firstcycle_V1.pdf']), '-dpdf','-fillpage')


       
%% find avg df/f and STD, plot by experiment
clc; clear all; close all;
ds = 'CrossOriRandPhase_15Hz_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = size(expt,2);
frame_rate = 15;
max_dist = 2;

V1 = [8 10 14 42 43 52 53 54 58];
AL = [29 31 7 49 50 55 57 60];
LM = [44 45 46 47 48];
axons = [74 76];
allexp = [8 10 14 42 43 52 53 54 58 29 31 7 49 50 55 57 60 44 45 46 47 48 74 76];


avg_dfof_t_all = [];

m_all = [];
s_all = [];

for iexp = allexp
    mouse = expt(iexp).mouse; date = expt(iexp).date; area = expt(iexp).img_loc{1}; ImgFolder = expt(iexp).coFolder; time = expt(iexp).coTime;
    nrun = length(ImgFolder); run_str = 'runs-002'; base = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara'];

    % Load timecourse data, stim data
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits.mat']))
    clear data_tc data_dfof_con_ph_tc_avg base_cell np_tc npSub_tc Val

    for im = 1:nMaskCon
        ind_mask = find(maskCon_all == maskCons(im));
        for it = 1:nStimCon
            ind_stim = find(stimCon_all == stimCons(it));
            if it>1 & im>1
                for ipp = 1:nMaskPhas
                    ind_phas = find(maskPhas_all == maskPhas(ipp));
                    ind = intersect(ind_phas, intersect(ind_stim,ind_mask));
                    trialInd{im,it,ipp} = ind;
                end
            end
        end
    end
    phaseresp_all = {};
    for ip = 1:8
        ptrials = trialInd{2,2,ip}; %set which specific plaid (1-8) you want to look at
        ptrials = setdiff(ptrials, find(centroid_dist>max_dist));
        avg_dfof = mean(data_dfof_tc(:,:,ptrials),3);
        avg_dfof_t = avg_dfof';
        
        phaseresp_all =[phaseresp_all; avg_dfof_t(resp_ind,:)];
    end
    phaseresp_mat = permute(cat(3,phaseresp_all{:}),[3 1 2]);

    respall = reshape(phaseresp_mat,[],1);
    m = mean(respall);
    
    stdresp = squeeze(mean(phaseresp_mat,1));
    s_cells = [];
    stand = [];
        for i = 1:size(stdresp,1)
            stand = std(stdresp(i,:));
            s_cells = [s_cells;stand];
        end

    s = mean(s_cells);
    m_all = [m_all;m];
    s_all = [s_all;s];
end

figure(1);
scatter(m_all, s_all,50)
xlabel('avg df/f')
ylabel('avg STD df/f (within trial)')
title('df/f noise')
set(gca,'TickDir','out');
text(m_all(1:(end)),s_all(1:(end)),sprintfc('  %d',(allexp')))

figure(1); print(fullfile(base, 'Analysis\2P\CrossOri\RandPhaseSummary', ['randPhase_Noise_allexperiments.pdf']), '-dpdf','-fillpage')


