%% BUTTERWORTH HIGH-PASS FILTER USED ON TC_AVG

clear; close all;
analysis_out = ('A:\home\carlo\analysis\2P\');

expType=input('monomodal(1) or dimodal(2) session: ');
seshType=input('interleaved (1) or blocked (2) session: ');
if seshType==1 && expType == 1
RCExptListCarlo_Inter;
elseif seshType==2 && expType == 1
RCExptListCarlo_Block;
elseif seshType==1 && expType == 2
RCExptListCarlo_VA_Inter;
end

pairings = loadRCList_Carlo(expt,seshType);
iexp = pairings{1,1}; isesh = pairings{3,1};
mouse = strtrim(expt{1,iexp}.mouse{1,isesh});
date = expt{1,iexp}.date{1,isesh};
run = expt{1,iexp}.run{1,isesh};
img_fn = [date '_img' mouse '\getTC_' run '\'];
load([analysis_out, img_fn, date '_img' mouse '_' run 'saveOutputs.mat']);
load([analysis_out, img_fn, date '_img' mouse '_' run '_nPCA' num2str(nPCA) '_mu' num2str(mu) '_nIC' num2str(nIC) '_thresh' num2str(cluster_threshold) '_TCave.mat']);

cutoffF = .035; FrameRate = 30;
%high-pass filter the time course before replacing laser off periods with NaN
[b,a]=butter(3,(cutoffF*2/FrameRate),'high');
tc_avg = filtfilt(b,a,raw_tc_avg); %zero-phase filtering using parameters from butterworth high-pass filter to output filtered TC
%to test the alignment by removing the session mean from the unfiltered TC
baseOrigTC = nanmean(nanmean(raw_tc_avg(1:1000,:),2),1); %find average of entire TC
origTC = nanmean(raw_tc_avg,2);  %call TC collapsed across neurons
normTC = tc_avg+baseOrigTC; %remove mean from collapsed TC
%plot raw and filtered TC
plot(origTC,'r');
hold on;
plot(nanmean(normTC,2),'k'); %collapse filtered TC across neurons to plot
title(['cutoff frequency at ' num2str(cutoffF) ' Hz']);
legend('raw TC','filtered TC','Location','northwest');
save([analysis_out, img_fn, date '_img' mouse '_' run '_nPCA' num2str(nPCA) '_mu' num2str(mu) '_nIC' num2str(nIC) '_thresh' num2str(cluster_threshold) '_TCave.mat'], 'normTCAvg','normOrigTC','-append');


%% META K-MEANS ANALYSIS VIA TEMPORAL ALIGNMENT OF CLIMBING FIBER ACTIVITY TO PC DENDRITES
clear; close all;
analysis_out = ('A:\home\carlo\analysis\2P\');
bdata_source = 'A:\home\carlo\rawData\behavioral\';

expType=input('monomodal(1) or dimodal(2): '); expIdx = {'monomodal', 'dimodal'};
seshType=input('interleaved (1), blocked (2), or naive (3) session: ');
if seshType==1
    if expType==1; RCExptListAll_Inter; elseif expType==2; RCExptListCarlo_VA_Inter; end
sesh = 'PL'; sesh2=expIdx{1,expType};
elseif seshType==2
    if expType==1; RCExptListCarlo_Block; elseif expType==2; RCExptListCarlo_VA_Block; end
sesh = 'Block'; sesh2=expIdx{1,expType};
elseif seshType==3
    if expType==1; RCExptListAll_Naive; elseif expType==2; RCExptListCarlo_VA_Inter; end
sesh = 'Naive'; sesh2=expIdx{1,expType};
end

pairings = loadRCList_Carlo(expt,seshType);
iexp = pairings{1,1}; isesh = pairings{3,1};
mouse = strtrim(expt{1,iexp}.mouse{1,isesh});
date = expt{1,iexp}.date{1,isesh};
run = expt{1,iexp}.run{1,isesh};
img_fn = [date '_img' mouse '\kMeans\'];
img_fn2 = [date '_img' mouse '\getTC_' run '\'];
if ~exist(fullfile(analysis_out,img_fn))
    mkdir(fullfile(analysis_out, img_fn))
end
thresh = input('deconvolution threshold [-x.x]: ');
load([analysis_out date '_img' mouse '\getTC_' run '\' date '_img' mouse '_' run 'saveOutputs.mat']);
load([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_deconvolution_thresh' num2str(thresh) '_TCave_cl.mat']);
load([analysis_out img_fn2 date '_img' mouse '_' run '_nPCA' num2str(nPCA) '_mu' num2str(mu) '_nIC' num2str(nIC) '_thresh' num2str(cluster_threshold) '_TCave.mat']);

%create correlation matrix and determine correct values for kMeans
filtCorr = corrcoef(TCave_cl); figure; imagesc(filtCorr);
clusterCount = linkage(filtCorr,'ward'); 
co4HC = median([clusterCount(end-2,3) clusterCount(end-1,3)]);
figure; dendrogram(clusterCount,'ColorThreshold',co4HC);
co4HC = median([clusterCount(end-5,3) clusterCount(end-4,3)]);
figure; dendrogram(clusterCount,'ColorThreshold',co4HC);

%are initial clusters grouped based on basseline amplitude 
%(spoilers, probably)
c = cluster(clusterCount,'Maxclust',3);
colours = {[1 0.2 0.2], [0.5 1 0.5], [0.2 0.2 1]}; %[ r g b ]%
figure;
for i = 1:max(c)
    plot(nanmean(TCave_cl(:,c==i),2),'color',colours{i});
    hold on;
end
%what do these initial clusters look like spatially
load([analysis_out, img_fn2, date '_img' mouse '_' run '_nPCA' num2str(nPCA) '_mu' num2str(mu) '_nIC' num2str(nIC) '_thresh' num2str(cluster_threshold) '_coor' num2str(threshold) '_mask3D_final.mat']);
figure;
imshow([analysis_out,img_fn2 'AVG_' date '_img' mouse '_' run '_rgstr_tiff_every50_ref' num2str(ref) '_jpeg.jpg']); hold on;
for m = 1:max(c)
    for n = 1:length(c)
    if c(n,1) == m
        cD{1,m}(1,n) = n;
    else
        cD{1,m}(1,n) = NaN;
    end
    end
end
for m = 1:max(c)
cD{1,m}(find(isnan(cD{1,m})))=[];
end
for m = 1:max(c)
    for n = 1:size(cD{1,m},2)
    bound = cell2mat(bwboundaries(mask3D_final(:,:,cD{1,m}(1,n))));
    plot(bound(:,2),bound(:,1),'.','color',colours{m}); hold on;
    text(mean(bound(:,2)),mean(bound(:,1)), ...
    num2str(cD{1,m}(1,n)), 'color', 'k', 'FontSize', 8);
    hold on;
    end
end
%{ 
%NEED TO SEGMENT NEURAL DATA TO SELECTIVELY OUTPUT FRAMES THAT EXIST BETWEEN LASERON AND CUE_ONSET
    bxStartMWorksTime  = round(double(cell2mat(mworks.tThisTrialStartTimeMs(1)))); %first frame in mworks-time
    laserOnTrial = double(cell2mat(mworks.holdStartsMs)) - bxStartMWorksTime;  %find laser on trial times
    preCueInterval  = double(cell2mat(mworks.holdTimesMs)); %find time between laser on and cue onset (~2667)
    rewardOnsetTrial = laserOnTrial + preCueInterval; %calculate reward onset for each trial from start of session
    postCueInterval = double(cell2mat(mworks.reactTimesMs)); %time between cue and rew delivery (~767)
    cuePresentation = rewardOnsetTrial-postCueInterval; %cue presentation

    cTargetOn = mworks.cTargetOn;

    frameStartLaser = []; %finding frame number for the start of each trial
    for s = 1:length(laseron)-1
        if laseron(1,s+1)-laseron(1,s)>100
            frameStartLaser = [frameStartLaser s+1];
        else
        end
    end
    
    cueOnLaserTime = []; whichTrial = 1; %finding frame number for cue pres. of each trial
    for ss = 1:length(laseron)
        if laseron(1,ss) == cTargetOn(1,whichTrial)
            cueOnLaserTime = [cueOnLaserTime ss];
            whichTrial=whichTrial+1;
        elseif isnan(cTargetOn(1,whichTrial))
            whichTrial=whichTrial+1;
        end
    end
    
    preCuePeriod = []; %finding all frame numbers between start and cue for each trial
    for sss=1:length(cueOnLaserTime)
    preCuePeriod = [preCuePeriod (laseron(1,frameStartLaser(1,sss):cueOnLaserTime(1,sss)))];
    end
    
for ss = 1:size(TCave_cl,2) %now take all those pre-cue times and pull out neural data for each dendrite
    preCueTC(:,ss) = all_TCave_cl(preCuePeriod,ss);
end
%}
expr = TCave_cl';%preCueTC';
    clusterVal = size(clusterCount,2)+2;
    numden = size(TCave_cl,2);    %
    numtimes = size(TCave_cl,1);  % these 4 were taken from meta_k_means.m
    counts = zeros(numden);       % as variables req for the following lines
    centroidmat = [];             %
%THE FOLLOWING WAS TAKEN FROM meta_k_means.m - credit: Ozden et al., 2008
%MINOR ALTERATIONS MADE THROUGHOUT WITH COMMENTS
%for k=3 >> removed for loop for k (clusterVal here -#of presumed clusters in data based on hierarchical clustering - 'linkage' above)
   % time elapsed for 28920x205 TC array - ~7:27 minutes
for it = 1:1000
    try 

        [IDX,C,SUMD,D] = kmeans(expr,clusterVal,'Distance','correlation');
%IDX  - cluster assignment for each dendrite
%C    - centroid location
%SUMD - within-cluster sums of distance to centroid
%D    - distance from each dendrite to each centroid
        for x = 1:numden
            for y = 1:numden
                if IDX(x)==IDX(y)
                counts(x,y) = counts(x,y) + 1; %returns
                end
            end
        end
    catch
    IDX=[];
    end
end


if ~isempty(IDX)
thr = 0.8; %correlation threshold - calculation of number of times two dendrites were clustered/total cluster attempts
clusters = [];
for xx = 1:numden
    for yy = 1:xx-1
        if counts(xx, yy) >= thr*it
            clusters = [clusters xx];
            clusters = [clusters yy];
        end
    end
end

% sort list of dendrites
clusters = sort(clusters);
% delete the repeated dendrites from the list
check = clusters(1);
list = check;
for c = 2:length(clusters)
    if clusters(c) ~= check
        list = [list clusters(c)];
        check = clusters(c);
    end
end

noThrCell = [];
for z = 1:size(TCave_cl,2)
    if ~ismember(z,list)
        noThrCell = [noThrCell z]; %find cells that didn't have +800 occurances with another dendrite in kmeans
    end
end

load([analysis_out, img_fn2, date '_img' mouse '_' run '_nPCA' num2str(nPCA) '_mu' num2str(mu) '_nIC' num2str(nIC) '_thresh' num2str(cluster_threshold) '_coor' num2str(threshold) '_mask3D_final.mat']);
%colours = {[1 0.2 0.2],[0.5 0.184 0.556],[0 0.5 1],[0 0.6 0.3],[0.2 0.8 0.8],[1 0.7 0.1250]}; % r m b g c y
figure;
imshow([analysis_out,img_fn2 'AVG_' date '_img' mouse '_' run '_rgstr_tiff_every50_ref' num2str(ref) '_jpeg.jpg']); hold on;
for m  = 1:length(noThrCell)
    bound = cell2mat(bwboundaries(mask3D_final(:,:,noThrCell(1,m))));
    colours = rand(1,4);
    plot(bound(:,2),bound(:,1),'.','color',colours); hold on;
    text(mean(bound(:,2)),mean(bound(:,1)), ...
        num2str(noThrCell(1,m)), 'color', 'k', 'FontSize', 8);
end

% vector of average point to centroid distances for each cluster
numvect = [];
cnt=1;
dendmem=[];
dendmem=dendmem';
% find the clusters and print out graded memberships for the dendrites
for n = 1:length(list)
    item = list(n);
    if item ~= 0
        members = [item];
        for nn = 1:numden
            if nn ~= item
                if counts(nn, item) >= thr*it
                    members = [members nn];
                    list(find(list == nn)) = 0;
                end
            end
        end
        eval(['cluster' num2str(cnt) '=members;']);
        cnt=cnt+1;
        % find centroid for cluster
        centroid = zeros(1, numtimes);
        for m = 1:length(members)
            centroid = centroid + expr(members(m), :);
        end
        centroid = centroid/length(members);
        centroidmat = [centroidmat; centroid];
        numerator = 0;
        % calculate the validity numerator
        for m = 1:length(members)
            numerator = numerator + (sum((expr(members(m),:)-centroid).^2))^0.5;
        end
        numerator = numerator/length(members);
        numvect = [numvect numerator];
        %graded membership using correlations
        grmem = zeros(numden, 1);
        for g = 1:numden
            grmem(g, 1) = corr2(centroid, expr(g, :));
        end
        dendmem=[dendmem grmem];
        grmem = [(1:numden); grmem'];

    end
end

allclusters=cell(cnt-1,2);
for i=1:cnt-1
    eval(['tempcluster=cluster' num2str(i) ';']);
    allclusters{i,1}=tempcluster;
    allclusters{i,2}=centroidmat(i,:);
end
centroidcorr = corrcoef(centroidmat');
combclusters=0;
[newclusters, dunnsinitial, centroidcorr, dendmem]=ioclustervalidity2(allclusters, expr, combclusters);
%ioclustervalidity2 IS A PARTNER SCRIPT WRITTEN BY THOSE INVOLVED IN OZDEN ET AL., 2008 - NO CHANGES WERE MADE TO THAT SCRIPT
% time elapsed on filtered 28920x205 pool after returning to each dendrites
% initial fluorescence
endclustering=0;
numclusters=clusterVal; %value changed from 10 to clusterVal
while endclustering==0
    numclusters=size(allclusters,1);
    combineclusters=[];
    for i=1:numclusters
        for j=(i+1):numclusters
            [newclusters, dunns, centroidcorr, dendmem]=ioclustervalidity2(allclusters, expr, [i j]);
            if dunns>=dunnsinitial
                combineclusters=[i j];
                dunnsinitial=dunns;
            end
        end
    end
    if ~isempty(combineclusters)
        [allclusters, dunnsinitial, centroidcorr, dendmem]=ioclustervalidity2(allclusters, expr, [combineclusters(1) combineclusters(2)]);
    else
        endclustering=1;
    end
    if numclusters==3
        endclustering=1;
    end
end
end
% numclusters=size(allclusters,1);
% for zz=1:numclusters;
%     clusters=allclusters{zz,1};
%     for zzz=1:length(clusters);
%         clusters(zzz)=clusters(zzz)+sum(badPCs(1:clusters(zzz)));
%     end
%     allclusters{zz,1}=clusters;
% end
% end


%needs to modify the line below due to different file names
load([analysis_out, img_fn2, date '_img' mouse '_' run '_nPCA' num2str(nPCA) '_mu' num2str(mu) '_nIC' num2str(nIC) '_thresh' num2str(cluster_threshold) '_coor' num2str(threshold) '_mask3D_final.mat']);
colours = {[1 0.2 0.2],[0.5 0.184 0.556],[0 0.5 1],[0 0.6 0.3],[0.2 0.8 0.8],[1 0.7 0.1250],[0.4 0.8 0.3]}; 
          %[r m b g c y]%

figure;
imshow([analysis_out,img_fn2 'AVG_' date '_img' mouse '_' run '_rgstr_tiff_every50_ref' num2str(ref) '_jpeg.jpg']); hold on;
for m  = 1:size(allclusters,1)
    for n = 1:size(allclusters{m,1},2)
    bound = cell2mat(bwboundaries(mask3D_final(:,:,allclusters{m,1}(1,n))));
    plot(bound(:,2),bound(:,1),'.','color',colours{m}); hold on;
    text(mean(bound(:,2)),mean(bound(:,1)), ...
        num2str(allclusters{m,1}(1,n)), 'color', 'k', 'FontSize', 8);
    end
end
hold on;
title([date ' img' mouse ' mask for k-means']);
savefig([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_mask_wdendrites_origPreCueKMeans.fig']);

save([analysis_out, img_fn, 'filtered_denormPreCueKMeans.mat'],'cluster*','centroid*','badPCs','*clusters','counts','dunnsinitial','dendmem','goodcells','grmem','mask3D_final','noThrCell','thr');

%% quick TC plot
    filename = dir([analysis_out,date '_img' mouse '\' date '_img' mouse '_' run  '*' '_TCave_cl.mat']);
    load([analysis_out,date '_img' mouse '\', filename.name]);
    %threshold = -3;
    spk = load([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_spk_deconvolve_threshold' num2str(thresholdDeco) '.mat' ]);
    spk_logic = spk.spk_logic_cl;
    nIC = size(spk_logic,2);
    %% fill laser off period with nans
    all_events = spk_logic;
    if expt{1,iexp}.ttl(1,isesh)
        all_events_temp = nan(size(tc_avg_all,1),size(spk_logic,2));
        all_events_temp(laseron(1,1:end),:) = all_events;
        all_events = all_events_temp;
    end
    sessions = [date '_img' mouse];% for behavior data - format YYMMDD[2/0/1]_imgMOUSE 2 = interleaved; 0 = CS- block; 1 = CS+ block
    if strlength(sessions) < 16
        bfile = dir([bdata_source 'data-i' mouse '-' date '*' ]);
    else
        bfile = dir([bdata_source 'data-i' mouse '-' date sessions(end-4:end) '.mat']);
    end
%load behavior file
    behav_dest = [bdata_source bfile.name];
    assert(length(bfile)) = 1;
    b_data = load(behav_dest);
    mworks = b_data.input;
    
    react_time = double(cell2mat(mworks.reactTimesMs(1,2:end)));
    cue_rew_int = mworks.RewardDelayDurationMs + round(mean(react_time),-1); %"total interval between cueonset and reward delivery"       
    frameRateHz = double(mworks.frameRateHz);
    cTargetOn = (mworks.cTargetOn);
    nTrials = length(cTargetOn);
    prewin_frames = round(1500./frameRateHz);
    postwin_frames = round(3000./frameRateHz);
    nFrames=cTargetOn(1,end);
    for itrial = 1:nTrials
        if cTargetOn(itrial)+postwin_frames-1 <= nFrames %& input.counterValues{itrial}(end)-cTargetOn(itrial) > postwin_frames
            targetAlign_tc(:,:,itrial) = all_TCave_cl(cTargetOn(itrial)-prewin_frames:cTargetOn(itrial)+postwin_frames-1,:);
            targetAlign_events(:,:,itrial) = all_events(cTargetOn(itrial)-prewin_frames:cTargetOn(itrial)+postwin_frames-1,:);
        end
    end
    
    targetAlignF = nanmean(targetAlign_tc(1:prewin_frames,:,:),3); % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.
    targetAligndFoverF = zeros(size(targetAlign_tc,1),size(targetAlign_tc,2),size(targetAlign_tc,3));%frame*cell*trials
    % calculate df/F
    for c = 1:size(TCave_cl,2)
        targetAligndFoverF(:,c,:) = (targetAlign_tc(:,c,:)-targetAlignF(c))./targetAlignF(c);
    end
    if mworks.doBlock2 %Specific CS+/CS- experiment added by Mike
        block2Trial = celleqel2mat_padded(mworks.tBlock2TrialNumber(1,2:end));
        ind_block2 = find(block2Trial); 
        ind_rew = find(~block2Trial); 
    else
        ind_rew = intersect(find(omitRewTrial == 0),find(unexpRewTrial == 0));
        if intersect(ind_omit, ind_unexp)
            x = ismember(ind_omit,intersect(ind_omit, ind_unexp));
            ind_omit(x) = [];
            x = ismember(ind_unexp,intersect(ind_omit, ind_unexp));
            ind_unexp(x) = [];
        end
    end
    s=1; if mworks.doBlock2; s = s+1; end
    tt = (-prewin_frames:postwin_frames-1).*(1000./frameRateHz);   
for i = 1:length(allclusters)
    clusTC = targetAlign_events(:,:,allclusters{i,1});
    clusDFoF = targetAligndFoverF(:,:,allclusters{i,1});
     figure; n = 1;
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,allclusters{i,1},ind_rew),3),2), nanstd(nanmean(targetAligndFoverF(:,allclusters{i,1},ind_rew),3),[],2)./sqrt(nIC), 'k');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        %ylim([-3 5]);
        title('CS+')
        n = n+1;
        if mworks.doBlock2
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,allclusters{i,1},ind_block2),3),2), nanstd(nanmean(targetAligndFoverF(:,allclusters{i,1},ind_block2),3),[],2)./sqrt(nIC),'r');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        %ylim([-3 5]);
        title('CS-')
        end
        sgtitle([date ' ' mouse  ' ' num2str(i)])
        savefig(fullfile(analysis_out,img_fn, [num2str(i) '_cueAlign_dFoverF.fig']))

      figure; n=1;
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,allclusters{i,1},ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,allclusters{i,1},ind_rew),3),[],2)./sqrt(nIC)).*(1000./frameRateHz));
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([0 2]);
        title('CS+')
        n = n+1;
        if mworks.doBlock2
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,allclusters{i,1},ind_block2),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,allclusters{i,1},ind_block2),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'r');
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([0 2]);
        title('CS-')
        end
        sgtitle([date ' ' mouse ' ' num2str(i)])
        savefig(fullfile(analysis_out,img_fn, [num2str(i) '_cueAlign_events_Hz.fig'])) 
end

for i = 1:length(allclusters)
    plotTC_Jin(TCave_cl(:,allclusters{i,1}), mask3D, 0, 1:size(allclusters{i,1},2), frameRateHz);

end

%%
%}
%{
%% PS BAND ANALYSIS     
% next next: take each cluster and run secondary cue align and lick align
% analysis on them to get a better picture of the activity driven
% differences between bands (likely need piezo data too)    

clear; close all;
analysis_out = 'A:\home\carlo\analysis\2P\';
bdata_source = 'A:\home\carlo\rawData\behavioral\';

removing_bad_trials = false; %bad_trials = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]; %added to remove shitty trials, turn off if not needed

RCExptListMike_Inter
%RCExptListCarlo_Inter
%RCExptListCarlo_Block
pairings = loadRCList_Mike;%Carlo;

%for  j = 1:size(expt.ttl,2)

    id = pairings{1,1};
    exp_subset = id;
    iexp = exp_subset;  %1:nexp
    mouse = strtrim(expt(iexp).mouse);
    date = expt(iexp).date;
    run = expt(iexp).run;
    fprintf([date ' ' mouse ' ' run '\n']);
    img_fn = [date '_img' mouse '\getTC_' run '\'];
    img_fn2 = [date '_img' mouse '\kMeans\'];
    load([analysis_out, img_fn2, 'filtered_denormKMeans.mat']);
    load([analysis_out, img_fn, date '_img' mouse '_' run 'saveOutputs.mat']);
    thresholdDeco = input('deconvolution threshold [-x.x]: ');
    load([analysis_out '\' date '_img' mouse '\' date '_img' mouse '_' run '_deconvolution_thresh' num2str(thresholdDeco) '_TCave_cl.mat']);
    spk = load([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_spk_deconvolve_threshold' num2str(thresholdDeco) '.mat' ]);
    load([analysis_out, img_fn, date '_img' mouse '_' run, '_nPCA', num2str(nPCA), '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), '_TCave.mat']);

%%
for c = 1:size(allclusters,1)
    tempClusName=sprintf('cluster_%d',c);
    spk_logic = spk.spk_logic_cl(:,allclusters{c,1});
    nIC = size(spk_logic,2);

    all_events = spk_logic;
    if expt(iexp).ttl
        all_events_temp = nan(size(all_TCave_cl,1),size(spk_logic,2));
        all_events_temp(laseron(1,1:end),:) = all_events;
        all_events = all_events_temp;
    end
    
    %save(fullfile([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_spk_deconvolve_threshold' num2str(thresholdDeco) '.mat' ]), 'all_events', '-append')
    
    sessions = [date '_img' mouse];% for behavior data - format YYMMDD[2/0/1]_imgMOUSE 2 = interleaved; 0 = CS- block; 1 = CS+ block
    if strlength(sessions) < 16
        bfile = dir([bdata_source 'data-i' mouse '-' date '*' ]);
    else
        bfile = dir([bdata_source 'data-i' mouse '-' date sessions(end-4:end) '.mat']);
    end
%load behavior file
    behav_dest = [bdata_source bfile.name];
    assert(length(bfile)) = 1;
    b_data = load(behav_dest);
    mworks = b_data.input;
    
    react_time = double(cell2mat(mworks.reactTimesMs(1,2:end)));
    cue_rew_int = mworks.RewardDelayDurationMs + round(mean(react_time),-1); %"total interval between cueonset and reward delivery"       
    img_fn2 = sessions;
    frameRateHz = double(mworks.frameRateHz);
    cTargetOn = (mworks.cTargetOn);
    nTrials = length(cTargetOn);
    prewin_frames = round(1500./frameRateHz);
    postwin_frames = round(3000./frameRateHz);
    targetAlign_tc = nan(prewin_frames+postwin_frames,nIC,nTrials);
    targetAlign_events = nan(prewin_frames+postwin_frames,nIC,nTrials);
    nFrames = size(tc_avg_all,1);

%%
    tempTCA = all_TCave_cl(:,allclusters{c,1});
    tempEVT = all_events;
    for itrial = 1:nTrials
        if cTargetOn(itrial)+postwin_frames-1 <= nFrames %& input.counterValues{itrial}(end)-cTargetOn(itrial) > postwin_frames
            targetAlign_tc(:,:,itrial) = tempTCA(cTargetOn(itrial)-prewin_frames:cTargetOn(itrial)+postwin_frames-1,:);
            targetAlign_events(:,:,itrial) = tempEVT(cTargetOn(itrial)-prewin_frames:cTargetOn(itrial)+postwin_frames-1,:);
        end
    end
    
    targetAlignF = nanmean(targetAlign_tc(1:prewin_frames,:,:),3); % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.
    targetAligndFoverF = zeros(size(targetAlign_tc,1),size(targetAlign_tc,2),size(targetAlign_tc,3));%frame*cell*trials
    % calculate df/F
    for d = 1:size(TCave_cl(allclusters{c,1}),2)
        targetAligndFoverF(:,d,:) = (targetAlign_tc(:,d,:)-targetAlignF(d))./targetAlignF(d);
    end
    
    
        omitRewTrial = celleqel2mat_padded(mworks.tRewardOmissionTrial(1,2:end));
        ind_omit = find(omitRewTrial);
        
        unexpRewTrial = celleqel2mat_padded(mworks.tDoNoStimulusChange(1,2:end));
        ind_unexp = find(unexpRewTrial);
        
        if mworks.doBlock2 %Specific CS+/CS- experiment added by Mike
            block2Trial = celleqel2mat_padded(mworks.tBlock2TrialNumber(1,2:end));
            ind_block2 = find(block2Trial); 
            ind_rew = find(~block2Trial); 
        else
            ind_rew = intersect(find(omitRewTrial == 0),find(unexpRewTrial == 0));
            if intersect(ind_omit, ind_unexp)
                x = ismember(ind_omit,intersect(ind_omit, ind_unexp));
                ind_omit(x) = [];
                x = ismember(ind_unexp,intersect(ind_omit, ind_unexp));
                ind_unexp(x) = [];
            end
        end
        
        
        if removing_bad_trials
            for b = 1:length(bad_trials)
                ind_omit = ind_omit(ind_omit ~= bad_trials(b));
                ind_unexp = ind_unexp(ind_unexp ~= bad_trials(b));
                ind_rew = ind_rew(ind_rew ~= bad_trials(b));
                if mworks.doBlock2
                    ind_block2 = ind_block2(ind_block2 ~= bad_trials(b));
                end
            end
        end
        
        tt = (-prewin_frames:postwin_frames-1).*(1000./frameRateHz);
        rewardDelayDurationMs = double(max(celleqel2mat_padded(mworks.tRewardDelayDurationMs(1,2:end)),[],2));
        reactTimesMs = double(mean(celleqel2mat_padded(mworks.reactTimesMs(1,2:end)),2));
        delayTimeMs = reactTimesMs+rewardDelayDurationMs;
        
        save([analysis_out, img_fn, 'figDataUntrimmed_' tempClusName '.mat']);
%% PLOTS    
        
        s=1;
        if length(ind_omit>5); s = s+1; end
        if length(ind_unexp>5); s = s+1; end
        if mworks.doBlock2; s = s+1; end
        
        if unique(expt(exp_subset).name == 'Block')
            indivFig = true;
        else
            indivFig = false;
        end
            
        if ~indivFig 
        figure; n = 1;
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_rew),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_rew),3),[],2)./sqrt(nIC), 'k');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        ylim([-10 10]);
        title('CS+')
        n = n+1;
        if length(ind_omit>5)
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_omit),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_omit),3),[],2)./sqrt(nIC),'r');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        title('Omit')
        n = n+1;
        end
        if length(ind_unexp>5)
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_unexp),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_unexp),3),[],2)./sqrt(nIC),'g');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        title('Unexpected reward')
        n = n+1;
        end
        if mworks.doBlock2
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_block2),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_block2),3),[],2)./sqrt(nIC),'r');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        ylim([-10 10]);
        title('CS-')
        end
        sgtitle([mouse ' ' tempClusName])
        savefig(fullfile(analysis_out,img_fn,[tempClusName '_cueAlign_dFoverF.fig']))

        figure; n=1;
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew),3),[],2)./sqrt(nIC)).*(1000./frameRateHz));
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([0 6]);
        title('CS+')
        n = n+1;
        if length(ind_omit>5)
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_omit),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'b');
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Omit')
        n = n+1;
        end
        if length(ind_unexp>5)
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_unexp),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_unexp),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'g');
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Unexpected reward')
        n = n+1;
        end
        if mworks.doBlock2
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_block2),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_block2),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'r');
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([0 6]);
        title('CS-')
        end
        sgtitle([mouse ' ' tempClusName])
        savefig(fullfile(analysis_out,img_fn, [tempClusName '_cueAlign_events_Hz.fig']))
        
        elseif indivFig
            
        figure; %n = 1;
        %subplot(s,1,n) 
        if expt.whichBlock(1,exp_subset) == 0
            shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_rew),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_rew),3),[],2)./sqrt(nIC),'red');
            vline(cue_rew_int)        
            xlabel('Time from cue')
            ylabel('Avg dF/F')
            title('CS-')
        elseif expt.whichBlock(1,exp_subset) == 1
            shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_rew),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_rew),3),[],2)./sqrt(nIC),'black');
            vline(cue_rew_int)        
            xlabel('Time from cue')
            ylabel('Avg dF/F')
            title('CS+')
        end
        sgtitle([mouse ' ' tempClusName])
        savefig(fullfile(analysis_out,img_fn,tempClusName, '_cueAlign_dFoverF.fig'))

        figure; %n=1;
        %subplot(s,1,n)
        if expt.whichBlock(1,exp_subset) == 0
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'red');
            vline(cue_rew_int)
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([0 6]);
            title('CS-')
        elseif expt.whichBlock(1,exp_subset) == 1
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'black');
            vline(cue_rew_int)
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([0 6]);
            title('CS+')
        end
        sgtitle([mouse ' ' tempClusName])
        savefig(fullfile(analysis_out,img_fn,tempClusName, '_cueAlign_events_Hz.fig'))
        end
        
        
        if length(ind_omit>10)
            figure;
            subplot(2,1,1)
            thresh = mean(diff(ind_omit));
            ind_omit_short = find(diff(ind_omit)<thresh)+1;
            ind_omit_long = find(diff(ind_omit)>thresh)+1;
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_omit(ind_omit_short)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit(ind_omit_short)),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'k');
            hold on
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_omit(ind_omit_long)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit(ind_omit_long)),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'b');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([0 6]);
            title(['Short (black n = ' num2str(length(ind_omit_short)) ') vs long (blue n = ' num2str(length(ind_omit_long)) ') interval between omits'])
            ind_rew_preomit = ind_omit-1;
            if find(ind_rew_preomit==0)
              ind_rew_preomit(find(ind_rew_preomit==0)) = [];
            end
            ind_rew_postomit = ind_omit+1;
            if find(ind_rew_postomit)
              ind_rew_postomit(find(ind_rew_postomit>nTrials)) = [];
            end
            subplot(2,1,2)
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew_preomit),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew_preomit),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'k');
            hold on
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew_postomit),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew_postomit),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'b');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Reward trial pre (black n = ' num2str(length(ind_rew_preomit)) ') vs post (blue n = ' num2str(length(ind_rew_postomit)) ') omit trials'])
            sgtitle([mouse ' ' tempClusName])
            savefig(fullfile(analysis_out,img_fn, [tempClusName '_omitByInterval.fig']))
        else
            ind_omit_short = [];
            ind_omit_long = [];
            ind_rew_preomit = [];
            ind_rew_postomit = [];
        end
        if length(ind_unexp>10)
            figure;
            thresh = mean(diff(ind_unexp));
            ind_unexp_short = find(diff(ind_unexp)<thresh)+1;
            ind_unexp_long = find(diff(ind_unexp)>thresh)+1;
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_unexp(ind_unexp_short)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_unexp(ind_unexp_short)),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'k');
            hold on
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_unexp(ind_unexp_long)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_unexp(ind_unexp_long)),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'b');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title([date ' ' mouse '- Short (black n = ' num2str(length(ind_unexp_short)) ') vs long (blue n = ' num2str(length(ind_unexp_long)) ') interval between unexpected reward'])
            savefig(fullfile(analysis_out,img_fn, [tempClusName '_unexpByInterval.fig']))
        else
            ind_unexp_short = [];
            ind_unexp_long = [];
        end
        
        if unique(expt(exp_subset).name == 'Inter')
         figure; n=1;
        %subplot(s,1,n)
        plot(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew),3),2).*(1000./frameRateHz),'black');
        hold on;
        plot(tt, nanmean(nanmean(targetAlign_events(:,:,ind_block2),3),2).*(1000./frameRateHz),'r');
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        %ylim([0 6]);
        title('Stacked spike rates separating CS+ and CS- trials');
        legend('CS+','CS-');
        sgtitle([mouse ' ' tempClusName])
        savefig(fullfile(analysis_out,img_fn, [tempClusName '_line_stackedSpike_events_Hz.fig']))
        
         figure; n=1;
        %subplot(s,1,n)
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'black');
        hold on;
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_block2),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_block2),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'r');
        vline(cue_rew_int)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        %ylim([0 6]);
        title([mouse ' ' tempClusName ' CS+/CS-'])
        %legend('CS+', 'sd CS+' ,'CS-', 'sd CS-');
        savefig(fullfile(analysis_out,img_fn, [tempClusName '_error_stackedSpike_events_Hz.fig']))
        
        end 
        
        if unique(expt(exp_subset).name == 'Inter')
        figure;
        n = floor(nTrials./3);
        start = 1;
        for i = 1:3
            subplot(3,2,start)
            ind_rew_temp = intersect(ind_rew,1+(i-1).*n:i*n);
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew_temp),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew_temp),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
            hold on
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([0 6]);
            title(['Trials ' num2str(1+(i-1).*n) ':' num2str(i*n)])
            vline(cue_rew_int)
            if length(ind_block2)>=10
                subplot(3,2,start+1)
                ind_block2_temp = intersect(ind_block2,1+(i-1).*n:i*n);
                shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_block2_temp),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_block2_temp),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'r');
                ylim([0 6])
                xlabel('Time from cue')
                ylabel('Spike rate (Hz)')
                title(['Trials ' num2str(1+(i-1).*n) ':' num2str(i*n)])
                vline(cue_rew_int)
            end
            start = start+2;
        end
        sgtitle([mouse ' ' tempClusName ' | CS+ (black), CS- (red)'])
        savefig(fullfile(analysis_out,img_fn, [tempClusName '_repsByTrial.fig']))
        
        elseif unique(expt(exp_subset).name == 'Block')
            if expt(exp_subset).whichBlock(1,exp_subset) == 0
        figure;
        n = floor(nTrials./3);
        start = 1;
        for i = 1:3
            subplot(3,2,start)
            ind_rew_temp = intersect(ind_rew,1+(i-1).*n:i*n);
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew_temp),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew_temp),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'r');
            hold on
            ylim([0 6])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Trials ' num2str(1+(i-1).*n) ':' num2str(i*n)])
            vline(cue_rew_int)
           
            start=start+2;
        end 
        sgtitle([mouse ' ' tempClusName ' | CS- (red)'])
        savefig(fullfile(analysis_out,img_fn, [tempClusName '_repsByTrial.fig']))
            elseif expt(exp_subset).whichBlock(1,exp_subset) == 1
        figure;
        n = floor(nTrials./3);
        start = 1;
        for i = 1:3
            subplot(3,2,start)
            ind_rew_temp = intersect(ind_rew,1+(i-1).*n:i*n);
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew_temp),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew_temp),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
            hold on
            ylim([0 6])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Trials ' num2str(1+(i-1).*n) ':' num2str(i*n)])
            vline(cue_rew_int)
        
            start=start+2;
        end
        sgtitle([mouse ' ' tempClusName ' | CS+ (black)'])
        savefig(fullfile(analysis_out,img_fn, [tempClusName '_repsByTrial.fig']))
            end
        end
        %{    
        if unique(expt.name == 'Crus')
            load(fullfile(analysis_out,img_fn, [img_fn '_splitImage.mat']))
            figure;
            indL = find(maskCat==1);
            indR = find(maskCat==2);
            subplot(2,2,1)
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,indL,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,indL,ind_rew),3),[],2).*(1000./frameRateHz))./sqrt(length(indL)),'k');
            ylim([0 5])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Reward- Left side- n=' num2str(length(indL))])
            vline(cue_rew_int)
            subplot(2,2,2)
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,indR,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,indR,ind_rew),3),[],2).*(1000./frameRateHz))./sqrt(length(indR)),'k');
            ylim([0 5])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            vline(cue_rew_int)
            title(['Reward- Right side- n=' num2str(length(indR))])
            subplot(2,2,3)
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,indL,ind_omit),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,indL,ind_omit),3),[],2).*(1000./frameRateHz))./sqrt(length(indL)),'r');
            ylim([0 5])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Omit- Left side- n=' num2str(length(indL))])
            vline(cue_rew_int)
            subplot(2,2,4)
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,indR,ind_omit),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,indR,ind_omit),3),[],2).*(1000./frameRateHz))./sqrt(length(indR)),'r');
            ylim([0 5])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['omit- Right side- n=' num2str(length(indR))])
            vline(cue_rew_int)
            sgtitle([date ' ' mouse '- Reward (black), Omit (red)'])
            savefig(fullfile(analysis_out,img_fn, [img_fn '_repsByCrus.fig']))
        end
        %}
        
        save(fullfile(analysis_out,img_fn, [tempClusName '_targetAlign.mat']), 'ind_rew', 'ind_block2', 'targetAlign_events', 'targetAligndFoverF', 'prewin_frames', 'postwin_frames', 'tt', 'frameRateHz')
        save(fullfile(analysis_out,img_fn, [tempClusName '_input.mat']), 'mworks')
end 

%%

clear; close all;
analysis_out = 'A:\home\carlo\analysis\2P\';
bdata_source = 'A:\home\carlo\rawData\behavioral\';

removing_bad_trials = false; %bad_trials = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]; %added to remove shitty trials, turn off if not needed

RCExptListCarlo_Inter
%RCExptListCarlo_Block
pairings = loadRCList_Carlo;
crpTrainDayList; %clearvars -except an* b* days
ii = length(days);
sessions=days{1,1};
%for  j = 1:size(expt.ttl,2)

    id = pairings{1,1};
    exp_subset = id;
    iexp = exp_subset;  %1:nexp
    mouse = strtrim(expt(iexp).mouse);
    date = expt(iexp).date;
    run = expt(iexp).run;
    fprintf([date ' ' mouse ' ' run '\n']);
    img_fn = [date '_img' mouse '\getTC_' run '\'];
    img_fn2 = [date '_img' mouse '\kMeans\'];
    load([analysis_out, img_fn2, 'filtered_denormKMeans.mat']);
    load([analysis_out, img_fn, date '_img' mouse '_' run 'saveOutputs.mat']);
    input = getBxData_Carlo(bdata_source, sessions);

for c = 1:size(allclusters,1)
    tempClusName=sprintf('cluster_%d',c);
    load(fullfile(analysis_out,img_fn, [tempClusName '_targetAlign.mat']));
    
    rewDelay_frames =  round((mean(celleqel2mat_padded(input.reactTimesMs))/1000).*frameRateHz);% there's 700ms between the cue and the reward delivery !!!!! if you change tooFastMs in the varibles in MWorks, this needs to be changed
    
    cTargetOn = input.cTargetOn;
    if iscell(cTargetOn) % if it is a cell, it means cTargetOn wasn't being over written in extract TC. If it's not a cell, it's already over written in extract TC. can be used directly
        cTargetOn = celleqel2mat_padded(input.cTargetOn);
        cTargetOn(1) = nan; % first trial doesn't have reward 
    end
 
    nTrials = size(cTargetOn(1:end),2);
    nIC = size(targetAlign_events,2); %you get targetalign_events from neural align, this is the neural data (in spikes) aligned to cue
    lickCueAlign =  nan(prewin_frames+postwin_frames,nTrials);
    lickCounterVals = cell(1,nTrials);
    nFrames = input.counterValues{nTrials}(end);
    lickDelay_frames =  round(0.1.*frameRateHz);
    lickSearch_frames =  round(0.3.*frameRateHz);
    
    postLick_frames = round(0.5.*frameRateHz);
    lickSearchRange = prewin_frames+lickDelay_frames:size(lickCueAlign,1)-lickSearch_frames;
    preRew_lickSearchRange_700ms = prewin_frames+lickDelay_frames:prewin_frames+lickDelay_frames+rewDelay_frames;
    lickBurstStart = nan(1,nTrials);
    postRew_lickBurstStart = nan(1,nTrials);
    lickBurstHz_all = nan(1,nTrials);
    preRew_lickBurstHz = nan(1,nTrials);
    postRew_lickBurstHz = nan(1,nTrials);
    postRew_lickAlignEvents = nan(2.*postLick_frames,nIC,nTrials);
    postRew_lickAlign = nan(2.*postLick_frames,nTrials);
    lastPreRew_lickAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    firstPostRew_lickAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    lastPreRew_lickAlign = nan(3.*postLick_frames,nTrials);
    firstPostRew_lickAlign = nan(3.*postLick_frames,nTrials);
    rewAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    preRewTrials = [];
    postRewTrials = [];
    postRew_lickSearchRange = prewin_frames+rewDelay_frames+lickDelay_frames:size(lickCueAlign,1)-lickSearch_frames-postLick_frames;% need to '-lickSearch_frames-postLick_frames' because later the inds needs to + lickSearch_frames or +postLick_frames, this is just for the following inds to be within matrix dimentions
    lastPreRewLickFrame = nan(1,nTrials);
    firstPostRewLickFrame = nan(1,nTrials);
    
    tTInx = 1;
    for itrial = 1:nTrials
        if ~isnan(cTargetOn(itrial))
            if cTargetOn(itrial)+postwin_frames-1 < nFrames
                lickTimes = input.lickometerTimesUs{itrial}; 
                counterTimes = input.counterTimesUs{itrial};
                counterVals = input.counterValues{itrial};
                lickCounterVals{itrial} = zeros(size(lickTimes));
                lickTC{itrial} = zeros(size(counterVals));% find how many licks for each frame
             %%%this for loop pulls out every lickTime that falls between the start of two frames
               %and if licks are found it outputs the frame licks occur to lickCounterVals
                for icount = 1:length(counterTimes)-1
                    ind = find(lickTimes>counterTimes(icount) & lickTimes<counterTimes(icount+1));
                    if ~isempty(ind)
                        lickCounterVals{itrial}(1,ind) = icount; % find which counter in this trial has licks
                    end
                end
             %%%this for loop counts down the number of frames in a trial and outputs a logical
               %outlining whether or not a lick occured between two counters (or frames)
                for ival = 1:length(counterTimes)
                    ind = find(lickCounterVals{itrial} == ival);
                    if ~isempty(ind)
                        lickTC{itrial}(1,ival) = length(ind);% if the mice doesn't lick during that counter than it's zero, if it does, it's one. length(ind) should always be 1.
                    end
                end
             %%%IF the last frame of the trial (from the start of the ITI preceeding the trial to the start of the ITI following the trial)
               %minus the frame where the cue is shown IS greater than the frame difference between   
                if input.counterValues{itrial}(end)-cTargetOn(itrial) > postwin_frames
                    %lickcuealign aligns licking of each frame to cue
                    lickCueAlign(:,itrial) = lickTC{itrial}(1,cTargetOn(itrial)-prewin_frames-counterVals(1):cTargetOn(itrial)+postwin_frames-1-counterVals(1));
                end
                ind = intersect(lickSearchRange,find(lickCueAlign(:,itrial)));
             %%%this for loop find frames where more than 3 licks occurred within a fixed window, and index that frame as the start of burst lick
                for i = 1:length(ind)
                    ilick = ind(i);
                    if sum(lickCueAlign(ilick:ilick+lickSearch_frames-1,itrial),1) >= 3 % more than 3 bursts within about 300ms
                        lickBurstStart(:,itrial) = ilick;% when licking burst happens
                        break
                    end
                end
                ind_pre = intersect(preRew_lickSearchRange_700ms,find(lickCueAlign(:,itrial))); %finds every instance of a lick that occurs within the search window
                ind_post = intersect(postRew_lickSearchRange,find(lickCueAlign(:,itrial))); %finds every instance of a lick that occurs [after the reward delivery]
                ind_all = intersect(lickSearchRange,find(lickCueAlign(:,itrial)));
                preRew_lickBurstHz(:,itrial) = length(ind_pre)./(mean(celleqel2mat_padded(input.reactTimesMs))/1000);
                postRew_lickBurstHz(:,itrial) = length(ind_post)./(length(postRew_lickSearchRange)./frameRateHz);
                lickBurstHz_all(:,itrial) = length(ind_all)./(length(lickSearchRange)./frameRateHz);
                ind = intersect(postRew_lickSearchRange,find(lickCueAlign(:,itrial)));
             %%%similar idea as the above for loop, but instead aligns neural data as POST-REWARD lick events 
                for i = 1:length(ind)
                    ilick = ind(i);
                    if sum(lickCueAlign(ilick:ilick+lickSearch_frames-1,itrial),1) >= 3
                        postRew_lickBurstStart(:,itrial) = ilick; %array of every lick within post-reward window if 3+ licks were recorded for that trial
                        postRew_lickAlignEvents(:,:,itrial) = targetAlign_events(ilick-postLick_frames:ilick+postLick_frames-1,:,itrial);% align neural data to first lick after reward
                        postRew_lickAlign(:,itrial) = lickCueAlign(ilick-postLick_frames:ilick+postLick_frames-1,itrial);% align licking data to first lick
                        break
                    end
                end
                ind_pre = find(lickCueAlign(prewin_frames:prewin_frames+rewDelay_frames,itrial),1,'last');
             %%%if there are licks recorded between cue and reward delivery (767 ms), align neural data to those licks 
                if ~isempty(ind_pre)
                    lastPreRewLickFrame(1,itrial) = ind_pre;
                    preRewTrials = [preRewTrials itrial];
                    % align neural and licking data to the last lick
                    lastPreRew_lickAlignEvents(:,:,itrial) = targetAlign_events(ind_pre-postLick_frames+prewin_frames:ind_pre+postLick_frames+postLick_frames-1+prewin_frames,:,itrial);
                    lastPreRew_lickAlign(:,itrial) = lickCueAlign(ind_pre-postLick_frames+prewin_frames:ind_pre+postLick_frames+postLick_frames-1+prewin_frames,itrial);
                end
                ind_post = find(lickCueAlign(prewin_frames+rewDelay_frames:prewin_frames+postwin_frames-postLick_frames-postLick_frames,itrial),1,'first');
             %%%   
                if ~isempty(ind_post)
                    firstPostRewLickFrame(1,itrial) = ind_post;
                    postRewTrials = [postRewTrials itrial];
                    % align neural and licking data to the first lick
                    firstPostRew_lickAlignEvents(:,:,itrial) = targetAlign_events(ind_post-postLick_frames+prewin_frames+rewDelay_frames:ind_post+postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,:,itrial);
                    firstPostRew_lickAlign(:,itrial) = lickCueAlign(ind_post-postLick_frames+prewin_frames+rewDelay_frames:ind_post+postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,itrial);
                end
                rewAlignEvents(:,:,itrial) =targetAlign_events(-postLick_frames+prewin_frames+rewDelay_frames:postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,:,itrial);
            
            tTInx = tTInx+1;
            end
        end
    end
    
    figure; %align licking data to cue
    shadedErrorBar(tt, nanmean(lickCueAlign,2).*(1000./frameRateHz), (nanstd(lickCueAlign,[],2)./sqrt(unique(sum(~isnan(lickCueAlign),2))).*(1000./frameRateHz)));
    hold on;
    scatter((lickBurstStart-prewin_frames).*(1000./frameRateHz), ones(size(lickBurstStart)).*(nanmean(nanmean(lickCueAlign,2).*(1000./frameRateHz),1)+5), 'mx'); %10.*ones(size(lickBurstStart)) replaced with nanmean(lickCueAlign,2)+5
    xlabel('Time from cue');
    ylabel('Lick rate (Hz)');
    ylim([0 10])
    vline(0,'k');
    vline(767,'r');
    sgtitle([tempClusName '- licks aligned to cue']);
    savefig(fullfile(analysis_out,img_fn, [tempClusName '_cueAlign_lickHz.fig']));
    
    nIC = size(targetAlign_events,2);
    if sum(~isnan(lickBurstStart))>6
        [sortlick sortlick_ind] = sort(lickBurstStart,'ascend');
        nburst = sum(~isnan(lickBurstStart));
        nnan = sum(isnan(lickBurstStart));
        ind_early_bst = sortlick_ind(1:floor(nburst/4));
        ind_late_bst = sortlick_ind(nburst-floor(nburst/4)+1:end-nnan);
        early_bst_time = mean((lickBurstStart(:,ind_early_bst)-prewin_frames).*(1000./frameRateHz),2);
        late_bst_time = mean((lickBurstStart(:,ind_late_bst)-prewin_frames).*(1000./frameRateHz),2);
    else
        ind_early_bst = [];
        ind_late_bst = [];
        early_bst_time = [];
        late_bst_time = [];
    end
   
    if sum(~isnan(lickBurstStart))>6
        figure; %plot neural data of trials of early vs. late bursts
        shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_early_bst),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_early_bst),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'm');
        hold on;
        shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_late_bst),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_late_bst),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
        scatter((lickBurstStart(:,ind_early_bst)-prewin_frames).*(1000./frameRateHz),(((nanmean(nanmean(nanmean(targetAlign_events,3),2),1).*(1000./frameRateHz))+3).*ones(1,length(ind_early_bst))),'mx');
        scatter((lickBurstStart(:,ind_late_bst)-prewin_frames).*(1000./frameRateHz),(((nanmean(nanmean(nanmean(targetAlign_events,3),2),1).*(1000./frameRateHz))+3).*ones(1,length(ind_late_bst))),'kx');
        xlabel('Time from cue');
        ylabel('Spike rate (Hz)');
        ylim([0 5]);
        vline(0,'k');
        vline(767,'r');
        title([(tempClusName) 'early lick bursts (n = ' num2str(length(ind_early_bst)) ')'; (tempClusName) 'later lick bursts (n = ' num2str(length(ind_late_bst)) ')']);
        savefig(fullfile(analysis_out,img_fn, [tempClusName '_cueAlignSpiking_byLickTime.fig']));
        hold off;
    end
    
    pct_precue_burst = length(find((lickBurstStart-prewin_frames).*(1000./frameRateHz)<600))./size(lickBurstStart,2);
    tl = (1-postLick_frames:postLick_frames).*(1000./frameRateHz);
    
    figure;
    subplot(2,1,1); % align neural activity to lick onset, only burst trials are included 
    shadedErrorBar(tl, nanmean(nanmean(postRew_lickAlignEvents,3),2).*(1000./frameRateHz), (nanstd(nanmean(postRew_lickAlignEvents,3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    hold on;
    xlabel('Time from lick');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    title([num2str(sum(~isnan(squeeze(postRew_lickAlignEvents(1,1,:))))) ' trials post-reward lick bursts (black)']);
    subplot(2,1,2); % seperate early burst trials vs. late burst trials
    [sortlick sortlick_ind] = sort(postRew_lickBurstStart,'ascend'); %sorts trial-based instances of lick and outputs [sorted instances, index of original position before sorting]
    nburst = sum(~isnan(postRew_lickBurstStart)); %total trials with lick burst
    nnan = sum(isnan(postRew_lickBurstStart)); %total trials without lick burst (some may be CS-; some may have no burst but still be CS+)
    ind_prerew_early_bst = sortlick_ind(1:floor(nburst/4)); %first 1/4 chosen as early trials
    ind_prerew_late_bst = sortlick_ind(nburst-floor(nburst/4)+1:end-nnan); %last 1/4 chosen as late trials 
    early_bst_time = nanmean((postRew_lickBurstStart(:,ind_prerew_early_bst)-prewin_frames-rewDelay_frames).*(1000./frameRateHz),2);
    late_bst_time = nanmean((postRew_lickBurstStart(:,ind_prerew_late_bst)-prewin_frames-rewDelay_frames).*(1000./frameRateHz),2);
    shadedErrorBar(tl, nanmean(nanmean(postRew_lickAlignEvents(:,:,ind_prerew_early_bst),3),2).*(1000./frameRateHz), (nanstd(nanmean(postRew_lickAlignEvents(:,:,ind_prerew_early_bst),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    hold on;
    shadedErrorBar(tl, nanmean(nanmean(postRew_lickAlignEvents(:,:,ind_prerew_late_bst),3),2).*(1000./frameRateHz), (nanstd(nanmean(postRew_lickAlignEvents(:,:,ind_prerew_late_bst),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'b');
    xlabel('Time from lick');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    title([num2str(floor(nburst/4)) ' earliest burst lick trials [blk]: avg = ' num2str(chop(early_bst_time,3)) ' ms';  num2str(floor(nburst/4)) ' latest burst lick trials [blu]: avg = ' num2str(chop(late_bst_time,3)) ' ms)'])
    sgtitle(tempClusName);
    savefig(fullfile(analysis_out,img_fn, [tempClusName '_postRew_lickBurstAlignSpikingEvents.fig']))

    figure;
    subplot(2,1,1); % align licking data to first lick
    shadedErrorBar(tl, nanmean(postRew_lickAlign,2).*(1000./frameRateHz), (nanstd(postRew_lickAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(postRew_lickAlign),2))),'k');
    hold on;
    xlabel('Time from lick');
    ylabel('Lick rate (Hz)');
    ylim([0 inf]);
    title('Neural data aligned to first lick post-reward per trial')
    subplot(2,1,2); % plot early burst trials and late burst trials separately
    shadedErrorBar(tl, nanmean(postRew_lickAlign(:,ind_prerew_early_bst),2).*(1000./frameRateHz), (nanstd(postRew_lickAlign(:,ind_prerew_early_bst),[],2).*(1000./frameRateHz))./sqrt(length(ind_prerew_early_bst)),'k');
    hold on;
    shadedErrorBar(tl, nanmean(postRew_lickAlign(:,ind_prerew_late_bst),2).*(1000./frameRateHz), (nanstd(postRew_lickAlign(:,ind_prerew_late_bst),[],2).*(1000./frameRateHz))./sqrt(length(ind_prerew_late_bst)),'b');
    xlabel('Time from lick');
    ylabel('Lick rate (Hz)');
    ylim([0 inf]);
    title('Early lick burst trials [blk] and late lick burst trials [blu]')
    sgtitle([tempClusName '- post reward lick burst aligned spiking']);
    savefig(fullfile(analysis_out,img_fn, [tempClusName '_postRew_lickBurstAlignSpiking.fig']));
    hold off;
    
    figure;
    colour = {'k', 'b', 'r', 'm'};
    for i = 1:4
        plot(tt,(cumsum(nansum(lickCueAlign(:,1+((i-1)*floor(nTrials/4)):floor(nTrials/4)+((i-1)*floor(nTrials/4))),2))),colour{1,i}, 'LineWidth',2.0);
        hold on;
    end
    title([tempClusName '- cumulative licking by quarter session']);
    xlabel('Time from cue');
    ylabel('Cumulative Licks');
    legend({'first quar', 'second quar', 'third quar', 'fourth quar'},'Location','northwest');
    hold off;
    savefig(fullfile(analysis_out,img_fn, [tempClusName '_cumulativeLicking.fig']));
    
    tl_rew = (1-postLick_frames:(postLick_frames*2)).*(1000./frameRateHz);
    figure; % align neural and licking data to the first and last lick, respectively
    subplot(2,2,1);
    shadedErrorBar(tl_rew, nanmean(nanmean(lastPreRew_lickAlignEvents,3),2).*(1000./frameRateHz), (nanstd(nanmean(lastPreRew_lickAlignEvents,3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    hold on;
    title('Aligned to last lick before reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    subplot(2,2,3);
    shadedErrorBar(tl_rew, nanmean(lastPreRew_lickAlign,2).*(1000./frameRateHz), (nanstd(lastPreRew_lickAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(lastPreRew_lickAlign)))),'k');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 30]);
    title([num2str(length(preRewTrials)) ' pre-reward lick frequency']);
    subplot(2,2,2);
    shadedErrorBar(tl_rew, nanmean(nanmean(firstPostRew_lickAlignEvents,3),2).*(1000./frameRateHz), (nanstd(nanmean(firstPostRew_lickAlignEvents,3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    title('Aligned to first lick after reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    subplot(2,2,4);
    shadedErrorBar(tl_rew, nanmean(firstPostRew_lickAlign,2).*(1000./frameRateHz), (nanstd(firstPostRew_lickAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(firstPostRew_lickAlign)))),'k');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 inf]);
    title([num2str(length(postRewTrials)) ' post-reward lick frequency']);
    sgtitle([tempClusName '- Licks relative to reward']);
    hold off;
    savefig(fullfile(analysis_out,img_fn, [tempClusName '_lastVsFirstLick.fig']));
    
    [sortlickHz sortlickHz_ind] = sort(lickBurstHz_all,'ascend');
    nburst = sum(~isnan(lickBurstHz_all));
    nnan = sum(isnan(lickBurstHz_all));
    ind_low_bst = sortlickHz_ind(1:floor(nburst/4));
    ind_high_bst = sortlickHz_ind(nburst-floor(nburst/4)+1:end-nnan);
    HL_lickrate.low_rew = mean(lickBurstHz_all(:,ind_low_bst),2);
    HL_lickrate.high_rew = mean(lickBurstHz_all(:,ind_high_bst),2);
    
    [sortlickHz sortlickHz_ind] = sort(preRew_lickBurstHz,'ascend');
    nnan = sum(isnan(preRew_lickBurstHz));
    nburst = sum(~isnan(preRew_lickBurstHz));
    ind_low_prerew = sortlickHz_ind(1:floor(nburst/4));
    ind_high_prerew = sortlickHz_ind(nburst-floor(nburst/4)+1:end-nnan);
    HL_lickrate.low_prerew = mean(preRew_lickBurstHz(:,ind_low_prerew),2);
    HL_lickrate.high_prerew = mean(preRew_lickBurstHz(:,ind_high_prerew),2);
    
    [sortlickHz sortlickHz_ind] = sort(postRew_lickBurstHz,'ascend');
    nburst = sum(~isnan(postRew_lickBurstHz));
    nnan = sum(isnan(postRew_lickBurstHz));
    ind_low_postrew = sortlickHz_ind(1:floor(nburst/4));
    ind_high_postrew = sortlickHz_ind(nburst-floor(nburst/4)+1:end-nnan);
    HL_lickrate.low_postrew = mean(postRew_lickBurstHz(:,ind_low_postrew),2);
    HL_lickrate.high_postrew = mean(postRew_lickBurstHz(:,ind_high_postrew),2);
    
    figure; % still plotting neural data, but seperate the trials based on licking rate of that trial
    subplot(3,1,1);
    shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_low_bst),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_low_bst),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
    hold on;
    shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_high_bst),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_high_bst),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([-1 2]);
    title(['Whole trial: ' num2str(chop(HL_lickrate.low_rew,2)) ' vs ' num2str(chop(HL_lickrate.high_rew,2)) ' Hz']);
    subplot(3,1,2);
    shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_low_prerew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_low_prerew),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');hold on;
    shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_high_prerew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_high_prerew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([-1 2]);
    title(['Pre-reward: ' num2str(chop(HL_lickrate.low_prerew,2)) ' vs ' num2str(chop(HL_lickrate.high_prerew,2)) ' Hz']);
    subplot(3,1,3);
    shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_low_postrew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_low_postrew),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');hold on;
    shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_high_postrew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_high_postrew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([-1 2]);
    title(['Post-reward: ' num2str(chop(HL_lickrate.low_postrew,2)) ' vs ' num2str(chop(HL_lickrate.high_postrew,2)) ' Hz']);
    sgtitle([tempClusName '- lick bursts by rate: low (blue) & high (black)']);
    hold off;
    savefig(fullfile(analysis_out,img_fn, [tempClusName '_cueAlignSpiking_byLickRate.fig']));
    
    save(fullfile(analysis_out,img_fn, [tempClusName '_cueAlignLick.mat']), 'firstPostRewLickFrame', ...
        'lastPreRewLickFrame', 'tl_rew', 'firstPostRew_lickAlignEvents', 'lastPreRew_lickAlignEvents', ...
        'lickCueAlign', 'lickBurstStart', 'lickCounterVals', 'lickSearch_frames', 'lickDelay_frames',...
        'lickTC', 'postwin_frames', 'prewin_frames', 'frameRateHz', 'tt', 'ind_early_bst', 'ind_late_bst', ...
        'early_bst_time', 'late_bst_time', 'pct_precue_burst', 'postRew_lickAlignEvents', 'postLick_frames', ...
        'postRew_lickBurstStart','tl', 'ind_prerew_early_bst', 'ind_prerew_late_bst','postRew_lickAlign',...
        'preRew_lickBurstHz','postRew_lickBurstHz','ind_low_prerew','ind_high_prerew','ind_low_postrew',...
        'ind_high_postrew','ind_high_bst','ind_low_bst','HL_lickrate');
end

%%

   
    %% stacked TC for each dendrite in a cluster
    colours = ['r','m','b','g','c','y'];%{[1 0.2 0.2],[0.5 0.184 0.556],[0 0.5 1],[0 0.6 0.3],[0.2 0.8 0.8],[1 0.7 0.1250]}; % r m b g c y
    for ii = 1:size(allclusters,1)
        shift=0; fig=figure; 
        if size(allclusters{ii},2)>10
        randDend=randperm(size(allclusters{ii},2),10);
        else
        randDend=1:size(allclusters{ii},2);
        end
        tempTC = TCave_cl(:,allclusters{ii,1}(1,randDend));
        for iii = 1:size(tempTC,2); v=1;
                for iv = 1:size(tempTC,1)
                if spk.spk_logic_cl(iv,allclusters{ii,1}(1,randDend(1,iii)))==1
                    spkTags(v,1) = iv; 
                    v=v+1;
                end
                end        
        plot(tempTC(:,iii)+shift,'k');
        hold on;
        plot(spkTags,nanmean(tempTC(:,iii)+shift+20000),'k.');%'color',[colours(ii) '*']);
        shift = shift+50000;  %clearvars spkTags;
        end
    set(gca,'YTick',50000:50000:shift);
    set(gca,'YTicklabel',(allclusters{ii,1}(1,randDend)));    
    end

   
     %}