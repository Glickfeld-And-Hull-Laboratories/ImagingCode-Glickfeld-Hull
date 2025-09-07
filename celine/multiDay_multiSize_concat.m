
clear all; close all; clc

ds = 'DART_V1_YM90K_Celine';
dataStructLabels = {'contrastxori'};
experimentFolder = 'SST_YM90K';
rc = behavConstsDART;
eval(ds);

sess_list = [8 10 20 22];
nSess = length(sess_list);
nd = 2;

% Session reference selection
x = input('Which session was used as reference for matching: 0- baseline, 1- post-DART: ');
if x == 0
    pre = 1; post = 2;
    disp('baseline used as reference')
else
    pre = 2; post = 1;
    disp('post-DART used as reference')
end

targetCon = [.25 .5 1];
frame_rate = 15;
nCon = length(targetCon);
nSize = 3;

% Create output directory
sess_title = num2str(sess_list(1));
for iSess = 2:nSess
    sess_title = [sess_title '_' num2str(sess_list(iSess))];
end
d = char(string(datetime('today')));

if nSess == 1
    if expt(sess_list(1)).multiday_timesincedrug_hours > 0
        dart_str = [expt(sess_list(1)).drug '_' num2str(expt(sess_list(1)).multiday_timesincedrug_hours) 'Hr'];
    else
        dart_str = 'control';
    end
    fnout = fullfile(rc.achAnalysis, experimentFolder, expt(sess_list(1)).mouse, ['multiday_' dart_str], d);
else
    fnout = fullfile(rc.achAnalysis, experimentFolder, ['concat' sess_title], d);
end
mkdir(fnout); cd(fnout)

%% Concatenating data
% Initialize concatenation variables
mice = []; red_concat = []; green_concat = []; nKeep_concat = [];
dirs_concat = []; cons_concat = []; 
% Initialize cell arrays for concatenation
cell_vars = {'tc_trial_avrg_stat', 'tc_trial_avrg_loc', 'pref_responses_stat', 'pref_responses_loc', ...
'h', 'data_resp', 'RIx', 'norm_dir_resp', 'pref_dir', 'noiseCorr', 'sigCorr', ...
'nonPref_trial_avrg_stat', 'nonPref_trial_avrg_loc', 'data_dfof_runOnset'};
for i = 1:length(cell_vars)
 eval([cell_vars{i} '_concat = cell(1,nd);']);
end
pref_allTrials_stat_concat = cell(nCon, nSize, nd);
pref_allTrials_loc_concat = cell(nCon, nSize, nd);
drug = cell(1, nSess);
% Main concatenation loop
for iSess = 1:nSess
     thisSess = sess_list(iSess);
     mouse = expt(thisSess).mouse;
     mice = [mice; mouse];
     drug{iSess} = expt(thisSess).drug;
    % Determine data path
    if expt(thisSess).multiday_timesincedrug_hours > 0
     dart_str = [expt(thisSess).drug '_' num2str(expt(thisSess).multiday_timesincedrug_hours) 'Hr'];
    else
     dart_str = 'control';
    end
     fn_multi = fullfile(rc.achAnalysis, experimentFolder, mouse, ['multiday_' dart_str]);
    % Load data
     load(fullfile(fn_multi, 'tc_keep.mat'));
     load(fullfile(fn_multi, 'resp_keep.mat'));
     load(fullfile(fn_multi, 'input.mat'));
     load(fullfile(fn_multi, 'behavioral_state.mat'));
     load(fullfile(fn_multi, 'cell_analysis.mat'));
     load(fullfile(fn_multi, 'HT_pyr_relationship.mat'));
     nKeep = size(tc_trial_avrg_stat{post}, 2);
    % Process trial conditions
     tCon_match = cell(1, nd);
     nTrials = [];
    for id = 1:nd
     nTrials = [nTrials, length(input(id).tBaseGratingContrast)];
     tCon_match{id} = celleqel2mat_padded(input(id).tGratingContrast(1:nTrials(id)));
    end
     dirs = unique(celleqel2mat_padded(input(post).tGratingDirectionDeg(1:nTrials(post))));
     cons = unique(tCon_match{post});
     sharedCon = find(ismember(cons, targetCon));
     sizes=unique(cell2mat(input(id).tGratingDiameterDeg));
    % Concatenate basic variables
     dirs_concat = [dirs_concat, dirs];
     cons_concat = [cons_concat, cons(sharedCon)];
     red_concat = [red_concat, red_cells_keep];
     green_concat = [green_concat, green_cells_keep];
     nKeep_concat = [nKeep_concat, nKeep];
    % Concatenate day-specific data
    for id = 1:nd
     tc_trial_avrg_stat_concat{id} = cat(2, tc_trial_avrg_stat_concat{id}, tc_trial_avrg_stat{id}(:,:,:,:));
     tc_trial_avrg_loc_concat{id} = cat(2, tc_trial_avrg_loc_concat{id}, tc_trial_avrg_loc{id}(:,:,sharedCon,:));
     pref_responses_loc_concat{id} = cat(1, pref_responses_loc_concat{id}, pref_responses_loc{id}(:,sharedCon,:));
     pref_responses_stat_concat{id} = cat(1, pref_responses_stat_concat{id}, pref_responses_stat{id}(:,sharedCon,:));
     RIx_concat{id} = cat(1, RIx_concat{id}, sum(RIx{id}));
     norm_dir_resp_concat{id} = cat(1, norm_dir_resp_concat{id}, norm_dir_resp{id});
     pref_dir_concat{id} = cat(2, pref_dir_concat{id}, prefDir_keep{id});
     noiseCorr_concat{id} = cat(2, noiseCorr_concat{id}, noiseCorr{id});
     sigCorr_concat{id} = cat(2, sigCorr_concat{id}, sigCorr{id});
     h_concat{id} = cat(1, h_concat{id}, h_keep{id});
     data_resp_concat{id} = cat(1, data_resp_concat{id}, data_resp_keep{id});
    end
    if iSess==1
        norm_diff_concat=norm_diff;
    else
        norm_diff_concat=cat(4,norm_diff_concat, norm_diff);
    end
     fprintf('Session %d of %d completed\n', iSess, nSess);
end
% Final variables
red_ind_concat = find(red_concat);
green_ind_concat = find(green_concat);
cons = targetCon;
nKeep_total = sum(nKeep_concat);

% Clear individual session data, keep only concatenated
clear tc_trial_avrg_stat tc_trial_avrg_loc pref_responses_stat pref_responses_loc
clear h data_resp RIx norm_dir_resp pref_dir noiseCorr sigCorr
clear nonPref_trial_avrg_stat nonPref_trial_avrg_loc data_dfof_runOnset norm_diff
clear red_cells_keep green_cells_keep prefDir_keep h_keep data_resp_keep
%%  Cell selection - find cells with running and stationary data for both days
haveRunning = cell(1,nd);
haveStat = cell(1,nd);
for id = 1:nd
    haveRunning{id} = sum(squeeze(sum(~isnan(pref_responses_loc_concat{id}),2)),2) == nSize*nCon;
    haveStat{id} = sum(squeeze(sum(~isnan(pref_responses_stat_concat{id}),2)),2) == nSize*nCon;
end

haveRunning_both = find(haveRunning{pre} .* haveRunning{post});
haveStat_both = find(haveStat{pre} .* haveStat{post});
runningCells = intersect(haveStat_both, haveRunning_both);

% Find cells responsive to small and large sizes
respToSmall = logical(sum(squeeze(sum(h_concat{pre}(:,:,:,1),2)),2) + ...
                     sum(squeeze(sum(h_concat{post}(:,:,:,1),2)),2));
respToLarge = logical(sum(squeeze(sum(h_concat{pre}(:,:,:,nSize),2)),2) + ...
                     sum(squeeze(sum(h_concat{post}(:,:,:,nSize),2)),2));


% Cell type assignments
runningGreen = intersect(runningCells, green_ind_concat);
runningRed = intersect(runningCells, red_ind_concat);
statGreen = green_ind_concat;
statRed = red_ind_concat;

% Mouse indices
mouseInds = cell(1,nSess);
start = 1;
for iMouse = 1:nSess
    mouseInds{iMouse} = start:(start-1)+nKeep_concat(iMouse);
    start = start + nKeep_concat(iMouse);
end

% Cell counts for each mouse
cellCounts = nan(nSess,2);
cellCountsGreen = nan(nSess,2);
mouseNames = [];

for iMouse = 1:nSess
    cellCounts(iMouse,1) = length(intersect(runningRed', mouseInds{iMouse}));
    cellCounts(iMouse,2) = length(intersect(statRed', mouseInds{iMouse}));
    cellCountsGreen(iMouse,1) = length(intersect(runningGreen', mouseInds{iMouse}));
    cellCountsGreen(iMouse,2) = length(intersect(statGreen', mouseInds{iMouse}));
    mouseNames = [mouseNames, string(mice(iMouse,:))];
end

% Save cell count tables
cellCountTableRed = table(cellCounts, RowNames=mouseNames)
cellCountTableGreen = table(cellCountsGreen, RowNames=mouseNames)
writetable(cellCountTableRed, fullfile(fnout,'cellCounts.csv'), 'WriteRowNames', true);

% Running by condition matrix
runningByCondition = nan(nKeep_total, nCon, nSize);
for iCon = 1:nCon
    for iSize = 1:nSize
        runningPre = ~isnan(pref_responses_loc_concat{pre}(:,iCon,iSize));
        runningPost = ~isnan(pref_responses_loc_concat{post}(:,iCon,iSize));
        runningByCondition(:,iCon,iSize) = runningPre .* runningPost;
    end
end


%% plot stationary timecourses for all cells
plotNeuralTimecourse(tc_trial_avrg_stat_concat, tc_trial_avrg_stat_concat, ...
    red_ind_concat, green_ind_concat, ...
    'UseDashedLines', [false, true], ...
    'Colors1', {'k', 'b'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'b'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'SST', 'Pyr'}, ...
    'StimStart', 31);

figs = findobj('Type', 'figure');
for i = 1:length(figs)
    figure(figs(i));
    saveas(gcf, sprintf('neural_timecourse_size_%d.pdf', i));
end


%% plot stationary timecourses for cells matched across behavioral state within each stim condition

% make figure with se shaded, one figure per contrast - stationary

tc_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_stat = cell(1,nd); %same for red
tc_green_se_stat = cell(1,nd); %this will be the se across all green cells
tc_red_se_stat = cell(1,nd); %same for red


nGreen = nan(nCon,nSize);
nRed = nan(nCon,nSize);

for id = 1:nd
  for iCon = 1:nCon
      for iSize = 1:nSize
        theseGreenCells = intersect(green_ind_concat, find(runningByCondition(:,iCon,iSize)));
        theseRedCells = intersect(red_ind_concat, find(runningByCondition(:,iCon,iSize)));
        nGreen(iCon, iSize)=length(theseGreenCells);
        nRed(iCon, iSize)=length(theseRedCells);

        tc_green_avrg_stat{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_stat_concat{id}(:,theseGreenCells,iCon,iSize),2);
        green_std=nanstd(tc_trial_avrg_stat_concat{id}(:,theseGreenCells,iCon,iSize),[],2);
        tc_green_se_stat{id}(:,iCon,iSize)=green_std/sqrt(length(theseGreenCells));
        
        tc_red_avrg_stat{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_stat_concat{id}(:,theseRedCells,iCon,iSize),2);
        red_std=nanstd(tc_trial_avrg_stat_concat{id}(:,theseRedCells,iCon,iSize),[],2);
        tc_red_se_stat{id}(:,iCon,iSize)=red_std/sqrt(length(theseRedCells));
        
        
        clear green_std red_std
      end 
    end
end

z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_stat{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
    for iSize = 1:nSize
    figure
    subplot(1,2,1) %for the first day
    
    
    
    ylim([-.02 .2]);
    hold on
    shadedErrorBar(t,tc_green_avrg_stat{pre}(:,iCon,iSize),tc_green_se_stat{pre}(:,iCon,iSize),'--k');
    hold on
    shadedErrorBar(t,tc_green_avrg_stat{post}(:,iCon,iSize),tc_green_se_stat{post}(:,iCon,iSize),'--b','transparent');
    hold on
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    title(['-HTP',' n = ', num2str(nGreen(iCon,iSize))])
    
    ylabel('dF/F') 
    xlabel('s') 
    set(gca,'XColor', 'none','YColor','none')
    
    
    subplot(1,2,2) %+HTP
    shadedErrorBar(t,tc_red_avrg_stat{pre}(:,iCon,iSize),tc_red_se_stat{pre}(:,iCon,iSize),'k');
    hold on
    shadedErrorBar(t,tc_red_avrg_stat{post}(:,iCon,iSize),tc_red_se_stat{post}(:,iCon,iSize),'b');
    ylim([-.02 .2]);
    hold on
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    ylabel('dF/F') 
    xlabel('s') 
    title(['+HTP',' n = ', num2str(nRed(iCon,iSize))])
    
    x0=5;
    y0=5;
    width=4;
    height=3;
    set(gcf,'units','inches','position',[x0,y0,width,height])
    set(gca,'XColor', 'none','YColor','none')
    
    sgtitle(['stationary, con ' num2str(cons(iCon)) ' size ' num2str(sizes(iSize))])
    
    print(fullfile(fnout,[num2str(cons(iCon)) '_' num2str(sizes(iSize)) 'matched_stat_cellType_timecourses.pdf']),'-dpdf');
    end
end 

% plots for running trials with cell matched across behavioral state

tc_green_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_loc = cell(1,nd); %same for red
tc_green_se_loc = cell(1,nd); %this will be the se across all green cells
tc_red_se_loc = cell(1,nd); %same for red


for id = 1:nd
  for iCon = 1:nCon
      for iSize = 1:nSize
        theseGreenCells = intersect(green_ind_concat, find(runningByCondition(:,iCon,iSize)));
        theseRedCells = intersect(red_ind_concat, find(runningByCondition(:,iCon,iSize)));

        tc_green_avrg_loc{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_loc_concat{id}(:,theseGreenCells,iCon,iSize),2);
        green_std=nanstd(tc_trial_avrg_loc_concat{id}(:,theseGreenCells,iCon,iSize),[],2);
        tc_green_se_loc{id}(:,iCon,iSize)=green_std/sqrt(length(theseGreenCells));
        
        tc_red_avrg_loc{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_loc_concat{id}(:,theseRedCells,iCon,iSize),2);
        red_std=nanstd(tc_trial_avrg_loc_concat{id}(:,theseRedCells,iCon,iSize),[],2);
        tc_red_se_loc{id}(:,iCon,iSize)=red_std/sqrt(length(theseRedCells));
        
        clear green_std red_std
      end 
  end
end

z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_loc{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
    for iSize = 1:nSize
    figure
    subplot(1,2,1) %for the first day
    
    
    
    ylim([-.02 .2]);
    hold on
    shadedErrorBar(t,tc_green_avrg_loc{pre}(:,iCon,iSize),tc_green_se_loc{pre}(:,iCon,iSize),'--k');
    hold on
    shadedErrorBar(t,tc_green_avrg_loc{post}(:,iCon,iSize),tc_green_se_loc{post}(:,iCon,iSize),'--b','transparent');
    hold on
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    title(['-HTP',' n = ', num2str(nGreen(iCon,iSize))])
    
    ylabel('dF/F') 
    xlabel('s') 
    set(gca,'XColor', 'none','YColor','none')
    
    
    subplot(1,2,2) %+HTP
    shadedErrorBar(t,tc_red_avrg_loc{pre}(:,iCon,iSize),tc_red_se_loc{pre}(:,iCon,iSize),'k');
    hold on
    shadedErrorBar(t,tc_red_avrg_loc{post}(:,iCon,iSize),tc_red_se_loc{post}(:,iCon,iSize),'b');
    ylim([-.02 .2]);
    hold on
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    ylabel('dF/F') 
    xlabel('s') 
    title(['+HTP',' n = ', num2str(nRed(iCon,iSize))])
    
    x0=5;
    y0=5;
    width=4;
    height=3;
    set(gcf,'units','inches','position',[x0,y0,width,height])
    set(gca,'XColor', 'none','YColor','none')
    
    sgtitle(['running, con ' num2str(cons(iCon)) ' size ' num2str(sizes(iSize))])
    
    print(fullfile(fnout,[num2str(cons(iCon)) '_' num2str(sizes(iSize)) 'matched_loc_cellType_timecourses.pdf']),'-dpdf');
    end
end 


%% population size tuning - averaged over contrast

%errorbar for stat resp and loc resp vs size, where error is across mice

%get average for each day

sizeResp_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
sizeResp_red_avrg_stat = cell(1,nd); %same for red
sizeResp_green_se_stat = cell(1,nd); %this will be the se across all green cells
sizeResp_red_se_stat = cell(1,nd); %same for red



for id = 1:nd
    green_data=squeeze(mean(pref_responses_stat_concat{id}(statGreen,:,:),2));%pulling the green cells and averaging over contrast
    sizeResp_green_avrg_stat{id}=nanmean(green_data,1);
    green_std=nanstd(green_data,1);
    sizeResp_green_se_stat{id}=green_std/sqrt(length(statGreen));
    
    
    red_data=squeeze(mean(pref_responses_stat_concat{id}(statRed,:,:),2));%pulling the red cells and averaging over contrast
    sizeResp_red_avrg_stat{id}=nanmean(red_data,1);
    red_std=nanstd(red_data,1);
    sizeResp_red_se_stat{id}=red_std/sqrt(length(statRed));
    
    clear green_std red_std green_data red_data
 
end



sizeResp_green_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
sizeResp_red_avrg_loc = cell(1,nd); %same for red
sizeResp_green_se_loc = cell(1,nd); %this will be the se across all green cells
sizeResp_red_se_loc = cell(1,nd); %same for red



for id = 1:nd
    green_data=squeeze(mean(pref_responses_loc_concat{id}(runningGreen,:,:),2));%pulling the green cells and averaging over contrast
    sizeResp_green_avrg_loc{id}=nanmean(green_data,1);
    green_std=nanstd(green_data,1);
    sizeResp_green_se_loc{id}=green_std/sqrt(length(runningGreen));
    
    red_data=squeeze(mean(pref_responses_loc_concat{id}(runningRed,:,:),2));%pulling the red cells and averaging over contrast
    sizeResp_red_avrg_loc{id}=nanmean(red_data,1);
    red_std=nanstd(red_data,1);
    sizeResp_red_se_loc{id}=red_std/sqrt(length(runningRed));
    
    clear green_std red_std green_data red_data
 
end



figure
subplot(2,2,1) %for the first day
errorbar(sizes,sizeResp_green_avrg_stat{pre},sizeResp_green_se_stat{pre},'k');
hold on
errorbar(sizes,sizeResp_green_avrg_stat{post},sizeResp_green_se_stat{post},'b');
title(['Stationary -HTP, DART',' n = ', num2str(length(statGreen))])
ylabel('dF/F, pref dir') 
xlabel('size (deg)') 
set(gca, 'TickDir', 'out')
box off
ylim([-.05 .2])

subplot(2,2,2) %for the second day
errorbar(sizes,sizeResp_red_avrg_stat{pre},sizeResp_red_se_stat{pre},'k');
hold on
errorbar(sizes,sizeResp_red_avrg_stat{post},sizeResp_red_se_stat{post},'b');
title(['Stationary +HTP',' n = ', num2str(length(statRed))])
ylabel('dF/F, pref dir') 
xlabel('size (deg)') 
set(gca, 'TickDir', 'out')
box off
ylim([-.05 .2])

subplot(2,2,3) %for the first day
errorbar(sizes,sizeResp_green_avrg_loc{pre},sizeResp_green_se_loc{pre},'k');
hold on
errorbar(sizes,sizeResp_green_avrg_loc{post},sizeResp_green_se_loc{post},'b');
title(['Running -HTP',' n = ', num2str(length(runningGreen))])
ylabel('dF/F, pref dir') 
xlabel('size (deg)') 
set(gca, 'TickDir', 'out')
box off
ylim([-.05 .3])

subplot(2,2,4) %for the second day
errorbar(sizes,sizeResp_red_avrg_loc{pre},sizeResp_red_se_loc{pre},'k');
hold on
errorbar(sizes,sizeResp_red_avrg_loc{post},sizeResp_red_se_loc{post},'b');
title(['Running +HTP',' n = ', num2str(length(runningRed))])
ylabel('dF/F, pref dir') 
xlabel('size (deg)') 
set(gca, 'TickDir', 'out')
box off
ylim([-.05 .3])

x0=5;
y0=5;
width=6;
height=6;
set(gcf,'units','inches','position',[x0,y0,width,height])



sgtitle(['population size tuning' ])

print(fullfile(fnout,['sizeTuning.pdf']),'-dpdf');


%% contrast response CONVERT TO MY FUNCTION
ymin=-0.01;
ymax=.15;
% errorbar for stat resp and loc resp vs size, where error is across mice
conResp_green_avrg_stat = cell(nSize,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_stat = cell(nSize,nd); %same for red
conResp_green_se_stat = cell(nSize,nd); %this will be the se across all green cells
conResp_red_se_stat = cell(nSize,nd); %same for red

consForPlotting = [25 50 100];
for id = 1:nd
    for iSize = 1:nSize
        green_data=pref_responses_stat_concat{id}(green_ind_concat,:,iSize);%pulling the green cells at this size
        conResp_green_avrg_stat{id}(iSize,:)=nanmean(green_data,1);
        green_std=nanstd(green_data,[],1);
        conResp_green_se_stat{id}(iSize,:)=green_std/sqrt(length(green_ind_concat));
        
        red_data=pref_responses_stat_concat{id}(red_ind_concat,:,iSize);%pulling the red cells at this size
        conResp_red_avrg_stat{id}(iSize,:)=nanmean(red_data,1);
        red_std=nanstd(red_data,[],1);
        conResp_red_se_stat{id}(iSize,:)=red_std/sqrt(length(red_ind_concat));
        
        clear green_std red_std green_data red_data
    end
end



figure
subplot(3,2,1) %for the first size, all contrasts
errorbar(consForPlotting,conResp_green_avrg_stat{pre}(1,:),conResp_green_se_stat{pre}(1,:),'--k');
hold on
errorbar(consForPlotting,conResp_green_avrg_stat{post}(1,:),conResp_green_se_stat{post}(1,:),'--b');
title(['Pyr n = ' , num2str(length(green_ind_concat))])
ylabel(['dF/F, ', num2str(sizes(1))]) 
set(gca, 'TickDir', 'out')
box off
ylim([ymin ymax-.12])
xlim([0 110])

subplot(3,2,2) 
errorbar(consForPlotting,conResp_red_avrg_stat{pre}(1,:),conResp_red_se_stat{pre}(1,:),'k');
hold on
errorbar(consForPlotting,conResp_red_avrg_stat{post}(1,:),conResp_red_se_stat{post}(1,:),'b');
title(['SST n = ' , num2str(length(red_ind_concat))])
set(gca, 'TickDir', 'out')
box off
ylim([ymin ymax-.12])
xlim([0 110])

subplot(3,2,3) %for the second size, all contrasts
errorbar(consForPlotting,conResp_green_avrg_stat{pre}(2,:),conResp_green_se_stat{pre}(2,:),'--k');
hold on
errorbar(consForPlotting,conResp_green_avrg_stat{post}(2,:),conResp_green_se_stat{post}(2,:),'--b');
ylabel(['dF/F, ', num2str(sizes(2))]) 
set(gca, 'TickDir', 'out')
box off
ylim([ymin ymax])
xlim([0 110])

subplot(3,2,4) 
errorbar(consForPlotting,conResp_red_avrg_stat{pre}(2,:),conResp_red_se_stat{pre}(2,:),'k');
hold on
errorbar(consForPlotting,conResp_red_avrg_stat{post}(2,:),conResp_red_se_stat{post}(2,:),'b');
set(gca, 'TickDir', 'out')
box off
ylim([ymin ymax])
xlim([0 110])

subplot(3,2,5) %for the second size, all contrasts
errorbar(consForPlotting,conResp_green_avrg_stat{pre}(3,:),conResp_green_se_stat{pre}(3,:),'--k');
hold on
errorbar(consForPlotting,conResp_green_avrg_stat{post}(3,:),conResp_green_se_stat{post}(3,:),'--b');
ylabel(['dF/F, ', num2str(sizes(3))]) 
set(gca, 'TickDir', 'out')
box off
ylim([ymin ymax])
xlim([0 110])

subplot(3,2,6) 
errorbar(consForPlotting,conResp_red_avrg_stat{pre}(3,:),conResp_red_se_stat{pre}(3,:),'k');
hold on
errorbar(consForPlotting,conResp_red_avrg_stat{post}(3,:),conResp_red_se_stat{post}(3,:),'b');
set(gca, 'TickDir', 'out')
box off
ylim([ymin ymax])
xlim([0 110])


han=axes('visible','off'); 
han.XLabel.Visible='on';
xlabel(han,'Contrast (%)');

x0=5;
y0=5;
width=4;
height=4;
set(gcf,'units','inches','position',[x0,y0,width,height])
sgtitle('Stationary')
print(fullfile(fnout,['contrastTuning.pdf']),'-dpdf');

%% for running 
ymin=-0.015;
ymax=.15;
% contrast response running
% errorbar for loc resp and loc resp vs size, where error is across mice
conResp_green_avrg_loc = cell(nSize,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_loc = cell(nSize,nd); %same for red
conResp_green_se_loc = cell(nSize,nd); %this will be the se across all green cells
conResp_red_se_loc = cell(nSize,nd); %same for red

consForPlotting = [12.5 25 50 100];
for id = 1:nd
    for iSize = 1:nSize
        green_data=pref_responses_loc_concat{id}(green_ind_concat,:,iSize);%pulling the green cells at this size
        conResp_green_avrg_loc{id}(iSize,:)=nanmean(green_data,1);
        green_std=nanstd(green_data,1);
        conResp_green_se_loc{id}(iSize,:)=green_std/sqrt(length(green_ind_concat));
        
        red_data=pref_responses_loc_concat{id}(red_ind_concat,:,iSize);%pulling the red cells at this size
        conResp_red_avrg_loc{id}(iSize,:)=nanmean(red_data,1);
        red_std=nanstd(red_data,1);
        conResp_red_se_loc{id}(iSize,:)=red_std/sqrt(length(red_ind_concat));
        
        clear green_std red_std green_data red_data
    end
end

conResp_green_avrg_loc = cell(nSize,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_loc = cell(nSize,nd); %same for red
conResp_green_se_loc = cell(nSize,nd); %this will be the se across all green cells
conResp_red_se_loc = cell(nSize,nd); %same for red

for id = 1:nd
    for iSize = 1:nSize
        green_data=pref_responses_loc_concat{id}(green_ind_concat,:,iSize);%pulling the green cells at this size
        conResp_green_avrg_loc{id}(iSize,:)=nanmean(green_data,1);
        green_std=nanstd(green_data,1);
        conResp_green_se_loc{id}(iSize,:)=green_std/sqrt(length(green_ind_concat));
        
        red_data=pref_responses_loc_concat{id}(red_ind_concat,:,iSize);%pulling the red cells at this size
        conResp_red_avrg_loc{id}(iSize,:)=nanmean(red_data,1);
        red_std=nanstd(red_data,1);
        conResp_red_se_loc{id}(iSize,:)=red_std/sqrt(length(red_ind_concat));
        
        clear green_std red_std green_data red_data
    end
end


figure
subplot(2,2,1) %for the first day
errorbar(consForPlotting,conResp_green_avrg_loc{pre}(1,:),conResp_green_se_loc{pre}(1,:),'--k');
hold on
errorbar(consForPlotting,conResp_green_avrg_loc{post}(1,:),conResp_green_se_loc{post}(1,:),'--b');
title(['Pyr n = ' , num2str(length(runningGreen))])
ylabel('dF/F, 20 deg') 
set(gca, 'TickDir', 'out')
box off
ylim([ymin ymax])
xlim([0 110])

subplot(2,2,3) %for the first day
errorbar(consForPlotting,conResp_green_avrg_loc{pre}(2,:),conResp_green_se_loc{pre}(2,:),'--k');
hold on
errorbar(consForPlotting,conResp_green_avrg_loc{post}(2,:),conResp_green_se_loc{post}(2,:),'--b');
ylabel('dF/F, Fullfield') 
set(gca, 'TickDir', 'out')
box off
ylim([ymin ymax])
xlim([0 110])

subplot(2,2,2) %for the first day
errorbar(consForPlotting,conResp_red_avrg_loc{pre}(1,:),conResp_red_se_loc{pre}(1,:),'k');
hold on
errorbar(consForPlotting,conResp_red_avrg_loc{post}(1,:),conResp_red_se_loc{post}(1,:),'b');
title(['SST n = ' , num2str(length(runningRed))])
set(gca, 'TickDir', 'out')
box off
ylim([ymin ymax])
xlim([0 110])

subplot(2,2,4) %for the first day
errorbar(consForPlotting,conResp_red_avrg_loc{pre}(2,:),conResp_red_se_loc{pre}(2,:),'k');
hold on
errorbar(consForPlotting,conResp_red_avrg_loc{post}(2,:),conResp_red_se_loc{post}(2,:),'b');
set(gca, 'TickDir', 'out')
box off
ylim([ymin ymax])
xlim([0 110])


han=axes('visible','off'); 
han.XLabel.Visible='on';
xlabel(han,'Contrast (%)');

x0=5;
y0=5;
width=6;
height=4;
set(gcf,'units','inches','position',[x0,y0,width,height])
sgtitle('Running')


print(fullfile(fnout,['contrastTuningRunning.pdf']),'-dpdf');


%% Normalized direction tuning

if length(dirs) ==8
    dirs_for_plotting=dirs-180;
elseif length(dirs) ==4
    dirs_for_plotting =dirs;
end

yMin = -.01;
yMax = .1;
yMax2 = .15;
iSize=3; %set to whichever size you want


green_dir_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
red_dir_avrg_stat = cell(1,nd); %same for red
green_dir_se_stat = cell(1,nd); %this will be the se across all green cells
red_dir_se_stat = cell(1,nd); %same for red

for iCon = 1:nCon
    for id = 1:nd
       
        green_dir_avrg_stat{id}=nanmean(nanmean(norm_dir_resp_concat{id}(statGreen,:,iCon,iSize),4),1);
        green_std=nanstd(nanmean(norm_dir_resp_concat{id}(statGreen,:,iCon,iSize),4),[],1);
        green_dir_se_stat{id}=green_std/sqrt(length(statGreen));
        green_dir_avrg_stat{id}=circshift(green_dir_avrg_stat{id},4);
        green_dir_se_stat{id}=circshift(green_dir_se_stat{id},4);
        
        red_dir_avrg_stat{id}=nanmean(nanmean(norm_dir_resp_concat{id}(statRed,:,iCon,iSize),4),1);
        red_std=nanstd(nanmean(norm_dir_resp_concat{id}(statRed,:,iCon,iSize),4),[],1);
        red_dir_se_stat{id}=red_std/sqrt(length(statRed));
        red_dir_avrg_stat{id}=circshift(red_dir_avrg_stat{id},4);
        red_dir_se_stat{id}=circshift(red_dir_se_stat{id},4);
        clear green_std red_std
        
    end
    
    
    
    
    figure
    subplot(1,2,1)
    errorbar(dirs_for_plotting,green_dir_avrg_stat{pre},green_dir_se_stat{pre},'--k')
    hold on
    errorbar(dirs_for_plotting,green_dir_avrg_stat{post},green_dir_se_stat{post},'--b')
    title('Stationary, Pyr')
    set(gca, 'TickDir', 'out')
    xticks(dirs)
    axis square
    box off
    ylabel('dF/F')
    ylim([yMin yMax])
    xlim([-10 140])
    
    subplot(1,2,2)
    errorbar(dirs_for_plotting,red_dir_avrg_stat{pre},red_dir_se_stat{pre},'k')
    hold on
    errorbar(dirs_for_plotting,red_dir_avrg_stat{post},red_dir_se_stat{post},'b')
   title('Stationary, SST')
    set(gca, 'TickDir', 'out')
    axis square
    box off
    xticks(dirs)
    ylim([yMin yMax])
    xlim([-10 140])
   
    sgtitle(['Normalized direction tuning ',num2str(cons(iCon))])
    
    
    print(fullfile(fnout,[num2str(cons(iCon)),'dirTuning.pdf']),'-dpdf','-bestfit')
end

%% change in pref direction
pref_dir_change=abs(pref_dir_concat{pre}-pref_dir_concat{post});
figure;subplot(1,2,1)
polarhistogram(pref_dir_change(green_ind_concat))
title('HTP-')
subplot(1,2,2)
polarhistogram(pref_dir_change(red_ind_concat))
title('HTP+')
print(fullfile(fnout,'prefDirChange.pdf'),'-dpdf','-bestfit')




%% plot fraction suppressed and facilitated


norm_diff_red = norm_diff_concat(:,:,:,red_ind_concat);
facil_red=norm_diff_red(:,:,:,:)>=1;
supp_red=norm_diff_red(:,:,:,:)<=-1;

N=length(red_ind_concat);



    facil_table_stat = squeeze(sum(facil_red(1,:,:,:),4)/N);
    supp_table_stat = squeeze(sum(supp_red(1,:,:,:),4)/N);
    
    figure;
    subplot(1,2,1)
    b=bar([1,2,3],[supp_table_stat(:,1),supp_table_stat(:,2),supp_table_stat(:,3)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
    b(1).FaceColor='k'
    b(2).FaceColor='r'
    b(3).FaceColor='b'
    xticklabels({'25','50','100'})
    hold on
    title('Suppressed')
    ylim([0 .4])
    ylabel(["Fraction SST cells"]) 
    xlabel(["Contrast"])
    set(gca,'TickDir','out')
    box off
    
    subplot(1,2,2)
    b=bar([1,2,3],[facil_table_stat(:,1),facil_table_stat(:,2),facil_table_stat(:,3)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
    b(1).FaceColor='k'
    b(2).FaceColor='r'
    b(3).FaceColor='b'
    xticklabels({,'25','50','100'})
    hold on
    title('Facilitated')
    ylim([0 .4])
    %ylabel(["Fraction SST cells"]) 
    xlabel(["Contrast"])
    set(gca,'TickDir','out')
    box off
    sgtitle('Stationary')
    
    x0=5;
    y0=5;
    width=3;
    height=1.75;
    set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'Facil_supp_stat.pdf'),'-dpdf');
%% 
norm_diff_red = norm_diff_concat(:,:,:,runningRed);
facil_red=norm_diff_red(:,:,:,:)>=1;
supp_red=norm_diff_red(:,:,:,:)<=-1;

N=length(runningRed);

 facil_table_loc = squeeze(sum(facil_red(2,:,:,:),4)/N);
 supp_table_loc = squeeze(sum(supp_red(2,:,:,:),4)/N);
    
    figure;
    subplot(1,2,1)
    b=bar([1,2,3],[supp_table_loc(:,1),supp_table_loc(:,2),supp_table_loc(:,3)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
    b(1).FaceColor='k'
    b(2).FaceColor='r'
    b(3).FaceColor='b'
    xticklabels({'25','50','100'})
    hold on
    title('Suppressed')
    ylim([0 .6])
    ylabel(["Fraction SST cells"]) 
    xlabel(["Contrast"])
    set(gca,'TickDir','out')
    box off
    
    subplot(1,2,2)
    b=bar([1,2,3],[facil_table_loc(:,1),facil_table_loc(:,2),facil_table_loc(:,3)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
    b(1).FaceColor='k'
    b(2).FaceColor='r'
    b(3).FaceColor='b'
    xticklabels({'25','50','100'})
    hold on
    title('Facilitated')
    ylim([0 .6])
    %ylabel(["Fraction SST cells"]) 
    xlabel(["Contrast"])
    set(gca,'TickDir','out')
    box off
    sgtitle('Running')
    x0=5;
    y0=5;
    width=3;
    height=1.75;
    set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'Facil_supp_loc.pdf'),'-dpdf');

%%
%make a table of suppresses and facilitated cells for the cells that have
%both stationary and running within each condition
facil=norm_diff(:,:,:,:)>=1;
supp=norm_diff(:,:,:,:)<=-1;

supp_table_stat=nan(nCon,nSize);
facil_table_stat=nan(nCon,nSize);
supp_table_loc=nan(nCon,nSize);
facil_table_loc=nan(nCon,nSize);

for iCon = 1:nCon
    for iSize=1:nSize
        theseRedCells = intersect(red_ind_concat, find(runningByCondition(:,iCon,iSize)));
        supp_table_stat(iCon,iSize)=sum(supp(1,iCon,iSize,theseRedCells),4)/nRed(iCon,iSize);
        supp_table_loc(iCon,iSize)=sum(supp(2,iCon,iSize,theseRedCells),4)/nRed(iCon,iSize);

        facil_table_stat(iCon,iSize)=sum(facil(1,iCon,iSize,theseRedCells),4)/nRed(iCon,iSize);
        facil_table_loc(iCon,iSize)=sum(facil(2,iCon,iSize,theseRedCells),4)/nRed(iCon,iSize);


    end
end


figure;
subplot(2,2,1)
b=bar([1,2,3],[supp_table_stat(:,1),supp_table_loc(:,1)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
b(1).FaceColor="#70D0F6"
b(2).FaceColor="#0C8ABB"
xticklabels({'12.5','25','50','100'})
hold on
title('Suppressed 20-deg')
ylim([0 .6])
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast"])
set(gca,'TickDir','out')
box off

subplot(2,2,2)
b=bar([1,2,3],[facil_table_stat(:,1),facil_table_loc(:,1)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
b(1).FaceColor="#C983B1"
b(2).FaceColor="#883367"
xticklabels({'12.5','25','50','100'})
hold on
title('Facilitated 20-deg')
ylim([0 .6])
%ylabel(["Fraction SST cells"]) 
xlabel(["Contrast"])
set(gca,'TickDir','out')
box off



subplot(2,2,3)
b=bar([1,2,3,4],[supp_table_stat(:,2),supp_table_loc(:,2)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
b(1).FaceColor="#70D0F6"
b(2).FaceColor="#0C8ABB"
xticklabels({'12.5','25','50','100'})
hold on
title('Suppressed fullfield')
ylim([0 .6])
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast"])
set(gca,'TickDir','out')
box off

subplot(2,2,4)
b=bar([1,2,3,4],[facil_table_stat(:,2),facil_table_loc(:,2)],'grouped','FaceColor',"#00AFEF",'EdgeColor', [1 1 1])
b(1).FaceColor="#C983B1"
b(2).FaceColor="#883367"
xticklabels({'12.5','25','50','100'})
hold on
title('Facilitated fullfield')
ylim([0 .6])
%ylabel(["Fraction SST cells"]) 
xlabel(["Contrast"])
set(gca,'TickDir','out')
box off
x0=5;
y0=5;
width=4;
height=4;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'matched_facil_supp_byState.pdf'),'-dpdf');

%% norm diff for stationary vs running
for iSize = 1:nSize
    figure;
    subplot(1,2,1)
    boxchart(squeeze(norm_diff(1,:,iSize,red_ind_concat))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75],BoxEdgeColor=[0 0 0]);
    hold on
    scatter([1, 2, 3],squeeze(norm_diff(1,:,iSize,red_ind_concat))',20,[.5 .15 .20], 'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)
    xticklabels({'25','50','100'})
    xlabel('Contrast(%)')
    ylabel('Normalized difference')
    ylim([-8 8])
    title('Stationary')
    hold off
    set(gca,'TickDir','out')
    box off
    x0=5;
    y0=5;
    width=1.25;
    height=2;
    set(gcf,'units','inches','position',[x0,y0,width,height])
    
    
    subplot(1,2,2)
    boxchart(squeeze(norm_diff(2,:,iSize,red_ind_concat))',MarkerStyle ="none",BoxFaceColor=	[.75 .75 .75],BoxEdgeColor=[0 0 0]);
    hold on
    scatter([1, 2, 3],squeeze(norm_diff(2,:,iSize,red_ind_concat))',20,[.5 .15 .20], 'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)
    xticklabels({'25','50','100'})
    xlabel('Contrast(%)')
    %ylabel('Normalized difference')
    ylim([-12 12])
    title('Running')
    hold off
    set(gca,'TickDir','out')
    box off
    x0=5;
    y0=5;
    width=6;
    height=4;
    set(gcf,'units','inches','position',[x0,y0,width,height])
    
    sgtitle(num2str(sizes(iSize)))

    print(fullfile(fnout,['normDiff_scatter_size', num2str(sizes(iSize)), '.pdf']),'-dpdf');
end

%% scatterplot and ttest for interneurons

%using cells matched across running and stationary for each stimulus
%condition

for iSize = 1:nSize
    for iCon = 1:nCon
        theseRedCells = intersect(red_ind_concat, find(runningByCondition(:,iCon,iSize)));
        figure; 
        subplot(1,2,1)
        scatter((pref_responses_stat_concat{pre}(theseRedCells,iCon,iSize)),(pref_responses_stat_concat{post}(theseRedCells,iCon,iSize)))
        set(gca, 'TickDir', 'out')
        axis square
        box off
        title(['Stationary con ', num2str(cons(iCon)), ' size ', num2str(sizes(iSize))])
        xlim([0 .2])
        ylim([0 .2])
        refline(1)
        [h,p]=ttest((pref_responses_stat_concat{pre}(theseRedCells,iCon,iSize)),(pref_responses_stat_concat{post}(theseRedCells,iCon,iSize)));
        if (p) < 0.05
            txt = {'p = ' num2str(round(p,3)),'n = ' num2str(length(theseRedCells))};
            text(.1,0.1,txt)
        end

        subplot(1,2,2)
        scatter((pref_responses_loc_concat{pre}(theseRedCells,iCon,iSize)),(pref_responses_loc_concat{post}(theseRedCells,iCon,iSize)))
        set(gca, 'TickDir', 'out')
        axis square
        box off
        title(['Running con ', num2str(cons(iCon)), ' size ', num2str(sizes(iSize))])
        xlim([0 .2])
        ylim([0 .2])
        refline(1)
        [h,p]=ttest((pref_responses_loc_concat{pre}(theseRedCells,iCon,iSize)),(pref_responses_loc_concat{post}(theseRedCells,iCon,iSize)));
        if (p) < 0.05
            txt = {'p = ' num2str(round(p,3)),'n = ' num2str(length(theseRedCells))};
            text(.15,0.15,txt)
        end

    end
end

%% plot large and small pupil stationary timecourses
%--> convert to using the function
% make figure with se shaded, one figure per contrast

tc_green_avrg_small = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_small = cell(1,nd); %same for red
tc_green_se_small = cell(1,nd); %this will be the se across all green cells
tc_red_se_small = cell(1,nd); %same for red


nGreen = nan(nCon,nSize);
nRed = nan(nCon,nSize);

for id = 1:nd
  for iCon = 1:nCon
      for iSize = 1:nSize
        tc_green_avrg_small{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_stat_smallPupil_concat{id}(:,green_ind_concat,iCon,iSize),2);
        green_std=nanstd(tc_trial_avrg_stat_smallPupil_concat{id}(:,green_ind_concat,iCon,iSize),[],2);
        tc_green_se_small{id}(:,iCon,iSize)=green_std/sqrt(length(green_ind_concat));
        
        tc_red_avrg_small{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_stat_smallPupil_concat{id}(:,red_ind_concat,iCon,iSize),2);
        red_std=nanstd(tc_trial_avrg_stat_smallPupil_concat{id}(:,red_ind_concat,iCon,iSize),[],2);
        tc_red_se_small{id}(:,iCon,iSize)=red_std/sqrt(length(red_ind_concat));
        
        
        clear green_std red_std
      end 
    end
end

z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_small{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
    for iSize = 1:nSize
    figure
    subplot(1,2,1) %for the first day
    
    
    
    ylim([-.05 .25]);
    hold on
    shadedErrorBar(t,tc_green_avrg_small{pre}(:,iCon,iSize),tc_green_se_small{pre}(:,iCon,iSize),'--k');
    hold on
    shadedErrorBar(t,tc_green_avrg_small{post}(:,iCon,iSize),tc_green_se_small{post}(:,iCon,iSize),'--b','transparent');
    hold on
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    title(['-HTP',' n = ', num2str(nGreen(iCon,iSize))])
    
    ylabel('dF/F') 
    xlabel('s') 
    set(gca,'XColor', 'none','YColor','none')
    
    
    subplot(1,2,2) %+HTP
    shadedErrorBar(t,tc_red_avrg_small{pre}(:,iCon,iSize),tc_red_se_small{pre}(:,iCon,iSize),'k');
    hold on
    shadedErrorBar(t,tc_red_avrg_small{post}(:,iCon,iSize),tc_red_se_small{post}(:,iCon,iSize),'b');
    ylim([-.05 .25]);
    hold on
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    ylabel('dF/F') 
    xlabel('s') 
    title(['+HTP',' n = ', num2str(nRed(iCon,iSize))])
    
    x0=5;
    y0=5;
    width=4;
    height=3;
    set(gcf,'units','inches','position',[x0,y0,width,height])
    set(gca,'XColor', 'none','YColor','none')
    
    sgtitle(['small pupil, con ' num2str(cons(iCon)) ' size ' num2str(sizes(iSize))])
    
    print(fullfile(fnout,[num2str(cons(iCon)) '_' num2str(sizes(iSize)) 'matched_small_cellType_timecourses.pdf']),'-dpdf');
    end
end 


tc_green_avrg_large = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_large = cell(1,nd); %same for red
tc_green_se_large = cell(1,nd); %this will be the se across all green cells
tc_red_se_large = cell(1,nd); %same for red


nGreen = nan(nCon,nSize);
nRed = nan(nCon,nSize);

for id = 1:nd
  for iCon = 1:nCon
      for iSize = 1:nSize
        tc_green_avrg_large{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_stat_largePupil_concat{id}(:,green_ind_concat,iCon,iSize),2);
        green_std=nanstd(tc_trial_avrg_stat_largePupil_concat{id}(:,green_ind_concat,iCon,iSize),[],2);
        tc_green_se_large{id}(:,iCon,iSize)=green_std/sqrt(length(green_ind_concat));
        
        tc_red_avrg_large{id}(:,iCon,iSize)=nanmean(tc_trial_avrg_stat_largePupil_concat{id}(:,red_ind_concat,iCon,iSize),2);
        red_std=nanstd(tc_trial_avrg_stat_largePupil_concat{id}(:,red_ind_concat,iCon,iSize),[],2);
        tc_red_se_large{id}(:,iCon,iSize)=red_std/sqrt(length(red_ind_concat));
        
        
        clear green_std red_std
      end 
    end
end

z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_large{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

for iCon = 1:nCon
    for iSize = 1:nSize
    figure
    subplot(1,2,1) %for the first day
    
    
    
    ylim([-.05 .25]);
    hold on
    shadedErrorBar(t,tc_green_avrg_large{pre}(:,iCon,iSize),tc_green_se_large{pre}(:,iCon,iSize),'--k');
    hold on
    shadedErrorBar(t,tc_green_avrg_large{post}(:,iCon,iSize),tc_green_se_large{post}(:,iCon,iSize),'--b','transparent');
    hold on
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    title(['-HTP',' n = ', num2str(nGreen(iCon,iSize))])
    
    ylabel('dF/F') 
    xlabel('s') 
    set(gca,'XColor', 'none','YColor','none')
    
    
    subplot(1,2,2) %+HTP
    shadedErrorBar(t,tc_red_avrg_large{pre}(:,iCon,iSize),tc_red_se_large{pre}(:,iCon,iSize),'k');
    hold on
    shadedErrorBar(t,tc_red_avrg_large{post}(:,iCon,iSize),tc_red_se_large{post}(:,iCon,iSize),'b');
    ylim([-.05 .25]);
    hold on
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    ylabel('dF/F') 
    xlabel('s') 
    title(['+HTP',' n = ', num2str(nRed(iCon,iSize))])
    
    x0=5;
    y0=5;
    width=4;
    height=3;
    set(gcf,'units','inches','position',[x0,y0,width,height])
    set(gca,'XColor', 'none','YColor','none')
    
    sgtitle(['large pupil, con ' num2str(cons(iCon)) ' size ' num2str(sizes(iSize))])
    
    print(fullfile(fnout,[num2str(cons(iCon)) '_' num2str(sizes(iSize)) 'matched_large_cellType_timecourses.pdf']),'-dpdf');
    end
end 
%% Split cells by noise correlation and examine DART effects

%within each cell of noiseCorr_concat is a 2 X nNeurons matrix, where on the
%first dimension 1=R and 2=p values.
figure;histogram(noiseCorr_concat{pre}(1,red_ind_concat));xlim([-.1 1])
figure;histogram(noiseCorr_concat{pre}(1,green_ind_concat));xlim([-.2 1.2])
figure;cdfplot(noiseCorr_concat{pre}(1,red_ind_concat))

median(noiseCorr_concat{pre}(1,red_ind_concat))

highNoiseCorr_red=find(noiseCorr_concat{pre}(1,red_ind_concat)>.3);
lowNoiseCorr_red=find(noiseCorr_concat{pre}(1,red_ind_concat)<=.3);


plotNeuralTimecourse(tc_trial_avrg_stat_concat, tc_trial_avrg_stat_concat, ...
    highNoiseCorr_red, highNoiseCorr_red, ...
    'UseDashedLines', [false, false], ...
    'Colors1', {'k', 'b'}, ...  % Black for pre, blue for post on left plots
    'Colors2', {'k', 'b'}, ...  % Black for pre, blue for post on right plots
    'Titles', {'lowCorr', 'highCorr'}, ...
    'StimStart', 31);
%% 
mean_pref_resp = mean(mean(pref_responses_stat_concat{pre},2),3);
mean_norm_diff = squeeze(squeeze(mean(mean(norm_diff_concat(1,:,:,:),2),3)));


figure;
subplot(1,3,1);scatter(mean_pref_resp(red_ind_concat),noiseCorr_concat{pre}(1,red_ind_concat));xlabel('mean resp');ylabel('noise corr');xlim([-.1 .6]);ylim([-.2 1])
subplot(1,3,2);scatter(mean_pref_resp(red_ind_concat),mean_norm_diff(red_ind_concat));xlabel('mean resp');ylabel('mean norm diff');xlim([-.1 .6]);ylim([-6 4])
subplot(1,3,3);scatter(noiseCorr_concat{pre}(1,red_ind_concat),mean_norm_diff(red_ind_concat));xlabel('noise corr');ylabel('mean norm diff');xlim([-.2 1]);ylim([-6 4])

x0=5;
y0=0;
width=10;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
sgtitle('LM')
%% 
histogram(mean_pref_resp(red_ind_concat))
