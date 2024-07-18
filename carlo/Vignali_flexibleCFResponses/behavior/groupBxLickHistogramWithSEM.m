%% groupWiseBxData

close all;
clear;

bxSourceBase = 'A:\home\carlo\RC\rawData\behavioral\'; %base folder for data from bx session
bxOutputBase = 'A:\home\carlo\RC\analysis\behavioral\'; %base folder for output for bx analysis
crpFigBase = 'A:\home\carlo\RC\analysis\crpFigures\'; %base folder for output for figures

% % only used when running bx analysis on mac (also need to disable saveas for some reason
% bxSourceBase = '/Volumes/All_staff/home/carlo/rawData/behavioral/';
% bxOutputBase = '/Volumes/All_staff/home/carlo/analysis/behavioral';
% crpFigBase = '/Volumes/All_staff/home/carlo/analysis/crpFigures/';

timeBeforeCue = 1000; %defines the window around the cue presentation which will be taken for plotting
timeAfterCue = 3000;

crpTrainDayList; %load blank variable arrays
%dateIdx = input('My mice :7: | Mikes mice :6: ');

%check and make sure the figure destinations exist
    %pull out mouse and date info from bx naming
[thisMouse, thisDate] = selectMouseInfo(days{1}); 
sessionName = input('which session; naive, PL, etc.?: ','s');
    %name folder for session figures
sessionFigs = [crpFigBase, sessionName '_session\'];
    %name folder for xSession figures
summaryFigs = [crpFigBase, sessionName, '_summary\'];
    %check if exists, if not, make file
if exist([crpFigBase, sessionName, '_session\'], 'file') == 0
    mkdir([crpFigBase, sessionName, '_session\']);
end
if exist([crpFigBase, sessionName, '_summary\'], 'file') == 0
    mkdir([crpFigBase, sessionName, '_summary\']);
end
    %determines the gap in days between training sessions (if a day is
        %skipped it returns (thatDate+0.5)
gapInx = gapDayFinder(days);

%% SECTION TWO - 

for ii = 1:length(days)
    %days{ii}
        %names output directory for bx data
    bxOutputDir  = [crpFigBase  'grouped_bxOutput'];
        %uses getBxData.m to find path to bx data matching thisMouse and
            %thisDate and loads the file that matches (or returns error if
            %multiple exist with the same name)
    bxData = getBxData_Carlo(bxSourceBase, days{ii});  %find the correct behavior file and loads it.
        %pulls out mouse ID and date of bx session
    [thisMouse, thisDate] = selectMouseInfo(days{ii});
        %finds each trial's start time in free floating MWorks time
    trialStart{ii,1} = round(double(cell2mat(bxData.tThisTrialStartTimeMs))); 
        %gets vectors for # of trials. Includes unimaged trials.
    numTrials{ii,1}  = length(bxData.tThisTrialStartTimeMs);
        %stores the time of the beginning of the first trial in MWorks time.
    bxStartMWorksTime{ii,1}  = trialStart{ii,1}(1,1);
    
        %use time of first frame to align licking times to start of imaging
    initLickTimes = zeros(size(bxData.lickometerTimesUs{1,length(bxData.lickometerTimesUs)}));
        %determines if 2nd var is in struct[1st var]
    if isfield(bxData, 'lickometerTimesUs') 
        for kk = 1:length(bxData.lickometerTimesUs)
                %concatenates licking times
            initLickTimes = [initLickTimes cell2mat(bxData.lickometerTimesUs(kk))/1000]; 
        end
            %aligns lick times to the start of imaging and converts from int64 to double array
        lickTimes{ii,1} = double(initLickTimes)-(bxStartMWorksTime{ii,1}); 
    end
    
    %%%Collects various events during session
        %aligns start of hold to start of imaging for each cue (cue onset)
    hold_start{ii,1} = double(cell2mat(bxData.holdStartsMs)) - bxStartMWorksTime{ii,1}; 
        %duration of the "lever hold" on that trial. meaningless here except to calculate cue onset
    hold_time{ii,1}  = double(cell2mat(bxData.holdTimesMs)); 
        %time between cue and reward [?]
    react_time{ii,1} = double(cell2mat(bxData.reactTimesMs)); 
        %fixed value throughout - required hold time for reward(?)
    req_hold{ii,1}   = double(cell2mat(bxData.tTotalReqHoldTimeMs)); 
        %fixed value throughout - rand add time req to hold (?)
    rnd_hold{ii,1}   = double(cell2mat(bxData.tRandReqHoldTimeMs)); 
        %sum of req hold times - useless in this experiment
    tot_req_hold{ii,1} = req_hold{ii,1} + rnd_hold{ii,1};
        %cue onset aligned to the start of imaging
    release_time{ii,1} = hold_start{ii,1} + hold_time{ii,1};
        %cumHistograms are aligned to cue - time from cue onset to licking
    cuePresentation{ii,1} = release_time{ii,1} - react_time{ii,1};  
        %this variable lengthens the time that the cue is on the screen
    TFT_ms{ii,1} = bxData.tooFastTimeMs; 
        %total interval between cue onset and reward delivery
    cue_rew_int{ii,1} = bxData.RewardDelayDurationMs + round(mean(react_time{ii,1}),-1); 
    
    block2_inx{ii,1} = find(cell2mat(bxData.tBlock2TrialNumber(1:end-1))); %isolates nonzero trial #s

    if exist(bxOutputDir)
            %updates saved bx data if additional data is analyized
        save(bxOutputDir, 'lickTimes', '-append');
    else
            %if no bx data has been saved, save as a .mat
        save(bxOutputDir, 'lickTimes');
    end
end
%% reorganize variables so they can be cat and finish any final math

for i = 1:length(days)
    trialStart{i,1} = trialStart{i,1} - bxStartMWorksTime{i,1};
    min_start_to_cue{i,1} = min([cuePresentation{i,1}-trialStart{i,1}]);
    min_cue_to_end{i,1} = min([trialStart{i,1}(2:end)-cuePresentation{i,1}(1:end-1)]);



%%
   bin_size = 100; %number of ms to bin licking. 
%     if min_start_to_cue > 20000
%         pre_cue_window_lick{i,1} = 20000;
%     else 
%         pre_cue_window_lick{i,1} = floor([min_start_to_cue/bin_size])*bin_size;
%     end
%     if min_cue_to_end > 6000
%         post_cue_window_lick{i,1} = 5999;
%     else 
%         post_cue_window_lick{i,1} = floor([min_cue_to_end/bin_size])*bin_size-1;
%     end
    pre_cue_window_lick = 18000;
    post_cue_window_lick = 5999;
    
    %get lick traces (1ms resolution) for rewarded and omission trials
    full_trial_licks{i,1} = zeros(length(cuePresentation{i,1})-1,(pre_cue_window_lick+post_cue_window_lick+1)); %dim1=trial# dim2=ms
    for kk = 1:length(cuePresentation{i,1})-1 %look at all trials except the last one.
        %find all the licking events Xms before and Yms after cue presentation
        licks_this_window{i,1} = lickTimes{i,1}(lickTimes{i,1}>cuePresentation{i,1}(kk)-pre_cue_window_lick & lickTimes{i,1}<cuePresentation{i,1}(kk)+post_cue_window_lick);
        alignment_this_trial{i,1} = cuePresentation{i,1}(kk)-(pre_cue_window_lick+1); %subtract off this time so that the lick times are converted to numbers which will correspond to there index in licks_by_trial
        licks_this_window{i,1} = licks_this_window{i,1} - alignment_this_trial{i,1};
        full_trial_licks{i,1}(kk, (licks_this_window{i,1})) = 1;
    end
    
    if bxData.doBlock2 == 0
        full_trial_licks_rewarded{i,1} = full_trial_licks{i,1};
%         full_trial_licks_rewarded{i,1}(sort([reward_omit_inx, unexp_rew_inx]) , :) = []; %exclude omission and unexpected trials from rewarded trials condition
    elseif bxData.doBlock2 ==1
        full_trial_licks_rewarded{i,1} = full_trial_licks{i,1};
        full_trial_licks_rewarded{i,1}(sort([block2_inx{i,1}]) , :) = []; %exclude omission and unexpected trials from rewarded trials condition
        full_trial_licks_block2{i,1} = full_trial_licks{i,1}(block2_inx{i,1}, :);
    end
%     full_trial_licks_omit = full_trial_licks(reward_omit_inx, :);
%     full_trial_licks_unexp = full_trial_licks(unexp_rew_inx, :);
     
    %stderr of trial licks
    full_trial_licks_rewarded_ste{i,1} = std(full_trial_licks_rewarded{i,1},1)./sqrt(length(full_trial_licks_rewarded{i,1}));
    full_trial_licks_rewarded_bin_ste{i,1} = zeros(1,(length(full_trial_licks_rewarded_ste{i,1})/bin_size));
    %bin licking by 50ms bins and convert to licks/sec
    full_trial_licks_rewarded_sum{i,1} = sum(full_trial_licks_rewarded{i,1},1);
%     full_trial_licks_omit_sum = sum(full_trial_licks_omit,1);
%     full_trial_licks_unexp_sum = sum(full_trial_licks_unexp,1);
    full_trial_licks_rewarded_bin{i,1} = zeros(1,(length(full_trial_licks_rewarded_sum{i,1})/bin_size));
%     full_trial_licks_omit_bin = zeros(1,(length(full_trial_licks_omit_sum)/bin_size));
%     full_trial_licks_unexp_bin = zeros(1,(length(full_trial_licks_unexp_sum)/bin_size));
    if bxData.doBlock2
        full_trial_licks_block2_sum{i,1} = sum(full_trial_licks_block2{i,1},1);
        full_trial_licks_block2_ste{i,1} = std(full_trial_licks_block2{i,1},1)./sqrt(length(full_trial_licks_block2{i,1}));
        full_trial_licks_block2_bin{i,1} = zeros(1,(length(full_trial_licks_block2_sum{i,1})/bin_size));
        full_trial_licks_block2_bin_ste{i,1} = zeros(1,(length(full_trial_licks_block2_ste{i,1})/bin_size));
    end
    cue_presentation_binned{i,1} = (pre_cue_window_lick/bin_size)+1;
    iii=1;
    for kk = 1:bin_size:length(full_trial_licks_rewarded_sum{i,1})
        iii=iii+1;
        full_trial_licks_rewarded_bin{i,1}(iii) = sum(full_trial_licks_rewarded_sum{i,1}(kk:[kk+bin_size-1]));
        full_trial_licks_rewarded_bin_ste{i,1}(iii) = sum(full_trial_licks_rewarded_ste{i,1}(kk:[kk+bin_size-1]));
%         full_trial_licks_omit_bin(iii) = sum(full_trial_licks_omit_sum(kk:[kk+bin_size-1]));
%         full_trial_licks_unexp_bin(iii) = sum(full_trial_licks_unexp_sum(kk:[kk+bin_size-1]));
        if bxData.doBlock2 ==1
            full_trial_licks_block2_bin{i,1}(iii) = sum(full_trial_licks_block2_sum{i,1}(kk:[kk+bin_size-1]));
            full_trial_licks_block2_bin_ste{i,1}(iii) = sum(full_trial_licks_block2_ste{i,1}(kk:[kk+bin_size-1]));
        end
    end
    
    %plot rewarded trials
    full_trial_licks_rewarded_bin{i,1} = (full_trial_licks_rewarded_bin{i,1}/size(full_trial_licks_rewarded{i,1},1))*(1000/bin_size); % convert to lick rate in Hz
    full_trial_licks_rewarded_bin_ste{i,1} = (full_trial_licks_rewarded_bin_ste{i,1}/size(full_trial_licks_rewarded_ste{i,1},1))*(1000/bin_size); % convert to lick rate in Hz
    full_trial_licks_block2_bin{i,1} = (full_trial_licks_block2_bin{i,1}/size(full_trial_licks_rewarded{i,1},1))*(1000/bin_size); % convert to lick rate in Hz
    full_trial_licks_block2_bin_ste{i,1} = (full_trial_licks_block2_bin_ste{i,1}/size(full_trial_licks_rewarded_ste{i,1},1))*(1000/bin_size); % convert to lick rate in Hz
    x_axis_bin{i,1} = ([1:length(full_trial_licks_rewarded_bin{i,1})]-cue_presentation_binned{i,1})*(bin_size/1000);
end

%% reformat each cell into nAnimals:nBin
x_axis_binGroup = nanmean(cell2mat(x_axis_bin),1);
full_trial_licks_rewarded_binGroup = nanmean(cell2mat(full_trial_licks_rewarded_bin),1);
full_trial_licks_block2_binGroup = nanmean(cell2mat(full_trial_licks_block2_bin),1);
full_trial_licks_rewarded_bin_steGroup = nanmean(cell2mat(full_trial_licks_rewarded_bin_ste),1);
full_trial_licks_block2_bin_steGroup = nanmean(cell2mat(full_trial_licks_block2_bin_ste),1);

for i=1:length(full_trial_licks_rewarded); nLickTrials{i,1} = size(full_trial_licks_rewarded{i,1},1); end
tLickTrials = sum(cell2mat(nLickTrials));

%% figure
     
   
    %if exist([sessionFigs thisDate '_img' thisMouse '_stackedLickHistSEM'], 'file') == 0
% lick histogram across CS+ and CS- trials in a session
tempFig = setFigure_Bx; hold on; 
bar(x_axis_binGroup, full_trial_licks_rewarded_binGroup, 'black'); hold on; 
erR = errorbar(x_axis_binGroup,full_trial_licks_rewarded_binGroup,(-full_trial_licks_rewarded_bin_steGroup(:,:)),full_trial_licks_rewarded_bin_steGroup);    
erR.Color = ('k'); erR.LineStyle = 'none'; 
% bar(x_axis_binGroup, full_trial_licks_block2_binGroup, 'red'); 
% erB2 = errorbar(x_axis_binGroup,full_trial_licks_block2_binGroup,(-full_trial_licks_block2_bin_steGroup(:,:)),full_trial_licks_block2_bin_steGroup);    
% erB2.Color = ('red'); erB2.LineStyle = 'none';  
title([sessionName ': rewarded trials']); xlabel('time (s) relative to cue'); ylabel('lick rate (Hz)'); 
xline(0,'k'); xline(770./1000,'k'); %rew_cue_int is an int64 variable so this simple division isn't working need to convert
axis([-1 3 0 10]); xticks([-1 0 1 2 3]); yticks([0 2 4 6 8 10]);

savefig([sessionFigs '_csplusstackedLickHistSEM']);
saveas(tempFig, [sessionFigs '_CSPlusStackedLickHistSEM.pdf']);

tempFig = setFigure_Bx; hold on; 
% bar(x_axis_binGroup, full_trial_licks_rewarded_binGroup, 'black'); hold on; 
% erR = errorbar(x_axis_binGroup,full_trial_licks_rewarded_binGroup,(-full_trial_licks_rewarded_bin_steGroup(:,:)),full_trial_licks_rewarded_bin_steGroup);    
% erR.Color = ('k'); erR.LineStyle = 'none'; 
bar(x_axis_binGroup, full_trial_licks_block2_binGroup, 'red'); 
erB2 = errorbar(x_axis_binGroup,full_trial_licks_block2_binGroup,(-full_trial_licks_block2_bin_steGroup(:,:)),full_trial_licks_block2_bin_steGroup);    
erB2.Color = ('red'); erB2.LineStyle = 'none';  
title([sessionName ': unrewarded trials']); xlabel('time (s) relative to cue'); ylabel('lick rate (Hz)'); 
xline(0,'k'); xline(770./1000,'k'); %rew_cue_int is an int64 variable so this simple division isn't working need to convert
axis([-1 3 0 10]); xticks([-1 0 1 2 3]); yticks([0 2 4 6 8 10]);

savefig([sessionFigs '_csminusstackedLickHistSEM']);
saveas(tempFig, [sessionFigs '_CSMinusStackedLickHistSEM.pdf']);

    
    
