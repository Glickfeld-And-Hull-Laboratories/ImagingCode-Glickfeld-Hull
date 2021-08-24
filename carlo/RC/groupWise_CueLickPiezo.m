clear
close all

analysis_out = 'A:\home\carlo\analysis\2P\';
% analysis_out{1,2} = 'A:\home\carlo\mikeAnalysis\2P\';

bdata_source = 'A:\home\carlo\rawData\behavioral\';
% bdata_source{1,2} = 'A:\home\carlo\mikeAnalysis\Behavior\';
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
% it's gonna have to be seperated again, so if animal ID > 1700 go to
% /analysis/2P if animal ID < 1700 go to /mikeAnalysis/2P
%% pull out all neural data across animals based on ID (<1600 vs. >1700)
output_fn = ([analysis_out 'groupwiseAlign' sesh '_' sesh2 '\']);
if ~exist(fullfile(output_fn))
 mkdir(fullfile(output_fn))
end

loadData = input('Reload dataset(1) or open current set(2): ');
if loadData == 2
    if seshType ~= 2
    load(fullfile(output_fn, '_groupwiseDataset.mat'));
    elseif seshType == 2
    blockName = input('which block is this? [type "plus/minus""First/Second"]: ','s');
    load(fullfile(output_fn, [blockName '_groupwiseDataset.mat']));
    end
else

%% call blanks for bx data concat
% gpreRew_lickBurstHz = []; gpostRew_lickBurstHz = []; glickBurstHz_all = []; gpostRew_lickBurstStart = []; gpostRew_lickAlignEvents = []; gpostRew_lickAlign = []; gpostRewTrials = []; gpreRewTrials = []; gfirstPostRew_lickAlignEvents = []; gfirstPostRew_lickAlign = []; glastPreRew_lickAlignEvents = []; glastPreRew_lickAlign = []; gTargetOn = []; gRctT = [];  gind_block2 = []; gind_rew = []; ggoodcells = [];
gind_block2 = []; gind_rew = []; gTargetOn = []; gRctT = []; initFrames = []; ggoodcells = [];gtt = [];
initLickTimes = []; initCounterTimes = []; initCounterVals = [];
ind_early_bst=[]; ind_late_bst=[]; early_bsttime=[]; late_bsttime=[];
earlyBurstEvents=[]; lateBurstEvents=[]; earlyBurstStart=[]; lateBurstStart=[];
%% 
if seshType == 2
allBlock = input('all block sessions (1-yes): ');
if allBlock == 1
    exptCount = 2.*length(expt);
else
    exptCount = length(expt);
end
else
    exptCount = length(expt);
end

for mouseIdx = 1:exptCount
[pairings] = loadRCList_Carlo(expt,seshType);
    %1512, 1513, 1516, 1520, 1081, 1702, 1705, 1706 
    % session number always 1 -- seperate PL and Naive lists
iexp = pairings{1,1}; isesh = pairings{3,1};
mouse = strtrim(expt{1,iexp}.mouse{1,isesh});
date = expt{1,iexp}.date{1,isesh};
run = expt{1,iexp}.run{1,isesh};
sessions = [date '_img' mouse];
%% this if statment is specific to this mouse as two lobules were imaged

img_fn = [date '_img' mouse '\getTC_' run '\'];

%%
if str2double(regexp(mouse,'\d*','match')) < 1700 && str2double(regexp(mouse,'\d*','match')) ~= 1081
analysis_out = 'A:\home\carlo\mikeAnalysis\2P\';
bdata_source = 'A:\home\mike\Data\Behavior\';  
dateIdx=6;
elseif str2double(regexp(mouse,'\d*','match')) > 1700 || str2double(regexp(mouse,'\d*','match')) == 1081 
analysis_out = 'A:\home\carlo\analysis\2P\';
bdata_source = 'A:\home\carlo\rawData\behavioral\';
dateIdx=7;
end
%% call data Bx and Neural
    mworks = getBxData_Carlo(bdata_source, sessions, dateIdx);
cTargetOn = mworks.cTargetOn;
    if iscell(cTargetOn) % if it is a cell, it means cTargetOn wasn't being over written in extract TC. If it's not a cell, it's already over written in extract TC. can be used directly
        cTargetOn = celleqel2mat_padded(mworks.cTargetOn);
        cTargetOn(1) = nan; % first trial doesn't have reward 
    end
    
gTargetOn{1,mouseIdx} = cTargetOn; %concat without adding prev session value
gRctT = [gRctT celleqel2mat_padded(mworks.reactTimesMs)];
initFrames{1,mouseIdx} = mworks.counterValues{length(cTargetOn)}(end);
animalTrials(1,mouseIdx) = length(cTargetOn);
initLickTimes{1,mouseIdx} = (mworks.lickometerTimesUs); 
initCounterTimes{1,mouseIdx} = (mworks.counterTimesUs);
initCounterVals{1,mouseIdx} = (mworks.counterValues);

%% slight deviation to calculate all the trial data for subplot figures
block2{1,mouseIdx} = celleqel2mat_padded(mworks.tBlock2TrialNumber);
rew{1,mouseIdx} = double(~cell2mat(mworks.tBlock2TrialNumber));
if seshType == 2
    pmpt = input('CS-? [1-yes|2-no]: ');
    if pmpt == 1
        block2{1,mouseIdx} = rew{1,mouseIdx}; rew{1,mouseIdx} = zeros([1 length(rew{1,mouseIdx})]);
    else
    end
end
       
for t = 1:length(block2{1,mouseIdx})
    if t == 1 
block2OneRew{1,mouseIdx}(1,t) = 0;
block2TwoRew{1,mouseIdx}(1,t) = 0;
block2ThreeRew{1,mouseIdx}(1,t) = 0;
block2OneB2{1,mouseIdx}(1,t) = 0;
block2TwoB2{1,mouseIdx}(1,t) = 0;
block2ThreeB2{1,mouseIdx}(1,t) = 0;
rewOneRew{1,mouseIdx}(1,t) = 0;
rewTwoRew{1,mouseIdx}(1,t) = 0;
rewThreeRew{1,mouseIdx}(1,t) = 0;
rewOneB2{1,mouseIdx}(1,t) = 0;
rewTwoB2{1,mouseIdx}(1,t) = 0;
rewThreeB2{1,mouseIdx}(1,t) = 0;
    elseif t == 2
        if block2{1,mouseIdx}(1,t-1) == 1
block2OneRew{1,mouseIdx}(1,t) = 0;
block2OneB2{1,mouseIdx}(1,t) = 1;
rewOneRew{1,mouseIdx}(1,t) = 0;
rewOneB2{1,mouseIdx}(1,t) = 1;
        elseif block2{1,mouseIdx}(1,t-1) == 0
block2OneRew{1,mouseIdx}(1,t) = 1;
block2OneB2{1,mouseIdx}(1,t) = 0;
rewOneRew{1,mouseIdx}(1,t) = 1;
rewOneB2{1,mouseIdx}(1,t) = 0;            
        end
block2TwoRew{1,mouseIdx}(1,t) = 0;
block2TwoB2{1,mouseIdx}(1,t) = 0;
rewTwoRew{1,mouseIdx}(1,t) = 0;
rewTwoB2{1,mouseIdx}(1,t) = 0;
block2ThreeRew{1,mouseIdx}(1,t) = 0;
block2ThreeB2{1,mouseIdx}(1,t) = 0;
rewThreeRew{1,mouseIdx}(1,t) = 0;
rewThreeB2{1,mouseIdx}(1,t) = 0;
block2FourRew{1,mouseIdx}(1,t) = 0;
block2FourB2{1,mouseIdx}(1,t) = 0;
rewFourRew{1,mouseIdx}(1,t) = 0;
rewFourB2{1,mouseIdx}(1,t) = 0;
    elseif t >= 3
        % trial prior is CS- and trial before is not CS-
        if block2{1,mouseIdx}(1,t-1) == 1 && block2{1,mouseIdx}(1,t-2) ~= 1
block2OneRew{1,mouseIdx}(1,t) = 0;
block2OneB2{1,mouseIdx}(1,t) = 1;
rewOneRew{1,mouseIdx}(1,t) = 0;
rewOneB2{1,mouseIdx}(1,t) = 1; 
        % trial prior is CS+ and trial before is not CS+
        elseif block2{1,mouseIdx}(1,t-1) == 0 && block2{1,mouseIdx}(1,t-2) ~= 0
block2OneRew{1,mouseIdx}(1,t) = 1;
block2OneB2{1,mouseIdx}(1,t) = 0;
rewOneRew{1,mouseIdx}(1,t) = 1;
rewOneB2{1,mouseIdx}(1,t) = 0;
        elseif sum(block2{1,mouseIdx}(1,t-2:t-1)) == 2 || sum(block2{1,mouseIdx}(1,t-2:t-1)) == 0
block2OneRew{1,mouseIdx}(1,t) = 0;
block2OneB2{1,mouseIdx}(1,t) = 0;
rewOneRew{1,mouseIdx}(1,t) = 0;
rewOneB2{1,mouseIdx}(1,t) = 0;            
        end
        
        % trial prior and before is CS- and trial number is 3        
        if t == 3 && sum(block2{1,mouseIdx}(1,t-2:t-1)) == 2
block2TwoRew{1,mouseIdx}(1,t) = 0;
block2TwoB2{1,mouseIdx}(1,t) = 1;
rewTwoRew{1,mouseIdx}(1,t) = 0;
rewTwoB2{1,mouseIdx}(1,t) = 1;
        % trial prior and before is CS-, trial number is greater than 3, and the third before is not CS-
        elseif  t > 3 && sum(block2{1,mouseIdx}(1,t-2:t-1)) == 2 && block2{1,mouseIdx}(1,t-3) ~= 1
block2TwoRew{1,mouseIdx}(1,t) = 0;
block2TwoB2{1,mouseIdx}(1,t) = 1;
rewTwoRew{1,mouseIdx}(1,t) = 0;
rewTwoB2{1,mouseIdx}(1,t) = 1;        
        % one of the trial prior and before is CS+ and the other is CS-
        elseif sum(block2{1,mouseIdx}(1,t-2:t-1)) == 1
block2TwoRew{1,mouseIdx}(1,t) = 0;
block2TwoB2{1,mouseIdx}(1,t) = 0;
rewTwoRew{1,mouseIdx}(1,t) = 0;
rewTwoB2{1,mouseIdx}(1,t) = 0;   
        % trial prior and 2 before is CS+ and trial three before is not CS+
        elseif t > 3 && sum(block2{1,mouseIdx}(1,t-2:t-1)) == 0 && block2{1,mouseIdx}(1,t-3) ~= 0
block2TwoRew{1,mouseIdx}(1,t) = 1;
block2TwoB2{1,mouseIdx}(1,t) = 0;
rewTwoRew{1,mouseIdx}(1,t) = 1;
rewTwoB2{1,mouseIdx}(1,t) = 0;  
        elseif t > 3 && sum(block2{1,mouseIdx}(1,t-3:t-1)) == 3 || t > 3 && sum(block2{1,mouseIdx}(1,t-3:t-1)) == 0
block2TwoRew{1,mouseIdx}(1,t) = 0;
block2TwoB2{1,mouseIdx}(1,t) = 0;
rewTwoRew{1,mouseIdx}(1,t) = 0;
rewTwoB2{1,mouseIdx}(1,t) = 0; 
        elseif t == 3
block2FourB2{1,mouseIdx}(1,t) = 0;
block2FourRew{1,mouseIdx}(1,t) = 0;
block2ThreeB2{1,mouseIdx}(1,t) = 0; 
block2ThreeRew{1,mouseIdx}(1,t) = 0;
        end
    end
    if t == 4
        % trial prior and two before are CS-        
        if sum(block2{1,mouseIdx}(1,t-3:t-1)) == 3 
block2ThreeRew{1,mouseIdx}(1,t) = 0;
block2ThreeB2{1,mouseIdx}(1,t) = 1;
rewThreeRew{1,mouseIdx}(1,t) = 0;
rewThreeB2{1,mouseIdx}(1,t) = 1;
        % trial prior and two before are CS+
        elseif sum(block2{1,mouseIdx}(1,t-3:t-1)) == 0 
block2ThreeRew{1,mouseIdx}(1,t) = 1;
block2ThreeB2{1,mouseIdx}(1,t) = 0;
rewThreeRew{1,mouseIdx}(1,t) = 1;
rewThreeB2{1,mouseIdx}(1,t) = 0; 
        else%if sum(block2{1,i}(1,t-3:t-1)) ~= 3 || sum(block2{1,i}(1,t-3:t-1)) ~= 0
block2ThreeRew{1,mouseIdx}(1,t) = 0;
block2ThreeB2{1,mouseIdx}(1,t) = 0;
rewThreeRew{1,mouseIdx}(1,t) = 0;
rewThreeB2{1,mouseIdx}(1,t) = 0; 
        end
block2FourB2{1,mouseIdx}(1,t) = 0;
block2FourRew{1,mouseIdx}(1,t) = 0;
rewFourB2{1,mouseIdx}(1,t) = 0;
rewFourB2{1,mouseIdx}(1,t) = 0;        
    end
    
    if t > 4
        if sum(block2{1,mouseIdx}(1,t-3:t-1)) == 3 && block2{1,mouseIdx}(1,t-4) == 0
block2ThreeRew{1,mouseIdx}(1,t) = 0;
block2ThreeB2{1,mouseIdx}(1,t) = 1;
rewThreeRew{1,mouseIdx}(1,t) = 0;
rewThreeB2{1,mouseIdx}(1,t) = 1;
        elseif sum(block2{1,mouseIdx}(1,t-3:t-1)) == 0 && block2{1,mouseIdx}(1,t-4) == 1
block2ThreeRew{1,mouseIdx}(1,t) = 1;
block2ThreeB2{1,mouseIdx}(1,t) = 0;
rewThreeRew{1,mouseIdx}(1,t) = 1;
rewThreeB2{1,mouseIdx}(1,t) = 0;  
        % trial prior and two before are not CS-
        else%if sum(block2{1,i}(1,t-3:t-1)) ~= 3 || sum(block2{1,i}(1,t-3:t-1)) ~= 0
block2ThreeRew{1,mouseIdx}(1,t) = 0;
block2ThreeB2{1,mouseIdx}(1,t) = 0;
rewThreeRew{1,mouseIdx}(1,t) = 0;
rewThreeB2{1,mouseIdx}(1,t) = 0;                    
        end
        
        if sum(block2{1,mouseIdx}(1,t-4:t-1)) == 4
block2FourRew{1,mouseIdx}(1,t) = 0;
block2FourB2{1,mouseIdx}(1,t) = 1;
rewFourRew{1,mouseIdx}(1,t) = 0;
rewFourB2{1,mouseIdx}(1,t) = 1;   
        elseif sum(block2{1,mouseIdx}(1,t-4:t-1)) == 0
block2FourRew{1,mouseIdx}(1,t) = 1;
block2FourB2{1,mouseIdx}(1,t) = 0;
rewFourRew{1,mouseIdx}(1,t) = 1;
rewFourB2{1,mouseIdx}(1,t) = 0;   
        else%if sum(block2{1,i}(1,t-4:t-1)) ~= 4 || sum(block2{1,i}(t-4:t-1)) ~= 0
block2FourRew{1,mouseIdx}(1,t) = 0;
block2FourB2{1,mouseIdx}(1,t) = 0;
rewFourRew{1,mouseIdx}(1,t) = 0;
rewFourB2{1,mouseIdx}(1,t) = 0;       
        end
    end
    
    if block2{1,mouseIdx}(1,t) == 0
block2OneRew{1,mouseIdx}(1,t) = 0;
block2TwoRew{1,mouseIdx}(1,t) = 0;
block2ThreeRew{1,mouseIdx}(1,t) = 0;
block2FourRew{1,mouseIdx}(1,t) = 0;
block2OneB2{1,mouseIdx}(1,t) = 0;
block2TwoB2{1,mouseIdx}(1,t) = 0;
block2ThreeB2{1,mouseIdx}(1,t) = 0;
block2FourB2{1,mouseIdx}(1,t) = 0;
    elseif block2{1,mouseIdx}(1,t) == 1
rewOneRew{1,mouseIdx}(1,t) = 0;
rewTwoRew{1,mouseIdx}(1,t) = 0;
rewThreeRew{1,mouseIdx}(1,t) = 0;
rewFourRew{1,mouseIdx}(1,t) = 0;
rewOneB2{1,mouseIdx}(1,t) = 0;
rewTwoB2{1,mouseIdx}(1,t) = 0;
rewThreeB2{1,mouseIdx}(1,t) = 0;
rewFourB2{1,mouseIdx}(1,t) = 0;
    end
end

%remove 'sum' from the following to get the initial arrays
b2{1,mouseIdx} = sum(([block2OneB2{1,mouseIdx}; block2OneRew{1,mouseIdx}; block2TwoB2{1,mouseIdx}; block2TwoRew{1,mouseIdx}; block2ThreeB2{1,mouseIdx}; block2ThreeRew{1,mouseIdx}]),1);
r{1,mouseIdx} = sum(([rewOneB2{1,mouseIdx}; rewOneRew{1,mouseIdx}; rewTwoB2{1,mouseIdx}; rewTwoRew{1,mouseIdx}; rewThreeB2{1,mouseIdx}; rewThreeRew{1,mouseIdx}]),1);
allT{1,mouseIdx} = sum(sum([b2{1,mouseIdx}; r{1,mouseIdx}],1));

%% 
load([analysis_out img_fn 'figDataUntrimmed_' sessions '.mat'],'-regexp','^(?!expt|img|ind).'); 

groupAlign.tc{1,mouseIdx} = targetAlign_tc;
groupAlign.events{1,mouseIdx} = targetAlign_events;
gtt = [gtt; tt];

%ggoodcells = [ggoodcells goodcells];
clearvars -except expt rew* seshType events loadData tc isesh exptCount piezo mouseIdx lick block2* mworks early* late* ind* iexp animalTrials init* tt target* g* date* day sessions mouse run i img_fn analysis_out bdata_source output_fn pairings

load([analysis_out img_fn '_cueAlignLick.mat'],'-regexp','^(?!expt|img|ind|early|late).');
  
groupAlign.licks{1,mouseIdx} = lickCueAlign; 
groupAlign.lickStart{1,mouseIdx} = lickBurstStart;
  if sum(~isnan(groupAlign.lickStart{1,mouseIdx}))>6
        [sortlick sortlick_ind] = sort(groupAlign.lickStart{1,mouseIdx},'ascend');
        nburst = sum(~isnan(groupAlign.lickStart{1,mouseIdx}));
        nnan = sum(isnan(groupAlign.lickStart{1,mouseIdx}));
%earliest and latest 25% lick - count uneven due to behavioral differences
        early_bst{1,mouseIdx} = sortlick_ind(1:floor(nburst/4)); %trials with the earliest 25% of licks
        late_bst{1,mouseIdx} = sortlick_ind(nburst-floor(nburst/4)+1:end-nnan); %trials with the latest 25% of licks
        early_time{1,mouseIdx} = mean((groupAlign.lickStart{1,mouseIdx}(:,early_bst{1,mouseIdx})-prewin_frames).*(1000./frameRateHz),2);
        late_time{1,mouseIdx} = mean((groupAlign.lickStart{1,mouseIdx}(:,late_bst{1,mouseIdx})-prewin_frames).*(1000./frameRateHz),2);
  else
        early_bst{1,mouseIdx} = [];
        late_bst{1,mouseIdx} = [];
        early_time{1,mouseIdx} = [];
        late_time{1,mouseIdx} = [];
  end
    ind_early_bst = [ind_early_bst early_bst{1,mouseIdx}];
    ind_late_bst = [ind_late_bst late_bst{1,mouseIdx}];
    early_bsttime = [early_bsttime early_time{1,mouseIdx}];
    late_bsttime = [late_bsttime late_time{1,mouseIdx}];
    
earlyBurstEvents = [earlyBurstEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,early_bst{1,mouseIdx}),3)];
lateBurstEvents = [lateBurstEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,late_bst{1,mouseIdx}),3)];
earlyBurstStart = [earlyBurstStart nanmean(groupAlign.lickStart{1,mouseIdx}(:,early_bst{1,mouseIdx}),3)];
lateBurstStart = [lateBurstStart nanmean(groupAlign.lickStart{1,mouseIdx}(:,late_bst{1,mouseIdx}),3)];
clearvars -except expt rew* seshType events tc loadData isesh exptCount piezo lick mouseIdx block2* mworks iexp ind* early* late* animalTrials init* target* g* date* day sessions mouse run i img_fn analysis_out bdata_source output_fn pairings
 

  load([analysis_out img_fn '_cueAlignPiezo.mat'],'-regexp','^(?!expt|img|ind).');
  
groupAlign.piezo{1,mouseIdx} = targetAlign_piezo;
clearvars -except expt rew* seshType events tc isesh loadData exptCount piezo lick mouseIdx block2* iexp ind* early* late* animalTrials init* target* g* Goodcells date* day sessions mouse run i img_fn analysis_out bdata_source output_fn pairings


tc.rewGroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rew{1,mouseIdx})),3);
tc.b2GroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2{1,mouseIdx})),3);
events.rewGroupAlign{1,mouseIdx} = nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(rew{1,mouseIdx})),3);
events.b2GroupAlign{1,mouseIdx} = nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(block2{1,mouseIdx})),3);
piezo.rewGroupAlign{1,mouseIdx} = groupAlign.piezo{1,mouseIdx}(:,logical(rew{1,mouseIdx}));
piezo.b2GroupAlign{1,mouseIdx} = groupAlign.piezo{1,mouseIdx}(:,logical(block2{1,mouseIdx}));
lick.rewGroupAlign{1,mouseIdx} = groupAlign.licks{1,mouseIdx}(:,logical(rew{1,mouseIdx}));
lick.b2GroupAlign{1,mouseIdx} = groupAlign.licks{1,mouseIdx}(:,logical(block2{1,mouseIdx}));

tc.rewOneRewGroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rewOneRew{1,mouseIdx})),3);
events.rewOneRewGroupAlign{1,mouseIdx} =  nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(rewOneRew{1,mouseIdx})),3);
piezo.rewOneRewGroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(rewOneRew{1,mouseIdx}));
lick.rewOneRewGroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(rewOneRew{1,mouseIdx}));
tc.b2OneRewGroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2OneRew{1,mouseIdx})),3);
events.b2OneRewGroupAlign{1,mouseIdx} =  nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(block2OneRew{1,mouseIdx})),3);
piezo.b2OneRewGroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(block2OneRew{1,mouseIdx}));
lick.b2OneRewGroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(block2OneRew{1,mouseIdx}));

tc.rewTwoRewGroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rewTwoRew{1,mouseIdx})),3);
events.rewTwoRewGroupAlign{1,mouseIdx} =  nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(rewTwoRew{1,mouseIdx})),3);
piezo.rewTwoRewGroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(rewTwoRew{1,mouseIdx}));
lick.rewTwoRewGroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(rewTwoRew{1,mouseIdx}));
tc.b2TwoRewGroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2TwoRew{1,mouseIdx})),3);
events.b2TwoRewGroupAlign{1,mouseIdx} =  nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(block2TwoRew{1,mouseIdx})),3);
piezo.b2TwoRewGroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(block2TwoRew{1,mouseIdx}));
lick.b2TwoRewGroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(block2TwoRew{1,mouseIdx}));

tc.rewThreeRewGroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rewThreeRew{1,mouseIdx})),3);
events.rewThreeRewGroupAlign{1,mouseIdx} =  nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(rewThreeRew{1,mouseIdx})),3);
piezo.rewThreeRewGroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(rewThreeRew{1,mouseIdx}));
lick.rewThreeRewGroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(rewThreeRew{1,mouseIdx}));
tc.b2ThreeRewGroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2ThreeRew{1,mouseIdx})),3);
events.b2ThreeRewGroupAlign{1,mouseIdx} =  nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(block2ThreeRew{1,mouseIdx})),3);
piezo.b2ThreeRewGroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(block2ThreeRew{1,mouseIdx}));
lick.b2ThreeRewGroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(block2ThreeRew{1,mouseIdx}));

tc.rewOneB2GroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rewOneB2{1,mouseIdx})),3);
events.rewOneB2GroupAlign{1,mouseIdx} =  nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(rewOneB2{1,mouseIdx})),3);
piezo.rewOneB2GroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(rewOneB2{1,mouseIdx}));
lick.rewOneB2GroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(rewOneB2{1,mouseIdx}));
tc.b2OneB2GroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2OneB2{1,mouseIdx})),3);
events.b2OneB2GroupAlign{1,mouseIdx} =  nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(block2OneB2{1,mouseIdx})),3);
piezo.b2OneB2GroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(block2OneB2{1,mouseIdx}));
lick.b2OneB2GroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(block2OneB2{1,mouseIdx}));

tc.rewTwoB2GroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rewTwoB2{1,mouseIdx})),3);
events.rewTwoB2GroupAlign{1,mouseIdx} =  nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(rewTwoB2{1,mouseIdx})),3);
piezo.rewTwoB2GroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(rewTwoB2{1,mouseIdx}));
lick.rewTwoB2GroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(rewTwoB2{1,mouseIdx}));
tc.b2TwoB2GroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2TwoB2{1,mouseIdx})),3);
events.b2TwoB2GroupAlign{1,mouseIdx} =  nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(block2TwoB2{1,mouseIdx})),3);
piezo.b2TwoB2GroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(block2TwoB2{1,mouseIdx}));
lick.b2TwoB2GroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(block2TwoB2{1,mouseIdx}));

tc.rewThreeB2GroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rewThreeB2{1,mouseIdx})),3);
events.rewThreeB2GroupAlign{1,mouseIdx} =  nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(rewThreeB2{1,mouseIdx})),3);
piezo.rewThreeB2GroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(rewThreeB2{1,mouseIdx}));
lick.rewThreeB2GroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(rewThreeB2{1,mouseIdx}));
tc.b2ThreeB2GroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2ThreeB2{1,mouseIdx})),3);
events.b2ThreeB2GroupAlign{1,mouseIdx} =  nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(block2ThreeB2{1,mouseIdx})),3);
piezo.b2ThreeB2GroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(block2ThreeB2{1,mouseIdx}));
lick.b2ThreeB2GroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(block2ThreeB2{1,mouseIdx}));

tc.rewFourB2GroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rewFourB2{1,mouseIdx})),3);
events.rewFourB2GroupAlign{1,mouseIdx} =  nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(rewFourB2{1,mouseIdx})),3);
piezo.rewFourB2GroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(rewFourB2{1,mouseIdx}));
lick.rewFourB2GroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(rewFourB2{1,mouseIdx}));
tc.b2FourB2GroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2FourB2{1,mouseIdx})),3);
events.b2FourB2GroupAlign{1,mouseIdx} =  nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(block2FourB2{1,mouseIdx})),3);
piezo.b2FourB2GroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(block2FourB2{1,mouseIdx}));
lick.b2FourB2GroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(block2FourB2{1,mouseIdx}));

tc.rewFourRewGroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rewFourRew{1,mouseIdx})),3);
events.rewFourRewGroupAlign{1,mouseIdx} =  nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(rewFourRew{1,mouseIdx})),3);
piezo.rewFourRewGroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(rewFourRew{1,mouseIdx}));
lick.rewFourRewGroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(rewFourRew{1,mouseIdx}));
tc.b2FourRewGroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2FourRew{1,mouseIdx})),3);
events.b2FourRewGroupAlign{1,mouseIdx} =  nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(block2FourRew{1,mouseIdx})),3);
piezo.b2FourRewGroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(block2FourRew{1,mouseIdx}));
lick.b2FourRewGroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(block2FourRew{1,mouseIdx}));

end
end
%
if seshType == 2 && loadData ~= 2
    blockName = input('which block is this? [type "plus/minus""First/Second"]: ','s');
end
%
if seshType ~= 2
save(fullfile(output_fn, '_groupwiseDataset.mat'));
else
    save(fullfile(output_fn, [blockName '_groupwiseDataset.mat']));
end

groupAlign_tc=[]; groupAlign_events=[]; groupAlign_piezo=[]; groupAlign_lick=[];
for a = 1:exptCount
nsize(1,a) = size(groupAlign.events{1,a},2);
end



for mouseIdx = 1:exptCount
  %rDend(i,:) = randperm(size(groupAlign.events{1,i},2),nsize);

  groupAlign_tc = cat(2,groupAlign_tc, nanmean(groupAlign.tc{1,mouseIdx},3));%(:,rDend(i,:)
  groupAlign_events = cat(2,groupAlign_events, nanmean(groupAlign.events{1,mouseIdx},3));
  groupAlign_piezo = [groupAlign_piezo groupAlign.piezo{1,mouseIdx}];
  groupAlign_lick = [groupAlign_lick groupAlign.licks{1,mouseIdx}];
    
end
%need to pull the same nNeurons from each mouse


%CHECK:: ggoodcells should equal groupAlign_tc(:,x) and
%   groupAlign_events(:,x)

cRctT = nanmean(gRctT);

frameRateHz = 30;
prewin_frames = round(1500./frameRateHz); 
postwin_frames = round(3000./frameRateHz);

cue = input('which session is this [CS+|CS-|Interleaved]?: '); % 1: CS+ | 2: CS- | 3: interleaved
if cue == 1
    figTitle = 'CS+';
elseif cue == 2
    figTitle = 'CS-';
elseif cue == 3
    figTitle = 'Interleaved';
end
%% lick align formatting

  nTrials = length(cell2mat(gTargetOn));
    rewDelay_frames = round((cRctT/1000).*frameRateHz);% there's 700ms between the cue and the reward delivery !!!!! if you change tooFastMs in the varibles in MWorks, this needs to be changed
    nIC = size(groupAlign_events,2); %you get targetalign_events from neural align, this is the neural data (in spikes) aligned to cue
    lickCueAlign =  nan(prewin_frames+postwin_frames,nTrials);
    lickCounterVals = cell(1,nTrials);
    for a = 1:exptCount
    nFrames{1,a} = initCounterVals{1,a}{length(gTargetOn{1,a})}(end);
    end
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
    preRewTrials = zeros(1,nTrials);
    postRewTrials = zeros(1,nTrials);
    postRew_lickSearchRange = prewin_frames+rewDelay_frames+lickDelay_frames:size(lickCueAlign,1)-lickSearch_frames-postLick_frames;% need to '-lickSearch_frames-postLick_frames' because later the inds needs to + lickSearch_frames or +postLick_frames, this is just for the following inds to be within matrix dimentions
    lastPreRewLickFrame = nan(1,nTrials);
    firstPostRewLickFrame = nan(1,nTrials);
    
    tTInx = 1; ntrials=0;
    for a=1:exptCount%lick align organization 
oDend = (1:nsize(1,a))+sum(nsize(1,1:a-1));
nDend = 1:nsize(1,a);       

    for itrial = 1:animalTrials(1,a) %nTrials%
        ind_post = [];
        ind_pre = [];
        ntrial=itrial+ntrials;
        if ~isnan(gTargetOn{1,a}(itrial))
            if gTargetOn{1,a}(itrial)+postwin_frames-1 < nFrames{1,a}
                lickTimes = initLickTimes{1,a}{itrial}; 
                counterTimes = initCounterTimes{1,a}{itrial};
                counterVals = initCounterVals{1,a}{itrial};
                lickCounterVals{ntrial} = zeros(size(lickTimes));
                lickTC{ntrial} = zeros(size(counterVals));% find how many licks for each frame
             %%%this for loop pulls out every lickTime that falls between the start of two frames
               %and if licks are found it outputs the frame licks occur to lickCounterVals
                for icount = 1:length(counterTimes)-1
                    ind = find(lickTimes>counterTimes(icount) & lickTimes<counterTimes(icount+1));
                    if ~isempty(ind)
                        lickCounterVals{ntrial}(1,ind) = icount; % find which counter in this trial has licks
                    end
                end
             %%%this for loop counts down the number of frames in a trial and outputs a logical
               %outlining whether or not a lick occured between two counters (or frames)
                for ival = 1:length(counterTimes)
                    ind = find(lickCounterVals{ntrial} == ival);
                    if ~isempty(ind)
                        lickTC{ntrial}(1,ival) = length(ind);% if the mice doesn't lick during that counter than it's zero, if it does, it's one. length(ind) should always be 1.
                    end
                end
             %%%IF the last frame of the trial (from the start of the ITI preceeding the trial to the start of the ITI following the trial)
               %minus the frame where the cue is shown IS greater than the frame difference between   
                if initCounterVals{1,a}{itrial}(end)-gTargetOn{1,a}(itrial) > postwin_frames
                    %lickcuealign aligns licking of each frame to cue
                    lickCueAlign(:,ntrial) = lickTC{ntrial}(1,gTargetOn{1,a}(itrial)-prewin_frames-counterVals(1):gTargetOn{1,a}(itrial)+postwin_frames-1-counterVals(1));
                end
                ind = intersect(lickSearchRange,find(lickCueAlign(:,ntrial)));
             %%%this for loop find frames where more than 3 licks occurred within a fixed window, and index that frame as the start of burst lick
                for mouseIdx = 1:length(ind)
                    ilick = ind(mouseIdx);
                    if sum(lickCueAlign(ilick:ilick+lickSearch_frames-1,ntrial),1) >= 3 % more than 3 bursts within about 300ms
                        lickBurstStart(:,ntrial) = ilick;% when licking burst happens
                        break
                    end
                end
                ind_pre = intersect(preRew_lickSearchRange_700ms,find(lickCueAlign(:,ntrial))); %finds every instance of a lick that occurs within the search window
                ind_post = intersect(postRew_lickSearchRange,find(lickCueAlign(:,ntrial))); %finds every instance of a lick that occurs [after the reward delivery]
                ind_all = intersect(lickSearchRange,find(lickCueAlign(:,ntrial)));
                preRew_lickBurstHz(:,ntrial) = length(ind_pre)./(cRctT/1000);
                postRew_lickBurstHz(:,ntrial) = length(ind_post)./(length(postRew_lickSearchRange)./frameRateHz);
                lickBurstHz_all(:,ntrial) = length(ind_all)./(length(lickSearchRange)./frameRateHz);
                ind = intersect(postRew_lickSearchRange,find(lickCueAlign(:,ntrial)));
             %%%similar idea as the above for loop, but instead aligns neural data as POST-REWARD lick events 
                for mouseIdx = 1:length(ind)
                    ilick = ind(mouseIdx);
                    if sum(lickCueAlign(ilick:ilick+lickSearch_frames-1,ntrial),1) >= 3
                        postRew_lickAlignEvents(:,oDend,ntrial) = groupAlign.events{1,a}(ilick-postLick_frames:ilick+postLick_frames-1,nDend,itrial);% align neural data to first lick after reward
                        postRew_lickBurstStart(:,ntrial) = ilick; %array of every lick within post-reward window if 3+ licks were recorded for that trial
                        postRew_lickAlign(:,ntrial) = lickCueAlign(ilick-postLick_frames:ilick+postLick_frames-1,ntrial);% align licking data to first lick
                        break
                    end
                end
                ind_pre = find(lickCueAlign(prewin_frames:prewin_frames+rewDelay_frames,ntrial),1,'last');
             %%%if there are licks recorded between cue and reward delivery (770 ms), align neural data to those licks 
                if ~isempty(ind_pre)
                    lastPreRewLickFrame(1,ntrial) = ind_pre;
                    preRewTrials(1,ntrial) =  1;
                    % align neural and licking data to the last lick
                    lastPreRew_lickAlignEvents(:,oDend,ntrial) = groupAlign.events{1,a}(ind_pre-postLick_frames+prewin_frames:ind_pre+postLick_frames+postLick_frames-1+prewin_frames,nDend,itrial);
                    lastPreRew_lickAlign(:,ntrial) = lickCueAlign(ind_pre-postLick_frames+prewin_frames:ind_pre+postLick_frames+postLick_frames-1+prewin_frames,ntrial);
                end
                ind_post = find(lickCueAlign(prewin_frames+rewDelay_frames:prewin_frames+postwin_frames-postLick_frames-postLick_frames,ntrial),1,'first');
             %%%   
                if ~isempty(ind_post)
                    firstPostRewLickFrame(1,ntrial) = ind_post;
                    postRewTrials(1,ntrial) = 1;
                    % align neural and licking data to the first lick
                    firstPostRew_lickAlignEvents(:,oDend,ntrial) = groupAlign.events{1,a}(ind_post-postLick_frames+prewin_frames+rewDelay_frames:ind_post+postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,nDend,itrial);
                    firstPostRew_lickAlign(:,ntrial) = lickCueAlign(ind_post-postLick_frames+prewin_frames+rewDelay_frames:ind_post+postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,ntrial);
                end
            rewAlignEvents(:,oDend,ntrial) = groupAlign.events{1,a}(-postLick_frames+prewin_frames+rewDelay_frames:postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,nDend,itrial);

            tTInx = tTInx+1;
            end
        end
    end
    ntrials=ntrial;
    end 
    
%% piezo align formatting

ind_nan = find(isnan(groupAlign_piezo(1,:)));

g= [0.4660 0.6740 0.1880];
r= [0.6350 0.0780 0.1840];


%%

%% lick align figures
tt=gtt(1,:);

% lickBurstStart(1,241:end) has no lick bursts recorded - double check to
% make sure thats actually the case and there isn't a coding error
  tempFig=figure; %align licking data to cue
    shadedErrorBar(tt, nanmean(lickCueAlign,2).*(1000./frameRateHz), (nanstd(lickCueAlign,[],2)./sqrt(unique(sum(~isnan(lickCueAlign),2))).*(1000./frameRateHz)));
    hold on;
    scatter((lickBurstStart-prewin_frames).*(1000./frameRateHz), 10.*ones(size(lickBurstStart)), 'x');
    xlabel('Time from cue');
    ylabel('Lick rate (Hz)');
    title(['lick events from ' num2str(size(expt,2)) ' animals | cue aligned licking']);
    vline(0,'k');
    vline(767,'k');
    if seshType ~= 2
    savefig(fullfile(output_fn, '_cueAlign_lickHz.fig'));
    saveas(tempFig, [output_fn  '_cueAlign_lickHz.jpeg']);
    else
    savefig(fullfile(output_fn, [blockName '_cueAlign_lickHz.fig']));
    saveas(tempFig, [output_fn blockName '_cueAlign_lickHz.jpeg']);
    end
    
    nIC = size(groupAlign_events,2);

    tempFig=figure; %plot neural data of trials of early vs. late bursts
    shadedErrorBar(tt,nanmean(nanmean(earlyBurstEvents,3),2).*(1000./frameRateHz), (nanstd(nanmean(earlyBurstEvents,3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'k');
    hold on;
    shadedErrorBar(tt,nanmean(nanmean(lateBurstEvents,3),2).*(1000./frameRateHz), (nanstd(nanmean(lateBurstEvents,3),[],2).*(1000./frameRateHz))./sqrt(nIC),'b');
    scatter((earlyBurstStart-prewin_frames).*(1000./frameRateHz),0.*ones(1,length(earlyBurstStart)),'xk');
    scatter((lateBurstStart-prewin_frames).*(1000./frameRateHz),0.*ones(1,length(lateBurstStart)),'xb');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    title([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | lick bursts: early (n = ' num2str(length(cell2mat(early_bst))) '| blk); late (n = ' num2str(length(cell2mat(late_bst))) '| blu)']);
    if seshType ~=2
    savefig(fullfile(output_fn, '_cueAlignSpiking_byLickTime.fig'));
    saveas(tempFig, [output_fn  '_cueAlignSpiking_byLickTime.jpeg']);
    else
    savefig(fullfile(output_fn, ['\' blockName '_cueAlignSpiking_byLickTime.fig']));
    saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byLickTime.jpeg']);
    end
    hold off;

    pct_precue_burst = length(find((lickBurstStart-prewin_frames).*(1000./frameRateHz)<600))./size(lickBurstStart,2);
    tl = (1-postLick_frames:postLick_frames).*(1000./frameRateHz);
   
    ntrials=0; ndends=0; earlyBurstAlignEvents=[]; lateBurstAlignEvents=[]; earlyPostLickAlign=[]; latePostLickAlign=[];
    for mouseIdx = 1:exptCount
        nDend = (1:nsize(1,mouseIdx))+ndends;
        nTrial = (1:animalTrials(1,mouseIdx))+ntrials; 
        perAnimal_postRewLickEvents{1,mouseIdx}=postRew_lickAlignEvents(:,nDend,nTrial);
        ntrials = max(nTrial);
        ndends = max(nDend);
        
    perAnimal_postRewBurstStart{1,mouseIdx} = postRew_lickBurstStart(:,nTrial);
    perAnimal_postRewLickAlign{1,mouseIdx} = postRew_lickAlign(:,nTrial);
    [sortlick sortlick_ind] = sort(perAnimal_postRewBurstStart{1,mouseIdx},'ascend'); %sorts trial-based instances of lick and outputs [sorted instances, index of original position before sorting]
    nburst{1,mouseIdx} = sum(~isnan(perAnimal_postRewBurstStart{1,mouseIdx})); %total trials with lick burst
    nnan = sum(isnan(perAnimal_postRewBurstStart{1,mouseIdx})); %total trials without lick burst (some may be CS-; some may have no burst but still be CS+)
    ind_prerew_early_bst{1,mouseIdx} = sortlick_ind(1:floor(nburst{1,mouseIdx}/4)); %first 1/4 chosen as early trials
    ind_prerew_late_bst{1,mouseIdx} = sortlick_ind(nburst{1,mouseIdx}-floor(nburst{1,mouseIdx}/4)+1:end-nnan); %last 1/4 chosen as late trials 
    early_bst_time{1,mouseIdx} = nanmean((perAnimal_postRewBurstStart{1,mouseIdx}(:,ind_prerew_early_bst{1,mouseIdx})-prewin_frames-rewDelay_frames).*(1000./frameRateHz),2);
    late_bst_time{1,mouseIdx} = nanmean((perAnimal_postRewBurstStart{1,mouseIdx}(:,ind_prerew_late_bst{1,mouseIdx})-prewin_frames-rewDelay_frames).*(1000./frameRateHz),2);
    
    earlyPostLickAlign = [earlyPostLickAlign perAnimal_postRewLickAlign{1,mouseIdx}(:,ind_prerew_early_bst{1,mouseIdx})];
    latePostLickAlign = [latePostLickAlign perAnimal_postRewLickAlign{1,mouseIdx}(:,ind_prerew_late_bst{1,mouseIdx})];
    earlyBurstAlignEvents = [earlyBurstAlignEvents nanmean(perAnimal_postRewLickEvents{1,mouseIdx}(:,:,ind_prerew_early_bst{1,mouseIdx}),3)];
    lateBurstAlignEvents = [lateBurstAlignEvents nanmean(perAnimal_postRewLickEvents{1,mouseIdx}(:,:,ind_prerew_late_bst{1,mouseIdx}),3)];
    end
    
    tempFig=figure;
    subplot(2,1,1); % align neural activity to lick onset, only burst trials are included 
    shadedErrorBar(tl, nanmean(nanmean(postRew_lickAlignEvents,3),2).*(1000./frameRateHz), (nanstd(nanmean(postRew_lickAlignEvents,3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    hold on;
    xlabel('Time from lick');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    title([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | ' num2str(sum(~isnan(postRew_lickBurstStart(1,:)))) ' trials post-reward lick bursts']);
    subplot(2,1,2); % seperate early burst trials vs. late burst trials
    shadedErrorBar(tl, nanmean(earlyBurstAlignEvents,2).*(1000./frameRateHz), (nanstd(earlyBurstAlignEvents,[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    hold on;
    shadedErrorBar(tl, nanmean(lateBurstAlignEvents,2).*(1000./frameRateHz), (nanstd(lateBurstAlignEvents,[],2).*(1000./frameRateHz))./sqrt(nIC),'b');
    xlabel('Time from lick');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    title([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | ' num2str(floor(sum(cell2mat(nburst))/4)) ' early [blk]: avg = ' num2str(nanmean(cell2mat(early_bst_time))) ' ms; ' num2str(floor(sum(cell2mat(nburst))/4)) ' late [blu]: avg = ' num2str(nanmean(cell2mat(late_bst_time))) ' ms'])
    if seshType ~= 2
    savefig(fullfile(output_fn, '_postRew_lickBurstAlignSpikingEvents.fig'))
    saveas(tempFig, [output_fn  '_postRew_lickBurstAlignSpikingEvents.jpeg']);
    else
    savefig(fullfile(output_fn, [blockName '_postRew_lickBurstAlignSpikingEvents.fig']));
    saveas(tempFig, [output_fn blockName '_postRew_lickBurstAlignSpikingEvents.jpeg']);
    end
    
    tempFig=figure;
    subplot(2,1,1); % align licking data to first lick
    shadedErrorBar(tl, nanmean(postRew_lickAlign,2).*(1000./frameRateHz), (nanstd(postRew_lickAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(postRew_lickAlign),2))),'k');
    hold on;
    xlabel('Time from lick');
    ylabel('Lick rate (Hz)');
    ylim([0 inf]);
    title('Neural align first lick post-reward')
    subplot(2,1,2); % plot early burst trials and late burst trials separately
    shadedErrorBar(tl, nanmean(earlyPostLickAlign,2).*(1000./frameRateHz), (nanstd(earlyPostLickAlign,[],2).*(1000./frameRateHz))./sqrt(length(cell2mat(ind_prerew_early_bst))),'k');
    hold on;
    shadedErrorBar(tl, nanmean(latePostLickAlign,2).*(1000./frameRateHz), (nanstd(latePostLickAlign,[],2).*(1000./frameRateHz))./sqrt(length(cell2mat(ind_prerew_late_bst))),'b');
    xlabel('Time from lick');
    ylabel('Lick rate (Hz)');
    ylim([0 inf]);
    title('Early burst trials [blk] and late burst trials [blu]')
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | post reward lick burst aligned spiking']);
    if seshType ~= 2
    savefig(fullfile(output_fn, '_postRew_lickBurstAlignSpiking.fig'));
    saveas(tempFig, [output_fn  '_postRew_lickBurstAlignSpiking.jpeg']);
    else
    savefig(fullfile(output_fn, [blockName '_postRew_lickBurstAlignSpiking.fig']));
    saveas(tempFig, [output_fn blockName '_postRew_lickBurstAlignSpiking.jpeg']);
    end
    hold off;
    
%     tempFig=figure;
%     colour = {'k', 'b', 'r', 'm'};
%     for i = 1:4
%         plot(tt,(cumsum(nansum(lickCueAlign(:,1+((i-1)*floor(nTrials/4)):floor(nTrials/4)+((i-1)*floor(nTrials/4))),2))),colour{1,i}, 'LineWidth',2.0);
%         hold on;
%     end
%     title([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | cumulative licking by quarter session']);
%     xlabel('Time from cue');
%     ylabel('Cumulative Licks');
%     legend({'first quar', 'second quar', 'third quar', 'fourth quar'},'Location','northwest');
%     hold off;
%     savefig(fullfile(output_fn, '_cumulativeLicking.fig'));

ndends=0; ntrials=0; lastPreLickEvent=[]; firstPostLickEvent=[];               
for mouseIdx = 1:exptCount
    nDend = (1:nsize(1,mouseIdx))+ndends;
    nTrial = (1:animalTrials(1,mouseIdx))+ntrials; 
               lastPreLickEvent = [lastPreLickEvent (nanmean(lastPreRew_lickAlignEvents(:,nDend,nTrial),3))];
               firstPostLickEvent = [firstPostLickEvent (nanmean(firstPostRew_lickAlignEvents(:,nDend,nTrial),3))];
    ntrials = max(nTrial);
    ndends = max(nDend);
end
                
    tl_rew = (1-postLick_frames:(postLick_frames*2)).*(1000./frameRateHz);
    tempFig=figure; % align neural and licking data to the first and last lick, respectively
    subplot(2,2,1);
    shadedErrorBar(tl_rew, nanmean(lastPreLickEvent,2).*(1000./frameRateHz), (nanstd(lastPreLickEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    hold on;
    title('Last lick before reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 2]);
    subplot(2,2,3);
    shadedErrorBar(tl_rew, nanmean(lastPreRew_lickAlign,2).*(1000./frameRateHz), (nanstd(lastPreRew_lickAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(lastPreRew_lickAlign)))),'k');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 35]);
    title([num2str(sum(preRewTrials)) 'trials with pre-reward licks']);
    subplot(2,2,2);
    shadedErrorBar(tl_rew, nanmean(firstPostLickEvent,2).*(1000./frameRateHz), (nanstd(firstPostLickEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    title('First lick after reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 2]);
    subplot(2,2,4);
    shadedErrorBar(tl_rew, nanmean(firstPostRew_lickAlign,2).*(1000./frameRateHz), (nanstd(firstPostRew_lickAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(firstPostRew_lickAlign)))),'k');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 35]);
    title([num2str(sum(postRewTrials)) 'trials with post-reward licks']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Licks relative to reward']);
    hold off;
    if seshType ~= 2
    savefig(fullfile(output_fn, '_lastVsFirstLick.fig'));
    saveas(tempFig, [output_fn  '_lastVsFirstLick.jpeg']);
    else
    savefig(fullfile(output_fn, [blockName '_lastVsFirstLick.fig']));
    saveas(tempFig, [output_fn blockName '_lastVsFirstLick.jpeg']);
    end
    
    ndends=0; ntrials=0; lowAllEvents=[]; highAllEvents=[]; lowPreEvents=[]; highPreEvents=[]; lowPostEvents=[]; highPostEvents=[];
    for mouseIdx = 1:exptCount
        nDend = (1:nsize(1,mouseIdx))+ndends;
        nTrial = (1:animalTrials(1,mouseIdx))+ntrials; 
       
        tempLickBurstHz = lickBurstHz_all(1,nTrial);
    [sortlickHz sortlickHz_ind] = sort(tempLickBurstHz,'ascend');
    nburst = sum(~isnan(tempLickBurstHz));
    nnan = sum(isnan(tempLickBurstHz));
    ind_low_bst = sortlickHz_ind(1:floor(nburst/4));
    ind_high_bst = sortlickHz_ind(nburst-floor(nburst/4)+1:end-nnan);
    HL_lickrate.low_rew{1,mouseIdx} = mean(tempLickBurstHz(:,ind_low_bst),2);
    HL_lickrate.high_rew{1,mouseIdx} = mean(tempLickBurstHz(:,ind_high_bst),2);
    
        tempLickBurstHz = preRew_lickBurstHz(1,nTrial);    
    [sortlickHz sortlickHz_ind] = sort(tempLickBurstHz,'ascend');
    nnan = sum(isnan(tempLickBurstHz));
    nburst = sum(~isnan(tempLickBurstHz));
    ind_low_prerew = sortlickHz_ind(1:floor(nburst/4));
    ind_high_prerew = sortlickHz_ind(nburst-floor(nburst/4)+1:end-nnan);
    HL_lickrate.low_prerew{1,mouseIdx} = mean(tempLickBurstHz(:,ind_low_prerew),2);
    HL_lickrate.high_prerew{1,mouseIdx} = mean(tempLickBurstHz(:,ind_high_prerew),2);

        tempLickBurstHz = postRew_lickBurstHz(1,nTrial);    
    [sortlickHz sortlickHz_ind] = sort(tempLickBurstHz,'ascend');
    nburst = sum(~isnan(tempLickBurstHz));
    nnan = sum(isnan(tempLickBurstHz));
    ind_low_postrew = sortlickHz_ind(1:floor(nburst/4));
    ind_high_postrew = sortlickHz_ind(nburst-floor(nburst/4)+1:end-nnan);
    HL_lickrate.low_postrew{1,mouseIdx} = mean(tempLickBurstHz(:,ind_low_postrew),2);
    HL_lickrate.high_postrew{1,mouseIdx} = mean(tempLickBurstHz(:,ind_high_postrew),2);

    lowAllEvents = [lowAllEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_low_bst),3)];
    highAllEvents = [highAllEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_high_bst),3)];
    lowPreEvents = [lowPreEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_low_prerew),3)];
    highPreEvents = [highPreEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_high_prerew),3)];
    lowPostEvents = [lowPostEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_low_postrew),3)];
    highPostEvents = [highPostEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_high_postrew),3)];

    ndends = max(nDend);
    ntrials = max(nTrial);
    end
    
    tempFig=figure; % still plotting neural data, but seperate the trials based on licking rate of that trial
    subplot(3,1,1);
    shadedErrorBar(tt,nanmean(lowAllEvents,2).*(1000./frameRateHz), (nanstd(lowAllEvents,[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
    hold on;
    shadedErrorBar(tt,nanmean(highAllEvents,2).*(1000./frameRateHz), (nanstd(highAllEvents,[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    title(['1s- Rew: ' num2str(chop(nanmean(cell2mat(HL_lickrate.low_rew)),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_lickrate.high_rew)),2)) ' Hz']);
    subplot(3,1,2);
    shadedErrorBar(tt,nanmean(lowPreEvents,2).*(1000./frameRateHz), (nanstd(lowPreEvents,[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
    hold on;
    shadedErrorBar(tt,nanmean(highPreEvents,2).*(1000./frameRateHz), (nanstd(highPreEvents,[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    title(['Pre- Rew: ' num2str(chop(nanmean(cell2mat(HL_lickrate.low_prerew)),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_lickrate.high_prerew)),2)) ' Hz']);
    subplot(3,1,3);
    shadedErrorBar(tt,nanmean(lowPostEvents,2).*(1000./frameRateHz), (nanstd(lowPostEvents,[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
    hold on;
    shadedErrorBar(tt,nanmean(highPostEvents,2).*(1000./frameRateHz), (nanstd(highPostEvents,[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    title(['Post- Rew: ' num2str(chop(nanmean(cell2mat(HL_lickrate.low_postrew)),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_lickrate.high_postrew)),2)) ' Hz']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | lick bursts by rate: low (blue) & high (black)']);
    hold off;
    if seshType ~= 2
    savefig(fullfile(output_fn, '_cueAlignSpiking_byLickRate.fig'));
    saveas(tempFig, [output_fn  '_cueAlignSpiking_byLickRate.jpeg']);
    else
    savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byLickRate.fig']));
    saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byLickRate.jpeg']);
    end
    
if seshType ~= 2
    save(fullfile(output_fn, '_cueAlignLick.mat'), 'firstPostRewLickFrame', ...
        'lastPreRewLickFrame', 'tl_rew', 'firstPostRew_lickAlignEvents', 'lastPreRew_lickAlignEvents', ...
        'lickCueAlign', 'lickBurstStart', 'lickCounterVals', 'lickSearch_frames', 'lickDelay_frames',...
        'lickTC', 'postwin_frames', 'prewin_frames', 'frameRateHz', 'tt', 'ind_early_bst', 'ind_late_bst', ...
        'early_bst_time', 'late_bst_time', 'pct_precue_burst', 'postRew_lickAlignEvents', 'postLick_frames', ...
        'postRew_lickBurstStart','tl', 'ind_prerew_early_bst', 'ind_prerew_late_bst','postRew_lickAlign',...
        'preRew_lickBurstHz','postRew_lickBurstHz','ind_low_prerew','ind_high_prerew','ind_low_postrew',...
        'ind_high_postrew','ind_high_bst','ind_low_bst','HL_lickrate','lowAllEvents','highAllEvents',...
        'lowPreEvents','highPreEvents','lowPostEvents','highPostEvents','lastPreLickEvent',...
        'firstPostLickEvent','earlyBurstAlignEvents','lateBurstAlignEvents','earlyPostLickAlign','latePostLickAlign');
else
    save(fullfile(output_fn, [blockName '_cueAlignLick.mat']), 'firstPostRewLickFrame', ...
        'lastPreRewLickFrame', 'tl_rew', 'firstPostRew_lickAlignEvents', 'lastPreRew_lickAlignEvents', ...
        'lickCueAlign', 'lickBurstStart', 'lickCounterVals', 'lickSearch_frames', 'lickDelay_frames',...
        'lickTC', 'postwin_frames', 'prewin_frames', 'frameRateHz', 'tt', 'ind_early_bst', 'ind_late_bst', ...
        'early_bst_time', 'late_bst_time', 'pct_precue_burst', 'postRew_lickAlignEvents', 'postLick_frames', ...
        'postRew_lickBurstStart','tl', 'ind_prerew_early_bst', 'ind_prerew_late_bst','postRew_lickAlign',...
        'preRew_lickBurstHz','postRew_lickBurstHz','ind_low_prerew','ind_high_prerew','ind_low_postrew',...
        'ind_high_postrew','ind_high_bst','ind_low_bst','HL_lickrate','lowAllEvents','highAllEvents',...
        'lowPreEvents','highPreEvents','lowPostEvents','highPostEvents','lastPreLickEvent',...
        'firstPostLickEvent','earlyBurstAlignEvents','lateBurstAlignEvents','earlyPostLickAlign','latePostLickAlign');
end

    % piezo align figures
if seshType == 2
tempFig=figure;
subplot(2,1,1);
shadedErrorBar(tt,nanmean(groupAlign_piezo,2),nanstd(groupAlign_piezo,[],2)./sqrt(nTrials),'k');
ylabel('Piezo voltage');
xlabel('Time from cue (ms)');
if cue == 1
    xline(0,'g');
elseif cue == 3
    xline(0,'k');
elseif cue == 2
    xline(0,'r');
end
vline(767,'k');
subplot(2,1,2);
shadedErrorBar(tt,nanmean(lickCueAlign,2),nanstd(lickCueAlign,[],2)./sqrt(nTrials),'k');
ylabel('Lick rate');
xlabel('Time from cue (ms)');
if cue == 1
    xline(0,'g');
elseif cue == 3
    xline(0,'k');
elseif cue == 2
    xline(0,'r');
end
vline(767,'k');
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | lickAlign and piezoAlign | whole trial']);
if seshType ~= 2
savefig(fullfile(output_fn, '_avgTrialPiezo_abs.fig'));
saveas(tempFig, [output_fn  '_avgTrialPiezo_abs.jpeg']);
else
savefig(fullfile(output_fn, [blockName '_avgTrialPiezo_abs.fig']));
saveas(tempFig, [output_fn blockName '_avgTrialPiezo_abs.jpeg']);
end
    
elseif seshType ~= 2
    intB2Piezo = []; intRewPiezo = []; intB2Lick = []; intRewLick = [];
    for mouseIdx = 1:length(expt)
       intB2Piezo = [intB2Piezo groupAlign.piezo{1,mouseIdx}(:,logical(cell2mat(block2(1,mouseIdx))))];
       intRewPiezo = [intRewPiezo groupAlign.piezo{1,mouseIdx}(:,logical(cell2mat(rew(1,mouseIdx))))];
       intB2Lick = [intB2Lick groupAlign.licks{1,mouseIdx}(:,logical(cell2mat(block2(1,mouseIdx))))];
       intRewLick = [intRewLick groupAlign.licks{1,mouseIdx}(:,logical(cell2mat(rew(1,mouseIdx))))]; 
    end
    
tempFig=figure;
subplot(2,2,1);
shadedErrorBar(tt,nanmean(intB2Piezo,2),nanstd(intB2Piezo,[],2)./sqrt(nTrials),'k');
ylabel('Piezo voltage');
xlabel('Time from cue (ms)');
ylim([0 .1]);
vline(0,'r');
vline(767,'k');
subplot(2,2,3);
shadedErrorBar(tt,nanmean(intB2Lick,2),nanstd(intB2Lick,[],2)./sqrt(nTrials),'k');
ylabel('Lick rate');
xlabel('Time from cue (ms)');
ylim([0 .3])
vline(0,'r');
vline(767,'k');

subplot(2,2,2);
shadedErrorBar(tt,nanmean(intRewPiezo,2),nanstd(intRewPiezo,[],2)./sqrt(nTrials),'k');
ylabel('Piezo voltage');
xlabel('Time from cue (ms)');
ylim([0 .1]);
vline(0,'g');
vline(767,'k');
subplot(2,2,4);
shadedErrorBar(tt,nanmean(intRewLick,2),nanstd(intRewLick,[],2)./sqrt(nTrials),'k');
ylabel('Lick rate');
xlabel('Time from cue (ms)');
ylim([0 .3])
vline(0,'g');
vline(767,'k');
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | lickAlign and piezoAlign | whole trial']);
if seshType ~= 2
savefig(fullfile(output_fn, '_avgTrialPiezo_abs.fig'));
saveas(tempFig, [output_fn  '_avgTrialPiezo_abs.jpeg']);
else
savefig(fullfile(output_fn, [blockName '_avgTrialPiezo_abs.fig']));
saveas(tempFig, [output_fn blockName '_avgTrialPiezo_abs.jpeg']);
end
end

preRew_lickSearchRange = prewin_frames+lickDelay_frames:prewin_frames+lickDelay_frames+rewDelay_frames;
postRew_lickSearchRange = prewin_frames+rewDelay_frames+lickDelay_frames:prewin_frames+rewDelay_frames+rewDelay_frames;
preRew_piezoAmp = mean(groupAlign_piezo(preRew_lickSearchRange,:),1);% average across frames for each trial
postRew_piezoAmp = mean(groupAlign_piezo(postRew_lickSearchRange,:),1);

tempFig=figure;
subplot(1,2,1);
scatter(preRew_lickBurstHz, preRew_piezoAmp,'ok'); 
xlabel('Lick rate');
ylabel('Piezo voltage');
title('Pre reward');
xlim([0 10]);
ylim([-0.2 0.2]);
axis square;
subplot(1,2,2);
scatter(postRew_lickBurstHz, postRew_piezoAmp,'ok');
xlabel('Lick rate');
ylabel('Piezo voltage');
title('Post reward');
xlim([0 10]);
ylim([-0.2 0.2]);
axis square;
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Lick vs. Piezo | pre- and post-reward']);
if seshType ~= 2
savefig(fullfile(output_fn, '_LickvsPiezo_abs.fig'));
saveas(tempFig, [output_fn  '_LickvsPiezo_abs.jpeg']);
else
savefig(fullfile(output_fn, [blockName '_LickvsPiezo_abs.fig']));
saveas(tempFig, [output_fn blockName '_LickvsPiezo_abs.jpeg']);
end

% plot neural activity during trials of top and bottom 10/25% of movement amplitude
ntrials=0; low25PreRewEvent=[];    high25PreRewEvent=[];    low10PreRewEvent=[];    high10PreRewEvent=[];    low25PostRewEvent=[];    high25PostRewEvent=[];    low10PostRewEvent=[];    high10PostRewEvent=[];
for mouseIdx = 1:exptCount
    nTrial = (1:animalTrials(1,mouseIdx))+ntrials;  
    tempPiezoAmp = preRew_piezoAmp(1,nTrial);
    
[sortPiezoAmp sortPiezoAmp_ind] = sort(tempPiezoAmp,'ascend');
nnan = sum(isnan(tempPiezoAmp));
ind_low25piezo_prerew{1,mouseIdx} = sortPiezoAmp_ind(1:floor(animalTrials(1,mouseIdx)/4));
ind_high25piezo_prerew{1,mouseIdx} = sortPiezoAmp_ind(animalTrials(1,mouseIdx)-floor(animalTrials(1,mouseIdx)/4)+1:end-nnan);
HL_piezo.low25_prerew{1,mouseIdx} = nanmean(tempPiezoAmp(:,ind_low25piezo_prerew{1,mouseIdx}),2);
HL_piezo.high25_prerew{1,mouseIdx} = nanmean(tempPiezoAmp(:,ind_high25piezo_prerew{1,mouseIdx}),2);

ind_low10piezo_prerew{1,mouseIdx} = sortPiezoAmp_ind(1:floor(animalTrials(1,mouseIdx)/10));
ind_high10piezo_prerew{1,mouseIdx} = sortPiezoAmp_ind(animalTrials(1,mouseIdx)-floor(animalTrials(1,mouseIdx)/10)+1:end-nnan);
HL_piezo.low10_prerew{1,mouseIdx} = nanmean(tempPiezoAmp(:,ind_low10piezo_prerew{1,mouseIdx}),2);
HL_piezo.high10_prerew{1,mouseIdx} = nanmean(tempPiezoAmp(:,ind_high10piezo_prerew{1,mouseIdx}),2);
    
    tempPiezoAmp = postRew_piezoAmp(1,nTrial);

[sortPiezoAmp sortPiezoAmp_ind] = sort(tempPiezoAmp,'ascend');
nnan = sum(isnan(tempPiezoAmp));
ind_low25piezo_postrew{1,mouseIdx} = sortPiezoAmp_ind(1:floor(animalTrials(1,mouseIdx)/4));
ind_high25piezo_postrew{1,mouseIdx} = sortPiezoAmp_ind(animalTrials(1,mouseIdx)-floor(animalTrials(1,mouseIdx)/4)+1:end-nnan);
HL_piezo.low25_postrew{1,mouseIdx} = nanmean(tempPiezoAmp(:,ind_low25piezo_postrew{1,mouseIdx}),2);
HL_piezo.high25_postrew{1,mouseIdx} = nanmean(tempPiezoAmp(:,ind_high25piezo_postrew{1,mouseIdx}),2);

ind_low10piezo_postrew{1,mouseIdx} = sortPiezoAmp_ind(1:floor(animalTrials(1,mouseIdx)/10));
ind_high10piezo_postrew{1,mouseIdx} = sortPiezoAmp_ind(animalTrials(1,mouseIdx) -floor(animalTrials(1,mouseIdx)/10)+1:end-nnan);
HL_piezo.low10_postrew{1,mouseIdx} = nanmean(tempPiezoAmp(:,ind_low10piezo_postrew{1,mouseIdx}),2);
HL_piezo.high10_postrew{1,mouseIdx} = nanmean(tempPiezoAmp(:,ind_high10piezo_postrew{1,mouseIdx}),2);

    ntrials = max(nTrial);
    
    low25PreRewEvent = [low25PreRewEvent nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_low25piezo_prerew{1,mouseIdx}),3)];
    high25PreRewEvent = [high25PreRewEvent nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_high25piezo_prerew{1,mouseIdx}),3)];
    low10PreRewEvent = [low10PreRewEvent nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_low10piezo_prerew{1,mouseIdx}),3)];
    high10PreRewEvent = [high10PreRewEvent nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_high10piezo_prerew{1,mouseIdx}),3)];
    low25PostRewEvent = [low25PostRewEvent nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_low25piezo_postrew{1,mouseIdx}),3)];
    high25PostRewEvent = [high25PostRewEvent nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_high25piezo_postrew{1,mouseIdx}),3)];
    low10PostRewEvent = [low10PostRewEvent nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_low10piezo_postrew{1,mouseIdx}),3)];
    high10PostRewEvent = [high10PostRewEvent nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_high10piezo_postrew{1,mouseIdx}),3)];
end

nIC = size(groupAlign_events,2);
tempFig=figure;
subplot(2,1,1);
shadedErrorBar(tt,nanmean(low25PreRewEvent,2).*(1000./frameRateHz), (nanstd(low25PreRewEvent,[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on;
shadedErrorBar(tt,nanmean(high25PreRewEvent,2).*(1000./frameRateHz), (nanstd(high25PreRewEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
if cue == 1
    vline(0,'g');
elseif cue == 3
    vline(0,'k');
elseif cue == 2
    vline(0,'r');
end
vline(767,'k');
ylim([0 inf]);
title(['Pre-reward: ' num2str(chop(nanmean(cell2mat(HL_piezo.low25_prerew)),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.high25_prerew)),2)) ' V']);
subplot(2,1,2);
shadedErrorBar(tt,nanmean(low25PostRewEvent,2).*(1000./frameRateHz), (nanstd(low25PostRewEvent,[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on;
shadedErrorBar(tt,nanmean(high25PostRewEvent,2).*(1000./frameRateHz), (nanstd(high25PostRewEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
xlabel('Time from cue (ms)');
ylabel('Spike rate (Hz)');
if cue == 1
    vline(0,'g');
elseif cue == 3
    vline(0,'k');
elseif cue == 2
    vline(0,'r');
end
vline(767,'k');
ylim([0 inf]);
title(['Post-reward: ' num2str(chop(nanmean(cell2mat(HL_piezo.low25_postrew)),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.high25_postrew)),2)) ' V']);
hold on;
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | cue aligned motion: low25(blu) high25(blk)']);
hold off;
if seshType ~= 2
savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoAmp25_abs.fig'))
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoAmp25_abs.jpeg']);
else
savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoAmp25_abs.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoAmp25_abs.jpeg']);
end

tempFig=figure;
subplot(2,1,1);
shadedErrorBar(tt,nanmean(low10PreRewEvent,2).*(1000./frameRateHz), (nanstd(low10PreRewEvent,[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on
shadedErrorBar(tt,nanmean(high10PreRewEvent,2).*(1000./frameRateHz), (nanstd(high10PreRewEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
xlabel('Time from cue (ms)');
ylabel('Spike rate (Hz)');
if cue == 1
    vline(0,'g');
elseif cue == 3
    vline(0,'k');
elseif cue == 2
    vline(0,'r');
end
vline(767,'k');
ylim([0 inf]);
title(['Pre-reward: ' num2str(chop(nanmean(cell2mat(HL_piezo.low10_prerew)),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.high10_prerew)),2)) ' V']);
subplot(2,1,2);
shadedErrorBar(tt,nanmean(low10PostRewEvent,2).*(1000./frameRateHz), (nanstd(low10PostRewEvent,[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on;
shadedErrorBar(tt,nanmean(high10PostRewEvent,2).*(1000./frameRateHz), (nanstd(high10PostRewEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
if cue == 1
    vline(0,'g');
elseif cue == 3
    vline(0,'k');
elseif cue == 2
    vline(0,'r');
end
vline(767,'k');
ylim([0 inf]);
title(['Post-reward: ' num2str(chop(nanmean(cell2mat(HL_piezo.low10_postrew)),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.high10_postrew)),2)) ' V']);
hold off;
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | piezo bursts by rate: low10(blu) & high10(blk)']);
if seshType ~= 2
savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoAmp10_abs.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoAmp10_abs.jpeg']);
else
savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoAmp10_abs.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoAmp10_abs.jpeg']);
end
%%
% sort out trials by looking at time of big movements: 25% and 10%
%nTrial=0;
for mouseIdx = 1:exptCount    
    ntrials = (animalTrials(1,mouseIdx));%+nTrial;
    
piezo_base_std{1,mouseIdx} = nanstd(reshape(groupAlign.piezo{1,mouseIdx}(1:prewin_frames,:),[prewin_frames.*ntrials 1]),[],1);
piezo_base_avg{1,mouseIdx} = nanmean(reshape(groupAlign.piezo{1,mouseIdx}(1:prewin_frames,:),[prewin_frames.*ntrials 1]),1);
groupAlign_piezo_thresh1{1,mouseIdx} = zeros(size(groupAlign.piezo{1,mouseIdx}));%>1std <2std
groupAlign_piezo_thresh1ALL{1,mouseIdx} = zeros(size(groupAlign.piezo{1,mouseIdx}));%>1std
groupAlign_piezo_thresh2ALL{1,mouseIdx} = zeros(size(groupAlign.piezo{1,mouseIdx})); %>2std
groupAlign_piezo_thresh2{1,mouseIdx} = zeros(size(groupAlign.piezo{1,mouseIdx}));%>2std <3std
groupAlign_piezo_thresh3{1,mouseIdx} = zeros(size(groupAlign.piezo{1,mouseIdx}));%>3std
piezoSearchRange = postwin_frames+lickDelay_frames:size(lickCueAlign,1)-lickSearch_frames;
piezoStart{1,mouseIdx} = nan(1,length(ntrials));

    %nTrial = max(ntrials);
end
% for each trial, frames that has a voltage>n*std = 1, others = 0
nTrial=0;
for mouseIdx = 1:exptCount
    ntrials = animalTrials(1,mouseIdx);%+nTrial;
for itrial = 1+nTrial:ntrials+nTrial
    tCo = 1:ntrials; tNo = tCo(1,itrial-nTrial); 
    groupAlign_piezo_thresh1{1,mouseIdx}(find(groupAlign.piezo{1,mouseIdx}(:,tNo)>=(piezo_base_avg{1,mouseIdx} + piezo_base_std{1,mouseIdx}.*1) & groupAlign.piezo{1,mouseIdx}(:,tNo)<(piezo_base_avg{1,mouseIdx} + piezo_base_std{1,mouseIdx}.*2)),tNo) = 1;
    groupAlign_piezo_thresh2{1,mouseIdx}(find(groupAlign.piezo{1,mouseIdx}(:,tNo)>=(piezo_base_avg{1,mouseIdx} + piezo_base_std{1,mouseIdx}.*2) & groupAlign.piezo{1,mouseIdx}(:,tNo)<(piezo_base_avg{1,mouseIdx} + piezo_base_std{1,mouseIdx}.*3)),tNo) = 1;
    groupAlign_piezo_thresh3{1,mouseIdx}(find(groupAlign.piezo{1,mouseIdx}(:,tNo)>=(piezo_base_avg{1,mouseIdx} + piezo_base_std{1,mouseIdx}.*3)),tNo) = 1;
    groupAlign_piezo_thresh1ALL{1,mouseIdx}(find(groupAlign.piezo{1,mouseIdx}(:,tNo)>=(piezo_base_avg{1,mouseIdx} + piezo_base_std{1,mouseIdx}.*1)),tNo) = 1;
    groupAlign_piezo_thresh2ALL{1,mouseIdx}(find(groupAlign.piezo{1,mouseIdx}(:,tNo)>=(piezo_base_avg{1,mouseIdx} + piezo_base_std{1,mouseIdx}.*2)),tNo) = 1;
    if  find(groupAlign_piezo_thresh3{1,mouseIdx}(piezoSearchRange,tNo),1,'first')
        piezoStart{1,mouseIdx}(:,tNo) = find(groupAlign_piezo_thresh3{1,mouseIdx}(piezoSearchRange,tNo),1,'first');% the start of very big movements in each trial (3std>baseline)
    else
        piezoStart{1,mouseIdx}(:,tNo) = 0;
    end
end
    nTrial=sum(animalTrials(1,1:mouseIdx));
end

                    %% CS+ and CS- piezo align with piezo trace!! %%
nTrial=0; ntrialsB2=[]; ntrialsRew=[]; earlyPiezoEventB2=[];    latePiezoEventB2=[];    early25PiezoEventB2=[];    late25PiezoEventB2=[];    early10PiezoEventB2=[];    late10PiezoEventB2=[];earlyPiezoEventRew=[];    latePiezoEventRew=[];    early25PiezoEventRew=[];    late25PiezoEventRew=[];    early10PiezoEventRew=[];    late10PiezoEventRew=[];
t3earlyPiezoB2=[]; t3earlyPiezoRew=[]; t3latePiezoB2=[]; t3latePiezoRew=[]; t3early25PiezoB2=[]; t3early25PiezoRew=[]; t3late25PiezoB2=[]; t3late25PiezoRew=[]; t3early10PiezoB2=[]; t3early10PiezoRew=[]; t3late10PiezoB2=[]; t3late10PiezoRew=[];
t2earlyPiezoB2=[]; t2earlyPiezoRew=[]; t2latePiezoB2=[]; t2latePiezoRew=[]; t2early25PiezoB2=[]; t2early25PiezoRew=[]; t2late25PiezoB2=[]; t2late25PiezoRew=[]; t2early10PiezoB2=[]; t2early10PiezoRew=[]; t2late10PiezoB2=[]; t2late10PiezoRew=[];
t1earlyPiezoB2=[]; t1earlyPiezoRew=[]; t1latePiezoB2=[]; t1latePiezoRew=[]; t1early25PiezoB2=[]; t1early25PiezoRew=[]; t1late25PiezoB2=[]; t1late25PiezoRew=[]; t1early10PiezoB2=[]; t1early10PiezoRew=[]; t1late10PiezoB2=[]; t1late10PiezoRew=[];
% % % % % POTENTIALLY 
t1earlyPiezoV2B2 = [];t1earlyPiezoV2Rew = [];t1latePiezoV2B2 = [];t1latePiezoV2Rew = [];t1early25PiezoV2B2 = [];t1early25PiezoV2Rew = [];t1late25PiezoV2B2 = [];t1late25PiezoV2Rew = [];t1early10PiezoV2B2 = [];t1early10PiezoV2Rew = [];t1late10PiezoV2B2 = [];t1late10PiezoV2Rew = [];
% % % % % TEMPORARY
for mouseIdx = 1:exptCount
    ntrials = (1:animalTrials(1,mouseIdx))+nTrial;
    for m=1:length(ntrials)
        if block2{1,mouseIdx}(1,m) == 1
    ntrialsB2 = [ntrialsB2 m];%+nTrial];
        elseif rew{1,mouseIdx}(1,m) == 1
    ntrialsRew = [ntrialsRew m];%+nTrial];
        end
    end
    
    tempPiezoB2{1,mouseIdx} = piezoStart{1,mouseIdx}(1,ntrialsB2);
    tempPiezoRew{1,mouseIdx} = piezoStart{1,mouseIdx}(1,ntrialsRew);
[sortPiezoStartB2 sortPiezoStart_indB2] = sort(tempPiezoB2{1,mouseIdx},'ascend');
nonzero = find(sortPiezoStartB2); sortPiezoStart_indB2=sortPiezoStart_indB2(1,nonzero);
[sortPiezoStartRew sortPiezoStart_indRew] = sort(tempPiezoRew{1,mouseIdx},'ascend');
nonzero = find(sortPiezoStartRew); sortPiezoStart_indRew=sortPiezoStart_indRew(1,nonzero);
nnanB2 = sum(isnan(tempPiezoB2{1,mouseIdx}));
nnanRew = sum(isnan(tempPiezoRew{1,mouseIdx}));
ind_earlypiezo_B2{1,mouseIdx} = sortPiezoStart_indB2(1:floor((length(sortPiezoStart_indB2)-nnanB2)/4));% top 25% of trials with early big movement
ind_earlypiezo_rew{1,mouseIdx} = sortPiezoStart_indRew(1:floor((length(sortPiezoStart_indRew)-nnanRew)/4));% top 25% of trials with early big movement
ind_latepiezo_B2{1,mouseIdx} = sortPiezoStart_indB2((length(sortPiezoStart_indB2)-nnanB2)-floor((length(sortPiezoStart_indB2)-nnanB2)/4)+1:end-nnanB2);
ind_latepiezo_rew{1,mouseIdx} = sortPiezoStart_indRew((length(sortPiezoStart_indRew)-nnanRew)-floor((length(sortPiezoStart_indRew)-nnanRew)/4)+1:end-nnanRew);
ind_allearlypiezo25_B2{1,mouseIdx} = sortPiezoStart_indB2(1:floor(length(sortPiezoStart_indB2)/4));
ind_allearlypiezo25_rew{1,mouseIdx} = sortPiezoStart_indRew(1:floor(length(sortPiezoStart_indRew)/4));
ind_alllatepiezo25_B2{1,mouseIdx} = sortPiezoStart_indB2(length(sortPiezoStart_indB2)-floor(length(sortPiezoStart_indB2)/4)+1:end);
ind_alllatepiezo25_rew{1,mouseIdx} = sortPiezoStart_indRew(length(sortPiezoStart_indRew)-floor(length(sortPiezoStart_indRew)/4)+1:end);
ind_allearlypiezo10_B2{1,mouseIdx} = sortPiezoStart_indB2(1:floor(length(sortPiezoStart_indB2)/10)); % top 10% of all trials, regardless of nans. does it make sense? 
ind_allearlypiezo10_rew{1,mouseIdx} = sortPiezoStart_indRew(1:floor(length(sortPiezoStart_indRew)/10)); % top 10% of all trials, regardless of nans. does it make sense? 
ind_alllatepiezo10_B2{1,mouseIdx} = sortPiezoStart_indB2(length(sortPiezoStart_indB2)-floor(length(sortPiezoStart_indB2)/10)+1:end);
ind_alllatepiezo10_rew{1,mouseIdx} = sortPiezoStart_indRew(length(sortPiezoStart_indRew)-floor(length(sortPiezoStart_indRew)/10)+1:end);
HL_piezo.early_B2{1,mouseIdx} = (tempPiezoB2{1,mouseIdx}(:,ind_earlypiezo_B2{1,mouseIdx}));
HL_piezo.early_rew{1,mouseIdx} = (tempPiezoRew{1,mouseIdx}(:,ind_earlypiezo_rew{1,mouseIdx}));
HL_piezo.late_B2{1,mouseIdx} = (tempPiezoB2{1,mouseIdx}(:,ind_latepiezo_B2{1,mouseIdx}));
HL_piezo.late_rew{1,mouseIdx} = (tempPiezoRew{1,mouseIdx}(:,ind_latepiezo_rew{1,mouseIdx}));
    nTrial=max(ntrials);
  
    earlyPiezoEventB2 = [earlyPiezoEventB2 nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_earlypiezo_B2{1,mouseIdx}),3)];
    earlyPiezoEventRew = [earlyPiezoEventRew nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_earlypiezo_rew{1,mouseIdx}),3)];
    latePiezoEventB2 = [latePiezoEventB2 nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_latepiezo_B2{1,mouseIdx}),3)];
    latePiezoEventRew = [latePiezoEventRew nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_latepiezo_rew{1,mouseIdx}),3)];
    early25PiezoEventB2 = [early25PiezoEventB2 nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_allearlypiezo25_B2{1,mouseIdx}),3)];
    early25PiezoEventRew = [early25PiezoEventRew nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_allearlypiezo25_rew{1,mouseIdx}),3)];
    late25PiezoEventB2 = [late25PiezoEventB2 nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_alllatepiezo25_B2{1,mouseIdx}),3)];
    late25PiezoEventRew = [late25PiezoEventRew nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_alllatepiezo25_rew{1,mouseIdx}),3)];
    early10PiezoEventB2 = [early10PiezoEventB2 nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_allearlypiezo10_B2{1,mouseIdx}),3)];
    early10PiezoEventRew = [early10PiezoEventRew nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_allearlypiezo10_rew{1,mouseIdx}),3)];
    late10PiezoEventB2 = [late10PiezoEventB2 nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_alllatepiezo10_B2{1,mouseIdx}),3)];
    late10PiezoEventRew = [late10PiezoEventRew nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_alllatepiezo10_rew{1,mouseIdx}),3)];
    % % % % % % % % % % % % % % % % % % TESTING SOMETHING BELOW
    t1earlyPiezoV2B2 = [t1earlyPiezoV2B2 groupAlign.piezo{1,mouseIdx}(:,ind_earlypiezo_B2{1,mouseIdx})];
    t1earlyPiezoV2Rew = [t1earlyPiezoV2Rew groupAlign.piezo{1,mouseIdx}(:,ind_earlypiezo_rew{1,mouseIdx})];
    t1latePiezoV2B2 = [t1latePiezoV2B2 groupAlign.piezo{1,mouseIdx}(:,ind_latepiezo_B2{1,mouseIdx})];
    t1latePiezoV2Rew = [t1latePiezoV2Rew groupAlign.piezo{1,mouseIdx}(:,ind_latepiezo_rew{1,mouseIdx})];
    t1early25PiezoV2B2 = [t1early25PiezoV2B2 groupAlign.piezo{1,mouseIdx}(:,ind_allearlypiezo25_B2{1,mouseIdx})];
    t1early25PiezoV2Rew = [t1early25PiezoV2Rew groupAlign.piezo{1,mouseIdx}(:,ind_allearlypiezo25_rew{1,mouseIdx})];
    t1late25PiezoV2B2 = [t1late25PiezoV2B2 groupAlign.piezo{1,mouseIdx}(:,ind_alllatepiezo25_B2{1,mouseIdx})];
    t1late25PiezoV2Rew = [t1late25PiezoV2Rew groupAlign.piezo{1,mouseIdx}(:,ind_alllatepiezo25_rew{1,mouseIdx})];
    t1early10PiezoV2B2 = [t1early10PiezoV2B2 groupAlign.piezo{1,mouseIdx}(:,ind_allearlypiezo10_B2{1,mouseIdx})];
    t1early10PiezoV2Rew = [t1early10PiezoV2Rew groupAlign.piezo{1,mouseIdx}(:,ind_allearlypiezo10_rew{1,mouseIdx})];
    t1late10PiezoV2B2 = [t1late10PiezoV2B2 groupAlign.piezo{1,mouseIdx}(:,ind_alllatepiezo10_B2{1,mouseIdx})];
    t1late10PiezoV2Rew = [t1late10PiezoV2Rew groupAlign.piezo{1,mouseIdx}(:,ind_alllatepiezo10_rew{1,mouseIdx})];
    % % % % % % % % % % % % % % % % % % TESTING SOMETHING ABOVE
    t1earlyPiezoB2 = [t1earlyPiezoB2 groupAlign_piezo_thresh1{1,mouseIdx}(:,ind_earlypiezo_B2{1,mouseIdx})];
    t1earlyPiezoRew = [t1earlyPiezoRew groupAlign_piezo_thresh1{1,mouseIdx}(:,ind_earlypiezo_rew{1,mouseIdx})];
    t1latePiezoB2 = [t1latePiezoB2 groupAlign_piezo_thresh1{1,mouseIdx}(:,ind_latepiezo_B2{1,mouseIdx})];
    t1latePiezoRew = [t1latePiezoRew groupAlign_piezo_thresh1{1,mouseIdx}(:,ind_latepiezo_rew{1,mouseIdx})];
    t1early25PiezoB2 = [t1early25PiezoB2 groupAlign_piezo_thresh1{1,mouseIdx}(:,ind_allearlypiezo25_B2{1,mouseIdx})];
    t1early25PiezoRew = [t1early25PiezoRew groupAlign_piezo_thresh1{1,mouseIdx}(:,ind_allearlypiezo25_rew{1,mouseIdx})];
    t1late25PiezoB2 = [t1late25PiezoB2 groupAlign_piezo_thresh1{1,mouseIdx}(:,ind_alllatepiezo25_B2{1,mouseIdx})];
    t1late25PiezoRew = [t1late25PiezoRew groupAlign_piezo_thresh1{1,mouseIdx}(:,ind_alllatepiezo25_rew{1,mouseIdx})];
    t1early10PiezoB2 = [t1early10PiezoB2 groupAlign_piezo_thresh1{1,mouseIdx}(:,ind_allearlypiezo10_B2{1,mouseIdx})];
    t1early10PiezoRew = [t1early10PiezoRew groupAlign_piezo_thresh1{1,mouseIdx}(:,ind_allearlypiezo10_rew{1,mouseIdx})];
    t1late10PiezoB2 = [t1late10PiezoB2 groupAlign_piezo_thresh1{1,mouseIdx}(:,ind_alllatepiezo10_B2{1,mouseIdx})];
    t1late10PiezoRew = [t1late10PiezoRew groupAlign_piezo_thresh1{1,mouseIdx}(:,ind_alllatepiezo10_rew{1,mouseIdx})];

    t2earlyPiezoB2 = [t2earlyPiezoB2 groupAlign_piezo_thresh2{1,mouseIdx}(:,ind_earlypiezo_B2{1,mouseIdx})];
    t2earlyPiezoRew = [t2earlyPiezoRew groupAlign_piezo_thresh2{1,mouseIdx}(:,ind_earlypiezo_rew{1,mouseIdx})];
    t2latePiezoB2 = [t2latePiezoB2 groupAlign_piezo_thresh2{1,mouseIdx}(:,ind_latepiezo_B2{1,mouseIdx})];
    t2latePiezoRew = [t2latePiezoRew groupAlign_piezo_thresh2{1,mouseIdx}(:,ind_latepiezo_rew{1,mouseIdx})];
    t2early25PiezoB2 = [t2early25PiezoB2 groupAlign_piezo_thresh2{1,mouseIdx}(:,ind_allearlypiezo25_B2{1,mouseIdx})];
    t2early25PiezoRew = [t2early25PiezoRew groupAlign_piezo_thresh2{1,mouseIdx}(:,ind_allearlypiezo25_rew{1,mouseIdx})];
    t2late25PiezoB2 = [t2late25PiezoB2 groupAlign_piezo_thresh2{1,mouseIdx}(:,ind_alllatepiezo25_B2{1,mouseIdx})];
    t2late25PiezoRew = [t2late25PiezoRew groupAlign_piezo_thresh2{1,mouseIdx}(:,ind_alllatepiezo25_rew{1,mouseIdx})];
    t2early10PiezoB2 = [t2early10PiezoB2 groupAlign_piezo_thresh2{1,mouseIdx}(:,ind_allearlypiezo10_B2{1,mouseIdx})];
    t2early10PiezoRew = [t2early10PiezoRew groupAlign_piezo_thresh2{1,mouseIdx}(:,ind_allearlypiezo10_rew{1,mouseIdx})];
    t2late10PiezoB2 = [t2late10PiezoB2 groupAlign_piezo_thresh2{1,mouseIdx}(:,ind_alllatepiezo10_B2{1,mouseIdx})];
    t2late10PiezoRew = [t2late10PiezoRew groupAlign_piezo_thresh2{1,mouseIdx}(:,ind_alllatepiezo10_rew{1,mouseIdx})];
    
    t3earlyPiezoB2 = [t3earlyPiezoB2 groupAlign_piezo_thresh3{1,mouseIdx}(:,ind_earlypiezo_B2{1,mouseIdx})];
    t3earlyPiezoRew = [t3earlyPiezoRew groupAlign_piezo_thresh3{1,mouseIdx}(:,ind_earlypiezo_rew{1,mouseIdx})];
    t3latePiezoB2 = [t3latePiezoB2 groupAlign_piezo_thresh3{1,mouseIdx}(:,ind_latepiezo_B2{1,mouseIdx})];
    t3latePiezoRew = [t3latePiezoRew groupAlign_piezo_thresh3{1,mouseIdx}(:,ind_latepiezo_rew{1,mouseIdx})];
    t3early25PiezoB2 = [t3early25PiezoB2 groupAlign_piezo_thresh3{1,mouseIdx}(:,ind_allearlypiezo25_B2{1,mouseIdx})];
    t3early25PiezoRew = [t3early25PiezoRew groupAlign_piezo_thresh3{1,mouseIdx}(:,ind_allearlypiezo25_rew{1,mouseIdx})];
    t3late25PiezoB2 = [t3late25PiezoB2 groupAlign_piezo_thresh3{1,mouseIdx}(:,ind_alllatepiezo25_B2{1,mouseIdx})];
    t3late25PiezoRew = [t3late25PiezoRew groupAlign_piezo_thresh3{1,mouseIdx}(:,ind_alllatepiezo25_rew{1,mouseIdx})];
    t3early10PiezoB2 = [t3early10PiezoB2 groupAlign_piezo_thresh3{1,mouseIdx}(:,ind_allearlypiezo10_B2{1,mouseIdx})];
    t3early10PiezoRew = [t3early10PiezoRew groupAlign_piezo_thresh3{1,mouseIdx}(:,ind_allearlypiezo10_rew{1,mouseIdx})];
    t3late10PiezoB2 = [t3late10PiezoB2 groupAlign_piezo_thresh3{1,mouseIdx}(:,ind_alllatepiezo10_B2{1,mouseIdx})];
    t3late10PiezoRew = [t3late10PiezoRew groupAlign_piezo_thresh3{1,mouseIdx}(:,ind_alllatepiezo10_rew{1,mouseIdx})];
    
ntrialsB2 = [];
ntrialsRew = [];
end

%% MOUSE EVENT and MOUSE PIEZO (MOST MOVEMENT, MED MOVEMENT, LEAST MOVEMENT)
%{

for i = 1:exptCount
for itrial = 1:animalTrials(1,i)
    for f = 3:size(lickCueAlign,1)-2        
        piezoWindow{1,i}(f,itrial) = nanmean(groupAlign.piezo{1,i}(f-2:f+2,itrial));       
    end
end 
end
for i = 1:length(piezoWindow)
[tempY(:,i),tempX(:,i)]=sort(nanmean(piezoWindow{1,i}(5:148,:),2),'descend');
yH(1,i) = tempY(10,i); xH(1,i) = tempX(10,i);
end
tempFig=figure; 
randcolor = {[0.2 0.8 0.8],[0.2 0 0],[1 0.5 0],[0 0.5 1],[0 0.6 0.3],[1 0.2 0.2],[0 0.8 1],[0.5 0.8 0]};
for i = 1:exptCount
    plot(tt(1,3:(size(lickCueAlign,1)-4)),nanmean(piezoWindow{1,i}(5:148,:),2)','color',randcolor{1,i});
    text(3500,yH(1,i),cell2mat(expt{1,i}.mouse),'color',randcolor{1,i});%
    hold on;
    %text((tt(1,sortPiezoRaw_ind{t}(1,1))), sortPiezoRaw{t}(1,1), num2str(t),'color','k');
end

[maxMotionMouse,maxMotionMouse_idx]=sort(yH,'descend');
highMoMo = [expt{1,maxMotionMouse_idx(1,1)}.mouse,expt{1,maxMotionMouse_idx(1,2)}.mouse];
medMoMo = [expt{1,maxMotionMouse_idx(1,4)}.mouse,expt{1,maxMotionMouse_idx(1,5)}.mouse];
lowMoMo = [expt{1,maxMotionMouse_idx(1,7)}.mouse,expt{1,maxMotionMouse_idx(1,8)}.mouse];

plot(tt,nanmean(nanmean(groupAlign.events{1,maxMotionMouse_idx(1,1:2)},3),2).*(1000./frameRateHz),(nanstd(nanmean(groupAlign.events{1,maxMotionMouse_idx(1,1:2)},3)).*(1000./frameRateHz))./sqrt(nIC));
plot(tt,nanmean(nanmean(groupAlign.events{1,maxMotionMouse_idx(1,4:5)},3),2).*(1000./frameRateHz),(nanstd(nanmean(groupAlign.events{1,maxMotionMouse_idx(1,4:5)},3)).*(1000./frameRateHz))./sqrt(nIC));
plot(tt,nanmean(nanmean(groupAlign.events{1,maxMotionMouse_idx(1,7:8)},3),2).*(1000./frameRateHz),(nanstd(nanmean(groupAlign.events{1,maxMotionMouse_idx(1,7:8)},3)).*(1000./frameRateHz))./sqrt(nIC));

randcolor = {[1 0.5 0],[1 0.2 0.2],[0 0.8 1],[0 0.5 1],[0 0.6 0.3],[0.5 0.8 0]};
tempFig = figure;
subplot(2,1,1);
plot(tt,nanmean(nanmean(groupAlign.events{1,maxMotionMouse_idx(1,1)},3),2).*(1000./frameRateHz),'color',randcolor{1,1}); %,(nanstd(nanmean(groupAlign.events{1,maxMotionMouse_idx(1,1)},3),[],2).*(1000./frameRateHz))./sqrt(nIC)
hold on;
plot(tt,nanmean(nanmean(groupAlign.events{1,maxMotionMouse_idx(1,2)},3),2).*(1000./frameRateHz),'color',randcolor{1,2});
plot(tt,nanmean(nanmean(groupAlign.events{1,maxMotionMouse_idx(1,4)},3),2).*(1000./frameRateHz),'color',randcolor{1,3});
plot(tt,nanmean(nanmean(groupAlign.events{1,maxMotionMouse_idx(1,5)},3),2).*(1000./frameRateHz),'color',randcolor{1,4});
plot(tt,nanmean(nanmean(groupAlign.events{1,maxMotionMouse_idx(1,7)},3),2).*(1000./frameRateHz),'color',randcolor{1,5});
plot(tt,nanmean(nanmean(groupAlign.events{1,maxMotionMouse_idx(1,8)},3),2).*(1000./frameRateHz),'color',randcolor{1,6});
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
vline(0,'k');
vline(767,'k');

subplot(2,1,2);
plot(tt(1,3:(size(lickCueAlign,1)-4)),nanmean(piezoWindow{1,maxMotionMouse_idx(1,1)}(5:148,:),2)','color',randcolor{1,1});
text(3500,yH(1,maxMotionMouse_idx(1,1)),cell2mat(expt{1,maxMotionMouse_idx(1,1)}.mouse),'color',randcolor{1,1});%
hold on;
plot(tt(1,3:(size(lickCueAlign,1)-4)),nanmean(piezoWindow{1,maxMotionMouse_idx(1,2)}(5:148,:),2)','color',randcolor{1,2});
text(3500,yH(1,maxMotionMouse_idx(1,2)),cell2mat(expt{1,maxMotionMouse_idx(1,2)}.mouse),'color',randcolor{1,2});%
plot(tt(1,3:(size(lickCueAlign,1)-4)),nanmean(piezoWindow{1,maxMotionMouse_idx(1,4)}(5:148,:),2)','color',randcolor{1,3});
text(3500,yH(1,maxMotionMouse_idx(1,4)),cell2mat(expt{1,maxMotionMouse_idx(1,4)}.mouse),'color',randcolor{1,3});%
plot(tt(1,3:(size(lickCueAlign,1)-4)),nanmean(piezoWindow{1,maxMotionMouse_idx(1,5)}(5:148,:),2)','color',randcolor{1,4});
text(3500,yH(1,maxMotionMouse_idx(1,5)),cell2mat(expt{1,maxMotionMouse_idx(1,5)}.mouse),'color',randcolor{1,4});%
plot(tt(1,3:(size(lickCueAlign,1)-4)),nanmean(piezoWindow{1,maxMotionMouse_idx(1,7)}(5:148,:),2)','color',randcolor{1,5});
text(3500,yH(1,maxMotionMouse_idx(1,7)),cell2mat(expt{1,maxMotionMouse_idx(1,7)}.mouse),'color',randcolor{1,5});%
plot(tt(1,3:(size(lickCueAlign,1)-4)),nanmean(piezoWindow{1,maxMotionMouse_idx(1,8)}(5:148,:),2)','color',randcolor{1,6});
text(3500,yH(1,maxMotionMouse_idx(1,8)),cell2mat(expt{1,maxMotionMouse_idx(1,8)}.mouse),'color',randcolor{1,6});%
xlabel('Time from cue');
ylabel('Piezo voltage');
vline(0,'k');
vline(767,'k');
%}
%% MOST VS. LEAST MOVEMENT - TOP AND BOTTOM 10% OF TRIALS FOR EACH MOUSE
mostMovementAlignedEvents=[]; leastMovementAlignedEvents=[]; mostTrialMovement=[]; leastTrialMovement=[];
for mouseIdx = 1:exptCount
[TrialMovement, TrialMovement_ind] = sort(nanmean(groupAlign.piezo{1,mouseIdx},1));
mostTrialMovement_ind{1,mouseIdx} = TrialMovement_ind(:,end-(floor(size(groupAlign.piezo{1,mouseIdx},2)./10)):end);
leastTrialMovement_ind{1,mouseIdx} = TrialMovement_ind(:,1:floor(size(groupAlign.piezo{1,mouseIdx},2)./10));

mostMovementAlignedEvents = [mostMovementAlignedEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,mostTrialMovement_ind{1,mouseIdx}),3)];
leastMovementAlignedEvents = [leastMovementAlignedEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,leastTrialMovement_ind{1,mouseIdx}),3)];
mostTrialMovement = [mostTrialMovement groupAlign.piezo{1,mouseIdx}(:,mostTrialMovement_ind{1,mouseIdx})];
leastTrialMovement = [leastTrialMovement groupAlign.piezo{1,mouseIdx}(:,leastTrialMovement_ind{1,mouseIdx})];
end

tempFig = figure;
subplot(2,1,1);
plot(tt,nanmean(mostMovementAlignedEvents,2).*(1000./frameRateHz),'color',[1 0.5 0]);
hold on;
plot(tt,nanmean(leastMovementAlignedEvents,2).*(1000./frameRateHz),'color',[0 0.5 1]);
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
vline(0,'k');
vline(767,'k');
legend('top 10%','bottom 10%','Location','NorthEast');
title('Neural data aligned to trials with top 10% (org) and bottom 10% (blu) of movement');
subplot(2,1,2);
plot(tt,nanmean(mostTrialMovement,2),'color',[1 0.5 0]);
hold on;
plot(tt,nanmean(leastTrialMovement,2),'color',[0 0.5 1]);
xlabel('Time from cue');
ylabel('Piezo voltage');
vline(0,'k');
vline(767,'k');
title('Piezo data for trials with top 10% (org) and bottom 10% (blu) of movement');
hold off;
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | neural aligned to trials with most and least movement']);
if seshType ~= 2
savefig(fullfile(output_fn, '_mostLeastMotionAlign.fig'));
saveas(tempFig, [output_fn  '_mostLeastMotionAlign.jpeg']);
else
savefig(fullfile(output_fn, [blockName '_mostLeastMotionAlign.fig']));
saveas(tempFig, [output_fn blockName '_mostLeastMotionAlign.jpeg']);
end



%% CS-
tempFig=figure;
subplot(2,1,1);
shadedErrorBar(tt,nanmean(earlyPiezoEventB2,2).*(1000./frameRateHz), (nanstd(earlyPiezoEventB2,[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on;
shadedErrorBar(tt,nanmean(latePiezoEventB2,2).*(1000./frameRateHz), (nanstd(latePiezoEventB2,[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
for mouseIdx=1:exptCount
scatter(tempPiezoB2{1,mouseIdx}(:,ind_earlypiezo_B2{1,mouseIdx})*(1000./frameRateHz),zeros(1,length(ind_earlypiezo_B2{1,mouseIdx}))+3,'ob');
scatter(tempPiezoB2{1,mouseIdx}(:,ind_latepiezo_B2{1,mouseIdx})*(1000./frameRateHz),zeros(1,length(ind_latepiezo_B2{1,mouseIdx}))+3,'ok');
end
vline(0,'r');
vline(767,'k');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 inf]);
title(['earliest 25% vs. latest 25% of motion of CS- trials aligned to cue: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_B2)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_B2)).*(1000./frameRateHz),2)) ' ms']);
subplot(2,1,2);
shadedErrorBar(tt,nanmean((t3earlyPiezoB2),2),nanstd((t3earlyPiezoB2),[],2)./sqrt(sum(cell2mat(block2))),'b');
hold on;
shadedErrorBar(tt,nanmean((t3latePiezoB2),2),nanstd((t3latePiezoB2),[],2)./sqrt(sum(cell2mat(block2))),'k');
vline(0,'r');
vline(767,'k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
title(['earliest 25% vs. latest 25% of motion of CS- trials: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_B2)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_B2)).*(1000./frameRateHz),2)) ' ms']);
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | motion bursts of CS- trials: early(blu) & late(blk)']);
if seshType ~= 2
savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatency_B2.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatency_B2.jpeg']);
else
savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatency_B2.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatency_B2.jpeg']);
end

% 1 SD over but under 2 SD of mean piezo voltage
tempFig=figure;
subplot(2,1,1);
shadedErrorBar(tt,nanmean(earlyPiezoEventB2,2).*(1000./frameRateHz), (nanstd(earlyPiezoEventB2,[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on;
shadedErrorBar(tt,nanmean(latePiezoEventB2,2).*(1000./frameRateHz), (nanstd(latePiezoEventB2,[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
for mouseIdx=1:exptCount
scatter(tempPiezoB2{1,mouseIdx}(:,ind_earlypiezo_B2{1,mouseIdx})*(1000./frameRateHz),zeros(1,length(ind_earlypiezo_B2{1,mouseIdx}))+3,'ob');
scatter(tempPiezoB2{1,mouseIdx}(:,ind_latepiezo_B2{1,mouseIdx})*(1000./frameRateHz),zeros(1,length(ind_latepiezo_B2{1,mouseIdx}))+3,'ok');
end
vline(0,'r');
vline(767,'k');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 inf]);
title(['earliest 25% vs. latest 25% of motion of CS- trials aligned to cue: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_B2)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_B2)).*(1000./frameRateHz),2)) ' ms']);
subplot(2,1,2);
shadedErrorBar(tt,nanmean((t1earlyPiezoB2),2),nanstd((t1earlyPiezoB2),[],2)./sqrt(sum(cell2mat(block2))),'b');
hold on;
shadedErrorBar(tt,nanmean((t1latePiezoB2),2),nanstd((t1latePiezoB2),[],2)./sqrt(sum(cell2mat(block2))),'k');
vline(0,'r');
vline(767,'k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
title(['earliest 25% vs. latest 25% of motion of CS- trials: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_B2)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_B2)).*(1000./frameRateHz),2)) ' ms']);
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | motion bursts of CS- trials: early(blu) & late(blk)']);
if seshType ~= 2
savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatency1SD_B2.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatency1SD_B2.jpeg']);
else
savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatency1SD_B2.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatency1SD_B2.jpeg']);
end
% 2 SD over but under 3 SD of mean piezo voltage
tempFig=figure;
subplot(2,1,1);
shadedErrorBar(tt,nanmean(earlyPiezoEventB2,2).*(1000./frameRateHz), (nanstd(earlyPiezoEventB2,[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on;
shadedErrorBar(tt,nanmean(latePiezoEventB2,2).*(1000./frameRateHz), (nanstd(latePiezoEventB2,[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
for mouseIdx=1:exptCount
scatter(tempPiezoB2{1,mouseIdx}(:,ind_earlypiezo_B2{1,mouseIdx})*(1000./frameRateHz),zeros(1,length(ind_earlypiezo_B2{1,mouseIdx}))+3,'ob');
scatter(tempPiezoB2{1,mouseIdx}(:,ind_latepiezo_B2{1,mouseIdx})*(1000./frameRateHz),zeros(1,length(ind_latepiezo_B2{1,mouseIdx}))+3,'ok');
end
vline(0,'r');
vline(767,'k');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 inf]);
title(['earliest 25% vs. latest 25% of motion of CS- trials aligned to cue: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_B2)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_B2)).*(1000./frameRateHz),2)) ' ms']);
subplot(2,1,2);
shadedErrorBar(tt,nanmean((t2earlyPiezoB2),2),nanstd((t2earlyPiezoB2),[],2)./sqrt(sum(cell2mat(block2))),'b');
hold on;
shadedErrorBar(tt,nanmean((t2latePiezoB2),2),nanstd((t2latePiezoB2),[],2)./sqrt(sum(cell2mat(block2))),'k');
vline(0,'r');
vline(767,'k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
title(['earliest 25% vs. latest 25% of motion of CS- trials: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_B2)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_B2)).*(1000./frameRateHz),2)) ' ms']);
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | motion bursts of CS- trials: early(blu) & late(blk)']);
if seshType ~= 2
savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatency2SD_B2.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatency2SD_B2.jpeg']);
else
savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatency2SD_B2.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatency2SD_B2.jpeg']);
end
%% pulling piezo data without considering frames with motion >=< mean motion
tempFig=figure;
subplot(2,1,1);
shadedErrorBar(tt,nanmean(earlyPiezoEventB2,2).*(1000./frameRateHz), (nanstd(earlyPiezoEventB2,[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on;
shadedErrorBar(tt,nanmean(latePiezoEventB2,2).*(1000./frameRateHz), (nanstd(latePiezoEventB2,[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
for mouseIdx=1:exptCount
scatter(tempPiezoB2{1,mouseIdx}(:,ind_earlypiezo_B2{1,mouseIdx})*(1000./frameRateHz),zeros(1,length(ind_earlypiezo_B2{1,mouseIdx}))+3,'ob');
scatter(tempPiezoB2{1,mouseIdx}(:,ind_latepiezo_B2{1,mouseIdx})*(1000./frameRateHz),zeros(1,length(ind_latepiezo_B2{1,mouseIdx}))+3,'ok');
end
vline(0,'r');
vline(767,'k');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 inf]);
title(['earliest 25% vs. latest 25% of motion of CS- trials aligned to cue: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_B2)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_B2)).*(1000./frameRateHz),2)) ' ms']);
subplot(2,1,2);
shadedErrorBar(tt,nanmean((t1earlyPiezoV2B2),2),nanstd((t1earlyPiezoV2B2),[],2)./sqrt(sum(cell2mat(block2))),'b');
hold on;
shadedErrorBar(tt,nanmean((t1latePiezoV2B2),2),nanstd((t1latePiezoV2B2),[],2)./sqrt(sum(cell2mat(block2))),'k');
vline(0,'r');
vline(767,'k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
title(['earliest 25% vs. latest 25% of motion of CS- trials: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_B2)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_B2)).*(1000./frameRateHz),2)) ' ms']);
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | motion bursts of CS- trials: early(blu) & late(blk)']);
if seshType ~= 2
savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatencyFromPiezoTrace_B2.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatencyFromPiezoTrace_B2.jpeg']);
else
savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatencyFromPiezoTrace_B2.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatencyFromPiezoTrace_B2.jpeg']);
end

tempFig=figure;
subplot(2,1,1);
shadedErrorBar(tt,nanmean(earlyPiezoEventRew,2).*(1000./frameRateHz), (nanstd(earlyPiezoEventRew,[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on;
shadedErrorBar(tt,nanmean(latePiezoEventRew,2).*(1000./frameRateHz), (nanstd(latePiezoEventRew,[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
for mouseIdx=1:exptCount
scatter(tempPiezoRew{1,mouseIdx}(:,ind_earlypiezo_rew{1,mouseIdx})*(1000./frameRateHz),zeros(1,length(ind_earlypiezo_rew{1,mouseIdx}))+3,'ob');
scatter(tempPiezoRew{1,mouseIdx}(:,ind_latepiezo_rew{1,mouseIdx})*(1000./frameRateHz),zeros(1,length(ind_latepiezo_rew{1,mouseIdx}))+3,'ok');
end
vline(0,'g');
vline(767,'k');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 inf]);
title(['earliest 25% vs. latest 25% of motion of CS+ trials aligned to cue: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_rew)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_rew)).*(1000./frameRateHz),2)) ' ms']);
subplot(2,1,2);
shadedErrorBar(tt,nanmean(t1earlyPiezoV2Rew,2),nanstd(t1earlyPiezoV2Rew,[],2)./sqrt(sum(cell2mat(rew))),'b');
hold on;
shadedErrorBar(tt,nanmean(t1latePiezoV2Rew,2),nanstd(t1latePiezoV2Rew,[],2)./sqrt(sum(cell2mat(rew))),'k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
vline(767,'k');
vline(0,'g');
title(['earliest 25% vs. latest 25% of motion of CS+ trials: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_rew)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_rew)).*(1000./frameRateHz),2)) ' ms']);
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | motion bursts of CS+ trials: early(blu) & late(blk)']);
if seshType ~= 2
savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatencyFromPiezoTrace_Rew.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatencyFromPiezoTrace_Rew.jpeg']);
else
savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatencyFromPiezoTrace_Rew.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatencyFromPiezoTrace_Rew.jpeg']);
end

%% CS+
tempFig=figure;
subplot(2,1,1);
shadedErrorBar(tt,nanmean(earlyPiezoEventRew,2).*(1000./frameRateHz), (nanstd(earlyPiezoEventRew,[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on;
shadedErrorBar(tt,nanmean(latePiezoEventRew,2).*(1000./frameRateHz), (nanstd(latePiezoEventRew,[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
for mouseIdx=1:exptCount
scatter(tempPiezoRew{1,mouseIdx}(:,ind_earlypiezo_rew{1,mouseIdx})*(1000./frameRateHz),zeros(1,length(ind_earlypiezo_rew{1,mouseIdx}))+3,'ob');
scatter(tempPiezoRew{1,mouseIdx}(:,ind_latepiezo_rew{1,mouseIdx})*(1000./frameRateHz),zeros(1,length(ind_latepiezo_rew{1,mouseIdx}))+3,'ok');
end
vline(0,'g');
vline(767,'k');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 inf]);
title(['earliest 25% vs. latest 25% of motion of CS+ trials aligned to cue: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_rew)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_rew)).*(1000./frameRateHz),2)) ' ms']);
subplot(2,1,2);
shadedErrorBar(tt,nanmean(t3earlyPiezoRew,2),nanstd(t3earlyPiezoRew,[],2)./sqrt(sum(cell2mat(rew))),'b');
hold on;
shadedErrorBar(tt,nanmean(t3latePiezoRew,2),nanstd(t3latePiezoRew,[],2)./sqrt(sum(cell2mat(rew))),'k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
vline(767,'k');
vline(0,'g');
title(['earliest 25% vs. latest 25% of motion of CS+ trials: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_rew)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_rew)).*(1000./frameRateHz),2)) ' ms']);
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | motion bursts of CS+ trials: early(blu) & late(blk)']);
if seshType ~= 2
savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatency_Rew.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatency_Rew.jpeg']);
else
savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatency_Rew.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatency_Rew.jpeg']);
end

% 1SD-2SD over mean piezo voltage
tempFig=figure;
subplot(2,1,1);
shadedErrorBar(tt,nanmean(earlyPiezoEventRew,2).*(1000./frameRateHz), (nanstd(earlyPiezoEventRew,[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on;
shadedErrorBar(tt,nanmean(latePiezoEventRew,2).*(1000./frameRateHz), (nanstd(latePiezoEventRew,[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
for mouseIdx=1:exptCount
scatter(tempPiezoRew{1,mouseIdx}(:,ind_earlypiezo_rew{1,mouseIdx})*(1000./frameRateHz),zeros(1,length(ind_earlypiezo_rew{1,mouseIdx}))+3,'ob');
scatter(tempPiezoRew{1,mouseIdx}(:,ind_latepiezo_rew{1,mouseIdx})*(1000./frameRateHz),zeros(1,length(ind_latepiezo_rew{1,mouseIdx}))+3,'ok');
end
vline(0,'g');
vline(767,'k');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 inf]);
title(['earliest 25% vs. latest 25% of motion of CS+ trials aligned to cue: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_rew)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_rew)).*(1000./frameRateHz),2)) ' ms']);
subplot(2,1,2);
shadedErrorBar(tt,nanmean(t1earlyPiezoRew,2),nanstd(t1earlyPiezoRew,[],2)./sqrt(sum(cell2mat(rew))),'b');
hold on;
shadedErrorBar(tt,nanmean(t1latePiezoRew,2),nanstd(t1latePiezoRew,[],2)./sqrt(sum(cell2mat(rew))),'k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
vline(767,'k');
vline(0,'g');
title(['earliest 25% vs. latest 25% of motion of CS+ trials: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_rew)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_rew)).*(1000./frameRateHz),2)) ' ms']);
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | motion bursts of CS+ trials: early(blu) & late(blk)']);
if seshType ~= 2
savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatency1SD_Rew.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatency1SD_Rew.jpeg']);
else
savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatency1SD_Rew.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatency1SD_Rew.jpeg']);
end

%2-3SD over mean piezo voltage
tempFig=figure;
subplot(2,1,1);
shadedErrorBar(tt,nanmean(earlyPiezoEventRew,2).*(1000./frameRateHz), (nanstd(earlyPiezoEventRew,[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on;
shadedErrorBar(tt,nanmean(latePiezoEventRew,2).*(1000./frameRateHz), (nanstd(latePiezoEventRew,[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
for mouseIdx=1:exptCount
scatter(tempPiezoRew{1,mouseIdx}(:,ind_earlypiezo_rew{1,mouseIdx})*(1000./frameRateHz),zeros(1,length(ind_earlypiezo_rew{1,mouseIdx}))+3,'ob');
scatter(tempPiezoRew{1,mouseIdx}(:,ind_latepiezo_rew{1,mouseIdx})*(1000./frameRateHz),zeros(1,length(ind_latepiezo_rew{1,mouseIdx}))+3,'ok');
end
vline(0,'g');
vline(767,'k');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 inf]);
title(['earliest 25% vs. latest 25% of motion of CS+ trials aligned to cue: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_rew)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_rew)).*(1000./frameRateHz),2)) ' ms']);
subplot(2,1,2);
shadedErrorBar(tt,nanmean(t2earlyPiezoRew,2),nanstd(t2earlyPiezoRew,[],2)./sqrt(sum(cell2mat(rew))),'b');
hold on;
shadedErrorBar(tt,nanmean(t2latePiezoRew,2),nanstd(t2latePiezoRew,[],2)./sqrt(sum(cell2mat(rew))),'k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
vline(767,'k');
vline(0,'g');
title(['earliest 25% vs. latest 25% of motion of CS+ trials: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_rew)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_rew)).*(1000./frameRateHz),2)) ' ms']);
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | motion bursts of CS+ trials: early(blu) & late(blk)']);
if seshType ~= 2
savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatency2SD_Rew.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatency2SD_Rew.jpeg']);
else
savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatency2SD_Rew.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatency2SD_Rew.jpeg']);
end

%% late B2 vs. Rew
tempFig=figure;
subplot(2,1,1);
shadedErrorBar(tt,nanmean(latePiezoEventRew,2).*(1000./frameRateHz), (nanstd(latePiezoEventRew,[],2).*(1000./frameRateHz))./sqrt(nIC), 'k');
hold on;
shadedErrorBar(tt,nanmean(latePiezoEventB2,2).*(1000./frameRateHz), (nanstd(latePiezoEventB2,[],2).*(1000./frameRateHz))./sqrt(nIC),'r');
for mouseIdx=1:exptCount
scatter(tempPiezoRew{1,mouseIdx}(:,ind_latepiezo_rew{1,mouseIdx})*(1000./frameRateHz),zeros(1,length(ind_latepiezo_rew{1,mouseIdx}))+3,'ok');
scatter(tempPiezoB2{1,mouseIdx}(:,ind_latepiezo_B2{1,mouseIdx})*(1000./frameRateHz),zeros(1,length(ind_latepiezo_B2{1,mouseIdx}))+3,'or');
end
vline(0,'b');
vline(767,'k');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 inf]);
title(['latest 25% of motion of CS+ vs CS- trials aligned to cue: ' num2str(chop(nanmean(cell2mat(HL_piezo.late_rew)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_B2)).*(1000./frameRateHz),2)) ' ms']);
subplot(2,1,2);
shadedErrorBar(tt,nanmean(t3latePiezoRew,2),nanstd(t3latePiezoRew,[],2)./sqrt(sum(cell2mat(rew))),'k');
hold on;
shadedErrorBar(tt,nanmean(t3latePiezoB2,2),nanstd(t3latePiezoB2,[],2)./sqrt(sum(cell2mat(block2))),'r');
vline(0,'b');
vline(767,'k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
title(['latest 25% of motion of CS+ vs CS- trials: ' num2str(chop(nanmean(cell2mat(HL_piezo.late_rew)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_B2)).*(1000./frameRateHz),2)) ' ms']);
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | motion bursts of late trials: CS+(blk) & CS-(red)']);
if seshType ~= 2
savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatency_Late.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatency_Late.jpeg']);
else
savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatency_Late.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatency_Late.jpeg']);
end

tempFig=figure;
subplot(2,1,1);
shadedErrorBar(tt,nanmean(earlyPiezoEventRew,2).*(1000./frameRateHz), (nanstd(earlyPiezoEventRew,[],2).*(1000./frameRateHz))./sqrt(nIC), 'k');
hold on;
shadedErrorBar(tt,nanmean(earlyPiezoEventB2,2).*(1000./frameRateHz), (nanstd(earlyPiezoEventB2,[],2).*(1000./frameRateHz))./sqrt(nIC),'r');
for mouseIdx=1:exptCount
scatter(tempPiezoRew{1,mouseIdx}(:,ind_earlypiezo_rew{1,mouseIdx})*(1000./frameRateHz),zeros(1,length(ind_earlypiezo_rew{1,mouseIdx}))+3,'ok');
scatter(tempPiezoB2{1,mouseIdx}(:,ind_earlypiezo_B2{1,mouseIdx})*(1000./frameRateHz),zeros(1,length(ind_earlypiezo_B2{1,mouseIdx}))+3,'or');
end
vline(0,'b');
vline(767,'k');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 inf]);
title(['earliest 25% of motion of CS+ vs CS- trials aligned to cue: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_rew)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.early_B2)).*(1000./frameRateHz),2)) ' ms']);
subplot(2,1,2);
shadedErrorBar(tt,nanmean(t3earlyPiezoRew,2),nanstd(t3earlyPiezoRew,[],2)./sqrt(sum(cell2mat(rew))),'k');
hold on;
shadedErrorBar(tt,nanmean(t3earlyPiezoB2,2),nanstd(t3earlyPiezoB2,[],2)./sqrt(sum(cell2mat(block2))),'r');
vline(0,'b');
vline(767,'k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
title(['earliest 25% of motion of CS+ vs. CS- trials: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_rew)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.early_B2)).*(1000./frameRateHz),2)) ' ms']);
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | motion bursts of early trials: CS+(blk) & CS-(red)']);

if seshType ~= 2
savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatency_Early.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatency_Early.jpeg']);
else
savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatency_Early.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatency_Early.jpeg']);
end
%% most and least motion
%{
for i = 1:exptCount
    for t = 1:size(groupAlign_piezo,2)
    [sortPiezoRaw{t} sortPiezoRaw_ind{t}] = sort(groupAlign_piezo(:,t),'descend');
    meanPiezoRaw(1,t) = nanmean(sortPiezoRaw{t});
    stdPiezoRaw(1,t) = nanstd(sortPiezoRaw{t});
    maxPiezoRaw(:,t) = sortPiezoRaw{t}(1:5,1);
    for f = 1:150 %always 150 frames per trial
    outliersPiezoRaw(f,t)=sortPiezoRaw{t}(f,1)>=meanPiezoRaw(1,t)+(1);
    end
    end
    c=1; extremeMotionTrial =[]; figure;
    randcolor = {[0.2 0.8 0.8],[0.2 0 0],[1 0.5 0],[0 0.5 1],[0 0.6 0.3],[1 0.2 0.2]};
    for t = 1:size(groupAlign_piezo,2)
        if outliersPiezoRaw(1,t) == 1
            plot(tt,groupAlign_piezo(:,t),'color',randcolor{1,c});
            extremeMotionTrial = [extremeMotionTrial t];
            text((tt(1,sortPiezoRaw_ind{t}(1,1))), sortPiezoRaw{t}(1,1), num2str(t),'color','k');
            hold on;
            c=c+1;
        end
    end
end
%}
%%
nTrial=0; ntrials=[]; earlyPiezoEvent=[]; latePiezoEvent=[]; early25PiezoEvent=[]; late25PiezoEvent=[]; early10PiezoEvent=[]; late10PiezoEvent=[];
t3earlyPiezo=[]; t3latePiezo=[]; t3early25Piezo=[]; t3late25Piezo=[]; t3early10Piezo=[]; t3late10Piezo=[];
t2earlyPiezo=[]; t2latePiezo=[]; t2early25Piezo=[]; t2late25Piezo=[]; t2early10Piezo=[]; t2late10Piezo=[];
t1earlyPiezo=[]; t1latePiezo=[]; t1early25Piezo=[]; t1late25Piezo=[]; t1early10Piezo=[]; t1late10Piezo=[];
for mouseIdx = 1:exptCount
    ntrials = (1:animalTrials(1,mouseIdx));%+nTrial;
    
    tempPiezo{1,mouseIdx} = piezoStart{1,mouseIdx}(1,ntrials);
[sortPiezoStart sortPiezoStart_ind] = sort(tempPiezo{1,mouseIdx},'ascend');
nonzero = find(sortPiezoStart); sortPiezoStart_ind=sortPiezoStart_ind(1,nonzero);
nnan = sum(isnan(tempPiezo{1,mouseIdx}));
ind_earlypiezo_{1,mouseIdx} = sortPiezoStart_ind(1:floor((length(sortPiezoStart_ind)-nnan)/4));% top 25% of trials with early big movement
ind_latepiezo_{1,mouseIdx} = sortPiezoStart_ind((length(sortPiezoStart_ind)-nnan)-floor((length(sortPiezoStart_ind)-nnan)/4)+1:end-nnan);
ind_allearlypiezo25_{1,mouseIdx} = sortPiezoStart_ind(1:floor(length(sortPiezoStart_ind)/4));
ind_alllatepiezo25_{1,mouseIdx} = sortPiezoStart_ind(length(sortPiezoStart_ind)-floor(length(sortPiezoStart_ind)/4)+1:end);
ind_allearlypiezo10_{1,mouseIdx} = sortPiezoStart_ind(1:floor(length(sortPiezoStart_ind)/10)); % top 10% of all trials, regardless of nans. does it make sense? 
ind_alllatepiezo10_{1,mouseIdx} = sortPiezoStart_ind(length(sortPiezoStart_ind)-floor(length(sortPiezoStart_ind)/10)+1:end);
HL_piezo.early_{1,mouseIdx} = (tempPiezo{1,mouseIdx}(:,ind_earlypiezo_{1,mouseIdx}));
HL_piezo.late_{1,mouseIdx} = (tempPiezo{1,mouseIdx}(:,ind_latepiezo_{1,mouseIdx}));
%     nTrial=max(ntrials);
  
    earlyPiezoEvent = [earlyPiezoEvent nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_earlypiezo_{1,mouseIdx}),3)];
    latePiezoEvent = [latePiezoEvent nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_latepiezo_{1,mouseIdx}),3)];
    early25PiezoEvent = [early25PiezoEvent nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_allearlypiezo25_{1,mouseIdx}),3)];
    late25PiezoEvent = [late25PiezoEvent nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_alllatepiezo25_{1,mouseIdx}),3)];
    early10PiezoEvent = [early10PiezoEvent nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_allearlypiezo10_{1,mouseIdx}),3)];
    late10PiezoEvent = [late10PiezoEvent nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_alllatepiezo10_{1,mouseIdx}),3)];

    t1earlyPiezo = [t1earlyPiezo groupAlign_piezo_thresh1{1,mouseIdx}(:,ind_earlypiezo_{1,mouseIdx})];
    t1latePiezo = [t1latePiezo groupAlign_piezo_thresh1{1,mouseIdx}(:,ind_latepiezo_{1,mouseIdx})];
    t1early25Piezo = [t1early25Piezo groupAlign_piezo_thresh1{1,mouseIdx}(:,ind_allearlypiezo25_{1,mouseIdx})];
    t1late25Piezo = [t1late25Piezo groupAlign_piezo_thresh1{1,mouseIdx}(:,ind_alllatepiezo25_{1,mouseIdx})];
    t1early10Piezo = [t1early10Piezo groupAlign_piezo_thresh1{1,mouseIdx}(:,ind_allearlypiezo10_{1,mouseIdx})];
    t1late10Piezo = [t1late10Piezo groupAlign_piezo_thresh1{1,mouseIdx}(:,ind_alllatepiezo10_{1,mouseIdx})];
    
    t2earlyPiezo = [t2earlyPiezo groupAlign_piezo_thresh2{1,mouseIdx}(:,ind_earlypiezo_{1,mouseIdx})];
    t2latePiezo = [t2latePiezo groupAlign_piezo_thresh2{1,mouseIdx}(:,ind_latepiezo_{1,mouseIdx})];
    t2early25Piezo = [t2early25Piezo groupAlign_piezo_thresh2{1,mouseIdx}(:,ind_allearlypiezo25_{1,mouseIdx})];
    t2late25Piezo = [t2late25Piezo groupAlign_piezo_thresh2{1,mouseIdx}(:,ind_alllatepiezo25_{1,mouseIdx})];
    t2early10Piezo = [t2early10Piezo groupAlign_piezo_thresh2{1,mouseIdx}(:,ind_allearlypiezo10_{1,mouseIdx})];
    t2late10Piezo = [t2late10Piezo groupAlign_piezo_thresh2{1,mouseIdx}(:,ind_alllatepiezo10_{1,mouseIdx})];
    
    t3earlyPiezo = [t3earlyPiezo groupAlign_piezo_thresh3{1,mouseIdx}(:,ind_earlypiezo_{1,mouseIdx})];
    t3latePiezo = [t3latePiezo groupAlign_piezo_thresh3{1,mouseIdx}(:,ind_latepiezo_{1,mouseIdx})];
    t3early25Piezo = [t3early25Piezo groupAlign_piezo_thresh3{1,mouseIdx}(:,ind_allearlypiezo25_{1,mouseIdx})];
    t3late25Piezo = [t3late25Piezo groupAlign_piezo_thresh3{1,mouseIdx}(:,ind_alllatepiezo25_{1,mouseIdx})];
    t3early10Piezo = [t3early10Piezo groupAlign_piezo_thresh3{1,mouseIdx}(:,ind_allearlypiezo10_{1,mouseIdx})];
    t3late10Piezo = [t3late10Piezo groupAlign_piezo_thresh3{1,mouseIdx}(:,ind_alllatepiezo10_{1,mouseIdx})];
    
% ntrials = [];

end

tempFig=figure;
subplot(2,1,1);
shadedErrorBar(tt,nanmean(earlyPiezoEvent,2)*(1000./frameRateHz), (nanstd(earlyPiezoEvent,[],2)*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on;
shadedErrorBar(tt,nanmean(latePiezoEvent,2)*(1000./frameRateHz), (nanstd(latePiezoEvent,[],2)*(1000./frameRateHz))./sqrt(nIC),'k');
vline(767,'k');
if cue == 1
    vline(0,'g');
elseif cue == 3
    vline(0,'k');
elseif cue == 2
    vline(0,'r');
end
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 inf]);
title('earliest 25% (blu) vs latest 25% (blk) motion of all trials');

subplot(2,1,2);
shadedErrorBar(tt,nanmean(t3earlyPiezo,2), (nanstd(t3earlyPiezo,[],2))./sqrt(nTrials), 'b');
hold on;
shadedErrorBar(tt,nanmean(t3latePiezo,2), (nanstd(t3latePiezo,[],2))./sqrt(nTrials),'k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
if cue == 1
    vline(0,'g');
elseif cue == 3
    vline(0,'k');
elseif cue == 2
    vline(0,'r');
end
vline(767,'k');
title('earliest 25% vs latest 25% motion of all trials');
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | motion bursts : early(blu) & late(blk)']);
if seshType ~= 2
savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatency_earlyLate.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatency_earlyLate.jpeg']);
else
savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatency_earlyLate.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatency_earlyLate.jpeg']);
end
%
tempFig=figure;
subplot(2,1,1);
shadedErrorBar(tt,nanmean(early25PiezoEvent,2)*(1000./frameRateHz), (nanstd(early25PiezoEvent,[],2)*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on;
shadedErrorBar(tt,nanmean(late25PiezoEvent,2)*(1000./frameRateHz), (nanstd(late25PiezoEvent,[],2)*(1000./frameRateHz))./sqrt(nIC),'k');
vline(767,'k');
if cue == 1
    vline(0,'g');
elseif cue == 3
    vline(0,'k');
elseif cue == 2
    vline(0,'r');
end
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 inf]);
title('earliest 25% (blu) vs latest 25% (blk) motion of all trials');

subplot(2,1,2);
shadedErrorBar(tt,nanmean(t3early25Piezo,2), (nanstd(t3early25Piezo,[],2))./sqrt(nTrials), 'b');
hold on;
shadedErrorBar(tt,nanmean(t3late25Piezo,2), (nanstd(t3late25Piezo,[],2))./sqrt(nTrials),'k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
if cue == 1
    vline(0,'g');
elseif cue == 3
    vline(0,'k');
elseif cue == 2
    vline(0,'r');
end
vline(767,'k');
title('earliest 25% vs latest 25% motion of all trials');
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | motion bursts : early(blu) & late(blk)']);
if seshType ~= 2
savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatency_25.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatency_25.jpeg']);
else
savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatency_25.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatency_25.jpeg']);
end
%{
%this figure plots all trial piezo traces
figure;
for p = 1:size(groupAlign_piezo,2)
    if ~isnan(groupAlign_piezo(1,p))
     randcolor = .38 + (.85-.38) .* rand(1,4);
     plot(tt,groupAlign_piezo(:,p),'color',randcolor);
     hold on;
    else
    end
        plot(tt,nanmean(nanmean(groupAlign_piezo,2),1),'r','LineWidth',1.0);
end
%}
tempFig=figure;
subplot(2,1,1);
shadedErrorBar(tt,nanmean(early10PiezoEvent,2)*(1000./frameRateHz), (nanstd(early10PiezoEvent,[],2)*(1000./frameRateHz))./sqrt(nIC), 'b');
hold on;
shadedErrorBar(tt,nanmean(late10PiezoEvent,2)*(1000./frameRateHz), (nanstd(late10PiezoEvent,[],2)*(1000./frameRateHz))./sqrt(nIC),'k');
vline(767,'k');
if cue == 1
    vline(0,'g');
elseif cue == 3
    vline(0,'k');
elseif cue == 2
    vline(0,'r');
end
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 inf]);
title('earliest 10% (blu) vs latest 10% (blk) motion of all trials');

subplot(2,1,2);
shadedErrorBar(tt,nanmean(t3early10Piezo,2), (nanstd(t3early10Piezo,[],2))./sqrt(nTrials), 'b');
hold on;
shadedErrorBar(tt,nanmean(t3late10Piezo,2), (nanstd(t3late10Piezo,[],2))./sqrt(nTrials),'k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
if cue == 1
    vline(0,'g');
elseif cue == 3
    vline(0,'k');
elseif cue == 2
    vline(0,'r');
end
vline(767,'k');
title('earliest 10% vs latest 10% motion of all trials');
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | motion bursts : early(blu) & late(blk)']);
if seshType ~= 2
savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatency_10.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatency_10.jpeg']);
else
savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatency_10.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatency_10.jpeg']);
end
% save piezo data

%% interleaved figuress
if cue == 3
%% SPIKE RATE SEQUENCE
 tempFig=figure; %reward sequences - spikes
    subplot(3,1,1);
    ylim([0 inf]);
    vline(0,'g');
    vline(767,'k');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(events.rewOneRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.rewOneRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(events.b2OneRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.b2OneRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
   title(['CS+ (blu) [' num2str(sum(cell2mat(rewOneRew))) '] and CS- (blk) [' num2str(sum(cell2mat(block2OneRew))) '] trials following one CS+ trial']);
    subplot(3,1,2);
    ylim([0 inf]);
    vline(0,'g');
    vline(767,'k');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(events.rewTwoRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.rewTwoRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(events.b2TwoRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.b2TwoRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
   title(['CS+ (blu) [' num2str(sum(cell2mat(rewTwoRew))) '] and CS- (blk) [' num2str(sum(cell2mat(block2TwoRew))) '] trials following two CS+ trial']);
    subplot(3,1,3);
    ylim([0 inf]);
    vline(0,'g');
    vline(767,'k');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(events.rewThreeRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.rewThreeRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(events.b2ThreeRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.b2ThreeRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
   title(['CS+ (blu) [' num2str(sum(cell2mat(rewThreeRew))) '] and CS- (blk) [' num2str(sum(cell2mat(block2ThreeRew))) '] trials following three CS+ trial']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | spike rate following repeated CS+ presentation']);
    hold off;
    if seshType ~= 2
    savefig(fullfile(output_fn, '_CSplusHabituation.fig'));
    saveas(tempFig, [output_fn  '_CSplusHabituation.jpeg']);
    else
    savefig(fullfile(output_fn, [blockName '_CSplusHabituation.fig']));
    saveas(tempFig, [output_fn blockName '_CSplusHabituation.jpeg']);
    end
 tempFig=figure; %block 2 sequences - spikes
    subplot(3,1,1);
    ylim([0 inf]);
    vline(0,'r');
    vline(767,'k');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(events.rewOneB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.rewOneB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(events.b2OneB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.b2OneB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
   title(['CS+ (blu) [' num2str(sum(cell2mat(rewOneB2))) '] and CS- (blk) [' num2str(sum(cell2mat(block2OneB2))) '] trials following one CS- trial']);
    subplot(3,1,2);
    ylim([0 inf]);
    vline(0,'r');
    vline(767,'k');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(events.rewTwoB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.rewTwoB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(events.b2TwoB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.b2TwoB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
   title(['CS+ (blu) [' num2str(sum(cell2mat(rewTwoB2))) '] and CS- (blk) [' num2str(sum(cell2mat(block2TwoB2))) '] trials following two CS- trial']);
    subplot(3,1,3);
    ylim([0 inf]);
    vline(0,'r');
    vline(767,'k');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(events.rewThreeB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.rewThreeB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(events.b2ThreeB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.b2ThreeB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
   title(['CS+ (blu) [' num2str(sum(cell2mat(rewThreeB2))) '] and CS- (blk) [' num2str(sum(cell2mat(block2ThreeB2))) '] trials following three CS- trial']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | spike rate following repeated CS- presentation']);
    hold off;
    if seshType ~= 2
    savefig(fullfile(output_fn, '_CSminusHabituation.fig'));
    saveas(tempFig, [output_fn  '_CSminusHabituation.jpeg']);
    else
    savefig(fullfile(output_fn, [blockName '_CSminusHabituation.fig']));
    saveas(tempFig, [output_fn blockName '_CSminusHabituation.jpeg']);
    end

   tempFig=figure; %cue sequences - four in a row - spikes
    subplot(2,1,1);
    ylim([0 inf]);
    vline(0,'r');
    vline(767,'k');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(events.rewFourB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.rewFourB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(events.b2FourB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.b2FourB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
   title(['CS+ (blu) [' num2str(sum(cell2mat(rewFourB2))) '] and CS- (blk) [' num2str(sum(cell2mat(block2FourB2))) '] trials following four CS- trial']);
    subplot(2,1,2);
    ylim([0 inf]);
    vline(0,'g');
    vline(767,'k');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(events.rewFourRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.rewFourRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC), 'b');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(events.b2FourRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.b2FourRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
   title(['CS+ (blu) [' num2str(sum(cell2mat(rewFourRew))) '] and CS- (blk) [' num2str(sum(cell2mat(block2FourRew))) '] trials following four CS+ trial']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | spike rate following repeated CS-/CS+ presentation']);
    hold off;
    if seshType ~= 2
    savefig(fullfile(output_fn, '_fourInARowSpikes.fig'));
    saveas(tempFig, [output_fn  '_fourInARowSpikes.jpeg']);
    else
    savefig(fullfile(output_fn, [blockName '_fourInARowSpikes.fig']));
    saveas(tempFig, [output_fn blockName '_fourInARowSpikes.jpeg']);
    end
%% PIEZO SEQUENCE   
 tempFig=figure; %block 2 sequences - piezo
    subplot(3,1,1);
    vline(0,'r');
    vline(767,'k');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(piezo.rewOneB2GroupAlign),2), (nanstd(cell2mat(events.rewOneB2GroupAlign),[],2))./sqrt(nTrials), 'b');    
    shadedErrorBar(tt,nanmean(cell2mat(piezo.b2OneB2GroupAlign),2), (nanstd(cell2mat(events.b2OneB2GroupAlign),[],2))./sqrt(nTrials),'k');
    xlabel('Time from cue');
    ylabel('Piezo voltage');
    ylim([0 .2]);   
   title(['CS+ (blu) [' num2str(sum(cell2mat(rewOneB2))) '] and CS- (blk) [' num2str(sum(cell2mat(block2OneB2))) '] trials following one CS- trial']);
    subplot(3,1,2);
    vline(0,'r');
    vline(767,'k');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(piezo.rewTwoB2GroupAlign),2), (nanstd(cell2mat(events.rewTwoB2GroupAlign),[],2))./sqrt(nTrials), 'b');
    shadedErrorBar(tt,nanmean(cell2mat(piezo.b2TwoB2GroupAlign),2), (nanstd(cell2mat(events.b2TwoB2GroupAlign),[],2))./sqrt(nTrials),'k');
    xlabel('Time from cue');
    ylabel('Piezo voltage');
    ylim([0 .2]);   
   title(['CS+ (blu) [' num2str(sum(cell2mat(rewTwoB2))) '] and CS- (blk) [' num2str(sum(cell2mat(block2TwoB2))) '] trials following two CS- trial']);
    subplot(3,1,3);
    vline(0,'r');
    vline(767,'k');
    hold on;    
    shadedErrorBar(tt,nanmean(cell2mat(piezo.rewThreeB2GroupAlign),2), (nanstd(cell2mat(events.rewThreeB2GroupAlign),[],2))./sqrt(nTrials), 'b');
    shadedErrorBar(tt,nanmean(cell2mat(piezo.b2ThreeB2GroupAlign),2), (nanstd(cell2mat(events.b2ThreeB2GroupAlign),[],2))./sqrt(nTrials),'k');
    xlabel('Time from cue');
    ylabel('Piezo voltage');
    ylim([0 .2]); 
   title(['CS+ (blu) [' num2str(sum(cell2mat(rewThreeB2))) '] and CS- (blk) [' num2str(sum(cell2mat(block2ThreeB2))) '] trials following three CS- trial']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | spike rate following repeated CS- presentation']);
    hold off;
    if seshType ~= 2
    savefig(fullfile(output_fn, '_CSminusPiezo.fig'));
    saveas(tempFig, [output_fn  '_CSminusPiezo.jpeg']);
    else
    savefig(fullfile(output_fn, [blockName '_CSminusPiezo.fig']));
    saveas(tempFig, [output_fn blockName '_CSminusPiezo.jpeg']);
    end
    
 tempFig=figure; %reward sequences - piezo
    subplot(3,1,1);
    vline(0,'g');
    vline(767,'k');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(piezo.rewOneRewGroupAlign),2), (nanstd(cell2mat(events.rewOneRewGroupAlign),[],2))./sqrt(nTrials), 'b');    
    shadedErrorBar(tt,nanmean(cell2mat(piezo.b2OneRewGroupAlign),2), (nanstd(cell2mat(events.b2OneRewGroupAlign),[],2))./sqrt(nTrials),'k');
    xlabel('Time from cue');
    ylabel('Piezo voltage');
    ylim([0 inf]);   
   title(['CS+ (blu) [' num2str(sum(cell2mat(rewOneRew))) '] and CS- (blk) [' num2str(sum(cell2mat(block2OneRew))) '] trials following one CS+ trial']);
    subplot(3,1,2);
    vline(0,'g');
    vline(767,'k');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(piezo.rewTwoRewGroupAlign),2), (nanstd(cell2mat(events.rewTwoRewGroupAlign),[],2))./sqrt(nTrials), 'b');
    shadedErrorBar(tt,nanmean(cell2mat(piezo.b2TwoRewGroupAlign),2), (nanstd(cell2mat(events.b2TwoRewGroupAlign),[],2))./sqrt(nTrials),'k');
    xlabel('Time from cue');
    ylabel('Piezo voltage');
    ylim([0 inf]);    
   title(['CS+ (blu) [' num2str(sum(cell2mat(rewTwoRew))) '] and CS- (blk) [' num2str(sum(cell2mat(block2TwoRew))) '] trials following two CS+ trial']);
    subplot(3,1,3);
    vline(0,'g');
    vline(767,'k');
    hold on;    
    shadedErrorBar(tt,nanmean(cell2mat(piezo.rewThreeRewGroupAlign),2), (nanstd(cell2mat(events.rewThreeRewGroupAlign),[],2))./sqrt(nTrials), 'b');
    shadedErrorBar(tt,nanmean(cell2mat(piezo.b2ThreeRewGroupAlign),2), (nanstd(cell2mat(events.b2ThreeRewGroupAlign),[],2))./sqrt(nTrials),'k');
    xlabel('Time from cue');
    ylabel('Piezo voltage');
    ylim([0 inf]);     
   title(['CS+ (blu) [' num2str(sum(cell2mat(rewThreeRew))) '] and CS- (blk) [' num2str(sum(cell2mat(block2ThreeRew))) '] trials following three CS+ trial']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | spike rate following repeated CS+ presentation']);
    hold off;
    if seshType ~= 2
    savefig(fullfile(output_fn, '_CSplusPiezo.fig'));
    saveas(tempFig, [output_fn  '_CSplusPiezo.jpeg']);
    else
    savefig(fullfile(output_fn, [blockName '_CSplusPiezo.fig']));
    saveas(tempFig, [output_fn blockName '_CSplusPiezo.jpeg']);
    end
%% LICK SEQUENCE

 tempFig=figure; % still plotting neural data, but seperate the trials based on licking rate of that trial
    subplot(3,1,1);
    shadedErrorBar(tt,nanmean(cell2mat(lick.rewOneRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(lick.rewOneRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(cell2mat(lick.rewOneRewGroupAlign))))), 'b');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(lick.b2OneRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(lick.b2OneRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(cell2mat(lick.b2OneRewGroupAlign))))),'k');
    xlabel('Time from cue');
    ylabel('Lick rate (Hz)');
    ylim([0 inf]);
    title(['CS+ (blu) [' num2str(sum(cell2mat(rewOneRew))) '] and CS- (blk) [' num2str(sum(cell2mat(block2OneRew))) '] trials following one CS+ trial']);
    subplot(3,1,2);
    shadedErrorBar(tt,nanmean(cell2mat(lick.rewTwoRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(lick.rewTwoRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(cell2mat(lick.rewTwoRewGroupAlign))))), 'b');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(lick.b2TwoRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(lick.b2TwoRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(cell2mat(lick.b2TwoRewGroupAlign))))),'k');
    xlabel('Time from cue');
    ylabel('Lick rate (Hz)');
    ylim([0 inf]);
    title(['CS+ (blu) [' num2str(sum(cell2mat(rewTwoRew))) '] and CS- (blk) [' num2str(sum(cell2mat(block2TwoRew))) '] trials following two CS+ trial']);
    subplot(3,1,3);
    shadedErrorBar(tt,nanmean(cell2mat(lick.rewThreeRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(lick.rewThreeRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(cell2mat(lick.rewThreeRewGroupAlign))))), 'b');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(lick.b2ThreeRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(lick.b2ThreeRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(cell2mat(lick.b2ThreeRewGroupAlign))))),'k');
    xlabel('Time from cue');
    ylabel('Lick rate (Hz)');
    ylim([0 inf]);
    title(['CS+ (blu) [' num2str(sum(cell2mat(rewThreeRew))) '] and CS- (blk) [' num2str(sum(cell2mat(block2ThreeRew))) '] trials following three CS+ trial']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | lick bursts based on prior stimuli presentation']);
    hold off;
    if seshType ~= 2
    savefig(fullfile(output_fn, '_LickRateRewRepetition.fig'));
    saveas(tempFig, [output_fn  '_LickRateRewRepetition.jpeg']);
    else
    savefig(fullfile(output_fn, [blockName '_LickRateRewRepetition.fig']));
    saveas(tempFig, [output_fn blockName '_LickRateRewRepetition.jpeg']);
    end
    
    tempFig=figure; % still plotting neural data, but seperate the trials based on licking rate of that trial
    subplot(3,1,1);
    shadedErrorBar(tt,nanmean(cell2mat(lick.rewOneB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(lick.rewOneB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(cell2mat(lick.rewOneB2GroupAlign))))), 'b');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(lick.b2OneB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(lick.b2OneB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(cell2mat(lick.b2OneB2GroupAlign))))),'k');
    xlabel('Time from cue');
    ylabel('Lick rate (Hz)');
    ylim([0 inf]);
    title(['CS+ (blu) [' num2str(sum(cell2mat(rewOneB2))) '] and CS- (blk) [' num2str(sum(cell2mat(block2OneB2))) '] trials following one CS- trial']);
    subplot(3,1,2);
    shadedErrorBar(tt,nanmean(cell2mat(lick.rewTwoB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(lick.rewTwoB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(cell2mat(lick.rewTwoB2GroupAlign))))), 'b');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(lick.b2TwoB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(lick.b2TwoB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(cell2mat(lick.b2TwoB2GroupAlign))))),'k');
    xlabel('Time from cue');
    ylabel('Lick rate (Hz)');
    ylim([0 inf]);
    title(['CS+ (blu) [' num2str(sum(cell2mat(rewTwoB2))) '] and CS- (blk) [' num2str(sum(cell2mat(block2TwoB2))) '] trials following two CS- trial']);
    subplot(3,1,3);
    shadedErrorBar(tt,nanmean(cell2mat(lick.rewThreeB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(lick.rewThreeB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(cell2mat(lick.rewThreeB2GroupAlign))))), 'b');
    hold on;
    shadedErrorBar(tt,nanmean(cell2mat(lick.b2ThreeB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(lick.b2ThreeB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(cell2mat(lick.b2ThreeB2GroupAlign))))),'k');
    xlabel('Time from cue');
    ylabel('Lick rate (Hz)');
    ylim([0 inf]);
    title(['CS+ (blu) [' num2str(sum(cell2mat(rewThreeB2))) '] and CS- (blk) [' num2str(sum(cell2mat(block2ThreeB2))) '] trials following three CS- trial']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | lick bursts based on prior stimuli presentation']);
    hold off;
    if seshType ~= 2
    savefig(fullfile(output_fn, '_LickRateB2Repetition.fig'));
    saveas(tempFig, [output_fn  '_LickRateB2Repetition.jpeg']);
    else
    savefig(fullfile(output_fn, [blockName '_LickRateB2Repetition.fig']));
    saveas(tempFig, [output_fn blockName '_LickRateB2Repetition.jpeg']);
    end
end
    %%
    
     s=1;
%         if length(ind_omit>5); s = s+1; end
%         if length(ind_unexp>5); s = s+1; end
%         if input.doBlock2
            s = s+1; 
%         end
        
        if unique(expt{1,1}.name == 'Block')
            indivFig = true;
        else
            indivFig = false;
        end
            
        b2GroupTC=[]; rewGroupTC=[]; targetAlignFB2=[]; targetAlignFRew=[];
        eventsB2=[]; eventsRew=[]; ntrials=0;
        for mouseIdx = 1:exptCount
               nTrial = (1:animalTrials(1,mouseIdx))+ntrials;  
b2GroupTC = [b2GroupTC nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2{1,mouseIdx})),3)];
rewGroupTC = [rewGroupTC nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rew{1,mouseIdx})),3)];

targetAlignFB2 = [targetAlignFB2 nanmean(groupAlign.tc{1,mouseIdx}(1:prewin_frames,:,logical(block2{1,mouseIdx})),3)]; % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.
targetAlignFRew = [targetAlignFRew nanmean(groupAlign.tc{1,mouseIdx}(1:prewin_frames,:,logical(rew{1,mouseIdx})),3)]; % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.

eventsB2 = [eventsB2 nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(block2{1,mouseIdx})),3)];
eventsRew = [eventsRew nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(rew{1,mouseIdx})),3)];
               ntrials=max(nTrial);
        end
        
    dFoFB2 = zeros(size(b2GroupTC,1),size(b2GroupTC,2),size(b2GroupTC,3));%frame*cell*trials
    dFoFRew = zeros(size(rewGroupTC,1),size(rewGroupTC,2),size(rewGroupTC,3));%frame*cell*trials
    % calculate df/F
    for c = 1:sum(nsize)
        dFoFB2(:,c) = (b2GroupTC(:,c)-targetAlignFB2(c))./targetAlignFB2(c);
        dFoFRew(:,c) = (rewGroupTC(:,c)-targetAlignFRew(c))./targetAlignFRew(c);
    end
    
        if ~indivFig 
        tempFig=figure; n = 1;
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(dFoFRew,2), nanstd(dFoFRew,[],2)./sqrt(nIC), 'k');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        %ylim([-3 5]);
        title('CS+')
        n = n+1;
%         if mworks.doBlock2
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(dFoFB2,2), nanstd(dFoFB2,[],2)./sqrt(nIC),'r');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        %ylim([-3 5]);
        title('CS-')
%         end
        sgtitle(['dFoF across animals | ' num2str(sum(nsize)) ' neurons from ' num2str(exptCount) ' animals']);
        if seshType ~= 2
        savefig(fullfile(output_fn, '_cueAlign_dFoverF.fig'))
        saveas(tempFig, [output_fn  '_cueAlign_dFoverF.jpeg']);
        else
        savefig(fullfile(output_fn, [blockName '_cueAlign_dFoverF.fig']));
        saveas(tempFig, [output_fn blockName '_cueAlign_dFoverF.jpeg']);
        end
        
        tempFig=figure; n=1;
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(eventsRew,2).*(1000./frameRateHz), (nanstd(eventsRew,[],2)./sqrt(nIC)).*(1000./frameRateHz));
        vline(767)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([0 inf]);
        title('CS+')
        n = n+1;
%         if mworks.doBlock2
        subplot(s,1,n)
        shadedErrorBar(tt, nanmean(eventsB2,2).*(1000./frameRateHz), (nanstd(eventsB2,[],2)./sqrt(nIC)).*(1000./frameRateHz),'r');
        vline(767)
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([0 inf]);
        title('CS-')
%         end
        sgtitle(['spike rate across animals | ' num2str(sum(nsize)) ' neurons from ' num2str(exptCount) ' animals'])
        if seshType ~= 2
        savefig(fullfile(output_fn, '_cueAlign_events_Hz.fig'))
        saveas(tempFig, [output_fn  '_cueAlign_events_Hz.jpeg']);
        else
        savefig(fullfile(output_fn, [blockName '_cueAlign_events_Hz.fig']));
        saveas(tempFig, [output_fn blockName '_cueAlign_events_Hz.jpeg']);
        end
      %%  
        elseif indivFig
            
        tempFig=figure; %n = 1;
        %subplot(s,1,n) 
        if expt{1,1}.whichBlock(1,isesh) == 0
            shadedErrorBar(tt, nanmean(dFoFB2,2), nanstd(dFoFB2,[],2)./sqrt(nIC),'red');
            vline(0,'r')
            vline(767,'k')        
            xlabel('Time from cue')
            ylabel('Avg dF/F')
        elseif expt{1,1}.whichBlock(1,isesh) == 1
            shadedErrorBar(tt, nanmean(dFoFRew,2), nanstd(dFoFRew,[],2)./sqrt(nIC),'black');
            vline(0, 'g')
            vline(767,'k')        
            xlabel('Time from cue')
            ylabel('Avg dF/F')
        end
        sgtitle(['dF/F of block sessions: ' blockName])
        if seshType ~= 2
        savefig(fullfile(output_fn, '_cueAlign_dFoverF.fig'))
        saveas(tempFig, [output_fn  '_cueAlign_dFoverF.jpeg']);
        else
        savefig(fullfile(output_fn, [blockName '_cueAlign_dFoverF.fig']));
        saveas(tempFig, [output_fn blockName '_cueAlign_dFoverF.jpeg']);
        end
        
        tempFig=figure; %n=1;
        %subplot(s,1,n)
        if expt{1,1}.whichBlock(1,isesh) == 0
            shadedErrorBar(tt, nanmean(eventsB2,2).*(1000./frameRateHz), (nanstd(eventsB2,[],2)./sqrt(nIC)).*(1000./frameRateHz),'red');
            ylim([0 3]);
            vline(0,'r')
            vline(767,'k')
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
        elseif expt{1,1}.whichBlock(1,isesh) == 1
            shadedErrorBar(tt, nanmean(eventsRew,2).*(1000./frameRateHz), (nanstd(eventsRew,[],2)./sqrt(nIC)).*(1000./frameRateHz),'black');
            ylim([0 3]);
            vline(0,'g')
            vline(767,'k')
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
        end
        sgtitle(['spike rate for blocked sessions: ' blockName])
        if seshType ~= 2
        savefig(fullfile(output_fn, '_cueAlign_events_Hz.fig'))
        saveas(tempFig, [output_fn  '_cueAlign_events_Hz.jpeg']);
        else
        savefig(fullfile(output_fn, [blockName '_cueAlign_events_Hz.fig']));
        saveas(tempFig, [output_fn blockName '_cueAlign_events_Hz.jpeg']);
        end
        end
       %% 
      
        
        if unique(expt{1,1}.name == 'Inter')
         tempFig=figure; n=1;
        %subplot(s,1,n)
        plot(tt, nanmean(eventsRew,2).*(1000./frameRateHz),'black');
        hold on;
        plot(tt, nanmean(eventsB2,2).*(1000./frameRateHz),'r');
        ylim([0 3]);
        vline(0,'k')
        vline(767,'k')
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        %ylim([0 6]);
        legend('CS+','CS-');
        sgtitle(['spike rate across animals | ' num2str(sum(nsize)) ' neurons from ' num2str(exptCount) ' animals'])
        if seshType ~= 2
        savefig(fullfile(output_fn, '_line_stackedSpike_events_Hz.fig'))
        saveas(tempFig, [output_fn  '_line_stackedSpike_events_Hz.jpeg']);
        else
        savefig(fullfile(output_fn, [blockName '_line_stackedSpike_events_Hz.fig']));
        saveas(tempFig, [output_fn blockName '_line_stackedSpike_events_Hz.jpeg']);
        end
        
         tempFig=figure; n=1;
        %subplot(s,1,n)
        shadedErrorBar(tt, nanmean(eventsRew,2).*(1000./frameRateHz), (nanstd(eventsRew,[],2)./sqrt(nIC)).*(1000./frameRateHz),'k');
        hold on;
        shadedErrorBar(tt, nanmean(eventsB2,2).*(1000./frameRateHz), (nanstd(eventsB2,[],2)./sqrt(nIC)).*(1000./frameRateHz),'r');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([0 3]);
        %legend('CS+','CS-');  
        vline(767,'k')
        vline(0,'k')
        sgtitle(['stacked spike rate across animals | ' num2str(sum(nsize)) ' neurons from ' num2str(exptCount) ' animals'])
        hold off;
        if seshType ~= 2
        savefig(fullfile(output_fn, '_error_stackedSpike_events_Hz.fig'))
        saveas(tempFig, [output_fn  '_error_stackedSpike_events_Hz.jpeg']);
        else
        savefig(fullfile(output_fn, [blockName '_error_stackedSpike_events_Hz.fig']));
        saveas(tempFig, [output_fn blockName '_error_stackedSpike_events_Hz.jpeg']);
        end
        end 
        %
                   
        if unique(expt{1,1}.name == 'Inter')
            z = [0 20 40]; for m = 1:3; indTemp{1,m} = (1:20)+z(1,m); end
        for mouseIdx=1:(exptCount)  
            logEventsB2{1,mouseIdx}(:,:,:) = groupAlign.events{1,mouseIdx}(:,:,logical(block2{1,mouseIdx}));
            logEventsRew{1,mouseIdx}(:,:,:) = groupAlign.events{1,mouseIdx}(:,:,logical(rew{1,mouseIdx}));
        end
        for mouseIdx=1:exptCount
tempEventsRew{1,mouseIdx} = nan(150,nsize(1,mouseIdx),60);
tempEventsB2{1,mouseIdx} = nan(150,nsize(1,mouseIdx),60); 
         if size(logEventsB2{1,mouseIdx},3) < 60
                tTr = size(logEventsB2{1,mouseIdx}(:,:,:),3);
                tempEventsB2{1,mouseIdx}(:,:,1:tTr) = logEventsB2{1,mouseIdx}(:,:,1:end);                
         else
                tempEventsB2{1,mouseIdx}(:,:,:) = logEventsB2{1,mouseIdx}(:,:,1:60);
         end 
         if size(logEventsRew{1,mouseIdx},3) < 60
                tTr = size(logEventsRew{1,mouseIdx}(:,:,:),3);
                tempEventsRew{1,mouseIdx}(:,:,1:tTr) = logEventsRew{1,mouseIdx}(:,:,1:end);                
         else
                tempEventsRew{1,mouseIdx}(:,:,:) = logEventsRew{1,mouseIdx}(:,:,1:60);
         end 
        end   
        
        fMLTrialRew{1,1} = []; fMLTrialRew{1,2} = []; fMLTrialRew{1,3} = [];
        fMLTrialB2{1,1} = []; fMLTrialB2{1,2} = []; fMLTrialB2{1,3} = [];
        for m = 1:3
            for i = 1:exptCount
            fMLTrialRew{1,m} = [fMLTrialRew{1,m} tempEventsRew{1,i}(:,:,indTemp{1,m})];
            fMLTrialB2{1,m} = [fMLTrialB2{1,m} tempEventsB2{1,i}(:,:,indTemp{1,m})];
            end
        end
        fMLTrialB2{1,1} = nanmean(fMLTrialB2{1,1},3);        fMLTrialRew{1,1} = nanmean(fMLTrialRew{1,1},3);
        fMLTrialB2{1,2} = nanmean(fMLTrialB2{1,2},3);        fMLTrialRew{1,2} = nanmean(fMLTrialRew{1,2},3);
        fMLTrialB2{1,3} = nanmean(fMLTrialB2{1,3},3);        fMLTrialRew{1,3} = nanmean(fMLTrialRew{1,3},3);
       
        tempFig=figure;
        n = floor(nTrials./3);
        start = 1;
        for i = 1:3
            subplot(3,2,start)
            ind_rew_temp = intersect(cell2mat(rew),1+(i-1).*n:i*n);
            shadedErrorBar(tt,nanmean(fMLTrialRew{1,i},2).*(1000./frameRateHz), (nanstd(fMLTrialRew{1,i},[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
            hold on
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([0 3]);
            title(['Trials ' num2str(indTemp{1,i}(1,1)) ':' num2str(indTemp{1,i}(1,end))])
            vline(767,'k')
            vline(0,'k')
            if length(cell2mat(block2))>=10
                subplot(3,2,start+1)
                ind_block2_temp = intersect(cell2mat(block2),1+(i-1).*n:i*n);
                shadedErrorBar(tt,nanmean(fMLTrialB2{1,i},2).*(1000./frameRateHz), (nanstd(fMLTrialB2{1,i},[],2).*(1000./frameRateHz))./sqrt(nIC),'r');
                ylim([0 3])
                xlabel('Time from cue')
                ylabel('Spike rate (Hz)')
                title(['Trials ' num2str(indTemp{1,i}(1,1)) ':' num2str(indTemp{1,i}(1,end))])
                vline(767,'k')
                vline(0,'k')
            end
            start = start+2;
        end
        sgtitle([num2str(sum(nsize)) 'neurons | CS+ (black), CS- (red)'])
        if seshType ~= 2
        savefig(fullfile(output_fn, '_repsByTrial.fig'))
        saveas(tempFig, [output_fn  '_repsByTrial.jpeg']);
        else
        savefig(fullfile(output_fn, [blockName '_repsByTrial.fig']));
        saveas(tempFig, [output_fn blockName '_repsByTrial.jpeg']);
        end
        
        elseif unique(expt{1,1}.name == 'Block')
            
z = [0 20 40]; for m = 1:3; indTemp{1,m} = (1:20)+z(1,m); end
        for mouseIdx=1:(exptCount)  
            tempEvents{1,mouseIdx} = nan(150,nsize(1,mouseIdx),60);
            if size(groupAlign.events{1,mouseIdx}(:,:,:),3) < 60
                tTr = size(groupAlign.events{1,mouseIdx}(:,:,:),3);
                tempEvents{1,mouseIdx}(:,:,1:tTr) = groupAlign.events{1,mouseIdx}(:,:,1:end);                
            elseif ~isnan(nanmean(nanmean(eventsB2(:,:),2),1))
                tempEvents{1,mouseIdx}(:,:,:) = groupAlign.events{1,mouseIdx}(:,:,1:60);
            elseif isnan(nanmean(nanmean(eventsB2(:,:),2),1))
                tempEvents{1,mouseIdx}(:,:,:) = groupAlign.events{1,mouseIdx}(:,:,1:60);
            end
        end
        fMLTrial{1,1} = []; fMLTrial{1,2} = []; fMLTrial{1,3} = [];
        for m = 1:3
            for i = 1:length(tempEvents)
            fMLTrial{1,m} = [fMLTrial{1,m} tempEvents{1,i}(:,:,indTemp{1,m})];
            end
        end
        fMLTrial{1,1} = nanmean(fMLTrial{1,1},3);
        fMLTrial{1,2} = nanmean(fMLTrial{1,2},3);
        fMLTrial{1,3} = nanmean(fMLTrial{1,3},3);
        %
            if ~isnan(nanmean(nanmean(eventsB2(:,:),2),1))
        tempFig=figure;
        n = floor(nTrials./3);
        start = 1;
        for i = 1:3
            subplot(3,1,start)
            ind_rew_temp = intersect(cell2mat(block2),1+(i-1).*n:i*n);
            shadedErrorBar(tt,nanmean(fMLTrial{1,i},2).*(1000./frameRateHz), (nanstd(fMLTrial{1,i},[],2).*(1000./frameRateHz))./sqrt(nIC),'r');
            hold on
            ylim([0 3])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Trials ' num2str(indTemp{1,i}(1,1)) ':' num2str(indTemp{1,i}(1,end))])
            vline(767,'k')
            vline(0,'r')
            start=start+1;
        end 
        sgtitle([num2str(sum(nsize)) 'neurons | CS- (red)'])
        if seshType ~= 2
        savefig(fullfile(output_fn, '_minusrepsByTrial.fig'))
        saveas(tempFig, [output_fn  '_minusrepsByTrial.jpeg']);
        else
        savefig(fullfile(output_fn, [blockName '_minusrepsByTrial.fig']));
        saveas(tempFig, [output_fn blockName '_minusrepsByTrial.jpeg']);
        end
            elseif isnan(nanmean(nanmean(eventsB2(:,:),2),1))
        tempFig=figure;
        n = floor(nTrials./3);
        start = 1;
        for i = 1:3
            subplot(3,1,start)
            ind_rew_temp = intersect(cell2mat(rew),1+(i-1).*n:i*n);
            shadedErrorBar(tt,nanmean(fMLTrial{1,i},2).*(1000./frameRateHz), (nanstd(fMLTrial{1,i},[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
            hold on
            ylim([0 3])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Trials ' num2str(indTemp{1,i}(1,1)) ':' num2str(indTemp{1,i}(1,end))])
            vline(767,'k')
            vline(0,'g')
            start=start+1;
        end
        sgtitle([num2str(sum(nsize)) 'neurons | CS+ (black)'])
        if seshType ~= 2
        savefig(fullfile(output_fn, '_plusrepsByTrial.fig'))
        saveas(tempFig, [output_fn  '_plusrepsByTrial.jpeg']);
        else
        savefig(fullfile(output_fn, [blockName '_plusrepsByTrial.fig']));
        saveas(tempFig, [output_fn blockName '_plusrepsByTrial.jpeg']);
        end
            end
        end
        
        
        
        