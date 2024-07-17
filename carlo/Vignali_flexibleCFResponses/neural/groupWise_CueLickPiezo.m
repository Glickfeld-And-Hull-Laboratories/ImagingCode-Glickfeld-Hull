clear
close all

analysis_out = 'A:\home\carlo\RC\analysis\2P\';
% analysis_out = '/Volumes/All_staff/home/carlo/analysis/2P/';

bdata_source = 'A:\home\carlo\RC\rawData\behavioral\';
% bdata_source = '/Volumes/All_staff/home/carlo/rawData/behavioral/';
seshTitle={'Naive';'PL';'RuleFlipNaive';'RuleFlipPL'};
expType=input('monomodal(1), dimodal(2), or CS+ only (3): '); expIdx = {'monomodal', 'dimodal', 'csPlus'};
seshType=input('naive (1), blocked (2), post-learning (3), or rule change naive(4)or PL(5) session: ');
if seshType==3
    if expType==1; RCExptListAll_Inter; elseif expType==2; RCExptListCarlo_VA_Inter; elseif expType==3; RCExptListCarlo_CSPlus_Inter; end
seshID=2;
if expType==3; whichPL = input('Which post-learning session, monostimulus(1) or interleaved(2)?: '); 
    if whichPL==1; sesh = 'PL_csPlusOnly'; elseif whichPL==2; sesh='PL_interleaved'; end
else
sesh = 'PL'; sesh2=expIdx{1,expType};
end
elseif seshType==2
    if expType==1; RCExptListCarlo_Block; elseif expType==2; RCExptListCarlo_VA_Block; elseif expType==3; RCExptListCarlo_CSPlus_Block; end
sesh = 'Block'; sesh2=expIdx{1,expType};
elseif seshType==1
seshID=1;
    if expType==1; RCExptListAll_Naive; elseif expType==2; RCExptListCarlo_VA_Inter; elseif expType==3; RCExptListCarlo_CSPlus_Inter; end
sesh = 'Naive'; sesh2=expIdx{1,expType};
elseif seshType==4
seshID=3;    
     if expType==1; disp('This does not exist'); return; elseif expType==2; RCExptListCarlo_VA_Inter; end
 sesh = 'RuleFlipNaive'; sesh2=expIdx{1,expType};
elseif seshType==5
seshID=4;    
     if expType==1; disp('This does not exist'); return; elseif expType==2; RCExptListCarlo_VA_Inter; end
 sesh = 'RuleFlipPL'; sesh2=expIdx{1,expType};
end
% it's gonna have to be seperated again, so if animal ID > 1700 go to
% /analysis/2P if animal ID < 1700 go to /mikeAnalysis/2P
%% pull out all neural data across animals based on ID (<1600 vs. >1700)
thresholdDeco = input('Which threshold do you wish to plot? ("-#" for decon | "#.#" for first der.): ');
analysis_out = 'A:\home\carlo\RC\analysis\2P\';
output_fn = ([analysis_out 'groupwiseAlign' num2str(thresholdDeco) '_' sesh '_' sesh2 '\']);
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
analysis_out = 'A:\home\carlo\RC\analysis\2P\';
output_fn = ([analysis_out 'groupwiseAlign' num2str(thresholdDeco) '_' sesh '_' sesh2 '\']);


%% call blanks for bx data concat
% gpreRew_lickBurstHz = []; gpostRew_lickBurstHz = []; glickBurstHz_all = []; gpostRew_lickBurstStart = []; gpostRew_lickAlignEvents = []; gpostRew_lickAlign = []; gpostRewTrials = []; gpreRewTrials = []; gfirstPostRew_lickAlignEvents = []; gfirstPostRew_lickAlign = []; glastPreRew_lickAlignEvents = []; glastPreRew_lickAlign = []; gTargetOn = []; gRctT = [];  gind_block2 = []; gind_rew = []; ggoodcells = [];
gind_block2 = []; gind_rew = []; gTargetOn = []; gRctT = []; initFrames = []; ggoodcells = [];gtt = [];
initLickTimes = []; initCounterTimes = []; initCounterVals = [];
ind_early_bst=[]; ind_late_bst=[]; early_bsttime=[]; late_bsttime=[];
earlyBurstEvents=[]; lateBurstEvents=[]; earlyBurstStart=[]; lateBurstStart=[];
earlyRewBurstEvents=[]; lateRewBurstEvents=[]; earlyRewBurstStart=[]; lateRewBurstStart=[];
earlyB2BurstEvents=[]; lateB2BurstEvents=[]; earlyB2BurstStart=[]; lateB2BurstStart=[];
tempearlyCSPLick_ind=[]; templateCSPLick_ind=[]; tempearlyCSMLick_ind=[]; templateCPMLick_ind=[];
%% 
if seshType == 2
allBlock = input('all block sessions (1-yes): ');
pmpt = input('CS-? [1-yes|2-no]: ');
if allBlock == 1
    exptCount = 2.*length(expt);
else
    exptCount = length(expt);
end
else
    exptCount = length(expt);
end


for mouseIdx = 1:exptCount %group all mouse data
analysis_out = 'A:\home\carlo\RC\analysis\2P\';
[pairings] = loadRCList_Carlo(expt,seshType);
    %1512, 1513, 1516, 1520, 1081, 1702, 1705, 1706 
    % session number always 1 -- seperate PL and Naive lists
iexp = pairings{1,1}; isesh = pairings{3,1};
mouse = strtrim(expt{1,iexp}.mouse{1,isesh});
date = expt{1,iexp}.date{1,isesh};
run = expt{1,iexp}.run{1,isesh};
sessions = [date '_img' mouse];
    fprintf([date ' ' mouse ' ' run '\n']);

if isempty(thresholdDeco)
img_fn = [date '_img' mouse '\getTC_' run '\'];
elseif thresholdDeco<0 %deconvolve data
img_fn = [date '_img' mouse '\getTC_' run '_' num2str(thresholdDeco) '\'];
elseif thresholdDeco>0 %first derivative data
img_fn = [date '_img' mouse '\getTC_' run '_FD' num2str(thresholdDeco) '\'];        
end
%%
if str2double(regexp(mouse,'\d*','match')) < 1700 && str2double(regexp(mouse,'\d*','match')) ~= 1081
analysis_out = 'A:\home\carlo\RC\mikeAnalysis\2P\';
bdata_source = 'A:\home\carlo\RC\rawData\behavioral\';%mike\Data\Behavior\';  
elseif str2double(regexp(mouse,'\d*','match')) > 1700 || str2double(regexp(mouse,'\d*','match')) == 1081 
analysis_out = 'A:\home\carlo\RC\analysis\2P\';
bdata_source = 'A:\home\carlo\RC\rawData\behavioral\';
end
%% call data Bx and Neural
    mworks = getBxData_Carlo(bdata_source, sessions);
    cTargetOn = mworks.cTargetOn;
    if iscell(cTargetOn) % if it is a cell, it means cTargetOn wasn't being over written in extract TC. If it's not a cell, it's already over written in extract TC. can be used directly
        cTargetOn = celleqel2mat_padded(mworks.cTargetOn);
    end
        cTargetOn(1) = nan; % first trial doesn't have reward 

load([analysis_out,date '_img' mouse '\getTC_' run '\' date '_img' mouse '_' run 'saveOutputs.mat']);   
load([analysis_out date '_img' mouse '\' date '_img' mouse '_' run '_deconvolution_thresh', num2str(thresholdDeco), '_TCave_cl.mat']);
load([analysis_out,date '_img' mouse '\getTC_' run '\' date '_img' mouse '_' run,...
    '_nPCA',num2str(nPCA), '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', ...
    num2str(cluster_threshold), '_TCave.mat'],'cTargetOn_cutted');
 
load([analysis_out img_fn 'figDataUntrimmed_' sessions '.mat'],'-regexp','^(?!expt|img|ind|seshType).'); 
animalTrials(1,mouseIdx) = size(targetAlign_events,3);
range{1,mouseIdx}=1:animalTrials(1,mouseIdx);
gTargetOn{1,mouseIdx} = cTargetOn;%(1,2:length(cTargetOn_cutted)+1); %concat without adding prev session value
gRctT = [gRctT celleqel2mat_padded(mworks.reactTimesMs)];
initFrames{1,mouseIdx} = mworks.counterValues{range{1,mouseIdx}(end)}(end);
initLickTimes{1,mouseIdx} = (mworks.lickometerTimesUs(1,range{1,mouseIdx})); 
initCounterTimes{1,mouseIdx} = (mworks.counterTimesUs(1,range{1,mouseIdx}));
initCounterVals{1,mouseIdx} = (mworks.counterValues(1,range{1,mouseIdx}));


%% calculate all the trial data for subplot figures identifying trials of a condition preceeded by n trials of the other/same condition
if size(celleqel2mat_padded(mworks.tBlock2TrialNumber),2) > size(targetAlign_events,3)
block2{1,mouseIdx} = celleqel2mat_padded(mworks.tBlock2TrialNumber(1,1:size(targetAlign_events,3)));
rew{1,mouseIdx} = double(~cell2mat(mworks.tBlock2TrialNumber(1,1:size(targetAlign_events,3))));
else
block2{1,mouseIdx} = celleqel2mat_padded(mworks.tBlock2TrialNumber(1,range{1,mouseIdx}));
rew{1,mouseIdx} = double(~cell2mat(mworks.tBlock2TrialNumber(1,range{1,mouseIdx})));    
end

if seshType == 2
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
groupAlign.goodcells{1,mouseIdx} = goodcells; 

if size(targetAlign_tc,2)==size(groupAlign.goodcells{1,mouseIdx},2) && size(targetAlign_tc,3)<size(range{1,mouseIdx},2)
groupAlign.tc{1,mouseIdx} = targetAlign_tc(:,:,:);
elseif size(targetAlign_tc,2)==size(groupAlign.goodcells{1,mouseIdx},2)
groupAlign.tc{1,mouseIdx} = targetAlign_tc(:,:,range{1,mouseIdx});
else
groupAlign.tc{1,mouseIdx} = targetAlign_tc(:,groupAlign.goodcells{1,mouseIdx},range{1,mouseIdx});
end
%adding segment to output dFoF as well
    if size(targetAlign_tc,2)==size(groupAlign.goodcells{1,mouseIdx},2)
    targetAlignF = mean(mean(targetAlign_tc(1:prewin_frames,:,:),3,'omitnan'),1,'omitnan'); % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.
    else
    targetAlignF = mean(mean(targetAlign_tc(1:prewin_frames,groupAlign.goodcells{1,mouseIdx},:),3,'omitnan'),1,'omitnan'); % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.
    end
groupAlign.dFoF{1,mouseIdx} = zeros(size(targetAlign_tc,1),size(targetAlignF,2),size(targetAlign_tc,3));%frame*cell*trials
    % calculate df/F
    for c = 1:size(targetAlignF,2)
        groupAlign.dFoF{1,mouseIdx}(:,c,:) = (targetAlign_tc(:,c,:)-targetAlignF(1,c))./targetAlignF(1,c);
    end
   
if size(targetAlign_events,2)==size(groupAlign.goodcells{1,mouseIdx},2) && size(targetAlign_events,3)<size(range{1,mouseIdx},2)
groupAlign.events{1,mouseIdx} = targetAlign_events(:,:,:);
elseif size(targetAlign_tc,2)==size(groupAlign.goodcells{1,mouseIdx},2)
groupAlign.events{1,mouseIdx} = targetAlign_events(:,:,range{1,mouseIdx});
else
groupAlign.events{1,mouseIdx} = targetAlign_events(:,groupAlign.goodcells{1,mouseIdx},range{1,mouseIdx});
end
gtt = [gtt; tt];
%ggoodcells = [ggoodcells goodcells];
clearvars -except sortlick* pmpt range thresholdDeco exp* rew* sesh* events loadData tc isesh exptCount piezo mouseIdx lick block2* mworks early* late* ind* iexp animalTrials init* tt target* g* date* day sessions mouse run i img_fn analysis_out bdata_source output_fn pairings

analysis_out = 'A:\home\carlo\RC\analysis\2P\';
load([analysis_out img_fn '_cueAlignLick.mat'],'-regexp','^(?!expt|img|ind|early|late|seshType).');

if size(lickCueAlign,2)<size(range{1,mouseIdx},2)
groupAlign.licks{1,mouseIdx} = lickCueAlign(:,:); 
groupAlign.lickStart{1,mouseIdx} = lickBurstStart(:,:);
else
groupAlign.licks{1,mouseIdx} = lickCueAlign(:,range{1,mouseIdx}); 
groupAlign.lickStart{1,mouseIdx} = lickBurstStart(:,range{1,mouseIdx});
end
  if sum(~isnan(groupAlign.lickStart{1,mouseIdx}))>6
        inImageWindow_lickStart{1,mouseIdx} = groupAlign.lickStart{1,mouseIdx};
        inImageWindow_lickStart{1,mouseIdx}(1,(groupAlign.lickStart{1,mouseIdx}>96)) = NaN;
        groupAlign.lickStart{1,mouseIdx} = inImageWindow_lickStart{1,mouseIdx};
        [sortlick, sortlick_ind] = sort(groupAlign.lickStart{1,mouseIdx},'ascend');
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

[sortlickB2{mouseIdx}, sortlick_indB2{mouseIdx}] = sort(groupAlign.lickStart{1,mouseIdx}(1,block2{1,mouseIdx} & ~isnan(groupAlign.lickStart{1,mouseIdx})),'ascend');
[sortlickRew{mouseIdx}, sortlick_indRew{mouseIdx}] = sort(groupAlign.lickStart{1,mouseIdx}(1,~block2{1,mouseIdx} & ~isnan(groupAlign.lickStart{1,mouseIdx})),'ascend');
early_bstB2{1,mouseIdx} = sortlick_indB2{mouseIdx}(1:floor(length(sortlickB2{mouseIdx})/4)); %trials with the earliest 25% of licks
late_bstB2{1,mouseIdx} = sortlick_indB2{mouseIdx}(length(sortlickB2{mouseIdx})-floor(length(sortlickB2{mouseIdx})/4)+1:length(sortlickB2{mouseIdx})); %trials with the latest 25% of licks
early_bstRew{1,mouseIdx} = sortlick_indRew{mouseIdx}(1:floor(length(sortlickRew{mouseIdx})/4)); %trials with the earliest 25% of licks
late_bstRew{1,mouseIdx} = sortlick_indRew{mouseIdx}(length(sortlickRew{mouseIdx})-floor(length(sortlickRew{mouseIdx})/4)+1:length(sortlickRew{mouseIdx})); %trials with the latest 25% of licks

%     early_bstOverLimInd{1,mouseIdx} = early_bst{1,mouseIdx}>96;
%     late_bstOverLimInd{1,mouseIdx} = late_bst{1,mouseIdx}>96;
%     early_bstOverLim{1,mouseIdx} = early_bst{1,mouseIdx};
%     early_bstOverLim{1,mouseIdx}(early_bstOverLimInd{1,mouseIdx}) = [];
%     late_bstOverLim{1,mouseIdx} = late_bst{1,mouseIdx};
%     late_bstOverLim{1,mouseIdx}(late_bstOverLimInd{1,mouseIdx}) = [];
    

earlyCSPLick_ind{1,mouseIdx} = early_bstRew{1,mouseIdx};
earlyCSMLick_ind{1,mouseIdx} = early_bstB2{1,mouseIdx};
lateCSPLick_ind{1,mouseIdx} = late_bstRew{1,mouseIdx};
lateCSMLick_ind{1,mouseIdx} = late_bstB2{1,mouseIdx};

tempB2Events = groupAlign.events{1,mouseIdx}(:,:,logical(block2{1,mouseIdx} & ~isnan(groupAlign.lickStart{1,mouseIdx})));
tempRewEvents = groupAlign.events{1,mouseIdx}(:,:,~block2{1,mouseIdx} & ~isnan(groupAlign.lickStart{1,mouseIdx}));
tempB2Start = groupAlign.lickStart{1,mouseIdx}(:,logical(block2{1,mouseIdx} & ~isnan(groupAlign.lickStart{1,mouseIdx})));
tempRewStart = groupAlign.lickStart{1,mouseIdx}(:,~block2{1,mouseIdx} & ~isnan(groupAlign.lickStart{1,mouseIdx}));

earlyBurstEvents = [earlyBurstEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,early_bst{1,mouseIdx}),3)];
lateBurstEvents = [lateBurstEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,late_bst{1,mouseIdx}),3)];
earlyBurstStart = [earlyBurstStart nanmean(groupAlign.lickStart{1,mouseIdx}(:,early_bst{1,mouseIdx}),3)];
lateBurstStart = [lateBurstStart nanmean(groupAlign.lickStart{1,mouseIdx}(:,late_bst{1,mouseIdx}),3)];

earlyRewBurstEvents = [earlyRewBurstEvents nanmean(tempRewEvents(:,:,earlyCSPLick_ind{1,mouseIdx}),3)];
lateRewBurstEvents = [lateRewBurstEvents nanmean(tempRewEvents(:,:,lateCSPLick_ind{1,mouseIdx}),3)];
earlyRewBurstStart = [earlyRewBurstStart nanmean(tempRewStart(:,earlyCSPLick_ind{1,mouseIdx}),3)];
lateRewBurstStart = [lateRewBurstStart nanmean(tempRewStart(:,lateCSPLick_ind{1,mouseIdx}),3)];
earlyB2BurstEvents = [earlyB2BurstEvents nanmean(tempB2Events(:,:,earlyCSMLick_ind{1,mouseIdx}),3)];
lateB2BurstEvents = [lateB2BurstEvents nanmean(tempB2Events(:,:,lateCSMLick_ind{1,mouseIdx}),3)];
earlyB2BurstStart = [earlyB2BurstStart nanmean(tempB2Start(:,earlyCSMLick_ind{1,mouseIdx}),3)];
lateB2BurstStart = [lateB2BurstStart nanmean(tempB2Start(:,lateCSMLick_ind{1,mouseIdx}),3)];
clearvars -except sortlick* pmpt range thresholdDeco exp* rew* sesh* events tc loadData isesh exptCount piezo lick mouseIdx block2* mworks iexp ind* early* late* animalTrials init* target* g* date* day sessions mouse run i img_fn analysis_out bdata_source output_fn pairings
 

  analysis_out = 'A:\home\carlo\RC\analysis\2P\';
load([analysis_out img_fn '_cueAlignPiezo.mat'],'-regexp','^(?!expt|img|ind|seshType).');
  
if size(targetAlign_piezo,2)<size(range{1,mouseIdx},2)
groupAlign.piezo{1,mouseIdx} = targetAlign_piezo(:,:);
else
groupAlign.piezo{1,mouseIdx} = targetAlign_piezo(:,range{1,mouseIdx});
end
clearvars -except sortlick* pmpt range thresholdDeco exp* rew* sesh* events tc isesh loadData exptCount piezo lick mouseIdx block2* iexp ind* early* late* animalTrials init* target* g* Goodcells date* day sessions mouse run i img_fn analysis_out bdata_source output_fn pairings


tc.rewGroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rew{1,mouseIdx})),3);
tc.b2GroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2{1,mouseIdx})),3);
events.rewGroupAlign{1,mouseIdx} = nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(rew{1,mouseIdx})),3);
events.b2GroupAlign{1,mouseIdx} = nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(block2{1,mouseIdx})),3);
piezo.rewGroupAlign{1,mouseIdx} = groupAlign.piezo{1,mouseIdx}(:,logical(rew{1,mouseIdx}));
piezo.b2GroupAlign{1,mouseIdx} = groupAlign.piezo{1,mouseIdx}(:,logical(block2{1,mouseIdx}));
lick.rewGroupAlign{1,mouseIdx} = groupAlign.licks{1,mouseIdx}(:,logical(rew{1,mouseIdx}));
lick.b2GroupAlign{1,mouseIdx} = groupAlign.licks{1,mouseIdx}(:,logical(block2{1,mouseIdx}));

tc.rewOneRewGroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rewOneRew{1,mouseIdx})),3);
events.rewOneRewGroupAlign{1,mouseIdx} = groupAlign.events{1,mouseIdx}(:,:,logical(rewOneRew{1,mouseIdx}));
piezo.rewOneRewGroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(rewOneRew{1,mouseIdx}));
lick.rewOneRewGroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(rewOneRew{1,mouseIdx}));
tc.b2OneRewGroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2OneRew{1,mouseIdx})),3);
events.b2OneRewGroupAlign{1,mouseIdx} = groupAlign.events{1,mouseIdx}(:,:,logical(block2OneRew{1,mouseIdx}));
piezo.b2OneRewGroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(block2OneRew{1,mouseIdx}));
lick.b2OneRewGroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(block2OneRew{1,mouseIdx}));

tc.rewTwoRewGroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rewTwoRew{1,mouseIdx})),3);
events.rewTwoRewGroupAlign{1,mouseIdx} = groupAlign.events{1,mouseIdx}(:,:,logical(rewTwoRew{1,mouseIdx}));
piezo.rewTwoRewGroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(rewTwoRew{1,mouseIdx}));
lick.rewTwoRewGroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(rewTwoRew{1,mouseIdx}));
tc.b2TwoRewGroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2TwoRew{1,mouseIdx})),3);
events.b2TwoRewGroupAlign{1,mouseIdx} = groupAlign.events{1,mouseIdx}(:,:,logical(block2TwoRew{1,mouseIdx}));
piezo.b2TwoRewGroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(block2TwoRew{1,mouseIdx}));
lick.b2TwoRewGroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(block2TwoRew{1,mouseIdx}));

tc.rewThreeRewGroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rewThreeRew{1,mouseIdx})),3);
events.rewThreeRewGroupAlign{1,mouseIdx} = groupAlign.events{1,mouseIdx}(:,:,logical(rewThreeRew{1,mouseIdx}));
piezo.rewThreeRewGroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(rewThreeRew{1,mouseIdx}));
lick.rewThreeRewGroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(rewThreeRew{1,mouseIdx}));
tc.b2ThreeRewGroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2ThreeRew{1,mouseIdx})),3);
events.b2ThreeRewGroupAlign{1,mouseIdx} = groupAlign.events{1,mouseIdx}(:,:,logical(block2ThreeRew{1,mouseIdx}));
piezo.b2ThreeRewGroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(block2ThreeRew{1,mouseIdx}));
lick.b2ThreeRewGroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(block2ThreeRew{1,mouseIdx}));

tc.rewOneB2GroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rewOneB2{1,mouseIdx})),3);
events.rewOneB2GroupAlign{1,mouseIdx} = groupAlign.events{1,mouseIdx}(:,:,logical(rewOneB2{1,mouseIdx}));
piezo.rewOneB2GroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(rewOneB2{1,mouseIdx}));
lick.rewOneB2GroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(rewOneB2{1,mouseIdx}));
tc.b2OneB2GroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2OneB2{1,mouseIdx})),3);
events.b2OneB2GroupAlign{1,mouseIdx} = groupAlign.events{1,mouseIdx}(:,:,logical(block2OneB2{1,mouseIdx}));
piezo.b2OneB2GroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(block2OneB2{1,mouseIdx}));
lick.b2OneB2GroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(block2OneB2{1,mouseIdx}));

tc.rewTwoB2GroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rewTwoB2{1,mouseIdx})),3);
events.rewTwoB2GroupAlign{1,mouseIdx} = groupAlign.events{1,mouseIdx}(:,:,logical(rewTwoB2{1,mouseIdx}));
piezo.rewTwoB2GroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(rewTwoB2{1,mouseIdx}));
lick.rewTwoB2GroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(rewTwoB2{1,mouseIdx}));
tc.b2TwoB2GroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2TwoB2{1,mouseIdx})),3);
events.b2TwoB2GroupAlign{1,mouseIdx} = groupAlign.events{1,mouseIdx}(:,:,logical(block2TwoB2{1,mouseIdx}));
piezo.b2TwoB2GroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(block2TwoB2{1,mouseIdx}));
lick.b2TwoB2GroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(block2TwoB2{1,mouseIdx}));

tc.rewThreeB2GroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rewThreeB2{1,mouseIdx})),3);
events.rewThreeB2GroupAlign{1,mouseIdx} = groupAlign.events{1,mouseIdx}(:,:,logical(rewThreeB2{1,mouseIdx}));
piezo.rewThreeB2GroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(rewThreeB2{1,mouseIdx}));
lick.rewThreeB2GroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(rewThreeB2{1,mouseIdx}));
tc.b2ThreeB2GroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2ThreeB2{1,mouseIdx})),3);
events.b2ThreeB2GroupAlign{1,mouseIdx} = groupAlign.events{1,mouseIdx}(:,:,logical(block2ThreeB2{1,mouseIdx}));
piezo.b2ThreeB2GroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(block2ThreeB2{1,mouseIdx}));
lick.b2ThreeB2GroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(block2ThreeB2{1,mouseIdx}));

tc.rewFourB2GroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rewFourB2{1,mouseIdx})),3);
events.rewFourB2GroupAlign{1,mouseIdx} = groupAlign.events{1,mouseIdx}(:,:,logical(rewFourB2{1,mouseIdx}));
piezo.rewFourB2GroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(rewFourB2{1,mouseIdx}));
lick.rewFourB2GroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(rewFourB2{1,mouseIdx}));
tc.b2FourB2GroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2FourB2{1,mouseIdx})),3);
events.b2FourB2GroupAlign{1,mouseIdx} = groupAlign.events{1,mouseIdx}(:,:,logical(block2FourB2{1,mouseIdx}));
piezo.b2FourB2GroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(block2FourB2{1,mouseIdx}));
lick.b2FourB2GroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(block2FourB2{1,mouseIdx}));

tc.rewFourRewGroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rewFourRew{1,mouseIdx})),3);
events.rewFourRewGroupAlign{1,mouseIdx} = groupAlign.events{1,mouseIdx}(:,:,logical(rewFourRew{1,mouseIdx}));
piezo.rewFourRewGroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(rewFourRew{1,mouseIdx}));
lick.rewFourRewGroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(rewFourRew{1,mouseIdx}));
tc.b2FourRewGroupAlign{1,mouseIdx} = nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2FourRew{1,mouseIdx})),3);
events.b2FourRewGroupAlign{1,mouseIdx} = groupAlign.events{1,mouseIdx}(:,:,logical(block2FourRew{1,mouseIdx}));
piezo.b2FourRewGroupAlign{1,mouseIdx} =  groupAlign.piezo{1,mouseIdx}(:,logical(block2FourRew{1,mouseIdx}));
lick.b2FourRewGroupAlign{1,mouseIdx} =  groupAlign.licks{1,mouseIdx}(:,logical(block2FourRew{1,mouseIdx}));

end
end
%
if seshType == 2 && loadData ~= 2
    blockName = input('which block is this? [type "plus/minus""First/Second"]: ','s');
end
%
analysis_out = 'A:\home\carlo\RC\analysis\2P\';
output_fn = ([analysis_out 'groupwiseAlign' num2str(thresholdDeco) '_' sesh '_' sesh2 '\']);

groupAlign_tc=[]; groupAlign_events=[]; groupAlign_piezo=[]; groupAlign_lick=[];
for a = 1:exptCount
nsize(1,a) = size(groupAlign.events{1,a},2);
end


for mouseIdx = 1:exptCount

groupAlign_tc = cat(2,groupAlign_tc, nanmean(groupAlign.tc{1,mouseIdx},3));%(:,rDend(i,:)
groupAlign_events = cat(2,groupAlign_events, nanmean(groupAlign.events{1,mouseIdx},3));
groupAlign_piezo = [groupAlign_piezo groupAlign.piezo{1,mouseIdx}];
groupAlign_lick = [groupAlign_lick groupAlign.licks{1,mouseIdx}];
 
end

% for i = 1:sum(~isnan(groupAlign_piezo(1,:)))
% if isnan(groupAlign_piezo(1,i)); groupAlign_piezo(:,i)=[]; end
% if isnan(groupAlign_lick(1,i)); groupAlign_lick(:,i)=[]; end
% end
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
%% piezo align formatting (to be shoved into additional formatting for lick align to pull out first motion 2o above u after reward

for imouse = 1:size(groupAlign.piezo,2)
    piezoWholeSession{1,imouse} = []; piezoWholeSessionPostReward{1,imouse} = []; piezoWholeSessionPreCue{1,imouse} = [];
for itrial = 1:size(groupAlign.piezo{1,imouse},2)
piezoWholeSession{1,imouse} = [piezoWholeSession{1,imouse}; groupAlign.piezo{1,imouse}(:,itrial)];
piezoWholeSessionPreCue{1,imouse} = [piezoWholeSessionPreCue{1,imouse}; groupAlign.piezo{1,imouse}(1:50,itrial)];
piezoWholeSessionPostReward{1,imouse} = [piezoWholeSessionPostReward{1,imouse}; groupAlign.piezo{1,imouse}(73:150,itrial)];
end

uWholeSession(1,imouse) = nanmean(piezoWholeSession{1,imouse});
oWholeSession(1,imouse) = nanstd(piezoWholeSession{1,imouse},[],1);
uWholeSessionPreCue(1,imouse) = nanmean(piezoWholeSessionPreCue{1,imouse});
oWholeSessionPreCue(1,imouse) = nanstd(piezoWholeSessionPreCue{1,imouse},[],1);
uWholeSessionPostReward(1,imouse) = nanmean(piezoWholeSessionPostReward{1,imouse});
oWholeSessionPostReward(1,imouse) = nanstd(piezoWholeSessionPostReward{1,imouse},[],1);
end

tempFig = setFigure; hold on;
plot(1:length(piezoWholeSession{1,imouse}),piezoWholeSession{1,imouse},'k');
hold on;
yline(nanmean(uWholeSession),'r--');
yline((nanmean(oWholeSession).*1)+(nanmean(uWholeSession)),'Color',[0.8500, 0.3250, 0.0980]);
yline((nanmean(oWholeSession).*2)+(nanmean(uWholeSession)),'Color',[0.9290, 0.6940, 0.1250]);
yline((nanmean(oWholeSession).*3)+(nanmean(uWholeSession)),'g');
yline((nanmean(oWholeSession).*4)+(nanmean(uWholeSession)),'b');
yline((nanmean(oWholeSession).*5)+(nanmean(uWholeSession)),'c');
yline((nanmean(oWholeSession).*6)+(nanmean(uWholeSession)),'m');
legend({'piezo','mean','1SD','2SD','3SD','4SD','5SD','6SD'},'Location','northeast');
ylim([0 1])
title('piezo transients across imaging session');
supertitle('mean and s.d. found using whole dataset');


tempFig = setFigure; hold on;
plot(1:length(piezoWholeSession{1,imouse}),piezoWholeSession{1,imouse},'k');
hold on;
yline(nanmean(uWholeSessionPreCue),'r--');
yline((nanmean(oWholeSessionPreCue).*1)+(nanmean(uWholeSessionPreCue)),'Color',[0.8500, 0.3250, 0.0980]);
yline((nanmean(oWholeSessionPreCue).*2)+(nanmean(uWholeSessionPreCue)),'Color',[0.9290, 0.6940, 0.1250]);
yline((nanmean(oWholeSessionPreCue).*3)+(nanmean(uWholeSessionPreCue)),'g');
yline((nanmean(oWholeSessionPreCue).*4)+(nanmean(uWholeSessionPreCue)),'b');
yline((nanmean(oWholeSessionPreCue).*5)+(nanmean(uWholeSessionPreCue)),'c');
yline((nanmean(oWholeSessionPreCue).*6)+(nanmean(uWholeSessionPreCue)),'m');
legend({'piezo','mean','1SD','2SD','3SD','4SD','5SD','6SD'},'Location','northeast');
ylim([0 1])
title('piezo transients across of imaging session');
supertitle('mean and s.d. found using only pre-cue frames');


tempFig = setFigure; hold on;
plot(1:length(piezoWholeSessionPostReward{1,imouse}),piezoWholeSessionPostReward{1,imouse},'k');
hold on;
yline(nanmean(uWholeSessionPreCue),'r--');
yline((nanmean(oWholeSessionPreCue).*1)+(nanmean(uWholeSessionPreCue)),'Color',[0.8500, 0.3250, 0.0980]);
yline((nanmean(oWholeSessionPreCue).*2)+(nanmean(uWholeSessionPreCue)),'Color',[0.9290, 0.6940, 0.1250]);
yline((nanmean(oWholeSessionPreCue).*3)+(nanmean(uWholeSessionPreCue)),'g');
yline((nanmean(oWholeSessionPreCue).*4)+(nanmean(uWholeSessionPreCue)),'b');
yline((nanmean(oWholeSessionPreCue).*5)+(nanmean(uWholeSessionPreCue)),'c');
yline((nanmean(oWholeSessionPreCue).*6)+(nanmean(uWholeSessionPreCue)),'m');
legend({'piezo','mean','1SD','2SD','3SD','4SD','5SD','6SD'},'Location','northeast');
ylim([0 1])
title('piezo transients across post-reward frames of imaging session');
supertitle('mean and s.d. found using pre-cue frames');

%create a logical (150XnTrial) indicating frames where value resides 2sd or
%1sd above the average for that session

for imouse = 1:size(groupAlign.piezo,2)
    motionPeaks{1,imouse} = zeros(150,size(groupAlign.piezo{1,imouse},2));
    motionPeaks3{1,imouse} = zeros(150,size(groupAlign.piezo{1,imouse},2));
    motionPeaks2{1,imouse} = zeros(150,size(groupAlign.piezo{1,imouse},2));
    motionPeaks1{1,imouse} = zeros(150,size(groupAlign.piezo{1,imouse},2));
    cutoffPiezo(1,imouse) = uWholeSessionPreCue(1,imouse)+(3.*oWholeSessionPreCue(1,imouse));
for itrial = 1:size(groupAlign.piezo{1,imouse},2)
    for iframe = 1:149
        changeThisTrial(1,iframe) = groupAlign.piezo{1,imouse}(iframe+1,itrial)-groupAlign.piezo{1,imouse}(iframe,itrial);
    end
    cutoffChange3(imouse,itrial) = mean(changeThisTrial,2,'omitnan')+(3.*std(changeThisTrial,[],2,'omitnan'));
    cutoffChange2(imouse,itrial) = mean(changeThisTrial,2,'omitnan')+(2.*std(changeThisTrial,[],2,'omitnan'));
    cutoffChange1(imouse,itrial) = mean(changeThisTrial,2,'omitnan')+(1.*std(changeThisTrial,[],2,'omitnan'));
    for iframe = 2:150
        if groupAlign.piezo{1,imouse}(iframe,itrial) > cutoffPiezo(1,imouse) && changeThisTrial(1,iframe-1)>cutoffChange3(imouse,itrial)
        motionPeaks3{1,imouse}(iframe,itrial) = 1;
        elseif groupAlign.piezo{1,imouse}(iframe,itrial) > cutoffPiezo(1,imouse) && changeThisTrial(1,iframe-1)>cutoffChange2(imouse,itrial)
        motionPeaks2{1,imouse}(iframe,itrial) = 1;            
        elseif groupAlign.piezo{1,imouse}(iframe,itrial) > cutoffPiezo(1,imouse) && changeThisTrial(1,iframe-1)>cutoffChange1(imouse,itrial)
        motionPeaks1{1,imouse}(iframe,itrial) = 1;    
        else
        end
        if groupAlign.piezo{1,imouse}(iframe,itrial) > cutoffPiezo(1,imouse)
        motionPeaks{1,imouse}(iframe,itrial) = 1;
        else
        end
    end
end
end

gMotionPeaks = cell2mat(motionPeaks); gBlock2 = cell2mat(block2); gBlock2 = logical(gBlock2);
gMotionPeaks3 = cell2mat(motionPeaks3); gMotionPeaks2 = cell2mat(motionPeaks2); gMotionPeaks1 = cell2mat(motionPeaks1);
%% lick align formatting

  nTrials = sum(animalTrials);
    rewDelay_frames = round((cRctT/1000).*frameRateHz);% there's 700ms between the cue and the reward delivery !!!!! if you change tooFastMs in the varibles in MWorks, this needs to be changed
    nIC = size(groupAlign_events,2); %you get targetalign_events from neural align, this is the neural data (in spikes) aligned to cue
    lickCueAlign =  nan(prewin_frames+postwin_frames,nTrials);
    lickCounterVals = cell(1,nTrials);
    for a = 1:exptCount
    nFrames{1,a} = initCounterVals{1,a}{animalTrials(1,a)}(end);
    end
    lickDelay_frames =  round(0.1.*frameRateHz);
    lickSearch_frames =  round(0.3.*frameRateHz);
    
    postLick_frames = round(0.5.*frameRateHz);
    lickSearchRange = prewin_frames+lickDelay_frames:size(lickCueAlign,1)-lickSearch_frames;
    preRew_lickSearchRange = prewin_frames+1:prewin_frames+rewDelay_frames+1;
    lickBurstStart_postReward = nan(1,nTrials);
    lickBurstStart_preReward = nan(1,nTrials);
    postRew_lickBurstStart = nan(1,nTrials);
    preRew_lickBurstStart = nan(1,nTrials);
    lickBurstHz_all = nan(1,nTrials);
    preRew_lickBurstHz = nan(1,nTrials);
    preRewCSP_lickBurstHz = nan(1,nTrials);
    preRewCSM_lickBurstHz = nan(1,nTrials);
    postRew_lickBurstHz = nan(1,nTrials);
    postRewCSP_lickBurstHz = nan(1,nTrials);
    postRewCSM_lickBurstHz = nan(1,nTrials);
    postRew_lickAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    postRew_lickAlign = nan(3.*postLick_frames,nTrials);
    preRew_lickAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    preRew_lickAlign = nan(3.*postLick_frames,nTrials);
    lastPreRew_lickAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    lastPreRew_lickAlignDFoF = nan(3.*postLick_frames,nIC,nTrials);
    firstPostRew_lickAlignEvents = nan(3.*postLick_frames,nIC,nTrials); 
    firstPostRew_lickAlignDFoF = nan(3.*postLick_frames,nIC,nTrials); 
    firstPreRew_lickAlignEvents = nan(3.*postLick_frames,nIC,nTrials); 
    firstPreRew_lickAlignDFoF = nan(3.*postLick_frames,nIC,nTrials); 
    lastPreRew_lickAlign = nan(3.*postLick_frames,nTrials); 
    firstPostRew_lickAlign = nan(3.*postLick_frames,nTrials);
    firstPreRew_lickAlign = nan(3.*postLick_frames,nTrials);
    
    firstLickAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    firstLickAlignDFoF = nan(3.*postLick_frames,nIC,nTrials);
    firstLickAlign = nan(3.*postLick_frames,nTrials);
    firstPiezoAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    firstPiezoAlignDFoF = nan(3.*postLick_frames,nIC,nTrials);
    firstPiezoAlign = nan(3.*postLick_frames,nTrials);
    % piezo align
    firstPostRew_piezoAlignEvents = nan(3.*postLick_frames,nIC,nTrials); 
    firstPreRew_piezoAlignEvents = nan(3.*postLick_frames,nIC,nTrials); 
    firstPostRew_piezoAlignDFoF = nan(3.*postLick_frames,nIC,nTrials); 
    firstPreRew_piezoAlignDFoF = nan(3.*postLick_frames,nIC,nTrials); 
    firstPostRew_piezoAlign = nan(3.*postLick_frames,nTrials);
    firstPreRew_piezoAlign = nan(3.*postLick_frames,nTrials);
    
    lastPreRew_piezoAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    lastPreRew_piezoAlignDFoF = nan(3.*postLick_frames,nIC,nTrials);
    lastPreRew_piezoAlign = nan(3.*postLick_frames,nTrials);    
    
    preRewPiezoTrials = zeros(1,nTrials);
    preFRewPiezoTrials = zeros(1,nTrials);
    postRewPiezoTrials = zeros(1,nTrials);
    
    lastPreRewPiezoFrame = nan(1,nTrials);
    firstPostRewPiezoFrame = nan(1,nTrials);
    %
    rewAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    rewAlignDFoF = nan(3.*postLick_frames,nIC,nTrials);
    preRewTrials = zeros(1,nTrials);
    preFRewTrials = zeros(1,nTrials);
    postRewTrials = zeros(1,nTrials);
    postRew_lickSearchRange = prewin_frames+rewDelay_frames+1:prewin_frames+rewDelay_frames+rewDelay_frames;% need to '-lickSearch_frames-postLick_frames' because later the inds needs to + lickSearch_frames or +postLick_frames, this is just for the following inds to be within matrix dimentions
    preRew_lickSearchRange = prewin_frames+1:prewin_frames+rewDelay_frames;% need to '-lickSearch_frames-postLick_frames' because later the inds needs to + lickSearch_frames or +postLick_frames, this is just for the following inds to be within matrix dimentions
    lastPreRewLickFrame = nan(1,nTrials);
    firstPostRewLickFrame = nan(1,nTrials);
    %
    firstPostRew_lickStart = NaN(1,nTrials);
    lastPreRew_lickStart = NaN(1,nTrials);
    firstPreRew_lickStart = NaN(1,nTrials);
    firstPostRew_piezoStart = NaN(1,nTrials);
    lastPreRew_piezoStart = NaN(1,nTrials);
    firstPreRew_piezoStart = NaN(1,nTrials);
    firstPostRew_lickTrial = NaN(1,nTrials);
    lastPreRew_lickTrial = NaN(1,nTrials);
    firstPreRew_lickTrial = NaN(1,nTrials);
    firstPostRew_piezoTrial = NaN(1,nTrials);
    lastPreRew_piezoTrial = NaN(1,nTrials);
    firstPreRew_piezoTrial = NaN(1,nTrials);
    %
    isolate_lastLickInd = NaN(1,nTrials);
    lastLick_lickAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    lastLick_lickAlignDFoF = nan(3.*postLick_frames,nIC,nTrials);
    lastLick_lickAlign = NaN(3.*postLick_frames,nTrials);
    isolate_lastMoveInd = NaN(1,nTrials);
    lastMove_piezoAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    lastMove_piezoAlignDFoF = nan(3.*postLick_frames,nIC,nTrials);
    lastMove_piezoAlign = NaN(3.*postLick_frames,nTrials);
    %
    firstLickStart = nan(1,nTrials); firstLickStartAll = cell(1,exptCount);
    firstLickHz = nan(1,nTrials);
    firstLickAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    firstLickAlignDFoF = nan(3.*postLick_frames,nIC,nTrials);
    firstLickAlign = nan(3.*postLick_frames,nTrials);
    firstLickTrial = nan(1,nTrials);    
    firstPiezoStart = nan(1,nTrials);
    firstPiezoHz = nan(1,nTrials);
    firstPiezoAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    firstPiezoAlignDFoF = nan(3.*postLick_frames,nIC,nTrials);
    firstPiezoAlign = nan(3.*postLick_frames,nTrials);
    firstPiezoTrial = nan(1,nTrials);
    %
    firstLickBurstStart = nan(1,nTrials);
    firstLickBurstHz = nan(1,nTrials);
    firstLickBurstAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    firstLickBurstAlignDFoF = nan(3.*postLick_frames,nIC,nTrials);
    firstLickBurstAlign = nan(3.*postLick_frames,nTrials);
    firstLickBurstTrial = nan(1,nTrials);
    firstSolitaryLickStart = nan(1,nTrials);
    firstSolitaryLickHz = nan(1,nTrials);
    firstSolitaryLickAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
    firstSolitaryLickAlignDFoF = nan(3.*postLick_frames,nIC,nTrials);
    firstSolitaryLickAlign = nan(3.*postLick_frames,nTrials);
    firstSolitaryLickTrial = nan(1,nTrials);
    for i = 1:4
    firstPiezoBurstStart{i,1} = nan(1,nTrials);
    firstPiezoBurstHz{i,1} = nan(1,nTrials);
    firstPiezoBurstAlignEvents{i,1} = nan(3.*postLick_frames,nIC,nTrials);
    firstPiezoBurstAlignDFoF{i,1} = nan(3.*postLick_frames,nIC,nTrials);
    firstPiezoBurstAlign{i,1} = nan(3.*postLick_frames,nTrials);
    firstPiezoBurstTrial{i,1} = nan(1,nTrials);
    firstSolitaryPiezoStart{i,1} = nan(1,nTrials);
    firstSolitaryPiezoHz{i,1} = nan(1,nTrials);
    firstSolitaryPiezoAlignEvents{i,1} = nan(3.*postLick_frames,nIC,nTrials);
    firstSolitaryPiezoAlignDFoF{i,1} = nan(3.*postLick_frames,nIC,nTrials);
    firstSolitaryPiezoAlign{i,1} = nan(3.*postLick_frames,nTrials);
    firstSolitaryPiezoTrial{i,1} = nan(1,nTrials);
    end
    %
    tTInx = 1; ntrials=0;
    for a=1:exptCount %lick align organization 
oDend = (1:nsize(1,a))+sum(nsize(1,1:a-1));
nDend = 1:nsize(1,a);       

    for itrial = 1:animalTrials(1,a) %nTrials
        ind_post = [];
        ind_pre = [];
        %for piezo
        ind_postP = [];
        ind_preP = [];
        %
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
                   ansTemp{1,a}(1,itrial)=initCounterVals{1,a}{itrial}(end)-gTargetOn{1,a}(itrial);
                if initCounterVals{1,a}{itrial}(end)-gTargetOn{1,a}(itrial) > postwin_frames
                    %lickcuealign aligns licking of each frame to cue
                    lickCueAlign(:,ntrial) = lickTC{ntrial}(1,gTargetOn{1,a}(itrial)-prewin_frames-counterVals(1):gTargetOn{1,a}(itrial)+postwin_frames-1-counterVals(1));
                end
                ind = intersect(lickSearchRange,find(lickCueAlign(:,ntrial)));
             %%%this for loop find frames where more than 3 licks occurred within a fixed window, and index that frame as the start of burst lick
                for i = 1:length(ind)
                    ilick = ind(i);
                    if sum(lickCueAlign(ilick:ilick+lickSearch_frames-1,ntrial),1) >= 3 % more than 3 bursts within about 300ms
                        lickBurstStart_postReward(:,ntrial) = ilick;% when licking burst happens
                        break
                    end
                end
             %%% simply, find the first lick with frame index
                ind = find(lickCueAlign(:,ntrial));
                for i = 1:length(ind)
                    ilick = ind(i);
                    if ~isnan(firstLickBurstStart(:,ntrial)) && ~isnan(firstSolitaryLickStart(:,ntrial))
                    elseif ilick+lickSearch_frames-1 > 150 || ilick+postLick_frames+postLick_frames-1 > 150
                    elseif ilick<=postLick_frames || ilick<=lickSearch_frames
                    else
                        if isnan(firstLickStart(:,ntrial)) && sum(lickCueAlign(ilick-postLick_frames:ilick-1,ntrial),1) == 0 && ilick+postLick_frames+postLick_frames <= 150
                        firstLickStart(:,ntrial) = ilick;
                        firstLickHz(:,ntrial) = length(ind)./(cRctT/1000);
                        firstLickAlignEvents(:,oDend,ntrial) = groupAlign.events{1,a}(ilick-postLick_frames+1:ilick+postLick_frames+postLick_frames,nDend,itrial);% align neural data to first lick
                        firstLickAlignDFoF(:,oDend,ntrial) = groupAlign.dFoF{1,a}(ilick-postLick_frames+1:ilick+postLick_frames+postLick_frames,nDend,itrial);% align neural data to first lick
                        firstLickAlign(:,ntrial) = lickCueAlign(ilick-postLick_frames+1:ilick+postLick_frames+postLick_frames,ntrial);
                        firstLickTrial(:,ntrial) = itrial;
                        end
                        
                        if isnan(firstPreRew_lickStart(:,ntrial)) && ilick>prewin_frames && ilick<prewin_frames+rewDelay_frames && sum(lickCueAlign(ilick-postLick_frames:ilick-1,ntrial),1) == 0
                        firstPreRew_lickStart(:,ntrial) = ilick;
                        firstPreRew_lickHz(:,ntrial) = length(ind)./(cRctT/1000);
                        firstPreRew_lickAlignEvents(:,oDend,ntrial) = groupAlign.events{1,a}(ilick-postLick_frames:ilick+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick
                        firstPreRew_lickAlignDFoF(:,oDend,ntrial) = groupAlign.dFoF{1,a}(ilick-postLick_frames:ilick+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick
                        firstPreRew_lickAlign(:,ntrial) = lickCueAlign(ilick-postLick_frames:ilick+postLick_frames+postLick_frames-1,ntrial);
                        firstPreRew_lickTrial(:,ntrial) = itrial;
                        end
                        
                        if isnan(firstPostRew_lickStart(:,ntrial)) && ilick>prewin_frames+rewDelay_frames && ilick<prewin_frames+(2*rewDelay_frames) && sum(lickCueAlign(ilick-postLick_frames:ilick-1,ntrial),1) == 0
                        firstPostRew_lickStart(:,ntrial) = ilick;
                        firstPostRew_lickHz(:,ntrial) = length(ind)./(cRctT/1000);
                        firstPostRew_lickAlignEvents(:,oDend,ntrial) = groupAlign.events{1,a}(ilick-postLick_frames:ilick+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick
                        firstPostRew_lickAlignDFoF(:,oDend,ntrial) = groupAlign.dFoF{1,a}(ilick-postLick_frames:ilick+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick
                        firstPostRew_lickAlign(:,ntrial) = lickCueAlign(ilick-postLick_frames:ilick+postLick_frames+postLick_frames-1,ntrial);
                        firstPostRew_lickTrial(:,ntrial) = itrial;
                        end
                        
                        if sum(lickCueAlign(ilick:ilick+lickSearch_frames-1,ntrial),1) >= 3 && isnan(firstLickBurstStart(:,ntrial)) && sum(lickCueAlign(ilick-postLick_frames:ilick-1,ntrial),1) == 0
                        firstLickBurstStart(:,ntrial) = ilick;
                        firstLickBurstHz(:,ntrial) = length(ind)./(cRctT/1000);
                        firstLickBurstAlignEvents(:,oDend,ntrial) = groupAlign.events{1,a}(ilick-postLick_frames:ilick+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick
                        firstLickBurstAlignDFoF(:,oDend,ntrial) = groupAlign.dFoF{1,a}(ilick-postLick_frames:ilick+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick
                        firstLickBurstAlign(:,ntrial) = lickCueAlign(ilick-postLick_frames:ilick+postLick_frames+postLick_frames-1,ntrial);% align licking data to first lick
                        firstLickBurstTrial(:,ntrial) = itrial;
                        elseif sum(lickCueAlign(ilick:ilick+lickSearch_frames-1,ntrial),1) == 1 && isnan(firstSolitaryLickStart(:,ntrial)) && sum(lickCueAlign(ilick-postLick_frames:ilick-1,ntrial),1) == 0
                        firstSolitaryLickStart(:,ntrial) = ilick;
                        firstSolitaryLickHz(:,ntrial) = length(ind)./(cRctT/1000);  
                        firstSolitaryLickAlignEvents(:,oDend,ntrial) = groupAlign.events{1,a}(ilick-postLick_frames:ilick+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick
                        firstSolitaryLickAlignDFoF(:,oDend,ntrial) = groupAlign.dFoF{1,a}(ilick-postLick_frames:ilick+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick
                        firstSolitaryLickAlign(:,ntrial) = lickCueAlign(ilick-postLick_frames:ilick+postLick_frames+postLick_frames-1,ntrial);% align licking data to first lick
                        firstSolitaryLickTrial(:,ntrial) = itrial;
                        end
                    end
                end
                %pull out the motion events array for different cutoff values
                for cutoffIdx = 1:4
                    if cutoffIdx < 4
                    eval(['tempMotionPeaks = gMotionPeaks' num2str(cutoffIdx) ';']);
                    elseif cutoffIdx == 4
                    eval('tempMotionPeaks = gMotionPeaks;');
                    end
                indP = find(tempMotionPeaks(:,ntrial));
                for i = 1:length(indP)
                    ipiezo = indP(i);
                    if ~isnan(firstPiezoBurstStart{cutoffIdx,1}(:,ntrial)) && ~isnan(firstSolitaryPiezoStart{cutoffIdx,1}(:,ntrial))
                    elseif ipiezo+lickSearch_frames-1 > 150 || ipiezo+postLick_frames+postLick_frames-1 > 150
                    elseif ipiezo<=lickSearch_frames || ipiezo<=postLick_frames
                    else
                        if isnan(firstPiezoStart(:,ntrial)) && sum(tempMotionPeaks(ipiezo-postLick_frames:ipiezo-1,ntrial),1) == 0 && cutoffIdx == 4
                        firstPiezoStart(:,ntrial) = ipiezo;
                        firstPiezoHz(:,ntrial) = length(indP)./(cRctT/1000);
                        firstPiezoAlignEvents(:,oDend,ntrial) = groupAlign.events{1,a}(ipiezo-postLick_frames:ipiezo+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick
                        firstPiezoAlignDFoF(:,oDend,ntrial) = groupAlign.dFoF{1,a}(ipiezo-postLick_frames:ipiezo+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick
                        firstPiezoAlign(:,ntrial) = tempMotionPeaks(ipiezo-postLick_frames:ipiezo+postLick_frames+postLick_frames-1,ntrial);
                        firstPiezoTrial(:,ntrial) = itrial;
                        end
                        
                        if isnan(firstPreRew_piezoStart(:,ntrial)) && ipiezo>prewin_frames && ipiezo<prewin_frames+rewDelay_frames && sum(tempMotionPeaks(ipiezo-postLick_frames:ipiezo-1,ntrial),1) == 0
                        firstPreRew_piezoStart(:,ntrial) = ipiezo;
                        firstPreRew_piezoHz(:,ntrial) = length(ind)./(cRctT/1000);
                        firstPreRew_piezoAlignEvents(:,oDend,ntrial) = groupAlign.events{1,a}(ipiezo-postLick_frames:ipiezo+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick
                        firstPreRew_piezoAlignDFoF(:,oDend,ntrial) = groupAlign.dFoF{1,a}(ipiezo-postLick_frames:ipiezo+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick
                        firstPreRew_piezoAlign(:,ntrial) = tempMotionPeaks(ipiezo-postLick_frames:ipiezo+postLick_frames+postLick_frames-1,ntrial);
                        firstPreRew_piezoTrial(:,ntrial) = itrial;
                        end
                        
                        if isnan(firstPostRew_piezoStart(:,ntrial)) && ipiezo>prewin_frames+rewDelay_frames && ipiezo<prewin_frames+(2*rewDelay_frames) && sum(tempMotionPeaks(ipiezo-postLick_frames:ipiezo-1,ntrial),1) == 0
                        firstPostRew_piezoStart(:,ntrial) = ipiezo;
                        firstPostRew_piezoHz(:,ntrial) = length(ind)./(cRctT/1000);
                        firstPostRew_piezoAlignEvents(:,oDend,ntrial) = groupAlign.events{1,a}(ipiezo-postLick_frames:ipiezo+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick
                        firstPostRew_piezoAlignDFoF(:,oDend,ntrial) = groupAlign.dFoF{1,a}(ipiezo-postLick_frames:ipiezo+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick
                        firstPostRew_piezoAlign(:,ntrial) = tempMotionPeaks(ipiezo-postLick_frames:ipiezo+postLick_frames+postLick_frames-1,ntrial);
                        firstPostRew_piezoTrial(:,ntrial) = itrial;
                        end
                        
                        if sum(tempMotionPeaks(ipiezo:ipiezo+lickSearch_frames-1,ntrial),1) >= 3 && isnan(firstPiezoBurstStart{cutoffIdx,1}(:,ntrial)) && sum(tempMotionPeaks(ipiezo-postLick_frames:ipiezo-1,ntrial),1) == 0
                        firstPiezoBurstStart{cutoffIdx,1}(:,ntrial) = ipiezo;
                        firstPiezoBurstHz{cutoffIdx,1}(:,ntrial) = length(indP)./(cRctT/1000);
                        
                        firstPiezoBurstAlignEvents{cutoffIdx,1}(:,oDend,ntrial) = groupAlign.events{1,a}(ipiezo-postLick_frames:ipiezo+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick
                        firstPiezoBurstAligndFoF{cutoffIdx,1}(:,oDend,ntrial) = groupAlign.dFoF{1,a}(ipiezo-postLick_frames:ipiezo+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick
                        firstPiezoBurstAlign{cutoffIdx,1}(:,ntrial) = tempMotionPeaks(ipiezo-postLick_frames:ipiezo+postLick_frames+postLick_frames-1,ntrial);% align licking data to first lick
                        firstPiezoBurstTrial{cutoffIdx,1}(:,ntrial) = itrial;
                        elseif sum(tempMotionPeaks(ipiezo:ipiezo+lickSearch_frames-1,ntrial),1) == 1 && isnan(firstSolitaryPiezoStart{cutoffIdx,1}(:,ntrial)) && sum(tempMotionPeaks(ipiezo-postLick_frames:ipiezo-1,ntrial),1) == 0
                        firstSolitaryPiezoStart{cutoffIdx,1}(:,ntrial) = ipiezo;
                        firstSolitaryPiezoHz{cutoffIdx,1}(:,ntrial) = length(indP)./(cRctT/1000);
                        
                        firstSolitaryPiezoAlignEvents{cutoffIdx,1}(:,oDend,ntrial) = groupAlign.events{1,a}(ipiezo-postLick_frames:ipiezo+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick
                        firstSolitaryPiezoAligndFoF{cutoffIdx,1}(:,oDend,ntrial) = groupAlign.dFoF{1,a}(ipiezo-postLick_frames:ipiezo+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick
                        firstSolitaryPiezoAlign{cutoffIdx,1}(:,ntrial) = tempMotionPeaks(ipiezo-postLick_frames:ipiezo+postLick_frames+postLick_frames-1,ntrial);% align licking data to first lick
                        firstSolitaryPiezoTrial{cutoffIdx,1}(:,ntrial) = itrial;
                        end
                    end
                end
                end
                %%%
                ind_pre = intersect(preRew_lickSearchRange,find(lickCueAlign(:,ntrial))); %finds every instance of a lick that occurs within the search window
                ind_post = intersect(postRew_lickSearchRange,find(lickCueAlign(:,ntrial))); %finds every instance of a lick that occurs [after the reward delivery]
                % piezo
                ind_preP = intersect(preRew_lickSearchRange,find(gMotionPeaks(:,ntrial))); %finds every instance of a lick that occurs within the search window
                ind_postP = intersect(postRew_lickSearchRange,find(gMotionPeaks(:,ntrial))); %finds every instance of a lick that occurs [after the reward delivery]
                %
                ind_all = intersect(lickSearchRange,find(lickCueAlign(:,ntrial)));
                preRew_lickBurstHz(:,ntrial) = length(ind_pre)./(cRctT/1000);
            if gBlock2(ntrial)==1
                preRewCSM_lickBurstHz(:,ntrial) = length(ind_pre)./(cRctT/1000);
                postRewCSM_lickBurstHz(:,ntrial) = length(ind_post)./(length(postRew_lickSearchRange)./frameRateHz);
            end
                postRew_lickBurstHz(:,ntrial) = length(ind_post)./(length(postRew_lickSearchRange)./frameRateHz);
            if gBlock2(ntrial)==0
                preRewCSP_lickBurstHz(:,ntrial) = length(ind_pre)./(cRctT/1000);
                postRewCSP_lickBurstHz(:,ntrial) = length(ind_post)./(length(postRew_lickSearchRange)./frameRateHz);
            end
            lickBurstHz_all(:,ntrial) = length(ind_all)./(length(lickSearchRange)./frameRateHz);
                ind = intersect(postRew_lickSearchRange,find(lickCueAlign(:,ntrial)));
             %%%similar idea as the above for loop, but instead aligns neural data as POST-REWARD lick events 
             %{
                for i = 1:length(ind)
                    ilick = ind(i);
                    if sum(lickCueAlign(ilick:ilick+lickSearch_frames-1,ntrial),1) >= 3
                        postRew_lickAlignEvents(:,oDend,ntrial) = groupAlign.events{1,a}(ilick-postLick_frames:ilick+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick after reward
                        postRew_lickAlignDFoF(:,oDend,ntrial) = groupAlign.dFoF{1,a}(ilick-postLick_frames:ilick+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick after reward
                        postRew_lickBurstStart(:,ntrial) = ilick; %array of every lick within post-reward window if 3+ licks were recorded for that trial
                        postRew_lickAlign(:,ntrial) = lickCueAlign(ilick-postLick_frames:ilick+postLick_frames+postLick_frames-1,ntrial);% align licking data to first lick
                        break
                    end
                end
                ind = intersect(preRew_lickSearchRange,find(lickCueAlign(:,ntrial)));
                for i = 1:length(ind)
                    ilick = ind(i);
                    if sum(lickCueAlign(ilick:ilick+lickSearch_frames-1,ntrial),1) >= 3
                        preRew_lickAlignEvents(:,oDend,ntrial) = groupAlign.events{1,a}(ilick-postLick_frames:ilick+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick after reward
                        preRew_lickAlignDFoF(:,oDend,ntrial) = groupAlign.dFoF{1,a}(ilick-postLick_frames:ilick+postLick_frames+postLick_frames-1,nDend,itrial);% align neural data to first lick after reward
                        preRew_lickBurstStart(:,ntrial) = ilick; %array of every lick within post-reward window if 3+ licks were recorded for that trial
                        preRew_lickAlign(:,ntrial) = lickCueAlign(ilick-postLick_frames:ilick+postLick_frames+postLick_frames-1,ntrial);% align licking data to first lick
                        break
                    end
                end
                %}
                ind_pre = find(lickCueAlign(prewin_frames:prewin_frames+rewDelay_frames,ntrial),1,'last'); %for orig: -500 to 1000ms around reward
                ind_preP = find(gMotionPeaks(prewin_frames:prewin_frames+rewDelay_frames,ntrial),1,'last');
         %ind_pre = find(lickCueAlign(prewin_frames-(rewDelay_frames+1):prewin_frames+rewDelay_frames,ntrial),1,'last'); %for novel: -1000 to 500ms around reward
%%%if there are licks recorded between cue and reward delivery (770 ms), align neural data to those licks 
                if ~isempty(ind_pre)
                    lastPreRewLickFrame(1,ntrial) = ind_pre;
                    preRewTrials(1,ntrial) =  1;
                    % align neural and licking data to the last lick
                    lastPreRew_lickAlignEvents(:,oDend,ntrial) = groupAlign.events{1,a}(ind_pre-postLick_frames+prewin_frames:ind_pre+postLick_frames+postLick_frames-1+prewin_frames,nDend,itrial); %for orig: -500 to 1000ms around reward
                    lastPreRew_lickAlignDFoF(:,oDend,ntrial) = groupAlign.dFoF{1,a}(ind_pre-postLick_frames+prewin_frames:ind_pre+postLick_frames+postLick_frames-1+prewin_frames,nDend,itrial); %for orig: -500 to 1000ms around reward
                    lastPreRew_lickAlign(:,ntrial) = lickCueAlign(ind_pre-postLick_frames+prewin_frames:ind_pre+postLick_frames+postLick_frames-1+prewin_frames,ntrial); %for orig: -500 to 1000ms around reward
                    lastPreRew_lickStart(:,ntrial) = ind_pre;

          %lastPreRew_lickAlignEvents(:,oDend,ntrial) = groupAlign.events{1,a}(ind_pre-postLick_frames+prewin_frames:ind_pre+postLick_frames+postLick_frames-1+prewin_frames,nDend,itrial); %for novel: -1000 to 500ms around reward
          %lastPreRew_lickAlign(:,ntrial) = lickCueAlign(ind_pre-postLick_frames+prewin_frames:ind_pre+postLick_frames+postLick_frames-1+prewin_frames,ntrial); %for novel: -1000 to 500ms around reward
                end
                
          % adding a section here to index trials with a final lick (for termination analysis
                ind_last = find(lickCueAlign(size(lickCueAlign,1)-6:size(lickCueAlign,1),ntrial),1);
                ind_lastP = find(gMotionPeaks(size(gMotionPeaks,1)-6:size(gMotionPeaks,1),ntrial),1);
                if isempty(ind_last)
                isolate_lastLick = find(lickCueAlign(prewin_frames+rewDelay_frames:end,ntrial),1,'last');
                if ~isempty(isolate_lastLick) && isolate_lastLick+prewin_frames+rewDelay_frames-1<size(lickCueAlign,1)-postLick_frames && isolate_lastLick>postLick_frames && sum(lickCueAlign(isolate_lastLick-lickSearch_frames:isolate_lastLick-1,ntrial),1) >= 3 %must be earlier than 135 otherwise I can't pull 15 frames after event
                isolate_lastLickInd(:,ntrial) = isolate_lastLick+rewDelay_frames;                 
                lastLick_lickAlignEvents(:,oDend,ntrial) = groupAlign.events{1,a}(isolate_lastLick-postLick_frames-postLick_frames+prewin_frames+rewDelay_frames:isolate_lastLick+postLick_frames-1+prewin_frames+rewDelay_frames,nDend,itrial); 
                lastLick_lickAlignDFoF(:,oDend,ntrial) = groupAlign.dFoF{1,a}(isolate_lastLick-postLick_frames-postLick_frames+prewin_frames+rewDelay_frames:isolate_lastLick+postLick_frames-1+prewin_frames+rewDelay_frames,nDend,itrial); 
                lastLick_lickAlign(:,ntrial) = lickCueAlign(isolate_lastLick-postLick_frames-postLick_frames+prewin_frames+rewDelay_frames:isolate_lastLick+postLick_frames-1+prewin_frames+rewDelay_frames,ntrial);
                end
                end
                if isempty(ind_lastP)
                isolate_lastMove = find(gMotionPeaks(prewin_frames+rewDelay_frames:end,ntrial),1,'last');
                if ~isempty(isolate_lastMove) && isolate_lastMove+prewin_frames+rewDelay_frames-1<size(gMotionPeaks,1)-postLick_frames
                isolate_lastMoveInd(:,ntrial) = isolate_lastMove;                 
                lastMove_piezoAlignEvents(:,oDend,ntrial) = groupAlign.events{1,a}(isolate_lastMove-postLick_frames-postLick_frames+prewin_frames:isolate_lastMove+postLick_frames-1+prewin_frames,nDend,itrial); 
                lastMove_piezoAlignDFoF(:,oDend,ntrial) = groupAlign.dFoF{1,a}(isolate_lastMove-postLick_frames-postLick_frames+prewin_frames:isolate_lastMove+postLick_frames-1+prewin_frames,nDend,itrial); 
                lastMove_piezoAlign(:,ntrial) = gMotionPeaks(isolate_lastMove-postLick_frames-postLick_frames+prewin_frames:isolate_lastMove+postLick_frames-1+prewin_frames,ntrial);
                end
                end
                
                %
                if ~isempty(ind_preP)
                    lastPreRewPiezoFrame(1,ntrial) = ind_preP;
                    preRewPiezoTrials(1,ntrial) = 1;
                    lastPreRew_piezoAlignEvents(:,oDend,ntrial) = groupAlign.events{1,a}(ind_preP-postLick_frames+prewin_frames:ind_preP+postLick_frames+postLick_frames-1+prewin_frames,nDend,itrial);
                    lastPreRew_piezoAlignDFoF(:,oDend,ntrial) = groupAlign.dFoF{1,a}(ind_preP-postLick_frames+prewin_frames:ind_preP+postLick_frames+postLick_frames-1+prewin_frames,nDend,itrial);
                    %lastPreRew_piezoAlign(:,ntrial) = gMotionPeaks(ind_preP-postLick_frames+prewin_frames:ind_preP+postLick_frames+postLick_frames-1+prewin_frames,ntrial);
                    lastPreRew_piezoAlign(:,ntrial) = groupAlign_piezo(ind_preP-postLick_frames+prewin_frames:ind_preP+postLick_frames+postLick_frames-1+prewin_frames,ntrial);
                    lastPreRew_piezoStart(:,ntrial) = ind_preP;

                end
                
             
            tTInx = tTInx+1;
            end
        end
    end
    ntrials=ntrial;
    end 
    
%% save

clearvars tempFig
if seshType ~= 2
save(fullfile(output_fn, '_groupwiseDataset.mat'));
else
save(fullfile(output_fn, [blockName '_groupwiseDataset.mat']));
end
%% piezo align formatting

ind_nan = find(isnan(groupAlign_piezo(1,:)));

g= [0.4660 0.6740 0.1880];
r= [0.6350 0.0780 0.1840];


ymax = (input('Set a y-axis limit for spike frequency graphs (5 is a good start): '));

%% %{
%% lick align figures%{
tt=gtt(1,:);
  %align licking data to cue
tempFig=setFigure; hold on;  
    shadedErrorBar_CV(tt, nanmean(lickCueAlign(:,~gBlock2),2).*(1000./frameRateHz), (nanstd(lickCueAlign(:,~gBlock2),[],2)./sqrt(unique(sum(~isnan(nanmean(lickCueAlign(:,~gBlock2),1)),2))).*(1000./frameRateHz)),'lineProps','k'); hold on;
    shadedErrorBar_CV(tt, nanmean(lickCueAlign(:,logical(gBlock2)),2).*(1000./frameRateHz), (nanstd(lickCueAlign(:,logical(gBlock2)),[],2)./sqrt(unique(sum(~isnan(nanmean(lickCueAlign(:,logical(gBlock2),1))),2))).*(1000./frameRateHz)),'lineProps','r');   
    scatter((lickBurstStart_postReward(:,~gBlock2)-prewin_frames).*(1000./frameRateHz), 10.*ones(size(lickBurstStart_postReward(:,~gBlock2))), 'k.');
    scatter((lickBurstStart_postReward(:,logical(gBlock2))-prewin_frames).*(1000./frameRateHz), 10.*ones(size(lickBurstStart_postReward(:,logical(gBlock2)))), 'r.');
    xlabel('Time from cue');
    ylabel('Lick rate (Hz)');
    title(['lick events from ' num2str(size(expt,2)) ' animals | cue aligned licking']);
    vline(0,'k');
    vline(767,'k');
    if seshType ~= 2
    %savefig(fullfile(output_fn, '_cueAlign_lickHz.fig'));
    saveas(tempFig, [output_fn  '_cueAlign_lickHz.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName '_cueAlign_lickHz.fig']));
    saveas(tempFig, [output_fn blockName '_cueAlign_lickHz.pdf'],'pdf');
    end
    
    nIC = size(groupAlign_events,2);
tt=gtt(1,:);
    tempFig=setFigure; hold on; %plot neural data of trials of early vs. late bursts
    shadedErrorBar_CV(tt,nanmean(nanmean(earlyBurstEvents,3),2).*(1000./frameRateHz), (nanstd(nanmean(earlyBurstEvents,3),[],2)./sqrt(nIC).*(1000./frameRateHz)),'lineProps','k');
    hold on;
    shadedErrorBar_CV(tt,nanmean(nanmean(lateBurstEvents,3),2).*(1000./frameRateHz), (nanstd(nanmean(lateBurstEvents,3),[],2)./sqrt(nIC).*(1000./frameRateHz)),'lineProps','b');
    errorbar(mean(((earlyBurstStart-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(earlyBurstEvents,3),2).*(1000./frameRateHz) nanmean(nanmean(lateBurstEvents,3),2).*(1000./frameRateHz)])),std(((earlyBurstStart-prewin_frames).*(1000./frameRateHz)./sqrt(size(gBlock2,2))),[],2,'omitnan'),'horizontal','k');
    errorbar(mean(((lateBurstStart-prewin_frames)),2,'omitnan').*(1000./frameRateHz),max(max([nanmean(nanmean(earlyBurstEvents,3),2).*(1000./frameRateHz) nanmean(nanmean(lateBurstEvents,3),2).*(1000./frameRateHz)])),std(((lateBurstStart-prewin_frames)./sqrt(size(gBlock2,2)).*(1000./frameRateHz)),[],2,'omitnan'),'horizontal','b');    
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    vline(0,'k');
    vline(767,'k');
    title([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | lick bursts: early (n = ' num2str(length(cell2mat(early_bst))) '| blk); late (n = ' num2str(length(cell2mat(late_bst))) '| blu)']);
    if seshType ~=2
    %savefig(fullfile(output_fn, 'earlyLate_10PercentOfLicks_allTrial.fig'));
    saveas(tempFig, [output_fn  'earlyLate_10PercentOfLicks_allTrial.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, ['\' blockName 'earlyLate_10PercentOfLicks_allTrial.fig']));
    saveas(tempFig, [output_fn blockName 'earlyLate_10PercentOfLicks_allTrial.pdf'],'pdf');
    end
    hold off;
    
     nIC = size(groupAlign_events,2);
tt=gtt(1,:); tt=tt(:,21:141);
    tempFig=setFigure; hold on; %plot neural data of trials of early vs. late bursts
    shadedErrorBar_CV(tt,nanmean(nanmean(earlyB2BurstEvents(21:141,:),3),2).*(1000./frameRateHz), (nanstd(nanmean(earlyB2BurstEvents(21:141,:),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    hold on;
    shadedErrorBar_CV(tt,nanmean(nanmean(lateB2BurstEvents(21:141,:),3),2).*(1000./frameRateHz), (nanstd(nanmean(lateB2BurstEvents(21:141,:),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
    errorbar(mean(((earlyB2BurstStart-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(earlyB2BurstEvents(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(lateB2BurstEvents(21:141,:),3),2).*(1000./frameRateHz)])),std(((earlyB2BurstStart-prewin_frames).*(1000./frameRateHz)./sqrt(size(earlyB2BurstStart,2))),[],2,'omitnan'),'horizontal','k');
    errorbar(mean(((lateB2BurstStart-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(earlyB2BurstEvents(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(lateB2BurstEvents(21:141,:),3),2).*(1000./frameRateHz)])),std(((lateB2BurstStart-prewin_frames).*(1000./frameRateHz)./sqrt(size(lateB2BurstStart,2))),[],2,'omitnan'),'horizontal','b');    
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    vline(0,'k');
    vline(767,'k');
    title([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | b2 lick bursts: early (n = ' num2str(sum(~isnan((earlyB2BurstStart)))) '| blk); late (n = ' num2str(sum(~isnan((lateB2BurstStart)))) '| blu)']);
    if seshType ~=2
    %savefig(fullfile(output_fn, 'earlyLate_25PercentOfLicks_CSMinus.fig'));
    saveas(tempFig, [output_fn  'earlyLate_25PercentOfLicks_CSMinus.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, ['\' blockName 'earlyLate_25PercentOfLicks_CSMinus.fig']));
    saveas(tempFig, [output_fn blockName 'earlyLate_25PercentOfLicks_CSMinus.pdf'],'pdf');
    end
    hold off;
    
     nIC = size(groupAlign_events,2);
tt=gtt(1,:); tt=tt(:,21:141);
    tempFig=setFigure; hold on; %plot neural data of trials of early vs. late bursts
    shadedErrorBar_CV(tt,nanmean(nanmean(earlyRewBurstEvents(21:141,:),3),2).*(1000./frameRateHz), (nanstd(nanmean(earlyRewBurstEvents(21:141,:),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    hold on;
    shadedErrorBar_CV(tt,nanmean(nanmean(lateRewBurstEvents(21:141,:),3),2).*(1000./frameRateHz), (nanstd(nanmean(lateRewBurstEvents(21:141,:),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
    errorbar(mean(((earlyRewBurstStart-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(earlyRewBurstEvents(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(lateRewBurstEvents(21:141,:),3),2).*(1000./frameRateHz)])),std(((earlyRewBurstStart-prewin_frames).*(1000./frameRateHz)./sqrt(size(earlyRewBurstStart,2))),[],2,'omitnan'),'horizontal','k');
    errorbar(mean(((lateRewBurstStart-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(earlyRewBurstEvents(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(lateRewBurstEvents(21:141,:),3),2).*(1000./frameRateHz)])),std(((lateRewBurstStart-prewin_frames).*(1000./frameRateHz)./sqrt(size(lateRewBurstStart,2))),[],2,'omitnan'),'horizontal','b');    
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    vline(0,'k');
    vline(767,'k');
    title([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | rew lick bursts: early (n = ' num2str(sum(~isnan((earlyRewBurstStart)))) '| blk); late (n = ' num2str(sum(~isnan((lateRewBurstStart)))) '| blu)']);
    if seshType ~=2
    %savefig(fullfile(output_fn, 'earlyLate_25PercentOfLicks_CSPlus.fig'));
    saveas(tempFig, [output_fn  'earlyLate_25PercentOfLicks_CSPlus.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, ['\' blockName 'earlyLate_25PercentOfLicks_CSPlus.fig']));
    saveas(tempFig, [output_fn blockName 'earlyLate_25PercentOfLicks_CSPlus.pdf'],'pdf');
    end
    hold off;
    
tt=gtt(1,:);
    pct_precue_burst = length(find((lickBurstStart_postReward-prewin_frames).*(1000./frameRateHz)<600))./size(lickBurstStart_postReward,2);
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
    
    tempFig=setFigure; hold on;
    subplot(2,1,1); hold on;% align neural activity to lick onset, only burst trials are included 
    shadedErrorBar_CV(tt(1,prewin_frames-postLick_frames+1:prewin_frames+postLick_frames+postLick_frames), nanmean(nanmean(postRew_lickAlignEvents,3),2).*(1000./frameRateHz), (nanstd(nanmean(postRew_lickAlignEvents,3),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    hold on;
    xlabel('Time from lick');
    ylabel('Spike rate (Hz)'); 
    ylim([0 inf]);
    title([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | ' num2str(sum(~isnan(postRew_lickBurstStart(1,:)))) ' trials post-reward lick bursts']);
    subplot(2,1,2); hold on; % seperate early burst trials vs. late burst trials
    shadedErrorBar_CV(tt(1,prewin_frames-postLick_frames+1:prewin_frames+postLick_frames+postLick_frames), nanmean(earlyBurstAlignEvents,2).*(1000./frameRateHz), (nanstd(earlyBurstAlignEvents,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    hold on;
    shadedErrorBar_CV(tt(1,prewin_frames-postLick_frames+1:prewin_frames+postLick_frames+postLick_frames), nanmean(lateBurstAlignEvents,2).*(1000./frameRateHz), (nanstd(lateBurstAlignEvents,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
    xlabel('Time from lick');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    title([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | ' num2str(floor(sum(cell2mat(nburst))/4)) ' early [blk]: avg = ' num2str(nanmean(cell2mat(early_bst_time))) ' ms; ' num2str(floor(sum(cell2mat(nburst))/4)) ' late [blu]: avg = ' num2str(nanmean(cell2mat(late_bst_time))) ' ms'])
    if seshType ~= 2
    %savefig(fullfile(output_fn, '_postRew_lickBurstAlignSpikingEvents.fig'))
    saveas(tempFig, [output_fn  '_postRew_lickBurstAlignSpikingEvents.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName '_postRew_lickBurstAlignSpikingEvents.fig']));
    saveas(tempFig, [output_fn blockName '_postRew_lickBurstAlignSpikingEvents.pdf'],'pdf');
    end
    
    tempFig=setFigure; hold on;
    subplot(2,1,1); hold on;% align licking data to first lick
    shadedErrorBar_CV(tt(1,prewin_frames-postLick_frames+1:prewin_frames+postLick_frames+postLick_frames), nanmean(postRew_lickAlign,2), (nanstd(postRew_lickAlign,[],2))./sqrt(unique(sum(~isnan(nanmean(postRew_lickAlign,1)),2))),'lineProps','k');
    hold on;
    xlabel('Time from lick');
    ylabel('Lick rate (Hz)');
    ylim([0 inf]);
    title('Neural align first lick post-reward')
    subplot(2,1,2); hold on;% plot early burst trials and late burst trials separately
    shadedErrorBar_CV(tt(1,prewin_frames-postLick_frames+1:prewin_frames+postLick_frames+postLick_frames), nanmean(earlyPostLickAlign,2), (nanstd(earlyPostLickAlign,[],2))./sqrt(length(cell2mat(ind_prerew_early_bst))),'lineProps','k');
    hold on;
    shadedErrorBar_CV(tt(1,prewin_frames-postLick_frames+1:prewin_frames+postLick_frames+postLick_frames), nanmean(latePostLickAlign,2), (nanstd(latePostLickAlign,[],2))./sqrt(length(cell2mat(ind_prerew_late_bst))),'lineProps','b');
    xlabel('Time from lick');
    ylabel('Lick rate (Hz)');
    ylim([0 inf]);
    title('Early burst trials [blk] and late burst trials [blu]')
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | post reward lick burst aligned spiking']);
    if seshType ~= 2
    %savefig(fullfile(output_fn, '_postRew_lickBurstAlignSpiking.fig'));
    saveas(tempFig, [output_fn  '_postRew_lickBurstAlignSpiking.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName '_postRew_lickBurstAlignSpiking.fig']));
    saveas(tempFig, [output_fn blockName '_postRew_lickBurstAlignSpiking.pdf'],'pdf');
    end
    hold off;
      
%%
%% MOST VS. LEAST MOVEMENT - TOP AND BOTTOM 25% OF TRIALS FOR EACH MOUSE
mostMovementAlignedEvents=[]; leastMovementAlignedEvents=[]; mostTrialMovement=[]; leastTrialMovement=[];
mostCSPMovementAlignedEvents=[]; leastCSPMovementAlignedEvents=[]; mostCSPTrialMovement=[]; leastCSPTrialMovement=[];
mostCSMMovementAlignedEvents=[]; leastCSMMovementAlignedEvents=[]; mostCSMTrialMovement=[]; leastCSMTrialMovement=[];

for mouseIdx = 1:exptCount
    block2Trials{mouseIdx}=[]; rewardTrials{mouseIdx}=[];
    CSPTrialMovement_ind=[]; CSMTrialMovement_ind=[];
    peakTrialChange=[]; preakTrialChange_ind=[]; trialChange=[];
for i = 1:size(groupAlign.piezo{mouseIdx},2)
for ii = 1:size(groupAlign.piezo{mouseIdx},1)-1
    trialChange(ii,i) = groupAlign.piezo{mouseIdx}(ii+1,i)-groupAlign.piezo{mouseIdx}(ii,i);
end
[peakTrialChange(1,i),peakTrialChange_ind(1,i)] = max(trialChange(50:81,i));
end
[TrialMovement, TrialMovement_ind] = sort(peakTrialChange);

for i = 1:size(block2{1,mouseIdx},2)
   if block2{1,mouseIdx}(1,i)==1; block2Trials{mouseIdx} = [block2Trials{mouseIdx} i]; end
   if block2{1,mouseIdx}(1,i)==0; rewardTrials{mouseIdx} = [rewardTrials{mouseIdx} i]; end
end
for i = 1:size(TrialMovement_ind,2)
    if ~isempty(find(rewardTrials{1,mouseIdx} == TrialMovement_ind(1,i)))
CSPTrialMovement_ind = [CSPTrialMovement_ind TrialMovement_ind(1,i)];
    elseif ~isempty(find(block2Trials{1,mouseIdx} == TrialMovement_ind(1,i)))
CSMTrialMovement_ind = [CSMTrialMovement_ind TrialMovement_ind(1,i)];
    end
end
mostCSPTrialMovement_ind{1,mouseIdx} = []; leastCSPTrialMovement_ind{1,mouseIdx} = [];
mostCSMTrialMovement_ind{1,mouseIdx} = []; leastCSMTrialMovement_ind{1,mouseIdx} = [];
mostTrialMovement_ind{1,mouseIdx} = TrialMovement_ind(:,size(groupAlign.piezo{1,mouseIdx},2)-(floor(size(groupAlign.piezo{1,mouseIdx},2)./4)):end);
leastTrialMovement_ind{1,mouseIdx} = TrialMovement_ind(:,1:round(size(groupAlign.piezo{1,mouseIdx},2)./4));

if seshType == 2 && sum(gBlock2) > 0
else
mostCSPTrialMovement_ind{1,mouseIdx} = CSPTrialMovement_ind(:,size(groupAlign.piezo{1,mouseIdx}(:,~block2{1,mouseIdx}),2)-(floor(size(groupAlign.piezo{1,mouseIdx}(:,~block2{1,mouseIdx}),2)./4)):size(groupAlign.piezo{1,mouseIdx}(:,~block2{1,mouseIdx}),2));
leastCSPTrialMovement_ind{1,mouseIdx} = CSPTrialMovement_ind(:,1:ceil(size(groupAlign.piezo{1,mouseIdx}(:,~block2{1,mouseIdx}),2)./4));
mostCSPMovementAlignedEvents = [mostCSPMovementAlignedEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,mostCSPTrialMovement_ind{1,mouseIdx}),3)];
leastCSPMovementAlignedEvents = [leastCSPMovementAlignedEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,leastCSPTrialMovement_ind{1,mouseIdx}),3)];
mostCSPTrialMovement = [mostCSPTrialMovement groupAlign.piezo{1,mouseIdx}(:,mostCSPTrialMovement_ind{1,mouseIdx})];
leastCSPTrialMovement = [leastCSPTrialMovement groupAlign.piezo{1,mouseIdx}(:,leastCSPTrialMovement_ind{1,mouseIdx})];

end
if seshType == 2 && sum(gBlock2) == 0
elseif expType == 3 && cue ~= 3
else
mostCSMTrialMovement_ind{1,mouseIdx} = CSMTrialMovement_ind(:,size(groupAlign.piezo{1,mouseIdx}(:,logical(block2{1,mouseIdx})),2)-(floor(size(groupAlign.piezo{1,mouseIdx}(:,logical(block2{1,mouseIdx})),2)./4)):size(groupAlign.piezo{1,mouseIdx}(:,logical(block2{1,mouseIdx})),2));
leastCSMTrialMovement_ind{1,mouseIdx} = CSMTrialMovement_ind(:,1:ceil(size(groupAlign.piezo{1,mouseIdx}(:,logical(block2{1,mouseIdx})),2)./4));
mostCSMMovementAlignedEvents = [mostCSMMovementAlignedEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,mostCSMTrialMovement_ind{1,mouseIdx}),3)];
leastCSMMovementAlignedEvents = [leastCSMMovementAlignedEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,leastCSMTrialMovement_ind{1,mouseIdx}),3)];
mostCSMTrialMovement = [mostCSMTrialMovement groupAlign.piezo{1,mouseIdx}(:,mostCSMTrialMovement_ind{1,mouseIdx})];
leastCSMTrialMovement = [leastCSMTrialMovement groupAlign.piezo{1,mouseIdx}(:,leastCSMTrialMovement_ind{1,mouseIdx})];

end
%
mostMovementAlignedEvents = [mostMovementAlignedEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,mostTrialMovement_ind{1,mouseIdx}),3)];
leastMovementAlignedEvents = [leastMovementAlignedEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,leastTrialMovement_ind{1,mouseIdx}),3)];
mostTrialMovement = [mostTrialMovement groupAlign.piezo{1,mouseIdx}(:,mostTrialMovement_ind{1,mouseIdx})];
leastTrialMovement = [leastTrialMovement groupAlign.piezo{1,mouseIdx}(:,leastTrialMovement_ind{1,mouseIdx})];
end

for i = 1:size(animalTrials,2)
mostMotion_csp_ind{1,i} = zeros(1,animalTrials(1,i)); mostMotion_csm_ind{1,i} = zeros(1,animalTrials(1,i));
leastMotion_csp_ind{1,i} = zeros(1,animalTrials(1,i)); leastMotion_csm_ind{1,i} = zeros(1,animalTrials(1,i));
earlyLick_csp_ind{1,i} = zeros(1,animalTrials(1,i)); earlyLick_csm_ind{1,i} = zeros(1,animalTrials(1,i));
lateLick_csp_ind{1,i} = zeros(1,animalTrials(1,i)); lateLick_csm_ind{1,i} = zeros(1,animalTrials(1,i));
    
mostMotion_csp_ind{1,i}(1,mostCSPTrialMovement_ind{1,i}) = 1; mostMotion_csm_ind{1,i}(1,mostCSMTrialMovement_ind{1,i}) = 1;
leastMotion_csp_ind{1,i}(1,leastCSPTrialMovement_ind{1,i}) = 1; leastMotion_csm_ind{1,i}(1,leastCSMTrialMovement_ind{1,i}) = 1;
earlyLick_csp_ind{1,i}(1,earlyCSPLick_ind{1,i}) = 1; earlyLick_csm_ind{1,i}(1,earlyCSMLick_ind{1,i}) = 1;
lateLick_csp_ind{1,i}(1,lateCSPLick_ind{1,i}) = 1; lateLick_csm_ind{1,i}(1,lateCSMLick_ind{1,i}) = 1;
end

earlyMostStart=[]; earlyLeastStart=[]; lateMostStart=[]; lateLeastStart=[];
for i = 1:size(animalTrials,2)
earlyLick_mostMotion_csp_Cspk{1,i} = []; lateLick_mostMotion_csp_Cspk{1,i} = []; earlyLick_leastMotion_csp_Cspk{1,i} = []; lateLick_leastMotion_csp_Cspk{1,i} = [];
earlyLick_mostMotion_csm_Cspk{1,i} = []; lateLick_mostMotion_csm_Cspk{1,i} = []; earlyLick_leastMotion_csm_Cspk{1,i} = []; lateLick_leastMotion_csm_Cspk{1,i} = [];

earlyLick_mostMotion_csp_Cspk{1,i} = [earlyLick_mostMotion_csp_Cspk{1,i} reshape(mean(groupAlign.events{1,i}(:,:,(mostMotion_csp_ind{1,i}==1 & earlyLick_csp_ind{1,i}==1)),2,'omitnan'),[150 sum(mostMotion_csp_ind{1,i}==1 & earlyLick_csp_ind{1,i}==1)])];
lateLick_mostMotion_csp_Cspk{1,i} = [lateLick_mostMotion_csp_Cspk{1,i} reshape(mean(groupAlign.events{1,i}(:,:,(mostMotion_csp_ind{1,i}==1 & lateLick_csp_ind{1,i}==1)),2,'omitnan'),[150 sum(mostMotion_csp_ind{1,i}==1 & lateLick_csp_ind{1,i}==1)])];
earlyLick_leastMotion_csp_Cspk{1,i} = [earlyLick_leastMotion_csp_Cspk{1,i} reshape(mean(groupAlign.events{1,i}(:,:,(leastMotion_csp_ind{1,i}==1 & earlyLick_csp_ind{1,i}==1)),2,'omitnan'),[150 sum(leastMotion_csp_ind{1,i}==1 & earlyLick_csp_ind{1,i}==1)])]; 
lateLick_leastMotion_csp_Cspk{1,i} = [lateLick_leastMotion_csp_Cspk{1,i} reshape(mean(groupAlign.events{1,i}(:,:,(leastMotion_csp_ind{1,i}==1 & lateLick_csp_ind{1,i}==1)),2,'omitnan'),[150 sum(leastMotion_csp_ind{1,i}==1 & lateLick_csp_ind{1,i}==1)])];    

earlyLick_mostMotion_csm_Cspk{1,i} = [earlyLick_mostMotion_csm_Cspk{1,i} reshape(mean(groupAlign.events{1,i}(:,:,(mostMotion_csm_ind{1,i}==1 & earlyLick_csm_ind{1,i}==1)),2,'omitnan'),[150 sum(mostMotion_csm_ind{1,i}==1 & earlyLick_csm_ind{1,i}==1)])];
lateLick_mostMotion_csm_Cspk{1,i} = [lateLick_mostMotion_csm_Cspk{1,i} reshape(mean(groupAlign.events{1,i}(:,:,(mostMotion_csm_ind{1,i}==1 & lateLick_csm_ind{1,i}==1)),2,'omitnan'),[150 sum(mostMotion_csm_ind{1,i}==1 & lateLick_csm_ind{1,i}==1)])];
earlyLick_leastMotion_csm_Cspk{1,i} = [earlyLick_leastMotion_csm_Cspk{1,i} reshape(mean(groupAlign.events{1,i}(:,:,(leastMotion_csm_ind{1,i}==1 & earlyLick_csm_ind{1,i}==1)),2,'omitnan'),[150 sum(leastMotion_csm_ind{1,i}==1 & earlyLick_csm_ind{1,i}==1)])]; 
lateLick_leastMotion_csm_Cspk{1,i} =[lateLick_leastMotion_csm_Cspk{1,i} reshape(mean(groupAlign.events{1,i}(:,:,(leastMotion_csm_ind{1,i}==1 & lateLick_csm_ind{1,i}==1)),2,'omitnan'),[150 sum(leastMotion_csm_ind{1,i}==1 & lateLick_csm_ind{1,i}==1)])];    

tempRewStart = groupAlign.lickStart{1,i}(:,~block2{1,i} & ~isnan(groupAlign.lickStart{1,i}));
earlyMostStart = [earlyMostStart (tempRewStart(:,mostMotion_csp_ind{1,i}==1 & earlyLick_csp_ind{1,i}==1))];
earlyLeastStart = [earlyLeastStart (tempRewStart(:,leastMotion_csp_ind{1,i}==1 & earlyLick_csp_ind{1,i}==1))];
lateMostStart = [lateMostStart (tempRewStart(:,mostMotion_csp_ind{1,i}==1 & lateLick_csp_ind{1,i}==1))];
lateLeastStart = [lateLeastStart (tempRewStart(:,leastMotion_csp_ind{1,i}==1 & lateLick_csp_ind{1,i}==1))];
end
%%
    
ndends=0; ntrials=0; lastPreLickEvent=[]; firstPostLickEvent=[];  
csPluslastPreLickEvent=[]; csPlusfirstPostLickEvent=[]; csMinuslastPreLickEvent=[]; csMinusfirstPostLickEvent=[];
csPlusfirstPreLickEvent=[]; csMinusfirstPreLickEvent=[]; csPlusLickTerminationEvent=[]; csMinusLickTerminationEvent=[];

%for most\least motion alt
mostMotion_csp_all = cell2mat(mostMotion_csp_ind); csPlusfirstPostLickEvent_leastMotion=[]; csPlusfirstPreLickEvent_leastMotion=[];
leastMotion_csp_all = cell2mat(leastMotion_csp_ind); csPlusfirstPostLickEvent_mostMotion=[]; csPlusfirstPreLickEvent_mostMotion=[];
for mouseIdx = 1:exptCount
    nDend = (1:nsize(1,mouseIdx))+ndends;
    nTrial = (1:animalTrials(1,mouseIdx))+ntrials; 
               lastPreLickEvent = [lastPreLickEvent (nanmean(lastPreRew_lickAlignEvents(:,nDend,nTrial),3))]; 
               firstPostLickEvent = [firstPostLickEvent (nanmean(firstPostRew_lickAlignEvents(:,nDend,nTrial),3))];
                csPluslastPreLickEvent = [csPluslastPreLickEvent (nanmean(lastPreRew_lickAlignEvents(:,nDend,nTrial(:,~gBlock2(:,nTrial))),3))];
                csPlusfirstPostLickEvent = [csPlusfirstPostLickEvent (nanmean(firstPostRew_lickAlignEvents(:,nDend,nTrial(:,~gBlock2(:,nTrial))),3))];
                csMinuslastPreLickEvent = [csMinuslastPreLickEvent (nanmean(lastPreRew_lickAlignEvents(:,nDend,nTrial(:,gBlock2(:,nTrial))),3))];
                csMinusfirstPostLickEvent = [csMinusfirstPostLickEvent (nanmean(firstPostRew_lickAlignEvents(:,nDend,nTrial(:,gBlock2(:,nTrial))),3))];
                csPlusfirstPreLickEvent = [csPlusfirstPreLickEvent (nanmean(firstPreRew_lickAlignEvents(:,nDend,nTrial(:,~gBlock2(:,nTrial))),3))];
                csMinusfirstPreLickEvent = [csMinusfirstPreLickEvent (nanmean(firstPreRew_lickAlignEvents(:,nDend,nTrial(:,gBlock2(:,nTrial))),3))];
                csPlusLickTerminationEvent = [csPlusLickTerminationEvent (nanmean(lastLick_lickAlignEvents(:,nDend,nTrial(:,~gBlock2(:,nTrial))),3))];
                csMinusLickTerminationEvent = [csMinusLickTerminationEvent (nanmean(lastLick_lickAlignEvents(:,nDend,nTrial(:,gBlock2(:,nTrial))),3))];
 
                %getting data for most\least motion
                csPlusfirstPostLickEvent_leastMotion = [csPlusfirstPostLickEvent_leastMotion (nanmean(firstPostRew_lickAlignEvents(:,nDend,nTrial(:,~gBlock2(:,nTrial) & leastMotion_csp_all(1,nTrial)==1)),3))];
                csPlusfirstPostLickEvent_mostMotion = [csPlusfirstPostLickEvent_mostMotion (nanmean(firstPostRew_lickAlignEvents(:,nDend,nTrial(:,~gBlock2(:,nTrial) & mostMotion_csp_all(1,nTrial)==1)),3))];
                csPlusfirstPreLickEvent_leastMotion = [csPlusfirstPreLickEvent_leastMotion (nanmean(firstPreRew_lickAlignEvents(:,nDend,nTrial(:,~gBlock2(:,nTrial) & leastMotion_csp_all(1,nTrial)==1)),3))];
                csPlusfirstPreLickEvent_mostMotion = [csPlusfirstPreLickEvent_mostMotion (nanmean(firstPreRew_lickAlignEvents(:,nDend,nTrial(:,~gBlock2(:,nTrial) & mostMotion_csp_all(1,nTrial)==1)),3))];
                
                
    ntrials = max(nTrial);
    ndends = max(nDend);
end

    tl_rew = (1-postLick_frames:(postLick_frames*2)).*(1000./frameRateHz); %orig :: -500ms to 1000 ms around reward
    %tl_rew = (1-postLick_frames*2:(postLick_frames)).*(1000./frameRateHz); %novel :: -1000ms to 500ms around reward
    tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
    subplot(2,2,1); hold on;
    shadedErrorBar_CV(tl_rew, nanmean(lastPreLickEvent,2).*(1000./frameRateHz), (nanstd(lastPreLickEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    hold on;
    title('Last lick before reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(lastPreRew_lickAlign,2), (nanstd(lastPreRew_lickAlign,[],2))./sqrt(unique(sum(~isnan(nanmean(lastPreRew_lickAlign,1))))),'lineProps','k');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title([num2str(sum(preRewTrials)) 'trials with pre-reward licks']);
    subplot(2,2,2);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(firstPostLickEvent,2).*(1000./frameRateHz), (nanstd(firstPostLickEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    title('First lick after reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(firstPostRew_lickAlign,2), (nanstd(firstPostRew_lickAlign,[],2))./sqrt(unique(sum(~isnan(nanmean(firstPostRew_lickAlign,1))))),'lineProps','k');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title([num2str(sum(postRewTrials)) 'trials with post-reward licks']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Licks relative to reward']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, '_lastVsFirstLick.fig'));
    saveas(tempFig, [output_fn  '_lastVsFirstLick.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName '_lastVsFirstLick.fig']));
    saveas(tempFig, [output_fn blockName '_lastVsFirstLick.pdf'],'pdf');
    end

 
    
ndends=0; ntrials=0; lastPrePiezoEvent=[]; firstPostPiezoEvent=[];  
csPluslastPrePiezoEvent=[]; csPlusfirstPostPiezoEvent=[]; csMinuslastPrePiezoEvent=[]; csMinusfirstPostPiezoEvent=[];
csPlusfirstPrePiezoEvent=[]; csMinusfirstPrePiezoEvent=[]; csPlusLastPiezoEvent=[]; csMinusLastPiezoEvent=[];
for mouseIdx = 1:exptCount
    nDend = (1:nsize(1,mouseIdx))+ndends;
    nTrial = (1:animalTrials(1,mouseIdx))+ntrials; 
               lastPrePiezoEvent = [lastPrePiezoEvent (nanmean(lastPreRew_piezoAlignEvents(:,nDend,nTrial),3))]; 
               firstPostPiezoEvent = [firstPostPiezoEvent (nanmean(firstPostRew_piezoAlignEvents(:,nDend,nTrial),3))];
                csPluslastPrePiezoEvent = [csPluslastPrePiezoEvent (nanmean(lastPreRew_piezoAlignEvents(:,nDend,nTrial(:,~gBlock2(:,nTrial))),3))];
                csPlusfirstPostPiezoEvent = [csPlusfirstPostPiezoEvent (nanmean(firstPostRew_piezoAlignEvents(:,nDend,nTrial(:,~gBlock2(:,nTrial))),3))];
                csMinuslastPrePiezoEvent = [csMinuslastPrePiezoEvent (nanmean(lastPreRew_piezoAlignEvents(:,nDend,nTrial(:,gBlock2(:,nTrial))),3))];
                csMinusfirstPostPiezoEvent = [csMinusfirstPostPiezoEvent (nanmean(firstPostRew_piezoAlignEvents(:,nDend,nTrial(:,gBlock2(:,nTrial))),3))];
                csPlusLastPiezoEvent = [csPlusLastPiezoEvent (nanmean(lastMove_piezoAlignEvents(:,nDend,nTrial(:,~gBlock2(:,nTrial))),3))];
                csMinusLastPiezoEvent = [csMinusLastPiezoEvent (nanmean(lastMove_piezoAlignEvents(:,nDend,nTrial(:,gBlock2(:,nTrial))),3))];
                
                csPlusfirstPrePiezoEvent = [csPlusfirstPrePiezoEvent (nanmean(firstPreRew_piezoAlignEvents(:,nDend,nTrial(:,~gBlock2(:,nTrial))),3))];
                csMinusfirstPrePiezoEvent = [csMinusfirstPrePiezoEvent (nanmean(firstPreRew_piezoAlignEvents(:,nDend,nTrial(:,gBlock2(:,nTrial))),3))];
ntrials = max(nTrial);
    ndends = max(nDend);
end    
    

    tempFig=setFigure; hold on; % align neural and Piezoing data to the first and last Piezo, respectively
    subplot(2,2,1); hold on;
    shadedErrorBar_CV(tl_rew, nanmean(lastPrePiezoEvent,2).*(1000./frameRateHz), (nanstd(lastPrePiezoEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    hold on;
    title('Last movement before reward');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(lastPreRew_piezoAlign,2).*(1000./frameRateHz), (nanstd(lastPreRew_piezoAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(nanmean(lastPreRew_piezoAlign,1))))),'lineProps','k');
    xlabel('Time from movement (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 10]);
    title([num2str(sum(preRewPiezoTrials)) 'trials with pre-reward movement']);
    subplot(2,2,2);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(firstPostPiezoEvent,2).*(1000./frameRateHz), (nanstd(firstPostPiezoEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    title('First movement after reward');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(firstPostRew_piezoAlign,2).*(1000./frameRateHz), (nanstd(firstPostRew_piezoAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(nanmean(firstPostRew_piezoAlign,1))))),'lineProps','k');
    xlabel('Time from movement (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 10]);
    title([num2str(sum(postRewPiezoTrials)) 'trials with post-reward movement']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Movement relative to reward']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, '_lastVsFirstPiezo.fig'));
    saveas(tempFig, [output_fn  '_lastVsFirstPiezo.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName '_lastVsFirstPiezo.fig']));
    saveas(tempFig, [output_fn blockName '_lastVsFirstPiezo.pdf'],'pdf');
    end

    %% rapid deviation to seperate first and last lick into CS+ and CS- trials
    csPlusFirstPostRew_lickAlign = firstPostRew_lickAlign(:,~gBlock2);
    csMinusFirstPostRew_lickAlign = firstPostRew_lickAlign(:,gBlock2);
    csPlusLastPreRew_lickAlign = lastPreRew_lickAlign(:,~gBlock2);
    csMinusLastPreRew_lickAlign = lastPreRew_lickAlign(:,gBlock2);

    csPlusFirstPreRew_lickAlign = firstPreRew_lickAlign(:,~gBlock2);
    csMinusFirstPreRew_lickAlign = firstPreRew_lickAlign(:,gBlock2);
    csPlusFpreRewTrials = preFRewTrials(:,~gBlock2);
    csMinusFpreRewTrials = preFRewTrials(:,gBlock2);

    csPluspostRewTrials = postRewTrials(:,~gBlock2);
    csMinuspostRewTrials = postRewTrials(:,gBlock2);
    csPluspreRewTrials = preRewTrials(:,~gBlock2);
    csMinuspreRewTrials = preRewTrials(:,gBlock2);
    %CS+
        tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
    subplot(2,2,1);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csPluslastPreLickEvent,2).*(1000./frameRateHz), (nanstd(csPluslastPreLickEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    hold on;
    title('Last lick before reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csPlusLastPreRew_lickAlign,2), (nanstd(csPlusLastPreRew_lickAlign,[],2))./sqrt(unique(sum(~isnan(nanmean(csPlusLastPreRew_lickAlign,1))))),'lineProps','k');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title([num2str(sum(csPluspreRewTrials)) 'CS+ trials with pre-reward licks']);
    subplot(2,2,2);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csPlusfirstPostLickEvent,2).*(1000./frameRateHz), (nanstd(csPlusfirstPostLickEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    title('First lick after reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csPlusFirstPostRew_lickAlign,2), (nanstd(csPlusFirstPostRew_lickAlign,[],2))./sqrt(unique(sum(~isnan(nanmean(csPlusFirstPostRew_lickAlign,1))))),'lineProps','k');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title([num2str(sum(csPluspostRewTrials)) 'CS+ trials with post-reward licks']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Licks relative to reward']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, '_csPluslastVsFirstLick.fig'));
    saveas(tempFig, [output_fn  '_csPluslastVsFirstLick.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName '_csPluslastVsFirstLick.fig']));
    saveas(tempFig, [output_fn blockName '_csPluslastVsFirstLick.pdf'],'pdf');
    end

    %CS-
        tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
    subplot(2,2,1);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csMinuslastPreLickEvent,2).*(1000./frameRateHz), (nanstd(csMinuslastPreLickEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    hold on;
    title('Last lick before reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csMinusLastPreRew_lickAlign,2), (nanstd(csMinusLastPreRew_lickAlign,[],2))./sqrt(unique(sum(~isnan(nanmean(lastPreRew_lickAlign,1))))),'lineProps','k');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title([num2str(sum(csMinuspreRewTrials)) 'CS- trials with pre-reward licks']);
    subplot(2,2,2);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csMinusfirstPostLickEvent,2).*(1000./frameRateHz), (nanstd(csMinusfirstPostLickEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    title('First lick after reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csMinusFirstPostRew_lickAlign,2), (nanstd(csMinusFirstPostRew_lickAlign,[],2))./sqrt(unique(sum(~isnan(nanmean(firstPostRew_lickAlign,1))))),'lineProps','k');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title([num2str(sum(csMinuspostRewTrials)) 'CS- trials with post-reward licks']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Licks relative to reward']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, '_csMinuslastVsFirstLick.fig'));
    saveas(tempFig, [output_fn  '_csMinuslastVsFirstLick.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName '_csMinuslastVsFirstLick.fig']));
    saveas(tempFig, [output_fn blockName '_csMinuslastVsFirstLick.pdf'],'pdf');
    end

    %% dead last licks - i.e. lick termination
    csPlusLastLick_lickAlign = lastLick_lickAlign(:,~gBlock2);
    csMinusLastLick_lickAlign = lastLick_lickAlign(:,gBlock2);
    csPlusLast_PiezoAlign = lastMove_piezoAlign(:,~gBlock2);
    csMinusLast_PiezoAlign = lastMove_piezoAlign(:,gBlock2);
      tl_last = (1-postLick_frames*2:(postLick_frames)).*(1000./frameRateHz); %orig :: -500ms to 1000 ms around reward

        tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
    subplot(2,2,1);hold on;
    shadedErrorBar_CV(tl_last, nanmean(csPlusLickTerminationEvent,2).*(1000./frameRateHz), (nanstd(csPlusLickTerminationEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    hold on;
    title('Last lick after reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    shadedErrorBar_CV(tl_last, nanmean(csPlusLastLick_lickAlign,2), (nanstd(csPlusLastLick_lickAlign,[],2))./sqrt(unique(sum(~isnan(nanmean(csPlusLastLick_lickAlign,1))))),'lineProps','k');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title([num2str(sum(~isnan(mean(csPlusLastLick_lickAlign,1,'omitnan')))) 'CS+ trials with terminating lick']);
    subplot(2,2,2);hold on;
    shadedErrorBar_CV(tl_last, nanmean(csMinusLickTerminationEvent,2).*(1000./frameRateHz), (nanstd(csMinusLickTerminationEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','r');
    title('Last lick after reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    shadedErrorBar_CV(tl_last, nanmean(csMinusLastLick_lickAlign,2), (nanstd(csMinusLastLick_lickAlign,[],2))./sqrt(unique(sum(~isnan(nanmean(csMinusLastLick_lickAlign,1))))),'lineProps','r');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title([num2str(sum(~isnan(mean(csMinusLastLick_lickAlign,1,'omitnan')))) 'CS- trials with terminating lick']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Licks relative to reward']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, '_lastLickCSMinusVsCSPlus.fig'));
    saveas(tempFig, [output_fn  '_lastLickCSMinusVsCSPlus.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName '_lastLickCSMinusVsCSPlus.fig']));
    saveas(tempFig, [output_fn blockName '_lastLickCSMinusVsCSPlus.pdf'],'pdf');
    end
    
    %motion
       tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
    subplot(2,2,1);hold on;
    shadedErrorBar_CV(tl_last, nanmean(csPlusLastPiezoEvent,2).*(1000./frameRateHz), (nanstd(csPlusLastPiezoEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    hold on;
    title('Last movement after reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    shadedErrorBar_CV(tl_last, nanmean(csPlusLast_PiezoAlign,2), (nanstd(csPlusLast_PiezoAlign,[],2))./sqrt(unique(sum(~isnan(nanmean(csPlusLast_PiezoAlign,1))))),'lineProps','k');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title([num2str(sum(~isnan(mean(csPlusLast_PiezoAlign,1,'omitnan')))) 'CS+ trials with terminating motion']);
    subplot(2,2,2);hold on;
    shadedErrorBar_CV(tl_last, nanmean(csMinusLastPiezoEvent,2).*(1000./frameRateHz), (nanstd(csMinusLastPiezoEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','r');
    title('Last movement after reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    shadedErrorBar_CV(tl_last, nanmean(csMinusLast_PiezoAlign,2), (nanstd(csMinusLast_PiezoAlign,[],2))./sqrt(unique(sum(~isnan(nanmean(csMinusLast_PiezoAlign,1))))),'lineProps','r');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title([num2str(sum(~isnan(mean(csMinusLast_PiezoAlign,1,'omitnan')))) 'CS- trials with terminating motion']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Movement relative to reward']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, '_lastMotionCSMinusVsCSPlus.fig'));
    saveas(tempFig, [output_fn  '_lastMotionCSMinusVsCSPlus.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName '_lastMotionCSMinusVsCSPlus.fig']));
    saveas(tempFig, [output_fn blockName '_lastMotionCSMinusVsCSPlus.pdf'],'pdf');
    end
    
    %%
    %quick addition here - lick align to trials with most and least movement
    if seshType~=2 || seshType==2 && cue==1
    csPlusFirstPostRew_lickAlign_mostMotion = firstPostRew_lickAlign(:,~gBlock2 & logical(mostMotion_csp_all));
    csPlusFirstPostRew_lickAlign_leastMotion = firstPostRew_lickAlign(:,~gBlock2 & logical(leastMotion_csp_all));
    
    csPlusFirstPreRew_lickAlign_mostMotion = firstPreRew_lickAlign(:,~gBlock2 & logical(mostMotion_csp_all));
    csPlusFirstPreRew_lickAlign_leastMotion = firstPreRew_lickAlign(:,~gBlock2 & logical(leastMotion_csp_all));
    
    
    
     tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
    subplot(2,2,1);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csPlusfirstPostLickEvent_mostMotion,2).*(1000./frameRateHz), (nanstd(csPlusfirstPostLickEvent_mostMotion,[],2).*(1000./frameRateHz))./sqrt(size(csPlusfirstPostLickEvent_mostMotion,2)),'lineProps','k');
    hold on;
    title('First lick after reward during trials with most movement');
    xlabel('Time from lick (ms)');
    ylabel('Cspk rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csPlusFirstPostRew_lickAlign_mostMotion,2), (nanstd(csPlusFirstPostRew_lickAlign_mostMotion,[],2))./sqrt(unique(sum(~isnan(nanmean(csPlusFirstPostRew_lickAlign_mostMotion,1))))),'lineProps','k');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title([num2str(size(csPlusFirstPostRew_lickAlign_mostMotion,2)) 'CS+ trials with large motion and post-reward licks']);
    subplot(2,2,2);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csPlusfirstPostLickEvent_leastMotion,2).*(1000./frameRateHz), (nanstd(csPlusfirstPostLickEvent_leastMotion,[],2).*(1000./frameRateHz))./sqrt(size(csPlusfirstPostLickEvent_leastMotion,2)),'lineProps','b');
    title('First lick after reward during trials with least movement');
    xlabel('Time from lick (ms)');
    ylabel('Cspk rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csPlusFirstPostRew_lickAlign_leastMotion,2), (nanstd(csPlusFirstPostRew_lickAlign_leastMotion,[],2))./sqrt(unique(sum(~isnan(nanmean(csPlusFirstPostRew_lickAlign_leastMotion,1))))),'lineProps','b');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title([num2str(size(csPlusFirstPostRew_lickAlign_leastMotion,2)) 'CS+ trials with small motion and post-reward licks']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Licks relative to reward']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, '_csPlus_lickAlign_postRew_mostAndLeastMovement.fig'));
    saveas(tempFig, [output_fn  '_csPlus_lickAlign_postRew_mostAndLeastMovement.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName '_csPlus_lickAlign_postRew_mostAndLeastMovement.fig']));
    saveas(tempFig, [output_fn blockName '_csPlus_lickAlign_postRew_mostAndLeastMovement.pdf'],'pdf');
    end
    
    
      
     tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
    subplot(2,2,1);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csPlusfirstPreLickEvent_mostMotion,2).*(1000./frameRateHz), (nanstd(csPlusfirstPreLickEvent_mostMotion,[],2).*(1000./frameRateHz))./sqrt(size(csPlusfirstPreLickEvent_mostMotion,2)),'lineProps','k');
    hold on;
    title('First lick after cue during trials with most movement');
    xlabel('Time from lick (ms)');
    ylabel('Cspk rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csPlusFirstPreRew_lickAlign_mostMotion,2), (nanstd(csPlusFirstPreRew_lickAlign_mostMotion,[],2))./sqrt(unique(sum(~isnan(nanmean(csPlusFirstPreRew_lickAlign_mostMotion,1))))),'lineProps','k');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title([num2str(size(csPlusFirstPreRew_lickAlign_mostMotion,2)) 'CS+ trials with large motion and pre-reward licks']);
    subplot(2,2,2);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csPlusfirstPreLickEvent_leastMotion,2).*(1000./frameRateHz), (nanstd(csPlusfirstPreLickEvent_leastMotion,[],2).*(1000./frameRateHz))./sqrt(size(csPlusfirstPreLickEvent_leastMotion,2)),'lineProps','b');
    title('First lick after cue during trials with least movement');
    xlabel('Time from lick (ms)');
    ylabel('Cspk rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csPlusFirstPreRew_lickAlign_leastMotion,2), (nanstd(csPlusFirstPreRew_lickAlign_leastMotion,[],2))./sqrt(unique(sum(~isnan(nanmean(csPlusFirstPreRew_lickAlign_leastMotion,1))))),'lineProps','b');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title([num2str(size(csPlusFirstPreRew_lickAlign_leastMotion,2)) 'CS+ trials with small motion and pre-reward licks']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Licks relative to reward']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, '_csPlus_lickAlign_preRew_mostAndLeastMovement.fig'));
    saveas(tempFig, [output_fn  '_csPlus_lickAlign_preRew_mostAndLeastMovement.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName '_csPlus_lickAlign_preRew_mostAndLeastMovement.fig']));
    saveas(tempFig, [output_fn blockName '_csPlus_lickAlign_preRew_mostAndLeastMovement.pdf'],'pdf');
    end
    end
    % second section moved to 2800
    %% Crazy shit real quick - this script is gonna need some work cause it's disgustingly long
    % Anyway, going to make align plots for CS+ and CS- 3rds to go with the
    % later 3rds figures, but will only 'currently' be making these for the
    % rule-flip naive day where changes appear over time in the neural data

    
    tt=gtt(1,:); 
    repTrials = [1:20; 21:40; 41:60];
    clearvars ndends ntrials
    ndends=[0 nsize];  ntrials=[0 animalTrials];  
    
for rep=1:3
csPluslastPreLickRep{rep}=[]; csPlusfirstPostLickRep{rep}=[]; 
csMinuslastPreLickRep{rep}=[]; csMinusfirstPostLickRep{rep}=[]; 
csMinusfirstPreLickRep{rep}=[]; csPlusfirstPreLickRep{rep}=[];
csPluslastPrePiezoEventRep{rep}=[]; csPlusfirstPostPiezoEventRep{rep}=[]; 
csMinuslastPrePiezoEventRep{rep}=[]; csMinusfirstPostPiezoEventRep{rep}=[]; 
csPlusfirstPrePiezoEventRep{rep}=[]; csMinusfirstPrePiezoEventRep{rep}=[]; 
piezoPlusRep{rep}=[]; piezoMinusRep{rep}=[]; lickPlusRep{rep}=[]; lickMinusRep{rep}=[];
remainderP = []; remainderM=[];
for mouseIdx = 1:exptCount
endL(mouseIdx,rep) = 0;
cspT{mouseIdx} = find(~block2{mouseIdx}); 
csmT{mouseIdx} = find(block2{mouseIdx});
    nDend = 1:nsize(1,mouseIdx)+ndends(1,mouseIdx);
    nTrial = 1:animalTrials(1,mouseIdx)+ntrials(1,mouseIdx); 
    if  seshType~=2 && length(cspT{mouseIdx})>=60 && length(csmT{mouseIdx})>=60 || seshType==2 && length(cspT{mouseIdx})>=60 || seshType==2 && length(cspT{mouseIdx})>=60 || rep<3
        thirdsP{mouseIdx} = [(1:floor(length(cspT{mouseIdx})./3)); (floor(length(cspT{mouseIdx})./3)+1:floor(length(cspT{mouseIdx})./3).*2); ((floor(length(cspT{mouseIdx})./3).*2)+1:floor(length(cspT{mouseIdx})./3).*3)];
        thirdsM{mouseIdx} = [(1:floor(length(csmT{mouseIdx})./3)); (floor(length(csmT{mouseIdx})./3)+1:floor(length(csmT{mouseIdx})./3).*2); ((floor(length(csmT{mouseIdx})./3).*2)+1:floor(length(csmT{mouseIdx})./3).*3)];
        if seshType~=2 || cue==1; remainderP{mouseIdx} = length(cspT{mouseIdx})-thirdsP{mouseIdx}(3,end); end
        if seshType~=2 || cue==2; remainderM{mouseIdx} = length(csmT{mouseIdx})-thirdsM{mouseIdx}(3,end); end 
                csPluslastPreLickRep{rep} = [csPluslastPreLickRep{rep} (nanmean(lastPreRew_lickAlignEvents(:,nDend,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,:))),3))];
                csPlusfirstPostLickRep{rep} = [csPlusfirstPostLickRep{rep} (nanmean(firstPostRew_lickAlignEvents(:,nDend,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,:))),3))];
                csMinuslastPreLickRep{rep} = [csMinuslastPreLickRep{rep} (nanmean(lastPreRew_lickAlignEvents(:,nDend,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,:))),3))];
                csMinusfirstPostLickRep{rep} = [csMinusfirstPostLickRep{rep} (nanmean(firstPostRew_lickAlignEvents(:,nDend,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,:))),3))];
                csPlusfirstPreLickRep{rep} = [csPlusfirstPreLickRep{rep} (nanmean(firstPreRew_lickAlignEvents(:,nDend,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,:))),3))];
                csMinusfirstPreLickRep{rep} = [csMinusfirstPreLickRep{rep} (nanmean(firstPreRew_lickAlignEvents(:,nDend,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,:))),3))];

                csPluslastPrePiezoEventRep{rep} = [csPluslastPrePiezoEventRep{rep} (nanmean(lastPreRew_piezoAlignEvents(:,nDend,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,:))),3))];
                csPlusfirstPostPiezoEventRep{rep} = [csPlusfirstPostPiezoEventRep{rep} (nanmean(firstPostRew_piezoAlignEvents(:,nDend,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,:))),3))];
                csMinuslastPrePiezoEventRep{rep} = [csMinuslastPrePiezoEventRep{rep} (nanmean(lastPreRew_piezoAlignEvents(:,nDend,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,:))),3))];
                csMinusfirstPostPiezoEventRep{rep} = [csMinusfirstPostPiezoEventRep{rep} (nanmean(firstPostRew_piezoAlignEvents(:,nDend,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,:))),3))];
                csPlusfirstPrePiezoEventRep{rep} = [csPlusfirstPrePiezoEventRep{rep} (nanmean(firstPreRew_piezoAlignEvents(:,nDend,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,:))),3))];
                csMinusfirstPrePiezoEventRep{rep} = [csMinusfirstPrePiezoEventRep{rep} (nanmean(firstPreRew_piezoAlignEvents(:,nDend,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,:))),3))];

    csPluspostRewTrialsRep{rep,mouseIdx} = postRewTrials(:,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,:)));
    csMinuspostRewTrialsRep{rep,mouseIdx} = postRewTrials(:,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,:)));
    csPluspreRewTrialsRep{rep,mouseIdx} = preRewTrials(:,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,:)));
    csMinuspreRewTrialsRep{rep,mouseIdx} = preRewTrials(:,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,:)));
    
    csPluspostRewPiezoTrialsRep{rep,mouseIdx} = postRewPiezoTrials(:,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,:)));
    csMinuspostRewPiezoTrialsRep{rep,mouseIdx} = postRewPiezoTrials(:,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,:)));
    csPluspreRewPiezoTrialsRep{rep,mouseIdx} = preRewPiezoTrials(:,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,:)));
    csMinuspreRewPiezoTrialsRep{rep,mouseIdx} = preRewPiezoTrials(:,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,:)));
    
    csPlusFpreRewTrialsRep{rep,mouseIdx} = preFRewTrials(:,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,:)));
    csMinusFpreRewTrialsRep{rep,mouseIdx} = preFRewTrials(:,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,:)));
    csPluspreFRewPiezoTrialsRep{rep,mouseIdx} = preFRewPiezoTrials(:,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,:)));
    csMinuspreFRewPiezoTrialsRep{rep,mouseIdx} = preFRewPiezoTrials(:,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,:)));
    else 
        if seshType~=2; endL(mouseIdx,rep) = min([length(cspT{mouseIdx}) length(csmT{mouseIdx})]); end
        if seshType==2 && cue==1; endL(mouseIdx,rep) = length(cspT{mouseIdx}); elseif seshType==2 && cue==2; endL(mouseIdx,rep) = length(csmT{mouseIdx}); end
        if seshType~=2 || seshType==2 && cue==1
                csPluslastPreLickRep{rep} = [csPluslastPreLickRep{rep} (nanmean(lastPreRew_lickAlignEvents(:,nDend,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,1):endL(mouseIdx,rep))),3))];
                csPlusfirstPostLickRep{rep} = [csPlusfirstPostLickRep{rep} (nanmean(firstPostRew_lickAlignEvents(:,nDend,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,1):endL(mouseIdx,rep))),3))];
                csPlusfirstPreLickRep{rep} = [csPlusfirstPreLickRep{rep} (nanmean(firstPreRew_lickAlignEvents(:,nDend,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,1):endL(mouseIdx,rep))),3))];
                csPluslastPrePiezoEventRep{rep} = [csPluslastPrePiezoEventRep{rep} (nanmean(lastPreRew_piezoAlignEvents(:,nDend,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,1):endL(mouseIdx,rep))),3))];
                csPlusfirstPostPiezoEventRep{rep} = [csPlusfirstPostPiezoEventRep{rep} (nanmean(firstPostRew_piezoAlignEvents(:,nDend,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,1):endL(mouseIdx,rep))),3))];
                csPlusfirstPrePiezoEventRep{rep} = [csPlusfirstPrePiezoEventRep{rep} (nanmean(firstPreRew_piezoAlignEvents(:,nDend,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,1):endL(mouseIdx,rep))),3))];
    csPluspostRewTrialsRep{rep,mouseIdx} = postRewTrials(:,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,1):endL(mouseIdx,rep)));
    csPluspreRewTrialsRep{rep,mouseIdx} = preRewTrials(:,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,1):endL(mouseIdx,rep)));
    csPluspostRewPiezoTrialsRep{rep,mouseIdx} = postRewPiezoTrials(:,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,1):endL(mouseIdx,rep)));
    csPluspreRewPiezoTrialsRep{rep,mouseIdx} = preRewPiezoTrials(:,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,1):endL(mouseIdx,rep)));
    csPlusFpreRewTrialsRep{rep,mouseIdx} = preFRewTrials(:,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,1):endL(mouseIdx,rep)));
    csPluspreFRewPiezoTrialsRep{rep,mouseIdx} = preFRewPiezoTrials(:,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,1):endL(mouseIdx,rep)));
        end
        if seshType~=2 || seshType==2 && cue==2
                csMinuslastPreLickRep{rep} = [csMinuslastPreLickRep{rep} (nanmean(lastPreRew_lickAlignEvents(:,nDend,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,1):endL(mouseIdx,rep))),3))];
                csMinusfirstPostLickRep{rep} = [csMinusfirstPostLickRep{rep} (nanmean(firstPostRew_lickAlignEvents(:,nDend,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,1):endL(mouseIdx,rep))),3))];
                csMinusfirstPreLickRep{rep} = [csMinusfirstPreLickRep{rep} (nanmean(firstPreRew_lickAlignEvents(:,nDend,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,1):endL(mouseIdx,rep))),3))];
                csMinuslastPrePiezoEventRep{rep} = [csMinuslastPrePiezoEventRep{rep} (nanmean(lastPreRew_piezoAlignEvents(:,nDend,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,1):endL(mouseIdx,rep))),3))];
                csMinusfirstPostPiezoEventRep{rep} = [csMinusfirstPostPiezoEventRep{rep} (nanmean(firstPostRew_piezoAlignEvents(:,nDend,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,1):endL(mouseIdx,rep))),3))];
                csMinusfirstPrePiezoEventRep{rep} = [csMinusfirstPrePiezoEventRep{rep} (nanmean(firstPreRew_piezoAlignEvents(:,nDend,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,1):endL(mouseIdx,rep))),3))];        
    csMinuspostRewPiezoTrialsRep{rep,mouseIdx} = postRewPiezoTrials(:,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,1):endL(mouseIdx,rep)));
    csMinuspreRewPiezoTrialsRep{rep,mouseIdx} = preRewPiezoTrials(:,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,1):endL(mouseIdx,rep)));
    csMinuspostRewTrialsRep{rep,mouseIdx} = postRewTrials(:,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,1):endL(mouseIdx,rep)));
    csMinuspreRewTrialsRep{rep,mouseIdx} = preRewTrials(:,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,1):endL(mouseIdx,rep)));
    csMinusFpreRewTrialsRep{rep,mouseIdx} = preFRewTrials(:,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,1):endL(mouseIdx,rep)));
    csMinuspreFRewPiezoTrialsRep{rep,mouseIdx} = preFRewPiezoTrials(:,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,1):endL(mouseIdx,rep)));
        end
    end
    
    if  seshType~=2 && length(cspT{mouseIdx})>=60 && length(csmT{mouseIdx})>=60 || seshType==2 && length(cspT{mouseIdx})>=60 || seshType==2 && length(cspT{mouseIdx})>=60 || rep<3
    piezoPlusRep{rep} = [piezoPlusRep{rep} groupAlign.piezo{1,mouseIdx}(:,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,:)))]; 
    piezoMinusRep{rep} = [piezoMinusRep{rep} groupAlign.piezo{1,mouseIdx}(:,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,:)))];
    lickPlusRep{rep} = [lickPlusRep{rep} groupAlign.licks{1,mouseIdx}(:,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,:)))];
    lickMinusRep{rep} = [lickMinusRep{rep} groupAlign.licks{1,mouseIdx}(:,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,:)))];
    else
    if seshType~=2 || seshType==2 && cue==1
    piezoPlusRep{rep} = [piezoPlusRep{rep} groupAlign.piezo{1,mouseIdx}(:,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,1):endL(mouseIdx,rep)))];
    lickPlusRep{rep} = [lickPlusRep{rep} groupAlign.licks{1,mouseIdx}(:,cspT{mouseIdx}(1,thirdsP{mouseIdx}(rep,1):endL(mouseIdx,rep)))];
    end
    if seshType~=2 || seshType==2 && cue==2
    piezoMinusRep{rep} = [piezoMinusRep{rep} groupAlign.piezo{1,mouseIdx}(:,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,1):endL(mouseIdx,rep)))];
    lickMinusRep{rep} = [lickMinusRep{rep} groupAlign.licks{1,mouseIdx}(:,csmT{mouseIdx}(1,thirdsM{mouseIdx}(rep,1):endL(mouseIdx,rep)))];
    end
    end
end
    
    if seshType~=2 || seshType==2 && cue==1
    csPluspostRewTrialsRepT{rep} = cell2mat(csPluspostRewTrialsRep(rep,:));
    csPluspreRewTrialsRepT{rep} = cell2mat(csPluspreRewTrialsRep(rep,:));
    csPluspostRewPiezoTrialsRepT{rep} = cell2mat(csPluspostRewPiezoTrialsRep(rep,:));
    csPluspreRewPiezoTrialsRepT{rep} = cell2mat(csPluspreRewPiezoTrialsRep(rep,:));   
    csPlusFpreRewTrialsRepT{rep} = cell2mat(csPlusFpreRewTrialsRep(rep,:));
    csPluspreFRewPiezoTrialsRepT{rep} = cell2mat(csPluspreFRewPiezoTrialsRep(rep,:));
    end
    if seshType~=2 || seshType==2 && cue==2
    csMinuspostRewTrialsRepT{rep} = cell2mat(csMinuspostRewTrialsRep(rep,:));
    csMinuspreRewTrialsRepT{rep} = cell2mat(csMinuspreRewTrialsRep(rep,:));
    csMinuspostRewPiezoTrialsRepT{rep} = cell2mat(csMinuspostRewPiezoTrialsRep(rep,:));
    csMinuspreRewPiezoTrialsRepT{rep} = cell2mat(csMinuspreRewPiezoTrialsRep(rep,:));
    csMinusFpreRewTrialsRepT{rep} = cell2mat(csMinusFpreRewTrialsRep(rep,:));
    csMinuspreFRewPiezoTrialsRepT{rep} = cell2mat(csMinuspreFRewPiezoTrialsRep(rep,:));
    end
    
end

    %temp
    csPlusfirstPostRew_piezoAlign = firstPostRew_piezoAlign(:,~gBlock2);
    csMinusfirstPostRew_piezoAlign = firstPostRew_piezoAlign(:,gBlock2);
    csPluslastPreRew_piezoAlign = lastPreRew_piezoAlign(:,~gBlock2);
    csMinuslastPreRew_piezoAlign = lastPreRew_piezoAlign(:,gBlock2);
    csPlusfirstPreRew_piezoAlign = firstPreRew_piezoAlign(:,~gBlock2);
    csMinusfirstPreRew_piezoAlign = firstPreRew_piezoAlign(:,gBlock2);
    %temp
colorsRepCSM={[63/255 0 0],[128/255 0 0],[255/255 0 0]};
colorsRepCSP={[0 0 0],[76.5/255 76.5/255 76.5/255],[153/255 153/255 153/255]};
    
    if seshType~=2 || seshType==2 && cue==1
         tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
    subplot(2,2,1);hold on;
    for rep=1:3
    shadedErrorBar_CV(tl_rew, nanmean(csPluslastPreLickRep{rep},2).*(1000./frameRateHz), (nanstd(csPluslastPreLickRep{rep},[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps',{'color',colorsRepCSP{rep}});
    hold on;
    end
    title('Last lick before reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    for rep=1:3
    shadedErrorBar_CV(tt(1,30:74), nanmean(lickPlusRep{rep}(30:74,:),2), (nanstd(lickPlusRep{rep}(30:74,:),[],2))./sqrt(size(lickPlusRep{rep}(37:81,:),2)),'lineProps',{'color',colorsRepCSP{rep}});
    hold on;
    end
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 .5]);
    title([num2str(sum(csPluspreRewTrialsRepT{rep})) 'CS+ trials with pre-reward licks' '-' num2str(rep)]);
    subplot(2,2,2);hold on;
    for rep = 1:3
    shadedErrorBar_CV(tl_rew, nanmean(csPlusfirstPostLickRep{rep},2).*(1000./frameRateHz), (nanstd(csPlusfirstPostLickRep{rep},[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps',{'color',colorsRepCSP{rep}});
    hold on;
    end
    title('First lick after reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    for rep=1:3
    shadedErrorBar_CV(tt(1,74:118), nanmean(lickPlusRep{rep}(74:118,:),2), (nanstd(lickPlusRep{rep}(74:118,:),[],2))./sqrt(size(lickPlusRep{rep}(37:81,:),2)),'lineProps',{'color',colorsRepCSP{rep}});
    hold on;
    end
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 .5]);
    title([num2str(sum(csPluspostRewTrialsRepT{rep}))  'CS+ trials with post-reward licks' '-' num2str(rep)]);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Licks relative to reward']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, [num2str(rep) '_csPluslastVsFirstLick.fig']));
    saveas(tempFig, [output_fn num2str(rep) '_csPluslastVsFirstLick.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [num2str(rep) blockName '_csPluslastVsFirstLick.fig']));
    saveas(tempFig, [output_fn num2str(rep) blockName '_csPluslastVsFirstLick.pdf'],'pdf');
    end

   tempFig=setFigure; hold on; % align neural and Piezoing data to the first and last Piezo, respectively
    subplot(2,2,1); hold on;
    for rep = 1:3
    shadedErrorBar_CV(tl_rew, nanmean(csPluslastPrePiezoEventRep{rep},2).*(1000./frameRateHz), (nanstd(csPluslastPrePiezoEventRep{rep},[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps',{'color',colorsRepCSP{rep}});
    hold on;
    end
    hold on;
    title('Last movement before reward');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    for rep = 1:3
    shadedErrorBar_CV(tt(1,30:74), nanmean(piezoPlusRep{rep}(30:74,:),2).*(1000./frameRateHz), (nanstd(piezoPlusRep{rep}(30:74,:),[],2).*(1000./frameRateHz))./sqrt(size(piezoPlusRep{rep}(37:81,:),2)),'lineProps',{'color',colorsRepCSP{rep}});
    hold on;
    end
    xlabel('Time from movement (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 10]);
    title([num2str(sum(csPluspreRewPiezoTrialsRepT{rep})) 'CS+ trials with pre-reward movement -' num2str(rep)]);
    subplot(2,2,2);hold on;
    for rep = 1:3
    shadedErrorBar_CV(tl_rew, nanmean(csPlusfirstPostPiezoEventRep{rep},2).*(1000./frameRateHz), (nanstd(csPlusfirstPostPiezoEventRep{rep},[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps',{'color',colorsRepCSP{rep}});
    hold on;
    end
    title('First movement after reward');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    for rep = 1:3
    shadedErrorBar_CV(tt(1,74:118), nanmean(piezoPlusRep{rep}(74:118,:),2).*(1000./frameRateHz), (nanstd(piezoPlusRep{rep}(74:118,:),[],2).*(1000./frameRateHz))./sqrt(size(piezoPlusRep{rep}(37:81,:),2)),'lineProps',{'color',colorsRepCSP{rep}});
    hold on;
    end
    xlabel('Time from movement (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 10]);
    title([num2str(sum(csPluspostRewPiezoTrialsRepT{rep})) 'CS+ trials with post-reward movement - ' num2str(rep)]);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Movement relative to reward']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn,  [num2str(rep) '_csPluslastVsFirstPiezo.fig']));
    saveas(tempFig, [output_fn num2str(rep) '_csPluslastVsFirstPiezo.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [num2str(rep) blockName '_csPluslastVsFirstPiezo.fig']));
    saveas(tempFig, [output_fn num2str(rep) blockName '_csPluslastVsFirstPiezo.pdf'],'pdf');
    end
    
        %CS+
    tempFig=setFigure; hold on; % align neural and Piezoing data to the first and last Piezo, respectively
    subplot(2,2,1); hold on;
    for rep = 1:3
    shadedErrorBar_CV(tl_rew, nanmean(csPlusfirstPrePiezoEventRep{rep},2).*(1000./frameRateHz), (nanstd(csPlusfirstPrePiezoEventRep{rep},[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps',{'color',colorsRepCSP{rep}});
    hold on;
    end
    hold on;
    title('First movement after cue');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    for rep = 1:3
    shadedErrorBar_CV(tl_rew, nanmean(piezoPlusRep{rep}(37:81,:),2).*(1000./frameRateHz), (nanstd(piezoPlusRep{rep}(37:81,:),[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(nanmean(piezoPlusRep{rep}(37:81,:),1))))),'lineProps',{'color',colorsRepCSP{rep}});
    hold on;
    end
    xlabel('Time from movement (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 10]);
    title([num2str(sum(csPluspreFRewPiezoTrialsRepT{rep})) 'CS+ trials with post-cue movement - ' num2str(rep)]);
    subplot(2,2,2);hold on;
    for rep = 1:3
    shadedErrorBar_CV(tl_rew, nanmean(csPlusfirstPreLickRep{rep},2).*(1000./frameRateHz), (nanstd(csPlusfirstPreLickRep{rep},[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps',{'color',colorsRepCSP{rep}});
    hold on;
    end
    title('First lick after cue');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    for rep = 1:3
    shadedErrorBar_CV(tl_rew, nanmean(csPlusFirstPreRew_lickAlign,2), (nanstd(csPlusFirstPreRew_lickAlign,[],2))./sqrt(unique(sum(~isnan(nanmean(csPlusFirstPreRew_lickAlign,1))))),'lineProps',{'color',colorsRepCSP{rep}});
    hold on;
    end
    xlabel('Time from lick (ms)');
    ylabel('Lick Rate (mV)');
    ylim([0 1]);
    title([num2str(sum(csPlusFpreRewTrialsRepT{rep})) 'CS+ trials with post-cue licks - ' num2str(rep)]);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Movement and Licks directly following cue']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, [num2str(rep) '_csPlusfirstPreRewPiezovsLick.fig']));
    saveas(tempFig, [output_fn num2str(rep) '_csPlusfirstPreRewPiezovsLick.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [num2str(rep) blockName '_csPlusFirstPreRewPiezovsLick.fig']));
    saveas(tempFig, [output_fn num2str(rep) blockName '_csPlusFirstPreRewPiezovsLick.pdf'],'pdf');
    end
    end
    
    %CS-
    if seshType~=2 || seshType==2 && cue==2
        tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
    subplot(2,2,1);hold on;
    for rep = 1:3
    shadedErrorBar_CV(tl_rew, nanmean(csMinuslastPreLickRep{rep},2).*(1000./frameRateHz), (nanstd(csMinuslastPreLickRep{rep},[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps',{'color',colorsRepCSM{rep}});
    hold on;
    end
    title('Last lick before reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    for rep = 1:3
    shadedErrorBar_CV(tt(1,30:74), nanmean(lickMinusRep{rep}(30:74,:),2), (nanstd(lickMinusRep{rep}(30:74,:),[],2))./sqrt(size(lickMinusRep{rep}(37:81,:),2)),'lineProps',{'color',colorsRepCSM{rep}});
    hold on;
    end
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title([num2str(sum(csMinuspreRewTrialsRepT{rep}))  'CS- trials with pre-reward licks' '-' num2str(rep)]);
    subplot(2,2,2);hold on;
    for rep = 1:3
    shadedErrorBar_CV(tl_rew, nanmean(csMinusfirstPostLickRep{rep},2).*(1000./frameRateHz), (nanstd(csMinusfirstPostLickRep{rep},[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps',{'color',colorsRepCSM{rep}});
    hold on;
    end
    title('First lick after reward');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    for rep = 1:3
    shadedErrorBar_CV(tt(1,74:118), nanmean(lickMinusRep{rep}(74:118,:),2), (nanstd(lickMinusRep{rep}(74:118,:),[],2))./sqrt(size(lickMinusRep{rep}(37:81,:),2)),'lineProps',{'color',colorsRepCSM{rep}});
    hold on;
    end
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title([num2str(sum(csMinuspostRewTrialsRepT{rep})) 'CS- trials with post-reward licks' '-' num2str(rep)]);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Licks relative to reward']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, [num2str(rep) '_csMinuslastVsFirstLick.fig']));
    saveas(tempFig, [output_fn num2str(rep) '_csMinuslastVsFirstLick.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [num2str(rep) blockName '_csMinuslastVsFirstLick.fig']));
    saveas(tempFig, [output_fn num2str(rep) blockName '_csMinuslastVsFirstLick.pdf'],'pdf');
    end
    
    %CS-
    tempFig=setFigure; hold on; % align neural and Piezoing data to the first and last Piezo, respectively
    subplot(2,2,1); hold on;
    for rep = 1:3
    shadedErrorBar_CV(tl_rew, nanmean(csMinuslastPrePiezoEventRep{rep},2).*(1000./frameRateHz), (nanstd(csMinuslastPrePiezoEventRep{rep},[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps',{'color',colorsRepCSM{rep}});
    hold on;
    end
    hold on;
    title('Last movement before reward');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    for rep = 1:3
    shadedErrorBar_CV(tt(1,30:74), nanmean(piezoMinusRep{rep}(30:74,:),2).*(1000./frameRateHz), (nanstd(piezoMinusRep{rep}(30:74,:),[],2).*(1000./frameRateHz))./sqrt(size(piezoMinusRep{rep}(37:81,:),2)),'lineProps',{'color',colorsRepCSM{rep}});
    hold on;
    end
    xlabel('Time from movement (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 10]);
    title([num2str(sum(csMinuspreRewPiezoTrialsRepT{rep})) 'CS- trials with pre-reward movement - ' num2str(rep)]);
    subplot(2,2,2);hold on;
    for rep = 1:3
    shadedErrorBar_CV(tl_rew, nanmean(csMinusfirstPostPiezoEventRep{rep},2).*(1000./frameRateHz), (nanstd(csMinusfirstPostPiezoEventRep{rep},[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps',{'color',colorsRepCSM{rep}});
    hold on;
    end
    title('First movement after reward');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    for rep = 1:3
    shadedErrorBar_CV(tt(1,74:118), nanmean(piezoMinusRep{rep}(74:118,:),2).*(1000./frameRateHz), (nanstd(piezoMinusRep{rep}(74:118,:),[],2).*(1000./frameRateHz))./sqrt(size(piezoMinusRep{rep}(37:81,:),2)),'lineProps',{'color',colorsRepCSM{rep}});
    hold on;
    end
    xlabel('Time from movement (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 10]);
    title([num2str(sum(csMinuspostRewPiezoTrialsRepT{rep})) 'CS- trials with post-reward movement - ' num2str(rep)]);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Movement relative to reward']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, [num2str(rep) '_csMinuslastVsFirstPiezo.fig']));
    saveas(tempFig, [output_fn num2str(rep) '_csMinuslastVsFirstPiezo.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [num2str(rep) blockName '_csMinuslastVsFirstPiezo.fig']));
    saveas(tempFig, [output_fn num2str(rep) blockName '_csMinuslastVsFirstPiezo.pdf'],'pdf');
    end
    %%%

    %CS-
    tempFig=setFigure; hold on; % align neural and Piezoing data to the first and last Piezo, respectively
    subplot(2,2,1); hold on;
    for rep = 1:3
    shadedErrorBar_CV(tl_rew, nanmean(csMinusfirstPrePiezoEventRep{rep},2).*(1000./frameRateHz), (nanstd(csMinusfirstPrePiezoEventRep{rep},[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps',{'color',colorsRepCSM{rep}});
    hold on;
    end
    title('First movement after cue');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    for rep = 1:3
    shadedErrorBar_CV(tl_rew, nanmean(csMinusfirstPreRew_piezoAlign,2).*(1000./frameRateHz), (nanstd(csMinusfirstPreRew_piezoAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(nanmean(csMinusfirstPreRew_piezoAlign,1))))),'lineProps',{'color',colorsRepCSM{rep}});
    hold on;
    end
    xlabel('Time from movement (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 10]);
    title([num2str(sum(csMinuspreFRewPiezoTrialsRepT{rep})) 'CS- trials with post-cue movement - ' num2str(rep)]);
    subplot(2,2,2);hold on;
    for rep = 1:3
    shadedErrorBar_CV(tl_rew, nanmean(csMinusfirstPreLickRep{rep},2).*(1000./frameRateHz), (nanstd(csMinusfirstPreLickRep{rep},[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps',{'color',colorsRepCSM{rep}});
    hold on;
    end
    title('First licks after cue');
    xlabel('Time from licks (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    for rep = 1:3
    shadedErrorBar_CV(tl_rew, nanmean(csMinusFirstPreRew_lickAlign,2), (nanstd(csMinusFirstPreRew_lickAlign,[],2))./sqrt(unique(sum(~isnan(nanmean(csMinusFirstPreRew_lickAlign,1))))),'lineProps',{'color',colorsRepCSM{rep}});
    hold on;
    end
    xlabel('Time from licks (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title([num2str(sum(csMinusFpreRewTrialsRepT{rep})) 'CS- trials with post-cue licks - ' num2str(rep)]);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Movement and Licks directly following cue']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, [num2str(rep) '_csMinusFirstPreRewPiezovsLick.fig']));
    saveas(tempFig, [output_fn num2str(rep) '_csMinusFirstPreRewPiezovsLick.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [num2str(rep) blockName '_csMinusFirstPreRewPiezovsLick.fig']));
    saveas(tempFig, [output_fn num2str(rep) blockName '_csMinusFirstPreRewPiezovsLick.pdf'],'pdf');
    end
    end
    
    %% CS+ and CS- for spike rate aligned to first and last movement (2o above u) in reference to reward delivery
    csPlusfirstPostRew_piezoAlign = firstPostRew_piezoAlign(:,~gBlock2);
    csMinusfirstPostRew_piezoAlign = firstPostRew_piezoAlign(:,gBlock2);
    csPluslastPreRew_piezoAlign = lastPreRew_piezoAlign(:,~gBlock2);
    csMinuslastPreRew_piezoAlign = lastPreRew_piezoAlign(:,gBlock2);
    
    csPlusfirstPreRew_piezoAlign = firstPreRew_piezoAlign(:,~gBlock2);
    csMinusfirstPreRew_piezoAlign = firstPreRew_piezoAlign(:,gBlock2);
    csPluspreFRewPiezoTrials = preFRewPiezoTrials(:,~gBlock2);
    csMinuspreFRewPiezoTrials = preFRewPiezoTrials(:,gBlock2);
    
    csPluspostRewPiezoTrials = postRewPiezoTrials(:,~gBlock2);
    csMinuspostRewPiezoTrials = postRewPiezoTrials(:,gBlock2);
    csPluspreRewPiezoTrials = preRewPiezoTrials(:,~gBlock2);
    csMinuspreRewPiezoTrials = preRewPiezoTrials(:,gBlock2);
    if cue==3 || cue==1
    %CS+
     tempFig=setFigure; hold on; % align neural and Piezoing data to the first and last Piezo, respectively
    subplot(2,2,1); hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csPluslastPrePiezoEvent,2).*(1000./frameRateHz), (nanstd(csPluslastPrePiezoEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    hold on;
    title('Last movement before reward');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csPluslastPreRew_piezoAlign,2).*(1000./frameRateHz), (nanstd(csPluslastPreRew_piezoAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(nanmean(csPluslastPreRew_piezoAlign,1))))),'lineProps','k');
    xlabel('Time from movement (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 10]);
    title([num2str(sum(csPluspreRewPiezoTrials)) 'CS+ trials with pre-reward movement']);
    subplot(2,2,2);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csPlusfirstPostPiezoEvent,2).*(1000./frameRateHz), (nanstd(csPlusfirstPostPiezoEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    title('First movement after reward');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csPlusfirstPostRew_piezoAlign,2).*(1000./frameRateHz), (nanstd(csPlusfirstPostRew_piezoAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(nanmean(csPlusfirstPostRew_piezoAlign,1))))),'lineProps','k');
    xlabel('Time from movement (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 10]);
    title([num2str(sum(csPluspostRewPiezoTrials)) 'CS+ trials with post-reward movement']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Movement relative to reward']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, '_csPluslastVsFirstPiezo.fig'));
    saveas(tempFig, [output_fn  '_csPluslastVsFirstPiezo.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName '_csPluslastVsFirstPiezo.fig']));
    saveas(tempFig, [output_fn blockName '_csPluslastVsFirstPiezo.pdf'],'pdf');
    end
    else
    end
    
    if cue==3 || cue==2    
    %CS-
    tempFig=setFigure; hold on; % align neural and Piezoing data to the first and last Piezo, respectively
    subplot(2,2,1); hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csMinuslastPrePiezoEvent,2).*(1000./frameRateHz), (nanstd(csMinuslastPrePiezoEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','r');
    hold on;
    title('Last movement before reward');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csMinuslastPreRew_piezoAlign,2).*(1000./frameRateHz), (nanstd(csMinuslastPreRew_piezoAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(nanmean(csMinuslastPreRew_piezoAlign,1))))),'lineProps','r');
    xlabel('Time from movement (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 10]);
    title([num2str(sum(csMinuspreRewPiezoTrials)) 'CS- trials with pre-reward movement']);
    subplot(2,2,2);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csMinusfirstPostPiezoEvent,2).*(1000./frameRateHz), (nanstd(csMinusfirstPostPiezoEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','r');
    title('First movement after reward');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csMinusfirstPostRew_piezoAlign,2).*(1000./frameRateHz), (nanstd(csMinusfirstPostRew_piezoAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(nanmean(csMinusfirstPostRew_piezoAlign,1))))),'lineProps','r');
    xlabel('Time from movement (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 10]);
    title([num2str(sum(csMinuspostRewPiezoTrials)) 'CS- trials with post-reward movement']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Movement relative to reward']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, '_csMinuslastVsFirstPiezo.fig'));
    saveas(tempFig, [output_fn  '_csMinuslastVsFirstPiezo.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName '_csMinuslastVsFirstPiezo.fig']));
    saveas(tempFig, [output_fn blockName '_csMinuslastVsFirstPiezo.pdf'],'pdf');
    end
    else
    end
    
    %% FIRST LICK/MOTION AFTER CUE
     if cue==3 || cue==1
    %CS+
     tempFig=setFigure; hold on; % align neural and Piezoing data to the first and last Piezo, respectively
    subplot(2,2,1); hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csPlusfirstPrePiezoEvent,2).*(1000./frameRateHz), (nanstd(csPlusfirstPrePiezoEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    hold on;
    title('First movement after cue');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csPlusfirstPreRew_piezoAlign,2).*(1000./frameRateHz), (nanstd(csPlusfirstPreRew_piezoAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(nanmean(csPlusfirstPreRew_piezoAlign,1))))),'lineProps','k');
    xlabel('Time from movement (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 10]);
    title([num2str(sum(csPluspreFRewPiezoTrials)) 'CS+ trials with post-cue movement']);
    subplot(2,2,2);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csPlusfirstPreLickEvent,2).*(1000./frameRateHz), (nanstd(csPlusfirstPreLickEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    title('First lick after cue');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csPlusFirstPreRew_lickAlign,2), (nanstd(csPlusFirstPreRew_lickAlign,[],2))./sqrt(unique(sum(~isnan(nanmean(csPlusFirstPreRew_lickAlign,1))))),'lineProps','k');
    xlabel('Time from lick (ms)');
    ylabel('Lick Rate (mV)');
    ylim([0 1]);
    title([num2str(sum(csPlusFpreRewTrials)) 'CS+ trials with post-cue licks']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Movement and Licks directly following cue']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, '_csPlusfirstPreRewPiezovsLick.fig'));
    saveas(tempFig, [output_fn  '_csPlusfirstPreRewPiezovsLick.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName '_csPlusFirstPreRewPiezovsLick.fig']));
    saveas(tempFig, [output_fn blockName '_csPlusFirstPreRewPiezovsLick.pdf'],'pdf');
    end
    else
    end
    
    if cue==3 || cue==2    
    %CS-
    tempFig=setFigure; hold on; % align neural and Piezoing data to the first and last Piezo, respectively
    subplot(2,2,1); hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csMinusfirstPrePiezoEvent,2).*(1000./frameRateHz), (nanstd(csMinusfirstPrePiezoEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','r');
    hold on;
    title('First movement after cue');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csMinusfirstPreRew_piezoAlign,2).*(1000./frameRateHz), (nanstd(csMinusfirstPreRew_piezoAlign,[],2).*(1000./frameRateHz))./sqrt(unique(sum(~isnan(nanmean(csMinusfirstPreRew_piezoAlign,1))))),'lineProps','r');
    xlabel('Time from movement (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 10]);
    title([num2str(sum(csMinuspreFRewPiezoTrials)) 'CS- trials with post-cue movement']);
    subplot(2,2,2);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csMinusfirstPreLickEvent,2).*(1000./frameRateHz), (nanstd(csMinusfirstPreLickEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','r');
    title('First licks before cue');
    xlabel('Time from licks (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(csMinusFirstPreRew_lickAlign,2), (nanstd(csMinusFirstPreRew_lickAlign,[],2))./sqrt(unique(sum(~isnan(nanmean(csMinusFirstPreRew_lickAlign,1))))),'lineProps','r');
    xlabel('Time from licks (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title([num2str(sum(csMinusFpreRewTrials)) 'CS- trials with post-cue licks']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Movement and Licks directly following cue']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, '_csMinusFirstPreRewPiezovsLick.fig'));
    saveas(tempFig, [output_fn  '_csMinusFirstPreRewPiezovsLick.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName '_csMinusFirstPreRewPiezovsLick.fig']));
    saveas(tempFig, [output_fn blockName '_csMinusFirstPreRewPiezovsLick.pdf'],'pdf');
    end
    else
    end
    %% another deviation slightly removed from the last deviation
    %Here I'm looking at a small 200ms window around reward delivery for
    %licking events and separating the data into 4 groups (CS+/CS- &
    %with/without licking events)
    
    csPlusNearRewardLick = [csPlusLastPreRew_lickAlign(end-4:end,:); csPlusFirstPostRew_lickAlign(1:25,:)];
    csMinusNearRewardLick = [csMinusLastPreRew_lickAlign(end-4:end,:); csMinusFirstPostRew_lickAlign(1:25,:)];
 
    csPlusPerTrialLickSum = (sum(csPlusNearRewardLick,1,'omitnan'))>0;
    csMinusPerTrialLickSum = (sum(csMinusNearRewardLick,1,'omitnan'))>0;
    
        eventsLickB2=[]; eventsLickRew=[]; ntrials=1; indexB2=[]; indexRew=[];
        
        for mouseIdx = 1:exptCount
        nTrial = (ntrials:animalTrials(1,mouseIdx)*mouseIdx);  
        ind_block2 = find(block2{1,mouseIdx});
        ind_rew = find(~block2{1,mouseIdx});        
        sizeB2=size(nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_block2),2));
        sizeRew=size(nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_rew),2));
        eventsLickB2 = [eventsLickB2 reshape(nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_block2),2),[sizeB2(1) sizeB2(3)])];
        eventsLickRew = [eventsLickRew reshape(nanmean(groupAlign.events{1,mouseIdx}(:,:,ind_rew),2),[sizeRew(1) sizeRew(3)])];
        indexB2=[indexB2 ind_block2];
        indexRew=[indexRew ind_rew];
ntrials=max(nTrial)+1;
        end
        
    eventsWithLickRew = eventsLickRew(:,csPlusPerTrialLickSum==1);
    eventsWithoutLickRew = eventsLickRew(:,csPlusPerTrialLickSum==0);
    eventsWithLickB2 = eventsLickB2(:,csMinusPerTrialLickSum==1);
    eventsWithoutLickB2 = eventsLickB2(:,csMinusPerTrialLickSum==0);
%     
    
    
    tempFig=setFigure; hold on; n=1;
        subplot(2,1,n);hold on;
        shadedErrorBar_CV(tt, nanmean(eventsWithLickRew,2).*(1000./frameRateHz), (nanstd(eventsWithLickRew,[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','k');
        ylim([0 3]);
        vline(767,'g')
        vline(0,'k')       
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title(['CS+ | ' num2str(sum(csPlusPerTrialLickSum==1)) ' trials with lick'])
        n = n+1;
%         if mworks.doBlock2
        subplot(2,1,n);hold on;
        shadedErrorBar_CV(tt, nanmean(eventsWithLickB2,2).*(1000./frameRateHz), (nanstd(eventsWithLickB2,[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','k');
        ylim([0 3]);
        vline(767,'r')
        vline(0,'k')
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title(['CS- | ' num2str(sum(csMinusPerTrialLickSum==1)) ' trials with lick'])
%         end
        sgtitle(['spike rate across animals during trials with lick around reward window [+/- 170ms]']); 
        supertitle([num2str(sum(nsize)) ' neurons from ' num2str(exptCount) ' animals - ' ]);
        if seshType ~= 2
        %savefig(fullfile(output_fn, '_withLickcueAlign_events_Hz.fig'))
        saveas(tempFig, [output_fn  '_withLickcueAlign_events_Hz.pdf'],'pdf');
        else
        %savefig(fullfile(output_fn, [blockName '_withLickcueAlign_events_Hz.fig']));
        saveas(tempFig, [output_fn blockName '_withLickcueAlign_events_Hz.pdf'],'pdf');
        end
        %%%%%%%%%%%%
         tempFig=setFigure; hold on; n=1;
        subplot(2,1,n);hold on;
        shadedErrorBar_CV(tt, nanmean(eventsWithoutLickRew,2).*(1000./frameRateHz), (nanstd(eventsWithoutLickRew,[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','k');
        ylim([0 3]);
        vline(767,'g')
        vline(0,'k')       
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title(['CS+ | ' num2str(sum(csPlusPerTrialLickSum==0)) ' trials without lick'])
        n = n+1;
%         if mworks.doBlock2
        subplot(2,1,n);hold on;
        shadedErrorBar_CV(tt, nanmean(eventsWithoutLickB2,2).*(1000./frameRateHz), (nanstd(eventsWithoutLickB2,[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','k');
        ylim([0 3]);
        vline(767,'r')
        vline(0,'k')
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title(['CS- | ' num2str(sum(csMinusPerTrialLickSum==0)) ' trials without lick'])
%         end
        sgtitle(['spike rate across animals during trials without lick around reward window [+/- 170ms]']);
        supertitle([num2str(sum(nsize)) ' neurons from ' num2str(exptCount) ' animals']);
        if seshType ~= 2
        %savefig(fullfile(output_fn, '_withoutLickcueAlign_events_Hz.fig'))
        saveas(tempFig, [output_fn  '_withoutLickcueAlign_events_Hz.pdf'],'pdf');
        else
        %savefig(fullfile(output_fn, [blockName '_withoutLickcueAlign_events_Hz.fig']));
        saveas(tempFig, [output_fn blockName '_withoutLickcueAlign_events_Hz.pdf'],'pdf');
        end
        
    
    
    %% return to the rest
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
    
   
    tempFig=setFigure; hold on; % still plotting neural data, but seperate the trials based on licking rate of that trial
    subplot(3,1,1);hold on;
    shadedErrorBar_CV(tt,nanmean(lowAllEvents,2).*(1000./frameRateHz), (nanstd(lowAllEvents,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
    hold on;
    shadedErrorBar_CV(tt,nanmean(highAllEvents,2).*(1000./frameRateHz), (nanstd(highAllEvents,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    title(['1s- Rew: ' num2str(chop(nanmean(cell2mat(HL_lickrate.low_rew)),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_lickrate.high_rew)),2)) ' Hz']);
    subplot(3,1,2);hold on;
    shadedErrorBar_CV(tt,nanmean(lowPreEvents,2).*(1000./frameRateHz), (nanstd(lowPreEvents,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
    hold on;
    shadedErrorBar_CV(tt,nanmean(highPreEvents,2).*(1000./frameRateHz), (nanstd(highPreEvents,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    title(['Pre- Rew: ' num2str(chop(nanmean(cell2mat(HL_lickrate.low_prerew)),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_lickrate.high_prerew)),2)) ' Hz']);
    subplot(3,1,3);hold on;
    shadedErrorBar_CV(tt,nanmean(lowPostEvents,2).*(1000./frameRateHz), (nanstd(lowPostEvents,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
    hold on;
    shadedErrorBar_CV(tt,nanmean(highPostEvents,2).*(1000./frameRateHz), (nanstd(highPostEvents,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    title(['Post- Rew: ' num2str(chop(nanmean(cell2mat(HL_lickrate.low_postrew)),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_lickrate.high_postrew)),2)) ' Hz']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | lick bursts by rate: low (blue) & high (black)']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, '_cueAlignSpiking_byLickRate.fig'));
    saveas(tempFig, [output_fn  '_cueAlignSpiking_byLickRate.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byLickRate.fig']));
    saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byLickRate.pdf'],'pdf');
    end
    %%
%Cspk PSTHs for CS+/- trials with low/high lick rate before/after reward%%   

%First, going to remove any 'licks occurring preReward in postRew array and
%vice versa + make a new variable with the lick frequency

preRew_lickAlignHz = sum(preRew_lickAlign(postLick_frames:3*postLick_frames-1,:),1)./length(preRew_lickSearchRange)/frameRateHz;
postRew_lickAlignHz = sum(postRew_lickAlign(postLick_frames+1:3*postLick_frames-1,:),1)./length(postRew_lickSearchRange)/frameRateHz;


ndends=0; ntrials=0; 
low_preRewHz_CSPEvents=[]; high_preRewHz_CSPEvents=[]; low_preRewHz_CSMEvents=[]; high_preRewHz_CSMEvents=[];
low_preRewHz_CSPLicks=[]; high_preRewHz_CSPLicks=[]; low_preRewHz_CSMLicks=[]; high_preRewHz_CSMLicks=[];
low_preRewHz_CSPLickStart = []; low_preRewHz_CSMLickStart = []; high_preRewHz_CSPLickStart = []; high_preRewHz_CSMLickStart = [];

low_postRewHz_CSPEvents=[]; high_postRewHz_CSPEvents=[]; low_postRewHz_CSMEvents=[]; high_postRewHz_CSMEvents=[];
low_postRewHz_CSPLicks=[]; high_postRewHz_CSPLicks=[]; low_postRewHz_CSMLicks=[]; high_postRewHz_CSMLicks=[];
low_postRewHz_CSPLickStart = []; low_postRewHz_CSMLickStart = []; high_postRewHz_CSPLickStart = []; high_postRewHz_CSMLickStart = [];

for mouseIdx = 1:exptCount
    %nDend = (1:nsize(1,mouseIdx))+ndends;
    nTrial = (1:animalTrials(1,mouseIdx))+ntrials; 

    temp_preRewHz_CSP = preRewCSP_lickBurstHz(1,(1:animalTrials(1,mouseIdx))+ntrials);
    temp_preRewHz_CSM = preRewCSM_lickBurstHz(1,(1:animalTrials(1,mouseIdx))+ntrials);
temp_postRewHz_CSP = postRewCSP_lickBurstHz(1,(1:animalTrials(1,mouseIdx))+ntrials);
temp_postRewHz_CSM = postRewCSM_lickBurstHz(1,(1:animalTrials(1,mouseIdx))+ntrials);

    %if lick rate is 0, turn it into a NaN, when I try to align to last lick
    %before/first lick after reward, having 'low lick rate' trials which
    %actually contain no licks is both frustrating and doesn't actually allow
    %me to look at the last lick before reward in instances where licks are
    %occurring infrequently, which is the actual intrigue here.
%temp_preRewHz_CSP(1,temp_preRewHz_CSP==0) = NaN;
%temp_preRewHz_CSM(1,temp_preRewHz_CSM==0) = NaN;
%temp_postRewHz_CSP(1,temp_postRewHz_CSP==0) = NaN;
%temp_postRewHz_CSP(1,temp_postRewHz_CSP==0) = NaN;
    %
    nanIdx_preRewHz_CSP = ~isnan(temp_preRewHz_CSP);
    nanIdx_preRewHz_CSM = ~isnan(temp_preRewHz_CSM);
nanIdx_postRewHz_CSP = ~isnan(temp_postRewHz_CSP);
nanIdx_postRewHz_CSM = ~isnan(temp_postRewHz_CSM);

    [sort_preRewHz_CSP, sort_preRewHz_CSP_ind{1,mouseIdx}] = sort(temp_preRewHz_CSP);
    [sort_preRewHz_CSM, sort_preRewHz_CSM_ind{1,mouseIdx}] = sort(temp_preRewHz_CSM);
[sort_postRewHz_CSP, sort_postRewHz_CSP_ind{1,mouseIdx}] = sort(temp_postRewHz_CSP);
[sort_postRewHz_CSM, sort_postRewHz_CSM_ind{1,mouseIdx}] = sort(temp_postRewHz_CSM);

if seshType ~= 2
       
sort_preRewHz_CSM_ind{1,mouseIdx}(isnan(sort_preRewHz_CSM)) = []; sort_preRewHz_CSM(isnan(sort_preRewHz_CSM)) = [];
sort_postRewHz_CSM_ind{1,mouseIdx}(isnan(sort_postRewHz_CSM)) = []; sort_postRewHz_CSM(isnan(sort_postRewHz_CSM)) = [];

    if size(sort_preRewHz_CSM,2)./4 == ceil(size(sort_preRewHz_CSM,2)./4)
    low_preRewHz_CSM{1,mouseIdx} = sort_preRewHz_CSM(1,1:ceil(size(sort_preRewHz_CSM,2)./4)+1); 
    low_preRewHz_CSM_ind{1,mouseIdx} = sort_preRewHz_CSM_ind{1,mouseIdx}(1,1:ceil(size(sort_preRewHz_CSM,2)./4)+1); 
    else
    low_preRewHz_CSM{1,mouseIdx} = sort_preRewHz_CSM(1,1:ceil(size(sort_preRewHz_CSM,2)./4)); 
    low_preRewHz_CSM_ind{1,mouseIdx} = sort_preRewHz_CSM_ind{1,mouseIdx}(1,1:ceil(size(sort_preRewHz_CSM,2)./4)); 
    end
    high_preRewHz_CSM{1,mouseIdx} = sort_preRewHz_CSM(1,size(sort_preRewHz_CSM,2)-floor(size(sort_preRewHz_CSM,2)./4):size(sort_preRewHz_CSM,2));
    high_preRewHz_CSM_ind{1,mouseIdx} = sort_preRewHz_CSM_ind{1,mouseIdx}(1,size(sort_preRewHz_CSM,2)-floor(size(sort_preRewHz_CSM,2)./4):size(sort_preRewHz_CSM,2));
    high_preRewHz_CSMEvents = [high_preRewHz_CSMEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,sort_preRewHz_CSM_ind{1,mouseIdx}(1,size(sort_preRewHz_CSM_ind{1,mouseIdx},2)-floor(size(sort_preRewHz_CSM_ind{1,mouseIdx},2)./4):size(sort_preRewHz_CSM_ind{1,mouseIdx},2))),3)];
    low_preRewHz_CSMEvents = [low_preRewHz_CSMEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,sort_preRewHz_CSM_ind{1,mouseIdx}(1,1:ceil(size(sort_preRewHz_CSM_ind{1,mouseIdx},2)./4))),3)];  
    high_preRewHz_CSMLicks = [high_preRewHz_CSMLicks groupAlign.licks{1,mouseIdx}(:,sort_preRewHz_CSM_ind{1,mouseIdx}(1,size(sort_preRewHz_CSM_ind{1,mouseIdx},2)-floor(size(sort_preRewHz_CSM_ind{1,mouseIdx},2)./4):size(sort_preRewHz_CSM_ind{1,mouseIdx},2)))];
    low_preRewHz_CSMLicks = [low_preRewHz_CSMLicks groupAlign.licks{1,mouseIdx}(:,sort_preRewHz_CSM_ind{1,mouseIdx}(1,1:ceil(size(sort_preRewHz_CSM_ind{1,mouseIdx},2)./4)))];  
    high_preRewHz_CSMLickStart = [high_preRewHz_CSMLickStart groupAlign.lickStart{1,mouseIdx}(:,sort_preRewHz_CSM_ind{1,mouseIdx}(1,size(sort_preRewHz_CSM_ind{1,mouseIdx},2)-floor(size(sort_preRewHz_CSM_ind{1,mouseIdx},2)./4):size(sort_preRewHz_CSM_ind{1,mouseIdx},2)))];
    low_preRewHz_CSMLickStart = [low_preRewHz_CSMLickStart groupAlign.lickStart{1,mouseIdx}(:,sort_preRewHz_CSM_ind{1,mouseIdx}(1,1:ceil(size(sort_preRewHz_CSM_ind{1,mouseIdx},2)./4)))];
if size(sort_postRewHz_CSM,2)./4 == ceil(size(sort_postRewHz_CSM,2)./4)
low_postRewHz_CSM{1,mouseIdx} = sort_postRewHz_CSM(1,1:ceil(size(sort_postRewHz_CSM,2)./4)+1); 
low_postRewHz_CSM_ind{1,mouseIdx} = sort_postRewHz_CSM_ind{1,mouseIdx}(1,1:ceil(size(sort_postRewHz_CSM,2)./4)+1); 
else
low_postRewHz_CSM{1,mouseIdx} = sort_postRewHz_CSM(1,1:ceil(size(sort_postRewHz_CSM,2)./4)); 
low_postRewHz_CSM_ind{1,mouseIdx} = sort_postRewHz_CSM_ind{1,mouseIdx}(1,1:ceil(size(sort_postRewHz_CSM,2)./4)); 
end
high_postRewHz_CSM{1,mouseIdx} = sort_postRewHz_CSM(1,size(sort_postRewHz_CSM,2)-floor(size(sort_postRewHz_CSM,2)./4):size(sort_postRewHz_CSM,2));
high_postRewHz_CSM_ind{1,mouseIdx} = sort_postRewHz_CSM_ind{1,mouseIdx}(1,size(sort_postRewHz_CSM,2)-floor(size(sort_postRewHz_CSM,2)./4):size(sort_postRewHz_CSM,2));
high_postRewHz_CSMEvents = [high_postRewHz_CSMEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,sort_postRewHz_CSM_ind{1,mouseIdx}(1,size(sort_postRewHz_CSM_ind{1,mouseIdx},2)-floor(size(sort_postRewHz_CSM_ind{1,mouseIdx},2)./4):size(sort_postRewHz_CSM_ind{1,mouseIdx},2))),3)];
low_postRewHz_CSMEvents = [low_postRewHz_CSMEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,sort_postRewHz_CSM_ind{1,mouseIdx}(1,1:ceil(size(sort_postRewHz_CSM_ind{1,mouseIdx},2)./4))),3)];  
high_postRewHz_CSMLicks = [high_postRewHz_CSMLicks groupAlign.licks{1,mouseIdx}(:,sort_postRewHz_CSM_ind{1,mouseIdx}(1,size(sort_postRewHz_CSM_ind{1,mouseIdx},2)-floor(size(sort_postRewHz_CSM_ind{1,mouseIdx},2)./4):size(sort_postRewHz_CSM_ind{1,mouseIdx},2)))];
low_postRewHz_CSMLicks = [low_postRewHz_CSMLicks groupAlign.licks{1,mouseIdx}(:,sort_postRewHz_CSM_ind{1,mouseIdx}(1,1:ceil(size(sort_postRewHz_CSM_ind{1,mouseIdx},2)./4)))];  
high_postRewHz_CSMLickStart = [high_postRewHz_CSMLickStart groupAlign.lickStart{1,mouseIdx}(:,sort_postRewHz_CSM_ind{1,mouseIdx}(1,size(sort_postRewHz_CSM_ind{1,mouseIdx},2)-floor(size(sort_postRewHz_CSM_ind{1,mouseIdx},2)./4):size(sort_postRewHz_CSM_ind{1,mouseIdx},2)))];
low_postRewHz_CSMLickStart = [low_postRewHz_CSMLickStart groupAlign.lickStart{1,mouseIdx}(:,sort_postRewHz_CSM_ind{1,mouseIdx}(1,1:ceil(size(sort_postRewHz_CSM_ind{1,mouseIdx},2)./4)))];

sort_preRewHz_CSP_ind{1,mouseIdx}(isnan(sort_preRewHz_CSP)) = []; sort_preRewHz_CSP(isnan(sort_preRewHz_CSP)) = []; 
sort_postRewHz_CSP_ind{1,mouseIdx}(isnan(sort_postRewHz_CSP)) = []; sort_postRewHz_CSP(isnan(sort_postRewHz_CSP)) = []; 

    if size(sort_preRewHz_CSP,2)./4 == ceil(size(sort_preRewHz_CSP,2)./4)
    low_preRewHz_CSP{1,mouseIdx} = sort_preRewHz_CSP(1,1:ceil(size(sort_preRewHz_CSP,2)./4)+1);
    low_preRewHz_CSP_ind{1,mouseIdx} = sort_preRewHz_CSP_ind{1,mouseIdx}(1,1:ceil(size(sort_preRewHz_CSP,2)./4)+1);
    else
    low_preRewHz_CSP{1,mouseIdx} = sort_preRewHz_CSP(1,1:ceil(size(sort_preRewHz_CSP,2)./4));
    low_preRewHz_CSP_ind{1,mouseIdx} = sort_preRewHz_CSP_ind{1,mouseIdx}(1,1:ceil(size(sort_preRewHz_CSP,2)./4));
    end
    high_preRewHz_CSP{1,mouseIdx} = sort_preRewHz_CSP(1,size(sort_preRewHz_CSP,2)-floor(size(sort_preRewHz_CSP,2)./4):size(sort_preRewHz_CSP,2));
    high_preRewHz_CSP_ind{1,mouseIdx} = sort_preRewHz_CSP_ind{1,mouseIdx}(1,size(sort_preRewHz_CSP,2)-floor(size(sort_preRewHz_CSP,2)./4):size(sort_preRewHz_CSP,2));
    high_preRewHz_CSPEvents = [high_preRewHz_CSPEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,sort_preRewHz_CSP_ind{1,mouseIdx}(1,size(sort_preRewHz_CSP_ind{1,mouseIdx},2)-floor(size(sort_preRewHz_CSP_ind{1,mouseIdx},2)./4):size(sort_preRewHz_CSP_ind{1,mouseIdx},2))),3)];
    low_preRewHz_CSPEvents = [low_preRewHz_CSPEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,sort_preRewHz_CSP_ind{1,mouseIdx}(1,1:ceil(size(sort_preRewHz_CSP_ind{1,mouseIdx},2)./4))),3)];
    high_preRewHz_CSPLicks = [high_preRewHz_CSPLicks groupAlign.licks{1,mouseIdx}(:,sort_preRewHz_CSP_ind{1,mouseIdx}(1,size(sort_preRewHz_CSP_ind{1,mouseIdx},2)-floor(size(sort_preRewHz_CSP_ind{1,mouseIdx},2)./4):size(sort_preRewHz_CSP_ind{1,mouseIdx},2)))];
    low_preRewHz_CSPLicks = [low_preRewHz_CSPLicks groupAlign.licks{1,mouseIdx}(:,sort_preRewHz_CSP_ind{1,mouseIdx}(1,1:ceil(size(sort_preRewHz_CSP_ind{1,mouseIdx},2)./4)))];
    high_preRewHz_CSPLickStart = [high_preRewHz_CSPLickStart groupAlign.lickStart{1,mouseIdx}(:,sort_preRewHz_CSP_ind{1,mouseIdx}(1,size(sort_preRewHz_CSP_ind{1,mouseIdx},2)-floor(size(sort_preRewHz_CSP_ind{1,mouseIdx},2)./4):size(sort_preRewHz_CSP_ind{1,mouseIdx},2)))];
    low_preRewHz_CSPLickStart = [low_preRewHz_CSPLickStart groupAlign.lickStart{1,mouseIdx}(:,sort_preRewHz_CSP_ind{1,mouseIdx}(1,1:ceil(size(sort_preRewHz_CSP_ind{1,mouseIdx},2)./4)))];
if size(sort_postRewHz_CSP,2)./4 == ceil(size(sort_postRewHz_CSP,2)./4)
low_postRewHz_CSP{1,mouseIdx} = sort_postRewHz_CSP(1,1:ceil(size(sort_postRewHz_CSP,2)./4)+1);
low_postRewHz_CSP_ind{1,mouseIdx} = sort_postRewHz_CSP_ind{1,mouseIdx}(1,1:ceil(size(sort_postRewHz_CSP,2)./4)+1);
else
low_postRewHz_CSP{1,mouseIdx} = sort_postRewHz_CSP(1,1:ceil(size(sort_postRewHz_CSP,2)./4));
low_postRewHz_CSP_ind{1,mouseIdx} = sort_postRewHz_CSP_ind{1,mouseIdx}(1,1:ceil(size(sort_postRewHz_CSP,2)./4));
end
high_postRewHz_CSP{1,mouseIdx} = sort_postRewHz_CSP(1,size(sort_postRewHz_CSP,2)-floor(size(sort_postRewHz_CSP,2)./4):size(sort_postRewHz_CSP,2));
high_postRewHz_CSP_ind{1,mouseIdx} = sort_postRewHz_CSP_ind{1,mouseIdx}(1,size(sort_postRewHz_CSP,2)-floor(size(sort_postRewHz_CSP,2)./4):size(sort_postRewHz_CSP,2));
high_postRewHz_CSPEvents = [high_postRewHz_CSPEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,sort_postRewHz_CSP_ind{1,mouseIdx}(1,size(sort_postRewHz_CSP_ind{1,mouseIdx},2)-floor(size(sort_postRewHz_CSP_ind{1,mouseIdx},2)./4):size(sort_postRewHz_CSP_ind{1,mouseIdx},2))),3)];
low_postRewHz_CSPEvents = [low_postRewHz_CSPEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,sort_postRewHz_CSP_ind{1,mouseIdx}(1,1:ceil(size(sort_postRewHz_CSP_ind{1,mouseIdx},2)./4))),3)];
high_postRewHz_CSPLicks = [high_postRewHz_CSPLicks groupAlign.licks{1,mouseIdx}(:,sort_postRewHz_CSP_ind{1,mouseIdx}(1,size(sort_postRewHz_CSP_ind{1,mouseIdx},2)-floor(size(sort_postRewHz_CSP_ind{1,mouseIdx},2)./4):size(sort_postRewHz_CSP_ind{1,mouseIdx},2)))];
low_postRewHz_CSPLicks = [low_postRewHz_CSPLicks groupAlign.licks{1,mouseIdx}(:,sort_postRewHz_CSP_ind{1,mouseIdx}(1,1:ceil(size(sort_postRewHz_CSP_ind{1,mouseIdx},2)./4)))];
high_postRewHz_CSPLickStart = [high_postRewHz_CSPLickStart groupAlign.lickStart{1,mouseIdx}(:,sort_postRewHz_CSP_ind{1,mouseIdx}(1,size(sort_postRewHz_CSP_ind{1,mouseIdx},2)-floor(size(sort_postRewHz_CSP_ind{1,mouseIdx},2)./4):size(sort_postRewHz_CSP_ind{1,mouseIdx},2)))];
low_postRewHz_CSPLickStart = [low_postRewHz_CSPLickStart groupAlign.lickStart{1,mouseIdx}(:,sort_postRewHz_CSP_ind{1,mouseIdx}(1,1:ceil(size(sort_postRewHz_CSP_ind{1,mouseIdx},2)./4)))];

elseif seshType==2    
    if sum(nanIdx_preRewHz_CSM)>0
          low_preRewHz_CSP{1,mouseIdx} = []; low_preRewHz_CSP_ind{1,mouseIdx} = []; low_postRewHz_CSP{1,mouseIdx} = []; low_postRewHz_CSP_ind{1,mouseIdx} = [];
          high_preRewHz_CSP{1,mouseIdx} = []; high_preRewHz_CSP_ind{1,mouseIdx} = []; high_postRewHz_CSP{1,mouseIdx} = []; high_postRewHz_CSP_ind{1,mouseIdx} = [];
sort_preRewHz_CSM_ind{1,mouseIdx}(isnan(sort_preRewHz_CSM)) = []; sort_preRewHz_CSM(isnan(sort_preRewHz_CSM)) = [];
sort_postRewHz_CSM_ind{1,mouseIdx}(isnan(sort_postRewHz_CSM)) = []; sort_postRewHz_CSM(isnan(sort_postRewHz_CSM)) = [];

    if size(sort_preRewHz_CSM,2)./4 == ceil(size(sort_preRewHz_CSM,2)./4)
    low_preRewHz_CSM{1,mouseIdx} = sort_preRewHz_CSM(1,1:ceil(size(sort_preRewHz_CSM,2)./4)+1); 
    low_preRewHz_CSM_ind{1,mouseIdx} = sort_preRewHz_CSM_ind{1,mouseIdx}(1,1:ceil(size(sort_preRewHz_CSM,2)./4)+1); 
    else
    low_preRewHz_CSM{1,mouseIdx} = sort_preRewHz_CSM(1,1:ceil(size(sort_preRewHz_CSM,2)./4)); 
    low_preRewHz_CSM_ind{1,mouseIdx} = sort_preRewHz_CSM_ind{1,mouseIdx}(1,1:ceil(size(sort_preRewHz_CSM,2)./4)); 
    end
    high_preRewHz_CSM{1,mouseIdx} = sort_preRewHz_CSM(1,size(sort_preRewHz_CSM,2)-floor(size(sort_preRewHz_CSM,2)./4):size(sort_preRewHz_CSM,2));
    high_preRewHz_CSM_ind{1,mouseIdx} = sort_preRewHz_CSM_ind{1,mouseIdx}(1,size(sort_preRewHz_CSM,2)-floor(size(sort_preRewHz_CSM,2)./4):size(sort_preRewHz_CSM,2));
    high_preRewHz_CSMEvents = [high_preRewHz_CSMEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,sort_preRewHz_CSM_ind{1,mouseIdx}(1,size(sort_preRewHz_CSM_ind{1,mouseIdx},2)-floor(size(sort_preRewHz_CSM_ind{1,mouseIdx},2)./4):size(sort_preRewHz_CSM_ind{1,mouseIdx},2))),3)];
    low_preRewHz_CSMEvents = [low_preRewHz_CSMEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,sort_preRewHz_CSM_ind{1,mouseIdx}(1,1:ceil(size(sort_preRewHz_CSM_ind{1,mouseIdx},2)./4))),3)];  
    high_preRewHz_CSMLicks = [high_preRewHz_CSMLicks groupAlign.licks{1,mouseIdx}(:,sort_preRewHz_CSM_ind{1,mouseIdx}(1,size(sort_preRewHz_CSM_ind{1,mouseIdx},2)-floor(size(sort_preRewHz_CSM_ind{1,mouseIdx},2)./4):size(sort_preRewHz_CSM_ind{1,mouseIdx},2)))];
    low_preRewHz_CSMLicks = [low_preRewHz_CSMLicks groupAlign.licks{1,mouseIdx}(:,sort_preRewHz_CSM_ind{1,mouseIdx}(1,1:ceil(size(sort_preRewHz_CSM_ind{1,mouseIdx},2)./4)))];  
    high_preRewHz_CSMLickStart = [high_preRewHz_CSMLickStart groupAlign.lickStart{1,mouseIdx}(:,sort_preRewHz_CSM_ind{1,mouseIdx}(1,size(sort_preRewHz_CSM_ind{1,mouseIdx},2)-floor(size(sort_preRewHz_CSM_ind{1,mouseIdx},2)./4):size(sort_preRewHz_CSM_ind{1,mouseIdx},2)))];
    low_preRewHz_CSMLickStart = [low_preRewHz_CSMLickStart groupAlign.lickStart{1,mouseIdx}(:,sort_preRewHz_CSM_ind{1,mouseIdx}(1,1:ceil(size(sort_preRewHz_CSM_ind{1,mouseIdx},2)./4)))];
if size(sort_postRewHz_CSM,2)./4 == ceil(size(sort_postRewHz_CSM,2)./4)
low_postRewHz_CSM{1,mouseIdx} = sort_postRewHz_CSM(1,1:ceil(size(sort_postRewHz_CSM,2)./4)+1); 
low_postRewHz_CSM_ind{1,mouseIdx} = sort_postRewHz_CSM_ind{1,mouseIdx}(1,1:ceil(size(sort_postRewHz_CSM,2)./4)+1); 
else
low_postRewHz_CSM{1,mouseIdx} = sort_postRewHz_CSM(1,1:ceil(size(sort_postRewHz_CSM,2)./4)); 
low_postRewHz_CSM_ind{1,mouseIdx} = sort_postRewHz_CSM_ind{1,mouseIdx}(1,1:ceil(size(sort_postRewHz_CSM,2)./4)); 
end
high_postRewHz_CSM{1,mouseIdx} = sort_postRewHz_CSM(1,size(sort_postRewHz_CSM,2)-floor(size(sort_postRewHz_CSM,2)./4):size(sort_postRewHz_CSM,2));
high_postRewHz_CSM_ind{1,mouseIdx} = sort_postRewHz_CSM_ind{1,mouseIdx}(1,size(sort_postRewHz_CSM,2)-floor(size(sort_postRewHz_CSM,2)./4):size(sort_postRewHz_CSM,2));
high_postRewHz_CSMEvents = [high_postRewHz_CSMEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,sort_postRewHz_CSM_ind{1,mouseIdx}(1,size(sort_postRewHz_CSM_ind{1,mouseIdx},2)-floor(size(sort_postRewHz_CSM_ind{1,mouseIdx},2)./4):size(sort_postRewHz_CSM_ind{1,mouseIdx},2))),3)];
low_postRewHz_CSMEvents = [low_postRewHz_CSMEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,sort_postRewHz_CSM_ind{1,mouseIdx}(1,1:ceil(size(sort_postRewHz_CSM_ind{1,mouseIdx},2)./4))),3)];  
high_postRewHz_CSMLicks = [high_postRewHz_CSMLicks groupAlign.licks{1,mouseIdx}(:,sort_postRewHz_CSM_ind{1,mouseIdx}(1,size(sort_postRewHz_CSM_ind{1,mouseIdx},2)-floor(size(sort_postRewHz_CSM_ind{1,mouseIdx},2)./4):size(sort_postRewHz_CSM_ind{1,mouseIdx},2)))];
low_postRewHz_CSMLicks = [low_postRewHz_CSMLicks groupAlign.licks{1,mouseIdx}(:,sort_postRewHz_CSM_ind{1,mouseIdx}(1,1:ceil(size(sort_postRewHz_CSM_ind{1,mouseIdx},2)./4)))];  
high_postRewHz_CSMLickStart = [high_postRewHz_CSMLickStart groupAlign.lickStart{1,mouseIdx}(:,sort_postRewHz_CSM_ind{1,mouseIdx}(1,size(sort_postRewHz_CSM_ind{1,mouseIdx},2)-floor(size(sort_postRewHz_CSM_ind{1,mouseIdx},2)./4):size(sort_postRewHz_CSM_ind{1,mouseIdx},2)))];
low_postRewHz_CSMLickStart = [low_postRewHz_CSMLickStart groupAlign.lickStart{1,mouseIdx}(:,sort_postRewHz_CSM_ind{1,mouseIdx}(1,1:ceil(size(sort_postRewHz_CSM_ind{1,mouseIdx},2)./4)))];
    
    elseif sum(nanIdx_preRewHz_CSP)>0
        low_preRewHz_CSM{1,mouseIdx} = []; low_preRewHz_CSM_ind{1,mouseIdx} = [];
        low_postRewHz_CSM{1,mouseIdx} = []; low_postRewHz_CSM_ind{1,mouseIdx} = [];
        high_preRewHz_CSM{1,mouseIdx} = []; high_preRewHz_CSM_ind{1,mouseIdx} = [];
        high_postRewHz_CSM{1,mouseIdx} = []; high_postRewHz_CSM_ind{1,mouseIdx} = [];
sort_preRewHz_CSP_ind{1,mouseIdx}(isnan(sort_preRewHz_CSP)) = []; sort_preRewHz_CSP(isnan(sort_preRewHz_CSP)) = []; 
sort_postRewHz_CSP_ind{1,mouseIdx}(isnan(sort_postRewHz_CSP)) = []; sort_postRewHz_CSP(isnan(sort_postRewHz_CSP)) = []; 

    if size(sort_preRewHz_CSP,2)./4 == ceil(size(sort_preRewHz_CSP,2)./4)
    low_preRewHz_CSP{1,mouseIdx} = sort_preRewHz_CSP(1,1:ceil(size(sort_preRewHz_CSP,2)./4)+1);
    low_preRewHz_CSP_ind{1,mouseIdx} = sort_preRewHz_CSP_ind{1,mouseIdx}(1,1:ceil(size(sort_preRewHz_CSP,2)./4)+1);
    else
    low_preRewHz_CSP{1,mouseIdx} = sort_preRewHz_CSP(1,1:ceil(size(sort_preRewHz_CSP,2)./4));
    low_preRewHz_CSP_ind{1,mouseIdx} = sort_preRewHz_CSP_ind{1,mouseIdx}(1,1:ceil(size(sort_preRewHz_CSP,2)./4));
    end
    high_preRewHz_CSP{1,mouseIdx} = sort_preRewHz_CSP(1,size(sort_preRewHz_CSP,2)-floor(size(sort_preRewHz_CSP,2)./4):size(sort_preRewHz_CSP,2));
    high_preRewHz_CSP_ind{1,mouseIdx} = sort_preRewHz_CSP_ind{1,mouseIdx}(1,size(sort_preRewHz_CSP,2)-floor(size(sort_preRewHz_CSP,2)./4):size(sort_preRewHz_CSP,2));
    high_preRewHz_CSPEvents = [high_preRewHz_CSPEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,sort_preRewHz_CSP_ind{1,mouseIdx}(1,size(sort_preRewHz_CSP_ind{1,mouseIdx},2)-floor(size(sort_preRewHz_CSP_ind{1,mouseIdx},2)./4):size(sort_preRewHz_CSP_ind{1,mouseIdx},2))),3)];
    low_preRewHz_CSPEvents = [low_preRewHz_CSPEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,sort_preRewHz_CSP_ind{1,mouseIdx}(1,1:ceil(size(sort_preRewHz_CSP_ind{1,mouseIdx},2)./4))),3)];
    high_preRewHz_CSPLicks = [high_preRewHz_CSPLicks groupAlign.licks{1,mouseIdx}(:,sort_preRewHz_CSP_ind{1,mouseIdx}(1,size(sort_preRewHz_CSP_ind{1,mouseIdx},2)-floor(size(sort_preRewHz_CSP_ind{1,mouseIdx},2)./4):size(sort_preRewHz_CSP_ind{1,mouseIdx},2)))];
    low_preRewHz_CSPLicks = [low_preRewHz_CSPLicks groupAlign.licks{1,mouseIdx}(:,sort_preRewHz_CSP_ind{1,mouseIdx}(1,1:ceil(size(sort_preRewHz_CSP_ind{1,mouseIdx},2)./4)))];
    high_preRewHz_CSPLickStart = [high_preRewHz_CSPLickStart groupAlign.lickStart{1,mouseIdx}(:,sort_preRewHz_CSP_ind{1,mouseIdx}(1,size(sort_preRewHz_CSP_ind{1,mouseIdx},2)-floor(size(sort_preRewHz_CSP_ind{1,mouseIdx},2)./4):size(sort_preRewHz_CSP_ind{1,mouseIdx},2)))];
    low_preRewHz_CSPLickStart = [low_preRewHz_CSPLickStart groupAlign.lickStart{1,mouseIdx}(:,sort_preRewHz_CSP_ind{1,mouseIdx}(1,1:ceil(size(sort_preRewHz_CSP_ind{1,mouseIdx},2)./4)))];
if size(sort_postRewHz_CSP,2)./4 == ceil(size(sort_postRewHz_CSP,2)./4)
low_postRewHz_CSP{1,mouseIdx} = sort_postRewHz_CSP(1,1:ceil(size(sort_postRewHz_CSP,2)./4)+1);
low_postRewHz_CSP_ind{1,mouseIdx} = sort_postRewHz_CSP_ind{1,mouseIdx}(1,1:ceil(size(sort_postRewHz_CSP,2)./4)+1);
else
low_postRewHz_CSP{1,mouseIdx} = sort_postRewHz_CSP(1,1:ceil(size(sort_postRewHz_CSP,2)./4));
low_postRewHz_CSP_ind{1,mouseIdx} = sort_postRewHz_CSP_ind{1,mouseIdx}(1,1:ceil(size(sort_postRewHz_CSP,2)./4));
end
high_postRewHz_CSP{1,mouseIdx} = sort_postRewHz_CSP(1,size(sort_postRewHz_CSP,2)-floor(size(sort_postRewHz_CSP,2)./4):size(sort_postRewHz_CSP,2));
high_postRewHz_CSP_ind{1,mouseIdx} = sort_postRewHz_CSP_ind{1,mouseIdx}(1,size(sort_postRewHz_CSP,2)-floor(size(sort_postRewHz_CSP,2)./4):size(sort_postRewHz_CSP,2));
high_postRewHz_CSPEvents = [high_postRewHz_CSPEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,sort_postRewHz_CSP_ind{1,mouseIdx}(1,size(sort_postRewHz_CSP_ind{1,mouseIdx},2)-floor(size(sort_postRewHz_CSP_ind{1,mouseIdx},2)./4):size(sort_postRewHz_CSP_ind{1,mouseIdx},2))),3)];
low_postRewHz_CSPEvents = [low_postRewHz_CSPEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,sort_postRewHz_CSP_ind{1,mouseIdx}(1,1:ceil(size(sort_postRewHz_CSP_ind{1,mouseIdx},2)./4))),3)];
high_postRewHz_CSPLicks = [high_postRewHz_CSPLicks groupAlign.licks{1,mouseIdx}(:,sort_postRewHz_CSP_ind{1,mouseIdx}(1,size(sort_postRewHz_CSP_ind{1,mouseIdx},2)-floor(size(sort_postRewHz_CSP_ind{1,mouseIdx},2)./4):size(sort_postRewHz_CSP_ind{1,mouseIdx},2)))];
low_postRewHz_CSPLicks = [low_postRewHz_CSPLicks groupAlign.licks{1,mouseIdx}(:,sort_postRewHz_CSP_ind{1,mouseIdx}(1,1:ceil(size(sort_postRewHz_CSP_ind{1,mouseIdx},2)./4)))];
high_postRewHz_CSPLickStart = [high_postRewHz_CSPLickStart groupAlign.lickStart{1,mouseIdx}(:,sort_postRewHz_CSP_ind{1,mouseIdx}(1,size(sort_postRewHz_CSP_ind{1,mouseIdx},2)-floor(size(sort_postRewHz_CSP_ind{1,mouseIdx},2)./4):size(sort_postRewHz_CSP_ind{1,mouseIdx},2)))];
low_postRewHz_CSPLickStart = [low_postRewHz_CSPLickStart groupAlign.lickStart{1,mouseIdx}(:,sort_postRewHz_CSP_ind{1,mouseIdx}(1,1:ceil(size(sort_postRewHz_CSP_ind{1,mouseIdx},2)./4)))];
   end
end

    ndends = max(nDend);
    ntrials = max(nTrial);
end

tt=gtt(1,:); tt=tt(:,21:141);
if seshType==2 && cue==2
else
    tempFig=setFigure; hold on; % still plotting neural data, but seperate the trials based on licking rate of that trial
subplot(2,1,1);    
    shadedErrorBar_CV(tt,mean(low_preRewHz_CSPEvents(21:141,:),2,'omitnan').*(1000./frameRateHz), (std(low_preRewHz_CSPEvents(21:141,:),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
    hold on;
    shadedErrorBar_CV(tt,mean(high_preRewHz_CSPEvents(21:141,:),2,'omitnan').*(1000./frameRateHz), (std(high_preRewHz_CSPEvents(21:141,:),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    errorbar(mean(((high_preRewHz_CSPLickStart-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(low_preRewHz_CSPEvents(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(high_preRewHz_CSPEvents(21:141,:),3),2).*(1000./frameRateHz)])),std(((high_preRewHz_CSPLickStart-prewin_frames).*(1000./frameRateHz)./sqrt(size(high_preRewHz_CSPLickStart,2))),[],2,'omitnan'),'horizontal','k');
    errorbar(mean(((low_preRewHz_CSPLickStart-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(low_preRewHz_CSPEvents(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(high_preRewHz_CSPEvents(21:141,:),3),2).*(1000./frameRateHz)])),std(((low_preRewHz_CSPLickStart-prewin_frames).*(1000./frameRateHz)./sqrt(size(low_preRewHz_CSPLickStart,2))),[],2,'omitnan'),'horizontal','b');    
    xlabel('Time from cue (s)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    title(['CS+ pre-cue: ' num2str(chop(nanmean(cell2mat(low_preRewHz_CSP)),2)) ' Hz (blue ' num2str(length(cell2mat(low_preRewHz_CSP))) ') vs ' num2str(chop(nanmean(cell2mat(high_preRewHz_CSP)),2)) ' Hz (black ' num2str(length(cell2mat(high_preRewHz_CSP))) ')']);
    sgtitle('Cspk rate of neurons from CS+ trials with the lowest and highest lick rate PRE-reward delivery');
subplot(2,1,2);
    shadedErrorBar_CV(tt,mean(low_preRewHz_CSPLicks(21:141,:),2,'omitnan'), (std(low_preRewHz_CSPLicks(21:141,:),[],2,'omitnan'))./sqrt(size(low_preRewHz_CSPLicks,2)),'lineProps','b');
    hold on;
    shadedErrorBar_CV(tt,mean(high_preRewHz_CSPLicks(21:141,:),2,'omitnan'), (std(high_preRewHz_CSPLicks(21:141,:),[],2,'omitnan'))./sqrt(size(high_preRewHz_CSPLicks,2)),'lineProps','k');
    ylim([0 0.5]);
    vline(767,'k');
    xlabel('Time from cue (s)');
    ylabel('Lick rate (Hz)');
    title('Raw lick rate trace for low frequency (blu) trials and high frequency (blk) trials');
    if seshType ~= 2
    %savefig(fullfile(output_fn, 'lowHighLickRate_PreCueHz_CSPlusCueAlign.fig'));
    saveas(tempFig, [output_fn  'lowHighLickRate_PreCueHz_CSPlusCueAlign.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName 'lowHighLickRate_PreCueHz_CSPlusCueAlign.fig']));
    saveas(tempFig, [output_fn blockName 'lowHighLickRate_PreCueHz_CSPlusCueAlign.pdf'],'pdf');
    end
    
    tempFig=setFigure; hold on; % still plotting neural data, but seperate the trials based on licking rate of that trial
subplot(2,1,1);
    shadedErrorBar_CV(tt,nanmean(low_postRewHz_CSPEvents(21:141,:),2).*(1000./frameRateHz), (nanstd(low_postRewHz_CSPEvents(21:141,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
    hold on;
    shadedErrorBar_CV(tt,nanmean(high_postRewHz_CSPEvents(21:141,:),2).*(1000./frameRateHz), (nanstd(high_postRewHz_CSPEvents(21:141,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    hold on;
    errorbar(mean(((high_postRewHz_CSPLickStart-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(low_postRewHz_CSPEvents(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(high_postRewHz_CSPEvents(21:141,:),3),2).*(1000./frameRateHz)])),std(((high_postRewHz_CSPLickStart-prewin_frames).*(1000./frameRateHz)./sqrt(size(high_postRewHz_CSPLickStart,2))),[],2,'omitnan'),'horizontal','k');
    errorbar(mean(((low_postRewHz_CSPLickStart-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(low_postRewHz_CSPEvents(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(high_postRewHz_CSPEvents(21:141,:),3),2).*(1000./frameRateHz)])),std(((low_postRewHz_CSPLickStart-prewin_frames).*(1000./frameRateHz)./sqrt(size(low_postRewHz_CSPLickStart,2))),[],2,'omitnan'),'horizontal','b');    
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    title(['CS+ post-cue: ' num2str(chop(nanmean(cell2mat(low_postRewHz_CSP)),2)) ' Hz (blue ' num2str(length(cell2mat(low_postRewHz_CSP))) ') vs ' num2str(chop(nanmean(cell2mat(high_postRewHz_CSP)),2)) ' Hz (black ' num2str(length(cell2mat(high_postRewHz_CSP))) ')']);
    sgtitle('Cspk rate of neurons from CS+ trials with the lowest and highest lick rate POST-reward delivery');
subplot(2,1,2);
    shadedErrorBar_CV(tt,mean(low_postRewHz_CSPLicks(21:141,:),2,'omitnan'), (std(low_postRewHz_CSPLicks(21:141,:),[],2,'omitnan'))./sqrt(size(low_postRewHz_CSPLicks,2)),'lineProps','b');
    hold on;
    shadedErrorBar_CV(tt,mean(high_postRewHz_CSPLicks(21:141,:),2,'omitnan'), (std(high_postRewHz_CSPLicks(21:141,:),[],2,'omitnan'))./sqrt(size(high_postRewHz_CSPLicks,2)),'lineProps','k');
    ylim([0 0.5]);
    vline(767,'k');
    xlabel('Time from cue (s)');
    ylabel('Lick rate (Hz)');
    title('Raw lick rate trace for low frequency (blu) trials and high frequency (blk) trials');
    if seshType ~= 2
    %savefig(fullfile(output_fn, 'lowHighLickRate_PostCueHz_CSPlusCueAlign.fig'));
    saveas(tempFig, [output_fn  'lowHighLickRate_PostCueHz_CSPlusCueAlign.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName 'lowHighLickRate_PostCueHz_CSPlusCueAlign.fig']));
    saveas(tempFig, [output_fn blockName 'lowHighLickRate_PostCueHz_CSPlusCueAlign.pdf'],'pdf');
    end
end
if seshType==2 && cue==1
else
    tempFig=setFigure; hold on; % still plotting neural data, but seperate the trials based on licking rate of that trial
subplot(2,1,1);
    shadedErrorBar_CV(tt,nanmean(low_preRewHz_CSMEvents(21:141,:),2).*(1000./frameRateHz), (nanstd(low_preRewHz_CSMEvents(21:141,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
    hold on;
    shadedErrorBar_CV(tt,nanmean(high_preRewHz_CSMEvents(21:141,:),2).*(1000./frameRateHz), (nanstd(high_preRewHz_CSMEvents(21:141,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    hold on;
    errorbar(mean(((high_preRewHz_CSMLickStart-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(low_preRewHz_CSMEvents(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(high_preRewHz_CSMEvents(21:141,:),3),2).*(1000./frameRateHz)])),std(((high_preRewHz_CSMLickStart-prewin_frames).*(1000./frameRateHz)./sqrt(size(high_preRewHz_CSMLickStart,2))),[],2,'omitnan'),'horizontal','k');
    errorbar(mean(((low_preRewHz_CSMLickStart-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(low_preRewHz_CSMEvents(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(high_preRewHz_CSMEvents(21:141,:),3),2).*(1000./frameRateHz)])),std(((low_preRewHz_CSMLickStart-prewin_frames).*(1000./frameRateHz)./sqrt(size(low_preRewHz_CSMLickStart,2))),[],2,'omitnan'),'horizontal','b');    
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    title(['CS- pre-cue: ' num2str(chop(nanmean(cell2mat(low_preRewHz_CSM)),2)) ' Hz (blue ' num2str(length(cell2mat(low_preRewHz_CSM))) ') vs ' num2str(chop(nanmean(cell2mat(high_preRewHz_CSM)),2)) ' Hz (black ' num2str(length(cell2mat(high_preRewHz_CSM))) ')']);
    sgtitle('Cspk rate of neurons from CS- trials with the lowest and highest lick rate PRE-reward delivery');
subplot(2,1,2);
    shadedErrorBar_CV(tt,mean(low_preRewHz_CSMLicks(21:141,:),2,'omitnan'), (std(low_preRewHz_CSMLicks(21:141,:),[],2,'omitnan'))./sqrt(size(low_preRewHz_CSMLicks,2)),'lineProps','b');
    hold on;
    shadedErrorBar_CV(tt,mean(high_preRewHz_CSMLicks(21:141,:),2,'omitnan'), (std(high_preRewHz_CSMLicks(21:141,:),[],2,'omitnan'))./sqrt(size(high_preRewHz_CSMLicks,2)),'lineProps','k');
    ylim([0 0.5]);
    vline(767,'k');
    xlabel('Time from cue (s)');
    ylabel('Lick rate (Hz)');
    title('Raw lick rate trace for low frequency (blu) trials and high frequency (blk) trials');
    if seshType ~= 2
    %savefig(fullfile(output_fn, 'lowHighLickRate_PreCueHz_CSMinusCueAlign.fig'));
    saveas(tempFig, [output_fn  'lowHighLickRate_PreCueHz_CSMinusCueAlign.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName 'lowHighLickRate_PreCueHz_CSMinusCueAlign.fig']));
    saveas(tempFig, [output_fn blockName 'lowHighLickRate_PreCueHz_CSMinusCueAlign.pdf'],'pdf');
    end
 
    tempFig=setFigure; hold on; % still plotting neural data, but seperate the trials based on licking rate of that trial
subplot(2,1,1);    
    shadedErrorBar_CV(tt,nanmean(low_postRewHz_CSMEvents(21:141,:),2).*(1000./frameRateHz), (nanstd(low_postRewHz_CSMEvents(21:141,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
    hold on;
    shadedErrorBar_CV(tt,nanmean(high_postRewHz_CSMEvents(21:141,:),2).*(1000./frameRateHz), (nanstd(high_postRewHz_CSMEvents(21:141,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    hold on;
    errorbar(mean(((high_postRewHz_CSMLickStart-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(low_postRewHz_CSMEvents(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(high_postRewHz_CSMEvents(21:141,:),3),2).*(1000./frameRateHz)])),std(((high_postRewHz_CSMLickStart-prewin_frames).*(1000./frameRateHz)./sqrt(size(high_postRewHz_CSMLickStart,2))),[],2,'omitnan'),'horizontal','k');
    errorbar(mean(((low_postRewHz_CSMLickStart-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(low_postRewHz_CSMEvents(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(high_postRewHz_CSMEvents(21:141,:),3),2).*(1000./frameRateHz)])),std(((low_postRewHz_CSMLickStart-prewin_frames).*(1000./frameRateHz)./sqrt(size(low_postRewHz_CSMLickStart,2))),[],2,'omitnan'),'horizontal','b');    
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    title(['CS- post-cue: ' num2str(chop(nanmean(cell2mat(low_postRewHz_CSM)),2)) ' Hz (blue ' num2str(length(cell2mat(low_postRewHz_CSM))) ') vs ' num2str(chop(nanmean(cell2mat(high_postRewHz_CSM)),2)) ' Hz (black ' num2str(length(cell2mat(high_postRewHz_CSM))) ')']);
    sgtitle('Cspk rate of neurons from CS- trials with the lowest and highest lick rate');
subplot(2,1,2);
    shadedErrorBar_CV(tt,mean(low_postRewHz_CSMLicks(21:141,:),2,'omitnan'), (std(low_postRewHz_CSMLicks(21:141,:),[],2,'omitnan'))./sqrt(size(low_postRewHz_CSMLicks,2)),'lineProps','b');
    hold on;
    shadedErrorBar_CV(tt,mean(high_postRewHz_CSMLicks(21:141,:),2,'omitnan'), (std(high_postRewHz_CSMLicks(21:141,:),[],2,'omitnan'))./sqrt(size(high_postRewHz_CSMLicks,2)),'lineProps','k');
    ylim([0 0.5]);
    vline(767,'k');
    xlabel('Time from cue (s)');
    ylabel('Lick rate (Hz)');
    title('Raw lick rate trace for low frequency (blu) trials and high frequency (blk) trials');
    if seshType ~= 2
    %savefig(fullfile(output_fn, 'lowHighLickRate_PostCueHz_CSMinusCueAlign.fig'));
    saveas(tempFig, [output_fn  'lowHighLickRate_PostCueHz_CSMinusCueAlign.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName 'lowHighLickRate_PostCueHz_CSMinusCueAlign.fig']));
    saveas(tempFig, [output_fn blockName 'lowHighLickRate_PostCueHz_CSMinusCueAlign.pdf'],'pdf');
    end
end   
 

%lick frequency-trial block density histogram%%
if seshType~=2
for i = 1:12; trialBlock_ind(i,:) = 2+((i-1)*10):11+((i-1)*10); end
elseif seshType==2
for i = 1:6; trialBlock_ind(i,:) = 2+((i-1)*10):11+((i-1)*10); end
end
for i = 1:size(trialBlock_ind,1)
    if seshType ~= 2 || seshType==2 && cue==1
    for ii = 1:size(animalTrials,2)
    countIn_trialBlock_lowPreCSP(ii,i)= sum(histcounts(low_preRewHz_CSP_ind{1,ii},trialBlock_ind(i,:)));
    countIn_trialBlock_lowPostCSP(ii,i)= sum(histcounts(low_postRewHz_CSP_ind{1,ii},trialBlock_ind(i,:)));
    countIn_trialBlock_highPreCSP(ii,i)= sum(histcounts(high_preRewHz_CSP_ind{1,ii},trialBlock_ind(i,:)));
    countIn_trialBlock_highPostCSP(ii,i)= sum(histcounts(high_postRewHz_CSP_ind{1,ii},trialBlock_ind(i,:)));
    end
    countIn_trialBlock_lowPreCSP(3,i) = mean(countIn_trialBlock_lowPreCSP(:,i),1,'omitnan');
    countIn_trialBlock_lowPostCSP(3,i) = mean(countIn_trialBlock_lowPostCSP(:,i),1,'omitnan');
    countIn_trialBlock_highPreCSP(3,i) = mean(countIn_trialBlock_highPreCSP(:,i),1,'omitnan');
    countIn_trialBlock_highPostCSP(3,i) = mean(countIn_trialBlock_highPostCSP(:,i),1,'omitnan');
    end
    if seshType~=2 || seshType==2 && cue==2
    for ii = 1:size(animalTrials,2)
    countIn_trialBlock_lowPreCSM(ii,i)= sum(histcounts(low_preRewHz_CSM_ind{1,ii},trialBlock_ind(i,:)));
    countIn_trialBlock_lowPostCSM(ii,i)= sum(histcounts(low_postRewHz_CSM_ind{1,ii},trialBlock_ind(i,:)));
    countIn_trialBlock_highPreCSM(ii,i)= sum(histcounts(high_preRewHz_CSM_ind{1,ii},trialBlock_ind(i,:)));
    countIn_trialBlock_highPostCSM(ii,i)= sum(histcounts(high_postRewHz_CSM_ind{1,ii},trialBlock_ind(i,:)));
    end
    countIn_trialBlock_lowPreCSM(3,i) = mean(countIn_trialBlock_lowPreCSM(:,i),1,'omitnan');
    countIn_trialBlock_lowPostCSM(3,i) = mean(countIn_trialBlock_lowPostCSM(:,i),1,'omitnan');
    countIn_trialBlock_highPreCSM(3,i) = mean(countIn_trialBlock_highPreCSM(:,i),1,'omitnan');
    countIn_trialBlock_highPostCSM(3,i) = mean(countIn_trialBlock_highPostCSM(:,i),1,'omitnan');
    end
end
if seshType~=2
countsToPlot = {'countIn_trialBlock_lowPreCSP','countIn_trialBlock_lowPostCSP','countIn_trialBlock_lowPreCSM','countIn_trialBlock_lowPostCSM';'countIn_trialBlock_highPreCSP','countIn_trialBlock_highPostCSP','countIn_trialBlock_highPreCSM','countIn_trialBlock_highPostCSM'};
namesToPlot = {'lowPreCSP','lowPostCSP','lowPreCSM','lowPostCSM';'highPreCSP','highPostCSP','highPreCSM','highPostCSM'};
elseif seshType==2 && cue==1
countsToPlot = {'countIn_trialBlock_lowPreCSP','countIn_trialBlock_lowPostCSP';'countIn_trialBlock_highPreCSP','countIn_trialBlock_highPostCSP'};
namesToPlot = {'lowPreCSP','lowPostCSP';'highPreCSP','highPostCSP'};
elseif seshType==2 && cue==2
countsToPlot = {'countIn_trialBlock_lowPreCSM','countIn_trialBlock_lowPostCSM';'countIn_trialBlock_highPreCSM','countIn_trialBlock_highPostCSM'};
namesToPlot = {'lowPreCSM','lowPostCSM';'highPreCSM','highPostCSM'};
end
for i = 1:size(countsToPlot,2)
    for ii = 1:2; eval(['plot_thisCondition(ii,:) = ' countsToPlot{ii,i} '(3,:);']); end
    tempFig=setFigure; subplot(1,2,1); hold on; 
bar(1:size(trialBlock_ind,1),[plot_thisCondition(1,:); plot_thisCondition(2,:)]);
ylim([0 ceil(max(max(plot_thisCondition)))+1]); ylabel('Count of trials in 10 trial bin'); xlabel('Bin number');
title([namesToPlot{1,i} ' and ' namesToPlot{2,i}]);
hold on;

    if seshType ~= 2
    %savefig(fullfile(output_fn, [countsToPlot{1,i} '.fig']));
    saveas(tempFig, [output_fn countsToPlot{1,i} '.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName countsToPlot{1,i} '.fig']));
    saveas(tempFig, [output_fn blockName countsToPlot{1,i} '.pdf'],'pdf');
    end

end

%making lick aligned now%%
glow_preRewHz_CSP_ind = []; glow_postRewHz_CSP_ind = []; ghigh_preRewHz_CSP_ind = []; ghigh_postRewHz_CSP_ind = []; tempAnimalTrials = [0 animalTrials];
glow_preRewHz_CSM_ind = []; glow_postRewHz_CSM_ind = []; ghigh_preRewHz_CSM_ind = []; ghigh_postRewHz_CSM_ind = [];
for i = 1:size(animalTrials,2)
    if seshType~=2 || seshType==2 && cue==1
glow_preRewHz_CSP_ind = [glow_preRewHz_CSP_ind low_preRewHz_CSP_ind{1,i}+sum(tempAnimalTrials(1,1:i))];
glow_postRewHz_CSP_ind = [glow_postRewHz_CSP_ind low_postRewHz_CSP_ind{1,i}+sum(tempAnimalTrials(1,1:i))];
ghigh_preRewHz_CSP_ind = [ghigh_preRewHz_CSP_ind high_preRewHz_CSP_ind{1,i}+sum(tempAnimalTrials(1,1:i))];
ghigh_postRewHz_CSP_ind = [ghigh_postRewHz_CSP_ind high_postRewHz_CSP_ind{1,i}+sum(tempAnimalTrials(1,1:i))];
    elseif seshType~=2 || seshType==2 && cue==2
glow_preRewHz_CSM_ind = [glow_preRewHz_CSM_ind low_preRewHz_CSM_ind{1,i}+sum(tempAnimalTrials(1,1:i))];
glow_postRewHz_CSM_ind = [glow_postRewHz_CSM_ind low_postRewHz_CSM_ind{1,i}+sum(tempAnimalTrials(1,1:i))];
ghigh_preRewHz_CSM_ind = [ghigh_preRewHz_CSM_ind high_preRewHz_CSM_ind{1,i}+sum(tempAnimalTrials(1,1:i))];
ghigh_postRewHz_CSM_ind = [ghigh_postRewHz_CSM_ind high_postRewHz_CSM_ind{1,i}+sum(tempAnimalTrials(1,1:i))];
    end
end
low_preRewHz_lickAlignEvents_CSP = mean(lastPreRew_lickAlignEvents(:,:,glow_preRewHz_CSP_ind),3,'omitnan'); low_preRew_lickAlign_CSP = lastPreRew_lickAlign(:,glow_preRewHz_CSP_ind);
high_preRewHz_lickAlignEvents_CSP = mean(lastPreRew_lickAlignEvents(:,:,ghigh_preRewHz_CSP_ind),3,'omitnan'); high_preRew_lickAlign_CSP = lastPreRew_lickAlign(:,ghigh_preRewHz_CSP_ind);
low_postRewHz_lickAlignEvents_CSP = mean(firstPostRew_lickAlignEvents(:,:,glow_postRewHz_CSP_ind),3,'omitnan'); low_postRew_lickAlign_CSP = firstPostRew_lickAlign(:,glow_postRewHz_CSP_ind);
high_postRewHz_lickAlignEvents_CSP = mean(firstPostRew_lickAlignEvents(:,:,ghigh_postRewHz_CSP_ind),3,'omitnan'); high_postRew_lickAlign_CSP = firstPostRew_lickAlign(:,ghigh_postRewHz_CSP_ind);

low_preRewHz_lickAlignEvents_CSM = mean(lastPreRew_lickAlignEvents(:,:,glow_preRewHz_CSM_ind),3,'omitnan'); low_preRew_lickAlign_CSM = lastPreRew_lickAlign(:,glow_preRewHz_CSM_ind);
high_preRewHz_lickAlignEvents_CSM = mean(lastPreRew_lickAlignEvents(:,:,ghigh_preRewHz_CSM_ind),3,'omitnan'); high_preRew_lickAlign_CSM = lastPreRew_lickAlign(:,ghigh_preRewHz_CSM_ind);
low_postRewHz_lickAlignEvents_CSM = mean(firstPostRew_lickAlignEvents(:,:,glow_postRewHz_CSM_ind),3,'omitnan'); low_postRew_lickAlign_CSM = firstPostRew_lickAlign(:,glow_postRewHz_CSM_ind);
high_postRewHz_lickAlignEvents_CSM = mean(firstPostRew_lickAlignEvents(:,:,ghigh_postRewHz_CSM_ind),3,'omitnan'); high_postRew_lickAlign_CSM = firstPostRew_lickAlign(:,ghigh_postRewHz_CSM_ind);



%low|high preRew CSP
%low|high postRew CSP
    tl_rew = (1-postLick_frames:(postLick_frames*2)).*(1000./frameRateHz); %orig :: -500ms to 1000 ms around reward
    tempFig=setFigure; hold on; % align neural and Piezoing data to the first and last Piezo, respectively
    subplot(2,2,1); hold on;
    shadedErrorBar_CV(tl_rew, mean(low_preRewHz_lickAlignEvents_CSP,2,'omitnan').*(1000./frameRateHz), (std(low_preRewHz_lickAlignEvents_CSP,[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
    shadedErrorBar_CV(tl_rew, mean(high_preRewHz_lickAlignEvents_CSP,2,'omitnan').*(1000./frameRateHz), (std(high_preRewHz_lickAlignEvents_CSP,[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    hold on;
    title('Last lick before reward');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    shadedErrorBar_CV(tl_rew, mean(low_preRew_lickAlign_CSP,2,'omitnan'), (std(low_preRew_lickAlign_CSP,[],2,'omitnan'))./sqrt(size(low_preRew_lickAlign_CSP,2)),'lineProps','b');
    shadedErrorBar_CV(tl_rew, mean(high_preRew_lickAlign_CSP,2,'omitnan'), (std(high_preRew_lickAlign_CSP,[],2,'omitnan'))./sqrt(size(high_preRew_lickAlign_CSP,2)),'lineProps','k');
    xlabel('Time from movement (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title(['CS+ trials with low(blu-' num2str(sum(~isnan(mean(low_preRew_lickAlign_CSP,1)))) ')|high(blk-' num2str(sum(~isnan(mean(high_preRew_lickAlign_CSP,1)))) ') pre-reward lick rate']);
    subplot(2,2,2);hold on;
    shadedErrorBar_CV(tl_rew, mean(low_postRewHz_lickAlignEvents_CSP,2,'omitnan').*(1000./frameRateHz), (std(low_postRewHz_lickAlignEvents_CSP,[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
    shadedErrorBar_CV(tl_rew, mean(high_postRewHz_lickAlignEvents_CSP,2,'omitnan').*(1000./frameRateHz), (std(high_postRewHz_lickAlignEvents_CSP,[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    title('First lick after reward');
    xlabel('Time from  (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    shadedErrorBar_CV(tl_rew, mean(low_postRew_lickAlign_CSP,2,'omitnan'), (std(low_postRew_lickAlign_CSP,[],2,'omitnan'))./sqrt(size(low_postRew_lickAlign_CSP,2)),'lineProps','b');
    shadedErrorBar_CV(tl_rew, mean(high_postRew_lickAlign_CSP,2,'omitnan'), (std(high_postRew_lickAlign_CSP,[],2,'omitnan'))./sqrt(size(high_postRew_lickAlign_CSP,2)),'lineProps','k');
    xlabel('Time from movement (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title(['CS+ trials with low(blu-' num2str(sum(~isnan(mean(low_postRew_lickAlign_CSP,1)))) ')|high(blk-' num2str(sum(~isnan(mean(high_postRew_lickAlign_CSP,1)))) ') post-reward lick rate']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Licks relative to reward']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, '_lowHigh_preRew_postRew_CS+.fig'));
    saveas(tempFig, [output_fn  '_lowHigh_preRew_postRew_CS+.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName '_lowHigh_preRew_postRew_CS+.fig']));
    saveas(tempFig, [output_fn blockName '_lowHigh_preRew_postRew_CS+.pdf'],'pdf');
    end
%low|high preRew CSM
%low|high postRew CSM
    tl_rew = (1-postLick_frames:(postLick_frames*2)).*(1000./frameRateHz); %orig :: -500ms to 1000 ms around reward
    tempFig=setFigure; hold on; % align neural and Piezoing data to the first and last Piezo, respectively
    subplot(2,2,1); hold on;
    shadedErrorBar_CV(tl_rew, mean(low_preRewHz_lickAlignEvents_CSM,2,'omitnan').*(1000./frameRateHz), (std(low_preRewHz_lickAlignEvents_CSM,[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
    shadedErrorBar_CV(tl_rew, mean(high_preRewHz_lickAlignEvents_CSM,2,'omitnan').*(1000./frameRateHz), (std(high_preRewHz_lickAlignEvents_CSM,[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    hold on;
    title('Last lick before reward');
    xlabel('Time from movement (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    shadedErrorBar_CV(tl_rew, mean(low_preRew_lickAlign_CSM,2,'omitnan'), (std(low_preRew_lickAlign_CSM,[],2,'omitnan'))./sqrt(size(low_preRew_lickAlign_CSM,2)),'lineProps','b');
    shadedErrorBar_CV(tl_rew, mean(high_preRew_lickAlign_CSM,2,'omitnan'), (std(high_preRew_lickAlign_CSM,[],2,'omitnan'))./sqrt(size(high_preRew_lickAlign_CSM,2)),'lineProps','k');
    xlabel('Time from movement (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title(['CS- trials with low(blu-' num2str(sum(~isnan(mean(low_preRew_lickAlign_CSM,1)))) ')|high(blk-' num2str(sum(~isnan(mean(high_preRew_lickAlign_CSM,1)))) ') pre-reward lick rate']);
    subplot(2,2,2);hold on;
    shadedErrorBar_CV(tl_rew, mean(low_postRewHz_lickAlignEvents_CSM,2,'omitnan').*(1000./frameRateHz), (std(low_postRewHz_lickAlignEvents_CSM,[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
    shadedErrorBar_CV(tl_rew, mean(high_postRewHz_lickAlignEvents_CSM,2,'omitnan').*(1000./frameRateHz), (std(high_postRewHz_lickAlignEvents_CSM,[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    title('First lick after reward');
    xlabel('Time from  (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    shadedErrorBar_CV(tl_rew, mean(low_postRew_lickAlign_CSM,2,'omitnan'), (std(low_postRew_lickAlign_CSM,[],2,'omitnan'))./sqrt(size(low_postRew_lickAlign_CSM,2)),'lineProps','b');
    shadedErrorBar_CV(tl_rew, mean(high_postRew_lickAlign_CSM,2,'omitnan'), (std(high_postRew_lickAlign_CSM,[],2,'omitnan'))./sqrt(size(high_postRew_lickAlign_CSM,2)),'lineProps','k');
    xlabel('Time from movement (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title(['CS- trials with low(blu-' num2str(sum(~isnan(mean(low_postRew_lickAlign_CSM,1)))) ')|high(blk-' num2str(sum(~isnan(mean(high_postRew_lickAlign_CSM,1)))) ') post-reward lick rate']);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Licks relative to reward']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, '_lowHigh_preRew_postRew_CS-.fig'));
    saveas(tempFig, [output_fn  '_lowHigh_preRew_postRew_CS-.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName '_lowHigh_preRew_postRew_CS-.fig']));
    saveas(tempFig, [output_fn blockName '_lowHigh_preRew_postRew_CS-.pdf'],'pdf');
    end
%%end - for now%
%% %it is now later%
%need to look at high vs. low lick rate trials across segments of a session, so 1,2,3rd phase of a session%%
for ii = 1:size(animalTrials,2)
for i = 1:size(repTrials,1)
    preRHz_CSM_indInd{ii,i} = sort_preRewHz_CSM_ind{1,ii}(ismember(sort_preRewHz_CSM_ind{1,ii},repTrials(i,:))); 
    postRHz_CSM_indInd{ii,i} = sort_postRewHz_CSM_ind{1,ii}(ismember(sort_postRewHz_CSM_ind{1,ii},repTrials(i,:)));
    preRHz_CSP_indInd{ii,i} = sort_preRewHz_CSP_ind{1,ii}(ismember(sort_preRewHz_CSP_ind{1,ii},repTrials(i,:)));
    postRHz_CSP_indInd{ii,i} = sort_postRewHz_CSP_ind{1,ii}(ismember(sort_postRewHz_CSP_ind{1,ii},repTrials(i,:)));
end
end
%%quick indexing for data organized across mice%%
tempAnimalTrials = [0 animalTrials];
for i = 1:size(repTrials,1)
for ii = 1:size(animalTrials,2)
glowRep_preRewHz_CSP_ind{ii,i}= []; glowRep_postRewHz_CSP_ind{ii,i} = []; ghighRep_preRewHz_CSP_ind{ii,i} = []; ghighRep_postRewHz_CSP_ind{ii,i} = []; 
glowRep_preRewHz_CSM_ind{ii,i} = []; glowRep_postRewHz_CSM_ind{ii,i} = []; ghighRep_preRewHz_CSM_ind{ii,i} = []; ghighRep_postRewHz_CSM_ind{ii,i} = [];

glowRep_preRewHz_CSP_ind{ii,i} = [glowRep_preRewHz_CSP_ind{ii,i} low_preRewHz_CSP_ind{1,ii}+sum(tempAnimalTrials(1,1:ii))];
glowRep_postRewHz_CSP_ind{ii,i} = [glowRep_postRewHz_CSP_ind{ii,i} low_postRewHz_CSP_ind{1,ii}+sum(tempAnimalTrials(1,1:ii))];
ghighRep_preRewHz_CSP_ind{ii,i} = [ghighRep_preRewHz_CSP_ind{ii,i} high_preRewHz_CSP_ind{1,ii}+sum(tempAnimalTrials(1,1:ii))];
ghighRep_postRewHz_CSP_ind{ii,i} = [ghighRep_postRewHz_CSP_ind{ii,i} high_postRewHz_CSP_ind{1,ii}+sum(tempAnimalTrials(1,1:ii))];
glowRep_preRewHz_CSM_ind{ii,i} = [glowRep_preRewHz_CSM_ind{ii,i} low_preRewHz_CSM_ind{1,ii}+sum(tempAnimalTrials(1,1:ii))];
glowRep_postRewHz_CSM_ind{ii,i} = [glowRep_postRewHz_CSM_ind{ii,i} low_postRewHz_CSM_ind{1,ii}+sum(tempAnimalTrials(1,1:ii))];
ghighRep_preRewHz_CSM_ind{ii,i} = [ghighRep_preRewHz_CSM_ind{ii,i} high_preRewHz_CSM_ind{1,ii}+sum(tempAnimalTrials(1,1:ii))];
ghighRep_postRewHz_CSM_ind{ii,i} = [ghighRep_postRewHz_CSM_ind{ii,i} high_postRewHz_CSM_ind{1,ii}+sum(tempAnimalTrials(1,1:ii))];
end    
end
%%end%%
for i = 1:size(repTrials,1)
low_preRHz_CSM_licks{1,i} = [];low_postRHz_CSM_licks{1,i} = [];low_preRHz_CSP_licks{1,i} = [];low_postRHz_CSP_licks{1,i} = [];high_preRHz_CSM_licks{1,i} = [];high_postRHz_CSM_licks{1,i} = [];high_preRHz_CSP_licks{1,i} = [];high_postRHz_CSP_licks{1,i} =[];
low_preRHz_CSM_lickStart{1,i} = [];low_postRHz_CSM_lickStart{1,i} = [];low_preRHz_CSP_lickStart{1,i} = [];low_postRHz_CSP_lickStart{1,i} = [];high_preRHz_CSM_lickStart{1,i} = [];high_postRHz_CSM_lickStart{1,i} = [];high_preRHz_CSP_lickStart{1,i} = [];high_postRHz_CSP_lickStart{1,i} =[];
low_preRHz_CSM_lickAlignEvents{1,i} = [];low_postRHz_CSM_lickAlignEvents{1,i} = [];low_preRHz_CSP_lickAlignEvents{1,i} = [];low_postRHz_CSP_lickAlignEvents{1,i} = [];  high_preRHz_CSM_lickAlignEvents{1,i} = [];high_postRHz_CSM_lickAlignEvents{1,i} = [];high_preRHz_CSP_lickAlignEvents{1,i} = [];high_postRHz_CSP_lickAlignEvents{1,i} =[];
low_preRHz_CSM_lickAlign{1,i} = [];low_postRHz_CSM_lickAlign{1,i} = [];low_preRHz_CSP_lickAlign{1,i} = [];low_postRHz_CSP_lickAlign{1,i} = [];high_preRHz_CSM_lickAlign{1,i} = [];high_postRHz_CSM_lickAlign{1,i} = [];high_preRHz_CSP_lickAlign{1,i} = [];high_postRHz_CSP_lickAlign{1,i} =[];
low_preRHz_CSM_events{1,i} = []; low_postRHz_CSM_events{1,i} = []; low_preRHz_CSP_events{1,i} = []; low_postRHz_CSP_events{1,i} = []; high_preRHz_CSM_events{1,i} = []; high_postRHz_CSM_events{1,i} = []; high_preRHz_CSP_events{1,i} = [];high_postRHz_CSP_events{1,i} = [];
for ii = 1:size(animalTrials,2)
    % cue align    
    % lick align    
if cue~=1
    if ~isempty(preRHz_CSM_indInd{ii,i})
    low_preRHz_CSM_events{1,i} = [low_preRHz_CSM_events{1,i} mean(groupAlign.events{1,ii}(:,:,preRHz_CSM_indInd{ii,i}(1,1:ceil(size(preRHz_CSM_indInd{ii,i},2)./4))),3,'omitnan')];
    low_preRHz_CSM_licks{1,i} = [low_preRHz_CSM_licks{1,i} groupAlign.licks{1,ii}(:,preRHz_CSM_indInd{ii,i}(1,1:ceil(size(preRHz_CSM_indInd{ii,i},2)./4)))];
    low_preRHz_CSM_lickStart{1,i} = [low_preRHz_CSM_lickStart{1,i} groupAlign.lickStart{1,ii}(:,preRHz_CSM_indInd{ii,i}(1,1:ceil(size(preRHz_CSM_indInd{ii,i},2)./4)))];
    low_preRHz_CSM_lickAlignEvents{1,i} = [low_preRHz_CSM_lickAlignEvents{1,i} mean(lastPreRew_lickAlignEvents(:,:,glowRep_preRewHz_CSM_ind{ii,i}),3,'omitnan')];
    low_preRHz_CSM_lickAlign{1,i} = [low_preRHz_CSM_lickAlign{1,i} lastPreRew_lickAlign(:,glowRep_preRewHz_CSM_ind{ii,i})];
    
    high_preRHz_CSM_events{1,i} = [high_preRHz_CSM_events{1,i} mean(groupAlign.events{1,ii}(:,:,preRHz_CSM_indInd{ii,i}(1,size(preRHz_CSM_indInd{ii,i},2)-floor(size(preRHz_CSM_indInd{ii,i},2)./4):size(preRHz_CSM_indInd{ii,i},2))),3,'omitnan')];
    high_preRHz_CSM_licks{1,i} = [high_preRHz_CSM_licks{1,i} groupAlign.licks{1,ii}(:,preRHz_CSM_indInd{ii,i}(1,size(preRHz_CSM_indInd{ii,i},2)-floor(size(preRHz_CSM_indInd{ii,i},2)./4):size(preRHz_CSM_indInd{ii,i},2)))];
    high_preRHz_CSM_lickStart{1,i} = [high_preRHz_CSM_lickStart{1,i} groupAlign.lickStart{1,ii}(:,preRHz_CSM_indInd{ii,i}(1,size(preRHz_CSM_indInd{ii,i},2)-floor(size(preRHz_CSM_indInd{ii,i},2)./4):size(preRHz_CSM_indInd{ii,i},2)))];
    high_preRHz_CSM_lickAlignEvents{1,i} = [high_preRHz_CSM_lickAlignEvents{1,i} mean(lastPreRew_lickAlignEvents(:,:,ghighRep_preRewHz_CSM_ind{ii,i}),3,'omitnan')];
    high_preRHz_CSM_lickAlign{1,i} =[high_preRHz_CSM_lickAlign{1,i}  lastPreRew_lickAlign(:,ghighRep_preRewHz_CSM_ind{ii,i})];
    end
    if ~isempty(postRHz_CSM_indInd{ii,i})
    low_postRHz_CSM_events{1,i} = [low_postRHz_CSM_events{1,i} mean(groupAlign.events{1,ii}(:,:,postRHz_CSM_indInd{ii,i}(1,1:ceil(size(postRHz_CSM_indInd{ii,i},2)./4))),3,'omitnan')];
    low_postRHz_CSM_licks{1,i} = [low_postRHz_CSM_licks{1,i} groupAlign.licks{1,ii}(:,postRHz_CSM_indInd{ii,i}(1,1:ceil(size(postRHz_CSM_indInd{ii,i},2)./4)))];
    low_postRHz_CSM_lickStart{1,i} = [low_postRHz_CSM_lickStart{1,i} groupAlign.lickStart{1,ii}(:,postRHz_CSM_indInd{ii,i}(1,1:ceil(size(postRHz_CSM_indInd{ii,i},2)./4)))];
    low_postRHz_CSM_lickAlignEvents{1,i} = [low_postRHz_CSM_lickAlignEvents{1,i} mean(firstPostRew_lickAlignEvents(:,:,glowRep_postRewHz_CSM_ind{ii,i}),3,'omitnan')];
    low_postRHz_CSM_lickAlign{1,i} = [low_postRHz_CSM_lickAlign{1,i} firstPostRew_lickAlign(:,glowRep_postRewHz_CSM_ind{ii,i})];
        
    high_postRHz_CSM_events{1,i} = [high_postRHz_CSM_events{1,i} mean(groupAlign.events{1,ii}(:,:,postRHz_CSM_indInd{ii,i}(1,size(postRHz_CSM_indInd{ii,i},2)-floor(size(postRHz_CSM_indInd{ii,i},2)./4):size(postRHz_CSM_indInd{ii,i},2))),3,'omitnan')];
    high_postRHz_CSM_licks{1,i} = [high_postRHz_CSM_licks{1,i} groupAlign.licks{1,ii}(:,postRHz_CSM_indInd{ii,i}(1,size(postRHz_CSM_indInd{ii,i},2)-floor(size(postRHz_CSM_indInd{ii,i},2)./4):size(postRHz_CSM_indInd{ii,i},2)))];
    high_postRHz_CSM_lickStart{1,i} = [high_postRHz_CSM_lickStart{1,i} groupAlign.lickStart{1,ii}(:,postRHz_CSM_indInd{ii,i}(1,size(postRHz_CSM_indInd{ii,i},2)-floor(size(postRHz_CSM_indInd{ii,i},2)./4):size(postRHz_CSM_indInd{ii,i},2)))];
    high_postRHz_CSM_lickAlignEvents{1,i} = [high_postRHz_CSM_lickAlignEvents{1,i} mean(firstPostRew_lickAlignEvents(:,:,ghighRep_postRewHz_CSM_ind{ii,i}),3,'omitnan')];
    high_postRHz_CSM_lickAlign{1,i} = [high_postRHz_CSM_lickAlign{1,i} firstPostRew_lickAlign(:,ghighRep_postRewHz_CSM_ind{ii,i})];
    end
end
if cue~=2
    if ~isempty(preRHz_CSP_indInd{ii,i})
    low_preRHz_CSP_events{1,i} = [low_preRHz_CSP_events{1,i} mean(groupAlign.events{1,ii}(:,:,preRHz_CSP_indInd{ii,i}(1,1:ceil(size(preRHz_CSP_indInd{ii,i},2)./4))),3,'omitnan')];
    low_preRHz_CSP_licks{1,i} = [low_preRHz_CSP_licks{1,i} groupAlign.licks{1,ii}(:,preRHz_CSP_indInd{ii,i}(1,1:ceil(size(preRHz_CSP_indInd{ii,i},2)./4)))];
    low_preRHz_CSP_lickStart{1,i} = [low_preRHz_CSP_lickStart{1,i} groupAlign.lickStart{1,ii}(:,preRHz_CSP_indInd{ii,i}(1,1:ceil(size(preRHz_CSP_indInd{ii,i},2)./4)))];
    low_preRHz_CSP_lickAlignEvents{1,i} = [low_preRHz_CSP_lickAlignEvents{1,i} mean(lastPreRew_lickAlignEvents(:,:,glowRep_preRewHz_CSP_ind{ii,i}),3,'omitnan')];
    low_preRHz_CSP_lickAlign{1,i} = [low_preRHz_CSP_lickAlign{1,i} lastPreRew_lickAlign(:,glowRep_preRewHz_CSP_ind{ii,i})];
    
    high_preRHz_CSP_events{1,i} = [high_preRHz_CSP_events{1,i} mean(groupAlign.events{1,ii}(:,:,preRHz_CSP_indInd{ii,i}(1,size(preRHz_CSP_indInd{ii,i},2)-floor(size(preRHz_CSP_indInd{ii,i},2)./4):size(preRHz_CSP_indInd{ii,i},2))),3,'omitnan')];
    high_preRHz_CSP_licks{1,i} = [high_preRHz_CSP_licks{1,i} groupAlign.licks{1,ii}(:,preRHz_CSP_indInd{ii,i}(1,size(preRHz_CSP_indInd{ii,i},2)-floor(size(preRHz_CSP_indInd{ii,i},2)./4):size(preRHz_CSP_indInd{ii,i},2)))];
    high_preRHz_CSP_lickStart{1,i} = [high_preRHz_CSP_lickStart{1,i} groupAlign.lickStart{1,ii}(:,preRHz_CSP_indInd{ii,i}(1,size(preRHz_CSP_indInd{ii,i},2)-floor(size(preRHz_CSP_indInd{ii,i},2)./4):size(preRHz_CSP_indInd{ii,i},2)))];
    high_preRHz_CSP_lickAlignEvents{1,i} = [high_preRHz_CSP_lickAlignEvents{1,i} mean(lastPreRew_lickAlignEvents(:,:,ghighRep_preRewHz_CSP_ind{ii,i}),3,'omitnan')];
    high_preRHz_CSP_lickAlign{1,i} = [high_preRHz_CSP_lickAlign{1,i} lastPreRew_lickAlign(:,ghighRep_preRewHz_CSP_ind{ii,i})];
    end
    if ~isempty(postRHz_CSP_indInd{ii,i})
    low_postRHz_CSP_events{1,i} = [low_postRHz_CSP_events{1,i} mean(groupAlign.events{1,ii}(:,:,postRHz_CSP_indInd{ii,i}(1,1:ceil(size(postRHz_CSP_indInd{ii,i},2)./4))),3,'omitnan')];
    low_postRHz_CSP_licks{1,i} = [low_postRHz_CSP_licks{1,i} groupAlign.licks{1,ii}(:,postRHz_CSP_indInd{ii,i}(1,1:ceil(size(postRHz_CSP_indInd{ii,i},2)./4)))];
    low_postRHz_CSP_lickStart{1,i} = [low_postRHz_CSP_lickStart{1,i} groupAlign.lickStart{1,ii}(:,postRHz_CSP_indInd{ii,i}(1,1:ceil(size(postRHz_CSP_indInd{ii,i},2)./4)))];
    low_postRHz_CSP_lickAlignEvents{1,i} = [low_postRHz_CSP_lickAlignEvents{1,i} mean(firstPostRew_lickAlignEvents(:,:,glowRep_postRewHz_CSP_ind{ii,i}),3,'omitnan')];
    low_postRHz_CSP_lickAlign{1,i} = [low_postRHz_CSP_lickAlign{1,i} firstPostRew_lickAlign(:,glowRep_postRewHz_CSP_ind{ii,i})];

    high_postRHz_CSP_events{1,i} = [high_postRHz_CSP_events{1,i} mean(groupAlign.events{1,ii}(:,:,postRHz_CSP_indInd{ii,i}(1,size(postRHz_CSP_indInd{ii,i},2)-floor(size(postRHz_CSP_indInd{ii,i},2)./4):size(postRHz_CSP_indInd{ii,i},2))),3,'omitnan')];
    high_postRHz_CSP_licks{1,i} = [high_postRHz_CSP_licks{1,i} groupAlign.licks{1,ii}(:,postRHz_CSP_indInd{ii,i}(1,size(postRHz_CSP_indInd{ii,i},2)-floor(size(postRHz_CSP_indInd{ii,i},2)./4):size(postRHz_CSP_indInd{ii,i},2)))];
    high_postRHz_CSP_lickStart{1,i} = [high_postRHz_CSP_lickStart{1,i} groupAlign.lickStart{1,ii}(:,postRHz_CSP_indInd{ii,i}(1,size(postRHz_CSP_indInd{ii,i},2)-floor(size(postRHz_CSP_indInd{ii,i},2)./4):size(postRHz_CSP_indInd{ii,i},2)))];
    high_postRHz_CSP_lickAlignEvents{1,i} = [high_postRHz_CSP_lickAlignEvents{1,i} mean(firstPostRew_lickAlignEvents(:,:,ghighRep_postRewHz_CSP_ind{ii,i}),3,'omitnan')];
    high_postRHz_CSP_lickAlign{1,i} = [high_postRHz_CSP_lickAlign{1,i} firstPostRew_lickAlign(:,ghighRep_postRewHz_CSP_ind{ii,i})];
    end
end
end
end
% colormap RGB

colorsRepCSP={[0 0 63/255],[0 0 128/255],[0 0 255/255];[0 0 0],[.4 .4 .4],[.7 .7 .7]};
% blue for low - black for high >> CS+
colorsRepCSM={[63/255 0 63/255],[128/255 0 128/255],[255/255 0 255/255];[63/255 0 0],[128/255 0 0],[255/255 0 0]};
% purple/magenta for low - red for high >> CS-

% cue align
tt=gtt(1,:); tt=tt(:,21:141);
if seshType==2 && cue==2
else
    %% low vs. high lick rate in CS+ trials before reward delivery (1,2,3rd)
    tempFig=setFigure; hold on; % still plotting neural data, but seperate the trials based on licking rate of that trial
subplot(2,2,1);    
for i=1:3
    shadedErrorBar_CV(tt,mean(low_preRHz_CSP_events{1,i}(21:141,:),2,'omitnan').*(1000./frameRateHz), (std(low_preRHz_CSP_events{1,i}(21:141,:),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps',{ 'color', colorsRepCSP{1,i}});
    hold on;
    errorbar(mean(((low_preRHz_CSP_lickStart{1,i}-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(low_preRHz_CSP_events{1,i}(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(high_preRHz_CSP_events{1,i}(21:141,:),3),2).*(1000./frameRateHz)])),std(((low_preRHz_CSP_lickStart{1,i}-prewin_frames).*(1000./frameRateHz)./sqrt(size(low_preRHz_CSP_lickStart{1,i},2))),[],2,'omitnan'),'horizontal','Color', colorsRepCSP{1,i});    
end
    xlabel('Time from cue (s)');
    ylabel('Cspk rate (Hz)');
    ylim([0 10]);
    title(['CS+ pre-cue: 1st>' num2str(round(nanmean(nanmean(low_preRHz_CSP_licks{1,1},2),1),2)) ' Hz (dark) | 2nd> ' num2str(round(nanmean(nanmean(low_preRHz_CSP_licks{1,2},2),1),2)) ' Hz (med) | 3rd> ' num2str(round(nanmean(nanmean(low_preRHz_CSP_licks{1,3},2),1),2)) ' Hz (light) ']);
subplot(2,2,3);
for i=1:3
    shadedErrorBar_CV(tt,mean(high_preRHz_CSP_events{1,i}(21:141,:),2,'omitnan').*(1000./frameRateHz), (std(high_preRHz_CSP_events{1,i}(21:141,:),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps',{'color', colorsRepCSP{2,i}});
    hold on;
    errorbar(mean(((high_preRHz_CSP_lickStart{1,i}-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(low_preRHz_CSP_events{1,i}(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(high_preRHz_CSP_events{1,i}(21:141,:),3),2).*(1000./frameRateHz)])),std(((high_preRHz_CSP_lickStart{1,i}-prewin_frames).*(1000./frameRateHz)./sqrt(size(high_preRHz_CSP_lickStart{1,i},2))),[],2,'omitnan'),'horizontal','Color', colorsRepCSP{2,i});
end 
    xlabel('Time from cue (s)');
    ylabel('Cspk rate (Hz)');
    ylim([0 10]);
    title(['CS+ pre-cue: 1st>' num2str(round(nanmean(nanmean(high_preRHz_CSP_licks{1,1},2),1),2)) ' Hz (dark) | 2nd> ' num2str(round(nanmean(nanmean(high_preRHz_CSP_licks{1,2},2),1),2)) ' Hz (med) | 3rd> ' num2str(round(nanmean(nanmean(high_preRHz_CSP_licks{1,3},2),1),2)) ' Hz (light) ']);
    sgtitle(['CS+ trials with the lowest (blu) and highest (blk) lick rate PRE-reward delivery | preLow: ' num2str(size(low_preRHz_CSP_licks{1,1},2)) ',' num2str(size(low_preRHz_CSP_licks{1,2},2)) ',' num2str(size(low_preRHz_CSP_licks{1,3},2)) ' - preHigh: '  num2str(size(high_preRHz_CSP_licks{1,1},2)) ',' num2str(size(high_preRHz_CSP_licks{1,2},2)) ',' num2str(size(high_preRHz_CSP_licks{1,3},2))]);
subplot(2,2,2);
for i=1:3
    shadedErrorBar_CV(tt,mean(low_preRHz_CSP_licks{1,i}(21:141,:),2,'omitnan'), (std(low_preRHz_CSP_licks{1,i}(21:141,:),[],2,'omitnan'))./sqrt(size(low_preRHz_CSP_licks{1,i},2)),'lineProps',{ 'color', colorsRepCSP{1,i}});
    hold on;
end
    vline(767,'k');
    xlabel('Time from cue (s)');
    ylabel('Lick rate (Hz)');
    title('Raw lick rate trace for low frequency (blu) pre-reward lick trials');
subplot(2,2,4);
for i=1:3
    shadedErrorBar_CV(tt,mean(high_preRHz_CSP_licks{1,i}(21:141,:),2,'omitnan'), (std(high_preRHz_CSP_licks{1,i}(21:141,:),[],2,'omitnan'))./sqrt(size(high_preRHz_CSP_licks{1,i},2)),'lineProps',{ 'color', colorsRepCSP{2,i}});
    hold on;
end
    vline(767,'k');
    xlabel('Time from cue (s)');
    ylabel('Lick rate (Hz)');
    title('Raw lick rate trace for high frequency (blk) pre-reward lick trials');

if seshType ~= 2
    %savefig(fullfile(output_fn, 'rep_lowHighLickRate_PreCueHz_CSPlusCueAlign.fig'));
    saveas(tempFig, [output_fn  'rep_lowHighLickRate_PreCueHz_CSPlusCueAlign.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName 'rep_lowHighLickRate_PreCueHz_CSPlusCueAlign.fig']));
    saveas(tempFig, [output_fn blockName 'rep_lowHighLickRate_PreCueHz_CSPlusCueAlign.pdf'],'pdf');
end

%% moved from 1350
    % piezo align to trials with most and least lick rate
    if seshType~=2 || seshType==2 && cue==1
    for i = 1:size(animalTrials,2)
        highPost_lickRate_CSPind{i,1} = zeros(1,animalTrials(1,i)); highPre_lickRate_CSPind{i,1} = zeros(1,animalTrials(1,i)); lowPost_lickRate_CSPind{i,1} = zeros(1,animalTrials(1,i)); lowPre_lickRate_CSPind{i,1} = zeros(1,animalTrials(1,i));
        
        highPost_lickRate_CSPind{i,1}(1,high_postRewHz_CSP_ind{1,i}) = 1;
        highPre_lickRate_CSPind{i,1}(1,high_preRewHz_CSP_ind{1,i}) = 1;
        lowPost_lickRate_CSPind{i,1}(1,low_postRewHz_CSP_ind{1,i}) = 1;
        lowPre_lickRate_CSPind{i,1}(1,low_preRewHz_CSP_ind{1,i}) = 1;    
    end
    highPost_lickRate_CSP = cell2mat(highPost_lickRate_CSPind');
    highPre_lickRate_CSP = cell2mat(highPre_lickRate_CSPind');
    lowPost_lickRate_CSP = cell2mat(lowPost_lickRate_CSPind');
    lowPre_lickRate_CSP = cell2mat(lowPre_lickRate_CSPind');
    
    csPlusFirstPostRew_piezoAlign_mostLickRate = firstPostRew_piezoAlignEvents(:,:,~gBlock2 & logical(highPost_lickRate_CSP));
    csPlusFirstPostRew_piezoAlign_leastLickRate = firstPostRew_piezoAlignEvents(:,:,~gBlock2 & logical(lowPost_lickRate_CSP));
    csPlusLastPreRew_piezoAlign_mostLickRate = lastPreRew_piezoAlignEvents(:,:,~gBlock2 & logical(highPre_lickRate_CSP));
    csPlusLastPreRew_piezoAlign_leastLickRate = lastPreRew_piezoAlignEvents(:,:,~gBlock2 & logical(lowPre_lickRate_CSP));
    
     tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
    subplot(2,2,1);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(nanmean(csPlusFirstPostRew_piezoAlign_mostLickRate,3),2).*(1000./frameRateHz), (nanstd(nanmean(csPlusFirstPostRew_piezoAlign_mostLickRate,3),[],2).*(1000./frameRateHz))./sqrt(size(csPlusFirstPostRew_piezoAlign_mostLickRate,2)),'lineProps','k');
    shadedErrorBar_CV(tl_rew, nanmean(nanmean(csPlusFirstPostRew_piezoAlign_leastLickRate,3),2).*(1000./frameRateHz), (nanstd(nanmean(csPlusFirstPostRew_piezoAlign_leastLickRate,3),[],2).*(1000./frameRateHz))./sqrt(size(csPlusFirstPostRew_piezoAlign_leastLickRate,2)),'lineProps','b');
    hold on;
    title(['Cspk first movement after reward from [highest-' num2str(size(csPlusFirstPostRew_piezoAlign_mostLickRate==1,3)) ' | lowest-' num2str(size(csPlusFirstPostRew_piezoAlign_leastLickRate==1,3)) '] lick rate trials']);
    xlabel('Time from lick (ms)');
    ylabel('Cspk rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,3);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(nanmean(csPlusLastPreRew_piezoAlign_mostLickRate,3),2).*(1000./frameRateHz), (nanstd(nanmean(csPlusLastPreRew_piezoAlign_mostLickRate,3),[],2).*(1000./frameRateHz))./sqrt(size(csPlusLastPreRew_piezoAlign_mostLickRate,2)),'lineProps','k');
    shadedErrorBar_CV(tl_rew, nanmean(nanmean(csPlusLastPreRew_piezoAlign_leastLickRate,3),2).*(1000./frameRateHz), (nanstd(nanmean(csPlusLastPreRew_piezoAlign_leastLickRate,3),[],2).*(1000./frameRateHz))./sqrt(size(csPlusLastPreRew_piezoAlign_leastLickRate,2)),'lineProps','b');
    title(['Cspk last movement before reward from [highest-' num2str(size(csPlusLastPreRew_piezoAlign_mostLickRate==1,3)) ' | lowest-' num2str(size(csPlusLastPreRew_piezoAlign_leastLickRate==1,3)) '] lick rate trials']);
    xlabel('Time from lick (ms)');
    ylabel('Cspk rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,2);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(nanmean(csPlusFirstPostRew_piezoAlign_mostLickRate,3),2).*(1000./frameRateHz), (nanstd(nanmean(csPlusFirstPostRew_piezoAlign_mostLickRate,3),[],2).*(1000./frameRateHz))./sqrt(size(csPlusFirstPostRew_piezoAlign_mostLickRate,2)),'lineProps','r');
    shadedErrorBar_CV(tl_rew, nanmean(nanmean(csPlusLastPreRew_piezoAlign_mostLickRate,3),2).*(1000./frameRateHz), (nanstd(nanmean(csPlusLastPreRew_piezoAlign_mostLickRate,3),[],2).*(1000./frameRateHz))./sqrt(size(csPlusLastPreRew_piezoAlign_mostLickRate,2)),'lineProps','g');
    title(['Cspk [first-' num2str(size(csPlusFirstPostRew_piezoAlign_mostLickRate==1,3)) ' | last-' num2str(size(csPlusLastPreRew_piezoAlign_mostLickRate==1,3)) ' movement] from highest lick rate trials']);
    xlabel('Time from lick (ms)');
    ylabel('Cspk rate (Hz)');
    ylim([0 ymax]);
    subplot(2,2,4);hold on;
    shadedErrorBar_CV(tl_rew, nanmean(nanmean(csPlusFirstPostRew_piezoAlign_leastLickRate,3),2).*(1000./frameRateHz), (nanstd(nanmean(csPlusFirstPostRew_piezoAlign_leastLickRate,3),[],2).*(1000./frameRateHz))./sqrt(size(csPlusFirstPostRew_piezoAlign_leastLickRate,2)),'lineProps','r');
    shadedErrorBar_CV(tl_rew, nanmean(nanmean(csPlusLastPreRew_piezoAlign_leastLickRate,3),2).*(1000./frameRateHz), (nanstd(nanmean(csPlusLastPreRew_piezoAlign_leastLickRate,3),[],2).*(1000./frameRateHz))./sqrt(size(csPlusLastPreRew_piezoAlign_leastLickRate,2)),'lineProps','g');
    title(['Cspk [first-' num2str(size(csPlusFirstPostRew_piezoAlign_leastLickRate==1,3)) ' | last-' num2str(size(csPlusLastPreRew_piezoAlign_leastLickRate==1,3)) ' movement] from lowest lick rate trials']);
    xlabel('Time from lick (ms)');
    ylabel('Cspk rate (Hz)');
    ylim([0 ymax]);
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Licks relative to reward']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, '_csPlus_motionAlign_highAndLowLickRate.fig'));
    saveas(tempFig, [output_fn  '_csPlus_motionAlign_highAndLowLickRate.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName '_csPlus_motionAlign_highAndLowLickRate.fig']));
    saveas(tempFig, [output_fn blockName '_csPlus_motionAlign_highAndLowLickRate.pdf'],'pdf');
    end
    end
    %% low vs. high lick rate in CS+ trials after reward delivery (1,2,3rd)
tempFig=setFigure; hold on; % still plotting neural data, but seperate the trials based on licking rate of that trial
subplot(2,2,1);
for i=1:3
    shadedErrorBar_CV(tt,nanmean(low_postRHz_CSP_events{1,i}(21:141,:),2).*(1000./frameRateHz), (nanstd(low_postRHz_CSP_events{1,i}(21:141,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps',{ 'color', colorsRepCSP{1,i}});
    hold on;
    errorbar(mean(((low_postRHz_CSP_lickStart{1,i}-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(low_postRHz_CSP_events{1,i}(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(high_postRHz_CSP_events{1,i}(21:141,:),3),2).*(1000./frameRateHz)])),std(((low_postRHz_CSP_lickStart{1,i}-prewin_frames).*(1000./frameRateHz)./sqrt(size(low_postRHz_CSP_lickStart{1,i},2))),[],2,'omitnan'),'horizontal','Color', colorsRepCSP{1,i});    
end
    xlabel('Time from cue (s)');
    ylabel('Cspk rate (Hz)');
    ylim([0 10]);
    title(['CS+ post-cue: 1st>' num2str(round(nanmean(nanmean(low_postRHz_CSP_licks{1,1},2),1),2)) ' Hz (dark) | 2nd> ' num2str(round(nanmean(nanmean(low_postRHz_CSP_licks{1,2},2),1),2)) ' Hz (med) | 3rd> ' num2str(round(nanmean(nanmean(low_postRHz_CSP_licks{1,3},2),1),2)) ' Hz (light) ']);
subplot(2,2,3);
for i=1:3
    shadedErrorBar_CV(tt,nanmean(high_postRHz_CSP_events{1,i}(21:141,:),2).*(1000./frameRateHz), (nanstd(high_postRHz_CSP_events{1,i}(21:141,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps',{ 'color', colorsRepCSP{2,i}});
    hold on;
    errorbar(mean(((high_postRHz_CSP_lickStart{1,i}-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(low_postRHz_CSP_events{1,i}(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(high_postRHz_CSP_events{1,i}(21:141,:),3),2).*(1000./frameRateHz)])),std(((high_postRHz_CSP_lickStart{1,i}-prewin_frames).*(1000./frameRateHz)./sqrt(size(high_postRHz_CSP_lickStart{1,i},2))),[],2,'omitnan'),'horizontal','Color', colorsRepCSP{2,i});
end
    xlabel('Time from cue (s)');
    ylabel('Cspk rate (Hz)');
    ylim([0 10]);
    title(['CS+ post-cue: 1st>' num2str(round(nanmean(nanmean(high_postRHz_CSP_licks{1,1},2),1),2)) ' Hz (dark) | 2nd> ' num2str(round(nanmean(nanmean(high_postRHz_CSP_licks{1,2},2),1),2)) ' Hz (med) | 3rd> ' num2str(round(nanmean(nanmean(high_postRHz_CSP_licks{1,3},2),1),2)) ' Hz (light) ']);
    sgtitle(['CS+ trials with the lowest (blu) and highest (blk) lick rate POST-reward delivery | postLow: ' num2str(size(low_postRHz_CSP_licks{1,1},2)) ',' num2str(size(low_postRHz_CSP_licks{1,2},2)) ',' num2str(size(low_postRHz_CSP_licks{1,3},2)) ' - postHigh: '  num2str(size(high_postRHz_CSP_licks{1,1},2)) ',' num2str(size(high_postRHz_CSP_licks{1,2},2)) ',' num2str(size(high_postRHz_CSP_licks{1,3},2))]);
subplot(2,2,2);
for i=1:3
    shadedErrorBar_CV(tt,mean(low_postRHz_CSP_licks{1,i}(21:141,:),2,'omitnan'), (std(low_postRHz_CSP_licks{1,i}(21:141,:),[],2,'omitnan'))./sqrt(size(low_postRHz_CSP_licks{1,i},2)),'lineProps',{ 'color', colorsRepCSP{1,i}});
    hold on;
end
    vline(767,'k');
    xlabel('Time from cue (s)');
    ylabel('Lick rate (Hz)');
    title('Raw lick rate trace for low frequency (blu) post-reward lick trials');
subplot(2,2,4);
for i=1:3
    shadedErrorBar_CV(tt,mean(high_postRHz_CSP_licks{1,i}(21:141,:),2,'omitnan'), (std(high_postRHz_CSP_licks{1,i}(21:141,:),[],2,'omitnan'))./sqrt(size(high_postRHz_CSP_licks{1,i},2)),'lineProps',{ 'color', colorsRepCSP{2,i}});
    hold on;
end
    vline(767,'k');
    xlabel('Time from cue (s)');
    ylabel('Lick rate (Hz)');
    title('Raw lick rate trace for high frequency (blk) post-reward lick trials');
    
    if seshType ~= 2
    %savefig(fullfile(output_fn, 'rep_lowHighLickRate_PostCueHz_CSPlusCueAlign.fig'));
    saveas(tempFig, [output_fn  'rep_lowHighLickRate_PostCueHz_CSPlusCueAlign.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName 'rep_lowHighLickRate_PostCueHz_CSPlusCueAlign.fig']));
    saveas(tempFig, [output_fn blockName 'rep_lowHighLickRate_PostCueHz_CSPlusCueAlign.pdf'],'pdf');
    end
end
if seshType==2 && cue==1
else
    
        %% low vs. high lick rate in CS- trials before reward delivery (1,2,3rd)
tempFig=setFigure; hold on; % still plotting neural data, but seperate the trials based on licking rate of that trial
subplot(2,2,1); hold on;
for i=1:3
    shadedErrorBar_CV(tt,nanmean(low_preRHz_CSM_events{1,i}(21:141,:),2).*(1000./frameRateHz), (nanstd(low_preRHz_CSM_events{1,i}(21:141,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps',{ 'color', colorsRepCSM{1,i}});
    hold on;
    errorbar(mean(((low_preRHz_CSM_lickStart{1,i}-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(low_preRHz_CSM_events{1,i}(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(high_preRHz_CSM_events{1,i}(21:141,:),3),2).*(1000./frameRateHz)])),std(((low_preRHz_CSM_lickStart{1,i}-prewin_frames).*(1000./frameRateHz)./sqrt(size(low_preRHz_CSM_lickStart{1,i},2))),[],2,'omitnan'),'horizontal','Color', colorsRepCSM{1,i});
end
    xlabel('Time from cue (s)');
    ylabel('Cspk rate (Hz)');
    ylim([0 10]);
    title(['CS- pre-cue: 1st>' num2str(round(nanmean(nanmean(low_preRHz_CSM_licks{1,1},2),1),2)) ' Hz (dark) | 2nd> ' num2str(round(nanmean(nanmean(low_preRHz_CSM_licks{1,2},2),1),2)) ' Hz (med) | 3rd> ' num2str(round(nanmean(nanmean(low_preRHz_CSM_licks{1,3},2),1),2)) ' Hz (light) ']);
subplot(2,2,3);
for i=1:3
    shadedErrorBar_CV(tt,nanmean(high_preRHz_CSM_events{1,i}(21:141,:),2).*(1000./frameRateHz), (nanstd(high_preRHz_CSM_events{1,i}(21:141,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps',{ 'color', colorsRepCSM{2,i}});
    hold on;
    errorbar(mean(((high_preRHz_CSM_lickStart{1,i}-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(low_preRHz_CSM_events{1,i}(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(high_preRHz_CSM_events{1,i}(21:141,:),3),2).*(1000./frameRateHz)])),std(((high_preRHz_CSM_lickStart{1,i}-prewin_frames).*(1000./frameRateHz)./sqrt(size(high_preRHz_CSM_lickStart{1,i},2))),[],2,'omitnan'),'horizontal','Color', colorsRepCSM{2,i});    
end
    xlabel('Time from cue (s)');
    ylabel('Cspk rate (Hz)');
    ylim([0 10]);
    title(['CS- pre-cue: 1st>' num2str(round(nanmean(nanmean(high_preRHz_CSM_licks{1,1},2),1),2)) ' Hz (dark) | 2nd> ' num2str(round(nanmean(nanmean(high_preRHz_CSM_licks{1,2},2),1),2)) ' Hz (med) | 3rd> ' num2str(round(nanmean(nanmean(high_preRHz_CSM_licks{1,3},2),1),2)) ' Hz (light) ']);
    sgtitle(['CS- trials with the lowest (mag) and highest (red) lick rate PRE-reward delivery | preLow: ' num2str(size(low_preRHz_CSM_licks{1,1},2)) ',' num2str(size(low_preRHz_CSM_licks{1,2},2)) ',' num2str(size(low_preRHz_CSM_licks{1,3},2)) ' - preHigh: '  num2str(size(high_preRHz_CSM_licks{1,1},2)) ',' num2str(size(high_preRHz_CSM_licks{1,2},2)) ',' num2str(size(high_preRHz_CSM_licks{1,3},2))]);
subplot(2,2,2);
for i=1:3
    shadedErrorBar_CV(tt,mean(low_preRHz_CSM_licks{1,i}(21:141,:),2,'omitnan'), (std(low_preRHz_CSM_licks{1,i}(21:141,:),[],2,'omitnan'))./sqrt(size(low_preRHz_CSM_licks{1,i},2)),'lineProps',{ 'color', colorsRepCSM{1,i}});
    hold on;
end
    vline(767,'k');
    xlabel('Time from cue (s)');
    ylabel('Lick rate (Hz)');
    title('Raw lick rate trace for low frequency (mag) pre-reward lick trials');   
subplot(2,2,4);
for i=1:3
    shadedErrorBar_CV(tt,mean(high_preRHz_CSM_licks{1,i}(21:141,:),2,'omitnan'), (std(high_preRHz_CSM_licks{1,i}(21:141,:),[],2,'omitnan'))./sqrt(size(high_preRHz_CSM_licks{1,i},2)),'lineProps',{ 'color', colorsRepCSM{2,i}});
    hold on;
end
    vline(767,'k');
    xlabel('Time from cue (s)');
    ylabel('Lick rate (Hz)');
    title('Raw lick rate trace for high frequency (red) pre-reward lick trials');
    
    if seshType ~= 2
    %savefig(fullfile(output_fn, 'rep_lowHighLickRate_PreCueHz_CSMinusCueAlign.fig'));
    saveas(tempFig, [output_fn  'rep_lowHighLickRate_PreCueHz_CSMinusCueAlign.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName 'rep_lowHighLickRate_PreCueHz_CSMinusCueAlign.fig']));
    saveas(tempFig, [output_fn blockName 'rep_lowHighLickRate_PreCueHz_CSMinusCueAlign.pdf'],'pdf');
    end
 
        %% low vs. high lick rate in CS- trials after reward delivery (1,2,3rd)
tempFig=setFigure; hold on; % still plotting neural data, but seperate the trials based on licking rate of that trial
subplot(2,2,1);    
for i=1:3
    shadedErrorBar_CV(tt,nanmean(low_postRHz_CSM_events{1,i}(21:141,:),2).*(1000./frameRateHz), (nanstd(low_postRHz_CSM_events{1,i}(21:141,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps',{ 'color', colorsRepCSM{1,i}});
    hold on;
    errorbar(mean(((low_postRHz_CSM_lickStart{1,i}-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(low_postRHz_CSM_events{1,i}(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(high_postRHz_CSM_events{1,i}(21:141,:),3),2).*(1000./frameRateHz)])),std(((low_postRHz_CSM_lickStart{1,i}-prewin_frames).*(1000./frameRateHz)./sqrt(size(low_postRHz_CSM_lickStart{1,i},2))),[],2,'omitnan'),'horizontal','Color', colorsRepCSM{1,i});    
end
    xlabel('Time from cue (s)');
    ylabel('Cspk rate (Hz)');
    ylim([0 10]);
    title(['CS- post-cue: 1st>' num2str(round(nanmean(nanmean(low_postRHz_CSM_licks{1,1},2),1),2)) ' Hz (dark) | 2nd> ' num2str(round(nanmean(nanmean(low_postRHz_CSM_licks{1,2},2),1),2)) ' Hz (med) | 3rd> ' num2str(round(nanmean(nanmean(low_postRHz_CSM_licks{1,3},2),1),2)) ' Hz (light) ']);
subplot(2,2,3);
for i=1:3
    shadedErrorBar_CV(tt,nanmean(high_postRHz_CSM_events{1,i}(21:141,:),2).*(1000./frameRateHz), (nanstd(high_postRHz_CSM_events{1,i}(21:141,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps',{ 'color', colorsRepCSM{2,i}});
    hold on;
    errorbar(mean(((high_postRHz_CSM_lickStart{1,i}-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(low_postRHz_CSM_events{1,i}(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(high_postRHz_CSM_events{1,i}(21:141,:),3),2).*(1000./frameRateHz)])),std(((high_postRHz_CSM_lickStart{1,i}-prewin_frames).*(1000./frameRateHz)./sqrt(size(high_postRHz_CSM_lickStart{1,i},2))),[],2,'omitnan'),'horizontal','Color', colorsRepCSM{2,i});
end
    xlabel('Time from cue (s)');
    ylabel('Cspk rate (Hz)');
    ylim([0 10]);
    title(['CS- post-cue: 1st>' num2str(round(nanmean(nanmean(high_postRHz_CSM_licks{1,1},2),1),2)) ' Hz (dark) | 2nd> ' num2str(round(nanmean(nanmean(high_postRHz_CSM_licks{1,2},2),1),2)) ' Hz (med) | 3rd> ' num2str(round(nanmean(nanmean(high_postRHz_CSM_licks{1,3},2),1),2)) ' Hz (light) ']);
    sgtitle(['CS- trials with the lowest (mag) and highest (red) lick rate POST-reward delivery | postLow: ' num2str(size(low_postRHz_CSM_licks{1,1},2)) ',' num2str(size(low_postRHz_CSM_licks{1,2},2)) ',' num2str(size(low_postRHz_CSM_licks{1,3},2)) ' - postHigh: '  num2str(size(high_postRHz_CSM_licks{1,1},2)) ',' num2str(size(high_postRHz_CSM_licks{1,2},2)) ',' num2str(size(high_postRHz_CSM_licks{1,3},2))]);
subplot(2,2,2);
for i=1:3
    shadedErrorBar_CV(tt,mean(low_postRHz_CSM_licks{1,i}(21:141,:),2,'omitnan'), (std(low_postRHz_CSM_licks{1,i}(21:141,:),[],2,'omitnan'))./sqrt(size(low_postRHz_CSM_licks{1,i},2)),'lineProps',{ 'color', colorsRepCSM{1,i}});
    hold on;
end
    vline(767,'k');
    xlabel('Time from cue (s)');
    ylabel('Lick rate (Hz)');
    title('Raw lick rate trace for low frequency (mag) post-reward lick trials');
subplot(2,2,4);
for i=1:3
    shadedErrorBar_CV(tt,mean(high_postRHz_CSM_licks{1,i}(21:141,:),2,'omitnan'), (std(high_postRHz_CSM_licks{1,i}(21:141,:),[],2,'omitnan'))./sqrt(size(high_postRHz_CSM_licks{1,i},2)),'lineProps',{ 'color', colorsRepCSM{2,i}});
    hold on;
end
    vline(767,'k');
    xlabel('Time from cue (s)');
    ylabel('Lick rate (Hz)');
    title('Raw lick rate trace for high frequency (red) post-reward lick trials');
    
    if seshType ~= 2
    %savefig(fullfile(output_fn, 'rep_lowHighLickRate_PostCueHz_CSMinusCueAlign.fig'));
    saveas(tempFig, [output_fn  'rep_lowHighLickRate_PostCueHz_CSMinusCueAlign.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName 'rep_lowHighLickRate_PostCueHz_CSMinusCueAlign.fig']));
    saveas(tempFig, [output_fn blockName 'rep_lowHighLickRate_PostCueHz_CSMinusCueAlign.pdf'],'pdf');
    end
end   
 
% counts

% lick align








%%end%
%%
tt=gtt(1,:);
if seshType ~= 2
    save(fullfile(output_fn, '_cueAlignLick.mat'), 'firstPostRewLickFrame', ...
        'lastPreRewLickFrame', 'tl_rew', 'firstPostRew_lickAlignEvents', 'lastPreRew_lickAlignEvents', ...
        'lickCueAlign', 'lickBurstStart_postReward', 'lickBurstStart_preReward', 'lickCounterVals', 'lickSearch_frames', 'lickDelay_frames',...
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
        'lickCueAlign', 'lickBurstStart_postReward', 'lickBurstStart_preReward', 'lickCounterVals', 'lickSearch_frames', 'lickDelay_frames',...
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
    tt=tt(1,36:96);
        yMaxPiezotemp = max(max([nanmean(groupAlign_piezo(36:96,:),2)]));
        yMaxLicktemp = max(max([nanmean(groupAlign_lick(36:96,:),2)]));
        yMaxPiezo = round(yMaxPiezotemp,1);
        yMaxLick = round(yMaxLicktemp,1);
        if yMaxPiezo < yMaxPiezotemp; yMaxPiezo = yMaxPiezo+.1; else; end
        if yMaxLick < yMaxLicktemp; yMaxLick=yMaxLick+.1; else; end
        
tempFig=setFigure; hold on;
subplot(2,1,1);hold on;
shadedErrorBar_CV(tt,nanmean(groupAlign_piezo(36:96,:),2),nanstd(groupAlign_piezo(36:96,:),[],2)./sqrt(nTrials),'lineProps','k');
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
ylim([0 yMaxPiezo]);
subplot(2,1,2);hold on;
shadedErrorBar_CV(tt,nanmean(lickCueAlign(36:96,:),2),nanstd(lickCueAlign(36:96,:),[],2)./sqrt(nTrials),'lineProps','k');
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
ylim([0 yMaxLick]);
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | lickAlign and piezoAlign | whole trial']);
if seshType ~= 2
%savefig(fullfile(output_fn, '_avgTrialPiezo_abs.fig'));
saveas(tempFig, [output_fn  '_avgTrialPiezo_abs.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_avgTrialPiezo_abs.fig']));
saveas(tempFig, [output_fn blockName '_avgTrialPiezo_abs.pdf'],'pdf');
end
    tt=gtt(1,:);
elseif seshType ~= 2
    intB2Piezo = []; intRewPiezo = []; intB2Lick = []; intRewLick = [];
    for mouseIdx = 1:length(expt)
       intB2Piezo = [intB2Piezo groupAlign.piezo{1,mouseIdx}(:,logical(cell2mat(block2(1,mouseIdx))))];
       intRewPiezo = [intRewPiezo groupAlign.piezo{1,mouseIdx}(:,logical(cell2mat(rew(1,mouseIdx))))];
       intB2Lick = [intB2Lick groupAlign.licks{1,mouseIdx}(:,logical(cell2mat(block2(1,mouseIdx))))];
       intRewLick = [intRewLick groupAlign.licks{1,mouseIdx}(:,logical(cell2mat(rew(1,mouseIdx))))]; 
    end
    
    tt=tt(1,36:96); 
        yMaxPiezotemp = max(max([nanmean(intB2Piezo(36:96,:),2) nanmean(intRewPiezo(36:96,:),2)]));
        yMaxLicktemp = max(max([nanmean(intB2Lick(36:96,:),2) nanmean(intRewLick(36:96,:),2)]));
        yMaxPiezo = round(yMaxPiezotemp,1);
        yMaxLick = round(yMaxLicktemp,1);
        if yMaxPiezo < yMaxPiezotemp; yMaxPiezo = yMaxPiezo+.1; else; end
        if yMaxLick < yMaxLicktemp; yMaxLick=yMaxLick+.1; else; end
        
tempFig=setFigure; hold on;
subplot(2,2,1);hold on;
shadedErrorBar_CV(tt,nanmean(intB2Piezo(36:96,:),2),nanstd(intB2Piezo(36:96,:),[],2)./sqrt(nTrials),'lineProps','k');
ylabel('Piezo voltage');
xlabel('Time from cue (ms)');
ylim([0 yMaxPiezo]);

vline(0,'r');
vline(767,'k');
subplot(2,2,3);hold on;
shadedErrorBar_CV(tt,nanmean(intB2Lick(36:96,:),2),nanstd(intB2Lick(36:96,:),[],2)./sqrt(nTrials),'lineProps','k');
ylabel('Lick rate');
xlabel('Time from cue (ms)');
ylim([0 yMaxLick])

vline(0,'r');
vline(767,'k');

subplot(2,2,2);hold on;
shadedErrorBar_CV(tt,nanmean(intRewPiezo(36:96,:),2),nanstd(intRewPiezo(36:96,:),[],2)./sqrt(nTrials),'lineProps','k');
ylabel('Piezo voltage');
xlabel('Time from cue (ms)');
ylim([0 yMaxPiezo]);

vline(0,'g');
vline(767,'k');
subplot(2,2,4);hold on;
shadedErrorBar_CV(tt,nanmean(intRewLick(36:96,:),2),nanstd(intRewLick(36:96,:),[],2)./sqrt(nTrials),'lineProps','k');
ylabel('Lick rate');
xlabel('Time from cue (ms)');
ylim([0 yMaxLick])

vline(0,'g');
vline(767,'k');
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | lickAlign and piezoAlign | whole trial']);
if seshType ~= 2
%savefig(fullfile(output_fn, '_avgTrialPiezo_abs.fig'));
saveas(tempFig, [output_fn  '_avgTrialPiezo_abs.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_avgTrialPiezo_abs.fig']));
saveas(tempFig, [output_fn blockName '_avgTrialPiezo_abs.pdf'],'pdf');
end

tt=gtt(1,:); tt=tt(1,21:141);

tempFig=setFigure; hold on;
subplot(2,2,1);hold on;
shadedErrorBar_CV(tt,nanmean(intB2Piezo(21:141,:),2),nanstd(intB2Piezo(21:141,:),[],2)./sqrt(nTrials),'lineProps','k');
ylabel('Piezo voltage');
xlabel('Time from cue (ms)');
ylim([0 yMaxPiezo]);

vline(0,'r');
vline(767,'k');
subplot(2,2,3);hold on;
shadedErrorBar_CV(tt,nanmean(intB2Lick(21:141,:),2),nanstd(intB2Lick(21:141,:),[],2)./sqrt(nTrials),'lineProps','k');
ylabel('Lick rate');
xlabel('Time from cue (ms)');
ylim([0 yMaxLick])

vline(0,'r');
vline(767,'k');

subplot(2,2,2);hold on;
shadedErrorBar_CV(tt,nanmean(intRewPiezo(21:141,:),2),nanstd(intRewPiezo(21:141,:),[],2)./sqrt(nTrials),'lineProps','k');
ylabel('Piezo voltage');
xlabel('Time from cue (ms)');
ylim([0 yMaxPiezo]);

vline(0,'g');
vline(767,'k');
subplot(2,2,4);hold on;
shadedErrorBar_CV(tt,nanmean(intRewLick(21:141,:),2),nanstd(intRewLick(21:141,:),[],2)./sqrt(nTrials),'lineProps','k');
ylabel('Lick rate');
xlabel('Time from cue (ms)');
ylim([0 yMaxLick])

vline(0,'g');
vline(767,'k');
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | lickAlign and piezoAlign | whole trial']);
if seshType ~= 2
%savefig(fullfile(output_fn, '_avgTrialPiezo_ablick.fig'));
saveas(tempFig, [output_fn  '_avgTrialPiezo_absext.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_avgTrialPiezo_absext.fig']));
saveas(tempFig, [output_fn blockName '_avgTrialPiezo_absext.pdf'],'pdf');
end
end

tt=gtt(1,:); 
 
preRew_lickSearchRange = prewin_frames+lickDelay_frames:prewin_frames+lickDelay_frames+rewDelay_frames;
postRew_lickSearchRange = prewin_frames+rewDelay_frames+lickDelay_frames:prewin_frames+rewDelay_frames+rewDelay_frames;
preRew_piezoAmp = mean(groupAlign_piezo(preRew_lickSearchRange,:),1);% average across frames for each trial
postRew_piezoAmp = mean(groupAlign_piezo(postRew_lickSearchRange,:),1);

tempFig=setFigure; hold on;
subplot(1,2,1);hold on;
scatter(preRew_lickBurstHz, preRew_piezoAmp,'ok'); 
xlabel('Lick rate');
ylabel('Piezo voltage');
title('Pre reward');
xlim([0 10]);
ylim([-0.2 0.2]);
axis square;
subplot(1,2,2);hold on;
scatter(postRew_lickBurstHz, postRew_piezoAmp,'ok');
xlabel('Lick rate');
ylabel('Piezo voltage');
title('Post reward');
xlim([0 10]);
ylim([-0.2 0.2]);
axis square;
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Lick vs. Piezo | pre- and post-reward']);
if seshType ~= 2
%savefig(fullfile(output_fn, '_LickvsPiezo_abs.fig'));
saveas(tempFig, [output_fn  '_LickvsPiezo_abs.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_LickvsPiezo_abs.fig']));
saveas(tempFig, [output_fn blockName '_LickvsPiezo_abs.pdf'],'pdf');
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
tempFig=setFigure; hold on;
subplot(2,1,1);hold on;
shadedErrorBar_CV(tt,nanmean(low25PreRewEvent,2).*(1000./frameRateHz), (nanstd(low25PreRewEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
hold on;
shadedErrorBar_CV(tt,nanmean(high25PreRewEvent,2).*(1000./frameRateHz), (nanstd(high25PreRewEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
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
subplot(2,1,2);hold on;
shadedErrorBar_CV(tt,nanmean(low25PostRewEvent,2).*(1000./frameRateHz), (nanstd(low25PostRewEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
hold on;
shadedErrorBar_CV(tt,nanmean(high25PostRewEvent,2).*(1000./frameRateHz), (nanstd(high25PostRewEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
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
%savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoAmp25_abs.fig'))
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoAmp25_abs.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoAmp25_abs.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoAmp25_abs.pdf'],'pdf');
end

tempFig=setFigure; hold on;
subplot(2,1,1);hold on;
shadedErrorBar_CV(tt,nanmean(low10PreRewEvent,2).*(1000./frameRateHz), (nanstd(low10PreRewEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
hold on
shadedErrorBar_CV(tt,nanmean(high10PreRewEvent,2).*(1000./frameRateHz), (nanstd(high10PreRewEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
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
subplot(2,1,2);hold on;
shadedErrorBar_CV(tt,nanmean(low10PostRewEvent,2).*(1000./frameRateHz), (nanstd(low10PostRewEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
hold on;
shadedErrorBar_CV(tt,nanmean(high10PostRewEvent,2).*(1000./frameRateHz), (nanstd(high10PostRewEvent,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
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
%savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoAmp10_abs.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoAmp10_abs.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoAmp10_abs.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoAmp10_abs.pdf'],'pdf');
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
groupAlign_piezo_thresh4{1,mouseIdx} = zeros(size(groupAlign.piezo{1,mouseIdx}));%>4std
piezoSearchRange = prewin_frames+lickDelay_frames:prewin_frames+rewDelay_frames+rewDelay_frames;
piezoStart{1,mouseIdx} = nan(1,length(ntrials));

    %nTrial = max(ntrials);
end
% for each trial, frames that has a voltage>n*std = 1, others = 0
nTrial=0;
for mouseIdx = 1:exptCount
    ntrials = animalTrials(1,mouseIdx);%+nTrial;
for itrial = 1+nTrial:ntrials+nTrial
    tCo = 1:ntrials; tNo = tCo(1,itrial-nTrial); 
    groupAlign_piezo_thresh1{1,mouseIdx}(find(groupAlign.piezo{1,mouseIdx}(prewin_frames+lickDelay_frames:end,tNo)>=(piezo_base_avg{1,mouseIdx} + piezo_base_std{1,mouseIdx}.*1) & groupAlign.piezo{1,mouseIdx}(prewin_frames+lickDelay_frames:end,tNo)<(piezo_base_avg{1,mouseIdx} + piezo_base_std{1,mouseIdx}.*2)),tNo) = 1;
    groupAlign_piezo_thresh2{1,mouseIdx}(find(groupAlign.piezo{1,mouseIdx}(prewin_frames+lickDelay_frames:end,tNo)>=(piezo_base_avg{1,mouseIdx} + piezo_base_std{1,mouseIdx}.*2)),tNo)=1;% & groupAlign.piezo{1,mouseIdx}(prewin_frames:end,tNo)<(piezo_base_avg{1,mouseIdx} + piezo_base_std{1,mouseIdx}.*3)),tNo) = 1;
    groupAlign_piezo_thresh3{1,mouseIdx}(find(groupAlign.piezo{1,mouseIdx}(prewin_frames+lickDelay_frames:end,tNo)>=(piezo_base_avg{1,mouseIdx} + piezo_base_std{1,mouseIdx}.*3)),tNo) = 1;
    groupAlign_piezo_thresh4{1,mouseIdx}(find(groupAlign.piezo{1,mouseIdx}(prewin_frames+lickDelay_frames:end,tNo)>=(piezo_base_avg{1,mouseIdx} + piezo_base_std{1,mouseIdx}.*4)),tNo) = 1;
    groupAlign_piezo_thresh1ALL{1,mouseIdx}(find(groupAlign.piezo{1,mouseIdx}(prewin_frames+lickDelay_frames:end,tNo)>=(piezo_base_avg{1,mouseIdx} + piezo_base_std{1,mouseIdx}.*1)),tNo) = 1;
    groupAlign_piezo_thresh2ALL{1,mouseIdx}(find(groupAlign.piezo{1,mouseIdx}(prewin_frames+lickDelay_frames:end,tNo)>=(piezo_base_avg{1,mouseIdx} + piezo_base_std{1,mouseIdx}.*2)),tNo) = 1;
    if  find(groupAlign_piezo_thresh2{1,mouseIdx}(piezoSearchRange,tNo),1,'first')
        piezoStart{1,mouseIdx}(:,tNo) = prewin_frames+lickDelay_frames+find(groupAlign_piezo_thresh2{1,mouseIdx}(:,tNo),1,'first');% the start of very big movements in each trial (3std>baseline)
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

%%
%Going to look at the Cspk rate for trials with 
%  most motion     least motion
%       |       X       |
%  early licks     late licks

%array building at line 989

earlyLick_mostMotion_csp_Cspk_fin = cell2mat(earlyLick_mostMotion_csp_Cspk); 
lateLick_mostMotion_csp_Cspk_fin = cell2mat(lateLick_mostMotion_csp_Cspk);
earlyLick_leastMotion_csp_Cspk_fin = cell2mat(earlyLick_leastMotion_csp_Cspk);
lateLick_leastMotion_csp_Cspk_fin = cell2mat(lateLick_leastMotion_csp_Cspk);

earlyLick_mostMotion_csm_Cspk_fin = cell2mat(earlyLick_mostMotion_csm_Cspk);
lateLick_mostMotion_csm_Cspk_fin = cell2mat(lateLick_mostMotion_csm_Cspk);
earlyLick_leastMotion_csm_Cspk_fin = cell2mat(earlyLick_leastMotion_csm_Cspk);
lateLick_leastMotion_csm_Cspk_fin = cell2mat(lateLick_leastMotion_csm_Cspk);
    
%plottime
tt=gtt(1,:);
tempFig=setFigure; hold on; subplot(2,1,1);
shadedErrorBar_CV(tt,mean(earlyLick_mostMotion_csp_Cspk_fin,2,'omitnan').*(1000./frameRateHz),(std(earlyLick_mostMotion_csp_Cspk_fin,[],2,'omitnan')./sqrt(size(earlyLick_mostMotion_csp_Cspk_fin,2))).*(1000./frameRateHz),'lineProps','k');
hold on;
shadedErrorBar_CV(tt,mean(lateLick_mostMotion_csp_Cspk_fin,2,'omitnan').*(1000./frameRateHz),(std(lateLick_mostMotion_csp_Cspk_fin,[],2,'omitnan')./sqrt(size(lateLick_mostMotion_csp_Cspk_fin,2))).*(1000./frameRateHz),'lineProps','b');
xlabel('Time from cue (s)');
ylabel('Cspk rate (Hz)');
vline(0,'k');
vline(767,'k');
supertitle(['early lick trials - ' num2str(size(earlyLick_mostMotion_csp_Cspk_fin,2)) ' and late lick trials - ' num2str(size(lateLick_mostMotion_csp_Cspk_fin,2))]);
title('Cspk rate for high motion trials with the earliest (blk) and latest (blu) licks');
subplot(2,1,2); hold on;
shadedErrorBar_CV(tt,mean(mostTrialMovement,2,'omitnan'),(std(mostTrialMovement,[],2,'omitnan')./sqrt(size(mostTrialMovement,2))),'lineProps','k');
hold on;
errorbar(mean(((earlyMostStart-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),round(max(mean(mostTrialMovement,2,'omitnan')),2)+.01,std(((earlyMostStart-prewin_frames).*(1000./frameRateHz)./sqrt(size(earlyMostStart,2))),[],2,'omitnan'),'horizontal','k');
errorbar(mean(((lateMostStart-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),round(max(mean(mostTrialMovement,2,'omitnan')),2)+.01,std(((lateMostStart-prewin_frames).*(1000./frameRateHz)./sqrt(size(lateMostStart,2))),[],2,'omitnan'),'horizontal','b');    
ylim([0 round(max(mean(mostTrialMovement,2,'omitnan')),2)+.05]);
title('piezo trace of most movement with errorbars for early and late licks');

if seshType ~= 2
%savefig(fullfile(output_fn, '_mostMotion_earlyLateLicks_Events.fig'));
saveas(tempFig, [output_fn  '_mostMotion_earlyLateLicks_Events.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_mostMotion_earlyLateLicks_Events.fig']));
saveas(tempFig, [output_fn blockName '_mostMotion_earlyLateLicks_Events.pdf'],'pdf');
end

tempFig=setFigure; hold on; subplot(2,1,1);
shadedErrorBar_CV(tt,mean(earlyLick_leastMotion_csp_Cspk_fin,2,'omitnan').*(1000./frameRateHz),(std(earlyLick_leastMotion_csp_Cspk_fin,[],2,'omitnan')./sqrt(size(earlyLick_leastMotion_csp_Cspk_fin,2))).*(1000./frameRateHz),'lineProps','k');
hold on;
shadedErrorBar_CV(tt,mean(lateLick_leastMotion_csp_Cspk_fin,2,'omitnan').*(1000./frameRateHz),(std(lateLick_leastMotion_csp_Cspk_fin,[],2,'omitnan')./sqrt(size(lateLick_leastMotion_csp_Cspk_fin,2))).*(1000./frameRateHz),'lineProps','b');
xlabel('Time from cue (s)');
ylabel('Cspk rate (Hz)');
vline(0,'k');
vline(767,'k');
supertitle(['early lick trials - ' num2str(size(earlyLick_leastMotion_csp_Cspk_fin,2)) ' and late lick trials - ' num2str(size(lateLick_leastMotion_csp_Cspk_fin,2))]);
title('Cspk rate for low motion trials with the earliest (blk) and latest (blu) licks');
subplot(2,1,2); hold on;
shadedErrorBar_CV(tt,mean(leastTrialMovement,2,'omitnan'),(std(leastTrialMovement,[],2,'omitnan')./sqrt(size(leastTrialMovement,2))),'lineProps','b');
hold on;
errorbar(mean(((earlyLeastStart-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),round(max(mean(mostTrialMovement,2,'omitnan')),2)+.01,std(((earlyLeastStart-prewin_frames).*(1000./frameRateHz)./sqrt(size(earlyLeastStart,2))),[],2,'omitnan'),'horizontal','k');
errorbar(mean(((lateLeastStart-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),round(max(mean(mostTrialMovement,2,'omitnan')),2)+.01,std(((lateLeastStart-prewin_frames).*(1000./frameRateHz)./sqrt(size(lateLeastStart,2))),[],2,'omitnan'),'horizontal','b');    
ylim([0 round(max(mean(mostTrialMovement,2,'omitnan')),2)+.05]);
title('piezo trace of least movement with errorbars for early and late licks');

if seshType ~= 2
%savefig(fullfile(output_fn, '_leastMotion_earlyLateLicks_Events.fig'));
saveas(tempFig, [output_fn  '_leastMotion_earlyLateLicks_Events.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_leastMotion_earlyLateLicks_Events.fig']));
saveas(tempFig, [output_fn blockName '_leastMotion_earlyLateLicks_Events.pdf'],'pdf');
end


for i = 1:size(animalTrials,2)
    licksMostCSM_t{i,1} = groupAlign.licks{1,i}(:,logical(mostMotion_csm_ind{1,i}));
    licksLeastCSM_t{i,1} = groupAlign.licks{1,i}(:,logical(leastMotion_csm_ind{1,i}));
    licksMostCSP_t{i,1} = groupAlign.licks{1,i}(:,logical(mostMotion_csp_ind{1,i}));
    licksLeastCSP_t{i,1} = groupAlign.licks{1,i}(:,logical(leastMotion_csp_ind{1,i}));

    licksMostCSM_start{i,1} = groupAlign.lickStart{1,i}(:,logical(mostMotion_csm_ind{1,i}));
    licksLeastCSM_start{i,1} = groupAlign.lickStart{1,i}(:,logical(leastMotion_csm_ind{1,i}));
    licksMostCSP_start{i,1} = groupAlign.lickStart{1,i}(:,logical(mostMotion_csp_ind{1,i}));
    licksLeastCSP_start{i,1} = groupAlign.lickStart{1,i}(:,logical(leastMotion_csp_ind{1,i}));
end

licksMostCSM = cell2mat(licksMostCSM_t');
licksMostCSP = cell2mat(licksMostCSP_t');
licksLeastCSM = cell2mat(licksLeastCSM_t');
licksLeastCSP = cell2mat(licksLeastCSP_t');

licksMostCSMStart = cell2mat(licksMostCSM_start');
licksMostCSPStart = cell2mat(licksMostCSP_start');
licksLeastCSMStart = cell2mat(licksLeastCSM_start');
licksLeastCSPStart = cell2mat(licksLeastCSP_start');


tempFig = setFigure;
subplot(2,1,1);hold on;
shadedErrorBar_CV(tt,mean(licksMostCSM,2,'omitnan'),(std(licksMostCSM,[],2,'omitnan')./sqrt(size(licksMostCSM,2))),'lineProps','k');
hold on;
shadedErrorBar_CV(tt,mean(licksLeastCSM,2,'omitnan'),(std(licksLeastCSM,[],2,'omitnan')./sqrt(size(licksLeastCSM,2))),'lineProps','b');
hold on;
errorbar(mean(((licksMostCSMStart-prewin_frames)),2,'omitnan').*(1000./frameRateHz),.5,std((licksMostCSMStart-prewin_frames),[],2,'omitnan')./sqrt(size(licksMostCSMStart,2)).*(1000./frameRateHz),'horizontal','k');
errorbar(mean(((licksLeastCSMStart-prewin_frames)),2,'omitnan').*(1000./frameRateHz),.5,std((licksLeastCSMStart-prewin_frames),[],2,'omitnan')./sqrt(size(licksLeastCSMStart,2)).*(1000./frameRateHz),'horizontal','b');    
xlabel('Time from cue (s)');
ylabel('Lick Rate (Hz)');
ylim([0 .6]);
vline(0,'k');
vline(767,'k');
title('Lick data for CS- trials with most (blk) and least (blu) movement');
%supertitle('Start of lick bursts for trials with most (blk) and least (blu) movement');
subplot(2,1,2);hold on;
shadedErrorBar_CV(tt,mean(licksMostCSP,2,'omitnan'),(std(licksMostCSP,[],2,'omitnan')./sqrt(size(licksMostCSP,2))),'lineProps','k');
hold on;
shadedErrorBar_CV(tt,mean(licksLeastCSP,2,'omitnan'),(std(licksLeastCSP,[],2,'omitnan')./sqrt(size(licksLeastCSP,2))),'lineProps','b');
hold on;
errorbar(mean(((licksMostCSPStart-prewin_frames)),2,'omitnan').*(1000./frameRateHz),.5,std(((licksMostCSPStart-prewin_frames)),[],2,'omitnan')./sqrt(size(licksMostCSPStart,2)).*(1000./frameRateHz),'horizontal','k');
errorbar(mean(((licksLeastCSPStart-prewin_frames)),2,'omitnan').*(1000./frameRateHz),.5,std(((licksLeastCSPStart-prewin_frames)),[],2,'omitnan')./sqrt(size(licksLeastCSPStart,2)).*(1000./frameRateHz),'horizontal','b');    
xlabel('Time from cue (s)');
ylabel('Lick Rate (Hz)');
ylim([0 .6]);
vline(0,'k');
vline(767,'k');
title('Lick data for CS+ trials with most (blk) and least (blu) movement');

if seshType ~= 2
%savefig(fullfile(output_fn, '_lickRateForMostLeastMotionAlign.fig'));
saveas(tempFig, [output_fn  '_lickRateForMostLeastMotionAlign.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_lickRateForMostLeastMotionAlign.fig']));
saveas(tempFig, [output_fn blockName '_lickRateForMostLeastMotionAlign.pdf'],'pdf');
end

%%
%time to plot most and least movement for CS+ trials and CS- trials
%seperately
tt=gtt(1,:);
tempFig = setFigure;
subplot(2,1,1);hold on;
shadedErrorBar_CV(tt,mean(mostMovementAlignedEvents,2,'omitnan').*(1000./frameRateHz),(std(mostMovementAlignedEvents,[],2,'omitnan')./sqrt(size(mostMovementAlignedEvents,2))).*(1000./frameRateHz),'lineProps','k');
hold on;
shadedErrorBar_CV(tt,mean(leastMovementAlignedEvents,2,'omitnan').*(1000./frameRateHz),(std(leastMovementAlignedEvents,[],2,'omitnan')./sqrt(size(mostMovementAlignedEvents,2))).*(1000./frameRateHz),'lineProps','b');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 inf]);
vline(0,'k');
vline(767,'k');
legend('SEM top 10%','mean top 10%','SEM bot 10%','mean bot 10%','Location','NorthEast','NumColumns',2);
title('Neural data aligned to trials with top 10% (black) and bottom 10% (blue) of movement');
subplot(2,1,2);hold on;
shadedErrorBar_CV(tt,mean(mostTrialMovement,2,'omitnan'),(std(mostTrialMovement,[],2,'omitnan')./sqrt(size(mostTrialMovement,2))),'lineProps','k');
hold on;
shadedErrorBar_CV(tt,mean(leastTrialMovement,2,'omitnan'),(std(leastTrialMovement,[],2,'omitnan')./sqrt(size(mostTrialMovement,2))),'lineProps','b');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
vline(0,'k');
vline(767,'k');
title('Piezo data for trials with top 10% (black) and bottom 10% (blue) of movement');
hold off;
sgtitle([num2str(size(mostTrialMovement,2)) ' | ' num2str(size(leastTrialMovement,2)) ' trials with most and least movement']);
if seshType ~= 2
%savefig(fullfile(output_fn, '_mostLeastMotionAlign.fig'));
saveas(tempFig, [output_fn  '_mostLeastMotionAlign.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_mostLeastMotionAlign.fig']));
saveas(tempFig, [output_fn blockName '_mostLeastMotionAlign.pdf'],'pdf');
end

% CS+
tt = gtt(1,21:141);
if seshType == 2 && sum(gBlock2) > 0
else
tempFig = setFigure; hold on;
subplot(2,1,1);hold on;
shadedErrorBar_CV(tt,mean(mostCSPMovementAlignedEvents(21:141,:),2,'omitnan').*(1000./frameRateHz),(std(mostCSPMovementAlignedEvents(21:141,:),[],2,'omitnan')./sqrt(size(mostMovementAlignedEvents,2))).*(1000./frameRateHz),'lineProps','k');
hold on;
shadedErrorBar_CV(tt,mean(leastCSPMovementAlignedEvents(21:141,:),2,'omitnan').*(1000./frameRateHz),(std(leastCSPMovementAlignedEvents(21:141,:),[],2,'omitnan')./sqrt(size(mostMovementAlignedEvents,2))).*(1000./frameRateHz),'lineProps','b');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 ymax]);
vline(0,'k');
vline(767,'k');
legend('SEM top 25%','mean top 25%','SEM bot 25%','mean bot 25%','Location','NorthEast','NumColumns',2);
title('Neural data aligned to trials with top 25% (black) and bottom 25% (blue) of movement');
subplot(2,1,2);hold on;
shadedErrorBar_CV(tt,mean(mostCSPTrialMovement(21:141,:),2,'omitnan'),(std(mostCSPTrialMovement(21:141,:),[],2,'omitnan')./sqrt(size(mostTrialMovement,2))),'lineProps','k');
hold on;
shadedErrorBar_CV(tt,mean(leastCSPTrialMovement(21:141,:),2,'omitnan'),(std(leastCSPTrialMovement(21:141,:),[],2,'omitnan')./sqrt(size(leastTrialMovement,2))),'lineProps','b');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
vline(0,'k');
vline(767,'k');
title(['Piezo data for CS+ trials with top 20% (' num2str(round(mean(mean(mostCSPTrialMovement(21:141,:),2,'omitnan'),1,'omitnan'),2)) ') and bottom 20% (' num2str(round(mean(mean(leastCSPTrialMovement(21:141,:),2,'omitnan'),1,'omitnan'),2)) ') of movement']);
hold off;
sgtitle([num2str(size(mostCSPTrialMovement,2)) ' | ' num2str(size(leastCSPTrialMovement,2)) ' CS+ trials with most and least movement']);
if seshType ~= 2
%savefig(fullfile(output_fn, '_mostCSPLeastMotionAlign.fig'));
saveas(tempFig, [output_fn  '_mostCSPLeastMotionAlign.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_mostCSPLeastMotionAlign.fig']));
saveas(tempFig, [output_fn blockName '_mostCSPLeastMotionAlign.pdf'],'pdf');
end
end
% CS-
if seshType == 2 && sum(gBlock2) == 0
elseif expType == 3 && cue ~= 3
else
tempFig = figure;
subplot(2,1,1);hold on;
shadedErrorBar_CV(tt,mean(mostCSMMovementAlignedEvents(21:141,:),2,'omitnan').*(1000./frameRateHz),(std(mostCSMMovementAlignedEvents(21:141,:),[],2,'omitnan')./sqrt(size(mostMovementAlignedEvents,2))).*(1000./frameRateHz),'lineProps','k');
hold on;
shadedErrorBar_CV(tt,mean(leastCSMMovementAlignedEvents(21:141,:),2,'omitnan').*(1000./frameRateHz),(std(leastCSMMovementAlignedEvents(21:141,:),[],2,'omitnan')./sqrt(size(mostMovementAlignedEvents,2))).*(1000./frameRateHz),'lineProps','b');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 ymax]);
vline(0,'k');
vline(767,'k');
legend('SEM top 20%','mean top 20%','SEM bot 20%','mean bot 20%','Location','NorthEast','NumColumns',2);
title('Neural data aligned to trials with top 20% (black) and bottom 20% (blue) of movement');
subplot(2,1,2);hold on;
shadedErrorBar_CV(tt,mean(mostCSMTrialMovement(21:141,:),2,'omitnan'),(std(mostCSMTrialMovement(21:141,:),[],2,'omitnan')./sqrt(size(mostTrialMovement,2))),'lineProps','k');
hold on;
shadedErrorBar_CV(tt,mean(leastCSMTrialMovement(21:141,:),2,'omitnan'),(std(leastCSMTrialMovement(21:141,:),[],2,'omitnan')./sqrt(size(leastTrialMovement,2))),'lineProps','b');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
vline(0,'k');
vline(767,'k');
title(['Piezo data for CS+ trials with top 20% (' num2str(round(mean(mean(mostCSMTrialMovement(21:141,:),2,'omitnan'),1,'omitnan'),2)) ') and bottom 20% (' num2str(round(mean(mean(leastCSMTrialMovement(21:141,:),2,'omitnan'),1,'omitnan'),2)) ') of movement']);
hold off;
sgtitle([num2str(size(mostCSMTrialMovement,2)) ' | ' num2str(size(leastCSMTrialMovement,2)) ' CS- trials with most and least movement']);
if seshType ~= 2
%savefig(fullfile(output_fn, '_mostCSMLeastMotionAlign.fig'));
saveas(tempFig, [output_fn  '_mostCSMLeastMotionAlign.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_mostCSMLeastMotionAlign.fig']));
saveas(tempFig, [output_fn blockName '_mostCSMLeastMotionAlign.pdf'],'pdf');
end
end
%least CS- and least CS+
tt = gtt(1,21:141);
if seshType == 2 
else
tempFig = figure;
subplot(2,1,1);hold on;
shadedErrorBar_CV(tt,mean(leastCSMMovementAlignedEvents(21:141,:),2,'omitnan').*(1000./frameRateHz),(std(leastCSMMovementAlignedEvents(21:141,:),[],2,'omitnan')./sqrt(size(leastMovementAlignedEvents,2))).*(1000./frameRateHz),'lineProps','r');
hold on;
shadedErrorBar_CV(tt,mean(leastCSPMovementAlignedEvents(21:141,:),2,'omitnan').*(1000./frameRateHz),(std(leastCSPMovementAlignedEvents(21:141,:),[],2,'omitnan')./sqrt(size(mostMovementAlignedEvents,2))).*(1000./frameRateHz),'lineProps','k');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 ymax]);
vline(0,'k');
vline(767,'k');
legend('SEM bot 25% CS-','mean bot 25% CS-','SEM bot 25% CS+','mean bot 25% CS+','Location','NorthEast','NumColumns',2);
title('Neural data aligned to trials with bot 25% CS- (red) and bot 25% CS+ (black) of movement');
subplot(2,1,2);hold on;
shadedErrorBar_CV(tt,mean(leastCSMTrialMovement(21:141,:),2,'omitnan'),(std(leastCSMTrialMovement(21:141,:),[],2,'omitnan')./sqrt(size(leastTrialMovement,2))),'lineProps','r');
hold on;
shadedErrorBar_CV(tt,mean(leastCSPTrialMovement(21:141,:),2,'omitnan'),(std(leastCSPTrialMovement(21:141,:),[],2,'omitnan')./sqrt(size(leastTrialMovement,2))),'lineProps','k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
vline(0,'k');
vline(767,'k');
title(['Piezo data for trials with bottom 25% CS- (' num2str(round(mean(mean(leastCSMTrialMovement(21:141,:),2,'omitnan'),1,'omitnan'),2)) ') and bottom 25% CS+ (' num2str(round(mean(mean(leastCSPTrialMovement(21:141,:),2,'omitnan'),1,'omitnan'),2)) ') of movement']);
hold off;
sgtitle([num2str(size(leastCSMTrialMovement,2)) ' | ' num2str(size(leastCSPTrialMovement,2)) ' CS-/CS+ trials with least movement']);
if seshType ~= 2
%savefig(fullfile(output_fn, '_leastCSPandCSMMotion.fig'));
saveas(tempFig, [output_fn  '_leastCSPandCSMMotion.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_leastCSPandCSMMotion.fig']));
saveas(tempFig, [output_fn blockName '_leastCSPandCSMMotion.pdf'],'pdf');
end
end
%
%most CS- and most CS+
tt = gtt(1,21:141);
if seshType == 2 
else
tempFig = figure;
subplot(2,1,1);hold on;
shadedErrorBar_CV(tt,mean(mostCSMMovementAlignedEvents(21:141,:),2,'omitnan').*(1000./frameRateHz),(std(mostCSMMovementAlignedEvents(21:141,:),[],2,'omitnan')./sqrt(size(mostMovementAlignedEvents,2))).*(1000./frameRateHz),'lineProps','r');
hold on;
shadedErrorBar_CV(tt,mean(mostCSPMovementAlignedEvents(21:141,:),2,'omitnan').*(1000./frameRateHz),(std(mostCSPMovementAlignedEvents(21:141,:),[],2,'omitnan')./sqrt(size(mostMovementAlignedEvents,2))).*(1000./frameRateHz),'lineProps','k');
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 ymax]);
vline(0,'k');
vline(767,'k');
legend('SEM top 25% CS-','mean top 25% CS-','SEM top 25% CS+','mean top 25% CS+','Location','NorthEast','NumColumns',2);
title('Neural data aligned to trials with top 25% CS- (red) and top 25% CS+ (black) of movement');
subplot(2,1,2);hold on;
shadedErrorBar_CV(tt,mean(mostCSMTrialMovement(21:141,:),2,'omitnan'),(std(mostCSMTrialMovement(21:141,:),[],2,'omitnan')./sqrt(size(mostTrialMovement,2))),'lineProps','r');
hold on;
shadedErrorBar_CV(tt,mean(mostCSPTrialMovement(21:141,:),2,'omitnan'),(std(mostCSPTrialMovement(21:141,:),[],2,'omitnan')./sqrt(size(mostTrialMovement,2))),'lineProps','k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
vline(0,'k');
vline(767,'k');
title(['Piezo data for trials with top 25% CS- (' num2str(round(mean(mean(mostCSMTrialMovement(21:141,:),2,'omitnan'),1,'omitnan'),2)) ') and top 25% CS+ (' num2str(round(mean(mean(mostCSPTrialMovement(21:141,:),2,'omitnan'),1,'omitnan'),2)) ') of movement']);
hold off;
sgtitle([num2str(size(mostCSMTrialMovement,2)) ' | ' num2str(size(mostCSPTrialMovement,2)) ' CS-/CS+ trials with most movement']);
if seshType ~= 2
%savefig(fullfile(output_fn, '_mostCSPandCSMMotion.fig'));
saveas(tempFig, [output_fn  '_mostCSPandCSMMotion.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_mostCSPandCSMMotion.fig']));
saveas(tempFig, [output_fn blockName '_mostCSPandCSMMotion.pdf'],'pdf');
end
end

tt=gtt(1,:);
%% early vs. late piezo events
latestMovementAlignedEvents=[]; earliestMovementAlignedEvents=[]; latestTrialMovementStart=[]; earliestTrialMovementStart=[];

latestCSPMovementAlignedEvents=[]; earliestCSPMovementAlignedEvents=[]; latestCSPTrialMovementStart=[]; earliestCSPTrialMovementStart=[];
latestCSMMovementAlignedEvents=[]; earliestCSMMovementAlignedEvents=[]; latestCSMTrialMovementStart=[]; earliestCSMTrialMovementStart=[];
latestCSPTrialMovement=[]; earliestCSPTrialMovement=[]; latestCSMTrialMovement=[]; earliestCSMTrialMovement=[];
latestTrialMovement=[]; earliestTrialMovement=[];
for mouseIdx = 1:exptCount
    CSPTimePiezo_ind{1,mouseIdx}=[]; CSMTimePiezo_ind{1,mouseIdx}=[]; sigMove_ind{1,mouseIdx}=zeros([1,size(piezoStart{1,mouseIdx},2)]);
    CSPTimePiezo{1,mouseIdx}=[]; CSMTimePiezo{1,mouseIdx}=[];
    %make an logical array to index trials with sig. movement
for i = 1:size(piezoStart{1,mouseIdx},2)
    if piezoStart{1,mouseIdx}(1,i)~=0; sigMove_ind{1,mouseIdx}(1,i) = 1; else; end
end

    %sort the movement data from earliest to latest and remove zeros so 
    %they don't throw off the top/bot 10% calculations
    [timePiezo{1,mouseIdx}, timePiezo_ind{1,mouseIdx}]=sort(piezoStart{1,mouseIdx});
    timePiezo_ind{1,mouseIdx}(1,timePiezo{1,mouseIdx}==0)=NaN; timePiezo_ind{1,mouseIdx}(isnan(timePiezo_ind{1,mouseIdx}))=[];
    timePiezo{1,mouseIdx}(1,timePiezo{1,mouseIdx}==0)=NaN; timePiezo{1,mouseIdx}(isnan(timePiezo{1,mouseIdx}))=[];

    if isempty(timePiezo_ind{1,mouseIdx})
        timePiezo_ind{1,mouseIdx}=NaN; timePiezo{1,mouseIdx}=NaN;
    end
for i = 1:size(timePiezo_ind{1,mouseIdx},2)
    if isnan(timePiezo{1,mouseIdx})
CSPTimePiezo_ind{1,mouseIdx} = [CSPTimePiezo_ind{1,mouseIdx} NaN];
CSPTimePiezo{1,mouseIdx} = [CSPTimePiezo{1,mouseIdx} NaN];
CSMTimePiezo_ind{1,mouseIdx} = [CSMTimePiezo_ind{1,mouseIdx} NaN];
CSMTimePiezo{1,mouseIdx} = [CSMTimePiezo{1,mouseIdx} NaN];        
    elseif ~isempty(find(rewardTrials{1,mouseIdx} == timePiezo_ind{1,mouseIdx}(1,i)))
CSPTimePiezo_ind{1,mouseIdx} = [CSPTimePiezo_ind{1,mouseIdx} timePiezo_ind{1,mouseIdx}(1,i)];
CSPTimePiezo{1,mouseIdx} = [CSPTimePiezo{1,mouseIdx} timePiezo{1,mouseIdx}(1,i)];
    elseif ~isempty(find(block2Trials{1,mouseIdx} == timePiezo_ind{1,mouseIdx}(1,i)))
CSMTimePiezo_ind{1,mouseIdx} = [CSMTimePiezo_ind{1,mouseIdx} timePiezo_ind{1,mouseIdx}(1,i)];
CSMTimePiezo{1,mouseIdx} = [CSMTimePiezo{1,mouseIdx} timePiezo{1,mouseIdx}(1,i)];
    end
end


latestTrialMovement_ind{1,mouseIdx} = timePiezo_ind{1,mouseIdx}(:,size(timePiezo_ind{1,mouseIdx},2)-(floor(size(timePiezo_ind{1,mouseIdx},2)./4)):end);
earliestTrialMovement_ind{1,mouseIdx} = timePiezo_ind{1,mouseIdx}(:,1:ceil(size(timePiezo_ind{1,mouseIdx},2)./4));

if seshType == 2 && sum(gBlock2) > 0
else
latestCSPTrialMovement_ind{1,mouseIdx} = CSPTimePiezo_ind{1,mouseIdx}(:,size(CSPTimePiezo_ind{1,mouseIdx},2)-(floor(size(CSPTimePiezo_ind{1,mouseIdx},2)./4)):end);
earliestCSPTrialMovement_ind{1,mouseIdx} = CSPTimePiezo_ind{1,mouseIdx}(:,1:ceil(size(CSPTimePiezo_ind{1,mouseIdx},2)./4));
if isnan(timePiezo{1,mouseIdx})
elseif ~isnan(timePiezo{1,mouseIdx})
latestCSPMovementAlignedEvents = [latestCSPMovementAlignedEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,latestCSPTrialMovement_ind{1,mouseIdx}),3)];
earliestCSPMovementAlignedEvents = [earliestCSPMovementAlignedEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,earliestCSPTrialMovement_ind{1,mouseIdx}),3)];
latestCSPTrialMovement = [latestCSPTrialMovement groupAlign.piezo{1,mouseIdx}(:,latestCSPTrialMovement_ind{1,mouseIdx})];
earliestCSPTrialMovement = [earliestCSPTrialMovement groupAlign.piezo{1,mouseIdx}(:,earliestCSPTrialMovement_ind{1,mouseIdx})];
end
latestCSPTrialMovementStart = [latestCSPTrialMovementStart CSPTimePiezo{1,mouseIdx}(:,size(CSPTimePiezo_ind{1,mouseIdx},2)-(floor(size(CSPTimePiezo_ind{1,mouseIdx},2)./4)):end)];
earliestCSPTrialMovementStart = [earliestCSPTrialMovementStart CSPTimePiezo{1,mouseIdx}(:,1:ceil(size(CSPTimePiezo_ind{1,mouseIdx},2)./4))];
end
if seshType == 2 && sum(gBlock2) == 0
elseif expType == 3 && cue ~= 3
else
latestCSMTrialMovement_ind{1,mouseIdx} = CSMTimePiezo_ind{1,mouseIdx}(:,size(CSMTimePiezo_ind{1,mouseIdx},2)-(floor(size(CSMTimePiezo_ind{1,mouseIdx},2)./4)):end);
earliestCSMTrialMovement_ind{1,mouseIdx} = CSMTimePiezo_ind{1,mouseIdx}(:,1:ceil(size(CSMTimePiezo_ind{1,mouseIdx},2)./4));
if isnan(timePiezo{1,mouseIdx})
elseif ~isnan(timePiezo{1,mouseIdx})
latestCSMMovementAlignedEvents = [latestCSMMovementAlignedEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,latestCSMTrialMovement_ind{1,mouseIdx}),3)];
earliestCSMMovementAlignedEvents = [earliestCSMMovementAlignedEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,earliestCSMTrialMovement_ind{1,mouseIdx}),3)];
latestCSMTrialMovement = [latestCSMTrialMovement groupAlign.piezo{1,mouseIdx}(:,latestCSMTrialMovement_ind{1,mouseIdx})];
earliestCSMTrialMovement = [earliestCSMTrialMovement groupAlign.piezo{1,mouseIdx}(:,earliestCSMTrialMovement_ind{1,mouseIdx})];
end
latestCSMTrialMovementStart = [latestCSMTrialMovementStart CSMTimePiezo{1,mouseIdx}(:,size(CSMTimePiezo_ind{1,mouseIdx},2)-(floor(size(CSMTimePiezo_ind{1,mouseIdx},2)./4)):end)];
earliestCSMTrialMovementStart = [earliestCSMTrialMovementStart CSMTimePiezo{1,mouseIdx}(:,1:ceil(size(CSMTimePiezo_ind{1,mouseIdx},2)./4))];
end
%
if isnan(timePiezo{1,mouseIdx})
elseif ~isnan(timePiezo{1,mouseIdx})
latestMovementAlignedEvents = [latestMovementAlignedEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,latestTrialMovement_ind{1,mouseIdx}),3)];
earliestMovementAlignedEvents = [earliestMovementAlignedEvents nanmean(groupAlign.events{1,mouseIdx}(:,:,earliestTrialMovement_ind{1,mouseIdx}),3)];
latestTrialMovement = [latestTrialMovement groupAlign.piezo{1,mouseIdx}(:,latestTrialMovement_ind{1,mouseIdx})];
earliestTrialMovement = [earliestTrialMovement groupAlign.piezo{1,mouseIdx}(:,earliestTrialMovement_ind{1,mouseIdx})];
end
latestTrialMovementStart = [latestTrialMovementStart timePiezo{1,mouseIdx}(:,size(timePiezo_ind{1,mouseIdx},2)-(floor(size(timePiezo_ind{1,mouseIdx},2)./10)):end)];
earliestTrialMovementStart = [earliestTrialMovementStart timePiezo{1,mouseIdx}(:,1:ceil(size(timePiezo_ind{1,mouseIdx},2)./10))];
end
%% early late movement    

nIC = size(groupAlign_events,2);
tt=gtt(1,:); tt=tt(:,21:141);
    tempFig=setFigure; hold on;
    subplot(2,1,1); hold on;
    shadedErrorBar_CV(tt,nanmean(earliestMovementAlignedEvents(21:141,:),2).*(1000./frameRateHz), (nanstd(earliestMovementAlignedEvents(21:141,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    hold on;
    shadedErrorBar_CV(tt,nanmean(latestMovementAlignedEvents(21:141,:),2).*(1000./frameRateHz), (nanstd(latestMovementAlignedEvents(21:141,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
    errorbar(mean(((earliestTrialMovementStart-prewin_frames-lickDelay_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(earliestMovementAlignedEvents(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(latestMovementAlignedEvents(21:141,:),3),2).*(1000./frameRateHz)])),std(((earliestTrialMovementStart-prewin_frames-lickDelay_frames).*(1000./frameRateHz)./sqrt(size(earliestTrialMovementStart,2))),[],2,'omitnan'),'horizontal','k');
    errorbar(mean(((latestTrialMovementStart-prewin_frames-lickDelay_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(earliestMovementAlignedEvents(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(latestMovementAlignedEvents(21:141,:),3),2).*(1000./frameRateHz)])),std(((latestTrialMovementStart-prewin_frames-lickDelay_frames).*(1000./frameRateHz)./sqrt(size(latestTrialMovementStart,2))),[],2,'omitnan'),'horizontal','b');    
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    vline(0,'k');
    vline(767,'k');
    title([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | movement bursts: early (n = ' num2str(length(cell2mat(earliestTrialMovement_ind))) '| blk); late (n = ' num2str(length(cell2mat(latestTrialMovement_ind))) '| blu)']); 
    subplot(2,1,2); hold on;
    shadedErrorBar_CV(tt,nanmean(earliestTrialMovement(21:141,:),2), (nanstd(earliestTrialMovement(21:141,:),[],2))./sqrt(size(earliestTrialMovement(21:141,:),2)),'lineProps','k');
    hold on;
    shadedErrorBar_CV(tt,nanmean(latestTrialMovement(21:141,:),2), (nanstd(latestTrialMovement(21:141,:),[],2))./sqrt(size(latestTrialMovement(21:141,:),2)),'lineProps','b');
    xlabel('Time from cue');
    ylabel('Piezo (V)');
    ylim([0 inf]);
    vline(0,'k');
    vline(767,'k');
    title('earliest and latest movement peaks - traces associated with whiskers above');
    if seshType ~=2
    %savefig(fullfile(output_fn, 'earlyLate_25PercentOfMovement_allTrial.fig'));
    saveas(tempFig, [output_fn  'earlyLate_25PercentOfMovement_allTrial.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, ['\' blockName 'earlyLate_25PercentOfMovement_allTrial.fig']));
    saveas(tempFig, [output_fn blockName 'earlyLate_25PercentOfMovement_allTrial.pdf'],'pdf');
    end
    hold off;
%% CS+ early late movement    
if seshType==2 & cue==2
else
nIC = size(groupAlign_events,2);
tt=gtt(1,:); tt=tt(1,21:141);
tempFig=setFigure; hold on;
subplot(2,1,1); hold on;
shadedErrorBar_CV(tt,nanmean(earliestCSPMovementAlignedEvents(21:141,:),2).*(1000./frameRateHz), (nanstd(earliestCSPMovementAlignedEvents(21:141,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
hold on;
shadedErrorBar_CV(tt,nanmean(latestCSPMovementAlignedEvents(21:141,:),2).*(1000./frameRateHz), (nanstd(latestCSPMovementAlignedEvents(21:141,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
errorbar(mean(((earliestCSPTrialMovementStart-prewin_frames-lickDelay_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(earliestCSPMovementAlignedEvents(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(latestCSPMovementAlignedEvents(21:141,:),3),2).*(1000./frameRateHz)])),std(((earliestCSPTrialMovementStart-prewin_frames-lickDelay_frames).*(1000./frameRateHz)./sqrt(size(earliestCSPTrialMovementStart,2))),[],2,'omitnan'),'horizontal','k');
errorbar(mean(((latestCSPTrialMovementStart-prewin_frames-lickDelay_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(earliestCSPMovementAlignedEvents(21:141,:),3),2).*(1000./frameRateHz) nanmean(nanmean(latestCSPMovementAlignedEvents(21:141,:),3),2).*(1000./frameRateHz)])),std(((latestCSPTrialMovementStart-prewin_frames-lickDelay_frames).*(1000./frameRateHz)./sqrt(size(latestCSPTrialMovementStart,2))),[],2,'omitnan'),'horizontal','b');    
xlabel('Time from cue');
ylabel('Spike rate (Hz)');
ylim([0 inf]);
vline(0,'k');
vline(767,'k');
title([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | movement bursts: early (n = ' num2str(length(cell2mat(earliestCSPTrialMovement_ind))) '| blk); late (n = ' num2str(length(cell2mat(latestCSPTrialMovement_ind))) '| blu)']);
subplot(2,1,2); hold on;
shadedErrorBar_CV(tt,nanmean(earliestCSPTrialMovement(21:141,:),2), (nanstd(earliestCSPTrialMovement(21:141,:),[],2))./sqrt(size(earliestCSPTrialMovement(21:141,:),2)),'lineProps','k');
hold on;
shadedErrorBar_CV(tt,nanmean(latestCSPTrialMovement(21:141,:),2), (nanstd(latestCSPTrialMovement(21:141,:),[],2))./sqrt(size(latestCSPTrialMovement(21:141,:),2)),'lineProps','b');
xlabel('Time from cue');
ylabel('Piezo (V)');
ylim([0 inf]);
vline(0,'k');
vline(767,'k');
title('earliest and latest movement peaks - traces associated with whiskers above');
if seshType ~=2
%savefig(fullfile(output_fn, 'earlyLate_25PercentOfMovement_CSPlusTrial.fig'));
saveas(tempFig, [output_fn  'earlyLate_25PercentOfMovement_CSPlusTrial.pdf'],'pdf');
else
%savefig(fullfile(output_fn, ['\' blockName 'earlyLate_25PercentOfMovement_CSPlusTrial.fig']));
saveas(tempFig, [output_fn blockName 'earlyLate_25PercentOfMovement_CSPlusTrial.pdf'],'pdf');
end
hold off;
end
%% CS+ early late movement    
if seshType==2 && cue==1
else
    nIC = size(groupAlign_events,2);
    tt=gtt(1,:); tt=tt(:,21:141);
    tempFig=setFigure; hold on;
    subplot(2,1,1); hold on;
    shadedErrorBar_CV(tt,nanmean(earliestCSMMovementAlignedEvents(21:141,:),2).*(1000./frameRateHz), (nanstd(earliestCSMMovementAlignedEvents(21:141,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    hold on;
    shadedErrorBar_CV(tt,nanmean(latestCSMMovementAlignedEvents(21:141,:),2).*(1000./frameRateHz), (nanstd(latestCSMMovementAlignedEvents(21:141,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
    errorbar(mean(((earliestCSMTrialMovementStart-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(earliestCSMMovementAlignedEvents,3),2).*(1000./frameRateHz) nanmean(nanmean(latestCSMMovementAlignedEvents,3),2).*(1000./frameRateHz)])),std(((earliestCSMTrialMovementStart-prewin_frames).*(1000./frameRateHz)./sqrt(size(earliestCSMTrialMovementStart,2))),[],2,'omitnan'),'horizontal','k');
    errorbar(mean(((latestCSMTrialMovementStart-prewin_frames).*(1000./frameRateHz)),2,'omitnan'),max(max([nanmean(nanmean(earliestCSMMovementAlignedEvents,3),2).*(1000./frameRateHz) nanmean(nanmean(latestCSMMovementAlignedEvents,3),2).*(1000./frameRateHz)])),std(((latestCSMTrialMovementStart-prewin_frames).*(1000./frameRateHz)./sqrt(size(latestCSMTrialMovementStart,2))),[],2,'omitnan'),'horizontal','b');    
    xlabel('Time from cue');
    ylabel('Spike rate (Hz)');
    ylim([0 inf]);
    vline(0,'k');
    vline(767,'k');
    title([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | movement bursts: early (n = ' num2str(length(cell2mat(earliestCSMTrialMovement_ind))) '| blk); late (n = ' num2str(length(cell2mat(latestCSMTrialMovement_ind))) '| blu)']);
    subplot(2,1,2); hold on;
    shadedErrorBar_CV(tt,nanmean(earliestCSMTrialMovement(21:141,:),2), (nanstd(earliestCSMTrialMovement(21:141,:),[],2))./sqrt(size(earliestCSMTrialMovement(21:141,:),2)),'lineProps','k');
    hold on;
    shadedErrorBar_CV(tt,nanmean(latestCSMTrialMovement(21:141,:),2), (nanstd(latestCSMTrialMovement(21:141,:),[],2))./sqrt(size(latestCSMTrialMovement(21:141,:),2)),'lineProps','b');
    xlabel('Time from cue');
    ylabel('Piezo (V)');
    ylim([0 inf]);
    vline(0,'k');
    vline(767,'k');
    title('earliest and latest movement peaks - traces associated with whiskers above');
    if seshType ~=2
    %savefig(fullfile(output_fn, 'earlyLate_25PercentOfMovement_CSMinusTrial.fig'));
    saveas(tempFig, [output_fn  'earlyLate_25PercentOfMovement_CSMinusTrial.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, ['\' blockName 'earlyLate_25PercentOfMovement_CSMinusTrial.fig']));
    saveas(tempFig, [output_fn blockName 'earlyLate_25PercentOfMovement_CSMinusTrial.pdf'],'pdf');
    end
    hold off;
end
tt=gtt(1,:);
%% pulling piezo data without considering frames with motion >=< mean motion
tempFig=setFigure; hold on;
subplot(2,1,1);hold on;
shadedErrorBar_CV(tt,nanmean(earlyPiezoEventB2,2).*(1000./frameRateHz), (nanstd(earlyPiezoEventB2,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
hold on;
shadedErrorBar_CV(tt,nanmean(latePiezoEventB2,2).*(1000./frameRateHz), (nanstd(latePiezoEventB2,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
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
subplot(2,1,2);hold on;
shadedErrorBar_CV(tt,nanmean((t1earlyPiezoV2B2),2),nanstd((t1earlyPiezoV2B2),[],2)./sqrt(sum(cell2mat(block2))),'lineProps','b');
hold on;
shadedErrorBar_CV(tt,nanmean((t1latePiezoV2B2),2),nanstd((t1latePiezoV2B2),[],2)./sqrt(sum(cell2mat(block2))),'lineProps','k');
vline(0,'r');
vline(767,'k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
title(['earliest 25% vs. latest 25% of motion of CS- trials: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_B2)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_B2)).*(1000./frameRateHz),2)) ' ms']);
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | motion bursts of CS- trials: early(blu) & late(blk)']);
if seshType ~= 2
%savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatencyFromPiezoTrace_B2.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatencyFromPiezoTrace_B2.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatencyFromPiezoTrace_B2.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatencyFromPiezoTrace_B2.pdf'],'pdf');
end

tempFig=setFigure; hold on;
subplot(2,1,1);hold on;
shadedErrorBar_CV(tt,nanmean(earlyPiezoEventRew,2).*(1000./frameRateHz), (nanstd(earlyPiezoEventRew,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
hold on;
shadedErrorBar_CV(tt,nanmean(latePiezoEventRew,2).*(1000./frameRateHz), (nanstd(latePiezoEventRew,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
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
subplot(2,1,2);hold on;
shadedErrorBar_CV(tt,nanmean(t1earlyPiezoV2Rew,2),nanstd(t1earlyPiezoV2Rew,[],2)./sqrt(sum(cell2mat(rew))),'lineProps','b');
hold on;
shadedErrorBar_CV(tt,nanmean(t1latePiezoV2Rew,2),nanstd(t1latePiezoV2Rew,[],2)./sqrt(sum(cell2mat(rew))),'lineProps','k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
vline(767,'k');
vline(0,'g');
title(['earliest 25% vs. latest 25% of motion of CS+ trials: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_rew)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_rew)).*(1000./frameRateHz),2)) ' ms']);
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | motion bursts of CS+ trials: early(blu) & late(blk)']);
if seshType ~= 2
%savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatencyFromPiezoTrace_Rew.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatencyFromPiezoTrace_Rew.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatencyFromPiezoTrace_Rew.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatencyFromPiezoTrace_Rew.pdf'],'pdf');
end

%% CS+
tempFig=setFigure; hold on;
subplot(2,1,1);hold on;
shadedErrorBar_CV(tt,nanmean(earlyPiezoEventRew,2).*(1000./frameRateHz), (nanstd(earlyPiezoEventRew,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
hold on;
shadedErrorBar_CV(tt,nanmean(latePiezoEventRew,2).*(1000./frameRateHz), (nanstd(latePiezoEventRew,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
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
subplot(2,1,2);hold on;
shadedErrorBar_CV(tt,nanmean(t3earlyPiezoRew,2),nanstd(t3earlyPiezoRew,[],2)./sqrt(sum(cell2mat(rew))),'lineProps','b');
hold on;
shadedErrorBar_CV(tt,nanmean(t3latePiezoRew,2),nanstd(t3latePiezoRew,[],2)./sqrt(sum(cell2mat(rew))),'lineProps','k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
vline(767,'k');
vline(0,'g');
title(['earliest 25% vs. latest 25% of motion of CS+ trials: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_rew)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_rew)).*(1000./frameRateHz),2)) ' ms']);
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | motion bursts of CS+ trials: early(blu) & late(blk)']);
if seshType ~= 2
%savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatency_Rew.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatency_Rew.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatency_Rew.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatency_Rew.pdf'],'pdf');
end

% 1SD-2SD over mean piezo voltage
tempFig=setFigure; hold on;
subplot(2,1,1);hold on;
shadedErrorBar_CV(tt,nanmean(earlyPiezoEventRew,2).*(1000./frameRateHz), (nanstd(earlyPiezoEventRew,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
hold on;
shadedErrorBar_CV(tt,nanmean(latePiezoEventRew,2).*(1000./frameRateHz), (nanstd(latePiezoEventRew,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
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
subplot(2,1,2);hold on;
shadedErrorBar_CV(tt,nanmean(t1earlyPiezoRew,2),nanstd(t1earlyPiezoRew,[],2)./sqrt(sum(cell2mat(rew))),'lineProps','b');
hold on;
shadedErrorBar_CV(tt,nanmean(t1latePiezoRew,2),nanstd(t1latePiezoRew,[],2)./sqrt(sum(cell2mat(rew))),'lineProps','k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
vline(767,'k');
vline(0,'g');
title(['earliest 25% vs. latest 25% of motion of CS+ trials: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_rew)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_rew)).*(1000./frameRateHz),2)) ' ms']);
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | motion bursts of CS+ trials: early(blu) & late(blk)']);
if seshType ~= 2
%savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatency1SD_Rew.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatency1SD_Rew.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatency1SD_Rew.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatency1SD_Rew.pdf'],'pdf');
end

%2-3SD over mean piezo voltage
tempFig=setFigure; hold on;
subplot(2,1,1);hold on;
shadedErrorBar_CV(tt,nanmean(earlyPiezoEventRew,2).*(1000./frameRateHz), (nanstd(earlyPiezoEventRew,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
hold on;
shadedErrorBar_CV(tt,nanmean(latePiezoEventRew,2).*(1000./frameRateHz), (nanstd(latePiezoEventRew,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
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
subplot(2,1,2);hold on;
shadedErrorBar_CV(tt,nanmean(t2earlyPiezoRew,2),nanstd(t2earlyPiezoRew,[],2)./sqrt(sum(cell2mat(rew))),'lineProps','b');
hold on;
shadedErrorBar_CV(tt,nanmean(t2latePiezoRew,2),nanstd(t2latePiezoRew,[],2)./sqrt(sum(cell2mat(rew))),'lineProps','k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
vline(767,'k');
vline(0,'g');
title(['earliest 25% vs. latest 25% of motion of CS+ trials: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_rew)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_rew)).*(1000./frameRateHz),2)) ' ms']);
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | motion bursts of CS+ trials: early(blu) & late(blk)']);
if seshType ~= 2
%savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatency2SD_Rew.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatency2SD_Rew.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatency2SD_Rew.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatency2SD_Rew.pdf'],'pdf');
end

%% late B2 vs. Rew
tempFig=setFigure; hold on;
subplot(2,1,1);hold on;
shadedErrorBar_CV(tt,nanmean(latePiezoEventRew,2).*(1000./frameRateHz), (nanstd(latePiezoEventRew,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
hold on;
shadedErrorBar_CV(tt,nanmean(latePiezoEventB2,2).*(1000./frameRateHz), (nanstd(latePiezoEventB2,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','r');
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
subplot(2,1,2);hold on;
shadedErrorBar_CV(tt,nanmean(t3latePiezoRew,2),nanstd(t3latePiezoRew,[],2)./sqrt(sum(cell2mat(rew))),'lineProps','k');
hold on;
shadedErrorBar_CV(tt,nanmean(t3latePiezoB2,2),nanstd(t3latePiezoB2,[],2)./sqrt(sum(cell2mat(block2))),'lineProps','r');
vline(0,'b');
vline(767,'k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
title(['latest 25% of motion of CS+ vs CS- trials: ' num2str(chop(nanmean(cell2mat(HL_piezo.late_rew)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.late_B2)).*(1000./frameRateHz),2)) ' ms']);
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | motion bursts of late trials: CS+(blk) & CS-(red)']);
if seshType ~= 2
%savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatency_Late.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatency_Late.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatency_Late.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatency_Late.pdf'],'pdf');
end

tempFig=setFigure; hold on;
subplot(2,1,1);hold on;
shadedErrorBar_CV(tt,nanmean(earlyPiezoEventRew,2).*(1000./frameRateHz), (nanstd(earlyPiezoEventRew,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
hold on;
shadedErrorBar_CV(tt,nanmean(earlyPiezoEventB2,2).*(1000./frameRateHz), (nanstd(earlyPiezoEventB2,[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','r');
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
subplot(2,1,2);hold on;
shadedErrorBar_CV(tt,nanmean(t3earlyPiezoRew,2),nanstd(t3earlyPiezoRew,[],2)./sqrt(sum(cell2mat(rew))),'lineProps','k');
hold on;
shadedErrorBar_CV(tt,nanmean(t3earlyPiezoB2,2),nanstd(t3earlyPiezoB2,[],2)./sqrt(sum(cell2mat(block2))),'lineProps','r');
vline(0,'b');
vline(767,'k');
xlabel('Time from cue');
ylabel('Piezo voltage');
ylim([0 inf]);
title(['earliest 25% of motion of CS+ vs. CS- trials: ' num2str(chop(nanmean(cell2mat(HL_piezo.early_rew)).*(1000./frameRateHz),2)) ' vs ' num2str(chop(nanmean(cell2mat(HL_piezo.early_B2)).*(1000./frameRateHz),2)) ' ms']);
sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | motion bursts of early trials: CS+(blk) & CS-(red)']);

if seshType ~= 2
%savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatency_Early.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatency_Early.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatency_Early.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatency_Early.pdf'],'pdf');
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

tempFig=setFigure; hold on;
subplot(2,1,1);hold on;
shadedErrorBar_CV(tt,nanmean(earlyPiezoEvent,2)*(1000./frameRateHz), (nanstd(earlyPiezoEvent,[],2)*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
hold on;
shadedErrorBar_CV(tt,nanmean(latePiezoEvent,2)*(1000./frameRateHz), (nanstd(latePiezoEvent,[],2)*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
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

subplot(2,1,2);hold on;
shadedErrorBar_CV(tt,nanmean(t3earlyPiezo,2), (nanstd(t3earlyPiezo,[],2))./sqrt(nTrials),'lineProps','b');
hold on;
shadedErrorBar_CV(tt,nanmean(t3latePiezo,2), (nanstd(t3latePiezo,[],2))./sqrt(nTrials),'lineProps','k');
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
%savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatency_earlyLate.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatency_earlyLate.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatency_earlyLate.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatency_earlyLate.pdf'],'pdf');
end
%
tempFig=setFigure; hold on;
subplot(2,1,1);hold on;
shadedErrorBar_CV(tt,nanmean(early25PiezoEvent,2)*(1000./frameRateHz), (nanstd(early25PiezoEvent,[],2)*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
hold on;
shadedErrorBar_CV(tt,nanmean(late25PiezoEvent,2)*(1000./frameRateHz), (nanstd(late25PiezoEvent,[],2)*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
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

subplot(2,1,2);hold on;
shadedErrorBar_CV(tt,nanmean(t3early25Piezo,2), (nanstd(t3early25Piezo,[],2))./sqrt(nTrials),'lineProps','b');
hold on;
shadedErrorBar_CV(tt,nanmean(t3late25Piezo,2), (nanstd(t3late25Piezo,[],2))./sqrt(nTrials),'lineProps','k');
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
%savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatency_25.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatency_25.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatency_25.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatency_25.pdf'],'pdf');
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
        plot(tt,nanmean(nanmean(groupAlign_piezo,2),1),'r',1,'LineWidth',1.0);
end
%}
tempFig=setFigure; hold on;
subplot(2,1,1);hold on;
shadedErrorBar_CV(tt,nanmean(early10PiezoEvent,2)*(1000./frameRateHz), (nanstd(early10PiezoEvent,[],2)*(1000./frameRateHz))./sqrt(nIC),'lineProps','b');
hold on;
shadedErrorBar_CV(tt,nanmean(late10PiezoEvent,2)*(1000./frameRateHz), (nanstd(late10PiezoEvent,[],2)*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
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

subplot(2,1,2);hold on;
shadedErrorBar_CV(tt,nanmean(t3early10Piezo,2), (nanstd(t3early10Piezo,[],2))./sqrt(nTrials),'lineProps','b');
hold on;
shadedErrorBar_CV(tt,nanmean(t3late10Piezo,2), (nanstd(t3late10Piezo,[],2))./sqrt(nTrials),'lineProps','k');
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
%savefig(fullfile(output_fn, '_cueAlignSpiking_byPiezoLatency_10.fig'));
saveas(tempFig, [output_fn  '_cueAlignSpiking_byPiezoLatency_10.pdf'],'pdf');
else
%savefig(fullfile(output_fn, [blockName '_cueAlignSpiking_byPiezoLatency_10.fig']));
saveas(tempFig, [output_fn blockName '_cueAlignSpiking_byPiezoLatency_10.pdf'],'pdf');
end
% response_frames piezo data
%}

%%%{    
s=1; s = s+1; 
        
        if unique(expt{1,1}.name == 'Block')
            indivFig = true;
        else
            indivFig = false;
        end
            
        b2GroupTC_acrossTrials=[]; rewGroupTC_acrossTrials=[]; targetAlignFB2_acrossTrials=[]; targetAlignFRew_acrossTrials=[]; eventsCSM_acrossTrials=[]; eventsCSP_acrossTrials=[]; 
        b2GroupTC_acrossNeurons=[]; rewGroupTC_acrossNeurons=[]; targetAlignFB2_acrossNeurons=[]; targetAlignFRew_acrossNeurons=[]; eventsCSM_acrossNeurons=[]; eventsCSP_acrossNeurons=[];
        ntrials=0;        
        for mouseIdx = 1:exptCount
               nTrial = (1:animalTrials(1,mouseIdx))+ntrials;  
b2GroupTC_acrossTrials = [b2GroupTC_acrossTrials nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2{1,mouseIdx})),3)];
rewGroupTC_acrossTrials = [rewGroupTC_acrossTrials nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rew{1,mouseIdx})),3)];

targetAlignFB2_acrossTrials = [targetAlignFB2_acrossTrials nanmean(groupAlign.tc{1,mouseIdx}(1:prewin_frames,:,logical(block2{1,mouseIdx})),3)]; % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.
targetAlignFRew_acrossTrials = [targetAlignFRew_acrossTrials nanmean(groupAlign.tc{1,mouseIdx}(1:prewin_frames,:,logical(rew{1,mouseIdx})),3)]; % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.

eventsCSM_acrossTrials = [eventsCSM_acrossTrials nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(block2{1,mouseIdx})),3)];
eventsCSP_acrossTrials = [eventsCSP_acrossTrials nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(rew{1,mouseIdx})),3)];

b2GroupTC_acrossNeurons = [b2GroupTC_acrossNeurons reshape(nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(block2{1,mouseIdx})),2),[150 size(groupAlign.tc{1,mouseIdx}(:,:,logical(block2{1,mouseIdx})),3)])];
rewGroupTC_acrossNeurons = [rewGroupTC_acrossNeurons reshape(nanmean(groupAlign.tc{1,mouseIdx}(:,:,logical(rew{1,mouseIdx})),2),[150 size(groupAlign.tc{1,mouseIdx}(:,:,logical(rew{1,mouseIdx})),3)])];

targetAlignFB2_acrossNeurons = [targetAlignFB2_acrossNeurons reshape(nanmean(groupAlign.tc{1,mouseIdx}(1:prewin_frames,:,logical(block2{1,mouseIdx})),2),[50 size(groupAlign.tc{1,mouseIdx}(1:prewin_frames,:,logical(block2{1,mouseIdx})),3)])]; % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.
targetAlignFRew_acrossNeurons = [targetAlignFRew_acrossNeurons reshape(nanmean(groupAlign.tc{1,mouseIdx}(1:prewin_frames,:,logical(rew{1,mouseIdx})),2),[50 size(groupAlign.tc{1,mouseIdx}(1:prewin_frames,:,logical(rew{1,mouseIdx})),3)])]; % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.

eventsCSM_acrossNeurons = [eventsCSM_acrossNeurons reshape(nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(block2{1,mouseIdx})),2),[150 size(groupAlign.events{1,mouseIdx}(:,:,logical(block2{1,mouseIdx})),3)])];
eventsCSP_acrossNeurons = [eventsCSP_acrossNeurons reshape(nanmean(groupAlign.events{1,mouseIdx}(:,:,logical(rew{1,mouseIdx})),2),[150 size(groupAlign.events{1,mouseIdx}(:,:,logical(rew{1,mouseIdx})),3)])];

ntrials=max(nTrial);
        end
        %
     for c = 1:sum(nsize)
        dFoFB2(:,c,:) = (b2GroupTC_acrossTrials(:,c,:)-targetAlignFB2_acrossTrials(1,c))./targetAlignFB2_acrossTrials(1,c);
        dFoFRew(:,c,:) = (rewGroupTC_acrossTrials(:,c,:)-targetAlignFRew_acrossTrials(1,c))./targetAlignFRew_acrossTrials(1,c);
     end
     
        if ~indivFig 
                tt=tt(1,36:96);
        tempFig=setFigure; hold on; n = 1;
        subplot(s,1,n);hold on;
        shadedErrorBar_CV(tt, nanmean(dFoFRew(36:96,:),2), nanstd(dFoFRew(36:96,:),[],2)./sqrt(nIC),'lineProps','k');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        %ylim([-3 5]);
        title('CS+')
        n = n+1;
%         if mworks.doBlock2
        subplot(s,1,n);hold on;
        shadedErrorBar_CV(tt, nanmean(dFoFB2(36:96,:),2), nanstd(dFoFB2(36:96,:),[],2)./sqrt(nIC),'lineProps','r');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        %ylim([-3 5]);
        title('CS-')
%         end
        sgtitle(['dFoF across animals | ' num2str(sum(nsize)) ' neurons from ' num2str(exptCount) ' animals']);
        if seshType ~= 2
        %savefig(fullfile(output_fn, '_cueAlign_dFoverF.fig'))
        saveas(tempFig, [output_fn  '_cueAlign_dFoverF.pdf'],'pdf');
        else
        %savefig(fullfile(output_fn, [blockName '_cueAlign_dFoverF.fig']));
        saveas(tempFig, [output_fn blockName '_cueAlign_dFoverF.pdf'],'pdf');
        end
        
        tempFig=setFigure; hold on; n=1;
        subplot(s,1,n);hold on;
        shadedErrorBar_CV(tt, nanmean(eventsCSP_acrossTrials(36:96,:),2).*(1000./frameRateHz), (nanstd(eventsCSP_acrossTrials(36:96,:),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','k');
        ylim([0 ymax]);
        vline(767,'k')
        vline(0,'k')       
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('CS+')
        n = n+1;
%         if mworks.doBlock2 (:,~gBlock2)
        subplot(s,1,n);hold on;
        shadedErrorBar_CV(tt, nanmean(eventsCSM_acrossTrials(36:96,:),2).*(1000./frameRateHz), (nanstd(eventsCSM_acrossTrials(36:96,:),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','r');
        ylim([0 ymax]);
        vline(767,'k')
        vline(0,'k')
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('CS-')
%         end
        sgtitle(['spike rate across animals | ' num2str(sum(nsize)) ' neurons from ' num2str(exptCount) ' animals'])
        if seshType ~= 2
        %savefig(fullfile(output_fn, '_cueAlign_events_Hz.fig'))
        saveas(tempFig, [output_fn  '_cueAlign_events_Hz.pdf'],'pdf');
        else
        %savefig(fullfile(output_fn, [blockName '_cueAlign_events_Hz.fig']));
        saveas(tempFig, [output_fn blockName '_cueAlign_events_Hz.pdf'],'pdf');
        end
        
                tempFig=setFigure; hold on; n=1;
        subplot(s,1,n);hold on;
        shadedErrorBar_CV(tt, nanmean(eventsCSP_acrossTrials(36:96,:),2).*(1000./frameRateHz), (nanstd(eventsCSP_acrossTrials(36:96,:),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','k');
        ylim([0 ymax]);
        errorbar((mean(firstLickStart(1,~gBlock2),2,'omitnan')-prewin_frames).*(1000./frameRateHz),max(max([nanmean(eventsCSP_acrossTrials,2).*(1000./frameRateHz) nanmean(eventsCSM_acrossTrials,2).*(1000./frameRateHz)]))+.5,std((firstLickStart(1,~gBlock2)),[],2,'omitnan')./sqrt(size(~gBlock2,2)).*(1000./frameRateHz),'horizontal','k');
%         errorbar(mean((firstPostRew_lickInd(1,~gBlock2)),2,'omitnan').*(1000./frameRateHz),max(max([nanmean(nanmean(eventsCSP,3),2).*(1000./frameRateHz) nanmean(nanmean(eventsCSM,3),2).*(1000./frameRateHz)]))+.5,std((firstPostRew_lickInd(1,~gBlock2)),[],2,'omitnan').*(1000./frameRateHz)./sqrt(size(~gBlock2,2)),'horizontal','k');
        vline(767,'k')
        vline(0,'k')       
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('CS+')
        n = n+1;
%         if mworks.doBlock2 (:,~gBlock2)
        subplot(s,1,n);hold on;
        shadedErrorBar_CV(tt, nanmean(eventsCSM_acrossTrials(36:96,:),2).*(1000./frameRateHz), (nanstd(eventsCSM_acrossTrials(36:96,:),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','r');
        ylim([0 ymax]);
        errorbar((mean(firstLickStart(1,gBlock2),2,'omitnan')-prewin_frames).*(1000./frameRateHz),max(max([nanmean(eventsCSP_acrossTrials,2).*(1000./frameRateHz) nanmean(eventsCSM_acrossTrials,2).*(1000./frameRateHz)]))+.5,std((firstLickStart(1,gBlock2)),[],2,'omitnan')./sqrt(size(gBlock2,2)).*(1000./frameRateHz),'horizontal','r');
%        errorbar(mean((firstPostRew_lickInd(1,gBlock2)),2,'omitnan').*(1000./frameRateHz),max(max([nanmean(nanmean(eventsCSP,3),2).*(1000./frameRateHz) nanmean(nanmean(eventsCSM,3),2).*(1000./frameRateHz)]))+.5,std((firstPostRew_lickInd(1,gBlock2)),[],2,'omitnan').*(1000./frameRateHz)./sqrt(size(gBlock2,2)),'horizontal','r');
        vline(767,'k')
        vline(0,'k')
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('CS-')
%         end
        sgtitle(['spike rate across animals | ' num2str(sum(nsize)) ' neurons from ' num2str(exptCount) ' animals'])
        if seshType ~= 2
        %savefig(fullfile(output_fn, '_cueAlign_events_lickBar_Hz.fig'))
        saveas(tempFig, [output_fn  '_cueAlign_events_lickBar_Hz.pdf'],'pdf');
        else
        %savefig(fullfile(output_fn, [blockName '_cueAlign_events_lickBar_Hz.fig']));
        saveas(tempFig, [output_fn blockName '_cueAlign_events_lickBar_Hz.pdf'],'pdf');
        end
        
                tempFig=setFigure; hold on; n=1;
        subplot(s,1,n);hold on;
        shadedErrorBar_CV(tt, nanmean(eventsCSP_acrossTrials(36:96,:),2).*(1000./frameRateHz), (nanstd(eventsCSP_acrossTrials(36:96,:),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','k');
        ylim([0 ymax]);
        errorbar((mean(firstPiezoStart(1,~gBlock2),2,'omitnan')-prewin_frames).*(1000./frameRateHz),max(max([nanmean(eventsCSP_acrossTrials,2).*(1000./frameRateHz) nanmean(eventsCSM_acrossTrials,2).*(1000./frameRateHz)]))+.5,std((firstPiezoStart(1,~gBlock2)),[],2,'omitnan')./sqrt(size(~gBlock2,2)).*(1000./frameRateHz),'horizontal','k');
%        errorbar(mean((firstPostRew_piezoInd(1,~gBlock2)),2,'omitnan').*(1000./frameRateHz),max(max([nanmean(nanmean(eventsCSP,3),2).*(1000./frameRateHz) nanmean(nanmean(eventsCSM,3),2).*(1000./frameRateHz)]))+.5,std((firstPostRew_piezoInd(1,~gBlock2)),[],2,'omitnan').*(1000./frameRateHz)./sqrt(size(~gBlock2,2)),'horizontal','k');
        vline(767,'k')
        vline(0,'k')       
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('CS+')
        n = n+1;
%         if mworks.doBlock2 (:,~gBlock2)
        subplot(s,1,n);hold on;
        shadedErrorBar_CV(tt, nanmean(eventsCSM_acrossTrials(36:96,:),2).*(1000./frameRateHz), (nanstd(eventsCSM_acrossTrials(36:96,:),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','r');
        errorbar((mean(firstPiezoStart(1,gBlock2),2,'omitnan')-prewin_frames).*(1000./frameRateHz),max(max([nanmean(eventsCSP_acrossTrials,2).*(1000./frameRateHz) nanmean(eventsCSM_acrossTrials,2).*(1000./frameRateHz)]))+.5,std((firstPiezoStart(1,gBlock2)),[],2,'omitnan')./sqrt(size(gBlock2,2)).*(1000./frameRateHz),'horizontal','r');
%        errorbar(mean((firstPostRew_piezoInd(1,gBlock2)),2,'omitnan').*(1000./frameRateHz),max(max([nanmean(nanmean(eventsCSP,3),2).*(1000./frameRateHz) nanmean(nanmean(eventsCSM,3),2).*(1000./frameRateHz)]))+.5,std((firstPostRew_piezoInd(1,gBlock2)),[],2,'omitnan').*(1000./frameRateHz)./sqrt(size(gBlock2,2)),'horizontal','r');
        ylim([0 ymax]);
        vline(767,'k')
        vline(0,'k')
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('CS-')
%         end
        sgtitle(['spike rate across animals | ' num2str(sum(nsize)) ' neurons from ' num2str(exptCount) ' animals'])
        if seshType ~= 2
        %savefig(fullfile(output_fn, '_cueAlign_events_piezoBar_Hz.fig'))
        saveas(tempFig, [output_fn  '_cueAlign_events_piezoBar_Hz.pdf'],'pdf');
        else
        %savefig(fullfile(output_fn, [blockName '_cueAlign_events_piezoBar_Hz.fig']));
        saveas(tempFig, [output_fn blockName '_cueAlign_events_piezoBar_Hz.pdf'],'pdf');
        end
        tt=gtt(1,:);
        %  
        elseif indivFig
        tt=tt(1,36:96);       
        tempFig=setFigure; hold on; %n = 1;
        %subplot(s,1,n) 
        if expt{1,1}.whichBlock(1,isesh) == 0
            shadedErrorBar_CV(tt, nanmean(dFoFB2(36:96,:),2), nanstd(dFoFB2(36:96,:),[],2)./sqrt(nIC),'lineProps','r');
            vline(0,'r')
            vline(767,'k')        
            xlabel('Time from cue')
            ylabel('Avg dF/F')
        elseif expt{1,1}.whichBlock(1,isesh) == 1
            shadedErrorBar_CV(tt, nanmean(dFoFRew(36:96,:),2), nanstd(dFoFRew(36:96,:),[],2)./sqrt(nIC),'lineProps','k');
            vline(0, 'g')
            vline(767,'k')        
            xlabel('Time from cue')
            ylabel('Avg dF/F')
        end
        sgtitle(['dF/F of block sessions: ' blockName])
        if seshType ~= 2
        %savefig(fullfile(output_fn, '_cueAlign_dFoverF.fig'))
        saveas(tempFig, [output_fn  '_cueAlign_dFoverF.pdf'],'pdf');
        else
        %savefig(fullfile(output_fn, [blockName '_cueAlign_dFoverF.fig']));
        saveas(tempFig, [output_fn blockName '_cueAlign_dFoverF.pdf'],'pdf');
        end
        
        tempFig=setFigure; hold on; %n=1;
        %subplot(s,1,n)
        if expt{1,1}.whichBlock(1,isesh) == 0
            shadedErrorBar_CV(tt, nanmean(eventsCSM_acrossTrials(36:96,:),2).*(1000./frameRateHz), (nanstd(eventsCSM_acrossTrials(36:96,:),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','r');
            ylim([0 ymax]);
            vline(0,'r')
            vline(767,'k')
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
        elseif expt{1,1}.whichBlock(1,isesh) == 1
            shadedErrorBar_CV(tt, nanmean(eventsCSP_acrossTrials(36:96,:),2).*(1000./frameRateHz), (nanstd(eventsCSP_acrossTrials(36:96,:),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','k');
            ylim([0 ymax]);
            vline(0,'g')
            vline(767,'k')
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
        end
        sgtitle(['spike rate for blocked sessions: ' blockName])
        if seshType ~= 2
        %savefig(fullfile(output_fn, '_cueAlign_events_Hz.fig'))
        saveas(tempFig, [output_fn  '_cueAlign_events_Hz.pdf'],'pdf');
        else
        %savefig(fullfile(output_fn, [blockName '_cueAlign_events_Hz.fig']));
        saveas(tempFig, [output_fn blockName '_cueAlign_events_Hz.pdf'],'pdf');
        end
        
          tempFig=setFigure; hold on; %n=1;
        %subplot(s,1,n)
        if expt{1,1}.whichBlock(1,isesh) == 0
            shadedErrorBar_CV(tt, nanmean(eventsCSM_acrossTrials(36:96,:),2).*(1000./frameRateHz), (nanstd(eventsCSM_acrossTrials(36:96,:),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','r');
            errorbar((mean(firstLickStart(1,gBlock2),2,'omitnan')-prewin_frames).*(1000./frameRateHz),max(max([nanmean(eventsCSP_acrossTrials,2).*(1000./frameRateHz) nanmean(eventsCSM_acrossTrials,2).*(1000./frameRateHz)]))+.5,std((firstLickStart(1,gBlock2)),[],2,'omitnan')./sqrt(size(gBlock2,2)).*(1000./frameRateHz),'horizontal','r');
            %errorbar(mean((firstPostRew_lickInd(1,gBlock2)),2,'omitnan').*(1000./frameRateHz),max(max([nanmean(nanmean(eventsCSP,3),2).*(1000./frameRateHz) nanmean(nanmean(eventsCSM,3),2).*(1000./frameRateHz)]))+.5,std((firstPostRew_lickInd(1,gBlock2)),[],2,'omitnan').*(1000./frameRateHz)./sqrt(size(gBlock2,2)),'horizontal','r');
            ylim([0 ymax]);
            vline(0,'r')
            vline(767,'k')
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
        elseif expt{1,1}.whichBlock(1,isesh) == 1
            shadedErrorBar_CV(tt, nanmean(eventsCSP_acrossTrials(36:96,:),2).*(1000./frameRateHz), (nanstd(eventsCSP_acrossTrials(36:96,:),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','k');
            errorbar((mean(firstLickStart(1,~gBlock2),2,'omitnan')-prewin_frames).*(1000./frameRateHz),max(max([nanmean(eventsCSP_acrossTrials,2).*(1000./frameRateHz) nanmean(eventsCSM_acrossTrials,2).*(1000./frameRateHz)]))+.5,std((firstLickStart(1,~gBlock2)),[],2,'omitnan')./sqrt(size(~gBlock2,2)).*(1000./frameRateHz),'horizontal','k');
            %errorbar(mean((firstPostRew_lickInd(1,~gBlock2)),2,'omitnan').*(1000./frameRateHz),max(max([nanmean(nanmean(eventsCSP,3),2).*(1000./frameRateHz) nanmean(nanmean(eventsCSM,3),2).*(1000./frameRateHz)]))+.5,std((firstPostRew_lickInd(1,~gBlock2)),[],2,'omitnan').*(1000./frameRateHz)./sqrt(size(gBlock2,2)),'horizontal','k');
            ylim([0 ymax]);
            vline(0,'g')
            vline(767,'k')
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
        end
        sgtitle(['spike rate for blocked sessions: ' blockName])
        if seshType ~= 2
        %savefig(fullfile(output_fn, '_cueAlign_events_lickBar_Hz.fig'))
        saveas(tempFig, [output_fn  '_cueAlign_events_lickBar_Hz.pdf'],'pdf');
        else
        %savefig(fullfile(output_fn, [blockName '_cueAlign_events_lickBar_Hz.fig']));
        saveas(tempFig, [output_fn blockName '_cueAlign_events_lickBar_Hz.pdf'],'pdf');
        end
        
          tempFig=setFigure; hold on; %n=1;
        %subplot(s,1,n)
        if expt{1,1}.whichBlock(1,isesh) == 0
            shadedErrorBar_CV(tt, nanmean(eventsCSM_acrossTrials(36:96,:),2).*(1000./frameRateHz), (nanstd(eventsCSM_acrossTrials(36:96,:),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','r');
            errorbar((mean(firstPiezoStart(1,gBlock2),2,'omitnan')-prewin_frames).*(1000./frameRateHz),max(max([nanmean(eventsCSP_acrossTrials,2).*(1000./frameRateHz) nanmean(eventsCSM_acrossTrials,2).*(1000./frameRateHz)]))+.5,std((firstPiezoStart(1,gBlock2)),[],2,'omitnan')./sqrt(size(gBlock2,2)).*(1000./frameRateHz),'horizontal','r');
            %errorbar(mean((firstPostRew_piezoInd(1,gBlock2)),2,'omitnan').*(1000./frameRateHz),max(max([nanmean(nanmean(eventsCSP,3),2).*(1000./frameRateHz) nanmean(nanmean(eventsCSM,3),2).*(1000./frameRateHz)]))+.5,std((firstPostRew_piezoInd(1,gBlock2)),[],2,'omitnan').*(1000./frameRateHz)./sqrt(size(gBlock2,2)),'horizontal','r');
            ylim([0 ymax]);
            vline(0,'r')
            vline(767,'k')
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
        elseif expt{1,1}.whichBlock(1,isesh) == 1
            shadedErrorBar_CV(tt, nanmean(eventsCSP_acrossTrials(36:96,:),2).*(1000./frameRateHz), (nanstd(eventsCSP_acrossTrials(36:96,:),[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','k');
            errorbar((mean(firstPiezoStart(1,~gBlock2),2,'omitnan')-prewin_frames).*(1000./frameRateHz),max(max([nanmean(eventsCSP_acrossTrials,2).*(1000./frameRateHz) nanmean(eventsCSM_acrossTrials,2).*(1000./frameRateHz)]))+.5,std((firstPiezoStart(1,~gBlock2)),[],2,'omitnan')./sqrt(size(~gBlock2,2)).*(1000./frameRateHz),'horizontal','k');
            %errorbar(mean((firstPostRew_piezoInd(1,~gBlock2)),2,'omitnan').*(1000./frameRateHz),max(max([nanmean(nanmean(eventsCSP,3),2).*(1000./frameRateHz) nanmean(nanmean(eventsCSM,3),2).*(1000./frameRateHz)]))+.5,std((firstPostRew_piezoInd(1,~gBlock2)),[],2,'omitnan').*(1000./frameRateHz)./sqrt(size(gBlock2,2)),'horizontal','k');
            ylim([0 ymax]);
            vline(0,'g')
            vline(767,'k')
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
        end
        sgtitle(['spike rate for blocked sessions: ' blockName])
        if seshType ~= 2
        %savefig(fullfile(output_fn, '_cueAlign_events_piezoBar_Hz.fig'))
        saveas(tempFig, [output_fn  '_cueAlign_events_piezoBar_Hz.pdf'],'pdf');
        else
        %savefig(fullfile(output_fn, [blockName '_cueAlign_events_piezoBar_Hz.fig']));
        saveas(tempFig, [output_fn blockName '_cueAlign_events_piezoBar_Hz.pdf'],'pdf');
        end
        end
       %
       eventsCSP = eventsCSP_acrossTrials;
       eventsCSM = eventsCSM_acrossTrials;
      tt=gtt(1,:);
        
        if unique(expt{1,1}.name == 'Inter')
         tempFig=setFigure; hold on; n=1;
        %subplot(s,1,n)
        plot(tt, nanmean(eventsCSP,2).*(1000./frameRateHz),'k');
        hold on;
        plot(tt, nanmean(eventsCSM,2).*(1000./frameRateHz),'r');
        ylim([0 5]);
        vline(0,'k')
        vline(767,'k')
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        %ylim([0 6]);
        legend('CS+','CS-');
        sgtitle(['spike rate across animals | ' num2str(sum(nsize)) ' neurons from ' num2str(exptCount) ' animals'])
        if seshType ~= 2
        %savefig(fullfile(output_fn, '_line_stackedSpike_events_Hz.fig'))
        saveas(tempFig, [output_fn  '_line_stackedSpike_events_Hz.pdf'],'pdf');
        else
        %savefig(fullfile(output_fn, [blockName '_line_stackedSpike_events_Hz.fig']));
        saveas(tempFig, [output_fn blockName '_line_stackedSpike_events_Hz.pdf'],'pdf');
        end
        
         tempFig=setFigure; hold on; n=1;
        %subplot(s,1,n)
        shadedErrorBar_CV(tt, nanmean(eventsCSP,2).*(1000./frameRateHz), (nanstd(eventsCSP,[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','k');
        hold on;
        shadedErrorBar_CV(tt, nanmean(eventsCSM,2).*(1000./frameRateHz), (nanstd(eventsCSM,[],2)./sqrt(nIC)).*(1000./frameRateHz),'lineProps','r');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([0 5]);
        %legend('CS+','CS-');  
        vline(767,'k')
        vline(0,'k')
        sgtitle(['stacked spike rate across animals | ' num2str(sum(nsize)) ' neurons from ' num2str(exptCount) ' animals'])
        hold off;
        if seshType ~= 2
        %savefig(fullfile(output_fn, '_error_stackedSpike_events_Hz.fig'))
        saveas(tempFig, [output_fn  '_error_stackedSpike_events_Hz.pdf'],'pdf');
        else
        %savefig(fullfile(output_fn, [blockName '_error_stackedSpike_events_Hz.fig']));
        saveas(tempFig, [output_fn blockName '_error_stackedSpike_events_Hz.pdf'],'pdf');
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
        
        fMLTrialCSP{1,1} = []; fMLTrialCSP{1,2} = []; fMLTrialCSP{1,3} = [];
        fMLTrialCSM{1,1} = []; fMLTrialCSM{1,2} = []; fMLTrialCSM{1,3} = [];
        for m = 1:3
            for i = 1:exptCount
            fMLTrialCSP{1,m} = [fMLTrialCSP{1,m} tempEventsRew{1,i}(:,:,indTemp{1,m})];
            fMLTrialCSM{1,m} = [fMLTrialCSM{1,m} tempEventsB2{1,i}(:,:,indTemp{1,m})];
            end
        end
        fMLTrialCSM{1,1} = nanmean(fMLTrialCSM{1,1},3);        fMLTrialCSP{1,1} = nanmean(fMLTrialCSP{1,1},3);
        fMLTrialCSM{1,2} = nanmean(fMLTrialCSM{1,2},3);        fMLTrialCSP{1,2} = nanmean(fMLTrialCSP{1,2},3);
        fMLTrialCSM{1,3} = nanmean(fMLTrialCSM{1,3},3);        fMLTrialCSP{1,3} = nanmean(fMLTrialCSP{1,3},3);

colorsRepCSP={[0 0 63/255],[0 0 128/255],[0 0 255/255];[0 0 0],[.4 .4 .4],[.7 .7 .7]};
% blue for low - black for high >> CS+
colorsRepCSM={[63/255 0 63/255],[128/255 0 128/255],[255/255 0 255/255];[63/255 0 0],[128/255 0 0],[255/255 0 0]};
% purple/magenta for low - red for high >> CS-

        tempFig=setFigure; hold on;
        n = floor(nTrials./3); nIC = sum(nsize);
        start = 1;
        for i = 1:3
            subplot(3,2,start); hold on;
            ind_rew_temp = intersect(cell2mat(rew),1+(i-1).*n:i*n);
            shadedErrorBar_CV(tt(1,36:96),nanmean(fMLTrialCSP{1,i}(36:96,:),2).*(1000./frameRateHz), (nanstd(fMLTrialCSP{1,i}(36:96,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps',{'color',colorsRepCSP{2,i}});
            hold on
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([0 ymax]);
            title(['Trials ' num2str(indTemp{1,i}(1,1)) ':' num2str(indTemp{1,i}(1,end))])
            vline(767,'k')
            vline(0,'k')
            if length(cell2mat(block2))>=10
                subplot(3,2,start+1); hold on;
                ind_block2_temp = intersect(cell2mat(block2),1+(i-1).*n:i*n);
                shadedErrorBar_CV(tt(1,36:96),nanmean(fMLTrialCSM{1,i}(36:96,:),2).*(1000./frameRateHz), (nanstd(fMLTrialCSM{1,i}(36:96,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps',{'color',colorsRepCSM{2,i}});
                ylim([0 ymax])
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
        %savefig(fullfile(output_fn, '_repsByTrial.fig'))
        saveas(tempFig, [output_fn  '_repsByTrial.pdf'],'pdf');
        else
        %savefig(fullfile(output_fn, [blockName '_repsByTrial.fig']));
        saveas(tempFig, [output_fn blockName '_repsByTrial.pdf'],'pdf');
        end
        
        
        
        elseif unique(expt{1,1}.name == 'Block')
            
z = [0 20 40]; for m = 1:3; indTemp{1,m} = (1:20)+z(1,m); end
        for mouseIdx=1:(exptCount)  
            tempEvents{1,mouseIdx} = nan(150,nsize(1,mouseIdx),60);
            if size(groupAlign.events{1,mouseIdx}(:,:,:),3) < 60
                tTr = size(groupAlign.events{1,mouseIdx}(:,:,:),3);
                tempEvents{1,mouseIdx}(:,:,1:tTr) = groupAlign.events{1,mouseIdx}(:,:,1:end);                
            elseif ~isnan(nanmean(nanmean(eventsCSM(:,:),2),1))
                tempEvents{1,mouseIdx}(:,:,:) = groupAlign.events{1,mouseIdx}(:,:,1:60);
            elseif isnan(nanmean(nanmean(eventsCSM(:,:),2),1))
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
            if ~isnan(nanmean(nanmean(eventsCSM(:,:),2),1))
        tempFig=setFigure; hold on;
        n = floor(nTrials./3);
        start = 1;
        for i = 1:3
            subplot(3,1,start);hold on;
            ind_rew_temp = intersect(cell2mat(block2),1+(i-1).*n:i*n);
            shadedErrorBar_CV(tt(1,36:96),nanmean(fMLTrial{1,i}(36:96,:),2).*(1000./frameRateHz), (nanstd(fMLTrial{1,i}(36:96,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','r');
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
        %savefig(fullfile(output_fn, '_minusrepsByTrial.fig'))
        saveas(tempFig, [output_fn  '_minusrepsByTrial.pdf'],'pdf');
        else
        %savefig(fullfile(output_fn, [blockName '_minusrepsByTrial.fig']));
        saveas(tempFig, [output_fn blockName '_minusrepsByTrial.pdf'],'pdf');
        end
            elseif isnan(nanmean(nanmean(eventsCSM(:,:),2),1))
        tempFig=setFigure; hold on;
        n = floor(nTrials./3);
        start = 1;
        for i = 1:3
            subplot(3,1,start);hold on;
            ind_rew_temp = intersect(cell2mat(rew),1+(i-1).*n:i*n);
            shadedErrorBar_CV(tt(1,36:96),nanmean(fMLTrial{1,i}(36:96,:),2).*(1000./frameRateHz), (nanstd(fMLTrial{1,i}(36:96,:),[],2).*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
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
        %savefig(fullfile(output_fn, '_plusrepsByTrial.fig'))
        saveas(tempFig, [output_fn  '_plusrepsByTrial.pdf'],'pdf');
        else
        %savefig(fullfile(output_fn, [blockName '_plusrepsByTrial.fig']));
        saveas(tempFig, [output_fn blockName '_plusrepsByTrial.pdf'],'pdf');
        end
            end
        end
        
sumPriorDend = 0; m=1;      
for i = 1:ndends
    idends = i - sumPriorDend; 
    block2PeakAmp(i,1) = max(nanmean(groupAlign.events{1,m}(:,idends,logical(block2{1,m}')),3).*(1000./frameRateHz));
    rewPeakAmp(i,1) = max(nanmean(groupAlign.events{1,m}(:,idends,~block2{1,m}'),3).*(1000./frameRateHz));
    
    currentDend = size(groupAlign.events{1,m},2);
    if i == currentDend+sumPriorDend
        sumPriorDend = sumPriorDend+currentDend;
        m=m+1;
    end    
end

colorsRepCSM={[63/255 0 0],[128/255 0 0],[255/255 0 0]};
colorsRepCSP={[0 0 0],[76.5/255 76.5/255 76.5/255],[153/255 153/255 153/255]};

if expType==2 && seshType==4
    %plot gui with 3 sections of activity from individual dendrites
    %averaged across ~20 trials
    q={'CSP','CSM'}; window={[75:87],[54:66]};
    for i=1:2
        eval(['currentArray = fMLTrial' q{i} ';']);
        eval(['colorMap=colorsRep' q{i} ';'])
           GUI_plotUniqueCell_sessionThirds(currentArray,colorMap,tt,frameRateHz,ymax,window{i});
    end

    %plot gui with average activity from individual dendrites across ~60
    %trials
    q={'CSP','CSM'}; window={[75:87],[54:66]}; colorMap={'k','r'};
    for i=1:2
        eval(['currentArray = events' q{i} ';']);
           GUI_plotUniqueCell(currentArray,colorMap{1,i},tt,frameRateHz,ymax,window{i});
    end
    %ID dendrites with larger dF/change in activity to CS+ vs. unexpected
    %reward
    for i = 1:2
         eval(['currentArray = events' q{i} ';']);
    avg_peak(i,1) = mean(max(currentArray(window{1,i},:),[],1),2,'omitnan');
    std_peak(i,1) = std(max(currentArray(window{1,i},:),[],1),[],2,'omitnan');
        ind_sigPosPeak(i,:) = zeros(1,size(eventsCSP,2)); ind_sigNegPeak(i,:) = zeros(1,size(eventsCSP,2));
    for ii = 1:size(eventsCSP,2)
        if max(eventsCSP(window{1,1},ii)) > avg_peak(i,1)+std_peak(i,1)
            ind_sigPosPeak(i,ii) = 1; %row 1 = CS+ | row 2 = CS-
        elseif max(eventsCSP(window{1,1},ii)) < avg_peak(i,1)-std_peak(i,1)
            ind_sigNegPeak(i,ii) = 1; %row 1 = CS+ | row 2 = CS-
        end
    end
    end
    
    sigPosCSP_sigPosCSM = zeros(1,size(eventsCSP,2)); sigPosCSP_sigNegCSM = zeros(1,size(eventsCSP,2)); sigPosCSP_nonSigCSM = zeros(1,size(eventsCSP,2)); sigNegCSP_sigPosCSM = zeros(1,size(eventsCSP,2)); sigNegCSP_sigNegCSM = zeros(1,size(eventsCSP,2)); sigNegCSP_nonSigCSM = zeros(1,size(eventsCSP,2)); nonSigCSP_sigPosCSM = zeros(1,size(eventsCSP,2)); nonSigCSP_sigNegCSM = zeros(1,size(eventsCSP,2)); nonSigCSP_nonSigCSM = zeros(1,size(eventsCSP,2));
                for i = 1:size(eventsCSP,2)
                    if ind_sigPosPeak(1,i)==1 && ind_sigPosPeak(2,i)==1
sigPosCSP_sigPosCSM(1,i) = 1;
                    elseif ind_sigPosPeak(1,i)==1 && ind_sigNegPeak(2,i)==1
sigPosCSP_sigNegCSM(1,i) = 1;
                    elseif ind_sigPosPeak(1,i)==1 && ind_sigNegPeak(2,i)==0 && ind_sigPosPeak(2,i)==0
sigPosCSP_nonSigCSM(1,i) = 1;
                    elseif ind_sigNegPeak(1,i)==1 && ind_sigPosPeak(2,i)==1
sigNegCSP_sigPosCSM(1,i) = 1;
                    elseif ind_sigNegPeak(1,i)==1 && ind_sigNegPeak(2,i)==1
sigNegCSP_sigNegCSM(1,i) = 1;
                    elseif ind_sigNegPeak(1,i)==1 && ind_sigNegPeak(2,i)==0 && ind_sigPosPeak(2,i)==0
sigNegCSP_nonSigCSM(1,i) = 1;
                    elseif ind_sigPosPeak(1,i)==0 && ind_sigNegPeak(1,i)==0 && ind_sigPosPeak(2,i)==1
nonSigCSP_sigPosCSM(1,i) = 1;
                    elseif ind_sigPosPeak(1,i)==0 && ind_sigNegPeak(1,i)==0 && ind_sigNegPeak(2,i)==1
nonSigCSP_sigNegCSM(1,i) = 1;
                    elseif ind_sigPosPeak(1,i)==0 && ind_sigNegPeak(1,i)==0 && ind_sigNegPeak(2,i)==0 && ind_sigPosPeak(2,i)==0
nonSigCSP_nonSigCSM(1,i) = 1;
                    end
                end

                arrayInQuestion = {'sigPosCSP_sigPosCSM', 'sigPosCSP_sigNegCSM', 'sigPosCSP_nonSigCSM', 'sigNegCSP_sigPosCSM', 'sigNegCSP_sigNegCSM', 'sigNegCSP_nonSigCSM', 'nonSigCSP_sigPosCSM', 'nonSigCSP_sigNegCSM', 'nonSigCSP_nonSigCSM'};
              tempFig=setFigure; hold on;
                for i = 1:size(arrayInQuestion,2)
                    eval(['thisCross=' arrayInQuestion{1,i} ';']);
                    subplot(3,3,i); hold on;
                    shadedErrorBar(tt(1,36:96),mean(eventsCSP(36:96,logical(thisCross)),2,'omitnan').*(1000./frameRateHz),std(eventsCSP(36:96,logical(thisCross)),[],2,'omitnan')./sqrt(sum(thisCross)).*(1000./frameRateHz),'k');
                    shadedErrorBar(tt(1,36:96),mean(eventsCSM(36:96,logical(thisCross)),2,'omitnan').*(1000./frameRateHz),std(eventsCSM(36:96,logical(thisCross)),[],2,'omitnan')./sqrt(sum(thisCross)).*(1000./frameRateHz),'r');
                    ylim([0 2*ymax]); vline(0,'k'); vline(767,'k'); title(arrayInQuestion{1,i});
                end
                    
                % plot in FOV
    iNeuron=1; arrayColors={[250/255 0 250/255],[175/255 0 175/255],[100/255 0 100/255],[250/255 12/5255 0],[175/255 98/255 0],[100/255 50/255 0],[0 250/255 250/255],[0 175/255 175/255],[0 100/255 100/255]};
for i = 1:size(expt,2)
    nNeuron = iNeuron:iNeuron+nsize(1,i)-1;
    
mouse = strtrim(expt{1,i}.mouse{1,seshID});
date = expt{1,i}.date{1,seshID};
run = expt{1,i}.run{1,seshID};    
img_fn2 = [date '_img' mouse '\getTC_' run '\'];
load([analysis_out, img_fn2, date '_img' mouse '_' run 'saveOutputs.mat']);
load([analysis_out, img_fn2, date '_img' mouse '_' run '_nPCA' num2str(nPCA) '_mu' num2str(mu) '_nIC' num2str(nIC) '_thresh' num2str(cluster_threshold) '_coor' num2str(threshold) '_mask3D_final.mat']);
tempFig=setFigure; hold on;
imshow([analysis_out,img_fn2 'AVG_' date '_img' mouse '_' run '_rgstr_tiff_every50_ref' num2str(ref) '_jpeg.jpg']); hold on;
for ii = 1:size(arrayInQuestion,2)
    eval(['thisCross=' arrayInQuestion{1,ii} ';']); cells=find(thisCross(1,nNeuron));
    for n = 1:size(cells,2)
    bound = cell2mat(bwboundaries(mask3D_final(:,:,cells(1,n))));
    plot(bound(:,2),bound(:,1),'.','color',arrayColors{1,ii});
    hold on;
    end
end  
    iNeuron=1+nNeuron(1,end);
end
          %output counts of cells with sig activity in each of the 3 sections and
    %count the # of sig activating dendrites


end
%close all;

 
%}
colorWheel ={[200/255 78/255 0/255],[153/255 51/255 153/255],[161/255 183/255 13/255],[255/255 217/255 96/255]};       
%% first event (pure and simple, no fuss : that's above)%{
    tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
    subplot(2,2,1); hold on;
    shadedErrorBar_CV(tt(1,prewin_frames-postLick_frames+1:prewin_frames+postLick_frames+postLick_frames), mean(mean(firstPiezoAlignEvents(:,:,~gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstPiezoAlignEvents(:,:,~gBlock2),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    shadedErrorBar_CV(tt(1,prewin_frames-postLick_frames+1:prewin_frames+postLick_frames+postLick_frames), mean(mean(firstPiezoAlignEvents(:,:,gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstPiezoAlignEvents(:,:,gBlock2),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','r');
    hold on;
    title('First motion');
    xlabel('Time from motion (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0+0.5 ymax+0.5]);
    subplot(2,2,3);hold on;
    shadedErrorBar_CV(tt(1,prewin_frames-postLick_frames+1:prewin_frames+postLick_frames+postLick_frames), mean(mean(firstPiezoAlign(:,~gBlock2),3,'omitnan'),2,'omitnan'), (std(mean(firstPiezoAlign(:,~gBlock2),3,'omitnan'),[],2,'omitnan'))./sqrt(unique(sum(~isnan(mean(firstPiezoAlign,1,'omitnan'))))),'lineProps','k');
    shadedErrorBar_CV(tt(1,prewin_frames-postLick_frames+1:prewin_frames+postLick_frames+postLick_frames), mean(mean(firstPiezoAlign(:,gBlock2),3,'omitnan'),2,'omitnan'), (std(mean(firstPiezoAlign(:,gBlock2),3,'omitnan'),[],2,'omitnan'))./sqrt(unique(sum(~isnan(mean(firstPiezoAlign,1,'omitnan'))))),'lineProps','r');
    xlabel('Time from motion (ms)');
    ylabel('Piezo voltage (mV)');
    ylim([0 1]);
    title([num2str(sum(~isnan(firstPiezoStart(:,~gBlock2)))) '+ ' num2str(sum(~isnan(firstPiezoStart(:,gBlock2)))) '- trials with motion']);
    
    subplot(2,2,2); hold on;
    shadedErrorBar_CV(tt(1,prewin_frames-postLick_frames+1:prewin_frames+postLick_frames+postLick_frames), mean(mean(firstLickAlignEvents(:,:,~gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,~gBlock2),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    shadedErrorBar_CV(tt(1,prewin_frames-postLick_frames+1:prewin_frames+postLick_frames+postLick_frames), mean(mean(firstLickAlignEvents(:,:,gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,gBlock2),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','r');
    hold on;
    title('First lick');
    xlabel('Time from lick (ms)');
    ylabel('Spike rate (Hz)');
    ylim([0+0.5 ymax+0.5]);
    subplot(2,2,4);hold on;
    shadedErrorBar_CV(tt(1,prewin_frames-postLick_frames+1:prewin_frames+postLick_frames+postLick_frames), mean(mean(firstLickAlign(:,~gBlock2),3,'omitnan'),2,'omitnan'), (std(mean(firstLickAlign(:,~gBlock2),3,'omitnan'),[],2,'omitnan'))./sqrt(unique(sum(~isnan(mean(firstLickAlign,1,'omitnan'))))),'lineProps','k');
    shadedErrorBar_CV(tt(1,prewin_frames-postLick_frames+1:prewin_frames+postLick_frames+postLick_frames), mean(mean(firstLickAlign(:,gBlock2),3,'omitnan'),2,'omitnan'), (std(mean(firstLickAlign(:,gBlock2),3,'omitnan'),[],2,'omitnan'))./sqrt(unique(sum(~isnan(mean(firstLickAlign,1,'omitnan'))))),'lineProps','r');
    xlabel('Time from lick (ms)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title([num2str(sum(~isnan(firstLickStart(:,~gBlock2)))) '+ ' num2str(sum(~isnan(firstLickStart(:,gBlock2)))) '- trials with lick']);
    
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Cspk relative to lick']);
    hold off;
    if seshType ~= 2
    %savefig(fullfile(output_fn, '_firstLick.fig'));
    saveas(tempFig, [output_fn  '_firstEvent.pdf'],'pdf');
    else
    %savefig(fullfile(output_fn, [blockName '_firstLick.fig']));
    saveas(tempFig, [output_fn blockName '_firstEvent.pdf'],'pdf');
    end
%}    
