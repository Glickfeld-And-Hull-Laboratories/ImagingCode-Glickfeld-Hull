%% group-wise analysis that outputs blocked figures for presentation

%% data org | load experiment data 
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
%% data org | load mouse data for given experiment 
thresholdDeco = input('Which threshold do you wish to plot? ("-#" for decon | "#.#" for first der.): ');
analysis_out = 'A:\home\carlo\RC\analysis\2P\';
output_fn = ([analysis_out 'groupwiseAlign' num2str(thresholdDeco) '_' sesh '_' sesh2 '\']);
if seshType ~= 2
load(fullfile(output_fn, '_groupwiseDataset.mat'));
elseif seshType == 2
blockName = input('which day to analyze? Call the CS+ session and it will load the companion CS- session. [type "plus + First/Second" or "allPlus"]: ','s');
    if strcmp(blockName,'plusFirst')
        plusData = load(fullfile(output_fn, [blockName '_groupwiseDataset.mat']));
        minusData = load(fullfile(output_fn, ['minusSecond_groupwiseDataset.mat']));
    elseif strcmp(blockName,'plusSecond')
        plusData = load(fullfile(output_fn, [blockName '_groupwiseDataset.mat']));
        minusData = load(fullfile(output_fn, ['minusFirst_groupwiseDataset.mat']));
    elseif strcmp(blockName,'allPlus')
        plusData = load(fullfile(output_fn, [blockName '_groupwiseDataset.mat']));
        minusData = load(fullfile(output_fn, ['allMinus_groupwiseDataset.mat']));
    end
end
analysis_out = 'A:\home\carlo\RC\analysis\2P\';
output_fn = ([analysis_out 'groupwiseAlign' num2str(thresholdDeco) '_' sesh '_' sesh2 '\']);

if ~exist(fullfile([output_fn 'stats\'])) 
 mkdir(fullfile([output_fn 'stats\']))
end
%% data org | fixed variables 
frameRateHz = 30;
prewin_frames = round(1500./frameRateHz);
postwin_frames = round(3000./frameRateHz);
response_delay = 1; %provides a window for the average delay between cue and cue-response 
tt = (-prewin_frames:postwin_frames-1).*(1000./frameRateHz);   
rewDelay_frames=23;
response_frames=11; %333 milliseconds
bx_response_frames=11; %333 milliseconds

respReq_cueAlign = 3; respReq_lickAlign = 3; respReq_piezoAlign = 3; 
suppReq_cueAlign = 2; suppReq_lickAlign = 1; suppReq_piezoAlign = 1;
chngReq_cueAlign = 3; chngReq_lickAlign = 2; chngReq_piezoAlign = 2;
% !! for lick and piezo align - 2 is actually u*1.5o
statsWindow = 1; %this determines the +/- window around the peak

baselineRange=prewin_frames-response_frames+response_delay:prewin_frames;
postCueRange=prewin_frames+response_delay:prewin_frames+response_frames;
postRewRange=prewin_frames+rewDelay_frames+response_delay:prewin_frames+rewDelay_frames+response_frames;

%% data org | cue-aligned 
for c = 1:size(plusData.groupAlign.goodcells,2)
    nsize(1,c) = size(plusData.groupAlign.tc{1,c},2);
    nsize(2,c) = size(minusData.groupAlign.tc{1,c},2);
end

CSMtimecourse_acrossTrials=[]; CSPtimecourse_acrossTrials=[]; CSMbaseline_acrossTrials=[]; CSPbaseline_acrossTrials=[]; CSMevents_acrossTrials=[]; CSPevents_acrossTrials=[]; 
CSMtimecourse_acrossNeurons=[]; CSPtimecourse_acrossNeurons=[]; CSMbaseline_acrossNeurons=[]; CSPbaseline_acrossNeurons=[]; eventsCSM_acrossNeurons=[]; eventsCSP_acrossNeurons=[];
first_lickAlignCSP_acrossTrials=[]; firstPostCue_lickAlignCSP_acrossTrials=[]; firstPostRew_lickAlignCSP_acrossTrials=[]; first_lickAlignCSM_acrossTrials=[]; firstPostCue_lickAlignCSM_acrossTrials=[]; firstPostRew_lickAlignCSM_acrossTrials=[];
first_piezoAlignCSP_acrossTrials=[]; firstPostCue_piezoAlignCSP_acrossTrials=[]; firstPostRew_piezoAlignCSP_acrossTrials=[]; first_piezoAlignCSM_acrossTrials=[]; firstPostCue_piezoAlignCSM_acrossTrials=[]; firstPostRew_piezoAlignCSM_acrossTrials=[];
CSMdFoF_acrossTrials=[]; CSPdFoF_acrossTrials=[]; ntrials=0;  
first_lickAlignCSPevents_acrossTrials=[]; first_lickAlignCSMevents_acrossTrials=[]; postCue_lickAlignCSPevents_acrossTrials=[]; postCue_lickAlignCSMevents_acrossTrials=[]; postRew_lickAlignCSPevents_acrossTrials=[]; postRew_lickAlignCSMevents_acrossTrials=[];
first_piezoAlignCSPevents_acrossTrials=[]; first_piezoAlignCSMevents_acrossTrials=[]; postCue_piezoAlignCSPevents_acrossTrials=[]; postCue_piezoAlignCSMevents_acrossTrials=[]; postRew_piezoAlignCSPevents_acrossTrials=[]; postRew_piezoAlignCSMevents_acrossTrials=[];
first_lickAlignDFoFCSP_acrossTrials=[]; first_lickAlignDFoFCSM_acrossTrials=[]; postCue_lickAlignDFoFCSP_acrossTrials=[]; postCue_lickAlignDFoFCSM_acrossTrials=[]; postRew_lickAlignDFoFCSP_acrossTrials=[]; postRew_lickAlignDFoFCSM_acrossTrials=[];
first_piezoAlignDFoFCSP_acrossTrials=[]; first_piezoAlignDFoFCSM_acrossTrials=[]; postCue_piezoAlignDFoFCSP_acrossTrials=[]; postCue_piezoAlignDFoFCSM_acrossTrials=[]; postRew_piezoAlignDFoFCSP_acrossTrials=[]; postRew_piezoAlignDFoFCSM_acrossTrials=[];

%for CS+ block
for mouse = 1:plusData.exptCount
nTrial = (1:plusData.animalTrials(1,mouse))+ntrials;  

for c = 1:nsize(1,mouse)
CSPdFoF_acrossTrials = [CSPdFoF_acrossTrials mean(plusData.groupAlign.dFoF{1,mouse}(:,c,logical(plusData.rew{1,mouse})),3,'omitnan')];%-mean(groupAlign.tc{1,mouseIdx}(1:prewin_frames,c,logical(rew{1,mouseIdx})),1,'omitnan'),3,'omitnan')];
CSPevents_acrossTrials = [CSPevents_acrossTrials mean(plusData.groupAlign.events{1,mouse}(:,c,logical(plusData.rew{1,mouse})),3,'omitnan')];
end
CSPtimecourse_acrossTrials = [CSPtimecourse_acrossTrials nanmean(plusData.groupAlign.tc{1,mouse}(:,:,logical(plusData.rew{1,mouse})),3)];
CSPbaseline_acrossTrials = [CSPbaseline_acrossTrials nanmean(plusData.groupAlign.tc{1,mouse}(1:prewin_frames,:,logical(plusData.rew{1,mouse})),3)]; % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.

CSPtimecourse_acrossNeurons = [CSPtimecourse_acrossNeurons reshape(nanmean(plusData.groupAlign.tc{1,mouse}(:,:,logical(plusData.rew{1,mouse})),2),[150 size(plusData.groupAlign.tc{1,mouse}(:,:,logical(plusData.rew{1,mouse})),3)])];
CSPbaseline_acrossNeurons = [CSPbaseline_acrossNeurons reshape(nanmean(plusData.groupAlign.tc{1,mouse}(1:prewin_frames,:,logical(plusData.rew{1,mouse})),2),[50 size(plusData.groupAlign.tc{1,mouse}(1:prewin_frames,:,logical(plusData.rew{1,mouse})),3)])]; % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.
eventsCSP_acrossNeurons = [eventsCSP_acrossNeurons reshape(nanmean(plusData.groupAlign.events{1,mouse}(:,:,logical(plusData.rew{1,mouse})),2),[150 size(plusData.groupAlign.events{1,mouse}(:,:,logical(plusData.rew{1,mouse})),3)])];

ntrials=max(nTrial);
end

first_lickAlignCSP_acrossTrials = [first_lickAlignCSP_acrossTrials sum(~isnan(plusData.firstLickStart(:,~plusData.gBlock2)))];
firstPostCue_lickAlignCSP_acrossTrials = [firstPostCue_lickAlignCSP_acrossTrials sum(~isnan(plusData.firstPreRew_lickStart(:,~plusData.gBlock2)))];
firstPostRew_lickAlignCSP_acrossTrials = [firstPostRew_lickAlignCSP_acrossTrials sum(~isnan(plusData.firstPostRew_lickStart(:,~plusData.gBlock2)))];

first_piezoAlignCSP_acrossTrials = [first_piezoAlignCSP_acrossTrials sum(~isnan(plusData.firstPiezoStart(:,~plusData.gBlock2)))];
firstPostCue_piezoAlignCSP_acrossTrials = [firstPostCue_piezoAlignCSP_acrossTrials sum(~isnan(plusData.firstPreRew_piezoStart(:,~plusData.gBlock2)))];
firstPostRew_piezoAlignCSP_acrossTrials = [firstPostRew_piezoAlignCSP_acrossTrials sum(~isnan(plusData.firstPostRew_piezoStart(:,~plusData.gBlock2)))];

first_lickAlignCSPevents_acrossTrials = [first_lickAlignCSPevents_acrossTrials mean(plusData.firstLickAlignEvents(:,:,~plusData.gBlock2),3,'omitnan')]; 
postCue_lickAlignCSPevents_acrossTrials = [postCue_lickAlignCSPevents_acrossTrials mean(plusData.firstPreRew_lickAlignEvents(:,:,~plusData.gBlock2),3,'omitnan')]; 
postRew_lickAlignCSPevents_acrossTrials = [postRew_lickAlignCSPevents_acrossTrials mean(plusData.firstPostRew_lickAlignEvents(:,:,~plusData.gBlock2),3,'omitnan')]; 

first_piezoAlignCSPevents_acrossTrials = [first_piezoAlignCSPevents_acrossTrials mean(plusData.firstPiezoAlignEvents(:,:,~plusData.gBlock2),3,'omitnan')]; 
postCue_piezoAlignCSPevents_acrossTrials = [postCue_piezoAlignCSPevents_acrossTrials mean(plusData.firstPreRew_piezoAlignEvents(:,:,~plusData.gBlock2),3,'omitnan')]; 
postRew_piezoAlignCSPevents_acrossTrials = [postRew_piezoAlignCSPevents_acrossTrials mean(plusData.firstPostRew_piezoAlignEvents(:,:,~plusData.gBlock2),3,'omitnan')]; 

first_lickAlignDFoFCSP_acrossTrials = [first_lickAlignDFoFCSP_acrossTrials mean(plusData.firstLickAlignDFoF(:,:,~plusData.gBlock2),3,'omitnan')]; 
postCue_lickAlignDFoFCSP_acrossTrials = [postCue_lickAlignDFoFCSP_acrossTrials mean(plusData.firstPreRew_lickAlignDFoF(:,:,~plusData.gBlock2),3,'omitnan')]; 
postRew_lickAlignDFoFCSP_acrossTrials = [postRew_lickAlignDFoFCSP_acrossTrials mean(plusData.firstPostRew_lickAlignDFoF(:,:,~plusData.gBlock2),3,'omitnan')]; 

first_piezoAlignDFoFCSP_acrossTrials = [first_piezoAlignDFoFCSP_acrossTrials mean(plusData.firstPiezoAlignDFoF(:,:,~plusData.gBlock2),3,'omitnan')]; 
postCue_piezoAlignDFoFCSP_acrossTrials = [postCue_piezoAlignDFoFCSP_acrossTrials mean(plusData.firstPreRew_piezoAlignDFoF(:,:,~plusData.gBlock2),3,'omitnan')]; 
postRew_piezoAlignDFoFCSP_acrossTrials = [postRew_piezoAlignDFoFCSP_acrossTrials mean(plusData.firstPostRew_piezoAlignDFoF(:,:,~plusData.gBlock2),3,'omitnan')]; 

%for CS- block
for mouse = 1:minusData.exptCount
nTrial = (1:minusData.animalTrials(1,mouse))+ntrials;  

for c = 1:nsize(2,mouse)
CSMdFoF_acrossTrials = [CSMdFoF_acrossTrials mean(minusData.groupAlign.dFoF{1,mouse}(:,c,logical(minusData.block2{1,mouse})),3,'omitnan')];%-mean(groupAlign.tc{1,mouseIdx}(1:prewin_frames,c,logical(block2{1,mouseIdx})),1,'omitnan'),3,'omitnan')];
CSMevents_acrossTrials = [CSMevents_acrossTrials mean(minusData.groupAlign.events{1,mouse}(:,c,logical(minusData.block2{1,mouse})),3,'omitnan')];
end

CSMtimecourse_acrossTrials = [CSMtimecourse_acrossTrials nanmean(minusData.groupAlign.tc{1,mouse}(:,:,logical(minusData.block2{1,mouse})),3)];
CSMbaseline_acrossTrials = [CSMbaseline_acrossTrials nanmean(minusData.groupAlign.tc{1,mouse}(1:prewin_frames,:,logical(minusData.block2{1,mouse})),3)]; % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.

CSMtimecourse_acrossNeurons = [CSMtimecourse_acrossNeurons reshape(nanmean(minusData.groupAlign.tc{1,mouse}(:,:,logical(minusData.block2{1,mouse})),2),[150 size(minusData.groupAlign.tc{1,mouse}(:,:,logical(minusData.block2{1,mouse})),3)])];
CSMbaseline_acrossNeurons = [CSMbaseline_acrossNeurons reshape(nanmean(minusData.groupAlign.tc{1,mouse}(1:prewin_frames,:,logical(minusData.block2{1,mouse})),2),[50 size(minusData.groupAlign.tc{1,mouse}(1:prewin_frames,:,logical(minusData.block2{1,mouse})),3)])]; % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.
eventsCSM_acrossNeurons = [eventsCSM_acrossNeurons reshape(nanmean(minusData.groupAlign.events{1,mouse}(:,:,logical(minusData.block2{1,mouse})),2),[150 size(minusData.groupAlign.events{1,mouse}(:,:,logical(minusData.block2{1,mouse})),3)])];

ntrials=max(nTrial);
end

first_lickAlignCSM_acrossTrials = [first_lickAlignCSM_acrossTrials sum(~isnan(minusData.firstLickStart(:,minusData.gBlock2)))];
firstPostCue_lickAlignCSM_acrossTrials = [firstPostCue_lickAlignCSM_acrossTrials sum(~isnan(minusData.firstPreRew_lickStart(:,minusData.gBlock2)))];
firstPostRew_lickAlignCSM_acrossTrials = [firstPostRew_lickAlignCSM_acrossTrials sum(~isnan(minusData.firstPostRew_lickStart(:,minusData.gBlock2)))];

first_piezoAlignCSM_acrossTrials = [first_piezoAlignCSM_acrossTrials sum(~isnan(minusData.firstPiezoStart(:,minusData.gBlock2)))];
firstPostCue_piezoAlignCSM_acrossTrials = [firstPostCue_piezoAlignCSM_acrossTrials sum(~isnan(minusData.firstPreRew_piezoStart(:,minusData.gBlock2)))];
firstPostRew_piezoAlignCSM_acrossTrials = [firstPostRew_piezoAlignCSM_acrossTrials sum(~isnan(minusData.firstPostRew_piezoStart(:,minusData.gBlock2)))];

first_lickAlignCSMevents_acrossTrials = [first_lickAlignCSMevents_acrossTrials mean(minusData.firstLickAlignEvents(:,:,minusData.gBlock2),3,'omitnan')];
postCue_lickAlignCSMevents_acrossTrials = [postCue_lickAlignCSMevents_acrossTrials mean(minusData.firstPreRew_lickAlignEvents(:,:,minusData.gBlock2),3,'omitnan')];
postRew_lickAlignCSMevents_acrossTrials = [postRew_lickAlignCSMevents_acrossTrials mean(minusData.firstPostRew_lickAlignEvents(:,:,minusData.gBlock2),3,'omitnan')];

first_piezoAlignCSMevents_acrossTrials = [first_piezoAlignCSMevents_acrossTrials mean(minusData.firstPiezoAlignEvents(:,:,minusData.gBlock2),3,'omitnan')];
postCue_piezoAlignCSMevents_acrossTrials = [postCue_piezoAlignCSMevents_acrossTrials mean(minusData.firstPreRew_piezoAlignEvents(:,:,minusData.gBlock2),3,'omitnan')];
postRew_piezoAlignCSMevents_acrossTrials = [postRew_piezoAlignCSMevents_acrossTrials mean(minusData.firstPostRew_piezoAlignEvents(:,:,minusData.gBlock2),3,'omitnan')];

first_lickAlignDFoFCSM_acrossTrials = [first_lickAlignDFoFCSM_acrossTrials mean(minusData.firstLickAlignDFoF(:,:,minusData.gBlock2),3,'omitnan')];
postCue_lickAlignDFoFCSM_acrossTrials = [postCue_lickAlignDFoFCSM_acrossTrials mean(minusData.firstPreRew_lickAlignDFoF(:,:,minusData.gBlock2),3,'omitnan')];
postRew_lickAlignDFoFCSM_acrossTrials = [postRew_lickAlignDFoFCSM_acrossTrials mean(minusData.firstPostRew_lickAlignDFoF(:,:,minusData.gBlock2),3,'omitnan')];

postCue_piezoAlignDFoFCSM_acrossTrials = [postCue_piezoAlignDFoFCSM_acrossTrials mean(minusData.firstPreRew_piezoAlignDFoF(:,:,minusData.gBlock2),3,'omitnan')];
first_piezoAlignDFoFCSM_acrossTrials = [first_piezoAlignDFoFCSM_acrossTrials mean(minusData.firstPiezoAlignDFoF(:,:,minusData.gBlock2),3,'omitnan')];
postRew_piezoAlignDFoFCSM_acrossTrials = [postRew_piezoAlignDFoFCSM_acrossTrials mean(minusData.firstPostRew_piezoAlignDFoF(:,:,minusData.gBlock2),3,'omitnan')];
%
colorWheel ={[200/255 78/255 0/255],[153/255 51/255 153/255],[161/255 183/255 13/255],[255/255 217/255 96/255]};
%
%% data org | segmentation of neurons into response groups by means AND change (+/-) 
eventLabel={'postCue','postRew'}; cueLabel={'CSP','CSM'}; 
  %CS+
  for event=1:size(eventLabel,2)
      eval(['cell.cueAlign.response_' eventLabel{1,event} 'CSP=[];']);
      eval(['cell.cueAlign.suppress_' eventLabel{1,event} 'CSP=[];']);
    for c = 1:sum(nsize(1,:))
%
 eval(['windowVar = CSPdFoF_acrossTrials(' eventLabel{1,event} 'Range,c);']);
 baselineVar = CSPdFoF_acrossTrials(baselineRange,c);
 pos.Diff.cueAlign{c,1} = ([0; diff(CSPdFoF_acrossTrials(:,c))]);

pos.Limit.cueAlign{c,1}.one = (mean(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),1,'omitnan')+(std(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),[],1,'omitnan').*1)); 
pos.Change.cueAlign{c,1}.one = find(pos.Diff.cueAlign{c,1}>pos.Limit.cueAlign{c,1}.one); 
pos.Limit.cueAlign{c,1}.two = (mean(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),1,'omitnan')+(std(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),[],1,'omitnan').*2)); 
pos.Change.cueAlign{c,1}.two = find(pos.Diff.cueAlign{c,1}>pos.Limit.cueAlign{c,1}.two); 
pos.Limit.cueAlign{c,1}.three = (mean(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),1,'omitnan')+(std(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),[],1,'omitnan').*3)); 
pos.Change.cueAlign{c,1}.three = find(pos.Diff.cueAlign{c,1}>pos.Limit.cueAlign{c,1}.three);

eval(['frameIdx.cueAlign.' eventLabel{1,event} 'CSP_posChangeDiff{1,c} = intersect(' eventLabel{1,event} 'Range, pos.Change.cueAlign{c,1}.one);']); 
eval(['frameIdx.cueAlign.' eventLabel{1,event} 'CSP_posChangeDiff{2,c} = intersect(' eventLabel{1,event} 'Range, pos.Change.cueAlign{c,1}.two);']); 
eval(['frameIdx.cueAlign.' eventLabel{1,event} 'CSP_posChangeDiff{3,c} = intersect(' eventLabel{1,event} 'Range, pos.Change.cueAlign{c,1}.three);']);

neg.Diff.cueAlign{c,1} = ([0; diff(CSPdFoF_acrossTrials(:,c))]);

neg.Limit.cueAlign{c,1}.one = (mean(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),1,'omitnan')-(std(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),[],1,'omitnan').*1)); 
neg.Change.cueAlign{c,1}.one = find(neg.Diff.cueAlign{c,1}<neg.Limit.cueAlign{c,1}.one); 
neg.Limit.cueAlign{c,1}.two = (mean(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),1,'omitnan')-(std(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),[],1,'omitnan').*2)); 
neg.Change.cueAlign{c,1}.two = find(neg.Diff.cueAlign{c,1}<neg.Limit.cueAlign{c,1}.two); 
neg.Limit.cueAlign{c,1}.three = (mean(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),1,'omitnan')-(std(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),[],1,'omitnan').*3)); 
neg.Change.cueAlign{c,1}.three = find(neg.Diff.cueAlign{c,1}<neg.Limit.cueAlign{c,1}.three);

eval(['frameIdx.cueAlign.' eventLabel{1,event} 'CSP_negChangeDiff{1,c} = intersect(' eventLabel{1,event} 'Range, neg.Change.cueAlign{c,1}.one);']); 
eval(['frameIdx.cueAlign.' eventLabel{1,event} 'CSP_negChangeDiff{2,c} = intersect(' eventLabel{1,event} 'Range, neg.Change.cueAlign{c,1}.two);']); 
eval(['frameIdx.cueAlign.' eventLabel{1,event} 'CSP_negChangeDiff{3,c} = intersect(' eventLabel{1,event} 'Range, neg.Change.cueAlign{c,1}.three);']);
%    

    tempResp_DiffIdx = diff(intersect(1:response_frames, find(windowVar>(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*respReq_cueAlign))))';
    tempResp_FrameIdx = intersect(1:response_frames, find(windowVar>(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*respReq_cueAlign)));
    if length(tempResp_DiffIdx)
      if strfind(tempResp_DiffIdx, [1 1])
        if eval(['frameIdx.cueAlign.' eventLabel{1,event} 'CSP_posChangeDiff{chngReq_cueAlign,c}'])
           eval(['cell.cueAlign.response_' eventLabel{1,event} 'CSP=[cell.cueAlign.response_' eventLabel{1,event} 'CSP c];']);
        end
      end
    end
    tempSupp_DiffIdx = diff(intersect(1:response_frames, find(windowVar<(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*suppReq_cueAlign))))';
    tempSupp_FrameIdx = intersect(1:response_frames, find(windowVar<(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*suppReq_cueAlign)));
    if length(tempSupp_DiffIdx)
      if strfind(tempSupp_DiffIdx, [1 1 1 1]) 
        if isempty(strfind(tempResp_DiffIdx, [1 1]))
          if eval(['frameIdx.cueAlign.' eventLabel{1,event} 'CSP_negChangeDiff{chngReq_cueAlign,c}'])
           eval(['cell.cueAlign.suppress_' eventLabel{1,event} 'CSP=[cell.cueAlign.suppress_' eventLabel{1,event} 'CSP c];']);
          end
        end
      end
    end
    end
  end
  
  %CS-
  for event=1:size(eventLabel,2)
      eval(['cell.cueAlign.response_' eventLabel{1,event} 'CSM=[];']);
      eval(['cell.cueAlign.suppress_' eventLabel{1,event} 'CSM=[];']);
    for c = 1:sum(nsize(2,:))
 eval(['windowVar = CSMdFoF_acrossTrials(' eventLabel{1,event} 'Range,c);']);
 baselineVar = CSMdFoF_acrossTrials(baselineRange,c);
 pos.Diff.cueAlign{c,1} = ([0; diff(CSMdFoF_acrossTrials(:,c))]);

pos.Limit.cueAlign{c,1}.one = (mean(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),1,'omitnan')+(std(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),[],1,'omitnan').*1)); 
pos.Change.cueAlign{c,1}.one = find(pos.Diff.cueAlign{c,1}>pos.Limit.cueAlign{c,1}.one); 
pos.Limit.cueAlign{c,1}.two = (mean(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),1,'omitnan')+(std(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),[],1,'omitnan').*2)); 
pos.Change.cueAlign{c,1}.two = find(pos.Diff.cueAlign{c,1}>pos.Limit.cueAlign{c,1}.two); 
pos.Limit.cueAlign{c,1}.three = (mean(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),1,'omitnan')+(std(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),[],1,'omitnan').*3)); 
pos.Change.cueAlign{c,1}.three = find(pos.Diff.cueAlign{c,1}>pos.Limit.cueAlign{c,1}.three);

eval(['frameIdx.cueAlign.' eventLabel{1,event} 'CSM_posChangeDiff{1,c} = intersect(' eventLabel{1,event} 'Range, pos.Change.cueAlign{c,1}.one);']); 
eval(['frameIdx.cueAlign.' eventLabel{1,event} 'CSM_posChangeDiff{2,c} = intersect(' eventLabel{1,event} 'Range, pos.Change.cueAlign{c,1}.two);']); 
eval(['frameIdx.cueAlign.' eventLabel{1,event} 'CSM_posChangeDiff{3,c} = intersect(' eventLabel{1,event} 'Range, pos.Change.cueAlign{c,1}.three);']);

neg.Diff.cueAlign{c,1} = ([0; diff(CSMdFoF_acrossTrials(:,c))]);

neg.Limit.cueAlign{c,1}.one = (mean(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),1,'omitnan')-(std(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),[],1,'omitnan').*1)); 
neg.Change.cueAlign{c,1}.one = find(neg.Diff.cueAlign{c,1}<neg.Limit.cueAlign{c,1}.one); 
neg.Limit.cueAlign{c,1}.two = (mean(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),1,'omitnan')-(std(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),[],1,'omitnan').*2)); 
neg.Change.cueAlign{c,1}.two = find(neg.Diff.cueAlign{c,1}<neg.Limit.cueAlign{c,1}.two); 
neg.Limit.cueAlign{c,1}.three = (mean(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),1,'omitnan')-(std(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),[],1,'omitnan').*3)); 
neg.Change.cueAlign{c,1}.three = find(neg.Diff.cueAlign{c,1}<neg.Limit.cueAlign{c,1}.three);

eval(['frameIdx.cueAlign.' eventLabel{1,event} 'CSM_negChangeDiff{1,c} = intersect(' eventLabel{1,event} 'Range, neg.Change.cueAlign{c,1}.one);']); 
eval(['frameIdx.cueAlign.' eventLabel{1,event} 'CSM_negChangeDiff{2,c} = intersect(' eventLabel{1,event} 'Range, neg.Change.cueAlign{c,1}.two);']); 
eval(['frameIdx.cueAlign.' eventLabel{1,event} 'CSM_negChangeDiff{3,c} = intersect(' eventLabel{1,event} 'Range, neg.Change.cueAlign{c,1}.three);']);
%    

    tempResp_DiffIdx = diff(intersect(1:response_frames, find(windowVar>(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*respReq_cueAlign))))';
    tempResp_FrameIdx = intersect(1:response_frames, find(windowVar>(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*respReq_cueAlign)));
    if length(tempResp_DiffIdx)
      if strfind(tempResp_DiffIdx, [1 1])
        if eval(['frameIdx.cueAlign.' eventLabel{1,event} 'CSM_posChangeDiff{chngReq_cueAlign,c}'])
           eval(['cell.cueAlign.response_' eventLabel{1,event} 'CSM=[cell.cueAlign.response_' eventLabel{1,event} 'CSM c];']);
        end
      end
    end
    tempSupp_DiffIdx = diff(intersect(1:response_frames, find(windowVar<(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*suppReq_cueAlign))))';
    tempSupp_FrameIdx = intersect(1:response_frames, find(windowVar<(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*suppReq_cueAlign)));
    if length(tempSupp_DiffIdx)
      if strfind(tempSupp_DiffIdx, [1 1 1 1]) 
        if isempty(strfind(tempResp_DiffIdx, [1 1]))
          if eval(['frameIdx.cueAlign.' eventLabel{1,event} 'CSM_negChangeDiff{chngReq_cueAlign,c}'])
           eval(['cell.cueAlign.suppress_' eventLabel{1,event} 'CSM=[cell.cueAlign.suppress_' eventLabel{1,event} 'CSM c];']);
          end
        end
      end
    end
    end
  end
%% data org | isolation of cells that fall into a single response group (e.g. removal of neurons which respond and suppress following cue) 
cellDiff.cueAlign.response_postCueCSP = setdiff(cell.cueAlign.response_postCueCSP,cell.cueAlign.suppress_postCueCSP);
cellDiff.cueAlign.suppress_postCueCSP = setdiff(cell.cueAlign.suppress_postCueCSP,cell.cueAlign.response_postCueCSP);
cellDiff.cueAlign.response_postRewCSP = setdiff(cell.cueAlign.response_postRewCSP,cell.cueAlign.suppress_postRewCSP);
cellDiff.cueAlign.suppress_postRewCSP = setdiff(cell.cueAlign.suppress_postRewCSP,cell.cueAlign.response_postRewCSP);
cellDiff.cueAlign.response_postCueCSM = setdiff(cell.cueAlign.response_postCueCSM,cell.cueAlign.suppress_postCueCSM);
cellDiff.cueAlign.suppress_postCueCSM = setdiff(cell.cueAlign.suppress_postCueCSM,cell.cueAlign.response_postCueCSM);
cellDiff.cueAlign.response_postRewCSM = setdiff(cell.cueAlign.response_postRewCSM,cell.cueAlign.suppress_postRewCSM);
cellDiff.cueAlign.suppress_postRewCSM = setdiff(cell.cueAlign.suppress_postRewCSM,cell.cueAlign.response_postRewCSM);

cellUnion.cueAlign.postCueCSP = union(cell.cueAlign.response_postCueCSP,cell.cueAlign.suppress_postCueCSP);
cellUnion.cueAlign.postRewCSP = union(cell.cueAlign.response_postRewCSP,cell.cueAlign.suppress_postRewCSP);
cellUnion.cueAlign.postCueCSM = union(cell.cueAlign.response_postCueCSM,cell.cueAlign.suppress_postCueCSM);
cellUnion.cueAlign.postRewCSM = union(cell.cueAlign.response_postRewCSM,cell.cueAlign.suppress_postRewCSM);

%this array is used later to contense text required for statistical tests
iterArray = {'cellDiff.cueAlign.response_postCueCSP','cellDiff.cueAlign.suppress_postCueCSP','cellDiff.cueAlign.response_postRewCSP','cellDiff.cueAlign.suppress_postRewCSP','cellDiff.cueAlign.response_postCueCSM','cellDiff.cueAlign.suppress_postCueCSM','cellDiff.cueAlign.response_postRewCSM','cellDiff.cueAlign.suppress_postRewCSM'};
%% plot     | cue-aligned representation of response groups 

respArray = {'response_postCueCSP','response_postCueCSM','suppress_postCueCSP','suppress_postCueCSM','response_postRewCSP','response_postRewCSM','suppress_postRewCSP','suppress_postRewCSM'};
titleArray = {'response postCue CS+','response postCue CS-','suppress postCue CS+','suppress postCue CS-','response postReward CS+','response postReward CS-','suppress postReward CS+','suppress postReward CS-'};
colorArray = {'k','r','b','m','k','r','b','m'};

ymax=input(['Set ymax: ']); 
ymin=input(['Set ymin: ']);
highlight = [255/255 217/255 96/255]; greylight = [.75 .75 .75]; redlight = [1 .75 .75]; windowPanes = [baselineRange(1), postCueRange(1), postRewRange(1)];
dataArray = {'CSPdFoF_acrossTrials','CSMdFoF_acrossTrials','CSPdFoF_acrossTrials','CSMdFoF_acrossTrials','CSPdFoF_acrossTrials','CSMdFoF_acrossTrials','CSPdFoF_acrossTrials','CSMdFoF_acrossTrials'};
tempFig = setFigure; hold on;
for c = 1:size(dataArray,2)  
    eval(['callTrue = ' dataArray{1,c} '(36:96,cellDiff.cueAlign.' respArray{1,c} ');']);
%     eval(['callFalse = ' dataArray{1,i} '(36:96,setdiff(sum(nsize(cue,:)),cellDiff.cueAlign.' respArray{1,i} '));']);
    eval(['dendrites = cellDiff.cueAlign.' respArray{1,c} ';']);
subplot(4,2,c); hold on;
ylim([ymin*2 ymax*2]); xline(0,'k'); xline(.767,'k'); title(titleArray{1,c});
if c<5; rectangle('Position',[tt(1,postCueRange(1)),ymin*2+.15,(1000./frameRateHz)*response_frames,ymax+5],'EdgeColor',highlight,'FaceColor',highlight); elseif c>4; rectangle('Position',[tt(1,postRewRange(1)),ymin*2+.15,(1000./frameRateHz)*response_frames,ymax+5],'EdgeColor',highlight,'FaceColor',highlight); end
shadedErrorBar_CV(tt(1,36:96),mean(callTrue,2,'omitnan').*(1000./frameRateHz),std(callTrue,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',colorArray{1,c});
%shadedErrorBar_CV(tt(1,36:96),mean(callFalse,2,'omitnan').*(1000./frameRateHz),std(callFalse,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',{'Color',[153/255 153/255 0/255]});
end
    saveas(tempFig, [output_fn 'stats\dFoF_PSTHs_' blockName '.pdf'],'pdf');

dataArray = {'CSPevents_acrossTrials','CSMevents_acrossTrials','CSPevents_acrossTrials','CSMevents_acrossTrials','CSPevents_acrossTrials','CSMevents_acrossTrials','CSPevents_acrossTrials','CSMevents_acrossTrials'};
tempFig = setFigure; hold on;
for c = 1:size(dataArray,2)  
    eval(['callTrue = ' dataArray{1,c} '(36:96,cellDiff.cueAlign.' respArray{1,c} ');']);
%     eval(['callFalse = ' dataArray{1,i} '(36:96,setdiff(sum(nsize(cue,:)),cellDiff.cueAlign.' respArray{1,i} '));']);
    eval(['dendrites = cellDiff.cueAlign.' respArray{1,c} ';']);
subplot(4,2,c); hold on;
ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title(titleArray{1,c});
if c<5; rectangle('Position',[tt(1,postCueRange(1)),0+.1,(1000./frameRateHz)*response_frames,ymax+1],'EdgeColor',highlight,'FaceColor',highlight); elseif c>4; rectangle('Position',[tt(1,postRewRange(1)),0+.1,(1000./frameRateHz)*response_frames,ymax+1],'EdgeColor',highlight,'FaceColor',highlight); end
shadedErrorBar_CV(tt(1,36:96),mean(callTrue,2,'omitnan').*(1000./frameRateHz),std(callTrue,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',colorArray{1,c});
%shadedErrorBar_CV(tt(1,36:96),mean(callFalse,2,'omitnan').*(1000./frameRateHz),std(callFalse,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',{'Color',[153/255 153/255 0/255]});
end
    saveas(tempFig, [output_fn 'stats\CspkEvents_PSTHs_' blockName '.pdf'],'pdf');

  tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
ylim([-1 ymax]);  title('Average dFoF for CS+/- trials with statistical windows');
rectangle('Position',[tt(1,baselineRange(1)),-1+.01,(1000./frameRateHz)*response_frames,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postCueRange(1)),-1+.01,(1000./frameRateHz)*response_frames,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postRewRange(1)),-1+.01,(1000./frameRateHz)*response_frames,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:96),mean(CSPdFoF_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSPdFoF_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSPdFoF_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:96),mean(CSMdFoF_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSMdFoF_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSMdFoF_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
saveas(tempFig, [output_fn 'stats\dFoF_' blockName '.pdf'],'pdf');

  tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
ylim([0 ymax]); title('Average Cspk for CS+/- trials with statistical windows');
rectangle('Position',[tt(1,baselineRange(1)),-1+.01,(1000./frameRateHz)*response_frames,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postCueRange(1)),-1+.01,(1000./frameRateHz)*response_frames,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postRewRange(1)),-1+.01,(1000./frameRateHz)*response_frames,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:96),mean(CSPevents_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSPevents_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSPevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:96),mean(CSMevents_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSMevents_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSMevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
saveas(tempFig, [output_fn 'stats\Cspk_' blockName '.pdf'],'pdf');

respArray = {'response_postCueCSP','response_postCueCSM','response_postRewCSP','response_postRewCSM';'suppress_postCueCSP','suppress_postCueCSM','suppress_postRewCSP','suppress_postRewCSM'};
dataArray = {'CSPevents_acrossTrials','CSMevents_acrossTrials','CSPevents_acrossTrials','CSMevents_acrossTrials'};

    tempFig = setFigure; hold on; c=1;
for event=1:size(eventLabel,2)
  for cue=1:size(cueLabel,2)
eval(['respN = size(' dataArray{1,c} '(36:96,cellDiff.cueAlign.' respArray{1,c} '),2);']);
eval(['suppN = size(' dataArray{1,c} '(36:96,cellDiff.cueAlign.' respArray{2,c} '),2);']);
       elseN = sum(nsize(cue,:))-(respN+suppN);

subplot(3,2,c); 
pie([respN, suppN, elseN]);
if c == sum(size(eventLabel,2)+size(cueLabel,2)); legend('resp','supp','neither'); end
c=c+1; title([eventLabel{1,event} '-' cueLabel{1,cue}]);
    
  end
end
    saveas(tempFig, [output_fn 'stats\pie_propCells_cueAlign_' blockName '.pdf'],'pdf');
    
respArray = {'response_postCueCSP','response_postRewCSP','suppress_postCueCSP','suppress_postRewCSP';'response_postCueCSM','response_postRewCSM','suppress_postCueCSM','suppress_postRewCSM'};
typeArray = {'response','suppress','resposne','suppress'};
%% ttest    | response-window avg vs. baseline avg 
  n.cueAlign.trialAvg_CSPtrials = sum(~plusData.gBlock2);
  n.cueAlign.trialAvg_CSMtrials = sum(minusData.gBlock2);
  n.cueAlign.trialAvg_CSPneurons = size(CSPevents_acrossTrials,2);
  n.cueAlign.trialAvg_CSMneurons = size(CSMevents_acrossTrials,2);
  n.cueAlign.trialAvg_mice = size(expt,2);

  [peak_cueActCSP, peakIdx_cueActCSP] = max(CSPevents_acrossTrials(postCueRange,:),[],1);
  [peak_rewActCSP, peakIdx_rewActCSP] = max(CSPevents_acrossTrials(postRewRange,:),[],1);
  [peak_baseCSP, peakIdx_baseCSP] = max(CSPevents_acrossTrials(baselineRange,:),[],1);
  [peak_cueActCSM, peakIdx_cueActCSM] = max(CSMevents_acrossTrials(postCueRange,:),[],1);
  [peak_rewActCSM, peakIdx_rewActCSM] = max(CSMevents_acrossTrials(postRewRange,:),[],1);
  [peak_baseCSM, peakIdx_baseCSM] = max(CSMevents_acrossTrials(baselineRange,:),[],1);

  [peak_cueActCSP_allN, peakIdx_cueActCSP_allN] = max(mean(CSPevents_acrossTrials(postCueRange,:),2,'omitnan'),[],1);
  [peak_rewActCSP_allN, peakIdx_rewActCSP_allN] = max(mean(CSPevents_acrossTrials(postRewRange,:),2,'omitnan'),[],1);
  [peak_baseCSP_allN, peakIdx_baseCSP_allN] = max(mean(CSPevents_acrossTrials(baselineRange,:),2,'omitnan'),[],1);
  [peak_cueActCSM_allN, peakIdx_cueActCSM_allN] = max(mean(CSMevents_acrossTrials(postCueRange,:),2,'omitnan'),[],1);
  [peak_rewActCSM_allN, peakIdx_rewActCSM_allN] = max(mean(CSMevents_acrossTrials(postRewRange,:),2,'omitnan'),[],1);
  [peak_baseCSM_allN, peakIdx_baseCSM_allN] = max(mean(CSMevents_acrossTrials(baselineRange,:),2,'omitnan'),[],1);

  m.cueAlign.CSP_baseline_trialAvg_Cspk = mean(mean(CSPevents_acrossTrials(baselineRange,:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.cueAlign.CSP_baseline_trialAvg_Cspk = std(mean(CSPevents_acrossTrials(baselineRange,:),1,'omitnan'),[],2,'omitnan')./sqrt(size(CSPevents_acrossTrials,2)).*(1000./frameRateHz);
  m.cueAlign.CSM_baseline_trialAvg_Cspk = mean(mean(CSMevents_acrossTrials(baselineRange,:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.cueAlign.CSM_baseline_trialAvg_Cspk = std(mean(CSMevents_acrossTrials(baselineRange,:),1,'omitnan'),[],2,'omitnan')./sqrt(size(CSMevents_acrossTrials,2)).*(1000./frameRateHz);
  
% difference in activity after cue (a.k.a. before reward - compared to pre-cue) - averaged trials :: SAA
for c = 1:sum(nsize(1,:))
postCue_baseline_CSP(:,c) = CSPevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSP_allN-statsWindow):baselineRange(1)-1+(peakIdx_baseCSP_allN+statsWindow),c);
postCue_activity_CSP(:,c) = CSPevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSP_allN+statsWindow),c);
end
[~, p.cueAlign.CSP_postCue_trialAvg_Cspk, ~, stats.cueAlign.CSP_postCue_trialAvg_Cspk] = ttest(mean(postCue_baseline_CSP,1,'omitnan'),mean(postCue_activity_CSP,1,'omitnan'),'tail','left');
  m.cueAlign.CSP_postCue_trialAvg_Cspk = mean(mean(CSPevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSP_allN+statsWindow),:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.cueAlign.CSP_postCue_trialAvg_Cspk = std(mean(CSPevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSP_allN+statsWindow),:),1,'omitnan'),[],2,'omitnan')./sqrt(size(CSPevents_acrossTrials,2)).*(1000./frameRateHz);
for c = 1:sum(nsize(2,:))
postCue_baseline_CSM(:,c) = CSMevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSM_allN-statsWindow):baselineRange(1)-1+(peakIdx_baseCSM_allN+statsWindow),c);
postCue_activity_CSM(:,c) = CSMevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSM_allN-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSM_allN+statsWindow),c); 
end
[~, p.cueAlign.CSM_postCue_trialAvg_Cspk, ~, stats.cueAlign.CSM_postCue_trialAvg_Cspk] = ttest(mean(postCue_baseline_CSM,1,'omitnan'),mean(postCue_activity_CSM,1,'omitnan'),'tail','left');
  m.cueAlign.CSM_postCue_trialAvg_Cspk = mean(mean(CSMevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSM_allN-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSM_allN+statsWindow),:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.cueAlign.CSM_postCue_trialAvg_Cspk = std(mean(CSMevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSM_allN-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSM_allN+statsWindow),:),1,'omitnan'),[],2,'omitnan')./sqrt(size(CSMevents_acrossTrials,2)).*(1000./frameRateHz);

% difference in activity after reward (compared to post-cue) - averaged trials :: comparing activity across neurons
for c = 1:sum(nsize(1,:))
postReward_baseline_CSP(:,c) = CSPevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSP_allN-statsWindow):baselineRange(1)-1+(peakIdx_baseCSP_allN+statsWindow),c);
postReward_activity_CSP(:,c) = CSPevents_acrossTrials(postRewRange(1)-1+(peakIdx_rewActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_rewActCSP_allN+statsWindow),c);
end
[~, p.cueAlign.CSP_postRew_trialAvg_Cspk, ~, stats.cueAlign.CSP_postRew_trialAvg_Cspk] = ttest(mean(postReward_baseline_CSP,1,'omitnan'),mean(postReward_activity_CSP,1,'omitnan'),'tail','left');
  m.cueAlign.CSP_postRew_trialAvg_Cspk = mean(mean(CSPevents_acrossTrials(postRewRange(1)-1+(peakIdx_rewActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_rewActCSP_allN+statsWindow),:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.cueAlign.CSP_postRew_trialAvg_Cspk = std(mean(CSPevents_acrossTrials(postRewRange(1)-1+(peakIdx_rewActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_rewActCSP_allN+statsWindow),:),1,'omitnan'),[],2,'omitnan')./sqrt(size(CSPevents_acrossTrials,2)).*(1000./frameRateHz);
for c = 1:sum(nsize(2,:))
postReward_baseline_CSM(:,c) = CSMevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSM_allN-statsWindow):baselineRange(1)-1+(peakIdx_baseCSM_allN+statsWindow),c);  
postReward_activity_CSM(:,c) = CSMevents_acrossTrials(postRewRange(1)-1+(peakIdx_rewActCSM_allN-statsWindow):postRewRange(1)-1+(peakIdx_rewActCSM_allN+statsWindow) ,c);   
end
[~, p.cueAlign.CSM_postRew_trialAvg_Cspk, ~, stats.cueAlign.CSM_postRew_trialAvg_Cspk] = ttest(mean(postReward_baseline_CSM,1,'omitnan'),mean(postReward_activity_CSM,1,'omitnan'),'tail','left');
  m.cueAlign.CSM_postRew_trialAvg_Cspk = mean(mean(CSMevents_acrossTrials(postRewRange(1)-1+(peakIdx_rewActCSM_allN-statsWindow):postRewRange(1)-1+(peakIdx_rewActCSM_allN+statsWindow),:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.cueAlign.CSM_postRew_trialAvg_Cspk = std(mean(CSMevents_acrossTrials(postRewRange(1)-1+(peakIdx_rewActCSM_allN-statsWindow):postRewRange(1)-1+(peakIdx_rewActCSM_allN+statsWindow),:),1,'omitnan'),[],2,'omitnan')./sqrt(size(CSPevents_acrossTrials,2)).*(1000./frameRateHz);

  %%    analysis for modulation (CS+ - CS- / CS+ + CS-)
    
  allCell_peakCSP = mean(mean(CSPevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSP-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSP+statsWindow),:,:),3,'omitnan'),1,'omitnan').*(1000./frameRateHz);
  allCell_peakCSM = mean(mean(CSMevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSM-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSM+statsWindow),:,:),3,'omitnan'),1,'omitnan').*(1000./frameRateHz);
  
  sig_allCell_peakCSP = sum(allCell_peakCSP >= 2*mean(postCue_baseline_CSP,1).*(1000./frameRateHz));
  sig_allCell_peakCSM = sum(allCell_peakCSM >= 2*mean(postCue_baseline_CSM,1).*(1000./frameRateHz));
  

%% plot     | statistical window for reward-window avg vs. baseline avg
peak_baselinePane(1,1) = round(mean(peakIdx_baseCSP_allN,2,'omitnan')); peak_baselinePane(1,2) = round(mean(peakIdx_baseCSM_allN,2,'omitnan'));
peak_postCuePane(1,1) = round(mean(peakIdx_cueActCSP_allN,2,'omitnan')); peak_postCuePane(1,2) = round(mean(peakIdx_cueActCSM_allN,2,'omitnan'));
peak_postRewPane(1,1) = round(mean(peakIdx_rewActCSP_allN,2,'omitnan')); peak_postRewPane(1,2) = round(mean(peakIdx_rewActCSM_allN,2,'omitnan'));


  tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('dFoF'); xlabel('Time from cue onset (s)')
subplot(2,1,1); hold on;
 ylim([ymin*2 ymax*2]); xline(0,'k'); xline(.767,'k'); title('Average dFoF for CS+ trials with statistical windows');
rectangle('Position',[tt(1,baselineRange(1)-1+peak_baselinePane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSPdFoF_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSPdFoF_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSPdFoF_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
subplot(2,1,2); hold on;
 ylim([ymin*2 ymax*2]); xline(0,'k'); xline(.767,'k'); title('Average dFoF for CS- trials with statistical windows');
rectangle('Position',[tt(1,baselineRange(1)-1+peak_baselinePane(1,2)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,2)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,2)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSMdFoF_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSMdFoF_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSMdFoF_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
saveas(tempFig, [output_fn 'stats\dFoFWindow_' blockName '.pdf'],'pdf');

  tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
subplot(2,1,1); hold on;
 ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('Average Cpsk for CS+ trials with statistical windows');
rectangle('Position',[tt(1,baselineRange(1)-1+peak_baselinePane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSPevents_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSPevents_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSPevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
subplot(2,1,2); hold on;
 ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('Average cSpk for CS- trials with statistical windows');
rectangle('Position',[tt(1,baselineRange(1)-1+peak_baselinePane(1,2)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,2)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,2)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSMevents_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSMevents_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSMevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
saveas(tempFig, [output_fn 'stats\CspkWindow_' blockName '.pdf'],'pdf');

  tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
 subplot(2,1,1); hold on;
ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('Average Cspk for CS+ trials with statistical windows');
rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
for i = 1:sum(nsize(1,:))
    plot(tt(1,36:96)/1000,CSPevents_acrossTrials(36:96,i).*(100./frameRateHz),'Color',greylight);
end
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSPevents_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSPevents_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSPevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');

subplot(2,1,2); hold on; 
ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('Average Cspk for CS- trials with statistical windows');
rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,2)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
for i = 1:sum(nsize(2,:))
    plot(tt(1,36:96)/1000,CSMevents_acrossTrials(36:96,i).*(100./frameRateHz),'Color',redlight);
end
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSMevents_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSMevents_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSMevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
saveas(tempFig, [output_fn 'stats\indivTraces_CspkWindow_' blockName '.pdf'],'pdf');

%% plot     | split plot cue-align 
 %plot cue-aligned
 gAnimalTrials = [(plusData.animalTrials),(minusData.animalTrials)];
 gNSize = [(plusData.nsize),(minusData.nsize)];
    

    cueAlignEvents_CSP = NaN(150,sum(plusData.nsize),sum((plusData.animalTrials))); nNP = 1; nTP = 1;
    cueAlignEvents_CSM = NaN(150,sum(minusData.nsize),sum((minusData.animalTrials))); nNM = 1; nTM = 1;



    for c = 1:2
    for i = 1:size(expt,2)
        if c==1
        cueAlignEvents_CSP(:,nNP:nNP+plusData.nsize(1,i)-1,nTP:nTP+plusData.animalTrials(1,i)-1) = plusData.groupAlign.events{1,i};
        nNP = nNP+plusData.nsize(1,i); nTP = nTP+plusData.animalTrials(1,i);
        elseif c==2
        cueAlignEvents_CSM(:,nNM:nNM+minusData.nsize(1,i)-1,nTM:nTM+minusData.animalTrials(1,i)-1) = minusData.groupAlign.events{1,i};
        nNM = nNM+minusData.nsize(1,i); nTM = nTM+minusData.animalTrials(1,i);
        end
    end
    end
    

    % split trial numbers into N phases
    splitBy = 3;
    for mus = 1:size(expt,2)
        b_temp = find(minusData.block2{1,mus}); b_split = floor(size(b_temp,2)./splitBy); b_range = 1:b_split:size(b_temp,2); 
        r_temp = find(~plusData.block2{1,mus}); r_split = floor(size(r_temp,2)./splitBy); r_range = 1:r_split:size(r_temp,2); 
    for i = 1:splitBy
        rewSplit{mus,i} = r_temp(1,r_range(1,i):(r_range(1,i)-1+r_split));
        b2Split{mus,i} = b_temp(1,b_range(1,i):(b_range(1,i)-1+b_split));
    end
    end

    
    gATP=0;gATM=0;
    for i = 1:size(plusData.animalTrials,2)
        gATP = [gATP gATP(1,i)-1+plusData.animalTrials(1,i)];
    end
    for i = 1:size(minusData.animalTrials,2)
        gATM = [gATM gATM(1,i)-1+minusData.animalTrials(1,i)];
    end
    for i = 1:splitBy
        gRewSplit{1,i}=[]; gB2Split{1,i}=[]; 
        for mus = 1:size(expt,2)
        gRewSplit{1,i} = [gRewSplit{1,i} (rewSplit{mus,i}+(gATP(1,mus)))];
        gB2Split{1,i} = [gB2Split{1,i} (b2Split{mus,i}+(gATM(1,mus)))];
        end
    end
    
    tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
    for i = 1:splitBy
        subplot(splitBy,1,i); hold on;
    shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents_CSP(36:96,:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(cueAlignEvents_CSP(36:96,:,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(sum(plusData.nsize)),'lineProps','k');
    shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents_CSM(36:96,:,gB2Split{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(cueAlignEvents_CSM(36:96,:,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(sum(minusData.nsize)),'lineProps','r');
    hold on;
    xlabel('Time from cue (s)');
    ylabel('Cspk rate (Hz)');
    ylim([0 ymax]);  vline(0,'k'); vline(.767,'k');
    end
    sgtitle([num2str(sum(gNSize)) ' neurons from ' num2str(size(expt,2)) ' animals | Cspk relative to cue']);
    saveas(tempFig, [output_fn 'stats\cueAlign_splits' num2str(splitBy) '_' blockName '.pdf'],'pdf'); 

    tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
    shadedErrorBar_CV(tt(1,36:96)/1000, mean(eventsCSP_acrossNeurons(36:96,:),2,'omitnan').*(1000./frameRateHz), (std(eventsCSP_acrossNeurons(36:96,:),[],2,'omitnan').*(1000./frameRateHz))./sqrt(sum(plusData.nsize)),'lineProps','k');
    shadedErrorBar_CV(tt(1,36:96)/1000, mean(eventsCSM_acrossNeurons(36:96,:),2,'omitnan').*(1000./frameRateHz), (std(eventsCSM_acrossNeurons(36:96,:),[],2,'omitnan').*(1000./frameRateHz))./sqrt(sum(minusData.nsize)),'lineProps','r');
    hold on;
    xlabel('Time from cue (s)');
    ylabel('Cspk rate (Hz)');
    ylim([0 ymax]);  vline(0,'k'); vline(.767,'k');
    sgtitle([num2str(sum(gNSize)) ' neurons from ' num2str(size(expt,2)) ' animals | Cspk relative to cue']);
 
 %% ANOVA for cue align comparisons across split
         for i = 1:splitBy
                   
         [peak_cueAlignSplit_postCueCSP(1,i), peakIdx_cueAlignSplit_postCueCSP(1,i)] = max(mean(mean(cueAlignEvents_CSP(postCueRange,:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan'),[],1);
         [peak_cueAlignSplit_postCueCSM(1,i), peakIdx_cueAlignSplit_postCueCSM(1,i)] = max(mean(mean(cueAlignEvents_CSM(postCueRange,:,gB2Split{1,i}),3,'omitnan'),2,'omitnan'),[],1);
         [peak_cueAlignSplit_postRewCSP(1,i), peakIdx_cueAlignSplit_postRewCSP(1,i)] = max(mean(mean(cueAlignEvents_CSP(postRewRange,:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan'),[],1);
         [peak_cueAlignSplit_postRewCSM(1,i), peakIdx_cueAlignSplit_postRewCSM(1,i)] = max(mean(mean(cueAlignEvents_CSM(postRewRange,:,gB2Split{1,i}),3,'omitnan'),2,'omitnan'),[],1);
         [peak_baselineSplitCSP(1,i), peakIdx_baselineSplitCSP(1,i)] = max(mean(mean(cueAlignEvents_CSP(baselineRange,:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan'),[],1);
         [peak_baselineSplitCSM(1,i), peakIdx_baselineSplitCSM(1,i)] = max(mean(mean(cueAlignEvents_CSM(baselineRange,:,gB2Split{1,i}),3,'omitnan'),2,'omitnan'),[],1);

         m.cueAlignSplit.CSP_postCue_trialAvg_Cspk(1,i) = mean(mean(mean(cueAlignEvents_CSP((postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz));
         s.cueAlignSplit.CSP_postCue_trialAvg_Cspk(1,i) = mean(std(mean(cueAlignEvents_CSP((postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(sum(plusData.nsize)).*(1000./frameRateHz));
         m.cueAlignSplit.CSP_postRew_trialAvg_Cspk(1,i) = mean(mean(mean(cueAlignEvents_CSP((postRewRange(peakIdx_cueAlignSplit_postRewCSP(1,i))-statsWindow):(postRewRange(peakIdx_cueAlignSplit_postRewCSP(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz));
         s.cueAlignSplit.CSP_postRew_trialAvg_Cspk(1,i) = mean(std(mean(cueAlignEvents_CSP((postRewRange(peakIdx_cueAlignSplit_postRewCSP(1,i))-statsWindow):(postRewRange(peakIdx_cueAlignSplit_postRewCSP(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(sum(plusData.nsize)).*(1000./frameRateHz));
         m.cueAlignSplit.CSP_baseline_trialAvg_Cspk(1,i) = mean(mean(mean(cueAlignEvents_CSP((baselineRange(peakIdx_baselineSplitCSP(1,i))-statsWindow):(baselineRange(peakIdx_baselineSplitCSP(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz));
         s.cueAlignSplit.CSP_baseline_trialAvg_Cspk(1,i) = mean(std(mean(cueAlignEvents_CSP((baselineRange(peakIdx_baselineSplitCSP(1,i))-statsWindow):(baselineRange(peakIdx_baselineSplitCSP(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(sum(plusData.nsize)).*(1000./frameRateHz));
         
         m.cueAlignSplit.CSM_postCue_trialAvg_Cspk(1,i) = mean(mean(mean(cueAlignEvents_CSM((postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))+statsWindow),:,gB2Split{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz));
         s.cueAlignSplit.CSM_postCue_trialAvg_Cspk(1,i) = mean(std(mean(cueAlignEvents_CSM((postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))+statsWindow),:,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(sum(minusData.nsize)).*(1000./frameRateHz));
         m.cueAlignSplit.CSM_postRew_trialAvg_Cspk(1,i) = mean(mean(mean(cueAlignEvents_CSM((postRewRange(peakIdx_cueAlignSplit_postRewCSM(1,i))-statsWindow):(postRewRange(peakIdx_cueAlignSplit_postRewCSM(1,i))+statsWindow),:,gB2Split{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz));
         s.cueAlignSplit.CSM_postRew_trialAvg_Cspk(1,i) = mean(std(mean(cueAlignEvents_CSM((postRewRange(peakIdx_cueAlignSplit_postRewCSM(1,i))-statsWindow):(postRewRange(peakIdx_cueAlignSplit_postRewCSM(1,i))+statsWindow),:,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(sum(minusData.nsize)).*(1000./frameRateHz));
         m.cueAlignSplit.CSM_baseline_trialAvg_Cspk(1,i) = mean(mean(mean(cueAlignEvents_CSM((baselineRange(peakIdx_baselineSplitCSM(1,i))-statsWindow):(baselineRange(peakIdx_baselineSplitCSM(1,i))+statsWindow),:,gB2Split{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz));
         s.cueAlignSplit.CSM_baseline_trialAvg_Cspk(1,i) = mean(std(mean(cueAlignEvents_CSM((baselineRange(peakIdx_baselineSplitCSM(1,i))-statsWindow):(baselineRange(peakIdx_baselineSplitCSM(1,i))+statsWindow),:,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(sum(minusData.nsize)).*(1000./frameRateHz));
         
             
         end

         tempFig=setFigure; hold on; 
        plot(1:splitBy, m.cueAlignSplit.CSP_postCue_trialAvg_Cspk, 'k', 'linewidth', 2);
        errorbar(1:splitBy, m.cueAlignSplit.CSP_postCue_trialAvg_Cspk, s.cueAlignSplit.CSP_postCue_trialAvg_Cspk, 'k', 'linewidth', 1.5);
        plot(1:splitBy, m.cueAlignSplit.CSM_postCue_trialAvg_Cspk, 'r', 'linewidth', 2);
        errorbar(1:splitBy, m.cueAlignSplit.CSM_postCue_trialAvg_Cspk, s.cueAlignSplit.CSM_postCue_trialAvg_Cspk, 'r', 'linewidth', 1.5);
        xlabel('Phase');
        ylabel('Cspk rate (Hz)');
        sgtitle([num2str(sum(gNSize)) ' neurons from ' num2str(size(expt,2)) ' animals | cue aligned peak Cspk rate after cue by phase']);
        ylim([0 ymax]); xticks([1 2 3]);
        saveas(tempFig, [output_fn 'stats\postCue_cueAlignPeak_splits' num2str(splitBy)  '_' blockName '_lineGraph.pdf'],'pdf'); 
 
       
            for i = 1:splitBy
            for c = 1:size(cueAlignEvents_CSP,2)
                % for f = 1:size(cueAlignEvents,1)
                cSpkRate_CSP_split{1,i}(:,c) = mean(cueAlignEvents_CSP(:,c,gRewSplit{1,i}),3,'omitnan');%==1)./sum(~isnan(cueAlignEvents(:,c,gRewSplit{1,i})));
                    cSpkRate_CSPtrials_split{1,i}(:,c,:) = cueAlignEvents_CSP(:,c,gRewSplit{1,i});
            end
            for c = 1:size(cueAlignEvents_CSM,2)
                cSpkRate_CSM_split{1,i}(:,c) = mean(cueAlignEvents_CSM(:,c,gB2Split{1,i}),3,'omitnan');%==1)./sum(~isnan(cueAlignEvents(:,c,gB2Split{1,i})));
                    cSpkRate_CSMtrials_split{1,i}(:,c,:) = cueAlignEvents_CSM(:,c,gB2Split{1,i});
                % end
            end
            end
        

        for i = 1:splitBy
        cueAlign_cSpk_CSP_split{1,i} = mean(cSpkRate_CSP_split{1,i}(postCueRange(1)-1+(peakIdx_cueAlignSplit_postCueCSP-statsWindow):postCueRange(1)-1+(peakIdx_cueAlignSplit_postCueCSP+statsWindow),:), 1,'omitnan');
        cueAlign_cSpk_CSP_splitData{1,i} = cSpkRate_CSP_split{1,i}(postCueRange(1)-1+(peakIdx_cueAlignSplit_postCueCSP-statsWindow):postCueRange(1)-1+(peakIdx_cueAlignSplit_postCueCSP+statsWindow),:);% mean(, 1,'omitnan');
            cueAlign_cSpk_CSP_splitAllData{1,i} = squeeze(cSpkRate_CSPtrials_split{1,i}(postCueRange(1)-1+(peakIdx_cueAlignSplit_postCueCSP-statsWindow):postCueRange(1)-1+(peakIdx_cueAlignSplit_postCueCSP+statsWindow),:,:));% mean(, 1,'omitnan');
        cueAlign_cSpk_CSM_split{1,i} = mean(cSpkRate_CSM_split{1,i}(postCueRange(1)-1+(peakIdx_cueAlignSplit_postCueCSM-statsWindow):postCueRange(1)-1+(peakIdx_cueAlignSplit_postCueCSM+statsWindow),:), 1,'omitnan');
        cueAlign_cSpk_CSM_splitData{1,i} = cSpkRate_CSM_split{1,i}(postCueRange(1)-1+(peakIdx_cueAlignSplit_postCueCSM-statsWindow):postCueRange(1)-1+(peakIdx_cueAlignSplit_postCueCSM+statsWindow),:);% mean(, 1,'omitnan');
            cueAlign_cSpk_CSM_splitAllData{1,i} = squeeze(cSpkRate_CSM_split{1,i}(postCueRange(1)-1+(peakIdx_cueAlignSplit_postCueCSM-statsWindow):postCueRange(1)-1+(peakIdx_cueAlignSplit_postCueCSM+statsWindow),:,:));% mean(, 1,'omitnan');
        end

        phaseID={'first','second','third'}; nPhase=size(phaseID,2); cueID={'CS+','CS-'}; nCue=size(cueID,2); if strcmp(blockName,'allPlus'); cells = gNSize; else; cells = sum(plusData.nsize+minusData.nsize); end
        cueAlign_cSpk_firstArray = [mean(cueAlign_cSpk_CSM_splitData{1,1}*1000/frameRateHz,1,'omitnan'),mean(cueAlign_cSpk_CSP_splitData{1,1}*1000/frameRateHz,1,'omitnan')];
        cueAlign_cSpk_secondArray = [mean(cueAlign_cSpk_CSM_splitData{1,2}*1000/frameRateHz,1,'omitnan'),mean(cueAlign_cSpk_CSP_splitData{1,2}*1000/frameRateHz,1,'omitnan')];
        cueAlign_cSpk_thirdArray = [mean(cueAlign_cSpk_CSM_splitData{1,3}*1000/frameRateHz,1,'omitnan'),mean(cueAlign_cSpk_CSP_splitData{1,3}*1000/frameRateHz,1,'omitnan')];
        cueAlign_cSpk_array{1,1} = cueAlign_cSpk_firstArray; cueAlign_cSpk_array{1,2} = cueAlign_cSpk_secondArray; cueAlign_cSpk_array{1,3} = cueAlign_cSpk_thirdArray;
        
        
        dataArray = NaN(sum(plusData.nsize+minusData.nsize),nPhase);
        
%        dataArray(2, :, :, :) = NaN(nCue, nPhase, nIC);% Initialize the 3D array
        for i = 1:length(phaseID)
          dataArray(:,i) = cueAlign_cSpk_array{1,i};
        end


%dataForTable=[cueAlign_wholeResponseData(:,1),cueAlign_wholeResponseData(:,2),cueAlign_wholeResponseData(:,3)];
% Create a table reflecting the within subject factors 
cue = categorical([repmat({'CS-'}, sum(minusData.nsize),1); repmat({'CS+'},sum(plusData.nsize),1)]);
cue = reordercats(cue, {'CS-', 'CS+'});  % Set CS- as the baseline
cueRep = repmat(cue,3,1);
phase = categorical([1 2 3]);
phases = categorical([ones(length(cue), 1); ones(length(cue), 1)*2; ones(length(cue), 1)*3]);
cellID = [repmat(1:length(cue),3,1)]'; cellID=categorical(cellID(:));
dataArray=dataArray(:);
t = table(cellID,cueRep,phases,dataArray,VariableNames=["cell","cue","phase","response"]);
% 3. Call fitrm with the modified within design.
lme_cueAlign_split = fitlme(t,'response ~ cue*phase + (1|cell)');
%resLME = compare(lme_cueAlign_split);

disp('postCueAlign LME')
lmeTable=dataset2table(lme_cueAlign_split.Coefficients);
disp(lmeTable)
writetable(lmeTable,[output_fn 'stats\' blockName '_lme_cueAlign_split.xlsx']);


%% data org | lick-aligned representation of response groups 
eventLabel={'first','postCue','postRew'}; cueLabel={'CSP','CSM'}; preLick_frames=15; preLickRange=(1:bx_response_frames)+(preLick_frames-bx_response_frames);
  %CS+
  for event=1:size(eventLabel,2)
      eval(['cell.lickAlign.response_' eventLabel{1,event} 'CSP=[];']);
      eval(['cell.lickAlign.suppress_' eventLabel{1,event} 'CSP=[];']);
    for c = 1:sum(nsize(1,:))
%
 eval(['preLickVar = ' eventLabel{1,event} '_lickAlignDFoFCSP_acrossTrials(1:preLick_frames,c);']);
 baselineVar = CSPdFoF_acrossTrials(baselineRange,c);
 eval(['pos.Diff.lickAlign{c,1} = ([0; diff(' eventLabel{1,event} '_lickAlignDFoFCSP_acrossTrials(:,c))]);']);
 
pos.Limit.lickAlign{c,1}.one = (mean(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),1,'omitnan')+(std(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),[],1,'omitnan').*1)); 
pos.Change.lickAlign{c,1}.one = find(pos.Diff.lickAlign{c,1}>pos.Limit.lickAlign{c,1}.one); 
pos.Limit.lickAlign{c,1}.two = (mean(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),1,'omitnan')+(std(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),[],1,'omitnan').*1.5)); 
pos.Change.lickAlign{c,1}.two = find(pos.Diff.lickAlign{c,1}>pos.Limit.lickAlign{c,1}.two); 
pos.Limit.lickAlign{c,1}.three = (mean(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),1,'omitnan')+(std(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),[],1,'omitnan').*2)); 
pos.Change.lickAlign{c,1}.three = find(pos.Diff.lickAlign{c,1}>pos.Limit.lickAlign{c,1}.three);

eval(['frameIdx.lickAlign.' eventLabel{1,event} 'CSP_posChangeDiff{1,c} = intersect(preLickRange, pos.Change.lickAlign{c,1}.one);']); 
eval(['frameIdx.lickAlign.' eventLabel{1,event} 'CSP_posChangeDiff{2,c} = intersect(preLickRange, pos.Change.lickAlign{c,1}.two);']); 
eval(['frameIdx.lickAlign.' eventLabel{1,event} 'CSP_posChangeDiff{3,c} = intersect(preLickRange, pos.Change.lickAlign{c,1}.three);']);

 eval(['neg.Diff.lickAlign{c,1} = ([0; diff(' eventLabel{1,event} '_lickAlignDFoFCSP_acrossTrials(:,c))]);']);
neg.Limit.lickAlign{c,1}.one = (mean(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),1,'omitnan')-(std(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),[],1,'omitnan').*1)); 
neg.Change.lickAlign{c,1}.one = find(neg.Diff.lickAlign{c,1}<neg.Limit.lickAlign{c,1}.one); 
neg.Limit.lickAlign{c,1}.two = (mean(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),1,'omitnan')-(std(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),[],1,'omitnan').*1.5)); 
neg.Change.lickAlign{c,1}.two = find(neg.Diff.lickAlign{c,1}<neg.Limit.lickAlign{c,1}.two); 
neg.Limit.lickAlign{c,1}.three = (mean(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),1,'omitnan')-(std(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),[],1,'omitnan').*2)); 
neg.Change.lickAlign{c,1}.three = find(neg.Diff.lickAlign{c,1}<neg.Limit.lickAlign{c,1}.three);

eval(['frameIdx.lickAlign.' eventLabel{1,event} 'CSP_negChangeDiff{1,c} = intersect(preLickRange, neg.Change.lickAlign{c,1}.one);']); 
eval(['frameIdx.lickAlign.' eventLabel{1,event} 'CSP_negChangeDiff{2,c} = intersect(preLickRange, neg.Change.lickAlign{c,1}.two);']); 
eval(['frameIdx.lickAlign.' eventLabel{1,event} 'CSP_negChangeDiff{3,c} = intersect(preLickRange, neg.Change.lickAlign{c,1}.three);']);

%    
    tempResp_DiffIdx = diff(intersect(preLickRange, find(preLickVar>(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*respReq_lickAlign))))';
    tempResp_FrameIdx = intersect(preLickRange, find(preLickVar>(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*respReq_lickAlign)));
    if length(tempResp_DiffIdx)
      if strfind(tempResp_DiffIdx, [1 1])
        if eval(['frameIdx.lickAlign.' eventLabel{1,event} 'CSP_posChangeDiff{chngReq_lickAlign,c}'])
           eval(['cell.lickAlign.response_' eventLabel{1,event} 'CSP=[cell.lickAlign.response_' eventLabel{1,event} 'CSP c];']);
        end
      end
    end
    tempSupp_DiffIdx = diff(intersect(preLickRange, find(preLickVar<(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*suppReq_lickAlign))))';
    tempSupp_FrameIdx = intersect(preLickRange, find(preLickVar<(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*suppReq_lickAlign)));
    if length(tempSupp_DiffIdx)
      if strfind(tempSupp_DiffIdx, [1 1 1 1]) 
        if isempty(strfind(tempResp_DiffIdx, [1 1]))
          if eval(['frameIdx.lickAlign.' eventLabel{1,event} 'CSP_negChangeDiff{chngReq_lickAlign,c}'])
           eval(['cell.lickAlign.suppress_' eventLabel{1,event} 'CSP=[cell.lickAlign.suppress_' eventLabel{1,event} 'CSP c];']);
          end
        end
      end
    end
    end
  end
  
  %CS-
  for event=1:size(eventLabel,2)
      eval(['cell.lickAlign.response_' eventLabel{1,event} 'CSM=[];']);
      eval(['cell.lickAlign.suppress_' eventLabel{1,event} 'CSM=[];']);
    for c = 1:sum(nsize(2,:))
%
 eval(['preLickVar = ' eventLabel{1,event} '_lickAlignDFoFCSM_acrossTrials(1:preLick_frames,c);']);
 baselineVar = CSMdFoF_acrossTrials(baselineRange,c);
 eval(['pos.Diff.lickAlign{c,1} = ([0; diff(' eventLabel{1,event} '_lickAlignDFoFCSM_acrossTrials(:,c))]);']);
 
pos.Limit.lickAlign{c,1}.one = (mean(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),1,'omitnan')+(std(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),[],1,'omitnan').*1)); 
pos.Change.lickAlign{c,1}.one = find(pos.Diff.lickAlign{c,1}>pos.Limit.lickAlign{c,1}.one); 
pos.Limit.lickAlign{c,1}.two = (mean(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),1,'omitnan')+(std(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),[],1,'omitnan').*1.5)); 
pos.Change.lickAlign{c,1}.two = find(pos.Diff.lickAlign{c,1}>pos.Limit.lickAlign{c,1}.two); 
pos.Limit.lickAlign{c,1}.three = (mean(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),1,'omitnan')+(std(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),[],1,'omitnan').*2)); 
pos.Change.lickAlign{c,1}.three = find(pos.Diff.lickAlign{c,1}>pos.Limit.lickAlign{c,1}.three);

eval(['frameIdx.lickAlign.' eventLabel{1,event} 'CSM_posChangeDiff{1,c} = intersect(preLickRange, pos.Change.lickAlign{c,1}.one);']); 
eval(['frameIdx.lickAlign.' eventLabel{1,event} 'CSM_posChangeDiff{2,c} = intersect(preLickRange, pos.Change.lickAlign{c,1}.two);']); 
eval(['frameIdx.lickAlign.' eventLabel{1,event} 'CSM_posChangeDiff{3,c} = intersect(preLickRange, pos.Change.lickAlign{c,1}.three);']);

 eval(['neg.Diff.lickAlign{c,1} = ([0; diff(' eventLabel{1,event} '_lickAlignDFoFCSM_acrossTrials(:,c))]);']);
neg.Limit.lickAlign{c,1}.one = (mean(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),1,'omitnan')-(std(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),[],1,'omitnan').*1)); 
neg.Change.lickAlign{c,1}.one = find(neg.Diff.lickAlign{c,1}<neg.Limit.lickAlign{c,1}.one); 
neg.Limit.lickAlign{c,1}.two = (mean(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),1,'omitnan')-(std(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),[],1,'omitnan').*1.5)); 
neg.Change.lickAlign{c,1}.two = find(neg.Diff.lickAlign{c,1}<neg.Limit.lickAlign{c,1}.two); 
neg.Limit.lickAlign{c,1}.three = (mean(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),1,'omitnan')-(std(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),[],1,'omitnan').*2)); 
neg.Change.lickAlign{c,1}.three = find(neg.Diff.lickAlign{c,1}<neg.Limit.lickAlign{c,1}.three);

eval(['frameIdx.lickAlign.' eventLabel{1,event} 'CSM_negChangeDiff{1,c} = intersect(preLickRange, neg.Change.lickAlign{c,1}.one);']); 
eval(['frameIdx.lickAlign.' eventLabel{1,event} 'CSM_negChangeDiff{2,c} = intersect(preLickRange, neg.Change.lickAlign{c,1}.two);']); 
eval(['frameIdx.lickAlign.' eventLabel{1,event} 'CSM_negChangeDiff{3,c} = intersect(preLickRange, neg.Change.lickAlign{c,1}.three);']);

%    
    tempResp_DiffIdx = diff(intersect(preLickRange, find(preLickVar>(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*respReq_lickAlign))))';
    tempResp_FrameIdx = intersect(preLickRange, find(preLickVar>(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*respReq_lickAlign)));
    if length(tempResp_DiffIdx)
      if strfind(tempResp_DiffIdx, [1 1])
        if eval(['frameIdx.lickAlign.' eventLabel{1,event} 'CSM_posChangeDiff{chngReq_lickAlign,c}'])
           eval(['cell.lickAlign.response_' eventLabel{1,event} 'CSM=[cell.lickAlign.response_' eventLabel{1,event} 'CSM c];']);
        end
      end
    end
    tempSupp_DiffIdx = diff(intersect(preLickRange, find(preLickVar<(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*suppReq_lickAlign))))';
    tempSupp_FrameIdx = intersect(preLickRange, find(preLickVar<(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*suppReq_lickAlign)));
    if length(tempSupp_DiffIdx)
      if strfind(tempSupp_DiffIdx, [1 1 1 1]) 
        if isempty(strfind(tempResp_DiffIdx, [1 1]))
          if eval(['frameIdx.lickAlign.' eventLabel{1,event} 'CSM_negChangeDiff{chngReq_lickAlign,c}'])
           eval(['cell.lickAlign.suppress_' eventLabel{1,event} 'CSM=[cell.lickAlign.suppress_' eventLabel{1,event} 'CSM c];']);
          end
        end
      end
    end
    end
  end

it=1; sideLabel={'response','suppress';'suppress','response'};
for cue=1:size(cueLabel,2)
  for event=1:size(eventLabel,2)
    for side=1:size(sideLabel,2)
eval(['cellDiff.lickAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} '= setdiff(cell.lickAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} ',cell.lickAlign.' sideLabel{2,side} '_'  eventLabel{1,event} cueLabel{1,cue} ');']);
it=it+1;
    end
  end
end
%
titleArray = {'CS+ reponsive firstLick','CS+ suppressive firstLick','CS+ reponsive postCueLick','CS+ suppressive postCueLick','CS+ reponsive postRewLick','CS+ suppressive postRewLick','CS- reponsive firstLick','CS- suppressive firstLick','CS- reponsive postCueLick','CS- suppressive postCueLick','CS- reponsive postRewLick','CS- suppressive postRewLick'};
colorArray = {'k','b','k','b','k','b','r','m','r','m','r','m'};
%% plot     | first-lick (response/suppression) 
tempFig = setFigure; hold on; c=1;
for cue=1:size(cueLabel,2)
  for event=1:size(eventLabel,2)
    for side=1:size(sideLabel,2)
    eval(['callTrue = ' eventLabel{1,event} '_lickAlignDFoF' cueLabel{1,cue} '_acrossTrials(:,cellDiff.lickAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} ');']);
    %eval(['callFalse = ' eventLabel{1,event} '_lickAlignDFoF' cueLabel{1,cue} '_acrossTrials(:,setdiff(sum(nsize(cue,:)),cellDiff.lickAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} '));']);
    eval(['dendrites = cellDiff.lickAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} ';']);
subplot(6,2,c); hold on;
xlim([-500 1000]); ylim([ymin*5 ymax*2]); xline(0,'k');  title(titleArray{1,c});
rectangle('Position',[tt(1,baselineRange(1))+10,ymin*2+.25,(1000./frameRateHz)*response_frames-10,ymax+5],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,37:81)/1000,mean(callTrue,2,'omitnan').*(1000./frameRateHz),std(callTrue,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',colorArray{1,c});
% shadedErrorBar_CV(tt(1,36:96)/1000,mean(callFalse,2,'omitnan'),std(callFalse,[],2,'omitnan')./sqrt(size(dendrites,2)),'lineProps',{'Color',[0/255 0/255 201/255]});
    c=c+1;
    end
  end
end
    saveas(tempFig, [output_fn 'stats\dFoF_lickAlign_PSTHs_' blockName '.pdf'],'pdf');

tempFig = setFigure; hold on; c=1;
for cue=1:size(cueLabel,2)
  for event=1:size(eventLabel,2)
    for side=1:size(sideLabel,2)
    eval(['callTrue = ' eventLabel{1,event} '_lickAlign' cueLabel{1,cue} 'events_acrossTrials(:,cellDiff.lickAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} ');']);
    %eval(['callFalse = ' eventLabel{1,event} '_lickAlign' cueLabel{1,cue} 'events_acrossTrials(:,setdiff(sum(nsize(cue,:)),cellDiff.lickAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} '));']);
    eval(['dendrites = cellDiff.lickAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} ';']);
subplot(6,2,c); hold on;
xlim([-500 1000]); ylim([ymin ymax]); xline(0,'k');  title(titleArray{1,c});
if c<5; rectangle('Position',[tt(1,baselineRange(1)),0+.15,(1000./frameRateHz)*response_frames,ymax+1],'EdgeColor',highlight,'FaceColor',highlight); elseif c>4; rectangle('Position',[tt(1,baselineRange(1)),0+.15,(1000./frameRateHz)*response_frames,ymax+1],'EdgeColor',highlight,'FaceColor',highlight); end
shadedErrorBar_CV(tt(1,37:81)/1000,mean(callTrue,2,'omitnan').*(1000./frameRateHz),std(callTrue,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',colorArray{1,c});
% shadedErrorBar_CV(tt(1,36:96)/1000,mean(callFalse,2,'omitnan').*(1000./frameRateHz),std(callFalse,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',{'Color',[0/255 0/255 201/255]});
    c=c+1;
    end
  end
end
    saveas(tempFig, [output_fn 'stats\CspkEvents_lickAlign_PSTHs_' blockName '.pdf'],'pdf');

    tempFig = setFigure; hold on; c=1;
for event=1:size(eventLabel,2)
  for cue=1:size(cueLabel,2)
eval(['respN = size(' eventLabel{1,event} '_lickAlign' cueLabel{1,cue} 'events_acrossTrials(:,cellDiff.lickAlign.response_' eventLabel{1,event} cueLabel{1,cue} '),2);']);
eval(['suppN = size(' eventLabel{1,event} '_lickAlign' cueLabel{1,cue} 'events_acrossTrials(:,cellDiff.lickAlign.suppress_' eventLabel{1,event} cueLabel{1,cue} '),2);']);
       elseN = sum(nsize(cue,:))-(respN+suppN);

subplot(3,2,c); 
pie([respN, suppN, elseN]);
c=c+1; title([eventLabel{1,event} '-' cueLabel{1,cue}]);
  end
end
legend('responsive','unresponsive','other','Position',[0.475,0.1,0.1,0.15])    

    saveas(tempFig, [output_fn 'stats\pie_propCells_lickAlign_' blockName '.pdf'],'pdf');

    
tempFig=setFigure; hold on;
subplot(3,1,1); hold on;
shadedErrorBar_CV(tt(1,37:81)/1000,mean(first_lickAlignCSPevents_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(first_lickAlignCSPevents_acrossTrials,[],2,'omitnan')./sqrt(size(first_lickAlignCSPevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,37:81)/1000,mean(first_lickAlignCSMevents_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(first_lickAlignCSMevents_acrossTrials,[],2,'omitnan')./sqrt(size(first_lickAlignCSMevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
xlim([-500 1000]); ylim([0 ymax]); xline(0,'k');  title('first lick');
subplot(3,1,2); hold on;
shadedErrorBar_CV(tt(1,37:81)/1000,mean(postCue_lickAlignCSPevents_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(postCue_lickAlignCSPevents_acrossTrials,[],2,'omitnan')./sqrt(size(postCue_lickAlignCSPevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,37:81)/1000,mean(postCue_lickAlignCSMevents_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(postCue_lickAlignCSMevents_acrossTrials,[],2,'omitnan')./sqrt(size(postCue_lickAlignCSMevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
xlim([-500 1000]); ylim([0 ymax]); xline(0,'k');  title('first postCue lick');
subplot(3,1,3); hold on;
shadedErrorBar_CV(tt(1,37:81)/1000,mean(postRew_lickAlignCSPevents_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(postRew_lickAlignCSPevents_acrossTrials,[],2,'omitnan')./sqrt(size(postRew_lickAlignCSPevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,37:81)/1000,mean(postRew_lickAlignCSMevents_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(postRew_lickAlignCSMevents_acrossTrials,[],2,'omitnan')./sqrt(size(postRew_lickAlignCSMevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
xlim([-500 1000]); ylim([0 ymax]); xline(0,'k');  title('first postRew lick');
%% ttest    | first-lick avg vs. baseline avg 
% offset start frame by 1, because we need to take a +/-1 frame window
% around the target; script breaks if the peak is on frame 1
  [peak_firstLickActCSP, peakIdx_firstLickActCSP] = max(first_lickAlignCSPevents_acrossTrials(preLickRange,:),[],1);
  [peak_postCueLickActCSP, peakIdx_postCueLickActCSP] = max(postCue_lickAlignCSPevents_acrossTrials(preLickRange,:),[],1);
  [peak_postRewLickActCSP, peakIdx_postRewLickActCSP] = max(postRew_lickAlignCSPevents_acrossTrials(preLickRange,:),[],1);
  [peak_baseCSP, peakIdx_baseCSP] = max(CSPevents_acrossTrials(baselineRange,:),[],1);
  [peak_firstLickActCSM, peakIdx_firstLickActCSM] = max(first_lickAlignCSMevents_acrossTrials(preLickRange,:),[],1);
  [peak_postCueLickActCSM, peakIdx_postCueLickActCSM] = max(postCue_lickAlignCSMevents_acrossTrials(preLickRange,:),[],1);
  [peak_postRewLickActCSM, peakIdx_postRewLickActCSM] = max(postRew_lickAlignCSMevents_acrossTrials(preLickRange,:),[],1);
  [peak_baseCSM, peakIdx_baseCSM] = max(CSMevents_acrossTrials(baselineRange,:),[],1);

  [~, peakIdx_firstLickActCSP_allN] = max(mean(first_lickAlignCSPevents_acrossTrials(preLickRange,:),2,'omitnan'),[],1);
  [~, peakIdx_postCueLickActCSP_allN] = max(mean(postCue_lickAlignCSPevents_acrossTrials(preLickRange,:),2,'omitnan'),[],1);
  [~, peakIdx_postRewLickActCSP_allN] = max(mean(postRew_lickAlignCSPevents_acrossTrials(preLickRange,:),2,'omitnan'),[],1);
  [~, peakIdx_baseCSP_allN] = max(mean(CSPevents_acrossTrials(baselineRange,:),2,'omitnan'),[],1);
  [~, peakIdx_firstLickActCSM_allN] = max(mean(first_lickAlignCSMevents_acrossTrials(preLickRange,:),2,'omitnan'),[],1);
  [~, peakIdx_postCueLickActCSM_allN] = max(mean(postCue_lickAlignCSMevents_acrossTrials(preLickRange,:),2,'omitnan'),[],1);
  [~, peakIdx_postRewLickActCSM_allN] = max(mean(postRew_lickAlignCSMevents_acrossTrials(preLickRange,:),2,'omitnan'),[],1);
  [~, peakIdx_baseCSM_allN] = max(mean(CSMevents_acrossTrials(baselineRange,:),2,'omitnan'),[],1);
% counts for lick align
n.lickAlign.CSP_first_trialAvg_count = sum(first_lickAlignCSP_acrossTrials); 
n.lickAlign.CSM_first_trialAvg_count = sum(first_lickAlignCSM_acrossTrials);
n.lickAlign.CSP_firstPostCue_trialAvg_count = sum(firstPostCue_lickAlignCSP_acrossTrials);
n.lickAlign.CSM_firstPostCue_trialAvg_count = sum(firstPostCue_lickAlignCSM_acrossTrials);
n.lickAlign.CSP_firstPostRew_trialAvg_count = sum(firstPostRew_lickAlignCSP_acrossTrials);
n.lickAlign.CSM_firstPostRew_trialAvg_count = sum(firstPostRew_lickAlignCSM_acrossTrials);
% difference in [] between first lick and baseline
for c = 1:sum(nsize(1,:))
firstLick_baseline_CSP(:,c) = CSPevents_acrossTrials(baselineRange(1)-1+peakIdx_baseCSP_allN-statsWindow:baselineRange(1)-1+peakIdx_baseCSP_allN+statsWindow,c);
firstLick_activity_CSP(:,c) = first_lickAlignCSPevents_acrossTrials(preLickRange(1)-1+peakIdx_firstLickActCSP_allN-statsWindow:preLickRange(1)-1+peakIdx_firstLickActCSP_allN+statsWindow,c);
end
[~, p.lickAlign.CSP_first_trialAvg_Cspk, ~, stats.lickAlign.CSP_first_trialAvg_Cspk] = ttest(mean(firstLick_baseline_CSP,1,'omitnan'),mean(firstLick_activity_CSP,1,'omitnan'),'tail','left');
  m.lickAlign.CSP_first_trialAvg_Cspk = mean(mean(first_lickAlignCSPevents_acrossTrials(preLickRange(1)-1+peakIdx_firstLickActCSP_allN-statsWindow:preLickRange(1)-1+peakIdx_firstLickActCSP_allN+statsWindow,:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickAlign.CSP_first_trialAvg_Cspk = std(mean(first_lickAlignCSPevents_acrossTrials(preLickRange(1)-1+peakIdx_firstLickActCSP_allN-statsWindow:preLickRange(1)-1+peakIdx_firstLickActCSP_allN+statsWindow,:),1,'omitnan'),[],2,'omitnan')./sqrt(size(first_lickAlignCSPevents_acrossTrials,2)).*(1000./frameRateHz);
for c = 1:sum(nsize(2,:))
firstLick_baseline_CSM(:,c) = CSMevents_acrossTrials(baselineRange(1)-1+peakIdx_baseCSM_allN-statsWindow:baselineRange(1)-1+peakIdx_baseCSM_allN+statsWindow,c);
firstLick_activity_CSM(:,c) = first_lickAlignCSMevents_acrossTrials(preLickRange(1)-1+peakIdx_firstLickActCSM_allN-statsWindow:preLickRange(1)-1+peakIdx_firstLickActCSM_allN+statsWindow,c);
end
[~, p.lickAlign.CSM_first_trialAvg_Cspk, ~, stats.lickAlign.CSM_first_trialAvg_Cspk] = ttest(mean(firstLick_baseline_CSM,1,'omitnan'),mean(firstLick_activity_CSM,1,'omitnan'),'tail','left');
  m.lickAlign.CSM_first_trialAvg_Cspk = mean(mean(first_lickAlignCSMevents_acrossTrials(preLickRange(1)-1+peakIdx_firstLickActCSM_allN-statsWindow:preLickRange(1)-1+peakIdx_firstLickActCSM_allN+statsWindow,:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickAlign.CSM_first_trialAvg_Cspk = std(mean(first_lickAlignCSMevents_acrossTrials(preLickRange(1)-1+peakIdx_firstLickActCSM_allN-statsWindow:preLickRange(1)-1+peakIdx_firstLickActCSM_allN+statsWindow,:),1,'omitnan'),[],2,'omitnan')./sqrt(size(first_lickAlignCSMevents_acrossTrials,2)).*(1000./frameRateHz);

% difference in activity between first lick [POST REWARD] and baseline
for c = 1:sum(nsize(1,:))
postRewLick_baseline_CSP(:,c) = CSPevents_acrossTrials(baselineRange(1)-1+peakIdx_baseCSP_allN-statsWindow:(baselineRange(1))+peakIdx_baseCSP_allN+statsWindow,c);
postRewLick_activity_CSP(:,c) = postRew_lickAlignCSPevents_acrossTrials(preLickRange(1)-1+peakIdx_postRewLickActCSP_allN-statsWindow:preLickRange(1)-1+peakIdx_postRewLickActCSP_allN+statsWindow,c);
end
[~, p.lickAlign.CSP_postRew_trialAvg_Cspk, ~, stats.lickAlign.CSP_postRew_trialAvg_Cspk] = ttest(mean(postRewLick_baseline_CSP,1,'omitnan'),mean(postRewLick_activity_CSP,1,'omitnan'),'tail','left');
  m.lickAlign.CSP_postRew_trialAvg_Cspk = mean(mean(postRew_lickAlignCSPevents_acrossTrials(preLickRange(1)-1+peakIdx_postRewLickActCSP_allN-statsWindow:preLickRange(1)-1+peakIdx_postRewLickActCSP_allN+statsWindow,:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickAlign.CSP_postRew_trialAvg_Cspk = std(mean(postRew_lickAlignCSPevents_acrossTrials(preLickRange(1)-1+peakIdx_postRewLickActCSP_allN-statsWindow:preLickRange(1)-1+peakIdx_postRewLickActCSP_allN+statsWindow,:),1,'omitnan'),[],2,'omitnan')./sqrt(size(postRew_lickAlignCSPevents_acrossTrials,2)).*(1000./frameRateHz);
for c = 1:sum(nsize(2,:))
postRewLick_baseline_CSM(:,c) = CSMevents_acrossTrials(baselineRange(1)-1+peakIdx_baseCSM_allN-statsWindow:baselineRange(1)-1+peakIdx_baseCSM_allN+statsWindow,c);
postRewLick_activity_CSM(:,c) = postRew_lickAlignCSMevents_acrossTrials(preLickRange(1)-1+peakIdx_postRewLickActCSM_allN-statsWindow:preLickRange(1)-1+peakIdx_postRewLickActCSM_allN+statsWindow,c); %preLickRange
end
[~, p.lickAlign.CSM_postRew_trialAvg_Cspk, ~, stats.lickAlign.CSM_postRew_trialAvg_Cspk] = ttest(mean(postRewLick_baseline_CSM,1,'omitnan'),mean(postRewLick_activity_CSM,1,'omitnan'),'tail','left');
  m.lickAlign.CSM_postRew_trialAvg_Cspk = mean(mean(postRew_lickAlignCSMevents_acrossTrials(preLickRange(1)-1+peakIdx_postRewLickActCSM_allN-statsWindow:preLickRange(1)-1+peakIdx_postRewLickActCSM_allN+statsWindow,:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickAlign.CSM_postRew_trialAvg_Cspk = std(mean(postRew_lickAlignCSMevents_acrossTrials(preLickRange(1)-1+peakIdx_postRewLickActCSM_allN-statsWindow:preLickRange(1)-1+peakIdx_postRewLickActCSM_allN+statsWindow,:),1,'omitnan'),[],2,'omitnan')./sqrt(size(postRew_lickAlignCSMevents_acrossTrials,2)).*(1000./frameRateHz);
  
% difference in activity between first lick [POST CUE] and baseline
for c = 1:sum(nsize(1,:))
postCueLick_baseline_CSP(:,c) = CSPevents_acrossTrials(baselineRange(1)-1+peakIdx_baseCSP_allN-statsWindow:baselineRange(1)-1+peakIdx_baseCSP_allN+statsWindow,c);
postCueLick_activity_CSP(:,c) = postCue_lickAlignCSPevents_acrossTrials(preLickRange(1)-1+peakIdx_postCueLickActCSP_allN-statsWindow:preLickRange(1)-1+peakIdx_postCueLickActCSP_allN+statsWindow,c);
end
[~, p.lickAlign.CSP_postCue_trialAvg_Cspk, ~, stats.lickAlign.CSP_postCue_trialAvg_Cspk] = ttest(mean(postCueLick_baseline_CSP,1,'omitnan'),mean(postCueLick_activity_CSP,1,'omitnan'),'tail','left');
  m.lickAlign.CSP_postCue_trialAvg_Cspk = mean(mean(postCue_lickAlignCSPevents_acrossTrials(preLickRange(1)-1+peakIdx_postCueLickActCSP_allN-statsWindow:preLickRange(1)-1+peakIdx_postCueLickActCSP_allN+statsWindow,:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickAlign.CSP_postCue_trialAvg_Cspk = std(mean(postCue_lickAlignCSPevents_acrossTrials(preLickRange(1)-1+peakIdx_postCueLickActCSP_allN-statsWindow:preLickRange(1)-1+peakIdx_postCueLickActCSP_allN+statsWindow,:),1,'omitnan'),[],2,'omitnan')./sqrt(size(postCue_lickAlignCSPevents_acrossTrials,2)).*(1000./frameRateHz);
for c = 1:sum(nsize(2,:))
postCueLick_baseline_CSM(:,c) = CSMevents_acrossTrials(baselineRange(1)-1+peakIdx_baseCSM_allN-statsWindow:baselineRange(1)-1+peakIdx_baseCSM_allN+statsWindow,c);
postCueLick_activity_CSM(:,c) = postCue_lickAlignCSMevents_acrossTrials(preLickRange(1)-1+peakIdx_postCueLickActCSM_allN-statsWindow:preLickRange(1)-1+peakIdx_postCueLickActCSM_allN+statsWindow,c);
end
[~, p.lickAlign.CSM_postCue_trialAvg_Cspk, ~, stats.lickAlign.CSM_postCue_trialAvg_Cspk] = ttest(mean(postCueLick_baseline_CSM,1,'omitnan'),mean(postCueLick_activity_CSM,1,'omitnan'),'tail','left');
  m.lickAlign.CSM_postCue_trialAvg_Cspk = mean(mean(postCue_lickAlignCSMevents_acrossTrials(preLickRange(1)-1+peakIdx_postCueLickActCSM_allN-statsWindow:preLickRange(1)-1+peakIdx_postCueLickActCSM_allN+statsWindow,:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickAlign.CSM_postCue_trialAvg_Cspk = std(mean(postCue_lickAlignCSMevents_acrossTrials(preLickRange(1)-1+peakIdx_postCueLickActCSM_allN-statsWindow:preLickRange(1)-1+peakIdx_postCueLickActCSM_allN+statsWindow,:),1,'omitnan'),[],2,'omitnan')./sqrt(size(postCue_lickAlignCSMevents_acrossTrials,2)).*(1000./frameRateHz);
%% plot     | statistical window for first-lick avg vs. baseline avg
avg_firstLickPaneCSP(1,1) = round(mean(peakIdx_firstLickActCSP_allN,2,'omitnan')); avg_firstLickPaneCSM(1,2) = round(mean(peakIdx_firstLickActCSM_allN,2,'omitnan'));

  tempFig=setFigure; hold on; xline(0,'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from lick onset (s)')
subplot(2,1,1); hold on;
xlim([-500 1000]); ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('Average Cspk rate preceeding first lick in CS+ trials with statistical windows');
rectangle('Position',[tt(1,35+preLickRange(1)-1+avg_firstLickPaneCSP(1,1)-statsWindow),-1+.025,baselineRange(1)-1+peak_baselinePane(1,1)+statsWindow,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(first_lickAlignCSPevents_acrossTrials,3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(first_lickAlignCSPevents_acrossTrials,3,'omitnan'),[],2,'omitnan')./sqrt(size(first_lickAlignCSPevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
subplot(2,1,2); hold on;
xlim([-500 1000]); ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('Average Cspk rate preceeding first lick in CS- trials with statistical windows');
rectangle('Position',[tt(1,35+preLickRange(1)-1+avg_firstLickPaneCSM(1,2)-statsWindow),-1+.025,baselineRange(1)-1+peak_baselinePane(1,1)+statsWindow,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(first_lickAlignCSMevents_acrossTrials,3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(first_lickAlignCSMevents_acrossTrials,3,'omitnan'),[],2,'omitnan')./sqrt(size(first_lickAlignCSMevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
saveas(tempFig, [output_fn 'stats\lickWindow_' blockName '.pdf'],'pdf');
%% plot     | first-piezo (response/suppression) 

  %we're gonna set an eval loop here because individual blocks for first-postCue-postRew for CS+&- is vile
eventLabel={'first','postCue','postRew'}; cueLabel={'CSP','CSM'}; prePiezo_frames=15; prePiezoRange=(1:bx_response_frames)+(prePiezo_frames-bx_response_frames);
  %CS+
  for event=1:size(eventLabel,2)
      eval(['cell.piezoAlign.response_' eventLabel{1,event} 'CSP=[];']);
      eval(['cell.piezoAlign.suppress_' eventLabel{1,event} 'CSP=[];']);
    for c = 1:sum(nsize(1,:))
%
 eval(['prePiezoVar = ' eventLabel{1,event} '_piezoAlignDFoFCSP_acrossTrials(prePiezoRange,c);']);
 baselineVar = CSPdFoF_acrossTrials(baselineRange,c);
 eval(['pos.Diff.piezoAlign{c,1} = ([0; diff(' eventLabel{1,event} '_piezoAlignDFoFCSP_acrossTrials(:,c))]);']);
 
pos.Limit.piezoAlign{c,1}.one = (mean(pos.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)>0])),1,'omitnan')+(std(pos.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)>0])),[],1,'omitnan').*1)); 
pos.Change.piezoAlign{c,1}.one = find(pos.Diff.piezoAlign{c,1}>pos.Limit.piezoAlign{c,1}.one); 
pos.Limit.piezoAlign{c,1}.two = (mean(pos.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)>0])),1,'omitnan')+(std(pos.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)>0])),[],1,'omitnan').*1.5)); 
pos.Change.piezoAlign{c,1}.two = find(pos.Diff.piezoAlign{c,1}>pos.Limit.piezoAlign{c,1}.two); 
pos.Limit.piezoAlign{c,1}.three = (mean(pos.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)>0])),1,'omitnan')+(std(pos.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)>0])),[],1,'omitnan').*2)); 
pos.Change.piezoAlign{c,1}.three = find(pos.Diff.piezoAlign{c,1}>pos.Limit.piezoAlign{c,1}.three);

eval(['frameIdx.piezoAlign.' eventLabel{1,event} 'CSP_posChangeDiff{1,c} = intersect(prePiezoRange, pos.Change.piezoAlign{c,1}.one);']); 
eval(['frameIdx.piezoAlign.' eventLabel{1,event} 'CSP_posChangeDiff{2,c} = intersect(prePiezoRange, pos.Change.piezoAlign{c,1}.two);']); 
eval(['frameIdx.piezoAlign.' eventLabel{1,event} 'CSP_posChangeDiff{3,c} = intersect(prePiezoRange, pos.Change.piezoAlign{c,1}.three);']);

 eval(['neg.Diff.piezoAlign{c,1} = ([0; diff(' eventLabel{1,event} '_piezoAlignDFoFCSP_acrossTrials(:,c))]);']);
neg.Limit.piezoAlign{c,1}.one = (mean(neg.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)<0])),1,'omitnan')-(std(neg.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)<0])),[],1,'omitnan').*1)); 
neg.Change.piezoAlign{c,1}.one = find(neg.Diff.piezoAlign{c,1}<neg.Limit.piezoAlign{c,1}.one); 
neg.Limit.piezoAlign{c,1}.two = (mean(neg.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)<0])),1,'omitnan')-(std(neg.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)<0])),[],1,'omitnan').*1.5)); 
neg.Change.piezoAlign{c,1}.two = find(neg.Diff.piezoAlign{c,1}<neg.Limit.piezoAlign{c,1}.two); 
neg.Limit.piezoAlign{c,1}.three = (mean(neg.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)<0])),1,'omitnan')-(std(neg.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)<0])),[],1,'omitnan').*2)); 
neg.Change.piezoAlign{c,1}.three = find(neg.Diff.piezoAlign{c,1}<neg.Limit.piezoAlign{c,1}.three);

eval(['frameIdx.piezoAlign.' eventLabel{1,event} 'CSP_negChangeDiff{1,c} = intersect(prePiezoRange, neg.Change.piezoAlign{c,1}.one);']); 
eval(['frameIdx.piezoAlign.' eventLabel{1,event} 'CSP_negChangeDiff{2,c} = intersect(prePiezoRange, neg.Change.piezoAlign{c,1}.two);']); 
eval(['frameIdx.piezoAlign.' eventLabel{1,event} 'CSP_negChangeDiff{3,c} = intersect(prePiezoRange, neg.Change.piezoAlign{c,1}.three);']);
%    
    tempResp_diffIdx = diff(intersect(prePiezoRange, find(prePiezoVar>(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*respReq_piezoAlign))))';
    tempResp_frameIdx = intersect(prePiezoRange, find(prePiezoVar>(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*respReq_piezoAlign)));
    if length(tempResp_diffIdx)
      if strfind(tempResp_diffIdx, [1 1])
        if eval(['frameIdx.piezoAlign.' eventLabel{1,event} 'CSP_posChangeDiff{chngReq_piezoAlign,c}'])
           eval(['cell.piezoAlign.response_' eventLabel{1,event} 'CSP=[cell.piezoAlign.response_' eventLabel{1,event} 'CSP c];']);
        end
      end
    end
    tempSupp_diffIdx = diff(intersect(prePiezoRange, find(prePiezoVar<(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*suppReq_piezoAlign))))';
    tempSupp_frameIdx = intersect(prePiezoRange, find(prePiezoVar<(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*suppReq_piezoAlign)));
    if length(tempSupp_diffIdx)
      if strfind(tempSupp_diffIdx, [1 1 1 1]) 
        if isempty(strfind(tempResp_diffIdx, [1 1]))
          if eval(['frameIdx.piezoAlign.' eventLabel{1,event} 'CSP_negChangeDiff{chngReq_piezoAlign,c}'])
           eval(['cell.piezoAlign.suppress_' eventLabel{1,event} 'CSP=[cell.piezoAlign.suppress_' eventLabel{1,event} 'CSP c];']);
          end
        end
      end
    end
    end
  end
  
  %CS-
   for event=1:size(eventLabel,2)
      eval(['cell.piezoAlign.response_' eventLabel{1,event} 'CSM=[];']);
      eval(['cell.piezoAlign.suppress_' eventLabel{1,event} 'CSM=[];']);
    for c = 1:sum(nsize(2,:))
%
 eval(['prePiezoVar = ' eventLabel{1,event} '_piezoAlignDFoFCSM_acrossTrials(prePiezoRange,c);']);
 baselineVar = CSMdFoF_acrossTrials(baselineRange,c);
 eval(['pos.Diff.piezoAlign{c,1} = ([0; diff(' eventLabel{1,event} '_piezoAlignDFoFCSM_acrossTrials(:,c))]);']);
 
pos.Limit.piezoAlign{c,1}.one = (mean(pos.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)>0])),1,'omitnan')+(std(pos.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)>0])),[],1,'omitnan').*1)); 
pos.Change.piezoAlign{c,1}.one = find(pos.Diff.piezoAlign{c,1}>pos.Limit.piezoAlign{c,1}.one); 
pos.Limit.piezoAlign{c,1}.two = (mean(pos.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)>0])),1,'omitnan')+(std(pos.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)>0])),[],1,'omitnan').*1.5)); 
pos.Change.piezoAlign{c,1}.two = find(pos.Diff.piezoAlign{c,1}>pos.Limit.piezoAlign{c,1}.two); 
pos.Limit.piezoAlign{c,1}.three = (mean(pos.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)>0])),1,'omitnan')+(std(pos.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)>0])),[],1,'omitnan').*2)); 
pos.Change.piezoAlign{c,1}.three = find(pos.Diff.piezoAlign{c,1}>pos.Limit.piezoAlign{c,1}.three);

eval(['frameIdx.piezoAlign.' eventLabel{1,event} 'CSM_posChangeDiff{1,c} = intersect(prePiezoRange, pos.Change.piezoAlign{c,1}.one);']); 
eval(['frameIdx.piezoAlign.' eventLabel{1,event} 'CSM_posChangeDiff{2,c} = intersect(prePiezoRange, pos.Change.piezoAlign{c,1}.two);']); 
eval(['frameIdx.piezoAlign.' eventLabel{1,event} 'CSM_posChangeDiff{3,c} = intersect(prePiezoRange, pos.Change.piezoAlign{c,1}.three);']);

 eval(['neg.Diff.piezoAlign{c,1} = ([0; diff(' eventLabel{1,event} '_piezoAlignDFoFCSM_acrossTrials(:,c))]);']);
neg.Limit.piezoAlign{c,1}.one = (mean(neg.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)<0])),1,'omitnan')-(std(neg.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)<0])),[],1,'omitnan').*1)); 
neg.Change.piezoAlign{c,1}.one = find(neg.Diff.piezoAlign{c,1}<neg.Limit.piezoAlign{c,1}.one); 
neg.Limit.piezoAlign{c,1}.two = (mean(neg.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)<0])),1,'omitnan')-(std(neg.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)<0])),[],1,'omitnan').*1.5)); 
neg.Change.piezoAlign{c,1}.two = find(neg.Diff.piezoAlign{c,1}<neg.Limit.piezoAlign{c,1}.two); 
neg.Limit.piezoAlign{c,1}.three = (mean(neg.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)<0])),1,'omitnan')-(std(neg.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)<0])),[],1,'omitnan').*2)); 
neg.Change.piezoAlign{c,1}.three = find(neg.Diff.piezoAlign{c,1}<neg.Limit.piezoAlign{c,1}.three);

eval(['frameIdx.piezoAlign.' eventLabel{1,event} 'CSM_negChangeDiff{1,c} = intersect(prePiezoRange, neg.Change.piezoAlign{c,1}.one);']); 
eval(['frameIdx.piezoAlign.' eventLabel{1,event} 'CSM_negChangeDiff{2,c} = intersect(prePiezoRange, neg.Change.piezoAlign{c,1}.two);']); 
eval(['frameIdx.piezoAlign.' eventLabel{1,event} 'CSM_negChangeDiff{3,c} = intersect(prePiezoRange, neg.Change.piezoAlign{c,1}.three);']);
%    
    tempResp_diffIdx = diff(intersect(prePiezoRange, find(prePiezoVar>(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*respReq_piezoAlign))))';
    tempResp_frameIdx = intersect(prePiezoRange, find(prePiezoVar>(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*respReq_piezoAlign)));
    if length(tempResp_diffIdx)
      if strfind(tempResp_diffIdx, [1 1])
        if eval(['frameIdx.piezoAlign.' eventLabel{1,event} 'CSM_posChangeDiff{chngReq_piezoAlign,c}'])
           eval(['cell.piezoAlign.response_' eventLabel{1,event} 'CSM=[cell.piezoAlign.response_' eventLabel{1,event} 'CSM c];']);
        end
      end
    end
    tempSupp_diffIdx = diff(intersect(prePiezoRange, find(prePiezoVar<(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*suppReq_piezoAlign))))';
    tempSupp_frameIdx = intersect(prePiezoRange, find(prePiezoVar<(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*suppReq_piezoAlign)));
    if length(tempSupp_diffIdx)
      if strfind(tempSupp_diffIdx, [1 1 1 1]) 
        if isempty(strfind(tempResp_diffIdx, [1 1]))
          if eval(['frameIdx.piezoAlign.' eventLabel{1,event} 'CSM_negChangeDiff{chngReq_piezoAlign,c}'])
           eval(['cell.piezoAlign.suppress_' eventLabel{1,event} 'CSM=[cell.piezoAlign.suppress_' eventLabel{1,event} 'CSM c];']);
          end
        end
      end
    end
    end
  end


it=1; sideLabel={'response','suppress';'suppress','response'};
for cue=1:size(cueLabel,2)
  for event=1:size(eventLabel,2)
    for side=1:size(sideLabel,2)
eval(['cellDiff.piezoAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} '= setdiff(cell.piezoAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} ',cell.piezoAlign.' sideLabel{2,side} '_'  eventLabel{1,event} cueLabel{1,cue} ');']);
it=it+1;
    end
  end
end
%
titleArray = {'CS+ reponsive firstPiezo','CS+ suppressive firstPiezo','CS+ reponsive postCuePiezo','CS+ suppressive postCuePiezo','CS+ reponsive postRewPiezo','CS+ suppressive postRewPiezo','CS- reponsive firstPiezo','CS- suppressive firstPiezo','CS- reponsive postCuePiezo','CS- suppressive postCuePiezo','CS- reponsive postRewPiezo','CS- suppressive postRewPiezo'};
colorArray = {'k','b','k','b','k','b','r','m','r','m','r','m'};

tempFig = setFigure; hold on; c=1;
for cue=1:size(cueLabel,2)
  for event=1:size(eventLabel,2)
    for side=1:size(sideLabel,2)
    eval(['callTrue = ' eventLabel{1,event} '_piezoAlignDFoF' cueLabel{1,cue} '_acrossTrials(:,cellDiff.piezoAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} ');']);
    %eval(['callFalse = ' eventLabel{1,event} '_piezoAlignDFoF' cueLabel{1,cue} '_acrossTrials(:,setdiff(sum(nsize(cue,:)),cellDiff.piezoAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} '));']);
    eval(['dendrites = cellDiff.piezoAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} ';']);
subplot(6,2,c); hold on;
xlim([-500 1000]); ylim([ymin*2 ymax*2]); xline(0,'k');  title(titleArray{1,c});
rectangle('Position',[tt(1,baselineRange(1))+10,ymin*2+.25,(1000./frameRateHz)*response_frames-10,ymax+5],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,37:81)/1000,mean(callTrue,2,'omitnan').*(1000./frameRateHz),std(callTrue,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',colorArray{1,c});
%shadedErrorBar_CV(tt(1,37:81)/1000,mean(callFalse,2,'omitnan').*(1000./frameRateHz),std(callFalse,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',{'Color',[0/255 0/255 201/255]});
    c=c+1;
    end
  end
end
    saveas(tempFig, [output_fn 'stats\dFoF_piezoAlign_PSTHs_' blockName '.pdf'],'pdf');

tempFig = setFigure; hold on; c=1;
for cue=1:size(cueLabel,2)
  for event=1:size(eventLabel,2)
    for side=1:size(sideLabel,2)
    eval(['callTrue = ' eventLabel{1,event} '_piezoAlign' cueLabel{1,cue} 'events_acrossTrials(:,cellDiff.piezoAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} ');']);
    %eval(['callFalse = ' eventLabel{1,event} '_piezoAlign' cueLabel{1,cue} 'events_acrossTrials(:,setdiff(sum(nsize(cue,:)),cellDiff.piezoAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} '));']);
    eval(['dendrites = cellDiff.piezoAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} ';']);
subplot(6,2,c); hold on;
xlim([-500 1000]); ylim([0 ymax]); xline(0,'k');  title(titleArray{1,c});
if c<5; rectangle('Position',[tt(1,baselineRange(1)),0+.15,(1000./frameRateHz)*response_frames,ymax+1],'EdgeColor',highlight,'FaceColor',highlight); elseif c>4; rectangle('Position',[tt(1,baselineRange(1)),0+.15,(1000./frameRateHz)*response_frames,ymax+1],'EdgeColor',highlight,'FaceColor',highlight); end
shadedErrorBar_CV(tt(1,37:81)/1000,mean(callTrue,2,'omitnan').*(1000./frameRateHz),std(callTrue,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',colorArray{1,c});
%shadedErrorBar_CV(tt(1,37:81)/1000,mean(callFalse,2,'omitnan').*(1000./frameRateHz),std(callFalse,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',{'Color',[0/255 0/255 201/255]});
    c=c+1;
    end
  end
end
    saveas(tempFig, [output_fn 'stats\CspkEvents_piezoAlign_PSTHs_' blockName '.pdf'],'pdf');

    tempFig = setFigure; hold on; c=1;
for event=1:size(eventLabel,2)
  for cue=1:size(cueLabel,2)
eval(['respN = size(' eventLabel{1,event} '_piezoAlign' cueLabel{1,cue} 'events_acrossTrials(:,cellDiff.piezoAlign.response_' eventLabel{1,event} cueLabel{1,cue} '),2);']);
eval(['suppN = size(' eventLabel{1,event} '_piezoAlign' cueLabel{1,cue} 'events_acrossTrials(:,cellDiff.piezoAlign.suppress_' eventLabel{1,event} cueLabel{1,cue} '),2);']);
       elseN = sum(nsize(cue,:))-(respN+suppN);

subplot(3,2,c); 
pie([respN, suppN, elseN]);
c=c+1; title([eventLabel{1,event} '-' cueLabel{1,cue}]);
    
  end
end
    saveas(tempFig, [output_fn 'stats\pie_propCells_piezoAlign_' blockName '.pdf'],'pdf');
    
tempFig=setFigure; hold on;
subplot(3,1,1); hold on;
shadedErrorBar_CV(tt(1,37:81)/1000,mean(first_piezoAlignCSPevents_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(first_piezoAlignCSPevents_acrossTrials,[],2,'omitnan')./sqrt(size(first_piezoAlignCSPevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,37:81)/1000,mean(first_piezoAlignCSMevents_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(first_piezoAlignCSMevents_acrossTrials,[],2,'omitnan')./sqrt(size(first_piezoAlignCSMevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
xlim([-500 1000]); ylim([0 ymax]); xline(0,'k');  title('first piezo');
subplot(3,1,2); hold on;
shadedErrorBar_CV(tt(1,37:81)/1000,mean(postCue_piezoAlignCSPevents_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(postCue_piezoAlignCSPevents_acrossTrials,[],2,'omitnan')./sqrt(size(postCue_piezoAlignCSPevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,37:81)/1000,mean(postCue_piezoAlignCSMevents_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(postCue_piezoAlignCSMevents_acrossTrials,[],2,'omitnan')./sqrt(size(postCue_piezoAlignCSMevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
xlim([-500 1000]); ylim([0 ymax]); xline(0,'k');  title('first postCue piezo');
subplot(3,1,3); hold on;
shadedErrorBar_CV(tt(1,37:81)/1000,mean(postRew_piezoAlignCSPevents_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(postRew_piezoAlignCSPevents_acrossTrials,[],2,'omitnan')./sqrt(size(postRew_piezoAlignCSPevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,37:81)/1000,mean(postRew_piezoAlignCSMevents_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(postRew_piezoAlignCSMevents_acrossTrials,[],2,'omitnan')./sqrt(size(postRew_piezoAlignCSMevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
xlim([-500 1000]); ylim([0 ymax]); xline(0,'k');  title('first postRew piezo');
%% ttest    | first-piezo avg vs. baseline avg 
% offset start frame by 1, because we need to take a +/-1 frame window
% around the target; script breaks if the peak is on frame 1
  [peak_firstPiezoActCSP, peakIdx_firstPiezoActCSP] = max(first_piezoAlignCSPevents_acrossTrials(prePiezoRange,:),[],1);
  [peak_postCuePiezoActCSP, peakIdx_postCuePiezoActCSP] = max(postCue_piezoAlignCSPevents_acrossTrials(prePiezoRange,:),[],1);
  [peak_postRewPiezoActCSP, peakIdx_postRewPiezoActCSP] = max(postRew_piezoAlignCSPevents_acrossTrials(prePiezoRange,:),[],1);
  [peak_baseCSP, peakIdx_baseCSP] = max(CSPevents_acrossTrials(baselineRange,:),[],1);
  [peak_firstPiezoActCSM, peakIdx_firstPiezoActCSM] = max(first_piezoAlignCSMevents_acrossTrials(prePiezoRange,:),[],1);
  [peak_postCuePiezoActCSM, peakIdx_postCuePiezoActCSM] = max(postCue_piezoAlignCSMevents_acrossTrials(prePiezoRange,:),[],1);
  [peak_postRewPiezoActCSM, peakIdx_postRewPiezoActCSM] = max(postRew_piezoAlignCSMevents_acrossTrials(prePiezoRange,:),[],1);
  [peak_baseCSM, peakIdx_baseCSM] = max(CSMevents_acrossTrials(baselineRange,:),[],1);

  [~, peakIdx_firstPiezoActCSP_allN] = max(mean(first_piezoAlignCSPevents_acrossTrials(prePiezoRange,:),2,'omitnan'),[],1);
  [~, peakIdx_postCuePiezoActCSP_allN] = max(mean(postCue_piezoAlignCSPevents_acrossTrials(prePiezoRange,:),2,'omitnan'),[],1);
  [~, peakIdx_postRewPiezoActCSP_allN] = max(mean(postRew_piezoAlignCSPevents_acrossTrials(prePiezoRange,:),2,'omitnan'),[],1);
  [~, peakIdx_baseCSP_allN] = max(mean(CSPevents_acrossTrials(baselineRange,:),2,'omitnan'),[],1);
  [~, peakIdx_firstPiezoActCSM_allN] = max(mean(first_piezoAlignCSMevents_acrossTrials(prePiezoRange,:),2,'omitnan'),[],1);
  [~, peakIdx_postCuePiezoActCSM_allN] = max(mean(postCue_piezoAlignCSMevents_acrossTrials(prePiezoRange,:),2,'omitnan'),[],1);
  [~, peakIdx_postRewPiezoActCSM_allN] = max(mean(postRew_piezoAlignCSMevents_acrossTrials(prePiezoRange,:),2,'omitnan'),[],1);
  [~, peakIdx_baseCSM_allN] = max(mean(CSMevents_acrossTrials(baselineRange,:),2,'omitnan'),[],1);
% counts for piezo align
n.piezoAlign.CSP_first_trialAvg_count = sum(first_piezoAlignCSP_acrossTrials); 
n.piezoAlign.CSM_first_trialAvg_count = sum(first_piezoAlignCSM_acrossTrials);
n.piezoAlign.CSP_firstPostCue_trialAvg_count = sum(firstPostCue_piezoAlignCSP_acrossTrials);
n.piezoAlign.CSM_firstPostCue_trialAvg_count = sum(firstPostCue_piezoAlignCSM_acrossTrials);
n.piezoAlign.CSP_firstPostRew_trialAvg_count = sum(firstPostRew_piezoAlignCSP_acrossTrials);
n.piezoAlign.CSM_firstPostRew_trialAvg_count = sum(firstPostRew_piezoAlignCSM_acrossTrials);
for c = 1:sum(nsize(1,:))
firstPiezo_baseline_CSP(:,c) = CSPevents_acrossTrials(baselineRange(1)-1+peakIdx_baseCSP_allN-statsWindow:baselineRange(1)-1+peakIdx_baseCSP_allN+statsWindow,c);
firstPiezo_activity_CSP(:,c) = first_piezoAlignCSPevents_acrossTrials(prePiezoRange(1)-1+peakIdx_firstPiezoActCSP_allN-statsWindow:prePiezoRange(1)-1+peakIdx_firstPiezoActCSP_allN+statsWindow,c);
end
[~, p.piezoAlign.CSP_first_trialAvg_Cspk, ~, stats.piezoAlign.CSP_first_trialAvg_Cspk] = ttest(mean(firstPiezo_baseline_CSP,1,'omitnan'),mean(firstPiezo_activity_CSP,1,'omitnan'),'tail','left');
  m.piezoAlign.CSP_first_trialAvg_Cspk = mean(mean(first_piezoAlignCSPevents_acrossTrials(prePiezoRange(1)-1+peakIdx_firstPiezoActCSP_allN-statsWindow:prePiezoRange(1)-1+peakIdx_firstPiezoActCSP_allN+statsWindow,:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoAlign.CSP_first_trialAvg_Cspk = std(mean(first_piezoAlignCSPevents_acrossTrials(prePiezoRange(1)-1+peakIdx_firstPiezoActCSP_allN-statsWindow:prePiezoRange(1)-1+peakIdx_firstPiezoActCSP_allN+statsWindow,:),1,'omitnan'),[],2,'omitnan')./sqrt(size(first_piezoAlignCSMevents_acrossTrials,2)).*(1000./frameRateHz);
for c = 1:sum(nsize(2,:))
firstPiezo_baseline_CSM(:,c) = CSMevents_acrossTrials(baselineRange(1)-1+peakIdx_baseCSM_allN-statsWindow:baselineRange(1)-1+peakIdx_baseCSM_allN+statsWindow,c);
firstPiezo_activity_CSM(:,c) = first_piezoAlignCSMevents_acrossTrials(prePiezoRange(1)-1+peakIdx_firstPiezoActCSM_allN-statsWindow:prePiezoRange(1)-1+peakIdx_firstPiezoActCSM_allN+statsWindow,c);
end
[~, p.piezoAlign.CSM_first_trialAvg_Cspk, ~, stats.piezoAlign.CSM_first_trialAvg_Cspk] = ttest(mean(firstPiezo_baseline_CSM,1,'omitnan'),mean(firstPiezo_activity_CSM,1,'omitnan'),'tail','left');
  m.piezoAlign.CSM_first_trialAvg_Cspk = mean(mean(first_piezoAlignCSMevents_acrossTrials(prePiezoRange(1)-1+peakIdx_firstPiezoActCSM_allN-statsWindow:prePiezoRange(1)-1+peakIdx_firstPiezoActCSM_allN+statsWindow,:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoAlign.CSM_first_trialAvg_Cspk = std(mean(first_piezoAlignCSMevents_acrossTrials(prePiezoRange(1)-1+peakIdx_firstPiezoActCSM_allN-statsWindow:prePiezoRange(1)-1+peakIdx_firstPiezoActCSM_allN+statsWindow,:),1,'omitnan'),[],2,'omitnan')./sqrt(size(first_piezoAlignCSMevents_acrossTrials,2)).*(1000./frameRateHz);

% difference in activity between first piezo [POST REWARD] and baseline
for c = 1:sum(nsize(1,:))
postRewPiezo_baseline_CSP(:,c) = CSPevents_acrossTrials(baselineRange(1)-1+peakIdx_baseCSP_allN-statsWindow:baselineRange(1)-1+peakIdx_baseCSP_allN+statsWindow,c);
postRewPiezo_activity_CSP(:,c) = postRew_piezoAlignCSPevents_acrossTrials(prePiezoRange(1)-1+peakIdx_postRewPiezoActCSP_allN-statsWindow:prePiezoRange(1)-1+peakIdx_postRewPiezoActCSP_allN+statsWindow,c);
end
[~, p.piezoAlign.CSP_postRew_trialAvg_Cspk, ~, stats.piezoAlign.CSP_postRew_trialAvg_Cspk] = ttest(mean(postRewPiezo_baseline_CSP,1,'omitnan'),mean(postRewPiezo_activity_CSP,1,'omitnan'),'tail','left');
  m.piezoAlign.CSP_postRew_trialAvg_Cspk = mean(mean(postRew_piezoAlignCSPevents_acrossTrials(prePiezoRange(1)-1+peakIdx_postRewPiezoActCSP_allN-statsWindow:prePiezoRange(1)-1+peakIdx_postRewPiezoActCSP_allN+statsWindow,:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoAlign.CSP_postRew_trialAvg_Cspk = std(mean(postRew_piezoAlignCSPevents_acrossTrials(prePiezoRange(1)-1+peakIdx_postRewPiezoActCSP_allN-statsWindow:prePiezoRange(1)-1+peakIdx_postRewPiezoActCSP_allN+statsWindow,:),1,'omitnan'),[],2,'omitnan')./sqrt(size(postRew_piezoAlignCSPevents_acrossTrials,2)).*(1000./frameRateHz);
for c = 1:sum(nsize(2,:))
postRewPiezo_baseline_CSM(:,c) = CSMevents_acrossTrials(baselineRange(1)-1+peakIdx_baseCSM_allN-statsWindow:baselineRange(1)-1+peakIdx_baseCSM_allN+statsWindow,c);
postRewPiezo_activity_CSM(:,c) = postRew_piezoAlignCSMevents_acrossTrials(prePiezoRange(1)-1+peakIdx_postRewPiezoActCSM_allN-statsWindow:prePiezoRange(1)-1+peakIdx_postRewPiezoActCSM_allN+statsWindow,c);
end
[~, p.piezoAlign.CSM_postRew_trialAvg_Cspk, ~, stats.piezoAlign.CSM_postRew_trialAvg_Cspk] = ttest(mean(postRewPiezo_baseline_CSM,1,'omitnan'),mean(postRewPiezo_activity_CSM,1,'omitnan'),'tail','left');
  m.piezoAlign.CSM_postRew_trialAvg_Cspk = mean(mean(postRew_piezoAlignCSMevents_acrossTrials(prePiezoRange(1)-1+peakIdx_postRewPiezoActCSM_allN-statsWindow:prePiezoRange(1)-1+peakIdx_postRewPiezoActCSM_allN+statsWindow,:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoAlign.CSM_postRew_trialAvg_Cspk = std(mean(postRew_piezoAlignCSMevents_acrossTrials(prePiezoRange(1)-1+peakIdx_postRewPiezoActCSM_allN-statsWindow:prePiezoRange(1)-1+peakIdx_postRewPiezoActCSM_allN+statsWindow,:),1,'omitnan'),[],2,'omitnan')./sqrt(size(postRew_piezoAlignCSMevents_acrossTrials,2)).*(1000./frameRateHz);
  
% difference in activity between first piezo [POST CUE] and baseline
for c = 1:sum(nsize(1,:))
postCuePiezo_baseline_CSP(:,c) = CSPevents_acrossTrials(baselineRange(1)-1+peakIdx_baseCSP_allN-statsWindow:baselineRange(1)-1+peakIdx_baseCSP_allN+statsWindow,c);
postCuePiezo_activity_CSP(:,c) = postCue_piezoAlignCSPevents_acrossTrials(prePiezoRange(1)-1+peakIdx_postCuePiezoActCSP_allN-statsWindow:prePiezoRange(1)-1+peakIdx_postCuePiezoActCSP_allN+statsWindow,c);
end
[~, p.piezoAlign.CSP_postCue_trialAvg_Cspk, ~, stats.piezoAlign.CSP_postCue_trialAvg_Cspk] = ttest(mean(postCuePiezo_baseline_CSP,1,'omitnan'),mean(postCuePiezo_activity_CSP,1,'omitnan'),'tail','left');
  m.piezoAlign.CSP_postCue_trialAvg_Cspk = mean(mean(postCue_piezoAlignCSPevents_acrossTrials(prePiezoRange(1)-1+peakIdx_postCuePiezoActCSP_allN-statsWindow:prePiezoRange(1)-1+peakIdx_postCuePiezoActCSP_allN+statsWindow,:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoAlign.CSP_postCue_trialAvg_Cspk = std(mean(postCue_piezoAlignCSPevents_acrossTrials(prePiezoRange(1)-1+peakIdx_postCuePiezoActCSP_allN-statsWindow:prePiezoRange(1)-1+peakIdx_postCuePiezoActCSP_allN+statsWindow,:),1,'omitnan'),[],2,'omitnan')./sqrt(size(postCue_piezoAlignCSPevents_acrossTrials,2)).*(1000./frameRateHz);
for c = 1:sum(nsize(2,:))
postCuePiezo_baseline_CSM(:,c) = CSMevents_acrossTrials(baselineRange(1)-1+peakIdx_baseCSM_allN-statsWindow:baselineRange(1)-1+peakIdx_baseCSM_allN+statsWindow,c);
postCuePiezo_activity_CSM(:,c) = postCue_piezoAlignCSMevents_acrossTrials(prePiezoRange(1)-1+peakIdx_postCuePiezoActCSM_allN-statsWindow:prePiezoRange(1)-1+peakIdx_postCuePiezoActCSM_allN+statsWindow,c);
end
[~, p.piezoAlign.CSM_postCue_trialAvg_Cspk, ~, stats.piezoAlign.CSM_postCue_trialAvg_Cspk] = ttest(mean(postCuePiezo_baseline_CSM,1,'omitnan'),mean(postCuePiezo_activity_CSM,1,'omitnan'),'tail','left');
  m.piezoAlign.CSM_postCue_trialAvg_Cspk = mean(mean(postCue_piezoAlignCSMevents_acrossTrials(prePiezoRange(1)-1+peakIdx_postCuePiezoActCSM_allN-statsWindow:prePiezoRange(1)-1+peakIdx_postCuePiezoActCSM_allN+statsWindow,:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoAlign.CSM_postCue_trialAvg_Cspk = std(mean(postCue_piezoAlignCSMevents_acrossTrials(prePiezoRange(1)-1+peakIdx_postCuePiezoActCSM_allN-statsWindow:prePiezoRange(1)-1+peakIdx_postCuePiezoActCSM_allN+statsWindow,:),1,'omitnan'),[],2,'omitnan')./sqrt(size(postCue_piezoAlignCSMevents_acrossTrials,2)).*(1000./frameRateHz);
%% plot     | statistical window for first-piezo avg vs. baseline avg
mouse(1,1) = round(mean(peakIdx_firstPiezoActCSP_allN,2,'omitnan')); avg_firstPiezoPane(1,2) = round(mean(peakIdx_firstPiezoActCSM_allN,2,'omitnan'));

  tempFig=setFigure; hold on; xline(0,'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from movement onset (s)')
subplot(2,1,1); hold on;
xlim([-500 1000]); ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('Average Cspk rate preceeding first movement event in CS+ trials with statistical windows');
rectangle('Position',[tt(1,35+preLickRange(1)-1+avg_firstPiezoPane(1,1)-statsWindow),-1+.025,baselineRange(1)-1+peak_baselinePane(1,1)+statsWindow,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(plusData.firstPiezoAlignEvents(:,:,:),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(plusData.firstPiezoAlignEvents(:,:,:),3,'omitnan'),[],2,'omitnan')./sqrt(size(plusData.firstPiezoAlignEvents,2)).*(1000./frameRateHz),'lineProps','k');
subplot(2,1,2); hold on;
xlim([-500 1000]); ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('Average Cspk rate preceeding first movement event in CS- trials with statistical windows');
rectangle('Position',[tt(1,35+preLickRange(1)-1+avg_firstPiezoPane(1,2)-statsWindow),-1+.025,baselineRange(1)-1+peak_baselinePane(1,1)+statsWindow,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(minusData.firstPiezoAlignEvents(:,:,:),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(minusData.firstPiezoAlignEvents(:,:,:),3,'omitnan'),[],2,'omitnan')./sqrt(size(minusData.firstPiezoAlignEvents,2)).*(1000./frameRateHz),'lineProps','r');
saveas(tempFig, [output_fn 'stats\peizoWindow_' blockName '.pdf'],'pdf');
%% data org | most vs least motion 
mostMovementAlignedEvents=[]; leastMovementAlignedEvents=[]; mostTrialMovement=[]; leastTrialMovement=[];
mostCSPMovementAlignedEvents=[]; leastCSPMovementAlignedEvents=[]; mostCSPTrialMovement=[]; leastCSPTrialMovement=[];
mostCSMMovementAlignedEvents=[]; leastCSMMovementAlignedEvents=[]; mostCSMTrialMovement=[]; leastCSMTrialMovement=[];
mostCSPTrialMovement_ind{1,mouse} = []; leastCSPTrialMovement_ind{1,mouse} = [];
mostCSMTrialMovement_ind{1,mouse} = []; leastCSMTrialMovement_ind{1,mouse} = [];
%CS+
for mouse = 1:plusData.exptCount
    CSPTrialMovement_ind=[]; rewardTrials{mouse}=[];
    peakTrialChange=[]; preakTrialChange_ind=[]; trialChange=[];
for c = 1:size(plusData.groupAlign.piezo{mouse},2) %trials
for ii = 1:size(plusData.groupAlign.piezo{mouse},1)-1 %frames
    trialChange(ii,c) = plusData.groupAlign.piezo{mouse}(ii+1,c)-plusData.groupAlign.piezo{mouse}(ii,c);
end
[peakTrialChange(1,c),peakTrialChange_ind(1,c)] = max(trialChange(50:96,c));
end
[TrialMovement, TrialMovement_ind] = sort(peakTrialChange);

for c = 1:size(plusData.block2{1,mouse},2)
   if plusData.block2{1,mouse}(1,c)==0; rewardTrials{mouse} = [rewardTrials{mouse} c]; end
end
for c = 1:size(TrialMovement_ind,2)
    if ~isempty(find(rewardTrials{1,mouse} == TrialMovement_ind(1,c)))
CSPTrialMovement_ind = [CSPTrialMovement_ind TrialMovement_ind(1,c)];
    end
end
mostTrialMovement_ind{1,mouse} = TrialMovement_ind(:,size(plusData.groupAlign.piezo{1,mouse},2)-(floor(size(plusData.groupAlign.piezo{1,mouse},2)./4)):end);
leastTrialMovement_ind{1,mouse} = TrialMovement_ind(:,1:round(size(plusData.groupAlign.piezo{1,mouse},2)./4));

mostCSPTrialMovement_ind{1,mouse} = CSPTrialMovement_ind(:,size(plusData.groupAlign.piezo{1,mouse}(:,~plusData.block2{1,mouse}),2)-(floor(size(plusData.groupAlign.piezo{1,mouse}(:,~plusData.block2{1,mouse}),2)./4)):size(plusData.groupAlign.piezo{1,mouse}(:,~plusData.block2{1,mouse}),2));
leastCSPTrialMovement_ind{1,mouse} = CSPTrialMovement_ind(:,1:ceil(size(plusData.groupAlign.piezo{1,mouse}(:,~plusData.block2{1,mouse}),2)./4));
mostCSPMovementAlignedEvents = [mostCSPMovementAlignedEvents nanmean(plusData.groupAlign.events{1,mouse}(:,:,mostCSPTrialMovement_ind{1,mouse}),3)];
leastCSPMovementAlignedEvents = [leastCSPMovementAlignedEvents nanmean(plusData.groupAlign.events{1,mouse}(:,:,leastCSPTrialMovement_ind{1,mouse}),3)];
end

%CS-
for mouse = 1:minusData.exptCount
    CSMTrialMovement_ind=[];  block2Trials{mouse}=[]; 
    peakTrialChange=[]; preakTrialChange_ind=[]; trialChange=[];
for c = 1:size(minusData.groupAlign.piezo{mouse},2) %trials
for ii = 1:size(minusData.groupAlign.piezo{mouse},1)-1 %frames
    trialChange(ii,c) = minusData.groupAlign.piezo{mouse}(ii+1,c)-minusData.groupAlign.piezo{mouse}(ii,c);
end
[peakTrialChange(1,c),peakTrialChange_ind(1,c)] = max(trialChange(50:96,c));
end
[TrialMovement, TrialMovement_ind] = sort(peakTrialChange);

for c = 1:size(minusData.block2{1,mouse},2)
   if minusData.block2{1,mouse}(1,c)==1; block2Trials{mouse} = [block2Trials{mouse} c]; end
end
for c = 1:size(TrialMovement_ind,2)
    if ~isempty(find(block2Trials{1,mouse} == TrialMovement_ind(1,c)))
CSMTrialMovement_ind = [CSMTrialMovement_ind TrialMovement_ind(1,c)];
    end
end
mostTrialMovement_ind{1,mouse} = TrialMovement_ind(:,size(minusData.groupAlign.piezo{1,mouse},2)-(floor(size(minusData.groupAlign.piezo{1,mouse},2)./4)):end);
leastTrialMovement_ind{1,mouse} = TrialMovement_ind(:,1:round(size(minusData.groupAlign.piezo{1,mouse},2)./4));

mostCSMTrialMovement_ind{1,mouse} = CSMTrialMovement_ind(:,size(minusData.groupAlign.piezo{1,mouse}(:,logical(minusData.block2{1,mouse})),2)-(floor(size(minusData.groupAlign.piezo{1,mouse}(:,logical(minusData.block2{1,mouse})),2)./4)):size(minusData.groupAlign.piezo{1,mouse}(:,logical(minusData.block2{1,mouse})),2));
leastCSMTrialMovement_ind{1,mouse} = CSMTrialMovement_ind(:,1:ceil(size(minusData.groupAlign.piezo{1,mouse}(:,logical(minusData.block2{1,mouse})),2)./4));
mostCSMMovementAlignedEvents = [mostCSMMovementAlignedEvents nanmean(minusData.groupAlign.events{1,mouse}(:,:,mostCSMTrialMovement_ind{1,mouse}),3)];
leastCSMMovementAlignedEvents = [leastCSMMovementAlignedEvents nanmean(minusData.groupAlign.events{1,mouse}(:,:,leastCSMTrialMovement_ind{1,mouse}),3)];
end
%% ttest    | most vs least motion avg vs. baseline avg 
  [peak_lowPostCuePiezoActCSP, peakIdx_lowPostCuePiezoActCSP] = max(leastCSPMovementAlignedEvents(postCueRange,:),[],1);
  [peak_lowPostRewPiezoActCSP, peakIdx_lowPostRewPiezoActCSP] = max(leastCSPMovementAlignedEvents(postRewRange,:),[],1);
  [peak_highPostCuePiezoActCSP, peakIdx_highPostCuePiezoActCSP] = max(mostCSPMovementAlignedEvents(postCueRange,:),[],1);
  [peak_highPostRewPiezoActCSP, peakIdx_highPostRewPiezoActCSP] = max(mostCSPMovementAlignedEvents(postRewRange,:),[],1);
  [peak_lowPostCuePiezoActCSM, peakIdx_lowPostCuePiezoActCSM] = max(leastCSMMovementAlignedEvents(postCueRange,:),[],1);
  [peak_lowPostRewPiezoActCSM, peakIdx_lowPostRewPiezoActCSM] = max(leastCSMMovementAlignedEvents(postRewRange,:),[],1);
  [peak_highPostCuePiezoActCSM, peakIdx_highPostCuePiezoActCSM] = max(mostCSMMovementAlignedEvents(postCueRange,:),[],1);
  [peak_highPostRewPiezoActCSM, peakIdx_highPostRewPiezoActCSM] = max(mostCSMMovementAlignedEvents(postRewRange,:),[],1);

  [peak_lowPostCuePiezoActCSP_allN, peakIdx_lowPostCuePiezoActCSP_allN] = max(mean(leastCSPMovementAlignedEvents(postCueRange,:),2,'omitnan'),[],1);
  [peak_lowPostRewPiezoActCSP_allN, peakIdx_lowPostRewPiezoActCSP_allN] = max(mean(leastCSPMovementAlignedEvents(postRewRange,:),2,'omitnan'),[],1);
  [peak_highPostCuePiezoActCSP_allN, peakIdx_highPostCuePiezoActCSP_allN] = max(mean(mostCSPMovementAlignedEvents(postCueRange,:),2,'omitnan'),[],1);
  [peak_highPostRewPiezoActCSP_allN, peakIdx_highPostRewPiezoActCSP_allN] = max(mean(mostCSPMovementAlignedEvents(postRewRange,:),2,'omitnan'),[],1);
  [peak_lowPostCuePiezoActCSM_allN, peakIdx_lowPostCuePiezoActCSM_allN] = max(mean(leastCSMMovementAlignedEvents(postCueRange,:),2,'omitnan'),[],1);
  [peak_lowPostRewPiezoActCSM_allN, peakIdx_lowPostRewPiezoActCSM_allN] = max(mean(leastCSMMovementAlignedEvents(postRewRange,:),2,'omitnan'),[],1);
  [peak_highPostCuePiezoActCSM_allN, peakIdx_highPostCuePiezoActCSM_allN] = max(mean(mostCSMMovementAlignedEvents(postCueRange,:),2,'omitnan'),[],1);
  [peak_highPostRewPiezoActCSM_allN, peakIdx_highPostRewPiezoActCSM_allN] = max(mean(mostCSMMovementAlignedEvents(postRewRange,:),2,'omitnan'),[],1);

% lowVhigh postCue - CS+
for c = 1:sum(nsize(1,:))
low_postCuePiezo_activity_CSP(:,c) = leastCSPMovementAlignedEvents(postCueRange(1)-1+(peakIdx_lowPostCuePiezoActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_lowPostCuePiezoActCSP_allN+statsWindow),c);  
high_postCuePiezo_activity_CSP(:,c) = mostCSPMovementAlignedEvents(postCueRange(1)-1+(peakIdx_highPostCuePiezoActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_highPostCuePiezoActCSP_allN+statsWindow),c);   
end
[~, p.piezoRate.CSP_postCue_trialAvg_Cspk, ~, stats.piezoRate.CSP_postCue_trialAvg_Cspk] = ttest(mean(low_postCuePiezo_activity_CSP,1,'omitnan'),mean(high_postCuePiezo_activity_CSP,1,'omitnan'));
  n.piezoRate.CSP_postCue_lowRate_Cspk = length(cell2mat(leastCSPTrialMovement_ind));
  m.piezoRate.CSP_postCue_lowRate_Cspk = mean(mean(leastCSPMovementAlignedEvents(postCueRange(1)-1+(peakIdx_lowPostCuePiezoActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_lowPostCuePiezoActCSP_allN+statsWindow),:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoRate.CSP_postCue_lowRate_Cspk = std(mean(leastCSPMovementAlignedEvents(postCueRange(1)-1+(peakIdx_lowPostCuePiezoActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_lowPostCuePiezoActCSP_allN+statsWindow),:),1,'omitnan'),[],2,'omitnan')./sqrt(size(leastCSPMovementAlignedEvents,2)).*(1000./frameRateHz);
  n.piezoRate.CSP_postCue_highRate_Cspk = length(cell2mat(mostCSPTrialMovement_ind));
  m.piezoRate.CSP_postCue_highRate_Cspk = mean(mean(mostCSPMovementAlignedEvents(postCueRange(1)-1+(peakIdx_highPostCuePiezoActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_highPostCuePiezoActCSP_allN+statsWindow),:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoRate.CSP_postCue_highRate_Cspk = std(mean(mostCSPMovementAlignedEvents(postCueRange(1)-1+(peakIdx_highPostCuePiezoActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_highPostCuePiezoActCSP_allN+statsWindow),:),1,'omitnan'),[],2,'omitnan')./sqrt(size(mostCSPMovementAlignedEvents,2)).*(1000./frameRateHz);

  tempFig=setFigure; hold on; xline(0,'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(leastCSPMovementAlignedEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(leastCSPMovementAlignedEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(leastCSPMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(mostCSPMovementAlignedEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(mostCSPMovementAlignedEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(mostCSPMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L movement - postCue - CS+'); ylim([0 ymax]);
      saveas(tempFig, [output_fn 'stats\HL_piezoRate_postCue_CS+_PSTHs_' blockName '.pdf'],'pdf');

  tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(leastCSPMovementAlignedEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(leastCSPMovementAlignedEvents(36:96,:),[],2,'omitnan')./sqrt(length(leastCSPMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(mostCSPMovementAlignedEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(mostCSPMovementAlignedEvents(36:96,:),[],2,'omitnan')./sqrt(length(mostCSPMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L movement - postCue - CS+'); ylim([0 ymax]); 
      saveas(tempFig, [output_fn 'stats\HL_piezoRate_postCue_wholeTC_CS+_PSTHs_' blockName '.pdf'],'pdf');
    
% lowVhigh postRew
for c = 1:sum(nsize(1,:))
low_postRewPiezo_activity_CSP(:,c) = leastCSPMovementAlignedEvents(postRewRange(1)-1+(peakIdx_lowPostRewPiezoActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_lowPostRewPiezoActCSP_allN+statsWindow),c);  
high_postRewPiezo_activity_CSP(:,c) = mostCSPMovementAlignedEvents(postRewRange(1)-1+(peakIdx_highPostRewPiezoActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_highPostRewPiezoActCSP_allN+statsWindow) ,c);   
end
[~, p.piezoRate.CSP_postRew_trialAvg_Cspk, ~, stats.piezoRate.CSP_postRew_trialAvg_Cspk] = ttest(mean(low_postRewPiezo_activity_CSP,1,'omitnan'),mean(high_postRewPiezo_activity_CSP,1,'omitnan'));
  n.piezoRate.CSP_postRew_lowRate_Cspk = length(cell2mat(leastCSPTrialMovement_ind));
  m.piezoRate.CSP_postRew_lowRate_Cspk = mean(mean(leastCSPMovementAlignedEvents(postRewRange(1)-1+(peakIdx_lowPostRewPiezoActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_lowPostRewPiezoActCSP_allN+statsWindow),:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoRate.CSP_postRew_lowRate_Cspk = std(mean(leastCSPMovementAlignedEvents(postRewRange(1)-1+(peakIdx_lowPostRewPiezoActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_lowPostRewPiezoActCSP_allN+statsWindow),:),1,'omitnan'),[],2,'omitnan')./sqrt(size(leastCSPMovementAlignedEvents,2)).*(1000./frameRateHz);
  n.piezoRate.CSP_postRew_highRate_Cspk = length(cell2mat(mostCSPTrialMovement_ind));
  m.piezoRate.CSP_postRew_highRate_Cspk = mean(mean(mostCSPMovementAlignedEvents(postRewRange(1)-1+(peakIdx_highPostRewPiezoActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_highPostRewPiezoActCSP_allN+statsWindow),:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoRate.CSP_postRew_highRate_Cspk = std(mean(mostCSPMovementAlignedEvents(postRewRange(1)-1+(peakIdx_highPostRewPiezoActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_highPostRewPiezoActCSP_allN+statsWindow),:),1,'omitnan'),[],2,'omitnan')./sqrt(size(mostCSPMovementAlignedEvents,2)).*(1000./frameRateHz);

  tempFig=setFigure; hold on; xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  shadedErrorBar_CV(tt(:,postRewRange)/1000,mean(leastCSPMovementAlignedEvents(postRewRange,:),2,'omitnan').*(1000./frameRateHz), std(leastCSPMovementAlignedEvents(postRewRange,:),[],2,'omitnan')./sqrt(length(leastCSPMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,postRewRange)/1000,mean(mostCSPMovementAlignedEvents(postRewRange,:),2,'omitnan').*(1000./frameRateHz), std(mostCSPMovementAlignedEvents(postRewRange,:),[],2,'omitnan')./sqrt(length(mostCSPMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L movement - postRew - CS+'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_piezoRate_postRew_CS+_PSTHs_' blockName '.pdf'],'pdf');

  tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(leastCSPMovementAlignedEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(leastCSPMovementAlignedEvents(36:96,:),[],2,'omitnan')./sqrt(length(leastCSPMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(mostCSPMovementAlignedEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(mostCSPMovementAlignedEvents(36:96,:),[],2,'omitnan')./sqrt(length(mostCSPMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L movement - postRew - CS+'); ylim([0 ymax]); 
      saveas(tempFig, [output_fn 'stats\HL_piezoRate_postRew_wholeTC_CS+_PSTHs_' blockName '.pdf'],'pdf');

% lowVhigh postCue - CS-
for c = 1:sum(nsize(2,:))
low_postCuePiezo_activity_CSM(:,c) = leastCSMMovementAlignedEvents(postCueRange(1)-1+(peakIdx_lowPostCuePiezoActCSM_allN-statsWindow):postCueRange(1)-1+(peakIdx_lowPostCuePiezoActCSM_allN+statsWindow),c);  
high_postCuePiezo_activity_CSM(:,c) = mostCSMMovementAlignedEvents(postCueRange(1)-1+(peakIdx_highPostCuePiezoActCSM_allN-statsWindow):postCueRange(1)-1+(peakIdx_highPostCuePiezoActCSM_allN+statsWindow),c);   
end
[~, p.piezoRate.CSM_postCue_trialAvg_Cspk, ~, stats.piezoRate.CSM_postCue_trialAvg_Cspk] = ttest(mean(low_postCuePiezo_activity_CSM,1,'omitnan'),mean(high_postCuePiezo_activity_CSM,1,'omitnan'));
  n.piezoRate.CSM_postCue_lowRate_Cspk = length(cell2mat(leastCSMTrialMovement_ind));
  m.piezoRate.CSM_postCue_lowRate_Cspk = mean(mean(leastCSMMovementAlignedEvents(postCueRange(1)-1+(peakIdx_lowPostCuePiezoActCSM_allN-statsWindow):postCueRange(1)-1+(peakIdx_lowPostCuePiezoActCSM_allN+statsWindow),:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoRate.CSM_postCue_lowRate_Cspk = std(mean(leastCSMMovementAlignedEvents(postCueRange(1)-1+(peakIdx_lowPostCuePiezoActCSM_allN-statsWindow):postCueRange(1)-1+(peakIdx_lowPostCuePiezoActCSM_allN+statsWindow),:),1,'omitnan'),[],2,'omitnan')./sqrt(size(leastCSMMovementAlignedEvents,2)).*(1000./frameRateHz);
  n.piezoRate.CSM_postCue_highRate_Cspk = length(cell2mat(mostCSMTrialMovement_ind));
  m.piezoRate.CSM_postCue_highRate_Cspk = mean(mean(mostCSMMovementAlignedEvents(postCueRange(1)-1+(peakIdx_highPostCuePiezoActCSM_allN-statsWindow):postCueRange(1)-1+(peakIdx_highPostCuePiezoActCSM_allN+statsWindow),:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoRate.CSM_postCue_highRate_Cspk = std(mean(mostCSMMovementAlignedEvents(postCueRange(1)-1+(peakIdx_highPostCuePiezoActCSM_allN-statsWindow):postCueRange(1)-1+(peakIdx_highPostCuePiezoActCSM_allN+statsWindow),:),1,'omitnan'),[],2,'omitnan')./sqrt(size(mostCSMMovementAlignedEvents,2)).*(1000./frameRateHz);
  
  tempFig=setFigure; hold on; xline(0,'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(leastCSMMovementAlignedEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(leastCSMMovementAlignedEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(leastCSMMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(mostCSMMovementAlignedEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(mostCSMMovementAlignedEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(mostCSMMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L movement - postCue - CS-'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_piezoRate_postCue_CS-_PSTHs_' blockName '.pdf'],'pdf');

  tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(leastCSMMovementAlignedEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(leastCSMMovementAlignedEvents(36:96,:),[],2,'omitnan')./sqrt(length(leastCSMMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(mostCSMMovementAlignedEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(mostCSMMovementAlignedEvents(36:96,:),[],2,'omitnan')./sqrt(length(mostCSMMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L movement - postCue - CS-'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_piezoRate_postCue_wholeTC_CS-_PSTHs_' blockName '.pdf'],'pdf');
      
% lowVhigh postRew
for c = 1:sum(nsize(2,:))
low_postRewPiezo_activity_CSM(:,c) = leastCSMMovementAlignedEvents(postRewRange(1)-1+(peakIdx_lowPostRewPiezoActCSM_allN-statsWindow):postRewRange(1)-1+(peakIdx_lowPostRewPiezoActCSM_allN+statsWindow),c);  
high_postRewPiezo_activity_CSM(:,c) = mostCSMMovementAlignedEvents(postRewRange(1)-1+(peakIdx_highPostRewPiezoActCSM_allN-statsWindow):postRewRange(1)-1+(peakIdx_highPostRewPiezoActCSM_allN+statsWindow) ,c);   
end
[~, p.piezoRate.CSM_postRew_trialAvg_Cspk, ~, stats.piezoRate.CSM_postRew_trialAvg_Cspk] = ttest(mean(low_postRewPiezo_activity_CSM,1,'omitnan'),mean(high_postRewPiezo_activity_CSM,1,'omitnan'));
  n.piezoRate.CSM_postRew_lowRate_Cspk = length(cell2mat(leastCSMTrialMovement_ind));
  m.piezoRate.CSM_postRew_lowRate_Cspk = mean(mean(leastCSMMovementAlignedEvents(postRewRange(1)-1+(peakIdx_lowPostRewPiezoActCSM_allN-statsWindow):postRewRange(1)-1+(peakIdx_lowPostRewPiezoActCSM_allN+statsWindow),:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoRate.CSM_postRew_lowRate_Cspk = std(mean(leastCSMMovementAlignedEvents(postRewRange(1)-1+(peakIdx_lowPostRewPiezoActCSM_allN-statsWindow):postRewRange(1)-1+(peakIdx_lowPostRewPiezoActCSM_allN+statsWindow),:),1,'omitnan'),[],2,'omitnan')./sqrt(size(leastCSMMovementAlignedEvents,2)).*(1000./frameRateHz);
  n.piezoRate.CSM_postRew_highRate_Cspk = length(cell2mat(mostCSMTrialMovement_ind));
  m.piezoRate.CSM_postRew_highRate_Cspk = mean(mean(mostCSMMovementAlignedEvents(postRewRange(1)-1+(peakIdx_highPostRewPiezoActCSM_allN-statsWindow):postRewRange(1)-1+(peakIdx_highPostRewPiezoActCSM_allN+statsWindow),:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoRate.CSM_postRew_highRate_Cspk = std(mean(mostCSMMovementAlignedEvents(postRewRange(1)-1+(peakIdx_highPostRewPiezoActCSM_allN-statsWindow):postRewRange(1)-1+(peakIdx_highPostRewPiezoActCSM_allN+statsWindow),:),1,'omitnan'),[],2,'omitnan')./sqrt(size(mostCSMMovementAlignedEvents,2)).*(1000./frameRateHz);

  tempFig=setFigure; hold on; xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  shadedErrorBar_CV(tt(:,postRewRange)/1000,mean(leastCSMMovementAlignedEvents(postRewRange,:),2,'omitnan').*(1000./frameRateHz), std(leastCSMMovementAlignedEvents(postRewRange,:),[],2,'omitnan')./sqrt(length(leastCSMMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,postRewRange)/1000,mean(mostCSMMovementAlignedEvents(postRewRange,:),2,'omitnan').*(1000./frameRateHz), std(mostCSMMovementAlignedEvents(postRewRange,:),[],2,'omitnan')./sqrt(length(mostCSMMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L movement - postRew - CS-'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_piezoRate_postRew_CS-_PSTHs_' blockName '.pdf'],'pdf');

  tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,2)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,2)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(leastCSMMovementAlignedEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(leastCSMMovementAlignedEvents(36:96,:),[],2,'omitnan')./sqrt(length(leastCSMMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(mostCSMMovementAlignedEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(mostCSMMovementAlignedEvents(36:96,:),[],2,'omitnan')./sqrt(length(mostCSMMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L movement - CS-'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_piezoRate_postRew_wholeTC_CS-_PSTHs_' blockName '.pdf'],'pdf');

%% data org | high vs low lick rate 
postLick_frames=15;

CSPpreRew_lickAlignHz = sum(plusData.preRew_lickAlign(postLick_frames:3*postLick_frames-1,:),1)./length(prewin_frames+1:prewin_frames+rewDelay_frames)/frameRateHz;
CSPpostRew_lickAlignHz = sum(plusData.postRew_lickAlign(postLick_frames+1:3*postLick_frames-1,:),1)./length(prewin_frames+rewDelay_frames+1:prewin_frames+rewDelay_frames+rewDelay_frames)/frameRateHz;

CSMpreRew_lickAlignHz = sum(plusData.preRew_lickAlign(postLick_frames:3*postLick_frames-1,:),1)./length(prewin_frames+1:prewin_frames+rewDelay_frames)/frameRateHz;
CSMpostRew_lickAlignHz = sum(plusData.postRew_lickAlign(postLick_frames+1:3*postLick_frames-1,:),1)./length(prewin_frames+rewDelay_frames+1:prewin_frames+rewDelay_frames+rewDelay_frames)/frameRateHz;

%
low_preRewHz_CSPEvents=[]; high_preRewHz_CSPEvents=[]; low_preRewHz_CSMEvents=[]; high_preRewHz_CSMEvents=[];
low_preRewHz_CSPLicks=[]; high_preRewHz_CSPLicks=[]; low_preRewHz_CSMLicks=[]; high_preRewHz_CSMLicks=[];
low_preRewHz_CSPLickStart = []; low_preRewHz_CSMLickStart = []; high_preRewHz_CSPLickStart = []; high_preRewHz_CSMLickStart = [];

low_postRewHz_CSPEvents=[]; high_postRewHz_CSPEvents=[]; low_postRewHz_CSMEvents=[]; high_postRewHz_CSMEvents=[];
low_postRewHz_CSPLicks=[]; high_postRewHz_CSPLicks=[]; low_postRewHz_CSMLicks=[]; high_postRewHz_CSMLicks=[];
low_postRewHz_CSPLickStart = []; low_postRewHz_CSMLickStart = []; high_postRewHz_CSPLickStart = []; high_postRewHz_CSMLickStart = [];


%CS+
ntrials=0;
for mouse = 1:plusData.exptCount
    nTrial = (1:plusData.animalTrials(1,mouse))+ntrials; 

temp_preRewHz_CSP = plusData.preRewCSP_lickBurstHz(1,(1:plusData.animalTrials(1,mouse))+ntrials);
temp_postRewHz_CSP = plusData.postRewCSP_lickBurstHz(1,(1:plusData.animalTrials(1,mouse))+ntrials);
nanIdx_preRewHz_CSP = ~isnan(temp_preRewHz_CSP);
nanIdx_postRewHz_CSP = ~isnan(temp_postRewHz_CSP);
[sort_preRewHz_CSP, sort_preRewHz_CSP_ind{1,mouse}] = sort(temp_preRewHz_CSP);
[sort_postRewHz_CSP, sort_postRewHz_CSP_ind{1,mouse}] = sort(temp_postRewHz_CSP);

sort_preRewHz_CSP_ind{1,mouse}(isnan(sort_preRewHz_CSP)) = []; sort_preRewHz_CSP(isnan(sort_preRewHz_CSP)) = []; 
sort_postRewHz_CSP_ind{1,mouse}(isnan(sort_postRewHz_CSP)) = []; sort_postRewHz_CSP(isnan(sort_postRewHz_CSP)) = []; 

    if size(sort_preRewHz_CSP,2)./4 == ceil(size(sort_preRewHz_CSP,2)./4)
    low_preRewHz_CSP{1,mouse} = sort_preRewHz_CSP(1,1:ceil(size(sort_preRewHz_CSP,2)./4)+1);
    low_preRewHz_CSP_ind{1,mouse} = sort_preRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSP,2)./4)+1);
    else
    low_preRewHz_CSP{1,mouse} = sort_preRewHz_CSP(1,1:ceil(size(sort_preRewHz_CSP,2)./4));
    low_preRewHz_CSP_ind{1,mouse} = sort_preRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSP,2)./4));
    end
    high_preRewHz_CSP{1,mouse} = sort_preRewHz_CSP(1,size(sort_preRewHz_CSP,2)-floor(size(sort_preRewHz_CSP,2)./4):size(sort_preRewHz_CSP,2));
    high_preRewHz_CSP_ind{1,mouse} = sort_preRewHz_CSP_ind{1,mouse}(1,size(sort_preRewHz_CSP,2)-floor(size(sort_preRewHz_CSP,2)./4):size(sort_preRewHz_CSP,2));
    high_preRewHz_CSPEvents = [high_preRewHz_CSPEvents nanmean(plusData.groupAlign.events{1,mouse}(:,:,sort_preRewHz_CSP_ind{1,mouse}(1,size(sort_preRewHz_CSP_ind{1,mouse},2)-floor(size(sort_preRewHz_CSP_ind{1,mouse},2)./4):size(sort_preRewHz_CSP_ind{1,mouse},2))),3)];
    low_preRewHz_CSPEvents = [low_preRewHz_CSPEvents nanmean(plusData.groupAlign.events{1,mouse}(:,:,sort_preRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSP_ind{1,mouse},2)./4))),3)];
    high_preRewHz_CSPLicks = [high_preRewHz_CSPLicks plusData.groupAlign.licks{1,mouse}(:,sort_preRewHz_CSP_ind{1,mouse}(1,size(sort_preRewHz_CSP_ind{1,mouse},2)-floor(size(sort_preRewHz_CSP_ind{1,mouse},2)./4):size(sort_preRewHz_CSP_ind{1,mouse},2)))];
    low_preRewHz_CSPLicks = [low_preRewHz_CSPLicks plusData.groupAlign.licks{1,mouse}(:,sort_preRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSP_ind{1,mouse},2)./4)))];
    high_preRewHz_CSPLickStart = [high_preRewHz_CSPLickStart plusData.groupAlign.lickStart{1,mouse}(:,sort_preRewHz_CSP_ind{1,mouse}(1,size(sort_preRewHz_CSP_ind{1,mouse},2)-floor(size(sort_preRewHz_CSP_ind{1,mouse},2)./4):size(sort_preRewHz_CSP_ind{1,mouse},2)))];
    low_preRewHz_CSPLickStart = [low_preRewHz_CSPLickStart plusData.groupAlign.lickStart{1,mouse}(:,sort_preRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSP_ind{1,mouse},2)./4)))];
if size(sort_postRewHz_CSP,2)./4 == ceil(size(sort_postRewHz_CSP,2)./4)
low_postRewHz_CSP{1,mouse} = sort_postRewHz_CSP(1,1:ceil(size(sort_postRewHz_CSP,2)./4)+1);
low_postRewHz_CSP_ind{1,mouse} = sort_postRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSP,2)./4)+1);
else
low_postRewHz_CSP{1,mouse} = sort_postRewHz_CSP(1,1:ceil(size(sort_postRewHz_CSP,2)./4));
low_postRewHz_CSP_ind{1,mouse} = sort_postRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSP,2)./4));
end
high_postRewHz_CSP{1,mouse} = sort_postRewHz_CSP(1,size(sort_postRewHz_CSP,2)-floor(size(sort_postRewHz_CSP,2)./4):size(sort_postRewHz_CSP,2));
high_postRewHz_CSP_ind{1,mouse} = sort_postRewHz_CSP_ind{1,mouse}(1,size(sort_postRewHz_CSP,2)-floor(size(sort_postRewHz_CSP,2)./4):size(sort_postRewHz_CSP,2));
high_postRewHz_CSPEvents = [high_postRewHz_CSPEvents nanmean(plusData.groupAlign.events{1,mouse}(:,:,sort_postRewHz_CSP_ind{1,mouse}(1,size(sort_postRewHz_CSP_ind{1,mouse},2)-floor(size(sort_postRewHz_CSP_ind{1,mouse},2)./4):size(sort_postRewHz_CSP_ind{1,mouse},2))),3)];
low_postRewHz_CSPEvents = [low_postRewHz_CSPEvents nanmean(plusData.groupAlign.events{1,mouse}(:,:,sort_postRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSP_ind{1,mouse},2)./4))),3)];
high_postRewHz_CSPLicks = [high_postRewHz_CSPLicks plusData.groupAlign.licks{1,mouse}(:,sort_postRewHz_CSP_ind{1,mouse}(1,size(sort_postRewHz_CSP_ind{1,mouse},2)-floor(size(sort_postRewHz_CSP_ind{1,mouse},2)./4):size(sort_postRewHz_CSP_ind{1,mouse},2)))];
low_postRewHz_CSPLicks = [low_postRewHz_CSPLicks plusData.groupAlign.licks{1,mouse}(:,sort_postRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSP_ind{1,mouse},2)./4)))];
high_postRewHz_CSPLickStart = [high_postRewHz_CSPLickStart plusData.groupAlign.lickStart{1,mouse}(:,sort_postRewHz_CSP_ind{1,mouse}(1,size(sort_postRewHz_CSP_ind{1,mouse},2)-floor(size(sort_postRewHz_CSP_ind{1,mouse},2)./4):size(sort_postRewHz_CSP_ind{1,mouse},2)))];
low_postRewHz_CSPLickStart = [low_postRewHz_CSPLickStart plusData.groupAlign.lickStart{1,mouse}(:,sort_postRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSP_ind{1,mouse},2)./4)))];

    ntrials = max(nTrial);
end

%CS-
ntrials=0;
for mouse = 1:minusData.exptCount
    nTrial = (1:minusData.animalTrials(1,mouse))+ntrials; 

temp_preRewHz_CSM = minusData.preRewCSM_lickBurstHz(1,(1:minusData.animalTrials(1,mouse))+ntrials);
temp_postRewHz_CSM = minusData.postRewCSM_lickBurstHz(1,(1:minusData.animalTrials(1,mouse))+ntrials);
nanIdx_preRewHz_CSM = ~isnan(temp_preRewHz_CSM);
nanIdx_postRewHz_CSM = ~isnan(temp_postRewHz_CSM);
[sort_preRewHz_CSM, sort_preRewHz_CSM_ind{1,mouse}] = sort(temp_preRewHz_CSM);
[sort_postRewHz_CSM, sort_postRewHz_CSM_ind{1,mouse}] = sort(temp_postRewHz_CSM);
       
sort_preRewHz_CSM_ind{1,mouse}(isnan(sort_preRewHz_CSM)) = []; sort_preRewHz_CSM(isnan(sort_preRewHz_CSM)) = [];
sort_postRewHz_CSM_ind{1,mouse}(isnan(sort_postRewHz_CSM)) = []; sort_postRewHz_CSM(isnan(sort_postRewHz_CSM)) = [];

    if size(sort_preRewHz_CSM,2)./4 == ceil(size(sort_preRewHz_CSM,2)./4)
    low_preRewHz_CSM{1,mouse} = sort_preRewHz_CSM(1,1:ceil(size(sort_preRewHz_CSM,2)./4)+1); 
    low_preRewHz_CSM_ind{1,mouse} = sort_preRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSM,2)./4)+1); 
    else
    low_preRewHz_CSM{1,mouse} = sort_preRewHz_CSM(1,1:ceil(size(sort_preRewHz_CSM,2)./4)); 
    low_preRewHz_CSM_ind{1,mouse} = sort_preRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSM,2)./4)); 
    end
    high_preRewHz_CSM{1,mouse} = sort_preRewHz_CSM(1,size(sort_preRewHz_CSM,2)-floor(size(sort_preRewHz_CSM,2)./4):size(sort_preRewHz_CSM,2));
    high_preRewHz_CSM_ind{1,mouse} = sort_preRewHz_CSM_ind{1,mouse}(1,size(sort_preRewHz_CSM,2)-floor(size(sort_preRewHz_CSM,2)./4):size(sort_preRewHz_CSM,2));
    high_preRewHz_CSMEvents = [high_preRewHz_CSMEvents nanmean(minusData.groupAlign.events{1,mouse}(:,:,sort_preRewHz_CSM_ind{1,mouse}(1,size(sort_preRewHz_CSM_ind{1,mouse},2)-floor(size(sort_preRewHz_CSM_ind{1,mouse},2)./4):size(sort_preRewHz_CSM_ind{1,mouse},2))),3)];
    low_preRewHz_CSMEvents = [low_preRewHz_CSMEvents nanmean(minusData.groupAlign.events{1,mouse}(:,:,sort_preRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSM_ind{1,mouse},2)./4))),3)];  
    high_preRewHz_CSMLicks = [high_preRewHz_CSMLicks minusData.groupAlign.licks{1,mouse}(:,sort_preRewHz_CSM_ind{1,mouse}(1,size(sort_preRewHz_CSM_ind{1,mouse},2)-floor(size(sort_preRewHz_CSM_ind{1,mouse},2)./4):size(sort_preRewHz_CSM_ind{1,mouse},2)))];
    low_preRewHz_CSMLicks = [low_preRewHz_CSMLicks minusData.groupAlign.licks{1,mouse}(:,sort_preRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSM_ind{1,mouse},2)./4)))];  
    high_preRewHz_CSMLickStart = [high_preRewHz_CSMLickStart minusData.groupAlign.lickStart{1,mouse}(:,sort_preRewHz_CSM_ind{1,mouse}(1,size(sort_preRewHz_CSM_ind{1,mouse},2)-floor(size(sort_preRewHz_CSM_ind{1,mouse},2)./4):size(sort_preRewHz_CSM_ind{1,mouse},2)))];
    low_preRewHz_CSMLickStart = [low_preRewHz_CSMLickStart minusData.groupAlign.lickStart{1,mouse}(:,sort_preRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSM_ind{1,mouse},2)./4)))];
if size(sort_postRewHz_CSM,2)./4 == ceil(size(sort_postRewHz_CSM,2)./4)
low_postRewHz_CSM{1,mouse} = sort_postRewHz_CSM(1,1:ceil(size(sort_postRewHz_CSM,2)./4)+1); 
low_postRewHz_CSM_ind{1,mouse} = sort_postRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSM,2)./4)+1); 
else
low_postRewHz_CSM{1,mouse} = sort_postRewHz_CSM(1,1:ceil(size(sort_postRewHz_CSM,2)./4)); 
low_postRewHz_CSM_ind{1,mouse} = sort_postRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSM,2)./4)); 
end
high_postRewHz_CSM{1,mouse} = sort_postRewHz_CSM(1,size(sort_postRewHz_CSM,2)-floor(size(sort_postRewHz_CSM,2)./4):size(sort_postRewHz_CSM,2));
high_postRewHz_CSM_ind{1,mouse} = sort_postRewHz_CSM_ind{1,mouse}(1,size(sort_postRewHz_CSM,2)-floor(size(sort_postRewHz_CSM,2)./4):size(sort_postRewHz_CSM,2));
high_postRewHz_CSMEvents = [high_postRewHz_CSMEvents nanmean(minusData.groupAlign.events{1,mouse}(:,:,sort_postRewHz_CSM_ind{1,mouse}(1,size(sort_postRewHz_CSM_ind{1,mouse},2)-floor(size(sort_postRewHz_CSM_ind{1,mouse},2)./4):size(sort_postRewHz_CSM_ind{1,mouse},2))),3)];
low_postRewHz_CSMEvents = [low_postRewHz_CSMEvents nanmean(minusData.groupAlign.events{1,mouse}(:,:,sort_postRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSM_ind{1,mouse},2)./4))),3)];  
high_postRewHz_CSMLicks = [high_postRewHz_CSMLicks minusData.groupAlign.licks{1,mouse}(:,sort_postRewHz_CSM_ind{1,mouse}(1,size(sort_postRewHz_CSM_ind{1,mouse},2)-floor(size(sort_postRewHz_CSM_ind{1,mouse},2)./4):size(sort_postRewHz_CSM_ind{1,mouse},2)))];
low_postRewHz_CSMLicks = [low_postRewHz_CSMLicks minusData.groupAlign.licks{1,mouse}(:,sort_postRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSM_ind{1,mouse},2)./4)))];  
high_postRewHz_CSMLickStart = [high_postRewHz_CSMLickStart minusData.groupAlign.lickStart{1,mouse}(:,sort_postRewHz_CSM_ind{1,mouse}(1,size(sort_postRewHz_CSM_ind{1,mouse},2)-floor(size(sort_postRewHz_CSM_ind{1,mouse},2)./4):size(sort_postRewHz_CSM_ind{1,mouse},2)))];
low_postRewHz_CSMLickStart = [low_postRewHz_CSMLickStart minusData.groupAlign.lickStart{1,mouse}(:,sort_postRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSM_ind{1,mouse},2)./4)))];

    ntrials = max(nTrial);
end
%% ttest    | high vs low lick rate avg vs. baseline avg 
  [peak_lowPostCueLickActCSP, peakIdx_lowPostCueLickActCSP] = max(low_preRewHz_CSPEvents(postCueRange,:),[],1);
  [peak_lowPostRewLickActCSP, peakIdx_lowPostRewLickActCSP] = max(low_postRewHz_CSPEvents(postRewRange,:),[],1);
  [peak_highPostCueLickActCSP, peakIdx_highPostCueLickActCSP] = max(high_preRewHz_CSPEvents(postCueRange,:),[],1);
  [peak_highPostRewLickActCSP, peakIdx_highPostRewLickActCSP] = max(high_postRewHz_CSPEvents(postRewRange,:),[],1);
  [peak_lowPostCueLickActCSM, peakIdx_lowPostCueLickActCSM] = max(low_preRewHz_CSMEvents(postCueRange,:),[],1);
  [peak_lowPostRewLickActCSM, peakIdx_lowPostRewLickActCSM] = max(low_postRewHz_CSMEvents(postRewRange,:),[],1);
  [peak_highPostCueLickActCSM, peakIdx_highPostCueLickActCSM] = max(high_preRewHz_CSMEvents(postCueRange,:),[],1);
  [peak_highPostRewLickActCSM, peakIdx_highPostRewLickActCSM] = max(high_postRewHz_CSMEvents(postRewRange,:),[],1);

  [peak_lowPostCueLickActCSP_allN, peakIdx_lowPostCueLickActCSP_allN] = max(mean(low_preRewHz_CSPEvents(postCueRange,:),2,'omitnan'),[],1);
  [peak_lowPostRewLickActCSP_allN, peakIdx_lowPostRewLickActCSP_allN] = max(mean(low_postRewHz_CSPEvents(postRewRange,:),2,'omitnan'),[],1);
  [peak_highPostCueLickActCSP_allN, peakIdx_highPostCueLickActCSP_allN] = max(mean(high_preRewHz_CSPEvents(postCueRange,:),2,'omitnan'),[],1);
  [peak_highPostRewLickActCSP_allN, peakIdx_highPostRewLickActCSP_allN] = max(mean(high_postRewHz_CSPEvents(postRewRange,:),2,'omitnan'),[],1);
  [peak_lowPostCueLickActCSM_allN, peakIdx_lowPostCueLickActCSM_allN] = max(mean(low_preRewHz_CSMEvents(postCueRange,:),2,'omitnan'),[],1);
  [peak_lowPostRewLickActCSM_allN, peakIdx_lowPostRewLickActCSM_allN] = max(mean(low_postRewHz_CSMEvents(postRewRange,:),2,'omitnan'),[],1);
  [peak_highPostCueLickActCSM_allN, peakIdx_highPostCueLickActCSM_allN] = max(mean(high_preRewHz_CSMEvents(postCueRange,:),2,'omitnan'),[],1);
  [peak_highPostRewLickActCSM_allN, peakIdx_highPostRewLickActCSM_allN] = max(mean(high_postRewHz_CSMEvents(postRewRange,:),2,'omitnan'),[],1);

% lowVhigh postCue - CS+
for c = 1:sum(nsize(1,:))
low_postCueLick_activity_CSP(:,c) = low_preRewHz_CSPEvents(postCueRange(1)-1+(peakIdx_lowPostCueLickActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_lowPostCueLickActCSP_allN+statsWindow),c);  
high_postCueLick_activity_CSP(:,c) = high_preRewHz_CSPEvents(postCueRange(1)-1+(peakIdx_highPostCueLickActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_highPostCueLickActCSP_allN+statsWindow) ,c);   
end
[~, p.lickRate.CSP_postCue_trialAvg_Cspk, ~, stats.lickRate.CSP_postCue_trialAvg_Cspk] = ttest(mean(low_postCueLick_activity_CSP,1,'omitnan'),mean(high_postCueLick_activity_CSP,1,'omitnan'));
  n.lickRate.CSP_postCue_lowRate_Cspk = length(cell2mat(low_preRewHz_CSP_ind));
  m.lickRate.CSP_postCue_lowRate_Cspk = mean(mean(low_preRewHz_CSPEvents(postCueRange(1)-1+(peakIdx_lowPostCueLickActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_lowPostCueLickActCSP_allN+statsWindow),:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickRate.CSP_postCue_lowRate_Cspk = std(mean(low_preRewHz_CSPEvents(postCueRange(1)-1+(peakIdx_lowPostCueLickActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_lowPostCueLickActCSP_allN+statsWindow),:),1,'omitnan'),[],2,'omitnan')./sqrt(size(low_preRewHz_CSPEvents,2)).*(1000./frameRateHz);
  n.lickRate.CSP_postCue_highRate_Cspk = length(cell2mat(high_preRewHz_CSP_ind));
  m.lickRate.CSP_postCue_highRate_Cspk = mean(mean(high_preRewHz_CSPEvents(postCueRange(1)-1+(peakIdx_highPostCueLickActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_highPostCueLickActCSP_allN+statsWindow),:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickRate.CSP_postCue_highRate_Cspk = std(mean(high_preRewHz_CSPEvents(postCueRange(1)-1+(peakIdx_highPostCueLickActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_highPostCueLickActCSP_allN+statsWindow),:),1,'omitnan'),[],2,'omitnan')./sqrt(size(high_preRewHz_CSPEvents,2)).*(1000./frameRateHz);

  tempFig=setFigure; hold on; xline(0,'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(low_preRewHz_CSPEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(low_preRewHz_CSPEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(low_preRewHz_CSPEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(high_preRewHz_CSPEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(high_preRewHz_CSPEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(high_preRewHz_CSPEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L lick rate - postCue - CS+'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_lickRate_postCue_CS+_PSTHs_' blockName '.pdf'],'pdf');

  tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(low_preRewHz_CSPEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(low_preRewHz_CSPEvents(36:96,:),[],2,'omitnan')./sqrt(length(low_preRewHz_CSPEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(high_preRewHz_CSPEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(high_preRewHz_CSPEvents(36:96,:),[],2,'omitnan')./sqrt(length(high_preRewHz_CSPEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L lick rate - postCue - CS+'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_lickRate_postCue_wholeTC_CS+_PSTHs_' blockName '.pdf'],'pdf');

% lowVhigh postRew
for c = 1:sum(nsize(1,:))
low_postRewLick_activity_CSP(:,c) = low_postRewHz_CSPEvents(postRewRange(1)-1+(peakIdx_lowPostRewLickActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_lowPostRewLickActCSP_allN+statsWindow),c);  
high_postRewLick_activity_CSP(:,c) = high_postRewHz_CSPEvents(postRewRange(1)-1+(peakIdx_highPostRewLickActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_highPostRewLickActCSP_allN+statsWindow) ,c);   
end
[~, p.lickRate.CSP_postRew_trialAvg_Cspk, ~, stats.lickRate.CSP_postRew_trialAvg_Cspk] = ttest(mean(low_postRewLick_activity_CSP,1,'omitnan'),mean(high_postRewLick_activity_CSP,1,'omitnan'));
  n.lickRate.CSP_postRew_lowRate_Cspk = length(cell2mat(low_postRewHz_CSP_ind));
  m.lickRate.CSP_postRew_lowRate_Cspk = mean(mean(low_postRewHz_CSPEvents(postRewRange(1)-1+(peakIdx_lowPostRewLickActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_lowPostRewLickActCSP_allN+statsWindow),:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickRate.CSP_postRew_lowRate_Cspk = std(mean(low_postRewHz_CSPEvents(postRewRange(1)-1+(peakIdx_lowPostRewLickActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_lowPostRewLickActCSP_allN+statsWindow),:),1,'omitnan'),[],2,'omitnan')./sqrt(size(low_postRewHz_CSPEvents,2)).*(1000./frameRateHz);
  n.lickRate.CSP_postRew_highRate_Cspk = length(cell2mat(high_postRewHz_CSP_ind));
  m.lickRate.CSP_postRew_highRate_Cspk = mean(mean(high_postRewHz_CSPEvents(postRewRange(1)-1+(peakIdx_highPostRewLickActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_highPostRewLickActCSP_allN+statsWindow),:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickRate.CSP_postRew_highRate_Cspk = std(mean(high_postRewHz_CSPEvents(postRewRange(1)-1+(peakIdx_highPostRewLickActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_highPostRewLickActCSP_allN+statsWindow),:),1,'omitnan'),[],2,'omitnan')./sqrt(size(high_postRewHz_CSPEvents,2)).*(1000./frameRateHz);

  tempFig=setFigure; hold on; xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  shadedErrorBar_CV(tt(:,postRewRange)/1000,mean(low_postRewHz_CSPEvents(postRewRange,:),2,'omitnan').*(1000./frameRateHz), std(low_postRewHz_CSPEvents(postRewRange,:),[],2,'omitnan')./sqrt(length(low_postRewHz_CSPEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,postRewRange)/1000,mean(high_postRewHz_CSPEvents(postRewRange,:),2,'omitnan').*(1000./frameRateHz), std(high_postRewHz_CSPEvents(postRewRange,:),[],2,'omitnan')./sqrt(length(high_postRewHz_CSPEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L lick rate - postRew - CS+'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_lickRate_postRew_CS+_PSTHs_' blockName '.pdf'],'pdf');

  tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(low_postRewHz_CSPEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(low_postRewHz_CSPEvents(36:96,:),[],2,'omitnan')./sqrt(length(low_postRewHz_CSPEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(high_postRewHz_CSPEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(high_postRewHz_CSPEvents(36:96,:),[],2,'omitnan')./sqrt(length(high_postRewHz_CSPEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L lick rate - postRew - CS+'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_lickRate_postRew_wholeTC_CS+_PSTHs_' blockName '.pdf'],'pdf');

% lowVhigh postCue - CS-
for c = 1:sum(nsize(2,:))
low_postCueLick_activity_CSM(:,c) = low_preRewHz_CSMEvents(postCueRange(1)-1+(peakIdx_lowPostCueLickActCSM_allN-statsWindow):postCueRange(1)-1+(peakIdx_lowPostCueLickActCSM_allN+statsWindow),c);  
high_postCueLick_activity_CSM(:,c) = high_preRewHz_CSMEvents(postCueRange(1)-1+(peakIdx_highPostCueLickActCSM_allN-statsWindow):postCueRange(1)-1+(peakIdx_highPostCueLickActCSM_allN+statsWindow) ,c);   
end
[~, p.lickRate.CSM_postCue_trialAvg_Cspk, ~, stats.lickRate.CSM_postCue_trialAvg_Cspk] = ttest(mean(low_postCueLick_activity_CSM,1,'omitnan'),mean(high_postCueLick_activity_CSM,1,'omitnan'));
  n.lickRate.CSM_postCue_lowRate_Cspk = length(cell2mat(low_preRewHz_CSM_ind));
  m.lickRate.CSM_postCue_lowRate_Cspk = mean(mean(low_preRewHz_CSMEvents(postCueRange(1)-1+(peakIdx_lowPostCueLickActCSM_allN-statsWindow):postCueRange(1)-1+(peakIdx_lowPostCueLickActCSM_allN+statsWindow),:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickRate.CSM_postCue_lowRate_Cspk = std(mean(low_preRewHz_CSMEvents(postCueRange(1)-1+(peakIdx_lowPostCueLickActCSM_allN-statsWindow):postCueRange(1)-1+(peakIdx_lowPostCueLickActCSM_allN+statsWindow),:),1,'omitnan'),[],2,'omitnan')./sqrt(size(low_preRewHz_CSMEvents,2)).*(1000./frameRateHz);
  n.lickRate.CSM_postCue_highRate_Cspk = length(cell2mat(high_preRewHz_CSM_ind));
  m.lickRate.CSM_postCue_highRate_Cspk = mean(mean(high_preRewHz_CSMEvents(postCueRange(1)-1+(peakIdx_highPostCueLickActCSM_allN-statsWindow):postCueRange(1)-1+(peakIdx_highPostCueLickActCSM_allN+statsWindow),:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickRate.CSM_postCue_highRate_Cspk = std(mean(high_preRewHz_CSMEvents(postCueRange(1)-1+(peakIdx_highPostCueLickActCSM_allN-statsWindow):postCueRange(1)-1+(peakIdx_highPostCueLickActCSM_allN+statsWindow),:),1,'omitnan'),[],2,'omitnan')./sqrt(size(high_preRewHz_CSMEvents,2)).*(1000./frameRateHz);
  
  tempFig=setFigure; hold on; xline(0,'k');  yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(low_preRewHz_CSMEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(low_preRewHz_CSMEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(low_preRewHz_CSMEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(high_preRewHz_CSMEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(high_preRewHz_CSMEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(high_preRewHz_CSMEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L lick rate - postCue - CS-'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_lickRate_postCue_CS-_PSTHs_' blockName '.pdf'],'pdf');
  
  tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(low_preRewHz_CSMEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(low_preRewHz_CSMEvents(36:96,:),[],2,'omitnan')./sqrt(length(low_preRewHz_CSMEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(high_preRewHz_CSMEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(high_preRewHz_CSMEvents(36:96,:),[],2,'omitnan')./sqrt(length(high_preRewHz_CSMEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L lick rate - postCue - CS-'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_lickRate_postCue_wholeTC_CS-_PSTHs_' blockName '.pdf'],'pdf');
      
% lowVhigh postRew
for c = 1:sum(nsize(2,:))
low_postRewLick_activity_CSM(:,c) = low_postRewHz_CSMEvents(postRewRange(1)-1+(peakIdx_lowPostRewLickActCSM_allN-statsWindow):postRewRange(1)-1+(peakIdx_lowPostRewLickActCSM_allN+statsWindow),c);  
high_postRewLick_activity_CSM(:,c) = high_postRewHz_CSMEvents(postRewRange(1)-1+(peakIdx_highPostRewLickActCSM_allN-statsWindow):postRewRange(1)-1+(peakIdx_highPostRewLickActCSM_allN+statsWindow),c);   
end
[~, p.lickRate.CSM_postRew_trialAvg_Cspk, ~, stats.lickRate.CSM_postRew_trialAvg_Cspk] = ttest(mean(low_postRewLick_activity_CSM,1,'omitnan'),mean(high_postRewLick_activity_CSM,1,'omitnan'));
  n.lickRate.CSM_postRew_lowRate_Cspk = length(cell2mat(low_postRewHz_CSM_ind));
  m.lickRate.CSM_postRew_lowRate_Cspk = mean(mean(low_postRewHz_CSMEvents(postRewRange(1)-1+(peakIdx_lowPostRewLickActCSM_allN-statsWindow):postRewRange(1)-1+(peakIdx_lowPostRewLickActCSM_allN+statsWindow),:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickRate.CSM_postRew_lowRate_Cspk = std(mean(low_postRewHz_CSMEvents(postRewRange(1)-1+(peakIdx_lowPostRewLickActCSM_allN-statsWindow):postRewRange(1)-1+(peakIdx_lowPostRewLickActCSM_allN+statsWindow),:),1,'omitnan'),[],2,'omitnan')./sqrt(size(low_postRewHz_CSMEvents,2)).*(1000./frameRateHz);
  n.lickRate.CSM_postRew_highRate_Cspk = length(cell2mat(high_postRewHz_CSM_ind));
  m.lickRate.CSM_postRew_highRate_Cspk = mean(mean(high_postRewHz_CSMEvents(postRewRange(1)-1+(peakIdx_highPostRewLickActCSM_allN-statsWindow):postRewRange(1)-1+(peakIdx_highPostRewLickActCSM_allN+statsWindow),:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickRate.CSM_postRew_highRate_Cspk = std(mean(high_postRewHz_CSMEvents(postRewRange(1)-1+(peakIdx_highPostRewLickActCSM_allN-statsWindow):postRewRange(1)-1+(peakIdx_highPostRewLickActCSM_allN+statsWindow),:),1,'omitnan'),[],2,'omitnan')./sqrt(size(high_postRewHz_CSMEvents,2)).*(1000./frameRateHz);

  tempFig=setFigure; hold on;  xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  shadedErrorBar_CV(tt(:,postRewRange)/1000,mean(low_postRewHz_CSMEvents(postRewRange,:),2,'omitnan').*(1000./frameRateHz), std(low_postRewHz_CSMEvents(postRewRange,:),[],2,'omitnan')./sqrt(length(low_postRewHz_CSMEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,postRewRange)/1000,mean(high_postRewHz_CSMEvents(postRewRange,:),2,'omitnan').*(1000./frameRateHz), std(high_postRewHz_CSMEvents(postRewRange,:),[],2,'omitnan')./sqrt(length(high_postRewHz_CSMEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L lick rate - postRew - CS-'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_lickRate_postCue_CS-_PSTHs_' blockName '.pdf'],'pdf');
  
  tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(low_preRewHz_CSMEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(low_preRewHz_CSMEvents(36:96,:),[],2,'omitnan')./sqrt(length(low_preRewHz_CSMEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(high_preRewHz_CSMEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(high_preRewHz_CSMEvents(36:96,:),[],2,'omitnan')./sqrt(length(high_preRewHz_CSMEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L lick rate - postRew - CS-'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_lickRate_postRew_wholeTC_CS-_PSTHs_' blockName '.pdf'],'pdf');
      
%% ttest    | earliest vs latest licks avg vs. baseline avg 
earlyRewBurstStart = plusData.earlyRewBurstStart;
lateRewBurstStart = plusData.lateRewBurstStart;
earlyRewBurstEvents = plusData.earlyRewBurstEvents;
lateRewBurstEvents = plusData.lateRewBurstEvents;
  [peak_earlyPostCueLickActCSP, peakIdx_earlyPostCueLickActCSP] = max(earlyRewBurstEvents(postCueRange,:),[],1);
  [peak_earlyPostRewLickActCSP, peakIdx_earlyPostRewLickActCSP] = max(earlyRewBurstEvents(postRewRange,:),[],1);
  [peak_latePostCueLickActCSP, peakIdx_latePostCueLickActCSP] = max(lateRewBurstEvents(postCueRange,:),[],1);
  [peak_latePostRewLickActCSP, peakIdx_latePostRewLickActCSP] = max(lateRewBurstEvents(postRewRange,:),[],1);

  [peak_earlyPostCueLickActCSP_allN, peakIdx_earlyPostCueLickActCSP_allN] = max(mean(earlyRewBurstEvents(postCueRange,:),2,'omitnan'),[],1);
  [peak_earlyPostRewLickActCSP_allN, peakIdx_earlyPostRewLickActCSP_allN] = max(mean(earlyRewBurstEvents(postRewRange,:),2,'omitnan'),[],1);
  [peak_latePostCueLickActCSP_allN, peakIdx_latePostCueLickActCSP_allN] = max(mean(lateRewBurstEvents(postCueRange,:),2,'omitnan'),[],1);
  [peak_latePostRewLickActCSP_allN, peakIdx_latePostRewLickActCSP_allN] = max(mean(lateRewBurstEvents(postRewRange,:),2,'omitnan'),[],1);

n.lickTiming.uEarlyTiming = (mean(earlyRewBurstStart,2,'omitnan')-50).*(1000./frameRateHz);
n.lickTiming.sEarlyTiming = (std(earlyRewBurstStart,[],2,'omitnan')./sqrt(sum(plusData.animalTrials))).*(1000./frameRateHz);
n.lickTiming.uLateTiming = (mean(lateRewBurstStart,2,'omitnan')-50).*(1000./frameRateHz);
n.lickTiming.sLateTiming = (std(lateRewBurstStart,[],2,'omitnan')./sqrt(sum(plusData.animalTrials))).*(1000./frameRateHz);
[~, p.lickTiming.diffInTiming] = ttest(mean(earlyRewBurstStart,1,'omitnan'),mean(lateRewBurstStart,1,'omitnan'));

% earlyVlate postCue - CS+
for c = 1:sum(nsize)
early_postCueLick_activity_CSP(:,c) = earlyRewBurstEvents(postCueRange(1)-1+(peakIdx_earlyPostCueLickActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_earlyPostCueLickActCSP_allN+statsWindow),c);  
late_postCueLick_activity_CSP(:,c) = lateRewBurstEvents(postCueRange(1)-1+(peakIdx_latePostCueLickActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_latePostCueLickActCSP_allN+statsWindow),c);   
end
[~, p.lickTiming.CSP_postCue_trialAvg_Cspk, ~, stats.lickTiming.CSP_postCue_trialAvg_Cspk] = ttest(mean(early_postCueLick_activity_CSP,1,'omitnan'),mean(late_postCueLick_activity_CSP,1,'omitnan'));
  n.lickTiming.CSP_postCue_early_Cspk = length((earlyRewBurstStart));
  m.lickTiming.CSP_postCue_early_Cspk = mean(mean(earlyRewBurstEvents(postCueRange,:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickTiming.CSP_postCue_early_Cspk = std(mean(earlyRewBurstEvents(postCueRange,:),1,'omitnan'),[],2,'omitnan')./sqrt(size(earlyRewBurstEvents,2)).*(1000./frameRateHz);
  n.lickTiming.CSP_postCue_late_Cspk = length((lateRewBurstStart));
  m.lickTiming.CSP_postCue_late_Cspk = mean(mean(lateRewBurstEvents(postCueRange,:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickTiming.CSP_postCue_late_Cspk = std(mean(lateRewBurstEvents(postCueRange,:),1,'omitnan'),[],2,'omitnan')./sqrt(size(lateRewBurstEvents,2)).*(1000./frameRateHz);

  tempFig=setFigure; hold on; xline(0,'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(earlyRewBurstEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(earlyRewBurstEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(earlyRewBurstEvents)).*(1000./frameRateHz),'lineProps','k');
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(lateRewBurstEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(lateRewBurstEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(lateRewBurstEvents)).*(1000./frameRateHz),'lineProps','b');
  title('PSTH of early vs. late licks - postCue - CS+'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\EL_lickTiming_postCue_CS+_PSTHs_' blockName '.pdf'],'pdf');

  tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(earlyRewBurstEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(earlyRewBurstEvents(36:96,:),[],2,'omitnan')./sqrt(length(earlyRewBurstEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(lateRewBurstEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(lateRewBurstEvents(36:96,:),[],2,'omitnan')./sqrt(length(lateRewBurstEvents)).*(1000./frameRateHz),'lineProps','k');
  plot(mean(((earlyRewBurstStart-prewin_frames)*(1000/frameRateHz)),2,'omitnan')/1000, 2, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5); 
  line([mean(((earlyRewBurstStart-prewin_frames)*(1000/frameRateHz)),2,'omitnan')/1000 - (std(((earlyRewBurstStart-prewin_frames)*(1000/frameRateHz)),[],2,'omitnan')/1000)/sqrt(size(expt,2)) ...
      , mean(((earlyRewBurstStart-prewin_frames)*(1000/frameRateHz)),2,'omitnan')/1000 + (std(((earlyRewBurstStart-prewin_frames)*(1000/frameRateHz)),[],2,'omitnan')/1000)/sqrt(size(expt,2))], [2,2], 'Color', 'k', 'LineWidth', 1); % Horizontal line for STD
  plot(mean(((lateRewBurstStart-prewin_frames)*(1000/frameRateHz)),2,'omitnan')/1000, 2, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5); 
  line([mean(((lateRewBurstStart-prewin_frames)*(1000/frameRateHz)),2,'omitnan')/1000 - (std(((lateRewBurstStart-prewin_frames)*(1000/frameRateHz)),[],2,'omitnan')/1000)/sqrt(size(expt,2)) ...
      , mean(((lateRewBurstStart-prewin_frames)*(1000/frameRateHz)),2,'omitnan')/1000 + (std(((lateRewBurstStart-prewin_frames)*(1000/frameRateHz)),[],2,'omitnan')/1000)/sqrt(size(expt,2))], [2,2], 'Color', 'b', 'LineWidth', 1); % Horizontal line for STD
  title('PSTH of early vs. late licks - postCue - CS+'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\EL_lickTiming_wholeTC_CS+_PSTHs_' blockName '.pdf'],'pdf');
   
% earlyVlate postRew
for c = 1:sum(nsize)
early_postRewLick_activity_CSP(:,c) = earlyRewBurstEvents(postRewRange(1)-1+(peakIdx_earlyPostRewLickActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_earlyPostRewLickActCSP_allN+statsWindow),c);  
late_postRewLick_activity_CSP(:,c) = lateRewBurstEvents(postRewRange(1)-1+(peakIdx_latePostRewLickActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_latePostRewLickActCSP_allN+statsWindow),c);   
end
[~, p.lickTiming.CSP_postRew_trialAvg_Cspk, ~, stats.lickTiming.CSP_postRew_trialAvg_Cspk] = ttest(mean(early_postRewLick_activity_CSP,1,'omitnan'),mean(late_postRewLick_activity_CSP,1,'omitnan'));
  n.lickTiming.CSP_postRew_early_Cspk = length((earlyRewBurstStart));
  m.lickTiming.CSP_postRew_early_Cspk = mean(mean(earlyRewBurstEvents(postRewRange,:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickTiming.CSP_postRew_early_Cspk = std(mean(earlyRewBurstEvents(postRewRange,:),1,'omitnan'),[],2,'omitnan')./sqrt(size(earlyRewBurstEvents,2)).*(1000./frameRateHz);
  n.lickTiming.CSP_postRew_late_Cspk = length((lateRewBurstStart));
  m.lickTiming.CSP_postRew_late_Cspk = mean(mean(lateRewBurstEvents(postRewRange,:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickTiming.CSP_postRew_late_Cspk = std(mean(lateRewBurstEvents(postRewRange,:),1,'omitnan'),[],2,'omitnan')./sqrt(size(lateRewBurstEvents,2)).*(1000./frameRateHz);

  tempFig=setFigure; hold on; xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  shadedErrorBar_CV(tt(:,postRewRange)/1000,mean(earlyRewBurstEvents(postRewRange,:),2,'omitnan').*(1000./frameRateHz), std(earlyRewBurstEvents(postRewRange,:),[],2,'omitnan')./sqrt(length(earlyRewBurstEvents)).*(1000./frameRateHz),'lineProps','k');
  shadedErrorBar_CV(tt(:,postRewRange)/1000,mean(lateRewBurstEvents(postRewRange,:),2,'omitnan').*(1000./frameRateHz), std(lateRewBurstEvents(postRewRange,:),[],2,'omitnan')./sqrt(length(lateRewBurstEvents)).*(1000./frameRateHz),'lineProps','b');
  title('PSTH of early vs. late licks - postRew - CS+'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\EL_lickTiming_postRew_CS+_PSTHs_' blockName '.pdf'],'pdf');
        
%% ttest    | lick-align avg vs. cue-align avg 
[~, p.lickCue.CSP_postCue_trialAvg_Cspk, ~, stats.lickCue.CSP_postCue_trialAvg_Cspk] = ttest(mean(CSPevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSP_allN+statsWindow),:),1,'omitnan'),mean(postCue_lickAlignCSPevents_acrossTrials(preLickRange(1)-1+peakIdx_postCueLickActCSP_allN-statsWindow:preLickRange(1)-1+peakIdx_postCueLickActCSP_allN+statsWindow,:),1,'omitnan'),'tail','left');
[~, p.lickCue.CSM_postCue_trialAvg_Cspk, ~, stats.lickCue.CSM_postCue_trialAvg_Cspk] = ttest(mean(CSMevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSM_allN-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSM_allN+statsWindow),:),1,'omitnan'),mean(postCue_lickAlignCSMevents_acrossTrials(preLickRange(1)-1+peakIdx_postCueLickActCSM_allN-statsWindow:preLickRange(1)-1+peakIdx_postCueLickActCSM_allN+statsWindow,:),1,'omitnan'),'tail','left');
[~, p.lickCue.CSP_postRew_trialAvg_Cspk, ~, stats.lickCue.CSP_postRew_trialAvg_Cspk] = ttest(mean(CSPevents_acrossTrials(postRewRange(1)-1+(peakIdx_rewActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_rewActCSP_allN+statsWindow),:),1,'omitnan'),mean(postRew_lickAlignCSPevents_acrossTrials(preLickRange(1)-1+peakIdx_postRewLickActCSP_allN-statsWindow:preLickRange(1)-1+peakIdx_postRewLickActCSP_allN+statsWindow,:),1,'omitnan'),'tail','left');
[~, p.lickCue.CSM_postRew_trialAvg_Cspk, ~, stats.lickCue.CSM_postRew_trialAvg_Cspk] = ttest(mean(CSMevents_acrossTrials(postRewRange(1)-1+(peakIdx_rewActCSM_allN-statsWindow):postRewRange(1)-1+(peakIdx_rewActCSM_allN+statsWindow),:),1,'omitnan'),mean(postRew_lickAlignCSMevents_acrossTrials(preLickRange(1)-1+peakIdx_postRewLickActCSM_allN-statsWindow:preLickRange(1)-1+peakIdx_postRewLickActCSM_allN+statsWindow,:),1,'omitnan'),'tail','left');
%% ttest    | piezo-align avg vs. cue-align avg 
[~, p.piezoCue.CSP_postCue_trialAvg_Cspk, ~, stats.piezoCue.CSP_postCue_trialAvg_Cspk] = ttest(mean(CSPevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSP_allN+statsWindow),:),1,'omitnan'),mean(postCue_piezoAlignCSPevents_acrossTrials(prePiezoRange(1)-1+peakIdx_postCuePiezoActCSP_allN-statsWindow:prePiezoRange(1)-1+peakIdx_postCuePiezoActCSP_allN+statsWindow,:),1,'omitnan'),'tail','left');
[~, p.piezoCue.CSM_postCue_trialAvg_Cspk, ~, stats.piezoCue.CSM_postCue_trialAvg_Cspk] = ttest(mean(CSMevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSM_allN-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSM_allN+statsWindow),:),1,'omitnan'),mean(postCue_piezoAlignCSMevents_acrossTrials(prePiezoRange(1)-1+peakIdx_postCuePiezoActCSM_allN-statsWindow:prePiezoRange(1)-1+peakIdx_postCuePiezoActCSM_allN+statsWindow,:),1,'omitnan'),'tail','left');
[~, p.piezoCue.CSP_postRew_trialAvg_Cspk, ~, stats.piezoCue.CSP_postRew_trialAvg_Cspk] = ttest(mean(CSPevents_acrossTrials(postRewRange(1)-1+(peakIdx_rewActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_rewActCSP_allN+statsWindow),:),1,'omitnan'),mean(postRew_piezoAlignCSPevents_acrossTrials(prePiezoRange(1)-1+peakIdx_postRewPiezoActCSP_allN-statsWindow:prePiezoRange(1)-1+peakIdx_postRewPiezoActCSP_allN+statsWindow,:),1,'omitnan'),'tail','left');
[~, p.piezoCue.CSM_postRew_trialAvg_Cspk, ~, stats.piezoCue.CSM_postRew_trialAvg_Cspk] = ttest(mean(CSMevents_acrossTrials(postRewRange(1)-1+(peakIdx_rewActCSM_allN-statsWindow):postRewRange(1)-1+(peakIdx_rewActCSM_allN+statsWindow),:),1,'omitnan'),mean(postRew_piezoAlignCSMevents_acrossTrials(prePiezoRange(1)-1+peakIdx_postRewPiezoActCSM_allN-statsWindow:prePiezoRange(1)-1+peakIdx_postRewPiezoActCSM_allN+statsWindow,:),1,'omitnan'),'tail','left');
%% ttest    | CS+ avg vs. CS- avg 
[p.cueCue.peak_postCue_trialAvg_Cspk, ~, stats.cueCue.peak_postCue_trialAvg_Cspk] = ranksum(mean(CSMevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSP_allN+statsWindow),:),1,'omitnan'),mean(CSPevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSP_allN-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSP_allN+statsWindow),:),1,'omitnan'));
[p.cueCue.peak_postRew_trialAvg_Cspk, ~, stats.cueCue.peak_postRew_trialAvg_Cspk] = ranksum(mean(CSMevents_acrossTrials(postRewRange(1)-1+(peakIdx_rewActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_rewActCSP_allN+statsWindow),:),1,'omitnan'),mean(CSPevents_acrossTrials(postRewRange(1)-1+(peakIdx_rewActCSP_allN-statsWindow):postRewRange(1)-1+(peakIdx_rewActCSP_allN+statsWindow),:),1,'omitnan'));

  tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('Average Cspk for CS+/- trials with statistical windows');
rectangle('Position',[tt(1,postCueRange(1)),-1+.01,(1000./frameRateHz)*response_frames,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postRewRange(1)),-1+.01,(1000./frameRateHz)*response_frames,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:96),mean(CSPevents_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSPevents_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSPevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:96),mean(CSMevents_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSMevents_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSMevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
%% save stats 
 save(fullfile([output_fn 'stats\simpleStatistics.mat']), 'p','m','s','stats','n','postCue*','postReward*','peak*','cell*','pos','neg','frameIdx','*Range');
%% Make tables for all t-tests

cueAlign_table = struct2table(p.cueAlign);
lickAlign_table = struct2table(p.lickAlign);
piezoAlign_table = struct2table(p.piezoAlign);
save(fullfile([output_fn 'stats\tableOfP-values.mat']), '*_table');
clearvars *_table

cueAlign_table = struct2table(m.cueAlign);
lickAlign_table = struct2table(m.lickAlign);
piezoAlign_table = struct2table(m.piezoAlign);
save(fullfile([output_fn 'stats\tableOfMeans.mat']), '*_table');
clearvars *_table

cueAlign_table = struct2table(s.cueAlign);
lickAlign_table = struct2table(s.lickAlign);
piezoAlign_table = struct2table(s.piezoAlign);
save(fullfile([output_fn 'stats\tableOfStdDev.mat']), '*_table');
clearvars *_table

cueAlign_table = struct2table(stats.cueAlign);
lickAlign_table = struct2table(stats.lickAlign);
piezoAlign_table = struct2table(stats.piezoAlign);
save(fullfile([output_fn 'stats\tableOfTStats.mat']), '*_table');
 %}