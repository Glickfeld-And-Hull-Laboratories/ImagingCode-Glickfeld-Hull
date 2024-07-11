%% group-wise analysis that outputs interleaved figures for presentation

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
load_fn = ([analysis_out 'groupwiseAlign' num2str(thresholdDeco) '_' sesh '_' sesh2 '\']);
output_fn = ([analysis_out 'groupwiseAlign_' sesh '_' sesh2 '_threshold' num2str(thresholdDeco) '_240317\']);
if seshType ~= 2
load(fullfile(load_fn, '_groupwiseDataset.mat'));
elseif seshType == 2
disp('There is a separate script to pull statistics from blocked sessions.');
return
end

analysis_out = 'A:\home\carlo\RC\analysis\2P\';
load_fn = ([analysis_out 'groupwiseAlign' num2str(thresholdDeco) '_' sesh '_' sesh2 '\']);
output_fn = ([analysis_out 'groupwiseAlign_' sesh '_' sesh2 '_threshold' num2str(thresholdDeco) '_240317\']);
if ~exist(fullfile([output_fn 'stats\'])) 
 mkdir(fullfile([output_fn 'stats\']))
end
%% data org | fixed variables 
frameRateHz = 30;
prewin_frames = round(1500./frameRateHz);
postwin_frames = round(3000./frameRateHz);
response_delay = 2; %provides a window for the average delay between cue and cue-response 
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

iDend=0; iTrial=0;
for mouse = 1:exptCount
 dends=(1:nsize(1,mouse))+iDend; trials=(1:animalTrials(1,mouse))+iTrial;
  cueAlignEvents(:,dends,trials)=groupAlign.events{1,mouse}(:,1:nsize(1,mouse),1:animalTrials(1,mouse));
 iDend=dends(end); iTrial=trials(end);
end

ymax=input(['Set ymax: ']); 
ymin=input(['Set ymin: ']);
%% data org | cue-aligned 
for c = 1:size(expt,2)
    nsize(1,c) = size(groupAlign.tc{1,c},2);
end

CSMtimecourse_acrossTrials=[]; CSPtimecourse_acrossTrials=[]; CSMbaseline_acrossTrials=[]; CSPbaseline_acrossTrials=[]; CSMevents_acrossTrials=[]; CSPevents_acrossTrials=[]; 
CSMtimecourse_acrossNeurons=[]; CSPtimecourse_acrossNeurons=[]; CSMbaseline_acrossNeurons=[]; CSPbaseline_acrossNeurons=[]; CSMevents_acrossNeurons=[]; CSPevents_acrossNeurons=[];
first_lickAlignCSP_acrossTrials=[]; firstPostCue_lickAlignCSP_acrossTrials=[]; firstPostRew_lickAlignCSP_acrossTrials=[]; first_lickAlignCSM_acrossTrials=[]; firstPostCue_lickAlignCSM_acrossTrials=[]; firstPostRew_lickAlignCSM_acrossTrials=[];
first_piezoAlignCSP_acrossTrials=[]; firstPostCue_piezoAlignCSP_acrossTrials=[]; firstPostRew_piezoAlignCSP_acrossTrials=[]; first_piezoAlignCSM_acrossTrials=[]; firstPostCue_piezoAlignCSM_acrossTrials=[]; firstPostRew_piezoAlignCSM_acrossTrials=[];
CSMdFoF_acrossTrials=[]; CSPdFoF_acrossTrials=[]; ntrials=0;        
for mouse = 1:exptCount
nTrial = (1:animalTrials(1,mouse))+ntrials;  

for c = 1:nsize(1,mouse)
CSMdFoF_acrossTrials = [CSMdFoF_acrossTrials mean(groupAlign.dFoF{1,mouse}(:,c,logical(block2{1,mouse})),3,'omitnan')];%-mean(groupAlign.tc{1,mouseIdx}(1:prewin_frames,c,logical(block2{1,mouseIdx})),1,'omitnan'),3,'omitnan')];
CSPdFoF_acrossTrials = [CSPdFoF_acrossTrials mean(groupAlign.dFoF{1,mouse}(:,c,logical(rew{1,mouse})),3,'omitnan')];%-mean(groupAlign.tc{1,mouseIdx}(1:prewin_frames,c,logical(rew{1,mouseIdx})),1,'omitnan'),3,'omitnan')];
CSMevents_acrossTrials = [CSMevents_acrossTrials nanmean(groupAlign.events{1,mouse}(:,c,logical(block2{1,mouse})),3)];
CSPevents_acrossTrials = [CSPevents_acrossTrials nanmean(groupAlign.events{1,mouse}(:,c,logical(rew{1,mouse})),3)];
end

CSMtimecourse_acrossTrials = [CSMtimecourse_acrossTrials nanmean(groupAlign.tc{1,mouse}(:,:,logical(block2{1,mouse})),3)];
CSPtimecourse_acrossTrials = [CSPtimecourse_acrossTrials nanmean(groupAlign.tc{1,mouse}(:,:,logical(rew{1,mouse})),3)];
CSMbaseline_acrossTrials = [CSMbaseline_acrossTrials nanmean(groupAlign.tc{1,mouse}(1:prewin_frames,:,logical(block2{1,mouse})),3)]; % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.
CSPbaseline_acrossTrials = [CSPbaseline_acrossTrials nanmean(groupAlign.tc{1,mouse}(1:prewin_frames,:,logical(rew{1,mouse})),3)]; % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.

CSMtimecourse_acrossNeurons = [CSMtimecourse_acrossNeurons reshape(nanmean(groupAlign.tc{1,mouse}(:,:,logical(block2{1,mouse})),2),[150 size(groupAlign.tc{1,mouse}(:,:,logical(block2{1,mouse})),3)])];
CSPtimecourse_acrossNeurons = [CSPtimecourse_acrossNeurons reshape(nanmean(groupAlign.tc{1,mouse}(:,:,logical(rew{1,mouse})),2),[150 size(groupAlign.tc{1,mouse}(:,:,logical(rew{1,mouse})),3)])];
CSMbaseline_acrossNeurons = [CSMbaseline_acrossNeurons reshape(nanmean(groupAlign.tc{1,mouse}(1:prewin_frames,:,logical(block2{1,mouse})),2),[50 size(groupAlign.tc{1,mouse}(1:prewin_frames,:,logical(block2{1,mouse})),3)])]; % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.
CSPbaseline_acrossNeurons = [CSPbaseline_acrossNeurons reshape(nanmean(groupAlign.tc{1,mouse}(1:prewin_frames,:,logical(rew{1,mouse})),2),[50 size(groupAlign.tc{1,mouse}(1:prewin_frames,:,logical(rew{1,mouse})),3)])]; % average across trials, use nanmean because if the imaging data is cutted due to z shift, then cTargetOn will indexing nans.
CSMevents_acrossNeurons = [CSMevents_acrossNeurons reshape(nanmean(groupAlign.events{1,mouse}(:,:,logical(block2{1,mouse})),2),[150 size(groupAlign.events{1,mouse}(:,:,logical(block2{1,mouse})),3)])];
CSPevents_acrossNeurons = [CSPevents_acrossNeurons reshape(nanmean(groupAlign.events{1,mouse}(:,:,logical(rew{1,mouse})),2),[150 size(groupAlign.events{1,mouse}(:,:,logical(rew{1,mouse})),3)])];

ntrials=max(nTrial);
end

first_lickAlignCSP_acrossTrials = [first_lickAlignCSP_acrossTrials sum(~isnan(firstLickStart(:,~gBlock2)))];
first_lickAlignCSM_acrossTrials = [first_lickAlignCSM_acrossTrials sum(~isnan(firstLickStart(:,gBlock2)))];
firstPostCue_lickAlignCSP_acrossTrials = [firstPostCue_lickAlignCSP_acrossTrials sum(~isnan(firstPreRew_lickStart(:,~gBlock2)))];
firstPostCue_lickAlignCSM_acrossTrials = [firstPostCue_lickAlignCSM_acrossTrials sum(~isnan(firstPreRew_lickStart(:,gBlock2)))];
firstPostRew_lickAlignCSP_acrossTrials = [firstPostRew_lickAlignCSP_acrossTrials sum(~isnan(firstPostRew_lickStart(:,~gBlock2)))];
firstPostRew_lickAlignCSM_acrossTrials = [firstPostRew_lickAlignCSM_acrossTrials sum(~isnan(firstPostRew_lickStart(:,gBlock2)))];

first_piezoAlignCSP_acrossTrials = [first_piezoAlignCSP_acrossTrials sum(~isnan(firstPiezoStart(:,~gBlock2)))];
first_piezoAlignCSM_acrossTrials = [first_piezoAlignCSM_acrossTrials sum(~isnan(firstPiezoStart(:,gBlock2)))];
firstPostCue_piezoAlignCSP_acrossTrials = [firstPostCue_piezoAlignCSP_acrossTrials sum(~isnan(firstPreRew_piezoStart(:,~gBlock2)))];
firstPostCue_piezoAlignCSM_acrossTrials = [firstPostCue_piezoAlignCSM_acrossTrials sum(~isnan(firstPreRew_piezoStart(:,gBlock2)))];
firstPostRew_piezoAlignCSP_acrossTrials = [firstPostRew_piezoAlignCSP_acrossTrials sum(~isnan(firstPostRew_piezoStart(:,~gBlock2)))];
firstPostRew_piezoAlignCSM_acrossTrials = [firstPostRew_piezoAlignCSM_acrossTrials sum(~isnan(firstPostRew_piezoStart(:,gBlock2)))];

first_lickAlignEventsCSP_acrossTrials=[]; first_lickAlignEventsCSM_acrossTrials=[]; postCue_lickAlignEventsCSP_acrossTrials=[]; postCue_lickAlignEventsCSM_acrossTrials=[]; postRew_lickAlignEventsCSP_acrossTrials=[]; postRew_lickAlignEventsCSM_acrossTrials=[];
first_piezoAlignEventsCSP_acrossTrials=[]; first_piezoAlignEventsCSM_acrossTrials=[]; postCue_piezoAlignEventsCSP_acrossTrials=[]; postCue_piezoAlignEventsCSM_acrossTrials=[]; postRew_piezoAlignEventsCSP_acrossTrials=[]; postRew_piezoAlignEventsCSM_acrossTrials=[];
first_lickAlignDFoFCSP_acrossTrials=[]; first_lickAlignDFoFCSM_acrossTrials=[]; postCue_lickAlignDFoFCSP_acrossTrials=[]; postCue_lickAlignDFoFCSM_acrossTrials=[]; postRew_lickAlignDFoFCSP_acrossTrials=[]; postRew_lickAlignDFoFCSM_acrossTrials=[];
first_piezoAlignDFoFCSP_acrossTrials=[]; first_piezoAlignDFoFCSM_acrossTrials=[]; postCue_piezoAlignDFoFCSP_acrossTrials=[]; postCue_piezoAlignDFoFCSM_acrossTrials=[]; postRew_piezoAlignDFoFCSP_acrossTrials=[]; postRew_piezoAlignDFoFCSM_acrossTrials=[];

first_lickAlignEventsCSP_acrossTrials = [first_lickAlignEventsCSP_acrossTrials mean(firstLickAlignEvents(:,:,~gBlock2),3,'omitnan')]; first_lickAlignEventsCSM_acrossTrials = [first_lickAlignEventsCSM_acrossTrials mean(firstLickAlignEvents(:,:,gBlock2),3,'omitnan')];
postCue_lickAlignEventsCSP_acrossTrials = [postCue_lickAlignEventsCSP_acrossTrials mean(firstPreRew_lickAlignEvents(:,:,~gBlock2),3,'omitnan')]; postCue_lickAlignEventsCSM_acrossTrials = [postCue_lickAlignEventsCSM_acrossTrials mean(firstPreRew_lickAlignEvents(:,:,gBlock2),3,'omitnan')];
postRew_lickAlignEventsCSP_acrossTrials = [postRew_lickAlignEventsCSP_acrossTrials mean(firstPostRew_lickAlignEvents(:,:,~gBlock2),3,'omitnan')]; postRew_lickAlignEventsCSM_acrossTrials = [postRew_lickAlignEventsCSM_acrossTrials mean(firstPostRew_lickAlignEvents(:,:,gBlock2),3,'omitnan')];

first_piezoAlignEventsCSP_acrossTrials = [first_piezoAlignEventsCSP_acrossTrials mean(firstPiezoAlignEvents(:,:,~gBlock2),3,'omitnan')]; first_piezoAlignEventsCSM_acrossTrials = [first_piezoAlignEventsCSM_acrossTrials mean(firstPiezoAlignEvents(:,:,gBlock2),3,'omitnan')];
postCue_piezoAlignEventsCSP_acrossTrials = [postCue_piezoAlignEventsCSP_acrossTrials mean(firstPreRew_piezoAlignEvents(:,:,~gBlock2),3,'omitnan')]; postCue_piezoAlignEventsCSM_acrossTrials = [postCue_piezoAlignEventsCSM_acrossTrials mean(firstPreRew_piezoAlignEvents(:,:,gBlock2),3,'omitnan')];
postRew_piezoAlignEventsCSP_acrossTrials = [postRew_piezoAlignEventsCSP_acrossTrials mean(firstPostRew_piezoAlignEvents(:,:,~gBlock2),3,'omitnan')]; postRew_piezoAlignEventsCSM_acrossTrials = [postRew_piezoAlignEventsCSM_acrossTrials mean(firstPostRew_piezoAlignEvents(:,:,gBlock2),3,'omitnan')];

first_lickAlignDFoFCSP_acrossTrials = [first_lickAlignDFoFCSP_acrossTrials mean(firstLickAlignDFoF(:,:,~gBlock2),3,'omitnan')]; first_lickAlignDFoFCSM_acrossTrials = [first_lickAlignDFoFCSM_acrossTrials mean(firstLickAlignDFoF(:,:,gBlock2),3,'omitnan')];
postCue_lickAlignDFoFCSP_acrossTrials = [postCue_lickAlignDFoFCSP_acrossTrials mean(firstPreRew_lickAlignDFoF(:,:,~gBlock2),3,'omitnan')]; postCue_lickAlignDFoFCSM_acrossTrials = [postCue_lickAlignDFoFCSM_acrossTrials mean(firstPreRew_lickAlignDFoF(:,:,gBlock2),3,'omitnan')];
postRew_lickAlignDFoFCSP_acrossTrials = [postRew_lickAlignDFoFCSP_acrossTrials mean(firstPostRew_lickAlignDFoF(:,:,~gBlock2),3,'omitnan')]; postRew_lickAlignDFoFCSM_acrossTrials = [postRew_lickAlignDFoFCSM_acrossTrials mean(firstPostRew_lickAlignDFoF(:,:,gBlock2),3,'omitnan')];

first_piezoAlignDFoFCSP_acrossTrials = [first_piezoAlignDFoFCSP_acrossTrials mean(firstPiezoAlignDFoF(:,:,~gBlock2),3,'omitnan')]; first_piezoAlignDFoFCSM_acrossTrials = [first_piezoAlignDFoFCSM_acrossTrials mean(firstPiezoAlignDFoF(:,:,gBlock2),3,'omitnan')];
postCue_piezoAlignDFoFCSP_acrossTrials = [postCue_piezoAlignDFoFCSP_acrossTrials mean(firstPreRew_piezoAlignDFoF(:,:,~gBlock2),3,'omitnan')]; postCue_piezoAlignDFoFCSM_acrossTrials = [postCue_piezoAlignDFoFCSM_acrossTrials mean(firstPreRew_piezoAlignDFoF(:,:,gBlock2),3,'omitnan')];
postRew_piezoAlignDFoFCSP_acrossTrials = [postRew_piezoAlignDFoFCSP_acrossTrials mean(firstPostRew_piezoAlignDFoF(:,:,~gBlock2),3,'omitnan')]; postRew_piezoAlignDFoFCSM_acrossTrials = [postRew_piezoAlignDFoFCSM_acrossTrials mean(firstPostRew_piezoAlignDFoF(:,:,gBlock2),3,'omitnan')];

colorWheel ={[200/255 78/255 0/255],[153/255 51/255 153/255],[161/255 183/255 13/255],[255/255 217/255 96/255]};
%
%% data org | segmentation of neurons into response groups by means AND change (+/-) 
eventLabel={'postCue','postRew'}; cueLabel={'CSP','CSM'}; 
for cue=1:size(cueLabel,2)
  for event=1:size(eventLabel,2)
      eval(['dendrite.cueAlign.response_' eventLabel{1,event} cueLabel{1,cue} '=[];']);
      eval(['dendrite.cueAlign.suppress_' eventLabel{1,event} cueLabel{1,cue} '=[];']);
    for c = 1:sum(nsize)
%
 eval(['windowVar = ' cueLabel{1,cue} 'dFoF_acrossTrials(' eventLabel{1,event} 'Range,c);']);
 eval(['baselineVar = ' cueLabel{1,cue} 'dFoF_acrossTrials(baselineRange,c);']);
 eval(['pos.Diff.cueAlign{c,1} = ([0; diff(' cueLabel{1,cue} 'dFoF_acrossTrials(:,c))]);']);

pos.Limit.cueAlign{c,1}.one = (mean(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),1,'omitnan')+(std(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),[],1,'omitnan').*1)); 
pos.Change.cueAlign{c,1}.one = find(pos.Diff.cueAlign{c,1}>pos.Limit.cueAlign{c,1}.one); 
pos.Limit.cueAlign{c,1}.two = (mean(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),1,'omitnan')+(std(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),[],1,'omitnan').*2)); 
pos.Change.cueAlign{c,1}.two = find(pos.Diff.cueAlign{c,1}>pos.Limit.cueAlign{c,1}.two); 
pos.Limit.cueAlign{c,1}.three = (mean(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),1,'omitnan')+(std(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),[],1,'omitnan').*3)); 
pos.Change.cueAlign{c,1}.three = find(pos.Diff.cueAlign{c,1}>pos.Limit.cueAlign{c,1}.three);

eval(['frameIdx.cueAlign.' eventLabel{1,event} cueLabel{1,cue} '_posChangeDiff{1,c} = intersect(' eventLabel{1,event} 'Range, pos.Change.cueAlign{c,1}.one);']); 
eval(['frameIdx.cueAlign.' eventLabel{1,event} cueLabel{1,cue} '_posChangeDiff{2,c} = intersect(' eventLabel{1,event} 'Range, pos.Change.cueAlign{c,1}.two);']); 
eval(['frameIdx.cueAlign.' eventLabel{1,event} cueLabel{1,cue} '_posChangeDiff{3,c} = intersect(' eventLabel{1,event} 'Range, pos.Change.cueAlign{c,1}.three);']);

 eval(['neg.Diff.cueAlign{c,1} = ([0; diff(' cueLabel{1,cue} 'dFoF_acrossTrials(:,c))]);']);

neg.Limit.cueAlign{c,1}.one = (mean(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),1,'omitnan')-(std(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),[],1,'omitnan').*1)); 
neg.Change.cueAlign{c,1}.one = find(neg.Diff.cueAlign{c,1}<neg.Limit.cueAlign{c,1}.one); 
neg.Limit.cueAlign{c,1}.two = (mean(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),1,'omitnan')-(std(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),[],1,'omitnan').*2)); 
neg.Change.cueAlign{c,1}.two = find(neg.Diff.cueAlign{c,1}<neg.Limit.cueAlign{c,1}.two); 
neg.Limit.cueAlign{c,1}.three = (mean(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),1,'omitnan')-(std(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),[],1,'omitnan').*3)); 
neg.Change.cueAlign{c,1}.three = find(neg.Diff.cueAlign{c,1}<neg.Limit.cueAlign{c,1}.three);

eval(['frameIdx.cueAlign.' eventLabel{1,event} cueLabel{1,cue} '_negChangeDiff{1,c} = intersect(' eventLabel{1,event} 'Range, neg.Change.cueAlign{c,1}.one);']); 
eval(['frameIdx.cueAlign.' eventLabel{1,event} cueLabel{1,cue} '_negChangeDiff{2,c} = intersect(' eventLabel{1,event} 'Range, neg.Change.cueAlign{c,1}.two);']); 
eval(['frameIdx.cueAlign.' eventLabel{1,event} cueLabel{1,cue} '_negChangeDiff{3,c} = intersect(' eventLabel{1,event} 'Range, neg.Change.cueAlign{c,1}.three);']);
%    

    tempResp_DiffIdx = diff(intersect(1:response_frames, find(windowVar>(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*respReq_cueAlign))))';
    tempResp_FrameIdx = intersect(1:response_frames, find(windowVar>(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*respReq_cueAlign)));
    if length(tempResp_DiffIdx)
      if strfind(tempResp_DiffIdx, [1 1])
        if eval(['frameIdx.cueAlign.' eventLabel{1,event} cueLabel{1,cue} '_posChangeDiff{chngReq_cueAlign,c}'])
           eval(['dendrite.cueAlign.response_' eventLabel{1,event} cueLabel{1,cue} '=[dendrite.cueAlign.response_' eventLabel{1,event} cueLabel{1,cue} ' c];']);
        end
      end
    end
    tempSupp_DiffIdx = diff(intersect(1:response_frames, find(windowVar<(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*suppReq_cueAlign))))';
    tempSupp_FrameIdx = intersect(1:response_frames, find(windowVar<(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*suppReq_cueAlign)));
    if length(tempSupp_DiffIdx)
      if strfind(tempSupp_DiffIdx, [1 1 1 1]) 
        if isempty(strfind(tempResp_DiffIdx, [1 1]))
          if eval(['frameIdx.cueAlign.' eventLabel{1,event} cueLabel{1,cue} '_negChangeDiff{chngReq_cueAlign,c}'])
           eval(['dendrite.cueAlign.suppress_' eventLabel{1,event} cueLabel{1,cue} '=[dendrite.cueAlign.suppress_' eventLabel{1,event} cueLabel{1,cue} ' c];']);
          end
        end
      end
    end
    end
  end
end
%% data org | isolation of cells that fall into a single response group (e.g. removal of neurons which respond and suppress following cue) 
dendriteDiff.cueAlign.response_postCueCSP = setdiff(dendrite.cueAlign.response_postCueCSP,dendrite.cueAlign.suppress_postCueCSP);
dendriteDiff.cueAlign.suppress_postCueCSP = setdiff(dendrite.cueAlign.suppress_postCueCSP,dendrite.cueAlign.response_postCueCSP);
dendriteDiff.cueAlign.response_postRewCSP = setdiff(dendrite.cueAlign.response_postRewCSP,dendrite.cueAlign.suppress_postRewCSP);
dendriteDiff.cueAlign.suppress_postRewCSP = setdiff(dendrite.cueAlign.suppress_postRewCSP,dendrite.cueAlign.response_postRewCSP);
dendriteDiff.cueAlign.response_postCueCSM = setdiff(dendrite.cueAlign.response_postCueCSM,dendrite.cueAlign.suppress_postCueCSM);
dendriteDiff.cueAlign.suppress_postCueCSM = setdiff(dendrite.cueAlign.suppress_postCueCSM,dendrite.cueAlign.response_postCueCSM);
dendriteDiff.cueAlign.response_postRewCSM = setdiff(dendrite.cueAlign.response_postRewCSM,dendrite.cueAlign.suppress_postRewCSM);
dendriteDiff.cueAlign.suppress_postRewCSM = setdiff(dendrite.cueAlign.suppress_postRewCSM,dendrite.cueAlign.response_postRewCSM);

dendriteOverlap.cueAlign.response_postCue = dendrite.cueAlign.response_postCueCSP(1,ismember(dendrite.cueAlign.response_postCueCSP, dendrite.cueAlign.response_postCueCSM));
dendriteOverlap.cueAlign.suppress_postCue = dendrite.cueAlign.suppress_postCueCSP(1,ismember(dendrite.cueAlign.suppress_postCueCSP, dendrite.cueAlign.suppress_postCueCSM));
dendriteOverlap.cueAlign.response_postRew = dendrite.cueAlign.response_postRewCSP(1,ismember(dendrite.cueAlign.response_postRewCSP, dendrite.cueAlign.response_postRewCSM));
dendriteOverlap.cueAlign.suppress_postRew = dendrite.cueAlign.suppress_postRewCSP(1,ismember(dendrite.cueAlign.suppress_postRewCSP, dendrite.cueAlign.suppress_postRewCSM));

dendriteUnion.cueAlign.postCueCSP = union(dendrite.cueAlign.response_postCueCSP,dendrite.cueAlign.suppress_postCueCSP);
dendriteUnion.cueAlign.postRewCSP = union(dendrite.cueAlign.response_postRewCSP,dendrite.cueAlign.suppress_postRewCSP);
dendriteUnion.cueAlign.postCueCSM = union(dendrite.cueAlign.response_postCueCSM,dendrite.cueAlign.suppress_postCueCSM);
dendriteUnion.cueAlign.postRewCSM = union(dendrite.cueAlign.response_postRewCSM,dendrite.cueAlign.suppress_postRewCSM);

dendriteUnique.cueAlign.postCueCSP = setdiff(dendriteDiff.cueAlign.response_postCueCSP,dendriteDiff.cueAlign.response_postCueCSM);
dendriteUnique.cueAlign.postRewCSP = setdiff(dendriteDiff.cueAlign.response_postRewCSP,dendriteDiff.cueAlign.response_postRewCSM);
dendriteUnique.cueAlign.postCueCSM = setdiff(dendriteDiff.cueAlign.response_postCueCSM,dendriteDiff.cueAlign.response_postCueCSP);
dendriteUnique.cueAlign.postRewCSM = setdiff(dendriteDiff.cueAlign.response_postRewCSM,dendriteDiff.cueAlign.response_postRewCSP);

dendriteUnresp.cueAlign.postCue = setdiff(1:sum(nsize),union(dendriteDiff.cueAlign.response_postCueCSP,dendriteDiff.cueAlign.response_postCueCSM));
dendriteUnresp.cueAlign.postRew = setdiff(1:sum(nsize),union(dendriteDiff.cueAlign.response_postRewCSP,dendriteDiff.cueAlign.response_postRewCSM));

iterArray = {'dendriteDiff.cueAlign.response_postCueCSP','dendriteDiff.cueAlign.suppress_postCueCSP','dendriteDiff.cueAlign.response_postRewCSP','dendriteDiff.cueAlign.suppress_postRewCSP','dendriteDiff.cueAlign.response_postCueCSM','dendriteDiff.cueAlign.suppress_postCueCSM','dendriteDiff.cueAlign.response_postRewCSM','dendriteDiff.cueAlign.suppress_postRewCSM'};
%% plot     | cue-aligned representation of response groups 
allcells = 1:sum(nsize);
%}
respArray = {'response_postCueCSP','response_postCueCSM','suppress_postCueCSP','suppress_postCueCSM','response_postRewCSP','response_postRewCSM','suppress_postRewCSP','suppress_postRewCSM'};
titleArray = {'response postCue CS+','response postCue CS-','suppress postCue CS+','suppress postCue CS-','response postReward CS+','response postReward CS-','suppress postReward CS+','suppress postReward CS-'};
colorArray = {'k','r','b','m','k','r','b','m'};

highlight = [255/255 217/255 96/255]; windowPanes = [baselineRange(1), postCueRange(1), postRewRange(1)];
dataArray = {'CSPdFoF_acrossTrials','CSMdFoF_acrossTrials','CSPdFoF_acrossTrials','CSMdFoF_acrossTrials','CSPdFoF_acrossTrials','CSMdFoF_acrossTrials','CSPdFoF_acrossTrials','CSMdFoF_acrossTrials'};
tempFig = setFigure; hold on;
for c = 1:size(dataArray,2)  
    eval(['callTrue = ' dataArray{1,c} '(36:96,dendriteDiff.cueAlign.' respArray{1,c} ');']);
%     eval(['callFalse = ' dataArray{1,i} '(36:96,setdiff(allcells,dendriteDiff.cueAlign.' respArray{1,i} '));']);
    eval(['dendrites = dendriteDiff.cueAlign.' respArray{1,c} ';']);
subplot(4,2,c); hold on;
 ylim([ymin*2 ymax*2]); xline(0,'k'); xline(.767,'k'); title(titleArray{1,c});
if c<5; rectangle('Position',[tt(1,postCueRange(1))/1000,ymin*2+.15,(1000./frameRateHz)*response_frames/1000,ymax+5],'EdgeColor',highlight,'FaceColor',highlight); elseif c>4; rectangle('Position',[tt(1,postRewRange(1))/1000,ymin*2+.15,(1000./frameRateHz)*response_frames/1000,ymax+5],'EdgeColor',highlight,'FaceColor',highlight); end
shadedErrorBar_CV(tt(1,36:96)/1000,mean(callTrue,2,'omitnan').*(1000./frameRateHz),std(callTrue,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',colorArray{1,c});
%shadedErrorBar_CV(tt(1,36:96)/1000,mean(callFalse,2,'omitnan').*(1000./frameRateHz),std(callFalse,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',{'Color',[153/255 153/255 0/255]});
end
    saveas(tempFig, [output_fn 'stats\dFoF_PSTHs.pdf'],'pdf');

dataArray = {'CSPevents_acrossTrials','CSMevents_acrossTrials','CSPevents_acrossTrials','CSMevents_acrossTrials','CSPevents_acrossTrials','CSMevents_acrossTrials','CSPevents_acrossTrials','CSMevents_acrossTrials'};
tempFig = setFigure; hold on;
for c = 1:size(dataArray,2)  
    eval(['callTrue = ' dataArray{1,c} '(36:96,dendriteDiff.cueAlign.' respArray{1,c} ');']);
%     eval(['callFalse = ' dataArray{1,i} '(36:96,setdiff(allcells,dendriteDiff.cueAlign.' respArray{1,i} '));']);
    eval(['dendrites = dendriteDiff.cueAlign.' respArray{1,c} ';']);
subplot(4,2,c); hold on;
 ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title(titleArray{1,c});
if c<5; rectangle('Position',[tt(1,postCueRange(1))/1000,0+.1,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight); elseif c>4; rectangle('Position',[tt(1,postRewRange(1))/1000,0+.1,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight); end
shadedErrorBar_CV(tt(1,36:96)/1000,mean(callTrue,2,'omitnan').*(1000./frameRateHz),std(callTrue,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',colorArray{1,c});
%shadedErrorBar_CV(tt(1,36:96)/1000,mean(callFalse,2,'omitnan').*(1000./frameRateHz),std(callFalse,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',{'Color',[153/255 153/255 0/255]});
end
    saveas(tempFig, [output_fn 'stats\CspkEvents_PSTHs.pdf'],'pdf');

tempFig=setFigure; hold on; 
 ylim([-1 ymax]); xline(0,'k'); xline(.767,'k'); title('Average dFoF for CS+/- trials with statistical windows');
rectangle('Position',[tt(1,baselineRange(1))/1000,-1+.01,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postCueRange(1))/1000,-1+.01,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postRewRange(1))/1000,-1+.01,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSPdFoF_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSPdFoF_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSPdFoF_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSMdFoF_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSMdFoF_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSMdFoF_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
saveas(tempFig, [output_fn 'stats\dFoF.pdf'],'pdf');

tempFig=setFigure; hold on; 
 ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('Average Cspk for CS+/- trials with statistical windows');
rectangle('Position',[tt(1,baselineRange(1))/1000,-1+.01,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postCueRange(1))/1000,-1+.01,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postRewRange(1))/1000,-1+.01,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSPevents_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSPevents_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSPevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSMevents_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSMevents_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSMevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
saveas(tempFig, [output_fn 'stats\Cspk.pdf'],'pdf');    
    
tempFig=setFigure; hold on; 
 subplot(2,2,1); hold on;
 ylim([0 ymax*2]); xline(0,'k'); xline(.767,'k'); title('Dendrites that respond to CS+/- presentation');
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSPevents_acrossTrials(36:96,dendriteOverlap.cueAlign.response_postCue),2,'omitnan').*(1000./frameRateHz),std(CSPevents_acrossTrials(36:96,dendriteOverlap.cueAlign.response_postCue),[],2,'omitnan')./sqrt(size(CSPevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSMevents_acrossTrials(36:96,dendriteOverlap.cueAlign.response_postCue),2,'omitnan').*(1000./frameRateHz),std(CSMevents_acrossTrials(36:96,dendriteOverlap.cueAlign.response_postCue),[],2,'omitnan')./sqrt(size(CSMevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
 subplot(2,2,2); hold on;
 ylim([0 ymax*2]); xline(0,'k'); xline(.767,'k'); title('Dendrites that respond to CS+ presentation');
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSPevents_acrossTrials(36:96,dendriteUnique.cueAlign.postCueCSP),2,'omitnan').*(1000./frameRateHz),std(CSPevents_acrossTrials(36:96,dendriteUnique.cueAlign.postCueCSP),[],2,'omitnan')./sqrt(size(CSPevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSMevents_acrossTrials(36:96,dendriteUnique.cueAlign.postCueCSP),2,'omitnan').*(1000./frameRateHz),std(CSMevents_acrossTrials(36:96,dendriteUnique.cueAlign.postCueCSP),[],2,'omitnan')./sqrt(size(CSMevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
 subplot(2,2,3); hold on;
 ylim([0 ymax*2]); xline(0,'k'); xline(.767,'k'); title('Dendrites that respond to CS- presentation');
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSPevents_acrossTrials(36:96,dendriteUnique.cueAlign.postCueCSM),2,'omitnan').*(1000./frameRateHz),std(CSPevents_acrossTrials(36:96,dendriteUnique.cueAlign.postCueCSM),[],2,'omitnan')./sqrt(size(CSPevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSMevents_acrossTrials(36:96,dendriteUnique.cueAlign.postCueCSM),2,'omitnan').*(1000./frameRateHz),std(CSMevents_acrossTrials(36:96,dendriteUnique.cueAlign.postCueCSM),[],2,'omitnan')./sqrt(size(CSMevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
 subplot(2,2,4); hold on;
 ylim([0 ymax*2]); xline(0,'k'); xline(.767,'k'); title('Dendrites that do not respond to cue presentation');
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSPevents_acrossTrials(36:96,dendriteUnresp.cueAlign.postCue),2,'omitnan').*(1000./frameRateHz),std(CSPevents_acrossTrials(36:96,dendriteUnresp.cueAlign.postCue),[],2,'omitnan')./sqrt(size(CSPevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSMevents_acrossTrials(36:96,dendriteUnresp.cueAlign.postCue),2,'omitnan').*(1000./frameRateHz),std(CSMevents_acrossTrials(36:96,dendriteUnresp.cueAlign.postCue),[],2,'omitnan')./sqrt(size(CSMevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
saveas(tempFig, [output_fn 'stats\Cspk_cueSpecificity.pdf'],'pdf'); 

    tempFig = setFigure; 
    bothN = size(dendriteOverlap.cueAlign.response_postCue,2);
    CSP_onlyN = size(dendriteUnique.cueAlign.postCueCSP,2);
    CSM_onlyN = size(dendriteUnique.cueAlign.postCueCSM,2);
    noneN = size(dendriteUnresp.cueAlign.postCue,2);
pie([bothN, CSP_onlyN, CSM_onlyN, noneN]); legend('cells responsive to both cues','cells only responsive to CS+','cells only responsive to CS-','cells unresponsive to either cue')
    saveas(tempFig, [output_fn 'stats\pie_propCells_postCueAlign.pdf'],'pdf');


    %
tempFig=setFigure; hold on;
 bv=bar(1,[bothN/sum(nsize), CSP_onlyN/sum(nsize), CSM_onlyN/sum(nsize), noneN/sum(nsize)],'stacked');
 legend('cells responsive to both cues','cells only responsive to CS+','cells only responsive to CS-','cells unresponsive to either cue')
 colors = {[153/255 51/255 153/255],[161/255 183/255 13/255],[200/255 78/255 0/255],[255/255 217/255 96/255]}; 
for i = 1:4
     bv(i).FaceColor=colors{i};
end
 saveas(tempFig, [output_fn 'stats\bar_propCells_postCueAlign.pdf'],'pdf');

 
    tempFig = setFigure; 
    bothN = size(dendriteOverlap.cueAlign.response_postRew,2);
    CSP_onlyN = size(dendrite.cueAlign.response_postRewCSP,2)-size(dendriteOverlap.cueAlign.response_postRew,2);
    CSM_onlyN = size(dendrite.cueAlign.response_postRewCSM,2)-size(dendriteOverlap.cueAlign.response_postRew,2);
    noneN = size(dendriteUnresp.cueAlign.postRew,2);
pie([bothN, CSP_onlyN, CSM_onlyN, noneN]); legend('cells responsive to reward delivery following both cues','cells only responsive to CS+','cells only responsive to CS-','cells unresponsive to reward delivery')
    saveas(tempFig, [output_fn 'stats\pie_propCells_postRewAlign.pdf'],'pdf');
%% ttest    | response-window avg vs. baseline avg  
  n.cueAlign.trialAvg_CSPtrials = sum(~gBlock2);
  n.cueAlign.trialAvg_CSMtrials = sum(gBlock2);
  n.cueAlign.trialAvg_neurons = size(CSPevents_acrossTrials,2);
  n.cueAlign.trialAvg_mice = size(expt,2);

  [peak_cueActCSP, peakIdx_cueActCSP] = max(mean(CSPevents_acrossTrials(postCueRange,:),2,'omitnan'),[],1);
  [peak_rewActCSP, peakIdx_rewActCSP] = max(mean(CSPevents_acrossTrials(postRewRange,:),2,'omitnan'),[],1);
  [peak_baseCSP, peakIdx_baseCSP] = max(mean(CSPevents_acrossTrials(baselineRange,:),2,'omitnan'),[],1);
  [peak_cueActCSM, peakIdx_cueActCSM] = max(mean(CSMevents_acrossTrials(postCueRange,:),2,'omitnan'),[],1);
  [peak_rewActCSM, peakIdx_rewActCSM] = max(mean(CSMevents_acrossTrials(postRewRange,:),2,'omitnan'),[],1);
  [peak_baseCSM, peakIdx_baseCSM] = max(mean(CSMevents_acrossTrials(baselineRange,:),2,'omitnan'),[],1);

  m.cueAlign.CSP_baseline_trialAvg_Cspk = mean(mean(CSPevents_acrossTrials(baselineRange,:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.cueAlign.CSP_baseline_trialAvg_Cspk = std(mean(CSPevents_acrossTrials(baselineRange,:),1,'omitnan'),[],2,'omitnan')./sqrt(size(CSPevents_acrossTrials,2)).*(1000./frameRateHz);
  m.cueAlign.CSM_baseline_trialAvg_Cspk = mean(mean(CSMevents_acrossTrials(baselineRange,:),1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.cueAlign.CSM_baseline_trialAvg_Cspk = std(mean(CSMevents_acrossTrials(baselineRange,:),1,'omitnan'),[],2,'omitnan')./sqrt(size(CSMevents_acrossTrials,2)).*(1000./frameRateHz);

% difference in activity after cue (a.k.a. before reward - compared to pre-cue) - averaged trials :: SAA
for c = 1:iDend
postCue_baseline_CSP(:,c) = CSPevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSP-statsWindow):baselineRange(1)-1+(peakIdx_baseCSP+statsWindow),c);
postCue_activity_CSP(:,c) = CSPevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSP-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSP+statsWindow),c);
end
[~, p.cueAlign.CSP_postCue_trialAvg_Cspk, ~, stats.cueAlign.CSP_postCue_trialAvg_Cspk] = ttest(mean(postCue_baseline_CSP,1,'omitnan'),mean(postCue_activity_CSP,1,'omitnan'));
  m.cueAlign.CSP_postCue_trialAvg_Cspk = mean(mean(postCue_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.cueAlign.CSP_postCue_trialAvg_Cspk = std(mean(postCue_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(postCue_activity_CSP,2)).*(1000./frameRateHz);
for c = 1:iDend
postCue_baseline_CSM(:,c) = CSMevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSM-statsWindow):baselineRange(1)-1+(peakIdx_baseCSM+statsWindow),c);
postCue_activity_CSM(:,c) = CSMevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSM-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSM+statsWindow),c); 
end
[~, p.cueAlign.CSM_postCue_trialAvg_Cspk, ~, stats.cueAlign.CSM_postCue_trialAvg_Cspk] = ttest(mean(postCue_baseline_CSM,1,'omitnan'),mean(postCue_activity_CSM,1,'omitnan'));
  m.cueAlign.CSM_postCue_trialAvg_Cspk = mean(mean(postCue_activity_CSM,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.cueAlign.CSM_postCue_trialAvg_Cspk = std(mean(postCue_activity_CSM,1,'omitnan'),[],2,'omitnan')./sqrt(size(postCue_activity_CSM,2)).*(1000./frameRateHz);

% difference in responsivity to CS+ and CS-
t=1;
for mus = 1:size(gTargetOn,2)
    trialCount_cue{mus}(1) = sum(~gBlock2(1,t:t+sum(~isnan(gTargetOn{mus}),2)-1));
    trialCount_cue{mus}(2) = sum(gBlock2(1,t:t+sum(~isnan(gTargetOn{mus}),2)-1));
    t=sum(~isnan(gTargetOn{mus}),2)+t-1;
end

% difference in activity after reward (compared to post-cue) - averaged trials :: comparing activity across neurons
for c = 1:iDend
postReward_baseline_CSP(:,c) = CSPevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSP-statsWindow):baselineRange(1)-1+(peakIdx_baseCSP+statsWindow),c);
postReward_activity_CSP(:,c) = CSPevents_acrossTrials(postRewRange(1)-1+(peakIdx_rewActCSP-statsWindow):postRewRange(1)-1+(peakIdx_rewActCSP+statsWindow),c);
end
[~, p.cueAlign.CSP_postRew_trialAvg_Cspk, ~, stats.cueAlign.CSP_postRew_trialAvg_Cspk] = ttest(mean(postReward_baseline_CSP,1,'omitnan'),mean(postReward_activity_CSP,1,'omitnan'));
  m.cueAlign.CSP_postRew_trialAvg_Cspk = mean(mean(postReward_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.cueAlign.CSP_postRew_trialAvg_Cspk = std(mean(postReward_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(postReward_activity_CSP,2)).*(1000./frameRateHz);
for c = 1:iDend
postReward_baseline_CSM(:,c) = CSMevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSM-statsWindow):baselineRange(1)-1+(peakIdx_baseCSM+statsWindow),c);  
postReward_activity_CSM(:,c) = CSMevents_acrossTrials(postRewRange(1)-1+(peakIdx_rewActCSM-statsWindow):postRewRange(1)-1+(peakIdx_rewActCSM+statsWindow),c);   
end
[~, p.cueAlign.CSM_postRew_trialAvg_Cspk, ~, stats.cueAlign.CSM_postRew_trialAvg_Cspk] = ttest(mean(postReward_baseline_CSM,1,'omitnan'),mean(postReward_activity_CSM,1,'omitnan'));
  m.cueAlign.CSM_postRew_trialAvg_Cspk = mean(mean(postReward_activity_CSM,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.cueAlign.CSM_postRew_trialAvg_Cspk = std(mean(postReward_activity_CSM,1,'omitnan'),[],2,'omitnan')./sqrt(size(postReward_activity_CSM,2)).*(1000./frameRateHz);

    %% analysis for modulation (CS+ - CS- / CS+ + CS-)
 if seshType == 4
  allCell_peakCSP = mean(mean(CSPevents_acrossTrials(postRewRange(1)-1+(peakIdx_cueActCSP-statsWindow):postRewRange(1)-1+(peakIdx_cueActCSP+statsWindow),:,:),3,'omitnan'),1,'omitnan').*(1000./frameRateHz);
  allCell_peakCSM = mean(mean(CSMevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSM-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSM+statsWindow),:,:),3,'omitnan'),1,'omitnan').*(1000./frameRateHz);
 else
  allCell_peakCSP = mean(mean(CSPevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSP-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSP+statsWindow),:,:),3,'omitnan'),1,'omitnan').*(1000./frameRateHz);
  allCell_peakCSM = mean(mean(CSMevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSM-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSM+statsWindow),:,:),3,'omitnan'),1,'omitnan').*(1000./frameRateHz);
 end
  % Adjust alpha for multiple comparisons using Bonferroni correction
  adjusted_alpha = 0.05 / size(animalTrials,2);
  % Revised loop to use adjusted alpha
  for i = 1:iDend
    [~, p.modulationIndex(i), ~] = ttest(mean(CSPevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSP-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSP+statsWindow),i,:),3,'omitnan'), mean(CSMevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSM-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSM+statsWindow),i,:),3,'omitnan'));
  end
  modulationIndex = (allCell_peakCSP-allCell_peakCSM)./(allCell_peakCSP+allCell_peakCSM);
  fav_CSP_indices = find(allCell_peakCSP > allCell_peakCSM & modulationIndex == 1);
  fav_CSM_indices = find(allCell_peakCSP < allCell_peakCSM & modulationIndex == -1);
  no_preference_indices = find(p.modulationIndex == 0 & modulationIndex ~= 1 | modulationIndex ~= -1);
  
  % Normalize modulation index for color mapping: range from 0 (red) to 1 (black)
  modIndexNormalized = (modulationIndex - min(modulationIndex)) / (max(modulationIndex) - min(modulationIndex));
   
    tempFig=setFigure;
    hold on; 
    % Loop through each cell to plot
    for i = 1:iDend
        % Determine color: Interpolate between red (1, 0, 0) and black (0, 0, 0)
        color = [1, 0, 0] * (1 - modIndexNormalized(i)); % Reduce red component
           % Plotting CSP vs. CSM with modulation-based color
        if ~isnan(modulationIndex(i))
            scatter(allCell_peakCSP(i), allCell_peakCSM(i), 40, 'filled', 'MarkerFaceColor', color);
        end
    end
    % Reference lines
    xLimits = [0 max(allCell_peakCSP)];
    yLimits = [0 max(allCell_peakCSM)];
    plotLimits = [0 max(xLimits(2), yLimits(2))]; % Ensuring lines extend to the edge of data

    % Calculate and plot 33-degree and 66-degree lines
    % 33 degrees
    theta33 = deg2rad(33);
    slope33 = tan(theta33);
    intercept33 = median(allCell_peakCSM) - slope33 * median(allCell_peakCSP);
    y33 = slope33 * plotLimits;
    
    % 66 degrees
    theta66 = deg2rad(66);
    slope66 = tan(theta66);
    intercept66 = median(allCell_peakCSM) - slope66 * median(allCell_peakCSP);
    y66 = slope66 * plotLimits;
    
    % Add the lines to the plot
    plot(plotLimits, y33, 'k--');
    plot(plotLimits, y66, 'r--');
    
    % Categorizing cells based on the zones they fall into
    thresholdCSP = 1;%2*m.cueAlign.CSP_baseline_trialAvg_Cspk;% + 10*s.cueAlign.CSP_baseline_trialAvg_Cspk;
    thresholdCSM = 1;%2*m.cueAlign.CSM_baseline_trialAvg_Cspk;% + 10*s.cueAlign.CSM_baseline_trialAvg_Cspk; 
    %rectangle('Position', [0, 0, thresholdCSP, thresholdCSM], 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
    zone = zeros(iDend, 1);  % 1 for below 33 degrees, 2 for between 33 and 66, 3 for above 66 degrees
    for i = 1:iDend
        if allCell_peakCSP(i) < thresholdCSP && allCell_peakCSM(i) < thresholdCSM
        zone(i) = 4; % Not significantly above baseline
       
        else
        yVal33 = slope33 * allCell_peakCSP(i);
        yVal66 = slope66 * allCell_peakCSP(i);
        
        if allCell_peakCSM(i) < yVal33
            zone(i) = 1;  % Below 33 degrees
        elseif allCell_peakCSM(i) < yVal66
            zone(i) = 2;  % Between 33 and 66 degrees
        else
            zone(i) = 3;  % Above 66 degrees
        end
        end
    end
    
    % Count cells in each zone
    cellPref_CSP = find(zone == 1);
    cellPref_none = find(zone == 2);
    cellPref_CSM = find(zone == 3);
    cellPref_insig = find(zone == 4);
   
    % Enhance plot
    xlabel('Peak CSP Response');
    ylabel('Peak CSM Response');
    title('CSP vs. CSM Responses with Angle Zones');
    grid on; ylim([0 ymax*3]); xlim([0 ymax*3])
    exportgraphics(tempFig, [output_fn 'stats/cuePreferenceApprox_peakSpkEvents_scatter.pdf'], 'ContentType', 'vector');

    tempFig = setFigure; hold on;
    titleArray={'CS+ preference','no preference','CS- preference','unresponsive'};
for c = 1:4     
subplot(4,1,c); hold on;
 ylim([0 ymax+2]); xline(0,'k'); xline(.767,'k'); title(titleArray{1,c});
%if c<5; rectangle('Position',[tt(1,postCueRange(1)-1)/1000,0+.1,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight); elseif c>4; rectangle('Position',[tt(1,postRewRange(1))/1000,0+.1,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight); end
shadedErrorBar_CV(tt(1,36:96)/1000,mean(mean(CSPevents_acrossTrials(36:96,zone==c,:),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(CSPevents_acrossTrials(36:96,zone==c,:),3,'omitnan'),[],2,'omitnan')./sqrt(iDend).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:96)/1000,mean(mean(CSMevents_acrossTrials(36:96,zone==c,:),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(CSMevents_acrossTrials(36:96,zone==c,:),3,'omitnan'),[],2,'omitnan')./sqrt(iDend).*(1000./frameRateHz),'lineProps','r');
end
    saveas(tempFig, [output_fn 'stats\cuePreferenceApprox_cSpkEvents_PSTHs.pdf'],'pdf');

%% plot     | statistical window for reward-window avg vs. baseline avg
peak_baselinePane(1,1) = round(mean(peakIdx_baseCSP,2,'omitnan')); peak_baselinePane(1,2) = round(mean(peakIdx_baseCSM,2,'omitnan'));
peak_postCuePane(1,1) = round(mean(peakIdx_cueActCSP,2,'omitnan')); peak_postCuePane(1,2) = round(mean(peakIdx_cueActCSM,2,'omitnan'));
peak_postRewPane(1,1) = round(mean(peakIdx_rewActCSP,2,'omitnan')); peak_postRewPane(1,2) = round(mean(peakIdx_rewActCSM,2,'omitnan'));

  tempFig=setFigure; hold on; 
subplot(2,1,1); hold on;
 ylim([ymin*2 ymax*2]); xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
 title('Average dFoF for CS+ trials with statistical windows');
rectangle('Position',[tt(1,baselineRange(1)-1+peak_baselinePane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSPdFoF_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSPdFoF_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSPdFoF_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
subplot(2,1,2); hold on;
 ylim([ymin*2 ymax*2]); xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
 title('Average dFoF for CS- trials with statistical windows');
rectangle('Position',[tt(1,baselineRange(1)-1+peak_baselinePane(1,2)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,2)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,2)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSMdFoF_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSMdFoF_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSMdFoF_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
saveas(tempFig, [output_fn 'stats\dFoFWindow.pdf'],'pdf');

tempFig=setFigure; hold on; 
subplot(2,1,1); hold on;
 ylim([0 ymax]); xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
 title('Average Cpsk for CS+ trials with statistical windows');
rectangle('Position',[tt(1,baselineRange(1)-1+peak_baselinePane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSPevents_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSPevents_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSPevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
subplot(2,1,2); hold on;
 ylim([0 ymax]); xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
title('Average cSpk for CS- trials with statistical windows');
rectangle('Position',[tt(1,baselineRange(1)-1+peak_baselinePane(1,2)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,2)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,2)-statsWindow)/1000,-1+.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSMevents_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSMevents_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSMevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
saveas(tempFig, [output_fn 'stats\CspkWindow.pdf'],'pdf');
%% ttest    | lick-present/lick-absent avg vs. baseline avg 

CSPevents_withLick_acrossTrials = []; CSPevents_woutLick_acrossTrials = []; 
CSMevents_withLick_acrossTrials = []; CSMevents_woutLick_acrossTrials = []; 
ntrials=0;  

for mouse = 1:exptCount
    nTrial = (1:animalTrials(1,mouse))+ntrials;  
    for c = 1:nsize(1,mouse)
CSPevents_withLick_acrossTrials = [CSPevents_withLick_acrossTrials mean(groupAlign.events{1,mouse}(:,c,(logical(rew{1,mouse}) & ~isnan(firstLickStart(1,nTrial)))),3,'omitnan')];
CSPevents_woutLick_acrossTrials = [CSPevents_woutLick_acrossTrials mean(groupAlign.events{1,mouse}(:,c,(logical(rew{1,mouse}) & isnan(firstLickStart(1,nTrial)))),3,'omitnan')];
CSMevents_withLick_acrossTrials = [CSMevents_withLick_acrossTrials mean(groupAlign.events{1,mouse}(:,c,(logical(block2{1,mouse}) & ~isnan(firstLickStart(1,nTrial)))),3,'omitnan')];
CSMevents_woutLick_acrossTrials = [CSMevents_woutLick_acrossTrials mean(groupAlign.events{1,mouse}(:,c,(logical(block2{1,mouse}) & isnan(firstLickStart(1,nTrial)))),3,'omitnan')];
    end
ntrials=max(nTrial);
end

  [peak_postCue_withLick_CSP, peakIdx_postCue_withLick_CSP] = max(mean(CSPevents_withLick_acrossTrials(postCueRange,:),2,'omitnan'),[],1);
  [peak_postCue_woutLick_CSP, peakIdx_postCue_woutLick_CSP] = max(mean(CSPevents_woutLick_acrossTrials(postCueRange,:),2,'omitnan'),[],1);
  [peak_postRew_withLick_CSP, peakIdx_postRew_withLick_CSP] = max(mean(CSPevents_withLick_acrossTrials(postRewRange,:),2,'omitnan'),[],1);
  [peak_postRew_woutLick_CSP, peakIdx_postRew_woutLick_CSP] = max(mean(CSPevents_woutLick_acrossTrials(postRewRange,:),2,'omitnan'),[],1);
  [peak_baseCSP, peakIdx_baseCSP] = max(mean(CSPevents_acrossTrials(baselineRange,:),2,'omitnan'),[],1);
  [peak_postCue_withLick_CSM, peakIdx_postCue_withLick_CSM] = max(mean(CSMevents_withLick_acrossTrials(postCueRange,:),2,'omitnan'),[],1);
  [peak_postCue_woutLick_CSM, peakIdx_postCue_woutLick_CSM] = max(mean(CSMevents_woutLick_acrossTrials(postCueRange,:),2,'omitnan'),[],1);
  [peak_postRew_withLick_CSM, peakIdx_postRew_withLick_CSM] = max(mean(CSMevents_withLick_acrossTrials(postRewRange,:),2,'omitnan'),[],1);
  [peak_postRew_woutLick_CSM, peakIdx_postRew_woutLick_CSM] = max(mean(CSMevents_woutLick_acrossTrials(postRewRange,:),2,'omitnan'),[],1);
  [peak_baseCSM, peakIdx_baseCSM] = max(mean(CSMevents_acrossTrials(baselineRange,:),2,'omitnan'),[],1);
%post-cue
for c = 1:iDend
postCue_withBaseline_CSP(:,c) = CSPevents_withLick_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSP-statsWindow):baselineRange(1)-1+(peakIdx_baseCSP+statsWindow),c);
postCue_woutBaseline_CSP(:,c) = CSPevents_woutLick_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSP-statsWindow):baselineRange(1)-1+(peakIdx_baseCSP+statsWindow),c);
postCue_withLick_CSP(:,c) = CSPevents_withLick_acrossTrials(postCueRange(1)-1+(peakIdx_postCue_withLick_CSP-statsWindow):postCueRange(1)-1+(peakIdx_postCue_withLick_CSP+statsWindow),c);
postCue_woutLick_CSP(:,c) = CSPevents_woutLick_acrossTrials(postCueRange(1)-1+(peakIdx_postCue_woutLick_CSP-statsWindow):postCueRange(1)-1+(peakIdx_postCue_woutLick_CSP+statsWindow),c);
end
[~, p.cueAlign.withLick.CSP_postCue_trialAvg_Cspk, ~, stats.cueAlign.withLick.CSP_postCue_trialAvg_Cspk] = ttest(mean(postCue_withBaseline_CSP,1,'omitnan'),mean(postCue_withLick_CSP,1,'omitnan'));
[~, p.cueAlign.woutLick.CSP_postCue_trialAvg_Cspk, ~, stats.cueAlign.woutLick.CSP_postCue_trialAvg_Cspk] = ttest(mean(postCue_woutBaseline_CSP,1,'omitnan'),mean(postCue_woutLick_CSP,1,'omitnan'));
[p.cueAlign.Lick.CSP_postCue_trialAvg_Cspk, ~, stats.cueAlign.Lick.CSP_postCue_trialAvg_Cspk] = ttest(mean(postCue_woutLick_CSP,1,'omitnan'),mean(postCue_withLick_CSP,1,'omitnan'));
  n.cueAlign.withLick.CSP_postCue_trialAvg_Cspk = sum(~gBlock2 & ~isnan(firstLickStart));
  m.cueAlign.withLick.CSP_postCue_trialAvg_Cspk = mean(mean(postCue_withLick_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.cueAlign.withLick.CSP_postCue_trialAvg_Cspk = (std(mean(postCue_withLick_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(postCue_withLick_CSP,2))).*(1000./frameRateHz);
  n.cueAlign.woutLick.CSP_postCue_trialAvg_Cspk = sum(~gBlock2 & isnan(firstLickStart));
  m.cueAlign.woutLick.CSP_postCue_trialAvg_Cspk = mean(mean(postCue_woutLick_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.cueAlign.woutLick.CSP_postCue_trialAvg_Cspk = (std(mean(postCue_woutLick_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(postCue_woutLick_CSP,2))).*(1000./frameRateHz);
for c = 1:iDend
postCue_withBaseline_CSM(:,c) = CSMevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSM-statsWindow):baselineRange(1)-1+(peakIdx_baseCSM+statsWindow),c);
postCue_woutBaseline_CSM(:,c) = CSMevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSM-statsWindow):baselineRange(1)-1+(peakIdx_baseCSM+statsWindow),c);
postCue_withLick_CSM(:,c) = CSMevents_withLick_acrossTrials(postCueRange(1)-1+(peakIdx_postCue_withLick_CSM-statsWindow):postCueRange(1)-1+(peakIdx_postCue_withLick_CSM+statsWindow),c);
postCue_woutLick_CSM(:,c) = CSMevents_woutLick_acrossTrials(postCueRange(1)-1+(peakIdx_postCue_woutLick_CSM-statsWindow):postCueRange(1)-1+(peakIdx_postCue_woutLick_CSM+statsWindow),c);
end
[~, p.cueAlign.withLick.CSM_postCue_trialAvg_Cspk, ~, stats.cueAlign.withLick.CSM_postCue_trialAvg_Cspk] = ttest(mean(postCue_withBaseline_CSM,1,'omitnan'),mean(postCue_withLick_CSM,1,'omitnan'));
[~, p.cueAlign.woutLick.CSM_postCue_trialAvg_Cspk, ~, stats.cueAlign.woutLick.CSM_postCue_trialAvg_Cspk] = ttest(mean(postCue_woutBaseline_CSM,1,'omitnan'),mean(postCue_woutLick_CSM,1,'omitnan'));
[p.cueAlign.Lick.CSM_postCue_trialAvg_Cspk, ~, stats.cueAlign.Lick.CSM_postCue_trialAvg_Cspk] = ranksum(mean(postCue_woutLick_CSM,1,'omitnan'),mean(postCue_withLick_CSM,1,'omitnan'));
  n.cueAlign.withLick.CSM_postCue_trialAvg_Cspk = sum(gBlock2 & ~isnan(firstLickStart));
  m.cueAlign.withLick.CSM_postCue_trialAvg_Cspk = mean(mean(postCue_withLick_CSM,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.cueAlign.withLick.CSM_postCue_trialAvg_Cspk = (std(mean(postCue_withLick_CSM,1,'omitnan'),[],2,'omitnan')./sqrt(size(postCue_withLick_CSM,2))).*(1000./frameRateHz);
  n.cueAlign.woutLick.CSM_postCue_trialAvg_Cspk = sum(gBlock2 & isnan(firstLickStart));
  m.cueAlign.woutLick.CSM_postCue_trialAvg_Cspk = mean(mean(postCue_woutLick_CSM,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.cueAlign.woutLick.CSM_postCue_trialAvg_Cspk = (std(mean(postCue_woutLick_CSM,1,'omitnan'),[],2,'omitnan')./sqrt(size(postCue_woutLick_CSM,2))).*(1000./frameRateHz);
%post-rew
for c = 1:iDend
postRew_withBaseline_CSP(:,c) = CSPevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSP-statsWindow):baselineRange(1)-1+(peakIdx_baseCSP+statsWindow),c);
postRew_woutBaseline_CSP(:,c) = CSPevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSP-statsWindow):baselineRange(1)-1+(peakIdx_baseCSP+statsWindow),c);
postRew_withLick_CSP(:,c) = CSPevents_withLick_acrossTrials(postRewRange(1)-1+(peakIdx_postRew_withLick_CSP-statsWindow):postRewRange(1)-1+(peakIdx_postRew_withLick_CSP+statsWindow),c);
postRew_woutLick_CSP(:,c) = CSPevents_woutLick_acrossTrials(postRewRange(1)-1+(peakIdx_postRew_woutLick_CSP-statsWindow):postRewRange(1)-1+(peakIdx_postRew_woutLick_CSP+statsWindow),c);
end
[~, p.cueAlign.withLick.CSP_postRew_trialAvg_Cspk, ~, stats.cueAlign.withLick.CSP_postRew_trialAvg_Cspk] = ttest(mean(postRew_withBaseline_CSP,1,'omitnan'),mean(postRew_withLick_CSP,1,'omitnan'));
[~, p.cueAlign.woutLick.CSP_postRew_trialAvg_Cspk, ~, stats.cueAlign.woutLick.CSP_postRew_trialAvg_Cspk] = ttest(mean(postRew_woutBaseline_CSP,1,'omitnan'),mean(postRew_woutLick_CSP,1,'omitnan'));
[p.cueAlign.Lick.CSP_postRew_trialAvg_Cspk, ~, stats.cueAlign.Lick.CSP_postRew_trialAvg_Cspk] = ttest(mean(postRew_woutLick_CSP,1,'omitnan'),mean(postRew_withLick_CSP,1,'omitnan'));
  n.cueAlign.withLick.CSP_postRew_trialAvg_Cspk = sum(~gBlock2 & ~isnan(firstLickStart));
  m.cueAlign.withLick.CSP_postRew_trialAvg_Cspk = mean(mean(postRew_withLick_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.cueAlign.withLick.CSP_postRew_trialAvg_Cspk = (std(mean(postRew_withLick_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(postRew_withLick_CSP,2))).*(1000./frameRateHz);
  n.cueAlign.woutLick.CSP_postRew_trialAvg_Cspk = sum(~gBlock2 & isnan(firstLickStart));
  m.cueAlign.woutLick.CSP_postRew_trialAvg_Cspk = mean(mean(postRew_woutLick_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.cueAlign.woutLick.CSP_postRew_trialAvg_Cspk = (std(mean(postRew_woutLick_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(postRew_woutLick_CSP,2))).*(1000./frameRateHz);
for c = 1:iDend
postRew_withBaseline_CSM(:,c) = CSMevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSM-statsWindow):baselineRange(1)-1+(peakIdx_baseCSM+statsWindow),c);
postRew_woutBaseline_CSM(:,c) = CSMevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSM-statsWindow):baselineRange(1)-1+(peakIdx_baseCSM+statsWindow),c);
postRew_withLick_CSM(:,c) = CSMevents_withLick_acrossTrials(postRewRange(1)-1+(peakIdx_postRew_withLick_CSM-statsWindow):postRewRange(1)-1+(peakIdx_postRew_withLick_CSM+statsWindow),c);
postRew_woutLick_CSM(:,c) = CSMevents_woutLick_acrossTrials(postRewRange(1)-1+(peakIdx_postRew_woutLick_CSM-statsWindow):postRewRange(1)-1+(peakIdx_postRew_woutLick_CSM+statsWindow),c);
end
[~, p.cueAlign.withLick.CSM_postRew_trialAvg_Cspk, ~, stats.cueAlign.withLick.CSM_postRew_trialAvg_Cspk] = ttest(mean(postRew_withBaseline_CSM,1,'omitnan'),mean(postRew_withLick_CSM,1,'omitnan'));
[~, p.cueAlign.woutLick.CSM_postRew_trialAvg_Cspk, ~, stats.cueAlign.woutLick.CSM_postRew_trialAvg_Cspk] = ttest(mean(postRew_woutBaseline_CSM,1,'omitnan'),mean(postRew_woutLick_CSM,1,'omitnan'));
[p.cueAlign.Lick.CSM_postRew_trialAvg_Cspk, ~, stats.cueAlign.Lick.CSM_postRew_trialAvg_Cspk] = ranksum(mean(postRew_woutLick_CSM,1,'omitnan'),mean(postRew_withLick_CSM,1,'omitnan'));
  n.cueAlign.withLick.CSM_postRew_trialAvg_Cspk = sum(gBlock2 & ~isnan(firstLickStart));
  m.cueAlign.withLick.CSM_postRew_trialAvg_Cspk = mean(mean(postRew_withLick_CSM,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.cueAlign.withLick.CSM_postRew_trialAvg_Cspk = (std(mean(postRew_withLick_CSM,1,'omitnan'),[],2,'omitnan')./sqrt(size(postRew_withLick_CSM,2))).*(1000./frameRateHz);
  n.cueAlign.woutLick.CSM_postRew_trialAvg_Cspk = sum(gBlock2 & isnan(firstLickStart));
  m.cueAlign.woutLick.CSM_postRew_trialAvg_Cspk = mean(mean(postRew_woutLick_CSM,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.cueAlign.woutLick.CSM_postRew_trialAvg_Cspk = (std(mean(postRew_woutLick_CSM,1,'omitnan'),[],2,'omitnan')./sqrt(size(postRew_woutLick_CSM,2))).*(1000./frameRateHz);

  
tempFig=setFigure; hold on; 
subplot(2,1,1); hold on; 
 ylim([0 ymax]); xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
title('Average Cspk for CS+ trials with/out licks');
rectangle('Position',[tt(1,baselineRange(1))/1000,-1+.1,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postCueRange(1))/1000,-1+.1,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postRewRange(1))/1000,-1+.1,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSPevents_withLick_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSPevents_withLick_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSPevents_withLick_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSPevents_woutLick_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSPevents_woutLick_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSPevents_woutLick_acrossTrials,2)).*(1000./frameRateHz),'lineProps','b');
 
subplot(2,1,2); hold on;
 ylim([0 ymax]); xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
 title('Average Cspk for CS- trials with/out licks');
rectangle('Position',[tt(1,baselineRange(1))/1000,-1+.1,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postCueRange(1))/1000,-1+.1,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
rectangle('Position',[tt(1,postRewRange(1))/1000,-1+.1,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSMevents_withLick_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSMevents_withLick_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSMevents_withLick_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSMevents_woutLick_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSMevents_woutLick_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSMevents_woutLick_acrossTrials,2)).*(1000./frameRateHz),'lineProps','m');
saveas(tempFig, [output_fn 'stats\withWoutLicking.pdf'],'pdf');
%% ttest    | first-lick avg vs. baseline avg 
preLick_frames=15; preLickRange=(1:bx_response_frames)+(preLick_frames-bx_response_frames);
baselineRange=prewin_frames-response_frames+1:prewin_frames;
postCueRange=prewin_frames+1:prewin_frames+response_frames;
postRewRange=prewin_frames+rewDelay_frames+1:prewin_frames+rewDelay_frames+response_frames;

% offset start frame by 1, because we need to take a +/-1 frame window
% around the target; script breaks if the peak is on frame 1
  [peak_firstLickActCSP, peakIdx_firstLickActCSP] = max(mean(first_lickAlignEventsCSP_acrossTrials(preLickRange,:),2,'omitnan'),[],1);
  [peak_postCueLickActCSP, peakIdx_postCueLickActCSP] = max(mean(postCue_lickAlignEventsCSP_acrossTrials(preLickRange,:),2,'omitnan'),[],1);
  [peak_postRewLickActCSP, peakIdx_postRewLickActCSP] = max(mean(postRew_lickAlignEventsCSP_acrossTrials(preLickRange,:),2,'omitnan'),[],1);
  [peak_baseCSP, peakIdx_baseCSP] = max(mean(CSPevents_acrossTrials(baselineRange,:),2,'omitnan'),[],1);
  [peak_firstLickActCSM, peakIdx_firstLickActCSM] = max(mean(first_lickAlignEventsCSM_acrossTrials(preLickRange,:),2,'omitnan'),[],1);
  [peak_postCueLickActCSM, peakIdx_postCueLickActCSM] = max(mean(postCue_lickAlignEventsCSM_acrossTrials(preLickRange,:),2,'omitnan'),[],1);
  [peak_postRewLickActCSM, peakIdx_postRewLickActCSM] = max(mean(postRew_lickAlignEventsCSM_acrossTrials(preLickRange,:),2,'omitnan'),[],1);
  [peak_baseCSM, peakIdx_baseCSM] = max(mean(CSMevents_acrossTrials(baselineRange,:),2,'omitnan'),[],1);

  % counts for lick align
n.lickAlign.CSP_first_trialAvg_count = sum(first_lickAlignCSP_acrossTrials); 
n.lickAlign.CSM_first_trialAvg_count = sum(first_lickAlignCSM_acrossTrials);
n.lickAlign.CSP_firstPostCue_trialAvg_count = sum(firstPostCue_lickAlignCSP_acrossTrials);
n.lickAlign.CSM_firstPostCue_trialAvg_count = sum(firstPostCue_lickAlignCSM_acrossTrials);
n.lickAlign.CSP_firstPostRew_trialAvg_count = sum(firstPostRew_lickAlignCSP_acrossTrials);
n.lickAlign.CSM_firstPostRew_trialAvg_count = sum(firstPostRew_lickAlignCSM_acrossTrials);
% difference in [] between first lick and baseline
for c = 1:iDend
lickAlign_baseline_CSP(:,c) = CSPevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSP-statsWindow):baselineRange(1)-1+(peakIdx_baseCSP+statsWindow),c);
lickAlign_activity_CSP(:,c) = first_lickAlignEventsCSP_acrossTrials(preLickRange(1)-1+(peakIdx_firstLickActCSP-statsWindow):preLickRange(1)-1+(peakIdx_firstLickActCSP+statsWindow),c);
end
[~, p.lickAlign.CSP_first_trialAvg_Cspk, ~, stats.lickAlign.CSP_first_trialAvg_Cspk] = ttest(mean(lickAlign_baseline_CSP,1,'omitnan'),mean(lickAlign_activity_CSP,1,'omitnan'));
  m.lickAlign.CSP_first_trialAvg_Cspk = mean(mean(lickAlign_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickAlign.CSP_first_trialAvg_Cspk = std(mean(lickAlign_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(lickAlign_activity_CSP,2)).*(1000./frameRateHz);
for c = 1:iDend % -statsWindow
lickAlign_baseline_CSM(:,c) = CSMevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSM-statsWindow):baselineRange(1)-1+(peakIdx_baseCSM+statsWindow),c);
lickAlign_activity_CSM(:,c) = first_lickAlignEventsCSM_acrossTrials(preLickRange(1)-1+(peakIdx_firstLickActCSM-statsWindow):preLickRange(1)-1+(peakIdx_firstLickActCSM+statsWindow),c);
end
[~, p.lickAlign.CSM_first_trialAvg_Cspk, ~, stats.lickAlign.CSM_first_trialAvg_Cspk] = ttest(mean(lickAlign_baseline_CSM,1,'omitnan'),mean(lickAlign_activity_CSM,1,'omitnan'));
  m.lickAlign.CSM_first_trialAvg_Cspk = mean(mean(lickAlign_activity_CSM,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickAlign.CSM_first_trialAvg_Cspk = std(mean(lickAlign_activity_CSM,1,'omitnan'),[],2,'omitnan')./sqrt(size(lickAlign_activity_CSM,2)).*(1000./frameRateHz);

% difference in activity between first lick [POST REWARD] and baseline
for c = 1:iDend
postRewLick_baseline_CSP(:,c) = CSPevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSP-statsWindow):baselineRange(1)-1+(peakIdx_baseCSP+statsWindow),c);
postRewLick_activity_CSP(:,c) = postRew_lickAlignEventsCSP_acrossTrials(preLickRange(1)-1+(peakIdx_postRewLickActCSP-statsWindow):preLickRange(1)-1+(peakIdx_postRewLickActCSP+statsWindow),c);
end
[~, p.lickAlign.CSP_postRew_trialAvg_Cspk, ~, stats.lickAlign.CSP_postRew_trialAvg_Cspk] = ttest(mean(postRewLick_baseline_CSP,1,'omitnan'),mean(postRewLick_activity_CSP,1,'omitnan'));
  m.lickAlign.CSP_postRew_trialAvg_Cspk = mean(mean(postRewLick_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickAlign.CSP_postRew_trialAvg_Cspk = std(mean(postRewLick_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(postRewLick_activity_CSP,2)).*(1000./frameRateHz);
for c = 1:iDend
postRewLick_baseline_CSM(:,c) = CSMevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSM-statsWindow):baselineRange(1)-1+(peakIdx_baseCSM+statsWindow),c);
postRewLick_activity_CSM(:,c) = postRew_lickAlignEventsCSM_acrossTrials(preLickRange(1)-1+(peakIdx_postRewLickActCSM-statsWindow):preLickRange(1)-1+(peakIdx_postRewLickActCSM+statsWindow),c);
end
[~, p.lickAlign.CSM_postRew_trialAvg_Cspk, ~, stats.lickAlign.CSM_postRew_trialAvg_Cspk] = ttest(mean(postRewLick_baseline_CSM,1,'omitnan'),mean(postRewLick_activity_CSM,1,'omitnan'));
  m.lickAlign.CSM_postRew_trialAvg_Cspk = mean(mean(postRewLick_activity_CSM,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickAlign.CSM_postRew_trialAvg_Cspk = std(mean(postRewLick_activity_CSM,1,'omitnan'),[],2,'omitnan')./sqrt(size(postRewLick_activity_CSM,2)).*(1000./frameRateHz);
  
% difference in activity between first lick [POST CUE] and baseline
for c = 1:iDend
postCueLick_baseline_CSP(:,c) = CSPevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSP-statsWindow):baselineRange(1)-1+(peakIdx_baseCSP+statsWindow),c);
postCueLick_activity_CSP(:,c) = postCue_lickAlignEventsCSP_acrossTrials(preLickRange(1)-1+(peakIdx_postCueLickActCSP-statsWindow):preLickRange(1)-1+(peakIdx_postCueLickActCSP+statsWindow),c);
end
[~, p.lickAlign.CSP_postCue_trialAvg_Cspk, ~, stats.lickAlign.CSP_postCue_trialAvg_Cspk] = ttest(mean(postCueLick_baseline_CSP,1,'omitnan'),mean(postCueLick_activity_CSP,1,'omitnan'));
  m.lickAlign.CSP_postCue_trialAvg_Cspk = mean(mean(postCueLick_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickAlign.CSP_postCue_trialAvg_Cspk = std(mean(postCueLick_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(postCueLick_activity_CSP,2)).*(1000./frameRateHz);
for c = 1:iDend
postCueLick_baseline_CSM(:,c) = CSMevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSM-statsWindow):baselineRange(1)-1+(peakIdx_baseCSM+statsWindow),c);
postCueLick_activity_CSM(:,c) = postCue_lickAlignEventsCSM_acrossTrials(preLickRange(1)-1+(peakIdx_postCueLickActCSM-statsWindow):preLickRange(1)-1+(peakIdx_postCueLickActCSM+statsWindow),c);
end
[~, p.lickAlign.CSM_postCue_trialAvg_Cspk, ~, stats.lickAlign.CSM_postCue_trialAvg_Cspk] = ttest(mean(postCueLick_baseline_CSM,1,'omitnan'),mean(postCueLick_activity_CSM,1,'omitnan'));
  m.lickAlign.CSM_postCue_trialAvg_Cspk = mean(mean(postCueLick_activity_CSM,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickAlign.CSM_postCue_trialAvg_Cspk = std(mean(postCueLick_activity_CSM,1,'omitnan'),[],2,'omitnan')./sqrt(size(postCueLick_activity_CSM,2)).*(1000./frameRateHz);

 save(fullfile([output_fn 'lickAlign.mat']), 'expt','nsize','lickAlign_baseline_CSP','lickAlign_activity_CSP','postCueLick_activity_CSP','postRewLick_activity_CSP','lickAlign_baseline_CSM','lickAlign_activity_CSM','postCueLick_activity_CSM','postRewLick_activity_CSM','peakIdx_firstLickActCSP','peakIdx_postCueLickActCSP','peakIdx_postRewLickActCSP','peakIdx_baseCSP','peakIdx_firstLickActCSM','peakIdx_postCueLickActCSM','peakIdx_postRewLickActCSM','peakIdx_baseCSM');

%%    analysis for lick modulation (CS+ - CS- / CS+ + CS-)
  allCell_peakLickCSP = mean(mean(first_lickAlignEventsCSP_acrossTrials(preLickRange(1)-1+(peakIdx_firstLickActCSP-statsWindow):preLickRange(1)-1+(peakIdx_firstLickActCSP+statsWindow),:,:),3,'omitnan'),1,'omitnan').*(1000./frameRateHz);
  allCell_peakLickCSM = mean(mean(first_lickAlignEventsCSM_acrossTrials(preLickRange(1)-1+(peakIdx_firstLickActCSM-statsWindow):preLickRange(1)-1+(peakIdx_firstLickActCSM+statsWindow),:,:),3,'omitnan'),1,'omitnan').*(1000./frameRateHz);
  
  % Adjust alpha for multiple comparisons using Bonferroni correction
  adjusted_alpha = 0.05 / size(animalTrials,2);
  % Revised loop to use adjusted alpha
  for i = 1:iDend
    [~, p.modulationLickIndex(i), ~] = ttest(mean(first_lickAlignEventsCSP_acrossTrials(preLickRange(1)-1+(peakIdx_firstLickActCSP-statsWindow):preLickRange(1)-1+(peakIdx_firstLickActCSP+statsWindow),i,:),3,'omitnan'), mean(first_lickAlignEventsCSM_acrossTrials(preLickRange(1)-1+(peakIdx_firstLickActCSM-statsWindow):preLickRange(1)-1+(peakIdx_firstLickActCSM+statsWindow),i,:),3,'omitnan'));
  end
  modulationIndex = (allCell_peakLickCSP-allCell_peakLickCSM)./(allCell_peakLickCSP+allCell_peakLickCSM);
  fav_LickCSP_indices = find(allCell_peakLickCSP > allCell_peakLickCSM & modulationIndex == 1);
  fav_LickCSM_indices = find(allCell_peakLickCSP < allCell_peakLickCSM & modulationIndex == -1);
  no_preference_indices = find(p.modulationLickIndex == 0 & modulationIndex ~= 1 | modulationIndex ~= -1);
  
  % Normalize modulation index for color mapping: range from 0 (red) to 1 (black)
  modIndexNormalized = (modulationIndex - min(modulationIndex)) / (max(modulationIndex) - min(modulationIndex));
   
    tempFig=setFigure;
    hold on; 
    % Loop through each cell to plot
    for i = 1:iDend
        % Determine color: Interpolate between red (1, 0, 0) and black (0, 0, 0)
        color = [1, 0, 0] * (1 - modIndexNormalized(i)); % Reduce red component
           % Plotting CSP vs. CSM with modulation-based color
        if ~isnan(modulationIndex(i))
            scatter(allCell_peakLickCSP(i), allCell_peakLickCSM(i), 40, 'filled', 'MarkerFaceColor', color);
        end
    end
    % Reference lines
    xLimits = [0 max(allCell_peakLickCSP)];
    yLimits = [0 max(allCell_peakLickCSM)];
    plotLimits = [0 max(xLimits(2), yLimits(2))]; % Ensuring lines extend to the edge of data

    % Calculate and plot 33-degree and 66-degree lines
    % 33 degrees
    theta33 = deg2rad(33);
    slope33 = tan(theta33);
    intercept33 = median(allCell_peakLickCSM) - slope33 * median(allCell_peakLickCSP);
    y33 = slope33 * plotLimits;
    
    % 66 degrees
    theta66 = deg2rad(66);
    slope66 = tan(theta66);
    intercept66 = median(allCell_peakLickCSM) - slope66 * median(allCell_peakLickCSP);
    y66 = slope66 * plotLimits;
    
    % Add the lines to the plot
    plot(plotLimits, y33, 'k--');
    plot(plotLimits, y66, 'r--');
    
    % Categorizing cells based on the zones they fall into
    thresholdCSP = 1;%m.cueAlign.CSP_baseline_trialAvg_Cspk + 10*s.cueAlign.CSP_baseline_trialAvg_Cspk;
    thresholdCSM = 1;%m.cueAlign.CSM_baseline_trialAvg_Cspk + 10*s.cueAlign.CSM_baseline_trialAvg_Cspk; 
    %rectangle('Position', [0, 0, thresholdCSP, thresholdCSM], 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
    zoneLick = zeros(iDend, 1);  % 1 for below 33 degrees, 2 for between 33 and 66, 3 for above 66 degrees
    for i = 1:iDend
        if allCell_peakLickCSP(i) < thresholdCSP && allCell_peakLickCSM(i) < thresholdCSM
        zoneLick(i) = 4; % Not significantly above baseline
       
        else
        yVal33 = slope33 * allCell_peakLickCSP(i);
        yVal66 = slope66 * allCell_peakLickCSP(i);
        
        if allCell_peakLickCSM(i) < yVal33
            zoneLick(i) = 1;  % Below 33 degrees
        elseif allCell_peakLickCSM(i) < yVal66
            zoneLick(i) = 2;  % Between 33 and 66 degrees
        else
            zoneLick(i) = 3;  % Above 66 degrees
        end
        end
    end
    
    % Count cells in each zone
    cellPref_lickCSP = find(zoneLick == 1);
    cellPref_lickNone = find(zoneLick == 2);
    cellPref_lickCSM = find(zoneLick == 3);
    cellPref_lickInsig = find(zoneLick == 4);
   
    % Enhance plot
    xlabel('Peak CSP Response');
    ylabel('Peak CSM Response');
    title('CSP vs. CSM Responses with Angle Zones');
    grid on; ylim([-.50 ymax*3]); xlim([-.50 ymax*3])
    hold off;
     saveas(tempFig, [output_fn 'stats\lickAlign_cuePreferenceApprox_peakSpkEvents_scatter.pdf'],'pdf');

    tempFig = setFigure; hold on;
    titleArray={'CS+ preference','no preference','CS- preference','unresponsive'};
for c = 1:4     
subplot(4,1,c); hold on;
 ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title(titleArray{1,c});
if c<5; rectangle('Position',[tt(1,36+preLickRange(1)-1)/1000,0+.1,(1000./frameRateHz)*response_frames/1000,ymax*3],'EdgeColor',highlight,'FaceColor',highlight); elseif c>4; rectangle('Position',[tt(1,postRewRange(1))/1000,0+.1,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight); end
shadedErrorBar_CV(tt(1,36:80)/1000,mean(first_lickAlignEventsCSP_acrossTrials(:,zoneLick==c,:),2,'omitnan').*(1000./frameRateHz),std(first_lickAlignEventsCSP_acrossTrials(:,zone==c,:),[],2,'omitnan')./sqrt(iDend).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:80)/1000,mean(first_lickAlignEventsCSM_acrossTrials(:,zoneLick==c,:),2,'omitnan').*(1000./frameRateHz),std(first_lickAlignEventsCSM_acrossTrials(:,zone==c,:),[],2,'omitnan')./sqrt(iDend).*(1000./frameRateHz),'lineProps','r');
end
    saveas(tempFig, [output_fn 'stats\lickAlign_cuePreferenceApprox_cSpkEvents_PSTHs.pdf'],'pdf');

    %% ANOVA for lick align comparisons

        cueID={'CS+','CS-'}; nCue=size(cueID,2); cells = 1:nIC; %nTrialsSplit=[size(firstPostRew_lick_spikeProbability_CSPtrials_split{1,1},3);size(firstLick_spikeProbability_CSMtrials_split{1,1},3)]; trialsSplit{1,1}=1:nTrialsSplit(1,1); trialsSplit{1,2}=1:nTrialsSplit(2,1);
        lickAlign_cSpk_array{1,1} = mean(lickAlign_activity_CSP,1,'omitnan'); lickAlign_cSpk_array{1,2} = mean(lickAlign_activity_CSM,1,'omitnan');
        dataArray = NaN(nCue, nIC);
%        dataArray(2, :, :, :) = NaN(nCue, nPhase, nIC);% Initialize the 3D array
        for i = 1:nCue
          dataArray(i, :) = lickAlign_cSpk_array{1,i};
        end

%      
numCues = 2; % Assuming there are 2 cues
numCells = nIC;
numMice = size(expt,2); % Number of mice

% Initialize cell arrays to store data for each factor
mouseData = cell(numMice, 1); for i = 1:numMice; mouseID{i} = sprintf('mouse %d',i); end
cueData = cell(numCues, 1);

cellIdx = []; 
for i = 1:size(expt,2)
cellIdx = [cellIdx; ones(nsize(1,i),1)*i];
end
wholeMouseData=[];wholeCueData=[];wholeCellData=[];
lickAlign_wholeResponseData=[];
for i = 1:size(cellIdx,1)
    for c=1:numCues
wholeMouseData(i,:) = cellIdx(i,1);
lickAlign_wholeResponseData(i,:) = [dataArray(1,i), dataArray(2,i)];
wholeCueData = [wholeCueData, c];
cellID(i,:) = i;
wholeCellData = [wholeCellData i];
        
    end
end

dataForTable=[lickAlign_wholeResponseData(:,1),lickAlign_wholeResponseData(:,2)];
t = table(categorical(cellID),lickAlign_wholeResponseData(:,1),lickAlign_wholeResponseData(:,2),VariableNames=["cell","response1","response2"]);
% Create a table reflecting the within subject factors 
cue = categorical(repmat([1 2]', 1, 1));
within = table(cue, 'VariableNames', {'cue'});
% 3. Call fitrm with the modified within design.
rm = fitrm(t,'response1-response2~1','WithinDesign',within);
[anovaTable_lickAlign] = ranova(rm, 'WithinModel','cue');

anovaTable_lickAlign = movevars(anovaTable_lickAlign, "DF", "Before", "SumSq");
writetable(anovaTable_lickAlign,[output_fn 'stats\lickAlign_two-wayANOVA.xlsx']);
disp('lickAlign')
disp(anovaTable_lickAlign)

%% ANOVA    | comparison between predictive and reactive licking
    predictiveReactiveID = {'predictive','reactive'}; cueID={'CS+','CS-'}; nCue=size(cueID,2); nSesh=size(predictiveReactiveID,2); cells = 1:nIC; %nTrialsSplit=[size(firstPostRew_lick_spikeProbability_CSPtrials_split{1,1},3);size(firstLick_spikeProbability_CSMtrials_split{1,1},3)]; trialsSplit{1,1}=1:nTrialsSplit(1,1); trialsSplit{1,2}=1:nTrialsSplit(2,1);
predictiveReactiveComp_lickAlign_cSpk_array{1,1} = mean(postCueLick_activity_CSP,1,'omitnan'); predictiveReactiveComp_lickAlign_cSpk_array{1,3} = mean(postRewLick_activity_CSP,1,'omitnan');
predictiveReactiveComp_lickAlign_cSpk_array{1,2} = mean(postCueLick_activity_CSM,1,'omitnan'); predictiveReactiveComp_lickAlign_cSpk_array{1,4} = mean(postRewLick_activity_CSM,1,'omitnan');
dataArray = NaN(nCue, (nIC));
for i = 1:nCue+nSesh
  dataArray(i, :) = predictiveReactiveComp_lickAlign_cSpk_array{1,i};
end

%      
numCues = 2; % Assuming there are 2 cues
numCells = nIC+nIC;
numMice = size(expt,2); % Number of mice

% Initialize cell arrays to store data for each factor
mouseData = cell(numMice, 1); for i = 1:numMice; mouseID{i} = sprintf('mouse %d',i); end
cueData = cell(numCues, 1);

cellIdx = []; 
for i = 1:size(expt,2)
cellIdx = [cellIdx; ones(nsize(1,i),1)*i];
end
wholeMouseData=[];wholeCueData=[];wholeCellData=[];
predictiveReactiveComp_lickAlign_wholeResponseData=[];
for i = 1:size(cellIdx,1)
    for c=1:numCues
wholeMouseData(i,:) = 1;%cellIdx(i,1);
predictiveReactiveComp_lickAlign_wholeResponseData(i,:) = [dataArray(:,i)];
wholeCueData = [wholeCueData, c];
cellID(i,:) = i;
wholeCellData = [wholeCellData i];   
    end
end

dataForTable=[predictiveReactiveComp_lickAlign_wholeResponseData(:,1),predictiveReactiveComp_lickAlign_wholeResponseData(:,2)];
% Create a table reflecting the within subject factors 
cue = categorical([1 2 1 2]');
time = categorical([1 1 2 2]');
% 3. Call fitrm with the modified within design.
t = table(cellID, dataArray(1,:)', dataArray(2,:)', dataArray(3,:)', dataArray(4,:)', ...
    'VariableNames', ["cell", "response1", "response2", "response3", "response4"]);
t = rmmissing(t);
within = table(cue, time, 'VariableNames', {'cue', 'time'});
rm = fitrm(t, 'response1-response4~1', 'WithinDesign', within);
anovaTable_predictiveReactiveComp_lickAlign = ranova(rm, 'WithinModel','cue*time');

anovaTable_predictiveReactiveComp_lickAlign = movevars(anovaTable_predictiveReactiveComp_lickAlign, "DF", "Before", "SumSq");
writetable(anovaTable_predictiveReactiveComp_lickAlign,[output_fn 'stats\predictiveReactiveComp_lickAlign_two-wayANOVA.xlsx']);
disp('predictive vs reactive comp_lickAlign')
disp(anovaTable_predictiveReactiveComp_lickAlign)

% CS+
[~, p.lickAlign.CSP_predVreac_trialAvg_Cspk, ~, stats.lickAlign.CSP_predVreac_trialAvg_Cspk] = ttest(mean(postCueLick_activity_CSP,1,'omitnan'),mean(postRewLick_activity_CSP,1,'omitnan'));
% CS-
[~, p.lickAlign.CSM_predVreac_trialAvg_Cspk, ~, stats.lickAlign.CSM_predVreac_trialAvg_Cspk] = ttest(mean(postCueLick_activity_CSM,1,'omitnan'),mean(postRewLick_activity_CSM,1,'omitnan'));

%% ANOVA    | if PL compare pre and postlearning responses
if seshType == 3 || seshType == 5
% load data
    if seshType == 3
    naive_fn = ['A:\home\carlo\RC\analysis\2P\groupwiseAlign_Naive_' sesh2 '_threshold' num2str(thresholdDeco) '_240317\'];
    naiveGroupData = load(fullfile(naive_fn, 'lickAlign.mat'));
    elseif seshType == 5
    naive_fn = ['A:\home\carlo\RC\analysis\2P\groupwiseAlign_RuleFlipNaive_' sesh2 '_threshold' num2str(thresholdDeco) '_240317\'];
    naiveGroupData = load(fullfile(naive_fn, 'lickAlign.mat'));
    end

% ttest for comp between pre and post learning first lick
seshComp_lickAlign_cSpk_array_parm{1,1} = mean(naiveGroupData.lickAlign_activity_CSP,1,'omitnan');
seshComp_lickAlign_cSpk_array_parm{1,2} = mean(naiveGroupData.lickAlign_activity_CSM,1,'omitnan');
seshComp_lickAlign_cSpk_array_parm{2,1} = mean(lickAlign_activity_CSP,1,'omitnan'); %second row
seshComp_lickAlign_cSpk_array_parm{2,2} = mean(lickAlign_activity_CSM,1,'omitnan');

[~,p.seshComp_lickAlign_cSpk_array_CSP,~,stats.seshComp_lickAlign_cSpk_array_CSP] = ttest2(seshComp_lickAlign_cSpk_array_parm{1,1},seshComp_lickAlign_cSpk_array_parm{2,1},'Vartype','equal');
[~,p.seshComp_lickAlign_cSpk_array_CSM,~,stats.seshComp_lickAlign_cSpk_array_CSM] = ttest2(seshComp_lickAlign_cSpk_array_parm{1,2},seshComp_lickAlign_cSpk_array_parm{2,2},'Vartype','equal');


% ANOVA for comparison between pre and post learning first lick  
seshID = {'naive','postlearning'}; cueID={'CS+','CS-'}; nCue=size(cueID,2); nSesh=size(seshID,2); cells = 1:nIC; %nTrialsSplit=[size(firstPostRew_lick_spikeProbability_CSPtrials_split{1,1},3);size(firstLick_spikeProbability_CSMtrials_split{1,1},3)]; trialsSplit{1,1}=1:nTrialsSplit(1,1); trialsSplit{1,2}=1:nTrialsSplit(2,1);
seshComp_lickAlign_cSpk_array{1,1} = [mean(naiveGroupData.lickAlign_activity_CSP,1,'omitnan'), mean(lickAlign_activity_CSP,1,'omitnan')];
seshComp_lickAlign_cSpk_array{1,2} = [mean(naiveGroupData.lickAlign_activity_CSM,1,'omitnan'), mean(lickAlign_activity_CSM,1,'omitnan')];
dataArray = NaN(nCue, sum([nIC,size(naiveGroupData.lickAlign_activity_CSP,2)]));
%        dataArray(2, :, :, :) = NaN(nCue, nPhase, nIC);% Initialize the 3D array
for i = 1:nCue
  dataArray(i, :) = seshComp_lickAlign_cSpk_array{1,i};
end

%      
numCues = 2; % Assuming there are 2 cues
numCells = sum([nIC,size(naiveGroupData.lickAlign_activity_CSP,2)]);
numMice = size(expt,2); % Number of mice

% Initialize cell arrays to store data for each factor
mouseData = cell(numMice, 1); for i = 1:numMice; mouseID{i} = sprintf('mouse %d',i); end
cueData = cell(numCues, 1);

cellIdx = []; 
for i = 1:size(expt,2)
cellIdx = [cellIdx; ones(nsize(1,i),1)*i];
end
for i = 1:size(naiveGroupData.expt,2)
cellIdx = [cellIdx; ones(naiveGroupData.nsize(1,i),1)*i+size(expt,2)];
end
wholeMouseData=[];wholeCueData=[];wholeCellData=[];
seshComp_lickAlign_wholeResponseData=[];
for i = 1:size(cellIdx,1)
    for c=1:numCues
wholeMouseData(i,:) = 1;%cellIdx(i,1);
seshComp_lickAlign_wholeResponseData(i,:) = [dataArray(1,i), dataArray(2,i)];
wholeCueData = [wholeCueData, c];
cellID(i,:) = i;
wholeCellData = [wholeCellData i];
        
    end
end

dataForTable=[seshComp_lickAlign_wholeResponseData(:,1),seshComp_lickAlign_wholeResponseData(:,2)];
% Create a table reflecting the within subject factors 
cue = categorical([ones(1,sum(size(naiveGroupData.lickAlign_activity_CSP,2)+nIC)) ,2*ones(1,sum(size(naiveGroupData.lickAlign_activity_CSP,2)+nIC))]');
sesh = categorical([ones(1, size(naiveGroupData.lickAlign_activity_CSP,2)), 2*ones(1, nIC)]');
% 3. Call fitrm with the modified within design.
t = table(categorical([cellID]),sesh,dataForTable(:,1),dataForTable(:,2),VariableNames=["cell","sesh","response1","response2"]);
within = table([1;2], 'VariableNames', {'cue'});
rm = fitrm(t,'response1-response2~sesh', 'WithinDesign', within);
[anovaTable_seshComp_lickAlign] = ranova(rm);

anovaTable_seshComp_lickAlign = movevars(anovaTable_seshComp_lickAlign, "DF", "Before", "SumSq");
writetable(anovaTable_seshComp_lickAlign,[output_fn 'stats\seshComp_lickAlign_two-wayANOVA.xlsx']);
disp('seshComp_lickAlign')
disp(anovaTable_seshComp_lickAlign)

% ANOVA for comparison between pre and post learning - PREDICTIVE LICK  
seshID = {'naive','postlearning'}; cueID={'CS+','CS-'}; nCue=size(cueID,2); nSesh=size(seshID,2); cells = 1:nIC; %nTrialsSplit=[size(firstPostRew_lick_spikeProbability_CSPtrials_split{1,1},3);size(firstLick_spikeProbability_CSMtrials_split{1,1},3)]; trialsSplit{1,1}=1:nTrialsSplit(1,1); trialsSplit{1,2}=1:nTrialsSplit(2,1);
seshComp_postCueLick_cSpk_array{1,1} = [mean(naiveGroupData.postCueLick_activity_CSP,1,'omitnan'), mean(postCueLick_activity_CSP,1,'omitnan')];
seshComp_postCueLick_cSpk_array{1,2} = [mean(naiveGroupData.postCueLick_activity_CSM,1,'omitnan'), mean(postCueLick_activity_CSM,1,'omitnan')];
dataArray = NaN(nCue, sum([nIC,size(naiveGroupData.postCueLick_activity_CSP,2)]));
%        dataArray(2, :, :, :) = NaN(nCue, nPhase, nIC);% Initialize the 3D array
for i = 1:nCue
  dataArray(i, :) = seshComp_postCueLick_cSpk_array{1,i};
end

%      
numCues = 2; % Assuming there are 2 cues
numCells = sum([nIC,size(naiveGroupData.postCueLick_activity_CSP,2)]);
numMice = size(expt,2); % Number of mice

% Initialize cell arrays to store data for each factor
mouseData = cell(numMice, 1); for i = 1:numMice; mouseID{i} = sprintf('mouse %d',i); end
cueData = cell(numCues, 1);

cellIdx = []; 
for i = 1:size(expt,2)
cellIdx = [cellIdx; ones(nsize(1,i),1)*i];
end
for i = 1:size(naiveGroupData.expt,2)
cellIdx = [cellIdx; ones(naiveGroupData.nsize(1,i),1)*i+size(expt,2)];
end
wholeMouseData=[];wholeCueData=[];wholeCellData=[];
seshComp_postCueLick_wholeResponseData=[];
for i = 1:size(cellIdx,1)
    for c=1:numCues
wholeMouseData(i,:) = 1;%cellIdx(i,1);
seshComp_postCueLick_wholeResponseData(i,:) = [dataArray(1,i), dataArray(2,i)];
wholeCueData = [wholeCueData, c];
cellID(i,:) = i;
wholeCellData = [wholeCellData i];
        
    end
end

dataForTable=[seshComp_postCueLick_wholeResponseData(:,1),seshComp_postCueLick_wholeResponseData(:,2)];
% Create a table reflecting the within subject factors 
cue = categorical([ones(1,sum(size(naiveGroupData.postCueLick_activity_CSP,2)+nIC)) ,2*ones(1,sum(size(naiveGroupData.postCueLick_activity_CSP,2)+nIC))]');
sesh = categorical([ones(1, size(naiveGroupData.postCueLick_activity_CSP,2)), 2*ones(1, nIC)]');
% 3. Call fitrm with the modified within design.
t = table(categorical([cellID]),sesh,dataForTable(:,1),dataForTable(:,2),VariableNames=["cell","sesh","response1","response2"]);
within = table([1;2], 'VariableNames', {'cue'});
rm = fitrm(t,'response1-response2~sesh', 'WithinDesign', within);
[anovaTable_seshComp_postCueLick] = ranova(rm);

anovaTable_seshComp_postCueLick = movevars(anovaTable_seshComp_postCueLick, "DF", "Before", "SumSq");
writetable(anovaTable_seshComp_postCueLick,[output_fn 'stats\seshComp_postCueLick_two-wayANOVA.xlsx']);
disp('seshComp_postCueLick')
disp(anovaTable_seshComp_postCueLick)

% ANOVA for comparison between pre and post learning - REACTIVE LICK 
seshID = {'naive','postlearning'}; cueID={'CS+','CS-'}; nCue=size(cueID,2); nSesh=size(seshID,2); cells = 1:nIC; %nTrialsSplit=[size(firstPostRew_lick_spikeProbability_CSPtrials_split{1,1},3);size(firstLick_spikeProbability_CSMtrials_split{1,1},3)]; trialsSplit{1,1}=1:nTrialsSplit(1,1); trialsSplit{1,2}=1:nTrialsSplit(2,1);
seshComp_postRewLick_cSpk_array{1,1} = [mean(naiveGroupData.postRewLick_activity_CSP,1,'omitnan'), mean(postRewLick_activity_CSP,1,'omitnan')];
seshComp_postRewLick_cSpk_array{1,2} = [mean(naiveGroupData.postRewLick_activity_CSM,1,'omitnan'), mean(postRewLick_activity_CSM,1,'omitnan')];
dataArray = NaN(nCue, sum([nIC,size(naiveGroupData.postRewLick_activity_CSP,2)]));
%        dataArray(2, :, :, :) = NaN(nCue, nPhase, nIC);% Initialize the 3D array
for i = 1:nCue
  dataArray(i, :) = seshComp_postRewLick_cSpk_array{1,i};
end

%      
numCues = 2; % Assuming there are 2 cues
numCells = sum([nIC,size(naiveGroupData.postRewLick_activity_CSP,2)]);
numMice = size(expt,2); % Number of mice

% Initialize cell arrays to store data for each factor
mouseData = cell(numMice, 1); for i = 1:numMice; mouseID{i} = sprintf('mouse %d',i); end
cueData = cell(numCues, 1);

cellIdx = []; 
for i = 1:size(expt,2)
cellIdx = [cellIdx; ones(nsize(1,i),1)*i];
end
for i = 1:size(naiveGroupData.expt,2)
cellIdx = [cellIdx; ones(naiveGroupData.nsize(1,i),1)*i+size(expt,2)];
end
wholeMouseData=[];wholeCueData=[];wholeCellData=[];
seshComp_postRewLick_wholeResponseData=[];
for i = 1:size(cellIdx,1)
    for c=1:numCues
wholeMouseData(i,:) = 1;%cellIdx(i,1);
seshComp_postRewLick_wholeResponseData(i,:) = [dataArray(1,i), dataArray(2,i)];
wholeCueData = [wholeCueData, c];
cellID(i,:) = i;
wholeCellData = [wholeCellData i];
        
    end
end

dataForTable=[seshComp_postRewLick_wholeResponseData(:,1),seshComp_postRewLick_wholeResponseData(:,2)];
% Create a table reflecting the within subject factors 
cue = categorical([ones(1,sum(size(naiveGroupData.postRewLick_activity_CSP,2)+nIC)) ,2*ones(1,sum(size(naiveGroupData.postRewLick_activity_CSP,2)+nIC))]');
sesh = categorical([ones(1, size(naiveGroupData.postRewLick_activity_CSP,2)), 2*ones(1, nIC)]');
% 3. Call fitrm with the modified within design.
t = table(categorical([cellID]),sesh,dataForTable(:,1),dataForTable(:,2),VariableNames=["cell","sesh","response1","response2"]);
within = table([1;2], 'VariableNames', {'cue'});
rm = fitrm(t,'response1-response2~sesh', 'WithinDesign', within);
[anovaTable_seshComp_postRewLick] = ranova(rm);

anovaTable_seshComp_postRewLick = movevars(anovaTable_seshComp_postRewLick, "DF", "Before", "SumSq");
writetable(anovaTable_seshComp_postRewLick,[output_fn 'stats\seshComp_postRewLick_two-wayANOVA.xlsx']);
disp('seshComp_postRewLick')
disp(anovaTable_seshComp_postRewLick)

end

%% plot     | statistical window for first-lick avg vs. baseline avg
%sample average
  [sampeak_firstLickActCSP, sampeakIdx_firstLickActCSP] = max(mean(first_lickAlignEventsCSP_acrossTrials(preLickRange,:),2),[],1);
  [sampeak_postCueLickActCSP, sampeakIdx_postCueLickActCSP] = max(mean(postCue_lickAlignEventsCSP_acrossTrials(preLickRange,:),2),[],1);
  [sampeak_postRewLickActCSP, sampeakIdx_postRewLickActCSP] = max(mean(postRew_lickAlignEventsCSP_acrossTrials(preLickRange,:),2),[],1);
  [sampeak_baseCSP, sampeakIdx_baseCSP] = max(mean(CSPevents_acrossTrials(baselineRange,:),2),[],1);
  [sampeak_firstLickActCSM, sampeakIdx_firstLickActCSM] = max(mean(first_lickAlignEventsCSM_acrossTrials(preLickRange,:),2),[],1);
  [sampeak_postCueLickActCSM, sampeakIdx_postCueLickActCSM] = max(mean(postCue_lickAlignEventsCSM_acrossTrials(preLickRange,:),2),[],1);
  [sampeak_postRewLickActCSM, sampeakIdx_postRewLickActCSM] = max(mean(postRew_lickAlignEventsCSM_acrossTrials(preLickRange,:),2),[],1);
  [sampeak_baseCSM, sampeakIdx_baseCSM] = max(mean(CSMevents_acrossTrials(baselineRange,:),2),[],1);
avg_firstLickPane(1,1) = round(mean(sampeakIdx_firstLickActCSP,2,'omitnan')); avg_firstLickPane(1,2) = round(mean(sampeakIdx_firstLickActCSM,2,'omitnan'));

%
tempFig=setFigure; hold on; 
subplot(2,1,1); hold on;
ylim([0 ymax]);  xline(0,'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)'); title('Average Cspk rate preceeding first lick in CS+ trials with statistical windows');
rectangle('Position',[tt(1,36-2+preLickRange(1)+avg_firstLickPane(1,1)-statsWindow)/1000,.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,~gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(firstLickAlignEvents(:,:,~gBlock2),3,'omitnan'),[],2,'omitnan')./sqrt(size(firstLickAlignEvents,2)).*(1000./frameRateHz),'lineProps','k');
subplot(2,1,2); hold on;
ylim([0 ymax]); xline(0,'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)'); title('Average Cspk rate preceeding first lick in CS- trials with statistical windows');
rectangle('Position',[tt(1,36-2+preLickRange(1)+avg_firstLickPane(1,2)-statsWindow)/1000,.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(firstLickAlignEvents(:,:,gBlock2),3,'omitnan'),[],2,'omitnan')./sqrt(size(firstLickAlignEvents,2)).*(1000./frameRateHz),'lineProps','r');
saveas(tempFig, [output_fn 'stats\lickWindow.pdf'],'pdf');
%% ttest    | cell lick-align avg vs. baseline avg
cueAlignEvents = NaN(size(groupAlign_events,1),sum(nsize),sum(animalTrials));
iDend=0; iTrial=0;
for mouse = 1:exptCount
 dends=(1:nsize(1,mouse))+iDend; trials=(1:animalTrials(1,mouse))+iTrial;
  cueAlignEvents(:,dends,trials)=groupAlign.events{1,mouse}(:,1:nsize(1,mouse),1:animalTrials(1,mouse));
 iDend=dends(end); iTrial=trials(end);
end
responseLogical_toLickAll=zeros(1,sum(nsize));
responseLogical_toLickOth=zeros(1,sum(nsize));
responseLogical_toLickCSP=zeros(1,sum(nsize));
responseLogical_toLickCSM=zeros(1,sum(nsize));
responseLogical_toLickNon=zeros(1,sum(nsize));

for c = 1:sum(nsize)
 [response_toLickCSP(1,c),~]=ttest(mean(cueAlignEvents(baselineRange(1)-1+(peakIdx_baseCSP-statsWindow):baselineRange(1)-1+(peakIdx_baseCSP+statsWindow),c,~gBlock2),1,'omitnan'),mean(firstLickAlignEvents(preLickRange(1)-1+(peakIdx_firstLickActCSP-statsWindow):preLickRange(1)-1+(peakIdx_firstLickActCSP+statsWindow),c,~gBlock2),1,'omitnan'));
 [response_toLickCSM(1,c),~]=ttest(mean(cueAlignEvents(baselineRange(1)-1+(peakIdx_baseCSM-statsWindow):baselineRange(1)-1+(peakIdx_baseCSM+statsWindow),c,gBlock2),1,'omitnan'),mean(firstLickAlignEvents(preLickRange(1)-1+(peakIdx_firstLickActCSM-statsWindow):preLickRange(1)-1+(peakIdx_firstLickActCSM+statsWindow),c,gBlock2),1,'omitnan'));
 if response_toLickCSP(1,c)==1 && response_toLickCSM(1,c)==1
  responseLogical_toLickAll(1,c)=1;
 elseif response_toLickCSP(1,c)==0 && response_toLickCSM(1,c)==0
  responseLogical_toLickNon(1,c)=1;
 elseif response_toLickCSP(1,c)==0 && response_toLickCSM(1,c)==1
  responseLogical_toLickCSM(1,c)=1;
 elseif response_toLickCSP(1,c)==1 && response_toLickCSM(1,c)==0
  responseLogical_toLickCSP(1,c)=1;
 else %if one of the responses of that cell is NaN
  responseLogical_toLickOth(1,c)=1;
 end
end

%% plot     | lick-aligned representaion of response groups
 tempFig = setFigure; 
n.lickResp.toAll = sum(responseLogical_toLickAll);
n.lickResp.toCSP = sum(responseLogical_toLickCSP);
n.lickResp.toCSM = sum(responseLogical_toLickCSM);
n.lickResp.toNon = sum(responseLogical_toLickNon);
 pie([n.lickResp.toAll, n.lickResp.toCSP, n.lickResp.toCSM, n.lickResp.toNon]); legend('cells responsive to licks in both CS+/CS- trials','cells only responsive to licks in CS+ trials','cells only responsive to licks in CS- trials','cells without responses to licking')
saveas(tempFig, [output_fn 'stats\pie_propCells_lickAlign_cueSpecificity.pdf'],'pdf');
 
tempFig=setFigure; hold on;
 bv=bar(1,[n.lickResp.toAll/sum(nsize), n.lickResp.toCSP/sum(nsize), n.lickResp.toCSM/sum(nsize), n.lickResp.toNon/sum(nsize)],'stacked');
 legend('cells responsive to licks in both CS+/CS- trials','cells only responsive to licks in CS+ trials','cells only responsive to licks in CS- trials','cells without responses to licking')
 colors = {[153/255 51/255 153/255],[161/255 183/255 13/255],[200/255 78/255 0/255],[255/255 217/255 96/255]}; 
for i = 1:4
     bv(i).FaceColor=colors{i};
end
 saveas(tempFig, [output_fn 'stats\bar_propCells_lickAlign_cueSpecificity.pdf'],'pdf');

 
 tempFig=setFigure; hold on; 
 subplot(2,2,1); hold on;
ylim([0 ymax*2]); xline(0,'k'); xline(.767,'k'); title('Dendrites that respond to licking in CS+/- trials');
shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,logical(responseLogical_toLickAll),~gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(firstLickAlignEvents(:,logical(responseLogical_toLickAll),~gBlock2),3,'omitnan'),[],2,'omitnan')./sqrt(size(firstLickAlignEvents,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,logical(responseLogical_toLickAll),gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(firstLickAlignEvents(:,logical(responseLogical_toLickAll),gBlock2),3,'omitnan'),[],2,'omitnan')./sqrt(size(firstLickAlignEvents,2)).*(1000./frameRateHz),'lineProps','r');
 subplot(2,2,2); hold on;
ylim([0 ymax*2]); xline(0,'k'); xline(.767,'k'); title('Dendrites that respond to licking in CS+ trials');
shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,logical(responseLogical_toLickCSP),~gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(firstLickAlignEvents(:,logical(responseLogical_toLickCSP),~gBlock2),3,'omitnan'),[],2,'omitnan')./sqrt(size(firstLickAlignEvents,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,logical(responseLogical_toLickCSP),gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(firstLickAlignEvents(:,logical(responseLogical_toLickCSP),gBlock2),3,'omitnan'),[],2,'omitnan')./sqrt(size(firstLickAlignEvents,2)).*(1000./frameRateHz),'lineProps','r');
 subplot(2,2,3); hold on;
ylim([0 ymax*2]); xline(0,'k'); xline(.767,'k'); title('Dendrites that respond to licking in CS- trials');
shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,logical(responseLogical_toLickCSM),~gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(firstLickAlignEvents(:,logical(responseLogical_toLickCSM),~gBlock2),3,'omitnan'),[],2,'omitnan')./sqrt(size(firstLickAlignEvents,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,logical(responseLogical_toLickCSM),gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(firstLickAlignEvents(:,logical(responseLogical_toLickCSM),gBlock2),3,'omitnan'),[],2,'omitnan')./sqrt(size(firstLickAlignEvents,2)).*(1000./frameRateHz),'lineProps','r');
 subplot(2,2,4); hold on;
ylim([0 ymax*2]); xline(0,'k'); xline(.767,'k'); title('Dendrites that do not respond to licking');
shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,logical(responseLogical_toLickNon),~gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(firstLickAlignEvents(:,logical(responseLogical_toLickNon),~gBlock2),3,'omitnan'),[],2,'omitnan')./sqrt(size(firstLickAlignEvents,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,logical(responseLogical_toLickNon),gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(firstLickAlignEvents(:,logical(responseLogical_toLickNon),gBlock2),3,'omitnan'),[],2,'omitnan')./sqrt(size(firstLickAlignEvents,2)).*(1000./frameRateHz),'lineProps','r');
saveas(tempFig, [output_fn 'stats\Cspk_lickRelated_cueSpecificity.pdf'],'pdf');
%% plot     | first-lick (response/suppression) 

  %we're gonna set an eval loop here because individual blocks for first-postCue-postRew for CS+&- is vile
eventLabel={'first','postCue','postRew'}; cueLabel={'CSP','CSM'}; 
for cue=1:size(cueLabel,2)
  for event=1:size(eventLabel,2)
      eval(['dendrite.lickAlign.response_' eventLabel{1,event} cueLabel{1,cue} '=[];']);
      eval(['dendrite.lickAlign.suppress_' eventLabel{1,event} cueLabel{1,cue} '=[];']);
    for c = 1:sum(nsize)
%
 eval(['preLickVar = ' eventLabel{1,event} '_lickAlignDFoF' cueLabel{1,cue} '_acrossTrials(1:preLick_frames,c);']);
 eval(['baselineVar = ' cueLabel{1,cue} 'dFoF_acrossTrials(baselineRange,c);']);
 eval(['pos.Diff.lickAlign{c,1} = ([0; diff(' eventLabel{1,event} '_lickAlignDFoF' cueLabel{1,cue} '_acrossTrials(:,c))]);']);
pos.Limit.lickAlign{c,1}.one = (mean(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),1,'omitnan')+(std(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),[],1,'omitnan').*1)); 
pos.Change.lickAlign{c,1}.one = find(pos.Diff.lickAlign{c,1}>pos.Limit.lickAlign{c,1}.one); 
pos.Limit.lickAlign{c,1}.two = (mean(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),1,'omitnan')+(std(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),[],1,'omitnan').*1.5)); 
pos.Change.lickAlign{c,1}.two = find(pos.Diff.lickAlign{c,1}>pos.Limit.lickAlign{c,1}.two); 
pos.Limit.lickAlign{c,1}.three = (mean(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),1,'omitnan')+(std(pos.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)>0])),[],1,'omitnan').*2)); 
pos.Change.lickAlign{c,1}.three = find(pos.Diff.lickAlign{c,1}>pos.Limit.lickAlign{c,1}.three);

eval(['frameIdx.lickAlign.' eventLabel{1,event} cueLabel{1,cue} '_posChangeDiff{1,c} = intersect(preLickRange, pos.Change.lickAlign{c,1}.one);']); 
eval(['frameIdx.lickAlign.' eventLabel{1,event} cueLabel{1,cue} '_posChangeDiff{2,c} = intersect(preLickRange, pos.Change.lickAlign{c,1}.two);']); 
eval(['frameIdx.lickAlign.' eventLabel{1,event} cueLabel{1,cue} '_posChangeDiff{3,c} = intersect(preLickRange, pos.Change.lickAlign{c,1}.three);']);

 eval(['neg.Diff.lickAlign{c,1} = ([0; diff(' eventLabel{1,event} '_lickAlignDFoF' cueLabel{1,cue} '_acrossTrials(:,c))]);']);
neg.Limit.lickAlign{c,1}.one = (mean(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),1,'omitnan')-(std(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),[],1,'omitnan').*1)); 
neg.Change.lickAlign{c,1}.one = find(neg.Diff.lickAlign{c,1}<neg.Limit.lickAlign{c,1}.one); 
neg.Limit.lickAlign{c,1}.two = (mean(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),1,'omitnan')-(std(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),[],1,'omitnan').*1.5)); 
neg.Change.lickAlign{c,1}.two = find(neg.Diff.lickAlign{c,1}<neg.Limit.lickAlign{c,1}.two); 
neg.Limit.lickAlign{c,1}.three = (mean(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),1,'omitnan')-(std(neg.Diff.cueAlign{c,1}(logical([0; diff(baselineVar)<0])),[],1,'omitnan').*2)); 
neg.Change.lickAlign{c,1}.three = find(neg.Diff.lickAlign{c,1}<neg.Limit.lickAlign{c,1}.three);

eval(['frameIdx.lickAlign.' eventLabel{1,event} cueLabel{1,cue} '_negChangeDiff{1,c} = intersect(preLickRange, neg.Change.lickAlign{c,1}.one);']); 
eval(['frameIdx.lickAlign.' eventLabel{1,event} cueLabel{1,cue} '_negChangeDiff{2,c} = intersect(preLickRange, neg.Change.lickAlign{c,1}.two);']); 
eval(['frameIdx.lickAlign.' eventLabel{1,event} cueLabel{1,cue} '_negChangeDiff{3,c} = intersect(preLickRange, neg.Change.lickAlign{c,1}.three);']);

%    
    tempResp_DiffIdx = diff(intersect(preLickRange, find(preLickVar>(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*respReq_lickAlign))))';
    tempResp_FrameIdx = intersect(preLickRange, find(preLickVar>(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*respReq_lickAlign)));
    if length(tempResp_DiffIdx)
      if strfind(tempResp_DiffIdx, [1 1])
        if eval(['frameIdx.lickAlign.' eventLabel{1,event} cueLabel{1,cue} '_posChangeDiff{chngReq_lickAlign,c}'])
           eval(['dendrite.lickAlign.response_' eventLabel{1,event} cueLabel{1,cue} '=[dendrite.lickAlign.response_' eventLabel{1,event} cueLabel{1,cue} ' c];']);
        end
      end
    end
    tempSupp_DiffIdx = diff(intersect(preLickRange, find(preLickVar<(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*suppReq_lickAlign))))';
    tempSupp_FrameIdx = intersect(preLickRange, find(preLickVar<(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*suppReq_lickAlign)));
    if length(tempSupp_DiffIdx)
      if strfind(tempSupp_DiffIdx, [1 1 1 1]) 
        if isempty(strfind(tempResp_DiffIdx, [1 1]))
          if eval(['frameIdx.lickAlign.' eventLabel{1,event} cueLabel{1,cue} '_negChangeDiff{chngReq_lickAlign,c}'])
           eval(['dendrite.lickAlign.suppress_' eventLabel{1,event} cueLabel{1,cue} '=[dendrite.lickAlign.suppress_' eventLabel{1,event} cueLabel{1,cue} ' c];']);
          end
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
eval(['dendriteDiff.lickAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} '= setdiff(dendrite.lickAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} ',dendrite.lickAlign.' sideLabel{2,side} '_'  eventLabel{1,event} cueLabel{1,cue} ');']);
it=it+1;
    end
  end
end
%
titleArray = {'CS+ reponsive firstLick','CS+ suppressive firstLick','CS+ reponsive postCueLick','CS+ suppressive postCueLick','CS+ reponsive postRewLick','CS+ suppressive postRewLick','CS- reponsive firstLick','CS- suppressive firstLick','CS- reponsive postCueLick','CS- suppressive postCueLick','CS- reponsive postRewLick','CS- suppressive postRewLick'};
colorArray = {'k','b','k','b','k','b','r','m','r','m','r','m'};

highlight = [255/255 217/255 96/255]; 
tempFig = setFigure; hold on; c=1;
for cue=1:size(cueLabel,2)
  for event=1:size(eventLabel,2)
    for side=1:size(sideLabel,2)
    eval(['callTrue = ' eventLabel{1,event} '_lickAlignDFoF' cueLabel{1,cue} '_acrossTrials(:,dendriteDiff.lickAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} ');']);
    eval(['callFalse = ' eventLabel{1,event} '_lickAlignDFoF' cueLabel{1,cue} '_acrossTrials(:,setdiff(allcells,dendriteDiff.lickAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} '));']);
    eval(['dendrites = dendriteDiff.lickAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} ';']);
subplot(6,2,c); hold on;
ylim([ymin*3 ymax*3]); xline(0,'k');  title(titleArray{1,c});
rectangle('Position',[(tt(1,baselineRange(1))+10)/1000,ymin*2+.25,((1000./frameRateHz)*response_frames-10)/1000,ymax+5],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:80)/1000,mean(callTrue,2,'omitnan').*(1000./frameRateHz),std(callTrue,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',colorArray{1,c});
% shadedErrorBar_CV(tt(1,36:96)/1000,mean(callFalse,2,'omitnan'),std(callFalse,[],2,'omitnan')./sqrt(size(dendrites,2)),'lineProps',{'Color',[0/255 0/255 201/255]});
    c=c+1;
    end
  end
end
    saveas(tempFig, [output_fn 'stats\dFoF_lickAlign_PSTHs.pdf'],'pdf');

tempFig = setFigure; hold on; c=1;
for cue=1:size(cueLabel,2)
  for event=1:size(eventLabel,2)
    for side=1:size(sideLabel,2)
    eval(['callTrue = ' eventLabel{1,event} '_lickAlignEvents' cueLabel{1,cue} '_acrossTrials(:,dendriteDiff.lickAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} ');']);
    eval(['callFalse = ' eventLabel{1,event} '_lickAlignEvents' cueLabel{1,cue} '_acrossTrials(:,setdiff(allcells,dendriteDiff.lickAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} '));']);
    eval(['dendrites = dendriteDiff.lickAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} ';']);
subplot(6,2,c); hold on;
ylim([0 ymax]); xline(0,'k');  title(titleArray{1,c});
if c<5; rectangle('Position',[tt(1,baselineRange(1))/1000,0+.15,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight); elseif c>4; rectangle('Position',[tt(1,baselineRange(1))/1000,0+.15,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight); end
shadedErrorBar_CV(tt(1,36:80)/1000,mean(callTrue,2,'omitnan').*(1000./frameRateHz),std(callTrue,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',colorArray{1,c});
% shadedErrorBar_CV(tt(1,36:96)/1000,mean(callFalse,2,'omitnan').*(1000./frameRateHz),std(callFalse,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',{'Color',[0/255 0/255 201/255]});
    c=c+1;
    end
  end
end
    saveas(tempFig, [output_fn 'stats\CspkEvents_lickAlign_PSTHs.pdf'],'pdf');

    tempFig = setFigure; hold on; c=1;
for event=1:size(eventLabel,2)
  for cue=1:size(cueLabel,2)
eval(['respN = size(' eventLabel{1,event} '_lickAlignEvents' cueLabel{1,cue} '_acrossTrials(:,dendriteDiff.lickAlign.response_' eventLabel{1,event} cueLabel{1,cue} '),2);']);
eval(['suppN = size(' eventLabel{1,event} '_lickAlignEvents' cueLabel{1,cue} '_acrossTrials(:,dendriteDiff.lickAlign.suppress_' eventLabel{1,event} cueLabel{1,cue} '),2);']);
       elseN = sum(nsize)-(respN+suppN);

subplot(3,2,c); 
pie([respN, suppN, elseN]);
c=c+1; title([eventLabel{1,event} '-' cueLabel{1,cue}]);
if c == (size(eventLabel,2)+size(cueLabel,2))
legend('responsive cells','suppressed cells','neutral cells');   
end
  end
end
    saveas(tempFig, [output_fn 'stats\pie_propCells_lickAlign.pdf'],'pdf');

    
tempFig=setFigure; hold on;
subplot(3,1,1); hold on;
shadedErrorBar_CV(tt(1,36:80)/1000,mean(first_lickAlignEventsCSP_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(first_lickAlignEventsCSP_acrossTrials,[],2,'omitnan')./sqrt(size(first_lickAlignEventsCSP_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:80)/1000,mean(first_lickAlignEventsCSM_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(first_lickAlignEventsCSM_acrossTrials,[],2,'omitnan')./sqrt(size(first_lickAlignEventsCSM_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
ylim([0 ymax]); xline(0,'k');  title('first lick');
subplot(3,1,2); hold on;
shadedErrorBar_CV(tt(1,36:80)/1000,mean(postCue_lickAlignEventsCSP_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(postCue_lickAlignEventsCSP_acrossTrials,[],2,'omitnan')./sqrt(size(postCue_lickAlignEventsCSP_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:80)/1000,mean(postCue_lickAlignEventsCSM_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(postCue_lickAlignEventsCSM_acrossTrials,[],2,'omitnan')./sqrt(size(postCue_lickAlignEventsCSM_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
ylim([0 ymax]); xline(0,'k');  title('first postCue lick');
subplot(3,1,3); hold on;
shadedErrorBar_CV(tt(1,36:80)/1000,mean(postRew_lickAlignEventsCSP_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(postRew_lickAlignEventsCSP_acrossTrials,[],2,'omitnan')./sqrt(size(postRew_lickAlignEventsCSP_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:80)/1000,mean(postRew_lickAlignEventsCSM_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(postRew_lickAlignEventsCSM_acrossTrials,[],2,'omitnan')./sqrt(size(postRew_lickAlignEventsCSM_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
ylim([0 ymax+1]); xline(0,'k');  title('first postRew lick');
saveas(tempFig, [output_fn 'stats\allEvents_firstLick.pdf'],'pdf');

tempFig=setFigure; hold on; 
 ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('Average Cspk for CS+/- trials with first lick range');
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSPevents_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSPevents_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSPevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:96)/1000,mean(CSMevents_acrossTrials(36:96,:),2,'omitnan').*(1000./frameRateHz),std(CSMevents_acrossTrials(36:96,:),[],2,'omitnan')./sqrt(size(CSMevents_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
errorbar((mean(firstLickStart(1,~gBlock2),2,'omitnan')-prewin_frames).*(1000./frameRateHz)/1000,.25,std((firstLickStart(1,~gBlock2)),[],2,'omitnan')./sqrt(size(~gBlock2,2)).*(1000./frameRateHz)/1000,'horizontal','k');
errorbar((mean(firstLickStart(1,gBlock2),2,'omitnan')-prewin_frames).*(1000./frameRateHz)/1000,.25,std((firstLickStart(1,gBlock2)),[],2,'omitnan')./sqrt(size(gBlock2,2)).*(1000./frameRateHz)/1000,'horizontal','r');
saveas(tempFig, [output_fn 'stats\Cspk_withFirstLick.pdf'],'pdf');   
%% plot     | first-lick and first-lick by triads
%lickFrameRange=prewin_frames-postLick_frames+2:prewin_frames+postLick_frames+postLick_frames+1;
lickFrameRange=prewin_frames-postLick_frames+1:prewin_frames+postLick_frames+postLick_frames;

 tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
    subplot(2,1,1); hold on;
    shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstLickAlignEvents(:,:,~gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,~gBlock2),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstLickAlignEvents(:,:,gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,gBlock2),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','r');
    hold on;
    title('First lick');
    xlabel('Time from lick (s)');
    ylabel('Cspk rate (Hz)');
    ylim([0 ymax]); vline(0,'k');
    
    subplot(2,1,2);hold on;
    shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstLickAlign(:,gBlock2),3,'omitnan'),2,'omitnan'), (std(mean(firstLickAlign(:,gBlock2),3,'omitnan'),[],2,'omitnan'))./sqrt((sum(~isnan(mean(firstLickAlign,1,'omitnan'))))),'lineProps','r');
    shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstLickAlign(:,~gBlock2),3,'omitnan'),2,'omitnan'), (std(mean(firstLickAlign(:,~gBlock2),3,'omitnan'),[],2,'omitnan'))./sqrt((sum(~isnan(mean(firstLickAlign,1,'omitnan'))))),'lineProps','k');
    xlabel('Time from lick (s)');
    ylabel('Lick rate (Hz)');
    ylim([0 1]);
    title([num2str(sum(~isnan(firstLickStart(:,~gBlock2)))) '+ ' num2str(sum(~isnan(firstLickStart(:,gBlock2)))) '- trials with lick']);
    
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Cspk relative to lick']);
    saveas(tempFig, [output_fn 'stats\firstLick.pdf'],'pdf'); 
    
    % find splits for each mouse
    splitBy = 3;
    for mus = 1:size(expt,2)
        b_temp = find(block2{1,mus}); b_split = floor(size(b_temp,2)./splitBy); b_range = 1:b_split:size(b_temp,2); 
        r_temp = find(~block2{1,mus}); r_split = floor(size(r_temp,2)./splitBy); r_range = 1:r_split:size(r_temp,2); 
    for i = 1:splitBy
        rewSplit{mus,i} = r_temp(1,r_range(1,i):(r_range(1,i)-1+r_split));
        b2Split{mus,i} = b_temp(1,b_range(1,i):(b_range(1,i)-1+b_split));
    end
    end
    % quick conversion to turn mus-specific trial no. into group trial no.
    gAT=0;
    for i = 1:size(animalTrials,2)
        gAT = [gAT gAT(1,i)+animalTrials(1,i)];
    end %the final value should equal the sum of all trials - as a check 
    for i = 1:splitBy
        gRewSplit{1,i}=[]; gB2Split{1,i}=[]; 
        for mus = 1:size(expt,2)
        gRewSplit{1,i} = [gRewSplit{1,i} (rewSplit{mus,i}+(gAT(1,mus)))];
        gB2Split{1,i} = [gB2Split{1,i} (b2Split{mus,i}+(gAT(1,mus)))];
        end
    end
    %plot
    tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
    for i = 1:splitBy
        subplot(2,splitBy,i); hold on;
    shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstLickAlignEvents(:,:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstLickAlignEvents(:,:,gB2Split{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','r');
    hold on;
    xlabel('Time from lick (s)');
    ylabel('Cspk rate (Hz)');
    ylim([0 ymax]); vline(0,'k');
        subplot(2,splitBy,i+splitBy); hold on;
    shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(firstLickAlign(:,gRewSplit{1,i}),2,'omitnan').*(1000./frameRateHz), (std(firstLickAlign(:,gRewSplit{1,i}),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(firstLickAlign(:,gB2Split{1,i}),2,'omitnan').*(1000./frameRateHz), (std(firstLickAlign(:,gB2Split{1,i}),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','r');
    hold on;
    xlabel('Time from lick (s)');
    ylabel('Lick rate (Hz)');
    ylim([0 35]); vline(0,'k');
    end
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Cspk relative to lick']);
    saveas(tempFig, [output_fn 'stats\firstLick_splitsBy' num2str(splitBy) '.pdf'],'pdf'); 
    
    %plot cue-aligned
    cueAlignEvents = NaN(150,sum(nsize),sum(animalTrials)); nN = 1; nT = 1;
    for i = 1:size(expt,2)
     cueAlignEvents(:,nN:nN+nsize(1,i)-1,nT:nT+animalTrials(1,i)-1) = groupAlign.events{1,i};
        nN = nN+nsize(1,i); nT = nT+animalTrials(1,i);
    end
    
    tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
    for i = 1:splitBy
        subplot(splitBy,1,i); hold on;
    shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents(36:96,:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(cueAlignEvents(36:96,:,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
    shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents(36:96,:,gB2Split{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(cueAlignEvents(36:96,:,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','r');
    hold on;
    xlabel('Time from cue (s)');
    ylabel('Cspk rate (Hz)');
    ylim([0 ymax]);  vline(0,'k'); vline(.767,'k');
    end
    sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Cspk relative to cue']);
    saveas(tempFig, [output_fn 'stats\cueAlign_splits' num2str(splitBy) '.pdf'],'pdf'); 
    
    %% embedding a section that will perform basic stats on the split
         for i = 1:splitBy
          for mus = 1:size(expt,2)
         n.lickAlignSplit.CSP_first_trialAvg_count(mus,i) = sum(~isnan(firstLickStart(:,gRewSplit{1,i}(1,gAT(mus)<gRewSplit{1,i}&gRewSplit{1,i}<gAT(mus+1))))); 
         n.lickAlignSplit.CSM_first_trialAvg_count(mus,i) = sum(~isnan(firstLickStart(:,gB2Split{1,i}(1,gAT(mus)<gB2Split{1,i}&gB2Split{1,i}<gAT(mus+1)))));
         n.lickAlignSplit.CSP_firstPostCue_trialAvg_count(mus,i) = sum(~isnan(firstPreRew_lickStart(:,gRewSplit{1,i}(1,gAT(mus)<gRewSplit{1,i}&gRewSplit{1,i}<gAT(mus+1)))));
         n.lickAlignSplit.CSM_firstPostCue_trialAvg_count(mus,i) = sum(~isnan(firstPreRew_lickStart(:,gB2Split{1,i}(1,gAT(mus)<gB2Split{1,i}&gB2Split{1,i}<gAT(mus+1)))));
         n.lickAlignSplit.CSP_firstPostRew_trialAvg_count(mus,i) = sum(~isnan(firstPostRew_lickStart(:,gRewSplit{1,i}(1,gAT(mus)<gRewSplit{1,i}&gRewSplit{1,i}<gAT(mus+1)))));
         n.lickAlignSplit.CSM_firstPostRew_trialAvg_count(mus,i) = sum(~isnan(firstPostRew_lickStart(:,gB2Split{1,i}(1,gAT(mus)<gB2Split{1,i}&gB2Split{1,i}<gAT(mus+1)))));
         n.cueAlignSplit.CSP_count(mus,i) = length(gRewSplit{1,i}(1,gAT(mus)<gRewSplit{1,i}&gRewSplit{1,i}<gAT(mus+1))); 
         n.cueAlignSplit.CSM_count(mus,i) = length(gB2Split{1,i}(1,gAT(mus)<gB2Split{1,i}&gB2Split{1,i}<gAT(mus+1)));
          end
         [peak_lickAlignSplitCSP(1,i), peakIdx_lickAlignSplitCSP(1,i)] = max(mean(mean(firstLickAlignEvents(preLickRange,:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan'),[],1);
         [peak_lickAlignSplitCSM(1,i), peakIdx_lickAlignSplitCSM(1,i)] = max(mean(mean(firstLickAlignEvents(preLickRange,:,gB2Split{1,i}),3,'omitnan'),2,'omitnan'),[],1);
         [peak_lickAlignSplit_postCueCSP(1,i), peakIdx_lickAlignSplit_postCueCSP(1,i)] = max(mean(mean(firstPreRew_lickAlignEvents(preLickRange,:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan'),[],1);
         [peak_lickAlignSplit_postCueCSM(1,i), peakIdx_lickAlignSplit_postCueCSM(1,i)] = max(mean(mean(firstPreRew_lickAlignEvents(preLickRange,:,gB2Split{1,i}),3,'omitnan'),2,'omitnan'),[],1);
         [peak_lickAlignSplit_postRewCSP(1,i), peakIdx_lickAlignSplit_postRewCSP(1,i)] = max(mean(mean(firstPostRew_lickAlignEvents(preLickRange,:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan'),[],1);
         [peak_lickAlignSplit_postRewCSM(1,i), peakIdx_lickAlignSplit_postRewCSM(1,i)] = max(mean(mean(firstPostRew_lickAlignEvents(preLickRange,:,gB2Split{1,i}),3,'omitnan'),2,'omitnan'),[],1);
         [peak_baselineSplitCSP(1,i), peakIdx_baselineSplitCSP(1,i)] = max(mean(mean(cueAlignEvents(baselineRange,:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan'),[],1);
         [peak_baselineSplitCSM(1,i), peakIdx_baselineSplitCSM(1,i)] = max(mean(mean(cueAlignEvents(baselineRange,:,gB2Split{1,i}),3,'omitnan'),2,'omitnan'),[],1);

         m.lickAlignSplit.CSP_trialAvg_Cspk(1,i) = mean(mean(mean(firstLickAlignEvents((preLickRange(peakIdx_lickAlignSplitCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplitCSP(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz));
         s.lickAlignSplit.CSP_trialAvg_Cspk(1,i) = mean(std(mean(firstLickAlignEvents((preLickRange(peakIdx_lickAlignSplitCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplitCSP(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz));
         m.lickAlignSplit.CSM_trialAvg_Cspk(1,i) = mean(mean(mean(firstLickAlignEvents((preLickRange(peakIdx_lickAlignSplitCSM(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplitCSM(1,i))+statsWindow),:,gB2Split{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz));
         s.lickAlignSplit.CSM_trialAvg_Cspk(1,i) = mean(std(mean(firstLickAlignEvents((preLickRange(peakIdx_lickAlignSplitCSM(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplitCSM(1,i))+statsWindow),:,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz));
         m.lickAlignSplit_postCue.CSP_trialAvg_Cspk(1,i) = mean(mean(mean(firstPreRew_lickAlignEvents((preLickRange(peakIdx_lickAlignSplit_postCueCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplit_postCueCSP(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz));
         s.lickAlignSplit_postCue.CSP_trialAvg_Cspk(1,i) = mean(std(mean(firstPreRew_lickAlignEvents((preLickRange(peakIdx_lickAlignSplit_postCueCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplit_postCueCSP(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz));
         m.lickAlignSplit_postCue.CSM_trialAvg_Cspk(1,i) = mean(mean(mean(firstPreRew_lickAlignEvents((preLickRange(peakIdx_lickAlignSplit_postCueCSM(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplit_postCueCSM(1,i))+statsWindow),:,gB2Split{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz));
         s.lickAlignSplit_postCue.CSM_trialAvg_Cspk(1,i) = mean(std(mean(firstPreRew_lickAlignEvents((preLickRange(peakIdx_lickAlignSplit_postCueCSM(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplit_postCueCSM(1,i))+statsWindow),:,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz));
         m.lickAlignSplit_postRew.CSP_trialAvg_Cspk(1,i) = mean(mean(mean(firstPostRew_lickAlignEvents((preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz));
         s.lickAlignSplit_postRew.CSP_trialAvg_Cspk(1,i) = mean(std(mean(firstPostRew_lickAlignEvents((preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz));
         m.lickAlignSplit_postRew.CSM_trialAvg_Cspk(1,i) = mean(mean(mean(firstPostRew_lickAlignEvents((preLickRange(peakIdx_lickAlignSplit_postRewCSM(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplit_postRewCSM(1,i))+statsWindow),:,gB2Split{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz));
         s.lickAlignSplit_postRew.CSM_trialAvg_Cspk(1,i) = mean(std(mean(firstPostRew_lickAlignEvents((preLickRange(peakIdx_lickAlignSplit_postRewCSM(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplit_postRewCSM(1,i))+statsWindow),:,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz));
         
         
         [peak_cueAlignSplit_postCueCSP(1,i), peakIdx_cueAlignSplit_postCueCSP(1,i)] = max(mean(mean(cueAlignEvents(postCueRange,:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan'),[],1);
         [peak_cueAlignSplit_postCueCSM(1,i), peakIdx_cueAlignSplit_postCueCSM(1,i)] = max(mean(mean(cueAlignEvents(postCueRange,:,gB2Split{1,i}),3,'omitnan'),2,'omitnan'),[],1);
         [peak_cueAlignSplit_postRewCSP(1,i), peakIdx_cueAlignSplit_postRewCSP(1,i)] = max(mean(mean(cueAlignEvents(postRewRange,:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan'),[],1);
         [peak_cueAlignSplit_postRewCSM(1,i), peakIdx_cueAlignSplit_postRewCSM(1,i)] = max(mean(mean(cueAlignEvents(postRewRange,:,gB2Split{1,i}),3,'omitnan'),2,'omitnan'),[],1);
         
         % if seshType == 4
         % %aligned to well-timed peak off CS-
         % m.cueAlignSplit.CSP_postCue_trialAvg_Cspk(1,i) = mean(mean(mean(cueAlignEvents((postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz));
         % s.cueAlignSplit.CSP_postCue_trialAvg_Cspk(1,i) = mean(std(mean(cueAlignEvents((postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz));
         % else
         m.cueAlignSplit.CSP_postCue_trialAvg_Cspk(1,i) = mean(mean(mean(cueAlignEvents((postCueRange(peakIdx_cueAlignSplit_postCueCSP(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSP(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz));
         s.cueAlignSplit.CSP_postCue_trialAvg_Cspk(1,i) = mean(std(mean(cueAlignEvents((postCueRange(peakIdx_cueAlignSplit_postCueCSP(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSP(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz));    
         % end
         m.cueAlignSplit.CSP_postRew_trialAvg_Cspk(1,i) = mean(mean(mean(cueAlignEvents((postRewRange(peakIdx_cueAlignSplit_postRewCSP(1,i))-statsWindow):(postRewRange(peakIdx_cueAlignSplit_postRewCSP(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz));
         s.cueAlignSplit.CSP_postRew_trialAvg_Cspk(1,i) = mean(std(mean(cueAlignEvents((postRewRange(peakIdx_cueAlignSplit_postRewCSP(1,i))-statsWindow):(postRewRange(peakIdx_cueAlignSplit_postRewCSP(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz));
         m.cueAlignSplit.CSP_baseline_trialAvg_Cspk(1,i) = mean(mean(mean(cueAlignEvents((baselineRange(peakIdx_baselineSplitCSP(1,i))-statsWindow):(baselineRange(peakIdx_baselineSplitCSP(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz));
         s.cueAlignSplit.CSP_baseline_trialAvg_Cspk(1,i) = mean(std(mean(cueAlignEvents((baselineRange(peakIdx_baselineSplitCSP(1,i))-statsWindow):(baselineRange(peakIdx_baselineSplitCSP(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz));
         
         m.cueAlignSplit.CSM_postCue_trialAvg_Cspk(1,i) = mean(mean(mean(cueAlignEvents((postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))+statsWindow),:,gB2Split{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz));
         s.cueAlignSplit.CSM_postCue_trialAvg_Cspk(1,i) = mean(std(mean(cueAlignEvents((postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))+statsWindow),:,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz));
         m.cueAlignSplit.CSM_postRew_trialAvg_Cspk(1,i) = mean(mean(mean(cueAlignEvents((postRewRange(peakIdx_cueAlignSplit_postRewCSM(1,i))-statsWindow):(postRewRange(peakIdx_cueAlignSplit_postRewCSM(1,i))+statsWindow),:,gB2Split{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz));
         s.cueAlignSplit.CSM_postRew_trialAvg_Cspk(1,i) = mean(std(mean(cueAlignEvents((postRewRange(peakIdx_cueAlignSplit_postRewCSM(1,i))-statsWindow):(postRewRange(peakIdx_cueAlignSplit_postRewCSM(1,i))+statsWindow),:,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz));
         m.cueAlignSplit.CSM_baseline_trialAvg_Cspk(1,i) = mean(mean(mean(cueAlignEvents((baselineRange(peakIdx_baselineSplitCSM(1,i))-statsWindow):(baselineRange(peakIdx_baselineSplitCSM(1,i))+statsWindow),:,gB2Split{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz));
         s.cueAlignSplit.CSM_baseline_trialAvg_Cspk(1,i) = mean(std(mean(cueAlignEvents((baselineRange(peakIdx_baselineSplitCSM(1,i))-statsWindow):(baselineRange(peakIdx_baselineSplitCSM(1,i))+statsWindow),:,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz));
         
         end
        
        tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
        for i = 1:splitBy
        subplot(1,splitBy,i); hold on;
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstLickAlignEvents(:,:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstLickAlignEvents(:,:,gB2Split{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','r');
        hold on;
        xlabel('Time from lick (s)');
        ylabel('Cspk rate (Hz)');
        ylim([0 ymax]); vline(0,'k');
        subtitle(['phase ' num2str(i) ': ' num2str(sum(n.lickAlignSplit.CSP_first_trialAvg_count(:,i))) ' CS+| ' num2str(sum(n.lickAlignSplit.CSM_first_trialAvg_count(:,i))) ' CS-'])
        end
        sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Cspk relative to first lick']);
        saveas(tempFig, [output_fn 'stats\firstLickResponseOnly_splitsBy' num2str(splitBy) '.pdf'],'pdf'); 
        
        tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
        for i = 1:splitBy
        subplot(1,splitBy,i); hold on;
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstPreRew_lickAlignEvents(:,:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstPreRew_lickAlignEvents(:,:,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstPreRew_lickAlignEvents(:,:,gB2Split{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstPreRew_lickAlignEvents(:,:,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','r');
        hold on;
        xlabel('Time from lick (s)');
        ylabel('Cspk rate (Hz)');
        ylim([0 ymax]); vline(0,'k');
        subtitle(['phase ' num2str(splitBy) ': ' num2str(sum(n.lickAlignSplit.CSP_firstPostCue_trialAvg_count(:,i))) ' CS+| ' num2str(sum(n.lickAlignSplit.CSM_firstPostCue_trialAvg_count(:,i))) ' CS-'])
        end
        sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Post cue Cspk relative to lick']);
        saveas(tempFig, [output_fn 'stats\firstLickPostCueResponseOnly_splitsBy' num2str(splitBy) '.pdf'],'pdf'); 
        
        tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
        for i = 1:splitBy
        subplot(1,splitBy,i); hold on;
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstPostRew_lickAlignEvents(:,:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstPostRew_lickAlignEvents(:,:,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstPostRew_lickAlignEvents(:,:,gB2Split{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstPostRew_lickAlignEvents(:,:,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','r');
        hold on;
        xlabel('Time from lick (s)');
        ylabel('Cspk rate (Hz)');
        ylim([0 ymax]); vline(0,'k');
        subtitle(['phase ' num2str(splitBy) ': ' num2str(sum(n.lickAlignSplit.CSP_firstPostRew_trialAvg_count(:,i))) ' CS+| ' num2str(sum(n.lickAlignSplit.CSM_firstPostRew_trialAvg_count(:,i))) ' CS-'])
        end
        sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Post reward Cspk relative to lick']);
        saveas(tempFig, [output_fn 'stats\firstLickPostRewardResponseOnly_splitsBy' num2str(splitBy) '.pdf'],'pdf'); 
        
        tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
        for i = 1:splitBy
        subplot(1,splitBy,i); hold on;
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents(36:96,:,gRewSplit{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(cueAlignEvents(36:96,:,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','k');
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents(36:96,:,gB2Split{1,i}),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(cueAlignEvents(36:96,:,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC),'lineProps','r');
        hold on;
        xlabel('Time from lick (s)');
        ylabel('Cspk rate (Hz)');
        ylim([0 ymax]); vline(0,'k'); vline(.767,'k'); xlim([tt(1,36)/1000 tt(1,96)/1000]);
        subtitle(['phase ' num2str(i) ': ' num2str(sum(n.cueAlignSplit.CSP_count(:,i))) ' CS+| ' num2str(sum(n.cueAlignSplit.CSM_count(:,i))) ' CS-'])
        end
        sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Cspk relative to first lick']);
        saveas(tempFig, [output_fn 'stats\cueAlignResponseOnly_splitsBy' num2str(splitBy) '.pdf'],'pdf'); 
        
%
        tempFig=setFigure; hold on; 
        shadedErrorBar_CV(1:splitBy, m.lickAlignSplit.CSP_trialAvg_Cspk, s.lickAlignSplit.CSP_trialAvg_Cspk, 'lineProps', 'k');
        shadedErrorBar_CV(1:splitBy, m.lickAlignSplit.CSM_trialAvg_Cspk, s.lickAlignSplit.CSM_trialAvg_Cspk, 'lineProps', 'r');
        xlabel('Phase');
        ylabel('Cspk rate (Hz)');
        sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | lick aligned peak Cspk rate by phase']);
        ylim([0 ymax+1]); xticks([1 2 3]);
        saveas(tempFig, [output_fn 'stats\lickAlignPeak_splits' num2str(splitBy) '.pdf'],'pdf'); 
     
        tempFig=setFigure; hold on; 
        plot(1:splitBy, m.lickAlignSplit.CSP_trialAvg_Cspk, 'k', 'linewidth', 2);
        errorbar(1:splitBy, m.lickAlignSplit.CSP_trialAvg_Cspk, s.lickAlignSplit.CSP_trialAvg_Cspk, 'k', 'linewidth', 1.5);
        plot(1:splitBy, m.lickAlignSplit.CSM_trialAvg_Cspk, 'r', 'linewidth', 2);
        errorbar(1:splitBy, m.lickAlignSplit.CSM_trialAvg_Cspk, s.lickAlignSplit.CSM_trialAvg_Cspk, 'r', 'linewidth', 1.5);
        xlabel('Phase');
        ylabel('Cspk rate (Hz)');
        sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | lick aligned peak Cspk rate by phase']);
        ylim([0 ymax+1]); xticks([1 2 3]);
        saveas(tempFig, [output_fn 'stats\lickAlignPeak_splits' num2str(splitBy) '_lineGraph.pdf'],'pdf'); 
     
        tempFig=setFigure; hold on; 
        plot(1:splitBy, m.lickAlignSplit_postCue.CSP_trialAvg_Cspk, 'k', 'linewidth', 2);
        errorbar(1:splitBy, m.lickAlignSplit_postCue.CSP_trialAvg_Cspk, s.lickAlignSplit_postCue.CSP_trialAvg_Cspk, 'k', 'linewidth', 1.5);
        plot(1:splitBy, m.lickAlignSplit_postCue.CSM_trialAvg_Cspk, 'r', 'linewidth', 2);
        errorbar(1:splitBy, m.lickAlignSplit_postCue.CSM_trialAvg_Cspk, s.lickAlignSplit_postCue.CSM_trialAvg_Cspk, 'r', 'linewidth', 1.5);
        xlabel('Phase');
        ylabel('Cspk rate (Hz)');
        sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | lick aligned peak Cspk rate after cue by phase']);
        ylim([0 ymax+1]); xticks([1 2 3]);
        saveas(tempFig, [output_fn 'stats\postCue_lickAlignPeak_splits' num2str(splitBy)  '_lineGraph.pdf'],'pdf'); 
     
        tempFig=setFigure; hold on; 
        plot(1:splitBy, m.lickAlignSplit_postRew.CSP_trialAvg_Cspk, 'k', 'linewidth', 2);
        errorbar(1:splitBy, m.lickAlignSplit_postRew.CSP_trialAvg_Cspk, s.lickAlignSplit_postRew.CSP_trialAvg_Cspk, 'k', 'linewidth', 1.5);
        plot(1:splitBy, m.lickAlignSplit_postRew.CSM_trialAvg_Cspk, 'r', 'linewidth', 2);
        errorbar(1:splitBy, m.lickAlignSplit_postRew.CSM_trialAvg_Cspk, s.lickAlignSplit_postRew.CSM_trialAvg_Cspk, 'r', 'linewidth', 1.5);
        xlabel('Phase');
        ylabel('Cspk rate (Hz)');
        sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | lick aligned peak Cspk rate after reward by phase']);
        ylim([0 ymax+1]); xticks([1 2 3]);
        saveas(tempFig, [output_fn 'stats\postReward_lickAlignPeak_splits' num2str(splitBy)  '_lineGraph.pdf'],'pdf'); 
 
         tempFig=setFigure; hold on; 
        plot(1:splitBy, m.cueAlignSplit.CSP_postCue_trialAvg_Cspk, 'k', 'linewidth', 2);
        errorbar(1:splitBy, m.cueAlignSplit.CSP_postCue_trialAvg_Cspk, s.cueAlignSplit.CSP_postCue_trialAvg_Cspk, 'k', 'linewidth', 1.5);
        plot(1:splitBy, m.cueAlignSplit.CSM_postCue_trialAvg_Cspk, 'r', 'linewidth', 2);
        errorbar(1:splitBy, m.cueAlignSplit.CSM_postCue_trialAvg_Cspk, s.cueAlignSplit.CSM_postCue_trialAvg_Cspk, 'r', 'linewidth', 1.5);
        xlabel('Phase');
        ylabel('Cspk rate (Hz)');
        sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | cue aligned peak Cspk rate after cue by phase']);
        ylim([0 ymax+1]); xticks([1 2 3]);
        saveas(tempFig, [output_fn 'stats\postCue_cueAlignPeak_splits' num2str(splitBy)  '_lineGraph.pdf'],'pdf'); 
 
         tempFig=setFigure; hold on; 
        plot(1:splitBy, m.cueAlignSplit.CSP_postRew_trialAvg_Cspk, 'k', 'linewidth', 2);
        errorbar(1:splitBy, m.cueAlignSplit.CSP_postRew_trialAvg_Cspk, s.cueAlignSplit.CSP_postRew_trialAvg_Cspk, 'k', 'linewidth', 1.5);
        plot(1:splitBy, m.cueAlignSplit.CSM_postRew_trialAvg_Cspk, 'r', 'linewidth', 2);
        errorbar(1:splitBy, m.cueAlignSplit.CSM_postRew_trialAvg_Cspk, s.cueAlignSplit.CSM_postRew_trialAvg_Cspk, 'r', 'linewidth', 1.5);
        xlabel('Phase');
        ylabel('Cspk rate (Hz)');
        sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | cue aligned peak Cspk rate after reward by phase']);
        ylim([0 ymax+1]); xticks([1 2 3]);
        saveas(tempFig, [output_fn 'stats\postReward_cueAlignPeak_splits' num2str(splitBy)  '_lineGraph.pdf'],'pdf'); 
 
% want to also plot the rate as a probability metric - #spikes/total trials
        for c = 1:size(firstLickAlignEvents,2)
            for f = 1:size(firstLickAlignEvents,1)
                firstLick_spikeProbability_CSP(f,c) = sum(firstLickAlignEvents(f,c,~gBlock2)==1)./sum(~isnan(firstLickAlignEvents(f,c,~gBlock2)));
                firstLick_spikeProbability_CSM(f,c) = sum(firstLickAlignEvents(f,c,gBlock2)==1)./sum(~isnan(firstLickAlignEvents(f,c,gBlock2)));
            end
        end

    % Here, I am plotting the Cspk probability of cells aligned to the first lick
        tempFig=setFigure; hold on;
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(firstLick_spikeProbability_CSP,2), std(firstLick_spikeProbability_CSP,[],2)./sqrt(nIC), 'lineProps', 'k')
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(firstLick_spikeProbability_CSM,2), std(firstLick_spikeProbability_CSM,[],2)./sqrt(nIC), 'lineProps', 'r')
        xlim([0 1]); ylim([0 .15]); vline(0,'k'); 
        xlabel('time from lick (s)');
        ylabel('Cspk probability (%)');
        sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | lick aligned peak Cspk probability']);
        saveas(tempFig, [output_fn 'stats\lickAlign_CspkProbability.pdf'],'pdf'); 

%% ANOVA    | for Split
if seshType == 4
% Now I want a probability histogram for Cspk probability at peak
    probabilityThreshold = 0.05;
    if find(max(mean(firstLick_spikeProbability_CSP,2,'omitnan')))>1
        tempFig=setFigure; hold on;
        histogram(firstLick_spikeProbability_CSP(find(mean(firstLick_spikeProbability_CSP,2,'omitnan') == max(mean(firstLick_spikeProbability_CSP,2,'omitnan')))-statsWindow:find(mean(firstLick_spikeProbability_CSP,2,'omitnan') == max(mean(firstLick_spikeProbability_CSP,2,'omitnan')))+statsWindow,:), 'BinWidth',0.05, 'Normalization','percentage', 'FaceColor','k' ,'EdgeColor','none');
        histogram(firstLick_spikeProbability_CSM(find(mean(firstLick_spikeProbability_CSM,2,'omitnan') == max(mean(firstLick_spikeProbability_CSM,2,'omitnan')))-statsWindow:find(mean(firstLick_spikeProbability_CSM,2,'omitnan') == max(mean(firstLick_spikeProbability_CSM,2,'omitnan')))+statsWindow,:), 'BinWidth',0.05, 'Normalization','percentage', 'FaceColor','r' ,'EdgeColor','none');
        xlabel('Cspk probability (%)');
        ylabel('Proportion of cells (%)');
        sgtitle(['lick aligned Cspk probability histogram at peak +/-' num2str(statsWindow) ': ' num2str(find(mean(firstLick_spikeProbability_CSP,2,'omitnan') == max(mean(firstLick_spikeProbability_CSP,2,'omitnan')))) ' CSP| ' num2str(find(mean(firstLick_spikeProbability_CSM,2,'omitnan') == max(mean(firstLick_spikeProbability_CSM,2,'omitnan')))) ' CSM']);
        saveas(tempFig, [output_fn 'stats\lickAligned_CspkProbabilityHistogram.pdf'],'pdf'); 
    
    % Normalize to log base 2 since spike probability is binary
        tempFig=setFigure; hold on;
        histogram(log2(firstLick_spikeProbability_CSP(mean(firstLick_spikeProbability_CSP,2,'omitnan') == max(mean(firstLick_spikeProbability_CSP,2,'omitnan')),:)), 'BinWidth',0.5, 'Normalization','percentage', 'FaceColor','k' ,'EdgeColor','none');
        histogram(log2(firstLick_spikeProbability_CSM(mean(firstLick_spikeProbability_CSM,2,'omitnan') == max(mean(firstLick_spikeProbability_CSM,2,'omitnan')),:)), 'BinWidth',0.5, 'Normalization','percentage', 'FaceColor','r' ,'EdgeColor','none');
        xlabel('Log 2 Cspk probability (%)');
        ylabel('Proportion of cells (%)');
        sgtitle(['lick aligned Cspk log(probability) histogram at peak: ' num2str(find(mean(firstLick_spikeProbability_CSP,2,'omitnan') == max(mean(firstLick_spikeProbability_CSP,2,'omitnan')))) ' CSP| ' num2str(find(mean(firstLick_spikeProbability_CSM,2,'omitnan') == max(mean(firstLick_spikeProbability_CSM,2,'omitnan')))) ' CSM']);
        saveas(tempFig, [output_fn 'stats\log_lickAligned_CspkProbabilityHistogram.pdf'],'pdf'); 
    
    % Plot Cspk rate of cells with Cspk probability > 0.05
        firstLick_cSpkProbThreshold_CSP = mean(firstLick_spikeProbability_CSP(find(mean(firstLick_spikeProbability_CSP,2,'omitnan') == max(mean(firstLick_spikeProbability_CSP,2,'omitnan')))-statsWindow:find(mean(firstLick_spikeProbability_CSP,2,'omitnan') == max(mean(firstLick_spikeProbability_CSP,2,'omitnan')))+statsWindow,:) > probabilityThreshold, 1,'omitnan') > 0;
        firstLick_cSpkProbThreshold_CSM = mean(firstLick_spikeProbability_CSM(find(mean(firstLick_spikeProbability_CSM,2,'omitnan') == max(mean(firstLick_spikeProbability_CSM,2,'omitnan')))-statsWindow:find(mean(firstLick_spikeProbability_CSM,2,'omitnan') == max(mean(firstLick_spikeProbability_CSM,2,'omitnan')))+statsWindow,:) > probabilityThreshold, 1,'omitnan') > 0;
        tempFig=setFigure; hold on;
        subplot(2,1,1); hold on;
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstLickAlignEvents(:,firstLick_cSpkProbThreshold_CSP,~gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), std(mean(firstLickAlignEvents(:,firstLick_cSpkProbThreshold_CSP,~gBlock2),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz)./sqrt(nIC), 'lineProps', 'k')
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstLickAlignEvents(:,firstLick_cSpkProbThreshold_CSM,gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), std(mean(firstLickAlignEvents(:,firstLick_cSpkProbThreshold_CSM,gBlock2),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz)./sqrt(nIC), 'lineProps', 'r')
        ylim([0 8]); vline(0,'k'); 
        xlabel('time from lick (s)');
        ylabel('Cspk probability (%)');
        title([num2str(sum(firstLick_cSpkProbThreshold_CSP)) ' CS+ and ' num2str(sum(firstLick_cSpkProbThreshold_CSM)) ' CS- responsive neurons from ' num2str(size(expt,2)) ' animals | lick aligned peak Cspk probability > ' num2str(probabilityThreshold)]);
        subplot(2,1,2); hold on;
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstLickAlignEvents(:,~firstLick_cSpkProbThreshold_CSP,~gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), std(mean(firstLickAlignEvents(:,~firstLick_cSpkProbThreshold_CSP,~gBlock2),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz)./sqrt(nIC), 'lineProps', 'k')
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstLickAlignEvents(:,~firstLick_cSpkProbThreshold_CSM,gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), std(mean(firstLickAlignEvents(:,~firstLick_cSpkProbThreshold_CSM,gBlock2),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz)./sqrt(nIC), 'lineProps', 'r')
        ylim([0 8]); vline(0,'k'); 
        xlabel('time from lick (s)');
        ylabel('Cspk probability (%)');
        title([num2str(sum(~firstLick_cSpkProbThreshold_CSP)) ' CS+ and ' num2str(sum(~firstLick_cSpkProbThreshold_CSM)) ' CS- nonresponsive neurons from ' num2str(size(expt,2)) ' animals | lick aligned peak Cspk probability < ' num2str(probabilityThreshold)]);
        saveas(tempFig, [output_fn 'stats\lickAlign_CspkProbability_' num2str(probabilityThreshold) '%threshold.pdf'],'pdf'); 
    end

        for c = 1:size(cueAlignEvents,2)
            for f = 1:size(cueAlignEvents,1)
                spikeProbability_CSP(f,c) = sum(cueAlignEvents(f,c,~gBlock2)==1)./sum(~isnan(cueAlignEvents(f,c,~gBlock2)));
                spikeProbability_CSM(f,c) = sum(cueAlignEvents(f,c,gBlock2)==1)./sum(~isnan(cueAlignEvents(f,c,gBlock2)));
            end
        end
    % Here, I am plotting the Cspk probability of cells aligned to cue
        tempFig=setFigure; hold on;
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(spikeProbability_CSP(36:96,:),2), std(spikeProbability_CSP(36:96,:),[],2)./sqrt(nIC), 'lineProps', 'k')
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(spikeProbability_CSM(36:96,:),2), std(spikeProbability_CSM(36:96,:),[],2)./sqrt(nIC), 'lineProps', 'r')
        xlim([0 1]); ylim([0 .15]); vline(0,'k'); vline(.767,'k')
        xlabel('time from cue (s)');
        ylabel('Cspk probability (%)');
        sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | cue aligned peak Cspk probability']);
        saveas(tempFig, [output_fn 'stats\cueAligned_CspkProbability.pdf'],'pdf'); 

    % Now I want a probability histogram for Cspk probability at peak        
        tempFig=setFigure; hold on;
        histogram(spikeProbability_CSP(mean(spikeProbability_CSP,2,'omitnan') == max(mean(spikeProbability_CSP,2,'omitnan')),:), 'BinWidth',0.05, 'Normalization','percentage', 'FaceColor','k' ,'EdgeColor','none');
        histogram(spikeProbability_CSM(mean(spikeProbability_CSM,2,'omitnan') == max(mean(spikeProbability_CSM,2,'omitnan')),:), 'BinWidth',0.05, 'Normalization','percentage', 'FaceColor','r' ,'EdgeColor','none');
        xlabel('Cspk probability (%)');
        ylabel('Proportion of cells (%)');
        sgtitle(['cue aligned Cspk probability histogram at peak: ' num2str(find(mean(spikeProbability_CSP,2,'omitnan') == max(mean(spikeProbability_CSP,2,'omitnan')))) ' CSP| ' num2str(find(mean(spikeProbability_CSM,2,'omitnan') == max(mean(spikeProbability_CSM,2,'omitnan')))) ' CSM']);
        saveas(tempFig, [output_fn 'stats\cueAligned_CspkProbabilityHistogram.pdf'],'pdf'); 
 
    % Normalize to log base 2 since spike probability is binary
        tempFig=setFigure; hold on;
        histogram(log2(spikeProbability_CSP(mean(spikeProbability_CSP,2,'omitnan') == max(mean(spikeProbability_CSP,2,'omitnan')),:)), 'BinWidth',0.5, 'Normalization','percentage', 'FaceColor','k' ,'EdgeColor','none');
        histogram(log2(spikeProbability_CSM(mean(spikeProbability_CSM,2,'omitnan') == max(mean(spikeProbability_CSM,2,'omitnan')),:)), 'BinWidth',0.5, 'Normalization','percentage', 'FaceColor','r' ,'EdgeColor','none');
        xlabel('Log 2 Cspk probability (%)');
        ylabel('Proportion of cells (%)');
        sgtitle(['lick aligned Cspk log(probability) histogram at peak: ' num2str(find(mean(spikeProbability_CSP,2,'omitnan') == max(mean(spikeProbability_CSP,2,'omitnan')))) ' CSP| ' num2str(find(mean(spikeProbability_CSM,2,'omitnan') == max(mean(spikeProbability_CSM,2,'omitnan')))) ' CSM']);
        saveas(tempFig, [output_fn 'stats\log_cueAligned_CspkProbabilityHistogram.pdf'],'pdf'); 
    
    % Plot Cspk rate of cells with Cspk probability > 0.10
        cSpkProbThreshold_CSP = mean(spikeProbability_CSP(find(mean(spikeProbability_CSP,2,'omitnan') == max(mean(spikeProbability_CSP,2,'omitnan')))-statsWindow:find(mean(spikeProbability_CSP,2,'omitnan') == max(mean(spikeProbability_CSP,2,'omitnan')))+statsWindow,:) > probabilityThreshold, 1,'omitnan') > 0;
        cSpkProbThreshold_CSM = mean(spikeProbability_CSM(find(mean(spikeProbability_CSM,2,'omitnan') == max(mean(spikeProbability_CSM,2,'omitnan')))-statsWindow:find(mean(spikeProbability_CSM,2,'omitnan') == max(mean(spikeProbability_CSM,2,'omitnan')))+statsWindow,:) > probabilityThreshold, 1,'omitnan') > 0;
        tempFig=setFigure; hold on;
        subplot(2,1,1); hold on;
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents(36:96,cSpkProbThreshold_CSP,~gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), std(mean(cueAlignEvents(36:96,cSpkProbThreshold_CSP,~gBlock2),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz)./sqrt(nIC), 'lineProps', 'k')
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents(36:96,cSpkProbThreshold_CSM,gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), std(mean(cueAlignEvents(36:96,cSpkProbThreshold_CSM,gBlock2),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz)./sqrt(nIC), 'lineProps', 'r')
        ylim([0 8]); vline(0,'k'); vline(.767,'k');
        xlabel('time from cue (s)');
        ylabel('Cspk probability (%)');
        title([num2str(sum(cSpkProbThreshold_CSP)) ' CS+ and ' num2str(sum(cSpkProbThreshold_CSM)) ' CS- responsive neurons from ' num2str(size(expt,2)) ' animals | cue aligned peak Cspk probability > ' num2str(probabilityThreshold)]);
        subplot(2,1,2); hold on;
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents(36:96,~cSpkProbThreshold_CSP,~gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), std(mean(cueAlignEvents(36:96,~cSpkProbThreshold_CSP,~gBlock2),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz)./sqrt(nIC), 'lineProps', 'k')
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents(36:96,~cSpkProbThreshold_CSM,gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), std(mean(cueAlignEvents(36:96,~cSpkProbThreshold_CSM,gBlock2),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz)./sqrt(nIC), 'lineProps', 'r')
        ylim([0 8]); vline(0,'k'); vline(.767,'k');
        xlabel('time from cue (s)');
        ylabel('Cspk probability (%)');
        title([num2str(sum(~cSpkProbThreshold_CSP)) ' CS+ and ' num2str(sum(~cSpkProbThreshold_CSM)) ' CS- nonresponsive neurons from ' num2str(size(expt,2)) ' animals | cue aligned peak Cspk probability < ' num2str(probabilityThreshold)]);
        saveas(tempFig, [output_fn 'stats\cueAlign_CspkProbability_' num2str(probabilityThreshold) '%threshold.pdf'],'pdf'); 
     
        for i = 1:splitBy
            for c = 1:size(cueAlignEvents,2)
                for f = 1:size(cueAlignEvents,1)
                spikeProbability_CSP_split{1,i}(f,c) = sum(cueAlignEvents(f,c,gRewSplit{1,i})==1)./sum(~isnan(cueAlignEvents(f,c,gRewSplit{1,i})));
                    spikeProbability_CSPtrials_split{1,i}(f,c,:) = cueAlignEvents(f,c,gRewSplit{1,i});
                spikeProbability_CSM_split{1,i}(f,c) = sum(cueAlignEvents(f,c,gB2Split{1,i})==1)./sum(~isnan(cueAlignEvents(f,c,gB2Split{1,i})));
                    spikeProbability_CSMtrials_split{1,i}(f,c,:) = cueAlignEvents(f,c,gB2Split{1,i});
                end
            end
            for c = 1:size(firstLickAlignEvents,2)
                for f = 1:size(firstLickAlignEvents,1)
                firstLick_spikeProbability_CSP_split{1,i}(f,c) = sum(firstLickAlignEvents(f,c,gRewSplit{1,i})==1)./sum(~isnan(firstLickAlignEvents(f,c,gRewSplit{1,i})));
                    firstLick_spikeProbability_CSPtrials_split{1,i}(f,c,:) = firstLickAlignEvents(f,c,gRewSplit{1,i});
                firstLick_spikeProbability_CSM_split{1,i}(f,c) = sum(firstLickAlignEvents(f,c,gB2Split{1,i})==1)./sum(~isnan(firstLickAlignEvents(f,c,gB2Split{1,i})));
                    firstLick_spikeProbability_CSMtrials_split{1,i}(f,c,:) = firstLickAlignEvents(f,c,gB2Split{1,i});
                end
            end
            for c = 1:size(firstPostRew_lickAlignEvents,2)
                for f = 1:size(firstPostRew_lickAlignEvents,1)
                firstPostRewLick_spikeProbability_CSP_split{1,i}(f,c) = sum(firstPostRew_lickAlignEvents(f,c,gRewSplit{1,i})==1)./sum(~isnan(firstPostRew_lickAlignEvents(f,c,gRewSplit{1,i})));
                    firstPostRewLick_spikeProbability_CSPtrials_split{1,i}(f,c,:) = firstPostRew_lickAlignEvents(f,c,gRewSplit{1,i});
                firstPostRewLick_spikeProbability_CSM_split{1,i}(f,c) = sum(firstPostRew_lickAlignEvents(f,c,gB2Split{1,i})==1)./sum(~isnan(firstPostRew_lickAlignEvents(f,c,gB2Split{1,i})));
                    firstPostRewLick_spikeProbability_CSMtrials_split{1,i}(f,c,:) = firstPostRew_lickAlignEvents(f,c,gB2Split{1,i});
                end
            end
        end
        
        for i = 1:splitBy
            for c = 1:size(cueAlignEvents,2)
                % for f = 1:size(cueAlignEvents,1)
                cSpkRate_CSP_split{1,i}(:,c) = mean(cueAlignEvents(:,c,gRewSplit{1,i}),3,'omitnan');%==1)./sum(~isnan(cueAlignEvents(:,c,gRewSplit{1,i})));
                    cSpkRate_CSPtrials_split{1,i}(:,c,:) = cueAlignEvents(:,c,gRewSplit{1,i});
                cSpkRate_CSM_split{1,i}(:,c) = mean(cueAlignEvents(:,c,gB2Split{1,i}),3,'omitnan');%==1)./sum(~isnan(cueAlignEvents(:,c,gB2Split{1,i})));
                    cSpkRate_CSMtrials_split{1,i}(:,c,:) = cueAlignEvents(:,c,gB2Split{1,i});
                % end
            end
            for c = 1:size(firstLickAlignEvents,2)
                %for f = 1:size(firstLickAlignEvents,1)
                firstLick_cSpkRate_CSP_split{1,i}(:,c) = mean(firstLickAlignEvents(:,c,gRewSplit{1,i}),3,'omitnan');%==1)./sum(~isnan(firstLickAlignEvents(:,c,gRewSplit{1,i})));
                    firstLick_cSpkRate_CSPtrials_split{1,i}(:,c,:) = firstLickAlignEvents(:,c,gRewSplit{1,i});
                firstLick_cSpkRate_CSM_split{1,i}(:,c) = mean(firstLickAlignEvents(:,c,gB2Split{1,i}),3,'omitnan');%==1)./sum(~isnan(firstLickAlignEvents(:,c,gB2Split{1,i})));
                    firstLick_cSpkRate_CSMtrials_split{1,i}(:,c,:) = firstLickAlignEvents(:,c,gB2Split{1,i});
                %end
            end
            for c = 1:size(firstPostRew_lickAlignEvents,2)
                %for f = 1:size(firstLickAlignEvents,1)
                firstPostRewLick_cSpkRate_CSP_split{1,i}(:,c) = mean(firstPostRew_lickAlignEvents(:,c,gRewSplit{1,i}),3,'omitnan');%==1)./sum(~isnan(firstLickAlignEvents(:,c,gRewSplit{1,i})));
                    firstPostRewLick_cSpkRate_CSPtrials_split{1,i}(:,c,:) = firstPostRew_lickAlignEvents(:,c,gRewSplit{1,i});
                firstPostRewLick_cSpkRate_CSM_split{1,i}(:,c) = mean(firstPostRew_lickAlignEvents(:,c,gB2Split{1,i}),3,'omitnan');%==1)./sum(~isnan(firstLickAlignEvents(:,c,gB2Split{1,i})));
                    firstPostRewLick_cSpkRate_CSMtrials_split{1,i}(:,c,:) = firstPostRew_lickAlignEvents(:,c,gB2Split{1,i});
                %end
            end
        end

        tempFig=setFigure; hold on;
        for i = 1:splitBy
        subplot(splitBy,1,i); hold on;
        histogram(spikeProbability_CSP_split{1,i}(mean(spikeProbability_CSP_split{1,i},2,'omitnan') == max(mean(spikeProbability_CSP_split{1,i},2,'omitnan')),:), 'BinWidth',0.05, 'Normalization','percentage', 'FaceColor','k' ,'EdgeColor','none');
        histogram(spikeProbability_CSM_split{1,i}(mean(spikeProbability_CSM_split{1,i},2,'omitnan') == max(mean(spikeProbability_CSM_split{1,i},2,'omitnan')),:), 'BinWidth',0.05, 'Normalization','percentage', 'FaceColor','r' ,'EdgeColor','none');
        xlabel('Cspk probability (%)');
        ylabel('Proportion of cells (%)');
        xlim([0 1.0]); ylim([0 50]);
        title([num2str(find(mean(spikeProbability_CSP_split{1,i},2,'omitnan') == max(mean(spikeProbability_CSP_split{1,i},2,'omitnan')))) ' CSP| ' num2str(find(mean(spikeProbability_CSM_split{1,i},2,'omitnan') == max(mean(spikeProbability_CSM_split{1,i},2,'omitnan')))) ' CSM']);
        end
        sgtitle('cue aligned Cspk probability histogram at peak: ');
        saveas(tempFig, [output_fn 'stats\cueAligned_CspkProbabilityHistogram.pdf'],'pdf'); 
 
        tempFig=setFigure; hold on; % align neural and licking data to the first and last lick, respectively
        for i = 1:splitBy
        subplot(splitBy,1,i); hold on;
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(cSpkRate_CSP_split{1,i}(36:96,:),2,'omitnan'), (std(cSpkRate_CSP_split{1,i}(36:96,:),[],2,'omitnan'))./sqrt(nIC),'lineProps','k');
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(cSpkRate_CSM_split{1,i}(36:96,:),2,'omitnan'), (std(cSpkRate_CSP_split{1,i}(36:96,:),[],2,'omitnan'))./sqrt(nIC),'lineProps','r');
        hold on;
        xlabel('time from cue (s)');
        ylabel('Cspk probability (%)');
        ylim([0 0.15]); vline(0,'k'); vline(.767,'k');
        title([]);
        end
        sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Cspk probability relative to cue']);
        saveas(tempFig, [output_fn 'stats\cueAligned_cSpkProbability_splitsBy' num2str(splitBy) '.pdf'],'pdf'); 
        
        for i = 1:splitBy
        cSpkProbThreshold_CSP_split{1,i} = mean(spikeProbability_CSP_split{1,i}(find(mean(spikeProbability_CSP_split{1,i},2,'omitnan') == max(mean(spikeProbability_CSP_split{1,i},2,'omitnan')))-statsWindow:find(mean(spikeProbability_CSP_split{1,i},2,'omitnan') == max(mean(spikeProbability_CSP_split{1,i},2,'omitnan')))+statsWindow,:) > probabilityThreshold, 1,'omitnan') > 0;
        cSpkProbThreshold_CSM_split{1,i} = mean(spikeProbability_CSM_split{1,i}(find(mean(spikeProbability_CSM_split{1,i},2,'omitnan') == max(mean(spikeProbability_CSM_split{1,i},2,'omitnan')))-statsWindow:find(mean(spikeProbability_CSM_split{1,i},2,'omitnan') == max(mean(spikeProbability_CSM_split{1,i},2,'omitnan')))+statsWindow,:) > probabilityThreshold, 1,'omitnan') > 0;
        end
        tempFig=setFigure; hold on; row=0; split=3;
        for i = 1:splitBy
        subplot(splitBy,2,i+row); hold on;
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents(36:96,cSpkProbThreshold_CSP_split{1,split},gRewSplit{1,i}),3,'omitnan'),2,'omitnan'), std(mean(cueAlignEvents(36:96,cSpkProbThreshold_CSP_split{1,split},gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC), 'lineProps', 'k')
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents(36:96,cSpkProbThreshold_CSM_split{1,split},gB2Split{1,i}),3,'omitnan'),2,'omitnan'), std(mean(cueAlignEvents(36:96,cSpkProbThreshold_CSM_split{1,split},gB2Split{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC), 'lineProps', 'r')
        ylim([0 .25]); vline(0,'k'); vline(.767,'k');
        xlabel('time from cue (s)');
        ylabel('Cspk probability (%)');
        title([num2str(sum(cSpkProbThreshold_CSP_split{1,split})) ' CS+ and ' num2str(sum(cSpkProbThreshold_CSM_split{1,split})) ' CS- responders']);
        subplot(splitBy,2,i+1+row); hold on;
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents(36:96,~cSpkProbThreshold_CSP_split{1,split},gRewSplit{1,i}),3,'omitnan'),2,'omitnan'), std(mean(cueAlignEvents(36:96,~cSpkProbThreshold_CSP_split{1,split},gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC), 'lineProps', 'k')
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents(36:96,~cSpkProbThreshold_CSM_split{1,split},gB2Split{1,i}),3,'omitnan'),2,'omitnan'), std(mean(cueAlignEvents(36:96,~cSpkProbThreshold_CSM_split{1,split},gB2Split{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC), 'lineProps', 'r')
        ylim([0 .25]); vline(0,'k'); vline(.767,'k');
        xlabel('time from cue (s)');
        ylabel('Cspk probability (%)');
        title([num2str(sum(~cSpkProbThreshold_CSP_split{1,split})) ' CS+ and ' num2str(sum(~cSpkProbThreshold_CSM_split{1,split})) ' CS- nonresponders']);
        row=i;
        end
        sgtitle(['Cspk probability across phases for cells with Cspk%>' num2str(probabilityThreshold) ' in phase ' num2str(split)])
        saveas(tempFig, [output_fn 'stats\cueAlign_CspkProbability_' num2str(probabilityThreshold) '%thresholdInSplit' num2str(split) '_splitBy' num2str(splitBy) '.pdf'],'pdf'); 
        
        % Here, I break lick align spike probability down by phase - first, a histogram
        tempFig=setFigure; hold on;
        for i = 1:splitBy
        subplot(splitBy,1,i); hold on;
        histogram(firstLick_spikeProbability_CSP_split{1,i}(mean(firstLick_spikeProbability_CSP_split{1,i},2,'omitnan') == max(mean(firstLick_spikeProbability_CSP_split{1,i},2,'omitnan')),:), 'BinWidth',0.05, 'Normalization','percentage', 'FaceColor','k' ,'EdgeColor','none');
        histogram(firstLick_spikeProbability_CSM_split{1,i}(mean(firstLick_spikeProbability_CSM_split{1,i},2,'omitnan') == max(mean(firstLick_spikeProbability_CSM_split{1,i},2,'omitnan')),:), 'BinWidth',0.05, 'Normalization','percentage', 'FaceColor','r' ,'EdgeColor','none');
        xlabel('Cspk probability (%)');
        ylabel('Proportion of cells (%)');
        xlim([0 0.8]); ylim([0 50]);
        title([num2str(find(mean(firstLick_spikeProbability_CSP_split{1,i},2,'omitnan') == max(mean(firstLick_spikeProbability_CSP_split{1,i},2,'omitnan')))) ' CSP| ' num2str(find(mean(firstLick_spikeProbability_CSM_split{1,i},2,'omitnan') == max(mean(firstLick_spikeProbability_CSM_split{1,i},2,'omitnan')))) ' CSM']);
        end
        sgtitle('lick aligned Cspk probability histogram at peak: ');
        saveas(tempFig, [output_fn 'stats\lickAligned_CspkProbabilityHistogram.pdf'],'pdf'); 
        % plot spike probability
        tempFig=setFigure; hold on; 
        for i = 1:splitBy
        subplot(splitBy,1,i); hold on;
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(firstLick_spikeProbability_CSP_split{1,i}(:,:),2,'omitnan'), (std(firstLick_spikeProbability_CSP_split{1,i}(:,:),[],2,'omitnan'))./sqrt(nIC),'lineProps','k');
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(firstLick_spikeProbability_CSM_split{1,i}(:,:),2,'omitnan'), (std(firstLick_spikeProbability_CSP_split{1,i}(:,:),[],2,'omitnan'))./sqrt(nIC),'lineProps','r');
        hold on;
        xlabel('time from lick (s)');
        ylabel('Cspk probability (%)');
        ylim([0 0.15]); vline(0,'k'); vline(.767,'k');
        title([]);
        end
        sgtitle([num2str(max(oDend)) ' neurons from ' num2str(size(expt,2)) ' animals | Cspk probability relative to lick']);
        saveas(tempFig, [output_fn 'stats\lickAligned_cSpkProbability_splitsBy' num2str(splitBy) '.pdf'],'pdf'); 
        
        for i = 1:splitBy
        firstLick_cSpkProbThreshold_CSP_split{1,i} = mean(firstLick_spikeProbability_CSP_split{1,i}((preLickRange(peakIdx_lickAlignSplitCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplitCSP(1,i))+statsWindow),:), 1,'omitnan');
        firstLick_cSpkProbThreshold_CSP_splitData{1,i} = firstLick_spikeProbability_CSP_split{1,i}((preLickRange(peakIdx_lickAlignSplitCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplitCSP(1,i))+statsWindow),:); % mean(, 1,'omitnan')
            firstLick_cSpkProbThreshold_CSP_splitAllData{1,i} = squeeze(firstLick_spikeProbability_CSPtrials_split{1,i}((preLickRange(peakIdx_lickAlignSplitCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplitCSP(1,i))+statsWindow),:,:));% mean(, 1,'omitnan');
        firstLick_cSpkProbThreshold_CSM_split{1,i} = mean(firstLick_spikeProbability_CSM_split{1,i}((preLickRange(peakIdx_lickAlignSplitCSM(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplitCSM(1,i))+statsWindow),:), 1,'omitnan');
        firstLick_cSpkProbThreshold_CSM_splitData{1,i} = firstLick_spikeProbability_CSM_split{1,i}((preLickRange(peakIdx_lickAlignSplitCSM(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplitCSM(1,i))+statsWindow),:);% mean(, 1,'omitnan');
            firstLick_cSpkProbThreshold_CSM_splitAllData{1,i} = squeeze(firstLick_spikeProbability_CSM_split{1,i}((preLickRange(peakIdx_lickAlignSplitCSM(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplitCSM(1,i))+statsWindow),:,:));% mean(, 1,'omitnan');
        end
        for i = 1:splitBy
        firstLick_cSpk_CSP_split{1,i} = mean(firstLick_cSpkRate_CSP_split{1,i}((preLickRange(peakIdx_lickAlignSplitCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplitCSP(1,i))+statsWindow),:), 1,'omitnan');
        firstLick_cSpk_CSP_splitData{1,i} = firstLick_cSpkRate_CSP_split{1,i}((preLickRange(peakIdx_lickAlignSplitCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplitCSP(1,i))+statsWindow),:); % mean(, 1,'omitnan')
            firstLick_cSpk_CSP_splitAllData{1,i} = squeeze(firstLick_cSpkRate_CSPtrials_split{1,i}((preLickRange(peakIdx_lickAlignSplitCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplitCSP(1,i))+statsWindow),:,:));% mean(, 1,'omitnan');
        firstLick_cSpk_CSM_split{1,i} = mean(firstLick_cSpkRate_CSM_split{1,i}((preLickRange(peakIdx_lickAlignSplitCSM(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplitCSM(1,i))+statsWindow),:), 1,'omitnan');
        firstLick_cSpk_CSM_splitData{1,i} = firstLick_cSpkRate_CSM_split{1,i}((preLickRange(peakIdx_lickAlignSplitCSM(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplitCSM(1,i))+statsWindow),:);% mean(, 1,'omitnan');
            firstLick_cSpk_CSM_splitAllData{1,i} = squeeze(firstLick_cSpkRate_CSM_split{1,i}((preLickRange(peakIdx_lickAlignSplitCSM(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplitCSM(1,i))+statsWindow),:,:));% mean(, 1,'omitnan');
        end

        tempFig=setFigure; hold on; row=0; split=1;
        for i = 1:splitBy
        subplot(splitBy,2,i+row); hold on;
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstLickAlignEvents(:,firstLick_cSpkProbThreshold_CSP_split{1,split}>probabilityThreshold,gRewSplit{1,i}),3,'omitnan'),2,'omitnan'), std(mean(firstLickAlignEvents(:,firstLick_cSpkProbThreshold_CSP_split{1,split}>probabilityThreshold,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC), 'lineProps', 'k')
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstLickAlignEvents(:,firstLick_cSpkProbThreshold_CSM_split{1,split}>probabilityThreshold,gB2Split{1,i}),3,'omitnan'),2,'omitnan'), std(mean(firstLickAlignEvents(:,firstLick_cSpkProbThreshold_CSM_split{1,split}>probabilityThreshold,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC), 'lineProps', 'r')
        ylim([0 .25]); vline(0,'k'); vline(.767,'k');
        xlabel('time from lick (s)');
        ylabel('Cspk probability (%)');
        title([num2str(sum(firstLick_cSpkProbThreshold_CSP_split{1,split}>probabilityThreshold)) ' CS+ and ' num2str(sum(firstLick_cSpkProbThreshold_CSM_split{1,split}>probabilityThreshold)) ' CS- responders']);
        subplot(splitBy,2,i+1+row); hold on;
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstLickAlignEvents(:,firstLick_cSpkProbThreshold_CSP_split{1,split}<=probabilityThreshold,gRewSplit{1,i}),3,'omitnan'),2,'omitnan'), std(mean(firstLickAlignEvents(:,firstLick_cSpkProbThreshold_CSP_split{1,split}<=probabilityThreshold,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC), 'lineProps', 'k')
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstLickAlignEvents(:,firstLick_cSpkProbThreshold_CSM_split{1,split}<=probabilityThreshold,gB2Split{1,i}),3,'omitnan'),2,'omitnan'), std(mean(firstLickAlignEvents(:,firstLick_cSpkProbThreshold_CSM_split{1,split}<=probabilityThreshold,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC), 'lineProps', 'r')
        ylim([0 .25]); vline(0,'k'); vline(.767,'k');
        xlabel('time from lick (s)');
        ylabel('Cspk probability (%)');
        title([num2str(sum(firstLick_cSpkProbThreshold_CSP_split{1,split}<=probabilityThreshold)) ' CS+ and ' num2str(sum(firstLick_cSpkProbThreshold_CSM_split{1,split}<=0.5)) ' CS- nonresponders']);
        row=i;
        end
        sgtitle(['Cspk probability across phases for cells with Cspk%>' num2str(probabilityThreshold) ' around first lick'])
        saveas(tempFig, [output_fn 'stats\lickAlign_CspkProbability_' num2str(probabilityThreshold) '%threshold_splitBy' num2str(splitBy) '.pdf'],'pdf'); 
     
        for i = 1:splitBy
        firstPostRewLick_cSpkProbThreshold_CSP_split{1,i} = mean(firstPostRewLick_spikeProbability_CSP_split{1,i}((preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))+statsWindow),:), 1,'omitnan');
        firstPostRewLick_cSpkProbThreshold_CSP_splitData{1,i} = firstPostRewLick_spikeProbability_CSP_split{1,i}((preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))+statsWindow),:); % mean(, 1,'omitnan')
            firstPostRewLick_cSpkProbThreshold_CSP_splitAllData{1,i} = squeeze(firstPostRewLick_spikeProbability_CSPtrials_split{1,i}((preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))+statsWindow),:,:));% mean(, 1,'omitnan');
        firstPostRewLick_cSpkProbThreshold_CSM_split{1,i} = mean(firstPostRewLick_spikeProbability_CSM_split{1,i}((preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))+statsWindow),:), 1,'omitnan');
        firstPostRewLick_cSpkProbThreshold_CSM_splitData{1,i} = firstPostRewLick_spikeProbability_CSM_split{1,i}((preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))+statsWindow),:);% mean(, 1,'omitnan');
            firstPostRewLick_cSpkProbThreshold_CSM_splitAllData{1,i} = squeeze(firstPostRewLick_spikeProbability_CSM_split{1,i}((preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))+statsWindow),:,:));% mean(, 1,'omitnan');
        end
        for i = 1:splitBy
        firstPostRewLick_cSpk_CSP_split{1,i} = mean(firstPostRewLick_cSpkRate_CSP_split{1,i}((preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))+statsWindow),:), 1,'omitnan');
        firstPostRewLick_cSpk_CSP_splitData{1,i} = firstPostRewLick_cSpkRate_CSP_split{1,i}((preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))+statsWindow),:); % mean(, 1,'omitnan')
            firstPostRewLick_cSpk_CSP_splitAllData{1,i} = squeeze(firstPostRewLick_cSpkRate_CSPtrials_split{1,i}((preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))+statsWindow),:,:));% mean(, 1,'omitnan');
        firstPostRewLick_cSpk_CSM_split{1,i} = mean(firstPostRewLick_cSpkRate_CSM_split{1,i}((preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))+statsWindow),:), 1,'omitnan');
        firstPostRewLick_cSpk_CSM_splitData{1,i} = firstPostRewLick_cSpkRate_CSM_split{1,i}((preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))+statsWindow),:);% mean(, 1,'omitnan');
            firstPostRewLick_cSpk_CSM_splitAllData{1,i} = squeeze(firstPostRewLick_cSpkRate_CSM_split{1,i}((preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))-statsWindow):(preLickRange(peakIdx_lickAlignSplit_postRewCSP(1,i))+statsWindow),:,:));% mean(, 1,'omitnan');
        end
        tempFig=setFigure; hold on; row=0; split=1;
        for i = 1:splitBy
        subplot(splitBy,2,i+row); hold on;
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstPostRew_lickAlignEvents(:,firstPostRewLick_cSpkProbThreshold_CSP_split{1,split}>probabilityThreshold,gRewSplit{1,i}),3,'omitnan'),2,'omitnan'), std(mean(firstPostRew_lickAlignEvents(:,firstPostRewLick_cSpkProbThreshold_CSP_split{1,split}>probabilityThreshold,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC), 'lineProps', 'k')
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstPostRew_lickAlignEvents(:,firstPostRewLick_cSpkProbThreshold_CSM_split{1,split}>probabilityThreshold,gB2Split{1,i}),3,'omitnan'),2,'omitnan'), std(mean(firstPostRew_lickAlignEvents(:,firstPostRewLick_cSpkProbThreshold_CSM_split{1,split}>probabilityThreshold,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC), 'lineProps', 'r')
        ylim([0 .25]); vline(0,'k'); vline(.767,'k');
        xlabel('time from lick (s)');
        ylabel('Cspk probability (%)');
        title([num2str(sum(firstPostRewLick_cSpkProbThreshold_CSP_split{1,split}>probabilityThreshold)) ' CS+ and ' num2str(sum(firstPostRewLick_cSpkProbThreshold_CSM_split{1,split}>probabilityThreshold)) ' CS- responders']);
        subplot(splitBy,2,i+1+row); hold on;
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstPostRew_lickAlignEvents(:,firstPostRewLick_cSpkProbThreshold_CSP_split{1,split}<=probabilityThreshold,gRewSplit{1,i}),3,'omitnan'),2,'omitnan'), std(mean(firstPostRew_lickAlignEvents(:,firstPostRewLick_cSpkProbThreshold_CSP_split{1,split}<=probabilityThreshold,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC), 'lineProps', 'k')
        shadedErrorBar_CV(tt(1,lickFrameRange)/1000, mean(mean(firstPostRew_lickAlignEvents(:,firstPostRewLick_cSpkProbThreshold_CSM_split{1,split}<=probabilityThreshold,gB2Split{1,i}),3,'omitnan'),2,'omitnan'), std(mean(firstPostRew_lickAlignEvents(:,firstPostRewLick_cSpkProbThreshold_CSM_split{1,split}<=probabilityThreshold,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC), 'lineProps', 'r')
        ylim([0 .25]); vline(0,'k'); vline(.767,'k');
        xlabel('time from lick (s)');
        ylabel('Cspk probability (%)');
        title([num2str(sum(firstPostRewLick_cSpkProbThreshold_CSP_split{1,split}<=probabilityThreshold)) ' CS+ and ' num2str(sum(firstPostRewLick_cSpkProbThreshold_CSM_split{1,split}<=0.5)) ' CS- nonresponders']);
        row=i;
        end
        sgtitle(['Cspk probability across phases for cells with Cspk%>' num2str(probabilityThreshold) ' around first lick post reward'])
        saveas(tempFig, [output_fn 'stats\postRewLickAlign_CspkProbability_' num2str(probabilityThreshold) '%threshold_splitBy' num2str(splitBy) '.pdf'],'pdf'); 
     
        for i = 1:splitBy
        cSpkProbThreshold_CSP_split{1,i} = mean(spikeProbability_CSP_split{1,i}((postCueRange(peakIdx_cueAlignSplit_postCueCSP(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSP(1,i))+statsWindow),:), 1,'omitnan');
        cSpkProbThreshold_CSP_splitData{1,i} = spikeProbability_CSP_split{1,i}((postCueRange(peakIdx_cueAlignSplit_postCueCSP(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSP(1,i))+statsWindow),:); % mean(, 1,'omitnan')
            cSpkProbThreshold_CSP_splitAllData{1,i} = squeeze(spikeProbability_CSPtrials_split{1,i}((postCueRange(peakIdx_cueAlignSplit_postCueCSP(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSP(1,i))+statsWindow),:,:));% mean(, 1,'omitnan');
        cSpkProbThreshold_CSM_split{1,i} = mean(spikeProbability_CSM_split{1,i}((postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))+statsWindow),:), 1,'omitnan');
        cSpkProbThreshold_CSM_splitData{1,i} = spikeProbability_CSM_split{1,i}((postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))+statsWindow),:);% mean(, 1,'omitnan');
            cSpkProbThreshold_CSM_splitAllData{1,i} = squeeze(spikeProbability_CSM_split{1,i}((postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))+statsWindow),:,:));% mean(, 1,'omitnan');
        end
        tempFig=setFigure; hold on; row=0; 
        for i = 1:splitBy
        subplot(splitBy,2,i+row); hold on;
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents(36:96,cSpkProbThreshold_CSP_split{1,i}>probabilityThreshold,gRewSplit{1,i}),3,'omitnan'),2,'omitnan'), std(mean(cueAlignEvents(36:96,cSpkProbThreshold_CSP_split{1,i}>probabilityThreshold,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC), 'lineProps', 'k')
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents(36:96,cSpkProbThreshold_CSM_split{1,i}>probabilityThreshold,gB2Split{1,i}),3,'omitnan'),2,'omitnan'), std(mean(cueAlignEvents(36:96,cSpkProbThreshold_CSM_split{1,i}>probabilityThreshold,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC), 'lineProps', 'r')
        ylim([0 .25]); vline(0,'k'); vline(.767,'k');
        xlabel('time from lick (s)');
        ylabel('Cspk probability (%)');
        title([num2str(sum(cSpkProbThreshold_CSP_split{1,i}>probabilityThreshold)) ' CS+ and ' num2str(sum(cSpkProbThreshold_CSM_split{1,i}>probabilityThreshold)) ' CS- responders']);
        subplot(splitBy,2,i+1+row); hold on;
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents(36:96,cSpkProbThreshold_CSP_split{1,i}<=probabilityThreshold,gRewSplit{1,i}),3,'omitnan'),2,'omitnan'), std(mean(cueAlignEvents(36:96,cSpkProbThreshold_CSP_split{1,i}<=probabilityThreshold,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC), 'lineProps', 'k')
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents(36:96,cSpkProbThreshold_CSM_split{1,i}<=probabilityThreshold,gB2Split{1,i}),3,'omitnan'),2,'omitnan'), std(mean(cueAlignEvents(36:96,cSpkProbThreshold_CSM_split{1,i}<=probabilityThreshold,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC), 'lineProps', 'r')
        ylim([0 .25]); vline(0,'k'); vline(.767,'k');
        xlabel('time from lick (s)');
        ylabel('Cspk probability (%)');
        title([num2str(sum(cSpkProbThreshold_CSP_split{1,i}<=probabilityThreshold)) ' CS+ and ' num2str(sum(cSpkProbThreshold_CSM_split{1,i}<=probabilityThreshold)) ' CS- nonresponders']);
        row=i;
        end
        sgtitle(['Cspk probability across phases for cells with Cspk%>' num2str(probabilityThreshold) ' after cue'])
        saveas(tempFig, [output_fn 'stats\cueAlign_CspkProbability_' num2str(probabilityThreshold) '%threshold_splitBy' num2str(splitBy) '.pdf'],'pdf'); 
     
        for i = 1:splitBy
        rewAlign_cSpkProbThreshold_CSP_split{1,i} = mean(spikeProbability_CSP_split{1,i}((postRewRange(peakIdx_cueAlignSplit_postRewCSP(1,i))-statsWindow):(postRewRange(peakIdx_cueAlignSplit_postRewCSP(1,i))+statsWindow),:), 1,'omitnan');
        rewAlign_cSpkProbThreshold_CSP_splitData{1,i} = spikeProbability_CSP_split{1,i}((postRewRange(peakIdx_cueAlignSplit_postRewCSP(1,i))-statsWindow):(postRewRange(peakIdx_cueAlignSplit_postRewCSP(1,i))+statsWindow),:); % mean(, 1,'omitnan')
            rewAlign_cSpkProbThreshold_CSP_splitAllData{1,i} = squeeze(spikeProbability_CSPtrials_split{1,i}((postRewRange(peakIdx_cueAlignSplit_postRewCSP(1,i))-statsWindow):(postRewRange(peakIdx_cueAlignSplit_postRewCSP(1,i))+statsWindow),:,:));% mean(, 1,'omitnan');
        rewAlign_cSpkProbThreshold_CSM_split{1,i} = mean(spikeProbability_CSM_split{1,i}((postRewRange(peakIdx_cueAlignSplit_postRewCSM(1,i))-statsWindow):(postRewRange(peakIdx_cueAlignSplit_postRewCSM(1,i))+statsWindow),:), 1,'omitnan');
        rewAlign_cSpkProbThreshold_CSM_splitData{1,i} = spikeProbability_CSM_split{1,i}((postRewRange(peakIdx_cueAlignSplit_postRewCSM(1,i))-statsWindow):(postRewRange(peakIdx_cueAlignSplit_postRewCSM(1,i))+statsWindow),:);% mean(, 1,'omitnan');
            rewAlign_cSpkProbThreshold_CSM_splitAllData{1,i} = squeeze(spikeProbability_CSM_split{1,i}((postRewRange(peakIdx_cueAlignSplit_postRewCSM(1,i))-statsWindow):(postRewRange(peakIdx_cueAlignSplit_postRewCSM(1,i))+statsWindow),:,:));% mean(, 1,'omitnan');
        end
        tempFig=setFigure; hold on; row=0; 
        for i = 1:splitBy
        subplot(splitBy,2,i+row); hold on;
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents(36:96,rewAlign_cSpkProbThreshold_CSP_split{1,i}>probabilityThreshold,gRewSplit{1,i}),3,'omitnan'),2,'omitnan'), std(mean(cueAlignEvents(36:96,rewAlign_cSpkProbThreshold_CSP_split{1,i}>probabilityThreshold,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC), 'lineProps', 'k')
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents(36:96,rewAlign_cSpkProbThreshold_CSM_split{1,i}>probabilityThreshold,gB2Split{1,i}),3,'omitnan'),2,'omitnan'), std(mean(cueAlignEvents(36:96,rewAlign_cSpkProbThreshold_CSM_split{1,i}>probabilityThreshold,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC), 'lineProps', 'r')
        ylim([0 .25]); vline(0,'k'); vline(.767,'k');
        xlabel('time from lick (s)');
        ylabel('Cspk probability (%)');
        title([num2str(sum(rewAlign_cSpkProbThreshold_CSP_split{1,i}>probabilityThreshold)) ' CS+ and ' num2str(sum(rewAlign_cSpkProbThreshold_CSM_split{1,i}>probabilityThreshold)) ' CS- responders']);
        subplot(splitBy,2,i+1+row); hold on;
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents(36:96,rewAlign_cSpkProbThreshold_CSP_split{1,i}<=probabilityThreshold,gRewSplit{1,i}),3,'omitnan'),2,'omitnan'), std(mean(cueAlignEvents(36:96,rewAlign_cSpkProbThreshold_CSP_split{1,i}<=probabilityThreshold,gRewSplit{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC), 'lineProps', 'k')
        shadedErrorBar_CV(tt(1,36:96)/1000, mean(mean(cueAlignEvents(36:96,rewAlign_cSpkProbThreshold_CSM_split{1,i}<=probabilityThreshold,gB2Split{1,i}),3,'omitnan'),2,'omitnan'), std(mean(cueAlignEvents(36:96,rewAlign_cSpkProbThreshold_CSM_split{1,i}<=probabilityThreshold,gB2Split{1,i}),3,'omitnan'),[],2,'omitnan')./sqrt(nIC), 'lineProps', 'r')
        ylim([0 .25]); vline(0,'k'); vline(.767,'k');
        xlabel('time from lick (s)');
        ylabel('Cspk probability (%)');
        title([num2str(sum(rewAlign_cSpkProbThreshold_CSP_split{1,i}<=probabilityThreshold)) ' CS+ and ' num2str(sum(rewAlign_cSpkProbThreshold_CSM_split{1,i}<=probabilityThreshold)) ' CS- nonresponders']);
        row=i;
        end
        sgtitle(['Cspk probability across phases for cells with Cspk%>' num2str(probabilityThreshold) ' after reward'])
        saveas(tempFig, [output_fn 'stats\rewAlign_CspkProbability_' num2str(probabilityThreshold) '%threshold_splitBy' num2str(splitBy) '.pdf'],'pdf'); 
     
    %% ANOVA for lick align comparisons across split

        phaseID={'first','second','third'}; nPhase=size(phaseID,2); cueID={'CS+','CS-'}; nCue=size(cueID,2); cells = 1:nIC; %nTrialsSplit=[size(firstPostRew_lick_spikeProbability_CSPtrials_split{1,1},3);size(firstLick_spikeProbability_CSMtrials_split{1,1},3)]; trialsSplit{1,1}=1:nTrialsSplit(1,1); trialsSplit{1,2}=1:nTrialsSplit(2,1);
        firstLick_cSpk_CSParray = [mean(firstLick_cSpk_CSP_splitData{1,1}*1000/frameRateHz,1,'omitnan');mean(firstLick_cSpk_CSP_splitData{1,2}*1000/frameRateHz,1,'omitnan');mean(firstLick_cSpk_CSP_splitData{1,3}*1000/frameRateHz,1,'omitnan')];
        firstLick_cSpk_CSMarray = [mean(firstLick_cSpk_CSM_splitData{1,1}*1000/frameRateHz,1,'omitnan');mean(firstLick_cSpk_CSM_splitData{1,2}*1000/frameRateHz,1,'omitnan');mean(firstLick_cSpk_CSM_splitData{1,3}*1000/frameRateHz,1,'omitnan')];
        firstLick_cSpk_array{1,1} = firstLick_cSpk_CSParray; firstLick_cSpk_array{1,2} = firstLick_cSpk_CSMarray;
        dataArray = NaN(nCue, nPhase, nIC);
%        dataArray(2, :, :, :) = NaN(nCue, nPhase, nIC);% Initialize the 3D array
        for i = 1:nCue
          dataArray(i, :, :) = firstLick_cSpk_array{1,i};
        end

%      
numCues = 2; % Assuming there are 2 cues
numPhases = 3;
numCells = nIC;
numMice = size(expt,2); % Number of mice

% Initialize cell arrays to store data for each factor
mouseData = cell(numMice, 1);
cueData = cell(numCues, 1);
phaseData = cell(numPhases, 1);

cellIdx = []; 
for i = 1:size(expt,2)
cellIdx = [cellIdx; ones(nsize(1,i),1)*i];
end
wholeMouseData=[];wholeCueData=[];wholeCellData=[];
firstLick_wholeResponseData=[];wholePhaseData=[];
for i = 1:size(cellIdx,1)
    for c=1:numCues
        for ph=1:numPhases
wholeMouseData(i,:) = cellIdx(i,1);
firstLick_wholeResponseData(i,:) = [dataArray(1,:,i), dataArray(2,:,i)];
wholeCueData = [wholeCueData, c];
cellID(i,:) = i;
wholeCellData = [wholeCellData i];
wholePhaseData = [wholePhaseData, ph];
        end
    end
end

dataForTable=[firstLick_wholeResponseData(:,1),firstLick_wholeResponseData(:,2),firstLick_wholeResponseData(:,3),firstLick_wholeResponseData(:,4),firstLick_wholeResponseData(:,5),firstLick_wholeResponseData(:,6)];
t = table(categorical(cellID),firstLick_wholeResponseData(:,1),firstLick_wholeResponseData(:,2),firstLick_wholeResponseData(:,3),firstLick_wholeResponseData(:,4),firstLick_wholeResponseData(:,5),firstLick_wholeResponseData(:,6),VariableNames=["cell","response1","response2","response3","response4","response5","response6"]);
% Create a table reflecting the within subject factors 
cue = categorical(repmat([1 1 1 2 2 2]', 1, 1));
phase = categorical(repmat([1 2 3 1 2 3]', 1, 1));
within = table(cue, phase, 'VariableNames', {'cue', 'phase'});
% 3. Call fitrm with the modified within design.
rm = fitrm(t,'response1-response6~1','WithinDesign',within);
[anovaTable_lickAlignSplit] = ranova(rm, 'WithinModel','cue*phase');

anovaTable_lickAlignSplit = movevars(anovaTable_lickAlignSplit, "DF", "Before", "SumSq");
writetable(anovaTable_lickAlignSplit,[output_fn 'stats\lickAlignSplit_two-wayANOVA.xlsx']);
disp('firstLickAlign')
disp(anovaTable_lickAlignSplit)

% !! CS+ only
t = table(categorical(cellID),firstLick_wholeResponseData(:,1),firstLick_wholeResponseData(:,2),firstLick_wholeResponseData(:,3),VariableNames=["cell","response1","response2","response3"]);
% Create a table reflecting the within subject factors 
cue = categorical(repmat([1 1 1]', 1, 1));
phase = categorical(repmat([1 2 3]', 1, 1));
within = table(phase, 'VariableNames', {'phase'});
% 3. Call fitrm with the modified within design.
rm = fitrm(t,'response1-response3~1','WithinDesign',within);
[anovaTable_lickAlignSplitCSP] = ranova(rm, 'WithinModel','phase');

anovaTable_lickAlignSplitCSP = movevars(anovaTable_lickAlignSplitCSP, "DF", "Before", "SumSq");
writetable(anovaTable_lickAlignSplitCSP,[output_fn 'stats\lickAlignSplit_CS+_two-wayANOVA.xlsx']);
disp(anovaTable_lickAlignSplitCSP)
% !! CS- only
t = table(categorical(cellID),firstLick_wholeResponseData(:,4),firstLick_wholeResponseData(:,5),firstLick_wholeResponseData(:,6),VariableNames=["cell","response1","response2","response3"]);
% Create a table reflecting the within subject factors 
cue = categorical(repmat([2 2 2]', 1, 1));
phase = categorical(repmat([1 2 3]', 1, 1));
within = table(phase, 'VariableNames', {'phase'});
% 3. Call fitrm with the modified within design.
rm = fitrm(t,'response1-response3~1','WithinDesign',within);
[anovaTable_lickAlignSplitCSM] = ranova(rm, 'WithinModel','phase');

anovaTable_lickAlignSplitCSM = movevars(anovaTable_lickAlignSplitCSM, "DF", "Before", "SumSq");
writetable(anovaTable_lickAlignSplitCSM,[output_fn 'stats\lickAlignSplit_CS-_two-wayANOVA.xlsx']);
disp(anovaTable_lickAlignSplitCSM)


    %% ANOVA for first lick after reward comparisons across split
        phaseID={'first','second','third'}; nPhase=size(phaseID,2); cueID={'CS+','CS-'}; nCue=size(cueID,2); cells = 1:nIC; %nTrialsSplit=[size(firstPostRewLick_spikeProbability_CSPtrials_split{1,1},3);size(firstLick_spikeProbability_CSMtrials_split{1,1},3)]; trialsSplit{1,1}=1:nTrialsSplit(1,1); trialsSplit{1,2}=1:nTrialsSplit(2,1);
        firstPostRewLick_cSpk_CSParray = [mean(firstPostRewLick_cSpk_CSP_splitData{1,1}*1000/frameRateHz,1,'omitnan');mean(firstPostRewLick_cSpk_CSP_splitData{1,2}*1000/frameRateHz,1,'omitnan');mean(firstPostRewLick_cSpk_CSP_splitData{1,3}*1000/frameRateHz,1,'omitnan')];
        firstPostRewLick_cSpk_CSMarray = [mean(firstPostRewLick_cSpk_CSM_splitData{1,1}*1000/frameRateHz,1,'omitnan');mean(firstPostRewLick_cSpk_CSM_splitData{1,2}*1000/frameRateHz,1,'omitnan');mean(firstPostRewLick_cSpk_CSM_splitData{1,3}*1000/frameRateHz,1,'omitnan')];
        firstPostRewLick_cSpk_array{1,1} = firstPostRewLick_cSpk_CSParray; firstPostRewLick_cSpk_array{1,2} = firstPostRewLick_cSpk_CSMarray;
        dataArray = NaN(nCue, nPhase, nIC);
%        dataArray(2, :, :, :) = NaN(nCue, nPhase, nIC);% Initialize the 3D array
        for i = 1:nCue
          dataArray(i, :, :) = firstPostRewLick_cSpk_array{1,i};
        end

%      
numCues = 2; % Assuming there are 2 cues
numPhases = 3;
numCells = nIC;
numMice = size(expt,2); % Number of mice

% Initialize cell arrays to store data for each factor
mouseData = cell(numMice, 1);
cueData = cell(numCues, 1);
phaseData = cell(numPhases, 1);

cellIdx = []; 
for i = 1:size(expt,2)
cellIdx = [cellIdx; ones(nsize(1,i),1)*i];
end
wholeMouseData=[];wholeCueData=[];wholeCellData=[];
firstPostRewLick_wholeResponseData=[];wholePhaseData=[];
for i = 1:size(cellIdx,1)
    for c=1:numCues
        for ph=1:numPhases
wholeMouseData(i,:) = cellIdx(i,1);
firstPostRewLick_wholeResponseData(i,:) = [dataArray(1,:,i), dataArray(2,:,i)];
wholeCueData = [wholeCueData, c];
cellID(i,:) = i;
wholeCellData = [wholeCellData i];
wholePhaseData = [wholePhaseData, ph];
        end
    end
end

dataForTable=[firstPostRewLick_wholeResponseData(:,1),firstPostRewLick_wholeResponseData(:,2),firstPostRewLick_wholeResponseData(:,3),firstPostRewLick_wholeResponseData(:,4),firstPostRewLick_wholeResponseData(:,5),firstPostRewLick_wholeResponseData(:,6)];
t = table(categorical(cellID),firstPostRewLick_wholeResponseData(:,1),firstPostRewLick_wholeResponseData(:,2),firstPostRewLick_wholeResponseData(:,3),firstPostRewLick_wholeResponseData(:,4),firstPostRewLick_wholeResponseData(:,5),firstPostRewLick_wholeResponseData(:,6),VariableNames=["cell","response1","response2","response3","response4","response5","response6"]);
% Create a table reflecting the within subject factors 
cue = categorical(repmat([1 1 1 2 2 2]', 1, 1));
phase = categorical(repmat([1 2 3 1 2 3]', 1, 1));
within = table(cue, phase, 'VariableNames', {'cue', 'phase'});
% 3. Call fitrm with the modified within design.
rm = fitrm(t,'response1-response6~1','WithinDesign',within);
[anovaTable_postRewLickAlignSplit] = ranova(rm, 'WithinModel','cue*phase');

anovaTable_postRewLickAlignSplit = movevars(anovaTable_postRewLickAlignSplit, "DF", "Before", "SumSq");
writetable(anovaTable_postRewLickAlignSplit,[output_fn 'stats\postRewLickAlignSplit_two-wayANOVA.xlsx']);
disp('postReward lickAlignSplit')
disp(anovaTable_postRewLickAlignSplit)

% !! CS+ only
t = table(categorical(cellID),firstPostRewLick_wholeResponseData(:,1),firstPostRewLick_wholeResponseData(:,2),firstPostRewLick_wholeResponseData(:,3),VariableNames=["cell","response1","response2","response3"]);
% Create a table reflecting the within subject factors 
cue = categorical(repmat([1 1 1]', 1, 1));
phase = categorical(repmat([1 2 3]', 1, 1));
within = table(phase, 'VariableNames', {'phase'});
% 3. Call fitrm with the modified within design.
rm = fitrm(t,'response1-response3~1','WithinDesign',within);
[anovaTable_postRewLickAlignSplitCSP] = ranova(rm, 'WithinModel','phase');

anovaTable_postRewLickAlignSplitCSP = movevars(anovaTable_postRewLickAlignSplitCSP, "DF", "Before", "SumSq");
writetable(anovaTable_postRewLickAlignSplitCSP,[output_fn 'stats\postRewLickAlignSplit_CS+_two-wayANOVA.xlsx']);
disp(anovaTable_postRewLickAlignSplitCSP)
% !! CS- only
t = table(categorical(cellID),firstPostRewLick_wholeResponseData(:,4),firstPostRewLick_wholeResponseData(:,5),firstPostRewLick_wholeResponseData(:,6),VariableNames=["cell","response1","response2","response3"]);
% Create a table reflecting the within subject factors 
cue = categorical(repmat([2 2 2]', 1, 1));
phase = categorical(repmat([1 2 3]', 1, 1));
within = table(phase, 'VariableNames', {'phase'});
% 3. Call fitrm with the modified within design.
rm = fitrm(t,'response1-response3~1','WithinDesign',within);
[anovaTable_postRewLickAlignSplitCSM] = ranova(rm, 'WithinModel','phase');

anovaTable_postRewLickAlignSplitCSM = movevars(anovaTable_postRewLickAlignSplitCSM, "DF", "Before", "SumSq");
writetable(anovaTable_postRewLickAlignSplitCSM,[output_fn 'stats\postRewLickAlignSplit_CS-_two-wayANOVA.xlsx']);
disp(anovaTable_postRewLickAlignSplitCSM)
    %% ANOVA for cue align comparisons across split

        for i = 1:splitBy
         if seshType==4
         cueAlign_cSpk_CSP_split{1,i} = mean(mean(cueAlignEvents((postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),1,'omitnan');
         cueAlign_cSpk_CSP_splitData{1,i} = mean(cueAlignEvents((postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan');
         else
         cueAlign_cSpk_CSP_split{1,i} = mean(mean(cueAlignEvents((postCueRange(peakIdx_cueAlignSplit_postCueCSP(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSP(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan'),1,'omitnan');
         cueAlign_cSpk_CSP_splitData{1,i} = mean(cueAlignEvents((postCueRange(peakIdx_cueAlignSplit_postCueCSP(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSP(1,i))+statsWindow),:,gRewSplit{1,i}),3,'omitnan');
         end
        cueAlign_cSpk_CSM_split{1,i} = mean(mean(cueAlignEvents((postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))+statsWindow),:,gB2Split{1,i}),3,'omitnan'),1,'omitnan');
        cueAlign_cSpk_CSM_splitData{1,i} = mean(cueAlignEvents((postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))-statsWindow):(postCueRange(peakIdx_cueAlignSplit_postCueCSM(1,i))+statsWindow),:,gB2Split{1,i}),3,'omitnan');
        end

        phaseID={'first','second','third'}; nPhase=size(phaseID,2); cueID={'CS+','CS-'}; nCue=size(cueID,2); cells = 1:nIC; %nTrialsSplit=[size(cueAlign_spikeProbability_CSPtrials_split{1,1},3);size(cueAlign_spikeProbability_CSMtrials_split{1,1},3)]; trialsSplit{1,1}=1:nTrialsSplit(1,1); trialsSplit{1,2}=1:nTrialsSplit(2,1);
        cueAlign_cSpk_CSParray = [mean(cueAlign_cSpk_CSP_splitData{1,1}*1000/frameRateHz,1,'omitnan');mean(cueAlign_cSpk_CSP_splitData{1,2}*1000/frameRateHz,1,'omitnan');mean(cueAlign_cSpk_CSP_splitData{1,3}*1000/frameRateHz,1,'omitnan')];
        cueAlign_cSpk_CSMarray = [mean(cueAlign_cSpk_CSM_splitData{1,1}*1000/frameRateHz,1,'omitnan');mean(cueAlign_cSpk_CSM_splitData{1,2}*1000/frameRateHz,1,'omitnan');mean(cueAlign_cSpk_CSM_splitData{1,3}*1000/frameRateHz,1,'omitnan')];
        cueAlign_cSpk_array{1,1} = cueAlign_cSpk_CSParray; cueAlign_cSpk_array{1,2} = cueAlign_cSpk_CSMarray;
        dataArray = NaN(nCue, nPhase, nIC);
%        dataArray(2, :, :, :) = NaN(nCue, nPhase, nIC);% Initialize the 3D array
        for i = 1:nCue
          dataArray(i, :, :) = cueAlign_cSpk_array{1,i};
        end

%      
numCues = 2; % Assuming there are 2 cues
numPhases = 3;
numCells = nIC;
numMice = size(expt,2); % Number of mice

% Initialize cell arrays to store data for each factor
mouseData = cell(numMice, 1); 
cueData = cell(numCues, 1);
phaseData = cell(numPhases, 1);

cellIdx = []; 
for i = 1:size(expt,2)
cellIdx = [cellIdx; ones(nsize(1,i),1)*i];
end
wholeMouseData=[];wholeCueData=[];wholeCellData=[];
cueAlign_wholeResponseData=[];wholePhaseData=[];
for i = 1:size(cellIdx,1)
    for c=1:numCues
        for ph=1:numPhases
wholeMouseData(i,:) = cellIdx(i,1);
cueAlign_wholeResponseData(i,:) = [dataArray(1,:,i), dataArray(2,:,i)];
wholeCueData = [wholeCueData, c];
cellID(i,:) = i;
wholeCellData = [wholeCellData i];
wholePhaseData = [wholePhaseData, ph];
        end
    end
end

dataForTable=[cueAlign_wholeResponseData(:,1),cueAlign_wholeResponseData(:,2),cueAlign_wholeResponseData(:,3),cueAlign_wholeResponseData(:,4),cueAlign_wholeResponseData(:,5),cueAlign_wholeResponseData(:,6)];
t = table(categorical(cellID),cueAlign_wholeResponseData(:,1),cueAlign_wholeResponseData(:,2),cueAlign_wholeResponseData(:,3),cueAlign_wholeResponseData(:,4),cueAlign_wholeResponseData(:,5),cueAlign_wholeResponseData(:,6),VariableNames=["cell","response1","response2","response3","response4","response5","response6"]);
% Create a table reflecting the within subject factors 
cue = categorical(repmat([1 1 1 2 2 2]', 1, 1));
phase = categorical(repmat([1 2 3 1 2 3]', 1, 1));
within = table(cue, phase, 'VariableNames', {'cue', 'phase'});
% 3. Call fitrm with the modified within design.
rm = fitrm(t,'response1-response6~1','WithinDesign',within);
[anovaTable_cueAlign] = ranova(rm, 'WithinModel','cue*phase');

anovaTable_cueAlign = movevars(anovaTable_cueAlign, "DF", "Before", "SumSq");
writetable(anovaTable_cueAlign,[output_fn 'stats\cueAlign_two-wayANOVA.xlsx']);
disp('postCueAlign')
disp(anovaTable_cueAlign)

% !! CS+ only
t = table(categorical(cellID),cueAlign_wholeResponseData(:,1),cueAlign_wholeResponseData(:,2),cueAlign_wholeResponseData(:,3),VariableNames=["cell","response1","response2","response3"]);
% Create a table reflecting the within subject factors 
cue = categorical(repmat([1 1 1]', 1, 1));
phase = categorical(repmat([1 2 3]', 1, 1));
within = table(phase, 'VariableNames', {'phase'});
% 3. Call fitrm with the modified within design.
rm = fitrm(t,'response1-response3~1','WithinDesign',within);
[anovaTable_cueAlignCSP] = ranova(rm, 'WithinModel','phase');

anovaTable_cueAlignCSP = movevars(anovaTable_cueAlignCSP, "DF", "Before", "SumSq");
writetable(anovaTable_cueAlignCSP,[output_fn 'stats\cueAlign_CS+_two-wayANOVA.xlsx']);
disp(anovaTable_cueAlignCSP)

% !! CS- only
t = table(categorical(cellID),cueAlign_wholeResponseData(:,4),cueAlign_wholeResponseData(:,5),cueAlign_wholeResponseData(:,6),VariableNames=["cell","response1","response2","response3"]);
% Create a table reflecting the within subject factors 
cue = categorical(repmat([2 2 2]', 1, 1));
phase = categorical(repmat([1 2 3]', 1, 1));
within = table(phase, 'VariableNames', {'phase'});
% 3. Call fitrm with the modified within design.
rm = fitrm(t,'response1-response3~1','WithinDesign',within);
[anovaTable_cueAlignCSM] = ranova(rm, 'WithinModel','phase');

anovaTable_cueAlignCSM = movevars(anovaTable_cueAlignCSM, "DF", "Before", "SumSq");
writetable(anovaTable_cueAlignCSM,[output_fn 'stats\cueAlign_CS-_two-wayANOVA.xlsx']);
disp(anovaTable_cueAlignCSM)

    %% ANOVA for reward align comparisons across split

     for i = 1:splitBy
        rewAlign_cSpk_CSP_split{1,i} = mean(cSpkRate_CSP_split{1,i}((postRewRange(peakIdx_cueAlignSplit_postRewCSP(1,i))-statsWindow):(postRewRange(peakIdx_cueAlignSplit_postRewCSP(1,i))+statsWindow),:), 1,'omitnan');
        rewAlign_cSpk_CSP_splitData{1,i} = cSpkRate_CSP_split{1,i}((postRewRange(peakIdx_cueAlignSplit_postRewCSP(1,i))-statsWindow):(postRewRange(peakIdx_cueAlignSplit_postRewCSP(1,i))+statsWindow),:);% mean(, 1,'omitnan');
            rewAlign_cSpk_CSP_splitAllData{1,i} = squeeze(cSpkRate_CSPtrials_split{1,i}((postRewRange(peakIdx_cueAlignSplit_postRewCSP(1,i))-statsWindow):(postRewRange(peakIdx_cueAlignSplit_postRewCSP(1,i))+statsWindow),:,:));% mean(, 1,'omitnan');
        rewAlign_cSpk_CSM_split{1,i} = mean(cSpkRate_CSM_split{1,i}((postRewRange(peakIdx_cueAlignSplit_postRewCSM(1,i))-statsWindow):(postRewRange(peakIdx_cueAlignSplit_postRewCSM(1,i))+statsWindow),:), 1,'omitnan');
        rewAlign_cSpk_CSM_splitData{1,i} = cSpkRate_CSM_split{1,i}((postRewRange(peakIdx_cueAlignSplit_postRewCSM(1,i))-statsWindow):(postRewRange(peakIdx_cueAlignSplit_postRewCSM(1,i))+statsWindow),:);% mean(, 1,'omitnan');
            rewAlign_cSpk_CSM_splitAllData{1,i} = squeeze(cSpkRate_CSM_split{1,i}((postRewRange(peakIdx_cueAlignSplit_postRewCSM(1,i))-statsWindow):(postRewRange(peakIdx_cueAlignSplit_postRewCSM(1,i))+statsWindow),:,:));% mean(, 1,'omitnan');
     end

        phaseID={'first','second','third'}; nPhase=size(phaseID,2); cueID={'CS+','CS-'}; nCue=size(cueID,2); cells = 1:nIC; %nTrialsSplit=[size(firstLick_spikeProbability_CSPtrials_split{1,1},3);size(firstLick_spikeProbability_CSMtrials_split{1,1},3)]; trialsSplit{1,1}=1:nTrialsSplit(1,1); trialsSplit{1,2}=1:nTrialsSplit(2,1);
        rewAlign_cSpk_CSParray = [mean(rewAlign_cSpk_CSP_splitData{1,1}*1000/frameRateHz,1,'omitnan');mean(rewAlign_cSpk_CSP_splitData{1,2}*1000/frameRateHz,1,'omitnan');mean(rewAlign_cSpk_CSP_splitData{1,3}*1000/frameRateHz,1,'omitnan')];
        rewAlign_cSpk_CSMarray = [mean(rewAlign_cSpk_CSM_splitData{1,1}*1000/frameRateHz,1,'omitnan');mean(rewAlign_cSpk_CSM_splitData{1,2}*1000/frameRateHz,1,'omitnan');mean(rewAlign_cSpk_CSM_splitData{1,3}*1000/frameRateHz,1,'omitnan')];
        rewAlign_cSpk_array{1,1} = rewAlign_cSpk_CSParray; rewAlign_cSpk_array{1,2} = rewAlign_cSpk_CSMarray;
        dataArray = NaN(nCue, nPhase, nIC);
%        dataArray(2, :, :, :) = NaN(nCue, nPhase, nIC);% Initialize the 3D array
        for i = 1:nCue
          dataArray(i, :, :) = rewAlign_cSpk_array{1,i};
        end

%      
numCues = 2; % Assuming there are 2 cues
numPhases = 3;
numCells = nIC;
numMice = size(expt,2); % Number of mice

% Initialize cell arrays to store data for each factor
mouseData = cell(numMice, 1); 
cueData = cell(numCues, 1);
phaseData = cell(numPhases, 1);

cellIdx = []; 
for i = 1:size(expt,2)
cellIdx = [cellIdx; ones(nsize(1,i),1)*i];
end
wholeMouseData=[];wholeCueData=[];wholeCellData=[];
rewAlign_wholeResponseData=[];wholePhaseData=[];
for i = 1:size(cellIdx,1)
    for c=1:numCues
        for ph=1:numPhases
wholeMouseData(i,:) = cellIdx(i,1);
rewAlign_wholeResponseData(i,:) = [dataArray(1,:,i), dataArray(2,:,i)];
wholeCueData = [wholeCueData, c];
cellID(i,:) = i;
wholeCellData = [wholeCellData i];
wholePhaseData = [wholePhaseData, ph];
        end
    end
end

dataForTable=[rewAlign_wholeResponseData(:,1),rewAlign_wholeResponseData(:,2),rewAlign_wholeResponseData(:,3),rewAlign_wholeResponseData(:,4),rewAlign_wholeResponseData(:,5),rewAlign_wholeResponseData(:,6)];
t = table(categorical(cellID),rewAlign_wholeResponseData(:,1),rewAlign_wholeResponseData(:,2),rewAlign_wholeResponseData(:,3),rewAlign_wholeResponseData(:,4),rewAlign_wholeResponseData(:,5),rewAlign_wholeResponseData(:,6),VariableNames=["cell","response1","response2","response3","response4","response5","response6"]);
% Create a table reflecting the within subject factors 
cue = categorical(repmat([1 1 1 2 2 2]', 1, 1));
phase = categorical(repmat([1 2 3 1 2 3]', 1, 1));
within = table(cue, phase, 'VariableNames', {'cue', 'phase'});
% 3. Call fitrm with the modified within design.
rm = fitrm(t,'response1-response6~1','WithinDesign',within);
[anovaTable_rewAlign] = ranova(rm, 'WithinModel','cue*phase');

anovaTable_rewAlign = movevars(anovaTable_rewAlign, "DF", "Before", "SumSq");
writetable(anovaTable_rewAlign,[output_fn 'stats\rewAlign_two-wayANOVA.xlsx']);
disp('postRewardAlign')
disp(anovaTable_rewAlign)

% !! CS+ only
t = table(categorical(cellID),rewAlign_wholeResponseData(:,1),rewAlign_wholeResponseData(:,2),rewAlign_wholeResponseData(:,3),VariableNames=["cell","response1","response2","response3"]);
% Create a table reflecting the within subject factors 
cue = categorical(repmat([1 1 1]', 1, 1));
phase = categorical(repmat([1 2 3]', 1, 1));
within = table(phase, 'VariableNames', {'phase'});
% 3. Call fitrm with the modified within design.
rm = fitrm(t,'response1-response3~1','WithinDesign',within);
[anovaTable_rewAlignCSP] = ranova(rm, 'WithinModel','phase');

anovaTable_rewAlignCSP = movevars(anovaTable_rewAlignCSP, "DF", "Before", "SumSq");
writetable(anovaTable_rewAlignCSP,[output_fn 'stats\rewAlign_CS+_two-wayANOVA.xlsx']);
disp(anovaTable_rewAlignCSP)
% !! CS- only
t = table(categorical(cellID),rewAlign_wholeResponseData(:,4),rewAlign_wholeResponseData(:,5),rewAlign_wholeResponseData(:,6),VariableNames=["cell","response1","response2","response3"]);
% Create a table reflecting the within subject factors 
cue = categorical(repmat([2 2 2]', 1, 1));
phase = categorical(repmat([1 2 3]', 1, 1));
within = table(phase, 'VariableNames', {'phase'});
% 3. Call fitrm with the modified within design.
rm = fitrm(t,'response1-response3~1','WithinDesign',within);
[anovaTable_rewAlignCSM] = ranova(rm, 'WithinModel','phase');

anovaTable_rewAlignCSM = movevars(anovaTable_rewAlignCSM, "DF", "Before", "SumSq");
writetable(anovaTable_rewAlignCSM,[output_fn 'stats\rewAlign_CS-_two-wayANOVA.xlsx']);
disp(anovaTable_rewAlignCSM)

% !! three-way ANOVA by condition alignment (lick, cue, rew)
dataForTable=[firstLick_wholeResponseData(:,1),firstLick_wholeResponseData(:,2),firstLick_wholeResponseData(:,3),firstLick_wholeResponseData(:,4),firstLick_wholeResponseData(:,5),firstLick_wholeResponseData(:,6),cueAlign_wholeResponseData(:,1),cueAlign_wholeResponseData(:,2),cueAlign_wholeResponseData(:,3),cueAlign_wholeResponseData(:,4),cueAlign_wholeResponseData(:,5),cueAlign_wholeResponseData(:,6),rewAlign_wholeResponseData(:,1),rewAlign_wholeResponseData(:,2),rewAlign_wholeResponseData(:,3),rewAlign_wholeResponseData(:,4),rewAlign_wholeResponseData(:,5),rewAlign_wholeResponseData(:,6)];
t = table(categorical(cellID),firstLick_wholeResponseData(:,1),firstLick_wholeResponseData(:,2),firstLick_wholeResponseData(:,3),firstLick_wholeResponseData(:,4),firstLick_wholeResponseData(:,5),firstLick_wholeResponseData(:,6),cueAlign_wholeResponseData(:,1),cueAlign_wholeResponseData(:,2),cueAlign_wholeResponseData(:,3),cueAlign_wholeResponseData(:,4),cueAlign_wholeResponseData(:,5),cueAlign_wholeResponseData(:,6),rewAlign_wholeResponseData(:,1),rewAlign_wholeResponseData(:,2),rewAlign_wholeResponseData(:,3),rewAlign_wholeResponseData(:,4),rewAlign_wholeResponseData(:,5),rewAlign_wholeResponseData(:,6),VariableNames=["cell","response1","response2","response3","response4","response5","response6","response7","response8","response9","response10","response11","response12","response13","response14","response15","response16","response17","response18"]);
% Create a table reflecting the within subject factors 
cue = categorical(repmat([1 1 1 2 2 2 1 1 1 2 2 2 1 1 1 2 2 2]', 1, 1));
phase = categorical(repmat([1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3]', 1, 1));
alignment = categorical(repmat([1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3]', 1, 1));
within = table(cue, phase, alignment, 'VariableNames', {'cue', 'phase', 'alignment'});
% 3. Call fitrm with the modified within design.
rm = fitrm(t,'response1-response18~1','WithinDesign',within);
[anovaTable_threeWay_alignment] = ranova(rm, 'WithinModel','cue*phase*alignment');

anovaTable_threeWay_alignment = movevars(anovaTable_threeWay_alignment, "DF", "Before", "SumSq");
writetable(anovaTable_threeWay_alignment,[output_fn 'stats\allAlignment_three-wayANOVA.xlsx']);
disp(anovaTable_threeWay_alignment)
end
%% plot     | first-piezo (response/suppression) 

  %we're gonna set an eval loop here because individual blocks for first-postCue-postRew for CS+&- is vile
eventLabel={'first','postCue','postRew'}; cueLabel={'CSP','CSM'}; prePiezo_frames=15; prePiezoRange=(1:bx_response_frames)+(prePiezo_frames-bx_response_frames);
for cue=1:size(cueLabel,2)
  for event=1:size(eventLabel,2)
      eval(['dendrite.piezoAlign.response_' eventLabel{1,event} cueLabel{1,cue} '=[];']);
      eval(['dendrite.piezoAlign.suppress_' eventLabel{1,event} cueLabel{1,cue} '=[];']);
    for c = 1:sum(nsize)
%
 eval(['prePiezoVar = ' eventLabel{1,event} '_piezoAlignDFoF' cueLabel{1,cue} '_acrossTrials(prePiezoRange,c);']);
 eval(['baselineVar = ' cueLabel{1,cue} 'dFoF_acrossTrials(baselineRange,c);']);
 eval(['pos.Diff.piezoAlign{c,1} = ([0; diff(' eventLabel{1,event} '_piezoAlignDFoF' cueLabel{1,cue} '_acrossTrials(:,c))]);']);
 
pos.Limit.piezoAlign{c,1}.one = (mean(pos.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)>0])),1,'omitnan')+(std(pos.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)>0])),[],1,'omitnan').*1)); 
pos.Change.piezoAlign{c,1}.one = find(pos.Diff.piezoAlign{c,1}>pos.Limit.piezoAlign{c,1}.one); 
pos.Limit.piezoAlign{c,1}.two = (mean(pos.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)>0])),1,'omitnan')+(std(pos.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)>0])),[],1,'omitnan').*1.5)); 
pos.Change.piezoAlign{c,1}.two = find(pos.Diff.piezoAlign{c,1}>pos.Limit.piezoAlign{c,1}.two); 
pos.Limit.piezoAlign{c,1}.three = (mean(pos.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)>0])),1,'omitnan')+(std(pos.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)>0])),[],1,'omitnan').*2)); 
pos.Change.piezoAlign{c,1}.three = find(pos.Diff.piezoAlign{c,1}>pos.Limit.piezoAlign{c,1}.three);

eval(['frameIdx.piezoAlign.' eventLabel{1,event} cueLabel{1,cue} '_posChangeDiff{1,c} = intersect(prePiezoRange, pos.Change.piezoAlign{c,1}.one);']); 
eval(['frameIdx.piezoAlign.' eventLabel{1,event} cueLabel{1,cue} '_posChangeDiff{2,c} = intersect(prePiezoRange, pos.Change.piezoAlign{c,1}.two);']); 
eval(['frameIdx.piezoAlign.' eventLabel{1,event} cueLabel{1,cue} '_posChangeDiff{3,c} = intersect(prePiezoRange, pos.Change.piezoAlign{c,1}.three);']);

 eval(['neg.Diff.piezoAlign{c,1} = ([0; diff(' eventLabel{1,event} '_piezoAlignDFoF' cueLabel{1,cue} '_acrossTrials(:,c))]);']);
neg.Limit.piezoAlign{c,1}.one = (mean(neg.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)<0])),1,'omitnan')-(std(neg.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)<0])),[],1,'omitnan').*1)); 
neg.Change.piezoAlign{c,1}.one = find(neg.Diff.piezoAlign{c,1}<neg.Limit.piezoAlign{c,1}.one); 
neg.Limit.piezoAlign{c,1}.two = (mean(neg.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)<0])),1,'omitnan')-(std(neg.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)<0])),[],1,'omitnan').*1.5)); 
neg.Change.piezoAlign{c,1}.two = find(neg.Diff.piezoAlign{c,1}<neg.Limit.piezoAlign{c,1}.two); 
neg.Limit.piezoAlign{c,1}.three = (mean(neg.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)<0])),1,'omitnan')-(std(neg.Diff.piezoAlign{c,1}(logical([0; diff(prePiezoVar)<0])),[],1,'omitnan').*2)); 
neg.Change.piezoAlign{c,1}.three = find(neg.Diff.piezoAlign{c,1}<neg.Limit.piezoAlign{c,1}.three);

eval(['frameIdx.piezoAlign.' eventLabel{1,event} cueLabel{1,cue} '_negChangeDiff{1,c} = intersect(prePiezoRange, neg.Change.piezoAlign{c,1}.one);']); 
eval(['frameIdx.piezoAlign.' eventLabel{1,event} cueLabel{1,cue} '_negChangeDiff{2,c} = intersect(prePiezoRange, neg.Change.piezoAlign{c,1}.two);']); 
eval(['frameIdx.piezoAlign.' eventLabel{1,event} cueLabel{1,cue} '_negChangeDiff{3,c} = intersect(prePiezoRange, neg.Change.piezoAlign{c,1}.three);']);
%    
    tempResp_diffIdx = diff(intersect(prePiezoRange, find(prePiezoVar>(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*respReq_piezoAlign))))';
    tempResp_frameIdx = intersect(prePiezoRange, find(prePiezoVar>(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*respReq_piezoAlign)));
    if length(tempResp_diffIdx)
      if strfind(tempResp_diffIdx, [1 1 1])
        if eval(['frameIdx.piezoAlign.' eventLabel{1,event} cueLabel{1,cue} '_posChangeDiff{chngReq_piezoAlign,c}'])
           eval(['dendrite.piezoAlign.response_' eventLabel{1,event} cueLabel{1,cue} '=[dendrite.piezoAlign.response_' eventLabel{1,event} cueLabel{1,cue} ' c];']);
        end
      end
    end
    tempSupp_diffIdx = diff(intersect(prePiezoRange, find(prePiezoVar<(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*suppReq_piezoAlign))))';
    tempSupp_frameIdx = intersect(prePiezoRange, find(prePiezoVar<(mean(baselineVar,1,'omitnan')+std(baselineVar,[],1,'omitnan').*suppReq_piezoAlign)));
    if length(tempSupp_diffIdx)
      if strfind(tempSupp_diffIdx, [1 1 1 1 1 1]) 
        if isempty(strfind(tempResp_diffIdx, [1 1 1]))
          if eval(['frameIdx.piezoAlign.' eventLabel{1,event} cueLabel{1,cue} '_negChangeDiff{chngReq_piezoAlign,c}'])
           eval(['dendrite.piezoAlign.suppress_' eventLabel{1,event} cueLabel{1,cue} '=[dendrite.piezoAlign.suppress_' eventLabel{1,event} cueLabel{1,cue} ' c];']);
          end
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
eval(['dendriteDiff.piezoAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} '= setdiff(dendrite.piezoAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} ',dendrite.piezoAlign.' sideLabel{2,side} '_'  eventLabel{1,event} cueLabel{1,cue} ');']);
it=it+1;
    end
  end
end
%
titleArray = {'CS+ reponsive firstPiezo','CS+ suppressive firstPiezo','CS+ reponsive postCuePiezo','CS+ suppressive postCuePiezo','CS+ reponsive postRewPiezo','CS+ suppressive postRewPiezo','CS- reponsive firstPiezo','CS- suppressive firstPiezo','CS- reponsive postCuePiezo','CS- suppressive postCuePiezo','CS- reponsive postRewPiezo','CS- suppressive postRewPiezo'};
colorArray = {'k','b','k','b','k','b','r','m','r','m','r','m'};

highlight = [255/255 217/255 96/255]; 
tempFig = setFigure; hold on; c=1;
for cue=1:size(cueLabel,2)
  for event=1:size(eventLabel,2)
    for side=1:size(sideLabel,2)
    eval(['callTrue = ' eventLabel{1,event} '_piezoAlignDFoF' cueLabel{1,cue} '_acrossTrials(:,dendriteDiff.piezoAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} ');']);
    eval(['callFalse = ' eventLabel{1,event} '_piezoAlignDFoF' cueLabel{1,cue} '_acrossTrials(:,setdiff(allcells,dendriteDiff.piezoAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} '));']);
    eval(['dendrites = dendriteDiff.piezoAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} ';']);
subplot(6,2,c); hold on;
ylim([ymin*2 ymax*2]); xline(0,'k');  title(titleArray{1,c});
rectangle('Position',[(tt(1,baselineRange(1))+10)/1000,ymin*2+.25,((1000./frameRateHz)*response_frames-10)/1000,ymax+5],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:80)/1000,mean(callTrue,2,'omitnan').*(1000./frameRateHz),std(callTrue,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',colorArray{1,c});
% shadedErrorBar_CV(tt(1,36:96)/1000,mean(callFalse,2,'omitnan'),std(callFalse,[],2,'omitnan')./sqrt(size(dendrites,2)),'lineProps',{'Color',[0/255 0/255 201/255]});
    c=c+1;
    end
  end
end
    saveas(tempFig, [output_fn 'stats\dFoF_piezoAlign_PSTHs.pdf'],'pdf');

tempFig = setFigure; hold on; c=1;
for cue=1:size(cueLabel,2)
  for event=1:size(eventLabel,2)
    for side=1:size(sideLabel,2)
    eval(['callTrue = ' eventLabel{1,event} '_piezoAlignEvents' cueLabel{1,cue} '_acrossTrials(:,dendriteDiff.piezoAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} ');']);
    eval(['callFalse = ' eventLabel{1,event} '_piezoAlignEvents' cueLabel{1,cue} '_acrossTrials(:,setdiff(allcells,dendriteDiff.piezoAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} '));']);
    eval(['dendrites = dendriteDiff.piezoAlign.' sideLabel{1,side} '_' eventLabel{1,event} cueLabel{1,cue} ';']);
subplot(6,2,c); hold on;
ylim([0 ymax]); xline(0,'k');  title(titleArray{1,c});
if c<5; rectangle('Position',[tt(1,baselineRange(1))/1000,0+.15,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight); elseif c>4; rectangle('Position',[tt(1,baselineRange(1))/1000,0+.15,(1000./frameRateHz)*response_frames/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight); end
shadedErrorBar_CV(tt(1,36:80)/1000,mean(callTrue,2,'omitnan').*(1000./frameRateHz),std(callTrue,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',colorArray{1,c});
% shadedErrorBar_CV(tt(1,36:96)/1000,mean(callFalse,2,'omitnan').*(1000./frameRateHz),std(callFalse,[],2,'omitnan')./sqrt(size(dendrites,2)).*(1000./frameRateHz),'lineProps',{'Color',[0/255 0/255 201/255]});
    c=c+1;
    end
  end
end
    saveas(tempFig, [output_fn 'stats\CspkEvents_piezoAlign_PSTHs.pdf'],'pdf');

    tempFig = setFigure; hold on; c=1;
for event=1:size(eventLabel,2)
  for cue=1:size(cueLabel,2)
eval(['respN = size(' eventLabel{1,event} '_piezoAlignEvents' cueLabel{1,cue} '_acrossTrials(:,dendriteDiff.piezoAlign.response_' eventLabel{1,event} cueLabel{1,cue} '),2);']);
eval(['suppN = size(' eventLabel{1,event} '_piezoAlignEvents' cueLabel{1,cue} '_acrossTrials(:,dendriteDiff.piezoAlign.suppress_' eventLabel{1,event} cueLabel{1,cue} '),2);']);
       elseN = sum(nsize)-(respN+suppN);

subplot(3,2,c); 
pie([respN, suppN, elseN]);
c=c+1; title([eventLabel{1,event} '-' cueLabel{1,cue}]);
    
  end
end
    saveas(tempFig, [output_fn 'stats\pie_propCells_piezoAlign.pdf'],'pdf');
    
tempFig=setFigure; hold on;
subplot(3,1,1); hold on;
shadedErrorBar_CV(tt(1,36:80)/1000,mean(first_piezoAlignEventsCSP_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(first_piezoAlignEventsCSP_acrossTrials,[],2,'omitnan')./sqrt(size(first_piezoAlignEventsCSP_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:80)/1000,mean(first_piezoAlignEventsCSM_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(first_piezoAlignEventsCSM_acrossTrials,[],2,'omitnan')./sqrt(size(first_piezoAlignEventsCSM_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
ylim([0 ymax]); xline(0,'k');  title('first piezo');
subplot(3,1,2); hold on;
shadedErrorBar_CV(tt(1,36:80)/1000,mean(postCue_piezoAlignEventsCSP_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(postCue_piezoAlignEventsCSP_acrossTrials,[],2,'omitnan')./sqrt(size(postCue_piezoAlignEventsCSP_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:80)/1000,mean(postCue_piezoAlignEventsCSM_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(postCue_piezoAlignEventsCSM_acrossTrials,[],2,'omitnan')./sqrt(size(postCue_piezoAlignEventsCSM_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
ylim([0 ymax]); xline(0,'k');  title('first postCue piezo');
subplot(3,1,3); hold on;
shadedErrorBar_CV(tt(1,36:80)/1000,mean(postRew_piezoAlignEventsCSP_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(postRew_piezoAlignEventsCSP_acrossTrials,[],2,'omitnan')./sqrt(size(postRew_piezoAlignEventsCSP_acrossTrials,2)).*(1000./frameRateHz),'lineProps','k');
shadedErrorBar_CV(tt(1,36:80)/1000,mean(postRew_piezoAlignEventsCSM_acrossTrials,2,'omitnan').*(1000./frameRateHz),std(postRew_piezoAlignEventsCSM_acrossTrials,[],2,'omitnan')./sqrt(size(postRew_piezoAlignEventsCSM_acrossTrials,2)).*(1000./frameRateHz),'lineProps','r');
ylim([0 ymax]); xline(0,'k');  title('first postRew piezo');
%% ttest    | first-piezo avg vs. baseline avg 
% offset start frame by 1, because we need to take a +/-1 frame window
% around the target; script breaks if the peak is on frame 1
  [peak_firstPiezoActCSP, peakIdx_firstPiezoActCSP] = max(mean(first_piezoAlignEventsCSP_acrossTrials(prePiezoRange,:),2,'omitnan'),[],1);
  [peak_postCuePiezoActCSP, peakIdx_postCuePiezoActCSP] = max(mean(postCue_piezoAlignEventsCSP_acrossTrials(prePiezoRange,:),2,'omitnan'),[],1);
  [peak_postRewPiezoActCSP, peakIdx_postRewPiezoActCSP] = max(mean(postRew_piezoAlignEventsCSP_acrossTrials(prePiezoRange,:),2,'omitnan'),[],1);
  [peak_baseCSP, peakIdx_baseCSP] = max(mean(CSPevents_acrossTrials(baselineRange,:),2,'omitnan'),[],1);
  [peak_firstPiezoActCSM, peakIdx_firstPiezoActCSM] = max(mean(first_piezoAlignEventsCSM_acrossTrials(prePiezoRange,:),2,'omitnan'),[],1);
  [peak_postCuePiezoActCSM, peakIdx_postCuePiezoActCSM] = max(mean(postCue_piezoAlignEventsCSM_acrossTrials(prePiezoRange,:),2,'omitnan'),[],1);
  [peak_postRewPiezoActCSM, peakIdx_postRewPiezoActCSM] = max(mean(postRew_piezoAlignEventsCSM_acrossTrials(prePiezoRange,:),2,'omitnan'),[],1);
  [peak_baseCSM, peakIdx_baseCSM] = max(mean(CSMevents_acrossTrials(baselineRange,:),2,'omitnan'),[],1);
% counts for piezo align
n.piezoAlign.CSP_first_trialAvg_count = sum(first_piezoAlignCSP_acrossTrials); 
n.piezoAlign.CSM_first_trialAvg_count = sum(first_piezoAlignCSM_acrossTrials);
n.piezoAlign.CSP_firstPostCue_trialAvg_count = sum(firstPostCue_piezoAlignCSP_acrossTrials);
n.piezoAlign.CSM_firstPostCue_trialAvg_count = sum(firstPostCue_piezoAlignCSM_acrossTrials);
n.piezoAlign.CSP_firstPostRew_trialAvg_count = sum(firstPostRew_piezoAlignCSP_acrossTrials);
n.piezoAlign.CSM_firstPostRew_trialAvg_count = sum(firstPostRew_piezoAlignCSM_acrossTrials);

for c = 1:iDend
firstPiezo_baseline_CSP(:,c) = CSPevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSP-statsWindow):baselineRange(1)-1+(peakIdx_baseCSP+statsWindow),c);
firstPiezo_activity_CSP(:,c) = first_piezoAlignEventsCSP_acrossTrials(prePiezoRange(1)-1+(peakIdx_firstPiezoActCSP-statsWindow):prePiezoRange(1)-1+(peakIdx_firstPiezoActCSP+statsWindow),c);
end
[~, p.piezoAlign.CSP_first_trialAvg_Cspk, ~, stats.piezoAlign.CSP_first_trialAvg_Cspk] = ttest(mean(firstPiezo_baseline_CSP,1,'omitnan'),mean(firstPiezo_activity_CSP,1,'omitnan'));
  m.piezoAlign.CSP_first_trialAvg_Cspk = mean(mean(firstPiezo_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoAlign.CSP_first_trialAvg_Cspk = std(mean(firstPiezo_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(firstPiezo_activity_CSP,2)).*(1000./frameRateHz);
for c = 1:iDend
firstPiezo_baseline_CSM(:,c) = CSMevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSM-statsWindow):baselineRange(1)-1+(peakIdx_baseCSM+statsWindow),c);
firstPiezo_activity_CSM(:,c) = first_piezoAlignEventsCSM_acrossTrials(prePiezoRange(1)-1+(peakIdx_firstPiezoActCSM-statsWindow):prePiezoRange(1)-1+(peakIdx_firstPiezoActCSM+statsWindow),c);
end
[~, p.piezoAlign.CSM_first_trialAvg_Cspk, ~, stats.piezoAlign.CSM_first_trialAvg_Cspk] = ttest(mean(firstPiezo_baseline_CSM,1,'omitnan'),mean(firstPiezo_activity_CSM,1,'omitnan'));
  m.piezoAlign.CSM_first_trialAvg_Cspk = mean(mean(firstPiezo_activity_CSM,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoAlign.CSM_first_trialAvg_Cspk = std(mean(firstPiezo_activity_CSM,1,'omitnan'),[],2,'omitnan')./sqrt(size(firstPiezo_activity_CSM,2)).*(1000./frameRateHz);

% difference in activity between first piezo [POST REWARD] and baseline
for c = 1:iDend
postRewPiezo_baseline_CSP(:,c) = CSPevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSP-statsWindow):baselineRange(1)-1+(peakIdx_baseCSP+statsWindow),c);
postRewPiezo_activity_CSP(:,c) = postRew_piezoAlignEventsCSP_acrossTrials(prePiezoRange(1)-1+(peakIdx_postRewPiezoActCSP-statsWindow):prePiezoRange(1)-1+(peakIdx_postRewPiezoActCSP+statsWindow),c);
end
[~, p.piezoAlign.CSP_postRew_trialAvg_Cspk, ~, stats.piezoAlign.CSP_postRew_trialAvg_Cspk] = ttest(mean(postRewPiezo_baseline_CSP,1,'omitnan'),mean(postRewPiezo_activity_CSP,1,'omitnan'));
  m.piezoAlign.CSP_postRew_trialAvg_Cspk = mean(mean(postRewPiezo_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoAlign.CSP_postRew_trialAvg_Cspk = std(mean(postRewPiezo_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(postRewPiezo_activity_CSP,2)).*(1000./frameRateHz);
for c = 1:iDend
postRewPiezo_baseline_CSM(:,c) = CSMevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSM-statsWindow):baselineRange(1)-1+(peakIdx_baseCSM+statsWindow),c);
postRewPiezo_activity_CSM(:,c) = postRew_piezoAlignEventsCSM_acrossTrials(prePiezoRange(1)-1+(peakIdx_postRewPiezoActCSM-statsWindow):prePiezoRange(1)-1+(peakIdx_postRewPiezoActCSM+statsWindow),c);
end
[~, p.piezoAlign.CSM_postRew_trialAvg_Cspk, ~, stats.piezoAlign.CSM_postRew_trialAvg_Cspk] = ttest(mean(postRewPiezo_baseline_CSM,1,'omitnan'),mean(postRewPiezo_activity_CSM,1,'omitnan'));
  m.piezoAlign.CSM_postRew_trialAvg_Cspk = mean(mean(postRewPiezo_activity_CSM,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoAlign.CSM_postRew_trialAvg_Cspk = std(mean(postRewPiezo_activity_CSM,1,'omitnan'),[],2,'omitnan')./sqrt(size(postRewPiezo_activity_CSM,2)).*(1000./frameRateHz);
  
% difference in activity between first piezo [POST CUE] and baseline
for c = 1:iDend
postCuePiezo_baseline_CSP(:,c) = CSPevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSP-statsWindow):baselineRange(1)-1+(peakIdx_baseCSP+statsWindow),c);
postCuePiezo_activity_CSP(:,c) = postCue_piezoAlignEventsCSP_acrossTrials(prePiezoRange(1)-1+(peakIdx_postCuePiezoActCSP-statsWindow):1+(peakIdx_postCuePiezoActCSP+statsWindow),c);
end
[~, p.piezoAlign.CSP_postCue_trialAvg_Cspk, ~, stats.piezoAlign.CSP_postCue_trialAvg_Cspk] = ttest(mean(postCuePiezo_baseline_CSP,1,'omitnan'),mean(postCuePiezo_activity_CSP,1,'omitnan'));
  m.piezoAlign.CSP_postCue_trialAvg_Cspk = mean(mean(postCuePiezo_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoAlign.CSP_postCue_trialAvg_Cspk = std(mean(postCuePiezo_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(postCuePiezo_activity_CSP,2)).*(1000./frameRateHz);
for c = 1:iDend
postCuePiezo_baseline_CSM(:,c) = CSMevents_acrossTrials(baselineRange(1)-1+(peakIdx_baseCSM-statsWindow):baselineRange(1)-1+(peakIdx_baseCSM+statsWindow),c);
postCuePiezo_activity_CSM(:,c) = postCue_piezoAlignEventsCSM_acrossTrials(prePiezoRange(1)-1+(peakIdx_postCuePiezoActCSM-statsWindow):prePiezoRange(1)-1+(peakIdx_postCuePiezoActCSM+statsWindow),c);
end
[~, p.piezoAlign.CSM_postCue_trialAvg_Cspk, ~, stats.piezoAlign.CSM_postCue_trialAvg_Cspk] = ttest(mean(postCuePiezo_baseline_CSM,1,'omitnan'),mean(postCuePiezo_activity_CSM,1,'omitnan'));
  m.piezoAlign.CSM_postCue_trialAvg_Cspk = mean(mean(postCuePiezo_activity_CSM,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoAlign.CSM_postCue_trialAvg_Cspk = std(mean(postCuePiezo_activity_CSM,1,'omitnan'),[],2,'omitnan')./sqrt(size(postCuePiezo_activity_CSM,2)).*(1000./frameRateHz);

   save(fullfile([output_fn 'piezoAlign.mat']), 'firstPiezo_baseline_CSP','firstPiezo_activity_CSP','postCuePiezo_activity_CSP','postRewPiezo_activity_CSP','firstPiezo_baseline_CSM','firstPiezo_activity_CSM','postCuePiezo_activity_CSM','postRewPiezo_activity_CSM','peakIdx_firstPiezoActCSP','peakIdx_postCuePiezoActCSP','peakIdx_postRewPiezoActCSP','peakIdx_baseCSP','peakIdx_firstPiezoActCSM','peakIdx_postCuePiezoActCSM','peakIdx_postRewPiezoActCSM','peakIdx_baseCSM');

%% plot     | statistical window for first-piezo avg vs. baseline avg
avg_firstPiezoPane(1,1) = round(mean(peakIdx_firstPiezoActCSP,2,'omitnan')); avg_firstPiezoPane(1,2) = round(mean(peakIdx_firstPiezoActCSM,2,'omitnan'));

tempFig=setFigure; hold on; 
subplot(2,1,1); hold on;
ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('Average Cspk rate preceeding first movement event in CS+ trials with statistical windows');
rectangle('Position',[tt(1,36-2+prePiezoRange(1)+avg_firstPiezoPane(1,1)-statsWindow)/1000,.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstPiezoAlignEvents(:,:,~gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(firstPiezoAlignEvents(:,:,~gBlock2),3,'omitnan'),[],2,'omitnan')./sqrt(size(firstPiezoAlignEvents,2)).*(1000./frameRateHz),'lineProps','k');
subplot(2,1,2); hold on;
ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('Average Cspk rate preceeding first movement event in CS- trials with statistical windows');
rectangle('Position',[tt(1,36-2+prePiezoRange(1)+avg_firstPiezoPane(1,2)-statsWindow)/1000,.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstPiezoAlignEvents(:,:,gBlock2),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(firstPiezoAlignEvents(:,:,gBlock2),3,'omitnan'),[],2,'omitnan')./sqrt(size(firstPiezoAlignEvents,2)).*(1000./frameRateHz),'lineProps','r');
saveas(tempFig, [output_fn 'stats\peizoWindow.pdf'],'pdf');

%% data org | most vs least motion 
mostMovementAlignedEvents=[]; leastMovementAlignedEvents=[]; mostTrialMovement=[]; leastTrialMovement=[];
mostCSPMovementAlignedEvents=[]; leastCSPMovementAlignedEvents=[]; mostCSPTrialMovement=[]; leastCSPTrialMovement=[];
mostCSMMovementAlignedEvents=[]; leastCSMMovementAlignedEvents=[]; mostCSMTrialMovement=[]; leastCSMTrialMovement=[];

for mouse = 1:size(expt,2)
    block2Trials{mouse}=[]; rewardTrials{mouse}=[];
    CSPTrialMovement_ind=[]; CSMTrialMovement_ind=[];
    peakTrialChange=[]; preakTrialChange_ind=[]; trialChange=[];
for c = 1:size(groupAlign.piezo{mouse},2)
for ii = 1:size(groupAlign.piezo{mouse},1)-1
    trialChange(ii,c) = groupAlign.piezo{mouse}(ii+1,c)-groupAlign.piezo{mouse}(ii,c);
end
[peakTrialChange(1,c),peakTrialChange_ind(1,c)] = max(trialChange(50:96,c));
end
[TrialMovement, TrialMovement_ind] = sort(peakTrialChange);

for c = 1:size(block2{1,mouse},2)
   if block2{1,mouse}(1,c)==1; block2Trials{mouse} = [block2Trials{mouse} c]; end
   if block2{1,mouse}(1,c)==0; rewardTrials{mouse} = [rewardTrials{mouse} c]; end
end
for c = 1:size(TrialMovement_ind,2)
    if ~isempty(find(rewardTrials{1,mouse} == TrialMovement_ind(1,c)))
CSPTrialMovement_ind = [CSPTrialMovement_ind TrialMovement_ind(1,c)];
    elseif ~isempty(find(block2Trials{1,mouse} == TrialMovement_ind(1,c)))
CSMTrialMovement_ind = [CSMTrialMovement_ind TrialMovement_ind(1,c)];
    end
end
mostCSPTrialMovement_ind{1,mouse} = []; leastCSPTrialMovement_ind{1,mouse} = [];
mostCSMTrialMovement_ind{1,mouse} = []; leastCSMTrialMovement_ind{1,mouse} = [];
mostTrialMovement_ind{1,mouse} = TrialMovement_ind(:,size(groupAlign.piezo{1,mouse},2)-(floor(size(groupAlign.piezo{1,mouse},2)./4)):end);
leastTrialMovement_ind{1,mouse} = TrialMovement_ind(:,1:round(size(groupAlign.piezo{1,mouse},2)./4));

mostCSPTrialMovement_ind{1,mouse} = CSPTrialMovement_ind(:,size(groupAlign.piezo{1,mouse}(:,~block2{1,mouse}),2)-(floor(size(groupAlign.piezo{1,mouse}(:,~block2{1,mouse}),2)./4)):size(groupAlign.piezo{1,mouse}(:,~block2{1,mouse}),2));
leastCSPTrialMovement_ind{1,mouse} = CSPTrialMovement_ind(:,1:ceil(size(groupAlign.piezo{1,mouse}(:,~block2{1,mouse}),2)./4));
mostCSPMovementAlignedEvents = [mostCSPMovementAlignedEvents nanmean(groupAlign.events{1,mouse}(:,:,mostCSPTrialMovement_ind{1,mouse}),3)];
leastCSPMovementAlignedEvents = [leastCSPMovementAlignedEvents nanmean(groupAlign.events{1,mouse}(:,:,leastCSPTrialMovement_ind{1,mouse}),3)];
mostCSMTrialMovement_ind{1,mouse} = CSMTrialMovement_ind(:,size(groupAlign.piezo{1,mouse}(:,logical(block2{1,mouse})),2)-(floor(size(groupAlign.piezo{1,mouse}(:,logical(block2{1,mouse})),2)./4)):size(groupAlign.piezo{1,mouse}(:,logical(block2{1,mouse})),2));
leastCSMTrialMovement_ind{1,mouse} = CSMTrialMovement_ind(:,1:ceil(size(groupAlign.piezo{1,mouse}(:,logical(block2{1,mouse})),2)./4));
mostCSMMovementAlignedEvents = [mostCSMMovementAlignedEvents nanmean(groupAlign.events{1,mouse}(:,:,mostCSMTrialMovement_ind{1,mouse}),3)];
leastCSMMovementAlignedEvents = [leastCSMMovementAlignedEvents nanmean(groupAlign.events{1,mouse}(:,:,leastCSMTrialMovement_ind{1,mouse}),3)];
end
%% ttest    | most vs least motion avg vs baseline avg 
peakStatsWindow=0;

for i = 1:iDend
    [cpeak_highMovement_postCue_CSP(i,1), cpeakIdx_highMovement_postCue_CSP(i,1)] = max(mostCSPMovementAlignedEvents(postCueRange,i),[],1);
    [cpeak_lowMovement_postCue_CSP(i,1), cpeakIdx_lowMovement_postCue_CSP(i,1)] = max(leastCSPMovementAlignedEvents(postCueRange,i),[],1);
    [cpeak_highMovement_postRew_CSP(i,1), cpeakIdx_highMovement_postRew_CSP(i,1)] = max(mostCSMMovementAlignedEvents(postRewRange,i),[],1);
    [cpeak_lowMovement_postRew_CSP(i,1), cpeakIdx_lowMovement_postRew_CSP(i,1)] = max(leastCSMMovementAlignedEvents(postRewRange,i),[],1);
end
[~, p.piezoRate.diffInPeak_postCue] = ttest((cpeakIdx_highMovement_postCue_CSP(~isnan(cpeakIdx_highMovement_postCue_CSP))),(cpeakIdx_lowMovement_postCue_CSP(~isnan(cpeakIdx_lowMovement_postCue_CSP))));
[~, p.piezoRate.diffInPeak_postRew] = ttest((cpeakIdx_highMovement_postRew_CSP(~isnan(cpeakIdx_highMovement_postRew_CSP))),(cpeakIdx_lowMovement_postRew_CSP(~isnan(cpeakIdx_lowMovement_postRew_CSP))));

  [peak_highMovement_postCue_CSP, peakIdx_highMovement_postCue_CSP] = max(mean(mostCSPMovementAlignedEvents(postCueRange,:),2,'omitnan'),[],1);
  [peak_lowMovement_postCue_CSP, peakIdx_lowMovement_postCue_CSP] = max(mean(leastCSPMovementAlignedEvents(postCueRange,:),2,'omitnan'),[],1);
  [peak_highMovement_postRew_CSP, peakIdx_highMovement_postRew_CSP] = max(mean(mostCSPMovementAlignedEvents(postRewRange,:),2,'omitnan'),[],1);
  [peak_lowMovement_postRew_CSP, peakIdx_lowMovement_postRew_CSP] = max(mean(leastCSPMovementAlignedEvents(postRewRange,:),2,'omitnan'),[],1);
  [peak_highMovement_postCue_CSM, peakIdx_highMovement_postCue_CSM] = max(mean(mostCSMMovementAlignedEvents(postCueRange,:),2,'omitnan'),[],1);
  [peak_lowMovement_postCue_CSM, peakIdx_lowMovement_postCue_CSM] = max(mean(leastCSMMovementAlignedEvents(postCueRange,:),2,'omitnan'),[],1);
  [peak_highMovement_postRew_CSM, peakIdx_highMovement_postRew_CSM] = max(mean(mostCSMMovementAlignedEvents(postRewRange,:),2,'omitnan'),[],1);
  [peak_lowMovement_postRew_CSM, peakIdx_lowMovement_postRew_CSM] = max(mean(leastCSMMovementAlignedEvents(postRewRange,:),2,'omitnan'),[],1);

  mostCSPTrialMovement=[]; leastCSPTrialMovement=[];
  mostCSMTrialMovement=[]; leastCSMTrialMovement=[];
  for mus = 1:size(expt,2)
      mostCSPTrialMovement = [mostCSPTrialMovement groupAlign.piezo{1,mus}(:,mostCSPTrialMovement_ind{1,mus})];
      leastCSPTrialMovement = [leastCSPTrialMovement groupAlign.piezo{1,mus}(:,leastCSPTrialMovement_ind{1,mus})];
      mostCSMTrialMovement = [mostCSMTrialMovement groupAlign.piezo{1,mus}(:,mostCSMTrialMovement_ind{1,mus})];
      leastCSMTrialMovement = [leastCSMTrialMovement groupAlign.piezo{1,mus}(:,leastCSMTrialMovement_ind{1,mus})];
  end

% lowVhigh postCue - CS+
for c = 1:sum(nsize)
low_postCuePiezo_activity_CSP(:,c) = leastCSPMovementAlignedEvents(postCueRange(1)-1+(peakIdx_lowMovement_postCue_CSP-peakStatsWindow):postCueRange(1)-1+(peakIdx_lowMovement_postCue_CSP+peakStatsWindow),c);  
high_postCuePiezo_activity_CSP(:,c) = mostCSPMovementAlignedEvents(postCueRange(1)-1+(peakIdx_highMovement_postCue_CSP-peakStatsWindow):postCueRange(1)-1+(peakIdx_highMovement_postCue_CSP+peakStatsWindow),c);   
end
[~, p.piezoRate.CSP_postCue_trialAvg_Cspk, ~, stats.piezoRate.CSP_postCue_trialAvg_Cspk] = ttest(mean(low_postCuePiezo_activity_CSP,1,'omitnan'),mean(high_postCuePiezo_activity_CSP,1,'omitnan'));
  n.piezoRate.CSP_postCue_lowRate_Cspk = length(cell2mat(leastCSPTrialMovement_ind));
  m.piezoRate.CSP_postCue_lowRate_Cspk = mean(mean(low_postCuePiezo_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoRate.CSP_postCue_lowRate_Cspk = std(mean(low_postCuePiezo_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(low_postCuePiezo_activity_CSP,2)).*(1000./frameRateHz);
  n.piezoRate.CSP_postCue_highRate_Cspk = length(cell2mat(mostCSPTrialMovement_ind));
  m.piezoRate.CSP_postCue_highRate_Cspk = mean(mean(high_postCuePiezo_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoRate.CSP_postCue_highRate_Cspk = std(mean(high_postCuePiezo_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(high_postCuePiezo_activity_CSP,2)).*(1000./frameRateHz);

  tempFig=setFigure; hold on; xline(0,'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  rectangle('Position',[tt(1,postCueRange(1)-1+peakIdx_lowMovement_postCue_CSP(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',[0.678,0.847,0.902],'FaceColor',[0.678,0.847,0.902]);
  rectangle('Position',[tt(1,postCueRange(1)-1+peakIdx_highMovement_postCue_CSP(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',[0.827,0.827,0.827],'FaceColor',[0.827,0.827,0.827]);
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(leastCSPMovementAlignedEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(leastCSPMovementAlignedEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(leastCSPMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(mostCSPMovementAlignedEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(mostCSPMovementAlignedEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(mostCSPMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L movement - postCue - CS+'); ylim([0 ymax]);
      saveas(tempFig, [output_fn 'stats\HL_piezoRate_postCue_CS+_PSTHs.pdf'],'pdf');

% lowVhigh postRew
for c = 1:sum(nsize)
low_postRewPiezo_activity_CSP(:,c) = leastCSPMovementAlignedEvents(postRewRange(1)-1+(peakIdx_lowMovement_postRew_CSP-peakStatsWindow):postRewRange(1)-1+(peakIdx_lowMovement_postRew_CSP+peakStatsWindow),c);  
high_postRewPiezo_activity_CSP(:,c) = mostCSPMovementAlignedEvents(postRewRange(1)-1+(peakIdx_highMovement_postRew_CSP-peakStatsWindow):postRewRange(1)-1+(peakIdx_highMovement_postRew_CSP+peakStatsWindow),c);   
end
[~, p.piezoRate.CSP_postRew_trialAvg_Cspk, ~, stats.piezoRate.CSP_postRew_trialAvg_Cspk] = ttest(mean(low_postRewPiezo_activity_CSP,1,'omitnan'),mean(high_postRewPiezo_activity_CSP,1,'omitnan'));
  n.piezoRate.CSP_postRew_lowRate_Cspk = length(cell2mat(leastCSPTrialMovement_ind));
  m.piezoRate.CSP_postRew_lowRate_Cspk = mean(mean(low_postRewPiezo_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoRate.CSP_postRew_lowRate_Cspk = std(mean(low_postRewPiezo_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(low_postRewPiezo_activity_CSP,2)).*(1000./frameRateHz);
  n.piezoRate.CSP_postRew_highRate_Cspk = length(cell2mat(mostCSPTrialMovement_ind));
  m.piezoRate.CSP_postRew_highRate_Cspk = mean(mean(high_postRewPiezo_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoRate.CSP_postRew_highRate_Cspk = std(mean(high_postRewPiezo_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(high_postRewPiezo_activity_CSP,2)).*(1000./frameRateHz);

tempFig=setFigure; hold on; xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  shadedErrorBar_CV(tt(:,postRewRange)/1000,mean(leastCSPMovementAlignedEvents(postRewRange,:),2,'omitnan').*(1000./frameRateHz), std(leastCSPMovementAlignedEvents(postRewRange,:),[],2,'omitnan')./sqrt(length(leastCSPMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,postRewRange)/1000,mean(mostCSPMovementAlignedEvents(postRewRange,:),2,'omitnan').*(1000./frameRateHz), std(mostCSPMovementAlignedEvents(postRewRange,:),[],2,'omitnan')./sqrt(length(mostCSPMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L movement - postRew - CS+'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_piezoRate_postRew_CS+_PSTHs.pdf'],'pdf');

tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(leastCSPMovementAlignedEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(leastCSPMovementAlignedEvents(36:96,:),[],2,'omitnan')./sqrt(length(leastCSPMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(mostCSPMovementAlignedEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(mostCSPMovementAlignedEvents(36:96,:),[],2,'omitnan')./sqrt(length(mostCSPMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L movement - CS+'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_piezoRate_wholeTC_CS+_PSTHs.pdf'],'pdf');

tempFig=setFigure; subplot(2,1,1); hold on; xline(0,'k'); xline(.767, 'k'); yticks(0:.1:.5); ylabel('piezo voltage (mV)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(leastCSPTrialMovement(36:96,:),2,'omitnan'), std(leastCSPTrialMovement(36:96,:),[],2,'omitnan')./sqrt(length(leastCSPTrialMovement)),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(mostCSPTrialMovement(36:96,:),2,'omitnan'), std(mostCSPTrialMovement(36:96,:),[],2,'omitnan')./sqrt(length(mostCSPTrialMovement)),'lineProps','k');
  ylim([0 .5]);  title('PSTH of H vs. L movement rate - CS+'); 
  subplot(2,1,2); hold on; xline(0,'k'); xline(.767, 'k'); yticks(0:.1:.5); ylabel('piezo voltage (mV)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(leastCSMTrialMovement(36:96,:),2,'omitnan'), std(leastCSMTrialMovement(36:96,:),[],2,'omitnan')./sqrt(length(leastCSMTrialMovement)),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(mostCSMTrialMovement(36:96,:),2,'omitnan'), std(mostCSMTrialMovement(36:96,:),[],2,'omitnan')./sqrt(length(mostCSMTrialMovement)),'lineProps','k');
  ylim([0 .5]);  title('PSTH of H vs. L movement rate - CS-');  
        saveas(tempFig, [output_fn 'stats\HL_piezoRate.pdf'],'pdf');
        
% lowVhigh postCue - CS-
for c = 1:sum(nsize)
low_postCuePiezo_activity_CSM(:,c) = leastCSMMovementAlignedEvents(postCueRange(1)-1+(peakIdx_lowMovement_postCue_CSM-peakStatsWindow):postCueRange(1)-1+(peakIdx_lowMovement_postCue_CSM+peakStatsWindow),c);  
high_postCuePiezo_activity_CSM(:,c) = mostCSMMovementAlignedEvents(postCueRange(1)-1+(peakIdx_highMovement_postCue_CSM-peakStatsWindow):postCueRange(1)-1+(peakIdx_highMovement_postCue_CSM+peakStatsWindow),c);   
end
[~, p.piezoRate.CSM_postCue_trialAvg_Cspk, ~, stats.piezoRate.CSM_postCue_trialAvg_Cspk] = ttest(mean(low_postCuePiezo_activity_CSM,1,'omitnan'),mean(high_postCuePiezo_activity_CSM,1,'omitnan'));
  n.piezoRate.CSM_postCue_lowRate_Cspk = length(cell2mat(leastCSMTrialMovement_ind));
  m.piezoRate.CSM_postCue_lowRate_Cspk = mean(mean(low_postCuePiezo_activity_CSM,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoRate.CSM_postCue_lowRate_Cspk = std(mean(low_postCuePiezo_activity_CSM,1,'omitnan'),[],2,'omitnan')./sqrt(size(low_postCuePiezo_activity_CSM,2)).*(1000./frameRateHz);
  n.piezoRate.CSM_postCue_highRate_Cspk = length(cell2mat(mostCSMTrialMovement_ind));
  m.piezoRate.CSM_postCue_highRate_Cspk = mean(mean(high_postCuePiezo_activity_CSM,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoRate.CSM_postCue_highRate_Cspk = std(mean(high_postCuePiezo_activity_CSM,1,'omitnan'),[],2,'omitnan')./sqrt(size(high_postCuePiezo_activity_CSM,2)).*(1000./frameRateHz);
  
  tempFig=setFigure; hold on; xline(0,'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(leastCSMMovementAlignedEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(leastCSMMovementAlignedEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(leastCSMMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(mostCSMMovementAlignedEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(mostCSMMovementAlignedEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(mostCSMMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L movement - postCue - CS-'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_piezoRate_postCue_CS-_PSTHs.pdf'],'pdf');

    
% lowVhigh postRew
for c = 1:sum(nsize)
low_postRewPiezo_activity_CSM(:,c) = leastCSMMovementAlignedEvents(postRewRange(1)-1+(peakIdx_lowMovement_postRew_CSM-peakStatsWindow):postRewRange(1)-1+(peakIdx_lowMovement_postRew_CSM+peakStatsWindow),c);  
high_postRewPiezo_activity_CSM(:,c) = mostCSMMovementAlignedEvents(postRewRange(1)-1+(peakIdx_highMovement_postRew_CSM-peakStatsWindow):postRewRange(1)-1+(peakIdx_highMovement_postRew_CSM+peakStatsWindow),c);   
end
[~, p.piezoRate.CSM_postRew_trialAvg_Cspk, ~, stats.piezoRate.CSM_postRew_trialAvg_Cspk] = ttest(mean(low_postRewPiezo_activity_CSM,1,'omitnan'),mean(high_postRewPiezo_activity_CSM,1,'omitnan'));
  n.piezoRate.CSM_postRew_lowRate_Cspk = length(cell2mat(leastCSMTrialMovement_ind));
  m.piezoRate.CSM_postRew_lowRate_Cspk = mean(mean(low_postRewPiezo_activity_CSM,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoRate.CSM_postRew_lowRate_Cspk = std(mean(low_postRewPiezo_activity_CSM,1,'omitnan'),[],2,'omitnan')./sqrt(size(low_postRewPiezo_activity_CSM,2)).*(1000./frameRateHz);
  n.piezoRate.CSM_postRew_highRate_Cspk = length(cell2mat(mostCSMTrialMovement_ind));
  m.piezoRate.CSM_postRew_highRate_Cspk = mean(mean(high_postRewPiezo_activity_CSM,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.piezoRate.CSM_postRew_highRate_Cspk = std(mean(high_postRewPiezo_activity_CSM,1,'omitnan'),[],2,'omitnan')./sqrt(size(high_postRewPiezo_activity_CSM,2)).*(1000./frameRateHz);

  tempFig=setFigure; hold on; xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  shadedErrorBar_CV(tt(:,postRewRange)/1000,mean(leastCSMMovementAlignedEvents(postRewRange,:),2,'omitnan').*(1000./frameRateHz), std(leastCSMMovementAlignedEvents(postRewRange,:),[],2,'omitnan')./sqrt(length(leastCSMMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,postRewRange)/1000,mean(mostCSMMovementAlignedEvents(postRewRange,:),2,'omitnan').*(1000./frameRateHz), std(mostCSMMovementAlignedEvents(postRewRange,:),[],2,'omitnan')./sqrt(length(mostCSMMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L movement - postRew - CS-'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_piezoRate_postRew_CS-_PSTHs.pdf'],'pdf');
tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,2)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,2)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(leastCSMMovementAlignedEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(leastCSMMovementAlignedEvents(36:96,:),[],2,'omitnan')./sqrt(length(leastCSMMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(mostCSMMovementAlignedEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(mostCSMMovementAlignedEvents(36:96,:),[],2,'omitnan')./sqrt(length(mostCSMMovementAlignedEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L movement - CS-'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_piezoRate_wholeTC_CS-_PSTHs.pdf'],'pdf');

%% data org | high vs low lick rate 
preRew_lickAlignHz = sum(preRew_lickAlign(postLick_frames:3*postLick_frames-1,:),1)./length(prewin_frames+1:prewin_frames+rewDelay_frames)/frameRateHz;
postRew_lickAlignHz = sum(postRew_lickAlign(postLick_frames+1:3*postLick_frames-1,:),1)./length(prewin_frames+rewDelay_frames+1:prewin_frames+rewDelay_frames+rewDelay_frames)/frameRateHz;


ndends=0; ntrials=0; 
low_preRewHz_CSPEvents=[]; high_preRewHz_CSPEvents=[]; low_preRewHz_CSMEvents=[]; high_preRewHz_CSMEvents=[];
low_preRewHz_CSPLicks=[]; high_preRewHz_CSPLicks=[]; low_preRewHz_CSMLicks=[]; high_preRewHz_CSMLicks=[];
low_preRewHz_CSPLickStart = []; low_preRewHz_CSMLickStart = []; high_preRewHz_CSPLickStart = []; high_preRewHz_CSMLickStart = [];

low_postRewHz_CSPEvents=[]; high_postRewHz_CSPEvents=[]; low_postRewHz_CSMEvents=[]; high_postRewHz_CSMEvents=[];
low_postRewHz_CSPLicks=[]; high_postRewHz_CSPLicks=[]; low_postRewHz_CSMLicks=[]; high_postRewHz_CSMLicks=[];
low_postRewHz_CSPLickStart = []; low_postRewHz_CSMLickStart = []; high_postRewHz_CSPLickStart = []; high_postRewHz_CSMLickStart = [];

for mouse = 1:exptCount
    %nDend = (1:nsize(1,mouseIdx))+ndends;
    nTrial = (1:animalTrials(1,mouse))+ntrials; 

    temp_preRewHz_CSP = preRewCSP_lickBurstHz(1,(1:animalTrials(1,mouse))+ntrials);
    temp_preRewHz_CSM = preRewCSM_lickBurstHz(1,(1:animalTrials(1,mouse))+ntrials);
temp_postRewHz_CSP = postRewCSP_lickBurstHz(1,(1:animalTrials(1,mouse))+ntrials);
temp_postRewHz_CSM = postRewCSM_lickBurstHz(1,(1:animalTrials(1,mouse))+ntrials);

    nanIdx_preRewHz_CSP = ~isnan(temp_preRewHz_CSP);
    nanIdx_preRewHz_CSM = ~isnan(temp_preRewHz_CSM);
nanIdx_postRewHz_CSP = ~isnan(temp_postRewHz_CSP);
nanIdx_postRewHz_CSM = ~isnan(temp_postRewHz_CSM);

    [sort_preRewHz_CSP, sort_preRewHz_CSP_ind{1,mouse}] = sort(temp_preRewHz_CSP);
    [sort_preRewHz_CSM, sort_preRewHz_CSM_ind{1,mouse}] = sort(temp_preRewHz_CSM);
[sort_postRewHz_CSP, sort_postRewHz_CSP_ind{1,mouse}] = sort(temp_postRewHz_CSP);
[sort_postRewHz_CSM, sort_postRewHz_CSM_ind{1,mouse}] = sort(temp_postRewHz_CSM);

if seshType ~= 2
       
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
    high_preRewHz_CSMEvents = [high_preRewHz_CSMEvents nanmean(groupAlign.events{1,mouse}(:,:,sort_preRewHz_CSM_ind{1,mouse}(1,size(sort_preRewHz_CSM_ind{1,mouse},2)-floor(size(sort_preRewHz_CSM_ind{1,mouse},2)./4):size(sort_preRewHz_CSM_ind{1,mouse},2))),3)];
    low_preRewHz_CSMEvents = [low_preRewHz_CSMEvents nanmean(groupAlign.events{1,mouse}(:,:,sort_preRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSM_ind{1,mouse},2)./4))),3)];  
    high_preRewHz_CSMLicks = [high_preRewHz_CSMLicks groupAlign.licks{1,mouse}(:,sort_preRewHz_CSM_ind{1,mouse}(1,size(sort_preRewHz_CSM_ind{1,mouse},2)-floor(size(sort_preRewHz_CSM_ind{1,mouse},2)./4):size(sort_preRewHz_CSM_ind{1,mouse},2)))];
    low_preRewHz_CSMLicks = [low_preRewHz_CSMLicks groupAlign.licks{1,mouse}(:,sort_preRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSM_ind{1,mouse},2)./4)))];  
    high_preRewHz_CSMLickStart = [high_preRewHz_CSMLickStart groupAlign.lickStart{1,mouse}(:,sort_preRewHz_CSM_ind{1,mouse}(1,size(sort_preRewHz_CSM_ind{1,mouse},2)-floor(size(sort_preRewHz_CSM_ind{1,mouse},2)./4):size(sort_preRewHz_CSM_ind{1,mouse},2)))];
    low_preRewHz_CSMLickStart = [low_preRewHz_CSMLickStart groupAlign.lickStart{1,mouse}(:,sort_preRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSM_ind{1,mouse},2)./4)))];
if size(sort_postRewHz_CSM,2)./4 == ceil(size(sort_postRewHz_CSM,2)./4)
low_postRewHz_CSM{1,mouse} = sort_postRewHz_CSM(1,1:ceil(size(sort_postRewHz_CSM,2)./4)+1); 
low_postRewHz_CSM_ind{1,mouse} = sort_postRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSM,2)./4)+1); 
else
low_postRewHz_CSM{1,mouse} = sort_postRewHz_CSM(1,1:ceil(size(sort_postRewHz_CSM,2)./4)); 
low_postRewHz_CSM_ind{1,mouse} = sort_postRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSM,2)./4)); 
end
high_postRewHz_CSM{1,mouse} = sort_postRewHz_CSM(1,size(sort_postRewHz_CSM,2)-floor(size(sort_postRewHz_CSM,2)./4):size(sort_postRewHz_CSM,2));
high_postRewHz_CSM_ind{1,mouse} = sort_postRewHz_CSM_ind{1,mouse}(1,size(sort_postRewHz_CSM,2)-floor(size(sort_postRewHz_CSM,2)./4):size(sort_postRewHz_CSM,2));
high_postRewHz_CSMEvents = [high_postRewHz_CSMEvents nanmean(groupAlign.events{1,mouse}(:,:,sort_postRewHz_CSM_ind{1,mouse}(1,size(sort_postRewHz_CSM_ind{1,mouse},2)-floor(size(sort_postRewHz_CSM_ind{1,mouse},2)./4):size(sort_postRewHz_CSM_ind{1,mouse},2))),3)];
low_postRewHz_CSMEvents = [low_postRewHz_CSMEvents nanmean(groupAlign.events{1,mouse}(:,:,sort_postRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSM_ind{1,mouse},2)./4))),3)];  
high_postRewHz_CSMLicks = [high_postRewHz_CSMLicks groupAlign.licks{1,mouse}(:,sort_postRewHz_CSM_ind{1,mouse}(1,size(sort_postRewHz_CSM_ind{1,mouse},2)-floor(size(sort_postRewHz_CSM_ind{1,mouse},2)./4):size(sort_postRewHz_CSM_ind{1,mouse},2)))];
low_postRewHz_CSMLicks = [low_postRewHz_CSMLicks groupAlign.licks{1,mouse}(:,sort_postRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSM_ind{1,mouse},2)./4)))];  
high_postRewHz_CSMLickStart = [high_postRewHz_CSMLickStart groupAlign.lickStart{1,mouse}(:,sort_postRewHz_CSM_ind{1,mouse}(1,size(sort_postRewHz_CSM_ind{1,mouse},2)-floor(size(sort_postRewHz_CSM_ind{1,mouse},2)./4):size(sort_postRewHz_CSM_ind{1,mouse},2)))];
low_postRewHz_CSMLickStart = [low_postRewHz_CSMLickStart groupAlign.lickStart{1,mouse}(:,sort_postRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSM_ind{1,mouse},2)./4)))];

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
    high_preRewHz_CSPEvents = [high_preRewHz_CSPEvents nanmean(groupAlign.events{1,mouse}(:,:,sort_preRewHz_CSP_ind{1,mouse}(1,size(sort_preRewHz_CSP_ind{1,mouse},2)-floor(size(sort_preRewHz_CSP_ind{1,mouse},2)./4):size(sort_preRewHz_CSP_ind{1,mouse},2))),3)];
    low_preRewHz_CSPEvents = [low_preRewHz_CSPEvents nanmean(groupAlign.events{1,mouse}(:,:,sort_preRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSP_ind{1,mouse},2)./4))),3)];
    high_preRewHz_CSPLicks = [high_preRewHz_CSPLicks groupAlign.licks{1,mouse}(:,sort_preRewHz_CSP_ind{1,mouse}(1,size(sort_preRewHz_CSP_ind{1,mouse},2)-floor(size(sort_preRewHz_CSP_ind{1,mouse},2)./4):size(sort_preRewHz_CSP_ind{1,mouse},2)))];
    low_preRewHz_CSPLicks = [low_preRewHz_CSPLicks groupAlign.licks{1,mouse}(:,sort_preRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSP_ind{1,mouse},2)./4)))];
    high_preRewHz_CSPLickStart = [high_preRewHz_CSPLickStart groupAlign.lickStart{1,mouse}(:,sort_preRewHz_CSP_ind{1,mouse}(1,size(sort_preRewHz_CSP_ind{1,mouse},2)-floor(size(sort_preRewHz_CSP_ind{1,mouse},2)./4):size(sort_preRewHz_CSP_ind{1,mouse},2)))];
    low_preRewHz_CSPLickStart = [low_preRewHz_CSPLickStart groupAlign.lickStart{1,mouse}(:,sort_preRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSP_ind{1,mouse},2)./4)))];
if size(sort_postRewHz_CSP,2)./4 == ceil(size(sort_postRewHz_CSP,2)./4)
low_postRewHz_CSP{1,mouse} = sort_postRewHz_CSP(1,1:ceil(size(sort_postRewHz_CSP,2)./4)+1);
low_postRewHz_CSP_ind{1,mouse} = sort_postRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSP,2)./4)+1);
else
low_postRewHz_CSP{1,mouse} = sort_postRewHz_CSP(1,1:ceil(size(sort_postRewHz_CSP,2)./4));
low_postRewHz_CSP_ind{1,mouse} = sort_postRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSP,2)./4));
end
high_postRewHz_CSP{1,mouse} = sort_postRewHz_CSP(1,size(sort_postRewHz_CSP,2)-floor(size(sort_postRewHz_CSP,2)./4):size(sort_postRewHz_CSP,2));
high_postRewHz_CSP_ind{1,mouse} = sort_postRewHz_CSP_ind{1,mouse}(1,size(sort_postRewHz_CSP,2)-floor(size(sort_postRewHz_CSP,2)./4):size(sort_postRewHz_CSP,2));
high_postRewHz_CSPEvents = [high_postRewHz_CSPEvents nanmean(groupAlign.events{1,mouse}(:,:,sort_postRewHz_CSP_ind{1,mouse}(1,size(sort_postRewHz_CSP_ind{1,mouse},2)-floor(size(sort_postRewHz_CSP_ind{1,mouse},2)./4):size(sort_postRewHz_CSP_ind{1,mouse},2))),3)];
low_postRewHz_CSPEvents = [low_postRewHz_CSPEvents nanmean(groupAlign.events{1,mouse}(:,:,sort_postRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSP_ind{1,mouse},2)./4))),3)];
high_postRewHz_CSPLicks = [high_postRewHz_CSPLicks groupAlign.licks{1,mouse}(:,sort_postRewHz_CSP_ind{1,mouse}(1,size(sort_postRewHz_CSP_ind{1,mouse},2)-floor(size(sort_postRewHz_CSP_ind{1,mouse},2)./4):size(sort_postRewHz_CSP_ind{1,mouse},2)))];
low_postRewHz_CSPLicks = [low_postRewHz_CSPLicks groupAlign.licks{1,mouse}(:,sort_postRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSP_ind{1,mouse},2)./4)))];
high_postRewHz_CSPLickStart = [high_postRewHz_CSPLickStart groupAlign.lickStart{1,mouse}(:,sort_postRewHz_CSP_ind{1,mouse}(1,size(sort_postRewHz_CSP_ind{1,mouse},2)-floor(size(sort_postRewHz_CSP_ind{1,mouse},2)./4):size(sort_postRewHz_CSP_ind{1,mouse},2)))];
low_postRewHz_CSPLickStart = [low_postRewHz_CSPLickStart groupAlign.lickStart{1,mouse}(:,sort_postRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSP_ind{1,mouse},2)./4)))];

elseif seshType==2    
    if sum(nanIdx_preRewHz_CSM)>0
          low_preRewHz_CSP{1,mouse} = []; low_preRewHz_CSP_ind{1,mouse} = []; low_postRewHz_CSP{1,mouse} = []; low_postRewHz_CSP_ind{1,mouse} = [];
          high_preRewHz_CSP{1,mouse} = []; high_preRewHz_CSP_ind{1,mouse} = []; high_postRewHz_CSP{1,mouse} = []; high_postRewHz_CSP_ind{1,mouse} = [];
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
    high_preRewHz_CSMEvents = [high_preRewHz_CSMEvents nanmean(groupAlign.events{1,mouse}(:,:,sort_preRewHz_CSM_ind{1,mouse}(1,size(sort_preRewHz_CSM_ind{1,mouse},2)-floor(size(sort_preRewHz_CSM_ind{1,mouse},2)./4):size(sort_preRewHz_CSM_ind{1,mouse},2))),3)];
    low_preRewHz_CSMEvents = [low_preRewHz_CSMEvents nanmean(groupAlign.events{1,mouse}(:,:,sort_preRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSM_ind{1,mouse},2)./4))),3)];  
    high_preRewHz_CSMLicks = [high_preRewHz_CSMLicks groupAlign.licks{1,mouse}(:,sort_preRewHz_CSM_ind{1,mouse}(1,size(sort_preRewHz_CSM_ind{1,mouse},2)-floor(size(sort_preRewHz_CSM_ind{1,mouse},2)./4):size(sort_preRewHz_CSM_ind{1,mouse},2)))];
    low_preRewHz_CSMLicks = [low_preRewHz_CSMLicks groupAlign.licks{1,mouse}(:,sort_preRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSM_ind{1,mouse},2)./4)))];  
    high_preRewHz_CSMLickStart = [high_preRewHz_CSMLickStart groupAlign.lickStart{1,mouse}(:,sort_preRewHz_CSM_ind{1,mouse}(1,size(sort_preRewHz_CSM_ind{1,mouse},2)-floor(size(sort_preRewHz_CSM_ind{1,mouse},2)./4):size(sort_preRewHz_CSM_ind{1,mouse},2)))];
    low_preRewHz_CSMLickStart = [low_preRewHz_CSMLickStart groupAlign.lickStart{1,mouse}(:,sort_preRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSM_ind{1,mouse},2)./4)))];
if size(sort_postRewHz_CSM,2)./4 == ceil(size(sort_postRewHz_CSM,2)./4)
low_postRewHz_CSM{1,mouse} = sort_postRewHz_CSM(1,1:ceil(size(sort_postRewHz_CSM,2)./4)+1); 
low_postRewHz_CSM_ind{1,mouse} = sort_postRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSM,2)./4)+1); 
else
low_postRewHz_CSM{1,mouse} = sort_postRewHz_CSM(1,1:ceil(size(sort_postRewHz_CSM,2)./4)); 
low_postRewHz_CSM_ind{1,mouse} = sort_postRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSM,2)./4)); 
end
high_postRewHz_CSM{1,mouse} = sort_postRewHz_CSM(1,size(sort_postRewHz_CSM,2)-floor(size(sort_postRewHz_CSM,2)./4):size(sort_postRewHz_CSM,2));
high_postRewHz_CSM_ind{1,mouse} = sort_postRewHz_CSM_ind{1,mouse}(1,size(sort_postRewHz_CSM,2)-floor(size(sort_postRewHz_CSM,2)./4):size(sort_postRewHz_CSM,2));
high_postRewHz_CSMEvents = [high_postRewHz_CSMEvents nanmean(groupAlign.events{1,mouse}(:,:,sort_postRewHz_CSM_ind{1,mouse}(1,size(sort_postRewHz_CSM_ind{1,mouse},2)-floor(size(sort_postRewHz_CSM_ind{1,mouse},2)./4):size(sort_postRewHz_CSM_ind{1,mouse},2))),3)];
low_postRewHz_CSMEvents = [low_postRewHz_CSMEvents nanmean(groupAlign.events{1,mouse}(:,:,sort_postRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSM_ind{1,mouse},2)./4))),3)];  
high_postRewHz_CSMLicks = [high_postRewHz_CSMLicks groupAlign.licks{1,mouse}(:,sort_postRewHz_CSM_ind{1,mouse}(1,size(sort_postRewHz_CSM_ind{1,mouse},2)-floor(size(sort_postRewHz_CSM_ind{1,mouse},2)./4):size(sort_postRewHz_CSM_ind{1,mouse},2)))];
low_postRewHz_CSMLicks = [low_postRewHz_CSMLicks groupAlign.licks{1,mouse}(:,sort_postRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSM_ind{1,mouse},2)./4)))];  
high_postRewHz_CSMLickStart = [high_postRewHz_CSMLickStart groupAlign.lickStart{1,mouse}(:,sort_postRewHz_CSM_ind{1,mouse}(1,size(sort_postRewHz_CSM_ind{1,mouse},2)-floor(size(sort_postRewHz_CSM_ind{1,mouse},2)./4):size(sort_postRewHz_CSM_ind{1,mouse},2)))];
low_postRewHz_CSMLickStart = [low_postRewHz_CSMLickStart groupAlign.lickStart{1,mouse}(:,sort_postRewHz_CSM_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSM_ind{1,mouse},2)./4)))];
    
    elseif sum(nanIdx_preRewHz_CSP)>0
        low_preRewHz_CSM{1,mouse} = []; low_preRewHz_CSM_ind{1,mouse} = [];
        low_postRewHz_CSM{1,mouse} = []; low_postRewHz_CSM_ind{1,mouse} = [];
        high_preRewHz_CSM{1,mouse} = []; high_preRewHz_CSM_ind{1,mouse} = [];
        high_postRewHz_CSM{1,mouse} = []; high_postRewHz_CSM_ind{1,mouse} = [];
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
    high_preRewHz_CSPEvents = [high_preRewHz_CSPEvents nanmean(groupAlign.events{1,mouse}(:,:,sort_preRewHz_CSP_ind{1,mouse}(1,size(sort_preRewHz_CSP_ind{1,mouse},2)-floor(size(sort_preRewHz_CSP_ind{1,mouse},2)./4):size(sort_preRewHz_CSP_ind{1,mouse},2))),3)];
    low_preRewHz_CSPEvents = [low_preRewHz_CSPEvents nanmean(groupAlign.events{1,mouse}(:,:,sort_preRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSP_ind{1,mouse},2)./4))),3)];
    high_preRewHz_CSPLicks = [high_preRewHz_CSPLicks groupAlign.licks{1,mouse}(:,sort_preRewHz_CSP_ind{1,mouse}(1,size(sort_preRewHz_CSP_ind{1,mouse},2)-floor(size(sort_preRewHz_CSP_ind{1,mouse},2)./4):size(sort_preRewHz_CSP_ind{1,mouse},2)))];
    low_preRewHz_CSPLicks = [low_preRewHz_CSPLicks groupAlign.licks{1,mouse}(:,sort_preRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSP_ind{1,mouse},2)./4)))];
    high_preRewHz_CSPLickStart = [high_preRewHz_CSPLickStart groupAlign.lickStart{1,mouse}(:,sort_preRewHz_CSP_ind{1,mouse}(1,size(sort_preRewHz_CSP_ind{1,mouse},2)-floor(size(sort_preRewHz_CSP_ind{1,mouse},2)./4):size(sort_preRewHz_CSP_ind{1,mouse},2)))];
    low_preRewHz_CSPLickStart = [low_preRewHz_CSPLickStart groupAlign.lickStart{1,mouse}(:,sort_preRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_preRewHz_CSP_ind{1,mouse},2)./4)))];
if size(sort_postRewHz_CSP,2)./4 == ceil(size(sort_postRewHz_CSP,2)./4)
low_postRewHz_CSP{1,mouse} = sort_postRewHz_CSP(1,1:ceil(size(sort_postRewHz_CSP,2)./4)+1);
low_postRewHz_CSP_ind{1,mouse} = sort_postRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSP,2)./4)+1);
else
low_postRewHz_CSP{1,mouse} = sort_postRewHz_CSP(1,1:ceil(size(sort_postRewHz_CSP,2)./4));
low_postRewHz_CSP_ind{1,mouse} = sort_postRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSP,2)./4));
end
high_postRewHz_CSP{1,mouse} = sort_postRewHz_CSP(1,size(sort_postRewHz_CSP,2)-floor(size(sort_postRewHz_CSP,2)./4):size(sort_postRewHz_CSP,2));
high_postRewHz_CSP_ind{1,mouse} = sort_postRewHz_CSP_ind{1,mouse}(1,size(sort_postRewHz_CSP,2)-floor(size(sort_postRewHz_CSP,2)./4):size(sort_postRewHz_CSP,2));
high_postRewHz_CSPEvents = [high_postRewHz_CSPEvents nanmean(groupAlign.events{1,mouse}(:,:,sort_postRewHz_CSP_ind{1,mouse}(1,size(sort_postRewHz_CSP_ind{1,mouse},2)-floor(size(sort_postRewHz_CSP_ind{1,mouse},2)./4):size(sort_postRewHz_CSP_ind{1,mouse},2))),3)];
low_postRewHz_CSPEvents = [low_postRewHz_CSPEvents nanmean(groupAlign.events{1,mouse}(:,:,sort_postRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSP_ind{1,mouse},2)./4))),3)];
high_postRewHz_CSPLicks = [high_postRewHz_CSPLicks groupAlign.licks{1,mouse}(:,sort_postRewHz_CSP_ind{1,mouse}(1,size(sort_postRewHz_CSP_ind{1,mouse},2)-floor(size(sort_postRewHz_CSP_ind{1,mouse},2)./4):size(sort_postRewHz_CSP_ind{1,mouse},2)))];
low_postRewHz_CSPLicks = [low_postRewHz_CSPLicks groupAlign.licks{1,mouse}(:,sort_postRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSP_ind{1,mouse},2)./4)))];
high_postRewHz_CSPLickStart = [high_postRewHz_CSPLickStart groupAlign.lickStart{1,mouse}(:,sort_postRewHz_CSP_ind{1,mouse}(1,size(sort_postRewHz_CSP_ind{1,mouse},2)-floor(size(sort_postRewHz_CSP_ind{1,mouse},2)./4):size(sort_postRewHz_CSP_ind{1,mouse},2)))];
low_postRewHz_CSPLickStart = [low_postRewHz_CSPLickStart groupAlign.lickStart{1,mouse}(:,sort_postRewHz_CSP_ind{1,mouse}(1,1:ceil(size(sort_postRewHz_CSP_ind{1,mouse},2)./4)))];
   end
end

    ndends = max(nDend);
    ntrials = max(nTrial);
end
%% ttest    | high vs low lick rate avg vs baseline avg 

for i = 1:iDend
    [cpeak_highLick_postCue_CSP(i,1), cpeakIdx_highLick_postCue_CSP(i,1)] = max(high_preRewHz_CSPEvents(postCueRange,i),[],1);
    [cpeak_lowLick_postCue_CSP(i,1), cpeakIdx_lowLick_postCue_CSP(i,1)] = max(low_preRewHz_CSPEvents(postCueRange,i),[],1);
    [cpeak_highLick_postRew_CSP(i,1), cpeakIdx_highLick_postRew_CSP(i,1)] = max(high_postRewHz_CSPEvents(postRewRange,i),[],1);
    [cpeak_lowLick_postRew_CSP(i,1), cpeakIdx_lowLick_postRew_CSP(i,1)] = max(low_postRewHz_CSPEvents(postRewRange,i),[],1);
     [cpeak_highLick_postCue_CSM(i,1), cpeakIdx_highLick_postCue_CSM(i,1)] = max(high_preRewHz_CSMEvents(postCueRange,i),[],1);
    [cpeak_lowLick_postCue_CSM(i,1), cpeakIdx_lowLick_postCue_CSM(i,1)] = max(low_preRewHz_CSMEvents(postCueRange,i),[],1);
    [cpeak_highLick_postRew_CSM(i,1), cpeakIdx_highLick_postRew_CSM(i,1)] = max(high_postRewHz_CSMEvents(postRewRange,i),[],1);
    [cpeak_lowLick_postRew_CSM(i,1), cpeakIdx_lowLick_postRew_CSM(i,1)] = max(low_postRewHz_CSMEvents(postRewRange,i),[],1);
end
[~, p.lickRate.diffInPeak_postCueCSP] = ttest((cpeakIdx_highLick_postCue_CSP(~isnan(cpeakIdx_highLick_postCue_CSP))),(cpeakIdx_lowLick_postCue_CSP(~isnan(cpeakIdx_lowLick_postCue_CSP))));
[~, p.lickRate.diffInPeak_postRewCSP] = ttest((cpeakIdx_highLick_postRew_CSP(~isnan(cpeakIdx_highLick_postRew_CSP))),(cpeakIdx_lowLick_postRew_CSP(~isnan(cpeakIdx_lowLick_postRew_CSP))));

[p.lickRate.rankSum_diffInPeak_postCueCSP, ~, stats.lickRate.rankSum_diffInPeak_postCueCSP] = ranksum(cpeakIdx_highLick_postCue_CSP(~isnan(cpeakIdx_highLick_postCue_CSP)),cpeakIdx_lowLick_postCue_CSP(~isnan(cpeakIdx_lowLick_postCue_CSP)));
[p.lickRate.rankSum_diffInPeak_postRewCSP,~,stats.lickRate.rankSum_diffInPeak_postRewCSP] = ranksum(cpeakIdx_highLick_postRew_CSP(~isnan(cpeakIdx_highLick_postRew_CSP)),cpeakIdx_lowLick_postRew_CSP(~isnan(cpeakIdx_lowLick_postRew_CSP)));

  [peak_highLick_postCue_CSP, peakIdx_highLick_postCue_CSP] = max(mean(high_preRewHz_CSPEvents(postCueRange,:),2,'omitnan'),[],1);
  [peak_lowLick_postCue_CSP, peakIdx_lowLick_postCue_CSP] = max(mean(low_preRewHz_CSPEvents(postCueRange,:),2,'omitnan'),[],1);
  [peak_highLick_postRew_CSP, peakIdx_highLick_postRew_CSP] = max(mean(high_postRewHz_CSPEvents(postRewRange,:),2,'omitnan'),[],1);
  [peak_lowLick_postRew_CSP, peakIdx_lowLick_postRew_CSP] = max(mean(low_postRewHz_CSPEvents(postRewRange,:),2,'omitnan'),[],1);
  [peak_highLick_postCue_CSM, peakIdx_highLick_postCue_CSM] = max(mean(high_preRewHz_CSMEvents(postCueRange,:),2,'omitnan'),[],1);
  [peak_lowLick_postCue_CSM, peakIdx_lowLick_postCue_CSM] = max(mean(low_preRewHz_CSMEvents(postCueRange,:),2,'omitnan'),[],1);
  [peak_highLick_postRew_CSM, peakIdx_highLick_postRew_CSM] = max(mean(high_postRewHz_CSMEvents(postRewRange,:),2,'omitnan'),[],1);
  [peak_lowLick_postRew_CSM, peakIdx_lowLick_postRew_CSM] = max(mean(low_postRewHz_CSMEvents(postRewRange,:),2,'omitnan'),[],1);

% lowVhigh postCue - CS+
for c = 1:sum(nsize)
low_postCueLick_activity_CSP(:,c) = low_preRewHz_CSPEvents(postCueRange(1)-1+(peakIdx_lowLick_postCue_CSP-peakStatsWindow):postCueRange(1)-1+(peakIdx_lowLick_postCue_CSP+peakStatsWindow),c);  
high_postCueLick_activity_CSP(:,c) = high_preRewHz_CSPEvents(postCueRange(1)-1+(peakIdx_highLick_postCue_CSP-peakStatsWindow):postCueRange(1)-1+(peakIdx_highLick_postCue_CSP+peakStatsWindow),c);   
end
[~, p.lickRate.CSP_postCue_trialAvg_Cspk, ~, stats.lickRate.CSP_postCue_trialAvg_Cspk] = ttest(mean(low_postCueLick_activity_CSP,1,'omitnan'),mean(high_postCueLick_activity_CSP,1,'omitnan'));
  n.lickRate.CSP_postCue_lowRate_Cspk = length(cell2mat(low_preRewHz_CSP_ind));
  m.lickRate.CSP_postCue_lowRate_Cspk = mean(mean(low_postCueLick_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickRate.CSP_postCue_lowRate_Cspk = std(mean(low_postCueLick_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(low_postCueLick_activity_CSP,2)).*(1000./frameRateHz);
  n.lickRate.CSP_postCue_highRate_Cspk = length(cell2mat(high_preRewHz_CSP_ind));
  m.lickRate.CSP_postCue_highRate_Cspk = mean(mean(high_postCueLick_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickRate.CSP_postCue_highRate_Cspk = std(mean(high_postCueLick_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(high_postCueLick_activity_CSP,2)).*(1000./frameRateHz);

  tempFig=setFigure; hold on; xline(0,'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  %rectangle('Position',[tt(1,postCueRange(1)-1+peakIdx_lowLick_postCue_CSP(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',[0.678,0.847,0.902],'FaceColor',[0.678,0.847,0.902]);
  %rectangle('Position',[tt(1,postCueRange(1)-1+peakIdx_highLick_postCue_CSP(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',[0.827,0.827,0.827],'FaceColor',[0.827,0.827,0.827]);
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(low_preRewHz_CSPEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(low_preRewHz_CSPEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(low_preRewHz_CSPEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(high_preRewHz_CSPEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(high_preRewHz_CSPEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(high_preRewHz_CSPEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L lick rate - postCue - CS+'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_lickRate_postCue_CS+_PSTHs.pdf'],'pdf');

% lowVhigh postRew
for c = 1:sum(nsize)
low_postRewLick_activity_CSP(:,c) = low_postRewHz_CSPEvents(postRewRange(1)-1+(peakIdx_lowLick_postRew_CSP-peakStatsWindow):postRewRange(1)-1+(peakIdx_lowLick_postRew_CSP),c);  
high_postRewLick_activity_CSP(:,c) = high_postRewHz_CSPEvents(postRewRange(1)-1+(peakIdx_highLick_postRew_CSP-peakStatsWindow):postRewRange(1)-1+(peakIdx_highLick_postRew_CSP),c);   
end
[~, p.lickRate.CSP_postRew_trialAvg_Cspk, ~, stats.lickRate.CSP_postRew_trialAvg_Cspk] = ttest(mean(low_postRewLick_activity_CSP,1,'omitnan'),mean(high_postRewLick_activity_CSP,1,'omitnan'));
  n.lickRate.CSP_postRew_lowRate_Cspk = length(cell2mat(low_postRewHz_CSP_ind));
  m.lickRate.CSP_postRew_lowRate_Cspk = mean(mean(low_postRewLick_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickRate.CSP_postRew_lowRate_Cspk = std(mean(low_postRewLick_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(low_postRewLick_activity_CSP,2)).*(1000./frameRateHz);
  n.lickRate.CSP_postRew_highRate_Cspk = length(cell2mat(high_postRewHz_CSP_ind));
  m.lickRate.CSP_postRew_highRate_Cspk = mean(mean(high_postRewLick_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickRate.CSP_postRew_highRate_Cspk = std(mean(high_postRewLick_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(high_postRewLick_activity_CSP,2)).*(1000./frameRateHz);

 tempFig=setFigure; hold on;xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postRewRange(1)-1+peakIdx_lowLick_postRew_CSP(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',[0.678,0.847,0.902],'FaceColor',[0.678,0.847,0.902]);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peakIdx_highLick_postRew_CSP(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',[0.827,0.827,0.827],'FaceColor',[0.827,0.827,0.827]);
  shadedErrorBar_CV(tt(:,postRewRange)/1000,mean(low_postRewHz_CSPEvents(postRewRange,:),2,'omitnan').*(1000./frameRateHz), std(low_postRewHz_CSPEvents(postRewRange,:),[],2,'omitnan')./sqrt(length(low_postRewHz_CSPEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,postRewRange)/1000,mean(high_postRewHz_CSPEvents(postRewRange,:),2,'omitnan').*(1000./frameRateHz), std(high_postRewHz_CSPEvents(postRewRange,:),[],2,'omitnan')./sqrt(length(high_postRewHz_CSPEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L lick rate - postRew - CS+'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_lickRate_postRew_CS+_PSTHs.pdf'],'pdf');
 tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax*1.5); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(low_postRewHz_CSPEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(low_postRewHz_CSPEvents(36:96,:),[],2,'omitnan')./sqrt(length(low_postRewHz_CSPEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(high_postRewHz_CSPEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(high_postRewHz_CSPEvents(36:96,:),[],2,'omitnan')./sqrt(length(high_postRewHz_CSPEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L lick rate - postRew - CS+'); ylim([0 ymax*1.5]);
        saveas(tempFig, [output_fn 'stats\HL_lickRate_wholeTC_rew_CS+_PSTHs.pdf'],'pdf');
 
        
  % lowVhigh postCue - CS+
for c = 1:sum(nsize)
low_postCueLick_activity_CSP(:,c) = low_preRewHz_CSPEvents(postCueRange(1)-1+(peakIdx_lowLick_postCue_CSP-peakStatsWindow):postCueRange(1)-1+(peakIdx_lowLick_postCue_CSP+peakStatsWindow),c);  
high_postCueLick_activity_CSP(:,c) = high_preRewHz_CSPEvents(postCueRange(1)-1+(peakIdx_highLick_postCue_CSP-peakStatsWindow):postCueRange(1)-1+(peakIdx_highLick_postCue_CSP+peakStatsWindow),c);   
end
[~, p.lickRate.CSP_postCue_trialAvg_Cspk, ~, stats.lickRate.CSP_postCue_trialAvg_Cspk] = ttest(mean(low_postCueLick_activity_CSP,1,'omitnan'),mean(high_postCueLick_activity_CSP,1,'omitnan'));
  n.lickRate.CSP_postCue_lowRate_Cspk = length(cell2mat(low_preRewHz_CSP_ind));
  m.lickRate.CSP_postCue_lowRate_Cspk = mean(mean(low_postCueLick_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickRate.CSP_postCue_lowRate_Cspk = std(mean(low_postCueLick_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(low_postCueLick_activity_CSP,2)).*(1000./frameRateHz);
  n.lickRate.CSP_postCue_highRate_Cspk = length(cell2mat(high_preRewHz_CSP_ind));
  m.lickRate.CSP_postCue_highRate_Cspk = mean(mean(high_postCueLick_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickRate.CSP_postCue_highRate_Cspk = std(mean(high_postCueLick_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(high_postCueLick_activity_CSP,2)).*(1000./frameRateHz);
  
  tempFig=setFigure; hold on; xline(0,'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(low_preRewHz_CSPEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(low_preRewHz_CSPEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(low_preRewHz_CSPEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(high_preRewHz_CSPEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(high_preRewHz_CSPEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(high_preRewHz_CSPEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L lick rate - postCue - CS-'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_lickRate_postCue_CS-_PSTHs.pdf'],'pdf');
      
        tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax*1.5); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(low_preRewHz_CSPEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(low_preRewHz_CSPEvents(36:96,:),[],2,'omitnan')./sqrt(length(low_preRewHz_CSPEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(high_preRewHz_CSPEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(high_preRewHz_CSPEvents(36:96,:),[],2,'omitnan')./sqrt(length(high_preRewHz_CSPEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L lick rate - postCue - CS+'); ylim([0 ymax*1.5]);
        saveas(tempFig, [output_fn 'stats\HL_lickRate_wholeTC_cue_CS+_PSTHs.pdf'],'pdf');

 tempFig=setFigure; subplot(2,1,1); hold on; xline(0,'k'); xline(.767, 'k'); yticks(0:.1:.5); ylabel('Lick rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(low_preRewHz_CSPLicks(36:96,:),2,'omitnan'), std(low_preRewHz_CSPLicks(36:96,:),[],2,'omitnan')./sqrt(length(low_preRewHz_CSPLicks)),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(high_preRewHz_CSPLicks(36:96,:),2,'omitnan'), std(high_preRewHz_CSPLicks(36:96,:),[],2,'omitnan')./sqrt(length(high_preRewHz_CSPLicks)),'lineProps','k');
  ylim([0 .5]);

  subplot(2,1,2); hold on; xline(0,'k'); xline(.767, 'k'); yticks(0:.1:.5); ylabel('Lick rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(low_postRewHz_CSPLicks(36:96,:),2,'omitnan'), std(low_postRewHz_CSPLicks(36:96,:),[],2,'omitnan')./sqrt(length(low_postRewHz_CSPLicks)),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(high_postRewHz_CSPLicks(36:96,:),2,'omitnan'), std(high_postRewHz_CSPLicks(36:96,:),[],2,'omitnan')./sqrt(length(low_postRewHz_CSPLicks)),'lineProps','k');
  ylim([0 .5]);
  title('PSTH of H vs. L lick rate - CS+'); 
        saveas(tempFig, [output_fn 'stats\HL_lickRate_CS+.pdf'],'pdf');

% lowVhigh postCue - CS-
for c = 1:sum(nsize)
low_postCueLick_activity_CSM(:,c) = low_preRewHz_CSMEvents(postCueRange(1)-1+(peakIdx_lowLick_postCue_CSM-peakStatsWindow):postCueRange(1)-1+(peakIdx_lowLick_postCue_CSM+peakStatsWindow),c);  
high_postCueLick_activity_CSM(:,c) = high_preRewHz_CSMEvents(postCueRange(1)-1+(peakIdx_highLick_postCue_CSM-peakStatsWindow):postCueRange(1)-1+(peakIdx_highLick_postCue_CSM+peakStatsWindow),c);   
end
[~, p.lickRate.CSM_postCue_trialAvg_Cspk, ~, stats.lickRate.CSM_postCue_trialAvg_Cspk] = ttest(mean(low_postCueLick_activity_CSM,1,'omitnan'),mean(high_postCueLick_activity_CSM,1,'omitnan'));
  n.lickRate.CSM_postCue_lowRate_Cspk = length(cell2mat(low_preRewHz_CSM_ind));
  m.lickRate.CSM_postCue_lowRate_Cspk = mean(mean(low_postCueLick_activity_CSM,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickRate.CSM_postCue_lowRate_Cspk = std(mean(low_postCueLick_activity_CSM,1,'omitnan'),[],2,'omitnan')./sqrt(size(low_postCueLick_activity_CSM,2)).*(1000./frameRateHz);
  n.lickRate.CSM_postCue_highRate_Cspk = length(cell2mat(high_preRewHz_CSM_ind));
  m.lickRate.CSM_postCue_highRate_Cspk = mean(mean(high_postCueLick_activity_CSM,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickRate.CSM_postCue_highRate_Cspk = std(mean(high_postCueLick_activity_CSM,1,'omitnan'),[],2,'omitnan')./sqrt(size(high_postCueLick_activity_CSM,2)).*(1000./frameRateHz);
  
  tempFig=setFigure; hold on; xline(0,'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(low_preRewHz_CSMEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(low_preRewHz_CSMEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(low_preRewHz_CSMEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(high_preRewHz_CSMEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(high_preRewHz_CSMEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(high_preRewHz_CSMEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L lick rate - postCue - CS-'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_lickRate_postCue_CS-_PSTHs.pdf'],'pdf');

% lowVhigh postRew
for c = 1:sum(nsize)
low_postRewLick_activity_CSM(:,c) = low_postRewHz_CSMEvents(postRewRange(1)-1+(peakIdx_lowLick_postRew_CSM-peakStatsWindow):postRewRange(1)-1+(peakIdx_lowLick_postRew_CSM+peakStatsWindow),c);  
high_postRewLick_activity_CSM(:,c) = high_postRewHz_CSMEvents(postRewRange(1)-1+(peakIdx_highLick_postRew_CSM-peakStatsWindow):postRewRange(1)-1+(peakIdx_highLick_postRew_CSM+peakStatsWindow),c);   
end
[~, p.lickRate.CSM_postRew_trialAvg_Cspk, ~, stats.lickRate.CSM_postRew_trialAvg_Cspk] = ttest(mean(low_postRewLick_activity_CSM,1,'omitnan'),mean(high_postRewLick_activity_CSM,1,'omitnan'));
  n.lickRate.CSM_postRew_lowRate_Cspk = length(cell2mat(low_postRewHz_CSM_ind));
  m.lickRate.CSM_postRew_lowRate_Cspk = mean(mean(low_postRewLick_activity_CSM,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickRate.CSM_postRew_lowRate_Cspk = std(mean(low_postRewLick_activity_CSM,1,'omitnan'),[],2,'omitnan')./sqrt(size(low_postRewLick_activity_CSM,2)).*(1000./frameRateHz);
  n.lickRate.CSM_postRew_highRate_Cspk = length(cell2mat(high_postRewHz_CSM_ind));
  m.lickRate.CSM_postRew_highRate_Cspk = mean(mean(high_postRewLick_activity_CSM,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickRate.CSM_postRew_highRate_Cspk = std(mean(high_postRewLick_activity_CSM,1,'omitnan'),[],2,'omitnan')./sqrt(size(high_postRewLick_activity_CSM,2)).*(1000./frameRateHz);

  tempFig=setFigure; hold on; xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  shadedErrorBar_CV(tt(:,postRewRange)/1000,mean(low_postRewHz_CSMEvents(postRewRange,:),2,'omitnan').*(1000./frameRateHz), std(low_postRewHz_CSMEvents(postRewRange,:),[],2,'omitnan')./sqrt(length(low_postRewHz_CSMEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,postRewRange)/1000,mean(high_postRewHz_CSMEvents(postRewRange,:),2,'omitnan').*(1000./frameRateHz), std(high_postRewHz_CSMEvents(postRewRange,:),[],2,'omitnan')./sqrt(length(high_postRewHz_CSMEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L lick rate - postRew - CS-'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_lickRate_postCue_CS-_PSTHs.pdf'],'pdf');
tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(low_preRewHz_CSMEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(low_preRewHz_CSMEvents(36:96,:),[],2,'omitnan')./sqrt(length(low_preRewHz_CSMEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(high_preRewHz_CSMEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(high_preRewHz_CSMEvents(36:96,:),[],2,'omitnan')./sqrt(length(high_preRewHz_CSMEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of H vs. L lick rate - postRew - CS-'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\HL_lickRate_wholeTC_CS-_PSTHs.pdf'],'pdf');
   
 tempFig=setFigure; subplot(2,1,1); hold on; xline(0,'k'); xline(.767, 'k'); yticks(0:.1:.5); ylabel('Lick rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(low_preRewHz_CSMLicks(36:96,:),2,'omitnan'), std(low_preRewHz_CSMLicks(36:96,:),[],2,'omitnan')./sqrt(length(low_preRewHz_CSMLicks)),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(high_preRewHz_CSMLicks(36:96,:),2,'omitnan'), std(high_preRewHz_CSMLicks(36:96,:),[],2,'omitnan')./sqrt(length(high_preRewHz_CSMLicks)),'lineProps','k');
  ylim([0 .5]);

  subplot(2,1,2); hold on; xline(0,'k'); xline(.767, 'k'); yticks(0:.1:.5); ylabel('Lick rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(low_postRewHz_CSMLicks(36:96,:),2,'omitnan'), std(low_postRewHz_CSMLicks(36:96,:),[],2,'omitnan')./sqrt(length(low_postRewHz_CSMLicks)),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(high_postRewHz_CSMLicks(36:96,:),2,'omitnan'), std(high_postRewHz_CSMLicks(36:96,:),[],2,'omitnan')./sqrt(length(low_postRewHz_CSMLicks)),'lineProps','k');
  ylim([0 .5]);
  title('PSTH of H vs. L lick rate - CS-'); 
        saveas(tempFig, [output_fn 'stats\HL_lickRate_CS-.pdf'],'pdf');
   
%% ttest    | earliest vs latest licks avg vs baseline avg 
n.lickTiming.uEarlyTiming = (mean(earlyRewBurstStart,2,'omitnan')-prewin_frames).*(1000./frameRateHz);
n.lickTiming.sEarlyTiming = (std(earlyRewBurstStart,[],2,'omitnan')./sqrt(sum(animalTrials))).*(1000./frameRateHz);
n.lickTiming.uLateTiming = (mean(lateRewBurstStart,2,'omitnan')-prewin_frames).*(1000./frameRateHz);
n.lickTiming.sLateTiming = (std(lateRewBurstStart,[],2,'omitnan')./sqrt(sum(animalTrials))).*(1000./frameRateHz);
[~, p.lickTiming.diffInTiming] = ttest(mean(earlyRewBurstStart,1,'omitnan'),mean(lateRewBurstStart,1,'omitnan'));

for i = 1:iDend
    [cpeak_earlyLick_first_CSP(i,1), cpeakIdx_earlyLick_first_CSP(i,1)] = max(earlyRewBurstEvents(preLickRange,i),[],1);
    [cpeak_lateLick_first_CSP(i,1), cpeakIdx_lateLick_first_CSP(i,1)] = max(lateRewBurstEvents(preLickRange,i),[],1);
    [cpeak_earlyLick_postCue_CSP(i,1), cpeakIdx_earlyLick_postCue_CSP(i,1)] = max(earlyRewBurstEvents(postCueRange,i),[],1);
    [cpeak_lateLick_postCue_CSP(i,1), cpeakIdx_lateLick_postCue_CSP(i,1)] = max(lateRewBurstEvents(postCueRange,i),[],1);
    [cpeak_earlyLick_postRew_CSP(i,1), cpeakIdx_earlyLick_postRew_CSP(i,1)] = max(earlyRewBurstEvents(postRewRange,i),[],1);
    [cpeak_lateLick_postRew_CSP(i,1), cpeakIdx_lateLick_postRew_CSP(i,1)] = max(lateRewBurstEvents(postRewRange,i),[],1);
end
[~, p.lickTiming.diffInPeak_first] = ttest((cpeakIdx_earlyLick_first_CSP(~isnan(cpeakIdx_earlyLick_first_CSP))),(cpeakIdx_lateLick_first_CSP(~isnan(cpeakIdx_lateLick_first_CSP))));
[~, p.lickTiming.diffInPeak_postCue] = ttest((cpeakIdx_earlyLick_postCue_CSP(~isnan(cpeakIdx_earlyLick_postCue_CSP))),(cpeakIdx_lateLick_postCue_CSP(~isnan(cpeakIdx_lateLick_postCue_CSP))));
[~, p.lickTiming.diffInPeak_postRew] = ttest((cpeakIdx_earlyLick_postRew_CSP(~isnan(cpeakIdx_earlyLick_postRew_CSP))),(cpeakIdx_lateLick_postRew_CSP(~isnan(cpeakIdx_lateLick_postRew_CSP))));

  [peak_earlyLick_first_CSP, peakIdx_earlyLick_first_CSP] = max(mean(earlyRewBurstEvents(preLickRange,:),2,'omitnan'),[],1);
  [peak_lateLick_first_CSP, peakIdx_lateLick_first_CSP] = max(mean(lateRewBurstEvents(preLickRange,:),2,'omitnan'),[],1);
  [peak_earlyLick_postCue_CSP, peakIdx_earlyLick_postCue_CSP] = max(mean(earlyRewBurstEvents(postCueRange,:),2,'omitnan'),[],1);
  [peak_lateLick_postCue_CSP, peakIdx_lateLick_postCue_CSP] = max(mean(lateRewBurstEvents(postCueRange,:),2,'omitnan'),[],1);
  [peak_earlyLick_postRew_CSP, peakIdx_earlyLick_postRew_CSP] = max(mean(earlyRewBurstEvents(postRewRange,:),2,'omitnan'),[],1);
  [peak_lateLick_postRew_CSP, peakIdx_lateLick_postRew_CSP] = max(mean(lateRewBurstEvents(postRewRange,:),2,'omitnan'),[],1);
 
 % earlyVlate first - CS+
for c = 1:sum(nsize)
early_firstLick_activity_CSP(:,c) = earlyRewBurstEvents(preLickRange(1)-1+(peakIdx_earlyLick_first_CSP-peakStatsWindow):preLickRange(1)-1+(peakIdx_earlyLick_first_CSP+peakStatsWindow),c);  
late_firstLick_activity_CSP(:,c) = lateRewBurstEvents(preLickRange(1)-1+(peakIdx_lateLick_first_CSP-peakStatsWindow):preLickRange(1)-1+(peakIdx_lateLick_first_CSP+peakStatsWindow),c);   
end
[~, p.lickTiming.CSP_first_trialAvg_Cspk, ~, stats.lickTiming.CSP_first_trialAvg_Cspk] = ttest(mean(early_firstLick_activity_CSP,1,'omitnan'),mean(late_firstLick_activity_CSP,1,'omitnan'));
  n.lickTiming.CSP_first_early_Cspk = length((earlyRewBurstStart));
  m.lickTiming.CSP_first_early_Cspk = mean(mean(early_firstLick_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickTiming.CSP_first_early_Cspk = std(mean(early_firstLick_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(early_firstLick_activity_CSP,2)).*(1000./frameRateHz);
      m.lickTiming.CSP_firstEarlyLickLat = mean((cpeakIdx_earlyLick_first_CSP(~isnan(cpeakIdx_earlyLick_first_CSP)))).*(1000./frameRateHz);
      s.lickTiming.CSP_firstEarlyLickLat = std((cpeakIdx_earlyLick_first_CSP(~isnan(cpeakIdx_earlyLick_first_CSP))))./sqrt(nIC).*(1000./frameRateHz);
  n.lickTiming.CSP_first_late_Cspk = length((lateRewBurstStart));
  m.lickTiming.CSP_first_late_Cspk = mean(mean(late_firstLick_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickTiming.CSP_first_late_Cspk = std(mean(late_firstLick_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(late_firstLick_activity_CSP,2)).*(1000./frameRateHz);
      m.lickTiming.CSP_firstLateLickLat = mean((cpeakIdx_lateLick_first_CSP(~isnan(cpeakIdx_lateLick_first_CSP)))).*(1000./frameRateHz);
      s.lickTiming.CSP_firstLateLickLat = std((cpeakIdx_lateLick_first_CSP(~isnan(cpeakIdx_lateLick_first_CSP))))./sqrt(nIC).*(1000./frameRateHz);
  
  tempFig=setFigure; hold on; xline(0, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,preLickRange(1)-1+peakIdx_lateLick_first_CSP(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',[0.678,0.847,0.902],'FaceColor',[0.678,0.847,0.902]);
  % rectangle('Position',[tt(1,preLickRange(1)-1+peakIdx_earlyLick_first_CSP(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',[0.827,0.827,0.827],'FaceColor',[0.827,0.827,0.827]);
  shadedErrorBar_CV(tt(:,preLickRange)/1000,mean(earlyRewBurstEvents(preLickRange,:),2,'omitnan').*(1000./frameRateHz), std(earlyRewBurstEvents(preLickRange,:),[],2,'omitnan')./sqrt(length(earlyRewBurstEvents)).*(1000./frameRateHz),'lineProps','k');
  shadedErrorBar_CV(tt(:,preLickRange)/1000,mean(lateRewBurstEvents(preLickRange,:),2,'omitnan').*(1000./frameRateHz), std(lateRewBurstEvents(preLickRange,:),[],2,'omitnan')./sqrt(length(lateRewBurstEvents)).*(1000./frameRateHz),'lineProps','b');
  title('PSTH of early vs. late licks - postCue - CS+'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\EL_lickTiming_first_CS+_PSTHs.pdf'],'pdf');
 
  
  tempFig=setFigure; hold on; xline(0, 'k'); xline(.767, 'k'); yticks(ymin:1:ymax*2); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,preLickRange(1)-1+peakIdx_lateLick_first_CSP(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',[0.678,0.847,0.902],'FaceColor',[0.678,0.847,0.902]);
  % rectangle('Position',[tt(1,preLickRange(1)-1+peakIdx_earlyLick_first_CSP(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',[0.827,0.827,0.827],'FaceColor',[0.827,0.827,0.827]);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(earlyRewBurstEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(earlyRewBurstEvents(36:96,:),[],2,'omitnan')./sqrt(length(earlyRewBurstEvents)).*(1000./frameRateHz),'lineProps','k');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(lateRewBurstEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(lateRewBurstEvents(36:96,:),[],2,'omitnan')./sqrt(length(lateRewBurstEvents)).*(1000./frameRateHz),'lineProps','b');
  
  plot(mean(((earlyRewBurstStart-prewin_frames)*(1000/frameRateHz)),2,'omitnan')/1000, 2, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5); 
  line([mean(((earlyRewBurstStart-prewin_frames)*(1000/frameRateHz)),2,'omitnan')/1000 - (std(((earlyRewBurstStart-prewin_frames)*(1000/frameRateHz)),[],2,'omitnan')/1000)/sqrt(size(expt,2)) ...
      , mean(((earlyRewBurstStart-prewin_frames)*(1000/frameRateHz)),2,'omitnan')/1000 + (std(((earlyRewBurstStart-prewin_frames)*(1000/frameRateHz)),[],2,'omitnan')/1000)/sqrt(size(expt,2))], [2,2], 'Color', 'k', 'LineWidth', 1); % Horizontal line for STD
  plot(mean(((lateRewBurstStart-prewin_frames)*(1000/frameRateHz)),2,'omitnan')/1000, 2, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5); 
  line([mean(((lateRewBurstStart-prewin_frames)*(1000/frameRateHz)),2,'omitnan')/1000 - (std(((lateRewBurstStart-prewin_frames)*(1000/frameRateHz)),[],2,'omitnan')/1000)/sqrt(size(expt,2)) ...
      , mean(((lateRewBurstStart-prewin_frames)*(1000/frameRateHz)),2,'omitnan')/1000 + (std(((lateRewBurstStart-prewin_frames)*(1000/frameRateHz)),[],2,'omitnan')/1000)/sqrt(size(expt,2))], [2,2], 'Color', 'b', 'LineWidth', 1); % Horizontal line for STD
  title('PSTH of early vs. late licks - postCue - CS+'); ylim([0 ymax*2]);
        saveas(tempFig, [output_fn 'stats\EL_lickTiming_full_CS+_PSTHs.pdf'],'pdf');

% earlyVlate postCue - CS+
for c = 1:sum(nsize)
early_postCueLick_activity_CSP(:,c) = earlyRewBurstEvents(postCueRange(1)-1+(peakIdx_earlyLick_postCue_CSP-peakStatsWindow):postCueRange(1)-1+(peakIdx_earlyLick_postCue_CSP+peakStatsWindow),c);  
late_postCueLick_activity_CSP(:,c) = lateRewBurstEvents(postCueRange(1)-1+(peakIdx_lateLick_postCue_CSP-peakStatsWindow):postCueRange(1)-1+(peakIdx_lateLick_postCue_CSP+peakStatsWindow),c);   
end
[~, p.lickTiming.CSP_postCue_trialAvg_Cspk, ~, stats.lickTiming.CSP_postCue_trialAvg_Cspk] = ttest(mean(early_postCueLick_activity_CSP,1,'omitnan'),mean(late_postCueLick_activity_CSP,1,'omitnan'));
  n.lickTiming.CSP_postCue_early_Cspk = length((earlyRewBurstStart));
  m.lickTiming.CSP_postCue_early_Cspk = mean(mean(early_postCueLick_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickTiming.CSP_postCue_early_Cspk = std(mean(early_postCueLick_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(early_postCueLick_activity_CSP,2)).*(1000./frameRateHz);
  n.lickTiming.CSP_postCue_late_Cspk = length((lateRewBurstStart));
  m.lickTiming.CSP_postCue_late_Cspk = mean(mean(late_postCueLick_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickTiming.CSP_postCue_late_Cspk = std(mean(late_postCueLick_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(late_postCueLick_activity_CSP,2)).*(1000./frameRateHz);

  tempFig=setFigure; hold on; xline(0, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peakIdx_lateLick_postCue_CSP(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',[0.678,0.847,0.902],'FaceColor',[0.678,0.847,0.902]);
  % rectangle('Position',[tt(1,postCueRange(1)-1+peakIdx_earlyLick_postCue_CSP(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',[0.827,0.827,0.827],'FaceColor',[0.827,0.827,0.827]);
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(earlyRewBurstEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(earlyRewBurstEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(earlyRewBurstEvents)).*(1000./frameRateHz),'lineProps','k');
  shadedErrorBar_CV(tt(:,postCueRange)/1000,mean(lateRewBurstEvents(postCueRange,:),2,'omitnan').*(1000./frameRateHz), std(lateRewBurstEvents(postCueRange,:),[],2,'omitnan')./sqrt(length(lateRewBurstEvents)).*(1000./frameRateHz),'lineProps','b');
  title('PSTH of early vs. late licks - postCue - CS+'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\EL_lickTiming_postCue_CS+_PSTHs.pdf'],'pdf');

% earlyVlate postRew
for c = 1:sum(nsize)
early_postRewLick_activity_CSP(:,c) = earlyRewBurstEvents(postRewRange(1)-1+(peakIdx_earlyLick_postRew_CSP-peakStatsWindow):postRewRange(1)-1+(peakIdx_earlyLick_postRew_CSP+peakStatsWindow),c);  
late_postRewLick_activity_CSP(:,c) = lateRewBurstEvents(postRewRange(1)-1+(peakIdx_lateLick_postRew_CSP-peakStatsWindow):postRewRange(1)-1+(peakIdx_lateLick_postRew_CSP+peakStatsWindow),c);   
end
[~, p.lickTiming.CSP_postRew_trialAvg_Cspk, ~, stats.lickTiming.CSP_postRew_trialAvg_Cspk] = ttest(mean(early_postRewLick_activity_CSP,1,'omitnan'),mean(late_postRewLick_activity_CSP,1,'omitnan'));
  n.lickTiming.CSP_postRew_early_Cspk = length((earlyRewBurstStart));
  m.lickTiming.CSP_postRew_early_Cspk = mean(mean(early_postRewLick_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickTiming.CSP_postRew_early_Cspk = std(mean(early_postRewLick_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(early_postRewLick_activity_CSP,2)).*(1000./frameRateHz);
  n.lickTiming.CSP_postRew_late_Cspk = length((lateRewBurstStart));
  m.lickTiming.CSP_postRew_late_Cspk = mean(mean(late_postRewLick_activity_CSP,1,'omitnan'),2,'omitnan').*(1000./frameRateHz);
  s.lickTiming.CSP_postRew_late_Cspk = std(mean(late_postRewLick_activity_CSP,1,'omitnan'),[],2,'omitnan')./sqrt(size(late_postRewLick_activity_CSP,2)).*(1000./frameRateHz);

  tempFig=setFigure; hold on; xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postRewRange(1)-1+peakIdx_lateLick_postRew_CSP(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',[0.678,0.847,0.902],'FaceColor',[0.678,0.847,0.902]);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peakIdx_earlyLick_postRew_CSP(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',[0.827,0.827,0.827],'FaceColor',[0.827,0.827,0.827]);
  shadedErrorBar_CV(tt(:,postRewRange)/1000,mean(earlyRewBurstEvents(postRewRange,:),2,'omitnan').*(1000./frameRateHz), std(earlyRewBurstEvents(postRewRange,:),[],2,'omitnan')./sqrt(length(earlyRewBurstEvents)).*(1000./frameRateHz),'lineProps','k');
  shadedErrorBar_CV(tt(:,postRewRange)/1000,mean(lateRewBurstEvents(postRewRange,:),2,'omitnan').*(1000./frameRateHz), std(lateRewBurstEvents(postRewRange,:),[],2,'omitnan')./sqrt(length(lateRewBurstEvents)).*(1000./frameRateHz),'lineProps','b');
  title('PSTH of early vs. late licks - postRew - CS+'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\EL_lickRate_postRew_CS+_PSTHs.pdf'],'pdf');

        tempFig=setFigure; hold on; xline(0,'k'); xline(.767, 'k'); yticks(ymin:1:ymax); ylabel('Cspk rate (Hz)'); xlabel('Time from cue onset (s)')
  % rectangle('Position',[tt(1,postCueRange(1)-1+peak_postCuePane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  % rectangle('Position',[tt(1,postRewRange(1)-1+peak_postRewPane(1,1)-peakStatsWindow)/1000,-1+.025,(1000./frameRateHz)*(peakStatsWindow*2)/1000,ymax*2+1],'EdgeColor',highlight,'FaceColor',highlight);
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(earlyRewBurstEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(earlyRewBurstEvents(36:96,:),[],2,'omitnan')./sqrt(length(earlyRewBurstEvents)).*(1000./frameRateHz),'lineProps','b');
  shadedErrorBar_CV(tt(:,36:96)/1000,mean(lateRewBurstEvents(36:96,:),2,'omitnan').*(1000./frameRateHz), std(lateRewBurstEvents(36:96,:),[],2,'omitnan')./sqrt(length(lateRewBurstEvents)).*(1000./frameRateHz),'lineProps','k');
  title('PSTH of early vs. late licks - postRew - CS+'); ylim([0 ymax]);
        saveas(tempFig, [output_fn 'stats\EL_lickTiming_postRew_wholeTC_CS+_PSTHs.pdf'],'pdf');

%% ttest    | lick-align avg vs cue-align avg 

[~, p.lickCue.CSP_postCue_trialAvg_Cspk, ~, stats.lickCue.CSP_postCue_trialAvg_Cspk] = ttest(mean(CSPevents_acrossTrials(postCueRange,:),1,'omitnan'),mean(postCue_lickAlignEventsCSP_acrossTrials(preLickRange,:),1,'omitnan'));
[~, p.lickCue.CSM_postCue_trialAvg_Cspk, ~, stats.lickCue.CSM_postCue_trialAvg_Cspk] = ttest(mean(CSMevents_acrossTrials(postCueRange,:),1,'omitnan'),mean(postCue_lickAlignEventsCSM_acrossTrials(preLickRange,:),1,'omitnan'));
[~, p.lickCue.CSP_postRew_trialAvg_Cspk, ~, stats.lickCue.CSP_postRew_trialAvg_Cspk] = ttest(mean(CSPevents_acrossTrials(postRewRange,:),1,'omitnan'),mean(postRew_lickAlignEventsCSP_acrossTrials(preLickRange,:),1,'omitnan'));
[~, p.lickCue.CSM_postRew_trialAvg_Cspk, ~, stats.lickCue.CSM_postRew_trialAvg_Cspk] = ttest(mean(CSMevents_acrossTrials(postRewRange,:),1,'omitnan'),mean(postRew_lickAlignEventsCSM_acrossTrials(preLickRange,:),1,'omitnan'));
%% ttest    | piezo-align avg vs cue-align avg 

[~, p.piezoCue.CSP_postCue_trialAvg_Cspk, ~, stats.piezoCue.CSP_postCue_trialAvg_Cspk] = ttest(mean(CSPevents_acrossTrials(postCueRange,:),1,'omitnan'),mean(postCue_piezoAlignEventsCSP_acrossTrials(prePiezoRange,:),1,'omitnan'));
[~, p.piezoCue.CSM_postCue_trialAvg_Cspk, ~, stats.piezoCue.CSM_postCue_trialAvg_Cspk] = ttest(mean(CSMevents_acrossTrials(postCueRange,:),1,'omitnan'),mean(postCue_piezoAlignEventsCSM_acrossTrials(prePiezoRange,:),1,'omitnan'));
[~, p.piezoCue.CSP_postRew_trialAvg_Cspk, ~, stats.piezoCue.CSP_postRew_trialAvg_Cspk] = ttest(mean(CSPevents_acrossTrials(postRewRange,:),1,'omitnan'),mean(postRew_piezoAlignEventsCSP_acrossTrials(prePiezoRange,:),1,'omitnan'));
[~, p.piezoCue.CSM_postRew_trialAvg_Cspk, ~, stats.piezoCue.CSM_postRew_trialAvg_Cspk] = ttest(mean(CSMevents_acrossTrials(postRewRange,:),1,'omitnan'),mean(postRew_piezoAlignEventsCSM_acrossTrials(prePiezoRange,:),1,'omitnan'));
%% ttest    | CS+ avg vs CS- avg 

[~, p.cueCue.CSMpostCue_trialAvg_Cspk] = ttest(mean(CSPevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSP-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSP+statsWindow),:),1,'omitnan'),mean(CSMevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSP-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSP+statsWindow),:),1,'omitnan'));
[~, p.cueCue.CSMpostRew_trialAvg_Cspk] = ttest(mean(CSPevents_acrossTrials(postRewRange(1)-1+(peakIdx_rewActCSP-statsWindow):postRewRange(1)-1+(peakIdx_rewActCSP+statsWindow),:),1,'omitnan'),mean(CSMevents_acrossTrials(postRewRange(1)-1+(peakIdx_rewActCSP-statsWindow):postRewRange(1)-1+(peakIdx_rewActCSP+statsWindow),:),1,'omitnan'));
[~, p.cueCue.CSPpostCue_trialAvg_Cspk] = ttest(mean(CSMevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSP-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSP+statsWindow),:),1,'omitnan'),mean(CSPevents_acrossTrials(postCueRange(1)-1+(peakIdx_cueActCSP-statsWindow):postCueRange(1)-1+(peakIdx_cueActCSP+statsWindow),:),1,'omitnan'));
[~, p.cueCue.CSPpostRew_trialAvg_Cspk] = ttest(mean(CSMevents_acrossTrials(postRewRange(1)-1+(peakIdx_rewActCSP-statsWindow):postRewRange(1)-1+(peakIdx_rewActCSP+statsWindow),:),1,'omitnan'),mean(CSPevents_acrossTrials(postRewRange(1)-1+(peakIdx_rewActCSP-statsWindow):postRewRange(1)-1+(peakIdx_rewActCSP+statsWindow),:),1,'omitnan'));
%% data org | n trials preceeding CS+ or CS- 

for mouse = 1:exptCount
tc.rewGroupAlign{1,mouse} = nanmean(groupAlign.tc{1,mouse}(:,:,logical(rew{1,mouse})),3);
tc.b2GroupAlign{1,mouse} = nanmean(groupAlign.tc{1,mouse}(:,:,logical(block2{1,mouse})),3);
events.rewGroupAlign{1,mouse} = nanmean(groupAlign.events{1,mouse}(:,:,logical(rew{1,mouse})),3);
events.b2GroupAlign{1,mouse} = nanmean(groupAlign.events{1,mouse}(:,:,logical(block2{1,mouse})),3);
piezo.rewGroupAlign{1,mouse} = groupAlign.piezo{1,mouse}(:,logical(rew{1,mouse}));
piezo.b2GroupAlign{1,mouse} = groupAlign.piezo{1,mouse}(:,logical(block2{1,mouse}));
lick.rewGroupAlign{1,mouse} = groupAlign.licks{1,mouse}(:,logical(rew{1,mouse}));
lick.b2GroupAlign{1,mouse} = groupAlign.licks{1,mouse}(:,logical(block2{1,mouse}));

tc.rewOneRewGroupAlign{1,mouse} = nanmean(groupAlign.tc{1,mouse}(:,:,logical(rewOneRew{1,mouse})),3);
events.rewOneRewGroupAlign{1,mouse} = nanmean(groupAlign.events{1,mouse}(:,:,logical(rewOneRew{1,mouse})),3);
piezo.rewOneRewGroupAlign{1,mouse} =  groupAlign.piezo{1,mouse}(:,logical(rewOneRew{1,mouse}));
lick.rewOneRewGroupAlign{1,mouse} =  groupAlign.licks{1,mouse}(:,logical(rewOneRew{1,mouse}));
tc.b2OneRewGroupAlign{1,mouse} = nanmean(groupAlign.tc{1,mouse}(:,:,logical(block2OneRew{1,mouse})),3);
events.b2OneRewGroupAlign{1,mouse} = nanmean(groupAlign.events{1,mouse}(:,:,logical(block2OneRew{1,mouse})),3);
piezo.b2OneRewGroupAlign{1,mouse} =  groupAlign.piezo{1,mouse}(:,logical(block2OneRew{1,mouse}));
lick.b2OneRewGroupAlign{1,mouse} =  groupAlign.licks{1,mouse}(:,logical(block2OneRew{1,mouse}));

tc.rewTwoRewGroupAlign{1,mouse} = nanmean(groupAlign.tc{1,mouse}(:,:,logical(rewTwoRew{1,mouse})),3);
events.rewTwoRewGroupAlign{1,mouse} = nanmean(groupAlign.events{1,mouse}(:,:,logical(rewTwoRew{1,mouse})),3);
piezo.rewTwoRewGroupAlign{1,mouse} =  groupAlign.piezo{1,mouse}(:,logical(rewTwoRew{1,mouse}));
lick.rewTwoRewGroupAlign{1,mouse} =  groupAlign.licks{1,mouse}(:,logical(rewTwoRew{1,mouse}));
tc.b2TwoRewGroupAlign{1,mouse} = nanmean(groupAlign.tc{1,mouse}(:,:,logical(block2TwoRew{1,mouse})),3);
events.b2TwoRewGroupAlign{1,mouse} = nanmean(groupAlign.events{1,mouse}(:,:,logical(block2TwoRew{1,mouse})),3);
piezo.b2TwoRewGroupAlign{1,mouse} =  groupAlign.piezo{1,mouse}(:,logical(block2TwoRew{1,mouse}));
lick.b2TwoRewGroupAlign{1,mouse} =  groupAlign.licks{1,mouse}(:,logical(block2TwoRew{1,mouse}));

tc.rewThreeRewGroupAlign{1,mouse} = nanmean(groupAlign.tc{1,mouse}(:,:,logical(rewThreeRew{1,mouse})),3);
events.rewThreeRewGroupAlign{1,mouse} = nanmean(groupAlign.events{1,mouse}(:,:,logical(rewThreeRew{1,mouse})),3);
piezo.rewThreeRewGroupAlign{1,mouse} =  groupAlign.piezo{1,mouse}(:,logical(rewThreeRew{1,mouse}));
lick.rewThreeRewGroupAlign{1,mouse} =  groupAlign.licks{1,mouse}(:,logical(rewThreeRew{1,mouse}));
tc.b2ThreeRewGroupAlign{1,mouse} = nanmean(groupAlign.tc{1,mouse}(:,:,logical(block2ThreeRew{1,mouse})),3);
events.b2ThreeRewGroupAlign{1,mouse} = nanmean(groupAlign.events{1,mouse}(:,:,logical(block2ThreeRew{1,mouse})),3);
piezo.b2ThreeRewGroupAlign{1,mouse} =  groupAlign.piezo{1,mouse}(:,logical(block2ThreeRew{1,mouse}));
lick.b2ThreeRewGroupAlign{1,mouse} =  groupAlign.licks{1,mouse}(:,logical(block2ThreeRew{1,mouse}));

tc.rewOneB2GroupAlign{1,mouse} = nanmean(groupAlign.tc{1,mouse}(:,:,logical(rewOneB2{1,mouse})),3);
events.rewOneB2GroupAlign{1,mouse} = nanmean(groupAlign.events{1,mouse}(:,:,logical(rewOneB2{1,mouse})),3);
piezo.rewOneB2GroupAlign{1,mouse} =  groupAlign.piezo{1,mouse}(:,logical(rewOneB2{1,mouse}));
lick.rewOneB2GroupAlign{1,mouse} =  groupAlign.licks{1,mouse}(:,logical(rewOneB2{1,mouse}));
tc.b2OneB2GroupAlign{1,mouse} = nanmean(groupAlign.tc{1,mouse}(:,:,logical(block2OneB2{1,mouse})),3);
events.b2OneB2GroupAlign{1,mouse} = nanmean(groupAlign.events{1,mouse}(:,:,logical(block2OneB2{1,mouse})),3);
piezo.b2OneB2GroupAlign{1,mouse} =  groupAlign.piezo{1,mouse}(:,logical(block2OneB2{1,mouse}));
lick.b2OneB2GroupAlign{1,mouse} =  groupAlign.licks{1,mouse}(:,logical(block2OneB2{1,mouse}));

tc.rewTwoB2GroupAlign{1,mouse} = nanmean(groupAlign.tc{1,mouse}(:,:,logical(rewTwoB2{1,mouse})),3);
events.rewTwoB2GroupAlign{1,mouse} = nanmean(groupAlign.events{1,mouse}(:,:,logical(rewTwoB2{1,mouse})),3);
piezo.rewTwoB2GroupAlign{1,mouse} =  groupAlign.piezo{1,mouse}(:,logical(rewTwoB2{1,mouse}));
lick.rewTwoB2GroupAlign{1,mouse} =  groupAlign.licks{1,mouse}(:,logical(rewTwoB2{1,mouse}));
tc.b2TwoB2GroupAlign{1,mouse} = nanmean(groupAlign.tc{1,mouse}(:,:,logical(block2TwoB2{1,mouse})),3);
events.b2TwoB2GroupAlign{1,mouse} = nanmean(groupAlign.events{1,mouse}(:,:,logical(block2TwoB2{1,mouse})),3);
piezo.b2TwoB2GroupAlign{1,mouse} =  groupAlign.piezo{1,mouse}(:,logical(block2TwoB2{1,mouse}));
lick.b2TwoB2GroupAlign{1,mouse} =  groupAlign.licks{1,mouse}(:,logical(block2TwoB2{1,mouse}));

tc.rewThreeB2GroupAlign{1,mouse} = nanmean(groupAlign.tc{1,mouse}(:,:,logical(rewThreeB2{1,mouse})),3);
events.rewThreeB2GroupAlign{1,mouse} = nanmean(groupAlign.events{1,mouse}(:,:,logical(rewThreeB2{1,mouse})),3);
piezo.rewThreeB2GroupAlign{1,mouse} =  groupAlign.piezo{1,mouse}(:,logical(rewThreeB2{1,mouse}));
lick.rewThreeB2GroupAlign{1,mouse} =  groupAlign.licks{1,mouse}(:,logical(rewThreeB2{1,mouse}));
tc.b2ThreeB2GroupAlign{1,mouse} = nanmean(groupAlign.tc{1,mouse}(:,:,logical(block2ThreeB2{1,mouse})),3);
events.b2ThreeB2GroupAlign{1,mouse} = nanmean(groupAlign.events{1,mouse}(:,:,logical(block2ThreeB2{1,mouse})),3);
piezo.b2ThreeB2GroupAlign{1,mouse} =  groupAlign.piezo{1,mouse}(:,logical(block2ThreeB2{1,mouse}));
lick.b2ThreeB2GroupAlign{1,mouse} =  groupAlign.licks{1,mouse}(:,logical(block2ThreeB2{1,mouse}));

tc.rewFourB2GroupAlign{1,mouse} = nanmean(groupAlign.tc{1,mouse}(:,:,logical(rewFourB2{1,mouse})),3);
events.rewFourB2GroupAlign{1,mouse} = nanmean(groupAlign.events{1,mouse}(:,:,logical(rewFourB2{1,mouse})),3);
piezo.rewFourB2GroupAlign{1,mouse} =  groupAlign.piezo{1,mouse}(:,logical(rewFourB2{1,mouse}));
lick.rewFourB2GroupAlign{1,mouse} =  groupAlign.licks{1,mouse}(:,logical(rewFourB2{1,mouse}));
tc.b2FourB2GroupAlign{1,mouse} = nanmean(groupAlign.tc{1,mouse}(:,:,logical(block2FourB2{1,mouse})),3);
events.b2FourB2GroupAlign{1,mouse} = nanmean(groupAlign.events{1,mouse}(:,:,logical(block2FourB2{1,mouse})),3);
piezo.b2FourB2GroupAlign{1,mouse} =  groupAlign.piezo{1,mouse}(:,logical(block2FourB2{1,mouse}));
lick.b2FourB2GroupAlign{1,mouse} =  groupAlign.licks{1,mouse}(:,logical(block2FourB2{1,mouse}));

tc.rewFourRewGroupAlign{1,mouse} = nanmean(groupAlign.tc{1,mouse}(:,:,logical(rewFourRew{1,mouse})),3);
events.rewFourRewGroupAlign{1,mouse} = nanmean(groupAlign.events{1,mouse}(:,:,logical(rewFourRew{1,mouse})),3);
piezo.rewFourRewGroupAlign{1,mouse} =  groupAlign.piezo{1,mouse}(:,logical(rewFourRew{1,mouse}));
lick.rewFourRewGroupAlign{1,mouse} =  groupAlign.licks{1,mouse}(:,logical(rewFourRew{1,mouse}));
tc.b2FourRewGroupAlign{1,mouse} = nanmean(groupAlign.tc{1,mouse}(:,:,logical(block2FourRew{1,mouse})),3);
events.b2FourRewGroupAlign{1,mouse} = nanmean(groupAlign.events{1,mouse}(:,:,logical(block2FourRew{1,mouse})),3);
piezo.b2FourRewGroupAlign{1,mouse} =  groupAlign.piezo{1,mouse}(:,logical(block2FourRew{1,mouse}));
lick.b2FourRewGroupAlign{1,mouse} =  groupAlign.licks{1,mouse}(:,logical(block2FourRew{1,mouse}));
end
% %% plot     | CS+\- preceeded by distinct CS+\- avg 

csPlusPlusSpectrum = {[.45 .75 0],[.15 .75 0],[.75 .95 0],[.5 1 .5]}; csPlusPlusIterations = [{'rewOneRewGroupAlign'},{'rewTwoRewGroupAlign'},{'rewThreeRewGroupAlign'},{'rewFourRewGroupAlign'}];
csPlusMinusSpectrum = {[0 0 1],[0 .45 .75],[0 .75 .75],[.3 .75 .95]}; csPlusMinusIterations = [{'rewOneB2GroupAlign'},{'rewTwoB2GroupAlign'},{'rewThreeB2GroupAlign'},{'rewFourB2GroupAlign'}];
csMinusPlusSpectrum = {[1 .25 .75],[1 .2 1],[.5 .1 .75],[.6 .2 1]}; csMinusPlusIterations = [{'b2OneRewGroupAlign'},{'b2TwoRewGroupAlign'},{'b2ThreeRewGroupAlign'},{'b2FourRewGroupAlign'}];
csMinusMinusSpectrum = {[1 0 0],[1 .5 .2],[.7 .4 .3],[1 .8 .7]}; csMinusMinusIterations = [{'b2OneB2GroupAlign'},{'b2TwoB2GroupAlign'},{'b2ThreeB2GroupAlign'},{'b2FourB2GroupAlign'}];


 tempFig=setFigure; hold on; %reward sequences - spikes
    subplot(2,2,1); hold on;
    ylim([0 ymax]);
    
    vline(0,'k');
    vline(.767,'k');
     hold on;
% CS+ preceeded by one CS+ trial   
     plotAlign=shadedErrorBar_CV(tt,nanmean(cell2mat(events.rewOneRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.rewOneRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csPlusPlusSpectrum{1}; plotAlign.patch.FaceColor = csPlusPlusSpectrum{1};
% CS+ preceeded by two CS+ trial   
     plotAlign=shadedErrorBar_CV(tt,nanmean(cell2mat(events.rewTwoRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.rewTwoRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csPlusPlusSpectrum{2}; plotAlign.patch.FaceColor = csPlusPlusSpectrum{2};
% CS+ preceeded by three CS+ trial   
     plotAlign=shadedErrorBar_CV(tt,nanmean(cell2mat(events.rewThreeRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.rewThreeRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csPlusPlusSpectrum{3}; plotAlign.patch.FaceColor = csPlusPlusSpectrum{3};
% CS+ preceeded by four CS+ trial   
     plotAlign=shadedErrorBar_CV(tt,nanmean(cell2mat(events.rewFourRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.rewFourRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csPlusPlusSpectrum{4}; plotAlign.patch.FaceColor = csPlusPlusSpectrum{4};
% set labels    
    xlabel('Time from cue');
    ylabel('Cspk rate (Hz)');
   title(['CS+ trials preceeded by one/two/three/four CS+ trial']);


 subplot(2,2,2); hold on;
    ylim([0 ymax]);
    
    vline(0,'k');
    vline(.767,'k');
     hold on;
% CS+ preceeded by one CS- trial   
     plotAlign=shadedErrorBar_CV(tt,nanmean(cell2mat(events.rewOneB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.rewOneB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csPlusMinusSpectrum{1}; plotAlign.patch.FaceColor = csPlusMinusSpectrum{1};
% CS+ preceeded by two CS- trial   
     plotAlign=shadedErrorBar_CV(tt,nanmean(cell2mat(events.rewTwoB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.rewTwoB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csPlusMinusSpectrum{2}; plotAlign.patch.FaceColor = csPlusMinusSpectrum{2};
% CS+ preceeded by three CS- trial   
     plotAlign=shadedErrorBar_CV(tt,nanmean(cell2mat(events.rewThreeB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.rewThreeB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csPlusMinusSpectrum{3}; plotAlign.patch.FaceColor = csPlusMinusSpectrum{3};
% CS+ preceeded by four CS- trial   
     plotAlign=shadedErrorBar_CV(tt,nanmean(cell2mat(events.rewFourB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.rewFourB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csPlusMinusSpectrum{4}; plotAlign.patch.FaceColor = csPlusMinusSpectrum{4};
% set labels    
    xlabel('Time from cue');
    ylabel('Cspk rate (Hz)');
   title(['CS+ trials preceeded by one/two/three/four CS- trial']);


    subplot(2,2,3); hold on;
    ylim([0 ymax]);
    
    vline(0,'k');
    vline(.767,'k');
     hold on;
% CS- preceeded by one CS+ trial   
     plotAlign=shadedErrorBar_CV(tt,nanmean(cell2mat(events.b2OneRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.b2OneRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csMinusPlusSpectrum{1}; plotAlign.patch.FaceColor = csMinusPlusSpectrum{1};
% CS- preceeded by two CS+ trial   
     plotAlign=shadedErrorBar_CV(tt,nanmean(cell2mat(events.b2TwoRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.b2TwoRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csMinusPlusSpectrum{2}; plotAlign.patch.FaceColor = csMinusPlusSpectrum{2};
% CS- preceeded by three CS+ trial   
     plotAlign=shadedErrorBar_CV(tt,nanmean(cell2mat(events.b2ThreeRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.b2ThreeRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csMinusPlusSpectrum{3}; plotAlign.patch.FaceColor = csMinusPlusSpectrum{3};
% CS- preceeded by four CS+ trial   
     plotAlign=shadedErrorBar_CV(tt,nanmean(cell2mat(events.b2FourRewGroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.b2FourRewGroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csMinusPlusSpectrum{4}; plotAlign.patch.FaceColor = csMinusPlusSpectrum{4};
% set labels    
    xlabel('Time from cue');
    ylabel('Cspk rate (Hz)');
   title(['CS- trials preceeded by one/two/three/four CS+ trial']);


 subplot(2,2,4); hold on;
    ylim([0 ymax]);
    
    vline(0,'k');
    vline(.767,'k');
     hold on;
% CS+ preceeded by one CS- trial   
     plotAlign=shadedErrorBar_CV(tt,nanmean(cell2mat(events.b2OneB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.b2OneB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csMinusMinusSpectrum{1}; plotAlign.patch.FaceColor = csMinusMinusSpectrum{1};
% CS+ preceeded by two CS- trial   
     plotAlign=shadedErrorBar_CV(tt,nanmean(cell2mat(events.b2TwoB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.b2TwoB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csMinusMinusSpectrum{2}; plotAlign.patch.FaceColor = csMinusMinusSpectrum{2};
% CS+ preceeded by three CS- trial   
     plotAlign=shadedErrorBar_CV(tt,nanmean(cell2mat(events.b2ThreeB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.b2ThreeB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csMinusMinusSpectrum{3}; plotAlign.patch.FaceColor = csMinusMinusSpectrum{3};
% CS+ preceeded by four CS- trial   
     plotAlign=shadedErrorBar_CV(tt,nanmean(cell2mat(events.b2FourB2GroupAlign),2).*(1000./frameRateHz), (nanstd(cell2mat(events.b2FourB2GroupAlign),[],2).*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csMinusMinusSpectrum{4}; plotAlign.patch.FaceColor = csMinusMinusSpectrum{4};
% set labels    
    xlabel('Time from cue');
    ylabel('Cspk rate (Hz)');
   title(['CS- trials preceeded by one/two/three/four CS- trial']);
%
    saveas(tempFig, [output_fn 'stats\CS_preceededByCS_PSTH.pdf'],'pdf');
%% plot     | CS+\- preceeded by distinct CS+\- maxima as bars 
maxWindow = postCueRange+5;

 tempFig=setFigure; hold on;
    subplot(2,2,1);hold on;
        ronr = cell2mat(events.rewOneRewGroupAlign);
    bar(1,max(nanmean(ronr(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csPlusPlusSpectrum{1});
        rtwr = cell2mat(events.rewTwoRewGroupAlign);
    bar(2,max(nanmean(rtwr(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csPlusPlusSpectrum{2});
        rthr = cell2mat(events.rewThreeRewGroupAlign);
    bar(3,max(nanmean(rthr(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csPlusPlusSpectrum{3});
        rfor = cell2mat(events.rewFourRewGroupAlign);
    bar(4,max(nanmean(rfor(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csPlusPlusSpectrum{4});
       ylim([0 ymax]);
       xlim([.5 4.5]);
       xlabel('# of preceeding cues');
       ylabel('Max Cspk rate (Hz)');
       title('CS+ preceeded by CS+');
   
     subplot(2,2,2);hold on;
        ronb = cell2mat(events.rewOneB2GroupAlign);
    bar(1,max(nanmean(ronb(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csPlusMinusSpectrum{1});
        rtwb = cell2mat(events.rewTwoB2GroupAlign);
    bar(2,max(nanmean(rtwb(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csPlusMinusSpectrum{2});
        rthb = cell2mat(events.rewThreeB2GroupAlign);
    bar(3,max(nanmean(rthb(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csPlusMinusSpectrum{3});
        rfob = cell2mat(events.rewFourB2GroupAlign);
    bar(4,max(nanmean(rfob(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csPlusMinusSpectrum{4});
       ylim([0 ymax]);
       xlim([.5 4.5]);
       xlabel('# of preceeding cues');
       ylabel('Max Cspk rate (Hz)');
       title('CS+ preceeded by CS-');
        
     subplot(2,2,3);hold on;
        bonr = cell2mat(events.b2OneRewGroupAlign);
    bar(1,max(nanmean(bonr(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csMinusPlusSpectrum{1});
        btwr = cell2mat(events.b2TwoRewGroupAlign);
    bar(2,max(nanmean(btwr(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csMinusPlusSpectrum{2});
        bthr = cell2mat(events.b2ThreeRewGroupAlign);
    bar(3,max(nanmean(bthr(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csMinusPlusSpectrum{3});
        bfor = cell2mat(events.b2FourRewGroupAlign);
    bar(4,max(nanmean(bfor(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csMinusPlusSpectrum{4});
       ylim([0 ymax]);
       xlim([.5 4.5]);
       xlabel('# of preceeding cues');
       ylabel('Max Cspk rate (Hz)');
       title('CS- preceeded by CS+');
        
     subplot(2,2,4);hold on;
        bonb = cell2mat(events.b2OneB2GroupAlign);
    bar(1,max(nanmean(bonb(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csMinusMinusSpectrum{1});
        btwb = cell2mat(events.b2TwoB2GroupAlign);
    bar(2,max(nanmean(btwb(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csMinusMinusSpectrum{2});
        bthb = cell2mat(events.b2ThreeB2GroupAlign);
    bar(3,max(nanmean(bthb(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csMinusMinusSpectrum{3});
        bfob = cell2mat(events.b2FourB2GroupAlign);
    bar(4,max(nanmean(bfob(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csMinusMinusSpectrum{4});
       ylim([0 ymax]);
       xlim([.5 4.5]);
       xlabel('# of preceeding cues');
       ylabel('Max Cspk rate (Hz)');
       title('CS- preceeded by CS-');
%   
    supertitle('Cue Window');
    saveas(tempFig, [output_fn 'stats\CS_preceededByCS_MaximaCue.pdf'],'pdf');


% plot maxima as bars rew
maxWindow = postRewRange+5;

 tempFig=setFigure; hold on;
    subplot(2,2,1);hold on;
        ronr = cell2mat(events.rewOneRewGroupAlign);
    bar(1,max(nanmean(ronr(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csPlusPlusSpectrum{1});
        rtwr = cell2mat(events.rewTwoRewGroupAlign);
    bar(2,max(nanmean(rtwr(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csPlusPlusSpectrum{2});
        rthr = cell2mat(events.rewThreeRewGroupAlign);
    bar(3,max(nanmean(rthr(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csPlusPlusSpectrum{3});
        rfor = cell2mat(events.rewFourRewGroupAlign);
    bar(4,max(nanmean(rfor(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csPlusPlusSpectrum{4});
       ylim([0 ymax]);
       xlim([.5 4.5]);
       xlabel('# of preceeding cues');
       ylabel('Max Cspk rate (Hz)');
       title('CS+ preceeded by CS+');
   
     subplot(2,2,2);hold on;
        ronb = cell2mat(events.rewOneB2GroupAlign);
    bar(1,max(nanmean(ronb(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csPlusMinusSpectrum{1});
        rtwb = cell2mat(events.rewTwoB2GroupAlign);
    bar(2,max(nanmean(rtwb(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csPlusMinusSpectrum{2});
        rthb = cell2mat(events.rewThreeB2GroupAlign);
    bar(3,max(nanmean(rthb(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csPlusMinusSpectrum{3});
        rfob = cell2mat(events.rewFourB2GroupAlign);
    bar(4,max(nanmean(rfob(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csPlusMinusSpectrum{4});
       ylim([0 ymax]);
       xlim([.5 4.5]);
       xlabel('# of preceeding cues');
       ylabel('Max Cspk rate (Hz)');
       title('CS+ preceeded by CS-');
        
     subplot(2,2,3);hold on;
        bonr = cell2mat(events.b2OneRewGroupAlign);
    bar(1,max(nanmean(bonr(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csMinusPlusSpectrum{1});
        btwr = cell2mat(events.b2TwoRewGroupAlign);
    bar(2,max(nanmean(btwr(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csMinusPlusSpectrum{2});
        bthr = cell2mat(events.b2ThreeRewGroupAlign);
    bar(3,max(nanmean(bthr(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csMinusPlusSpectrum{3});
        bfor = cell2mat(events.b2FourRewGroupAlign);
    bar(4,max(nanmean(bfor(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csMinusPlusSpectrum{4});
       ylim([0 ymax]);
       xlim([.5 4.5]);
       xlabel('# of preceeding cues');
       ylabel('Max Cspk rate (Hz)');
       title('CS- preceeded by CS+');
        
     subplot(2,2,4);hold on;
        bonb = cell2mat(events.b2OneB2GroupAlign);
    bar(1,max(nanmean(bonb(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csMinusMinusSpectrum{1});
        btwb = cell2mat(events.b2TwoB2GroupAlign);
    bar(2,max(nanmean(btwb(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csMinusMinusSpectrum{2});
        bthb = cell2mat(events.b2ThreeB2GroupAlign);
    bar(3,max(nanmean(bthb(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csMinusMinusSpectrum{3});
        bfob = cell2mat(events.b2FourB2GroupAlign);
    bar(4,max(nanmean(bfob(maxWindow,:),2).*(1000./frameRateHz)),'FaceColor',csMinusMinusSpectrum{4});
       ylim([0 ymax]);
       xlim([.5 4.5]);
       xlabel('# of preceeding cues');
       ylabel('Max Cspk rate (Hz)');
       title('CS- preceeded by CS-');
%   
    supertitle('Reward Window');
    saveas(tempFig, [output_fn 'stats\CS_preceededByCS_MaximaRew.pdf'],'pdf');  
%% plot     | CS+\- preceeded by distinct CS+\- aligned to lick avg 

tempFig=setFigure; hold on; 
    subplot(2,2,1); hold on;
ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('First lick CS+ trials following CS+ trials');
rectangle('Position',[tt(1,36-2+preLickRange(1)+avg_firstLickPane(1,1)-statsWindow)/1000,.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
% CS+ preceeded by one CS+ trial   
     plotAlign=shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,logical(cell2mat(rewOneRew))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,logical(cell2mat(rewOneRew))),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csPlusPlusSpectrum{1}; plotAlign.patch.FaceColor = csPlusPlusSpectrum{1};
% CS+ preceeded by two CS+ trial  
     plotAlign=shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,logical(cell2mat(rewTwoRew))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,logical(cell2mat(rewTwoRew))),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csPlusPlusSpectrum{2}; plotAlign.patch.FaceColor = csPlusPlusSpectrum{2};
% CS+ preceeded by three CS+ trial   
     plotAlign=shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,logical(cell2mat(rewThreeRew))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,logical(cell2mat(rewThreeRew))),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csPlusPlusSpectrum{3}; plotAlign.patch.FaceColor = csPlusPlusSpectrum{3};
% CS+ preceeded by four CS+ trial   
     plotAlign=shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,logical(cell2mat(rewFourRew))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,logical(cell2mat(rewFourRew))),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csPlusPlusSpectrum{4}; plotAlign.patch.FaceColor = csPlusPlusSpectrum{4};

    subplot(2,2,2); hold on;
ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('First lick CS+ trials following CS- trials');
rectangle('Position',[tt(1,36-2+preLickRange(1)+avg_firstLickPane(1,2)-statsWindow)/1000,.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
% CS+ preceeded by one CS+ trial   
     plotAlign=shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,logical(cell2mat(rewOneB2))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,logical(cell2mat(rewOneB2))),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csPlusMinusSpectrum{1}; plotAlign.patch.FaceColor = csPlusMinusSpectrum{1};
% CS+ preceeded by two CS+ trial  
     plotAlign=shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,logical(cell2mat(rewTwoB2))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,logical(cell2mat(rewTwoB2))),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csPlusMinusSpectrum{2}; plotAlign.patch.FaceColor = csPlusMinusSpectrum{2};
% CS+ preceeded by three CS+ trial   
     plotAlign=shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,logical(cell2mat(rewThreeB2))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,logical(cell2mat(rewThreeB2))),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csPlusMinusSpectrum{3}; plotAlign.patch.FaceColor = csPlusMinusSpectrum{3};
% CS+ preceeded by four CS+ trial   
     plotAlign=shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,logical(cell2mat(rewFourB2))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,logical(cell2mat(rewFourB2))),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csPlusMinusSpectrum{4}; plotAlign.patch.FaceColor = csPlusMinusSpectrum{4};
    
    subplot(2,2,3); hold on;
ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('First lick in CS- trials following CS+ trials');
rectangle('Position',[tt(1,36-2+preLickRange(1)+avg_firstLickPane(1,1)-statsWindow)/1000,.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
% CS+ preceeded by one CS+ trial   
     plotAlign=shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,logical(cell2mat(block2OneRew))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,logical(cell2mat(block2OneRew))),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csMinusPlusSpectrum{1}; plotAlign.patch.FaceColor = csMinusPlusSpectrum{1};
% CS+ preceeded by two CS+ trial  
     plotAlign=shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,logical(cell2mat(block2TwoRew))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,logical(cell2mat(block2TwoRew))),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csMinusPlusSpectrum{2}; plotAlign.patch.FaceColor = csMinusPlusSpectrum{2};
% CS+ preceeded by three CS+ trial   
     plotAlign=shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,logical(cell2mat(block2ThreeRew))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,logical(cell2mat(block2ThreeRew))),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csMinusPlusSpectrum{3}; plotAlign.patch.FaceColor = csMinusPlusSpectrum{3};
% CS+ preceeded by four CS+ trial   
     plotAlign=shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,logical(cell2mat(block2FourRew))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,logical(cell2mat(block2FourRew))),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csMinusPlusSpectrum{4}; plotAlign.patch.FaceColor = csMinusPlusSpectrum{4};

    subplot(2,2,4); hold on;
ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('First lick in CS- trials following CS- trials');
rectangle('Position',[tt(1,36-2+preLickRange(1)+avg_firstLickPane(1,2)-statsWindow)/1000,.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
% CS= preceeded by one CS- trial   
     plotAlign=shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,logical(cell2mat(block2OneB2))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,logical(cell2mat(block2OneB2))),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csMinusMinusSpectrum{1}; plotAlign.patch.FaceColor = csMinusMinusSpectrum{1};
% CS- preceeded by two CS- trial  
     plotAlign=shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,logical(cell2mat(block2TwoB2))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,logical(cell2mat(block2TwoB2))),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csMinusMinusSpectrum{2}; plotAlign.patch.FaceColor = csMinusMinusSpectrum{2};
% CS- preceeded by three CS- trial   
     plotAlign=shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,logical(cell2mat(block2ThreeB2))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,logical(cell2mat(block2ThreeB2))),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csMinusMinusSpectrum{3}; plotAlign.patch.FaceColor = csMinusMinusSpectrum{3};
% CS- preceeded by four CS- trial   
     plotAlign=shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,logical(cell2mat(block2FourB2))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz), (std(mean(firstLickAlignEvents(:,:,logical(cell2mat(block2FourB2))),3,'omitnan'),[],2,'omitnan').*(1000./frameRateHz))./sqrt(nIC));
     plotAlign.mainLine.Color = csMinusMinusSpectrum{4}; plotAlign.patch.FaceColor = csMinusMinusSpectrum{4};

    saveas(tempFig, [output_fn 'stats\CS_preceededByCS_lickAlignPSTH.pdf'],'pdf');
%% plot     | CS+\- preceeded by distinct CS+\- maxima as bars aligned to lick 
maxWindow = preLickRange;

 tempFig=setFigure; hold on;
    subplot(2,2,1);hold on;
    bar(1,max(mean(mean(firstLickAlignEvents(maxWindow,:,logical(cell2mat(rewOneRew))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz)),'FaceColor',csPlusPlusSpectrum{1});
    bar(2,max(mean(mean(firstLickAlignEvents(maxWindow,:,logical(cell2mat(rewTwoRew))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz)),'FaceColor',csPlusPlusSpectrum{2});
    bar(3,max(mean(mean(firstLickAlignEvents(maxWindow,:,logical(cell2mat(rewThreeRew))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz)),'FaceColor',csPlusPlusSpectrum{3});
    bar(4,max(mean(mean(firstLickAlignEvents(maxWindow,:,logical(cell2mat(rewFourRew))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz)),'FaceColor',csPlusPlusSpectrum{4});
       ylim([0 ymax]);
       xlim([.5 4.5]);
       xlabel('# of preceeding cues');
       ylabel('Max Cspk rate (Hz)');
       title('CS+ preceeded by CS+');
   
    subplot(2,2,2);hold on;
    bar(1,max(mean(mean(firstLickAlignEvents(maxWindow,:,logical(cell2mat(rewOneB2))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz)),'FaceColor',csPlusMinusSpectrum{1});
    bar(2,max(mean(mean(firstLickAlignEvents(maxWindow,:,logical(cell2mat(rewTwoB2))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz)),'FaceColor',csPlusMinusSpectrum{2});
    bar(3,max(mean(mean(firstLickAlignEvents(maxWindow,:,logical(cell2mat(rewThreeB2))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz)),'FaceColor',csPlusMinusSpectrum{3});
    bar(4,max(mean(mean(firstLickAlignEvents(maxWindow,:,logical(cell2mat(rewFourB2))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz)),'FaceColor',csPlusMinusSpectrum{4});
       ylim([0 ymax]);
       xlim([.5 4.5]);
       xlabel('# of preceeding cues');
       ylabel('Max Cspk rate (Hz)');
       title('CS+ preceeded by CS-');
        
    subplot(2,2,3);hold on;
    bar(1,max(mean(mean(firstLickAlignEvents(maxWindow,:,logical(cell2mat(block2OneRew))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz)),'FaceColor',csMinusPlusSpectrum{1});
    bar(2,max(mean(mean(firstLickAlignEvents(maxWindow,:,logical(cell2mat(block2TwoRew))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz)),'FaceColor',csMinusPlusSpectrum{2});
    bar(3,max(mean(mean(firstLickAlignEvents(maxWindow,:,logical(cell2mat(block2ThreeRew))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz)),'FaceColor',csMinusPlusSpectrum{3});
    bar(4,max(mean(mean(firstLickAlignEvents(maxWindow,:,logical(cell2mat(block2FourRew))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz)),'FaceColor',csMinusPlusSpectrum{4});
       ylim([0 ymax]);
       xlim([.5 4.5]);
       xlabel('# of preceeding cues');
       ylabel('Max Cspk rate (Hz)');
       title('CS- preceeded by CS+');
        
    subplot(2,2,4);hold on;
    bar(1,max(mean(mean(firstLickAlignEvents(maxWindow,:,logical(cell2mat(block2OneB2))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz)),'FaceColor',csMinusMinusSpectrum{1});
    bar(2,max(mean(mean(firstLickAlignEvents(maxWindow,:,logical(cell2mat(block2TwoB2))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz)),'FaceColor',csMinusMinusSpectrum{2});
    bar(3,max(mean(mean(firstLickAlignEvents(maxWindow,:,logical(cell2mat(block2ThreeB2))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz)),'FaceColor',csMinusMinusSpectrum{3});
    bar(4,max(mean(mean(firstLickAlignEvents(maxWindow,:,logical(cell2mat(block2FourB2))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz)),'FaceColor',csMinusMinusSpectrum{4});
       ylim([0 ymax]);
       xlim([.5 4.5]);
       xlabel('# of preceeding cues');
       ylabel('Max Cspk rate (Hz)');
       title('CS- preceeded by CS-');
%   
    supertitle('Cue Window');
    saveas(tempFig, [output_fn 'stats\CS_preceededByCS_MaximaLick.pdf'],'pdf');
%% data org | CS+\- preceeded by CS+\- 
for i = 1:size(expt,2)
rewRew{1,i} = sum([rewOneRew{1,i}; rewTwoRew{1,i}; rewThreeRew{1,i}; rewFourRew{1,i}]);
rewB2{1,i} = sum([rewOneB2{1,i}; rewTwoB2{1,i}; rewThreeB2{1,i}; rewFourB2{1,i}]);
b2Rew{1,i} = sum([block2OneRew{1,i}; block2TwoRew{1,i}; block2ThreeRew{1,i}; block2FourRew{1,i}]);
b2B2{1,i} = sum([block2OneB2{1,i}; block2TwoB2{1,i}; block2ThreeB2{1,i}; block2FourB2{1,i}]);
end
%% plot     | CS+\- preceeded by CS+\- aligned to lick avg

tempFig=setFigure; hold on; 
    subplot(2,2,1); hold on;
ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('First lick CS+ trials following CS+ trials');
rectangle('Position',[tt(1,36-2+preLickRange(1)+avg_firstLickPane(1,1)-statsWindow)/1000,.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
plotAlign=shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,logical(cell2mat(rewRew))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(firstLickAlignEvents(:,:,logical(cell2mat(rewRew))),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz),'lineProps','k');
plotAlign.mainLine.Color = csPlusPlusSpectrum{1}; plotAlign.patch.FaceColor = csPlusPlusSpectrum{1};

    subplot(2,2,2); hold on;
ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('First lick CS+ trials following CS- trials');
rectangle('Position',[tt(1,36-2+preLickRange(1)+avg_firstLickPane(1,2)-statsWindow)/1000,.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
plotAlign=shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,logical(cell2mat(rewB2))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(firstLickAlignEvents(:,:,logical(cell2mat(rewB2))),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz),'lineProps','r');
plotAlign.mainLine.Color = csPlusMinusSpectrum{1}; plotAlign.patch.FaceColor = csPlusMinusSpectrum{1};
    
    subplot(2,2,3); hold on;
ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('First lick in CS- trials following CS+ trials');
rectangle('Position',[tt(1,36-2+preLickRange(1)+avg_firstLickPane(1,1)-statsWindow)/1000,.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
plotAlign=shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,logical(cell2mat(b2Rew))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(firstLickAlignEvents(:,:,logical(cell2mat(b2Rew))),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz),'lineProps','k');
plotAlign.mainLine.Color = csMinusPlusSpectrum{1}; plotAlign.patch.FaceColor = csMinusPlusSpectrum{1};

    subplot(2,2,4); hold on;
ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('First lick in CS- trials following CS- trials');
rectangle('Position',[tt(1,36-2+preLickRange(1)+avg_firstLickPane(1,2)-statsWindow)/1000,.025,(1000./frameRateHz)*(statsWindow*2)/1000,ymax+1],'EdgeColor',highlight,'FaceColor',highlight);
plotAlign=shadedErrorBar_CV(tt(1,36:80)/1000,mean(mean(firstLickAlignEvents(:,:,logical(cell2mat(b2B2))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(firstLickAlignEvents(:,:,logical(cell2mat(b2B2))),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz),'lineProps','r');
plotAlign.mainLine.Color = csMinusMinusSpectrum{1}; plotAlign.patch.FaceColor = csMinusMinusSpectrum{1};

    saveas(tempFig, [output_fn 'stats\CS_preceededByCS_lickAlignPSTH_summary.pdf'],'pdf');
%% plot     | CS+\- preceeded by CS+\- avg
tempFig=setFigure; hold on; 
    subplot(2,2,1); hold on;
 ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('First CS+ trial following CS+ trials');
plotAlign=shadedErrorBar_CV(tt/1000,mean(mean(cueAlignEvents(:,:,logical(cell2mat(rewRew))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(cueAlignEvents(:,:,logical(cell2mat(rewRew))),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz),'lineProps','k');
plotAlign.mainLine.Color = csPlusPlusSpectrum{1}; plotAlign.patch.FaceColor = csPlusPlusSpectrum{1};

    subplot(2,2,2); hold on;
 ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('First CS+ trial following CS- trials');
plotAlign=shadedErrorBar_CV(tt/1000,mean(mean(cueAlignEvents(:,:,logical(cell2mat(rewB2))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(cueAlignEvents(:,:,logical(cell2mat(rewB2))),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz),'lineProps','r');
plotAlign.mainLine.Color = csPlusMinusSpectrum{1}; plotAlign.patch.FaceColor = csPlusMinusSpectrum{1};
    
    subplot(2,2,3); hold on;
 ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('First in CS- trial following CS+ trials');
plotAlign=shadedErrorBar_CV(tt/1000,mean(mean(cueAlignEvents(:,:,logical(cell2mat(b2Rew))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(cueAlignEvents(:,:,logical(cell2mat(b2Rew))),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz),'lineProps','k');
plotAlign.mainLine.Color = csMinusPlusSpectrum{1}; plotAlign.patch.FaceColor = csMinusPlusSpectrum{1};

    subplot(2,2,4); hold on;
 ylim([0 ymax]); xline(0,'k'); xline(.767,'k'); title('First in CS- trial following CS- trials');
plotAlign=shadedErrorBar_CV(tt/1000,mean(mean(cueAlignEvents(:,:,logical(cell2mat(b2B2))),3,'omitnan'),2,'omitnan').*(1000./frameRateHz),std(mean(cueAlignEvents(:,:,logical(cell2mat(b2B2))),3,'omitnan'),[],2,'omitnan')./sqrt(nIC).*(1000./frameRateHz),'lineProps','r');
plotAlign.mainLine.Color = csMinusMinusSpectrum{1}; plotAlign.patch.FaceColor = csMinusMinusSpectrum{1};

    saveas(tempFig, [output_fn 'stats\CS_preceededByCS_PSTH_summary.pdf'],'pdf');
%% save stats 
 save(fullfile([output_fn 'stats\simpleStatistics.mat']), 'p','m','s','stats','n','postCue*','postReward*','peak*','dendrite*','pos','neg','frameIdx','*Range');
%% Make tables for all t-tests

cueAlign_table = struct2table(p.cueAlign);
cueAlignWithLick_table = struct2table(p.cueAlign.withLick);
cueAlignWoutLick_table = struct2table(p.cueAlign.woutLick);
lickAlign_table = struct2table(p.lickAlign);
piezoAlign_table = struct2table(p.piezoAlign);
save(fullfile([output_fn 'stats\tableOfP-values.mat']), '*_table');
clearvars *_table

cueAlign_table = struct2table(m.cueAlign);
cueAlignWithLick_table = struct2table(m.cueAlign.withLick);
cueAlignWoutLick_table = struct2table(m.cueAlign.woutLick);
lickAlign_table = struct2table(m.lickAlign);
piezoAlign_table = struct2table(m.piezoAlign);
save(fullfile([output_fn 'stats\tableOfMeans.mat']), '*_table');
clearvars *_table

cueAlign_table = struct2table(s.cueAlign);
cueAlignWithLick_table = struct2table(s.cueAlign.withLick);
cueAlignWoutLick_table = struct2table(s.cueAlign.woutLick);
lickAlign_table = struct2table(s.lickAlign);
piezoAlign_table = struct2table(s.piezoAlign);
save(fullfile([output_fn 'stats\tableOfStdDev.mat']), '*_table');
clearvars *_table

cueAlign_table = struct2table(stats.cueAlign);
cueAlignWithLick_table = struct2table(stats.cueAlign.withLick);
cueAlignWoutLick_table = struct2table(stats.cueAlign.woutLick);
lickAlign_table = struct2table(stats.lickAlign);
piezoAlign_table = struct2table(stats.piezoAlign);
save(fullfile([output_fn 'stats\tableOfTStats.mat']), '*_table');
 %}