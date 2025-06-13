clear all; clear global;  close all
clc

ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories
eval(ds);

doMWCmPD = true; % generate the MW counter - photodiode counter plot or not

day_id = 27;

mouse = expt(day_id).mouse;
expDate = expt(day_id).date;
ExperimentFolder = expt(day_id).exptType;
dat = 'data-';
fn = fullfile(rc.achAnalysis,ExperimentFolder,mouse,expDate); %can make this flexible if folder structure is different
mkdir(fn)

runs = eval(['expt(day_id).' cell2mat(dataStructLabels) '_runs']);
times = eval(['expt(day_id).' cell2mat(dataStructLabels) '_time']);

fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\' dat mouse '-' expDate '-' times '.mat'];
load(cell2mat(fName));

%%

nOn = 30;
nOff = 60;

[stimOns stimOffs] = photoFrameFinder_Sanworks(info.frame);
nTrials = length(stimOns);
trialLength = nOn+nOff;
nFrames = trialLength*nTrials;
AllTrialsNFrames = nan(nTrials,1);
l1 = length(stimOns);
l2 = length(stimOffs);
if l1 > l2
    stimOffs(end+1) = 86400;
end
%MAXnFrames
% for itrial = 1:nTrials
%     if itrial == 1
%         startFrame = 0;
%     else
%         startFrame = stimOffs(itrial-1);
%     end
%     iOn_pd = stimOns(itrial);
%     nOff_pd = iOn_pd - startFrame; 
%     nOn_pd = stimOffs(itrial) - stimOns(itrial);
%     thisTrialNFrames = nOn_pd+nOff_pd;
%     AllTrialsNFrames(itrial) = thisTrialNFrames;
% end
% MAXnFrames = max(AllTrialsNFrames);
% data_tc = nan(MAXnFrames,nCells,nTrials);
% data_f_trial = nan(nCells,nTrials);
% %stimOffs(end+1) = nFrames;
% for itrial = 1:nTrials
%     if itrial == 1
%         firstFrame = 0;
%         startFrame = 1; % for indexing later
%     else
%         firstFrame = stimOffs(itrial-1);
%         startFrame = firstFrame;
%     end
%     iOn_pd = stimOns(itrial);
%     nOff_pd = iOn_pd - startFrame; 
%     nOn_pd = stimOffs(itrial) - stimOns(itrial);
%     thisTrialNFrames = nOn_pd+nOff_pd;
%     trialFrames = nan(MAXnFrames,nCells);
%     trialFrames(1:nOn_pd+nOff_pd,:) = npSub_tc(startFrame:stimOffs(itrial)-1,:);
%     %if ~isnan(stimOns(itrial)) & (stimOns(itrial)+nOn+nOff/2)<nFrames
%         data_tc(:,:,itrial) = trialFrames;
%     %end
% end
% 
% elements = zeros(size(data_tc,3),1);
% for m = 1:nTrials
%     elements(m,1) = sum(~isnan(data_tc(:,1,m)));
% end
% % plot 
% figure
% [counts centers] = hist(elements);
% histogram(elements);
% set(gca,'YScale','log')
% sgtitle('Log Trial Length Distribution')
% ylim([0.1 max(n)+max(n)*1/10])
% xlabel('nFrames')
% ylabel('log number of trials')
% 
% print(fullfile(fnout,'LogTrialLengthDistribution.pdf'),'-dpdf','-bestfit');
% 
% figure
% histogram(elements);
% sgtitle('Trial Length Distribution')
% ylim([0.1 max(n)+max(n)*1/10])
% xlabel('nFrames')
% ylabel('number of trials')
% 
% print(fullfile(fnout,'TrialLengthDistribution.pdf'),'-dpdf','-bestfit');

if doMWCmPD == true
    figure
    ptdcounter = stimOns;
    mwcounter = [60:90:86400];
    plot(ptdcounter-mwcounter);
    sgtitle('Counter Value Diff')
    xlabel('trial number')
    ylabel('Photodiode Counter - mWorks Counter')
    % print(fullfile(fnout,'ptd-mw.pdf'),'-dpdf','-bestfit');
else
    fprintf('No MWorks Counter - Photodiode Counter plot generated');
end