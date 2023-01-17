% to analyze individual day's data
clear all; clear global; close all

%identifying animal and run


mouse = 'WK24';
date = '220721';
time = char('1106');
ImgFolder = char('002');
RetImgFolder = char('001');

frame_rate = 30; %enter the frame rate, or I can edit this to enter the stimulus duration


%setting my paths
fn_base = 'Z:\home\Celine\Analysis\2p_analysis\';
fn = fullfile(fn_base,mouse,date,ImgFolder);
mkdir(fn);
cd(fn);

beh_prefix = strcat('Z:\Behavior\Data\data-');
beh_file = [beh_prefix mouse '-' date '-' time '.mat'];
load(beh_file); %load the mworks behavioral file

run_str = ['runs-' ImgFolder];
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];

% load data
load(fullfile([datemouserun '_TCs.mat']));
% loads 'mask_cell', 'mask_np'
load(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse,date,RetImgFolder, [date '_' mouse '_runs-' RetImgFolder '_lbub_fits.mat']))

%NEED TO LOAD MASKS TO GET THE INDEXES OF RED CELLS
ref = char('001');
load(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse, date, ref, [date '_' mouse '_runs-' ref '_mask_cell.mat']))
%% convert to trials
nOn = input.nScansOn;
nOff=input.nScansOff;
%change this to use a padded array, where I add zeros at the end. test=padarray(cellTCs_match{1},30,0,'post');
stimStart = (nOff/2)+1; %this indicates both the perdiod to trim off the start and the stim on period after trimming
stimEnd=stimStart+nOn-1;
nTrialsOG=length(input.tGratingDiameterDeg);
nTrials= nTrialsOG-1;

cellTCs_OG = npSub_tc;

%use the list of direction by trial to figure out how many trials there are
%currently the way I center the stim on period requires me to cut out
%one trial, hence the -1

nCells = size(npSub_tc,2);

npSub_tc = npSub_tc(stimStart:size(npSub_tc,1),:);
npSub_tc = padarray(npSub_tc,double(stimStart)-1,0,'post');
nFrames = size(npSub_tc,1);
data_trial = reshape(npSub_tc,[nOn+nOff nTrialsOG nCells]);
data_f= mean(data_trial(1: (nOff/2),:,:),1);
data_tc_trial = bsxfun(@rdivide,bsxfun(@minus,data_trial,data_f),data_f);
data_tc_trial = data_tc_trial(:,1:nTrials-1,:);
clear data_trial data_f
plot(squeeze(mean(data_tc_trial, 2)))
%% find the stimulus conditions
tCons = celleqel2mat_padded(input.tGratingContrast(1:nTrials-1)); %transforms cell array into matrix (1 x ntrials)
Cons = unique(tCons);
nCons = length(Cons);

tSize = celleqel2mat_padded(input.tGratingDiameterDeg(1:nTrials-1)); %transforms cell array into matrix (1 x ntrials)
Sizes = unique(tSize);
nSizes = length(Sizes);

% find preferred size and contrast, find cells that are responsive 

data_resp = zeros(nCells,nSizes,nCons,2);
h = zeros(nCells, nSizes,nCons);
p = zeros(nCells, nSizes,nCons);

resp_win = stimStart+2:stimStart+8;
base_win = (stimStart - nOff/2):stimStart-1;



for iSize = 1:nSizes
    ind_size = find(tSize == Sizes(iSize));
    for iCon = 1:nCons
        ind_con = find(tCons == Cons(iCon));
        ind = intersect(ind_size,ind_con); %for every orientation and then every contrast, find trials with that con/ori combination
        data_resp(:,iSize,iCon,1) = squeeze(mean(mean(data_tc_trial(resp_win,ind,:),1),2));
        data_resp(:,iSize,iCon,2) = squeeze(std(mean(data_tc_trial(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
        [h(:,iSize,iCon), p(:,iSize,iCon)] = ttest(mean(data_tc_trial(resp_win,ind,:),1), mean(data_tc_trial(base_win,ind,:),1),'dim',2,'tail','right','alpha',0.05./((nSizes*nCons)-1));
    end
end
% sizeOnly = squeeze(mean(data_resp(:,:,:,1),3));
% for iCell = 1:nCells
%    [~,prefSize(1,iCell)] = max(sizeOnly(iCell,:,:,1),[],2); 
% end

h_all = sum(sum(h,2),3);
resp_all=logical(h_all);
% resp_small_highCon = logical(h(:,1,4));

load([date '_' mouse '_' run_str '_centerCells.mat']); %load the centerCells matrix
% centeredResp=intersect(find(resp_all),centerCells);
% smallResp=intersect(centeredResp,find(prefSize<4));
% length(centeredResp)
goodFitResp=intersect(find(resp_all),goodfit_ind);
length(goodFitResp)
keepDists = cellDists(goodFitResp);

% clear sizeOnly
%% plot timecourses at different sizes and contrasts
meanTC_byCondition=nan(nOn+nOff,nSizes,nCons,length(goodFitResp));
%frame_rate=double(input.frameImagingRateMs);
frame_rate = 30;
t = 1:(size(data_tc_trial,1));
t=(t-double(stimStart))/frame_rate;

[n n2] = subplotn(nSizes*nCons);
x=1;
figure;
    for iSize = 1:nSizes %loop through the sizes
        inds1 = find(tSize == Sizes(iSize)); %find trials with that size
        for iCon = 1:nCons
        inds2 = find(tCons == Cons(iCon));
        inds = intersect(inds1,inds2);
        temp_trials = squeeze(nanmean(data_tc_trial(:,inds,goodFitResp),2));
        meanTC_byCondition(:,iSize,iCon,:)=temp_trials;
        temp_mean = nanmean(temp_trials,2);
        temp_se = std(temp_trials,[],2)/sqrt(length(goodFitResp));

        subplot(n,n2,x)

        shadedErrorBar(t,temp_mean,temp_se);
        hold on
        fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
        hold on
        fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
        hold on
        ylim([-.03 .1])
        xlim([-2 2])
        
        hline(0)
        hold off
        title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
        x=x+1;
        end
        

    end
  
sgtitle([mouse, ', ', num2str(length(goodFitResp)),' cells'])
print(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse, date, ImgFolder, ['tc_matrix.pdf']), '-dpdf')
%% same as above split by cell type
pyrCells = intersect(goodFitResp, find(~mask_label));
interNrns = intersect(goodFitResp, find(mask_label));


meanTC_byCondition_Pyr=nan(nOn+nOff,nSizes,nCons,length(pyrCells));
meanTC_byCondition_IN=nan(nOn+nOff,nSizes,nCons,length(interNrns));
%frame_rate=double(input.frameImagingRateMs);
frame_rate = 30;
t = 1:(size(data_tc_trial,1));
t=(t-double(stimStart))/frame_rate;

[n n2] = subplotn(nSizes*nCons);
x=1;
figure;
    for iSize = 1:nSizes %loop through the sizes
        inds1 = find(tSize == Sizes(iSize)); %find trials with that size
        for iCon = 1:nCons
        inds2 = find(tCons == Cons(iCon));
        inds = intersect(inds1,inds2);
        temp_trials1 = squeeze(nanmean(data_tc_trial(:,inds,pyrCells),2));
        meanTC_byCondition_Pyr(:,iSize,iCon,:)=temp_trials1;
        temp_mean1 = nanmean(temp_trials1,2);
        temp_se1 = std(temp_trials1,[],2)/sqrt(length(pyrCells));

        temp_trials2 = squeeze(nanmean(data_tc_trial(:,inds,interNrns),2));
        meanTC_byCondition_IN(:,iSize,iCon,:)=temp_trials2;
        temp_mean2 = nanmean(temp_trials2,2);
        temp_se2 = std(temp_trials2,[],2)/sqrt(length(interNrns));

        subplot(n,n2,x)

        shadedErrorBar(t,temp_mean1,temp_se1);
        hold on
        shadedErrorBar(t,temp_mean2,temp_se2,'r');
        hold on
        fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
        hold on
        fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
        hold on
        ylim([-.03 .1])
        xlim([-2 2])
        
        hline(0)
        hold off
        title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
        x=x+1;
        end
        

    end
  
sgtitle([mouse, ', ', num2str(length(goodFitResp)),' cells'])
%% for checking out individual size/con combinations if something looks unusual
        inds1 = find(tSize == 30);
        inds2 = find(tCons == .4);
        inds = intersect(inds1,inds2);
         temp_trials = squeeze(nanmean(data_tc_trial(:,inds,resp_all),2));
        temp_mean = nanmean(temp_trials,2);
        temp_se = std(temp_trials,[],2)/sqrt(sum(resp_all));
        figure
        shadedErrorBar(t(49:82),temp_mean(49:82,:),temp_se(49:82,:));
        hold on
        ylim([-.05 .13])
        vline(0)

%% looking at cells individually


for iCell = 1:length(goodFitResp)

figure;
[n n2] = subplotn(nSizes*nCons);
x=1;
cellInd = goodFitResp(iCell)
    for iSize = 1:nSizes %loop through the sizes
        inds1 = find(tSize == Sizes(iSize)); %find trials with that size
        for iCon = 1:nCons
        inds2 = find(tCons == Cons(iCon));
        inds = intersect(inds1,inds2);
        temp_trials = data_tc_trial(:,inds,cellInd);
        temp_mean = nanmean(temp_trials,2);
        temp_se = std(temp_trials,[],2)/sqrt(size(temp_trials,2));

        subplot(n,n2,x)

        
        fill([.2 .2 .4 .4],[-.1 .5 .5 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
        hold on
        fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
        hold on
        shadedErrorBar(t,temp_mean,temp_se);
        hold on
        ylim([-.03 .2])
        xlim([-2 2])
        
        hline(0)
        hold off
        title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
        x=x+1;
        end
        

    end
    sgtitle(num2str(cellInd));
    print(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse, date, ImgFolder, [num2str(cellInd),'_matrix.pdf']), '-dpdf')
end


%% looking at wheel speed
wheel_speed = wheelSpeedCalc(input,32,'orange'); 
nanmean(wheel_speed)



wheel_tc = zeros(nOn+nOff, nTrials);

for iTrial = 1:nTrials
    wheel_tc(:,iTrial) = wheel_speed(1+((iTrial-1).*(nOn+nOff)):iTrial.*(nOn+nOff));
end
wheel_trial_avg = mean(wheel_tc(nOff:nOn+nOff,:),1);
RIx = wheel_trial_avg>2; %.55 is the noise level in the wheel movement
mean(RIx)

%%
TC_byConditionLoc=nan(nOn+nOff,nSizes,nCons,length(interNrns));
TC_byConditionStat=nan(nOn+nOff,nSizes,nCons,length(interNrns));

[n n2] = subplotn(nSizes*nCons);
x=1;
figure;
    for iSize = 1:nSizes %loop through the sizes
        inds1 = find(tSize == Sizes(iSize)); %find trials with that size
        for iCon = 1:nCons
        inds2 = find(tCons == Cons(iCon));
        inds3 = intersect(inds1,inds2);
        inds4=intersect(inds3,find(RIx));
        inds5=intersect(inds3,setdiff(1:nTrials,find(RIx)));
        temp_trials = squeeze(nanmean(data_tc_trial(:,inds4,interNrns),2));
        TC_byConditionLoc(:,iSize,iCon,:)=temp_trials;
        temp_mean = nanmean(temp_trials,2);
        temp_se = std(temp_trials,[],2)/sqrt(length(interNrns));
        
        temp_trials2 = squeeze(nanmean(data_tc_trial(:,inds5,interNrns),2));
        TC_byConditionStat(:,iSize,iCon,:)=temp_trials2;
        temp_mean2 = nanmean(temp_trials2,2);
        temp_se2 = std(temp_trials2,[],2)/sqrt(length(interNrns));

        subplot(n,n2,x)

        shadedErrorBar(t,temp_mean,temp_se,'b');
        hold on
        shadedErrorBar(t,temp_mean2,temp_se2);
%         hold on
%         fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
        hold on
        fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
        hold on
        ylim([-.03 .1])
        xlim([-2 2])
        
        hline(0)
        hold off
        title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
        x=x+1;
        end
        

    end
  
sgtitle([mouse, ', ', num2str(length(interNrns)),' INs'])

%save(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse, date, ImgFolder,'TCs.mat'),'TC_byConditionStat','TC_byConditionLoc','meanTC_byCondition')
%%  ratio of dF/F in response window vs. late window(s)
responseWins = nan(3,nSizes,nCons,length(goodFitResp));

for iSize = 1:nSizes
    for iCon = 1:nCons
        for iCell = 1:length(goodFitResp)
            responseWins(1,iSize,iCon,iCell) = nanmean(meanTC_byCondition(62:65,iSize,iCon,iCell),1);
            responseWins(2,iSize,iCon,iCell) = nanmean(meanTC_byCondition(67:70,iSize,iCon,iCell),1);
            responseWins(3,iSize,iCon,iCell) = nanmean(meanTC_byCondition(76:106,iSize,iCon,iCell),1);
  
        end

    end
end
clear resp1 resp2 resp3
respWins_rect = responseWins;
respWins_rect( respWins_rect <= 0 ) = nan;

save(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse, date, ImgFolder,'responses.mat'),'respWins_rect','responseWins')

%%

[n n2] = subplotn(nSizes*nCons);
x=1;
figure;
    for iSize = 1:nSizes %loop through the sizes
        inds1 = find(tSize == Sizes(iSize)); %find trials with that size
        for iCon = 1:nCons

        subplot(n,n2,x)
        
        x_axis = [nanmean(responseWins(1,iSize,iCon,:)),nanmean(responseWins(2,iSize,iCon,:)),nanmean(responseWins(3,iSize,iCon,:))];
        y = [std(responseWins(1,iSize,iCon,:))/sqrt(size(responseWins,4)),std(responseWins(2,iSize,iCon,:))/sqrt(size(responseWins,4)),std(responseWins(3,iSize,iCon,:))/sqrt(size(responseWins,4))];

        labs =categorical({'Win1','Win2','Win3'});
        bar(labs,x_axis)                
        hold on
        er = errorbar(labs,x_axis,-y,y,LineStyle = 'none');    
        er.Color = [0 0 0]; 
        title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
        x=x+1;
        end
        

    end
  
sgtitle([mouse, ', ', num2str(length(goodFitResp)),' cells'])

%% 
DistCutoffs=[0,5,10,20];

[n n2] = subplotn(1*(length(DistCutoffs)-1));
x=1;
figure;
    for iDist = 1:(length(DistCutoffs)-1) %loop through the sizes
          minDist=DistCutoffs(iDist);
          maxDist=DistCutoffs(iDist+1);
         
          cellIndsTemp=find(keepDists>minDist & keepDists<maxDist)
        inds1 = find(tSize == Sizes(5)); %find trials with that size
        inds2 = find(tCons == Cons(4));
        inds = intersect(inds1,inds2);
        temp_trials = squeeze(nanmean(data_tc_trial(:,inds,cellIndsTemp),2));
        temp_mean = nanmean(temp_trials,2);
        temp_se = std(temp_trials,[],2)/sqrt(length(cellIndsTemp));

        subplot(n,n2,x)

        shadedErrorBar(t(60:90),temp_mean(60:90),temp_se(60:90));
        hold on
        fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
        hold on
        fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
        hold on
        ylim([-.03 .1])
        xlim([-2 2])
        
        hline(0)
        hold off
        title([num2str(minDist) ' to ' num2str(maxDist) ', n = ' num2str(length(cellIndsTemp))] )        
        x=x+1;
 
        

    end
