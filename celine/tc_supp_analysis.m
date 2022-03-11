% to analyze individual day's data
clear all; clear global; close all

%identifying animal and run


mouse = 'WK11';
date = '220114';
time = char('1410');
imgFolder = char('002');
RetImgFolder = char('001');
frame_rate = 30; %enter the frame rate, or I can edit this to enter the stimulus duration


%setting my paths
fn_base = 'Z:\home\Celine\Analysis\2p_analysis\';
fn = fullfile(fn_base,mouse,date,imgFolder);
mkdir(fn);
cd(fn);

beh_prefix = strcat('Z:\Behavior\Data\data-');
beh_file = [beh_prefix mouse '-' date '-' time '.mat'];
load(beh_file); %load the mworks behavioral file

run_str = catRunName(imgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];

% load data
load(fullfile([datemouserun '_TCs.mat']));

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

%for matched cell
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

%% find preferred size and contrast, find cells that are responsive 

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

h_all = sum(sum(h,2),3);
resp_all=logical(h_all);
resp_small_highCon = logical(h(:,1,4));

load([date '_' mouse '_' run_str '_centerCells.mat']); %load the centerCells matrix
centeredResp=intersect(find(resp_small_highCon),centerCells);
length(centeredResp)
%% plot timecourses at different sizes and contrasts

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
        temp_trials = squeeze(nanmean(data_tc_trial(:,inds,centeredResp),2));
        temp_mean = nanmean(temp_trials,2);
        temp_se = std(temp_trials,[],2)/sqrt(length(centeredResp));

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
print(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse, date, imgFolder, ['tc_matrix.pdf']), '-dpdf')
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




for iCell = 1:length(centeredResp)
figure;
[n n2] = subplotn(nSizes*nCons);
x=1;
cellInd = centeredResp(iCell)
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
    print(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse, date, imgFolder, [num2str(cellInd),'_matrix.pdf']), '-dpdf')
end