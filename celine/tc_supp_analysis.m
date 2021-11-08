% to analyze individual day's data
clear all; clear global; close all

%identifying animal and run
date = '211004';
imgFolder = '003';
time = '1457';
mouse = 'i2015';

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
nTrials=360;

cellTCs_OG = npSub_tc;

%use the list of direction by trial to figure out how many trials there are
%currently the way I center the stim on period requires me to cut out
%one trial, hence the -1

%for matched cell
nCells = size(npSub_tc,2);

npSub_tc = npSub_tc(stimStart:size(npSub_tc,1),:);
npSub_tc = padarray(npSub_tc,double(stimStart)-1,0,'post');
nFrames = size(npSub_tc,1);
data_trial = reshape(npSub_tc,[nOn+nOff nTrials nCells]);
data_f= mean(data_trial(1: (nOff/2),:,:),1);
data_tc_trial = bsxfun(@rdivide,bsxfun(@minus,data_trial,data_f),data_f);

clear data_trial data_f
 plot(squeeze(mean(data_tc_trial, 2)))
%% find the stimulus conditions
tCons = celleqel2mat_padded(input.tGratingContrast); %transforms cell array into matrix (1 x ntrials)
Cons = unique(tCons);
nCons = length(Cons);

tSize = celleqel2mat_padded(input.tGratingDiameterDeg); %transforms cell array into matrix (1 x ntrials)
Sizes = unique(tSize);
nSizes = length(Sizes);

%% find preferred size and contrast, find cells that are responsive 

data_resp = zeros(nCells,nSizes,nCons,2);
h = zeros(nCells, nSizes,nCons);
p = zeros(nCells, nSizes,nCons);

resp_win = stimStart:stimEnd;
base_win = (stimStart - input.frameImagingRateMs/2):stimStart-1;



for iSize = 1:nSizes
    ind_size = find(tSize == Sizes(iSize));
    for iCon = 1:nCons
        ind_con = find(tCons == Cons(iCon));
        ind = intersect(ind_size,ind_con); %for every orientation and then every contrast, find trials with that con/ori combination
        data_resp(:,iSize,iCon,1) = squeeze(mean(mean(data_tc_trial(resp_win,ind,:),1),2));
        data_resp(:,iSize,iCon,2) = squeeze(std(mean(data_tc_trial(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
        [h(:,iSize,iCon), p(:,iSize,iCon)] = ttest(mean(data_tc_trial(resp_win,ind,:),1), mean(data_tc_trial(base_win,ind,:),1),'dim',2,'tail','right','alpha',0.05./(nSizes.*3-1));
    end
end

h_all = sum(sum(h,2),3);
resp=logical(h_all);
sum(resp)
respInds = find(resp);
%% plot timecourses at different sizes and contrasts
frame_rate=double(input.frameImagingRateMs);
t = 1:(size(data_tc_trial,1));
t=(t-double(stimStart))/frame_rate;

[n n2] = subplotn(nSizes*nCons);
x=1;
figure;
    for iSize = 1:nSizes %loop through the sizes
        inds1 = find(tSize == Sizes(5)); %find trials with that size
        for iCon = 1:nCons
        inds2 = find(tCons == Cons(iCon));
        inds = intersect(inds1,inds2);
        temp_trials = squeeze(nanmean(data_tc_trial(:,inds,resp),2));
        temp_mean = nanmean(temp_trials,2);
        temp_se = std(temp_trials,[],2)/sqrt(sum(resp));

        subplot(n,n2,x)
        shadedErrorBar(t,temp_mean,temp_se);
        hold on
        ylim([-.05 .1])
%         x = [0 nOn/frame_rate nOn/frame_rate 0];
%         y = [.08 .08 .082 .082 ];
%         patch(x,y,'b')
        hold off
        title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
        x=x+1
        end
        

    end
    
%% 
conInds = find(tCons == .8);
sizeInds = find(tSize==120);
temp_trials = squeeze(nanmean(data_tc_trial(:,sizeInds,resp),2));
temp_mean = nanmean(temp_trials,2);
temp_se = std(temp_trials,[],2)/sqrt(sum(resp));
shadedErrorBar(t,temp_mean,temp_se);
