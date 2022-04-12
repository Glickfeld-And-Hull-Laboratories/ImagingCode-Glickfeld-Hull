clear all; clear global; close all


%% WK12
mouse = 'WK12';
date = '220107';
time = char('1532');
ImgFolder = char('003');
RetImgFolder = char('002');

frame_rate = 30; %enter the frame rate, or I can edit this to enter the stimulus duration

%setting my paths
fn_base = 'Z:\home\Celine\Analysis\2p_analysis\';
fn = fullfile(fn_base,mouse,date,ImgFolder);
mkdir(fn);
cd(fn);

beh_prefix = strcat('Z:\Behavior\Data\data-');
beh_file = [beh_prefix mouse '-' date '-' time '.mat'];
load(beh_file); %load the mworks behavioral file

run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];


load('TCs.mat');
temp_meanTC1 = meanTC_byCondition;
tempTC_loc1=TC_byConditionLoc;
tempTC_stat1=TC_byConditionStat;
load('responses.mat');
temp_responseWins1 = responseWins;
%% WK16

mouse = 'WK16';
date = '220315';
time = char('1654');
ImgFolder = char('004');
RetImgFolder = char('003');


frame_rate = 30; %enter the frame rate, or I can edit this to enter the stimulus duration

%setting my paths
fn_base = 'Z:\home\Celine\Analysis\2p_analysis\';
fn = fullfile(fn_base,mouse,date,ImgFolder);
mkdir(fn);
cd(fn);

beh_prefix = strcat('Z:\Behavior\Data\data-');
beh_file = [beh_prefix mouse '-' date '-' time '.mat'];
load(beh_file); %load the mworks behavioral file

run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];


load('TCs.mat');
temp_meanTC2 = meanTC_byCondition;
tempTC_loc2=TC_byConditionLoc;
tempTC_stat2=TC_byConditionStat;
load('responses.mat');
temp_responseWins2 = responseWins;
%% WK17
mouse = 'WK17';
date = '220327';
time = char('1504');
ImgFolder = char('004');
RetImgFolder = char('003');

frame_rate = 30; %enter the frame rate, or I can edit this to enter the stimulus duration

%setting my paths
fn_base = 'Z:\home\Celine\Analysis\2p_analysis\';
fn = fullfile(fn_base,mouse,date,ImgFolder);
mkdir(fn);
cd(fn);

beh_prefix = strcat('Z:\Behavior\Data\data-');
beh_file = [beh_prefix mouse '-' date '-' time '.mat'];
load(beh_file); %load the mworks behavioral file

run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];


load('TCs.mat');
temp_meanTC3 = meanTC_byCondition;
tempTC_loc3=TC_byConditionLoc;
tempTC_stat3=TC_byConditionStat;
load('responses.mat');
temp_responseWins3 = responseWins;

%% WK11 
mouse = 'WK11';
date = '220114';
time = char('1410');
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

run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];


load('TCs.mat');
temp_meanTC4 = meanTC_byCondition;
tempTC_loc4=TC_byConditionLoc;
tempTC_stat4=TC_byConditionStat;
load('responses.mat');
temp_responseWins4 = responseWins;
%% WK15
mouse = 'WK15';
date = '220309';
time = char('1653');
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

run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];


load('TCs.mat');
temp_meanTC5 = meanTC_byCondition;
tempTC_loc5=TC_byConditionLoc;
tempTC_stat5=TC_byConditionStat;
load('responses.mat');
temp_responseWins5 = responseWins;
%%
meanTC_byCondition = cat(4,temp_meanTC1,temp_meanTC2,temp_meanTC3,temp_meanTC4,temp_meanTC5);
TC_byConditionLoc = cat(4,tempTC_loc1,tempTC_loc2,tempTC_loc3,tempTC_loc4,tempTC_loc5);
TC_byConditionStat = cat(4,tempTC_stat1,tempTC_stat2,tempTC_stat3,tempTC_stat4,tempTC_stat5);
responseWins = cat(4,temp_responseWins1,temp_responseWins2,temp_responseWins3,temp_responseWins4,temp_responseWins5);

clear temp_meanTC1 temp_meanTC2 temp_meanTC3 temp_meanTC4 temp_meanTC5 tempTC_loc1 tempTC_loc2 tempTC_loc3 tempTC_loc4 tempTC_loc5 tempTC_stat1 tempTC_stat2 tempTC_stat3 tempTC_stat4 tempTC_stat5 temp_responseWins1 temp_responseWins2 temp_responseWins3 temp_responseWins4 temp_responseWins5
%%
nTrialsOG=length(input.tGratingDiameterDeg);
nTrials= nTrialsOG-1;
nOn = input.nScansOn;
nOff=input.nScansOff;
tCons = celleqel2mat_padded(input.tGratingContrast(1:nTrials-1)); %transforms cell array into matrix (1 x ntrials)
Cons = unique(tCons);
nCons = length(Cons);

tSize = celleqel2mat_padded(input.tGratingDiameterDeg(1:nTrials-1)); %transforms cell array into matrix (1 x ntrials)
Sizes = unique(tSize);
nSizes = length(Sizes);
%%
stimStart = (nOff/2)+1; %this indicates both the perdiod to trim off the start and the stim on period after trimming
stimEnd=stimStart+nOn-1;

t = 1:(size(meanTC_byCondition,1));
t=(t-double(stimStart))/frame_rate;
nCells = size(meanTC_byCondition,4);
[n n2] = subplotn(nSizes*nCons);
x=1;
figure;
    for iSize = 1:nSizes %loop through the sizes
        
        for iCon = 1:nCons

        temp_mean = nanmean(meanTC_byCondition(:,iSize,iCon,:),4);
        temp_se = std(meanTC_byCondition(:,iSize,iCon,:),[],4)/sqrt(nCells);

        subplot(n,n2,x)

        shadedErrorBar(t,temp_mean,temp_se);
        hold on
        fill([t(62) t(62) t(65) t(65)],[-.1 .15 .15 -.1],'g',FaceAlpha = 0.25,LineStyle='none')
        hold on
        fill([t(67) t(67) t(70) t(70)],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
        hold on
        fill([t(76) t(76) t(106) t(106)],[-.1 .15 .15 -.1],'p',FaceAlpha = 0.25,LineStyle='none')
        hold on
        fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'k',LineStyle='none')
        hold on
        ylim([-.03 .1])
        xlim([-2 2])
        
        hline(0)
        hold off
        title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
        x=x+1;
        end
        

    end
  
sgtitle([num2str(nCells),' cells'])

%% 


[n n2] = subplotn(nSizes*nCons);
x=1;
figure;
    for iSize = 1:nSizes %loop through the sizes
        
        for iCon = 1:nCons
        
        temp_mean = nanmean(TC_byConditionLoc(:,iSize,iCon,:),4);
        temp_se = nanstd(TC_byConditionLoc(:,iSize,iCon,:),[],4)/sqrt(nCells);
        
        temp_mean2 = nanmean(TC_byConditionStat(:,iSize,iCon,:),4);
        temp_se2 = nanstd(TC_byConditionStat(:,iSize,iCon,:),[],4)/sqrt(nCells);

        subplot(n,n2,x)

        shadedErrorBar(t,temp_mean,temp_se,'b');
        hold on
        shadedErrorBar(t,temp_mean2,temp_se2);
        hold on

        fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'k',LineStyle='none')
        hold on
        ylim([-.03 .1])
        xlim([-2 2])
        
        hline(0)
        hold off
        title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
        x=x+1;
        end
        

    end
    sgtitle([num2str(nCells),' cells'])
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
  
sgtitle([ num2str(nCells),' cells'])
%% 
respDiff = squeeze(responseWins(1,:,:,:)-responseWins(2,:,:,:));
figure;
scatter(squeeze(respDiff(1,4,:)),squeeze(respDiff(5,4,:)));
xlabel('Response difference 7.5 X 0.8')
ylabel('Response difference 120 X 0.8')
refline(1)
xlim([-.16 .1])
ylim([-.16 .1])
axis square
%% chaning shapes for LG
responseWinsLG = responseWins;
responseWinsLG=permute(responseWinsLG, [4 2 3 1]);
responseWinsLG = reshape(responseWinsLG, [162 20 3]);

meanTC_byConditionLG=meanTC_byCondition;
meanTC_byConditionLG=permute(meanTC_byConditionLG,[4,2,3,1]);
meanTC_byConditionLG = reshape(meanTC_byConditionLG, [162 20 123]);

TC_byConditionLocLG = TC_byConditionLoc;
TC_byConditionLocLG=permute(TC_byConditionLocLG,[4,2,3,1]);
TC_byConditionLocLG = reshape(TC_byConditionLocLG, [162 20 123]);

TC_byConditionStatLG = TC_byConditionStat;
TC_byConditionStatLG=permute(TC_byConditionStatLG,[4,2,3,1]);
TC_byConditionStatLG = reshape(TC_byConditionStatLG, [162 20 123]);

conditions(:,1)=repmat(Sizes,1,4);
conditions(:,2)=repelem(Cons,5);

save(fullfile('Z:\home\celine\Analysis\2p_analysis\suppressionDataforLG_033122','supressionData.mat'),'TC_byConditionStat','TC_byConditionLoc','meanTC_byCondition','responseWins')
save(fullfile('Z:\home\celine\Analysis\2p_analysis\suppressionDataforLG_033122','supressionDataforLG.mat'),'TC_byConditionStatLG','TC_byConditionLocLG','meanTC_byConditionLG','responseWinsLG','conditions')