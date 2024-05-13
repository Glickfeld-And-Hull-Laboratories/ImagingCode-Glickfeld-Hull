% to analyze individual day's data
clear all; clear global; close all

%identifying animal and run
netSupp_expt

expt_num = 25;

mouse = expt(expt_num).mouse
date = expt(expt_num).date;
time = expt(expt_num).contrastxsize_time;
ImgFolder = expt(expt_num).contrastxsize_runs;
RetImgFolder = expt(expt_num).retRun;
ori_run=expt(expt_num).oriRun;
depth = expt(expt_num).z;
frame_rate = expt(expt_num).frame_rate; 

if computer == 'GLNXA64'
    isilonName =  '/home/cc735@dhe.duke.edu/GlickfeldLabShare';
    base = fullfile('/All_Staff/home/ACh/Analysis/2p_analysis');
    beh_prefix = strcat(isilonName,'/All_Staff/Behavior/Data/data-');
else
    isilonName = 'duhs-user-nc1.dhe.duke.edu/';
    base = fullfile('/home/ACh/Analysis/2p_analysis');
   
   beh_prefix = strcat('Z:\Behavior\Data\data-');
end
%setting my paths
fn = fullfile(base,mouse,date,ImgFolder);
mkdir(fn);
cd(fn);

mkdir(fullfile(base, mouse, date, depth));

beh_file = [beh_prefix mouse '-' date '-' time '.mat'];
load(beh_file); %load the mworks behavioral file

run_str = ['runs-' ImgFolder];
ori_run_str = ['runs-' ori_run];
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];

% load data
load(fullfile([datemouserun '_TCs.mat']));
% loads 'mask_cell', 'mask_np'
load(fullfile(base, mouse,date,RetImgFolder, [date '_' mouse '_runs-' RetImgFolder '_lbub_fits.mat']))
load([date '_' mouse '_' run_str '_centerCells.mat']); %load the centerCells matrix
load(fullfile(base, mouse,date,ori_run, [date '_' mouse '_runs-' ori_run '_prefOris.mat']))
% load in the ori preferences

%NEED TO LOAD MASKS TO GET THE INDEXES OF RED CELLS

load(fullfile(base, mouse, date, RetImgFolder, [date '_' mouse '_runs-' RetImgFolder '_mask_cell.mat']))
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
%plot(squeeze(mean(data_tc_trial, 2)))
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
        ind = intersect(ind_size,ind_con); %for every size and then every contrast, find trials with that con/ori combination
        data_resp(:,iSize,iCon,1) = squeeze(mean(mean(data_tc_trial(resp_win,ind,:),1),2));
        data_resp(:,iSize,iCon,2) = squeeze(std(mean(data_tc_trial(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
        [h(:,iSize,iCon), p(:,iSize,iCon)] = ttest(mean(data_tc_trial(resp_win,ind,:),1), mean(data_tc_trial(base_win,ind,:),1),'dim',2,'tail','right','alpha',0.05);
        %[h(:,iSize,iCon), p(:,iSize,iCon)] = ttest(mean(data_tc_trial(resp_win,ind,:),1), mean(data_tc_trial(base_win,ind,:),1),'dim',2,'tail','right','alpha',0.05./((nSizes*nCons)-1));
    end
end
% sizeOnly = squeeze(mean(data_resp(:,:,:,1),3));
% for iCell = 1:nCells
%    [~,prefSize(1,iCell)] = max(sizeOnly(iCell,:,:,1),[],2); 
% end

h_all = sum(sum(h,2),3);
resp_any=logical(h_all);
% resp_small_highCon = logical(h(:,1,4));


centeredResp=intersect(find(resp_any),centerCells);
% smallResp=intersect(centeredResp,find(prefSize<4));
% length(centeredResp)
goodFitResp=intersect(find(resp_any),goodfit_ind);
length(goodFitResp)
keepDists = cellDists(goodFitResp);
h_keep=h(goodFitResp,:,:);%subset the h matrix to only the cells that were responsive to at least one stimulus
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
%print(fullfile(base, mouse, date, depth, ['tc_matrix.pdf']), '-dpdf')
%% same as above split by cell type
interNrns = mask_label(goodFitResp);
interNrns_inds= find(interNrns);
pyrCells_inds = find(~interNrns_inds); %DOUBLE CHECK

%


meanTC_byCondition_Pyr=nan(nOn+nOff,nSizes,nCons,length(pyrCells_inds));
meanTC_byCondition_IN=nan(nOn+nOff,nSizes,nCons,length(interNrns_inds));
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
        temp_trials1 = squeeze(nanmean(data_tc_trial(:,inds,pyrCells_inds),2));
        meanTC_byCondition_Pyr(:,iSize,iCon,:)=temp_trials1;
        temp_mean1 = nanmean(temp_trials1,2);
        temp_se1 = std(temp_trials1,[],2)/sqrt(length(pyrCells_inds));

        temp_trials2 = squeeze(nanmean(data_tc_trial(:,inds,interNrns_inds),2));
        meanTC_byCondition_IN(:,iSize,iCon,:)=temp_trials2;
        temp_mean2 = nanmean(temp_trials2,2);
        temp_se2 = std(temp_trials2,[],2)/sqrt(length(interNrns_inds));

        subplot(n,n2,x)

        shadedErrorBar(t,temp_mean1,temp_se1);
        hold on
        shadedErrorBar(t,temp_mean2,temp_se2,'r');
        hold on
        %fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
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
        
clear temp_mean1 temp_trials1 temp_se1 temp_mean2 temp_trials2 temp_se2
    end
  
sgtitle([mouse, ', ', num2str(length(goodFitResp)),' cells'])
% %% for checking out individual size/con combinations if something looks unusual
%         inds1 = find(tSize == 30);
%         inds2 = find(tCons == .4);
%         inds = intersect(inds1,inds2);
%          temp_trials = squeeze(nanmean(data_tc_trial(:,inds,resp_all),2));
%         temp_mean = nanmean(temp_trials,2);
%         temp_se = std(temp_trials,[],2)/sqrt(sum(resp_all));
%         figure
%         shadedErrorBar(t(49:82),temp_mean(49:82,:),temp_se(49:82,:));
%         hold on
%         ylim([-.05 .13])
%         vline(0)

% %% looking at cells individually
% 
% 
% for iCell = 1:length(goodFitResp)
% 
% figure;
% [n n2] = subplotn(nSizes*nCons);
% x=1;
% cellInd = goodFitResp(iCell)
%     for iSize = 1:nSizes %loop through the sizes
%         inds1 = find(tSize == Sizes(iSize)); %find trials with that size
%         for iCon = 1:nCons
%         inds2 = find(tCons == Cons(iCon));
%         inds = intersect(inds1,inds2);
%         temp_trials = data_tc_trial(:,inds,cellInd);
%         temp_mean = nanmean(temp_trials,2);
%         temp_se = std(temp_trials,[],2)/sqrt(size(temp_trials,2));
% 
%         subplot(n,n2,x)
% 
%         
%         fill([.2 .2 .4 .4],[-.1 .5 .5 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
%         hold on
%         fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
%         hold on
%         shadedErrorBar(t,temp_mean,temp_se);
%         hold on
%         ylim([-.03 .2])
%         xlim([-2 2])
%         
%         hline(0)
%         hold off
%         title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
%         x=x+1;
%         end
%         
% 
%     end
%     sgtitle(num2str(cellInd));
%     print(fullfile(base, mouse, date, ImgFolder, [num2str(cellInd),'_matrix.pdf']), '-dpdf')
% end


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
sum(RIx)
%%

pyrCells_inds = intersect(goodFitResp, find(~mask_label));
interNrns_inds = intersect(goodFitResp, find(mask_label));

TC_byConditionStat=nan(nOn+nOff,nSizes,nCons,length(goodFitResp));
TC_byConditionLoc=nan(nOn+nOff,nSizes,nCons,length(goodFitResp));
TC_byConditionLoc_INs=nan(nOn+nOff,nSizes,nCons,length(interNrns_inds));
TC_byConditionStat_INs=nan(nOn+nOff,nSizes,nCons,length(interNrns_inds));

[n n2] = subplotn(nSizes*nCons);
x=1;
figure;
    for iSize = 1:nSizes %loop through the sizes
        inds1 = find(tSize == Sizes(iSize)); %find trials with that size
        for iCon = 1:nCons
        inds2 = find(tCons == Cons(iCon));
        inds3 = intersect(inds1,inds2);
        inds4=intersect(inds3,find(RIx));
        inds5=intersect(inds3,setdiff(1:nTrials,find(RIx))); %find the trials NOT in RIx

        temp_trials = squeeze(nanmean(data_tc_trial(:,inds4,interNrns_inds),2));
        TC_byConditionLoc_INs(:,iSize,iCon,:)=temp_trials;
        temp_mean = nanmean(temp_trials,2);
        temp_se = std(temp_trials,[],2)/sqrt(length(interNrns_inds));
        
        temp_trials2 = squeeze(nanmean(data_tc_trial(:,inds5,interNrns_inds),2));
        TC_byConditionStat_INs(:,iSize,iCon,:)=temp_trials2;
        temp_mean2 = nanmean(temp_trials2,2);
        temp_se2 = std(temp_trials2,[],2)/sqrt(length(interNrns_inds));

        TC_byConditionStat(:,iSize,iCon,:)=squeeze(nanmean(data_tc_trial(:,inds5,goodFitResp),2));
        TC_byConditionLoc(:,iSize,iCon,:)=squeeze(nanmean(data_tc_trial(:,inds4,goodFitResp),2));
        
        subplot(n,n2,x)

       % shadedErrorBar(t,temp_mean,temp_se,'b');
        hold on
        shadedErrorBar(t,temp_mean2,temp_se2);
        hold on
        %fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
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
  
sgtitle([mouse, ', ', num2str(length(interNrns_inds)),' INs'])
print(fullfile(base, mouse, date, depth, ['IntrNrn_tc_matrix.pdf']), '-dpdf')


TC_byConditionLoc_Pyr=nan(nOn+nOff,nSizes,nCons,length(pyrCells_inds));
TC_byConditionStat_Pyr=nan(nOn+nOff,nSizes,nCons,length(pyrCells_inds));

[n n2] = subplotn(nSizes*nCons);
x=1;
figure;
    for iSize = 1:nSizes %loop through the sizes
        inds1 = find(tSize == Sizes(iSize)); %find trials with that size
        for iCon = 1:nCons
        inds2 = find(tCons == Cons(iCon));
        inds3 = intersect(inds1,inds2);
        inds4=intersect(inds3,find(RIx));
        length(inds4)
        inds5=intersect(inds3,setdiff(1:nTrials,find(RIx)));
        temp_trials = squeeze(nanmean(data_tc_trial(:,inds4,pyrCells_inds),2));
        TC_byConditionLoc_Pyr(:,iSize,iCon,:)=temp_trials;
        temp_mean = nanmean(temp_trials,2);
        temp_se = std(temp_trials,[],2)/sqrt(length(pyrCells_inds));
        
        temp_trials2 = squeeze(nanmean(data_tc_trial(:,inds5,pyrCells_inds),2));
        TC_byConditionStat_Pyr(:,iSize,iCon,:)=temp_trials2;
        temp_mean2 = nanmean(temp_trials2,2);
        temp_se2 = std(temp_trials2,[],2)/sqrt(length(pyrCells_inds));

        subplot(n,n2,x)

        %shadedErrorBar(t,temp_mean,temp_se,'b');
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
  
sgtitle([mouse, ', ', num2str(length(pyrCells_inds)),' Pyr'])
print(fullfile(base, mouse, date, depth, ['Pyr_tc_matrix.pdf']), '-dpdf')

save(fullfile(base, mouse, date, ImgFolder,'TCs.mat'),'TC_byConditionStat_INs','TC_byConditionLoc_INs','meanTC_byCondition')
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

save(fullfile(base, mouse, date, ImgFolder,'responses.mat'),'respWins_rect','responseWins')

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

%% making seperate matrices for cells with different RF center distances
DistCutoffs=[0,5,10,20];

    for iDist = 1:(length(DistCutoffs)-1) %loop through the sizes
          minDist=DistCutoffs(iDist);
          maxDist=DistCutoffs(iDist+1);
         
        cellIndsTemp=find(keepDists>minDist & keepDists<maxDist);
        meanTC_byCondition=nan(nOn+nOff,nSizes,nCons,length(cellIndsTemp));
        figure;
        [n n2] = subplotn(nSizes*nCons);
        x=1;
        for iSize = 1:nSizes %loop through the sizes
            inds1 = find(tSize == Sizes(iSize)); %find trials with that size
            for iCon = 1:nCons
                inds2 = find(tCons == Cons(iCon));
                inds = intersect(inds1,inds2);
                temp_trials1 = squeeze(nanmean(data_tc_trial(:,inds,cellIndsTemp),2));
                meanTC_byCondition(:,iSize,iCon,:)=temp_trials1;
                temp_mean1 = nanmean(temp_trials1,2);
                temp_se1 = std(temp_trials1,[],2)/sqrt(length(cellIndsTemp));
        
                subplot(n,n2,x)
        
                shadedErrorBar(t,temp_mean1,temp_se1);

                hold on
                %fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
                hold on
                fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
                hold on
                ylim([-.03 .1])
                xlim([-2 2])
                
                hline(0)
                hold off
                title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
                x=x+1;
               sgtitle([num2str(minDist) ' to ' num2str(maxDist) ', n = ' num2str(length(cellIndsTemp))] ) 
              
            end
        

    end
 
    print(fullfile(base, mouse, date, depth, ['distance',num2str(iDist),'_matrix.pdf']), '-dpdf');

    end

    

%% making seperate matrices for cells with different pref ori
keepPrefOri = prefOri(goodFitResp); %make a subset of preOri for the responsive cells

 oris=unique(prefOri);

    for iOri = 1:length(oris) %loop through the sizes

        oriInds = find(prefOri == oris(iOri));
        cellIndsTemp=find(intersect(goodFitResp,oriInds));
        meanTC_byCondition=nan(nOn+nOff,nSizes,nCons,length(cellIndsTemp));

        figure;
        [n n2] = subplotn(nSizes*nCons);
        
        x=1;
        for iSize = 1:nSizes %loop through the sizes
                inds1 = find(tSize == Sizes(iSize)); %find trials with that size
                for iCon = 1:nCons
                    inds2 = find(tCons == Cons(iCon));
                    inds = intersect(inds1,inds2);
                    temp_trials1 = squeeze(nanmean(data_tc_trial(:,inds,cellIndsTemp),2)); %one average per cell
                    meanTC_byCondition(:,iSize,iCon,:)=temp_trials1;
                    temp_mean1 = nanmean(temp_trials1,2);
                    temp_se1 = std(temp_trials1,[],2)/sqrt(length(cellIndsTemp));
            
                    subplot(n,n2,x)
            
                    shadedErrorBar(t,temp_mean1,temp_se1);
    
                    hold on
                    %fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
                    hold on
                    fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
                    hold on
                    ylim([-.03 .1])
                    xlim([-2 2])
                    
                    hline(0)
                    hold off
                    title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
                    x=x+1;
                   sgtitle(['pref ori ', num2str(oris(iOri)) ', n = ' num2str(length(cellIndsTemp))] ) 
              
                end
        

    end
 
        print(fullfile(base, mouse, date, depth, ['prefOri',num2str(iOri),'_matrix.pdf']), '-dpdf');

    end
%% save data summary
save(fullfile(base, mouse, date,ImgFolder,'goodFitResp_summary_updated.mat'), ...
    'TC_byConditionLoc','TC_byConditionStat','h_keep','interNrns','keepDists','keepPrefOri')