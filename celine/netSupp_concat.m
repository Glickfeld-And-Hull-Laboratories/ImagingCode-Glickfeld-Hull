clear all; clear global; close all

%list the experiments I want
expt_list = [1,2,3,4,11,12,13,14,16,17,18,19,20,21,22,23,24,25];
nExpt = length(expt_list);

netSupp_expt

if computer == 'GLNXA64'
    isilonName =  '/home/cc735@dhe.duke.edu/GlickfeldLabShare';
    base = fullfile('/All_Staff/home/ACh/Analysis/2p_analysis');
    beh_prefix = strcat(isilonName,'/All_Staff/Behavior/Data/data-');
else
    isilonName = 'duhs-user-nc1.dhe.duke.edu/';
    base = fullfile('/home/ACh/Analysis/2p_analysis/NetworkSuppression/');
   
   beh_prefix = strcat('Z:\Behavior\Data\data-');
end

fn = fullfile(base,['netSupp_analysis-',date]);
mkdir(fn);
cd(fn);

%% collect and combine data

TCs_stat = []; 
TCs_loc = [];
interNrns_concat=[];
dists_concat=[];
prefOri_concat=[];
h_concat=[];
h_short_concat = [];
mouseID=[];
depth=[];
exptNumber = [];

for iExpt = 1:nExpt
    expt_num=expt_list(iExpt);
    %navigate to the right folder
    mouse = expt(expt_num).mouse;
    date = expt(expt_num).date;
    time = expt(expt_num).contrastxsize_time;
    frame_rate = expt(expt_num).frame_rate; 
    ImgFolder = expt(expt_num).contrastxsize_runs;
    load(fullfile(base, mouse, date,ImgFolder,'goodFitResp_summary_updated.mat'));
    beh_file = [beh_prefix mouse '-' date '-' time '.mat'];
    load(beh_file); %load the mworks behavioral file


    TCs_stat=cat(4,TCs_stat,TC_byConditionStat);
    TCs_loc=cat(4,TCs_loc,TC_byConditionLoc);
    interNrns_concat=cat(2,interNrns_concat,interNrns);
    dists_concat=cat(2,dists_concat,keepDists');
    prefOri_concat=cat(2,prefOri_concat,keepPrefOri);
    h_concat = cat(1,h_concat,h_keep);
    h_short_concat = cat(1,h_short_concat,h_short_keep);
    
    tCons = celleqel2mat_padded(input.tGratingContrast); %transforms cell array into matrix (1 x ntrials)
    Cons = unique(tCons);
    nCons = length(Cons);
    
    tSize = celleqel2mat_padded(input.tGratingDiameterDeg); %transforms cell array into matrix (1 x ntrials)
    Sizes = unique(tSize);
    nSizes = length(Sizes);

    nOn = input.nScansOn;
    nOff=input.nScansOff;
    stimStart = (nOff/2)+1; 
    stimEnd=stimStart+nOn-1;

    depth_temp = repmat(expt(expt_num).z,size(TC_byConditionStat,4),1);
    mouse_name = pad(mouse,5);
    mouseID_temp = repmat(mouse_name,size(TC_byConditionStat,4),1);
    exptNumber_temp = repmat(iExpt,size(TC_byConditionStat,4),1);


    depth=[depth;depth_temp];
    mouseID = [mouseID;mouseID_temp];
    exptNumber = vertcat(exptNumber, exptNumber_temp);


    clear TC_byConditionStat TC_byConditionLoc interNrns keepDists keepPrefOri
    clear mouse date time iExpt expt_num tCons tSize input depth_temp mouseID_temp h_keep

end

interNrns = find(interNrns_concat);
pyrCells = find(~interNrns_concat); 

nCells = size(TCs_stat,4);

TCs_stat_smooth=nan(((nOn+nOff)),nSizes,nCons,nCells);

for iSize = 1:nSizes %loop through the sizes
    for iCon = 1:nCons
        for iCell = 1:nCells
            TCs_stat_smooth(:,iSize,iCon,iCell)=smoothdata(TCs_stat(:,iSize,iCon,iCell),'movmean',3);
    
        end

    end
end

%create axis in seconds
t = 1:double(nOn+nOff);
t=(t-double(stimStart))/frame_rate;

%find cells that are well centered
centered = find(dists_concat<10);  
centerIN = intersect(centered,interNrns);
centerPyr = intersect(centered,pyrCells);

%find cells at various cortical depths
depth1 = find(depth<=175);
depth2= find(depth>175);
depth2 = intersect(depth2,find(depth<=225));
depth3 = find(depth>225);
depth3 = intersect(depth3,find(depth<=275));
depth4=find(depth>275);
depth4=intersect(depth4, find(depth<350));

PyrDepth1=(intersect(centerPyr,depth1));
INdepth1=(intersect(centerIN,depth1));

PyrDepth2=(intersect(centerPyr,depth2));
INdepth2=(intersect(centerIN,depth2));

PyrDepth3=(intersect(centerPyr,depth3));
INdepth3=(intersect(centerIN,depth3));

PyrDepth4=(intersect(centerPyr,depth4));
INdepth4=(intersect(centerIN,depth4));


PyrDepths=cell(1,4);
PyrDepths{1}=PyrDepth1;
PyrDepths{2}=PyrDepth2;
PyrDepths{3}=PyrDepth3;
PyrDepths{4}=PyrDepth4;

INdepths=cell(1,4);
INdepths{1}=INdepth1;
INdepths{2}=INdepth2;
INdepths{3}=INdepth3;
INdepths{4}=INdepth4;
%% find peak, trough and similar metrics

%make an empty matrix for values
peak_trough =nan(nSizes,nCons,nCells,2);
dip =nan(nSizes,nCons,nCells);
dip_peak_ratio =nan(nSizes,nCons,nCells);
peak_time=nan(nSizes,nCons,nCells,4);
resp_means = nan(nSizes,nCons,nCells,2);
fwhm=nan(nSizes,nCons,nCells); %only for 80% contrast
%first value will be the eary timepoint, second value will be the late
%timepoint

for iSize = 1:nSizes %loop through the sizes    
    for iCon = 1:nCons
        for iCell = 1:nCells
           

           thisTrace = TCs_stat(:,iSize,iCon,iCell);
           %interpolate data
           traceInterp=interp1(t,thisTrace,(t(1):0.01:t(123)));
           %traceInterp=smoothdata(traceInterp,'movmean',3);
           t2=t(1):0.01:t(123); % to get temporal values with the interpolated data
           stimStart_interp=find(t2==0)+5; %don't look in the first 50 ms, as this is too early to be a true peak
           %find(t2==.2)

           peak=max(traceInterp(stimStart_interp:stimStart_interp+20)); %find the max value within a set window
           %currently set to 61 (stim onset) through 67, 200ms after stim
           trough=min(traceInterp(stimStart_interp+20:stimStart_interp+30));
           %search for the trough between 200 and 300 ms after 
           
%           (thisTrace(stimStart+6:stimStart+9));

           
           peak_time_temp = find(traceInterp==peak);

           peakBin = [peak_time_temp-5,peak_time_temp+5]; %a ~100 ms bin around the peak
           
           trough_time_temp = find(traceInterp==trough);
           troughBin = [trough_time_temp-7,trough_time_temp+7]; 

           peak_time(iSize,iCon, iCell,1)=t2(peak_time_temp);
           peak_time(iSize,iCon, iCell,2)= t2(trough_time_temp);

                if peak > 0
                   if iCon > 0 %previously I used this to only take the half peak measurement for 
                       
                       halfPeak = peak/2;
                       %find the frame of the half peak
                       half_peak_frame = find(traceInterp(stimStart_interp:peak_time_temp)>halfPeak,1,'first');
                       half_peak_temp = stimStart_interp+(half_peak_frame-1); %adjust this for the frame we started on
                        %convert this to time
                       peak_time(iSize,iCon, iCell,3) = t2(half_peak_temp);
                        %find the frame of the equivalent point on the decay
                       if min(traceInterp(peak_time_temp:length(t2)))<halfPeak
                           half_dacay_frame=find(traceInterp(peak_time_temp:length(t2))<halfPeak,1,'first'); 
                           %look for the half decay anywhere after the peak, until the end of the trace
                           half_dacay_temp=peak_time_temp+half_dacay_frame-1; %adjust this for the frame we started on
                           peak_time(iSize,iCon, iCell,4) = t2(half_dacay_temp);
                            %find the difference, in time, between the half decay
                            %and the half peak
                           fwhm(iSize,iCon,iCell)=t2(half_dacay_temp)-t2(half_peak_temp);
                       end
        
                       % [minVal, endWin]=min(abs(t2-t(peak_time_temp))); %find the value in the interpoalted t vector that is closest to the time in the oringal t vector when the peak occurs
                       % half_peak_temp = find(traceInterp(stimStart_interp:endWin)>halfPeak,1,'first');
                       % peak_time(iSize,iCon, iCell,3) = t2(stimStart_interp+(half_peak_temp-1));
        
                    %dip is finding the local decrease within the general vicinty
                    %of through
                   end
                end
            dip_temp=trough-((traceInterp(troughBin(1))+traceInterp(troughBin(2)))/2);

            dip(iSize,iCon, iCell)=dip_temp;            
            dip_peak_ratio(iSize,iCon, iCell)=dip_temp/peak;
           
           %figure;plot(thisTrace);hold on; vline(peak_time_temp);vline(trough_time_temp);xlim([stimStart stimStart+15]);hold off
          
           peak_trough(iSize,iCon, iCell,1)=peak;
           %peak is the identified max
           peak_trough(iSize,iCon, iCell,2)=trough;
            %rough is the minimum within the identified timebin, which is
            %locked to the time of the peak 

            %resp_mean is the average (:,:,:,1) or cumulative integral (:,:,:,2) response from 0 to 400 ms
           resp_means(iSize,iCon, iCell,1)=mean(thisTrace(stimStart:stimStart+12));
           resp_means(iSize,iCon, iCell,2)=max(cumtrapz(t(stimStart:stimStart+12),thisTrace(stimStart:stimStart+12)));


           clear thisTrace peak peakBin troughBin dip_temp trough trough_time_temp peak_time_temp

        end
    end
end


%% look at the mean peak for each
mean_peak=(mean(peak_time(:,:,:,1),3,"omitmissing"));
mean_halfPeak=mean(peak_time(:,:,:,3),3,"omitmissing");
mean_trough=(mean(peak_time(:,:,:,2),3,"omitmissing"));
mean_halfDecay=(mean(peak_time(:,:,:,4),3,"omitmissing"));
std_peak=std(peak_time(:,:,:,1),[],3,"omitmissing");
std_trough=std(peak_time(:,:,:,2),[],3,"omitmissing");

[n n2] = subplotn(nSizes*nCons);
x=1;

figure;
    for iSize = 1:nSizes %loop through the sizes
        
        for iCon = 1:nCons

        temp_peak=mean_peak(iSize,iCon);
        temp_halfPeak=mean_halfPeak(iSize,iCon);
        temp_trough=mean_trough(iSize,iCon);
        temp_halfDecay=mean_halfDecay(iSize,iCon);
        
        temp_mean1 = mean(TCs_stat(:,iSize,iCon,centerPyr),4,"omitnan");
        temp_se1 = std(TCs_stat(:,iSize,iCon,centerPyr),[],4,"omitnan")/sqrt(length(centerPyr));

        temp_mean2 = mean(TCs_stat(:,iSize,iCon,centerIN),4,"omitnan");
        temp_se2 = std(TCs_stat(:,iSize,iCon,centerIN),[],4,"omitnan")/sqrt(length(centerIN));

        subplot(n,n2,x)
        vline(temp_halfPeak,'g')
        vline(temp_peak,'b')
        vline(temp_trough,'k')
        vline(temp_halfDecay,'r')
        hold on
        shadedErrorBar(t(:),temp_mean2,temp_se2,'r');
        hold on
        shadedErrorBar(t,temp_mean1,temp_se1);
        hold on
        alpha(.5)
        %fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
        hold on
        %fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
        hold on
        ylim([-.03 .1])
        xlim([-.25 .5])
        
        box off
        set(gca, 'TickDir', 'out')
        hline(0)
        hold off
        title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
        x=x+1;
        end
        
clear temp_mean1 temp_trials1 temp_se1 temp_mean2 temp_trials2 temp_se2
    end
  
x0=1;
y0=1;
width=5;
height=8;
set(gcf,'units','inches','position',[x0,y0,width,height])
sgtitle('Stationary')

sgtitle(['Stationary, ', num2str(length(centerPyr)),' Pyr cells, ',num2str(length(centerIN)),' SST cells'])

%% peak time traces

figure; 
subplot(1,2,1);
%peak_time(:,:,centerIN,1)
for iCon = 1:nCons
        temp_mean1 = mean(peak_time(:,iCon,centerPyr,1),3,"omitnan");
        temp_se1 = std(peak_time(:,iCon,centerPyr,1),[],3,"omitnan")/sqrt(length(centerPyr));

        errorbar(Sizes,temp_mean1,temp_se1);
        hold on
        
end

title('Centered Pyr cells')
ylabel('Peak time')
xlabel('Size')
xticks(Sizes)
ylim([.1 .17])
box off
set(gca, 'TickDir', 'out')
lgd=legend(string(Cons))
title(lgd,'Contrast')

subplot(1,2,2)
%peak_time(:,:,centerIN,1)
for iCon = 1:nCons
        temp_mean1 = mean(peak_time(:,iCon,centerIN,1),3,"omitnan");
        temp_se1 = std(peak_time(:,iCon,centerIN,1),[],3,"omitnan")/sqrt(length(centerIN));

        errorbar(Sizes,temp_mean1,temp_se1);
        hold on
end

title('Centered SST cells')
ylabel('Peak time')
xlabel('Size')
xticks(Sizes)
ylim([.1 .17])
box off
set(gca, 'TickDir', 'out')


x0=5;
y0=5;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height]);
sgtitle('Stationary')

print('Peak_time_bySizeandConstrast.pdf', '-dpdf');

%% Trough time traces

figure; 
subplot(1,2,1);
%peak_time(:,:,centerIN,1)
for iCon = 1:nCons
        temp_mean1 = mean(peak_time(:,iCon,centerPyr,2),3,"omitnan");
        temp_se1 = std(peak_time(:,iCon,centerPyr,2),[],3,"omitnan")/sqrt(length(centerPyr));

        errorbar(Sizes,temp_mean1,temp_se1);
        hold on
        
end

title('Centered Pyr cells')
ylabel('Trough time')
xlabel('Size')
xticks(Sizes)
%ylim([.1 .17])
box off
set(gca, 'TickDir', 'out')
lgd=legend(string(Cons))
title(lgd,'Contrast')

subplot(1,2,2)
%peak_time(:,:,centerIN,1)
for iCon = 1:nCons
        temp_mean1 = mean(peak_time(:,iCon,centerIN,2),3,"omitnan");
        temp_se1 = std(peak_time(:,iCon,centerIN,2),[],3,"omitnan")/sqrt(length(centerIN));

        errorbar(Sizes,temp_mean1,temp_se1);
        hold on
end

title('Centered SST cells')
ylabel('Trough time')
xlabel('Size')
xticks(Sizes)
%ylim([.1 .17])
box off
set(gca, 'TickDir', 'out')


x0=5;
y0=5;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height]);
sgtitle('Stationary')


print('Trough_time_bySizeandConstrast.pdf', '-dpdf');



% %% trough time heatmap
% 
% figure; 
% subplot(2,1,2)
% heatmap(mean(peak_time(:,:,centerIN,2),3),'Colormap',winter,'ColorLimits',[.2 .3])
% title('Centered SST cells trough time')
% ylabel('Size')
% xlabel('Contrast')
% 
% 
% subplot(2,1,1) 
% heatmap(mean(peak_time(:,:,centerPyr,2),3),'Colormap',winter,'ColorLimits',[.2 .3])
% title('Centered Pyr cells trough time')
% ylabel('Size')
% xlabel('Contrast')
% print('Trough_time_heatmap.pdf', '-dpdf');


%% dip heatmap - the local decrease to the trough
% 
% figure; 
% subplot(2,1,2) 
% heatmap(mean(dip(:,:,centerIN),3),'Colormap',winter,'ColorLimits',[-.026 -.016])
% title('Centered SST cells dip (df/f)')
% ylabel('Size')
% xlabel('Contrast')
% 
% 
% subplot(2,1,1)  
% heatmap(mean(dip(:,:,centerPyr),3),'Colormap',winter,'ColorLimits',[-.026 -.016])
% title('Centered Pyr cells dip (df/f)')
% ylabel('Size')
% xlabel('Contrast')
% 
% print('Dip_magnitude_heatmap.pdf', '-dpdf');

%% analysis on peak
PyrPeak1=peak_trough(:,:,centerPyr,1);

%flatten my matrix
PyrPeak = reshape(PyrPeak1,[20*size(PyrPeak1,3),1]);

%the flattened vector goes through all contrasts for eac size, then all
%contrasts for the next size, etc., for each cell, then starts over for the
%next cell.

%make vectors for the size and contrast factors
size_factor = repelem(Sizes,1,4);
size_factor  = repmat(size_factor,1,size(PyrPeak1,3));
con_factor = repmat(Cons,1,5);
con_factor=repmat(con_factor,1,size(PyrPeak1,3));
%make a cell identifier
cellID = 1:size(PyrPeak1,3);
cellID = repelem(cellID,1,20);


INPeak1=peak_trough(:,:,centerIN,1);

%flatten my matrix
INPeak = reshape(INPeak1,[20*size(INPeak1,3),1]);

%the flattened vector goes through all contrasts for eac size, then all
%contrasts for the next size, etc., for each cell, then starts over for the
%next cell.

%make vectors for the size and contrast factors
size_factor_IN = repelem(Sizes,1,4);
size_factor_IN  = repmat(size_factor_IN,1,size(INPeak1,3));
con_factor_IN = repmat(Cons,1,5);
con_factor_IN=repmat(con_factor_IN,1,size(INPeak1,3));
%make a cell identifier
cellID_IN = 1:size(INPeak1,3);
cellID_IN = repelem(cellID_IN,1,20);




%% plot stationary


[n n2] = subplotn(nSizes*nCons);
x=1;
figure;
    for iSize = 1:nSizes %loop through the sizes
        
        for iCon = 1:nCons
        
        temp_mean1 = mean(TCs_stat(:,iSize,iCon,centerPyr),4,"omitnan");
        temp_se1 = std(TCs_stat(:,iSize,iCon,centerPyr),[],4,"omitnan")/sqrt(length(centerPyr));


        temp_mean2 = mean(TCs_stat(:,iSize,iCon,centerIN),4,"omitnan");
        temp_se2 = std(TCs_stat(:,iSize,iCon,centerIN),[],4,"omitnan")/sqrt(length(centerIN));

        subplot(n,n2,x)

        shadedErrorBar(t(:),temp_mean2,temp_se2,'r');
        hold on
        shadedErrorBar(t(:),temp_mean1,temp_se1);
        hold on
        
        %fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
        hold on
        fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
        hold on
        ylim([-.03 .15])
        xlim([-.15 .5])
        box off
        set(gca, 'TickDir', 'out')
        hline(0)
        hold off
        title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
        x=x+1;
        end
        
clear temp_mean1 temp_trials1 temp_se1 temp_mean2 temp_trials2 temp_se2
    end
  
x0=1;
y0=1;
width=5;
height=8;
set(gcf,'units','inches','position',[x0,y0,width,height])
sgtitle(['Stationary, ', num2str(length(centerPyr)),' Pyr cells, ',num2str((length(centerIN))),' SST cells'])
print('Stationary_center_matrix.pdf', '-dpdf');

%% centered cells running trials 
%find out how many cells have running data in each condition
locResp = squeeze(mean(TCs_loc(61:68,:,:,:),1));
locCountPyr = sum(~isnan(locResp(:,:,centerPyr)),3)
locCountSST = sum(~isnan(locResp(:,:,centerIN)),3)
%% plot running


[n n2] = subplotn(nSizes*nCons);
x=1;
figure;
    for iSize = 1:nSizes %loop through the sizes
        
        for iCon = 1:nCons
         
        
        temp_mean1 = mean(TCs_loc(:,iSize,iCon,centerPyr),4,"omitnan");
        temp_se1 = std(TCs_loc(:,iSize,iCon,centerPyr),[],4,"omitnan")/sqrt(locCountPyr(iSize,iCon));


        temp_mean2 = mean(TCs_loc(:,iSize,iCon,centerIN),4,"omitnan");
        temp_se2 = std(TCs_loc(:,iSize,iCon,centerIN),[],4,"omitnan")/sqrt(locCountSST(iSize,iCon));

        subplot(n,n2,x)

        shadedErrorBar(t(:),temp_mean2,temp_se2,'r');
        hold on
        shadedErrorBar(t(:),temp_mean1,temp_se1);
        hold on
        
        %fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
        hold on
        fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
        hold on
        ylim([-.03 .15])
        xlim([-.15 .5])
        box off
        set(gca, 'TickDir', 'out')
        hline(0)
        hold off
        title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
        x=x+1;
        end
        
clear temp_mean1 temp_trials1 temp_se1 temp_mean2 temp_trials2 temp_se2
    end
  
sgtitle(['Running, ~', num2str(round(mean(mean(locCountPyr)))),' Pyr cells, ~',num2str(round(mean(mean(locCountSST)))),' SST cells'])
print('Running_center_matrix.pdf', '-dpdf');

%% plot running averaging over contrast

x=1;
y=nSizes+1;
figure;
    for iSize = 1:nSizes %loop through the sizes
        ContrastMeansLoc=nan(nOn+nOff,nCells);
        ContrastMeansStat=nan(nOn+nOff,nCells);

        for iCell = 1:nCells
        
            ContrastMeansLoc(:,iCell) = mean(TCs_loc(:,iSize,iCon,iCell),3,"omitnan");
            ContrastMeansStat(:,iCell) = mean(TCs_stat(:,iSize,iCon,iCell),3,"omitnan");

        end

    
        temp_mean1=mean(ContrastMeansLoc(:,centerPyr),2,"omitnan");
        nPyrLoc = length(find(~isnan(ContrastMeansLoc(1,centerPyr))));
        temp_se1 = std(ContrastMeansLoc(:,centerPyr),[],2,"omitnan")/sqrt(nPyrLoc);


        temp_mean2=mean(ContrastMeansStat(:,centerPyr),2,"omitnan");
        nPyrStat = length(find(~isnan(ContrastMeansStat(1,centerPyr))));
        temp_se2 = std(ContrastMeansStat(:,centerPyr),[],2,"omitnan")/sqrt(nPyrStat);


        subplot(2,nSizes,x)
        shadedErrorBar(t(:),temp_mean1,temp_se1,'m');
        hold on
        shadedErrorBar(t(:),temp_mean2,temp_se2,'k');
        %fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
        hold on
        fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
        hold on
        ylim([-.02 .16])
        xlim([-.15 .5])
        box off
        set(gca, 'TickDir', 'out')
        hline(0)
        hold off
        title([num2str(Sizes(iSize)) ] )        
        x=x+1;
        
        temp_mean1=mean(ContrastMeansLoc(:,centerIN),2,"omitnan");
        nINLoc = length(find(~isnan(ContrastMeansLoc(1,centerIN))));
        temp_se1 = std(ContrastMeansLoc(:,centerIN),[],2,"omitnan")/sqrt(nINLoc);


        temp_mean2=mean(ContrastMeansStat(:,centerIN),2,"omitnan");
        nINStat = length(find(~isnan(ContrastMeansStat(1,centerIN))));
        temp_se2 = std(ContrastMeansStat(:,centerIN),[],2,"omitnan")/sqrt(nINStat);
        
        subplot(2,nSizes,y)
        shadedErrorBar(t(:),temp_mean1,temp_se1,'m');
        hold on
        shadedErrorBar(t(:),temp_mean2,temp_se2,'k');
        %fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
        hold on
        fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
        hold on
        ylim([-.02 .16])
        xlim([-.15 .5])
        box off
        set(gca, 'TickDir', 'out')
        hline(0)
        hold off
        title([num2str(Sizes(iSize)) ] )        
        y=y+1;

        
clear temp_mean1 temp_trials1 temp_se1 temp_mean2 temp_trials2 temp_se2
    end

x0=1;
y0=1;
width=8;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
  
%sgtitle(['Running, ~', num2str(round(mean(mean(locCountPyr)))),' Pyr cells, ~',num2str(round(mean(mean(locCountSST)))),' SST cells, averaged over contrast'])
print('Running_avrgOverContrast.pdf', '-dpdf');

%% plot running averaging over size


x=1;
y=nCons+1;
figure;
    for iCon = 1:nCons %loop through the Cons
        SizeMeansLoc=nan(nOn+nOff,nCells);
        SizeMeansStat=nan(nOn+nOff,nCells);

        for iCell = 1:nCells
        
            SizeMeansLoc(:,iCell) = mean(TCs_loc(:,iSize,iCon,iCell),2,"omitnan");
            SizeMeansStat(:,iCell) = mean(TCs_stat(:,iSize,iCon,iCell),2,"omitnan");

        end

    
        temp_mean1=mean(SizeMeansLoc(:,centerPyr),2,"omitnan");
        nPyrLoc = length(find(~isnan(SizeMeansLoc(1,centerPyr))));
        temp_se1 = std(SizeMeansLoc(:,centerPyr),[],2,"omitnan")/sqrt(nPyrLoc);


        temp_mean2=mean(SizeMeansStat(:,centerPyr),2,"omitnan");
        nPyrStat = length(find(~isnan(SizeMeansStat(1,centerPyr))));
        temp_se2 = std(SizeMeansStat(:,centerPyr),[],2,"omitnan")/sqrt(nPyrStat);


        subplot(2,nCons,x)
        shadedErrorBar(t(:),temp_mean1,temp_se1,'m');
        hold on
        shadedErrorBar(t(:),temp_mean2,temp_se2,'k');
        %fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
        hold on
        fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
        hold on
        ylim([-.02 .16])
        xlim([-.15 .5])
        box off
        set(gca, 'TickDir', 'out')
        hline(0)
        hold off
        title([num2str(Cons(iCon)) ] )        
        x=x+1;
        


        temp_mean1=mean(SizeMeansLoc(:,centerIN),2,"omitnan");
        nINLoc = length(find(~isnan(SizeMeansLoc(1,centerIN))));
        temp_se1 = std(SizeMeansLoc(:,centerIN),[],2,"omitnan")/sqrt(nINLoc);


        temp_mean2=mean(SizeMeansStat(:,centerIN),2,"omitnan");
        nINStat = length(find(~isnan(SizeMeansStat(1,centerIN))));
        temp_se2 = std(SizeMeansStat(:,centerIN),[],2,"omitnan")/sqrt(nINStat);
        
        subplot(2,nCons,y)
        shadedErrorBar(t(:),temp_mean1,temp_se1,'m');
        hold on
        shadedErrorBar(t(:),temp_mean2,temp_se2,'k');
        %fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
        hold on
        fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
        hold on
        ylim([-.02 .16])
        xlim([-.15 .5])
        box off
        set(gca, 'TickDir', 'out')
        hline(0)
        hold off
        title([num2str(Cons(iCon)) ] )        
        y=y+1;

        
clear temp_mean1 temp_trials1 temp_se1 temp_mean2 temp_trials2 temp_se2
    end

x0=1;
y0=1;
width=8;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
  


print('Running_avrgOverSize.pdf', '-dpdf');
%% size curve at each contrast based on peak

figure;
subplot(1,2,1)
for iCon = 1:nCons
        temp_mean1 = mean(peak_trough(:,iCon,centerPyr,1),3,"omitnan");
        temp_se1 = std(peak_trough(:,iCon,centerPyr,1),[],3,"omitnan")/sqrt(length(centerPyr));

        errorbar(Sizes,temp_mean1,temp_se1);
        hold on
end

lgd=legend(string(Cons))
title(lgd,'Contrast')
title("Size tuning by contrast, Pyr cells")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('df/f peak')
ylim([0 .13])




subplot(1,2,2)
for iCon = 1:nCons
        temp_mean1 = mean(peak_trough(:,iCon,centerIN,1),3,"omitnan");
        temp_se1 = std(peak_trough(:,iCon,centerIN,1),[],3,"omitnan")/sqrt(length(centerIN));

        errorbar(Sizes,temp_mean1,temp_se1);
        hold on
end

title("Size tuning by contrast, SST cells")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('df/f peak')
ylim([0 .13])


x0=5;
y0=5;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
sgtitle('Stationary')


print('SizeTuning_byContrast_peak.pdf', '-dpdf');

%% size curve at each contrast based on mean
% cm = winter;
% "color",cm(iCon,:)


figure;
subplot(1,2,1)
for iCon = 1:nCons
        temp_mean1 = mean(resp_means(:,iCon,centerPyr,1),3,"omitnan");
        temp_se1 = std(resp_means(:,iCon,centerPyr,1),[],3,"omitnan")/sqrt(length(centerPyr));

        errorbar(Sizes,temp_mean1,temp_se1);
        hold on
end

lgd=legend(string(Cons))
title(lgd,'Contrast')
title("Size tuning by contrast, Pyr cells")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('df/f mean')
ylim([0 .06])




subplot(1,2,2)
for iCon = 1:nCons
        temp_mean1 = mean(resp_means(:,iCon,centerIN,1),3,"omitnan");
        temp_se1 = std(resp_means(:,iCon,centerIN,1),[],3,"omitnan")/sqrt(length(centerIN));

        errorbar(Sizes,temp_mean1,temp_se1);
        hold on
end
lgd=legend(string(Cons))
title(lgd,'Contrast')
title("Size tuning by contrast, SST cells")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('df/f mean')
ylim([0 .06])

sgtitle('Stationary')
x0=5;
y0=5;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
print('SizeTuning_byContrast_means.pdf', '-dpdf');

%% Contrast tuning by size, based on peak
figure;
subplot(1,2,1)
for iSize = 1:nSizes
        temp_mean1 = mean(peak_trough(iSize,:,centerPyr,1),3,"omitnan");
        temp_se1 = std(peak_trough(iSize,:,centerPyr,1),[],3,"omitnan")/sqrt(length(centerPyr));
        errorbar(Cons,temp_mean1,temp_se1);
        hold on
end
lgd=legend(string(Sizes),'location','best');
title(lgd,'Size')
title("Contrast tuning by contrast, Pyr cells")
box off
set(gca, 'TickDir', 'out')
xticks(Cons)
xlabel('Contrast')
ylabel('df/f peak')
ylim([0 .13])

subplot(1,2,2)
for iSize = 1:nSizes
        temp_mean1 = mean(peak_trough(iSize,:,centerIN,1),3,"omitnan");
        temp_se1 = std(peak_trough(iSize,:,centerIN,1),[],3,"omitnan")/sqrt(length(centerIN));

        errorbar(Cons,temp_mean1,temp_se1);
        hold on
end
lgd=legend(string(Sizes),'location','best');
title(lgd,'Size')
title("Contrast tuning by contrast, SST cells")
box off
set(gca, 'TickDir', 'out')
xticks(Cons)
xlabel('Contrast')
ylabel('df/f peak')
ylim([0 .13])


x0=5;
y0=5;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
print('ContrastTuning_bySize.pdf', '-dpdf');

%% TC matrix separated by depth



TCs_stat_normalized = nan(size(TCs_stat));
for iCell = 1:nCells
    normMax = peak_trough(1,4,iCell,1);
    TCs_stat_normalized(:,:,:,iCell)=TCs_stat(:,:,:,iCell)/normMax;
end

[n n2] = subplotn(nSizes*nCons);
x=1;
figure;
    for iSize = 1:nSizes %loop through the sizes
        
        for iCon = 1:nCons
        
        temp_mean1 = mean(TCs_stat_normalized(:,iSize,iCon,PyrDepth1),4,"omitnan");
        temp_se1 = std(TCs_stat_normalized(:,iSize,iCon,PyrDepth1),[],4,"omitnan")/sqrt(length(PyrDepth1));


        temp_mean2 = mean(TCs_stat_normalized(:,iSize,iCon,PyrDepth2),4,"omitnan");
        temp_se2 = std(TCs_stat_normalized(:,iSize,iCon,PyrDepth2),[],4,"omitnan")/sqrt(length(PyrDepth2));


        temp_mean3 = mean(TCs_stat_normalized(:,iSize,iCon,PyrDepth3),4,"omitnan");
        temp_se3 = std(TCs_stat_normalized(:,iSize,iCon,PyrDepth3),[],4,"omitnan")/sqrt(length(PyrDepth3));

        temp_mean4 = mean(TCs_stat_normalized(:,iSize,iCon,PyrDepth4),4,"omitnan");
        temp_se4 = std(TCs_stat_normalized(:,iSize,iCon,PyrDepth4),[],4,"omitnan")/sqrt(length(PyrDepth4));

        subplot(n,n2,x)

        shadedErrorBar(t(:),temp_mean1,temp_se1);
        hold on
        shadedErrorBar(t(:),temp_mean2,temp_se2,'b');
        hold on
        shadedErrorBar(t(:),temp_mean3,temp_se3,'m');
        hold on
        shadedErrorBar(t(:),temp_mean4,temp_se4,'r');
        hold on
        alpha(.5)
        %fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
        hold on
        fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
        hold on
        ylim([-.2 1])
        xlim([-.15 .5])  
        box off
        set(gca, 'TickDir', 'out')
        hline(0)
        hold off
        
        title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
        x=x+1;
        end
        
        clear temp_mean1 temp_trials1 temp_se1 temp_mean2 temp_trials2 temp_se2 temp_mean3 temsp_se3 temp_mean4 temp_se4

    end


x0=1;
y0=1;
width=7;
height=11;
set(gcf,'units','inches','position',[x0,y0,width,height])

%sgtitle(['Pyr by depth ', num2str(length(PyrDepth1)),' deep cells, ',num2str(length(PyrDepth2)),' shallow cells'])

 print('Pyr_depth_matrix2.pdf', '-dpdf');
%% peak time by size and depth, at 40% and 80% contrast

% 
% peak_time = nan(size(peak_time(:,4,:,1)));
% for iCell = 1:nCells
%     normMax = peak_time(1,4,iCell);
%     peak_time(:,4,iCell)=peak_time(:,4,iCell,1)/normMax;
% end


for iCon = 3:4
figure;

subplot(1,2,1)
     
temp_mean1 = mean(peak_time(:,iCon,PyrDepth1,1),3,"omitnan");
temp_se1 = std(peak_time(:,iCon,PyrDepth1,1),[],3,"omitnan")/sqrt(length(PyrDepth1));

temp_mean2 = mean(peak_time(:,iCon,PyrDepth2,1),3,"omitnan");
temp_se2 = std(peak_time(:,iCon,PyrDepth2,1),[],3,"omitnan")/sqrt(length(PyrDepth2));

temp_mean3 = mean(peak_time(:,iCon,PyrDepth3,1),3,"omitnan");
temp_se3 = std(peak_time(:,iCon,PyrDepth3,1),[],3,"omitnan")/sqrt(length(PyrDepth3));

temp_mean4 = mean(peak_time(:,iCon,PyrDepth4,1),3,"omitnan");
temp_se4 = std(peak_time(:,iCon,PyrDepth4,1),[],3,"omitnan")/sqrt(length(PyrDepth4));

errorbar(Sizes,temp_mean1,temp_se1,'k');
hold on
errorbar(Sizes,temp_mean2,temp_se2,'b');
hold on
errorbar(Sizes,temp_mean3,temp_se3,'m');
hold on
errorbar(Sizes,temp_mean4,temp_se4,'r');
legend('144-175','175-225','225-275','275-305')
title("Pyr")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('df/f mean')
%ylim([0 .1])


subplot(1,2,2)

temp_mean1 = mean(peak_time(:,iCon,INdepth1,1),3,"omitnan");
temp_se1 = std(peak_time(:,iCon,INdepth1,1),[],3,"omitnan")/sqrt(length(INdepth1));

temp_mean2 = mean(peak_time(:,iCon,INdepth2,1),3,"omitnan");
temp_se2 = std(peak_time(:,iCon,INdepth2,1),[],3,"omitnan")/sqrt(length(INdepth2));

temp_mean3 = mean(peak_time(:,iCon,INdepth3,1),3,"omitnan");
temp_se3 = std(peak_time(:,iCon,INdepth3,1),[],3,"omitnan")/sqrt(length(INdepth3));

temp_mean4 = mean(peak_time(:,iCon,INdepth4,1),3,"omitnan");
temp_se4 = std(peak_time(:,iCon,INdepth4,1),[],3,"omitnan")/sqrt(length(INdepth4));

errorbar(Sizes,temp_mean1,temp_se1,'k');
hold on
errorbar(Sizes,temp_mean2,temp_se2,'b');
hold on
errorbar(Sizes,temp_mean3,temp_se3,'m');
hold on
errorbar(Sizes,temp_mean4,temp_se4,'r');
title("SST")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('s')
%ylim([0 .1])

sgtitle(['Peak time by size and depth at ',num2str(Cons(iCon)), ' stationary'])

x0=1;
y0=1;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(['PeakTime_byDepth_',num2str(Cons(iCon)), '.pdf'], '-dpdf');

end
%% Trough time by size and depth at 40% and 80% contrast

% 
% peak_time = nan(size(peak_time(:,4,:,2)));
% for iCell = 1:nCells
%     normMax = peak_time(1,4,iCell);
%     peak_time(:,4,iCell)=peak_time(:,4,iCell,2)/normMax;
% end


for iCon = 3:4
figure;

subplot(1,2,1)
     
temp_mean1 = mean(peak_time(:,iCon,PyrDepth1,2),3,"omitnan");
temp_se1 = std(peak_time(:,iCon,PyrDepth1,2),[],3,"omitnan")/sqrt(length(PyrDepth1));

temp_mean2 = mean(peak_time(:,iCon,PyrDepth2,2),3,"omitnan");
temp_se2 = std(peak_time(:,iCon,PyrDepth2,2),[],3,"omitnan")/sqrt(length(PyrDepth2));

temp_mean3 = mean(peak_time(:,iCon,PyrDepth3,2),3,"omitnan");
temp_se3 = std(peak_time(:,iCon,PyrDepth3,2),[],3,"omitnan")/sqrt(length(PyrDepth3));

temp_mean4 = mean(peak_time(:,iCon,PyrDepth4,2),3,"omitnan");
temp_se4 = std(peak_time(:,iCon,PyrDepth4,2),[],3,"omitnan")/sqrt(length(PyrDepth4));

errorbar(Sizes,temp_mean1,temp_se1,'k');
hold on
errorbar(Sizes,temp_mean2,temp_se2,'b');
hold on
errorbar(Sizes,temp_mean3,temp_se3,'m');
hold on
errorbar(Sizes,temp_mean4,temp_se4,'r');
legend('144-175','175-225','225-275','275-305')
title("Pyr")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('df/f mean')
ylim([0.22 .3])


subplot(1,2,2)

temp_mean1 = mean(peak_time(:,iCon,INdepth1,2),3,"omitnan");
temp_se1 = std(peak_time(:,iCon,INdepth1,2),[],3,"omitnan")/sqrt(length(INdepth1));

temp_mean2 = mean(peak_time(:,iCon,INdepth2,2),3,"omitnan");
temp_se2 = std(peak_time(:,iCon,INdepth2,2),[],3,"omitnan")/sqrt(length(INdepth2));

temp_mean3 = mean(peak_time(:,iCon,INdepth3,2),3,"omitnan");
temp_se3 = std(peak_time(:,iCon,INdepth3,2),[],3,"omitnan")/sqrt(length(INdepth3));

temp_mean4 = mean(peak_time(:,iCon,INdepth4,2),3,"omitnan");
temp_se4 = std(peak_time(:,iCon,INdepth4,2),[],3,"omitnan")/sqrt(length(INdepth4));

errorbar(Sizes,temp_mean1,temp_se1,'k');
hold on
errorbar(Sizes,temp_mean2,temp_se2,'b');
hold on
errorbar(Sizes,temp_mean3,temp_se3,'m');
hold on
errorbar(Sizes,temp_mean4,temp_se4,'r');
title("SST")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('s')
ylim([0.22 .3])

sgtitle(['Trough time by size and depth at ',num2str(Cons(iCon)), ' stationary'])

x0=1;
y0=1;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(['TroughTime_byDepth_',num2str(Cons(iCon)), '.pdf'], '-dpdf');

end

%% Size tuning based on peak, separated by depth, at 80% contrast


% 
% peak_trough = nan(size(peak_trough(:,:,:,1)));
% for iCell = 1:nCells
%     normMax = peak_trough(1,4,iCell);
%     peak_trough(:,4,iCell)=peak_trough(:,4,iCell,1)/normMax;
% end


for iCon = 3:4
figure;

subplot(1,2,1)
     
temp_mean1 = mean(resp_means(:,iCon,PyrDepth1,1),3,"omitnan");
temp_se1 = std(resp_means(:,iCon,PyrDepth1,1),[],3,"omitnan")/sqrt(length(PyrDepth1));

temp_mean2 = mean(resp_means(:,iCon,PyrDepth2,1),3,"omitnan");
temp_se2 = std(resp_means(:,iCon,PyrDepth2,1),[],3,"omitnan")/sqrt(length(PyrDepth2));

temp_mean3 = mean(resp_means(:,iCon,PyrDepth3,1),3,"omitnan");
temp_se3 = std(resp_means(:,iCon,PyrDepth3,1),[],3,"omitnan")/sqrt(length(PyrDepth3));

temp_mean4 = mean(resp_means(:,iCon,PyrDepth4,1),3,"omitnan");
temp_se4 = std(resp_means(:,iCon,PyrDepth4,1),[],3,"omitnan")/sqrt(length(PyrDepth4));

errorbar(Sizes,temp_mean1,temp_se1,'k');
hold on
errorbar(Sizes,temp_mean2,temp_se2,'b');
hold on
errorbar(Sizes,temp_mean3,temp_se3,'m');
hold on
errorbar(Sizes,temp_mean4,temp_se4,'r');
legend('144-175','175-225','225-275','275-305')
title("Pyr")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('df/f peak')
ylim([0 .09])


subplot(1,2,2)

temp_mean1 = mean(resp_means(:,iCon,INdepth1,1),3,"omitnan");
temp_se1 = std(resp_means(:,iCon,INdepth1,1),[],3,"omitnan")/sqrt(length(INdepth1));

temp_mean2 = mean(resp_means(:,iCon,INdepth2,1),3,"omitnan");
temp_se2 = std(resp_means(:,iCon,INdepth2,1),[],3,"omitnan")/sqrt(length(INdepth2));

temp_mean3 = mean(resp_means(:,iCon,INdepth3,1),3,"omitnan");
temp_se3 = std(resp_means(:,iCon,INdepth3,1),[],3,"omitnan")/sqrt(length(INdepth3));

temp_mean4 = mean(resp_means(:,iCon,INdepth4,1),3,"omitnan");
temp_se4 = std(resp_means(:,iCon,INdepth4,1),[],3,"omitnan")/sqrt(length(INdepth4));

errorbar(Sizes,temp_mean1,temp_se1,'k');
hold on
errorbar(Sizes,temp_mean2,temp_se2,'b');
hold on
errorbar(Sizes,temp_mean3,temp_se3,'m');
hold on
errorbar(Sizes,temp_mean4,temp_se4,'r');
title("SST")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')

ylim([0 .09])

sgtitle(['Size tuning by depth, peak at ',num2str(Cons(iCon)), ' stationary'])

x0=1;
y0=1;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(['SizeTuning_byDepth_peak_',num2str(Cons(iCon)), '.pdf'], '-dpdf');

end



%% Size tuning based on mean, separated by depth, at 40% and 80% contrast

% 
% resp_means = nan(size(resp_means(:,:,:,1)));
% for iCell = 1:nCells
%     normMax = resp_means(1,4,iCell);
%     resp_means(:,4,iCell)=resp_means(:,4,iCell,1)/normMax;
% end

for iCon = 3:4
figure;

subplot(1,2,1)
temp_mean1 = mean(resp_means(:,iCon,PyrDepth1,1),3,"omitnan");
temp_se1 = std(resp_means(:,iCon,PyrDepth1,1),[],3,"omitnan")/sqrt(length(PyrDepth1));

temp_mean2 = mean(resp_means(:,iCon,PyrDepth2,1),3,"omitnan");
temp_se2 = std(resp_means(:,iCon,PyrDepth2,1),[],3,"omitnan")/sqrt(length(PyrDepth2));

temp_mean3 = mean(resp_means(:,iCon,PyrDepth3,1),3,"omitnan");
temp_se3 = std(resp_means(:,iCon,PyrDepth3,1),[],3,"omitnan")/sqrt(length(PyrDepth3));

temp_mean4 = mean(resp_means(:,iCon,PyrDepth4,1),3,"omitnan");
temp_se4 = std(resp_means(:,iCon,PyrDepth4,1),[],3,"omitnan")/sqrt(length(PyrDepth4));

errorbar(Sizes,temp_mean1,temp_se1,'k');
hold on
errorbar(Sizes,temp_mean2,temp_se2,'b');
hold on
errorbar(Sizes,temp_mean3,temp_se3,'m');
hold on
errorbar(Sizes,temp_mean4,temp_se4,'r');
legend('144-175','175-225','225-275','275-305')
title("Pyr")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('df/f mean')
ylim([0 .09])


subplot(1,2,2)

temp_mean1 = mean(resp_means(:,iCon,INdepth1,1),3,"omitnan");
temp_se1 = std(resp_means(:,iCon,INdepth1,1),[],3,"omitnan")/sqrt(length(INdepth1));

temp_mean2 = mean(resp_means(:,iCon,INdepth2,1),3,"omitnan");
temp_se2 = std(resp_means(:,iCon,INdepth2,1),[],3,"omitnan")/sqrt(length(INdepth2));

temp_mean3 = mean(resp_means(:,iCon,INdepth3,1),3,"omitnan");
temp_se3 = std(resp_means(:,iCon,INdepth3,1),[],3,"omitnan")/sqrt(length(INdepth3));

temp_mean4 = mean(resp_means(:,iCon,INdepth4,1),3,"omitnan");
temp_se4 = std(resp_means(:,iCon,INdepth4,1),[],3,"omitnan")/sqrt(length(INdepth4));

errorbar(Sizes,temp_mean1,temp_se1,'k');
hold on
errorbar(Sizes,temp_mean2,temp_se2,'b');
hold on
errorbar(Sizes,temp_mean3,temp_se3,'m');
hold on
errorbar(Sizes,temp_mean4,temp_se4,'r');
title("SST")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylim([0 .09])

sgtitle(['Size tuning by depth, mean at ',num2str(Cons(iCon)), ' stationary'])

x0=1;
y0=1;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(['SizeTuning_byDepth_means_',num2str(Cons(iCon)), '.pdf'], '-dpdf');

end
%% Size tuning based on integral, separated by depth, at 40% and 80% contrast

for iCon = 3:4
figure;

subplot(1,2,1)
     
temp_mean1 = mean(resp_means(:,iCon,PyrDepth1,2),3,"omitnan");
temp_se1 = std(resp_means(:,iCon,PyrDepth1,2),[],3,"omitnan")/sqrt(length(PyrDepth1));

temp_mean2 = mean(resp_means(:,iCon,PyrDepth2,2),3,"omitnan");
temp_se2 = std(resp_means(:,iCon,PyrDepth2,2),[],3,"omitnan")/sqrt(length(PyrDepth2));

temp_mean3 = mean(resp_means(:,iCon,PyrDepth3,2),3,"omitnan");
temp_se3 = std(resp_means(:,iCon,PyrDepth3,2),[],3,"omitnan")/sqrt(length(PyrDepth3));

temp_mean4 = mean(resp_means(:,iCon,PyrDepth4,2),3,"omitnan");
temp_se4 = std(resp_means(:,iCon,PyrDepth4,2),[],3,"omitnan")/sqrt(length(PyrDepth4));

errorbar(Sizes,temp_mean1,temp_se1,'k');
hold on
errorbar(Sizes,temp_mean2,temp_se2,'b');
hold on
errorbar(Sizes,temp_mean3,temp_se3,'m');
hold on
errorbar(Sizes,temp_mean4,temp_se4,'r');
legend('144-175','175-225','225-275','275-305')
title("Pyr")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('df/f integral')
%ylim([0 .1])


subplot(1,2,2)

temp_mean1 = mean(resp_means(:,iCon,INdepth1,2),3,"omitnan");
temp_se1 = std(resp_means(:,iCon,INdepth1,2),[],3,"omitnan")/sqrt(length(INdepth1));

temp_mean2 = mean(resp_means(:,iCon,INdepth2,2),3,"omitnan");
temp_se2 = std(resp_means(:,iCon,INdepth2,2),[],3,"omitnan")/sqrt(length(INdepth2));

temp_mean3 = mean(resp_means(:,iCon,INdepth3,2),3,"omitnan");
temp_se3 = std(resp_means(:,iCon,INdepth3,2),[],3,"omitnan")/sqrt(length(INdepth3));

temp_mean4 = mean(resp_means(:,iCon,INdepth4,2),3,"omitnan");
temp_se4 = std(resp_means(:,iCon,INdepth4,2),[],3,"omitnan")/sqrt(length(INdepth4));

errorbar(Sizes,temp_mean1,temp_se1,'k');
hold on
errorbar(Sizes,temp_mean2,temp_se2,'b');
hold on
errorbar(Sizes,temp_mean3,temp_se3,'m');
hold on
errorbar(Sizes,temp_mean4,temp_se4,'r');
title("SST")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
%ylim([0 .1])

sgtitle(['Size tuning by depth, integral at ',num2str(Cons(iCon)), ' stationary'])

x0=1;
y0=1;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(['SizeTuning_byDepth_integral_',num2str(Cons(iCon)), '.pdf'], '-dpdf');

end


%% full width at half max at 40% and 80% contrast 

figure;
for iCon=3:4

subplot(1,2,iCon-2)
temp_mean1 = mean(fwhm(:,iCon,centerPyr),3,"omitnan");
temp_se1 = std(fwhm(:,iCon,centerPyr),[],3,"omitnan")/sqrt(length(centerPyr));

temp_mean2 = mean(fwhm(:,iCon,centerIN),3,"omitnan");
temp_se2 = std(fwhm(:,iCon,centerIN),[],3,"omitnan")/sqrt(length(centerIN));

errorbar(Sizes,temp_mean1,temp_se1);
hold on
errorbar(Sizes,temp_mean2,temp_se2);
legend('Pyr','SST')
title(num2str(Cons(iCon)))
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('seconds')
ylim([0.06 .2])

end
sgtitle("FWHM by size, stationary")

x0=1;
y0=1;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
print('FWHM_by_size.pdf', '-dpdf');


%% full width at half max by depth

% fwhm = nan(size(fwhm));
% for iCell = 1:nCells
%     normMax = fwhm(1,1,iCell);
%     fwhm(:,:,iCell)=fwhm(:,:,iCell)/normMax;
% end

for iCon=3:4
figure;

subplot(1,2,1)
temp_mean1 = mean(fwhm(:,iCon,PyrDepth1),3,"omitnan");
temp_se1 = std(fwhm(:,iCon,PyrDepth1),[],3,"omitnan")/sqrt(length(PyrDepth1));

temp_mean2 = mean(fwhm(:,iCon,PyrDepth2),3,"omitnan");
temp_se2 = std(fwhm(:,iCon,PyrDepth2),[],3,"omitnan")/sqrt(length(PyrDepth2));

temp_mean3 = mean(fwhm(:,iCon,PyrDepth3),3,"omitnan");
temp_se3 = std(fwhm(:,iCon,PyrDepth3),[],3,"omitnan")/sqrt(length(PyrDepth3));

temp_mean4 = mean(fwhm(:,iCon,PyrDepth4),3,"omitnan");
temp_se4 = std(fwhm(:,iCon,PyrDepth4),[],3,"omitnan")/sqrt(length(PyrDepth4));

errorbar(Sizes,temp_mean1,temp_se1,'k');
hold on
errorbar(Sizes,temp_mean2,temp_se2,'b');
hold on
errorbar(Sizes,temp_mean3,temp_se3,'m');
hold on
errorbar(Sizes,temp_mean4,temp_se4,'r');
legend('144-175','175-225','225-275','275-305')
title("Pyr")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('seconds')
%ylim([0.4 1.1])


subplot(1,2,2)

temp_mean1 = mean(fwhm(:,iCon,INdepth1),3,"omitnan");
temp_se1 = std(fwhm(:,iCon,INdepth1),[],3,"omitnan")/sqrt(length(INdepth1));

temp_mean2 = mean(fwhm(:,iCon,INdepth2),3,"omitnan");
temp_se2 = std(fwhm(:,iCon,INdepth2),[],3,"omitnan")/sqrt(length(INdepth2));

temp_mean3 = mean(fwhm(:,iCon,INdepth3),3,"omitnan");
temp_se3 = std(fwhm(:,iCon,INdepth3),[],3,"omitnan")/sqrt(length(INdepth3));

temp_mean4 = mean(fwhm(:,iCon,INdepth4),3,"omitnan");
temp_se4 = std(fwhm(:,iCon,INdepth4),[],3,"omitnan")/sqrt(length(INdepth4));

errorbar(Sizes,temp_mean1,temp_se1,'k');
hold on
errorbar(Sizes,temp_mean2,temp_se2,'b');
hold on
errorbar(Sizes,temp_mean3,temp_se3,'m');
hold on
errorbar(Sizes,temp_mean4,temp_se4,'r');
title("SST")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('seconds')
%ylim([0 1.6])

sgtitle(['FWHM by size and depth, ',num2str(Cons(iCon)), ' stationary'])

x0=1;
y0=1;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(['FWHM_bySizeandDepth_',num2str(Cons(iCon)), '.pdf'], '-dpdf');


end

%% LINE PLOTS FOR RUNNING DATA
%% find peak, trough and similar metrics



%make an empty matrix for values
peak_trough_loc =nan(nSizes,nCons,nCells,2);
dip_loc =nan(nSizes,nCons,nCells);
dip_peak_ratio_loc =nan(nSizes,nCons,nCells);
peak_time_loc=nan(nSizes,nCons,nCells,4);
resp_means_loc = nan(nSizes,nCons,nCells,2);
fwhm_loc=nan(nSizes,nCons,nCells); 

for iSize = 1:nSizes %loop through the sizes    
    for iCon = 1:nCons
        for iCell = 1:nCells
           
        if ~isnan(locResp(iSize,iCon,iCell))
           thisTrace = TCs_loc(:,iSize,iCon,iCell);
           %interpolate data
           traceInterp=interp1(t,thisTrace,(t(1):0.01:t(123)));
           %traceInterp=smoothdata(traceInterp,'movmean',3);
           t2=t(1):0.01:t(123); % to get temporal values with the interpolated data
           stimStart_interp=find(t2==0);
           %find(t2==.2)

           peak=max(traceInterp(stimStart_interp:stimStart_interp+20)); %find the max value within a set window
           %currently set to 61 (stim onset) through 67, 200ms after stim
           trough=min(traceInterp(stimStart_interp+20:stimStart_interp+30));
           %search for the trough between 200 and 300 ms
           
%           (thisTrace(stimStart+6:stimStart+9));

           
           peak_time_temp = find(traceInterp==peak);

           peakBin = [peak_time_temp-5,peak_time_temp+5]; %a ~100 ms bin around the peak
           
           trough_time_temp = find(traceInterp==trough);
           troughBin = [trough_time_temp-7,trough_time_temp+7]; 

           peak_time_loc(iSize,iCon, iCell,1)=t2(peak_time_temp);
           peak_time_loc(iSize,iCon, iCell,2)= t2(trough_time_temp);

                if peak > 0
                   if iCon > 2
                       
                       halfPeak = peak/2;
                       %find the frame of the half peak
                       half_peak_frame = find(traceInterp(stimStart_interp:peak_time_temp)>halfPeak,1,'first');
                       half_peak_temp = stimStart_interp+(half_peak_frame-1); %adjust this for the frame we started on
                        %convert this to time
                       peak_time_loc(iSize,iCon, iCell,3) = t2(half_peak_temp);
                        %find the frame of the equivalent point on the decay
                       if min(traceInterp(peak_time_temp:length(t2)))<halfPeak
                           half_dacay_frame=find(traceInterp(peak_time_temp:length(t2))<halfPeak,1,'first'); 
                           %look for the half decay anywhere after the peak, until the end of the trace
                           half_dacay_temp=peak_time_temp+half_dacay_frame-1; %adjust this for the frame we started on
                           peak_time_loc(iSize,iCon, iCell,4) = t2(half_dacay_temp);
                            %find the difference, in time, between the half decay
                            %and the half peak
                           fwhm_loc(iSize,iCon,iCell)=t2(half_dacay_temp)-t2(half_peak_temp);
                       end
        
                       % [minVal, endWin]=min(abs(t2-t(peak_time_temp))); %find the value in the interpoalted t vector that is closest to the time in the oringal t vector when the peak occurs
                       % half_peak_temp = find(traceInterp(stimStart_interp:endWin)>halfPeak,1,'first');
                       % peak_time(iSize,iCon, iCell,3) = t2(stimStart_interp+(half_peak_temp-1));
        
                    %dip is finding the local decrease within the general vicinty
                    %of through
                   end
                end
            dip_temp=trough-((traceInterp(troughBin(1))+traceInterp(troughBin(2)))/2);

            dip_loc(iSize,iCon, iCell)=dip_temp;            
            dip_peak_ratio_loc(iSize,iCon, iCell)=dip_temp/peak;
           
           %figure;plot(thisTrace);hold on; vline(peak_time_temp);vline(trough_time_temp);xlim([stimStart stimStart+15]);hold off
          
           peak_trough_loc(iSize,iCon, iCell,1)=peak;
           %peak is the identified max
           peak_trough_loc(iSize,iCon, iCell,2)=trough;
            %rough is the minimum within the identified timebin, which is
            %locked to the time of the peak 

            %resp_mean is the average (:,:,:,1) or cumulative integral (:,:,:,2) response from 0 to 400 ms
           resp_means_loc(iSize,iCon, iCell,1)=mean(thisTrace(stimStart:stimStart+12));
           resp_means_loc(iSize,iCon, iCell,2)=max(cumtrapz(t(stimStart:stimStart+12),thisTrace(stimStart:stimStart+12)));


           clear thisTrace peak peakBin troughBin dip_temp trough trough_time_temp peak_time_temp
        end
        end
    end
end


%% look at the mean peak for each
mean_peak_loc=(mean(peak_time_loc(:,:,centered,1),3,"omitmissing"));
mean_halfPeak_loc=mean(peak_time_loc(:,:,centered,3),3,"omitmissing");
mean_trough_loc=(mean(peak_time_loc(:,:,centered,2),3,"omitmissing"));
mean_halfDecay_loc=(mean(peak_time_loc(:,:,centered,4),3,"omitmissing"));

[n n2] = subplotn(nSizes*nCons);
x=1;

figure;
    for iSize = 1:nSizes %loop through the sizes
        
        for iCon = 1:nCons

        temp_peak=mean_peak_loc(iSize,iCon);
        temp_halfPeak=mean_halfPeak_loc(iSize,iCon);
        temp_trough=mean_trough_loc(iSize,iCon);
        temp_halfDecay=mean_halfDecay_loc(iSize,iCon);
        
        temp_mean1 = mean(TCs_loc(:,iSize,iCon,centerPyr),4,"omitnan");
        temp_se1 = std(TCs_loc(:,iSize,iCon,centerPyr),[],4,"omitnan")/sqrt(length(centerPyr));

        temp_mean2 = mean(TCs_loc(:,iSize,iCon,centerIN),4,"omitnan");
        temp_se2 = std(TCs_loc(:,iSize,iCon,centerIN),[],4,"omitnan")/sqrt(length(centerIN));

        subplot(n,n2,x)
  
        hold on
        shadedErrorBar(t(:),temp_mean2,temp_se2,'r');
        hold on
        shadedErrorBar(t,temp_mean1,temp_se1);
        hold on
        alpha(.5)
        %fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
        hold on
        %fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
        hold on
        ylim([-.03 .14])
        xlim([-.25 .5])
        vline(temp_halfPeak,'g')
        vline(temp_peak,'b')
        vline(temp_trough,'k')
        vline(temp_halfDecay,'r')
        box off
        set(gca, 'TickDir', 'out')
        hline(0)
        hold off
        title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
        x=x+1;
        end
        
clear temp_mean1 temp_trials1 temp_se1 temp_mean2 temp_trials2 temp_se2
    end
  
x0=1;
y0=1;
width=5;
height=8;
set(gcf,'units','inches','position',[x0,y0,width,height])

sgtitle(['Running, ~', num2str(round(mean(mean(locCountPyr)))),' Pyr cells, ~',num2str(round(mean(mean(locCountSST)))),' SST cells'])
print('Running_peaks_labelled.pdf', '-dpdf');

%% size curve at each contrast based on peak
figure;
subplot(1,2,1)
for iCon = 1:nCons
        temp_mean1 = mean(peak_trough_loc(:,iCon,centerPyr,1),3,"omitnan");
        temp_se1 = std(peak_trough_loc(:,iCon,centerPyr,1),[],3,"omitnan")/sqrt(length(centerPyr));

        errorbar(Sizes,temp_mean1,temp_se1);
        hold on
end

lgd=legend(string(Cons))
title(lgd,'Contrast')
title("Size tuning by contrast, Pyr cells")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('df/f peak')
ylim([0 .16])




subplot(1,2,2)
for iCon = 1:nCons
        temp_mean1 = mean(peak_trough_loc(:,iCon,centerIN,1),3,"omitnan");
        temp_se1 = std(peak_trough_loc(:,iCon,centerIN,1),[],3,"omitnan")/sqrt(length(centerIN));

        errorbar(Sizes,temp_mean1,temp_se1);
        hold on
end

title("Size tuning by contrast, SST cells")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('df/f peak')
ylim([0 .13])


x0=5;
y0=5;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
sgtitle('Running')

print('SizeTuning_peak_Running.pdf', '-dpdf');

%% size curve at each contrast based on mean
figure;
subplot(1,2,1)
for iCon = 1:nCons
        temp_mean1 = mean(resp_means_loc(:,iCon,centerPyr,1),3,"omitnan");
        temp_se1 = std(resp_means_loc(:,iCon,centerPyr,1),[],3,"omitnan")/sqrt(length(centerPyr));

        errorbar(Sizes,temp_mean1,temp_se1);
        hold on
end

lgd=legend(string(Cons))
title(lgd,'Contrast')
title("Size tuning by contrast, Pyr cells")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('df/f mean')
ylim([-.005 .08])




subplot(1,2,2)
for iCon = 1:nCons
        temp_mean1 = mean(resp_means_loc(:,iCon,centerIN,1),3,"omitnan");
        temp_se1 = std(resp_means_loc(:,iCon,centerIN,1),[],3,"omitnan")/sqrt(length(centerIN));

        errorbar(Sizes,temp_mean1,temp_se1);
        hold on
end

title("Size tuning by contrast, SST cells")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('df/f mean')
ylim([-.005 .08])


x0=5;
y0=5;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
sgtitle('Running')

print('SizeTuning_mean_Running.pdf', '-dpdf');

%%
figure;
for iCon=3:4

subplot(1,2,iCon-2)
temp_mean1 = mean(fwhm_loc(:,iCon,centerPyr),3,"omitnan");
temp_se1 = std(fwhm_loc(:,iCon,centerPyr),[],3,"omitnan")/sqrt(length(centerPyr));

temp_mean2 = mean(fwhm_loc(:,iCon,centerIN),3,"omitnan");
temp_se2 = std(fwhm_loc(:,iCon,centerIN),[],3,"omitnan")/sqrt(length(centerIN));

errorbar(Sizes,temp_mean1,temp_se1);
hold on
errorbar(Sizes,temp_mean2,temp_se2);
legend('Pyr','SST')
title(num2str(Cons(iCon)))
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('seconds')
ylim([0.06 .25])

end
sgtitle("FWHM by size, running")

x0=1;
y0=1;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
print('FWHM_by_size_running.pdf', '-dpdf');

%% peak time traces

figure; 
subplot(1,2,1);
%peak_time(:,:,centerIN,1)
for iCon = 1:nCons
        temp_mean1 = mean(peak_time_loc(:,iCon,centerPyr,1),3,"omitnan");
        temp_se1 = std(peak_time_loc(:,iCon,centerPyr,1),[],3,"omitnan")/sqrt(length(centerPyr));

        errorbar(Sizes,temp_mean1,temp_se1);
        hold on
        
end

title('Centered Pyr cells')
ylabel('Peak time')
xlabel('Size')
xticks(Sizes)
ylim([.1 .18])
box off
set(gca, 'TickDir', 'out')
lgd=legend(string(Cons))
title(lgd,'Contrast')

subplot(1,2,2)
%peak_time(:,:,centerIN,1)
for iCon = 1:nCons
        temp_mean1 = mean(peak_time_loc(:,iCon,centerIN,1),3,"omitnan");
        temp_se1 = std(peak_time_loc(:,iCon,centerIN,1),[],3,"omitnan")/sqrt(length(centerIN));

        errorbar(Sizes,temp_mean1,temp_se1);
        hold on
end

title('Centered SST cells')
ylabel('Peak time')
xlabel('Size')
xticks(Sizes)
ylim([.1 .18])
box off
set(gca, 'TickDir', 'out')


x0=5;
y0=5;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height]);
sgtitle('Running')


print('Peak_time_bySizeandConstrast_running.pdf', '-dpdf');

%% peak time traces

figure; 
subplot(1,2,1);
%peak_time(:,:,centerIN,1)
for iCon = 1:nCons
        temp_mean1 = mean(peak_time_loc(:,iCon,centerPyr,2),3,"omitnan");
        temp_se1 = std(peak_time_loc(:,iCon,centerPyr,2),[],3,"omitnan")/sqrt(length(centerPyr));

        errorbar(Sizes,temp_mean1,temp_se1);
        hold on
        
end

title('Centered Pyr cells')
ylabel('Trough time')
xlabel('Size')
xticks(Sizes)
ylim([.23 .3])
box off
set(gca, 'TickDir', 'out')
lgd=legend(string(Cons))
title(lgd,'Contrast')

subplot(1,2,2)
%peak_time(:,:,centerIN,1)
for iCon = 1:nCons
        temp_mean1 = mean(peak_time_loc(:,iCon,centerIN,2),3,"omitnan");
        temp_se1 = std(peak_time_loc(:,iCon,centerIN,2),[],3,"omitnan")/sqrt(length(centerIN));

        errorbar(Sizes,temp_mean1,temp_se1);
        hold on
end

title('Centered SST cells')
ylabel('Trough time')
xlabel('Size')
xticks(Sizes)
ylim([.23 .3])
box off
set(gca, 'TickDir', 'out')


x0=5;
y0=5;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height]);
sgtitle('Running')


print('Trough_time_bySizeandConstrast_running.pdf', '-dpdf');


%% EXPLORATORY VISUALIZATIONS BELOW

%% make size X contrast matrix of average timecouses, seperated by cell type and condition




[n n2] = subplotn(nSizes*nCons);
x=1;
figure;
    for iSize = 1:nSizes %loop through the sizes
        
        for iCon = 1:nCons
        
        temp_mean1 = mean(TCs_stat(:,iSize,iCon,pyrCells),4,"omitnan");
        temp_se1 = std(TCs_stat(:,iSize,iCon,pyrCells),[],4,"omitnan")/sqrt(length(pyrCells));


        temp_mean2 = mean(TCs_stat(:,iSize,iCon,interNrns),4,"omitnan");
        temp_se2 = std(TCs_stat(:,iSize,iCon,interNrns),[],4,"omitnan")/sqrt(length(interNrns));

        subplot(n,n2,x)

        shadedErrorBar(t(:),temp_mean2,temp_se2,'r');
        hold on
        shadedErrorBar(t(:),temp_mean1,temp_se1);
        hold on
        alpha(.5)
        %fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
        hold on
        fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
        hold on
        ylim([-.03 .1])
        xlim([-1 1])
        box off
        set(gca, 'TickDir', 'out')
        hline(0)
        hold off
        title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
        x=x+1;
        end
        
clear temp_mean1 temp_trials1 temp_se1 temp_mean2 temp_trials2 temp_se2
    end
  
sgtitle(['Stationary, ', num2str(length(pyrCells)),' Pyr cells, ',num2str(length(interNrns)),' SST cells'])

 print('Stationary_matrix.pdf', '-dpdf');

%% for running trials
[n n2] = subplotn(nSizes*nCons);
x=1;
figure;
    for iSize = 1:nSizes %loop through the sizes
        
        for iCon = 1:nCons
        
        temp_mean1 = mean(TCs_loc(:,iSize,iCon,pyrCells),4,"omitnan");
        temp_se1 = std(TCs_loc(:,iSize,iCon,pyrCells),[],4,"omitnan")/sqrt(length(pyrCells));


        temp_mean2 = mean(TCs_loc(:,iSize,iCon,interNrns),4,"omitnan");
        temp_se2 = std(TCs_loc(:,iSize,iCon,interNrns),[],4,"omitnan")/sqrt(length(interNrns));

        subplot(n,n2,x)

        shadedErrorBar(t(:),temp_mean2,temp_se2,'r');
        hold on
        shadedErrorBar(t(:),temp_mean1,temp_se1);
        hold on
        
        %fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
        hold on
        fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
        hold on
        ylim([-.03 .15])
        xlim([-1 1])
        box off
        set(gca, 'TickDir', 'out')
        hline(0)
        hold off
        title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
        x=x+1;
        end
        
clear temp_mean1 temp_trials1 temp_se1 temp_mean2 temp_trials2 temp_se2
    end
  
sgtitle(['Running, ', num2str(length(pyrCells)),' Pyr cells, ',num2str(length(interNrns)),' SST cells'])
print('Running_matrix.pdf', '-dpdf');




% 
% 
% [n n2] = subplotn(nSizes*nCons);
% x=1;
% figure;
%     for iSize = 1:nSizes %loop through the sizes
% 
%         for iCon = 1:nCons
% 
%         temp_mean1 = mean(TCs_stat(:,iSize,iCon,INdepth1),4,"omitnan");
%         temp_se1 = std(TCs_stat(:,iSize,iCon,INdepth1),[],4,"omitnan")/sqrt(length(INdepth1));
% 
% 
%         temp_mean2 = mean(TCs_stat(:,iSize,iCon,INdepth2),4,"omitnan");
%         temp_se2 = std(TCs_stat(:,iSize,iCon,INdepth2),[],4,"omitnan")/sqrt(length(INdepth2));
% 
%         subplot(n,n2,x)
% 
%         shadedErrorBar(t(:),temp_mean2,temp_se2,'b');
%         hold on
%         shadedErrorBar(t(:),temp_mean1,temp_se1);
%         hold on
%         alpha(.5)
%         %fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
%         hold on
%         fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
%         hold on
%         ylim([-.03 .125])
%         xlim([-.15 .5])  
%         box off
%         set(gca, 'TickDir', 'out')
%         hline(0)
%         hold off
%         title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
%         x=x+1;
%         end
% 
% clear temp_mean1 temp_trials1 temp_se1 temp_mean2 temp_trials2 temp_se2
%     end
% 
% sgtitle(['SST by depth, ', num2str(length(INdepth1)),' deep cells, ',num2str(length(INdepth2)),' shallow cells (blue)'])
% 
% x0=1;
% y0=1;
% width=7;
% height=11;
% set(gcf,'units','inches','position',[x0,y0,width,height])
% 
%  print('SST_depth_matrix.pdf', '-dpdf');
 
 %% peak time traces by depth +/- 200 microns

PyrDeep=(intersect(centerPyr,deepCells));
INdeep=(intersect(centerIN,deepCells));

shallowCells = find(depth<200);
PyrShallow=(intersect(centerPyr,shallowCells));
INshallow=(intersect(centerIN,shallowCells));


figure; 
subplot(2,2,1)
%peak_time(:,:,centerIN,1)
for iCon = 1:nCons
        temp_mean1 = mean(peak_time(:,iCon,PyrDeep,1),3,"omitnan");
        temp_se1 = std(peak_time(:,iCon,PyrDeep,1),[],3,"omitnan")/sqrt(length(PyrDeep));

        errorbar(Sizes,temp_mean1,temp_se1);
        hold on
        
end

title('Deep Pyr cells')
ylabel('Peak time')
xlabel('Size')
xticks(Sizes)
ylim([.08 .18])
box off
set(gca, 'TickDir', 'out')
lgd=legend(string(Cons))
title(lgd,'Contrast')

subplot(2,2,2)
%peak_time(:,:,centerIN,1)
for iCon = 1:nCons
        temp_mean1 = mean(peak_time(:,iCon,INdeep,1),3,"omitnan");
        temp_se1 = std(peak_time(:,iCon,INdeep,1),[],3,"omitnan")/sqrt(length(INdeep));

        errorbar(Sizes,temp_mean1,temp_se1);
        hold on
end

title('Deep SST cells')
ylabel('Peak time')
xlabel('Size')
xticks(Sizes)
ylim([.08 .18])
box off
set(gca, 'TickDir', 'out')

subplot(2,2,3)
%peak_time(:,:,centerIN,1)
for iCon = 1:nCons
        temp_mean1 = mean(peak_time(:,iCon,PyrShallow,1),3,"omitnan");
        temp_se1 = std(peak_time(:,iCon,PyrShallow,1),[],3,"omitnan")/sqrt(length(PyrShallow));

        errorbar(Sizes,temp_mean1,temp_se1);
        hold on
        
end

title('Shallow Pyr cells')
ylabel('Peak time')
xlabel('Size')
xticks(Sizes)
ylim([.08 .18])
box off
set(gca, 'TickDir', 'out')

subplot(2,2,4)
%peak_time(:,:,centerIN,1)
for iCon = 1:nCons
        temp_mean1 = mean(peak_time(:,iCon,INshallow,1),3,"omitnan");
        temp_se1 = std(peak_time(:,iCon,INshallow,1),[],3,"omitnan")/sqrt(length(INshallow));

        errorbar(Sizes,temp_mean1,temp_se1);
        hold on
end

title('Shallow SST cells')
ylabel('Peak time')
xlabel('Size')
xticks(Sizes)
ylim([.08 .18])
box off
set(gca, 'TickDir', 'out')


x0=0;
y0=0;
width=7;
height=6;
set(gcf,'units','inches','position',[x0,y0,width,height])

print('Peak_time_by_depth.pdf', '-dpdf');


%% TC matrix separated distance of RF from stim
DistCutoffs=[0,10,20,100];

    for iDist = 1:(length(DistCutoffs)-1) %loop through the sizes
          minDist=DistCutoffs(iDist);
          maxDist=DistCutoffs(iDist+1);
         
        PyrIndsTemp=intersect(pyrCells,find(dists_concat>minDist & dists_concat<maxDist));
        INIndsTemp=intersect(interNrns,find(dists_concat>minDist & dists_concat<maxDist));

        figure;
        [n n2] = subplotn(nSizes*nCons);
        x=1;
        for iSize = 1:nSizes %loop through the sizes

                for iCon = 1:nCons

                    temp_mean1 = mean(TCs_stat(:,iSize,iCon,PyrIndsTemp),4,"omitnan");
                    temp_se1 = std(TCs_stat(:,iSize,iCon,PyrIndsTemp),[],4,"omitnan")/sqrt(length(PyrIndsTemp));
                
                    temp_mean2 = mean(TCs_stat(:,iSize,iCon,INIndsTemp),4,"omitnan");
                    temp_se2 = std(TCs_stat(:,iSize,iCon,INIndsTemp),[],4,"omitnan")/sqrt(length(INIndsTemp));
        
                    subplot(n,n2,x)
            
                    shadedErrorBar(t(:),temp_mean2,temp_se2,'r');
                    hold on
                    shadedErrorBar(t(:),temp_mean1,temp_se1);
                    hold on
                    alpha(.5)
                    %fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
                    hold on
                    fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
                    hold on
                    ylim([-.03 .1])
                    xlim([-1 1])
                    
                    hline(0)
                    hold off
                    title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
                    x=x+1;
                    sgtitle([num2str(minDist) ' to ' num2str(maxDist), ', ', num2str(length(PyrIndsTemp)),' Pyr cells, ',num2str(length(INIndsTemp)),' SST cells'])
                    clear temp_se2 temp_mean2 temp_se1 temp_mean1
            end
        

    end
 
   print(fullfile(['distance',num2str(iDist),'_matrix.pdf']), '-dpdf');

    end

%% TC matrix sparated by preferred Ori

 oris=unique(prefOri_concat);

    for iOri = 1:length(oris) %loop through the sizes

        oriInds = find(prefOri_concat == oris(iOri));
        PyrIndsTemp=find(intersect(pyrCells,oriInds));
        INIndsTemp=find(intersect(interNrns,oriInds));

        figure;
        [n n2] = subplotn(nSizes*nCons);
        
        x=1;
        for iSize = 1:nSizes %loop through the sizes
                
                for iCon = 1:nCons
                                        temp_mean1 = mean(TCs_stat(:,iSize,iCon,PyrIndsTemp),4,"omitnan");
                    temp_se1 = std(TCs_stat(:,iSize,iCon,PyrIndsTemp),[],4,"omitnan")/sqrt(length(PyrIndsTemp));
                
                    temp_mean2 = mean(TCs_stat(:,iSize,iCon,INIndsTemp),4,"omitnan");
                    temp_se2 = std(TCs_stat(:,iSize,iCon,INIndsTemp),[],4,"omitnan")/sqrt(length(INIndsTemp));
        
                    subplot(n,n2,x)
            
                    shadedErrorBar(t(:),temp_mean2,temp_se2,'r');
                    hold on
                    shadedErrorBar(t(:),temp_mean1,temp_se1);
                    hold on
                    
                    alpha(.5)
                    %fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
                    hold on
                    fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
                    hold on
                    ylim([-.03 .15])
                    xlim([-1 1])

                    
                    hline(0)
                    hold off
                    title([num2str(Sizes(iSize)) ' X ' num2str(Cons(iCon))] )        
                    x=x+1;
                   sgtitle(['pref ori ', num2str(oris(iOri)) ', ', num2str(length(PyrIndsTemp)),' Pyr cells, ',num2str(length(INIndsTemp)),' SST cells'] ) 
              clear temp_se2 temp_mean2 temp_se1 temp_mean1
              
                end

        end
        
print(['prefOri',num2str(iOri),'_matrix.pdf'], '-dpdf');
    end

%% responsive at each stim condition, not necc. centered

nSST_resp_by_cond = nan(nCons,nSizes);
h_any = logical(sum(sum(h_concat,3),2));
nSST_resp_any = sum(h_any(interNrns));


t = 1:double(nOn+nOff);
t=(t-double(stimStart))/frame_rate;
[n n2] = subplotn(nSizes*nCons);
x=1;
fig=figure;
for iCon = 1:nCons
    for iSize = 1:nSizes %loop through the sizes
        
        
        responsiveTheseTrials = find(h_concat(:,iSize,iCon));
        thesePyr=intersect(pyrCells,responsiveTheseTrials);
        theseIN=intersect(interNrns,responsiveTheseTrials);
        nSST_resp_by_cond(iCon,iSize)=length(theseIN);

        temp_mean1 = mean(TCs_stat(:,iSize,iCon,thesePyr),4,"omitnan");
        temp_se1 = std(TCs_stat(:,iSize,iCon,thesePyr),[],4,"omitnan")/sqrt(length(thesePyr));


        temp_mean2 = mean(TCs_stat(:,iSize,iCon,theseIN),4,"omitnan");
        temp_se2 = std(TCs_stat(:,iSize,iCon,theseIN),[],4,"omitnan")/sqrt(length(theseIN));



        subplot(n2,n,x)

        shadedErrorBar(t(:),temp_mean2,temp_se2,'b');
        hold on
        shadedErrorBar(t(:),temp_mean1,temp_se1);
        hold on
        alpha(.5)
        %fill([.2 .2 .4 .4],[-.1 .15 .15 -.1],'b',FaceAlpha = 0.25,LineStyle='none')
        hold on
        fill([0 0 .1 .1],[-.015 -.01 -.01 -.015],'r',FaceAlpha = 0.25,LineStyle='none')
        hold on
        ylim([-.03 .1])
        xlim([-.15 .5])
        box off
        set(gca, 'TickDir', 'out')
        hline(0)
        hold off
        title([{num2str(length(thesePyr)) ' Pyr cells'} {num2str(length(theseIN)) ' SST cells'}] )
        %title([num2str(Sizes(iSize)),' contrast ',num2str(Cons(iCon))])
        x=x+1;
        end
        
clear temp_mean1 temp_trials1 temp_se1 temp_mean2 temp_trials2 temp_se2
    end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Contrast');
xlabel(han,'Size');
x0=1;
y0=1;
width=7;
height=11;
set(gcf,'units','inches','position',[x0,y0,width,height])

sgtitle(['Stationary'])

print('TCs_respAtEach', '-dpdf');

nSST_PCT_by_cond = nSST_resp_by_cond./nSST_resp_any;
figure; heatmap(nSST_PCT_by_cond);
xlabel('Size')
ylabel('Contrast')
ax = gca;
ax.XData = [7.5, 15, 30, 60, 120]
ax.YData =[0.1, 0.2, 0.4, 0.8]
print('Prct_respAtEach', '-dpdf');
%% rise time, decay time, and FWHM for cells responsive at each size, 80% contrast
temp_mean1 = nan(nCons,nSizes,3);
temp_se1 = nan(nCons,nSizes,3);
temp_mean2 = nan(nCons,nSizes,3);
temp_se2 = nan(nCons,nSizes,3);

for iCon = 1:nCons
    for iSize = 1:nSizes %loop through the sizes
        
        
        responsiveTheseTrials = find(h_concat(:,iSize,iCon));
        thesePyr=intersect(pyrCells,responsiveTheseTrials);
        theseIN=intersect(interNrns,responsiveTheseTrials);

        temp_mean1(iCon,iSize,1:2) = mean(peak_time(iSize,iCon,thesePyr,3:4),3,"omitnan");
        temp_se1(iCon,iSize,1:2) = (std(peak_time(iSize,iCon,thesePyr,3:4),[],3,"omitnan"))./length(thesePyr);
        temp_mean1(iCon,iSize,3)=mean(fwhm(iSize,iCon,thesePyr),3,"omitnan");
        temp_se1(iCon,iSize,3)=(std(fwhm(iSize,iCon,thesePyr),[],3,"omitnan"))./length(thesePyr);


        temp_mean2(iCon,iSize,1:2) =  mean(peak_time(iSize,iCon,theseIN,3:4),3,"omitnan");
        temp_se2(iCon,iSize,1:2) = (std(peak_time(iSize,iCon,theseIN,3:4),[],3,"omitnan"))./length(theseIN);
        temp_mean2(iCon,iSize,3)=mean(fwhm(iSize,iCon,theseIN),3,"omitnan");
        temp_se2(iCon,iSize,3)=(std(fwhm(iSize,iCon,theseIN),[],3,"omitnan"))./length(theseIN);


    end
end

figure;
subplot(1,3,1)
errorbar(Sizes,temp_mean2(4,:,1),temp_se2(4,:,1),'Color',	"#4e701f");
hold on
errorbar(Sizes,temp_mean1(4,:,1),temp_se1(4,:,1),'k');
%ylim([0,.11])
xticks(Sizes)
xlabel('Size')
ylabel('Seconds')
title('Rise time')
set(gca, 'TickDir', 'out')
box off

subplot(1,3,2)
errorbar(Sizes,temp_mean2(4,:,2),temp_se2(4,:,2),'Color',	"#4e701f");
hold on
errorbar(Sizes,temp_mean1(4,:,2),temp_se1(4,:,2),'k');
%ylim([0,.28])
xticks(Sizes)
xlabel('Size')
title('Decay time')
set(gca, 'TickDir', 'out')
box off

subplot(1,3,3)
errorbar(Sizes,temp_mean2(4,:,3),temp_se2(4,:,3),'Color',	"#4e701f");
hold on
errorbar(Sizes,temp_mean1(4,:,3),temp_se1(4,:,3),'k');
%ylim([0,.28])
xticks(Sizes)
xlabel('Size')
title('FWHM')
set(gca, 'TickDir', 'out')
box off



x0=5;
y0=5;
width=7;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])

print('dynamics_respAtEach', '-dpdf');
%% size tuning for only responsive cells are each condition
figure;
subplot(1,2,1)
for iCon = 1:nCons
    sizeMeans=[];
    sizeSE=[];
    for iSize = 1:nSizes
        responsiveTheseTrials = find(h_concat(:,iSize,iCon));
        thesePyr=intersect(pyrCells,responsiveTheseTrials);

        temp_mean = mean(peak_trough(iSize,iCon,thesePyr,1),3,"omitnan");
        temp_se = std(peak_trough(iSize,iCon,thesePyr,1),[],3,"omitnan")/sqrt(length(thesePyr));

        sizeMeans=[sizeMeans,temp_mean];
        sizeSE=[sizeSE,temp_se];
    end
    errorbar(Sizes,sizeMeans,sizeSE);
    hold on
end

lgd=legend(string(Cons))
title(lgd,'Contrast')
title("Size tuning by contrast, Pyr cells")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('df/f peak')
ylim([0 .13])

hold off


subplot(1,2,2)
for iCon = 1:nCons
    sizeMeans=[];
    sizeSE=[];
    for iSize = 1:nSizes
        responsiveTheseTrials = find(h_concat(:,iSize,iCon));
        theseIN=intersect(interNrns,responsiveTheseTrials);

        temp_mean = mean(peak_trough(iSize,iCon,theseIN,1),3,"omitnan");
        temp_se = std(peak_trough(iSize,iCon,theseIN,1),[],3,"omitnan")/sqrt(length(theseIN));

        sizeMeans=[sizeMeans,temp_mean];
        sizeSE=[sizeSE,temp_se];
    end
    errorbar(Sizes,sizeMeans,sizeSE);
    hold on
end


lgd=legend(string(Cons))
title(lgd,'Contrast')
title("Size tuning by contrast, SST cells")
box off
set(gca, 'TickDir', 'out')
xticks(Sizes)
xlabel('Size')
ylabel('df/f peak')
ylim([0 .13])


x0=5;
y0=5;
width=6;
height=3;
set(gcf,'units','inches','position',[x0,y0,width,height])
print('SizeTuning_byContrast.pdf', '-dpdf');
%% Relationship to distance

responsiveTheseTrials = find(h_concat(:,nSizes,nCons)); %getting the cells taht are responsive at the largest size and highest contrast
thesePyr=intersect(pyrCells,responsiveTheseTrials);
theseIN=intersect(interNrns,responsiveTheseTrials);

these_peaks_sst =  squeeze(peak_time(nSizes,nCons,theseIN,3));
these_FWHM_sst =  squeeze(fwhm(nSizes,nCons,theseIN));
these_dists_sst = dists_concat(:,theseIN);

these_peaks_pyr =  squeeze(peak_time(nSizes,nCons,thesePyr,3));
these_FWHM_pyr =  squeeze(fwhm(nSizes,nCons,thesePyr));
these_dists_pyr= dists_concat(:,thesePyr);

figure;
subplot(1,2,1)
swarmchart(these_dists_sst, these_peaks_sst)
set(gca, 'TickDir', 'out')
xlabel('RF-stim distance')
ylabel('Rise time')
title('SST cells')
subplot(1,2,2)
swarmchart(these_dists_pyr, these_peaks_pyr)
set(gca, 'TickDir', 'out')
xlabel('RF-stim distance')
%ylabel('Rise time')
title('Pyr cells')

x0=1;
y0=1;
width=6;
height=2.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
print('Rise_timeByDist', '-dpdf');


figure;
subplot(1,2,1)
swarmchart(these_dists_sst, these_FWHM_sst)
set(gca, 'TickDir', 'out')
xlabel('RF-stim distance')
ylabel('FWHM')
title('SST cells')
subplot(1,2,2)
swarmchart(these_dists_pyr, these_FWHM_pyr)
set(gca, 'TickDir', 'out')
xlabel('RF-stim distance')
%ylabel('Rise time')
title('Pyr cells')

x0=1;
y0=1;
width=6;
height=2.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
print('FWHM_ByDist', '-dpdf');


%% SSI
SSI_peak = squeeze((max(peak_trough(:,3,center,1)))-peak_trough(5,3,:,1));

SSI_mean = squeeze((max(resp_means(:,3,:,1)))-resp_means(5,3,:,1));

scatter(SSI_peak,SSI_mean)
axis square
ylim([0 .5])
xlabel('SSI based on peak')
ylabel('SSI based on mean')

%% splitting SST cells into two groups based on response to small size or not at 80% contrast

%group 1: responsive at ALL sizes. Sum h_concat(iCell,:,4) = nSizes
%group 2: responsive at all sizes ABOVE 15. h_concat(iCell,:,4) == [0 0 1 1
%1]
group1Cells=zeros(1,nCells);
group2Cells=zeros(1,nCells);

group1Vec=[1 1 1 1 1];
group2Vec=[0 0 1 1 1];

for iCell = 1:nCells
    if all(h_concat(iCell,:,4)==group1Vec)
        group1Cells(iCell)=1;
    elseif all(h_concat(iCell,:,4)==group2Vec)
        group2Cells(iCell)=1;
    end
end


sum(group1Cells)
sum(group2Cells)

SST_group1=intersect(interNrns,find(group1Cells));
SST_group2=intersect(interNrns,find(group2Cells));

%% rise time, decay time, and FWHM for SST cells in the two groups above
temp_mean1 = nan(nCons,nSizes,3);
temp_se1 = nan(nCons,nSizes,3);
temp_mean2 = nan(nCons,nSizes,3);
temp_se2 = nan(nCons,nSizes,3);

for iCon = 1:nCons
    for iSize = 1:nSizes %loop through the sizes

        temp_mean1(iCon,iSize,1:2) = mean(peak_time(iSize,iCon,SST_group1,3:4),3,"omitnan");
        temp_se1(iCon,iSize,1:2) = (std(peak_time(iSize,iCon,SST_group1,3:4),[],3,"omitnan"))./length(SST_group1);
        temp_mean1(iCon,iSize,3)=mean(fwhm(iSize,iCon,SST_group1),3,"omitnan");
        temp_se1(iCon,iSize,3)=(std(fwhm(iSize,iCon,SST_group1),[],3,"omitnan"))./length(SST_group1);

        temp_mean2(iCon,iSize,1:2) =  mean(peak_time(iSize,iCon,SST_group2,3:4),3,"omitnan");
        temp_se2(iCon,iSize,1:2) = (std(peak_time(iSize,iCon,SST_group2,3:4),[],3,"omitnan"))./length(SST_group2);
        temp_mean2(iCon,iSize,3)=mean(fwhm(iSize,iCon,SST_group2),3,"omitnan");
        temp_se2(iCon,iSize,3)=(std(fwhm(iSize,iCon,SST_group2),[],3,"omitnan"))./length(SST_group2);


    end
end
temp_mean2(iCon,1:2,:)=NaN;

figure;
subplot(1,3,1)
errorbar(Sizes,temp_mean2(4,:,1),temp_se2(4,:,1),'Color',	"#4e701f");
hold on
errorbar(Sizes,temp_mean1(4,:,1),temp_se1(4,:,1),'Color',	"#3BB806");
%ylim([0,.11])
xticks(Sizes)
xlabel('Size')
ylabel('Seconds')
title('Rise time')
set(gca, 'TickDir', 'out')
box off

subplot(1,3,2)
errorbar(Sizes,temp_mean2(4,:,2),temp_se2(4,:,2),'Color',	"#4e701f");
hold on
errorbar(Sizes,temp_mean1(4,:,2),temp_se1(4,:,2),'Color',	"#3BB806");
%ylim([0,.28])
xticks(Sizes)
xlabel('Size')
title('Decay time')
set(gca, 'TickDir', 'out')
box off

subplot(1,3,3)
errorbar(Sizes,temp_mean2(4,:,3),temp_se2(4,:,3),'Color',	"#4e701f");
hold on
errorbar(Sizes,temp_mean1(4,:,3),temp_se1(4,:,3),'Color',	"#3BB806");
%ylim([0,.28])
xticks(Sizes)
xlabel('Size')
title('FWHM')
set(gca, 'TickDir', 'out')
box off



x0=5;
y0=5;
width=7;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])

print('dynamics_byGroup', '-dpdf');

%% New way of splitting groups

group1Cells=zeros(1,nCells);
group2Cells=zeros(1,nCells);


for iCell = 1:nCells
    if (h_concat(iCell,1,4)==1) & (h_concat(iCell,5,4)==0)
        group1Cells(iCell)=1;
    elseif (h_concat(iCell,1,4)==0) & (h_concat(iCell,5,4)==1)
        group2Cells(iCell)=1;
    end
end
%define a subset of SST cells with med-range retinotopic distances
%SST_midDistance = intersect(interNrns,find(dists_concat >6 & dists_concat < 18));

SST_group1=intersect(interNrns,find(group1Cells));
length(SST_group1)
SST_group2=intersect(interNrns,find(group2Cells));
length(SST_group2)

mean(dists_concat(SST_group1))
mean(dists_concat(SST_group2))

figure;
histogram(dists_concat(SST_group1))
hold on
histogram(dists_concat(SST_group2))
legend('Responsive at 7.5, NOT 120','Responsive at 120, NOT 7.5')
box off
print('RF_distance_byGroup', '-dpdf');

%% how are the groups distributed across depths
mean(depth(SST_group1))
mean(depth(SST_group2))

figure;
histogram(depth(SST_group1))
hold on
histogram(depth(SST_group2))
legend('Responsive at 7.5, NOT 120','Responsive at 120, NOT 7.5')
xlabel('Cortical depth')
box off
print('corticalDepth_byGroup', '-dpdf');


%% rise time, decay time, and FWHM for cells responsive at each size WITHIN groups 1 and 2, 80% contrast
temp_mean1 = nan(nCons,nSizes,3);
temp_se1 = nan(nCons,nSizes,3);
temp_mean2 = nan(nCons,nSizes,3);
temp_se2 = nan(nCons,nSizes,3);

for iCon = 1:nCons
    for iSize = 1:nSizes %loop through the sizes
        
        responsiveTheseTrials = find(h_concat(:,iSize,iCon));
        theseGroup1=intersect(SST_group1,responsiveTheseTrials);
        theseGroup2=intersect(SST_group2,responsiveTheseTrials);

        temp_mean1(iCon,iSize,1:2) = mean(peak_time(iSize,iCon,theseGroup1,3:4),3,"omitnan");
        temp_se1(iCon,iSize,1:2) = (std(peak_time(iSize,iCon,theseGroup1,3:4),[],3,"omitnan"))./length(theseGroup1);
        temp_mean1(iCon,iSize,3)=mean(fwhm(iSize,iCon,theseGroup1),3,"omitnan");
        temp_se1(iCon,iSize,3)=(std(fwhm(iSize,iCon,theseGroup1),[],3,"omitnan"))./length(theseGroup1);


        temp_mean2(iCon,iSize,1:2) =  mean(peak_time(iSize,iCon,theseGroup2,3:4),3,"omitnan");
        temp_se2(iCon,iSize,1:2) = (std(peak_time(iSize,iCon,theseGroup2,3:4),[],3,"omitnan"))./length(theseGroup2);
        temp_mean2(iCon,iSize,3)=mean(fwhm(iSize,iCon,theseGroup2),3,"omitnan");
        temp_se2(iCon,iSize,3)=(std(fwhm(iSize,iCon,theseGroup2),[],3,"omitnan"))./length(theseGroup2);


    end
end

figure;
subplot(1,3,1)
errorbar(Sizes,temp_mean2(4,:,1),temp_se2(4,:,1),'Color',	"#4e701f");
hold on
errorbar(Sizes,temp_mean1(4,:,1),temp_se1(4,:,1),'k');
%ylim([0,.11])
xticks(Sizes)
xlabel('Size')
ylabel('Seconds')
title('Rise time')
set(gca, 'TickDir', 'out')
box off

subplot(1,3,2)
errorbar(Sizes,temp_mean2(4,:,2),temp_se2(4,:,2),'Color',	"#4e701f");
hold on
errorbar(Sizes,temp_mean1(4,:,2),temp_se1(4,:,2),'k');
%ylim([0,.28])
xticks(Sizes)
xlabel('Size')
title('Decay time')
set(gca, 'TickDir', 'out')
box off

subplot(1,3,3)
errorbar(Sizes,temp_mean2(4,:,3),temp_se2(4,:,3),'Color',	"#4e701f");
hold on
errorbar(Sizes,temp_mean1(4,:,3),temp_se1(4,:,3),'k');
%ylim([0,.28])
xticks(Sizes)
xlabel('Size')
title('FWHM')
set(gca, 'TickDir', 'out')
box off



x0=5;
y0=5;
width=7;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])

print('dynamics_respAtEach_byGroup', '-dpdf');

%%
%group 1: responsive at ALL sizes. Sum h_concat(iCell,:,4) = nSizes
%group 2: responsive at all sizes ABOVE 15. h_concat(iCell,:,4) == [0 0 1 1
%1]
group1Cells=zeros(1,nCells);
group2Cells=zeros(1,nCells);

%check whether responsive to the smallest size at the highest contrast

for iCell = 1:nCells
    if h_concat(iCell,1,4)==1
        group1Cells(iCell)=1;
    elseif h_concat(iCell,1,4)~=1
        group2Cells(iCell)=1;
    end
end


SST_group1=intersect(interNrns,find(group1Cells));
length(SST_group1)
SST_group2=intersect(interNrns,find(group2Cells));
length(SST_group2)
%% rise time, decay time, and FWHM for SST cells in the two groups above
temp_mean1 = nan(nCons,nSizes,3);
temp_se1 = nan(nCons,nSizes,3);
temp_mean2 = nan(nCons,nSizes,3);
temp_se2 = nan(nCons,nSizes,3);

for iCon = 1:nCons
    for iSize = 1:nSizes %loop through the sizes

        temp_mean1(iCon,iSize,1:2) = mean(peak_time(iSize,iCon,SST_group1,3:4),3,"omitnan");
        temp_se1(iCon,iSize,1:2) = (std(peak_time(iSize,iCon,SST_group1,3:4),[],3,"omitnan"))./length(SST_group1);
        temp_mean1(iCon,iSize,3)=mean(fwhm(iSize,iCon,SST_group1),3,"omitnan");
        temp_se1(iCon,iSize,3)=(std(fwhm(iSize,iCon,SST_group1),[],3,"omitnan"))./length(SST_group1);

        temp_mean2(iCon,iSize,1:2) =  mean(peak_time(iSize,iCon,SST_group2,3:4),3,"omitnan");
        temp_se2(iCon,iSize,1:2) = (std(peak_time(iSize,iCon,SST_group2,3:4),[],3,"omitnan"))./length(SST_group2);
        temp_mean2(iCon,iSize,3)=mean(fwhm(iSize,iCon,SST_group2),3,"omitnan");
        temp_se2(iCon,iSize,3)=(std(fwhm(iSize,iCon,SST_group2),[],3,"omitnan"))./length(SST_group2);


    end
end
temp_mean2(iCon,1,:)=NaN;

figure;
subplot(1,3,1)
errorbar(Sizes,temp_mean2(4,:,1),temp_se2(4,:,1),'.','MarkerSize',20,'Color',"#4e701f");
hold on
errorbar(Sizes,temp_mean1(4,:,1),temp_se1(4,:,1),'.','MarkerSize',20,'Color',	"#3BB806");
%ylim([0,.11])
xticks(Sizes)
xlabel('Size')
ylabel('Seconds')
title('Rise time')
set(gca, 'TickDir', 'out')
box off

subplot(1,3,2)
errorbar(Sizes,temp_mean2(4,:,2),temp_se2(4,:,2),'.','MarkerSize',20,'Color',	"#4e701f");
hold on
errorbar(Sizes,temp_mean1(4,:,2),temp_se1(4,:,2),'.','MarkerSize',20,'Color',	"#3BB806");
%ylim([0,.28])
xticks(Sizes)
xlabel('Size')
title('Decay time')
set(gca, 'TickDir', 'out')
box off

subplot(1,3,3)
errorbar(Sizes,temp_mean2(4,:,3),temp_se2(4,:,3),'.','MarkerSize',20,'Color',	"#4e701f");
hold on
errorbar(Sizes,temp_mean1(4,:,3),temp_se1(4,:,3),'.','MarkerSize',20,'Color',	"#3BB806");
%ylim([0,.28])
xticks(Sizes)
xlabel('Size')
title('FWHM')
set(gca, 'TickDir', 'out')
box off
sgtitle('Responsive at 7.5 VS NOT responsive at 7.5')


x0=5;
y0=5;
width=7;
height=2;
set(gcf,'units','inches','position',[x0,y0,width,height])

print('dynamicsByGroup', '-dpdf','-bestfit');

%% dynamics vs pref size

%find preferred size for each cell, at 80% contrast, based on peak response
prefSize=NaN(1,nCells);
for iCell = 1:nCells
   [~,prefSize(1,iCell)] = max(peak_trough(:,nCons,iCell,1),[],1); 
end

%plot FWHM vs pref size at 30 degrees for SST cells responsive at 30 degrees
responsiveSST = intersect(find(h_concat(:,3,4)),interNrns); %30 deg, 80% contrast

means_peak=NaN(1,nSizes);
se_peak=NaN(1,nSizes);
means_decay=NaN(1,nSizes);
se_decay=NaN(1,nSizes);
means_FWHM=NaN(1,nSizes);
se_FWHM=NaN(1,nSizes);


for iSize = 1:nSizes
    prefThisSize=find(prefSize==iSize);
    theseSST=intersect(prefThisSize,responsiveSST);

    means_peak(iSize)=mean(squeeze(peak_time(3,4, theseSST,3)),"omitmissing");
    se_peak(iSize) = (std(squeeze(peak_time(3,4, theseSST,3)),[],"omitmissing"))/length(theseSST);

    means_decay(iSize)=mean(squeeze(peak_time(3,4, theseSST,4)),"omitmissing");
    se_decay(iSize) = (std(squeeze(peak_time(3,4, theseSST,4)),[],"omitmissing"))/length(theseSST);

    means_FWHM(iSize)=mean(squeeze(fwhm(iSize,iCon,theseSST)),"omitmissing");
    se_FWHM(iSize) = (std(squeeze(fwhm(iSize,iCon,theseSST)),[],"omitmissing"))/length(theseSST);

end


figure
subplot(1,3,1)
errorbar(means_peak,se_peak,'.','MarkerSize',20,'Color',	"#4e701f")
ylabel('Rise')
xlabel('Pref size')
xlim([.5 5.5])
xticks([1,2,3,4,5])
xticklabels(Sizes)
box off
set(gca, 'TickDir', 'out')

subplot(1,3,2)
errorbar(means_decay,se_decay,'.','MarkerSize',20,'Color',	"#4e701f")
ylabel('Decay')
xlabel('Pref size')
xlim([.5 5.5])
xticks([1,2,3,4,5])
xticklabels(Sizes)
box off
set(gca, 'TickDir', 'out')

subplot(1,3,3)
errorbar(means_FWHM,se_FWHM,'.','MarkerSize',20,'Color',	"#4e701f")
ylabel('FWHM')
xlabel('Pref size')
xlim([.5 5.5])
xticks([1,2,3,4,5])
xticklabels(Sizes)
box off
set(gca, 'TickDir', 'out')

sgtitle('30 degree size, 80% contrast')

x0=5;
y0=5;
width=7;
height=2;
set(gcf,'units','inches','position',[x0,y0,width,height])


print('dynamicsAt30_byPrefSize', '-dpdf','-bestfit');
%% seeing how groups of cells are clustered among/within mice
SST_group1_logical = group1Cells.*interNrns_concat;
SST_group2_logical = group2Cells.*interNrns_concat;

mouseIndx = NaN(2,nExpt);
for iExpt = 1:nExpt
    mouseIndx(1,iExpt)=min(find(exptNumber==iExpt));
    mouseIndx(2,iExpt)=max(find(exptNumber==iExpt));
end

figure
for iExpt = 1:nExpt
    subplot(3,6,iExpt)
    plot(SST_group1_logical(mouseIndx(1,iExpt):mouseIndx(2,iExpt)),'k')
    hold on
    plot(SST_group2_logical(mouseIndx(1,iExpt):mouseIndx(2,iExpt)),'color',"#3BB806")
    yticks([0,1])
end

x0=5;   
y0=5;
width=9;
height=6;
set(gcf,'units','inches','position',[x0,y0,width,height])
sgtitle('Responsive at 7.5, NOT 120 VS Responsive at 120, NOT 7.5')


print('SST_group_distribution2', '-dpdf','-bestfit');